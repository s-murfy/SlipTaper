"""Put a quadrangle mesh into the cell order that taper.x needs.

taper.x only works if the cells in the mesh file are in "raster" order and use a
consistent node convention:

  * nodes 1,4 are on the deep edge of a cell and nodes 2,3 on the shallow edge,
  * the cell that follows in the file is the next cell along strike, i.e. its
    nodes 1,2 are nodes 4,3 of the previous cell,
  * the file lists one complete row of cells (along strike) after another,
    starting with the row at the surface and working down dip.

Shuffling the lines of a mesh that works makes taper.x fail for many cells
(index 0 out of bounds in LAT_mesh.f90 / slip.f90), so the file order matters
even when every cell on its own is fine. 

Usage:  python order_mesh.py mesh_in.dat mesh_out.dat
"""
import collections
import sys

import numpy as np


def raster_order(quads, decimals=3, return_source=False):
    """Re-order a mesh into the order described above.

    quads : (ncell, 4, 3) array of [lon, lat, z] for the nodes of each quadrangle
    Returns the re-ordered (ncell, 4, 3) array and the number of cells in each row.
    With return_source=True it also returns, for each output cell, the index of the
    input cell it came from (used by describe_changes).
    """
    key = lambda p: tuple(np.round(p[:2], decimals))
    n = len(quads)

    # the two highest nodes (z closest to the surface) form the shallow edge
    shallow, deep = [], []
    for q in quads:
        by_z = np.argsort(-q[:, 2])
        shallow.append({key(q[k]) for k in by_z[:2]})
        deep.append({key(q[k]) for k in by_z[2:]})

    node_cells = collections.defaultdict(list)
    for i, q in enumerate(quads):
        for p in q:
            node_cells[key(p)].append(i)

    # along-strike neighbours share one shallow node and one deep node
    neighbours = [[] for _ in range(n)]
    for i in range(n):
        around = {c for k in shallow[i] | deep[i] for c in node_cells[k]} - {i}
        neighbours[i] = [j for j in around
                         if len(shallow[i] & shallow[j]) == 1 and len(deep[i] & deep[j]) == 1]
    if any(len(nb) > 2 for nb in neighbours):
        raise ValueError('a cell has more than 2 along-strike neighbours; is the mesh a single grid?')

    # a row is a chain of cells, walk it from one of its two ends
    seen = [False] * n
    rows = []
    for start in range(n):
        if seen[start] or len(neighbours[start]) == 2:
            continue
        row = [start]
        seen[start] = True
        while True:
            nxt = [j for j in neighbours[row[-1]] if not seen[j]]
            if not nxt:
                break
            row.append(nxt[0])
            seen[nxt[0]] = True
        rows.append(row)
    if not all(seen):
        raise ValueError('some cells are not part of an open row of cells')
    if any(len(r) < 2 for r in rows):
        raise ValueError('found a row with a single cell; cannot tell its strike direction')

    # shallowest row first, and every row running in the same direction as the first
    centre = lambda i: quads[i][:, :2].mean(axis=0)
    rows.sort(key=lambda r: -np.mean([quads[i][:, 2].max() for i in r]))
    reference = centre(rows[0][-1]) - centre(rows[0][0])
    if reference[1] < 0:                      # first row runs south -> north
        rows[0].reverse()
        reference = -reference
    for r in rows[1:]:
        if np.dot(centre(r[-1]) - centre(r[0]), reference) < 0:
            r.reverse()

    ordered, source = [], []
    for row in rows:
        for pos, i in enumerate(row):
            pts = {key(p): p for p in quads[i]}
            if pos > 0:                       # near side is shared with the previous cell
                shared = {key(p) for p in quads[row[pos - 1]]}
                near_s, near_d = (shallow[i] & shared).pop(), (deep[i] & shared).pop()
                far_s, far_d = (shallow[i] - {near_s}).pop(), (deep[i] - {near_d}).pop()
            else:                             # first cell: far side is shared with the next cell
                shared = {key(p) for p in quads[row[pos + 1]]}
                far_s, far_d = (shallow[i] & shared).pop(), (deep[i] & shared).pop()
                near_s, near_d = (shallow[i] - {far_s}).pop(), (deep[i] - {far_d}).pop()
            ordered.append([pts[near_d], pts[near_s], pts[far_s], pts[far_d]])
            source.append(i)

    if return_source:
        return np.array(ordered), [len(r) for r in rows], np.array(source)
    return np.array(ordered), [len(r) for r in rows]


def describe_changes(quads_in, quads_out, source, row_lengths, max_listed=10):
    """Text report of what raster_order changed. Cell numbers are 1-based line numbers."""
    n = len(quads_out)
    positions = np.arange(n)
    moved = positions[source != positions]
    # a cell's nodes were re-ordered if its lon/lat corners differ, node by node, from the input cell
    renumbered = positions[np.any(quads_in[source][:, :, :2] != quads_out[:, :, :2], axis=(1, 2))]

    def listed(idx):
        txt = ', '.join(str(i + 1) for i in idx[:max_listed])
        return txt + (f', ... ({len(idx) - max_listed} more)' if len(idx) > max_listed else '')

    lines = [f'{n} cells in {len(row_lengths)} rows: {row_lengths}']
    if len(moved) == 0 and len(renumbered) == 0:
        lines.append('The mesh was already in the correct order: nothing was changed.')
    else:
        lines.append(f'Cells moved to a different position in the file: {len(moved)} of {n}'
                     + (f' (new positions {listed(moved)})' if len(moved) else ''))
        lines.append(f'Cells whose node order was changed: {len(renumbered)} of {n}'
                     + (f' (new positions {listed(renumbered)})' if len(renumbered) else ''))
        lines.append('Cell numbers in input.txt (cell=) refer to the new order.')
    return '\n'.join(lines)


if __name__ == '__main__':
    quads_in = np.loadtxt(sys.argv[1]).reshape(-1, 4, 3)
    quads_out, row_lengths, source = raster_order(quads_in, return_source=True)
    print(describe_changes(quads_in, quads_out, source, row_lengths))
    np.savetxt(sys.argv[2], quads_out.reshape(-1, 12), fmt='%.6f')
