# SlipTaper
Code to discretise a quadrangle into triangles producing a smooth slip distribution that tapers to zero at boundary.

## Description 
``SlipTaper`` take a coarse quadrangles mesh, descritizes it to a higher specified resolution of triangles. A smooth slip distribution is then appiled to the specified quadrangle which is smoothly tapered into adjoining quadrangles which have also been meshed at the given higher resolution. Distances across the mesh are calculated using [Trilateration] while the coversion between geographical coorindates and UTM are performed with [GeoREF]


## How to compile and run
Compliation requires a fortran compiler. If using `gfortran` open the `makefile` comment line 4 (FC=ifx) and uncomment line 3 (FC=gfortran). 
To compile by typing ``make `` to create the execuitable `taper.x`. To run type `.\taper.x` in the folder containing the `input.txt` and quadrangle mesh file.  

### Input files 
Two input files are required by the code. The first is a file called ``input.txt`` in this file general parameters for SlipTaper are provided (e.g. number of subcells within quadrangles, name of quadrangle file, length of taper, etc).  The second is the name of the coarse quadrangle mesh file that will be used to calculate the finer triangular grid onto. The mesh file will contain the nodes of the quadrangles for each cell on each line in geographical coordinates. 

### Output file 
Output files contain a triangular mesh in geographical coordinates in vtk format which can be opened using *Paraview* and as txt file that can be read by *HySEA*. Setting the verbose in the input.txt to 1 will generate vtk files containing the intermediate meshes in UTM coordinates. All output files contain a number which corresponds to the quadrangle from the original mesh which the slip is centred on. 

## Repository layout
- `src/` : Fortran source and `makefile`.
- `example/` : an example `input.txt`, quadrangle mesh files (`southamerica_quadrilateral_coords_order.dat`) and a notebook (`plot_slip_results.ipynb`) for checking meshes.

## Quick start
```bash
cd src
make                          # builds taper.x (edit FC in the makefile if not using gfortran)
cd ../example
../src/taper.x                # reads ./input.txt and the mesh file it names
```
`taper.x` reads `input.txt` from the **current directory**, and the mesh file named in it is also looked up relative to the current directory, so run it from the folder holding both (e.g. `example/`).

## Input file (`input.txt`)
The input uses Fortran namelists; lines starting with `#` are comments.

```
&MESHFILE   name='southamerica_quadrilateral_coords_order.dat' /
&MESHGRID   nx_cells=10, nz_cells=10, taper=15.d3 /
&CELLCHOICE cell=144, surface_depth=-100.d0, rake=90.d0 /
&OUTPUT     verbose=1 /
```

| Namelist | Parameter | Meaning |
|---|---|---|
| `MESHFILE` | `name` | Quadrangle mesh file (see below). |
| `MESHGRID` | `nx_cells`, `nz_cells` | Number of subdivisions along strike and down dip within each quadrangle. |
| `MESHGRID` | `taper` | Distance (m) over which slip tapers to zero. |
| `CELLCHOICE` | `cell` | Quadrangle (line number in the mesh file) the slip is centred on. |
| `CELLCHOICE` | `surface_depth` | If a node on a quadrangle edge is shallower than this value (m), that edge is treated as the surface and no taper is applied there. |
| `CELLCHOICE` | `rake` | Rake (degrees) written to the HySEA file. |
| `OUTPUT` | `verbose` | 0: final lon/lat outputs only; 1: also intermediate UTM vtk files; 2: also print progress to screen. |

The four namelists must appear in this order.

## Mesh file format
One quadrangle per line, with its four corner nodes given in order as
`lon1 lat1 depth1  lon2 lat2 depth2  lon3 lat3 depth3  lon4 lat4 depth4`
(longitude/latitude in degrees, depth in metres, negative downwards).

```
   shallow (surface) edge
   node 2 ----------- node 3
     |                   |        along strike  ->
     |      cell         |
     |                   |
   node 1 ----------- node 4
   deep edge
```

### Requirements on the mesh
The quadrangles form a regular grid (rows along strike, columns down dip) and the triangle mesh is built cell by cell, reusing the nodes already created on shared edges. This only works if the mesh file is laid out as follows; `taper.x` will otherwise stop with an error.

1. **Neighbouring cells share nodes exactly.** Shared corners have identical coordinates, so every quadrangle edge is shared by at most two cells.
2. **Every cell uses the same node convention.** Nodes 2-3 are the shallow edge and nodes 1-4 the deep edge, so 1→4 runs along strike and 1→2 runs down dip in all cells. An edge shared by two cells must be a strike edge (1-4 or 2-3) in both or a dip edge (1-2 or 3-4) in both. The surface test (`surface_depth`) is only applied to the 2-3 edge.
3. **Cells are listed row by row.** The file lists a complete row of cells along strike, then the next row, starting with the row at the surface. Within a row each cell is followed by its along-strike neighbour, so its nodes 1,2 are nodes 4,3 of the previous cell.

Putting the cells in the wrong order is enough to break the mesh even if every cell is individually correct: shuffling the lines of a working mesh will make `taper.x` fail.

### Preparing a mesh
`example/order_mesh.py` takes any mesh of this kind and rewrites it to meet the requirements above:
```bash
python example/order_mesh.py my_mesh.dat my_mesh_ordered.dat
```
It puts the shallow edge on nodes 2-3, makes the strike direction the same in every cell, and orders the cells row by row. It prints the number of cells in each row; check that this is what you expect (for example 4 rows of 165 cells). It stops with an error if the mesh is not a single grid. Cell numbers in `input.txt` (`cell=`) refer to the order in the file you give to `taper.x`, which is the re-ordered file. 

## Output files
Written to the current directory; `NNN` is the chosen cell number (3 digits, zero padded).
- `slip_taper_geo_NNN.vtk` : triangular mesh with slip in lon/lat, for Paraview.
- `slip_taper_geo_NNN_hysea.dat` : the same slip distribution (with the rake) in the format read by HySEA.
- With `verbose=1`: `slip_taper_utm_NNN.vtk` (triangular mesh in UTM) and a `*_rec.vtk` file with the original quadrangles.

To produce a slip distribution for a different quadrangle, change `cell` in `input.txt` and rerun.

[Trilateration]: https://github.com/andherit/trilateration
[GeoREF]:https://github.com/andherit/georef

## Reference 
If you use this code please reference the following paper which it is based on: 

> Herrero, A. and Murphy, S., 2018. Self-similar slip distributions on irregular shaped faults. *Geophysical Journal International*, 213(3), pp.2060–2070, https://doi.org/10.1093/gji/ggy104

**Acknowledgements:** The plotting notebook, order_mesh.py, the mesh error checks in LAT_mesh.f90 and parts of this README were drafted with the assistance of Claude.

