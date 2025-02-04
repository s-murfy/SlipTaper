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

[Trilateration]: https://github.com/andherit/trilateration
[GeoREF]:https://github.com/andherit/georef
