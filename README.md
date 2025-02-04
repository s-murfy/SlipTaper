# SlipTaper
Code to discretise a quadrangle into triangles producing a smooth slip distribution that tapers to zero at boundary.

##Description 
This code 

## How to compile 

### Input files 
Two input files are required by the code. The first is a file called ``input.txt`` in this file general parameters for SlipTaper are provided (e.g. number of subcells within quadrangles, name of quadrangle file, length of taper, etc).  The second is the name of the coarse quadrangle mesh file that will be used to calculate the finer triangular grid onto. 

### Output file 

Distances across the mesh are calculated using [Trilateration] while the coversion between geographical coorindates and UTM are performed with [GeoREF]

[Trilateration]: https://github.com/andherit/trilateration
[GeoREF]:https://github.com/andherit/georef
