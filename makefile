# Makefile
#
FC = gfortran 
#FC = ifx 

OBJ = generic.o georef.o LAT_mesh.o lists.o distance.o io.o slip.o main.o 


# debugging options 
#
#FCFLAGS = -O -g -fbounds-check -ffpe-trap=invalid -fcheck=all -fbacktrace
# FCFLAGS = -O -g -traceback -implicitnone #-check all 
FCFLAGS = -O2 

.SUFFIXES: .f90


all: taper.x

# objects
#
#
generic.o : generic.f90
	$(FC) $(FCFLAGS) -c $<

LAT_mesh.o : LAT_mesh.f90
	$(FC) $(FCFLAGS) -c $<

georef.o : georef.f90
	$(FC) $(FCFLAGS) -c $<

lists.o : lists.f90	
	$(FC) $(FCFLAGS) -c $<

distance.o : distance.f90 
	$(FC) $(FCFLAGS) -c $<

io.o : io.f90
	$(FC) $(FCFLAGS) -c $<

slip.o : slip.f90
	$(FC) $(FCFLAGS) -c $<

main.o : main.f90
	$(FC) $(FCFLAGS) -c $<

taper.x :$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod

