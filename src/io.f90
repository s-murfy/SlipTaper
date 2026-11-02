module io
use generic 

implicit none

contains 
!##############################################################
subroutine read_input(fname,nx,nz,slipping_cell,taper_dist,slip_rake,surface,output_level)
character(:),allocatable,intent(out) :: fname
character(len=256) :: name
integer,intent(out) :: nx,nz,slipping_cell
real(pr)::rake ,surface_depth,taper
real(pr),intent(out) :: slip_rake,surface,taper_dist
integer,intent(out) :: output_level
integer :: nx_cells,nz_cells,cell,verbose 


NAMELIST / MESHFILE / name
NAMELIST / MESHGRID / nx_cells,nz_cells,taper
NAMELIST / CELLCHOICE / cell,surface_depth, rake
NAMELIST / OUTPUT / verbose

! set default value
verbose = 1

open(11,file='input.txt',status='old',action='read')
read(11,nml=MESHFILE)
read(11,nml=MESHGRID)
read(11,nml=CELLCHOICE)
read(11,nml=OUTPUT)

allocate(character(len=len_trim(name)) :: fname)
fname = trim(name) 
nx = nx_cells+1  ! an additional point is need for defining nodes associated with cells
nz = nz_cells+1  ! an additional point is need for defining nodes associated with cells
! nx = 11
! nz = 11

slipping_cell = cell
slip_rake = rake 
surface = surface_depth
output_level = verbose 
taper_dist = taper 

close(11)
return 
end subroutine read_input
!##############################################################
subroutine read_quadrangles(fname,lon,lat,z,nquads,output_level)
integer :: cell_id, iostat, line_count
integer :: i
integer, intent(out) :: nquads
real(pr) :: lon1,lat1,d1,lon2,lat2,d2,lon3,lat3,d3,lon4,lat4,d4
real(pr),allocatable,dimension(:,:),intent(out) :: lon,lat,z
character(len=256) :: line
character(:),ALLOCATABLE,intent(in) :: fname
character(len=*), parameter :: fmt ="(12E14.6)"
integer,intent(in) :: output_level

open(12,file=fname, status='old', action='read', iostat=iostat)

line_count = 0
if (iostat /= 0) then 
    write(*,*) 'error reading file'
    stop
endif

do 
    read(12,'(A)',iostat=iostat) line
    if (iostat /= 0) exit
    line_count = line_count+1
enddo
close(12)

! number of quadrangles on fault
if(output_level == 2) write(*,*) 'number of quadrangles in file...',line_count
nquads = line_count
allocate(lon(nquads,4),lat(nquads,4),z(nquads,4))

open(12,file=fname)
do i = 1,nquads
    read(12,*) lon1,lat1,d1,lon2,lat2,d2,lon3,lat3,d3,lon4,lat4,d4

    lon(i,1:4) = (/ lon1,lon2,lon3,lon4/)
    lat(i,1:4) = (/ lat1,lat2,lat3,lat4/)
    z(i,1:4) =   (/ d1,d2,d3,d4 /)
enddo

return 
end subroutine read_quadrangles 
!##############################################################
subroutine convert2geo(amesh,latref,lonref,geomesh)
use georef
use LAT_mesh
type(mesh),intent(in) :: amesh
real(pr),intent(in) :: latref,lonref
type(mesh),intent(out) :: geomesh
real(pr) :: lat, lon 
integer :: i 

geomesh%Nnodes = amesh%Nnodes
geomesh%Ncells = amesh%Ncells

allocate(geomesh%cell(geomesh%Ncells,3))
do i = 1,geomesh%Ncells
    geomesh%cell(i,:) = amesh%cell(i,:)
enddo 

allocate(geomesh%px(geomesh%Nnodes),geomesh%py(geomesh%Nnodes),geomesh%pz(geomesh%Nnodes))
geomesh%pz = amesh%pz
do i = 1,geomesh%Nnodes
    lat = 0._pr
    lon = 0._pr
    call convgeoutm(lat,lon,latref,lonref,amesh%px(i),amesh%py(i),-1)
    geomesh%px(i) = lat
    geomesh%py(i) = lon
enddo 


return 
end subroutine convert2geo
!##############################################################
end module io 
