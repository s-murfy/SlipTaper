module io
use generic 

implicit none

contains 
!##############################################################
subroutine read_all_nodes(lon,lat,z,nquads)
integer :: cell_id, iostat, line_count
integer :: i,nquads
real(pr) :: lon1,lat1,d1,lon2,lat2,d2,lon3,lat3,d3
real(pr) :: lon4,lat4,d4,lon5,lat5,d5,lon6,lat6,d6
real(pr),allocatable,dimension(:,:) :: lon,lat,z
character(len=256) :: line
character(len=*), parameter :: fmt ="(I2,F8.4,F8.4,F8.0,F8.4,F8.4,F8.0,F8.4,F8.4,F8.0)"

open(12,file='puysegur2_25x25_lon_lat_depth.txt', status='old', action='read', iostat=iostat)

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
nquads = int(line_count/2.)
allocate(lon(nquads,4),lat(nquads,4),z(nquads,4))

open(12,file='puysegur2_25x25_lon_lat_depth.txt')
do i = 1,nquads
    read(12,*) cell_id,lon1,lat1,d1,lon2,lat2,d2,lon3,lat3,d3
    read(12,*) cell_id,lon4,lat4,d4,lon5,lat5,d5,lon6,lat6,d6

    lon(i,1:4) = (/ lon1,lon2,lon3,lon5/)
    lat(i,1:4) = (/ lat1,lat2,lat3,lat5/)
    z(i,1:4) =   (/ d1,d2,d3,d5 /)
enddo

return
end subroutine read_all_nodes
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

! open(12,file='puysegur2_quadrilateral_coords.dat', status='old', action='read', iostat=iostat)
! open(12,file='southamerica_mod_quadrilateral_coords.dat', status='old', action='read', iostat=iostat)
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
! !##############################################################
! subroutine read_nodes(lon,lat,z)
! integer :: int
! real(pr) :: lon1,lat1,d1,lon2,lat2,d2,lon3,lat3,d3
! real(pr) :: lon4,lat4,d4,lon5,lat5,d5,lon6,lat6,d6
! real(pr),dimension(4) :: lon,lat,z
! character(len=*), parameter :: fmt ="(I2,F8.4,F8.4,F8.0,F8.4,F8.4,F8.0,F8.4,F8.4,F8.0)"

! open(12,file='puysegur2_25x25_lon_lat_depth.txt')
! read(12,*) int,lon1,lat1,d1,lon2,lat2,d2,lon3,lat3,d3
! read(12,*) int,lon4,lat4,d4,lon5,lat5,d5,lon6,lat6,d6

! ! next cell
! ! read(12,*) int,lon1,lat1,d1,lon2,lat2,d2,lon3,lat3,d3
! ! read(12,*) int,lon4,lat4,d4,lon5,lat5,d5,lon6,lat6,d6

! close(12)
! ! lon(1:4) = (/ 165.7354,  166.0931,  166.4346,  166.2403/)
! ! lat(1:4) = (/ -46.6238, -46.2828,  -46.3874,  -46.7742/)
! lon(1:4) = (/ 166.0931,   165.7354,  166.4346,  166.2403/)
! lat(1:4) = (/ -46.2828, -46.6238,  -46.3874,  -46.7742/)
! z(1:4) = (/-1.4726e4, -1.3175e4, -40000.0, -40000.0/)
! ! lon(1:4) = (/ lon1,lon2,lon3,lon5/)
! ! lat(1:4) = (/ lat1,lat2,lat3,lat5/)
! ! z(1:4) =   (/ d1,d2,d3,d5 /)
! end subroutine read_nodes
! !##############################################################
!##############################################################
function get_lengths(x,y,z) 
real(pr),dimension(4),intent(in) :: x,y,z
real(pr),dimension(2) :: get_lengths
real(pr) :: len_up,len_down

len_up = sqrt((x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2)
len_down = sqrt((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)

!write(*,*) len_up/1000.,len_down/1000.

get_lengths = (/ len_up,len_down /) 
end function get_lengths
!##############################################################
pure function get_widths(x,y,z)
real(pr),dimension(4),intent(in) :: x,y,z
real(pr),dimension(2):: get_widths
real(pr) :: width_s,width_e

! width_s = sqrt((x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2)
! width_e = sqrt((x(4)-x(1))**2+(y(4)-y(1))**2+(z(4)-z(1))**2)

width_s = sqrt((x(3)-x(1))**2+(y(3)-y(1))**2+(z(3)-z(1))**2)
width_e = sqrt((x(4)-x(2))**2+(y(4)-y(2))**2+(z(4)-z(2))**2)

!write(*,*) 'width...',width_s/1000,width_e/1000.

get_widths = (/ width_s,width_e /) 

end function get_widths
!##############################################################
subroutine tapered_slip(amesh,dist,central_dist,slip)
use LAT_mesh
real(pr),allocatable, dimension(:,:),intent(in) :: dist
type(mesh),intent(in) :: amesh
real(pr), allocatable, dimension(:),intent(out) :: slip 
real(pr), allocatable, dimension(:),intent(out) :: central_dist
integer :: i, central_node
real(pr) :: d1,d2,d3,taper_s,taper_e,taper_length,mean_d
real(pr),parameter :: kilo = 1000._pr 

! Calculate tapered slip from centre of cell 
central_node = anint(real(amesh%Nnodes)/2.)
allocate(central_dist(amesh%Nnodes))
central_dist = dist(central_node,:)

allocate(slip(amesh%Ncells))
slip  = 0._pr
taper_s = 5._pr*kilo
taper_e = 12.5_pr*kilo 
taper_length = taper_e-taper_s
do i = 1,amesh%Ncells
    d1 = central_dist(amesh%cell(i,1))
    d2 = central_dist(amesh%cell(i,2))
    d3 = central_dist(amesh%cell(i,3))

    mean_d = (d1+d2+d3)/3._pr

    if (mean_d <= taper_s) then
       slip(i)  = 1._pr 
    elseif ((mean_d > taper_s).and.(mean_d <= taper_e)) then 
        slip(i) = 1._pr - (mean_d-taper_s)/taper_length
    endif 
enddo 

return 
end subroutine tapered_slip
!##############################################################
subroutine tapered_slip_rectangle(amesh,dist,central_dist,rect_mesh,slip)
use LAT_mesh
real(pr),allocatable, dimension(:,:),intent(in) :: dist
type(mesh),intent(in) :: amesh
type(mesh),intent(out) :: rect_mesh

real(pr), allocatable, dimension(:),intent(out) :: slip 
real(pr), allocatable, dimension(:),intent(out) :: central_dist
integer :: i, inx, central_node
real(pr) :: d1,d2,d3,d4, taper_s,taper_e,taper_length,mean_d
real(pr),parameter :: kilo = 1000._pr 

! Calculate tapered slip from centre of cell 
central_node = anint(real(amesh%Nnodes)/2.)
allocate(central_dist(amesh%Nnodes))
central_dist = dist(central_node,:)

rect_mesh%Ncells = int(real(amesh%Ncells,pr)/2.) 
allocate(rect_mesh%cell(rect_mesh%Ncells,4))
allocate(slip(rect_mesh%Ncells))
rect_mesh%cell = 0._pr

rect_mesh%Nnodes = amesh%Nnodes 
allocate(rect_mesh%px(rect_mesh%Nnodes),rect_mesh%py(rect_mesh%Nnodes),rect_mesh%pz(rect_mesh%Nnodes))
rect_mesh%px = amesh%px
rect_mesh%py = amesh%py
rect_mesh%pz = amesh%pz

slip  = 0._pr
taper_s = 5._pr*kilo
taper_e = 12.5_pr*kilo 
taper_length = taper_e-taper_s
i = 1
inx = 1
do while(i <= amesh%Ncells)
    d1 = central_dist(amesh%cell(i,1))
    d2 = central_dist(amesh%cell(i,2))
    d3 = central_dist(amesh%cell(i,3))
    d4 = central_dist(amesh%cell(i+1,3))

    rect_mesh%cell(inx,:) =  (/amesh%cell(i,1),amesh%cell(i,2),amesh%cell(i,3),amesh%cell(i+1,3) /)

    mean_d = (d1+d2+d3+d4)/4._pr
    if (mean_d <= taper_s) then
       slip(inx)  = 1._pr 
    elseif ((mean_d > taper_s).and.(mean_d <= taper_e)) then 
        slip(inx) = 1._pr - (mean_d-taper_s)/taper_length
    endif 
    inx = inx+1
    i = i+2
enddo 



return 
end subroutine tapered_slip_rectangle
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
