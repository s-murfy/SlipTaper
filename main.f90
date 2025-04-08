program taper_slip
use generic
use io 
use georef
use LAT_mesh
use mod_slip
implicit none 
integer :: i,j,nx,nz 
real(pr) :: latref,lonref
real(pr),dimension(4) :: z,x,y
real(pr),allocatable,dimension(:,:) :: fault_lon,fault_lat,fault_z
real(pr),allocatable,dimension(:,:) :: quad_lon,quad_lat,quad_z,quad_x,quad_y
integer,allocatable,dimension(:,:,:) :: edge_x,edge_z
integer,allocatable,dimension(:):: tri2quad
! real(pr),allocatable,dimension(:):: bc_dist,edge_dist,corner_dist
real(pr),allocatable, dimension(:) :: slipping_dist
real(pr),allocatable,dimension(:,:):: bc_dist
real(pr) :: accuracy
real(pr), dimension(2) :: flenghts,fwidths
real(pr), allocatable, dimension(:) :: central_dist,slip,rect_slip
real(pr),dimension(2) :: cell_lengths,cell_widths
character(30) :: out_file
integer :: central_node
integer:: nfault,nquad,slipping_cell,original_slipping_cell
type(mesh) :: amesh,geomesh,rect_mesh,qmesh,fault_mesh
integer :: anticlockwise
character(:),allocatable :: fname 
real(pr) :: slip_rake,surface,taper_dist
integer :: output_level
logical :: edge_exist

! reaad input file 
call read_input(fname,nx,nz,slipping_cell,taper_dist,slip_rake,surface,output_level)

original_slipping_cell = slipping_cell

! read in coordinates from file
call read_quadrangles(fname, fault_lon,fault_lat,fault_z,nfault,output_level)

! create mesh for whole fault 
accuracy = 0.01_pr
call create_quad_mesh(fault_lon,fault_lat,fault_z,nfault,accuracy,fault_mesh,output_level)

! choose cells connected to slipping cell
!slipping_cell = 6 !2 ! this is the slipping cell in the original mesh
                   ! it will be updated in find_slipping_cells


call find_slipping_cells(slipping_cell,fault_mesh,quad_lon,quad_lat,quad_z,nquad,output_level)
if (output_level >= 2) write(*,*) 'new slipping cell id............',slipping_cell
allocate(quad_x(nquad,4),quad_y(nquad,4))

! set first point as the reference point 
latref=quad_lat(1,1) 
lonref=quad_lon(1,1)
! transfrom all points to utm
subfault: do i = 1,nquad
    inside_cell: do j = 1,4
                call convgeoutm(quad_lat(i,j),quad_lon(i,j),latref,lonref,quad_x(i,j),quad_y(i,j),1)
        enddo inside_cell
        ! write(*,*) 'quadrangle cell no', i
        ! anticlockwise = check_anticlockwise(quad_x(i,:),quad_y(i,:),4)
        ! if (anticlockwise > 0) write(*,*) 'anti_clockwise'
        ! if (anticlockwise < 0) write(*,*) 'clockwise'
enddo subfault
 
 
! create quadrangle mesh
accuracy = 1._pr
call create_quad_mesh(quad_x,quad_y,quad_z,nquad,accuracy,qmesh,output_level)



! create triangle mesh from quadrangle mesh
call create_tri_mesh(qmesh,nx,nz,amesh,edge_x,edge_z,tri2quad,output_level)



! ! calculate distance between all points on mesh
if (output_level >= 2) write(*,*) 'Calculate border distances '
! call boundary_dist(slipping_cell,qmesh,amesh,edge_x,edge_z,surface,bc_dist,edge_dist,corner_dist,edge_exist)
call boundary_distV2(qmesh,amesh,edge_x,edge_z,surface,bc_dist)

! ! calc tapered slip 
!call tapered_slip(amesh,dist,central_dist,slip)
if (output_level >= 2) write(*,*) 'calculate slip on fault...'
! call add_slip(slipping_cell,qmesh,amesh,edge_x,edge_z,tri2quad,bc_dist,edge_dist,corner_dist,edge_exist,taper_dist,slip)
call add_slipV2(slipping_cell,qmesh,amesh,edge_x,edge_z,bc_dist,taper_dist,slip)


! ! calculate tapered slip and transform mesh to rectilinear cells 
! call tapered_slip_rectangle(amesh,dist,central_dist,rect_mesh,rect_slip)

! ! output slip in utm coordinates 
if (output_level == 2) write(*,*) 'output slip'

allocate(slipping_dist(amesh%Ncells))
slipping_dist = bc_dist(slipping_cell,:)
if (output_level >= 1) then 
        write(out_file,'(a15,i3.3)') 'slip_taper_utm_',original_slipping_cell
        ! call write_out(amesh,bc_dist,slip,out_file)     ! writing out distance from boundary 
        

        call write_out(amesh,slipping_dist,slip,out_file)   ! writing out corner distance 




        if (output_level == 2) write(*,*) 'output original quadrangles'
!call write_out(amesh,central_dist,slip,out_file)
        call write_out_rect(qmesh,out_file) ! possible to activate nodal and cell values outputs if needed
! call write_out_rect(rect_mesh,central_dist,rect_slip,out_file)  
endif 


! ! convert mesht to lat/lon
if (output_level >= 2) write(*,*) 'convert to lat/lon'
call convert2geo(amesh,latref,lonref,geomesh)


if (output_level >= 2) write(*,*) 'output slip in lat/lon coordinates'

write(out_file,'(a15,i3.3)') 'slip_taper_geo_',original_slipping_cell
call write_out(geomesh,slipping_dist,slip,out_file)


call hysea_output(geomesh,slip,slip_rake,out_file)

stop 
end program taper_slip
