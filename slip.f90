module mod_slip
use generic 
use LAT_mesh
implicit none 

contains

!====================================================================
subroutine boundary_dist(slipping_cell,qmesh,amesh,edge_x,edge_z,surface,bc_dist)
use distance
integer,allocatable,dimension(:,:,:),intent(in) :: edge_x,edge_z
type(mesh),intent(in) :: qmesh,amesh
real(pr),allocatable,dimension(:) :: bc_dist
real(pr),intent(in) :: surface
real(pr) :: i_value
integer,intent(in) :: slipping_cell
integer :: qnode_1,qnode_2,qnode_3,qnode_4
logical :: fast ,not_surface
integer :: i

qnode_1 = qmesh%cell(slipping_cell,1)
qnode_2 = qmesh%cell(slipping_cell,2)
qnode_3 = qmesh%cell(slipping_cell,3)
qnode_4 = qmesh%cell(slipping_cell,4)



if (allocated(bc_dist)) deallocate(bc_dist)
allocate(bc_dist(amesh%Nnodes))
bc_dist = infinity
i_value = 10000000._pr  ! use a large number to generate for the generation of plane wavefront
!set sides of quadrangle to zero 
bc_dist(edge_z(qnode_1,qnode_2,:)) = i_value
bc_dist(edge_z(qnode_3,qnode_4,:)) = i_value 
bc_dist(edge_x(qnode_1,qnode_4,:)) = i_value  

!only possible edge connected with surface is qnode_2,qnode_3
! if any node is equal to or above surface leave out of distance calculation meaing taper is not applied
not_surface = .true. 
do i = 1,size(edge_x,3)
    if (amesh%pz(edge_x(qnode_2,qnode_3,i)) >= surface-0.001) then   ! the 0.001 is to account for possible rounding error 
        not_surface = .false.
        exit 
    endif 
enddo 

if (not_surface) then 
    bc_dist(edge_x(qnode_2,qnode_3,:)) = i_value !if not on surface include in border
endif 

!bc_dist(edge_x(qnode_2,qnode_3,:)) =  i_value   

fast = .true.
call onevsall2d(amesh,bc_dist,fast)
bc_dist = bc_dist-i_value  !remove large number to get actual distance


return 
end subroutine boundary_dist
!====================================================================
subroutine add_slip(slipping_cell,qmesh,amesh,edge_x,edge_z,tri2quad,bc_dist,taper_dist,slip)
type(mesh),intent(in) :: qmesh,amesh
integer,dimension(:,:,:),intent(in) :: edge_x,edge_z
real(pr), dimension(:),intent(in) :: bc_dist
integer, dimension(:),intent(in) :: tri2quad
integer,intent(in) :: slipping_cell
integer,dimension(3) :: query_cell
real(pr), allocatable, dimension(:),intent(out) :: slip 
real(pr) :: d1,d2,d3,mean_d
real(pr),intent(in) :: taper_dist
integer:: qnode_1,qnode_2,qnode_3,qnode_4
integer :: istart,iend,i,j,nz
logical :: inside
! real(pr),parameter :: kilo = 1000._pr 

qnode_1 = qmesh%cell(slipping_cell,1)
qnode_2 = qmesh%cell(slipping_cell,2)
qnode_3 = qmesh%cell(slipping_cell,3)
qnode_4 = qmesh%cell(slipping_cell,4)

!taper_dist = 25._pr*kilo

nz = size(edge_z,3)
istart = edge_z(qnode_1,qnode_2,1)
iend = edge_z(qnode_3,qnode_4,nz)
! write(*,*) 'node range....',istart,iend

allocate(slip(amesh%Ncells))
slip = 0._pr

!apply uniform slip 
! do i = 1,amesh%ncells
!     inside = .true. 
!     query_cell = tri2quad(amesh%cell(i,1:3)) 
!     if (check_cell_quad(slipping_cell,query_cell)) then   ! inside slipping cell 
!         slip(i) = 1._pr 
!     endif 
! enddo 

! apply tapered slip 
do i = 1,amesh%ncells
    inside = .true. 
    query_cell = tri2quad(amesh%cell(i,1:3)) 
    d1 = bc_dist(amesh%cell(i,1))
    d2 = bc_dist(amesh%cell(i,2))
    d3 = bc_dist(amesh%cell(i,3))
    mean_d = (d1+d2+d3)/3._pr

    if (check_cell_quad(slipping_cell,query_cell)) then  ! inside slipping cell 
        slip(i) = mean_d/taper_dist+0.5
        if (slip(i) > 1._pr) slip(i) = 1._pr 
    else ! outside slipping cell 
        slip(i) = -mean_d/taper_dist+0.5
        if (slip(i) < 0._pr) slip(i) = 0._pr 
    endif 

enddo 



return
end subroutine add_slip
! !====================================================================
function check_cell_quad(slipping_cell,query_cell)    result(match)
integer,intent(in) :: slipping_cell
integer,dimension(3) :: query_cell
integer,dimension(7,3) :: test
integer :: i 
logical :: match

match = .false.
test(1,:) = (/slipping_cell, slipping_cell, slipping_cell /)
test(2,:) = (/0, 0, slipping_cell /)
test(3,:) = (/0, slipping_cell, 0 /)
test(4,:) = (/slipping_cell, 0, 0 /)
test(5,:) = (/slipping_cell, slipping_cell, 0 /)
test(6,:) = (/slipping_cell, 0, slipping_cell /)
test(7,:) = (/0, slipping_cell, slipping_cell /)

do i = 1,7
    if (all(query_cell == test(i,:))) then 
        match = .true.
        exit 
    endif
enddo 

return
end function check_cell_quad
! !====================================================================

end module mod_slip
