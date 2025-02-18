module mod_slip
use generic 
use LAT_mesh
implicit none 

contains

!====================================================================
subroutine boundary_dist(slipping_cell,qmesh,amesh,edge_x,edge_z,surface,bc_dist,edge_dist,corner_dist,edge_exist)
use distance
integer,allocatable,dimension(:,:,:),intent(in) :: edge_x,edge_z
type(mesh),intent(in) :: qmesh,amesh
real(pr),allocatable,dimension(:) :: bc_dist,edge_dist,corner_dist
real(pr),intent(in) :: surface
real(pr) :: i_value
integer,intent(in) :: slipping_cell
integer :: qnode_1,qnode_2,qnode_3,qnode_4
logical :: fast ,not_surface
integer :: i,nx,nz
logical :: edge_12,edge_34,edge_14,edge_23,edge_exist


qnode_1 = qmesh%cell(slipping_cell,1)
qnode_2 = qmesh%cell(slipping_cell,2)
qnode_3 = qmesh%cell(slipping_cell,3)
qnode_4 = qmesh%cell(slipping_cell,4)



if (allocated(bc_dist)) deallocate(bc_dist)
if (allocated(edge_dist)) deallocate(edge_dist)

allocate(bc_dist(amesh%Nnodes))
bc_dist = infinity
allocate(edge_dist(amesh%Nnodes))
edge_dist = infinity

i_value = 10000000._pr  ! use a large number to generate for the generation of plane wavefront
!set sides of quadrangle to zero 
! bc_dist(edge_z(qnode_1,qnode_2,:)) = i_value
! bc_dist(edge_z(qnode_3,qnode_4,:)) = i_value 
! bc_dist(edge_x(qnode_1,qnode_4,:)) = i_value  

! check if any of the edges are interior edges 
edge_23 = interior_edge(qmesh,qnode_2,qnode_3,slipping_cell)

edge_12 = interior_edge(qmesh,qnode_1,qnode_2,slipping_cell)
edge_34 = interior_edge(qmesh,qnode_3,qnode_4,slipping_cell)
edge_14 = interior_edge(qmesh,qnode_1,qnode_4,slipping_cell)

! if (.not.edge_12) write(*,*) 'edge 12 is on edge of mesh'
! if (.not.edge_34) write(*,*) 'edge 34 is on edge of mesh'
! if (.not.edge_14) write(*,*) 'edge 14 is on edge of mesh'

if((.not.edge_12).or.(.not.edge_34).or.(.not.edge_14)) then 
    ! write(*,*) 'we have an external edge !'
    edge_exist = .true.
endif 

if (edge_12) then 
    bc_dist(edge_z(qnode_1,qnode_2,:)) = i_value
else    ! edge of slab , using different taper condition 
    edge_dist(edge_z(qnode_1,qnode_2,:)) = i_value
    do i =1,qmesh%Nnodes ! also need to set sides on adjoining quadrangle cells also 
        edge_dist(edge_z(qnode_1,i,:)) = i_value
        edge_dist(edge_z(i,qnode_2,:)) = i_value
    enddo 
endif 

if (edge_34) then 
    bc_dist(edge_z(qnode_3,qnode_4,:)) = i_value
else     ! edge of slab , using different taper condition 
    edge_dist(edge_z(qnode_3,qnode_4,:)) = i_value
    do i =1,qmesh%Nnodes ! also need to set sides on adjoining quadrangle cells also 
        edge_dist(edge_z(qnode_3,i,:)) = i_value
        edge_dist(edge_z(i,qnode_4,:)) = i_value
    enddo
endif

if (edge_14) then 
    bc_dist(edge_x(qnode_1,qnode_4,:)) = i_value
else 
    edge_dist(edge_x(qnode_1,qnode_4,:)) = i_value
    do i =1,qmesh%Nnodes ! also need to set sides on adjoining quadrangle cells also 
        edge_dist(edge_x(qnode_1,i,:)) = i_value
        edge_dist(edge_x(i,qnode_4,:)) = i_value
    enddo
endif

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

    if(edge_23) then 
        bc_dist(edge_x(qnode_2,qnode_3,:)) = i_value !if not on surface include in border
    else
        edge_dist(edge_x(qnode_2,qnode_3,:)) = i_value
        do i =1,qmesh%Nnodes ! also need to set sides on adjoining quadrangle cells also 
            edge_dist(edge_x(i,qnode_3,:)) = i_value
            edge_dist(edge_x(qnode_2,i,:)) = i_value
        enddo
    endif 
endif 

!bc_dist(edge_x(qnode_2,qnode_3,:)) =  i_value   

fast = .true.
call onevsall2d(amesh,bc_dist,fast)
bc_dist = bc_dist-i_value  !remove large number to get actual distance

if (edge_exist) then 
    ! write(*,*) 'edge exist and distance from boundary being calculated'
    call onevsall2d(amesh,edge_dist,fast)
    edge_dist = edge_dist-i_value 
endif 

allocate(corner_dist(amesh%Nnodes))
corner_dist = infinity
nx = size(edge_x,3)
nz = size(edge_z,3)


! set just the point in the corner to zero 
! corner_dist(edge_x(qnode_1,qnode_4,1)) = i_value 
! corner_dist(edge_x(qnode_1,qnode_4,nx)) = i_value
! corner_dist(edge_x(qnode_2,qnode_3,1)) = i_value 
! corner_dist(edge_x(qnode_2,qnode_3,nx)) = i_value


! set the 3 points in the corner to zero
corner_dist(edge_x(qnode_1,qnode_4,1:2)) = i_value 
corner_dist(edge_x(qnode_1,qnode_4,nx-1:nx)) = i_value
corner_dist(edge_x(qnode_2,qnode_3,1:2)) = i_value 
corner_dist(edge_x(qnode_2,qnode_3,nx-1:nx)) = i_value 

corner_dist(edge_z(qnode_3,qnode_4,1:2)) = i_value  !0._pr
corner_dist(edge_z(qnode_3,qnode_4,nz-1:nz)) = i_value !0._pr
corner_dist(edge_z(qnode_1,qnode_2,1:2)) = i_value  !0._pr
corner_dist(edge_z(qnode_1,qnode_2,nz-1:nz)) = i_value !0._pr

call onevsall2d(amesh,corner_dist,fast)
corner_dist = corner_dist-i_value 



return 
end subroutine boundary_dist
!====================================================================
function interior_edge(qmesh,n1,n2,current_cell) result(is_interior_edge)
type(mesh),intent(in) :: qmesh
integer,intent(in) :: n1,n2,current_cell
integer :: i
logical :: is_interior_edge

is_interior_edge = .false.

do i = 1,qmesh%Ncells
    if (i == current_cell) cycle

    if (any(qmesh%cell(i,:)==n1).and.any(qmesh%cell(i,:)==n2)) then 
        is_interior_edge = .true.
        exit 
    endif
enddo 

return
end function interior_edge
!====================================================================
subroutine add_slip(slipping_cell,qmesh,amesh,edge_x,edge_z,tri2quad,bc_dist,edge_dist,corner_dist,edge_exist,taper_dist,slip)
type(mesh),intent(in) :: qmesh,amesh
integer,dimension(:,:,:),intent(in) :: edge_x,edge_z
real(pr), dimension(:),intent(in) :: bc_dist,edge_dist,corner_dist
integer, dimension(:),intent(in) :: tri2quad
integer,intent(in) :: slipping_cell
integer,dimension(3) :: query_cell
real(pr), allocatable, dimension(:),intent(out) :: slip 
real(pr) :: d1,d2,d3,mean_d,edge_d,bary_x,bary_y
real(pr),intent(in) :: taper_dist
integer:: qnode_1,qnode_2,qnode_3,qnode_4
integer,dimension(3) :: query_nodes
integer :: istart,iend,i,j,nz,nx
logical :: inside,edge_exist,is_edge_node
real(pr),dimension(4) :: quad_x,quad_y
real(pr) :: slip_edge
real(pr) :: c1,c2,c3,c4,dx
! real(pr),parameter :: kilo = 1000._pr 

qnode_1 = qmesh%cell(slipping_cell,1)
qnode_2 = qmesh%cell(slipping_cell,2)
qnode_3 = qmesh%cell(slipping_cell,3)
qnode_4 = qmesh%cell(slipping_cell,4)

quad_x = (/ qmesh%px(qnode_1),qmesh%px(qnode_2),qmesh%px(qnode_3),qmesh%px(qnode_4)  /)
quad_y = (/ qmesh%py(qnode_1),qmesh%py(qnode_2),qmesh%py(qnode_3),qmesh%py(qnode_4) /)

! nz = size(edge_z,3)
! istart = edge_z(qnode_1,qnode_2,1)
! iend = edge_z(qnode_3,qnode_4,nz)
! write(*,*) 'node range....',istart,iend
nx = size(edge_x,3)
c1 = edge_x(qnode_1,qnode_4,1) 
c2 = edge_x(qnode_2,qnode_3,nx)
c3 = edge_x(qnode_2,qnode_3,1)
c4 = edge_x(qnode_1,qnode_4,nx)
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
! estimate the lenght of a cell 
dx = max(abs(bc_dist(amesh%cell(1,1))-bc_dist(amesh%cell(1,2))),      &
         abs(bc_dist(amesh%cell(1,1))-bc_dist(amesh%cell(1,3))),      &
         abs(bc_dist(amesh%cell(1,2))-bc_dist(amesh%cell(1,3)))  )

! apply tapered slip 
do i = 1,amesh%ncells

    bary_x = (amesh%px(amesh%cell(i,1))+amesh%px(amesh%cell(i,2))+amesh%px(amesh%cell(i,3)))/3._pr
    bary_y = (amesh%py(amesh%cell(i,1))+amesh%py(amesh%cell(i,2))+amesh%py(amesh%cell(i,3)))/3._pr

    ! inside = .true.
    inside = pt_inside_quad(bary_x,bary_y,quad_x,quad_y)   ! check if cell is inside the slipping quadrangle

    query_cell = tri2quad(amesh%cell(i,1:3)) 
    d1 = bc_dist(amesh%cell(i,1))
    d2 = bc_dist(amesh%cell(i,2))
    d3 = bc_dist(amesh%cell(i,3))
    mean_d = (d1+d2+d3)/3._pr
    slip_edge = 0._pr


    d1 = corner_dist(amesh%cell(i,1))
    d2 = corner_dist(amesh%cell(i,2))
    d3 = corner_dist(amesh%cell(i,3))
    edge_d = (d1+d2+d3)/3._pr

    
    query_nodes = amesh%cell(i,:)
    is_edge_node = check_edge_node(query_nodes,c1,c2,c3,c4) ! check if cell has a node on the edge 

    if (inside) then  ! inside slipping cell 

    ! if (check_cell_quad(slipping_cell,query_cell)) then  ! inside slipping cell 
        ! if (is_edge_node) then 
        if (edge_d <= dx) then   ! apply a different taper condition for edges
        ! if (edge_d <= taper_dist*0.25) then
            slip(i) = edge_d/taper_dist+0.25 !
        else
            slip(i) = mean_d/taper_dist+0.5
        endif 

        if (edge_exist) then ! case when on the edge of larger mesh
            d1 = edge_dist(amesh%cell(i,1))
            d2 = edge_dist(amesh%cell(i,2))
            d3 = edge_dist(amesh%cell(i,3))
            mean_d = (d1+d2+d3)/3._pr
            slip_edge = edge_d/taper_dist+0.1
        endif 
        slip(i) = max(slip(i),slip_edge)

    else ! outside slipping cell 

        ! if (is_edge_node) then      
        if (edge_d <= dx) then     ! apply a different taper condition for edges
        ! if (edge_d <= taper_dist*0.25) then
            slip(i) = -edge_d/taper_dist+0.25   !0.35
        else
             slip(i) = -mean_d/taper_dist+0.5
        endif 
        ! slip(i) = 0._pr
        if (edge_exist) then ! case when on the edge of larger mesh
            d1 = edge_dist(amesh%cell(i,1))
            d2 = edge_dist(amesh%cell(i,2))
            d3 = edge_dist(amesh%cell(i,3))
            mean_d = (d1+d2+d3)/3._pr
            slip_edge = -edge_d/taper_dist+0.1
        endif 
        slip(i) = max(slip(i),slip_edge)

    endif 

    if (slip(i) > 1._pr) slip(i) = 1._pr 
    if (slip(i) < 0._pr) slip(i) = 0._pr 

enddo 



return
end subroutine add_slip
!====================================================================
function pt_inside_quad(x,y,quad_x,quad_y) result(inside)
real(pr),intent(in) :: x,y
real(pr),dimension(4),intent(in) :: quad_x,quad_y
logical :: inside
real(pr) :: x1,x2,y1,y2,x_intersect
integer :: i,j

inside = .false.

! if ((x+value >= minval(bc_x)).and.(x+value <= maxval(bc_x)).and.       &
!     (y+value >= minval(bc_y)).and.(y+value <= maxval(bc_y))) then 
!     inside = .true.
! endif 



j = 4  ! Last vertex index (since Fortran arrays start at 1)

do i = 1, 4
    x1 = quad_x(i)
    y1 = quad_y(i)
    x2 = quad_x(j)
    y2 = quad_y(j)

    ! Check if the point is between the y bounds of the edge
    if ((y1 > y) .neqv. (y2 > y)) then
    ! Compute intersection of polygon edge with horizontal ray
        x_intersect = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
        if (x < x_intersect) then
            inside = .not. inside
        end if
    end if
    j = i  ! Move to next edge
end do



return
end function pt_inside_quad
!====================================================================
function check_edge_node(nodes,c1,c2,c3,c4) result(is_edge_node)
integer,dimension(3),intent(in) :: nodes
real(pr),intent(in) :: c1,c2,c3,c4
logical :: is_edge_node

is_edge_node = .false.
if (any(nodes == c1).or.any(nodes == c2).or.any(nodes == c3).or.any(nodes == c4)) then 
    is_edge_node = .true.
endif

return
end function check_edge_node
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
