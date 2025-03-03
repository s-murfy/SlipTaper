module mod_slip
use generic 
use LAT_mesh
implicit none 

contains

!====================================================================
subroutine boundary_distV2(qmesh,amesh,edge_x,edge_z,surface,bc_dist)
use distance
integer,allocatable,dimension(:,:,:),intent(in) :: edge_x,edge_z
type(mesh),intent(in) :: qmesh,amesh
real(pr),allocatable,dimension(:,:) :: bc_dist
real(pr),allocatable,dimension(:) :: quad_dist

real(pr),intent(in) :: surface
real(pr) :: i_value
integer :: qnode_1,qnode_2,qnode_3,qnode_4
logical :: fast ,not_surface
integer :: i,j


if (allocated(bc_dist)) deallocate(bc_dist)
allocate(bc_dist(qmesh%Nnodes,amesh%Nnodes))
allocate(quad_dist(amesh%Nnodes))
i_value = 10000000._pr  ! use a large number to generate for the generation of plane wavefront
fast = .true.

do i = 1,qmesh%Ncells
    quad_dist = infinity
    qnode_1 = qmesh%cell(i,1)
    qnode_2 = qmesh%cell(i,2)
    qnode_3 = qmesh%cell(i,3)
    qnode_4 = qmesh%cell(i,4)
!set sides of quadrangle to zero 
    quad_dist(edge_z(qnode_1,qnode_2,:)) = i_value
    quad_dist(edge_z(qnode_3,qnode_4,:)) = i_value
    quad_dist(edge_x(qnode_1,qnode_4,:)) = i_value
!only possible edge connected with surface is qnode_2,qnode_3
! if any node is equal to or above surface leave out of distance calculation meaing taper is not applied
    not_surface = .true. 
    do j = 1,size(edge_x,3)
        if (amesh%pz(edge_x(qnode_2,qnode_3,j)) >= surface-0.001) then   ! the 0.001 is to account for possible rounding error 
            not_surface = .false.
            exit 
        endif 
    enddo 
    if (not_surface) then 
        quad_dist(edge_x(qnode_2,qnode_3,:)) = i_value !if not on surface include in border
    endif 
    call onevsall2d(amesh,quad_dist,fast)
    quad_dist = quad_dist-i_value  !remove large number to get actual distance
    bc_dist(i,:) = quad_dist
enddo
return 
end subroutine boundary_distV2
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
i_value = 4000._pr
corner_dist(edge_x(qnode_1,qnode_4,1)) =  0._pr !i_value 
corner_dist(edge_x(qnode_1,qnode_4,nx)) = 0._pr  !i_value
corner_dist(edge_x(qnode_2,qnode_3,1)) =  0._pr !i_value 
corner_dist(edge_x(qnode_2,qnode_3,nx)) = 0._pr !i_value


! corner_dist(edge_x(qnode_1,qnode_4,1)) =  i_value 
! corner_dist(edge_x(qnode_1,qnode_4,nx)) = i_value
! corner_dist(edge_x(qnode_2,qnode_3,1)) =  i_value 
! corner_dist(edge_x(qnode_2,qnode_3,nx)) = i_value

! set the 3 points in the corner to zero
! corner_dist(edge_x(qnode_1,qnode_4,1:2)) = i_value 
! corner_dist(edge_x(qnode_1,qnode_4,nx-1:nx)) = i_value
! corner_dist(edge_x(qnode_2,qnode_3,1:2)) = i_value 
! corner_dist(edge_x(qnode_2,qnode_3,nx-1:nx)) = i_value 

! corner_dist(edge_z(qnode_3,qnode_4,1:2)) = i_value  
! corner_dist(edge_z(qnode_3,qnode_4,nz-1:nz)) = i_value 
! corner_dist(edge_z(qnode_1,qnode_2,1:2)) = i_value  
! corner_dist(edge_z(qnode_1,qnode_2,nz-1:nz)) = i_value 
fast =.true.
call onevsall2d(amesh,corner_dist,fast)
! corner_dist = corner_dist-i_value 



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
subroutine add_slipV2(slipping_cell,qmesh,amesh,edge_x,edge_z,bc_dist,taper_dist,slip)
type(mesh),intent(in) :: qmesh,amesh
integer,dimension(:,:,:),intent(in) :: edge_x,edge_z
real(pr), dimension(:,:),intent(in) :: bc_dist
integer,intent(in) :: slipping_cell
integer,dimension(3) :: query_cell
real(pr), allocatable, dimension(:),intent(out) :: slip 
real(pr) :: d1,d2,d3,d4,mean_d,edge_d,bary_x,bary_y
real(pr),intent(in) :: taper_dist
integer:: qnode_1,qnode_2,qnode_3,qnode_4
integer,dimension(3) :: query_nodes
integer :: i,j,nz,nx
logical :: inside,edge_exist,is_edge_node
real(pr),dimension(4) :: quad_x,quad_y
real(pr),allocatable,dimension(:,:) :: all_slip
real(pr) :: slip_edge
real(pr) :: const,total_slip,correction
! real(pr),parameter :: kilo = 1000._pr 
integer :: count 
real(pr) :: c1,c2,c3,c4


allocate(all_slip(qmesh%Ncells,amesh%Ncells))
allocate(slip(amesh%Ncells))
slip = 0._pr
all_slip = 0._pr

do j = 1,qmesh%Ncells
    qnode_1 = qmesh%cell(j,1)
    qnode_2 = qmesh%cell(j,2)
    qnode_3 = qmesh%cell(j,3)
    qnode_4 = qmesh%cell(j,4)

    quad_x = (/ qmesh%px(qnode_1),qmesh%px(qnode_2),qmesh%px(qnode_3),qmesh%px(qnode_4)  /)
    quad_y = (/ qmesh%py(qnode_1),qmesh%py(qnode_2),qmesh%py(qnode_3),qmesh%py(qnode_4) /)

    i = 1
    do while (i<amesh%ncells)
        bary_x = (amesh%px(amesh%cell(i,1))+amesh%px(amesh%cell(i,2))+amesh%px(amesh%cell(i,3))+amesh%px(amesh%cell(i+1,3)))/4._pr
        bary_y = (amesh%py(amesh%cell(i,1))+amesh%py(amesh%cell(i,2))+amesh%py(amesh%cell(i,3))+amesh%py(amesh%cell(i+1,3)))/4._pr

    ! inside = .true.
        inside = pt_inside_quad(bary_x,bary_y,quad_x,quad_y)   ! check if cell is inside the slipping quadrangle

        d1 = bc_dist(j,amesh%cell(i,1))
        d2 = bc_dist(j,amesh%cell(i,2))
        d3 = bc_dist(j,amesh%cell(i,3))
        d4 = bc_dist(j,amesh%cell(i+1,3))
        mean_d = (d1+d2+d3+d4)/4._pr
        ! mean_d = (d1+d2+d3)/3._pr
    const =0.5_pr
    if (inside) then  ! inside slipping cell 
        all_slip(j,i) = const+mean_d/taper_dist
    else ! outside slipping cell 
        all_slip(j,i) = const-mean_d/taper_dist  
    endif

    if (all_slip(j,i) > 1._pr) all_slip(j,i) = 1._pr 
    if (all_slip(j,i) < 0._pr) all_slip(j,i) = 0._pr 
    all_slip(j,i+1) = all_slip(j,i) 
    i = i+2
    enddo 
enddo 

! 
! qnode_1 = qmesh%cell(slipping_cell,1)
! qnode_2 = qmesh%cell(slipping_cell,2)
! qnode_3 = qmesh%cell(slipping_cell,3)
! qnode_4 = qmesh%cell(slipping_cell,4)
! nz = size(edge_z,3)
! nx = size(edge_x,3)
! c1 = edge_x(qnode_1,qnode_4,1) 
! c2 = edge_x(qnode_2,qnode_3,nx)
! c3 = edge_x(qnode_2,qnode_3,1)
! c4 = edge_x(qnode_1,qnode_4,nx)


i = 1
do while (i<amesh%ncells)

    query_nodes = amesh%cell(i,:)
    is_edge_node = check_slab_edge(query_nodes,slipping_cell,edge_x,edge_z,qmesh)! check if cell has a node on the edge 

    if (is_edge_node) then  ! at the edge of slab and we do not  apply a global balanace of slip
    ! write(*,*) i,'on an edge'  
        slip(i) = all_slip(slipping_cell,i)
    else 
        total_slip = sum(all_slip(:,i))
        correction = 1._pr/total_slip
        slip(i) = all_slip(slipping_cell,i)*correction
    endif 


    ! count =  0
    ! do j = 1,qmesh%Ncells
    !     if (all_slip(j,i) > epsilon(1._pr)) then
    !       count = count +1
    !     endif 
    ! enddo 

    ! if (count > 1) then 
    !     total_slip = sum(all_slip(:,i))
    !     ! correction = 1._pr
    !     ! if (total_slip > 1) correction = 1._pr/total_slip
    !     correction = 1._pr/total_slip
    !     slip(i) = all_slip(slipping_cell,i)*correction
    ! else    ! at the edge of slab and we do not need to apply a global balanace 
    !     slip(i) = all_slip(slipping_cell,i)
    ! endif 

    slip(i+1) = slip(i)
    i = i+2
enddo 

return
end subroutine add_slipV2
!====================================================================



!====================================================================
subroutine add_slip(slipping_cell,qmesh,amesh,edge_x,edge_z,tri2quad,bc_dist,edge_dist,corner_dist,edge_exist,taper_dist,slip)
type(mesh),intent(in) :: qmesh,amesh
integer,dimension(:,:,:),intent(in) :: edge_x,edge_z
real(pr), dimension(:),intent(in) :: bc_dist,edge_dist,corner_dist
integer, dimension(:),intent(in) :: tri2quad
integer,intent(in) :: slipping_cell
integer,dimension(3) :: query_cell
real(pr), allocatable, dimension(:),intent(out) :: slip 
real(pr) :: d1,d2,d3,d4,mean_d,edge_d,bary_x,bary_y
real(pr),intent(in) :: taper_dist
integer:: qnode_1,qnode_2,qnode_3,qnode_4
integer,dimension(3) :: query_nodes
integer :: istart,iend,i,j,nz,nx
logical :: inside,edge_exist,is_edge_node
real(pr),dimension(4) :: quad_x,quad_y
real(pr) :: slip_edge
real(pr) :: c1,c2,c3,c4,dx,const,alpha
! real(pr),parameter :: kilo = 1000._pr 

qnode_1 = qmesh%cell(slipping_cell,1)
qnode_2 = qmesh%cell(slipping_cell,2)
qnode_3 = qmesh%cell(slipping_cell,3)
qnode_4 = qmesh%cell(slipping_cell,4)

quad_x = (/ qmesh%px(qnode_1),qmesh%px(qnode_2),qmesh%px(qnode_3),qmesh%px(qnode_4)  /)
quad_y = (/ qmesh%py(qnode_1),qmesh%py(qnode_2),qmesh%py(qnode_3),qmesh%py(qnode_4) /)

! istart = edge_z(qnode_1,qnode_2,1)
! iend = edge_z(qnode_3,qnode_4,nz)
! write(*,*) 'node range....',istart,iend
nz = size(edge_z,3)
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
! dx = sum(abs(bc_dist(amesh%cell(i,1))-bc_dist(amesh%cell(i,2))),      &
!          abs(bc_dist(amesh%cell(i,1))-bc_dist(amesh%cell(i,3))),      &
!          abs(bc_dist(amesh%cell(i,2))-bc_dist(amesh%cell(i,3)))  )/3._pr
!     write(*,*) 'dx',dx !this give a rough estimate of cell length
!     alpha = dx/taper_dist  
! apply tapered slip 
! do i = 1,amesh%ncells
i = 1
do while (i<amesh%ncells)
    bary_x = (amesh%px(amesh%cell(i,1))+amesh%px(amesh%cell(i,2))+amesh%px(amesh%cell(i,3))+amesh%px(amesh%cell(i+1,3)))/4._pr
    bary_y = (amesh%py(amesh%cell(i,1))+amesh%py(amesh%cell(i,2))+amesh%py(amesh%cell(i,3))+amesh%py(amesh%cell(i+1,3)))/4._pr

! estimate the lenght of a cell 
    dx = max(abs(bc_dist(amesh%cell(i,1))-bc_dist(amesh%cell(i,2))),      &
         abs(bc_dist(amesh%cell(i,1))-bc_dist(amesh%cell(i,3))),      &
         abs(bc_dist(amesh%cell(i,2))-bc_dist(amesh%cell(i,3)))  )
    dx = 6000.0_pr
    ! write(*,*) 'dx',dx !this give a rough estimate of cell length
    alpha = dx/taper_dist  

    ! inside = .true.
    inside = pt_inside_quad(bary_x,bary_y,quad_x,quad_y)   ! check if cell is inside the slipping quadrangle

    query_cell = tri2quad(amesh%cell(i,1:3)) 
    d1 = bc_dist(amesh%cell(i,1))
    d2 = bc_dist(amesh%cell(i,2))
    d3 = bc_dist(amesh%cell(i,3))
    ! mean_d = (d1+d2+d3)/3._pr

    d4 = bc_dist(amesh%cell(i+1,3))
    mean_d = (d1+d2+d3+d4)/4._pr

    slip(i) = 0._pr
    slip_edge = 0._pr


    d1 = corner_dist(amesh%cell(i,1))
    d2 = corner_dist(amesh%cell(i,2))
    d3 = corner_dist(amesh%cell(i,3))
    ! edge_d = (d1+d2+d3)/3._pr

    d4 = bc_dist(amesh%cell(i+1,3))
    edge_d = (d1+d2+d3+d4)/4._pr

    query_nodes = amesh%cell(i,:)
    is_edge_node = check_edge_node(query_nodes,c1,c2,c3,c4) ! check if cell has a node on the edge 

    const = intercept_value(i,nx,nz,alpha,inside) ! check if cell is on the edge of the mesh
    const =0.5_pr
    if (inside) then  ! inside slipping cell 

    ! if (check_cell_quad(slipping_cell,query_cell)) then  ! inside slipping cell 
        ! if (is_edge_node) then 
        if (abs(const-0.5_pr) <= 0.0001) then 
            slip(i) = mean_d/taper_dist+const !0.4
        else
            slip(i) = edge_d/taper_dist+const !+0.25 !
        endif

        ! if (edge_d <= dx*3) then   ! apply a different taper condition for edges
        !     slip(i) = edge_d/(taper_dist*4)+intercept_value !+0.25 !
        ! else
        !     slip(i) = mean_d/taper_dist+intercept_value !0.4
        ! endif 

        if (edge_exist) then ! case when on the edge of larger mesh
            d1 = edge_dist(amesh%cell(i,1))
            d2 = edge_dist(amesh%cell(i,2))
            d3 = edge_dist(amesh%cell(i,3))
            ! mean_d = (d1+d2+d3)/3._pr
            d4 = edge_dist(amesh%cell(i+1,3))
            mean_d = (d1+d2+d3+d4)/4._pr

            slip_edge = mean_d/taper_dist+0.1
        endif 
        slip(i) = max(slip(i),slip_edge)
        
    else ! outside slipping cell 

        if (abs(const-0.5_pr) <= 0.0001) then 
            slip(i) = -mean_d/taper_dist+const !0.4     ! using side distance
        else
            slip(i) = -edge_d/taper_dist+const !+0.25 !   ! using corner distance
        endif


        ! if (is_edge_node) then      
        ! if (edge_d <= dx*2) then     ! apply a different taper condition for edges
        !     slip(i) = -edge_d/(taper_dist)+0.42 !25   !0.35
        ! else
        !      slip(i) = -mean_d/taper_dist+0.54 !0.4
        ! endif 

        ! if (i==401) write(*,*) 'slip value',i, slip(i), -edge_d/taper_dist+0.25,edge_d,edge_exist
        ! if (i==997) write(*,*) 'slip value',i, slip(i), -edge_d/taper_dist+0.25,edge_d,edge_exist

        ! if (i==19) write(*,*) 'slip value',i, slip(i),-edge_d/taper_dist+0.25,edge_d,edge_exist
        ! if (i==981) write(*,*) 'slip value',i, slip(i), -edge_d/taper_dist+0.25,edge_d,edge_exist


        ! slip(i) = 0._pr
        if (edge_exist) then ! case when on the edge of larger mesh
            d1 = edge_dist(amesh%cell(i,1))
            d2 = edge_dist(amesh%cell(i,2))
            d3 = edge_dist(amesh%cell(i,3))
            ! d3 = edge_dist(amesh%cell(i,3))
            d4 = edge_dist(amesh%cell(i+1,3))
            mean_d = (d1+d2+d3+d4)/4._pr

            ! mean_d = (d1+d2+d3)/3._pr
            slip_edge = -mean_d/taper_dist+0.1
        endif 
        slip(i) = max(slip(i),slip_edge)

    endif 

        ! if (i==401) write(*,*) 'slip value',i, slip(i), -edge_d/taper_dist+0.25,edge_d,edge_exist
        ! if (i==997) write(*,*) 'slip value',i, slip(i), -edge_d/taper_dist+0.25,edge_d,edge_exist

        ! if (i==19) write(*,*) 'slip value',i, slip(i),-edge_d/taper_dist+0.25,edge_d,edge_exist
        ! if (i==981) write(*,*) 'slip value',i, slip(i), -edge_d/taper_dist+0.25,edge_d,edge_exist

    ! if (const == 10000) slip(i) = 0.25_pr

    if (slip(i) > 1._pr) slip(i) = 1._pr 
    if (slip(i) < 0._pr) slip(i) = 0._pr 

    slip(i+1) = slip(i) 
    i = i+2
enddo 

! correct for slip in corners 
call smoothen_slip(qmesh,slipping_cell,slip,nx,nz)

return
end subroutine add_slip

!====================================================================
function find_adjacent_cells(qmesh,slipping_cell) result(list)
type(mesh),intent(in) :: qmesh
integer,intent(in) :: slipping_cell
integer :: qnode_1,qnode_2,qnode_3,qnode_4,i,j,inx
integer,dimension(4) :: nodes
integer,dimension(4,4) :: list 

qnode_1 = qmesh%cell(slipping_cell,1)
qnode_2 = qmesh%cell(slipping_cell,2)
qnode_3 = qmesh%cell(slipping_cell,3)
qnode_4 = qmesh%cell(slipping_cell,4)
nodes = (/qnode_1,qnode_2,qnode_3,qnode_4/)
list = 0 

do j = 1,4
    inx = 0
    do i = 1,qmesh%Ncells
        if ( any( qmesh%cell(i,:) == nodes(j) )) then 
            inx = inx+1
            list(j,inx) = i
        endif 
    enddo
enddo 

call sort_rows(list)  ! order the cells in the list starting from the smallest

return 
end function find_adjacent_cells
!====================================================================
! Subroutine to sort each row of a 4x4 array
subroutine sort_rows(a)
integer, intent(inout) :: a(4,4)
integer :: i, j, k, temp

    do i = 1, 4
        do j = 1, 3
            do k = 1, 4 - j
                if (a(i, k) > a(i, k+1)) then
                    temp = a(i, k)
                    a(i, k) = a(i, k+1)
                    a(i, k+1) = temp
                end if
            end do
        end do
    end do
return 
end subroutine sort_rows
!====================================================================
subroutine smoothen_slip(qmesh,slipping_cell,slip,nx,nz)
type(mesh),intent(in) :: qmesh
integer,intent(in) :: slipping_cell
real(pr),intent(inout) :: slip(:)
integer,intent(in) :: nx,nz
integer :: i,j,k,ii,nx_cells,nz_cells,i0,ncells_per_quad
integer :: i1,i2,i3,i4
integer,dimension(4,4) :: cell_list
integer,dimension(4) :: indices,list 
real(pr) :: total_slip

nx_cells = 2*(nx-1)
nz_cells = (nz-1)
! i0 = int(ii/nx_cells/nz_cells)
ncells_per_quad = nx_cells*nz_cells
! write(*,*)'total no. of cells, ncells per quad', size(slip),ncells_per_quad 

cell_list = find_adjacent_cells(qmesh,slipping_cell)
cell_list = cell_list-1
do i = 1, 4
    ! print '(4I5)', (cell_list(i, j), j=1, 4)
! condition to correct the cells right in the corner 
! conditions for point in the corner cell (1,1),(nx,nz),(1,nx),(1,nz)
    if (cell_list(i, 1) < 0) cycle
    do k = 1,4  ! loop over corner cell combinations 
        list = cell_list(i,:)
        indices = corner_indices(k, nx_cells, nz_cells,list)
        i1 = indices(1)
        i2 = indices(2) 
        i3 = indices(3)
        i4 = indices(4)

        total_slip = slip(i1)+slip(i2)+slip(i3)+slip(i4)

        slip(i1) = slip(i1)/total_slip   
        slip(i2) = slip(i2)/total_slip 
        slip(i3) = slip(i3)/total_slip   
        slip(i4) = slip(i4)/total_slip    
        !make sure you correct both cells! 
        slip(i1+1) = slip(i1+1)/total_slip   
        slip(i2+1) = slip(i2+1)/total_slip 
        slip(i3+1) = slip(i3+1)/total_slip   
        slip(i4+1) = slip(i4+1)/total_slip    
    enddo 
enddo 
return 

end subroutine smoothen_slip
!====================================================================
function corner_indices(i, nx_cells, nz_cells,list) result(indices)
integer,intent(in) :: i,nx_cells,nz_cells
integer,dimension(4),intent(in) :: list
integer :: i1,i2,i3,i4
integer :: ncells_per_quad
integer,dimension(4) :: indices

ncells_per_quad = nx_cells*nz_cells

if (i == 1) then  ! this is the index for the cells in the corner
!  if ((i==nx_cells*nz_cells-1).or.(i==1).or.(i==nx_cells-1).or.(i==nx_cells*(nz_cells-1)+1)) then 
    i1  = nx_cells*nz_cells-1 +list(3)*ncells_per_quad    ! in the quadrangle below the first 
    i2 = 1+list(2)*ncells_per_quad                            !in the second quadrangle (right, slipping cell)
    i3 = nx_cells-1+list(1)*ncells_per_quad               ! in the first quadrangle (left) 
    i4 = nx_cells*(nz_cells-1)+1+list(4)*ncells_per_quad  ! in the quadrangle below second

elseif(i==2) then ! this the index for the cells diagonal from the corner 
!  ((i==nx_cells+3).or.(i==2*nx_cells-3).or.                                &
!         (i==(nz_cells-1)*nx_cells-3).or.(i==nx_cells*(nz_cells-2)+3)) then
    i1 = (nz_cells-1)*nx_cells-3+list(3)*ncells_per_quad  
    i2 = nx_cells+3+list(2)*ncells_per_quad 
    i3 = 2*nx_cells-3+list(1)*ncells_per_quad    
    i4 = nx_cells*(nz_cells-2)+3+list(4)*ncells_per_quad
elseif(i==3) then ! this the index for the cells with sides on the boundary

    i1 = nx_cells*nz_cells-3+list(3)*ncells_per_quad
    i2 = 3+list(2)*ncells_per_quad
    i3 = nx_cells-3+list(1)*ncells_per_quad
    i4 = nx_cells*(nz_cells-1)+3+list(4)*ncells_per_quad

elseif(i==4) then ! this the index for the cells that touch other boundary

    i1  = nx_cells*(nz_cells-1)-1 +list(3)*ncells_per_quad    ! in the quadrangle below the first 
    i2 = nx_cells+1+list(2)*ncells_per_quad                            !in the second quadrangle (right, slipping cell)
    i3 = 2*nx_cells-1+list(1)*ncells_per_quad               ! in the first quadrangle (left) 
    i4 = nx_cells*(nz_cells-2)+1+list(4)*ncells_per_quad  ! in the quadrangle below second

! elseif (   (i==3).or.(i==nx_cells-3).or.(i==nx_cells+1).or.(i==2*nx_cells-1).or.     &
!            (i==nx_cells*nz_cells-3).or.(i==nx_cells*(nz_cells-1)+3).or.              &
!            (i==nx_cells*(nz_cells-1)-1).or.(i==nx_cells*(nz_cells-2)+1)  ) then       

endif 

indices = (/i1,i2,i3,i4/)
return 
end function corner_indices
!====================================================================

function intercept_value(ii,nx,nz,alpha,inside) result(constant)
integer,intent(in) :: ii,nx,nz
real(pr),intent(in) :: alpha
logical,intent(in) :: inside
real(pr) :: constant,c_inside
integer :: nx_cells,nz_cells
integer :: i,i0

!i = ii+1  ! add one more to get cell that is right on the boundary 
nx_cells = 2*(nx-1)
nz_cells = (nz-1)
i0 = int(ii/nx_cells/nz_cells)
i = ii - i0*nx_cells*nz_cells

! conditions for point in the corner cell (1,1),(nx,nz),(1,nx),(1,nz)
 if ((i==nx_cells*nz_cells-1).or.(i==1).or.(i==nx_cells-1).or.(i==nx_cells*(nz_cells-1)+1)) then 
    ! constant = 0.25_pr*(1._pr-3._pr/sqrt(2._pr)*alpha) ! cell is in the corner
    c_inside = 0.25_pr

    if (inside) then 
        constant = c_inside
    else
        constant = 1._pr/3._pr*(1._pr-c_inside+sqrt(2._pr)*alpha) ! cell is in the corner
    endif

    ! write(*,*) 'corner cell intercep',i,i0,constant,alpha

!cell is one cell away from either boundary
elseif ((i==nx_cells+3).or.(i==2*nx_cells-3).or.                                &
        (i==(nz_cells-1)*nx_cells-3).or.(i==nx_cells*(nz_cells-2)+3)) then

    ! constant = 0.25_pr*(1._pr-9./sqrt(2._pr)*alpha) ! cell is one cell away from the boundary
    ! constant = 0.25_pr*(1._pr+3.*sqrt(2._pr)*alpha) ! cell is one cell away from the boundary

    c_inside = 0.2_pr
    if (inside) then 
        constant = c_inside
    else
        constant = 1._pr/3._pr*(1._pr-c_inside+3.*sqrt(2._pr)*alpha) ! cell is one cell away from the boundary
    endif
    ! write(*,*) 'corner cell intercep',i,i0,constant,alpha
    ! write(*,*) 'away from edges cell intercep',i,constant,alpha


! cell has one size on a boundary 
elseif (   (i==3).or.(i==nx_cells-3).or.(i==nx_cells+1).or.(i==2*nx_cells-1).or.     &
           (i==nx_cells*nz_cells-3).or.(i==nx_cells*(nz_cells-1)+3).or.              &
           (i==nx_cells*(nz_cells-1)-1).or.(i==nx_cells*(nz_cells-2)+1)  ) then       
   
    ! constant = 0.25_pr*(1._pr-sqrt(10._pr)*3._pr/2._pr*alpha) ! cell is one cell away from the boundary
    ! constant = 0.25_pr*(1._pr+sqrt(10._pr)*alpha) ! cell is one cell away from the boundary

    ! constant = 10000.0

    c_inside = 0.05_pr
    if (inside) then 
        constant = c_inside 
    else
        constant = 1._pr/3._pr*(1._pr-c_inside+sqrt(10._pr)*alpha) ! cell is one cell away from the boundary
    endif

    ! write(*,*) 'touching one boundary intercep',i,constant,alpha
    ! write(*,*) 'corner cell intercep',i,i0,constant,alpha


else 
    constant = 0.5_pr
endif

return 
end function intercept_value
!====================================================================
function pt_inside_quad(x,y,quad_x,quad_y) result(inside)
real(pr),intent(in) :: x,y
real(pr),dimension(4),intent(in) :: quad_x,quad_y
logical :: inside
real(pr) :: x1,x2,y1,y2,x_intersect
integer :: i,j

inside = .false.

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
function check_slab_edge(nodes,slipping_cell,edge_x,edge_z,qmesh) result(is_edge_node)
type(mesh),intent(in) :: qmesh
integer,dimension(3),intent(in) :: nodes
integer,dimension(:,:,:),intent(in) :: edge_x,edge_z
integer,intent(in) :: slipping_cell
integer :: q1,q2,q3,q4
integer :: i,j,count
logical :: is_edge_node

q1 = qmesh%cell(slipping_cell,1)
q2 = qmesh%cell(slipping_cell,2)
q3 = qmesh%cell(slipping_cell,3)
q4 = qmesh%cell(slipping_cell,4)

is_edge_node = .false.


! check one vertical side of quadrangle 
count = 0
do i = 1,qmesh%Ncells
    if (i == slipping_cell) cycle
    do j = 1,4 
        if (q1==qmesh%cell(i,j) .or. q2==qmesh%cell(i,j)) then
            count = count +1
        endif 
    enddo 
enddo 
if (count < 2) then 
    !write(*,*) 'on right slab edge'
    do j = 1,3
        if (any(nodes(j)==edge_z(q1,q2,:))) then 
                is_edge_node = .true.
                return
        endif 
    enddo 
endif 


! check other vertical side of quadrangle 

count = 0
do i = 1,qmesh%Ncells
    if (i == slipping_cell) cycle
    do j = 1,4 
        if (q3 == qmesh%cell(i,j) .or. q4 == qmesh%cell(i,j)) then
            count = count +1
        endif 
    enddo 
enddo 
if (count < 2) then 
    !write(*,*) 'on left slab edge'
        do j = 1,3
            if (any(nodes(j)==edge_z(q3,q4,:))) then 
                is_edge_node = .true.
                return 
            endif 
        enddo 
     
endif 


! check bottom side of quadrangle 
count = 0
do i = 1,qmesh%Ncells
    if (i == slipping_cell) cycle
    do j = 1,4 
        if (q1 == qmesh%cell(i,j) .or. q4 == qmesh%cell(i,j)) then
            count = count +1
        endif 
    enddo 
enddo 

if (count < 2) then 
        write(*,*) 'on bottom slab edge'
        do j = 1,3
            if (any(nodes(j)==edge_x(q1,q2,:))) then 
                is_edge_node = .true.
                return
            endif 
        enddo 
endif 

return
end function check_slab_edge
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
