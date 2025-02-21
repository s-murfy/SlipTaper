module LAT_mesh
use generic
implicit none

type mesh
  integer(pin) :: Nnodes,Ncells! Nnodes : no of nodes ; Ncells : number of cells 
  real(pr), dimension(:), allocatable :: px,py,pz  !          Px,Py,Pz : x,y,z location of nodes (Nnodes)
  integer(pin), dimension(:,:), allocatable :: cell !             cell : relationship between cell and nodes (Ncells,3)
                                                    !                    The row gives the cell and the three columns refer 
                                                    !                    to node indexes
  integer(pin),allocatable,dimension(:,:) ::  EToE  !             EToE : element to element connectivity     (Ncells,3)
                                                    !                    This array gives the cell indexes of cells with conjoining faces  
                                                    !                    Row is the query cell, and the three columns are 
                                                    !                    the indexes to the cells that share a face with the query cell
                                                    !                    If column index = row there is no adjointing cell (i.e. EToE(i,3) = i)
  integer(pin), allocatable, dimension(:,:,:) ::   FVsNodes ! FVsNodes : array linking faces (pairs of nodes) with cells (Nnodes,Nodes,2)
                                                    !                    Sparse matrix where the first two indexes give the node index corresponding
                                                    !                    to one face between cells. The two values provided in the 3rd index are the 
                                                    !                    cell indexes corresponding to the face.
                                                    !                    Null values are set to 0 
end type mesh     


interface dumpnodeattributevtk
    module procedure idumpnodeattributevtk
    module procedure rdumpnodeattributevtk
end interface 


contains
!###############################################################################
subroutine find_slipping_cells(slipping_cell,fault_mesh,quad_lon,quad_lat,quad_z,nquad,output_level)
type(mesh),intent(in) :: fault_mesh
integer,intent(inout) :: slipping_cell
real(pr),allocatable,dimension(:,:) :: quad_lon,quad_lat,quad_z
integer,intent(out) :: nquad  
integer,dimension(4) :: nodes
integer,dimension(9) ::  choosen_cells
integer :: i,j,new_slipping_cell
integer,intent(in) :: output_level
nquad = 0   
choosen_cells = 0 
!choosen_cells(1) = slipping_cell
nodes = fault_mesh%cell(slipping_cell,:)

full_fault: do i = 1,fault_mesh%Ncells
      cell_nodes: do j = 1,4
            if (any(fault_mesh%cell(i,j) == nodes)) then 
               nquad = nquad+1
               choosen_cells(nquad) = i
               if (i == slipping_cell) new_slipping_cell = nquad

               exit cell_nodes
            endif
      enddo cell_nodes
enddo full_fault

slipping_cell = new_slipping_cell

if(output_level== 2) write(*,*) 'choosen quadrangles ................', choosen_cells
allocate(quad_lon(nquad,4),quad_lat(nquad,4),quad_z(nquad,4))

do i = 1,nquad
   quad_lon(i,:) = fault_mesh%px(fault_mesh%cell(choosen_cells(i),:))
   quad_lat(i,:) = fault_mesh%py(fault_mesh%cell(choosen_cells(i),:))
   quad_z(i,:) = fault_mesh%pz(fault_mesh%cell(choosen_cells(i),:))
enddo

return 
end subroutine find_slipping_cells
!###############################################################################
pure function check_anticlockwise(x,y,n) result(anticlockwise)
real(pr),dimension(:),intent(in) :: x,y
integer,intent(in) :: n
integer :: anticlockwise
real(pr) :: sum 
integer :: i

sum = 0._pr
do i = 1,n-1
     sum = sum + y(i+1)*x(i)-x(i+1)*y(i)
enddo 
sum = sum + y(1)*x(n)-x(1)*y(n)
anticlockwise = 0
if (sum > 0) then ! anti_clockwise
   anticlockwise = 1
elseif (sum < 0) then   !clockwise
   anticlockwise = -1
endif 


return
end function check_anticlockwise
!###############################################################################
subroutine create_tri_mesh(qmesh,nx,nz,amesh,edge_x,edge_z,tri2quad,output_level)
type(mesh),intent(in) :: qmesh 
integer,intent(in) :: nx,nz
type(mesh),intent(out) :: amesh
integer,allocatable,dimension(:,:,:),intent(out) :: edge_x,edge_z
integer, allocatable,dimension(:) :: nodes,old_to_new
integer,allocatable,dimension(:) :: tmp_tri2quad,tri2quad
logical,allocatable,dimension(:,:,:) :: edge_to_do_z,edge_to_do_x
integer :: nnodes
integer :: i,j,k,icount,c_count
integer :: qnode_1,qnode_2,qnode_3,qnode_4
integer :: ex_up_count,ex_down_count,ez_left_count,ez_right_count 
integer :: existing_node1,existing_node2,existing_node3,existing_node4
logical :: bc,cell_updated 
integer,dimension(2,3) :: cell 
real(pr),dimension(4) :: z,x,y
real(pr) :: dw,dl,xs,ys,zs,xe,ye,ze
real(pr) :: num ! used for random seed 
type(mesh) :: tmp
integer :: id_start, id_end
integer :: new_index,old_index 
integer :: node_exists
integer,intent(in) :: output_level
!
! Scheme for discretising quadrangles into triangles
!
!  ^ increaseing depth
!  |
!     quad_node_1.--------------------------. quad_node_4
!                |                          |
!                |      i.--------.i+1      |
!                |       |        |         |
!                |       |        |         |
!                |   i+nx.--------.i+nx+1   |
!                |                          |
!     quad_node_2.--------------------------. quad_node_3      


! ----------------surface ----------------------
!     quad_node_2.--------------------------. quad_node_3
!                |                          |
!                |   i+nx.--------.i+nx+1   |
!                |       |        |         |
!                |       |        |         |
!                |      i.--------.i+1      |
!                |                          |
!     quad_node_1.--------------------------. quad_node_4      


! this array defines nodes on the edge of the quadranles 
allocate(edge_x(qmesh%Nnodes,qmesh%Nnodes,nx),edge_z(qmesh%Nnodes,qmesh%Nnodes,nz))
edge_x = 0
edge_z = 0

tmp%Nnodes = nx*nz*qmesh%Ncells
allocate(tmp%px(tmp%Nnodes),tmp%py(tmp%Nnodes),tmp%pz(tmp%Nnodes))
tmp%px = 0.d0
tmp%py = 0.d0
tmp%pz = 0.d0

allocate(tmp_tri2quad(tmp%Nnodes))
tmp_tri2quad = 0

! discretise large quadrangles into smaller triangles: generate nodes 
icount = 0
large_cells: do k = 1,qmesh%Ncells 
      x = qmesh%px(qmesh%cell(k,:))
      y = qmesh%py(qmesh%cell(k,:))
      z = qmesh%pz(qmesh%cell(k,:))

   downdip: do i = 1,nz
      dw = real(i-1,pr)/real(nz-1,pr)  !downdip step 
! ! start point for line along strike
      xs = x(1)+(x(2)-x(1))*dw
      ys = y(1)+(y(2)-y(1))*dw
      zs = z(1)+(z(2)-z(1))*dw
! end point for line along strike 
      xe = x(4)+(x(3)-x(4))*dw
      ye = y(4)+(y(3)-y(4))*dw
      ze = z(4)+(z(3)-z(4))*dw

      alongstrike: do j = 1,nx !
        icount = icount+1
        dl = real(j-1,pr)/real(nx-1,pr) !step along strike
        tmp%px(icount) = xs+(xe-xs)*dl
        tmp%py(icount) = ys+(ye-ys)*dl
        tmp%pz(icount) = zs+(ze-zs)*dl     

        tmp_tri2quad(icount) = k 

      enddo alongstrike  
   enddo downdip

   if (output_level >= 2)  write(*,*)'est. cell, width, length',k,               &
                   sqrt((xe-xs)**2+(ye-ys)**2+(ze-zs)**2),                       &
                   sqrt((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
enddo large_cells


call random_seed ! intialise random seed 

! now link the nodes to cells and across the large cells. 
amesh%Ncells = (nx-1)*(nz-1)*2*qmesh%Ncells   
allocate(amesh%cell(amesh%Ncells,3))
amesh%cell = 0 
icount = 0 
c_count = 1
full_fault: do k = 1,qmesh%Ncells 

  !          write(*,*) 'k....',k
            qnode_1 = qmesh%cell(k,1)
            qnode_2 = qmesh%cell(k,2)
            qnode_3 = qmesh%cell(k,3)
            qnode_4 = qmesh%cell(k,4)
   !         write(*,*) 'qnodes....',qnode_1,qnode_2,qnode_3,qnode_4
            id_start = (k-1)*nx*nz+1  !
            id_end = (k-1)*nz*nx+nx*(nz-1) 
 !           write(*,*) 'start/end indices...',id_start,id_end,id_end-id_start
      subfault: do i = id_start, id_end     ! loop over all nodes in the quadrangle
                  cell_updated = .false.
                  if (mod(i,nx) == 0) cycle  ! don't take last point along strike in quadrangles as we deal with 2 nodes at a time for cells 
                  bc = .false.   
                  if (i <= id_start+nx-1)  then  !internal cells are along bottom boundary of quadrangle
 !                    write(*,*) 'i....',i,id_start,id_start+nx-1
                     bc = .true.
            ! check if edge point is already in the list
                     if (i==id_start) then ! checking bottom left corner
                           call bottom_left_corner(qnode_1,qnode_2,qnode_4,i,i+nx,i+nx+1,i+1,tmp,edge_x,edge_z,cell)
  !                         write(*,*) 'applying bottom_left_corner test using bottom edge........'
                           cell_updated = .true.
                           if ((any(cell(1,:) == 0)).or.(any(cell(2,:) == 0))) then 
                               write(*,*) 'problem bottom_left_corner'
                               stop 
                           endif 
                           
                     elseif (i == id_start+nx-2 ) then  ! check bottom right corner
                           call bottom_right_corner(qnode_1,qnode_3,qnode_4,i,i+nx,i+nx+1,i+1,tmp,edge_x,edge_z,cell)
   !                        write(*,*) 'applying bottom right corner test using bottom edge........'
                           cell_updated = .true.
                           if ((any(cell(1,:) == 0)).or.(any(cell(2,:) == 0))) then 
                               write(*,*)  cell(1,:)
                               write(*,*)  cell(2,:)
                               write(*,*) 'problem bottom_right_corner'
                               stop 
                           endif 
                     else  ! check the rest of the bottom edge 
                           call bottom_edge(qnode_1,qnode_4,i,i+nx,i+nx+1,i+1,tmp,edge_x,cell)
                           cell_updated = .true.
                       
                           if ((any(cell(1,:) == 0)).or.(any(cell(2,:) == 0))) then 
                               write(*,*) 'problem bottom_edge'
                               stop 
                           endif 
                     endif 
                     amesh%cell(c_count,:) = cell(1,:)
                     amesh%cell(c_count+1,:) = cell(2,:)
                     cell_updated = .true.
                  elseif (i > id_end-nx)  then  !along top boundary of quadrangle 
                     bc = .true.                  
                     if (i == id_end-nx+1 ) then ! check top left corner
                           call top_left_corner(qnode_1,qnode_2,qnode_3,i,i+nx,i+nx+1,i+1,tmp,edge_x,edge_z,cell) 
                           cell_updated = .true.
!                           write(*,*) 'applying top_left_corner'
                           if ((any(cell(1,:) == 0)).or.(any(cell(2,:) == 0))) then 
                              write(*,*) 'problem top_left_corner'
                              stop
                           endif 
                     elseif (i==id_end-1) then         !check top right corner 
                           call top_right_corner(qnode_2,qnode_3,qnode_4,i,i+nx,i+nx+1,i+1,tmp,edge_x,edge_z,cell)
                           cell_updated = .true.
                           if ((any(cell(1,:) == 0)).or.(any(cell(2,:) == 0))) then 
                              write(*,*) 'problem top_right_corner'
                              stop
                           endif 
                     else 
                           call top_edge(qnode_2,qnode_3,i,i+nx,i+nx+1,i+1,tmp,edge_x,cell) 
                           cell_updated = .true.
                           if ((any(cell(1,:) == 0)).or.(any(cell(2,:) == 0))) then 
                              write(*,*) 'problem top_edge'
                              stop
                           endif 
                     endif 
                     amesh%cell(c_count,:) = cell(1,:)
                     amesh%cell(c_count+1,:) = cell(2,:)
!                     cell_updated = .true.
                  elseif (mod(i,nx) == 1) then !cell is along left boundary 
                     ! write(*,*) 'left boundary start'
                     bc = .true.
                     call left_edge(qnode_1,qnode_2,i,i+nx,i+nx+1,i+1,tmp,edge_z,cell) 
                      if ((any(cell(1,:) == 0)).or.(any(cell(2,:) == 0))) then 
                              write(*,*) 'problem left_edge'
                              stop
                     endif 
                     amesh%cell(c_count,:) = cell(1,:)
                     amesh%cell(c_count+1,:) = cell(2,:)
   

                  elseif (mod(i+1,nx) == 0) then !cell is along right boundary
                     bc = .true.
                     ! write(*,*) 'right boundary start'
                      call  right_edge(qnode_3,qnode_4,i,i+nx,i+nx+1,i+1,tmp,edge_z,cell)
                     if ((any(cell(1,:) == 0)).or.(any(cell(2,:) == 0))) then 
                           write(*,*) 'problem right_edge'
                           stop 
                     endif 
                                           
                     amesh%cell(c_count,:) = cell(1,:)
                     amesh%cell(c_count+1,:) = cell(2,:)
                      
                  endif

                  if(.not.bc)  then ! cell is inside subfault 

                     ! call random_number(num)        ! randomly split of cell
                     ! if (num <0.5) then
                        amesh%cell(c_count,:)   =   (/ i, i+1, nx+i+1/)
                        amesh%cell(c_count+1,:) =   (/ i, nx+i+1, nx+i/)
                     ! else
                     !    amesh%cell(c_count,:)   =   (/ nx+i, i, i+1 /)
                     !    amesh%cell(c_count+1,:) =   (/ nx+i, i+1, i+1+nx /)
                     ! endif

                     ! amesh%cell(c_count,:)   =   (/ i, i+1, nx+i+1/)
                     ! amesh%cell(c_count+1,:) =   (/ i, nx+i+1, nx+i/)
                  endif
 
                  c_count = c_count+2
               enddo subfault
   enddo full_fault


  if (output_level >= 2)  write(*,*) 'number of triangular cells', amesh%Ncells 

! find and remove orphan nodes 
   nnodes = 0
   allocate(nodes(tmp%Nnodes))
   nodes = 0
   allocate(old_to_new(tmp%Nnodes))
   old_to_new = -1 


   !write(*,*) 'maxval in cell....', maxval(amesh%cell)
! loop over all nodes and check if they are needed
   node_loop: do i =1,tmp%Nnodes
         if ((any(amesh%cell(:,1)==i)).or.      &
             (any(amesh%cell(:,2)==i)).or.      &
             (any(amesh%cell(:,3)==i)) ) then !node is needed
            nnodes = nnodes+1
            nodes(nnodes) = i
            old_to_new(i) = nnodes
         else   ! node can be removed 
            old_to_new(i) = -1
         endif 
   enddo node_loop   
   
   ! correct the cell indexes
   do j = 1,amesh%Ncells
      do k = 1,3 
         if(old_to_new(amesh%cell(j,k))>0) then 
            amesh%cell(j,k) = old_to_new(amesh%cell(j,k))
         else 
            write(*,*) 'node not detected in cell',j,k
         endif 
      enddo

   enddo
 

! copy the nodes to the new mesh
   amesh%Nnodes = nnodes 
   allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
   amesh%px = tmp%px(nodes(1:nnodes))
   amesh%py = tmp%py(nodes(1:nnodes))
   amesh%pz = tmp%pz(nodes(1:nnodes))


! correct edge array indices 
   allocate(edge_to_do_x(qmesh%Nnodes,qmesh%Nnodes,nx),edge_to_do_z(qmesh%Nnodes,qmesh%Nnodes,nz))
   edge_to_do_x =.true.
   edge_to_do_z =.true.

   do i = 1,qmesh%Ncells
      qnode_1 = qmesh%cell(i,1)
      qnode_2 = qmesh%cell(i,2)
      qnode_3 = qmesh%cell(i,3)
      qnode_4 = qmesh%cell(i,4)

      do k = 1,nx
         ! if (edge_x(qnode_1,qnode_4 ,k) == 0 ) cycle   ! edge array is empty 
         ! if (edge_x(qnode_2,qnode_3 ,k) == 0 ) cycle   ! edge array is empty 
         if ((edge_x(qnode_1,qnode_4 ,k) /= 0 ).and.(edge_to_do_x(qnode_1,qnode_4,k)))  then
            if (old_to_new(edge_x(qnode_1,qnode_4 ,k)) > 0 ) then 
               edge_x(qnode_1,qnode_4 ,k) = old_to_new(edge_x(qnode_1,qnode_4 ,k))
               edge_x(qnode_4,qnode_1 ,k) = old_to_new(edge_x(qnode_4,qnode_1 ,k))
            else 
               old_index = edge_x(qnode_1,qnode_4 ,k)
               new_index = find_replacement_node(old_index,tmp,amesh)
               if (new_index == -1) then 
                  write(*,*) 'Error: no replacement node found (qnode_1,qnode_4)'
               else    ! replacement node found and applied 
                  edge_x(qnode_1,qnode_4 ,k) = new_index
                  edge_x(qnode_4,qnode_1 ,k) = new_index
               endif 
            endif 
            edge_to_do_x(qnode_1,qnode_4,k) = .false.
            edge_to_do_x(qnode_4,qnode_1,k) = .false.
         endif

         if ((edge_x(qnode_2,qnode_3 ,k) /= 0 ).and.(edge_to_do_x(qnode_2,qnode_3,k)) )then          
            if(old_to_new(edge_x(qnode_2,qnode_3 ,k)) > 0 ) then
               edge_x(qnode_2,qnode_3 ,k) = old_to_new(edge_x(qnode_2,qnode_3 ,k))
               edge_x(qnode_3,qnode_2 ,k) = old_to_new(edge_x(qnode_3,qnode_2 ,k))
            else  
               old_index = edge_x(qnode_2,qnode_3,k)
               new_index = find_replacement_node(old_index,tmp,amesh)
               if (new_index == -1) then 
                  write(*,*) 'Error: no replacement node found (qnode_2,qnode_3)'
                  write(*,*) 'index in question:    ',old_index
               else    
                  edge_x(qnode_2,qnode_3 ,k) = new_index
                  edge_x(qnode_3,qnode_2 ,k) = new_index
               endif 
            endif 
            edge_to_do_x(qnode_2,qnode_3,k) = .false.
            edge_to_do_x(qnode_3,qnode_2,k) = .false.
         endif

         if ((edge_z(qnode_1,qnode_2 ,k) /= 0 ).and.(edge_to_do_z(qnode_1,qnode_2,k))) then    ! edge array is empty 
            if(old_to_new(edge_z(qnode_1,qnode_2 ,k)) > 0 ) then 
               edge_z(qnode_1,qnode_2 ,k) = old_to_new(edge_z(qnode_1,qnode_2 ,k))
               edge_z(qnode_2,qnode_1 ,k) = old_to_new(edge_z(qnode_2,qnode_1 ,k))
            else 
               old_index = edge_z(qnode_1,qnode_2 ,k)
               new_index = find_replacement_node(old_index,tmp,amesh)
               if (new_index == -1) then 
                  write(*,*) 'Error: no replacement node found (qnode_1,qnode_2)'
                  write(*,*) 'index in question:    ',old_index
               else    
                  edge_z(qnode_1,qnode_2 ,k) = new_index
                  edge_z(qnode_2,qnode_1 ,k) = new_index
               endif 
            endif
            edge_to_do_z(qnode_1,qnode_2,k ) = .false.
            edge_to_do_z(qnode_2,qnode_1,k ) = .false.
         endif 



         if ((edge_z(qnode_3,qnode_4 ,k ) /= 0 ).and.(edge_to_do_z(qnode_3,qnode_4,k))) then    ! edge array is empty 
            if(old_to_new(edge_z(qnode_3,qnode_4 ,k)) > 0 ) then
                edge_z(qnode_4,qnode_3 ,k) = old_to_new(edge_z(qnode_4,qnode_3 ,k))
               edge_z(qnode_3,qnode_4 ,k) = old_to_new(edge_z(qnode_3,qnode_4 ,k)) 
            else 
               old_index = edge_z(qnode_3,qnode_4 ,k)
               new_index = find_replacement_node(old_index,tmp,amesh)
               if (new_index == -1) then 
                  write(*,*) 'Error: no replacement node found (qnode_3,qnode_4)'
                  write(*,*) 'index in question:    ',old_index
               else    
                  edge_z(qnode_3,qnode_4 ,k) = new_index
                  edge_z(qnode_4,qnode_3 ,k) = new_index
               endif 
            endif
            edge_to_do_z(qnode_3,qnode_4,k) = .false.
            edge_to_do_z(qnode_4,qnode_3,k) = .false.
         endif 
      enddo


   enddo

   ! write(*,*) 'number of original nodes', tmp%Nnodes
   ! write(*,*) 'number of nodes', nnodes
   ! write(*,*) 'new maxval in cell....', maxval(amesh%cell)



   allocate(tri2quad(amesh%Nnodes))
   tri2quad = tmp_tri2quad(nodes(1:nnodes))
   ! set the values on the edge to zero  
   do i = 1,qmesh%Nnodes
      do k = i,qmesh%Nnodes
         if (edge_x(i,k ,1) /= 0 ) then  ! edge array is not empty   
            tri2quad(edge_x(i,k,:)) = 0
         endif  
         if (edge_z(i,k ,1) /= 0 ) then    ! edge array is not empty   
            tri2quad(edge_z(i,k,:)) = 0
         endif 
      enddo 
   enddo

return 
end subroutine create_tri_mesh
!###############################################################################
subroutine bottom_left_corner(q1,q2,q4,n1,n2,n3,n4,tmp,edge_x,edge_z,cell) 
type(mesh),intent(in) :: tmp
integer, intent(in) :: q1,q2,q4,n1,n2,n3,n4
integer,dimension(:,:,:),intent(inout) :: edge_x,edge_z
integer,dimension(2,3),intent(out) :: cell
integer :: pre_n1,pre_n2,pre_n4,id 
!
!
!        pre_n2
!      .n2------------.n3
!      |              |
!      |              |
! -----.n1------------.n4
!      | pre_n1
!      |
pre_n1 = check_edge(q1,q2,edge_z,tmp,n1) !left edge
pre_n2 = check_edge(q1,q2,edge_z,tmp,n2) !left edge
pre_n4 = check_edge(q1,q4,edge_x,tmp,n4) !bottom edge 
! write(*,* ) 'bottom_left_corner'
! write(*,*) 'n1,n2,n3,n4',n1,n2,n3,n4
! write (*,*) 'n1:     ',tmp%px(n1),tmp%py(n1),tmp%pz(n1)
! write (*,*) 'n2:     ',tmp%px(n2),tmp%py(n2),tmp%pz(n2)
! write(*,*)'pre_n1,pre_n2,pre_n4',pre_n1,pre_n2,pre_n4
! write(*,*) 'edge_z....',edge_z(q1,q2,:)

if ((pre_n1==0).and.(pre_n2==0).and.(pre_n4==0)) then ! all new nodes 
   cell(1,:) = (/ n1, n4, n3/)
   cell(2,:) = (/ n1, n3, n2/)
! update both edges 
   id = find_empty_position(edge_z(q1,q2,:)) ! left edge 
   edge_z(q2,q1,id) = n1
   edge_z(q1,q2,id) = n1

   edge_z(q2,q1,id+1) = n2
   edge_z(q1,q2,id+1) = n2

   id = find_empty_position(edge_x(q1,q4,:))  ! bottom edge 
   edge_x(q4,q1,id) = n1
   edge_x(q1,q4,id) = n1

   edge_x(q4,q1,id+1) = n4
   edge_x(q1,q4,id+1) = n4

elseif ((pre_n1/=0).and.(pre_n2/=0).and.(pre_n4==0)) then ! left side already exists 
   cell(1,:) = (/ pre_n1, n4, n3/)
   cell(2,:) = (/ pre_n1, n3, pre_n2/)
!update bottom edge 
   id = find_empty_position(edge_x(q1,q4,:))
   edge_x(q4,q1,id) = n1
   edge_x(q1,q4,id) = n1

   edge_x(q4,q1,id+1) = n4
   edge_x(q1,q4,id+1) = n4

elseif ((pre_n1 /= 0).and.(pre_n2 == 0).and.(pre_n4 /=0)) then ! bottom already exists 
   cell(1,:) = (/ pre_n1, pre_n4, n3/)
   cell(2,:) = (/ pre_n1, n3, n2/)

! update left edge 
   id = find_empty_position(edge_z(q1,q2,:))
   edge_z(q2,q1,id) = n1
   edge_z(q1,q2,id) = n1

   edge_z(q2,q1,id+1) = n2
   edge_z(q1,q2,id+1) = n2

elseif ((pre_n1/=0).and.(pre_n2==0).and.(pre_n4 ==0)) then ! bottom corner pre-exists (I don't think we ever reach this condition)
   cell(1,:) = (/ pre_n1, n4, n3/)
   cell(2,:) = (/ pre_n1, n3, n2/)

! add corner to both edges 
   id = find_empty_position(edge_z(q1,q2,:))
   edge_z(q2,q1,id) = n1
   edge_z(q1,q2,id) = n1

   ! id = find_empty_position(edge_x)  !remove as we haven't tested this 
   ! edge_z(q4,q1,id) = n1
   ! edge_z(q1,q4,id) = n1

elseif ((pre_n1 /= 0).and.(pre_n2/=0).and.(pre_n4 /=0)) then ! all sides defined    
   cell(1,:) = (/ pre_n1, pre_n4, n3/)
   cell(2,:) = (/ pre_n1, n3, pre_n2/)   
endif
return
end subroutine  bottom_left_corner
!###############################################################################
subroutine bottom_right_corner(q1,q3,q4,n1,n2,n3,n4,tmp,edge_x,edge_z,cell) 
type(mesh),intent(in) :: tmp
integer, intent(in) :: q1,q3,q4,n1,n2,n3,n4
integer,dimension(:,:,:),intent(inout) :: edge_x,edge_z
integer,dimension(2,3),intent(out) :: cell
integer :: pre_n1,pre_n3,pre_n4,pre_n4_a, id
!      
!                     |
!                     |pre_n3
!      .n2------------.n3 
!      |              |
!      |              |
!      .n1------------.n4 -------
!      pre_n1         |pre_n4
!                     |


pre_n1 = check_edge(q1,q4,edge_x,tmp,n1) !bottom edge 
pre_n3 = check_edge(q3,q4,edge_z,tmp,n3) ! right edge
pre_n4 = check_edge(q3,q4,edge_z,tmp,n4) !right edge  
! write(*,* ) 'bottom_right_corner'
! write(*,*) 'n1,n2,n3,n4',n1,n2,n3,n4
! write (*,*) 'n4:     ',tmp%px(n4),tmp%py(n4),tmp%pz(n4)
! write (*,*) 'n3:     ',tmp%px(n3),tmp%py(n3),tmp%pz(n3)
! write(*,*)'pre_n1,pre_n2,pre_n4',pre_n1,pre_n3,pre_n4
! write(*,*) 'edge_z....',edge_z(q3,q4,:)

cell = 0 
if ((pre_n1==0).and.(pre_n3==0).and.(pre_n4==0)) then ! all new nodes 
   cell(1,:) = (/ n1, n4, n3/)
   cell(2,:) = (/ n1, n3, n2/)

   ! cell(1,:) = (/ n2, n1, n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ n2, n4, n3/)

! update both edges 
   id = find_empty_position(edge_z(q3,q4,:))  ! right edge 
   edge_z(q3,q4,id) = n4
   edge_z(q4,q3,id) = n4

   edge_z(q3,q4,id+1) = n3
   edge_z(q4,q3,id+1) = n3

   id = find_empty_position(edge_x(q1,q4,:))   ! bottom edge 
   edge_x(q4,q1,id) = n1
   edge_x(q1,q4,id) = n1

   edge_x(q4,q1,id+1) = n4  ! we haven't checked this value but if n1 is not in the list then neither should n4
   edge_x(q1,q4,id+1) = n4

elseif ((pre_n1==0).and.(pre_n3/=0).and.(pre_n4/=0)) then ! right side already logged 
   cell(1,:) = (/ n1, pre_n4, pre_n3/)
   cell(2,:) = (/ n1, pre_n3, n2/)

   ! cell(1,:) = (/ n2, n1, pre_n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ n2, pre_n4, pre_n3/)


   !add nodes to bottom edge 
   id = find_empty_position(edge_x(q1,q4,:))   ! bottom edge 
   edge_x(q4,q1,id) = n1
   edge_x(q1,q4,id) = n1

   pre_n4_a = check_edge(q1,q4,edge_x,tmp,n4) ! we havent checked this point on other edge 
   if (pre_n4_a == 0) then 
      edge_x(q4,q1,id+1) = n4  
      edge_x(q1,q4,id+1) = n4
   endif 
elseif ((pre_n1 /= 0).and.(pre_n3 == 0).and.(pre_n4 ==0)) then ! bottom already exists 
   cell(1,:) = (/ pre_n1, n4, n3/)
   cell(2,:) = (/ pre_n1, n3, n2/)

   ! cell(1,:) = (/ n2, pre_n1, n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ n2, n4, n3/)

!add nodes to the right edge
   id = find_empty_position(edge_z(q3,q4,:))  ! right edge 
   edge_z(q3,q4,id) = n4
   edge_z(q4,q3,id) = n4

   edge_z(q3,q4,id+1) = n3
   edge_z(q4,q3,id+1) = n3

  pre_n4_a = check_edge(q1,q4,edge_x,tmp,n4) ! we havent checked this point on other edge 
   if (pre_n4_a == 0) then 
      id = find_empty_position(edge_x(q1,q4,:))   ! bottom edge 
      edge_x(q4,q1,id) = n4  
      edge_x(q1,q4,id) = n4
   endif 

elseif ((pre_n1==0).and.(pre_n3==0).and.(pre_n4 /= 0)) then ! bottom right corner pre-exists 
   cell(1,:) = (/ n1, pre_n4, n3/)
   cell(2,:) = (/ n1, n3, n2/)

   ! cell(1,:) = (/ n2, n1, pre_n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ n2, pre_n4, n3/)


   id = find_empty_position(edge_z(q3,q4,:))  ! right edge 
   edge_z(q3,q4,id) = n4
   edge_z(q4,q3,id) = n4

elseif ((pre_n1 /= 0).and.(pre_n3/=0).and.(pre_n4 /=0)) then ! all sides defined    
   cell(1,:) = (/ pre_n1, pre_n4, pre_n3/)
   cell(2,:) = (/ pre_n1, pre_n3, n2/)   

   ! cell(1,:) = (/ n2, pre_n1, pre_n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ n2, pre_n4, pre_n3/)


endif

! write(*,*) 'edge_z(end of bottom_right_corner)....',q3,q4,edge_z(q3,q4,:)
! write(*,*) 'edge_x(end of bottom_right_corner)....',q1,q4,edge_x(q1,q4,:)

return
end subroutine bottom_right_corner
!###############################################################################
subroutine top_left_corner(q1,q2,q3,n1,n2,n3,n4,tmp,edge_x,edge_z,cell) 
type(mesh),intent(in) :: tmp
integer, intent(in) :: q1,q2,q3,n1,n2,n3,n4
integer,dimension(:,:,:),intent(inout) :: edge_x,edge_z
integer,dimension(2,3),intent(out) :: cell
integer :: pre_n1,pre_n2,pre_n3,pre_n2_a,pre_n2_top, id

!
!      |
!      |
!      |pre_n2       pre_n3
!------.n2------------.n3 -------
!      |              |
!      |              |
!------.n1------------.n4 
!      |pre_n1       
!      |

pre_n1 = check_edge(q1,q2,edge_z,tmp,n1) !left edge 
pre_n2 = check_edge(q1,q2,edge_z,tmp,n2) ! left edge
pre_n3 = check_edge(q2,q3,edge_x,tmp,n3) ! top edge 

pre_n2_top = check_edge(q2,q3,edge_x,tmp,n2) ! top edge 


! write(*,* ) 'top_left_corner'
! write(*,*) 'n1,n2,n3,n4',n1,n2,n3,n4
! write (*,*) 'n4:     ',tmp%px(n2),tmp%py(n2),tmp%pz(n2)
! write (*,*) 'n3:     ',tmp%px(n3),tmp%py(n3),tmp%pz(n3)
! write(*,*)'pre_n1,pre_n2,pre_n3',pre_n1,pre_n2,pre_n3
! write(*,*) 'edge_x....',q2,q3,edge_x(q2,q3,:)
! write(*,*) 'edge_z....',q1,q2,edge_z(q1,q2,:)


cell = 0
if ((pre_n1==0).and.(pre_n2==0).and.(pre_n3==0)) then ! all new nodes 
   cell(1,:) = (/ n1, n4, n3/)
   cell(2,:) = (/ n1, n3, n2/)

   ! cell(1,:) = (/ n2, n1, n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ n2, n4, n3/)

! update both edge 
   id = find_empty_position(edge_z(q1,q2,:))  ! left edge 
   edge_z(q1,q2,id) = n1
   edge_z(q2,q1,id) = n1

   edge_z(q1,q2,id+1) = n2
   edge_z(q2,q1,id+1) = n2

   id = find_empty_position(edge_x(q2,q3,:))   ! top edge 
   edge_x(q2,q3,id) = n2 ! we haven't checked this value but if n1 is not in the list then neither should n3
   edge_x(q3,q2,id) = n2

   edge_x(q2,q3,id+1) = n3  
   edge_x(q3,q2,id+1) = n3


elseif ((pre_n1/=0).and.(pre_n2/=0).and.(pre_n3==0)) then ! left side already exists 
   cell(1,:) = (/ pre_n1, n4, n3/)
   cell(2,:) = (/ pre_n1, n3, pre_n2/)

   ! cell(1,:) = (/ pre_n2, pre_n1, n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ pre_n2, n4, n3/)

   id = find_empty_position(edge_x(q2,q3,:))   ! update top edge 
   edge_x(q2,q3,id) = n2 ! we haven't checked this value but if n1 is not in the list then neither should n3
   edge_x(q3,q2,id) = n2

   edge_x(q2,q3,id+1) = n3  
   edge_x(q3,q2,id+1) = n3

elseif ((pre_n1 == 0).and.(pre_n2 /= 0).and.(pre_n3 /=0)) then ! top already exists 
   cell(1,:) = (/ n1, n4, pre_n3/)
   cell(2,:) = (/ n1, pre_n3, pre_n2/)

   ! cell(1,:) = (/ pre_n2, n1, n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ pre_n2, n4, pre_n3/)


! update both edge 
   id = find_empty_position(edge_z(q1,q2,:))  ! left edge 
   edge_z(q1,q2,id) = n1
   edge_z(q2,q1,id) = n1

   edge_z(q1,q2,id+1) = n2
   edge_z(q2,q1,id+1) = n2


elseif ((pre_n1 ==0).and.(pre_n2/=0).and.(pre_n3 ==0)) then ! top corner pre-exists 
   cell(1,:) = (/ n1, n4, n3/)
   cell(2,:) = (/ n1, n3, pre_n2/)

   ! cell(1,:) = (/ pre_n2, n1, n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ pre_n2, n4, n3/)

   id = find_empty_position(edge_z(q1,q2,:))  ! left edge 
   edge_z(q1,q2,id) = n1
   edge_z(q2,q1,id) = n1

   id = find_empty_position(edge_x(q2,q3,:))   ! top edge 
   edge_x(q2,q3,id) = n3  
   edge_x(q3,q2,id) = n3

elseif((pre_n1 /=0).and.(pre_n2/=0).and.(pre_n3 /=0)) then ! all sides already defined  
   cell(1,:) = (/ pre_n1, n4, pre_n3/)
   cell(2,:) = (/ pre_n1, pre_n3, pre_n2/)   

   ! cell(1,:) = (/ pre_n2, pre_n1, n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ pre_n2, n4, pre_n3/)

elseif((pre_n1 /=0).and.(pre_n2==0).and.(pre_n3 ==0)) then ! left bottom point pre-exists 
   cell(1,:) = (/ pre_n1, n4, n3/)
   cell(2,:) = (/ pre_n1, n3, n2/)

   ! cell(1,:) = (/ n2, pre_n1, n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ n2, n4, n3/)


! update both edge 
   id = find_empty_position(edge_z(q1,q2,:))  ! left edge corner value added  
   edge_z(q1,q2,id) = n2
   edge_z(q2,q1,id) = n2

   id = find_empty_position(edge_x(q2,q3,:))   ! top edge 
   edge_x(q2,q3,id) = n2 ! we haven't checked this value but if n1 is not in the list then neither should n3
   edge_x(q3,q2,id) = n2

   edge_x(q2,q3,id+1) = n3  
   edge_x(q3,q2,id+1) = n3


elseif((pre_n1 ==0).and.(pre_n2==0).and.(pre_n3 /=0)) then ! top right point pre-exists 
   cell(1,:) = (/ n1, n4, pre_n3/)
   cell(2,:) = (/ n1, pre_n3, n2/)

   ! cell(1,:) = (/ n2, n1, n4/)  !use this orientation so that slipping area can be defined correctly
   ! cell(2,:) = (/ n2, n4, pre_n3/)


   id = find_empty_position(edge_x(q2,q3,:))   ! top edge 

   pre_n2_a = check_edge(q2,q3,edge_x,tmp,n2) ! we havent checked this point on other edge 
   if (pre_n2_a == 0) then 
    
      edge_x(q2,q3,id) = n2  
      edge_x(q3,q2,id) = n2

      edge_x(q2,q3,id+1) = n3  
      edge_x(q3,q2,id+1) = n3

   else 
      edge_x(q2,q3,id) = n3  
      edge_x(q3,q2,id) = n3
   endif 


elseif((pre_n1 /=0).and.(pre_n2==0).and.(pre_n3 /=0)) then ! corner point not defined pre-exists 

   if (pre_n2_top /= 0) then ! check if other edge has a node 
     cell(1,:) = (/ pre_n1, n4, pre_n3/)
     cell(2,:) = (/ pre_n1, pre_n3, pre_n2_top/)  

   !   cell(1,:) = (/ pre_n2_top, pre_n1, n4/)  !use this orientation so that slipping area can be defined correctly
   !   cell(2,:) = (/ pre_n2_top, n4, pre_n3/)

      id = find_empty_position(edge_z(q1,q2,:))  ! left edge, assume we don't need to update top edge as pre_n3 is present
      edge_z(q1,q2,id) = pre_n2_top
      edge_z(q2,q1,id) = pre_n2_top
   else 
     cell(1,:) = (/ pre_n1, n4, pre_n3/)
     cell(2,:) = (/ pre_n1, pre_n3, n2/)  

   !   cell(1,:) = (/ n2, pre_n1, n4/)  !use this orientation so that slipping area can be defined correctly
   !   cell(2,:) = (/ n2, n4, pre_n3/)

      id = find_empty_position(edge_z(q1,q2,:))  ! left edge, assume we don't need to update top edge as pre_n3 is present
      edge_z(q1,q2,id) = n2
      edge_z(q2,q1,id) = n2
   endif 

endif

if (cell(1,1) == 0)  then 
   write(*,*) 'problem top_left_corner'
   write(*,*) n1,n2,n3,n4
   write(*,*) pre_n1,pre_n2,pre_n3
   stop
endif 

return
end subroutine top_left_corner
!###############################################################################
subroutine top_right_corner(q2,q3,q4,n1,n2,n3,n4,tmp,edge_x,edge_z,cell)
type(mesh),intent(in) :: tmp
integer, intent(in) :: n1,n2,n3,n4,q2,q3,q4
integer,dimension(:,:,:),intent(inout) :: edge_x,edge_z
integer,dimension(2,3),intent(out) :: cell
integer :: pre_n4,pre_n2,pre_n3,id,pre_n3_a,pre_n3_top
!
!      
!                     |      
!       pre_n2        | pre_n3
!      .n2------------.n3 -------
!      |              |
!      |              |
!      .n1------------.n4 
!                     | pre_n4       
!                     |

pre_n2 = check_edge(q2,q3,edge_x,tmp,n2) !top edge 
pre_n3_top = check_edge(q2,q3,edge_x,tmp,n3) !top edge 

pre_n3 = check_edge(q3,q4,edge_z,tmp,n3) ! right edge
pre_n4 = check_edge(q3,q4,edge_z,tmp,n4) ! right edge 

! write(*,* ) 'top_right_corner'
! write(*,*) 'n1,n2,n3,n4',n1,n2,n3,n4
! write (*,*) 'n4:     ',tmp%px(n2),tmp%py(n2),tmp%pz(n2)
! write (*,*) 'n3:     ',tmp%px(n3),tmp%py(n3),tmp%pz(n3)
! write(*,*)'pre_n2,pre_n3,pre_n4',pre_n2,pre_n3,pre_n4
! write(*,*) 'edge_x....',edge_x(q2,q3,:)
! write(*,*) 'edge_z(top_right_corner)....',q3,q4,edge_z(q3,q4,:)

cell = 0
if ((pre_n2==0).and.(pre_n3==0).and.(pre_n4==0)) then ! all new nodes 
   cell(1,:) = (/ n1, n4, n3/)
   cell(2,:) = (/ n1, n3, n2/)

! update both edge 
   id = find_empty_position(edge_z(q3,q4,:))  ! update right edge 
   edge_z(q3,q4,id) = n4
   edge_z(q4,q3,id) = n4

   edge_z(q3,q4,id+1) = n3
   edge_z(q4,q3,id+1) = n3

   id = find_empty_position(edge_x(q2,q3,:))   ! update top edge 
   edge_x(q2,q3,id) = n2 
   edge_x(q3,q2,id) = n2

   ! pre_n3_a = check_edge(q2,q3,edge_x,tmp,n3) ! we havent checked this point on other edge 
   ! if (pre_n3_a == 0) then 
   if (pre_n3_top == 0) then 
      edge_x(q2,q3,id+1) = n3  
      edge_x(q3,q2,id+1) = n3
   endif 

elseif ((pre_n2 ==0).and.(pre_n3/=0).and.(pre_n4/=0)) then ! right side already exists 

   cell(1,:) = (/ n1, pre_n4, pre_n3/)
   cell(2,:) = (/ n1, pre_n3, n2/)

   id = find_empty_position(edge_x(q2,q3,:))   ! update top edge 
   edge_x(q2,q3,id) = n2 
   edge_x(q3,q2,id) = n2

   if (pre_n3_top == 0) then 
      edge_x(q2,q3,id+1) = pre_n3 
      edge_x(q3,q2,id+1) = pre_n3
   endif 

   ! id = find_empty_position(edge_z(q3,q4,:))  ! update right edge 
   ! edge_z(q3,q4,id) = n4
   ! edge_z(q4,q3,id) = n4

   ! edge_z(q3,q4,id+1) = n3
   ! edge_z(q4,q3,id+1) = n3


elseif ((pre_n2 == 0).and.(pre_n3 /= 0).and.(pre_n4 ==0)) then ! top already exists 
   write(*,*) 'got here',q2,q3
   cell(1,:) = (/ n1, n4, pre_n3/)
   cell(2,:) = (/ n1, pre_n3, pre_n2/)

   id = find_empty_position(edge_x(q2,q3,:))   ! update top edge 
   edge_x(q2,q3,id) = n2 
   edge_x(q3,q2,id) = n2

   ! pre_n3_a = check_edge(q2,q3,edge_x,tmp,n3) ! we havent checked this point on other edge 
   ! if (pre_n3_a == 0) then 
   if (pre_n3_top == 0) then 
      edge_x(q2,q3,id+1) = n3  
      edge_x(q3,q2,id+1) = n3
   endif 

   id = find_empty_position(edge_z(q3,q4,:))  ! update right edge 
   edge_z(q3,q4,id) = n4
   edge_z(q4,q3,id) = n4
! elseif ((pre_n2 == 0).and.(pre_n3/=0).and.(pre_n4 ==0)) then ! right corner pre-exists 
!    cell(1,:) = (/ n1, n4, pre_n3/)
!    cell(2,:) = (/ n1, pre_n3, n2/)

!    id = find_empty_position(edge_x(q2,q3,:))   ! update top edge 
!    edge_x(q2,q3,id) = n2 
!    edge_x(q3,q2,id) = n2

!    pre_n3_a = check_edge(q2,q3,edge_x,tmp,n3) ! we havent checked this point on other edge 
!    if (pre_n3_a == 0) then 
!       edge_x(q2,q3,id+1) = n3  
!       edge_x(q3,q2,id+1) = n3
!    endif 

   ! id = find_empty_position(edge_z(q3,q4,:))  ! update right edge 
   ! edge_z(q3,q4,id) = n4
   ! edge_z(q4,q3,id) = n4

elseif ((pre_n2 /= 0).and.(pre_n3/=0).and.(pre_n4 /=0)) then ! all sides defined    
   cell(1,:) = (/ n1, pre_n4, pre_n3/)
   cell(2,:) = (/ n1, pre_n3, pre_n2/)


elseif ((pre_n2 /= 0).and.(pre_n3==0).and.(pre_n4 /=0)) then ! corner not saved   

   if(pre_n3_top /= 0) then ! use value from other edge 
      cell(1,:) = (/ n1, pre_n4, pre_n3_top/)
      cell(2,:) = (/ n1, pre_n3_top, pre_n2/)

      id = find_empty_position(edge_z(q3,q4,:))  ! update right edge 
      edge_z(q3,q4,id) = pre_n3_top
      edge_z(q4,q3,id) = pre_n3_top
   
   else   !if there is no value on the other edge take new node 
      cell(1,:) = (/ n1, pre_n4, n3/)
      cell(2,:) = (/ n1, n3, pre_n2/)
   
      id = find_empty_position(edge_z(q3,q4,:))  ! update right edge 
      edge_z(q3,q4,id) = n3
      edge_z(q4,q3,id) = n3

      id = find_empty_position(edge_x(q2,q3,:))   ! update top edge corner
      edge_x(q2,q3,id) = n3 
      edge_x(q3,q2,id) = n3

   endif 
endif

if (cell(1,1) == 0) then 
   write(*,*) 'problem top_right_corner'
   write(*,*) n1,n2,n3,n4
   write(*,*) pre_n2,pre_n3,pre_n4   
   stop
endif 

return 
end subroutine top_right_corner
!###############################################################################
subroutine bottom_edge(q1,q4,n1,n2,n3,n4,tmp,edge,cell) 
type(mesh),intent(in) :: tmp
integer, intent(in) :: n1,n2,n3,n4,q1,q4
integer,dimension(:,:,:),intent(inout) :: edge
integer :: pre_n4,pre_n1
integer,dimension(2,3) :: cell
real(pr) :: num 
integer :: id 
!
!      
!      .n2------------.n3 
!      |              |
!      |              |
! -----.n1------------.n4 --------------
!      pre_n1         pre_n4       
!                     
pre_n1 = check_edge(q1,q4,edge,tmp,n1) ! bottom of quadrangle
pre_n4 = check_edge(q1,q4,edge,tmp,n4) !      

cell = 0
if ((pre_n1==0).and.(pre_n4==0)) then ! all new nodes 
   ! cell(1,:) = (/ n1, n4, n3/)
   ! cell(2,:) = (/ n1, n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, n4, n3/)
      cell(2,:) = (/ n1, n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, n1, n4/)
   !    cell(2,:) = (/ n2, n4, n3/)
   ! endif

   id = find_empty_position(edge(q1,q4,:))
   edge(q4,q1,id) = n1
   edge(q1,q4,id) = n1
  
   edge(q4,q1,id+1) = n4
   edge(q1,q4,id+1) = n4

elseif ((pre_n1 /=0).and.(pre_n4 /=0)) then !both nodes are already logged
   ! cell(1,:) = (/ pre_n1, pre_n4, n3/)
   ! cell(2,:) = (/ pre_n1, n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ pre_n1, pre_n4, n3/)
      cell(2,:) = (/ pre_n1, n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, pre_n1, pre_n4/)
   !    cell(2,:) = (/ n2, pre_n4, n3/)
   ! endif


elseif ((pre_n1 /=0).and.(pre_n4 ==0)) then ! nodes at n1 is already logged
   ! cell(1,:) = (/ pre_n1, n4, n3/)
   ! cell(2,:) = (/ pre_n1, n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ pre_n1, n4, n3/)
      cell(2,:) = (/ pre_n1, n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, pre_n1, n4/)
   !    cell(2,:) = (/ n2, n4, n3/)
   ! endif

   id = find_empty_position(edge(q1,q4,:))  
   edge(q4,q1,id) = n4
   edge(q1,q4,id) = n4

elseif ((pre_n1 == 0).and.(pre_n4 /=0)) then !nodes at np is already logged
   ! cell(1,:) = (/ n1, pre_n4, n3/)
   ! cell(2,:) = (/ n1, n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, pre_n4, n3/)
      cell(2,:) = (/ n1, n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, n1, pre_n4/)
   !    cell(2,:) = (/ n2, pre_n4, n3/)
   ! endif

   id = find_empty_position(edge(q1,q4,:))
   edge(q4,q1,id) = n1
   edge(q1,q4,id) = n1

endif 

return
end subroutine bottom_edge 
!###############################################################################
subroutine top_edge(q2,q3,n1,n2,n3,n4,tmp,edge,cell) 
type(mesh),intent(in) :: tmp
integer, intent(in) :: n1,n2,n3,n4,q2,q3
integer,dimension(:,:,:),intent(inout) :: edge
integer :: pre_n2,pre_n3
integer,dimension(2,3) :: cell
real(pr) :: num
integer :: id 
!
!       pre_n2         pre_n3        
! -----.n2------------.n3 -------------- 
!      |              |
!      |              |
!      .n1------------.n4
!    
!

pre_n2 = check_edge(q2,q3,edge,tmp,n2) !top of quadrangle 
pre_n3 = check_edge(q2,q3,edge,tmp,n3) 


cell = 0
if ((pre_n2==0).and.(pre_n3==0)) then ! all new nodes 
   ! cell(1,:) = (/ n1, n4, n3/)
   ! cell(2,:) = (/ n1, n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, n4, n3/)
      cell(2,:) = (/ n1, n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, n1, n4/)
   !    cell(2,:) = (/ n2, n4, n3/)
   ! endif


! add all nodes to edge 
   id = find_empty_position(edge(q2,q3,:))
   edge(q2,q3,id) = n2
   edge(q3,q2,id) = n2
  
   edge(q2,q3,id+1) = n3
   edge(q3,q2,id+1) = n3

elseif ((pre_n2/=0).and.(pre_n3/=0)) then!both nodes are already logged 
   ! cell(1,:) = (/ n1, n4, pre_n3/)
   ! cell(2,:) = (/ n1, pre_n3, pre_n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, n4, pre_n3/)
      cell(2,:) = (/ n1, pre_n3, pre_n2/)
   ! else
   !    cell(1,:) = (/ pre_n2, n1, n4/)
   !    cell(2,:) = (/ pre_n2, n4, pre_n3/)
   ! endif


elseif ((pre_n2 /=0).and.(pre_n3 == 0)) then ! nodes at n2 is already logged
   ! cell(1,:) = (/ n1, n4, n3/)
   ! cell(2,:) = (/ n1, n3, pre_n2/)
   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, n4, n3/)
      cell(2,:) = (/ n1, n3, pre_n2/)
   ! else
   !    cell(1,:) = (/ pre_n2, n1, n4/)
   !    cell(2,:) = (/ pre_n2, n4, n3/)
   ! endif


   id = find_empty_position(edge(q2,q3,:))
   edge(q2,q3,id) = n3
   edge(q3,q2,id) = n3

elseif ((pre_n2 == 0).and.(pre_n3 /=0)) then !nodes at n3 is already logged
   ! cell(1,:) = (/ n1, n4, pre_n3/)
   ! cell(2,:) = (/ n1, pre_n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, n4, pre_n3/)
      cell(2,:) = (/ n1, pre_n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, n1, n4/)
   !    cell(2,:) = (/ n2, n4, pre_n3/)
   ! endif

   id = find_empty_position(edge(q2,q3,:))
   edge(q2,q3,id) = n2
   edge(q3,q2,id) = n2

endif 

return
end subroutine top_edge 
!###############################################################################
subroutine left_edge(q1,q2,n1,n2,n3,n4,tmp,edge,cell) 
type(mesh),intent(in) :: tmp
integer, intent(in) :: n1,n2,n3,n4,q1,q2
integer,dimension(:,:,:),intent(inout) :: edge
integer :: pre_n1,pre_n2
integer,dimension(2,3) :: cell
integer :: id 
real(pr) :: num 
!      |
!      |pre_n2                
! -----.n2------------.n3 
!      |              |
!      |              |
!      .n1------------.n4
!      |pre_n1     
!      |
pre_n1 = check_edge(q1,q2,edge,tmp,n1) 
pre_n2 = check_edge(q1,q2,edge,tmp,n2) 

cell = 0
if ((pre_n1==0).and.(pre_n2==0)) then ! all new nodes 
   ! cell(1,:) = (/ n1, n4, n3/)
   ! cell(2,:) = (/ n1, n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, n4, n3/)
      cell(2,:) = (/ n1, n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, n1, n4/)
   !    cell(2,:) = (/ n2, n4, n3/)
   ! endif

! add all nodes to edge 
   id = find_empty_position(edge(q1,q2,:))
   edge(q1,q2,id) = n1
   edge(q2,q1,id) = n1
  
   edge(q1,q2,id+1) = n2
   edge(q2,q1,id+1) = n2

elseif ((pre_n1/=0).and.(pre_n2/=0)) then!both nodes are already logged 
   ! cell(1,:) = (/ pre_n1, n4, n3/)
   ! cell(2,:) = (/ pre_n1, n3, pre_n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ pre_n1, n4, n3/)
      cell(2,:) = (/ pre_n1, n3, pre_n2/)
   ! else
   !    cell(1,:) = (/ pre_n2, pre_n1, n4/)
   !    cell(2,:) = (/ pre_n2, n4, n3/)
   ! endif


elseif ((pre_n1 /=0).and.(pre_n2 == 0)) then ! node at n1 is already logged
   ! cell(1,:) = (/ pre_n1, n4, n3/)
   ! cell(2,:) = (/ pre_n1, n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ pre_n1, n4, n3/)
      cell(2,:) = (/ pre_n1, n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, pre_n1, n4/)
   !    cell(2,:) = (/ n2, n4, n3/)
   ! endif

   id = find_empty_position(edge(q1,q2,:))
   edge(q1,q2,id) = n2
   edge(q2,q1,id) = n2

elseif ((pre_n1 == 0).and.(pre_n2 /=0)) then !nodes at n2 is already logged
   ! cell(1,:) = (/ n1, n4, n3/)
   ! cell(2,:) = (/ n1, n3, pre_n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, n4, n3/)
      cell(2,:) = (/ n1, n3, pre_n2/)
   ! else
   !    cell(1,:) = (/ pre_n2, n1, n4/)
   !    cell(2,:) = (/ pre_n2, n4, n3/)
   ! endif

   id = find_empty_position(edge(q1,q2,:))
   edge(q1,q2,id) = n1
   edge(q2,q1,id) = n1

endif 


return 
end subroutine left_edge
!###############################################################################
subroutine right_edge(q3,q4,n1,n2,n3,n4,tmp,edge,cell) 
type(mesh),intent(in) :: tmp
integer, intent(in) :: n1,n2,n3,n4,q3,q4
integer,dimension(:,:,:),intent(inout) :: edge
integer :: pre_n3,pre_n4
integer,dimension(2,3) :: cell
integer :: id 
real(pr) :: num 
!                     |
!                     |pre_n3                
!      .n2------------.n3------------ 
!      |              |
!      |              |
!      .n1------------.n4
!                     | pre_n4     
!                     |
pre_n3 = check_edge(q3,q4,edge,tmp,n3) !left side of quadrangle 
pre_n4 = check_edge(q3,q4,edge,tmp,n4) 

cell = 0
if ((pre_n3==0).and.(pre_n4==0)) then ! all new nodes 
   ! cell(1,:) = (/ n1, n4, n3/)
   ! cell(2,:) = (/ n1, n3, n2/)
   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, n4, n3/)
      cell(2,:) = (/ n1, n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, n1, n4/)
   !    cell(2,:) = (/ n2, n4, n3/)
   ! endif

! add all nodes to edge 
   id = find_empty_position(edge(q3,q4,:))

   edge(q3,q4,id) = n4
   edge(q4,q3,id) = n4

   edge(q3,q4,id+1) = n3
   edge(q4,q3,id+1) = n3
  
elseif ((pre_n3/=0).and.(pre_n4/=0)) then!both nodes are already logged 
!   cell(1,:) = (/ n1, pre_n4, pre_n3/)
!   cell(2,:) = (/ n1, pre_n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
   cell(1,:) = (/ n1, pre_n4, pre_n3/)
   cell(2,:) = (/ n1, pre_n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, n1, pre_n4/)
   !    cell(2,:) = (/ n2, pre_n4, pre_n3/)
   ! endif


elseif ((pre_n3 /=0).and.(pre_n4 == 0)) then ! node at n3 is already logged
   ! cell(1,:) = (/ n1, n4, pre_n3/)
   ! cell(2,:) = (/ n1, pre_n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, n4, pre_n3/)
      cell(2,:) = (/ n1, pre_n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, n1, n4/)
   !    cell(2,:) = (/ n2, n4, pre_n3/)
   ! endif

   id = find_empty_position(edge(q3,q4,:))
   edge(q3,q4,id) = n4
   edge(q4,q3,id) = n4

elseif ((pre_n3 == 0).and.(pre_n4 /=0)) then !nodes at n4 is already logged
   ! cell(1,:) = (/ n1, pre_n4, n3/)
   ! cell(2,:) = (/ n1, n3, n2/)

   ! call random_number(num)        ! randomly split of cell
   ! if (num <0.5) then
      cell(1,:) = (/ n1, pre_n4, n3/)
      cell(2,:) = (/ n1, n3, n2/)
   ! else
   !    cell(1,:) = (/ n2, n1, pre_n4/)
   !    cell(2,:) = (/ n2, pre_n4, n3/)
   ! endif
   
   id = find_empty_position(edge(q3,q4,:))
   edge(q3,q4,id) = n3
   edge(q4,q3,id) = n3
endif 

return 
end subroutine right_edge 
!###############################################################################

function find_empty_position(arr) result(id)
integer,dimension(:), intent(in) :: arr
integer :: nn,i,id

id = -1 
nn = size(arr)
do i = 1,nn
   if (arr(i) == 0) then
      id = i 
      return 
   endif 
enddo 

if (id ==-1) then 
   write(*,*) 'error: array is full!!'
   stop 
endif 

return 
end function find_empty_position
!###############################################################################
function find_replacement_node(old_id,oldmesh,newmesh)  result(found_index)
type(mesh),intent(in) :: newmesh,oldmesh 
integer,intent(in) :: old_id
integer :: found_index
integer ::i 

found_index = -1 
do i = 1,newmesh%Nnodes
      if ((abs(oldmesh%px(old_id)-newmesh%px(i)) <= 0.1).and.     &
          (abs(oldmesh%py(old_id)-newmesh%py(i)) <= 0.1).and.     &
          (abs(oldmesh%pz(old_id)-newmesh%pz(i)) <= 0.1)) then
         found_index = i
         exit
      endif 
enddo 

return 
end function find_replacement_node
!###############################################################################
function check_edge(n1,n2,edge_list,tmesh,query_node) result(isfound)
! this function returns true if the node is not in the list 
integer,intent(in) :: n1,n2,query_node
integer, dimension(:,:,:),intent(in) :: edge_list
integer :: i,node_id
integer :: isfound
type(mesh),intent(in) :: tmesh

isfound = 0
do i = 1,size(edge_list,3) 
   node_id = edge_list(n1,n2,i)
   if (node_id == 0 ) then !if the first element in the array is 0 all the elements are zero 
       exit
   elseif (abs(tmesh%px(node_id)-tmesh%px(query_node)) <= 0.1 .and. &
           abs(tmesh%py(node_id)-tmesh%py(query_node)) <= 0.1 .and. &
           abs(tmesh%pz(node_id)-tmesh%pz(query_node)) <= 0.1 ) then 
         isfound = node_id
         exit
   endif
enddo

return
end function check_edge
!###############################################################################
subroutine create_quad_mesh(x,y,z,nquad,accuracy,qmesh,output_level)
integer :: i,j,k,nquad,icount
real(pr),allocatable, dimension(:,:) :: x,y,z
real(pr),intent(in) :: accuracy
logical :: not_in_list
type(mesh) :: tmp
type(mesh),intent(out) :: qmesh 
integer,intent(in) :: output_level


allocate(tmp%px(nquad*4),tmp%py(nquad*4),tmp%pz(nquad*4))
tmp%px = 0.d0
tmp%py = 0.d0
tmp%pz = 0.d0

tmp%px(1) = x(1,1)
tmp%py(1) = y(1,1)
tmp%pz(1) = z(1,1)

! make list of individual nodes 
icount = 1
do i = 1,nquad
    do j = 1,4
        not_in_list = .true.
        do k = 1,icount  !assume if the nodes are close : note accuracy is dependent on whether utm or lat/lon is being used
            if ((abs(tmp%px(k)-x(i,j)) <= accuracy ).and.       &
                (abs(tmp%py(k)-y(i,j)) <= accuracy  ).and.       &
                (abs(tmp%pz(k)-z(i,j)) <= accuracy  ) ) then ! 
                 not_in_list = .false.
                 exit   !exit the loop if the node is found
            endif 
        enddo
        if (not_in_list) then 
            icount = icount+1
            tmp%px(icount) = x(i,j)
            tmp%py(icount) = y(i,j)
            tmp%pz(icount) = z(i,j)
        endif 
    enddo 
enddo

if (output_level >= 2) write(*,*) 'number of individual nodes in quadrangle mesh',icount
qmesh%Nnodes = icount 

allocate(qmesh%px(qmesh%Nnodes),qmesh%py(qmesh%Nnodes),qmesh%pz(qmesh%Nnodes))
qmesh%px = tmp%px(1:icount)
qmesh%py = tmp%py(1:icount)
qmesh%pz = tmp%pz(1:icount)

deallocate(tmp%px,tmp%py,tmp%pz)

!link nodes to quadilateral cells
qmesh%ncells = nquad
allocate(qmesh%cell(nquad,4))

do i = 1,nquad
    do j = 1,4
        qmesh%cell(i,j) = i
        do k = 1,qmesh%Nnodes
            if (x(i,j)==qmesh%px(k).and.     &
                y(i,j)==qmesh%py(k).and.     &
                z(i,j)==qmesh%pz(k)) then 
                    qmesh%cell(i,j) = k
                    exit !exit the loop if the node is found
            endif
        enddo 
    enddo
enddo 

return
end subroutine create_quad_mesh
!###############################################################################
subroutine hysea_output(amesh,slip,rake,out_file)
implicit none 
real(pr), allocatable,dimension(:)  :: slip
real(pr),intent(in) :: rake
type(mesh),intent(in) :: amesh
character(30) :: out_file  ! name of output file
character(40) :: filename  ! name of output file with extension
real(pr), parameter :: kilo = 1000.d0

character(len=*), parameter :: fmt ="(12E14.6)"
! character(len=*), parameter :: fmt ="(I3,12E14.6)"
!character(len=*), parameter :: fmt ="(I3,F8.4,F8.4,F8.0,F8.4,F8.4,F8.0,F8.4,F8.4,F8.0,F8.4,F8.4)"
integer :: i,n1,n2,n3

! output for HySEA:  lon, lat,depth,x3,rake,slip 
write(filename,'(a,a)') trim(adjustl(out_file)),'_hysea.dat'

open(11,file=trim(adjustl(filename)),form = 'formatted')
! lat = px, lon = py
!output is: lon, lat, depth, lon,lat,depth,lon,lat,depth,rake,slip
do i =1,amesh%Ncells
   n1 = amesh%cell(i,1)
   n2 = amesh%cell(i,2)
   n3 = amesh%cell(i,3)

   if (slip(i)> 0.001) then
     write(11,fmt) amesh%py(n1),amesh%px(n1),abs(amesh%pz(n1))/kilo,         &
                 amesh%py(n2),amesh%px(n2),abs(amesh%pz(n2))/kilo,         &
                 amesh%py(n3),amesh%px(n3),abs(amesh%pz(n3))/kilo,         &
                 rake,slip(i)
   endif 
enddo 
close(11)


return
end subroutine hysea_output
!###############################################################################
subroutine write_out_rect(out_mesh,out_file)
!subroutine write_out_rect(out_mesh,central_dist,slip,out_file)

  implicit none
  type(mesh),intent(in) :: out_mesh        ! complete mesh  
  character(30) :: out_file  ! name of output file
  character(34) :: filename  ! name of output file with extension
  character (len=20) :: fnamef
  real(pr), allocatable,dimension(:)  :: central_dist,slip
  
write(filename,'(a,a)') trim(adjustl(out_file)),'_rec.vtk'
open(11,file=trim(adjustl(filename)),form = 'formatted')
call dumpmeshvtk_jup_rect(11,out_mesh)  ! changed to this as it is read by meshio (other doesn't recognise POLYDATA)
! save distance from centre of patch 
! call dumpnodeattributevtk(11,out_mesh,central_dist,'Central_Distance',.true.)
! save slip on cells 
!call dumpcellattributevtk(11,out_mesh,slip,'slip',.true.)

return
end subroutine write_out_rect
!###############################################################################
subroutine dumpmeshvtk_jup_rect(dev,out_mesh)

  type(mesh),intent(in) :: out_mesh
  integer(pin) :: dev

  integer(pin) :: i,j

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a25:)') 'DATASET UNSTRUCTURED_GRID'
  if( pr == 4) then
     write(dev,'(a7,i10,a6)') 'POINTS ',out_mesh%Nnodes,' float'
  else
     write(dev,'(a7,i10,a7)') 'POINTS ',out_mesh%Nnodes,' double'
  endif

  do i=1,out_mesh%Nnodes
     write(dev,*) out_mesh%px(i),out_mesh%py(i),out_mesh%pz(i)
  enddo
  write(dev,'(a6,2i10)') 'CELLS ',out_mesh%Ncells,out_mesh%Ncells*5
  do i=1,out_mesh%Ncells
     write(dev,'(a2,$)') '4 '
     write(dev,*) (out_mesh%cell(i,j)-1,j=1,4)
  enddo

  write(dev,*)
  write(dev,'(a11,i10)') 'CELL_TYPES ',out_mesh%Ncells
  do i = 1,out_mesh%Ncells
       write(dev,*) 7
  enddo

  return
end subroutine dumpmeshvtk_jup_rect
!###############################################################################
subroutine write_out(out_mesh,central_dist,slip,out_file)
  implicit none
  type(mesh),intent(in) :: out_mesh        ! complete mesh  
  character(30) :: out_file  ! name of output file
  character(34) :: filename  ! name of output file with extension
  character (len=20) :: fnamef
  real(pr), allocatable,dimension(:)  :: central_dist,slip
  
write(filename,'(a,a)') trim(adjustl(out_file)),'.vtk'
open(11,file=trim(adjustl(filename)),form = 'formatted')
call dumpmeshvtk_jup(11,out_mesh)  ! changed to this as it is read by meshio (other doesn't recognise POLYDATA)
! save distance from centre of patch 
call dumpnodeattributevtk(11,out_mesh,central_dist,'Central_Distance',.true.)
! save slip on cells 
call dumpcellattributevtk(11,out_mesh,slip,'slip',.true.)

return
end subroutine write_out

!###############################################################################
subroutine dumpmeshvtk_jup(dev,out_mesh)

  type(mesh),intent(in) :: out_mesh
  integer(pin) :: dev

  integer(pin) :: i,j

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a25:)') 'DATASET UNSTRUCTURED_GRID'
  if( pr == 4) then
     write(dev,'(a7,i10,a6)') 'POINTS ',out_mesh%Nnodes,' float'
  else
     write(dev,'(a7,i10,a7)') 'POINTS ',out_mesh%Nnodes,' double'
  endif

  do i=1,out_mesh%Nnodes
     write(dev,*) out_mesh%px(i),out_mesh%py(i),out_mesh%pz(i)
  enddo
  write(dev,'(a6,2i10)') 'CELLS ',out_mesh%Ncells,out_mesh%Ncells*4
  do i=1,out_mesh%Ncells
     write(dev,'(a2,$)') '3 '
     write(dev,*) (out_mesh%cell(i,j)-1,j=1,3)
  enddo

  write(dev,*)
  write(dev,'(a11,i10)') 'CELL_TYPES ',out_mesh%Ncells
  do i = 1,out_mesh%Ncells
       write(dev,*) 5
  enddo

  return
end subroutine dumpmeshvtk_jup
!###############################################################################!###############################################################################
subroutine dumpmeshvtk(dev,out_mesh)

  type(mesh) :: out_mesh
  integer(pin) :: dev

  integer(pin) :: i,j

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  if( pr == 4) then
     write(dev,'(a7,i10,a6)') 'POINTS ',out_mesh%Nnodes,' float'
  else
     write(dev,'(a7,i10,a7)') 'POINTS ',out_mesh%Nnodes,' double'
  endif

  do i=1,out_mesh%Nnodes
     write(dev,*) out_mesh%px(i),out_mesh%py(i),out_mesh%pz(i)
  enddo
  write(dev,'(a9,2i10)') 'POLYGONS ',out_mesh%Ncells,out_mesh%Ncells*4
  do i=1,out_mesh%Ncells
     write(dev,'(a2,$)') '3 '
     write(dev,*) (out_mesh%cell(i,j)-1,j=1,3)
  enddo
  return
end subroutine dumpmeshvtk
!###############################################################################
subroutine idumpnodeattributevtk(dev,out_mesh,field,attname,init)
  type(mesh) :: out_mesh
  integer(pin), dimension(out_mesh%Nnodes) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  ! if (init) write(dev,'(a11,i10)') 'POINT_DATA ',out_mesh%Nnodes
  ! if( pr == 4) then
    write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' int 1'
  ! else
  !    write(dev,'(a8,a,a10)') 'SCALARS ',trim(attname),' double 1'
  ! endif
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,out_mesh%Nnodes
     write(dev,*) field(i)
  enddo
  return
end subroutine idumpnodeattributevtk
!###############################################################################
subroutine rdumpnodeattributevtk(dev,out_mesh,field,attname,init)
  type(mesh) :: out_mesh
  real(pr), dimension(out_mesh%Nnodes) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  if (init) write(dev,'(a11,i10)') 'POINT_DATA ',out_mesh%Nnodes
  if( pr == 4) then
    write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  else
     write(dev,'(a8,a,a10)') 'SCALARS ',trim(attname),' double 1'
  endif
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,out_mesh%Nnodes
     write(dev,*) field(i)
  enddo
  return
end subroutine rdumpnodeattributevtk
!###############################################################################
subroutine dumpcellattributevtk(dev,out_mesh,field,attname,init)
  type(mesh) :: out_mesh
  real(pr), dimension(out_mesh%Ncells) :: field
  integer(pin) :: dev
  character*(*) :: attname
  logical :: init

  integer(pin) :: i

  if (init) write(dev,'(a10,i10)') 'CELL_DATA ',out_mesh%Ncells
  if( pr == 4) then
     write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  else
     write(dev,'(a8,a,a10)') 'SCALARS ',trim(attname),' double 1'
  endif

  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,out_mesh%Ncells
     write(dev,*) field(i)
  enddo
  return
end subroutine dumpcellattributevtk
!###############################################################################

end module LAT_mesh
