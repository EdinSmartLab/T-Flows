!==============================================================================!
  subroutine Grid_Mod_Find_Nodes_Cells(grid)
!------------------------------------------------------------------------------!
!   Cells around each node are needed for Lagrangian particle tracking.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, ln1, ln2, n1, n2, s, n, per   ! counters
  integer              :: max_n_cells       ! max number of cells around nodes
  real                 :: lx, ly, lz, x1, y1, z1, x2, y2, z2
  logical              :: px, py, pz
  integer, allocatable :: tmp_nodes_c(:,:)

! Just for checking, erase this later
! integer :: lc
! real    :: xn, xc, yn, yc, zn, zc, max_del, min_del, del
!==============================================================================!

  ! Allocate memory for node coordinates
  n = grid % n_nodes
  allocate(grid % nodes_n_cells(1:n))
  grid % nodes_n_cells(:) = 0

  ! Take aliases
  px = .false.
  py = .false.
  pz = .false.
  lx = grid % per_x;  if(abs(lx) > TINY) px = .true.
  ly = grid % per_y;  if(abs(ly) > TINY) py = .true.
  lz = grid % per_z;  if(abs(lz) > TINY) pz = .true.

  ! This information is missing :-(
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      grid % cells_n_nodes(c2) = grid % faces_n_nodes(s)
      do n = 1, grid % cells_n_nodes(c2)
        grid % cells_n(n,c2) = grid % faces_n(n,s)
      end do
    end if
  end do

  !--------------------------------------------------------!
  !   Find maximum number of cells surrounding each node   !
  !     (use only inside cells at this point in time)      !
  !--------------------------------------------------------!

  ! All faces
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    do ln1 = 1, grid % cells_n_nodes(c1)    ! local node number 1
      n1 = grid % cells_n(ln1, c1)          ! global node number 1
      x1 = grid % xn(n1)
      y1 = grid % yn(n1)
      z1 = grid % zn(n1)
      do ln2 = 1, grid % cells_n_nodes(c2)  ! local node number 2
        n2 = grid % cells_n(ln2, c2)        ! global node number 2
        x2 = grid % xn(n2)
        y2 = grid % yn(n2)
        z2 = grid % zn(n2)
        if(Math_Mod_Distance_Squared(x1, y1, z1, x2, y2, z2) < PICO) then
          grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 2
          grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 2
        end if
        if( px .neqv. py .neqv. pz ) then  ! only one periodic direction
          if( ( px .and.                                                      &
                (Math_Mod_Distance_Squared(x1,y1,z1,x2+lx,y2,z2) < PICO .or.  &
                 Math_Mod_Distance_Squared(x1,y1,z1,x2-lx,y2,z2) < PICO)      &
              ) .or.                                                          &
              ( py .and.                                                      &
                (Math_Mod_Distance_Squared(x1,y1,z1, x2,y2+ly,z2) < PICO .or. &
                 Math_Mod_Distance_Squared(x1,y1,z1, x2,y2-ly,z2) < PICO)     &
              ) .or.                                                          &
              ( pz .and.                                                      &
                (Math_Mod_Distance_Squared(x1,y1,z1, x2,y2,z2+lz) < PICO .or. &
                 Math_Mod_Distance_Squared(x1,y1,z1, x2,y2,z2-lz) < PICO)     &
              ) ) then
            grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 2
            grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 2
          end if
        end if
        if( px .and. py .and. .not. pz ) then  ! x and y are periodic
          if(Math_Mod_Distance_Squared(x1,y1,z1,x2+lx,y2,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2-lx,y2,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2,y2+ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2,y2-ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2+lx,y2+ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2+lx,y2-ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2-lx,y2+ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2-lx,y2-ly,z2) < PICO) then
            grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 2
            grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 2
          end if
        end if
      end do
    end do
  end do

! ! Inside cells
! do c = 1, grid % n_cells
!   do ln = 1, grid % cells_n_nodes(c)  ! local node number
!     n = grid % cells_n(ln, c)         ! global node number
!
!     ! Increase number of cells surrounding the this node by one
!     grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1
!   end do
! end do
!
! ! Boundary cells
! do s = 1, grid % n_bnd_cells
!   c2 = grid % faces_c(2, s)
!   if(c2 >= 0) then
!     print *, 'PANIC!  Something is very wrong in Find_Nodes_Cells'
!   end if
!   do ln = 1, grid % faces_n_nodes(s)  ! local face number
!     n = grid % faces_n(ln, s)         ! global node number
!
!     ! Increase number of cells surrounding the this node by one
!     grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1
!   end do
! end do

  max_n_cells = maxval(grid % nodes_n_cells)
  write(100, *) 'max_n_cells(1) = ', max_n_cells

  ! Allocate memory for cells surrounding each node
  allocate(tmp_nodes_c(1:max_n_cells, 1:grid % n_nodes))

  !----------------------------------------------------------!
  !   Now you can really store the cells surrounding nodes   !
  !----------------------------------------------------------!
  grid % nodes_n_cells(:) = 0  ! re-initialize the cell count

  ! All faces
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    do ln1 = 1, grid % cells_n_nodes(c1)    ! local node number 1
      n1 = grid % cells_n(ln1, c1)          ! global node number 1
      x1 = grid % xn(n1)
      y1 = grid % yn(n1)
      z1 = grid % zn(n1)

      do ln2 = 1, grid % cells_n_nodes(c2)  ! local node number 1
        n2 = grid % cells_n(ln2, c2)        ! global node number 1
        x2 = grid % xn(n2)
        y2 = grid % yn(n2)
        z2 = grid % zn(n2)

        if(Math_Mod_Distance_Squared(x1, y1, z1, x2, y2, z2) < PICO) then
          grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
          tmp_nodes_c(grid % nodes_n_cells(n1), n1) = c1
          grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
          tmp_nodes_c(grid % nodes_n_cells(n1), n1) = c2
          grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
          tmp_nodes_c(grid % nodes_n_cells(n2), n2) = c1
          grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
          tmp_nodes_c(grid % nodes_n_cells(n2), n2) = c2
        end if

        if( px .neqv. py .neqv. pz ) then  ! only one periodic direction
          if( ( px .and.                                                      &
                (Math_Mod_Distance_Squared(x1,y1,z1,x2+lx,y2,z2) < PICO .or.  &
                 Math_Mod_Distance_Squared(x1,y1,z1,x2-lx,y2,z2) < PICO)      &
              ) .or.                                                          &
              ( py .and.                                                      &
                (Math_Mod_Distance_Squared(x1,y1,z1, x2,y2+ly,z2) < PICO .or. &
                 Math_Mod_Distance_Squared(x1,y1,z1, x2,y2-ly,z2) < PICO)     &
              ) .or.                                                          &
              ( pz .and.                                                      &
                (Math_Mod_Distance_Squared(x1,y1,z1, x2,y2,z2+lz) < PICO .or. &
                 Math_Mod_Distance_Squared(x1,y1,z1, x2,y2,z2-lz) < PICO)     &
              ) ) then
            grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
            tmp_nodes_c(grid % nodes_n_cells(n1), n1) = c1
            grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
            tmp_nodes_c(grid % nodes_n_cells(n1), n1) = c2
            grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
            tmp_nodes_c(grid % nodes_n_cells(n2), n2) = c1
            grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
            tmp_nodes_c(grid % nodes_n_cells(n2), n2) = c2
          end if
        end if
        if( px .and. py .and. .not. pz ) then  ! x and y are periodic
          if(Math_Mod_Distance_Squared(x1,y1,z1,x2+lx,y2,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2-lx,y2,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2,y2+ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2,y2-ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2+lx,y2+ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2+lx,y2-ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2-lx,y2+ly,z2) < PICO .or.  &
             Math_Mod_Distance_Squared(x1,y1,z1,x2-lx,y2-ly,z2) < PICO) then
            grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
            tmp_nodes_c(grid % nodes_n_cells(n1), n1) = c1
            grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
            tmp_nodes_c(grid % nodes_n_cells(n1), n1) = c2
            grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
            tmp_nodes_c(grid % nodes_n_cells(n2), n2) = c1
            grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
            tmp_nodes_c(grid % nodes_n_cells(n2), n2) = c2
          end if
        end if
      end do
    end do
  end do

! ! Inside cells
! do c = 1, grid % n_cells
!   do ln = 1, grid % cells_n_nodes(c)  ! local node number
!     n = grid % cells_n(ln, c)         ! global node number
!
!     ! Increase number of cells surrounding the this node by one ...
!     grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1
!
!     ! ... and store the current cell
!     tmp_nodes_c(grid % nodes_n_cells(n), n) = c
!   end do
! end do
!
! ! Boundary cells
! do s = 1, grid % n_bnd_cells
!   c2 = grid % faces_c(2,s)
!   do ln = 1, grid % faces_n_nodes(s)  ! local face number
!     n = grid % faces_n(ln, s)         ! global node number
!
!     ! Increase number of cells surrounding the this node by one ...
!     grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1
!
!     ! ... and store the current cell
!     tmp_nodes_c(grid % nodes_n_cells(n), n) = c2
!
!     ! Also store boundary face for boundary cell
!     grid % cells_bnd_face(c2) = s
!   end do
! end do

  do n = 1, grid % n_nodes
    call Sort_Mod_Unique_Int( grid % nodes_n_cells(n),  &
                tmp_nodes_c(1:grid % nodes_n_cells(n), n) )
  end do

  max_n_cells = maxval(grid % nodes_n_cells)
  write(100, *) 'max_n_cells(2) = ', max_n_cells

  ! Allocate memory for cells surrounding each node
  allocate(grid % nodes_c(1:max_n_cells, 1:grid % n_nodes))
  grid % nodes_c(:,:) = 0

  do n = 1, grid % n_nodes
    grid % nodes_c(1:grid % nodes_n_cells(n), n) = &
       tmp_nodes_c(1:grid % nodes_n_cells(n), n)
  end do

  do n = 1, grid % n_nodes
    write(100, '(3f12.5, i8, 8i8)')                               &
          grid % xn(n),                                           &
          grid % yn(n),                                           &
          grid % zn(n),                                           &
          grid % nodes_n_cells(n), grid % nodes_c(1:max_n_cells, n)
  end do
! ! Just for checking, erase this later
! max_del = -HUGE
! min_del = +HUGE
!
! do n = 1, grid % n_nodes
!   do lc = 1, grid % nodes_n_cells(n)  ! local cell number
!     c = tmp_nodes_c(lc, n)         ! global cell number
!
!     xn = grid % xn(n)
!     yn = grid % yn(n)
!     zn = grid % zn(n)
!
!     xc = grid % xc(c)
!     yc = grid % yc(c)
!     zc = grid % zc(c)
!
!     del = sqrt( (xn-xc)**2 + (yn-yc)**2 + (zn-zc)**2 )
!
!     min_del = min(del, min_del)
!     max_del = max(del, max_del)
!   end do
! end do
!
! print *, '# Checking: min and max del: ', min_del, max_del
! stop

  end subroutine
