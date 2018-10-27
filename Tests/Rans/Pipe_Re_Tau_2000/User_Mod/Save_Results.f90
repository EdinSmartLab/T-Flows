!==============================================================================!
  subroutine User_Mod_Save_Results(grid, save_name) 
!------------------------------------------------------------------------------!
!   This subroutine reads name.1r file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
!------------------------------------------------------------------------------!
  use Grid_Mod
  use Flow_Mod
  use Rans_Mod
  use Comm_Mod                       ! parallel stuff
  use Name_Mod,  only: problem_name
  use Const_Mod                      ! constants
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: save_name
!-----------------------------------[Locals]-----------------------------------!
  integer             :: n_prob, pl, c, i, count, s, c1, c2, n_points
  character(len=80)   :: coord_name, res_name, res_name_plus
  character(len=80)   :: store_name
  real,allocatable    :: z_p(:),                              &
                         u_p(:), v_p(:), w_p(:), y_plus_p(:), &
                         kin_p(:), eps_p(:), f22_p(:),        & 
                         vis_t_p(:), uw_p(:), zeta_p(:),      &  
                         t_p(:), tt_p(:),                    &   
                         ut_p(:), vt_p(:), wt_p(:), &  
                         ind(:),                    &
                         wall_p(:)
  integer,allocatable :: n_p(:), n_count(:)
  real                :: t_wall, t_tau, d_wall, nu_max, b11, b12, Rad 
  real                :: Ubulk, Error, Re, Cf_Dean, Cf, pr, u_tau_p, R
  logical             :: there
!==============================================================================!

  ! Set the name for coordinate file
  call Name_File(0, coord_name, ".1r")

  ! Store the name
  store_name = problem_name
  problem_name = save_name

  call Name_File(0, res_name,      "-res.dat")
  call Name_File(0, res_name_plus, "-res-plus.dat")

  !------------------!
  !   Read 1r file   !
  !------------------!
  inquire(file=coord_name, exist=there)
  if(.not. there) then
    if(this_proc < 2) then
      print *, '==============================================================='
      print *, 'In order to extract profiles and write them in ascii files'
      print *, 'the code has to read cell-faces coordinates '
      print *, 'in wall-normal direction in the ascii file ''case_name.1r.'''
      print *, 'The file format should be as follows:'
      print *, '10  ! number of cells + 1'
      print *, '1 0.0'
      print *, '2 0.1'
      print *, '3 0.2'
      print *, '... '
      print *, '==============================================================='
    end if

    ! Restore the name and return
    problem_name = store_name
    return
  end if

  Ubulk = bulk % flux_z / (density*bulk % area_z)
  t_wall = 0.0
  nu_max = 0.0
  n_points = 0

  open(9, file=coord_name)

  ! Write the number of searching intervals
  read(9,*) n_prob
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! Read the intervals positions
  do pl=1,n_prob
    read(9,*) ind(pl), z_p(pl)
  end do
  close(9)

  allocate(n_p     (n_prob));  n_p      = 0
  allocate(wall_p  (n_prob));  wall_p   = 0.0
  allocate(u_p    (n_prob));  u_p     = 0.0
  allocate(v_p    (n_prob));  v_p     = 0.0
  allocate(w_p    (n_prob));  w_p     = 0.0
  allocate(kin_p   (n_prob));  kin_p    = 0.0
  allocate(eps_p   (n_prob));  eps_p    = 0.0
  allocate(uw_p    (n_prob));  uw_p     = 0.0
  allocate(vis_t_p (n_prob));  vis_t_p  = 0.0
  allocate(f22_p   (n_prob));  f22_p    = 0.0
  allocate(zeta_p  (n_prob));  zeta_p   = 0.0
  allocate(y_plus_p(n_prob));  y_plus_p = 0.0

  allocate(n_count(n_prob)); n_count=0
  count = 0
  if(heat_transfer) then
    allocate(t_p(n_prob));  t_p = 0.0
    allocate(tt_p(n_prob));  tt_p = 0.0
    allocate(ut_p(n_prob));  ut_p = 0.0
    allocate(vt_p(n_prob));  vt_p = 0.0
    allocate(wt_p(n_prob));  wt_p = 0.0
  end if

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c = 1, grid % n_cells - grid % comm % n_buff_cells 
      Rad = 1.0 - grid % wall_dist(c)
      if( Rad < (z_p(i)) .and.  &
          Rad > (z_p(i+1))) then
        R = sqrt(grid % xc(c)**2 + grid % yc(c)**2)
        b11 = grid % xc(c)/R
        b12 = grid % yc(c)/R

        wall_p(i) = wall_p(i) + grid % wall_dist(c)
        u_p(i)   = u_p(i) + u % n(c)
        v_p(i)   = v_p(i) + v % n(c)
        w_p(i)   = w_p(i) + w % n(c)

        kin_p(i)   = kin_p(i) + kin % n(c)
        eps_p(i)   = eps_p(i) + eps % n(c)
        uw_p(i)    = uw_p(i) &
                   + b11 * vis_t(c) *(u % z(c) + w % x(c)) &
                   + b12 * vis_t(c) *(v % z(c) + w % y(c)) 
        vis_t_p(i) = vis_t_p(i) + vis_t(c) / viscosity
        y_plus_p(i)= y_plus_p(i) + y_plus(c)

        if(turbulence_model .eq. K_EPS_ZETA_F) then
          f22_p(i)  = f22_p(i) + f22  % n(c)
          zeta_p(i) = zeta_p(i) + zeta % n(c)
        end if

        if(heat_transfer) then
          t_p(i)   = t_p(i) + t % n(c)
          ut_p(i)   = ut_p(i) + ut % n(c)
          vt_p(i)   = vt_p(i) + vt % n(c)
          wt_p(i)   = wt_p(i) + wt % n(c)
        end if
        n_count(i) = n_count(i) + 1
      end if 
    end do
  end do

  ! Average over all processors
  do pl=1, n_prob-1
    call Comm_Mod_Global_Sum_Int(n_count(pl))

    call Comm_Mod_Global_Sum_Real(wall_p(pl))

    call Comm_Mod_Global_Sum_Real(u_p(pl))
    call Comm_Mod_Global_Sum_Real(v_p(pl))
    call Comm_Mod_Global_Sum_Real(w_p(pl))

    call Comm_Mod_Global_Sum_Real(kin_p(pl))
    call Comm_Mod_Global_Sum_Real(eps_p(pl))
    call Comm_Mod_Global_Sum_Real(uw_p(pl))
    call Comm_Mod_Global_Sum_Real(vis_t_p(pl))
    call Comm_Mod_Global_Sum_Real(y_plus_p(pl))

    call Comm_Mod_Global_Sum_Real(f22_p(pl))
    call Comm_Mod_Global_Sum_Real(zeta_p(pl))

    count =  count + n_count(pl)

    if(heat_transfer) then
      call Comm_Mod_Global_Sum_Real(t_p(pl))
      call Comm_Mod_Global_Sum_Real(tt_p(pl))
      call Comm_Mod_Global_Sum_Real(ut_p(pl))
      call Comm_Mod_Global_Sum_Real(vt_p(pl))
      call Comm_Mod_Global_Sum_Real(wt_p(pl))
    end if
  end do


  call Comm_Mod_Wait

  do i = 1, n_prob-1
    if(n_count(i) .ne. 0) then
      wall_p(i)  = wall_p(i)/n_count(i)
      u_p(i)    = u_p(i)/n_count(i)
      v_p(i)    = v_p(i)/n_count(i)
      w_p(i)    = w_p(i)/n_count(i)

      kin_p(i)   = kin_p(i)/n_count(i)
      eps_p(i)   = eps_p(i)/n_count(i)
      uw_p(i)    = uw_p(i)/n_count(i)
      vis_t_p(i) = vis_t_p(i)/n_count(i)
      f22_p(i)   = f22_p(i)/n_count(i)
      zeta_p(i)  = zeta_p(i)/n_count(i)
      y_plus_p(i)= y_plus_p(i)/n_count(i)
      if(heat_transfer) then
        t_p(i)  = t_p(i)/n_count(i)
        tt_p(i) = tt_p(i)/n_count(i)
        ut_p(i) = ut_p(i)/n_count(i)
        vt_p(i) = vt_p(i)/n_count(i)
        wt_p(i) = wt_p(i)/n_count(i)
      end if
    end if
  end do

  ! Calculating friction velocity and friction temperature
  if(y_plus_p(1) > 5.0) then
    u_tau_p = sqrt(max(abs(bulk % p_drop_x),  &
                       abs(bulk % p_drop_y),  &
                       abs(bulk % p_drop_z))/density)
  else  
    u_tau_p =  sqrt( (viscosity*sqrt(u_p(1)**2 +        &
                                     v_p(1)**2 +        &
                                     w_p(1)**2)         &
                                     / wall_p(1))       &
                                     / density)
  end if

  if(u_tau_p .eq. 0.0) then
    if(this_proc < 2) then
      write(*,*) '# Friction velocity is zero in Save_Results.f90!'
    end if

    ! Restore the name and return
    problem_name = store_name
    return
  end if

  if(heat_transfer) then 
    d_wall = 0.0
    do c = 1, grid % n_cells
      if(grid % wall_dist(c) > d_wall) then
        d_wall = grid % wall_dist(c)
        t_inf  = t % n(c)
      end if
    end do

    call Comm_Mod_Wait

    if(heat_flux> 0.0) then
      call Comm_Mod_Global_Min_Real(t_inf)
    else
      call Comm_Mod_Global_Max_Real(t_inf)
    end if
    
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2  < 0) then
        if( Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. WALL .or.  &
            Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. WALLFL) then

          t_wall = t_wall + t % n(c2)
          nu_max = nu_max + t % q(c2)/(conductivity*(t % n(c2) - t_inf))
          n_points = n_points + 1
        end if
      end if
    end do

    call Comm_Mod_Global_Sum_Real(t_wall)
    call Comm_Mod_Global_Sum_Real(nu_max)

    call Comm_Mod_Wait

    t_wall = t_wall / n_points
    nu_max = nu_max / n_points
    t_tau  = heat_flux / (density * capacity * u_tau_p)
  end if

  open(3, file = res_name)
  open(4, file = res_name_plus)
  
  do i = 3, 4
    pr = viscosity * capacity / conductivity
    Re = density * Ubulk * 2.0/viscosity
    Cf_Dean = 0.0791*(Re)**(-0.25)
    Cf      = u_tau_p**2/(0.5*Ubulk**2)
    Error   = abs(Cf_Dean - Cf)/Cf_Dean * 100.0
    write(i,'(A1,(A12,E12.6))')  &
    '#', 'Ubulk    = ', Ubulk 
    write(i,'(A1,(A12,E12.6))')  &
    '#', 'Re       = ', density * Ubulk * 2.0/viscosity
    write(i,'(A1,(A12,E12.6))')  &
    '#', 'Re_tau   = ', density*u_tau_p/viscosity
    write(i,'(A1,(A12,E12.6))')  &
    '#', 'Cf       = ', 2.0*(u_tau_p/Ubulk)**2
    write(i,'(A1,(A12,F12.6))')  &
    '#', 'Utau     = ', u_tau_p 
    write(i,'(A1,(A12,F12.6,A2,A22))') & 
    '#', 'Cf_error = ', Error, ' %', 'Dean formula is used.'
    if(heat_transfer) then
      write(i,'(A1,(A12, F12.6))')'#', 'Nu number =', nu_max 
      write(i,'(A1,(A12, F12.6,A2,A39))')'#', 'Nu_error  =', &
      abs(0.023*0.5*Re**0.8*pr**0.4 - & 
      nu_max)/(0.023*0.5*Re**0.8*pr**0.4) * 100.0, ' %',&
      'correlation of Dittus-Boelter is used.' 
    end if

    if(turbulence_model .eq. K_EPS) then
      if(heat_transfer) then
        write(i,'(A1,2X,A60)') '#',  ' r,'                    //  &
                                     ' w,'                    //  &
                                     ' kin, eps, uw,'         //  &
                                     ' vis_t/viscosity,'      //  &
                                     ' t, ut, vt, wt,'   
      else
        write(i,'(A1,2X,A60)')  '#', ' r,'                    //  &
                                     ' w,'                    //  &
                                     ' kin, eps, uw, vis_t/viscosity'
      end if
    else if(turbulence_model .eq. K_EPS_ZETA_F) then
      if(heat_transfer) then
        write(i,'(A1,2X,A60)') '#',  ' r,'                    //  &
                                     ' w,'                    //  &
                                     ' kin, eps, uw,'         //  &
                                     ' f22, zeta,'            //  &
                                     ' vis_t/viscosity,'      //  &
                                     ' t, ut, vt, wt'
      else
        write(i,'(A1,2X,A50)') '#', ' r,'                     //  &
                                    ' w,'                     //  &
                                    ' kin, eps, uw,'          //  &
                                    ' f22, zeta'              //  & 
                                    ' vis_t/viscosity,'       
      end if
    end if
  end do

  if(heat_transfer) then
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(3,'(12e15.7,I5)') wall_p(i),                       &
                             w_p(i),                          &
                             kin_p(i), eps_p(i), uw_p(i),     &
                             f22_p(i), zeta_p(i), vis_t_p(i), &
                             t_p(i), ut_p(i), vt_p(i), wt_p(i),&
                             n_count(i)
      end if
    end do
  else
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(3,'(8e15.7,I5)')  wall_p(i),                        &
                             w_p(i),                          &
                             kin_p(i), eps_p(i), uw_p(i),     &
                             f22_p(i), zeta_p(i), vis_t_p(i), &
                             n_count(i)
      end if
    end do
  end if

  close(3)

  do i = 1, n_prob-1
    wall_p(i)= density * wall_p(i)*u_tau_p/viscosity
    w_p(i) = w_p(i)/u_tau_p
    kin_p(i) = kin_p(i)/u_tau_p**2                      ! kin%n(c)
    eps_p(i) = eps_p(i)*viscosity/(u_tau_p**4.0*density)! eps%n(c)
    uw_p(i) = uw_p(i)/(u_tau_p**2*density)            ! vis_t(c)*(u%z(c)+w%x(c))

    if(turbulence_model .eq. K_EPS_ZETA_F) then
      f22_p(i) = f22_p(i)*viscosity/u_tau_p**2.0   ! f22%n(c)
    end if
 
    if(heat_transfer) then
      t_p(i) = (t_wall - t_p(i))/t_tau   ! t % n(c)
      ut_p(i) = ut_p(i)/(u_tau_p*t_tau)  ! ut % n(c)
      vt_p(i) = vt_p(i)/(u_tau_p*t_tau)  ! vt % n(c)
      wt_p(i) = wt_p(i)/(u_tau_p*t_tau)  ! wt % n(c)
    end if
  end do

  if(heat_transfer) then
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(4,'(12e15.7)') wall_p(i),                       &
                             w_p(i),                          &
                             kin_p(i), eps_p(i), uw_p(i),     &
                             f22_p(i), zeta_p(i), vis_t_p(i), &
                             t_p(i), ut_p(i), vt_p(i), wt_p(i)
      end if
    end do
  else
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(4,'(8e15.7)')  wall_p(i),                        &
                             w_p(i),                          &
                             kin_p(i), eps_p(i), uw_p(i),     &
                             f22_p(i), zeta_p(i), vis_t_p(i)
      end if
    end do
  end if

  close(4)

  deallocate(n_p)
  deallocate(z_p)
  deallocate(u_p)
  deallocate(v_p)
  deallocate(w_p)
  deallocate(kin_p)
  deallocate(eps_p)
  deallocate(uw_p)
  deallocate(vis_t_p)
  deallocate(f22_p)
  deallocate(zeta_p)
  deallocate(y_plus_p)
  if(heat_transfer) then
    deallocate(t_p)
    deallocate(tt_p)
    deallocate(ut_p)
    deallocate(vt_p)
    deallocate(wt_p)
  end if

  if(this_proc < 2)  write(6, *) '# Finished with User_Mod_Save_Results.f90.'

  ! Restore the name
  problem_name = store_name

  end subroutine
