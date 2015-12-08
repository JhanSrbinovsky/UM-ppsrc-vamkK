! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine Set_var_grid --------------------------------------------
!
! Purpose : to derive variable grid parameters at each point
!           lambda's and phi's no-longer have fixed intervals
!
! Method:
!    1. See paper
!                Setting and testing variable resolution grids
!                   inside LAM domains.              Terry Davies
!    2. Needs to be called separately for each direction
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Top Level

      SUBROUTINE Set_var_grid (                                         &
     &                        Lambda_p_in, Lambda_u_in,                 &
     &                        Phi_p_in, Phi_v_in,                       &
     &                        global_row_length, global_rows,           &
     &                        row_length, rows,  n_rows,                &
     &                        halo_i, halo_j, L_varres,                 &
     &                        delta_lambda_in, delta_phi_in,            &
     &                        base_lambda_in, base_phi_in,              &
     &                        lam_var, phi_var,                         &
     &                        var_ratio, lam_ratio, phi_ratio,          &
     &                        lam_frac, phi_frac,                       &
     &                        glambda_p, glambda_u, phi_p2d, phi_v2d,   &
     &                        gdlambda_p, gdlambda_u,                   &
     &                        dphi_p2d, dphi_v2d,                       &
     &                        grecip_dlamp, grecip_dlamu,               &
     &                        recip_dphip2d, recip_dphiv2d,             &
     &                        wt_lambda_p, wt_lambda_u,                 &
     &                        wt_phi_p2d, wt_phi_v2d,                   &
     &                        lambda_p_rm, lambda_p_rp,                 &
     &                        lambda_u_rm, lambda_u_rp,                 &
     &                        phi_p_rm2d, phi_p_rp2d,                   &
     &                        phi_v_rm2d, phi_v_rp2d,                   &
     &                        lambda_p_end, lambda_u_end,               &
     &                        phi_p_end, phi_v_end, dlambda_p_end,      &
     &                        dlambda_u_end, dphi_p_end, dphi_v_end,    &
     &                        recip_lambda_p_m, recip_lambda_p,         &
     &                        recip_lambda_p_p, recip_lambda_p_p2,      &
     &                        recip_lambda_u_m, recip_lambda_u,         &
     &                        recip_lambda_u_p, recip_lambda_u_p2,      &
     &                        recip_phi_p_m2d, recip_phi_p2d,           &
     &                        recip_phi_p_p2d, recip_phi_p_p22d,        &
     &                        recip_phi_v_m2d, recip_phi_v2d,           &
     &                        recip_phi_v_p2d, recip_phi_v_p22d,        &
     &                        max_pts, recip_dlam, recip_dphi,          &
     &                        halo_lam, halo_phi, look_lam, look_phi,   &
     &                        model_domain, datastart )


      USE conversions_mod, ONLY: pi_over_180

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE PrintStatus_mod
      IMPLICIT NONE


      INTEGER                                                           &
     &  row_length                                                      &
                    ! number of points on a row this pe
     &, rows                                                            &
                    ! number of points in a column this pe
     &, n_rows                                                          &
                    ! number of points in a column this pe
     &, halo_i                                                          &
                    ! NS halo size
     &, halo_j                                                          &
                    ! EW halo size
     &, global_row_length                                               &
                    ! global number of points on a row
     &, global_rows                                                     &
                    ! global number of points in a column
     &, lam_var                                                         &
                    ! number of variable res. lambda intervals
     &, phi_var                                                         &
                    ! number of variable res. phi intervals
     &, datastart(3)                                                    &
                    ! start pointers for this pe
     &, model_domain                                                    &
                    ! holds integer code for model domain
     &, max_pts
                    ! max size of look-up array

      REAL                                                              &
     &  Delta_lambda_in                                                 &
                          ! EW (x) grid spacing in degrees
     &, Delta_phi_in                                                    &
                          ! NS (y)grid spacing in degrees
     &, Base_lambda_in                                                  &
                          ! Longitude of first theta point in degs
     &, Base_phi_in                                                     &
                          ! Latitude of first theta point in degrees
     &, lambda_p_end                                                    &
     &, phi_p_end                                                       &
     &, lambda_u_end                                                    &
     &, phi_v_end                                                       &
     &, dlambda_p_end                                                   &
     &, dlambda_u_end                                                   &
     &, dphi_p_end                                                      &
     &, dphi_v_end                                                      &
     &, var_ratio                                                       &
                   ! grid-stretcing ratio for var grid
     &, lam_ratio                                                       &
                   ! scaling of original grid to high res grid
     &, phi_ratio                                                       &
                   ! scaling of original grid to high res grid
     &, lam_frac                                                        &
                   ! proportion of reg. points in West of domain
     &, phi_frac  
                   ! proportion of reg. points in South of domain

! Input VarRes grid info in degrees
      REAL                                                              &
     &  Lambda_p_in(global_row_length)                                  &
                                      !IN EW and NS VarRes gri
     &, Lambda_u_in(global_row_length)                                  &
                                      !IN EW and NS u,v grid
     &, Phi_p_in(global_rows)                                           &
                                      !location in degrees
     &, Phi_v_in(global_rows)       
                                      !location in degrees

      LOGICAL L_varres
! Output 

! VarRes Array co-ordinates in radians
      REAL:: glambda_p         ( 1-halo_i : global_row_length+halo_i )
      REAL:: glambda_u         ( 1-halo_i : global_row_length+halo_i )
      REAL:: Phi_p2d           ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : rows + halo_j )
      REAL:: Phi_v2d           ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : n_rows+halo_j )
      REAL:: gdlambda_p         ( 1-halo_i : global_row_length+halo_i )
      REAL:: gdlambda_u         ( 1-halo_i : global_row_length+halo_i )
      REAL:: dPhi_p2d          ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : rows+halo_j )
      REAL:: dPhi_v2d          ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : n_rows+halo_j )
      REAL:: grecip_dlamp      ( 1-halo_i : global_row_length+halo_i )
      REAL:: grecip_dlamu      ( 1-halo_i : global_row_length+halo_i )
      REAL:: recip_dPhip2d     ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : rows+halo_j )
      REAL:: recip_dPhiv2d     ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : n_rows+halo_j )
      REAL:: wt_lambda_p         ( 1-halo_i : row_length+halo_i )
      REAL:: wt_lambda_u         ( 1-halo_i : row_length+halo_i )
      REAL:: wt_Phi_p2d          ( 1-halo_i : row_length + halo_i,      &
     &                             1-halo_j : rows+halo_j )
      REAL:: wt_Phi_v2d          ( 1-halo_i : row_length + halo_i,      &
     &                             1-halo_j : n_rows+halo_j )
      REAL:: lambda_p_rm       ( 1-halo_i : row_length+halo_i )
      REAL:: lambda_p_rp       ( 1-halo_i : row_length+halo_i )
      REAL:: lambda_u_rm       ( 1-halo_i : row_length+halo_i )
      REAL:: lambda_u_rp       ( 1-halo_i : row_length+halo_i )
      REAL:: phi_p_rm2d        ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : rows+halo_j )
      REAL:: phi_p_rp2d        ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : rows+halo_j )
      REAL:: phi_v_rm2d        ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : n_rows+halo_j )
      REAL:: phi_v_rp2d        ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : n_rows+halo_j )
      REAL:: Recip_lambda_p_m  ( 1-halo_i : row_length+halo_i )
      REAL:: Recip_lambda_p    ( 1-halo_i : row_length+halo_i )
      REAL:: Recip_lambda_p_p  ( 1-halo_i : row_length+halo_i )
      REAL:: Recip_lambda_p_p2 ( 1-halo_i : row_length+halo_i )
      REAL:: Recip_lambda_u_m  ( 1-halo_i : row_length+halo_i )
      REAL:: Recip_lambda_u    ( 1-halo_i : row_length+halo_i )
      REAL:: Recip_lambda_u_p  ( 1-halo_i : row_length+halo_i )
      REAL:: Recip_lambda_u_p2 ( 1-halo_i : row_length+halo_i )
      REAL:: Recip_phi_p_m2d   ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : rows+halo_j )
      REAL:: Recip_phi_p2d    ( 1-halo_i : row_length + halo_i,         &
     &                           1-halo_j : rows+halo_j )
      REAL:: Recip_phi_p_p2d   ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : rows+halo_j )
      REAL:: Recip_phi_p_p22d  ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : rows+halo_j )
      REAL:: Recip_phi_v_m2d   ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : n_rows+halo_j )
      REAL:: Recip_phi_v2d    ( 1-halo_i : row_length + halo_i,         &
     &                           1-halo_j : n_rows+halo_j )
      REAL:: Recip_phi_v_p2d   ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : n_rows+halo_j )
      REAL:: Recip_phi_v_p22d  ( 1-halo_i : row_length + halo_i,        &
     &                           1-halo_j : n_rows+halo_j )

      REAL ::  recip_dlam   ! smallest delta_lambda p
      REAL ::  recip_dphi   ! smallest delta_phi p

      INTEGER :: halo_lam  ! halo for lamp look-up table
      INTEGER :: halo_phi  ! halo for phip look-up table

! Set of look-up tables for searching on variable grids
      INTEGER :: look_lam(max_pts)  ! lambda p search
      INTEGER :: look_phi(max_pts)  ! phi p search
         
! local  variables       

      REAL                                                              &
     &  delta_p, delta_u, delta_v

! loop counters
      INTEGER                                                           &
     &  i, j                                                            &
     &, info                                                            &
     &, gi, gj 

! local arrays       
      REAL                                                              &
     &  wk_lambda_p ( 1-halo_i : global_row_length + halo_i )           &
     &, wk_lambda_u ( 1-halo_i : global_row_length + halo_i )           &
     &, wk_dlambda_p( 1-halo_i : global_row_length + halo_i )           &
     &, wk_dlambda_u( 1-halo_i : global_row_length + halo_i )           &
     &, wk_phi_p  ( 1-halo_j : global_rows + halo_j )                   &
     &, wk_phi_v  ( 1-halo_j : global_rows + halo_j )                   &
     &, wk_dphi_p ( 1-halo_j : global_rows + halo_j )                   &
     &, wk_dphi_v ( 1-halo_j : global_rows + halo_j )                   &
     &, lambda_p ( 1-halo_i : row_length+halo_i )                       &
     &, lambda_u ( 1-halo_i : row_length+halo_i )                       &
     &, dlambda_p ( 1-halo_i : row_length+halo_i )                      &
     &, dlambda_u ( 1-halo_i : row_length+halo_i )

      REAL:: Phi_p ( 1-halo_j : rows+halo_j )
      REAL:: Phi_v( 1-halo_j : n_rows+halo_j )
      REAL:: dPhi_p ( 1-halo_j : rows+halo_j )
      REAL:: dPhi_v( 1-halo_j : n_rows+halo_j )
      REAL:: recip_dPhip ( 1-halo_j : rows+halo_j )
      REAL:: recip_dPhiv ( 1-halo_j : n_rows+halo_j )
      REAL:: phi_p_rm ( 1-halo_j : rows+halo_j )
      REAL:: phi_p_rp ( 1-halo_j : rows+halo_j )
      REAL:: phi_v_rm ( 1-halo_j : n_rows+halo_j )
      REAL:: phi_v_rp ( 1-halo_j : n_rows+halo_j )
      REAL:: Recip_phi_p_m ( 1-halo_j : rows+halo_j )
      REAL:: Recip_phi_p ( 1-halo_j : rows+halo_j )
      REAL:: Recip_phi_p_p ( 1-halo_j : rows+halo_j )
      REAL:: Recip_phi_p_p2( 1-halo_j : rows+halo_j )
      REAL:: Recip_phi_v_m ( 1-halo_j : n_rows+halo_j )
      REAL:: Recip_phi_v ( 1-halo_j : n_rows+halo_j )
      REAL:: Recip_phi_v_p ( 1-halo_j : n_rows+halo_j )
      REAL:: Recip_phi_v_p2( 1-halo_j : n_rows+halo_j )

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('SET_VAR_GRID',zhook_in,zhook_handle)
      If( .not. L_varres ) Then !  compute variable grid parameters
        
        If (PrintStatus >= PrStatus_Normal ) Then
          write(6,*)
          write(6,*) ' calling calc_var_grid for lambda and phi '
        End If  

! DEPENDS ON: calc_var_grid
        Call calc_var_grid (                                            &
     &                   lambda_p_in, lambda_u_in, model_domain, .true.,&
     &                   .true., global_row_length, datastart, 1,       &
     &                   global_row_length, global_row_length, lam_var, &
     &                   delta_lambda_in, base_lambda_in,               &
     &                   lam_ratio, var_ratio, lam_frac )

! DEPENDS ON: calc_var_grid
        Call calc_var_grid (                                            &
     &                   phi_p_in, phi_v_in, model_domain, .false.,     &
     &                   .true., global_rows, datastart, 2,             &
     &                   global_rows, global_rows - 1, phi_var,         &
     &                   delta_phi_in, base_phi_in,                     &
     &                   phi_ratio, var_ratio, phi_frac )

      End If

!      write(*,*) 'lambda_p_in =', lambda_p_in
!      write(*,*) 'lambda_u_in =', lambda_u_in
!      write(*,*) 'phi_p_in =', phi_p_in
!      write(*,*) 'phi_v_in =', phi_v_in
       
! Initialise wk arrays
      do i = 1-halo_i, global_row_length + halo_i
        wk_lambda_p(i) = 0.0
        wk_lambda_u(i) = 0.0
      end do
      do j = 1-halo_j, global_rows + halo_j
        wk_phi_p(j) = 0.0
        wk_phi_v(j) = 0.0
      end do

      do i = 1, global_row_length  
        wk_lambda_p(i) = lambda_p_in(i) 
        wk_lambda_u(i) = lambda_u_in(i)  
      end do 
      
      do i = 1, global_row_length - 1 
        wk_dlambda_p(i) = lambda_p_in(i+1) -  lambda_p_in(i)
        wk_dlambda_u(i) = lambda_u_in(i+1) -  lambda_u_in(i)  
      end do
      wk_dlambda_p(global_row_length)=wk_dlambda_p(global_row_length-1)
      wk_dlambda_u(global_row_length)=wk_dlambda_p(global_row_length-1)
      
      do j = 1, global_rows 
        wk_phi_p(j) = phi_p_in(j)  
      end do 
      do j = 1, global_rows -1
        wk_phi_v(j) = phi_v_in(j)
      end do 
      wk_phi_v(global_rows) = 2 * phi_v_in(global_rows-1) -             &
     &                            phi_v_in(global_rows-2)
      
      do j = 1, global_rows - 1
        wk_dphi_p(j) = phi_p_in(j+1) - phi_p_in(j) 
      end do
      do j = 1, global_rows - 2
        wk_dphi_v(j) = phi_v_in(j+1) - phi_v_in(j) 
      end do
      wk_dphi_p(global_rows ) =  wk_dphi_p(global_rows -1) 
      wk_dphi_v(global_rows-1 ) =  wk_dphi_v(global_rows -2)  
      wk_dphi_v(global_rows ) =  wk_dphi_v(global_rows -1)
         
!==========================
! fill in the halos for lambda and and phi
!==========================      
      delta_p = wk_dlambda_p(1)
      delta_u = wk_dlambda_u(1)
      do i =  0, 1 - halo_i, -1
        wk_lambda_p(i) = wk_lambda_p(i+1) - delta_p
        wk_lambda_u(i) = wk_lambda_u(i+1) - delta_u
        wk_dlambda_p(i) = delta_p
        wk_dlambda_u(i) = delta_u
      end do
      delta_p = wk_dlambda_p(global_row_length)
      delta_u = wk_dlambda_u(global_row_length)
      do i =  global_row_length + 1, global_row_length + halo_i
        wk_lambda_p(i) = wk_lambda_p(i-1) + delta_p
        wk_lambda_u(i) = wk_lambda_u(i-1) + delta_u
        wk_dlambda_p(i) = delta_p
        wk_dlambda_u(i) = delta_u
      end do

      delta_p = wk_dphi_p(1)
      delta_v = wk_dphi_v(1)
      do j =  0, 1 - halo_j, -1
        wk_phi_p(j) = wk_phi_p(j+1) - delta_p
        wk_phi_v(j) = wk_phi_v(j+1) - delta_v
        wk_dphi_p(j) = delta_p
        wk_dphi_v(j) = delta_v
      end do
      delta_p = wk_dphi_p(global_rows)
      delta_v = wk_dphi_v(global_rows)
      do j =  global_rows + 1, global_rows + halo_j
        wk_phi_p(j) = wk_phi_p(j-1) + delta_p
        wk_phi_v(j) = wk_phi_v(j-1) + delta_v
        wk_dphi_p(j) = delta_p
        wk_dphi_v(j) = delta_v
      end do
      
!==========================
! delta_lambda / phi and lambda /phi  
!==========================  

      do i = 1-halo_i , global_row_length + halo_i
        wk_lambda_p(i) = wk_lambda_p(i) * Pi_over_180
        wk_lambda_u(i) = wk_lambda_u(i) * Pi_over_180
        wk_dlambda_p(i) = wk_dlambda_p(i) * Pi_over_180
        wk_dlambda_u(i) = wk_dlambda_u(i) * Pi_over_180
      end do
      do j = 1-halo_j , global_rows + halo_j
        wk_phi_p(j) = wk_phi_p(j) * Pi_over_180
        wk_phi_v(j) = wk_phi_v(j) * Pi_over_180
        wk_dphi_p(j) = wk_dphi_p(j) * Pi_over_180
        wk_dphi_v(j) = wk_dphi_v(j) * Pi_over_180
      end do

      lambda_p_end = wk_lambda_p( global_row_length )
      lambda_u_end = wk_lambda_u( global_row_length )
      phi_p_end = wk_phi_p( global_rows )
      phi_v_end = wk_phi_v( global_rows-1 )
      dlambda_p_end = wk_dlambda_p( global_row_length )
      dlambda_u_end = wk_dlambda_u( global_row_length )
      dphi_p_end = wk_dphi_p( global_rows )
      dphi_v_end = wk_dphi_v( global_rows-1 )
 
!==========================
! delta_lambda / phi and lambda /phi  
!==========================  
     
      do gi = 1-halo_i , global_row_length + halo_i 
        glambda_p(gi) =  wk_lambda_p(gi)
        glambda_u(gi) =  wk_lambda_u(gi)
        gdlambda_p(gi) =  wk_dlambda_p(gi)
        gdlambda_u(gi) =  wk_dlambda_u(gi)
        grecip_dlamp(gi) =  1.0 / wk_dlambda_p(gi)
        grecip_dlamu(gi) =  1.0 / wk_dlambda_u(gi)
      end do
      do i = 1-halo_i , row_length + halo_i 
        gi = datastart(1) + i - 1
        lambda_p(i) =  wk_lambda_p(gi)
        lambda_u(i) =  wk_lambda_u(gi)
        dlambda_p(i) =  wk_dlambda_p(gi)
        dlambda_u(i) =  wk_dlambda_u(gi)
      end do
      
      do j = 1-halo_j , rows + halo_j 
        gj = datastart(2) + j - 1
          phi_p(j) =  wk_phi_p(gj)
          dphi_p(j) =  wk_dphi_p(gj)
          recip_dphip(j) =  1.0 / dphi_p(j)
      end do 
 
      do j = 1-halo_j , n_rows + halo_j 
        gj = datastart(2) + j - 1
          phi_v(j) =  wk_phi_v(gj)
          dphi_v(j) =  wk_dphi_v(gj)
          recip_dphiv(j) =  1.0 / dphi_v(j)
      end do 

!  Set denominators for cubic-Lagrange interpolation

      If (PrintStatus >= PrStatus_Normal ) Then
        write(6,*)
        write(6,*) ' calling Set_coeff_lagrange for lambda_p/u'
      End If 
 
! DEPENDS ON: Set_coeff_lagrange
       Call Set_coeff_lagrange(                                         &
     &                        lambda_p, row_length, halo_i, .false.,    &
     &                        lambda_p_rm,  lambda_p_rp,                &
     &                        recip_lambda_p_m,  recip_lambda_p,        &
     &                        recip_lambda_p_p, recip_lambda_p_p2 )
     
! DEPENDS ON: Set_coeff_lagrange
       Call Set_coeff_lagrange(                                         &
     &                        lambda_u, row_length, halo_i, .false.,    &
     &                        lambda_u_rm,  lambda_u_rp,                &
     &                        recip_lambda_u_m,  recip_lambda_u,        &
     &                        recip_lambda_u_p, recip_lambda_u_p2 )

      If (PrintStatus >= PrStatus_Normal ) Then
        write(6,*)
        write(6,*) ' calling Set_coeff_lagrange for phi_p/v'
      End If
      
! DEPENDS ON: Set_coeff_lagrange
       Call Set_coeff_lagrange(                                         &
     &                        phi_p, rows, halo_j, .false.,             &
     &                        phi_p_rm,  phi_p_rp,                      &
     &                        recip_phi_p_m,  recip_phi_p,              &
     &                        recip_phi_p_p, recip_phi_p_p2 )


! DEPENDS ON: Set_coeff_lagrange
       Call Set_coeff_lagrange(                                         &
     &                        phi_v, n_rows, halo_j, .false.,           &
     &                        phi_v_rm,  phi_v_rp,                      &
     &                        recip_phi_v_m,  recip_phi_v,              &
     &                        recip_phi_v_p, recip_phi_v_p2 )
                       
!  Set 2d phi arrays
      do j = 1-halo_j , rows + halo_j 
        do i = 1-halo_i , row_length + halo_i 
          phi_p2d(i,j) = phi_p(j)  
          dphi_p2d(i,j) = dphi_p(j)
          recip_dphip2d(i,j) = recip_dphip(j)
          phi_p_rm2d(i,j) = phi_p_rm(j) 
          phi_p_rp2d(i,j) = phi_p_rp(j) 
          recip_phi_p_m2d(i,j) = recip_phi_p_m(j) 
          recip_phi_p2d(i,j) = recip_phi_p(j) 
          recip_phi_p_p2d(i,j) = recip_phi_p_p(j) 
          recip_phi_p_p22d(i,j) = recip_phi_p_p2(j) 
        end do 
      end do 
 
      do j = 1-halo_j , n_rows + halo_j 
        do i = 1-halo_i , row_length + halo_i 
          phi_v2d(i,j) = phi_v(j) 
          dphi_v2d(i,j) =  dphi_v(j)
          recip_dphiv2d(i,j) = recip_dphiv(j) 
          phi_v_rm2d(i,j) = phi_v_rm(j)
          phi_v_rp2d(i,j) = phi_v_rp(j)
          recip_phi_v_m2d(i,j) = recip_phi_v_m(j)
          recip_phi_v2d(i,j) = recip_phi_v(j)
          recip_phi_v_p2d(i,j) = recip_phi_v_p(j)
          recip_phi_v_p22d(i,j) = recip_phi_v_p2(j)
        end do 
      end do 

!  Set weight arrays (only small halos are needed
!    so 1-halo values are not set since -halo values appear on RHS)
        wt_lambda_p(1-halo_i) = 0.0
        wt_lambda_u(1-halo_i) = 0.0
        Do i = 1-halo_i, row_length+halo_i
          wt_phi_p2d(i,1-halo_j) = 0.0
          wt_phi_v2d(i,1-halo_j) = 0.0
        enddo
        do i = 2-halo_i, row_length+halo_i
          wt_lambda_p(i) = (lambda_p(i) - lambda_u(i-1)) /              &
     &                                                   dlambda_p(i-1)
          wt_lambda_u(i) = (lambda_u(i) - lambda_p(i)) /                &
     &                                                   dlambda_u(i-1)
        enddo
        do j = 2-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            wt_phi_p2d(i,j) = (phi_p(j) - phi_v(j-1)) / dphi_p(j-1)
          enddo
        enddo
        do j = 2-halo_j, n_rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            wt_phi_v2d(i,j) = (phi_v(j) - phi_p(j)) / dphi_v(j-1)
          enddo
        enddo
     
! DEPENDS ON: set_var_look
      call set_var_look(                                                 &
     &                  wk_lambda_p, global_row_length, halo_i,          &
     &                  recip_dlam, look_lam, max_pts, halo_lam,         &
     &                  .false. )
     
! DEPENDS ON: set_var_look
      call set_var_look(                                                 &
     &                  wk_phi_p, global_rows, halo_j,                   &
     &                  recip_dphi, look_phi, max_pts, halo_phi,         &
     &                  .false. )

      IF (lhook) CALL dr_hook('SET_VAR_GRID',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Set_var_grid

