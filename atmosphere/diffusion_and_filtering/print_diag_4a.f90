! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine Print_diag
MODULE print_diag_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE Print_diag_4A(                                         &
     &                     u, v, theta, rho_r2, w, q, qcl, qcf,         &
     &                     rows, n_rows, row_length, model_levels,      &
     &                     off_x, off_y,timestep_number, print_step,    &
     &                     diag_interval,                               &
     &                     rpemax, rpemin, ipesum, rpesum, w_limit,     &
     &                     L_print_pe, L_print_w,                       &
     &                     L_print_wmax, L_print_max_wind,              &
     &                     L_print_div, L_print_lapse, L_print_theta1,  &
     &                     L_print_shear, L_diag_wind, L_diag_noise,    &
     &                     max_w_run, max_wind_run, min_theta1_run,     &
     &                     dtheta1_run, max_div_run, min_div_run,       &
     &                     min_lapse_run, max_shear_run, time_max_shear,&
     &                     time_div_max, time_div_min, time_lapse_min,  &
     &                     time_w_max, time_max_wind, time_theta1_min,  &
     &                     max_KE_run, min_KE_run, max_noise_run,       &
     &                     time_KE_max, time_KE_min, time_noise_max )
 
      USE trignometric_mod, ONLY: sec_theta_latitude,                   &
                            cos_theta_latitude, cos_v_latitude,         &
                            cos_theta_longitude, sin_theta_longitude
      USE proc_info_mod
      USE horiz_grid_mod, ONLY: delta_xi1, delta_xi2

      USE level_heights_mod

! Purpose:
!          Diagnostic routine for w, divergence, lapse_rates
!             and level 1 theta
!          v-at-poles version
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      USE atm_fields_bounds_mod, ONLY:                                  &
          udims, vdims, pdims, udims_s, vdims_s, tdims_s, pdims_s,      &
          qdims_l,wdims_s


      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      USE um_parparams, ONLY: pnorth, peast, psouth, pwest
      USE domain_params 
      USE u_to_p_mod, ONLY: u_to_p
      USE v_to_p_mod, ONLY: v_to_p
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer , Intent(IN) ::                                           &
     &                rows, n_rows, row_length,off_x,off_y, model_levels

      Integer , Intent(IN) ::                                           &
     &  timestep_number                                                 &
     &, print_step                                                      &
     &, diag_interval                                                   &
! the following 4 sizes are set in SETCONA
     &, rpemax                                                          & 
                  ! total size of array needed to find max over pe's
     &, rpemin                                                          & 
                  ! total size of array needed to find min over pe's
     &, ipesum                                                          & 
                  ! total size of integer array needed to sum over pe's
     &, rpesum    ! total size of real array needed to sum over pe's

      Real , Intent(IN) ::                                              &
     &  w_limit                                                         
                   ! w threshold for print
       

      Integer , Intent(INOUT) ::                                        &
     &  time_theta1_min                                                 &
     &, time_ke_max(model_levels + 1)                                   &
     &, time_ke_min(model_levels + 1)                                   &
     &, time_w_max(model_levels)                                        &
     &, time_max_shear(model_levels)                                    &
     &, time_max_wind(model_levels)                                     &
     &, time_div_max(model_levels)                                      &
     &, time_div_min(model_levels)                                      &
     &, time_lapse_min(model_levels)                                    &
     &, time_noise_max(model_levels)

      Real , Intent(INOUT) ::                                           &
     &  min_theta1_run                                                  &
     &, dtheta1_run                                                     &
     &, max_ke_run(model_levels + 1)                                    &
     &, min_ke_run(model_levels + 1)                                    &
     &, max_w_run(0:model_levels)                                         &
     &, max_shear_run(model_levels)                                     &
     &, max_wind_run(model_levels)                                      &
     &, max_div_run(model_levels)                                       &
     &, min_div_run(model_levels)                                       &
     &, min_lapse_run(model_levels)                                     &
     &, max_noise_run(model_levels)


      Real , Intent(IN) ::                                              &
     &  u(udims_s%i_start:udims_s%i_end,                                & 
     &    udims_s%j_start:udims_s%j_end,                                & 
     &    udims_s%k_start:udims_s%k_end)                                &
     &, v(vdims_s%i_start:vdims_s%i_end,                                & 
     &    vdims_s%j_start:vdims_s%j_end,                                & 
     &    vdims_s%k_start:vdims_s%k_end)                                &
     &, w(wdims_s%i_start:wdims_s%i_end,                                & 
     &    wdims_s%j_start:wdims_s%j_end,                                & 
     &    wdims_s%k_start:wdims_s%k_end)                                &
     &, theta(tdims_s%i_start:tdims_s%i_end,                            & 
     &        tdims_s%j_start:tdims_s%j_end,                            & 
     &        tdims_s%k_start:tdims_s%k_end)                            &
     &, rho_r2(pdims_s%i_start:pdims_s%i_end,                           & 
     &         pdims_s%j_start:pdims_s%j_end,                           & 
     &         pdims_s%k_start:pdims_s%k_end)                           &
     &,      q(qdims_l%i_start:qdims_l%i_end,                           & 
     &         qdims_l%j_start:qdims_l%j_end,                           & 
     &         qdims_l%k_start:qdims_l%k_end)                           &
     &,    qcl(qdims_l%i_start:qdims_l%i_end,                           & 
     &         qdims_l%j_start:qdims_l%j_end,                           & 
     &         qdims_l%k_start:qdims_l%k_end)                           &
     &,    qcf(qdims_l%i_start:qdims_l%i_end,                           & 
     &         qdims_l%j_start:qdims_l%j_end,                           & 
     &         qdims_l%k_start:qdims_l%k_end)

      Logical , Intent(IN) ::                                           &
     &  L_print_w                                                       &
                         ! Print control
     &, L_print_wmax                                                    &
                         ! Print control
     &, L_print_shear                                                   &
                          ! Print control
     &, L_print_max_wind                                                &
                             ! Print control
     &, L_print_div                                                     &
                         ! Print control
     &, L_print_lapse                                                   &
                         ! Print control
     &, L_print_theta1                                                  &
                         ! Print control
     &, L_diag_wind                                                     &
                         ! KE diagnostics switch
     &, L_diag_noise                                                    &
                         ! noise diagnostics switch
     &, L_print_pe      ! print diagnostics on all pe's if true
! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                     ! Loop counters
     &, gi, gj                                                          &
                    ! global pointers
     &, ki, kr                                                          & 
                    ! pointers for arrays summed over pe's
     &, kminr                                                           & 
                   ! pointers for arrays summed over pe's
     &, kmaxr                                                           & 
                   ! pointers for arrays summed over pe's
     &, j_start_l, j_stop_l                                            &
                                ! Loop indices
     &, info                                                            &
     &, size                                                            &
                 !  size needed for length for array sums in stats
     &, level_max                                                       &
     &, pe                                                              &
     &, pe_max                                                          &
     &, level_max_run                                                   &
     &, time_max_run                                                    &
     &, i_min_theta                                                     &
                        ! pointers for theta prints level 1
     &, j_min_theta     ! pointers for theta prints level 1

      Real                                                              &
     &  recip_delta_xi1                                                 &
     &, recip_delta_xi2                                                 &
     &, sum_n                                                           &
     &, sum_s                                                           &
     &, min_theta_pe                                                    &
                     ! min theta level 1 on this pe
     &, lambda_max                                                      &
     &, phi_max                                                         &
     &, max_w                                                           &
     &, max_run                                                         &
     &, wind                                                            &
     &, dvdz_at_u                                                       &
     &, dudz_at_v                                                       &
     &, shear                                                           &
     &, dtheta1                                                         &
     &, dtheta1n, dtheta1s, dtheta1e, dtheta1w

! Local arrays

      Real, Dimension (:,:), Allocatable ::                             &
     &  work1                                                           &
     &, dudz                                                            &
     &, dvdz

      Real, Dimension (:,:,:), Allocatable ::                           &
     &  u_rho                                                           &
     &, v_rho                                                           &
     &, w_rho                                                           &
     &, rho_dry   

      Real                                                              &
     &  sumr(rpesum)                                                    & 
                     ! array  for summing integers over pe's
     &, l_n_poles(row_length)                                           &
     &, l_s_poles(row_length)                                           &
     &, mag_vector_np (model_levels)                                    &
     &, dir_vector_np (model_levels)                                    &
     &, mag_vector_sp (model_levels)                                    &
     &, dir_vector_sp (model_levels)                                    &
     &, max_wind_pe(model_levels)                                       &
                                  ! max wind each level on this pe
     &, max_w_pe(0:model_levels - 1)                                      &
                                   ! max w each level on this pe
     &, max_shear_pe(model_levels - 1)                                  &
                                       ! max w each level on this pe
     &, min_lap_pe(model_levels)                                        &
                                 ! min lapse rate on this pe
     &, max_div_pe(model_levels)                                        &
                                 ! max div rate on this pe
     &, min_div_pe(model_levels)                                        &
                                 ! min div rate on this pe
     &, max_real(rpemax)                                                & 
                         ! array  for finding max over pe's
     &, min_real(rpemin)                                                & 
                         ! array  for finding min over pe's
     &, w_mean(model_levels - 1)                                        &
     &, w_variance(model_levels - 1)                                    &
     &, w_std_dev(model_levels - 1)                                     &
     &, w_skew(model_levels - 1)                                        &
     &, w_kurtosis(model_levels - 1)

      Integer                                                           &
     &  sumi(ipesum)                                                    & 
                     ! array  for summing integers over pe's
     &, i_min_lap(model_levels)                                         &
                                ! i pointers for prints
     &, j_min_lap(model_levels)                                         &
                                ! j pointers for prints
     &, i_max_div(model_levels)                                         &
                                ! i pointers for prints
     &, j_max_div(model_levels)                                         &
                                ! j pointers for prints
     &, i_min_div(model_levels)                                         &
                                ! i pointers for prints
     &, j_min_div(model_levels)                                         &
                                ! j pointers for prints
     &, i_max_wind(model_levels)                                        &
                                 ! i pointers for prints
     &, j_max_wind(model_levels)                                        &
                                 ! j pointers for prints
     &, i_max_shear(model_levels - 1)                                   &
                                      ! i pointers for prints
     &, j_max_shear(model_levels - 1)                                   &
                                      ! j pointers for prints
     &, i_max_w(0:model_levels - 1)                                       &
                                  ! i pointers for prints
     &, j_max_w(0:model_levels - 1)                                       &
                                  ! j pointers for prints
     &, w_count(0:model_levels - 1)                                       &
                                  ! counter for w > w_limit
     &, lap_count(model_levels) ! counter for lapse rate < 0

      integer,dimension(2)                  :: min_indices,max_indices


      Real                                                              &
     &  rad_to_deg                                                      &
     &, recip_row_length                                                &
     &, recip_rows                                                      &
     &, weight1                                                         &
     &, lambda                                                          &
                        ! longitude of max
     &, phi                                                             &
                        ! latitude of  max
     &, sum_levels_ke

      CHARACTER(LEN=8)                                                       &
     &  l_string                                                        &
     &, p_string                                                        &
     &, string_lon                                                      &
     &, string_lat

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

      IF (lhook) CALL dr_hook('PRINT_DIAG_4A',zhook_in,zhook_handle)
      If( mod( timestep_number , diag_interval )  ==  0) then

! ----------------------------------------------------------------------
! Section 0. Initialise values
! ----------------------------------------------------------------------
        rad_to_deg = 180.0 / Pi
        recip_row_length = 1.0 / real(global_row_length)
        recip_rows = 1.0 / real(global_rows)
        recip_delta_xi1 = 1. / delta_xi1
        recip_delta_xi2 = 1. / delta_xi2

        if ( L_print_div ) then
! DEPENDS ON: swap_bounds
          CALL swap_bounds(rho_r2, row_length,rows, model_levels,        &
     &                       off_x, off_y, fld_type_p,.false.)
         end if

        if ( L_print_div .OR. L_print_shear .OR. L_diag_wind          &
     &       .OR. L_print_max_wind) then
! DEPENDS ON: swap_bounds
          CALL swap_bounds(v, row_length,n_rows, model_levels,        &
     &                       off_x, off_y, fld_type_v,.true.)
! DEPENDS ON: swap_bounds
          CALL swap_bounds(u, row_length, rows, model_levels,          &
     &                       off_x,off_y,fld_type_u,.true.)
         end if

! ----------------------------------------------------------------------
! Section 1  !  Initialise counters
! ----------------------------------------------------------------------
        kminr = 0
        kmaxr = 0

! ----------------------------------------------------------------------
! Section 2  Find theta min on level 1
! ----------------------------------------------------------------------
        if ( L_print_theta1 ) then

          min_indices = minloc(theta(1:row_length,1:rows,1))
          i_min_theta = min_indices(1)
          j_min_theta = min_indices(2)
          min_theta_pe = theta(i_min_theta,j_min_theta,1)
! copy into min_real(1) which will hold global min later
          kminr = kminr + 1
          min_real(kminr) = min_theta_pe

! find dtheta's surrounding minimum (all will be negative)
          dtheta1e = theta(i_min_theta,j_min_theta,1) -                 &
     &               theta(i_min_theta+1,j_min_theta,1)
          dtheta1w = theta(i_min_theta,j_min_theta,1) -                 &
     &               theta(i_min_theta-1,j_min_theta,1)
          dtheta1n = theta(i_min_theta,j_min_theta,1) -                 &
     &               theta(i_min_theta,j_min_theta+1,1)
          dtheta1s = theta(i_min_theta,j_min_theta,1) -                 &
     &               theta(i_min_theta,j_min_theta-1,1)

        endif ! L_print_theta1

! ----------------------------------------------------------------------
! Section 3. Loop over levels > 2 for static stability
! ----------------------------------------------------------------------

        if ( L_print_lapse ) then

          allocate (work1(row_length,rows))

          lap_count(1) = 0
          Do k = 2, model_levels

            lap_count(k) = 0

            Do j = 1, rows
              Do i = 1, row_length
                work1(i,j) = (theta(i,j,k) - theta(i,j,k-1)) /          &
     &               (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
                if( work1(i,j) < 0.0 ) then
                  lap_count(k) = lap_count(k) + 1
                endif ! work1(i,j) < 0.0
              End Do
            End Do

!Find min lapse_rates on each processor
            min_indices  = minloc(work1(1:row_length,1:rows))
            i_min_lap(k) = min_indices(1)
            j_min_lap(k) = min_indices(2)
            min_lap_pe(k) = work1(i_min_lap(k),j_min_lap(k))

! Copy local min lap into global max/min_real which will hold
!          global min after calling gc_rmin
            kminr = kminr + 1
            min_real(kminr) = min_lap_pe(k)

          End Do  ! k = 2, model_levels

          deallocate ( work1 )

        endif ! L_print_lapse

! ----------------------------------------------------------------------
! Section 4. Calculate divergence
! ----------------------------------------------------------------------
        if ( L_print_div ) then

          allocate (work1(row_length, rows))

          Do  k = 1, model_levels

! Calculate u derivative at rho points
            Do j = 1, rows
              Do i = 1, row_length
                work1(i,j) = 0.5 * recip_delta_xi1 * ( u(i,j,k) *    &
     &                        ( rho_r2(i,j,k) + rho_r2(i+1,j,k) ) -     &
     &                                                 u(i-1,j,k) *     &
     &                       ( rho_r2(i-1,j,k) + rho_r2(i,j,k) ) )*     &
     &                            sec_theta_latitude(i,j)
              End Do
            End Do

! Calculate v derivative at rho points

            Do j = 1, rows           
              Do i = 1, row_length
                work1(i,j) = work1(i,j) + 0.5 * recip_delta_xi2 *       &
     &                                                     ( v(i,j,k) * &
     &                            ( rho_r2(i,j,k) + rho_r2(i,j+1,k) ) * &
     &                               cos_v_latitude(i,j) - v(i,j-1,k) * &
     &                            ( rho_r2(i,j,k) + rho_r2(i,j-1,k) ) * &
     &                                        cos_v_latitude(i,j-1) ) * &
     &                                 sec_theta_latitude(i,j)
              End Do
            End Do

            Do j = 1, rows
              Do i = 1, row_length
                work1(i,j) = work1(i,j)                                 &
     &                        / (rho_r2(i,j,k) * r_rho_levels(i,j,k))
              End Do
            End Do

! Find maximum absolute divergence on each processor
            min_indices  = minloc(work1(1:row_length,1:rows))
            i_min_div(k) = min_indices(1)
            j_min_div(k) = min_indices(2)
            min_div_pe(k) = work1(i_min_div(k),j_min_div(k))

            max_indices  = maxloc(work1(1:row_length,1:rows))
            i_max_div(k) = max_indices(1)
            j_max_div(k) = max_indices(2)
            max_div_pe(k) = work1(i_max_div(k),j_max_div(k))

! Copy local max/min div into global_max/min_real which will hold
!    global max/min_div after calling gc_rmax
            kminr = kminr + 1
            kmaxr = kmaxr + 1
            max_real(kmaxr) = max_div_pe(k)
            min_real(kminr) = min_div_pe(k)

          end Do !  k = 1, model_levels

          deallocate ( work1 )

        endif ! L_print_div

! ----------------------------------------------------------------------
! Section 5   Now find max of w for each level
!               w(model_levels) = 0.0 everywhere so omit
! ----------------------------------------------------------------------
        if (L_print_w .or. L_print_wmax ) then

          Do   k = 0, model_levels - 1

            w_count(k) = 0
            Do j = 1, rows
              Do i = 1, row_length
                if (w(i,j,k) > w_limit) then
                  w_count(k) = w_count(k) + 1
                end if
              End Do
            End Do

            max_indices = maxloc(w(1:row_length,1:rows,k))
            i_max_w(k) = max_indices(1)
            j_max_w(k) = max_indices(2)
            max_w_pe(k) = w(i_max_w(k),j_max_w(k),k)

! Copy local max w into w_max which will hold global w_max after gc_rmax
            kmaxr = kmaxr + 1
            max_real(kmaxr) = max_w_pe(k)

          End Do  ! k = 1, model_levels - 1

        endif ! L_print_w .or. L_print_wmax

! ----------------------------------------------------------------------
! Section 6   Wind shear diagnostics
! ----------------------------------------------------------------------
        if( L_print_shear ) then

        j_start_l = vdims%j_start
        j_stop_l = vdims%j_end
        If (model_domain  ==  mt_global ) Then
          If (at_extremity(PSouth)) j_start_l = vdims%j_start +1
          If (at_extremity(PNorth)) j_stop_l = vdims%j_end -1
        endIf ! model_domain  ==  mt_global

          allocate ( dudz(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end) )
          allocate ( dvdz(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end) )

          Do k = 1, model_levels - 1

            max_shear_pe(k) = 0.0
            i_max_shear(k) = 0
            j_max_shear(k) = 0

            Do j = udims_s%j_start,udims_s%j_end
              Do i = udims_s%i_start,udims_s%i_end
                dudz(i,j) = ( u(i,j,k+1) - u(i,j,k) ) /                 &
     &                     ( r_at_u(i,j,k+1) - r_at_u(i,j,k) )
              End Do
            End Do

            Do j = vdims_s%j_start,vdims_s%j_end
              Do i = vdims_s%i_start,vdims_s%i_end
                dvdz(i,j) = ( v(i,j,k+1) - v(i,j,k) ) /                 &
     &                     ( r_at_v(i,j,k+1) - r_at_v(i,j,k) )
              End Do
            End Do

            Do j = udims%j_start,udims%j_end
              Do i = udims%i_start,udims%i_end
                dvdz_at_u = 0.25 * ( dvdz(i,j) + dvdz(i+1,j) +          &
     &                               dvdz(i,j-1) + dvdz(i+1,j-1) )
                shear = sqrt( dudz(i,j) * dudz(i,j) +                   &
     &                       dvdz_at_u * dvdz_at_u )
                if( shear > max_shear_pe(k) ) then
                  max_shear_pe(k) = shear
                  i_max_shear(k) = i
                  j_max_shear(k) = j
                end if ! shear > max_shear_pe(k)
              End Do
            End Do

            Do j = j_start_l, j_stop_l
              Do i = vdims%i_start,vdims%i_end
                dudz_at_v = 0.25 * ( dudz(i,j  ) + dudz(i-1,j  ) +      &
     &                               dudz(i,j+1) + dudz(i-1,j+1) )
                shear = sqrt( dvdz(i,j) * dvdz(i,j) +                   &
     &                        dudz_at_v * dudz_at_v )
                if( shear > max_shear_pe(k) ) then
                  max_shear_pe(k) = shear
                  i_max_shear(k) = i
                  j_max_shear(k) = j
                end if ! shear > max_shear_pe(k)
              End Do
            End Do

! Copy local max shear into global_max_real which will hold
!    global max_shear after calling gc_rmax
            kmaxr = kmaxr + 1
            max_real(kmaxr) = max_shear_pe(k)

          EndDo  ! k = 1, model_levels - 1

          deallocate ( dudz )
          deallocate ( dvdz )

        endif ! L_print_shear

! ----------------------------------------------------------------------
! Section 7.  KE diagnostics and max_wind
! ----------------------------------------------------------------------

!  Initialise pointers in arrays for summing over pe's
        ki = 0
        kr = 0

        if ( L_diag_wind .or. L_print_max_wind ) then

          allocate ( u_rho(row_length, rows, model_levels) )
          allocate ( v_rho(row_length, rows, model_levels) )

!-------------------------------------------------------------------
! 7.1 Interpolate winds to rho points (already on same vertical level)
!-------------------------------------------------------------------

      CALL v_to_p(v,                                                   &
                        vdims_s%i_start,vdims_s%i_end,                  & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity,v_rho)

      CALL u_to_p(u,                                                   &
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity,u_rho)

!-------------------------------------------------------------------
! 7.2  Find max wind on this pe
!-------------------------------------------------------------------
          if ( L_print_max_wind ) then

          allocate (work1(1-off_x:row_length+off_x, 1-off_y:rows+off_y))

            DO k = 1, model_levels
              max_wind_pe(k) = 0.0
              i_max_wind(k) = 0
              j_max_wind(k) = 0
              do j = 1, rows         
                do i = 1, row_length
                  work1(i,j) = sqrt( u_rho(i,j,k) * u_rho(i,j,k) +      &
     &                             v_rho(i,j,k) * v_rho(i,j,k) )
                End do
              End do

              max_indices = maxloc(work1(1:row_length,1:rows))
              i_max_wind(k) = max_indices(1)
              j_max_wind(k) = max_indices(2)   !! + j_start - 1
              max_wind_pe(k) = work1(i_max_wind(k),j_max_wind(k))

! Copy max wind_local into max_real which for global gc_rmax later
              kmaxr = kmaxr + 1
              max_real(kmaxr) = max_wind_pe(k)

            End do  !  k = 1, model_levels

            deallocate (work1)

          end if ! L_print_max_wind

!-------------------------------------------------------------------
! 7.3  Polar processing for max_wind and wind_diag
!-------------------------------------------------------------------
! Uses v on row next to poles to calculate wind at the polar point

          if ( .not. L_diag_wind ) then
            deallocate ( u_rho )
            deallocate ( v_rho )
          endif ! .not. L_diag_wind

        endif ! L_diag_wind .or. L_print_max_wind

!-------------------------------------------------------------------
! 7.4  Now complete wind_diag
!-------------------------------------------------------------------
        if ( L_diag_wind ) then

          allocate ( w_rho(row_length, rows, model_levels) )
          allocate ( rho_dry(row_length, rows, model_levels) )

          k = 1
          kr = kr + 1
!  Initialise sumr for ke summing later
          sumr(kr) = 0.0
          Do j = 1, rows
            Do i = 1, row_length
              weight1 = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k)) / &
     &                 (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
              rho_dry(i,j,k) = rho_r2(i,j,k) *                          &
     &                (1. - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))
              w_rho(i,j,k) = ( 1.0 -  weight1 ) * w(i,j,k)
            End Do
          End Do

          Do k = 2, model_levels
            kr = kr + 1
            sumr(kr) = 0.0
            Do j = 1, rows
              Do i = 1, row_length
              weight1 = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k)) / &
     &                 (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
              rho_dry(i,j,k) = rho_r2(i,j,k) * ( (1. - weight1) *       &
     &                      (1. - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k) ) +&
     &                                                         weight1 *&
     &                (1. - q(i,j,k-1) - qcl(i,j,k-1) - qcf(i,j,k-1) ) )
              w_rho(i,j,k) = (1.0 -  weight1) * w(i,j,k) +              &
     &                                 weight1 * w(i,j,k-1)
              End Do
            End Do
          End Do  !  k = 2, model_levels

! Problem of u & v values at poles
! Uses v on row next to poles to calculate wind at the polar point
! reset sum array counter
        kr = kr - model_levels

!----------------------------------------------------------------------
! Integrals over fields using values interpolated to rho grid
! At present all integrals done after interpolation to rho points.
! Integrals for energy involve using dry mass.
!----------------------------------------------------------------------
! KE terms for u, v & w
        DO k = 1, model_levels

         kr = kr + 1
         do j = 1, rows        
            do i = 1, row_length
              sumr(kr) = sumr(kr) + 0.5 * rho_dry(i,j,k) *              &
     &                                       delta_xi1 * delta_xi2 * &
     &              (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)) * &
     &                                  ( u_rho(i,j,k) * u_rho(i,j,k) + &
     &                                    v_rho(i,j,k) * v_rho(i,j,k) + &
     &                                  w_rho(i,j,k) * w_rho(i,j,k) ) * &
     &                                 cos_theta_latitude(i,j)
            End do
          End do

        End do  !  k = 1, model_levels

        deallocate ( u_rho )
        deallocate ( v_rho )
        deallocate ( w_rho )
        deallocate ( rho_dry )

      endif ! L_diag_wind

! Global sum for KE done with all the other sums over processors

! ----------------------------------------------------------------------
! Section 8. Now find max/mins over all processors
!            All required fields done at same time
! ----------------------------------------------------------------------

      if ( kmaxr > 0 ) then
        call gc_rmax(kmaxr, n_proc, info, max_real)
      endif !  kmaxr > 0
      if ( kminr > 0 ) then
        call gc_rmin(kminr, n_proc, info, min_real)
      endif !  kminr > 0

!  Re-initialise pointers in max/min arrays
      kminr = 0
      kmaxr = 0

! ----------------------------------------------------------------------
! Section 9. Now find time and place of max/mins
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 9.1  Obtain max and mins from max_real, min_real
!              and fill sumi, sumr arrays for summing over pe's
! ----------------------------------------------------------------------

      if ( L_print_theta1 ) then
        kminr = kminr + 1
        dtheta1 = min(dtheta1e, dtheta1w, dtheta1n, dtheta1s)
        ki = ki + 1
        kr = kr + 3
        sumi(ki) = 0
        sumr(kr-2) = 0.0
        sumr(kr-1) = 0.0
        sumr(kr) = 0.0
        if( min_theta_pe <= min_real(kminr) )then
          sumi(ki) = me ! min is on this processor
! store largest negative gradient at min theta1 point
          sumr(kr-2) = dtheta1
          gi = l_datastart(1) + i_min_theta - 1
          gj = l_datastart(2) + j_min_theta - 1
          If ( model_domain  ==  mt_Global) then
            sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
            sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
          else ! for LAM output fraction of point
            sumr(kr-1) = (gi-1) * recip_row_length * 100.0
            sumr(kr) = (gj-1) * recip_rows * 100.0
          endif ! model_domain  ==  mt_Global
        endif ! min_theta_pe <= min_real(kminr)

      endif ! L_print_theta1

      if ( L_print_lapse ) then
        Do k = 2, model_levels
          kminr = kminr + 1
          if( min_real(kminr) < min_lapse_run(k) )then
              min_lapse_run(k) = min_real(kminr)
              time_lapse_min(k) = timestep_number
          endif ! min_real(kminr) > min_lapse_run(k)
          ki = ki + 2
          kr = kr + 2
          sumi(ki - 1) = lap_count(k)
          sumi(ki) = 0
          sumr(kr-1) = 0.0
          sumr(kr) = 0.0
          if ( min_lap_pe(k) <= min_real(kminr)) then
            sumi(ki) = me ! min is on this processor for this level
            gi = l_datastart(1) + i_min_lap(k) - 1
            gj = l_datastart(2) + j_min_lap(k) - 1
            If ( model_domain  ==  mt_Global) then
              sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
              sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
            else ! for LAM output fraction of point
              sumr(kr-1) = (gi-1) * recip_row_length * 100.0
              sumr(kr) = (gj-1) * recip_rows * 100.0
            endif ! model_domain  ==  mt_Global
          endif  ! min_lap_pe(k) <= min_real(kminr)
        end Do ! k = 2, model_levels

      endif ! L_print_lapse

      if ( L_print_div ) then

        Do k = 1, model_levels
          kmaxr = kmaxr + 1
          if( max_real(kmaxr) > max_div_run(k) )then
            max_div_run(k) =  max_real(kmaxr)
            time_div_max(k) = timestep_number
          endif ! max_real(kminr)  > max_div_run(k)
          ki = ki + 1
          kr = kr + 2
          sumi(ki) = 0
          sumr(kr-1) = 0.0
          sumr(kr) = 0.0
          if ( max_div_pe(k) >= max_real(kmaxr)) then
            sumi(ki) = me ! max is on this processor for this level
            gi = l_datastart(1) + i_max_div(k) - 1
            gj = l_datastart(2) + j_max_div(k) - 1
            If ( model_domain  ==  mt_Global) then
              sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
              sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
            else ! for LAM output fraction of point
              sumr(kr-1) = (gi-1) * recip_row_length * 100.0
              sumr(kr) = (gj-1) * recip_rows * 100.0
            endif ! model_domain  ==  mt_Global
          endif  ! max_div_pe(k) >= max_real(kmaxr)
        end Do  !  k = 1, model_levels

        Do k = 1, model_levels
          kminr = kminr + 1
          if( min_real(kminr) < min_div_run(k) )then
            min_div_run(k) =  min_real(kminr)
            time_div_min(k) = timestep_number
          endif ! min_real(kminr) > min_div_run(k)
          ki = ki + 1
          kr = kr + 2
          sumi(ki) = 0
          sumr(kr-1) = 0.0
          sumr(kr) = 0.0
          if ( min_div_pe(k) <= min_real(kminr)) then
            sumi(ki) = me ! min is on this processor for this level
            gi = l_datastart(1) + i_min_div(k) - 1
            gj = l_datastart(2) + j_min_div(k) - 1
            If ( model_domain  ==  mt_Global) then
              sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
              sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
            else ! for LAM output fraction of point
              sumr(kr-1) = (gi-1) * recip_row_length * 100.0
              sumr(kr) = (gj-1) * recip_rows * 100.0
            endif ! model_domain  ==  mt_Global
          endif  ! min_div_pe(k) <= min_real(kminr)
        end Do  !  k = 1, model_levels

      endif ! L_print_div

      if ( L_print_w .or. L_print_wmax ) then
        Do   k = 0, model_levels - 1
          kmaxr = kmaxr + 1
          if( max_real(kmaxr) > max_w_run(k) )then
             max_w_run(k) =  max_real(kmaxr)
             time_w_max(k+1) = timestep_number
          endif ! max_real(kmaxr) > max_w_run(k)
          ki = ki + 2
          kr = kr + 2
          sumi(ki - 1) = w_count(k)
          sumi(ki) = 0
          sumr(kr-1) = 0.0
          sumr(kr) = 0.0
          if ( max_w_pe(k) >= max_real(kmaxr)) then
            sumi(ki) = me ! max is on this processor for this level
            gi = l_datastart(1) + i_max_w(k) - 1
            gj = l_datastart(2) + j_max_w(k) - 1
            If ( model_domain  ==  mt_Global) then
              sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
              sumr(kr) = ( (gj-0.5) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
            else ! for LAM output fraction of point
              sumr(kr-1) = (gi-1) * recip_row_length * 100.0
              sumr(kr) = (gj-1) * recip_rows * 100.0
            endif ! model_domain  ==  mt_Global
          endif  ! max_w_pe(k) >= max_real(kmaxr)
        End do  ! k = 1, model_levels - 1
      endif ! L_print_w .or. L_print_wmax

      if ( L_print_shear ) then
        Do   k = 1, model_levels - 1
          kmaxr = kmaxr + 1
          if( max_real(kmaxr) > max_shear_run(k) )then
             max_shear_run(k) =  max_real(kmaxr)
             time_max_shear(k) = timestep_number
          endif ! max_real(kmaxr) > max_shear_run(k)
          ki = ki + 1
          kr = kr + 2
          sumi(ki) = 0
          sumr(kr-1) = 0.0
          sumr(kr) = 0.0
          if ( max_shear_pe(k) >= max_real(kmaxr)) then
            sumi(ki) = me ! max is on this processor for this level
            gi = l_datastart(1) + i_max_shear(k) - 1
            gj = l_datastart(2) + j_max_shear(k) - 1
            If ( model_domain  ==  mt_Global) then
              sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
              sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
            else ! for LAM output fraction of point
              sumr(kr-1) = (gi-1) * recip_row_length * 100.0
              sumr(kr) = (gj-1) * recip_rows * 100.0
            endif ! model_domain  ==  mt_Global
          endif  ! max_shear_pe(k) >= max_real(kmaxr)
        End do  ! k = 1, model_levels - 1
      endif ! L_print_shear

! ----------------------------------------------------------------------

      if ( L_print_max_wind ) then
        Do   k = 1, model_levels
          kmaxr = kmaxr + 1
          if( max_real(kmaxr) > max_wind_run(k) )then
             max_wind_run(k) =  max_real(kmaxr)
             time_max_wind(k) = timestep_number
          endif ! max_real(kmaxr) > max_wind_run(k)
          ki = ki + 1
          kr = kr + 2
          sumi(ki) = 0
          sumr(kr-1) = 0.0
          sumr(kr) = 0.0
          if ( max_wind_pe(k) >= max_real(kmaxr)) then
            sumi(ki) = me ! max is on this processor for this level
            gi = l_datastart(1) + i_max_wind(k) - 1
            gj = l_datastart(2) + j_max_wind(k) - 1
            If ( model_domain  ==  mt_Global) then
              sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
              sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
            else ! for LAM output fraction of point
              sumr(kr-1) = (gi-1) * recip_row_length * 100.0
              sumr(kr) = (gj-1) * recip_rows * 100.0
            endif ! model_domain  ==  mt_Global
          endif  ! max_wind_pe(k) <= max_real(kmaxr)
        End do  ! k = 1, model_levels
      endif ! L_print_max_wind

! ----------------------------------------------------------------------
!  Printing will be done in section 12
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 10  Summing over pe's to obtain sums and location of max/mins
! ----------------------------------------------------------------------

! We do not need a reproducible sum, as each element of sumr consists of
! one non-zero entry and all other entries as zero.
! Adding n_proc-1 zero's to X give X whether you use a reproducible sum
! or not.

      if( ki > 0 ) call gc_isum(ki, n_proc, info, sumi)
      if( kr > 0 ) call gc_rsum(kr, n_proc, info, sumr)

! ----------------------------------------------------------------------
! Section 10.1  Obtain max/mins for KE
! ----------------------------------------------------------------------
        if ( L_diag_wind ) then

! Re-initialise pointer
          kr = 0
          sum_levels_ke = 0.0
          Do   k = 1, model_levels
            kr = kr + 1
            if( sumr(kr) > max_ke_run(k) ) then
              max_ke_run(k) = sumr(kr)
              time_ke_max(k) = timestep_number
            endif ! sumr(kr) > max_ke_run(k)
            if( sumr(kr) < min_ke_run(k) ) then
              min_ke_run(k) =  sumr(kr)
              time_ke_min(k) = timestep_number
            endif ! sumr(kr) > max_ke_run(k)
            sum_levels_ke = sum_levels_ke + sumr(kr)
          enddo   !  k = 1, model_levels

          k = model_levels + 1  ! Put vertical sum in  model_levels + 1
          if( sum_levels_ke > max_ke_run(k) ) then
            max_ke_run(k) =  sum_levels_ke
            time_ke_max(k) = timestep_number
          endif ! sum_levels_ke > max_ke_run(k)
          if( sum_levels_ke < min_ke_run(k) ) then
            min_ke_run(k) =  sum_levels_ke
            time_ke_min(k) = timestep_number
          endif ! sum_levels_ke < min_ke_run(k)

        endif ! L_diag_wind

! ----------------------------------------------------------------------
! Section 11. calculate stats for w field to assess noise
! ----------------------------------------------------------------------

      if ( L_diag_noise ) then

        size = 4 * (model_levels - 1)

! DEPENDS ON: calc_stats
        call calc_stats(                                                &
     &                     w(1-off_x, 1-off_y, 1),                      &
     &                     row_length, rows, model_levels - 1,          &
     &                     off_x, off_y, n_proc, size,                  &
     &                     w_mean, w_variance, w_std_dev,               &
     &                     w_skew, w_kurtosis)

      endif ! L_diag_noise

      EndIf ! mod( timestep_number , diag_interval )  ==  0

! ----------------------------------------------------------------------
! Section 12. Print diagnostic information
! ----------------------------------------------------------------------

      If( mod( timestep_number , print_step ) == 0                      &
     &   .and.   (L_print_pe .or. me == 0)       ) then

! Re-initialise pointers
        kminr = 0
        kmaxr = 0
        ki = 0
        kr = 0

        if ( L_diag_noise ) then

          write(6,*) '  '
          write(6,*) ' Vertical velocity characteristics at '           &
     &              ,timestep_number,' time steps'
          write(6,*) '  '
          write(6,*) 'Level   Mean      Variance    std_dev   ',        &
     &             ' Skew      Kurtosis'

          max_w = -100.0
          max_run = -100.0
          Do k = 1, model_levels - 1

            if( w_variance(k) > max_noise_run(k) ) then
              max_noise_run(k) = w_variance(k)
              time_noise_max(k) = timestep_number
            endif ! w_variance(k) > max_noise_run(k)
            if( w_variance(k) > max_w ) then
              max_w = w_variance(k)
              level_max = k
            endif ! w_variance(k) > max_noise
            if( max_noise_run(k) > max_run ) then
              max_run = max_noise_run(k)
              level_max_run = k
              time_max_run = time_noise_max(k)
            endif ! max_noise_run(k) > max_run

            write(6,994) k, w_mean(k), w_variance(k), w_std_dev(k),     &
     &                   w_skew(k), w_kurtosis(k)
          end Do !   k = 1, model_levels - 1

          write(6,*) '    '
          write(6,*) ' Maximum noise (w variance) at timestep ',        &
     &                                      timestep_number
          write(6,*) '   Max noise      level',                         &
     &             '   max this run  level   timestep'
          write(6,992) max_w, level_max,                                &
     &               max_run, level_max_run, time_max_run

        endif ! L_diag_noise

        if ( L_diag_wind ) then
          write(6,*) '  '
          write(6,*) ' KE at ',timestep_number,' time steps'
          write(6,*) 'level   This step     min this run (time step)',  &
     &                                '  max this run  (time step)'
          Do   k = 1, model_levels
            kr = kr + 1
            write(6,993) k, sumr(kr), min_ke_run(k), time_ke_min(k),    &
     &                              max_ke_run(k), time_ke_max(k)
          enddo   !  k = 1, model_levels

          k = model_levels + 1  ! Put vertical sum in  model_levels + 1
          write(6,*)' Total KE    Min KE this run at timestep ',        &
     &                        ' Max KE this run at timestep '
          write(6,991) sum_levels_ke, min_ke_run(k), time_ke_min(k),    &
     &                              max_ke_run(k), time_ke_max(k)
        endif ! L_diag_wind

        if ( L_print_theta1 ) then
          kminr = kminr + 1
          ki = ki + 1
          kr = kr + 3
          pe = sumi(ki)
          dtheta1 = sumr(kr-2)
          lambda = sumr(kr-1)
          phi = sumr(kr)
          if( min_real(kminr) < min_theta1_run )then
            min_theta1_run = min_real(kminr)
            time_theta1_min = timestep_number
            dtheta1_run = dtheta1
          endif ! min_real(kminr) < min_theta1_run
          If ( model_domain  ==  mt_Global) then
            if(lambda > 180.0) then
              lambda = 360.0 - lambda
              l_string ='deg W'
            else
              l_string ='deg E'
            endif      !  lambda(k) > 180.0
            if(phi > 0.0) then
              p_string ='deg N'
            else
              p_string ='deg S'
            endif      !  phi(k) > 0.0
          else ! LAM domains
            p_string ='% North'
            l_string ='% East'
          endIf ! model_domain  ==  mt_Global
          write(6,*) '  '
         write(6,*)'Minimum theta level 1 for timestep ',timestep_number
          write(6,*) '               This timestep',                    &
     &             '                         This run'
          write(6,*) '  Min theta1     proc          position',         &
     &            '            Min theta1 timestep'
          write(6,996)min_real(kminr), pe, lambda, l_string,            &
     &           phi, p_string, min_theta1_run, time_theta1_min
          write(6,*) ' Largest negative delta theta1 at minimum theta1 '
          write(6,989) dtheta1, dtheta1_run

        endif ! L_print_theta1

        if ( L_print_lapse ) then
          write(6,*) '  '
          write(6,*) ' *****  Minimum lapse_rate at timestep ',         &
     &        timestep_number, '     ****** '
          write(6,*) '      number          this timestep',             &
     &             '                      this run'
          write(6,*) ' Level <0   dtheta/dz  proc       position',      &
     &           '            run min    timestep'
          Do k = 2, model_levels
            kminr = kminr + 1
            ki = ki + 2
            kr = kr + 2
            lap_count(1) = sumi(ki - 1)
            pe = sumi(ki)
            lambda = sumr(kr-1)
            phi = sumr(kr)
            If ( model_domain  ==  mt_Global) then
              if(lambda > 180.0) then
                lambda = 360.0 - lambda
                l_string ='deg W'
              else
                l_string ='deg E'
              endif      !  lambda_min(k) > 180.0
              if(phi > 0.0) then
                p_string ='deg N'
              else
                p_string ='deg S'
              endif      !  phi_min(k) > 0.0
            else ! LAM domains
              p_string ='% North'
              l_string ='% East'
            endIf ! model_domain  ==  mt_Global

            write(6,997) k, lap_count(1), min_real(kminr), pe,          &
     &                      lambda, l_string, phi, p_string,            &
     &                   min_lapse_run(k), time_lapse_min(k)

          end Do ! k = 2, model_levels

        endif ! L_print_lapse

        if ( L_print_div ) then

          write(6,*) '  '
          write(6,*) ' ***** Maximum divergence at timestep ',          &
     &                 timestep_number,' ****** '
          write(6,*) '         this timestep',                          &
     &        '                           this run'
          write(6,*) 'level  max div   proc     and   position',        &
     &             '      max div at timestep'

          Do k = 1, model_levels
            kmaxr = kmaxr + 1
            ki = ki + 1
            kr = kr + 2
            pe = sumi(ki)
            lambda = sumr(kr-1)
            phi = sumr(kr)
            If ( model_domain  ==  mt_Global) then
              if(lambda > 180.0) then
                lambda = 360.0 - lambda
                l_string ='deg W'
              else
                l_string ='deg E'
              endif      !  lambda > 180.0
              if(phi > 0.0) then
                p_string ='deg N'
              else
                p_string ='deg S'
              endif      !  phi > 0.0
            else ! LAM domains
              p_string ='% North'
              l_string ='% East'
            endIf ! model_domain  ==  mt_Global

            write(6,999)  k, max_real(kmaxr), pe,                       &
     &              lambda, l_string, phi, p_string,                    &
     &              max_div_run(k), time_div_max(k)

          end Do  !  k = 1, model_levels

          write(6,*) '  '
         write(6,*)' **** Minimum divergence (CONVERGENCE) at timestep '&
     &                ,timestep_number,' ***** '
          write(6,*) '         this timestep',                          &
     &        '                           this run'
          write(6,*) 'level  min div   proc     and   position',        &
     &                           '      min div at timestep'
          Do k = 1, model_levels
            kminr = kminr + 1
            ki = ki + 1
            kr = kr + 2
            pe = sumi(ki)
            lambda = sumr(kr-1)
            phi = sumr(kr)
            If ( model_domain  ==  mt_Global) then
              if(lambda > 180.0) then
                lambda = 360.0 - lambda
                l_string ='deg W'
              else
                l_string ='deg E'
              endif      !  lambda > 180.0
              if(phi > 0.0) then
                p_string ='deg N'
              else
                p_string ='deg S'
              endif      !  phi > 0.0
            else ! LAM domains
              p_string ='% North'
              l_string ='% East'
            endIf ! model_domain  ==  mt_Global
            write(6,999) k, min_real(kminr), pe,                        &
     &                lambda, l_string, phi, p_string,                  &
     &                min_div_run(k), time_div_min(k)
          end Do  !  k = 1, model_levels

        endif ! L_print_div

        if ( L_print_w ) then
          write(6,*) '  '
          write(6,*) ' ***** Maximum vertical velocity w at timestep ', &
     &                       timestep_number,' ****** '
          write(6,990) w_limit
          write(6,*) '         pts >        this timestep',             &
     &            '                     this run'
          write(6,*) 'level w_limit  w_max   proc     and   position',  &
     &          '        w_max at timestep'
        endif ! L_print_w
        if ( L_print_w .or. L_print_wmax ) then
          max_w = -100.0
          max_run = -100.0
          Do   k = 0, model_levels - 1
            kmaxr = kmaxr + 1
            ki = ki + 2
            kr = kr + 2
            w_count(1) = sumi(ki - 1)
            pe = sumi(ki)
            lambda = sumr(kr-1)
            phi = sumr(kr)
            If ( model_domain  ==  mt_Global) then
              if(lambda > 180.0) then
                lambda = 360.0 - lambda
                l_string ='deg W'
              else
                l_string ='deg E'
              endif      !  lambda > 180.0
              if(phi > 0.0) then
                p_string ='deg N'
              else
                p_string ='deg S'
              endif      !  phi > 0.0
            else ! LAM domains
              p_string ='% North'
              l_string ='% East'
            endIf ! model_domain  ==  mt_Global
            if ( max_w < max_real(kmaxr) ) then
              max_w = max_real(kmaxr)
              level_max = k
              pe_max = pe
              lambda_max = lambda
              phi_max = phi
              string_lon = l_string
              string_lat = p_string
            endif ! max_w < max_real(kmaxr)
            if ( max_run < max_w_run(k) ) then
              max_run = max_w_run(k)
              level_max_run = k
              time_max_run = time_w_max(k+1)
            endif ! max_run < max_w_run(k)
            if ( L_print_w ) then
            write(6,997)  k, w_count(1), max_real(kmaxr), pe,           &
     &                    lambda, l_string, phi, p_string,              &
     &                        max_w_run(k), time_w_max(k+1)
            endif ! L_print_w
          endDo   !  k = 1, model_levels - 1
        endif ! L_print_w .or. L_print_wmax

        if ( L_print_wmax ) then
          write(6,*) '  '
          write(6,*) ' Maximum vertical velocity at timestep ',         &
     &            timestep_number,'      Max w this run '
          write(6,*) '   w_max   level  proc         position        ', &
     &           '     run w_max level timestep'
          write(6,995) max_w, level_max, pe_max,                        &
     &            lambda_max, string_lon, phi_max, string_lat,          &
     &            max_run, level_max_run, time_max_run
        endif ! L_print_wmax

        if ( L_print_shear ) then
          write(6,*) '  '
          write(6,*) ' ***** Maximum vertical wind shear at timestep ', &
     &                       timestep_number,' ****** '
          write(6,*) 'level limit  max   proc     and   position',      &
     &          '        max shear at timestep'
          max_w = -100.0
          max_run = -100.0
          Do   k = 1, model_levels - 1
            kmaxr = kmaxr + 1
            ki = ki + 1
            kr = kr + 2
            pe = sumi(ki)
            lambda = sumr(kr-1)
            phi = sumr(kr)
            If ( model_domain  ==  mt_Global) then
              if(lambda > 180.0) then
                lambda = 360.0 - lambda
                l_string ='deg W'
              else
                l_string ='deg E'
              endif      !  lambda > 180.0
              if(phi > 0.0) then
                p_string ='deg N'
              else
                p_string ='deg S'
              endif      !  phi > 0.0
            else ! LAM domains
              p_string ='% North'
              l_string ='% East'
            endIf ! model_domain  ==  mt_Global
            write(6,999) k, max_real(kmaxr), pe,                        &
     &                    lambda, l_string, phi, p_string,              &
     &                        max_shear_run(k), time_max_shear(k)
          endDo   !  k = 1, model_levels - 1
        endif ! L_print_shear

        if ( L_print_max_wind ) then
          write(6,*) '  '
          write(6,*) ' ***** Maximum wind speed at timestep ',          &
     &                       timestep_number,' ****** '
          write(6,*) '                this timestep'                    &
     &           ,'                     this run'
          write(6,*) 'level  max wind   proc     and   position',       &
     &          '      max wind at timestep'
          max_w = 0.0
          max_run = 0.0
          Do   k = 1, model_levels
            kmaxr = kmaxr + 1
            ki = ki + 1
            kr = kr + 2
            pe = sumi(ki)
            lambda = sumr(kr-1)
            phi = sumr(kr)
            If ( model_domain  ==  mt_Global) then
              if(lambda > 180.0) then
                lambda = 360.0 - lambda
                l_string ='deg W'
              else
                l_string ='deg E'
              endif      !  lambda > 180.0
              if(phi > 0.0) then
                p_string ='deg N'
              else
                p_string ='deg S'
              endif      !  phi > 0.0
            else ! LAM domains
              p_string ='% North'
              l_string ='% East'
            endIf ! model_domain  ==  mt_Global
            if ( max_w < max_real(kmaxr) ) then
              max_w = max_real(kmaxr)
              level_max = k
              pe_max = pe
              lambda_max = lambda
              phi_max = phi
              string_lon = l_string
              string_lat = p_string
            endif ! max_w < max_real(kmaxr)
            if ( max_run < max_wind_run(k) ) then
              max_run = max_wind_run(k)
              level_max_run = k
              time_max_run = time_max_wind(k)
            endif ! max_run < max_wind_run(k)
            write(6,999) k,  max_real(kmaxr), pe,                       &
     &                    lambda, l_string, phi, p_string,              &
     &                        max_wind_run(k), time_max_wind(k)
          endDo   !  k = 1, model_levels - 1

          write(6,*) '  '
          write(6,*) ' Maximum horizontal wind at timestep ',           &
     &            timestep_number,'      Max wind this run '
          write(6,*) '   max_wind   level  proc         position   ',   &
     &           '     run max_wind level timestep'
          write(6,995) max_w, level_max, pe_max,                        &
     &            lambda_max, string_lon, phi_max, string_lat,          &
     &            max_run, level_max_run, time_max_run
        endif ! L_print_max_wind

      endif !  mod(timestep_number, print_step) == 0

! ----------------------------------------------------------------------
! Section 13. Print formats
! ----------------------------------------------------------------------

 989   FORMAT(' This timestep = ',F8.2,'K. At min for run = ',F8.2,'K')
 990   FORMAT(' The column below w_limit shows number of points > ',    &
     &        E9.3,' m/s')
 991   FORMAT(1X, 2E15.3,I9,E15.3,I9)
 992   FORMAT(1X, E12.3,I9,E14.3,2I9)
 993   FORMAT(1X, I4, 2E15.3,I9,E15.3,I9)
 994   FORMAT(1X, I4, 5E11.3)
 995   FORMAT(1X, E11.3, I4, I7, F7.1, A7, F7.1, A7, E11.3, I5, I6)
 996   FORMAT(1X, F11.2, I8, F8.1, A7, F8.1, A7, F11.2, I6)
 997   FORMAT(1X, I4, I6, E11.3, I7, F6.1, A7, F6.1, A7, E11.3, I6)
 998   FORMAT(1X, I4, I6, E11.3, I5, F6.1, A7, F6.1, A7, E11.3, I6)
 999   FORMAT(1X, I4, E11.3, I7, F6.1, A7, F6.1, A7, E11.3, I6)

      IF (lhook) CALL dr_hook('PRINT_DIAG_4A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Print_diag_4A
END MODULE
