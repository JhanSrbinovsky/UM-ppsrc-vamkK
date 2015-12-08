! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Setdiff_old
!
! Purpose : to set up polar filtering and diffusion
!
! Language: FORTRAN 77 + common extensions also in Fortran 90.
! Programming standard; Unified Model Documentation Paper No. 3
! version 7.2, dated 5/2/98
!
! Documentation : Unified Model Documentation Paper No P0
!
!----------------------------------------------------------------
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level

      SUBROUTINE Setdiff_old(                                           &
! Input size and control variables
     &      global_rows, row_length, rows, n_rows,                      &
     &      model_domain, at_extremity, datastart,                      &
     &      offx, offy, mype, nproc, max_121_rows,                      &
! other constants
     &      delta_lambda, delta_phi,                                    &
     &      polar_cap, scale_ratio, diff_coeff_ref,                     &
     &      cos_theta_latitude, sin_theta_latitude,                     &
     &      cos_v_latitude, sin_v_latitude,                             &
! Output data
     &       global_u_filter, global_v_filter,                          &
     &       u_sweeps, v_sweeps,                                        &
     &       u_begin, u_end, v_begin, v_end,                            &
     &       diff_coeff_u, diff_coeff_v,                                &
     &       diff_coeff_thermo, diff_coeff_wind,                        &
     &       diff_order_thermo, diff_order_wind,                        &
     &       diff_timescale_thermo, diff_timescale_wind,                &
     &       L_filter, L_pfcomb, L_pftheta, L_pfuv, L_pfw, L_pfincs,    &
     &       L_diff_thermo, L_diff_wind, L_diff_w, L_diff_incs)


      USE conversions_mod, ONLY: pi_over_180, pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

      INTEGER, Intent(In) ::                                            &
     &       global_rows,                                               &
                                 ! IN total number of rows in model
     &       row_length,                                                &
     &       rows,                                                      &
     &       n_rows,                                                    &
     &       offx,                                                      &
     &       offy,                                                      &
     &       mype,                                                      &
     &       nproc,                                                     &
     &       model_domain                                               &
     &, max_121_rows                                                    &
                      ! Max no. of rows 1-2-1 filtered in a hemisphere
     &, diff_order_thermo                                               &
     &, diff_order_wind                                                 &
     &, datastart(3)       ! First gridpoints held by this processor

      REAL, Intent(In) ::                                               &
     &  delta_lambda                                                    &
                     ! EW (x) grid spacing
     &, delta_phi                                                       &
                     ! NS (y) grid spacing
     &, scale_ratio

      LOGICAL, Intent(In) ::                                            &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_pfincs                                                        &
                       ! switch for combined polar filter
     &, L_pftheta                                                       &
                       ! switch for polar filter for theta
     &, L_pfw                                                           &
                       ! switch for polar filter for w
     &, L_pfuv                                                          &
                       ! switch for polar filter for horizontal winds
     &, L_diff_w                                                        &
                       ! switch for horiz. diffusion of w
     &, L_diff_incs    ! switch for horiz. diffusion of incs

      Real, Intent(In) ::                                               &
     &  sin_theta_latitude (row_length, rows)                           &
     &, sin_v_latitude (row_length, n_rows)                             &
     &, cos_theta_latitude(1-offx:row_length+offx,                      &
     &                     1-offy:rows+offy )                           &
     &, cos_v_latitude(1-offx:row_length+offx, 1-offy:n_rows+offy )

  ! Output arguments for diffusion/filtering control

      LOGICAL, Intent(InOut) ::                                         &
     &  L_filter                                                        &
                       ! general switch filter/diffusion
     &, L_pfcomb                                                        &
                       ! general switch for combined polar filter
     &, L_diff_thermo                                                   &
                       ! switch for horiz. diffusion of theta
     &, L_diff_wind    ! switch for horiz. diffusion of u,v

      INTEGER, Intent(InOut) ::                                         &
     &  diff_timescale_thermo                                           &
     &, diff_timescale_wind                                             &
     &, global_u_filter                                                 &
                         ! number of u rows filtered in a hemisphere
     &, global_v_filter                                                 &
                         ! number of v rows filtered in a hemisphere
     &, u_sweeps(max_121_rows)                                          &
                                ! sweeps for 1-2-1 filter
     &, v_sweeps(max_121_rows)                                          &
                                ! sweeps for 1-2-1 filter
     &, u_begin(0:max_121_rows)                                         &
                                 ! row pointers for 1-2-1 filter
     &, u_end(0:max_121_rows)                                           &
                                 ! row pointers for 1-2-1 filter
     &, v_begin(0:max_121_rows)                                         &
                                 ! row pointers for 1-2-1 filter
     &, v_end(0:max_121_rows)    ! row pointers for 1-2-1 filter

      REAL, Intent(InOut) ::                                            &
     &  diff_coeff_ref                                                  &
     &, polar_cap                                                       &
     &, diff_coeff_thermo                                               &
     &, diff_coeff_wind                                                 &
!    diffusion coefficient for u rows
     &, diff_coeff_u(1-offx:row_length+offx, 1-offy:rows+offy)          &
!    diffusion coefficient for v rows
     &, diff_coeff_v(1-offx:row_length+offx, 1-offy:n_rows+offy)


! Local variables
      REAL                                                              &
     &  cos_lat_new                                                     &
     &, cos_polar_cap                                                   &
     &, d_ratio                                                         &
     &, damp_ratio                                                      &
     &, diff_time_equ                                                   &
     &, filter_lat                                                      &
     &, lat_new                                                         &
     &, lat_row                                                         &
     &, Pi_over_2                                                       &
     &, recip_diffusion_order                                           &
     &, scale                                                           &
     &, scale_r2                                                        &
     &, scale_test                                                      &
     &, sin_ratio                                                       &
     &, tiny                                                            &
     &, diff_coeff(global_rows)

      LOGICAL                                                           &
     &  L_rows_even

      INTEGER                                                           &
              ! Mostly loop counters, but j also used for interp points
     &  I, J, GJ                                                        &
     &, j_start_u, j_stop_u                                             &
     &, j_start_v, j_stop_v                                             &
     &, n_sweeps                                                        &
     &, max_filter_rows                                                 &
     &, diff_order                                                      &
     &, diff_timescale                                                  &
     &, u_filt_rows                                                     &
     &, v_filt_rows                                                     &
     &, half_grows                                                      &
     &, info

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! 1.0 Set filtering/diffusion parameters
! ----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('SETDIFF_OLD',zhook_in,zhook_handle)
        tiny = 0.000001

      if ( L_pfincs .or. L_pftheta .or. L_pfuv .or. L_pfw )             &
     &    L_pfcomb = .true.

      if ( L_diff_thermo .or. L_diff_wind ) then
        L_filter = .true.
!  Only 1st order diffusion allowed for combi diffusion
        diff_order = 1
        diff_timescale = max(diff_timescale_thermo,diff_timescale_wind)
        If ( diff_timescale < 1 ) Then
          diff_timescale = 1
          write(6,*)' Diffusion timescale was set to < 1 timestep'      &
     &    ,' Did you mean to do this? '
          write(6,*)' Diffusion timescale has been reset to 1 timestep' &
     &    ,' but you might want to change this in the UMUI'
        endIf ! diff_timescale < 1
        If ( diff_coeff_ref < tiny ) Then
          write(6,*)' Reference diffusion coefficient is too small or 0'&
     &    ,' Reset this in the UMUI or switch off diffusion '
        endIf ! diff_coeff_ref < tiny
      endif ! L_diff_thermo .or. L_diff_wind

! ----------------------------------------------------------------------
! 1.2 Calculate diffusion coefficients given e-folding timesteps, order
! ----------------------------------------------------------------------
      Pi_over_2 = 0.5 * Pi
      damp_ratio = 1.0  ! hardwired for current HadGEM1 usage
      d_ratio = (delta_phi * delta_phi) / (delta_lambda * delta_lambda)
      filter_lat = Pi_over_180 * polar_cap
      If ( model_domain  /=  mt_global ) Then
        filter_lat = 0.0
      endIf !  model_domain  /=  mt_global
      cos_polar_cap = cos(filter_lat)
      if( cos_polar_cap > 0.0 ) then
        scale_test = damp_ratio / cos_polar_cap
      else ! cos_polar_cap > 0.0
        scale_test = 1000.0
      endif ! cos_polar_cap > 0.0

! If damp_ratio = 1, then the damping at the polar-cap of the physical
! scale equivalent to 2-grid-length waves equatorwards is the same.
! Increasing damp_ratio reduces the physical scale at the polar-cap
! chosen to match the 2 grid waves and thus increases the diffusion
      sin_ratio =   sin (Pi_over_2 * cos_polar_cap * damp_ratio) *      &
     &                sin (Pi_over_2 * cos_polar_cap * damp_ratio)
      If ( L_diff_thermo ) Then
        If ( diff_order_thermo > 0 ) Then
          recip_diffusion_order = 1.0 / real(diff_order_thermo)
          If ( diff_timescale_thermo > 0 ) Then

            diff_coeff_thermo = diff_coeff_ref *                        &
     &                       (1.0 - EXP( -1.0 / diff_timescale_thermo) )&
     &                                        ** recip_diffusion_order
          else !  diff_timescale_thermo = 0
            diff_coeff_thermo = diff_coeff_ref
          endIf  !  diff_timescale_thermo > 0

! determine EW diffusion coefficients will scaling relative to
!  adjacent row polewards
          half_grows = global_rows / 2
          if( 2 * half_grows < global_rows )then
            L_rows_even = .false.
            half_grows = half_grows + 1
          else !  2 * half_grows = global_rows
            L_rows_even = .true.
          endif !  2 * half_grows < global_rows
          diff_coeff(1) = 1.0
! even though doing S Hem, can use +ve lats since taking cosine values

          lat_new = Pi_over_2
          cos_lat_new = 0.0
          Do j = 2, half_grows
            lat_new = lat_new - delta_phi
            cos_lat_new = cos(lat_new)
            if ( lat_new < filter_lat ) then
              scale = sin ( Pi_over_2 * cos_lat_new )
              diff_coeff(j) = sin_ratio / ( scale * scale )
            else
              diff_coeff(j) = 1.0
            endif  !  lat_new < filter_lat
          EndDo ! j = 2, half_grows
          if( L_rows_even )then
            j_stop_u = half_grows
          else ! L_rows_even
            j_stop_u = half_grows - 1
          endif ! L_rows_even
          Do j = 1, j_stop_u
           diff_coeff(global_rows - j + 1 ) = diff_coeff(j)
          endDo !  j = 1, j_stop_u
!  Now copy appropriate values for this processor
          Do j = 1, rows
            gj = datastart(2) + j - 1
            Do i = 1, row_length
              diff_coeff_u(i,j) = diff_coeff(gj) * diff_coeff_thermo
            End Do
          End Do
!  Diffusion coeff NOT a vector field
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
     &                     diff_coeff_u, row_length, rows, 1,           &
     &                     offx, offy, fld_type_p, .false.)
!  NS Diffusion coeff scaled from equator EW coefficient
          diff_coeff_thermo = diff_coeff(half_grows) * diff_coeff_thermo&
     &                                                * d_ratio
          if ( cos_polar_cap > 0.0 .and. diff_coeff_thermo > 0.0) then
            diff_time_equ =  -1.0 / log( 1.0 -                          &
     &                                     (8.0 * diff_coeff_thermo)    &
     &                                         ** diff_order )
          else
            diff_time_equ =  diff_timescale_thermo
          endif ! cos_polar_cap > 0.0
       if ( diff_coeff_thermo > 0.0) then
       Write ( Unit = 6, fmt=*) 'Diffusion order for theta is '         &
     &                          ,diff_order_thermo,                     &
     &          ' with diffusion coefficient = ',diff_coeff_ref
       Write ( Unit = 6, fmt=*) 'Diffusion timescale for theta is '     &
     &                          , diff_timescale_thermo,' timesteps '
       If ( model_domain == mt_global ) Then
         Write ( Unit = 6, fmt=*) 'This is the time taken to damp the ',&
     &    ' 2 gridlength wave at ',polar_cap,' degrees by 1/e '
         Write ( Unit = 6, fmt=*) 'For theta, North-South diffusion'    &
     &           ,' coefficient and East-West diffusion coefficient'    &
     &           ,' at the equator = ', diff_coeff_thermo
         Write ( Unit = 6, fmt=*) 'Diffusion timescale for theta at ',  &
     &                ' the equator is ',diff_time_equ, ' timesteps '
         i = 2 * nint( scale_test )
        Write ( Unit = 6, fmt=*) 'The damping of the 2 gridlength wave',&
     & ' at the equator is the same as the ',i,'th gridlength wave at'  &
     &                      ,polar_cap,' degrees '
       else!  model_domain  /=  mt_global
         Write ( Unit = 6, fmt=*) 'This is the time taken to damp the ',&
     &    ' 2 gridlength wave by 1/e '
         Write ( Unit = 6, fmt=*) 'For theta, North-South diffusion'    &
     &           ,' coefficient and East-West diffusion coefficient = ' &
     &           , diff_coeff_thermo
       endIf !  model_domain == mt_global
       else ! diff_coeff_thermo = 0.0
       Write ( Unit = 6, fmt=*) 'Diffusion coefficient for theta is '   &
     &                          ,diff_coeff_thermo
       Write ( Unit = 6, fmt=*) 'The reference diffusion coefficient '  &
     &   , ' on the COMBI diffusion page is probably 0.0'
       Write ( Unit = 6, fmt=*) 'You may want to reset to 0.25 '
       endif ! diff_coeff_thermo > 0.0
        Else
          diff_coeff_thermo = 0.0
          diff_timescale_thermo = 0
       Write ( Unit = 6, fmt=*) 'Horizontal diffusion for theta is '    &
     &                          ,'  NOT ACTIVE '
        EndIf ! diff_order_thermo > 0
      Else ! L_diff_thermo = .false.
       Write ( Unit = 6, fmt=*) 'Horizontal diffusion for theta is '    &
     &                          ,'  NOT ACTIVE '
      EndIf ! L_diff_thermo

      If (L_diff_wind) then
        If ( diff_order_wind > 0 ) Then
          recip_diffusion_order = 1.0 / real(diff_order_wind)
          If ( diff_timescale_wind > 0 ) Then
            diff_coeff_wind = diff_coeff_ref *                          &
     &                       (1.0 - EXP( -1.0 / diff_timescale_wind) )  &
     &                                       ** recip_diffusion_order
          else  ! diff_timescale_wind = 0
            diff_coeff_wind = diff_coeff_ref
          endIf ! diff_timescale_wind > 0
! determine where diffusion will begin by finding
!                         how many v-rows are filtered
          half_grows = (global_rows - 1) / 2
          if( 2 * half_grows < global_rows - 1 )then
            L_rows_even = .false.
            half_grows = half_grows + 1
          else !  2 * half_grows = global_rows
            L_rows_even = .true.
          endif !  2 * half_grows < global_rows
          diff_coeff(1) = 1.0
! even though doing S Hem, can use +ve lats since taking cosine values
          lat_new = 0.5 * ( Pi - delta_phi)
          cos_lat_new = cos(lat_new)
          Do j = 2, half_grows
            lat_new = lat_new - delta_phi
            cos_lat_new = cos(lat_new)
            if ( lat_new < filter_lat ) then
              scale = sin ( Pi_over_2 * cos_lat_new )
              diff_coeff(j) = sin_ratio / ( scale * scale )
            else
              diff_coeff(j) = 1.0
            endif  !  lat_new > filter_lat
          EndDo ! j = 2, half_grows
          if( L_rows_even )then
            j_stop_v = half_grows
          else ! L_rows_even
            j_stop_v = half_grows - 1
          endif ! L_rows_even
          Do j = 1, j_stop_v
           diff_coeff(global_rows - j) = diff_coeff(j)
          endDo !  j = 1, j_stop_v
!  Now copy appropriate values for this processor
          Do j = 1, n_rows
            gj = datastart(2) + j - 1
            Do i = 1, row_length
              diff_coeff_v(i,j) = diff_coeff(gj) * diff_coeff_wind
            End Do
          End Do
!  Diffusion coeff NOT a vector field
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
     &                     diff_coeff_v, row_length, n_rows, 1,         &
     &                     offx, offy, fld_type_v, .false.)
!  NS Diffusion coeff scaling of equator EW coefficient
          diff_coeff_wind = diff_coeff(half_grows) * diff_coeff_wind    &
     &                                              * d_ratio
          if ( cos_polar_cap > 0.0 .and. diff_coeff_wind > 0.0) then
            diff_time_equ =  -1.0 / log( 1.0 -                          &
     &                                          (8.0 * diff_coeff_wind) &
     &                                              ** diff_order)
          else
            diff_time_equ =  diff_timescale_thermo
          endif ! cos_polar_cap > 0.0
       if ( diff_coeff_wind > 0.0 ) then
       Write ( Unit = 6, fmt=*) 'Diffusion order for wind is '          &
     &                          ,diff_order_wind,                       &
     &          ' with diffusion coefficient = ',diff_coeff_ref
       Write ( Unit = 6, fmt=*) 'Diffusion timescale for wind is '      &
     &                          , diff_timescale_wind,' timesteps '
       If ( model_domain == mt_global ) Then
         Write ( Unit = 6, fmt=*) 'This is the time taken to damp the ',&
     &    ' 2 gridlength wave at ',polar_cap,' degrees by 1/e '
         Write ( Unit = 6, fmt=*) 'For u and v, North-South diffusion'  &
     &   ,' coefficient and East-West diffusion coefficient '           &
     &   ,' at the equator = ', diff_coeff_wind
        Write ( Unit = 6, fmt=*) 'Diffusion timescale for wind at the ',&
     & ' equator is ',diff_time_equ, ' timesteps '
       Else !  model_domain  /=  mt_global
         Write ( Unit = 6, fmt=*) 'This is the time taken to damp the ',&
     &    ' 2 gridlength wave by 1/e '
         Write ( Unit = 6, fmt=*) 'For u and v, North-South diffusion ' &
     &             ,' coefficient and East-West diffusion coefficient ' &
     &             , diff_coeff_wind
       EndIf ! model_domain == mt_global
       else ! diff_coeff_wind = 0.0
       Write ( Unit = 6, fmt=*) 'Diffusion coefficient for wind is '    &
     &                          ,diff_coeff_wind
       Write ( Unit = 6, fmt=*) 'The reference diffusion coefficient '  &
     &   , ' on the COMBI diffusion page is probably 0.0'
       Write ( Unit = 6, fmt=*) 'You may want to reset to 0.25 '
       endif ! diff_coeff_wind > 0.0
        Else
          diff_coeff_wind = 0.0
          diff_timescale_wind = 0
       Write ( Unit = 6, fmt=*) 'Horizontal diffusion for winds is '    &
     &                          ,'  NOT ACTIVE '
        EndIf ! diff_order_wind > 0
      Else ! L_diff__wind = .false
       Write ( Unit = 6, fmt=*) 'Horizontal diffusion for winds is '    &
     &                          ,'  NOT ACTIVE '
      EndIf ! L_diff__wind
      if ( L_diff_incs ) then
        Write ( Unit = 6, fmt=*) 'L_diff_incs = ',L_diff_incs,          &
     &       ' so increments will be horizontally diffused.'
      endif ! L_diff_incs
      if ( L_diff_w ) then
        Write ( Unit = 6, fmt=*) 'L_diff_w = ',L_diff_w,                &
     &       ' so w will be horizontally diffused like theta'
      endif ! L_diff_w

! ----------------------------------------------------------------------
! 1.3 Initialise scanning structure for 1-2-1 filter
! ----------------------------------------------------------------------

      max_filter_rows = max_121_rows

!  Initialise sweep boundaries for 1-2-1 filter
        Do j = 0, max_filter_rows
          u_begin(j) = -1
          v_begin(j) = -1
          u_end(j) = -2
          v_end(j) = -2
        End Do ! j = 0, max_filter_rows
        u_filt_rows = 0
        v_filt_rows = 0
        Do j = 1, max_filter_rows
          u_sweeps(j) = 0
          v_sweeps(j) = 0
        End Do ! j = 0, max_filter_rows

      If (model_domain  ==  mt_global ) then
        If ( L_pftheta .or. L_pfuv .or. L_pfw ) then

        scale_r2 = scale_ratio * scale_ratio

!  Calculate sweeping structure required near poles

        j_start_u = 1
        j_stop_u = rows
        If (at_extremity(PSouth)) j_start_u = 2
        If (at_extremity(PNorth)) j_stop_u = rows - 1

        if( sin_theta_latitude(1,j_start_u) < -0.5 )then
!  This processor starts in S. Hemisphere and only extends to N. Hem
!   filtering area when nprocy=1.  Find filtering area (if it exists)
          u_end(0) = j_stop_u
          j = j_start_u
          lat_row = acos( cos_theta_latitude(1,j) )
          u_filt_rows = int( (lat_row - filter_lat) / delta_phi )
          if(u_filt_rows < 1 )then
            u_filt_rows = 0
            u_begin(0) = j
          else
            u_end(u_filt_rows) = j
            Do  i = 1 ,u_filt_rows
              u_begin(i) = j
            EndDo ! i = 1 ,max_filter_rows
            Do i = u_filt_rows-1, 1, -1
              j = j + 1
              u_end(i) = j
            EndDo  ! i = u_filt_rows-1, 1, -1
            u_begin(0) = u_end(1) + 1
! If 1-2-1 is active, always do at least 1 sweep
            u_sweeps(1) = 1
            Do i = 2, u_filt_rows
              scale_test = cos_theta_latitude(1,u_end(i)) /             &
     &                      cos_theta_latitude(1,u_end(i-1))
              if (scale_test < scale_ratio * scale_r2 ) then
                u_sweeps(i) = 3
              elseif (scale_test < scale_r2 ) then
                u_sweeps(i) = 2
              else
                u_sweeps(i) = 1
              endif ! scale_test < scale_ratio * scale_r2
            EndDo  ! i = 2, u_filt_rows
          endif ! u_filt_rows < 1
        else  !  This processor may need filtering in N. Hemisphere
          u_begin(0) = j_start_u
          j = j_stop_u
          lat_row = acos( cos_theta_latitude(1,j) )
          u_filt_rows = int( (lat_row - filter_lat) / delta_phi )
          if(u_filt_rows < 1 )then
            u_filt_rows = 0
            u_end(0) = j
          else
            u_begin(u_filt_rows) = j
            Do  i = 1 ,u_filt_rows
              u_end(i) = j
            EndDo ! i = 1 ,max_filter_rows
            Do i = u_filt_rows-1, 1, -1
              j = j - 1
              u_begin(i) = j
            EndDo  ! i = u_filt_rows-1, 1, -1
            u_end(0) = u_begin(1) - 1
! If 1-2-1 is active, always do at least 1 sweep
            u_sweeps(1) = 1
            Do i = 2, u_filt_rows
              scale_test = cos_theta_latitude(1,u_begin(i)) /           &
     &                      cos_theta_latitude(1,u_begin(i-1))
              if (scale_test < scale_ratio * scale_r2 ) then
                u_sweeps(i) = 3
              elseif (scale_test < scale_r2 ) then
                u_sweeps(i) = 2
              else
                u_sweeps(i) = 1
              endif ! scale_test < scale_ratio * scale_r2
            EndDo  ! i = 2, u_filt_rows
          endif ! u_filt_rows < 1
!  End for processors which may need filtering in N. Hemisphere
        Endif ! sin_theta_latitude(1,j_start_u) < -0.5

!  Do same as above for v rows
        if( sin_v_latitude(1,j_start_u) < -0.5 )then
!  This processor starts in S. Hemisphere and only extends to N. Hem
!   filtering area when nprocy=1.  Find filtering area (if it exists)
          v_end(0) = n_rows
          j = 1
          lat_row = acos( cos_v_latitude(1,j) )
          v_filt_rows = int( (lat_row - filter_lat) / delta_phi )
          if(v_filt_rows < 1 )then
            v_filt_rows = 0
            v_begin(0) = j
          else
            v_end(v_filt_rows) = j
            Do  i = 1 ,v_filt_rows
              v_begin(i) = j
            EndDo ! i = 1 ,max_filter_rows
            Do i = v_filt_rows-1, 1, -1
              j = j + 1
              v_end(i) = j
            EndDo  ! i = v_filt_rows-1, 1, -1
            v_begin(0) = v_end(1) + 1
! If 1-2-1 is active, always do at least 1 sweep
            v_sweeps(1) = 1
            Do i = 2, v_filt_rows
              scale_test = cos_v_latitude(1,v_end(i)) /                 &
     &                      cos_v_latitude(1,v_end(i-1))
              if (scale_test < scale_ratio * scale_r2 ) then
                v_sweeps(i) = 3
              elseif (scale_test < scale_r2 ) then
                v_sweeps(i) = 2
              else
                v_sweeps(i) = 1
              endif ! scale_test < scale_ratio * scale_r2
            EndDo  ! i = 2, v_filt_rows
          endif ! v_filt_rows < 1
        else  !  This processor may need filtering in N. Hemisphere
          v_begin(0) = 1
          j = n_rows
          lat_row = acos( cos_v_latitude(1,j) )
          v_filt_rows = int( (lat_row - filter_lat) / delta_phi )
          if(v_filt_rows < 1 )then
            v_filt_rows = 0
            v_end(0) = j
          else
            v_begin(v_filt_rows) = j
            Do  i = 1 ,v_filt_rows
              v_end(i) = j
            EndDo ! i = 1 ,max_filter_rows
            Do i = v_filt_rows-1, 1, -1
              j = j - 1
              v_begin(i) = j
            EndDo  ! i = v_filt_rows-1, 1, -1
            v_end(0) = v_begin(1) - 1
! If 1-2-1 is active, always do at least 1 sweep
            v_sweeps(1) = 1
            Do i = 2, v_filt_rows
              scale_test = cos_v_latitude(1,v_begin(i)) /               &
     &                      cos_v_latitude(1,v_begin(i-1))
              if (scale_test < scale_ratio * scale_r2 ) then
                v_sweeps(i) = 3
              elseif (scale_test < scale_r2 ) then
                v_sweeps(i) = 2
              else
                v_sweeps(i) = 1
              endif ! scale_test < scale_ratio * scale_r2
            EndDo  ! i = 2, v_filt_rows
          endif ! v_filt_rows < 1
!  End for processors which may need filtering in N. Hemisphere
        Endif ! sin_v_latitude(1,j) < -0.5

!  max values for filter sweeps need to be given to all processors
! to synchronize swap bounds. If sweep is not needed on a processor
! then the active loops will no-op due to loop bounds settings
        global_u_filter = u_filt_rows
        global_v_filter = v_filt_rows
        call gc_imax(1, nproc, info, global_u_filter)
        call gc_imax(1, nproc, info, global_v_filter)
        call gc_imax(global_u_filter, nproc, info, u_sweeps)
        call gc_imax(global_v_filter, nproc, info, v_sweeps)

        if( global_v_filter > max_filter_rows ) then
          Write ( Unit = 6, fmt=*) '*** WARNING *** '
          Write ( Unit = 6, fmt=*) 'You need to increase '              &
     &     ,' the PARAMETER max_121_rows in include file CMAXSIZE'      &
     &     ,' to ',global_v_filter
          Write ( Unit = 6, fmt=*)                                      &
     &         'OR set polar_cap closer to 90 degrees.'
        else !
          Write ( Unit = 6, fmt=*) 'This is processor ',mype
          Write ( Unit = 6, fmt=*) 'u_filt_rows = ',u_filt_rows         &
     &                          , ' v_filt_rows = ',v_filt_rows
          Write ( Unit = 6, fmt=*) 'global_u_filter = ',global_u_filter &
     &                          , ' global_v_filter = ',global_v_filter
          i = 0
          n_sweeps = 0
          Do
            if ( i > 0) then
              if ( u_sweeps(i) > 2 ) then
                n_sweeps = n_sweeps + 3
              elseif ( u_sweeps(i) > 1 ) then
                n_sweeps = n_sweeps + 2
              else
                n_sweeps = n_sweeps + 1
              endif  ! u_sweeps(i) > 2
            endif ! i > 0
            Write ( Unit = 6, fmt=*) 'Sweep ',n_sweeps,                 &
     &             ' u_begin = ', u_begin(i),' u_end = ',u_end(i)
            i = i + 1
            if ( i > u_filt_rows) EXIT
          CYCLE
          endDo
          i = 0
          n_sweeps = 0
          Do
            if ( i > 0) then
              if ( v_sweeps(i) > 2 ) then
                n_sweeps = n_sweeps + 3
              elseif ( v_sweeps(i) > 1 ) then
                n_sweeps = n_sweeps + 2
              else
                n_sweeps = n_sweeps + 1
              endif  ! v_sweeps(i) > 2
            endif ! i > 0
            Write ( Unit = 6, fmt=*) 'Sweep ',n_sweeps,                 &
     &             ' v_begin = ', v_begin(i),' v_end = ',v_end(i)
            i = i + 1
            if ( i > v_filt_rows) EXIT
          CYCLE
          endDo
          Write ( Unit = 6, fmt=*) 'Polar filter activates polewards '  &
     &   , ' of ',polar_cap, ' degrees and needs a maximum of '         &
     &   , n_sweeps, ' sweeps of the 1-2-1 filter '
        endif !  global_v_filter > max_filter_rows

        EndIf  !   L_pftheta .or. L_pfuv .or. L_pfw
      EndIf  !   model_domain  ==  mt_global

!  LAM domains can have diffusion so ...
      If  (model_domain  /=  mt_global) then

        j_start_u = 1
        j_stop_u = rows
        j_start_v = 1
        j_stop_v = n_rows
        If (at_extremity(PSouth)) then
         j_start_u = 2
         j_start_v = 2
        EndIf ! at_extremity(PSouth)
        If (at_extremity(PNorth)) then
          j_stop_u = rows - 1
          j_stop_v = n_rows - 1
        EndIf ! at_extremity(PNorth)
! Now set loop bounds for one sweep of diffusion
        u_begin(0) = j_start_u
        u_end(0) = j_stop_u
        v_begin(0) = j_start_v
        v_end(0) = j_stop_v
       Write ( Unit = 6, fmt=*) 'u_begin(0) = ',u_begin(0)              &
     &                         , ' u_end(0) = ',u_end(0)
       Write ( Unit = 6, fmt=*) 'v_begin(0) = ',v_begin(0)              &
     &                         , ' v_end(0) = ',v_end(0)

      EndIf ! model_domain  /=  mt_global

      IF (lhook) CALL dr_hook('SETDIFF_OLD',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE Setdiff_old
