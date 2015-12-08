! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Setdiff
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
!
      SUBROUTINE eg_setdiff(                                            &
! Input size and control variables
     &      global_rows, row_length, rows, n_rows, model_levels,        &
     &      model_domain, at_extremity, datastart,                      &
     &      offx, offy, mype, nproc_y,                                  &
     &      max_model_levels, max_121_rows, max_sweeps,                 &
! other constants
     &      delta_lambda, delta_phi,                                    &
     &      polar_cap, scale_ratio,                                     &
     &      ref_lat_deg, diff_coeff_ref,                                &
     &      cos_theta_latitude, sin_theta_latitude,                     &
     &      cos_v_latitude, sin_v_latitude,                             &
! Output data
     &       global_u_filter, global_v_filter,                          &
     &       u_sweeps, v_sweeps,                                        &
     &       u_begin, u_end, v_begin, v_end,                            &
     &       diff_coeff_u, diff_coeff_v, diff_coeff_phi,                &
     &       diffusion_coefficient_thermo, diffusion_order_thermo,      &
     &       diffusion_coefficient_wind, diffusion_order_wind,          &
     &       diff_order_thermo, diff_order_wind,                        &
     &       diff_timescale_thermo, diff_timescale_wind,                &
     &       L_sponge, sponge_power, sponge_ew, sponge_ns,              &
     &       sponge_wts_ew, sponge_wts_ns, L_subfilter_horiz,           &
     &       L_diffusion, L_cdiffusion, L_filter, L_diff_auto,          &
     &       L_pfcomb, L_pfexner, L_pftheta, L_pfuv, L_pfw, L_pfincs,   &
     &       L_diff_thermo, L_diff_wind, L_diff_w, L_diff_incs,         &
     &       halo_i, halo_j,xi2_p, xi2_v                              )


      USE turb_diff_ctl_mod, ONLY: turb_ustart, turb_uend,              &
                                   turb_vstart, turb_vend,              & 
                                   turb_phi_st, turb_phi_end
      USE conversions_mod, ONLY: pi_over_180,pi
      USE ereport_mod, ONLY: ereport
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
     &       model_levels,                                              &
     &       offx,                                                      &
     &       offy,                                                      &
     &       mype,                                                      &
     &       nproc_y,                                                   &
     &       model_domain                                               &
     &, max_model_levels                                                &
                          ! Max no. model levels
     &, max_121_rows                                                    &
                      ! Max no. of rows 1-2-1 filtered in a hemisphere
     &, max_sweeps                                                      &
                      ! Max no. of sweeps of 1-2-1 filter
     &, diff_timescale_thermo                                           &
     &, diff_timescale_wind                                             &
     &, diff_order_thermo                                               &
     &, diff_order_wind                                                 &
     &, datastart(3)                                                    &
                           ! First gridpoints held by this processor
     &, diffusion_order_thermo(max_model_levels)                        &
     &, diffusion_order_wind(max_model_levels)                          &
     &, sponge_power                                                    &
     &, sponge_ew                                                       &
     &, sponge_ns

      REAL, Intent(In) ::                                               &
     &  delta_lambda                                                    &
                     ! EW (x) grid spacing
     &, delta_phi                                                       &
                     ! NS (y) grid spacing
     &, scale_ratio                                                     &
     &, diffusion_coefficient_thermo(max_model_levels)                  &
     &, diffusion_coefficient_wind(max_model_levels)

      LOGICAL, Intent(In) ::                                            &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_subfilter_horiz                                               &
                    ! switch for horiz turb diffusion
     &, L_diffusion                                                     &
                    ! switch for old diffusion
     &, L_cdiffusion                                                    &
                     ! switch for conservative diffusion
     &, L_pfincs                                                        &
                 ! switch for polar filter
     &, L_pftheta                                                       &
                       ! switch for polar filter for theta
     &, L_pfw                                                           &
                       ! switch for polar filter for w
     &, L_pfuv                                                          &
                       ! switch for polar filter for horizontal winds
     &, L_pfexner                                                       &
                       ! switch for polar filter for Exner pressure
     &, L_diff_w                                                        &
                       ! switch for horiz. diffusion of w
     &, L_diff_auto                                                     &
                       ! T if automatic calculation of diffusion parms
     &, L_sponge       ! T if sponge zone active at lateral boundaries

      Real, Intent(In) ::                                               &
     &  sin_theta_latitude (row_length, rows)                           &
     &, sin_v_latitude (row_length, n_rows)                             &
     &, cos_theta_latitude(1-offx:row_length+offx,                      &
     &                     1-offy:rows+offy )                           &
     &, cos_v_latitude(1-offx:row_length+offx, 1-offy:n_rows+offy )

      INTEGER, INTENT(IN) :: halo_i, halo_j

      Real, Intent(In) ::                                               &
        xi2_p(1-halo_j:rows+halo_j),                                    &
        xi2_v(-halo_j:n_rows-1+halo_j)
     

  ! Output arguments for diffusion/filtering control

      LOGICAL, Intent(InOut) ::                                         &
     &  L_filter                                                        &
                 ! switch for polar filter
     &, L_pfcomb                                                        &
                       ! switch for combined polar filter/diffusion
     &, L_diff_thermo                                                   &
                       ! switch for horiz. diffusion of theta
     &, L_diff_wind                                                     &  
                       ! switch for horiz. diffusion of u,v
     &, L_diff_incs                                                     
                       ! switch for horiz. diffusion of incs

      INTEGER, Intent(InOut) ::                                         &
     &  global_u_filter                                                 &
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
     &, ref_lat_deg                                                     &
     &, diff_coeff_phi                                                  &
                         ! NS diffusion coeff
!    diffusion coefficient for u rows
     &, diff_coeff_u(1-offx:row_length+offx, 1-offy:rows+offy)          &
!    diffusion coefficient for v rows
     &, diff_coeff_v(1-offx:row_length+offx, 1-offy:n_rows+offy)        &
     &, sponge_wts_ew(sponge_ew)                                        &
     &, sponge_wts_ns(sponge_ns)


! Local variables
      REAL                                                              &
     &  angle                                                           &
     &, c1, c2                                                          &
     &, cos_lat1                                                        &
     &, cos_lat2                                                        &
     &, cos_ref_lat                                                     &
     &, denom                                                           &
     &, diff_new                                                        &
     &, diff_time_equ                                                   &
     &, diff_time_print                                                 &
     &, diff_coeff_thprint                                              &
     &, diff_coeff_uvprint                                              &
     &, diff_coeff_temp                                                 &
     &, dphi_half                                                       &
     &, exp1                                                            &
     &, filter_lat                                                      &
     &, first_lat                                                       &
     &, lat_new                                                         &
     &, Pi_over_2                                                       &
     &, pioverr                                                         &
     &, ref_lat_rad                                                     &
     &, rpower                                                          &
     &, scale_test                                                      &
     &, small                                                           &
     &, tiny                                                            &
     &, weight

! Local arrays
      REAL                                                              &
     & diff_scale(0:global_rows)

      INTEGER                                                           &
     &  diff_order                                                      &
     &, diff_start_row                                                  &
     &, diff_timescale                                                  &
     &, filt_row                                                        &
     &, filt_NP                                                         &
     &, filt_SP                                                         &
     &, first_row                                                       &
     &, last_row                                                        &
     &, global_filter                                                   &
     &, half_grows                                                      &
     &, n_sweeps                                                        &
     &, power

      INTEGER                                                           &
     &  sweeps(2 * max_121_rows)

      LOGICAL                                                           &
     &  L_rows_even                                                     &
     &, L_diff_message                                                  &
     &, L_test

      INTEGER                                                           &
              ! Mostly loop counters, but j also used for interp points
     &  I, J, GJ, k                                                     &
     &, j_start_u, j_stop_u                                             &
     &, j_start_v, j_stop_v                                             &
     &, j_begin, j_end, j_stop                                          &
      , js, je, icode
     
      CHARACTER(LEN=*),PARAMETER :: RoutineName = 'setdiff_4A'
      CHARACTER(LEN=80)          :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! 1.0 Set filtering/diffusion parameters
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('EG_SETDIFF',zhook_in,zhook_handle)
      L_filter = .false.
 
      if( L_pfincs .or. L_pftheta .or. L_pfuv .or. L_pfw .or. L_pfexner)&
     &    L_pfcomb = .true.

      If ( model_domain  /=  mt_global ) Then
        if ( L_pfcomb) then
          write(6,*) '  '
          write(6,*) ' POLAR filter is not allowed. Switch OFF in UMUI '
        endif ! L_pfcomb
        filter_lat = 0.0
      endIf !  model_domain  /=  mt_global

      turb_ustart = 1
      turb_uend = rows
      turb_vstart = 1
      turb_vend = n_rows
      turb_phi_st = 1
      turb_phi_end = rows

      Do j = 0, max_121_rows
        u_begin(j) = -1
        v_begin(j) = -1
        u_end(j) = -3
        v_end(j) = -3
      End Do ! j = 0, max_121_rows
      Do j = 1, max_121_rows
        u_sweeps(j) = 0
        v_sweeps(j) = 0
      End Do ! j = 0, max_121_rows

!  Main block for setting filtering and diffusion for global domain
      If ( model_domain == mt_global ) Then
        tiny = 0.000001
        small = 0.001
        Pi_over_2 = 0.5 * Pi
        pioverr = pi / scale_ratio
        ref_lat_rad = ref_lat_deg * Pi_over_180
        cos_ref_lat = cos(ref_lat_rad)
        dphi_half = 0.5 * delta_phi
        filter_lat = Pi_over_180 * polar_cap
        global_filter = 0

!  Initialise sweep boundaries for filtering and diffusion
        j_start_u = 1
        j_stop_u = rows
        j_start_v = 1
        j_stop_v = n_rows
        If (at_extremity(PSouth)) then
          j_start_v = 2
        endIf ! at_extremity(PSouth))
        If (at_extremity(PNorth)) then
          j_stop_v = n_rows - 1
        endIf ! at_extremity(PNorth)

!  Diffusion coefficients calculated for 1 hemisphere only, so ....
!  .... set pointers needed for switching hemispheres
        half_grows = global_rows / 2
        if( 2 * half_grows < global_rows )then
          L_rows_even = .false.
        else !  2 * half_grows = global_rows
          L_rows_even = .true.
        endif !  2 * half_grows < global_rows

! Set diffusion coefficient = 0.25 on all rows . Loop limits
!   in filtering routines will determine whether they are used or not.
!   Values for diffusion will be calculated and overwritten later
          Do j = 1, rows + offy ! extend to halo - TomM
            Do i = 1-offx, row_length
              diff_coeff_u(i,j) = 0.25
            End Do
          EndDo !    j = 1, rows
          Do j = 1, n_rows + offy ! extend to halo - TomM
            Do i = 1-offx, row_length
              diff_coeff_v(i,j) = 0.25
            End Do
          EndDo !    j = 1, n_rows

! ----------------------------------------------------------------------
! 1.2 Initialise scanning structure for 1-2-1 filter
!  Do the calculations on all processors to avoid communication
!  NB ... thus cannot use pre-calculated trigs since polar values
!  might not be on this processor
! ----------------------------------------------------------------------
        if ( L_pfcomb ) then
          L_filter = .true.
          diff_coeff_temp = 0.25
        else ! no polar filter wanted
! calculate diffusion coefficients for all rows if no polar filter
          diff_coeff_temp = diff_coeff_ref
        endif  !   L_pfcomb

        if (L_filter) then
!  Use positive latitudes (i.e. N Pole, S Pole is mirror image)
          first_lat = Pi_over_2 - dphi_half
          filter_lat = Pi_over_180 * polar_cap
! If polar filter has been switched on make sure the end latitude
!   after the row next to the pole
          if( filter_lat > first_lat ) then
            cos_lat1  = cos(first_lat)
            c1 = cos(pioverr * cos_lat1 / cos_ref_lat)
            lat_new = first_lat - dphi_half
            cos_lat2  = cos(lat_new)
            c2 = cos(pioverr * cos_lat2 / cos_ref_lat)
            diff_new = c2 - c1*c1
            sweeps(1) = 1

            if ( diff_new > 0 ) then
              j = 2
              filter_lat  = first_lat
            else ! diff_new < 0
              sweeps(2) = 2
              j = 3
              do    !  cycle loop
                j = j + 1
                c1 = c2
                lat_new = lat_new - dphi_half
                cos_lat2  = cos(lat_new)
                c2 = cos(pioverr * cos_lat2 / cos_ref_lat)
                diff_new = c2 - c1 * c1
                if( diff_new > 0 ) exit
              end do  !  cycle loop
            endif ! diff_new > 0

            filter_lat  =  lat_new
            angle = filter_lat / Pi_over_180
            global_filter = j - 1

            write(6,*) 'global_filter = ', global_filter
            if ( global_filter > 2 * max_121_rows) then
              cmessage = 'max_121_rows must be increased in CMAXSIZE '
              icode = 1
              CALL ereport(RoutineName, icode, cmessage)
            endif ! global_filter > 2 * max_121_rows

            polar_cap = filter_lat / Pi_over_180
            write(6,912) polar_cap

! Determine the number of additional sweeps for each row
            power = 2
! First 2 rows done
            lat_new = filter_lat + dphi_half
            Do i = 3, global_filter
              c2 = c1
              lat_new = lat_new + dphi_half
              cos_lat1  = cos(lat_new)
              c1 = cos(pioverr * cos_lat1 / cos_ref_lat)
              do   !  cycle loop
                power = power + 1
                if( power >= max_sweeps ) exit
                scale_test = c1 ** power
                if( c2 > scale_test ) exit
              enddo !   cycle loop
              sweeps(i) = power
              c1 = scale_test
            EndDo  ! i = 3, global_filter

          else  ! filter_lat set by user

            lat_new = first_lat
            global_filter  = 0
! Set number of sweeps for each row
            write(6,*) ' Setting polar filtering and diffusion settings'
            write(6,*) '     ****  ASSUMING V-AT-THE-POLES   ****'
            write(6,920) polar_cap, first_lat / Pi_over_180
            write(6,921) polar_cap
            write(6,*) ' Set polar_cap = 90.0 in UMUI'                  &
     &             , ' if filtering is to be set automatically.'
            power = 0
            do
              if ( lat_new < filter_lat ) EXIT
              global_filter = global_filter + 1
              if ( global_filter > 2 * max_121_rows) then
                cmessage = 'max_121_rows must be increased in CMAXSIZE '
                icode = 1
                CALL ereport(RoutineName, icode, cmessage)
              else
                power = power + 1
                sweeps(global_filter) = power
                lat_new  = lat_new - dphi_half
              endif ! global_filter > 2 * max_121_rows
            enddo !   while ( diff_new < 0.0 )
            write(6,*) 'global_filter = ', global_filter
            filter_lat  = first_lat - (global_filter - 1.0) * dphi_half
            polar_cap = filter_lat / Pi_over_180
            write(6,912) polar_cap

          endif  ! filter_lat > first_lat

! Now set required values for u and v rows
          global_v_filter = global_filter / 2
          global_u_filter = global_filter - global_v_filter

! sweeps are cumulative to set as additional sweeps
          if ( global_filter == 2 * global_v_filter ) then
            v_sweeps(1) = sweeps(1)
            u_sweeps(1) = sweeps(2)
            Do  i = 2, global_v_filter
              v_sweeps(i) = sweeps(2*i - 1) - v_sweeps(i-1)
            EndDo ! i = 2, global_v_filter
            Do  i = 2, global_u_filter
             u_sweeps(i) = sweeps(2*i) - u_sweeps(i-1)
            EndDo ! i = 2, global_u_filter
          elseif  ( global_v_filter > 0 ) then
            u_sweeps(1) = sweeps(1 )
            v_sweeps(1) = sweeps(2)
            Do  i = 2, global_u_filter
              u_sweeps(i) = sweeps(2*i - 1) - u_sweeps(i-1)
            EndDo ! i = 1, global_u_filter
            Do  i = 2, global_v_filter
              v_sweeps(i) = sweeps(2*i) - v_sweeps(i-1)
            EndDo ! i = 1, global_v_filter
          else  !
            u_sweeps(1) = sweeps(1)
          endif ! global_filter == 2 * global_v_filter

!   Set loop limits on active processors.

          first_row = datastart(2)
          last_row = first_row + rows - 1
          filt_NP = global_rows - global_u_filter + 1

          If ( first_row <= global_u_filter ) Then
! S Hem rows need filtering
 
! Are we filtering across multiple processors ?
            If ( last_row >= global_u_filter ) Then
! northern-most processor row being filtered, set loop limits
              je = global_u_filter - first_row + 2
              j = 1
            Else ! last_row < global_u_filter
! other processor rows being filtered
              je = rows
              Do  j = 1, global_u_filter - last_row + 1
                u_begin(j) = 1
                u_end(j) = je
              End Do ! j = 1, global_u_filter - last_row + 1
              j = global_u_filter - last_row + 2
            EndIf ! last_row > global_u_filter

            Do  ! cycle over rest of filter rows
              je = je - 1
              If (je < 1 ) EXIT
              u_begin(j) = 1
              u_end(j) = je
              j = j + 1
            End Do  ! cycle over filter rows

          Else If ( last_row >= filt_NP) Then
! N Hem rows need filtering
 
! Are we filtering across multiple processors ?
            If ( first_row <= filt_NP ) Then
! southern-most processor row being filtered, set loop limits
              js = filt_NP - first_row
              j = 1
            Else ! first_row > filt_NP
! other processor rows being filtered
              js = 1
              Do  j = 1, first_row - filt_NP + 1
                u_begin(j) = js
                u_end(j) = rows
              End Do ! j = 1, first_row - filt_NP + 1
              j = first_row - filt_NP + 2

            EndIf ! last_row > global_u_filter
 
            Do  ! cycle over rest of filter rows
              js = js + 1
              If (js > rows ) EXIT
              u_begin(j) = js
              u_end(j) = rows
              j = j + 1
            End Do  ! cycle over filter rows

          End If ! first_row <= global_u_filter

!  Repeat for v rows

          last_row = first_row + rows - 1
          filt_SP = global_v_filter + 1
          filt_NP = global_rows - global_v_filter + 1

          If ( first_row <= filt_SP ) Then
! S Hem rows need filtering (but not at pole)
            js = 1
            If ( first_row == 1 ) js = 2

! Are we filtering across multiple processors ?
            If ( last_row >= filt_SP ) Then
! northern-most processor row being filtered, set loop limits

              je = filt_SP - first_row + 2
              j = 1

            Else ! last_row < filt_SP
! other processor rows being filtered

              je = n_rows
              Do  j = 1, filt_SP - last_row + 1
                v_begin(j) = js
                v_end(j) = je
              End Do ! j = 1, filt_SP - last_row + 1
              j = filt_SP - last_row + 2

            EndIf ! last_row > filt_SP

            Do  ! cycle over rest of filter rows
              je = je - 1
              If (je < js ) EXIT
              v_begin(j) = js
              v_end(j) = je
              j = j + 1
            End Do  ! cycle over filter rows

          Else If ( last_row >= filt_NP) Then
! N Hem rows need filtering
            je = rows
 
! Are we filtering across multiple processors ?
            If ( first_row <= filt_NP ) Then
! southern-most processor row being filtered, set loop limits

              js = filt_NP - first_row
              j = 1
              
            Else ! first_row > filt_NP
! other processor rows being filtered

              js = 1
              Do  j = 1, first_row - filt_NP + 1
                v_begin(j) = js
                v_end(j) = je
              End Do ! j = 1, first_row - filt_NP + 1
              j = first_row - filt_NP + 2

            EndIf ! last_row > global_u_filter
 
            Do  ! cycle over rest of filter rows
              js = js + 1
              If (js > je ) EXIT
              v_begin(j) = js
              v_end(j) = je
              j = j + 1
            End Do  ! cycle over filter rows

          End If ! first_row <= global_u_filter

          write(6,*) '  '
          if( global_u_filter > max_121_rows ) then
            write(6,*) '  ******    WARNING      ******* '
            write(6,*) ' You need to increase the PARAMETER '           &
     &             ,'max_121_rows in include file CMAXSIZE for array '  &
     &             , global_u_filter
            write(6,*)' OR set polar_cap closer to 90 degrees.'
            icode = 1
            cmessage = "max_121_rows needs to be increased or "//       &
                         "polar_cap set closer to 90 degrees"
            CALL ereport(RoutineName, icode, cmessage)

          else !   global_u_filter =< max_121_rows

            write(6,*) ' Polar filter details for processor ', mype
            write(6,*) ' global_u_filter = ', global_u_filter,          &
     &               ' global_v_filter = ', global_v_filter
            Do i = 1, global_u_filter
              write(6,*) ' Pass ', i, ' sweeps = ', u_sweeps(i),        &
     &               ' u_begin = ', u_begin(i),' u_end = ', u_end(i)
            end do
            Do i = 1, global_v_filter
              write(6,*) ' Pass ', i, ' sweeps = ', v_sweeps(i),        &
     &               ' v_begin = ', v_begin(i),' v_end = ', v_end(i)
            end do
            write(6,910) polar_cap
            write(6,911) global_v_filter

          endif !  global_u_filter > max_121_rows

          IF ( L_subfilter_horiz ) THEN
!  Set turb_diff start and end rows if there is polar filtering
            IF ( u_end(1) > 0 .OR. v_end(1) > 0 ) THEN
!  polar filtering on this processor

              IF ( v_begin(1) < 3 ) THEN
!  S. Hem processor
                turb_ustart = v_end(1) + 1
                turb_uend = n_rows
              ELSE   ! v_begin(1) >= 3
!  N. Hem processor
                turb_vstart = 1
                turb_vend = v_begin(1) - 1
              END IF ! v_begin(1) < 3
              IF ( u_begin(1) < 2 ) THEN
!  S. Hem processor
                turb_ustart = u_end(1) + 1
                turb_uend = rows
                turb_phi_st = 2
              ELSE   ! u_begin(1) >= 2
!  N. Hem processor
                turb_ustart = 1
                turb_uend = u_begin(1) - 1
                turb_phi_end = n_rows - 1
              END IF ! u_begin(1) < 2

            END IF ! u_end(1) > 0 .OR. v_end(1) > 0

          END IF ! L_subfilter_horiz

        endif  !   L_filter

      endif !  model_domain == mt_global

! ----------------------------------------------------------------------
! 1.4 Calculate diffusion coefficients given e-folding timesteps, order
! ----------------------------------------------------------------------
      L_diff_message = .true.
      L_filter = .false.
      diff_order = max( diff_order_thermo, diff_order_wind )
      If ( diff_order < 1 ) Then
        diff_timescale = 0
        L_diff_thermo = .false.
        L_diff_wind = .false.
        L_diff_incs = .false.
      endIf ! diff_order < 1

      if ( L_diff_thermo .or. L_diff_wind .or. L_diff_incs  ) then
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
          write(6,*)' Reference diffusion coefficient has been reset to'&
     &    ,' 0.000001 to avoid error message in deriving damping time.'
          diff_coeff_ref = tiny
        endIf ! diff_coeff_ref < tiny
      endif ! L_diff_thermo .or. L_diff_wind

! Determine diffusion coefficients
      if ( L_filter .and. model_domain == mt_global) then

! Initialise values for horizontal diffusion
        do j = 1, global_filter
          diff_scale(j) = 0.25
        end do ! j = 1, global_filter
! Set starting row for diffusion (1st after polar filter)
        diff_start_row = global_filter + 1

        If ( L_diff_auto ) Then

          If ( global_filter > 0 ) then
            diff_start_row = global_filter
           write(6,*)'Horizontal diffusion will be applied equatorwards'&
     &            ,' of ',filter_lat / Pi_over_180,' degrees'
          else ! no polar filter
            diff_start_row = 1
            lat_new = Pi_over_2 - dphi_half
            do while ( lat_new > filter_lat )
              diff_scale(diff_start_row) = 0.0
              lat_new = lat_new - dphi_half
              diff_start_row = diff_start_row + 1
            enddo !  while  lat_new > filter_lat
            filter_lat = lat_new
           write(6,*)'Horizontal diffusion will be applied equatorwards'&
     &            ,' of ',filter_lat / Pi_over_180,' degrees'
            diff_scale(diff_start_row) = diff_coeff_temp
          endIf ! global_filter > 0

          cos_lat1 = cos(filter_lat)
          c1 = sin( pioverr * cos_lat1 / cos_ref_lat )
          lat_new = filter_lat - dphi_half
          j = diff_start_row
! This loop scales relative to ref_lat
          if ( lat_new < ref_lat_rad ) then
            diff_scale(j) = diff_coeff_temp
            cos_lat2 = cos_lat1
        write(6,*)'  row    lat    ref_lat_rad    diffusion coefficient'
          write(6,*)j, lat_new, ref_lat_rad, diff_scale(j)
          else
            write(6,*)'row      lat             diffusion coefficient '
          endif ! lat_new < ref_lat_rad
          do while ( lat_new > ref_lat_rad - small )
            j = j + 1
            if( lat_new < 0.0 ) lat_new = 0.0
            cos_lat2 = cos(lat_new)
            c2 = sin( pioverr * cos_lat2 / cos_ref_lat )
            diff_scale(j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
            write(6,*)j, lat_new, diff_scale(j)
            lat_new = lat_new - dphi_half
          enddo !  while ( lat_new > ref_lat_rad - small )

          c2 = sin(pioverr)
          denom = 1.0 / ( c2 * c2 )
! This loop scales r-grid wave to the same physical scale
! as next row polewards (if there are any rows left)
          if  ( lat_new > 0.0 - small ) then
            write(6,*)'  Now match successive rows '
            write(6,*)'  row    lat      diff_scale '
            do while ( lat_new > 0.0 - small )
              j = j + 1
              cos_lat1 = cos_lat2
              if( lat_new < 0.0 ) lat_new = 0.0
              cos_lat2 = cos(lat_new)
              c1 = sin(pioverr * cos_lat1 /cos_lat2)
              diff_scale(j) = diff_scale(j - 1) * c1 * c1 * denom
              write(6,*)j, lat_new, diff_scale(j)
              lat_new = lat_new - dphi_half
            enddo ! while ( lat_new > 0.0 - small )
          endif  ! lat_new > 0.0 - small

!  NS diffusion coefficient is the scaled relative to
!EW diffusion at equator
          diff_coeff_phi = diff_scale(j) * delta_phi * delta_phi /      &
     &                                     (delta_lambda * delta_lambda)
!  Now copy appropriate values for this processor
          v_begin(0) = -1
          v_end(0) = -2
          if( sin_v_latitude(1,1) + tiny < 0.0 )then
!  First row is in Southern  hemisphere
            if( v_begin(1) > 0 )then
              j_begin = v_begin(1)
            else
              j_begin = j_start_v
            endif ! v_begin(1) > 0
            j_end = j_stop_v
            j_stop = j_end

            IF ( nproc_y == 1 ) j_stop = half_grows + 1

          else ! sin_v_latitude(1,1) + tiny > 0.0
!  First row is in Northern  hemisphere
            if( v_end(1) > 0 )then
              j_end = v_end(1)
            else
              j_end = j_stop_v
            endif ! v_begin(1) > 0
            j_begin = j_start_v
            j_stop = j_end
          endif ! sin_v_latitude(1,1) + tiny < 0.0

          L_test = .false.
          j = j_begin - 1
          Do !  j = j_begin, j_stop
            j = j + 1
            if ( j > j_stop ) exit
            if ( sin_v_latitude(1,j) > 0.0 ) exit
            gj = 2 * ( datastart(2) + j - 2)
            if ( v_begin(0) < 1 .and. diff_scale(gj) > tiny ) then
              v_begin(0) = j
              v_end(0) = j_stop
              L_test = .true.
            endif ! v_begin(0) < 1 .and. diff_scale(gj) > tiny
            if ( L_test .and. diff_scale(gj) < tiny ) then
              v_end(0) = j - 1
              L_test = .false.
            endif ! L_test .and. diff_scale(gj) < tiny
            Do i = 1-offx, row_length

              ! Overwrite using ENDGame xi_2 array
              If ( ABS(xi2_v(j-1)) <  filter_lat ) THEN
                c2 = sin( pioverr * cos(xi2_v(j-1)) / cos_ref_lat )
                diff_coeff_v(i,j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
              ELSE ! Polar filter
                diff_coeff_v(i,j) = 0.25
              END IF

            End Do
          End Do  !  j cycle
!  now do N. Hem values if needed (loop will exit if not)
          j = j - 1  ! reset row counter
          Do !  j = j_stop + 1, j_end
            j = j + 1
            if ( j > j_end ) exit
            gj = 2 * (global_rows - datastart(2) - j + 1)
            if ( v_begin(0) < 1 .and. diff_scale(gj) > tiny ) then
              v_begin(0) = j
              v_end(0) = j_stop_v
              L_test = .true.
            endif ! v_begin(0) < 1 .and. diff_scale(gj) > tiny
            if ( L_test .and. diff_scale(gj) < tiny ) then
              v_end(0) = j - 1
              L_test = .false.
            endif ! L_test .and. diff_scale(gj) < tiny
            Do i = 1-offx, row_length

              If ( ABS(xi2_v(j-1)) <  filter_lat ) THEN
              ! Overwrite using ENDGame xi_2 array
              c2 = sin( pioverr * cos(xi2_v(j-1)) / cos_ref_lat )
              diff_coeff_v(i,j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
              ELSE ! Polar filter
                diff_coeff_v(i,j) = 0.25
              END IF
            End Do
          End Do  !  j cycle

! Reset v_begin(0) and v_end(0) to determine rows for NS diffusion only
          v_begin(0) = j_start_v
          v_end(0) = j_stop_v
          write(6,*)'On processor ', mype
          write(6,*)'v_begin(0) = ', v_begin(0),' v_end(0) = ', v_end(0)
          write(6,*)'v_begin(1) = ', v_begin(1),' v_end(1) = ', v_end(1)

! EW diffusion rows are contained within NS diffusion rows and row
! limits are combined with polar filtering rows

!  when  nproc_y == 1  S hem filtered/diffused first then N.Hem
! but   diff_coeff_v  contains coeff for all rows
          if( nproc_y == 1) then
            v_end(0) = j_stop
          endif ! nproc_y == 1

          u_begin(0) = -1
          u_end(0) = -2
          if( sin_theta_latitude(1,1) + tiny < 0.0 )then
!  First row is in Southern  hemisphere
            if( u_begin(1) > 0 )then
              j_begin = u_begin(1)
            else
              j_begin = 1
            endif ! u_begin(1) > 0
            j_end = rows
            j_stop = j_end
            if( nproc_y == 1 ) then
              if( L_rows_even) then
                j_stop = half_grows 
              else ! global_rows is odd
                j_stop = half_grows + 1
              endif !  L_rows_even
            endif ! nproc_y == 1

          else   !( sin_theta_latitude(1,1) + tiny > 0.0
!  First row is in Northern  hemisphere
            if( u_end(1) > 0 )then
              j_end = u_end(1)
            else
              j_end = rows
            endif ! u_end(1) > 0
            j_begin = 1
            j_stop = j_end
          endif !  sin_theta_latitude(1,1) + tiny < 0.0

          L_test = .false.
          j = j_begin - 1
          Do !  j = j_begin, j_stop
            j = j + 1
            if ( j > j_stop ) exit
            if ( sin_theta_latitude(1,j) > 0.0 ) exit
            gj = 2 * ( datastart(2) + j - 1) - 1
            if ( u_begin(0) < 1 .and. diff_scale(gj) > tiny ) then
              u_begin(0) = j
              u_end(0) = j_stop
              L_test = .true.
            endif ! u_begin(0) < 1 .and. diff_scale(gj) > tiny
            if ( L_test .and. diff_scale(gj) < tiny ) then
              u_end(0) = j - 1
              L_test = .false.
            endif ! L_test .and. diff_scale(gj) < tiny
            Do i = 1-offx, row_length

              If ( ABS(xi2_p(j)) <  filter_lat ) THEN
              ! Overwrite using ENDGame xi_2 array
                c2 = sin( pioverr * cos(xi2_p(j)) / cos_ref_lat )
                diff_coeff_u(i,j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
              ELSE ! Polar filter
                diff_coeff_u(i,j) = 0.25
              END IF

            End Do
          End Do !  j cycle
!  now do N. Hem values if needed (loop will exit if not)
          j = j - 1  ! reset row counter
          Do !  j = j_stop + 1, j_end
            j = j + 1
            if ( j > j_end ) exit
            gj = 2 * (global_rows - datastart(2) - j) + 1
            if ( u_begin(0) < 1 .and. diff_scale(gj) > tiny ) then
              u_begin(0) = j
              u_end(0) = rows
              L_test = .true.
            endif ! u_begin(0) < 1 .and. diff_scale(gj) > tiny
            if ( L_test .and. diff_scale(gj) < tiny ) then
              u_end(0) = j - 1
              L_test = .false.
            endif ! L_test .and. diff_scale(gj) < tiny
            Do i = 1-offx, row_length

              If ( ABS(xi2_p(j)) <  filter_lat ) THEN
              ! Overwrite using ENDGame xi_2 array
                c2 = sin( pioverr * cos(xi2_p(j)) / cos_ref_lat )
                diff_coeff_u(i,j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
              ELSE ! Polar filter
                diff_coeff_u(i,j) = 0.25
              END IF

            End Do
          End Do  !  j cycle

!           WRITE(6,*) ' ' 
!           Do j = 1-offy,rows+offy
!             WRITE(6,*) 'j, diff_u = ', j,diff_coeff_u(11-offx,j)
!           End Do
!           WRITE(6,*) ' ' 
!           Do j = 1-offy,n_rows+offy
!             WRITE(6,*) 'j, diff_v = ', j,diff_coeff_v(11-offx,j)
!           End Do
!           WRITE(6,*) ' ' 
!           Do j = 1,global_rows
!             WRITE(6,*) 'j, diff_scale = ', j,diff_scale(j)
!           End Do
!           WRITE(6,*) ' ' 
!           gj = 1
!           Do j = 1,rows
!             WRITE(6,*) ' v ',j,diff_coeff_v(11-offx,j),diff_scale(gj)
!             gj=gj+1
!             WRITE(6,*) ' u ',j,diff_coeff_u(11-offx,j),diff_scale(gj)
!             gj=gj+1
!           End Do
!           WRITE(6,*) ' v ',j,diff_coeff_v(11-offx,n_rows),diff_scale(gj)
!           WRITE(6,*) ' ' 
!           flush(6)

! Reset u_begin(0) and u_end(0) to determine rows for NS diffusion only
          u_begin(0) = j_start_u
          u_end(0) = j_stop_u
          write(6,*)'u_begin(0) = ', u_begin(0),' u_end(0) = ', u_end(0)
          write(6,*)'u_begin(1) = ', u_begin(1),' u_end(1) = ', u_end(1)

! EW diffusion rows are contained within NS diffusion rows and row
! limits are combined with polar filtering rows

!  when  nproc_y == 1  S hem filtered/diffused first then N.Hem
! but   diff_coeff_u  contains coeff for all rows
          if( nproc_y == 1) then
            u_end(0) = j_stop
          endif ! nproc_y == 1

          diff_time_equ = -1.0 / log ( 1.0 -                            &
     &                          4.0 * diff_scale(global_rows - 1))
          diff_time_print = -1.0 / log ( 1.0 -                          &
     &                        4.0 * diff_scale(diff_start_row + 1))
          if ( mype == 0 ) then
            write(6,*) ' '
            write(6,*) '  ****   Horizontal diffusion is ACTIVE ****'
            if ( 2 * (diff_start_row/2) == diff_start_row ) then
              write(6,918) diff_start_row/2 + 1, diff_time_print
              write(6,919) diff_time_equ
            else !  diffusion starts on v row
              write(6,917) diff_start_row/2 + 1, diff_time_print
              write(6,919) diff_time_equ
            endif ! 2 * (diff_start_row/2) == diff_start_row
          endif ! mype == 0

        else ! L_diff_auto false, using user defined parameters

          L_test = .false.
          exp1 = exp(1.0)
          rpower = -1.0 / diff_timescale
          diff_coeff_temp = 0.25 * (1.0 - exp1**rpower)
          v_begin(0) = j_start_v
          v_end(0) = j_stop_v
          Do j = j_start_v, j_stop_v
            Do i = 1-offx, row_length
              diff_coeff_v(i,j) = diff_coeff_temp
            End Do
          End Do

          u_begin(0) = 1
          u_end(0) = rows
          Do j = 1, rows
            Do i = 1-offx, row_length
              diff_coeff_u(i,j) = diff_coeff_temp
            End Do
          End Do

!  NS Diffusion coeff scaled from equator EW coefficient
          diff_coeff_thprint = diff_coeff_ref
          diff_coeff_uvprint = diff_coeff_ref

          If( L_diff_thermo) then
            diff_time_equ =  -1.0 / log( 1.0 -                          &
     &                                   (8.0 * diff_coeff_ref)         &
     &                                       ** diff_order)
            write(6,*) ' '
            write(6,*) '  ****  Horizontal diffusion of theta ',        &
     &               '    is ACTIVE **** '
            write(6,903) diff_order, diff_coeff_thprint
            write(6,904) diff_timescale
            If ( model_domain == mt_global ) Then
              write(6,905) polar_cap
              write(6,906) diff_coeff_ref
              write(6,907) diff_time_equ
            else!  model_domain  /=  mt_global
              write(6,*)' This is the time taken to damp',              &
     &                 ' the 2 gridlength wave by 1/e '
            endIf !  model_domain == mt_global
          Else ! L_diff_thermo = .false.
            L_diff_message = .true.
            If ( L_diffusion .or. L_cdiffusion ) Then
              k = 1
              do ! k = 1, model_levels but end on EXIT
                if ( diffusion_coefficient_thermo(k) > tiny .and.       &
     &                diffusion_order_thermo(k) > 0 ) then
                  L_diff_message = .false.
                  Exit
                endif
                k = k + 1
                if( k > model_levels) Exit
              end do
            endIf ! L_diffusion .or. L_cdiffusion
            If ( L_diff_message ) Then
              write(6,*) '  '
              write(6,*) ' *** Horizontal diffusion is  NOT ACTIVE ***'
            endIf ! L_diff_message
          EndIf ! L_diff_thermo

          If (L_diff_wind ) then
            diff_time_equ =  -1.0 / log( 1.0 -                          &
     &                                      (8.0 * diff_coeff_ref)      &
     &                                          ** diff_order)
            write(6,*) ' '
           write(6,*)'  ****   Horizontal diffusion of horizontal wind',&
     &                 ' is ACTIVE ****'
            write(6,903) diff_order, diff_coeff_uvprint
            write(6,901) diff_timescale
            If ( model_domain == mt_global ) Then
              write(6,905) polar_cap
              write(6,906) diff_coeff_ref
              write(6,907) diff_time_equ
            else !  model_domain  /=  mt_global
              write(6,*) ' This is the time taken to damp',             &
     &                 ' the 2 gridlength wave by 1/e '
            EndIf ! model_domain == mt_global
          Else ! L_diff_wind = .false.
            L_diff_message = .true.
            If ( L_diffusion .or. L_cdiffusion ) Then
              k = 1
              do ! k = 1, model_levels but end on EXIT
                if ( diffusion_coefficient_wind(k) > tiny .and.         &
     &                diffusion_order_wind(k) > 0 ) then
                  L_diff_message = .false.
                  Exit
                endif
                k = k + 1
                if( k > model_levels) Exit
              end do
            endIf ! L_diffusion .or. L_cdiffusion
            If ( L_diff_message ) Then
              write(6,*) ' *** Horizontal diffusion for winds is ',     &
     &                     ' NOT ACTIVE ***'
            endIf ! L_diff_message
          EndIf ! L_diff__wind

          if ( L_diff_incs ) then
            write(6,*) 'L_diff_incs = ',L_diff_incs,                    &
     &              ' so increments will be horizontally diffused.'
          endif ! L_diff_incs
          if ( L_diff_w ) then
            write(6,*)'L_diff_w = ',L_diff_w,                           &
     &             ' so w will be horizontally diffused like theta'
          endif ! L_diff_w

        endIf  !  L_diff_auto

 !  LAM domains can have diffusion so ...
      elseif (L_filter) then
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

        exp1 = exp(1.0)
        rpower = -1.0 / diff_timescale
        diff_coeff_temp = 0.25 * (1.0 - exp1**rpower)
        diff_coeff_phi = diff_coeff_temp

        Do j = j_start_u, j_stop_u
          Do i = 1-offx, row_length
            diff_coeff_u(i,j) = diff_coeff_temp
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1-offx, row_length
            diff_coeff_v(i,j) = diff_coeff_temp
          End Do
        End Do

        write(6,*) ' Check filter sweeping set for non-global domains'
        write(6,*)' u_begin(0) = ', u_begin(0), ' u_end(0) = ', u_end(0)
        write(6,*)' v_begin(0) = ', v_begin(0), ' v_end(0) = ', v_end(0)
        write(6,*) ' Only 1st order diffusion allowed - '               &
     &         , ' (overwriting UMUI value supplied). '
        write(6,*) ' 1st order diffusion e-folding time scale = '       &
     &         , diff_timescale, ' time steps'
        write(6,*) ' diffusion coefficient = ', diff_coeff_temp

      EndIf ! L_filter .and. model_domain  /=  mt_global

      if( L_sponge ) then
        Write ( Unit = 6, fmt=*) '  Sponge zone is active '
        Write ( Unit = 6, fmt=*) '  Sides sponge width  = ', sponge_ew
        Write ( Unit = 6, fmt=*) '  Northern/Southern sponge widths = ' &
     &                             , sponge_ns
        if( sponge_power == 1 ) then
          weight = 1.0 / (sponge_ew + 1)
          sponge_wts_ew(1) = 1.0 - weight
          do i = 2, sponge_ew
           sponge_wts_ew(i) = sponge_wts_ew(i-1) - weight
          enddo
          weight = 1.0 / (sponge_ns + 1)
           sponge_wts_ns(1) = 1.0 - weight
          do i = 2, sponge_ns
           sponge_wts_ns(i) = sponge_wts_ns(i-1) - weight
          enddo
          Write ( Unit = 6, fmt=*) ' sponge_power = ', sponge_power,    &
     &      '  therefore sponge blending is linear'
          Write ( Unit = 6, fmt=*) ' East/West sponge weights '
          Write ( Unit = 6, fmt=*) sponge_wts_ew
          Write ( Unit = 6, fmt=*) ' North/South sponge weights '
          Write ( Unit = 6, fmt=*) sponge_wts_ns
        elseif( sponge_power == 2 ) then
          weight = 1.0 / (sponge_ew + 1)
          sponge_wts_ew(1) = 1.0 - weight
          do i = 2, sponge_ew
           sponge_wts_ew(i) = sponge_wts_ew(i-1) - weight
           sponge_wts_ew(i-1) = sponge_wts_ew(i-1) * sponge_wts_ew(i-1)
          enddo
           sponge_wts_ew(sponge_ew) = sponge_wts_ew(sponge_ew) *        &
     &                                sponge_wts_ew(sponge_ew)
          weight = 1.0 / (sponge_ns + 1)
           sponge_wts_ns(1) = 1.0 - weight
          do i = 2, sponge_ns
           sponge_wts_ns(i) = sponge_wts_ns(i-1) - weight
           sponge_wts_ns(i-1) = sponge_wts_ns(i-1) * sponge_wts_ns(i-1)
          enddo
           sponge_wts_ns(sponge_ns) = sponge_wts_ns(sponge_ns) *        &
     &                                sponge_wts_ns(sponge_ns)
          Write ( Unit = 6, fmt=*) ' sponge_power = ', sponge_power,    &
     &      '  therefore sponge blending is quadratic'
          Write ( Unit = 6, fmt=*) ' East/West sponge weights '
          Write ( Unit = 6, fmt=*) sponge_wts_ew
          Write ( Unit = 6, fmt=*) ' North/South sponge weights '
          Write ( Unit = 6, fmt=*) sponge_wts_ns
        else
          Write ( Unit = 6, fmt=*) ' sponge_power = ', sponge_power,    &
     &      '  is NOT PERMITTED'
        endif ! sponge_power == 1
      endif ! L_sponge

!   quick fix to fix bit-reprod problem with filtering increments
! DEPENDS ON: swap_bounds
        call Swap_Bounds(                                               &
     &                   diff_coeff_u, row_length, rows, 1,             &
     &                   offx, offy, fld_type_u, .false.)
! DEPENDS ON: swap_bounds
        call Swap_Bounds(                                               &
     &                   diff_coeff_v, row_length, n_rows, 1,           &
     &                   offx, offy, fld_type_u, .false.)

  WRITE(6,*) ' ENDGame diffusion settings:'
  WRITE(6,*) '        L_diffusion', L_diffusion
  WRITE(6,*) '        L_cdiffusion', L_cdiffusion
  WRITE(6,*) '        L_filter', L_filter
  WRITE(6,*) '        L_diff_auto', L_diff_auto
  WRITE(6,*) '        L_pfcomb', L_pfcomb
  WRITE(6,*) '        L_pfexner', L_pfexner
  WRITE(6,*) '        L_pftheta', L_pftheta
  WRITE(6,*) '        L_pfuv', L_pfuv
  WRITE(6,*) '        L_pfw', L_pfw
  WRITE(6,*) '        L_pfincs',L_pfincs
  WRITE(6,*) '        L_diff_thermo', L_diff_thermo
  WRITE(6,*) '        L_diff_wind', L_diff_wind
  WRITE(6,*) '        L_diff_w', L_diff_w
  WRITE(6,*) '        L_diff_incs',L_diff_incs

! ----------------------------------------------------------------------
! Section 2. Print FORMATTING
! ----------------------------------------------------------------------

 901  FORMAT(' Diffusion timescale for wind is ',I4,' timesteps ')
 903  FORMAT(' Diffusion order is ',I2,                                 &
     &       ' with diffusion coefficient = ',F6.4)
 904  FORMAT(' Diffusion timescale for theta is ',I4,' timesteps ')
 905  FORMAT( ' This is the time taken to damp the 2 gridlength wave',  &
     &        ' at ',F6.2,' degrees by 1/e ')
 906  FORMAT( ' The diffusion coefficient at the equator = ',E8.2)
 907  FORMAT( ' This is a damping timescale of ',F7.1,' timesteps')
 910  FORMAT( '  Polar filter activates polewards of ',F6.2,' degrees')
 911  FORMAT( '  and needs a maximum of ',I3,                           &
     &        ' pass(es) of the 1-2-1 filter')
 912  FORMAT( ' polar_cap reset to ',F6.2,' degrees')
 917  FORMAT('Diffusion timescale for v-row ',I4                        &
     &                                   ,'  is ',F6.1,' timesteps')
 918  FORMAT('Diffusion timescale for u-row ',I4                        &
     &                                   ,'  is ',F6.1,' timesteps')
 919  FORMAT('Diffusion timescale at Equator is ',F6.1,' timesteps')
 920  FORMAT(' By setting polar_cap = ',F6.2,' <',F6.2,                 &
     &       ' (row next to pole), ')
 921  FORMAT(' polar filtering will stop at row before ',F6.2)

      IF (lhook) CALL dr_hook('EG_SETDIFF',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE eg_Setdiff
