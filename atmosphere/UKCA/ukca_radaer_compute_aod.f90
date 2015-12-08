! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Obtain UKCA-MODE aerosol optical properties at optical depth
!+ wavelengths using UKCA_RADAER look-up tables, and compute the
!+ UKCA-MODE modal optical depth for the mode requested.
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
SUBROUTINE ukca_radaer_compute_aod(                                     &
      ! UKCA_RADAER model switch
      l_ukca_radaer,                                                    &
      ! Actual array dimensions
      n_profile, n_layer, n_ukca_mode, n_ukca_cpnt, n_aod_wavel,        &
      ! UKCA_RADAER structure
      ukca_radaer,                                                      &
      ! Modal diameters from UKCA module
      ukca_dry, ukca_wet,                                               &
      ! Mass thickness of layers
      d_mass,                                                           &
      ! Component volumes
      ukca_cpnt_volume,                                                 &
      ! Modal volumes, densities, and water content
      ukca_modal_volume, ukca_modal_density, ukca_water_volume,         &
      ! Modal mass-mixing ratios
      ukca_modal_mixr,                                                  &
      ! Modal number concentrations
      ukca_modal_number,                                                &
      ! Type selection
      type_wanted, soluble_wanted,                                      &
      ! Modal extinction aerosol optical depths
      ukca_modal_aod,                                                   &
      ! Fixed array dimensions
      npd_profile,  npd_layer, npd_aod_wavel                            &
   )

USE parkind1, ONLY: jpim, jprb 
USE yomhook,  ONLY: lhook, dr_hook
USE conversions_mod, ONLY: pi

! UKCA look-up tables
!
USE ukca_radaer_lut

!
! UKCA pre-computed values
USE ukca_radaer_precalc

USE ukca_radaer_struct_mod

IMPLICIT NONE

!
! Arguments with intent in
!

! UKCA_RADAER model switch
LOGICAL l_ukca_radaer

! Fixed array dimensions
INTEGER npd_profile,  &
        npd_layer,    &
        npd_aod_wavel

! Actual array dimensions
INTEGER n_profile,   &
        n_layer,     &
        n_ukca_mode, &
        n_ukca_cpnt, &
        n_aod_wavel

! Structure for UKCA/radiation interaction
TYPE (ukca_radaer_struct) :: ukca_radaer

! Modal dry and wet diameters
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_dry
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_wet

! Mass-thickness of vertical layers
REAL, DIMENSION(npd_profile, npd_layer) :: d_mass

! Component volumes
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_cpnt) :: &
                                                        ukca_cpnt_volume

! Modal volumes
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: &
                                                       ukca_modal_volume

! Modal density
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: &
                                                      ukca_modal_density

! Volume of water in each mode
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: &
                                                       ukca_water_volume


! Modal mass-mixing ratios
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_modal_mixr

! Modal number concentrations
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_modal_number

! Type selection
INTEGER type_wanted
LOGICAL soluble_wanted

!
! Arguments with intent out
!

! Modal aerosol optical depths
REAL, DIMENSION(npd_profile, npd_aod_wavel) :: ukca_modal_aod
  
!
! Local variables
!
INTEGER i, &
        j, &
        k, &
        l, &
        n

!
! Values at the AOD wavelengths:
!      Mie parameter for the wet and dry diameters and 
!      the indices of their nearest neighbour
!      Complex refractive index and the index of its nearest neighbour
!
REAL x
INTEGER n_x
REAL x_dry
INTEGER n_x_dry
REAL re_m 
INTEGER n_nr
REAL im_m
INTEGER n_ni

!
! Values at given wavelength
!
REAL :: loc_abs
REAL :: loc_sca
REAL :: loc_vol
REAL :: factor

!
! Local copies of typedef members
!
INTEGER nx
REAL logxmin         ! log(xmin)
REAL logxmaxmlogxmin ! log(xmax) - log(xmin)
INTEGER nnr
REAL nrmin
REAL incr_nr
INTEGER nni
REAL nimin
REAL incr_ni

!
! Local copies of mode type, component index and component type
!
INTEGER this_mode_type
INTEGER this_cpnt
INTEGER this_cpnt_type

!
! Thresholds on the modal mass-mixing ratio and modal number
! concentrations above which aerosol optical properties are to be
! computed.
!
REAL, PARAMETER :: threshold_mmr = 0.0E+00 ! kg/kg
REAL, PARAMETER :: threshold_nbr = 1.0E+04 ! m-3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!
!
!

IF (lhook) CALL dr_hook('UKCA_RADAER_COMPUTE_AOD', zhook_in, zhook_handle)

ukca_modal_aod(:, :)     = 0.0E+00

IF (l_ukca_radaer) THEN

  DO i = 1, n_ukca_mode

    IF ((ukca_radaer%i_mode_type(i) == type_wanted) .AND. &
        (ukca_radaer%l_soluble(i).EQV.soluble_wanted)) THEN
    
      !
      ! Get the right look-up table. It depends on the type of mode
      ! (Aitken/accumulation or coarse). For the aerosol optical depth,
      ! the shortwave spectrum LUTs are used.
      !
    
      IF (ukca_radaer%i_mode_type(i) == ip_ukca_mode_aitken .OR. &
          ukca_radaer%i_mode_type(i) == ip_ukca_mode_accum) THEN
        this_mode_type = ip_ukca_lut_accum
      ELSE ! IP_UKCA_MODE_COARSE
        this_mode_type = ip_ukca_lut_coarse
      END IF
    
      nx      = ukca_lut(this_mode_type, ip_ukca_lut_sw)%n_x
      logxmin = log(ukca_lut(this_mode_type, ip_ukca_lut_sw)%x_min)
      logxmaxmlogxmin = &
              log(ukca_lut(this_mode_type, ip_ukca_lut_sw)%x_max) - logxmin
    
      nnr     = ukca_lut(this_mode_type, ip_ukca_lut_sw)%n_nr
      nrmin   = ukca_lut(this_mode_type, ip_ukca_lut_sw)%nr_min
      incr_nr = ukca_lut(this_mode_type, ip_ukca_lut_sw)%incr_nr
    
      nni     = ukca_lut(this_mode_type, ip_ukca_lut_sw)%n_ni
      nimin   = ukca_lut(this_mode_type, ip_ukca_lut_sw)%ni_min
      incr_ni = ukca_lut(this_mode_type, ip_ukca_lut_sw)%incr_ni
    
      DO k = 1, n_layer
    
        DO l = 1, n_profile
        
        IF (ukca_modal_mixr(l, k, i) > threshold_mmr .AND. &
            ukca_modal_number(l, k, i) > threshold_nbr) THEN
            
            DO n = 1, n_aod_wavel
              
              x = Pi * ukca_wet(l, k, i) / precalc%aod_wavel(n)        
              n_x = NINT( (log(x)    - logxmin) / &
                           logxmaxmlogxmin * (nx-1) ) + 1
              n_x = MIN(nx, MAX(1, n_x))
              
              x_dry = Pi * ukca_dry(l, k, i) / precalc%aod_wavel(n)
              n_x_dry = NINT( (log(x_dry) - logxmin) / &
                               logxmaxmlogxmin * (nx-1) ) + 1
              n_x_dry = MIN(nx, MAX(1, n_x_dry))
              
              re_m = 0.0E+00
              im_m = 0.0E+00
              
              DO j = 1, ukca_radaer%n_cpnt_in_mode(i)
                
                this_cpnt = ukca_radaer%i_cpnt_index(j, i)
                this_cpnt_type = ukca_radaer%i_cpnt_type(this_cpnt)
          
                re_m = re_m + ukca_cpnt_volume(l, k, this_cpnt) * &
                              precalc%aod_realrefr(this_cpnt_type, n)
                im_m = im_m +  ukca_cpnt_volume(l, k, this_cpnt) * &
                              precalc%aod_imagrefr(this_cpnt_type, n)
                
              END DO ! j
              
              IF (ukca_radaer%l_soluble(i)) THEN
              
                re_m = re_m + ukca_water_volume(l, k, i) * &
                              precalc%aod_realrefr(ip_ukca_water, n)
                im_m = im_m + ukca_water_volume(l, k, i) * &
                              precalc%aod_imagrefr(ip_ukca_water, n)
              
              END IF ! l_soluble
              
              re_m = re_m / ukca_modal_volume(l, k, i)
              im_m = im_m / ukca_modal_volume(l, k, i)
        
              n_nr = NINT( (re_m - nrmin) / incr_nr ) + 1
              n_nr = MIN(nnr, MAX(1, n_nr))
          
              n_ni = NINT( (im_m - nimin) / incr_ni ) + 1
              n_ni = MIN(nni, MAX(1, n_ni))
              
              loc_abs = ukca_lut(this_mode_type, ip_ukca_lut_sw)% &
                        ukca_absorption(n_x, n_ni, n_nr)

              loc_sca = ukca_lut(this_mode_type, ip_ukca_lut_sw)% &
                        ukca_scattering(n_x, n_ni, n_nr)
                        
              loc_vol = ukca_lut(this_mode_type, ip_ukca_lut_sw)% &
                        volume_fraction(n_x_dry)
                        
              factor = (ukca_modal_density(l, k, i) * loc_vol * &
                        precalc%aod_wavel(n))
              
              loc_abs = loc_abs / factor
              loc_sca = loc_sca / factor
              
              !
              ! aerosol optical depth for extinction
              !
              ukca_modal_aod(l, n) = &
              ukca_modal_aod(l, n) + &
              ukca_modal_mixr(l, k, i) * d_mass(l, k) * (loc_abs + loc_sca)
              
            END DO ! n
        
          END IF
        
        END DO ! l
      
      END DO ! k
    
    END IF

  END DO ! i

END IF ! l_ukca_radaer

IF (lhook) CALL dr_hook('UKCA_RADAER_COMPUTE_AOD', zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_compute_aod
