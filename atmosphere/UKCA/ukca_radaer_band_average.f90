! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Average optical properties of UKCA-MODE aerosols, as obtained from
!+ look-up tables, over spectral wavebands.
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
SUBROUTINE ukca_radaer_band_average(                                    &
      ! Spectral information
      n_band, isolir, l_exclude, n_band_exclude, index_exclude          &
      ! Actual array dimensions
   ,  n_profile, n_layer, n_ukca_mode, n_ukca_cpnt                      &
      ! structure for UKCA/radiation interaction
   ,  ukca_radaer                                                       &
      ! Modal mass-mixing ratios from UKCA module
   ,  ukca_modal_mmr                                                    &
      ! Modal number concentrations from UKCA module
   ,  ukca_modal_number                                                 &
      ! Modal diameters from UKCA module
   ,  ukca_dry_diam, ukca_wet_diam                                      &
      ! Other inputs from UKCA module
   ,  ukca_cpnt_volume, ukca_modal_volume, ukca_modal_density           &
   ,  ukca_water_volume                                                 &
      ! Band-averaged optical properties (outputs)
   ,  ukca_absorption, ukca_scattering, ukca_asymmetry                  &
      ! Fixed array dimensions 
   ,  npd_profile, npd_layer, npd_band, npd_exclude                     &
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

USE spcrg3a_mod, ONLY: ip_solar, ip_infra_red
IMPLICIT NONE

!
! Arguments with intent(in)
!
! Current spectrum
!
INTEGER :: isolir
!
! Fixed array dimensions
!
INTEGER :: npd_profile, &
           npd_layer,   &
           npd_band,    &
           npd_exclude
!
! Actual array dimensions
!
INTEGER :: n_profile,   &
           n_layer,     &
           n_band,      &
           n_ukca_mode, &
           n_ukca_cpnt
!
! Variables related to waveband exclusion
!
LOGICAL l_exclude
INTEGER, DIMENSION(npd_band) :: n_band_exclude
INTEGER, DIMENSION(npd_exclude, npd_band) :: index_exclude
!
! Structure for UKCA/radiation interaction
!
TYPE (ukca_radaer_struct) :: ukca_radaer

!
! Modal mass-mixing ratios
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_modal_mmr

!
! Modal number concentrations (m-3)
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_modal_number

!
! Dry and wet modal diameters
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_dry_diam
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_wet_diam

!
! Component volumes
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_cpnt) :: &
                                                     ukca_cpnt_volume

!
! Modal volumes and densities
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: &
                                                    ukca_modal_volume
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: &
                                                    ukca_modal_density

!
! Volume of water in modes
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: &
                                                    ukca_water_volume

!
! Arguments with intent(out)
!
! Band-averaged modal optical properties
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode, npd_band) :: &
  ukca_absorption
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode, npd_band) :: &
  ukca_scattering
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode, npd_band) :: &
  ukca_asymmetry

!
!
! Local variables
!
!
!
! Spectrum definitions
!

!
! Values at the point of integration:
!      Mie parameter for the wet and dry diameters and the indices of 
!      their nearest neighbour
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
! Integrals
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode, npd_band) :: &
  integrated_abs
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode, npd_band) :: &
  integrated_sca
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode, npd_band) :: &
  integrated_asy
REAL :: loc_abs(precalc%n_integ_pts)
REAL :: loc_sca(precalc%n_integ_pts)
REAL :: loc_asy(precalc%n_integ_pts)
REAL :: loc_vol
REAL :: factor

!
! Waveband-integrated flux corrected for exclusions
!
REAL exclflux

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
! Loop variables
!
INTEGER i_band, & ! loop on wavebands
        i_mode, & ! loop on aerosol modes
        i_cmpt, & ! loop on aerosol components
        i_layr, & ! loop on vertical dimension
        i_prof, & ! loop on horizontal dimension
        i_intg    ! loop on integration points and excluded bands

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

IF (lhook) CALL dr_hook('UKCA_RADAER_BAND_AVERAGE', zhook_in, zhook_handle)

!
! To band-average modal optical properties, we need first to compute
! adequate indices in the look-up tables. For that, we need:
! *** the modal dry radius (we've got the diameter as input)
! *** the modal wet radius (we've got the diameter as input)
! *** the modal refractive index (computed as volume-weighted component
!     refractive indices)
! In addition, in order to output specific coefficients for absorption
! and scattering (in m2/kg from m-1), we need the modal density.
!

DO i_band = 1, n_band

  DO i_mode = 1, n_ukca_mode
    
    !
    ! Mode type. From a look-up table point of view, Aitken and
    ! accumulation types are treated in the same way.
    ! Once we know which look-up table to select, make local copies
    ! of info needed for nearest-neighbour calculations.
    !
    IF (ukca_radaer%i_mode_type(i_mode) == IP_UKCA_MODE_AITKEN .OR. &
        ukca_radaer%i_mode_type(i_mode) == IP_UKCA_MODE_ACCUM) THEN
      this_mode_type = IP_UKCA_LUT_ACCUM
    ELSE ! IP_UKCA_MODE_COARSE
      this_mode_type = IP_UKCA_LUT_COARSE
    END IF

    nx      = ukca_lut(this_mode_type, isolir)%n_x
    logxmin = log(ukca_lut(this_mode_type, isolir)%x_min)
    logxmaxmlogxmin = &
              log(ukca_lut(this_mode_type, isolir)%x_max) - logxmin
    
    nnr     = ukca_lut(this_mode_type, isolir)%n_nr
    nrmin   = ukca_lut(this_mode_type, isolir)%nr_min
    incr_nr = ukca_lut(this_mode_type, isolir)%incr_nr
    
    nni     = ukca_lut(this_mode_type, isolir)%n_ni
    nimin   = ukca_lut(this_mode_type, isolir)%ni_min
    incr_ni = ukca_lut(this_mode_type, isolir)%incr_ni
  
    !
    ! Wavelength-dependent calculations.
    ! Waveband-integration is done at the same time, so most computed
    ! items are not stored in arrays.
    !
    DO i_layr = 1, n_layer
  
      DO i_prof = 1, n_profile
        
        !
        ! Only make calculations if there are some aerosols, and
        ! if the number concentration is large enough.
        ! This test is especially important for the first timestep,
        ! as UKCA has not run yet and its output is therefore
        ! not guaranteed to be valid. Mass mixing ratios and numbers 
        ! are initialised to zero as prognostics.
        ! Also, at low number concentrations, the size informations
        ! given by UKCA are unreliable and might produce erroneous
        ! optical properties.
        !

        IF (ukca_modal_mmr(i_prof, i_layr, i_mode) > threshold_mmr .AND. &
            ukca_modal_number(i_prof, i_layr, i_mode) > threshold_nbr) THEN
      
          DO i_intg = 1, precalc%n_integ_pts
        
            !
            ! Compute the Mie parameter from the wet diameter
            ! and get the LUT-array index of its nearest neighbour.
            !
            x = Pi * ukca_wet_diam(i_prof, i_layr, i_mode) / &
                precalc%wavelength(i_intg, i_band, isolir)        
            n_x = NINT( (log(x)    - logxmin) / &
                         logxmaxmlogxmin * (nx-1) ) + 1
            n_x = MIN(nx, MAX(1, n_x))
                        
            !
            ! Same for the dry diameter (needed to access the volume
            ! fraction)
            !
            x_dry = Pi * ukca_dry_diam(i_prof, i_layr, i_mode) / &
                    precalc%wavelength(i_intg, i_band, isolir)
            n_x_dry = NINT( (log(x_dry) - logxmin) / &
                             logxmaxmlogxmin * (nx-1) ) + 1
            n_x_dry = MIN(nx, MAX(1, n_x_dry))
               
            !
            ! Compute the modal complex refractive index as 
            ! volume-weighted component refractive indices.
            ! Get the LUT-array index of their nearest neighbours.
            !
            re_m = 0.0E+00
            im_m = 0.0E+00
        
            DO i_cmpt = 1, ukca_radaer%n_cpnt_in_mode(i_mode)
        
              this_cpnt = ukca_radaer%i_cpnt_index(i_cmpt, i_mode)
              this_cpnt_type = ukca_radaer%i_cpnt_type(this_cpnt)
          
              re_m = re_m + ukca_cpnt_volume(i_prof, i_layr, this_cpnt) * &
                            precalc%realrefr(this_cpnt_type, &
                                             i_intg, i_band, isolir)
              im_m = im_m + ukca_cpnt_volume(i_prof, i_layr, this_cpnt) * &
                            precalc%imagrefr(this_cpnt_type, &
                                             i_intg, i_band, isolir)
        
            END DO ! i_cmpt
          
            IF (ukca_radaer%l_soluble(i_mode)) THEN
            
              !
              ! Account for refractive index of water
              !
              re_m = re_m + ukca_water_volume(i_prof, i_layr, i_mode) * &
                            precalc%realrefr(IP_UKCA_WATER, &
                                             i_intg, i_band, isolir)
              im_m = im_m + ukca_water_volume(i_prof, i_layr, i_mode) * &
                            precalc%imagrefr(IP_UKCA_WATER, &
                                             i_intg, i_band, isolir)
            
            END IF ! l_soluble
            
            re_m = re_m / ukca_modal_volume(i_prof, i_layr, i_mode)
            im_m = im_m / ukca_modal_volume(i_prof, i_layr, i_mode)
            
            n_nr = NINT( (re_m - nrmin) / incr_nr ) + 1
            n_nr = MIN(nnr, MAX(1, n_nr))
          
            n_ni = NINT( (im_m - nimin) / incr_ni ) + 1
            n_ni = MIN(nni, MAX(1, n_ni))

            !
            ! Get local copies of the relevant look-up table entries.
            !
            loc_abs(i_intg) = ukca_lut(this_mode_type, isolir)% &
                         ukca_absorption(n_x, n_ni, n_nr)

            loc_sca(i_intg) = ukca_lut(this_mode_type, isolir)% &
                         ukca_scattering(n_x, n_ni, n_nr)
                  
            loc_asy(i_intg) = ukca_lut(this_mode_type, isolir)% &
                         ukca_asymmetry(n_x, n_ni, n_nr)
                  
            loc_vol = ukca_lut(this_mode_type, isolir)% &
                      volume_fraction(n_x_dry)

            !
            ! Offline Mie calculations were integrated using the Mie
            ! parameter. Compared to an integration using the particle
            ! radius, extra factors are introduced. Absorption and
            ! scattering efficiencies must be multiplied by the squared
            ! wavelength, and the volume fraction by the cubed wavelength.
            ! Consequently, ratios abs/volfrac and sca/volfrac have then
            ! to be divided by the wavelength.
            ! We also weight by the solar irradiance or Planckian 
            ! irradiance.
            !
            factor = precalc%irrad(i_intg, i_band, isolir) / &
               (ukca_modal_density(i_prof, i_layr, i_mode) * loc_vol * &
                precalc%wavelength(i_intg, i_band, isolir))
            loc_abs(i_intg) = loc_abs(i_intg) * factor
            loc_sca(i_intg) = loc_sca(i_intg) * factor
            loc_asy(i_intg) = loc_asy(i_intg) * loc_sca(i_intg)
            
          END DO ! i_intg
          
          !
          ! Trapezoidal integration 
          !
          integrated_abs(i_prof, i_layr, i_mode, i_band) = 0.0E+00
          integrated_sca(i_prof, i_layr, i_mode, i_band) = 0.0E+00
          integrated_asy(i_prof, i_layr, i_mode, i_band) = 0.0E+00
          DO i_intg = 1, precalc%n_integ_pts - 1
            integrated_abs(i_prof, i_layr, i_mode, i_band) = &
               integrated_abs(i_prof, i_layr, i_mode, i_band) + &
               (precalc%wavelength(i_intg+1, i_band, isolir) - &
                precalc%wavelength(i_intg, i_band, isolir)) * &
               (loc_abs(i_intg+1) + loc_abs(i_intg))
            integrated_sca(i_prof, i_layr, i_mode, i_band) = &
               integrated_sca(i_prof, i_layr, i_mode, i_band) + &
               (precalc%wavelength(i_intg+1, i_band, isolir) - &
                precalc%wavelength(i_intg, i_band, isolir)) * &
               (loc_sca(i_intg+1) + loc_sca(i_intg))
            integrated_asy(i_prof, i_layr, i_mode, i_band) = &
               integrated_asy(i_prof, i_layr, i_mode, i_band) + &
               (precalc%wavelength(i_intg+1, i_band, isolir) - &
                precalc%wavelength(i_intg, i_band, isolir)) * &
               (loc_asy(i_intg+1) + loc_asy(i_intg))
          END DO ! i_intg
          integrated_abs(i_prof, i_layr, i_mode, i_band) = &
             integrated_abs(i_prof, i_layr, i_mode, i_band) * 0.5
          integrated_sca(i_prof, i_layr, i_mode, i_band) = &
             integrated_sca(i_prof, i_layr, i_mode, i_band) * 0.5
          integrated_asy(i_prof, i_layr, i_mode, i_band) = &
             integrated_asy(i_prof, i_layr, i_mode, i_band) * 0.5
          
        ELSE
        
          integrated_abs(i_prof, i_layr, i_mode, i_band) = 0.0E+00
          integrated_sca(i_prof, i_layr, i_mode, i_band) = 0.0E+00
          integrated_asy(i_prof, i_layr, i_mode, i_band) = 0.0E+00
        
        END IF

      END DO ! i_prof
    
    END DO ! i_layr

  END DO ! i_mode

END DO ! i_band

!
! Final integrals. Depend on excluded bands.
!
DO i_band = 1, n_band
  
  IF (l_exclude) THEN
    
    IF (n_band_exclude(i_band) > 0) THEN
      
      !
      ! Remove contribution from excluded bands.
      !
      DO i_intg = 1, n_band_exclude(i_band)
      
        DO i_mode = 1, n_ukca_mode
          
          DO i_layr = 1, n_layer
          
            DO i_prof = 1, n_profile
          
              integrated_abs(i_prof, i_layr, i_mode, i_band) = &
                        integrated_abs(i_prof, i_layr, i_mode, i_band) - &
                        integrated_abs(i_prof, i_layr, i_mode, &
                                       index_exclude(i_intg, i_band))
              integrated_sca(i_prof, i_layr, i_mode, i_band) = &
                        integrated_sca(i_prof, i_layr, i_mode, i_band) - &
                        integrated_sca(i_prof, i_layr, i_mode, &
                                       index_exclude(i_intg, i_band))
              integrated_asy(i_prof, i_layr, i_mode, i_band) = &
                        integrated_asy(i_prof, i_layr, i_mode, i_band) - &
                        integrated_asy(i_prof, i_layr, i_mode, &
                                       index_exclude(i_intg, i_band))
      
            END DO ! i_prof
          
          END DO ! i_layr
        
        END DO ! i_mode
      
        exclflux = precalc%flux(i_band, isolir) - &
                   precalc%flux(index_exclude(i_intg, i_band), isolir)
      
      END DO ! i_intg
      
    ELSE
    
      exclflux = precalc%flux(i_band, isolir)
    
    END IF
    
  ELSE
  
    exclflux = precalc%flux(i_band, isolir)
  
  END IF
  
  DO i_mode = 1, n_ukca_mode

    DO i_layr = 1, n_layer
          
      DO i_prof = 1, n_profile
    
        ukca_absorption(i_prof, i_layr, i_mode, i_band) = &
          integrated_abs(i_prof, i_layr, i_mode, i_band) / exclflux
        ukca_scattering(i_prof, i_layr, i_mode, i_band) = &
          integrated_sca(i_prof, i_layr, i_mode, i_band) / exclflux
        IF (integrated_sca(i_prof, i_layr, i_mode, i_band) > 0.0E+00) THEN
          ukca_asymmetry(i_prof, i_layr, i_mode, i_band)  = &
            integrated_asy(i_prof, i_layr, i_mode, i_band) / &
            integrated_sca(i_prof, i_layr, i_mode, i_band)
        ELSE
          ukca_asymmetry(i_prof, i_layr, i_mode, i_band) = 0.0E+00
        END IF
      
      END DO ! i_prof
    
    END DO ! i_layr
    
  END DO  ! i_mode

END DO ! i_band

IF (lhook) CALL dr_hook('UKCA_RADAER_BAND_AVERAGE', zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_band_average
