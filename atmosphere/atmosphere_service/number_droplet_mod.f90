! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Function to calculate cloud droplet number concentration.

! Purpose:
!   Cloud droplet number concentration is calculated from aerosol
!   concentration or else fixed values are assigned.

! Method:
!   Sulphate aerosol mass concentration is converted to number
!   concentration by assuming a log-normal size distribution.
!   Sea-salt and/or biomass-burning aerosols may then
!   be added if required. The total is then converted to cloud
!   droplet concentration following the parametrization of
!   Jones et al. (1994) and lower limits are imposed.
!   Alternatively, fixed droplet values are assigned if the
!   parametrization is not required.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Service

! Description of Code:
!   FORTRAN 90

!- ---------------------------------------------------------------------
MODULE number_droplet_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: number_droplet

! Variables for this module only
LOGICAL :: nd_init = .FALSE.         ! Are particle volumes initialised?
REAL    :: particle_volume_nh42so4   ! Mean volume of (NH4)2SO4
REAL    :: particle_volume_biomass   ! Mean volume of biomass smoke
REAL    :: particle_volume_biogenic  ! Mean volume of biogenic aerosol
REAL    :: particle_volume_ocff      ! Mean volume of OCFF aerosol
REAL    :: particle_volume_nitrate   ! Mean volume of nitrate aerosol


!             Median radius of log-normal distribution for (NH4)2SO4
REAL, PARAMETER :: radius_0_nh42so4  = 9.5e-8
!             Geometric standard deviation of same
REAL, PARAMETER :: sigma_0_nh42so4   = 1.4
!             Density of ammonium sulphate aerosol
REAL, PARAMETER :: density_nh42so4   = 1.769e+03
!             Median radius of log-normal distribution for biomass smoke
REAL, PARAMETER :: radius_0_biomass  = 1.2e-07
!             Geometric standard deviation of same
REAL, PARAMETER :: sigma_0_biomass   = 1.30
!             Density of biomass smoke aerosol
REAL, PARAMETER :: density_biomass   = 1.35e+03
!             Median radius of log-normal dist. for biogenic aerosol
REAL, PARAMETER :: radius_0_biogenic = 9.5e-08
!             Geometric standard deviation of same
REAL, PARAMETER :: sigma_0_biogenic  = 1.50
!             Density of biogenic aerosol
REAL, PARAMETER :: density_biogenic  = 1.3e+03
!             Median radius of log-normal dist. for OCFF aerosol
REAL, PARAMETER :: radius_0_ocff     = 0.12e-06
!             Geometric standard deviation of same
REAL, PARAMETER :: sigma_0_ocff      = 1.30
!             Density of OCFF aerosol
REAL, PARAMETER :: density_ocff      = 1350.0
!             Median radius of log-normal dist. for nitrate aerosol
REAL, PARAMETER :: radius_0_nitrate  = 9.5e-8
!             Geometric standard deviation of same
REAL, PARAMETER :: sigma_0_nitrate   = 1.4
!             Density of nitrate aerosol
REAL, PARAMETER :: density_nitrate   = 1.725e+03


CONTAINS 
!- ---------------------------------------------------------------------
  SUBROUTINE number_droplet_init()
  USE conversions_mod, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('NUMBER_DROPLET_INIT',zhook_in,zhook_handle)

    ! Initialize particle volume constants
  particle_volume_nh42so4 =                                           &
      (4.0e+00*pi/3.0e+00)*radius_0_nh42so4**3                        &
      *EXP(4.5e+00*(LOG(sigma_0_nh42so4))**2)

  particle_volume_biomass =                                           &
      (4.0e+00*pi/3.0e+00)*radius_0_biomass**3                        &
      *EXP(4.5e+00*(LOG(sigma_0_biomass))**2)

  particle_volume_biogenic =                                          &
      (4.0e+00*pi/3.0e+00)*radius_0_biogenic**3                       &
      *EXP(4.5e+00*(LOG(sigma_0_biogenic))**2)

  particle_volume_ocff =                                              &
      (4.0e+00*pi/3.0e+00)*radius_0_ocff**3                           &
      *EXP(4.5e+00*(LOG(sigma_0_ocff))**2)

  particle_volume_nitrate =                                           &
      (4.0e+00*pi/3.0e+00)*radius_0_nitrate**3                        &
      *EXP(4.5e+00*(LOG(sigma_0_nitrate))**2)

  ! Now volumes are initialised
  nd_init = .TRUE.

  IF (lhook) CALL dr_hook('NUMBER_DROPLET_INIT',zhook_out,zhook_handle)

  RETURN
  END SUBROUTINE number_droplet_init

!- ---------------------------------------------------------------------

  SUBROUTINE number_droplet(                                        &
       i_start, i_end, j_start, j_end,                              &
       k_start, k_end, k_loop_start, k_loop_end,                    &
       l_aerosol_droplet, l_nh42so4,                                &
       accum_sulphate, diss_sulphate,                               &
       l_seasalt_ccn, sea_salt_film, sea_salt_jet,                  &
       l_biogenic_ccn, biogenic,                                    &
       l_biomass_ccn, biomass_aged, biomass_cloud,                  &
       l_ocff_ccn, ocff_aged, ocff_cloud,                           &
       l_nitrate_ccn, nitr_acc, nitr_diss,                          &
       density_air,                                                 &
       snow_depth,                                                  &
       land_fract,                                                  &
       ntot_land, ntot_sea,                                         &
       n_drop                                                       &
     )

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i_start      ! i indexing start
  INTEGER, INTENT(IN) :: i_end        ! i indexing end
  INTEGER, INTENT(IN) :: j_start      ! j indexing start
  INTEGER, INTENT(IN) :: j_end        ! j indexing end
  INTEGER, INTENT(IN) :: k_start      ! k indexing start
  INTEGER, INTENT(IN) :: k_end        ! k indexing end
  INTEGER, INTENT(IN) :: k_loop_start ! k looping start
  INTEGER, INTENT(IN) :: k_loop_end   ! k looping end

  !             Flag to use aerosols to find droplet number
  LOGICAL, INTENT(IN) :: l_aerosol_droplet
  !             Is the input "sulphate" aerosol in the form of
  !             ammonium sulphate (T) or just sulphur (F)? Used in the
  !             same way for the nitrate aerosol, in form of ammonium
  !             nitrate (T) or just nitrogen (F).
  LOGICAL, INTENT(IN) :: l_nh42so4
  !             Is sea-salt aerosol to be used?
  LOGICAL, INTENT(IN) :: l_seasalt_ccn
  !             Is biomass smoke aerosol to be used?
  LOGICAL, INTENT(IN) :: l_biomass_ccn
  !             Is biogenic aerosol to be used?
  LOGICAL, INTENT(IN) :: l_biogenic_ccn 
  !             Is fossil-fuel organic carbon aerosol to be used?
  LOGICAL, INTENT(IN) :: l_ocff_ccn
  !             Is ammonium nitrate aerosol to be used?
  LOGICAL, INTENT(IN) :: l_nitrate_ccn

  !             Mixing ratio of accumulation-mode sulphate aerosol
  REAL, INTENT(IN) :: accum_sulphate ( i_start : i_end,     &
                                       j_start : j_end,     &
                                       k_start : k_end )

!             Mixing ratio of dissolved sulphate aerosol
  REAL, INTENT(IN) :: diss_sulphate ( i_start : i_end,      &
                                      j_start : j_end,      &
                                      k_start : k_end )

  !             Mixing ratio of aged biomass smoke
  REAL, INTENT(IN) :: biomass_aged ( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Mixing ratio of in-cloud biomass smoke
  REAL, INTENT(IN) :: biomass_cloud( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Number concentration of film-mode sea salt aerosol (m-3)
  REAL, INTENT(IN) :: sea_salt_film( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Number concentration of jet-mode sea salt aerosol (m-3)
  REAL, INTENT(IN) :: sea_salt_jet ( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Mixing ratio of biogenic aerosol
  REAL, INTENT(IN) :: biogenic     ( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Mixing ratio of aged fossil-fuel organic carbon
  REAL, INTENT(IN) :: ocff_aged    ( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Mixing ratio of in-cloud fossil-fuel organic carbon
  REAL, INTENT(IN) :: ocff_cloud   ( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Mixing ratio of accumulation-mode nitrate aerosol
  REAL, INTENT(IN) :: nitr_acc     ( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Mixing ratio of dissolved-mode nitrate aerosol
  REAL, INTENT(IN) :: nitr_diss    ( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Density of air (kg m-3)
  REAL, INTENT(IN) :: density_air  ( i_start : i_end,      &
                                     j_start : j_end,      &
                                     k_start : k_end ) 

  !             Snow depth (m; >5000 is flag for ice-sheets)
  REAL, INTENT(IN) :: snow_depth   ( i_start : i_end,      &
                                     j_start : j_end)

  !             Land fraction
  REAL, INTENT(IN) :: land_fract   ( i_start : i_end,      &
                                     j_start : j_end)
 
  !             Droplet number over land if parameterization is off (m-3)
  REAL, INTENT(IN) :: ntot_land
  !             Droplet number over sea if parameterization is off (m-3)
  REAL, INTENT(IN) :: ntot_sea



  !             Returned number concentration of cloud droplets (m-3)
  REAL, INTENT(OUT) :: n_drop( i_start : i_end,      &
                               j_start : j_end,      &
                               k_start : k_end ) 

!     Local variables:

  INTEGER :: i       ! Looper
  INTEGER :: j       ! Looper
  INTEGER :: k       ! Looper

  !  Number density of CCN
  REAL :: n_ccn ( i_start : i_end,      &
                  j_start : j_end,      &
                  k_start : k_end )


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!---------------------------------------------------------------------------

  IF (lhook) CALL dr_hook('NUMBER_DROPLET',zhook_in,zhook_handle)

! Note the following call needs to be run by a single thread at a time
! to avoid race/deadlock conditions as this routine is sometimes called in
! a parallel context. 
!$OMP CRITICAL
  IF (.NOT. nd_init) CALL number_droplet_init()
!$OMP END CRITICAL

  IF (l_aerosol_droplet) THEN
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC)                      &
!$OMP&          PRIVATE(i, j, k)
    DO k = k_loop_start, k_loop_end
      DO j = j_start, j_end
        DO i = i_start, i_end

!        If active, aerosol concentrations are used to calculate the
!        number of CCN, which is then used to determine the number
!        concentration of cloud droplets (m-3).

        IF (l_nh42so4) THEN
!           Input data have already been converted to ammonium sulphate.
          n_ccn(i,j,k)=(accum_sulphate(i,j,k)+diss_sulphate(i,j,k))            &
             *density_air(i,j,k)/(density_nh42so4*particle_volume_nh42so4)
        ELSE
!           Convert m.m.r. of sulphur to ammonium sulphate by
!           multiplying by ratio of molecular weights:
          n_ccn(i,j,k)=(accum_sulphate(i,j,k)+diss_sulphate(i,j,k))*4.125      &
             *density_air(i,j,k)/(density_nh42so4*particle_volume_nh42so4)
        END IF

        IF (l_seasalt_ccn) THEN
          n_ccn(i,j,k)=n_ccn(i,j,k)+sea_salt_film(i,j,k)+sea_salt_jet(i,j,k)
        END IF

        IF (l_biomass_ccn) THEN
          n_ccn(i,j,k)=n_ccn(i,j,k)+((biomass_aged(i,j,k)+biomass_cloud(i,j,k))&
                       *density_air(i,j,k)                                     &
                       /(density_biomass*particle_volume_biomass))
        END IF

        IF (l_biogenic_ccn) THEN
          n_ccn(i,j,k)=n_ccn(i,j,k)+(biogenic(i,j,k)*density_air(i,j,k)        &
                       /(density_biogenic*particle_volume_biogenic))
        END IF

        IF (l_ocff_ccn) THEN
          n_ccn(i,j,k)=n_ccn(i,j,k)+((ocff_aged(i,j,k)+ocff_cloud(i,j,k))      &
                      *density_air(i,j,k)/(density_ocff*particle_volume_ocff))
        END IF

        IF (l_nitrate_ccn) THEN
          IF (l_nh42so4) THEN
!             Input data have already been converted to ammonium nitrate
            n_ccn(i,j,k)=n_ccn(i,j,k) + ((nitr_acc(i,j,k) + nitr_diss(i,j,k))  &
              * density_air(i,j,k)/(density_nitrate*particle_volume_nitrate))
          ELSE
!             Convert m.m.r. of nitrogen to ammonium nitrate by
!             multiplying by ratio of molecular weights
            n_ccn(i,j,k)=n_ccn(i,j,k) + ((nitr_acc(i,j,k) + nitr_diss(i,j,k))  &
                         * 5.714 * density_air(i,j,k)                          &
                         / (density_nitrate*particle_volume_nitrate))
          END IF
        END IF

!        Apply relation of Jones et al. (1994) to get droplet number
!        and apply minimum value (equivalent to 5 cm-3):

        n_drop(i,j,k)=3.75e+08*(1.0e+00-EXP(-2.5e-9*n_ccn(i,j,k)))

        IF (n_drop(i,j,k)  <   5.0e+06) THEN
          n_drop(i,j,k)=5.0e+06
        END IF

!        If gridbox is more than 20% land AND this land is not covered
!        by an ice-sheet, use larger minimum droplet number (=35 cm-3):

        IF (land_fract(i,j)  >   0.2 .AND. snow_depth(i,j)  <   5000.0         &
                            .AND. n_drop(i,j,k)  <   35.0e+06) THEN
          n_drop(i,j,k)=35.0e+06
        END IF

      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END PARALLEL DO

  ELSE  ! l_aerosol_droplet

!        Without aerosols, the number of droplets is fixed; a simple
!        50% criterion is used for land or sea in this case.
    DO k = k_loop_start, k_loop_end
      DO j = j_start, j_end
        DO i = i_start, i_end

         IF (land_fract(i,j)  >=  0.5) THEN
            n_drop(i,j,k)=ntot_land
         ELSE
            n_drop(i,j,k)=ntot_sea
         END IF
        END DO
      END DO
    END DO

  END IF ! l_aerosol_droplet



  IF (lhook) CALL dr_hook('NUMBER_DROPLET',zhook_out,zhook_handle)

  END SUBROUTINE number_droplet

END MODULE number_droplet_mod
