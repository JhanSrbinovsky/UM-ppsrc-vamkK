! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Cloud droplet number calculator
! Subroutine Interface:
MODULE lsp_taper_ndrop_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE lsp_taper_ndrop(                                             &
                   ! (Full) Aerosol tracers
                            so4_ait, so4_acc, so4_dis, sea_salt_film,   &
                            biogenic, sea_salt_jet, bmass_agd,          &
                            bmass_cld, ocff_agd, ocff_cld,              &
                            nitr_acc, nitr_diss, arcl,                  &
                   ! Murk aerosol
                            aerosol,                                    &
                   ! CDNC from UKCA 
                            ukca_cdnc,                                  &
                            cdnc_dim1, cdnc_dim2, cdnc_dim3,            &
                   ! Other parameters
                            rhodz_dry, rhodz_moist,                     &
                            deltaz,                                     &
                            snow_depth, land_fract,                     &
                   ! Output parameter of n_drop_tpr
                            n_drop_tpr                                  &
                                 )
! Microphysics modules

  USE mphys_inputs_mod,      ONLY: z_peak_nd, ndrop_surf,               &
                                   l_autoconv_murk, l_taper_new,        &
                                   max_drop_surf, l_droplet_tpr,        &
                                   l_mcr_arcl, l_clark_aero,            &
                                   arcl_inhom_sc

  USE mphys_constants_mod,   ONLY: ntot_land, ntot_sea, max_drop,       &
                                   n0_murk, m0_murk

  USE mphys_bypass_mod,      ONLY: mphys_mod_top
 
  USE lsp_autoc_consts_mod,  ONLY: power_murk, z_low_nd, min_drop_alt

! Dynamics and grid bounds modules

  USE level_heights_mod,     ONLY: eta_theta_levels
  USE atm_fields_bounds_mod, ONLY: qdims

! General constants modules
  USE arcl_mod,              ONLY: npd_arcl_compnts, ip_arcl_sulp_ak,   &
                                   ip_arcl_sulp_ac,  ip_arcl_sulp_di,   &
                                   ip_arcl_sslt_fi,  ip_arcl_sslt_jt,   &
                                   ip_arcl_biom_ag, ip_arcl_biom_ic,    &
                                   ip_arcl_ocff_ag, ip_arcl_ocff_ic

  USE conversions_mod,       ONLY: pi
  USE um_input_control_mod,  ONLY: l_mr_physics1,                       &
                                   l_use_sulphate_autoconv,             &
                                   l_use_bmass_autoconv,                &
                                   l_use_ocff_autoconv,                 &
                                   l_use_nitrate_autoconv,              &
                                   l_use_seasalt_autoconv,              &
                                   l_use_biogenic, l_use_arclbiom,      &
                                   l_use_arclsslt, l_use_arclsulp,      &
                                   l_use_arclocff 
  USE ukca_option_mod, ONLY:       l_ukca_aie2
  USE murk_inputs_mod,       ONLY: l_murk
  USE vectlib_mod,           ONLY: exp_v, powr_v
  USE number_droplet_mod,    ONLY: number_droplet

! Dr Hook modules
!--------------------------------------
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
!--------------------------------------

  IMPLICIT NONE

!---------------------------------------------------------------------
! Purpose:
!   Calculates cloud drop number concentration for the microphysics
!   scheme.
 
!   By using the tapering method, we can also reduce droplet number 
!   in the atmospheric boundary layer, from a peak at a given altitude 
!   (z_peak_nd) to a user defined value (ndrop_surf) or a variable
!   value dependent on aerosol amounts.

! Method:
!  1) Determine the height (specifically eta value) below which droplet
!     number is to taper. When the taper curve is inactive, this value
!     is not used.

!  2) Above this height, calculate droplet number using the
!     Haywood-Jones formulae for Murk aerosol
!     (as in the autoconversion routine).

!     Haywood et al (2008): (aerosol number) :
!              n_aer = n0_murk * ( Aerosol / m0_murk ) ** power_murk

!     Jones  et al (1994): Aerosol number to droplet number:
!              n_d   = 3.75e8 * ( 1 - exp ( -2.5e-9 * n_aer ) )

!     If full prognostic aerosols are used or aerosol climatologies
!     are used, then the routine calls the number_droplet function
!     to generate cloud droplet number concentration.

!     If no aerosol species are available and droplet tapering is
!     requested, then the routine uses a simple profile of
!     pre-determined values for droplet number.

!  3) Using the first height above z_peak_nd, determine the droplet
!     profile at each level towards the surface

!    n_drop = ndrop_surf + vala * log ( z / z(1) )

!    where vala = ( nd(z_peak_nd) - ndrop_surf ) / log ( z_peak_nd / z(1))

!    ndrop_surf may be fixed (from user input) or be calculated as a
!    variable value if l_taper_new is set to .TRUE.. In this latter case,
!    n_drop = ndrop_surf2 + vala * log ( z / z(1) )

!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

!------------------------------------------------
! Subroutine arguments
!------------------------------------------------

  INTEGER :: cdnc_dim1, cdnc_dim2, cdnc_dim3
 
  REAL, INTENT(IN) :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3) 
! CDNC from UKCA

  REAL, INTENT(IN) ::                                                   &
!-----------
! Aerosols
!-----------
  so4_ait      ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! SO4
! Aitken
  so4_acc      ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! SO4
! acc
  so4_dis      ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! SO4
! dis
  sea_salt_film( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! Sea salt
! film
  biogenic     ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! Biogenic
! aerosol
  sea_salt_jet ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! Sea salt
! Jet
  bmass_agd    ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! Aged
! Biomass
  bmass_cld    ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! Cloudy
! Biomass
  ocff_agd     ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! ocff
! aged
  ocff_cld     ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! ocff
! cloud
  nitr_acc     ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! Nitr
! aged
  nitr_diss    ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
!nitrate
!diss
  aerosol      ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
! Murk aerosol for
! droplet number
! calculations
  arcl         ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end,                           &
                            npd_arcl_compnts ),                         &
! Aerosol climatologies for droplet number calculations

  rhodz_dry    ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
  rhodz_moist  ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
  deltaz       ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                             1 : qdims%k_end ),                         &
!Grid box information

  land_fract   ( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end ),                         &
! Land fraction
!(for Jones et al,
! 1994 method)
  snow_depth   ( qdims%i_start : qdims%i_end ,                          &
                 qdims%j_start : qdims%j_end )
! Depth of snow
!(for Jones et al,
! 1994 method)

  REAL, INTENT(OUT) ::                                                  &
! intended output
    n_drop_tpr( qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end, 1:qdims%k_end )
! Droplet number after tapering has taken place

!------------------------------------------------
! Local variables
!------------------------------------------------

  REAL ::                                                               &
        eta_peak,                                                       &
           ! Eta value below which to taper the droplet number
        eta_low_nd,                                                     &
           ! Eta value above which to use low droplet number
        eta_before_taper,                                               &
           ! Eta value at the grid box immediately above the tapering
           ! level
        n_aer( qdims%i_start : qdims%i_end,                             &
               qdims%j_start : qdims%j_end  ),                          &
           ! Aerosol number from
           ! Haywood et al (2008)
        vala(        qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end ),                     &
        vala2,                                                          &
           ! Droplet number determined constants
        ndrop_surf2( qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end ),                     &
           ! Variable surface droplet number        
        half_range
           ! For cosine function, this is equivalent to
           ! ( 375.0E6 m 3 - min_drop_alt ) / 2.0
           ! Where 375.0E6 / m 3 is the maximum droplet number assigned
           ! by the Jones et al (1994) paper.

! 3D version of Air density in kg/m3
  REAL :: rho3d(  qdims%i_start : qdims%i_end,                          &
                  qdims%j_start : qdims%j_end,                          &
                  1 : qdims%k_end )

! A dummy array
  REAL :: dummy( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end,                           &
                 1 : qdims%k_end)

  INTEGER :: i, j, k              ! Loop counters
  INTEGER :: vec_length 
  INTEGER :: level_peak           ! level of eta_peak
  

! Logical to set nitrate climatology. Currently hardwired to .false.
! as a nitrate climatology is not yet available. 
  LOGICAL, PARAMETER :: l_use_arclnitr = .FALSE.

  REAL :: temp1, temp2
! temp for plane invariant reciprocals 
  REAL ::  tempv( qdims%i_start : qdims%i_end,                          &
       qdims%j_start : qdims%j_end )                          
! temp for output of vector computations

!------------------------------------------------
! Dr Hook subroutine timer details
!------------------------------------------------
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('LSP_TAPER_NDROP',zhook_in,zhook_handle)
!---------------------------------------------------------------------
! Start of physics
!---------------------------------------------------------------------

!----------------------------------------------------
! 1) Determine the eta level below which droplet number
!    is to taper
!----------------------------------------------------

  eta_peak    = z_peak_nd / mphys_mod_top
  eta_low_nd  = z_low_nd  / mphys_mod_top
  
  vec_length  = (qdims%i_end - qdims%i_start + 1)*(qdims%j_end -        &
       qdims%j_start + 1)

  IF (l_droplet_tpr) THEN

    !------------------------------------------------
    ! Droplet tapering is on, so we need to calculate
    ! droplet number concentration and add in a taper
    !------------------------------------------------

    ! First find index of level of eta_peak
    DO k = 1, qdims%k_end
      IF (eta_theta_levels(k) >= eta_peak) THEN
        level_peak = k
        eta_before_taper = eta_theta_levels(k)
        EXIT
      END IF
    END DO

    IF (l_murk .AND. l_autoconv_murk) THEN

      !----------------------------------------------
      ! Use murk aerosol to calculate droplet number
      ! and taper this profile
      !----------------------------------------------
     
      ! Calculate variable surface droplet number
      ! if required

      IF (l_taper_new) THEN

       DO j = qdims%j_start, qdims%j_end
         DO i = qdims%i_start, qdims%i_end

             ! Use Haywood et al (2008) formulae to get droplet number

             n_aer(i,j) = MAX( aerosol(i,j,1) / m0_murk*1.0e-9, 0.0001)
             ! 1.0E-9 converts from ug/kg to kg/kg

             !-----------------------------------------------
             ! Calculation of the aerosol number
             !-----------------------------------------------
             n_aer(i,j) = n0_murk * n_aer(i,j) ** power_murk

             !-----------------------------------------------
             ! Convert to CCN using a Jones et al (1994)
             ! relationship (as modified by Jonathan Wilkinson)
             !-----------------------------------------------

             ndrop_surf2(i,j) =                                        &
                max_drop_surf * (1.0 - EXP( - 1.5e-9 * n_aer(i,j) ) )

             !------------------------------------------------
             ! Ensure the surface droplet number doesn't get
             ! below the minimum value
             !------------------------------------------------

             ndrop_surf2(i,j) = MAX( ndrop_surf2(i,j), ndrop_surf)


         END DO ! qdims%i
       END DO   ! qdims%j

     END IF ! l_taper_new (new version of taper code)


     DO k = level_peak,qdims%k_end

       DO j = qdims%j_start, qdims%j_end
         DO i = qdims%i_start, qdims%i_end

          ! Above the taper level, so use Haywood-Jones formulae:

          n_aer(i,j) = MAX( aerosol(i,j, k) / m0_murk * 1.0e-9, 0.0001)
          ! 1.0E-9 converts from ug/kg to kg/kg

          END DO 
       END DO

       !-----------------------------------------------
       ! Calculation of the aerosol number
       !-----------------------------------------------
          
       CALL POWR_V(vec_length, n_aer, power_murk, tempv)

       !-----------------------------------------------
       ! Convert to CCN using the Jones et al (1994)
       ! relationship (follows number_droplet routine
       ! but reproduced here to avoid compiler issues).
       !-----------------------------------------------

       DO j = qdims%j_start, qdims%j_end
         DO i = qdims%i_start, qdims%i_end 
           n_aer(i,j) = -2.5e-9 *  n0_murk  * tempv(i,j)
         END DO
       END DO

       CALL EXP_V(vec_length, n_aer, tempv)

       DO j = qdims%j_start, qdims%j_end
         DO i = qdims%i_start, qdims%i_end 
           n_drop_tpr(i,j,k) = max_drop * ( 1.0e+00 - tempv(i,j)  )
         END DO
       END DO
       
       DO j = qdims%j_start, qdims%j_end
         DO i = qdims%i_start, qdims%i_end
              
           IF ( land_fract(i,j) >= 0.2 .AND.                          &
                snow_depth(i,j) < 5000.0    )    THEN
                
             n_drop_tpr(i,j,k) = MAX( n_drop_tpr(i,j,k), 35.0e+06 )

           ELSE
               
             n_drop_tpr(i,j,k) = MAX( n_drop_tpr(i,j,k), 5.0e+06  )
                
           END IF ! land fract > 0.2 etc
              
         END DO
       END DO
     END DO


     ! Start tapering
     DO k=1, level_peak-1
          
       temp1 = LOG( eta_before_taper / eta_theta_levels( 1 ) )
       temp1 = 1. / temp1
       temp2 =  LOG( eta_theta_levels(k) / eta_theta_levels(1) )
          
       IF (l_taper_new) THEN
         DO j = qdims%j_start, qdims%j_end      
           DO i = qdims%i_start, qdims%i_end             
                
             ! Use variable droplet number at surface
                
             vala( i, j ) = (n_drop_tpr(i,j,level_peak)              &
                  - ndrop_surf2(i,j)) /                              &
                  LOG( eta_before_taper /                            &
                  eta_theta_levels( 1 ) )
                
             n_drop_tpr( i, j, k ) =  ndrop_surf2(i,j)               &
                  + vala( i, j ) *                                   &
                  LOG( eta_theta_levels(k) /                         &
                  eta_theta_levels(1) )
                
                
           END DO ! qdims%i
         END DO   ! qdims%j
            
       ELSE  ! l_taper_new 
            
            ! Used fixed droplet number at surface
            
         DO j = qdims%j_start, qdims%j_end      
           DO i = qdims%i_start, qdims%i_end   
                
             vala( i, j ) = (n_drop_tpr(i,j,level_peak) - ndrop_surf) * &
                  temp1
                
             n_drop_tpr( i, j, k ) =  ndrop_surf + vala( i, j ) *    &
                  temp2

           END DO ! qdims%i
         END DO   ! qdims%j
            
       END IF

     END DO ! eta values below peak
        

   ELSE IF (l_mcr_arcl) THEN

   !------------------------------------------------------------
   ! Full droplet number calculation with aerosol climatologies
   !------------------------------------------------------------
      DO k = level_peak, qdims%k_end
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end

! Step 1: Calculate air density, rho

              IF (l_mr_physics1) THEN

                ! rho is the dry density
                rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)

              ELSE ! l_mr_physics1

                ! rho is the moist density
                rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)

              END IF  ! l_mr_physics1
          END DO
        END DO
      END DO

! Step 2: Call number_droplet routine
      CALL number_droplet(                                              &
                           qdims%i_start, qdims%i_end,                  &
                           qdims%j_start, qdims%j_end,                  &
                           1, qdims%k_end,                              &
                           level_peak, qdims%k_end,                     &
                           l_mcr_arcl,                                  &
                           .FALSE.,                                     &
                           arcl(:,:,:, ip_arcl_sulp_ac),                &
                           arcl(:,:,:, ip_arcl_sulp_di),                &
                           l_use_arclsslt,                              &
                           arcl(:,:,:, ip_arcl_sslt_fi),                &
                           arcl(:,:,:, ip_arcl_sslt_jt),                &
                           l_use_biogenic, biogenic,                    &
                           l_use_arclbiom,                              &
                           arcl(:,:,:, ip_arcl_biom_ag),                &
                           arcl(:,:,:, ip_arcl_biom_ic),                &
                           l_use_arclocff,                              &
                           arcl(:,:,:, ip_arcl_ocff_ag),                &
                           arcl(:,:,:, ip_arcl_ocff_ic),                &
                           l_use_arclnitr,                              &
                           dummy,                                       &
                           dummy,                                       &
                           rho3d,                                       &
                           snow_depth, land_fract,                      &
                           ntot_land, ntot_sea,                         &
                           n_drop_tpr)

      DO k = level_peak, qdims%k_end
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
! scaling the n_drop by the inhomogeneity scaling:
              n_drop_tpr(i,j,k) = n_drop_tpr(i,j,k) * arcl_inhom_sc

          END DO
        END DO
      END DO

      DO k = 1, level_peak-1
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
              vala( i, j ) = (n_drop_tpr(i,j,level_peak ) - ndrop_surf ) /&
                   LOG( eta_before_taper / eta_theta_levels( 1 ) )

              n_drop_tpr( i, j, k ) =  ndrop_surf + vala( i, j ) *      &
                                   LOG( eta_theta_levels(k) /           &
                                        eta_theta_levels(1) )

          END DO ! qdims%i
        END DO ! qdims%j
      END DO ! qdims%k

    ELSE IF (l_use_sulphate_autoconv .OR. l_ukca_aie2) THEN 
 
   !------------------------------------------------------------ 
   ! Full droplet number calculation based on prognostic aerosol 
   !------------------------------------------------------------

       DO k = level_peak, qdims%k_end
         DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end

              IF (l_ukca_aie2) THEN 
 
                n_drop_tpr(i,j,k) = ukca_cdnc(i,j,k) 
 
              ELSE  ! (not l_ukca_aie2)

! Step 1: Calculate air density, rho

                IF (l_mr_physics1) THEN 
 
                  ! rho is the dry density 
   
                  rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k) 
    
                ELSE ! l_mr_physics1 
    
                  ! rho is the moist density 
 
                  rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k) 
  
                END IF  ! l_mr_physics1 
              END IF  ! l_ukca_aie2
            END DO
          END DO
        END DO

! Step 2: Call number_droplet routine
        IF (.NOT. l_ukca_aie2) THEN

          CALL number_droplet(                                          &
                           qdims%i_start, qdims%i_end,                  &
                           qdims%j_start, qdims%j_end,                  &
                           1, qdims%k_end,                              &
                           level_peak, qdims%k_end,                     &
                           l_use_sulphate_autoconv,                     &
                           .FALSE.,                                     &
                           so4_acc,                                     &
                           so4_dis,                                     &
                           l_use_seasalt_autoconv,                      &
                           sea_salt_film,                               &
                           sea_salt_jet,                                &
                           l_use_biogenic, biogenic,                    &
                           l_use_bmass_autoconv,                        &
                           bmass_agd,                                   &
                           bmass_cld,                                   &
                           l_use_ocff_autoconv,                         &
                           ocff_agd, ocff_cld,                          &
                           l_use_nitrate_autoconv,                      &
                           nitr_acc,                                    &
                           nitr_diss,                                   &
                           rho3d,                                       &
                           snow_depth, land_fract,                      &
                           ntot_land, ntot_sea,                         &
                           n_drop_tpr)
      END IF


      DO k = 1, level_peak-1
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end

              ! start tapering
              vala( i, j ) = (n_drop_tpr( i, j, level_peak ) - ndrop_surf ) /&
                   LOG( eta_before_taper / eta_theta_levels( 1 ) )

              n_drop_tpr( i, j, k ) =  ndrop_surf + vala( i, j ) *      &
                                   LOG( eta_theta_levels(k) /           &
                                        eta_theta_levels(1) )

          END DO ! qdims%i
        END DO ! qdims%j
      END DO ! qdims%k

   ELSE   ! not l_murk, l_mcr_arcl, l_use_sulphate_autoconv or l_ukca_aie2

     !-----------------------------------------
     ! Use a simple profile for tapering
     !-----------------------------------------

     ! Do not allow a peak value of eta to exceed that of the
     ! pre-determined low droplet number (usually around 2km)

     IF (eta_peak > eta_low_nd) eta_peak = eta_low_nd

       DO k = qdims%k_end, 1, -1

         IF ( eta_theta_levels(k) >= eta_low_nd .AND.                   &
              eta_theta_levels(k) >= eta_peak         ) THEN

           ! Above taper level and level set for minimum droplet
           ! number, hence set n_drop at this level to low value

           n_drop_tpr(:,:,k) = 100.0e6


         ELSE IF ( eta_theta_levels(k) >= eta_peak .AND.                &
                   eta_theta_levels(k) < eta_low_nd     ) THEN

          ! Above taper level yet below level set for minimum
          ! droplet number, so use a simple function of eta to
          ! determine the droplet number. This is intended to
          ! be independent of model levels

           half_range = (max_drop - min_drop_alt) / 2.0

           n_drop_tpr(:,:,k) = half_range * (- COS (pi + pi *           &
                             ( ( eta_theta_levels (k) - eta_peak ) /    &
                             ( eta_low_nd - eta_peak ) ) ) )            &
                             + half_range + min_drop_alt

          ! No need to set up peak droplet as this should be max_drop

        ELSE ! eta_theta_levels


          vala2  =  ( max_drop - ndrop_surf ) /                         &
                    ( LOG( eta_peak / eta_theta_levels( 1 ) ) )

          n_drop_tpr( :,:, k ) = ndrop_surf + vala2 *                   &
                                 LOG( eta_theta_levels(k) /             &
                                      eta_theta_levels(1) )

        END IF !eta_levels above or below peak values

      END DO  ! qdims%k

    END IF ! l_murk / l_use_sulphate_autoconv

  ELSE ! l_droplet_tpr

    !-------------------------------------------------
    ! Droplet tapering is not active, but we still
    ! need to calculate potential cloud drop number
    ! concentration for use in autoconversion
    !-------------------------------------------------

    IF (l_murk .AND. l_autoconv_murk) THEN

      ! Calculate using MURK aerosol
      DO k = 1,qdims%k_end
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end

            ! Use Jones-Haywood Formulae
            n_aer(i,j) =                                                &
                       MAX( aerosol(i,j, k) / m0_murk * 1.0e-9, 0.0001)
            ! 1.0E-9 converts from ug/kg to kg/kg

            !-----------------------------------------------
            ! Calculation of the aerosol number
            !-----------------------------------------------
            n_aer(i,j) = n0_murk * n_aer(i,j) ** power_murk
            n_drop_tpr(i,j,k) =                                         &
                max_drop * ( 1.0e+00-EXP( -2.5e-9 * n_aer(i,j) ) )

            IF ( land_fract(i,j) >= 0.2 .AND.                           &
                 snow_depth(i,j) < 5000.0    )    THEN

              n_drop_tpr(i,j,k) = MAX( n_drop_tpr(i,j,k), 35.0e+06 )

            ELSE

              n_drop_tpr(i,j,k) = MAX( n_drop_tpr(i,j,k), 5.0e+06  )

            END IF ! land fract > 0.2 etc

          END DO ! i (qdims%i)
        END DO ! j (qdims%j)
      END DO ! k (qdims%k)

    ELSE IF (l_mcr_arcl) THEN

      ! Calculate using aerosol climatologies

      DO k = 1,qdims%k_end
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end

! Step 1: Calculate air density, rho

              IF (l_mr_physics1) THEN

                ! rho is the dry density
                rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)

              ELSE ! l_mr_physics1

                ! rho is the moist density
                rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)

              END IF  ! l_mr_physics1
            END DO
          END DO
        END DO

! Step 2: Call number_droplet routine
      CALL number_droplet(                                              &
                           qdims%i_start, qdims%i_end,                  &
                           qdims%j_start, qdims%j_end,                  &
                           1, qdims%k_end,                              &
                           1, qdims%k_end,                              &
                           l_mcr_arcl,                                  &
                           .FALSE.,                                     &
                           arcl(:,:,:, ip_arcl_sulp_ac),                &
                           arcl(:,:,:, ip_arcl_sulp_di),                &
                           l_use_arclsslt,                              &
                           arcl(:,:,:, ip_arcl_sslt_fi),                &
                           arcl(:,:,:, ip_arcl_sslt_jt),                &
                           l_use_biogenic, biogenic,                    &
                           l_use_arclbiom,                              &
                           arcl(:,:,:, ip_arcl_biom_ag),                &
                           arcl(:,:,:, ip_arcl_biom_ic),                &
                           l_use_arclocff,                              &
                           arcl(:,:,:, ip_arcl_ocff_ag),                &
                           arcl(:,:,:, ip_arcl_ocff_ic),                &
                           l_use_arclnitr,                              &
                           dummy,                                       &
                           dummy,                                       &
                           rho3d,                                       &
                           snow_depth, land_fract,                      &
                           ntot_land, ntot_sea,                         &
                           n_drop_tpr)

      DO k = 1,qdims%k_end
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
! scaling the n_drop by the inhomogeneity scaling:
              n_drop_tpr(i,j,k) = n_drop_tpr(i,j,k) * arcl_inhom_sc

          END DO
        END DO
      END DO

    ELSE IF (l_use_sulphate_autoconv .OR. l_ukca_aie2) THEN 
 
   !------------------------------------------------------------ 
   ! Full droplet number calculation based on prognostic aerosol 
   !------------------------------------------------------------

      DO k = 1, qdims%k_end
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end

              IF (l_ukca_aie2) THEN 
 
                n_drop_tpr(i,j,k) = ukca_cdnc(i,j,k) 
 
              ELSE  ! (not l_ukca_aie2)

! Step 1: Calculate air density, rho

                IF (l_mr_physics1) THEN 
 
                  ! rho is the dry density 
 
                  rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k) 
 
                ELSE ! l_mr_physics1 
 
                  ! rho is the moist density 
 
                  rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k) 
 
                END IF  ! l_mr_physics1 

              END IF  ! l_ukca_aie2

          END DO
        END DO
      END DO

    IF (.NOT. l_ukca_aie2) THEN

! Step 2: Call number_droplet routine
        CALL number_droplet(                                            &
                           qdims % i_start, qdims % i_end,              &
                           qdims % j_start, qdims % j_end,              &
                           1, qdims % k_end,                            &
                           1, qdims % k_end,                            &
                           l_use_sulphate_autoconv,                     &
                           .FALSE.,                                     &
                           so4_acc,                                     &
                           so4_dis,                                     &
                           l_use_seasalt_autoconv,                      &
                           sea_salt_film,                               &
                           sea_salt_jet,                                &
                           l_use_biogenic, biogenic,                    &
                           l_use_bmass_autoconv,                        &
                           bmass_agd,                                   &
                           bmass_cld,                                   &
                           l_use_ocff_autoconv,                         &
                           ocff_agd, ocff_cld,                          &
                           l_use_nitrate_autoconv,                      &
                           nitr_acc,                                    &
                           nitr_diss,                                   &
                           rho3d,                                       &
                           snow_depth, land_fract,                      &
                           ntot_land, ntot_sea,                         &
                           n_drop_tpr)
    END IF ! l_ukca_aie2

    ELSE ! (not l_murk, l_mcr_arcl, l_use_sulphate_autoconv or l_ukca_aie2)

      !-----------------------------------------------------------
      ! In this case, STASH 4/211 will be a simple land-sea split
      ! dependant on ntot_land and ntot_sea
      !-----------------------------------------------------------

      ! Set the first level
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end

          IF (land_fract(i,j) >= 0.5) THEN
            n_drop_tpr(i,j,1) = ntot_land
          ELSE ! land_fract
            n_drop_tpr(i,j,1) = ntot_sea
          END IF ! land fract

        END DO ! i (qdims%i)
      END DO ! j (qdims%j)

      ! copy up to higher levels - thus avoiding too many branches
      ! in loops
      DO k = 2, qdims%k_end
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
            n_drop_tpr(i,j,k) = n_drop_tpr(i,j,1)
          END DO
        END DO
      END DO

    END IF ! l_autoconv_murk, l_use_sulphate_autoconv
   
  END IF ! l_droplet_tpr

!---------------------------------------------------------------------
! End of physics
!---------------------------------------------------------------------
  IF (lhook) CALL dr_hook('LSP_TAPER_NDROP',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_taper_ndrop
END MODULE lsp_taper_ndrop_mod
