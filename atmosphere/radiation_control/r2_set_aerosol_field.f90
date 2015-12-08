! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_AEROSOL_FIELD(IERR                              &
     &   , N_PROFILE, NLEVS, N_LAYER, N_AEROSOL, N_AEROSOL_MR           &
     &   , TYPE_AEROSOL, I_GATHER, L_EXTRA_TOP                          &
     &   , L_CLIMAT_AEROSOL, L_CLIM_AERO_HGT, L_HadGEM1_Clim_Aero       &
     &   , BL_DEPTH, T, N_LEVELS_BL, L_MURK_RAD, AERO_MESO              &
     &   , L_DUST, L_USE_DUST, DUST_DIM1, DUST_DIM2                     &
     &   , DUST_1, DUST_2, DUST_3, DUST_4, DUST_5, DUST_6               &
     &   , L_USE_BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2, BIOGENIC       &
     &   , L_SULPC_SO2, L_USE_SULPC_DIRECT                              &
     &   , SULP_DIM1, SULP_DIM2                                         &
     &   , ACCUM_SULPHATE, AITKEN_SULPHATE                              &
     &   , L_USE_SEASALT_DIRECT, SALT_DIM_A, SALT_DIM_B                 &
     &   , SEA_SALT_FILM, SEA_SALT_JET, P                               &
     &   , L_SOOT, L_USE_SOOT_DIRECT, SOOT_DIM1, SOOT_DIM2              &
     &   , FRESH_SOOT, AGED_SOOT                                        &
     &   , L_BIOMASS, L_USE_BMASS_DIRECT, BMASS_DIM1, BMASS_DIM2        &
     &   , FRESH_BMASS, AGED_BMASS                                      &
     &   , L_OCFF, L_USE_OCFF_DIRECT, OCFF_DIM1, OCFF_DIM2              &
     &   , FRESH_OCFF, AGED_OCFF                                        &
     &   , L_NITRATE, L_USE_NITRATE_DIRECT, NITRATE_DIM1, NITRATE_DIM2  &
     &   , ACCUM_NITRATE                                                &
     &   , N_ARCL_SPECIES, N_ARCL_COMPNTS, I_ARCL_COMPNTS               &
     &   , L_USE_ARCL, ARCL_DIM1, ARCL_DIM2, ARCL                       &
     &   , LAND, LYING_SNOW, PSTAR                                      &
     &   , P_LAYER_BOUNDARIES                                           &
     &   , TRINDX, alat, PREVIOUS_TIME                                  &
     &   , AEROSOL_MIX_RATIO, AEROSOL_MR_SOURCE, AEROSOL_MR_TYPE_INDEX  &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES       &
     &   , NPD_AEROSOL_MIXRATIO, FIRST_LAYER                            &
     &   )
!
!
!
      USE atmos_constants_mod, ONLY: r
     
      USE conversions_mod, ONLY: pi
      USE dust_parameters_mod, ONLY: l_twobin_dust
      USE arcl_mod,           ONLY: npd_arcl_compnts, npd_arcl_species, &
                              ip_arcl_sulp, ip_arcl_dust, ip_arcl_sslt, &
                              ip_arcl_blck, ip_arcl_biom, ip_arcl_ocff, &
                              ip_arcl_dlta, &
                              ip_arcl_sulp_ac, ip_arcl_sulp_ak, &
                              ip_arcl_sulp_di, ip_arcl_dust_b1, &
                              ip_arcl_dust_b2, ip_arcl_dust_b3, &
                              ip_arcl_dust_b4, ip_arcl_dust_b5, &
                              ip_arcl_dust_b6, ip_arcl_sslt_fi, &
                              ip_arcl_sslt_jt, ip_arcl_blck_fr, &
                              ip_arcl_blck_ag, ip_arcl_biom_fr, &
                              ip_arcl_biom_ag, ip_arcl_biom_ic, &
                              ip_arcl_ocff_fr, ip_arcl_ocff_ag, &
                              ip_arcl_ocff_ic, ip_arcl_dlta_dl
      USE rad_input_mod, ONLY: aeroscl_csk_clim
      USE rad_pcf, ONLY: i_err_fatal
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY: ereport
      USE aercmp3a_mod, ONLY: npd_aerosol_component, ip_water_soluble,  &
                              ip_dust_like, ip_oceanic, ip_soot, ip_ash,&
                              ip_sulphuric, ip_accum_sulphate,          &
                              ip_aitken_sulphate, ip_fresh_soot,        &
                              ip_aged_soot, ip_seasalt_film,            &
                              ip_seasalt_jet, ip_dust_1, ip_dust_2,     &
                              ip_dust_3, ip_dust_4, ip_dust_5,          &
                              ip_dust_6, ip_biomass_1, ip_biomass_2,    &
                              ip_biogenic, ip_ocff_fresh, ip_ocff_aged, &
                              ip_delta, ip_nitrate, ip_twobindust_1,    &
                              ip_twobindust_2, ip_aersrc_cusack_ron,    &
                              ip_aersrc_cusack_roff,                    &
                              ip_aersrc_classic_ron,                    &
                              ip_aersrc_classic_roff,                   &
                              ip_aersrc_arcl_ron, ip_aersrc_arcl_roff
      USE stdio3a_mod, ONLY: iu_stdin, iu_stdout, iu_err
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
!
!     DUMMY ARGUMENTS.
!
!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE                                                  &
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES                                          &
!             MAXIMUM NUMBER OF AEROSOL SPECIES IN SPECTRAL INFORMATION
     &   , NPD_AEROSOL_MIXRATIO                                         &
!             MAXIMUM NUMBER OF AEROSOL SPECIES MIXING RATIO INFORMATION
     &   , FIRST_LAYER
!             First layer for some variables
!             0 for flux_calc(3A/C), 1 for radiance_calc(3Z)
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , NLEVS                                                        &
!             Number of layers used outside the radiation scheme
     &   , N_LAYER                                                      &
!             Number of layers seen by radiation
     &   , N_LEVELS_BL                                                  &
!             Number of layers occupied by boundary-layer aerosol
!             if L_CLIM_AERO_HGT is false.
     &   , N_AEROSOL                                                    &
!             NUMBER OF AEROSOLS IN SPECTRAL FILE
     &   , N_AEROSOL_MR                                                 &
!             NUMBER OF AEROSOLS IN AEROSOL_MIXING_RATIO ARRAY
     &   , TYPE_AEROSOL(NPD_AEROSOL_SPECIES)
!             ACTUAL TYPES OF AEROSOLS
!
!     GATHERING ARRAY:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO GATHER
!
!     FLAG FOR THE CLIMATOLOGICAL AEROSOL DISTRIBUTION.
      LOGICAL                                                           &
                   !, INTENT(IN)
     &     L_CLIMAT_AEROSOL                                             &
!             FLAG FOR CLIMATOLOGICAL AEROSOL DISTRIBUTION
     &   , L_CLIM_AERO_HGT                                              &
!             Flag to use the boundary layer depth in setting the
!             climatological aerosol
     &   , L_HadGEM1_Clim_Aero                                          &
!             Flag to use HadGEM1 setting for climatological aerosols
     &   , L_MURK_RAD                                                   &
!             FLAG FOR MESOSCALE MODEL AEROSOL
     &   , L_EXTRA_TOP
!             Flag to include an extra top layer in the radiation scheme
!
! Declare mineral dust aerosol variables:
      LOGICAL                                                           &
     &     L_DUST                                                       &
!             Mineral dust is available for radn. effects or diagnositcs
     &   , L_USE_DUST
!             Use direct radiative effect of mineral dust
      INTEGER DUST_DIM1,DUST_DIM2
!        Dimensions for mineral dust arrays (P_FIELD,P_LEVELS or 1,1)
      REAL DUST_1(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div1 dust
     &   , DUST_2(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div2 dust
     &   , DUST_3(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div3 dust
     &   , DUST_4(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div4 dust
     &   , DUST_5(DUST_DIM1, DUST_DIM2)                                 &
                                        !Mass mixing ratio of div5 dust
     &   , DUST_6(DUST_DIM1, DUST_DIM2) !Mass mixing ratio of div6 dust
!     VARIABLES FOR THE BIOGENIC AEROSOL
      LOGICAL L_USE_BIOGENIC
      INTEGER BIOGENIC_DIM1, BIOGENIC_DIM2
      REAL BIOGENIC(BIOGENIC_DIM1, BIOGENIC_DIM2)

!     VARIABLES FOR THE SULPHUR CYCLE:
      LOGICAL                                                           &
                   !, INTENT(IN)
     &     L_SULPC_SO2                                                  &
!             SULPHUR CYCLE AVAILABLE FOR RADN. EFFECTS OR DIAGNOSTICS
     &   , L_USE_SULPC_DIRECT
!             FLAG TO USE SULPHUR CYCLE FOR DIRECT EFFECT
      INTEGER                                                           &
                   !, INTENT(IN)
     &     SULP_DIM1,SULP_DIM2
!             DIMENSIONS FOR _SULPHATE ARRAYS, (P_FIELD,P_LEVELS or 1,1)
      REAL                                                              &
                !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)                         &
!             MASS MIXING RATIOS OF ACCUMULATION MODE AEROSOL
     &   , AITKEN_SULPHATE(SULP_DIM1, SULP_DIM2)
!             MASS MIXING RATIOS OF AITKEN MODE AEROSOL
!
! Declare soot variables:
      LOGICAL                                                           &
           L_SOOT                                                       &
!             SOOT AEROSOL AVAILABLE FOR RADN. EFFECTS OR DIAGNOSTICS
         , L_USE_SOOT_DIRECT
!             USE DIRECT RAD. EFFECT OF SOOT AEROSOL
      INTEGER SOOT_DIM1,SOOT_DIM2
                !DIMENSIONS FOR SOOT ARRAYS, (P_FIELD,P_LEVELS or 1,1)
      REAL FRESH_SOOT(SOOT_DIM1, SOOT_DIM2)                             &
                                                 ! MMR OF FRESH SOOT
     &   , AGED_SOOT(SOOT_DIM1, SOOT_DIM2)       ! MMR OF AGED SOOT
!
! Declare biomass smoke aerosol variables:
      LOGICAL                                                           &
           L_BIOMASS                                                    &
!             Biomass smoke available for radn. effects or diagnostics
         , L_USE_BMASS_DIRECT
!              Use direct radiative effect of biomass smoke
      INTEGER BMASS_DIM1,BMASS_DIM2
!              Dimensions for biomass arrays (P_FIELD,P_LEVELS or 1,1)
      REAL FRESH_BMASS(BMASS_DIM1, BMASS_DIM2)                          &
                                                    ! MMR OF FRESH SMOKE
     &   , AGED_BMASS(BMASS_DIM1, BMASS_DIM2)       ! MMR OF AGED SMOKE
!
! Declare fossil-fuel organic carbon aerosol variables:
      LOGICAL                                                           &
           L_OCFF                                                       &
!              OCFF aerosol available for radn. effects or diagnostics
         , L_USE_OCFF_DIRECT
!              Use direct radiative effect of OCFF aerosol
      INTEGER OCFF_DIM1, OCFF_DIM2
!              Dimensions for OCFF arrays (P_FIELD,P_LEVELS or 1,1)
      REAL FRESH_OCFF(OCFF_DIM1, OCFF_DIM2)                             &
                                                    ! MMR OF FRESH OCFF
     &   , AGED_OCFF(OCFF_DIM1, OCFF_DIM2)          ! MMR OF AGED OCFF
!
      LOGICAL                                                           &
     &     L_USE_SEASALT_DIRECT
!             Flag for direct effect of interactive sea-salt aerosol
!
      INTEGER                                                           &
     &     SALT_DIM_A, SALT_DIM_B
!             Array sizes of sea-salt aerosols
!                (either P_FIELD,P_LEVELS or 1,1)
!
      REAL                                                              &
     &     SEA_SALT_FILM(SALT_DIM_A, SALT_DIM_B)                        &
!             On input, number concentration (m-3) of film-mode
!             sea-salt aerosols; converted to mass mixing ratio.
     &   , SEA_SALT_JET(SALT_DIM_A, SALT_DIM_B)
!             On input, number concentration (m-3) of jet-mode
!             sea-salt aerosols; converted to mass mixing ratio.
! Declare nitrate aerosol variables: 
      LOGICAL                                                           &
           L_NITRATE                                                    &
!              Nitrate aerosol available for radn. effects or diagnostcs
         , L_USE_NITRATE_DIRECT 
!              Use direct radiative effect of nitrate aerosols 
      INTEGER NITRATE_DIM1, NITRATE_DIM2 
!              Dimensions for nitrate arrays (P_FIELD,P_LEVELS or 1,1) 
      REAL ACCUM_NITRATE(NITRATE_DIM1, NITRATE_DIM2) 
!              Mass-mixing ratio of accumulation nitrate 
! 
!
! Aerosol climatology for NWP:
!

      ! Number of requested species within the climatology
      Integer N_ARCL_SPECIES
      
      ! Corresponding number of requested components
      Integer N_ARCL_COMPNTS
      
      ! Model switches for each species
      Logical, dimension(NPD_ARCL_SPECIES) :: L_USE_ARCL
      
      ! Array index of components
      Integer, dimension(NPD_ARCL_COMPNTS) :: I_ARCL_COMPNTS
      
      ! Array dimensions
      Integer ARCL_DIM1, ARCL_DIM2
      
      ! Mass-mixing ratios
      Real                                                              &
     &    ARCL(ARCL_DIM1, ARCL_DIM2, N_ARCL_COMPNTS)
      
      REAL                                                              &
     &     AERO_MESO(NPD_FIELD, NLEVS)
!             MIXING RATIO OF 'URBAN' AEROSOL OF MESOSCALE MODEL
!     GENERAL ATMOSPHERIC PROPERTIES:
      INTEGER                                                           &
                !, INTENT(IN)
     &     TRINDX(NPD_FIELD)                                            &
!             LAYER BOUNDARY OF TROPOPAUSE
     &,    PREVIOUS_TIME(7)
!             Time
      REAL                                                              &
                !, INTENT(IN)
     &     PSTAR(NPD_FIELD)                                             &
!             SURFACE PRESSURES
     &,    P_LAYER_BOUNDARIES(NPD_FIELD,0:NLEVS)                        &
!             PRESSURE AT BOUNDARIES OF LAYERS
     &,    alat(npd_profile)
!             Latitudes in degrees
!
!     SURFACE FIELDS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     LAND(NPD_FIELD)
!             LAND SEA MASK
      REAL                                                              &
                !, INTENT(IN)
     &     LYING_SNOW(NPD_FIELD)                                        &
!             DEPTH OF LYING SNOW
     &   , BL_DEPTH(NPD_FIELD)                                          &
!             Depth of the boundary layer
     &   , T(NPD_PROFILE, FIRST_LAYER:NPD_LAYER)                        &
!             Temperatures of atmospheric layers
     &   , P(NPD_PROFILE, FIRST_LAYER:NPD_LAYER)
!             Pressures at layer centres
!
      REAL                                                              &
                !, INTENT(OUT)
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, FIRST_LAYER:NPD_LAYER         &
     &        , NPD_AEROSOL_MIXRATIO)
!             MIXING RATIOS OF AEROSOLS
      INTEGER                                                           &
                   !, INTENT(OUT)
     &     AEROSOL_MR_TYPE_INDEX(NPD_AEROSOL_MIXRATIO)                  &
!             Index relating aerosol_mix_ratio aerosols to aerosols in
!             the spectral information
     &   , AEROSOL_MR_SOURCE(NPD_AEROSOL_MIXRATIO)
!             Source of the aerosol data (compared to IP_AERSRC_... )
!
!
!     LOCAL VARIABLES:
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , J_PLUS                                                       &
!             ADDITION TO LOOP VARIABLE, FOR INDEXING AEROSOL_MIX_RATIO
     &   , J_MR                                                         &
!             LOOP VARIABLE FOR LOOPING OVER J+J_PLUS
     &   , L                                                            &
!             LOOP VARIABLE
     &   , LG                                                           &
!             INDEX FOR GATHERING
     &   , BLTOP                                                        &
!             INDEX OF UPPER BOUNDARY OF PLANETARY BOUNDARY LAYER
     &   , I_AEROSOL                                                    &
!             ACTUAL TYPE OF AEROSOL BEING CONSIDERED
     &   , I_TOP_COPY
!             Topmost radiative layer to be set directly from the
!             profile supplied
!
!
!     ARRAYS FOR THE CLIMATOLOGICAL AEROSOL MODEL
      Integer :: ii
!       Variable used in initialization of arrays.
      Logical, Dimension(NPD_AEROSOL_COMPONENT) :: L_IN_CLIMAT

!       Flags to indicate which aerosols are included in the
!       climatology: this may be used to enable various components
!       to be replaced by fully prognostic schemes.
      Integer, Dimension(NPD_AEROSOL_COMPONENT) :: I_CLIM_POINTER =     &
     &    (/ 1, 2, 3, 4, 0, 5,                                          &
     &    (0, ii=7, NPD_AEROSOL_COMPONENT) /)

!       Pointers to the indices of the original climatological
!       aerosol model.
      REAL                                                              &
     &     AEROSOL_MIX_RATIO_CLIM(NPD_PROFILE, FIRST_LAYER:NPD_LAYER, 5)
!             MIXING RATIOS OF THE CLIMATOLOGICAL AEROSOLS
!
      REAL                                                              &
     &     MODE_RADIUS_SS_FILM                                          &
!            Mode radius of film-mode sea-salt aerosol (m)
     &   , MODE_RADIUS_SS_JET                                           &
!            Mode radius of jet-mode sea-salt aerosol (m)
     &   , SIGMA_SS_FILM                                                &
!            Geometric standard deviation of film-mode sea-salt aerosol
     &   , SIGMA_SS_JET                                                 &
!            Geometric standard deviation of jet-mode sea-salt aerosol
     &   , DENSITY_SEA_SALT
!            Bulk density of sodium chloride, taken as representative
!            of sea-salt aerosol (kg m-3); taken from Pruppacher & Klett
!            2nd edition (1997), p. 244.
!
      Integer im
      Real, Dimension(NPD_AEROSOL_COMPONENT) :: Meso_frac =             &
     &  (/ 0.61, 0.17, 0.0, 0.22, (0.0, im=5,NPD_AEROSOL_COMPONENT) /)
      REAL INV17

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*),PARAMETER   :: RoutineName = 'r2_set_aerosol_field'
      CHARACTER(LEN=256) :: cmessage

      PARAMETER(                                                        &
     &     MODE_RADIUS_SS_FILM=0.1E-06                                  &
     &   , MODE_RADIUS_SS_JET=1.0E-06                                   &
     &   , SIGMA_SS_FILM=1.9                                            &
     &   , SIGMA_SS_JET=2.0                                             &
     &   , DENSITY_SEA_SALT=2165.0                                      &
     &   )
!
!

!
      IF (lhook) CALL dr_hook('R2_SET_AEROSOL_FIELD',zhook_in,zhook_handle)
!
!     INITIALISE SOME VALUES:
      J_PLUS=0
!
!     If an extra layer is used in the radiation scheme, any
!     aerosol profiles supplied will be copied into the profiles
!     seen by the radiation scheme starting with the second layer
!     down, rather than the first.
      I_TOP_COPY=1
      IF (L_EXTRA_TOP) I_TOP_COPY=2


      L_IN_CLIMAT(1:4) = .TRUE.
      L_IN_CLIMAT(5)   = .FALSE.
      L_IN_CLIMAT(6)   = .TRUE.
      L_IN_CLIMAT(7:NPD_AEROSOL_COMPONENT) = .FALSE.

!  Use HadGEM1 settings to remove boundary layer climatological aerosols.
      IF (L_HadGEM1_Clim_Aero) THEN
        DO ii=1,N_AEROSOL
          IF ((TYPE_AEROSOL(ii) == IP_WATER_SOLUBLE).OR.                &
     &        (TYPE_AEROSOL(ii) == IP_DUST_LIKE) .OR.                   &
     &        (TYPE_AEROSOL(ii) == IP_OCEANIC).OR.                      &
     &        (TYPE_AEROSOL(ii) == IP_SOOT)) THEN
            L_IN_CLIMAT(ii) = .FALSE.
          END IF
        END DO
!  Check to see if MURK is on in radiation
        IF (L_MURK_RAD) THEN
          WRITE(cmessage, '(A)') '*** ERROR: MURK AEROSOL SWITCHED'      & 
            // 'ON IN RADIATION BUT THIS IS CURRENTLY INCONSISTENT WITH' & 
            // 'ONLY HAVING THE STRATOSPHERIC CUSACK CLIMATLOGY - THE'   & 
            // 'MURK SPLITS ITS MIXING RATIO AMONGST THE BOUNDARY LAYER' & 
            // 'CUSACK AEROSOLS'
          ierr=i_err_fatal 
          CALL ereport(RoutineName,ierr,cmessage)
          IF (lhook) CALL dr_hook('R2_SET_AEROSOL_FIELD',zhook_out,zhook_handle)
          RETURN
        ENDIF
      END IF
!
!  Use climatological soot if climatological aerosols are on and not
!  using interactive soot.
       L_IN_CLIMAT(IP_SOOT) = L_IN_CLIMAT(IP_SOOT)                      &
     &                     .AND.(.NOT.L_USE_SOOT_DIRECT)
!
!  Use climatological "oceanic" aerosol if climatological aerosols have
!  been selected and not using sea-salt parametrization.
      L_IN_CLIMAT(IP_OCEANIC) = L_IN_CLIMAT(IP_OCEANIC)                 &
     &                          .AND.(.NOT. L_USE_SEASALT_DIRECT)
!
!  The climatological water-soluble aerosol should not be used
!  if the direct effects of sulphate aerosols are included.
!  (Note that this differs the situation applying in earlier
!  versions of the model, such as HadAM3).
      L_IN_CLIMAT(IP_WATER_SOLUBLE) = L_IN_CLIMAT(IP_WATER_SOLUBLE)     &
     &                     .AND.(.NOT.L_USE_SULPC_DIRECT)
!
!
!
      IF (L_CLIMAT_AEROSOL) THEN
!
!        SET THE MIXING RATIOS OF THE CLIMATOLOGICAL AEROSOLS
!        USED IN THE CLIMATOLOGY OF HADCM3. A SEPARATE SUBROUTINE
!        IS USED TO ENSURE BIT-REPRODUCIBLE RESULTS BY USING
!        EARLIER CODE. THIS COULD BE ALTERED IF A NEW CLIMATOLOGY WERE
!        USED.
!
! DEPENDS ON: r2_set_aero_clim_hadcm3
         CALL R2_SET_AERO_CLIM_HADCM3(N_PROFILE, NLEVS, N_LAYER         &
     &      , I_GATHER, L_EXTRA_TOP                                     &
     &      , L_CLIM_AERO_HGT, BL_DEPTH, T, N_LEVELS_BL                 &
     &      , LAND, LYING_SNOW, PSTAR, P_LAYER_BOUNDARIES, TRINDX       &
     &      , alat, PREVIOUS_TIME                                       &
     &      , AEROSOL_MIX_RATIO_CLIM                                    &
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER, FIRST_LAYER            &
     &      )
!
      ENDIF
!
!
!     THE AEROSOLS REQUIRED BY FOR THE CALCULATION SHOULD HAVE BEEN
!     SELECTED WHEN THE SPECTRAL FILE WAS READ IN. EACH TYPE SHOULD
!     BE SET APPROPRIATELY.
!
      DO J=1, N_AEROSOL
!
         I_AEROSOL=TYPE_AEROSOL(J)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF (L_CLIMAT_AEROSOL.AND.L_IN_CLIMAT(I_AEROSOL)) THEN
            AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
            AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CUSACK_RON
!
!           Here mxing ratios in all layers are set because the
!           possible presence of an extra layer has been allowed
!           for in the subroutine.
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =AEROSOL_MIX_RATIO_CLIM(L, I                       &
                     , I_CLIM_POINTER(I_AEROSOL))                       & 
                     * aeroscl_csk_clim(i_clim_pointer(i_aerosol))
               ENDDO
            ENDDO
!
! For the other aerosol species, they may be supplied from dedicated
! schemes (e.g. climate model) or climatologies (e.g. NWP model).
! The former are indicated by L_USE_<SPECIES>, the latter by
! L_USE_ARCL(IP_ARCL_<SPECIES>). Should both logical be set to true
! for a given aerosol species (i.e. an aerosol is provided by a 
! dedicated scheme AND a climatology), the climatology wins and will 
! be seen by radiation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Mineral dust aerosols:
!
!        Six-bin dust, bin 1, clim or prognostic
         ELSE IF (i_aerosol == ip_dust_1 .AND.                          & 
                  (l_use_arcl(ip_arcl_dust).or.( l_dust .AND.           &
                   (.NOT.l_twobin_dust) ) ) ) THEN
           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              ! RECORD THE J FOR THE AEROSOL MIXING RATIO AS THIS IS THE
              ! INDEX THAT RELATES IT TO SPECTRAL INFORMATION
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              ! AEROSOL CLIM ALWAYS SEEN BY RADIATION
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              ! NOW SET THE AEROSOL DATA VALUES
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B1))
                END DO
              END DO           
              IF ( l_dust .AND. .NOT.l_twobin_dust ) THEN
                 ! MAKE ROOM FOR THE 6 BIN PROGNOSTIC DUST ASWELL
                 J_PLUS=J_PLUS+1
              ENDIF
           END IF
           IF (L_DUST .AND. (.NOT.l_twobin_dust) ) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_DUST) ) THEN
                 ! PROG. NOT USED TO EFFECT RAD FLUXES IF CLIM IS ON
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_DUST) THEN
                 ! DUST DIRECT EFFECT IS ON, SO USED IN RADIATION:
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =DUST_1(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
!        Two-bin dust, bin 1. prognostic only, 6 bin clim used in radn
         ELSE IF (i_aerosol == ip_twobindust_1 .AND.                    & 
                  l_dust.AND.l_twobin_dust ) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( l_use_arcl(ip_arcl_dust) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_DUST) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              END IF                 
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =DUST_1(LG, N_LAYER+1-I)
                END DO
              END DO
           
!
!        Six-bin dust, bin 2, clim or prognostic
         ELSE IF (i_aerosol == ip_dust_2 .AND.                          & 
                  (l_use_arcl(ip_arcl_dust).or.( l_dust .AND.           &
                   (.NOT.l_twobin_dust) ) ) ) THEN
           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B2))
                END DO
              END DO           
              IF ( l_dust .AND. .NOT.l_twobin_dust ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           END IF
           IF (L_DUST .AND. (.NOT.l_twobin_dust) ) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_DUST) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_DUST) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =DUST_2(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
!        Two-bin dust, bin 2. prognostic only, 6 bin clim used in radn
         ELSE IF (i_aerosol == ip_twobindust_2 .AND.                    & 
                  l_dust.AND.l_twobin_dust ) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( l_use_arcl(ip_arcl_dust) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_DUST) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              END IF                 
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =DUST_2(LG, N_LAYER+1-I)
                END DO
              END DO
!
!        Six-bin dust, bin 3, clim or prognostic
         ELSE IF (i_aerosol == ip_dust_3 .AND.                          & 
                  (l_use_arcl(ip_arcl_dust).or.( l_dust .AND.           &
                   (.NOT.l_twobin_dust) ) ) ) THEN
           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B3))
                END DO
              END DO           
              IF ( l_dust .AND. .NOT.l_twobin_dust ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           END IF
           IF (L_DUST .AND. (.NOT.l_twobin_dust) ) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_DUST) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_DUST) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =DUST_3(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF           
!
!        Six-bin dust, bin 4, clim or prognostic
         ELSE IF (i_aerosol == ip_dust_4 .AND.                          & 
                  (l_use_arcl(ip_arcl_dust).or.( l_dust .AND.           &
                   (.NOT.l_twobin_dust) ) ) ) THEN
           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B4))
                END DO
              END DO           
              IF ( l_dust .AND. .NOT.l_twobin_dust ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           END IF
           IF (L_DUST .AND. (.NOT.l_twobin_dust) ) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_DUST) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_DUST) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =DUST_4(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF           
!
!        Six-bin dust, bin 5, clim or prognostic
         ELSE IF (i_aerosol == ip_dust_5 .AND.                          & 
                  (l_use_arcl(ip_arcl_dust).or.( l_dust .AND.           &
                   (.NOT.l_twobin_dust) ) ) ) THEN
           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B5))
                END DO
              END DO           
              IF ( l_dust .AND. .NOT.l_twobin_dust ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           END IF
           IF (L_DUST .AND. (.NOT.l_twobin_dust) ) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_DUST) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_DUST) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =DUST_5(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF           
!
!        Six-bin dust, bin 6, clim or prognostic
         ELSE IF (i_aerosol == ip_dust_6 .AND.                          & 
                  (l_use_arcl(ip_arcl_dust).or.( l_dust .AND.           &
                   (.NOT.l_twobin_dust) ) ) ) THEN
           IF (L_USE_ARCL(IP_ARCL_DUST)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_DUST_B6))
                END DO
              END DO           
              IF ( l_dust .AND. .NOT.l_twobin_dust ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           END IF
           IF (L_DUST .AND. (.NOT.l_twobin_dust) ) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_DUST) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_DUST) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                     &
     &               =DUST_6(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Aerosols related to the sulphur cycle (note that dissolved
!        sulphate does not contribute to the direct effect):
         ELSE IF (i_aerosol == ip_accum_sulphate .AND.                  & 
                  (l_use_arcl(ip_arcl_sulp).or.l_sulpc_so2)) THEN
           IF (L_USE_ARCL(IP_ARCL_SULP)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_SULP_AC))
                END DO
              END DO
              IF ( L_SULPC_SO2 ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           IF (L_SULPC_SO2) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_SULP) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (l_use_sulpc_direct) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &            
     &               =ACCUM_SULPHATE(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (i_aerosol == ip_aitken_sulphate .AND.                 & 
                  (l_use_arcl(ip_arcl_sulp).or.l_sulpc_so2)) THEN
           
           IF (L_USE_ARCL(IP_ARCL_SULP)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_SULP_AK))
                END DO
              END DO
              IF ( L_SULPC_SO2 ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           IF (L_SULPC_SO2) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_SULP) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (l_use_sulpc_direct) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &            
     &               =AITKEN_SULPHATE(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Black Carbon/Soot aerosols
         ELSE IF (i_aerosol == ip_fresh_soot .AND.                      & 
                  (l_use_arcl(ip_arcl_blck).or.l_soot)) THEN 
           IF (L_USE_ARCL(IP_ARCL_BLCK)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_BLCK_FR))
                END DO
              END DO
              IF ( L_SOOT ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           IF (L_SOOT) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_BLCK) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_SOOT_DIRECT) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &            
     &               =FRESH_SOOT(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (i_aerosol == ip_aged_soot .AND.                       & 
                  (l_use_arcl(ip_arcl_blck).or.l_soot)) THEN 
           IF (L_USE_ARCL(IP_ARCL_BLCK)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_BLCK_AG))
                END DO
              END DO
              IF ( L_SOOT ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           IF (L_SOOT) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_BLCK) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (L_USE_SOOT_DIRECT) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &            
     &               =AGED_SOOT(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Biomass burning aerosols
         ELSE IF (i_aerosol == ip_biomass_1 .AND.                       & 
                  (l_use_arcl(ip_arcl_biom).or.l_biomass)) THEN 

           IF (L_USE_ARCL(IP_ARCL_BIOM)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_BIOM_FR))
                END DO
              END DO
              IF ( L_BIOMASS ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           IF (L_BIOMASS) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_BIOM) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (l_use_bmass_direct) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &            
     &               =FRESH_BMASS(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (i_aerosol == ip_biomass_2 .AND.                       & 
                  (l_use_arcl(ip_arcl_biom).or.l_biomass)) THEN 

           IF (L_USE_ARCL(IP_ARCL_BIOM)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &
     &               =ARCL(LG, N_LAYER+1-I,                             &
     &                     I_ARCL_COMPNTS(IP_ARCL_BIOM_AG))
                END DO
              END DO
              IF ( L_BIOMASS ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           IF (L_BIOMASS) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_BIOM) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (l_use_bmass_direct) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)                    &            
     &               =AGED_BMASS(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Sea salt aerosols
         ELSE IF (i_aerosol == ip_seasalt_film .AND.                    & 
                (l_use_arcl(ip_arcl_sslt).or.l_use_seasalt_direct)) THEN     

           IF (L_USE_ARCL(IP_ARCL_SSLT)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                   &
     &                    (ARCL(LG, N_LAYER+1-I,                        &
     &                     I_ARCL_COMPNTS(IP_ARCL_SSLT_FI))*4.0*PI      &
     &                   *DENSITY_SEA_SALT*(MODE_RADIUS_SS_FILM**3)   &
     &                   *EXP(4.5*(ALOG(SIGMA_SS_FILM))**2))          &
     &                   /(3.0*P(L, I)/(R*T(L, I)))
                END DO
              END DO
              IF ( L_USE_SEASALT_DIRECT ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           IF (L_USE_SEASALT_DIRECT) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_SSLT) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                   &
     &                    (SEA_SALT_FILM(LG, N_LAYER+1-I)*4.0*PI        &
     &                   *DENSITY_SEA_SALT*(MODE_RADIUS_SS_FILM**3)   &
     &                   *EXP(4.5*(ALOG(SIGMA_SS_FILM))**2))          &
     &                   /(3.0*P(L, I)/(R*T(L, I)))
                END DO
              END DO
           END IF
!
         ELSE IF (i_aerosol == ip_seasalt_jet .AND.                     & 
                  (l_use_arcl(ip_arcl_sslt).or.l_use_seasalt_direct)) THEN     

           IF (L_USE_ARCL(IP_ARCL_SSLT)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                   &
     &                    (ARCL(LG, N_LAYER+1-I,                        &
     &                     I_ARCL_COMPNTS(IP_ARCL_SSLT_JT))*4.0*PI      &
     &                   *DENSITY_SEA_SALT*(MODE_RADIUS_SS_JET**3)    &
     &                   *EXP(4.5*(ALOG(SIGMA_SS_JET))**2))           &
     &                   /(3.0*P(L, I)/(R*T(L, I)))
                END DO
              END DO
              IF ( L_USE_SEASALT_DIRECT ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           iF (L_USE_SEASALT_DIRECT) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_SSLT) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                   LG=I_GATHER(L)
                   AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                   &
     &                    (SEA_SALT_JET(LG, N_LAYER+1-I)*4.0*PI         &
     &                   *DENSITY_SEA_SALT*(MODE_RADIUS_SS_JET**3)    &
     &                   *EXP(4.5*(ALOG(SIGMA_SS_JET))**2))           &
     &                   /(3.0*P(L, I)/(R*T(L, I)))
                END DO
              END DO
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Organic Carbon from Fossil Fuel aerosol:
         ELSE IF (i_aerosol == ip_ocff_fresh .AND.                      & 
     &            (l_use_arcl(ip_arcl_ocff).or.l_ocff)) THEN 
         
           IF (L_USE_ARCL(IP_ARCL_OCFF)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                    &
                         ARCL(LG, N_LAYER+1-I,                          &
     &                     I_ARCL_COMPNTS(IP_ARCL_OCFF_FR))
                END DO
              END DO
              IF ( L_OCFF ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           IF (L_OCFF) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_OCFF) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (l_use_ocff_direct) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                     &
     &                   FRESH_OCFF(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!
         ELSE IF (i_aerosol == ip_ocff_aged .AND.                       & 
     &            (l_use_arcl(ip_arcl_ocff).or.l_ocff)) THEN 
         
           IF (L_USE_ARCL(IP_ARCL_OCFF)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                    &
                         ARCL(LG, N_LAYER+1-I,                          &
     &                     I_ARCL_COMPNTS(IP_ARCL_OCFF_AG))
                END DO
              END DO
              IF ( L_OCFF ) THEN
                 J_PLUS=J_PLUS+1
              ENDIF
           ENDIF
           IF (L_OCFF) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF ( L_USE_ARCL(IP_ARCL_OCFF) ) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ELSE IF (l_use_ocff_direct) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              ENDIF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                     &
     &                    AGED_OCFF(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Nitrate aerosol: only available as prognostic
         ELSE IF (I_AEROSOL == IP_NITRATE) THEN
           IF (L_NITRATE) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              IF (L_USE_NITRATE_DIRECT) THEN
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_RON
              ELSE
                 AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CLASSIC_ROFF
              END IF
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                    &
     &                    ACCUM_NITRATE(LG, N_LAYER+1-I)
                END DO
              END DO
           END IF         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Biogenic aerosols, only a climatology:
         ELSE IF ((I_AEROSOL == IP_BIOGENIC) .AND. L_USE_BIOGENIC) THEN
           AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
           AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
           DO I=I_TOP_COPY, N_LAYER
             DO L=1, N_PROFILE
                LG=I_GATHER(L)
                AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                      &
     &                   BIOGENIC(LG, N_LAYER+1-I)
             END DO
           END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        'Delta' aerosol: only ever a climatology
         ELSE IF (I_AEROSOL == IP_DELTA) THEN
           IF (L_USE_ARCL(IP_ARCL_DLTA)) THEN
              AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
              AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_ARCL_RON
              DO I=I_TOP_COPY, N_LAYER
                DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=                    &
                         ARCL(LG, N_LAYER+1-I,                          &
      &                    I_ARCL_COMPNTS(IP_ARCL_DLTA_DL))
                END DO
              END DO
           END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ELSE
!           The options to the radiation code do not require this
!           aerosol to be considered: its mixing ratio is set to 0.
!           If the spectral information has already removed all the 
!           redundany aerosols then this block of code should not 
!           normally be executed, but may be required for ease of 
!           including modifications.
!
            AEROSOL_MR_TYPE_INDEX(J+J_PLUS)=J
            AEROSOL_MR_SOURCE(J+J_PLUS) = IP_AERSRC_CUSACK_ROFF
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  AEROSOL_MIX_RATIO(L, I, J+J_PLUS)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
!
!        If using an extra top layer extrapolate the mixing ratio
!        from the adjacent layer, except in the case of Cusack 
!        climatological aerosols which are specifically set.
         IF (L_EXTRA_TOP .AND.                                          &
            AEROSOL_MR_SOURCE(J+J_PLUS) /= IP_AERSRC_CUSACK_RON .AND.   &
            AEROSOL_MR_SOURCE(J+J_PLUS) /= IP_AERSRC_CUSACK_ROFF ) THEN
           DO L=1, N_PROFILE
             AEROSOL_MIX_RATIO(L, 1, J+J_PLUS)                          &
               =AEROSOL_MIX_RATIO(L, 2, J+J_PLUS)
           ENDDO
         ENDIF
!
      ENDDO
      
      ! Check the number of aerosols actually set matches the number 
      ! expected and report diagnostic info if the numbers mismatch
      IF (N_AEROSOL_MR /= N_AEROSOL+J_PLUS) THEN
        ! writing details of the error to stdout:
        WRITE(IU_STDOUT, '(A)') '*** ERROR: '                           &
          // 'AN UNEXPECTED NUMBER OF AEROSOL SPECIES HAS BEEN SET.***' &
          // '(details below)'
        WRITE(IU_STDOUT, '(A,I3)') 'NUMBER OF AEROSOLS IN SPEC FILE:',  &
          N_AEROSOL
        WRITE(IU_STDOUT, '(A,I3)') 'NUMBER OF AEROSOLS EXPECTED (CLIMS' &
          // ' AND PROGNSOTICS):' ,  N_AEROSOL_MR 
        WRITE(IU_STDOUT, '(A,I3)') 'NUMBER OF AEROSOLS ACTUALLY SET:',  &
          N_AEROSOL+J_PLUS
        WRITE (IU_STDOUT,'(A)')                                         &
                  ' OUTPUT AEROSOL ARRAYS FROM R2_SET_AEROSOL_FIELD: '
        WRITE (IU_STDOUT,'(A)') ' (SEE rad_pcf.F90 TO INTERPRET THE'    &
          // ' AEROSOL TYPES AND SOURCES)'
        WRITE (IU_STDOUT,'(A)') ' J_MR, AEROSOL_MR_TYPE_INDEX, '        &
          // 'TYPE_AEROSOL(AEROSOL_MR_TYPE_INDEX), AEROSOL_MR_SOURCE'
        DO J_MR=1, N_AEROSOL_MR
           IF (AEROSOL_MR_TYPE_INDEX(J_MR) >= 1 .AND.                   &
               AEROSOL_MR_TYPE_INDEX(J_MR) <= N_AEROSOL) THEN
             WRITE (IU_STDOUT,'(I3,I3,I3,I3)') J_MR,                       &
               AEROSOL_MR_TYPE_INDEX(J_MR),                             &
               TYPE_AEROSOL(AEROSOL_MR_TYPE_INDEX(J_MR)),               &
               AEROSOL_MR_SOURCE(J_MR)
           ELSE
             WRITE (IU_STDOUT,'(I3,I3,A3,I3)') J_MR,                       &
               AEROSOL_MR_TYPE_INDEX(J_MR),' - ', AEROSOL_MR_SOURCE(J_MR)
           END IF
        ENDDO
        ! now writing out the error message itself
        WRITE(cmessage, '(A)') '*** ERROR: '                &
          // 'AN UNEXPECTED NUMBER OF AEROSOL SPECIES HAS BEEN SET.*** '&
          // ' (see stdout for details)'
        ierr=i_err_fatal 
        CALL ereport(RoutineName,ierr,cmessage)
        IF (lhook) CALL dr_hook('R2_SET_AEROSOL_FIELD',zhook_out,zhook_handle)
        RETURN
      ENDIF

!
! Add the MURK aerosol onto the Cusack climatology (via meso_frac), 
! but only in the boundary layer levels
      IF (L_MURK_RAD) THEN
        INV17=1./1.7
        DO J_MR=1, N_AEROSOL_MR
          DO I=(NLEVS-N_LEVELS_BL+1),NLEVS
            DO L=1, N_PROFILE
              LG=I_GATHER(L)
! Note that aerosol mass mixing ratios are in units of 1E9 kg/kg
              AEROSOL_MIX_RATIO(L,I,J_MR) = AEROSOL_MIX_RATIO(L,I,J_MR)   &
     &          + AERO_MESO(LG,NLEVS+1-I)*1E-9*                           &
     &            MESO_FRAC(AEROSOL_MR_TYPE_INDEX(J_MR))*INV17
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
!

      IF (lhook) CALL dr_hook('R2_SET_AEROSOL_FIELD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE R2_SET_AEROSOL_FIELD
