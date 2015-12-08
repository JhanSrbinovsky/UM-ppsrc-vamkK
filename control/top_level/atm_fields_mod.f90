! DECLARE_ATM_FIELDS_MOD
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module atm_fields_mod

  implicit none

!     TRACERS are a special case
      ! dummy array used to prevent crashing when there are no tracers
      real, target  :: dummy_field(1) = (/-1.0/) 
      real, pointer :: TRACER(:)
      real, pointer :: TRACER_UKCA(:)
      real, pointer :: TRACER_UKCA_LBC(:)
      real, pointer :: TRACER_UKCA_LBC_TEND(:)
      real, pointer :: TRACER_LBC(:)
      real, pointer :: TRACER_LBC_TEND(:)
!     Add a specific pointer for ozone tracer and cariolle parameters
      real, pointer :: OZONE_TRACER(:)
      real, pointer :: O3_PROD_LOSS(:)
      real, pointer :: O3_P_L_VMR(:)
      real, pointer :: O3_VMR(:)
      real, pointer :: O3_P_L_TEMP(:)
      real, pointer :: O3_TEMP(:)
      real, pointer :: O3_P_L_COLO3(:)
      real, pointer :: O3_COLO3(:)
      ! 1: Array Variables (dimensions are resolution dependent.)

! Oxidant concentrations from UKCA for use in HadGEM sulphur cycle
      real, pointer :: OH_UKCA(:)
      real, pointer :: H2O2_UKCA(:)
      real, pointer :: HO2_UKCA(:)
      real, pointer :: O3_UKCA(:)
      real, pointer :: HNO3_UKCA(:)

      ! 1.1: Data variables stored in primary space.

      REAL, POINTER :: DryRho (:)
      REAL, POINTER :: EtaDot(:)
      REAL, POINTER :: ThetaV(:)
      REAL, POINTER :: psi_w_surf(:)
      REAL, POINTER :: psi_w_lid(:)
      REAL, POINTER :: m_v(:)
      REAL, POINTER :: m_cl(:)
      REAL, POINTER :: m_cf(:)
      REAL, POINTER :: m_cf2(:)
      REAL, POINTER :: m_gr(:)
      REAL, POINTER :: m_r(:)
      REAL, POINTER :: exner_surf(:)
      REAL, POINTER :: exner(:)

      real, pointer :: U(:)      ! u component of wind
      real, pointer :: V(:)      ! v component of wind
      real, pointer :: W(:)      ! w component of wind
      real, pointer :: RHO(:)    ! Density
      real, pointer :: THETA(:)  ! Potential temperature
      real, pointer :: Q(:)        ! Specific humidity
      real, pointer :: QCL(:)      ! qcl
      real, pointer :: QCF(:)      ! qcf
      real, pointer :: QCF2(:)     ! second ice
      real, pointer :: QRAIN(:)    ! rain
      real, pointer :: QGRAUP(:)   ! graupel

      ! TKE based turbulence scheme
      real, pointer :: E_TRB(:)     ! TKE
      real, pointer :: TSQ_TRB(:)   ! Self covariance of thetal'
      real, pointer :: QSQ_TRB(:)   ! Self coveriance of qw'
      real, pointer :: COV_TRB(:)   ! Correlation between thetal' and qw'
      real, pointer :: ZHPAR_SHCU(:)! Height of mixed layer used to
                                    ! evaluate the non-grad buoy flux

      ! Exner pressure on rho levels
      real, pointer :: EXNER_RHO_LEVELS(:)

      real, pointer :: U_ADV(:) ! Advective u component of wind
      real, pointer :: V_ADV(:) ! Advective v component of wind
      real, pointer :: W_ADV(:) ! Advective w component of wind

      ! 1.2: Data variables stored in secondary space.
      real, pointer :: P(:)                  ! Pressure on rho le

      ! Pressure on theta levels
      real, pointer :: P_THETA_LEVELS(:)

      ! Exner pressure on theta levels
      real, pointer :: EXNER_THETA_LEVELS(:)

      ! 1.3: Cloud Fields
      real, pointer :: CCA(:)                ! Convective cloud amoun
      real, pointer :: CF_AREA(:)            ! Area Cloud Fraction
      real, pointer :: CF_BULK(:)            ! Bulk Cloud Fraction
      real, pointer :: CF_LIQUID(:)          ! Liquid cloud fraction
      real, pointer :: CF_FROZEN(:)          ! Frozen cloud fraction

      ! 1.4: Soil Ancillary fields
      real, pointer :: DEEP_SOIL_TEMP(:)   ! Deep soil temperature

      real, pointer :: SMCL(:)   ! Soil moisture content in layers
      real, pointer :: STHU(:)   ! Unfrozen soil moisture fraction
      real, pointer :: STHF(:)   ! Frozen soil moisture fraction

      ! 1.5: Radiation Increments
      real, pointer :: SW_INCS(:)    ! SW radiation increments
      real, pointer :: LW_INCS(:)    ! LW radiation increments
! PAR radiation increment
      real, pointer :: DIRPAR(:)

      ! 1.6: Ozone
      real, pointer :: O3(:)          ! Ozone
!  tropopause-based ozone
      real, pointer :: TPPSOZONE(:)
      ! 1.7: Tracer and aerosol fields
      real, pointer :: MURK_SOURCE(:)    ! multi-level murk source
      real, pointer :: MURK(:)           ! multi-level murk concent
      real, pointer :: DUST_DIV1(:)      ! dust mmr, division 1
      real, pointer :: DUST_DIV2(:)      ! dust mmr, division 2
      real, pointer :: DUST_DIV3(:)      ! dust mmr, division 3
      real, pointer :: DUST_DIV4(:)      ! dust mmr, division 4
      real, pointer :: DUST_DIV5(:)      ! dust mmr, division 5
      real, pointer :: DUST_DIV6(:)      ! dust mmr, division 6
      real, pointer :: SO2(:)            ! sulphur dioxide gas
      real, pointer :: DMS(:)            ! dimethyl sulphide gas
      real, pointer :: SO4_AITKEN(:)     ! Aitken mode sulphate aer
      real, pointer :: SO4_ACCU(:)       ! accumulation mode sulpha
      real, pointer :: SO4_DISS(:)       ! dissloved  sulphate aero
      real, pointer :: H2O2(:)           ! hydrogen peroxide mmr
      real, pointer :: NH3(:)            ! ammonia gas mmr
      real, pointer :: SOOT_NEW(:)       ! fresh soot mmr
      real, pointer :: SOOT_AGD(:)       ! aged soot mmr
      real, pointer :: SOOT_CLD(:)       ! soot in cloud mmr
      real, pointer :: BMASS_NEW(:)      ! fresh biomass mmr
      real, pointer :: BMASS_AGD(:)      ! aged biomass mmr
      real, pointer :: BMASS_CLD(:)      ! cloud biomass mmr
      real, pointer :: OCFF_NEW(:)       ! fresh OCFF mmr
      real, pointer :: OCFF_AGD(:)       ! aged OCFF mmr
      real, pointer :: OCFF_CLD(:)       ! OCFF in cloud mmr
      real, pointer :: SO2_NATEM(:)      ! natural SO2 emissions
      real, pointer :: OH(:)             ! hydroxyl radical ancilla
      real, pointer :: HO2(:)            ! hydrogen dioxide ancilla
      real, pointer :: H2O2_LIMIT(:)     ! limiting H2O2 ancillary
      real, pointer :: O3_CHEM(:)        ! ozone for chemistry anci
      real, pointer :: NITR_ACC(:)       ! accumulation ammonium nitrate
      real, pointer :: NITR_DISS(:)      ! dissolved ammonium nitrate
      real, pointer :: CO2(:)            ! 3D CO2 FIELD
! 1.8: Multi-level user ancillary fields
      real, pointer :: USER_MULT1(:)     ! multi-level user ancilla
      real, pointer :: USER_MULT2(:)     ! multi-level user ancilla
      real, pointer :: USER_MULT3(:)     ! multi-level user ancilla
      real, pointer :: USER_MULT4(:)     ! multi-level user ancilla
      real, pointer :: USER_MULT5(:)     ! multi-level user ancilla
      real, pointer :: USER_MULT6(:)     ! multi-level user ancilla
      real, pointer :: USER_MULT7(:)     ! multi-level user ancilla
      real, pointer :: USER_MULT8(:)     ! multi-level user ancilla
      real, pointer :: USER_MULT9(:)     ! multi-level user ancilla
      real, pointer :: USER_MULT10(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT11(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT12(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT13(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT14(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT15(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT16(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT17(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT18(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT19(:)    ! multi-level user ancilla
      real, pointer :: USER_MULT20(:)    ! multi-level user ancilla

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
      real, pointer :: OROG_LBC(:)                       ! Orography LBC
      real, pointer :: U_LBC(:)                          ! U LBC
      real, pointer :: V_LBC(:)                          ! V LBC
      real, pointer :: W_LBC(:)                          ! W LBC
      real, pointer :: RHO_LBC(:)                        ! RHO LBC
      real, pointer :: THETA_LBC(:)                      ! Theta LBC
      real, pointer :: Q_LBC(:)                          ! Q LBC
      real, pointer :: QCL_LBC(:)                        ! QCL LBC
      real, pointer :: QCF_LBC(:)                        ! QCF LBC
      real, pointer :: QCF2_LBC(:)                       ! 2nd Ice LBC
      real, pointer :: QRAIN_LBC(:)                      ! Rain LBC
      real, pointer :: QGRAUP_LBC(:)                     ! Graupel LBC
      real, pointer :: CF_BULK_LBC(:)                    ! CF_BULK LBC
      real, pointer :: CF_LIQUID_LBC(:)                  ! CF_LIQUID_LBC
      real, pointer :: CF_FROZEN_LBC(:)                  ! CF_FROZEN_LBC
      real, pointer :: EXNER_LBC(:)                      ! Exner LBC
      real, pointer :: U_ADV_LBC(:)                      ! U_ADV LBC
      real, pointer :: V_ADV_LBC(:)                      ! V_ADV LBC
      real, pointer :: W_ADV_LBC(:)                      ! W_ADV LBC
      real, pointer :: MURK_LBC(:)                       ! Murk aerosol LBC
      real, pointer :: DUST_DIV1_LBC(:)      ! dust mmr, division 1 LBC
      real, pointer :: DUST_DIV2_LBC(:)      ! dust mmr, division 2 LBC
      real, pointer :: DUST_DIV3_LBC(:)      ! dust mmr, division 3 LBC
      real, pointer :: DUST_DIV4_LBC(:)      ! dust mmr, division 4 LBC
      real, pointer :: DUST_DIV5_LBC(:)      ! dust mmr, division 5 LBC
      real, pointer :: DUST_DIV6_LBC(:)      ! dust mmr, division 6 LBC
      real, pointer :: SO2_LBC(:)            ! sulphur dioxide gas LBC
      real, pointer :: DMS_LBC(:)            ! dimethyl sulphide gas LBC
      real, pointer :: SO4_AITKEN_LBC(:)     ! Aitken mode sulphate aerosol LBC
      real, pointer :: SO4_ACCU_LBC(:)       ! accumulation mode sulphate LBC
      real, pointer :: SO4_DISS_LBC(:)       ! dissolved sulphate aero LBC
      real, pointer :: NH3_LBC(:)            ! ammonia gas LBC
      real, pointer :: SOOT_NEW_LBC(:)       ! fresh soot LBC
      real, pointer :: SOOT_AGD_LBC(:)       ! aged soot LBC
      real, pointer :: SOOT_CLD_LBC(:)       ! soot in cloud LBC
      real, pointer :: BMASS_NEW_LBC(:)      ! fresh biomass LBC
      real, pointer :: BMASS_AGD_LBC(:)      ! aged biomass LBC
      real, pointer :: BMASS_CLD_LBC(:)      ! cloud biomass LBC
      real, pointer :: OCFF_NEW_LBC(:)       ! fresh fossil fuel LBC
      real, pointer :: OCFF_AGD_LBC(:)       ! aged fossil fuel LBC
      real, pointer :: OCFF_CLD_LBC(:)       ! cloud fossil fuel LBC
      real, pointer :: NITR_ACC_LBC(:)       ! accumulation ammonium nitrate LBC
      real, pointer :: NITR_DISS_LBC(:)      ! dissolved ammonium nitrate LBC

      ! Lateral Boundary Condition tendencies
      real, pointer :: U_LBC_TEND(:)                     ! U LBC  tendencies
      real, pointer :: V_LBC_TEND(:)                     ! V LBC tendencies
      real, pointer :: W_LBC_TEND(:)                     ! W LBC tendencies
      real, pointer :: RHO_LBC_TEND(:)                   ! RHO LBC tendencies
      real, pointer :: THETA_LBC_TEND(:)                 ! Theta LBC tendencies
      real, pointer :: Q_LBC_TEND(:)                     ! Q LBC tendencies
      real, pointer :: QCL_LBC_TEND(:)                   ! QCL LBC tendencies
      real, pointer :: QCF_LBC_TEND(:)                   ! QCF LBC tendencies
      real, pointer :: QCF2_LBC_TEND(:)                  ! 2nd Ice
      real, pointer :: QRAIN_LBC_TEND(:)                 ! Rain LBC tendencies
      real, pointer :: QGRAUP_LBC_TEND(:)                ! Graupel
      real, pointer :: CF_BULK_LBC_TEND(:)               ! CF_BULK LBC tend'cies
      real, pointer :: CF_LIQUID_LBC_TEND(:)             ! CF_LIQUID_LBC t'cies
      real, pointer :: CF_FROZEN_LBC_TEND(:)             ! CF_FROZEN_LBC t'cies
      real, pointer :: EXNER_LBC_TEND(:)                 ! Exner LBC tendencies
      real, pointer :: U_ADV_LBC_TEND(:)                 ! U_ADV LBC tendencies
      real, pointer :: V_ADV_LBC_TEND(:)                 ! V_ADV LBC tendencies
      real, pointer :: W_ADV_LBC_TEND(:)                 ! W_ADV LBCtendencies
      real, pointer :: MURK_LBC_TEND(:)                  ! Murk aerosol LBC tend
      real, pointer :: DUST_DIV1_LBC_TEND(:)      ! dust mmr, division 1 LBC tend
      real, pointer :: DUST_DIV2_LBC_TEND(:)      ! dust mmr, division 2 LBC tend
      real, pointer :: DUST_DIV3_LBC_TEND(:)      ! dust mmr, division 3 LBC tend
      real, pointer :: DUST_DIV4_LBC_TEND(:)      ! dust mmr, division 4 LBC tend
      real, pointer :: DUST_DIV5_LBC_TEND(:)      ! dust mmr, division 5 LBC tend
      real, pointer :: DUST_DIV6_LBC_TEND(:)      ! dust mmr, division 6 LBC tend
      real, pointer :: SO2_LBC_TEND(:)            ! sulphur dioxide gas LBC tend
      real, pointer :: DMS_LBC_TEND(:)            ! dimethyl sulphide gas LBC tend
      real, pointer :: SO4_AITKEN_LBC_TEND(:)     ! Aitken mode sulphate aerosol LBC
      real, pointer :: SO4_ACCU_LBC_TEND(:)       ! accumulation mode sulphate LBC
      real, pointer :: SO4_DISS_LBC_TEND(:)       ! dissolved sulphate aero LBC
      real, pointer :: NH3_LBC_TEND(:)            ! ammonia gas LBC tend
      real, pointer :: SOOT_NEW_LBC_TEND(:)       ! fresh soot LBC tend
      real, pointer :: SOOT_AGD_LBC_TEND(:)       ! aged soot LBC tend
      real, pointer :: SOOT_CLD_LBC_TEND(:)       ! soot in cloud LBC tend
      real, pointer :: BMASS_NEW_LBC_TEND(:)      ! fresh biomass LBC tend
      real, pointer :: BMASS_AGD_LBC_TEND(:)      ! aged biomass LBC tend
      real, pointer :: BMASS_CLD_LBC_TEND(:)      ! cloud biomass LBC tend
      real, pointer :: OCFF_NEW_LBC_TEND(:)       ! fresh fossil fuel LBC tend
      real, pointer :: OCFF_AGD_LBC_TEND(:)       ! aged fossil fuel LBC tend
      real, pointer :: OCFF_CLD_LBC_TEND(:)       ! cloud fossil fuel LBC tend
      real, pointer :: NITR_ACC_LBC_TEND(:)       ! accumulation ammonium nitrate LBC tend
      real, pointer :: NITR_DISS_LBC_TEND(:)      ! dissolved ammonium nitrate LBC tend


      ! 2: Scalar Variables

      ! 2.1: Data variables stored in primary space.
      real, pointer :: TSTAR(:)          ! Surface temperature
      real, pointer :: LAND(:)           ! Land sea mask
      real, pointer :: TSTAR_ANOM(:)     ! Surface temperature anomaly
!   2.15: Fields for coastal tiling
      real, pointer :: FRAC_LAND(:)      ! Land fraction in grid box
      real, pointer :: TSTAR_LAND(:)     ! Land surface temperature
      real, pointer :: TSTAR_SEA(:)      ! Sea surface temperature
      real, pointer :: TSTAR_SICE(:)     ! Sea-ice surface temperature
      real, pointer :: TSTAR_SICE_CAT(:) ! Sea-ice category surface temp
! Set pointers for sea-ice and land albedos
      real, pointer :: SICE_ALB(:)       ! Sea-ice albedo
      real, pointer :: LAND_ALB(:)       ! Mean land albedo

      ! 2.2: Data variables stored in secondary space.

      real, pointer :: PSTAR(:)          ! Surface pressure

      ! 2.3: Cloud fields

      real, pointer :: CCB(:)            ! Convective cloud base
      real, pointer :: CCT(:)            ! Convective cloud top

      real, pointer :: CCLWP(:)          ! Convective cloud liquid water path
      real, pointer :: deep_flag(:)      ! Deep convection flag
      real, pointer :: past_precip(:)    ! Past convective precipitation
      real, pointer :: past_conv_ht(:)   ! Past convective height

      ! 2.4: Boundary layer fields

      real, pointer :: ZH(:)             ! Boundary layer depth    
 
      REAL, POINTER :: ddmfx(:)          ! Convective dowdraught
                                         ! mass-flux at cloud base

      ! Standard deviation of turbulent fluctuations of layer 1
                                         ! temperature
      real, pointer :: T1_SD(:)

      ! Standard deviation of turbulent fluctuations of layer 1 humidity
      real, pointer :: Q1_SD(:)

      ! Decoupled screen-level temperature
      REAL, POINTER :: TScrnDcl_SSI(:)
      REAL, POINTER :: TScrnDcl_TILE(:)
      REAL, POINTER :: tStbTrans(:)

      ! Number of model levels in the  turbulently mixed layer
      real, pointer :: NTML(:)

      ! Top level for turb mixing in any decoupled Sc layer
      real, pointer :: NTDSC(:)

      ! Bottom level for turb mixing in any decoupled Sc layer
      real, pointer :: NBDSC(:)

      real, pointer :: CUMULUS(:)      ! Boundary layer convection flag

      ! 2.4: Soil Ancillary fields

      real, pointer :: SAT_SOILW_SUCTION(:) ! Saturated soil water suction
      real, pointer :: THERM_CAP(:)     ! Thermal capacity
      real, pointer :: THERM_COND(:)    ! Thermal conductivity
      real, pointer :: VOL_SMC_CRIT(:)  ! Vol smc at critical point
      real, pointer :: VOL_SMC_WILT(:)  ! Vol smc at wilting point
      real, pointer :: VOL_SMC_SAT(:)   ! Vol smc at saturation
      real, pointer :: SAT_SOIL_COND(:) ! Saturated soil conductivity
      real, pointer :: CLAPP_HORN(:)    ! Clapp-Hornberger B coefficient

      ! 2.5: Other surface fields
      REAL, POINTER :: canopy_water(:)  ! Canopy Water
      REAL, POINTER :: z0(:)            ! Roughness length, used for sea points
      REAL, POINTER :: gs(:)            ! Gridbox mean canopy conductance

      ! 2.6: Orographic Ancillary fields

      ! (changed from "OROG" due to naming conflicts)
      real, pointer :: OROGRAPHY(:)     ! Orographic height
      real, pointer :: OROG_SD(:)       ! Standard Deviation of orography
      real, pointer :: OROG_SIL(:)      ! Silhouette area of orography
      real, pointer :: OROG_HO2(:)      ! Peak to trough height/(2*sqrt2)
      real, pointer :: OROG_GRAD_X(:)  
      real, pointer :: OROG_GRAD_Y(:)  
      real, pointer :: OROG_UNFILT(:)  
      real, pointer :: OROG_GRAD_XX(:)  ! Orographic gradient xx
      real, pointer :: OROG_GRAD_XY(:)  ! Orographic gradient xy
      real, pointer :: OROG_GRAD_YY(:)  ! Orographic gradient yy

      ! 2.7: Sea/Sea Ice

      real, pointer :: U_SEA(:)         ! Surface current (u component)
      real, pointer :: V_SEA(:)         ! Surface current (v component)
      real, pointer :: U_0_P(:)         ! Surace  current (u) on p-grid
      real, pointer :: V_0_P(:)         ! Surface current (v) on p-grid
      real, pointer :: ICE_FRACTION(:)  ! Sea ice fraction
      real, pointer :: ICE_THICKNESS(:) ! Sea ice depth
      real, pointer :: TI(:)            ! Sea ice temperature
      real, pointer :: ICE_FRACT_CAT(:) ! Sea ice fraction on categories
      real, pointer :: ICE_THICK_CAT(:) ! Sea ice thickness on categories
      real, pointer :: TI_CAT(:)        ! Sea ice temperature on categories
      real, pointer :: ICE_K_CAT(:)     ! Sea ice effect cond on categories

      ! 2.8: Snow

      real, pointer :: SNODEP(:)        ! Snow depth on land
      real, pointer :: SNODEP_SEA(:)    ! Snow depth on sea ice
      real, pointer :: SNODEP_SEA_CAT(:) ! Snow depth on sea ice catagories
      real, pointer :: CATCH_SNOW(:)    ! Coniferous canopy
!                               ! snow capacity
      real, pointer :: SNOW_GRND(:)     ! Snow below canopy
      real, pointer :: SNSOOT(:)        ! Snow soot content

! 2.9: aerosol emission fields, including mineral dust parent soil props

      real, pointer :: SOIL_CLAY(:)                    ! soil clay fraction
      real, pointer :: SOIL_SILT(:)                    ! soil silt fraction
      real, pointer :: SOIL_SAND(:)                    ! soil sand fraction
      real, pointer :: DUST_MREL1(:)                   ! soil rel mass, div 1
      real, pointer :: DUST_MREL2(:)                   ! soil rel mass, div 2
      real, pointer :: DUST_MREL3(:)                   ! soil rel mass, div 3
      real, pointer :: DUST_MREL4(:)                   ! soil rel mass, div 4
      real, pointer :: DUST_MREL5(:)                   ! soil rel mass, div 5
      real, pointer :: DUST_MREL6(:)                   ! soil rel mass, div 6


      real, pointer :: SO2_EM(:)        ! sulphur dioxide emission
      real, pointer :: DMS_EM(:)        ! dimethyl sulphide emission
      real, pointer :: SO2_HILEM(:)     ! high level SO2 emissions
      real, pointer :: NH3_EM(:)        ! ammonia gas surface emiss
      real, pointer :: SOOT_EM(:)       ! fresh soot surface emissions
      real, pointer :: SOOT_HILEM(:)    ! fresh soot high lev emissions
      real, pointer :: BMASS_EM(:)      ! fresh bmass surface emissions
      real, pointer :: BMASS_HILEM(:)   ! fresh bmass high lev emissions
      real, pointer :: OCFF_EM(:)       ! fresh OCFF surface emissions
      real, pointer :: OCFF_HILEM(:)    ! fresh OCFF high lev emissions
      real, pointer :: DMS_CONC(:)      ! seawater dimethyl sulphide conc.
      real, pointer :: DMS_OFLUX(:)     ! DMS flux from ocean model

      ! 2.10: User ancillary fields
      real, pointer :: USER_ANC1(:)         ! user ancillary field 1
      real, pointer :: USER_ANC2(:)         ! user ancillary field 2
      real, pointer :: USER_ANC3(:)         ! user ancillary field 3
      real, pointer :: USER_ANC4(:)         ! user ancillary field 4
      real, pointer :: USER_ANC5(:)         ! user ancillary field 5
      real, pointer :: USER_ANC6(:)         ! user ancillary field 6
      real, pointer :: USER_ANC7(:)         ! user ancillary field 7
      real, pointer :: USER_ANC8(:)         ! user ancillary field 8
      real, pointer :: USER_ANC9(:)         ! user ancillary field 9
      real, pointer :: USER_ANC10(:)        ! user ancillary field 10
      real, pointer :: USER_ANC11(:)        ! user ancillary field 11
      real, pointer :: USER_ANC12(:)        ! user ancillary field 12
      real, pointer :: USER_ANC13(:)        ! user ancillary field 13
      real, pointer :: USER_ANC14(:)        ! user ancillary field 14
      real, pointer :: USER_ANC15(:)        ! user ancillary field 15
      real, pointer :: USER_ANC16(:)        ! user ancillary field 16
      real, pointer :: USER_ANC17(:)        ! user ancillary field 17
      real, pointer :: USER_ANC18(:)        ! user ancillary field 18
      real, pointer :: USER_ANC19(:)        ! user ancillary field 19
      real, pointer :: USER_ANC20(:)        ! user ancillary field 20

      !   2.11: Store arrays for energy correction calculation
      real, pointer :: NET_FLUX(:)                   ! Net energy flux
      real, pointer :: NET_MFLUX(:)                  ! Net moisture flux

      !   2.12: Tiled Vegetation and Triffid fields
      real, pointer :: FRAC_TYP(:)        ! Fractions of surface type
      real, pointer :: FRAC_CON1(:)       ! Fractions of surface type
      real, pointer :: FRAC_CON2(:)       ! Fractions of surface type
      real, pointer :: FRAC_CON3(:)       ! Fractions of surface type
      real, pointer :: FRAC_CON4(:)       ! Fractions of surface type
      real, pointer :: FRAC_CON5(:)       ! Fractions of surface type
      real, pointer :: FRAC_CON6(:)       ! Fractions of surface type
      real, pointer :: FRAC_CON7(:)       ! Fractions of surface type
      real, pointer :: FRAC_CON8(:)       ! Fractions of surface type
      real, pointer :: FRAC_CON9(:)       ! Fractions of surface type
      real, pointer :: LAI_PFT(:)         ! LAI of plant functional types
      real, pointer :: CANHT_PFT(:)       ! Canopy hght of plant func types
      ! (changed from "DISTURB" due to naming conflicts)
      real, pointer :: DISTURB_VEG(:)     ! Disturbed fraction of vegetation 
      real, pointer :: SOIL_ALB(:)        ! Snow-free albedo of bare soil
      REAL, POINTER :: obs_alb_sw(:)      ! Observed snow-free SW albedo 
      REAL, POINTER :: obs_alb_vis(:)     ! Observed snow-free VIS albedo 
      REAL, POINTER :: obs_alb_nir(:)     ! Observed snow-free NIR albedo 
      real, pointer :: SOIL_CARB(:)       ! Soil carbon content
      real, pointer :: SOIL_CARB1(:)      ! Soil carbon content DPM
      real, pointer :: SOIL_CARB2(:)      ! Soil carbon content RPM
      real, pointer :: SOIL_CARB3(:)      ! Soil carbon content BIO
      real, pointer :: SOIL_CARB4(:)      ! Soil carbon content HUM
      real, pointer :: NPP_PFT_ACC(:)     ! Accumulated NPP on PFTs
      real, pointer :: G_LF_PFT_ACC(:)    ! Accum. leaf turnover rate PFTs
      real, pointer :: G_PHLF_PFT_ACC(:)  ! Accumulated phenological leaf
                                    ! turnover rate on PFTs
      real, pointer :: RSP_W_PFT_ACC(:)   ! Accum. wood respiration on PFTs
      real, pointer :: RSP_S_ACC(:)       ! Accumulated soil respiration
      real, pointer :: RSP_S_ACC1(:)      ! Accumulated soil respiration DPM
      real, pointer :: RSP_S_ACC2(:)      ! Accumulated soil respiration RPM
      real, pointer :: RSP_S_ACC3(:)      ! Accumulated soil respiration BIO
      real, pointer :: RSP_S_ACC4(:)      ! Accumulated soil respiration HUM
      real, pointer :: CAN_WATER_TILE(:)  ! Canopy water content on tiles
      real, pointer :: CATCH_TILE(:)      ! Canopy capacity on tiles
      real, pointer :: INFIL_TILE(:)      ! Max infiltration rate on tiles
      real, pointer :: RGRAIN_TILE(:)     ! Snow grain size on tiles
      real, pointer :: SNODEP_TILE(:)     ! Snow depth on tiles
      real, pointer :: TSTAR_TILE(:)      ! Surface temperature on tiles
      real, pointer :: Z0_TILE(:)         ! Surface roughness on tiles
      REAL, POINTER :: z0h_tile(:)        ! Snow-free thermal roughness 
                                          ! on tiles
      REAL, POINTER :: dolr_field(:)      ! TOA - surface upward LW at
                                  ! radiation timestep
      real, pointer :: LW_DOWN(:)         ! Surface downward LW at
                                  ! radiation timestep
      ! (changed from SW_TILE b/c of naming conflicts)
      real, pointer :: SW_TILE_RTS(:)     ! Surface net SW on land tiles at
                                  ! radiation timestep

! MORUSES - new urban two-tile scheme
      real, pointer :: hgt(:)      ! Building height (m)
      real, pointer :: hwr(:)      ! Height to width
      real, pointer :: wrr(:)      ! Width ratio
      real, pointer :: disp(:)     ! Displacement height (m)
      real, pointer :: ztm(:)      ! Eff roughness length of momentum (m)
      real, pointer :: albwl(:)    ! Wall albedo
      real, pointer :: albrd(:)    ! Road albedo
      real, pointer :: emisw(:)    ! Wall emmissivity
      real, pointer :: emisr(:)    ! Road emmissivity

!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!      real, pointer :: TSLAB(:)     ! Temperature of slab ocean.
!      real, pointer :: TCLIM(:)     ! Climatological SST's
!      real, pointer :: HCLIM(:)     ! Climatological ice depth
!      real, pointer :: CHEAT(:)     ! Caryheat (from ice model to ocn)
!      real, pointer :: OIFLX(:)     ! Ocean-Ice heat flux
!      real, pointer :: UICE(:)      ! Zonal comp of ice velocity
!      real, pointer :: VICE(:)      ! Meridional comp of ice velocity
!      real, pointer :: SIG11NE(:)   !  Internal stresses for
!      real, pointer :: SIG11SE(:)   !  EVP ice dynamics
!      real, pointer :: SIG11SW(:)
!      real, pointer :: SIG11NW(:)
!      real, pointer :: SIG12NE(:)
!      real, pointer :: SIG12SE(:)
!      real, pointer :: SIG12SW(:)
!      real, pointer :: SIG12NW(:)
!      real, pointer :: SIG22NE(:)
!      real, pointer :: SIG22SE(:)
!      real, pointer :: SIG22SW(:)
!      real, pointer :: SIG22NW(:)
!! END MODIFIED BY AT

!   2.14: Carbon cycle fields
      real, pointer :: CO2FLUX(:)   ! Ocean CO2 flux (Kg CO2/m2/s1)
      real, pointer :: CO2_EMITS(:) ! Surface CO2 emissions (Kg CO2/m2/s1)

!   2.15: Fields carried forward from previous version.
!         May not be required
!!      real, pointer :: SURF_RESIST_NIT(:)  ! Surface resistance on
!!                                    ! non-ice tiles
!!      real, pointer :: ROOT_DEPTH_NIT(:)   ! Root depth on non-ice tiles
!!      real, pointer :: Z0V_TYP(:)          ! Vegetative Roughness length on
!!                                    ! tiles
!!      real, pointer :: TSNOW(:)            ! Snow surface layer temperature
!!      real, pointer :: ICE_EDGE(:)
!!      real, pointer :: OROG_TENDENCY(:)    ! Orographic tendencies
!!      real, pointer :: OROG_SD_TENDENCY(:) ! Orographic variable tendency

      ! Pointers for ATMOSPHERE model constants. Scalars only.
      ! Addresses in level dependent constants array.
!!      real, pointer :: ETATHETA(:)
!!      real, pointer :: ETARHO(:)
!!      real, pointer :: RHCRIT(:)
!!      real, pointer :: SOIL_THICKNESS(:)
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      real, pointer :: zseak_theta(:) ! zsea(k) on theta levels
      real, pointer :: Ck_theta(:)    ! C(k)    on theta levels
      real, pointer :: zseak_rho(:)   ! zsea(k) on rho levels
      real, pointer :: Ck_rho(:)      ! C(k)    on rho levels

!   2.16: Fields for large-scale hydrology scheme.

      real, pointer :: TI_MEAN(:)          !Mean topographic index
      real, pointer :: TI_SIG(:)           !Standard dev. in topographic index
      real, pointer :: FEXP(:)             !Exponential decay in soil
!                                  ! saturated conductivity
      real, pointer :: GAMMA_INT(:)        !Integrated gamma function
      real, pointer :: WATER_TABLE(:)      !Water table depth
      real, pointer :: FSFC_SAT(:)         !Surface saturation fraction
      real, pointer :: F_WETLAND(:)        !Wetland fraction

      real, pointer :: STHZW(:)
      real, pointer :: A_FSAT(:)
      real, pointer :: C_FSAT(:)
      real, pointer :: A_FWET(:)
      real, pointer :: C_FWET(:)

!   2.17: Fields for River routing.
      real, pointer :: RIV_SEQUENCE(:)   ! River sequence
      real, pointer :: RIV_DIRECTION(:)  ! River direction
      real, pointer :: RIV_STORAGE(:)    ! River water storage
      real, pointer :: TOT_SURFROFF(:)   ! Accumulated surface runoff
      real, pointer :: TOT_SUBROFF(:)    !     "       sub-surface runoff
      real, pointer :: RIV_INLANDATM(:)       ! inland basin outflow
 
! Fields for water conservation correction due to lake evaporation: 
      real, pointer :: ACC_LAKE_EVAP(:)  ! Acc lake evaporation 
 
! Fields for grid-to-grid river routing (river routing 2A)
      real, pointer :: RIV_IAREA(:)      ! Drainage area
      real, pointer :: RIV_SLOPE(:)      ! Grid-cell slope
      real, pointer :: RIV_FLOWOBS1(:)   ! Initial values of flow
      real, pointer :: RIV_INEXT(:)      ! Flow direction (x)
      real, pointer :: RIV_JNEXT(:)      ! Flow direction (y)
      real, pointer :: RIV_LAND(:)       ! Land-type (land/river/sea)
      real, pointer :: RIV_SUBSTORE(:)   ! Subsurface storage
      real, pointer :: RIV_SURFSTORE(:)  ! Surface storage
      real, pointer :: RIV_FLOWIN(:)     ! Surface inflow
      real, pointer :: RIV_BFLOWIN(:)    ! Subsurface inflow

! Fields used when coupling using OASIS.
      real, pointer :: C_SOLAR(:)       ! CPL solar radn
      real, pointer :: C_BLUE(:)        ! CPL blue radn
      real, pointer :: C_LONGWAVE(:)    ! CPL lw radn
      real, pointer :: C_TAUX(:)        ! CPL taux
      real, pointer :: C_TAUY(:)        ! CPL tauy
      real, pointer :: C_W10(:)         ! CPL 10m wind
      real, pointer :: C_SENSIBLE(:)    ! CPL sensible ht flx
      real, pointer :: C_SUBLIM(:)      ! CPL Sublimation rate
      real, pointer :: C_EVAP(:)        ! CPL Evaporation
      real, pointer :: C_FCONDTOPN(:)   ! CPL Fcondtop
      real, pointer :: C_TOPMELTN(:)    ! CPL Topmelt
      real, pointer :: C_LSRAIN(:)      
      real, pointer :: C_LSSNOW(:)
      real, pointer :: C_CVRAIN(:)
      real, pointer :: C_CVSNOW(:)
      real, pointer :: C_RIVEROUT(:)
      REAL, POINTER :: c_calving(:)

!   2.18: JULES variables
      REAL, POINTER :: SNOWDEPTH(:)      ! Snow depth on ground on tiles (m)
      REAL, POINTER :: RHO_SNOW_GRND(:)  ! Snowpack bulk density (kg/m3)
      REAL, POINTER :: NSNOW(:)          ! Number of snow layers on ground on tiles
      REAL, POINTER :: DS(:)             ! Snow layer thickness (m)
      REAL, POINTER :: SICE(:)           ! Snow layer ice mass on tiles (Kg/m2)
      REAL, POINTER :: SLIQ(:)           ! Snow layer liquid mass on tiles (Kg/m2)
      REAL, POINTER :: TSNOWLAYER(:)     ! Snow layer temperature (K)
      REAL, POINTER :: RHO_SNOW(:)       ! Snow layer densities (kg/m3)
      REAL, POINTER :: RGRAINL(:)        ! Snow layer grain size on tiles (microns)
!  FLake lake scheme
      REAL, POINTER :: lake_depth(:)
      REAL, POINTER :: lake_fetch(:)
      REAL, POINTER :: lake_t_mean(:)
      REAL, POINTER :: lake_t_mxl(:)
      REAL, POINTER :: lake_t_ice(:)
      REAL, POINTER :: lake_h_mxl(:)
      REAL, POINTER :: lake_h_ice(:)
      REAL, POINTER :: lake_shape(:)
      REAL, POINTER :: lake_g_dt(:)

! Aerosol climatologies
      real, pointer :: ARCLBIOG_BG(:)    ! Biogenic aerosol climatology
      real, pointer :: ARCLBIOM_FR(:)    ! Biomass burning (fresh) aerosol clim
      real, pointer :: ARCLBIOM_AG(:)    ! Biomass burning (aged) aerosol clim
      real, pointer :: ARCLBIOM_IC(:)    ! Biomass burning (in-cloud) aerosol clim
      real, pointer :: ARCLBLCK_FR(:)    ! Black carbon (fresh) aerosol clim
      real, pointer :: ARCLBLCK_AG(:)    ! Black carbon (aged) aerosol clim
      real, pointer :: ARCLSSLT_FI(:)    ! Sea salt (film mode) aerosol clim 
      real, pointer :: ARCLSSLT_JT(:)    ! Sea salt (jet mode) aerosol clim
      real, pointer :: ARCLSULP_AC(:)    ! Sulphate (accumulation mode) aero clim
      real, pointer :: ARCLSULP_AK(:)    ! Sulphate (Aitken mode) aerosol clim 
      real, pointer :: ARCLSULP_DI(:)    ! Sulphate (dissolved) aerosol clim
      real, pointer :: ARCLDUST_B1(:)    ! Dust (bin 1) aerosol climatology 
      real, pointer :: ARCLDUST_B2(:)    ! Dust (bin 2) aerosol climatology 
      real, pointer :: ARCLDUST_B3(:)    ! Dust (bin 3) aerosol climatology 
      real, pointer :: ARCLDUST_B4(:)    ! Dust (bin 4) aerosol climatology 
      real, pointer :: ARCLDUST_B5(:)    ! Dust (bin 5) aerosol climatology 
      real, pointer :: ARCLDUST_B6(:)    ! Dust (bin 6) aerosol climatology 
      real, pointer :: ARCLOCFF_FR(:)    ! Org carbon fossil fuel (fresh) aero clim
      real, pointer :: ARCLOCFF_AG(:)    ! Org carbon fossil fuel (aged) aero clim
      real, pointer :: ARCLOCFF_IC(:)    ! Org carbon fossil fuel (in-cloud) aero clim
      real, pointer :: ARCLDLTA_DL(:)    ! Delta aerosol climatology

! Convective Cloud Fields
      real, pointer :: LCBASE(:)
      real, pointer :: CCW_RAD(:)

end module atm_fields_mod

! DECLARE_ATM_FIELDS_MOD END
