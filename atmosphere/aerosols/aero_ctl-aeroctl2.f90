! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE aero_ctl_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE aero_ctl(                                                           &
  ! Parallel variables
  halo_i, halo_j, off_x, off_y, global_row_length, global_rows,                &
  proc_row_group, proc_col_group, at_extremity, n_proc, n_procx,               &
  n_procy, neighbour, g_rows, g_row_length, me,                                &
  ! model dimensions
  row_length, rows, n_rows, land_points,                                       &
  model_levels, wet_model_levels, bl_levels, n_cca_levels,                     &
  theta_field_size,                                                            &
  salt_dim1, salt_dim2, salt_dim3,                                             &
  aero_dim1, aero_dim2, aero_dim3,                                             &
  ! Model switches
  model_domain, l_cal360, l_sec_var, l_eqt, ltimer,                            &
  ! Model parameters
  ntot_land, ntot_sea,                                                         &
  ! Co-ordinate information
  delta_lambda, delta_phi,                                                     &
  lat_rot_np, long_rot_np,                                                     &
  ! Time stepping information
  timestep,                                                                    &
  val_year, val_day_number, val_hour, val_minute,                              &
  val_second, timestep_number,                                                 &
  previous_time,                                                               &
  call_chem_freq,                                                              &
  ! Trig arrays
  sin_theta_longitude, cos_theta_longitude,                                    &
  fv_cos_theta_latitude,                                                       &
  ! Grid-dependent arrays
  f3_at_u, true_longitude, true_latitude,                                      &
  ! Data fields IN
  u, v, tstar, tstar_sea,                                                      &
  theta, q, qcl, qcf,                                                          &
  rho, land_mask, fland_ctile,                                                 &
  p_theta_levels, exner_rho_levels, exner_theta_levels,                        &
  ice_fract, snow_depth,                                                       &
  cloudf_halos,                                                                &
  oh_conc, h2o2_lmt, ho2_conc, o3, hno3,                                       &
  so2_surfem, so2_hilem, so2_natem,                                            &
  dms_em_ancil, dms_conc, nh3_em,                                              &
  dms_ointer,                                                                  &
  soot_surem, soot_hilem, bmass_surem, bmass_hilem, ocff_surem,                &
  ocff_hilem, land_index,                                                      &
  ! Logicals IN
  l_sulpc_so2, l_sulpc_dms, l_sulpc_ozone,                                     &
  l_sulpc_so2_o3_nonbuffered, l_sulpc_nh3,                                     &
  l_sulpc_online_oxidants,                                                     &
  l_sulpc_2_way_coupling,                                                      &
  l_sulphate_ccn, l_seasalt_ccn, l_soot, l_soot_ccn,                           &
  l_biomass, l_biomass_ccn, l_ocff, l_ocff_ccn,                                &
  l_nitrate, l_nitrate_ccn,                                                    &
  l_so2_surfem, l_so2_hilem, l_so2_natem, l_dms_em,                            &
  l_dms_em_inter, l_dms_liss_merlivat,                                         &
  l_dms_wanninkhof, l_dms_nightingale,                                         &
  l_dms_ointer,                                                                &
  l_nh3_em, l_ctile,                                                           &
  l_soot_surem, l_soot_hilem, l_bmass_surem, l_bmass_hilem,                    &
  l_ocff_surem, l_ocff_hilem, l_use_biogenic,                                  &
  l_use_seasalt_direct, l_use_seasalt_indirect,                                &
  l_use_seasalt_autoconv, l_use_seasalt_pm, l_dust,                            &
  ! Data fields IN/OUT
  so2, dms,                                                                    &
  so4_ait, so4_acc, so4_dis,                                                   &
  h2o2_mxr, nh3,                                                               &
  soot_new, soot_agd, soot_cld,                                                &
  bmass_new, bmass_agd, bmass_cld,                                             &
  ocff_new, ocff_agd, ocff_cld,                                                &
  biogenic, nitr_acc, nitr_diss,                                               &
  ! Data fields IN
  dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6,            &
  ! Data fields OUT
  ! Diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
  stashwork17,                                                                 &
  error_code)

!----------------------------------------------------------------------
! Purpose: Interface for Aerosol Modelling, to include Sulphur Cycle,
!          soot, biomass burning aerosol, OCFF and nitrate modelling
!          Note that dust is simply passed in for use in PM diagnostics
!          Main mineral dust routines are called from boundary layer
!          routines

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols

! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards

! System components covered:

! System task:

! Documentation: Not yet available

!-----------------------------------------------------------------
USE earth_constants_mod, ONLY: g , two_omega
USE atmos_constants_mod, ONLY: r

USE level_heights_mod, ONLY:                          &
    r_theta_levels, r_rho_levels

USE run_aerosol_mod, ONLY: so2_high_level, soot_high_level, &
    bmass_high_level_1,  bmass_high_level_2, ocff_high_level
USE calc_pm_diags_mod, ONLY: calc_pm_diags
USE global_2d_sums_mod, ONLY: global_2d_sums

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE agebmass_mod,         ONLY: agebmass
USE ageocff_mod,          ONLY: ageocff
USE agesoot_mod,          ONLY: agesoot
USE bmassnuclscav_mod,    ONLY: bmassnuclscav
USE diagnostics_aero_mod, ONLY: diagnostics_aero
USE dms_flux_mod,         ONLY: dms_flux
USE nitrate_mod,          ONLY: nitrate
USE ocffnuclscav_mod,     ONLY: ocffnuclscav
USE sootdiffscav_mod,     ONLY: sootdiffscav
USE sulphr_mod,           ONLY: sulphr
USE trsrce_mod,           ONLY: trsrce
USE number_droplet_mod,   ONLY: number_droplet
USE Submodel_Mod

IMPLICIT NONE

! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_S,                                                       &
                                 ! relative molecular mass S kg/mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_HO2,                                                     &
                                 ! relative molecular mass HO2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_S    = 3.20E-2,                                    &
                                          ! kg/mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &           RMM_HO2  = 3.30E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------

!
! Arguments with intent IN:

! Parallel setup variables
INTEGER :: halo_i                     ! Size of halo in i direction.
INTEGER :: halo_j                     ! Size of halo in j direction.
INTEGER :: off_x                      ! Size of small halo in i
INTEGER :: off_y                      ! Size of small halo in j.
INTEGER :: global_row_length          ! number of points on a row
INTEGER :: proc_row_group             ! Group id for processors on the same row
INTEGER :: proc_col_group             ! Group id for processors on the same col
INTEGER :: global_rows                ! NUMBER OF global rows
INTEGER :: n_proc                     ! Total number of processors
INTEGER :: n_procx                    ! Number of processors in longitude
INTEGER :: n_procy                    ! Number of processors in latitude
! Array with the Ids of the four neighbours in the horizontal plane
INTEGER :: neighbour(4)
INTEGER :: g_rows (0:n_proc-1)
INTEGER :: g_row_length (0:n_proc-1)
INTEGER :: me                         ! My processor number

! Indicates if this processor is at north, south east or west of the processor grid
LOGICAL  ::  at_extremity(4)


! Model dimensions
INTEGER :: row_length
INTEGER :: rows
INTEGER :: n_rows
INTEGER :: land_points
INTEGER :: model_levels
INTEGER :: wet_model_levels
INTEGER :: bl_levels
INTEGER :: n_cca_levels        ! No. conv cloud levels (1 if 2D, nlevs if 3D)
INTEGER :: theta_field_size
! dimensions of seasalt array
INTEGER :: salt_dim1
INTEGER :: salt_dim2
INTEGER :: salt_dim3
INTEGER :: aero_dim1           ! row_length or 1
INTEGER :: aero_dim2           ! rows or 1
INTEGER :: aero_dim3           ! model_levels or 1
! Model switches
INTEGER :: model_domain
LOGICAL :: l_cal360             ! T if using 360 day calendar
LOGICAL :: l_sec_var            ! if T include secular varn of earth's orbit,
LOGICAL :: l_eqt                ! if T include equn of time in SOLPOS
LOGICAL :: ltimer               ! if T then output some timing information

! Co-ordinate arrays
REAL :: delta_lambda
REAL :: delta_phi
REAL :: lat_rot_np
REAL :: long_rot_np

! Trig arrays
REAL :: cos_theta_longitude (row_length, rows)
REAL :: sin_theta_longitude (row_length, rows)
! Grid-dependent arrays
REAL :: f3_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)
REAL :: true_longitude(row_length, rows)
REAL :: true_latitude(row_length, rows)
REAL :: fv_cos_theta_latitude (1-off_x:row_length+off_x,1-off_y:rows+off_y)

! Time stepping information
REAL    :: timestep           !atmosphere model timestep
INTEGER :: val_year
INTEGER :: val_day_number
INTEGER :: val_hour
INTEGER :: val_minute
INTEGER :: val_second
INTEGER :: timestep_number
INTEGER :: previous_time(7)
! frequency of calling chemistry per atmos phys timestep
INTEGER :: call_chem_freq


! Diagnostics info
REAL :: stashwork17(*)  ! STASH workspace for section 17 (Aero_Ctl)

REAL :: u    (1-off_x:row_length+off_x, 1-off_y:rows+off_y,    model_levels)
REAL :: v    (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,  model_levels)
REAL :: theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    model_levels)
REAL :: q    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,wet_model_levels)
REAL :: qcl  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,wet_model_levels)
REAL :: qcf  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,wet_model_levels)
REAL :: rho  (1-off_x:row_length+off_x,   1-off_y:rows+off_y,  model_levels)
REAL :: p_theta_levels(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)
REAL :: exner_rho_levels  (1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels+1)
REAL :: exner_theta_levels(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)
REAL :: ice_fract (row_length, rows)
REAL :: snow_depth(row_length, rows)
REAL :: land_fract(row_length, rows)
REAL :: tstar(row_length, rows)
REAL :: tstar_sea(row_length, rows)
REAL :: fland_ctile(land_points)
REAL :: cloudf_halos(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,wet_model_levels)
REAL :: oh_conc(row_length,rows,model_levels)
REAL :: ho2_conc(row_length,rows,model_levels)
REAL :: h2o2_lmt(row_length,rows,model_levels)
REAL :: o3(row_length,rows,model_levels)
REAL :: so2_surfem(row_length,rows)                    ! SO2 emiss at surface
REAL :: so2_hilem(row_length,rows)                     ! SO2 emiss at chimney level
REAL :: so2_natem(row_length,rows,model_levels)        ! Volcanic SO2 emiss
REAL :: dms_em_ancil(row_length,rows)                  ! Ancillary DMS emiss (surf)
REAL :: dms_conc(row_length,rows)                      ! Seawater DMS conc (n mol l-1)
REAL :: nh3_em(row_length,rows)                        ! NH3 emiss (surface)
REAL :: soot_surem(row_length,rows)                    ! SOOT emiss at surface
REAL :: soot_hilem(row_length,rows)                    ! SOOT emiss at chimney lev
REAL :: bmass_surem(row_length,rows)                   ! Biomass surface emissions
REAL :: bmass_hilem(row_length,rows)                   ! Biomass high level emissions
REAL :: ocff_surem(row_length,rows)                    ! OCFF emiss at surface
REAL :: ocff_hilem(row_length,rows)                    ! OCFF emiss at chimney lev

REAL ::  ntot_land                      ! Number of droplets over land / m-3
REAL ::  ntot_sea                       ! Number of droplets over sea / m-3

INTEGER :: land_index(land_points)

LOGICAL :: l_sulpc_so2                  ! T if Sulphur Cycle required
LOGICAL :: l_sulpc_dms                  ! T if DMS chemistry required
LOGICAL :: l_sulpc_ozone                ! T if O3 oxidn required

! l_sulpc_so2_o3_nonbuffered T if SO2+O3 reaction is NOT to be buffered by NH3.
LOGICAL :: l_sulpc_so2_o3_nonbuffered

! l_sulpc_nh3 T if NH3 buffering required (always T if L_SULPC_OZONE is T)
LOGICAL :: l_sulpc_nh3
LOGICAL :: l_sulpc_online_oxidants      ! T if oxidants from UKCA are used.

! l_sulpc_2_way_coupling T if coupling to UKCA is 2-way,
! i.e. if depleted oxidants are passed back to UKCA
LOGICAL :: l_sulpc_2_way_coupling
LOGICAL :: l_sulphate_ccn               !T if sulphate used for CCN
LOGICAL :: l_seasalt_ccn                !T if sea-salt used for CCN
LOGICAL :: l_soot_ccn                   !T if soot used for CCN
LOGICAL :: l_biomass_ccn                !T if biomass used for CCN
LOGICAL :: l_ocff_ccn                   !T if OCFF used for CCN
LOGICAL :: l_nitrate_ccn                !T if nitrate used for CCN
LOGICAL :: l_ctile                      !T if coastal tiling is on
LOGICAL :: l_soot                       !T if SOOT modelling required
LOGICAL :: l_biomass                    !T if biomass modelling reqd
LOGICAL :: l_ocff                       !T if OCFF modelling required
LOGICAL :: l_nitrate                    !T if nitrate modelling required
LOGICAL :: l_so2_surfem                 !T if surface SO2 ems present
LOGICAL :: l_so2_hilem                  !T if high lev SO2 em present
LOGICAL :: l_so2_natem                  !T if volcanic SO2 ems present
LOGICAL :: l_dms_em                     !T if DMS emiss present (surf)
LOGICAL :: l_dms_em_inter               !T if interactive DMS emiss
LOGICAL :: l_dms_ointer                 !T if ocean DMS emiss model

! Switches to determine which scheme to use for interactive DMS emissions
LOGICAL :: l_dms_liss_merlivat
LOGICAL :: l_dms_wanninkhof
LOGICAL :: l_dms_nightingale

LOGICAL :: l_nh3_em                     ! T if NH3 emiss present (surf)
LOGICAL :: l_soot_surem                 ! T if surface SOOT ems present
LOGICAL :: l_soot_hilem                 ! T if high lev SOOT ems presnt
LOGICAL :: l_bmass_surem                ! T if sfc biomass ems present
LOGICAL :: l_bmass_hilem                ! T if hi lev bmass ems present
LOGICAL :: l_ocff_surem                 ! T if surface OCFF ems present
LOGICAL :: l_ocff_hilem                 ! T if high lev OCFF ems presnt
LOGICAL :: l_use_biogenic               ! T if using biogenics for CCN
LOGICAL :: land_mask(row_length,rows)   ! T IF LAND, F IF SEA
LOGICAL :: l_dust                       ! T if mineral dust used
LOGICAL :: l_use_seasalt_direct         ! T if SS dir. rad. effect.
LOGICAL :: l_use_seasalt_indirect       ! T if SS 1st indir. effect
LOGICAL :: l_use_seasalt_autoconv       ! T if SS 2nd indir. effect
LOGICAL :: l_use_seasalt_pm

! Arguments with intent IN/OUT:
! mmr S in SO2
REAL :: so2      (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
! mmr S in DMS
REAL :: dms      (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
! mmr S in AIT
REAL, TARGET :: so4_ait  (1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          model_levels)
! mmr S in ACC
REAL, TARGET :: so4_acc  (1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          model_levels)
! mmr S in DIS
REAL, TARGET :: so4_dis  (1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          model_levels)
! mmr N in NH3
REAL :: nh3      (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
! mmr H2O2
REAL :: h2o2_mxr (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)
! mmr fresh soot
REAL :: soot_new (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
! mmr aged soot
REAL :: soot_agd (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
! mmr soot in cloud
REAL :: soot_cld (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
! mmr fresh smoke
REAL :: bmass_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
! mmr aged smoke
REAL, TARGET :: bmass_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          model_levels)
! mmr cloud smoke
REAL, TARGET :: bmass_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          model_levels)
! mmr fresh OCFF
REAL :: ocff_new (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
! mmr aged OCFF
REAL, TARGET :: ocff_agd (1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          model_levels)
! mmr OCFF in cloud
REAL, TARGET :: ocff_cld (1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          model_levels)

REAL, TARGET :: biogenic(row_length, rows, model_levels)    ! mmr biogenics
REAL :: hno3(row_length,rows,model_levels)                  ! mmr HNO3

! mmr accum mode nitrate
REAL, TARGET :: nitr_acc (1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          model_levels)
! mmr dissolved nitrate
REAL, TARGET :: nitr_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          model_levels)

! Arguments with intent IN:
! Dust in 6 divisions
REAL, INTENT(IN) :: dust_div1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div4(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div5(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div6(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)

! arguments with intent OUT:
INTEGER ::   error_code

! Variables with intent OUT (diagnostics):
REAL :: msa(row_length,rows,model_levels)             ! mmr S in MSA
REAL :: f_dms_to_so2(row_length,rows,model_levels)    ! frac oxid DMS to SO2
REAL :: f_dms_to_so4(row_length,rows,model_levels)    ! frac oxid DMS to SO4
REAL :: f_dms_to_msa(row_length,rows,model_levels)    ! frac oxid DMS to MSA
REAL :: nh3_dep(row_length,rows,model_levels)         ! NH3 depleted
REAL :: deltas_dry(row_length,rows,model_levels)      ! SO2 dry ox per ts
REAL :: deltas_wet(row_length,rows,model_levels)      ! SO2 wet ox by H2O2
REAL :: deltas_wet_o3(row_length,rows,model_levels)   ! SO2 wet ox by O3
REAL :: deltas_tot(row_length,rows,model_levels)      ! total SO2 ox per ts
REAL :: deltas_dms(row_length,rows,model_levels)      !DMS dry ox per ts
! SO4_DIS released by evapn of cloud droplets to SO4_ACC per ts
REAL :: deltas_evap(row_length,rows,model_levels)
! SO4_ACC transfd by nucleation to SO4_DIS per ts
REAL :: deltas_nucl(row_length,rows,model_levels)
!SO4_AIT transfd to SO4_DIS by diffusion per ts
REAL :: deltas_diffuse(row_length,rows,model_levels)
REAL :: deltas_merge(row_length,rows,model_levels)
REAL :: deltas_coag(row_length,rows,model_levels)
REAL :: psi(row_length,rows,model_levels)
! Flux thro' reaction HNO3 + NH3 --> NH4NO3
REAL :: delta_n_chem(row_length,rows,model_levels)
! nitr_diss released by evaporation of cloud droplets to nitr_diss per ts
REAL :: delta_n_evap(row_length,rows,model_levels)
! nitr_acc transferred by nucleation to nitr_diss by nucleation per ts.
REAL :: delta_n_nuc(row_length,rows,model_levels)

! Local variables

INTEGER :: i,j,k,n
INTEGER :: i_start              ! Row start point for polar row tidying
INTEGER :: istat                ! Status (error code) indicator

REAL    :: chemstep             ! chemistry timestep

! Diagnostic increments generated by each call of SULPHR
REAL :: msa_inc    (row_length,rows,model_levels)
REAL :: nh3_dep_inc(row_length,rows,model_levels)

REAL :: t     (row_length,rows,model_levels)   ! Temp (K) on theta levels
REAL :: t_surf(row_length,rows)                ! Surface temperature (K)

! For calculation of cos zenith angle
REAL :: cosza2d(row_length,rows)      !cos zenith angle
REAL :: sindec                        ! sin(solar declination)
REAL :: hour_angle(row_length,rows)
REAL :: cosz_beg(row_length,rows)
REAL :: cosz_end(row_length,rows)
REAL :: scs                           ! solar constant scaling factor
REAL :: seconds_since_midnight
REAL :: eq_time                       ! The Equation of time
REAL :: day_fraction(row_length,rows)

! For sea-salt aerosol and DMS emissions calculations
REAL :: u_1(aero_dim1,aero_dim2)                      ! surface wind u
REAL :: v_1(aero_dim1,aero_dim2)                      ! surface wind v
REAL :: u_1_mean(1)                                   ! mean of u_1 for N Pole
REAL :: v_1_mean(1)                                   ! mean of v_1 for N Pole
! Sea-salt film aerosol.
REAL, TARGET :: sea_salt_film(salt_dim1,salt_dim2,salt_dim3)  
! Sea-salt jet aerosol.
REAL, TARGET :: sea_salt_jet(salt_dim1,salt_dim2,salt_dim3)   
REAL :: height(aero_dim1,aero_dim2,aero_dim3)         ! Layer centre heights
REAL :: dms_em_inter(aero_dim1,aero_dim2)             ! Interactive DMS emissns
REAL :: dms_ointer(aero_dim1,aero_dim2)               ! Interactive DMS emissns
REAL :: dms_ointer_mean(1)                            ! N Pole mean DMS_Ointer
REAL :: dms_em_combined(row_length,rows)              ! Combined DMS emissns

! Increments from soot ageing and diffusional scavenging
! Increment from soot ageing
REAL :: delta_agesoot(row_length,rows,model_levels)
! Increment from soot diff. scav.
REAL :: delta_sootdiffscav(row_length,rows,model_levels)

! Increments from biomass smoke ageing and nucleation scavenging
REAL :: delta_agebmass(row_length,rows,model_levels)      ! Increment from smoke ageing
REAL :: delta_bmassnuclscav(row_length,rows,model_levels) ! Increment from smoke nucl. scav.
! For emissions of biomass smoke
INTEGER ::  biomass_level                                  ! Model level of biomass emissions
REAL    ::  biomass_per_level(row_length,rows)
! Quantity of biomass smoke emitted
! per model level

! Increments from fossil-fuel organic carbon ageing and nucl. scavenging
REAL :: delta_ageocff(row_length, rows, model_levels)      ! Increment from OCFF ageing
REAL :: delta_ocffnuclscav(row_length, rows, model_levels) ! Increment from OCFF nucl. scav.


! For cloud fraction without halos
REAL :: cloudf(row_length,rows,wet_model_levels)

LOGICAL,PARAMETER :: l_sulpc_newdms= .TRUE.       ! TRUE if new DMS scheme required

! For droplet number calculation
REAL, TARGET :: dummy(row_length, rows, wet_model_levels)
REAL :: n_droplet(row_length, rows, wet_model_levels)
REAL :: rho_air(row_length, rows, wet_model_levels)
REAL, POINTER :: ss_film(:,:,:)
REAL, POINTER :: ss_jet(:,:,:)
REAL, POINTER :: sulp_acc(:,:,:)
REAL, POINTER :: sulp_dis(:,:,:)
REAL, POINTER :: agd_bmass(:,:,:)
REAL, POINTER :: cld_bmass(:,:,:)
REAL, POINTER :: agd_ocff(:,:,:)
REAL, POINTER :: cld_ocff(:,:,:)
REAL, POINTER :: acc_nitrate(:,:,:)
REAL, POINTER :: dis_nitrate(:,:,:)
REAL, POINTER :: biogenic_ccn(:,:,:)

! For PM10 & PM2.5 calculation
REAL :: pm10       (row_length, rows, model_levels)
REAL :: pm2p5      (row_length, rows, model_levels)
REAL :: pm10_so4   (row_length, rows, model_levels)
REAL :: pm2p5_so4  (row_length, rows, model_levels)
REAL :: pm10_bc    (row_length, rows, model_levels)
REAL :: pm2p5_bc   (row_length, rows, model_levels)
REAL :: pm10_bb    (row_length, rows, model_levels)
REAL :: pm2p5_bb   (row_length, rows, model_levels)
REAL :: pm10_ocff  (row_length, rows, model_levels)
REAL :: pm2p5_ocff (row_length, rows, model_levels)
REAL :: pm10_soa   (row_length, rows, model_levels)
REAL :: pm2p5_soa  (row_length, rows, model_levels)
REAL :: pm10_ss    (row_length, rows, model_levels)
REAL :: pm2p5_ss   (row_length, rows, model_levels)
REAL :: conc_dust  (row_length, rows, model_levels)
REAL :: pm10_dust  (row_length, rows, model_levels)
REAL :: pm2p5_dust (row_length, rows, model_levels)
REAL :: pm10_nitr  (row_length, rows, model_levels)
REAL :: pm2p5_nitr (row_length, rows, model_levels)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! Set up global land fraction field

IF (lhook) CALL dr_hook('AERO_CTL',zhook_in,zhook_handle)

!$OMP  PARALLEL DEFAULT(SHARED)                                                &
!$OMP& PRIVATE(k, j, i,  ss_film, ss_jet, sulp_acc,                            &
!$OMP& agd_bmass, cld_bmass, agd_ocff, dis_nitrate, biogenic_ccn,              &
!$OMP& acc_nitrate, sulp_dis, cld_ocff)

IF (l_ctile) THEN

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      land_fract(i,j) = 0.0
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k = 1,land_points
    j = (land_index(k)-1)/row_length + 1
    i = land_index(k) - (j-1)*row_length
    land_fract(i,j) = fland_ctile(k)
  END DO
!$OMP END DO NOWAIT

ELSE

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      IF (land_mask(i,j)) THEN
        land_fract(i,j) = 1.0
      ELSE
        land_fract(i,j) = 0.0
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT

END IF


! Copy cloud fraction into an array without halos and set to zero if
! smaller than machine precision

!$OMP DO SCHEDULE(STATIC)
DO k = 1,wet_model_levels
  DO j = 1,rows
    DO i = 1,row_length
      cloudf(i,j,k)=cloudf_halos(i,j,k)
      IF (cloudf(i,j,k) <  EPSILON(cloudf(i,j,k)))                             &
        cloudf(i,j,k)=0.0
      IF (cloudf(i,j,k) >  1.0) THEN
        cloudf(i,j,k) = 1.0
      END IF
    END DO
  END DO
END DO
!$OMP END DO NOWAIT


! If sea-salt aerosol or interactive DMS emissions are to be
! used then calculate surface windspeed and height information

IF (l_seasalt_ccn .OR. l_use_seasalt_pm .OR. l_dms_em_inter) THEN

!$OMP DO SCHEDULE(STATIC)
  DO j = 1,rows
    DO i = 1,row_length
      u_1(i,j) = u(i,j,1)
      v_1(i,j) = v(i,j,1)
    END DO
  END DO
!$OMP END DO



!$OMP MASTER

!
! Tidy up at North Pole; not done at South Pole as it's land.
! Polar values are mean of 1st interior row:
        If (at_extremity(PNorth)) Then
! Start point of first interior (i.e. non-polar) row:
          i_start = rows-1
! Sum over points on PEs in order along first interior row:
          CALL global_2d_sums(u_1(:,i_start:i_start), row_length, 1,   &
                              0, 0, 1, u_1_mean, proc_row_group)
          CALL global_2d_sums(v_1(:,i_start:i_start), row_length, 1,   &
                              0, 0, 1, v_1_mean, proc_row_group)

          u_1_mean(1)=u_1_mean(1)/global_row_length
          v_1_mean(1)=v_1_mean(1)/global_row_length
          Do i = 1, row_length
            u_1(i, rows) = u_1_mean(1)
            v_1(i, rows) = v_1_mean(1)
          End Do
        Endif
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, aero_dim3
    DO j = 1, aero_dim2
      DO i = 1, aero_dim1
        height(i,j,k) = r_theta_levels(i,j,k)                                  &
          -r_theta_levels(i,j,0)
      END DO
    END DO
  END DO
!$OMP END DO


END IF ! L_seasalt_CCN .OR. L_DMS_em_inter


!$OMP MASTER

      If (L_DMS_Ointer) Then
!
! Tidy up ocean DMS flux at North Pole, as done above for u_1,v_1.
! Polar values are mean of 1st interior row:
        If (at_extremity(PNorth)) Then
! Start point of first interior (i.e. non-polar) row:
          i_start = rows-1
! Sum over points on PEs in order along first interior row:
          CALL global_2d_sums(dms_ointer(:,i_start:i_start), row_length,&
                              1, 0, 0, 1, dms_ointer_mean,              &
                              proc_row_group)

          dms_Ointer_mean(1)=dms_Ointer_mean(1)/global_row_length
          Do i = 1, row_length
            dms_Ointer(i, rows) = dms_Ointer_mean(1)
          End Do
        Endif
      Endif ! L_DMS_Ointer
!
   

IF (l_seasalt_ccn .OR. l_use_seasalt_pm) THEN
  ! DEPENDS ON: set_seasalt
  CALL set_seasalt(u_1, v_1, height, land_fract, ice_fract,                    &
    row_length, rows, model_levels,                                            &
    salt_dim1, salt_dim2, salt_dim3,                                           &
    bl_levels, sea_salt_film, sea_salt_jet)
ELSE
  sea_salt_film(1,1,1) = 0.0
  sea_salt_jet(1,1,1) = 0.0
END IF

!$OMP END MASTER
!$OMP BARRIER

! Calculate T from theta if soot, S cycle, biomass or OCFF switched on
! and droplet numbers for diffusional scavenging.

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      t(i,j,k)=theta(i,j,k)*exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO

IF (l_sulpc_so2 .OR. l_soot .OR. l_biomass .OR. l_ocff) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k=1,wet_model_levels
    DO j=1,rows
      DO i=1,row_length
        rho_air(i,j,k)=p_theta_levels(i,j,k)/                                  &
          (r * t(i,j,k))
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

END IF ! L_SULPC_SO2 .OR. L_SOOT .OR. L_BIOMASS .OR. L_OCFF

!$OMP END PARALLEL

IF (l_sulpc_so2 .OR. l_soot .OR. l_biomass .OR. l_ocff) THEN

  IF (l_seasalt_ccn) THEN
    ss_film=>sea_salt_film(1:row_length,1:rows,:)
    ss_jet=>sea_salt_jet(1:row_length,1:rows,:)
  ELSE
    ss_film=>dummy
    ss_jet=>dummy
  END IF

  IF (l_sulphate_ccn) THEN
    sulp_acc=>so4_acc(1:row_length,1:rows,:)
    sulp_dis=>so4_dis(1:row_length,1:rows,:)
  ELSE
    sulp_acc=>dummy
    sulp_dis=>dummy
  END IF

  IF (l_biomass_ccn) THEN
    agd_bmass=>bmass_agd(1:row_length,1:rows,:)
    cld_bmass=>bmass_cld(1:row_length,1:rows,:)
  ELSE
    agd_bmass=>dummy
    cld_bmass=>dummy
  END IF

  IF (l_ocff_ccn) THEN
    agd_ocff=>ocff_agd(1:row_length,1:rows,:)
    cld_ocff=>ocff_cld(1:row_length,1:rows,:)
  ELSE
    agd_ocff=>dummy
    cld_ocff=>dummy
  END IF

  IF (l_nitrate_ccn) THEN
    acc_nitrate=>nitr_acc(1:row_length,1:rows,:)
    dis_nitrate=>nitr_diss(1:row_length,1:rows,:)
  ELSE
    acc_nitrate=>dummy
    dis_nitrate=>dummy
  END IF

  IF (l_use_biogenic) THEN
    biogenic_ccn=>biogenic(1:row_length,1:rows,:)
  ELSE
    biogenic_ccn=>dummy
  END IF

  CALL number_droplet(                                                   &
    1, row_length,                                                       &
    1, rows,                                                             &
    1, model_levels,                                                     &
    1, model_levels,                                                     &
    l_sulphate_ccn,.FALSE.,                                              &
    sulp_acc,sulp_dis,                                                   &
    l_seasalt_ccn,ss_film,ss_jet,                                        &
    l_use_biogenic,biogenic_ccn,                                         &
    l_biomass_ccn,agd_bmass,cld_bmass,                                   &
    l_ocff_ccn,agd_ocff,cld_ocff,                                        &
    l_nitrate_ccn,acc_nitrate,dis_nitrate,                               &
    rho_air,snow_depth,land_fract,                                       &
    ntot_land,ntot_sea, n_droplet )

END IF

! Soot cycle code.
! Calculate emissions of soot from surface and chimney level, ageing
! of soot and diffusional scavenging of aged soot to soot-in-cloud.

IF (l_soot) THEN  ! If soot modelling is included

  IF (l_soot_surem) THEN  ! If surface soot emissions are included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      model_levels, wet_model_levels,                                          &
      halo_i, halo_j,                                                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      soot_new(:,:,1), soot_surem, 1,                                          &
      timestep, 1, 1, 0.0 )
  END IF  ! L_SOOT_SUREM

  IF (l_soot_hilem) THEN  ! If chimney level soot emissions
    ! are included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      model_levels, wet_model_levels,                                          &
      halo_i, halo_j,                                                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      soot_new(:,:,soot_high_level), soot_hilem,                               &
      soot_high_level, timestep, 1, 1, 0.0 )
  END IF  ! L_SOOT_HILEM

  ! Calculate quantity of fresh soot converted to aged soot

  CALL agesoot(                                                                &
    row_length, rows, off_x, off_y,                                            &
    model_levels, timestep,                                                    &
    soot_new,                                                                  &
    delta_agesoot                                                              &
    )

  ! Calculate quantity of aged soot scavenged to cloud soot

  CALL sootdiffscav(                                                           &
    rows, row_length, off_x, off_y, halo_i, halo_j,                            &
    model_levels, wet_model_levels, timestep,                                  &
    cloudf, qcl, qcf, p_theta_levels, t,                                       &
    n_droplet,                                                                 &
    soot_agd, soot_cld,                                                        &
    delta_sootdiffscav                                                         &
    )

  ! Update soot arrays with increments from ageing and
  ! diffusional scavenging

  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        soot_new(i,j,k)=soot_new(i,j,k) - delta_agesoot(i,j,k)
        soot_agd(i,j,k)=soot_agd(i,j,k) + delta_agesoot(i,j,k)                 &
          - delta_sootdiffscav(i,j,k)
        soot_cld(i,j,k)=soot_cld(i,j,k)                                        &
          + delta_sootdiffscav(i,j,k)
      END DO
    END DO
  END DO

END IF  ! L_SOOT

! End of soot cycle code

! Biomass aerosol code.
! Calculate emissions of smoke from surface and high level, ageing
! of smoke and nucleation scavenging of aged smoke to smoke-in-cloud.

IF (l_biomass) THEN  ! If biomass smoke modelling is included

  IF (l_bmass_surem) THEN  ! If surface smoke emissions included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      model_levels, wet_model_levels,                                          &
      halo_i, halo_j,                                                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      bmass_new(:,:,1), bmass_surem, 1,                                        &
      timestep, 1, 1, 0.0  )
  END IF  ! L_BMASS_SUREM

  IF (l_bmass_hilem) THEN  ! If high level smoke emissions
    ! are included

    ! Check that the range of emission levels is correctly specified
    ! and, if so, emit equal quantities of smoke on all model levels
    ! between bmass_high_level_1 and bmass_high_level_2.

    IF (bmass_high_level_1  >   bmass_high_level_2 .OR.                        &
      bmass_high_level_1  >   model_levels .OR.                                &
      bmass_high_level_2  >   model_levels) THEN
      WRITE(6,*) 'Aero_Ctl: Invalid range of biomass emission '
      WRITE(6,*) 'levels specified.'
      WRITE(6,*) 'Lowest level: ',bmass_high_level_1
      WRITE(6,*) 'Highest level: ',bmass_high_level_2
    ELSE

      DO j=1,rows
        DO i=1,row_length
          biomass_per_level(i,j) = bmass_hilem(i,j)/                           &
            ((bmass_high_level_2 - bmass_high_level_1) + 1)
        END DO
      END DO

      DO biomass_level = bmass_high_level_1, bmass_high_level_2

        CALL trsrce(                                                           &
          rows, row_length, off_x, off_y, halo_i, halo_j,                      &
          model_levels, wet_model_levels,                                      &
          halo_i, halo_j,                                                      &
          theta, q , qcl , qcf , exner_rho_levels, rho,                        &
          bmass_new(:,:,biomass_level), biomass_per_level,                     &
          biomass_level, timestep, 1, 1, 0.0 )

      END DO  ! biomass_level

    END IF  ! Test on emissions levels

  END IF  ! L_BMASS_HILEM

  ! Calculate quantity of fresh smoke converted to aged smoke

  CALL agebmass(                                                               &
    row_length, rows, off_x, off_y,                                            &
    model_levels, timestep,                                                    &
    bmass_new,                                                                 &
    delta_agebmass                                                             &
    )

  ! Calculate quantity of aged smoke scavenged to cloud smoke

  CALL bmassnuclscav(                                                          &
    rows, row_length, off_x, off_y, halo_i, halo_j,                            &
    model_levels, wet_model_levels, timestep,                                  &
    cloudf, qcl, qcf,                                                          &
    bmass_agd, bmass_cld,                                                      &
    delta_bmassnuclscav                                                        &
    )

  ! Update smoke arrays with increments from ageing and
  ! nucleation scavenging

  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        bmass_new(i,j,k)=bmass_new(i,j,k) - delta_agebmass(i,j,k)
        !    Simulate the condensation of VOCs onto aged biomass by
        !    increasing the mass transferred upon ageing by (8.75/5.4),
        !    the ratio of the fraction of BC in fresh (8.75%) and aged
        !    (5.4%) biomass aerosol:
        bmass_agd(i,j,k)=bmass_agd(i,j,k)                                      &
          + ((8.75/5.4)*delta_agebmass(i,j,k))                                 &
          - delta_bmassnuclscav(i,j,k)
        bmass_cld(i,j,k)=bmass_cld(i,j,k)                                      &
          + delta_bmassnuclscav(i,j,k)
      END DO
    END DO
  END DO

END IF  ! L_BIOMASS

! End of biomass aerosol code


IF ( l_sulpc_so2 ) THEN

  ! Calculate cos zenith angle for SULPH2 by calling SOLPOS and SOLANG
  ! Calculate number of seconds since midnight to the beginning of
  ! the timetsep.

  seconds_since_midnight = REAL( previous_time(4) * 3600                       &
    + previous_time(5) * 60  + previous_time(6))

  ! DEPENDS ON: solpos
  CALL solpos (previous_time(7), previous_time(1),                             &
    l_cal360, l_sec_var, l_eqt, eq_time, sindec, scs)

  ! DEPENDS ON: solang
  CALL solang(                                                                 &
                                ! input constants
    sindec, seconds_since_midnight,                                            &
    timestep, eq_time,                                                         &
                                ! row and column dependent constants
    true_latitude,                                                             &
    true_longitude,                                                            &
                                ! size variables
    row_length*rows,                                                           &
                                ! output fields
    day_fraction, cosza2d, hour_angle, cosz_beg, cosz_end )


  ! Calculate interactive DMS emissions if required

  IF (l_dms_em_inter) THEN

    DO j = 1, rows
      DO i = 1, row_length
        IF (l_ctile) THEN
          t_surf(i, j) = tstar_sea(i, j)
        ELSE
          t_surf(i, j) = tstar(i, j)
        END IF
      END DO
    END DO

    CALL dms_flux (                                                            &
      ! size variables
      row_length, rows,                                                        &
      ! input fields
      u_1, v_1,                                                                &
      height(:,:,1),                                                           &
      t_surf,                                                                  &
      land_fract,                                                              &
      dms_conc,                                                                &
      ! logical switches
      l_dms_liss_merlivat, l_dms_wanninkhof,  l_dms_nightingale,               &
      ! output field
      dms_em_inter )

  END IF


  ! Zero output diagnostics for full model timestep
  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        msa(i,j,k) = 0.0
        nh3_dep(i,j,k) = 0.0
      END DO
    END DO
  END DO

  ! Calculate length of chemistry timestep and use it to control input
  ! of emissions and S Chemistry

  chemstep=timestep/call_chem_freq

  DO n=1,call_chem_freq

    ! Call TRSRCE to insert emissions

    IF (l_so2_surfem) THEN         ! Insert surface SO2 emiss
      CALL trsrce(                                                             &
        rows, row_length, off_x, off_y, halo_i, halo_j,                        &
        model_levels, wet_model_levels,                                        &
        halo_i, halo_j,                                                        &
        theta, q , qcl , qcf , exner_rho_levels, rho,                          &
        so2(:,:,1), so2_surfem, 1,                                             &
        chemstep, 1, 1, 0.0 )
    END IF

    IF (l_so2_hilem) THEN          ! Insert chimney SO2 emiss
      CALL trsrce(                                                             &
        rows, row_length, off_x, off_y, halo_i, halo_j,                        &
        model_levels, wet_model_levels,                                        &
        halo_i, halo_j,                                                        &
        theta, q , qcl , qcf , exner_rho_levels, rho,                          &
        so2(:,:,so2_high_level), so2_hilem, so2_high_level,                    &
        chemstep, 1, 1, 0.0 )
    END IF

    IF (l_so2_natem) THEN          ! Insert volcanic SO2 emiss
      DO k= 1,model_levels
        CALL trsrce(                                                           &
          rows, row_length, off_x, off_y, halo_i, halo_j,                      &
          model_levels, wet_model_levels,                                      &
          halo_i, halo_j,                                                      &
          theta, q , qcl , qcf , exner_rho_levels, rho,                        &
          so2(:,:,k), so2_natem(:,:,k), k,                                     &
          chemstep, 1, 1, 0.0 )
      END DO
    END IF

    ! Merge ocean DMS emissions with land DMS emissions
    IF (l_dms_ointer) THEN
      dms_em_inter(:,:) = dms_ointer(:,:)
    END IF

    IF (l_dms_em_inter .OR. l_dms_ointer) THEN
      DO j = 1, rows
        DO i = 1, row_length
          dms_em_combined(i, j)                                                &
            = (land_fract(i, j) * dms_em_ancil(i, j))                          &
            + ( ((1.0 - land_fract(i, j)) * dms_em_inter(i, j))                &
            * (1.0 - ice_fract(i, j)) )
        END DO
      END DO
    ELSE
      IF (l_dms_em) THEN     ! Just copy over the standard ancil
        DO j = 1, rows
          DO i = 1, row_length
            dms_em_combined(i, j) = dms_em_ancil(i, j)
          END DO
        END DO
      END IF
    END IF

    IF (l_dms_em) THEN             ! Insert DMS emiss (surface)
      CALL trsrce(                                                             &
        rows, row_length, off_x, off_y, halo_i, halo_j,                        &
        model_levels, wet_model_levels,                                        &
        halo_i, halo_j,                                                        &
        theta, q , qcl , qcf , exner_rho_levels, rho,                          &
        dms(:,:,1), dms_em_combined, 1,                                        &
        chemstep, 1, 1, 0.0 )
    END IF

    IF (l_nh3_em) THEN             ! Insert NH3 emiss (surface)
      CALL trsrce(                                                             &
        rows, row_length, off_x, off_y, halo_i, halo_j,                        &
        model_levels, wet_model_levels,                                        &
        halo_i, halo_j,                                                        &
        theta, q , qcl , qcf , exner_rho_levels, rho,                          &
        nh3(:,:,1), nh3_em, 1,                                                 &
        chemstep, 1, 1, 0.0 )
    END IF


    CALL sulphr(                                                               &
                                ! Arguments IN
      halo_i, halo_j, off_x, off_y,                                            &
      row_length, rows,                                                        &
      model_levels, wet_model_levels ,                                         &
      theta_field_size,                                                        &
      chemstep,                                                                &
      cloudf, cosza2d,                                                         &
      p_theta_levels, t, q, qcl, qcf,                                          &
      oh_conc, h2o2_lmt, ho2_conc, o3,                                         &
      n_droplet,                                                               &
      l_sulpc_dms, l_sulpc_newdms,                                             &
      l_sulpc_ozone,                                                           &
      l_sulpc_so2_o3_nonbuffered,                                              &
      l_sulpc_nh3,                                                             &
      l_sulpc_online_oxidants,                                                 &
      l_sulpc_2_way_coupling,                                                  &
                                ! Arguments IN/OUT
      so2, dms, so4_ait, so4_acc, so4_dis,                                     &
      nh3, h2o2_mxr,                                                           &
                                ! Arguments OUT (diagnostics)
      msa_inc,                                                                 &
      nh3_dep_inc,                                                             &
      f_dms_to_so2, f_dms_to_so4, f_dms_to_msa,                                &
      deltas_dry, deltas_wet, deltas_wet_o3,                                   &
      deltas_tot, deltas_dms,                                                  &
      deltas_evap, deltas_nucl, deltas_diffuse,                                &
      deltas_merge,                                                            &
      deltas_coag, psi )

    ! Add diagnostic increments from SULPHR to total for model timestep

    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          msa(i,j,k)=msa(i,j,k)+msa_inc(i,j,k)
          nh3_dep(i,j,k)=nh3_dep(i,j,k)+nh3_dep_inc(i,j,k)
        END DO
      END DO
    END DO

    IF (l_nitrate) THEN
      ! Call nitrate chemistry routine
      CALL nitrate(                                                            &
                                ! Arguments IN
        halo_i, halo_j, off_x, off_y,                                          &
        row_length, rows,                                                      &
        model_levels, wet_model_levels,                                        &
        chemstep, cloudf, p_theta_levels, t, q, qcl, qcf,                      &
                                ! Arguments IN/OUT
        hno3, nh3, nitr_acc, nitr_diss,                                        &
                                ! Arguments OUT
        delta_n_chem, delta_n_evap, delta_n_nuc )
      !
    END IF           ! End L_NITRATE condition

  END DO             ! End CALL_CHEM_FREQ loop

END IF               ! End L_SULPC_SO2 test


! OCFF cycle code.
! Calculate emissions of ocff from surface and chimney level, ageing
! of ocff and nucleation scavenging of aged ocff to ocff-in-cloud.

IF (l_ocff) THEN  ! If ocff modelling is included

  IF (l_ocff_surem) THEN  ! If surface ocff emissions are included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      model_levels, wet_model_levels,                                          &
      halo_i, halo_j,                                                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      ocff_new(:,:,1), ocff_surem, 1,                                          &
      timestep, 1, 1, 0.0 )
  END IF  ! L_OCFF_SUREM

  IF (l_ocff_hilem) THEN  ! If chimney level ocff emissions
    ! are included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      model_levels, wet_model_levels,                                          &
      halo_i, halo_j,                                                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      ocff_new(:,:,ocff_high_level), ocff_hilem,                               &
      ocff_high_level, timestep, 1, 1, 0.0 )
  END IF  ! L_OCFF_HILEM

  ! Calculate quantity of fresh ocff converted to aged ocff

  CALL ageocff(                                                                &
    row_length, rows, off_x, off_y,                                            &
    model_levels, timestep,                                                    &
    ocff_new,                                                                  &
    delta_ageocff                                                              &
    )

  ! Calculate quantity of aged ocff scavenged to cloud ocff

  CALL ocffnuclscav(                                                           &
    rows, row_length, off_x, off_y, halo_i, halo_j,                            &
    model_levels, wet_model_levels, timestep,                                  &
    cloudf, qcl, qcf,                                                          &
    ocff_agd, ocff_cld,                                                        &
    delta_ocffnuclscav                                                         &
    )

  ! Update ocff arrays with increments from ageing and
  ! nucleation scavenging

  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        ocff_new(i,j,k)=ocff_new(i,j,k) - delta_ageocff(i,j,k)
        ocff_agd(i,j,k)=ocff_agd(i,j,k) + delta_ageocff(i,j,k)                 &
          - delta_ocffnuclscav(i,j,k)
        ocff_cld(i,j,k)=ocff_cld(i,j,k)                                        &
          + delta_ocffnuclscav(i,j,k)
      END DO
    END DO
  END DO

END IF  ! L_OCFF

! End of ocff cycle code



! PM concentration code: If any of the diagnostics on PM10 and PM2.5
! or total mass concentrations (item numbers 220-235, Sect. 17) is requested
! then call the subroutine that considers all aerosol species and modes

IF (sf(220,17) .OR. sf(221,17) .OR. sf(222,17) .OR.                            &
    sf(223,17) .OR. sf(224,17) .OR. sf(225,17) .OR.                            &
    sf(226,17) .OR. sf(227,17) .OR. sf(228,17) .OR.                            &
    sf(229,17) .OR. sf(230,17) .OR. sf(231,17) .OR.                            &
    sf(232,17) .OR. sf(233,17) .OR. sf(234,17) .OR.                            &
    sf(235,17) .OR. sf(236,17) .OR. sf(237,17) .OR.                            &
    sf(257,17) ) THEN

  CALL calc_pm_diags (                                                         &
    off_x, off_y,                                                              &
    row_length, rows,                                                          &
    model_levels,                                                              &
    salt_dim1, salt_dim2, salt_dim3,                                           &
    l_sulpc_so2, l_soot, l_biomass,                                            &
    l_ocff, l_use_biogenic, l_dust, l_nitrate,                                 &
    l_use_seasalt_pm, p_theta_levels, t,                                       &
    so4_ait, so4_acc, soot_new, soot_agd, bmass_new, bmass_agd,                &
    ocff_new, ocff_agd, biogenic, sea_salt_film, sea_salt_jet,                 &
    dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6,          &
    nitr_acc, pm10, pm2p5,                                                     &
    pm10_so4, pm2p5_so4,                                                       &
    pm10_bc, pm2p5_bc,                                                         &
    pm10_bb, pm2p5_bb,                                                         &
    pm10_ocff, pm2p5_ocff,                                                     &
    pm10_soa, pm2p5_soa,                                                       &
    pm10_ss, pm2p5_ss,                                                         &
    conc_dust, pm10_dust, pm2p5_dust,                                          &
    pm10_nitr, pm2p5_nitr)

END IF


! Call diagnostic routine

CALL diagnostics_aero(                                                         &
  row_length, rows, model_levels,                                              &
  wet_model_levels,                                                            &
  n_rows, global_row_length, global_rows,                                      &
  halo_i, halo_j, off_x, off_y, me,                                            &
  n_proc, n_procx, n_procy,                                                    &
  g_rows, g_row_length,                                                        &
  timestep,                                                                    &
  at_extremity,                                                                &
  l_sulpc_so2, l_dms_em,                                                       &
  l_sulpc_dms, l_sulpc_newdms,                                                 &
  l_sulpc_ozone, l_sulpc_nh3,                                                  &
  l_soot,                                                                      &
  msa, nh3_dep,                                                                &
  dms_em_combined,                                                             &
  deltas_dms,                                                                  &
  f_dms_to_so2,                                                                &
  f_dms_to_so4,                                                                &
  f_dms_to_msa,                                                                &
  deltas_dry,                                                                  &
  deltas_wet,                                                                  &
  deltas_wet_o3,                                                               &
  deltas_evap,                                                                 &
  deltas_nucl,                                                                 &
  deltas_diffuse,                                                              &
  deltas_coag,                                                                 &
  deltas_merge,                                                                &
  delta_n_chem,                                                                &
  delta_n_evap,                                                                &
  delta_n_nuc,                                                                 &
  psi,                                                                         &
  pm10, pm2p5,                                                                 &
  pm10_so4, pm2p5_so4,                                                         &
  pm10_bc, pm2p5_bc,                                                           &
  pm10_bb, pm2p5_bb,                                                           &
  pm10_ocff, pm2p5_ocff,                                                       &
  pm10_soa, pm2p5_soa,                                                         &
  pm10_ss, pm2p5_ss,                                                           &
  conc_dust, pm10_dust, pm2p5_dust,                                            &
  pm10_nitr, pm2p5_nitr,                                                       &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
  stashwork17)
IF (lhook) CALL dr_hook('AERO_CTL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE aero_ctl
END MODULE aero_ctl_mod
