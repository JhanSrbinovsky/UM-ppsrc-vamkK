! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


MODULE cosp_init_mod
  USE cosp_constants_mod, ONLY: n_hydro
  USE cosp_types_mod, ONLY: cosp_config, cosp_gridbox, construct_cosp_gridbox
  USE cosp_input_mod
  USE earth_constants_mod, ONLY: earth_radius
  USE atm_fields_bounds_mod, ONLY: pdims,pdims_s,pdims_l,tdims, &
                                  tdims_s,tdims_l,qdims
  USE ereport_mod
  USE domain_params
  USE Submodel_Mod
  IMPLICIT NONE

! Description:
!   Routine that allocates and initialises the derived types with
!   the gridbox-mean inputs and configuration information for COSP.
!
! Method:
!   Sets the logical flags in cosp_cfg that control the COSP outputs.
!   It also deals with the dependencies between diagnostics and instrument
!   simulators.
!   It allocates the arrays in cosp_gbx and fills them in with relevant model
!   variables.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CONTAINS
  SUBROUTINE cosp_init(L_cosp,L_Rad_Step_prog,L_radiation,row_length,rows,     &
      model_levels,cosp_crain_3d,cosp_csnow_3d,p,q_n,T_n,                      &
      t_surf,p_star,p_theta_levels,r_theta_levels,r_rho_levels,                &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
      L_cosp_call,cosp_npoints,cosp_cfg,cosp_gbx)

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

!-----Input arguments
! Switch that tells if COSP is requested
  LOGICAL,INTENT(IN) :: L_cosp
! Is this a prognostic radiation timestep?
  LOGICAL,INTENT(IN) :: L_Rad_Step_prog
! Is this a radiation timestep
  LOGICAL,INTENT(IN) :: L_radiation
! Domain dimensions
  INTEGER,INTENT(IN) :: row_length,rows,model_levels
! Convective rainfall and snowfall 3D fields
  REAL,INTENT(IN) :: cosp_crain_3d(row_length,rows,model_levels)
  REAL,INTENT(IN) :: cosp_csnow_3d(row_length,rows,model_levels)
! Other model fields needed by COSP
  REAL,INTENT(IN) :: p(pdims_s%i_start:pdims_s%i_end, &
                       pdims_s%j_start:pdims_s%j_end, &
                       pdims_s%k_start:pdims_s%k_end)
  REAL,INTENT(IN) :: q_n(qdims%i_start:qdims%i_end, &
                          qdims%j_start:qdims%j_end,&
                          1:qdims%k_end)
  REAL,INTENT(IN) :: p_theta_levels(tdims_s%i_start:tdims_s%i_end, &
                                    tdims_s%j_start:tdims_s%j_end, &
                                    tdims_s%k_start:tdims_s%k_end)
  REAL,INTENT(IN) :: r_theta_levels(tdims_l%i_start:tdims_l%i_end, &
                                    tdims_l%j_start:tdims_l%j_end, &
                                    tdims_l%k_start:tdims_l%k_end)
  REAL,INTENT(IN) :: r_rho_levels(pdims_l%i_start:pdims_l%i_end, &
                                  pdims_l%j_start:pdims_l%j_end, &
                                  pdims_l%k_start:pdims_l%k_end)
  REAL,INTENT(IN) :: T_n(tdims%i_start:tdims%i_end, &
                         tdims%j_start:tdims%j_end, &
                         1:tdims%k_end)
  REAL,INTENT(IN) :: T_surf(row_length, rows)
  REAL,INTENT(IN) :: p_star(pdims%i_start:pdims%i_end, &
                            pdims%j_start:pdims%j_end)

!-----Output arguments
! COSP needs to be called in this timestep
  LOGICAL,INTENT(OUT) :: L_cosp_call
! Number of horizontal points in the COSP domain
  INTEGER,INTENT(OUT) :: cosp_npoints
! COSP Configuration options
  TYPE(cosp_config),INTENT(OUT)  :: cosp_cfg
! Gridbox information. Input for COSP
  TYPE(cosp_gridbox),INTENT(OUT) :: cosp_gbx

!-----Local variables
  INTEGER :: i, k, icode
  REAL :: cosp_time,cosp_time_bnds(2)
  CHARACTER(LEN=100) :: cmessage ! Error/warning message
  CHARACTER(LEN=9) :: routine_name='COSP_INIT'


   icode = 9
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Set diagnostic flags
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF (L_cosp.AND.L_Rad_Step_prog) THEN
!   Copy instrument flags to cfg structure
    cosp_cfg%Lradar_sim = cosp_cloudsat_sim
    cosp_cfg%Llidar_sim = cosp_lidar_sim
    cosp_cfg%Lisccp_sim = cosp_isccp_sim
    cosp_cfg%Lmisr_sim  = cosp_misr_sim
    cosp_cfg%Lmodis_sim = cosp_modis_sim
    cosp_cfg%Lrttov_sim = cosp_rttov_sim
!   Flag to control output to file. Always false in in-line version
    cosp_cfg%Lwrite_output = .FALSE.
!-------Copy diagnostic flags to cfg structure
!   ISCCP
    cosp_cfg%Lalbisccp       = sf(331,2)
    cosp_cfg%Ltauisccp       = sf(332,2)
    cosp_cfg%Lpctisccp       = sf(333,2)
    cosp_cfg%Lcltisccp       = sf(334,2)
    cosp_cfg%Lmeantbisccp    = sf(335,2)
    cosp_cfg%Lmeantbclrisccp = sf(336,2)
    cosp_cfg%Lclisccp        = sf(337,2)
    cosp_cfg%Lboxptopisccp   = .FALSE.
    cosp_cfg%Lboxtauisccp    = .FALSE.
!   CALIPSO
    cosp_cfg%LlidarBetaMol532 = sf(340,2)
    cosp_cfg%Latb532          = sf(341,2)
    cosp_cfg%Latb532gbx       = (sf(355,2).OR.sf(356,2))
    cosp_cfg%LcfadLidarsr532  = (sf(342,2).OR.sf(370,2))
    cosp_cfg%Lclcalipso       = (sf(343,2).OR.sf(371,2))
    cosp_cfg%Lcllcalipso      = sf(344,2)
    cosp_cfg%Lclmcalipso      = sf(345,2)
    cosp_cfg%Lclhcalipso      = sf(346,2)
    cosp_cfg%Lcltcalipso      = sf(347,2)
    cosp_cfg%LparasolRefl     = sf(348,2)
!   CloudSat
    cosp_cfg%Lcfaddbze94  = (sf(350,2).OR.sf(372,2))
    cosp_cfg%Ldbze94      = sf(351,2)
    cosp_cfg%Ldbze94gbx   = (sf(353,2).OR.sf(354,2))
!   MISR
    cosp_cfg%LclMISR = sf(360,2)
!   MODIS
    cosp_cfg%Lclmodis      = sf(450,2)
    cosp_cfg%Lcltmodis     = sf(451,2)
    cosp_cfg%Lclwmodis     = sf(452,2)
    cosp_cfg%Lclimodis     = sf(453,2)
    cosp_cfg%Lclhmodis     = sf(454,2)
    cosp_cfg%Lclmmodis     = sf(455,2)
    cosp_cfg%Lcllmodis     = sf(456,2)
    cosp_cfg%Ltautmodis    = sf(457,2)
    cosp_cfg%Ltauwmodis    = sf(458,2)
    cosp_cfg%Ltauimodis    = sf(459,2)
    cosp_cfg%Ltautlogmodis = sf(460,2)
    cosp_cfg%Ltauwlogmodis = sf(461,2)
    cosp_cfg%Ltauilogmodis = sf(462,2)
    cosp_cfg%Lreffclwmodis = sf(463,2)
    cosp_cfg%Lreffclimodis = sf(464,2)
    cosp_cfg%Lpctmodis     = sf(465,2)
    cosp_cfg%Llwpmodis     = sf(466,2)
    cosp_cfg%Liwpmodis     = sf(467,2)
!   CloudSat and CALIPSO
    cosp_cfg%Lclcalipso2    = (sf(349,2).OR.sf(374,2))
    cosp_cfg%Lcltlidarradar = .FALSE.
!   RTTOV
    cosp_cfg%Ltbrttov = .FALSE.
!   Other
    cosp_cfg%Lfracout = sf(352,2)

!   Deal with dependencies
    IF (.NOT.cosp_cfg%Lradar_sim) THEN
      cosp_cfg%Lcfaddbze94    = .FALSE.
      cosp_cfg%Lclcalipso2    = .FALSE.
      cosp_cfg%Lcltlidarradar = .FALSE.
      cosp_cfg%Ldbze94        = .FALSE.
      cosp_cfg%Ldbze94gbx     = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Llidar_sim) THEN
      cosp_cfg%Latb532 = .FALSE.
      cosp_cfg%LcfadLidarsr532  = .FALSE.
      cosp_cfg%Lclcalipso2      = .FALSE.
      cosp_cfg%Lclcalipso       = .FALSE.
      cosp_cfg%Lclhcalipso      = .FALSE.
      cosp_cfg%Lcllcalipso      = .FALSE.
      cosp_cfg%Lclmcalipso      = .FALSE.
      cosp_cfg%Lcltcalipso      = .FALSE.
      cosp_cfg%Lcltlidarradar   = .FALSE.
      cosp_cfg%LparasolRefl     = .FALSE.
      cosp_cfg%LlidarBetaMol532 = .FALSE.
      cosp_cfg%Latb532gbx       = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Lisccp_sim) THEN
      cosp_cfg%Lalbisccp       = .FALSE.
      cosp_cfg%Lboxptopisccp   = .FALSE.
      cosp_cfg%Lboxtauisccp    = .FALSE.
      cosp_cfg%Lclisccp        = .FALSE.
      cosp_cfg%Lpctisccp       = .FALSE.
      cosp_cfg%Ltauisccp       = .FALSE.
      cosp_cfg%Lcltisccp       = .FALSE.
      cosp_cfg%Lmeantbisccp    = .FALSE.
      cosp_cfg%Lmeantbclrisccp = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Lmisr_sim) THEN
      cosp_cfg%LclMISR = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Lmodis_sim) THEN
      cosp_cfg%Lclmodis      = .FALSE.
      cosp_cfg%Lcltmodis     = .FALSE.
      cosp_cfg%Lclwmodis     = .FALSE.
      cosp_cfg%Lclimodis     = .FALSE.
      cosp_cfg%Lclhmodis     = .FALSE.
      cosp_cfg%Lclmmodis     = .FALSE.
      cosp_cfg%Lcllmodis     = .FALSE.
      cosp_cfg%Ltautmodis    = .FALSE.
      cosp_cfg%Ltauwmodis    = .FALSE.
      cosp_cfg%Ltauimodis    = .FALSE.
      cosp_cfg%Ltautlogmodis = .FALSE.
      cosp_cfg%Ltauwlogmodis = .FALSE.
      cosp_cfg%Ltauilogmodis = .FALSE.
      cosp_cfg%Lreffclwmodis = .FALSE.
      cosp_cfg%Lreffclimodis = .FALSE.
      cosp_cfg%Lpctmodis     = .FALSE.
      cosp_cfg%Llwpmodis     = .FALSE.
      cosp_cfg%Liwpmodis     = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Lrttov_sim) THEN
      cosp_cfg%Ltbrttov = .FALSE.
    END IF
    IF ((.NOT.cosp_cfg%Lradar_sim).AND.(.NOT.cosp_cfg%Llidar_sim).AND.     &
        (.NOT.cosp_cfg%Lisccp_sim).AND.(.NOT.cosp_cfg%Lmisr_sim).AND.      &
        (.NOT.cosp_cfg%Lmodis_sim)) THEN
      cosp_cfg%Lfracout = .FALSE.
    END IF

!       Are the appropriate masks requested?
    IF (sf(343,2).AND.(.NOT.sf(320,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.320 along 2.343")
    IF (sf(344,2).AND.(.NOT.sf(321,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.321 along 2.344")
    IF (sf(345,2).AND.(.NOT.sf(322,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.322 along 2.345")
    IF (sf(346,2).AND.(.NOT.sf(323,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.323 along 2.346")
    IF (sf(347,2).AND.(.NOT.sf(324,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.324 along 2.347")
    IF (sf(371,2).AND.(.NOT.sf(325,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.325 along 2.371")
    IF ((sf(331,2).OR.sf(332,2).OR.sf(333,2).OR. &
      sf(334,2).OR.sf(337,2).OR.sf(360,2)).AND.(.NOT.sf(330,2))) THEN
      cmessage = " You need stash 2.325 along any of "// &
                 "the following: 2.331-2.334,3.337,2.360"
      CALL ereport(routine_name,icode,cmessage)
    END IF
!   Does the stats routines need to be called?
    cosp_cfg%Lstats = .FALSE.
    IF ((cosp_cfg%Lradar_sim).OR.(cosp_cfg%Llidar_sim).OR. &
      (cosp_cfg%Lisccp_sim)) cosp_cfg%Lstats = .TRUE.
!   Check that vertical grid and diagnostics requested match
    IF (cosp_use_vgrid) THEN
      IF (sf(342,2).OR.sf(343,2).OR.sf(350,2).OR. &
        sf(353,2).OR.sf(355,2)) THEN
        cmessage = " COSP: Stash 2.342,2.343,2.350,2.353,"// &
                   " 2.355 cannot be selected with use_vgrid==.true."
        CALL Ereport(routine_name,icode,cmessage)
      END IF
    ELSE
      IF (sf(354,2).OR.sf(356,2).OR.sf(370,2).OR. &
        sf(371,2).OR.sf(372,2)) THEN
        cmessage = " COSP: Stash 2.354,2.356,2.370,2.371,"// &
                         " 2.372 cannot be selected with use_vgrid==.FALSE."
        CALL Ereport(routine_name,icode,cmessage)
      END IF
    END IF
!   Use of mixing ratios for precip not yet supported in in-line version
    IF (.NOT.cosp_use_precipitation_fluxes) THEN
      cmessage = " COSP_USE_PRECIPITATION_FLUXES must be"// &
                 " .TRUE. in in-line version."
      CALL ereport(routine_name,icode,cmessage)
    END IF
  END IF ! L_cosp.and.L_Rad_Step_prog

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Allocate and fill in input variables on full radiation timesteps
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  L_cosp_call = .FALSE.
  IF (L_radiation.AND.L_Rad_Step_prog.AND.L_cosp) THEN
    L_cosp_call = .TRUE.
    cosp_npoints=row_length*rows
!   cosp_time* are irrelevant in in-line version. Set to sensible values.
    cosp_time=1.0
    cosp_time_bnds(1)=0.0
    cosp_time_bnds(2)=2.0
!-------Allocate memory for gridbox type
    CALL construct_cosp_gridbox(cosp_time, cosp_time_bnds,                     &
          cosp_radar_freq, cosp_surface_radar, cosp_use_mie_tables,            &
          cosp_use_gas_abs, cosp_do_ray, cosp_melt_lay, cosp_k2,               &
          cosp_npoints, model_levels, cosp_Ncolumns, N_HYDRO,                  &
          cosp_Nprmts_max_hydro, cosp_Naero, cosp_Nprmts_max_aero,             &
          cosp_Npoints_it, cosp_lidar_ice_type, cosp_isccp_topheight,          &
          cosp_isccp_topheight_direction, cosp_overlap, cosp_emsfc_lw,         &
          cosp_use_precipitation_fluxes,cosp_use_reff, cosp_platform,          &
          cosp_satellite, cosp_instrument, cosp_Nchannels, cosp_ZenAng,        &
          cosp_channels(1:cosp_nchannels), cosp_surfem(1:cosp_nchannels),      &
          cosp_co2,cosp_ch4,cosp_n2o,cosp_co,cosp_gbx)
!-------Here code to populate input structure
! Following EG dimension changes assume that 
! COSP is not interested in the tdims zeroth level. 
    cosp_gbx%sh        = RESHAPE(q_n,(/cosp_npoints,model_levels/))
    cosp_gbx%p         = RESHAPE(p_theta_levels(1:row_length,1:rows,           &
                                                1:tdims_s%k_end),              &
                                  (/cosp_npoints,model_levels/))
    cosp_gbx%ph        = RESHAPE(p(1:row_length,1:rows,:),                     &
                                  (/cosp_npoints,model_levels/))
    cosp_gbx%psfc      = RESHAPE(p_star,(/cosp_npoints/))
    cosp_gbx%zlev      = RESHAPE(r_theta_levels(1:row_length,                  &
                                  1:rows,1:tdims_l%k_end),                     &
                                  (/cosp_npoints,model_levels/))
    cosp_gbx%zlev      = cosp_gbx%zlev - Earth_radius
    cosp_gbx%zlev_half = RESHAPE(r_rho_levels(1:row_length,1:rows,:),          &
                                  (/cosp_npoints,model_levels/))
    cosp_gbx%zlev_half = cosp_gbx%zlev_half - Earth_radius
    cosp_gbx%T         = RESHAPE(T_n,(/cosp_npoints,model_levels/))
    cosp_gbx%skt       = RESHAPE(t_surf,(/cosp_npoints/))
!-------Compute relative humidity. Using saturated mixing ratio in cosp_gbx%q.
    DO k=1,model_levels
! DEPENDS ON: qsat
      CALL qsat(cosp_gbx%q(:,k),cosp_gbx%T(:,k), &
                cosp_gbx%p(:,k),cosp_npoints)
      DO i=1,cosp_npoints
        IF (cosp_gbx%q(i,k)/=0.0) &
              cosp_gbx%q(i,k) = cosp_gbx%sh(i,k)/cosp_gbx%q(i,k)*100.0
      END DO
    END DO ! k, model levels
!   Sunlit is filled in glue_rad
!   LS rain is filled in mycrophys_ctl. The mass of LS snow is already
!   contained in the LS ice water variable (aggr + crystals +[graupel]),
!   so it does not need to be considered.
    cosp_gbx%rain_cv = RESHAPE(cosp_crain_3d,(/cosp_npoints,model_levels/))
    cosp_gbx%snow_cv = RESHAPE(cosp_csnow_3d,(/cosp_npoints,model_levels/))
  END IF !L_radiation.and.L_Rad_Step_prog.and.L_cosp

  RETURN
  END SUBROUTINE cosp_init
END MODULE cosp_init_mod
