! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE cosp_main_mod
  USE atmos_constants_mod, ONLY: r, repsilon
  USE cosp_constants_mod, ONLY: i_lscliq, i_lscice, i_lsrain, i_lssnow, &
                                i_cvcliq, i_cvcice, i_cvrain, i_cvsnow, &
                                i_lsgrpl, n_hydro,  parasol_nrefl
  USE cosp_types_mod, cosp_rttov_type => cosp_rttov
  USE cosp_input_mod, ONLY: cosp_ncolumns, cosp_csat_vgrid, cosp_nchannels, &
                            cosp_nlr, cosp_use_vgrid, cosp_overlap

  USE cosp_diagnostics_mod

  USE mod_cosp
  USE cosp_reff_mod
  USE ereport_mod
  USE mphys_psd_mod,    ONLY: cr, dr, x2i, ci0, di0, x3i, x4i
  USE mphys_inputs_mod, ONLY: x1r, x2r, ai, bi
  USE mphys_constants_mod, ONLY: x1i, x4r
  IMPLICIT NONE

! Description:
!   Routine that calls the main COSP routine and passes the outputs 
!   to the diagnostics routine/
!
! Method:
!   Control routine that calls the routines that computes effective radii
!   of the requested hydrometeors, allocate the output types, 
!   call the COSP routines, and copy the diagnotics to STASH.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CONTAINS
  SUBROUTINE cosp_main(Ltimer,model_domain,at_extremity,row_length,rows,       &
                  n_rows,model_levels,cosp_cfg,cosp_gbx,                       &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
      STASHwork2)
  
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

!-Input arguments
  LOGICAL,INTENT(IN) :: Ltimer ! If true then output some timing information
  LOGICAL,INTENT(IN) :: at_extremity(4) ! PE at N,S,E or W of the grid
  INTEGER,INTENT(IN) :: model_domain,row_length,rows,n_rows,model_levels
  TYPE(cosp_config),INTENT(IN)     :: cosp_cfg  ! Configuration options
!-Input/output arguments
  TYPE(cosp_gridbox),INTENT(INOUT)  :: cosp_gbx      ! Gridbox-mean inputs
  REAL,INTENT(INOUT)                :: STASHwork2(*) ! STASH workspace
!-Local variables
  TYPE(cosp_subgrid)    :: cosp_sgx     ! Subgrid outputs
  TYPE(cosp_sgradar)    :: cosp_sgrad   ! Output from radar simulator
  TYPE(cosp_sglidar)    :: cosp_sglid   ! Output from lidar simulator
  TYPE(cosp_isccp)      :: cosp_is      ! Output from ISCCP simulator
  TYPE(cosp_misr)       :: cosp_ms      ! Output from MISR simulator
  TYPE(cosp_modis)      :: cosp_mds     ! Output from MODIS simulator
  TYPE(cosp_rttov_type) :: cosp_rttov   ! Output from RTTOV
  TYPE(cosp_vgrid)      :: cosp_vg      ! Vertical grid of stats
  TYPE(cosp_radarstats) :: cosp_stradar ! Summary statistics from radar
  TYPE(cosp_lidarstats) :: cosp_stlidar ! Summary statistics from lidar
  INTEGER :: icode
  REAL, ALLOCATABLE :: aux(:,:)

! Exponent that controls the temperature dependence of the intercept 
! rainfall of the PSD (0.0 -> no dependence)
  REAL, PARAMETER :: X3R = 0.0
! Exponent of the normalised density in the terminal fall speed (UMDP26)
  REAL, PARAMETER :: GX = 0.4
! Coefficients for the density distribution of rainfall (Homogeneous liquid
! spheres)
  REAL, PARAMETER :: AR = 523.6
  REAL, PARAMETER :: BR = 3.0

  LOGICAL :: Lcosp_run
  CHARACTER(LEN=9), PARAMETER :: RoutineName='COSP_MAIN'
  CHARACTER(LEN=200) :: cmessage

!     Check if COSP needs to be run or is configured in diagnostic mode
  IF (cosp_cfg%Lradar_sim .OR. cosp_cfg%Llidar_sim .OR. &
      cosp_cfg%Lisccp_sim .OR. cosp_cfg%Lmodis_sim .OR. &
      cosp_cfg%Lmisr_sim  .OR. cosp_cfg%Lrttov_sim) Lcosp_run = .TRUE.

  IF (cosp_cfg%Lmodis_sim) THEN
     icode = 1
     cmessage = " COSP MODIS simulator not functional at this version."// &
                " Contact the code owner for details."
     CALL Ereport(RoutineName,icode,cmessage)
  END IF

  IF (cosp_cfg%Lrttov_sim) THEN
     icode = 1
     cmessage = " COSP RTTOV simulator not functional at this version."// &
                " Contact the code owner for details."
     CALL Ereport(RoutineName,icode,cmessage)
  END IF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Populate input structure
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute effective radius consistent with microphysics
! The mixing ratios are gathered from several places:
! r2_lwrad for cloud water mixing ratios
! microphys_ctl for LS precip fluxes
! ni_conv_ctl for convective precip
! The LS and CONV liquid effective radius are gathered from r2_lwrad3c
  ALLOCATE(aux(cosp_gbx%npoints,cosp_gbx%nlevels))
! Density
  aux = cosp_gbx%p/(r*cosp_gbx%T*(1.0+((1.0 - repsilon)/repsilon)*cosp_gbx%sh -&
        cosp_gbx%mr_hydro(:,:,I_LSCLIQ) - cosp_gbx%mr_hydro(:,:,I_LSCICE) -    &
        cosp_gbx%mr_hydro(:,:,I_CVCLIQ) - cosp_gbx%mr_hydro(:,:,I_CVCICE)))
! Ice aggregates
  CALL cosp_reff(.FALSE.,X1I,X2I,X3I,X4I,AI,BI,CI0,DI0,GX,                     &
                cosp_gbx%npoints,model_levels,cosp_gbx%T,aux,                  &
                cosp_gbx%mr_hydro(:,:,I_LSCICE),cosp_gbx%Reff(:,:,I_LSCICE))
  CALL cosp_reff(.FALSE.,X1I,X2I,X3I,X4I,AI,BI,CI0,DI0,GX,                     &
                cosp_gbx%npoints,model_levels,cosp_gbx%T,aux,                  &
                cosp_gbx%mr_hydro(:,:,I_CVCICE),cosp_gbx%Reff(:,:,I_CVCICE))
! Rain
  CALL cosp_reff(.TRUE.,X1R,X2R,X3R,X4R,AR,BR,CR,DR,GX,                        &
                cosp_gbx%npoints,model_levels,cosp_gbx%T,aux,                  &
                cosp_gbx%rain_ls,cosp_gbx%Reff(:,:,I_LSRAIN))
  CALL cosp_reff(.TRUE.,X1R,X2R,X3R,X4R,AR,BR,CR,DR,GX,                        &
                    cosp_gbx%npoints,model_levels,cosp_gbx%T,aux,              &
                    cosp_gbx%rain_cv,cosp_gbx%Reff(:,:,I_CVRAIN))
! Convective Snow (LS is included in the ice content)
  CALL cosp_reff(.TRUE.,X1I,X2I,X3I,X4I,AI,BI,CI0,DI0,GX,         &
                    cosp_gbx%npoints,model_levels,cosp_gbx%T,aux, &
                    cosp_gbx%snow_cv,cosp_gbx%Reff(:,:,I_CVSNOW))
  DEALLOCATE(aux)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Define new vertical grid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL construct_cosp_vgrid(cosp_gbx,cosp_nlr,cosp_use_vgrid,                  &
                            cosp_csat_vgrid,cosp_vg)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Allocate memory for other types
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF (Lcosp_run) THEN
    CALL construct_cosp_subgrid(cosp_gbx%npoints, cosp_Ncolumns,               &
                    model_levels, cosp_sgx)
    CALL construct_cosp_sgradar(cosp_cfg,cosp_gbx%npoints,cosp_ncolumns,       &
                    model_levels,N_HYDRO,cosp_sgrad)
    CALL construct_cosp_radarstats(cosp_cfg,cosp_gbx%npoints,cosp_ncolumns,    &
                    cosp_vg%Nlvgrid,N_HYDRO,cosp_stradar)
    CALL construct_cosp_sglidar(cosp_cfg,cosp_gbx%npoints,cosp_ncolumns,       &
                    model_levels,N_HYDRO,PARASOL_NREFL,cosp_sglid)
    CALL construct_cosp_lidarstats(cosp_cfg,cosp_gbx%npoints,cosp_ncolumns,    &
                    cosp_vg%Nlvgrid,N_HYDRO,PARASOL_NREFL,cosp_stlidar)
    CALL construct_cosp_isccp(cosp_cfg,cosp_gbx%npoints,cosp_ncolumns,         &
                    model_levels,cosp_is)
    CALL construct_cosp_misr(cosp_cfg,cosp_gbx%npoints,cosp_ms)
    CALL construct_cosp_modis(cosp_cfg,cosp_gbx%npoints,cosp_mds)
    CALL construct_cosp_rttov(cosp_cfg,cosp_gbx%npoints,cosp_nchannels,        &
                    cosp_rttov)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Call simulator
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     IF (Ltimer) CALL timer ('COSP',5)
     CALL cosp(cosp_overlap,cosp_ncolumns,cosp_cfg,cosp_vg,cosp_gbx,           &
                  cosp_sgx,cosp_sgrad,cosp_sglid,cosp_is,cosp_ms,cosp_mds,     &
                  cosp_rttov,cosp_stradar,cosp_stlidar)
     IF (Ltimer) CALL timer ('COSP',6)
  ELSE ! Minumum allocation
    CALL construct_cosp_subgrid(1, 1, 1, cosp_sgx)
    CALL construct_cosp_sgradar(cosp_cfg,1, 1, 1,N_HYDRO,cosp_sgrad)
    CALL construct_cosp_radarstats(cosp_cfg,1, 1, 1, 1,cosp_stradar)
    CALL construct_cosp_sglidar(cosp_cfg,1, 1, 1, 1, 1,cosp_sglid)
    CALL construct_cosp_lidarstats(cosp_cfg,1, 1, 1, 1, 1,cosp_stlidar)
    CALL construct_cosp_isccp(cosp_cfg,1, 1, 1,cosp_is)
    CALL construct_cosp_misr(cosp_cfg, 1,cosp_ms)
    CALL construct_cosp_modis(cosp_cfg, 1,cosp_mds)
    CALL construct_cosp_rttov(cosp_cfg, 1, 1,cosp_rttov)
  ENDIF ! Lcosp_run

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Copy COSP diagnostics to stash
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL cosp_diagnostics(model_domain,at_extremity,row_length,rows,             &
         n_rows,model_levels,Lcosp_run,cosp_cfg,cosp_vg,cosp_gbx,              &
         cosp_sgx,cosp_sgrad,cosp_sglid,cosp_is,cosp_ms,cosp_mds,cosp_rttov,   &
         cosp_stradar,cosp_stlidar,                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
              STASHwork2)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Deallocate memory in derived types.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL free_cosp_rttov(cosp_rttov)
  CALL free_cosp_modis(cosp_mds)
  CALL free_cosp_misr(cosp_ms)
  CALL free_cosp_isccp(cosp_is)
  CALL free_cosp_lidarstats(cosp_stlidar)
  CALL free_cosp_sglidar(cosp_sglid)
  CALL free_cosp_radarstats(cosp_stradar)
  CALL free_cosp_sgradar(cosp_sgrad)
  CALL free_cosp_subgrid(cosp_sgx)
  CALL free_cosp_vgrid(cosp_vg)

  RETURN
  END SUBROUTINE cosp_main
END MODULE cosp_main_mod
