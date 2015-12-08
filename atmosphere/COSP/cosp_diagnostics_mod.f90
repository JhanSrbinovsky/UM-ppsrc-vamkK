! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE cosp_diagnostics_mod
  USE cosp_constants_mod, ONLY: r_undef, parasol_nrefl, misr_n_cth, sr_bins,   &
                        dbze_bins, n_hydro
  USE cosp_types_mod, ONLY: cosp_config, cosp_vgrid, cosp_gridbox,             &
                        cosp_subgrid, cosp_sgradar, cosp_sglidar, cosp_isccp,  &
                        cosp_misr, cosp_rttov, cosp_radarstats, cosp_lidarstats
  USE cosp_input_mod, ONLY: COSP_NCOLUMNS_MAX
  USE mod_cosp_modis_simulator, only: cosp_modis
  USE UM_ParParams
  USE ereport_mod
  USE domain_params
  USE Submodel_Mod
  IMPLICIT NONE
  PRIVATE create_mask, change_units, undef_to_zero

! Description:
!   Routine that receives the COSP derived types with the outputs and write
!   them to STASH.
!
! Method:
!   Calls to the appropriate STASH routines.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CONTAINS
  SUBROUTINE cosp_diagnostics(model_domain,at_extremity,row_length,rows,       &
      n_rows,model_levels,Lcosp_run,cfg,vgrid,gbx,                             &
      sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar,              &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
      STASHwork)
 
  IMPLICIT NONE
!----Input arguments
  INTEGER,INTENT(IN) :: model_domain
! Indicates if this processor is at north,south, east or west of the PE grid
  LOGICAL,INTENT(IN) :: at_extremity(4)
! Number of points on a row
  INTEGER,INTENT(IN) :: row_length
! Number of rows in a theta field
  INTEGER,INTENT(IN) :: rows
! Number of rows in a v field
  INTEGER,INTENT(IN) :: n_rows
! Number of model levels
  INTEGER,INTENT(IN) :: model_levels
! True if cosp is run
  LOGICAL :: Lcosp_run
! Configuration options
  TYPE(cosp_config),INTENT(IN) :: cfg
! Information on vertical grid of stats
  TYPE(cosp_vgrid),INTENT(IN) :: vgrid
! Gridbox-mean COSP inputs
  TYPE(cosp_gridbox),INTENT(IN) :: gbx
! Subgrid info
  TYPE(cosp_subgrid),INTENT(IN) :: sgx
! Output from radar simulator
  TYPE(cosp_sgradar),INTENT(IN) :: sgradar
! Output from lidar simulator
  TYPE(cosp_sglidar),INTENT(IN) :: sglidar
! Output from ISCCP simulator
  TYPE(cosp_isccp),INTENT(INOUT):: isccp
! Output from MISR simulator
  TYPE(cosp_misr),INTENT(IN)    :: misr
! Output from MODIS simulator
  TYPE(cosp_modis),INTENT(INOUT)    :: modis
! Output from RTTOV
  TYPE(cosp_rttov),INTENT(IN)   :: rttov
! Summary statistics from radar simulator
  TYPE(cosp_radarstats),INTENT(IN) :: stradar
! Summary statistics from lidar simulator
  TYPE(cosp_lidarstats),INTENT(IN) :: stlidar
!----Input/Output arguments
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
! STASH workspace
  REAL,INTENT(INOUT) :: STASHwork(*)

!----Local variables
! Number of pressure levels for ISCCP histograms (2.337)
  INTEGER, PARAMETER :: N_ISCCP_PRESSURE_LEVELS = 7
! Number of tau levels for stash 2.337
  INTEGER, PARAMETER :: N_TAU_LEVELS = 7
! Section 2, LW
  INTEGER,PARAMETER :: sect=2
! Stash code
  INTEGER :: sc_code
! Error status
  INTEGER :: icode
! Cumulative error status
  INTEGER :: errorsum
! Stash-related indices
  INTEGER ::si_tmp,i,pslevel_out,pslevel,si_a,si_b,st_length,il
! Internal model index
  INTEGER :: im_index
! Temporary arrays
  REAL :: aux1D_1(row_length*rows)
  REAL :: aux2D_1(row_length, rows)
  REAL :: aux3D_1(row_length, rows, N_ISCCP_PRESSURE_LEVELS)
  REAL :: aux3D_2(row_length, rows, model_levels)
! Dummy array passed to rad3d_inp, set to 0 because interpolation is not needed.
  REAL :: ES_space_interp(4, row_length, rows)
! Logical flags for tau plevels
  LOGICAL :: l_tau_levels(N_TAU_LEVELS)
! Logical flags for scattering ratio plevels
  LOGICAL :: l_sr_levels(SR_BINS)
! Logical flags for reflectivity plevels
  LOGICAL :: l_dbze_levels(DBZE_BINS)
! Logical flags for hydrometeor plevels
  LOGICAL :: l_hydro_levels(N_HYDRO)
! Logical flags for subcolumns plevels
  LOGICAL :: l_ncol_levels(COSP_NCOLUMNS_MAX)
! Flag to complete field on the northern polar row
  LOGICAL :: L_complete_North
! Flag to complete field on the southern polar row
  LOGICAL :: L_complete_South
! Error message
  CHARACTER(len=80) :: cmessage
! Routine name
  CHARACTER(LEN=16),PARAMETER :: RoutineName='COSP_DIAGNOSTICS'

  icode = 0 ! Initialise error status
  errorsum = 0
  im_index = internal_model_index(atmos_im)

! Initialise local variables
  l_tau_levels(:)   = .FALSE.
  l_sr_levels(:)    = .FALSE.
  l_dbze_levels(:)  = .FALSE.
  l_hydro_levels(:) = .FALSE.
  l_ncol_levels(:)  = .FALSE.
  ES_space_interp = 0.0

! Flags to check whether North/South Poles should be filled in
  L_complete_North=(model_domain==mt_global).AND.(at_extremity(PNorth))
  L_complete_South=(model_domain==mt_global).AND.(at_extremity(PSouth))

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BMDI Masks (Heaviside functions) for CALIPSO diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! 2.320: Mask for 2.343 CALIPSO cloud area on levels
!
  IF (sf(320,2).AND.cfg%Lclcalipso.AND.Lcosp_run) THEN
    sc_code = 320
    aux3D_2 = RESHAPE(stlidar%lidarcld,(/row_length, rows, model_levels/))
    CALL create_mask(R_UNDEF,aux3D_2)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(sc_code,sect,im_index)),aux3D_2,            &
        row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
        stlist(1,stindex(1,sc_code,sect,im_index)),len_stlist,stash_levels,    &
        num_stash_levels+1,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 320)"//cmessage
      GO TO 9999
    END IF
  END IF
!
! 2.321: Mask for 2.344 CALIPSO low-level cloud
!
  IF (sf(321,2).AND.cfg%Lcllcalipso.AND.Lcosp_run) THEN
    sc_code = 321
    aux2D_1 = RESHAPE(stlidar%cldlayer(:,1),(/row_length, rows/))
    CALL create_mask(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 321)"
      GO TO 9999
    END IF
  END IF
!
! 2.322: Mask for 2.345 CALIPSO mid-level cloud
!
  IF (sf(322,2).AND.cfg%Lclmcalipso.AND.Lcosp_run) THEN
    sc_code = 322
    aux2D_1 = RESHAPE(stlidar%cldlayer(:,2),(/row_length, rows/))
    CALL create_mask(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 322)"
      GO TO 9999
    END IF
  END IF
!
! 2.323: Mask for 2.346 CALIPSO high-level cloud
!
  IF (sf(323,2).AND.cfg%Lclmcalipso.AND.Lcosp_run) THEN
    sc_code = 323
    aux2D_1 = RESHAPE(stlidar%cldlayer(:,3),(/row_length, rows/))
    CALL create_mask(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 323)"
      GO TO 9999
    END IF
  END IF
!
! 2.324: Mask for 2.347 CALIPSO total cloud fraction
!
  IF (sf(324,2).AND.cfg%Lcltcalipso.AND.Lcosp_run) THEN
    sc_code = 324
    aux2D_1 = RESHAPE(stlidar%cldlayer(:,4),(/row_length, rows/))
    CALL create_mask(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 324)"
      GO TO 9999
    END IF
  END IF
!
! 2.325: Mask for 2.371 CALIPSO CLOUD AREA ON 40 LEVELS
!
  IF (sf(325,2).AND.sf(371,2).AND.cfg%Lclcalipso.AND.Lcosp_run) THEN
    sc_code = 325
    errorsum = 0
    DO i=1,stlidar%Nlevels
!     Copy diagnostics by looping over levels in call to copydiag.
      si_tmp = si(sc_code,sect,im_index) + (i-1)*row_length*rows
      aux2D_1 = RESHAPE(stlidar%lidarcld(:,stlidar%Nlevels-i+1),               &
                        (/row_length, rows/))
      CALL create_mask(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
      CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,&
                     0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
      CALL copydiag(STASHwork(si_tmp),aux2D_1,row_length,rows,0,0,0,0,         &
                     at_extremity,atmos_im,sect,sc_code,icode,cmessage)
      errorsum = errorsum + icode
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 325)"
      GO TO 9999
    END IF
  END IF


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ISCCP diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The ISCCP cloud fraction (2.334) will be output as fraction, not in %
! as is output by COSP. We do the conversion here because it is used 
! to weight tau, albedo and CTP. Time-averages of these diagnostics
! will then need to be unweighted by the average 2.334 to produce
! mean in-cloud values.
! Hence, isccp%totalcldarea needs to be converted to fraction here.
  IF (cfg%Lcltisccp.AND.Lcosp_run) THEN
      CALL change_units(R_UNDEF,isccp%totalcldarea,0.01) ! % to 1
  END IF

!
! Weights: Stash: sf(330,2)
!
  sc_code = 330
  IF (sf(sc_code,sect).AND.Lcosp_run) THEN
    aux2D_1 = RESHAPE(gbx%sunlit,(/row_length, rows/))
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 330)"
      GO TO 9999
    END IF
  END IF
!
! Weighted Cloud albedo: Stash: sf(331,2)
!
  IF (cfg%Lalbisccp .AND. cfg%Lcltisccp.AND.Lcosp_run) THEN
    sc_code = 331
    aux1D_1 = isccp%totalcldarea*isccp%meanalbedocld
    aux2D_1 = RESHAPE(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
        cmessage="Error in copydiag( item 331)"
        GO TO 9999
    END IF
  END IF
!
! Weighted Cloud optical depth: Stash: sf(332,2)
!
  IF (cfg%Ltauisccp .AND. cfg%Lcltisccp.AND.Lcosp_run) THEN
    sc_code = 332
    aux1D_1 = isccp%totalcldarea*isccp%meantaucld
    aux2D_1 = RESHAPE(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 332)"
      GO TO 9999
    END IF
  END IF
!
! Weighted Mean CTP: Stash: sf(333,2)
!
  IF (cfg%Lpctisccp .AND. cfg%Lcltisccp.AND.Lcosp_run) THEN
    sc_code = 333
    aux1D_1 = isccp%totalcldarea*isccp%meanptop
    aux2D_1 = RESHAPE(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 333)"
      GO TO 9999
    END IF
  END IF
!
! Weighted cloud area: Stash: sf(334,2)
!
  IF (cfg%Lcltisccp.AND.Lcosp_run) THEN
    sc_code = 334
    aux2D_1 = RESHAPE(isccp%totalcldarea,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 334)"
      GO TO 9999
    END IF
  END IF
!
! Cloud brightness temperature: Stash: sf(335,2)
!
  IF (cfg%Lmeantbisccp.AND.Lcosp_run) THEN
    sc_code = 335
    aux2D_1 = RESHAPE(isccp%meantb,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 335)"
      GO TO 9999
    END IF
  END IF
!
! Clear-sky 10.5 brightness temperature: Stash: sf(336,2)
!
  IF (cfg%Lmeantbclrisccp.AND.Lcosp_run) THEN
    sc_code = 336
    aux2D_1 = RESHAPE(isccp%meantbclr,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 336)"
      GO TO 9999
    END IF
  END IF
!
! Cloud top pressure versus optical depth histograms
!
  IF (cfg%Lclisccp.AND.Lcosp_run) THEN
    sc_code = 337
    errorsum = 0
! Vector of logical flags set to .TRUE. for requested plevels
! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(N_TAU_LEVELS,len_stlist,                              &
        stlist(1,stindex(1,sc_code,sect,im_index)),l_tau_levels,               &
        stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
    IF (icode >  0) THEN
      cmessage="Error in set_pseudo_list( item 337)"
      GO TO 9999
    END IF
! Copy the levels requested to STASHwork
    pslevel_out=0
    DO pslevel=1,N_TAU_LEVELS
      IF (l_tau_levels(pslevel)) THEN
        pslevel_out=pslevel_out+1
        DO i=1,N_ISCCP_PRESSURE_LEVELS
! Copy ISCCP diagnostics by looping over 7 levels in call to copydiag.
! This is because copydiag_3d cannot handle ISCCP levels.
          si_tmp = si(sc_code,sect,im_index) +                                 &
                   (pslevel-1)*row_length*rows*N_ISCCP_PRESSURE_LEVELS
          aux2D_1 = RESHAPE(isccp%fq_isccp(:,pslevel,i),(/row_length, rows/))
          CALL change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,     &
                   row_length,rows,0,0,0,0, at_extremity,atmos_im,sect,        &
                   sc_code,icode,cmessage)
          errorsum = errorsum + icode
        END DO
      END IF
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 337)"
      GO TO 9999
    END IF
  END IF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CALIPSO diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! 2.340 CALIPSO MOLECULAR BACKSCATTER
!
  IF (sf(340,2).AND.cfg%LlidarBetaMol532.AND.Lcosp_run) THEN
    sc_code = 340
    aux3D_2 = RESHAPE(sglidar%beta_mol,(/row_length, rows, model_levels/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length,  &
            rows,model_levels,0,0,0,0, at_extremity,                           &
            stlist(1,stindex(1,sc_code,sect,im_index)),len_stlist,stash_levels,&
            num_stash_levels+1,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 340)"//cmessage
    END IF
  END IF

!
! 2.341 Calipso attenuated backscatter
!
  IF (sf(341,2).and.cfg%Latb532) THEN
      IF (sglidar%Ncolumns > COSP_NCOLUMNS_MAX) THEN
        cmessage="Too many subcolumns. Max=30. ( item 341)"
        GOTO 9999
      END IF
      sc_code = 341
      ! Vector of logical flags set to .true. for requested plevels
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(sglidar%Ncolumns,len_stlist,                        &
          stlist(1,stindex(1,sc_code,sect,im_index)),                          &
          l_ncol_levels,stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
      IF (icode >  0) THEN
        cmessage="Error in set_pseudo_list( item 341)"
        GOTO 9999
      END IF
      ! Copy the levels requested to STASHwork
      pslevel_out=0
      DO pslevel=1,sglidar%Ncolumns
      IF (l_ncol_levels(pslevel)) THEN
        pslevel_out=pslevel_out+1
        DO i=1,sglidar%Nlevels
          ! Copy diagnostic by looping over levels in call to copydiag.
          si_tmp = si(sc_code,sect,im_index) +                                 &
                   (pslevel-1)*row_length*rows*sglidar%Nlevels
          aux2D_1 = reshape(sglidar%beta_tot(:,pslevel,i),(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(L_complete_North, L_complete_South,.false.,  &
                          row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
          CALL copydiag (STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,    &
                   row_length,rows,0,0,0,0, at_extremity,atmos_im,sect,sc_code,&
                   icode,cmessage)
          IF (icode  >   0) then
              cmessage="Error in copydiag( item 341)"
              goto 9999
          END IF
        END DO
      END IF
      END DO
  END IF

!
! 2.342 CALIPSO CFAD SCATTERING RATIO
!
  IF (sf(342,2).AND.cfg%LcfadLidarsr532.AND.Lcosp_run) THEN
    sc_code = 342
    errorsum = 0
! Vector of logical flags set to .TRUE. for requested plevels
! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(SR_BINS,len_stlist,                                   &
           stlist(1,stindex(1,sc_code,sect,im_index)),l_sr_levels,             &
           stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
    IF (icode >  0) THEN
      cmessage="Error in set_pseudo_list( item 342)"
      GO TO 9999
    END IF
!   Copy the levels requested to STASHwork
    pslevel_out=0
    DO pslevel=1,SR_BINS
      IF (l_sr_levels(pslevel)) THEN
        pslevel_out=pslevel_out+1
        DO i=1,stlidar%Nlevels
          ! Copy Calipso diagnostics by looping over levels in call to copydiag.
          si_tmp  = si(sc_code,sect,im_index) +                                &
                    (pslevel-1)*row_length*rows*stlidar%Nlevels
          aux2D_1 = RESHAPE(stlidar%cfad_sr(:,pslevel,i),(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,           &
                          row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,     &
                     row_length,rows,0,0,0,0, at_extremity,atmos_im,sect,      &
                     sc_code,icode,cmessage)
          errorsum = errorsum + icode
        END DO
      END IF
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 342)"
      GO TO 9999
    END IF
  END IF
!
! 2.343 CALIPSO CLOUD AREA ON LEVELS
!
  IF (sf(343,2).AND.cfg%Lclcalipso.AND.Lcosp_run) THEN
    sc_code = 343
    aux3D_2 = RESHAPE(stlidar%lidarcld,(/row_length, rows, model_levels/))
    CALL change_units(R_UNDEF,aux3D_2,0.01) ! % to 1
    CALL undef_to_zero(R_UNDEF,aux3D_2)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length, &
         rows,model_levels,0,0,0,0, at_extremity,                              &
         stlist(1,stindex(1,sc_code,sect,im_index)),len_stlist,stash_levels,   &
         num_stash_levels+1,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 343)"//cmessage
      GO TO 9999
    END IF
  END IF
!
! 2.344 CALIPSO LOW-LEVEL CLOUD
!
  IF (cfg%Lcllcalipso.AND.Lcosp_run) THEN
    sc_code = 344
    aux2D_1 = RESHAPE(stlidar%cldlayer(:,1),(/row_length, rows/))
    CALL change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
    CALL undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 344)"
      GO TO 9999
    END IF
  END IF
!
! 2.345 CALIPSO MID-LEVEL CLOUD
!
  IF (cfg%Lclmcalipso.AND.Lcosp_run) THEN
    sc_code = 345
    aux2D_1 = RESHAPE(stlidar%cldlayer(:,2),(/row_length, rows/))
    CALL change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
    CALL undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 345)"
      GO TO 9999
    END IF
  END IF
!
! 2.346 CALIPSO HIGH-LEVEL CLOUD
!
  IF (cfg%Lclhcalipso.AND.Lcosp_run) THEN
    sc_code = 346
    aux2D_1 = RESHAPE(stlidar%cldlayer(:,3),(/row_length, rows/))
    CALL change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
    CALL undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 346)"
      GO TO 9999
    END IF
  END IF
!
! 2.347 CALIPSO TOTAL CLOUD FRACTION
!
  IF (cfg%Lcltcalipso.AND.Lcosp_run) THEN
    sc_code = 347
    aux2D_1 = RESHAPE(stlidar%cldlayer(:,4),(/row_length, rows/))
    CALL change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
    CALL undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
    CALL copydiag(STASHwork(si(sc_code,sect,im_index)),aux2D_1,row_length,rows,&
                  0,0,0,0, at_extremity,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode  >   0) THEN
      cmessage="Error in copydiag( item 347)"
      GO TO 9999
    END IF
  END IF
!
! 2.348 PARASOL TOA REFLECTANCE
!
  IF (cfg%LparasolRefl.AND.Lcosp_run) THEN
    sc_code = 348
    DO i=1,PARASOL_NREFL
      si_tmp = si(sc_code,sect,im_index) + row_length*rows*(i-1)
      aux2D_1 = RESHAPE(stlidar%parasolrefl(:,i),(/row_length, rows/))
! DEPENDS ON: rad3d_inp
      CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,&
                     0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
      CALL copydiag(STASHwork(si_tmp),aux2D_1,row_length,rows,0,0,0,0,         &
                    at_extremity,atmos_im,sect,sc_code,icode,cmessage)
      IF (icode >  0) THEN
        cmessage=": error in copydiag_3d(item 348)"//cmessage
        GO TO 9999
      END IF
    END DO
  END IF
!
! 2.349 CALIPSO CLOUD AREA NOT DETECTED BY CLOUDSAT ON MODEL LEVELS
!
  IF (sf(349,2).AND.cfg%Lclcalipso2.AND.Lcosp_run) THEN
    sc_code = 349
    aux3D_2 = RESHAPE(stradar%lidar_only_freq_cloud,                           &
                      (/row_length, rows, model_levels/))
    CALL change_units(R_UNDEF,aux3D_2,0.01) ! % to 1
    CALL undef_to_zero(R_UNDEF,aux3D_2)
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length,  &
           rows,model_levels,0,0,0,0, at_extremity,                            &
           stlist(1,stindex(1,sc_code,sect,im_index)),len_stlist,stash_levels, &
           num_stash_levels+1,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 349)"//cmessage
      GO TO 9999
    END IF
  END IF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CloudSat diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! 2.350/2.372 CLOUDSAT CFAD ON MODEL LEVELS/40 LEVELS
!
  IF ((sf(350,2).OR.sf(372,2)).AND.cfg%Lcfaddbze94.AND.Lcosp_run) THEN
    IF (sf(350,2)) sc_code = 350
    IF (sf(372,2)) sc_code = 372
    errorsum = 0
! Vector of logical flags set to .TRUE. for requested plevels
! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(DBZE_BINS,len_stlist,                                 &
              stlist(1,stindex(1,sc_code,sect,im_index)),l_dbze_levels,        &
              stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
    IF (icode >  0) THEN
      cmessage="Error in set_pseudo_list( item 350 or 372)"
      GO TO 9999
    END IF
!   Copy the levels requested to STASHwork
    pslevel_out=0
    DO pslevel=1,DBZE_BINS
      IF (l_dbze_levels(pslevel)) THEN
        pslevel_out=pslevel_out+1
        DO i=1,stradar%Nlevels
! Copy CloudSat diagnostics by looping over levels in call to copydiag.
          si_tmp = si(sc_code,sect,im_index) +                                 &
                   (pslevel_out-1)*row_length*rows*stradar%Nlevels
          IF (sf(350,2)) il = i
! NOTE: for stash in cloudsat grid. Levels are written in reverse order.
! I don't understand why, but otherwise the diagnostics are upside down.
          IF (sf(372,2)) il = stradar%Nlevels-i+1
          aux2D_1 = RESHAPE(stradar%cfad_ze(:,pslevel,il),(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,&
                         rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,     &
                        row_length,rows,0,0,0,0, at_extremity,atmos_im,sect,   &
                        sc_code,icode,cmessage)
          errorsum = errorsum + icode
        END DO
      END IF
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 350/372)"
      GO TO 9999
    END IF
  END IF

!
! 2.351 CloudSat reflectivities
!
  IF (sf(351,2).AND.cfg%Ldbze94.AND.Lcosp_run) THEN
    IF (sgradar%Ncolumns > COSP_NCOLUMNS_MAX) THEN
      icode = 9
      cmessage="Too many subcolumns for radar reflectivity diagnostic. "//     &
               "Max=100. ( item 351)"
      GO TO 9999
    END IF
    sc_code = 351
    errorsum = 0
! Vector of logical flags set to .TRUE. for requested plevels
! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(sgradar%Ncolumns,len_stlist,                          &
             stlist(1,stindex(1,sc_code,sect,im_index)),                       &
             l_sr_levels,stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
    IF (icode >  0) THEN
      cmessage="Error in set_pseudo_list( item 351)"
      GO TO 9999
    END IF
! Copy the levels requested to STASHwork
    pslevel_out=0
    DO pslevel=1,sgradar%Ncolumns
      IF (l_sr_levels(pslevel)) THEN
        pslevel_out=pslevel_out+1
        DO i=1,sgradar%Nlevels
! Copy diagnostic by looping over levels in call to copydiag.
          si_tmp = si(sc_code,sect,im_index) +                                 &
                   (pslevel-1)*row_length*rows*sgradar%Nlevels
          aux2D_1 = RESHAPE(sgradar%Ze_tot(:,pslevel,i),(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,&
                         rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,     &
                        row_length,rows,0,0,0,0, at_extremity,atmos_im,sect,   &
                        sc_code,icode,cmessage)
          errorsum = errorsum + icode
        END DO
      END IF
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 351)"
      GO TO 9999
    END IF
  END IF
!
! 2.352 Subcolumns cloud array
!
  IF (sf(352,2).AND.cfg%Lfracout.AND.Lcosp_run) THEN
    IF (sgx%Ncolumns > COSP_NCOLUMNS_MAX) THEN
      icode = 9
      cmessage="Too many subcolumns. Max=100. ( item 352)"
      GO TO 9999
    END IF
    sc_code = 352
    errorsum = 0
! Vector of logical flags set to .TRUE. for requested plevels
! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(sgx%Ncolumns,len_stlist,                              &
          stlist(1,stindex(1,sc_code,sect,im_index)),                          &
          l_sr_levels,stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
    IF (icode >  0) THEN
      cmessage="Error in set_pseudo_list( item 352)"
      GO TO 9999
    END IF
    ! Copy the levels requested to STASHwork
    pslevel_out=0
    DO pslevel=1,sgx%Ncolumns
      IF (l_sr_levels(pslevel)) THEN
        pslevel_out=pslevel_out+1
        DO i=1,sgx%Nlevels
          ! Copy diagnostic by looping over levels in call to copydiag.
          si_tmp = si(sc_code,sect,im_index) +                                 &
                   (pslevel-1)*row_length*rows*sgx%Nlevels
          aux2D_1 = RESHAPE(sgx%frac_out(:,pslevel,i),(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,&
                         rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,     &
                       row_length,rows,0,0,0,0, at_extremity,atmos_im,sect,    &
                       sc_code,icode,cmessage)
          errorsum = errorsum + icode
        END DO
      END IF
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 352)"
      GO TO 9999
    END IF
  END IF
!
! 2.353/2.354 CLOUDSAT GBX MEAN REFLECTIVITY ON MODEL LVLS/40 LVLS
!
  IF ((sf(353,2).OR.sf(354,2)).AND.cfg%Ldbze94gbx.AND.Lcosp_run) THEN
    IF (sf(353,2)) sc_code = 353
    IF (sf(354,2)) sc_code = 354
    errorsum = 0
    DO i=1,stradar%Nlevels
! Copy diagnostics by looping over levels in call to copydiag.
      si_tmp = si(sc_code,sect,im_index) + (i-1)*row_length*rows
      IF (sf(353,2)) il = i
! NOTE: for stash in cloudsat grid. Levels are written in reverse order.
! I don't understand why, but otherwise the diagnostics are upside down.
      IF (sf(354,2)) il = stradar%Nlevels-i+1
      aux2D_1 = RESHAPE(stradar%dbze94gbx(:,il),(/row_length, rows/))
! DEPENDS ON: rad3d_inp
      CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,&
                     0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
      CALL copydiag (STASHwork(si_tmp),aux2D_1,row_length,rows,0,0,0,0,        &
                     at_extremity,atmos_im,sect,sc_code,icode,cmessage)
      errorsum = errorsum + icode
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 353/354)"
      GO TO 9999
    END IF
  END IF

!
! 2.355/2.356 CALIPSO GBX MEAN ATB ON MODEL LVLS/40 LVLS
!
  IF ((sf(355,2).OR.sf(356,2)).AND.cfg%Latb532gbx.AND.Lcosp_run) THEN
    IF (sf(355,2)) sc_code = 355
    IF (sf(356,2)) sc_code = 356
    errorsum = 0
    DO i=1,stlidar%Nlevels
! Copy diagnostics by looping over levels in call to copydiag.
      si_tmp = si(sc_code,sect,im_index) + (i-1)*row_length*rows
      IF (sf(355,2)) il = i
! NOTE: for stash in cloudsat grid. Levels are written in reverse order.
! I don't understand why, but otherwise the diagnostics are upside down.
      IF (sf(356,2)) il = stlidar%Nlevels-i+1
      aux2D_1 = RESHAPE(stlidar%atb532gbx(:,il),(/row_length, rows/))
! DEPENDS ON: rad3d_inp
      CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,&
                     0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
      CALL copydiag(STASHwork(si_tmp),aux2D_1,row_length,rows,0,0,0,0,         &
                    at_extremity,atmos_im,sect,sc_code,icode,cmessage)
      errorsum = errorsum + icode
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 355/356)"
      GO TO 9999
    END IF
  END IF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! MISR diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! 2.360 MISR CTH-TAU HISTOGRAM
!
  IF (cfg%LclMISR.AND.Lcosp_run) THEN
    sc_code = 360
    errorsum = 0
! Vector of logical flags set to .TRUE. for requested plevels
! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(N_TAU_LEVELS,len_stlist,                              &
            stlist(1,stindex(1,sc_code,sect,im_index)),                        &
            l_tau_levels,stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
    IF (icode >  0) THEN
      cmessage="Error in set_pseudo_list( item 360)"
      GO TO 9999
    END IF
! Copy the levels requested to STASHwork
    pslevel_out=0
    DO pslevel=1,N_TAU_LEVELS
      IF (l_tau_levels(pslevel)) THEN
        pslevel_out=pslevel_out+1
        DO i=1,MISR_N_CTH
! Copy MISR diagnostics by looping over 16 levels in call to copydiag.
          si_tmp = si(sc_code,sect,im_index) +                                 &
                   (pslevel-1)*row_length*rows*MISR_N_CTH
          aux2D_1 = RESHAPE(misr%fq_MISR(:,pslevel,misr%Nlevels-i+1),          &
                            (/row_length, rows/))
          CALL change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
          CALL undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length, &
                         rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,     &
                      row_length,rows,0,0,0,0, at_extremity,atmos_im,sect,     &
                      sc_code,icode,cmessage)
          errorsum = errorsum + icode
        END DO
      END IF
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 360)"
      GO TO 9999
    END IF
  END IF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Diagnostics on CloudSat vertical grid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! 2.370 CALIPSO CFAD SCATTERING RATIO 
!
  IF (sf(370,2).AND.cfg%LcfadLidarsr532.AND.Lcosp_run) THEN
    sc_code = 370
    errorsum = 0
! Vector of logical flags set to .TRUE. for requested plevels
! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(SR_BINS,len_stlist,                                   &
            stlist(1,stindex(1,sc_code,sect,im_index)),                        &
            l_sr_levels,stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
    IF (icode >  0) THEN
      cmessage="Error in set_pseudo_list( item 370)"
      GO TO 9999
    END IF
! Copy the levels requested to STASHwork
    pslevel_out=0
    DO pslevel=1,SR_BINS
      IF (l_sr_levels(pslevel)) THEN
        pslevel_out=pslevel_out+1
        DO i=1,stlidar%Nlevels
! Copy Calipso diagnostics by looping over levels in call to copydiag.
          si_tmp = si(sc_code,sect,im_index) +                                 &
          (pslevel_out-1)*row_length*rows*stlidar%Nlevels
          aux2D_1 = RESHAPE(stlidar%cfad_sr(:,pslevel,stlidar%Nlevels-i+1),    &
                            (/row_length, rows/))
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length, &
                        rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,     &
                row_length,rows,0,0,0,0,at_extremity,atmos_im,sect,sc_code,    &
                icode,cmessage)
          errorsum = errorsum + icode
        END DO
      END IF
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 370)"
      GO TO 9999
    END IF
  END IF
!
! 2.371 CALIPSO CLOUD AREA ON 40 LEVELS
!
  IF (sf(371,2).AND.cfg%Lclcalipso.AND.Lcosp_run) THEN
    sc_code = 371
    errorsum = 0
    DO i=1,stlidar%Nlevels
! Copy diagnostics by looping over levels in call to copydiag.
      si_tmp = si(sc_code,sect,im_index) + (i-1)*row_length*rows
      aux2D_1 = RESHAPE(stlidar%lidarcld(:,stlidar%Nlevels-i+1),               &
                        (/row_length, rows/))
      CALL change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
      CALL undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
      CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,     &
                     rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
      CALL copydiag(STASHwork(si_tmp),aux2D_1,row_length,rows,0,0,0,0,         &
                    at_extremity,atmos_im,sect,sc_code,icode,cmessage)
      errorsum = errorsum + icode
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 371)"
      GO TO 9999
    END IF
  END IF
!
! 2.374 CALIPSO CLOUD AREA NOT DETECTED BY CLOUDSAT ON 40 LEVELS
!
  IF (sf(374,2).AND.cfg%Lclcalipso2.AND.Lcosp_run) THEN
    sc_code = 374
    errorsum = 0
    DO i=1,stradar%Nlevels
! Copy diagnostics by looping over levels in call to copydiag.
      si_tmp = si(sc_code,sect,im_index) + (i-1)*row_length*rows
      aux2D_1 = RESHAPE(stradar%lidar_only_freq_cloud(:,stradar%Nlevels-i+1),  &
                        (/row_length, rows/))
      CALL change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
      CALL undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
      CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,&
                     0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
      CALL copydiag(STASHwork(si_tmp),aux2D_1,row_length,rows,0,0,0,0,         &
                    at_extremity,atmos_im,sect,sc_code,icode,cmessage)
      errorsum = errorsum + icode
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 374)"
      GO TO 9999
    END IF
  END IF


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! MODIS diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! The MODIS cloud fractions (2.451-3) will be output as fraction, not in %
     ! as is output by COSP. We do the conversion here because it is used 
     ! to weight tau, albedo and CTP. Time-averages of these diagnostics
     ! will then need to be unweighted by the corresponding cloud fraction
     ! to produce mean in-cloud values.
     ! Hence, we do the conversions to fractions here.
     ! For MODIS, we also replace R_UNDEFs with 0.0.
     ! This is not done for ISCCP because
     ! the missing value in ICARUS is already set to 0.0.
     If (Lcosp_run) then
        If (cfg%Lcltmodis) then
            call change_units(R_UNDEF,modis%Cloud_Fraction_Total_Mean,0.01)
            call undef_to_zero(R_UNDEF,modis%Cloud_Fraction_Total_Mean)
        End if
        If (cfg%Lclwmodis) then
            call change_units(R_UNDEF,modis%Cloud_Fraction_Water_Mean,0.01)
            call undef_to_zero(R_UNDEF,modis%Cloud_Fraction_Water_Mean)
        End if
        If (cfg%Lclimodis) then
            call change_units(R_UNDEF,modis%Cloud_Fraction_Ice_Mean,0.01)
            call undef_to_zero(R_UNDEF,modis%Cloud_Fraction_Ice_Mean)
        End if
        If (cfg%Ltautmodis)                                                    &
           call undef_to_zero(R_UNDEF,modis%Optical_Thickness_Total_Mean)
        If (cfg%Ltauwmodis)                                                    &
           call undef_to_zero(R_UNDEF,modis%Optical_Thickness_Water_Mean)
        If (cfg%Ltauimodis)                                                    &
           call undef_to_zero(R_UNDEF,modis%Optical_Thickness_Ice_Mean)
        If (cfg%Ltautlogmodis)                                                 &
           call undef_to_zero(R_UNDEF,modis%Optical_Thickness_Total_LogMean)
        If (cfg%Ltauwlogmodis)                                                 &
           call undef_to_zero(R_UNDEF,modis%Optical_Thickness_Water_LogMean)
        If (cfg%Ltauilogmodis)                                                 &
           call undef_to_zero(R_UNDEF,modis%Optical_Thickness_Ice_LogMean)
        If (cfg%Lreffclwmodis)                                                 &
           call undef_to_zero(R_UNDEF,modis%Cloud_Particle_Size_Water_Mean)
        If (cfg%Lreffclimodis)                                                 &
           call undef_to_zero(R_UNDEF,modis%Cloud_Particle_Size_Ice_Mean)
        If (cfg%Lpctmodis)                                                     &
           call undef_to_zero(R_UNDEF,modis%Cloud_Top_Pressure_Total_Mean)
        If (cfg%Llwpmodis)                                                     &
           call undef_to_zero(R_UNDEF,modis%Liquid_Water_Path_Mean)
        If (cfg%Liwpmodis)                                                     &
           call undef_to_zero(R_UNDEF,modis%Ice_Water_Path_Mean)
     Endif

!
! 2.450 MODIS CTP-TAU HISTOGRAM
!
     If (sf(450,2).and.cfg%Lclmodis.and.Lcosp_run) then
         sc_code = 450
         ! Vector of logical flags set to .true. for requested plevels
! DEPENDS ON: set_pseudo_list
         Call set_pseudo_list(n_tau_levels,len_stlist,                         &
              stlist(1,stindex(1,sc_code,sect,im_index)),                      &
              l_tau_levels,stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
         If (icode >  0) then
            cmessage="Error in set_pseudo_list( item 450)"
            goto 9999
         End if
         ! Copy the levels requested to STASHwork
         pslevel_out=0
         Do pslevel=1,n_tau_levels
          If (l_tau_levels(pslevel)) then
            pslevel_out=pslevel_out+1
            Do i=1,n_isccp_pressure_levels
              ! Copy diagnostics by looping over 7 levels in call to copydiag.
              ! This is because copydiag_3d cannot handle ISCCP levels.
              si_tmp = si(sc_code,sect,im_index) +                             &
                       (pslevel-1)*row_length*rows*n_isccp_pressure_levels
              aux2D_1 = reshape(                                               &
                  modis%Optical_Thickness_vs_Cloud_Top_Pressure(:,pslevel,i),  &
                  (/row_length, rows/))
              call change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
              call undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(L_complete_North, L_complete_South,.false.,       &
                             row_length,rows,0,0,1,1,1.0,ES_space_interp,      &
                             1,aux2D_1)
! DEPENDS ON: copydiag
              Call copydiag (STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,&
                      row_length,rows,0,0,0,0, at_extremity,atmos_im,sect,     &
                      sc_code,icode,cmessage)
              If (icode  >   0) then
                  cmessage="Error in copydiag( item 450)"
                  goto 9999
              End if
            End do
          End if
         End do
      End if

!
! 2.451 MODIS WEIGHTED TOT. CLOUD AREA
!
     If (sf(451,2).and.cfg%Lcltmodis.and.Lcosp_run) then
         sc_code = 451
         aux2D_1 = reshape(modis%Cloud_Fraction_Total_Mean,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
         Call rad3d_inp(L_complete_North, L_complete_South,.false.,            &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 451)"
            goto 9999
         End if
      End if

!
! 2.452 MODIS LIQUID CLOUD FRACTION 
!
     If (sf(452,2).and.cfg%Lclwmodis.and.Lcosp_run) then
         sc_code = 452
         aux2D_1 = reshape(modis%Cloud_Fraction_Water_Mean,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 452)"
            goto 9999
         End if
      End if

!
! 2.453 MODIS ICE CLOUD FRACTION
!
     If (sf(453,2).and.cfg%Lclimodis.and.Lcosp_run) then
         sc_code = 453
         aux2D_1 = reshape(modis%Cloud_Fraction_Ice_Mean,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 453)"
            goto 9999
         End if
      End if

!
! 2.454 MODIS HIGH-LVL CLOUD FRACTION
!
     If (sf(454,2).and.cfg%Lclhmodis.and.Lcosp_run) then
         sc_code = 454
         aux2D_1 = reshape(modis%Cloud_Fraction_High_Mean,(/row_length, rows/))
         call change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
         call undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,  &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 454)"
            goto 9999
         End if
      End if

!
! 2.455 MODIS MID-LVL CLOUD FRACTION
!
     If (sf(455,2).and.cfg%Lclmmodis.and.Lcosp_run) then
         sc_code = 455
         aux2D_1 = reshape(modis%Cloud_Fraction_Mid_Mean,(/row_length, rows/))
         call change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
         call undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 455)"
            goto 9999
         End if
      End if

!
! 2.456 MODIS LOW-LVL CLOUD FRACTION
!
     If (sf(456,2).and.cfg%Lcllmodis.and.Lcosp_run) then
         sc_code = 456
         aux2D_1 = reshape(modis%Cloud_Fraction_Low_Mean,(/row_length, rows/))
         call change_units(R_UNDEF,aux2D_1,0.01) ! % to 1
         call undef_to_zero(R_UNDEF,aux2D_1)
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
            icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 456)"
            goto 9999
         End if
      End if

!
! 2.457 MODIS WEIGHTED TAU - TOTAL
!
     If (cfg%Ltautmodis .AND. cfg%Lcltmodis.and.Lcosp_run) then
         sc_code = 457
         aux1D_1 = modis%Cloud_Fraction_Total_Mean *                           &
                   modis%Optical_Thickness_Total_Mean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 457)"
            goto 9999
         End if
      End if

!
! 2.458 MODIS WEIGHTED TAU - LIQUID
!
     If (cfg%Ltauwmodis .AND. cfg%Lclwmodis.and.Lcosp_run) then
         sc_code = 458
         aux1D_1 = modis%Cloud_Fraction_Water_Mean *                           &
                   modis%Optical_Thickness_Water_Mean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 458)"
            goto 9999
         End if
      End if

!
! 2.459 MODIS WEIGHTED TAU - ICE
!
     If (cfg%Ltauimodis .AND. cfg%Lclimodis.and.Lcosp_run) then
         sc_code = 459
         aux1D_1 = modis%Cloud_Fraction_Ice_Mean *                             &
                   modis%Optical_Thickness_Ice_Mean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 459)"
            goto 9999
         End if
      End if

!
! 2.460 MODIS WEIGHTED LOG(TAU) -TOTAL
!
     If (cfg%Ltautlogmodis .AND. cfg%Lcltmodis.and.Lcosp_run) then
         sc_code = 460
         aux1D_1 = modis%Cloud_Fraction_Total_Mean *                           &
                   modis%Optical_Thickness_Total_LogMean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 460)"
            goto 9999
         End if
      End if

!
! 2.461 MODIS WEIGHTED LOG(TAU) - LIQ
!
     If (cfg%Ltauwlogmodis .AND. cfg%Lclwmodis.and.Lcosp_run) then
         sc_code = 461
         aux1D_1 = modis%Cloud_Fraction_Water_Mean *                           &
                   modis%Optical_Thickness_Water_LogMean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 461)"
            goto 9999
         End if
      End if

!
! 2.462 MODIS WEIGHTED LOG(TAU) - ICE
!
     If (cfg%Ltauilogmodis .AND. cfg%Lclimodis.and.Lcosp_run) then
         sc_code = 462
         aux1D_1 = modis%Cloud_Fraction_Ice_Mean *                             &
                   modis%Optical_Thickness_Ice_LogMean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 462)"
            goto 9999
         End if
      End if

!
! 2.463 MODIS WEIGHTED Reff LIQUID
!
     If (cfg%Lreffclwmodis .AND. cfg%Lclwmodis.and.Lcosp_run) then
         sc_code = 463
         aux1D_1 = modis%Cloud_Fraction_Water_Mean *                           &
                   modis%Cloud_Particle_Size_Water_Mean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 463)"
            goto 9999
         End if
      End if

!
! 2.464 MODIS WEIGHTED Reff ICE
!
     If (cfg%Lreffclimodis .AND. cfg%Lclimodis.and.Lcosp_run) then
         sc_code = 464
         aux1D_1 = modis%Cloud_Fraction_Ice_Mean *                             &
                   modis%Cloud_Particle_Size_Ice_Mean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 464)"
            goto 9999
         End if
      End if

!
! 2.465 MODIS WEIGHTED CLOUD TOP PRES.
!
     If (cfg%Lpctmodis .AND. cfg%Lcltmodis.and.Lcosp_run) then
         sc_code = 465
         aux1D_1 = modis%Cloud_Fraction_Total_Mean *                           &
                   modis%Cloud_Top_Pressure_Total_Mean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 465)"
            goto 9999
         End if
      End if

!
! 2.466 MODIS WEIGHTED LIQUID W. PATH
!
     If (cfg%Llwpmodis .AND. cfg%Lclwmodis.and.Lcosp_run) then
         sc_code = 466
         aux1D_1 = modis%Cloud_Fraction_Water_Mean*modis%Liquid_Water_Path_Mean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 466)"
            goto 9999
         End if
      End if

!
! 2.467 MODIS WEIGHTED ICE WATER PATH
!
     If (cfg%Liwpmodis .AND. cfg%Lclimodis.and.Lcosp_run) then
         sc_code = 467
         aux1D_1 = modis%Cloud_Fraction_Ice_Mean*modis%Ice_Water_Path_Mean
         aux2D_1 = reshape(aux1D_1,(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(L_complete_North, L_complete_South,.false.,           &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(sc_code,sect,im_index)),                  &
             aux2D_1,row_length,rows,0,0,0,0, at_extremity,                    &
             atmos_im,sect,sc_code,                                            &
             icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 467)"
            goto 9999
         End if
      End if


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! COSP inputs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! 2.373 HYDROMETEOR EFFECTIVE RADIUS
!
  IF (sf(373,2)) THEN
    sc_code = 373
    errorsum = 0
! Vector of logical flags set to .TRUE. for requested plevels
! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(N_HYDRO,len_stlist,                                   &
            stlist(1,stindex(1,sc_code,sect,im_index)),                        &
            l_hydro_levels,stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
    IF (icode >  0) THEN
      cmessage="Error in set_pseudo_list( item 373)"
      GO TO 9999
    END IF
!   Copy the levels requested to STASHwork
    pslevel_out=0
    DO pslevel=1,N_HYDRO
      IF (l_hydro_levels(pslevel)) THEN
        pslevel_out=pslevel_out+1
        DO i=1,gbx%Nlevels
!         Copy CloudSat diagnostics by looping over levels in call to copydiag.
          si_tmp = si(sc_code,sect,im_index) +                                 &
                   (pslevel-1)*row_length*rows*gbx%Nlevels
          aux2D_1 = RESHAPE(gbx%reff(:,i,pslevel),(/row_length, rows/))
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,  &
                         row_length,rows,0,0,1,1,1.0,ES_space_interp,1,aux2D_1)
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si_tmp+(row_length*rows*(i-1))),aux2D_1,     &
                     row_length,rows,0,0,0,0, at_extremity,atmos_im,sect,      &
                     sc_code,icode,cmessage)
          errorsum = errorsum + icode
        END DO
      END IF
    END DO
    IF (errorsum  >   0) THEN
      icode = errorsum
      cmessage="Error in copydiag( item 373)"
      GO TO 9999
    END IF
  END IF
!
! 2.375 LARGE-SCALE CLOUD OPT. DEPTH
!
  IF (sf(375,2)) THEN
    sc_code = 375
    aux3D_2 = RESHAPE(gbx%dtau_s,(/row_length, rows, model_levels/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North, L_complete_South,.FALSE.,row_length,rows, &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length,  &
            rows,model_levels,0,0,0,0, at_extremity,                           &
            stlist(1,stindex(1,sc_code,sect,im_index)),                        &
            len_stlist,stash_levels,num_stash_levels+1,atmos_im,sect,sc_code,  &
            icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 375)"//cmessage
      GO TO 9999
    END IF
  END IF
!
! 2.376 LARGE-SCALE CLOUD EMISSIVITY
!
  IF (sf(376,2)) THEN
    sc_code = 376
    aux3D_2 = RESHAPE(gbx%dem_s,(/row_length, rows, model_levels/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,  &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length,  &
         rows,model_levels,0,0,0,0, at_extremity,                              &
         stlist(1,stindex(1,sc_code,sect,im_index)),len_stlist,stash_levels,   &
         num_stash_levels+1,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 376)"//cmessage
      GO TO 9999
    END IF
  END IF
!
! 2.377 CONVECTIVE CLOUD OPT. DEPTH
!
  IF (sf(377,2)) THEN
    sc_code = 377
    aux3D_2 = RESHAPE(gbx%dtau_c,(/row_length, rows, model_levels/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,  &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length,  &
           rows,model_levels,0,0,0,0, at_extremity,                            &
           stlist(1,stindex(1,sc_code,sect,im_index)),len_stlist,stash_levels, &
           num_stash_levels+1,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 377)"//cmessage
      GO TO 9999
    END IF
  END IF
!
! 2.378 CONVECTIVE CLOUD EMISSIVITY
!
  IF (sf(378,2)) THEN
    sc_code = 378
    aux3D_2 = RESHAPE(gbx%dem_c,(/row_length, rows, model_levels/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,  &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length,  &
            rows,model_levels,0,0,0,0, at_extremity,                           &
            stlist(1,stindex(1,sc_code,sect,im_index)),len_stlist,stash_levels,&
            num_stash_levels+1,atmos_im,sect,sc_code,icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 378)"//cmessage
      GO TO 9999
    END IF
  END IF
!
! 2.380-388 HYDROMETEOR EFFECTIVE RADIUS 
!
  DO pslevel=1,N_HYDRO
    IF (sf(379+pslevel,2)) THEN
      sc_code = 379+pslevel
      aux3D_2 = RESHAPE(gbx%reff(:,:,pslevel),                                 &
                        (/row_length, rows, model_levels/))
! DEPENDS ON: rad3d_inp
      CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,&
                     0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
      CALL copydiag_3d(stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length,&
              rows,model_levels,0,0,0,0, at_extremity,                         &
              stlist(1,stindex(1,sc_code,sect,im_index)),                      &
              len_stlist,stash_levels,num_stash_levels+1,atmos_im,sect,        &
              sc_code,icode,cmessage)
      IF (icode >  0) THEN
        cmessage=": error in copydiag_3d(item 380-88)"//cmessage
        GO TO 9999
      END IF
    END IF
  END DO
!
! 2.389 COSP: 3D CONVECTIVE RAINFALL RATE
!
  IF (sf(389,2)) THEN
    sc_code = 389
    aux3D_2 = RESHAPE(gbx%rain_cv,(/row_length, rows, model_levels/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,  &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length,  &
            rows,model_levels,0,0,0,0, at_extremity,                           &
            stlist(1,stindex(1,sc_code,sect,im_index)),                        &
            len_stlist,stash_levels,num_stash_levels+1,atmos_im,sect,sc_code,  &
            icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 389)"//cmessage
      GO TO 9999
    END IF
  END IF
!
! 2.390 COSP: 3D CONVECTIVE SNOWFALL RATE
!
  IF (sf(390,2)) THEN
    sc_code = 390
    aux3D_2 = RESHAPE(gbx%snow_cv,(/row_length, rows, model_levels/))
! DEPENDS ON: rad3d_inp
    CALL rad3d_inp(L_complete_North,L_complete_South,.FALSE.,row_length,rows,  &
                   0,0,1,1,1.0,ES_space_interp,model_levels,aux3D_2)
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(sc_code,sect,im_index)),aux3D_2,row_length,  &
            rows,model_levels,0,0,0,0, at_extremity,                           &
            stlist(1,stindex(1,sc_code,sect,im_index)),                        &
            len_stlist,stash_levels,num_stash_levels+1,atmos_im,sect,sc_code,  &
            icode,cmessage)
    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 390)"//cmessage
      GO TO 9999
    END IF
  END IF

 9999 CONTINUE  ! exit point on error
  IF(icode /= 0) THEN
    CALL Ereport(RoutineName,icode,cmessage)
  END IF

  RETURN
  END SUBROUTINE cosp_diagnostics

! Subroutine Interface:
  ELEMENTAL SUBROUTINE create_mask(mdi,x)
    IMPLICIT NONE
    REAL,INTENT(IN) :: mdi
    REAL,INTENT(INOUT) :: x

    IF (x == mdi) THEN
      x = 0.0
    else
      x = 1.0
    END IF
  END SUBROUTINE create_mask
! Subroutine Interface:
  ELEMENTAL SUBROUTINE change_units(mdi,x,factor)
    IMPLICIT NONE
    REAL,INTENT(IN) :: mdi,factor
    REAL,INTENT(INOUT) :: x

    IF (x /= mdi) THEN
      x = factor*x
    END IF
  END SUBROUTINE change_units
! Subroutine Interface:
  ELEMENTAL SUBROUTINE undef_to_zero(mdi,x)
    IMPLICIT NONE
    REAL,INTENT(IN) :: mdi
    REAL,INTENT(INOUT) :: x

    IF (x == mdi) THEN
      x = 0.0
    END IF
  END SUBROUTINE undef_to_zero
END MODULE cosp_diagnostics_mod
