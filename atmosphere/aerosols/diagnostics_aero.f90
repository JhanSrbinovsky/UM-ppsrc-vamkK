! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE diagnostics_aero_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE diagnostics_aero(                             &
  row_length, rows, model_levels,                        &
  wet_model_levels,                                      &
  n_rows, global_row_length, global_rows,                &
  halo_i, halo_j, off_x, off_y, me,                      &
  n_proc, n_procx, n_procy,                              &
  g_rows, g_row_length,                                  &
  timestep,                                              &
  at_extremity,                                          &
  l_sulpc_so2, l_dms_em,                                 &
  l_sulpc_dms, l_sulpc_newdms,                           &
  l_sulpc_ozone, l_sulpc_nh3,                            &
  l_soot,                                                &
  msa, nh3_dep,                                          &
  dms_emiss,                                             &
  deltas_dms,                                            &
  f_dms_to_so2,                                          &
  f_dms_to_so4,                                          &
  f_dms_to_msa,                                          &
  deltas_dry,                                            &
  deltas_wet,                                            &
  deltas_wet_o3,                                         &
  deltas_evap,                                           &
  deltas_nucl,                                           &
  deltas_diffuse,                                        &
  deltas_coag,                                           &
  deltas_merge,                                          &
  delta_n_chem,                                          &
  delta_n_evap,                                          &
  delta_n_nuc,                                           &
  psi,                                                   &
  pm10,      pm2p5,                                      &
  pm10_so4,  pm2p5_so4,                                  &
  pm10_bc,   pm2p5_bc,                                   &
  pm10_bb,   pm2p5_bb,                                   &
  pm10_ocff, pm2p5_ocff,                                 &
  pm10_soa,  pm2p5_soa,                                  &
  pm10_ss,   pm2p5_ss,                                   &
  conc_dust, pm10_dust, pm2p5_dust,                      &
  pm10_nitr, pm2p5_nitr,                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
  stashwork )

!----------------------------------------------------------------------
! Purpose:  Calculates diagnostics from section 17 AERO_CTL2 routine
!           and outputs them.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols

! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards

! Documentation:  UMDP 20

!-----------------------------------------------------------------------


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
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

!
! Arguments with Intent IN. ie: Input variables.

LOGICAL ::  at_extremity(4)  ! Indicates if this processor is at north,
                             ! south, east or west of the processor grid

! Parameters

INTEGER :: row_length           ! number of points on a row
INTEGER :: rows                 ! number of rows in a theta field
INTEGER :: n_rows               ! number of rows in a v field
INTEGER :: model_levels         ! number of model levels
INTEGER :: wet_model_levels     ! number of model levels where moisture
INTEGER :: global_row_length    ! NUMBER OF points on a global row
INTEGER :: global_rows          ! NUMBER OF global rows
INTEGER :: me                   ! Processor number
INTEGER :: halo_i               ! size of large halo in x direction
INTEGER :: halo_j               ! size of large halo in y direction
INTEGER :: off_x                ! size of small halo in x direction
INTEGER :: off_y                ! size of small halo in y direction
INTEGER :: n_proc
INTEGER :: n_procx
INTEGER :: n_procy
INTEGER :: g_rows (0:n_proc-1)
INTEGER :: g_row_length (0:n_proc-1)

REAL    :: timestep

LOGICAL :: l_sulpc_so2          ! T if S Cycle on
LOGICAL :: l_sulpc_dms          ! T if DMS included
LOGICAL :: l_dms_em             ! T if DMS emissions used
LOGICAL :: l_sulpc_newdms       ! T if new DMS scheme used (requires OZONE)
LOGICAL :: l_sulpc_ozone        ! T if OZONE field present
LOGICAL :: l_sulpc_nh3          ! T if NH3 field present
LOGICAL :: l_soot               ! T if SOOT modelling on

! Arguments with intent IN/OUT (diagnostics):
REAL    :: msa(row_length,rows,model_levels)             ! mmr S in MSA
REAL    :: nh3_dep(row_length,rows,model_levels)         ! NH3 depleted
REAL    :: dms_emiss(row_length,rows)                    ! DMS emiss (kgSm-2s-1)
REAL    :: deltas_dms(row_length,rows,model_levels)
REAL    :: f_dms_to_so2(row_length,rows,model_levels)
REAL    :: f_dms_to_so4(row_length,rows,model_levels)
REAL    :: f_dms_to_msa(row_length,rows,model_levels)
REAL    :: deltas_dry(row_length,rows,model_levels)
REAL    :: deltas_wet(row_length,rows,model_levels)
REAL    :: deltas_wet_o3(row_length,rows,model_levels)
REAL    :: deltas_evap(row_length,rows,model_levels)
REAL    :: deltas_nucl(row_length,rows,model_levels)
REAL    :: deltas_diffuse(row_length,rows,model_levels)
REAL    :: deltas_coag(row_length,rows,model_levels)
REAL    :: deltas_merge(row_length,rows,model_levels)
REAL    :: delta_n_chem(row_length,rows,model_levels)
REAL    :: delta_n_evap(row_length,rows,model_levels)
REAL    :: delta_n_nuc(row_length,rows,model_levels)
REAL    :: psi(row_length,rows,model_levels)
REAL    :: pm10(row_length,rows,model_levels)            ! PM10 (ug m-3)
REAL    :: pm2p5(row_length,rows,model_levels)           ! PM2.5 (ug m-3)
! Sulphate contributions to PM concs.
REAL    :: pm10_so4 (row_length, rows, model_levels)
REAL    :: pm2p5_so4(row_length, rows, model_levels)
! Black carbon contrib. to PM concs.
REAL    :: pm10_bc (row_length, rows, model_levels)
REAL    :: pm2p5_bc(row_length, rows, model_levels)
! Biomass aerosol contrib to PM concs.
REAL    :: pm10_bb (row_length, rows, model_levels)
REAL    :: pm2p5_bb(row_length, rows, model_levels)
! OCFF contributions to PM concs.
REAL    :: pm10_ocff (row_length, rows, model_levels)
REAL    :: pm2p5_ocff(row_length, rows, model_levels)
! SOA contributions to PM concs.
REAL    :: pm10_soa (row_length, rows, model_levels)
REAL    :: pm2p5_soa(row_length, rows, model_levels)
! Sea-salt contributions to PM concs.
REAL    :: pm10_ss (row_length, rows, model_levels)
REAL    :: pm2p5_ss(row_length, rows, model_levels)
! Dust contributions to PM concs.
REAL    :: conc_dust (row_length, rows, model_levels)  ! Total dust concentration
REAL    :: pm10_dust (row_length, rows, model_levels)
REAL    :: pm2p5_dust(row_length, rows, model_levels)
! Nitrate contributions to PM concs.
REAL    :: pm10_nitr (row_length, rows, model_levels)
REAL    :: pm2p5_nitr(row_length, rows, model_levels)

! Diagnostic variables
REAL    :: stashwork(*)         ! STASH workspace for section 17
! Local variables
INTEGER ::  i, j, level, k      ! loop counters
INTEGER ::  icode               ! Return code  =0 Normal exit  >1 Error

! STASH section, item no.s
INTEGER,PARAMETER :: sect = 17 !  for aero_ctl (S Cycle or soot)
INTEGER           :: item

REAL :: work_3d(row_length,rows,model_levels) ! work space

CHARACTER (LEN=80) :: cmessage
CHARACTER (LEN=* ), PARAMETER ::  routinename='diagnostics_aero'

INTEGER :: im_index        ! internal model index

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(routinename,zhook_in,zhook_handle)
icode = 0 ! Initialise error status
im_index = internal_model_index(atmos_im)

! Copy diagnostic information to STASHwork for STASH processing

! Write MSA to STASH if DMS included

IF (l_sulpc_dms) THEN

  item = 203                          !MSA
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! Convert to flux per sec
    DO level=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          work_3d(i,j,level) = msa(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                       &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0, at_extremity,                      &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 203)"//cmessage
    END IF

  END IF

END IF


! Write NH3_DEP to STASH if NH3 included

IF (l_sulpc_nh3) THEN

  item = 204                          !NH3_DEP
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! Convert to flux per sec
    DO level=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          work_3d(i,j,level) = nh3_dep(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                       &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0, at_extremity,                      &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 204)"//cmessage
    END IF

  END IF

END IF         ! End L_SULPC_NH3 condn



! Diagnose DMS emissions if requested

IF (l_dms_em) THEN

  item = 205                          !DMS emissions
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! DEPENDS ON: copydiag
    CALL copydiag (stashwork(si(item,sect,im_index)),                          &
      dms_emiss,                                                               &
      row_length,rows,0,0,0,0, at_extremity,                                   &
      atmos_im,sect,item,                                                      &
      icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": Error in copydiag (item 205)"//cmessage
    END IF

  END IF

END IF         ! End L_DMS_em condn

IF (l_sulpc_so2 .AND. l_sulpc_dms) THEN

  item = 206
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_dms(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 206)"//cmessage
    END IF

  END IF

  item = 207
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! get the fraction and convert to flux per second
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_dms(i,j,level) *                         &
            f_dms_to_so2(i,j,level) / timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 207)"//cmessage
    END IF

  END IF

  item = 208
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! get the fraction and convert to flux per second
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_dms(i,j,level) *                         &
            f_dms_to_so4(i,j,level) *                                          &
            (1.0e00 - psi(i,j,level))/ timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 208)"//cmessage
    END IF

  END IF

  item = 209
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! get the fraction and convert to flux per second
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_dms(i,j,level) *                         &
            f_dms_to_so4(i,j,level) *                                          &
            psi(i,j,level) / timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 209)"//cmessage
    END IF

  END IF

  item = 210
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_dry(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 210)"//cmessage
    END IF

  END IF

  item = 211
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_wet(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 211)"//cmessage
    END IF

  END IF

  item = 212
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_wet_o3(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 212)"//cmessage
    END IF

  END IF

  item = 213
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_evap(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 213)"//cmessage
    END IF

  END IF

  item = 214
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_nucl(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 214)"//cmessage
    END IF

  END IF

  item = 215
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_diffuse(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 215)"//cmessage
    END IF

  END IF

  item = 216
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_coag(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 216)"//cmessage
    END IF

  END IF

  item = 217
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! get the fraction and convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_dry(i,j,level)*                          &
            (1.0e00-psi(i,j,level))/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 217)"//cmessage
    END IF

  END IF

  item = 218
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! get the fraction and convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_dry(i,j,level)*                          &
            psi(i,j,level)/timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 218)"//cmessage
    END IF

  END IF

  item = 219
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    ! convert to flux per sec
    DO level=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,level) = deltas_merge(i,j,level)                         &
            /timestep
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0,at_extremity,                       &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                      &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,sect,item,                                                      &
      icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 219)"//cmessage
    END IF

  END IF


END IF

! ---------------------------------------------------------------------------

!     Diagnose PM10, PM2.5 and the contributions of the different
!     aerosols species to them if requested. Note that PM10 & PM2.5
!     can be calculated as long as any of the aerosol species is used.

! Diagnose PM10 if requested:

item = 220
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm10,                                                                      &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 220)"//cmessage
  END IF
END IF

! Diagnose PM2.5 if requested:

item = 221
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm2p5,                                                                     &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 221)"//cmessage
  END IF
END IF

! Diagnose PM10 & PM2.5 concs. due to different aerosol species
! if requested:

item = 222
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm10_so4,                                                                  &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 222)"//cmessage
  END IF
END IF

item = 223
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm2p5_so4,                                                                 &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 223)"//cmessage
  END IF
END IF

item = 224
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm10_bc,                                                                   &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 224)"//cmessage
  END IF
END IF

item = 225
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm2p5_bc,                                                                  &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 225)"//cmessage
  END IF
END IF

item = 226
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm10_bb,                                                                   &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 226)"//cmessage
  END IF
END IF

item = 227
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm2p5_bb,                                                                  &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 227)"//cmessage
  END IF
END IF

item = 228
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm10_ocff,                                                                 &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 228)"//cmessage
  END IF
END IF

item = 229
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm2p5_ocff,                                                                &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 229)"//cmessage
  END IF
END IF

item = 230
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm10_soa,                                                                  &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 230)"//cmessage
  END IF
END IF

item = 231
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm2p5_soa,                                                                 &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 231)"//cmessage
  END IF
END IF

item = 232
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm10_ss,                                                                   &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 232)"//cmessage
  END IF
END IF

item = 233
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm2p5_ss,                                                                  &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 233)"//cmessage
  END IF
END IF

item = 234
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm10_dust,                                                                 &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 234)"//cmessage
  END IF
END IF

item = 235
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm2p5_dust,                                                                &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 235)"//cmessage
  END IF
END IF

item = 236
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm10_nitr,                                                                 &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 236)"//cmessage
  END IF
END IF

item = 237
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                         &
    pm2p5_nitr,                                                                &
    row_length,rows,model_levels,0,0,0,0, at_extremity,                        &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 237)"//cmessage
  END IF
END IF
!
! Fluxes associated with nitrate chemistry:
!
item = 240
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! convert to flux per sec
  DO level=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d(i,j,level) = delta_n_chem(i,j,level)                           &
          /timestep
      END DO
    END DO
  END DO
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,sect,im_index)),                          &
    work_3d,                                                                   &
    row_length,rows,model_levels,0,0,0,0,at_extremity,                         &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item 240)"//cmessage
  END IF
END IF
!
item = 241
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! convert to flux per sec
  DO level=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d(i,j,level) = delta_n_nuc(i,j,level)                            &
          /timestep
      END DO
    END DO
  END DO
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,sect,im_index)),                          &
    work_3d,                                                                   &
    row_length,rows,model_levels,0,0,0,0,at_extremity,                         &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item 241)"//cmessage
  END IF
END IF
!
item = 242
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! convert to flux per sec
  DO level=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d(i,j,level) = delta_n_evap(i,j,level)                           &
          /timestep
      END DO
    END DO
  END DO
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,sect,im_index)),                          &
    work_3d,                                                                   &
    row_length,rows,model_levels,0,0,0,0,at_extremity,                         &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,                        &
    stash_levels,num_stash_levels+1,                                           &
    atmos_im,sect,item,                                                        &
    icode,cmessage)
  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item 242)"//cmessage
  END IF
END IF
! Diagnose total dust concentration if requested:
!
item = 257
IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),        &
    conc_dust,                                                &
    row_length,rows,model_levels,0,0,0,0, at_extremity,       &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
    stash_levels,num_stash_levels+1,                          &
    atmos_im,sect,item,                                       &
    icode,cmessage)
  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d (item 257)"//cmessage
  END IF
END IF
!
! ---------------------------------------------------------------------------

IF (icode > 0) THEN
  CALL ereport (routinename, icode, cmessage)
END IF

IF (lhook) CALL dr_hook(routinename,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_aero
END MODULE diagnostics_aero_mod
