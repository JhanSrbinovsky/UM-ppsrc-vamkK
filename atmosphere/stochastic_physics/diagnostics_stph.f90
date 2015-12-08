! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: stochastic physics
MODULE diagnostics_stph_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE diagnostics_stph( row_length, rows, model_levels,            &
                             n_rows, at_extremity, stph_diag,           &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                             stashwork35)

! Purpose:
!  Calculates diagnostics generated from stochastic physics routines
!  (UM section 35).

! Method:
! Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the stochastic
! physics routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing.

!  Diagnostics currently available:

USE stph_diag_mod,  ONLY:                                               &
    strstphdiag

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE Field_Types
USE Submodel_Mod

IMPLICIT NONE


! Arguments with Intent IN. ie: Input variables.

LOGICAL                                                                 &
  at_extremity(4)
                   ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

INTEGER                                                                 &
  row_length                                                            &
                   ! number of points on a row
, rows                                                                  &
                   ! number of rows in a theta field
, n_rows                                                                &
                   ! number of rows in a v field
, model_levels
                   ! number of model levels

!     Declaration of Stochastic Physics diagnostics.
TYPE (strstphdiag) :: stph_diag

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
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!  Global Variables:----------------------------------------------------

! Diagnostics info
REAL                                                                    &
 stashwork35(*)
                  ! STASH workspace
INTEGER                                                                 &
  im_index        ! internal model index

! Local variables

INTEGER                                                                 &
  icode                                                                 &
                  ! Return code  =0 Normal exit  >1 Error
 ,item                                                                  &
                  ! STASH item
 ,sect            ! STASH section
PARAMETER( sect = 35 ) ! for stochastic physics

CHARACTER(LEN=80)  cmessage

CHARACTER(LEN=*) routinename
PARAMETER ( routinename='diagnostics_stph')

INTEGER, PARAMETER :: pnorth=1
INTEGER, PARAMETER :: peast =2
INTEGER, PARAMETER :: psouth=3
INTEGER, PARAMETER :: pwest =4
INTEGER, PARAMETER :: nodomain = -1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------

! Initialise error status
IF (lhook) CALL dr_hook('DIAGNOSTICS_STPH',zhook_in,zhook_handle)
icode = 0
im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! DIAG.35001 Copy U wind after SKEB2 to stashwork
! ----------------------------------------------------------------------
item = 1  ! U wind after SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_u,                                              &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 1)"//cmessage
   END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35002 Copy V wind after SKEB2 to stashwork
! ----------------------------------------------------------------------

item = 2  ! V wind after SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_v,                                              &
        row_length,n_rows,model_levels,0,0,0,0, at_extremity,           &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 2)"//cmessage
   END IF

END IF  !  sf(item,sect)
! ----------------------------------------------------------------------
! DIAG.35003 Copy SKEB2: Full u increment to stashwork
! ----------------------------------------------------------------------
item = 3  ! skeb2 full u wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_u_incr,                                         &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 3)"//cmessage
   END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35004 Copy SKEB2: full V INCR to stashwork
! ----------------------------------------------------------------------

item = 4  ! skeb2 full v wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_v_incr,                                         &
        row_length,n_rows,model_levels,0,0,0,0, at_extremity,           &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 4)"//cmessage
   END IF

END IF  !  sf(item,sect)
! ----------------------------------------------------------------------
! DIAG.35005 Copy SKEB2: Rotational U INCR to stashwork
! ----------------------------------------------------------------------
item = 5  ! skeb2 rotational u wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_u_rot,                                          &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 1)"//cmessage
   END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35006 Copy SKEB2: Rotational V INCR to stashwork
! ----------------------------------------------------------------------

item = 6  ! skeb2 rotational  v wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_v_rot,                                          &
        row_length,n_rows,model_levels,0,0,0,0, at_extremity,           &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 2)"//cmessage
   END IF

END IF  !  sf(item,sect)
! ----------------------------------------------------------------------
! DIAG.35007 Copy SKEB2: Divergent U INCR to stashwork
! ----------------------------------------------------------------------
item = 7  ! skeb2 divergent u wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_u_div,                                          &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 3)"//cmessage
   END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35008 Copy SKEB2: Divergent V INCR to stashwork
! ----------------------------------------------------------------------

item = 8  ! skeb2 divergent v wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_v_div,                                          &
        row_length,n_rows,model_levels,0,0,0,0, at_extremity,           &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 4)"//cmessage
   END IF
     END IF

 ! ----------------------------------------------------------------------
! DIAG.35009 Copy SKEB2: dissipation field from smagorinsky code to stashwork
! ----------------------------------------------------------------------
item = 9  ! skeb2 dissipation field from smagorinsky
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_disp_smag,                                      &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 1)"//cmessage
   END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35010 Copy SKEB2: dissipation field from convection to stashwork
! ----------------------------------------------------------------------

item = 10  ! skeb2 dissipation field from convection
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_disp_conv,                                      &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 2)"//cmessage
   END IF

END IF  !  sf(item,sect)
! ----------------------------------------------------------------------
! DIAG.35011 Copy SKEB2: dissipation field from SKEB1 to stashwork
! ----------------------------------------------------------------------
item = 11  ! skeb2 dissipation field from SKEB1-type
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_disp_skeb1,                                     &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 3)"//cmessage
   END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35012 Copy SKEB2: smoothed modulating field to stashwork
! ----------------------------------------------------------------------

item = 12  ! skeb2 smoothed modulating field
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_smodfield,                                      &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35013 Copy SKEB2: final stream function field
! ----------------------------------------------------------------------
item = 13  ! skeb2 raw modulating field
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_streamfunction,                                 &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 3)"//cmessage
   END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35014 Copy SKEB2: intial random pattern
! ----------------------------------------------------------------------

item = 14  ! skeb2 smoothed modulating field
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork35(si(item,sect,im_index)),               &
        stph_diag%skeb2_random_pattern,                                 &
        row_length,rows,model_levels,0,0,0,0, at_extremity,             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35015 Copy SKEB2: Vert Integ. KE of initial SF forcing
! ----------------------------------------------------------------------

item = 15  ! Vert Integ. KE of initial SF forcing
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag (stashwork35(si(item,sect,im_index)),                  &
        stph_diag%skeb2_ke_psif,                                        &
        row_length,rows,0,0,0,0, at_extremity,                          &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35016 Copy SKEB2: Vert Integ. KE of numerical diss
! ----------------------------------------------------------------------

item = 16  ! Vert Integ. KE of numerical diss
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag (stashwork35(si(item,sect,im_index)),                  &
        stph_diag%skeb2_ke_sdisp,                                       &
        row_length,rows,0,0,0,0, at_extremity,                          &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35017 Copy SKEB2: Vert Integ. KE of convection diss
! ----------------------------------------------------------------------

item = 17  ! Vert Integ. KE of convection diss
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag (stashwork35(si(item,sect,im_index)),                  &
        stph_diag%skeb2_ke_cdisp,                                       &
        row_length,rows,0,0,0,0, at_extremity,                          &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35018 Copy SKEB2: Vert Integ. KE of 2nd convection diss
! ----------------------------------------------------------------------

item = 18  ! Vert Integ. KE of mflx-based "w" convection diss
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag (stashwork35(si(item,sect,im_index)),                  &
        stph_diag%skeb2_ke_kdisp,                                       &
        row_length,rows,0,0,0,0, at_extremity,                          &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35019 Copy SKEB2: Vert Integ. KE of modulated SF forcing
! ----------------------------------------------------------------------

item = 19  ! Vert Integ. KE of modulated SF forcing
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag (stashwork35(si(item,sect,im_index)),                  &
        stph_diag%skeb2_ke_m_psif,                                      &
        row_length,rows,0,0,0,0, at_extremity,                          &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35020 Copy SKEB2: Vert Integ. KE of total wind incr before SKEB2
! ----------------------------------------------------------------------

item = 20  ! Vert Integ. KE of total wind incr before SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag (stashwork35(si(item,sect,im_index)),                  &
        stph_diag%skeb2_ke_prewindincr,                                 &
        row_length,rows,0,0,0,0, at_extremity,                          &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35021 Copy SKEB2: Vert Integ. KE of wind incr from SKEB2
! ----------------------------------------------------------------------

item = 21  ! Vert Integ. KE of wind incr from SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag (stashwork35(si(item,sect,im_index)),                  &
        stph_diag%skeb2_ke_windincr,                                    &
        row_length,rows,0,0,0,0, at_extremity,                          &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35022 Copy SKEB2: Vert Integ. KE of total wind incr after SKEB2
! ----------------------------------------------------------------------

item = 22  ! Vert Integ. KE of total wind incr after SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag (stashwork35(si(item,sect,im_index)),                  &
        stph_diag%skeb2_ke_postwindincr,                                &
        row_length,rows,0,0,0,0, at_extremity,                          &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
      cmessage=": error in copydiag(item 4)"//cmessage
   END IF
 END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
!  single point error handling.
! ----------------------------------------------------------------------

 IF (icode /= 0) THEN

   CALL ereport(routinename,icode,cmessage)
 END IF

 IF (lhook) CALL dr_hook('DIAGNOSTICS_STPH',zhook_out,zhook_handle)
 RETURN
 END SUBROUTINE diagnostics_stph
END MODULE diagnostics_stph_mod
