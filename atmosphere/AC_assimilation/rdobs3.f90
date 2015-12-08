! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  4 Subroutines in deck : RDOBS, RDOBS2, RDOBS3 and DAYS -----
!LL
!LL  Purpose : Read from ACOBS Files,reformat and place OBS header
!LL            details in COMOBS. The bulk of the required OBS data
!LL            is put into dynamic work array OBS for transmission via
!LL            argument list to GETOBS. OBS is written out to a cache
!LL            file for subsequent reading at later timesteps.
!LL            Thus reread of ACOBS files only required intermittently
!LL            (The routine DAYS does a dd/mm/yy to dayno)
!LL
!LL
!LL S.Bell      <- programmer of some or all of previous code or changes
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical components covered:
!LL
!LL  Project Task : P3
!LL
!LL  External documentation:
!LL
!LLEND------------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: AC Assimilation
MODULE rdobs3_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE RDOBS3 (UNIT_NO,NOBTYPMX,NOBTYP,OBSTYP,NOBS,NDATAV,    &
     &                   NOBLEV,OBLEVTYP,OBS_LEVELS,OBS_DATA,           &
! DUMP_AR2 array dimensions
     &  LEN_FIXHD, LEN_INTHD, LEN_REALHD,                               &
     &  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,         &
     &  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC,LEN2_FLDDEPC,          &
     &  LEN_EXTCNST,  LEN_DUMPHIST,                                     &
     &  LEN_CFI1, LEN_CFI2, LEN_CFI3,                                   &
     &  LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS,                               &
! DUMP_AR2 end
     &                   LEN_DATA,MAX_NDV,TNOBS,MISSD,P_LEVELS,         &
     &                   ICODE,CMESSAGE                                 &
     &                         ,IPT                                     &
     &                             )

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE IO
      IMPLICIT NONE

!L----------------------------------------------------------------------
      INTEGER IPT
      INTEGER                                                           &
     &   UNIT_NO                                                        &
                          ! IN  Unit no of observation file
     &  ,NOBTYPMX                                                       &
                          ! IN  Max no of observation types
     &  ,NOBTYP                                                         &
                          ! OUT No of observation types
     &  ,OBSTYP(NOBTYPMX)                                               &
                          ! OUT Observation types
     &  ,NOBS(NOBTYPMX)                                                 &
                          ! OUT No of Observations
     &  ,NDATAV(NOBTYPMX)                                               &
                          ! OUT No of data values
     &  ,NOBLEV(NOBTYPMX)                                               &
                          ! OUT No of observation levels
     &  ,OBLEVTYP(NOBTYPMX)                                             &
                           !OUT Observation level type
     &  ,LEN_DATA                                                       &
                          ! IN  Dimension of data section
     &  ,TNOBS                                                          &
                          ! OUT Total no of observations
     &  ,P_LEVELS                                                       &
                          ! IN  No of model levels
     &  ,ICODE            ! OUT Return Code

      REAL                                                              &
     &   MISSD                                                          &
                                   ! OUT Real missing data indicator
     &  ,OBS_LEVELS(P_LEVELS+1,*)                                       &
                                   ! OUT Observation levels
     &  ,OBS_DATA(LEN_DATA)        ! OUT Observation data

!DR   REAL DATALEVS(LEN1_LEVDEPC-2,*)     ! OUT Obs levels

      CHARACTER(LEN=*) CMESSAGE ! OUT Error message if ICODE > 0
!L----------------------------------------------------------------------
!L
!*L----------------- COMDECK DUMP_LEN --------------------------------

      INTEGER LEN_FIXHD
      INTEGER LEN_INTHD
      INTEGER LEN_REALHD
      INTEGER LEN1_LEVDEPC, LEN2_LEVDEPC
      INTEGER LEN1_ROWDEPC, LEN2_ROWDEPC
      INTEGER LEN1_COLDEPC, LEN2_COLDEPC
      INTEGER LEN1_FLDDEPC, LEN2_FLDDEPC
      INTEGER LEN_EXTCNST
      INTEGER LEN_DUMPHIST
      INTEGER LEN_CFI1, LEN_CFI2, LEN_CFI3
      INTEGER LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS
!*----------------------------------------------------------------------
!     Dynamic allocate arrays for READDUMP
!*L----------------- COMDECK DUMP_DIM ----------------------------------

      INTEGER FIXHD(LEN_FIXHD)
      INTEGER INTHD(LEN_INTHD)
      REAL    REALHD(LEN_REALHD)
      REAL    LEVDEPC(LEN1_LEVDEPC,LEN2_LEVDEPC)
      REAL    ROWDEPC(LEN1_ROWDEPC,LEN2_ROWDEPC)
      REAL    COLDEPC(LEN1_COLDEPC,LEN2_COLDEPC)
      REAL    FLDDEPC(LEN1_FLDDEPC,LEN2_FLDDEPC)
      REAL    EXTCNST(LEN_EXTCNST)
      REAL    DUMPHIST(LEN_DUMPHIST)
      INTEGER CFI1(LEN_CFI1), CFI2(LEN_CFI2), CFI3(LEN_CFI3)
      INTEGER LOOKUP(LEN1_LOOKUP_OBS,LEN2_LOOKUP_OBS)
!*----------------------------------------------------------------------

!
! EXTERNAL SUBROUTINE CALLS
      EXTERNAL READACOBS
!     LOCAL VARIABLES
      INTEGER JOBT,JLEV
      INTEGER MAX_NDV   !  Maximum no of data values (for obs type)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('RDOBS3',zhook_in,zhook_handle)

!     Go to start of obs file

      CALL SETPOS (UNIT_NO,0,ICODE)

!     Read in the observation file
! DEPENDS ON: readacobs
      CALL READACOBS (UNIT_NO,                                          &
     &  FIXHD,    LEN_FIXHD,                                            &
     &  INTHD,    LEN_INTHD,                                            &
     &  REALHD,   LEN_REALHD,                                           &
     &  LEVDEPC,  LEN1_LEVDEPC,LEN2_LEVDEPC,                            &
     &  ROWDEPC,  LEN1_ROWDEPC,LEN2_ROWDEPC,                            &
     &  COLDEPC,  LEN1_COLDEPC,LEN2_COLDEPC,                            &
     &  FLDDEPC,  LEN1_FLDDEPC,LEN2_FLDDEPC,                            &
     &  EXTCNST,  LEN_EXTCNST,                                          &
     &  DUMPHIST, LEN_DUMPHIST,                                         &
     &  CFI1,     LEN_CFI1,                                             &
     &  CFI2,     LEN_CFI2,                                             &
     &  CFI3,     LEN_CFI3,                                             &
     &  LOOKUP,                                                         &
     &  LEN1_LOOKUP_OBS,LEN2_LOOKUP_OBS,                                &
     &               LEN_DATA,OBS_DATA,                                 &
     &               ICODE,CMESSAGE                                     &
     &                         ,IPT                                     &
     &                             )

      IF (ICODE >  0) GOTO 9999

      TNOBS  = INTHD(28)   !  Total no of observations
      NOBTYP = INTHD(32)   !  No of observation types
      MISSD  = REALHD(29)  !  Real MDI used in observations

      IF (NOBTYP >  0) THEN

        DO JOBT = 1,NOBTYP
          OBSTYP(JOBT) = LOOKUP(65,JOBT)   ! Observation Types
          NOBS  (JOBT) = LOOKUP(66,JOBT)   ! No of obs
          NDATAV(JOBT) = LOOKUP(67,JOBT)   ! No of data values
          NOBLEV(JOBT) = LEVDEPC(2,JOBT)   ! No of obs levels
          OBLEVTYP(JOBT) = LEVDEPC(1,JOBT) ! Level type
          DO JLEV =1,LEN1_LEVDEPC-2
            OBS_LEVELS(JLEV,JOBT) = LEVDEPC(2+JLEV,JOBT) ! Obs levels
          ENDDO
        ENDDO

!       MAXNLEV1 = LEN1_LEVDEPC-2  ! Max no of levels + 1

      ENDIF
9999  CONTINUE
      IF (lhook) CALL dr_hook('RDOBS3',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE RDOBS3
END MODULE rdobs3_mod
