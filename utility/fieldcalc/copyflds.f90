! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines for copying fields from one UM fieldsfile to another

SUBROUTINE CopyFlds( NumFlds_In,  &  ! in
                     MaxFlds,     &  ! in
                     STCode,      &  ! in
                     MO8Level,    &  ! in
                     MO8LevelLwr, &  ! in
                     MO8LevelUpr, &  ! in
                     FCTime,      &  ! in
                     FCTimeStart, &  ! in
                     FCTimeFreq,  &  ! in
                     FCTimeEnd,   &  ! in
                     LBProc,      &  ! in
                     MinsPastHr,  &  ! in
                     MinsDif,     &  ! in
                     Factor,      &  ! in
                     Store,       &  ! in
                     PPHdrMod,    &  ! in
                     UMHdr_in,    &  ! in
                     UMHdr_out,   &  ! inout
                     Fields,      &  ! inout
                     ErrorStatus )   ! inout

! Description:
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Field_type,          &
  PP_Header_type,         &
  UM_Header_type,         &
  LSMField
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds_In
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: STCode     (MaxFlds)
INTEGER, INTENT(INOUT) :: MO8Level    (MaxFlds)
INTEGER, INTENT(IN)    :: MO8LevelLwr (MaxFlds)
INTEGER, INTENT(IN)    :: MO8LevelUpr (MaxFlds)
INTEGER, INTENT(INOUT) :: FCTime      (MaxFlds)
INTEGER, INTENT(IN)    :: FCTimeStart (MaxFlds)
INTEGER, INTENT(IN)    :: FCTimeFreq  (MaxFlds)
INTEGER, INTENT(IN)    :: FCTimeEnd   (MaxFlds)
INTEGER, INTENT(IN) :: LBProc     (MaxFlds)
INTEGER, INTENT(IN) :: Store      (MaxFlds)
INTEGER, INTENT(IN) :: MinsPastHr (MaxFlds)
INTEGER, INTENT(IN) :: MinsDif    (MaxFlds)
REAL,    INTENT(IN) :: Factor     (MaxFlds)
LOGICAL, INTENT(IN) :: PPHdrMod
TYPE(UM_Header_type), INTENT(IN) :: UMHdr_in

TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr_out
TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "CopyFlds"
REAL, PARAMETER :: VSmall = 1.e-9
CHARACTER(LEN=*), PARAMETER :: CPackType = "WGDOS"

! Local Variables:
INTEGER :: NumLookups
INTEGER :: MaxFldSize
INTEGER :: i
INTEGER :: ihdr
INTEGER :: NumScl
INTEGER :: CFldCount
INTEGER :: WFldCount
INTEGER :: SFldCount
INTEGER :: NumFlds
INTEGER :: dt_days, dt_secs
INTEGER :: vt_days, vt_secs
INTEGER :: LocMinsDif
REAL    :: PackAcc
CHARACTER(LEN=80) :: ErrMessage
LOGICAL :: LDecode
LOGICAL :: LCopyall ! True if we are really copying all fields.
TYPE(PP_Field_type) :: TempField
TYPE(PP_Field_type) :: XpndField
TYPE(PP_Field_type) :: CompField
LOGICAL, ALLOCATABLE :: FoundField(:)
LOGICAL              :: l_return_warn

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Default return a warning to false
l_return_warn = .FALSE.

! Default specifying Numflds to be the input.
LCopyAll = .FALSE.
NumFlds = NumFlds_In
NumLookups = UMHdr_in % NumFlds

! Deal with whether we are being called by Copyall.  Numflds could be 0 so we
! have added passed down numflds + maxflds + 1 to signify copyall.
IF ( NumFlds_In > MaxFlds ) THEN
  ! This signals to copy all fields
  LCopyall = .TRUE.
  NumFlds = NumFlds_In - Maxflds - 1
END IF

! Allocate array to store whether field is found.
ALLOCATE(FoundField(Numflds))
FoundField(:) = .FALSE.

! Allocate temporary field
MaxFldSize = MAXVAL( UMHdr_in % Lookup(1:NumLookups) % NumCols ) * &
             MAXVAL( UMHdr_in % Lookup(1:NumLookups) % NumRows )
ALLOCATE( TempField % RData( MaxFldSize, 1) )
NumScl    = 0
CFldCount = 0
SFldCount = 0

! Lets store all the fields requested in store array.

CALL GetFlds(NumFlds,MaxFlds,StCode,MO8Level,FCTime,LBProc,MinsPastHr,  &
             MinsDif, Store,PPHdrMod,UMHdr_in,Fields,ErrorStatus)
! Add check here... if problem
IF (ErrorStatus /= StatusOk) THEN
  errmessage = 'Problem storing field.  Trying to continue copying.'
  CALL ereport(routinename, errorstatus, errmessage)
! Temporarily remove error but set it again at the end.
  errorstatus = statusok
  l_return_warn = .TRUE.
END IF

!--------------------
! Loop through list
!--------------------
DO ihdr = 1,NumLookups
! If we are copying all we can read the field here.
  IF (LCopyall) THEN
    CFldCount = CFldCount + 1
    !---------------------------------------
    ! Found matching field.  Retrieve data
    !---------------------------------------
    !CFldCount = CFldCount + 1
    ! We may have already retrieved data in getflds above.
    TempField % LookupPos = ihdr
    TempField % Hdr = UMHdr_in % Lookup( TempField % LookupPos )
    ! Read field into temporary space.
    LDecode = .FALSE.
! DEPENDS ON: readfld
    CALL ReadFld( UMHdr_in, LDecode, PPHdrMod, TempField,   &
                    ErrorStatus )
    IF ( ErrorStatus /= StatusOK ) THEN
      WRITE( ErrMessage, '(A36,I4)' ) &
             "Could not read input field position ", ihdr
      Errorstatus = StatusWarning

      CALL EReport( RoutineName, ErrorStatus, ErrMessage )
      l_return_warn = .TRUE.
      CYCLE
    END IF
  END IF

  !----------------------------------------
  !  Search lookup for the required FIELD
  !----------------------------------------
  ! -1 is a wildcard for MO8Level and FCRange. STCode must be fixed.
  ! There is a similar check in getflds.
  DO i = 1, NumFlds

! See if we are using range and set MO8Level to correct value
    IF ( UMHdr_in % Lookup(ihdr) % MO8Level >= MO8LevelLwr(i) .AND. &
         UMHdr_in % Lookup(ihdr) % MO8Level <= MO8LevelUpr(i)) THEN
      MO8Level(i) = UMHdr_in % Lookup(ihdr) % MO8Level
    END IF

! See if we are using range to specify time, set FCTime to correct value.
! If MinsPastHr is required specify each MinsPastHr separately in namelist.
    IF ( UMHdr_in % Lookup(ihdr) % FCRange >= FCTimeStart(i) .AND. &
         UMHdr_in % Lookup(ihdr) % FCRange <= FCTimeEnd(i)   .AND. &
         MODULO(UMHdr_in % Lookup(ihdr) % FCRange - FCTimeStart(i),&
                FCTimeFreq(i)) == 0) THEN
      FCTime(i) = UMHdr_in % Lookup(ihdr) % FCRange
    END IF

    IF ( STCode(i)      == UMHdr_in % Lookup(ihdr) % STCode       .AND. &
         (MO8Level(i)   == UMHdr_in % Lookup(ihdr) % MO8Level     .OR.  &
          MO8Level(i)   == -1                                 )   .AND. &
         (FCTime(i)     == UMHdr_in % Lookup(ihdr) % FCRange      .OR.  &
          FCTime(i)     == -1                                 )   .AND. &
         (MinsPastHr(i) == UMHdr_in % Lookup(ihdr) % ValidMin     .OR.  &
          MinsPastHr(i) == -1                                 )   .AND. &
         (LBProc(i)     == UMHdr_in % Lookup(ihdr) % LBProc       .OR.  &
          LBProc(i)     == -1                                 ) ) THEN

      ! If we have a valid minsdif value then check it.
      IF ( MinsDif(i)    /= -1 ) THEN
        ! Calculate data time as reference time.
! DEPENDS ON: time2sec
        CALL time2sec(UMHdr_in % Lookup(ihdr) % DataYear,  &
                      UMHdr_in % Lookup(ihdr) % DataMonth, &
                      UMHdr_in % Lookup(ihdr) % DataDate,  &
                      UMHdr_in % Lookup(ihdr) % DataHour,  &
                      UMHdr_in % Lookup(ihdr) % DataMin,   &
                      0,                                   &
                      0, 0,                                &
                      dt_days, dt_secs,                    &
                      MOD(UMHdr_in % Lookup(ihdr) % LBTIM,10) == 2)
        ! Calculate validity time using data time as basis.
! DEPENDS ON: time2sec
        CALL time2sec(UMHdr_in % Lookup(ihdr) % ValidYear,  &
                      UMHdr_in % Lookup(ihdr) % ValidMonth, &
                      UMHdr_in % Lookup(ihdr) % ValidDate,  &
                      UMHdr_in % Lookup(ihdr) % ValidHour,  &
                      UMHdr_in % Lookup(ihdr) % ValidMin,   &
                      0,                                    &
                      0, 0,                                 &
                      vt_days, vt_secs,                     &
                      MOD(UMHdr_in % Lookup(ihdr) % LBTIM,10) == 2)
        ! Calculate difference in minutes between data time and validity time.
        IF (UMHdr_in % Lookup(ihdr) % LBTim /= 11) THEN
          LocMinsDif = ((dt_days-vt_days)*24*60) + ((dt_secs-vt_secs)/60)
        ELSE
          LocMinsDif = ((vt_days-dt_days)*24*60) + ((vt_secs-dt_secs)/60)
        END IF
      ELSE
        LocMinsDif = -1
      END IF

        ! Check to see if time difference is the correct one
      IF ( MinsDif(i) /= LocMinsDif ) THEN
        ! Move to next field since this isnt the correct field
          CYCLE
      END IF
           
      IF ( ABS(Factor(i)) > VSmall ) THEN
        NumScl = NumScl + 1
      END IF
      !---------------------------------------
      ! Found matching field.  Retrieve data
      !---------------------------------------
      FoundField(i) = .TRUE.

      IF (.NOT. LCopyall) THEN
        CFldCount = CFldCount + 1
        TempField % LookupPos = ihdr
        TempField % Hdr = UMHdr_in % Lookup( TempField % LookupPos )
        ! Read field into temporary space.
        LDecode = .FALSE.
! DEPENDS ON: readfld
        CALL ReadFld( UMHdr_in, LDecode, PPHdrMod, TempField,   &
                      ErrorStatus )
        IF ( ErrorStatus /= StatusOK ) THEN
          WRITE( ErrMessage, '(A36,I4)' ) &
                 "Could not read input field position ", ihdr
          ErrorStatus = StatusWarning

          CALL EReport( RoutineName, ErrorStatus, ErrMessage )
          l_return_warn = .TRUE.
          CYCLE
        END IF
      END IF

      IF ( ABS(Factor(i)) > VSmall ) THEN ! Unpack
        ! Could be improved since we might already have this stored.
        PackAcc = TempField % Hdr % BAcc
        CALL GetFlds(1,1, &
                     TempField % Hdr % STCode,       &
                     TempField % Hdr % MO8Level,     &
                     TempField % Hdr % FCRange,      &
                     TempField % Hdr % LBProc,       &
                     TempField % Hdr % ValidMin,     &
                     LocMinsDif,                     &
                     1,PPHdrMod,UMHdr_in,XpndField,ErrorStatus)

        ! Scaling of output required
        SFldCount = SFldCount + 1
        XpndField % RData = Factor(i) * XpndField % RData
! DEPENDS ON: pack_single
        CALL Pack_Single( PackAcc, CPackType, XpndField, &
                          CompField, ErrorStatus )
        TempField % Hdr = CompField % Hdr
        TempField % RData( 1:SIZE(CompField % RData), :) = &
                                                CompField % RData
        DEALLOCATE( CompField % RData )
        NULLIFY( CompField % RData )
        DEALLOCATE( XpndField % RData )
        NULLIFY( XpndField % RData )
      END IF
! Write it out if we are just copying what is specified in NumFlds.
      IF (.NOT. LCopyall) THEN
! DEPENDS ON: fldout
        CALL FldOut( TempField, UMHdr_out, ErrorStatus )
      END IF
      EXIT
    END IF

  END DO
  ! Write every field out.
  IF (LCopyall) THEN
    CALL FldOut( TempField, UMHdr_out, ErrorStatus )
  END IF

END DO

! Now print out diagnostic info if we didnt find it.  For entries with
! Store(i) > 0 we should have reported error in getflds instead.
DO i = 1, NumFlds
  IF ( (.NOT. FoundField(i)) .AND. Store(i) < 1) THEN
    WRITE( ErrMessage, '(A46,4I5)' )                         &
        "Field not found for STCode, MO8Level, FCRange, MinsPastHr: ",   &
        STCode(i), MO8Level(i), FCTime(i), MinsPastHr(i)
    WRITE(6,*) ErrMessage
  END IF
END DO


WRITE(6,*) CFldCount, " fields copied."
WRITE(6,*) SFldCount, " scaled, out of ", NumScl, " requested."

9999 CONTINUE

IF(ALLOCATED(FoundField)) DEALLOCATE( FoundField )
IF(ASSOCIATED(TempField % Rdata)) THEN
  DEALLOCATE( TempField % RData )
  NULLIFY( TempField % Rdata )
END IF

! If error with storing fields then we at least try and complete action but
! return an error so further actions require intervention of a error clearing
! command.
IF (l_return_warn) THEN
  errorstatus = statuswarning
END IF

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE CopyFlds

