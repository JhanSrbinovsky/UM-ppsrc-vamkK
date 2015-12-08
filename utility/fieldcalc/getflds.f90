! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to read requested fields from a UM fieldsfile.

SUBROUTINE GetFlds ( NumFlds,     &  ! in
                     MaxFlds,     &  ! in
                     STCode,      &  ! in
                     MO8Level,    &  ! in
                     FCTime,      &  ! in
                     LBProc,      &  ! in
                     MinsPastHr,  &  ! in
                     MinsDif,     &  ! in
                     Store,       &  ! in
                     PPHdrMod,    &  ! in
                     UMHdr,       &  ! in
                     Fields,      &  ! inout
                     ErrorStatus)    ! inout

! Description:
!   This subroutine reads a given number of fields, matching given
!   criteria from an open UM fieldsfile.  The fields are then store in
!   given locations in the 'Fields' array.
!
! Method:
!   See online documentation.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type,         &
  LSMField
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: STCode  (NumFlds)
INTEGER, INTENT(IN) :: MO8Level(NumFlds)
INTEGER, INTENT(IN) :: FCTime  (NumFlds)
INTEGER, INTENT(IN) :: LBProc  (NumFlds)
INTEGER, INTENT(IN) :: Store   (NumFlds)
INTEGER, INTENT(IN) :: MinsPastHr (NumFlds)
INTEGER, INTENT(IN) :: MinsDif    (NumFlds)
LOGICAL, INTENT(IN) :: PPHdrMod
TYPE(UM_Header_type), INTENT(IN) :: UMHdr

TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "GetFlds"

! Local Variables:
INTEGER :: i                      ! local counter
INTEGER :: ifld
INTEGER :: ihdr
CHARACTER(LEN=80) :: ErrMessage
LOGICAL :: LDecode = .true.
INTEGER :: NumCols
INTEGER :: NumRows
REAL, POINTER :: TempArray(:,:)
INTEGER :: DFldCount
INTEGER :: WFldCount
INTEGER :: PFldCount
LOGICAL :: Field_Not_Found
INTEGER :: NumStore
INTEGER :: vt_days, vt_secs
INTEGER :: dt_days, dt_secs
INTEGER :: LocMinsDif
! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( .NOT. ASSOCIATED(UMHdr % Lookup) ) THEN
  WRITE(6,'(A)') "Cannot get field since file uninitialised - skipping"
  errorstatus = statuswarning
  GOTO 9999
END IF

! Set defaults
Field_Not_Found = .FALSE.
DfldCount = 0
WFldCount = 0
PFldCount = 0
NumStore  = 0
LocMinsDif = 0

!-------------------------------
! Loop through required fields
!-------------------------------
DO i = 1,NumFlds


  ifld = Store(i)     ! Array index where field should be stored


  IF ( ifld > 0 ) THEN
    IF (MO8Level(i) == -1 .OR.  &
        FCTime(i)   == -1) THEN

      WRITE(6,'(A)') "Cannot store multiple fields in one position"
      WRITE(6,'(4(A,I7))') "Store = ", Store(i), " Level = ", MO8Level(i), &
                           " FCTime = ", FCTime(i), " LBProc = ", LBProc(i)
      ifld = 0
    END IF
  END IF

  ! If Store is less than 1 then assume we can just skip
  IF (ifld < 1) THEN
    CYCLE
  END IF
  ! Search to make sure Store is not set
  IF (COUNT(Store(i) == Store(1:i)) > 1) THEN
    DFldCount = DFldCount + 1
    WRITE(6,'(A,I7,2A)') "Field ", ifld, " is already found - ",  &
                         "storing first match"
    CYCLE
  END IF

  NumStore = NumStore + 1

  ! Check to see if existing field has been read into header.
  IF ( ASSOCIATED( Fields(ifld) % RData ) ) THEN
  
    ! See if correct field is already there (to reuse fields - similar check in
    ! copyflds
    IF ( (Fields(ifld) % Hdr % STCode   == STCode(i)  )   .AND. &
         (Fields(ifld) % Hdr % MO8Level == MO8Level(i))   .AND. &
         (Fields(ifld) % Hdr % FCRange  == FCTime(i)  )   .AND. &
         (Fields(ifld) % Hdr % ValidMin == MinsPastHr(i)) .AND. &
         (Fields(ifld) % Hdr % LBProc   == LBProc(i)      .OR.  &
          LBProc(i)                     == -1         ) ) THEN
               ! Calculate data time as reference time.
      IF (MinsDif(i) /= -1) THEN
! DEPENDS ON: time2sec
        CALL time2sec(Fields(ifld) % Hdr % DataYear,  &
                      Fields(ifld) % Hdr % DataMonth, &
                      Fields(ifld) % Hdr % DataDate,  &
                      Fields(ifld) % Hdr % DataHour,  &
                      Fields(ifld) % Hdr % DataMin,   &
                      0,                              &
                      0, 0,                           &
                      dt_days, dt_secs,               &
                      MOD(Fields(ifld) % Hdr % LBTIM,10) == 2)
        ! Calculate validity time using data time as basis.
! DEPENDS ON: time2sec
        CALL time2sec(Fields(ifld) % Hdr % ValidYear,  &
                      Fields(ifld) % Hdr % ValidMonth, &
                      Fields(ifld) % Hdr % ValidDate,  &
                      Fields(ifld) % Hdr % ValidHour,  &
                      Fields(ifld) % Hdr % ValidMin,   &
                      0,                               &
                      0, 0,                            &
                      vt_days, vt_secs,                &
                      MOD(Fields(ifld) % Hdr % LBTIM,10) == 2)
        ! Calculate difference in minutes between data time and validity time.
        ! Requires change in calculation depending on whether its a mean/accum.
        IF (Fields(ifld) % Hdr % LBTim /= 11) THEN
          LocMinsDif = ((dt_days-vt_days)*24*60) + ((dt_secs-vt_secs)/60)
        ELSE
          LocMinsDif = ((vt_days-dt_days)*24*60) + ((vt_secs-dt_secs)/60)
        END IF
      ELSE
        LocMinsDif = -1
      END IF

      IF (MinsDif(i) == LocMinsDif) THEN
        PFldCount = PFldCount + 1
        CYCLE      ! Already exists - move on to next field
      END IF
    END IF
  END IF

  ! Clear space for field
  IF ( ASSOCIATED( Fields(ifld) % RData ) ) THEN
    DEALLOCATE( Fields(ifld) % RData )
    NULLIFY( Fields(ifld) % RData )
  END IF

  !----------------------------------------
  !  Search lookup for the required FIELD
  !----------------------------------------
  ! Search on LBTYP, LBLEV, LBFT - similar check in copyflds
  Fields(ifld) % LookupPos = 0
  DO ihdr = 1,UMHdr % NumFlds
    ! Need to use ABS since a negative is used by getb_flds to signify
    ! field has been interpolated to B grid.
    IF ( (ABS(STCode(i)) == UMHdr % Lookup(ihdr) % STCode)   .AND. &
         (   MO8Level(i) == UMHdr % Lookup(ihdr) % MO8Level) .AND. &
         (     FCTime(i) == UMHdr % Lookup(ihdr) % FCRange)  .AND. &
         ( MinsPastHr(i) == UMHdr % Lookup(ihdr) % ValidMin) .AND. & 
         (     LBProc(i) == UMHdr % Lookup(ihdr) % LBProc    .OR.  &
               LBProc(i) == -1 ) ) THEN
      IF (MinsDif(i) /= -1) THEN
        ! Calculate data time as reference time.
        CALL time2sec(UMHdr % Lookup(ihdr) % DataYear,  &
                      UMHdr % Lookup(ihdr) % DataMonth, &
                      UMHdr % Lookup(ihdr) % DataDate,  &
                      UMHdr % Lookup(ihdr) % DataHour,  &
                      UMHdr % Lookup(ihdr) % DataMin,   &
                      0,                                   &
                      0, 0,                                &
                      dt_days, dt_secs,                    &
                      MOD(UMHdr % Lookup(ihdr) % LBTIM,10) == 2)
        ! Calculate validity time using data time as basis.
        CALL time2sec(UMHdr % Lookup(ihdr) % ValidYear,  &
                      UMHdr % Lookup(ihdr) % ValidMonth, &
                      UMHdr % Lookup(ihdr) % ValidDate,  &
                      UMHdr % Lookup(ihdr) % ValidHour,  &
                      UMHdr % Lookup(ihdr) % ValidMin,   &
                      0,                                    &
                      0, 0,                                 &
                      vt_days, vt_secs,                     &
                      MOD(UMHdr % Lookup(ihdr) % LBTIM,10) == 2)
        ! Calculate difference in minutes between data time and validity time.
        ! Requires change in calculation depending on whether its a mean/accum.
        IF (UMHdr % Lookup(ihdr) % LBTim /= 11) THEN
          LocMinsDif = ((dt_days-vt_days)*24*60) + ((dt_secs-vt_secs)/60)
        ELSE
          LocMinsDif = ((vt_days-dt_days)*24*60) + ((vt_secs-dt_secs)/60)
        END IF
      ELSE
        LocMinsDif = -1
      END IF
      
      IF (MinsDif(i) == LocMinsDif) THEN
        Fields(ifld) % LookupPos = ihdr
        WFldCount = WFldCount + 1
        ! Take first field it finds (LBProc is allowed to be -1 here...)
        EXIT
      END IF
    END IF
  END DO

  IF ( Fields(ifld) % LookupPos == 0 ) THEN
    WRITE( ErrMessage, '(A55,4I5)' )                                &
        "Field not found (STCode,MO8Level,FCRange,MinsPastHr): ",   &
        STCode(i), MO8Level(i), FCTime(i), MinsPastHr(i)
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus, ErrMessage )
! Reset Errorstatus but return warning back to calling routine at end.
    ErrorStatus = StatusOk
    Field_Not_Found=.TRUE.
    CYCLE
  END IF

  Fields(ifld) % Hdr = UMHdr % Lookup( Fields(ifld) % LookupPos )

  ! Allocate space for unpacked field
  ! Need to check if land packed.
  IF ( MOD(Fields(ifld) % Hdr % LBPack / 10, 10) == 2 ) THEN
    IF(ASSOCIATED(LSMField % RData)) THEN
      Fields(ifld) % Hdr % NumCols = LSMField % Hdr % NumCols
      Fields(ifld) % Hdr % NumRows = LSMField % Hdr % NumRows
    ELSE
      ErrorStatus = StatusWarning

      CALL EReport( RoutineName, ErrorStatus, &
                  "Cannot store land/sea packed field without LSM - skipping.")
! Reset Errorstatus but return warning back to calling routine at end.
      ErrorStatus = StatusOk
      Field_Not_Found = .TRUE.
      CYCLE
    END IF
  END IF

  NumCols = Fields(ifld) % Hdr % NumCols
  NumRows = Fields(ifld) % Hdr % NumRows

  ALLOCATE ( Fields(ifld) % RData( NumCols*NumRows,1 ) )

  !----------------------------------------
  ! End of search.  Start retrieving data
  !----------------------------------------

! DEPENDS ON: readfld
  CALL ReadFld( UMHdr, LDecode, PPHdrMod, Fields(ifld),  &
                ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    GO TO 9999
  END IF

  ! Now we need to reshape the read in data to be the correct shape of the
  ! unpacked field.
  TempArray => Fields(ifld) % RData
  NULLIFY(Fields(ifld) % RData)
  ALLOCATE( Fields(ifld) % RData( NumCols,NumRows ) )
  Fields(ifld) % RData = RESHAPE(SOURCE = TempArray,                   &
                                 SHAPE  = (/NumCols,NumRows/))
  DEALLOCATE(TempArray)
  NULLIFY(TempArray)

END DO

WRITE(6,'(2(I7,A))') WFldCount, " stored, out of ", NumStore, " requested."
IF( DFldCount /= 0) THEN
  WRITE(6,'(A,I7,A)') " Warning: ", DFldCount, " duplicate fields found."
END IF
IF( PFldCount /= 0) THEN
  WRITE(6,'(A,I7,A)') " Info: ", PFldCount, " previous fields re-used."
END IF


9999 CONTINUE

IF (Field_Not_Found .AND. ErrorStatus == StatusOK) THEN
  ! Reset errorstatus to pass warning back that a store didnt succeed.
  ErrorStatus = StatusWarning
END IF

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE GetFlds

