! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to translate stash codes to GRIB codes.

!=======================================================================

! Description:
!   Functions to query the translation table.  This was an attempt to restrict
!   just reading and modifying values directly.  We have tried to reuse
!   existing translation table but have slightly modified it.
!   The translation table is structured like:
! stash, fc, disc, cat, item, level type, level, PDT#+octs, paramid, manip
! 
! Where:
!   stash = superstash code ([model][section][item])
!   Example: U COMPT OF WIND      01       00   002 = 0100002
!   fc    = fieldcode (old pp field code, structured like [lvc][ppcode])
!   Example: 10 METRE WIND U-COMP B GRID                      1   00057
!   disc  = GRIB2 discipline (unused)
!   cat   = GRIB2 category (unused)
!   item  = GRIB2 item (unused)
!   level type = GRIB2 level type (unused)
!   level = GRIB2 level (unused)
!   PDT   = GRIB2 product template (unused)
!   paramid = GRIB-API parameter, this results in not needing all unused above.
!   manip   = modification to data (e.g. multiply by factor) (optional)
!
!   The list of functions are:
!     find_stcode_to_paramid
!     find_ppcode_to_paramid
!     find_paramid_to_stcode
!     find_stcode_to_level
!     find_stcode_to_ppcode
!     find_stcode_to_lbvc
!     find_paramid_to_manip_val
!     find_paramid_to_manip_op
! 
! Method:
!   We first read in the translation table.  We then use functions to query it.
!  
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6


MODULE ff_to_grib2_mod

USE grib_api

IMPLICIT NONE

PRIVATE

! Number of bits of info
INTEGER, PARAMETER :: info_size = 11
! STASH code
INTEGER, PARAMETER :: stcode    = 1
! PP fieldcode
INTEGER, PARAMETER :: fieldcode = 2
! GRIB2 discipline number
INTEGER, PARAMETER :: disc      = 3
! GRIB2 category number 
INTEGER, PARAMETER :: cat       = 4
! GRIB2 item number
INTEGER, PARAMETER :: item      = 5
! required for pp files
INTEGER, PARAMETER :: lbtyp     = 6
! for fields where level is implied (e.g. 1.5m TEMPERATURE)
INTEGER, PARAMETER :: level     = 7
! template number
INTEGER, PARAMETER :: pdt       = 8
! ECMWF parameter number
INTEGER, PARAMETER :: paramid   = 9
! Scaling from STASH to GRIB
INTEGER, PARAMETER :: manip_val = 10
! Processing number
INTEGER, PARAMETER :: lbproc    = 11

TYPE stash_to_grib2
! Array holding info
  INTEGER   :: info(info_size)
  CHARACTER :: manip_op
! Lets setup a linked list - performance shouldnt be a problem.
  TYPE(stash_to_grib2), POINTER :: next_stash => NULL()
  TYPE(stash_to_grib2), POINTER :: next_grib2 => NULL()
END TYPE stash_to_grib2

TYPE lbvc_to_grib2
  INTEGER :: lbvc
  INTEGER :: grib2_lt
  TYPE(lbvc_to_grib2), POINTER :: next_lbvc => NULL()
END TYPE

TYPE(stash_to_grib2), POINTER, PUBLIC :: first_st_to_grib2   => NULL()
TYPE(lbvc_to_grib2),  POINTER, PUBLIC :: first_lbvc_to_grib2 => NULL()
TYPE(stash_to_grib2), POINTER         :: curr_st_to_grib2    => NULL()
TYPE(lbvc_to_grib2),  POINTER         :: curr_lbvc_to_grib2  => NULL()

! Interface

Interface ff_to_grib_data_manip
  Module Procedure ff_to_grib_data_manip_1D
End Interface

Interface grib_to_ff_data_manip
  Module Procedure grib_to_ff_data_manip_2D
End Interface


! Make types publicly available
PUBLIC stash_to_grib2
PUBLIC lbvc_to_grib2

! Make routines publicly available
PUBLIC read_st_to_grib2

! Make functions publicly available
PUBLIC find_stcode_to_paramid
PUBLIC find_superppcode_to_paramid
PUBLIC find_stcode_to_level
PUBLIC find_stcode_to_ppcode
PUBLIC find_stcode_to_lbvc
PUBLIC find_paramid_to_stcode
PUBLIC ff_to_grib_data_manip
PUBLIC grib_to_ff_data_manip

CONTAINS
SUBROUTINE read_st_to_grib2(Errorstatus)

USE Err_mod, ONLY: &
  statusOK,        &
  statusFatal,     &
  endoffile
  
USE ereport_mod, ONLY: ereport  

USE filenamelength_mod, ONLY: filenamelength

INTEGER, INTENT(INOUT)        :: errorstatus
CHARACTER(LEN=filenamelength) :: c_st_to_grib2
CHARACTER(LEN=80)             :: errmessage
INTEGER, PARAMETER            :: unitno_start = 100
INTEGER                       :: unitno
INTEGER                       :: readstatus
CHARACTER(LEN=512)            :: aline
CHARACTER(LEN=100)            :: string_val
INTEGER                       :: cmtpos
INTEGER                       :: delimpos
INTEGER                       :: delimpos_old
INTEGER                       :: i
INTEGER                       :: ind_loc
CHARACTER(LEN=*), PARAMETER   :: routinename = 'READ_ST_TO_GRIB2'
REAL                          :: temp_real
LOGICAL                       :: unitno_open
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

IF ( Errorstatus /= StatusOk ) THEN
! Previous error - do not proceed
  GOTO 9999
END IF

IF ( ASSOCIATED(first_st_to_grib2) ) THEN
! Already ran this - lets skip.
  GOTO 9999
ELSE
  ALLOCATE(first_st_to_grib2)
END IF

curr_st_to_grib2 => first_st_to_grib2

Call FORT_GET_ENV('ST_TO_GRIB2',12, c_st_to_grib2,filenamelength, ErrorStatus)
IF ( Errorstatus /= StatusOk ) THEN
  ErrMessage = 'Cannot find translation table from STASH to GRIB2'
  errorstatus = statusfatal
  ! Ereport will reset this to statusOK (i.e. 0)
  CALL ereport( routinename, errorstatus, errmessage)
END IF

! Find available unitno - fixes issues on systems where stdout/err/in are not
! in the locations assumed to be used.
unitno = unitno_start
DO
  INQUIRE(unitno, OPENED=unitno_open)
  IF (.NOT. unitno_open) EXIT
  unitno = unitno + 1
END DO


OPEN( UNIT   = unitno,        &
      ACCESS = "SEQUENTIAL",  &
      ACTION = "READ",        &
      FILE   = c_st_to_grib2, &
      FORM   = "FORMATTED",   &
      IOSTAT = ErrorStatus,   &
      STATUS = "OLD" )

IF ( Errorstatus /= StatusOk ) THEN
  ErrMessage = 'Cannot open translation table from STASH to GRIB2'
  errorstatus = statusfatal
  ! Ereport will reset this to statusOK (i.e. 0)
  CALL ereport( routinename, errorstatus, errmessage)
END IF


DO
  READ(unitno, '(a)', IOSTAT=readstatus) aline
  IF (readstatus == EndofFile) THEN
    EXIT
  END IF
! Find if a comment exists
  cmtpos   = SCAN(aline, '#')
! Delete any comments from input line.
  IF (cmtpos > 0) THEN
    aline(cmtpos:) = ""
  END IF

  IF (LEN_TRIM(aline) == 0) THEN
    ! This was an comment or empty line move to next line.
    CYCLE
  END IF

! Allocate space for record.
  IF (.NOT. ASSOCIATED(curr_st_to_grib2)) THEN
    ALLOCATE(curr_st_to_grib2)
  END IF
 
  delimpos_old = 0
! Set defaults for information.
  curr_st_to_grib2 % info = (/ (-1,i=1,info_size) /)

! Set defaults for real numbers.
  curr_st_to_grib2 % info(manip_val) = &
    TRANSFER(1.0,curr_st_to_grib2 % info(manip_val))
  curr_st_to_grib2 % manip_op = "*"
  curr_st_to_grib2 % info(level) = &
    TRANSFER(-1.0,curr_st_to_grib2 % info(level))

! Now process line - assume it will have info_size fields.
  DO i = 1, info_size
   
! Check that we previous delimiter is within string.  Okay if we are the last
! column since it is optional (manip operation).
    IF (delimpos_old >= LEN_TRIM(aline)) THEN
      IF (i < manip_val) THEN
        WRITE(0,*) "INFO: Ran out of columns - using default for ", &
                    curr_st_to_grib2 % info(1)
      END IF
      EXIT
    END IF

    ! Find new delimpos.
    delimpos = delimpos_old + SCAN(aline(delimpos_old+1:), ',')
    
    ! Only last field should have no comma at the end of the line.
    IF (delimpos == delimpos_old) THEN
    ! Assume this is the last field so set to one past end of line.
      delimpos = len_trim(aline) + 1
    END IF

    ! Read into segment of line into another string.
    READ(aline(delimpos_old+1:delimpos-1), '(a)') string_val
    IF (i == manip_val) THEN
      ! Default is to multiply by 1.0
      IF (SCAN(string_val,'+-*/') > 0) THEN
        ! We have a mathematical operation and a number to use.
        ind_loc = SCAN(string_val,'+-*/')
        curr_st_to_grib2 % manip_op = string_val(ind_loc:ind_loc)
        READ(string_val(ind_loc+1:),'(f10.3)') temp_real 
        curr_st_to_grib2 % info(i) = TRANSFER(temp_real, &
                                     curr_st_to_grib2 % info(i))
      ELSE IF(SCAN(string_val,'EL') > 0) THEN
        ind_loc = SCAN(string_val,'EL')
        ! We have a description of an operation (E for Exp, L for Log)
        curr_st_to_grib2 % info(i) = 0.0
        curr_st_to_grib2 % manip_op = string_val(ind_loc:ind_loc)
      ELSE
        WRITE(6,*) "WARNING: Unrecognised manipulation option for ", &
                   curr_st_to_grib2 % info(stcode)
      END IF
    ! If item has no characters  x;. then treat as integer.
    ELSE IF (SCAN(string_val,'x;.') == 0) THEN
      READ(string_val,'(i10)') curr_st_to_grib2 % info(i)
    ! If item has a period treat as real.
    ELSE IF (SCAN(string_val,'.') /= 0) THEN
      READ(string_val,'(f10.3)') temp_real
      curr_st_to_grib2 % info(i) = TRANSFER(temp_real, &
                                     curr_st_to_grib2 % info(i))
    
    END IF
    delimpos_old = delimpos
  END DO
  ALLOCATE(curr_st_to_grib2 % next_stash)
  curr_st_to_grib2 => curr_st_to_grib2 % next_stash
END DO
! Tidy up after ourselves (we have to allocate next before pointing to it
! at the end of the loop).
IF (ASSOCIATED(curr_st_to_grib2 % next_stash)) THEN
  DEALLOCATE(curr_st_to_grib2 % next_stash)
END IF

NULLIFY(curr_st_to_grib2 % next_stash)

CLOSE(unitno)

9999 CONTINUE

END SUBROUTINE read_st_to_grib2

SUBROUTINE read_lbvc_to_grib2

END SUBROUTINE read_lbvc_to_grib2

FUNCTION find_stcode_to_paramid (stcode_in,lbproc_in_arg)

INTEGER :: find_stcode_to_paramid
INTEGER :: stcode_in
INTEGER, OPTIONAL :: lbproc_in_arg
INTEGER(KIND=KindOfInt) :: paramid_out

INTEGER :: lbproc_in

paramid_out = -1
! If missing in argument call the translation table will have -1 too,
lbproc_in = -1

IF (PRESENT(lbproc_in_arg)) THEN
  lbproc_in = lbproc_in_arg
END IF

curr_st_to_grib2 => first_st_to_grib2

DO WHILE ( ASSOCIATED(curr_st_to_grib2 ) )
  IF (curr_st_to_grib2 % info(stcode) == stcode_in) THEN
    IF (curr_st_to_grib2 % info(lbproc) == -1 .OR. &
        curr_st_to_grib2 % info(lbproc) == lbproc_in) THEN
      paramid_out = curr_st_to_grib2 % info(paramid)
      EXIT
    END IF
  END IF
  curr_st_to_grib2 => curr_st_to_grib2 % next_stash
END DO

find_stcode_to_paramid = paramid_out

END FUNCTION find_stcode_to_paramid

FUNCTION find_superppcode_to_paramid (superppcode_in)

INTEGER :: find_superppcode_to_paramid
INTEGER :: superppcode_in
INTEGER(KIND=KindOfInt) :: paramid_out

paramid_out = -1

curr_st_to_grib2 => first_st_to_grib2

DO WHILE ( ASSOCIATED(curr_st_to_grib2 ) )
  IF (curr_st_to_grib2 % info(fieldcode) == superppcode_in) THEN
    paramid_out = curr_st_to_grib2 % info(paramid)
    EXIT
  END IF
  curr_st_to_grib2 => curr_st_to_grib2 % next_stash
END DO

find_superppcode_to_paramid = paramid_out

END FUNCTION find_superppcode_to_paramid

FUNCTION find_paramid_to_stcode (paramid_in)

INTEGER :: find_paramid_to_stcode
INTEGER(KIND=KindOfInt) :: paramid_in
INTEGER :: stcode_out

stcode_out = -1

curr_st_to_grib2 => first_st_to_grib2

DO WHILE ( ASSOCIATED(curr_st_to_grib2 ) )
  IF (curr_st_to_grib2 % info(paramid) == paramid_in) THEN
    stcode_out = curr_st_to_grib2 % info(stcode)
    EXIT
  END IF
  curr_st_to_grib2 => curr_st_to_grib2 % next_stash
END DO

find_paramid_to_stcode = stcode_out

END FUNCTION find_paramid_to_stcode

FUNCTION find_stcode_to_level(stcode_in)

REAL :: find_stcode_to_level
INTEGER :: stcode_in
INTEGER :: level_out

level_out = -1

curr_st_to_grib2 => first_st_to_grib2

DO WHILE ( ASSOCIATED(curr_st_to_grib2 ) )
  IF (curr_st_to_grib2 % info(stcode) == stcode_in) THEN
    level_out = curr_st_to_grib2 % info(level)
    EXIT
  END IF
  curr_st_to_grib2 => curr_st_to_grib2 % next_stash
END DO

find_stcode_to_level = TRANSFER(level_out,find_stcode_to_level)

END FUNCTION find_stcode_to_level

FUNCTION find_stcode_to_ppcode(stcode_in)

INTEGER :: find_stcode_to_ppcode
INTEGER :: stcode_in
INTEGER :: ppcode_out

ppcode_out = 0

curr_st_to_grib2 => first_st_to_grib2

DO WHILE ( ASSOCIATED(curr_st_to_grib2 ) )
  IF (curr_st_to_grib2 % info(stcode) == stcode_in) THEN
    ppcode_out = MOD(curr_st_to_grib2 % info(fieldcode),10000)
    ppcode_out = MAX(ppcode_out,0)
    EXIT
  END IF
  curr_st_to_grib2 => curr_st_to_grib2 % next_stash
END DO

find_stcode_to_ppcode = ppcode_out

END FUNCTION find_stcode_to_ppcode

FUNCTION find_stcode_to_lbvc(stcode_in)

INTEGER :: find_stcode_to_lbvc
INTEGER :: stcode_in
INTEGER :: lbvc_out

lbvc_out = 0

curr_st_to_grib2 => first_st_to_grib2

DO WHILE ( ASSOCIATED(curr_st_to_grib2 ) )
  IF (curr_st_to_grib2 % info(stcode) == stcode_in) THEN
    lbvc_out = ABS(curr_st_to_grib2 % info(fieldcode))/10000
    EXIT
  END IF
  curr_st_to_grib2 => curr_st_to_grib2 % next_stash
END DO

find_stcode_to_lbvc = lbvc_out

END FUNCTION find_stcode_to_lbvc

FUNCTION find_paramid_to_manip_val(paramid_in)

REAL :: find_paramid_to_manip_val
INTEGER :: paramid_in
INTEGER :: manip_val_out

manip_val_out = 1.0

curr_st_to_grib2 => first_st_to_grib2

DO WHILE ( ASSOCIATED(curr_st_to_grib2 ) )
  IF (curr_st_to_grib2 % info(paramid) == paramid_in) THEN
    manip_val_out = curr_st_to_grib2 % info(manip_val)
    EXIT
  END IF
  curr_st_to_grib2 => curr_st_to_grib2 % next_stash
END DO

find_paramid_to_manip_val = TRANSFER(manip_val_out,find_paramid_to_manip_val)

END FUNCTION find_paramid_to_manip_val

FUNCTION find_paramid_to_manip_op(paramid_in)

CHARACTER :: find_paramid_to_manip_op
INTEGER :: paramid_in

curr_st_to_grib2 => first_st_to_grib2

DO WHILE ( ASSOCIATED(curr_st_to_grib2 ) )
  IF (curr_st_to_grib2 % info(paramid) == paramid_in) THEN
    find_paramid_to_manip_op = curr_st_to_grib2 % manip_op
    EXIT
  END IF
  curr_st_to_grib2 => curr_st_to_grib2 % next_stash
END DO

END FUNCTION find_paramid_to_manip_op

SUBROUTINE ff_to_grib_data_manip_1D(paramid_in, values)
INTEGER :: paramid_in
REAL    :: values(:)

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

REAL      :: manip_val_out
CHARACTER :: manip_op_out

manip_val_out = find_paramid_to_manip_val(paramid_in)
manip_op_out = find_paramid_to_manip_op(paramid_in)

! Only perfom not unity operations (this is the default)
IF (.NOT. (SCAN(manip_op_out,"*") > 0 .AND. manip_val_out == 1.0) ) THEN 
  ! We now apply the operation with the manip. val
  SELECT CASE(manip_op_out)
  CASE ('*')
    WHERE (values /= RMDI) values = manip_val_out*values
  CASE ('/')
    WHERE (values /= RMDI) values = values/manip_val_out
  CASE ('+')
    WHERE (values /= RMDI) values = values+manip_val_out
  CASE ('-')
    WHERE (values /= RMDI) values = values-manip_val_out
  CASE ('E')
    WHERE (values /= RMDI) values = EXP(values)
  CASE ('L')
    WHERE (values /= RMDI) values = LOG(values)
  END SELECT
END IF

END SUBROUTINE ff_to_grib_data_manip_1D

SUBROUTINE grib_to_ff_data_manip_2D(paramid_in, values)
INTEGER :: paramid_in
REAL :: values(:,:)

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

REAL      :: manip_val_out
CHARACTER :: manip_op_out

manip_val_out = find_paramid_to_manip_val(paramid_in)
manip_op_out = find_paramid_to_manip_op(paramid_in)

! Only perfom not unity operations
IF (.NOT. (SCAN(manip_op_out,"*/") > 0 .AND. manip_val_out == 1.0) ) THEN 
  ! Since we are going from grib to ff we need to reverse the operation.
  
  SELECT CASE(manip_op_out)
  CASE ('*')
    WHERE (values /= RMDI) values = values/manip_val_out
  CASE ('/')
    WHERE (values /= RMDI) values = values*manip_val_out
  CASE ('+')
    WHERE (values /= RMDI) values = values-manip_val_out
  CASE ('-')
    WHERE (values /= RMDI) values = values+manip_val_out
  CASE ('E')
    WHERE (values /= RMDI) values = LOG(values)
  CASE ('L')
    WHERE (values /= RMDI) values = EXP(values)
  END SELECT
END IF

END SUBROUTINE grib_to_ff_data_manip_2D

END MODULE ff_to_grib2_mod
