! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets the data_source array corresponding to the fields array.

MODULE rcf_locate_alt_field_mod
!  Subroutine rcf_locate_alt_field - sets field to use alternative.
!
! Description:
!   When the field is not in the input dump sometimes we can fallback to using
!   an alternative.  This is similar to fieldcalcs but we may want to
!   interpolate and perform other types of processing on some fields which
!   fieldcalcs does not do.
!
! Method:
!   Using predefined list we can fallback to using other fields if the field
!   is not available in the input dump.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


CONTAINS
SUBROUTINE rcf_locate_alt_field( field, fields, fields_count, pos)

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Field_Type_Mod, Only : &
    Field_type

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Submodel_Mod, Only :   &
    Atmos_IM

Use Rcf_Stashcodes_Mod

Use UM_Parvars, Only : &
    mype

Use PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_Ppx_Info_Mod, Only : &
    STM_record_type

IMPLICIT NONE

! Arguments
TYPE (field_type)                         :: field
TYPE (field_type), POINTER                :: fields(:)
INTEGER                                   :: fields_count
INTEGER                                   :: pos
CHARACTER(LEN=80)                         :: cmessage

! Local variables
CHARACTER(LEN=*), PARAMETER :: routinename = "rcf_read_other_field"
INTEGER :: other_sec
INTEGER :: other_item
INTEGER :: errorstatus
TYPE(STM_record_type), POINTER :: other_stm

other_sec = -1
other_item = -1
NULLIFY(other_stm)

! Lets see what field we want
SELECT CASE( field % stashmaster % model )
CASE ( Atmos_IM )
  SELECT CASE( field % stashmaster % section )
  CASE( stashcode_prog_sec )
    SELECT CASE( field % stashmaster % item )
    CASE( stashcode_theta )
      ! We keep the STM as stashcode_t.
      other_sec  = stashcode_t/1000
      other_item = MOD(stashcode_t, 1000)
      If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
        Write (6,*) "Using ", stashcode_t, " for ", field % stashmaster % name
      End If
      
    END SELECT
  END SELECT
END SELECT

IF (other_sec == -1 .OR. other_item == -1) THEN
  WRITE(cmessage,*) "Alternative field not found for ", &
                    field % stashmaster % name
  errorstatus = 10
  CALL ereport(RoutineName, errorstatus, cmessage)
END IF

! Now lets find the field.
CALL rcf_locate(other_sec, other_item, fields, fields_count, pos)

! Reset the STASHmaster if defined.
IF (ASSOCIATED(other_stm)) THEN
  fields(pos) % stashmaster => other_stm
END IF

RETURN
END SUBROUTINE rcf_locate_alt_field
END MODULE rcf_locate_alt_field_mod
