! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Wrapper for the READFLDS routine.

MODULE Rcf_Read_Field_Mod
IMPLICIT NONE

!  Subroutine Rcf_Read_Field - read in a field
!
! Description:
!   A wrappr for the READFLDS routine - thus reads a field into the
!   field data-type
!
! Method:
!   Selects real,int or logical data pointer into which to read data.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

CONTAINS

SUBROUTINE Rcf_Read_Field( field, hdr, decomp )

USE Rcf_Field_Type_Mod, ONLY : &
    field_type

USE Ereport_mod, ONLY : &
    Ereport

USE Rcf_UMhead_Mod, ONLY : &
    um_header_type

USE PrintStatus_mod, ONLY : &
    LTimer

USE UM_ParVars, ONLY :   &
    current_decomp_type, &
    change_decomposition

USE cppxref_mod, ONLY: &
    ppx_type_real,     &
    ppx_type_int,      &
    ppx_type_log

IMPLICIT NONE

! Arguments
TYPE( field_type ), INTENT( InOut )    :: field
TYPE( um_header_type), INTENT( In )    :: hdr
INTEGER, INTENT( In )                  :: decomp

! Local vars
CHARACTER (Len=*), PARAMETER :: RoutineName='Rcf_Read_Field'
CHARACTER (Len=80)           :: Cmessage
INTEGER                      :: ErrorStatus
INTEGER                      :: orig_decomp     ! decomposition on
                                                ! entry to routine

EXTERNAL Timer
!----------------------------------------------------------------

! DEPENDS ON: timer
IF (LTimer) CALL Timer( RoutineName, 3)
!------------------------------------------------------------------
! Change decomposition if required not same as current
!------------------------------------------------------------------
orig_decomp = current_decomp_type
IF (orig_decomp /= decomp) THEN
  CALL Change_Decomposition( decomp )
END IF

!------------------------------------------------------------------
! Read in the data - need 3 cases real, logical, integer
!------------------------------------------------------------------
SELECT CASE( field % stashmaster % data_type )
  CASE (ppx_type_real)
! DEPENDS ON: rcf_readflds
    CALL Rcf_ReadFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                       hdr % Lookup, hdr % Len1Lookup, field % DATA,   &
                       field % level_size, hdr % FixHd,                &
                       ErrorStatus, Cmessage )

  CASE (ppx_type_int)
! DEPENDS ON: rcf_readflds
    CALL Rcf_ReadFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                      hdr % Lookup, hdr % Len1Lookup, field % Data_Int,&
                       field % level_size, hdr % FixHd,                &
                       ErrorStatus, Cmessage )

  CASE (ppx_type_log)
! DEPENDS ON: rcf_readflds
    CALL Rcf_ReadFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                      hdr % Lookup, hdr % Len1Lookup, field % Data_Log,&
                       field % level_size, hdr % FixHd,                &
                       ErrorStatus, Cmessage )

  CASE Default
    ErrorStatus = -10
    Cmessage = 'Unable to read in field as datatype cannot be '//&
               'determined'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

!-------------------------------------------------------------------
! Check the returned error conditions
!-------------------------------------------------------------------
IF ( ErrorStatus /= 0 ) THEN
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!--------------------------------------------------------------------
! Change the decomposition back if required
!---------------------------------------------------------------------
IF (current_decomp_type /= orig_decomp) THEN
  CALL Change_Decomposition( orig_decomp )
END IF

! DEPENDS ON: timer
IF (LTimer) CALL Timer( RoutineName, 4)

RETURN

END SUBROUTINE Rcf_Read_Field
END MODULE Rcf_Read_Field_Mod
