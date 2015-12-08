! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Wrapper for WRITFLDS

MODULE Rcf_Write_Field_Mod
IMPLICIT NONE

!  Subroutine Rcf_Write_Field
!
! Description:
!    Wrapper for writflds - utilising fields data types
!
! Method:
!    Deals seperately with real, integer and logical data.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!

CONTAINS

SUBROUTINE Rcf_Write_Field( field, hdr, decomp, free )


USE Rcf_Alloc_Field_mod, ONLY : &
    Rcf_Dealloc_Field

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

USE cppxref_mod, ONLY : &
    ppx_type_real,      &
    ppx_type_int,       &
    ppx_type_log

IMPLICIT NONE

! Arguments
TYPE( field_type ), INTENT( InOut )    :: field
TYPE( um_header_type), INTENT( In )    :: hdr
INTEGER, INTENT( In )                  :: decomp
LOGICAL, OPTIONAL, INTENT( In )        :: free

! Local vars
CHARACTER (Len=*), PARAMETER :: RoutineName='Rcf_Write_Field'
CHARACTER (Len=80)           :: Cmessage
INTEGER                      :: ErrorStatus
LOGICAL                      :: free_data
INTEGER                      :: orig_decomp ! decomposition when
                                            ! routine is entered

EXTERNAL Rcf_WritFlds, Timer

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
! First set the free_data logical appropriately
!------------------------------------------------------------------
IF ( PRESENT( free ) ) THEN
  free_data = free
ELSE
  free_data = .FALSE.
END IF

!------------------------------------------------------------------
! Now write out the data - need 3 cases real, logical, integer
!------------------------------------------------------------------
SELECT CASE( field % stashmaster % data_type )
  CASE (ppx_type_real)
! DEPENDS ON: rcf_writflds
    CALL Rcf_WritFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                   hdr % Lookup, hdr % Len1Lookup, field % DATA,       &
                   field % level_size, hdr % FixHd,                    &
                   ErrorStatus, Cmessage , .TRUE., .FALSE.)

  CASE (ppx_type_int)
! DEPENDS ON: rcf_writflds
    CALL Rcf_WritFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                   hdr % Lookup, hdr % Len1Lookup, field % Data_Int,   &
                   field % level_size, hdr % FixHd,                    &
                   ErrorStatus, Cmessage , .TRUE., .FALSE.)

  CASE (ppx_type_log)
! DEPENDS ON: rcf_writflds
    CALL Rcf_WritFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                   hdr % Lookup, hdr % Len1Lookup, field % Data_Log,   &
                   field % level_size, hdr % FixHd,                    &
                   ErrorStatus, Cmessage , .TRUE., .FALSE.)

  CASE Default
    ErrorStatus = -10
    Cmessage = 'Unable to write out field as datatype cannot be '//&
        'determined'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

!-------------------------------------------------------------------
! Check the returned error conditions
!-------------------------------------------------------------------
IF ( ErrorStatus /= 0 ) THEN
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!-------------------------------------------------------------------
! Free the allocated data-space if required to
!-------------------------------------------------------------------
IF (free_data) THEN
  CALL Rcf_Dealloc_Field( field )
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

END SUBROUTINE Rcf_Write_Field
END MODULE Rcf_Write_Field_Mod
