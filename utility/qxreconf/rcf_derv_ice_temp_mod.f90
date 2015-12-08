! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reinitialises sea ice surface temperature

MODULE Rcf_Derv_Ice_Temp_Mod

! Description:
!     This subroutine is the initialises the sea ice temperature
!     on ice catagories when these are called for in stash
!     It sets all ice category pseudo levels to have temperatures
!     equal to the single level sea ice temperature (stash=49)
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

CONTAINS

SUBROUTINE Rcf_Derv_Ice_Temp( fields_out, field_count_out, hdr_out,  &
                            ice_temp_cat )

USE Rcf_Locate_Mod, ONLY : &
    Rcf_Locate

USE Rcf_Alloc_Field_Mod, ONLY : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Stashcodes_Mod, ONLY : &
    stashcode_sea_ice_temp,    &
    stashcode_prog_sec

USE Rcf_Field_Type_Mod, ONLY : &
    Field_Type

USE Rcf_UMhead_Mod, ONLY : &
    um_header_type

USE PrintStatus_mod, ONLY : &
    PrintStatus,                &
    PrStatus_Normal

USE UM_ParVars, ONLY : &
    mype

USE decomp_params, ONLY : &
    decomp_rcf_output

USE Rcf_Read_Field_Mod, ONLY : &
    Rcf_Read_Field


IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( um_header_type), INTENT(In) :: hdr_out
TYPE( field_type ), INTENT(InOut) :: ice_temp_cat
INTEGER, INTENT(In)               :: field_count_out

! Internal variables
TYPE( field_type ), POINTER       ::  ice_temp


INTEGER                           ::  pos   ! position in array
INTEGER                           ::  i     ! loop index



!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
    WRITE (6,*) 'Reinitialising ice surface temperature'
  END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in
!----------------------------------------------------------------------
! GBM ice temperature; will abort if ice_temp not found

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_sea_ice_temp,       &
                   fields_out,field_count_out,pos)

  ice_temp => fields_out(pos)
  CALL Rcf_Alloc_Field( ice_temp )
  CALL Rcf_Read_Field( ice_temp, hdr_out, decomp_rcf_output )


!----------------------------------------------------------------------
! Loop through ice_temp_cat
!----------------------------------------------------------------------

  DO i = 1,ice_temp_cat % levels
    ice_temp_cat % DATA(:,i) = ice_temp % DATA(:,1)
  ENDDO

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------
  CALL Rcf_Dealloc_Field( ice_temp )

RETURN
END SUBROUTINE Rcf_Derv_Ice_Temp
END MODULE Rcf_Derv_Ice_Temp_Mod
