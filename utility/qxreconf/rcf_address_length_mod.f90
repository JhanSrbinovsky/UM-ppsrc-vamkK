! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+Determine STASH input length per vertical level for prog var

Module Rcf_Address_Length_Mod

!  Subroutine Rcf_Address_Length - determines field size
!
! Description:
!    Calculates size of field for output dump addressing (single level).
!
! Method:
!    Calculates field sizes based on Grid_Type from stashmaster
!    Generally single level, but old LBC needs level info.
!
!    Based on UM 4.5/5.0 code.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains

SUBROUTINE Rcf_Address_Length( IGPL, halo_type, LEN )

Use rimtypes
Use lbc_mod
Use Submodel_Mod

Use Rcf_CntlAtm_Mod

Use Rcf_Lsm_Mod, Only : &
    glob_land_out

Use Rcf_Model_Mod, Only :  &
    ZonAvOzone

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

USE UM_ParVars

USE nlsizes_namelist_mod, ONLY: &
    tr_vars

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

IMPLICIT NONE

! Subroutine arguments:
Integer, Intent(In)   :: IGPL      ! Grid code
Integer, Intent(In)   :: halo_type ! code from stashmaster
Integer, Intent(Out)  :: LEN       ! Length

! Local scalars:
INTEGER               :: IP   !  pressure levels
INTEGER               :: IQ   !  wet levels
INTEGER               :: IT   !  tracer levels
INTEGER               :: IX1  !  leftmost point
INTEGER               :: IX2  !  rightmost point
INTEGER               :: IY1  !  lower point
INTEGER               :: IY2  ! upper point


! Function and subroutine calls:
EXTERNAL LLTORC

!- End of Header ---------------------------------------------------

! Determine row/column nos. for global domain on output grid
! DEPENDS ON: lltorc
CALL LLTORC(IGPL,90,-90,0,360,IY1,IY2,IX1,IX2, Output_Grid)


Select Case (IGPL)
  Case ( 21 )
    LEN= glob_land_out     ! Land compressed

  Case ( 22 )
    IF (ZonAvOzone) THEN   !  Zonal
      LEN = IY2-IY1+1
    ELSE                   !  Full fields
      LEN = (IX2-IX1+1)*(IY2-IY1+1)
    ENDIF

  Case (23)
    LEN = Output_Grid % glob_r_field

  Case ( 25 )              ! Old LBC
    CALL Rcf_Level_Code( 2, IP, Output_Grid)   ! pressure levels
    CALL Rcf_Level_Code( 3, IQ, Output_Grid)   ! wet levels
    CALL Rcf_Level_Code(11, IT, Output_Grid)   ! tracer levels
    LEN=( Output_Grid % GLOB_P_ROWS +                                &
          Output_Grid % GLOB_P_ROW_LENGTH                            &
          -2*RIMWIDTHA(rima_type_norm))*2*RIMWIDTHA(rima_type_norm)* &
        (1+3*IP+2*IQ+TR_VARS*IT)-2*IP*4*RIMWIDTHA(rima_type_norm)


  Case Default
    LEN =(IX2-IX1+1)*(IY2-IY1+1)

End Select

RETURN
END Subroutine Rcf_Address_Length

End Module Rcf_Address_Length_Mod
