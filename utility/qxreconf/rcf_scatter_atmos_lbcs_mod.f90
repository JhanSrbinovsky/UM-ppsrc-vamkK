! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Passes out atmosphere LBCs equally among all processors

MODULE Rcf_Scatter_Atmos_Lbcs_Mod

!  Subroutine Rcf_Scatter_Atmos_LBCs - scatters atmos lbcs
!
! Description:
!   Scatters atmos LBCs evenly to all pes from a single PE which
!   holds whole field at routine start - excess goes to last PE.
!
! Method:
!   Local sizes are calculated and data is scattered using GC_RSEND.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

CONTAINS

SUBROUTINE Rcf_Scatter_Atmos_Lbcs(Full_LBC,     Full_LBC_Size,   &
                                  Part_LBC,     Part_LBC_Size,   &
                                  Stash_Record, Scatter_pe)

USE UM_ParVars, ONLY : &
    mype,               &
    nproc,              &
    glsize

USE Ereport_Mod, ONLY : &
    Ereport

USE Rcf_Ppx_Info_Mod, ONLY : &
    STM_Record_Type

USE Rcf_Level_Code_Mod, ONLY : &
    Rcf_Level_Code

USE Rcf_Address_Length_Mod, ONLY : &
    Rcf_Address_length

USE Rcf_Grid_Type_Mod, ONLY : &
    Input_Grid

USE gcom_mod, ONLY : gc_ok

IMPLICIT NONE

! Arguments
INTEGER, INTENT(In)                 :: Full_LBC_Size
INTEGER, INTENT(Out)                :: Part_LBC_Size
INTEGER, INTENT(In)                 :: Scatter_pe

REAL, INTENT(In)                    :: Full_LBC( Full_LBC_Size )
REAL, INTENT(Out)                   :: Part_LBC( * )

TYPE( STM_Record_Type ), INTENT(In) :: Stash_Record

! Local variables
INTEGER                      :: i
INTEGER                      :: global_size
INTEGER                      :: local_size
INTEGER                      :: local_extra
INTEGER                      :: start_level
INTEGER                      :: end_level
INTEGER                      :: istat
INTEGER                      :: ErrorStatus
CHARACTER (Len=*), PARAMETER :: RoutineName='Rcf_Scatter_Atmos_LBCs'
CHARACTER (Len=80)           :: Cmessage


! Global size of LBC - found from Rcf_Address_Length( input and
! output grids same size )

CALL Rcf_Address_Length( Stash_Record % Grid_Type,             &
                         Stash_Record % Halo_Type, global_size )

! Need to multiply this size by the number of levels
CALL Rcf_Level_Code( Stash_Record % lb_code, start_level,  Input_Grid )
CALL Rcf_Level_Code( Stash_Record % lt_code, end_level,    Input_Grid )

global_size = global_size * ( (end_level - start_level) + 1)


! Work out how much data per PE now.

local_size = global_size / nproc
local_extra = MOD( global_size, nproc )

! Use gc_rsend and gc_rrecv - could use gcg_ralltoall but this is
! simpler to work out

! First send from processor scatter_pe
IF (mype == scatter_pe) THEN
  DO i = 0, nproc - 2
    CALL Gc_rsend( 100, local_size, i, istat, Part_LBC, &
                                       Full_LBC( i * local_size + 1 ) )
    IF (istat /= GC_OK) THEN
      ErrorStatus = 10
      Cmessage = 'Failure in Gcom - gc_rsend'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF
  END DO

  ! Treat nproc - 1 seperately
  i = nproc - 1
  CALL Gc_rsend( 100, local_size + local_extra, i, istat, Part_LBC, &
                                       Full_LBC( i * local_size + 1 ) )
  IF (istat /= GC_OK) THEN
    ErrorStatus = 20
    Cmessage = 'Failure in Gcom - gc_rsend'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF
END IF

IF (mype /= nproc - 1) THEN
  Part_LBC_size = local_size
ELSE
  Part_LBC_size = local_size + local_extra
END IF

! Now receive from Processor scatter_pe on all pes
CALL Gc_rrecv( 100, Part_LBC_size, scatter_pe, istat, Part_LBC,    &
                                      Full_LBC( mype * local_size + 1) )
IF (istat /= GC_OK) THEN
  ErrorStatus = 30
  Cmessage = 'Failure in Gcom - gc_rrecv'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

RETURN
END SUBROUTINE Rcf_Scatter_Atmos_Lbcs

END MODULE Rcf_Scatter_Atmos_Lbcs_Mod


