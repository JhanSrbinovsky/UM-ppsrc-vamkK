! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers atmosphere LBCs from all processors.

MODULE Rcf_Gather_Atmos_Lbcs_Mod

!  Subroutine Rcf_Gather_Atmos_LBCs
!
! Description:
!   Gathers a distributed atmosphere LBC field onto a single PE.
!
! Method:
!   Calculates remote sizes and gathers using GC_RSEND/RRECV.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

CONTAINS

SUBROUTINE Rcf_Gather_Atmos_Lbcs(Full_LBC,     Full_LBC_Size,   &
                                 Part_LBC,     Part_LBC_Size,   &
                                 Stash_Record, Gather_pe)

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
    Rcf_Address_Length

USE Rcf_Grid_Type_Mod, ONLY : &
    Output_Grid

USE gcom_mod, ONLY : gc_ok

IMPLICIT NONE

! Arguments
INTEGER, INTENT(In)                 :: Full_LBC_Size
INTEGER, INTENT(Out)                :: Part_LBC_Size
INTEGER, INTENT(In)                 :: Gather_pe

REAL, INTENT(Out)                   :: Full_LBC( Full_LBC_Size )
REAL, INTENT(In)                    :: Part_LBC( * )

TYPE( STM_Record_Type ), INTENT(In) :: Stash_Record

! Local variables
INTEGER                      :: i
INTEGER                      :: global_size
INTEGER                      :: local_size
INTEGER                      :: local_extra
INTEGER                      :: part_size
INTEGER                      :: start_level
INTEGER                      :: end_level
INTEGER                      :: istat
INTEGER                      :: ErrorStatus
CHARACTER (Len=*), PARAMETER :: RoutineName='Rcf_Gather_Atmos_LBCs'
CHARACTER (Len=80)           :: Cmessage


! Global size of LBC - found from Rcf_Address_Length( input and
! output grids same size )

CALL Rcf_Address_Length( Stash_Record % Grid_Type,             &
                         Stash_Record % Halo_Type, global_size )

! Need to multiply this size by the number of levels
CALL Rcf_Level_Code( Stash_Record % lb_code, start_level,  Output_Grid )
CALL Rcf_Level_Code( Stash_Record % lt_code, end_level,    Output_Grid )

global_size = global_size * ( (end_level - start_level) + 1)


! Work out how much data per PE now.

local_size = global_size / nproc
local_extra = MOD( global_size, nproc )

! Use gc_rsend and gc_rrecv - could use gcg_ralltoall but this is
! simpler to work out

IF (mype /= nproc - 1) THEN
  Part_size = local_size
ELSE
  Part_size = local_size + local_extra
END IF

! Now send to Processor gather_pe on all pes
CALL Gc_rsend( 200, Part_size, gather_pe, istat, &
               Full_LBC( mype * local_size + 1), Part_LBC )
IF (istat /= GC_OK) THEN
  ErrorStatus = 30
  Cmessage = 'Failure in Gcom - gc_rsend'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Now recv on processor gather_pe
IF (mype == gather_pe) THEN
  DO i = 0, nproc - 2
    CALL Gc_rrecv( 200, local_size, i, istat, &
                   Full_LBC( i * local_size + 1 ), Part_LBC )
    IF (istat /= GC_OK) THEN
      ErrorStatus = 10
      Cmessage = 'Failure in Gcom - gc_rrecv'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF
  END DO

  ! Treat nproc - 1 seperately
  i = nproc - 1
  CALL Gc_rrecv( 200, local_size + local_extra, i, istat, &
                 Full_LBC( i * local_size + 1 ), Part_LBC )
  IF (istat /= GC_OK) THEN
    ErrorStatus = 20
    Cmessage = 'Failure in Gcom - gc_rrecv'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

END IF

! Set the return size
part_LBC_size = local_size
IF (mype == nproc - 1) part_LBC_size = part_LBC_size + local_extra

RETURN
END SUBROUTINE Rcf_Gather_Atmos_Lbcs

END MODULE Rcf_Gather_Atmos_Lbcs_Mod


