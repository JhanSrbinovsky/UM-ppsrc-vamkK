! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Passes out Ocean LBCs equally among all processors

Module Rcf_Scatter_Ocean_Lbcs_Mod

!  Subroutine Rcf_Scatter_Ocean_LBCS - scatters ocean lbc fields
!
! Description:
!   Scatters ocean lbc field from 1 processor to all - equal amounts
!   go to each processor with any excess going to the last PE.
!
! Method:
!   Uses GC_RSEND and GC_RRECV. Assumes if input or output data will
!   be of same size.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Scatter_Ocean_Lbcs(Full_LBC,     Full_LBC_Size,   &
                                  Part_LBC,     Part_LBC_Size,   &
                                  Stash_Record, Scatter_pe)

USE UM_ParVars, Only : &
    mype,               &
    nproc,              &
    glsize

Use Rcf_Address_Length_Mod, Only : &
    Rcf_Address_Length

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Ppx_Info_Mod, Only : &
    STM_Record_Type

USE gcom_mod, ONLY : gc_ok

IMPLICIT NONE

! Arguments
Integer, Intent(In)                 :: Full_LBC_Size
Integer, Intent(Out)                :: Part_LBC_Size
Integer, Intent(In)                 :: Scatter_pe

Real, Intent(In)                    :: Full_LBC( Full_LBC_Size )
Real, Intent(Out)                   :: Part_LBC( * )

Type( STM_Record_Type ), Intent(In) :: Stash_Record

! Local variables
Integer                      :: i
Integer                      :: global_size
Integer                      :: local_size
Integer                      :: local_extra
Integer                      :: istat
Integer                      :: ErrorStatus
Character (Len=*), Parameter :: RoutineName='Rcf_Scatter_Ocean_LBCs'
Character (Len=80)           :: Cmessage



! Global size of LBC - found from Rcf_Address_Length( input and
! output grids same size )

Call Rcf_Address_Length( Stash_Record % grid_type,             &
                         Stash_Record % halo_type, global_size )

local_size = global_size / nproc
local_extra = Mod( global_size, nproc )

! Use gc_rsend and gc_rrecv - could use gcg_ralltoall but this is
! simpler to work out

! First send from processor scatter_pe
If (mype == scatter_pe) Then
  Do i = 0, nproc - 2
    Call Gc_rsend( 100, local_size, i, istat, Part_LBC, &
                                       Full_LBC( i * local_size + 1 ) )
    If (istat /= GC_OK) Then
      ErrorStatus = 10
      Cmessage = 'Failure in Gcom - gc_rsend'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If
  End Do

  ! Treat nproc - 1 seperately
  i = nproc - 1
  Call Gc_rsend( 100, local_size + local_extra, i, istat, Part_LBC, &
                                       Full_LBC( i * local_size + 1 ) )
  If (istat /= GC_OK) Then
    ErrorStatus = 20
    Cmessage = 'Failure in Gcom - gc_rsend'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If
End If

If (mype /= nproc - 1) Then
  Part_LBC_size = local_size
Else
  Part_LBC_size = local_size + local_extra
End If

! Now receive from Processor scatter_pe on all pes
Call Gc_rrecv( 100, Part_LBC_size, scatter_pe, istat, Part_LBC,    &
                                      Full_LBC( mype * local_size + 1) )
If (istat /= GC_OK) Then
  ErrorStatus = 30
  Cmessage = 'Failure in Gcom - gc_rrecv'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If


Return
End Subroutine Rcf_Scatter_Ocean_Lbcs

End Module Rcf_Scatter_Ocean_Lbcs_Mod


