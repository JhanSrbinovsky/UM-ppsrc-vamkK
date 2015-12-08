! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ broadcasts a STASHmaster record from 1 to many PEs

Module Rcf_bcastSTM_Mod

!  Subroutine Rcf_bcastSTM - broadcasts a STM record
!
! Description:
!   Used to broadcast a STASHmaster record to all PEs for PE source.
!   An assumption is made that the STASHmaster records are
!   contiguous in memory (ie sequenced or COMMON).
!
! Method:
!   Two calls to GCOM are used - 1 for integer and 1 for character
!   data. Sizes need to change if STASHmaster format is changed
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.5   10/02/03   Extension of option code to 30 digits. T. White
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains

Subroutine Rcf_BcastSTM( STASHmaster , source )

Use Rcf_Ppx_Info_Mod, Only : &
    Ppxref_pack_profs,   &
    STM_IntDataLen, &
    STM_CharDataLen, &
    STM_record_type

USE UM_ParVars, Only :  &
    mype,                &
    nproc

Implicit None

! Arguments
Integer                :: source      ! PE that holds source record
Type (STM_record_type) :: STASHmaster ! The record to broadcast

! Local variables
Integer                :: msg         ! message tag
Integer                :: info        ! dummy

!------------------------------------------------------------------
! Need 2 broadcasts to do whole thing!
!------------------------------------------------------------------

msg = 8001
Call gc_ibcast( msg, STM_IntDataLen, source, nproc, info, &
                STASHmaster % model )

msg = 8002
Call gc_cbcast( msg, STM_CharDataLen, source, nproc, info, &
                STASHmaster % name )

Return
End Subroutine Rcf_BcastSTM

End Module Rcf_BcastSTM_Mod
