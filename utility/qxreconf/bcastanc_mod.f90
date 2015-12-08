! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
Module bcastANC_Mod

! This is used to broadcast a ANCILmaster record
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

Contains

Subroutine BcastANC_Fields( ANCILmaster , source )

Use Ancil_Mod, Only : &
    ANC_record_type,  &
    AncRec_IntDataLen,   &
    AncRec_CharDataLen

USE UM_ParVars, Only :  &
    mype,                &
    nproc

Implicit None

! Arguments
Integer                :: source      ! PE that holds source record
Type (ANC_record_type) :: ANCILmaster ! The record to broadcast

! Local variables
Integer                :: msg         ! message tag
Integer                :: info        ! dummy

!------------------------------------------------------------------
! Need 2 broadcasts : First for Integers, Second for Characters
!------------------------------------------------------------------

msg = 8051
Call gc_ibcast( msg, AncRec_IntDataLen, source, nproc, info,  &
                ANCILmaster % ancil_ref_number )

msg = 8052
Call gc_cbcast( msg, AncRec_CharDataLen, source, nproc, info, &
                ANCILmaster % anc_name )

Return
End Subroutine BcastANC_Fields

Subroutine BcastANC_Files ( ANCILmaster , source )

Use Ancil_Mod, Only :   &
    ANC_file_type,      &
    AncFile_IntDataLen, &
    AncFile_CharDataLen

USE UM_ParVars, Only :  &
    mype,                    &
    nproc

Implicit None

! Arguments
Integer              :: source      ! PE that holds source record
Type (ANC_file_type) :: ANCILmaster ! The record to broadcast

! Local variables
Integer              :: msg         ! message tag
Integer              :: info        ! dummy

!------------------------------------------------------------------
! Need 2 broadcasts : First for Integers, Second for Characters
!------------------------------------------------------------------

msg = 8053
Call gc_ibcast( msg, AncFile_IntDataLen, source, nproc, info,  &
                ANCILmaster % anc_file_number )

msg = 8054
Call gc_cbcast( msg, AncFile_CharDataLen, source, nproc, info, &
                ANCILmaster % anc_env_var )

Return
End Subroutine BcastANC_Files


End Module BcastANC_Mod
