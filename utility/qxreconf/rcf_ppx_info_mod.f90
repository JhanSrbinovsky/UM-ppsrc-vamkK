! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Defines the STASHmaster record

Module Rcf_Ppx_Info_Mod

! Description:
!   Defines the STASHmaster record format, the USTSNUM namelist
!   and other related variables
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Use Submodel_Mod, Only  :     &
    N_Internal_Model_Max

Use ppxlook_mod, Only       : &
    ppxref_sections,          &
    ppxref_items

Implicit None

! A parameter needed for the following definition.
! Number of packing codes in stashmaster record.
Integer, Parameter        :: Ppxref_pack_profs = 10

! Length of Option and Version Codes.
Integer, Parameter        :: STM_OptCodeLen = 30 !Must be multiple of 5
Integer, Parameter        :: STM_VerCodeLen = 20

! For broadcast: quantity of integer/character data in each STMrecord.
! Must be changed if STMrecord format changes.
Integer, Parameter        :: STM_IntDataLen  = 28 + Ppxref_pack_profs &
                                             + STM_OptCodeLen / 5
Integer, Parameter        :: STM_CharDataLen = 37

! Define STMrecord type - to hold STASHmaster records
! Note that as these are in a sequence and are integer followed by
! character, we can do some fairly efficient comms with them.
! Thus order *is* vital!
Type STM_record_type
  Sequence
  Integer                 :: model
  Integer                 :: section
  Integer                 :: item
  Integer                 :: space_code
  Integer                 :: ptr_code
  Integer                 :: timavail_code
  Integer                 :: grid_type
  Integer                 :: lv_code
  Integer                 :: lb_code
  Integer                 :: lt_code
  Integer                 :: pt_code
  Integer                 :: pf_code
  Integer                 :: pl_code
  Integer                 :: lev_flag
  Integer                 :: opt_code( STM_OptCodeLen / 5 )
  Integer                 :: version_mask
  Integer                 :: halo_type
  Integer                 :: data_type
  Integer                 :: dump_packing
  Integer                 :: packing_acc(Ppxref_pack_profs)
  Integer                 :: rotate_code
  Integer                 :: field_code
  Integer                 :: user_code
  Integer                 :: lbvc_code
  Integer                 :: base_level
  Integer                 :: top_level
  Integer                 :: ref_lbvc_code
  Integer                 :: cf_levelcode
  Integer                 :: cf_fieldcode
  Integer                 :: RowIndex

  Character (Len=36)      :: name
  Character (Len=1)       :: OriginFlag
End Type STM_record_type

!-------------------------------------------------------------------
! Guts for the storage of STASHmaster information
!-------------------------------------------------------------------
Integer, Target, Save    :: ppxRecs = 0     ! No. of stash records
Integer, Target, Save    :: ppxRecs_in = 0     ! No. of stash records input

Type (STM_record_type), Pointer, Save   :: STM_record(:)
Type (STM_record_type), Pointer, Save   :: STM_record_in(:)

! Referencing for the above
Integer, Target, Save    :: ppxptr( N_Internal_Model_Max,              &
                                 0:Ppxref_Sections , Ppxref_Items ) = 0
Integer, Target, Save    :: ppxptr_in( N_Internal_Model_Max,           &
                                 0:Ppxref_Sections , Ppxref_Items ) = 0

! Info for user STASHmaster file
Integer, Target, Save    :: n_ustash        = 0
Integer, Target, Save    :: n_ustash_in     = 0
Integer, Target, Save    :: nrecs_ustash    = 0
Integer, Target, Save    :: nrecs_ustash_in = 0
Character (Len=80), Target, Save :: ustsfils(20)
Character (Len=80), Target, Save :: ustsfils_in(20)

Namelist/ustsnum /n_ustash, nrecs_ustash, ustsfils
Namelist/ustsnum_in /n_ustash_in, nrecs_ustash_in, ustsfils_in

End Module Rcf_Ppx_Info_mod
