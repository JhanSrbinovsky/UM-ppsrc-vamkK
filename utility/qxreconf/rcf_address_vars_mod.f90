! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Storage of Addressing variables

Module Rcf_Address_Vars_Mod

! Description:
!   Module specifically for variable storage. Variables used in
!   RCF_Address and TSTMSK.
!
! Method:
!   Declare a number of variables and USE them from the appropriate
!   locations
!
!   Original variables stripped from rcf_address.F90 to prevent circular
!   usage of modules between rcf_address and tstmsk. FCM didn't like
!   the pre-exisiting arrangement.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Use Rcf_Ppx_Info_Mod, Only :    &
    STM_OptCodeLen

IMPLICIT NONE

Integer, Save   ::    ISPACE         ! Space code
Integer, Save   ::    IGP            ! Grid of data code
Integer, Save   ::    ILEV           ! Level type code
Integer, Save   ::    IBOT           ! First level code
Integer, Save   ::    ITOP           ! Last level code
Integer, Save   ::    IFLAG          ! Level compression flag
Integer, Save   ::    IOPN(STM_OptCodeLen/5) ! Sectional option code
Integer, Save   ::    VMSK           ! Integer equiv of bin vers mask
Integer, Save   ::    IPSEUDO        ! Pseudo dimension type
Integer, Save   ::    IPFIRST        ! First pseudo dim code
Integer, Save   ::    IPLAST         ! Last pseudo dim code
Integer, Save   ::    HALO           ! Halo type code

End Module Rcf_Address_Vars_Mod
