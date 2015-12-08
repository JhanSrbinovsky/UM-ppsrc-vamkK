! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set the precision of real variables.

MODULE realtype_rd

! Description:
!   There are two versions of this routine. Here in the UM version we
!   set the precision used in the radiance_core routines to be real64.
!   There is a corresponding version for the offline radiation code 
!   where RealK is set explicitly.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

  USE um_types, ONLY: real64

  INTEGER, PARAMETER :: RealK=real64

END MODULE realtype_rd
