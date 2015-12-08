! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ SCM module for making SCMop accessible to SCMoutput from SCM_main.

MODULE global_scmop
  USE UM_types
! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE
!
! Description:
!   All the SCM diagnostic information is stored in a single
!   structure, SCMop, first declared in SCM_main. This module acts
!   like a common block to allow both SCM_main and SCMoutput access
!   to this structure without it having to be passed down the call tree.
!   A module has to be used because a derived-type structure cannot be
!   put in a common.
!
! Original Owner: Luke Jones
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
! Language: Fortran 90.

! Declare SCMop
  TYPE(SCMop_type) :: SCMop

END MODULE global_scmop
