! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE sw_diag_mod

  USE max_calls
  USE swrdiag_mod

  IMPLICIT NONE

!   Declaration of SW diagnostics.

!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Radiation Control

!   Those calculated within the radiation code and not used to
!   advance the integration in any way are contained within
!   a structure defined in swrdiag_mod.

    TYPE (strswdiag), SAVE :: sw_diag(npd_swcall)

END MODULE sw_diag_mod
