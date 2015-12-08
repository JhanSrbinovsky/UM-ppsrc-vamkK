! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lw_diag_mod

  USE max_calls
  USE lwrdiag_mod

  IMPLICIT NONE

!     Declaration of LW diagnostics.
!     Those quantities which are purely diagnostic are bundled into
!     a structure.

!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: Radiation Control

      TYPE (strlwdiag), SAVE :: lw_diag(npd_lwcall)

END MODULE lw_diag_mod
