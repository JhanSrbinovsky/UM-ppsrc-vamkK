! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
MODULE o3intp_mod

IMPLICIT NONE
!+ ---------------------------------------------------------------------
!  Module to specify allowed methods of interpolating from the
!  ancillary file.
!- ---------------------------------------------------------------------
!
      Integer, parameter :: IO3_3DSPEC = 1
!       Ozone is provided as a full 3D field.
      Integer, parameter :: IO3_2DSPEC = 2
!       Ozone is expanded from a 2D field by direct copying.
      Integer, parameter :: IO3_2DMASSCON = 3
!       Ozone is expanded from a 2D field with conservation of mass.
      Integer, parameter :: IO3_TROP_MAP = 4
!       Ozone mixing ratios are set by mapping each height at each
!       grid-point in the real profile to a height in the ancillary
!       profile and using the mixing ratio there. The mapping is
!       set using the height of the tropopause
      Integer, parameter :: IO3_TROP_MAP_MASSCON = 5
!       Ozone mixing ratios are set as above, but scaled so as to
!       preserve the vertically integrated column ozone
!
! ----------------------------------------------------------------------

END MODULE o3intp_mod
