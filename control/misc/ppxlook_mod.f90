! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Misc

! Purpose: Declares ppxref look-up arrays used by the UM and associated
!           arrays and parameters.

MODULE ppxlook_mod

USE version_mod, ONLY: nitemp, nsectp, ndiagp
USE cppxref_mod, ONLY: ppxref_codelen, ppxref_charlen
USE submodel_mod, ONLY: n_internal_model_max

IMPLICIT NONE

PRIVATE n_internal_model_max

! No. of STASH items per section
INTEGER, PARAMETER :: ppxref_items = nitemp

! No. of STASH sections per internal model
INTEGER, PARAMETER :: ppxref_sections = nsectp-49

! Max. number of non-null records in ppxref file (>1200)
INTEGER, PARAMETER :: num_diag_max = ndiagp

! Max. number of user-defined ppxref records allowed
INTEGER, PARAMETER :: num_usr_diag_max = 850

! No. of ppxref records read into PPXI,PPXC (for dyn. allocation)

INTEGER :: ppxrecs


! Global arrays:
! ppxref look-up arrays



INTEGER, DIMENSION(:,:), ALLOCATABLE :: ppxi


CHARACTER :: ppxc(num_diag_max,ppxref_charlen)


! Arrays for temporary storage of user-ppxref records -
!   used to transfer these records from STASH_PROC into U_MODEL
INTEGER :: ppxi_u(num_usr_diag_max,ppxref_codelen)
CHARACTER :: ppxc_u(num_usr_diag_max,ppxref_codelen)

! Array of flags to indicate origin of ppxref record
! 'P' for ppxref file; 'U' for user-stash master file
CHARACTER :: OriginFlag(num_diag_max)

! Array of indices to identify which ppxref record corresponds to
! any given row of PPXI, PPXC
INTEGER :: RowIndex(num_diag_max)

! Pointer array for PPXI, PPXC arrays



INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ppxptr


END MODULE ppxlook_mod
