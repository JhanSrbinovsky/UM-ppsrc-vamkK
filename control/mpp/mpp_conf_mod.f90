! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE mpp_conf_mod

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Purpose:
!    This module contains configuration variables for MPP code.
!
! Code Description:
!   Language: FORTRAN 90

USE global_2d_sums_mod, ONLY : global_sum_method

IMPLICIT NONE

! Variables declared with a sensible default
INTEGER :: extended_halo_size_ew    = 4
INTEGER :: extended_halo_size_ns    = 4
INTEGER :: gcom_coll_limit          = 1   

NAMELIST / NLST_MPP / extended_halo_size_ew,            &
                      extended_halo_size_ns,            &
                      gcom_coll_limit,                  &
                      global_sum_method

END MODULE mpp_conf_mod
