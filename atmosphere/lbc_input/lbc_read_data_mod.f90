! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for LBC reading

! Description:
!   Module containing data concerning LBC reading
!   
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC input

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE lbc_read_data_mod

IMPLICIT NONE

! Items from the retired cbound.h include file

INTEGER, PARAMETER :: max_rimwidth = 30    ! Maximum rim width for LBCs
REAL :: rimweightsa(max_rimwidth)          ! Rim weights

INTEGER :: current_lbc_step         ! Timestep at which LBCs were last updated
INTEGER :: albc_num                 ! Number of atmos boundary file currently
                                    !  in use
INTEGER :: albc2_starttime_steps    ! VT of first block of data in 2nd atmos
                                    !  boundary file, in steps from start of run

INTEGER :: albc_swapstep            ! Step on which to swap to second atmos
                                    ! boundary_file

NAMELIST /bouncnst/ rimweightsa

END MODULE lbc_read_data_mod
