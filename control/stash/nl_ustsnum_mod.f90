! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! User STASHmaster information

MODULE nl_ustsnum_mod

IMPLICIT NONE
!
! Description:
!   Contains the namelist information for reading userSTASHmaster data.
!   For the Unified Model this is held in a single PRESM_A file.
!   For def(UTILIO) up to 20 named userSTASHmaster files may be used.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!

      INTEGER      n_uSTASH        ! Number of user STASHmaster files




      INTEGER      nrecs_uSTASH    ! Total no. of user stash records
      CHARACTER(LEN=8)  ustsfils(20)    ! Names of user ppxref file (PRESM_A)

      NAMELIST/USTSNUM /n_uSTASH,nrecs_uSTASH,ustsfils

END MODULE nl_ustsnum_mod
