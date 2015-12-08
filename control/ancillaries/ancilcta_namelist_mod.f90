! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Ancil-related logicals read in by namelist
!
MODULE ancilcta_namelist_mod

! Description:
!   Module containing the ANCILCTA namelist and associated declarations.
!
! Method:
!  Contents of the ANCILCTA namelist are read in and stored here.
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Ancillaries
IMPLICIT NONE

LOGICAL :: l_sstanom = .FALSE. ! Indicator if SST anom to be formed
                               !  (RECON=T) or used (-DEF,RECON)
LOGICAL :: lamipii   = .FALSE. ! True if special AMIP II updating


NAMELIST/ANCILCTA/                                                 &
  l_sstanom, lamipii 

END MODULE ancilcta_namelist_mod
