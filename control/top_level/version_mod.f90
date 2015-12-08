! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! Purpose: STASH parameter definitions

MODULE version_mod

IMPLICIT NONE

! Max. no. of STASH sections per internal model
INTEGER,PARAMETER :: nsectp=99
! Max. no. of STASH items per section
INTEGER,PARAMETER :: nitemp=999
! Max. no. of STASH list records (prognostic + diagnostic)
INTEGER,PARAMETER :: nrecdp=1500
! Max. no. of output times tables in STASHC
INTEGER,PARAMETER :: ntimep=100
! Max. no. of time profiles in STASHC
INTEGER,PARAMETER :: nproftp=100
! Max. no. of domain profiles/levels lists in STASHC (used for both)
INTEGER,PARAMETER :: nprofdp=100
! Max. total no. of time series in STASHC
INTEGER,PARAMETER :: NTimSerP=1500
! Max. no. time series per domain profile
INTEGER,PARAMETER :: tsdp=250
! Max. no. of useage profiles in STASHC
INTEGER,PARAMETER :: nprofup=40
! Max. no. of levels in a levels list
INTEGER,PARAMETER :: nlevp=50
! Max. no. of pseudo levels in a  pseudo levels list
INTEGER,PARAMETER :: npslevp=100
! Max. no. of pseudo levels lists in STASHC
INTEGER,PARAMETER :: npslistp=100
! Max. no. non-blank records in PPXREF file
INTEGER,PARAMETER :: ndiagp=4000
INTEGER,PARAMETER :: ndiagpm=nrecdp
INTEGER,PARAMETER :: nelemp=33
INTEGER,PARAMETER :: nlevp_S=nlevp*6+1
INTEGER,PARAMETER :: nlevlstsp=nprofdp

! OUTFILE_E must be consistent with NUNITS in CHSUNITS.h
! Ranges of output file numbers
INTEGER,PARAMETER :: outfile_s=20
INTEGER,PARAMETER :: outfile_e=179

END MODULE version_mod
