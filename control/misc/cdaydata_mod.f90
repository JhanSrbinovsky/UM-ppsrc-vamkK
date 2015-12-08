! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc
MODULE cdaydata_mod

IMPLICIT NONE
! CDATDATA start
!
! Constants needed by routines to calculate day/month/year from
! incremental day number since calendar zero point, and vice-versa,
! when using Gregorian calendar

      INTEGER,PARAMETER:: DAYS_PER_4C = 146097
      INTEGER,PARAMETER:: DAYS_PER_C  = 36524
      INTEGER,PARAMETER:: DAYS_PER_4Y = 1461
      INTEGER,PARAMETER:: DAYS_PER_Y  = 365

      INTEGER,PARAMETER:: DAYS_IN_MONTH(12) =                           &
     &  (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

      INTEGER,PARAMETER:: DAYS_TO_MONTH(12) =                           &
     &  (/0, 31, 59, 90,120,151,181,212,243,273,304,334/)

! CDAYDATA end

END MODULE cdaydata_mod
