! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Initailise a stash code to new grib code table
!
! Subroutine Interface:
      SUBROUTINE GRIB_TABLE_INIT1

      IMPLICIT NONE
!
! Description:
!   Initialises array GRIB_TABLE used to map stash codes to
!   standard table 2 grib grib edition 1 field codes
!
! Method:
!   See WMO documentation on grib for list of codes.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.0     31/03/95  Original code. R.A.Stratton
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:
!
! -------------------------------------------------------------------
! Declarations:
!
! CGRIBTAB start
!
! Description:
!  Holds a lookup table for converting between Unified Model
! stash codes and some other grib table of codes
!
! Declarations:
      INTEGER, PARAMETER:: MAX_SECT_GRBTAB=16
      INTEGER, PARAMETER:: MAX_ITEM_GRBTAB=300

      INTEGER :: GRIB_TABLE(0:MAX_SECT_GRBTAB,MAX_ITEM_GRBTAB)

      COMMON /GRIBTAB/ GRIB_TABLE

! CGRIBTAB end

! Subroutine arguments - none

! Local varaibles:
      INTEGER                                                           &
     &      I,J        ! loop counters

! Function & Subroutine calls: none
! --------------------------------------------------------------------
! Initilise table with -99 - indicate no matching code

      DO I=0,MAX_SECT_GRBTAB,1
        DO J=1,MAX_ITEM_GRBTAB
          GRIB_TABLE(I,J)=-99
        ENDDO
      ENDDO
! --------------------------------------------------------------------
!  Change entries for field with matching grib codes
      GRIB_TABLE(16,222)=2     ! mslp
      GRIB_TABLE(15,201)=3     ! pressure tendency
      GRIB_TABLE(10,206)=6     ! geopotential
      GRIB_TABLE(10,206)=6     ! geopotential
      GRIB_TABLE(16,201)=7     ! geopotential height
      GRIB_TABLE(16,202)=7     ! geopotential height
      GRIB_TABLE(1,  4)=11    ! temperature
      GRIB_TABLE(2,  4)=11    ! temperature
      GRIB_TABLE(3,  4)=11    ! temperature
      GRIB_TABLE(4,  4)=11    ! temperature
      GRIB_TABLE(5,209)=11    ! temperature
      GRIB_TABLE(9,209)=11    ! temperature
      GRIB_TABLE(15,216)=11    ! temperature
      GRIB_TABLE(16,203)=11    ! temperature
      GRIB_TABLE(0,4)=13       ! potential temperature
      GRIB_TABLE(5,4)=13       ! potential temperature
      GRIB_TABLE(10,4)=13      ! potential temperature
      GRIB_TABLE(0,2)=33      ! u_component of wind speed
      GRIB_TABLE(3,2)=33      ! u_component of wind speed
      GRIB_TABLE(6,2)=33      ! u_component of wind speed
      GRIB_TABLE(7,2)=33      ! u_component of wind speed
      GRIB_TABLE(10,2)=33      ! u_component of wind speed
      GRIB_TABLE(12,2)=33      ! u_component of wind speed
      GRIB_TABLE(13,2)=33      ! u_component of wind speed
      GRIB_TABLE(15,201)=33      ! u_component of wind speed
      GRIB_TABLE(0,3)=34      ! v_component of wind speed
      GRIB_TABLE(3,3)=34      ! v_component of wind speed
      GRIB_TABLE(6,3)=34      ! v_component of wind speed
      GRIB_TABLE(7,3)=34      ! v_component of wind speed
      GRIB_TABLE(10,3)=34      ! v_component of wind speed
      GRIB_TABLE(12,3)=34      ! v_component of wind speed
      GRIB_TABLE(13,3)=34      ! v_component of wind speed
      GRIB_TABLE(15,202)=34      ! v_component of wind speed
      GRIB_TABLE(12,201)=39      ! vertical velocity
      GRIB_TABLE(12,202)=39      ! vertical velocity
      GRIB_TABLE(15,222)=39      ! vertical velocity
      GRIB_TABLE(0,10)=51      ! specific humidity
      GRIB_TABLE(3,10)=51      ! specific humidity
      GRIB_TABLE(4,10)=51      ! specific humidity
      GRIB_TABLE(5,10)=51      ! specific humidity
      GRIB_TABLE(9,10)=51      ! specific humidity
      GRIB_TABLE(15,226)=51      ! specific humidity
      GRIB_TABLE(16,204)=52      ! relative humidity
      GRIB_TABLE(3,223)=57      ! evaporation (units not correct?)
      GRIB_TABLE(5,216)=59      ! precipitation rate
      GRIB_TABLE(2,204)=71      ! total cloud cover (should be %)
      GRIB_TABLE(0,13)=72      ! total convective cloud
      GRIB_TABLE(5,13)=72      ! total convective cloud
      GRIB_TABLE(9,203)=73      ! low cloud
      GRIB_TABLE(9,204)=74      ! medium cloud
      GRIB_TABLE(9,205)=75      ! high cloud
      GRIB_TABLE(3,238)=85      ! soil temperature
      GRIB_TABLE(8,208)=86      ! soil moisture
      GRIB_TABLE(2,201)=112     ! net long-wave radiation surface
      GRIB_TABLE(2,205)=114     ! net long-wave radiation toa
      GRIB_TABLE(3,234)=121     ! latent heat flux
      GRIB_TABLE(3,217)=122     ! sensible heat flux
! --------------------------------------------------------------------
      RETURN
      END SUBROUTINE GRIB_TABLE_INIT1
! ======================================================================
!+ Initailise a stash code to new user grib code table
!
! Subroutine Interface:
