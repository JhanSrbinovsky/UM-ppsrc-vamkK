! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Initailise a stash code to new grib code table
!
! Subroutine Interface:
! ======================================================================
!+ Initailise a stash code to new user grib code table
!
! Subroutine Interface:
      SUBROUTINE GRIB_TABLE_INIT2

      IMPLICIT NONE
!
! Description:
!   Initialises array GRIB_TABLE used to map stash codes to
!   a user defined  grib table 2 set of field codes
!
! Method:
!   Uses special set of mapping for AMIP highres experiment.
! NOTE : This subroutine is present  to allow a user to
!        edit it for their own set of mappings.
!        Codes must be between 129 and 255 to conform to user
!        extensions to grib edition 1 table 2.
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
      GRIB_TABLE(9,206)=129   ! cloud liquid water vertical mean
      GRIB_TABLE(9,207)=130   ! cloud ice content vertical mean
      GRIB_TABLE(1,201)=131   ! net surface downward SW radiation
      GRIB_TABLE(1,207)=132   ! Incoming SW TOA
      GRIB_TABLE(1,208)=133   ! outgoing SW TOA
      GRIB_TABLE(1,209)=134   ! outgoing SW TOA clear sky
      GRIB_TABLE(1,210)=135   ! downward surface SW clear sky
      GRIB_TABLE(1,235)=136   ! total downward surface SW
      GRIB_TABLE(2,201)=137   ! net surface downward LW radiation
      GRIB_TABLE(2,204)=138   ! Total cloud amount
      GRIB_TABLE(2,205)=139   ! OLR TOA
      GRIB_TABLE(2,206)=140   ! OLR TOA  clear sky
      GRIB_TABLE(2,207)=141   ! downward surface LW
      GRIB_TABLE(2,208)=142   ! downward surface LW clear sky
      GRIB_TABLE(3,217)=143   ! surface sensible heat flux
      GRIB_TABLE(3,219)=144   ! u - comp of surface wind stress
      GRIB_TABLE(3,220)=145   ! v - comp of surface wind stress
      GRIB_TABLE(3,223)=146   ! surface  evaportaion rate
      GRIB_TABLE(3,225)=147   ! u - comp of wind 10m
      GRIB_TABLE(3,226)=148   ! v - comp of wind 10m
      GRIB_TABLE(3,234)=149   ! latent heat flux
      GRIB_TABLE(3,236)=150   ! 1.5m temperature
      GRIB_TABLE(3,237)=151   ! 1.5m specific humidity
      GRIB_TABLE(3,238)=152   ! deep soil temperature level 1
      GRIB_TABLE(3,245)=153   ! 1.5m relative humidity
      GRIB_TABLE(4,203)=154   ! surface rainfall rate (L_S)
      GRIB_TABLE(4,204)=155   ! surface snowfall rate (L_S)
      GRIB_TABLE(5,205)=156   ! surface snowfall rate (CONV)
      GRIB_TABLE(5,206)=157   ! surface rainfall rate (CONV)
      GRIB_TABLE(8,23)=158   !snow mass
      GRIB_TABLE(8,204)=159   ! surface runoff
      GRIB_TABLE(8,205)=160   ! sub-surface runoff
      GRIB_TABLE(8,208)=161   ! soil moisture
      GRIB_TABLE(8,209)=162   ! canopy water content
      GRIB_TABLE(6,201)=163   ! u -comp of GWD stress
      GRIB_TABLE(6,202)=164   ! v -comp of GWD stress
      GRIB_TABLE(16,222)=165   ! mslp
      GRIB_TABLE(0,1)=166      ! surface pressure
      GRIB_TABLE(0,10)=167      ! specific humidity vertical mean
      GRIB_TABLE(0,24)=168      ! surface temperature
! --------------------------------------------------------------------
      RETURN
      END SUBROUTINE GRIB_TABLE_INIT2
