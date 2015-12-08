! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: write_one_field,write_pp
!
! Purpose: Flux processing routine.
!          write_one_field:
!          This routine writes out one pp-field to the required output
!          pp file. The pp-header is atered for the correct stashcode
!          and the field written to the file using the routine write_pp.
!          The routine also sets land points to missing data and checks
!          for consistency in the field dimensions.
!          write_pp: Writes out a pp header then a pp field
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine write_pp ( IOutUnit, Int_Head, Real_Head,              &
     &                      ncolsOut, nrowsOut, field_out, icode)



      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn, cerr,&
                           cstd, csub
      USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,         &
                              len_realhd, itgrid, iugrid,                &
                              max_num_fc_times, max_num_clim_times,      &
                              max_num_in_flux, len2_lookuppreferred,     &
                              len2_lookupprevious, len2_lookupclimate
      implicit none

! declaration of parameters

! declaration of argument list
      integer IOutUnit  ! IN output unit number
      integer Int_Head(Len_IntHd)  ! integer part of lookup table
      real Real_Head(Len_RealHd)   ! real part of lookup table
      integer ncolsOut             ! # of columns in output field
      integer nrowsOut             ! # of rows in output field
      real field_out( ncolsOut, nrowsOut ) ! field output
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of globals used
! 5.3      22/10/01     Update input file units. A. Hines
! 6.1      30/07/04     Update input file units for UM6.0. A. Hines
! CUNITNOS defines parameters defining unit numbers and logicals for
! units of fields which are open

      ! for input control files
      INTEGER, PARAMETER :: UnitDbg = 7
      INTEGER, PARAMETER :: UnitHK  = 8
      INTEGER, PARAMETER :: UnitVT  = 9
      INTEGER, PARAMETER :: UnitSlt = 94

      ! for flux fields
      INTEGER, PARAMETER :: UnitPreferred = 10
      INTEGER, PARAMETER :: UnitPrevious  = 11
      INTEGER, PARAMETER :: UnitClimate   = 12

      ! for land / sea masks
      INTEGER, PARAMETER :: UnitNWPlsmt  = 13 ! NWP tracer grid
      INTEGER, PARAMETER :: UnitFOAMlsmt = 15 ! FOAM tracer grid
      INTEGER, PARAMETER :: UnitFOAMlsmu = 16 ! FOAM velocity grid

      ! for output of ancillary file validity times namelist
      ! output ancillary file validity times
      INTEGER, PARAMETER :: UnitVTOut = 18

      ! for main set of output pp files
      ! wind stress & mixing energy
      INTEGER,PARAMETER:: UnitWindsOut     = 20

      INTEGER,PARAMETER:: UnitHeatOut      = 21 ! heat fluxes
      INTEGER,PARAMETER:: UnitMoistureOut  = 25 ! moisture fluxes
      INTEGER,PARAMETER:: UnitSeaIceOut    = 23 ! sea-ice field
      INTEGER,PARAMETER:: UnitReferencesOut= 24 ! reference fields
      INTEGER,PARAMETER:: UnitPressureOut  = 26 ! pressure fields
      INTEGER,PARAMETER:: UnitWindspdOut   = 27 ! wind speed

! for output units
      INTEGER,PARAMETER:: IUnOutLow        = 20
      INTEGER,PARAMETER:: IUnOutHi         = 44

      COMMON / OutUnits / LUnOutOpen, LPreferred, LPrevious, LClimate

      !  logical T => output ! unit is open
      LOGICAL :: LUnOutOpen(IUnOutLow:IUnOutHi)

      LOGICAL :: LPreferred ! T => preferred NWP file is available
      LOGICAL :: LPrevious  ! T => previous NWP file is available
      LOGICAL :: LClimate   ! T => climate file is available
! CUNITNOS end

! no local arrays

! declaration of local scalars
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'write_pp'  ! subroutine name for error messages

! 1. Write out header
      write (IOutUnit, IOStat = icode) Int_Head, Real_Head
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1. error writing out lookup table  '
        icode = 47
        go to 9999
      end if

! 2. Write out data
      write (IOutUnit, IOStat = icode) field_out
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2. error writing out data field  '
        icode = 48
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE write_pp
!----------------------------------------------------------------------
