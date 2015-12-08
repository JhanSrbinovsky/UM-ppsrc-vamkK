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
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine write_one_field (                                      &
!----------------------------------------------------------------------
! comdeck: AFIELDS
! Purpose: argument list for all dynamically allocated arrays used by
!          FOAM_Flux_Process.
!          This deck is linked to CFIELDS.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AFLDDIMS
! Purpose: argument list for dimensions of fields.
!          This deck is linked to FLDDIMS.
!----------------------------------------------------------------------
     & ncols, nrowst, nrowsu, nrowsv, nrowsuv, ncolsO, nrowstO, nrowsuO,&
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: ALSMS
! Purpose: argument list for land / sea masks.
!          This deck is linked to CLSMS.
!----------------------------------------------------------------------
     &  lsmt, lsmu, lsmtO, lsmuO,                                       &
!----------------------------------------------------------------------
! ACOORDS argument list for latitude and longitude grid coordinates.
! This file is linked to CCOORDS.
     & lambda_t, phi_t, lambda_u, phi_u,                                &
     & lambda_tO, phi_tO, lambda_uO, phi_uO,                            &
! ACOORDS end
!----------------------------------------------------------------------
! comdeck: AINTERP
! Purpose: argument list for interpolation coefficients for
!          interpolation from atmosphere to ocean grid.
!          This deck is linked to CINTERP.
!----------------------------------------------------------------------
     & index_bl_t, index_br_t, index_bl_u, index_br_u,                  &
     & weight_tr_t, weight_bl_t, weight_br_t, weight_tl_t,              &
     & weight_tr_u, weight_bl_u, weight_br_u, weight_tl_u,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AFILLIN
! Purpose: argument list for control of filling in using spiral fill
!          This deck is linked to CFILLIN.
!----------------------------------------------------------------------
     & n_pts_unres_t, index_unres_t, n_pts_unres_u, index_unres_u,      &
     & n_calls_spiral_t,n_pts_spiral_t,n_calls_spiral_u,n_pts_spiral_u, &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AROTGRD
! Purpose: argument list for variables associted with the
!          definition of the rotated grid.
!          This deck is linked to CROTGRID
!----------------------------------------------------------------------
     & rotg, rotgO, pole_lat, pole_lon, poleO_lat, poleO_lon,           &
     & coef_angle1, coef_angle2, coef_angle3, coef_angle4,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
     &       StCode, FFCode, PPCode, IVTOffHr,                          &
     &       IGridtype, nrows,                                          &
     &       Int_Head, Real_Head, IOutUnit, ldebug,                     &
     &       field_atm, icode )

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn, cerr,&
                           cstd, csub
      USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,         &
                              len_realhd, itgrid, iugrid,                &
                              max_num_fc_times, max_num_clim_times,      &
                              max_num_in_flux, len2_lookuppreferred,     &
                              len2_lookupprevious, len2_lookupclimate
      implicit none

! declaration of parameters used in argument list

! declaration of argument list

! array dimensions, lsms, interpolation coeffs etc. : all intent IN
!----------------------------------------------------------------------
! comdeck: CFIELDS
! Purpose: local declaration of all dynamically allocated arrays used by
!          FOAM_Flux_Process.
!          This deck is linked to AFIELDS.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: CFLDDIMS
! Purpose: declares dimensions of fields.
!          This deck is linked to AFLDDIMS.
!----------------------------------------------------------------------
! declarations:

!atmosphere grid
      integer ncols     !  number of columns
      integer nrowst    !  number of rows (tracer grid)
      integer nrowsu    !  number of rows (u grid)
      integer nrowsv    !  number of rows (v grid)
      integer nrowsuv   !  number of rows on grid with coincident velys
                        !  (at B grid locations both for B & C grids)
!ocean grid
      integer ncolsO    !  number of columns
      integer nrowstO   !  number of rows (tracer grid)
      integer nrowsuO   !  number of rows (velocity grid)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: CLSMS
! Purpose: local declarations of land / sea masks.
!          This deck is linked to ALSMS.
!----------------------------------------------------------------------
! declarations:

!land/sea masks for atmosphere grid
      integer lsmt(ncols,nrowst)        ! t grid
      integer lsmu(ncols,nrowsuv)      ! uv grid (u v at coincident pts)
!land/sea masks for ocean grid
      integer lsmtO(ncolsO,nrowstO)     ! t grid
      integer lsmuO(ncolsO,nrowsuO)     ! u grid
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: CCOORDS
! Purpose: local declarations of latitude and longitude grid
!          coordinates.
!          This deck is linked to ACOORDS.
!----------------------------------------------------------------------
! declarations:
!coordinates of grids
      real lambda_t(ncols)   ! coords of longitudes: atmosphere tracer
      real phi_t(nrowst)     ! coords of latitudes : atmosphere tracer

      real lambda_u(ncols)  ! coords of longitudes: { atmos vely, u & v
      real phi_u(nrowsuv)   ! coords of latitudes:  { at coincident pts

      real lambda_tO(ncolsO) ! coords of longitudes: ocean tracer
      real phi_tO(nrowstO)   ! coords of latitudes : ocean tracer

      real lambda_uO(ncolsO) ! coords of longitudes: ocean velocity
      real phi_uO(nrowsuO)   ! coords of latitudes:  ocean velocity
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: CINTERP
! Purpose: declares interpolation coefficients for interpolation from
!          atmosphere to ocean grid.
!          This deck is linked to AINTERP.
!----------------------------------------------------------------------
! declarations:

!indices of  ## corners of source gridbox
!for tracer grid interpolation
      integer index_bl_t(ncolsO*nrowstO)  ! bottom lefthand  tracer
      integer index_br_t(ncolsO*nrowstO)  ! bottom righthand tracer

!Weights applied to value at ## corners of source gridbox
!for tracer grid interpolation
      real weight_tr_t(ncolsO*nrowstO)  ! top right    tracer
      real weight_bl_t(ncolsO*nrowstO)  ! bottom left  tracer
      real weight_br_t(ncolsO*nrowstO)  ! bottom right tracer
      real weight_tl_t(ncolsO*nrowstO)  ! top left     tracer

!indices of  ## corners of source gridbox
!for velocity grid interpolation
      integer index_bl_u(ncolsO*nrowsuO)  ! bottom lefthand  velocity
      integer index_br_u(ncolsO*nrowsuO)  ! bottom righthand velocity

!Weight applied to value at ## corner of source gridbox
!for velocity grid interpolation
      real weight_tr_u(ncolsO*nrowsuO)  ! top right    velocity
      real weight_bl_u(ncolsO*nrowsuO)  ! bottom left  velocity
      real weight_br_u(ncolsO*nrowsuO)  ! bottom right velocity
      real weight_tl_u(ncolsO*nrowsuO)  ! top left     velocity
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: CFILLIN
! Purpose: declares indices for points on ocean grid which need
!          filling in and values controlling the use of spiral fill.
!          This deck is linked to AFILLIN.
!----------------------------------------------------------------------
! declarations:

! parameters
      integer max_no_searches
      parameter ( max_no_searches = 10)

! indices to points on ocean grid which need filling in (i.e.
! seapoints on ocean grid which are not fully surrounded by
! seapoints on the atmosphere grid)

! for tracer grid
      integer n_pts_unres_t                  ! number of unresolved pts
      integer index_unres_t(ncolsO*nrowstO)  ! indices for each pt

! for velocity grid
      integer n_pts_unres_u                  ! number of unresolved pts
      integer index_unres_u(ncolsO*nrowsuO)  ! indices for each pt

!control of first calls to spiral search: tracer grid
      integer n_calls_spiral_t   ! number of times to call spiral search
      integer n_pts_spiral_t(max_no_searches) ! # of pts (nsearch)
                                                ! for each call

!control of calls of spiral search: velocity grid
      integer n_calls_spiral_u   ! number of times to call spiral search
      integer n_pts_spiral_u(max_no_searches)   ! # of pts (nsearch)
                                                ! for each call
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: CROTGRD
! Purpose: declares variables associted with the
!          definition of the rotated grid.
!          This deck is linked to AROTGRID
!----------------------------------------------------------------------
! declarations:

!logicals to define type of grid, T => rotated grid
      logical rotg    ! T => atmosphere on a rotated grid
      logical rotgO   ! T => ocean on a rotated grid

!positions of the pole on the rotated grid
      real pole_lat   ! latitude of pole for atmosphere
      real pole_lon   ! longitude of pole for atmosphere
      real poleO_lat  ! latitude of pole for ocean
      real poleO_lon  ! longitude of pole for ocean

!coefficients for converting the wind vectors to the correct direction
! for use with wind vector components defined at coincident points

      real coef_angle1(ncols,nrowsuv) ! atmosphere to
                                      ! standard lat-long grid
      real coef_angle2(ncols,nrowsuv) ! atmosphere to
                                      ! standard lat-long grid
      real coef_angle3(ncols,nrowsuv) ! standard lat-long to ocean grid
      real coef_angle4(ncols,nrowsuv) ! standard lat-long to ocean grid

!----------------------------------------------------------------------
!----------------------------------------------------------------------

! field codes to insert in integer header that is output
      integer StCode   ! IN stash code
      integer FFCode   ! IN Met O 8 field code
      integer PPCode   ! IN PP package code
      integer IVTOffHr ! IN offset of validity time from reference

! other input
      integer IGridtype  ! IN  grid type (0 = tracer, 1 = velocity)
      integer nrows      ! IN  number of rows in input field

      integer Int_Head(Len_IntHd)   ! IN integer part of lookup table
      real Real_Head(Len_RealHd)    ! IN real part of lookup table
      integer IOutUnit   ! IN  output unit
      logical ldebug     ! IN  T => output debugging info

      real field_atm( ncols, nrows ) ! IN  field on NWP grid
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected


! declaration of parameters
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

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
!----------------------------------------------------------------------
! comdeck: CLOOKUPS
! Purpose: declares and stores lookup tables for flux input files.
!          Also stores LCyclic, ISeaPt and ILandPt.
!          Some arrays are dimensioned using parameters in PLOOKUPS
!          so *CALL CLOOKUPS must always be preceded by *CALL PLOOKUPS
!----------------------------------------------------------------------
! parameters
      integer ISeaPt    ! value of land / sea mask at a sea point
                        ! ISeaPt = 0 because land / sea mask is
                        ! zero at sea points and 1 at land points
      integer ILandPt   ! ILandPt = 1
      parameter ( ISeaPt = 0, ILandPt = 1 )

! common block:
      common / Lookups /                                                &
     &   Len2_ActualPreferred, Len2_ActualPrevious, Len2_ActualClimate, &
     &   LCyclic, LCyclicO,                                             &
     &   FixHdPreferred, FixHdPrevious, FixHdClimate,                   &
     &   FixHdlsmt,FixHdlsmtO, FixHdlsmuO,                              &
     &   LookupPreferred, LookupPrevious, LookupClimate,                &
     &   LookFldNoPreferred, LookFldNoPrevious, LookFldNoClimate,       &
     &   Lookuplsmt, LookuplsmtO, LookuplsmuO


! actual numbers of fields in files
      integer Len2_ActualPreferred ! for NWP file
      integer Len2_ActualPrevious  ! for NWP file
      integer Len2_ActualClimate   ! for climate file

! additional information on fields read in
      logical LCyclic   ! T => atmosphere grid is cyclic
      logical LCyclicO  ! T => ocean grid is cyclic

! fixed headers
      integer FixHdPreferred(Len_FixHd)  ! for preferred NWP file
      integer FixHdPrevious(Len_FixHd)   ! for previous NWP file
      integer FixHdClimate(Len_FixHd)    ! for climate file
      integer FixHdlsmt(Len_FixHd)       ! for atmos. tracer lsm
      integer FixHdlsmtO(Len_FixHd)      ! for ocean tracer lsm
      integer FixHdlsmuO(Len_FixHd)      ! for ocean velocity lsm

! lookup tables for NWP preferred and previous files and climate files
      integer LookupPreferred(Len1_Lookup, Len2_LookupPreferred)
      integer LookupPrevious(Len1_Lookup, Len2_LookupPrevious)
      integer LookupClimate(Len1_Lookup, Len2_LookupClimate)

! additional lookup entry indicating the field (and lookup table)
! to use to access data with this validity time and stash code.
! This allows the user to specify copies of data with amended dates to
! be used. (See subroutine add_lookups)
      integer LookFldNoPreferred(Len2_LookupPreferred)
      integer LookFldNoPrevious(Len2_LookupPrevious)
      integer LookFldNoClimate(Len2_LookupClimate)

! lookup tables for land / sea masks for 4 grids:
      integer Lookuplsmt(Len1_Lookup)   ! atmosphere tracer
!  Lookuplsmu(Len1_Lookup)   is not used from vn 5.3
      integer LookuplsmtO(Len1_Lookup)  ! FOAM ocean tracer
      integer LookuplsmuO(Len1_Lookup)  ! FOAM ocean velocity
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck CVALOFF
! Purpose: stores user inputs relating to validity times and units to
!          use for output of FOAM fluxes
!----------------------------------------------------------------------

      ! maximum number of validity times
      INTEGER,PARAMETER:: MaxTimes = 50

      ! period in hours covered by each of output flux fields
                              ! (6 hours for operational FOAM)
      INTEGER :: ValidityPeriod

      ! number of (validity) times to process
      INTEGER :: NoValidTimes

      ! offset from reference time of validity time of end of flux
      ! period (in hours)
      INTEGER :: IValidOffHr(MaxTimes)

      ! offset from main output unit. Used to output fields for last
                                    ! validity time to separate files
      INTEGER :: IOutUnitOff(MaxTimes)

! control for inserting additional "copies" of lookup tables / fields

      ! For preferred file
      INTEGER:: NoAddTimesPreferred ! number of lookup tables to insert
      INTEGER:: ISrchOffHrPreferred(MaxTimes) ! offset hours to look for
      INTEGER:: INewOffHrPreferred(MaxTimes) ! new offset hours

      ! For previous file
      INTEGER:: NoAddTimesPrevious ! number of lookup tables to insert
      INTEGER:: ISrchOffHrPrevious(MaxTimes) ! offset hour to look for
      INTEGER:: INewOffHrPrevious(MaxTimes) ! new offset hour

      ! For climate file
      INTEGER:: NoAddTimesClimate ! number of lookup tables to insert
      INTEGER:: ISrchOffHrClimate(MaxTimes) ! offset hour to look for
      INTEGER:: INewOffHrClimate(MaxTimes) ! new offset hour

      ! Value to use at land points in final output file (the default
      ! value is rmdi; for testing it is sometimes useful  to set this
      ! value to zero).
      REAL :: output_land_value   ! value at land points

      COMMON / ValOff / ValidityPeriod,                                 &
     &  NoValidTimes, IValidOffHr, IOutUnitOff,                         &
     &  NoAddTimesPreferred, ISrchOffHrPreferred, INewOffHrPreferred,   &
     &  NoAddTimesPrevious, ISrchOffHrPrevious, INewOffHrPrevious,      &
     &  NoAddTimesClimate, ISrchOffHrClimate, INewOffHrClimate,         &
     &  output_land_value

! CVALOFF end

! declaration of local arrays
      real field_ocean( ncolsO, nrowstO ) ! used for both t and u cases
      integer index_unres ( ncolsO * nrowstO ) ! indices to unresolved
                                               ! points on ocean grid
!----------------------------------------------------------------------

! declaration of local scalars

      integer ipts          ! loop index over unresolved points
      integer isearch       ! loop index over calls to spiral_s
      integer nsearch       ! # of pts in search "radius"
      integer n_pts_unres   ! local counter of # of unresolved points
      integer ncolsOut      ! # of columns in output field
      integer nrowsOut      ! # of rows in output field

      external lsm_set,h_int_lsm,spiral_s,write_pp
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'write_one_field'  ! subroutine name for error messages

! 0.1 check that nrows and IGridtype are consistent
      if (  IGridtype  ==  0 ) then

        if ( nrows  /=  nrowst ) then
          icode = 44
          write(UnErr,*)CErr,CSub,                                      &
     &       ' 0.1.1 nrows and IGridtype inconsistent: ',               &
     &       ' nrows, nrowst, IGridtype =', nrows, nrowst, IGridtype
          go to 9999
        end if

      else if ( IGridtype  ==  1 ) then

        if ( nrows  /=  nrowsuv ) then
          icode = 45
          write(UnErr,*)CErr,CSub,                                      &
     &       ' 0.1.2 nrows and IGridtype inconsistent: ',               &
     &       ' nrows, nrowsuv, IGridtype =', nrows, nrowsuv, IGridtype
          go to 9999
        end if

      else

        icode = 46
        write(UnErr,*)CErr,CSub,                                        &
     &       ' 0.1.3 not coded for IGridtype =', IGridtype
          go to 9999

      end if ! IGridtype

! 1. Set land points to missing data (use atmosphere grids)

! for tracer grid
      if ( IGridtype  ==  0 ) then
! DEPENDS ON: lsm_set
        call lsm_set( ncols, nrows, lsmt, ILandPt,                      &
     &       rmdi, ldebug, field_atm )

      else if ( IGridtype  ==  1 ) then
! DEPENDS ON: lsm_set
        call lsm_set( ncols, nrows, lsmu, ILandPt,                      &
     &       rmdi, ldebug, field_atm )

      end if

! 2. Interpolate to ocean grid

      if ( IGridtype  ==  0) then

        ncolsOut = ncolsO
        nrowsOut = nrowstO
! DEPENDS ON: h_int_lsm
        call h_int_lsm(nrowst,ncols,ncolsOut*nrowsOut, rmdi,            &
     &     index_bl_t,index_br_t, field_atm,                            &
     &     weight_bl_t,weight_br_t,weight_tl_t,weight_tr_t,             &
     &     lsmtO,                                                       &
     &     field_ocean)

      else if ( IGridtype  ==  1) then

        ncolsOut = ncolsO
        nrowsOut = nrowsuO
! DEPENDS ON: h_int_lsm
        call h_int_lsm(nrowsu,ncols,ncolsOut*nrowsOut, rmdi,            &
     &     index_bl_u,index_br_u, field_atm,                            &
     &     weight_bl_u,weight_br_u,weight_tl_u,weight_tr_u,             &
     &     lsmuO,                                                       &
     &     field_ocean)

      end if

! 3. fill in coastal values

! 3.1 for a tracer grid
      if ( IGridtype  ==  0) then

! 3.1.1 copy unresolved points into a local array (which is
!       updated by each call to spiral_s)

        n_pts_unres = n_pts_unres_t
        do ipts = 1, n_pts_unres
          index_unres(ipts) = index_unres_t(ipts)
        end do

! 3.1.2 do spiral searches

        do isearch = 1, n_calls_spiral_t

          nsearch = n_pts_spiral_t(isearch)

! DEPENDS ON: spiral_s
          call spiral_s(lsmtO,index_unres,n_pts_unres,                  &
     &      nrowsOut,ncolsOut,field_ocean,nsearch,ISeaPt,LCyclicO)

        end do ! isearch

! 3.2 for a velocity grid
      else if ( IGridtype  ==  1) then

! 3.2.1 copy unresolved points into a local array (which is
!       updated by each call to spiral_s)

        n_pts_unres = n_pts_unres_u
        do ipts = 1, n_pts_unres
          index_unres(ipts) = index_unres_u(ipts)
        end do

! 3.2.2 do spiral searches

        do isearch = 1, n_calls_spiral_u

          nsearch = n_pts_spiral_u(isearch)

! DEPENDS ON: spiral_s
          call spiral_s(lsmuO,index_unres,n_pts_unres,                  &
     &      nrowsOut,ncolsOut,field_ocean,nsearch,ISeaPt,LCyclicO)

        end do ! isearch


      end if  ! IGridtype

! 4. Reset missing data values at land points if user has
!    chosen to do so

      if ( output_land_value   /=  rmdi ) then
        if ( IGridtype  ==  0) then
! DEPENDS ON: lsm_set
          call lsm_set( ncolsOut, nrowsOut, lsmtO, ILandPt,             &
     &                  output_land_value, ldebug, field_ocean )
        else if ( IGridtype  ==  1 ) then
! DEPENDS ON: lsm_set
          call lsm_set( ncolsOut, nrowsOut, lsmuO, ILandPt,             &
     &                  output_land_value, ldebug, field_ocean )
        end if
      end if

! 5. Amend grid information in lookup table
      if ( IGridtype  ==  0) then
! DEPENDS ON: amend_lookup
        call amend_lookup (  LookuplsmtO, Int_Head, Real_Head,          &
     &                       output_land_value,                         &
     &                       StCode, FFCode, PPCode, IVTOffHr )

      else if ( IGridtype  ==  1 ) then
! DEPENDS ON: amend_lookup
        call amend_lookup (  LookuplsmuO, Int_Head, Real_Head,          &
     &                       output_land_value,                         &
     &                       StCode, FFCode, PPCode, IVTOffHr )

      end if

! 6. write out filled pp field on ocean grid
! DEPENDS ON: write_pp
      call write_pp(IOutUnit, Int_Head, Real_Head,                      &
     &              ncolsOut, nrowsOut, field_ocean, icode)
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 5. error writing out a pp header and field  '
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE write_one_field
!----------------------------------------------------------------------
!----------------------------------------------------------------------
