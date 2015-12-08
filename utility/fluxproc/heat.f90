! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: heat
!
! Purpose: Flux processing routine.
!          To produce a pp file containing:
!          Net Penetraing Solar Radiation (SOL)
!          Net non Penetraing Heat        (HTN)
!          for the times required.
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine heat(                                                  &
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
     &                 icode )

      USE cfdcodes_mod, ONLY: stcwindstressu, stcwindstressv,           &
                              stcwindmixeng, stcsw, stcsw1, stclongwave,&
                              stcsensibleheat, stcsublim, stctopmelt,   &
                              stcbotmelt, stcevaporation, stcdrain,     &
                              stcconvrain, stcdsnow, stcconvsnow,       &
                              stcsst, stcsss, stchice, stcaice, stcssp, &
                              stcwindspeedu, stcwindspeedv, outstctaux, &
                              outstctauy, outstcwme, outstcsol,         &
                              outstchtn, outstcple, outstcsno,          &
                              outstcsub, outstctop, outstcbot,          &
                              outstcsst, outstcsss, outstchice,         &
                              outstcssp, outstcwspx, outstcwspy, fftaux,&
                              fftauy, ffwme, ffsol, ffhtn, ffple, ffsno,&
                              ffsub, fftop, ffbot, ffsst, ffsss, ffhice,&
                              ffssp, ffwspx, ffwspy, pptaux, pptauy,    &
                              ppwme, ppsol, pphtn, ppple, ppsno, ppsub, &
                              pptop, ppbot, ppsst, ppsss, pphice, ppssp,&
                              ppwspx, ppwspy
      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,        &
                              len_realhd, itgrid, iugrid,               &
                              max_num_fc_times, max_num_clim_times,     &
                              max_num_in_flux, len2_lookuppreferred,    &
                              len2_lookupprevious, len2_lookupclimate
      implicit none

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

      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

      real lhevap               ! latent heat of evaporation
      parameter ( lhevap    = 2.25E6)

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
!----------------------------------------------------------------------
! comdeck: CDEBUG
! Purpose: declares and stores information needed to output debugging
!          diagnostics at user selected points as fields are processed.
!----------------------------------------------------------------------
! common block:
      COMMON / CDebug /                                                 &
     &  NoDbgPts,IColDbg,JRowDbg,l_winds_dbg,l_heat_dbg,l_moisture_dbg, &
     &  l_sea_ice_dbg, l_references_dbg, l_pressure_dbg,l_windspd_dbg

      common / CCDebug / CValues

! declarations:
      integer MaxNoDbgPts   ! parameter: max. number of points to output
      parameter ( MaxNoDbgPts = 20)

! points to output
      integer NoDbgPts      ! actual number of points to output
      integer IColDbg(MaxNoDbgPts)   ! column of each point
      integer JRowDbg(MaxNoDbgPts)   ! row    of each point

! character array for output
      CHARACTER(LEN=11) CValues(MaxNoDbgPts) ! values to write out

! debug logical for each output file
      LOGICAL :: l_winds_dbg
      LOGICAL :: l_heat_dbg
      LOGICAL :: l_moisture_dbg
      LOGICAL :: l_sea_ice_dbg
      LOGICAL :: l_references_dbg
      LOGICAL :: l_pressure_dbg
      LOGICAL :: l_windspd_dbg

! CDEBUG end

! declaration of local arrays
      integer Int_Head_SW1(Len_IntHd)   ! integer part of lookup table
      integer Int_Head_SW(Len_IntHd)    ! integer part of lookup table
      integer Int_Head_LW(Len_IntHd)    ! integer part of lookup table
      integer Int_Head_EVAP(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_SH(Len_IntHd)    ! integer part of lookup table
      real Real_Head_SW1(Len_RealHd)    ! real part of lookup table
      real Real_Head_SW(Len_RealHd)     ! real part of lookup table
      real Real_Head_LW(Len_RealHd)     ! real part of lookup table
      real Real_Head_EVAP(Len_RealHd)   ! real part of lookup table
      real Real_Head_SH(Len_RealHd)     ! real part of lookup table
      real SW_radiation_band1(ncols, nrowst)! short wave flux (band 1)
      real SW_radiation(ncols, nrowst)      ! short wave flux
      real LW_radiation(ncols, nrowst)      ! long_wave_radiation
      real evaporation(ncols, nrowst)       ! evaporation
      real sensible_heat(ncols, nrowst)     ! sensible heat
      real latent_heat(ncols,nrowst)        ! latent heat
      real non_pen_heat(ncols,nrowst)       ! net non-penetrating heat
      real fieldint(ncols,nrowst)           ! intermediate field
! declaration of local scalars

      integer ivt           ! loop index over validity times
      integer IVTOffHr      ! offset of validity time from reference
      integer IOutUnit      ! output unit

      logical ldebug        ! T => output debugging info (set in 0.)
      logical l_leads       ! T => using minleadsfrac
                            ! F => using minicefrac

      CHARACTER(LEN=256) cmessage   ! error message

! declaration of externals
      external read_leads_flds, write_one_field,                        &
     &         ScalarMult,FieldSub,FieldAdd

!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'heat'  ! subroutine name for error messages

      ldebug = l_heat_dbg      ! set by debug input control file

!----------------------------------------------------------------------
! 1. start loop over validity times
!----------------------------------------------------------------------
      do ivt = 1, NoValidTimes

        IVTOffHr = IValidOffHr(ivt)
        IOutUnit = IOutUnitOff(ivt) + UnitHeatOut

!----------------------------------------------------------------------
! 2. Read in net down short wave flux over open sea (band 1)
!----------------------------------------------------------------------
        l_leads = .true.                ! set to use minleadsfrac
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCSW1,StCAICE,                            &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_SW1,                         &
     &                    Real_Head_SW1, ncols, nrowst,                 &
     &                    SW_radiation_band1,                           &
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2. unable to read SW Radiation Flux (band 1)'
          icode = 1001
          go to 9999
        end if

! 2.2 Write out solar radiation
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
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
     &       OutStCSOL, FFSOL, PPSOL, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_SW1, Real_Head_SW1, IOutUnit,ldebug,              &
     &       SW_radiation_band1, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2. unable to write penetrating '                    &
     &       ,'solar radiation (SOL)'
          icode = 1101
          go to 9999
        end if

!----------------------------------------------------------------------
! 3. Read in fields to calculate net non penetrating heat
!----------------------------------------------------------------------
! 3.1 Read net down short wave readiation
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCSW,StCAICE,                             &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_SW,                          &
     &                    Real_Head_SW, ncols, nrowst,                  &
     &                    SW_radiation,                                 &
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read SW Radiation Flux'
          icode = 1002
          go to 9999
        end if
! 3.2 Use Field Sub to work out first component of HTN
! DEPENDS ON: fieldsub
        call FieldSub (ncols, nrowst, rmdi,                             &
     &            SW_radiation, SW_radiation_band1,                     &
     &            fieldint,                                             &
     &            icode, cmessage)
! 3.3 Read net down long wave flux
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCLongWave,StCAICE,                       &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_LW,                          &
     &                    Real_Head_LW, ncols, nrowst,                  &
     &                    LW_radiation,                                 &
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read LW Radiation Flux'
          icode = 1003
          go to 9999
        end if

! 3.4 Use FieldAdd to do HTN = fieldint + LW_radiation
! DEPENDS ON: fieldadd
        call FieldAdd (ncols, nrowst, rmdi,                             &
     &            fieldint, LW_radiation,                               &
     &            non_pen_heat,                                         &
     &            icode, cmessage)

! 3.5 Read evaporation from sea
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCEvaporation,StCAICE,                    &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_EVAP,                        &
     &                    Real_Head_EVAP, ncols, nrowst,                &
     &                    evaporation,                                  &
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read evaporation from sea'
          icode = 1004
          go to 9999
        end if

! 3.6 Use ScalarMult and FieldSub to work out latent heat
!     and subtract it from HTN
! DEPENDS ON: scalarmult
        call ScalarMult (ncols, nrowst, rmdi,                           &
     &            lhevap, evaporation,                                  &
     &            latent_heat,                                          &
     &            icode, cmessage)
! DEPENDS ON: fieldsub
        call FieldSub (ncols, nrowst, rmdi,                             &
     &            non_pen_heat, latent_heat,                            &
     &            fieldint,                                             &
     &            icode, cmessage)

! 3.7 Read Sensible Heat Flux
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCSensibleHeat,StCAICE,                   &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_SH,                          &
     &                    Real_Head_SH, ncols, nrowst,                  &
     &                    sensible_heat,                                &
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read sensible heat flux'
          icode = 1005
          go to 9999
        end if

! 3.8 Use FieldSub to calculate final HTN
! DEPENDS ON: fieldsub
        call FieldSub (ncols, nrowst, rmdi,                             &
     &            fieldint, sensible_heat,                              &
     &            non_pen_heat,                                         &
     &            icode, cmessage)

! 3.9 Write out net non penetrating heat
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
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
     &       OutStCHTN, FFHTN, PPHTN, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_SH, Real_Head_SH, IOutUnit,ldebug,                &

     &       non_pen_heat, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to write net non penetrating heat'
          icode = 1102
          go to 9999
        end if

!----------------------------------------------------------------------
! 4. end loop over validity times
!----------------------------------------------------------------------
        enddo    !  ivt

9999  continue
      return
      END SUBROUTINE heat
!----------------------------------------------------------------------
