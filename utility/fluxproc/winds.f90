! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: winds
!
! Purpose: To produce a pp field containing:
!            wind stress (x-direction)
!            wind stress (y-direction)
!            wind mixing energy
!          for all the fields required.
!          Addition of rotated grid (S. Spall)
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine winds(                                                 &
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
      USE lookup_addresses

      USE w_eqtoll_mod, ONLY: w_eqtoll
      USE w_lltoeq_mod, ONLY: w_lltoeq
      IMPLICIT NONE

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
!----------------------------------------------------------------------
! comdeck: CSELCT
! Purpose: declares and stores information needed for flux selectiong
!----------------------------------------------------------------------
! common block:
      common / CSelect /                                                &
     &   l_B_grid,                                                      &
     &   l_winds_slt,   l_heat_slt,       l_moisture_slt,               &
     &   l_sea_ice_slt, l_references_slt, l_pressure_slt,               &
     &   l_windspd_slt

! debug logical for each selected flux
      logical l_B_grid,                                                 &
     &        l_winds_slt,   l_heat_slt,      l_moisture_slt,           &
     &        l_sea_ice_slt, l_references_slt,l_pressure_slt,           &
     &        l_windspd_slt
!----------------------------------------------------------------------
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


! declaration of local arrays
      integer Int_Head_taux(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_tauy(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_WME(Len_IntHd)  ! integer part of lookup table
      real Real_Head_taux(Len_RealHd)   ! real part of lookup table
      real Real_Head_tauy(Len_RealHd)   ! real part of lookup table
      real Real_Head_WME(Len_RealHd)   ! real part of lookup table
! windstresses on true atmosphere grids
      real winds_u(ncols,nrowsu)
      real winds_v(ncols,nrowsv)
! windstresses with components calculated at coincident points
      real windstressu(ncols, nrowsuv)   ! wind stress fields (taux)
      real windstressv(ncols, nrowsuv)   ! wind stress fields (tauy)
      real wind_mixing_energy(ncols,nrowst) ! wind mixing energy
      real wndu_tmp(ncols, nrowsuv)   ! wind str. on reg. lat-long gri
      real wndv_tmp(ncols, nrowsuv)   ! wind str. on reg. lat-long gri

! declaration of local scalars

      integer ivt           ! loop index over validity times
      integer IVTOffHr      ! offset of validity time from reference
      integer IOutUnit      ! output unit

      integer i       ! loop index for columns
      integer j       ! loop index for rows

      logical ldebug        ! T => output debugging info (set in 0.)
      logical l_leads       ! T => using minleadsfrac
                            ! F => using minicefrac

! declaration of externals
      external read_fields, write_one_field, read_leads_flds

!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'winds'  ! subroutine name for error messages

      ldebug = l_winds_dbg     ! set by debug input control file

!----------------------------------------------------------------------
! 1. start loop over validity times
!----------------------------------------------------------------------
      do ivt = 1, NoValidTimes

        IVTOffHr = IValidOffHr(ivt)
        IOutUnit = IOutUnitOff(ivt) + UnitWindsOut

!----------------------------------------------------------------------
! 2.1 read in windstresses and convert to B grid values if necessary
!----------------------------------------------------------------------

! DEPENDS ON: read_fields
        call read_fields(StCWindStressU, IVTOffHr,                      &
     &                 ldebug, Int_Head_taux, Real_Head_taux,           &
     &                 ncols, nrowsu,                                   &
     &                 winds_u,                                         &
     &                 icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.1 unable to read u windstresses'
          icode = 1006
          go to 9999
        end if

! DEPENDS ON: read_fields
        call read_fields(StCWindStressV, IVTOffHr,                      &
     &                 ldebug, Int_Head_tauy, Real_Head_tauy,           &
     &                 ncols, nrowsv,                                   &
     &                 winds_v,                                         &
     &                 icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.1 unable to read v windstresses'
          icode = 1006
          go to 9999
        end if

        if ( l_b_grid ) then
          do j=1, nrowsuv
            do i = 1, ncols
              windstressu(i,j) = winds_u(i,j)
              windstressv(i,j) = winds_v(i,j)
            end do
          end do

        else   ! l_b_grid
          do j=1, nrowsuv
            do i = 1, ncols-1
              windstressu(i,j) = 0.5*(winds_u(i,j)+winds_u(i,j+1))
              windstressv(i,j) = 0.5*(winds_v(i,j)+winds_v(i+1,j))
            end do
            if (LCyclic) then
              windstressu(ncols,j)=windstressu(1,j)
              windstressv(ncols,j)=windstressv(1,j)
            else
              windstressu(ncols,j)=0.0 ! missing data
              windstressv(ncols,j)=0.0 ! missing data
            end if
          end do

! Update Lookup tables to reflect changes to grid
          Real_Head_taux(BZY-Len_IntHd)=                                &
     &          Real_Head_tauy(BZY-Len_IntHd)
          Real_Head_tauy(BZX-Len_IntHd)=                                &
     &          Real_Head_taux(BZX-Len_IntHd)

        end if  ! l_b_grid


!----------------------------------------------------------------------
! 2.2 Rotate wind vectors if rotated grids are used
!----------------------------------------------------------------------

      if (rotg) then
        call w_eqtoll(coef_angle1, coef_angle2, windstressu,            &
     &           windstressv, wndu_tmp, wndv_tmp, ncols*nrowsuv,        &
     &           .TRUE.)
      else
        do j = 1, nrowsuv
          do i = 1, ncols
            wndu_tmp(i,j)=windstressu(i,j)
            wndv_tmp(i,j)=windstressv(i,j)
          enddo
        enddo
      endif

      if (rotgO) then
        call w_lltoeq(coef_angle3, coef_angle4, wndu_tmp,               &
     &           wndv_tmp, windstressu, windstressv, ncols*nrowsuv,     &
     &           .TRUE. )
      else
        do j = 1, nrowsuv
          do i = 1, ncols
            windstressu(i,j)=wndu_tmp(i,j)
            windstressv(i,j)=wndv_tmp(i,j)
          enddo
        enddo
      endif

!----------------------------------------------------------------------
! 2.3 write out U + V component of windstress
!----------------------------------------------------------------------

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
     &       OutStCTAUX, FFTAUX, PPTAUX, IVTOffHr,                      &
     &          IUGrid, nrowsuv,                                        &
     &          Int_Head_taux, Real_Head_taux, IOutUnit, ldebug,        &
     &       windstressu, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.2 unable to write U component of windstresses'
          icode = 1103
          go to 9999
        end if
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
     &       OutStCTAUY, FFTAUY, PPTAUY, IVTOffHr,                      &
     &          IUGrid, nrowsuv,                                        &
     &          Int_Head_tauy, Real_Head_tauy, IOutUnit, ldebug,        &
     &       windstressv, icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.2 unable to write V component of windstresses'
          icode = 1104
          go to 9999
        end if

!----------------------------------------------------------------------
! 3. Read in wind mixing energy
!----------------------------------------------------------------------
        l_leads = .true.                ! set to use minleadsfrac
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCWindMixEng, StCAICE,                    &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_WME,                         &
     &                    Real_Head_WME, ncols, nrowst,                 &
     &                    wind_mixing_energy,                           &
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read Wind Mixing Energy'
          icode = 1007
          go to 9999
        end if

!----------------------------------------------------------------------
! 3.1 Write out wind mixing energy
!----------------------------------------------------------------------
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
     &       OutStCWME, FFWME, PPWME, IVTOffHr,                         &
     &          ITGrid, nrowst,                                         &
     &          Int_Head_WME, Real_Head_WME, IOutUnit, ldebug,          &
     &       wind_mixing_energy, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3.1. unable to write wind mixing energy'
          icode = 1105
          go to 9999
        end if

!----------------------------------------------------------------------
! Last. end loop over validity times
!----------------------------------------------------------------------
        enddo    !  ivt

9999  continue
      return
      END SUBROUTINE winds
!----------------------------------------------------------------------
