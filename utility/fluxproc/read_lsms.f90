! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! deck: RDLSMS
!
! contains routines: read_lsms
!
! Purpose: reads land / sea masks, calculates coefficients for
!          interpolation between grids and indices for "unresolved"
!          seapoints on ocean grid which are not surrounded by seapoints
!          on the atmosphere grid.
!          Addition to handle rotated grids (S. Spall)
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_lsms(                                             &
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
     &           icode)

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,        &
                             len_realhd, itgrid, iugrid,                &
                             max_num_fc_times, max_num_clim_times,      &
                             max_num_in_flux, len2_lookuppreferred,     &
                             len2_lookupprevious, len2_lookupclimate
      USE lookup_addresses

      USE coast_aj_mod, ONLY: coast_aj
      USE eqtoll_mod,   ONLY: eqtoll
      USE h_int_co_mod, ONLY: h_int_co
      USE lltoeq_mod,   ONLY: lltoeq
      USE w_coeff_mod,  ONLY: w_coeff
      IMPLICIT NONE

! declaration of arguments

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

! Globals
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

! declaration of local arrays

! arrays required as input by h_int_co
       real phi_full(ncolsO*nrowstO)      ! { allocated with largest
       real lambda_full(ncolsO*nrowstO)   ! { row dimension

! temporary arrays of positions on atmosphere grid required for w_coeff
       real lambda_eqA(ncols*nrowsuv) ! Longitude on atmos  ==  grid
       real phi_eqA(ncols*nrowsuv)    ! Latitude on atmos  ==  grid
       real lambda_tmp1(ncols*nrowsuv) ! Long. on reg. lat-long grid
       real phi_tmp1(ncols*nrowsuv)    ! Latitude on reg. lat-long grid
       real lambda_tmp2(ncols*nrowsuv) ! Longitude on ocean eq. grid
       real phi_tmp2(ncols*nrowsuv)    ! Latitude on ocean eq. grid

! arrays required in the conversion of lat-long for ocean points
       real phi_eq(ncolsO*nrowstO)      ! {
       real lambda_eq(ncolsO*nrowstO)   ! { allocated with largest
       real phi_ll(ncolsO*nrowstO)      ! { row dimension
       real lambda_ll(ncolsO*nrowstO)   ! {

! arrays output by coast_aj which are not subsequently used
! they are allocated with largest ocean grid dimensions (i.e. tracer)
       integer index_targ(ncolsO*nrowstO)
       integer index_srce(ncols*nrowst)
!       integer coastal_points(ncolsO*nrowstO)
       integer coastal_points
       integer index_land_unres(ncolsO*nrowstO)

! declaration of local scalars
       logical mask    ! T => land sea mask is provided
       integer i       ! loop index for columns
       integer j       ! loop index for rows
       integer ij      ! loop index for points in 2D field

! scalar output by coast_aj which is not subsequently used
       integer n_pts_unres_land   ! number of unresolved land points

!----------------------------------------------------------------------

! 0. Preliminaries
      CSub = 'read_lsms'  ! subroutine name for error messages

!----------------------------------------------------------------------
! 1. Read land / sea masks
!----------------------------------------------------------------------

! 1.1 read atmosphere tracer land / sea mask and calculate
!    grid coordinates
! DEPENDS ON: read_lsm_anc
      call read_lsm_anc(UnitNWPlsmt, Len_FixHd, Len1_Lookup, FixHdlsmt, &
     &       Lookuplsmt, ncols, nrowst, lsmt, lambda_t, phi_t,          &
     &       icode)

! check icode
      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.1 failed reading NWP tracer land / sea mask'
        go to 9999
      end if

! 1.2 set atmosphere velocity land / sea mask from tracer land / sea
!     mask and calculate grid coordinates

! DEPENDS ON: set_lsmu
      call set_lsmu (  ncols, nrowst, nrowsuv, LCyclic,                 &
     &                 lambda_t, phi_t, lsmt, ILandPt,                  &
     &                 lsmu, lambda_u, phi_u )



! 1.3 read ocean tracer land / sea mask and calculate
!    grid coordinates
! DEPENDS ON: read_lsm_anc
      call read_lsm_anc(UnitFOAMlsmt, Len_FixHd, Len1_Lookup,           &
     &       FixHdlsmtO, LookuplsmtO, ncolsO, nrowstO, lsmtO,           &
     &       lambda_tO, phi_tO,                                         &
     &       icode)

! check icode
      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.3 failed reading ocean tracer land / sea mask'
        icode = icode + 2000
        go to 9999
      end if

! 1.4 read ocean velocity land / sea mask and calculate
!    grid coordinates
! DEPENDS ON: read_lsm_anc
      call read_lsm_anc(UnitFOAMlsmu, Len_FixHd, Len1_Lookup,           &
     &       FixHdlsmuO, LookuplsmuO, ncolsO, nrowsuO, lsmuO,           &
     &       lambda_uO, phi_uO,                                         &
     &       icode)

! check icode
      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.4 failed reading ocean velocity land / sea mask'
        icode = icode + 2000
        go to 9999
      end if

!----------------------------------------------------------------------
! 2. Find if the atmosphere grid is rotated
!---------------------------------------------------------------------

! 2.1 Get the position of the poles from the
!     atmosphere lsm header

! DEPENDS ON: copy_to_real
      call copy_to_real ( Lookuplsmt(BPLAT), pole_lat )
! DEPENDS ON: copy_to_real
      call copy_to_real ( Lookuplsmt(BPLON), pole_lon )

! 2.2 Find if a rotated grid is being used

      rotg=.true.
      if ( pole_lat  >   89.99 .and. pole_lat  <   90.01) then
        rotg=.false.
      end if

      if (pole_lat  <   -1.0e5) then
        rotg=.false.
      end if

! 2.3 Do error checking on the positions of the poles

      if (rotg) then

        if ( pole_lat  >   90.0 .or. pole_lat  <   -90.0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.3 incorrect latitude of pole in atmos lsm header'
          icode = icode + 2000
          go to 9999
        end if

        if ( pole_lon  >   360.0 .or. pole_lon  <   -360.0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.3 incorrect longitude of pole in atmos lsm header'
          icode = icode + 2000
          go to 9999
        end if

      end if

! 2.4 Write out details of the type of grid

      if (rotg) then
        write(UnStd,*)CStd//CSub//'Atmosphere on a rotated grid;',      &
     &    ' BPLAT = ', pole_lat, '; BPLON = ', pole_lon
      else
        write(UnStd,*)CStd//CSub//'Atmosphere on a non-rotated grid'
      end if

!----------------------------------------------------------------------
! 3. Find if the ocean grid is rotated
!---------------------------------------------------------------------

! 3.1 Get the position of the poles from the
!     ocean tracer grid lsm header

! DEPENDS ON: copy_to_real
      call copy_to_real ( LookuplsmtO(BPLAT), poleO_lat )
! DEPENDS ON: copy_to_real
      call copy_to_real ( LookuplsmtO(BPLON), poleO_lon )

! 3.2 Find if a rotated grid is being used

      rotgO=.true.
      if ( poleO_lat  >   89.99 .and. poleO_lat  <   90.01) then
        rotgO=.false.
      end if

      if (poleO_lat  <   -1.0e5) then
        rotgO=.false.
      end if

! 3.3 Do error checking on the positions of the poles

      if (rotgO) then

        if ( poleO_lat  >   90.0 .or. poleO_lat  <   -90.0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3.3 incorrect latitude of pole in ocean lsm header'
          icode = icode + 2000
          go to 9999
        end if

        if ( poleO_lon  >   360.0 .or. poleO_lon  <   -360.0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3.3 incorrect longitude of pole in ocean lsm header'
          icode = icode + 2000
          go to 9999
        end if

      end if

! 3.4 Write out details of the type of grid

      if (rotgO) then
        write(UnStd,*)CStd//CSub//'Ocean on a rotated grid;',           &
     &    ' BPLAT = ', poleO_lat, '; BPLON = ', poleO_lon
      else
        write(UnStd,*)CStd//CSub//'Ocean on a non-rotated grid'
      end if

!----------------------------------------------------------------------
! 4. Calculate interpolation coefficients for interpolation
!    from atmosphere to ocean tracer grids
!----------------------------------------------------------------------

! 4.1 prepare target grid coordinates required
!     by h_int_co for tracer grid

      do j = 1, nrowstO
        do i = 1, ncolsO
          ij = i + (j-1) * ncolsO
          phi_eq( ij ) = phi_tO ( j )
          lambda_eq( ij ) =  lambda_tO ( i )
       end do
      end do

! 4.2 If the ocean uses a rotated grid, convert the ocean lat-long
!     vector to a standard lat-long grid

      if (rotgO) then
        call eqtoll(phi_eq, lambda_eq, phi_ll, lambda_ll,               &
     &                      poleO_lat, poleO_lon, ncolsO*nrowstO)
      else
        do i = 1, ncolsO*nrowstO
          phi_ll( i ) = phi_eq( i )
          lambda_ll( i ) =  lambda_eq ( i )
        end do
      end if

! 4.3 If the atmosphere uses a rotated grid, convert the ocean standard
!     lat-long to the atmosphere rotated grid

      if (rotg) then
        call lltoeq(phi_ll, lambda_ll, phi_full, lambda_full,           &
     &                      pole_lat, pole_lon, ncolsO*nrowstO)
      else
        do i = 1, ncolsO*nrowstO
          phi_full( i ) = phi_ll( i )
          lambda_full( i ) =  lambda_ll ( i )
        end do
      end if

! 4.4 Convert target longitude to correct range

      do i = 1, ncolsO*nrowstO
        lambda_full(i)=mod(lambda_full(i)-lambda_t(1)+720.,360.)        &
     &                         +lambda_t(1)
      end do

! 4.5 Calculate interpolation coefficients for tracer grids

      call h_int_co(index_bl_t,index_br_t,                              &
     & weight_tr_t,weight_br_t,weight_tl_t,weight_bl_t,                 &
     & lambda_t, phi_t, lambda_full, phi_full,                          &
     & ncols,nrowst,ncolsO*nrowstO,LCyclic)

!----------------------------------------------------------------------
! 5. Calculate interpolation coefficients for interpolation
!    from atmosphere to ocean velocity grids. Atmosphere velocity
!    components must be defined at coincident points.
!----------------------------------------------------------------------

! 5.1 prepare target grid coordinates required
!     by h_int_co for velocity grid

      do j = 1, nrowsuO
        do i = 1, ncolsO
          ij = i + (j-1) * ncolsO
          phi_eq( ij ) = phi_uO ( j )
          lambda_eq (ij ) =  lambda_uO ( i )
       end do
      end do

! 5.2 If the ocean uses a rotated grid, convert the ocean lat-long
!     vector to a standard lat-long grid

      if (rotgO) then
        call eqtoll(phi_eq, lambda_eq, phi_ll, lambda_ll,               &
     &                      poleO_lat, poleO_lon, ncolsO*nrowsuO)
      else
        do i = 1, ncolsO*nrowsuO
          phi_ll( i ) = phi_eq( i )
          lambda_ll( i ) =  lambda_eq ( i )
        end do
      end if

! 5.3 If the atmosphere uses a rotated grid, convert the ocean standard
!     lat-long to the atmosphere rotated grid

      if (rotg) then
        call lltoeq(phi_ll, lambda_ll, phi_full, lambda_full,           &
     &                      pole_lat, pole_lon, ncolsO*nrowsuO)
      else
        do i = 1, ncolsO*nrowsuO
          phi_full( i ) = phi_ll( i )
          lambda_full( i ) =  lambda_ll ( i )
        end do
      end if

! 5.4 Convert target longitude to correct range

      do i = 1, ncolsO*nrowsuO
        lambda_full(i)=mod(lambda_full(i)-lambda_u(1)+720.,360.)        &
     &                         +lambda_u(1)
      end do

! 5.5 Calculate interpolation coefficients for velocity grids
      call h_int_co(index_bl_u,index_br_u,                              &
     & weight_tr_u,weight_br_u,weight_tl_u,weight_bl_u,                 &
     & lambda_u, phi_u, lambda_full, phi_full,                          &
     & ncols,nrowsuv,ncolsO*nrowsuO,LCyclic)

!----------------------------------------------------------------------
! 6. Calculate the coefficients for rotating wind vectors to align
!    with the ocean grid
!----------------------------------------------------------------------

! 6.1 Set up the coefficients for atmosphere to reg. lat-long

      do j = 1, nrowsuv
        do i = 1, ncols
          ij = i + (j-1) * ncols
          phi_eqA( ij ) = phi_u ( j )
          lambda_eqA ( ij ) =  lambda_u ( i )
       end do
      end do

      if (rotg) then
        call eqtoll(phi_eqA, lambda_eqA, phi_tmp1, lambda_tmp1,         &
     &                      pole_lat, pole_lon, ncols*nrowsuv)
        call w_coeff(coef_angle1, coef_angle2, lambda_tmp1,             &
     &             lambda_eqA, pole_lat, pole_lon, ncols*nrowsuv)
      else
        do ij = 1, ncols*nrowsuv
          lambda_tmp1(ij)=lambda_eqA(ij)
          phi_tmp1(ij)=phi_eqA(ij)
        enddo
      endif

! 6.2 Set up the coefficients for reg. lat-long to ocean

      if (rotgO) then
        call lltoeq(phi_tmp1, lambda_tmp1, phi_tmp2, lambda_tmp2,       &
     &                      poleO_lat, poleO_lon, ncols*nrowsuv)
        call w_coeff(coef_angle3, coef_angle4, lambda_tmp1,             &
     &             lambda_tmp2, poleO_lat, poleO_lon, ncols*nrowsuv)
      endif

!----------------------------------------------------------------------
! 7. Calculate indices for unresolved points i.e. seapoints on ocean
!    grid which are not surrounded by seapoints on the atmosphere grid
!----------------------------------------------------------------------

! 7.1 Calculate indices for unresolved points for tracer grids

      mask = .true.  ! land / sea mask for target grid is to be used

      call coast_aj (index_bl_t,index_br_t,                             &
     & weight_tr_t,weight_br_t,weight_tl_t,weight_bl_t,                 &
     & ncols,nrowst,ncolsO*nrowstO,                                     &
     & lsmt,lsmtO,                                                      &
     & index_targ,index_srce,coastal_points,mask,                       &
     & index_unres_t,n_pts_unres_t,                                     &
     & index_land_unres,n_pts_unres_land)

! 7.2 Calculate indices for unresolved points for velocity grids

      call coast_aj (index_bl_u,index_br_u,                             &
     & weight_tr_u,weight_br_u,weight_tl_u,weight_bl_u,                 &
     & ncols,nrowsuv,ncolsO*nrowsuO,                                    &
     & lsmu,lsmuO,                                                      &
     & index_targ,index_srce,coastal_points,mask,                       &
     & index_unres_u,n_pts_unres_u,                                     &
     & index_land_unres,n_pts_unres_land)

!----------------------------------------------------------------------
! 8. Determine number of searchs needed to fill in unresolved points
!----------------------------------------------------------------------

! 8.1 on tracer grid
! DEPENDS ON: set_searches
      call set_searches ( ncolsO, nrowstO, LCyclicO, ISeaPt,            &
     &     lsmtO, n_pts_unres_t, index_unres_t, max_no_searches,        &
     &     n_calls_spiral_t, n_pts_spiral_t, icode)

      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 4.1 unable to fill in all unresolved points '
        icode = icode + 2000
        go to 9999
      end if


! 8.2 on velocity grid
! DEPENDS ON: set_searches
      call set_searches ( ncolsO, nrowsuO, LCyclicO, ISeaPt,            &
     &     lsmuO, n_pts_unres_u, index_unres_u, max_no_searches,        &
     &     n_calls_spiral_u, n_pts_spiral_u, icode)

      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 4.2 unable to fill in all unresolved points '
        icode = icode + 2000
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE read_lsms
!----------------------------------------------------------------------
