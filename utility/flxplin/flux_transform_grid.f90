! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! contains routines: Flux_Transform_Grid
!
! Purpose: Flux processing routine.
!          Controls main processing for Flux_Transform_Main
!
!    Programming standard : UMDP 3
!
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine Flux_Transform_Grid (                                  &
!----------------------------------------------------------------------
! comdeck: AFLDDIMS
! Purpose: argument list for dimensions of fields.
!          This deck is linked to FLDDIMS.
!----------------------------------------------------------------------
     & ncols, nrowst, nrowsu, nrowsv, nrowsuv, ncolsO, nrowstO, nrowsuO,&
!----------------------------------------------------------------------
                                      NoInFiles,icode)

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
                             len_realhd, itgrid, iugrid,                &
                             max_num_fc_times, max_num_clim_times,      &
                             max_num_in_flux, len2_lookuppreferred,     &
                             len2_lookupprevious, len2_lookupclimate
      USE lookup_addresses
      USE w_eqtoll_mod, ONLY: w_eqtoll
      USE w_lltoeq_mod, ONLY: w_lltoeq
 
     IMPLICIT NONE

! declaration of argument list
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
      integer NoInFiles ! IN number of input files
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! no local parameters

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

! declaration of local arrays (all arrays in COMDECK CFIELDS)
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

! local fields
      integer Int_Head(Len_IntHd)  ! integer header
      real Real_Head(Len_RealHd)   ! real header

      integer Int_Head_y(Len_IntHd) ! headers for y component
      real Real_Head_y(Len_RealHd)   ! of wind type fields

      real scalar(ncols, nrowst) ! scalar field

      real windx(ncols, nrowsuv)  ! wind field - x component
      real windy(ncols, nrowsuv)  ! wind field - x component
      real windx_tmp(ncols, nrowsuv) ! partially rotated
      real windy_tmp(ncols, nrowsuv) ! wind fields

! local scalars
      integer IROW_NUMBER
      CHARACTER(LEN=80) cmessage
      integer IStC   ! stash code
      integer i, j, input_field, ifile   ! loop counters
      logical ldebug
      logical L_more_fields

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'Flux_Transform_Grid'  ! subroutine name for err messages

! 0.1 Read StashMaster files

      IROW_NUMBER=0
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_A',IROW_NUMBER,                     &
     &  ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_O',IROW_NUMBER,                     &
     &  ICODE,CMESSAGE)

! 1. Read in land sea masks and calculate grid coordinates and
!    coefficients for interpolation from atmosphere to ocean grids.

! DEPENDS ON: read_lsms
      call read_lsms (                                                  &
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
     &                icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1. error in read_lsms '
        go to 9999
      end if

! 2. Loop over fields to be read and written

! need to avoid unit 22 which is assigned to STASHMaster
      do ifile=IUnOutLow+10, IUnOutLow+NoInFiles+9

       print*,'processing file number ',ifile
       input_field = 0
       L_more_fields = .true.

       do while ( L_more_fields )
         input_field = input_field + 1

! 2.1 read the next field header

        read (ifile, IOStat = icode) Int_Head, Real_Head


        if ( icode  <   0 ) then
          write(UnWarn,*)CWarn,CSub,                                    &
     &       ' step 2.1 End of file ',ifile,' reached after',           &
     &       input_field-1,' fields'
          icode=0
          L_more_fields = .false.
        end if


        if ( L_more_fields ) then
          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &         ' step 2.1 Error reading ', input_field,'th field'
            go to 9999
          end if


          IStC = Int_Head(ITEM_CODE)

          print*,'STASH code for field ',input_field, ' is ',IStC
          print*,'scalar dimensions are ', ncols, nrowst

! 2.2 IF it is not a wind type field THEN

        if ( IStC  /=  StCWindSpeedU  .and. IStC  /=  OutStCTAUX        &
     & .and. IStC  /=  OutStCTAUY .and. IStC  /=  StCWindSpeedV )       &
     &  then

! 2.2.1 read in the scalar field

          read (ifile, IOStat = icode) scalar

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &         ' step 2.2.1 Error reading in ',                         &
     &         input_field, 'th field '
            go to 9999
          end if

! 2.2.2 call routine to interpolate to new grid and write it out

! DEPENDS ON: field_interpolate_write
          call field_interpolate_write(                                 &
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
     &        Int_Head, Real_Head, ldebug, ITGrid, nrowst,              &
     &        IUnOutLow+NoInFiles+10, scalar, icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' Step 2.2.2 error outputting scalar field '
            go to 9999
          end if

! 2.3  ELSE

        else ! not a scalar field

! 2.3.1 read in x component wind field

          read (ifile, IOStat = icode) windx

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &         ' step 2.3.1 Error reading in ',                         &
     &         input_field, 'th field - which is a wind field '
            go to 9999
          end if

! 2.3.2 read in the next field header

          read (ifile, IOStat = icode) Int_Head_y, Real_Head_y

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 2.3.2 Error reading header of ',                    &
     &       input_field, 'th field header '
            go to 9999
          end if

! 2.3.3 check that it is y component of a velocity type field

          IStC=Int_Head_y(ITEM_CODE)
          if (      IStC  /=  OutStCTAUY                                &
     &        .and. IStC  /=  StCWindSpeedV ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 2.3.3 Header must be y component wind field ',      &
     &       'STASH code is ',IStC
            go to 9999
          end if

! 2.3.4 read in the y component wind field

          read (ifile, IOStat = icode) windy

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 2.3.4 Error reading in ',                           &
     &       input_field, 'th field - which is a wind field '
            go to 9999
          end if

! 2.3.5 rotate wind components if necessary
          if (rotg) then
           call w_eqtoll(coef_angle1, coef_angle2, windx,               &
     &           windy, windx_tmp, windy_tmp, ncols*nrowsuv,            &
     &           .TRUE.)
          else
           do j = 1, nrowsuv
             do i = 1, ncols
               windx_tmp(i,j)=windx(i,j)
               windy_tmp(i,j)=windy(i,j)
             enddo
           enddo
          endif

          if (rotgO) then
           call w_lltoeq(coef_angle3, coef_angle4, windx_tmp,           &
     &           windy_tmp, windx, windy, ncols*nrowsuv,                &
     &           .TRUE.)
          else
           do j = 1, nrowsuv
             do i = 1, ncols
               windx(i,j)=windx_tmp(i,j)
               windy(i,j)=windy_tmp(i,j)
             enddo
           enddo
          endif

! 2.3.6 Output both fields

! DEPENDS ON: field_interpolate_write
          call field_interpolate_write(                                 &
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
     &        Int_Head, Real_Head, ldebug, IUGrid, nrowsuv,             &
     &        IUnOutLow+NoInFiles+10, windx, icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' Step 2.3.6 error outputting wind-x field '
            go to 9999
          end if

! DEPENDS ON: field_interpolate_write
          call field_interpolate_write(                                 &
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
     &        Int_Head_y, Real_Head_y, ldebug, IUGrid, nrowsuv,         &
     &        IUnOutLow+NoInFiles+10, windy, icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' Step 2.3.6 error outputting wind-y field '
            go to 9999
          end if

        end if !   IStC  /=  wind codes

        end if ! L_more_fields

       end do
       print*,'finished processing file ',ifile

      end do

9999  continue
      return
      END SUBROUTINE Flux_Transform_Grid
!----------------------------------------------------------------------
