! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!-----SUBROUTINE PRINT_INITDATA---------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!     **************************************************************
!     *                                                            *
!     *  Single Column Unified Model version.                      *
!     *                                                            *
!     **************************************************************
!
!
!     Purpose: To print out initial run data for SCM integration
!
!---------------------------------------------------------------------

SUBROUTINE print_initdata                                                     &
  ( row_length,rows, land_pts, nlevs, nwet, nozone, nfor, nbl_levs            &
  , nsoilt_levs, nsoilm_levs, ntrop, dayno_init, a_sw_radstep_diag            &
  , a_sw_radstep_prog, tls, qls, uls, vls, wls, ichgf, ilscnt, ti, time_init  &
  , nout, numnout )


  USE scm_utils, ONLY:                                                        &
    zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_main_force, ONLY:                                                     &
    year_init, ndayin, nminin, nsecin, timestep, lat, long, ancyc, local_time &
  , change_clim, exname_in, runno_in, exname_out, runno_out, soil_type        &
  , land_sea_mask, obs, geoforce, geoinit, stats, noforce, tapein, tapeout    &
  , smi_opt, smcli, sth, fsmc, ug, vg, conv_mode, altdat, t_inc               &
  , tstar_forcing, q_star, u_inc, v_inc, w_inc, flux_h, flux_e, ui, vi, wi    &
  , qi, ccai, iccbi, iccti, p_in, canopy_gbi, smci, snodepi, t_deep_soili     &
  , tstari, z0mseai, ntrad1, radcloud_fixed, cca_rad, iccb_rad, icct_rad      &
  , ccwpin_rad, layer_cloud_rad, qcl_rad, ozone, ug_opt, vg_opt


  IMPLICIT NONE

!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------


  INTEGER ::                    &
    row_length                  &! In leading dimension of SCM arrays.
  , rows                        &
  , nlevs                       &! In no of levels.
  , nwet                        &! In no of model levels in which Q is set.
  , nfor                        &! In Number terms for observational forcing
  , nbl_levs                    &! In Number of Boundary layer levels
  , nsoilt_levs                 &! In Number of soil temperature levels
  , nsoilm_levs                 &! In Number of soil moisture levels
  , ntrop                       &! In Max number of levels in the troposphere
  , nozone                      &! In Model levels in which ozone is set.
  , land_pts                     ! In no of land points


  INTEGER ::                    &
    dayno_init                  &! Initial day (assume 360 day year)
  , ichgf                       &! No. of timesteps between change in
                                 ! observational forcing
  , ilscnt                      &! Counts for change in observational forcing
  , a_sw_radstep_diag           &! No. of diagnostic timesteps
  , a_sw_radstep_prog           &! No. of prognostic timesteps
                                 ! between calls to radiation
  , numnout                     &! No. units for output
  , nout(numnout)                ! Units for output


! Increments due to large-scale horizontal and vertical advection
  REAL ::                           &
    tls(row_length,rows,nfor,nlevs) &! Temperature
  , qls(row_length,rows,nfor,nwet)  &! Specific humidity
  , uls(row_length,rows,nfor,nlevs) &! Zonal wind
  , vls(row_length,rows,nfor,nlevs) &! Meridional wind
  , wls(row_length,rows,nfor,0:nlevs)  ! Vertical velocity


  REAL ::                           &
    ti(row_length,rows,nlevs)        ! Initial temp. profile  (K)



  REAL ::      &
    time_init   ! Start time of run(secs)


! Local Variables

  INTEGER ::                 &
    i, j, k, l, m, land_cnt  &! Counters
  , multnl,resnl             &! Counters for printing nlevs variables
  , multnf,resnf             &! Counter for printing nfor variables
  , multsl,ressl             &! Counter for printing nsoil variables
  , multoz,resoz             &! Counter for printing nozone variables
  , multcl,rescl

  CHARACTER (Len=100) :: c1fmt, c2fmt, c3fmt, c4fmt, c5fmt
  CHARACTER (Len=100) :: c6fmt, c7fmt, c8fmt, c10fmt
  CHARACTER (Len=100) :: c11fmt1, c11fmt2
  CHARACTER (Len=100) :: c12fmt, c13fmt, c14fmt, c15fmt
  CHARACTER (Len=100) :: c16fmt, c17fmt, c18fmt, c19fmt, c20fmt
  CHARACTER (Len=100) :: c31fmt, c32fmt, c33fmt, c34fmt, c35fmt
  CHARACTER (Len=100) :: c36fmt, c37fmt, c38fmt, c39fmt, c40fmt
  CHARACTER (Len=100) :: c41fmt, c42fmt, c43fmt, c44fmt, c45fmt
  CHARACTER (Len=100) :: c46fmt, c47fmt, c48fmt, c49fmt, c50fmt
  CHARACTER (Len=100) :: c51fmt, c52fmt, c53fmt, c54fmt, c55fmt
  CHARACTER (Len=100) :: c56fmt, c57fmt, c58fmt, c59fmt, c60fmt
  CHARACTER (Len=100) :: c61fmt, c62fmt, c63fmt, c64fmt, c65fmt
  CHARACTER (Len=100) :: c66fmt, c67fmt, c68fmt, c69fmt, c70fmt
  CHARACTER (Len=100) :: c71fmt, c72fmt, c73fmt, c74fmt, c75fmt
  CHARACTER (Len=100) :: c76fmt, c77fmt, c78fmt, c79fmt, c80fmt
  CHARACTER (Len=100) :: c81fmt, c82fmt, c83fmt, c84fmt, c85fmt
  CHARACTER (Len=100) :: c86fmt, c87fmt, c88fmt, c89fmt, c90fmt
  CHARACTER (Len=100) :: c91fmt, c92fmt, c93fmt, c94fmt, c95fmt
  CHARACTER (Len=100) :: c96fmt, c97fmt, c98fmt, c99fmt, c100fmt
  CHARACTER (Len=100) :: c101fmt, c102fmt, c103fmt

  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('PRINT_INITDATA',zhook_in,zhook_handle)

!---------------------------------------------------------------------
!     Set up format statements
!---------------------------------------------------------------------

  c1fmt = '('                                                       &
    //'" * Initial data used for scm integration      ",27x,"*")'
  c2fmt = '('                                                       &
   //'"******************************************************"'     &
   //',"*******************")'
  c3fmt = '('                                                       &
   //'" * Model run is for a land point              ",26x," *")'
  c4fmt = '('                                                       &
   //'" * Model run is for a sea point               ",26x," *")'
  c5fmt = '('                                                       &
   //'" * Run started on year                        ",i26," *")'
  c6fmt = '('                                                       &
   //'" * Run started on day                         ",i26," *")'
  c7fmt = '('                                                       &
   //'" * Run started at time (secs)               ",f28.2," *")'
  c8fmt = '('                                                       &
   //'" * length of run (days:mins:secs)",13x,i8,":",i8,":",'       &
   //'i8, " *")'
  c10fmt = '('                                                      &
   //'" * Timestep (seconds)                       ",f28.1," *")'
  c11fmt1 = '('                                                     &
   //'" * Prognostic Timestep for radiation (secs) ",f28.1," *")'
  c11fmt2 = '('                                                     &
   //'" * Diagnostic Timestep for radiation (secs) ",f28.1," *")'
  c12fmt = '('                                                      &
   //'" * 1st timestep for radiation (seconds)     ",f28.1," *")'
  c13fmt = '('                                                      &
   //'" * Latitude of gridpoint                    ",f28.1," *")'
  c14fmt = '('                                                      &
   //'" * Longitude of gridpoint                   ",f28.1," *")'
  c15fmt = '('                                                      &
   //'" * Number of model levels                     ",i26," *")'
  c16fmt = '('                                                      &
   //'" * Number of atmospheric moist levels         ",i26," *")'
  c17fmt = '('                                                      &
   //'" * Number of levels in boundary layer         ",i26," *")'
  c18fmt = '('                                                      &
   //'" * Number of deep soil temperatures           ",i26," *")'
  c19fmt = '('                                                      &
   //'" * Convection switched off                    ",26x," *")'
  c20fmt = '('                                                      &
   //'" * Convection switched on for diagnostics only",26x," *")'
  c31fmt = '('                                                      &
   //'" * Annual cycle included                      ",26x," *")'
  c32fmt = '('                                                      &
   //'" * Annual cycle not included                  ",26x," *")'
  c33fmt = '('                                                      &
   //'" * All diagnostics refer to local time        ",26x," *")'
  c34fmt = '('                                                      &
   //'" * All diagnostics refer to gmt               ",26x," *")'
  c35fmt = '('                                                      &
   //'" * Large-scale observational forcing selected ",26x," *")'
  c36fmt = '('                                                      &
   //'" * Number of timesteps between change in forcing",i25,"*")'
  c37fmt = '('                                                      &
   //'" * Counter for forcing                        ",i26," *")'
  c38fmt = '('                                                      &
   //'" * Number of terms for observational forcing  ",i26," *")'
  c39fmt = '('                                                      &
   //'" * Temperature increment tls (k)              ",26x," *")'
  c40fmt = '('                                                      &
   //'" * Zonal wind increment uls     (m/s)         ",26x," *")'
  c41fmt = '('                                                      &
   //'" * Meridional wind increment vls  (m/s)       ",26x," *")'
  c42fmt = '('                                                      &
   //'" * Specific humidity increment qls (kg/kg)    ",26x," *")'
  c43fmt = '('                                                      &
   //'" * Flux_h - observational forcing             ",26x," *")'
  c44fmt = '('                                                      &
   //'" * Flux_e - observational forcing             ",26x," *")'
  c45fmt = '('                                                      &
   //'" * Large-scale statistical forcing selected   ",26x," *")'
  c46fmt = '('                                                      &
   //'" * Number of days between change in forcing   ",i26," *")'
  c47fmt = '('                                                      &
   //'" * No large-scale forcing selected            ",26x," *")'
  c48fmt = '('                                                      &
   //'" * Geostrophic forcing selected               ",26x," *")'
  c49fmt = '('                                                      &
   //'" * Geostrophic initialisation selected        ",26x," *")'
  c50fmt = '('                                                      &
   //'" * Initial zonal wind profile ui  (m/s)       ",26x," *")'
  c51fmt = '('                                                      &
   //'" * Initial meridional wind profile vi  (m/s)  ",26x," *")'
  c52fmt = '('                                                      &
   //'" * Initial geostrophic zonal wind ug (m/s)"'                 &
   //',20x,f10.3," *")'
  c53fmt = '('                                                      &
   //'" * Initial geostrophic meridional wind vg (m/s)",15x,'       &
   //' f10.3," *")'
  c54fmt = '('                                                      &
   //'" * Initial temperature profile ti  (k)        ",26x," *")'
  c55fmt = '('                                                      &
   //'" * Initial specific humidity profile qi (kg/kg)",26x,"*")'
  c56fmt = '('                                                      &
   //'" * Ozone                                      ",26x," *")'
  c57fmt = '('                                                      &
   //'" * Initial convective cloud amount ccai     ",f28.4," *")'
  c58fmt = '('                                                      &
   //'" * Initial convective cloud base iccbi        ",i26," *")'
  c59fmt = '('                                                      &
   //'" * Initial convective cloud top iccti         ",i26," *")'
  c60fmt = '('                                                      &
   //'" * Initial pressure p_in     ",10f7.5," *")'
  c61fmt = '('                                                      &
   //'" * Fixed radiation conv. cloud amount cca_rad",f27.4," *")'
  c62fmt = '('                                                      &
   //'" * Fixed radiation conv. cloud base iccb_rad  ",i26," *")'
  c63fmt = '('                                                      &
   //'" * Fixed radiation conv. cloud top icct_rad   ",i26," *")'
  c64fmt = '('                                                      &
   //'" * Fixed radiation conv. cloud ccwpin_rad    ",f27.4," *")'
  c65fmt = '('                                                      &
   //'" * Fixed radiation layer cloud amount          ",25x," *")'
  c66fmt = '('                                                      &
   //'" * Fixed radiation total water and ice in cloud",25x," *")'
  c67fmt = '('                                                      &
   //'" * Initial canopy water content kg m-2      ",f28.3," *")'
  c68fmt = '('                                                      &
   //'" * Initial soil moisture content kg m-2     ",f28.3," *")'
  c69fmt = '('                                                      &
   //'" * Initial snow depth kg m-2                ",f28.3," *")'
  c70fmt = '('                                                      &
   //'" * Initial deep soil temperatures (k)",6x,"*")'
  c71fmt = '('                                                      &
   //'" * Initial surface temperature (K)          ",f28.3," *")'
  c72fmt = '('                                                      &
   //'" * JULES initialised by input of actual soil "'              &
   //'"moisture in ,"layers        *")'
  c73fmt = '('                                                      &
   //'" * Initial smcl  (Kg m-2)                      ",26X,"*")'
  c74fmt = '('                                                      &
   //'"* JULES initialisation by input of soil moisture ",'         &
   //'"stress factor    *")'
  c75fmt = '('                                                      &
   //'" * Initial soil moisture stress factor      ",f28.3," *")'
  c76fmt = '('                                                      &
   //'" * JULES initialised by input of smcl as a fraction of ",'   &
   //'"saturation        *")'
  c77fmt = '('                                                      &
   //'" * Initial sth                                 ",26x,"*")'
  c78fmt = '('                                                      &
   //'" * Initial sea surface roughness length     ",f28.3," *")'
  c79fmt = '('                                                      &
   //'" * Initial data read from tape                ",26x," *")'
  c80fmt = '('                                                      &
   //'" * With experment name                        ",a26," *")'
  c81fmt = '('                                                      &
   //'" * With run number                            ",i26," *")'
  c82fmt = '('                                                      &
   //'" * No data read from tape                     ",26x," *")'
  c83fmt = '('                                                      &
   //'" * Data written to tape                       ",26x," *")'
  c84fmt = '('                                                      &
   //'" * With experment name                        ",a26," *")'
  c85fmt = '('                                                      &
   //'" * With run number                            ",i26," *")'
  c86fmt = '('                                                      &
   //'" * No data written to tape                    ",26x," *")'
  c96fmt = '('                                                      &
   //'" * Soil type ice                              ",26x," *")'
  c97fmt = '('                                                      &
   //'" * Soil type fine                             ",26x," *")'
  c98fmt = '('                                                      &
   //'" * Soil type medium                           ",26x," *")'
  c99fmt = '('                                                      &
   //'" * Soil type coarse                           ",26x," *")'
  c100fmt = '('                                                     &
   //'" * Initial soil temperature profile derived   ",26x," *")'
  c101fmt = '('                                                     &
   //'" * From soil temp. cycle with mean annual   ",f28.1," *")'
  c102fmt = '('                                                     &
   //'" * Daily amplitude of soil temperature cycle",f28.1," *")'
  c103fmt = '('                                                     &
   //'" * Annual amplitude of soil temperature cycle",f27.1," *")'

!---------------------------------------------------------------------
!     Output the initial data to all the unit nos. where diagnostics
!     are required
!---------------------------------------------------------------------

  land_cnt = 0

  DO l=1, rows          ! loop on points
    DO k=1, row_length

      IF (land_sea_mask(k,l)) THEN
        land_cnt = land_cnt + 1
      END IF

      IF (row_length*rows  >  1 )                                   &
        WRITE (nout(numnout),*) " COLUMN NUMBER : (",k,l, ")"
      IF (nout(numnout)  /=  0) THEN

        WRITE (nout(numnout),c1fmt)
        WRITE (nout(numnout),c2fmt)
        IF  (land_sea_mask(k,l)) THEN
          WRITE (nout(numnout),c3fmt)
        ELSE
          WRITE (nout(numnout),c4fmt)
        END IF
        WRITE (nout(numnout),c5fmt) year_init
        WRITE (nout(numnout),c6fmt) dayno_init
        WRITE (nout(numnout),c7fmt) time_init
        WRITE (nout(numnout),c8fmt) ndayin, nminin, nsecin
        WRITE (nout(numnout),c10fmt) timestep
        WRITE (nout(numnout),c11fmt1) timestep*a_sw_radstep_prog
        WRITE (nout(numnout),c11fmt2) timestep*a_sw_radstep_diag
        WRITE (nout(numnout),c12fmt) timestep*ntrad1
        WRITE (nout(numnout),c13fmt) lat(k,l)
        WRITE (nout(numnout),c14fmt) long(k,l)
        WRITE (nout(numnout),c15fmt) nlevs
        WRITE (nout(numnout),c16fmt) nwet
        WRITE (nout(numnout),c17fmt) nbl_levs
        IF (land_sea_mask(k,l)) THEN
          WRITE (nout(numnout),c18fmt) nsoilt_levs
        END IF

! If any physics routines are switched off or switched to
! diagnostics only print out here
        IF (conv_mode  ==  2) THEN
          WRITE (nout(numnout),c19fmt)
        ELSE IF (conv_mode  ==  1) THEN
          WRITE (nout(numnout),c20fmt)
        END IF

! ancyc
        IF (ancyc) THEN
          WRITE (nout(numnout),c31fmt)
        ELSE
          WRITE (nout(numnout),c32fmt)
        END IF

! local time
        IF (local_time) THEN
          WRITE (nout(numnout),c33fmt)
        ELSE
          WRITE (nout(numnout),c34fmt)
        END IF

! obs
        IF (obs) THEN
          WRITE (nout(numnout),c35fmt)
          WRITE (nout(numnout),c36fmt) ichgf
          WRITE (nout(numnout),c37fmt) ilscnt
          WRITE (nout(numnout),c38fmt) nfor

          multnf = INT(nfor/10)
          resnf = MOD(nfor,10)

! tls
          WRITE (nout(numnout),c39fmt)
          IF (multnf  /=  0) THEN
            DO j=1, nlevs
              DO m=1, multnf
                WRITE (nout(numnout),'(10f7.2)')                    &
                  (tls(k,l,i,j), i = (m-1)*10+1, 10*m)
              END DO
              WRITE (nout(numnout),'(10f7.2)')                      &
                (tls(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
            END DO
          ELSE
            DO j=1, nlevs
              WRITE (nout(numnout),'(10f7.2)')                      &
                (tls(k,l,i,j), i = 1, resnf)
            END DO
          END IF               ! multnf

! tstar_forcing
          WRITE (nout(numnout),c39fmt)
          IF (multnf  /=  0) THEN
            DO m=1, multnf
              WRITE (nout(numnout),'(10f7.2)')                      &
                (tstar_forcing(k,l,i), i = (m-1)*10+1, 10*m)
            END DO
            WRITE (nout(numnout),'(10f7.2)')                        &
              (tstar_forcing(k,l,i), i = (m-1)*10+1,                &
                                               (m-1)*10+resnf)
          ELSE
            WRITE (nout(numnout),'(10f7.2)')                        &
              (tstar_forcing(k,l,i), i = 1, resnf)
          END IF               ! multnf

! uls
          WRITE (nout(numnout),c40fmt)
          IF (multnf  /=  0) THEN
            DO j=1, nlevs
              DO m=1, multnf
                WRITE (nout(numnout),'(10f7.2)')                    &
                  (uls(k,l,i,j), i = (m-1)*10+1, 10*m)
              END DO
              WRITE (nout(numnout),'(10f7.2)')                      &
                (uls(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
            END DO
          ELSE
            DO j=1, nlevs
              WRITE (nout(numnout),'(10f7.2)')                      &
                (uls(k,l,i,j), i = 1, resnf)
            END DO
          END IF               ! multnf
! vls
          WRITE (nout(numnout),c41fmt)
          IF (multnf  /=  0) THEN
            DO j=1, nlevs
              DO m=1, multnf
                WRITE (nout(numnout),'(10f7.2)')                    &
                  (vls(k,l,i,j),i=(m-1)*10+1,10*m)
              END DO
              WRITE (nout(numnout),'(10f7.2)')                      &
                (vls(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
            END DO
          ELSE
            DO j=1, nlevs
              WRITE (nout(numnout),'(10f7.2)')                      &
                (vls(k,l,i,j),i=1,resnf)
            END DO
          END IF               !multnf

! wls
          WRITE (nout(numnout),c41fmt)
          IF (multnf  /=  0) THEN
            DO j=0, nlevs
              DO m=1, multnf
                WRITE (nout(numnout),'(10f7.2)')                    &
                  (wls(k,l,i,j),i=(m-1)*10+1,10*m)
              END DO
              WRITE (nout(numnout),'(10f7.2)')                      &
                (wls(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
            END DO
          ELSE
            DO j=0, nlevs
              WRITE (nout(numnout),'(10f7.2)')                      &
                (wls(k,l,i,j),i=1,resnf)
            END DO
          END IF               !multnf

! qls
          WRITE (nout(numnout),c42fmt)
          multnf = INT(nfor/5)
          resnf = MOD(nfor,5)
          IF (multnf  /=  0) THEN
            DO j=1, nwet
              DO m=1, multnf
                WRITE (nout(numnout),'(5e10.2)')                    &
                (qls(k,l,i,j),i=(m-1)*5+1,5*m)
              END DO
              WRITE (nout(numnout),'(5e10.2)')                      &
              (qls(k,l,i,j), i = (m-1)*5+1, (m-1)*5+resnf)
            END DO
          ELSE
            DO j=1, nlevs
              WRITE (nout(numnout),'(5e10.2)')                      &
                (qls(k,l,i,j),i=1,resnf)
            END DO
          END IF               ! multnf

! flux_h
          WRITE (nout(numnout),c43fmt)
          multnf = INT(nfor/5)
          resnf = MOD(nfor,5)
          IF (multnf  /=  0) THEN
            DO m=1, multnf
              WRITE (nout(numnout),'(5f10.5)')                      &
                (flux_h(k,l,i),i=(m-1)*5+1,5*m)
            END DO
            WRITE (nout(numnout),'(5f10.5)')                        &
              (flux_h(k,l,i), i=(m-1)*5+1,(m-1)*5+resnf)
          ELSE
            WRITE (nout(numnout),'(5f10.5)')                        &
              (flux_h(k,l,i),i=1,resnf)
          END IF               !multnf

! flux_e
          WRITE (nout(numnout),c44fmt)
          IF (multnf  /=  0) THEN
            DO m=1, multnf
              WRITE (nout(numnout),'(5f10.5)')                      &
                (flux_e(k,l,i),i=(m-1)*5+1,5*m)
            END DO
            WRITE (nout(numnout),'(5f10.5)')                        &
              (flux_e(k,l,i), i = (m-1)*5+1, (m-1)*5+resnf)
          ELSE
            WRITE (nout(numnout),'(5f10.5)')                        &
             (flux_e(k,l,i),i=1,resnf)
          END IF               ! multnf
        END IF                 ! obs

! STATS
        IF (stats) THEN
          WRITE (nout(numnout),c45fmt)
          WRITE (nout(numnout),c46fmt) change_clim
        END IF                 ! stats

! noforce
        IF (noforce) THEN
          WRITE (nout(numnout),c47fmt)
        END IF                 ! noforce

! geoforce
        IF (geoforce) THEN
          WRITE (nout(numnout),c48fmt)
          IF (geoinit) THEN
            WRITE (nout(numnout),c49fmt)
          END IF
        END IF                 ! goforce

! obs, noforce, altdat
        IF (noforce .OR.  obs .OR. altdat .OR. geoforce) THEN
          multnl = INT((nlevs+1)/10)
          resnl = MOD((nlevs+1),10)
          IF (noforce .OR. obs                                      &
           .OR. (geoforce .AND. .NOT. geoinit)) THEN

! wi
            WRITE (nout(numnout),c50fmt)
            IF (multnl /= 0) THEN
              DO j=1,multnl
                WRITE (nout(numnout),'("  ",10f7.2)')               &
                  (wi(k,l,i), i = 10*(j-1)+1, 10*j)
              END DO
              WRITE (nout(numnout),'("  ",10f7.2)')                 &
                (wi(k,l,i), i = 10*(j-1)+1, 10*(j-1)+resnl)
            ELSE
              WRITE (nout(numnout),'("  ",10f7.2)')                 &
                (wi(k,l,i),i=1,resnl)
            END IF             ! multnl
          END IF

          multnl = INT(nlevs/10)
          resnl = MOD(nlevs,10)

          IF (noforce .OR. obs                                      &
           .OR. (geoforce .AND. .NOT. geoinit)) THEN
! ui
            WRITE (nout(numnout),c50fmt)
            IF (multnl /= 0) THEN
              DO j=1,multnl
                WRITE (nout(numnout),'("  ",10f7.2)')               &
                  (ui(k,l,i), i = 10*(j-1)+1, 10*j)
              END DO
              WRITE (nout(numnout),'("  ",10f7.2)')                 &
                (ui(k,l,i), i = 10*(j-1)+1, 10*(j-1)+resnl)
            ELSE
              WRITE (nout(numnout),'("  ",10f7.2)')                 &
                (ui(k,l,i),i=1,resnl)
            END IF             ! multnl

! vi
            WRITE (nout(numnout),c51fmt)
            IF (multnl  /=  0) THEN
              DO j=1, multnl
                WRITE (nout(numnout),'("  ",10f7.2)')               &
                  (vi(k,l,i), i = 10*(j-1)+1, 10*j)
              END DO
              WRITE (nout(numnout),'("  ",10f7.2)')                 &
                (vi(k,l,i), i = 10*(j-1)+1, 10*(j-1)+resnl)
            ELSE
              WRITE (nout(numnout),'("  ",10f7.2)')                 &
                (vi(k,l,i), i = 1, resnl)
            END IF             ! multnl

          ELSE IF (geoforce .AND. geoinit) THEN
! ug, vg
            multnf = INT(nfor/10)
            resnf = MOD(nfor,10)

            IF (ug_opt == 1) THEN
              WRITE (nout(numnout),c52fmt) ug(k,l,1,1)
            ELSE IF (ug_opt == 2) THEN
              IF (multnf  /=  0) THEN
                DO j=1, nlevs
                  DO m=1, multnf
                    WRITE (nout(numnout),'(10f7.2)')                &
                       (ug(k,l,i,j), i = (m-1)*10+1, 10*m)
                  END DO
                  WRITE (nout(numnout),'(10f7.2)')                  &
                     (ug(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
                END DO
              ELSE
                DO j=1, nlevs
                  WRITE (nout(numnout),'(10f7.2)')                  &
                     (ug(k,l,i,j), i = 1, resnf)
                END DO
              END IF               ! multnf
            END IF  ! if ug_opt

            IF (vg_opt == 1) THEN
              WRITE (nout(numnout),c53fmt) vg(k,l,1,1)
            ELSE IF (vg_opt == 2) THEN
              IF (multnf  /=  0) THEN
                DO j=1, nlevs
                  DO m=1, multnf
                    WRITE (nout(numnout),'(10f7.2)')                &
                       (vg(k,l,i,j), i = (m-1)*10+1, 10*m)
                  END DO
                  WRITE (nout(numnout),'(10f7.2)')                  &
                     (vg(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
                END DO
              ELSE
                DO j=1, nlevs
                  WRITE (nout(numnout),'(10f7.2)')                  &
                     (vg(k,l,i,j), i = 1, resnf)
                END DO
              END IF               ! multnf
            END IF  ! if vg_opt
          END IF               ! obs or noforce or geoforce

! ti
          WRITE (nout(numnout),c54fmt)
          IF (multnl  /=  0) THEN
            DO j=1,multnl
              WRITE (nout(numnout),'("  ",10f7.2)')                 &
                (ti(k,l,i),i=10*(j-1)+1,10*j)
            END DO
            WRITE (nout(numnout),'("  ",10f7.2)')                   &
              (ti(k,l,i),i=10*(j-1)+1,10*(j-1)+resnl)
          ELSE
            WRITE (nout(numnout),'("  ",10f7.2)')                   &
              (ti(k,l,i),i=1,resnl)
          END IF               !multnl

! qi
          WRITE (nout(numnout),c55fmt)
          multnl = INT(nwet/5)
          resnl = MOD(nwet,5)
          IF (multnl  /=  0) THEN
            DO j=1, multnl
              WRITE (nout(numnout),'("  ",5e11.5)')                 &
                (qi(k,l,i),i=5*(j-1)+1,5*j)
            END DO
            WRITE (nout(numnout),'("  ",5e11.5)')                   &
              (qi(k,l,i),i=5*(j-1)+1,5*(j-1)+resnl)
          ELSE
            WRITE (nout(numnout),'("  ",5e11.5)')                   &
              (qi(k,l,i),i=1,resnl)
          END IF               ! multnl
        END IF                 ! end if altdat,obs,noforce

! Ozone
        WRITE (nout(numnout),c56fmt)
        multoz = INT(nozone/5)
        resoz = MOD(nozone,5)
        IF (multoz  /=  0) THEN
          DO j=1, multoz
            WRITE (nout(numnout),'("  ",5e11.5)')                   &
              (ozone(k,l,i),i=5*(j-1)+1,5*j)
          END DO
          WRITE (nout(numnout),'("  ",5e11.5)')                     &
            (ozone(k,l,i),i=5*(j-1)+1,5*(j-1)+resoz)
        ELSE
          WRITE (nout(numnout),'("  ",5e11.5)')                     &
            (ozone(k,l,i),i=1,resoz)
        END IF                 ! multoz

        WRITE (nout(numnout),c57fmt) ccai(k,l)
        WRITE (nout(numnout),c58fmt) iccbi(k,l)
        WRITE (nout(numnout),c59fmt) iccti(k,l)
        multnl = INT(nlevs+1/10)
        resnl = MOD(nlevs+1,10)
        DO j=1, multnl
          WRITE (nout(numnout),c60fmt)                              &
            (p_in(k,l,i), i = 10*(j-1)+1, 10*j)
        END DO
        multnl = INT(nlevs/10)
        resnl = MOD(nlevs,10)

! radcloud_fixed
        IF (radcloud_fixed) THEN
          WRITE (nout(numnout),c61fmt) cca_rad(k,l)
          WRITE (nout(numnout),c62fmt) iccb_rad(k,l)
          WRITE (nout(numnout),c63fmt) icct_rad(k,l)
          WRITE (nout(numnout),c64fmt) ccwpin_rad(k,l)

          multcl = INT(nwet/5)
          rescl = MOD(nwet,5)

! layer_cloud_rad
          WRITE (nout(numnout),c65fmt)
          IF (multcl  /=  0) THEN
            DO j=1, multcl
              WRITE (nout(numnout),'("  ",5e11.5)')                 &
                (layer_cloud_rad(k,l,i), i=5*(j-1)+1, 5*j)
            END DO
            WRITE (nout(numnout),'("  ",5e11.5)')                   &
              (layer_cloud_rad(k,l,i), i = 5*(j-1)+1, 5*(j-1)+rescl)
          ELSE
            WRITE (nout(numnout),'("  ",5e11.5)')                   &
              (layer_cloud_rad(k,l,i), i = 1, rescl)
          END IF

! qcl_rad
          WRITE (nout(numnout),c66fmt)
          IF (multcl  /=  0) THEN
            DO j=1, multcl
              WRITE (nout(numnout),'("  ",5e11.5)')                 &
                (qcl_rad(k,l,i), i = 5*(j-1)+1, 5*j)
            END DO
            WRITE (nout(numnout),'("  ",5e11.5)')                   &
              (qcl_rad(k,l,i),i=5*(j-1)+1,5*(j-1)+rescl)
          ELSE
            WRITE (nout(numnout),'("  ",5e11.5)')                   &
              (qcl_rad(k,l,i), i = 1, rescl)
          END IF
        END IF                 ! radcloud fixed

! landmask
        IF (land_sea_mask(k,l)) THEN

          WRITE (nout(numnout),c67fmt) canopy_gbi(land_cnt)
          WRITE (nout(numnout),c68fmt) smci(land_cnt)
          WRITE (nout(numnout),c69fmt) snodepi(k,l)

          multsl = INT((nsoilt_levs)/10)
          ressl  = MOD((nsoilt_levs),10)

          ! t_deep_soil
          WRITE (nout(numnout),c70fmt)
          IF (multsl  /=  0) THEN

            DO j=1, multsl
              WRITE (nout(numnout),'("  ",10f7.2)')                       &
                (t_deep_soili(land_cnt,i),i = 10*(j-1)+1, 10*j)
            END DO

            WRITE (nout(numnout),'("  ",10f7.2)')                         &
              (t_deep_soili(land_cnt,i), i = 10*(j-1)+1,10*(j-1)+ressl)
          ELSE
            WRITE (nout(numnout),'("  ",10f7.2)')                         &
              (t_deep_soili(land_cnt,i), i = 1, ressl)
          END IF

          WRITE (nout(numnout),c71fmt) tstari(k,l)

        ELSE

          WRITE (nout(numnout),c78fmt) z0mseai

        END IF      ! land point

        IF (tapein) THEN
          WRITE (nout(numnout),c79fmt)
          WRITE (nout(numnout),c80fmt) exname_in
          WRITE (nout(numnout),c81fmt) runno_in
        ELSE
          WRITE (nout(numnout),c82fmt)
        END IF

        IF (tapeout) THEN
          WRITE (nout(numnout),c83fmt)
          WRITE (nout(numnout),c84fmt) exname_out
          WRITE (nout(numnout),c85fmt) runno_out
        ELSE
          WRITE (nout(numnout),c86fmt)
        END IF

        IF (land_sea_mask(k,l)) THEN
          IF (soil_type(land_cnt) ==  1) THEN
            WRITE (nout(numnout),c96fmt)
          ELSE IF (soil_type(land_cnt)  ==  2) THEN
            WRITE (nout(numnout),c97fmt)
          ELSE IF (soil_type(land_cnt)  ==  3) THEN
            WRITE (nout(numnout),c98fmt)
          ELSE IF (soil_type(land_cnt)  ==  4) THEN
            WRITE (nout(numnout),c99fmt)
          END IF
        END IF

        WRITE (nout(numnout),c2fmt)
      END IF                   ! nout
    END DO                     ! k
  END DO                      ! l

  IF (lhook) CALL dr_hook('PRINT_INITDATA',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE print_initdata

