! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     Subroutine forcing
!
!     Purpose: Called by scm_main (Single Column Model main routine)
!              to apply the appropriate Forcing (code was previously
!              in the main Calling routine scm_main ).
!
!     Code Description:
!     Language - FORTRAN 90
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!     Documentation: Single Column Model Guide - J. Lean
!=====================================================================
! Options to set initial profiles
!=====================================================================
! (i)   Observational large scale forcing (OBS=TRUE of namelist LOGIC)
!         Initial data is then from namelist INPROF
! (ii)  Statistical large scale forcing (STATS=TRUE of namelist LOGIC)
!         Initial data can either be derived from climate datasets
!         using subroutine INITSTAT or set from namelist
!         INPROF (set ALTDAT=TRUE in namelist LOGIC)
! (iii) No large-scale forcing initial data is set fron namelist
!         INPROF
! (iv)  Continuation from previous run stored on tape
!         (Set TAPEIN=TRUE in namelist LOGIC).  All other initial data
!         is overwritten
!=====================================================================
!---------------------------------------------------------------------

SUBROUTINE forcing                                                            &

  ! (In)
  ( row_length, rows, nlevs, nwet, nfor, nbl_levs, nsoilt_levs, nsoilm_levs   &
  , ntrop, sec_day, stepcount, daycount, dayno_wint, daysteps, nSCMdpkgs      &
  , ntab, ichgf, t, q, qcl, qcf, u, v, w, l_scmdiags, p, exner_theta_levels   &
  , rp, r_theta_levs, ch_tstar_forcing, ch_flux_h, ch_flux_e                  &
  , ch_tls,  ch_qls,  ch_uls,  ch_vls,  ch_wls                                &
  , ch_t_bg, ch_q_bg, ch_u_bg, ch_v_bg, ch_w_bg                               &
  , ad, at, avn, aw, cdbar, ctbar, cvnbar, cwbar, dbar, tbar, vnbar, wbar     &
  , vpbar, cdsd, ctsd, cvnsd, cwsd, dsd, tsd, vnsd, wsd, tdash, ddash         &
  , deltan, px, py                                                            &

  ! (InOut)
  , ilscnt, flux_h_scm, flux_e_scm, ti_scm, qi_scm, t_inc_scm, q_star_scm     &
  , qcl_inc, qcf_inc, u_inc_scm, v_inc_scm, w_inc_scm                         &
  , t_bg_scm, q_bg_scm, u_bg_scm, v_bg_scm, w_bg_scm                          &
  , tls, qls, uls, vls, wls, tr, qr, vnr, vpr, wr                             &

  ! (Out)
  , iv, iy, idum, tstar, factor_rhokh, rhokh, dab1, dap1 )

  USE scm_utils, ONLY:                                                        &
    old_rlx, old_vertadv, old_nml, zhook_in, zhook_out, jprb, lhook, dr_hook  &
  , rlx_none, rlx_init, rlx_bgrd, rlx_inst_init, rlx_inst_bgrd

  USE s_main_force, ONLY:                                                     &
    timestep, ui, vi, wi, obs, obs_surf, stats, prindump_obs, prinstat        &
  , l_vertadv, rlx_t, rlx_q, rlx_u, rlx_v, rlx_w                              &
  , plev_t, plev_q, plev_u, plev_v, plev_w, tau_t, tau_q, tau_u, tau_v, tau_w


  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE

!---------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------

!==========================
! Arguments with Intent(In)
!==========================
  INTEGER, INTENT(In) ::  &
    row_length            &! Leading dimension of SCM
  , rows                  &
  , nlevs                 &! No of levels.
  , nwet                  &! No of model levels in which Q is
  , nfor                  &! Number terms for observational forcing
  , nbl_levs              &! Number of Boundary layer levels
  , nsoilt_levs           &! Number of soil temp levels
  , nsoilm_levs           &! Number of soil moisture levels
  , ntrop                  ! Max number of levels in the troposphere

  INTEGER, INTENT(In) ::  &
    sec_day               &
  , stepcount             &! Timestep counter
  , daycount              &! Daynumber (1 represents 1st january)
  , dayno_wint            &! Daynumber relative to winter solstice
  , daysteps              &! No. of timesteps in 1 day
  , nSCMdpkgs             &! No of SCM diagnostics packages
  , ntab                  &! Dimension of array used in random generator
  , ichgf                  ! No. of timesteps between change
                           ! in observational forcing

  REAL, INTENT(In) ::               &
    t(row_length,rows,nlevs)        &! Temperature (K)
  , q(row_length,rows,nwet)         &! Specific humidity (kg/kg)
  , qcl(row_length,rows,nwet)       &! Cloud water content (kg/kg)
  , qcf(row_length,rows,nwet)       &! Cloud ice content (kg/kg)
  , u(row_length,rows,nlevs)        &! Zonal wind (m/s)
  , v(row_length,rows,nlevs)        &! Meridional wind (m/s)
  , w(row_length,rows,0:nlevs)       ! Vertical velocity (m/s)

  LOGICAL, INTENT(In) ::  &
    l_scmdiags(nscmdpkgs)  ! Logicals for SCM diagnostics packages


  REAL, INTENT(In) ::                         &
    p(row_length,rows,nlevs)                  &! Pressure coordinates (Pa)
  , exner_theta_levels(row_length,rows,nlevs) &
  , rp(row_length,rows,nlevs)                 &! 1/p for rho levels
                                               ! (1/HPa or 1/mb)
  , r_theta_levs(row_length,rows,0:nlevs)      ! Distance of theta levels
                                               ! from centre of earth (m)


  ! Rate of change of surface forcing for ...
  REAL, INTENT(In) ::                         &
    ch_tstar_forcing(row_length,rows,nfor-1)  &! ... tstar,  K/s
  , ch_flux_h(row_length,rows,nfor-1)         &! ... flux_h, W/(m^2*s)
  , ch_flux_e(row_length,rows,nfor-1)          ! ... flux_e, W/(m^2*s)

  ! Rate of change of atmospheric forcing tendencies for ...
  REAL, INTENT(In) ::                       &
    ch_tls(row_length,rows,nfor-1,nlevs)    &! ... Temperature K/(day*s)
  , ch_qls(row_length,rows,nfor-1,nwet)     &! ... Spec. humid (kg/kg)/(day*s)
  , ch_uls(row_length,rows,nfor-1,nlevs)    &! ... Zonal  wind m/(day*s^2)
  , ch_vls(row_length,rows,nfor-1,nlevs)    &! ... Merid. wind m/(day*s^2)
  , ch_wls(row_length,rows,nfor-1,0:nlevs)   ! ... Vert.  wind m/(day*s^2)

  ! Rate of change of background fields for ...
  REAL, INTENT(In) ::                       &
    ch_t_bg(row_length,rows,nfor-1,nlevs)   &! ... Temperature  K/s
  , ch_q_bg(row_length,rows,nfor-1,nwet)    &! ... Spec. humid  (kg/kg)/s
  , ch_u_bg(row_length,rows,nfor-1,nlevs)   &! ... Zonal  wind  m/(s^2)
  , ch_v_bg(row_length,rows,nfor-1,nlevs)   &! ... Merid. wind  m/(s^2)
  , ch_w_bg(row_length,rows,nfor-1,0:nlevs)  ! ... Vert.  wind  m/(s^2)


  ! Declarations for Statistical forcing
  !--------------------------------------
  REAL, INTENT(In) ::               &
    ad(row_length,rows,nwet-1)      &! Term a of eqn. 2.22 for
  , at(row_length,rows,nlevs-1)     &!   dew pt depression and temperature
  , avn(row_length,rows,nlevs-1)    &! Term a of equ. 2.22 for
  , aw(row_length,rows,ntrop-1)      !   horiz. and vertical velocities

  ! Mean of random variable for ...
  REAL, INTENT(In) ::               &
    cdbar(row_length,rows,nwet)     &!   ... Dew Point Depression (K)
  , ctbar(row_length,rows,nlevs)    &!   ... Temperature (K)
  , cvnbar(row_length,rows,nlevs)   &!   ... Velocity, VN (m/s)
  , cwbar(row_length,rows,ntrop)     !   ... Vertical Velocity (mb/s)

  ! Mean at daycount days from winter solstice of ...
  REAL, INTENT(In) ::               &
    dbar(row_length,rows,nwet)      &!   ... Dew Point Depression (K)
  , tbar(row_length,rows,nlevs)     &!   ... Temperature (K)
  , vnbar(row_length, rows,nlevs)   &!   ... Velocity, VN (m/s)
  , wbar(row_length,rows,ntrop)     &!   ... Vertical velocity (mb/s)
  , vpbar(row_length,rows,nlevs)     !   ... Velocity, VP (m/s)

  ! Standard Deviation of random variable for ...
  REAL, INTENT(In) ::               &
    cdsd(row_length,rows,nwet)      &!   ... Dew Point Depression (K)
  , ctsd(row_length,rows,nlevs)     &!   ... Temperature (K)
  , cvnsd(row_length,rows,nlevs)    &!   ... Velocity, VN (m/s)
  , cwsd(row_length,rows,ntrop)      !   ... Vertical velocity (mb/s)

  ! Standard Deviation at daycount days from winter solstice of ...
  REAL, INTENT(In) ::               &
    dsd(row_length,rows,nwet)       &!   ... Dew Point Depression (K)
  , tsd(row_length,rows,nlevs)      &!   ... Temperature (K)
  , vnsd(row_length,rows,nlevs)     &!   ... Velocity, VN (m/s)
  , wsd(row_length,rows,ntrop)       !   ... Vertical Velocity (mb/s)


  REAL, INTENT(In) ::               &
    tdash(row_length,rows,nlevs)    &! Temp. corrections (K)
  , ddash(row_length,rows,nwet)     &! Dew pt. corrections
  , deltan(row_length,rows)         &! Radius of area (m)
  , px(row_length,rows,ntrop)       &! Reciprocal log functions for calc. of
  , py(row_length,rows,ntrop-1)      !   advection used in eqns 2.12 and 2.13


!==============================
! Arguments with intent (InOut)
!==============================
  INTEGER, INTENT(InOut) ::             &
    ilscnt                               ! Count for observational forcing

  REAL, INTENT(InOut) ::                &
    flux_h_scm(row_length,rows)         &! Srf sensible heat flux (W/m2)
  , flux_e_scm(row_length,rows)          ! Srf latent heat flux   (W/m2)

  REAL, INTENT(InOut) ::                &
    ti_scm(row_length,rows,nlevs)       &! Initial temperature  (K)
  , qi_scm(row_length,rows,nwet)         ! Initial spec. humid. (kg/kg)

  ! Forcing increments of ...
  REAL, INTENT(InOut) ::                &
    t_inc_scm(row_length,rows,nlevs)    &! ... Temperature  (K)
  , q_star_scm(row_length,rows,nwet)    &! ... Spec. humid. (kg/kg)
  , qcl_inc(row_length,rows,nwet)       &! ... Cloud liquid (kg/kg)
  , qcf_inc(row_length,rows,nwet)       &! ... Cloud ice    (kg/kg)
  , u_inc_scm(row_length,rows,nlevs)    &! ... Zonal  wind  (m/s)
  , v_inc_scm(row_length,rows,nlevs)    &! ... Merid. wind  (m/s)
  , w_inc_scm(row_length,rows,0:nlevs)   ! ... Vert.  wind  (m/s)

  ! Obs. background field of ...
  REAL, INTENT(InOut) ::                &
    t_bg_scm(row_length,rows,nlevs)     &! ... Temperature  (K)
  , q_bg_scm(row_length,rows,nwet)      &! ... Spec. humid. (kg/kg)
  , u_bg_scm(row_length,rows,nlevs)     &! ... Zonal  wind  (m/s)
  , v_bg_scm(row_length,rows,nlevs)     &! ... Merid. wind  (m/s)
  , w_bg_scm(row_length,rows,0:nlevs)    ! ... Vert.  wind  (m/s)


  ! Forcing tendency due to large-scale horizontal/vertical advection of ...
  REAL, INTENT(InOut) ::                &
    tls(row_length,rows,nlevs)          &! ... Temperature (K)/day
  , qls(row_length,rows,nwet)           &! ... Spec. Humid (kg/kg)/day
  , uls(row_length,rows,nlevs)          &! ... Zonal wind  (m/s)/day
  , vls(row_length,rows,nlevs)          &! ... Merid wind  (m/s)/day
  , wls(row_length,rows,0:nlevs)         ! ... Vert. wind  (m/s)/day


  ! Statistical forcing
  !---------------------
  REAL, INTENT(InOut) ::                &! Randomly sampled ...
    tr(row_length,rows,nlevs,2)         &!   ... Temperature (K)
  , qr(row_length,rows,nwet ,2)         &!   ... Specific humidity (kg/kg)
  , vnr(row_length,rows,nlevs,2)        &!   ... Horizontal velocity (m/s)
  , vpr(row_length,rows,nlevs,2)        &!   ... Horizontal velocity (m/s)
  , wr(row_length,rows,ntrop,2)          !   ... Vertical velocity (mb/s)


!==============================
! Arguments with intent (Out)
!==============================
  INTEGER, INTENT(Out) ::               &
    iv(ntab)                            &! State of number generator saved
  , iy                                  &!    on tape for continuation run
  , idum

  REAL, INTENT(Out) ::                  &
    tstar(row_length,rows)              &! Surface Temperature
  , factor_rhokh(row_length,rows)       &
  , rhokh(row_length,rows,nbl_levs)     &
  , dab1(row_length,rows,44)            &! Observational diagnostics
  , dap1(row_length,rows,36,nlevs)       ! Observational diagnostics


!---------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------
  CHARACTER(LEN=80) ::                     &
    Cmessage                             ! Error message if ErrorStatus > 0

  CHARACTER (len=*), PARAMETER ::       &
    RoutineName = 'forcing'

  INTEGER ::                            &
    i, j, k, l

  INTEGER ::                            &
    ErrorStatus

  REAL ::                               &
    tstpfd                              &
  , alpha                               &
  , a2out_all(nlevs)                    &
  , a2out_wet(nwet)

  REAL ::                               &
    qclls(row_length,rows,nwet)         &
  , qcfls(row_length,rows,nwet)

! Include parameters necessary for calls to SCMoutput...
! Start of include file: s_scmop.h
! Description:
!  Declares and defines some parameters necessary for calling SCMoutput
!
!

! Integers to represent the different time profiles. All must
! be non-negative and less than "only_radsteps".

  INTEGER, PARAMETER :: &
    t_inst        = 1   &! Give the instantaneous value
  , t_avg         = 2   &! Construct the average value
  , t_max         = 3   &! " maximum value
  , t_min         = 4   &! " minimum value
  , t_acc         = 5   &! " accumulated value
  , t_div         = 7   &! " average value divided
                         !   by another diagnostic
  , t_mult        = 8   &! " average value multiplied
                         !   by another diagnostic
  , t_acc_div     = 9   &! " accumulated value divided
                         !   by another diagnostic
  , t_acc_mult    = 10  &! " accumulated value multiplied
                         !   by another diagnostic
  , t_const       = 11  &! The value is constant.
  , only_radsteps = 100  ! When added to one of the above parameters,
                         ! flags that the diagnostic is only available
                         ! on radiation timesteps

! Integers to represent the different domain profiles
  INTEGER, PARAMETER :: &
    d_sl      = 1       &
  , d_soilt   = 2       &
  , d_bl      = 3       &
  , d_wet     = 4       &
  , d_all     = 5       &
  , d_soilm   = 6       &
  , d_tile    = 7       &
  , d_vis     = 9       &
  , d_point   = 13      &
  , d_allxtra = 14      &
  , d_land    = 15      &
  , d_cloud   = 16

! Statement function to encode a stream number into an integer
  INTEGER :: &
    Stream   &
  , strm

  Stream(strm) = 2**(strm-1)

! The default streams for diagnostics to go to will be 1,2,3,4,5 and 6.
! The following should thus be equal to:
!
! Stream(1) [2^0=1] + Stream(2) [2^1=2]  + Stream(3) [2^2=4]
! Stream(4) [2^3=8] + Stream(5) [2^4=16] + Stream(6) [2^5=32]
! Total = 63
!
! where Stream() is the statement function defined above.
! default is 63 (all)

  INTEGER, PARAMETER :: &
    default_streams = 63

! Integers to represent the different diagnostics packages
  INTEGER, PARAMETER :: &
    SCMDiag_gen   = 1   & ! General diagnostics
  , SCMDiag_rad   = 2   & ! Radiation
  , SCMDiag_bl    = 3   & ! Boundary layer
  , SCMDiag_surf  = 4   & ! Surface
  , SCMDiag_land  = 5   & ! Land points only
  , SCMDiag_sea   = 6   & ! Sea points only
  , SCMDiag_lsp   = 7   & ! Large scale precip
  , SCMDiag_conv  = 8   & ! Convection
  , SCMDiag_lscld = 9   & ! Large scale cloud
  , SCMDiag_pc2   = 10  & ! PC2
  , SCMDiag_forc  = 11  & ! Forcing
  , SCMDiag_incs  = 12  & ! Increments
  , SCMDiag_gwd   = 13    ! Gravity Wave Drag

! End of include file: s_scmop.h

! The increments to t, u, v, q, qcl, qcf on input
  REAL ::                               &
    t_inc_in(row_length,rows,nlevs)     &
  , u_inc_in(row_length,rows,nlevs)     &
  , v_inc_in(row_length,rows,nlevs)     &
  , w_inc_in(row_length,rows,0:nlevs)   &
  , q_inc_in(row_length,rows,nlevs)     &
  , qcl_inc_in(row_length,rows,nlevs)   &
  , qcf_inc_in(row_length,rows,nlevs)

! The increments to T, u, v, q due to observational forcing
  REAL ::                               &
    t_inc_obs(row_length,rows,nlevs)    &
  , u_inc_obs(row_length,rows,nlevs)    &
  , v_inc_obs(row_length,rows,nlevs)    &
  , w_inc_obs(row_length,rows,0:nlevs)  &
  , q_inc_obs(row_length,rows,nlevs)

! The increments to T, q, qcl, qcf due to interactive vertical
! advection
  REAL ::                               &
    t_vertadv(row_length,rows,nlevs)    &
  , q_vertadv(row_length,rows,nlevs)    &
  , qcl_vertadv(row_length,rows,nlevs)  &
  , qcf_vertadv(row_length,rows,nlevs)

  REAL ::                               &
    th_p                                &
  , th_m

  REAL ::                               &
    w_vertadv(row_length,rows,0:nlevs)

  REAL ::  &
    factor &! Holding variable for calculations
  , dth    &! Change in potential temperature
  , dz      ! Layer thickness

  ! Dr Hook
  !=============================================================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('FORCING',zhook_in,zhook_handle)

!---------------------------------------------------------------------
! Control variable
!---------------------------------------------------------------------
  tstpfd = timestep / sec_day
  alpha  = timestep / 3600.0

!---------------------------------------------------------------------
! Make copies of the increments as they are on input
!---------------------------------------------------------------------

  t_inc_in   (:,:,:) = t_inc_scm  (:,:,:)
  u_inc_in   (:,:,:) = u_inc_scm  (:,:,:)
  v_inc_in   (:,:,:) = v_inc_scm  (:,:,:)
  w_inc_in   (:,:,:) = w_inc_scm  (:,:,:)
  q_inc_in   (:,:,:) = q_star_scm (:,:,:)
  qcl_inc_in (:,:,:) = qcl_inc    (:,:,:)
  qcf_inc_in (:,:,:) = qcf_inc    (:,:,:)

!---------------------------------------------------------------------
! Set instantaneous profiles and budgets to zero for OBS forcing
!---------------------------------------------------------------------

  IF (obs) THEN
    IF (prindump_obs) THEN
      dap1(:,:,:,:) = 0.0
      dab1(:,:,:)   = 0.0
    END IF ! prindump_obs
  END IF  ! obs

!---------------------------------------------------------------------
! If statistical forcing required:-
!   Set up 2 profiles. 1 for start of day plus 1 for start of
!   following day and linearly interpolate between 2 values for
!   all forcing variables.  Increments to T and Q added and U
!   and V calculated
!---------------------------------------------------------------------


! Initialise qcl and qcf keep
  qclls(:,:,:) = 0.0
  qcfls(:,:,:) = 0.0

  IF (stats) THEN
! _ls incs = Atmos_physics1 ; _inc = forced incs
!  We need to keep hold of the increments calculated from physics1 and
!  then add them on to the ones calculated from the stats forcing

    tls(:,:,:) = t_inc_scm(:,:,:)
    uls(:,:,:) = u_inc_scm(:,:,:)
    vls(:,:,:) = v_inc_scm(:,:,:)
    wls(:,:,:) = w_inc_scm(:,:,:)

    qls  (:,:,:) = q_star_scm(:,:,:)
    qclls(:,:,:) = qcl_inc(:,:,:)
    qcfls(:,:,:) = qcf_inc(:,:,:)

! DEPENDS ON: statstep
    CALL statstep                                                             &
      ! In
      ( row_length, rows, nlevs, nwet, ntrop, deltan, px, py, daysteps        &
      , stepcount, dayno_wint, tr, vnr, vpr, qr, wr, tbar, tsd, tdash, dbar   &
      , dsd, ddash, vnbar, vpbar, vnsd, wbar, wsd, ctbar, ctsd, at, cdbar     &
      , cdsd, ad, cvnbar, cvnsd, avn, cwbar, cwsd, aw, p, rp, u, v, w, t, q   &
      , prinstat, u_inc_scm, v_inc_scm, w_inc_scm, t_inc_scm, q_star_scm      &
      , daycount, timestep, iv, ntab, iy, idum )

! qcl and qcf temporarily forced using q/100.0
    qcl_inc(:,:,:) = q_star_scm(:,:,:)/100.0
    qcf_inc(:,:,:) = q_star_scm(:,:,:)/100.0

!  Need to add on increments calculated in physics1
    t_inc_scm(:,:,:) = t_inc_scm(:,:,:) + tls(:,:,:)
    u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uls(:,:,:)
    v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vls(:,:,:)

    q_star_scm (:,:,:) = q_star_scm (:,:,:) + qls   (:,:,:)
    qcl_inc    (:,:,:) = qcl_inc    (:,:,:) + qclls (:,:,:)
    qcf_inc    (:,:,:) = qcf_inc    (:,:,:) + qcfls (:,:,:)


! W not used in 4.5, added here, but will need someone to use
! stats forcing in anger to modify.
    DO k=1, ntrop
      DO j=1, rows
        DO i=1, row_length
          w_inc_scm(i,j,k) = w_inc_scm(i,j,k) + wls(i,j,k)
        END DO
      END DO
    END DO


  ELSE IF (obs) THEN
! _ls = forced increments ; _inc = Atmos_Physics1 incs


!---------------------------------------------------------------------
!       Select forcing value for time of day
!---------------------------------------------------------------------
!
! NOTE: ilscnt is used to count the number of observational profiles
!       as the run proceeds. Among other things, it is used to select
!       the appropriate change/gradient to the forcing tendencies when
!       observational intervals greater than the model timestep.
!       THIS FUNCTIONALITY IS FLAWED! It is recommended that users
!       have ilscnt set to 0 in their SCM namelists and have obs
!       forcings commence at the same time they wish to begin their
!       SCM run.                                    Wong/Kerr-Munslow
    IF (MOD((daycount-1) * INT(sec_day)                             &
      + (stepcount-1) * INT(timestep)                               &
      , ichgf * INT(timestep))  ==  0) THEN
      ilscnt = ilscnt + 1
    END IF

    IF (ilscnt  ==  0) THEN
      ilscnt = 1
    ELSE IF (ilscnt  >=  nfor) THEN
      WRITE(Cmessage,*) 'time exceeds forcing period'
      ErrorStatus = 1

      CALL ereport(RoutineName,ErrorStatus,Cmessage)
    END IF

    IF (obs_surf) THEN
      ! Update surface forcings
      ! -----------------------
      flux_h_scm(:,:) = flux_h_scm(:,:) + timestep*ch_flux_h(:,:,ilscnt)
      flux_e_scm(:,:) = flux_e_scm(:,:) + timestep*ch_flux_e(:,:,ilscnt)
      tstar(:,:)      = tstar(:,:)      + timestep*ch_tstar_forcing(:,:,ilscnt)

      rhokh(:,:,1)      = flux_h_scm(:,:)
      factor_rhokh(:,:) = flux_e_scm(:,:)
    END IF

    ! Update current large-scale forcings
    ! -----------------------------------
    tls(:,:,:) = tls(:,:,:) + timestep*ch_tls(:,:,ilscnt,:)
    uls(:,:,:) = uls(:,:,:) + timestep*ch_uls(:,:,ilscnt,:)
    vls(:,:,:) = vls(:,:,:) + timestep*ch_vls(:,:,ilscnt,:)
    wls(:,:,:) = wls(:,:,:) + timestep*ch_wls(:,:,ilscnt,:)
    qls(:,:,:) = qls(:,:,:) + timestep*ch_qls(:,:,ilscnt,:)


    ! Update observed background state (T,q,u,v,w)
    !---------------------------------------------
    t_bg_scm (:,:,:)   = t_bg_scm (:,:,:) + timestep*ch_t_bg(:,:,ilscnt,:)
    q_bg_scm (:,:,:)   = q_bg_scm (:,:,:) + timestep*ch_q_bg(:,:,ilscnt,:)
    u_bg_scm (:,:,:)   = u_bg_scm (:,:,:) + timestep*ch_u_bg(:,:,ilscnt,:)
    v_bg_scm (:,:,:)   = v_bg_scm (:,:,:) + timestep*ch_v_bg(:,:,ilscnt,:)
    w_bg_scm (:,:,:)   = w_bg_scm (:,:,:) + timestep*ch_w_bg(:,:,ilscnt,:)


    ! Add increment due to large-scale forcings
    !------------------------------------------
    t_inc_scm  (:,:,:) = t_inc_scm(:,:,:)  + tstpfd*tls(:,:,:)
    q_star_scm (:,:,:) = q_star_scm(:,:,:) + tstpfd*qls(:,:,:)

    IF (.NOT. old_vertadv) THEN
      w_inc_scm(:,:,:) = w_inc_scm(:,:,:)  + tstpfd*wls(:,:,:)
    END IF

    ! DEPENDS ON: relax_forcing
    CALL relax_forcing                                                        &
      ! In
      ( row_length, rows, nlevs, t, ti_scm, t_bg_scm                          &
      , p, plev_t, tau_t, rlx_t                                               &

      ! InOut
      , t_inc_scm )

    ! DEPENDS ON: relax_forcing
    CALL relax_forcing                                                        &
      ! In
      ( row_length, rows, nwet, q, qi_scm, q_bg_scm                           &
      , p, plev_q, tau_q, rlx_q                                               &

      ! InOut
      , q_star_scm )

    IF (old_rlx) THEN
      !------------------------------------------------------------------------
      ! OLD BEHAVIOUR, NOT RECOMMENDED
      ! (would like to get rid of this when possible)
      !------------------------------------------------------------------------
      SELECT CASE (rlx_u)

      CASE (rlx_none)
        u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uls(:,:,:)*tstpfd

      CASE (rlx_init, rlx_inst_init)
        ! Old EUROCS code
        ! This relaxs back to initial conditions.  This is a mess because the
        ! original eurocs logical to relax uv forcings was not affected by the
        ! old l_windrlx logical. The old eurocs relaxation was also done
        ! BEFORE the winds were updated, while the original non-eurocs
        ! wind relaxtion was done AFTER the winds were updated.

        ! DEPENDS ON: relax_forcing
        CALL relax_forcing                                                    &
          ! In
          ( row_length, rows, nlevs, u, ui, u_bg_scm, p, plev_u, tau_u, rlx_u &

          ! InOut
          , u_inc_scm )

        u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uls(:,:,:)*tstpfd

      CASE (rlx_bgrd, rlx_inst_bgrd)
        ! This relaxs back to obs increments.
        ! DEPENDS ON: relax_forcing
        CALL relax_forcing                                                    &
          ! In
          ( row_length, rows, nlevs, u, ui, u_bg_scm, p, plev_u, tau_u, rlx_u &

          ! InOut
          , u_inc_scm )
      END SELECT


      SELECT CASE (rlx_v)

      CASE (rlx_none)
        v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vls(:,:,:)*tstpfd

      CASE (rlx_init, rlx_inst_init)
        ! Old EUROCS code
        ! This relaxs back to initial conditions.  This is a mess because the
        ! original eurocs logical to relax uv forcings was not affected by the
        ! old l_windrlx logical. The old eurocs relaxation was also done
        ! BEFORE the winds were updated, while the original non-eurocs
        ! wind relaxtion was done AFTER the winds were updated.

        ! DEPENDS ON: relax_forcing
        CALL relax_forcing                                                    &
          ! In
          ( row_length, rows, nlevs, v, vi, v_bg_scm, p, plev_v, tau_v, rlx_v &

          ! InOut
          , v_inc_scm )

        v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vls(:,:,:)*tstpfd

      CASE (rlx_bgrd, rlx_inst_bgrd)
        ! This relaxs back to obs increments
        ! DEPENDS ON: relax_forcing
        CALL relax_forcing                                                    &
          ! In
          ( row_length, rows, nlevs, v, vi, v_bg_scm, p, plev_v, tau_v, rlx_v &

          ! InOut
          , v_inc_scm )

      END SELECT


    ELSE ! Use new revised wind relaxation code

      ! Update
      u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uls(:,:,:)*tstpfd
      v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vls(:,:,:)*tstpfd

      ! DEPENDS ON: relax_forcing
      CALL relax_forcing                                                      &
        ! In
        ( row_length, rows, nlevs, u, ui, u_bg_scm, p, plev_u, tau_u, rlx_u   &

        ! InOut
        , u_inc_scm )

      ! DEPENDS ON: relax_forcing
      CALL relax_forcing                                                      &
        ! In
        ( row_length, rows, nlevs, v, vi, v_bg_scm, p, plev_v, tau_v, rlx_v   &

        ! InOut
        , v_inc_scm )

    END IF ! (old_rlx)


    ! DEPENDS ON: relax_forcing
    CALL relax_forcing                                                        &
      ! In
      ( row_length, rows, nlevs, w(:,:,1:nlevs), wi(:,:,1:nlevs)              &
      , w_bg_scm(:,:,1:nlevs), p, plev_w, tau_w, rlx_w                        &

      ! InOut
      , w_inc_scm(:,:,1:nlevs) )

!-----------------------------------------------------------------------
! SCM Forcing OR Increments Diagnostics Package
!-----------------------------------------------------------------------

    IF (l_SCMdiags(SCMdiag_forc)                                              &
         .OR. l_SCMdiags(SCMdiag_incs)) THEN

      ! These are large-scale forcings and INCLUDE any relaxation effects
      t_inc_obs(:,:,:) = t_inc_scm(:,:,:)  - t_inc_in(:,:,:)
      q_inc_obs(:,:,:) = q_star_scm(:,:,:) - q_inc_in(:,:,:)
      u_inc_obs(:,:,:) = u_inc_scm(:,:,:)  - u_inc_in(:,:,:)
      v_inc_obs(:,:,:) = v_inc_scm(:,:,:)  - v_inc_in(:,:,:)
      w_inc_obs(:,:,:) = w_inc_scm(:,:,:)  - w_inc_in(:,:,:)

    END IF


!-----------------------------------------------------------------------
! Interactive vertical advection
!-----------------------------------------------------------------------
! NOTE: This assumes that the specified forcing does NOT contain
!       increments resulting from LS vertical advection.
!       If this is not the case, using this option will result in
!       double counting.
!-----------------------------------------------------------------------


    IF (l_vertadv) THEN

      IF (old_vertadv) THEN
        w_vertadv(:,:,:) = wls(:,:,:)
      ELSE
        ! w_inc_scm at this point has allowed for increments due to
        ! LS-tendency and relaxation. The original code calculates
        ! w*d(x)/dz, though uses w = large scale w tendency. Shouldn't
        ! it be using w + w increment from specified LS forcing routine.

        ! Interactive vertical advection should be using the
        ! current wind profile not the specified LS tendency
        w_vertadv(:,:,:) = w(:,:,:) + w_inc_scm(:,:,:)
      END IF

      t_vertadv(:,:,:)   = 0.0
      q_vertadv(:,:,:)   = 0.0
      qcl_vertadv(:,:,:) = 0.0
      qcf_vertadv(:,:,:) = 0.0

      !-----------------------------------------------------------------
      ! Vertical temperature advection
      !-----------------------------------------------------------------
      ! Specified large scale forcing does not include contribution
      ! from vertical advection.  Approximate contributions due to
      ! vertical advection using upstream approximation.
      ! NOTE: tls holds horizontal advection plus other fixed forcing.

      DO k=1, nlevs
        DO j=1, rows
          DO i=1, row_length

            factor = - timestep*w_vertadv(i,j,k)*exner_theta_levels(i,j,k)

            IF (w_vertadv(i,j,k) < 0.0) THEN

              ! Subsidence
              IF (k < nlevs) THEN
                th_p = t(i,j,k+1) / exner_theta_levels(i,j,k+1)
                th_m = t(i,j,k)   / exner_theta_levels(i,j,k)

                dth  = th_p - th_m
                dz   = r_theta_levs(i,j,k+1) - r_theta_levs(i,j,k)

                t_vertadv(i,j,k)  = factor*dth / dz
              END IF

            ELSE

              ! Ascent
              IF (k > 1) THEN
                th_p = t(i,j,k)   / exner_theta_levels(i,j,k)
                th_m = t(i,j,k-1) / exner_theta_levels(i,j,k-1)

                dth  = th_p - th_m
                dz   = r_theta_levs(i,j,k) - r_theta_levs(i,j,k-1)

                t_vertadv(i,j,k) = factor*dth / dz
              END IF

            END IF ! w_vertadv < 0.0

          END DO
        END DO
      END DO

      !-----------------------------------------------------------------
      ! Vertical moisture advection
      !-----------------------------------------------------------------
      ! Specified large scale forcing does not include contribution
      ! from vertical advection.  Approximate contributions due to
      ! vertical advection using upstream approximation.
      ! NOTE: qls holds horizontal advection plus other fixed forcing.

      DO k=1, nwet
        DO j=1, rows
          DO i=1, row_length

            factor = - timestep*w_vertadv(i,j,k)

            IF (w_vertadv(i,j,k) < 0.0) THEN

              ! Subsidence
              IF (k < nlevs) THEN
                dz = r_theta_levs(i,j,k+1) - r_theta_levs(i,j,k)

                q_vertadv(i,j,k)   = factor*(q(i,j,k+1)   - q(i,j,k))   / dz
                qcl_vertadv(i,j,k) = factor*(qcl(i,j,k+1) - qcl(i,j,k)) / dz
                qcf_vertadv(i,j,k) = factor*(qcf(i,j,k+1) - qcf(i,j,k)) / dz
              END IF

            ELSE

              ! Ascent
              IF (k > 1) THEN
                dz = r_theta_levs(i,j,k) - r_theta_levs(i,j,k-1)

                q_vertadv(i,j,k)   = factor*(q(i,j,k)   - q(i,j,k-1))   / dz
                qcl_vertadv(i,j,k) = factor*(qcl(i,j,k) - qcl(i,j,k-1)) / dz
                qcf_vertadv(i,j,k) = factor*(qcf(i,j,k) - qcf(i,j,k-1)) / dz
              END IF

            END IF ! w_vertadv < 0.0

          END DO
        END DO
      END DO

      t_inc_scm(:,:,:)  = t_inc_scm(:,:,:)  + t_vertadv(:,:,:)
      q_star_scm(:,:,:) = q_star_scm(:,:,:) + q_vertadv(:,:,:)
      qcl_inc(:,:,:)    = qcl_inc(:,:,:)    + qcl_vertadv(:,:,:)
      qcf_inc(:,:,:)    = qcf_inc(:,:,:)    + qcf_vertadv(:,:,:)

    END IF ! vert_adv

    IF (prindump_obs) THEN
      dap1(:,:,10,:) = t_inc_scm(:,:,:) / sec_day
      dap1(:,:,20,:) = q_star_scm(:,:,:) * 1000.0 / sec_day
    END IF

  END IF                     ! stats or obs

!-----------------------------------------------------------------------
!     SCM Forcing OR Increments Diagnostics Packages
!-----------------------------------------------------------------------
  IF (l_SCMdiags(SCMdiag_forc)                                                &
      .OR. l_SCMdiags(SCMdiag_incs)) THEN

    CALL scmoutput(uls, 'ls_u_inc'                                            &
      , 'Large-scale u-wind forcing tendency', '(m/s)/day'                    &
      , t_inst, d_all, default_streams, '', routinename)

    CALL scmoutput(vls, 'ls_v_inc'                                            &
      , 'Large-scale v-wind forcing tendency', '(m/s)/day'                    &
      , t_inst, d_all, default_streams, '', routinename)

    DO k=1, nlevs
      a2out_all(k) = wls(1,1,k)
    END DO

    CALL scmoutput(a2out_all, 'ls_w_inc'                                      &
      , 'Large-scale w-wind forcing tendency', '(m/s)/day'                    &
      , t_inst, d_all, default_streams, '', routinename)

    If (.NOT. old_nml) THEN

      CALL scmoutput(u_bg_scm, 'bg_u'                                         &
        , 'Background u-wind', 'm/s'                                          &
        , t_inst, d_all, default_streams, '', routinename)

      CALL scmoutput(v_bg_scm, 'bg_v'                                         &
        , 'Background v-wind', 'm/s'                                          &
        , t_inst, d_all, default_streams, '', routinename)

      DO k=1, nlevs
        a2out_all(k) = w_bg_scm(1,1,k)
      END DO

      CALL scmoutput(a2out_all, 'bg_w'                                        &
        , 'Background w-wind', 'm/s'                                          &
        , t_inst, d_all, default_streams, '', routinename)

      CALL scmoutput(q_bg_scm, 'bg_q'                                         &
        , 'Background specific humidity', 'kg/kg'                             &
        , t_inst, d_all, default_streams, '', routinename)

      CALL scmoutput(t_bg_scm, 'bg_t'                                         &
        , 'Background temperature', 'K'                                       &
        , t_inst, d_all, default_streams, '', routinename)

    END IF ! old_nml

    DO k=1, nlevs
      a2out_all(k) = t_inc_scm(1,1,k) - t_inc_in(1,1,k)
    END DO

! DEPENDS ON: scmoutput
    CALL scmoutput(a2out_all, 'dt_totforc'                                    &
      , 'Total temperature increment from s_forcng','K'                       &
      , t_avg, d_all, default_streams, '', RoutineName)

    DO k=1, nlevs
      a2out_all(k) = u_inc_scm(1,1,k) - u_inc_in(1,1,k)
    END DO

! DEPENDS ON: scmoutput
    CALL scmoutput(a2out_all, 'du_totforc'                                    &
      , 'Total u increment from s_forcng', 'm/s'                              &
      , t_avg, d_all, default_streams, '', RoutineName)

    DO k=1, nlevs
      a2out_all(k) = v_inc_scm(1,1,k) - v_inc_in(1,1,k)
    END DO

! DEPENDS ON: scmoutput
    CALL scmoutput(a2out_all, 'dv_totforc'                                    &
      , 'Total v increment from s_forcng', 'm/s'                              &
      , t_avg, d_all, default_streams, '', RoutineName)

    DO k=1, nwet
      a2out_wet(k) = q_star_scm(1,1,k) - q_inc_in(1,1,k)
    END DO

! DEPENDS ON: scmoutput
    CALL scmoutput(a2out_wet, 'dq_totforc'                                    &
      , 'Total humidity increment from s_forcng', 'kg/kg'                     &
      , t_avg, d_wet, default_streams, '', RoutineName)

    DO k=1, nwet
      a2out_wet(k) = qcl_inc(1,1,k) - qcl_inc_in(1,1,k)
    END DO

! DEPENDS ON: scmoutput
    CALL scmoutput(a2out_wet, 'dqcl_totforc'                                  &
      , 'Total QCL increment from s_forcng', 'kg/kg'                          &
      , t_avg, d_wet, default_streams, '', RoutineName)

    DO k=1, nwet
      a2out_wet(k) = qcf_inc(1,1,k) - qcf_inc_in(1,1,k)
    END DO

! DEPENDS ON: scmoutput
    CALL scmoutput(a2out_wet, 'dqcf_totforc'                                  &
      , 'Total QCF increment from s_forcng', 'kg/kg'                          &
      , t_avg, d_wet, default_streams, '', RoutineName)

!-----------------------------------------------------------------------
!       SCM Forcing OR Increments Diagnostics Packages
!       when observational based forcing specified
!-----------------------------------------------------------------------
    IF (obs) THEN

! DEPENDS ON: scmoutput
      CALL scmoutput(t_inc_obs, 'dt_obsforc'                                  &
        , 'Temperature increment from observational forcing', 'K'             &
        , t_avg, d_all, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
      CALL scmoutput(u_inc_obs, 'du_obsforc'                                  &
        , 'U increment from observational forcing', 'm/s'                     &
        , t_avg, d_all, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
      CALL scmoutput(v_inc_obs, 'dv_obsforc'                                  &
        , 'V increment from observational forcing', 'm/s'                     &
        , t_avg, d_all, default_streams, '', RoutineName)

      DO k=1, nlevs
        a2out_all(k) = w_inc_obs(1,1,k)
      END DO

! DEPENDS ON: scmoutput
      CALL scmoutput(a2out_all, 'dw_obsforc'                                  &
        , 'W increment from observational forcing', 'm/s'                     &
        , t_avg, d_all, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
      CALL scmoutput(q_inc_obs, 'dq_obsforc'                                  &
        , 'Humidity increment from observational forcing', 'kg/kg'            &
        , t_avg, d_wet, default_streams, '', RoutineName)

    END IF ! obs

!-----------------------------------------------------------------------
!       SCM Forcing OR Increments Diagnostics Packages
!       when vertical advection specified
!-----------------------------------------------------------------------
    IF (l_vertadv) THEN

! DEPENDS ON: scmoutput
      CALL scmoutput(t_vertadv, 'dt_vertadv'                                  &
        , 'Temperature increment from vertical advection', 'K'                &
        , t_avg, d_all, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
      CALL scmoutput(q_vertadv, 'dq_vertadv'                                  &
        , 'Humidity increment from vertical advection', 'kg/kg'               &
        , t_avg, d_wet, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
      CALL scmoutput(qcl_vertadv, 'dqcl_vertadv'                              &
        , 'QCL increment from vertical advection', 'kg/kg'                    &
        , t_avg, d_wet, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
      CALL scmoutput(qcf_vertadv, 'dqcf_vertadv'                              &
        , 'QCF increment from vertical advection', 'kg/kg'                    &
        , t_avg, d_wet, default_streams, '', RoutineName)

    END IF ! l_vertadv

  END IF ! l_SCMdiags(SCMdiag_forc).OR l_SCMdiags(SCMdiag_incs)

  IF (lhook) CALL dr_hook('FORCING',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE forcing
