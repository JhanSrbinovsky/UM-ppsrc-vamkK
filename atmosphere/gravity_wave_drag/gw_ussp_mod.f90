! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gravity Wave (Ultra-Simple Spectral Parametrization) Scheme.
! Subroutine Interface:

MODULE GW_USSP_MOD

IMPLICIT NONE
  
CONTAINS

      SUBROUTINE GW_USSP(LEVELS, MODEL_DOMAIN, ROWS, NROWS,             &
     &  OFF_X, OFF_Y, HALO_I, HALO_J, ROW_LENGTH,                       &
     &  global_row_length,n_proc, n_procy, proc_row_group,at_extremity, &
     &  R_RHO_LEVELS, R_THETA_LEVELS, P_LAYER_BOUNDARIES,               &
     &  R_U, R_V, T_inc,                                                &
     &  SIN_THETA_LONGITUDE, SIN_THETA_LATITUDE,                        &
     &  THETA, RHO, TIMESTEP, U, V,                                     &
     &  L_ussp_heating,                                                 &
     &  GWSPEC_EFLUX,GWSPEC_SFLUX,GWSPEC_WFLUX,GWSPEC_NFLUX,            &
     &  GWSPEC_EWACC,GWSPEC_NSACC,                                      &
     &  GWSPEC_EFLUX_ON, GWSPEC_EFLUX_P_ON, GWSPEC_SFLUX_ON,            &
     &  GWSPEC_WFLUX_ON, GWSPEC_WFLUX_P_ON, GWSPEC_NFLUX_ON,            &
     &  GWSPEC_EWACC_ON, GWSPEC_EWACC_P_ON, GWSPEC_NSACC_ON)
!
! purpose:    This subroutine calculates the vertical momentum flux
!             divergence due to gravity waves as parametrised by the
!             Warner and McIntyre Ultra Simple Spectral gravity wave
!             Parametrization adapted for use in the UM.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Gravity Wave Drag
!
! code description:
!   language: fortran 90
!   this code is written to umdp3 programming standards.
!   documentation: Unified Model Documentation Paper 34 (Non-Orog GW)

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE g_wave_input_mod, ONLY: L_USSP_OPAQUE, USSP_LAUNCH_FACTOR
      USE earth_constants_mod, ONLY: g, earth_radius, two_omega
      USE atmos_constants_mod, ONLY: R, cp ! Gas constant and
                                           ! heat capacity for dry air
      
! Definitions of prognostic variable array sizes
      USE atm_fields_bounds_mod, ONLY:                                  &
         udims, vdims, pdims, tdims, pdims,                             &
         udims_s, vdims_s, pdims_s, tdims_s, pdims_l, tdims_l                    
     

! Model level values
      USE level_heights_mod, ONLY:                                      & 
         eta_theta_levels           ! Eta values of theta levels

      USE eg_v_at_poles_mod, ONLY: eg_v_at_poles  
          
      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
!$    USE omp_lib
      USE domain_params
      USE p_to_t_mod, ONLY: p_to_t
      USE p_to_u_mod, ONLY: p_to_u
      USE p_to_v_mod, ONLY: p_to_v
      USE polar_row_mean_mod, ONLY: polar_row_mean
      USE u_to_p_mod, ONLY: u_to_p
      USE v_to_p_mod, ONLY: v_to_p
      IMPLICIT NONE

! ----------------------------------------------------------------------+-------
!     Subroutine arguments of GW_USSP
! ----------------------------------------------------------------------+-------
!     Fixed starting rho-level (so Work_PsmallHALO exactly replicated?)
      Integer, Parameter :: pkfix0start    = 0
!     Fixed starting theta-level (e.g. for P_layer_boundaries)
      Integer, Parameter :: tkfix0start    = 0
!
!     Fixed starting theta-level (e.g. for STASH diagnostic arrays)
      Integer, Parameter :: tkfix1start    = 1
!
      INTEGER                                                           &
     &  LEVELS                                                          &
                             !IN Number of model levels
     &, ROWS                                                            &
                             !IN Number of rows for u field
     &, NROWS                                                           &
                             !IN Number of rows for v field
     &, OFF_X                                                           &
                             !IN offset longitude
     &, OFF_Y                                                           &
                             !IN offset latitude
     &, HALO_I                                                          &
                             !IN Halo in longitude
     &, HALO_J                                                          &
                             !IN Halo in latitude
     &, global_row_length                                               &
                             ! number of points on a row
     &, proc_row_group                                                  &
                             ! Group id for processors on the same row
     &, n_proc                                                          &
                             ! Total number of processors
     &, n_procy                                                         &
                             ! Number of processors in latitude                             
     &, ROW_LENGTH                                                      &
                             !IN Number of grid points in row
     &, MODEL_DOMAIN         !IN Model type (global, LAM etc)
      REAL                                                              &
     &  TIMESTEP           & !IN Timestep
     &, SIN_THETA_LONGITUDE(tdims%i_start:tdims%i_end,                  &
     &                      tdims%j_start:tdims%j_end)                  &
                                              !IN Grid point longitudes
     &, SIN_THETA_LATITUDE(tdims%i_start:tdims%i_end,                   &
     &                     tdims%j_start:tdims%j_end)                   &
                                              !IN P-GRID Latitudes
     &, R_THETA_LEVELS(tdims_l%i_start:tdims_l%i_end,                   &
     &                 tdims_l%j_start:tdims_l%j_end,                   &
     &                 0:tdims_l%k_end)                                 &
     &, R_RHO_LEVELS(pdims_l%i_start:pdims_l%i_end,                     &
     &               pdims_l%j_start:pdims_l%j_end,                     &
     &               pdims_l%k_end)                                     &
     &, GWSPEC_EFLUX(udims%i_start:udims%i_end,                         &
     &               udims%j_start:udims%j_end,                         &
     &                 tkfix1start:tdims%k_end)  & ! OUT Fp in E azimuth
     &, GWSPEC_SFLUX(vdims%i_start:vdims%i_end,                         &
     &               vdims%j_start:vdims%j_end,                         &
     &                 tkfix1start:tdims%k_end)  & ! OUT Fp in S azimuth
     &, GWSPEC_WFLUX(udims%i_start:udims%i_end,                         &
     &               udims%j_start:udims%j_end,                         &
     &                 tkfix1start:tdims%k_end)  & ! OUT Fp in W azimuth
     &, GWSPEC_NFLUX(vdims%i_start:vdims%i_end,                         &
     &               vdims%j_start:vdims%j_end,                         &
     &                 tkfix1start:tdims%k_end)  & ! OUT Fp in N azimuth
     &, GWSPEC_EWACC(udims%i_start:udims%i_end,                         &
     &               udims%j_start:udims%j_end,                         &
     &               udims%k_start:udims%k_end)  & ! OUT Accel of U wind
     &, GWSPEC_NSACC(vdims%i_start:vdims%i_end,                         &
     &               vdims%j_start:vdims%j_end,                         &
     &               vdims%k_start:vdims%k_end)  & ! OUT Accel of V wind
     &, THETA(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
     &        tkfix1start:tdims%k_end)                                  &
                             !IN    Primary model array for theta
     &, RHO(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,&
     &      pdims_s%k_start:pdims_s%k_end)    & ! IN Primary model array
                             ! for Density x (radius earth)^2.
     &, U(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,  &
     &    udims_s%k_start:udims_s%k_end)                                &
                             !INOUT Primary model array for U field
     &, V(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,  &
     &    vdims_s%k_start:vdims_s%k_end)                                &
                             !INOUT Primary model array for V field
     &, T_inc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
     &        tkfix1start:tdims%k_end)                                  &
                             !INOUT Temperature increment
     &, P_LAYER_BOUNDARIES(tdims%i_start:tdims%i_end,                   &
     &                     tdims%j_start:tdims%j_end,                   &
     &                     tkfix0start:tdims%k_end)                     &
     &, R_U(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,&
     &      udims_s%k_start:udims_s%k_end)                              &
                             !INOUT U wind increment diagnostic
     &, R_V(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,&
     &      vdims_s%k_start:vdims_s%k_end)
                             !INOUT V wind increment diagnostic
      LOGICAL                                                           &
     &  L_ussp_heating    & !IN Switch to calculate heating tendency
     &, GWSPEC_EFLUX_ON   & !IN 
     &, GWSPEC_EFLUX_P_ON & !IN Switches for diagnostics of
     &, GWSPEC_SFLUX_ON   & !IN Fp in each of 4 azimuths
     &, GWSPEC_WFLUX_ON   & !IN 
     &, GWSPEC_WFLUX_P_ON & !IN 
     &, GWSPEC_NFLUX_ON   & !IN 
     &, GWSPEC_EWACC_ON   & !IN  Switches for acceleration diagnostics
     &, GWSPEC_EWACC_P_ON & !IN 
     &, GWSPEC_NSACC_ON   & !IN 
     &, AT_EXTREMITY(4)     !IN Edge of domain indicator
!
! ----------------------------------------------------------------------+-------
!     Local parameters
! ----------------------------------------------------------------------+-------
!
!     Max number of directions, typically four.
      Integer, Parameter :: IDIR           = 4
!
!     N processor address
      Integer, Parameter :: PNORTH         = 1
!
!     E processor address
      Integer, Parameter :: PEAST          = 2
!
!     S processor address
      Integer, Parameter :: PSOUTH         = 3
!
!     W processor address
      Integer, Parameter :: PWEST          = 4
!
!     Maximum number of iterations of Newton Raphson DO (While) loop
      Integer, Parameter :: MAXWHILE       = 9

!
      Real, Parameter :: R_EARTH_RADIUS_SQ =                            &
     &                                  1. / (Earth_Radius*Earth_Radius)

!
!     Reciprocal of middle atmosphere mean scale height for pressure,
!     normally assumed to be around 7km
      Real, Parameter :: RSCALE_H          = G / (R * 239.145)
!
!     Eta (hybrid height) level for model launch
      Real, Parameter :: ETALAUNCH         = 0.045
!
!     Parameter beta in the launch spectrum total energy equation
      Real, Parameter :: BETA_E0           = 1.0227987125E-1
!
!     Azimuthal sector for launch spectrum integral Delta Phi / 2
      Real, Parameter :: DDPHIR2           = PI / IDIR
!
!     Parameter p in B_0(p) for launch spectrum intrinsic frequency
!     NOTE: This parameter determines the intrinsic frequency spectrum
!           shape and hence the integral form in 4.1, which is strictly
!           valid only for p > 1. !!IF contemplating changes BE WARNED!!
      Real, Parameter :: PSAT              = 5.0 / 3.0
!
!     Psat - 1
      Real, Parameter :: PSATM1            = PSAT - 1.0
!
!     2 - Psat
      Real, Parameter :: TWOMPSAT          = 2.0 - PSAT
!
!     Reciprocal of max wavelength at launch (/m)
      Real, Parameter :: LMINL             = 1./20000
!
!     Reciprocal of characteristic (spectrum peak) wavelength (/m)
      Real, Parameter :: LSTAR             = 1./4300
!
!     Wavenumber at peak in spectrum
      Real, Parameter ::  MSTAR            = 2.0 * PI * LSTAR
!
!     Reciprocal of mstar (m) and mstar^2 (m^2)
      Real, Parameter ::  RMSTAR           = 1. / MSTAR
      Real, Parameter ::  RMSTARSQ         = RMSTAR * RMSTAR
!
!     Minimum vertical wavenumber at launch (/m )
      Real, Parameter ::  MMINL            = 2.0 * PI * LMINL
!
!     Normalised minimum vertical wavenumber at launch
      Real, Parameter ::  MNLMIN           = LMINL / LSTAR
!
!     Equatorial planetary vorticity gradient parameter B_eq (/m /s )
      Real, Parameter :: BETA_EQ_RMSTAR    = 2.3E-11 * RMSTAR
!
!     Power s of vertical wavenumber spectrum A_0(s,t) at low m
      Real, Parameter :: SS                = 1.0
!
!     s + 1, s - 1
      Real, Parameter :: SSP1              = SS + 1.0
!
!     Power t=t_sat of vertical wavenumber spectrum at large m due to
!     saturation by gravity wave breaking (and shape of chopping fn)
      Real, Parameter :: TT                = 3.0
!
!     t - 1, t - 2, 1 / (t-2), (t-3) / (t-2), 2 - t
      Real, Parameter :: TTM1              = TT - 1.0
      Real, Parameter :: TTM2              = TT - 2.0
      Real, Parameter :: RTTM2             = 1.0 / TTM2
      Real, Parameter :: TTRAT             = (TT - 3.0) * RTTM2
      Real, Parameter :: TWOMTT            = 2.0 - TT
!
!     s + t, 1 / (s+t)
      Real, Parameter :: SSPTT             = SS + TT
      Real, Parameter :: RSSPTT            = 1.0 / SSPTT
!
!     Weight for (n+1)th guess at mNlX in iteration solution
      Real, Parameter :: MWEIGHT           = 0.8
!
!     Strength coefficient constant for Launch spectrum (CCL / A0)
      Real, Parameter :: CCL0 = 3.41910625e-9
!
! ----------------------------------------------------------------------+-------
!     Security parameters
! ----------------------------------------------------------------------+-------
!
!     Minimum allowed value of buoyancy frequency squared
      Real, Parameter ::  SQNMIN           = 1.0E-4
!
!     Minimum allowed non-zero value of Curvature Coefficient A
      Real, Parameter ::  ASECP            =  1.0E-20
      Real, Parameter ::  ASECN            = -(ASECP)
!
! ----------------------------------------------------------------------+-------
!     Local Constants (a) Physical
! ----------------------------------------------------------------------+-------
      REAL                                                              &
     &  A0_R_SP1TM1       & ! A0/(s+1)(t-1) normalisation factor for the
!                                    launch spectrum vertical wavenumber
     &, A0_R_1MT          & !-A0/(t-1) normalisation factor for the
!                                    launch spectrum vertical wavenumber
     &, GLOB_LAUNCH_FLUX  & ! mStar**(-2) * [CCL0 * USSP_LAUNCH_FACTOR]
     &, TAIL_CHOP2B       & ! Integral segment [(s+1) * mNLmin**(1-t)]
     &, HEAD_CHOP2A       & ! Integral segment [(t-1) * mNLmin**(s+1)]
     &, COSPHI(4)         & ! Cos(phi_j)
     &, SINPHI(4)           ! Sin(phi_j)
!
! ----------------------------------------------------------------------+-------
!     Local variables (scalars) used in GW_USSP
!     Some effectively expanded to workspace (using vector registers)
! ----------------------------------------------------------------------+-------
      INTEGER                                                           &
     &  I                                                               &
                            !longitude
     &, J                                                               &
                            !latitude
     &, K                                                               &
                            !level
     &, KBASE             & !Base level (1 + minimum theta level)
     &, ILAUNCH                                                         &
                            !Minimum level for launch level (all points)
     &, LAUNCHLEV                                                       &
                            !Launch level at specific point
     &, NNI, NNJ, NNJD                                                  &
                            !Index values
     &, NCHOP2                                                          &
                            !Number of spectra with low m intersect
!    &, NSPIRA              !Number spiral A solution points
     &, JDIR                                                            &
                            !direction
     &, JLEV                                                            &
                            !level
     &, JWHILE                                                          &
                            !Counter for while loop
     &, JJ                                                              &
                            ! omp block iterator
     &, OMP_BLOCK 
                            ! blocking size for omp_block
      REAL                                                              &
     &  F_F                                                             &
                            !Inertial frequency at current latitude
                            !(rad s^-1)
     &, MNLY                                                            &
                            ! High wavenumber intersection point
     &, MKILL                                                           &
                            ! Wavenumber where Doppler transformed
!                             spectrum is reduced to zero
     &, OMINRNBV                                                        &
                            ! omega_min(launch) / N (k)
     &, CCS0RMSTARSQ                                                    &
                            ! Constant component of saturation spectrum
     &, FMINUS                                                          &
                            ! Minimum range value of function f
     &, FTERM                                                           &
                            ! Intermediate  value of function f
     &, FPLUS                                                           &
                            ! Maximum range value of function f
     &, GMINUS                                                          &
                            ! Minimum range value of function g
     &, GTERM                                                           &
                            ! Intermediate  value of function g
     &, GPLUS               ! Maximum range value of function g

      REAL :: G_G           ! Wave-induced force per
                            ! unit mass due to azimuthal sectors (m s^-2)

      REAL :: dzb,dzu,dzl   ! Layer depths

      REAL :: uhat,vhat     ! Velocities on theta levels

      REAL :: ududt,vdvdt   ! Kinetic energy tendencies on theta levels

      LOGICAL                                                           &
     &  L_CHOP2             ! Indicates spectra with low m intersect
! ----------------------------------------------------------------------+-------
!     Local variables (dynamic arrays) used in GW_USSP
! ----------------------------------------------------------------------+-------
!     INTEGER
!    &  NVIEW(LEVELS)       ! check output array
!    &, NVIEW2(LEVELS)      ! check output array
      REAL                                                              &
     &  DDU_a(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
     &        tkfix1start:tdims%k_end,IDIR)                             &
!                            Delta U=udotk(launch)-udotk(jlev)
     &, FPFAC(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
     &        tkfix1start:tdims%k_end,IDIR) & ! Record of total flux
     &, FPTOT(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
     &        tkfix1start:tdims%k_end,IDIR) & ! Pseudomomentum flux
!                         integrated over azimuthal sector (kg m^-1 s^-2)
     &, G_X(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
     &      pdims%k_start:pdims%k_end)      & ! Zonal component of
                           !wave-induced force on density level (m s^-2)
                           !Note that in our notation G_X equates to
                           ! DU_DT, the zonal wind tendency (m s^-2)
     &, G_Y(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
     &      pdims%k_start:pdims%k_end)      & ! Meridional component of 
                           !wave-induced force on density level (m s^-2)
                           !Note that in our notation G_Y equates to
                           ! DV_DT, meridional wind tendency (m s^-2)
     &, ACOEFF(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
     &         tkfix1start:tdims%k_end,IDIR)                            &
!                            Coefficient A in intersect point equation
     &, CURVATURE(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
     &            tkfix1start:tdims%k_end,IDIR)                         &
!                            Term (A / B) in intersect point equation
     &, ATTENUATION(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
     &              tkfix1start:tdims%k_end,IDIR)                       &
!                            Coefficient B in intersect point equation
     &, INTERCEPT1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
     &             tkfix1start:tdims%k_end,IDIR)                        &
!                            Chop function B*[1 + (A/B)]^(t-2)
     &, MINTERCEPT(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
     &             tkfix1start:tdims%k_end,IDIR)                        &
!                            Chop function B*[1 + (A/B)*mNlmin]^(t-2)
     &, MGUESS_a(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
     &           tkfix1start:tdims%k_end+1,IDIR)                        &
!                            Starting value of vertical wavenumber for
!                            crossing point search
     &, MNLX(IDIR*ROWS*ROW_LENGTH,0:MAXWHILE)                           &
                                              ! Intersect mNlX estimates
     &, INDEXI(IDIR*ROWS*ROW_LENGTH)                                    &
                                      ! I location of chop type points
     &, INDEXJ(IDIR*ROWS*ROW_LENGTH)                                    &
                                      ! J location of chop type points
     &, INDEXJD(IDIR*ROWS*ROW_LENGTH)                                   &
                                      ! JDIR location of chop type pnts
     &, ATTE_C(IDIR*ROWS*ROW_LENGTH)                                    &
                                      ! Compressed attenuation array
     &, CURV_C(IDIR*ROWS*ROW_LENGTH)                                    &
                                      ! Compressed curvature array
     &, WGTN(IDIR*ROWS*ROW_LENGTH)                                      &
                                      ! Weighting of n term in iter
     &, OMIN(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                               ! Either f_f or the equatorial minimum
!                            frequency, whichever is less  (rad s^-1)
     &, NBV(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
     &      tkfix1start:tdims%k_end)     & ! Buoyancy [Brunt Vaisala]
!                              frequency on half-levels (rad s^-1)
     &, RHOCL(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,IDIR) &
                                           ! [Rho . Cl]_klaunch
     &, FSATK_SCALE(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
     &         tkfix1start:tdims%k_end)  & ! [Rho(z) . Csat(z) / m*^2]_k
     &, RHO_TH(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
     &         tkfix1start:tdims%k_end)  & ! Rho on theta levels
     &, UDOTK(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
     &        tkfix1start:tdims%k_end,IDIR)  & ! Component of wind in
!                                            phi_jdir direction (m s^-1)
     &, UONP(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
     &       pdims%k_start:pdims%k_end)    & ! LOCAL U on Rho grid
     &, VONP(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
     &       pdims%k_start:pdims%k_end)    & ! LOCAL V on Rho grid
     &, UONP_SMALLHALO(pdims_s%i_start:pdims_s%i_end,                   &
     &            pdims_s%j_start:pdims_s%j_end,                        &
     &            pdims_s%k_start:pdims_s%k_end) & ! LOCAL U with halo
     &, VONP_SMALLHALO(pdims_s%i_start:pdims_s%i_end,                   &
     &            pdims_s%j_start:pdims_s%j_end,                        &
     &            pdims_s%k_start:pdims_s%k_end) & ! LOCAL V with halo
     &, RHONT_SMALLHALO(tdims_s%i_start:tdims_s%i_end,                  &
     &            tdims_s%j_start:tdims_s%j_end,                        &
     &            tkfix1start:tdims_s%k_end)     & !Rho*(radius_earth)^2 
                                                   !on Theta grid.
     &, UINC(udims%i_start:udims%i_end,udims%j_start:udims%j_end,       &
     &       udims%k_start:udims%k_end)  &  ! U Increments on rho grid.
     &, VINC(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,       &
     &       vdims%k_start:vdims%k_end)  &  ! V Increments on rho grid.
     &, Work_PsmallHALO(pdims_s%i_start:pdims_s%i_end,                  &
     &                  pdims_s%j_start:pdims_s%j_end,                  &
     &                  pkfix0start:pdims_s%k_end)                      &
     &, Work_TsmallHALO(tdims_s%i_start:tdims_s%i_end,                  &
     &                  tdims_s%j_start:tdims_s%j_end,                  &
     &                  tkfix0start:tdims_s%k_end)
      LOGICAL                                                           &
     &  L_FTHENG(IDIR*ROWS*ROW_LENGTH) ! Indicate dir of spiral solution

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!- End of Header
!
! ==Main Block==--------------------------------------------------------+-------
! --------------------------------
!     Local Constants (a) Physical
! --------------------------------
      DATA COSPHI/0,-1,0,1/
      DATA SINPHI/1,0,-1,0/
!
      IF (lhook) CALL dr_hook('GW_USSP',zhook_in,zhook_handle)

      GLOB_LAUNCH_FLUX = RMSTARSQ * CCL0 * USSP_LAUNCH_FACTOR
      FMINUS = MNLMIN**SSPTT
      TAIL_CHOP2B = SSP1 / (MNLMIN**TTM1)
      HEAD_CHOP2A = TTM1 * (MNLMIN**SSP1)
      A0_R_SP1TM1 = 1.0 / ( SS + TT - HEAD_CHOP2A )
      A0_R_1MT    = (-(SSP1) ) * A0_R_SP1TM1
      CCS0RMSTARSQ = RMSTARSQ * BETA_E0 * SIN(DDPHIR2) * PSATM1 /       &
     &              (PI * TWOMPSAT)
      KBASE  = tkfix1start + 1
!
! ----------------------------------------------------------------------+-------
! Find model level for launch
! ----------------------------------------------------------------------+-------
!     Minimum value if variable height sources are used
!     WARNING: Not defined relative to tdims%k_start as vatpoles alters it
      ILAUNCH=tkfix1start
      DO K=2,tdims%k_end
        IF (ETA_THETA_LEVELS(K) >  ETALAUNCH.AND.                       &
     &      ETA_THETA_LEVELS(K-1) <  ETALAUNCH) THEN
          IF ((ETA_THETA_LEVELS(K)-ETALAUNCH) <                         &
     &        (ETALAUNCH-ETA_THETA_LEVELS(K-1))) THEN
            ILAUNCH=K
          ELSE
            ILAUNCH=K-1
          ENDIF
        ENDIF
      ENDDO
      LAUNCHLEV = ILAUNCH
!
! ----------------------------------------------------------------------+-------
!     Interpolate : [Vertical]   RHO onto T grid 
!             and : [Horizontal] U,V onto P grid
! ----------------------------------------------------------------------+-------
      CALL P_TO_T (ROW_LENGTH, ROWS, HALO_I, HALO_J,                    &
     &             OFF_X,OFF_Y,LEVELS-1,R_THETA_LEVELS,                 &
     &             R_RHO_LEVELS,RHO,RHONT_SMALLHALO)
!     Extrapolate topmost level to be a scale height from level below
      Rows_do_init: DO J=tdims%j_start,tdims%j_end
        Row_length_do_init: DO I=tdims%i_start,tdims%i_end
          RHONT_SMALLHALO(I,J,tdims%k_end)=                             &
     &      RHONT_SMALLHALO(I,J,tdims%k_end-1) *                        &
     &      EXP(-(R_THETA_LEVELS(I,J,tdims%k_end)-                      &
     &            R_THETA_LEVELS(I,J,tdims%k_end-1)) * RSCALE_H)
        END DO  Row_length_do_init
      END DO  Rows_do_init
!
      CALL u_to_p(u,                                                    &
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        levels,                                         &
                        model_domain,at_extremity,uonp)

      CALL v_to_p(v,                                                    &
                        vdims_s%i_start,vdims_s%i_end,                  & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        levels,                                         &
                        model_domain,at_extremity,vonp)
IF (.NOT. l_vatpoles) THEN 
! ---------------------------
!     Set polar winds to zero
! ---------------------------
      If(model_domain  ==  mt_global) Then
        If (at_extremity(psouth) ) Then
          VONP(:,pdims%j_start,:) = 0.0
          UONP(:,pdims%j_start,:) = 0.0
        Endif
        If (at_extremity(pnorth) ) Then
          VONP(:,pdims%j_end,:) = 0.0
          UONP(:,pdims%j_end,:) = 0.0
        Endif
      Endif
END IF ! vatpoles
!
! ----------------------------------------------------------------------+-------
!     Initialize local arrays to zero
! ----------------------------------------------------------------------+-------
! ---------------------------------
!     Set winds with haloes to zero
! ---------------------------------


      Levels_do1: DO K=pdims_s%k_start,pdims_s%k_end
        Rows_do1: DO J=pdims_s%j_start,pdims_s%j_end
          Row_length_do1: DO I=pdims_s%i_start,pdims_s%i_end
            UONP_SMALLHALO(I,J,K) = 0.
            VONP_SMALLHALO(I,J,K) = 0.
          END DO  Row_length_do1
        END DO  Rows_do1
      END DO  Levels_do1
!
      Levels_do1a: DO K=pdims%k_start,pdims%k_end
        Rows_do1a: DO J=pdims%j_start,pdims%j_end
          Row_length_do1a: DO I=pdims%i_start,pdims%i_end
! ----------------------------------------------------------------------+-------
!           Zero vertical divergence of pseudomomentum flux
! ----------------------------------------------------------------------+-------
            G_X(I,J,K) = 0.0
            G_Y(I,J,K) = 0.0
          END DO  Row_length_do1a
        END DO  Rows_do1a
!       NVIEW(K) = 0
!       NVIEW2(K) = 0
      END DO  Levels_do1a

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, jdir, ominrnbv, f_f,   &
!$OMP& jlev)

!
! ----------------------------------------------------------------------+-------
! 1.0   Set variables that are to be defined on all model levels
! ----------------------------------------------------------------------+-------
!

!$OMP DO SCHEDULE(STATIC) 
      Levels_do2: DO JLEV=KBASE,(tdims%k_end - 1)
! ----------------------------------------------------------------------+-------
! 1.1   Density, buoyancy frequency and altitude for middle levels
! ----------------------------------------------------------------------+-------
        Rows_do2: DO J=tdims%j_start,tdims%j_end
          Row_length_do2: DO I=tdims%i_start,tdims%i_end
            RHO_TH(I,J,JLEV) = RHONT_SMALLHALO(I,J,JLEV) * r_earth_radius_sq
!           Buoyancy (Brunt-Vaisala) frequency calculation
            NBV(I,J,JLEV) = ( g*(THETA(I,J,JLEV+1)-THETA(I,J,JLEV-1))/  &
     &                          (THETA(I,J,JLEV) *                      &
     &       (R_THETA_LEVELS(I,J,JLEV+1)-R_THETA_LEVELS(I,J,JLEV-1))) )
            NBV(I,J,JLEV) = MAX(NBV(I,J,JLEV), SQNMIN)
            NBV(I,J,JLEV) = SQRT(NBV(I,J,JLEV))
          END DO  Row_length_do2
        END DO  Rows_do2
      END DO  Levels_do2
!$OMP END DO

!
! ----------------------------------------------------------------------+-------
!     Density, Buoyancy (Brunt-Vaisala) frequency at top and bottom
! ----------------------------------------------------------------------+-------

!$OMP DO SCHEDULE(STATIC)
      Rows_do3: DO J=tdims%j_start,tdims%j_end
        Row_length_do3: DO I=tdims%i_start,tdims%i_end
!       Set rho and theta level value of Z_TH.
          RHO_TH(I,J,tkfix1start) = RHONT_SMALLHALO(I,J,tkfix1start) *  &
     &                              r_earth_radius_sq
          NBV(I,J,tkfix1start)    = NBV(I,J,KBASE)
          RHO_TH(I,J,tdims%k_end) = RHONT_SMALLHALO(I,J,tdims%k_end) *  &
     &                              r_earth_radius_sq
          NBV(I,J,tdims%k_end)    = NBV(I,J,tdims%k_end-1)
        END DO  Row_length_do3
      END DO  Rows_do3
!$OMP END DO
!

!$OMP SINGLE
      Levels_do4: DO JLEV=tdims%k_end,KBASE,-1
! ----------------------------------------------------------------------+-------
! 1.2   Set buoyancy frequency constant up to 1km altitude
! ----------------------------------------------------------------------+-------
        Rows_do4: DO J=tdims%j_start,tdims%j_end
          Row_length_do4: DO I=tdims%i_start,tdims%i_end
            IF ( (R_THETA_LEVELS(I,J,JLEV) - EARTH_RADIUS) <  1.0E3)    &
     &      NBV(I,J,JLEV-1) = NBV(I,J,JLEV)
          END DO  Row_length_do4
        END DO  Rows_do4
      END DO  Levels_do4
!$OMP END SINGLE
 
!
! ----------------------------------------------------------------------+-------
! 1.3  Compute component of wind U in each wave-propagation direction.
!      U is the dot product of (u,v) with k_0 but n.b. UDOTK is half-way
!      between rho levels.
!      Interpolate : [Vertical]   RHO onto T grid 
! ----------------------------------------------------------------------+-------

!$OMP DO SCHEDULE(STATIC) 
      IDir_do1: DO JDIR=1,IDIR
!
        Levels_do5: DO JLEV=tkfix1start,(tdims%k_end - 1)
          Rows_do5: DO J=tdims%j_start,tdims%j_end
            Row_length_do5: DO I=tdims%i_start,tdims%i_end
!           Assume theta levels are half way between rho levels.
!     WARNING: This is not very robust, preferable to recode with P_TO_T call,
!              but won't be a bit comparable change of course.
              UDOTK(I,J,JLEV,JDIR) =                                    &
     &        0.5*(UONP(I,J,JLEV) + UONP(I,J,JLEV+1))*COSPHI(JDIR) +    &
     &        0.5*(VONP(I,J,JLEV) + VONP(I,J,JLEV+1))*SINPHI(JDIR)
            END DO  Row_length_do5
          END DO  Rows_do5
        END DO  Levels_do5 
!
! ----------------------------------------------------------------------+-------
!      Set wind component for top level, to be equal to that on the top
!      Rho level, and total flux of horizontal pseudomomentum at bottom
! ----------------------------------------------------------------------+-------

        Rows_do5a: DO J=tdims%j_start,tdims%j_end
          Row_length_do5a: DO I=tdims%i_start,tdims%i_end
            UDOTK(I,J,tdims%k_end,JDIR) = UONP(I,J,pdims%k_end)         &
     &       *COSPHI(JDIR) + VONP(I,J,pdims%k_end)*SINPHI(JDIR)
            FPTOT(I,J,tkfix1start,JDIR) = 0.0
          END DO  Row_length_do5a
        END DO  Rows_do5a

!
! ----------------------------------------------------------------------+-------
! 2.0  Initialize variables that need to be defined up to launch level
! ----------------------------------------------------------------------+-------
!
!
! ----------------------------------------------------------------------+-------
! 3.0  Initialize gravity wave spectrum variables for the launch level
! ----------------------------------------------------------------------+-------
!
        Rows_do7: DO J=tdims%j_start,tdims%j_end
          Row_length_do7: DO I=tdims%i_start,tdims%i_end
! ----------------------------------------------------------------------+-------
!         Globally invariant value at Launch of Total vertical Flux of
!         horizontal pseudomomentum ... UMDP 34: 1.14
! ----------------------------------------------------------------------+-------
            RHOCL(I,J,JDIR) = RHO_TH(I,J,ILAUNCH) * GLOB_LAUNCH_FLUX
          END DO  Row_length_do7
        END DO  Rows_do7
      END DO  IDir_do1
!$OMP END DO

!
! ----------------------------------------------------------------------+-------
! 3.1 Compute minimum intrinsic frequency OMIN from inertial frequency
!     squared F_F. See UMDP 34, eqn. (A.7).
!-----------------------------------------------------------------------+-------

!$OMP DO SCHEDULE(STATIC)
      Rows_do7a: DO J=tdims%j_start,tdims%j_end
        Row_length_do7a: DO I=tdims%i_start,tdims%i_end
          F_F = (two_omega * SIN_THETA_LATITUDE(I,J))**2
          OMIN(I,J) = NBV(I,J,ILAUNCH) * BETA_EQ_RMSTAR
          OMIN(I,J) = MAX( OMIN(I,J), F_F )
          OMIN(I,J) = SQRT(OMIN(I,J))
        END DO  Row_length_do7a
      END DO  Rows_do7a
!$OMP END DO

!
! ----------------------------------------------------------------------+-------
! 4.0 Calculations carried out at levels from Launch to Lid
! ----------------------------------------------------------------------+-------
!

!$OMP DO SCHEDULE(STATIC)
      IDir_do1a: DO JDIR=1,IDIR
        Levels_do8: DO JLEV=KBASE,tdims%k_end
          Rows_do8: DO J=tdims%j_start,tdims%j_end
            Row_length_do8: DO I=tdims%i_start,tdims%i_end
! ----------------------------------------------------------------------+-------
!           Total vertical flux of horizontal pseudomomentum at launch level, 
!           analytic integral under curve = rho_l * C_l / m*^2 ... UMDP 34: 1.14
!           Note: the total flux at levels above is initialised to the sum of 
!           all sources propagated up through the column.
! ----------------------------------------------------------------------+-------
              IF (JLEV == LAUNCHLEV)  THEN
                FPTOT(I,J,JLEV,JDIR) = RHOCL(I,J,JDIR)
              ELSE
                FPTOT(I,J,JLEV,JDIR) = 0.0
              END IF
              FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV-1,JDIR) +           &
     &                               FPTOT(I,J,JLEV,JDIR)
            END DO  Row_length_do8
          END DO  Rows_do8
        END DO  Levels_do8
      END DO  IDir_do1a
!$OMP END DO

!
! ----------------------------------------------------------------------+-------
! 4.1 Compute [rho(z) . C(z)]_k / m*^2 scaling factor for quasi-saturated
!     spectrum (as per ... UMDP 34: 1.15,1.16)
! ----------------------------------------------------------------------+-------
!
!     IF (ABS(PSATM1) >= 0.1) THEN
!     For current setting of parameter psat this test is always true

!$OMP DO SCHEDULE(STATIC)
        Levels_do8a: DO JLEV=ILAUNCH,tdims%k_end
          Rows_do8a: DO J=tdims%j_start,tdims%j_end
            Row_length_do8a: DO I=tdims%i_start,tdims%i_end
              OMINRNBV = OMIN(I,J) / NBV(I,J,JLEV)
              FSATK_SCALE(I,J,JLEV) = RHO_TH(I,J,JLEV) * CCS0RMSTARSQ * &
     &         (NBV(I,J,JLEV))**2 * (OMINRNBV**PSATM1) *                &
     &         (1.0 - (OMINRNBV**TWOMPSAT)) / (1.0 - (OMINRNBV**PSATM1))
            END DO  Row_length_do8a
          END DO  Rows_do8a
        END DO  Levels_do8a
!$OMP END DO 

!$OMP END PARALLEL

!     ELSE
!     Require a different functional form for normalisation factor B0
!             BBS = 1.0 / ALOG(NBV(I,J,JLEV) / OMIN(I,J))
!     END IF
! ----------------------------------------------------------------------+-------
!     Loop over directions and levels and calculate horizontal component
!     of the vertical flux of pseudomomentum for each azimuthal
!     direction and for each altitude
! ----------------------------------------------------------------------+-------

! gives each thread the largest block possible to execute


! Parameters used: psat, mnlmin, maxwhile, mstar, twomtt, ttm2, rttm2, 
! ttm1, rssptt, ssptt
!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(jlev, jdir, jj, i,              & 
!$OMP& j, mkill, mnly, indexjd, indexj, indexi, nni, fplus, gplus,      & 
!$OMP& gminus, atte_c, curv_c, fterm, gterm, wgtn, mnlx, l_chop2,       &
!$OMP& jwhile, nchop2, l_ftheng, nnjd, nnj, omp_block)

        OMP_BLOCK = IDIR
!$      OMP_BLOCK = CEILING(REAL(IDIR)/omp_get_num_threads())

!$OMP DO SCHEDULE(STATIC)
        omp_blocking1: DO JJ=1, IDIR, OMP_BLOCK 

           Levels_do92: DO JLEV=ILAUNCH+1,tdims%k_end
              L_CHOP2 = .FALSE.
              IDir_do2a: DO JDIR=JJ,MIN(JJ+OMP_BLOCK-1, IDIR)
                 Rows_do92: DO J=tdims%j_start,tdims%j_end
                    Row_length_do92: DO I=tdims%i_start,tdims%i_end

! ----------------------------------------------------------------------+-------
!     Initialise MGUESS (start point for iterative searches if needed)
! ----------------------------------------------------------------------+-------
                       MGUESS_a(I,J,JLEV,JDIR) = 0.0
                       Fptot_if1: IF (FPTOT(I,J,JLEV-1,JDIR) >  0.0) THEN
! ----------------------------------------------------------------------+-------
! 4.2       Calculate variables that define the Chop Type Cases.
! ----------------------------------------------------------------------+-------
                          DDU_a(I,J,JLEV,JDIR) = UDOTK(I,J,ILAUNCH,JDIR)&
                               &                  - UDOTK(I,J,JLEV,JDIR)
! ----------------------------------------------------------------------+-------
!             UMDP 34: 1.23 coefficient B
!             Using ratio of flux scalings rather than densities is more 
!             robust because total flux imported from other schemes may 
!             have different relationship to launch density than 1.14
!             Note: the total flux is initialised to the sum of all sources
!             propagated up through the column.
! ----------------------------------------------------------------------+-------
              ATTENUATION(I,J,JLEV,JDIR) =                              &
     &           ( FSATK_SCALE(I,J,JLEV) / FPTOT(I,J,JLEV,JDIR) ) *     &
     &           ( NBV(I,J,ILAUNCH) / NBV(I,J,JLEV) )**TTM1
! ----------------------------------------------------------------------+-------
!             UMDP 34: 1.22 coefficient A = (A/B) * B
! ----------------------------------------------------------------------+-------
              CURVATURE(I,J,JLEV,JDIR)  =  DDU_a(I,J,JLEV,JDIR) * MSTAR &
     &                                     / NBV(I,J,ILAUNCH)
!
              ACOEFF(I,J,JLEV,JDIR)     =  CURVATURE(I,J,JLEV,JDIR) *   &
     &                                   ATTENUATION(I,J,JLEV,JDIR)
!
              MINTERCEPT(I,J,JLEV,JDIR) = ATTENUATION(I,J,JLEV,JDIR)    &
     &         * ( (MNLMIN * CURVATURE(I,J,JLEV,JDIR)) + 1.0 )**TTM2
!
              Curv_if1: IF (CURVATURE(I,J,JLEV,JDIR) <  ASECN)  THEN
! ----------------------------------------------------------------------+-------
!             Negative Doppler Shift : factor will hit zero (kill point)
! ----------------------------------------------------------------------+-------
                MKILL = 1.0 / ABS(CURVATURE(I,J,JLEV,JDIR))
!
                Mkill_if1: IF (MKILL <= MNLMIN )  THEN
! ----------------------------------------------------------------------+-------
!               Chop Type IV : No flux propagates
! ----------------------------------------------------------------------+-------
                  FPTOT(I,J,JLEV,JDIR) = 0.0
                ELSE
                  IF (MKILL >  1.0)  THEN
                  INTERCEPT1(I,J,JLEV,JDIR) = ATTENUATION(I,J,JLEV,JDIR)&
     &              * ( 1.0 + CURVATURE(I,J,JLEV,JDIR) )**TTM2
                  ELSE
! ----------------------------------------------------------------------+-------
!                 Doppler factor minimum (kill point) situated below
!                 mstar in the low-m part of the launch spectrum
! ----------------------------------------------------------------------+-------
                    INTERCEPT1(I,J,JLEV,JDIR) = 0.0
                  END IF
!
                  Lowend_if1: IF (INTERCEPT1(I,J,JLEV,JDIR) >= 1.0) THEN
! ----------------------------------------------------------------------+-------
!                 Chop Type I: Intersection in high wavenumber part only
! ----------------------------------------------------------------------+-------
                    MNLY = ( ATTENUATION(I,J,JLEV,JDIR)**TTRAT -        &
     &              ATTENUATION(I,J,JLEV,JDIR)) / ACOEFF(I,J,JLEV,JDIR)
!
                    FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) *       &
     &     (1.0 - (A0_R_1MT * CURVATURE(I,J,JLEV,JDIR) * MNLY**TWOMTT))
                  ELSE
                    IF (MINTERCEPT(I,J,JLEV,JDIR) <= FMINUS)  THEN
! ----------------------------------------------------------------------+-------
!                 Chop Type IIb: Low wavenumber intersect only below min
! ----------------------------------------------------------------------+-------
                      FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) *     &
     &            A0_R_SP1TM1 * TAIL_CHOP2B * MINTERCEPT(I,J,JLEV,JDIR) &
     &                 * ( (MNLMIN * CURVATURE(I,J,JLEV,JDIR)) + 1.0 )
                    ELSE
! ----------------------------------------------------------------------+-------
!                 Chop Type IIa: Low wavenumber intersect only
! ----------------------------------------------------------------------+-------
                      L_CHOP2 = .TRUE.
                      MGUESS_a(I,J,JLEV,JDIR) = MIN(MKILL, 1.0)
                      FPFAC(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) * A0_R_SP1TM1
                      FPTOT(I,J,JLEV,JDIR) = 0.0
                    END IF
                  END IF  Lowend_if1
!
                END IF  Mkill_if1
!
              ELSE IF (CURVATURE(I,J,JLEV,JDIR) >  ASECP)  THEN
! ----------------------------------------------------------------------+-------
!             Positive Doppler Shift : non-zero factor (no kill point)
! ----------------------------------------------------------------------+-------
                INTERCEPT1(I,J,JLEV,JDIR) = ATTENUATION(I,J,JLEV,JDIR)  &
     &            * ( 1.0 + CURVATURE(I,J,JLEV,JDIR) )**TTM2
!
                Chop3_if1: IF (INTERCEPT1(I,J,JLEV,JDIR) <  1.0)  THEN
! ----------------------------------------------------------------------+-------
!               Chop Type III: Intersection in both wavenumber parts
! ----------------------------------------------------------------------+-------
                  FPFAC(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) * A0_R_SP1TM1
!
! ----------------------------------------------------------------------+-------
!                 First find intersect in high wavenumber part
!                 UMDP 34: 1.25
! ----------------------------------------------------------------------+-------
                  MNLY = ( ATTENUATION(I,J,JLEV,JDIR)**TTRAT -          &
     &              ATTENUATION(I,J,JLEV,JDIR)) / ACOEFF(I,J,JLEV,JDIR)
!
                  FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) *         &
     &              A0_R_1MT * CURVATURE(I,J,JLEV,JDIR) * MNLY**TWOMTT
!
! ----------------------------------------------------------------------+-------
!                 Then find intersect in low wavenumber part to reckon
!                 its flux contribution for addition when available
! ----------------------------------------------------------------------+-------
                  IF (MINTERCEPT(I,J,JLEV,JDIR) <= FMINUS)  THEN
! ----------------------------------------------------------------------+-------
!                 Chop Type IIIb: Low wavenumber intersect below min
! ----------------------------------------------------------------------+-------
                    FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) +       &
     &                                   ( FPFAC(I,J,JLEV,JDIR) *       &
     &                 TAIL_CHOP2B *  MINTERCEPT(I,J,JLEV,JDIR) *       &
     &                ( (MNLMIN * CURVATURE(I,J,JLEV,JDIR)) + 1.0 ) )
                  ELSE
! ----------------------------------------------------------------------+-------
!                 Chop Type IIIa: Low wavenumber intersect
! ----------------------------------------------------------------------+-------
                    L_CHOP2 = .TRUE.
                    MGUESS_a(I,J,JLEV,JDIR) = 1.0
                  END IF
!
! ----------------------------------------------------------------------+-------
!               ELSE Chop Type 0: No intersection (spectrum unaltered)
! ----------------------------------------------------------------------+-------
                END IF  Chop3_if1
              ELSE
! ----------------------------------------------------------------------+-------
!             Negligible Doppler shift
! ----------------------------------------------------------------------+-------
!               Strictly this is analytic solution mNLX.  UMDP 34: 1.27
                MNLY = ATTENUATION(I,J,JLEV,JDIR)**RSSPTT
                IF (MNLY <= MNLMIN)  THEN
! ----------------------------------------------------------------------+-------
!               Chop Type IIb: Low wavenumber intersect only below min
! ----------------------------------------------------------------------+-------
                  FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) *         &
     &           A0_R_SP1TM1 * TAIL_CHOP2B * ATTENUATION(I,J,JLEV,JDIR)
                ELSE
                  IF (MNLY <  1.0)  FPTOT(I,J,JLEV,JDIR) =              &
! ----------------------------------------------------------------------+-------
!               Chop Type IIc: Low wavenumber intersect only (analytic)
! ----------------------------------------------------------------------+-------
     &              FPTOT(I,J,JLEV,JDIR) * A0_R_SP1TM1 *                &
     &               ( (SSPTT * (MNLY**SSP1)) - HEAD_CHOP2A )
! ----------------------------------------------------------------------+-------
!                 ELSE Chop Type 0: No intersection (spectrum unaltered)
! ----------------------------------------------------------------------+-------
                END IF
              END IF  Curv_if1
!
              END IF  Fptot_if1
            END DO  Row_length_do92
          END DO  Rows_do92
        END DO  IDir_do2a
!
        Lchop2_if1: IF (L_CHOP2)  THEN
! ----------------------------------------------------------------------+-------
!       Process low wavenumber contribution: evaluate intersect mNX
! ----------------------------------------------------------------------+-------
          NCHOP2 = 0
!
          IDir_do2b: DO JDIR=JJ,MIN(JJ+OMP_BLOCK-1, IDIR)
            Rows_do93: DO J=tdims%j_start,tdims%j_end
              Row_length_do93: DO I=tdims%i_start,tdims%i_end
                IF (MGUESS_a(I,J,JLEV,JDIR) >  0.0)  THEN
                  NCHOP2 = NCHOP2 + 1
!
                  INDEXJD(NCHOP2) = JDIR
                  INDEXJ(NCHOP2)  = J
                  INDEXI(NCHOP2)  = I
                END IF
              END DO  Row_length_do93
            END DO  Rows_do93
          END DO  IDir_do2b
!         NVIEW(JLEV) = NCHOP2
!         NSPIRA = 0
!
          Nchop2_do1: DO I=1,NCHOP2
! ----------------------------------------------------------------------+-------
!               Chop Type IIa : / Full solution required for mNlX
!          or   Chop Type IIIa: ! ----------------------------------------------------------------------+-------

            NNJD = INDEXJD(I)
            NNJ  = INDEXJ(I)
            NNI  = INDEXI(I)
!
            FPLUS  = MGUESS_a(NNI,NNJ,JLEV,NNJD)**SSPTT
            GPLUS  = INTERCEPT1(NNI,NNJ,JLEV,NNJD)
!           FMINUS = MNLMIN**SSPTT    Defined as a constant
            GMINUS = MINTERCEPT(NNI,NNJ,JLEV,NNJD)
            ATTE_C(I) = ATTENUATION(NNI,NNJ,JLEV,NNJD)
            CURV_C(I) = CURVATURE(NNI,NNJ,JLEV,NNJD)
!
            FTERM = ( ((FMINUS / ATTE_C(I))**RTTM2) - 1.0 ) / CURV_C(I)
            GTERM = GMINUS**RSSPTT
            L_FTHENG(I) = .FALSE.
!
            Curv_if2: IF (CURVATURE(NNI,NNJ,JLEV,NNJD) >  ASECP)  THEN
! ----------------------------------------------------------------------+-------
!           Positive Doppler Shift
! ----------------------------------------------------------------------+-------
              WGTN(I) = 0.0
            ELSE
! ----------------------------------------------------------------------+-------
!           Negative Doppler Shift
! ----------------------------------------------------------------------+-------
              WGTN(I) = 1.0 - MWEIGHT
!
              IF (FPLUS <= GMINUS  .AND.  GPLUS >  FMINUS)  THEN
                FTERM = (((FPLUS / ATTE_C(I))**RTTM2) - 1.0)/ CURV_C(I)
                GTERM = GPLUS**RSSPTT
                L_FTHENG(I) = (GTERM  <   FTERM)
!
              ELSE IF (FPLUS >  GMINUS  .AND.  GPLUS <= FMINUS)  THEN
                L_FTHENG(I) = (GTERM >= FTERM)
!
              ELSE IF (FPLUS <= GMINUS  .AND.  GPLUS <= FMINUS)  THEN
                L_FTHENG(I) = .TRUE.
!
!             ELSE Use default settings
              END IF
            END IF  Curv_if2
!
            IF (L_FTHENG(I))  THEN
!             NSPIRA = NSPIRA + 1
              MNLX(I,0) = FTERM
            ELSE
              MNLX(I,0) = GTERM
            END IF
          END DO  Nchop2_do1
!         NVIEW2(JLEV) = NSPIRA
!
          Jwhile_do2: DO JWHILE=0,MAXWHILE-1
            Nchop2_do2: DO I=1,NCHOP2
!
              IF (L_FTHENG(I))  THEN
! ----------------------------------------------------------------------+-------
!           Obtain m_n+1 from g_n+1  = f_n (m_n)
! ----------------------------------------------------------------------+-------
                MNLX(I,JWHILE+1) = (                                    &
     &           (((MNLX(I,JWHILE)**SSPTT) / ATTE_C(I))**RTTM2) - 1.0 ) &
     &           / CURV_C(I)
              ELSE
! ----------------------------------------------------------------------+-------
!           Obtain m_n+1 from f_n+1  = g_n (m_n)
! ----------------------------------------------------------------------+-------
                MNLX(I,JWHILE+1) = ( (ATTE_C(I) *                       &
     &          ((1.0 + (CURV_C(I) * MNLX(I,JWHILE)))**TTM2))**RSSPTT )
              END IF
!
              MNLX(I,JWHILE+1) = ((1.0 - WGTN(I)) * MNLX(I,JWHILE+1)) + &
     &                                  (WGTN(I)  * MNLX(I,JWHILE))
!
            END DO  Nchop2_do2
          END DO  Jwhile_do2
!
!CDIR NODEP
          Nchop2_do3: DO I=1,NCHOP2
            NNJD = INDEXJD(I)
            NNJ  = INDEXJ(I)
            NNI  = INDEXI(I)
!
            FPTOT(NNI,NNJ,JLEV,NNJD) = FPTOT(NNI,NNJ,JLEV,NNJD) +       &
     &     (FPFAC(NNI,NNJ,JLEV,NNJD) * ( ((MNLX(I,MAXWHILE)**SSP1) *    &
     &      ( SSPTT + (SSP1 * MNLX(I,MAXWHILE) * CURV_C(I)) )) - HEAD_CHOP2A ))
!
          END DO  Nchop2_do3
        END IF  Lchop2_if1
!
        IDir_do2c: DO JDIR=JJ,MIN(JJ+OMP_BLOCK-1, IDIR)
          Rows_do10: DO J=tdims%j_start,tdims%j_end
            Row_length_do10: DO I=tdims%i_start,tdims%i_end
!-----------------------------------------------------------------------+-------
!         Now correct pseudomomentum flux in the evolved spectrum if the
!         new value is non-physical (pseudomomentum flux cannot increase
!         with altitude)
!-----------------------------------------------------------------------+-------
              FPTOT(I,J,JLEV,JDIR) =                                    &
     &          MIN(FPTOT(I,J,JLEV,JDIR), FPTOT(I,J,JLEV-1,JDIR))
!
            END DO  Row_length_do10
          END DO  Rows_do10
        END DO  IDir_do2c
!
      END DO  Levels_do92
   END DO omp_blocking1
!$OMP END DO

!$OMP END PARALLEL


!
! ----------------------------------------------------------------------+-------
! 4.5   If choosing Opaque Upper Boundary set fluxes to zero at top
! ----------------------------------------------------------------------+-------
      IF (L_USSP_OPAQUE) THEN
        IDir_do3: DO JDIR=1,IDIR
          Rows_do12: DO J=tdims%j_start,tdims%j_end
            Row_length_do12: DO I=tdims%i_start,tdims%i_end
              FPTOT(I,J,tdims%k_end,JDIR) =  0.0
            END DO  Row_length_do12
          END DO  Rows_do12
        END DO  IDir_do3
      ENDIF


!$OMP  PARALLEL DEFAULT(NONE) SHARED(fptot, g_x, g_y, cosphi,           &
!$OMP& p_layer_boundaries, uonp_smallhalo, vonp_smallhalo, timestep,    &
!$OMP& sinphi, pdims, ilaunch, L_ussp_heating, r_theta_levels,          &
!$OMP& r_rho_levels, uonp, vonp, T_inc, tdims)                          &
!$OMP& PRIVATE(g_g, i, j, k, jlev, jdir, jj, dzb, dzu, dzl, uhat, vhat, &
!$OMP& ududt, vdvdt, omp_block)

!
! ----------------------------------------------------------------------+-------
! 5.0   Compute vertical divergence of pseudomomentum flux.
!       Column Integral Converts : [Vertical]   T onto RHO grid 
! ----------------------------------------------------------------------+-------

! gives each thread the largest block possible to execute

   OMP_BLOCK = pdims%k_end-(ILAUNCH+1)+1
!$  OMP_BLOCK = CEILING(REAL(((pdims%k_end - (ILAUNCH+1)) + 1))/        &
!$& omp_get_num_threads())

!$OMP DO SCHEDULE(STATIC)
      omp_blocking2: DO JJ=ILAUNCH+1, pdims%k_end, OMP_BLOCK
        Levels_do14: DO JLEV=JJ, MIN(JJ+OMP_BLOCK-1,pdims%k_end)
          Rows_do14: DO J=pdims%j_start,pdims%j_end
            Row_length_do14: DO I=pdims%i_start,pdims%i_end
              IDir_do4: DO JDIR=1,IDIR
!             Pseudomomentum flux
                G_G =                                                   &
     &           g * (FPTOT(I,J,JLEV,JDIR) - FPTOT(I,J,JLEV-1,JDIR)) /  &
     &            (P_LAYER_BOUNDARIES(I,J,JLEV) -                       &
     &            P_LAYER_BOUNDARIES(I,J,JLEV-1))
                G_X(I,J,JLEV)= G_X(I,J,JLEV)+G_G* COSPHI(JDIR)
                G_Y(I,J,JLEV)= G_Y(I,J,JLEV)+G_G* SINPHI(JDIR)
              END DO  IDir_do4
            END DO  Row_length_do14
          END DO  Rows_do14
        END DO  Levels_do14
      END DO  omp_blocking2
!$OMP END DO

!
! ----------------------------------------------------------------------+-------
!     WRITE(6,*) ' Level of Launching : ',ILAUNCH
!     WRITE(6,*) ' Number of Specials : ',NVIEW
!     WRITE(6,*) ' Number of SpiralAs : ',NVIEW2
! ----------------------------------------------------------------------+-------
!
! ----------------------------------------------------------------------+-------
! 5.1   Wind and temperature increments from wave dissipation
! ----------------------------------------------------------------------+-------

!$OMP DO SCHEDULE(STATIC)
      Levels_do15: DO JLEV=ILAUNCH,pdims%k_end
        Rows_do15: DO J=pdims%j_start,pdims%j_end
          Row_length_do15: DO I=pdims%i_start,pdims%i_end
            UONP_SMALLHALO(I,J,JLEV) = TIMESTEP * G_X(I,J,JLEV)
            VONP_SMALLHALO(I,J,JLEV) = TIMESTEP * G_Y(I,J,JLEV)
          END DO  Row_length_do15
        END DO  Rows_do15
      END DO  Levels_do15
!$OMP END DO



!-----------------------------------------------------------------
! Calculate heating due to gravity wave dissipation
!-----------------------------------------------------------------
      GW_heating: IF ( L_ussp_heating ) THEN

!$OMP DO SCHEDULE(STATIC)
       Levels_do16: DO k = tkfix1start,tdims%k_end-1
         Rows_do16: DO j = tdims%j_start,tdims%j_end
           Row_length_do16: DO i = tdims%i_start,tdims%i_end

           dzb = r_theta_levels(i,j,k) -  r_rho_levels(i,j,k)
           dzu = r_rho_levels(i,j,k+1) -  r_theta_levels(i,j,k)
           dzl = r_rho_levels(i,j,k+1) -  r_rho_levels(i,j,k)

!          u and v on theta_level(k)
           uhat   = dzu * uonp(i,j,k) + dzb * uonp(i,j,k+1)
           vhat   = dzu * vonp(i,j,k) + dzb * vonp(i,j,k+1)

!          u*du and v*dv on theta_level(k)
           ududt  = uhat *( dzu * g_x(i,j,k) + dzb * g_x(i,j,k+1) )
           vdvdt  = vhat *( dzu * g_y(i,j,k) + dzb * g_y(i,j,k+1) )

!          dT/dt on theta_level(k)
           T_inc(i,j,k)=T_inc(i,j,k) - timestep * ( ududt + vdvdt ) /   &
                       ( cp * dzl * dzl )

         END DO  Row_length_do16
        END DO  Rows_do16
       END DO  Levels_do16
!$OMP END DO

      END IF GW_heating

!$OMP END PARALLEL

!
! ----------------------------------------------------------------------+-------
!     Put U,V increments onto U,V grid after initialising increments.
!     Interpolate : [Horizontal] P onto U,V grid
! ----------------------------------------------------------------------+-------
      UINC(:,:,:)=0.
      VINC(:,:,:)=0.
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(UONP_SMALLHALO(pdims_s%i_start,pdims_s%j_start,1)&
     &                 ,ROW_LENGTH,ROWS,LEVELS,                         &
     &                OFF_X,OFF_Y,FLD_TYPE_P,.FALSE.)

      CALL p_to_u(uonp_smallhalo,                                       &
                   pdims_s%i_start,pdims_s%i_end,                       &
                   pdims_s%j_start,pdims_s%j_end,                       &
                   udims%i_start,udims%i_end,                           &
                   udims%j_start,udims%j_end,                           &
                   1,levels,uinc)      
!
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(VONP_SMALLHALO(pdims_s%i_start,pdims_s%j_start,1)&
     &                ,ROW_LENGTH,ROWS,LEVELS,                          &
     &                OFF_X,OFF_Y,FLD_TYPE_P,.FALSE.)

      CALL p_to_v(vonp_smallhalo,                                       &
                   pdims_s%i_start,pdims_s%i_end,                       &
                   pdims_s%j_start,pdims_s%j_end,                       &
                   vdims%i_start,vdims%i_end,                           &
                   vdims%j_start,vdims%j_end,                           &
                   1,levels,vinc) 
                   
IF (l_vatpoles) THEN
! correct stress at poles
      IF( model_domain == mt_global ) THEN
        IF( at_extremity(psouth) ) THEN


          CALL eg_v_at_poles(uinc, vinc, 1.0,                     &
                             udims%j_start, vdims%j_start,        &
                             udims,vdims)
        
        END IF
        
        IF( at_extremity(pnorth) ) THEN


          CALL eg_v_at_poles(uinc, vinc, -1.0,                    &
                             udims%j_end, vdims%j_end,            &
                             udims,vdims)
        
        END IF           
      END IF
END IF ! vatpoles        
!
! ----------------------------------------------------------------------+-------
! Add increments to wind and temperature
! ----------------------------------------------------------------------+-------
      DO K=ILAUNCH,udims%k_end
!     NB: assumes that U and V winds will always be on same vertical levels
        DO J=udims%j_start,udims%j_end
          DO I=udims%i_start,udims%i_end
            R_U(I,J,K)=R_U(I,J,K)+UINC(I,J,K)
          ENDDO
        ENDDO
        DO J=vdims%j_start,vdims%j_end
          DO I=vdims%i_start,vdims%i_end
            R_V(I,J,K)=R_V(I,J,K)+VINC(I,J,K)
          ENDDO
        ENDDO
      ENDDO

!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Northward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
      GWspec_Flux_if1n: IF(GWSPEC_NFLUX_ON) THEN
        Levels_do20n: DO K=tkfix1start,tdims%k_end
          Rows_do20n: DO J=tdims%j_start,tdims%j_end
            Row_length_do20n: DO I=tdims%i_start,tdims%i_end
              Work_TsmallHALO(I,J,K) = FPTOT(I,J,K,1)
            END DO  Row_length_do20n
          END DO  Rows_do20n
        END DO  Levels_do20n
!
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &       Work_TsmallHALO(tdims_s%i_start,tdims_s%j_start,1),        &
     &       ROW_LENGTH, ROWS,                                          &
     &    LEVELS,OFF_X,OFF_Y,FLD_TYPE_P,.FALSE.)
!

!     Interpolate : [Horizontal] P onto V grid
!
      CALL p_to_v(work_TsmallHALO(tdims_s%i_start,                      & 
                                  tdims_s%j_start,1),                   &
                   pdims_s%i_start,pdims_s%i_end,                       &
                   pdims_s%j_start,pdims_s%j_end,                       &
                   vdims%i_start,vdims%i_end,                           &
                   vdims%j_start,vdims%j_end,                           &
                   1,levels,gwspec_nflux(vdims%i_start,vdims%j_start,1))     
IF (l_vatpoles) THEN
! set polar rows to common mean value. 
      CALL polar_row_mean(                                        & 
                      gwspec_nflux(vdims%i_start,vdims%j_start,1),&
                      vdims%i_start,vdims%i_end,                  &
                      vdims%j_start,vdims%j_end,                  &
                      1,levels,                                   &
                      global_row_length,                          &
                      n_proc, n_procy, proc_row_group,            &
                      at_extremity)
END IF ! vatpoles
     
      END IF  GWspec_Flux_if1n
!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Westward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
      GWspec_Flux_if1w: IF(GWSPEC_WFLUX_ON .or. GWSPEC_WFLUX_P_ON) THEN
        Levels_do20w: DO K=tkfix1start,tdims%k_end
          Rows_do20w: DO J=tdims%j_start,tdims%j_end
            Row_length_do20w: DO I=tdims%i_start,tdims%i_end
              Work_TsmallHALO(I,J,K) = FPTOT(I,J,K,2)
            END DO  Row_length_do20w
          END DO  Rows_do20w
        END DO  Levels_do20w
!
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &       Work_TsmallHALO(tdims_s%i_start,tdims_s%j_start,1),        &
     &       ROW_LENGTH, ROWS,                                          &
     &    LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
!
!     Interpolate : [Horizontal] P onto U grid
!
      CALL p_to_u(                                                      &
             Work_TsmallHALO(tdims_s%i_start,tdims_s%j_start,1),        &
                   pdims_s%i_start,pdims_s%i_end,                       &
                   pdims_s%j_start,pdims_s%j_end,                       &
                   udims%i_start,udims%i_end,                           &
                   udims%j_start,udims%j_end,                           & 
                   1,levels,gwspec_wflux(udims%i_start,udims%j_start,1))     
      END IF  GWspec_Flux_if1w
!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Southward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
      GWspec_Flux_if1s: IF(GWSPEC_SFLUX_ON) THEN
        Levels_do20s: DO K=tkfix1start,tdims%k_end
          Rows_do20s: DO J=tdims%j_start,tdims%j_end
            Row_length_do20s: DO I=tdims%i_start,tdims%i_end
              Work_TsmallHALO(I,J,K) = FPTOT(I,J,K,3)
            END DO  Row_length_do20s
          END DO  Rows_do20s
        END DO  Levels_do20s
!
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &       Work_TsmallHALO(tdims_s%i_start,tdims_s%j_start,1),        &
     &       ROW_LENGTH, ROWS,                                          &
     &    LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
!
!     Interpolate : [Horizontal] P onto V grid
!
      CALL p_to_v(                                                      &
             Work_TsmallHALO(tdims_s%i_start,tdims_s%j_start,1),        &
                   pdims_s%i_start,pdims_s%i_end,                       &
                   pdims_s%j_start,pdims_s%j_end,                       &
                   vdims%i_start,vdims%i_end,                           &
                   vdims%j_start,vdims%j_end,                           &
                   1,levels,gwspec_sflux(vdims%i_start,vdims%j_start,1))
IF (l_vatpoles) THEN
! set polar rows to common mean value. 
      CALL polar_row_mean(                                        & 
                      gwspec_sflux(vdims%i_start,vdims%j_start,1),&
                      vdims%i_start,vdims%i_end,                  &
                      vdims%j_start,vdims%j_end,                  &
                      1,levels,                                   &
                      global_row_length,                          &
                      n_proc, n_procy, proc_row_group,            &
                      at_extremity)
END IF ! vatpoles        
                   
      END IF  GWspec_Flux_if1s
!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Eastward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
      GWspec_Flux_if1e: IF(GWSPEC_EFLUX_ON .or. GWSPEC_EFLUX_P_ON) THEN
        Levels_do20e: DO K=tkfix1start,tdims%k_end
          Rows_do20e: DO J=tdims%j_start,tdims%j_end
            Row_length_do20e: DO I=tdims%i_start,tdims%i_end
              Work_TsmallHALO(I,J,K) = FPTOT(I,J,K,4)
            END DO  Row_length_do20e
          END DO  Rows_do20e
        END DO  Levels_do20e
!
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &       Work_TsmallHALO(tdims_s%i_start,tdims_s%j_start,1),        &
     &       ROW_LENGTH, ROWS,                                          &
     &    LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
!
!     Interpolate : [Horizontal] P onto U grid
!
      CALL p_to_u(                                                      &
             Work_TsmallHALO(tdims_s%i_start,tdims_s%j_start, 1),       &
                   pdims_s%i_start,pdims_s%i_end,                       &
                   pdims_s%j_start,pdims_s%j_end,                       &
                   udims%i_start,udims%i_end,                           &
                   udims%j_start,udims%j_end,                           & 
                   1,levels,gwspec_eflux(udims%i_start,udims%j_start,1))     
      END IF  GWspec_Flux_if1e
!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Acceleration of Zonal Wind (on rho levels)
! ----------------------------------------------------------------------+-------
      GWspec_Acc_if1ew: IF(GWSPEC_EWACC_ON .or. GWSPEC_EWACC_P_ON) THEN
        Levels_do20ew: DO K=pdims%k_start,pdims%k_end
          Rows_do20ew: DO J=pdims%j_start,pdims%j_end
            Row_length_do20ew: DO I=pdims%i_start,pdims%i_end
              Work_PsmallHALO(I,J,K) = G_X(I,J,K)
            END DO  Row_length_do20ew
          END DO  Rows_do20ew
        END DO  Levels_do20ew
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &       Work_PsmallHALO(pdims_s%i_start,pdims_s%j_start,1),        &
     &       ROW_LENGTH, ROWS,                                          &
     &      LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
!
!     Interpolate : [Horizontal] P onto U grid
!
      CALL p_to_u(                                                      &
             Work_PsmallHALO(pdims_s%i_start,pdims_s%j_start,1),        &
                   pdims_s%i_start,pdims_s%i_end,                       &
                   pdims_s%j_start,pdims_s%j_end,                       &
                   udims%i_start,udims%i_end,                           &
                   udims%j_start,udims%j_end,                           &
                   1,levels,gwspec_ewacc)     
      END IF GWspec_Acc_if1ew
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Acceleration of Meridional Wind (on rho levels)
! ----------------------------------------------------------------------+-------
      GWspec_Acc_if1ns: IF(GWSPEC_NSACC_ON) THEN
        Levels_do20ns: DO K=pdims%k_start,pdims%k_end
          Rows_do20ns: DO J=pdims%j_start,pdims%j_end
            Row_length_do20ns: DO I=pdims%i_start,pdims%i_end
              Work_PsmallHALO(I,J,K) = G_Y(I,J,K)
            END DO  Row_length_do20ns
          END DO  Rows_do20ns
        END DO  Levels_do20ns
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &       Work_PsmallHALO(pdims_s%i_start,pdims_s%j_start,1),        &
     &       ROW_LENGTH, ROWS,                                          &
     &      LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
!
!     Interpolate : [Horizontal] P onto V grid
!
      CALL p_to_v(                                                      &
             Work_PsmallHALO(pdims_s%i_start,pdims_s%j_start,1),        &
                   pdims_s%i_start,pdims_s%i_end,                       &
                   pdims_s%j_start,pdims_s%j_end,                       &
                   vdims%i_start,vdims%i_end,                           &
                   vdims%j_start,vdims%j_end,                           &
                   1,levels,gwspec_nsacc)
                   
IF (l_vatpoles) THEN
! correct stress at poles
      IF( model_domain == mt_global ) THEN
        IF( at_extremity(psouth) ) THEN


          CALL eg_v_at_poles(gwspec_ewacc, gwspec_nsacc, 1.0,     &
                             udims%j_start, vdims%j_start,        &
                             udims,vdims)
        
        END IF
        
        IF( at_extremity(pnorth) ) THEN


          CALL eg_v_at_poles(gwspec_ewacc, gwspec_nsacc, -1.0,    &
                             udims%j_end, vdims%j_end,            &
                             udims,vdims)
        
        END IF           
      END IF
END IF ! vatpoles        
                   
      END IF  GWspec_Acc_if1ns

      IF (lhook) CALL dr_hook('GW_USSP',zhook_out,zhook_handle)
      RETURN
!
      END SUBROUTINE GW_USSP

END MODULE GW_USSP_MOD
