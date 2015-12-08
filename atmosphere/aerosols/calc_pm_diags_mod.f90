! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the total PM10 and PM2.5 concentrations and for each aerosol type

MODULE calc_pm_diags_mod

CONTAINS

SUBROUTINE calc_pm_diags (                                                     &
  ! Arguments IN
  ! Parallel variables
  off_x, off_y,                                                                &
  !   Model dimensions
  row_length, rows,                                                            &
  model_levels,                                                                &
  salt_dim1, salt_dim2, salt_dim3,                                             &
  !   Logicals IN
  l_sulpc_so2, l_soot, l_biomass, l_ocff, l_use_biogenic, l_dust,              &
  l_nitrate, l_use_seasalt_pm,                                  &
  !   Data fields IN
  p_theta_levels, t,                                                           &
  so4_ait, so4_acc, soot_new, soot_agd, bmass_new, bmass_agd,                  &
  ocff_new, ocff_agd, biogenic, sea_salt_film, sea_salt_jet,                   &
  dust_1, dust_2, dust_3, dust_4, dust_5, dust_6, nitr_acc,                    &
  ! Arguments OUT
  pm10, pm2p5, pm10_so4, pm2p5_so4, pm10_bc, pm2p5_bc, pm10_bb, pm2p5_bb,      &
  pm10_ocff, pm2p5_ocff, pm10_soa, pm2p5_soa,  pm10_ss, pm2p5_ss,              &
  conc_dust, pm10_dust, pm2p5_dust, pm10_nitr, pm2p5_nitr)

USE atmos_constants_mod, ONLY: r
USE conversions_mod, ONLY: pi
USE dust_parameters_mod, ONLY: l_twobin_dust

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


!--------------------------------------------------------------------------
! Description : 
! Calculates the PM10 and PM2.5 mass concentrations (in microgram/m3)
! Method : 
!     Calculation of the PM10 and PM2.5 mass concentrations (in
!     microgram/m3) by using the cumulative volume (from 0 to
!     cut-off diameter) of a log-normal, with median radius
!     (or gemetric mean radius) r_bar, geometric standard deviation
!     sigma and total number Ntot.
!     The contributions by the different aerosol species and modes
!     are added up to get the total PM10 & PM2.5 concs (microgram/m3).
!
!     Formulation derived from eqs. (7.51) & (7.51a) of Seinfeld
!     and Pandis (1998) using method of (7.39). Also, checked with
!     "The Mechanics of Inhaled Pharmaceutical Aerosols: An
!     Introduction", by Warren A. Finlay, Academic Press, 2001
!
!     (1) Cumulative_volume for particles with diameter between
!     0 and d_cut (assuming spherical particles):
!            V(d_cut) = (A*Ntot/2.0) + (A*Ntot/2.0) *
!                   erf((alog(d_cut)-B)/(sqrt(2.0)*alog(sigma)))
!     with
!            A = (pi/6.0)*exp(3.0*alog(dbar)+4.5*(alog(sigma))**2)
!            B = alog (volume median diameter) =
!              = alog(d_bar)+3.0*(alog(sigma))**2
!     See also note (3) for the calculation of Ntot

!     (2) Known V(d_cut) in m3, the mass concentration (ug / m3) of
!         particles with diameter lower than d_cut is:
!            V(d_cut) * rho_particle * 1e9
!            (need to multiply by 1e9 because rho is in kg/m3)

!     (3) Sea salt diagnostics are aerosol number (N) while the
!     diagnostics for the other aerosol species are m.m.r.
!     In those cases one can calculate Ntot with the relationship:
!
!                            rho_air           1
!      N (m-3) = m.m.r. * ------------  * -------------, with
!                         rho_particle     V_particle

!      V_particle = (4*PI/3) * (r_bar**3) * exp (4.5*(ln (sigma))**2),
!      because 2*r_bar is the geometric diameter

!                                          3 * rho_air
!      We finally use:  N (m-3) = m.m.r. * -------------
!                                            denom
!                      (where denom = 3*rho_particle*V_particle)

!     (4) The size distribution of mineral dust is currently not lognormal.
!         Six size bins are used instead, covering
!         3.162e-8 to 1.000e-7 m (bin 1), 1.000e-8 to 3.162e-7 m (bin 2),
!         3.162e-7 to 1.000e-6 m (bin 3), 1.000e-6 to 3.162e-6 m (bin 4),
!         3.162e-6 to 1.000e-5 (bin 5), and 1.000e-5 to 3.162e-5 m (bin 6).
!         Consequently, PM2.5 concentrations are derived from m.m.r. of
!         bins 1 to 3 and 19.3941% of bin 4, and PM10 concentrations are
!         derived from m.m.r. of bins 1 to 4 and 39.8317% of bin 5:
!            m.m.r. * rho_air (kg/m3) * 1e9 --> ug / m3
!         We also calculate a total concentration by using the total
!         m.m.r. across all 6 bins.  

!     (5) All contributions are finally added to get the total
!         PM10 and PM2.5 mass concentrations

! Limitations:
!      This code currently ignores hygroscopic growth of aerosol

!  Calculation of total concentration is also possible, but
!  much simpler since no cut-off diameter is required. This is
!  currently only done for mineral dust.
!
!  Called by Aero_Ctl 

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols

! Code Description:
!   Language:  Fortran 90
!   This code is written to UMDP3 version 8.3 programming standards.
!--------------------------------------------------------------------------

! Model dimensions
INTEGER, INTENT(IN) :: off_x      ! Size of small halo in i
INTEGER, INTENT(IN) :: off_y      ! Size of small halo in j.
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
!dimensions of seasalt array
INTEGER, INTENT(IN) :: salt_dim1
INTEGER, INTENT(IN) :: salt_dim2
INTEGER, INTENT(IN) :: salt_dim3
!
LOGICAL, INTENT(IN) :: l_sulpc_so2            !T if Sulphur Cycle required
LOGICAL, INTENT(IN) :: l_soot                 !T if SOOT modelling required
LOGICAL, INTENT(IN) :: l_biomass              !T if biomass modelling reqd
LOGICAL, INTENT(IN) :: l_ocff                 !T if OCFF modelling required
LOGICAL, INTENT(IN) :: l_use_biogenic         !T if clim. of biog. aeros. used
LOGICAL, INTENT(IN) :: l_dust                 !T if mineral dust used
LOGICAL, INTENT(IN) :: l_nitrate              !T if nitrate modelling required
LOGICAL, INTENT(IN) :: l_use_seasalt_pm       !T if want to include SS in PM

! Data fields IN
! pressure on theta levels
REAL, INTENT(IN) :: p_theta_levels(1-off_x:row_length+off_x,                    &
                                   1-off_y:rows+off_y, model_levels)
! Temp (K) on theta levels
REAL, INTENT(IN) ::  t(row_length,rows,model_levels)

! mmr S in AITKEN and ACCUMULATION modes
REAL, INTENT(IN) :: so4_ait(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
                               model_levels)
REAL, INTENT(IN) :: so4_acc(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
                            model_levels)

! mmr fresh and aged soot
REAL, INTENT(IN) :: soot_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
                                model_levels)
REAL, INTENT(IN) :: soot_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
                                model_levels)

!mmr fresh and aged smoke
REAL, INTENT(IN) ::    bmass_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                                 model_levels)
REAL, INTENT(IN) ::    bmass_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                                 model_levels)

! mmr fresh and aged OCFF
REAL, INTENT(in) ::    ocff_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
                                model_levels)
REAL, INTENT(in) ::    ocff_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
                                model_levels)
! mmr biogenics
REAL, INTENT(in) ::    biogenic(row_length, rows, model_levels)

! Sea-salt aerosol number densities
REAL, INTENT(in) ::    sea_salt_film(salt_dim1,salt_dim2,salt_dim3)
REAL, INTENT(in) ::    sea_salt_jet (salt_dim1,salt_dim2,salt_dim3)

! mmr dust in 5 divisions
REAL, INTENT(in) ::    dust_1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
                              model_levels)
REAL, INTENT(in) ::    dust_2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
                              model_levels)
REAL, INTENT(in) ::    dust_3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
                              model_levels)
REAL, INTENT(in) ::    dust_4(1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
                              model_levels)
REAL, INTENT(in) ::    dust_5(1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
                              model_levels)
REAL, INTENT(in) ::    dust_6(1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
                              model_levels)

! mmr Nitrate in ACC
REAL, INTENT(in) ::    nitr_acc(1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
                                model_levels)

! Output arguments
!PM10  (micrograms m-3)
REAL, INTENT(OUT) :: pm10       (row_length, rows, model_levels)
!PM2.5 (micrograms m-3)
REAL, INTENT(OUT) :: pm2p5      (row_length, rows, model_levels)
!PM concentrations due to sulphate
REAL, INTENT(OUT) :: pm10_so4   (row_length, rows, model_levels)
REAL, INTENT(OUT) :: pm2p5_so4  (row_length, rows, model_levels)
!PM concs. due to black carbon
REAL, INTENT(OUT) :: pm10_bc    (row_length, rows, model_levels)
REAL, INTENT(OUT) :: pm2p5_bc   (row_length, rows, model_levels)
!PM concs. due to biomass burning aerosol
REAL, INTENT(OUT) :: pm10_bb    (row_length, rows, model_levels)
REAL, INTENT(OUT) :: pm2p5_bb   (row_length, rows, model_levels)
!PM concs. due to OCFF
REAL, INTENT(OUT) :: pm10_ocff  (row_length, rows, model_levels)
REAL, INTENT(OUT) :: pm2p5_ocff (row_length, rows, model_levels)
!PM concs. due to SOA
REAL, INTENT(OUT) :: pm10_soa   (row_length, rows, model_levels)
REAL, INTENT(OUT) :: pm2p5_soa  (row_length, rows, model_levels)
!PM concs. due to sea-salt aerosol
REAL, INTENT(OUT) :: pm10_ss   (row_length, rows, model_levels)
REAL, INTENT(OUT) :: pm2p5_ss  (row_length, rows, model_levels)
! concs. due to mineral dust
REAL, INTENT(OUT) :: conc_dust (row_length, rows, model_levels)
REAL, INTENT(OUT) :: pm10_dust (row_length, rows, model_levels)
REAL, INTENT(OUT) :: pm2p5_dust(row_length, rows, model_levels)
!PM concs. due to nitrate
REAL, INTENT(OUT) :: pm10_nitr (row_length, rows, model_levels)
REAL, INTENT(OUT) :: pm2p5_nitr(row_length, rows, model_levels)


! Local parameters for PM10 and PM2.5 calculations
! Parameters related to the volume distribution were
! calculated offline by the program calc_pm_params.f90
! All densities are in kg/m3
! Parameters for sulphate
REAL, PARAMETER :: rho_su          = 1769.0     ! density 
!What is advected is sulphur. Therefore,
!necessary to convert m.m.r. of sulphur to ammonium
!sulphate by multiplying by ratio of molecular weights:
!(NH4)2SO4 / S = 132 / 32 =  4.125
REAL, PARAMETER :: s_conv_fac      = 4.125
REAL, PARAMETER :: denom_su_ait    = 0.832154093803605184E-20
REAL, PARAMETER :: a_param_su_ait  = 0.156803107933598664E-23
REAL, PARAMETER :: erf_su_ait_d1   = 1.00000000000000000
REAL, PARAMETER :: erf_su_ait_d2   = 1.00000000000000000
REAL, PARAMETER :: denom_su_acc    = 0.317222975403968212E-16
REAL, PARAMETER :: a_param_su_acc  = 0.597744442065136312E-20
REAL, PARAMETER :: erf_su_acc_d1   = 1.00000000000000000
REAL, PARAMETER :: erf_su_acc_d2   = 0.999999999970596409
! Parameters for Black Carbon (soot)
REAL, PARAMETER :: rho_bc          = 1900.0     ! density 
REAL, PARAMETER :: denom_bc_fr     = 0.132771498349419138E-16
REAL, PARAMETER :: a_param_bc_fr   = 0.232932453244595951E-20
REAL, PARAMETER :: erf_bc_fr_d1    = 0.999998972736938496
REAL, PARAMETER :: erf_bc_fr_d2    = 0.996102525359877422
REAL, PARAMETER :: denom_bc_ag     = 0.132771498349419138E-16
REAL, PARAMETER :: a_param_bc_ag   = 0.232932453244595951E-20
REAL, PARAMETER :: erf_bc_ag_d1    = 0.999998972736938496
REAL, PARAMETER :: erf_bc_ag_d2    = 0.996102525359877422
! Parameters for biomass burning aerosol
REAL, PARAMETER :: rho_bb          = 1350.0     ! density 
REAL, PARAMETER :: denom_bb_fr     = 0.231243545609857872E-16
REAL, PARAMETER :: a_param_bb_fr   = 0.570971717555206423E-20
REAL, PARAMETER :: erf_bb_fr_d1    = 1.00000000000000000
REAL, PARAMETER :: erf_bb_fr_d2    = 1.00000000000000000
REAL, PARAMETER :: denom_bb_ag     = 0.399588846813834411E-16
REAL, PARAMETER :: a_param_bb_ag   = 0.986639127935394024E-20
REAL, PARAMETER :: erf_bb_ag_d1    = 1.00000000000000000
REAL, PARAMETER :: erf_bb_ag_d2    = 0.999999999999999667
! Parameters for OCFF
REAL, PARAMETER :: rho_ocff        = 1350.0     ! density 
REAL, PARAMETER :: denom_ocff_fr   = 0.231243545609857872E-16
REAL, PARAMETER :: a_param_ocff_fr = 0.570971717555206423E-20
REAL, PARAMETER :: erf_ocff_fr_d1  = 1.00000000000000000
REAL, PARAMETER :: erf_ocff_fr_d2  = 1.00000000000000000
REAL, PARAMETER :: denom_ocff_ag   = 0.399588846813834411E-16
REAL, PARAMETER :: a_param_ocff_ag = 0.986639127935394024E-20
REAL, PARAMETER :: erf_ocff_ag_d1  = 1.00000000000000000
REAL, PARAMETER :: erf_ocff_ag_d2  = 0.999999999999999667
! Parameters for biogenic aerosol (SOA)
REAL, PARAMETER :: rho_soa         = 1300.0     ! density 
REAL, PARAMETER :: denom_soa       = 0.293507300762206618E-16
REAL, PARAMETER :: a_param_soa     = 0.752582822467193671E-20
REAL, PARAMETER :: erf_soa_d1      = 1.00000000000000000
REAL, PARAMETER :: erf_soa_d2      = 0.999999724269635237
! Parameters for sea-salt
REAL, PARAMETER :: rho_ss          = 2165.0     ! density 
REAL, PARAMETER :: a_param_ss_fi   = 0.267438839195005737E-19
REAL, PARAMETER :: erf_ss_fi_d1    = 0.999969448922318982
REAL, PARAMETER :: erf_ss_fi_d2    = 0.955514879308947629
REAL, PARAMETER :: a_param_ss_je   = 0.363956958194679524E-16
REAL, PARAMETER :: erf_ss_je_d1    = 0.191596825012604527
REAL, PARAMETER :: erf_ss_je_d2    = -0.921169668112190365
! Parameters for nitrate
REAL, PARAMETER :: rho_ni          = 1725.0     ! density 
!What is advected is nitrogen. Therefore,
!necessary to convert m.m.r. of N to ammonium
!nitrate by multiplying by ratio of molecular weights:
!NH4NO3 / N =  80 / 14 =  5.714
REAL, PARAMETER :: n_conv_fac      = 5.714
REAL, PARAMETER :: denom_ni_acc    = 0.309332748768708366E-16
REAL, PARAMETER :: a_param_ni_acc  = 0.597744442065136312E-20
REAL, PARAMETER :: erf_ni_acc_d1   = 1.00000000000000000
REAL, PARAMETER :: erf_ni_acc_d2   = 0.999999999970596409
! Parameters for six-bin dust
REAL, PARAMETER :: dust5_pm10_frac  = 0.398317 ! Fraction of bin 5 <10 microns
REAL, PARAMETER :: dust4_pm2p5_frac = 0.193941 ! Fraction of bin 4 <2.5 microns
! Parameters for two-bin dust
REAL, PARAMETER :: dust1_pm2p5_frac_twobin = 0.405559 
                                               ! Fraction of bin 1 < 2.5 microns
REAL, PARAMETER :: dust2_pm10_frac_twobin  = 0.507086
                                               ! Fraction of bin 2 < 10 microns

REAL :: rho_air    ! air density
REAL :: n_tot      ! total nr. aerosol conc. 

INTEGER :: i, j         ! loop counters (horizontal field indices)
INTEGER :: k            ! loop counter  (vertical index)
INTEGER :: icode        ! passed to UM routine ereport
CHARACTER(len=80)             :: cmessage
CHARACTER (LEN=* ), PARAMETER :: routinename='PM10_PM2P5'

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(routinename,zhook_in,zhook_handle)

! Add contributions from the different aerosol species and modes
! to calculate pm10 & pm2.5 concentrations

!initialise pm10 & pm2.5 concentrations
pm10       (:,:,:) = 0.0
pm2p5      (:,:,:) = 0.0
pm10_so4   (:,:,:) = 0.0
pm2p5_so4  (:,:,:) = 0.0
pm10_bc    (:,:,:) = 0.0
pm2p5_bc   (:,:,:) = 0.0
pm10_bb    (:,:,:) = 0.0
pm2p5_bb   (:,:,:) = 0.0
pm10_ocff  (:,:,:) = 0.0
pm2p5_ocff (:,:,:) = 0.0
pm10_soa   (:,:,:) = 0.0
pm2p5_soa  (:,:,:) = 0.0
pm10_ss    (:,:,:) = 0.0
pm2p5_ss   (:,:,:) = 0.0
pm10_dust  (:,:,:) = 0.0
conc_dust  (:,:,:) = 0.0
pm2p5_dust (:,:,:) = 0.0
pm10_nitr  (:,:,:) = 0.0
pm2p5_nitr (:,:,:) = 0.0

DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length

      rho_air = p_theta_levels(i, j, k)/(r * t(i, j, k))

      IF (l_sulpc_so2) THEN
        ! includes conversion of m.m.r. of sulphur to ammonium
        n_tot = 3.0*(so4_ait(i, j, k)*s_conv_fac)*rho_air / denom_su_ait

        pm10_so4(i, j, k) =                                                    &
          ((a_param_su_ait * n_tot / 2.0)                                      &
          +(a_param_su_ait * n_tot / 2.0)                                      &
          * erf_su_ait_d1) * rho_su * 1.0e9

        pm2p5_so4(i, j, k) =                                                   &
          ((a_param_su_ait * n_tot/2.0)                                        &
          +(a_param_su_ait * n_tot/2.0)                                        &
          * erf_su_ait_d2) * rho_su * 1.0e9

        n_tot = 3.0*(so4_acc(i, j, k)*s_conv_fac)*rho_air / denom_su_acc

        pm10_so4(i, j, k) = pm10_so4(i, j, k) +                                &
          ((a_param_su_acc * n_tot / 2.0)                                      &
          +(a_param_su_acc * n_tot / 2.0)                                      &
          * erf_su_acc_d1) * rho_su * 1.0e9

        pm2p5_so4(i, j, k) = pm2p5_so4(i, j, k) +                              &
          ((a_param_su_acc * n_tot/2.0)                                        &
          +(a_param_su_acc * n_tot/2.0)                                        &
          * erf_su_acc_d2) * rho_su * 1.0e9

        pm10  (i, j, k) = pm10_so4 (i, j, k)
        pm2p5 (i, j, k) = pm2p5_so4(i, j, k)
      END IF


      IF (l_soot) THEN
        n_tot = 3.0 * soot_new(i, j, k) * rho_air / denom_bc_fr

        pm10_bc(i, j, k) =                                                     &
          ((a_param_bc_fr * n_tot / 2.0)                                       &
          +(a_param_bc_fr * n_tot / 2.0)                                       &
          * erf_bc_fr_d1) * rho_bc * 1.0e9

        pm2p5_bc(i, j, k) =                                                    &
          ((a_param_bc_fr * n_tot / 2.0)                                       &
          +(a_param_bc_fr * n_tot / 2.0)                                       &
          * erf_bc_fr_d2) * rho_bc * 1.0e9

        n_tot = 3.0 * soot_agd (i, j, k) * rho_air / denom_bc_ag

        pm10_bc(i, j, k) =  pm10_bc(i, j, k) +                                 &
          ((a_param_bc_ag * n_tot / 2.0)                                       &
          +(a_param_bc_ag * n_tot / 2.0)                                       &
          * erf_bc_ag_d1) * rho_bc * 1.0e9

        pm2p5_bc(i, j, k) = pm2p5_bc(i, j, k) +                                &
          ((a_param_bc_ag * n_tot / 2.0)                                       &
          +(a_param_bc_ag * n_tot / 2.0)                                       &
          * erf_bc_ag_d2) * rho_bc * 1.0e9

        pm10  (i, j, k) = pm10  (i, j, k) + pm10_bc (i, j, k)
        pm2p5 (i, j, k) = pm2p5 (i, j, k) + pm2p5_bc(i, j, k)
      END IF


      IF (l_biomass) THEN
        n_tot = 3.0 * bmass_new(i, j, k) * rho_air / denom_bb_fr

        pm10_bb(i, j, k) =                                                     &
          ((a_param_bb_fr * n_tot / 2.0)                                       &
          +(a_param_bb_fr * n_tot / 2.0)                                       &
          * erf_bb_fr_d1) * rho_bb * 1.0e9

        pm2p5_bb(i, j, k) =                                                    &
          ((a_param_bb_fr * n_tot / 2.0)                                       &
          +(a_param_bb_fr * n_tot / 2.0)                                       &
          * erf_bb_fr_d2) * rho_bb * 1.0e9

        n_tot = 3.0 * bmass_agd(i, j, k) * rho_air / denom_bb_ag

        pm10_bb(i, j, k) =  pm10_bb(i, j, k) +                                 &
          ((a_param_bb_ag * n_tot / 2.0)                                       &
          +(a_param_bb_ag * n_tot / 2.0)                                       &
          * erf_bb_ag_d1) * rho_bb * 1.0e9

        pm2p5_bb(i, j, k) = pm2p5_bb(i, j, k) +                                &
          ((a_param_bb_ag * n_tot / 2.0)                                       &
          +(a_param_bb_ag * n_tot / 2.0)                                       &
          * erf_bb_ag_d2) * rho_bb * 1.0e9

        pm10  (i, j, k) = pm10  (i, j, k) + pm10_bb (i, j, k)
        pm2p5 (i, j, k) = pm2p5 (i, j, k) + pm2p5_bb(i, j, k)
      END IF


      IF (l_ocff) THEN
        n_tot = 3.0 * ocff_new(i, j, k) * rho_air / denom_ocff_fr

        pm10_ocff(i, j, k) =                                                   &
          ((a_param_ocff_fr * n_tot / 2.0)                                     &
          +(a_param_ocff_fr * n_tot / 2.0)                                     &
          * erf_ocff_fr_d1) * rho_ocff * 1.0e9

        pm2p5_ocff(i, j, k) =                                                  &
          ((a_param_ocff_fr * n_tot/2.0)                                       &
          +(a_param_ocff_fr * n_tot/2.0)                                       &
          * erf_ocff_fr_d2) * rho_ocff * 1.0e9

        n_tot = 3.0 * ocff_agd(i, j, k) * rho_air / denom_ocff_ag

        pm10_ocff(i, j, k) =  pm10_ocff(i, j, k) +                             &
          ((a_param_ocff_ag * n_tot / 2.0)                                     &
          +(a_param_ocff_ag * n_tot / 2.0)                                     &
          * erf_ocff_ag_d1) * rho_ocff * 1.0e9

        pm2p5_ocff(i, j, k) = pm2p5_ocff(i, j, k) +                            &
          ((a_param_ocff_ag * n_tot/2.0)                                       &
          +(a_param_ocff_ag * n_tot/2.0)                                       &
          * erf_ocff_ag_d2) * rho_ocff * 1.0e9

        pm10  (i, j, k) = pm10  (i, j, k) + pm10_ocff (i, j, k)
        pm2p5 (i, j, k) = pm2p5 (i, j, k) + pm2p5_ocff(i, j, k)
      END IF

      IF (l_use_biogenic) THEN
        n_tot = 3.0 * biogenic(i, j, k) * rho_air / denom_soa

        pm10_soa(i, j, k) =                                                    &
          ((a_param_soa * n_tot / 2.0)                                         &
          +(a_param_soa * n_tot / 2.0)                                         &
          * erf_soa_d1) * rho_soa * 1.0e9

        pm2p5_soa(i, j, k) =                                                   &
          ((a_param_soa * n_tot / 2.0)                                         &
          +(a_param_soa * n_tot / 2.0)                                         &
          * erf_soa_d2) * rho_soa * 1.0e9

        pm10  (i, j, k) = pm10  (i, j, k) + pm10_soa (i, j, k)
        pm2p5 (i, j, k) = pm2p5 (i, j, k) + pm2p5_soa(i, j, k)
      END IF

      IF (l_use_seasalt_pm) THEN
        n_tot = sea_salt_film(i, j, k)

        pm10_ss(i, j, k) =                                                     &
          ((a_param_ss_fi * n_tot / 2.0)                                       &
          +(a_param_ss_fi * n_tot / 2.0)                                       &
          * erf_ss_fi_d1) * rho_ss * 1.0e9

        pm2p5_ss(i, j, k) =                                                    &
          ((a_param_ss_fi * n_tot / 2.0)                                       &
          +(a_param_ss_fi * n_tot / 2.0)                                       &
          * erf_ss_fi_d2) * rho_ss * 1.0e9

        n_tot = sea_salt_jet(i, j, k)

        pm10_ss(i, j, k) =  pm10_ss(i, j, k) +                                 &
          ((a_param_ss_je * n_tot / 2.0)                                       &
          +(a_param_ss_je * n_tot / 2.0)                                       &
          * erf_ss_je_d1) * rho_ss * 1.0e9

        pm2p5_ss(i, j, k) = pm2p5_ss(i, j, k) +                                &
          ((a_param_ss_je * n_tot / 2.0)                                       &
          +(a_param_ss_je * n_tot / 2.0)                                       &
          * erf_ss_je_d2) * rho_ss * 1.0e9

        pm10  (i, j, k) = pm10  (i, j, k) + pm10_ss (i, j, k)
        pm2p5 (i, j, k) = pm2p5 (i, j, k) + pm2p5_ss(i, j, k)
      END IF

      IF (l_dust) THEN
        IF (l_twobin_dust) THEN   ! two-bin dust      
          conc_dust(i, j, k) =                                                 &
             ((dust_1(i, j, k)                                                 &
             + dust_2(i, j, k))                                                &
             * rho_air * 1.0e9)
             
          pm10_dust(i, j, k) =                                                 &
             ((dust_1(i, j, k)                                                 &          
             + (dust2_pm10_frac_twobin * dust_2(i, j, k)))                     &
             * rho_air * 1.0e9)
             
          pm2p5_dust(i,j, k) =                                                 &
             ((dust1_pm2p5_frac_twobin  * dust_1(i, j, k))                     &
             * rho_air * 1.0e9)
        ELSE ! Six-bin dust
          conc_dust(i, j, k) =                                                 &
             ((DUST_1(i, j, k)                                                 &
             + DUST_2(i, j, k)                                                 &
             + DUST_3(i, j, k)                                                 &
             + DUST_4(i, j, k)                                                 &
             + DUST_5(i, j, k)                                                 &
             + DUST_6(i, j, k))                                                &
             * rho_air * 1.0e9)

         pm10_dust(i, j, k) =                                                  &
            ((dust_1(i, j, k)                                                  &
            + dust_2(i, j, k)                                                  &
            + dust_3(i, j, k)                                                  &
            + dust_4(i, j, k)                                                  &
            + (dust5_pm10_frac * dust_5(i, j, k)))                             &
            * rho_air * 1.0e9)

          pm2p5_dust(i, j, k) =                                                &
            ((dust_1(i, j, k)                                                  &
            + dust_2(i, j, k)                                                  &
            + dust_3(i, j, k)                                                  &
            + (dust4_pm2p5_frac * dust_4(i, j, k)))                            &
            * rho_air * 1.0e9)

          pm10  (i, j, k) = pm10  (i, j, k) + pm10_dust (i, j, k)
          pm2p5 (i, j, k) = pm2p5 (i, j, k) + pm2p5_dust(i, j, k)
        END IF ! l_twobin_dust
      END IF

      IF (l_nitrate) THEN
        !includes  conversion of  m.m.r. of N to ammonium nitrate 
        n_tot = 3.0*(nitr_acc(i, j, k)*n_conv_fac)*rho_air / denom_ni_acc

        pm10_nitr(i, j, k) =                                                   &
          ((a_param_ni_acc * n_tot / 2.0)                                      &
          +(a_param_ni_acc * n_tot / 2.0)                                      &
          * erf_ni_acc_d1) * rho_ni * 1.0e9

        pm2p5_nitr(i, j, k) =                                                  &
          ((a_param_ni_acc * n_tot/2.0)                                        &
          +(a_param_ni_acc * n_tot/2.0)                                        &
          * erf_ni_acc_d2) * rho_ni * 1.0e9

        pm10  (i, j, k) = pm10  (i, j, k) + pm10_nitr (i, j, k)
        pm2p5 (i, j, k) = pm2p5 (i, j, k) + pm2p5_nitr(i, j, k)
      END IF

    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(routinename,zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_pm_diags

END MODULE calc_pm_diags_mod

