! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE conv_diag_0a_mod

IMPLICIT NONE

CONTAINS

!  Diagnosis of convective occurrence and type

SUBROUTINE conv_diag_0a(                                          &

! IN values defining field dimensions and subset to be processed :
  row_length, rows                                                &

! IN values defining vertical grid of model atmosphere :
, bl_levels, model_levels, wet_model_levels                       &
, land_points                                                     &
, p, p_theta_lev,exner_rho                                        &
, rho_only, rho_theta, z_full, z_half                             &

! IN Model switches
, l_mixing_ratio,l_ctile, l_extra_call                            &
, no_cumulus                                                      &

! IN Cloud data :
, qcf, qcl, cloud_fraction                                        &

! IN everything not covered so far :
, pstar, q, theta, exner_theta_levels, u_p, v_p, u_0_p, v_0_p     &
, tstar_land, tstar_sea, tstar_sice, z0msea                       &
, l_flux_bc, flux_e, flux_h, l_spec_z0, z0m_scm, z0h_scm          &
, tstar, land_mask, flandg, ice_fract, timestep                   &
, w_copy, w_max                                                   &
, deep_flag, past_precip, past_conv_ht                            &

! SCM Diagnostics (dummy values in full UM)
, nscmdpkgs,l_scmdiags                                            &

! OUT data required elsewhere in UM system :
, zh,zhpar,dzh,z_lcl,z_lcl_uv,delthvu,ql_ad,ntml,ntpar,nlcl       &
, cumulus, l_shallow,l_congestus, l_congestus2, conv_type, cin    &
, cape, wstar, wthvs, entrain_coef, qsat_lcl                      &
, error                                                           &
 )

!-----------------------------------------------------------------------
! Purpose:
!   Diagnosis of convective occurrence - version for use with no
!   convection scheme.
!
!   Called by Atmos_physics2
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3  programming standards v8.3
!
!-----------------------------------------------------------------------

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                    &
  pdims, pdims_s, tdims_s, tdims, qdims, wdims

! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels             &  ! Radii on theta levels (m)
 ,r_rho_levels                  ! Radii on rho levels (m)

USE atmos_constants_mod, ONLY:                                      &
  cp, kappa, r, repsilon, c_virtual

USE water_constants_mod, ONLY: lc, lf, tm

USE earth_constants_mod, ONLY: g, earth_radius

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


! subroutines
USE conv_surf_flux_mod, ONLY: conv_surf_flux

USE Field_Types
IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

! Arguments with intent IN:


INTEGER, INTENT(IN) ::  &
  row_length            & ! Local number of points on a row
 ,rows                    ! Local number of rows in a theta field

! (b) Defining vertical grid of model atmosphere.

INTEGER, INTENT(IN) :: &
  bl_levels            & ! Max. no. of "boundary" levels allowed.
 ,model_levels         & ! number of model levels
 ,wet_model_levels     & ! number of wet model levels
 ,land_points            ! number of land points

REAL, INTENT(IN) ::                          &
  p(pdims_s%i_start:pdims_s%i_end,           & ! pressure  on rho levels (Pa)
    pdims_s%j_start:pdims_s%j_end,           &
    pdims_s%k_start:pdims_s%k_end)           &
 ,p_theta_lev(tdims%i_start:tdims%i_end,     & ! P on theta lev (Pa)
              tdims%j_start:tdims%j_end,     &
                          1:tdims%k_end)     &
 ,exner_rho(pdims_s%i_start:pdims_s%i_end,   & ! Exner on rho level
            pdims_s%j_start:pdims_s%j_end,   & !
            pdims_s%k_start:pdims_s%k_end)   &
 ,rho_only(row_length,rows,1:pdims%k_end)    & ! density (kg/m3)
 ,rho_theta(row_length,rows,1:tdims%k_end-1) & ! rho th lev (kg/m3)
 ,z_full(row_length,rows,1:tdims%k_end)      & ! height th lev (m)
 ,z_half(row_length,rows,1:pdims%k_end)        ! height rho lev (m)

LOGICAL,INTENT(IN) ::   &
  l_mixing_ratio        & ! true moisture input as mixing ratios
                          ! false moisture input as specific humidity
 ,l_ctile               & ! true if coastal tiling
 ,l_flux_bc             & ! true if SCM using specified surface fluxes
 ,l_spec_z0             & ! true if roughness length has been specified
 ,l_extra_call            ! true this is an additional call to conv_diag
                          ! within a timestep

LOGICAL,INTENT(IN) ::   &
  no_cumulus(row_length,rows)   ! Points overruled by BL

! (c) Cloud data.

REAL, INTENT(IN) ::                          &
  qcf(qdims%i_start:qdims%i_end,             & ! Cloud ice (kg per kg air)
      qdims%j_start:qdims%j_end,             &
                  1:qdims%k_end)             &
 ,qcl(qdims%i_start:qdims%i_end,             & !Cloud liquid water (kg/kg air)
      qdims%j_start:qdims%j_end,             &
                  1:qdims%k_end)             &
 ,cloud_fraction(qdims%i_start:qdims%i_end,  & !  Cloud fraction
                 qdims%j_start:qdims%j_end,  &
                             1:qdims%k_end)

! (d) Atmospheric + any other data not covered so far, incl control.

REAL, INTENT(IN) ::                             &
  pstar(row_length, rows)                       & ! Surface pressure (Pa)
 ,q(qdims%i_start:qdims%i_end,                  & ! water vapour (kg/kg)
    qdims%j_start:qdims%j_end,                  &
                1:qdims%k_end)                  &
 ,theta(tdims%i_start:tdims%i_end,              & ! Theta (Kelvin)
        tdims%j_start:tdims%j_end,              &
                    1:tdims%k_end)              &
 ,exner_theta_levels(tdims%i_start:tdims%i_end, & ! exner pressure theta lev
                     tdims%j_start:tdims%j_end, & !  (Pa)
                                 1:tdims%k_end)

REAL, INTENT(IN) ::            &
  u_p(row_length, rows)        & ! U(1) on P-grid.
 ,v_p(row_length, rows)        & ! V(1) on P-grid.
 ,u_0_p(row_length,rows)       & ! W'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,v_0_p(row_length,rows)       & ! S'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,flux_e(row_length,rows)      & ! Specified surface
                                 !    latent heat flux (W/m^2)
 ,flux_h(row_length,rows)      & ! Specified surface
                                 !    sensible heat fluxes (in W/m2)
 ,z0msea(row_length,rows)      & ! Sea roughness length for momentum (m)
 ,z0m_scm(row_length,rows)     & ! Namelist input z0m (if >0)
 ,z0h_scm(row_length,rows)       ! Namelist input z0h (if >0)

REAL, INTENT(IN) :: &
  tstar_land(row_length, rows) & ! Surface T on land
 ,tstar_sea(row_length, rows)  & ! Surface T on sea
 ,tstar_sice(row_length, rows)   ! Surface T on sea-ice 

REAL, INTENT(INOUT) ::    &
  tstar(row_length,rows)    ! Surface temperature (= top soil layer temperature)
                            ! (K). NOTE only inout for SCM

LOGICAL,INTENT(IN) ::         &
  land_mask(row_length, rows)   ! T If land, F Elsewhere.

REAL, INTENT(IN) ::            &
  flandg(pdims_s%i_start:pdims_s%i_end, & ! Land fraction of gridbox
         pdims_s%j_start:pdims_s%j_end) & ! on all points
 ,ice_fract(row_length,rows)              ! fraction of sea that has ice

REAL, INTENT(IN) ::                  &
  timestep                           & ! timestep (seconds).
 ,w_copy(wdims%i_start:wdims%i_end,  & ! vertical velocity W (m/s)
         wdims%j_start:wdims%j_end,  &
                     0:wdims%k_end)  & ! Not exact match to module values
 ,w_max(row_length,rows)               ! Column max vertical velocity (m/s)

REAL, INTENT(IN) ::              &
  deep_flag(row_length,rows)     & ! 0-1.0, 1 if deep last time step
 ,past_precip(row_length,rows)   & ! convective precip rate last step
                                   ! or a decayed value.
 ,past_conv_ht(row_length,rows)    ! convective height (m)


! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::  &
  nscmdpkgs               ! No of diagnostics packages

LOGICAL,INTENT(IN) ::   &
  l_scmdiags(nscmdpkgs)   ! Logicals for diagnostics packages

REAL, INTENT(INOUT) ::         &
  zh(row_length,rows)            ! Height above surface of top
                                 !  of boundary layer (metres).

REAL, INTENT(OUT) ::             &
  zhpar(row_length,rows)         & ! Height of max parcel ascent (m)
 ,dzh(row_length,rows)           & ! Height of inversion top (m)
 ,z_lcl(row_length,rows)         & ! Height of lifting condensation
                                   ! level  (m)
 ,z_lcl_uv(row_length,rows)      & ! Height of lifting condensation
                                   ! level on uv grid (m) (Not set 0A scheme)
 ,delthvu(row_length,rows)       & ! Integral of undilute parcel buoyancy
                                   ! over convective cloud layer
                                   ! (for convection scheme)
 ,ql_ad(row_length,rows)         & ! adiabatic liquid water content at
                                   ! inversion or cloud top (kg/kg) (Not 0A)
 ,entrain_coef(row_length,rows)  & ! Entrainment coefficient        (Not 0A)
 ,qsat_lcl(row_length,rows)      & ! qsat at cloud base (kg/kg)     (Not 0A)
 ,cape(row_length, rows)         & ! CAPE from parcel ascent (m2/s2)
 ,cin(row_length, rows)            ! CIN from parcel ascent (m2/s2)

INTEGER, INTENT(OUT) ::      &
  ntml(row_length,rows)      & ! Number of model levels in the
                               ! turbulently mixed layer.
 ,ntpar(row_length,rows)     & ! Max levels for parcel ascent
 ,nlcl(row_length,rows)        ! No. of model layers below the
                               ! lifting condensation level.

LOGICAL,INTENT(OUT) ::         &
  cumulus(row_length,rows)     & ! Logical indicator for convection
 ,l_shallow(row_length,rows)   & ! Logical indicator for shallow Cu   (Not 0A)
 ,l_congestus(row_length,rows) & ! Logical indicator for congestus Cu (Not 0A)
 ,l_congestus2(row_length,rows)  ! Logical ind 2 for congestus Cu     (Not 0A)

! Required for extra call to conv_diag as values from original BL call
! may return zero values based on initial timestep profiles. (Not 0A)
REAL, INTENT(INOUT) ::    &
  wstar(row_length,rows)  & ! Convective sub-cloud velocity scale (m/s)
 ,wthvs(row_length,rows)    ! surface flux of wthv (K m/s)
INTEGER, INTENT(INOUT) :: &
  error                     ! 0 - no error in this routine


! Redundant arguments
! -------------------
! Convective type array ::
INTEGER, INTENT(IN) :: &
  conv_type(row_length, rows)

!-----------------------------------------------------------------------
! Variables defined locally
! NOTE - compression arrays are defined for row_length*rows though only
!        nunstable points are used. Memory could be reduced by using
!        ALLOCATE. Less memory may be required if we did not use compression
!        in the first place as some arrays are copies.
!        Compression is being done as it speeds up some of the code.
!-----------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER ::  routinename = 'conv_diag_0a'

INTEGER ::   &
  i,j,ii     & ! Local Loop counter (horizontal field index).
 ,k          & ! Locla Loop counter (vertical level index).
 ,nunstable    ! total number of unstable points

INTEGER :: &
  index_i(row_length*rows)   & ! column number of unstable points
 ,index_j(row_length*rows)     ! row number of unstable points

! uncompressed arrays - all points

REAL ::                         &
  fb_surf(row_length, rows)       ! Change in theta_v from surface
                                  ! to layer 1 (note diff from BL)

REAL  ::                        &
  tv1_sd( row_length*rows)      & ! Approx to standard dev of level
                                  ! 1 virtual temperature (K).
 ,bl_vscale2(row_length*rows)   & ! Velocity scale squared for 
                                  ! boundary layer eddies (m2/s2)
 ,flandg_c(row_length*rows)       ! fraction of land compressed 



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Mixing ratio, r,  versus specific humidity, q
!
! In most cases the expression to first order are the same
!
!  Tl = T - (lc/cp)qcl - [(lc+lf)/cp]qcf
!  Tl = T - (lc/cp)rcl - [(lc+lf)/cp]rcf  - equally correct definition
!
! thetav = theta(1+cvq)         accurate
!        = theta(1+r/repsilon)/(1+r) ~ theta(1+cvr) approximate
!
! svl = (Tl+gz/cp)*(1+(1/repsilon-1)qt)
!     ~ (Tl+gz/cp)*(1+(1/repsilon-1)rt)
!
! dqsat/dT = repsilon*Lc*qsat/(R*T*T)
! drsat/dT = repsilon*Lc*rsat/(R*T*T)  equally approximate
!
! Only altering the expression for vapour pressure
!
!  e = qp/repsilon       - approximation
!  e = rp/(repsilon+r)   - accurate
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CONV_DIAG_0A',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

error = 0

!-----------------------------------------------------------------------
! 1.1 Verify grid/subset definitions.
!-----------------------------------------------------------------------

IF ( bl_levels <  1 .OR. rows <  1 .OR. tdims%k_end <  1 .OR.    &
     qdims%k_end <  1 ) THEN
  error = 1
  GO TO 9999

END IF

!-----------------------------------------------------------------------
! 1.2 initialisation of output arrays
!-----------------------------------------------------------------------
DO j=1,rows
  DO i=1,row_length

    cumulus(i,j)     = .FALSE.
    l_shallow(i,j)   = .FALSE.

    ! set false as not used if no convection

    l_congestus(i,j) = .FALSE.
    l_congestus2(i,j) = .FALSE.

    ! Requried by Boundary layer scheme

    ntml(i,j) = 1
    nlcl(i,j) = 1

    ! Only used in this routine as no convection call and no convection
    ! and no call to stash for section 5.

    delthvu(i,j) = 0.0
    ql_ad(i,j)    = 0.0
    cape(i,j)     = 0.0
    cin(i,j)      = 0.0

    ! required by Boundary layer scheme

    z_lcl(i,j) = z_half(i,j,nlcl(i,j)+1)


    ! set to zero as not used

    z_lcl_uv(i,j) = 0.0
    zhpar(i,j)    = 0.0          ! set to zero 
    dzh(i,j)      =-9.9E9        ! initialise to large and negative 
    entrain_coef(i,j) = -99.0    ! indicates not set
    ntpar(i,j) = 1

  END DO
END DO
!-----------------------------------------------------------------------
! Calculation of surface buoyancy flux and level one standard deviation
! of virtual temperature.
!-----------------------------------------------------------------------

CALL conv_surf_flux(                                              &
  row_length, rows                                                &
, model_levels, wet_model_levels, land_points                     &
, l_mixing_ratio,l_ctile, l_flux_bc, l_spec_z0                    &
, land_mask                                                       &
, pstar, tstar_land, tstar_sea, tstar_sice, zh, flandg            &
, ice_fract, u_p, v_p, u_0_p, v_0_p                               &
, flux_e, flux_h,  z0msea, z0m_scm, z0h_scm                       &
, z_full, q, theta, exner_theta_levels                            &
, tstar, fb_surf, tv1_sd, bl_vscale2 )


!-----------------------------------------------------------------------
! 2.0 Decide on unstable points  ( fb_surf > 0.0 )
!     Only work on these points for the rest of the calculations.
!-----------------------------------------------------------------------

nunstable = 0           ! total number of unstable points
DO j=1,rows
  DO i=1,row_length
    IF ( fb_surf(i,j)  >   0.0 ) THEN
      nunstable = nunstable + 1
      index_i(nunstable) = i
      index_j(nunstable) = j
    END IF
  END DO ! I
END DO ! J

! land fraction on just unstable points
IF (nunstable > 0) THEN 
  DO ii=1,nunstable
    i = index_i(ii)   
    j = index_j(ii)   
    flandg_c(ii) = flandg(i,j)
  END DO
!-----------------------------------------------------------------------
! Work on just unstable points to diagnose whether convective.
!-----------------------------------------------------------------------
! DEPENDS ON: conv_diag_comp_0a
  CALL conv_diag_comp_0a(                                               &
          row_length, rows                                              &
        , bl_levels, tdims%k_end, qdims%k_end, nunstable                &
        , index_i, index_j                                              &
        , l_mixing_ratio                                                &
        , p, P_theta_lev, exner_rho                                     &
        , rho_only, rho_theta, z_full, z_half                           &
        , qcf, qcl, cloud_fraction                                      &
        , pstar, q, theta, exner_theta_levels                           &
        , land_mask, flandg_c, timestep                                 &
        , w_copy ,w_max, tv1_sd, bl_vscale2                             &
        , deep_flag, past_precip, past_conv_ht                          &
        , nSCMDpkgs, L_SCMDiags                                         &
        , ntml, ntpar, nlcl                                             &
        , cumulus, L_shallow, l_congestus, l_congestus2, conv_type      &
        , zh,zhpar,dzh,z_lcl,z_lcl_uv,delthvu,ql_ad                     &
        , cin, CAPE,entrain_coef                                        &
        , qsat_lcl                                                      &
          )

ELSE

  ! Reset zh as done in conv_diag_comp for all points
  DO j=1,rows
    DO i=1,row_length
      zh(i,j) = z_half(i,j,2)
    END DO
  END DO
        
END IF        ! unstable points only

!-----------------------------------------------------------------------
! SCM output
!-----------------------------------------------------------------------
!----------------------------------------------------------------------------

9999   CONTINUE  ! Branch for error exit.

IF (lhook) CALL dr_hook('CONV_DIAG_0A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE conv_diag_0a

END MODULE conv_diag_0a_mod
