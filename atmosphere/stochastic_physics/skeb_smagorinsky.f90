! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE skeb_smagorinsky_mod

IMPLICIT NONE

CONTAINS

    SUBROUTINE skeb_smagorinsky(                                        &
                     rows, row_length, n_rows                           &
,                    model_levels                                       &
,                    u, v                                               &
,                    delta_lambda, delta_phi                            &
,                    disp                                               &
,                    timestep_number                                    &
                     )

! Calculates energy dissipation based on the Smagorinksy diffusion
! scheme, where the dissipation is proportional to D^3=SQRT(Ds+Dt)
! the shearing and tension stresses. This routine is derived from
! the subroutine "turb_smagorinky" used for subgrid turbulence.

USE dynamics_input_mod, ONLY: l_endgame

USE stochastic_physics_run_mod, ONLY: skeb2_toplev, skeb2_botlev

USE earth_constants_mod, ONLY: earth_radius

! Bounds of arrays
USE atm_fields_bounds_mod, ONLY:                                        &
    pdims, pdims_s, udims_s, vdims_s

! Model grid trigonometry
USE trignometric_mod, ONLY:                                             &
    cos_theta_latitude, cos_v_latitude

! Model level-height modules
USE level_heights_mod,     ONLY:                                        &
    r_rho_levels

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE PrintStatus_mod
USE UM_ParVars
IMPLICIT NONE
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


! Variables with Intent (In)

INTEGER, INTENT(IN) ::                                                  &
  row_length                                                            &
           ! number of points on a row.
, rows                                                                  &
           ! number of rows.
, n_rows                                                                &
           ! number of v rows.
, model_levels                                                          &
           ! number of model levels.
, timestep_number
           ! Forecast time-step

REAL, INTENT(IN) ::                                                     &
  u            (udims_s%i_start:udims_s%i_end,                          &
                udims_s%j_start:udims_s%j_end,                          &
                udims_s%k_start:udims_s%k_end)                          &
           ! main prog. u field
, v            (vdims_s%i_start:vdims_s%i_end,                          &
                vdims_s%j_start:vdims_s%j_end,                          &
                vdims_s%k_start:vdims_s%k_end)                          &
           ! main prog. v field
, delta_lambda                                                          &
, delta_phi

! Variables with Intent (Out)
 REAL, INTENT(OUT) ::                                                   &
  disp         (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)
           ! Numerical Dissipation (m**2.s**-3)

! Local parameters
INTEGER ::                                                              &
  i                                                                     &
           ! loop over row_length
, j                                                                     &
           ! loop over rows
, k                                                                     &
           ! loop over model levels
, ju                                                                    &
           ! special case in j-loop for P,U at pole
, jv                                                                    &
           ! special case in j-loop for V at pole
, im1                                                                   &
           ! i minus one
, jm1                                                                   &
           ! j minus one
, ip1                                                                   &
           ! i plus one
, jp1                                                                   &
           ! j plus one
, kp1

REAL ::                                                                 &
  delta_x, delta_y                                                      &
           ! grid-spacing in m
, smallp = 1.e-14                                                       &
           ! A small number within available precision
, kh = 0.28                                                             &
           ! Numerical dissipation tuning factor (used to get desired
           ! dissipation rate of 0.7 W/m^2
, rmlmax                                                                &
           ! horizontal gridlength
, shear        (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)                              &
           ! horizontal strain rate
, ssq          (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)                              &
           ! squared horizontal strain rate
, visc         (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)
           ! Viscosity=MAX(dx,dy)^2*Shear

REAL ::                                                                 &
  dx_rho       (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)                              &
           ! dx on rho points
, dx_rho_centre(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)                              &
           ! dx in centre of grid box on rho levels
, dy_rho       (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)                              &
           ! dy on rho points
, dy_rho_centre(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)                              &
           ! dy in centre of grid box on rho levels
, cx_rho       (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)                              &
           ! reciprocal of dx_rho
, cy_rho       (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)
           ! reciprocal of dy_rho

!  These values are calculated once at timestep=1 then saved
!  for further calls
!  NB****************************************************
!  NB** ASSUMING THAT RHO & THETA HEIGHTS ARE CONSTANT **
!  NB****************************************************
!  ------------------------------------------------------
!  :: Note about the grid points ::
!  -- U is on U-points and V on V-points (Arakawa C)
!  -- cf. pp 4.5-4.8 in http://www-nwp/~frax/ND3D_UM6p3_doc.html
!  -- Ouput DISP is on Pi-points
!  ------------------------------------------------------
REAL, ALLOCATABLE, SAVE ::                                              &
  cx_rho_centre(:,:,:), cy_rho_centre(:,:,:)                            &
           ! reciprocal of d[x|y]_rho_centre
, cx2_rho(:,:,:), cy2_rho(:,:,:)                                        &
           ! c[x|y]_rho^2
, latitude_mask(:,:)
           ! Theta-level thickness
           ! Mask-out high latitudes

REAL ::                                                                 &
  s11                                                                   &
           ! ii component of shear stress
, s22                                                                   &
           ! jj component of shear stress
, s12                                                                   &
           ! ij component of shear stress
, piover2                                                               &
           ! half of PI
, maxlat                                                                &
           ! Lat threshold for truncation of shear to zero at pole
, mylat
           ! Work value of Lat in Radians

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Initialise variables
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('SKEB_SMAGORINSKY',zhook_in,zhook_handle)

! Reciprocal of dx_rho_centre (with halo)
IF (.NOT.ALLOCATED(cx_rho_centre)) THEN
     ALLOCATE(cx_rho_centre(pdims_s%i_start:pdims_s%i_end,              &
                            pdims_s%j_start:pdims_s%j_end,              &
                            pdims_s%k_start:pdims_s%k_end))
END IF
! Reciprocal of dy_rho_centre (with halo)
IF (.NOT.ALLOCATED(cy_rho_centre)) THEN
     ALLOCATE(cy_rho_centre(pdims_s%i_start:pdims_s%i_end,              &
                            pdims_s%j_start:pdims_s%j_end,              &
                            pdims_s%k_start:pdims_s%k_end))
END IF
! Square of cx_rho (with halo)
IF (.NOT.ALLOCATED(cx2_rho)) THEN
     ALLOCATE(cx2_rho      (pdims_s%i_start:pdims_s%i_end,              &
                            pdims_s%j_start:pdims_s%j_end,              &
                            pdims_s%k_start:pdims_s%k_end))
END IF
! Square of cy_rho (with halo)
IF (.NOT.ALLOCATED(cy2_rho)) THEN
     ALLOCATE(cy2_rho      (pdims_s%i_start:pdims_s%i_end,              &
                            pdims_s%j_start:pdims_s%j_end,              &
                            pdims_s%k_start:pdims_s%k_end))
END IF

! Initialise dissipation array in parts not written later
DO k = 1, skeb2_botlev - 1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      disp(i,j,k) = 0.0
    END DO
  END DO
END DO
DO k = skeb2_toplev + 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      disp(i,j,k) = 0.0
    END DO
  END DO
END DO
delta_x = earth_radius * delta_lambda
delta_y = earth_radius * delta_phi

!----------------------------------------------------------------------
!   Hybrid-height model :: rho and theta-levels constant throughout run
!   Calculated variables saved for next call
!----------------------------------------------------------------------

IF  (timestep_number == 1) THEN

! Initialise parts of arrays outside SKEB levels
  DO k = pdims%k_start, skeb2_botlev - 1
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        cx_rho_centre(i,j,k) = 0.0
        cy_rho_centre(i,j,k) = 0.0
      END DO
    END DO
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        cy2_rho(i,j,k) = 0.0
      END DO
    END DO
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        cx2_rho(i,j,k) = 0.0
      END DO
    END DO
  END DO
  DO k = skeb2_toplev + 1, pdims%k_end
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        cx_rho_centre(i,j,k) = 0.0
        cy_rho_centre(i,j,k) = 0.0
      END DO
    END DO
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        cy2_rho(i,j,k) = 0.0
      END DO
    END DO
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        cx2_rho(i,j,k) = 0.0
      END DO
    END DO
  END DO

! Mask for latitudes > MAXLAT deg
! Decreases field with a COSINE-shaped drop off across a half-wave
!  cycle from 1 at maxlat to 0 at pole
  piover2=ACOS(0.)
  maxlat=1.46   ! ~84 deg = Arbitrary choice based on energy plots
  IF (.NOT.ALLOCATED(latitude_mask)) THEN
       ALLOCATE(latitude_mask(pdims%i_start:pdims%i_end,                &
                            pdims%j_start:pdims%j_end))
! Default global value of one
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      latitude_mask(i,j) = 1.
    END DO
  END DO
! Apply ramp near poles
  DO j = pdims%j_start, pdims%j_end
    mylat = (piover2 - MAX(maxlat,                                      &
          ACOS(cos_theta_latitude(pdims%i_start,j))))/(piover2 - maxlat)
    mylat = 2. * piover2 * mylat     ! convert span over 0 <-> 180 deg
    DO i = pdims%i_start, pdims%i_end
      latitude_mask(i,j) = 0.5 - 0.5 * COS(mylat)
    END DO
  END DO
! Set value to zero at pole (ND grid) or last theta-row (EG grid)
! to avoid spurious turbulence values close to pole
  IF (at_extremity(psouth)) THEN
    DO i = pdims%i_start, pdims%i_end
      latitude_mask(i, pdims%j_start) = 0.
    END DO
  END IF
  IF (at_extremity(pnorth)) THEN
    DO i = pdims%i_start, pdims%i_end
      latitude_mask(i, pdims%j_end) = 0.
    END DO
  END IF
END IF

!----------------------------------------------------------------
! Calculate grid spacings
!----------------------------------------------------------------
! Horizontal grid spacings used in calculation of rate of strain
! term delta_lambda and delta_phi are in radians
!

IF (.NOT. l_endgame) THEN
DO k = skeb2_botlev, skeb2_toplev + 1

  DO j = pdims%j_start, pdims%j_end
    ju = j
    jv = j
    jp1 = j + 1
!   For Pole: use cos_theta_latitude value at adjacent row
    IF (at_extremity(psouth)) ju = MAX(ju, pdims%j_start + 1)
    IF (at_extremity(pnorth)) ju = MIN(ju, pdims%j_end - 1)

    DO i = pdims%i_start, pdims%i_end
      ip1 = i + 1

      dx_rho(i,j,k) =                                                   &
               r_rho_levels(i,j,k)                                      &
             * delta_lambda * cos_theta_latitude(i,ju)

      dy_rho(i,j,k) =                                                   &
               r_rho_levels(i,j,k) * delta_phi

      dy_rho_centre(i,j,k) = 0.25*(                                     &
         r_rho_levels(i,jp1,k) + r_rho_levels(ip1,jp1,k)                &
        +r_rho_levels(i,j,k) + r_rho_levels(ip1,j,k) )                  &
        * delta_phi

      cx_rho(i,j,k) = 1./dx_rho(i,j,k)
      cy_rho(i,j,k) =  1./dy_rho(i,j,k)
      cy_rho_centre(i,j,k) = 1./dy_rho_centre(i,j,k)
      cx2_rho(i,j,k) = cx_rho(i,j,k)*cx_rho(i,j,k)
      cy2_rho(i,j,k) = cy_rho(i,j,k)*cy_rho(i,j,k)

      dx_rho_centre(i,j,k) = 0.25*(                                     &
         r_rho_levels(i,jp1,k) + r_rho_levels(ip1,jp1,k)                &
        +r_rho_levels(i,j,k) + r_rho_levels(ip1,j,k) )                  &
                * delta_lambda * cos_v_latitude(i,jv)
      cx_rho_centre(i,j,k) =  1./dx_rho_centre(i,j,k)

    END DO
  END DO

  IF (at_extremity(pnorth)) THEN
    DO i = pdims%i_start, pdims%i_end
      cx_rho_centre(i,pdims%j_end,k) =                                  &
               1./dx_rho_centre(i,pdims%j_end-1,k)
    END DO
  END IF
END DO
ELSE
  DO k = skeb2_botlev, skeb2_toplev + 1

    DO j = pdims%j_start, pdims%j_end
      ju = j
      jv = j
      jp1 = j + 1
  !   Trap in case cos_theta_latitude=0, use value at adjacent row
      IF (at_extremity(psouth).AND.cos_theta_latitude(1,ju).LT.1.0e-3)  &
          ju = j + 1
      IF (at_extremity(pnorth).AND.cos_theta_latitude(1,ju).LT.1.0e-3)  &
          ju = j - 1
  !   Trap in case cos_v_latitude~0, use value at adjacent row
      IF (at_extremity(psouth)) jv = MAX(jv, pdims%j_start + 1)
      IF (at_extremity(pnorth)) jv = MIN(jv, pdims%j_end - 1)

      DO i = pdims%i_start, pdims%i_end
        ip1 = i + 1

        dx_rho(i,j,k) =                                                 &
                 r_rho_levels(i,j,k)                                    &
               * delta_lambda * cos_theta_latitude(i,ju)

        dy_rho(i,j,k) =                                                 &
                 r_rho_levels(i,j,k) * delta_phi

        dy_rho_centre(i,j,k) = 0.25*(                                   &
           r_rho_levels(i,jp1,k) + r_rho_levels(ip1,jp1,k)              &
          +r_rho_levels(i,j,k) + r_rho_levels(ip1,j,k) )                &
          * delta_phi

        cx_rho(i,j,k) = 1./dx_rho(i,j,k)
        cy_rho(i,j,k) =  1./dy_rho(i,j,k)
        cy_rho_centre(i,j,k) = 1./dy_rho_centre(i,j,k)
        cx2_rho(i,j,k) = cx_rho(i,j,k)*cx_rho(i,j,k)
        cy2_rho(i,j,k) = cy_rho(i,j,k)*cy_rho(i,j,k)

        dx_rho_centre(i,j,k) = 0.25*(                                   &
           r_rho_levels(i,jp1,k) + r_rho_levels(ip1,jp1,k)              &
          +r_rho_levels(i,j,k) + r_rho_levels(ip1,j,k) )                &
                  * delta_lambda * cos_v_latitude(i,jv)
        cx_rho_centre(i,j,k) =  1./dx_rho_centre(i,j,k)

      END DO
    END DO

  END DO
END IF

! DEPENDS ON: swap_bounds
CALL swap_bounds(cx2_rho, row_length, rows, model_levels,               &
                 offx, offy, fld_type_p, .FALSE.)

! DEPENDS ON: swap_bounds
CALL swap_bounds(cx_rho_centre, row_length, rows, model_levels,         &
                 offx, offy, fld_type_v, .FALSE.)

! DEPENDS ON: swap_bounds
CALL swap_bounds(cy2_rho, row_length, rows, model_levels,               &
                 offx, offy, fld_type_p, .FALSE.)

! DEPENDS ON: swap_bounds
CALL swap_bounds(cy_rho_centre, row_length, rows, model_levels,         &
                 offx, offy, fld_type_v, .FALSE.)

IF  (printstatus  >  prstatus_normal) THEN
  WRITE(6,*) '*********************************************'
  WRITE(6,*) 'Smagorinsky TimeStep=1'
  WRITE(6,*) 'max val cx2_rho', MAXVAL(cx2_rho)
  WRITE(6,*) 'max val cy2_rho', MAXVAL(cy2_rho)
  WRITE(6,*) 'max val cx_rho_centre', MAXVAL(cx_rho_centre)
  WRITE(6,*) 'max val cy_rho_centre', MAXVAL(cy_rho_centre)
  WRITE(6,*)
END IF

END IF    ! timestep_number == 1
!----------------------------------------------------------------------


IF  (printstatus  >  prstatus_oper) THEN
  WRITE(6,*) '*********************************************'
  WRITE(6,*) 'STPH: skeb_smagorinsky TimeStep=',timestep_number
  WRITE(6,*) 'max val cx2_rho', MAXVAL(cx2_rho)
  WRITE(6,*) 'max val cy2_rho', MAXVAL(cy2_rho)
  WRITE(6,*) 'max val cx_rho_centre', MAXVAL(cx_rho_centre)
  WRITE(6,*) 'max val cy_rho_centre', MAXVAL(cy_rho_centre)
  WRITE(6,*) 'max loc cx2_rho', MAXLOC(cx2_rho)
  WRITE(6,*) 'max loc cy2_rho', MAXLOC(cy2_rho)
  WRITE(6,*) 'max loc cx_rho_centre', MAXLOC(cx_rho_centre)
  WRITE(6,*) 'max loc cy_rho_centre', MAXLOC(cy_rho_centre)
  WRITE(6,*)
END IF


! Now calculate 3-components of the deviatoric stress
DO k = skeb2_botlev, skeb2_toplev
  kp1 = k + 1

  DO j = pdims%j_start, pdims%j_end
    jm1 = j - 1
    jp1 = j + 1
    ! Prevent array access at pole 
    IF (at_extremity(psouth)) jm1 = MAX(jm1, pdims%j_start + 1) 
    IF (at_extremity(pnorth)) jp1 = MIN(jp1, pdims%j_end - 1)

    DO i = pdims%i_start, pdims%i_end
      im1 = i - 1
      ip1 = i + 1

      s11 =                                                             &
       cx2_rho(i,j,kp1)*(u(i,j,kp1)-u(im1,j,kp1))**2 +                  &
       cx2_rho(i,j,k)*(u(i,j,k)-u(im1,j,k))**2

      s22 =                                                             &
       cy2_rho(i,j,kp1)*(v(i,j,kp1)-v(i,jm1,kp1))**2 +                  &
       cy2_rho(i,j,k)*(v(i,j,k)-v(i,jm1,k))**2

      s12 = (  (  (                                                     &
        ((u(im1,j,k)-u(im1,jm1,k))*cy_rho_centre(im1,jm1,k)             &
       +(v(i,jm1,k)-v(im1,jm1,k))*cx_rho_centre(im1,jm1,k))**2 +        &
        ((u(im1,jp1,k)-u(im1,j,k))*cy_rho_centre(im1,j,k)               &
       +(v(i,j,k)-v(im1,j,k))*cx_rho_centre(im1,j,k))**2                &
              )  +  (                                                   &
        ((u(i,j,k)-u(i,jm1,k))*cy_rho_centre(i,jm1,k)                   &
       +(v(ip1,jm1,k)-v(i,jm1,k))*cx_rho_centre(i,jm1,k))**2 +          &
        ((u(i,jp1,k)-u(i,j,k))*cy_rho_centre(i,j,k)                     &
       +(v(ip1,j,k)-v(i,j,k))*cx_rho_centre(i,j,k))**2                  &
         ) )  +  ( (                                                    &
        ((u(im1,j,kp1)-u(im1,jm1,kp1))*cy_rho_centre(im1,jm1,kp1)       &
       +(v(i,jm1,kp1)-v(im1,jm1,kp1))                                   &
         *cx_rho_centre(im1,jm1,kp1))**2 +                              &
        ((u(im1,jp1,kp1)-u(im1,j,kp1))*cy_rho_centre(im1,j,kp1)         &
       +(v(i,j,kp1)-v(im1,j,kp1))*cx_rho_centre(im1,j,kp1))**2          &
              )  +  (                                                   &
        ((u(i,j,kp1)-u(i,jm1,kp1))*cy_rho_centre(i,jm1,kp1)             &
       +(v(ip1,jm1,kp1)-v(i,jm1,kp1))                                   &
        *cx_rho_centre(i,jm1,kp1))**2 +                                 &
        ((u(i,jp1,kp1)-u(i,j,kp1))*cy_rho_centre(i,j,kp1)               &
       +(v(ip1,j,kp1)-v(i,j,kp1))*cx_rho_centre(i,j,kp1))**2            &
          ) ) )*0.125      ! _averaging s12 over 8 points

      ssq(i,j,k) = s11 + s22 + s12 + smallp
      shear(i,j,k) = SQRT(ssq(i,j,k))

    END DO   ! on i
  END DO   ! on j
END DO   ! on k

! Numerical tuning factor * horizontal gridlength
rmlmax = kh * MAX(delta_x, delta_y)

DO k = skeb2_botlev, skeb2_toplev
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      visc(i,j,k) = (rmlmax**2)*shear(i,j,k)
!     Reduce values to zero from latitudes from maxlat up to pole
!     to avoid spurious values as dx -> 0
      disp(i,j,k) = ssq(i,j,k) * visc(i,j,k) * latitude_mask(i,j)

    END DO
  END DO
END DO

IF (lhook) CALL dr_hook('SKEB_SMAGORINSKY',zhook_out,zhook_handle)
RETURN
END SUBROUTINE skeb_smagorinsky
END MODULE skeb_smagorinsky_mod
