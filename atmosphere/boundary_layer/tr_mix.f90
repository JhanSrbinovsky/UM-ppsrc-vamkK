! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE TR_MIX ------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Boundary Layer
!
!    Purpose: Calculate tracer flux and pass through to IMP_MIX to solve
!
!    Programming standard: UM Documentation Paper No 3
!
!    Documentation: UM Documentation Paper No 24.
!
!*----------------------------------------------------------------------
MODULE tr_mix_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE tr_mix (                                                     &
! IN fields
 bl_levels,GAMMA,rhokh_rdz,rhokh_1,dtrdz,surf_em,res_factor,            &
 kent, we_lim, t_frac, zrzi,                                            &
 kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,zh ,zhsc, z_uv,             &
! INOUT / OUT fields
 field,f_field,surf_dep_flux                                            &
)

  USE atm_fields_bounds_mod, ONLY: pdims, tdims
  USE imp_mix_mod, ONLY: imp_mix
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

! arguments passed IN
  INTEGER, INTENT(IN) ::                                                &
    bl_levels   ! IN No. of atmospheric levels for boundary layer

  REAL, INTENT(IN) ::                                                   &
   GAMMA(bl_levels),                                                    &
   rhokh_rdz(pdims%i_start:,pdims%j_start:, 2:),                        &
                                   ! IN Mixing coeff. above surface
                                   !    = RHOKH(,K)*RDZ(K)
                                   !    for K>=2 (from IMP_SOLVER).
   rhokh_1(pdims%i_start:,pdims%j_start:),                              &
                                   ! IN  Surface exchange coeff.
                                   !     from P243 (SF_EXCH)
   dtrdz(pdims%i_start:,pdims%j_start:, :),                             &
                                   ! IN  dt/(rho*r*r*dz) for scalar
                                   !     flux divergence
   surf_em(pdims%i_start:,pdims%j_start:),                              &
                                   ! IN, Surface emissions in kg/m2/s
   res_factor(pdims%i_start:,pdims%j_start:),                           &
                                   ! IN, dry dep coeff=Ra/(Ra+Rb+Rc)
    we_lim(pdims%i_start:,pdims%j_start:,:),                            &
                                   ! IN rho*entrainment rate implied by
                                   !     placing of subsidence
    zrzi(pdims%i_start:,pdims%j_start:,:),                              &
                                   ! IN (z-z_base)/(z_i-z_base)
    t_frac(pdims%i_start:,pdims%j_start:,:),                            &
                                   ! IN a fraction of the timestep
    we_lim_dsc(pdims%i_start:,pdims%j_start:,:),                        &
                                   ! IN rho*entrainment rate implied by
                                   !     placing of subsidence
    zrzi_dsc(pdims%i_start:,pdims%j_start:,:),                          &
                                   ! IN (z-z_base)/(z_i-z_base)
    t_frac_dsc(pdims%i_start:,pdims%j_start:,:),                        &
                                   ! IN a fraction of the timestep
    z_uv(pdims%i_start:,pdims%j_start:,:),                              &
                                   ! IN Z_uv(*,K) is height of half
                                   !    level k-1/2.
    zhsc(pdims%i_start:,pdims%j_start:),                                &
                                   ! IN Top of decoupled layer
    zh(pdims%i_start:,pdims%j_start:)
                                   ! IN Top of surface mixed layer

  INTEGER, INTENT(IN) ::                                                &
    kent(pdims%i_start:,pdims%j_start:),                                &
                                   ! IN grid-level of SML inversion
    kent_dsc(pdims%i_start:,pdims%j_start:)
                                   ! IN grid-level of DSC inversion

! INOUT arguments
  REAL, INTENT(INOUT) ::                                                &
   field(pdims%i_start:,pdims%j_start:,:)
                                   ! INOUT Tracer amount in kg/kg.

! OUT arguments
  REAL, INTENT(OUT) ::                                                  &
   f_field(pdims%i_start:,pdims%j_start:,:),                            &
                                   ! OUT Flux of tracer in kg/m2/s.
   surf_dep_flux(pdims%i_start:,pdims%j_start:)
                                   ! OUT, surface deposn flux (kg/m2/s)

!    Local and other symbolic constants :-

  REAL,PARAMETER:: one=1.0
  REAL,PARAMETER:: smallp=TINY(one)

!   Workspace :-

  REAL ::                                                               &
   rhok_dep(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                    ! Surface deposition coeficient
   gamma_rhokh_rdz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, &
                   2:bl_levels),                                        &
                             ! gamma*RHOKH_RDZ
   dz_disc,                                                             &
                             ! Temporary in subgrid zi calculation
   dzlkp1,                                                              &
                             ! Thickness of theta-level K+1
   f_field_ent,                                                         &
                             ! Time-level n entrainment flux
   dfield_inv            ! inversion jump

! Arrays for vectorisation
  REAL, DIMENSION(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)::&
   dfield_sml,                                                          &
                             ! Jump in field across SML inversion
   dfield_dsc            ! Jump in field across DSC inversion
!*
!  Local scalars :-
  INTEGER ::                                                            &
   i,j,                                                                 &
                ! Loop counter (horizontal field index).
   k,                                                                   &
                ! Loop counter (vertical index).
   km1,                                                                 &
                ! Max(K-1,2)
   ient     ! Loop counter for entrainment

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
!   0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------
  IF (lhook) CALL dr_hook('TR_MIX',zhook_in,zhook_handle)


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,j,i,dzlkp1,dz_disc)

!-----------------------------------------------------------------------
!  1.  Calculate flux of tracer:
!-----------------------------------------------------------------------
!  1.1 Above the surface
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        gamma_rhokh_rdz(i,j,k) = GAMMA(k) * rhokh_rdz(i,j,k)
        f_field(i,j,k) = - rhokh_rdz(i,j,k) *                           &
                           (field(i,j,k) - field(i,j,k-1))
      END DO
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!  1.2 At the surface: (i) set surface flux equal to input emissions
!                    (should be passed in as ancillary file, else ZERO)
!                      (ii) Use input resistance factors to calculate
!                    surface deposition (if ZERO then no dry deposition)
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      f_field(i,j,1) = surf_em(i,j)  ! Inject surface emissions
      rhok_dep(i,j) = res_factor(i,j) * rhokh_1(i,j)
      surf_dep_flux(i,j) = -rhok_dep(i,j) * field(i,j,1)
      rhok_dep(i,j) = GAMMA(1) * rhok_dep(i,j)
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Add on explicit entrainment fluxes
! These are calculated from the time-level n profile and parametrized
! entrainment rate to give a time-level n flux (F_FIELD_ENT).  This is
! then implemented using an equivalent RHOKH as the diagnosed inversion
! jump can change significantly across the timestep and so, therefore,
! should the parametrized flux.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k = kent(i,j)-1     ! equal to originally diagnosed NTML
        !----------------------------
        ! diagnose SML inversion jump
        !----------------------------
      IF (k  ==  bl_levels-1) THEN
        dfield_sml(i,j) = field(i,j,k+1) - field(i,j,k)
      ELSE
        dzlkp1  = z_uv(i,j,k+2) - z_uv(i,j,k+1)
        dz_disc = z_uv(i,j,k+2) - zh(i,j)
        IF (dz_disc/dzlkp1  >   0.1) THEN
          dfield_sml(i,j) = (field(i,j,k+1)-field(i,j,k))               &
                              * dzlkp1 /dz_disc

          IF ( field(i,j,k+2)  >   field(i,j,k+1) .AND.                 &
                     field(i,j,k+1)  >   field(i,j,k) ) THEN
            dfield_sml(i,j) = MIN( field(i,j,k+2)-field(i,j,k),         &
                                   dfield_sml(i,j) )
          ELSE IF ( field(i,j,k+2)  <   field(i,j,k+1) .AND.            &
                    field(i,j,k+1)  <   field(i,j,k) ) THEN
            dfield_sml(i,j) = MAX( field(i,j,k+2)-field(i,j,k),         &
                                   dfield_sml(i,j) )
          ELSE  ! FIELD non-monotonic
            dfield_sml(i,j) = field(i,j,k+1)-field(i,j,k)
          END IF
        ELSE
          dfield_sml(i,j) = field(i,j,k+2) - field(i,j,k)
        END IF
      END IF
    END DO
  END DO
!OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      k = kent_dsc(i,j)-1     ! equal to originally diagnosed NTDSC

        !----------------------------
        ! diagnose DSC inversion jump
        !----------------------------
      IF (k  ==  bl_levels-1) THEN
        dfield_dsc(i,j) = field(i,j,k+1) - field(i,j,k)
      ELSE IF ( kent_dsc(i,j)  >   2 ) THEN

        dzlkp1  = z_uv(i,j,k+2) - z_uv(i,j,k+1)
        dz_disc = z_uv(i,j,k+2) - zhsc(i,j)
        IF (dz_disc/dzlkp1  >   0.1) THEN
          dfield_dsc(i,j) = (field(i,j,k+1)-field(i,j,k))               &
                              * dzlkp1 /dz_disc

          IF ( field(i,j,k+2)  >   field(i,j,k+1) .AND.                 &
               field(i,j,k+1)  >   field(i,j,k) ) THEN
            dfield_dsc(i,j) = MIN( field(i,j,k+2)-field(i,j,k),         &
                                   dfield_dsc(i,j) )
          ELSE IF ( field(i,j,k+2)  <   field(i,j,k+1) .AND.            &
                    field(i,j,k+1)  <   field(i,j,k) ) THEN
            dfield_dsc(i,j) = MAX( field(i,j,k+2)-field(i,j,k),         &
                                   dfield_dsc(i,j) )
          ELSE  ! FIELD non-monotonic
            dfield_dsc(i,j) = field(i,j,k+1)-field(i,j,k)
          END IF
        ELSE
          dfield_dsc(i,j) = field(i,j,k+2) - field(i,j,k)
        END IF
      END IF
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

      !------------------------------------
      ! calculate entrainment fluxes and KH
      !------------------------------------
  DO ient = 1, 3

!OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k,km1,    &
!OMP& f_field_ent, dfield_inv)
    DO j = pdims%j_start, pdims%j_end
!CDIR nodep
      DO i = pdims%i_start, pdims%i_end

        k = kent(i,j)-2+ient
        IF ( k > 1 .AND. k <= bl_levels .AND.                           &
                         t_frac(i,j,ient) > 0.0) THEN
          IF ( ABS(field(i,j,k)-field(i,j,k-1)) >= smallp ) THEN
            dfield_inv = dfield_sml(i,j)
            ! DFIELD_INV must have same sign as
            ! local gradient to get right sign of flux
            IF ( dfield_sml(i,j) / (field(i,j,k)-field(i,j,k-1))        &
                  <   0.0 ) dfield_inv = field(i,j,k)-field(i,j,k-1)
            f_field_ent = t_frac(i,j,ient) * ( f_field(i,j,1)           &
              - ( we_lim(i,j,ient) * dfield_inv + f_field(i,j,1) )      &
                * zrzi(i,j,ient) )
            ! interpolation to surface flux must not change the sign
            ! of the entrainment flux otherwise KH will be <0!
            IF ( f_field_ent / (field(i,j,k)-field(i,j,k-1))  >   0.0 ) &
              f_field_ent = - t_frac(i,j,ient) *                        &
                        we_lim(i,j,ient) * dfield_inv * zrzi(i,j,ient)
            ! Restrict size of RHOKH for numerical safety
            km1 = MAX( k-1, 2 )
            gamma_rhokh_rdz(i,j,k) = gamma_rhokh_rdz(i,j,k) +           &
                    MIN( - GAMMA(k)*f_field_ent/                        &
                           ( field(i,j,k)-field(i,j,k-1) ),             &
                         10.0*gamma_rhokh_rdz(i,j,km1) )
            ! Recalculate explicit flux using entrainment KH
            f_field(i,j,k) = - (gamma_rhokh_rdz(i,j,k) / GAMMA(k)) *    &
                               (field(i,j,k) - field(i,j,k-1))
          END IF
        END IF
      END DO
    END DO
!OMP END PARALLEL DO
  END DO   ! IENT

  DO ient = 1, 3

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k,km1,    &
!$OMP& f_field_ent, dfield_inv)
    DO j = pdims%j_start, pdims%j_end
!CDIR nodep
      DO i = pdims%i_start, pdims%i_end

        k = kent_dsc(i,j)-2+ient
        IF ( kent_dsc(i,j) >= 3 .AND. k <= bl_levels ) THEN
          IF ( t_frac_dsc(i,j,ient) > 0.0                               &
               .AND. ABS(field(i,j,k)-field(i,j,k-1)) >= smallp ) THEN
            dfield_inv = dfield_dsc(i,j)
              ! DFIELD_INV must have same sign as
              ! local gradient to get right sign of flux
            IF ( dfield_dsc(i,j) / (field(i,j,k)-field(i,j,k-1))        &
                  <   0.0 ) dfield_inv = field(i,j,k)-field(i,j,k-1)
            f_field_ent = - t_frac_dsc(i,j,ient) *                      &
               we_lim_dsc(i,j,ient) * dfield_inv * zrzi_dsc(i,j,ient)
              ! Restrict size of RHOKH for numerical safety
            km1 = MAX( k-1, 2 )
            gamma_rhokh_rdz(i,j,k) = gamma_rhokh_rdz(i,j,k) +           &
                     MIN( - GAMMA(k)*f_field_ent/                       &
                                      ( field(i,j,k)-field(i,j,k-1) ),  &
                          10.0*gamma_rhokh_rdz(i,j,km1) )
              ! Recalculate explicit flux using entrainment KH
            f_field(i,j,k) = - (gamma_rhokh_rdz(i,j,k) / GAMMA(k)) *    &
                             (field(i,j,k) - field(i,j,k-1))
          END IF
        END IF

      END DO
    END DO
!$OMP END PARALLEL DO
  END DO ! IENT

!-----------------------------------------------------------------------
!  2.  Call routine IMPL_CAL to calculate incrememnts to tracer field
!      and suface deposition flux for output
!-----------------------------------------------------------------------

  CALL imp_mix (                                                        &
   bl_levels,dtrdz,                                                     &
   gamma_rhokh_rdz,rhok_dep,f_field,surf_dep_flux,field                 &
   )

  IF (lhook) CALL dr_hook('TR_MIX',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE tr_mix
END MODULE tr_mix_mod
