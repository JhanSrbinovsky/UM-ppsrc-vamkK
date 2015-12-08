! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine IDL_ideal_set_fields

      Subroutine IDL_ideal_set_fields_3a(                               &
     &                      model_domain, row_length, rows, n_rows      &
     &,                     model_levels, wet_model_levels              &
     &,                     off_x, off_y, halo_i, halo_j                &
     &,                     me, n_proc, at_extremity                    &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     theta_eh, rho_eh, exner_rho_levels_eh       &
     &,                     u, v, w, u_adv, v_adv, w_adv                &
     &,                     theta, q, qcl, qcf, rho                     &
     &,                     exner_rho_levels, exner_theta_levels        &
     &,                     p, p_theta_levels, p_star                   &
     &,                     L_include_halos)

! Purpose:
!      1. Copy large halo initial data from work arrays into data arrays
!      2. Set work arrays (qcl, qcf, w_adv) and w to zero
!      3. Initialise pressure arrays!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      USE earth_constants_mod, ONLY: g, earth_radius


      USE atmos_constants_mod, ONLY:                                    &
          r, p_zero, recip_kappa

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE


!  starting level for theta/q variables (= 0/1 for VATPOLES/ND) 
INTEGER, PARAMETER :: stlev = 0
!INTEGER, PARAMETER :: stlev = 1

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                             ! Local number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j                                                          &
                             ! Size of halo in j direction.
     &, off_x                                                           &
                             ! Size of small halo in i
     &, off_y                ! Size of small halo in j.

      Integer                                                           &
     &  me                                                              &
                   ! My processor number
     &, n_proc     ! Total number of processors

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_include_halos


      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)

      ! Work arrays with extended halos (_eh)
      ! Needed so that external halo values in LAMS can be set correctly
      ! for the lateral boundary arrays.
      REAL, Intent (InOut) ::                                           &
     &  theta_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           stlev:model_levels)                                          &
     &, rho_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &           model_levels)                                          &
     &, exner_rho_levels_eh(1-halo_i:row_length+halo_i,                 &
     &           1-halo_j:rows+halo_j, model_levels+1)

! Primary Arrays
      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &    model_levels)                                                 &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &    model_levels)                                                 &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &    0:model_levels)                                               &
     &, u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        model_levels)                                             &
     &, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,       &
     &        model_levels)                                             &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        0:model_levels)                                           &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &      model_levels)                                               &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &    model_levels+1)                                               &
     &, p_star(row_length, rows)                                        &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        stlev:model_levels)                                             &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &    stlev:wet_model_levels)                                             &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      stlev:wet_model_levels)                                           &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      stlev:wet_model_levels)                                           &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, model_levels+1)            &
     &, exner_theta_levels(1-off_x:row_length+off_x,1-off_y:rows+off_y, &
     &        model_levels)                                             &
     &, p_theta_levels(1-off_x:row_length+off_x,  1-off_y:rows+off_y,   &
     &        model_levels)

! local variables
      Integer                                                           &
     & i, j, k   ! loop counters

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('IDL_IDEAL_SET_FIELDS_3A',zhook_in,zhook_handle)

! First copy from big haloed work arrays into normal arrays
! u,v fields are in u_adv, v_adv
! DEPENDS ON: copy_field
      CALL COPY_FIELD(THETA_EH, THETA                                   &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels+1, model_levels+1, 1, model_levels+1 &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_p, .true., .false., .false.)
! DEPENDS ON: copy_field
      CALL COPY_FIELD(EXNER_RHO_LEVELS_EH, EXNER_RHO_LEVELS             &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels+1, model_levels+1, 1, model_levels+1 &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_p, .true., .false., .false.)
! DEPENDS ON: copy_field
      CALL COPY_FIELD(RHO_EH, RHO                                       &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels, model_levels, 1, model_levels       &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_p, .true., .false., .false.)
! DEPENDS ON: copy_field
      CALL COPY_FIELD(U_ADV, U                                          &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels, model_levels, 1, model_levels       &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_u, .true., .false., .true.)
! DEPENDS ON: copy_field
      CALL COPY_FIELD(V_ADV, V                                          &
     &,               row_length, row_length, n_rows, n_rows            &
     &,               model_levels, model_levels, 1, model_levels       &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_v, .true., .false., .true.)

! Set required arrays to zero, ( w_adv, qcl, qcf)
      w_adv(:,:,:) = 0.0
      qcl(:,:,:) = 0.0
      qcf(:,:,:) = 0.0

! Initialise remaining field w=0
! DEPENDS ON: copy_field
      CALL COPY_FIELD(W_ADV, W                                          &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels+1, model_levels+1, 1, model_levels+1 &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_p, .true., .false., .false.)
!
! Must now re-calculate the pressure-based variables, namely pressure
! on both rho and theta levels, exner on theta levels and p_star so that
! they are all consistent with the new LBC-updated values of exner on
! rho levels.
!
! DEPENDS ON: calc_exner_at_theta
      Call Calc_Exner_at_theta(                                         &
     &                      r_theta_levels, r_rho_levels,               &
     &                      Exner_rho_levels,                           &
     &                      row_length, rows, model_levels,             &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      Exner_theta_levels, L_include_halos)

! Calculate pressure from Exner at rho levels.
      Do k = 1, model_levels+1
        Do j = 1, rows
          Do i = 1, row_length
            p(i,j,k) = p_zero * exner_rho_levels(i,j,k)**recip_kappa
          End Do
        End Do
      End Do
! Calculate p_theta_levels from Exner at theta levels.
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            p_theta_levels(i,j,k) = p_zero *                            &
     &                      exner_theta_levels(i,j,k)**recip_kappa
          End Do
        End Do
      End Do
! Halos updated
! DEPENDS ON: swap_bounds
      call Swap_Bounds(p,                                               &
     &                   row_length, rows, model_levels+1,              &
     &                   off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
      call Swap_Bounds(p_theta_levels,                                  &
     &                   row_length, rows, model_levels,                &
     &                   off_x, off_y, fld_type_p, .false.)

! DEPENDS ON: calc_p_star
      Call Calc_P_star(                                                 &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   p, rho, row_length, rows, model_levels,        &
     &                   off_x, off_y, halo_i, halo_j,                  &
     &                   p_star)

      IF (lhook) CALL dr_hook('IDL_IDEAL_SET_FIELDS_3A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_ideal_set_fields_3a
