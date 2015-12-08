! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

   SUBROUTINE init_lbc_dynamics(u,v,w, exner, rho, theta, etadot,              &
                                m_v, m_cl, m_cf, m_r, m_cf2, m_gr,             &
                                u_lbc, v_lbc, w_lbc, rho_lbc, exner_lbc,       &
                                theta_lbc, mv_lbc, mcl_lbc, mcf_lbc,           &
                                mr_lbc, mcf2_lbc, mgr_lbc,                     &
                                u_lbc_tnd, v_lbc_tnd, w_lbc_tnd,               &
                                rho_lbc_tnd, exner_lbc_tnd,                    &
                                theta_lbc_tnd, mv_lbc_tnd, mcl_lbc_tnd,        &
                                mcf_lbc_tnd, mr_lbc_tnd,                       &
                                mcf2_lbc_tnd, mgr_lbc_tnd,                     &
                                inc_fact, rim_stepsa, timestep_number,         &
                                row_length, rows, n_rows, model_levels,        &
                                offx, offy, halo_i, halo_j,                    &
                                l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,         &
                                lenrim, rimwidth, rimweights,                  &
                                lbc_size, lbc_start,                           &
                                at_extremity,                                  &
                                L_do_boundaries, L_do_halos                    &
                               )

   USE um_parparams
   USE rimtypes
   USE field_types
   USE atm_fields_bounds_mod
   USE metric_terms_mod
   USE horiz_grid_mod, ONLY : intw_rho2w, intw_u2p, intw_v2p

   USE yomhook, ONLY: lhook, dr_hook
   USE parkind1, ONLY: jprb, jpim

   IMPLICIT NONE

!
! Description: Update tendancies for LBC's and copy into field
!              (used by ENDGame to simplify atm_step).
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

   INTEGER, INTENT(IN)    :: row_length, rows, n_rows, model_levels
   INTEGER, INTENT(IN)    :: offx, offy, halo_i, halo_j
   INTEGER, INTENT(IN)    :: rim_stepsa, timestep_number
   INTEGER, INTENT(IN)    :: rimwidth, lenrim(Nfld_max,NHalo_max),             &
                                   lbc_size(4,Nfld_max,NHalo_max),             &
                                  lbc_start(4,Nfld_max,NHalo_max)

   REAL,    INTENT(IN)    :: rimweights(rimwidth)

   LOGICAL, INTENT(IN)    :: l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup
   REAL,    INTENT(INOUT) :: inc_fact
   LOGICAL, INTENT(INOUT) :: at_extremity(4)
   LOGICAL, INTENT(INOUT) :: L_do_boundaries, L_do_halos

   REAL,    INTENT(INOUT) :: u(udims_s%i_start:udims_s%i_end,                  &
                               udims_s%j_start:udims_s%j_end,                  &
                               udims_s%k_start:udims_s%k_end)
   REAL,    INTENT(INOUT) :: v(vdims_s%i_start:vdims_s%i_end,                  &
                               vdims_s%j_start:vdims_s%j_end,                  &
                               vdims_s%k_start:vdims_s%k_end)
   REAL,    INTENT(INOUT) :: w(wdims_s%i_start:wdims_s%i_end,                  &
                               wdims_s%j_start:wdims_s%j_end,                  &
                               wdims_s%k_start:wdims_s%k_end),                 &
                        etadot(wdims_s%i_start:wdims_s%i_end,                  &
                               wdims_s%j_start:wdims_s%j_end,                  &
                               wdims_s%k_start:wdims_s%k_end)
   REAL,    INTENT(INOUT) :: rho(pdims_s%i_start:pdims_s%i_end,                &
                                 pdims_s%j_start:pdims_s%j_end,                &
                                 pdims_s%k_start:pdims_s%k_end),               &
                             exner(pdims_s%i_start:pdims_s%i_end,              &
                                   pdims_s%j_start:pdims_s%j_end,              &
                                   pdims_s%k_start:pdims_s%k_end+1)
   REAL,    INTENT(INOUT) :: theta(tdims_s%i_start:tdims_s%i_end,              &
                                   tdims_s%j_start:tdims_s%j_end,              &
                                   tdims_s%k_start:tdims_s%k_end)
   REAL,    INTENT(INOUT) :: m_v(qdims_s%i_start:qdims_s%i_end,                &
                                 qdims_s%j_start:qdims_s%j_end,                &
                                 qdims_s%k_start:qdims_s%k_end),               &
                             m_r(qdims_s%i_start:qdims_s%i_end,                &
                                 qdims_s%j_start:qdims_s%j_end,                &
                                 qdims_s%k_start:qdims_s%k_end),               &
                            m_cl(qdims_s%i_start:qdims_s%i_end,                &
                                 qdims_s%j_start:qdims_s%j_end,                &
                                 qdims_s%k_start:qdims_s%k_end),               &
                            m_cf(qdims_s%i_start:qdims_s%i_end,                &
                                 qdims_s%j_start:qdims_s%j_end,                &
                                 qdims_s%k_start:qdims_s%k_end),               &
                            m_gr(qdims_s%i_start:qdims_s%i_end,                &
                                 qdims_s%j_start:qdims_s%j_end,                &
                                 qdims_s%k_start:qdims_s%k_end),               &
                           m_cf2(qdims_s%i_start:qdims_s%i_end,                &
                                 qdims_s%j_start:qdims_s%j_end,                &
                                 qdims_s%k_start:qdims_s%k_end)

   REAL,    INTENT (IN)   :: u_lbc(lenrim(fld_type_u,halo_type_extended),      &
                                   udims_s%k_start:udims_s%k_end),             &
                         u_lbc_tnd(lenrim(fld_type_u,halo_type_extended),      &
                                   udims_s%k_start:udims_s%k_end)
   REAL,    INTENT (IN)   :: v_lbc(lenrim(fld_type_v,halo_type_extended),      &
                                   vdims_s%k_start:vdims_s%k_end),             &
                         v_lbc_tnd(lenrim(fld_type_v,halo_type_extended),      &
                                   vdims_s%k_start:vdims_s%k_end)
   REAL,    INTENT (IN)   :: w_lbc(lenrim(fld_type_p,halo_type_extended),      &
                                   wdims_s%k_start:wdims_s%k_end),             &
                         w_lbc_tnd(lenrim(fld_type_p,halo_type_extended),      &
                                   wdims_s%k_start:wdims_s%k_end)

   REAL,    INTENT (IN)   :: rho_lbc(lenrim(fld_type_p,halo_type_extended),    &
                                     pdims_s%k_start:pdims_s%k_end),           &
                         rho_lbc_tnd(lenrim(fld_type_p,halo_type_extended),    &
                                     pdims_s%k_start:pdims_s%k_end),           &
                           exner_lbc(lenrim(fld_type_p,halo_type_extended),    &
                                     pdims_s%k_start:pdims_s%k_end+1),         &
                       exner_lbc_tnd(lenrim(fld_type_p,halo_type_extended),    &
                                     pdims_s%k_start:pdims_s%k_end+1)

   REAL,    INTENT (IN)   :: theta_lbc(lenrim(fld_type_p,halo_type_extended),  &
                                       tdims_s%k_start:tdims_s%k_end),         &
                         theta_lbc_tnd(lenrim(fld_type_p,halo_type_extended),  &
                                       tdims_s%k_start:tdims_s%k_end)

   REAL,    INTENT (IN)   :: mv_lbc(lenrim(fld_type_p,halo_type_extended),     &
                                    qdims_s%k_start:qdims_s%k_end),            &
                             mr_lbc(lenrim(fld_type_p,halo_type_extended),     &
                                    qdims_s%k_start:qdims_s%k_end),            &
                            mcf_lbc(lenrim(fld_type_p,halo_type_extended),     &
                                    qdims_s%k_start:qdims_s%k_end),            &
                            mcl_lbc(lenrim(fld_type_p,halo_type_extended),     &
                                    qdims_s%k_start:qdims_s%k_end),            &
                           mcf2_lbc(lenrim(fld_type_p,halo_type_extended),     &
                                    qdims_s%k_start:qdims_s%k_end),            &
                            mgr_lbc(lenrim(fld_type_p,halo_type_extended),     &
                                    qdims_s%k_start:qdims_s%k_end)

   REAL,    INTENT (IN)   :: mv_lbc_tnd(lenrim(fld_type_p,halo_type_extended), &
                                        qdims_s%k_start:qdims_s%k_end),        &
                             mr_lbc_tnd(lenrim(fld_type_p,halo_type_extended), &
                                        qdims_s%k_start:qdims_s%k_end),        &
                            mcf_lbc_tnd(lenrim(fld_type_p,halo_type_extended), &
                                        qdims_s%k_start:qdims_s%k_end),        &
                            mcl_lbc_tnd(lenrim(fld_type_p,halo_type_extended), &
                                        qdims_s%k_start:qdims_s%k_end),        &
                           mcf2_lbc_tnd(lenrim(fld_type_p,halo_type_extended), &
                                        qdims_s%k_start:qdims_s%k_end),        &
                            mgr_lbc_tnd(lenrim(fld_type_p,halo_type_extended), &
                                        qdims_s%k_start:qdims_s%k_end)
!
   INTEGER                :: i, j, k, lbc_len

   REAL                   :: u_at_w, v_at_w

   REAL                   :: u_tend(lenrim(fld_type_u,halo_type_extended),     &
                                    udims_s%k_start:udims_s%k_end),            &
                             v_tend(lenrim(fld_type_v,halo_type_extended),     &
                                    vdims_s%k_start:vdims_s%k_end),            &
                             w_tend(lenrim(fld_type_p,halo_type_extended),     &
                                    wdims_s%k_start:wdims_s%k_end),            &
                             p_tend(lenrim(fld_type_p,halo_type_extended),     &
                                    pdims_s%k_start:pdims_s%k_end+1)
!
   INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
   INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
   REAL(KIND=jprb)               :: zhook_handle

   IF (lhook) CALL dr_hook('init_lbc_dynamics',zhook_in,zhook_handle)

   l_do_halos      = .TRUE.
   l_do_boundaries = .TRUE.

   IF (rim_stepsa == 0) THEN
      inc_fact = 0.0
   ELSE
      inc_fact = 1.0/ (rim_stepsa-MOD(timestep_number-1,rim_stepsa))
   END IF

! update u

   lbc_len = lenrim(fld_type_u,halo_type_extended)

   DO k = 1, model_levels
      DO i = 1, lbc_len
         u_tend(i,k) = u_lbc(i,k) + inc_fact*(u_lbc_tnd(i,k) - u_lbc(i,k))
      END DO
   END DO

! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(row_length, rows, offx, offy,                  &
         model_levels, fld_type_u, u,                                          &
         lenrim(fld_type_u,halo_type_extended),                                &
         lbc_size(1,fld_type_u,halo_type_extended),                            &
         lbc_start(1,fld_type_u,halo_type_extended),                           &
         halo_i, halo_j,                                                       &
         u_tend, rimwidth, rimwidth,  rimweights, at_extremity,                &
         L_do_boundaries, L_do_halos)

! update v

   lbc_len = lenrim(fld_type_v,halo_type_extended)

   DO k = 1, model_levels
      DO i = 1, lbc_len
         v_tend(i,k) = v_lbc(i,k) + inc_fact*(v_lbc_tnd(i,k) - v_lbc(i,k))
      END DO
   END DO

! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(row_length, n_rows, offx, offy,                &
         model_levels, fld_type_v, v,                                          &
         lenrim(fld_type_v,halo_type_extended),                                &
         lbc_size(1,fld_type_v,halo_type_extended),                            &
         lbc_start(1,fld_type_v,halo_type_extended),                           &
         halo_i, halo_j,                                                       &
         v_tend, rimwidth, rimwidth,  rimweights, at_extremity,                &
         L_do_boundaries, L_do_halos)

! update w

   lbc_len = lenrim(fld_type_p,halo_type_extended)

   DO k = 0, model_levels
      DO i = 1, lbc_len
         w_tend(i,k) = w_lbc(i,k) + inc_fact*(w_lbc_tnd(i,k) - w_lbc(i,k))
      END DO
   END DO

! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(row_length, rows, offx, offy,                  &
         model_levels+1, fld_type_p, w,                                        &
         lenrim(fld_type_p,halo_type_extended),                 &
         lbc_size(1,fld_type_p,halo_type_extended),             &
         lbc_start(1,fld_type_p,halo_type_extended),            &
         halo_i, halo_j,                                                       &
         w_tend, rimwidth, rimwidth,  rimweights, at_extremity,                &
         L_do_boundaries, L_do_halos)


! update rho
! cannot do rho here because of r^2 tern in tendancy array
!  DO k = 1, model_levels
!     DO i = 1, lbc_len
!        p_tend(i,k) = rho_lbc(i,k) + inc_fact*(rho_lbc_tnd(i,k) - rho_lbc(i,k))
!     END DO
!  END DO

! DEPENDS ON: set_lateral_boundaries
!   CALL set_lateral_boundaries(row_length, rows, offx, offy,                  &
!        model_levels, fld_type_p, rho,                                        &
!        lenrim(fld_type_p,halo_type_extended),                                &
!        lbc_size(1,fld_type_p,halo_type_extended),                            &
!        lbc_start(1,fld_type_p,halo_type_extended),                           &
!        halo_i, halo_j,                                                       &
!        p_tend, rimwidth, rimwidth,  rimweights, at_extremity,                &
!        L_do_boundaries, L_do_halos)

! update Exner

   DO k = 1, model_levels+1
      DO i = 1, lbc_len
         p_tend(i,k) = exner_lbc(i,k)                                          &
                        + inc_fact*(exner_lbc_tnd(i,k) - exner_lbc(i,k))
      END DO
   END DO

! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(row_length, rows, offx, offy,                  &
         model_levels+1, fld_type_p, exner,                                    &
         lenrim(fld_type_p,halo_type_extended),                                &
         lbc_size(1,fld_type_p,halo_type_extended),                            &
         lbc_start(1,fld_type_p,halo_type_extended),                           &
         halo_i, halo_j,                                                       &
         p_tend, rimwidth, rimwidth,  rimweights, at_extremity,                &
         L_do_boundaries, L_do_halos)

! update thetav

   DO k = 1, model_levels+1
      DO i = 1, lbc_len
         w_tend(i,k) = theta_lbc(i,k)                                          &
                        + inc_fact*(theta_lbc_tnd(i,k) - theta_lbc(i,k))
      END DO
   END DO

! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(row_length, rows, offx, offy,                  &
         model_levels+1, fld_type_p, theta,                                    &
         lenrim(fld_type_p,halo_type_extended),                                &
         lbc_size(1,fld_type_p,halo_type_extended),                            &
         lbc_start(1,fld_type_p,halo_type_extended),                           &
         halo_i, halo_j,                                                       &
         w_tend, rimwidth, rimwidth,  rimweights, at_extremity,                &
         L_do_boundaries, L_do_halos)

! 3 main Moisture fields

   DO k = 1, model_levels+1
      DO i = 1, lbc_len
         w_tend(i,k) = mv_lbc(i,k) + inc_fact*(mv_lbc_tnd(i,k) - mv_lbc(i,k))
      END DO
   END DO

! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(row_length, rows, offx, offy,                  &
         model_levels+1, fld_type_p, m_v,                                      &
         lenrim(fld_type_p,halo_type_extended),                                &
         lbc_size(1,fld_type_p,halo_type_extended),                            &
         lbc_start(1,fld_type_p,halo_type_extended),                           &
         halo_i, halo_j,                                                       &
         w_tend, rimwidth, rimwidth,  rimweights, at_extremity,                &
         L_do_boundaries, L_do_halos)

   DO k = 1, model_levels+1
      DO i = 1, lbc_len
         w_tend(i,k) = mcl_lbc(i,k) + inc_fact*(mcl_lbc_tnd(i,k) - mcl_lbc(i,k))
      END DO
   END DO

! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(row_length, rows, offx, offy,                  &
         model_levels+1, fld_type_p, m_cl,                                     &
         lenrim(fld_type_p,halo_type_extended),                                &
         lbc_size(1,fld_type_p,halo_type_extended),                            &
         lbc_start(1,fld_type_p,halo_type_extended),                           &
         halo_i, halo_j,                                                       &
         w_tend, rimwidth, rimwidth,  rimweights, at_extremity,                &
         L_do_boundaries, L_do_halos)

   DO k = 1, model_levels+1
      DO i = 1, lbc_len
         w_tend(i,k) = mcf_lbc(i,k) + inc_fact*(mcf_lbc_tnd(i,k) - mcf_lbc(i,k))
      END DO
   END DO

! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(row_length, rows, offx, offy,                  &
         model_levels+1, fld_type_p, m_cf,                                     &
         lenrim(fld_type_p,halo_type_extended),                                &
         lbc_size(1,fld_type_p,halo_type_extended),                            &
         lbc_start(1,fld_type_p,halo_type_extended),                           &
         halo_i, halo_j,                                                       &
         w_tend, rimwidth, rimwidth,  rimweights, at_extremity,                &
         L_do_boundaries, L_do_halos)

! rain

   IF( l_mcr_qrain ) THEN
      DO k = 1, model_levels+1
         DO i = 1, lbc_len
            w_tend(i,k) = mr_lbc(i,k)                                          &
                          + inc_fact*(mr_lbc_tnd(i,k) - mr_lbc(i,k))
         END DO
      END DO

! DEPENDS ON: set_lateral_boundaries
       CALL set_lateral_boundaries(row_length, rows, offx, offy,               &
            model_levels+1, fld_type_p, m_r,                                   &
            lenrim(fld_type_p,halo_type_extended),                             &
            lbc_size(1,fld_type_p,halo_type_extended),                         &
            lbc_start(1,fld_type_p,halo_type_extended),                        &
            halo_i, halo_j,                                                    &
            w_tend, rimwidth, rimwidth,  rimweights, at_extremity,             &
            L_do_boundaries, L_do_halos)
   END IF

! graup

   IF( l_mcr_qgraup ) THEN
      DO k = 1, model_levels+1
         DO i = 1, lbc_len
            w_tend(i,k) = mgr_lbc(i,k)                                         &
                          + inc_fact*(mgr_lbc_tnd(i,k) - mgr_lbc(i,k))
         END DO
      END DO

! DEPENDS ON: set_lateral_boundaries
       CALL set_lateral_boundaries(row_length, rows, offx, offy,               &
            model_levels+1, fld_type_p, m_gr,                                  &
            lenrim(fld_type_p,halo_type_extended),                             &
            lbc_size(1,fld_type_p,halo_type_extended),                         &
            lbc_start(1,fld_type_p,halo_type_extended),                        &
            halo_i, halo_j,                                                    &
            w_tend, rimwidth, rimwidth,  rimweights, at_extremity,             &
            L_do_boundaries, L_do_halos)
   END IF

! cf2

   IF( l_mcr_qcf2 ) THEN
      DO k = 1, model_levels+1
         DO i = 1, lbc_len
            w_tend(i,k) = mcf2_lbc(i,k)                                        &
                          + inc_fact*(mcf2_lbc_tnd(i,k) - mcf2_lbc(i,k))
         END DO
      END DO

! DEPENDS ON: set_lateral_boundaries
       CALL set_lateral_boundaries(row_length, rows, offx, offy,               &
            model_levels+1, fld_type_p, m_cf2,                                 &
            lenrim(fld_type_p,halo_type_extended),                             &
            lbc_size(1,fld_type_p,halo_type_extended),                         &
            lbc_start(1,fld_type_p,halo_type_extended),                        &
            halo_i, halo_j,                                                    &
            w_tend, rimwidth, rimwidth,  rimweights, at_extremity,             &
            L_do_boundaries, L_do_halos)
   END IF

! Fix etadot

  DO k = 1, model_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

          u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +              &
                                     intw_u2p(i,2)*u(i,j,k+1) ) +              &
                   intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k)   +              &
                                     intw_u2p(i,2)*u(i,j,k) )

          v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +              &
                                     intw_v2p(j,2)*v(i,j,k+1) ) +              &
                   intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +                &
                                     intw_v2p(j,2)*v(i,j,k) )

          etadot(i,j,k) = ( w(i,j,k)/h3_p_eta(i,j,k) -                         &
                             u_at_w*dxi1_xi3(i,j,k)/                           &
                                           h1_p_eta(i,j,k) -                   &
                             v_at_w*dxi2_xi3(i,j,k)/                           &
                                           h2_p_eta(i,j,k) ) /                 &
                                             deta_xi3_theta(i,j,k)
      END DO
    END DO
  END DO

  etadot(:,:,0) = 0.0
  etadot(:,:,model_levels) = 0.0

! DEPENDS ON: swap_bounds
   CALL swap_bounds(etadot, row_length, rows, model_levels+1,                  &
                    offx, offy, fld_type_p, .FALSE.)

   IF (lhook) CALL dr_hook('init_lbc_dynamics',zhook_out,zhook_handle)

   END SUBROUTINE init_lbc_dynamics
