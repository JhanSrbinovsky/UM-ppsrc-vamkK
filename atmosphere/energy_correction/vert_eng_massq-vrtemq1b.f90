! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE VERT_ENG_MASSQ
!  
!    PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!              Given the prognostic input fields it integrates
!              vertically various fields required by the energy
!              correction.
!              Also if call from section 30  (climate diagnostics)
!              integrates qu fluxes etc
!  
!  
!    PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!  
!    DOCUMENTATION :  Energy correction documentation
!  
!----------------------------------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Energy Correction

      MODULE vert_eng_massq_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE VERT_ENG_MASSQ(                                        &
                            wet_model_levels,theta ,u,v,w, rho_r2, q,   &
                            qcl,qcf,                                    &
                            wet_to_dry_n,                               &
                            exner_theta_levels,                         &
                            Lqflux,                                     &
! output fields
                            rho_dry,dry_to_wet,                         &
                            vert_int_array,vert_qflux)

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE atm_fields_bounds_mod, ONLY: udims, vdims, wdims, tdims,      &
                                       pdims, udims_s, vdims_s, tdims_s,&
                                       pdims_s, wdims_s, qdims_l,       &
                                       qdims_s
          
      USE level_heights_mod,     ONLY: r_theta_levels, r_rho_levels     
          
      USE swapable_field_mod,    ONLY: swapable_field_pointer_type

      USE earth_constants_mod,   ONLY: g, earth_radius
      USE atmos_constants_mod,   ONLY: r, cp

      USE yomhook,               ONLY: lhook, dr_hook
      USE parkind1,              ONLY: jprb, jpim
      USE Field_Types,           ONLY: fld_type_u,fld_type_v

      USE eng_mass_param_mod,    ONLY: ip_dry_mass, ip_wet_mass, ip_cvT,&
                                       ip_gr,       ip_keu,      ip_kev,&
                                       ip_kew,      ip_q,        ip_qcl,&
                                       ip_qcf,      ip_qu,       ip_qv, &
                                       ip_qw,       n_sums,      n_flux

      USE proc_info_mod,         ONLY: at_extremity,model_domain,       &
                                       proc_row_group=>                 &
                                          gc_proc_row_group,            &
                                       global_row_length

      USE trignometric_mod,      ONLY : cos_theta_longitude,            &
                                        sin_theta_longitude
      USE um_parparams,          ONLY : psouth, pnorth

      USE u_to_p_mod, ONLY: u_to_p
      USE v_to_p_mod, ONLY: v_to_p
      IMPLICIT NONE

! constants required for calculations
! Required for interpolation

!----------------------------------------------------------------------
! INPUT variables
!----------------------------------------------------------------------

      REAL, INTENT(IN) ::                                               &
        theta(tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end)
                                              
      REAL, target ::                                                   &
            U(udims_s%i_start:udims_s%i_end,                            &
              udims_s%j_start:udims_s%j_end,                            &
              udims_s%k_start:udims_s%k_end),                           &
                                                
            V(vdims_s%i_start:vdims_s%i_end,                            &
              vdims_s%j_start:vdims_s%j_end,                            &
              vdims_s%k_start:vdims_s%k_end)
                                                
      REAL, INTENT(IN) ::                                               &
            w(wdims_s%i_start:wdims_s%i_end,                            &
              wdims_s%j_start:wdims_s%j_end,                            &
              wdims_s%k_start:wdims_s%k_end),                           &
            q(qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            & 
              qdims_l%k_start:qdims_l%k_end),                           &
          qcl(qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end),                           &
          qcf(qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end),                           &
 wet_to_dry_n(tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end),                           &
        exner_theta_levels(tdims_s%i_start:tdims_s%i_end,               &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end),                           &
       rho_r2(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end)

      LOGICAL, INTENT (IN) ::                                           &
        lqflux            ! true if q flux calculations required.

!   Vertical integrals plus info on rho_dry

      REAL, INTENT(INOUT)  ::                                           &
        vert_int_array(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end,n_sums),               &
                                              ! vertical integrals
        vert_qflux(pdims%i_start:pdims%i_end,                           &
                       pdims%j_start:pdims%j_end,n_flux),               &
                                              ! q flux integrals
        rho_dry(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,pdims%k_start:pdims%k_end),   &
                                                ! rho dry x r^2
        dry_to_wet(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end)

!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------
! pointers for sum_array

      REAL :: weight1, weight2, weight3,                                &
              tempd, tempw, ww2, ww1,                                   &
              cv                        ! specific heat at constant vol

      REAL ::                                                           &
        rho(       pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
                                                ! rho only
        delr_rho(  pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
                                                ! dr
        T(         pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
                                                ! TEMPERATURE
        rho_theta( pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
                                                 ! not used
             t_rho(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
             u_rho(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
             v_rho(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
             w_rho(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
             q_rho(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
           qcl_rho(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end),                          &
           qcf_rho(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end)

      REAL ::                                                           &
        mag_vector_np(pdims%k_end),                                     &
                                     ! magnitude of the vector wind NP
        dir_vector_np(pdims%k_end),                                     &
                                     ! direction of the vector wind NP
        mag_vector_sp(pdims%k_end),                                     &
                                     ! magnitude of the vector wind SP
        dir_vector_sp(pdims%k_end)   ! direction of the vector wind SP

      INTEGER :: i, j, k, n             ! LOOP COUNTER

      INTEGER :: i_field                ! counter for swapbounds

      TYPE(swapable_field_pointer_type) :: fields_to_swap(2)

      INTEGER :: wet_model_levels

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('VERT_ENG_MASSQ',zhook_in,zhook_handle)

      cv=cp-R       ! value for dry air

!----------------------------------------------------------------------
! zero output arrays
!----------------------------------------------------------------------
      DO n=1,n_sums         ! energy integrals
        DO j=pdims%j_start,pdims%j_end
          DO i=pdims%i_start,pdims%i_end
            vert_int_array(i,j,n) = 0.0
          END DO
        END DO
      END DO
      IF (lqflux) THEN
        DO n=1,n_flux          ! q flux integrals
          DO j=pdims%j_start,pdims%j_end
            DO i=pdims%i_start,pdims%i_end
              vert_qflux(i,j,n) = 0.0
            END DO
          END DO
        END DO
      END IF

!----------------------------------------------------------------------
! First convert theta to temperature
!----------------------------------------------------------------------
      DO k = pdims%k_start,pdims%k_end
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start,pdims%i_end
            T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
          END DO
        END DO
      END DO

!----------------------------------------------------------------------
! CALCULATE layer thickness for rho layers and calculate rho
!----------------------------------------------------------------------
      DO K=pdims%k_start,pdims%k_end
        DO j =pdims%j_start,pdims%j_end
          DO I=pdims%i_start,pdims%i_end

            DELR_RHO(I,J,K) = R_THETA_LEVELS(I,J,K) -                   &
     &                        R_THETA_LEVELS(I,J,K-1)
            rho(i,j,k) = rho_r2(i,j,k)/                                 &
     &           (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))

          END DO
        END DO
      END DO

!----------------------------------------------------------------------
! Convert rho to rho dry in the same way as done in dynamics - Flux_rho
! Note now uses linear vertical intepolation to be consistent with
! new dynamics code at UM 5.1.
! Now using the wet_to_dry_n array from atm_step. Doing this ensures
! consistency with dynamics conversions and means the code is also correct
! if more moist prognostics e.g. rain and graupel are added.
!----------------------------------------------------------------------

      DO k = pdims%k_start,pdims%k_end
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start,pdims%i_end
            dry_to_wet(i,j,k)= 1./wet_to_dry_n(i,j,k)
            rho_dry(i,j,k) = rho_r2(i,j,k) * wet_to_dry_n(i,j,k)
          END DO
        END DO
      END DO

!-------------------------------------------------------------------
! Intepolate T, w & q to rho points
! Note needs to remain consistent with intepolation used in dynamics.
! Using linear interpolation
!
!                      K               for rho
!      K                          K-1  for theta
!      X<--- w1------->X<-- w2--->X
!       <----------w3------------>
!-------------------------------------------------------------------
      k=1
      DO j = pdims%j_start,pdims%j_end
        DO i=pdims%i_start,pdims%i_end
! assume bottom rho level value equal to bottom theta level value
          q_rho(i,j,k)   = q(i,j,k)
          qcl_rho(i,j,k) = qcl(i,j,k)
          qcf_rho(i,j,k) = qcf(i,j,k)
          T_rho(i,j,k)   = T(i,j,k)
! only w has a value at surface
          weight1 = r_theta_levels(i,j,k) - r_rho_levels  (i,j,k)
          weight2 = r_rho_levels  (i,j,k) - r_theta_levels(i,j,k-1)
          weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
          ww1 = weight1/weight3
          ww2 = weight2/weight3
          w_rho(i,j, K) =  ww2 * w(i,j,K)  + ww1 * w(i,j,K-1)
        END DO
      END DO

      DO K=2,wet_model_levels
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start,pdims%i_end
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            ww2 = weight2/weight3
            t_rho(i,j, K) = ww2 * t(i,j,K) + ww1 * t(i,j,K-1)
            w_rho(i,j, K) = ww2 * w(i,j,K) + ww1 * w(i,j,K-1)
            q_rho (i,j, K)  = ww2 * q(i,j,K)   + ww1 * q(i,j,K-1)
            qcl_rho(i,j, K) = ww2 * qcl(i,j,K) + ww1 * qcl(i,j,K-1)
            qcf_rho(i,j, K) = ww2 * qcf(i,j,K) + ww1 * qcf(i,j,K-1)
          END DO
        END DO
      END DO

! Special case of wet_model_levels < model_levels
! wet to dry conversions imply q_rho has a value on wet_model_levels+1
! Assume also true for qcl and qcf ?

      IF (wet_model_levels < pdims%k_end) THEN
        k=wet_model_levels+1
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start,pdims%i_end
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            q_rho (i,j, K)  = ww1 * q(i,j,K-1)
            qcl_rho(i,j, K) = ww1 * qcl(i,j,K-1)
            qcf_rho(i,j, K) = ww1 * qcf(i,j,K-1)
          END DO
        END DO
      END IF

! Any dry levels at top of model
      DO K=wet_model_levels+1,pdims%k_end
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start,pdims%i_end
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            ww2 = weight2/weight3
            t_rho(i,j, K) = ww2 * t(i,j,K) + ww1 * t(i,j,K-1)
            w_rho(i,j, K) = ww2 * w(i,j,K) + ww1 * w(i,j,K-1)
          END DO
        END DO
      END DO
!-------------------------------------------------------------------
! Interpolate winds to rho points (already on same vertical level)
!-------------------------------------------------------------------
! Need to call swap bounds as halo points not set
      i_field = 0
      i_field = i_field + 1
      fields_to_swap(i_field) % field      => u(:,:,:)
      fields_to_swap(i_field) % field_type = fld_type_u
      fields_to_swap(i_field) % levels     = udims%k_end
      fields_to_swap(i_field) % rows       = udims%j_end-udims%j_start+1
      fields_to_swap(i_field) % vector     = .TRUE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field      => v(:,:,:)
      fields_to_swap(i_field) % field_type = fld_type_v
      fields_to_swap(i_field) % levels     = vdims%k_end
      fields_to_swap(i_field) % rows       = vdims%j_end-vdims%j_start+1
      fields_to_swap(i_field) % vector     = .TRUE.

! DEPENDS ON: swap_bounds_mv
      CALL swap_bounds_mv(fields_to_swap, i_field,                      &
                          udims%i_end-udims%i_start+1,                  &
                          udims_s%halo_i,udims_s%halo_j)

      CALL v_to_p(v,                                                    &
                        vdims_s%i_start,vdims_s%i_end,                  & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        pdims%k_end,                                    &
                        model_domain,at_extremity,v_rho)

      CALL u_to_p(u,                                                    &
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        pdims%k_end,                                    &
                        model_domain,at_extremity,u_rho)

! Problem of u & v values at poles
! Uses v on row next to poles to calculate wind at the polar point

IF (.NOT. l_vatpoles) THEN

        IF (model_domain  ==  1 ) THEN    ! global model

! DEPENDS ON: polar_vector_wind_n
          Call Polar_vector_wind_n                                      &
                            (v,                                         &
                             sin_theta_longitude,cos_theta_longitude,   &
                             vdims%i_end-vdims%i_start+1,               &
                             vdims%j_end-vdims%j_start+1, pdims%k_end,  &
                             mag_vector_np,dir_vector_np,               &
                             mag_vector_sp,dir_vector_sp,               &
                             vdims_s%halo_i, vdims_s%halo_j,            &
                             global_row_length,                         &
                             proc_row_group, at_extremity)

          IF (at_extremity(PSouth) ) THEN
            DO k = pdims%k_start,pdims%k_end
              DO i = pdims%i_start,pdims%i_end
                u_rho(i,pdims%j_start,k) = 0.0
                v_rho(i,pdims%j_start,k) = mag_vector_sp(k)
              END DO
            END DO
          END IF
         IF (at_extremity(PNorth) ) THEN
            DO k = pdims%k_start,pdims%k_end
              DO i = pdims%i_start,pdims%i_end
                u_rho(i,pdims%j_end,k) = 0.0
                v_rho(i,pdims%j_end,k) = mag_vector_np(k)
              END DO
            END DO
          END IF
        END IF              ! end test on global model
END IF ! vatpoles
!----------------------------------------------------------------------
! Integrals over fields using values interpolated to rho grid
! At present all integrals done after interpolation to rho points.
! Integrals for energy involve using dry mass. This may need to be
! altered.
!----------------------------------------------------------------------
! (a) fields on all model levels

      DO k = pdims%k_start,pdims%k_end
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start,pdims%i_end

            tempw = RHO_r2(I,j,k)*DELR_RHO(I,j,k)   ! wet mass

            tempd = RHO_dry(I,j,k)*DELR_RHO(I,j,k)  ! dry mass

! dry mass
            vert_int_array(I,J,ip_dry_mass) =                           &
                          vert_int_array(i,j,ip_dry_mass) + tempd

! total mass wet + dry
            vert_int_array(I,J,ip_wet_mass) =                           &
                          vert_int_array(i,j,ip_wet_mass) + tempw

! gr term
            vert_int_array(I,J,ip_gr) = vert_int_array(i,j,ip_gr) +     &
               g*(r_rho_levels(i,j,k)-earth_radius)*tempd

! KE terms for u, v & w
            vert_int_array(I,J,ip_keu) = vert_int_array(i,j,ip_keu) +   &
                       0.5*u_rho(i,j,k)*u_rho(i,j,k)*tempd
            vert_int_array(I,J,ip_kev) = vert_int_array(i,j,ip_kev) +   &
                       0.5*v_rho(i,j,k)*v_rho(i,j,k)*tempd
            vert_int_array(I,J,ip_kew) = vert_int_array(i,j,ip_kew) +   &
                       0.5*w_rho(i,j,k)*w_rho(i,j,k)*tempd

! cvT term
            vert_int_array(I,J,ip_cvt) = vert_int_array(i,j,ip_cvt) +   &
                                  cv*T_rho(i,j,k)*tempd
          END DO
        END DO
      END DO


      DO k=pdims%k_start,wet_model_levels
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start,pdims%i_end
            tempw = RHO_r2(I,j,k)*DELR_RHO(I,j,k)  ! wet mass
! q*rho
            vert_int_array(I,J,ip_q) = vert_int_array(i,j,ip_q) +       &
                                            q_rho(i,j,k)*tempw
! qcl*rho
            vert_int_array(I,J,ip_qcl) = vert_int_array(i,j,ip_qcl) +   &
                                            qcl_rho(i,j,k)*tempw
! qcf*rho
            vert_int_array(I,J,ip_qcf) = vert_int_array(i,j,ip_qcf) +   &
                                            qcf_rho(i,j,k)*tempw

          END DO
        END DO
      END DO

! Special case of wet_model_levels < model_levels
      IF (wet_model_levels < pdims%k_end) THEN
        k=wet_model_levels+1
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start,pdims%i_end
            tempw = RHO_r2(I,j,k)*DELR_RHO(I,j,k)  ! wet mass
            vert_int_array(I,J,ip_q) = vert_int_array(i,j,ip_q) +       &
                                            q_rho(i,j,k)*tempw
            vert_int_array(I,J,ip_qcl) = vert_int_array(i,j,ip_qcl) +   &
                                            qcl_rho(i,j,k)*tempw
            vert_int_array(I,J,ip_qcf) = vert_int_array(i,j,ip_qcf) +   &
                                            qcf_rho(i,j,k)*tempw
          END DO
        END DO
      END IF   !         wet_model_levels < model_levels

!----------------------------------------------------------------------
! Additional q flux integrals not required by energy correction
! At present all integrals done after interpolation to rho points.
!----------------------------------------------------------------------


      IF (Lqflux) THEN
        DO K = pdims%k_start,wet_model_levels
          DO j = pdims%j_start,pdims%j_end
            DO I = pdims%i_start,pdims%i_end
              tempw = RHO_r2(I,j,k)*DELR_RHO(I,j,k)
              vert_qflux(I,J,ip_qu) = vert_qflux(i,j,ip_qu) +           &
                                     u_rho(i,j,k)*q_rho(i,j,k)*tempw
              vert_qflux(I,J,ip_qv) = vert_qflux(i,j,ip_qv) +           &
                                     v_rho(i,j,k)*q_rho(i,j,k)*tempw
              vert_qflux(I,J,ip_qw) = vert_qflux(i,j,ip_qw) +           &
                                     w_rho(i,j,k)*q_rho(i,j,k)*tempw

            END DO
          END DO
        END DO

! Special case of wet_model_levels < model_levels

        IF (wet_model_levels < pdims%k_end) THEN
          k=wet_model_levels+1
          DO j = pdims%j_start,pdims%j_end
            DO i = pdims%i_start,pdims%i_end
              vert_qflux(I,J,ip_qu) = vert_qflux(i,j,ip_qu) +           &
                                     u_rho(i,j,k)*q_rho(i,j,k)*tempw
              vert_qflux(I,J,ip_qv) = vert_qflux(i,j,ip_qv) +           &
                                     v_rho(i,j,k)*q_rho(i,j,k)*tempw
              vert_qflux(I,J,ip_qw) = vert_qflux(i,j,ip_qw) +           &
                                     w_rho(i,j,k)*q_rho(i,j,k)*tempw
            END DO
          END DO
        END IF  !           wet_model_levels < model_levels
      END IF
!----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('VERT_ENG_MASSQ',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE VERT_ENG_MASSQ
      END MODULE
