! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE ENG_MASS_DIAG------------------------------------------
!  
!    PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!              - TO GLOBALLY INTERGATE TOTAL ENERGY AMD MASS OF
!                THE ATMOSPHERE
!  
!    NOT SUITABLE FOR SINGLE COLUMN MODEL USE
!  
!    PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!  
!    DOCUMENTATION :
!  
!----------------------------------------------------------------------
!
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Energy Correction

      MODULE eng_mass_diag_mod

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE eng_mass_diag(                                         &
                            wet_model_levels,                           &
                            delta_lambda,delta_phi,                     & 
                            theta ,u,v,w, rho_r2,q,qcl,qcf,             &
                            wet_to_dry_n,                               &
                            exner_theta_levels,                         &
                            sum_moist_flux,                             &
                            tot_mass_init,tot_m_init,                   &
                            Lmass_corr,lqt_corr,lemq_print,             &
                            a_energysteps,timestep,tot_energy_final,    &
                            tot_mass_final,tot_m_final)

      USE dynamics_input_mod, ONLY: l_endgame

      USE water_constants_mod, ONLY: lc, lf
      USE global_2d_sums_mod, ONLY: global_2d_sums
      
      USE level_heights_mod, ONLY:  r_theta_levels, r_rho_levels 
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE um_parparams,          ONLY: NoDomain,PNorth,PEast,PSouth,    &
                                       PWest
      USE proc_info_mod,       ONLY: at_extremity,global_row_length,    &     
                                     gc_proc_row_group,                 &
                                     gc_proc_col_group,                 &
                                     model_domain,mype=>me , n_proc,    &
                                     n_procx, n_procy

      USE atm_fields_bounds_mod, ONLY: tdims_s,udims_s,vdims_s,         &
                                       wdims_s,qdims_l,tdims_s,         &
                                       pdims_s,pdims  ,pdims_s
      USE trignometric_mod,    ONLY: FV_cos_theta_latitude,             &
                                     cos_v_latitude,                    &
                                     cos_theta_longitude,               &
                                     sin_theta_longitude
 
      USE vert_eng_massq_mod,  ONLY: vert_eng_massq
      USE eng_mass_param_mod,  ONLY: ip_dry_mass, ip_wet_mass, ip_cvT,  &
                                     ip_gr,       ip_keu,      ip_kev,  &
                                     ip_kew,      ip_q,        ip_qcl,  &
                                     ip_qcf,      ip_qu,       ip_qv,   &
                                     ip_qw,       n_sums,      n_flux

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: wet_model_levels
      REAL delta_lambda,delta_phi 

!Input model data
      REAL, INTENT (IN) ::                                              &
        theta(tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end),                           &
            U(udims_s%i_start:udims_s%i_end,                            &
              udims_s%j_start:udims_s%j_end,                            &
              udims_s%k_start:udims_s%k_end),                           &
            V(vdims_s%i_start:vdims_s%i_end,                            &
              vdims_s%j_start:vdims_s%j_end,                            &
              vdims_s%k_start:vdims_s%k_end),                           &
            w(wdims_s%i_start:wdims_s%i_end,                            &
              wdims_s%j_start:wdims_s%j_end,                            &
              wdims_s%k_start:wdims_s%k_end),                           &
        exner_theta_levels(tdims_s%i_start:tdims_s%i_end,               &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end)


!     these cannot be INTENT IN because they are passed on to
!     subroutines, however they are input fields
      REAL, INTENT (INOUT) ::                                           &
  wet_to_dry_n(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
            q(qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            & 
              qdims_l%k_start:qdims_l%k_end),                           &
          qcl(qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end),                           &
          qcf(qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end),                           &
        sum_moist_flux(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end)

! In/OUT  model data  Altered in Lmass_corr set
      REAL :: rho_r2(pdims_s%i_start:pdims_s%i_end,                     &
                   pdims_s%j_start:pdims_s%j_end,                       &
                   pdims_s%k_start:pdims_s%k_end)


      LOGICAL ::                                                        &
        Lmass_corr,                                                     &
                        ! if true do a mass correction
        lqt_corr,                                                       &
                        ! if true do a moisture correction
        lemq_print      ! if true print addtional info


      REAL ::                                                           &
           tot_energy_final,                                            &
                             ! TOTAL ENERGY OF ATMOSPHERE
           tot_mass_final,                                              &
                             ! total dry mass of atmosphere
           tot_m_final       ! total moisture

      INTEGER, INTENT (IN) ::                                           &
           a_energysteps      ! number of timesteps between energy cor
      REAL, INTENT (IN) ::                                              &
           tot_mass_init,                                               &
                              ! total dry mass of atmosphere at T+0
           tot_m_init,                                                  &
                              ! total moisture at previous time
           timestep           ! tiemstep of model in seconds.

!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------

      REAL ::                                                           &
       sum_array(pdims%i_start:pdims%i_end,                             &
                 pdims%j_start:pdims%j_end,n_sums),                     &
                                            ! array to be summed
       dummy_array(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,n_flux),                   &
                                            ! dummy array
       sum_results(n_sums),                                             &
                                            ! sum of SUM_ARRAY
       moist_in(1)                     ! Array for global_sums

      REAL ::                                                           &
        cv,                                                             &
                                       ! specific heat at constant vol
        tot_dry_energy,                                                 &
        tot_wet_mass,                                                   &
        tot_dry_mass,                                                   &
        tot_ke,                                                         &
        err_mass,                                                       &
                                       ! error in dry mass
        mass_corr,                                                      &
                                       ! correction factor dry mass
        factor,                                                         &
                                       ! grid resolution factor
        tot_moist_in,                                                   &
        err_moist,                                                      &
        moist_corr,                                                     &
                                       ! moist correction
        change_moist,                                                   &
                                       ! change in qt over period
        per_err,                                                        &
        q_cld,                                                          &
                                       ! total cloud water
        weight1, weight2 , weight3,                                     &
                                       ! weights
        tempw, tempd

      REAL ::                                                           &
        rho_dry(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end),                             &
        dry_to_wet(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   pdims%k_start:pdims%k_end)

      INTEGER :: i, j, k , n               ! LOOP COUNTER

      LOGICAL :: lqflux         ! switch for calculations
      INTEGER :: istat
      REAL    :: sum_rows(pdims%j_end)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!----------------------------------------------------------------------
! zero summation arrays
!----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('ENG_MASS_DIAG',zhook_in,zhook_handle)

      DO n=1,n_sums
        sum_results(n)=0.0
      END DO

! set logical so unnecessary calculation not done

      Lqflux = .false.

!----------------------------------------------------------------------
!   call vert_eng_massq to do vertical integrations
!----------------------------------------------------------------------

      CALL vert_eng_massq(                                              &
                            wet_model_levels,                           &
                            theta ,u,v,w, rho_r2, q,qcl,qcf             &
                           ,wet_to_dry_n                                &
                           ,exner_theta_levels                          &
                           ,Lqflux                                      &
                           ,rho_dry,dry_to_wet                          &
                           ,sum_array, dummy_array )

!----------------------------------------------------------------------
!  multiply all vertical integrals by cos(latitude)
!----------------------------------------------------------------------


      DO n = 1,n_sums
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start,pdims%i_end
             sum_array(i,j,n) = sum_array(i,j,n)*                       &
                                FV_cos_theta_latitude(i,j)
          END DO
        END DO
      END DO

      IF ( .NOT. l_endgame ) THEN
!----------------------------------------------------------------------
! reset sum_array so only one contribution from each pole
! NOTE : FV_cos_theta_latitude at poles holds area for whole pole
!----------------------------------------------------------------------
      IF (at_extremity(PSouth)) THEN
        IF (at_extremity(PWest)) THEN
          DO k=1,n_sums
            DO i=2,pdims%i_end
              sum_array(i,1,k)=0.0
            END DO
          END DO
        else
          DO k=1,n_sums
            DO i= pdims%i_start,pdims%i_end
              sum_array(i,1,k)=0.0
            END DO
          END DO
        END IF
      END IF
      IF (at_extremity(PNorth)) THEN
        IF (at_extremity(PWest)) THEN
          DO k=1,n_sums
            DO i=2,pdims%i_end
              sum_array(i,pdims%j_end,k)=0.0
            END DO
          END DO
        else
          DO k=1,n_sums
            DO i=pdims%i_start,pdims%i_end
              sum_array(i,pdims%j_end,k)=0.0
            END DO
          END DO
        END IF
      END IF
      END IF

!----------------------------------------------------------------------
! Section 3 - Calculation of energy
! Work out the global sums and multiply by correct factor for grid res.
!----------------------------------------------------------------------
      factor=delta_lambda*delta_phi   ! already in radians

! Do the global sums
! Note no halo on input array, grid_type 1, halo_type 3

      CALL global_2d_sums(sum_array,pdims%i_end,pdims%j_end,0,0,n_sums, &
                                     sum_results)


! total dry air energy minus latent heat  held by cloud

      tot_dry_energy =( sum_results(ip_cvT)+ sum_results(ip_gr)         &
                 + sum_results(ip_keu) + sum_results(ip_kev)            &
                 + sum_results(ip_kew) )*factor

      tot_energy_final  = tot_dry_energy -                              &
       (lc*sum_results(ip_qcl)+(lc+lf)*sum_results(ip_qcf) )*factor

      tot_dry_mass =sum_results(ip_dry_mass)*factor
      tot_mass_final = tot_dry_mass
      tot_wet_mass =sum_results(ip_wet_mass)*factor

! Print out info if required

      IF (mype == 0.AND.lemq_print) THEN

        tot_ke= (sum_results(ip_keu) + sum_results(ip_kev)              &
                 + sum_results(ip_kew) )*factor

        WRITE(6,fmt='(A,E15.5)')'Tot dry mass       ',                  &
               tot_dry_mass
        WRITE(6,fmt='(A,E15.5)')'Tot mass           ',                  &
               tot_wet_mass
        WRITE(6,fmt='(A,E15.5)')'Tot energy         ',                  &
               tot_energy_final
        WRITE(6,fmt='(A,E15.5)')'tot dry energy     ',                  &
               tot_dry_energy

        WRITE(6,fmt='(A,E15.5)')'gr( rho cal)       ',                  &
               sum_results(ip_gr)*factor
        WRITE(6,fmt='(A,E15.5)')'KE( rho cal)       ',                  &
               tot_ke
        WRITE(6,fmt='(A,E15.5)')'KEu(rho cal)       ',                  &
               sum_results(ip_keu)*factor
        WRITE(6,fmt='(A,E15.5)')'KEv(rho cal)       ',                  &
               sum_results(ip_kev)*factor
        WRITE(6,fmt='(A,E15.5)')'KEw(rho cal)       ',                  &
               sum_results(ip_kew)*factor

        WRITE(6,fmt='(A,E15.5)')'cvT( rho cal)      ',                  &
               sum_results(ip_cvt)*factor
        WRITE(6,fmt='(A,E15.5)')'lq ( rho cal)      ',                  &
               sum_results(ip_q)*factor*lc
        WRITE(6,fmt='(A,E15.5)')'lqcf( rho cal)     ',                  &
               sum_results(ip_qcl)*factor*lc
        WRITE(6,fmt='(A,E15.5)')'lqcl( rho cal)     ',                  &
               sum_results(ip_qcf)*factor*(lc+lf)

      END IF

!----------------------------------------------------------------------
! Section 4. Mass conservation & correction
! Note this is the best place to do this as rho_dry has already been
! evaluated earlier in subroutine
!----------------------------------------------------------------------
      IF (Lemq_print.or.Lmass_corr) THEN
!----------------------------------------------------------------------
! Work out drift in dry mass and printout results
!----------------------------------------------------------------------

        err_mass  = tot_mass_final - tot_mass_init
        mass_corr = tot_mass_init/tot_mass_final

        IF (mype == 0) THEN
          WRITE(6,fmt='(A,E20.5,A)') 'Final dry mass of atmosphere  = ',&
                      tot_mass_final, ' KG'
          WRITE(6,fmt='(A,E20.5,A)') 'Initial dry mass of atmosphere= ',&
                      tot_mass_init, ' KG'
          WRITE(6,fmt='(A,E20.5  )') 'Correction factor for rho_dry = ',&
                      mass_corr
        END IF

      END IF

!----------------------------------------------------------------------
! Apply dry mass correction
!----------------------------------------------------------------------
      IF (lmass_corr) THEN

        DO K = pdims%k_start,pdims%k_end
          DO j = pdims%j_start,pdims%j_end
            DO I = pdims%i_start,pdims%i_end

! adjust rho_dry in case required for q correction
              rho_dry(i,j,k) = rho_dry(i,j,k)*mass_corr
! use new rho dry to recalculate rho
              rho_r2(i,j,k) = rho_dry(i,j,k)*dry_to_wet(i,j,k)

            END DO
          END DO
        END DO

      END IF        ! lmass_corr

!----------------------------------------------------------------------
! Section 5. Moisture conservation
!----------------------------------------------------------------------
! total moisture at this time

        tot_m_final = (sum_results(ip_q)+sum_results(ip_qcl)            &
                     +sum_results(ip_qcf)) *factor

      IF (lqt_corr.or.lemq_print) THEN

! total moisture into atmosphere over energy correction period

        tot_moist_in=0.0
        moist_in(1)=0.0

        CALL global_2d_sums(sum_moist_flux,pdims%i_end,pdims%j_end,0,0, &
                            1, moist_in)

        tot_moist_in = moist_in(1)*factor

!----------------------------------------------------------------------
! Work out drift in moisture
!----------------------------------------------------------------------

        change_moist = tot_m_final - tot_m_init
        err_moist    = change_moist - tot_moist_in
        per_err      = 100.*err_moist/change_moist

        IF (mype == 0) THEN
          WRITE(6,fmt='(A,E15.5,A)') 'Final   moisture              = ',&
                      tot_m_final, ' KG'
          WRITE(6,fmt='(A,E15.5,A)') 'Initial moisture              = ',&
                      tot_m_init, ' KG'
          WRITE(6,fmt='(A,E15.5,A)') 'change in moisture            = ',&
                     change_moist, ' KG'
          WRITE(6,fmt='(A,E15.5,A)') 'Moisture added E-P in period  = ',&
                      tot_moist_in, ' KG'
          WRITE(6,fmt='(A,E15.5,A)') 'Error in moisture             = ',&
                      err_moist, ' KG'
          WRITE(6,fmt='(A,E15.5  )') 'Error as % of change          = ',&
                      per_err

          WRITE(6,fmt='(A,E15.5)') 'q ( rho cal)      ',                &
                      sum_results(ip_q)*factor
          WRITE(6,fmt='(A,E15.5)') 'qcf( rho cal)     ',                &
                      sum_results(ip_qcl)*factor
          WRITE(6,fmt='(A,E15.5)') 'qcl( rho cal)     ',                &
                      sum_results(ip_qcf)*factor

        END IF


!----------------------------------------------------------------------
! Apply moisture correction - highly experimental and not simple.
!----------------------------------------------------------------------
! Not tested as basic conservation of moist very poor when running
! using UM5.0
!----------------------------------------------------------------------
! Note problems as rho and q are related, therefore must alter
! both fields
! At present all error in qt applied to q (no correction of qcl or qcf)
!
!   rho = rho_dry / (1-q-qcl-qcf)
!
!  total global qt = sum qt*rho_dry*r^2 dr dlambda dphi cos(phi)
!                    -------------------------------------------
!                       (1-q-qcl-qcf)
!
!Let (qnew+qcl+qcf)/(1-qnew-qcl-qcf)= f (qold+qcl+qcf)/(1-qold-qcl-qcf)
! then
!     qnew =[ A/(1+A) -qcl -qcf ]
!  where
!     A= f ((qold+qcl+qcf)/(1-qold-qcl-qcf))
!
! where f = 1 - (error in qt)/(total qt old)
!
!  rho new = rho_dry / (1 - qnew-qcl-qcf)
!----------------------------------------------------------------------

        IF (lqt_corr) THEN   ! correct water vapour and wet_rho

! correction factor f
          moist_corr = 1. - err_moist/tot_m_init

! Firstly correct q (on q levels or should this be q_rho ?)

          DO k = pdims%k_start,pdims%k_end
            DO j = pdims%j_start,pdims%j_end
              DO i= pdims%i_start,pdims%i_end
                q_cld = qcl(i,j,k)+qcf(i,j,k)
!               tempw holds A
                tempw  = moist_corr * (q(i,j,k) + q_cld)/               &
                                     (1. - q(i,j,k) - q_cld)
                q(i,j,k) = ( tempw/(1.+tempw)  ) - q_cld

              END DO
            END DO
          END DO

! Work out for this new q the dry to wet conversion this order to the
! calculations should ensure dry mass is conserved.
! Using linear interpolation weights.

          IF ( l_endgame ) THEN
          k = 1
          DO j = pdims%j_start,pdims%j_end
            DO i = pdims%i_start,pdims%i_end

              dry_to_wet(i,j,k)= 1./(1. - q(i,j,k) - qcl(i,j,k)         &
                                          - qcf(i,j,k))

            END DO
          END DO

          DO k = 2, pdims%k_end
            DO j = pdims%j_start,pdims%j_end
              DO i = pdims%i_start,pdims%i_end

                weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
                weight2 = r_rho_levels(i,j,k)-r_theta_levels(i,j,k-1)
                weight3 = r_theta_levels(i,j,k)-r_theta_levels(i,j,k-1)

                tempd = ( weight2 *                                     &
                    (1. - q(i,j,k)- qcl(i,j,k) - qcf(i,j,k) ) +         &
                      weight1 *                                         &
                    (1. - q(i,j,k-1)- qcl(i,j,k-1) - qcf(i,j,k-1) ) )   &
                     / weight3

                dry_to_wet(i,j,k)= 1./tempd

              END DO
            END DO
          END DO


! convert rho dry back to rho wet
          DO k = pdims%k_start,pdims%k_end
            DO j = pdims%j_start,pdims%j_end
              DO i = pdims%i_start,pdims%i_end

                rho_r2(i,j,k) = rho_dry(i,j,k)*dry_to_wet(i,j,k)

              END DO
            END DO
          END DO

          ELSE

          k = 1
          DO j = pdims%j_start,pdims%j_end
            DO i = pdims%i_start,pdims%i_end

              dry_to_wet(i,j,k)= 1./(1. - q(i,j,k) - qcl(i,j,k)         &
                                       - qcf(i,j,k))

            END DO
          END DO

          DO k = 2, wet_model_levels
            DO j = pdims%j_start,pdims%j_end
              DO i = pdims%i_start,pdims%i_end
                 weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
                 weight2 = r_rho_levels(i,j,k)-r_theta_levels(i,j,k-1)
                 weight3 = r_theta_levels(i,j,k)                        &
                                            - r_theta_levels(i,j,k-1)
                 tempd = ( weight2 *                                    &
                    (1. - q(i,j,k)- qcl(i,j,k) - qcf(i,j,k) ) +         &
                      weight1 *                                         &
                    (1. - q(i,j,k-1)- qcl(i,j,k-1) - qcf(i,j,k-1) ) )   &
                     / weight3
                 dry_to_wet(i,j,k)= 1./tempd
              END DO
            END DO
          END DO

! Special case of wet_model_levels < model_levels

          IF (wet_model_levels <  pdims%k_end) THEN
            k=wet_model_levels+1
            DO j = pdims%j_start,pdims%j_end
              DO i = pdims%i_start,pdims%i_end
                weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
                weight2 = r_rho_levels(i,j,k)-r_theta_levels(i,j,k-1)
                weight3 = r_theta_levels(i,j,k)                         &
                                            - r_theta_levels(i,j,k-1)
                tempd = ( weight2 +                                     &
                          weight1 *                                     &
                         (1. - q  (i,j,k-1)- qcl(i,j,k-1)               &
                             - qcf(i,j,k-1)) ) / weight3

                dry_to_wet(i,j,k)= 1./tempd
              END DO
            END DO

          END IF   ! wet_model_levels < model_levels

! convert rho dry back to rho wet
          DO k = 1, wet_model_levels
            DO j = pdims%j_start,pdims%j_end
              DO i = pdims%i_start,pdims%i_end
                rho_r2(i,j,k) = rho_dry(i,j,k)*dry_to_wet(i,j,k)
              END DO
            END DO
          END DO

! Case of wet_model_levels < model_levels

          IF (wet_model_levels <  pdims%k_end ) THEN
            k=wet_model_levels+1
            DO j = pdims%j_start,pdims%j_end
              DO i = pdims%i_start,pdims%i_end
                rho_r2(i,j,k) = rho_dry(i,j,k)*dry_to_wet(i,j,k)
              END DO
            END DO
          END IF
          END IF

        END IF ! lqt_corr

      END IF ! lqt_corr or print
!----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('ENG_MASS_DIAG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ENG_MASS_DIAG
      END MODULE eng_mass_diag_mod
