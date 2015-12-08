! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eg_pofil_vatp
!
! Purpose:
!          Filter/diffusion based on first order
!          conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.
!          Applies z-level diffusion
!
! Method:
!          For any field. Pointers used to input grid
!          This version tested for v-at-the-poles/ENDGAME
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      MODULE  eg_pofil_zlevel_mod
      CONTAINS
      SUBROUTINE eg_pofil_zlevel(                                       &
                            field, fld_type, field_dims,                & 
                            field_dims_s, field_dims_l,                 &
                            field_dims_offj, field_dims_offj_s,         &
                            field_dims_offj_l, field_dims_offi,         &
                            field_dims_offi_s, field_dims_offi_l,       &
                            row_length, in_rows, levels, offx, offy,    &
                            delta_lambda, delta_phi, active_levels,     &
                            metric_level,                               &
                            r_levels, r_levels_offk,                    &
                            r_levels_offi, r_levels_offj,               &
                            eta_levels, eta_levels_offk,                &
                            lambda, lambda_offi, phi, phi_offj,         &
                            n_procy, max_filter_rows, global_filter,    &
                            sweeps, sweep_begin, sweep_end,             &
                            diff_coeff_phi, diff_coeff, L_diff,         &
                            L_vector, L_pofil_hadgem2,                  &
                            csphi_on, csphi_off,L_u_ns_diff)

      USE atm_fields_bounds_mod
      USE field_types
      USE ereport_mod, ONLY: ereport
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      TYPE (array_dims) :: field_dims, field_dims_s, field_dims_l,      &
                           field_dims_offi, field_dims_offi_s,          &
                           field_dims_offi_l,                           &
                           field_dims_offj, field_dims_offj_s,          &
                           field_dims_offj_l                            

      LOGICAL, Intent(In) ::                                            &
        L_diff,                                                         &
                         ! true IF diffusion active
        L_vector         ! true IF vector field

      LOGICAL, Intent(In) :: L_pofil_hadgem2       ! setting for hadgem2

      LOGICAL, Intent(In) :: L_u_ns_diff ! true if diffusing u

                             
      INTEGER, Intent(In) ::                                            &
        max_filter_rows,                                                &
                         ! max dimension for sweeping arrays
        global_filter,                                                  &
                         ! max number of filter sweeps 0=diffusion only 
        fld_type,                                                       &
                      ! field type (p=1, u=2 or v=3)
        n_procy,                                                        &
                       ! Number of processors in latitude
        active_levels,                                                  &
                          ! number of levels to be filtered
        metric_level                                                    
                        ! uppermost non-constant layer for metric terms

      REAL, Intent(In) :: delta_lambda,  delta_phi

! field dimensions for swap bounds
       INTEGER, Intent(In) ::     row_length, in_rows, levels, offx, offy


! Number of sweeps and NS sweeping limits
      INTEGER, Intent(In) ::                                            &
        sweeps   (  max_filter_rows),                                   &
        sweep_begin(0:max_filter_rows),                                 &
        sweep_end  (0:max_filter_rows)                                      

      REAL, Intent(In) ::                                               &
! EW diffusion coefficient
                 diff_coeff(field_dims_s%i_start:field_dims_s%i_end,    &
                            field_dims_s%j_start:field_dims_s%j_end),   &
! NS diffusion coefficient
             diff_coeff_phi(field_dims_s%i_start:field_dims_s%i_end,    &
                            field_dims_s%j_start:field_dims_s%j_end)

! r level arrays
      REAL, Intent(In) ::                                               &
        r_levels (field_dims_l%i_start:field_dims_l%i_end,              &
                  field_dims_l%j_start:field_dims_l%j_end,              &
                  field_dims_l%k_start:field_dims_l%k_end),             &
        r_levels_offk(field_dims_l%i_start:field_dims_l%i_end,          &
                      field_dims_l%j_start:field_dims_l%j_end,          &
                    1-field_dims_l%k_start:field_dims_l%k_end),         &
        r_levels_offi(field_dims_offi_l%i_start:field_dims_offi_l%i_end,&
                      field_dims_l%j_start     :field_dims_l%j_end,     &
                      field_dims_l%k_start:field_dims_l%k_end),         &
        r_levels_offj(field_dims_l%i_start     :field_dims_l%i_end,     &
                      field_dims_offj_l%j_start:field_dims_offj_l%j_end,&
                      field_dims_l%k_start:field_dims_l%k_end)
!         r_levels_offi(field_dims_offi_l%i_start:field_dims_offi_l%i_end,&
!                       field_dims_offi_l%j_start:field_dims_offi_l%j_end,&
!                       field_dims_l%k_start:field_dims_l%k_end),         &
!         r_levels_offj(field_dims_offj_l%i_start:field_dims_offj_l%i_end,&
!                       field_dims_offj_l%j_start:field_dims_offj_l%j_end,&
!                       field_dims_l%k_start:field_dims_l%k_end)

! eta level arrays
      REAL, Intent(In) ::                                               &
        eta_levels     (  field_dims_l%k_start:field_dims_l%k_end),     &
        eta_levels_offk(1-field_dims_l%k_start:field_dims_l%k_end)

! Horizontal grid arrays
      REAL, Intent(In) ::                                               &
        lambda     (field_dims_l%i_start:field_dims_l%i_end),           &
        phi        (field_dims_l%j_start:field_dims_l%j_end),           &
        lambda_offi(field_dims_offi_l%i_start:field_dims_offi_l%i_end), &
        phi_offj   (field_dims_offj_l%j_start:field_dims_offj_l%j_end)  


!       REAL, Intent(In) ::                                              &
!      &  r_levels (1-halo_ri:row_length+halo_ri,                        &
!      &            1-halo_rj:in_rows+halo_rj, base_f:model_levels)      &
!      &, r_half_levels (1-halo_hi:row_length+halo_hi,                   &
!      &                 1-halo_hj:in_rows+halo_hj, base_h:model_levels) &
!      &, r_midi (1-halo_ii:row_length+halo_ii,                          &
!      &          1-halo_ij:in_rows+halo_ij, base_i:model_levels)        &
!      &, r_midj (1-halo_ji:row_length+halo_ji,                          &
!      &          1-halo_jj:cos_rows+halo_jj, base_j:model_levels)


     REAL, Intent(In)    ::                                             &
           csphi_on(field_dims_l%j_start:field_dims_l%j_end),           &
           csphi_off(field_dims_offj_l%j_start:field_dims_offj_l%j_end)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      REAL, Intent(InOut) ::                                            &
           field (field_dims_s%i_start:field_dims_s%i_end,              &
                  field_dims_s%j_start:field_dims_s%j_end,              &
                  field_dims_s%k_start:field_dims_s%k_end)


! Local Variables.
      INTEGER :: i, j, k, i_filter, i_sweep, num_pass,                  &
                 j_sweep_begin, j_sweep_end, j_sweep_store

      INTEGER :: base_k, base_i, base_j

      LOGICAL :: L_cycle, L_combine

! Local copies of NS sweeping limits
      INTEGER             ::                                            &
        idx_begin(0:max_filter_rows),                                   &
        idx_end  (0:max_filter_rows)  

! Local copies of cos(phi) for NS U diffusion
     REAL  ::                                                           &
           csphi_on_u(field_dims_l%j_start:field_dims_l%j_end),         &
           csphi_off_u(field_dims_offj_l%j_start:                       &
                       field_dims_offj_l%j_end),                        &
           u_fac(field_dims_l%j_start:field_dims_l%j_end)

      REAL ::                                                           &
        delta_z(field_dims_s%i_start:field_dims_s%i_end,                &
                field_dims_s%j_start:field_dims_s%j_end,                &
                field_dims_s%k_start:metric_level ),                    &
        recip_r_squared_delz(field_dims_s%i_start:field_dims_s%i_end,   &
                             field_dims_s%j_start:field_dims_s%j_end,   &
                             field_dims_s%k_start:metric_level)

      REAL ::                                                           &
        wt_z_off2on(field_dims_s%i_start:field_dims_s%i_end,            &
                    field_dims_s%j_start:field_dims_s%j_end,            &
                    1:field_dims_s%k_end),                              &
        wt_z_on2off(field_dims_s%i_start:field_dims_s%i_end,            &
                    field_dims_s%j_start:field_dims_s%j_end,            &
                    1:field_dims_s%k_end),                              &
        wt_x_off2on(field_dims_s%i_start:field_dims_s%i_end,            &
                    field_dims_s%j_start:field_dims_s%j_end),           &
        wt_y_off2on(field_dims_s%i_start:field_dims_s%i_end,            &
                    field_dims_s%j_start:field_dims_s%j_end),           &
!         wt_x_on2off(field_dims_offi_s%i_start:field_dims_offi_s%i_end,  &
!                     field_dims_offi_s%j_start:field_dims_offi_s%j_end), &
!         wt_y_on2off(field_dims_offj_s%i_start:field_dims_offj_s%i_end,  &
!                     field_dims_offj_s%j_start:field_dims_offj_s%j_end) 
        wt_x_on2off(field_dims_offi_s%i_start:field_dims_offi_s%i_end,  &
                    field_dims_s%j_start     :field_dims_s%j_end),      &
        wt_y_on2off(field_dims_s%i_start     :field_dims_s%i_end,       &
                    field_dims_offj_s%j_start:field_dims_offj_s%j_end) 

! Arrays used for trisolve
      REAL ::    a(field_dims%i_start:field_dims%i_end,                 &
                  field_dims%j_start:field_dims%j_end,                  &
                  field_dims%k_start:field_dims%k_end),                 &
                ap(field_dims%i_start:field_dims%i_end,                 &
                  field_dims%j_start:field_dims%j_end,                  &
                  field_dims%k_start:field_dims%k_end),                 &
                am(field_dims%i_start:field_dims%i_end,                 &
                  field_dims%j_start:field_dims%j_end,                  &
                  field_dims%k_start:field_dims%k_end),                 &
                 b(field_dims%i_start:field_dims%i_end,                 &
                  field_dims%j_start:field_dims%j_end,                  &
                  field_dims%k_start:field_dims%k_end) 

      REAL :: G_lambda(field_dims_offi_s%i_start:field_dims_offi_s%i_end, &
                       field_dims_s%j_start     :field_dims_s%j_end,      &
                       field_dims_s%k_start     :field_dims_s%k_end),     &
              G_lambda_eta(field_dims_s%i_start :field_dims_s%i_end,      &
                           field_dims_s%j_start :field_dims_s%j_end,      &
                           0                    :metric_level+1),         &
              G_phi   (field_dims_s%i_start     :field_dims_s%i_end,      &
                       field_dims_offj_s%j_start:field_dims_offj_s%j_end, &
                       field_dims_s%k_start     :field_dims_s%k_end),     &
              G_phi_eta(field_dims_s%i_start    :field_dims_s%i_end,      &
                        field_dims_s%j_start    :field_dims_s%j_end,      &
                        0                       :metric_level+1)

      REAL :: f_av_xz_k, f_av_xz_km, f_av_xz_i, f_av_xz_im,               &
              f_av_yz_k, f_av_yz_km, f_av_yz_j, f_av_yz_jm, del_r_av_z

      REAL :: flux_EW

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! External functions
INTERFACE
  ! For tri-diagonal solver
  FUNCTION eg_vert_trisolve(a, b, c, d, KK, tri_idx) 
    IMPLICIT NONE
INTEGER, INTENT(IN)     :: KK, tri_idx
REAL                    :: eg_vert_trisolve(tri_idx:KK) 
REAL, INTENT(IN)        :: b(tri_idx:KK), c(tri_idx:KK-1), d(tri_idx:KK)
REAL, INTENT(IN)        :: a(tri_idx+1:KK) 
  END FUNCTION eg_vert_trisolve
END INTERFACE

! ----------------------------------------------------------------------
! Section 1.   Set metric terms and calculate D(r)/D(eta)
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('EG_POFIL_ZLEVEL',zhook_in,zhook_handle)


      base_k = 1 - field_dims%k_start
      base_i = 1 - field_dims%i_start
      base_j = 1 - field_dims%j_start

! WRITE(6,*) 'At top of pofil_zlevel'
! WRITE(6,*) base_i,base_j,base_k
! i = 10
! WRITE(6,*) ' If i = ',i
! WRITE(6,*) ' then lambda(i) = ',lambda(i)
! WRITE(6,*) ' this lies between lambda_offi(i-1+base_i,i+base_i) = ', &
!              lambda_offi(i-1+base_i),lambda_offi(i+base_i)
! j = 25
! WRITE(6,*) ' If j = ',j
! WRITE(6,*) ' then phi(j) = ',phi(j)
! WRITE(6,*) ' this lies between phi_offj(j-1+base_j,j+base_j) = ', &
!              phi_offj(j-1+base_j),phi_offj(j+base_j)
! flush(6)

!   No metric terms needed above first constant rho-level
!   since they cancel on constant r-surfaces
      DO k = 1,metric_level
        DO j = field_dims_s%j_start,field_dims_s%j_end
          DO i = field_dims_s%i_start,field_dims_s%i_end
            delta_z(i,j,k) = r_levels_offk(i,j,k+base_k) -               &
                             r_levels_offk(i,j,k+base_k-1)
            recip_r_squared_delz(i,j,k) = 1.0 / ( r_levels(i,j,k) *      &
                                                  r_levels(i,j,k) *      &
                                                   delta_z(i,j,k) )            
          END DO
        END DO
      END DO

      IF ( field_dims%k_start == 0 ) THEN
        DO j = field_dims_s%j_start,field_dims_s%j_end
          DO i = field_dims_s%i_start,field_dims_s%i_end
            delta_z(i,j,0) = r_levels_offk(i,j,1) -                      &
                             r_levels(i,j,0)
            recip_r_squared_delz(i,j,0) = 1.0 / ( r_levels(i,j,0) *      &
                                                  r_levels(i,j,0) *      &
                                                   delta_z(i,j,0) )
          END DO
        END DO
      END IF

! ----------------------------------------------------------------------
! Section 1.1   Set up averaging arrays 
! ----------------------------------------------------------------------
! Vertical Weights
      DO k = field_dims_s%k_start+base_k,field_dims_s%k_end-field_dims_s%k_start
        DO j = field_dims_s%j_start,field_dims_s%j_end
          DO i = field_dims_s%i_start,field_dims_s%i_end

! Endgame does averaging in eta
            ! Vertical weight from  k+/- 1/2 to  field level         
!             wt_z_off2on(i,j,k) = (eta_levels(k) -                        & 
!                                         eta_levels_offk(k-1+base_k))     &
!                                 /(eta_levels_offk(k+base_k) -            &
!                                         eta_levels_offk(k-1+base_k))
            ! Vertical weight from field level to k+/- 1/2      
            wt_z_on2off(i,j,k) = (eta_levels_offk(k) -                   & 
                                        eta_levels(k-base_k))            &
                                /(eta_levels(k+1-base_k) -               &
                                        eta_levels(k-base_k))  
          END DO
        END DO
      END DO


! horizontal weights
      DO j = field_dims_s%j_start,field_dims_s%j_end
        DO i = field_dims_s%i_start,field_dims_s%i_end
          ! Lambda weight from  i+/- 1/2 to field point 
          wt_x_off2on(i,j) = (lambda(i) - lambda_offi(i-1+base_i))        &
                            /(lambda_offi(i+base_i) -                     &
                                        lambda_offi(i-1+base_i))
          ! phi weight from  j+/- 1/2 to field point 
          wt_y_off2on(i,j) = (phi(j)             - phi_offj(j-1+base_j))  &
                            /(phi_offj(j+base_j) - phi_offj(j-1+base_j))

        END DO
      END DO

!       DO j = field_dims_offi_s%j_start,field_dims_offi_s%j_end
      DO j = field_dims_s%j_start,field_dims_s%j_end
        DO i = field_dims_offi_s%i_start,field_dims_offi_s%i_end
          ! Lambda weight from  field point to i+/- 1/2 
          wt_x_on2off(i,j) = (lambda_offi(i)     - lambda(i-base_i))      &
                            /(lambda(i+1-base_i) - lambda(i-base_i))
        END DO
      END DO 

      DO j = field_dims_offj_s%j_start,field_dims_offj_s%j_end
!         DO i = field_dims_offj_s%i_start,field_dims_offj_s%i_end
        DO i = field_dims_s%i_start,field_dims_s%i_end
          ! Phi weight from  field point to j+/- 1/2 
          wt_y_on2off(i,j) = (phi_offj(j)     - phi(j-base_j))            &
                            /(phi(j+1-base_j) - phi(j-base_j)) 
        END DO
      END DO


      wt_x_on2off = 0.5
      wt_y_on2off = 0.5
!       wt_z_on2off = 0.5

! ----------------------------------------------------------------------
! Section 2.0  Filtering and diffusion
! ----------------------------------------------------------------------
      DO i= 0, global_filter
        idx_begin(i) = sweep_begin(i)
        idx_end(i) = sweep_end(i)
      END DO
! Shift in sweeping index for Endgame v points
      IF ( fld_type == fld_type_v ) THEN
        DO i= 0, global_filter
          idx_begin(i) = sweep_begin(i)-1
          idx_end(i) = sweep_end(i)-1
        END DO
      END IF

      L_combine = .false.
      num_pass = global_filter
      j_sweep_begin = idx_begin(1)
      j_sweep_end = idx_end(1)
! Combine EW diffusion with first sweep of polar filter
      IF ( L_diff ) THEN
        L_combine = .true.
        IF( global_filter < 1 ) num_pass = 1
        IF ( j_sweep_begin < 0 ) THEN
          j_sweep_begin = idx_begin(0)
          j_sweep_end = idx_end(0)
        ELSEif ( j_sweep_end < idx_end(0) ) THEN
          j_sweep_end = idx_end(0)
        ELSEif ( j_sweep_begin > idx_begin(0) ) THEN
          j_sweep_begin = idx_begin(0)
        END IF !  j_sweep_begin < 0
      END IF ! L_diff

!       WRITE(6,*) 'Altered J sweeps'
!       DO i= 0, num_pass
!         WRITE(6,*) i,idx_begin(i),idx_end(i)
!       END DO

      DO i_filter = 1, num_pass
!         WRITE(6,*) 'Starting pass :',i_filter
!         flush(6)
        IF( i_filter > 1 ) THEN
          L_combine = .false.
          j_sweep_begin = idx_begin(i_filter)
          j_sweep_end = idx_end(i_filter)
        END IF ! i_filter > 1

       L_cycle = .false.
!  IF global and n_procy = 1, need to filter hemispheres separately
        IF ( n_procy == 1) L_cycle = .TRUE.

        i_sweep = 1
        DO  ! Sweeping loop  i_sweep = 1, sweeps(i_filter)

        G_lambda = 0.0
        G_lambda_eta = 0.0

!           WRITE(6,*) 'EW J Sweeps ',i_filter,i_sweep,j_sweep_begin,j_sweep_end
! ----------------------------------------------------------------------
! Section 4.1   EW filtering and diffusion
! ----------------------------------------------------------------------
          SELECT CASE ( fld_type ) 
! -----------------------------------------------------------------------------
! --- U fields ---
          CASE ( fld_type_u )

            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims_offi_s%i_start,field_dims_offi%i_end

! This is G_lambda stored on p points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                  IF ( k == field_dims%k_end ) THEN
                    f_av_xz_k =       wt_x_on2off(i,j) *field(i-1,j,k) +    &
                                 (1.0-wt_x_on2off(i,j))*field(i,j,k) 
                  ELSE
                    f_av_xz_k = wt_z_on2off(i,j,k)*(                        &
                                      wt_x_on2off(i,j) *field(i-1,j,k+1)  + &
                                 (1.0-wt_x_on2off(i,j))*field(i,j,k+1)) +   &
                           (1.0-wt_z_on2off(i,j,k))*(                       &
                                      wt_x_on2off(i,j) *field(i-1,j,k) +    &
                                 (1.0-wt_x_on2off(i,j))*field(i,j,k) )
                  END IF
                  IF ( k == 1 ) THEN
                    f_av_xz_km=  wt_x_on2off(i,j) *field(i-1,j,k) +         &
                            (1.0-wt_x_on2off(i,j))*field(i,j,k) 

                  ELSE
                    f_av_xz_km= wt_z_on2off(i,j,k-1)*(                      &
                                      wt_x_on2off(i,j) *field(i-1,j,k) +    &
                                 (1.0-wt_x_on2off(i,j))*field(i,j,k) ) +    &
                           (1.0-wt_z_on2off(i,j,k-1))*(                     &
                                      wt_x_on2off(i,j)* field(i-1,j,k-1) +  &
                                 (1.0-wt_x_on2off(i,j))*field(i,j,k-1) )
                  END IF

                  G_lambda(i,j,k) = (( field(i,j,k) - field(i-1,j,k))         &
                                     *delta_z(i,j,k)                          &
                                   - (f_av_xz_k - f_av_xz_km   )              &
                                    *(r_levels(i,j,k) - r_levels(i-1,j,k)) )  &
                                   /( lambda(i) - lambda(i-1))                &
                                   *delta_lambda*delta_lambda*diff_coeff(i,j) &
                                   *r_levels_offi(i,j,k) * r_levels_offi(i,j,k)
                END DO
              END DO
            END DO

            DO k = 1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_lambda_eta on (i,j-1/2,k) points                    
                  f_av_xz_i = wt_z_on2off(i,j,k)*(                           &
                                wt_x_on2off(i,j)*field(i,j,k+1) +            &
                           (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1)) +        &
                         (1.0-wt_z_on2off(i,j,k))*(                          &
                                wt_x_on2off(i,j)*field(i,j,k) +              &
                           (1.0-wt_x_on2off(i,j))*field(i+1,j,k) ) 
                 f_av_xz_im= wt_z_on2off(i,j,k)*(                            &
                                wt_x_on2off(i,j)*field(i-1,j,k+1)  +         &
                           (1.0-wt_x_on2off(i,j))*field(i,j,k+1) ) +         &
                         (1.0-wt_z_on2off(i,j,k))*(                          &
                                wt_x_on2off(i,j)*field(i-1,j,k) +            &
                           (1.0-wt_x_on2off(i,j))*field(i,j,k) )

                  G_lambda_eta(i,j,k) = (( f_av_xz_i - f_av_xz_im ))         &
                                   *0.5*((r_levels_offi(i+1,j,k+1)           &
                                                    -r_levels_offi(i,j,k+1)) &
                                       + (r_levels_offi(i+1,j,k)             &
                                                     -r_levels_offi(i,j,k))) &
                              /( lambda_offi(i+1) - lambda_offi(i))**2       &
                              *delta_lambda*delta_lambda*diff_coeff(i,j)     &
                              *r_levels_offk(i,j,k) * r_levels_offk(i,j,k)
                END DO
              END DO               
            END DO
            ! Top & Bottom boundary condition
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                G_lambda_eta(i,j,metric_level) = 0.0
                G_lambda_eta(i,j,0) = 0.0
              END DO
            END DO

            k = 1
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                del_r_av_z = 0.5*((r_levels_offi(i+1,j,k+1)           &
                                             -r_levels_offi(i,j,k+1)) &
                                + (r_levels_offi(i+1,j,k)             &
                                              -r_levels_offi(i,j,k)))

                ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &        
                            /((lambda_offi(i+1) - lambda_offi(i))**2  &
                            *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                            *delta_lambda**2*diff_coeff(i,j)          &
                            *r_levels_offk(i,j,k)**2
                am(i,j,k) = 0.0
              END DO
            END DO
            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  ! LHS terms

                  del_r_av_z = 0.5*((r_levels_offi(i+1,j,k+1)           &
                                               -r_levels_offi(i,j,k+1)) &
                                  + (r_levels_offi(i+1,j,k)             &
                                                -r_levels_offi(i,j,k)))

                  ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &        
                              /((lambda_offi(i+1) - lambda_offi(i))**2  &
                              *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                              *delta_lambda**2*diff_coeff(i,j)          &
                              *r_levels_offk(i,j,k)**2

                  del_r_av_z = 0.5*((r_levels_offi(i+1,j,k)             &
                                               -r_levels_offi(i,j,k))   &
                                  + (r_levels_offi(i+1,j,k-1)           &
                                               -r_levels_offi(i,j,k-1)))

                  am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &        
                              /((lambda_offi(i+1) - lambda_offi(i))**2  &
                              *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                              *delta_lambda**2*diff_coeff(i,j)          &
                              *r_levels_offk(i,j,k-1)**2
                END DO
              END DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                ap(i,j,k) = 0.0

                del_r_av_z = 0.5*((r_levels_offi(i+1,j,k)             &
                                             -r_levels_offi(i,j,k))   &
                                + (r_levels_offi(i+1,j,k-1)           &
                                             -r_levels_offi(i,j,k-1)))

                am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &        
                            /((lambda_offi(i+1) - lambda_offi(i))**2  &
                            *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                            *delta_lambda**2*diff_coeff(i,j)          &
                            *r_levels_offk(i,j,k-1)**2
              END DO
            END DO
            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  a(i,j,k) = 1.0 - ap(i,j,k) - am(i,j,k)
                    ! RHS term
                  b(i,j,k) = field(i,j,k) +                                &
                            ( (G_lambda(i+1,j,k) - G_lambda(i,j,k))        &
                             /(lambda_offi(i+1)  - lambda_offi(i))         &
                            - (G_lambda_eta(i,j,k)                         &
                                                 - G_lambda_eta(i,j,k-1))) &
                             * recip_r_squared_delz(i,j,k) 
                END DO
              END DO
            END DO

            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end 

!DEPENDS ON: eg_vert_trisolve
                field(i,j,1:metric_level) =                                  &
                       eg_vert_trisolve(                                     &
                                        am(i,j,2:metric_level),              &
                                        a(i,j,1:metric_level),               &
                                        ap(i,j,1:metric_level-1),            &
                                        b(i,j,1:metric_level),               &
                                        metric_level,1 ) 

              END DO
            END DO

          DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
            DO j = j_sweep_begin, j_sweep_end

!               DO i = field_dims%i_start-1,field_dims%i_end 
              DO i = field_dims_offi_s%i_start,field_dims_offi%i_end
                G_lambda(i,j,k) = ( field(i,j,k) - field(i-1,j,k) )          &
                                 /( lambda(i)    - lambda(i-1)    )          & 
                           * delta_lambda**2*diff_coeff(i,j)                     
              END DO

              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = field(i,j,k)                                  &
                            + (G_lambda(i+1,j,k) - G_lambda(i,j,k))          &
                              /(lambda_offi(i+1)  - lambda_offi(i))         

              END DO
            END DO
          END DO

! -----------------------------------------------------------------------------
! --- V fields ---
          CASE ( fld_type_v ) 

            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims_offi%i_start,field_dims_offi_s%i_end

! This is G_lambda stored on points (i,j,k-1/2)
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  IF ( k == field_dims%k_end ) THEN
                    f_av_xz_k = wt_x_on2off(i,j) *field(i,j,k) +            &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) 
                  ELSE
                    f_av_xz_k = wt_z_on2off(i,j,k)*(                        &
                                      wt_x_on2off(i,j) *field(i,j,k+1)    + &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1)) + &
                           (1.0-wt_z_on2off(i,j,k))*(                       &
                                      wt_x_on2off(i,j) *field(i,j,k) +      &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )
                  END IF

                  IF ( k == 1 ) THEN
                    f_av_xz_km= wt_x_on2off(i,j) *field(i,j,k) +             &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k)
                  ELSE
                    f_av_xz_km= wt_z_on2off(i,j,k-1)*(                      &
                                      wt_x_on2off(i,j) *field(i,j,k) +      &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) ) +  &
                           (1.0-wt_z_on2off(i,j,k-1))*(                     &
                                      wt_x_on2off(i,j)* field(i,j,k-1) +    &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k-1) )
                  END IF

                  G_lambda(i,j,k) = (( field(i+1,j,k) - field(i,j,k)    )     &
                                   *delta_z(i,j,k)                            &
                                   - ( f_av_xz_k - f_av_xz_km   )             &
                                   *(r_levels(i+1,j,k) - r_levels(i,j,k)))    & 
                                   /( lambda(i+1) - lambda(i))                &
                                   *delta_lambda*delta_lambda*diff_coeff(i,j) &
                                   *r_levels_offi(i,j,k) * r_levels_offi(i,j,k)

                END DO
              END DO
            END DO

            DO k = 1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_lambda_eta on (i-1/2,j,k) points
                  f_av_xz_i = wt_z_on2off(i,j,k)*(                           &
                                wt_x_on2off(i,j)*field(i,j,k+1-base_k)  +    &
                           (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k)) + &
                         (1.0-wt_z_on2off(i,j,k))*(                          &
                                wt_x_on2off(i,j)*field(i,j,k-base_k) +       &
                           (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )
                  f_av_xz_im= wt_z_on2off(i,j,k)*(                           &
                                wt_x_on2off(i,j)*field(i-1,j,k+1-base_k)  +  &
                           (1.0-wt_x_on2off(i,j))*field(i,j,k+1-base_k) ) +  &
                         (1.0-wt_z_on2off(i,j,k))*(                          &
                                wt_x_on2off(i,j)*field(i-1,j,k-base_k) +     &
                           (1.0-wt_x_on2off(i,j))*field(i,j,k-base_k) )

                  G_lambda_eta(i,j,k) = (( f_av_xz_i - f_av_xz_im ))         &
                                   *0.5*((r_levels_offi(i,j,k+1)             &
                                                  -r_levels_offi(i-1,j,k+1)) &
                                       + (r_levels_offi(i,j,k)               &
                                                   -r_levels_offi(i-1,j,k))) &
                              /( lambda_offi(i) - lambda_offi(i-1))**2       &
                              *delta_lambda*delta_lambda*diff_coeff(i,j)     &
                              *r_levels_offk(i,j,k) * r_levels_offk(i,j,k)

                END DO
              END DO
            END DO
            ! Top & Bottom boundary condition
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                G_lambda_eta(i,j,metric_level) = 0.0
                G_lambda_eta(i,j,0) = 0.0
              END DO
            END DO

            k = 1
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
               del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)             &
                                          -r_levels_offi(i-1,j,k+1)) &
                               + (r_levels_offi(i,j,k)               &
                                           -r_levels_offi(i-1,j,k))) 

                ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((lambda_offi(i) - lambda_offi(i-1))**2  &
                            *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                            *delta_lambda**2*diff_coeff(i,j)          &
                            *r_levels_offk(i,j,k)**2

                am(i,j,k) = 0.0
              END DO
            END DO
            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                 ! LHS terms

                 del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)             &
                                            -r_levels_offi(i-1,j,k+1)) &
                                 + (r_levels_offi(i,j,k)               &
                                             -r_levels_offi(i-1,j,k))) 

                  ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((lambda_offi(i) - lambda_offi(i-1))**2  &
                              *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                              *delta_lambda**2*diff_coeff(i,j)          &
                              *r_levels_offk(i,j,k)**2

                 del_r_av_z = 0.5*((r_levels_offi(i,j,k)               &
                                            -r_levels_offi(i-1,j,k))   &
                                 + (r_levels_offi(i,j,k-1)             &
                                             -r_levels_offi(i-1,j,k-1))) 

                  am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((lambda_offi(i) - lambda_offi(i-1))**2  &
                              *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                              *delta_lambda**2*diff_coeff(i,j)          &
                              *r_levels_offk(i,j,k-1)**2
                END DO
              END DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                ap(i,j,k) = 0.0

               del_r_av_z = 0.5*((r_levels_offi(i,j,k)               &
                                          -r_levels_offi(i-1,j,k))   &
                               + (r_levels_offi(i,j,k-1)             &
                                           -r_levels_offi(i-1,j,k-1))) 

                am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((lambda_offi(i) - lambda_offi(i-1))**2  &
                            *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                            *delta_lambda**2*diff_coeff(i,j)          &
                            *r_levels_offk(i,j,k-1)**2
              END DO
            END DO
            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                 ! LHS terms
                  a(i,j,k) = 1.0 - ap(i,j,k) - am(i,j,k)
                 ! RHS term
                  b(i,j,k) = field(i,j,k) +                                &
                            ( (G_lambda(i,j,k)   - G_lambda(i-1,j,k))      &
                             /(lambda_offi(i)    - lambda_offi(i-1))       &
                            - (G_lambda_eta(i,j,k)                         &
                                                 - G_lambda_eta(i,j,k-1))) &
                               * recip_r_squared_delz(i,j,k) 
                END DO
              END DO
            END DO

            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end 
!DEPENDS ON: eg_vert_trisolve
                field(i,j,1:metric_level) =                                  &
                       eg_vert_trisolve(                                     &
                                        am(i,j,2:metric_level),              &
                                        a(i,j,1:metric_level),               &
                                        ap(i,j,1:metric_level-1),            &
                                        b(i,j,1:metric_level),               &
                                        metric_level,1 ) 
              END DO
            END DO

          DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
            DO j = j_sweep_begin, j_sweep_end

              DO i = field_dims_offi%i_start,field_dims_offi_s%i_end 
                G_lambda(i,j,k) = ( field(i+1,j,k) - field(i,j,k) )          &
                                 /( lambda(i+1)    - lambda(i)    )          & 
                                 * delta_lambda**2*diff_coeff(i,j)            
              END DO

              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = field(i,j,k)                                   &
                            + (G_lambda(i,j,k) - G_lambda(i-1,j,k))           &
                             /(lambda_offi(i)  - lambda_offi(i-1))
              END DO
            END DO
          END DO

! -----------------------------------------------------------------------------
! --- P fields ---
          CASE ( fld_type_p )   

            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims_offi%i_start,field_dims_offi_s%i_end

! This is G_lambda stored on u points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  IF ( k == field_dims%k_end ) THEN
                    f_av_xz_k = wt_x_on2off(i,j) *field(i,j,k) +            &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) 
                  ELSE
                    f_av_xz_k = wt_z_on2off(i,j,k)*(                        &
                                      wt_x_on2off(i,j) *field(i,j,k+1)    + &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1)) + &
                           (1.0-wt_z_on2off(i,j,k))*(                       &
                                      wt_x_on2off(i,j) *field(i,j,k) +      &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )
                  END IF
                  IF ( k == 1 ) THEN
                    f_av_xz_km= wt_x_on2off(i,j) *field(i,j,k) +            &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k)
                  ELSE
                    f_av_xz_km= wt_z_on2off(i,j,k-1)*(                      &
                                      wt_x_on2off(i,j) *field(i,j,k) +      &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) ) +  &
                           (1.0-wt_z_on2off(i,j,k-1))*(                     &
                                      wt_x_on2off(i,j)* field(i,j,k-1) +    &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k-1) )
                  END IF
                    G_lambda(i,j,k) = (( field(i+1,j,k) - field(i,j,k)    )     &
                                     *delta_z(i,j,k)                            &
                                     - ( f_av_xz_k - f_av_xz_km   )             &
                                     *(r_levels(i+1,j,k) - r_levels(i,j,k)))    &
                                     /( lambda(i+1) - lambda(i))                &
                                     *delta_lambda*delta_lambda*diff_coeff(i,j) &
                                     *r_levels_offi(i,j,k) * r_levels_offi(i,j,k)
                END DO
              END DO
            END DO

            DO k = 1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_lambda_eta on w points
                  f_av_xz_i = wt_z_on2off(i,j,k)*(                           &
                                wt_x_on2off(i,j)*field(i,j,k+1-base_k)  +    &
                           (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k)) + &
                         (1.0-wt_z_on2off(i,j,k))*(                          &
                                wt_x_on2off(i,j)*field(i,j,k-base_k) +       &
                           (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )
                  f_av_xz_im= wt_z_on2off(i,j,k)*(                           &
                                wt_x_on2off(i,j)*field(i-1,j,k+1-base_k)  +  &
                           (1.0-wt_x_on2off(i,j))*field(i,j,k+1-base_k) ) +  &
                         (1.0-wt_z_on2off(i,j,k))*(                          &
                                wt_x_on2off(i,j)*field(i-1,j,k-base_k) +     &
                           (1.0-wt_x_on2off(i,j))*field(i,j,k-base_k) )


                  G_lambda_eta(i,j,k) = (( f_av_xz_i - f_av_xz_im ))         &
                                   *0.5*((r_levels_offi(i,j,k+1)             &
                                                  -r_levels_offi(i-1,j,k+1)) &
                                       + (r_levels_offi(i,j,k)               &
                                                   -r_levels_offi(i-1,j,k))) &
                              /( lambda_offi(i) - lambda_offi(i-1))**2       &
                              *delta_lambda*delta_lambda*diff_coeff(i,j)     &
                              *r_levels_offk(i,j,k) * r_levels_offk(i,j,k)
                END DO
              END DO
            END DO
            !Top & Bottom boundary condition
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                G_lambda_eta(i,j,metric_level) = 0.0
                G_lambda_eta(i,j,0) = 0.0
              END DO
            END DO

            k = 1
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms

                del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)             &
                                          -r_levels_offi(i-1,j,k+1))  &
                                + (r_levels_offi(i,j,k)               &
                                          -r_levels_offi(i-1,j,k))) 

                ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((lambda_offi(i) - lambda_offi(i-1))**2  &
                            *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                            *delta_lambda**2*diff_coeff(i,j)          &
                            *r_levels_offk(i,j,k)**2


                am(i,j,k) = 0.0
              END DO
            END DO
            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  ! LHS terms

                  del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)             &
                                            -r_levels_offi(i-1,j,k+1))  &
                                  + (r_levels_offi(i,j,k)               &
                                            -r_levels_offi(i-1,j,k))) 

                  ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((lambda_offi(i) - lambda_offi(i-1))**2  &
                              *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                              *delta_lambda**2*diff_coeff(i,j)          &
                              *r_levels_offk(i,j,k)**2

                  del_r_av_z = 0.5*((r_levels_offi(i,j,k)               &
                                            -r_levels_offi(i-1,j,k))    &
                                  + (r_levels_offi(i,j,k-1)             &
                                            -r_levels_offi(i-1,j,k-1))) 

                  am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((lambda_offi(i) - lambda_offi(i-1))**2  &
                              *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                              *delta_lambda**2*diff_coeff(i,j)          &
                              *r_levels_offk(i,j,k-1)**2
                END DO
              END DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                ap(i,j,k) = 0.0 

                del_r_av_z = 0.5*((r_levels_offi(i,j,k)               &
                                          -r_levels_offi(i-1,j,k))    &
                                + (r_levels_offi(i,j,k-1)             &
                                          -r_levels_offi(i-1,j,k-1))) 

                am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((lambda_offi(i) - lambda_offi(i-1))**2  &
                            *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                            *delta_lambda**2*diff_coeff(i,j)          &
                            *r_levels_offk(i,j,k-1)**2
              END DO
            END DO
            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                ! LHS terms
                  a(i,j,k) = 1.0 - ap(i,j,k) - am(i,j,k)
                  ! RHS term
                  b(i,j,k) = field(i,j,k) +                                &
                            ( (G_lambda(i,j,k)   - G_lambda(i-1,j,k))      &
                            /( lambda_offi(i)    - lambda_offi(i-1) )      & 
                            - (G_lambda_eta(i,j,k)                         &
                                                 - G_lambda_eta(i,j,k-1))) &
                            * recip_r_squared_delz(i,j,k) 
                END DO
              END DO
            END DO

            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end 
!DEPENDS ON: eg_vert_trisolve
                field(i,j,1:metric_level) =                                  &
                       eg_vert_trisolve(                                     &
                                        am(i,j,2:metric_level),              &
                                        a(i,j,1:metric_level),               &
                                        ap(i,j,1:metric_level-1),            &
                                        b(i,j,1:metric_level),               &
                                        metric_level,1 ) 
              END DO
            END DO

          DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
            DO j = j_sweep_begin, j_sweep_end

              DO i = field_dims_offi%i_start,field_dims_offi_s%i_end 
                G_lambda(i,j,k) = ( field(i+1,j,k) - field(i,j,k) )          &
                                 /( lambda(i+1)    - lambda(i)    )          &
                                 * delta_lambda**2*diff_coeff(i,j)               
              END DO

              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = field(i,j,k)                                   &
                            + (G_lambda(i,j,k) - G_lambda(i-1,j,k))           &
                             /(lambda_offi(i)  - lambda_offi(i-1))
              END DO
            END DO
          END DO

! -----------------------------------------------------------------------------
! --- w fields ---
! --- As for p fields except bottom boundary condition is more complex
          CASE ( fld_type_w ) 
 
            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims_offi%i_start,field_dims_offi_s%i_end

! This is G_lambda stored on (i,j-1/2,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                  IF ( k == field_dims%k_end ) THEN
                    f_av_xz_k = wt_x_on2off(i,j) *field(i,j,k) +            &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) 
                  ELSE
                    f_av_xz_k = wt_z_on2off(i,j,k+1)*(                        &
                                      wt_x_on2off(i,j) *field(i,j,k+1)    + &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1)) + &
                           (1.0-wt_z_on2off(i,j,k+1))*(                       &
                                      wt_x_on2off(i,j) *field(i,j,k) +      &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )
                  END IF

                  f_av_xz_km= wt_z_on2off(i,j,k)*(                        &
                                    wt_x_on2off(i,j) *field(i,j,k) +      &
                               (1.0-wt_x_on2off(i,j))*field(i+1,j,k) ) +  &
                         (1.0-wt_z_on2off(i,j,k))*(                       &
                                    wt_x_on2off(i,j)* field(i,j,k-1) +    &
                               (1.0-wt_x_on2off(i,j))*field(i+1,j,k-1) )

                  G_lambda(i,j,k) = (( field(i+1,j,k) - field(i,j,k)    )     &
                                    *delta_z(i,j,k)                           &
                                   - ( f_av_xz_k - f_av_xz_km   )             &
                                     *(r_levels(i+1,j,k) - r_levels(i,j,k)))  &
                                   /( lambda(i+1) - lambda(i))                &
                                   *delta_lambda*delta_lambda*diff_coeff(i,j) &
                                   *r_levels_offi(i,j,k) * r_levels_offi(i,j,k)
                END DO
              END DO
            END DO
! Bottom boundary
            k = 0
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims_offi%i_start,field_dims_offi_s%i_end
! This is G_lambda stored on (i,j-1/2,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                f_av_xz_k = wt_z_on2off(i,j,k+1)*(                              &
                                  wt_x_on2off(i,j) *field(i,j,k+1)    +         &
                             (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1)) +         &
                       (1.0-wt_z_on2off(i,j,k+1))*(                             &
                                  wt_x_on2off(i,j) *field(i,j,k) +              &
                             (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                f_av_xz_km= wt_x_on2off(i,j) *field(i,j,k) +                    &
                             (1.0-wt_x_on2off(i,j))*field(i+1,j,k)  


                G_lambda(i,j,k) = (( field(i+1,j,k) - field(i,j,k)    )         &
                                 - ( f_av_xz_k - f_av_xz_km   ))                &
                                 *(delta_z(i+1,j,k) + delta_z(i,j,k) )*0.5      &
                                 /( lambda(i+1) - lambda(i))                    &
                                 *delta_lambda*delta_lambda*diff_coeff(i,j)     &
                                 *r_levels_offi(i,j,k) * r_levels_offi(i,j,k)
             END DO
           END DO  

            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_lambda_eta on p points
                  f_av_xz_i = wt_z_on2off(i,j,k)*(                            &
                                wt_x_on2off(i,j)*field(i,j,k)  +              &
                           (1.0-wt_x_on2off(i,j))*field(i+1,j,k)) +           &
                         (1.0-wt_z_on2off(i,j,k))*(                           &
                                wt_x_on2off(i,j)*field(i,j,k-1) +             &
                           (1.0-wt_x_on2off(i,j))*field(i+1,j,k-1) )
                  f_av_xz_im= wt_z_on2off(i,j,k)*(                            &
                                wt_x_on2off(i,j)*field(i-1,j,k)  +            &
                           (1.0-wt_x_on2off(i,j))*field(i,j,k) ) +            &
                         (1.0-wt_z_on2off(i,j,k))*(                           &
                                wt_x_on2off(i,j)*field(i-1,j,k-1) +           &
                           (1.0-wt_x_on2off(i,j))*field(i,j,k-1) )


                  G_lambda_eta(i,j,k) = (( f_av_xz_i - f_av_xz_im ))          &
                                   *0.5*((r_levels_offi(i,j,k)                &
                                                  -r_levels_offi(i-1,j,k))    &
                                       + (r_levels_offi(i,j,k-1)              &
                                                   -r_levels_offi(i-1,j,k-1)))&
                              /( lambda_offi(i) - lambda_offi(i-1))**2        &
                              *delta_lambda*delta_lambda*diff_coeff(i,j)      &
                              *r_levels_offk(i,j,k) * r_levels_offk(i,j,k)
                END DO
              END DO

            END DO
            ! Top & Bottom boundary condition
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                G_lambda_eta(i,j,metric_level+1) = G_lambda_eta(i,j,metric_level)
                G_lambda_eta(i,j,0) = G_lambda_eta(i,j,1)
              END DO
            END DO
            
            k = 0
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)               &
                                            -r_levels_offi(i-1,j,k+1))  &
                                 + (r_levels_offi(i,j,k)                &
                                             -r_levels_offi(i-1,j,k)))  

                ap(i,j,k) = -recip_r_squared_delz(i,j,k)                &
                            *del_r_av_z**2                              &
                            /((lambda_offi(i) - lambda_offi(i-1))**2    &
                            *(r_levels(i,j,k+1)-r_levels(i,j,k)))       &
                            *delta_lambda**2*diff_coeff(i,j)            &
                            *r_levels_offk(i,j,k+1)**2
                am(i,j,k) = 0.0
              END DO
            END DO
            DO k = 1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  ! LHS terms

                  del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)               &
                                              -r_levels_offi(i-1,j,k+1))  &
                                   + (r_levels_offi(i,j,k)                &
                                               -r_levels_offi(i-1,j,k)))  

                  ap(i,j,k) = -recip_r_squared_delz(i,j,k)                &
                              *del_r_av_z**2                              &
                              /((lambda_offi(i) - lambda_offi(i-1))**2    &
                              *(r_levels(i,j,k+1)-r_levels(i,j,k)))       &
                              *delta_lambda**2*diff_coeff(i,j)            &
                              *r_levels_offk(i,j,k+1)**2

                  del_r_av_z = 0.5*((r_levels_offi(i,j,k)                 &
                                              -r_levels_offi(i-1,j,k))    &
                                   + (r_levels_offi(i,j,k-1)              &
                                               -r_levels_offi(i-1,j,k-1)))

                  am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((lambda_offi(i) - lambda_offi(i-1))**2  &
                              *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                              *delta_lambda**2*diff_coeff(i,j)          &
                              *r_levels_offk(i,j,k)**2
                END DO
              END DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                ap(i,j,k) = 0.0

                del_r_av_z = 0.5*((r_levels_offi(i,j,k)                 &
                                            -r_levels_offi(i-1,j,k))    &
                                 + (r_levels_offi(i,j,k-1)              &
                                             -r_levels_offi(i-1,j,k-1)))

                am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((lambda_offi(i) - lambda_offi(i-1))**2  &
                            *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                            *delta_lambda**2*diff_coeff(i,j)          &
                            *r_levels_offk(i,j,k)**2
              END DO
            END DO
            DO k = 0, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  ! LHS terms
                  a(i,j,k) = 1.0 - ap(i,j,k) - am(i,j,k)
                  ! RHS term
                  b(i,j,k) = field(i,j,k) +                                &
                            ( (G_lambda(i,j,k)   - G_lambda(i-1,j,k))      &
                            /( lambda(i+1)       - lambda(i)        )      & 
                            - (G_lambda_eta(i,j,k+1)                       &
                                                 - G_lambda_eta(i,j,k)))   &
                               * recip_r_squared_delz(i,j,k) 
                END DO
              END DO
            END DO
 
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end 
!DEPENDS ON: eg_vert_trisolve
                field(i,j,0:metric_level) =                                  &
                       eg_vert_trisolve(                                     &
                                        am(i,j,1:metric_level),              &
                                        a(i,j,0:metric_level),               &
                                        ap(i,j,0:metric_level-1),            &
                                        b(i,j,0:metric_level),               &
                                        metric_level,0 ) 
              END DO
            END DO

          DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
            DO j = j_sweep_begin, j_sweep_end

              DO i = field_dims_offi%i_start,field_dims_offi_s%i_end 
                G_lambda(i,j,k) = ( field(i+1,j,k) - field(i,j,k) )          &
                                 /( lambda(i+1)    - lambda(i)    )          &
                                 * delta_lambda**2*diff_coeff(i,j)         
              END DO

              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = field(i,j,k)                                  &
                            + (G_lambda(i,j,k) - G_lambda(i-1,j,k))          &
                             /(lambda_offi(i)  - lambda_offi(i-1))
              END DO
            END DO
          END DO
! -----------------------------------------------------------------------------
! --- Default ---
      CASE DEFAULT
        CALL ereport('eg_pofil_zlevel',fld_type, 'Invalid field type')
      END SELECT
! -----------------------------------------------------------------------------
! Flux computations  
! !       flux_EW = 0.0
! !       DO k = field_dims%k_start,metric_level
! !         DO j = j_sweep_begin,j_sweep_end
! ! !           flux_EW = 0.0
! !           DO i = field_dims%i_start,field_dims%i_end
! !             flux_EW = flux_EW + (G_lambda(i+base_i,j,k) -                    &
! !                                            G_lambda(i-1+base_i,j,k))*        &
! !                                    recip_r_squared_delz(1,1,1)         
! ! !                               / (lambda_offi(i+base_i) -                     &
! ! !                                            lambda_offi(i-1+base_i))          & 
! !           END DO
! !         END DO
! !       END DO
! !       WRITE(6,*) 'EW G_lambda flux = ',i_filter,flux_EW
! !       flux_EW = 0.0     
! !       DO j = j_sweep_begin,j_sweep_end
! !         DO i = field_dims%i_start,field_dims%i_end
! ! !           flux_EW = 0.0   
! !           DO k = field_dims%k_start,metric_level
! !             flux_EW = flux_EW - (G_lambda_eta(i,j,k+base_k) -                &
! !                                            G_lambda_eta(i,j,k-1+base_k))*    &
! !                                    recip_r_squared_delz(1,1,1)      
! !           END DO
! !         END DO
! !       END DO
! !       WRITE(6,*) 'EW G_lambda_eta k inner flux = ',i_filter,flux_EW


!------------------------------------------------------------------------
! IF n_procy=1, Northern hemisphere needs to be done after S. Hem
          IF ( L_cycle ) THEN
            j_sweep_store = j_sweep_begin
            IF ( i_sweep == 1 .AND. L_combine ) THEN
              j_sweep_begin = j_sweep_end + 1
            ELSE
              j_sweep_begin = field_dims%j_end - j_sweep_end + 1

              IF ( fld_type == fld_type_v ) &
                  j_sweep_begin = field_dims%j_end - j_sweep_end

            END IF ! i_sweep == 1 .AND. L_combine
            j_sweep_end = field_dims%j_end - j_sweep_store + 1

              IF ( fld_type == fld_type_v ) &
                  j_sweep_end = field_dims%j_end - j_sweep_store

            L_cycle = .false.
        CYCLE
          END IF ! L_cycle)

! Reset pointers for next sweep
! either because N Hem has just been done or
! 1st sweep was combined filter and diffusion
          j_sweep_begin = idx_begin(i_filter)
          j_sweep_end = idx_end(i_filter)
          IF( n_procy == 1 ) L_cycle = .true.

          IF ( fld_type == fld_type_w ) THEN
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
                           field, row_length, in_rows, levels,          &
                           offx, offy, fld_type_p, L_vector) 
          ELSE
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
                           field, row_length, in_rows, levels,          &
                           offx, offy, fld_type, L_vector)
          END IF
          i_sweep = i_sweep + 1
          IF ( i_sweep > sweeps(i_filter) ) EXIT
        CYCLE
          END DO ! sweeping  loop
      END DO   !  i_filter = 1, num_pass


! ----------------------------------------------------------------------
! Section 4.2   NS diffusion
! ----------------------------------------------------------------------

      IF( L_diff ) THEN
        j_sweep_begin = idx_begin(0)
        j_sweep_end = idx_end(0)
        IF ( n_procy == 1 ) THEN
          j_sweep_end = field_dims%j_end - j_sweep_begin + 1

          IF ( fld_type == fld_type_v ) &
            j_sweep_end = field_dims%j_end - j_sweep_begin

        END IF

        G_phi = 0.0
        G_phi_eta = 0.0

!         WRITE(6,*) 'NS Sweep limits ',j_sweep_begin,j_sweep_end
!         flush(6)

        SELECT CASE ( fld_type ) 
! -----------------------------------------------------------------------------
! --- U fields ---
        CASE ( fld_type_u )
        ! settings for diffusion of u
        IF ( L_u_ns_diff ) THEN
          DO k = 1, active_levels
            DO j = j_sweep_begin-1, j_sweep_end+1
              DO i = field_dims%i_start, field_dims%i_end
                field(i,j,k) = field(i,j,k)/csphi_on(j)
              END DO
            END DO
          END DO

          DO j = field_dims_l%j_start,field_dims_l%j_end
            u_fac(j) = csphi_on(j)
            csphi_on_u(j) = csphi_on(j)*csphi_on(j)
          ENDDO
          DO j = field_dims_offj_l%j_start,field_dims_offj_l%j_end
            csphi_off_u(j) = csphi_off(j)*csphi_off(j)*csphi_off(j)
          ENDDO 
        ELSE
          u_fac = 1.0    
          csphi_off_u = csphi_off
          csphi_on_u = csphi_on
        END IF
          DO k = 1, metric_level
            DO j = j_sweep_begin-1, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on (i,j,k-1/2) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                IF ( k == field_dims%k_end ) THEN
                  f_av_yz_k =       wt_y_on2off(i,j) *field(i,j,k) +      &
                               (1.0-wt_y_on2off(i,j))*field(i,j+1,k) 
                ELSE
                  f_av_yz_k = wt_z_on2off(i,j,k)*(                        &
                                    wt_y_on2off(i,j) *field(i,j,k+1)  +   &
                               (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1)) + &
                         (1.0-wt_z_on2off(i,j,k))*(                       &
                                    wt_y_on2off(i,j) *field(i,j,k) +      &
                               (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )
                END IF
                IF ( k == 1 ) THEN
                  f_av_yz_km=  wt_y_on2off(i,j) *field(i,j,k) +         &
                          (1.0-wt_y_on2off(i,j))*field(i,j+1,k) 

                ELSE
                  f_av_yz_km= wt_z_on2off(i,j,k-1)*(                      &
                                    wt_y_on2off(i,j) *field(i,j,k) +      &
                               (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) +  &
                         (1.0-wt_z_on2off(i,j,k-1))*(                     &
                                    wt_y_on2off(i,j)* field(i,j,k-1) +    &
                               (1.0-wt_y_on2off(i,j))*field(i,j+1,k-1) )
                END IF

                G_phi(i,j,k) = (( field(i,j+1,k) - field(i,j,k)    )        &
                                 *delta_z(i,j,k)                            &
                                 - ( f_av_yz_k - f_av_yz_km   )             &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k)))    &
                                 /( phi(j+1) - phi(j))                      &
                                 *delta_phi*delta_phi*diff_coeff_phi(i,j)   &
                                 *r_levels_offj(i,j,k)*r_levels_offj(i,j,k) &
                                 *csphi_off_u(j)
              END DO
            END DO
          END DO
          DO k = 1, metric_level-1
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
! This is G_phi_eta on (i,j-1/2,k) points                    
                f_av_yz_j = wt_z_on2off(i,j,k)*(                           &
                              wt_y_on2off(i,j)*field(i,j,k+1) +            &
                         (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1)) +        &
                       (1.0-wt_z_on2off(i,j,k))*(                          &
                              wt_y_on2off(i,j)*field(i,j,k) +              &
                         (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) 
               f_av_yz_jm= wt_z_on2off(i,j,k)*(                            &
                              wt_y_on2off(i,j)*field(i,j-1,k+1)  +         &
                         (1.0-wt_y_on2off(i,j))*field(i,j,k+1) ) +         &
                       (1.0-wt_z_on2off(i,j,k))*(                          &
                              wt_y_on2off(i,j)*field(i,j-1,k) +            &
                         (1.0-wt_y_on2off(i,j))*field(i,j,k) )

                G_phi_eta(i,j,k) = (( f_av_yz_j - f_av_yz_jm ))            &
                                   *0.5*((r_levels_offj(i,j,k+1)           &
                                                -r_levels_offj(i,j-1,k+1)) &
                                       + (r_levels_offj(i,j,k)             &
                                                 -r_levels_offj(i,j-1,k))) &
                            /( phi_offj(j) - phi_offj(j-1))**2             &
                            *delta_phi*delta_phi*diff_coeff_phi(i,j)       &
                            *r_levels_offk(i,j,k) * r_levels_offk(i,j,k)   &
                            *csphi_on_u(j)
              END DO
            END DO               
          END DO
          ! Top & Bottom boundary condition
          DO j = j_sweep_begin, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end
              G_phi_eta(i,j,metric_level) = 0.0
              G_phi_eta(i,j,0) = 0.0
            END DO
          END DO
          k = 1
          DO j = j_sweep_begin, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end
            ! LHS terms
              del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)             &
                                         -r_levels_offj(i,j-1,k+1)) &
                                + (r_levels_offj(i,j,k)             &
                                          -r_levels_offj(i,j-1,k))) 

              ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                          *del_r_av_z**2                            &
                          /((phi_offj(j) - phi_offj(j-1))**2        &
                          *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                          *delta_phi**2*diff_coeff_phi(i,j)         &
                          *r_levels_offk(i,j,k)**2

              am(i,j,k) = 0.0
            END DO
          END DO
          DO k = 2, metric_level-1
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                ! LHS terms

                del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)             &
                                           -r_levels_offj(i,j-1,k+1)) &
                                  + (r_levels_offj(i,j,k)             &
                                            -r_levels_offj(i,j-1,k))) 

                ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((phi_offj(j) - phi_offj(j-1))**2        &
                            *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                            *delta_phi**2*diff_coeff_phi(i,j)         &
                              *r_levels_offk(i,j,k)**2

                del_r_av_z = 0.5*((r_levels_offj(i,j,k)               &
                                           -r_levels_offj(i,j-1,k))   &
                                  + (r_levels_offj(i,j,k-1)           &
                                            -r_levels_offj(i,j-1,k-1))) 

                am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((phi_offj(j) - phi_offj(j-1))**2        &
                            *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                            *delta_phi**2*diff_coeff_phi(i,j)         &
                            *r_levels_offk(i,j,k-1)**2
              END DO
            END DO
          END DO
          k = metric_level
          DO j = j_sweep_begin, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end
            ! LHS terms
              ap(i,j,k) = 0.0

              del_r_av_z = 0.5*((r_levels_offj(i,j,k)               &
                                         -r_levels_offj(i,j-1,k))   &
                                + (r_levels_offj(i,j,k-1)           &
                                          -r_levels_offj(i,j-1,k-1))) 

              am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                          *del_r_av_z**2                            &
                          /((phi_offj(j) - phi_offj(j-1))**2        &
                          *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                          *delta_phi**2*diff_coeff_phi(i,j)         &
                          *r_levels_offk(i,j,k-1)**2

            END DO
          END DO
          DO k = 1, metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                a(i,j,k) = 1.0*u_fac(j) - ap(i,j,k) - am(i,j,k)
                  ! RHS term
                b(i,j,k) = field(i,j,k)*u_fac(j) +                       &
                          ( (G_phi(i,j,k) - G_phi(i,j-1,k))              &
                           /(phi_offj(j)  - phi_offj(j-1))               &
                          - (G_phi_eta(i,j,k)                            &
                                               - G_phi_eta(i,j,k-1)))    &
                           * recip_r_squared_delz(i,j,k)/csphi_on_u(j) 
              END DO
            END DO
          END DO
          DO j = j_sweep_begin, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end 

!DEPENDS ON: eg_vert_trisolve
              field(i,j,1:metric_level) =                                  &
                     eg_vert_trisolve(                                     &
                                      am(i,j,2:metric_level),              &
                                      a(i,j,1:metric_level),               &
                                      ap(i,j,1:metric_level-1),            &
                                      b(i,j,1:metric_level),               &
                                      metric_level,1 ) 

            END DO
          END DO
        DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
          DO j = j_sweep_begin-1, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end
              G_phi(i,j,k) = ( field(i,j+1,k) - field(i,j,k) )            &
                               /( phi(j+1)    - phi(j)    )               & 
                         * delta_phi**2*diff_coeff_phi(i,j)*csphi_off_u(j)     
            END DO
          END DO

          DO j = j_sweep_begin, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end
              field(i,j,k) = field(i,j,k)                                  &
                          + (G_phi(i,j,k) - G_phi(i,j-1,k))                &
                            /((phi_offj(j)  - phi_offj(j-1))*csphi_on_u(j))  

            END DO
          END DO
        END DO

! -----------------------------------------------------------------------------
! --- V fields ---
          CASE ( fld_type_v ) 
            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end+1
                DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on p points 
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  IF ( k == field_dims%k_end ) THEN
                    f_av_yz_k = wt_y_on2off(i,j) *field(i,j-1,k) +          &
                                 (1.0-wt_y_on2off(i,j))*field(i,j,k) 
                  ELSE
                    f_av_yz_k = wt_z_on2off(i,j,k)*(                        &
                                      wt_y_on2off(i,j) *field(i,j-1,k+1)  + &
                                 (1.0-wt_y_on2off(i,j))*field(i,j,k+1)) +   &
                           (1.0-wt_z_on2off(i,j,k))*(                       &
                                      wt_y_on2off(i,j) *field(i,j-1,k) +    &
                                 (1.0-wt_y_on2off(i,j))*field(i,j,k) )
                  END IF

                  IF ( k == 1 ) THEN
                    f_av_yz_km= wt_y_on2off(i,j) *field(i,j-1,k) +          &
                                 (1.0-wt_y_on2off(i,j))*field(i,j,k)
                  ELSE
                    f_av_yz_km= wt_z_on2off(i,j,k-1)*(                      &
                                      wt_y_on2off(i,j) *field(i,j-1,k) +    &
                                 (1.0-wt_y_on2off(i,j))*field(i,j,k) ) +    &
                           (1.0-wt_z_on2off(i,j,k-1))*(                     &
                                      wt_y_on2off(i,j)* field(i,j-1,k-1) +  &
                                 (1.0-wt_y_on2off(i,j))*field(i,j,k-1) )
                  END IF

                  G_phi(i,j,k) = (( field(i,j,k) - field(i,j-1,k)    )      &
                                   *delta_z(i,j,k)                          &
                                   - ( f_av_yz_k - f_av_yz_km   )           &
                                   *(r_levels(i,j,k)-r_levels(i,j-1,k)))    &
                                   /( phi(j) - phi(j-1))                    &
                                  *delta_phi*delta_phi*diff_coeff_phi(i,j)  &
                                   *r_levels_offj(i,j,k)*r_levels_offj(i,j,k)&
                                  *csphi_off(j)
                END DO
              END DO              
            END DO
            DO k = 1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_phi_eta on (i-1/2,j,k) points
                  f_av_yz_j = wt_z_on2off(i,j,k)*(                          &
                                wt_y_on2off(i,j)*field(i,j,k+1-base_k)  +   &
                           (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1-base_k)) +&
                         (1.0-wt_z_on2off(i,j,k))*(                         &
                                wt_y_on2off(i,j)*field(i,j,k-base_k) +      &
                           (1.0-wt_y_on2off(i,j))*field(i,j+1,k-base_k) )
                  f_av_yz_jm= wt_z_on2off(i,j,k)*(                          &
                                wt_y_on2off(i,j)*field(i,j-1,k+1-base_k)  + &
                           (1.0-wt_y_on2off(i,j))*field(i,j,k+1-base_k) ) + &
                         (1.0-wt_z_on2off(i,j,k))*(                         &
                                wt_y_on2off(i,j)*field(i,j-1,k-base_k) +    &
                           (1.0-wt_y_on2off(i,j))*field(i,j,k-base_k) )

                  G_phi_eta(i,j,k) = (( f_av_yz_j - f_av_yz_jm ))           &
                                  *0.5*((r_levels_offj(i,j+1,k+1)           &
                                                   -r_levels_offj(i,j,k+1)) &
                                      + (r_levels_offj(i,j+1,k)             &
                                                   -r_levels_offj(i,j,k)))  &
                              /( phi_offj(j+1) - phi_offj(j))**2            &
                              *delta_phi*delta_phi*diff_coeff_phi(i,j)      &
                              *r_levels_offk(i,j,k) * r_levels_offk(i,j,k)  &
                              *csphi_on(j)
                END DO
              END DO
            END DO
            ! Top & Bottom boundary condition
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                G_phi_eta(i,j,metric_level) = 0.0
                G_phi_eta(i,j,0) = 0.0
              END DO
            END DO
            k = 1
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                del_r_av_z = 0.5*((r_levels_offj(i,j+1,k+1)           &
                                             -r_levels_offj(i,j,k+1)) &
                                + (r_levels_offj(i,j+1,k)             &
                                             -r_levels_offj(i,j,k)))  

                ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((phi_offj(j+1) - phi_offj(j))**2        &
                            *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                            *delta_phi**2*diff_coeff_phi(i,j)         &
                             *r_levels_offk(i,j,k)**2

                am(i,j,k) = 0.0
              END DO
            END DO
            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                 ! LHS terms

                  del_r_av_z = 0.5*((r_levels_offj(i,j+1,k+1)           &
                                               -r_levels_offj(i,j,k+1)) &
                                  + (r_levels_offj(i,j+1,k)             &
                                               -r_levels_offj(i,j,k)))  

                  ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((phi_offj(j+1) - phi_offj(j))**2        &
                              *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                              *delta_phi**2*diff_coeff_phi(i,j)         &
                               *r_levels_offk(i,j,k)**2

                  del_r_av_z = 0.5*((r_levels_offj(i,j+1,k)             &
                                               -r_levels_offj(i,j,k))   &
                                  + (r_levels_offj(i,j+1,k-1)           &
                                               -r_levels_offj(i,j,k-1)))  

                  am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((phi_offj(j+1) - phi_offj(j))**2        &
                              *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                              *delta_phi**2*diff_coeff_phi(i,j)         &
                              *r_levels_offk(i,j,k-1)**2
                END DO
              END DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                ap(i,j,k) = 0.0

                del_r_av_z = 0.5*((r_levels_offj(i,j+1,k)             &
                                             -r_levels_offj(i,j,k))   &
                                + (r_levels_offj(i,j+1,k-1)           &
                                             -r_levels_offj(i,j,k-1)))  

                am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((phi_offj(j+1) - phi_offj(j))**2        &
                            *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                            *delta_phi**2*diff_coeff_phi(i,j)         &
                            *r_levels_offk(i,j,k-1)**2

              END DO
            END DO
            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                 ! LHS terms
                  a(i,j,k) = 1.0 - ap(i,j,k) - am(i,j,k)
                 ! RHS term
                  b(i,j,k) = field(i,j,k) +                                &
                            ( (G_phi(i,j+1,k)   - G_phi(i,j,k))            &
                             /(phi_offj(j+1)    - phi_offj(j))             &
                            - (G_phi_eta(i,j,k)                            &
                                                 - G_phi_eta(i,j,k-1)))    &
                               * recip_r_squared_delz(i,j,k)/csphi_on(j)
                END DO
              END DO
            END DO
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end 
!DEPENDS ON: eg_vert_trisolve
                field(i,j,1:metric_level) =                                  &
                       eg_vert_trisolve(                                     &
                                        am(i,j,2:metric_level),              &
                                        a(i,j,1:metric_level),               &
                                        ap(i,j,1:metric_level-1),            &
                                        b(i,j,1:metric_level),               &
                                        metric_level,1 ) 
              END DO
            END DO
          DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
            DO j = j_sweep_begin, j_sweep_end+1

              DO i = field_dims%i_start,field_dims%i_end 
                G_phi(i,j,k) = ( field(i,j,k) - field(i,j-1,k) )          &
                                 /( phi(j)    - phi(j-1)    )             & 
                                 * delta_phi**2*diff_coeff_phi(i,j)       &
                                 *csphi_off(j)
              END DO
            END DO

            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = field(i,j,k)                               &
                            + (G_phi(i,j+1,k) - G_phi(i,j,k))             &
                             /((phi_offj(j+1)  - phi_offj(j))/csphi_on(j))
              END DO
            END DO
          END DO

! -----------------------------------------------------------------------------
! --- P fields ---
          CASE ( fld_type_p )   
            DO k = 1, metric_level
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on v points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  IF ( k == field_dims%k_end ) THEN
                    f_av_yz_k = wt_y_on2off(i,j) *field(i,j,k) +            &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k) 
                  ELSE
                    f_av_yz_k = wt_z_on2off(i,j,k)*(                        &
                                      wt_y_on2off(i,j) *field(i,j,k+1)    + &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1)) + &
                           (1.0-wt_z_on2off(i,j,k))*(                       &
                                      wt_y_on2off(i,j) *field(i,j,k) +      &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )
                  END IF
                  IF ( k == 1 ) THEN
                    f_av_yz_km= wt_y_on2off(i,j) *field(i,j,k) +            &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k)
                  ELSE
                    f_av_yz_km= wt_z_on2off(i,j,k-1)*(                      &
                                      wt_y_on2off(i,j) *field(i,j,k) +      &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) +  &
                           (1.0-wt_z_on2off(i,j,k-1))*(                     &
                                      wt_y_on2off(i,j)* field(i,j,k-1) +    &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k-1) )
                  END IF
                  G_phi(i,j,k) = (( field(i,j+1,k) - field(i,j,k) )         &
                                *delta_z(i,j,k)                             &
                                - ( f_av_yz_k      - f_av_yz_km   )         &
                                *(r_levels(i,j+1,k) - r_levels(i,j,k)))     &
                                /( phi(j+1) - phi(j))                       &
                                *delta_phi*delta_phi*diff_coeff_phi(i,j)    &
                                *r_levels_offj(i,j,k)*r_levels_offj(i,j,k)  &
                                *csphi_off(j)
                END DO
              END DO
            END DO

            DO k = 1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_phi_eta on w points
                  f_av_yz_j = wt_z_on2off(i,j,k)*(                          &
                                wt_y_on2off(i,j)*field(i,j,k+1-base_k)  +   &
                           (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1-base_k)) +&
                         (1.0-wt_z_on2off(i,j,k))*(                         &
                                wt_y_on2off(i,j)*field(i,j,k-base_k) +      &
                           (1.0-wt_y_on2off(i,j))*field(i,j+1,k-base_k) )
                  f_av_yz_jm= wt_z_on2off(i,j,k)*(                          &
                                wt_y_on2off(i,j)*field(i,j-1,k+1-base_k)  + &
                           (1.0-wt_y_on2off(i,j))*field(i,j,k+1-base_k) ) + &
                         (1.0-wt_z_on2off(i,j,k))*(                         &
                                wt_y_on2off(i,j)*field(i,j-1,k-base_k) +    &
                           (1.0-wt_y_on2off(i,j))*field(i,j,k-base_k) )


                  G_phi_eta(i,j,k) = (( f_av_yz_j - f_av_yz_jm ))           &
                                  *0.5*((r_levels_offj(i,j,k+1)             &
                                                 -r_levels_offj(i,j-1,k+1)) &
                                      + (r_levels_offj(i,j,k)               &
                                                 -r_levels_offj(i,j-1,k)))  &
                              /( phi_offj(j) - phi_offj(j-1))**2            &
                              *delta_phi*delta_phi*diff_coeff_phi(i,j)      &
                              *r_levels_offk(i,j,k) * r_levels_offk(i,j,k)  &
                              *csphi_on(j)
                END DO
              END DO
            END DO
            !Top & Bottom boundary condition
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                G_phi_eta(i,j,metric_level) = 0.0
                G_phi_eta(i,j,0) = 0.0
              END DO
            END DO

            k = 1
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)             &
                                           -r_levels_offj(i,j-1,k+1)) &
                                + (r_levels_offj(i,j,k)               &
                                           -r_levels_offj(i,j-1,k)))  

                ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((phi_offj(j) - phi_offj(j-1))**2        &
                            *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                            *delta_phi**2*diff_coeff_phi(i,j)         &
                            *r_levels_offk(i,j,k)**2

                am(i,j,k) = 0.0
              END DO
            END DO
            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  ! LHS terms

                  del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)             &
                                             -r_levels_offj(i,j-1,k+1)) &
                                  + (r_levels_offj(i,j,k)               &
                                             -r_levels_offj(i,j-1,k)))  

                  ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((phi_offj(j) - phi_offj(j-1))**2        &
                              *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                              *delta_phi**2*diff_coeff_phi(i,j)         &
                              *r_levels_offk(i,j,k)**2

                  del_r_av_z = 0.5*((r_levels_offj(i,j,k)               &
                                             -r_levels_offj(i,j-1,k))   &
                                  + (r_levels_offj(i,j,k-1)             &
                                             -r_levels_offj(i,j-1,k-1)))  

                  am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((phi_offj(j) - phi_offj(j-1))**2        &
                              *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                              *delta_phi**2*diff_coeff_phi(i,j)         &
                              *r_levels_offk(i,j,k-1)**2
                END DO
              END DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                ap(i,j,k) = 0.0 

                del_r_av_z = 0.5*((r_levels_offj(i,j,k)               &
                                           -r_levels_offj(i,j-1,k))   &
                                + (r_levels_offj(i,j,k-1)             &
                                           -r_levels_offj(i,j-1,k-1)))  

                am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((phi_offj(j) - phi_offj(j-1))**2        &
                            *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                            *delta_phi**2*diff_coeff_phi(i,j)         &
                            *r_levels_offk(i,j,k-1)**2              
              END DO
            END DO
            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                ! LHS terms
                  a(i,j,k) = 1.0 - ap(i,j,k) - am(i,j,k)
                  ! RHS term
                  b(i,j,k) = field(i,j,k) +                                &
                            ( (G_phi(i,j,k)   - G_phi(i,j-1,k))            &
                            /( phi_offj(j)    - phi_offj(j-1) )            & 
                            - (G_phi_eta(i,j,k)                            &
                                                 - G_phi_eta(i,j,k-1)))    &
                            * recip_r_squared_delz(i,j,k)/csphi_on(j) 
                END DO
              END DO
            END DO

            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end 
!DEPENDS ON: eg_vert_trisolve
                field(i,j,1:metric_level) =                                  &
                       eg_vert_trisolve(                                     &
                                        am(i,j,2:metric_level),              &
                                        a(i,j,1:metric_level),               &
                                        ap(i,j,1:metric_level-1),            &
                                        b(i,j,1:metric_level),               &
                                        metric_level,1 ) 
              END DO
            END DO

          DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
            DO j = j_sweep_begin-1, j_sweep_end

              DO i = field_dims%i_start,field_dims%i_end 
                G_phi(i,j,k) = ( field(i,j+1,k) - field(i,j,k) )          &
                                 /( phi(j+1)    - phi(j)    )             &
                                 * delta_phi**2*diff_coeff_phi(i,j)
              END DO
            END DO

            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = field(i,j,k)                               &
                            + (G_phi(i,j,k) - G_phi(i,j-1,k))             &
                             /(phi_offj(j)  - phi_offj(j-1))
              END DO
            END DO
          END DO

! -----------------------------------------------------------------------------
! --- w fields ---
! --- As for p fields except bottom boundary condition is more complex
          CASE ( fld_type_w ) 
            DO k = 1, metric_level
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on (i-1/2,j,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                  IF ( k == field_dims%k_end ) THEN
                    f_av_yz_k = wt_y_on2off(i,j) *field(i,j,k) +            &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k) 
                  ELSE
                    f_av_yz_k = wt_z_on2off(i,j,k+1)*(                      &
                                      wt_y_on2off(i,j) *field(i,j,k+1)    + &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1)) + &
                           (1.0-wt_z_on2off(i,j,k+1))*(                     &
                                      wt_y_on2off(i,j) *field(i,j,k) +      &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )
                  END IF

                  f_av_yz_km= wt_z_on2off(i,j,k)*(                        &
                                    wt_y_on2off(i,j) *field(i,j,k) +      &
                               (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) +  &
                         (1.0-wt_z_on2off(i,j,k))*(                       &
                                    wt_y_on2off(i,j)* field(i,j,k-1) +    &
                               (1.0-wt_y_on2off(i,j))*field(i,j+1,k-1) )

                  G_phi(i,j,k) = (( field(i,j+1,k) - field(i,j,k)    )      &
                                 *delta_z(i,j,k)                            &
                                 - ( f_av_yz_k - f_av_yz_km   )             &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k)))    &
                                 /( phi(j+1) - phi(j))                      &
                                 *delta_phi*delta_phi*diff_coeff_phi(i,j)   &
                                 *r_levels_offj(i,j,k)*r_levels_offj(i,j,k) &
                                 *csphi_off(j)
                END DO
              END DO
            END DO
! Bottom boundary
            k = 0
            DO j = j_sweep_begin-1, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
! This is G_phi stored on (i-1/2,j,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                f_av_yz_k = wt_z_on2off(i,j,k+1)*(                          &
                                  wt_y_on2off(i,j) *field(i,j,k+1)    +     &
                             (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1)) +     &
                       (1.0-wt_z_on2off(i,j,k+1))*(                         &
                                  wt_y_on2off(i,j) *field(i,j,k) +          &
                             (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                f_av_yz_km= wt_y_on2off(i,j) *field(i,j,k) +                &
                             (1.0-wt_x_on2off(i,j))*field(i,j+1,k)  


                G_phi(i,j,k) = (( field(i+1,j,k) - field(i,j,k)    )        &
                               *delta_z(i,j,k)                              &
                               - ( f_av_yz_k - f_av_yz_km   )               &
                               *(r_levels(i,j+1,k) - r_levels(i,j,k)))      &
                               /( phi(j+1) - phi(j) )                       &
                               *delta_phi*delta_phi*diff_coeff_phi(i,j)     &
                               *r_levels_offj(i,j,k) * r_levels_offj(i,j,k) &
                               *csphi_off(j)
             END DO
           END DO  

            DO k = 1, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_phi_eta on p points
                  f_av_yz_j = wt_z_on2off(i,j,k)*(                          &
                                wt_y_on2off(i,j)*field(i,j,k)  +            &
                           (1.0-wt_y_on2off(i,j))*field(i,j+1,k)) +         &
                         (1.0-wt_z_on2off(i,j,k))*(                         &
                                wt_y_on2off(i,j)*field(i,j,k-1) +           &
                           (1.0-wt_y_on2off(i,j))*field(i,j+1,k-1) )
                  f_av_yz_jm= wt_z_on2off(i,j,k)*(                          &
                                wt_y_on2off(i,j)*field(i,j-1,k)  +          &
                           (1.0-wt_y_on2off(i,j))*field(i,j,k) ) +          &
                         (1.0-wt_z_on2off(i,j,k))*(                         &
                                wt_y_on2off(i,j)*field(i,j-1,k-1) +         &
                           (1.0-wt_y_on2off(i,j))*field(i,j,k-1) )


                  G_phi_eta(i,j,k) = (( f_av_yz_j - f_av_yz_jm ))           &
                                  *0.5*((r_levels_offj(i,j,k)               &
                                                -r_levels_offj(i,j-1,k))    &
                                      + (r_levels_offj(i,j,k-1)             &
                                                -r_levels_offj(i,j-1,k-1))) &
                              /( phi_offj(j) - phi_offj(j-1))**2            &
                              *delta_phi*delta_phi*diff_coeff_phi(i,j)      &
                              *r_levels_offk(i,j,k) * r_levels_offk(i,j,k)  &
                              *csphi_on(j)
                END DO
              END DO

            END DO
            ! Top & Bottom boundary condition
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                G_phi_eta(i,j,metric_level+1) = G_phi_eta(i,j,metric_level)
                G_phi_eta(i,j,0) = G_phi_eta(i,j,1)
              END DO
            END DO
           
            k = 0
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)             &
                                          -r_levels_offj(i,j-1,k+1))  &
                                + (r_levels_offj(i,j,k)               &
                                          -r_levels_offj(i,j-1,k)))   

                ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((phi_offj(j) - phi_offj(j-1))**2        &
                            *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                            *delta_phi**2*diff_coeff_phi(i,j)         &
                            *r_levels_offk(i,j,k+1)**2


                am(i,j,k) = 0.0
              END DO
            END DO
            DO k = 1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  ! LHS terms
                  del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)             &
                                            -r_levels_offj(i,j-1,k+1))  &
                                  + (r_levels_offj(i,j,k)               &
                                            -r_levels_offj(i,j-1,k)))   

                  ap(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((phi_offj(j) - phi_offj(j-1))**2        &
                              *(r_levels(i,j,k+1)-r_levels(i,j,k)))     &
                              *delta_phi**2*diff_coeff_phi(i,j)         &
                              *r_levels_offk(i,j,k+1)**2

                  del_r_av_z = 0.5*((r_levels_offj(i,j,k)               &
                                            -r_levels_offj(i,j-1,k))    &
                                  + (r_levels_offj(i,j,k-1)             &
                                            -r_levels_offj(i,j-1,k-1)))   

                  am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                              *del_r_av_z**2                            &
                              /((phi_offj(j) - phi_offj(J-1))**2        &
                              *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                              *delta_phi**2*diff_coeff_phi(i,j)         &
                              *r_levels_offk(i,j,k)**2
                END DO
              END DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
              ! LHS terms
                ap(i,j,k) = 0.0

                del_r_av_z = 0.5*((r_levels_offj(i,j,k)               &
                                          -r_levels_offj(i,j-1,k))    &
                                + (r_levels_offj(i,j,k-1)             &
                                          -r_levels_offj(i,j-1,k-1)))   

                am(i,j,k) = -recip_r_squared_delz(i,j,k)              &
                            *del_r_av_z**2                            &
                            /((phi_offj(j) - phi_offj(J-1))**2        &
                            *(r_levels(i,j,k)-r_levels(i,j,k-1)))     &
                            *delta_phi**2*diff_coeff_phi(i,j)         &
                            *r_levels_offk(i,j,k)**2

              END DO
            END DO
            DO k = 0, metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  ! LHS terms
                  a(i,j,k) = 1.0 - ap(i,j,k) - am(i,j,k)
                  ! RHS term
                  b(i,j,k) = field(i,j,k) +                                &
                            ( (G_phi(i,j,k)   - G_phi(i,j-1,k))            &
                            /( phi(j+1)       - phi(j)        )            & 
                            - (G_phi_eta(i,j,k+1)                          &
                                                 - G_phi_eta(i,j,k)))      &
                               * recip_r_squared_delz(i,j,k)/csphi_on(j) 
                END DO
              END DO
            END DO

            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end 
!DEPENDS ON: eg_vert_trisolve
                field(i,j,0:metric_level) =                                 &
                       eg_vert_trisolve(                                    &
                                        am(i,j,1:metric_level),             &
                                        a(i,j,0:metric_level),              &
                                        ap(i,j,0:metric_level-1),           &
                                        b(i,j,0:metric_level),              &
                                        metric_level,0 ) 
              END DO
            END DO

          DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
            DO j = j_sweep_begin-1, j_sweep_end

              DO i = field_dims%i_start,field_dims%i_end 
                G_phi(i,j,k) = ( field(i,j+1,k) - field(i,j,k) )            &
                                 /( phi(j+1)    - phi(j)    )               &
                                 * delta_phi**2*diff_coeff_phi(i,j)         &
                                 *csphi_off(j)
              END DO
            END DO
         
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = field(i,j,k)                                 &
                            + (G_phi(i,j,k) - G_phi(i,j-1,k))               &
                             /((phi_offj(j)  - phi_offj(j-1))*csphi_on(j))
              END DO
            END DO
          END DO
! -----------------------------------------------------------------------------
! --- Default ---
      CASE DEFAULT
        CALL ereport('eg_pofil_zlevel',fld_type, 'Invalid field type')
      END SELECT

      IF ( fld_type == fld_type_w ) THEN
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                             &
                         field, row_length, in_rows, levels,          &
                         offx, offy, fld_type_p, L_vector) 
      ELSE
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                             &
                         field, row_length, in_rows, levels,          &
                         offx, offy, fld_type, L_vector)
      END IF

! -----------------------------------------------------------------------------
! Flux computations  
! !       flux_EW = 0.0
! !       DO k = field_dims%k_start,metric_level
! !         DO j = j_sweep_begin,j_sweep_end
! ! !           flux_EW = 0.0
! !           DO i = field_dims%i_start,field_dims%i_end
! !             flux_EW = flux_EW + (G_phi(i,j+base_j,k) -                    &
! !                                            G_phi(i,j-1+base_j,k))*        &
! !                                    recip_r_squared_delz(1,1,1)         
! ! !                               / (lambda_offi(i+base_i) -                     &
! ! !                                            lambda_offi(i-1+base_i))          & 
! !           END DO
! !         END DO
! !       END DO
! !       WRITE(6,*) 'NS G_phi flux = ',flux_EW
! !       flux_EW = 0.0     
! !       DO j = j_sweep_begin,j_sweep_end
! !         DO i = field_dims%i_start,field_dims%i_end
! ! !           flux_EW = 0.0   
! !           DO k = field_dims%k_start,metric_level
! !             flux_EW = flux_EW - (G_phi_eta(i,j,k+base_k) -                &
! !                                            G_phi_eta(i,j,k-1+base_k))*    &
! !                                    recip_r_squared_delz(1,1,1)      
! !           END DO
! !         END DO
! !       END DO
! !       WRITE(6,*) 'NS G_phi_eta k inner flux = ',flux_EW

END IF ! L_diff


! !       DO k = 1, active_levels
! ! 
! !         IF ( k <= metric_level ) THEN
! ! !   No metric terms since levels are horizontal
! ! 
! !           DO j = j_sweep_begin-1, j_sweep_end
! !             DO i = 1, row_length
! !               temp(i,j) = (field(i,j+1,k) * delta_z(i,j+1,k) -          &
! !                            field(i,j,  k) * delta_z(i,j, k) )           &
! !                                           * diff_coeff_phi              &
! !                            *r_midj(i,j+off_v,k) * r_midj(i,j+off_v,k)   &
! !                            *csphi_off(j)
! !             END DO
! !           END DO
! ! 
! !           DO j = j_sweep_begin, j_sweep_end            
! !             DO i = 1, row_length              
! !               field(i,j,k) = field(i,j,k)*u_fac(j+v_shift-1)            &
! !                            + (temp(i,j) - temp(i,j-1))                  &
! !                            *1.0/csphi_on(j+v_shift-1)                   &
! !                            *recip_r_squared_delz(i,j,k)
! !             END DO
! !           END DO
! ! 
! !         ELSE  ! k > metric_level
! ! !   No metric terms since levels are horizontal
! ! 
! !           DO j = j_sweep_begin-1, j_sweep_end
! !             DO i = 1, row_length
! !               temp(i,j) = (field(i,j+1,k) - field(i,j,k) ) *            &
! !                            diff_coeff_phi * csphi_off(j)                
! !             END DO
! !           END DO
! ! 
! !           DO j = j_sweep_begin, j_sweep_end
! !             DO i = 1, row_length
! !               field(i,j,k) = field(i,j,k)*u_fac(j+v_shift-1)            & 
! !                            + (temp(i,j) - temp(i,j-1))                  &
! !                              *1.0/csphi_on(j+v_shift-1)
! !             END DO
! !           END DO
! ! 
! !         END IF ! k <= metric_level
! ! 
! !       END DO  ! k = 1, active_levels

! ! ! DEPENDS ON: swap_bounds
! !         CALL Swap_Bounds(                                               &
! !                          field, row_length, in_rows, levels,            &
! !                          halo_x, halo_y, fld_type, L_vector)
! ! 
! !       END IF ! L_diff

   
      IF (lhook) CALL dr_hook('EG_POFIL_ZLEVEL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eg_pofil_zlevel
      END MODULE  eg_pofil_zlevel_mod
