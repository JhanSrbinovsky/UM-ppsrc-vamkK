! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine pofil_zlevel
!
! Purpose:
!          Filter/diffusion based on first order
!          conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.
!          Applies z-level diffusion (truly horiztonal
!          diffusion)
!
! Method:
!          For any field. Pointers used to input grid
!          This version tested for u-at-the-poles new dynamics
!          (a13_2a) and has extensions for future ENDGame 
!          implentation (A13_4A), (not currently tested). 
!          There is currently no support for v at poles ND (A13_3A)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
      MODULE  pofil_zlevel_mod
        IMPLICIT NONE
      CONTAINS
      SUBROUTINE pofil_zlevel(                                           &
                            field, fld_type, field_dims,                 & 
                            field_dims_s, field_dims_l,                  &
                            field_dims_offj, field_dims_offj_s,          &
                            field_dims_offj_l, field_dims_offi,          &
                            field_dims_offi_s, field_dims_offi_l,        &
                            row_length, in_rows, levels, offx, offy,     &
                            delta_lambda, delta_phi, active_levels,      &
                            metric_level,                                &
                            r_levels, r_levels_offk,                     &
                            r_levels_offi, r_levels_offj,                &
                            eta_levels, eta_levels_offk,                 &
                            lambda, lambda_offi, phi, phi_offj,          &
                            n_procy, max_filter_rows, global_filter,     &
                            sweeps, sweep_begin, sweep_end,              &
                            diff_coeff_phi, diff_coeff, L_diff,          &
                            L_vector,                                    &
                            csphi_on, csphi_off,L_u_ns_diff,kl,          &
                            at_extremity, r_dim_k_start)

      USE atm_fields_bounds_mod
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      USE ereport_mod, ONLY : ereport

      IMPLICIT NONE


! Arguments with Intent IN. ie: Input variables.

      TYPE (array_dims) :: field_dims, field_dims_s, field_dims_l,       &
                           field_dims_offi, field_dims_offi_s,           &
                           field_dims_offi_l,                            &
                           field_dims_offj, field_dims_offj_s,           &
                           field_dims_offj_l                            

      LOGICAL, Intent(In) ::                                             &
        L_diff,                                                          &
                         ! true IF diffusion active
        L_vector         ! true IF vector field


      LOGICAL, Intent(In) :: at_extremity(4)  
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
      LOGICAL, Intent(In) :: L_u_ns_diff ! true if diffusing u

                             
      INTEGER, Intent(In) ::                                             &
        max_filter_rows,                                                 &
                         ! max dimension for sweeping arrays
        global_filter,                                                   &
                         ! max number of filter sweeps 0=diffusion only 
        fld_type,                                                        &
                      ! field type (p=1, u=2 or v=3)
        n_procy,                                                         &
                       ! Number of processors in latitude
        active_levels,                                                   &
                          ! number of levels to be filtered
        metric_level,                                                    &
                        ! uppermost non-constant layer for metric terms
        kl,                                                              &
                        ! Start level for w(=1) and theta (=0) diffusion
        r_dim_k_start
                        ! k_start index for r dimensions
                        ! r_dim_k_start = field_dims%k_start for EG
                        ! r_dim_k_start /= field_dims%k_start for w/theta ND fields

      REAL, Intent(In) :: delta_lambda,  delta_phi

! field dimensions for swap bounds
      INTEGER, Intent(In) :: row_length, in_rows, levels, offx, offy


! Number of sweeps and NS sweeping limits
      INTEGER, Intent(In) ::                                             &
        sweeps   (  max_filter_rows),                                    &
        sweep_begin(0:max_filter_rows),                                  &
        sweep_end  (0:max_filter_rows)                       

      REAL, Intent(In) ::                                                &
! EW diffusion coefficient
                 diff_coeff(field_dims_s%i_start:field_dims_s%i_end,     &
                            field_dims_s%j_start:field_dims_s%j_end),    &
! NS diffusion coefficient
             diff_coeff_phi(field_dims_s%i_start:field_dims_s%i_end,     &
                            field_dims_s%j_start:field_dims_s%j_end)

! r level arrays
      REAL, Intent(In) ::                                                &
        r_levels (field_dims_l%i_start:field_dims_l%i_end,               &
                  field_dims_l%j_start:field_dims_l%j_end,               &
                  r_dim_k_start       :field_dims_l%k_end),              &
        r_levels_offk(field_dims_l%i_start:field_dims_l%i_end,           &
                      field_dims_l%j_start:field_dims_l%j_end,           &
                    1-r_dim_k_start       :field_dims_l%k_end),          &
        r_levels_offi(field_dims_offi_l%i_start:field_dims_offi_l%i_end, &
                      field_dims_l%j_start     :field_dims_l%j_end,      &
                      r_dim_k_start            :field_dims_l%k_end),     &
        r_levels_offj(field_dims_l%i_start     :field_dims_l%i_end,      &
                      field_dims_offj_l%j_start:field_dims_offj_l%j_end, &
                      r_dim_k_start            :field_dims_l%k_end)
! eta level arrays
      REAL, Intent(In) ::                                                &
        eta_levels     (  r_dim_k_start:field_dims_l%k_end),             &
        eta_levels_offk(1-r_dim_k_start:field_dims_l%k_end)

! Horizontal grid arrays
      REAL, Intent(In) ::                                                &
        lambda     (field_dims_l%i_start:field_dims_l%i_end),            &
        phi        (field_dims_l%j_start:field_dims_l%j_end),            &
        lambda_offi(field_dims_offi_l%i_start:field_dims_offi_l%i_end),  &
        phi_offj   (field_dims_offj_l%j_start:field_dims_offj_l%j_end)  

      REAL, Intent(In)    ::                                             &
           csphi_on(field_dims_l%j_start:field_dims_l%j_end),            &
           csphi_off(field_dims_offj_l%j_start:field_dims_offj_l%j_end)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
      REAL, Intent(InOut) ::                                             &
           field (field_dims_s%i_start:field_dims_s%i_end,               &
                  field_dims_s%j_start:field_dims_s%j_end,               &
                  field_dims_s%k_start:field_dims_s%k_end)


! Local Variables.
      INTEGER :: i, j, k, kp, i_filter, i_sweep, num_pass,               &
                 j_sweep_begin, j_sweep_end, j_sweep_store

      INTEGER :: base_k, base_i, base_j

      LOGICAL :: L_cycle, L_combine

! Local copies of NS sweeping limits
      INTEGER             ::                                             &
        idx_begin(0:max_filter_rows),                                    &
        idx_end  (0:max_filter_rows)  

! Local copies of cos(phi) for NS U diffusion
      REAL  ::                                                           &
           csphi_on_u(field_dims_l%j_start:field_dims_l%j_end),          &
           csphi_off_u(field_dims_offj_l%j_start:                        &
                       field_dims_offj_l%j_end),                         &
           u_fac(field_dims_l%j_start:field_dims_l%j_end)

      REAL ::                                                            &
        recip_r_squared_delz(field_dims_s%i_start:field_dims_s%i_end,    & 
                             field_dims_s%j_start:field_dims_s%j_end,    &
                             field_dims_s%k_start:metric_level)

      REAL :: del_z, diff_coeff_av

      REAL ::                                                            &
        wt_z_on2off(field_dims_s%i_start:field_dims_s%i_end,             &
                    field_dims_s%j_start:field_dims_s%j_end,             &
                    0:field_dims_s%k_end),                               &
        wt_x_on2off(field_dims_offi_s%i_start:field_dims_offi_s%i_end,   &
                    field_dims_s%j_start     :field_dims_s%j_end),       &
        wt_y_on2off(field_dims_s%i_start     :field_dims_s%i_end,        &
                    field_dims_offj_s%j_start:field_dims_offj_s%j_end) 

      REAL :: G_lambda(field_dims_offi_s%i_start:                        &
                                             field_dims_offi_s%i_end,    &
                       field_dims_s%j_start     :field_dims_s%j_end,     &
                       field_dims_s%k_start     :field_dims_s%k_end),    &
              G_phi   (field_dims_s%i_start     :field_dims_s%i_end,     &
                       field_dims_offj_s%j_start:                        &
                                             field_dims_offj_s%j_end,    &
                       field_dims_s%k_start     :field_dims_s%k_end)

      REAL ::  G_lambda_eta_k_2D(field_dims%i_start:field_dims%i_end,    &
                                 field_dims%j_start:field_dims%j_end),   &
               G_phi_eta_k_2D(field_dims%i_start:field_dims%i_end,       &
                              field_dims%j_start:field_dims%j_end)

      REAL :: G_lambda_eta_km, G_phi_eta_km

      REAL :: f_av_xz_k_2D(field_dims_offi_s%i_start:                    & 
                                             field_dims_offi_s%i_end,    &
                           field_dims%j_start:field_dims%j_end),         &
              f_av_yz_k_2D(field_dims%i_start:field_dims%i_end,          &
                           field_dims_s%j_start:field_dims_s%j_end),     &
              f_av_yz_j_1D(field_dims%i_start:field_dims%i_end)

      REAL :: f_av_xz_km, f_av_xz_i, f_av_xz_im,                         &
              f_av_yz_km, f_av_yz_jm, del_r_av_z

      INTEGER :: ND_shift ! Shift in index for New dynamics

      INTEGER :: j_v_b, j_v_e

! Off-centring for implicit terms
      REAL :: alpha, beta

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Extra variables for optimisation
! Arrays used for trisolve
      REAL ::    b_3D(field_dims%i_start:field_dims%i_end,               &
                      field_dims%j_start:field_dims%j_end,               &
                      field_dims%k_start:metric_level),                  &
                ap_3D(field_dims%i_start:field_dims%i_end,               &
                      field_dims%j_start:field_dims%j_end,               & 
                      1:metric_level-1),                                 &
                am_3D(field_dims%i_start:field_dims%i_end,               &
                      field_dims%j_start:field_dims%j_end,               &
                      2:metric_level),                                   &
                denom(field_dims%i_start:field_dims%i_end,               &
                      field_dims%j_start:field_dims%j_end,               &
                      1:metric_level),                                   &
                c_new(field_dims%i_start:field_dims%i_end,               &
                      field_dims%j_start:field_dims%j_end,               &
                      1:metric_level-1)  

      LOGICAL :: L_init_matrix

      REAL ::  recip_del_lambda(field_dims%i_start:field_dims%i_end),    &
               diff_coeff_av_2D(field_dims%i_start:field_dims%i_end,     &
                                field_dims%j_start:field_dims%j_end),    &
               diff_coeff_2D(field_dims_offi_s%i_start-1:                & 
                                             field_dims_offi_s%i_end,    &
                             field_dims%j_start-1:field_dims%j_end+1),   & 
               recip_del_phi(field_dims_s%j_start:field_dims_s%j_end),   &
               recip_csphi(field_dims_s%j_start:field_dims_s%j_end) 

      INTEGER :: swp_bnds_fld_type

      INTEGER :: js, je, is, ie, ke, is2, ie2, is3, ie3
              


! ! External functions
! INTERFACE
!   ! For tri-diagonal solver
!   FUNCTION vert_trisolve(a, b, c, d, KK, tri_idx) 
!     IMPLICIT NONE
! INTEGER, INTENT(IN)     :: KK, tri_idx
! REAL                    :: vert_trisolve(tri_idx:KK) 
! REAL, INTENT(IN)        :: b(tri_idx:KK), c(tri_idx:KK-1), d(tri_idx:KK)
! REAL, INTENT(IN)        :: a(tri_idx+1:KK) 
!   END FUNCTION vert_trisolve
! END INTERFACE

! ----------------------------------------------------------------------
! Section 1.   Set metric terms and calculate D(r)/D(eta)
! ----------------------------------------------------------------------      

      IF (lhook) CALL dr_hook('POFIL_ZLEVEL',zhook_in,zhook_handle)

      L_init_matrix = .TRUE.

      swp_bnds_fld_type = fld_type
      IF ( fld_type == fld_type_w ) swp_bnds_fld_type = fld_type_p

      ND_shift = 0

      alpha = 2.0

      beta  = (1.0 - alpha)/alpha
      
      base_k = 1 - field_dims%k_start
      base_i = 1 - field_dims%i_start
      base_j = 1 - field_dims%j_start


! New dynamics grid shifting indices
      IF (fld_type == fld_type_u) base_i = 1
      IF (fld_type == fld_type_v) base_j = 1
      IF (fld_type == fld_type_w) base_k = 1


      js = field_dims_s%j_start
      je = field_dims_s%j_end
      is = field_dims_s%i_start
      ie = field_dims_s%i_end   
      ke = field_dims_s%k_end

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
      DO k = 1,metric_level
!   No metric terms needed above first constant rho-level
!   since they cancel on constant r-surfaces
        DO j = js,je
          DO i = is,ie
            recip_r_squared_delz(i,j,k) = 1.0 / ( r_levels(i,j,k) *      &
                                                  r_levels(i,j,k) *      &
                                  (r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1)) )   
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO


! ----------------------------------------------------------------------
! Section 1.1   Set up averaging arrays 
! ----------------------------------------------------------------------
! Vertical Weights
! Initialise bottom and top levels
      DO j = js,je
        DO i = is,ie
          wt_z_on2off(i,j,0) = 0.5
        END DO
      END DO
      DO j = js,je
        DO i = is,ie
          wt_z_on2off(i,j,ke) = 0.5
        END DO
      END DO

! ND does averaging in r
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
      DO k = 1,ke-r_dim_k_start
        DO j = js,je
          DO i =is,ie
            ! Vertical weight from field level to k+/- 1/2      
            wt_z_on2off(i,j,k) = (r_levels_offk(i,j,k) -                 &
                                        r_levels(i,j,k-base_k))          &
                                /(r_levels(i,j,k+1-base_k) -             &
                                        r_levels(i,j,k-base_k))       
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

! horizontal weights - Currently hard wired to 1/2
      DO j = field_dims_s%j_start,field_dims_s%j_end
        DO i = field_dims_offi_s%i_start,field_dims_offi_s%i_end
          ! Lambda weight from  field point to i+/- 1/2 
!           wt_x_on2off(i,j) = (lambda_offi(i)     - lambda(i-base_i))     &
!                             /(lambda(i+1-base_i) - lambda(i-base_i))
          wt_x_on2off(i,j) = 0.5
        END DO
      END DO 

      DO j = field_dims_offj_s%j_start,field_dims_offj_s%j_end
        DO i = field_dims_s%i_start,field_dims_s%i_end
          ! Phi weight from  field point to j+/- 1/2 
!           wt_y_on2off(i,j) = (phi_offj(j)     - phi(j-base_j))           &
!                             /(phi(j+1-base_j) - phi(j-base_j)) 
          wt_y_on2off(i,j) = 0.5
        END DO
      END DO


! ------------------------------------------------------------------------
! Section 2.0  Filtering and diffusion
! ------------------------------------------------------------------------
      DO i= 0, max_filter_rows !global_filter
        idx_begin(i) = sweep_begin(i)
        idx_end(i) = sweep_end(i)
      END DO

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
        ELSEIF ( j_sweep_end < idx_end(0) ) THEN
          j_sweep_end = idx_end(0)
        ELSEIF ( j_sweep_begin > idx_begin(0) ) THEN
          j_sweep_begin = idx_begin(0)
        END IF !  j_sweep_begin < 0
      END IF ! L_diff

      DO i_filter = 1, num_pass
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
! ------------------------------------------------------------------------
! Section 4.1   EW filtering and diffusion
! ------------------------------------------------------------------------
          SELECT CASE ( fld_type ) 
! ------------------------------------------------------------------------
! --- U fields ---
          CASE ( fld_type_u )
            ND_shift = 1
            is = field_dims%i_start
            ie = field_dims%i_end
! === Compute Matrix a ===================================================
            IF ( L_init_matrix ) THEN   
              DO i = is,ie
                recip_del_lambda(i) = 1.0/( lambda_offi(i+1)             &
                                          - lambda_offi(i)  )
              END DO

              k = 1
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie            
                  diff_coeff_av_2D(i,j) = 0.5*(diff_coeff(i-ND_shift+1,j)&
                                             + diff_coeff(i-ND_shift,j)) &
                                             *delta_lambda**2            &
                                             *recip_del_lambda(i)**2
                  del_r_av_z = 0.5*((r_levels_offi(i+1,j,k+1)            &
                                               -r_levels_offi(i,j,k+1))  &
                                  + (r_levels_offi(i+1,j,k)              &
                                                -r_levels_offi(i,j,k)))

                  ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)      &
                                 *del_r_av_z**2                          &
                                 /(r_levels(i,j,k+1)-r_levels(i,j,k))    &
                                 *diff_coeff_av_2D(i,j)                  &
                                 *r_levels_offk(i,j,k)**2
                END DO
              END DO

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(del_r_av_z,  &
!$OMP& i, j, k)
              DO k = 2, metric_level-1
                DO j = j_sweep_begin, j_sweep_end
                  DO i = is,ie
                    del_r_av_z = 0.5*((r_levels_offi(i+1,j,k+1)          &
                                                 -r_levels_offi(i,j,k+1))&
                                    + (r_levels_offi(i+1,j,k)            &
                                                  -r_levels_offi(i,j,k)))

                    ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)    &
                                   *del_r_av_z**2                        &
                                   /(r_levels(i,j,k+1)-r_levels(i,j,k))  &
                                   *diff_coeff_av_2D(i,j)                &
                                   *r_levels_offk(i,j,k)**2

                    del_r_av_z = 0.5*((r_levels_offi(i+1,j,k)            &
                                                 -r_levels_offi(i,j,k))  &
                                    + (r_levels_offi(i+1,j,k-1)          &
                                                 -r_levels_offi(i,j,k-1)))

                    am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)    &
                                   *del_r_av_z**2                        &
                                   /(r_levels(i,j,k)-r_levels(i,j,k-1))  &
                                   *diff_coeff_av_2D(i,j)                &
                                   *r_levels_offk(i,j,k-1)**2
                  END DO 
                END DO
              END DO
!$OMP END PARALLEL DO
              k = metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie 

                  del_r_av_z = 0.5*((r_levels_offi(i+1,j,k)              &
                                               -r_levels_offi(i,j,k))    &
                                  + (r_levels_offi(i+1,j,k-1)            &
                                               -r_levels_offi(i,j,k-1)))

                  am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)      &
                                 *del_r_av_z**2                          &
                                 /(r_levels(i,j,k)-r_levels(i,j,k-1))    &
                                 *diff_coeff_av_2D(i,j)                  &
                                 *r_levels_offk(i,j,k-1)**2
                END DO
              END DO

! Inversion of tridiagonal matrix
              k = 1
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/(1.0 - ap_3D(i,j,k))
                  c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                END DO
              END DO
              DO k = 2, metric_level-1
                DO j =  j_sweep_begin, j_sweep_end
                  DO i = is,ie
                    denom(i,j,k) = 1.0/(                                 &
                                  (1.0 - ap_3D(i,j,k) - am_3D(i,j,k))    &
                                       - c_new(i,j,k-1)*am_3D(i,j,k))
                    c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                  END DO
                END DO
              END DO
              k = metric_level
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/((1.0 - am_3D(i,j,k))               &
                                           - c_new(i,j,k-1)*am_3D(i,j,k))
                END DO
              END DO
              L_init_matrix = .FALSE.
            ENDIF !L_init_matrix
! ========================================================================
! === Compute G_lambda ===================================================
            k = 1
            DO j = j_sweep_begin, j_sweep_end
              is3 = field_dims_offi%i_start 
              ie3 = field_dims_offi_s%i_end
              DO i = is3,ie3

! This is G_lambda stored on p points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                diff_coeff_2D(i-ND_shift,j) = diff_coeff(i-ND_shift,j)   &
                                             *delta_lambda**2            &
                                             /( lambda(i) - lambda(i-1)) 

                f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                 &
                                    wt_x_on2off(i,j) *field(i-1,j,k+1)   &
                             + (1.0-wt_x_on2off(i,j))*field(i,j,k+1))    &
                             + (1.0-wt_z_on2off(i,j,k))*(                &
                                    wt_x_on2off(i,j) *field(i-1,j,k)     &
                             + (1.0-wt_x_on2off(i,j))*field(i,j,k) )

                f_av_xz_km =  wt_x_on2off(i,j) *field(i-1,j,k) +         &
                         (1.0-wt_x_on2off(i,j))*field(i,j,k) 

                G_lambda(i,j,k) = ( (field(i,j,k) - field(i-1,j,k))      &
                                *(r_levels_offk(i,j,k+base_k)            &
                                - r_levels_offk(i,j,k+base_k-1))         &
                                - (f_av_xz_k_2D(i,j) - f_av_xz_km)       &
                                *(r_levels(i,j,k) - r_levels(i-1,j,k)) ) &
                                *diff_coeff_2D(i-ND_shift,j)             &
                                *r_levels_offi(i,j,k)**2
              END DO
            END DO

            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                is3 = field_dims_offi%i_start 
                ie3 = field_dims_offi_s%i_end
                DO i = is3,ie3
! This is G_lambda stored on p points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_xz_km = f_av_xz_k_2D(i,j)

                  f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k)*(               &
                                      wt_x_on2off(i,j) *field(i-1,j,k+1) &
                               + (1.0-wt_x_on2off(i,j))*field(i,j,k+1))  &
                               + (1.0-wt_z_on2off(i,j,k))*(              &
                                      wt_x_on2off(i,j) *field(i-1,j,k)   &
                               + (1.0-wt_x_on2off(i,j))*field(i,j,k) )

                  G_lambda(i,j,k) = ( (field(i,j,k) - field(i-1,j,k))    &
                                  *(r_levels_offk(i,j,k+base_k)          &
                                  - r_levels_offk(i,j,k+base_k-1))       &
                                  - (f_av_xz_k_2D(i,j) - f_av_xz_km)     &
                                  *(r_levels(i,j,k) - r_levels(i-1,j,k)))&
                                  *diff_coeff_2D(i-ND_shift,j)           &
                                  *r_levels_offi(i,j,k)**2
                END DO
              END DO
            END DO

            k = metric_level
            IF ( metric_level == field_dims%k_end ) THEN
              DO j = j_sweep_begin, j_sweep_end

                is3 = field_dims_offi%i_start 
                ie3 = field_dims_offi_s%i_end
                DO i = is3,ie3

! This is G_lambda stored on p points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_xz_km = f_av_xz_k_2D(i,j)

                  f_av_xz_k_2D(i,j) = wt_x_on2off(i,j) *field(i-1,j,k)   &
                               + (1.0-wt_x_on2off(i,j))*field(i,j,k) 

                  G_lambda(i,j,k) = ( (field(i,j,k) - field(i-1,j,k))    &
                                  *(r_levels_offk(i,j,k+base_k)          &
                                  - r_levels_offk(i,j,k+base_k-1))       &
                                  - (f_av_xz_k_2D(i,j) - f_av_xz_km)     &
                                  *(r_levels(i,j,k) - r_levels(i-1,j,k)))&
                                  *diff_coeff_2D(i-ND_shift,j)           &
                                  *r_levels_offi(i,j,k)**2
                END DO
              END DO
            ELSE 
              DO j = j_sweep_begin, j_sweep_end
                is3 = field_dims_offi%i_start 
                ie3 = field_dims_offi_s%i_end
                DO i = is3,ie3

! This is G_lambda stored on p points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_xz_km = f_av_xz_k_2D(i,j)

                  f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k)*(               &
                                      wt_x_on2off(i,j) *field(i-1,j,k+1) &
                               + (1.0-wt_x_on2off(i,j))*field(i,j,k+1))  &
                               + (1.0-wt_z_on2off(i,j,k))*(              &
                                      wt_x_on2off(i,j) *field(i-1,j,k)   &
                               + (1.0-wt_x_on2off(i,j))*field(i,j,k) )

                  G_lambda(i,j,k) = ( (field(i,j,k) - field(i-1,j,k))    &
                                  *(r_levels_offk(i,j,k+base_k)          &
                                  - r_levels_offk(i,j,k+base_k-1))       &
                                  - (f_av_xz_k_2D(i,j) - f_av_xz_km)     &
                                  *(r_levels(i,j,k) - r_levels(i-1,j,k)))&
                                  *diff_coeff_2D(i-ND_shift,j)           &
                                  *r_levels_offi(i,j,k)**2
                END DO
              END DO
            END IF
! end of G_lambda
! ========================================================================
! === Compute RHS b ======================================================
            k = 1
            DO j = j_sweep_begin, j_sweep_end
              i = field_dims%i_start - 1
              f_av_xz_i = wt_z_on2off(i,j,k)*(                           &
                          wt_x_on2off(i,j)*field(i,j,k+1)                &
                   + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))            &
                   + (1.0-wt_z_on2off(i,j,k))*(                          &
                          wt_x_on2off(i,j)*field(i,j,k)                  &
                   + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

              DO i = field_dims%i_start,field_dims%i_end               
                ! This is G_lambda_eta on (i,j-1/2,k) points  
                f_av_xz_im = f_av_xz_i
                  
                f_av_xz_i = wt_z_on2off(i,j,k)*(                         &
                            wt_x_on2off(i,j)*field(i,j,k+1)              &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))          &
                     + (1.0-wt_z_on2off(i,j,k))*(                        &
                            wt_x_on2off(i,j)*field(i,j,k)                &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                G_lambda_eta_km = 0.0

                G_lambda_eta_k_2D(i,j) = ( f_av_xz_i - f_av_xz_im )      &
                                        *0.5*((r_levels_offi(i+1,j,k+1)  &
                                                -r_levels_offi(i,j,k+1)) &
                                            + (r_levels_offi(i+1,j,k)    &
                                                -r_levels_offi(i,j,k)))  &
                                       *diff_coeff_av_2D(i,j)            &
                                       *r_levels_offk(i,j,k)**2

                b_3D(i,j,k) = field(i,j,k) +                             &
                              ((G_lambda(i+1,j,k) - G_lambda(i,j,k))     &
                             *recip_del_lambda(i)                        &
                            - (G_lambda_eta_k_2D(i,j) - G_lambda_eta_km))&
                             *recip_r_squared_delz(i,j,k)                &
                            - beta*ap_3D(i,j,k)                          &
                             *(field(i,j,k+1) - field(i,j,k))
              END DO
            END DO
            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                i = field_dims%i_start - 1
                f_av_xz_i = wt_z_on2off(i,j,k)*(                         &
                            wt_x_on2off(i,j)*field(i,j,k+1)              &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))          &
                     + (1.0-wt_z_on2off(i,j,k))*(                        &
                            wt_x_on2off(i,j)*field(i,j,k)                &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                DO i = field_dims%i_start,field_dims%i_end
                  ! This is G_lambda_eta on (i,j-1/2,k) points  
                  f_av_xz_im = f_av_xz_i
                  
                  f_av_xz_i = wt_z_on2off(i,j,k)*(                       &
                              wt_x_on2off(i,j)*field(i,j,k+1)            &
                       + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))        &
                       + (1.0-wt_z_on2off(i,j,k))*(                      &
                              wt_x_on2off(i,j)*field(i,j,k)              &
                       + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) ) 

                  G_lambda_eta_km = G_lambda_eta_k_2D(i,j)

                  G_lambda_eta_k_2D(i,j) = (f_av_xz_i - f_av_xz_im)      &
                                          *0.5*((r_levels_offi(i+1,j,k+1)&
                                                 -r_levels_offi(i,j,k+1))&
                                              + (r_levels_offi(i+1,j,k)  &
                                                 -r_levels_offi(i,j,k))) &
                                          *diff_coeff_av_2D(i,j)         &
                                          *r_levels_offk(i,j,k)**2       

                  b_3D(i,j,k) = field(i,j,k) +                           &
                                ((G_lambda(i+1,j,k) - G_lambda(i,j,k))   &
                               *recip_del_lambda(i)                      &
                              - (G_lambda_eta_k_2D(i,j)                  &
                                                    - G_lambda_eta_km))  &
                               *recip_r_squared_delz(i,j,k)              &
                              - beta*ap_3D(i,j,k)                        &
                               *(field(i,j,k+1) - field(i,j,k))          &
                              + beta*am_3D(i,j,k)                        &
                               *(field(i,j,k)   - field(i,j,k-1))
                END DO
              END DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end              
                G_lambda_eta_km = G_lambda_eta_k_2D(i,j)

                G_lambda_eta_k_2D(i,j) = 0.0

                b_3D(i,j,k) = field(i,j,k) +                             &
                              ((G_lambda(i+1,j,k) - G_lambda(i,j,k))     &
                             *recip_del_lambda(i)                        &
                            - (G_lambda_eta_k_2D(i,j) - G_lambda_eta_km))&
                             *recip_r_squared_delz(i,j,k)                &
                            + beta*am_3D(i,j,k)                          &
                             *(field(i,j,k) - field(i,j,k-1))

              END DO
            END DO
! ========================================================================
! === trisolve ===========================================================
! Inlining of trisolve
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,1) = b_3D(i,j,1)*denom(i,j,1)
              END DO
            END DO

            DO k = 2,metric_level-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  field(i,j,k) = (b_3D(i,j,k)-field(i,j,k-1)             &
                                 *am_3D(i,j,k))*denom(i,j,k)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = (b_3D(i,j,k) - field(i,j,k-1)             &
                               *am_3D(i,j,k))*denom(i,j,k)
              END DO
            END DO
! Backward solve sweep
           DO k = metric_level-1,1,-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
             DO j = j_sweep_begin, j_sweep_end
               DO i = field_dims%i_start,field_dims%i_end
                  field(i,j,k) = field(i,j,k)-c_new(i,j,k)*field(i,j,k+1)
               END DO
             END DO
!$OMP END PARALLEL DO
           END DO

! Start of upper level
           is = field_dims_offi%i_start
           ie = field_dims_offi_s%i_end
           is2 = field_dims%i_start
           ie2 = field_dims%i_end

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
           DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
             DO j = j_sweep_begin, j_sweep_end
               DO i = is,ie
                 G_lambda(i,j,k) = ( field(i,j,k) - field(i-1,j,k) )     &
                                  /( lambda(i)    - lambda(i-1)    )     &
                                  *delta_lambda**2                       &
                                  *diff_coeff(i-ND_shift,j)
               END DO
               DO i = is2,ie2
                 field(i,j,k) = field(i,j,k)                             &
                             + (G_lambda(i+1,j,k) - G_lambda(i,j,k))     &
                              * recip_del_lambda(i)

               END DO
             END DO
           END DO
!$OMP END PARALLEL DO
! ========================================================================
! --- V field ---
         CASE ( fld_type_v ) 
! === Compute Matrix a ===================================================
           is = field_dims%i_start
           ie = field_dims%i_end

           IF ( L_init_matrix ) THEN
             DO i = is,ie
               recip_del_lambda(i) = 1.0/                                &
                                      (lambda_offi(i) - lambda_offi(i-1))
             END DO
             k = 1
             DO j = j_sweep_begin, j_sweep_end
               DO i = is,ie
                 diff_coeff_av_2D(i,j) = 0.5*(diff_coeff(i,j)            &
                                             +diff_coeff(i-1,j))         &
                                         *delta_lambda**2                &
                                         *recip_del_lambda(i)**2

                 del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)               &
                                            -r_levels_offi(i-1,j,k+1))   &
                                 + (r_levels_offi(i,j,k)                 &
                                            -r_levels_offi(i-1,j,k))) 

                 ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)       &
                                *del_r_av_z**2                           &
                                /(r_levels(i,j,k+1)-r_levels(i,j,k))     &
                                *diff_coeff_av_2D(i,j)                   &
                                *r_levels_offk(i,j,k)**2

               END DO
             END DO
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(del_r_av_z,  &
!$OMP& i, j, k)
             DO k = 2, metric_level-1
               DO j = j_sweep_begin, j_sweep_end
                 DO i = is,ie
                   del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)             &
                                             -r_levels_offi(i-1,j,k+1))  &
                                  + (r_levels_offi(i,j,k)                &
                                             -r_levels_offi(i-1,j,k))) 

                   ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)     &
                                  *del_r_av_z**2                         &
                                  /(r_levels(i,j,k+1)-r_levels(i,j,k))   &
                                  *diff_coeff_av_2D(i,j)                 &
                                  *r_levels_offk(i,j,k)**2

                   del_r_av_z = 0.5*((r_levels_offi(i,j,k)               &
                                             -r_levels_offi(i-1,j,k))    &
                                  + (r_levels_offi(i,j,k-1)              &
                                             -r_levels_offi(i-1,j,k-1))) 

                   am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)     &
                                  *del_r_av_z**2                         &
                                  /(r_levels(i,j,k)-r_levels(i,j,k-1))   &
                                  *diff_coeff_av_2D(i,j)                 &
                                  *r_levels_offk(i,j,k-1)**2
                 END DO
               END DO
             END DO
!$OMP END PARALLEL DO
             k = metric_level
             DO j = j_sweep_begin, j_sweep_end
               DO i = is,ie

                 del_r_av_z = 0.5*((r_levels_offi(i,j,k)                 &
                                           -r_levels_offi(i-1,j,k))      &
                                + (r_levels_offi(i,j,k-1)                &
                                           -r_levels_offi(i-1,j,k-1))) 

                 am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)       &
                                *del_r_av_z**2                           &
                                /(r_levels(i,j,k)-r_levels(i,j,k-1))     &
                                *diff_coeff_av_2D(i,j)                   &
                                *r_levels_offk(i,j,k-1)**2
               END DO
             END DO
! Inversion of tridiagonal matrix
             k = 1
             DO j =  j_sweep_begin, j_sweep_end
               DO i = is,ie
                 denom(i,j,k) = 1.0/(1.0 - ap_3D(i,j,k) )
                 c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
               END DO
             END DO
             DO k = 2, metric_level-1
               DO j =  j_sweep_begin, j_sweep_end
                 DO i = is,ie
                   denom(i,j,k) = 1.0/(                                  &
                                    (1.0 - ap_3D(i,j,k) - am_3D(i,j,k))  &
                                         - c_new(i,j,k-1)*am_3D(i,j,k) )
                   c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                 END DO
               END DO
             END DO
             k = metric_level
             DO j =  j_sweep_begin, j_sweep_end
               DO i = is,ie
                   denom(i,j,k) = 1.0/( (1.0 - am_3D(i,j,k))             &
                                - c_new(i,j,k-1)*am_3D(i,j,k))
               END DO
             END DO
             L_init_matrix = .FALSE.
           END IF !L_init_matrix
! ========================================================================
! === G_lambda ===========================================================
           k = 1
           DO j = j_sweep_begin, j_sweep_end
             is3 = field_dims_offi_s%i_start 
             ie3 = field_dims_offi%i_end
             DO i = is3,ie3
! This is G_lambda stored on points (i,j,k-1/2)
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
               diff_coeff_2D(i,j) = diff_coeff(i,j)*delta_lambda**2      &
                                    /( lambda(i+1) - lambda(i)) 

               f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                  &
                                   wt_x_on2off(i,j) *field(i,j,k+1)      &
                            + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))   &
                            + (1.0-wt_z_on2off(i,j,k))*(                 &
                                   wt_x_on2off(i,j) *field(i,j,k)        &
                            + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

               f_av_xz_km = wt_x_on2off(i,j) *field(i,j,k) +             &
                       (1.0-wt_x_on2off(i,j))*field(i+1,j,k)

               G_lambda(i,j,k) = ( (field(i+1,j,k) - field(i,j,k))       &
                                  *(r_levels_offk(i,j,k+base_k) -        &
                                    r_levels_offk(i,j,k+base_k-1))       &
                                 - (f_av_xz_k_2D(i,j) - f_av_xz_km)      &
                                  *(r_levels(i+1,j,k) - r_levels(i,j,k)))&
                                  *diff_coeff_2D(i,j)                    &
                                  *r_levels_offi(i,j,k)**2
             END DO
           END DO
           DO k = 2, metric_level-1
             DO j = j_sweep_begin, j_sweep_end
               is3 = field_dims_offi_s%i_start 
               ie3 = field_dims_offi%i_end
               DO i = is3,ie3
! This is G_lambda stored on points (i,j,k-1/2)
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                 f_av_xz_km = f_av_xz_k_2D(i,j)

                 f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                &
                                     wt_x_on2off(i,j) *field(i,j,k+1)    &
                              + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1)) &
                              + (1.0-wt_z_on2off(i,j,k))*(               &
                                     wt_x_on2off(i,j) *field(i,j,k)      &
                              + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                 G_lambda(i,j,k) = ( (field(i+1,j,k) - field(i,j,k))     &
                                    *(r_levels_offk(i,j,k+base_k) -      &
                                      r_levels_offk(i,j,k+base_k-1))     &
                                   - (f_av_xz_k_2D(i,j) - f_av_xz_km)    &
                                    *(r_levels(i+1,j,k)-r_levels(i,j,k)))&
                                    *diff_coeff_2D(i,j)                  &
                                    *r_levels_offi(i,j,k)**2
               END DO
             END DO
           END DO
           k = metric_level
           IF ( metric_level == field_dims%k_end ) THEN
             DO j = j_sweep_begin, j_sweep_end
               is3 = field_dims_offi_s%i_start 
               ie3 = field_dims_offi%i_end
               DO i = is3,ie3
! This is G_lambda stored on points (i,j,k-1/2)
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                 f_av_xz_km = f_av_xz_k_2D(i,j)

                 f_av_xz_k_2D(i,j) = wt_x_on2off(i,j) *field(i,j,k) +    &
                                (1.0-wt_x_on2off(i,j))*field(i+1,j,k) 

                 G_lambda(i,j,k) = ( (field(i+1,j,k) - field(i,j,k))     &
                                    *(r_levels_offk(i,j,k+base_k) -      &
                                      r_levels_offk(i,j,k+base_k-1))     &
                                   - (f_av_xz_k_2D(i,j) - f_av_xz_km)    &
                                    *(r_levels(i+1,j,k)-r_levels(i,j,k)))&
                                    *diff_coeff_2D(i,j)                  &
                                    *r_levels_offi(i,j,k)**2
               END DO
             END DO
           ELSE
             DO j = j_sweep_begin, j_sweep_end
               is3 = field_dims_offi_s%i_start 
               ie3 = field_dims_offi%i_end
               DO i = is3,ie3
! This is G_lambda stored on points (i,j,k-1/2)
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                 f_av_xz_km = f_av_xz_k_2D(i,j)

                 f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                &
                                     wt_x_on2off(i,j) *field(i,j,k+1)    &
                              + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1)) &
                              + (1.0-wt_z_on2off(i,j,k))*(               &
                                     wt_x_on2off(i,j) *field(i,j,k)      &
                              + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                 G_lambda(i,j,k) = ( (field(i+1,j,k) - field(i,j,k))     &
                                    *(r_levels_offk(i,j,k+base_k) -      &
                                      r_levels_offk(i,j,k+base_k-1))     &
                                   - (f_av_xz_k_2D(i,j) - f_av_xz_km)    &
                                    *(r_levels(i+1,j,k)-r_levels(i,j,k)))&
                                    *diff_coeff_2D(i,j)                  &
                                    *r_levels_offi(i,j,k)**2
               END DO
             END DO
           ENDIF
! ========================================================================
! === Compute RHS b ======================================================
           k = 1
           DO j = j_sweep_begin, j_sweep_end
             i = field_dims%i_start - 1
             f_av_xz_i = wt_z_on2off(i,j,k)*(                            &
                         wt_x_on2off(i,j)*field(i,j,k+1-base_k)          &
                  + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k))      &
                  + (1.0-wt_z_on2off(i,j,k))*(                           &
                         wt_x_on2off(i,j)*field(i,j,k-base_k)            &
                  + (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )

             DO i = field_dims%i_start,field_dims%i_end
               ! This is G_lambda_eta on (i-1/2,j,k) points
               f_av_xz_im = f_av_xz_i

               f_av_xz_i = wt_z_on2off(i,j,k)*(                          &
                           wt_x_on2off(i,j)*field(i,j,k+1-base_k)        &
                    + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k))    &
                    + (1.0-wt_z_on2off(i,j,k))*(                         &
                           wt_x_on2off(i,j)*field(i,j,k-base_k)          &
                    + (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )
               
               G_lambda_eta_km = 0.0

               G_lambda_eta_k_2D(i,j) = (( f_av_xz_i - f_av_xz_im ))     &
                                        *0.5*((r_levels_offi(i,j,k+1)    &
                                              -r_levels_offi(i-1,j,k+1)) &
                                            + (r_levels_offi(i,j,k)      &
                                              -r_levels_offi(i-1,j,k)))  &
                                        *diff_coeff_av_2D(i,j)           &
                                        *r_levels_offk(i,j,k)**2

               b_3D(i,j,k) = field(i,j,k) +                              &
                            ( (G_lambda(i,j,k) - G_lambda(i-1,j,k))      &
                             *recip_del_lambda(i)                        &
                            - (G_lambda_eta_k_2D(i,j) - G_lambda_eta_km))&
                             * recip_r_squared_delz(i,j,k)               &
                            - beta*ap_3D(i,j,k)                          &
                             *(field(i,j,k+1) - field(i,j,k)) 
             END DO
           END DO

           DO k = 2, metric_level-1
             DO j = j_sweep_begin, j_sweep_end
               i = field_dims%i_start - 1

               f_av_xz_i = wt_z_on2off(i,j,k)*(                          &
                             wt_x_on2off(i,j)*field(i,j,k+1-base_k)      &
                      + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k))  &
                    + (1.0-wt_z_on2off(i,j,k))*(                         &
                             wt_x_on2off(i,j)*field(i,j,k-base_k)        &
                      + (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )

               DO i = field_dims%i_start,field_dims%i_end
                 ! This is G_lambda_eta on (i-1/2,j,k) points
                 f_av_xz_im = f_av_xz_i

                 f_av_xz_i = wt_z_on2off(i,j,k)*(                        &
                              wt_x_on2off(i,j)*field(i,j,k+1-base_k)     &
                       + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k)) &
                      + (1.0-wt_z_on2off(i,j,k))*(                       &
                              wt_x_on2off(i,j)*field(i,j,k-base_k)       &
                       + (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )
               
                 G_lambda_eta_km = G_lambda_eta_k_2D(i,j)
 
                 G_lambda_eta_k_2D(i,j) = (( f_av_xz_i - f_av_xz_im ))   &
                                          *0.5*((r_levels_offi(i,j,k+1)  &
                                             -r_levels_offi(i-1,j,k+1))  &
                                          + (r_levels_offi(i,j,k)        &
                                            -r_levels_offi(i-1,j,k)))    &
                                          *diff_coeff_av_2D(i,j)         &
                                          *r_levels_offk(i,j,k)**2

                 b_3D(i,j,k) = field(i,j,k) +                            &
                              ( (G_lambda(i,j,k) - G_lambda(i-1,j,k))    &
                               *recip_del_lambda(i)                      &
                              - (G_lambda_eta_k_2D(i,j)                  &
                                                   - G_lambda_eta_km))   &
                               * recip_r_squared_delz(i,j,k)             &
                              - beta*ap_3D(i,j,k)                        &
                               *(field(i,j,k+1) - field(i,j,k))          &
                              + beta*am_3D(i,j,k)                        &
                               *(field(i,j,k)   - field(i,j,k-1))
               END DO
             END DO
           END DO
           k = metric_level
             DO j = j_sweep_begin, j_sweep_end
               DO i = field_dims%i_start,field_dims%i_end

                 G_lambda_eta_km = G_lambda_eta_k_2D(i,j)
 
                 G_lambda_eta_k_2D(i,j) = 0.0

                 b_3D(i,j,k) = field(i,j,k) +                            &
                            ( (G_lambda(i,j,k)   - G_lambda(i-1,j,k))    &
                              *recip_del_lambda(i)                       &
                            - (G_lambda_eta_k_2D(i,j) - G_lambda_eta_km))&
                              *recip_r_squared_delz(i,j,k)               &
                             + beta*am_3D(i,j,k)                         &
                              *(field(i,j,k) - field(i,j,k-1)) 
               END DO
             END DO
! ========================================================================
! === trisolve ===========================================================
! Inlining of trisolve
           DO j = j_sweep_begin, j_sweep_end
             DO i = field_dims%i_start,field_dims%i_end
               field(i,j,1) = b_3D(i,j,1)*denom(i,j,1)
             END DO
           END DO

           DO k = 2,metric_level-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
             DO j = j_sweep_begin, j_sweep_end
               DO i = field_dims%i_start,field_dims%i_end
                 field(i,j,k) = (b_3D(i,j,k)-field(i,j,k-1)              &
                                *am_3D(i,j,k))*denom(i,j,k)
               END DO
             END DO
!$OMP END PARALLEL DO
           END DO
           k = metric_level
             DO j = j_sweep_begin, j_sweep_end
               DO i = field_dims%i_start,field_dims%i_end
                 field(i,j,k) = (b_3D(i,j,k) - field(i,j,k-1)            &
                                *am_3D(i,j,k))*denom(i,j,k)
               END DO
             END DO
! Backward solve sweep
           DO k = metric_level-1,1,-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
             DO j = j_sweep_begin, j_sweep_end
               DO i = field_dims%i_start,field_dims%i_end
                 field(i,j,k) = field(i,j,k)-c_new(i,j,k)*field(i,j,k+1)
               END DO
             END DO
!$OMP END PARALLEL DO
           END DO
! ========================================================================
! === Upper level ========================================================
           is = field_dims_offi_s%i_start
           ie = field_dims_offi%i_end 
           is2 = field_dims%i_start
           ie2 = field_dims%i_end
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
           DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
             DO j = j_sweep_begin, j_sweep_end
               DO i = is,ie 
                 G_lambda(i,j,k) = ( field(i+1,j,k) - field(i,j,k) )     &
                                  /( lambda(i+1)    - lambda(i)    )     &
                                  * delta_lambda**2*diff_coeff(i,j)
               END DO

               DO i = is2,ie2
                 field(i,j,k) = field(i,j,k)                             &
                             + (G_lambda(i,j,k) - G_lambda(i-1,j,k))     &
                              *recip_del_lambda(i)
               END DO
             END DO
           END DO
!$OMP END PARALLEL DO
! ------------------------------------------------------------------------
! --- P fields ---
         CASE ( fld_type_p )   
           WRITE(6,*) 'Diffusion of fields on p points not supported'
! === Compute Matrix a ===================================================
           is = field_dims%i_start
           ie = field_dims%i_end

           IF ( L_init_matrix ) THEN
             DO i = is,ie
               recip_del_lambda(i) = 1.0/                                &
                                 (lambda_offi(i) - lambda_offi(i-1))
             END DO

             k = 1
             DO j = j_sweep_begin, j_sweep_end
               DO i = is,ie
                 diff_coeff_av_2D(i,j) = 0.5*(diff_coeff(i,j)            &
                                                     +diff_coeff(i-1,j)) &
                                         *delta_lambda**2                &
                                         *recip_del_lambda(i)**2

                 del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)               &
                                           -r_levels_offi(i-1,j,k+1))    &
                                 + (r_levels_offi(i,j,k)                 &
                                           -r_levels_offi(i-1,j,k))) 

                 ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)       &
                                *del_r_av_z**2                           &
                                /(r_levels(i,j,k+1)-r_levels(i,j,k))     &
                                *diff_coeff_av_2D(i,j)                   &
                                *r_levels_offk(i,j,k)**2

               END DO
             END DO
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(del_r_av_z, &
!$OMP& i, j, k)
             DO k = 2, metric_level-1
               DO j = j_sweep_begin, j_sweep_end
                 DO i = is,ie
                   del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)             &
                                             -r_levels_offi(i-1,j,k+1))  &
                                   + (r_levels_offi(i,j,k)               &
                                             -r_levels_offi(i-1,j,k))) 

                   ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)     &
                                  *del_r_av_z**2                         &
                                  /(r_levels(i,j,k+1)-r_levels(i,j,k))   &
                                  *diff_coeff_av_2D(i,j)                 &
                                  *r_levels_offk(i,j,k)**2

                   del_r_av_z = 0.5*((r_levels_offi(i,j,k)               &
                                             -r_levels_offi(i-1,j,k))    &
                                   + (r_levels_offi(i,j,k-1)             &
                                             -r_levels_offi(i-1,j,k-1))) 

                   am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)     &
                                  *del_r_av_z**2                         &
                                  /(r_levels(i,j,k)-r_levels(i,j,k-1))   &
                                  *diff_coeff_av_2D(i,j)                 &
                                  *r_levels_offk(i,j,k-1)**2
                 END DO
               END DO
             END DO
!$OMP END PARALLEL DO
             k = metric_level
             DO j = j_sweep_begin, j_sweep_end
               DO i = is,ie

                 del_r_av_z = 0.5*((r_levels_offi(i,j,k)                 &
                                           -r_levels_offi(i-1,j,k))      &
                                 + (r_levels_offi(i,j,k-1)               &
                                           -r_levels_offi(i-1,j,k-1))) 

                 am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)       &
                                *del_r_av_z**2                           &
                                /(r_levels(i,j,k)-r_levels(i,j,k-1))     &
                                *diff_coeff_av_2D(i,j)                   &
                                *r_levels_offk(i,j,k-1)**2
              END DO
            END DO

! Inversion of tridiagonal matrix
            k = 1
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/(1.0 - ap_3D(i,j,k))
                  c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                END DO
              END DO
            DO k = 2, metric_level-1
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/((1.0 - ap_3D(i,j,k) - am_3D(i,j,k))&
                                 -c_new(i,j,k-1)*am_3D(i,j,k))
                  c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                END DO
              END DO
            END DO
            k = metric_level
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/((1.0 - am_3D(i,j,k))               &             
                                 -c_new(i,j,k-1)*am_3D(i,j,k))
                END DO
              END DO
              L_init_matrix = .FALSE.
            END IF
! ========================================================================
! === G_lambda ===========================================================
            k = 1
              DO j = j_sweep_begin, j_sweep_end
               is3 = field_dims_offi_s%i_start 
               ie3 = field_dims_offi%i_end
               DO i = is3,ie3
! This is G_lambda stored on u points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  diff_coeff_2D(i,j) = diff_coeff(i,j)*delta_lambda**2   &
                                       /( lambda(i+1) - lambda(i)) 

                  f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k)*(               &
                                      wt_x_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))&
                               + (1.0-wt_z_on2off(i,j,k))*(              &
                                      wt_x_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                  f_av_xz_km= wt_x_on2off(i,j) *field(i,j,k) +           &
                               (1.0-wt_x_on2off(i,j))*field(i+1,j,k)

                  G_lambda(i,j,k) = (( field(i+1,j,k) - field(i,j,k))    &
                                    *(r_levels_offk(i,j,k+base_k) -      &
                                      r_levels_offk(i,j,k+base_k-1))     &
                                   - (f_av_xz_k_2D(i,j) - f_av_xz_km)    &
                                    *(r_levels(i+1,j,k)                  &
                                                    - r_levels(i,j,k)))  &
                                    *diff_coeff_2D(i,j)                  &
                                    *r_levels_offi(i,j,k)**2
                END DO
              END DO
            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                is3 = field_dims_offi_s%i_start 
                ie3 = field_dims_offi%i_end
               DO i = is3,ie3
! This is G_lambda stored on u points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_xz_km = f_av_xz_k_2D(i,j)

                  f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k)*(               &
                                      wt_x_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))&
                               + (1.0-wt_z_on2off(i,j,k))*(              &
                                      wt_x_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                  G_lambda(i,j,k) = ((field(i+1,j,k) - field(i,j,k))     &
                                    *(r_levels_offk(i,j,k+base_k)        &
                                    - r_levels_offk(i,j,k+base_k-1))     &
                                    - (f_av_xz_k_2D(i,j) - f_av_xz_km )  &
                                    *(r_levels(i+1,j,k)                  &
                                                  - r_levels(i,j,k)))    &
                                    *diff_coeff_2D(i,j)                  &
                                    *r_levels_offi(i,j,k)**2
                END DO
              END DO
            END DO
            k = metric_level
            IF ( metric_level == field_dims%k_end ) THEN
              DO j = j_sweep_begin, j_sweep_end
                is3 = field_dims_offi_s%i_start 
                ie3 = field_dims_offi%i_end
                DO i = is3,ie3
! This is G_lambda stored on u points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_xz_km = f_av_xz_k_2D(i,j)

                  f_av_xz_k_2D(i,j) = wt_x_on2off(i,j) *field(i,j,k) +   &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) 

                  G_lambda(i,j,k) = ((field(i+1,j,k) - field(i,j,k))     &
                                    *(r_levels_offk(i,j,k+base_k)        &
                                    - r_levels_offk(i,j,k+base_k-1))     &
                                    - (f_av_xz_k_2D(i,j) - f_av_xz_km)   &
                                    *(r_levels(i+1,j,k)                  &
                                                    - r_levels(i,j,k)))  &
                                    *diff_coeff_2D(i,j)                  &
                                    *r_levels_offi(i,j,k)**2
                END DO
              END DO
            ELSE
              DO j = j_sweep_begin, j_sweep_end
                is3 = field_dims_offi_s%i_start 
                ie3 = field_dims_offi%i_end
                DO i = is3,ie3
! This is G_lambda stored on u points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_xz_km = f_av_xz_k_2D(i,j)

                  f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k)*(               &
                                      wt_x_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))&
                               + (1.0-wt_z_on2off(i,j,k))*(              &
                                      wt_x_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                  G_lambda(i,j,k) = ((field(i+1,j,k) - field(i,j,k))     &
                                    *(r_levels_offk(i,j,k+base_k)        &
                                    - r_levels_offk(i,j,k+base_k-1))     &
                                    - (f_av_xz_k_2D(i,j) - f_av_xz_km)   &
                                    *(r_levels(i+1,j,k)                  &
                                                  - r_levels(i,j,k)))    &
                                    *diff_coeff_2D(i,j)                  &
                                    *r_levels_offi(i,j,k)**2
                END DO
              END DO
            END IF
! ========================================================================
! === Compute RHS b ======================================================
            k = 1
            DO j = j_sweep_begin, j_sweep_end
              i = is - 1
              f_av_xz_i = wt_z_on2off(i,j,k)*(                           &
                          wt_x_on2off(i,j)*field(i,j,k+1-base_k)         &
                   + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k))     &
                   + (1.0-wt_z_on2off(i,j,k))*(                          &
                          wt_x_on2off(i,j)*field(i,j,k-base_k)           &
                  + (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )
              DO i = is,ie
                ! RHS term
                ! This is G_lambda_eta on w points
                f_av_xz_im = f_av_xz_i

                f_av_xz_i = wt_z_on2off(i,j,k)*(                         &
                            wt_x_on2off(i,j)*field(i,j,k+1-base_k)       &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k))   &
                     + (1.0-wt_z_on2off(i,j,k))*(                        &
                            wt_x_on2off(i,j)*field(i,j,k-base_k)         &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )

                G_lambda_eta_km = 0.0

                G_lambda_eta_k_2D(i,j) = (( f_av_xz_i - f_av_xz_im ))    &
                                         *0.5*((r_levels_offi(i,j,k+1)   &
                                              - r_levels_offi(i-1,j,k+1))&
                                       + (r_levels_offi(i,j,k)           &
                                              - r_levels_offi(i-1,j,k))) &
                                         *diff_coeff_av_2D(i,j)          &
                                         *r_levels_offk(i,j,k)**2

                b_3D(i,j,k) = field(i,j,k) +                             &
                            ((G_lambda(i,j,k) - G_lambda(i-1,j,k))       &
                             *recip_del_lambda(i)                        &
                            - (G_lambda_eta_k_2D(i,j) - G_lambda_eta_km))&
                             *recip_r_squared_delz(i,j,k)                &
                            - beta*ap_3D(i,j,k)                          &
                             *(field(i,j,k+1) - field(i,j,k))
              END DO
            END DO
            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                i = is - 1
                f_av_xz_i = wt_z_on2off(i,j,k)*(                         &
                            wt_x_on2off(i,j)*field(i,j,k+1-base_k)       &
                      + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k))  &
                      + (1.0-wt_z_on2off(i,j,k))*(                       &
                             wt_x_on2off(i,j)*field(i,j,k-base_k)        &
                      + (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )
                DO i = is,ie
                ! RHS term
                ! This is G_lambda_eta on w points
                  f_av_xz_im = f_av_xz_i

                  f_av_xz_i = wt_z_on2off(i,j,k)*(                       &
                              wt_x_on2off(i,j)*field(i,j,k+1-base_k)     &
                       + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1-base_k)) &
                       + (1.0-wt_z_on2off(i,j,k))*(                      &
                              wt_x_on2off(i,j)*field(i,j,k-base_k)       &
                       + (1.0-wt_x_on2off(i,j))*field(i+1,j,k-base_k) )

                  G_lambda_eta_km = G_lambda_eta_k_2D(i,j)

                  G_lambda_eta_k_2D(i,j) = (( f_av_xz_i - f_av_xz_im ))  &
                                          *0.5*((r_levels_offi(i,j,k+1)  &
                                              - r_levels_offi(i-1,j,k+1))&
                                          + (r_levels_offi(i,j,k)        &
                                               - r_levels_offi(i-1,j,k)))&
                                          *diff_coeff_av_2D(i,j)         &
                                          *r_levels_offk(i,j,k)**2

                  b_3D(i,j,k) = field(i,j,k) +                           &
                              ((G_lambda(i,j,k) - G_lambda(i-1,j,k))     &
                               *recip_del_lambda(i)                      &
                             - (G_lambda_eta_k_2D(i,j)                   &
                                                  - G_lambda_eta_km))    &
                               *recip_r_squared_delz(i,j,k)              &
                              - beta*ap_3D(i,j,k)                        &
                               *(field(i,j,k+1) - field(i,j,k))          &
                              + beta*am_3D(i,j,k)                        &
                               *(field(i,j,k)   - field(i,j,k-1)) 

                END DO
              END DO
            END DO

            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie
                G_lambda_eta_km = G_lambda_eta_k_2D(i,j)

                G_lambda_eta_k_2D(i,j) = 0.0

                b_3D(i,j,k) = field(i,j,k) +                             &
                            ((G_lambda(i,j,k) - G_lambda(i-1,j,k))       &
                             *recip_del_lambda(i)                        &
                           - (G_lambda_eta_k_2D(i,j) - G_lambda_eta_km)) &
                             *recip_r_squared_delz(i,j,k)                &
                            + beta*am_3D(i,j,k)                          &
                            *(field(i,j,k) - field(i,j,k-1)) 
              END DO
            END DO
! === trisolve ===========================================================
! Inlining of trisolve
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie
                field(i,j,1) = b_3D(i,j,1)*denom(i,j,1)
              END DO
            END DO

            DO k = 2,metric_level-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  field(i,j,k) = (b_3D(i,j,k)-field(i,j,k-1)             &
                                 *am_3D(i,j,k))*denom(i,j,k)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
            k = metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  field(i,j,k) = (b_3D(i,j,k) - field(i,j,k-1)           &
                                 *am_3D(i,j,k))*denom(i,j,k)
                END DO
              END DO
! Backward solve sweep
            DO k = metric_level-1,1,-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  field(i,j,k) = field(i,j,k)-c_new(i,j,k)*field(i,j,k+1)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
! ========================================================================
! === Upper level ========================================================
            is = field_dims_offi_s %i_start
            ie = field_dims_offi   %i_end 
            is2 = field_dims%i_start
            ie2 = field_dims%i_end

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
            DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie 

                  G_lambda(i,j,k) = ( field(i+1,j,k) - field(i,j,k) )    &
                                   /( lambda(i+1)    - lambda(i)    )    &
                                   * delta_lambda**2*diff_coeff(i,j) 
                END DO
                DO i = is2,ie2
                  field(i,j,k) = field(i,j,k)                            &
                              + (G_lambda(i,j,k) - G_lambda(i-1,j,k))    &
                                *recip_del_lambda(i)
                END DO
              END DO
            END DO
!$OMP END PARALLEL DO
! ------------------------------------------------------------------------
! --- w fields ---
! --- As for p fields except bottom boundary condition is more complex
          CASE ( fld_type_w )          
            is = field_dims%i_start 
            ie = field_dims%i_end    
! === Compute Matrix a ===================================================
            IF ( L_init_matrix ) THEN
              DO i = is,ie
                recip_del_lambda(i) = 1.0                                &
                                      /(lambda_offi(i) - lambda_offi(i-1))
              END DO

              k = kl
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie        
                  diff_coeff_av_2D(i,j) = 0.5*(diff_coeff(i,j)           &
                                             + diff_coeff(i-1,j))        &
                                          *delta_lambda**2               &
                                          *recip_del_lambda(i)**2

                  del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)              &
                                            -r_levels_offi(i-1,j,k+1))   &
                                  + (r_levels_offi(i,j,k)                &
                                            -r_levels_offi(i-1,j,k)))  

                  ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)      &
                                 *del_r_av_z**2                          &
                                 /(r_levels(i,j,k+1)-r_levels(i,j,k))    &
                                 *diff_coeff_av_2D(i,j)                  &
                                 *r_levels_offk(i,j,k+1)**2
                END DO
              END DO
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(del_r_av_z, &
!$OMP& i, j, k)
              DO k = kl+1, metric_level-1
                DO j = j_sweep_begin, j_sweep_end
                  DO i = is,ie    
                    del_r_av_z = 0.5*((r_levels_offi(i,j,k+1)            &
                                              -r_levels_offi(i-1,j,k+1)) &
                                    + (r_levels_offi(i,j,k)              &
                                              -r_levels_offi(i-1,j,k)))  

                    ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)    &
                                   *del_r_av_z**2                        &
                                   /(r_levels(i,j,k+1)-r_levels(i,j,k))  &
                                   *diff_coeff_av_2D(i,j)                &
                                   *r_levels_offk(i,j,k+1)**2

                    del_r_av_z = 0.5*((r_levels_offi(i,j,k)              &
                                              -r_levels_offi(i-1,j,k))   &
                                    + (r_levels_offi(i,j,k-1)            &
                                              -r_levels_offi(i-1,j,k-1)))

                    am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)    &
                                   *del_r_av_z**2                        &
                                   /(r_levels(i,j,k)-r_levels(i,j,k-1))  &
                                   *diff_coeff_av_2D(i,j)                &
                                   *r_levels_offk(i,j,k)**2
                  END DO 
                END DO
              END DO
!$OMP END PARALLEL DO
              k = metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie 
                
                  del_r_av_z = 0.5*((r_levels_offi(i,j,k)                &
                                              -r_levels_offi(i-1,j,k))   &
                                   + (r_levels_offi(i,j,k-1)             &
                                              -r_levels_offi(i-1,j,k-1)))

                  am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)      &
                                 *del_r_av_z**2                          &
                                 /(r_levels(i,j,k)-r_levels(i,j,k-1))    &
                                 *diff_coeff_av_2D(i,j)                  &
                                 *r_levels_offk(i,j,k)**2

                END DO
              END DO
! Inversion of tridiagonal matrix
              k = kl
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/(1.0 - ap_3D(i,j,k))
                  c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                END DO
              END DO
              DO k = kl+1, metric_level-1
                DO j =  j_sweep_begin, j_sweep_end
                  DO i = is,ie
                    denom(i,j,k) = 1.0/((1.0 - ap_3D(i,j,k)              &
                                   - am_3D(i,j,k))                       &
                                   - c_new(i,j,k-1)*am_3D(i,j,k))
                    c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                  END DO
                END DO
              END DO
              k = metric_level
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/((1.0 - am_3D(i,j,k))               &
                                 - c_new(i,j,k-1)*am_3D(i,j,k))
                END DO
              END DO
              L_init_matrix = .FALSE.
            ENDIF
! ========================================================================
! === G_Lambda ===========================================================
! Bottom boundary
            k = kl
            DO j = j_sweep_begin, j_sweep_end
              is3 = field_dims_offi_s%i_start 
              ie3 = field_dims_offi%i_end
              DO i = is3,ie3

! This is G_lambda stored on (i,j-1/2,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                diff_coeff_2D(i,j) = diff_coeff(i,j)*delta_lambda**2     &
                                     /( lambda(i+1) - lambda(i)) 

! New dynamics never diffuses on surface
                  del_z =  r_levels_offk(i,j,k+base_k) -                 &
                           r_levels_offk(i,j,k+base_k-1)    

                f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k+1)*(               &
                                    wt_x_on2off(i,j) *field(i,j,k+1)     &
                             + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))  &
                             + (1.0-wt_z_on2off(i,j,k+1))*(              &
                                    wt_x_on2off(i,j) *field(i,j,k)       &
                             + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                f_av_xz_km= wt_x_on2off(i,j) *field(i,j,k) +             &
                             (1.0-wt_x_on2off(i,j))*field(i+1,j,k)  

                G_lambda(i,j,k) = ((field(i+1,j,k) - field(i,j,k))*del_z &
                                - (f_av_xz_k_2D(i,j) - f_av_xz_km)       &
                                 *(r_levels(i+1,j,k) - r_levels(i,j,k))) &
                                 *diff_coeff_2D(i,j)                     &
                                 *r_levels_offi(i,j,k)**2
              END DO
            END DO 
 
            DO k = kl+1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                is3 = field_dims_offi_s%i_start 
                ie3 = field_dims_offi%i_end
                DO i = is3,ie3
! This is G_lambda stored on (i,j-1/2,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                  f_av_xz_km = f_av_xz_k_2D(i,j)

                  f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k+1)*(             &
                                      wt_x_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))&
                               + (1.0-wt_z_on2off(i,j,k+1))*(            &
                                      wt_x_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                  G_lambda(i,j,k) = ((field(i+1,j,k) - field(i,j,k))     &
                                    *(r_levels_offk(i,j,k+base_k)        &
                                    - r_levels_offk(i,j,k+base_k-1))     &
                                    - (f_av_xz_k_2D(i,j) - f_av_xz_km)   &
                                    *(r_levels(i+1,j,k)                  &
                                                    - r_levels(i,j,k)))  &
                                   *diff_coeff_2D(i,j)                   &
                                   *r_levels_offi(i,j,k)**2
                END DO
              END DO
            END DO
 
            k = metric_level
            IF ( metric_level == field_dims%k_end ) THEN 
              DO j = j_sweep_begin, j_sweep_end
                is3 = field_dims_offi_s%i_start 
                ie3 = field_dims_offi%i_end
                DO i = is3,ie3
! This is G_lambda stored on (i,j-1/2,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                  f_av_xz_km = f_av_xz_k_2D(i,j)

                  f_av_xz_k_2D(i,j) = wt_x_on2off(i,j) *field(i,j,k) +   &
                                 (1.0-wt_x_on2off(i,j))*field(i+1,j,k) 

                  G_lambda(i,j,k) = ((field(i+1,j,k) - field(i,j,k))     &
                                    *(r_levels_offk(i,j,k+base_k)        &
                                    - r_levels_offk(i,j,k+base_k-1))     &
                                   - (f_av_xz_k_2D(i,j) - f_av_xz_km)    &
                                    *(r_levels(i+1,j,k)                  &
                                                    - r_levels(i,j,k)))  &
                                    *diff_coeff_2D(i,j)                  &
                                    *r_levels_offi(i,j,k)**2
                END DO
              END DO
            ELSE
              DO j = j_sweep_begin, j_sweep_end
                is3 = field_dims_offi_s%i_start 
                ie3 = field_dims_offi%i_end
                DO i = is3,ie3

! This is G_lambda stored on (i,j-1/2,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)

                  f_av_xz_km = f_av_xz_k_2D(i,j)

                  f_av_xz_k_2D(i,j) = wt_z_on2off(i,j,k+1)*(             &
                                      wt_x_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k+1))&
                               + (1.0-wt_z_on2off(i,j,k+1))*(            &
                                      wt_x_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                  G_lambda(i,j,k) = ((field(i+1,j,k) - field(i,j,k))     &
                                    *(r_levels_offk(i,j,k+base_k)        &
                                    - r_levels_offk(i,j,k+base_k-1))     &
                                    - (f_av_xz_k_2D(i,j) - f_av_xz_km)   &
                                    *(r_levels(i+1,j,k)                  &
                                                    - r_levels(i,j,k)))  &
                                    *diff_coeff_2D(i,j)                  &
                                    *r_levels_offi(i,j,k)**2
                END DO
              END DO
            END IF
! ========================================================================
! === Compute RHS b ======================================================

            k = kl  
            kp = k + 1
            DO j = j_sweep_begin, j_sweep_end
              i = field_dims%i_start - 1
              f_av_xz_i = wt_z_on2off(i,j,kp)*(                          &
                          wt_x_on2off(i,j)*field(i,j,kp)                 &
                   + (1.0-wt_x_on2off(i,j))*field(i+1,j,kp))             &
                   + (1.0-wt_z_on2off(i,j,kp))*(                         &
                          wt_x_on2off(i,j)*field(i,j,k)                  &
                   + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

              DO i = field_dims%i_start,field_dims%i_end
                ! RHS term
                ! This is G_lambda_eta on p points
                f_av_xz_im = f_av_xz_i

                f_av_xz_i = wt_z_on2off(i,j,kp)*(                        &
                            wt_x_on2off(i,j)*field(i,j,kp)               &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,kp))           &
                     + (1.0-wt_z_on2off(i,j,kp))*(                       &
                            wt_x_on2off(i,j)*field(i,j,k)                &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                G_lambda_eta_km = 0.0

                G_lambda_eta_k_2D(i,j) = (( f_av_xz_i - f_av_xz_im ))    &
                                         *0.5*((r_levels_offi(i,j,kp)    &
                                              - r_levels_offi(i-1,j,kp)) &
                                         + (r_levels_offi(i,j,k)         &
                                               - r_levels_offi(i-1,j,k)))&
                                         *diff_coeff_av_2D(i,j)          &
                                         *r_levels_offk(i,j,kp)**2

                b_3D(i,j,k) = field(i,j,k) +                             &
                             ((G_lambda(i,j,k) - G_lambda(i-1,j,k))      &
                             *recip_del_lambda(i)                        &
                           - (G_lambda_eta_k_2D(i,j) - G_lambda_eta_km)) &
                             *recip_r_squared_delz(i,j,k)                &
                           - beta*ap_3D(i,j,k)                           &
                             *(field(i,j,k+1) - field(i,j,k))
              END DO
            END DO

            DO k = kl+1, metric_level-1
              kp = k + 1
              DO j = j_sweep_begin, j_sweep_end
                i = field_dims%i_start - 1
                f_av_xz_i = wt_z_on2off(i,j,kp)*(                        &
                            wt_x_on2off(i,j)*field(i,j,kp)               &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,kp))           &
                     + (1.0-wt_z_on2off(i,j,kp))*(                       &
                            wt_x_on2off(i,j)*field(i,j,k)                &
                     + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                DO i = field_dims%i_start,field_dims%i_end
                  ! This is G_lambda_eta on p points
                  f_av_xz_im = f_av_xz_i

                  f_av_xz_i = wt_z_on2off(i,j,kp)*(                      &
                              wt_x_on2off(i,j)*field(i,j,kp)             &
                       + (1.0-wt_x_on2off(i,j))*field(i+1,j,kp))         &
                       + (1.0-wt_z_on2off(i,j,kp))*(                     &
                              wt_x_on2off(i,j)*field(i,j,k)              &
                       + (1.0-wt_x_on2off(i,j))*field(i+1,j,k) )

                  G_lambda_eta_km = G_lambda_eta_k_2D(i,j)

                  G_lambda_eta_k_2D(i,j) = (( f_av_xz_i - f_av_xz_im ))  &
                                          *0.5*((r_levels_offi(i,j,kp)   &
                                               - r_levels_offi(i-1,j,kp))&
                                          + (r_levels_offi(i,j,k)        &
                                             - r_levels_offi(i-1,j,k)))  &
                                          *diff_coeff_av_2D(i,j)         &
                                          *r_levels_offk(i,j,kp)**2

                  b_3D(i,j,k) = field(i,j,k) +                           &
                              ((G_lambda(i,j,k) - G_lambda(i-1,j,k))     &
                              *recip_del_lambda(i)                       &
                             - (G_lambda_eta_k_2D(i,j)                   &
                                              - G_lambda_eta_km))        &
                              *recip_r_squared_delz(i,j,k)               &
                             - beta*ap_3D(i,j,k)                         &
                              *(field(i,j,k+1) - field(i,j,k))           &
                             + beta*am_3D(i,j,k)                         &
                              *(field(i,j,k)   - field(i,j,k-1)) 
                END DO
              END DO
            END DO

            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                ! This is G_lambda_eta on p points
                ! RHS term
                G_lambda_eta_km = G_lambda_eta_k_2D(i,j)

                G_lambda_eta_k_2D(i,j) = 0.0

                b_3D(i,j,k) = field(i,j,k) +                             &
                            ((G_lambda(i,j,k) - G_lambda(i-1,j,k))       &
                             *recip_del_lambda(i)                        &
                           - (G_lambda_eta_k_2D(i,j) - G_lambda_eta_km)) &
                             *recip_r_squared_delz(i,j,k)                &
                            + beta*am_3D(i,j,k)*                         &
                             (field(i,j,k)   - field(i,j,k-1)) 
              END DO
            END DO
!=========================================================================
! === trisolve ===========================================================
! Inlining of trisolve
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,kl) = b_3D(i,j,kl)*denom(i,j,kl)
              END DO
            END DO

            DO k = kl+1,metric_level-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  field(i,j,k) = (b_3D(i,j,k)-field(i,j,k-1)             &
                                 *am_3D(i,j,k))*denom(i,j,k)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = (b_3D(i,j,k) - field(i,j,k-1)             &
                               *am_3D(i,j,k))*denom(i,j,k)
              END DO
            END DO
! Backward solve sweep
            DO k = metric_level-1,kl,-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                   field(i,j,k) = field(i,j,k)-c_new(i,j,k)*field(i,j,k+1)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
! ========================================================================
! === Upper Level ========================================================
            is = field_dims_offi_s%i_start
            ie = field_dims_offi%i_end 
            is2 = field_dims%i_start
            ie2 = field_dims%i_end
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
            DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  G_lambda(i,j,k) = ( field(i+1,j,k) - field(i,j,k) )    &
                                   /( lambda(i+1)    - lambda(i)    )    &
                                   * delta_lambda**2*diff_coeff(i,j) 
                END DO
                DO i = is2,ie2
                  field(i,j,k) = field(i,j,k)                            &
                              + (G_lambda(i,j,k) - G_lambda(i-1,j,k))    &
                                *recip_del_lambda(i)                   
                END DO
              END DO
            END DO
!$OMP END PARALLEL DO
! ------------------------------------------------------------------------
! --- Default ---
          CASE DEFAULT
            CALL ereport('pofil_zlevel',fld_type, 'Invalid field type')
          END SELECT
!-------------------------------------------------------------------------
! IF n_procy=1, Northern hemisphere needs to be done after S. Hem
          IF ( L_cycle ) THEN
            L_init_matrix = .TRUE. ! Recompute tri-diagonal matrix
            j_sweep_store = j_sweep_begin
            IF ( i_sweep == 1 .AND. L_combine ) THEN
              j_sweep_begin = j_sweep_end + 1
            ELSE
              j_sweep_begin = field_dims%j_end - j_sweep_end + 1
            END IF ! i_sweep == 1 .AND. L_combine
            j_sweep_end = field_dims%j_end - j_sweep_store + 1
            L_cycle = .false.
        CYCLE
          END IF ! L_cycle
! Reset pointers for next sweep
! either because N Hem has just been done or
! 1st sweep was combined filter and diffusion
          j_sweep_begin = idx_begin(i_filter)
          j_sweep_end = idx_end(i_filter)
          IF( n_procy == 1 ) L_cycle = .true.

          i_sweep = i_sweep + 1

          IF ( i_sweep > sweeps(i_filter) ) THEN
! Full Swap bounds and Exit
! DEPENDS ON: swap_bounds
            CALL Swap_Bounds(                                            &
                             field, row_length, in_rows, levels,         &
                             offx, offy, swp_bnds_fld_type, L_vector) 


            EXIT
          ELSE
! Swap bounds East-West
! DEPENDS ON: swap_bounds
            CALL Swap_Bounds(                                            &
                             field, row_length, in_rows+2*offy, levels,  &
                             offx, 0, swp_bnds_fld_type, L_vector) 
          END IF
        CYCLE
        END DO ! sweeping  loop
      END DO   !  i_filter = 1, num_pass
! ------------------------------------------------------------------------
! Section 4.2   NS diffusion
! ------------------------------------------------------------------------
      IF( L_diff) THEN
        j_sweep_begin = idx_begin(0)
        j_sweep_end = idx_end(0)
        IF ( n_procy == 1 ) THEN
          j_sweep_end = field_dims%j_end - j_sweep_begin + 1
        END IF

        SELECT CASE ( fld_type ) 
! ------------------------------------------------------------------------
! --- U fields ---
        CASE ( fld_type_u )

          is = field_dims%i_start
          ie = field_dims%i_end

        ! settings for diffusion of u
          IF ( L_u_ns_diff ) THEN
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
            DO k = 1, active_levels
              DO j = j_sweep_begin-1, j_sweep_end+1
                DO i = is,ie
                  field(i,j,k) = field(i,j,k)/csphi_on(j)
                END DO
              END DO
            END DO
!$OMP END PARALLEL DO

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

!===== Compute a =========================================================
          DO j = j_sweep_begin, j_sweep_end
            recip_del_phi(j) = 1.0/(phi_offj(j) - phi_offj(j-1))
            recip_csphi(j) = 1.0/csphi_on_u(j)
          END DO

          k = 1
          DO j = j_sweep_begin, j_sweep_end
            DO i = is,ie
              diff_coeff_av_2D(i,j) = 0.5*(diff_coeff_phi(i,j)           &
                                          +diff_coeff_phi(i,j-1))        &
                                      *delta_phi**2*recip_del_phi(j)**2

              del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)                  &
                                         -r_levels_offj(i,j-1,k+1))      &
                                + (r_levels_offj(i,j,k)                  &
                                          -r_levels_offj(i,j-1,k))) 

              ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)          &
                             *del_r_av_z**2                              &
                             /(r_levels(i,j,k+1)-r_levels(i,j,k))        &
                             *diff_coeff_av_2D(i,j)                      &
                             *r_levels_offk(i,j,k)**2

            END DO
          END DO
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(del_r_av_z, &
!$OMP& i, j, k)
          DO k = 2, metric_level-1
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie

                del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)                &
                                           -r_levels_offj(i,j-1,k+1))    &
                                  + (r_levels_offj(i,j,k)                &
                                            -r_levels_offj(i,j-1,k))) 

                ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)        &
                               *del_r_av_z**2                            &
                               /(r_levels(i,j,k+1)-r_levels(i,j,k))      &
                               *diff_coeff_av_2D(i,j)                    &
                               *r_levels_offk(i,j,k)**2

                del_r_av_z = 0.5*((r_levels_offj(i,j,k)                  &
                                           -r_levels_offj(i,j-1,k))      &
                                  + (r_levels_offj(i,j,k-1)              &
                                            -r_levels_offj(i,j-1,k-1))) 

                am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)        &
                               *del_r_av_z**2                            &
                               /(r_levels(i,j,k)-r_levels(i,j,k-1))      &
                               *diff_coeff_av_2D(i,j)                    &
                               *r_levels_offk(i,j,k-1)**2

              END DO
            END DO
          END DO
!$OMP END PARALLEL DO
          k = metric_level
          DO j = j_sweep_begin, j_sweep_end
            DO i = is,ie

              del_r_av_z = 0.5*((r_levels_offj(i,j,k)                    &
                                         -r_levels_offj(i,j-1,k))        &
                                + (r_levels_offj(i,j,k-1)                &
                                          -r_levels_offj(i,j-1,k-1))) 

              am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)          &
                             *del_r_av_z**2                              &
                             /(r_levels(i,j,k)-r_levels(i,j,k-1))        &
                             *diff_coeff_av_2D(i,j)                      &
                             *r_levels_offk(i,j,k-1)**2
            END DO
          END DO

! Inversion of tridiagonal matrix
          k = 1
          DO j =  j_sweep_begin, j_sweep_end
            DO i = is,ie
              denom(i,j,k) = 1.0/(u_fac(j) - ap_3D(i,j,k))
              c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
            END DO
          END DO
          DO k = 2, metric_level-1
            DO j =  j_sweep_begin, j_sweep_end
              DO i = is,ie
                denom(i,j,k) = 1.0/((u_fac(j) - ap_3D(i,j,k)             &
                               - am_3D(i,j,k))                           &
                               - c_new(i,j,k-1)*am_3D(i,j,k))
                c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
              END DO
            END DO
          END DO
          k = metric_level
          DO j =  j_sweep_begin, j_sweep_end
            DO i = is,ie
              denom(i,j,k) = 1.0/((u_fac(j)                              &
                             - am_3D(i,j,k))                             &
                             - c_new(i,j,k-1)*am_3D(i,j,k))
            END DO
          END DO
!=========================================================================
!===== G_Phi =============================================================
          k = 1
          DO j = j_sweep_begin-1, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end
! This is G_phi stored on (i,j,k-1/2) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
              diff_coeff_2D(i,j) = diff_coeff_phi(i,j)*delta_phi**2      &
                                   /( phi(j+1) - phi(j))*csphi_off_u(j) 

              f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                   &
                                  wt_y_on2off(i,j) *field(i,j,k+1)       &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))    &
                           + (1.0-wt_z_on2off(i,j,k))*(                  &
                                  wt_y_on2off(i,j) *field(i,j,k)         &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

              f_av_yz_km=  wt_y_on2off(i,j) *field(i,j,k)                &
                    + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) 


              G_phi(i,j,k) = (( field(i,j+1,k) - field(i,j,k))           &
                               *(r_levels_offk(i,j,k+base_k) -           &
                                 r_levels_offk(i,j,k+base_k-1))          &
                               - (f_av_yz_k_2D(i,j) - f_av_yz_km)        &
                               *(r_levels(i,j+1,k) - r_levels(i,j,k)))   &
                               *diff_coeff_2D(i,j)                       &
                               *r_levels_offj(i,j,k)**2
                                 
            END DO
          END DO
          DO k = 2, metric_level-1
            DO j = j_sweep_begin-1, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on (i,j,k-1/2) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                f_av_yz_km = f_av_yz_k_2D(i,j)

                f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                 &
                                    wt_y_on2off(i,j) *field(i,j,k+1)     &
                             + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))  &
                             + (1.0-wt_z_on2off(i,j,k))*(                &
                                    wt_y_on2off(i,j) *field(i,j,k)       &
                             + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                G_phi(i,j,k) = (( field(i,j+1,k) - field(i,j,k))         &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                 - ( f_av_yz_k_2D(i,j) - f_av_yz_km   )  &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                 *diff_coeff_2D(i,j)                     &
                                 *r_levels_offj(i,j,k)**2
              END DO
            END DO
          END DO
          k = metric_level
          IF (metric_level == field_dims%k_end  ) THEN
            DO j = j_sweep_begin-1, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on (i,j,k-1/2) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                f_av_yz_km = f_av_yz_k_2D(i,j)

                f_av_yz_k_2D(i,j) = wt_y_on2off(i,j) *field(i,j,k) +     &
                               (1.0-wt_y_on2off(i,j))*field(i,j+1,k) 

                G_phi(i,j,k) = (( field(i,j+1,k) - field(i,j,k))         &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                 - ( f_av_yz_k_2D(i,j) - f_av_yz_km   )  &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                 *diff_coeff_2D(i,j)                     &
                                 *r_levels_offj(i,j,k)**2
              END DO
            END DO
          ELSE
            DO j = j_sweep_begin-1, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on (i,j,k-1/2) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                f_av_yz_km = f_av_yz_k_2D(i,j)

                f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                 &
                                    wt_y_on2off(i,j) *field(i,j,k+1)     &
                             + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))  &
                             + (1.0-wt_z_on2off(i,j,k))*(                &
                                    wt_y_on2off(i,j) *field(i,j,k)       &
                             + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                G_phi(i,j,k) = (( field(i,j+1,k) - field(i,j,k))         &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                 - ( f_av_yz_k_2D(i,j) - f_av_yz_km)     &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                 *diff_coeff_2D(i,j)                     &
                                 *r_levels_offj(i,j,k)**2
              END DO
            END DO
          END IF
          IF( at_extremity(PSouth) ) THEN
            DO k = 1, metric_level
              DO i = field_dims%i_start,field_dims%i_end
                G_phi(i,j_sweep_begin-1,k) = 0.0 
              END DO              
            END DO
          END IF 
          IF( at_extremity(PNorth) ) THEN
            DO k = 1, metric_level
              DO i = field_dims%i_start,field_dims%i_end
                G_phi(i,j_sweep_end,k) = 0.0 
              END DO              
            END DO
          END IF            
!=========================================================================
!===== Compute RHS b =====================================================

          k = 1
          j = j_sweep_begin - 1
          DO i = field_dims%i_start,field_dims%i_end
            f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                       &
                              wt_y_on2off(i,j)*field(i,j,k+1)            &
                       + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))        &
                       + (1.0-wt_z_on2off(i,j,k))*(                      &
                              wt_y_on2off(i,j)*field(i,j,k)              &
                       + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) 
          END DO
          DO j = j_sweep_begin, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end

              ! RHS term
              ! This is G_phi_eta on (i,j-1/2,k) points   
              f_av_yz_jm = f_av_yz_j_1D(i)
                 
              f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                     &
                                wt_y_on2off(i,j)*field(i,j,k+1)          &
                         + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))      &
                         + (1.0-wt_z_on2off(i,j,k))*(                    &
                                wt_y_on2off(i,j)*field(i,j,k)            &
                         + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) 
          
              G_phi_eta_km = 0.0

              G_phi_eta_k_2D(i,j) = (( f_av_yz_j_1D(i) - f_av_yz_jm ))   &
                                    *0.5*((r_levels_offj(i,j,k+1)        &
                                              -r_levels_offj(i,j-1,k+1)) &
                                     + (r_levels_offj(i,j,k)             &
                                               -r_levels_offj(i,j-1,k))) &
                                    *diff_coeff_av_2D(i,j)               &
                                    *r_levels_offk(i,j,k)**2*csphi_on_u(j)

              b_3D(i,j,k) = field(i,j,k)*u_fac(j) +                      &
                          ((G_phi(i,j,k) - G_phi(i,j-1,k))               &
                           *recip_del_phi(j)                             &
                         - (G_phi_eta_k_2D(i,j)  - G_phi_eta_km))        &
                          * recip_r_squared_delz(i,j,k)*recip_csphi(j)   &
                          - beta*ap_3D(i,j,k)                            &
                          *(field(i,j,k+1) - field(i,j,k)) 

            END DO
          END DO

          DO k = 2, metric_level-1
            j = j_sweep_begin - 1
            DO i = field_dims%i_start,field_dims%i_end
              f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                     &
                                wt_y_on2off(i,j)*field(i,j,k+1)          &
                         + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))      &
                         + (1.0-wt_z_on2off(i,j,k))*(                    &
                                wt_y_on2off(i,j)*field(i,j,k)            &
                         + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) 
            END DO
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
               ! RHS term
               ! This is G_phi_eta on (i,j-1/2,k) points    
                f_av_yz_jm = f_av_yz_j_1D(i)
                
                f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                   &
                                  wt_y_on2off(i,j)*field(i,j,k+1)        &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))    &
                           + (1.0-wt_z_on2off(i,j,k))*(                  &
                                  wt_y_on2off(i,j)*field(i,j,k)          &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) 

                G_phi_eta_km = G_phi_eta_k_2D(i,j)       

                G_phi_eta_k_2D(i,j) = ((f_av_yz_j_1D(i) - f_av_yz_jm ))  &
                                       *0.5*((r_levels_offj(i,j,k+1)     &
                                             -r_levels_offj(i,j-1,k+1))  &
                                       + (r_levels_offj(i,j,k)           &
                                             -r_levels_offj(i,j-1,k)))   &
                                        *diff_coeff_av_2D(i,j)           &
                                        *r_levels_offk(i,j,k)**2         &
                                        *csphi_on_u(j)

                b_3D(i,j,k) = field(i,j,k)*u_fac(j) +                    &
                            ((G_phi(i,j,k) - G_phi(i,j-1,k))             &
                             *recip_del_phi(j)                           &
                           - (G_phi_eta_k_2D(i,j)  - G_phi_eta_km))      &
                            * recip_r_squared_delz(i,j,k)*recip_csphi(j) & 
                            - beta*ap_3D(i,j,k)                          &
                            *(field(i,j,k+1) - field(i,j,k))             &
                            + beta*am_3D(i,j,k)                          &
                            *(field(i,j,k)   - field(i,j,k-1))
              END DO
            END DO
          END DO

          k = metric_level
          DO j = j_sweep_begin, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end

              ! RHS term
              G_phi_eta_km = G_phi_eta_k_2D(i,j)       

              G_phi_eta_k_2D(i,j) = 0.0

              b_3D(i,j,k) = field(i,j,k)*u_fac(j) +                      &
                          ((G_phi(i,j,k) - G_phi(i,j-1,k))               &
                           *recip_del_phi(j)                             &
                         - (G_phi_eta_k_2D(i,j)  - G_phi_eta_km))        &
                          * recip_r_squared_delz(i,j,k)*recip_csphi(j)   &
                          + beta*am_3D(i,j,k)                            &
                          *(field(i,j,k)- field(i,j,k-1)) 
            END DO
          END DO
!=========================================================================
!===== Trisolve ==========================================================
! Inlining of trisolve
          DO j = j_sweep_begin, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end
              field(i,j,1) = b_3D(i,j,1)*denom(i,j,1)
            END DO
          END DO

          DO k = 2,metric_level-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = (b_3D(i,j,k)-field(i,j,k-1)*am_3D(i,j,k)) &
                               *denom(i,j,k)
              END DO
            END DO
!$OMP END PARALLEL DO
          END DO
          k = metric_level
          DO j = j_sweep_begin, j_sweep_end
            DO i = field_dims%i_start,field_dims%i_end
              field(i,j,k) = (b_3D(i,j,k) - field(i,j,k-1)*am_3D(i,j,k)) &
                             *denom(i,j,k)
            END DO
          END DO

! Backward solve sweep
          DO k = metric_level-1,1,-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,k) = field(i,j,k)-c_new(i,j,k)*field(i,j,k+1)
              END DO
            END DO
!$OMP END PARALLEL DO
          END DO
!=========================================================================
!===== Upper level =======================================================
          is = field_dims%i_start
          ie = field_dims%i_end
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
          DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
            DO j = j_sweep_begin-1, j_sweep_end
              DO i = is,ie
                G_phi(i,j,k) = ( field(i,j+1,k) - field(i,j,k) )         &
                              /( phi(j+1)    - phi(j)    )               & 
                               * delta_phi**2*diff_coeff_phi(i,j)        &
                               * csphi_off_u(j)     
              END DO
            END DO
            IF( at_extremity(PSouth) ) THEN
              DO i = is,ie
                G_phi(i,j_sweep_begin-1,k) = 0.0 
              END DO
            END IF 
            IF( at_extremity(PNorth) ) THEN
              DO i = is,ie
                G_phi(i,j_sweep_end,k) = 0.0             
              END DO
            END IF            
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie
                field(i,j,k) = field(i,j,k)                              &
                            + (G_phi(i,j,k) - G_phi(i,j-1,k))            &
                              *recip_del_phi(j)*recip_csphi(j)

              END DO
            END DO
          END DO
!$OMP END PARALLEL DO
!=========================================================================
! ------------------------------------------------------------------------
! --- V fields ---
        CASE ( fld_type_v )
          j_v_b = j_sweep_begin
          j_v_e = j_sweep_end 
          IF( at_extremity(PSouth) ) j_v_b = j_sweep_begin - 1
          IF( at_extremity(PNorth) ) j_v_e = j_sweep_end   + 1
          is = field_dims%i_start
          ie = field_dims%i_end
!===== Compute a =========================================================
          k = 1
          DO j = j_v_b, j_v_e
            recip_del_phi(j) =  1.0/(phi_offj(j+1) - phi_offj(j))
            recip_csphi(j) = 1.0/csphi_on(j)
            DO i = is,ie
              diff_coeff_av_2D(i,j) = 0.5*(diff_coeff_phi(i,j+1)         &
                                         + diff_coeff_phi(i,j))          &
                                         *delta_phi**2                   &
                                         *recip_del_phi(j)**2

              del_r_av_z = 0.5*((r_levels_offj(i,j+1,k+1)                &
                                           -r_levels_offj(i,j,k+1))      &
                              + (r_levels_offj(i,j+1,k)                  &
                                           -r_levels_offj(i,j,k)))  

              ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)          &
                             *del_r_av_z**2                              &
                            /(r_levels(i,j,k+1)-r_levels(i,j,k))         &
                             *diff_coeff_av_2D(i,j)                      &
                             *r_levels_offk(i,j,k)**2

!               am_3D(i,j,k) = 0.0
            END DO
          END DO

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(del_r_av_z, &
!$OMP& i, j, k)
          DO k = 2, metric_level-1
            DO j = j_v_b, j_v_e
              DO i = is,ie

                del_r_av_z = 0.5*((r_levels_offj(i,j+1,k+1)              &
                                                -r_levels_offj(i,j,k+1)) &
                                 + (r_levels_offj(i,j+1,k)               &
                                                -r_levels_offj(i,j,k)))  

                ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)        &
                               *del_r_av_z**2                            &
                              /(r_levels(i,j,k+1)-r_levels(i,j,k))       &
                               *diff_coeff_av_2D(i,j)                    &
                               *r_levels_offk(i,j,k)**2

                del_r_av_z = 0.5*((r_levels_offj(i,j+1,k)                &
                                               -r_levels_offj(i,j,k))    &
                                + (r_levels_offj(i,j+1,k-1)              &
                                               -r_levels_offj(i,j,k-1)))  

                am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)        &
                               *del_r_av_z**2                            &
                              /(r_levels(i,j,k)-r_levels(i,j,k-1))       &
                               *diff_coeff_av_2D(i,j)                    &
                               *r_levels_offk(i,j,k-1)**2
              END DO
            END DO
          END DO
!$OMP END PARALLEL DO
          k = metric_level
          DO j = j_v_b, j_v_e
            DO i = is,ie

              del_r_av_z = 0.5*((r_levels_offj(i,j+1,k)                  &
                                           -r_levels_offj(i,j,k))        &
                              + (r_levels_offj(i,j+1,k-1)                &
                                           -r_levels_offj(i,j,k-1)))  

              am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)          &
                             *del_r_av_z**2                              &
                            /(r_levels(i,j,k)-r_levels(i,j,k-1))         &
                             *diff_coeff_av_2D(i,j)                      &
                             *r_levels_offk(i,j,k-1)**2
            END DO
          END DO

! Inversion of tridiagonal matrix
          k = 1
          DO j = j_v_b, j_v_e
            DO i = is,ie
              denom(i,j,k) = 1.0/(1.0 - ap_3D(i,j,k))
              c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
            END DO
          END DO
          DO k = 2, metric_level-1
            DO j = j_v_b, j_v_e
              DO i = is,ie
                denom(i,j,k) = 1.0/((1.0 - ap_3D(i,j,k) - am_3D(i,j,k))  &
                              -c_new(i,j,k-1)*am_3D(i,j,k))
                c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
              END DO
            END DO
          END DO
          k = metric_level
          DO j =  j_v_b, j_v_e
            DO i = is,ie
              denom(i,j,k) = 1.0/((1.0 - am_3D(i,j,k))                   &
                            -c_new(i,j,k-1)*am_3D(i,j,k))
            END DO
          END DO
!=========================================================================
!===== G_Phi =============================================================
          k = 1
          DO j = j_sweep_begin, j_sweep_end+1
            DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on p points 
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
              diff_coeff_2D(i,j) = diff_coeff_phi(i,j)*delta_phi**2      &
                                 /(phi(j) - phi(j-1))*csphi_off(j) 

              f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                   &
                                  wt_y_on2off(i,j) *field(i,j-1,k+1)     &
                           + (1.0-wt_y_on2off(i,j))*field(i,j,k+1))      &
                           + (1.0-wt_z_on2off(i,j,k))*(                  &
                                  wt_y_on2off(i,j) *field(i,j-1,k)       &
                           + (1.0-wt_y_on2off(i,j))*field(i,j,k) )

              f_av_yz_km= wt_y_on2off(i,j) *field(i,j-1,k) +             &
                     (1.0-wt_y_on2off(i,j))*field(i,j,k)

              G_phi(i,j,k) = ((field(i,j,k) - field(i,j-1,k) )           &
                             *(r_levels_offk(i,j,k+base_k)               &
                             - r_levels_offk(i,j,k+base_k-1))            &
                            - (f_av_yz_k_2D(i,j) - f_av_yz_km   )        &
                             *(r_levels(i,j,k)-r_levels(i,j-1,k)))       &
                              *diff_coeff_2D(i,j)                        &
                              *r_levels_offj(i,j,k)**2
            END DO
          END DO              
          DO k = 2, metric_level-1
            DO j = j_sweep_begin, j_sweep_end+1
              DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on p points 
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                f_av_yz_km = f_av_yz_k_2D(i,j)

                f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                 &
                                    wt_y_on2off(i,j) *field(i,j-1,k+1)   &
                             + (1.0-wt_y_on2off(i,j))*field(i,j,k+1))    &
                             + (1.0-wt_z_on2off(i,j,k))*(                &
                                    wt_y_on2off(i,j) *field(i,j-1,k)     &
                             + (1.0-wt_y_on2off(i,j))*field(i,j,k) )

                G_phi(i,j,k) = ((field(i,j,k) - field(i,j-1,k))          &
                               *(r_levels_offk(i,j,k+base_k)             &
                               - r_levels_offk(i,j,k+base_k-1))          &
                              - (f_av_yz_k_2D(i,j) - f_av_yz_km)         &
                               *(r_levels(i,j,k)-r_levels(i,j-1,k)))     &
                                *diff_coeff_2D(i,j)                      &
                                *r_levels_offj(i,j,k)**2
              END DO
            END DO              
          END DO
          k = metric_level
          IF ( metric_level == field_dims%k_end ) THEN
            DO j = j_sweep_begin, j_sweep_end+1
              DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on p points 
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                f_av_yz_km = f_av_yz_k_2D(i,j)

                f_av_yz_k_2D(i,j) = wt_y_on2off(i,j) *field(i,j-1,k) +   &
                               (1.0-wt_y_on2off(i,j))*field(i,j,k) 

                G_phi(i,j,k) = ((field(i,j,k) - field(i,j-1,k))          &
                               *(r_levels_offk(i,j,k+base_k)             &
                               - r_levels_offk(i,j,k+base_k-1))          &
                              - (f_av_yz_k_2D(i,j) - f_av_yz_km   )      &
                               *(r_levels(i,j,k)-r_levels(i,j-1,k)))     &
                                *diff_coeff_2D(i,j)                      &
                                *r_levels_offj(i,j,k)**2
              END DO
            END DO 
          ELSE
            DO j = j_sweep_begin, j_sweep_end+1
              DO i = field_dims%i_start,field_dims%i_end

! This is G_phi stored on p points 
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                f_av_yz_km = f_av_yz_k_2D(i,j)

                f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k)*(                 &
                                    wt_y_on2off(i,j) *field(i,j-1,k+1)   &
                             + (1.0-wt_y_on2off(i,j))*field(i,j,k+1))    &
                             + (1.0-wt_z_on2off(i,j,k))*(                &
                                    wt_y_on2off(i,j) *field(i,j-1,k)     &
                             + (1.0-wt_y_on2off(i,j))*field(i,j,k) )

                G_phi(i,j,k) = ((field(i,j,k) - field(i,j-1,k))          &
                               *(r_levels_offk(i,j,k+base_k) -           &
                                 r_levels_offk(i,j,k+base_k-1))          &
                              - (f_av_yz_k_2D(i,j) - f_av_yz_km   )      &
                               *(r_levels(i,j,k)-r_levels(i,j-1,k)))     &
                                *diff_coeff_2D(i,j)                      &
                                *r_levels_offj(i,j,k)**2
              END DO
            END DO 
          END IF             
          j = j_sweep_begin - 1
          IF( at_extremity(PSouth) ) THEN
            DO k = 1, metric_level
              DO i = field_dims%i_start,field_dims%i_end
                G_phi(i,j,k) = 0.0 
              END DO              
            END DO 
          END IF
          j = j_sweep_end + 2
          IF( at_extremity(PNorth) ) THEN
            DO k = 1, metric_level
              DO i = field_dims%i_start,field_dims%i_end
                G_phi(i,j,k) = 0.0 
              END DO              
            END DO
          END IF  
!=========================================================================
!===== Compute RHS b =====================================================       

          k = 1
          j = j_v_b - 1
          DO i = field_dims%i_start,field_dims%i_end
            f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                       &
                              wt_y_on2off(i,j)*field(i,j,k+1)            &
                       + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))        &
                       + (1.0-wt_z_on2off(i,j,k))*(                      &
                              wt_y_on2off(i,j)*field(i,j,k)              &
                       + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) 
          END DO
          DO j = j_v_b, j_v_e
            DO i = field_dims%i_start,field_dims%i_end
                ! RHS term
                ! This is G_phi_eta on (i-1/2,j,k) points
              f_av_yz_jm = f_av_yz_j_1D(i)

              f_av_yz_j_1D(i)= wt_z_on2off(i,j,k)*(                      &
                               wt_y_on2off(i,j)*field(i,j,k+1-base_k)    &
                        + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1-base_k))&
                        + (1.0-wt_z_on2off(i,j,k))*(                     &
                               wt_y_on2off(i,j)*field(i,j,k-base_k)      &
                        + (1.0-wt_y_on2off(i,j))*field(i,j+1,k-base_k) )

              G_phi_eta_km = 0.0

              G_phi_eta_k_2D(i,j) = ((f_av_yz_j_1D(i) - f_av_yz_jm ))    &
                                     *0.5*((r_levels_offj(i,j+1,k+1)     &
                                           -r_levels_offj(i,j,k+1))      &
                                   + (r_levels_offj(i,j+1,k)             &
                                           -r_levels_offj(i,j,k)))       &
                                     *diff_coeff_av_2D(i,j)              &
                                     *r_levels_offk(i,j,k)**2*csphi_on(j)

              b_3D(i,j,k) = field(i,j,k) +                               &
                          ((G_phi(i,j+1,k) - G_phi(i,j,k))               &
                           *recip_del_phi(j)                             &
                         - (G_phi_eta_k_2D(i,j) - G_phi_eta_km))         &
                           *recip_r_squared_delz(i,j,k)*recip_csphi(j)   &
                          - beta*ap_3D(i,j,k)                            &
                          *(field(i,j,k+1) - field(i,j,k))
            END DO 
          END DO

          DO k = 2, metric_level-1
            j = j_v_b - 1
            DO i = field_dims%i_start,field_dims%i_end
              f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                     &
                                wt_y_on2off(i,j)*field(i,j,k+1)          &
                         + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))      &
                         + (1.0-wt_z_on2off(i,j,k))*(                    &
                                wt_y_on2off(i,j)*field(i,j,k)            &
                         + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) 
            END DO
            DO j = j_v_b, j_v_e
              DO i = field_dims%i_start,field_dims%i_end
                 ! RHS term
                 ! This is G_phi_eta on (i-1/2,j,k) points
                f_av_yz_jm = f_av_yz_j_1D(i)

                f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                   &
                                  wt_y_on2off(i,j)*field(i,j,k+1-base_k) &
                           + (1.0-wt_y_on2off(i,j))                      &
                                               *field(i,j+1,k+1-base_k)) &
                           + (1.0-wt_z_on2off(i,j,k))*(                  &
                                  wt_y_on2off(i,j)*field(i,j,k-base_k)   &
                           + (1.0-wt_y_on2off(i,j))                      &
                                               *field(i,j+1,k-base_k))

                G_phi_eta_km = G_phi_eta_k_2D(i,j)

                G_phi_eta_k_2D(i,j) = ((f_av_yz_j_1D(i) - f_av_yz_jm ))  &
                                       *0.5*((r_levels_offj(i,j+1,k+1)   &
                                             -r_levels_offj(i,j,k+1))    &
                                     + (r_levels_offj(i,j+1,k)           &
                                             -r_levels_offj(i,j,k)))     &
                                       *diff_coeff_av_2D(i,j)            &
                                       *r_levels_offk(i,j,k)**2          &
                                       *csphi_on(j)

                b_3D(i,j,k) = field(i,j,k) +                             &
                            ((G_phi(i,j+1,k) - G_phi(i,j,k))             &
                             *recip_del_phi(j)                           &
                           - (G_phi_eta_k_2D(i,j) - G_phi_eta_km))       &
                             *recip_r_squared_delz(i,j,k)*recip_csphi(j) & 
                            - beta*ap_3D(i,j,k)                          &
                            *(field(i,j,k+1) - field(i,j,k))             &
                            + beta*am_3D(i,j,k)                          &
                            *(field(i,j,k)   - field(i,j,k-1))
                END DO
              END DO
            END DO

            k = metric_level
            DO j = j_v_b, j_v_e
              DO i = field_dims%i_start,field_dims%i_end
                ! RHS term
                G_phi_eta_km = G_phi_eta_k_2D(i,j)

                G_phi_eta_k_2D(i,j) = 0.0

                b_3D(i,j,k) = field(i,j,k) +                             &
                            ((G_phi(i,j+1,k) - G_phi(i,j,k))             &
                             *recip_del_phi(j)                           &
                           - (G_phi_eta_k_2D(i,j) - G_phi_eta_km))       &
                             *recip_r_squared_delz(i,j,k)*recip_csphi(j) &
                            + beta*am_3D(i,j,k)                          &
                            *(field(i,j,k)   - field(i,j,k-1))
              END DO
            END DO

!=========================================================================
!===== Trisolve ==========================================================
! Inlining of trisolve
            DO j = j_v_b, j_v_e
              DO i = is,ie
                field(i,j,1) = b_3D(i,j,1)*denom(i,j,1)
              END DO
            END DO

            DO k = 2,metric_level-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_v_b, j_v_e
                DO i = is,ie
                  field(i,j,k) = (b_3D(i,j,k)-field(i,j,k-1)             &
                                 *am_3D(i,j,k))*denom(i,j,k)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
            k = metric_level
            DO j = j_v_b, j_v_e
              DO i = is,ie
                field(i,j,k) = (b_3D(i,j,k) - field(i,j,k-1)             &
                               *am_3D(i,j,k))*denom(i,j,k)
              END DO
            END DO

! Backward solve sweep
            DO k = metric_level-1,1,-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_v_b, j_v_e
                DO i = is,ie
                  field(i,j,k) = field(i,j,k)-c_new(i,j,k)               &
                                *field(i,j,k+1)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
!=========================================================================
!===== Upper level =======================================================

            is = field_dims%i_start
            ie = field_dims%i_end 

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
            DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
              DO j = j_sweep_begin, j_sweep_end+1
                DO i = is,ie
                  G_phi(i,j,k) = ( field(i,j,k) - field(i,j-1,k) )       &
                                /( phi(j)    - phi(j-1)    )             & 
                                 * delta_phi**2*diff_coeff_phi(i,j)      &
                                 * csphi_off(j)
                END DO
              END DO
              j = j_sweep_begin - 1
              IF( at_extremity(PSouth) ) THEN
                DO i = is,ie
                  G_phi(i,j,k) = 0.0 
                END DO              
              END IF
              j = j_sweep_end + 2
              IF( at_extremity(PNorth) ) THEN
                DO i = is,ie
                  G_phi(i,j,k) = 0.0 
                END DO              
              END IF            

              DO j = j_v_b, j_v_e
                DO i = is,ie
                  field(i,j,k) = field(i,j,k)                            &
                              + (G_phi(i,j+1,k) - G_phi(i,j,k))          &
                               *recip_del_phi(j)*recip_csphi(j)
                END DO
              END DO
            END DO
!$OMP END PARALLEL DO
! ------------------------------------------------------------------------
! --- P fields ---
          CASE ( fld_type_p )   
            WRITE(6,*) ' Diffusion of fields on p points not supported'
            is = field_dims%i_start
            ie = field_dims%i_end
!===== Compute a =========================================================
            DO j = j_sweep_begin, j_sweep_end
              recip_del_phi(j) = 1.0/(phi_offj(j) - phi_offj(j-1))
              recip_csphi(j) = 1.0/csphi_on(j)
            END DO
            k = 1
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie
                diff_coeff_av_2D(i,j) = 0.5*(diff_coeff_phi(i,j)         &
                                            +diff_coeff_phi(i,j-1))      &
                                       *delta_phi**2*recip_del_phi(j)**2

                del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)                &
                                           -r_levels_offj(i,j-1,k+1))    &
                                + (r_levels_offj(i,j,k)                  &
                                           -r_levels_offj(i,j-1,k)))  

                ap_3d(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)        &
                               *del_r_av_z**2                            &
                              /(r_levels(i,j,k+1)-r_levels(i,j,k))       &
                               *diff_coeff_av_2D(i,j)                    &
                               *r_levels_offk(i,j,k)**2

              END DO
            END DO
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(del_r_av_z, &
!$OMP& i, j, k)
            DO k = 2, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie

                  del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)              &
                                             -r_levels_offj(i,j-1,k+1))  &
                                  + (r_levels_offj(i,j,k)                &
                                             -r_levels_offj(i,j-1,k)))  

                  ap_3d(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)      &
                                 *del_r_av_z**2                          &
                                /(r_levels(i,j,k+1)-r_levels(i,j,k))     &
                                 *diff_coeff_av_2D(i,j)                  &
                                 *r_levels_offk(i,j,k)**2

                  del_r_av_z = 0.5*((r_levels_offj(i,j,k)                &
                                             -r_levels_offj(i,j-1,k))    &
                                  + (r_levels_offj(i,j,k-1)              &
                                             -r_levels_offj(i,j-1,k-1)))  

                  am_3d(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)      &
                                 *del_r_av_z**2                          &
                                /(r_levels(i,j,k)-r_levels(i,j,k-1))     &
                                 *diff_coeff_av_2D(i,j)                  &
                                 *r_levels_offk(i,j,k-1)**2
                 
                END DO
              END DO
            END DO
!$OMP END PARALLEL DO
            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie

                del_r_av_z = 0.5*((r_levels_offj(i,j,k)                  &
                                           -r_levels_offj(i,j-1,k))      &
                                + (r_levels_offj(i,j,k-1)                &
                                           -r_levels_offj(i,j-1,k-1)))  

                am_3d(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)        &
                               *del_r_av_z**2                            &
                              /(r_levels(i,j,k)-r_levels(i,j,k-1))       &
                               *diff_coeff_av_2D(i,j)                    &
                               *r_levels_offk(i,j,k-1)**2     
         
              END DO
            END DO

! Inversion of tridiagonal matrix
            k = 1
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/(1.0 - ap_3D(i,j,k))
                  c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                END DO
              END DO
            DO k = 2, metric_level-1
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/((1.0 - ap_3D(i,j,k) - am_3D(i,j,k))&
                                 -c_new(i,j,k-1)*am_3D(i,j,k))
                  c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                END DO
              END DO
            END DO
            k = metric_level
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/((1.0 - am_3D(i,j,k))               &
                                 -c_new(i,j,k-1)*am_3D(i,j,k))
                END DO
              END DO
!=========================================================================
!===== G_Phi =============================================================

            k = 1
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = is,ie
                  diff_coeff_2D(i,j) = diff_coeff_phi(i,j)*delta_phi**2  &
                                       /( phi(j+1) - phi(j))*csphi_off(j)

! This is G_phi stored on v points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k)*(               &
                                      wt_y_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))&
                               + (1.0-wt_z_on2off(i,j,k))*(              &
                                      wt_y_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                  f_av_yz_km= wt_y_on2off(i,j) *field(i,j,k) +           &
                         (1.0-wt_y_on2off(i,j))*field(i,j+1,k)

                  G_phi(i,j,k) = ((field(i,j+1,k) - field(i,j,k) )       &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                - (f_av_yz_k_2D(i,j) - f_av_yz_km)       &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                 *diff_coeff_2D(i,j)                     &
                                 *r_levels_offj(i,j,k)**2
                END DO
              END DO

            DO k = 2, metric_level-1
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = is,ie
! This is G_phi stored on v points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_yz_km = f_av_yz_k_2D(i,j)

                  f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k)*(               &
                                      wt_y_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))&
                               + (1.0-wt_z_on2off(i,j,k))*(              &
                                      wt_y_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                  G_phi(i,j,k) = ((field(i,j+1,k) - field(i,j,k) )       &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                - (f_av_yz_k_2D(i,j) - f_av_yz_km)       &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                 *diff_coeff_2D(i,j)                     &
                                 *r_levels_offj(i,j,k)**2
                END DO
              END DO
            END DO

            k = metric_level
            IF (metric_level == field_dims%k_end  ) THEN
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_phi stored on v points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_yz_km = f_av_yz_k_2D(i,j)

                  f_av_yz_k_2D(i,j) = wt_y_on2off(i,j) *field(i,j,k) +   &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k) 

                  G_phi(i,j,k) = ((field(i,j+1,k) - field(i,j,k) )       &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                - (f_av_yz_k_2D(i,j) - f_av_yz_km)       &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                  *diff_coeff_2D(i,j)                    &
                                  *r_levels_offj(i,j,k)**2
                END DO
              END DO
            ELSE
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_phi stored on v points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_yz_km = f_av_yz_k_2D(i,j)

                  f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k)*(               &
                                      wt_y_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))&
                               + (1.0-wt_z_on2off(i,j,k))*(              &
                                      wt_y_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                  G_phi(i,j,k) = ((field(i,j+1,k) - field(i,j,k) )       &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                - (f_av_yz_k_2D(i,j) - f_av_yz_km   )    &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                  *diff_coeff_2D(i,j)                    &
                                  *r_levels_offj(i,j,k)**2
                END DO
              END DO
            END IF

            IF( at_extremity(PSouth) ) THEN
              DO k = 1, metric_level
                DO i = is,ie
                  G_phi(i,j_sweep_begin-1,k) = 0.0
                END DO              
              END DO 
            END IF
            IF( at_extremity(PNorth) ) THEN
              DO k = 1, metric_level
                DO i = is,ie
                  G_phi(i,j_sweep_end,k) = 0.0 
                END DO              
              END DO
            ENDIF            
!=========================================================================
!===== Compute RHS b =====================================================
            k = 1
            j = j_sweep_begin - 1
            DO i = is,ie
              f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                     &
                                wt_y_on2off(i,j)*field(i,j,k+1-base_k)   &
                         + (1.0-wt_y_on2off(i,j))                        &
                                             *field(i,j+1,k+1-base_k))   &
                         + (1.0-wt_z_on2off(i,j,k))*(                    &
                                wt_y_on2off(i,j)*field(i,j,k-base_k)     &
                         + (1.0-wt_y_on2off(i,j))*field(i,j+1,k-base_k) )
            END DO
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie
                ! RHS term
                ! This is G_phi_eta on w points
                f_av_yz_jm = f_av_yz_j_1D(i)

                f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                   &
                                  wt_y_on2off(i,j)*field(i,j,k+1-base_k) &
                           + (1.0-wt_y_on2off(i,j))                      &
                                               *field(i,j+1,k+1-base_k)) &
                           + (1.0-wt_z_on2off(i,j,k))*(                  &
                                  wt_y_on2off(i,j)*field(i,j,k-base_k)   &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,k-base_k))

                G_phi_eta_km = 0.0

                G_phi_eta_k_2D(i,j) = ((f_av_yz_j_1D(i) - f_av_yz_jm))   &
                                       *0.5*((r_levels_offj(i,j,k+1)     &
                                               -r_levels_offj(i,j-1,k+1))&
                                     + (r_levels_offj(i,j,k)             &
                                               -r_levels_offj(i,j-1,k))) &
                                       *diff_coeff_av_2D(i,j)            &
                                       *r_levels_offk(i,j,k)**2          &
                                       *csphi_on(j)

                b_3d(i,j,k) = field(i,j,k) +                             &
                            ((G_phi(i,j,k) - G_phi(i,j-1,k))             &
                             *recip_del_phi(j)                           &
                           - (G_phi_eta_k_2D(i,j)- G_phi_eta_km))        &
                            * recip_r_squared_delz(i,j,k)*recip_csphi(j) &
                            - beta*ap_3d(i,j,k)                          &
                            *(field(i,j,k+1) - field(i,j,k)) 
              END DO
            END DO 

            DO k = 2, metric_level-1
              j = j_sweep_begin - 1
              DO i = is,ie
                f_av_yz_j_1D(i) = wt_z_on2off(i,j,k)*(                   &
                                  wt_y_on2off(i,j)*field(i,j,k+1-base_k) &
                           + (1.0-wt_y_on2off(i,j))                      &
                                               *field(i,j+1,k+1-base_k)) &
                           + (1.0-wt_z_on2off(i,j,k))*(                  &
                                  wt_y_on2off(i,j)*field(i,j,k-base_k)   &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,k-base_k))
              END DO
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  ! RHS term
                  ! This is G_phi_eta on w points
                  f_av_yz_jm = f_av_yz_j_1D(i)

                  f_av_yz_j_1D(i)= wt_z_on2off(i,j,k)*(                  &
                                   wt_y_on2off(i,j)*field(i,j,k+1-base_k)&
                            + (1.0-wt_y_on2off(i,j))                     &
                                              *field(i,j+1,k+1-base_k))  &
                            + (1.0-wt_z_on2off(i,j,k))*(                 &
                                   wt_y_on2off(i,j)*field(i,j,k-base_k)  &
                            + (1.0-wt_y_on2off(i,j))                     &
                                                  *field(i,j+1,k-base_k))

                  G_phi_eta_km = G_phi_eta_k_2D(i,j)

                  G_phi_eta_k_2D(i,j)= ((f_av_yz_j_1D(i) - f_av_yz_jm )) &
                                        *0.5*((r_levels_offj(i,j,k+1)    &
                                               -r_levels_offj(i,j-1,k+1))&
                                      + (r_levels_offj(i,j,k)            &
                                               -r_levels_offj(i,j-1,k))) &
                                        *diff_coeff_av_2D(i,j)           &
                                        *r_levels_offk(i,j,k)**2         &
                                        *csphi_on(j)

                  b_3D(i,j,k) = field(i,j,k) +                           &
                              ((G_phi(i,j,k) - G_phi(i,j-1,k))           &
                               *recip_del_phi(j)                         &
                             - (G_phi_eta_k_2D(i,j) - G_phi_eta_km))     &
                              * recip_r_squared_delz(i,j,k)              &
                              * recip_csphi(j)                           &
                              - beta*ap_3D(i,j,k)                        &
                              *(field(i,j,k+1) - field(i,j,k))           &
                              + beta*am_3D(i,j,k)                        &
                              *(field(i,j,k)   - field(i,j,k-1)) 
                END DO
              END DO
            END DO

            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie
                ! RHS term
                G_phi_eta_km = G_phi_eta_k_2D(i,j)

                G_phi_eta_k_2D(i,j) = 0.0

                b_3D(i,j,k) = field(i,j,k) +                             &
                            ((G_phi(i,j,k) - G_phi(i,j-1,k))             &
                             *recip_del_phi(j)                           &
                           - (G_phi_eta_k_2D(i,j) - G_phi_eta_km))       &
                            * recip_r_squared_delz(i,j,k)*recip_csphi(j) &
                            + beta*am_3D(i,j,k)                          &
                            *(field(i,j,k)   - field(i,j,k-1))

              END DO
            END DO
!=========================================================================
!===== Trisolve ==========================================================
! Inlining of trisolve
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie
                field(i,j,1) = b_3D(i,j,1)*denom(i,j,1)
              END DO
            END DO
            DO k = 2,metric_level-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  field(i,j,k) = (b_3D(i,j,k)-field(i,j,k-1)             &
                                 *am_3D(i,j,k))*denom(i,j,k)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
            k = metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  field(i,j,k) = (b_3D(i,j,k) - field(i,j,k-1)           &
                                 *am_3D(i,j,k))*denom(i,j,k)
                END DO
              END DO
! Backward solve sweep
            DO k = metric_level-1,1,-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  field(i,j,k) = field(i,j,k)-c_new(i,j,k)*field(i,j,k+1)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
!==========================================================================
!===== Upper level ========================================================
            is = field_dims%i_start
            ie = field_dims%i_end

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
            DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = is,ie
                  G_phi(i,j,k) = ( field(i,j+1,k) - field(i,j,k) )       &
                                /( phi(j+1)    - phi(j)    )             &
                                 * delta_phi**2*diff_coeff_phi(i,j)      &
                                 * csphi_off(j)

                END DO
              END DO
              IF( at_extremity(PSouth) ) THEN
                DO i = is,ie
                  G_phi(i,j_sweep_begin-1,k) = 0.0             
                END DO 
              END IF
              IF( at_extremity(PNorth) ) THEN
                DO i = is,ie
                  G_phi(i,j_sweep_end,k) = 0.0              
                END DO
              ENDIF            
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  field(i,j,k) = field(i,j,k)                            &
                              + (G_phi(i,j,k) - G_phi(i,j-1,k))          &
                                *recip_del_phi(j)*recip_csphi(j)
                END DO
              END DO
            END DO
!$OMP END PARALLEL DO

! -------------------------------------------------------------------------
! --- w fields ---
! --- As for p fields except bottom boundary condition is more complex
          CASE ( fld_type_w ) 

            is = field_dims%i_start
            ie = field_dims%i_end
!===== Compute a ==========================================================
            DO j = j_sweep_begin, j_sweep_end
              recip_del_phi(j) = 1.0/(phi_offj(j) - phi_offj(j-1))
              recip_csphi(j) = 1.0/csphi_on(j)
            END DO

            k = kl
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie

                diff_coeff_av_2D(i,j) = 0.5*(diff_coeff_phi(i,j)         &
                                            +diff_coeff_phi(i,j-1))      &
                                       *delta_phi**2*recip_del_phi(j)**2

                del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)                &
                                          -r_levels_offj(i,j-1,k+1))     &
                                + (r_levels_offj(i,j,k)                  &
                                          -r_levels_offj(i,j-1,k)))   

                ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)        &
                               *del_r_av_z**2                            &
                              /(r_levels(i,j,k+1)-r_levels(i,j,k))       &
                               *diff_coeff_av_2D(i,j)                    &
                               *r_levels_offk(i,j,k+1)**2

              END DO
            END DO

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(del_r_av_z, &
!$OMP& i, j, k)
            DO k = kl+1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  del_r_av_z = 0.5*((r_levels_offj(i,j,k+1)              &
                                            -r_levels_offj(i,j-1,k+1))   &
                                  + (r_levels_offj(i,j,k)                &
                                            -r_levels_offj(i,j-1,k)))   

                  ap_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)      &
                                 *del_r_av_z**2                          &
                                /(r_levels(i,j,k+1)-r_levels(i,j,k))     &
                                 *diff_coeff_av_2D(i,j)                  &
                                 *r_levels_offk(i,j,k+1)**2

                  del_r_av_z = 0.5*((r_levels_offj(i,j,k)                &
                                            -r_levels_offj(i,j-1,k))     &
                                  + (r_levels_offj(i,j,k-1)              &
                                            -r_levels_offj(i,j-1,k-1)))   

                  am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)      &
                                 *del_r_av_z**2                          &
                                /(r_levels(i,j,k)-r_levels(i,j,k-1))     &
                                 *diff_coeff_av_2D(i,j)                  &
                                 *r_levels_offk(i,j,k)**2
                END DO
              END DO
            END DO
!$OMP END PARALLEL DO

            k = metric_level
            DO j = j_sweep_begin, j_sweep_end
              DO i = is,ie

                del_r_av_z = 0.5*((r_levels_offj(i,j,k)                  &
                                          -r_levels_offj(i,j-1,k))       &
                                + (r_levels_offj(i,j,k-1)                &
                                          -r_levels_offj(i,j-1,k-1)))   

                am_3D(i,j,k) = -alpha*recip_r_squared_delz(i,j,k)        &
                               *del_r_av_z**2                            &
                              /(r_levels(i,j,k)-r_levels(i,j,k-1))       &
                               *diff_coeff_av_2D(i,j)                    &
                               *r_levels_offk(i,j,k)**2
              END DO
            END DO
! Inversion of tridiagonal matrix
            k = kl
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/(1.0 - ap_3D(i,j,k))
                  c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                END DO
              END DO
            DO k = kl+1, metric_level-1
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/((1.0 - ap_3D(i,j,k) - am_3D(i,j,k))&
                                 -c_new(i,j,k-1)*am_3D(i,j,k))
                  c_new(i,j,k) = ap_3D(i,j,k)*denom(i,j,k)
                END DO
              END DO
            END DO
            k = metric_level
              DO j =  j_sweep_begin, j_sweep_end
                DO i = is,ie
                  denom(i,j,k) = 1.0/((1.0 - am_3D(i,j,k))               &
                               - c_new(i,j,k-1)*am_3D(i,j,k))
                END DO
              END DO
!=========================================================================
!===== G_phi ============================================================= 
! Bottom boundary
            k = kl
            DO j = j_sweep_begin-1, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
! This is G_phi stored on (i-1/2,j,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                diff_coeff_2D(i,j) = diff_coeff_phi(i,j)*delta_phi**2    &
                                   /(phi(j+1) - phi(j))*csphi_off(j) 


                f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k+1)*(               &
                                    wt_y_on2off(i,j) *field(i,j,k+1)     &
                             + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))  &
                             + (1.0-wt_z_on2off(i,j,k+1))*(              &
                                    wt_y_on2off(i,j) *field(i,j,k)       &
                             + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                f_av_yz_km= wt_y_on2off(i,j) *field(i,j,k) +             &
                       (1.0-wt_y_on2off(i,j))*field(i,j+1,k) 

                del_z =  r_levels_offk(i,j,k+base_k) -                   &
                         r_levels_offk(i,j,k+base_k-1)  


                G_phi(i,j,k) = ((field(i,j+1,k) - field(i,j,k))*del_z    &
                              - (f_av_yz_k_2D(i,j) - f_av_yz_km)         &
                               *(r_levels(i,j+1,k) - r_levels(i,j,k)))   &
                                *diff_coeff_2D(i,j)                      &
                                *r_levels_offj(i,j,k)**2
              END DO
            END DO 

            DO k = kl+1, metric_level-1
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_phi stored on (i-1/2,j,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_yz_km = f_av_yz_k_2D(i,j)

                  f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k+1)*(             &
                                      wt_y_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))&
                               + (1.0-wt_z_on2off(i,j,k+1))*(            &
                                      wt_y_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k))

                  G_phi(i,j,k) = ((field(i,j+1,k) - field(i,j,k))        &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                - (f_av_yz_k_2D(i,j) - f_av_yz_km   )    &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                  *diff_coeff_2D(i,j)                    &
                                  *r_levels_offj(i,j,k)**2
                END DO
              END DO
            END DO

            k = metric_level
            IF ( metric_level == field_dims%k_end ) THEN
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_phi stored on (i-1/2,j,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_yz_km = f_av_yz_k_2D(i,j)

                  f_av_yz_k_2D(i,j) = wt_y_on2off(i,j) *field(i,j,k) +   &
                                 (1.0-wt_y_on2off(i,j))*field(i,j+1,k) 

                  G_phi(i,j,k) = ((field(i,j+1,k) - field(i,j,k))        &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                - (f_av_yz_k_2D(i,j) - f_av_yz_km)       &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                 *diff_coeff_2D(i,j)                     &
                                 *r_levels_offj(i,j,k)**2
                END DO
              END DO
            ELSE
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
! This is G_phi stored on (i-1/2,j,k) points
! 1/d(eta) can be taken outside and canceled with d(eta)/d(r)
                  f_av_yz_km = f_av_yz_k_2D(i,j)

                  f_av_yz_k_2D(i,j) = wt_z_on2off(i,j,k+1)*(             &
                                      wt_y_on2off(i,j) *field(i,j,k+1)   &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))&
                               + (1.0-wt_z_on2off(i,j,k+1))*(            &
                                      wt_y_on2off(i,j) *field(i,j,k)     &
                               + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                  G_phi(i,j,k) = ((field(i,j+1,k) - field(i,j,k))        &
                                 *(r_levels_offk(i,j,k+base_k) -         &
                                   r_levels_offk(i,j,k+base_k-1))        &
                                - (f_av_yz_k_2D(i,j) - f_av_yz_km)       &
                                 *(r_levels(i,j+1,k) - r_levels(i,j,k))) &
                                 *diff_coeff_2D(i,j)                     &
                                 *r_levels_offj(i,j,k)**2
                END DO
              END DO
            END IF
            IF( at_extremity(PSouth) ) THEN
              DO k = kl, metric_level
                DO i = field_dims%i_start,field_dims%i_end
                  G_phi(i,j_sweep_begin-1,k) = 0.0 
                END DO              
              END DO 
            END IF
            IF( at_extremity(PNorth) ) THEN
              DO k = kl, metric_level
                DO i = field_dims%i_start,field_dims%i_end
                  G_phi(i,j_sweep_end,k) = 0.0
                END DO              
              END DO
            ENDIF            
!==========================================================================
!===== Compute RHS b ====================================================== 
            k = kl
            kp = k + 1
            j = j_sweep_begin - 1 
            DO i = field_dims%i_start,field_dims%i_end
              f_av_yz_j_1D(i) = wt_z_on2off(i,j,kp)*(                    &
                                wt_y_on2off(i,j)*field(i,j,k+1)          &
                         + (1.0-wt_y_on2off(i,j))*field(i,j+1,k+1))      &
                         + (1.0-wt_z_on2off(i,j,kp))*(                   &
                                wt_y_on2off(i,j)*field(i,j,k)            &
                         + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) 
            END DO           
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                ! This is G_phi_eta on p points
                f_av_yz_jm = f_av_yz_j_1D(i)

                f_av_yz_j_1D(i) = wt_z_on2off(i,j,kp)*(                  &
                                  wt_y_on2off(i,j)*field(i,j,kp)         &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,kp))     &
                           + (1.0-wt_z_on2off(i,j,kp))*(                 &
                                  wt_y_on2off(i,j)*field(i,j,k)          &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                G_phi_eta_km = 0.0

                G_phi_eta_k_2D(i,j) = ((f_av_yz_j_1D(i) - f_av_yz_jm))   &
                                       *0.5*((r_levels_offj(i,j,kp)      &
                                             -r_levels_offj(i,j-1,kp))   &
                                     + (r_levels_offj(i,j,k)             &
                                             -r_levels_offj(i,j-1,k)))   &
                                       *diff_coeff_av_2D(i,j)            &
                                       *r_levels_offk(i,j,kp)**2         &
                                       *csphi_on(j)

                b_3D(i,j,k) = field(i,j,k) +                             &
                            ((G_phi(i,j,k)   - G_phi(i,j-1,k))           &
                             *recip_del_phi(j)                           &
                           - (G_phi_eta_k_2D(i,j)   - G_phi_eta_km))     &
                            * recip_r_squared_delz(i,j,k)*recip_csphi(j) &
                            - beta*ap_3D(i,j,k)                          &
                            *(field(i,j,k+1) - field(i,j,k))
              END DO
            END DO

            DO k = kl+1, metric_level-1
              kp = k + 1
              j = j_sweep_begin - 1 
              DO i = field_dims%i_start,field_dims%i_end
                f_av_yz_j_1D(i) = wt_z_on2off(i,j,kp)*(                  &
                                  wt_y_on2off(i,j)*field(i,j,kp)         &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,kp))     &
                           + (1.0-wt_z_on2off(i,j,kp))*(                 &
                               +  wt_y_on2off(i,j)*field(i,j,k)          &
                           + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) ) 
              END DO           
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  ! This is G_phi_eta on p points
                  f_av_yz_jm = f_av_yz_j_1D(i)

                  f_av_yz_j_1D(i) = wt_z_on2off(i,j,kp)*(                &
                                   wt_y_on2off(i,j)*field(i,j,kp)        &
                            + (1.0-wt_y_on2off(i,j))*field(i,j+1,kp))    &
                            + (1.0-wt_z_on2off(i,j,kp))*(                &
                                   wt_y_on2off(i,j)*field(i,j,k)         &
                            + (1.0-wt_y_on2off(i,j))*field(i,j+1,k) )

                  G_phi_eta_km = G_phi_eta_k_2D(i,j)

                  G_phi_eta_k_2D(i,j)= ((f_av_yz_j_1D(i) - f_av_yz_jm))  &
                                        *0.5*((r_levels_offj(i,j,kp)     &
                                              -r_levels_offj(i,j-1,kp))  &
                                      + (r_levels_offj(i,j,k)            &
                                              -r_levels_offj(i,j-1,k)))  &
                                        *diff_coeff_av_2D(i,j)           &
                                        *r_levels_offk(i,j,kp)**2        &
                                        *csphi_on(j)

                  b_3D(i,j,k) = field(i,j,k) +                           &
                              ((G_phi(i,j,k)   - G_phi(i,j-1,k))         &
                               *recip_del_phi(j)                         &
                             - (G_phi_eta_k_2D(i,j) - G_phi_eta_km))     &
                              * recip_r_squared_delz(i,j,k)              &
                              * recip_csphi(j)                           &
                              - beta*ap_3D(i,j,k)                        &
                              *(field(i,j,k+1) - field(i,j,k))           &
                              + beta*am_3D(i,j,k)                        &
                              *(field(i,j,k)   - field(i,j,k-1))
                END DO
              END DO
            END DO
            k = metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  G_phi_eta_km = G_phi_eta_k_2D(i,j)
 
                  G_phi_eta_k_2D(i,j) = 0.0

                  b_3D(i,j,k) = field(i,j,k) +                           &
                              ((G_phi(i,j,k)   - G_phi(i,j-1,k))         &
                               *recip_del_phi(j)                         &
                             - (G_phi_eta_k_2D(i,j) - G_phi_eta_km))     &
                              * recip_r_squared_delz(i,j,k)              &
                              * recip_csphi(j)                           &
                              + beta*am_3D(i,j,k)                        &
                              *(field(i,j,k)   - field(i,j,k-1))
                END DO
              END DO
!==========================================================================
!===== Trisolve =========================================================== 
! Inlining of trisolve
            DO j = j_sweep_begin, j_sweep_end
              DO i = field_dims%i_start,field_dims%i_end
                field(i,j,kl) = b_3D(i,j,kl)*denom(i,j,kl)
              END DO
            END DO
            DO k = kl+1,metric_level-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  field(i,j,k) = (b_3D(i,j,k)-field(i,j,k-1)             &
                                 *am_3D(i,j,k))*denom(i,j,k)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
            k = metric_level
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  field(i,j,k) = (b_3D(i,j,k) - field(i,j,k-1)           &
                                 *am_3D(i,j,k))*denom(i,j,k)
                END DO
               END DO
! Backward solve sweep
            DO k = metric_level-1,kl,-1
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j)
              DO j = j_sweep_begin, j_sweep_end
                DO i = field_dims%i_start,field_dims%i_end
                  field(i,j,k) = field(i,j,k)-c_new(i,j,k)*field(i,j,k+1)
                END DO
              END DO
!$OMP END PARALLEL DO
            END DO
!=========================================================================
!===== Upper level ======================================================= 
            is = field_dims%i_start
            ie = field_dims%i_end
!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k)
            DO k = metric_level+1,active_levels
! Grid is flat so diffusion operator simplifies to standard 1-2-1 filter
              DO j = j_sweep_begin-1, j_sweep_end
                DO i = is,ie 
                  G_phi(i,j,k) = (field(i,j+1,k) - field(i,j,k) )        &
                                /(phi(j+1)    - phi(j)    )              &
                                * delta_phi**2*diff_coeff_phi(i,j)       &
                                * csphi_off(j)
                END DO
              END DO
              IF( at_extremity(PSouth) ) THEN
                DO i = is,ie
                  G_phi(i,j_sweep_begin-1,k) = 0.0              
                END DO 
              END IF
              IF( at_extremity(PNorth) ) THEN
                DO i = is,ie
                  G_phi(i,j_sweep_end,k) = 0.0           
                END DO
              ENDIF            
              DO j = j_sweep_begin, j_sweep_end
                DO i = is,ie
                  field(i,j,k) = field(i,j,k)                            &
                              + (G_phi(i,j,k) - G_phi(i,j-1,k))          &
                                *recip_del_phi(j)*recip_csphi(j)
                END DO
              END DO
            END DO
!$OMP END PARALLEL DO
! ------------------------------------------------------------------------
! --- Default ---
          CASE DEFAULT
            CALL ereport('pofil_zlevel',fld_type, 'Invalid field type')
          END SELECT

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                                &
                         field, row_length, in_rows, levels,             &
                         offx, offy, swp_bnds_fld_type, L_vector) 

! ------------------------------------------------------------------------
      END IF ! L_diff
   
      IF (lhook) CALL dr_hook('POFIL_ZLEVEL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE pofil_zlevel
      END MODULE  pofil_zlevel_mod

