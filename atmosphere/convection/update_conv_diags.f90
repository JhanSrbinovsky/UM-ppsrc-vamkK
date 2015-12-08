! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! Update convection diagnostics
!
SUBROUTINE  update_conv_diags(rows, row_length, model_levels, wet_levels,    &
            n_conv_levels, call_number, n_conv_calls,                        & 
          ntml,ntpar, freeze_lev, it_kterm_deep, it_kterm_shall, it_cg_term, &
          cumulus, l_shallow, l_congestus, l_congestus2,it_dp_cfl_limited,   &
          it_md_cfl_limited, it_mid_level,                                   &
          one_over_conv_calls, timestep_conv , wstar,                        &
          exner_theta_levels, z_rho,                                         &

         ! Input (values output by latest convection call)
         it_cape_out, it_conv_rain, it_conv_snow,                            &
         it_precip_dp, it_precip_sh, it_precip_md, it_precip_cg,             &
         ind_cape_reduced, cape_ts_used, it_ind_deep, it_ind_shall,          &
         it_wstar_up, it_mb1, it_mb2,                                        &
         dthbydt, dqbydt, dqclbydt, dqcfbydt, dcflbydt, dcffbydt, dbcfbydt,  &
         dubydt_p, dvbydt_p,                                                 &
         it_up_flux, it_up_flux_half, it_dwn_flux, it_entrain_up,            &
         it_entrain_dwn, it_detrain_up, it_detrain_dwn,                      &
         it_conv_rain_3d, it_conv_snow_3d,                                   &
         it_wqt_flux,it_wthetal_flux,it_wthetav_flux,it_wql_flux,            &
         it_uw_dp, it_vw_dp, it_uw_shall, it_vw_shall, it_uw_mid, it_vw_mid, &
         it_mf_deep, it_mf_congest, it_mf_shall, it_mf_midlev,               &
         it_dt_deep, it_dt_congest, it_dt_shall, it_dt_midlev,               &
         it_dq_deep, it_dq_congest, it_dq_shall, it_dq_midlev,               &
         it_du_deep, it_du_congest, it_du_shall, it_du_midlev,               &
         it_dv_deep, it_dv_congest, it_dv_shall, it_dv_midlev,               &
         ! in/out diagnostics
         conv_rain, conv_snow )


! Modules

USE atm_fields_bounds_mod, ONLY:                                            &
   tdims_s

USE cv_run_mod,  ONLY:                                                      &
   l_mom,                                                                   &
   i_convection_vn,                                              &
   i_convection_vn_4a,                                           &
   i_convection_vn_5a,                                           &
   i_convection_vn_6a

USE cv_stash_flg_mod, ONLY:                                                 &
  flg_up_flx, flg_up_flx_half, flg_dwn_flx,                                 &
  flg_entr_up, flg_entr_dwn, flg_detr_up, flg_detr_dwn,                     &
  flg_conv_rain_3d, flg_conv_snow_3d,                                       &
  flg_uw_dp, flg_vw_dp, flg_uw_shall, flg_vw_shall, flg_uw_mid, flg_vw_mid, &
  flg_wqt_flux, flg_wql_flux, flg_wthetal_flux, flg_wthetav_flux,           &
  flg_mf_deep, flg_mf_congest, flg_mf_shall, flg_mf_midlev,                 &
  flg_dt_deep, flg_dt_congest, flg_dt_shall, flg_dt_midlev,                 &
  flg_dq_deep, flg_dq_congest, flg_dq_shall, flg_dq_midlev,                 &
  flg_du_deep, flg_du_congest, flg_du_shall, flg_du_midlev,                 &
  flg_dv_deep, flg_dv_congest, flg_dv_shall, flg_dv_midlev,                 &
  flg_deep_tops,                                                            &
  l_qcl_incr_cinh, l_qcf_incr_cinh, l_cfl_incr_cinh,                        &
  l_cff_incr_cinh, l_bcf_incr_cinh,                                         &
  l_t_incr_conv, l_q_incr_conv  ,l_qcl_incr_conv, l_qcf_incr_conv,          &
  l_cfl_incr_conv, l_cff_incr_conv, l_bcf_incr_conv

! Convective diagnostic output arrays
USE cv_diagnostic_array_mod, ONLY:                                      &
        precip_deep, precip_shall ,precip_mid ,precip_cong ,cape_out    &
       ,deep_ind ,shallow_ind ,congestus_ind ,congestus_ind2 ,mid_ind   &
       ,ntml_diag ,ntpar_diag ,freeze_diag ,kterm_diag                  &
       ,wstar_up_diag ,wstar_dn_diag ,mb1_diag ,mb2_diag ,wqt_cb        &
       ,wthetal_cb ,wqt_inv ,wthetal_inv ,sh_top ,sh_base ,cg_top       &
       ,cg_base ,cg_term ,cape_ts_diag ,ind_cape_reduced_diag           &
       ,uw_dp, vw_dp ,uw_shall ,vw_shall ,uw_mid ,vw_mid ,up_flux       &
       ,dwn_flux ,entrain_up ,detrain_up ,entrain_dwn ,detrain_dwn      &
       ,up_flux_half ,T_incr_diag_conv ,q_incr_diag_conv                &
       ,qcl_incr_diag_conv ,qcf_incr_diag_conv                          & 
       ,cf_liquid_incr_diag_conv ,cf_frozen_incr_diag_conv              &
       ,bulk_cf_incr_diag_conv,u_incr_diag_conv ,v_incr_diag_conv       &
       ,theta_diag, q_diag                                              &
       ,mf_deep ,mf_congest ,mf_shall ,mf_midlev                        &
       ,dt_deep ,dt_congest ,dt_shall ,dt_midlev                        &
       ,dq_deep ,dq_congest ,dq_shall ,dq_midlev                        &
       ,du_deep ,du_congest ,du_shall ,du_midlev                        &
       ,dv_deep ,dv_congest ,dv_shall ,dv_midlev                        &
       ,wqt_flux_sh ,wql_flux_sh ,wthetal_flux_sh ,wthetav_flux_sh      &
       ,dubydt_pout ,dvbydt_pout ,conv_rain_3d ,conv_snow_3d            &
       ,deep_cfl_limited, mid_cfl_limited, deep_tops

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Updates convection diagnostics after call to glue each substep of
!   convection.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection

! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.

! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) ::  &
  row_length            & ! row length
 ,rows                  & ! rows
 ,model_levels          & ! model levels
 ,wet_levels            & ! wet levels
 ,n_conv_levels         & ! Number of convection levels
 ,call_number           & ! call number for convection sub-steps
 ,n_conv_calls            ! number of convection calls

INTEGER, INTENT(IN) ::           &
  ntml(row_length, rows)         & ! Top level of surface mixed layer
 ,ntpar(row_length, rows)        & ! Top level of initial parcel ascent
 ,freeze_lev(row_length, rows)   & ! freezing level
 ,it_kterm_deep(row_length,rows) & ! lev no for terminating deep convection
 ,it_kterm_shall(row_length,rows) & ! lev no for terminating shallow convection
 ,it_cg_term(row_length,rows)      ! lev no for terminating congestus

LOGICAL, INTENT(IN) ::            &
  cumulus(row_length, rows)       & ! Cumulus present if true (from BL)
 ,l_shallow(row_length, rows)     & ! Logical switch for shallow Cu
 ,l_congestus(row_length, rows)   & ! Logical switch for congestus
 ,l_congestus2(row_length, rows)  & ! congestus in descending air
 ,it_mid_level(row_length, rows)    ! Mid level present in column
 
REAL, INTENT(IN) ::      &
  one_over_conv_calls    & ! 1/n_conv_calls
 ,timestep_conv            ! convection tiemstep (s)

REAL, INTENT(IN) ::          &
  wstar(row_length, rows)    & ! Sub-cloud convective velocity scale (m/s)
 ,exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                     tdims_s%j_start:tdims_s%j_end,  & !
                     tdims_s%k_start:tdims_s%k_end)  &
 ,z_rho(row_length, rows, model_levels)                ! rho levels (m)

REAL, INTENT(IN) ::                 &
  it_cape_out(row_length, rows)     &
 ,it_conv_rain(row_length, rows)    &
 ,it_conv_snow(row_length, rows)    &
 ,it_precip_dp(row_length, rows)    &  ! deep precip
 ,it_precip_sh(row_length, rows)    &  ! shallow precip
 ,it_precip_md(row_length, rows)    &  ! mid-level precip
 ,it_precip_cg(row_length, rows)    &  ! congestus precip
 ,it_wstar_up(row_length, rows)     &
 ,it_mb1(row_length, rows)          &
 ,it_mb2(row_length, rows)          &
 ,ind_cape_reduced(row_length,rows) &  ! indicator of reduced cape timescale
 ,cape_ts_used(row_length, rows)    &  ! Actual cape timescale used for deep(s)
 ,it_dp_cfl_limited(row_length, rows) &! indicator for CFL limited deep conv
 ,it_md_cfl_limited(row_length, rows) &! indicator for CFL limited mid conv
 ,it_ind_deep(row_length, rows)       &! indicator of real deep
 ,it_ind_shall(row_length, rows)       ! indicator of real shallow

REAL, INTENT(IN) ::                                &
  it_up_flux(row_length, rows, model_levels)       &
 ,it_up_flux_half(row_length, rows, model_levels)  &  !up flux on half levs.
 ,it_dwn_flux(row_length, rows, model_levels)      &
 ,it_entrain_up(row_length, rows, model_levels)    &
 ,it_detrain_up(row_length, rows, model_levels)    &
 ,it_entrain_dwn(row_length, rows, model_levels)   &
 ,it_detrain_dwn(row_length, rows, model_levels)   &
 ,it_conv_rain_3d (row_length, rows, wet_levels)   &
 ,it_conv_snow_3d (row_length, rows, wet_levels)   &
 ,it_wqt_flux(row_length, rows, wet_levels)        &
 ,it_wql_flux(row_length, rows, wet_levels)        &
 ,it_wthetal_flux(row_length, rows, wet_levels)    &
 ,it_wthetav_flux(row_length, rows, wet_levels)    &
 ,it_uw_dp(row_length, rows, model_levels)         & 
 ,it_vw_dp(row_length, rows, model_levels)         &
 ,it_uw_shall(row_length, rows, model_levels)      &
 ,it_vw_shall(row_length, rows, model_levels)      &
 ,it_uw_mid(row_length, rows, model_levels)        &
 ,it_vw_mid(row_length, rows, model_levels)        &
 ,it_mf_deep(row_length, rows, wet_levels)                  &
 ,it_mf_congest (row_length, rows, wet_levels)              &
 ,it_mf_shall(row_length, rows, wet_levels)                 &
 ,it_mf_midlev(row_length, rows, wet_levels)                &
 ,it_dt_deep(row_length, rows, wet_levels)                  &
 ,it_dt_congest(row_length, rows, wet_levels)               &
 ,it_dt_shall(row_length, rows, wet_levels)                 &
 ,it_dt_midlev(row_length, rows, wet_levels)                &
 ,it_dq_deep(row_length, rows, wet_levels)                  &
 ,it_dq_congest(row_length, rows, wet_levels)               &
 ,it_dq_shall(row_length, rows, wet_levels)                 &
 ,it_dq_midlev (row_length, rows, wet_levels)               &
 ,it_du_deep (row_length, rows, wet_levels)                 &
 ,it_du_congest(row_length, rows, wet_levels)               &
 ,it_du_shall(row_length, rows, wet_levels)                 &
 ,it_du_midlev (row_length, rows, wet_levels)               &
 ,it_dv_deep (row_length, rows, wet_levels)                 &
 ,it_dv_congest(row_length, rows, wet_levels)               &
 ,it_dv_shall (row_length, rows, wet_levels)                &
 ,it_dv_midlev (row_length, rows, wet_levels)

REAL, INTENT(IN) ::                         &
  dthbydt(row_length, rows, wet_levels)     &
 ,dqbydt(row_length, rows, wet_levels)      &
 ,dqclbydt(row_length, rows, wet_levels)    & ! Q4 Increment qcl
 ,dqcfbydt(row_length, rows, wet_levels)    & ! Q4 Increment qcf
 ,dcflbydt(row_length, rows, wet_levels)    & ! Cloud Increment
 ,dcffbydt(row_length, rows, wet_levels)    & ! Cloud Increment
 ,dbcfbydt(row_length, rows, wet_levels)    & ! Cloud Increment
 ,dubydt_p(row_length, rows, model_levels)  &
 ,dvbydt_p(row_length, rows, model_levels) 

REAL,INTENT(INOUT) ::                     &
  conv_rain(row_length,rows)              & ! convective rainfall
, conv_snow(row_length,rows)                ! convective snowfall


! Local declarations:
INTEGER  ::         &
  i,j,k                ! loop counters

REAL     ::         &
  rsteps               ! 1./(Number of real shallow sub-steps)

! Required by Dr hook 
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('UPDATE_CONV_DIAGS',zhook_in,zhook_handle)

!------------------------------------------------------------------------------
! Diagnostics not on switches
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Single level fields
!------------------------------------------------------------------------------


DO j = 1, rows
  DO i = 1, row_length     
    ntml_diag(i,j)  = ntml_diag(i,j)                                         &
                         + FLOAT(ntml(i,j))*one_over_conv_calls
    ntpar_diag(i,j) = ntpar_diag(i,j)                                        &
                         + FLOAT(ntpar(i,j))*one_over_conv_calls
    freeze_diag(i,j) = freeze_diag(i,j)                                      &
                         + FLOAT(freeze_lev(i,j)) *one_over_conv_calls
  END DO
END DO

DO j = 1, rows
  DO i = 1, row_length

    wstar_dn_diag(i,j) = wstar_dn_diag(i,j) + wstar(i,j)*one_over_conv_calls

  END DO
END DO

! dependent on type of convection

DO j = 1, rows
  DO i = 1, row_length

    IF (it_mid_level(i,j)) THEN
      mid_ind(i,j) = mid_ind(i,j)+one_over_conv_calls
    END IF

    ! deep_ind - now set dependent on whether deep convection really occurs
    ! in the deep convection scheme. Failed deep attempts will be set to zero
    ! and not counted.
    deep_ind(i,j) = deep_ind(i,j) + it_ind_deep(i,j)*one_over_conv_calls

    ! Dependent on whether shallow really occurs
    shallow_ind(i,j) = shallow_ind(i,j) + it_ind_shall(i,j)*one_over_conv_calls

    ! Only work out top and bottom for points where real shallow
    ! convection occurred
    IF (it_ind_shall(i,j) == 1.0) THEN
      ! Use actual termination level rather than ntpar (5A & 6A)
      ! (4A value set to ntpar)
      k=it_kterm_shall(i,j)+1
      sh_top(i,j) =sh_top(i,j) + z_rho(i,j,k)

      k=ntml(i,j)+1
      sh_base(i,j)=sh_base(i,j)+ z_rho(i,j,k)

    END IF
  END DO
END DO

! Final sub-step
IF (call_number == n_conv_calls) THEN
  ! correctly scale heights by number of sub-steps with real shallow

  DO j = 1, rows
    DO i = 1, row_length
      IF (shallow_ind(i,j) > 0.0) THEN   ! there are shallow sub-steps
        rsteps=1./(shallow_ind(i,j)*n_conv_calls)
        sh_top(i,j) =sh_top(i,j)*rsteps
        sh_base(i,j)=sh_base(i,j)*rsteps     
      END IF  
    END DO
  END DO
END IF



! Congestus type only defined for 5A & 6A 
IF (i_convection_vn == i_convection_vn_5a .OR.    &
    i_convection_vn == i_convection_vn_6a ) THEN
DO j = 1, rows
  DO i = 1, row_length
    IF (l_congestus(i,j)) THEN
      congestus_ind(i,j) = congestus_ind(i,j) + one_over_conv_calls

      k=it_cg_term(i,j)+1
      cg_top(i,j) =cg_top(i,j) + z_rho(i,j,k) *one_over_conv_calls

      k=ntml(i,j)+1
      cg_base(i,j)=cg_base(i,j)+ z_rho(i,j,k) *one_over_conv_calls

    END IF
    IF (l_congestus2(i,j)) THEN
      congestus_ind2(i,j) = congestus_ind2(i,j)  +one_over_conv_calls
    END IF
  END DO
END DO
END IF


DO j = 1, rows
  DO i = 1, row_length

    conv_rain(i,j) = conv_rain(i,j) + it_conv_rain(i,j) * one_over_conv_calls

    conv_snow(i,j) = conv_snow(i,j) + it_conv_snow(i,j) * one_over_conv_calls

    precip_deep(i,j) = precip_deep(i,j) + it_precip_dp(i,j)                  &
                                                  * one_over_conv_calls
    precip_shall(i,j)= precip_shall(i,j) + it_precip_sh(i,j)                 &
                                                  * one_over_conv_calls
    precip_mid(i,j)  = precip_mid(i,j) + it_precip_md(i,j)                    &
                                                  * one_over_conv_calls
    cape_out(i,j)   = cape_out(i,j) +it_cape_out(i,j)*one_over_conv_calls

    kterm_diag(i,j) = kterm_diag(i,j)+                        &
                        FLOAT(it_kterm_deep(i,j)) *one_over_conv_calls

    ind_cape_reduced_diag(i,j) = ind_cape_reduced_diag(i,j) + &
                         ind_cape_reduced(i,j) *one_over_conv_calls
    cape_ts_diag(i,j) = cape_ts_diag(i,j) +                   &
                         cape_ts_used(i,j) *one_over_conv_calls

    deep_cfl_limited(i,j) = deep_cfl_limited(i,j) +                         &
                                 it_dp_cfl_limited(i,j) *one_over_conv_calls
    mid_cfl_limited(i,j)  = mid_cfl_limited(i,j) +                          &
                                 it_md_cfl_limited(i,j) *one_over_conv_calls

  END DO
END DO

! 5A & 6A only
IF (i_convection_vn == i_convection_vn_5a .OR.    &
    i_convection_vn == i_convection_vn_6a ) THEN
DO j = 1, rows
  DO i = 1, row_length
    wstar_up_diag(i,j) = wstar_up_diag(i,j)                                 &
                           +it_wstar_up(i,j)*one_over_conv_calls
    mb1_diag(i,j) = mb1_diag(i,j)                                           &
                       +it_mb1(i,j)*one_over_conv_calls
    mb2_diag(i,j) = mb2_diag(i,j)                                           &
                           +it_mb2(i,j)*one_over_conv_calls
    cg_term(i,j) = cg_term(i,j)+                                            &
                      FLOAT(it_cg_term(i,j)) *one_over_conv_calls

    precip_cong(i,j) = precip_cong(i,j) + it_precip_cg(i,j)                 &
                                                  * one_over_conv_calls
  END DO
END DO
END IF
! ------------------------------------------------------------------------------
! Diagnostics on switches
! ------------------------------------------------------------------------------

IF (flg_conv_rain_3d) THEN
  DO k=1, wet_levels
    DO j=1, rows
      DO i=1, row_length
        conv_rain_3d(i,j,k) = conv_rain_3d(i,j,k)                             &
                               + it_conv_rain_3d(i,j,k) * one_over_conv_calls
      END DO      ! i (row_length)
    END DO      ! j (rows)
  END DO      ! k (wet_levels)
END IF

IF (flg_conv_snow_3d) THEN
  DO k=1, wet_levels
    DO j=1, rows
      DO i=1, row_length
        conv_snow_3d(i,j,k) = conv_snow_3d(i,j,k)                             &
                               + it_conv_snow_3d(i,j,k) * one_over_conv_calls
      END DO      ! i (row_length)
    END DO      ! j (rows)
  END DO      ! k (wet_levels)
END IF

IF (flg_deep_tops) THEN

  DO j = 1, rows
    DO i = 1, row_length
      ! Deep convection
      ! Altered condition as can get failed deep cases 
      IF (it_ind_deep(i,j) == 1.0) THEN
        k=it_kterm_deep(i,j)
        IF (k > 0) THEN  ! in case still get a zero value
          deep_tops(i,j,k) = deep_tops(i,j,k) + one_over_conv_calls 
        END IF
      END IF
    END DO
  END DO

END IF

IF (flg_mf_deep) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        mf_deep(i,j,k) = mf_deep(i,j,k)+it_mf_deep(i,j,k)*one_over_conv_calls 
      END DO
    END DO
  END DO
END IF
IF (flg_mf_congest) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        mf_congest(i,j,k) = mf_congest(i,j,k)                                 &
                                    +it_mf_congest(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_mf_shall) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        mf_shall(i,j,k) = mf_shall(i,j,k)                                     &
                                      +it_mf_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_mf_midlev) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        mf_midlev(i,j,k) = mf_midlev(i,j,k)                                   &
                                      +it_mf_midlev(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_dt_deep) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dt_deep(i,j,k) = dt_deep(i,j,k) + it_dt_deep(i,j,k)*one_over_conv_calls 
      END DO
    END DO
  END DO
END IF
IF (flg_dt_congest) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dt_congest(i,j,k) = dt_congest(i,j,k)                                 &
                                    +it_dt_congest(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_dt_shall) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dt_shall(i,j,k) = dt_shall(i,j,k)                                     &
                                      +it_dt_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_dt_midlev) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dt_midlev(i,j,k) = dt_midlev(i,j,k)                                   &
                                     +it_dt_midlev(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_dq_deep) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dq_deep(i,j,k) = dq_deep(i,j,k) +it_dq_deep(i,j,k)*one_over_conv_calls 
      END DO
    END DO
  END DO
END IF
IF (flg_dq_congest) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dq_congest(i,j,k) = dq_congest(i,j,k)                                 &
                                    +it_dq_congest(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_dq_shall) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dq_shall(i,j,k) = dq_shall(i,j,k)                                     &
                                      +it_dq_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_dq_midlev) THEN
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dq_midlev(i,j,k) = dq_midlev(i,j,k)                                   &
                                     +it_dq_midlev(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_du_deep) THEN
  DO k=1,n_conv_levels+1
    DO j = 1, rows
      DO i = 1, row_length
        du_deep(i,j,k) =  du_deep(i,j,k) +it_du_deep(i,j,k)*one_over_conv_calls 
      END DO
    END DO
  END DO
END IF
IF (flg_du_congest) THEN
  DO k=1,n_conv_levels+1
    DO j = 1, rows
      DO i = 1, row_length
        du_congest(i,j,k) = du_congest(i,j,k)                                 &
                                    +it_du_congest(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_du_shall) THEN
  DO k=1,wet_levels
    DO j = 1, rows
      DO i = 1, row_length
        du_shall(i,j,k) = du_shall(i,j,k)                                     &
                                      +it_du_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_du_midlev) THEN
  DO k=1,wet_levels
    DO j = 1, rows
      DO i = 1, row_length
        du_midlev(i,j,k) = du_midlev(i,j,k)                                   &
                                      +it_du_midlev(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
IF (flg_dv_deep) THEN
  DO k=1,wet_levels
    DO j = 1, rows
      DO i = 1, row_length
        dv_deep(i,j,k) = dv_deep(i,j,k) +it_dv_deep(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
END IF
    IF (flg_dv_congest) THEN
      DO k=1,wet_levels
        DO j = 1, rows
          DO i = 1, row_length
            dv_congest(i,j,k) = dv_congest(i,j,k)                 &
                     +it_dv_congest(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF
    IF (flg_dv_shall) THEN
      DO k=1,wet_levels
        DO j = 1, rows
          DO i = 1, row_length
            dv_shall(i,j,k) = dv_shall(i,j,k)                     &
                      +it_dv_shall(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF
    IF (flg_dv_midlev) THEN
      DO k=1,wet_levels
        DO j = 1, rows
          DO i = 1, row_length
            dv_midlev(i,j,k) = dv_midlev(i,j,k)                   &
                    +it_dv_midlev(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF



! Update T, q diagnostic incerements

IF ( l_t_incr_conv ) THEN
  DO k=1,n_conv_levels
    DO j=1,rows
      DO i=1,row_length
        t_incr_diag_conv(i,j,k) = t_incr_diag_conv(i,j,k)                   &
                                       + dthbydt(i,j,k) * timestep_conv       &
                                         * exner_theta_levels(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
END IF                   ! on STASHflag

    IF ( l_q_incr_conv ) THEN
      DO k=1,n_conv_levels
        DO j=1,rows
          DO i=1,row_length
            q_incr_diag_conv(i,j,k) = q_incr_diag_conv(i,j,k)   &
                          + dqbydt(i,j,k) * timestep_conv
          END DO ! i
        END DO ! j
      END DO ! k
    END IF                   ! on STASHflag

    IF ( l_qcl_incr_conv ) THEN
      DO k=1,n_conv_levels
        DO j=1,rows
          DO i=1,row_length
            qcl_incr_diag_conv(i,j,k) =                          &
                     qcl_incr_diag_conv(i,j,k) +                 &
                         dqclbydt(i,j,k) * timestep_conv
          END DO ! i
        END DO ! j
      END DO
    END IF

    IF ( l_qcf_incr_conv ) THEN
      DO k=1,n_conv_levels
        DO j=1,rows
          DO i=1,row_length
            qcf_incr_diag_conv(i,j,k) =                          &
                     qcf_incr_diag_conv(i,j,k) +                 &
                            dqcfbydt(i,j,k) * timestep_conv
          END DO ! i
        END DO ! j
      END DO
    END IF

    IF ( l_cfl_incr_conv ) THEN
      DO k=1,n_conv_levels
        DO j=1,rows
          DO i=1,row_length
            cf_liquid_incr_diag_conv(i,j,k) =                    &
                          cf_liquid_incr_diag_conv(i,j,k) +      &
                          dcflbydt(i,j,k) * timestep_conv
          END DO ! i
        END DO ! j
      END DO
    END IF

    IF ( l_cff_incr_conv ) THEN
      DO k=1,n_conv_levels
        DO j=1,rows
          DO i=1,row_length
            cf_frozen_incr_diag_conv(i,j,k) =                    &
                          cf_frozen_incr_diag_conv(i,j,k) +      &
                          dcffbydt(i,j,k) * timestep_conv
          END DO ! i
        END DO ! j
      END DO ! k
    END IF

    IF ( l_bcf_incr_conv ) THEN
      DO k=1,n_conv_levels
        DO j=1,rows
          DO i=1,row_length
            bulk_cf_incr_diag_conv(i,j,k) =                      &
                       bulk_cf_incr_diag_conv(i,j,k)   +         &
                         (dbcfbydt(i,j,k) * timestep_conv)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF

    IF (l_mom) THEN
      ! Require fields whether diagnostics required or not now as 
      ! dubydt_pout and dvbydt_pout hold total increments to U and V from
      ! convection and are used to update r_u and r_v in atmos_physics2
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            dubydt_pout(i,j,k) =dubydt_pout(i,j,k)                &
                                  +dubydt_p(i,j,k)*one_over_conv_calls
            dvbydt_pout(i,j,k) =dvbydt_pout(i,j,k)                &
                                  +dvbydt_p(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_up_flx) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            up_flux(i,j,k)=up_flux(i,j,k)+it_up_flux(i,j,k)*      &
                                       one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_up_flx_half) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            up_flux_half(i,j,k)=up_flux_half(i,j,k)               &
                +it_up_flux_half(i,j,k)* one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_dwn_flx) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            dwn_flux(i,j,k)=dwn_flux(i,j,k)+it_dwn_flux(i,j,k)*   &
                                       one_over_conv_calls
          END DO
        END DO
      END DO
    END IF
    IF (flg_entr_up) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            entrain_up(i,j,k)=entrain_up(i,j,k)+                  &
                   it_entrain_up(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_entr_dwn) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            entrain_dwn(i,j,k)=entrain_dwn(i,j,k)+                &
                   it_entrain_dwn(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_detr_up) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            detrain_up(i,j,k)=detrain_up(i,j,k)+                  &
                   it_detrain_up(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_detr_dwn) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            detrain_dwn(i,j,k)=detrain_dwn(i,j,k)+                &
                   it_detrain_dwn(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_uw_dp) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            uw_dp(i,j,k)=uw_dp(i,j,k)+                            &
                   it_uw_dp(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_vw_dp) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            vw_dp(i,j,k)=vw_dp(i,j,k)+                            &
                   it_vw_dp(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_uw_shall) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            uw_shall(i,j,k)=uw_shall(i,j,k)+                      &
                   it_uw_shall(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_vw_shall) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            vw_shall(i,j,k)=vw_shall(i,j,k)+                      &
                   it_vw_shall(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_uw_mid) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            uw_mid(i,j,k)=uw_mid(i,j,k)+                          &
                               it_uw_mid(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_vw_mid) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            vw_mid(i,j,k)=vw_mid(i,j,k)+                      &
                               it_vw_mid(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

    IF (flg_wqt_flux) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            wqt_flux_sh(i,j,k)=wqt_flux_sh(i,j,k)+                &
                   it_wqt_flux(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
        DO j = 1, rows
          DO i = 1, row_length
            k=ntml(i,j)+1
            wqt_cb(i,j)=wqt_cb(i,j)+                              &
                   it_wqt_flux(i,j,k)*one_over_conv_calls
          END DO
        END DO
        DO j = 1, rows
          DO i = 1, row_length
            k=ntpar(i,j)+1
            wqt_inv(i,j)=wqt_inv(i,j)+                            &
                   it_wqt_flux(i,j,k)*one_over_conv_calls
          END DO
        END DO
    END IF
    IF (flg_wql_flux) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            wql_flux_sh(i,j,k)=wql_flux_sh(i,j,k)+                &
                   it_wql_flux(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF
    IF (flg_wthetal_flux) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            wthetal_flux_sh(i,j,k)=wthetal_flux_sh(i,j,k)+        &
                   it_wthetal_flux(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
        DO j = 1, rows
          DO i = 1, row_length
            k=ntml(i,j)+1
            wthetal_cb(i,j)=wthetal_cb(i,j)+                      &
                   it_wthetal_flux(i,j,k)*one_over_conv_calls
          END DO
        END DO
        DO j = 1, rows
          DO i = 1, row_length
            k=ntpar(i,j)+1
            wthetal_inv(i,j)=wthetal_inv(i,j)+                    &
                   it_wthetal_flux(i,j,k)*one_over_conv_calls
          END DO
        END DO
    END IF
    IF (flg_wthetav_flux) THEN
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            wthetav_flux_sh(i,j,k)=wthetav_flux_sh(i,j,k)+         &
                   it_wthetav_flux(i,j,k)*one_over_conv_calls
          END DO
        END DO
      END DO
    END IF

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('UPDATE_CONV_DIAGS',zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
END SUBROUTINE
