! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Restore convection variables.
!
! Subroutine Interface:
SUBROUTINE restore_conv                                                       &
! Input data
  ( resdump, nvars, row_length, rows, model_levels, wet_levels, tr_levels     &
  , tr_vars, tr_ukca, offx, offy                                              &
! Output data
  , theta_conv, q_conv, qcl_conv, qcf_conv, cf_liquid_conv, cf_frozen_conv    &
  , bulk_cf_conv, theta_inc, q_inc, qcl_inc, qcf_inc, cf_liquid_inc           &
  , cf_frozen_inc, bulk_cf_inc, r_u, r_v, aerosol, dust_div1, dust_div2       &
  , dust_div3, dust_div4, dust_div5, dust_div6, so2, so4_aitken, so4_accu     &
  , so4_diss, dms, nh3, soot_new, soot_agd, soot_cld, ocff_new, ocff_agd      &
  , ocff_cld, nitr_acc, nitr_diss, co2, free_tracers, ukca_tracers            &
  , ozone_tracer, cclwp, conv_rain, conv_snow )

  USE dust_parameters_mod, ONLY: l_twobin_dust

  IMPLICIT NONE

!
! Description: To restore convection variables after the call to
!              convection in the case where convection is called
!              to get the diagnostics out only.
!
! Method:

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!

! Code description:
!   FORTRAN 77 + common extensions also in fortran 90.
!   This code is written to UM programming standards version 7.4

!     INCLUDED COMDECKS

!     Inputs

  INTEGER, INTENT(In) ::                                            &
    row_length                                                      &
  , rows                                                            &
  , model_levels                                                    &
  , wet_levels                                                      &
  , tr_levels                                                       &
  , tr_vars                                                         &
  , tr_ukca                                                         &
  , offx                                                            &
  , offy                                                            &
  , nvars

   REAL, INTENT(Out) ::                                             &
    theta_conv(row_length, rows, model_levels)                      &
  , q_conv(row_length, rows, wet_levels)                            &
  , qcl_conv(row_length, rows, wet_levels)                          &
  , qcf_conv(row_length, rows, wet_levels)                          &
  , cf_liquid_conv(row_length, rows, wet_levels)                    &
  , cf_frozen_conv(row_length, rows, wet_levels)                    &
  , bulk_cf_conv(row_length, rows, wet_levels)                      &
  , theta_inc(row_length, rows, model_levels)                       &
  , q_inc(row_length, rows, wet_levels)                             &
  , qcl_inc(row_length, rows, wet_levels)                           &
  , qcf_inc(row_length, rows, wet_levels)                           &
  , cf_liquid_inc(row_length, rows, wet_levels)                     &
  , cf_frozen_inc(row_length, rows, wet_levels)                     &
  , bulk_cf_inc(row_length, rows, wet_levels)                       &
  , r_u(row_length, rows, model_levels)                             &
  , r_v(row_length, rows, model_levels)                             &
  , aerosol(row_length, rows, model_levels)                         &
  , dust_div1(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       !dust mmr in div1
  , dust_div2(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       !dust mmr in div2
  , dust_div3(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       !dust mmr in div3
  , dust_div4(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       !dust mmr in div4
  , dust_div5(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       !dust mmr in div5
  , dust_div6(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       !dust mmr in div6
  , so2(row_length, rows, model_levels)                             &
  , so4_aitken(row_length, rows, model_levels)                      &
  , so4_accu(row_length, rows, model_levels)                        &
  , so4_diss(row_length, rows, model_levels)                        &
  , dms(row_length, rows, model_levels)                             &
  , nh3(row_length, rows, model_levels)                             &
  , soot_new(row_length, rows, model_levels)                        &
  , soot_agd(row_length, rows, model_levels)                        &
  , soot_cld(row_length, rows, model_levels)                        &
  , ocff_new(row_length, rows, model_levels)                        &
  , ocff_agd(row_length, rows, model_levels)                        &
  , ocff_cld(row_length, rows, model_levels)                        &
  , nitr_acc(row_length, rows, model_levels)                        &
  , nitr_diss(row_length, rows, model_levels)                       &
  , co2(row_length, rows, model_levels)                             &
  , free_tracers(row_length, rows, tr_levels, tr_vars)              &
  , ukca_tracers(row_length, rows, tr_levels, tr_ukca)              &
  , cclwp(row_length, rows)                                         &
  , conv_rain(row_length, rows)                                     &
  , conv_snow(row_length, rows)                                     &
  , ozone_tracer(1-offx:row_length+offx, 1-offy:rows+offy,          &
           model_levels)
                       ! Cariolle ozone tracer

  REAL, INTENT(In) ::                                               &
    resdump(row_length, rows, nvars)

! Local Variables

  INTEGER ::                                                        &
    i, j, k, l, kcount

!----------------------------------------------------------------------

  DO i=1, row_length
    DO j=1, rows

      DO k=1, model_levels
        theta_conv(i,j,k) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        q_conv(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        qcl_conv(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        qcf_conv(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        cf_liquid_conv(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        cf_frozen_conv(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        bulk_cf_conv(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        theta_inc(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        q_inc(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        qcl_inc(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        qcf_inc(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        cf_liquid_inc(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        cf_frozen_inc(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        bulk_cf_inc(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        R_u(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        R_v(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        aerosol(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        dust_div1(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        dust_div2(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      IF (.NOT.l_twobin_dust) THEN
        kcount = k
        DO k=kcount, kcount + model_levels - 1
          dust_div3(i,j,k-kcount+1) = resdump(i,j,k)
        END DO

        kcount = k
        DO k=kcount, kcount + model_levels - 1
          dust_div4(i,j,k-kcount+1) = resdump(i,j,k)
        END DO

        kcount = k
        DO k=kcount, kcount + model_levels - 1
          dust_div5(i,j,k-kcount+1) = resdump(i,j,k)
        END DO

        kcount = k
        DO k=kcount, kcount + model_levels - 1
          dust_div6(i,j,k-kcount+1) = resdump(i,j,k)
        END DO
      END IF

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        so2(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        so4_aitken(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        so4_accu(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        so4_diss(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        dms(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        nh3(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        soot_new(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        soot_agd(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        soot_cld(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        ocff_new(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        ocff_agd(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        ocff_cld(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        nitr_acc(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        nitr_diss(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        co2(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        ozone_tracer(i,j,k-kcount+1) = resdump(i,j,k)
      END DO

      kcount = k
      cclwp(i,j) = resdump(i,j,k)

      k = k+1
      kcount = kcount+1
      conv_snow(i,j) = resdump(i,j,k)
      
      k = k+1
      kcount = kcount+1
      conv_rain(i,j) = resdump(i,j,k)
      
      kcount = kcount+1

      DO l=1, tr_vars
        DO k=kcount, kcount + tr_levels - 1
          free_tracers(i,j,k-kcount+1,l) = resdump(i,j,k)
        END DO
        kcount = k
      END DO

      DO l = 1, tr_ukca
        DO k = kcount, kcount + tr_levels - 1
          ukca_tracers(i,j,k-kcount+1,l) = resdump(i,j,k)
        END DO
        kcount = k
      END DO


    END DO
  END DO

  RETURN

END SUBROUTINE restore_conv

!----------------------------------------------------------------------
