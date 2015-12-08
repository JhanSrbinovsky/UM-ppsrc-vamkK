! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Save convection variables.
!
! Subroutine Interface:
! Input data

SUBROUTINE save_conv                                                          &
  ( row_length, rows, model_levels, wet_levels, tr_levels, tr_vars, tr_ukca   &
  , offx, offy, theta_conv, q_conv, qcl_conv, qcf_conv, cf_liquid_conv        &
  , cf_frozen_conv, bulk_cf_conv, theta_inc, q_inc, qcl_inc, qcf_inc          &
  , cf_liquid_inc, cf_frozen_inc, bulk_cf_inc, r_u, r_v, aerosol, dust_div1   &
  , dust_div2, dust_div3, dust_div4, dust_div5, dust_div6, so2, so4_aitken    &
  , so4_accu, so4_diss, dms, nh3, soot_new, soot_agd, soot_cld, ocff_new      &
  , ocff_agd, ocff_cld, nitr_acc, nitr_diss, co2, free_tracers, ukca_tracers  &
  , ozone_tracer, cclwp, conv_rain, conv_snow, nvars                          &
! Output data
  , resdump )

  USE dust_parameters_mod, ONLY: l_twobin_dust

  IMPLICIT NONE

!
! Description:  To save convection variables before the call to
!               convection in the case where convection is called
!               to get the diagnostics out only.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code description:
!    Programming standards :
!      Fortran 90, Written to UM coding standards
!      as specified in UMDP 3, vn8.2

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
  , nvars                                                           &
  , offx                                                            &
  , offy

   REAL, INTENT(In) ::                                              &
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
                       ! dust mmr in div1
  , dust_div2(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       ! dust mmr in div2
  , dust_div3(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       ! dust mmr in div3
  , dust_div4(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       ! dust mmr in div4
  , dust_div5(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       ! dust mmr in div5
  , dust_div6(1-offx:row_length+offx, 1-offy:rows+offy,             &
         model_levels)                                              &
                       ! dust mmr in div6
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
                       ! Ozone tracer for cariolle scheme

  REAL, INTENT(Out) ::                                              &
    resdump(row_length, rows, nvars)

! Local Variables

  INTEGER ::                                                        &
    i,j,k,l, kcount

!----------------------------------------------------------------------

  resdump(:,:,:) = 0.0

  DO i=1, row_length
    DO j=1, rows

      DO k=1, model_levels
        resdump(i,j,k) = theta_conv(i,j,k)
      END DO

      kcount=k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = q_conv(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = qcl_conv(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = qcf_conv(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = cf_liquid_conv(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = cf_frozen_conv(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = bulk_cf_conv(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = theta_inc(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = q_inc(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = qcl_inc(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = qcf_inc(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = cf_liquid_inc(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = cf_frozen_inc(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + wet_levels - 1
        resdump(i,j,k) = bulk_cf_inc(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = R_u(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = R_v(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = aerosol(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = dust_div1(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = dust_div2(i,j,k)
      END DO

      IF (.NOT. l_twobin_dust) THEN
        kcount = k
        DO k=kcount, kcount + model_levels - 1
          resdump(i,j,k) = dust_div3(i,j,k)
        END DO

        kcount = k
        DO k=kcount, kcount + model_levels - 1
          resdump(i,j,k) = dust_div4(i,j,k)
        END DO

        kcount = k
        DO k=kcount, kcount + model_levels - 1
          resdump(i,j,k) = dust_div5(i,j,k)
        END DO

        kcount = k
        DO k=kcount, kcount + model_levels - 1
          resdump(i,j,k) = dust_div6(i,j,k)
        END DO
      END IF ! l_twobin_dust

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = so2(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = so4_aitken(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = so4_accu(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = so4_diss(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = dms(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = nh3(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = soot_new(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = soot_agd(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = soot_cld(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = ocff_new(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = ocff_agd(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = ocff_cld(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = nitr_acc(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = nitr_diss(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = co2(i,j,k)
      END DO

      kcount = k
      DO k=kcount, kcount + model_levels - 1
        resdump(i,j,k) = ozone_tracer(i,j,k)
      END DO

      kcount = k
      resdump(i,j,k) = cclwp(i,j)
      
      k = k+1
      kcount = kcount+1
      resdump(i,j,k) = conv_snow(i,j)
      
      k = k+1
      kcount = kcount+1
      resdump(i,j,k) = conv_rain(i,j)

      kcount = kcount+1
      
      DO l=1, tr_vars
        DO k=kcount, kcount + tr_levels - 1
          resdump(i,j,k) = free_tracers(i,j,k,l)
        END DO
        kcount = k
      END DO

      DO l = 1, tr_ukca
        DO k = kcount, kcount + tr_levels - 1
          resdump(i,j,k) = ukca_tracers(i,j,k,l)
        END DO
        kcount = k
      END DO


    END DO
  END DO

  RETURN

END SUBROUTINE save_conv

!----------------------------------------------------------------------
