! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine qt_bal_cld
SUBROUTINE qt_bal_cld                                                   &
    (p_star,p_theta_levels,p,                                           &
     theta,exner_theta_levels,                                          &
     q,qcl,qcf,qcf2,                                                    &
     rhcpt, rhc_row_length, rhc_rows, bl_levels,                        &
     delta_lambda, delta_phi,                                           &
     fv_cos_theta_latitude,                                             &
     l_cld_area, l_acf_cusack, l_acf_brooks,                            &
     l_mcr_qcf2, l_mixing_ratio, ntml, cumulus,                         &
     area_cloud_fraction,  bulk_cloud_fraction,                         &
     cloud_fraction_liquid,  cloud_fraction_frozen,                     &
     mype)

! Purpose:
!        reset q, t and the cloud fields to be consistent at the
!        end of the timestep
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

  USE atmos_constants_mod, ONLY: cp
  USE atm_fields_bounds_mod, ONLY: pdims, pdims_s,          &
                                   tdims, tdims_s, tdims_l, &
                                   qdims, qdims_l
  USE water_constants_mod, ONLY: lc
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

  INTEGER ::                                                            &
    mype,               &  ! My processor number
    rhc_row_length,     &  ! Array size for RHcrit array
    rhc_rows,           &  ! Array size for RHcrit array
    bl_levels


  LOGICAL ::                                                            &
    l_cld_area,                                                         &
                           ! true if using area cloud fraction (ACF)
    l_acf_cusack,                                                       &
                           ! ... and selected Cusack and PC2 off
    l_acf_brooks,                                                       &
                           ! ... and selected Brooks
    l_mcr_qcf2,                                                         &
                           ! true if second cloud ice variable in use
    l_mixing_ratio         ! true if using mixing ratio formulation

  REAL ::                                                               &
    p(                     pdims_s%i_start:pdims_s%i_end,               & 
                           pdims_s%j_start:pdims_s%j_end,               &
                           pdims_s%k_start:pdims_s%k_end+1),            &
    p_theta_levels(        pdims_s%i_start:pdims_s%i_end,               & 
                           pdims_s%j_start:pdims_s%j_end,               &
                           pdims_s%k_start:pdims_s%k_end),              &
    p_star(                  pdims%i_start:pdims%i_end,                 & 
                             pdims%j_start:pdims%j_end),                &
    theta(                 tdims_s%i_start:tdims_s%i_end,               & 
                           tdims_s%j_start:tdims_s%j_end,               &
                           tdims_s%k_start:tdims_s%k_end),              &
    exner_theta_levels(    tdims_s%i_start:tdims_s%i_end,               & 
                           tdims_s%j_start:tdims_s%j_end,               &
                           tdims_s%k_start:tdims_s%k_end),              &
    q(                     qdims_l%i_start:qdims_l%i_end,               & 
                           qdims_l%j_start:qdims_l%j_end,               &
                           qdims_l%k_start:qdims_l%k_end),              &
    qcl(                   qdims_l%i_start:qdims_l%i_end,               & 
                           qdims_l%j_start:qdims_l%j_end,               &
                           qdims_l%k_start:qdims_l%k_end),              &
    qcf(                   qdims_l%i_start:qdims_l%i_end,               & 
                           qdims_l%j_start:qdims_l%j_end,               &
                           qdims_l%k_start:qdims_l%k_end),              &
    qcf2(                  qdims_l%i_start:qdims_l%i_end,               & 
                           qdims_l%j_start:qdims_l%j_end,               &
                           qdims_l%k_start:qdims_l%k_end),              &
    rhcpt(rhc_row_length, rhc_rows, qdims%k_end),                       &
! coordinate arrays
    delta_lambda,                                                       &
    delta_phi,                                                          &
! trig arrays
    fv_cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,               & 
                           tdims_s%j_start:tdims_s%j_end)

! Diagnostic variables

  REAL ::                                                               &
    area_cloud_fraction(     qdims%i_start:qdims%i_end,                 & 
                             qdims%j_start:qdims%j_end,                 &
                                         1:qdims%k_end),                &
    bulk_cloud_fraction(   qdims_l%i_start:qdims_l%i_end,               & 
                           qdims_l%j_start:qdims_l%j_end,               &
                           qdims_l%k_start:qdims_l%k_end),              &
    cloud_fraction_liquid( qdims_l%i_start:qdims_l%i_end,               & 
                           qdims_l%j_start:qdims_l%j_end,               &
                           qdims_l%k_start:qdims_l%k_end),              &
    cloud_fraction_frozen( qdims_l%i_start:qdims_l%i_end,               & 
                           qdims_l%j_start:qdims_l%j_end,               &
                           qdims_l%k_start:qdims_l%k_end)

  INTEGER ::                                                            &
    ntml (                   pdims%i_start:pdims%i_end,                 & 
                             pdims%j_start:pdims%j_end)

  LOGICAL ::                                                            &
    cumulus (                qdims%i_start:qdims%i_end,                 & 
                             qdims%j_start:qdims%j_end) 
                           ! bl convection flag

! Local variables
  REAL ::                                                               &
   p_layer_centres(          pdims%i_start:pdims%i_end,                 & 
                             pdims%j_start:pdims%j_end,                 &
                                         0:pdims%k_end),                &
              ! pressure at layer centres. Same as p_theta_levels
              ! except bottom level = p_star, and at top = 0.
   p_layer_boundaries(       pdims%i_start:pdims%i_end,                 & 
                             pdims%j_start:pdims%j_end,                 &
                                         0:pdims%k_end),                &
   tl(                       tdims%i_start:tdims%i_end,                 & 
                             tdims%j_start:tdims%j_end,                 &
                                         1:tdims%k_end),                &
   qt(                       qdims%i_start:qdims%i_end,                 & 
                             qdims%j_start:qdims%j_end,                 &
                                         1:qdims%k_end),                &
   qcf_in(                   qdims%i_start:qdims%i_end,                 & 
                             qdims%j_start:qdims%j_end,                 &
                                         1:qdims%k_end),                &
   qcl_out(                  qdims%i_start:qdims%i_end,                 & 
                             qdims%j_start:qdims%j_end,                 &
                                         1:qdims%k_end),                &
   cf_inout(                 qdims%i_start:qdims%i_end,                 & 
                             qdims%j_start:qdims%j_end,                 &
                                         1:qdims%k_end),                &
   cfl_inout(                qdims%i_start:qdims%i_end,                 & 
                             qdims%j_start:qdims%j_end,                 &
                                         1:qdims%k_end),                &
   cff_inout(                qdims%i_start:qdims%i_end,                 & 
                             qdims%j_start:qdims%j_end,                 &
                                         1:qdims%k_end)

  INTEGER ::                                                            &
    large_levels,                                                       &
    levels_per_level

  INTEGER ::                                                            &
    i,j,k,errorstatus     ! loop variables

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('QT_BAL_CLD',zhook_in,zhook_handle)

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      p_layer_centres(i,j,0) = p_star(i,j)
      p_layer_boundaries(i,j,0) = p_star(i,j)
    END DO
  END DO
  DO k = pdims%k_start, pdims%k_end-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
        p_layer_boundaries(i,j,k) = p(i,j,k+1)
      END DO
    END DO
  END DO
  k=pdims%k_end
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
      p_layer_boundaries(i,j,k) = 0.0
    END DO
  END DO

! ----------------------------------------------------------------------
! Section  Convert qT and Tl for input to cloud scheme.
! ----------------------------------------------------------------------

! Create Tl and qT
  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        tl(i,j,k) =                                                     &
               (theta(i,j,k) *exner_theta_levels(i,j,k))                &
                       - (lc * qcl(i,j,k)) / cp
        qt(i,j,k) = q(i,j,k) + qcl(i,j,k)
        qcf_in(i,j,k)=qcf(i,j,k)
        cf_inout(i,j,k)= bulk_cloud_fraction(i,j,k)
        cfl_inout(i,j,k)= cloud_fraction_liquid(i,j,k)
        cff_inout(i,j,k)= cloud_fraction_frozen(i,j,k)
      END DO
    END DO
  END DO

      ! If second cloud ice variable in use then add to qcf
      ! for the cloud scheme call
  IF (l_mcr_qcf2)                                                       &
    qcf_in(:,:,:) = qcf_in(:,:,:) + qcf2(qdims%i_start:qdims%i_end,     &
                                         qdims%j_start:qdims%j_end,     &
                                                     1:qdims%k_end)

! ----------------------------------------------------------------------
! Section BL.4b Call cloud scheme to convert Tl and qT to T, q and qcl
!              calculate bulk_cloud fields from qT and qcf
!               and calculate area_cloud fields.
! ----------------------------------------------------------------------


! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
  levels_per_level = 3
  large_levels = ((qdims%k_end - 2)*levels_per_level) + 2

! DEPENDS ON: ls_arcld
  CALL ls_arcld( p_layer_centres, rhcpt, p_layer_boundaries,            &
                   rhc_row_length, rhc_rows, bl_levels,                 &
                   levels_per_level, large_levels,                      &
                   l_cld_area, l_acf_cusack,l_acf_brooks,               &
                   delta_lambda, delta_phi,                             &
                   fv_cos_theta_latitude,                               &
                   ntml, cumulus, l_mixing_ratio, qcf_in,               &
                   tl, qt, qcl_out,                                     &
                   area_cloud_fraction,  cf_inout,                      &
                   cfl_inout, cff_inout ,                               &
                   errorstatus, mype)

! qt holds q (no halos), tl holds t(no halos),
! qcl_out holds qcl(no halos)
  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        theta(i,j,k) = tl(i,j,k)/exner_theta_levels(i,j,k)
        q(i,j,k) = qt(i,j,k)
        qcl(i,j,k) =qcl_out(i,j,k)
        bulk_cloud_fraction(i,j,k) = cf_inout(i,j,k)
        cloud_fraction_liquid(i,j,k) = cfl_inout(i,j,k)
        cloud_fraction_frozen(i,j,k) = cff_inout(i,j,k)
      END DO
    END DO
  END DO

  IF (lhook) CALL dr_hook('QT_BAL_CLD',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE qt_bal_cld
