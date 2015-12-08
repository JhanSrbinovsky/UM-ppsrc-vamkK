! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE KMKH---------------------------------------------------
!
!  Purpose: To set the turbulent mixing coefficients KM and KH
!           (Note: should be used after any vertical interpolation
!                  but before any horizontal interpolation.)
!
!  Programming standard:
!
!  Documentation: UMDP No.24
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE kmkh (                                                       &
! IN data
 bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                                &
 ntml,cumulus,ntdsc,dsc,sml_disc_inv,dsc_disc_inv,                      &
 weight_1dbl,weight_1dbl_rho,                                           &
! INOUT data
 rhokm,rhokh,rhokmz,rhokhz,rhokm_top,rhokh_top                          &
 )

  USE bl_option_mod, ONLY:                                              &
      Keep_Ri_FA, off, on, Kprof_cu, except_disc_inv
  USE bl_diags_mod, ONLY: strnewbldiag
  USE cv_run_mod, ONLY:                                                 &
      l_param_conv
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: pdims,pdims_s,tdims

  IMPLICIT NONE

! IN arguments
  INTEGER, INTENT(IN) ::                                                &
   bl_levels

  LOGICAL, INTENT(IN) ::                                                &
   cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                ! IN flag for Cu in the bl
   dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! IN Flag set if decoupled stratocumulus layer found

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER,INTENT(IN) ::                                                 &
   nSCMDpkgs              ! No of SCM diagnostics packages

  LOGICAL,INTENT(IN) ::                                                 &
   L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages

! Declaration of new BL diagnostics.
  TYPE (strnewbldiag) :: BL_diag

  INTEGER, INTENT(IN) ::                                                &
   ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                ! IN Number of model levels in the
!                                   turbulently mixed layer.
   ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                ! IN Top level for turb mixing in
!                                   cloud layer
   sml_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                ! IN Flags for whether discontinuous
   dsc_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! IN inversions are diagnosed

  REAL, INTENT(IN) ::                                                   &
   weight_1dbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               bl_levels),                                              &
   weight_1dbl_rho(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, &
               bl_levels)
                              ! IN Weighting applied to 1D BL scheme, 
                              !    to blend with Smagorinsky scheme,
                              !    on theta and rho levels
! INOUT arguments
  REAL, INTENT(INOUT) ::                                                &
   rhokmz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          2:bl_levels),                                                 &
                              ! INOUT Non-local turbulent mixing
!                                  coefficient for momentum.
   rhokhz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:bl_levels),                                                 &
                              ! INOUT Non-local turbulent mixing
!                                  coefficient for heat and moisture
   rhokm_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             2:bl_levels),                                              &
                              ! INOUT Non-local top-down turbulent
                              !    mixing coefficient for momentum.
   rhokh_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             2:bl_levels)
                              ! INOUT Non-local top-down turbulent
                              !    mixing coefficient for heat
                              !    and moisture.
  REAL,INTENT(INOUT) ::                                                 &
   rhokm(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,   &
         bl_levels),                                                    &
                              ! INOUT Layer K-1 - to - layer K
!                                     turbulent mixing coefficient
!                                     for momentum.
   rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                              ! INOUT Layer K-1 - to - layer K
!                                     turbulent mixing coefficient
!                                     for heat and moisture.
!----------------------------------------------------------------------
!  Define local storage.

  INTEGER ::                                                            &
   i,j,                                                                 &
                       ! Loop counter (horizontal field index).
   k             ! Loop counter (vertical level index).

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('KMKH',zhook_in,zhook_handle)


  IF (Keep_Ri_FA == on) THEN
!-----------------------------------------------------------------------
! Set local K's to zero at the LCL in cumulus and at the
! top of a turbulent layer with a well-defined inversion
!-----------------------------------------------------------------------
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF ( ( cumulus(i,j) .OR. sml_disc_inv(i,j) == 1) .AND.        &
               (k == ntml(i,j)+1 .OR. k == ntml(i,j)+2) ) THEN
            rhokh(i,j,k) = 0.0
            rhokm(i,j,k) = 0.0
          END IF

          IF ( dsc_disc_inv(i,j)  ==  1 .AND.                           &
               (k == ntdsc(i,j)+1 .OR. k == ntdsc(i,j)+2) ) THEN
            rhokh(i,j,k) = 0.0
            rhokm(i,j,k) = 0.0
          END IF

        END DO ! P_POINTS,i
      END DO ! P_POINTS,j
    END DO ! BL_LEVELS

  ELSE IF (Keep_Ri_FA == except_disc_inv) THEN
!-----------------------------------------------------------------------
        ! Reduce local K's only at the top of a turbulent
        ! layer with a well-defined inversion
!-----------------------------------------------------------------------
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF ( sml_disc_inv(i,j)  ==  1 .AND.                           &
               (k == ntml(i,j)+1 .OR. k == ntml(i,j)+2) ) THEN
            rhokh(i,j,k) = (1.0-weight_1dbl(i,j,k))*rhokh(i,j,k) 
            rhokm(i,j,k) = (1.0-weight_1dbl(i,j,k))*rhokm(i,j,k) 
          END IF

          IF ( dsc_disc_inv(i,j)  ==  1 .AND.                           &
               (k == ntdsc(i,j)+1 .OR. k == ntdsc(i,j)+2) ) THEN
            rhokh(i,j,k) = (1.0-weight_1dbl(i,j,k))*rhokh(i,j,k) 
            rhokm(i,j,k) = (1.0-weight_1dbl(i,j,k))*rhokm(i,j,k) 
          END IF

        END DO ! P_POINTS,i
      END DO ! P_POINTS,j
    END DO ! BL_LEVELS

  ELSE
!-----------------------------------------------------------------------
! Set local K's to zero from the LCL in cumulus and from the
! top of a turbulent layer with a well-defined inversion
!-----------------------------------------------------------------------
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF(cumulus(i,j) .AND. ( (l_param_conv .AND. k >  ntml(i,j))       &
             .OR. (.NOT. l_param_conv .AND. k >= ntml(i,j)) )) THEN
            rhokh(i,j,k)=0.0
            rhokm(i,j,k)=0.0
          END IF

          IF ( dsc_disc_inv(i,j)  ==  1 .AND. k  >   ntdsc(i,j) ) THEN
            rhokh(i,j,k) = 0.0
            rhokm(i,j,k) = 0.0
          END IF

          IF ( sml_disc_inv(i,j)  ==  1 .AND. k  >   ntml(i,j) ) THEN
              !   This also means no local mixing within any DSC layer
            rhokh(i,j,k) = 0.0
            rhokm(i,j,k) = 0.0
          END IF

        END DO ! P_POINTS,i
      END DO ! P_POINTS,j
    END DO ! BL_LEVELS

  END IF ! test on Keep_Ri_FA

  IF (Kprof_cu == off) THEN
!-----------------------------------------------------------------------
! Set non-local K's to zero at the LCL in cumulus layers,
! including level NTML if not l_param_conv convection scheme
!-----------------------------------------------------------------------
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        IF(cumulus(i,j) .AND. ( (l_param_conv .AND. k == ntml(i,j)+1)       &
             .OR. (.NOT. l_param_conv .AND.                                 &
                         k >= ntml(i,j).AND.k <  ntml(i,j)+2) )) THEN
          rhokhz(i,j,k)=0.0
          rhokmz(i,j,k)=0.0
          rhokh_top(i,j,k)=0.0
          rhokm_top(i,j,k)=0.0
        END IF
      END DO ! P_POINTS,i
    END DO ! P_POINTS,j
  END DO ! BL_LEVELS

  END IF  ! test on Kprof_cu
!-----------------------------------------------------------------------
! Save diffusion coefficients for diagnostics
!-----------------------------------------------------------------------
  IF (BL_diag%L_rhoKmloc) THEN
    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%rhoKmloc(i,j,k)=rhokm(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%L_rhoKhloc) THEN
    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%rhoKhloc(i,j,k)=rhokh(i,j,k)
        END DO
      END DO
    END DO
  END IF      

  IF (BL_diag%L_rhoKmsurf) THEN
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%rhoKmsurf(i,j,k)=weight_1dbl(i,j,k)*rhokmz(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%L_rhoKhsurf) THEN
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%rhoKhsurf(i,j,k)=weight_1dbl_rho(i,j,k)*rhokhz(i,j,k)
        END DO
      END DO
    END DO
  END IF      

  IF (BL_diag%L_rhoKmSc) THEN
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%rhoKmSc(i,j,k)=weight_1dbl(i,j,k)*rhokm_top(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%L_rhoKhSc) THEN
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%rhoKhSc(i,j,k)=weight_1dbl_rho(i,j,k)*rhokh_top(i,j,k)
        END DO
      END DO
    END DO
  END IF
!-----------------------------------------------------------------------
! Set KM and KH to be the maximum of the local and non-local
! values andstore RHO_KM on P-grid for output.
!-----------------------------------------------------------------------
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        rhokh(i,j,k) = MAX( rhokh(i,j,k) ,                              &
               weight_1dbl_rho(i,j,k)*(rhokhz(i,j,k)+rhokh_top(i,j,k)) )
        rhokm(i,j,k) = MAX( rhokm(i,j,k) ,                              &
               weight_1dbl(i,j,k)*(rhokmz(i,j,k)+rhokm_top(i,j,k)) )

      END DO ! P_POINTS,i
    END DO ! P_POINTS,j
  END DO ! BL_LEVELS


  IF (lhook) CALL dr_hook('KMKH',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE kmkh
