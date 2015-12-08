! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------

! Checks on condensate at the beginning of convection if PC2 being used.

SUBROUTINE  conv_pc2_init(rows, row_length, model_levels, wet_levels,      &
                          exner_theta_levels,                              &
                          theta_conv, q_conv, qcl_conv, qcf_conv,          &
                          cf_liquid_conv, cf_frozen_conv, bulk_cf_conv)

USE atm_fields_bounds_mod, ONLY:                                          &
  tdims_s
  
USE atmos_constants_mod, ONLY: cp 
USE water_constants_mod, ONLY: lc, lf 
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Updates convection diagnostics after call to glue each substep of
!   convection.

!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Convection

! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.

! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN)  ::    &
  row_length               & ! Row length
 ,rows                     & ! Number of rows
 ,model_levels             & ! Number of model levels
 ,wet_levels                 ! Number of model wet levels

REAL, INTENT(IN)  ::                                 &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                     tdims_s%j_start:tdims_s%j_end,  & !
                     tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::                           &
  theta_conv(row_length, rows, model_levels)     & ! theta (K)
 ,q_conv(row_length, rows, wet_levels)           & ! water vapour  (kq/kg)
 ,qcl_conv(row_length, rows, wet_levels)         & ! cloud liquid water (kq/kg)
 ,qcf_conv(row_length, rows, wet_levels)         & ! cloud ice water (kq/kg)
 ,cf_liquid_conv(row_length, rows, wet_levels)   & ! liquid cloud fraction
 ,cf_frozen_conv(row_length, rows, wet_levels)   & ! ice cloud fraction
 ,bulk_cf_conv(row_length, rows, wet_levels)       ! total cloud fraction

! ------------------------------------------------------------------------------
! Local declarations:
! ------------------------------------------------------------------------------
INTEGER  ::         &
  i,j,k             &  ! loop counters
 ,n_test               ! counter

INTEGER  ::              &
  n_qcx(1+wet_levels,2)    ! Record  out-of-range input events


! Required by Dr hook 
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('CONV_PC2_INIT',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Prevent negative condensate problems by condensing vapour if needed.
! Input fields are updated without updating increments so change will
! only affect PC2 increments and not the top-level QX_STAR fields.
! ----------------------------------------------------------------------

  n_qcx(1+wet_levels,1) = 0
  n_qcx(1+wet_levels,2) = 0

  DO k=1, wet_levels

    n_test = 0
    DO j=1, rows
      DO i=1, row_length
        IF(qcf_conv(i,j,k)  <   0.0) THEN
!             Freeze some liquid to zero ice content
          qcl_conv(i,j,k) = qcl_conv(i,j,k) + qcf_conv(i,j,k)
          theta_conv(i,j,k) = theta_conv(i,j,k) -                       &
                  ( (qcf_conv(i,j,k) * lf) / (cp * exner_theta_levels(i,j,k)) )
          qcf_conv(i,j,k) = 0.0
          cf_frozen_conv(i,j,k) = 0.0
          bulk_cf_conv(i,j,k) = cf_liquid_conv(i,j,k)
          n_test = n_test + 1
        END IF
      END DO  ! I
    END DO  ! J
    n_qcx(k,1) = n_test
    n_qcx(1+wet_levels,1) = n_qcx(1+wet_levels,1) + n_qcx(k,1)

    n_test = 0
    DO j=1, rows
      DO i=1, row_length
        IF(qcl_conv(i,j,k)  <   0.0) THEN
!             Condense some vapour to zero liquid content
          q_conv(i,j,k) = q_conv(i,j,k) + qcl_conv(i,j,k)
          theta_conv(i,j,k) = theta_conv(i,j,k) -                 &
                   ( (qcl_conv(i,j,k) * lc) / (cp * exner_theta_levels(i,j,k)) )
          qcl_conv(i,j,k) = 0.0
          cf_liquid_conv(i,j,k) = 0.0
          bulk_cf_conv(i,j,k) = cf_frozen_conv(i,j,k)
          n_test = n_test + 1
        END IF
      END DO  ! I
    END DO  ! J
    n_qcx(k,2) = n_test
    n_qcx(1+wet_levels,2) = n_qcx(1+wet_levels,2) + n_qcx(k,2)


!     Might also be necessary to place a limit on supersaturation.
!     Convection copes but cloud scheme response is less predictable.
!     Would need eg. LS_CLD_C to reduce supersaturation consistently.

    DO j=1, rows
      DO i=1, row_length

!           Ensure that input values of cloud fields lie within the
!           bounds of physical possibility (should do nothing).

        bulk_cf_conv(i,j,k)   = MAX(0., ( MIN(1., bulk_cf_conv(i,j,k)) ))
        cf_liquid_conv(i,j,k) = MAX(0., ( MIN(1., cf_liquid_conv(i,j,k)) ))
        cf_frozen_conv(i,j,k) = MAX(0., ( MIN(1., cf_frozen_conv(i,j,k)) ))
      END DO  ! I
    END DO  ! J

  END DO  ! K

  IF(n_qcx(1+wet_levels,1)  >   0)                                   &
    WRITE(6,*) 'Qcf < 0 fixed by PC2 ',n_qcx(1+wet_levels,1)
  IF(n_qcx(1+wet_levels,2)  >   0)                                   &
    WRITE(6,*) 'Qcl < 0 fixed by PC2 ',n_qcx(1+wet_levels,2)

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('CONV_PC2_INIT',zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
END SUBROUTINE
