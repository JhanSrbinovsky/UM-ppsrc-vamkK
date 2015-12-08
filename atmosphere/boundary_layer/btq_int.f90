! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE BTQ_INT
!  Purpose: To interpolate buoyancy parameters BT and BQ from full
!  levels to half levels
!
! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!   THIS CODE IS WRITTEN TO UMDP 3 PROGRAMMING STANDARDS.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE btq_int (                                                    &
! IN levels
 bl_levels,                                                             &
! IN fields
 z_tq,z_uv,bq,bt,bq_cld,bt_cld,a_qs,a_dqsdt,cf,                         &
! OUT fields
 bqm,btm,bqm_cld,btm_cld,a_qsm,a_dqsdtm,cfm                             &
  )

  USE atm_fields_bounds_mod, ONLY: pdims,tdims
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.
  INTEGER, INTENT(IN) ::                                                &
   bl_levels              ! IN No. of atmospheric levels for which

  REAL, INTENT(IN) ::                                                   &
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels), &
   bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                              ! IN A buoyancy parameter for clear air
                              !    on p,T,q-levels (full levels).
   bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                              ! IN A buoyancy parameter for clear air
                              !    on p,T,q-levels (full levels).
   bq_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                              ! IN A buoyancy parameter for cloudy air
                              !    on p,T,q-levels (full levels).
   bt_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                              ! IN A buoyancy parameter for cloudy air
                              !    on p,T,q-levels (full levels).
   a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                              ! IN Saturated lapse rate factor
                              !    on p,T,q-levels (full levels).
   a_dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           bl_levels),                                                  &
                              ! IN Saturated lapse rate factor
                              !    on p,T,q-levels (full levels).
   cf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                              ! IN Cloud fraction (decimal).

! OUT fields
  REAL, INTENT(OUT) ::                                                  &
   bqm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                              ! OUT A buoyancy parameter for clear air
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   btm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                              ! OUT A buoyancy parameter for clear air
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   bqm_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels),                                                  &
                              ! OUT A buoyancy parameter for cloudy air
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   btm_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels),                                                  &
                              ! OUT A buoyancy parameter for cloudy air
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   a_qsm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                              ! OUT Saturated lapse rate factor
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   a_dqsdtm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
                              ! OUT Saturated lapse rate factor
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   cfm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels)
                              ! OUT Estimate of cloud fraction
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
!-----------------------------------------------------------------------
!    Local and other symbolic constants :-

!  Define local storage.

!  (b) Scalars.

  REAL ::                                                               &
   wk,                                                                  &
               ! Temporary in weighting factor.
   wkm1,                                                                &
               ! Temporary in weighting factor.
   weight1,weight2,weight3

  INTEGER ::                                                            &
   i,j,                                                                 &
               ! Loop counter (horizontal field index).
   k       ! Loop counter (vertical level index).

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('BTQ_INT',zhook_in,zhook_handle)
!-----------------------------------------------------------------------
! 1.  Loop round levels.
!-----------------------------------------------------------------------

  DO k = 2, bl_levels
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end

!-----------------------------------------------------------------------
! 1.1 Calculate buoyancy parameters at half levels,
!     i.e. at level K-1/2, if current level is level K.
!-----------------------------------------------------------------------

        weight1 = 1. / ( z_tq(i,j,k) -                                  &
                         z_tq(i,j,k-1))
        weight2 = z_tq(i,j,k) -                                         &
                  z_uv(i,j,k)
        weight3 = z_uv(i,j,k) -                                         &
                  z_tq(i,j,k-1)
        wkm1 = weight3 * weight1         ! P243.C5 (2nd eqn)
        wk = weight2 * weight1             ! P243.C5 (1st eqn)

        btm(i,j,k-1) = wkm1*bt(i,j,k) + wk*bt(i,j,k-1)
        bqm(i,j,k-1) = wkm1*bq(i,j,k) + wk*bq(i,j,k-1)
        btm_cld(i,j,k-1) = wkm1*bt_cld(i,j,k) + wk*bt_cld(i,j,k-1)
        bqm_cld(i,j,k-1) = wkm1*bq_cld(i,j,k) + wk*bq_cld(i,j,k-1)
        a_qsm(i,j,k-1) = wkm1*a_qs(i,j,k) + wk*a_qs(i,j,k-1)
        a_dqsdtm(i,j,k-1) = wkm1*a_dqsdt(i,j,k) + wk*a_dqsdt(i,j,k-1)
!          CFM(I,j,K-1) = WKM1*CF(I,j,K) + WK*CF(I,j,K-1)
          ! to compensate for lack of resolution of cloud layers
        cfm(i,j,k-1) = MAX( cf(i,j,k), cf(i,j,k-1) )

      END DO !
    END DO !
  END DO ! bl_levels

  IF (lhook) CALL dr_hook('BTQ_INT',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE btq_int
