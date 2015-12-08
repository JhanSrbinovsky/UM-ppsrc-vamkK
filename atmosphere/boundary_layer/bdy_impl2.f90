! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE BDY_IMPL2----------------------------------------------
!!!
!!!  Purpose: Calculate implicit correction to boundary layer fluxes of
!!!           heat, moisture and momentum.
!!!
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version ?, dated ?.
!!!
!!!  System component covered: P24.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Boundary Layer
SUBROUTINE bdy_impl2 (                                                  &

! IN values defining field dimensions and subset to be processed :
 offx ,offy ,row_length,rows,n_rows,                                    &

! IN values defining vertical grid of model atmosphere :
 bl_levels,                                                             &

! IN Substepping information
 l_ftl, l_fqw, l_taux, l_tauy,                                          &

! IN data :
 GAMMA,                                                                 &
 rhokh,rhokm_u,rhokm_v,rdz_charney_grid,rdz_u,rdz_v,                    &

! INOUT data :
 qw,tl,fqw,ftl,taux,tauy,                                               &
 du,dv,ct_ctq,dqw,dtl,cq_cm_u,cq_cm_v,                                  &

! OUT data :
 t_latest,q_latest,rhokh_mix,                                           &

! LOGICAL LTIMER
 ltimer                                                                 &
 )

  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, udims_s, vdims_s, udims_l, vdims_l
  USE atmos_constants_mod, ONLY: cp

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

  INTEGER ::                                                            &
    row_length,                                                         &
                   ! Local number of points on a row
    rows,                                                               &
                   ! Local number of rows in a theta field
    n_rows,                                                             &
                   ! Local number of rows in a v field
    offx ,                                                              &
                   ! Size of small halo in i
    offy       ! Size of small halo in j.

! (b) Defining vertical grid of model atmosphere.

  INTEGER ::                                                            &
   bl_levels                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH

!  In :-

  LOGICAL ::                                                            &
   l_ftl,                                                               &
   l_fqw,                                                               &
   l_taux,                                                              &
   l_tauy


  REAL ::                                                               &
   rhokh(row_length,rows,bl_levels),                                    &
                                   ! IN Exchange coeffs for moisture.
   rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
            bl_levels),                                                 &
                                   ! IN Exchange coefficients for U
   rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
            bl_levels),                                                 &
                                   ! IN Exchange coefficients for V
   rdz_charney_grid(row_length,rows,bl_levels),                         &
                                   ! IN RDZ(,1) is the reciprocal of the
                                   ! height of level 1, i.e. of the
                                   ! middle of layer 1.  For K > 1,
                                   ! RDZ(,K) is the reciprocal
                                   ! of the vertical distance
                                   ! from level K-1 to level K.
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                   ! IN  RDZ (K > 1) on U-grid.
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                   ! IN  RDZ (K > 1) on V-grid.
   GAMMA(bl_levels)         ! IN implicit weights

  LOGICAL :: ltimer               ! Logical switch for TIMER diags

!  In/outs :-

  REAL ::                                                               &
   qw(row_length,rows,bl_levels),                                       &
                                   ! INOUT Total water content, but
                                   !       replaced by specific
                                   !       humidity in LS_CLD.
   tl(row_length,rows,bl_levels),                                       &
                                   ! INOUT Ice/liquid water temperature,
                                   !       but replaced by T in LS_CLD.
   fqw(row_length,rows,bl_levels),                                      &
                                   ! INOUT Moisture flux between layers
                                   !       (kg per square metre per sec)
                                   !       FQW(,1) is total water flux
                                   !       from surface, 'E'.
   ftl(row_length,rows,bl_levels),                                      &
                                   ! INOUT FTL(,K) contains net
                                   !       turbulent sensible heat flux
                                   !       into layer K from below; so
                                   !       FTL(,1) is the surface
                                   !       sensible heat, H. (W/m2)
   taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,            &
         bl_levels),                                                    &
                                   ! INOUT W'ly component of surface
                                   !       wind stress (N/sq m).(On
                                   !       UV-grid with first and last
                                   !       rows undefined or at present,
                                   !       set to  missing data
   tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,            &
         bl_levels),                                                    &
                                   ! INOUT S'ly component of surface
                                   !       wind stress (N/sq m).  On
                                   !       UV-grid; comments as per TAUX
   ct_ctq(row_length,rows,bl_levels),                                   &
                                   ! INOUT Coefficient in T and q
                                   !       tri-diagonal implicit matrix
   cq_cm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                   ! INOUT Coefficient in U tri-diagonal
                                   !       implicit matrix
   cq_cm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                   ! INOUT Coefficient in V tri-diagonal
                                   !       implicit matrix
   dqw(row_length,rows,bl_levels),                                      &
                                   ! INOUT BL increment to q field
   dtl(row_length,rows,bl_levels),                                      &
                                   ! INOUT BL increment to T field
   du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
        bl_levels),                                                     &
                                   ! INOUT BL increment to u wind field
   dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,      &
        bl_levels)
                                   ! INOUT BL increment to v wind field

  REAL ::                                                               &
   t_latest(row_length,rows,bl_levels),                                 &
                                   ! OUT Temperature (K)
   q_latest(row_length,rows,bl_levels),                                 &
                                   ! OUT Specific Humidity (Kg/Kg)
   rhokh_mix(row_length,rows,bl_levels)
                                   ! OUT Exchange coeffs for moisture.


!-----------------------------------------------------------------------

!  Local scalars :-

  INTEGER ::                                                            &
   i,j,                                                                 &
                  ! LOCAL Loop counter (horizontal field index).
   k          ! LOCAL Loop counter (vertical level index).

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('BDY_IMPL2',zhook_in,zhook_handle)

! DEPENDS ON: im_bl_pt2
  CALL im_bl_pt2 (                                                      &
   offx ,offy ,row_length,rows,n_rows,bl_levels,                        &
   l_ftl,l_fqw,l_taux,l_tauy,                                           &
   rhokh(1,1,2),rhokm_u(udims%i_start,udims%j_start,2),                 &
   rhokm_v(vdims%i_start,vdims%j_start,2),                              &
   rdz_charney_grid,rdz_u,rdz_v,GAMMA,                                  &
   ct_ctq,dqw,dtl,cq_cm_u,cq_cm_v,du,dv,                                &
   fqw,ftl,taux,tauy,qw,tl,                                             &
   ltimer                                                               &
  )


!-----------------------------------------------------------------------
!! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
!      Also, IMPL_CAL only updates FTL_TILE(*,1) and FQW_TILE(*,1)
!      over sea points, so copy this to remaining tiles
!-----------------------------------------------------------------------

  DO k = 2, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        ftl(i,j,k) = ftl(i,j,k)*cp
      END DO
    END DO
  END DO

!7.1 Copy T and Q from workspace to INOUT space.

  DO k = 1, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        t_latest(i,j,k)=tl(i,j,k)
        q_latest(i,j,k)=qw(i,j,k)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
! 10 Set RHOKH, the coefficients required for tracer mixing.
!    Required 5B and after due to change in contents of RHOKH in rest
!    of routine.
!-----------------------------------------------------------------------

  DO k = 2, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        rhokh_mix(i,j,k) = rhokh(i,j,k)*                                &
                           rdz_charney_grid(i,j,k)
      END DO
    END DO
  END DO

  IF (lhook) CALL dr_hook('BDY_IMPL2',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bdy_impl2
