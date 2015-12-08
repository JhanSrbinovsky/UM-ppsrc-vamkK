! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE BDY_IMPL1 ----------------------------------------------
!!!
!!!  Purpose: Calculate downward sweep of matrix for increments to
!!!           T and Q in the boundary layer, by calling subroutine
!!!           IM_BL_PT1.
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------
!!  Arguments :-
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
SUBROUTINE bdy_impl1 (                                                  &
 halo_i,halo_j,offx ,offy ,row_length,rows,n_rows,bl_levels,            &
 global_row_length,n_proc, n_procy, proc_row_group,at_extremity,        &
 q,qcl,qcf,q_latest,qcl_latest,qcf_latest,                              &
 t,t_latest,                                                            &
 dtrdz_charney_grid,dtrdz_u,dtrdz_v,                                    &
 rhokh,rhokm_u,rhokm_v,                                                 &
 rdz_charney_grid,rdz_u,rdz_v,GAMMA,                                    &
 du_nt,dv_nt,                                                           &
 fqw,ftl,tau_x,tau_y,                                                   &
 qw,tl,                                                                 &
 ct_ctq,dqw,dtl,cq_cm_u,cq_cm_v,du,dv,                                  &
 ltimer                                                                 &
)

  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, udims_s, vdims_s, udims_l, vdims_l
  USE atmos_constants_mod, ONLY: cp
  USE water_constants_mod, ONLY: lc, lf

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

  LOGICAL ::                                                            &
   ltimer

  INTEGER ::                                                            &
    row_length,                                                         &
                                 ! IN Number of points on a row
    rows,                                                               &
                                 ! IN Number of rows in a theta field
    n_rows,                                                             &
                                 ! IN Number of rows in a v field
    halo_i,                                                             &
                                 ! IN Size of halo in i direction.
    halo_j,                                                             &
                                 ! IN Size of halo in j direction.
    offx ,                                                              &
                                 ! IN Size of small halo in i
    offy ,                                                              &
                                 ! IN Size of small halo in j.
    n_proc,                                                             &
                                ! Total number of processors
    n_procy,                                                            &
                                ! Number of procs N/S
    global_row_length,                                                  &
                                ! number of points per global row
    proc_row_group,                                                     &
                                ! Group id for processors on the same row
                                 
    bl_levels                 ! IN No. of atmospheric levels for
                                 !    which boundary layer fluxes are
                                 !    calculated.

  LOGICAL ::                                                            &
    at_extremity(4)  ! Indicates if this processor is at north,
                         !  south,east or west of the processor grid 
 

  REAL ::                                                               &
   q(row_length,rows,bl_levels),                                        &
                                   ! IN specific humidity
   qcl(row_length,rows,bl_levels),                                      &
                                   ! IN Cloud liquid water
   qcf(row_length,rows,bl_levels),                                      &
                                   ! IN Cloud ice (kg per kg air)
   q_latest(row_length,rows,bl_levels),                                 &
                                   ! IN specific humidity
   qcl_latest(row_length,rows,bl_levels),                               &
                                   ! IN Cloud liquid water
   qcf_latest(row_length,rows,bl_levels),                               &
                                   ! IN Cloud ice (kg per kg air)
   t(row_length,rows,bl_levels),                                        &
                                   ! IN temperature
                                   !    Latest estimates to time
                                   !    level n+1 values
   t_latest(row_length,rows,bl_levels),                                 &
                                   ! IN temperature
   dtrdz_charney_grid(row_length,rows,bl_levels),                       &
                                   ! IN dz for bottom BL_LEVELS
   dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                   ! IN -g.dt/dp for model wind layers
   dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                   ! IN -g.dt/dp for model wind layers
   rhokh(row_length,rows,bl_levels),                                    &
                                   ! IN Exchange coeff for FTL above
                                   !    surface.
   rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
            bl_levels),                                                 &
                                   ! IN Exchange coefficients for
                                   !    momentum, on U-grid with
                                   !    first and last rows ignored.
                                   !    for K>=2 (from KMKH).
   rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
            bl_levels),                                                 &
                                   ! IN Exchange coefficients for
                                   !    momentum, on V-grid with
                                   !    first and last rows ignored.
                                   !    for K>=2 (from KMKH).
   rdz_charney_grid(row_length,rows,bl_levels),                         &
                                   ! IN RDZ(,1) is the reciprocal of the
                                   ! height of level 1, i.e. of the
                                   ! middle of layer 1.  For K > 1,
                                   ! RDZ(,K) is the reciprocal
                                   ! of the vertical distance
                                   ! from level K-1 to level K.
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                   ! IN Reciprocal of the vertical
                                   !    distance from level K-1 to
                                   !    level K. (K > 1) on wind levels
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                   ! IN Reciprocal of the vertical
                                   !    distance from level K-1 to
                                   !    level K. (K > 1) on wind levels
   GAMMA(bl_levels),                                                    &
                                   ! IN Implicit weighting.
   du_nt(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,   &
          bl_levels),                                                   &
                                   ! IN u non-turbulent increments.
   dv_nt(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,   &
          bl_levels),                                                   &
                                   ! IN v non-turbulent increments.
   fqw(row_length,rows,bl_levels),                                      &
                                   ! IN Flux of QW (ie., for surface,
                                   !    total evaporation). Kg/sq m/s
   ftl(row_length,rows,bl_levels),                                      &
                                   ! IN Flux of TL (ie., for surface,
                                   !    H/Cp where H is sensible heat
                                   !    in W per sq m).
   tau_x(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
         bl_levels),                                                    &
                                   ! IN x-component of turbulent
                                   !    stress at levels k-1/2;
                                   !    eg. TAUX(,1) is surface stress.
                                   !    U-grid, 1st and last rows set
                                   !    to "missing data". (N/sq m)
                                   !    IN as "explicit" fluxes from
                                   !    ex_flux_uv, OUT as "implicit
   tau_y(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
         bl_levels)
                                   ! IN y-component of turbulent
                                   !    stress at levels k-1/2;
                                   !    eg. TAUX(,1) is surface stress.
                                   !    V-grid, 1st and last rows set
                                   !    to "missing data". (N/sq m)
                                   !    IN as "explicit" fluxes from
                                   !    ex_flux_uv, OUT as "implicit


  REAL ::                                                               &
   qw(row_length, rows, bl_levels),                                     &
                                   ! OUT total water
   tl(row_length, rows, bl_levels),                                     &
                                   ! OUT liquid water temperature
   ct_ctq(row_length,rows,bl_levels),                                   &
                                   ! OUT Coefficient in T and q
                                   !     tri-diagonal implicit matrix
   cq_cm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                   ! OUT Coefficient in U and V
                                   !     tri-diagonal implicit matrix
   cq_cm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                   ! OUT Coefficient in U and V
                                   !     tri-diagonal implicit matrix
   dqw(row_length,rows,bl_levels),                                      &
                                   ! OUT BL increment to q field
   dtl(row_length,rows,bl_levels),                                      &
                                   ! OUT BL increment to T field
   du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
        bl_levels),                                                     &
                                   ! INOUT BL increment to u wind field
   dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,      &
        bl_levels) 
                                   ! INOUT BL increment to v wind field

! Workspace
  REAL ::                                                               &
   dqw_nt(row_length,rows,bl_levels),                                   &
   dtl_nt(row_length,rows,bl_levels)


! Derived local parameters.

  REAL :: lcrcp,ls,lsrcp

  PARAMETER (                                                           &
   lcrcp=lc/cp,                                                         &
                             ! Evaporation-to-dT conversion factor.
   ls=lf+lc,                                                            &
                             ! Latent heat of sublimation.
   lsrcp=ls/cp                                                          &
                             ! Sublimation-to-dT conversion factor.
    )


!  Local scalars :-

  INTEGER ::                                                            &
   i,j,                                                                 &
                ! Loop counter (horizontal field index).
   k        ! Loop counter (vertical index).

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle




  IF (lhook) CALL dr_hook('BDY_IMPL1',zhook_in,zhook_handle)

  DO k = 1, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        qw(i,j,k) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
        tl(i,j,k) = t(i,j,k) - lcrcp*qcl(i,j,k) - lsrcp*qcf(i,j,k)
        dqw_nt(i,j,k) = q_latest(i,j,k) + qcl_latest(i,j,k)             &
                        + qcf_latest(i,j,k) - qw(i,j,k)
        dtl_nt(i,j,k) = t_latest(i,j,k)                                 &
                        - lcrcp * qcl_latest(i,j,k)                     &
                        - lsrcp * qcf_latest(i,j,k)                     &
                        - tl(i,j,k)
      END DO
    END DO
  END DO

! DEPENDS ON: im_bl_pt1
  CALL im_bl_pt1 (                                                      &
   offx ,offy ,row_length,rows,n_rows,bl_levels,                        &
   global_row_length,n_proc, n_procy, proc_row_group,at_extremity,      &
   halo_i,halo_j,                                                       &
   dtrdz_charney_grid,dtrdz_u,dtrdz_v,                                  &
   rhokh(1,1,2),rhokm_u(udims%i_start,udims%j_start,2),                 &
   rhokm_v(vdims%i_start,vdims%j_start,2),                              &
   rdz_charney_grid,rdz_u,rdz_v,GAMMA,                                  &
   dqw_nt,dtl_nt,du_nt,dv_nt,                                           &
   fqw,ftl,tau_x,tau_y,                                                 &
   ct_ctq,dqw,dtl,cq_cm_u,cq_cm_v,du,dv,                                &
   ltimer                                                               &
  )

  IF (lhook) CALL dr_hook('BDY_IMPL1',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bdy_impl1
