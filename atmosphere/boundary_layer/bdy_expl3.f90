! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE BDY_EXPL3----------------------------------------------
!
!  Purpose: Calculate explicit fluxes taux and tauy on their
!           native (u or v) grids
!
! Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 24.
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE bdy_expl3 (                                                  &
! IN grid related variables 
 bl_levels,                                                             &
! IN SCM diags
 nSCMDpkgs,L_SCMDiags,                                                  &
! IN variables used in flux calculations
 u, v, u_0, v_0, rhokm_u_land, rhokm_v_land, flandfac_u, flandfac_v,    &
 rhokm_u_ssi, rhokm_v_ssi, fseafac_u, fseafac_v, flandg_u, flandg_v,    &
 zh, rdz_u, rdz_v, rhokm_u, rhokm_v, taux_fd_u, tauy_fd_v,              &
 rhogamu_u, rhogamv_v, f_ngstress_u, f_ngstress_v,                      &
! OUT explicit momentum fluxes
 taux_land, tauy_land, taux_ssi, tauy_ssi, taux, tauy                   &
 )

  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE atm_fields_bounds_mod, ONLY:   udims, vdims, udims_s, vdims_s,    &
                                     pdims
  USE switches, ONLY: l_ctile
  USE bl_option_mod, ONLY: buddy_sea, on, ng_stress, formdrag,          &
       explicit_stress
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

!  Inputs :-

  INTEGER, INTENT(IN) ::                                                &
   bl_levels
                   ! IN Max. no. of "boundary" levels

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN) ::                                                &
   nSCMDpkgs             ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
   L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

  REAL, INTENT(IN) ::                                                   &
   u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,       &
     bl_levels),                                                        &
   v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,       &
     bl_levels),                                                        &
                   ! horizontal winds
   u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end),            &
   v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),            &
                   ! surface currents
   rhokm_u_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end),   &
   rhokm_u_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end),    &
   rhokm_v_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),   &
   rhokm_v_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),    &
   rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
           bl_levels),                                                  &
   rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
           bl_levels),                                                  &
                   ! rho * Km terms on u/v grids
   taux_fd_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,       &
         bl_levels),                                                    &
   tauy_fd_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,       &
         bl_levels),                                                    &
                   ! explicit form drag
   flandfac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),     &
   flandfac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),     &
   fseafac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
   fseafac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
                   ! land and sea scaling factors for coastal tiling
   flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
   flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),       &
                   ! land fraction on u/v grids
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
         2:bl_levels),                                                  &
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
         2:bl_levels),                                                  &
                   ! 1 / distance between levels
! Counter gradient stress terms for 1A Bl scheme
   rhogamu_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,       &
                 2:bl_levels),                                          &
   rhogamv_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,       &
                 2:bl_levels),                                          &
! Counter gradient stress terms for other BL schemes
   f_ngstress_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,    &
                 2:bl_levels),                                          &
   f_ngstress_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,    &
                 2:bl_levels),                                          &
   zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! boundary layer depth

! Outputs :-
  REAL, INTENT(OUT) ::                                                  &
   taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
                   ! Taux over land part of grid box.
   taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                   ! Taux over sea part of grid box.
   tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
                   ! Tauy over land part of grid box.
   tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),       &
                   ! Tauy over sea part of grid box.
   taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,            &
        bl_levels),                                                     &
   tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,            &
        bl_levels)
                   ! explicit momentum fluxes
!-----------------------------------------------------------------------

! local variables :-







  REAL ::                                                               &
   tau_grad_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,      &
              bl_levels),                                               &
   tau_non_grad_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,  &
              bl_levels),                                               &
   tau_grad_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,      &
              bl_levels),                                               &
   tau_non_grad_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,  &
              bl_levels)
  REAL, DIMENSION (:,:), ALLOCATABLE ::                                 &
   zh_uv            ! zh on u or v points for ex_flux_uv call
  Integer ::                                                            &
   i, j

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
  IF (lhook) CALL dr_hook('BDY_EXPL3',zhook_in,zhook_handle)
!-----------------------------------------------------------------------

! surface fluxes
  IF (l_ctile .AND. buddy_sea == on) THEN

    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        taux_land(i,j) = rhokm_u_land(i,j)* u(i,j,1) * flandfac_u(i,j)
        taux_ssi(i,j)  = rhokm_u_ssi(i,j) * ( u(i,j,1) - u_0(i,j) )     &
                                          * fseafac_u(i,j)
        taux(i,j,1) = flandg_u(i,j)*taux_land(i,j)                      &
                      + (1.-flandg_u(i,j))*taux_ssi(i,j)
      END DO
    END DO

    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        tauy_land(i,j) = rhokm_v_land(i,j)* v(i,j,1) * flandfac_v(i,j)
        tauy_ssi(i,j)  = rhokm_v_ssi(i,j) * ( v(i,j,1) - v_0(i,j) )     &
                                          * fseafac_v(i,j)
        tauy(i,j,1) = flandg_v(i,j)*tauy_land(i,j)                      &
                      + (1.-flandg_v(i,j))*tauy_ssi(i,j)
      END DO
    END DO

  ELSE   ! Standard code

    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        taux_land(i,j) = rhokm_u_land(i,j) * u(i,j,1)
        taux_ssi(i,j) = rhokm_u_ssi(i,j) * ( u(i,j,1) - u_0(i,j) )

        taux(i,j,1) = flandg_u(i,j)*taux_land(i,j)                      &
                      + (1.-flandg_u(i,j))*taux_ssi(i,j)
      END DO
    END DO

    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        tauy_land(i,j) = rhokm_v_land(i,j) * v(i,j,1)
        tauy_ssi(i,j) = rhokm_v_ssi(i,j) * ( v(i,j,1) - v_0(i,j) )

        tauy(i,j,1) = flandg_v(i,j)*tauy_land(i,j)                      &
                      + (1.-flandg_v(i,j))*tauy_ssi(i,j)
      END DO
    END DO

  END IF

! above surface fluxes
!-----------------------------------------------------------------------
! 5.6 Calculation of explicit fluxes of U and V, other BL schemes method
!-----------------------------------------------------------------------

! Copy ZH into an array defined on u-points
! In principle, ZH should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
  ALLOCATE (zh_uv(udims%i_start:udims%i_end,udims%j_start:udims%j_end))
IF (l_vatpoles) THEN
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
       zh_uv(i,j) = zh(i+1,j)
    END DO
  END DO
ELSE
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
       zh_uv(i,j) = zh(i,j)
    END DO
  END DO
END IF ! vatpoles

! DEPENDS ON: ex_flux_uv
  CALL ex_flux_uv (                                                     &
                        ! For U
    udims%i_start, udims%i_end, udims%j_start, udims%j_end,             &
    udims_s%halo_i , udims_s%halo_j, bl_levels,                         &
    u,zh_uv,rdz_u,rhokm_u,f_ngstress_u,taux_fd_u,taux,                  &
    tau_grad_u,tau_non_grad_u                                           &
    )
  DEALLOCATE (zh_uv)

! Copy ZH into an array defined on v-points
! In principle, ZH should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
  ALLOCATE (zh_uv(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end))
IF (l_vatpoles) THEN
  DO i = vdims%i_start, vdims%i_end
    DO j = vdims%j_start, vdims%j_end-1
       zh_uv(i,j) = zh(i,j+1)
    END DO
    zh_uv(i,vdims%j_end) = zh(i,pdims%j_end)
  END DO
ELSE
  DO i = vdims%i_start, vdims%i_end
    DO j = vdims%j_start, vdims%j_end
       zh_uv(i,j) = zh(i,j)
    END DO
  END DO
END IF ! vatpoles

! DEPENDS ON: ex_flux_uv
  CALL ex_flux_uv (                                                     &
                        ! For V
    vdims%i_start, vdims%i_end, vdims%j_start, vdims%j_end,             &
    vdims_s%halo_i , vdims_s%halo_j, bl_levels,                         &
    v,zh_uv,rdz_v,rhokm_v,f_ngstress_v,tauy_fd_v,tauy,                  &
    tau_grad_v,tau_non_grad_v                                           &
    )
  DEALLOCATE (zh_uv)


! Add additional orographic stress to surface stress over land

  IF(formdrag  ==  explicit_stress) THEN

    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        IF(flandg_u(i,j) >  0.0)THEN
          taux_land(i,j) = taux_land(i,j) + taux_fd_u(i,j,1)
        END IF
        IF(flandg_u(i,j) <  1.0)THEN
          taux_ssi(i,j) = taux_ssi(i,j) + taux_fd_u(i,j,1)
        END IF
      END DO
    END DO

    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        IF(flandg_v(i,j) >  0.0)THEN
          tauy_land(i,j) = tauy_land(i,j) + tauy_fd_v(i,j,1)
        END IF
        IF(flandg_v(i,j) <  1.0)THEN
          tauy_ssi(i,j) = tauy_ssi(i,j) + tauy_fd_v(i,j,1)
        END IF
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
  IF (lhook) CALL dr_hook('BDY_EXPL3',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bdy_expl3
