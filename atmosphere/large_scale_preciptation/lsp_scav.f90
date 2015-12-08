! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE LSP_SCAV-----------------------------------------------
!    Purpose: Scavenge aerosol by large scale precipitation.
!
!    Programming standard: Unified Model Documentation Paper No 3,
!                          Version 7, dated 11/3/93.
!
!    Logical component covered: Part of P26.
!
!    Documentation: Unified Model Documentation Paper No 26.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation
MODULE lsp_scav_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_scav( points, rain, snow, droplet_flux, aerosol )

  ! General modules
  USE timestep_mod, ONLY: timestep

  ! Dr Hook Modules
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE
  INTEGER ::                                                            &
                      ! Input integer scalar :-
   points         ! IN Number of points to be processed.

  REAL ::                                                               &
                      ! Input real arrays :-
   rain(points),                                                        &
                      ! IN Rate of rainfall in this layer from
!                     !       above
!*                    !       (kg per sq m per s).
   snow(points),                                                        &
                      ! IN Rate of snowfall in this layer from
!                     !       above
!*                    !       (kg per sq m per s).
   droplet_flux(points)
                      ! In Rate of droplet settling in this layer
                      !       from above
                      !       (kg per sq m per s).
  REAL ::                                                               &
                      ! Updated real arrays :-
   aerosol(points) ! INOUT Aerosol mixing ratio


!-----------------------------------------------------------------------
!  Define local scalars.
!-----------------------------------------------------------------------
!  (a) Reals effectively expanded to workspace by the Cray (using
!      vector registers).
  REAL ::                                                               &
                      ! Real workspace.
   krain,ksnow
  PARAMETER(krain=1.0e-4,ksnow=1.0e-4)
  REAL ::                                                               &
                      ! Real workspace.
   rrain,rsnow
!  (b) Others.
  INTEGER :: i       ! Loop counter (horizontal field index).

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


! Overall rate = KRAIN*(R) where R is in mm/hr=kg/m2/s*3600.0
  IF (lhook) CALL dr_hook('LSP_SCAV',zhook_in,zhook_handle)
  rrain=krain*timestep*3600.0
  rsnow=ksnow*timestep*3600.0
  DO i = 1, points
    aerosol(i)=aerosol(i)/                                              &
          (1.0+rrain*rain(i)+rrain*droplet_flux(i)+rsnow*snow(i))
  END DO
  IF (lhook) CALL dr_hook('LSP_SCAV',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_scav
END MODULE lsp_scav_mod
