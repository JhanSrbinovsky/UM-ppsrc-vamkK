! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_helmholtz_mod

REAL, SAVE, ALLOCATABLE ::  ec_vol(:,:,:)
REAL, SAVE, ALLOCATABLE ::  ec_area(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_u(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_v(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_w(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_theta(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_etadot(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_vol(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_rhox(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_rhoy(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_rhoz(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_p(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_b(:,:,:)
REAL, SAVE              ::  hm_pp

CONTAINS


SUBROUTINE eg_init_helmholtz(row_length,rows,n_rows,model_levels,offx,offy)
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod, ONLY: ereport

IMPLICIT NONE
!
! Description: allocates constants required for the
!              solution to the  Helmholtz problem.
!  
!
! Method: Appendix F, ENDGame formulation version 1.01
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


INTEGER row_length,rows,n_rows,model_levels,offx,offy

INTEGER alloc_stat

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('eg_init_helmholtz',zhook_in,zhook_handle)


ALLOCATE ( ec_vol(row_length,rows,model_levels), STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  ec_area(row_length,rows,model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_u(-offx:row_length-1+offx,1-offy:rows+offy,            &
                  model_levels), STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_v(1-offx:row_length+offx,-offy:n_rows-1+offy,          &
                 model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_w(1-offx:row_length+offx,1-offy:rows+offy,             &
                 0:model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_b(1-offx:row_length+offx,1-offy:rows+offy,             &
                 0:model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")


ALLOCATE (  hm_theta(1-offx:row_length+offx,                          &
                1-offy:rows+offy,0:model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")      

ALLOCATE (  hm_etadot(1-offx:row_length+offx,1-offy:rows+offy,        &
            0:model_levels), STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")                          

ALLOCATE (  hm_vol(1-offx:row_length+offx,1-offy:rows+offy,           &
                       model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_rhox(-offx:row_length-1+offx,1-offy:rows+offy,         &
                         model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_rhoy(1-offx:row_length+offx,-offy:n_rows-1+offy,       &
          model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")                    

ALLOCATE (  hm_rhoz(1-offx:row_length+offx, 1-offy:rows+offy,         &
          0:model_levels),          STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_p(1-offx:row_length+offx,                              &
                1-offy:rows+offy,0:model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

IF (lhook) CALL dr_hook('eg_init_helmholtz',zhook_out,zhook_handle)

RETURN

END SUBROUTINE eg_init_helmholtz

SUBROUTINE eg_destroy_helmholtz()
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
!
! Description: frees constants required for the
!              solution to the  Helmholtz problem.
!  
!
! Method: Appendix F, ENDGame formulation version 1.01
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('eg_destroy_helmholtz',zhook_in,zhook_handle)

DEALLOCATE ( hm_p )
DEALLOCATE ( hm_rhoz )
DEALLOCATE ( hm_rhoy )
DEALLOCATE ( hm_rhox )
DEALLOCATE ( hm_vol )
DEALLOCATE ( hm_etadot )
DEALLOCATE ( hm_theta )
DEALLOCATE ( hm_w )
DEALLOCATE ( hm_v )
DEALLOCATE ( hm_u )
DEALLOCATE ( ec_area )
DEALLOCATE ( ec_vol )

IF (lhook) CALL dr_hook('eg_destroy_helmholtz',zhook_out,zhook_handle)

RETURN

END SUBROUTINE eg_destroy_helmholtz




SUBROUTINE eg_helmholtz_update_integrity()

USE integrity_mod

IMPLICIT NONE

IF (integrity_test)                                                  &
  CALL update_hash_m(hm_p,            SIZE(hm_p),            'hm_p_',&
                     hm_rhoz,         SIZE(hm_rhoz),         'hmroz',&
                     hm_rhoy,         SIZE(hm_rhoy),         'hmroy',&
                     hm_rhox,         SIZE(hm_rhox),         'hmrox',&
                     hm_vol,          SIZE(hm_vol),          'hmvol',&
                     hm_etadot,       SIZE(hm_etadot),       'hm_et',&
                     hm_theta,        SIZE(hm_theta),        'hm_th',&
                     hm_w,            SIZE(hm_w),            'hm_w_',&
                     hm_v,            SIZE(hm_v),            'hm_v_',&
                     hm_u,            SIZE(hm_u),            'h_mu_',&
                     ec_area,         SIZE(ec_area),         'ecare',&
                     ec_vol,          SIZE(ec_vol),          'ecvol')

END SUBROUTINE eg_helmholtz_update_integrity




SUBROUTINE eg_helmholtz_check_integrity()

USE integrity_mod

IMPLICIT NONE

IF (integrity_test)                                                  &
  CALL  check_hash_m(hm_p,            SIZE(hm_p),            'hm_p_',&
                     hm_rhoz,         SIZE(hm_rhoz),         'hmroz',&
                     hm_rhoy,         SIZE(hm_rhoy),         'hmroy',&
                     hm_rhox,         SIZE(hm_rhox),         'hmrox',&
                     hm_vol,          SIZE(hm_vol),          'hmvol',&
                     hm_etadot,       SIZE(hm_etadot),       'hm_et',&
                     hm_theta,        SIZE(hm_theta),        'hm_th',&
                     hm_w,            SIZE(hm_w),            'hm_w_',&
                     hm_v,            SIZE(hm_v),            'hm_v_',&
                     hm_u,            SIZE(hm_u),            'h_mu_',&
                     ec_area,         SIZE(ec_area),         'ecare',&
                     ec_vol,          SIZE(ec_vol),          'ecvol')

END SUBROUTINE eg_helmholtz_check_integrity

END MODULE eg_helmholtz_mod
