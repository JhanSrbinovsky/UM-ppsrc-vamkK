! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing tcs warm rain subroutine for calculating 
! scalings for tcs scheme
!
MODULE tcs_calc_scales


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  !
  ! Description:
  !   This routine calculates scalings for tcs warm rain calculations
  !
  ! Method:
  !   <reference to documentation to go here, once available>
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.2 programming standards.

CONTAINS

  SUBROUTINE calc_scales_warm(scales, conv_type, ntml, ntpar,      &
     wstar_dn, delthvu, theta, q_mix, wthvs, dthetav_cb, ztop_uv,  &
     zlcl_uv, z_rho)
!-----------------------------------------------------------------------
! Description:
!   Calculate scales for convection scheme - warm rain version.
!
!   Allocation for the scales array is done in the calling procedure
!-----------------------------------------------------------------------

    USE tcs_parameters_warm, only:                                 &
       mb_a1, mb_a2, wup_a1

    USE tcs_constants, only:                                       &
       g, c_virtual

    USE tcs_class_scales, only :                                   &
       scales_conv

    IMPLICIT NONE
    TYPE(scales_conv), INTENT(inout) :: scales

    INTEGER, INTENT(in) :: conv_type(:)
                  ! Indicator of type of convection 
                  !    1=non-precipitating shallow
                  !    2=drizzling shallow
                  !    3=warm congestus

    INTEGER, INTENT(in) ::  ntml(:)

    INTEGER, INTENT(in) ::  ntpar(:)                                                       

    REAL, INTENT(in) :: wstar_dn(:)
    
    REAL, INTENT(in) :: delthvu(:)
    
    REAL, INTENT(in) :: theta(:,:)
    
    REAL, INTENT(in) :: q_mix(:,:)
    
    REAL, INTENT(in) :: wthvs(:)
    
    REAL, INTENT(in) :: dthetav_cb(:)
    
    REAL, INTENT(in) :: ztop_uv(:)
    
    REAL, INTENT(in) :: zlcl_uv(:)
    
    REAL, INTENT(in) :: z_rho(:,:)
    
    !-------------------------------------------------------------------
    ! Loop counters
    !-------------------------------------------------------------------
    INTEGER :: i

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('TCS_CALC_SCALES:CALC_SCALES_WARM',zhook_in,zhook_handle)
    
    !-------------------------------------------------------------------
    ! Put wstar_dn (from the BL) into the scales structure.
    !-------------------------------------------------------------------
    scales%wstar_dn(:)=wstar_dn(:) 
    
    !-------------------------------------------------------------------
    ! Cloud base mass flux closure
    !-------------------------------------------------------------------
    
       scales%mb(:)  = mb_a1 * wstar_dn(:)
       !
       ! The offset for congestus needs to be more carefully 
       ! considered and related to evaporation of precip/Ri
       !

       DO i=1,SIZE(conv_type)
         IF (conv_type(i)>=3)THEN
           scales%mb(i)=MAX(scales%mb(i) - 0.016, .005)
         END IF
       END DO

    !-------------------------------------------------------------------
    ! Alternative cloud base mass flux
    !-------------------------------------------------------------------
       scales%mb_new(:)=scales%mb(:)

!  Not currently doing this...      
!     WHERE (ABS(wthvs) < TINY(wthvs(1))) ! Since wthvs can be zero from 
!                                         ! BL in full model
!        scales%mb_new(:) = mb_a2*wstar_dn(:)
!     ELSEWHERE
!        ! added condition on dthetav_cb to avoid problems
!        scales%mb_new(:) = mb_a2*wthvs(:)/               &
!             MERGE(dthetav_cb(:), wthvs(:)/wstar_dn(:),  &
!             dthetav_cb(:) > wthvs(:)/wstar_dn(:))
!     END WHERE

    !-------------------------------------------------------------------
    ! Ensemble vertical velocity at cloud base
    !-------------------------------------------------------------------
    scales%wup2_cb(:) = wup_a1 * wstar_dn(:) * wstar_dn(:)

    !-------------------------------------------------------------------
    ! Working note: Check this may wish to use rho and dz rather than
    ! delthuv 
    !-------------------------------------------------------------------
    DO i=1,SIZE(ntml)
       scales%wstar_up(i) = (delthvu(i) * scales%mb(i) * G &
            / (theta(i,ntml(i))                     &
            * (1.0 + C_VIRTUAL * q_mix(i,ntml(i)))))**0.3333
       scales%zlcl(i) = z_rho(i,ntml(i)+1)
       scales%zcld(i) = z_rho(i,ntpar(i)+1) - z_rho(i,ntml(i)+1)
    END DO

    scales%zcld_uv(:) = ztop_uv(:) - zlcl_uv(:)

    !-------------------------------------------------------------------
    ! pre-calculate frequently used expressions
    !-------------------------------------------------------------------
    scales%wsc_o_mb(:) = scales%wstar_up(:)/scales%mb(:)
    scales%mb_o_wsc(:) = scales%mb(:)/scales%wstar_up(:)
    scales%root_mb_o_wsc(:) = SQRT(scales%mb_o_wsc(:))
    scales%mb_new_o_wsc(:) = scales%mb_new(:)/scales%wstar_up(:)
    scales%root_mb_new_o_wsc(:) = SQRT(scales%mb_new_o_wsc(:))
    scales%wstar_up3(:)=scales%root_mb_new_o_wsc(:) *(scales%wstar_up(:)**3)
    IF (lhook) CALL dr_hook('TCS_CALC_SCALES:CALC_SCALES_WARM',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE calc_scales_warm
  
END MODULE tcs_calc_scales
