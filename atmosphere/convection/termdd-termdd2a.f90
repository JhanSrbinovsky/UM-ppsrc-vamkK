! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate whether downdraught is able to continue
!
! Subroutine Interface: 
!
SUBROUTINE termdd (npnts, k, bdd_start                                     &
                   , thdd_k, qdd_k, the_k, qe_k, ppn_mix_dd                &
                   , b_dd_end, bdd_on)

USE cv_run_mod, ONLY:                                                      &
    dd_opt

USE atmos_constants_mod, ONLY: c_virtual

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
! 
! Description: Calculate whether downdraught is able to continue
!              Calculate buoyancy
!
! Method: UM documentataion paper 27
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER, INTENT(IN) :: &
  npnts                &  ! Vector length
 ,k                       ! Present model layer

LOGICAL, INTENT(IN) :: &
  bdd_start(npnts)        ! Mask for those points where downdraught may occur
                          ! in layer k-1

REAL, INTENT(IN) :: &
  thdd_k(npnts)        & ! Potential temperature of downdraught in layer k (K)

 ,qdd_k(npnts)         & ! Mixing ratio of downdraught in layer k  (kg/kg) 

 ,the_k(npnts)         & ! Potential temperature of environment in layer k (K)

 ,qe_k(npnts)          & ! Mixing ratio of environment in layer k   (kg/kg) 

 ,ppn_mix_dd(npnts)      ! precipitation mixing ratio (kg/kg)      
 
LOGICAL, INTENT(OUT) :: &
  b_dd_end(npnts)       & ! Mask for those points where downdraught is 
                          ! terminating
 ,bdd_on(npnts)           ! Mask for those points where downdraught continues
                          ! layer k-1 (as bdd_start here)

!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::   &
  i             ! Loop counters

REAL ::      &
  buoy1      &  ! Buoyancy of parcel

 ,thdd_v     &  ! Used in calculation of buoyancy

 ,the_v         ! Used in calculation of buoyancy


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('TERMDD',zhook_in,zhook_handle)

IF (dd_opt == 1) THEN
  !-----------------------------------------------------------------------
  ! Check if parcel still negatively buoyant such that downdraught
  ! can continue to next layer
  !-----------------------------------------------------------------------

  DO i=1,npnts
    thdd_v = thdd_k(i)*(1.0+c_virtual*qdd_k(i)) /(1.+ppn_mix_dd(i))
    the_v  = the_k(i) *(1.0+c_virtual*qe_k(i))
    buoy1  = thdd_v - the_v

   !-----------------------------------------------------------------------
   ! Calculate state of downdraught
   !-----------------------------------------------------------------------

    IF (bdd_start(i) .AND. buoy1 >  0.5) THEN
      bdd_on(i) = .FALSE.
    ELSE IF (buoy1 >  0.5 .OR. k == 2) THEN
      b_dd_end(i) = .TRUE.
    END IF
  END DO

ELSE
  !-----------------------------------------------------------------------
  ! Check if parcel still negatively buoyant such that downdraught
  ! can continue to next layer
  !-----------------------------------------------------------------------

  DO i=1,npnts
    thdd_v = thdd_k(i)*(1.0+c_virtual*qdd_k(i))
    the_v  = the_k(i)* (1.0+c_virtual*qe_k(i))
    buoy1  = thdd_v - the_v

    !-----------------------------------------------------------------------
    ! Calculate state of downdraught
    !-----------------------------------------------------------------------

    IF (bdd_start(i) .AND. buoy1 >  0.5) THEN
      bdd_on(i) = .FALSE.
    ELSE IF (buoy1 >  0.5 .OR. k == 2) THEN
      b_dd_end(i) = .TRUE.
    END IF
  END DO

END IF         ! dd_opt test

IF (lhook) CALL dr_hook('TERMDD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE termdd

