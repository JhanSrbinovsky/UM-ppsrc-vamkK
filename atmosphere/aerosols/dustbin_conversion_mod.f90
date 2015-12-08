! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE dustbin_conversion_mod

!
! Description:
!   A module containing constants and subroutines involving converting
!   input LBC fields between 2- and 6-bin
!  
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Aerosols
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards version 8.4
!

IMPLICIT NONE

! Extinction proportions
REAL, PARAMETER :: kext6(1:6) = (/ 652.797, 3626.70, 976.290, 260.675,  &
                                   77.6338, 23.9075 /)

REAL, PARAMETER :: kext2(1:2) = (/ 700.367, 141.453 /)



  CONTAINS


! Convert six-bin dust to two-bin dust
  SUBROUTINE convert_dust_six_to_two(                                   &
        data_size,                                                      &
        dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6 )

    USE parkind1, ONLY: jpim, jprb 
    USE yomhook,  ONLY: lhook, dr_hook


    IMPLICIT NONE
    
! Array bounds
    INTEGER, INTENT(IN) :: data_size

! Dust data fields
    REAL, INTENT(INOUT) ::                                              &
        dust_div1(data_size),                                           &
        dust_div2(data_size),                                           &
        dust_div3(data_size),                                           &
        dust_div4(data_size),                                           &
        dust_div5(data_size),                                           &
        dust_div6(data_size)                                            

! Internal variables - mass proportions
    REAL :: p1(1:6)
    REAL :: p2(1:6)

! Temporary variables
    REAL :: dust_div1_t
    REAL :: dust_div2_t
    INTEGER :: i

! Dr Hook variables
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('CONVERT_DUST_SIX_TO_TWO',zhook_in,zhook_handle)

    ! Set p1,p2 here as setting these in the declaration enables the 'save'
    ! property by default. Which isn't desirable.
    p1(1:6) = (/ 0., 1., 1., 0.43, 0., 0. /)
    p2(1:6) = (/ 0., 0., 0., 0.57, 1., 1. /)

! Because we want to preserve AOD at 550nm rather than total mass
! the proportions must be adjusted to account for extinction
! efficiency
    p1(:) = p1(:) * kext6(:) / kext2(1)
    p2(:) = p2(:) * kext6(:) / kext2(2)
    
! Generate 2-bin dust into temporary arrays from 6-bin arrays
    DO i = 1, data_size
      dust_div1_t = dust_div1(i) * p1(1) +                         &
                    dust_div2(i) * p1(2) +                         &
                    dust_div3(i) * p1(3) +                         &
                    dust_div4(i) * p1(4) +                         &
                    dust_div5(i) * p1(5) +                         &  
                    dust_div6(i) * p1(6) 

      dust_div2_t = dust_div1(i) * p2(1) +                         &
                    dust_div2(i) * p2(2) +                         &
                    dust_div3(i) * p2(3) +                         &
                    dust_div4(i) * p2(4) +                         &
                    dust_div5(i) * p2(5) +                         &
                    dust_div6(i) * p2(6)

! Copy temporary arrays and set other bins to zero
      dust_div1(i) = dust_div1_t
      dust_div2(i) = dust_div2_t
      dust_div3(i) = 0.0
      dust_div4(i) = 0.0
      dust_div5(i) = 0.0
      dust_div6(i) = 0.0
    END DO
  IF (lhook) CALL dr_hook('CONVERT_DUST_SIX_TO_TWO',zhook_out,zhook_handle)
  RETURN
  
  END SUBROUTINE convert_dust_six_to_two


! Convert two-bin dust to six-bin dust                         
  SUBROUTINE convert_dust_two_to_six(                                   &
        data_size,                                                      &
        dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6 )

    USE parkind1, ONLY: jpim, jprb 
    USE yomhook,  ONLY: lhook, dr_hook

    IMPLICIT NONE

! Array bounds
    INTEGER, INTENT(IN) :: data_size

! Dust data fields
    REAL, INTENT(INOUT) ::                                              &
        dust_div1(data_size),                                           &
        dust_div2(data_size),                                           &
        dust_div3(data_size),                                           &
        dust_div4(data_size),                                           &
        dust_div5(data_size),                                           &
        dust_div6(data_size)                                            

! Internal variables - mass proportions
    REAL :: p1(1:6)
    REAL :: p2(1:6)

! Temporary variables
    REAL :: dust_div1_t
    REAL :: dust_div2_t
    INTEGER :: i
    
! Dr Hook variables
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('CONVERT_DUST_TWO_TO_SIX',zhook_in,zhook_handle)

    ! Set p1,p2 here as setting these in the declaration enables the 'save'
    ! property by default. Which isn't desirable.
    p1(1:6) = (/ 0.0035, 0.035, 0.220, 0.745, 0.,     0. /)
    p2(1:6) = (/ 0.,     0.,    0.,    0.219, 0.781,  0. /)

! Because we want to preserve AOD at 550nm rather than total mass
! the proportions must be adjusted to account for extinction
! efficiency
    p1(:) = p1(:) * kext2(1) / SUM( kext6(:) * p1(:) )
    p2(:) = p2(:) * kext2(2) / SUM( kext6(:) * p2(:) )

    DO i = 1, data_size
! Copy original 2-bin dust into temporary arrays
      dust_div1_t = dust_div1(i)
      dust_div2_t = dust_div2(i)
    
! Interpolate new 6-bin from 2-bin scheme
      dust_div1(i) = p1(1) * dust_div1_t + p2(1) * dust_div2_t
      dust_div2(i) = p1(2) * dust_div1_t + p2(2) * dust_div2_t
      dust_div3(i) = p1(3) * dust_div1_t + p2(3) * dust_div2_t
      dust_div4(i) = p1(4) * dust_div1_t + p2(4) * dust_div2_t
      dust_div5(i) = p1(5) * dust_div1_t + p2(5) * dust_div2_t
      dust_div6(i) = p1(6) * dust_div1_t + p2(6) * dust_div2_t
    END DO
      
    IF (lhook) CALL dr_hook('CONVERT_DUST_TWO_TO_SIX',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE convert_dust_two_to_six


END MODULE dustbin_conversion_mod
