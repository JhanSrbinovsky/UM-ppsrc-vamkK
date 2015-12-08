! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE update_dpsidt_mod

IMPLICIT NONE

CONTAINS


 SUBROUTINE update_dpsidt(alpha, dpsidtc, dpsidts, n1, n2, nlim,        &
                          gspect, zero, icode, info, timestep_number)

! This routine updates the spherical harmonic coefficients using
! a first-order auto-regressive process. gspect is the wave-number
! dependent power spectrum generated at the start of the run in
! subroutine "backscatter_spectrum"

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim
 USE UM_ParVars
 IMPLICIT NONE

! Parameters and common blocks required by the MPP-UM

 INTEGER, INTENT(IN) :: n1, n2, nlim, icode, info, zero, timestep_number
 REAL, INTENT(IN) :: alpha, gspect(nlim)
 REAL, INTENT(OUT) :: dpsidtc(0:n2,n1:n2), dpsidts(0:n2,n1:n2)

 !     Local variables
 REAL, SAVE :: oneminusalpha, squarerootalpha, sdx2
 INTEGER :: m, n
 ! MPP quantities
 INTEGER :: sndcount

 ! Random numbers generated on processor=0 (only needed if REPROD = T)
 REAL, DIMENSION(:,:), ALLOCATABLE :: z_rand_nums, z_rand_nums2
 ! Random numbers on each processor
 REAL :: rand_nums(0:n2,n1:n2), rand_nums2(0:n2,n1:n2)

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('UPDATE_DPSIDT',zhook_in,zhook_handle)

 IF ( mype == zero ) THEN
   ! Generate a separate set of random numbers for sin and cosine coeffs
   CALL random_number(rand_nums)
   CALL random_number(rand_nums2)
 END IF

 ! Scatter Random numbers to all processors
 sndcount = (n2+1)*(n2-n1+1)
 CALL gc_rbcast(3245, sndcount, zero, nproc, icode, rand_nums)
 CALL gc_rbcast(3246, sndcount, zero, nproc, icode, rand_nums2)

 IF (timestep_number == 1) THEN

! Set up AR1-process parameters and save
   oneminusalpha=1-alpha
! SQRT(alpha) ensures that the noise (random nums) decorrelates
! faster than the timestep
   squarerootalpha=SQRT(alpha)
! Twice standard deviation of random numbers
! These have a variance of (1/12)
   sdx2 = 2.0 * SQRT(0.083333333333333)

   DO n = n1, n2
     DO m = 0, n
! Initialise with random value around expected mean
       dpsidtc(m,n) = sdx2 * (rand_nums(m,n)-0.5)*gspect(n)
       dpsidts(m,n) = sdx2 * (rand_nums2(m,n)-0.5)*gspect(n)
     END DO
   END DO
 ELSE
   DO n = n1, n2
     DO m = 0, n
! (RANDOM-0.5) provides white noise with mean of 0 and var. of (1/12)
       dpsidtc(m,n) = oneminusalpha*dpsidtc(m,n) + squarerootalpha      &
                      *(rand_nums(m,n)-0.5)*gspect(n)
       dpsidts(m,n) = oneminusalpha*dpsidts(m,n) + squarerootalpha      &
                      *(rand_nums2(m,n)-0.5)*gspect(n)
     END DO
   END DO
 ENDIF

 IF (lhook) CALL dr_hook('UPDATE_DPSIDT',zhook_out,zhook_handle)
 RETURN

END SUBROUTINE update_dpsidt

END MODULE update_dpsidt_mod
