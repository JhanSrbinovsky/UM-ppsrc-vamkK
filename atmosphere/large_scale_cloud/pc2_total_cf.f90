! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ PC2 Cloud Scheme: Change in total cloud fraction
! Subroutine Interface:
SUBROUTINE pc2_total_cf(                                                &
!      Number of points
 points,                                                                &
!      Input fields
 cfl, cff, deltacl, deltacf,                                            &
!      Input Output fields
 cf)

  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim

  IMPLICIT NONE

! Purpose:
!   Update the total cloud fraction due to changes in
!   liquid and ice cloud fractions.
!
! Method:
!   Assumes a random overlap of the forced cloud with already existing
!   conditions in the gridbox. See Annex D of the PC2 cloud scheme
!   documentation.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation

!  Subroutine Arguments:------------------------------------------------
  INTEGER ::                                                            &
                        !, INTENT(IN)
   points
!       No. of points being processed.

  REAL ::                                                               &
                        !, INTENT(IN)
   cfl(points),                                                         &
!       Liquid cloud fraction
   cff(points),                                                         &
!       Ice cloud fraction
   deltacl(points),                                                     &
!       Change in liquid cloud fraction
   deltacf(points)
!       Change in ice cloud fraction

  REAL ::                                                               &
                        !, INTENT(INOUT)
   cf(points)
!       Total cloud fraction

!  External functions:

!  Local scalars--------------------------------------------------------

!  (a)  Scalars effectively expanded to workspace by the Cray (using
!       vector registers).
  REAL ::                                                               &
   deltac_1, & ! Change in total cloud fraction due to liquid cloud
   deltac_2, & ! Change in total cloud fraction due to ice cloud
   cf_1,     & ! Input cf + deltac_1 (Used to avoid rounding error)
   cf_2        ! Input cf + deltac_2 (Used to avoid rounding error)

!  (b)  Others
  INTEGER :: i ! Loop counter

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

!- End of Header

  IF (lhook) CALL dr_hook('PC2_TOTAL_CF',zhook_in,zhook_handle)

! ==Main Block==--------------------------------------------------------
! Points_do1:
  DO i = 1, points

! ----------------------------------------------------------------------
! 1. Update total cloud fraction.
! ----------------------------------------------------------------------

! Calculate change in total cloud fraction due to a change in liquid
! cloud fraction. This depends upon the sign of the change of liquid
! cloud fraction.

    IF (deltacl(i) > 0.0 .AND. cfl(i) < 1.0) THEN

!     Random overlap
!     deltac_1 = deltacl(i) *(1.0 - cf(i))/(1.0 - cfl(i))

!     Minimum overlap
      deltac_1 = MIN(         deltacl(i), (1.0-cf(i)) )
      cf_1     = MIN( cf(i) + deltacl(i),  1.0        )

    ELSE IF (deltacl(i) < 0.0 .AND. cfl(i) > 0.0) THEN

!     Random overlap
!     deltac_1 = deltacl(i) * (cf(i)-cff(i)) / cfl(i)

!     Minimum overlap
      deltac_1 = MAX(         deltacl(i), (cff(i)-cf(i)) )
      cf_1     = MAX( cf(i) + deltacl(i),  cff(i)        )

    ELSE
      deltac_1 = 0.0
      cf_1     = cf(i)
    END IF

! Calculate change in total cloud fraction due to a change in ice
! cloud fraction. This depends upon the sign of the change of ice
! cloud fraction.

    IF (deltacf(i) > 0.0 .AND. cff(i) < 1.0) THEN

!     Random overlap
!     deltac_2 = deltacf(i) *(1.0 - cf(i))/(1.0 - cff(i))

!     Minimum overlap
      deltac_2 = MIN(         deltacf(i), (1.0-cf(i)) )
      cf_2     = MIN( cf(i) + deltacf(i),  1.0        )

    ELSE IF (deltacf(i)  <   0.0 .AND. cff(i)  >   0.0) THEN

!     Random overlap
!     deltac_2 = deltacf(i) * (cf(i)-cfl(i)) / cff(i)

!     Minimum overlap
      deltac_2 = MAX(         deltacf(i), (cfl(i)-cf(i)) )
      cf_2     = MAX( cf(i) + deltacf(i),  cfl(i)        )

    ELSE
      deltac_2 = 0.0
      cf_2     = cf(i)
    END IF

! Sum the two changes

    cf(i) = cf(i) + deltac_1 + deltac_2

! Avoid rounding error in limiting cases

    IF (deltac_1 == 0.0) cf(i) = cf_2
    IF (deltac_2 == 0.0) cf(i) = cf_1

! For minimum overlap we need to check that the total cloud
! fraction is constrained within 0 and 1

    cf(i) = MAX(  MIN( cf(i),1.0 )  , 0.0)

! Points_do1:
  END DO

! End of the subroutine

  IF (lhook) CALL dr_hook('PC2_TOTAL_CF',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_total_cf
