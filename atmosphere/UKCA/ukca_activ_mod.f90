  MODULE UKCA_ACTIV_MOD
! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:  Contains subroutines for generating arrays for pdf of updraught
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge
!  University of Leeds, University of Oxford, and The Met Office.
!  See: www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
!  Author: Rosalind West, AOPP, Oxford, 2010
!  -------

    
    USE yomhook,  ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    IMPLICIT NONE

  CONTAINS
  
    SUBROUTINE ACTIVMKLIN(kbdim, klev, arrmin, arrmax, nbins, binwidth, &
                          linarr)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN)  :: kbdim
      INTEGER, INTENT(IN)  :: klev
      INTEGER, INTENT(IN)  :: nbins
      REAL,    INTENT(IN)  :: arrmin(kbdim, klev)
      REAL,    INTENT(IN)  :: arrmax(kbdim, klev)
      REAL,    INTENT(OUT) :: binwidth(kbdim, klev, nbins)
      REAL,    INTENT(OUT) :: linarr(kbdim, klev, nbins)

      INTEGER :: i,j,k
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_ACTIV_MOD:ACTIVMKLIN',zhook_in,   &
                              zhook_handle)
      
      DO k = 1, nbins
        DO i=1, klev
          DO j=1, kbdim
               binwidth(j,i,k)=(arrmax(j,i)-arrmin(j,i))/REAL(nbins)
               linarr(j,i,k)= arrmin(j,i) + (REAL(k)-0.5)*binwidth(j,i,k)
          END DO
        END DO
      END DO

      IF (lhook) CALL dr_hook('UKCA_ACTIV_MOD:ACTIVMKLIN',zhook_out,    &
                              zhook_handle)
      RETURN
      END SUBROUTINE ACTIVMKLIN
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      SUBROUTINE ACTIVMKPDF (kbdim, klev, array, nbins, sigma, mean,    &
                             pdfarr)
 
      USE UKCA_CONSTANTS, ONLY: PPI
    
      IMPLICIT NONE
    
      INTEGER, INTENT(IN) :: kbdim
      INTEGER, INTENT(IN) :: klev
      INTEGER, INTENT(IN) :: nbins
      REAL,    INTENT(IN) :: array(kbdim, klev, nbins)
      REAL,    INTENT(IN) :: sigma(kbdim, klev)
      REAL,    INTENT(IN) :: mean(kbdim, klev)
      REAL,    INTENT(OUT) :: pdfarr(kbdim, klev, nbins)

      INTEGER :: i, j, k
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_ACTIV_MOD:ACTIVMKPDF',zhook_in,     &
                               zhook_handle)
 
      DO k=1, nbins
        DO i=1, klev
          DO j=1, kbdim
               ! normalised gaussian
            pdfarr(j,i,k) = (1.0/((2.0*ppi)**0.5))*(1.0/sigma(j,i))*    &
                  EXP(-((array(j,i,k)-mean(j,i))**2/(2*sigma(j,i)**2)))
            END DO
         END DO
      END DO

      IF (lhook) CALL dr_hook('UKCA_ACTIV_MOD:ACTIVMKPDF',zhook_out,    &
                               zhook_handle)

      RETURN
      END SUBROUTINE ACTIVMKPDF
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      SUBROUTINE ACTIVMKSKEW (kbdim, klev, array, nbins, sigma, mean,   &
                              alpha, pdfarr)
 
      USE UKCA_CONSTANTS, ONLY: PPI
    
      IMPLICIT NONE
    
      INTEGER, INTENT(IN) :: kbdim
      INTEGER, INTENT(IN) :: klev
      INTEGER, INTENT(IN) :: nbins
      REAL,    INTENT(IN) :: array(kbdim, klev, nbins)
      REAL,    INTENT(IN) :: sigma(kbdim, klev)
      REAL,    INTENT(IN) :: mean(kbdim, klev)
      REAL,    INTENT(IN) :: alpha(kbdim, klev)
      REAL,    INTENT(OUT) :: pdfarr(kbdim, klev, nbins)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      INTEGER :: i, j, k
 
      IF (lhook) CALL dr_hook('UKCA_ACTIV_MOD:ACTIVMKSKEW',zhook_in,   &
                               zhook_handle)
      DO i=1, klev
        DO j=1, kbdim
          DO k=1, nbins
           ! skew-normal distribution
            pdfarr(j,i,k)=(1.0/((2.0*ppi)**0.5))*(1.0/sigma(j,i))*     &
                EXP(-((array(j,i,k)-mean(j,i))**2/(2*sigma(j,i)**2)))* &
                (1+ERF(alpha(j,i)*(array(j,i,k)-mean(j,i))*            &
                (1.0/((2.0)**0.5))*(1.0/sigma(j,i))))
            END DO
         END DO
      END DO

      IF (lhook) CALL dr_hook('UKCA_ACTIV_MOD:ACTIVMKSKEW',zhook_out,   &
                               zhook_handle)

      RETURN
      END SUBROUTINE activmkskew
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      SUBROUTINE ACTIVCLOSEST (arrayin, nbins, value, closeval,         &
                               closeind)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: nbins
      REAL,    INTENT(IN)  :: arrayin(nbins)
      REAL,    INTENT(IN)  :: value
      INTEGER, INTENT(OUT) :: closeind
      REAL,    INTENT(OUT) :: closeval 

      INTEGER :: min_subs(1)
      REAL    :: subarr(nbins)
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_ACTIV_MOD:ACTIVCLOSEST',zhook_in,   &
                              zhook_handle)

    ! Subtract the value you wish to compare from
    ! the whole array you are comparing it to and 
    ! create a new array of these differences.
      subarr(:) = arrayin(:) - value

    ! Find the smallest absolute difference between value and array.
    ! *min_subs contains the location of the smallest absolute value of subarr
    !  NB the value of MINLOC must be returned to an array. 
    ! E.g. "icmin=minloc(c(1:3))" will only work if "icmin" has been 
    ! declared to be an integer array, and the value "3" will be placed 
    ! in the first element of "icmin".
    
      min_subs = MINLOC(ABS(subarr(1:nbins)))

    ! *closeind is just an integer of the 1st element of min_subs
      closeind = min_subs(1)
    
    ! *closeval contains the value of the closest element of arrayin to value
    ! *closest contains the location of that value
      closeval = arrayin(closeind)

      IF (lhook) CALL dr_hook('UKCA_ACTIV_MOD:ACTIVCLOSEST',zhook_out,   &
                              zhook_handle)
  
      RETURN
      END SUBROUTINE ACTIVCLOSEST


      END MODULE UKCA_ACTIV_MOD
