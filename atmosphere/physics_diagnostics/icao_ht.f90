! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE ICAO_HT-----------------------------------------
!LL
!LL  PURPOSE:   Performs an interpolation from pressure levels to the
!LL             I.C.A.O standard atmosphere heights.
!LL             e.g. tropopause heights or max wind heights in
!LL             thousands of feet
!LL
!LL  Logical components covered D4
!LL
!LL  Project TASK: D4
!LL
!LL  Programming standard: U M DOC  Paper NO. 4,
!LL
!LL
!LLEND------------------------------------------------------------------
!
!*L  ARGUMENTS:---------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Physics Diagnostics
MODULE icao_ht_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE ICAO_HT(                                               &
     & press_in, row_length, rows, icao_height)
!-----------------------------------------------------------------------

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
! Subroutine arguments
       INTEGER, INTENT(IN) ::                                           &
     &  row_length                                                      &
                                     ! number of points on a row
     &, rows                         ! number of rows in a theta field

      REAL, INTENT(IN) ::                                               &
     & press_in(row_length,rows)     ! Pressures to be converted to
                                     ! ICAO heights

      REAL, INTENT(OUT) ::                                              &
     & icao_height(row_length,rows)  ! ICAO height in thousands of feet

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! External routines
! None

! Local parameters:

      INTEGER, PARAMETER ::                                             &
     &      phst_len = 15

      REAL, PARAMETER ::                                                &
     &      p0 = 101325.0,                                              &
     &      pr = 1000.0

! Local variables:

      INTEGER ::  i,j,k              ! Loop Counters

      REAL ::     term1, term2

      REAL ::                                                           &
     & phst(phst_len)                                                   &
     &,hst(phst_len)                                                    &
     &,pressure(row_length,rows)     ! Local pressure array

      LOGICAL ::                                                        &
     & found(row_length,rows)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Data for local arrays
      DATA phst/101325.0, 95000.0, 85000.0, 70000.0, 50000.0, 40000.0,  &
     &           30000.0, 25000.0, 20000.0, 15000.0, 10000.0,  7000.0,  &
     &            5000.0,  3000.0,   999.0/
      DATA hst/      0.0,  2000.0,  5000.0, 10000.0, 18000.0, 24000.0,  &
     &           30000.0, 34000.0, 39000.0, 45000.0, 53000.0, 61000.0,  &
     &           68000.0, 78000.0, 99000.0/

!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('ICAO_HT',zhook_in,zhook_handle)

      DO j=1,rows
        DO i=1,row_length
          pressure(i,j)=press_in(i,j)

! 1. If the pressure is less than 10mb it is assumed to be 10mb
!    If pressure is greater than 1013.25mb it's assumed to be 1013.25mb

          IF (pressure(i,j) <= pr .AND. pressure(i,j) >= 0.) THEN
            pressure(i,j) = pr
          END IF
          IF (pressure(i,j)  >   p0) THEN
            pressure(i,j) = p0
          END IF
          found(i,j)=.false.
        END DO
      END DO

!  2. Find the first value of PHST that is less than the pressure

      DO k=1,phst_len-1
        DO j=1,rows
          DO i=1,row_length

! 3. If the pressure is set to missing data indicator, set ICAO height
!    to missing data indicator

            IF (pressure(i,j) == RMDI .AND. (.NOT.found(i,j))) THEN
              icao_height(i,j) = RMDI
              found(i,j) = .true.
            ELSE IF ((phst(k+1) <  pressure(i,j)).AND.(.NOT.found(i,j)))&
     &        THEN

! 4. Calculate the ICAO height

              term1=(3*phst(k)-pressure(i,j))*(pressure(i,j)-phst(k))
              term2=(3*phst(k)-phst(k+1))*(phst(k+1)-phst(k))
              icao_height(i,j)=((hst(k+1)-hst(k))*term1/term2           &
     &                         + hst(k))*0.001
              found(i,j)=.true.
            END IF
          END DO
        END DO
      END DO

      IF (lhook) CALL dr_hook('ICAO_HT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ICAO_HT
END MODULE icao_ht_mod
