! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate Clear Air Turbulence diagnostics.

SUBROUTINE DiffP( NumLevs,      &  ! in
                  PRef,         &  ! in
                  FFields,      &  ! in
                  PFields,      &  ! in
                  dFdP,         &  ! inout
                  ErrorStatus )    ! inout

! Description:
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs                      ! No. levels
REAL, INTENT(IN) :: PRef                            ! P where deriv reqd
TYPE(PP_Field_type), INTENT(IN) :: FFields(NumLevs) ! Input variable, F
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs) ! Corresp. P vals

TYPE(PP_Field_type), INTENT(INOUT) :: dFdP          ! Deriv wrt p
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DiffP"
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
INTEGER, PARAMETER :: NPts = 6                      ! No. spline coeffs

! Local Variables:
INTEGER :: i, j, jlwr, jupr, jmid, k, halfn
REAL :: dx, t
REAL :: b(NPts), c(NPts), d(NPts)                       ! Spline coeffs
REAL :: deltay(NPts-1), deltax(NPts-1), dydx(NPts-1)    ! Spline inputs

INTEGER, ALLOCATABLE :: Lower (:,:)    ! Lower model level
INTEGER, ALLOCATABLE :: Upper (:,:)    ! Upper model level

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( dFdP % RData ) ) THEN
  DEALLOCATE( dFdP % RData )
END IF
dFdP % Hdr = FFields(1) % Hdr
ALLOCATE( dFdP % RData(dFdP % Hdr % NumCols, &
                       dFdP % Hdr % NumRows) )
dFdP % Hdr % BMDI = RMDI
dFdP % Hdr % STCode = IMDI

ALLOCATE ( Lower ( FFields(1) % Hdr % NumCols,   &
                   FFields(1) % Hdr % NumRows ) )

ALLOCATE ( Upper ( FFields(1) % Hdr % NumCols,   &
                   FFields(1) % Hdr % NumRows ) )

halfn = NPts/2
DO j = 1,FFields(1) % Hdr % NumRows
  DO i = 1,FFields(1) % Hdr % NumCols
    ! Divide and conquer to find level containing PRef
    jlwr = halfn
    jupr = NumLevs+1-halfn
    DO WHILE ( (jupr-jlwr) > 1 )
      jmid = (jlwr+jupr)/2
      IF ( PRef >= PFields(jmid) % RData(i,j) ) THEN
        jupr = jmid
      END IF
      IF ( PRef <  PFields(jmid) % RData(i,j) ) THEN
        jlwr = jmid
      END IF
    END DO

    Lower (i,j) = jlwr
    Upper (i,j) = jupr

  ENDDO
ENDDO

DO j = 1,FFields(1) % Hdr % NumRows
   DO i = 1,FFields(1) % Hdr % NumCols
    jlwr = Lower (i,j)
    jupr = Upper (i,j)

    ! Use fields in reverse order so that pressure is increasing
    dx = PRef - PFields(jupr) % RData(i,j)
    DO k = 1,NPts-1
      deltax(k) = PFields(jupr+halfn  -k) % RData(i,j)-  &
                  PFields(jupr+halfn-1-k) % RData(i,j)
      deltay(k) = FFields(jupr+halfn  -k) % RData(i,j)-  &
                  FFields(jupr+halfn-1-k) % RData(i,j)
    END DO

!-----------------------------------------------------------------------
! the coefficients b(i), c(i), and d(i), i=1,2,...,NPts are computed
! for a cubic interpolating spline
! s(u) = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!   derivative= b(i) + 2.*c(i)*(u-x(i)) + 3.*d(i)*(u-x(i))**2
!
! where  x(i) < u < x(i+1), using horner's rule
! if  u < x(1) then  i = 1  is used.
! if  u >= x(n) then  i = n  is used.
!
! A binary search is performed to determine the proper interval.
! input..
!   NPts = the number of data points or knots (NPts.ge.2)
!   x = the abscissas of the knots in strictly increasing order
!   y = the ordinates of the knots
! output..
!   b, c, d  = arrays of spline coefficients as defined above.
! using  p  to denote differentiation,
!   y(i) = s(x(i))
!   b(i) = sp(x(i))
!   c(i) = spp(x(i))/2
!   d(i) = sppp(x(i))/6  (derivative from the right)
!-----------------------------------------------------------------------

! Tridiagonal system: b = diagonal, d = offdiagonal, c = RH side.
! Third derivs at  x(1)  and  x(n) obtained from divided differences
    ! Derivative to the right
    dydx(1:NPts-1) = deltay(1:NPts-1) / deltax(1:NPts-1)
    ! Difference in right and left derivative
    c(2:NPts-1) = dydx  (2:NPts-1) - dydx  (1:NPts-2)
    ! b(k)=x(k+1)-x(k-1)
    b(1)        = -deltax(1)
    b(NPts)     = -deltax(NPts-1)
    b(2:NPts-1) = deltax(1:NPts-2) + deltax(2:NPts-1)

    c(1)    = c(3)     /b(3)      - c(2)     /b(2)
    c(NPts) = c(NPts-2)/b(NPts-2) - c(NPts-1)/b(NPts-1)

    c(1)    = c(1)   *b(1)   **2/( b(3)     -b(1) )
    c(NPts) = c(NPts)*b(NPts)**2/( b(NPts-2)-b(NPts) )

    b(2:NPts-1) = 2.*b(2:NPts-1)
!----------------------------------------------  forward elimination
    DO k = 1, NPts-1
      t = deltax(k)/b(k)
      b(k+1) = b(k+1) - t*deltax(k)
      c(k+1) = c(k+1) - t*c(k)
    ENDDO
    c(NPts) = c(NPts)/b(NPts)
!------------------------------------------------  back substitution
    DO k = NPts-1, halfn, -1
      c(k) = (c(k) - deltax(k)*c(k+1))/b(k) !2nd derivative / 6
    ENDDO

!--------------------------------------------------------------------
    b(halfn) = dydx(halfn) - deltax(halfn)*(c(halfn+1) + 2.*c(halfn))
    d(halfn) = (c(halfn+1) - c(halfn)) /deltax(halfn) ! = 3rd deriv/6
    c(halfn) = 3.*c(halfn)                            ! = 2nd deriv/2

!-------------------------------------------------------------------
!--- evaluate the derivative of the cubic spline function
!-----------------------------------------------------------------

    dFdP % RData(i,j) = b(halfn) + dx*(2.0*c(halfn) + 3.0*dx*d(halfn))

  ENDDO
ENDDO

DEALLOCATE ( Lower )
DEALLOCATE ( Upper )

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE DiffP

