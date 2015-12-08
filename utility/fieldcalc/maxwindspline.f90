! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate Max Wind value and height
! ====================================================================

SUBROUTINE MaxWindSpline (NumLevs, i, j,                       &
                          UFields, VFields, PFields, MaxWLev,  &
                          Uinc, Vinc, Pinc )

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


USE IO_Mod,  ONLY:     &
  PP_Field_type

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
INTEGER, INTENT(IN) :: i, j       ! Row, Col on grid

TYPE(PP_Field_type), INTENT(IN) :: UFields(NumLevs) ! U-wind on B-grid
TYPE(PP_Field_type), INTENT(IN) :: VFields(NumLevs) ! V-wind on B-grid
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs) ! P on wind B-grid

INTEGER, INTENT(IN) :: MaxWLev(PFields(1)%Hdr%NumCols,  &
                               PFields(1)%Hdr%NumRows)

INTEGER, PARAMETER :: ninc=16           ! number of increments used

REAL, INTENT(OUT) :: Uinc(2*ninc)  ! u at increment points
REAL, INTENT(OUT) :: Vinc(2*ninc)  ! v at increment points
REAL, INTENT(OUT) :: Pinc(2*ninc)  ! p at increment points

! Local constants

INTEGER, PARAMETER :: kmax=5            ! total levels needed for spline
INTEGER, PARAMETER :: khalf=(kmax+1)/2  ! number of levels each side

! Local variables:
INTEGER :: k
REAL :: P_Lwr, P_Mid, P_Upr

! Variables involving spline interpolation
INTEGER :: spl(2*ninc)
REAL :: xsp(kmax), ysp(kmax)
REAL :: bsp(kmax), csp(kmax), dsp(kmax)
REAL :: dx  (2*ninc)


! Use levels above and below to create splines of U and V
! Use the splines to evaluate U and V at intervals between
! level-1 and level+1      OLD : level-1/2 and level+1/2

! P_Lwr : P on Rho Level below maximum
! P_Mid : P on Rho Level at    maximum
! P_Upr : P on Rho Level above maximum

P_Lwr = PFields(MaxWLev(i,j)-1) % RData(i,j)
P_Mid = PFields(MaxWLev(i,j)  ) % RData(i,j)
P_Upr = PFields(MaxWLev(i,j)+1) % RData(i,j)

DO k = 1, ninc                               ! Calc p increments
   Pinc(k)      = P_Lwr + (P_Mid-P_Lwr)*FLOAT(k)/FLOAT(ninc)
   Pinc(k+ninc) = P_Mid + (P_Upr-P_Mid)*FLOAT(k)/FLOAT(ninc)
END DO

!-- Pressure is x -- going down through levels so it is increasing
DO k = 1,kmax                      ! P on surrounding levels
  xsp(k) = PFields(MaxWLev(i,j)+khalf-k) % RData(i,j)
END DO

spl (1     :  ninc) = khalf
spl (1+ninc:2*ninc) = khalf-1
dx(:) = Pinc(:) - xsp(spl(:))

!-- Get U values --
DO k = 1,kmax                      ! U on surrounding levels
  ysp(k) = UFields(MaxWLev(i,j)+khalf-k) % RData(i,j)
END DO

! DEPENDS ON: splinesetup
CALL SplineSetup( kmax, xsp, ysp, bsp, csp, dsp )
Uinc(:) = ysp(spl) + dx*( bsp(spl) + dx*( csp(spl) + dx*dsp(spl)))

!-- Get V values --
DO k = 1,kmax                      ! V on surrounding levels
  ysp(k) = VFields(MaxWLev(i,j)+khalf-k) % RData(i,j)
END DO

! DEPENDS ON: splinesetup
CALL SplineSetup( kmax, xsp, ysp, bsp, csp, dsp )
Vinc(:) = ysp(spl) + dx*( bsp(spl) + dx*( csp(spl) + dx*dsp(spl)))

END SUBROUTINE MaxWindSpline


! ====================================================================
