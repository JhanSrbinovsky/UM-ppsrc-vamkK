! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE INITQLCF------------------------------------------------
!
!    Purpose: Diagnose (initialise) cloud water (frozen + liquid
!             separately) and cloud amount, from standard prognostic
!             model variables (Q, T, P*).
!
!        For initialisation of RHC see the documentation
!        of sub-component P292 (large-scale cloud). SRHC1 and SRHC2
!        contain the critical relative humidity discussed in the
!        paragraph incorporating equations P292.11 - P292.14.  The
!        values below are based on those used in the old GCM.
!
!   Logical components covered : P292
!
!    Arguments:---------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Single Column Model

SUBROUTINE initqlcf                                                           &
  ( p, rhcrit, q, t, p_levels, row_length, rows, ocf, oqcf, oqcl, nbl, jlev )

USE atmos_constants_mod, ONLY: repsilon, r, cp


  USE water_constants_mod, ONLY: lc, lf, tm
USE conversions_mod, ONLY: pi

  IMPLICIT NONE

  INTEGER ::    &
    nbl         &! In No of levels in b.l.
  , jlev        &! In Model level to be processed
  , row_length  &! In No. of points in x direction.
  , rows        &! In No. of points in y direction.
  , p_levels     ! In No. of hybrid levels in the model.

  REAL ::                  &
    rhcrit(p_levels)       &!   Critical relative humidities
  , p(row_length,rows)     &! In Pressure (Pa).
  , q(row_length,rows)     &! In Sp humidity (kg water per kg air).
  , t(row_length,rows)     &! In Temperature (K).
  , ocf(row_length,rows)   &! Out Cloud fraction (decimal fraction).
  , oqcf(row_length,rows)  &! Out Cloud ice content (kg per kg air).
  , oqcl(row_length,rows)   ! Out Cloud liquid water (kg per kg air).

! Workspace usage:-----------------------------------------------------
  REAL :: &         ! Workspace (see later comments for usage):-
    wqsat(row_length,rows)  ! WORK

! ---------------------------------------------------------------------
! Local varables:------------------------------------------------------
  REAL ::      &
    x          &! Local (see later comments for usage):-
  , wac        &! Local
  , wrh        &! Local
  , alphalrcp  &! Local
  , qc          ! Local Cloud water/ice content


! Derived + local constants.
  REAL, PARAMETER ::                    &
    pfiddle  = 0.25                     &! Used in 8/8 cloud case.
  , palcon2e = lc*lc*repsilon           &
  , paldep2e = (lc+lf)*(lc+lf)*repsilon &
  , prcp     = r*cp                     &
  , pc1      = 1.060660172              &! 3/sqrt(8).
  , pc2      = 2.0*pc1                  &
  , pc3      = pi/3.0                    ! pi/3

  REAL :: &
    srhc1 &! Critical relative humidity consts
  , srhc2

  DATA            &
    srhc1/0.925/  &  ! In b.l.
  , srhc2/0.85/      ! Above b.l.

  REAL :: &
   rhc     ! Critical relative humidity

  INTEGER :: &
    i, j      ! Do loop index

!-----------------------------------------------------------------------
!   0. Initialise critical relative humidity according to level
!-----------------------------------------------------------------------

  IF (jlev <  nbl) THEN
    rhc = srhc1
  ELSE
    rhc = srhc2
  END IF

  rhc = rhcrit(jlev)

!-----------------------------------------------------------------------
!   1. Calculate pressure (in array w), hence qsat, hence relative
!      humidity in WRH.
!-----------------------------------------------------------------------

! DEPENDS ON: qsat
  CALL qsat(wqsat,t(1,1),p,row_length*rows)

  DO j=1, rows
    DO i=1, row_length
      wrh = q(i,j)/wqsat(i,j)
      IF (wrh >  1.0) THEN
        wrh = 1.0
      END IF

!-----------------------------------------------------------------------
!   2. Calculate cloud fraction ocf.  (Known as C in formulae).
!-----------------------------------------------------------------------
      ocf(i,j) = 0.0

      IF (rhc < 1.0) THEN

        IF (wrh > rhc .AND. wrh < (5.0+rhc)/6.0) THEN
          ocf(i,j) = 2.0*COS(pc3+ACOS(pc1*(wrh-rhc)/(1.0-rhc) )/3.0)
          ocf(i,j) = ocf(i,j)*ocf(i,j)
        END IF

        IF (wrh >= (5.0+rhc)/6.0) THEN
          ocf(i,j) = pc2*(1.0-wrh)/(1.-rhc)
          ocf(i,j) = 1.0-ocf(i,j)**(2.0/3.0)
        END IF

        IF (ocf(i,j) < 0.0) ocf(i,j) = 0.0
        IF (ocf(i,j) > 1.0) ocf(i,j) = 1.0

      ELSE ! rhc = 1. so set cloud fraction to 0 or 1

        IF (wrh < 1.0) THEN
          ocf(i,j) = 0.0
        ELSE
          ocf(i,j) = 1.0
        END IF

      END IF

!-----------------------------------------------------------------------
!   3. Calculate F(C) - store in wac.
!-----------------------------------------------------------------------
      wac = 0.0

      IF (ocf(i,j) <= 0.5 .AND. ocf(i,j) > 0.0)THEN
        wac = 2.0*ocf(i,j)
        wac = SQRT(wac*wac*wac)/6.0
      END IF

      IF (ocf(i,j) > 0.5)THEN
        wac = 2.0*(1.0-ocf(i,j))
        wac = 1.0+                                                     &
               SQRT(wac*wac*wac)/6.0-                                  &
               SQRT(wac)
      END IF

!-----------------------------------------------------------------------
!   4. Calculate A(C) - store in wac.
!-----------------------------------------------------------------------
      wac = wac*(1.0-rhc)

!-----------------------------------------------------------------------
!   5. Calculate total cloud water - store in qc
!-----------------------------------------------------------------------
      IF (t(i,j) > tm) THEN
        alphalrcp = palcon2e*wqsat(i,j)/(prcp*t(i,j)*t(i,j))
      ELSE
        alphalrcp = paldep2e*wqsat(i,j)/(prcp*t(i,j)*t(i,j))
      END IF

! Special treatment of full liquid cloud cover case.
      IF (ocf(i,j) >= 1.0 .AND. t(i,j) >= tm) THEN
        qc = pfiddle/(1.0+alphalrcp*(1.0+pfiddle))
      ELSE
        qc = wac/(1.0+alphalrcp*(wrh+wac))
      END IF
      qc = qc*wqsat(i,j)
!-----------------------------------------------------------------------
!   6. Partition cloud water into liquid and ice, and store in output
!      arrays OQCL, OQCF.
!-----------------------------------------------------------------------
      oqcl(i,j) = 0.0
      oqcf(i,j) = 0.0

      IF (t(i,j) > tm) THEN
        oqcl(i,j) = qc
      ELSE
        oqcf(i,j) = qc
      END IF

    END DO ! i
  END DO ! j


  RETURN
END SUBROUTINE initqlcf

