! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Diagnostic RHcrit Calculation Scheme.
! Subroutine Interface:
SUBROUTINE ls_calc_rhcrit(                                              &
!      Pressure related fields
  p_layer_centres,                                                      &
!      Array dimensions
  levels, row_length, rows, global_row_length,                          &
!      Prognostic Fields
  t, q, qcf, land_frac, ice_frac,                                       &
!      Logical control
  l_mixing_ratio,                                                       &
!      Output
  rhcpt,                                                                &
! Input for calc height and gridbox scale
  delta_lambda)
      USE water_constants_mod, ONLY: lc, lf, tm
  USE global_2d_sums_mod, ONLY: global_2d_sums
  USE earth_constants_mod, ONLY: earth_radius
  USE atmos_constants_mod, ONLY: cp, r, repsilon
  
  USE level_heights_mod, ONLY:  r_theta_levels
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  USE UM_ParVars
  USE rhccon2b_mod, ONLY: rhc_con1, rhc_con2, rhc_con3, rhc_con4,       &
                          rhc_min, rhc_max
  IMPLICIT NONE

! Purpose: To calculate the critical relative humidity in every
!          grid-box.

! Method : The critical relative humidity of a certain grid-box is
!          determined from the variance in a 3*3 set of boxes centred
!          on it. A fit, dependent on pressure, relates the variance
!          of the 3*3 region to the variance within the one grid-box.
!          This variance is converted to a critical relative humidity
!          in a straightforward fashion.
!          Some points in the 3*3 region may be excluded from the
!          variance calculation in the BL layers, if their
!          surfaces do not 'match'. The criterion for matching is that
!          land and sea-ice match, but that open ocean does not match
!          with either of these.
!          In all layers, points in the 3*3 region which lie outside 2
!          std devs of the mean are excluded in a second iteration of
!          the main calculation.
!          The best estimate of the standard deviation in a sample
!          size of n elements in which the mean is calculated from the
!          sample uses (n-1) in the denominator. This subroutine uses
!          n in the denominator because errors in the estimate of the
!          std dev from using such small sample sizes gave rise to
!          concerns over the std dev being too large occasionally,
!          and using n in the denominator rather than (n-1) addresses
!          this issue indirectly. Addressing the issue of large
!          estimates of std dev was considered more important than
!          using the 'best' estimate of the std dev.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Code Description:
!   Language: FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

! Declarations:

!  Global Variables:----------------------------------------------------


!  Subroutine arguments
!-----------------------------------------------------------------------
! IN variables
!-----------------------------------------------------------------------
  INTEGER ::                                                            &
                        !, INTENT(IN)
    levels,                                                             &
!       No. of levels being processed.
    row_length, rows,                                                   &
!       Horizontal dimensions of arrays being processed by scheme.
    global_row_length
!       Length of full global row (ie. counting all relevant PEs).

  LOGICAL ::                                                            &
    l_mixing_ratio                    ! Use mixing ratios

  REAL ::                                                               &
                        !, INTENT(IN)
    ice_frac(0:row_length+1,0:rows+1),                                  &
!       Ice fraction.
    land_frac(0:row_length+1,0:rows+1),                                 &
    p_layer_centres(0:row_length+1,0:rows+1,0:levels),                  &
!       pressure at all points, on theta levels (Pa).
!       NB: Zero level of array is surface pressure.
    q(0:row_length+1,0:rows+1,levels),                                  &
!       Total water content, Q+QCL (kg per kg)
    qcf(0:row_length+1,0:rows+1,levels),                                &
!       Cloud ice content at processed levels (kg water per kg air).
    t(0:row_length+1,0:rows+1,levels),                                  &
!       Liquid/frozen water temperature (K)
    delta_lambda
!       EW (x) grid spacing in radians 

!-----------------------------------------------------------------------
! OUT variables
!-----------------------------------------------------------------------
  REAL ::                                                               &
                        !, INTENT(OUT)
    rhcpt(row_length,rows,levels)
!       Critical relative humidity at every gridpoint.

!  Local Parameters and other physical constants------------------------


  REAL :: inv9, ls, lsrcp, ercpr
  PARAMETER ( inv9 = 1./9.,                                             &
              ls = lc+lf,                                               &
              lsrcp = (lc+lf)/cp,                                       &
              ercpr = repsilon/(cp*r))

!  Local scalars--------------------------------------------------------

  INTEGER ::                                                            &
      i, j, k, j8,                                                      &
                            ! Simple loop variables
      COUNT,                                                            &
                            ! Counter of points
      im1,ip1, jm1,jp1,                                                 &
                            ! (I-1),(I+1),(J-1),(J+1)
      ij_field,                                                         &
                            ! Number of non-halo points in arrays
      i_length,                                                         &
                            ! Row length for polar row adjustments
      j_start,                                                          &
                            ! Row number for polar row adjustments
      istat             ! Status (error code) indicator

  REAL ::                                                               &
      mean_supsat,                                                      &
                            ! MEAN RH OF 3*3 REGION
      tot_var,                                                          &
                            ! TOTAL VARIANCE OF 3*3 REGION
      supsat_sd_1,                                                      &
                            ! STANDARD DEVIATION OF 'R.H.' IN GRID-BOX
      supsat_sd_3,                                                      &
                            ! RESOLVED STD DEV OF 'R.H.' IN 3*3 REGION
      root_6,                                                           &
                            ! =sqrt(6.)
      latht,                                                            &
                            ! =Lc if T>Tm, ELSE = Lc+Lf
      two_sigma,                                                        &
                          ! Two times sigma
      supsat1,                                                          &
                            ! Temporary variable
      supsat2,                                                          &
                            ! Temporary variable
      supsat3,                                                          &
                            ! Temporary variable
      supsat4,                                                          &
                            ! Temporary variable
      supsat5,                                                          &
                            ! Temporary variable
      supsat6,                                                          &
                            ! Temporary variable
      supsat7,                                                          &
                            ! Temporary variable
      supsat8,                                                          &
                            ! Temporary variable
      r_row_length      ! Reciprocal of number of points on FI row)

!  Local dynamic arrays-------------------------------------------------
  INTEGER ::                                                            &
     icount(row_length,rows)       ! Counter of points

  LOGICAL ::                                                            &
     ocean(0:row_length+1,0:rows+1)! Those points which are not
!                     land, and have a sea-ice fraction less than 0.25.

  REAL ::                                                               &
     tl(0:row_length+1,0:rows+1,levels),                                &
                                            ! Conserved temperature
!                                                      (P292.1, UMDP29)
     qt(0:row_length+1,0:rows+1,levels),                                &
                                            ! Conserved WATER
!                                                      (P292.2, UMDP29)
     qst(0:row_length+1,0:rows+1),                                      &
                                        ! SATURATION VAPOUR PRESSURE
     p_grad(row_length,rows,levels),                                    &
                                        ! TERM WHICH RELATES
!                                         RH_SD_3 TO RH_SD_1
     supsat(0:row_length+1,0:rows+1),                                   &
                                        ! 'RELATIVE HUMIDITY' OF GRIDBOX
     al(0:row_length+1,0:rows+1),                                       &
                                        ! Defined in P292.6 in UMDP 29
     surf_mult(row_length,rows,8),                                      &
                                        ! Multiplier to take into
!                                         account surface matching.
     rhcpt_mean(levels)                                               
                                        ! Mean of first interior row.

  REAL                                                                  & 
     height,                                                            & 
         ! Height above the surface (m) 
     h_scale,                                                           & 
         ! Horizontal scale size of the grid box (m) 
     rh_crit_factor 
         ! Resolution-dependent Correction 

     ! Parameters for resolution-dependent factor calculation
  REAL, PARAMETER ::   cubefit_a = -8.37e-18
  REAL, PARAMETER ::   cubefit_b =  1.09e-11
  REAL, PARAMETER ::   cubefit_c = -4.46e-6
  REAL, PARAMETER ::   cubefit_d =  1.571354
  REAL, PARAMETER :: height_crit =  3380.0 
                        ! Designed to mimic BL_LEVEL=13 
                        ! in L50 set-up.

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!- End of Header

! ==Main Block==--------------------------------------------------------

  IF (lhook) CALL dr_hook('LS_CALC_RHCRIT',zhook_in,zhook_handle)

  ROOT_6=SQRT(6.)

  h_scale = Earth_Radius * delta_lambda 

  rh_crit_factor = (cubefit_a * h_scale * h_scale * h_scale)           &
                 + (cubefit_b * h_scale * h_scale)                     &
                 + (cubefit_c * h_scale) + cubefit_d 

! Levels_do1:
  DO k = 1, levels
! Rows_do1:
    DO j = 1, rows
! Rowlen_do1:
      DO i = 1, row_length
        p_grad(i,j,k) = p_layer_centres(i,j,k) - rhc_con4
        p_grad(i,j,k) = rh_crit_factor * (rhc_con1 +                                      &
       (rhc_con2 * p_grad(i,j,k) / (rhc_con3 + ABS(p_grad(i,j,k))) ) )

      END DO ! Rowlen_do1
    END DO ! Rows_do1
  END DO ! Levels_do1

! Ocean points defined now as not land and where ice fraction LT 0.25
! Rows_do2:
  DO j = 0, rows+1
! Rowlen_do2:
    DO i = 0, row_length+1
      ocean(i,j) = (land_frac(i,j) <  0.5) .AND.                        &
                                      (ice_frac(i,j) <  2.5e-1)
    END DO ! Rowlen_do2
  END DO ! Rows_do2

! A real no. is now assigned to every neighbouring point of every point
! on the grid, if their surfaces match it has the value one, else it is
! zero.
! Eight_do1:
  DO j8 = 1, 8
! Rows_do3:
    DO j = 1, rows
! Rowlen_do3:
      DO i = 1, row_length
        surf_mult(i,j,j8)=0.
      END DO ! Rowlen_do3
    END DO ! Rows_do3
  END DO ! Eight_do1

! Rows_do4:
  DO j = 1, rows
    jm1 = j - 1
    jp1 = j + 1
! Rowlen_do4:
    DO i = 1, row_length
      icount(i,j)=1
      im1 = i - 1
      ip1 = i + 1

      IF ( (ocean(i,j) .AND. ocean(im1,jm1))  .OR.                      &
           (.NOT.ocean(i,j) .AND. .NOT.ocean(im1,jm1)) ) THEN
        surf_mult(i,j,1) = 1.
        icount(i,j) = icount(i,j) + 1
      END IF

      IF ( (ocean(i,j) .AND. ocean(i,jm1))  .OR.                        &
           (.NOT.ocean(i,j) .AND. .NOT.ocean(i,jm1)) ) THEN
        surf_mult(i,j,2) = 1.
        icount(i,j) = icount(i,j) + 1
      END IF

      IF ( (ocean(i,j) .AND. ocean(ip1,jm1))  .OR.                      &
           (.NOT.ocean(i,j) .AND. .NOT.ocean(ip1,jm1)) ) THEN
        surf_mult(i,j,3) = 1.
        icount(i,j) = icount(i,j) + 1
      END IF

      IF ( (ocean(i,j) .AND. ocean(im1,j))  .OR.                        &
           (.NOT.ocean(i,j) .AND. .NOT.ocean(im1,j)) ) THEN
        surf_mult(i,j,4) = 1.
        icount(i,j) = icount(i,j) + 1
      END IF

      IF ( (ocean(i,j) .AND. ocean(ip1,j))  .OR.                        &
           (.NOT.ocean(i,j) .AND. .NOT.ocean(ip1,j)) ) THEN
        surf_mult(i,j,5) = 1.
        icount(i,j) = icount(i,j) + 1
      END IF

      IF ( (ocean(i,j) .AND. ocean(im1,jp1))  .OR.                      &
           (.NOT.ocean(i,j) .AND. .NOT.ocean(im1,jp1)) ) THEN
        surf_mult(i,j,6) = 1.
        icount(i,j) = icount(i,j) + 1
      END IF

      IF ( (ocean(i,j) .AND. ocean(i,jp1))  .OR.                        &
           (.NOT.ocean(i,j) .AND. .NOT.ocean(i,jp1)) ) THEN
        surf_mult(i,j,7) = 1.
        icount(i,j) = icount(i,j) + 1
      END IF

      IF ( (ocean(i,j) .AND. ocean(ip1,jp1))  .OR.                      &
           (.NOT.ocean(i,j) .AND. .NOT.ocean(ip1,jp1)) ) THEN
        surf_mult(i,j,8) = 1.
        icount(i,j) = icount(i,j) + 1
      END IF

    END DO ! Rowlen_do4
  END DO ! Rows_do4

! An initial sweep is done for all grid-cells, obtaining an initial
! estimate of the variance of the 3*3 grid.

! Levels_do5:
  DO k = 1, levels
! Rows_do5a:
    DO j = 0, rows+1
! Rowlen_do5a:
      DO i = 0, row_length+1
!   Calculate Tl and QT as in P292.1, P292.2 in UMDP 29.
!  (Assumes version 3A onwards of Section 4)
        tl(i,j,k) = t(i,j,k) - lsrcp*qcf(i,j,k)
        qt(i,j,k) = q(i,j,k) + qcf(i,j,k)

      END DO ! Rowlen_do5a
    END DO ! Rows_do5a

! DEPENDS ON: qsat_mix
    CALL qsat_mix(qst(0,0),tl(0,0,k),p_layer_centres(0,0,k),            &
              (row_length+2)*(rows+2),l_mixing_ratio)
! Rows_do5b:
    DO j = 0, rows+1
! Rowlen_do5b:
      DO i = 0, row_length+1
        IF (tl(i,j,k)  >   tm) THEN
          latht = lc/tl(i,j,k)
        ELSE
          latht = ls/tl(i,j,k)
        END IF
        al(i,j) = 1./ (1.+latht*latht*ercpr*qst(i,j))
! SUPSAT given by P292.3 of UMDP 29.
        supsat(i,j) = al(i,j)*(qt(i,j,k) - qst(i,j))
      END DO ! Rowlen_do5b
    END DO ! Rows_do5b

! Rows_do5c:
    DO j = 1, rows
      jm1 = j - 1
      jp1 = j + 1
! Rowlen_do5c:
      DO i = 1, row_length

! If the altitude is greater than height_crit (3km) then ensure all 9 
! neighbouring gridboxes are assumed to be similar. 
! This mimics a previous version of code where two different 
! calculations were done for k< and > BL_LEVELS. 
         height=r_theta_levels(i,j,k) - r_theta_levels(i,j,0)

         IF (height > height_crit) THEN 
           DO J8 = 1, 8 
             surf_mult(i,j,j8)=1.0 
             icount(i,j)=9.0 
           END DO 
         END IF 

        im1 = i - 1
        ip1 = i + 1
        supsat1 = surf_mult(i,j,1) * supsat(im1,jm1)
        supsat2 = surf_mult(i,j,2) * supsat(i,jm1)
        supsat3 = surf_mult(i,j,3) * supsat(ip1,jm1)
        supsat4 = surf_mult(i,j,4) * supsat(im1,j)
        supsat5 = surf_mult(i,j,5) * supsat(ip1,j)
        supsat6 = surf_mult(i,j,6) * supsat(im1,jp1)
        supsat7 = surf_mult(i,j,7) * supsat(i,jp1)
        supsat8 = surf_mult(i,j,8) * supsat(ip1,jp1)
        mean_supsat=(supsat1 + supsat2 + supsat3 + supsat4              &
                   + supsat5 + supsat6 + supsat7 + supsat8              &
                   + supsat(i,j)) / icount(i,j)
        tot_var = supsat1*supsat1 + supsat2*supsat2                     &
                + supsat3*supsat3 + supsat4*supsat4                     &
                + supsat5*supsat5 + supsat6*supsat6                     &
                + supsat7*supsat7 + supsat8*supsat8                     &
                + supsat(i,j)*supsat(i,j)                               &
                - icount(i,j)*mean_supsat*mean_supsat
        tot_var = ABS(tot_var)

!  Now remove the statistical outliers from the 3*3 region, so that
!  sigma, and hence RHcrit, is not biased by extreme values.
!  Points outside 2*sigma of the mean are considered outliers and are
!  rejected.
        IF (icount(i,j)  >   1) THEN
          two_sigma = 2. * SQRT(tot_var/icount(i,j))
        ELSE
          two_sigma = qst(i,j) * 0.01
        END IF
        COUNT=1

        IF (ABS(supsat(im1,jm1)-mean_supsat)  >   two_sigma) THEN
          supsat1 = 0.
        ELSE IF (surf_mult(i,j,1)  >   0.5) THEN
          COUNT = COUNT + 1
        END IF

        IF (ABS(supsat(i,jm1)-mean_supsat)  >   two_sigma) THEN
          supsat2 = 0.
        ELSE IF (surf_mult(i,j,2)  >   0.5) THEN
          COUNT = COUNT + 1
        END IF

        IF (ABS(supsat(ip1,jm1)-mean_supsat)  >   two_sigma) THEN
          supsat3 = 0.
        ELSE IF (surf_mult(i,j,3)  >   0.5) THEN
          COUNT = COUNT + 1
        END IF

        IF (ABS(supsat(im1,j)-mean_supsat)  >   two_sigma) THEN
          supsat4 = 0.
        ELSE IF (surf_mult(i,j,4)  >   0.5) THEN
          COUNT = COUNT + 1
        END IF

        IF (ABS(supsat(ip1,j)-mean_supsat)  >   two_sigma) THEN
          supsat5 = 0.
        ELSE IF (surf_mult(i,j,5)  >   0.5) THEN
          COUNT = COUNT + 1
        END IF

        IF (ABS(supsat(im1,jp1)-mean_supsat)  >   two_sigma) THEN
          supsat6 = 0.
        ELSE IF (surf_mult(i,j,6)  >   0.5) THEN
          COUNT = COUNT + 1
        END IF

        IF (ABS(supsat(i,jp1)-mean_supsat)  >   two_sigma) THEN
          supsat7 = 0.
        ELSE IF (surf_mult(i,j,7)  >   0.5) THEN
          COUNT = COUNT + 1
        END IF

        IF (ABS(supsat(ip1,jp1)-mean_supsat)  >   two_sigma) THEN
          supsat8 = 0.
        ELSE IF (surf_mult(i,j,8)  >   0.5) THEN
          COUNT = COUNT + 1
        END IF

        IF (COUNT >  1) THEN
          mean_supsat=(supsat1 + supsat2 + supsat3 + supsat4            &
                     + supsat5 + supsat6 + supsat7 + supsat8            &
                     + supsat(i,j)) / COUNT
          tot_var = supsat1*supsat1 + supsat2*supsat2                   &
                  + supsat3*supsat3 + supsat4*supsat4                   &
                  + supsat5*supsat5 + supsat6*supsat6                   &
                  + supsat7*supsat7 + supsat8*supsat8                   &
                  + supsat(i,j)*supsat(i,j)                             &
                  - COUNT*mean_supsat*mean_supsat
          tot_var = ABS(tot_var)
          supsat_sd_3 = SQRT(tot_var/COUNT)
        ELSE
!           Limit the 3*3 grid variance when scatter is large.
          supsat_sd_3 = qst(i,j)*0.01
        END IF

! Try to detect if the central point (i,j) is an outlier, as can happen
! in a GPS situation. If so, set the variance to a small value.
        IF (COUNT >  2) THEN
          mean_supsat=(supsat1 + supsat2 + supsat3 + supsat4            &
           + supsat5 + supsat6 + supsat7 + supsat8) / (COUNT-1.0)
          tot_var = supsat1*supsat1 + supsat2*supsat2                   &
                  + supsat3*supsat3 + supsat4*supsat4                   &
                  + supsat5*supsat5 + supsat6*supsat6                   &
                  + supsat7*supsat7 + supsat8*supsat8                   &
                  - COUNT*mean_supsat*mean_supsat
          IF (ABS(supsat(i,j)-mean_supsat) >                            &
                      (2.0*SQRT(ABS(tot_var)/COUNT))) THEN
            supsat_sd_3 = qst(i,j)*0.01
          END IF
        END IF

! P_GRAD determines the relation between 3*3 and sub-grid variance.
        supsat_sd_1 = p_grad(i,j,k) * supsat_sd_3
! RHCPT defined from P292.14 in UMDP 29
        rhcpt(i,j,k) = 1. - (root_6*supsat_sd_1 /(al(i,j)*qst(i,j)))
! RHcrit is now limited to lie between a range defined in RHCCON2B
        rhcpt(i,j,k) = MAX(rhcpt(i,j,k),rhc_min)
        rhcpt(i,j,k) = MIN(rhcpt(i,j,k),rhc_max)
      END DO ! Rowlen_do5c
    END DO ! Rows_do5c

  END DO ! Levels_do5

! Tidy up at South Pole : Pole is mean of first interior row.
! SouthPole_if1:
  IF (at_extremity(PSouth)) THEN

    ij_field = row_length * rows
    i_length = row_length

! Row number of first interior (ie. non-polar) row
    j_start  = 2

    r_row_length = 1. / global_row_length
!       Number of points in sum should be global_row_length: might
!       need to modify r_row_length if i_length lt row_length.
!       r_row_length = r_row_length * row_length / i_length

! Sum over points in PEs in order along first interior row
! (gc_proc_row_group is group ID for rows of PEs, here only PSouth).
    CALL global_2d_sums(rhcpt(:,j_start:j_start,:), i_length, 1, 0, 0,  &
                        levels, rhcpt_mean, gc_proc_row_group)

! Levels_do7:
    DO k = 1, levels
      rhcpt_mean(k) = rhcpt_mean(k) * r_row_length
! Rowlen_do7:
      DO i = 1, row_length
        rhcpt(i,1,k) = rhcpt_mean(k)
      END DO ! Rowlen_do7
    END DO  ! Levels_do7

  END IF  ! SouthPole_if1

! Tidy up at North Pole : Pole is mean of first interior row.
! NorthPole_if1:
  IF (at_extremity(PNorth)) THEN

    ij_field = row_length * rows
    i_length = row_length

! Row number of first interior (ie. non-polar) row
    j_start  = rows - 1

    r_row_length = 1. / global_row_length
!       Number of points in sum should be global_row_length: might
!       need to modify r_row_length if i_length lt row_length.
!       r_row_length = r_row_length * row_length / i_length

! Sum over points in PEs in order along first interior row
! (gc_proc_row_group is group ID for rows of PEs, here only PNorth).
    CALL global_2d_sums(rhcpt(:,j_start:j_start,:), i_length, 1, 0, 0,  &
                        levels, rhcpt_mean, gc_proc_row_group)

! Levels_do8:
    DO k = 1, levels
      rhcpt_mean(k) = rhcpt_mean(k) * r_row_length
! Rowlen_do8:
      DO i = 1, row_length
        rhcpt(i,rows,k) = rhcpt_mean(k)
      END DO ! Rowlen_do8
    END DO  ! Levels_do8

  END IF  ! NorthPole_if1

  IF (lhook) CALL dr_hook('LS_CALC_RHCRIT',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE ls_calc_rhcrit
