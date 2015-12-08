! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_calc_eta_levels_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_calc_eta_levels(                                    &
                      model_levels,me,first_constant_rho_level,       &
                      eta_theta_levels, eta_rho_levels,               &
!  Grid information
                      grid_number, height_domain,                     &
                      first_theta_height, thin_theta_height,          &
                      big_layers, transit_layers, big_factor,         &
                      vert_grid_ratio)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE integrity_mod

IMPLICIT NONE
!
! Description:
!  
!          Sets up vertical grid (eta_levels) 
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



INTEGER                                                                &
  me
INTEGER                                                                &
  grid_number,                                                         &
  big_layers,                                                          &
  transit_layers

REAL                                                                   &
  height_domain,                                                       &
  first_theta_height,                                                  &
  thin_theta_height,                                                   &
  big_factor,                                                          &
  vert_grid_ratio

INTEGER                                                                &
  model_levels,                                                        &
                   ! number of model levels
  first_constant_rho_level

REAL                                                                   &
     ! vertical co-ordinate information
  eta_theta_levels(0:model_levels),                                    &
  eta_rho_levels(model_levels)

! local variables

REAL                                                                   &
  delta_z,                                                             &
  delta_eta,                                                           &
  mag_value1,                                                          &
  mag_value2,                                                          &
  mag_value3,                                                          &
  mag_factor1,                                                         &
  mag_factor2,                                                         &
  mag_factor3,                                                         &
  top_mag_value,                                                       &
  mag,                                                                 &
  mag3

INTEGER                                                                &
 k,                                                                    &
 fcrl,                                                                 &
 level,                                                                &
 mod_layers

 PARAMETER( mag_value1 = 1.04 )    ! value for 38 levels
 PARAMETER( mag_value2 = 1.1001 )  ! value for 38 levels
 PARAMETER( mag_value3 = 1.1)      ! value for 38 levels
 PARAMETER( top_mag_value = 1.5 )  ! value for 38 levels
 PARAMETER( mag_factor1 = 0.01 )   ! value for 38 levels
 PARAMETER( mag_factor2 = 0.02 )   ! value for 38 levels
 PARAMETER( mag_factor3 = 0.02)    ! value for 38 levels

! Description: COMDECK containing vertical grid types
!  for use in idealised problems
!
      INTEGER, PARAMETER :: vert_regular=1
      INTEGER, PARAMETER :: vert_quadratic_theta=21
      INTEGER, PARAMETER :: vert_bi_quadratic=22
      INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
      INTEGER, PARAMETER :: vert_schar=3
      INTEGER, PARAMETER :: vert_dwd=4
      INTEGER, PARAMETER :: vert_stretch_plus_regular=5
      INTEGER, PARAMETER :: vert_quad_stretch_thin=6
      INTEGER, PARAMETER :: vert_regular_thin=7
      INTEGER, PARAMETER :: vert_geometric_theta=8
      INTEGER, PARAMETER :: vert_dump=10
      INTEGER, PARAMETER :: vert_idl_um_grid=11

! Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_IDL_CALC_ETA_LEVELS',zhook_in,zhook_handle)


! ---------------------------------------------------------------------
! Section 1.  Initialise Data fields.
!             Set reference vertical grid
!             Set eta_theta_levels
! ---------------------------------------------------------------------

!   Change first_constant_rho_level to ensure levels are flat in
!   any stretching region
IF ( grid_number  ==  vert_stretch_plus_regular .and.                 &
     big_layers  >   0) THEN
  fcrl = model_levels - big_layers - transit_layers
  IF(fcrl  <   first_constant_rho_level)THEN
    first_constant_rho_level = fcrl
  END IF  !fcrl  <   first_constant_rho_level
END IF   !big_layers  >   0
!   fcrl shorthand for first_constant_rho_level
fcrl = first_constant_rho_level

eta_theta_levels(0) = 0.0

IF(grid_number  ==  vert_regular )THEN

!!!!!!!!!!      Regular grid start   !!!!!!!!!!!!!!!!!!!!!!!!!

  delta_eta = 1.0/REAL(model_levels)
  eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to regular grid
  DO k=1, model_levels
    eta_theta_levels(k) = REAL(k) * delta_eta
  END DO
  IF(me  ==  0)THEN
    WRITE(UNIT=6,FMT='(A)') ' Regular vertical grid selected'
    WRITE(UNIT=6,FMT='(A,E16.8)') '   eta_theta_levels(top) = ',          &
                         eta_theta_levels(model_levels)
    WRITE(UNIT=6,FMT='(A,E16.8,A)') '   height_domain =',height_domain,   &
                        ' metres'
    WRITE(UNIT=6,FMT='(A,E16.8)') '   each layer thickness =',            &
                         height_domain*delta_eta
  END IF    !(me  ==  0)
! Set eta_rho_levels
  DO k=1,model_levels
    eta_rho_levels(k) =                                               &
    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
  END DO
!!!!!!!!!!      Regular grid end  !!!!!!!!!!!!!!!!!!!!!!!!!

ELSE IF(grid_number  ==  vert_quadratic_theta )THEN

!!!!!!!!!!     Quadratic grid start  !!!!!!!!!!!!!!!!!!!!!!!!!

  delta_eta = 1.0/REAL(model_levels)
  eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to quadratic grid
  DO k=1, model_levels
    eta_theta_levels(k) = (REAL(k) * delta_eta)                       &
                        * (REAL(k) * delta_eta)
  END DO
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A)')'***** Quadratic grid for theta levels *****'
    WRITE(6,fmt='(A,E16.8)')'eta_theta_levels(top) =',                    &
    eta_theta_levels(model_levels)
    WRITE(6,fmt='(A,E16.8,A)')'height_domain =',height_domain,' metres'
  END IF      !(me  ==  0)
! Set eta_rho_levels
  DO k=1,model_levels
    eta_rho_levels(k) =                                               &
    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
  END DO
!!!!!!!!!!      Quadratic grid for theta levels end  !!!!!!!!

ELSE IF(grid_number  ==  vert_bi_quadratic)THEN

!!!!!!!!!!     Quadratic u & theta grids start  !!!!!!!!!!!!!!!!!    
  delta_eta = 1.0/REAL(model_levels)
  eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to quadratic grid
  DO k=1, model_levels
    eta_theta_levels(k) = (REAL(k) * delta_eta)                       &
                        * (REAL(k) * delta_eta)
  END DO
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A)')'***** Quadratic grid for u & theta levels *****'
    WRITE(6,fmt='(A,E16.8)')'eta_theta_levels(top) ='                     &
                        ,eta_theta_levels(model_levels)    
    WRITE(6,fmt='(A,E16.8,A)')'height_domain =',height_domain,' metres'
  END IF    !(me  ==  0)
!!!!!!!!!!      Quadratic grids for u & theta levels end  !!!!!!!!

ELSE IF(grid_number  ==  vert_quadratic_uvtheta)THEN

!!!!!!!!!!    Common Quadratic u & theta grid start  !!!!!!!!!!

  delta_eta = 1.0/REAL(2*model_levels)
  eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to quadratic grid
  DO k=1, model_levels
    eta_rho_levels(k) = (REAL(2*k-1) * delta_eta)                     &
                      * (REAL(2*k-1) * delta_eta)
    eta_theta_levels(k) = (REAL(2*k) * delta_eta)                     &
                      * (REAL(2*k) * delta_eta)    
  END DO
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A)')'** Common Quadratic grid for u & theta levels **'
    WRITE(6,fmt='(A,E16.8)')'eta_theta_levels(top) =',                    &
                         eta_theta_levels(model_levels)
    WRITE(6,fmt='(A,E16.8,A)')'height_domain =',height_domain,' metres'
  END IF    !(me  ==  0)

!!!!!!!!  Common Quadratic grid for u & theta levels end  !!!!!!!!

ELSE IF(grid_number  ==  vert_schar)THEN

!!!!!!!!!!     Schar grid start  !!!!!!!!!!!!!!!!!!!!!!!!!

  DO k=0, model_levels
! These values for eta_theta_levels equate approximately to those used
! by Schar.
    eta_theta_levels(k) = REAL(k)**1.2/                               &
                        REAL(model_levels)**1.2
  END DO
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A)')'***** Schar grid selected *****'
    WRITE(6,fmt='(A,E16.8)')'eta_theta_levels(top) =',                    &
                         eta_theta_levels(model_levels)
    WRITE(6,fmt='(A,E16.8,A)')'height_domain =',height_domain,' metres'
  END IF    !(me  ==  0)

! Set eta_rho_levels
  DO k=1,model_levels
    eta_rho_levels(k) =                                               &
    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
  END DO
!!!!!!!!!!     Schar grid end  !!!!!!!!!!!!!!!!!!!!!!!!!

ELSE IF(grid_number  ==  vert_dwd )THEN

!!!!!!!!!!      DWD stretched start   !!!!!!!!!!!!!!!!!!!!!!!!!

! These values for eta_theta_levels similar to DWD
! Set eta_theta_levels to actual height values
  eta_theta_levels(0) = 0.0
  eta_theta_levels(1) = 40.0
  eta_theta_levels(2) = 100.0
  delta_z = 100.0
  DO k = 3, model_levels
! These values for eta_theta_levels equate to stretched  grid
    eta_theta_levels(k) = eta_theta_levels(k-1) + delta_z
    delta_z = delta_z + 40.0
  END DO
  height_domain = eta_theta_levels( model_levels)
! Now normalise eta_theta_levels
  DO k=1,model_levels
    eta_theta_levels(k) =                                             &
    eta_theta_levels(k)/eta_theta_levels( model_levels)
  END DO
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A)')'***** DWD stretched grid selected *****'
    WRITE(6,fmt='(A,E16.8,A)')'eta_theta_levels(top) =',                  &
                           eta_theta_levels(model_levels)
    WRITE(6,fmt='(A,E16.8,A)')'height_domain =',height_domain,' metres'
  END IF    !(me  ==  0)
! Set eta_rho_levels
  DO k=1,model_levels
    eta_rho_levels(k) =                                               &
    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
  END DO
!!!!!!!!!!      DWD stretched end  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ELSE IF(grid_number  ==  vert_stretch_plus_regular)THEN

!!!!!!!!!!     regular grid followed by stretched (big_layers)!!!!!

! mod_layers is the number of regular layers
  mod_layers=model_levels-transit_layers-big_layers
! mag is the magnification factor applied over transit_layers
  mag=big_factor**(1.0/REAL(MAX(1,transit_layers)))
! interval chosen to scale to 1 over domain
!  NB   mag * (1.0 - big_factor)/(1.0 - mag) +
!    if first transit layer magnified
  delta_eta = 1.0/ (   REAL(mod_layers) +                             &
             (1.0 - big_factor)/(1.0 - mag) +                         &
             big_factor * REAL(big_layers) )
! Make height domain for problem the full height domain
! For this grid the INPUT height domain is for regular levels
  height_domain = height_domain/(REAL(mod_layers)*delta_eta)
  eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to regular grid
  DO k=1, mod_layers
    eta_theta_levels(k) = REAL(k) * delta_eta
  END DO
  DO k = mod_layers + 1 ,mod_layers + transit_layers
! These values for eta_theta_levels equate to stretched  grid
    eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
    delta_eta = delta_eta * mag
  END DO
! These values for eta_theta_levels equate to stretched  grid
  DO k = mod_layers + transit_layers + 1 , model_levels
    eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
  END DO
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A)')                                                &
                '** regular grid followed by stretched grid selected'
    WRITE(6,fmt='(A,E16.8,A,A,I5,A)')'     regular grid =',                &
                               height_domain*delta_eta,' metres'      &
                               ,' over ',mod_layers,' layers'
    WRITE(6,fmt='(A,I5,A,A,I5,A)')'     stretching over ',              &
                               transit_layers,' layers',              &
                         ' with magnification factor =',mag,' times'
    WRITE(6,fmt='(A,I5,A,I5,A)')'     thick ',big_factor,               &
                         '* regular over ',                           &
                         big_layers,' layers'
    WRITE(6,fmt='(A,E16.8)')'eta_theta_levels(top) =',                    &
    eta_theta_levels(model_levels)
    WRITE(6,fmt='(A,E16.8,A)')'height_domain =',height_domain,' metres'
  END IF    !(me  ==  0)
! Set eta_rho_levels
  DO k = 1, model_levels
    eta_rho_levels(k) =                                               &
    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
  END DO

!!!!!!! end of regular grid followed by stretched (big_layers)!!!

ELSE IF(grid_number  ==  vert_quad_stretch_thin)THEN
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A)')'********grid option 6  ****'
    WRITE(6,fmt='(A)')'*quadratic grid from level 2 *'
    WRITE(6,fmt='(A)')'*        followed by expanding layers    *  '
    WRITE(6,fmt='(A)')'*  thin layer near surface    *  '
  END IF    !(me  ==  0)

! For this grid need to work with model_levels - 1 to generate
!  quadratic part
  delta_eta = 1.0/REAL(model_levels - 1)
  eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to quadratic grid
  DO k = 1, model_levels - 1
    eta_theta_levels(k) = (REAL(k) * delta_eta)                       &
                      * (REAL(k) * delta_eta)
  END DO
! Use eta_rho_levels to hold ratio of succesive eta differences
  level = 0
  DO k = 1, model_levels - 2
   eta_rho_levels(k)=(eta_theta_levels(k+1)- eta_theta_levels(k))     &
                    /(eta_theta_levels(k)- eta_theta_levels(k-1))
    IF (level == 0) THEN
      IF (eta_rho_levels(k) <   mag_value1) level = k
    END IF
  END DO

  IF (level  /=  0) THEN
    delta_eta = eta_theta_levels(level ) -                            &
                eta_theta_levels(level - 1)
    mag = mag_value1
    mag3 = mag_factor3
    DO k = level + 1 , model_levels - 1
      delta_eta = mag * delta_eta
      eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
      mag = mag + mag_factor1
      IF(mag >  mag_value3)THEN
        mag = mag - mag_factor1 + mag3
        mag3 = 2.0 * mag3
        IF(mag >  2.0)mag = 1.3
      ELSE IF(mag >  mag_value2)THEN
        mag = mag - mag_factor1 + mag_factor2
      END IF
    END DO
  END IF  !  level  /=  0

! Re-normalise eta_theta_levels and set eta_rho_levels
! Store current top of quadratic grid in eta_theta_levels(model_levels)
  eta_theta_levels(model_levels) =                                    &
                   eta_theta_levels(model_levels - 1)    
! Work from top to allow level 1 values to be inserted
  DO k= model_levels - 1, 2, -1
    eta_theta_levels(k) = eta_theta_levels(k-1) /                     &
                        eta_theta_levels(model_levels)
  END DO
!!!!!  reset height domain to give r_theta = first_theta_height
!!!!!                                    at level 1
  height_domain = first_theta_height / eta_theta_levels(2)
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A,E16.8,A)')'height_domain =',height_domain,' metres'
  END IF    !(me  ==  0)
! Re-normalise eta_theta_levels(model_levels)
  eta_theta_levels(model_levels) = 1.0
! Now add in thin bottom layer
  eta_theta_levels(1) = thin_theta_height / height_domain
  eta_rho_levels(1) = 0.5 * eta_theta_levels(1)
  DO k= 2, model_levels
    eta_rho_levels(k) =                                               &
        0.5 * (eta_theta_levels(k-1) + eta_theta_levels(k) )
  END DO

!!!!!! end of quadratic grid followed by expanding layers!!!

ELSE IF(grid_number  ==  vert_geometric_theta )THEN

!!!!!!!!!! Geometric spacing of theta levels   !!!!!!!!!!!!!!!!!!!!!!!!!

  delta_eta = 1.0/REAL(model_levels)
  eta_theta_levels(0) = 0.0

  DO k=1, model_levels
    eta_theta_levels(k)=eta_theta_levels(k-1)+delta_eta
    eta_rho_levels(k)=0.5*(eta_theta_levels(k-1)+eta_theta_levels(k))
    delta_eta=delta_eta*vert_grid_ratio
  END DO
  eta_rho_levels(:)=eta_rho_levels(:)/eta_theta_levels(model_levels)
  eta_theta_levels(:)=eta_theta_levels(:)/eta_theta_levels(model_levels)

  IF(me  ==  0)THEN
    WRITE(UNIT=6,FMT='(A)') ' Semi-Geometric vertical grid selected'
    WRITE(UNIT=6,FMT='(A)') '   eta_theta_levels = '
    DO k=1, model_levels
     WRITE(UNIT=6,FMT='(E16.8)')      eta_theta_levels(k)
    END DO
    WRITE(UNIT=6,FMT='(A)') '   eta_rho_levels   = '
    DO k=1, model_levels
      WRITE(UNIT=6,FMT='(E16.8)')  eta_rho_levels(k)
    END DO  
  END IF    !(me  ==  0)
!!!!!!!!!! Geometric spacing of theta levels end !!!!!!!!!!!!!!!!!!!!!!!


ELSE IF(grid_number  >   vert_dump )THEN
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A)')'********grid option not supported****'
  END IF    !(me  ==  0)

END IF  ! grid_number options

IF (integrity_test)                                                   &
  CALL update_hash_m(                                                 &
                  eta_theta_levels,SIZE(eta_theta_levels),'etatl',    &
                  eta_rho_levels,  SIZE(eta_rho_levels),  'etarl')

IF (lhook) CALL dr_hook('EG_IDL_CALC_ETA_LEVELS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_calc_eta_levels
END MODULE eg_idl_calc_eta_levels_mod
