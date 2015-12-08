! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
SUBROUTINE eg_idl_surface_setup(                                      &
                      row_length, rows, model_levels                  &
,                     halo_i, halo_j                                  &
,                     me, n_proc, at_extremity, model_domain          &
,                     n_rows                                          &
,                     orography, orog_haloes, orog_u, orog_v          &
,                     surface_type                                    &
,                     h_o, lambda_fraction, phi_fraction              &
,                     half_width_x, half_width_y)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE earth_constants_mod, ONLY: earth_radius, g
USE conversions_mod, ONLY : pi
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE ereport_mod, ONLY : ereport
USE PrintStatus_mod
USE UM_ParParams
IMPLICIT NONE
!
! Description:
!          Sets up surface
!          surface_number defines options

!  
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


!*** Input arguments/declarations need tidying up, many are unused now!

INTEGER ::                                                            &
  row_length                                                          &
                   ! number of points on a processor row
, rows                                                                &
                   ! number of rows in a processor theta field
, model_levels                                                        &
                   ! number of model levels
, halo_i                                                              &
                       ! Size of halo in i direction.
, halo_j                                                              &
                       ! Size of halo in j direction.
, n_rows                                                              
                   ! number of rows in a processor v field


! Orography arrays
!    orog_haloes is r_theta_levels(0)
REAL :: orog_haloes(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j),&
        orog_u(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j),     &
        orog_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

!    orography has no haloes - for D1 array
REAL :: orography(row_length,rows) 

! Mountain dimensions
REAL :: h_o, lambda_fraction, phi_fraction, half_width_x, half_width_y     

INTEGER :: model_domain, surface_type

! this processor and Total number of processors
INTEGER :: me, n_proc                                          

! Indicates if this processor is at north,
! south, east or west of the processor grid
LOGICAL :: at_extremity(4)                           

! local variables
! loop counters and local halos
INTEGER :: i, j, haloi, haloj, info                                                                
! Stop-gap for introducing ENDGame indexing convention. Should be
! removed when whole subroutine is recoded under ENDGame convention
INTEGER :: i_shift, j_shift
! REMOVE ABOVE LINES (when ENDGame indexing is applied)

REAL ::x, y, h_max

! Description: COMDECK containing surface types
!  for use in idealised problems
!

      INTEGER, PARAMETER :: surface_zero=0
      INTEGER, PARAMETER :: surface_ellipse=1
      INTEGER, PARAMETER :: surface_ridge=2
      INTEGER, PARAMETER :: surface_plateau=3
      INTEGER, PARAMETER :: surface_massif=4
      INTEGER, PARAMETER :: surface_mask=5
      INTEGER, PARAMETER :: surface_gauss=6
      INTEGER, PARAMETER :: surface_ridge_series=7
! ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_schar_ridge=8
      INTEGER, PARAMETER :: surface_baroclinic=9
! End of ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_dump=10

! External routines
REAL, EXTERNAL :: eg_baro_geo_p

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_IDL_SURFACE_SETUP',zhook_in,zhook_handle)

IF (me  ==  0 .AND. PrintStatus >= PrStatus_Normal) THEN
  WRITE(6,fmt='(A,I2)')' Model_domain = ',model_domain
  IF ( model_domain  ==  mt_global ) THEN
    WRITE(6,fmt='(A,I2)')' Global model_domain  =  ',model_domain
  ELSE IF ( model_domain  ==  mt_lam ) THEN
    WRITE(6,fmt='(A,I2)')' Limited-area model_domain  =  ',model_domain
  ELSE IF ( model_domain  ==  mt_cyclic_lam ) THEN
    WRITE(6,fmt='(A,I2)')' Cyclic (East-West) Limited-area model_domain  =  '  &
          ,model_domain
  ELSE IF ( model_domain  ==  mt_bi_cyclic_lam ) THEN
    WRITE(6,fmt='(A,I2)')' Bi-Cyclic Limited-area model_domain  =  '           &
          ,model_domain
  END IF ! model_domain  ==  mt_global
END IF ! (me == 0)

  IF ( model_domain  /=  mt_global         .AND.                               &
       model_domain  /=  mt_lam            .AND.                               &
       model_domain  /=  mt_cyclic_lam     .AND.                               &
       model_domain  /=  mt_bi_cyclic_lam                                      &
     ) THEN
       CALL ereport( 'EG_IDL_SURFACE_SETUP',model_domain,                      &
                     'model domain NOT SUPPORTED'  )
  END IF ! model_domain  ==  mt_global

! ---------------------------------------------------------------------
! Section 2.  Insert mountain and set up orography arrays
! ---------------------------------------------------------------------


!-----------------------------------------------------------------------
! flat surface
!-----------------------------------------------------------------------
IF(surface_type  ==  surface_zero )THEN

  IF(me  ==  0 .AND. PrintStatus >= PrStatus_Normal)THEN
    WRITE(6,fmt='(A,I2)')' Flat surface_type = ',surface_type
    WRITE(6,fmt='(A)')' Orography set to 0.0 everywhere '
  END IF ! (me == 0)

!  Make sure h_o = 0 for diagnostic printing
  h_o = 0.0

  orog_haloes = 0.0
  orog_u = 0.0
  orog_v = 0.0

!-----------------------------------------------------------------------
! Barolcinic wave surface
!-----------------------------------------------------------------------
ELSE IF(surface_type  ==  surface_baroclinic) THEN

! Surface orography to give lower boundary surface pressure of
! 1000hPa in baroclinic wave test (QJRMS 132, 2943--2975).

  IF(me  ==  0 .AND. PrintStatus >= PrStatus_Normal)THEN
    WRITE(6,fmt='(A,I2)')' surface_type = ',surface_type
    WRITE(6,fmt='(A)')' Orography for Jablonowski & Williamson baroclinic'
    WRITE(6,fmt='(A)')' wave test (QJRMS 132, 2943--2975).'
    WRITE(6,fmt='(A)')
    WRITE(6,fmt='(A,2I4)') 'pdims%i_start, pdims%j_start: ', pdims%i_start, pdims%j_start 
    WRITE(6,fmt='(A,2I4)') 'udims%i_start, vdims%j_start: ', udims%i_start, vdims%j_start
    WRITE(6,fmt='(A,2I4)') 'haloi, haloj: ', haloi, haloj
  END IF ! (me == 0)

! theta-point orog
  DO j=pdims%j_start,pdims%j_end
    DO i=pdims%i_start,pdims%i_end
! DEPENDS ON: eg_baro_geo_p
      orog_haloes(i,j) = eg_baro_geo_p(xi1_p(i), xi2_p(j), 1.0)/g
    END DO
  END DO

! u-point orog
  DO j=pdims%j_start,pdims%j_end
    DO i=udims%i_start,udims%i_end
      i_shift=i+1 ! Stop-gap to convert to ENDGame indexing convention
!       i_shift=i
! DEPENDS ON: eg_baro_geo_p
      orog_u(i_shift,j) = eg_baro_geo_p(xi1_u(i), xi2_p(j), 1.0)/g
    END DO
  END DO

! v-point orog
  DO j=vdims%j_start,vdims%j_end
    j_shift=j+1 ! Stop-gap to convert to ENDGame indexing convention
!     j_shift=j
    DO i=pdims%i_start,pdims%i_end
! DEPENDS ON: eg_baro_geo_p
      orog_v(i,j_shift) = eg_baro_geo_p(xi1_p(i), xi2_v(j), 1.0)/g
    END DO
  END DO

!-----------------------------------------------------------------------
! Schar ridge surface
!-----------------------------------------------------------------------

ELSE IF(surface_type  ==  surface_schar_ridge) THEN

  IF(me  ==  0 .AND. PrintStatus >= PrStatus_Normal)THEN
    WRITE(6,fmt='(A,I2)')' surface_type = ',surface_type
    WRITE(6,fmt='(A)')' North-South ridge'
    WRITE(6,fmt='(A,E16.8,A,E16.8,A)')' ridge height = ',h_o,         &
                              ' metres. Centre at '                   &
          ,360.0*lambda_fraction,' degrees'
    WRITE(6,fmt='(A,E16.8,A)')' ridge half-width = ',half_width_x,' radians'
    WRITE(6,*)
    WRITE(6,fmt='(A,2I4)') 'pdims%i_start, pdims%j_start: ',          &
                            pdims%i_start, pdims%j_start 
    WRITE(6,fmt='(A,2I4)') 'udims%i_start, vdims%j_start: ',          &
                            udims%i_start, vdims%j_start
    WRITE(6,fmt='(A,2I4)') 'haloi, haloj: ', haloi, haloj
  END IF ! (me == 0)

! theta-point orog
  DO j=pdims%j_start,pdims%j_end
    DO i=pdims%i_start,pdims%i_end      
      x = xi1_p(i) - 2.0*pi*lambda_fraction
      orog_haloes(i,j) = h_o * COS(pi*x/(half_width_x))**2            &
                           *EXP(-(x/half_width_y)**2) 
    END DO
  END DO

! u-point orog
  DO j=pdims%j_start,pdims%j_end
    DO i=udims%i_start,udims%i_end
      i_shift=i+1 ! Stop-gap to convert to ENDGame indexing convention
      x = xi1_u(i) - 2.0*pi*lambda_fraction
      orog_u(i_shift,j) = h_o * COS(pi*x/(half_width_x))**2           &
                           *EXP(-(x/half_width_y)**2)
    END DO
  END DO

! v-point orog
  DO j=vdims%j_start,vdims%j_end
    j_shift=j+1 ! Stop-gap to convert to ENDGame indexing convention
    DO i=pdims%i_start,pdims%i_end
      x = xi1_p(i) - 2.0*pi*lambda_fraction
      orog_v(i,j_shift) = h_o * COS(pi*x/(half_width_x))**2           &
                           *EXP(-(x/half_width_y)**2)
    END DO
  END DO

!-----------------------------------------------------------------------
! Gaussian hill surface
!-----------------------------------------------------------------------

ELSE IF(surface_type  ==  surface_gauss)THEN

  IF(me  ==  0.AND. PrintStatus >= PrStatus_Normal)THEN
    WRITE(6,fmt='(A,I2,A)')' surface_type = ',surface_type,           &
                          ' ** Gaussian **'
    WRITE(6,fmt='(A,E16.8,A,E16.8,A,E16.8,A)')' hill height = ',h_o,  &
                                  ' metres. East-West centre at ',    &
                                  360.0*lambda_fraction,              &
                                  ' degrees, north-south centre at ', &
                                  180.0*phi_fraction-90.0,            &
                                  ' degrees'
    IF(half_width_x  ==  half_width_y)THEN
      WRITE(6,fmt='(A,E16.8,A)')' circular hill half-width = ',       &
                            half_width_x,' radians'
    ELSE
      WRITE(6,fmt='(A,E16.8,A,A,E16.8,A)')'elliptical hill half_width_x = ', &
                                   half_width_x,' radians',           &
                             ' half_width_y = ',half_width_y,' radians'
    END IF     !(half_width_x  ==  half_width_y)
    WRITE(6,*)
    WRITE(6,fmt='(A,2I4)') 'pdims%i_start, pdims%j_start: ',          &
                            pdims%i_start, pdims%j_start 
    WRITE(6,fmt='(A,2I4)') 'udims%i_start, vdims%j_start: ',          &
                            udims%i_start, vdims%j_start
    WRITE(6,fmt='(A,2I4)') 'haloi, haloj: ', haloi, haloj
  END IF ! (me == 0)

! theta-point orog
  DO j=pdims%j_start,pdims%j_end
    y = xi2_p(j) - pi*(phi_fraction-0.5)
    DO i=pdims%i_start,pdims%i_end
      x = xi1_p(i) - 2.0*pi*lambda_fraction
      orog_haloes(i,j) = h_o * EXP( -1.0 * (                          &
                         ( x/half_width_x)**2. +                      &
                         ( y/half_width_y)**2.))    
    END DO
  END DO
  IF (PrintStatus >= PrStatus_Normal)THEN
    WRITE(6,fmt='(A,I5)') 'Me = ',me
    WRITE(6,fmt='(2I4,E16.8)') pdims%j_start,pdims%j_end,phi_fraction
    WRITE(6,fmt='(2E16.8)') xi2_p(pdims%j_start),xi2_p(pdims%j_end)
    WRITE(6,*)' '
  END IF
! u-point orog
  DO j=pdims%j_start,pdims%j_end
    y = xi2_p(j) - pi*(phi_fraction-0.5)
    DO i=udims%i_start,udims%i_end
      i_shift=i+1 ! Stop-gap to convert to ENDGame indexing convention
      x = xi1_u(i) - 2.0*pi*lambda_fraction
      orog_u(i_shift,j) = h_o * EXP( -1.0 * (                         &
                    ( x/half_width_x)**2. +                           &
                    ( y/half_width_y)**2.)) 
    END DO
  END DO

! v-point orog
  DO j=vdims%j_start,vdims%j_end
    j_shift=j+1 ! Stop-gap to convert to ENDGame indexing convention
    y = xi2_v(j) - pi*(phi_fraction-0.5)
    DO i=pdims%i_start,pdims%i_end
      x = xi1_p(i) - 2.0*pi*lambda_fraction
      orog_v(i,j_shift) = h_o * EXP( -1.0 * (                         &
                    ( x/half_width_x)**2. +                           &
                    ( y/half_width_y)**2.)) 
    END DO
  END DO


!-----------------------------------------------------------------------
! Surface from start dump
!-----------------------------------------------------------------------
ELSE IF(surface_type  ==  surface_dump)THEN

!  Zero r_theta_levels(i,j,0)
  DO j = 1-halo_j, rows+halo_j
    DO i = 1-halo_i, row_length+halo_i
      orog_haloes(i,j)  = 0.0
    END DO
  END DO
!  Copy orography into r_theta_levels(i,j,0)
  DO j = 1, rows
    DO i = 1, row_length
      orog_haloes(i,j)  = orography(i,j)
    END DO
  END DO

! DEPENDS ON: swap_bounds
  CALL swap_bounds(orog_haloes, row_length, rows, 1,                  &
                 halo_i, halo_j, fld_type_p,.FALSE.)

 DO j = 1, n_rows
    DO i = 1, row_length
!     NOTE: these arrays are _not_ declared in ENDGame notation,
!           hence the -1
      orog_v(i,j) = 0.5*(orog_haloes(i,j-1)+orog_haloes(i,j))
    END DO
  END DO

  DO j = 1, rows
    DO i = 1, row_length
!     NOTE: these arrays are _not_ declared in ENDGame notation,
!           hence the -1
      orog_u(i,j) = 0.5*(orog_haloes(i-1,j)+orog_haloes(i,j))
    END DO
  END DO

!   Reset h_o and use to store max orography
  h_max = 0.0

  DO j=1,rows
    DO i=1,row_length
!  Set h_o to max orography found
      IF (orography(i,j)  >   h_max)THEN
        h_max = orography(i,j)
      END IF
    END DO
  END DO
! find max  over all processors
  CALL gc_rmax(1, n_proc, info, h_max)
! Put max in h_o so that z_orog_print can calculate grid over max h_o
  h_o = h_max
  IF(me  ==  0.AND. PrintStatus >= PrStatus_Normal)THEN
    WRITE(6,fmt='(A,I2)')'surface_type = ',surface_type
    WRITE(6,fmt='(A)')' ** Orography in input dump used everywhere **'
    WRITE(6,fmt='(A,E16.8,A)')'maximum orographic height = ',h_o,' metres'
  END IF ! (me == 0)


!-----------------------------------------------------------------------
! Unsurported option
!-----------------------------------------------------------------------
ELSE

  IF(me  ==  0.AND. PrintStatus >= PrStatus_Normal)THEN
    WRITE(6,fmt='(A,I2,A)')' surface_type = ',surface_type,               &
                          ' option unavailable '
    CALL ereport('eg_idl_surface_setup',surface_type,                    &
                 'Invalid surface type')
  END IF ! (me == 0)

END IF        ! on surface_type


!-----------------------------------------------------------------------
! Swap bounds etc
!-----------------------------------------------------------------------

! DEPENDS ON: swap_bounds
CALL swap_bounds(orog_haloes,row_length,rows,1,                        &
                 halo_i,halo_j,fld_type_p,.FALSE.)

! DEPENDS ON: swap_bounds
CALL swap_bounds(orog_u,row_length,rows,1,                             &
                 halo_i,halo_j,fld_type_u,.FALSE.)

! DEPENDS ON: swap_bounds
CALL swap_bounds(orog_v,row_length,n_rows,1,                           &
                 halo_i,halo_j,fld_type_v,.FALSE.)

!  Copy orog_haloes into D1 array
DO j = 1, rows
  DO i = 1, row_length
    orography(i,j) = orog_haloes(i,j)
  END DO
END DO

IF (lhook) CALL dr_hook('EG_IDL_SURFACE_SETUP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_surface_setup
