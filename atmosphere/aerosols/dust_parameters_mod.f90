! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   A module containing constants/parameters used in the dust scheme
!
MODULE dust_parameters_mod

!
! Description:
!   This module contains declarations for constants and tunable 
!   parameters used to diagnose the emission and deposition of 
!   mineral dust
!  
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards version 8.2.
!

  USE missing_data_mod, ONLY: RMDI, IMDI
  USE ereport_mod, ONLY : ereport
IMPLICIT NONE


!
! Parameters fundamental to the dust scheme 
! (Details of dust properties/size divisions etc.)
! =======================================================================

! Prognostic mineral dust aerosol
LOGICAL :: l_dust  
! Diagnostic mineral dust aerosol lifting
LOGICAL :: l_dust_diag

! Logical to contain whether using two-bin (true) or six-bin (false) dust 
LOGICAL :: l_twobin_dust = .FALSE.
!   (default to six-bin dust)

! Define which dust bins are active
LOGICAL :: l_dust_div1, l_dust_div2, l_dust_div3, l_dust_div4,        &
           l_dust_div5, l_dust_div6

! Number of discrete particle size divisions (bins)
INTEGER :: ndiv              ! number of divisions that can be lifted 
                             !  from the surface
INTEGER, PARAMETER :: ndivh=9! number of divisions that can be blown 
                             !  horizontally along the surface and contribute 
                             !  to the lifting of the 1st NDIV divisions
INTEGER, PARAMETER :: ndivl=6! number of divisions in dust soil ancillaries

INTEGER :: i_dust = IMDI                    ! dust scheme setting in namelist
INTEGER, PARAMETER :: i_dust_off = 0        ! No dust
INTEGER, PARAMETER :: i_dust_prognostic = 1 ! Prognostic dust       
INTEGER, PARAMETER :: i_dust_diagnostic = 2 ! Diagnostic dust

! The size of particles included in each division 
! Note that by using two arrays for the max/min diameter here we can set
! up overlapping divisions. We must take care to make them consistent

REAL, ALLOCATABLE ::  dmax(:)
        ! max diameter of particles in each div.
REAL, ALLOCATABLE ::  dmin(:)
        ! min diameter of particles in each div.
REAL, ALLOCATABLE ::  drep(:)
        ! representative particle diameter
!
! Physical properties of the dust
REAL, PARAMETER :: rhop = 2.65E+3  ! density of a dust particle (quartz)

!
! Parameters used during dust emissions calculations
! =======================================================================
!
! Parameters based on observations/published research
REAL, PARAMETER, DIMENSION(ndivh) :: ustd_bas =                       &
       (/ 0.85, 0.72, 0.59, 0.46, 0.33, 0.16, 0.14, 0.18, 0.28 /)
                                      ! impact U*t derived from Bagnold (1941)
REAL, PARAMETER :: horiz_c = 2.61     ! C in horizontal flux calc (White 1979)
REAL, PARAMETER :: vert_a = 13.4      ! A in vertical flux calc(Gillette 1979)
REAL, PARAMETER :: vert_b = -6.       ! B in vertical flux calc(Gillette 1979)
REAL, PARAMETER :: vert_c = 0.01      ! cgs to si conversion

!
! Input variables needed when using the dust lifting scheme with 1 tile
!
REAL            :: z0_soil            ! Roughness length over bare soil

!
! Tuning parameters - set by namelist
REAL            :: us_am = RMDI       ! ustar correction (multiplic.)
REAL            :: sm_corr = RMDI     ! soil moist. correction factor
REAL            :: horiz_d = RMDI     ! Global tuning param for horizontal
                                      !  (and hence vertical) flux
! Tuning parameters - defined here
REAL, PARAMETER :: us_aa  = 0.0       ! ustar correction (additional)
REAL, PARAMETER :: ust_aa = 0.0       ! ustar_t correction (add)
REAL, PARAMETER :: ust_am = 1.0       ! ustar_t correction (multi.)
!
! Limits used in dust flux calculations
REAL            :: u_s_min = 1.e-5    ! Minimum val of u_s_std,
                                      !  below which no dust is produced
                                      !  (avoids divide by zero problem)
REAL, PARAMETER :: clay_max = 0.1     ! Max clay fraction.
REAL            :: snowmin = 1.e-6    ! Min snow depth for dust
REAL            :: h_orog_limit = 150.! 1/2pk-trough height above which
                                      !  no dust is produced
REAL, PARAMETER :: fland_lim =0.99999 ! No ems if fland<lim as windspeed
                                      !  too high at coastal points
!
! Switch to diagnose vertical flux using a fixed (user definable) size
!  distribution. If set to false, the vertical flux in each bin at a point 
!  is poportional to the horizontal flux in that bin at that point.
LOGICAL         :: l_fix_size_dist = .false.

! Namelist components can't be allocatable, so set to six, but in two-bin
! case should only use the first two array positions.
REAL :: size_dist(6) = RMDI

!
! Switch to allow emission of dust on tiles other than the bare soil tile
! 0 allows emission only on bare soil, 1 uses a bare soil radiative frac
! from the LAI, as for surface albedo (scaling by dust_veg_sc_nojules),
! and other methods (2+ could be added later).
INTEGER         :: dust_veg_emiss = IMDI
!
! Parameters used during the gravitational settling of dust
! =======================================================================
!

REAL, PARAMETER :: accf = 1.257 ! Cunningham correction factor term A
REAL, PARAMETER :: bccf = 0.4   ! Cunningham correction factor term B
REAL, PARAMETER :: cccf = -1.1  ! Cunningham correction factor term C

!
! Parameters used during the scavenging of dust
! =======================================================================
!
REAL, ALLOCATABLE :: krain_dust(:)
REAL, ALLOCATABLE :: ksnow_dust(:)
!
! RUN_Dust namelist via which non-parameter values herein can be set
! =======================================================================
! 
  NAMELIST /RUN_Dust/                                                   &
       us_am, sm_corr, horiz_d, l_fix_size_dist, dust_veg_emiss,                      &
       i_dust, l_twobin_dust

  CONTAINS
!
! Internal subroutine to check that entries to RUN_dust are consistent
! =======================================================================
! 
  SUBROUTINE dust_parameters_check( )
!
!   Check that a user-defined emissions size distribution is normalised
!
    IMPLICIT NONE
    
    IF (l_fix_size_dist) THEN
      size_dist(:)=size_dist(:)/SUM(size_dist)
    END IF
    
  END SUBROUTINE dust_parameters_check


! Internal subroutine to set the default value of size_dist. This may be
! overwritten by setting size_dist in the namelist (hence it's in a separate
! subroutine) or by setting l_twobin_dust. If l_fix_size_dist is set it is
! possible to override size_dist in the namelist as that was the previous
! behaviour, and it's mandatory to use it in the case of l_twobin_dust. If 
! l_fix_size_dist is not set the value here is irrelevant as it's not used

  SUBROUTINE dust_size_dist_initialise( )

      IMPLICIT NONE

      size_dist(1:6) =                                                   &
            (/ 0.0005, 0.0049, 0.0299, 0.2329, 0.4839, 0.2479 /)

  END SUBROUTINE dust_size_dist_initialise

! Internal subroutine to load the correct parameters into the variables in this
! module depending on whether two- or six- bin dust is being used
  SUBROUTINE dust_parameters_load( )

    IMPLICIT NONE

    INTEGER           :: icode           ! Error code
    CHARACTER(LEN=80) :: cmessage        ! Error message

! Convert i_dust into appropriate logical settings
    SELECT CASE(i_dust)
      CASE(i_dust_off)
        l_dust      = .FALSE.
        l_dust_diag = .FALSE.
      CASE(i_dust_prognostic)
        l_dust      = .TRUE.
        l_dust_diag = .FALSE.
      CASE(i_dust_diagnostic)
        l_dust      = .FALSE.
        l_dust_diag = .TRUE.
      CASE DEFAULT
        icode = 2
        WRITE(cmessage, '(A,I2,A)') 'i_dust = ', i_dust, ' is not a valid'&
                                  // ' setting'
        CALL ereport('DUST_PARAMETERS_LOAD', icode, cmessage)
    END SELECT

! Ascertain which dust bins are active    
    IF (i_dust == i_dust_prognostic) THEN
      l_dust_div1 = .TRUE.
      l_dust_div2 = .TRUE.
      IF (l_twobin_dust) THEN
        l_dust_div3 = .FALSE.
        l_dust_div4 = .FALSE.
        l_dust_div5 = .FALSE.
        l_dust_div6 = .FALSE.
      ELSE
        l_dust_div3 = .TRUE.
        l_dust_div4 = .TRUE.
        l_dust_div5 = .TRUE.
        l_dust_div6 = .TRUE.
      END IF    
    ELSE
      l_dust_div1 = .FALSE.
      l_dust_div2 = .FALSE.
      l_dust_div3 = .FALSE.
      l_dust_div4 = .FALSE.
      l_dust_div5 = .FALSE.
      l_dust_div6 = .FALSE.
    END IF

! Mandatory use of fixed size distribution when using two-bin code
    IF (l_twobin_dust) l_fix_size_dist = .TRUE.

! Define ndiv
    IF (l_twobin_dust) THEN
      ndiv = 2
    ELSE
      ndiv = 6
    END IF      

! Allocate scavenging co-efficients
    ALLOCATE(krain_dust(ndiv))
    ALLOCATE(ksnow_dust(ndiv))
    ALLOCATE(dmax(ndiv))
    ALLOCATE(dmin(ndiv))
    ALLOCATE(drep(ndiv))


! Set variables depending on how many bins
    IF (l_twobin_dust) THEN
      
 ! scav. coeff. for rain      
      krain_dust(1:ndiv) = (/ 4.0e-5, 3.0E-4 /)
  ! scav. coeff. for snow
      ksnow_dust(1:ndiv) = (/ 4.0E-5, 3.0E-4 /)

      dmax(1:ndiv) =  (/ 4.0E-6, 2.0E-5 /)
                ! max diameter of particles in each div.
      dmin(1:ndiv) =  (/ 2.0E-7, 4.0E-6 /)
                ! min diameter of particles in each div.
      drep(1:ndiv) =  (/ 2.0E-6, 8.0E-6 /)
                ! representative particle diameter

! Two-bin dust must have this size distribution
      size_dist(1:ndivl) =                                                &
          (/ 0.1800, 0.8200, 0.0000, 0.0000, 0.0000, 0.0000 /)

    ELSE ! Six-bin dust

 ! scav. coeff. for rain
      krain_dust(1:ndiv) =                                                &
       (/ 2.e-5, 2.e-5, 3.e-5, 6.e-5, 4.e-4, 4.e-4 /)
  ! scav. coeff. for snow
      ksnow_dust(1:ndiv) =                                                &
       (/ 2.e-5, 2.e-5, 3.e-5, 6.e-5, 4.e-4, 4.e-4 /)


      dmax(1:ndiv) =                                                      &
       (/ 2.0E-7, 6.32456E-7, 2.0E-6, 6.32456E-6, 2.0E-5, 6.32456E-5 /)
                ! max diameter of particles in each div.
      dmin(1:ndiv) =                                                      &
       (/ 6.32456E-8, 2.0E-7, 6.32456E-7, 2.0E-6, 6.32456E-6, 2.0E-5 /)
                ! min diameter of particles in each div.
      drep(1:ndiv) =                                                      &
                     (/ 0.112468E-06, 0.355656E-06, 0.112468E-05,         &   
                        0.355656E-05, 0.112468E-04, 0.355656E-04 /)
                ! representative particle diameter

      
    ENDIF ! l_twobin_dust
  
  END SUBROUTINE dust_parameters_load  


! internal subroutine to deallocate the dust arrays - provided for consistency,
! not actually called at present
  SUBROUTINE dust_parameters_unload( )

    IMPLICIT NONE
    DEALLOCATE (drep)
    DEALLOCATE (dmin)
    DEALLOCATE (dmax)
    DEALLOCATE (ksnow_dust)
    DEALLOCATE (krain_dust)
  
  END SUBROUTINE dust_parameters_unload

END MODULE dust_parameters_mod

