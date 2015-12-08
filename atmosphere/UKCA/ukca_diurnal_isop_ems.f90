! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Description:
!  Main driver routine for chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_MAIN1.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Adapted from original code written by Guang Zeng
!
!------------------------------------------------------------------
!
      SUBROUTINE UKCA_DIURNAL_ISOP_EMS(row_length, &
                                       rows,       &
                                       emi_in,     &
                                       cosza_in,   &
                                       intza_in,   &
                                       sinlat,     &
                                       coslat,     &
                                       tanlat,     &
                                       timestep,   &
                                       emi_out,    &
                                       testdcycl)

      USE UKCA_CONSTANTS
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Control_Max_Sizes
      USE PrintStatus_mod
      USE Submodel_Mod

      IMPLICIT NONE

! CTIME ----------------------------------------------------
!
!  Purpose: Derived model time/step information including start/end
!           step numbers and frequencies (in steps) of interface field
!           generation, boundary field updating, ancillary field
!           updating; and assimilation start/end times.
!           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!           Also contains current time/date information, current
!           step number (echoed in history file) and steps-per-group.
!
!END -----------------------------------------------------------------

      INTEGER :: I_YEAR               ! Current model time (years)
      INTEGER :: I_MONTH              ! Current model time (months)
      INTEGER :: I_DAY                ! Current model time (days)
      INTEGER :: I_HOUR               ! Current model time (hours)
      INTEGER :: I_MINUTE             ! Current model time (minutes)
      INTEGER :: I_SECOND             ! Current model time (seconds)
      INTEGER :: I_DAY_NUMBER         ! Current model time (day no)
      INTEGER :: PREVIOUS_TIME(7)     ! Model time at previous step
      INTEGER :: IAU_DTResetStep      ! Data time reset step for IAU run

      INTEGER :: BASIS_TIME_DAYS  ! Integral no of days to basis time
      INTEGER :: BASIS_TIME_SECS  ! No of seconds-in-day at basis time

      LOGICAL :: L_C360DY

! UM6.5MODEL_ANALYSIS_HRS changed to REAL - 
!   requires FORECAST_HRS and DATA_MINUS_BASIS_HRS to REAL also 
      REAL    :: FORECAST_HRS     ! Hours since Data Time (ie T+nn)
      REAL    :: DATA_MINUS_BASIS_HRS ! Data time - basis time (hours)

      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
        I_DAY_NUMBER,PREVIOUS_TIME,                                     &
        BASIS_TIME_DAYS,BASIS_TIME_SECS,                                &
        FORECAST_HRS,DATA_MINUS_BASIS_HRS,                              &
        IAU_DTResetStep, L_C360DY

      INTEGER :: STEPim(INTERNAL_ID_MAX)  ! Step no since basis time
      INTEGER :: GROUPim(INTERNAL_ID_MAX) ! Number of steps per group

      ! Finish step number this run
      INTEGER :: TARGET_END_STEPim(INTERNAL_ID_MAX)

      REAL :: SECS_PER_STEPim(INTERNAL_ID_MAX) ! Timestep length in secs

      ! Frequency of interface field generation in steps
      INTEGER :: INTERFACE_STEPSim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Start steps for interface field generation
      INTEGER :: INTERFACE_FSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! End steps for interface field generation
      INTEGER :: INTERFACE_LSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Frequency of  updating boundary fields in steps
      INTEGER :: BOUNDARY_STEPSim(INTERNAL_ID_MAX)

      ! No of steps from boundary data prior to basis time to model
      ! basis time
      INTEGER :: BNDARY_OFFSETim(INTERNAL_ID_MAX)

      ! Lowest frequency for updating of ancillary fields in steps
      INTEGER :: ANCILLARY_STEPSim(INTERNAL_ID_MAX)

      ! Start steps for assimilation
      INTEGER :: ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX)

      ! Number of assimilation steps to analysis
      INTEGER :: ASSIM_STEPSim(INTERNAL_ID_MAX)

      ! Number of assimilation steps after analysis
      INTEGER :: ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)

      COMMON/CTIMEE/                                                    &
     &  STEPim,GROUPim,TARGET_END_STEPim,INTERFACE_STEPSim,             &
     &  INTERFACE_FSTEPim,INTERFACE_LSTEPim,BOUNDARY_STEPSim,           &
     &  BNDARY_OFFSETim,ANCILLARY_STEPSim,ASSIM_FIRSTSTEPim,            &
     &  ASSIM_STEPSim,ASSIM_EXTRASTEPSim,SECS_PER_STEPim

! CTIME end

      INTEGER,                              INTENT(IN)  :: row_length
      INTEGER,                              INTENT(IN)  :: rows
      LOGICAL,                              INTENT(IN)  :: testdcycl
      
      REAL,                                 INTENT(IN)  :: timestep
      REAL, DIMENSION(1:row_length,1:rows), INTENT(IN)  :: cosza_in    ! COS (zenith angle)
      REAL, DIMENSION(1:row_length,1:rows), INTENT(IN)  :: intza_in    ! INT(COS(sza))
      REAL, DIMENSION(1:row_length,1:rows), INTENT(IN)  :: sinlat      ! sin (latitude)
      REAL, DIMENSION(1:row_length,1:rows), INTENT(IN)  :: coslat      ! cos (latitude)
      REAL, DIMENSION(1:row_length,1:rows), INTENT(IN)  :: tanlat      ! tan (latitude)
      REAL, DIMENSION(1:row_length,1:rows), INTENT(IN)  :: emi_in      ! IN isoprene emission
      REAL, DIMENSION(1:row_length,1:rows), INTENT(OUT) :: emi_out     ! OUT diurnally varying isoprene emission

!     Local variables

      INTEGER :: i,j

      REAL, PARAMETER :: zerocheck       = 1.0E-6
      REAL, PARAMETER :: secs_per_hour   = 3600.0
      REAL, PARAMETER :: hours_per_day   = 24.0

      REAL :: declin                   ! Solar declination angle
      REAL :: int_a,b,int_h,sza_int    ! SZA integration variables
      REAL :: emit_day                 ! C5H8 emission in one day
      REAL :: trise                    ! Time of the sunrise
      REAL :: fxa,fxb,fxc
      REAL :: mpd, mod_mpd

      REAL, DIMENSION(1:row_length,1:rows) :: cosza       ! COS (zenith angle)
      REAL, DIMENSION(1:row_length,1:rows) :: intza       ! INT(COS(sza))
      REAL, DIMENSION(1:row_length,1:rows) :: daylen      ! day length for curent day
      REAL, DIMENSION(1:row_length,1:rows) :: cs_hour_ang ! cosine hour angle


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_DIURNAL_ISOP_EMS.F90',zhook_in,zhook_handle)

!     Calculate Declination Angle and Daylength for each row for
!     curent day of the year.
!     Ensure COS of HOUR ANGLE does not exceed + or - 1, and set DAY
!     LENGTH of 1st & last rows to be the same as that of adjacent rows
!     to avoid possible problems at the poles (tan(90)=infinity).

      fxa = Pi_Over_180
      fxb = 23.45/Recip_Pi_Over_180
      fxc = hours_per_day/pi

      daylen(:,:)      = 0.0
      cs_hour_ang(:,:) = 0.0
      emi_out(:,:)     = 0.0
     
      cosza(:,:)    = cosza_in(:,:)
      intza(:,:)    = intza_in(:,:)
           
      declin = fxb * SIN(fxa*(266.0+i_day_number))
      DO j = 1,rows
        DO i = 1,row_length
     
          cs_hour_ang(i,j) = -1.0*tanlat(i,j) * TAN(declin)
          IF (cs_hour_ang(i,j) < -1.0) cs_hour_ang(i,j)=-1.0
          IF (cs_hour_ang(i,j) > 1.0) cs_hour_ang(i,j)=1.0
          
          mpd     = (fxc * ACOS(cs_hour_ang(i,j)))*60.0 ! compute minutes of sunshine per day
          mod_mpd = MOD(mpd,(timestep/60.0)) ! compute residual mins (in excess to mins per timestep)
          daylen(i,j) = (mpd - mod_mpd)/60.0 ! compute sunshine hours a multiple of timestep
          
          trise = 12.0 - 0.5*daylen(i,j)
     
          int_a   = sinlat(i,j)*SIN(declin)*daylen(i,j)
          b       = coslat(i,j)*COS(declin)
          int_h   = (24.0/pi)*SIN(pi*trise/12.0)
          sza_int = int_a + b*int_h

!         adjust factor to model time step units

          sza_int = sza_int*secs_per_hour
          
!         calculate emission(day)

          emit_day = emi_in(i,j)*secs_per_hour*hours_per_day
    
!         now scale the emissions

          IF((cosza(i,j) > 0.0) .AND. (sza_int > 1.0E-1)) THEN
              emi_out(i,j) = emit_day*(cosza(i,j)/intza(i,j))
              IF (emi_out(i,j) < 0.0) THEN
                emi_out(i,j) = 0.0
              END IF
          ELSE
            emi_out(i,j) = 0.0
          END IF
     
          IF ((j == 1) .AND. (i == 1) .AND.               &
              (emi_out(i,j) > 0.0)     .AND.              &
              PrintStatus >= PrStatus_Diag .AND.          & 
              (testdcycl)) THEN
            
            WRITE(6,'(A8,2A5,6A15,3A12)')                 &
                    'UKCA_DIURNAL_ISOP_EMS:',             &
                    'i', 'j',                             &
                    'cont_daylen',                        &
                    'disc_daylen',                        &
                    'deg_lat',                            &
                    'SZA',                                &
                    'INT(SZA)',                           &
                    'INT(SZA) prec',                      &
                    'emi_in',                             &
                    'emi_out',                            &
                    'emi_day'
           
            WRITE(6,'(A8,2I5,6F15.1,3E12.4)')             &
                    i, j,                                 &
                    (fxc * ACOS(cs_hour_ang(i,j))),       &
                    daylen(i,j),                          &
                    ACOS(coslat(i,j))*Recip_Pi_Over_180,  &
                    ACOS(cosza(i,j))*Recip_Pi_Over_180,   &
                    sza_int,                              &
                    intza(i,j),                           &
                    emi_in(i,j),                          &
                    emi_out(i,j),                         &
                    emit_day
          END IF
        END DO                  !end domain loop
      END DO
      
      IF (lhook) CALL dr_hook('UKCA_DIURNAL_ISOP_EMS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ukca_diurnal_isop_ems
