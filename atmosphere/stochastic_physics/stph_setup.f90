! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE stph_setup_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE stph_setup( )

! Generate and store or read Random Seed for STOCHASTIC PHYSICS
!  Driven by logicals: l_rp2, l_skeb2, 
!                    : l_stphseed_write, l_stphseed_read
!
! --------------------------------------------------------------------
! Sets up a seed that should vary from run to run, i.e. to ensure that
!  each ensemble member gets a different set of random numbers
!  Every RANDOM call after this changes the seed in a predictable way.
!
! Notes on the use of the random seed:
! If L_STPHSEED_WRITE=.TRUE., generate a new seed based on the actual
!  date and write the seed to a file for future reference.
!  At every dump step the random seed at that point and the random
!  streamfunction pattern coefficients overwrite this file.
!
! If L_STPHSEED_READ=.TRUE., read the random seed from the named file.
!  If CRUN is detected, then the streamfunction coefficients are also
!  read in from the seed file.
!
! If both L_STPHSEED_READ + L_STPHSEED_WRITE=.FALSE., use the model
!  date to generate the random seed.
!
! In all instances, the SEED value is written to standard output.
!
! Creates masks used for smoothing and damping SKEB2 fields if set
!  by logicals. These are saved in memory to reduce runtime.

! Main SKEB2 switch
 USE stochastic_physics_run_mod,  ONLY:                                 &
     l_skeb2, nsmooth, offx_stph, offy_stph, mask_smooth,               &
     mask_pdamp, l_skebsmooth_adv, l_rp2, l_stphseed_read,              &
     l_stphseed_write
! Value of pi
 USE conversions_mod, ONLY: pi

! ENDGame compatible array-bounds pointers
 USE atm_fields_bounds_mod, ONLY:                                       &
          pdims, stphdims_l

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim
 USE ereport_mod, ONLY : ereport
 USE PrintStatus_mod
 USE UM_ParVars
 USE Control_Max_Sizes
 USE domain_params
 USE stph_closeinput_mod,  ONLY: stph_closeinput
 USE stph_closeoutput_mod, ONLY: stph_closeoutput
 USE stph_openinput_mod,   ONLY: stph_openinput
 USE stph_openoutput_mod,  ONLY: stph_openoutput
 USE stph_readentry_mod,   ONLY: stph_readentry
 USE stph_writeentry_mod,  ONLY: stph_writeentry

 USE Submodel_Mod
 IMPLICIT NONE

! This include contains the current model date
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

! parameter for finding PE0
 INTEGER, PARAMETER ::  stph_zero = 0   ! Zero PE

! local Arrays for Stochastic Physics Random Number generation
 INTEGER :: iarg     ! random seed input argument
 INTEGER :: max_iarg ! maximum integer size of seed argument 
 INTEGER :: tam      ! size of the random seed
 INTEGER :: nens = 0 ! ensemble member number 
 INTEGER :: icode    ! error code from fort_get_env
 INTEGER :: dt(8)    ! date/time info for the random seed
 INTEGER :: i, j, n, ip1, im1, jp1, jm1
              
! random numbers for calculating seed
 REAL, ALLOCATABLE :: rnum(:)

! local arrays for smoothing mask
 REAL              :: r_dist   ! Radial distance for mask_pdamp
 REAL, ALLOCATABLE :: mask_smooth_tmp(:,:)

! local temporary arrays
 CHARACTER(LEN=256)       :: cmessage      ! out error message
 CHARACTER(LEN=8  )       :: ens_member    ! ensemble member (string)
 CHARACTER(LEN=*), PARAMETER  :: routinename='stph_setup'

! random seeds
 INTEGER, DIMENSION(:), ALLOCATABLE :: iranseed
 INTEGER, DIMENSION(:), ALLOCATABLE :: prevseed

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

! --------------------------------------------------------------------

 IF (lhook) CALL dr_hook('STPH_SETUP',zhook_in,zhook_handle)

! Keep a record of the settings for this run in the PE output files

 IF (printstatus  >=  prstatus_normal) THEN 
   WRITE(6,'(A)') ' '
   WRITE(6,'(A)') '*** STPH SETUP ***'
   WRITE(6,'(A,L1)') 'L_RP2 = ', l_rp2
   WRITE(6,'(A,L1)') 'L_SKEB2 = ', l_skeb2
   WRITE(6,'(A,L1)') 'L_STPHSEED_READ = ', l_stphseed_read
   WRITE(6,'(A,L1)') 'L_STPHSEED_WRITE = ', l_stphseed_write
 END IF

 IF (l_skebsmooth_adv) THEN
! Set Halo size for spatial smoothing to number of smoothing iterations
   offx_stph = nsmooth
   offy_stph = nsmooth
   stphdims_l%i_start = pdims%i_start - offx_stph
   stphdims_l%i_end   = pdims%i_end   + offx_stph
   stphdims_l%i_len   = stphdims_l%i_end - stphdims_l%i_start + 1
   stphdims_l%j_start = pdims%j_start - offy_stph
   stphdims_l%j_end   = pdims%j_end   + offy_stph
   stphdims_l%j_len   = stphdims_l%j_end - stphdims_l%j_start + 1
   stphdims_l%k_start = pdims%k_start
   stphdims_l%k_end   = pdims%k_end
   stphdims_l%k_len   = stphdims_l%k_end - stphdims_l%k_start + 1
   stphdims_l%halo_i  = offx_stph
   stphdims_l%halo_j  = offy_stph

! Setup pattern damping mask ( 0 in middle -> 1 at [offy_stph+1] )
   IF (.NOT.ALLOCATED(mask_pdamp))                                      &
     ALLOCATE(mask_pdamp(-offx_stph:offx_stph, -offy_stph:offy_stph))
   DO j = -offy_stph, offy_stph
     DO i = -offx_stph, offx_stph
       r_dist = pi * SQRT(REAL(i**2 + j**2))/REAL(offx_stph + 1)
       r_dist = MIN(r_dist, pi)   ! Maximum radius
       mask_pdamp(i,j) = 0.5 - 0.5 * COS(r_dist)
     END DO
   END DO
   IF (printstatus  ==  prstatus_diag) THEN 
     WRITE(6,'(11ES11.4)') mask_pdamp
   END IF

! Setup spatial smoothing mask
   IF (.NOT.ALLOCATED(mask_smooth))                                     &
     ALLOCATE(mask_smooth(-offx_stph-1:offx_stph+1,                     &
                          -offy_stph-1:offy_stph+1))
   IF (.NOT.ALLOCATED(mask_smooth_tmp))                                 &
     ALLOCATE(mask_smooth_tmp(-offx_stph-1:offx_stph+1,                 &
                            -offy_stph-1:offy_stph+1))
   DO j = -offy_stph-1, offy_stph+1
     DO i = -offx_stph-1, offx_stph+1
       mask_smooth(i,j) = 0.
       mask_smooth_tmp(i,j) = 0.
     END DO
   END DO
   mask_smooth(0,0) = 1.
   DO n = 1, nsmooth
     IF (printstatus  ==  prstatus_diag) THEN 
       WRITE(6,'(A,ES22.15)') 'SUM(mask_smooth)= ', SUM(mask_smooth)
       WRITE(6,'(13ES11.4)') mask_smooth
     END IF
! Keep smoothed values in temporary array
     DO j = -offy_stph, offy_stph
       jm1 = j - 1
       jp1 = j + 1
       DO i = -offx_stph, offx_stph
         im1 = i - 1
         ip1 = i + 1
         mask_smooth_tmp(i,j) = 0.0625*(mask_smooth(im1,jp1) +          &
                                        mask_smooth(ip1,jp1) +          &
                                        mask_smooth(im1,jm1) +          &
                                        mask_smooth(ip1,jm1) +          &
                                        2*( mask_smooth(i,jp1) +        &
                                            mask_smooth(im1,j) +        &
                                            mask_smooth(ip1,j) +        &
                                            mask_smooth(i,jm1) ) +      &
                                          4*mask_smooth(i,j) )

       END DO
     END DO
! Update main array
     DO j = -offy_stph-1, offy_stph+1
       DO i = -offx_stph-1, offx_stph+1
         mask_smooth(i,j) = mask_smooth_tmp(i,j)
       END DO
     END DO
   END DO

   IF (printstatus  ==  prstatus_diag) THEN 
     WRITE(6,'(A,ES22.15)') 'SUM(mask_smooth)= ', SUM(mask_smooth)
     WRITE(6,'(13ES11.4)') mask_smooth
   END IF

   DEALLOCATE(mask_smooth_tmp)

 END IF ! l_skebsmooth_adv


! IBM documentation suggests using generator=2 to improve
!  random cycling period and calculation speed
! (http://publib.boulder.ibm.com/infocenter/comphelp/v8v101/
!  .......index.jsp?topic=/com.ibm.xlf101a.doc/xlfcr/runpgms.htm)
! Get size of random seed for read/write/allocate
 CALL random_seed(SIZE=tam)

 IF  (mype == stph_zero) THEN
   IF (l_stphseed_read) THEN
   ! ---------------------------------------
   ! read in random seed from file.
   ! ---------------------------------------

     IF (.NOT.ALLOCATED(iranseed)) ALLOCATE(iranseed(tam))

     CALL stph_openinput(.false.)

     CALL stph_readentry(iranseed,tam)

     CALL stph_closeinput()

     CALL random_seed(PUT=iranseed(1:tam))
     IF  (printstatus  >=  prstatus_normal) THEN
       WRITE(6,'(A)') 'STPH SETUP: random seed read from file'
       WRITE(6,'(I30)') iranseed
     ENDIF

     DEALLOCATE(iranseed)
   ELSE IF (l_stphseed_write) THEN
     ! -----------------------------------------------------------
     ! generate new random seed using actual date YYMMDDHHmmSSsss
     ! -----------------------------------------------------------
     CALL date_and_time(VALUES=dt)

     IF (.NOT.ALLOCATED(prevseed)) ALLOCATE(prevseed(tam))
     IF (.NOT.ALLOCATED(iranseed)) ALLOCATE(iranseed(tam))
     IF (.NOT.ALLOCATED(rnum)) ALLOCATE(rnum(tam))

     ! Use full date in randomising seed to reduce chance of
     !  recycling from one run to the next.
     ! The formula calculates the days since ~2000AD and adds
     !  in time suitably inflated to fully change the seed.
     ! Only use last two digits of year to prevent numerical
     !  overflow at some date in the future.
     ! A random number generated from this seed is used to
     !  multiply the seed again
     dt(1) = dt(1) - 100*INT(0.01*dt(1))
     iarg = (dt(3) - 32075 +                                            &
              1461*(dt(1) + 4800 + (dt(2) - 14)/12)/4 +                 &
              367*(dt(2) - 2 - (dt(2)-14)/12*12)/12 -                   &
              3*((dt(1)+4900+(dt(2)-14)/12)/100)/4)*1000 +              &
              dt(8)**2.86 + dt(7)**3.79 + dt(5)**5.12 +                 &
              dt(6)**3.24
     ! Constrain iarg in a range to prevent numerical overflow
     ! 2**8 < iarg < SQRT(HUGE(INT)), since iarg is raised to n<2
     max_iarg = FLOOR(SQRT(1.0*HUGE(iarg)))
     iarg = MOD(iarg, max_iarg)
     iarg = MAX(iarg, 256)
     ! Generate initial random number set to use as power
     prevseed(:) = iarg
     CALL random_seed(PUT=prevseed(1:tam))
     CALL random_number(rnum)

     ! Range of seed from 0 to 2**31 (32-bit Int)
     iranseed(:)=iarg*rnum(:)
     ! Set final seed
     CALL random_seed(PUT=iranseed(1:tam))

     IF  (printstatus  >=  prstatus_normal) THEN
       WRITE(6,'(A)') 'STPH SETUP: Size of Random Seed'
       WRITE(6,'(I10)') tam
       WRITE(6,'(A)') 'STPH SETUP: Random seed using computer date/time'
       WRITE(6,'(I30)') iranseed
     END IF

     ! Keep record of seed used for failure re-runs
     CALL stph_openoutput("rewind    ",.false.)

     CALL stph_writeentry(iranseed,tam)

     CALL stph_closeoutput()

     DEALLOCATE(prevseed) 
     DEALLOCATE(iranseed) 
     DEALLOCATE(rnum) 
   ELSE        ! No seed file
     ! -----------------------------------------------------------
     ! Standard MOGREPS operational/trial setup using model date and
     ! ENVIRONMENT variable ENS_MEMBER
     ! -----------------------------------------------------------
     ! Get ENS_MEMBER from ENV
     ens_member = '00000000'
     CALL fort_get_env('ENS_MEMBER',10,ens_member,8,icode)
     ens_member = TRIM(ens_member)
     IF (icode == 0) THEN
       WRITE(6,'(A,A8)') 'Successfully retrieved ens_member = ', ens_member
     ELSE
       ! Force abort error if ENS_MEMBER cannot be read
       WRITE(6,'(A)') '**ERROR**: SKEB2 Initial Random Seed'
       WRITE(6,'(A)') '  Problem retrieving ENS_MEMBER from ENVIRONMENT'
       WRITE(6,'(A)') '  Section 35: check UMUI settings'
       WRITE(6,'(A)') '  Use Panel: Script Inserts or Modifications'
       WRITE (cmessage,'(A)') 'STPH_SETUP: Cannot retrieve environment' &
                              //' variable ENS_MEMBER.'
       
       CALL ereport(routinename, ABS(icode), cmessage)
     END IF

     ! Convert STRING to INTEGER value
     READ(ens_member,'(I8)') nens

     IF (nens < 0) THEN
       ! Force abort error if ENS_MEMBER is negative
       WRITE(6,'(A)') '**ERROR**: SKEB2 Initial Random Seed'
       WRITE(6,'(A)') '  User provided ENS_MEMBER cannot be negative'
       WRITE(6,'(A)') '  Section 35: check UMUI settings'
       WRITE(6,'(A)') '  Use Panel: Script Inserts or Modifications'
       WRITE (cmessage,'(A)') 'STPH_SETUP: Environment variable'        &
                              //' ENS_MEMBER cannot be negative.'
       
       CALL ereport(routinename, 1, cmessage)

     ELSE IF (nens > 0) THEN
       ! Stochastic physics is not used if ens_member is set to zero
       !  in the case of no seed file being used. This is set later
       !  on all PEs, so need to bcast ens_member (nens) here from PE0

       dt(1) = i_year
       dt(2) = i_month
       dt(3) = i_day
       dt(4) = 0                ! Shift from UTC
       dt(5) = i_hour + 1       ! Set range 1 - 24
       dt(6) = i_minute
       dt(7) = nens             ! Ens mem 0-24 into second dimension
       dt(8) = nens + 100       ! Ens mem 100-124 into millisecond dim

       IF (.NOT.ALLOCATED(prevseed)) ALLOCATE(prevseed(tam))
       IF (.NOT.ALLOCATED(iranseed)) ALLOCATE(iranseed(tam))
       IF (.NOT.ALLOCATED(rnum)) ALLOCATE(rnum(tam))

       ! Use full date in randomising seed to reduce chance of
       !  recycling from one run to the next.
       ! The formula calculates the days since ~2000AD and adds
       !  in time suitably inflated to fully change the seed.
       ! Only use last two digits of year to prevent numerical
       !  overflow at some date in the future.
       ! A random number generated from this seed is used to
       !  multiply the seed again
       dt(1) = dt(1) - 100*INT(0.01*dt(1))
       iarg = (dt(3) - 32075 +                                          &
              1461*(dt(1) + 4800 + (dt(2) - 14)/12)/4 +                 &
              367*(dt(2) - 2 - (dt(2)-14)/12*12)/12 -                   &
              3*((dt(1)+4900+(dt(2)-14)/12)/100)/4)*1000 +              &
              dt(8)**2.86 + dt(7)**3.79 + dt(5)**5.12 +                 &
              dt(6)**3.24
       ! Constrain iarg in a range to prevent numerical overflow
       ! 2**8 < iarg < SQRT(HUGE(INT)), since iarg is raised to n<2
       max_iarg = FLOOR(SQRT(1.0*HUGE(iarg)))
       iarg = MOD(iarg, max_iarg)
       iarg = MAX(iarg, 256)
       ! Generate initial random number set to use as power
       prevseed(:) = iarg
       CALL random_seed(PUT=prevseed(1:tam))
       CALL random_number(rnum)
       ! Range of seed from 0 to 2**31 (32-bit Int)
       iranseed(:)=iarg*rnum(:)
       ! Set final seed
       CALL random_seed(PUT=iranseed(1:tam))

       IF  (printstatus  >=  prstatus_normal) THEN
         WRITE(6,'(A)') 'STPH SETUP: Size of Random Seed, INT(mem)'     &
                        //', CHAR(mem)'
         WRITE(6,'(2I10,2X,A8)') tam, nens, ens_member
         WRITE(6,'(A,6I4,2I10)') 'STPH SETUP: model date = ',dt
         WRITE(6,'(A)') 'STPH SETUP: Random seed generated using model date'
         WRITE(6,'(I30)') iranseed
       END IF

     END IF    ! nens > 0
   END IF ! .NOT.l_stphseed_write .NOR.l_stphseed_read
 END IF ! mype == stph_zero  

 ! Override stochastic physics logicals to FALSE when control ens_member
 !  is set and the job is not using any seed files
 CALL gc_ibcast(1,1,stph_zero,nproc,icode,nens)
 IF (icode /= 0) THEN
   WRITE (cmessage,'(A)') 'STPH_SETUP: ens_member not bcast correctly'
   CALL ereport(routinename, ABS(icode), cmessage)
 END IF

 IF (.NOT.l_stphseed_write .AND. .NOT.l_stphseed_read .AND. nens == 0) THEN
       l_rp2   = .FALSE.
       l_skeb2 = .FALSE. 
       IF (printstatus  >=  prstatus_min) THEN 
         WRITE(6,'(A)')    '**********************************'
         WRITE(6,'(A)')    '***** STPH SETUP :: OVERRIDE *****'
         WRITE(6,'(A,I4)') '**     Ensemble Member = ', nens
         WRITE(6,'(A,L1)') '**     L_RP2   =  ', l_rp2
         WRITE(6,'(A,L1)') '**     L_SKEB2 =  ', l_skeb2
         WRITE(6,'(A)')    '**********************************'
       END IF
 END IF

 
 IF (lhook) CALL dr_hook('STPH_SETUP',zhook_out,zhook_handle)
 RETURN

END SUBROUTINE stph_setup
END MODULE stph_setup_mod
