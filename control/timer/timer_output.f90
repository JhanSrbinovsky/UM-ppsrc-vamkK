! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!   SUBROUTINE TIMER_OUTPUT -----------------------------------------
!
!                      Purpose:
!
!   Reports timing information from timer.

!*********************************************************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Timer

SUBROUTINE timer_output(                                                       &
  in_number_of_timers, ni_number_of_timers,                                    &
  in_cpu_time_elapsed, ni_cpu_time_elapsed,                                    &
  in_wallclock_time_elapsed, ni_wallclock_time_elapsed,                        &
  in_number_of_times_timed, ni_number_of_times_timed,                          &
  in_timer_name, ni_timer_name,                                                &
  action,message)

USE IOS_Common, ONLY :                                                         &
    L_IOS_Active,                                                              &
    io_servers,                                                                &
    l_io_server,                                                               &
    global_procs,                                                               &
    IOS_Server_Groups,                                                         &
    IOS_tasks_per_server,                                                      &
    model_procs
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
!$ USE omp_lib                 ! Note OpenMP sentinel

USE ereport_mod, ONLY : ereport
USE UM_ParVars
IMPLICIT NONE

! Arguments

INTEGER                                                                        &
  in_number_of_timers                                                          &
                       ! IN number of inclusive timers
, ni_number_of_timers                                                          &
                       ! IN number of non-inclusive timers
, in_number_of_times_timed(in_number_of_timers)                                &
!                            ! IN number of times timed - inclusive
, ni_number_of_times_timed(ni_number_of_timers)                                &
!                            ! IN number of times timed - non-incl.
, action  ! final output or intermediate

REAL                                                                           &
  in_cpu_time_elapsed(in_number_of_timers)                                     &
!                            ! IN elapsed inclusive CPU time
, ni_cpu_time_elapsed(ni_number_of_timers)                                     &
!                            ! IN elapsed non-inclusive CPU time
, in_wallclock_time_elapsed(in_number_of_timers)                               &
!                            ! IN elapsed inclusive wallclock time
, ni_wallclock_time_elapsed(ni_number_of_timers)

CHARACTER(LEN=20)                                                              &
  in_timer_name(in_number_of_timers)                                           &
!                            ! IN name of timed section - inclusive
, ni_timer_name(ni_number_of_timers)
!                            ! IN name of timed section - non-incl.

CHARACTER(LEN=*)                                                               &
   message                   ! IN message to print


! Local variables

INTEGER max_timers
PARAMETER(max_timers=300)

INTEGER last_call_to_timer,intermediate_output
PARAMETER(last_call_to_timer=2,                                                &
          intermediate_output=7)

INTEGER                                                                        &
  number_of_timers                                                             &
, local_number_of_times_timed(max_timers)

REAL                                                                           &
  local_cpu_time_elapsed(max_timers)                                           &
, local_wallclock_time_elapsed(max_timers)

CHARACTER(LEN=20)                                                                   &
  local_timer_name(max_timers)

! Variables required for using intermediate timers
! They record the values on the last call to this routine
INTEGER                                                                        &
  last_in_number_of_times_timed(max_timers)                                    &
, last_ni_number_of_times_timed(max_timers)

REAL                                                                           &
  last_in_cpu_time_elapsed(max_timers)                                         &
, last_ni_cpu_time_elapsed(max_timers)                                         &
, last_in_wallclock_time_elapsed(max_timers)                                   &
, last_ni_wallclock_time_elapsed(max_timers)

LOGICAL                                                                        &
  first_intermediate_timer_call

DATA first_intermediate_timer_call /.TRUE./
DATA last_in_number_of_times_timed /max_timers*0/
DATA last_ni_number_of_times_timed /max_timers*0/
DATA last_in_cpu_time_elapsed /max_timers*0.0/
DATA last_ni_cpu_time_elapsed /max_timers*0.0/
DATA last_in_wallclock_time_elapsed /max_timers*0.0/
DATA last_ni_wallclock_time_elapsed /max_timers*0.0/

SAVE                                                                           &
  last_in_number_of_times_timed, last_ni_number_of_times_timed,                &
  last_in_cpu_time_elapsed, last_ni_cpu_time_elapsed,                          &
  last_in_wallclock_time_elapsed, last_ni_wallclock_time_elapsed,              &
  first_intermediate_timer_call



INTEGER sortwork_int    ! work variable for sort
REAL    sortwork_real   ! work variable for sort
CHARACTER(LEN=20) sortwork_char ! work variable for sort

REAL total_cpu_time,                                                           &
                             ! total cpu time spent in program
     total_wallclock_time,                                                     &
                             ! total wallclock time spent in
!                                  ! program
     average_cpu_elapsed,                                                      &
                             ! average cpu elapsed time
     average_wallclock_elapsed,                                                &
                                ! average wallclock elapsed time
     percent_of_cpu_total,                                                     &
                             ! % of cpu time spent in a section
     percent_of_wallclock_total,                                               &
                            ! % of wallclock time spent in a
!                                 ! section
     speed_up               ! speed_up=cpu/wallclock

! These are the declarations for MPP timer

INTEGER info,                                                                  &
  wallclock_max_pe(max_timers),wallclock_min_pe(max_timers),                   &
  cpu_max_pe(max_timers),cpu_min_pe(max_timers)

REAL wallclock_mean(max_timers),cpu_mean(max_timers),                          &
     wallclock_median(max_timers),cpu_median(max_timers),                      &
     wallclock_sd(max_timers),cpu_sd(max_timers),                              &
     wallclock_max(max_timers),wallclock_min(max_timers),                      &
     cpu_max(max_timers),cpu_min(max_timers),                                  &
     cpu_total(max_timers)

INTEGER                                                                        &
  summ_n_timers                                                                &
                   ! number of routines for ni summary
, routine_id          ! routine id on this processor

REAL, ALLOCATABLE :: wallclock_times(:)! wallclock time from each proc
REAL, ALLOCATABLE :: cpu_times(:)              ! cpu time from each proc
REAL              :: total_cpu,max_wall        ! total cpu, maxumum wallclock times

CHARACTER(LEN=20) summ_section(max_timers)  ! names of sections
COMMON /mpp_timer/ summ_n_timers,total_cpu,max_wall,   &
                   summ_section

! Variables for loops etc.
INTEGER i,j,k,timer_kind

CHARACTER(LEN=8) env_num_threads ! to get environment variable OMP_NUM_THREADS
CHARACTER(LEN=80) cmessage    ! Internal error message
INTEGER err                 ! err return of fort_env_get subroutine
INTEGER atm_num_threads, omp_max_threads     ! to store value read from env_num_threads
INTEGER :: ICODE
CHARACTER(LEN=*) routinename
PARAMETER (routinename = 'TIMER_OUTPUT')

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Check to see if this is an intermediate output, and the first
! time it has been called
IF (lhook) CALL dr_hook('TIMER_OUTPUT',zhook_in,zhook_handle)

ALLOCATE(wallclock_times(0:nproc_max-1))
ALLOCATE(cpu_times(0:nproc_max-1))      


IF (L_IO_Server) THEN 
  ! IOS makes no use of timer.
  WRITE(6,'(A,A)')'timer_output: This IO Server is not ', &
      'reporting timing information'
ELSE
IF ((action  ==  intermediate_output) .AND.                                    &
    (first_intermediate_timer_call  ) ) THEN
! Copy the arguments into the last_* arrays
  first_intermediate_timer_call=.FALSE.

  DO i=1,in_number_of_timers
    last_in_number_of_times_timed(i)=in_number_of_times_timed(i)
    last_in_cpu_time_elapsed(i)=in_cpu_time_elapsed(i)
    last_in_wallclock_time_elapsed(i)=                                         &
      in_wallclock_time_elapsed(i)
  END DO

  DO i=1,ni_number_of_timers
    last_ni_number_of_times_timed(i)=ni_number_of_times_timed(i)
    last_ni_cpu_time_elapsed(i)=ni_cpu_time_elapsed(i)
    last_ni_wallclock_time_elapsed(i)=                                         &
      ni_wallclock_time_elapsed(i)
  END DO

  GO TO 9999  ! jump to end - no output on first call
END IF

WRITE(6,*)
WRITE(6,*) '******************************************'
WRITE(6,*)

DO timer_kind=1,2  ! 1 is non-inclusive and 2 is inclusive
! Copy arguments into local arrays
  IF (action  ==  last_call_to_timer) THEN
    WRITE(6,*) 'END OF RUN - TIMER OUTPUT'
    WRITE(6,*) 'Timer information is for whole run'
    IF (timer_kind  ==  1) THEN  ! non-inclusive timer
      number_of_timers=ni_number_of_timers
      DO i=1,number_of_timers
        local_timer_name(i)=ni_timer_name(i)
        local_cpu_time_elapsed(i)=ni_cpu_time_elapsed(i)
        local_wallclock_time_elapsed(i)=                                       &
          ni_wallclock_time_elapsed(i)
        local_number_of_times_timed(i)=                                        &
          ni_number_of_times_timed(i)
      END DO
    ELSE ! timer_kind  ==  2 - inclusive timer
      number_of_timers=in_number_of_timers
      DO i=1,number_of_timers
        local_timer_name(i)=in_timer_name(i)
        local_cpu_time_elapsed(i)=in_cpu_time_elapsed(i)
        local_wallclock_time_elapsed(i)=                                       &
          in_wallclock_time_elapsed(i)
        local_number_of_times_timed(i)=                                        &
          in_number_of_times_timed(i)
      END DO
    END IF ! which timer kind this was
  ELSE  ! this is an intermediate output call
    WRITE(6,*) 'INTERMEDIATE TIMER OUTPUT :',message
    WRITE(6,*) 'Timer information is only for code executed ',                 &
               'since last intermediate timer output.'
    IF (timer_kind  ==  1) THEN  ! non-inclusive timer
      number_of_timers=ni_number_of_timers
      DO i=1,number_of_timers
        local_timer_name(i)=ni_timer_name(i)
        local_cpu_time_elapsed(i)=ni_cpu_time_elapsed(i)-                      &
                                  last_ni_cpu_time_elapsed(i)
        local_wallclock_time_elapsed(i)=                                       &
          ni_wallclock_time_elapsed(i)-                                        &
          last_ni_wallclock_time_elapsed(i)
        local_number_of_times_timed(i)=                                        &
          ni_number_of_times_timed(i)-                                         &
          last_ni_number_of_times_timed(i)
      END DO
    ELSE ! timer kind  ==  2 - inclusive timer
      number_of_timers=in_number_of_timers
      DO i=1,number_of_timers
        local_timer_name(i)=in_timer_name(i)
        local_cpu_time_elapsed(i)=in_cpu_time_elapsed(i)-                      &
                                  last_in_cpu_time_elapsed(i)
        local_wallclock_time_elapsed(i)=                                       &
          in_wallclock_time_elapsed(i)-                                        &
          last_in_wallclock_time_elapsed(i)
        local_number_of_times_timed(i)=                                        &
          in_number_of_times_timed(i)-                                         &
          last_in_number_of_times_timed(i)
      END DO
    END IF  ! what timer type
  END IF  ! what action to perform

! Do work for non-inclusive timers

! Calculate the total time in the program (based on non-inclusive
! timers)
  IF (timer_kind  ==  1) THEN
    total_cpu_time = 0.0
    total_wallclock_time = 0.0
    DO i=1,number_of_timers
    total_cpu_time = total_cpu_time + local_cpu_time_elapsed(i)
    total_wallclock_time =                                                     &
      total_wallclock_time + local_wallclock_time_elapsed(i)
    END DO

    WRITE(6,*) 'PE ',mype,' Elapsed CPU Time: ',                               &
               total_cpu_time
    WRITE(6,*) 'PE ',mype,'  Elapsed Wallclock Time: ',                        &
                total_wallclock_time

! Calculate the total cpu time over all processors and the
! maximum elapsed time - so allowing a speedup to be caclulated

    total_cpu=total_cpu_time
    max_wall=total_wallclock_time

    CALL gc_rsum(1,nproc,info,total_cpu)
    CALL gc_rmax(1,nproc,info,max_wall)

    max_wall=MAX(max_wall,0.000001)
    WRITE(6,*)
    WRITE(6,*) 'Total Elapsed CPU Time: ',                                     &
               total_cpu
    WRITE(6,*) 'Maximum Elapsed Wallclock Time: ',                             &
               max_wall
    WRITE(6,*) 'Speedup: ',total_cpu/max_wall
    WRITE(6,*) '--------------------------------------------'

  END IF

! Sort subroutines into time order (based on wallclock time)

  DO i=1,number_of_timers-1
    DO j=(i+1),number_of_timers
      IF (local_wallclock_time_elapsed(j)  >                                   &
          local_wallclock_time_elapsed(i)) THEN

!             Swap the two entries
        sortwork_real = local_cpu_time_elapsed(i)
        local_cpu_time_elapsed(i) = local_cpu_time_elapsed(j)
        local_cpu_time_elapsed(j) = sortwork_real

        sortwork_real = local_wallclock_time_elapsed(i)
        local_wallclock_time_elapsed(i) =                                      &
          local_wallclock_time_elapsed(j)
        local_wallclock_time_elapsed(j) = sortwork_real

        sortwork_int = local_number_of_times_timed(i)
        local_number_of_times_timed(i) =                                       &
          local_number_of_times_timed(j)
        local_number_of_times_timed(j) = sortwork_int

        sortwork_char = local_timer_name(i)
        local_timer_name(i) = local_timer_name(j)
        local_timer_name(j) = sortwork_char

      END IF
    END DO
  END DO

  IF (timer_kind  ==  1) THEN
    WRITE(6,'(20x,a,i6)') 'Non-Inclusive Timer Summary for PE ',mype
    WRITE(6,'(3x,a,14x,a,2x,a,4x,a,3x,a,2x,a,2x,a,4x,a,4x,a)')                 &
      'ROUTINE','CALLS','TOT CPU','AVERAGE','TOT WALL','AVERAGE',              &
      '% CPU','% WALL','SPEED-UP'
  ELSE
    WRITE(6,'(20x,a,i6)') 'Inclusive Timer Summary for PE ',mype
    WRITE(6,'(3x,a,14x,a,2x,a,4x,a,3x,a,2x,a,2x,a)')                           &
      'ROUTINE','CALLS','TOT CPU','AVERAGE','TOT WALL',                        &
                                  'AVERAGE','SPEED-UP'
  END IF

  DO i=1,number_of_timers
    IF (local_number_of_times_timed(i)  /=  0) THEN
      average_cpu_elapsed =  local_cpu_time_elapsed(i)/                        &
                             local_number_of_times_timed(i)
      average_wallclock_elapsed = local_wallclock_time_elapsed(i)/             &
                                  local_number_of_times_timed(i)
    ELSE
       average_cpu_elapsed = 0.0
       average_wallclock_elapsed = 0.0
    END IF

    IF (local_wallclock_time_elapsed(i)  >   0) THEN
      speed_up=local_cpu_time_elapsed(i)/                                      &
               local_wallclock_time_elapsed(i)
    ELSE
      speed_up=1.0
    END IF

    IF (timer_kind  ==  1) THEN  ! non-inclusive timer has some
!                                      ! extra output

      percent_of_cpu_total = 100.0*local_cpu_time_elapsed(i)/                  &
                             total_cpu_time
      percent_of_wallclock_total =                                             &
        100.0*local_wallclock_time_elapsed(i)/                                 &
        total_wallclock_time

      WRITE(6,'(i3,1x,a20,1x,i4,4(2x,f8.2),2(2x,f6.2),4x,f6.2)')               &
                  i,local_timer_name(i),                                       &
                  local_number_of_times_timed(i),                              &
                  local_cpu_time_elapsed(i),average_cpu_elapsed,               &
                  local_wallclock_time_elapsed(i),                             &
                  average_wallclock_elapsed,                                   &
                  percent_of_cpu_total,                                        &
                  percent_of_wallclock_total,speed_up

    ELSE ! inclusive timer has slightly less to output

      WRITE(6,'(i3,1x,a20,1x,i4,4(2x,f8.2),4x,f6.2)')                          &
                  i,local_timer_name(i),                                       &
                  local_number_of_times_timed(i),                              &
                  local_cpu_time_elapsed(i),average_cpu_elapsed,               &
                  local_wallclock_time_elapsed(i),                             &
                  average_wallclock_elapsed,speed_up

    END IF

  END DO





! We only want to process the statistics of the timings for real mpp
! jobs where nproc > 1.  Any utilities where UTILIO, FLDIO, UTILHIST and
! FLUXPROC are defined for the cpp use nproc = 1 as parameter in parvars.h.

! And now to assemble an overall timing assesment on PE0
! Each PE sends it total wallclock and cpu time spent in each routine
! to PE0, which calculates the average, s.d., max and min, and
! sorts on the basis of the average wallclock time


! We'll use the list of routines that PE0 already has as the master
! list.

  IF (mype  ==  0) THEN
    WRITE(6,*)
    WRITE(6,'(A)') 'MPP Timing information : '
    WRITE(6,'(I6,A,I4,A,I4)')nproc,' processors in atmosphere configuration ' &
        ,nproc_x,' x ',nproc_y
    IF (L_IOS_Active()) THEN
      WRITE(6,'(I4,A,I2,A,I2,A)')SIZE(io_servers),                              & 
        ' IO servers in configuration '                                         &
        ,IOS_Server_Groups,' x ',IOS_tasks_per_server,                          &
        ' are not included in stats'
    ELSE IF (global_procs /= model_procs) THEN
      WRITE(6,'(I4,A)')global_procs-model_procs, &
          ' processes are not included in stats (misconfigured IO servers?)'
    END IF
! write out number of threads, warn user if thread number has been 
! changed during run of um
    CALL FORT_GET_ENV('OMP_NUM_THREADS',15,env_num_threads,8,err)

! check for intel's KMP_NUM_THREADS
    IF (err /= 0) THEN 
      CALL FORT_GET_ENV('KMP_NUM_THREADS',15,env_num_threads,8,err)
    ENDIF

    IF (err  ==  0) THEN

      READ(env_num_threads,'(I4)') atm_num_threads
        omp_max_threads = atm_num_threads

! Only want OpenMP section executing if OpenMP is compiled in,
! so protect by sentinal
!$ omp_max_threads = omp_get_max_threads()
      
      IF( omp_max_threads /= atm_num_threads) THEN
          
!$      cmessage = 'Environment variable '//                                   &
!$        'OMP_NUM_THREADS does not match current number of threads'
!$        ICODE = -100

!$        CALL ereport(routinename,ICODE,cmessage)

!$        WRITE(6,*) 'Number of threads : ', omp_max_threads
!$        WRITE(6,*) 'Environment variable OMP_NUM_THREADS : ',                &
!$               atm_num_threads
      ELSE
!$        WRITE(6,*) 'Number of OMP threads : ', atm_num_threads
      END IF

   END IF


    summ_n_timers=number_of_timers
    DO i=1,summ_n_timers
      summ_section(i)=local_timer_name(i)
    END DO
  END IF

! tell everyone else how many routines to do summary on - and which
! routines they are
  CALL gc_ibcast(3213,1,0,nproc,info,summ_n_timers)
  CALL gc_cbcast(3214,20*summ_n_timers,0,nproc,info,                           &
                 summ_section)


  DO i=1,summ_n_timers

    CALL gc_gsync (nproc,info)

! which section_ref is this for me?

    routine_id=0
    DO j=1,number_of_timers
      IF (local_timer_name(j)  ==  summ_section(i))                            &
        routine_id=j
    END DO

    IF (routine_id  >   0) THEN
      wallclock_times(mype)=                                                   &
        local_wallclock_time_elapsed(routine_id)
      cpu_times(mype)=local_cpu_time_elapsed(routine_id)
    ELSE
      wallclock_times(mype)=0.0
      cpu_times(mype)=0.0
    END IF

! send my information to PE 0.
    CALL gc_rsend(1000+mype,1,0,info,wallclock_times(mype),                    &
                wallclock_times(mype))
    CALL gc_gsync(nproc,info)

    IF (mype  ==  0) THEN
      DO j=0,nproc-1
        CALL gc_rrecv(1000+j,1,j,info,wallclock_times(j),                      &
                      wallclock_times(j))
      END DO
    END IF
    CALL gc_gsync(nproc,info)

    CALL gc_rsend(10000+mype,1,0,info,cpu_times(mype),                         &
                cpu_times(mype))
    CALL gc_gsync(nproc,info)

    IF (mype  ==  0) THEN
      DO j=0,nproc-1
        CALL gc_rrecv(10000+j,1,j,info,cpu_times(j),                           &
                      cpu_times(j))
      END DO
    END IF

    IF (mype  ==  0) THEN
! collect all the information - and start calculating the statistics
      wallclock_mean(i)=0.0
      cpu_total(i)=0.0
      wallclock_max(i)=-1.0e30
      wallclock_min(i)=1.0e30
      cpu_max(i)=-1.0e30
      cpu_min(i)=1.0e30

      DO j=0,nproc-1

        wallclock_mean(i)=wallclock_mean(i)+wallclock_times(j)
        cpu_total(i)=cpu_total(i)+cpu_times(j)

        IF (wallclock_times(j) >  wallclock_max(i)) THEN
          wallclock_max(i)=wallclock_times(j)
          wallclock_max_pe(i)=j
        END IF
        IF (wallclock_times(j) <  wallclock_min(i)) THEN
          wallclock_min(i)=wallclock_times(j)
          wallclock_min_pe(i)=j
        END IF
        IF (cpu_times(j) >  cpu_max(i)) THEN
          cpu_max(i)=cpu_times(j)
          cpu_max_pe(i)=j
        END IF
        IF (cpu_times(j) <  cpu_min(i)) THEN
          cpu_min(i)=cpu_times(j)
          cpu_min_pe(i)=j
        END IF

      END DO ! loop over processors

! and calculate the statistics
! first calculate the means
      wallclock_mean(i)=wallclock_mean(i)/nproc
      cpu_mean(i)=cpu_total(i)/nproc
! To stop a divide by zero later:
      IF (wallclock_mean(i)  ==  0.0) wallclock_mean(i)=1.0e-20
      IF (cpu_mean(i)  ==  0.0) cpu_mean(i)=1.0e-20
! and now the standard deviation
      wallclock_sd(i)=0.0
      cpu_sd(i)=0.0
      DO j=0,nproc-1
        wallclock_sd(i)=wallclock_sd(i)+                                       &
          (wallclock_times(j)-wallclock_mean(i))*                              &
          (wallclock_times(j)-wallclock_mean(i))
        cpu_sd(i)=cpu_sd(i)+(cpu_times(j)-cpu_mean(i))*                        &
                      (cpu_times(j)-cpu_mean(i))
      END DO
      wallclock_sd(i)=SQRT(wallclock_sd(i)/nproc)
      cpu_sd(i)=SQRT(cpu_sd(i)/nproc)

! Calculate the median
      DO j=0,nproc-2
        DO k=j+1,nproc-1
          IF (wallclock_times(k)  >   wallclock_times(j)) THEN
            sortwork_real=wallclock_times(j)
            wallclock_times(j)=wallclock_times(k)
            wallclock_times(k)=sortwork_real
          END IF
          IF (cpu_times(k)  >   cpu_times(j)) THEN
            sortwork_real=cpu_times(j)
            cpu_times(j)=cpu_times(k)
            cpu_times(k)=sortwork_real
          END IF
        END DO
      END DO

      IF (MOD(nproc,2)  ==  0) THEN
        wallclock_median(i)=(wallclock_times((nproc/2)-1)+                     &
                             wallclock_times(nproc/2))*0.5
        cpu_median(i)=(cpu_times((nproc/2)-1)+                                 &
                       cpu_times(nproc/2))*0.5
      ELSE
        wallclock_median(i)=wallclock_times(nproc/2)
        cpu_median(i)=cpu_times(nproc/2)
      END IF

    END IF ! am I PE 0?

  END DO ! loop over sections

! Sort and output the information on PE 0

  IF (mype  ==  0) THEN

    DO i=1,summ_n_timers-1
      DO j=(i+1),summ_n_timers
        IF (wallclock_max(j)  >   wallclock_max(i)) THEN

! Swap the entries I and J

        sortwork_char=summ_section(i)
        summ_section(i)=summ_section(j)
        summ_section(j)=sortwork_char

        sortwork_real=wallclock_mean(i)
        wallclock_mean(i)=wallclock_mean(j)
        wallclock_mean(j)=sortwork_real

        sortwork_real=wallclock_median(i)
        wallclock_median(i)=wallclock_median(j)
        wallclock_median(j)=sortwork_real

        sortwork_real=wallclock_sd(i)
        wallclock_sd(i)=wallclock_sd(j)
        wallclock_sd(j)=sortwork_real

        sortwork_real=wallclock_max(i)
        wallclock_max(i)=wallclock_max(j)
        wallclock_max(j)=sortwork_real

        sortwork_real=wallclock_min(i)
        wallclock_min(i)=wallclock_min(j)
        wallclock_min(j)=sortwork_real

        sortwork_int=wallclock_min_pe(i)
        wallclock_min_pe(i)=wallclock_min_pe(j)
        wallclock_min_pe(j)=sortwork_int

        sortwork_int=wallclock_max_pe(i)
        wallclock_max_pe(i)=wallclock_max_pe(j)
        wallclock_max_pe(j)=sortwork_int

        sortwork_real=cpu_mean(i)
        cpu_mean(i)=cpu_mean(j)
        cpu_mean(j)=sortwork_real

        sortwork_real=cpu_median(i)
        cpu_median(i)=cpu_median(j)
        cpu_median(j)=sortwork_real

        sortwork_real=cpu_sd(i)
        cpu_sd(i)=cpu_sd(j)
        cpu_sd(j)=sortwork_real

        sortwork_real=cpu_max(i)
        cpu_max(i)=cpu_max(j)
        cpu_max(j)=sortwork_real

        sortwork_real=cpu_min(i)
        cpu_min(i)=cpu_min(j)
        cpu_min(j)=sortwork_real

        sortwork_int=cpu_min_pe(i)
        cpu_min_pe(i)=cpu_min_pe(j)
        cpu_min_pe(j)=sortwork_int

        sortwork_int=cpu_max_pe(i)
        cpu_max_pe(i)=cpu_max_pe(j)
        cpu_max_pe(j)=sortwork_int

        END IF
      END DO
    END DO

! and write out the information
    WRITE(6,*)
    IF (timer_kind  ==  1) THEN
      WRITE(6,*) 'MPP : Non Inclusive timer summary'
    ELSE
      WRITE(6,*) 'MPP : Inclusive timer summary'
    END IF

    WRITE(6,*)
    WRITE(6,*)  'WALLCLOCK  TIMES'

    WRITE(6,'(4X,a,19X,a,3X,a,7X,a,3X,a,6X,a,3X,a,8X,a,3X,a)')                 &
      'ROUTINE','MEAN','MEDIAN','SD','% of mean','MAX','(PE)','MIN','(PE)'

    DO i=1,summ_n_timers

      WRITE(6,                                                                 &
      '(I3,1X,A20,1X,3(1X,F8.2),5X,F6.2,a,2(1X,F8.2,1X,a,I6,a) )')             &
                 i,summ_section(i),                                            &
                 wallclock_mean(i),wallclock_median(i),                        &
                 wallclock_sd(i),                                              &
                 (wallclock_sd(i)/wallclock_mean(i))*100.0,'%',                &
                 wallclock_max(i),'(',wallclock_max_pe(i),')',                 &
                 wallclock_min(i),'(',wallclock_min_pe(i),')'
    END DO

    WRITE(6,*)
    WRITE(6,*)  'CPU TIMES (sorted by wallclock times)'

    WRITE(6,'(4X,a,19X,a,3X,a,7X,a,3X,a,6X,a,3X,a,8X,a,3X,a)')                 &
     'ROUTINE','MEAN','MEDIAN','SD','% of mean','MAX','(PE)','MIN','(PE)' 

    DO i=1,summ_n_timers
      WRITE(6,                                                                 &
      '(I3,1X,A20,1X,3(1X,F8.2),5X,F6.2,a,2(1X,F8.2,1X,a,I6,a) )')             &
                 i,summ_section(i),                                            &
                 cpu_mean(i),cpu_median(i),                                    &
                 cpu_sd(i),                                                    &
                 (cpu_sd(i)/cpu_mean(i))*100.0,'%',                            &
                 cpu_max(i),'(',cpu_max_pe(i),')',                             &
                 cpu_min(i),'(',cpu_min_pe(i),')'
    END DO

  END IF
  WRITE(6,*)



END DO ! loop over timer kind

! Finally copy the timer info into the last_* arrays so that the
! intermediate timer can calculate the timings since this point

DO i=1,in_number_of_timers
  last_in_number_of_times_timed(i)=in_number_of_times_timed(i)
  last_in_cpu_time_elapsed(i)=in_cpu_time_elapsed(i)
  last_in_wallclock_time_elapsed(i)=                                           &
    in_wallclock_time_elapsed(i)
END DO

DO i=1,ni_number_of_timers
  last_ni_number_of_times_timed(i)=ni_number_of_times_timed(i)
  last_ni_cpu_time_elapsed(i)=ni_cpu_time_elapsed(i)
  last_ni_wallclock_time_elapsed(i)=                                           &
    ni_wallclock_time_elapsed(i)
END DO


 9999 CONTINUE
END IF ! L_IO_Server

DEALLOCATE(cpu_times)      
DEALLOCATE(wallclock_times)

IF (lhook) CALL dr_hook('TIMER_OUTPUT',zhook_out,zhook_handle)
RETURN

END SUBROUTINE timer_output
