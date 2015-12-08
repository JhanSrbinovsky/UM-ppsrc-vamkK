! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Control variables

!  Description:  

! Migrated from include file cntlall.h

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90

MODULE nlstcall_mod

  USE chsunits_mod, ONLY : nunits
  USE missing_data_mod, ONLY : imdi, rmdi


  IMPLICIT NONE

  ! Array holding original data time (prior to assimilation)
  INTEGER:: model_basis_time(6) = imdi

  ! Model analysis time in hours since Basis Time
  ! UM6.5 - Replace model_analysis_hrs by model_analysis_mins 
  !         model_analysis_hrs changed to REAL
  REAL :: model_analysis_hrs     = rmdi
  INTEGER :: model_analysis_mins = imdi

  INTEGER:: ncpu              = imdi  ! No of CPUs assigned to the program
  INTEGER:: ancil_reftime(6)  = imdi  ! Ref. time for updating ancillaries
  INTEGER:: run_target_end(6) = imdi  ! Target end time for this run

  INTEGER:: Num_ALBCs            = imdi ! Number of atmos boundary files
  INTEGER:: ALBC2_StartTime_mins = imdi ! VT of first block of data in 2nd
                                        ! atmos boundary file, in minutes
                                        ! from start of run
  
  ! Increment to be added on each resubmission of the job.
  INTEGER:: run_resubmit_inc(6) = imdi

  ! Number of field headers reserved for non-mean PPfiles on each
  ! unit
  INTEGER:: pp_len2_look(20:NUNITS) = imdi

  ! Internally defined PP packing code
  INTEGER:: pp_pack_code(20:NUNITS) = imdi

  ! Frequency of initialisation of FTunit
  INTEGER:: ft_steps(20:NUNITS)      = imdi
  INTEGER :: ft_firststep(20:NUNITS) = imdi    ! ... starting at step number .
  INTEGER :: ft_laststep(20:NUNITS)  = imdi    ! ... ending at step number ..
  LOGICAL:: latmosnext               = .FALSE.
  LOGICAL :: loceannext              = .FALSE.
  LOGICAL:: LPP                      = .FALSE. ! Activate PPCTL
  LOGICAL:: lpp_select(20:NUNITS)    = .FALSE. ! Activate PP init on unit
  LOGICAL:: ldump                    = .FALSE. ! Activate DUMPCTL
  LOGICAL:: lmean                    = .FALSE. ! Activate MEANCTL
  LOGICAL:: lhistory                 = .FALSE. ! Update TEMP history file
  LOGICAL:: lprint                   = .FALSE. ! Activate PRINTCTL
  LOGICAL:: linterface               = .FALSE. ! Activate GEN_INTF
  LOGICAL:: lexit                    = .FALSE. ! Activate EXITCHEK
  LOGICAL:: ljobrelease              = .FALSE. ! Activate JOBCTL

  ! Select printed diags from means
  LOGICAL:: lmeanPR(4)               = .FALSE.

  LOGICAL:: lancillary               = .FALSE. ! Activate UP_ANCIL
  LOGICAL:: lboundary                = .FALSE. ! Activate UP_BOUND
  LOGICAL:: lassimilation            = .FALSE. ! Activate assimilation

  !Activate detailed TIMER routine
  LOGICAL:: ltimer                   = .FALSE.
  
  ! T : D1 copied to memory for AO coupling
  LOGICAL:: l_ao_d1_memory           = .FALSE.

  ! Real-period climate means
  LOGICAL:: lclimrealyr             = .FALSE.


  CHARACTER(LEN=4) :: expt_id           ! Unique alphanumeric serial number
!                                       ! associated with model
!                                       ! (Non-Operational expts)
!                                       !
!                                       ! Operational run name
!                                       ! (Operational expts)
  CHARACTER(LEN=1) ::  job_id           ! Unique alphanumeric job identifier
!                                       ! used for networking
  CHARACTER(LEN=14) :: model_status     ! Operational or NonOperational
  CHARACTER(LEN=14) :: model_assim_mode ! Atmosphere,Ocean,Coupled or None
  CHARACTER(LEN=17) :: time_convention  ! Relative, Timestep, Absolute_long,
!                                         Absolute_standard or Absolute_short
!
  CHARACTER(LEN=1) :: type_letter_1(20:NUNITS) ! File type letter #1
  CHARACTER(LEN=1) :: type_letter_2(20:NUNITS) ! File type letter #2
  CHARACTER(LEN=1) :: type_letter_3(20:NUNITS) ! File type letter #3
!
  CHARACTER(LEN=1) :: ft_output(20:NUNITS)  ! "Y" if output file on unit
  CHARACTER(LEN=1) :: ft_select(20:NUNITS)  ! "Y" if file selected for post
!                                             processing request.
  CHARACTER(LEN=1) :: ft_archsel(20:NUNITS) ! "Y" if file to be archived.
!
  CHARACTER(LEN=10) ::run_assim_mode      ! cf model_assim_mode (Oper use)
  CHARACTER(LEN=1) :: control_resubmit    ! User flag for auto resubmit
  
  NAMELIST / NLSTCALL /                                              &
     model_basis_time, model_analysis_mins,                          &
     ncpu, ancil_reftime, run_target_end,                            &
     Num_ALBCs, ALBC2_StartTime_mins,                                &
     run_resubmit_inc, pp_len2_look, pp_pack_code,                   &
     ft_steps, ft_firststep, ft_laststep,                            &
     latmosnext, loceannext, LPP, lpp_select, ldump, lmean,          &
     lhistory, lprint, linterface, lexit, ljobrelease,               &
     lmeanPR, lancillary, lboundary, lassimilation,                  &
     ltimer, l_ao_d1_memory,                                         &
     lclimrealyr,                                                    &
     expt_id, job_id,                                                &
     model_status, model_assim_mode,                                 &
     time_convention,                                                &
     type_letter_2,                                                  &
     ft_output, ft_select, ft_archsel,                               &
     run_assim_mode, control_resubmit

END MODULE nlstcall_mod
