! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Control variables

!  Description:  Control variables for generic aspects of internal models.
!                Generic means values likely to be common to the control
!                of any sub-model/internal model.

! Migrated from include file cntlgen.h

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90

MODULE nlstgen_mod

  USE Submodel_Mod, ONLY: N_Internal_Model_Max
  USE missing_data_mod, ONLY:  rmdi, imdi
  IMPLICIT NONE
      
! Max no. of irregular times for dumps
  INTEGER, PARAMETER :: Dumptimes_Len1 = 40

! No. of time intervals for climate meaning
  INTEGER, PARAMETER :: Meanfreq_Len1 = 4

! Max no. of irregular times for job release
  INTEGER, PARAMETER :: Jobrel_Len1 = 10

  INTEGER       :: steps_per_periodim(N_Internal_Model_Max) = imdi
  INTEGER       :: secs_per_periodim(N_Internal_Model_Max)  = imdi

! Number of steps between atmosphere restart dumps
  INTEGER       :: dumpfreqim(N_Internal_Model_Max) = imdi

! Archiving frequency  for atmos dumps
  INTEGER       :: archdump_freqim(N_Internal_Model_Max) = imdi

! Timesteps (from start of run) at which restart dumps are written
  INTEGER       :: dumptimesim(Dumptimes_Len1,N_Internal_Model_Max) = imdi

! Indicators for mean dump frequency
  INTEGER       :: meanfreqim(Meanfreq_Len1,N_Internal_Model_Max) = imdi

! PP field selectors
  INTEGER       :: ppselectim(Meanfreq_Len1,N_Internal_Model_Max) = imdi

! Switches for pp field archive
  INTEGER       :: archppselim(Meanfreq_Len1,N_Internal_Model_Max) = imdi

! Number of field headers to reserve for internal model mean  PPfiles
  INTEGER       :: pp_len2_meanim(Meanfreq_Len1,N_Internal_Model_Max) = 30000

! Reference time for production of means
  INTEGER       :: mean_reftimeim(6,N_Internal_Model_Max) = imdi

! Step numbers  at which to release user-specified scripts
  INTEGER       :: jobrel_stepim(Jobrel_Len1,N_Internal_Model_Max) = imdi
  
! Offset in timesteps to release user-specified scripts
  INTEGER       :: jobrel_offsetim = 0
  
! Offset for dump archiving
  INTEGER       :: archdump_offsetim(N_Internal_Model_Max) = imdi

! Unit reserved for mean PPs
  INTEGER       :: ft_meanim(N_Internal_Model_Max) = 27

! Packing indicator for dumps
  INTEGER       :: dump_packim(N_Internal_Model_Max) = imdi

  LOGICAL   :: llboutim(N_Internal_Model_Max) = .FALSE.  ! Lateral b.c.'s

  NAMELIST / NLSTCGEN /                                                 &
        steps_per_periodim, secs_per_periodim, dumpfreqim,              &
        archdump_freqim, dumptimesim, ppselectim,                       &
        archppselim,  meanfreqim, mean_reftimeim,                       &
        jobrel_stepim, archdump_offsetim,                               &
        dump_packim,                                                    &
        llboutim, jobrel_offsetim


END MODULE nlstgen_mod 
