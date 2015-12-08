! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation


MODULE comobs_mod

USE atmos_max_sizes, ONLY: model_levels_max
IMPLICIT NONE

! This module contains parameters relating to observations in
!   atmospheric AC assimilation 

!    The AC Observation files are read in at the start of an
!    assimilation (in routine RDOBS2). The observations are then
!    merged and any observations not required by the assimilation
!    to follow will be compressed out. The required observations are
!    then written to a temporary file attached to Unit 15.
!
!    On each timestep the observations are read in from Unit 15 and
!    assimilated. This is done until the AC Observation files
!    need to be read in again which will normally be every 6 (3) hours
!    for a global (elf) assimilation. A new set of observations will
!    then be written to unit 15.
!
!    COMOBS contains informtion on the observations in store.
!
!    Observations are assigned to types. Each type provides one sort
!    of information and has its own detailed format. Note that one
!    observation report may be in several types, for example, a sonde
!    observation may report temperature, wind and Relative Humidity and
!    the data will come under types 201, 301 and 401 respectively.
!
!    There are NOBTYP types which are listed in OBSTYP. The
!    The observation types are listed in decreasing order of the
!    number of data values in NDATAV. Each type OBSTYP(J) has :-
!
!        1. NOBS(J) observations stored for it.
!        2. NDATAV(J) data values in each obs.
!        3. A total of NDATAV(J)*NOBS(J) data values.
!
!    The array MDISPOBT gives the displacement to the first
!    observation of each type.
!
!    ie. MDISPOBT(1)=0 and MDISPOBT(J+1)=MDISPOBT(J)+NOBS(J)
!
!    The Total number of observations and data values after merging
!    the observation files are TNOBS and TNDV respectively.
!
!    Maximum allowed values of NOBTYP is NOBTYPMX
!    Maximum allowed values of NDATAV is NDATAVMX
!    Maximum allowed values of TNDV   is TNDVMAX
!
!  MAXNLEV1    - Maximum no of levels for any obs type in OBSTYP + 1
!  MDISPOBT    - Offset to first obs for each type in OBS.
!  MISSD       - Value used for missing data in obs files (=5555.0)
!  NDATAV      - No of data values for each type in OBSTYP.
!  NDVHDR      - First NDVHDR data values common to all observations.
!  NERLEV1     - No of data value corresponding to first error ratio
!              - for each type in OBSTYP.
!  NOBLEV      - No of levels for each type in OBSTYP.
!  NOBS        - No of obs    for each type in OBSTYP.
!  NOBTYP      - No of observation types in OBSTYP.
!  OBLAYERB    - Pressure Levels of layer boundaries for each type.
!              - Set to MISSD for obs types not concerned.
!  OBLEVELS    - Obs Pressure Levels for each type.
!              - Set to MISSD for obs types not concerned.
!  OBLEVTYP    - Level Type obs data in on for each type in OBSTYP.
!  OBSTYP      - List of observation types in observation files.
!  OBS_LAT_N   - N ) Boundaries of area for
!  OBS_LAT_S   - S ) which obs are to fetched for.
!  OBS_LONG_W  - W )
!  OBS_LONG_E  - E )
!  OBS_REF_YY  - Year  ) Reference
!  OBS_REF_MM  - Month ) Time/Date for observations
!  OBS_REF_DD  - Day   ) = Start of assimilation.
!  OBS_REF_HH  - Hour  )
!  OBS_REF_MIN - Mins  )
!  OBS_NO_ST   - Offset to first observation number for each type
!  TIMEINT     - Interval in minutes between reading obs files.
!  TIMENEXT    - Time in minutes relative to start of assimilation
!              - when obs files are next to be read.
! ---------------------------------------------------------------------

! Note that some variables are stored as an instance of a sequenced type.
! This allows communication via MPI in a fairly safe manner.
! Any changes to the sequenced type will need to ensure that INTEGERS are
! seperate from REALs and that the obs_info_int_len and obs_info_real_len
! variables are updated appropriately.

INTEGER, PARAMETER :: ndatavmx           = 6+3*model_levels_max
INTEGER, PARAMETER :: noblevmx           = model_levels_max+1
INTEGER, PARAMETER :: max_num_acob_files = 100
INTEGER, PARAMETER :: num_ob_file_types  = 11
INTEGER, PARAMETER :: nobtypmx           = 182


TYPE obs_info_type
  SEQUENCE
  INTEGER :: maxnlev1
  INTEGER :: nobtyp
  INTEGER :: ndvhdr
  INTEGER :: obs_ref_yy
  INTEGER :: obs_ref_mm
  INTEGER :: obs_ref_dd
  INTEGER :: obs_ref_hh
  INTEGER :: obs_ref_min
  INTEGER :: num_used_files
  INTEGER :: obstyp(nobtypmx)
  INTEGER :: noblev(nobtypmx)
  INTEGER :: ndatav(nobtypmx)
  INTEGER :: nerlev1(nobtypmx)
  INTEGER :: nobs(nobtypmx)
  INTEGER :: oblevtyp(nobtypmx)
  INTEGER :: mdispobt(nobtypmx)
  INTEGER :: obs_no_st(nobtypmx)
  INTEGER :: per_file_tndvmax(max_num_acob_files)
  INTEGER :: filename_len(max_num_acob_files)


  REAL    :: missd
  REAL    :: timeint
  REAL    :: timenext
  REAL    :: obs_lat_n
  REAL    :: obs_lat_s
  REAL    :: obs_long_e
  REAL    :: obs_long_w
  REAL    :: oblevels(noblevmx,nobtypmx)
  REAL    :: oblayerb(noblevmx+1,nobtypmx)

END TYPE obs_info_type

TYPE (obs_info_type) :: obs_info

! Amount of integer data in type - for MPI comms
INTEGER, PARAMETER  :: obs_info_int_len = (8 * nobtypmx) +         &
                                          (2 * max_num_acob_files) &
                                             + 9 

! Amount of real data in type - for MPI comms
INTEGER, PARAMETER  :: obs_info_real_len = (2 * noblevmx + 1) * nobtypmx + 7



CHARACTER (LEN=256) :: used_files(max_num_acob_files)
CHARACTER (LEN=30)  :: ob_file_type(num_ob_file_types) =                &
                       (/'Surface                       ',              &
                         'Sonde                         ',              &
                         'Aircraft                      ',              &
                         'Sat120                        ',              &
                         'Sat500                        ',              &
                         'GLOSS                         ',              &
                         'Satwind                       ',              &
                         'Scatwind                      ',              &
                         'MOPS                          ',              &
                         'Precip                        ',              &
                         'Test                          '/)


END MODULE comobs_mod
