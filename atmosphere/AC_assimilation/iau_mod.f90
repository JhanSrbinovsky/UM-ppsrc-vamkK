! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parameters and variables for Incremental Analysis Update (IAU) scheme

MODULE IAU_mod

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE um_types
USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,        &
                       len_realhd, itgrid, iugrid,                &
                       max_num_fc_times, max_num_clim_times,      &
                       max_num_in_flux, len2_lookuppreferred,     &
                       len2_lookupprevious, len2_lookupclimate

USE missing_data_mod, ONLY: imdi, rmdi

IMPLICIT NONE

! Common blocks:

! For Len_FixHd:

!-------------------------------------------------------------------------------
! [1]: Global constants.
!-------------------------------------------------------------------------------

INTEGER, PARAMETER :: MaxFileNameLen = 255 ! (As used in GEN)

INTEGER, PARAMETER :: MaxNum_IAU_incs = 99 ! Maximum number of increment files

INTEGER, PARAMETER :: IAU_unit = 108 ! Unit number for IAU increments

! Data for fields that may be read from the increment files:
INTEGER, PARAMETER :: IAU_NumFldCodes = 22
INTEGER, PARAMETER :: IAU_FldCodes(IAU_NumFldCodes) =       &
                        (/  2       ,  3       ,  4       , &
                            9       ,  10      ,  12      , &
                            24      ,  90      ,  150     , &
                            253     ,  254     ,  255     , &
                            407     ,  431     ,  432     , &
                            433     ,  434     ,  435     , &
                            436     ,  480     ,  16207   , &
                            18001    /)
CHARACTER(LEN=7), PARAMETER :: &
                      IAU_FldDescs(IAU_NumFldCodes) =       &
                        (/ 'u      ', 'v      ', 'theta  ', &
                           'soilm  ', 'q      ', 'qCF    ', &
                           'sfctemp', 'aerosol', 'w      ', &
                           'rho    ', 'qCL    ', 'exner  ', &
                           'p      ', 'dust1  ', 'dust2  ', &
                           'dust3  ', 'dust4  ', 'dust5  ', &
                           'dust6  ', 'ozone  ', 'qT     ', &
                           'qT     ' /)

! TimeID choices, copied from VAR:
INTEGER, PARAMETER :: TimeID_Once     = 0 ! Instantaneous increment
INTEGER, PARAMETER :: TimeID_Constant = 3 ! Constant error tendency

! Call frequency codes:
INTEGER, PARAMETER :: CallFreq_Never              = 0
INTEGER, PARAMETER :: CallFreq_EveryCallIncsToAdd = 1
INTEGER, PARAMETER :: CallFreq_EndInc1Window      = 2
INTEGER, PARAMETER :: CallFreq_LastIAUCall        = 3

!-------------------------------------------------------------------------------
! [2]: Global type definitions.
!-------------------------------------------------------------------------------

! Structure to hold IAU increment data:
TYPE IAU_inc_type

  CHARACTER(MaxFileNameLen) :: FilePath 

  LOGICAL :: FieldsToAdd    ! Any fields to be added to the model?

  LOGICAL :: TendencyIncs   ! Are increments (per second) tendencies?

  LOGICAL :: Contains_q_qCL_or_qCF ! Are there any increments to q, qCL or qCF?
  LOGICAL :: Contains_qT           ! Are there any increments to qT?
  LOGICAL :: Contains_sfctemp      ! Are there any increments to sfctemp?
  LOGICAL :: Contains_soilm        ! Are there any increments to soilm?

  INTEGER       :: InsertionStartTS ! Start timestep for insertion
  INTEGER       :: InsertionEndTS   ! End   timestep for insertion
  REAL, POINTER :: Weights(:)       ! Insertion weights

  INTEGER          :: FixHd(Len_FixHd)
  INTEGER          :: Len1Lookup
  INTEGER          :: Len2Lookup
  INTEGER, POINTER :: Lookup(:,:)

  LOGICAL :: StoreFields  ! Store the increment fields between IAU calls?
  LOGICAL :: StorePacked  ! Store them 32-bit-packed?
  LOGICAL :: FieldsStored ! Are the fields currently stored?

  INTEGER                    :: flds_len     ! Required length of storage array
  REAL,              POINTER :: flds     (:) ! Storage array for increments
  REAL(KIND=real32), POINTER :: flds32bit(:) ! 32-bit alternative

END TYPE IAU_inc_type

!-------------------------------------------------------------------------------
! [3]: Global variables.
!-------------------------------------------------------------------------------

TYPE (IAU_inc_type), ALLOCATABLE :: IAU_incs(:) ! IAU increments

INTEGER :: IAU_LocFldLens(IAU_NumFldCodes) ! Local lengths of fields that may
                                           ! be read from the increment files

REAL :: q_min ! Minimum value allowed for q after addition of q increments

INTEGER :: IAU_FirstCallTS = -1 ! Timestep number for first call to IAU
INTEGER :: IAU_LastCallTS  = -1 ! Timestep number for last  call to IAU

!-------------------------------------------------------------------------------
! [3.1]: Namelist variables.
!-------------------------------------------------------------------------------

LOGICAL :: L_IAU = .FALSE.  ! Activate IAU scheme?

INTEGER :: Num_IAU_incs = imdi ! Number of increment files

LOGICAL :: L_IAU_SpecifyInc1Filter = .FALSE. ! Specify filter details for first
                                            ! increment? Was TRUE

! Filter details for first increment:
INTEGER :: IAU_FilterType    = imdi   ! Filter type 1- uniform, 2- triangular,
                                      !             3- Lanczos, 4- Dolph
INTEGER :: IAU_StartMin      = imdi   ! Start minute of filtering period
                                      ! Was 0
INTEGER :: IAU_EndMin        = imdi   ! End   minute of filtering period
                                      ! Was 0
INTEGER :: IAU_ApexMin       = imdi   ! Apex minute for triangular filter
                                      ! Was 0
REAL    :: IAU_Cutoff_period = rmdi   ! Filter cutoff period in hours - was 0.0
REAL    :: IAU_SBE_period    = rmdi   ! Stop band edge period in hours - was 0.0
      
! Switches:
LOGICAL :: L_IAU_IncDiags       = .FALSE. ! Write increment diagnostics?
LOGICAL :: L_IAU_CalcExnerIncs  = .FALSE. ! Use p increments to calculate
                                          ! exner increments?
LOGICAL :: L_IAU_CalcThetaIncs  = .FALSE. ! Calculate theta increments using
                                          ! exner and q increments?
LOGICAL :: L_IAU_CalcRhoIncs    = .FALSE. ! Calculate rho increments using
                                          ! exner, theta and (possibly) q
                                          ! increments?
LOGICAL :: L_IAU_IncTStar       = .FALSE. ! Add level-one temperature
                                          ! increments to surface temperature
                                          ! and top-level soil temperature?
LOGICAL :: L_IAU_IncTStar_tile  = .FALSE. ! Add level-one temperature
                                          ! increments to surface temperature
                                          ! on land-surface tiles
                                          ! Was TRUE
LOGICAL :: L_IAU_IncTSurf_Snow  = .FALSE. ! Add level-one temperature
                                          ! increments to snow temperature
LOGICAL :: L_IAU_IncTLake       = .FALSE. ! Add level-one temperature
                                          ! increments to lake temperature
LOGICAL :: L_IAU_UseSfctempPerts= .FALSE. ! Include surface temperature
                                          ! perturbations and add to TStar
                                          ! and top-level soil temperature?
LOGICAL :: L_IAU_UseSoilmPerts  = .FALSE. ! Add soil moisture perturbations
                                          ! to soil-levels?
LOGICAL :: L_IAU_ResetPoles     = .FALSE. ! Reset polar rows of relevant
                                          ! increment fields to their mean
                                          ! values? Was TRUE
LOGICAL :: L_IAU_RemoveSS       = .TRUE.  ! Remove supersaturation wrt water?
                                          ! Not in name list
LOGICAL :: L_IAU_DumpTS0State   = .FALSE. ! Write out model state
                                          ! immediately after timestep-zero
                                          ! call to IAU scheme
LOGICAL :: L_IAU_ApplyQTCorrections = .FALSE.     ! Apply corrections for
                                                  ! processing of qT incs.
LOGICAL :: L_IAU_Add_qT_prime_to_q = .FALSE. ! Replace q with qT_plus,
                                                  ! bypassing Var_DiagCloud.
LOGICAL :: L_IAU_IncrementIce   = .FALSE. ! Diagnose ice cloud increments as
                                          ! well as humidity and liquid cloud
                                          ! increments from qT increments?
LOGICAL :: L_IAU_ScaleCloud     = .FALSE. ! If diagnosing humidity and cloud
                                          ! increments from qT increments, 
                                          ! scale increments to be within
                                          ! physical bounds?
LOGICAL :: L_IAU_LimitUpperThetaIncs = .FALSE. ! Apply limits to size of 
                                               ! upper-level theta increments?
LOGICAL :: L_IAU_SetOzoneMin    = .FALSE. ! Reset ozone to oz_min if ozone was
                                          ! negative?

! Parameters for use with L_IAU_LimitUpperThetaIncs:
REAL :: IAU_LimitUpperThetaIncs_pBound = rmdi ! Pressure boundary - was 200.0
REAL :: IAU_LimitUpperThetaIncs_maxInc = rmdi ! Maximum absolute value of
                                               ! increment after multiplication
                                               ! by IAU weight - was 100.0

! Parameters and switches for use with QLimits:
INTEGER :: IAU_QLimitsCallFreq  = imdi    ! Call frequency 
                                          !  - was CallFreq_EndInc1Window
LOGICAL :: L_IAU_QLimitsDiags   = .FALSE. ! Write diagnostics?
LOGICAL :: L_IAU_RmNonTropQIncs = .FALSE. ! Remove non-trop q increments?
REAL    :: IAU_trop_min_RH      = rmdi    ! Lower limit to apply to trop RH
                                          ! Was 0.0
REAL    :: IAU_nonTrop_max_q    = rmdi    ! Upper limit to apply to non-trop q
                                          ! Was 3.0E-06
REAL    :: IAU_nonTrop_min_q    = rmdi    ! Lower limit to apply to non-trop q
                                          ! Was 1.0E-06
REAL    :: IAU_nonTrop_max_RH   = rmdi    ! Upper limit to apply to non-trop RH
                                          ! Was 0.1
REAL    :: IAU_trop_min_p       = rmdi    ! Minimum tropospheric pressure
                                          ! Was 1.0E+04
REAL    :: IAU_trop_max_PV      = rmdi    ! Maximum tropospheric ABS(PV)
                                          ! Was 2.5E-06
REAL    :: IAU_nonTrop_max_p    = rmdi    ! Maximum non-trop     pressure
                                          ! Was 4.0E+04

! Want to make this switch obsolete ASAP:
LOGICAL :: L_IAU_IgnoreTopLevThetaIncs = .FALSE.

! Modifications to limits applied when L_IAU_RemoveSS is true:
REAL    :: IAU_qsatScale        = 1.0   ! Fraction of qsat to be used as upper
                                        ! bound for q - not in name list
INTEGER :: IAU_qsatScale_maxLev = 0     ! Maximum level at which to apply
                                        ! IAU_qsatScale - not in name list 

! Parameters to modify output from Var_DiagCloud in IAU routine:
LOGICAL :: L_IAU_LimitIncOp   = .FALSE. ! Limit qCL increments generated by
                                        ! Var_DiagCloud?
REAL    :: IAU_qCLThreshScale = 1.0     ! Scaling factor for limit to qCL
                                        ! increments - not in name list 

! Parameters and switches for use within Var_DiagCloud itself.
! (Currently Var_DiagCloud is called only from IAU, so it makes sense to
! control them via the IAU namelist.)

REAL :: DiagCloud_Tol_q = rmdi  ! Was 0.0
  ! Tolerance for qcl==0 and qcf==0 tests in VarDiag_Cloud.
  ! (Tests suggest that a value of 1.0e-14 is more-or-less optimal.)

REAL :: DiagCloud_Tol_FM = rmdi  ! Was 0.0
  ! Tolerance for (1-FM)==0  tests in VarDiag_Cloud.
  ! (Tests suggest that a value of 1.0e-28 is more-or-less optimal.)

LOGICAL :: DiagCloud_ApplyCompLimits = .FALSE. ! Was TRUE
  ! If true, take measures to avoid machine precision issues in Var_DiagCloud.
  ! (Once we're confident that this option has no nasty side-effects, we should
  ! remove the switch and always use it.)   

INTEGER :: DiagCloud_NumMaxLoops = imdi  ! Was 150
  ! Maximum number of loops for convergence in Var_DiagCloud.

REAL    :: DiagCloud_QN_CompRegimeLimit = rmdi ! Was 20.0  
  ! When the magnitude of QN gets much beyond this parameter, the standard
  ! calculations for qCL_calc in Var_DiagCloud give deviations from its upper
  ! or lower bound that are highly sensitive to rounding. In the case where QN
  ! is negative, this makes the iterative algorithm used to calculate cloud
  ! quantities extremely sensitive to small changes in the input data, and when
  ! the IncrementIce option is used often leads to non-convergence. To get
  ! around this, when ApplyCompBounds is set to true, and the magnitude of QN
  ! exceeds the above limit, we set qCL_calc directly to the appropriate
  ! bounding value.
  !
  ! QN < -QN_CompRegimeLimit occurs at points very far from any saturation.
  ! QN >  QN_CompRegimeLimit occurs at points very far from any unsaturation.
  !
  ! The parameter value was determined experimentally by determining when the
  ! path through the code starts becoming affected by least-significant-bit
  ! changes to the LS states. The current value is appropriate for 64-bit
  ! reals.

!-------------------------------------------------------------------------------
! [4]: IAU namelist.
!-------------------------------------------------------------------------------

NAMELIST / IAU_nl /             &
L_IAU,                          &
Num_IAU_incs,                   &
L_IAU_SpecifyInc1Filter,        &
IAU_FilterType,                 &
IAU_StartMin,                   &
IAU_EndMin,                     &
IAU_ApexMin,                    &
IAU_Cutoff_period,              &
IAU_SBE_period,                 &
L_IAU_IncDiags,                 &
L_IAU_CalcExnerIncs,            &
L_IAU_CalcThetaIncs,            &
L_IAU_CalcRhoIncs,              &
L_IAU_IncTStar,                 &
L_IAU_IncTStar_tile,            &
L_IAU_IncTSurf_Snow,            &
L_IAU_UseSoilmPerts,            &
L_IAU_UseSfctempPerts,          &
L_IAU_ResetPoles,               &
L_IAU_DumpTS0State,             &
L_IAU_ApplyQTCorrections,       &
L_IAU_Add_qT_prime_to_q,        &
L_IAU_IncrementIce,             &
L_IAU_ScaleCloud,               &
L_IAU_LimitUpperThetaIncs,      &
L_IAU_SetOzoneMin,              &
IAU_LimitUpperThetaIncs_pBound, &
IAU_LimitUpperThetaIncs_maxInc, &
IAU_QLimitsCallFreq,            &
L_IAU_QLimitsDiags,             &
L_IAU_RmNonTropQIncs,           &
IAU_trop_min_RH,                &
IAU_nonTrop_max_q,              &
IAU_nonTrop_min_q,              &
IAU_nonTrop_max_RH,             &
IAU_trop_min_p,                 &
IAU_trop_max_PV,                &
IAU_nonTrop_max_p,              &
L_IAU_IgnoreTopLevThetaIncs,    &
DiagCloud_Tol_q,                &
DiagCloud_Tol_FM,               &
DiagCloud_ApplyCompLimits,      &
DiagCloud_NumMaxLoops,          & ! must be greater than 50
DiagCloud_QN_CompRegimeLimit

END MODULE IAU_mod
