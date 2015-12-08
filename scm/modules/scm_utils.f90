!
! MODULE SCM_UTILS--------------------------------------------------------------
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! ******************************COPYRIGHT***************************************


MODULE scm_utils

!-------------------------------------------------------------------------------
! Description:
!   SCM UTILITY MODULE; For use with SCM in routines with "USE scm_utils"
!
!   Contains subroutines:
!   scm_message: To send messages to stdout or unit number with timestep
!                information. Calls with no arguments will output the current
!                model timestep.
!
!   scm_trap_nan:
!   alloc_common_scm:
!   dealloc_common_scm:
!-------------------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
  USE um_types
  USE missing_data_mod
  USE nstypes, ONLY: ntype, npft, nnvg, soil
  USE parkind1, ONLY: jpim, jprb
  USE yomhook,  ONLY: lhook, dr_hook

  USE soil_param, ONLY: dzsoil
  IMPLICIT NONE

! Missing data indicators

  ! Precision of integers passed in/out out netcdf routines
  INTEGER, PARAMETER, PRIVATE :: incdf = integer32
  INTEGER           , PRIVATE :: nan_list_count = 0 ! Updated by scm trap nan
  CHARACTER(LEN=50) , PRIVATE :: nan_list(50)       ! Updated by scm trap nan

  INTEGER :: time_info(7)

  ! Integer id for time and date variables in netcdf files
  INTEGER(incdf) :: scm_nc_time_id
  INTEGER(incdf) :: scm_nc_date_id
  INTEGER(incdf) :: scm_nc_local_id
  INTEGER(incdf) :: scm_nc_hour_id
  INTEGER        :: scm_timestep_count ! Number of current timestep
  REAL           :: scm_timestep       ! SCM run timestep

  INTEGER ::    &
    rw          &! Rows
  , rw_lng      &! Row length
  , nmod_lv     &! No. of model levels
  , nwet_lv     &! No. of wet levels
  , nbl_lv      &! No. of boundary layer levels
  , o3_lv       &! No. of ozone levels
  , nlnd        &! No. of land points
  , ntile       &! No. of tiles per land point
  , ntyp        &! No. of plant profiles
  , st_lv       &! No. of soil temperature levels
  , sm_lv       &! No. of soil moisture levels
  , ntr_lv      &! No. of tracer levels
  , ntr_var      ! No. of tracers


  INTEGER ::    &
    nml_nmod_lv &! No. of model levels specified in obs forcing
  , nobs        &! No. of observational forcing profiles in namelist
  , obs_t0       ! Observational profile at which to begin run

  REAL, ALLOCATABLE :: &
    z_th(:)     &! Height of theta levels
  , z_rh(:)     &! Height of rho   levels
  , nml_z_th(:) &
  , nml_z_rh(:)

  LOGICAL ::              &
    cv_run_mess = .FALSE. &! Display runtime convective messages
  , nc_obsfor   = .FALSE. &! 
  , old_nml     = .FALSE. &! Old namelist format (pre vn7.7)
  , old_rlx     = .FALSE. &! Old forcing relaxation code (pre vn7.7)
  , old_vertadv = .FALSE. &! Old interactive vertical advection code (pre vn7.7)
  , interpolate = .FALSE.

  ! Parameters for forcing relaxation case options
  INTEGER, PARAMETER ::   &
    rlx_none      = 0     &! No relaxation
  , rlx_init      = 1     &! Relax to initial profile across tau
  , rlx_bgrd      = 2     &! Relax to background profile across tau
  , rlx_inst_init = 3     &! Reset to initial profile
  , rlx_inst_bgrd = 4      ! Reset to background profile


  ! Dr Hook Parameters
  !=================================
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1


CONTAINS

  SUBROUTINE scm_message(string, unit)

    IMPLICIT NONE

    CHARACTER(LEN=*), OPTIONAL, INTENT(In) :: string
    INTEGER,      OPTIONAL, INTENT(In) :: unit

    CHARACTER(LEN=10) :: ts_str
    INTEGER       :: unitno

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('SCM_MESSAGE',zhook_in,zhook_handle)

    IF (PRESENT(unit)) THEN
      unitno = unit
    ELSE
      unitno = 6
    END IF

    WRITE(ts_str,'(A4,I4,A2)') '[TS ', scm_timestep_count,'] '

    IF (PRESENT(string)) THEN
      WRITE(unitno,*) ts_str, TRIM(ADJUSTL(string))
    ELSE
      WRITE(unitno,*) ts_str
    END IF

    IF (lhook) CALL dr_hook('SCM_MESSAGE',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE scm_message

!-----------------------------------------------------------------------------

  SUBROUTINE scm_trap_nan(string, routine)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(In) :: string  ! Diagnostic with NaN
    CHARACTER(LEN=*), INTENT(In) :: routine ! Calling routine
    INTEGER                  :: i

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('SCM_TRAP_NAN',zhook_in,zhook_handle)

    IF (nan_list_count <= 50) THEN

      DO i=1, nan_list_count
        IF (nan_list(i) == string) THEN
          RETURN
        END IF
      END DO

      IF (nan_list_count == 50) THEN
        CALL scm_message ('50+ Diagnostics contain NaNs')
        nan_list_count = nan_list_count + 1
      ELSE
        CALL scm_message ('NaN detected: '//string//' from '//routine)
        nan_list(nan_list_count+1) = string
        nan_list_count = nan_list_count + 1
      END IF
    END IF

    IF (lhook) CALL dr_hook('SCM_TRAP_NAN',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE scm_trap_nan

!=============================================================================

  SUBROUTINE alloc_common_scm

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('ALLOC_COMMON_SCM',zhook_in,zhook_handle)

    ALLOCATE                    &
      ( z_th(nmod_lv)           &
      , z_rh(nmod_lv+1)         &
      , nml_z_th(nml_nmod_lv)   &
      , nml_z_rh(nml_nmod_lv+1) )

    IF (lhook) CALL dr_hook('ALLOC_COMMON_SCM',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE alloc_common_scm

!-------------------------------------------------------------------------------

  SUBROUTINE dealloc_common_scm

    IMPLICIT NONE

    ! Dr Hook
    !=============================================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('DEALLOC_COMMON_SCM',zhook_in,zhook_handle)

    DEALLOCATE   &
      ( z_th     &
      , z_rh     &
      , nml_z_th &
      , nml_z_rh )

    IF (lhook) CALL dr_hook('DEALLOC_COMMON_SCM',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE dealloc_common_scm


!-------------------------------------------------------------------------------
END MODULE scm_utils
