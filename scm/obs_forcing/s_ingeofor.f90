! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : INGEOFOR

MODULE s_ingeofor

  USE scm_cntl_mod, ONLY: scm_nml

  USE s_maxdim, ONLY:                                   &
    mx_rw_lng, mx_rw, mx_nobs, mx_mod_lv

  USE scm_utils, ONLY:                                  &
    rmdi, rw_lng, rw, nobs, nmod_lv                     &
  , zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_main_force, ONLY:                               &
    ug_opt, vg_opt                                      &
  , scm_ug => ug                                        &
  , scm_vg => vg

  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in INGEOFOR namelist from forcing file, scm_nml.
!   INGEOFOR contains information related to geostrophic wind forcing.
!
! Method:
!   Namelist INGEOFOR is defined in this module and read in by contained
!   subroutine read_ingeofor.  Scalar variables are read directly to
!   s_main_force.  Array variables are declared as private single rank arrays
!   of fixed length (arrays cannot be varying length in namelists) and
!   reshaped before being transferred to arrays of the correct size/shape in
!   s_main_force.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================
!
! Declarations:

  !--------------------------------------------------
  ! Fixed arrays to read in namelist before reshaping
  !--------------------------------------------------
  REAL, PRIVATE ::                                &
    ug (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv) = rmdi &! Geostrophic u-wind (m/s)
  , vg (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv) = rmdi  ! Geostrophic v-wind (m/s)


  !---------------------------------------------------------------------------
  ! Define namelist
  !---------------------------------------------------------------------------
  NAMELIST/ingeofor/         &
    ug, vg, ug_opt, vg_opt

  PRIVATE :: ingeofor

!=============================================================================
CONTAINS

  SUBROUTINE read_ingeofor

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    INTEGER :: istatus
    INTEGER :: icode

    INTEGER :: i, j, j2, k

    INTEGER, PARAMETER ::       &
      constant_value       = 0  &
    , constant_profile     = 1  &
    , time_varying_profile = 2

    CHARACTER(LEN=11), PARAMETER :: routinename='read_ingeofor'

    IF (lhook) CALL dr_hook('READ_INGEOFOR',zhook_in,zhook_handle)

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, ingeofor)
    CLOSE (10)


    ! Geostrophic - Zonal wind
    !---------------------------------------
    IF (ug(1) /= rmdi) THEN
      SELECT CASE (ug_opt)

      CASE (constant_value)
        ! Single value ug applied everywhere
        DO k=1, nmod_lv
          DO j2=1, nobs
            DO j=1, rw
              DO i=1, rw_lng
                scm_ug(i,j,j2,k) = ug(1)
              END DO
            END DO
          END DO
        END DO

      CASE (constant_profile)
        ! Single profile, non-time varying
        DO k=1, nmod_lv
          DO j2=1, nobs
            DO j=1, rw
              DO i=1, rw_lng
                scm_ug(i,j,j2,k) = ug(k)
              END DO
            END DO
          END DO
        END DO

      CASE (time_varying_profile)
        ! Time varying profile
        scm_ug = RESHAPE(ug ,(/rw_lng, rw, nobs, nmod_lv/))
      END SELECT
    END IF


    ! Geostrophic - Meridional wind
    !---------------------------------------
    IF (vg(1) /= rmdi) THEN
      SELECT CASE (vg_opt)

      CASE (constant_value)
        ! Single value vg applied everywhere
        DO k=1, nmod_lv
          DO j2=1, nobs
            DO j=1, rw
              DO i=1, rw_lng
                scm_vg(i,j,j2,k) = vg(1)
              END DO
            END DO
          END DO
        END DO

      CASE (constant_profile)
        ! Single profile, non-time varying
        DO k=1, nmod_lv
          DO j2=1, nobs
            DO j=1, rw
              DO i=1, rw_lng
                scm_vg(i,j,j2,k) = vg(k)
              END DO
            END DO
          END DO
        END DO

      CASE (time_varying_profile)
        ! Time varying profile
        scm_vg = RESHAPE(vg ,(/rw_lng, rw, nobs, nmod_lv/))
      END SELECT
    END IF


    IF (lhook) CALL dr_hook('READ_INGEOFOR',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_ingeofor

!=============================================================================
END MODULE s_ingeofor
