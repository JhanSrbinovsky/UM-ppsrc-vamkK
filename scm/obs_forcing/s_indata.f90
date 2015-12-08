! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : INDATA

MODULE s_indata

  USE scm_cntl_mod, ONLY: scm_nml

  USE s_maxdim, ONLY:                                                &
    mx_rw_lng, mx_rw, mx_nlnd

  USE scm_utils, ONLY:                                               &
    imdi, rmdi, rw_lng, rw, nlnd, zhook_in, zhook_out                &
  , jprb, lhook, dr_hook

  USE s_main_force, ONLY:                                            &
    salt_dim1, salt_dim2, salt_dim3                                  &
  , year_init, month_init, day_init, hour_init, min_init, sec_init   &
  , tapeyear_init, tapemonth_init, tapeday_init                      &
  , tapehour_init, tapemin_init,   tapesec_init, gather              &
  , scm_soil_type    => soil_type                                    &
  , scm_gridbox_area => gridbox_area                                 &
  , scm_lat          => lat                                          &
  , scm_long         => long

  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in INDATA namelist from forcing file, scm_nml. INDATA
!   contains initial data required to run the model.
!
! Method:
!   Namelist INDATA is defined in this module and read in by contained
!   subroutine read_indata.  Scalar variables are read directly to
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
  INTEGER, PRIVATE ::                           &
    soil_type    (mx_nlnd) = imdi

  REAL, PRIVATE ::                              &
    gridbox_area (mx_rw_lng*mx_rw) = rmdi       &
  , lat          (mx_rw_lng*mx_rw) = rmdi       &
  , long         (mx_rw_lng*mx_rw) = rmdi


  !---------------------------------------------------------------------------
  ! Define namelist
  !---------------------------------------------------------------------------
  NAMELIST/indata/                                                            &
    soil_type, tapeyear_init, tapemonth_init, tapeday_init, tapehour_init     &
  , tapemin_init, tapesec_init, year_init, month_init, day_init, hour_init    &
  , min_init, sec_init, lat, long, gridbox_area, salt_dim1, salt_dim2         &
  , salt_dim3, gather

  PRIVATE :: indata

!=============================================================================
CONTAINS

  SUBROUTINE read_indata

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    ! Local variables
    INTEGER :: istatus
    INTEGER :: icode
    CHARACTER(LEN=11), PARAMETER :: routinename='read_indata'

    IF (lhook) CALL dr_hook('READ_INDATA',zhook_in,zhook_handle)

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, indata)
    CLOSE (10)

    IF (lat(1)          /= rmdi) scm_lat                    &
      = RESHAPE( lat,          (/rw_lng,rw/))
    IF (long(1)         /= rmdi) scm_long                   &
      = RESHAPE( long,         (/rw_lng,rw/))
    IF (gridbox_area(1) /= rmdi) scm_gridbox_area           &
      = RESHAPE( gridbox_area, (/rw_lng,rw/))

    IF (nlnd > 0) THEN
      IF (soil_type(1)  /= imdi) scm_soil_type              &
        = RESHAPE( soil_type, (/nlnd/))
    END IF

    IF (lhook) CALL dr_hook('READ_INDATA',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_indata

!=============================================================================
END MODULE s_indata
