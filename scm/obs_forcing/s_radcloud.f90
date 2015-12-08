! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : RACLOUD

MODULE s_radcloud

  USE scm_cntl_mod, ONLY: scm_nml

  USE s_maxdim,  ONLY:                                       &
    mx_rw_lng, mx_rw, mx_wet_lv

  USE scm_utils, ONLY:                                       &
    rmdi, imdi, rw_lng, rw, nwet_lv, nml_nmod_lv, z_th, z_rh &
  , nml_z_th, nml_z_rh, interpolate                          &
  , zhook_in, zhook_out, jpim, jprb, lhook, dr_hook


  USE s_interp_mod, ONLY: interp1d
  USE s_main_force, ONLY:                                    &
    scm_iccb_rad        => iccb_rad                          &
  , scm_icct_rad        => icct_rad                          &
  , scm_cca_rad         => cca_rad                           &
  , scm_layer_cloud_rad => layer_cloud_rad                   &
  , scm_qcl_rad         => qcl_rad                           &
  , scm_qcf_rad         => qcf_rad                           &
  , scm_ccwpin_rad      => ccwpin_rad

  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in RADCLOUD namelist from forcing file, <scm_nml>.
!   RADCLOUD contains information related to fixed radiative clouds, it is
!   only use when the radcloud_fixed logical is enabled.
!
! Method:
!   Namelist RADCLOUD is defined in this module and read in by contained
!   subroutine read_radcloud.  Scalar variables are read directly to
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

  INTEGER, PRIVATE ::                &
    iccb_rad(mx_rw_lng*mx_rw) = imdi &
  , icct_rad(mx_rw_lng*mx_rw) = imdi

  REAL, PRIVATE ::                                     &
    cca_rad         (mx_rw_lng*mx_rw) = rmdi           &
  , ccwpin_rad      (mx_rw_lng*mx_rw) = rmdi           &
  , qcl_rad         (mx_rw_lng*mx_rw*mx_wet_lv) = rmdi &
  , qcf_rad         (mx_rw_lng*mx_rw*mx_wet_lv) = rmdi &
  , layer_cloud_rad (mx_rw_lng*mx_rw*mx_wet_lv) = rmdi

  REAL, ALLOCATABLE ::          &
    nml_layer_cloud_rad (:,:,:) &
  , nml_qcl_rad         (:,:,:) &
  , nml_qcf_rad         (:,:,:)


  NAMELIST/radcloud/                                                          &
    cca_rad, iccb_rad, icct_rad, layer_cloud_rad, qcl_rad, qcf_rad,           &
    ccwpin_rad

  PRIVATE :: radcloud

!=============================================================================
CONTAINS

  SUBROUTINE read_radcloud

    IMPLICIT NONE

    INTEGER :: istatus
    INTEGER :: icode

    CHARACTER(LEN=11), PARAMETER :: routinename='read_radcloud'

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('READ_RADCLOUD',zhook_in,zhook_handle)

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, radcloud)
    CLOSE (10)


    IF (icct_rad(1) /= imdi)                                                  &
      scm_icct_rad  = RESHAPE( icct_rad,  (/rw_lng,rw/) )

    IF (iccb_rad(1) /= imdi)                                                  &
      scm_iccb_rad  = RESHAPE( iccb_rad,  (/rw_lng,rw/) )

    IF (cca_rad(1)   /= rmdi)                                                 &
      scm_cca_rad    = RESHAPE( cca_rad,    (/rw_lng,rw/) )

    IF (ccwpin_rad(1) /= rmdi)                                                &
      scm_ccwpin_rad = RESHAPE( ccwpin_rad, (/rw_lng,rw/) )

    IF (interpolate) THEN

      CALL alloc_radcloud
      CALL interp_radcloud
      CALL dealloc_radcloud

    ELSE

      IF (layer_cloud_rad(1) /= rmdi)                                         &
        scm_layer_cloud_rad = RESHAPE( layer_cloud_rad                        &
                                     , (/rw_lng,rw,nwet_lv/) )

      IF (qcl_rad(1) /= rmdi)                                                 &
        scm_qcl_rad = RESHAPE( qcl_rad, (/rw_lng,rw,nwet_lv/) )

      IF (qcl_rad(1) /= rmdi)                                                 &
        scm_qcf_rad = RESHAPE( qcf_rad, (/rw_lng,rw,nwet_lv/) )

    END IF

    IF (lhook) CALL dr_hook('READ_RADCLOUD',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_radcloud

!-----------------------------------------------------------------------------

  SUBROUTINE interp_radcloud

    IMPLICIT NONE

    ! Dr Hook
    !=============================================================
    REAL(KIND=jprb) :: zhook_handle

    ! Local variables
    INTEGER :: i,j

    IF (lhook) CALL dr_hook('INTERP_RADCLOUD',zhook_in,zhook_handle)

    IF (layer_cloud_rad(1) /= rmdi) THEN

      nml_layer_cloud_rad = RESHAPE( layer_cloud_rad                          &
                                   , (/rw_lng,rw,nml_nmod_lv/) )
      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th(:), nml_layer_cloud_rad(i,j,:)                         &
            , z_th,        scm_layer_cloud_rad(i,j,:) )
        END DO
      END DO

    END IF


    IF (qcl_rad(1) /= rmdi) THEN

      nml_qcl_rad = RESHAPE( qcl_rad, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th(:), nml_qcl_rad(i,j,:)                                 &
            , z_th,        scm_qcl_rad(i,j,:) )
        END DO
      END DO

    END IF


    IF (qcf_rad(1) /= rmdi) THEN

      nml_qcf_rad = RESHAPE( qcf_rad, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th(:), nml_qcf_rad(i,j,:)                                 &
            , z_th,        scm_qcf_rad(i,j,:) )
        END DO
      END DO

    END IF

    IF (lhook) CALL dr_hook('INTERP_RADCLOUD',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE interp_radcloud

!-----------------------------------------------------------------------------

  SUBROUTINE alloc_radcloud

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('ALLOC_RADCLOUD',zhook_in,zhook_handle)

    ALLOCATE                                          &
      ( nml_layer_cloud_rad (rw_lng,rw,nml_nmod_lv)   &
      , nml_qcl_rad         (rw_lng,rw,nml_nmod_lv)   &
      , nml_qcf_rad         (rw_lng,rw,nml_nmod_lv) )

    IF (lhook) CALL dr_hook('ALLOC_RADCLOUD',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE alloc_radcloud

!-----------------------------------------------------------------------------

  SUBROUTINE dealloc_radcloud

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('DEALLOC_RADCLOUD',zhook_in,zhook_handle)

    DEALLOCATE               &
      ( nml_layer_cloud_rad  &
      , nml_qcl_rad          &
      , nml_qcf_rad )

    IF (lhook) CALL dr_hook('DEALLOC_RADCLOUD',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE dealloc_radcloud

!-----------------------------------------------------------------------------
END MODULE s_radcloud
