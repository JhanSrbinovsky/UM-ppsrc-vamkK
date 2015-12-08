! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to aggregate clouds into regions.
!
! Method:
!       The clouds in a layer are combined in groups to form regions
!       which will be considered as bulk entities in the solution of the
!       equation of transfer. The extents of these regions are also
!       determined.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
SUBROUTINE aggregate_cloud(ierr                                         &
     , n_profile, n_layer, n_cloud_top                                  &
     , i_cloud, i_cloud_representation, n_cloud_type                    &
     , frac_cloud                                                       &
     , i_region_cloud, frac_region                                      &
     , nd_profile, nd_layer, nd_cloud_type, nd_region                   &
     , id_ct                                                            &
     )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport

  IMPLICIT NONE


! Dummy array sizes
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , nd_cloud_type                                                     &
!       Maximum number of types of cloud
    , nd_region                                                         &
!       Maximum number of cloudy regions
    , id_ct
!       Topmost declared cloudy layer

! Dummy variables.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top
!       Topmost cloudy layer
  INTEGER, INTENT(IN) ::                                                &
      i_cloud                                                           &
!       Cloud scheme used
    , i_cloud_representation                                            &
!       Representation of clouds used
    , n_cloud_type
!       Number of types of cloud

  REAL (RealK), INTENT(IN) ::                                           &
      frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)
!       Fractions of each type of cloud

  INTEGER, INTENT(OUT) ::                                               &
      i_region_cloud(nd_cloud_type)
!       Regions in which particular types of cloud fall
  REAL (RealK), INTENT(OUT) ::                                          &
      frac_region(nd_profile, id_ct: nd_layer, nd_region)
!       Fractions of total cloud occupied by each region


! Local variables
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , k
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=256)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'aggregate_cloud'


  IF (lhook) CALL dr_hook('AGGREGATE_CLOUD',zhook_in,zhook_handle)

  IF ( (i_cloud == ip_cloud_triple).OR.                                 &
       (i_cloud == ip_cloud_part_corr_cnv) ) THEN

    IF (i_cloud_representation == ip_cloud_csiw) THEN

      DO k=1, n_cloud_type
        IF (k == ip_cloud_type_sw) THEN
          i_region_cloud(k)=ip_region_strat
        ELSE IF (k == ip_cloud_type_si) THEN
          i_region_cloud(k)=ip_region_strat
        ELSE IF (k == ip_cloud_type_cw) THEN
          i_region_cloud(k)=ip_region_conv
        ELSE IF (k == ip_cloud_type_ci) THEN
          i_region_cloud(k)=ip_region_conv
        END IF
      END DO

      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          frac_region(l, i, ip_region_strat)                            &
            =frac_cloud(l, i, ip_cloud_type_sw)                         &
            +frac_cloud(l, i, ip_cloud_type_si)
          frac_region(l, i, ip_region_conv)                             &
            =frac_cloud(l, i, ip_cloud_type_cw)                         &
            +frac_cloud(l, i, ip_cloud_type_ci)
        END DO
      END DO

    ELSE IF (i_cloud_representation == ip_cloud_conv_strat) THEN

      DO k=1, n_cloud_type
        IF (k == ip_cloud_type_strat) THEN
          i_region_cloud(k)=ip_region_strat
        ELSE IF (k == ip_cloud_type_conv) THEN
          i_region_cloud(k)=ip_region_conv
        END IF
      END DO

      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          frac_region(l, i, ip_region_strat)                            &
            =frac_cloud(l, i, ip_cloud_type_strat)
          frac_region(l, i, ip_region_conv)                             &
            =frac_cloud(l, i, ip_cloud_type_conv)
         END DO
      END DO

    ELSE
      cmessage = '*** Error: This representation of clouds is not '//   &
        'compatible with separate convective and stratiform overlap.'
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF

  END IF


  IF (lhook) CALL dr_hook('AGGREGATE_CLOUD',zhook_out,zhook_handle)

END SUBROUTINE aggregate_cloud
