! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************

!+ Subroutine to set the parametrization schemes for clouds.

! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of code:
!   fortran 77  with extensions listed in documentation.

!- ---------------------------------------------------------------------
      SUBROUTINE r2_set_cloud_parametrization(ierr, n_band              &
         , i_st_water, i_cnv_water, i_st_ice, i_cnv_ice                 &
         , l_drop_type, i_drop_parametrization                          &
         , n_drop_phf_term, drop_parameter_list                         &
         , drop_parm_min_dim, drop_parm_max_dim                         &
         , l_ice_type, i_ice_parametrization                            &
         , n_ice_phf_term, ice_parameter_list                           &
         , ice_parm_min_dim, ice_parm_max_dim                           &
         , i_condensed_param, condensed_n_phf, condensed_param_list     &
         , condensed_min_dim, condensed_max_dim                         &
         , nd_band, nd_drop_type, nd_ice_type, nd_cloud_parameter       &
         , nd_cloud_component                                           &
         )


      USE rad_pcf
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE


!     dummy arguments:

      INTEGER                                                           &
                !, intent(out)
           ierr
!             error flag

!     sizes of arrays:
      INTEGER                                                           &
                !, intent(in)
           nd_band                                                      &
!             maximum number of spectral bands
         , nd_drop_type                                                 &
!             maximum number of types of droplets
         , nd_ice_type                                                  &
!             maximum number of types of ice crystals
         , nd_cloud_parameter                                           &
!             maximum number of parameters for clouds
         , nd_cloud_component
!             size allocated for components in clouds

      INTEGER                                                           &
                !, intent(in)
           n_band
!             number of spectral bands

!     types of droplets and crystals:
      INTEGER                                                           &
                !, intent(in)
           i_st_water                                                   &
!             type of water droplets in stratiform clouds
         , i_cnv_water                                                  &
!             type of water droplets in convective clouds
         , i_st_ice                                                     &
!             type of ice crystals in stratiform clouds
         , i_cnv_ice
!             type of ice crystals in convective clouds

      LOGICAL                                                           &
                !, intent(in)
           l_drop_type(nd_drop_type)                                    &
!             flags for types of droplet present
         , l_ice_type(nd_ice_type)
!             flags for types of ice crystal present
      INTEGER                                                           &
                !, intent(in)
           i_drop_parametrization(nd_drop_type)                         &
!             parametrizations of types of droplets
         , n_drop_phf_term(nd_drop_type)                                &
!             Number of terms in the phase function for
!             droplets
         , i_ice_parametrization(nd_ice_type)                           &
!             parametrizations of types of ice crystals
         , n_ice_phf_term(nd_ice_type)
!             Number of terms in the phase function for
!             ice crystals
      REAL                                                              &
                !, intent(in)
           drop_parameter_list(nd_cloud_parameter                       &
              , nd_band, nd_drop_type)                                  &
!             parameters for optical parametrizations of droplets
         , drop_parm_min_dim(nd_drop_type)                              &
!             minimum size of droplets permitted in parametrizations
         , drop_parm_max_dim(nd_drop_type)                              &
!             maximum size of droplets permitted in parametrizations
         , ice_parameter_list(nd_cloud_parameter                        &
              , nd_band, nd_ice_type)                                   &
!             parameters for optical parametrizations of ice crystals
         , ice_parm_min_dim(nd_ice_type)                                &
!             minimum size of ice crystals permitted in parametrizations
         , ice_parm_max_dim(nd_ice_type)
!             maximum size of ice crystals permitted in parametrizations

      INTEGER                                                           &
                !, intent(out)
           i_condensed_param(nd_cloud_component)                        &
!             types of parametrization used for condensed
!             components in clouds
         , condensed_n_phf(nd_cloud_component)
!             Number of terms in the phase function
      REAL                                                              &
                !, intent(out)
           condensed_param_list(nd_cloud_parameter                      &
              , nd_cloud_component, nd_band)                            &
!             coefficients for parametrization of condensed phases
         , condensed_min_dim(nd_cloud_component)                        &
!             minimum dimension of each condensed component
         , condensed_max_dim(nd_cloud_component)
!             maximum dimension of each condensed component


!     local variables:
      INTEGER                                                           &
           i                                                            &
!             loop variable
         , j                                                            &
!             loop variable
         , i_scheme
!             parametrization scheme

!     Functions called:
      INTEGER, EXTERNAL :: set_n_cloud_parameter
!       Function to find number of parameters for clouds


      CHARACTER (LEN=*), PARAMETER :: RoutineName =                     &
        'r2_set_cloud_parametrization'
      CHARACTER (LEN=240) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook(                                          &
        'R2_SET_CLOUD_PARAMETRIZATION', zhook_in, zhook_handle)


!     Select parametrization for water in stratiform clouds:
      IF ( (i_st_water <= nd_drop_type).AND.                            &
           (l_drop_type(i_st_water)) ) THEN
         i_scheme=i_drop_parametrization(i_st_water)
         i_condensed_param(ip_clcmp_st_water)=i_scheme
         condensed_n_phf(ip_clcmp_st_water)=n_drop_phf_term(i_st_water)
         condensed_min_dim(ip_clcmp_st_water)                           &
            =drop_parm_min_dim(i_st_water)
         condensed_max_dim(ip_clcmp_st_water)                           &
            =drop_parm_max_dim(i_st_water)
      ELSE
         cmessage = '*** error: no data exist for type ' //             &
           'of droplet selected in stratiform water clouds.'
         ierr=i_err_fatal
         GO TO 9999
      END IF

      DO i=1, n_band
! DEPENDS ON: set_n_cloud_parameter
         DO j=1, set_n_cloud_parameter(i_scheme                         &
            , ip_clcmp_st_water, condensed_n_phf(ip_clcmp_st_water))
            condensed_param_list(j, ip_clcmp_st_water, i)               &
               =drop_parameter_list(j, i, i_st_water)
         END DO
      END DO


!     Select parametrization for water in convective clouds:
      IF ( (i_cnv_water <= nd_drop_type).AND.                           &
           (l_drop_type(i_cnv_water)) ) THEN
         i_scheme=i_drop_parametrization(i_cnv_water)
         i_condensed_param(ip_clcmp_cnv_water)=i_scheme
         condensed_n_phf(ip_clcmp_cnv_water)                            &
            =n_drop_phf_term(i_cnv_water)
         condensed_min_dim(ip_clcmp_cnv_water)                          &
            =drop_parm_min_dim(i_cnv_water)
         condensed_max_dim(ip_clcmp_cnv_water)                          &
            =drop_parm_max_dim(i_cnv_water)
      ELSE
         cmessage = '*** error: no data exist for type ' //             &
           'of crystal selected in convective water clouds.'
         ierr=i_err_fatal
         GO TO 9999
      END IF

      DO i=1, n_band
         DO j=1, set_n_cloud_parameter(i_scheme                         &
            , ip_clcmp_cnv_water, condensed_n_phf(ip_clcmp_cnv_water))
            condensed_param_list(j, ip_clcmp_cnv_water, i)              &
               =drop_parameter_list(j, i, i_cnv_water)
         END DO
      END DO


!     Select parametrization for ice in stratiform clouds:
      IF ( (i_st_ice <= nd_ice_type).AND.                               &
           (l_ice_type(i_st_ice)) ) THEN
         i_scheme=i_ice_parametrization(i_st_ice)
         i_condensed_param(ip_clcmp_st_ice)=i_scheme
         condensed_n_phf(ip_clcmp_st_ice)=n_ice_phf_term(i_st_ice)
         condensed_min_dim(ip_clcmp_st_ice)                             &
            =ice_parm_min_dim(i_st_ice)
         condensed_max_dim(ip_clcmp_st_ice)                             &
            =ice_parm_max_dim(i_st_ice)
      ELSE
         cmessage = '*** error: no data exist for type ' //             &
           'of crystal selected in stratiform ice clouds.'
         ierr=i_err_fatal
         GO TO 9999
      END IF

      DO i=1, n_band
         DO j=1, set_n_cloud_parameter(i_scheme                         &
            , ip_clcmp_st_ice, condensed_n_phf(ip_clcmp_st_ice))
            condensed_param_list(j, ip_clcmp_st_ice, i)                 &
               =ice_parameter_list(j, i, i_st_ice)
         END DO
      END DO


!     Select parametrization for ice in convective clouds:
      IF ( (i_cnv_ice <= nd_ice_type).AND.                              &
           (l_ice_type(i_cnv_ice)) ) THEN
         i_scheme=i_ice_parametrization(i_cnv_ice)
         i_condensed_param(ip_clcmp_cnv_ice)=i_scheme
         condensed_n_phf(ip_clcmp_cnv_ice)=n_ice_phf_term(i_cnv_ice)
         condensed_min_dim(ip_clcmp_cnv_ice)                            &
            =ice_parm_min_dim(i_cnv_ice)
         condensed_max_dim(ip_clcmp_cnv_ice)                            &
            =ice_parm_max_dim(i_cnv_ice)
      ELSE
         cmessage = '*** error: no data exist for type ' //             &
           'of crystal selected in convective ice clouds.'
         ierr=i_err_fatal
         GO TO 9999
      END IF

      DO i=1, n_band
         DO j=1, set_n_cloud_parameter(i_scheme                         &
            , ip_clcmp_cnv_ice, condensed_n_phf(ip_clcmp_cnv_ice))
            condensed_param_list(j, ip_clcmp_cnv_ice, i)                &
               =ice_parameter_list(j, i, i_cnv_ice)
         END DO
      END DO


 9999 CONTINUE
! Check error condition
      IF (ierr /= i_normal) THEN

        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook(                                          &
        'R2_SET_CLOUD_PARAMETRIZATION', zhook_out, zhook_handle)
      END SUBROUTINE r2_set_cloud_parametrization
