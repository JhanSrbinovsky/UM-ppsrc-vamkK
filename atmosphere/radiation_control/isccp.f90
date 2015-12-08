! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Routine to diagnose ISCCP cloud types

! Purpose:
!   This subroutine diagnoses ISCCP cloud types

! Method:
!   Argument list mirrors that of R2_SWRAD

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.

!- --------------------------------------------------------------------
      SUBROUTINE isccp(ierr, i_segment                                  &
         , h2o                                                          &
!                       Pressure Fields
         , tstar, pstar, p_layer_boundaries, p_layer_centres            &
!                       Temperatures
         , tac                                                          &
!                       Options for treating clouds
         , lca_area, lca_bulk                                           &
!                       Convective Cloud Fields
         , cca, cccwp, ccb, cct                                         &
!                       Solar Fields
         , lit, list, scs                                               &
         , trindx_fld                                                   &
!                       General Diagnostics
         , lw_diag, sw_diag, row_list, col_list                         &
!                       Physical Dimensions
         , nlit                                                         &
         , n_points, nlevs, n_layer, nclds                              &
!          Number of points in segment, number of levels,
!          cloud levels
         , nwet, nozone                                                 &
         , npd_field, npd_profile, npd_layer, npd_column                &
         , n_cca_lev                                                    &
           ! Variables needed to calculate layer masses
         , rho_r2, r_rho_levels, r_theta_levels                         &
         , q, qcl, qcf, qcf2, qrain, qgraup                             &
         , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
         )


      USE rad_pcf
      USE rad_com_mod, ONLY: is_ncol
      USE swrdiag_mod, ONLY:                                            &
          strswdiag

      USE lwrdiag_mod, ONLY:                                            &
          strlwdiag

      USE cv_run_mod, ONLY:                                             &
          l_3d_cca

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE



!     DIMENSIONS OF ARRAYS:
      INTEGER ierr   ! error code

      INTEGER                                                           &
                !, INTENT(IN)
           npd_field                                                    &
!             FIELD SIZE IN CALLING PROGRAM
         , npd_profile                                                  &
!             SIZE OF ARRAY OF PROFILES
         , npd_layer                                                    &
!             ARRAY SIZES FOR LAYERS
         , npd_column
!             NUMBER OF COLUMNS PER POINT

      INTEGER                                                           &
           top_height                                                   &
!             FLAG TO SPECIFY METHOD FOR ISCCP CLOUD TOP HEIGHT
         , overlap
!             FLAG TO SPECIFY METHOD FOR ISCCP CLOUD OVERLAP



!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
           n_points                                                     &
!             Total number of points, including unlit ones
         , nwet                                                         &
!             NUMBER OF WET LEVELS
         , nozone                                                       &
!             NUMBER OF LEVELS WITH OZONE
         , nlevs                                                        &
!             NUMBER OF ATMOSPHERIC LAYERS
         , n_layer                                                      &
!             Number of layers seen in the radiation scheme
         , nclds                                                        &
!             NUMBER OF CLOUDY LEVELS
         , n_cca_lev
!             NUMBER OF CONVECTIVE CLOUD LEVELS




!     GASEOUS MIXING RATIOS
      REAL                                                              &
                !, INTENT(IN)
           h2o(npd_field, nwet)
!             MASS MIXING RATIO OF WATER

!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL                                                              &
                !, INTENT(IN)
           tstar(npd_field)                                             &
!             SURFACE TEMPERATURES
         , pstar(npd_field)                                             &
!             SURFACE PRESSURES
         , p_layer_boundaries(npd_field,0:nlevs)                        &
!            PRESSURE AT BOUNDARIES OF LAYERS
         , p_layer_centres(npd_field,0:nlevs)                           &
!            PRESSURE AT CENTRES OF LAYERS
         , tac(npd_field, nlevs)
!             TEMPERATURES AT CENTRES OF LAYERS

!     INCIDENT SOLAR RADIATION:
      INTEGER                                                           &
                !, INTENT(IN)
           nlit                                                         &
!             NUMBER OF LIT POINTS
         , list(npd_field)
!             LIST OF LIT POINTS
      REAL                                                              &
                !, INTENT(IN)
           scs                                                          &
!             SCALING OF SOLAR INCIDENT FIELD
         , lit(npd_field)
!             FRACTION OF TIME POINT IS LIT

      REAL                                                              &
                !, INTENT(IN)
           lca_area(npd_field, nclds+1/(nclds+1))                       &
!             AREA FRACTIONS OF LAYER CLOUDS OUTSIDE CONVECTIVE TOWERS
         , lca_bulk(npd_field, nclds+1/(nclds+1))
!             BULK FRACTIONS OF LAYER CLOUDS OUTSIDE CONVECTIVE TOWERS

!     PROPERTIES OF CONVECTIVE CLOUDS:
      INTEGER                                                           &
                !, INTENT(IN)
           ccb(npd_field)                                               &
!             BASE OF CONVECTIVE CLOUD
         , cct(npd_field)
!             TOP OF CONVECTIVE CLOUD
      REAL                                                              &
                !, INTENT(IN)
           cccwp(npd_field)                                             &
!             WATER PATH OF CONVECTIVE CLOUD
         , cca(npd_field,n_cca_lev)
!             FRACTION OF CONVECTIVE CLOUD
      LOGICAL                                                           &
                !, INTENT(IN)
           l_mcr_qcf2                                                   &
                          ! Use second ice category
      ,    l_mcr_qrain                                                  &
                          ! Use prognostic rain
      ,    l_mcr_qgraup                                                 &
                          ! Use graupel
      ,    l_mixing_ratio ! Use mixing ratios in layer mass calculation

!                       Level of tropopause
      INTEGER                                                           &
           trindx_fld(npd_field)
!             THE LAYER BOUNDARY OF THE TROPOPAUSE
!     Information for the calculation of layer masses
      REAL, INTENT(IN)::                                                &
        rho_r2(npd_field,nlevs)                                         &
                                ! Air density*radius of earth**2 / kg m-1
      , r_rho_levels(npd_field,nlevs)                                   &
                                      ! Height of rho levels / m
      , r_theta_levels(npd_field,0:nlevs)                               &
                                           ! Height of theta levels / m
      , q(npd_field,nwet)                                               &
                                ! Water vapour mixing ratio / kg kg-1
      , qcl(npd_field,nwet)                                             &
                                ! Liquid water mixing ratio / kg kg-1
      , qcf(npd_field,nwet)                                             &
                                ! Ice mixing ratio / kg kg-1
      , qcf2(npd_field,nwet)                                            &
                                ! Second ice category mr / kg kg-1
      , qrain(npd_field,nwet)                                           &
                                ! Rain mixing ratio / kg kg-1
      , qgraup(npd_field,nwet)  ! Graupel mixing ratio / kg kg-1

!     CALCULATED FLUXES:
      REAL                                                              &
                !, INTENT(OUT)
           swsea(npd_field)
!             SEA-SURFACE COMPONENTS OF FLUX

!     DIAGNOSTICS:

!     Definition of the diagnostic structure

      TYPE (strlwdiag) :: lw_diag
      TYPE (strswdiag) :: sw_diag

      INTEGER, INTENT(IN) :: row_list(npd_field)
!                              List of row indices of lit points
      INTEGER, INTENT(IN) :: col_list(npd_field)
!                              List of column indices of lit points


!     CALCULATED LOCAL DIAGNOSTICS:
      REAL                                                              &
                !, INTENT(OUT)
           pfull_fld(npd_field, nclds)                                  &
         , phalf_fld(npd_field, nclds+1)                                &
         , qv_fld(npd_field, nclds)                                     &
         , cc_fld(npd_field, nclds)                                     &
         , conv_fld(npd_field, nclds)                                   &
         , dtau_s_fld(npd_field, nclds)                                 &
         , dtau_c_fld(npd_field, nclds)                                 &
         , skt_fld(npd_field)                                           &
         , emsfc_lw_fld(npd_field)                                      &
         , at_fld(npd_field, nclds)                                     &
         , dem_s_fld(npd_field, nclds)                                  &
         , dem_c_fld(npd_field, nclds)                                  &
         , fq_isccp_fld(npd_field, 7, 7)                                &
         , meanalbedocld(npd_field)                                     &
         , meantaucld(npd_field)                                        &
         , meanptop(npd_field)                                          &
         , totalcldarea(npd_field)

!     LOCAL VARIABLES.

      INTEGER                                                           &
           i                                                            &
!             LOOP VARIABLE
         , j                                                            &
!             LOOP VARIABLE
         , l
!             LOOP VARIABLE

!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL                                                              &
           d_mass(npd_profile, npd_layer)                               &
!             MASS THICKNESSES OF LAYERS
         , p(npd_profile, 0: npd_layer)                                 &
!             PRESSURE FIELD
         , t(npd_profile, 0: npd_layer)
!             TEMPERATURE FIELD

      REAL :: layer_heat_capacity(npd_profile, npd_layer)
!             Specific heat capacity of layer * d_mass

      REAL                                                              &
           dummy2d(1,1)
!     ADDITIONAL ARGUMENTS
      INTEGER                                                           &
           i_segment

      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'isccp'
      CHARACTER (LEN=240) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('ISCCP',zhook_in,zhook_handle)

      !Initialize error code
      ierr=0

      IF (is_ncol  <   10) THEN
       WRITE(0,*) 'WARNING: No. of ISCCP sampling columns less than 10'
       WRITE(0,*) 'This may introduce _systematic_ biases'
        IF (is_ncol  <   1) THEN
          cmessage =                                                    &
            'ERROR: Number of ISCCP sampling columns not specified'
          ierr=i_err_fatal
          GO TO 9999
        END IF
      END IF

      overlap=3
      top_height=1

!     SET THE THERMODYNAMIC PROPERTIES OF THE ATMOSPHERE.
! DEPENDS ON: r2_set_thermodynamic
      CALL r2_set_thermodynamic(nlit, nlevs, n_layer, nwet, list        &
         , .FALSE., .FALSE.                                             &
         , pstar                                                        &
         , p_layer_boundaries                                           &
         , p_layer_centres                                              &
         , dummy2d                                                      &
         , dummy2d, tac                                                 &
         , rho_r2, r_rho_levels, r_theta_levels                         &
         , q, qcl, qcf, qcf2, qrain, qgraup                             &
         , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
         , p(:,1:), t(:,1:), dummy2d, d_mass                            &
         , layer_heat_capacity                                          &
         , npd_field, npd_profile, npd_layer                            &
         )

      IF (lw_diag%l_isccp_cf .OR.                                       &
          lw_diag%l_isccp_cf_tau_0_to_p3 .OR.                           &
          lw_diag%l_isccp_cf_tau_p3_to_1p3 .OR.                         &
          lw_diag%l_isccp_cf_tau_1p3_to_3p6 .OR.                        &
          lw_diag%l_isccp_cf_tau_3p6_to_9p4 .OR.                        &
          lw_diag%l_isccp_cf_tau_9p4_to_23 .OR.                         &
          lw_diag%l_isccp_cf_tau_23_to_60 .OR.                          &
          lw_diag%l_isccp_cf_tau_ge_60) THEN
          IF (.NOT.lw_diag%l_isccp_weights) THEN
            cmessage =                                                  &
              '*** ERROR: ISCCP CLOUD FRACTIONS MAY BE DIAGNOSED ' //   &
              'ONLY IN CONJUNCTION WITH THE CORREPONDING WEIGHT.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF


!     Set up arrays for ISCCP_CLOUDTYPES_FLD

      DO i=1,nclds+1
        DO l=1,n_points
          phalf_fld(l,i)=p_layer_boundaries(l,(i-1))
        END DO
      END DO

       IF (l_3d_cca) THEN
         DO i=1, nclds
            DO l=1, n_points
               conv_fld(l,i)=cca(l,i)
            END DO
         END DO
       ELSE
         DO i=1, nclds
            DO l=1, n_points
               IF ( (i <= cct(l)-1).AND.(i >= ccb(l)) ) THEN
                  conv_fld(l,i)=cca(l,1)
               ELSE
                  conv_fld(l,i)=0.0e+00
               END IF
            END DO
         END DO
       END IF

      DO i=1,nclds
        DO l=1,n_points
          pfull_fld(l,i)=p_layer_boundaries(l,(i-1))
          qv_fld(l,i)=h2o(l,i)
          cc_fld(l,i)=conv_fld(l,i)+(1-conv_fld(l,i))*lca_area(l,i)
          dtau_s_fld(l,i)=0        ! lit points filled in below
          dtau_c_fld(l,i)=0        ! lit points filled in below
          skt_fld(l)=tstar(l)
          emsfc_lw_fld(l)=0.99
          at_fld(l,i)=tac(l,i)
          dem_s_fld(l,i)=0         ! lit points filled in below
          dem_c_fld(l,i)=0         ! lit points filled in below
        END DO
      END DO

      DO i=1,7
        DO j=1,7
          DO l=1,n_points
            fq_isccp_fld(l,i,j)=0
          END DO
        END DO
      END DO
      DO l=1,n_points
        meanalbedocld(l)=0.0
        meantaucld(l)=0.0
        meanptop(l)=0.0
        totalcldarea(l)=0.0
      END DO

!     Fill in optical thickness and emissivity for lit points.
!     1.0E-300 is used rather than 0.0 on the IF test to stop very
!     small weights from causing an overflow failure

      DO i=1,nclds
        DO l=1,nlit

          IF (sw_diag%ls_cloud_weight_extinction                        &
          (col_list(l), row_list(l),i) >  1.0e-300) THEN
            dtau_s_fld(list(l),i)=d_mass(l,nlevs+1-i)                   &
             *sw_diag%ls_cloud_extinction                               &
             (col_list(l), row_list(l),i)                               &
             /sw_diag%ls_cloud_weight_extinction                        &
             (col_list(l), row_list(l),i)
          END IF


          IF (sw_diag%cnv_cloud_weight_extinction                       &
          (col_list(l), row_list(l),i) >  1.0e-300) THEN
            dtau_c_fld(list(l),i)=d_mass(l,nlevs+1-i)                   &
             *sw_diag%cnv_cloud_extinction                              &
             (col_list(l), row_list(l),i)                               &
             /sw_diag%cnv_cloud_weight_extinction                       &
             (col_list(l), row_list(l),i)
          END IF

          IF (lw_diag%ls_cloud_weight_absorptivity                      &
          (col_list(l), row_list(l),i) >  1.0e-300) THEN
            dem_s_fld(list(l),i)=1.0-EXP(-1.666*d_mass(l,nlevs+1-i)     &
             *lw_diag%ls_cloud_absorptivity                             &
             (col_list(l), row_list(l),i)                               &
             /lw_diag%ls_cloud_weight_absorptivity                      &
             (col_list(l), row_list(l),i))
          END IF

          IF (lw_diag%cnv_cloud_weight_absorptivity                     &
          (col_list(l), row_list(l),i) >  1.0e-300) THEN
            dem_c_fld(list(l),i)=1.0-EXP(-1.666*d_mass(l,nlevs+1-i)     &
             *lw_diag%cnv_cloud_absorptivity                            &
             (col_list(l), row_list(l),i)                               &
             /lw_diag%cnv_cloud_weight_absorptivity                     &
             (col_list(l), row_list(l),i))
          END IF
        END DO
      END DO



! DEPENDS ON: isccp_cloudtypes_fld
      CALL isccp_cloudtypes_fld(i_segment                               &
!          Input
         , npd_field                                                    &
                        ! Number of points in whole field
         , n_points                                                     &
                        ! Number of points in segment
         , nlit                                                         &
                        ! Number of lit points
         , list                                                         &
                        ! Indexes of lit points
         , nclds                                                        &
                        ! Number of levels
         , is_ncol                                                      &
                        ! Number of columns for gridbox decomposition
         , pfull_fld                                                    &
                        ! Pressure on full levels
         , phalf_fld                                                    &
                        ! Pressure on full levels
         , qv_fld                                                       &
                        ! Water vapour
         , cc_fld                                                       &
                        ! Total cloud amount (cc+(1-cc)*ls)
         , conv_fld                                                     &
                        ! Convective cloud amount
         , dtau_s_fld                                                   &
                        ! Large-scale optical thickness
         , dtau_c_fld                                                   &
                        ! Convective optical thickness
         , top_height                                                   &
                        ! Flag to specify cloud top method
         , overlap                                                      &
                        ! Flag to specify overlap method
         , skt_fld                                                      &
                        ! Surface temperature
         , emsfc_lw_fld                                                 &
                        ! Surface emissivity
         , at_fld                                                       &
                        ! Atmospheric temperature
         , dem_s_fld                                                    &
                        ! Large-scale cloud emissivities
         , dem_c_fld                                                    &
                        ! Convective cloud emissivities
         , trindx_fld                                                   &
                        ! tropopause index
!          Output
         , fq_isccp_fld                                                 &
                        ! Cloud amounts in various ISCCP classes
         , meanalbedocld                                                &
                        ! Weighted mean cloud albedo
         , meantaucld                                                   &
                        ! Weighted mean cloud optical depth
         , meanptop                                                     &
                        ! Weighted mean cloud top pressure
         , totalcldarea                                                 &
                        ! Total cloud area
         )


!     ISCCP diagnostics

        DO i=1,7
          DO l=1,nlit

            IF (lw_diag%l_isccp_cf_tau_0_to_p3)                         &
                 lw_diag%isccp_cf_tau_0_to_p3                           &
                   (col_list(l),row_list(l),i)                          &
                   =fq_isccp_fld(list(l),1,i)

            IF (lw_diag%l_isccp_cf_tau_p3_to_1p3)                       &
                 lw_diag%isccp_cf_tau_p3_to_1p3                         &
                   (col_list(l),row_list(l),i)                          &
                   =fq_isccp_fld(list(l),2,i)

            IF (lw_diag%l_isccp_cf_tau_1p3_to_3p6)                      &
                 lw_diag%isccp_cf_tau_1p3_to_3p6                        &
                   (col_list(l),row_list(l),i)                          &
                   =fq_isccp_fld(list(l),3,i)

            IF (lw_diag%l_isccp_cf_tau_3p6_to_9p4)                      &
                 lw_diag%isccp_cf_tau_3p6_to_9p4                        &
                   (col_list(l),row_list(l),i)                          &
                   =fq_isccp_fld(list(l),4,i)

            IF (lw_diag%l_isccp_cf_tau_9p4_to_23)                       &
                 lw_diag%isccp_cf_tau_9p4_to_23                         &
                   (col_list(l),row_list(l),i)                          &
                   =fq_isccp_fld(list(l),5,i)

            IF (lw_diag%l_isccp_cf_tau_23_to_60)                        &
                 lw_diag%isccp_cf_tau_23_to_60                          &
                   (col_list(l),row_list(l),i)                          &
                   =fq_isccp_fld(list(l),6,i)

            IF (lw_diag%l_isccp_cf_tau_ge_60)                           &
                 lw_diag%isccp_cf_tau_ge_60                             &
                   (col_list(l),row_list(l),i)                          &
                   =fq_isccp_fld(list(l),7,i)

            IF (lw_diag%l_isccp_cf) THEN
              lw_diag%isccp_cf(col_list(l),row_list(l),i)=0.0
              DO j=1,7
                lw_diag%isccp_cf(col_list(l),row_list(l),i)=            &
                lw_diag%isccp_cf(col_list(l),row_list(l),i)+            &
                   fq_isccp_fld(list(l),j,i)
              END DO
            END IF
          END DO
        END DO

        IF (lw_diag%l_meanalbedocld .AND. lw_diag%l_totalcldarea) THEN
          DO l=1,nlit
            lw_diag%meanalbedocld(col_list(l), row_list(l))             &
                   =meanalbedocld(list(l)) * totalcldarea(list(l))
          END DO
        END IF

        IF (lw_diag%l_meantaucld .AND. lw_diag%l_totalcldarea) THEN
          DO l=1,nlit
            lw_diag%meantaucld(col_list(l), row_list(l))                &
                   =meantaucld(list(l)) * totalcldarea(list(l))
          END DO
        END IF

        IF (lw_diag%l_meanptop .AND. lw_diag%l_totalcldarea) THEN
          DO l=1,nlit
            lw_diag%meanptop(col_list(l), row_list(l))                  &
                   =meanptop(list(l)) * totalcldarea(list(l))
          END DO
        END IF

        IF (lw_diag%l_totalcldarea) THEN
          DO l=1,nlit
            lw_diag%totalcldarea(col_list(l), row_list(l))              &
                   =totalcldarea(list(l))
          END DO
        END IF

        IF (lw_diag%l_isccp_weights) THEN
          DO l=1,nlit
            lw_diag%isccp_weights(col_list(l), row_list(l))=1.0
          END DO
        END IF

 9999 CONTINUE
! Check error condition
      IF (ierr /= i_normal) THEN

        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook('ISCCP',zhook_out,zhook_handle)
      END SUBROUTINE isccp
