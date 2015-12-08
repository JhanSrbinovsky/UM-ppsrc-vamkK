! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to set thermodynamic properties

! Purpose:
!   Pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.

! Method:
!   Straightforward.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
      SUBROUTINE r2_set_thermodynamic(n_profile, nlevs, n_layer,q_levels&
         , i_gather, l_extra_top, l_boundary_temperature                &
         , pstar                                                        &
         , p_layer_boundaries                                           &
         , p_layer_centres                                              &
         , height_theta                                                 &
         , height_rho                                                   &
         , tac                                                          &
           ! IN. Height and moisture information for the calculation
           ! of layer masses
         , rho_r2, r_rho_levels, r_theta_levels                         &
         , q, qcl, qcf, qcf2, qrain, qgraup                             &
         , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
           ! OUT
         , p, t, t_bdy, d_mass                                          &
         , layer_heat_capacity                                          &
         , npd_field, npd_profile, npd_layer                            &
         )


      USE earth_constants_mod, ONLY: g
      USE atmos_constants_mod, ONLY: cp

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE


!     INCLUDED COMDECKS
! C_PERMA start

      ! Specific heat capacity of water vapour (J/kg/K)
      REAL,PARAMETER:: HCAPV=1850.0

      ! Specific heat capacity of water (J/kg/K)
      REAL,PARAMETER:: HCAPW=4180.0

      ! Specific heat capacity of ice (J/kg/K)
      REAL,PARAMETER:: HCAPI=2100.0

      ! Density of ice (kg/m3)
      REAL,PARAMETER:: RHO_ICE=917

      ! Rate of change of ice potential with temperature
      ! RHO_ICE*LF/ZERODEGC*1/(RHO_WATER*G) (m/K)
      REAL,PARAMETER:: DPSIDT=114.3

! C_PERMA end

!     DUMMY ARGUMENTS.
!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
           npd_field                                                    &
!             SIZE OF ARRAY FROM UM
         , npd_profile                                                  &
!             MAXIMUM NUMBER OF PROFILES
         , npd_layer
!             MAXIMUM NUMBER OF LAYERS

!     SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
           n_profile                                                    &
!             NUMBER OF PROFILES
         , nlevs                                                        &
!             Number of levels in the main model
         , q_levels                                                     &
                     ! Number of moist model levels
         , n_layer                                                      &
!             Number of layers seen by radiation
!             NUMBER OF LEVELS
         , i_gather(npd_field)
!             LIST OF POINTS TO GATHER

      REAL                                                              &
                !, INTENT(IN)
           pstar(npd_field)                                             &
!             SURFACE PRESSURES
         ,p_layer_boundaries(npd_field,0:nlevs)                         &
!             PRESSURE AT EDGES OF LAYERS
         ,p_layer_centres(npd_field,0:nlevs)                            &
!             PRESSURE AT CENTRES OF LAYERS
         ,height_theta(npd_field,0:nlevs)                               &
!             HEIGHT AT CENTRES OF LAYERS
         ,height_rho(npd_field,nlevs)                                   &
!             HEIGHT AT EDGES OF LAYERS
         , tac(npd_field, nlevs)
!             TEMPERATURES AT CENTRES OF LAYERS
       REAL, INTENT(IN) ::                                              &
        rho_r2(npd_field,nlevs)                                         &
                                ! Air density*radius of earth**2 / kg m-1
      , r_rho_levels(npd_field,nlevs)                                   &
                                      ! Height of rho levels / m
      , r_theta_levels(npd_field,0:nlevs)                               &
                                           ! Height of theta levels / m
      , q(npd_field,q_levels)                                           &
                                      ! Vapour content / kg kg-1
      , qcl(npd_field,q_levels)                                         &
                                      ! Liquid water content / kg kg-1
      , qcf(npd_field,q_levels)                                         &
                                      ! Ice content / kg kg-1
      , qcf2(npd_field,q_levels)                                        &
                                      ! Second ice content / kg kg-1
      , qrain(npd_field,q_levels)                                       &
                                      ! Rain water content / kg kg-1
      , qgraup(npd_field,q_levels)    ! Graupel content / kg kg-1

       LOGICAL, INTENT(IN) ::                                           &
        l_mcr_qcf2                                                      &
                                      ! Use second prognostic ice
      , l_mcr_qrain                                                     &
                                      ! Use prognostic rain
      , l_mcr_qgraup                                                    &
                                      ! Use prognostic graupel
      , l_mixing_ratio                ! Use mixing ratios

      LOGICAL                                                           &
                !, INTENT(IN)
           l_extra_top                                                  &
!             Flag to create an extra top layer for radiation
         , l_boundary_temperature
!             FLAG TO CALCULATE TEMPERATURES AT BOUNADRIES OF LAYERS.


      REAL                                                              &
                !, INTENT(OUT)
           d_mass(npd_profile, npd_layer)                               &
!             MASS THICKNESSES OF LAYERS
         , p(npd_profile, npd_layer)                                    &
!             PRESSURE FIELD
         , t(npd_profile, npd_layer)                                    &
!             TEMPERATURE FIELD
         , t_bdy(npd_profile, 0: npd_layer)
!             TEMPERATURES AT EDGES OF LAYERS

      REAL, INTENT(OUT) :: layer_heat_capacity(npd_profile, npd_layer)
!             Specific heat capacity of layer * d_mass

!     LOCAL VARIABLES.

      INTEGER                                                           &
           i                                                            &
!             LOOP VARIABLE
         , ii                                                           &
!             LOOP VARIABLE
         , l                                                            &
!             LOOP VARIABLE
         , lg                                                           &
!             INDEX TO GATHER
         , i_top_copy                                                   &
!             Topmost layer where properties are set by copying the
!             input fields.
         , level_number  ! Model level in the upward counting
                         ! coordinate system

      REAL                                                              &
           pu                                                           &
!             PRESSURE FOR UPPER LAYER
         , pl                                                           &
!             PRESSURE FOR LOWER LAYER
         , pml1                                                         &
!             PRESSURE FOR INTERPOLATION
         , wtl                                                          &
!             WEIGHT FOR LOWER LAYER
         , wtu                                                          &
!             WEIGHT FOR UPPER LAYER
         , rho1(npd_profile)                                            &
                                         ! Air density / kg m-3
         , rho2(npd_profile)                                            &
                                         ! Air density / kg m-3
         , deltaz(npd_profile, npd_layer)                               &
                                         ! Layer thickness / m
         , rhodz_moist(npd_profile,npd_layer)                           &
                                              ! Thick*density / kg m-2
         , q_total(npd_profile)          ! Total moisture mr / kg kg-1

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('R2_SET_THERMODYNAMIC',zhook_in,zhook_handle)

!     Calculate properties at the centres of layers: a top layer is
!     artificially created, if required, reaching from the top of
!     the model's actual atmosphere to zero pressure, assuming an
!     isothermal profile.
      IF (l_extra_top) THEN
!       Note that we will continue to use this calculation for the
!       artificial top layer when the more accurate mixing ratio code
!       is used for the lower levels because we do not have the model
!       information to make a better estimate.
        DO l=1, n_profile
          lg=i_gather(l)
          p(l, 1)=0.5e+00*p_layer_boundaries(lg, nlevs)
          t(l, 1)=tac(lg, nlevs)
          d_mass(l, 1)=p_layer_boundaries(lg, nlevs)/g
          layer_heat_capacity(l, 1)=cp
        END DO
!       The second radiative layer will be the first to have properties
!       set by copying input fields.
        i_top_copy=2
      ELSE
!       The second radiative layer will be the first to have properties
!       set by copying input fields.
        i_top_copy=1
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          p(l, i)=p_layer_centres(lg, n_layer+1-i)
          t(l, i)=tac(lg, n_layer+1-i)
          IF (.NOT. l_mixing_ratio) THEN
            d_mass(l, i)                                                &
              =ABS(p_layer_boundaries(lg, n_layer-i)                    &
              -p_layer_boundaries(lg, n_layer+1-i))/g
          END IF
        END DO
      END DO
      IF (l_mixing_ratio) THEN
        ! Perform more accurate calculation of the mass of the
        ! layer

        DO i=i_top_copy, n_layer      ! This is the level in the
                                      ! downward counting coord. system

          level_number = n_layer+1-i  ! This is the level in the
                                      ! upward counting coord. system

          ! We note that rho_r2, r_rho_levels, q, qcl, qcf, qrain
          ! qcf2 and qgraup are in the upward
          ! counting label system and defined on all points. rho1,
          ! rho2, deltaz, rhodz_moist, rhodz_dry and q_total are in the
          ! downward counting label system and defined on gathered
          ! points.

          DO l=1, n_profile
            lg=i_gather(l)

            ! Calculate densities at the boundaries of the layer
            ! by removing the r**2 term from rho_r2
            ! Rho1 is the density at the lower boundary
            rho1(l)= rho_r2(lg,level_number)                            &
                             /( r_rho_levels(lg,level_number) *         &
                                r_rho_levels(lg,level_number) )


            ! Check whether there is a rho level above the current
            ! moist level.
            IF (level_number  <   nlevs) THEN
              ! Rho2 is the density at the upper boundary.
              rho2(l)= rho_r2(lg,level_number+1)                        &
                               /( r_rho_levels(lg,level_number+1) *     &
                                  r_rho_levels(lg,level_number+1) )

              ! Calculate the average value of rho across the layer
              ! multiplied by the layer thickness and the layer
              ! thickness.
              rhodz_moist(l,i) =                                        &
                          rho2(l) * ( r_theta_levels(lg,level_number)   &
                                    - r_rho_levels(lg,level_number) )   &
                       +  rho1(l) * ( r_rho_levels(lg,level_number+1)   &
                                    - r_theta_levels(lg,level_number))
              deltaz(l,i) = r_rho_levels(lg,level_number+1)             &
                                - r_rho_levels(lg,level_number)

              IF (level_number  ==  1) THEN
                ! For the lowest layer we need to extend the lower
                ! boundary from the lowest (in height) rho level
                ! to the surface. The surface is the 0'th theta level
                ! when counting upwards.
                deltaz(l,i) = r_rho_levels(lg,2)                        &
                                - r_theta_levels(lg,0)
                rhodz_moist(l,i) = rhodz_moist(l,i)*deltaz(l,i)         &
                       / (r_rho_levels(lg,2)-r_rho_levels(lg,1))
              END IF  ! level_number  ==  1

            ELSE
              ! For a top layer higher than the highest rho level
              ! we can calculate a pseudo rho level. We will assume
              ! it has a similar density to the rho level below
              ! and that the intervening theta level is in the centre
              ! of the layer.
              deltaz(l,i) = 2.0*(r_theta_levels(lg,level_number)        &
                                  -r_rho_levels(lg,level_number))
              rhodz_moist(l,i) = rho1(l) * deltaz(l,i)

            END IF  ! level_number  <   nlevs

            IF (level_number  <=  q_levels) THEN
              ! Calculate total moisture
              q_total(l) = q(lg,level_number) + qcl(lg,level_number)    &
                                              + qcf(lg,level_number)

              ! Calculate the specific heat capacity of the moist air
              layer_heat_capacity(l,i)=cp + q(lg,level_number)*hcapv +  &
                qcl(lg,level_number)*hcapw + qcf(lg,level_number)*hcapi

              ! Add on contributions from optional prognostics
              IF (l_mcr_qcf2) THEN
                q_total(l) = q_total(l) + qcf2(lg,level_number)
                layer_heat_capacity(l,i) = layer_heat_capacity(l,i) +   &
                  qcf2(lg,level_number)*hcapi
              END IF  ! l_mcr_qcf2

              IF (l_mcr_qrain) THEN
                q_total(l) = q_total(l) + qrain(lg,level_number)
                layer_heat_capacity(l,i) = layer_heat_capacity(l,i) +   &
                  qrain(lg,level_number)*hcapw
              END IF  ! l_mcr_qrain

              IF (l_mcr_qgraup) THEN
                q_total(l) = q_total(l) + qgraup(lg,level_number)
                layer_heat_capacity(l,i) = layer_heat_capacity(l,i) +   &
                  qgraup(lg,level_number)*hcapi
              END IF  ! l_mcr_qgraup

            ELSE

              ! There is no moisture in the model
              q_total(l) = 0.0
              layer_heat_capacity(l,i) = cp

            END IF  ! level_number le q_levels

            ! Convert from the moist density of air to the dry
            ! density of air
            d_mass(l,i) = rhodz_moist(l,i) / (1.0 + q_total(l))

          END DO  ! l
        END DO  ! i

        ! If l_mixing ratio is true, layer_heat_capacity is the true
        ! heat capacity of the moist air including condensate.
        DO i=1,n_layer
          DO l=1,n_profile
            layer_heat_capacity(l,i)=d_mass(l,i)*                       &
                                       layer_heat_capacity(l,i)
          END DO
        END DO

      ELSE

        ! If l_mixing_ratio is false, layer_heat_capacity is based
        ! on the specific heat of dry air and a mass derived from the
        ! hydrostatic approximation.
        DO i=1,n_layer
          DO l=1,n_profile
            layer_heat_capacity(l,i)=d_mass(l,i)*cp
          END DO
        END DO

      END IF  ! L_mixing_ratio


      IF (l_boundary_temperature) THEN

!        INTERPOLATE TEMPERATURES AT THE BOUNDARIES OF LAYERS
!        FROM THE EXNER FUNCTION.
         DO l=1, n_profile
            lg=i_gather(l)

!           TAKE THE TEMPERATURE OF THE AIR JUST ABOVE THE SURFACE AS
!           THE TEMPERATURE AT THE MIDDLE OF THE BOTTOM LAYER.
            t_bdy(l, n_layer)=tac(lg, 1)
!        Set the temperature at the top of the model by isothermal
!        extrapolation of the temperature in the middle of the topmost
!        layer used in the full atmospheric model.
!        otherwise.
         t_bdy(l, 0)=tac(lg, nlevs)

         END DO

!        If an extra layer is to be added, the temperature at the
!        bottom of this layer must be set by isothermal extrapolation.
         IF (l_extra_top) THEN
           DO l=1, n_profile
             t_bdy(l, 1)=tac(i_gather(l), nlevs)
           END DO
         END IF


!        Now set the temperatures at interior boundaries of the
!        profile supplied using interpolation with the EXNER function.
         DO i=i_top_copy, n_layer-1
            ii=n_layer-i
            DO l=1, n_profile
               lg=i_gather(l)
               wtu=height_rho(lg,ii+1)-height_theta(lg,ii)
               wtl=height_theta(lg,ii+1)-height_rho(lg,ii+1)
               t_bdy(l, i)=(wtu*tac(lg, n_layer+1-i)                    &
                  +wtl*tac(lg, n_layer-i))/(wtl+wtu)
            END DO
         END DO

      END IF


      IF (lhook) CALL dr_hook('R2_SET_THERMODYNAMIC',zhook_out,zhook_handle)
      END SUBROUTINE r2_set_thermodynamic
