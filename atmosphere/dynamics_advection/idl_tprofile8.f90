! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine IDL_tprofile8

      SUBROUTINE IDL_tprofile8(                                         &
                            row_length, rows, model_levels,             &
                            me, halo_i, halo_j,                         &
                            r_theta_levels, r_rho_levels,               &
                            eta_theta_levels, eta_rho_levels,           &
                            theta, exner_rho_levels, exner_theta_levels,&
                            theta_ref, exner_ref,                       &
                            height_domain,                              &
                            p_surface, theta_surface)

! Purpose:
!          Sets up initial data for idealised problems.
!          Initial isothermal temperature profile at all points
!          is set to the same everywhere
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE earth_constants_mod, ONLY: g, earth_radius
      USE atmos_constants_mod, ONLY:                                    &
          r, cp, kappa, p_zero 

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER                                                           &
        row_length                                                      &
                         ! number of points on a row
      , rows                                                            &
                         ! number of rows in a theta field
      , model_levels                                                    &
                         ! number of model levels
      , halo_i                                                          &
                             ! Size of halo in i direction.
      , halo_j               ! Size of halo in j direction.

      REAL                                                              &
        theta_surface                                                   &
      , p_surface                                                       &
      , height_domain

      INTEGER                                                           &
        me         ! My processor number



      REAL                                                              &
           ! vertical co-ordinate information
        r_theta_levels(1-halo_i:row_length+halo_i,                      &
                       1-halo_j:rows+halo_j,0:model_levels)             &
      , r_rho_levels(1-halo_i:row_length+halo_i,                        &
                     1-halo_j:rows+halo_j, model_levels)                &
      , eta_theta_levels(0:model_levels)                                &
      , eta_rho_levels(model_levels)

! Output Arrays from this routine
      REAL                                                              &
        theta(1-halo_i:row_length+halo_i                                &
      ,       1-halo_j:rows+halo_j, model_levels)                       &
      , exner_rho_levels(1-halo_i:row_length+halo_i                     &
      ,                  1-halo_j:rows+halo_j, model_levels+1)          &
      , exner_theta_levels(1-halo_i:row_length+halo_i                   &
      ,                    1-halo_j:rows+halo_j, model_levels)          &
      , theta_ref(model_levels)                                         &
                                 !theta profile for use in sponge & lbcs
      , exner_ref(model_levels + 1)  ! Exner profile for use in lbcs

! local variables
      INTEGER                                                           &
        i, j, k

      REAL                                                              &
        weight                                                          &
      , temp                                                            &
      , goverCpT                                                        &
      , delta_z                                                         &
      , exner_surface
     
      
     REAL   :: exner_ref_pro(0:model_levels+1) 
     REAL   :: thetav_ref_pro(0:model_levels)
     REAL   :: rho_ref_pro(model_levels)
     REAL   :: T0_ref(0:model_levels)
     REAL   :: geopotential_factor(0:model_levels)
     REAL   :: levs(0:model_levels)
     REAL   :: u_term(1:model_levels)
     REAL   :: intw_rho2w(model_levels,2)
     REAL   :: intw_w2rho(model_levels,2)
     
     REAL :: grad
      INTEGER :: prof_type

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------
! Section 1  ISOTHERMAL   temperature constant
!        (constant) temperature =  theta_surface
!-------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_TPROFILE8',zhook_in,zhook_handle)
!
! Set initial guess
!     
      exner_ref_pro(0)  = (p_surface / p_zero)**kappa
      exner_ref_pro(:)  = exner_ref_pro(0)
      thetav_ref_pro(:) = theta_surface
      rho_ref_pro(:)    = 1.0
      T0_ref(:)         = theta_surface*exner_ref_pro(0)
      geopotential_factor(:) = g
! U_term not needed
      u_term(:) = 0
      
      DO k=1, model_levels
        intw_w2rho(k,1) = ( eta_rho_levels(k)-eta_theta_levels(k-1) ) / &
                          ( eta_theta_levels(k)-eta_theta_levels(k-1) )
        intw_w2rho(k,2) = 1.0 - intw_w2rho(k,1)
      END DO

      DO k=1, model_levels-1
        intw_rho2w(k,1) = ( eta_theta_levels(k)-eta_rho_levels(k) ) /   &
                          ( eta_rho_levels(k+1)-eta_rho_levels(k) )
        intw_rho2w(k,2) = 1.0 - intw_rho2w(k,1)
      END DO

      intw_rho2w(model_levels,1) = 0.5
      intw_rho2w(model_levels,2) = 0.5  ! as in New dynamics 
      
      prof_type = 2     ! T specified
      T0_ref(0) = theta_surface*exner_ref_pro(0)
      
      DO j = 1-halo_j, rows+halo_j
        DO i = 1-halo_i, row_length+halo_i 
          T0_ref(0) = theta_surface*exner_ref_pro(0)                          &
                      -6.0/1000.0*(r_theta_levels(i,j,0)-earth_radius)
        
          levs(0) = r_theta_levels(i,j,0) - earth_radius
          DO k = 1,model_levels
            grad = -6.0/1000.0
            IF(  r_rho_levels(i,j,k) > 13000.0 ) grad = -2.0/1000.0
            IF(  r_rho_levels(i,j,k) > 15000.0 ) grad = 0.0
            IF(  r_rho_levels(i,j,k) > 16000.0 ) grad = 5.0/8000.0
            IF(  r_rho_levels(i,j,k) > 25000.0 ) grad = 0.0
            T0_ref(k) = T0_ref(k-1) + grad * (r_theta_levels(i,j,k)    &
                                            - r_theta_levels(i,j,k-1))
                                            
            levs(k) = r_rho_levels(i,j,k) - earth_radius
          END DO
! DEPENDS ON: eg_Newton
          CALL eg_newton(exner_ref_pro, thetav_ref_pro, rho_ref_pro, T0_ref, &
                 geopotential_factor, levs, intw_rho2w, intw_w2rho,          &
                 p_surface, model_levels, prof_type, u_term)
                 
! Copy data back
          DO k=1,model_levels
            theta(i,j,k) = thetav_ref_pro(k)
            exner_rho_levels(i,j,k) = exner_ref_pro(k)
          END DO
          exner_rho_levels(i,j,model_levels+1) = exner_ref_pro(model_levels+1)
        END DO
      END DO


      DO k = 1, model_levels - 1
        DO j = 1-halo_j, rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            weight = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/   &
                      (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
            exner_theta_levels(i,j,k) =   weight *                      &
                                         exner_rho_levels(i,j,k) +      &
                        (1.0 - weight) * exner_rho_levels(i,j,k+1)
          END DO
        END DO
      END DO

      k = model_levels
      DO j = 1-halo_j, rows+halo_j
        DO i = 1-halo_i, row_length+halo_i
          exner_theta_levels(i,j,k) =  0.5 * exner_rho_levels(i,j,k) +  &
                                       0.5 * exner_rho_levels(i,j,k+1)
        END DO
      END DO

!  exner_ref required to set lbcs
      DO k = 1, model_levels+1
        exner_ref(k) = exner_rho_levels(1,1,k)
        theta_ref(k) = theta(1,1,k)
      END DO
      exner_ref(model_levels+1) = exner_rho_levels(1,1,model_levels+1)

      IF (lhook) CALL dr_hook('IDL_TPROFILE8',zhook_out,zhook_handle)

      END SUBROUTINE IDL_tprofile8
!  End subroutine IDL_tprofile8
