! *****************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*********************************

!+ Balance Exner in LBCs and set rho using equation of state.
!  Reset w=0
      MODULE eg_balance_lbc_values_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_balance_lbc_values(                                 &
        EXNER_LBC, RHO_LBC, THETA_LBC, Q_LBC, W_LBC, W_ADV_LBC,         &
        U_LBC, V_LBC,                                                   &
        QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,              &
        L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &  
        R_RHO_LEVELS, R_THETA_LEVELS,                                   &
        ROW_LENGTH, ROWS, model_levels, wet_levels,                     &
        HALO_I, HALO_J, LENRIM,                                         &
        LENRIMU, LENRIMV,                                               &
        LBC_START, LBC_START_U, LBC_START_V,                            &
        RIMWIDTH, N_RIMS_TO_DO, RIMWEIGHTS, AT_EXTREMITY,               &
        DELTA_PHI, DELTA_LAMBDA,                                        &
        BASE_PHI, BASE_LAMBDA,                                          &
        DATASTART, LAT_ROT_NP,                                          &  
        GLOBAL_ROW_LENGTH, GLOBAL_ROWS)

USE atmos_constants_mod, ONLY: cp, kappa, p_zero, r, c_virtual

USE earth_constants_mod, ONLY : g, two_omega

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE metric_terms_mod
USE gravity_mod
USE eg_alpha_mod,      ONLY : alpha_w
USE timestep_mod,      ONLY : timestep
USE level_heights_mod, ONLY : eta_rho_levels, eta_theta_levels
USE coriolis_mod,      ONLY : f1_comp, f2_comp
USE horiz_grid_mod,    ONLY : xi1_u, xi2_v, intw_p2u,intw_p2v,          &
                              intw_w2rho,intw_rho2w 
USE atm_fields_bounds_mod
USE UM_ParParams
USE eg_parameters_mod, ONLY : l_cartesian
USE lam_lbc_weights_mod

      IMPLICIT NONE

! routine of the same name with changes to take into
! account the different EG disretisation of the vertical 
! momentum equation - Coriolis terms and non-constant g(z)
!
! Since EG coriolis discretisation has a weak dependence on 
! rho, there is a small inconsistency in the rho used to compute
! exner and the subsequent resetting of rho to satisfy the 
! equation of state.
!
! Exner: reset in the LBC region of a LAM model
! such that the vertical pressure gradient balances
! gravity, the vertical coriolis and metric terms.
! This should substantially reduce the vertical
! accelerations in the LBC region, thereby reducing
! problems occuring when vertical levels are different in
! the driving and driven models.
!
!  NB. vertical accelerations resultsing from parameterised
!      processes are not considered here.
!
!        2    2
!  Dw = u  + v                              d Exner    w
!  __   _______ +  f u - f v - g - C theta  _______ + S
!                   2     1         p     v
!  Dt      r                                  dr
!
!  where f1 = 2 Omega sin(lambda) cos(phi )
!                                        0
!
!        f2 = 2 Omega ( cos(phi)sin(phi ) - sin(phi)cos(lambda)cos(phi )
!                                      0                              0
!
! Rho: reset using the equation of state.
!
! W wind: set to zero.
!
! The N_RIMS_TO_DO variable allows only a certain depth of
! rim width to be done.
! The RIMWEIGHTS array contains a weighting factor that
! controls what proportion of FIELD is used and what
! proportion of LBC. In the halo area only LBC data is used.
!
! If RIMWIDTH is zero then the routine assumes there are no
! LBCs available and exits without updating any fields
!
! The logicals L_DO_BOUNDARIES and L_DO_HALOS have both
! been hard-wired to .TRUE.
!
! Code structure:
!       Section 1. Setup constants
!       Section 2. Setup blending weight arrays
!       Section 3. Calculate coriolis terms
!       Section 4. Impose boundary conditions for each boundary in turn
!                       -reset exner
!                       -reset rho
!                       -reset w
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input
!
! Arguments:
! Parameters required for dimensioning some of the arguments

      LOGICAL, INTENT(IN) ::                                            &
       AT_EXTREMITY(4)    ! Indicates if this processor is at
                          ! the edge (North,East,South,West)
                          ! of the processor grid

      LOGICAL, INTENT(IN) :: L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup

      INTEGER, INTENT(IN) ::                                            &       
         ROW_LENGTH    &  ! number of points in row             
       , ROWS          &  ! number of rows                      
       , HALO_I        &  ! size of FIELD halo in EW direction
       , HALO_J        &  ! size of FIELD halo in NS direction
       , MODEL_LEVELS  &  ! number of vertical levels
       , WET_LEVELS    &  ! number of moist levels
       , LENRIM        &  ! size of one level of LBC exner data
       , LENRIMU       &  ! size of one level of LBC u data
       , LENRIMV       &  ! size of one level of LBC v data
       , LBC_START(4)  &  ! offset of each side of exner LBC data
       , LBC_START_U(4)&  ! offset of each side of u LBC data
       , LBC_START_V(4)&  ! offset of each side of v LBC data
       , RIMWIDTH      &  ! size (width) of the boundary area
       , N_RIMS_TO_DO  &  ! number of rims to do (counting from
                          !      the outside in)
       , DATASTART(3)                                                   &  
               ! index of first grid point held by this processor
       , GLOBAL_ROW_LENGTH &! no. of columns across LAM domain
       , GLOBAL_ROWS      ! no. of rows in the LAM domain

      REAL, INTENT(IN) ::                                               &  
         delta_phi      &  ! grid length in latitude  (radians)
       , delta_lambda   &  ! grid length in longitude (radians)
       , base_phi       &  ! latitude of first grid pt (radians)
       , base_lambda    &  ! longitude of first grid pt (radians)
       , lat_rot_NP     &  ! lat of rotated pole
       , rimweights(rimwidth) &  !Weight to apply to each successive rim

        ! vertical co-ordinate information
       , r_theta_levels(1-halo_i:row_length+halo_i,                     &
                     1-halo_j:rows+halo_j,0:model_levels)               &
       , r_rho_levels(1-halo_i:row_length+halo_i,                       &
                   1-halo_j:rows+halo_j, model_levels)                  &

        ! LBC data arrays 
       , THETA_LBC(LENRIM, 0:MODEL_LEVELS)                              & 
       , Q_LBC(LENRIM, 0:MODEL_LEVELS)                                  & 
       , QCL_LBC(LENRIM, 0:MODEL_LEVELS)                                & 
       , QCF_LBC(LENRIM, 0:MODEL_LEVELS)                                & 
       , QCF2_LBC(LENRIM, 0:MODEL_LEVELS)                               & 
       , QRAIN_LBC(LENRIM, 0:MODEL_LEVELS)                              & 
       , QGRAUP_LBC(LENRIM, 0:MODEL_LEVELS)                             & 
       , U_LBC(LENRIMU, MODEL_LEVELS)                                   & 
       , V_LBC(LENRIMV, MODEL_LEVELS)

      REAL, INTENT(INOUT) ::                                            &
         EXNER_LBC(LENRIM, MODEL_LEVELS+1)                              &
       , RHO_LBC(LENRIM, MODEL_LEVELS)                                  &
       , W_LBC(LENRIM, 0:MODEL_LEVELS)                                  &
       , W_ADV_LBC(LENRIM, 0:MODEL_LEVELS) 

! Local variables
      INTEGER :: I          ! loop counters 
      INTEGER :: J
      INTEGER :: K 
      INTEGER :: kp
      INTEGER :: gj         ! j position in global field
      INTEGER :: gi         ! i position in global field
      INTEGER :: weights_i  ! modified I to point to ew_weights array
      INTEGER :: weights_j  ! modified J to point to ns_weights array
      INTEGER :: LBCarea    ! marker to designate which LBC region
                            !   a point lies within
      INTEGER :: row_start_pt      ! first point along row to update
      INTEGER :: row_end_pt        ! last point along row to update
      INTEGER :: first_row         ! first row to update
      INTEGER :: last_row          ! last row to update
      INTEGER :: lbc_row_len       ! Length of a row of LBC data
      INTEGER :: lbc_row_len_u     ! Length of a row of LBC u data
      INTEGER :: first_pt_of_lbc   
                            ! First point on row contained in LBC data
      INTEGER :: first_row_of_lbc  ! First row contained in LBC data
      INTEGER :: LBC_address(1-halo_i:row_length+halo_i,                &
                    1-halo_j:rows+halo_j)
                            ! address in Exner LBC array at i,j
      INTEGER :: LBC_address_jm1(1-halo_i:row_length+halo_i,            &
                   1-halo_j:rows+halo_j)
                            ! address in V LBC array at i,j-1
      INTEGER :: LBC_address_im1(1-halo_i:row_length+halo_i,            &
                   1-halo_j:rows+halo_j)
                            ! address in U LBC array at i-1,j
      INTEGER :: LBC_address_u(1-halo_i:row_length+halo_i,              &
                   1-halo_j:rows+halo_j)
                            ! address in U LBC field at i,j
      INTEGER :: LBC_address_v(1-halo_i:row_length+halo_i,              &
                   1-halo_j:rows+halo_j)
                            ! address in V LBC field at i,j
      INTEGER :: LBC_address_u_jm1(1-halo_i:row_length+halo_i,          &
                   1-halo_j:rows+halo_j) 
      INTEGER :: LBC_address_u_im1(1-halo_i:row_length+halo_i,          &
                   1-halo_j:rows+halo_j)
      INTEGER :: LBC_address_ip1(1-halo_i:row_length+halo_i,            &
                   1-halo_j:rows+halo_j)
      INTEGER :: LBC_address_jp1(1-halo_i:row_length+halo_i,            &
                   1-halo_j:rows+halo_j)

      REAL :: temp, denom1, dpdz, u_av, v_av, metric_term 
      REAL :: theta_wet_lbc(lenrim, 0:model_levels),                     &
              dry2wet(lenrim, 0:model_levels)

      REAL :: work1(1-halo_i:row_length+halo_i,                          &
                    1-halo_j:rows+halo_j,model_levels)

      REAL :: f1_comp_eh(pdims_l%i_start:pdims_l%i_end,                  &
                         pdims_l%j_start:pdims_l%j_end),                 &
              f2_comp_eh(pdims_l%i_start:pdims_l%i_end,                  &
                         pdims_l%j_start:pdims_l%j_end),                 &
              g_theta_eh(tdims_l%i_start:tdims_l%i_end,                  &
                         tdims_l%j_start:tdims_l%j_end,                  &
                         tdims_l%k_start:tdims_l%k_end)

      INTEGER :: iloc, ii, jj, ext_row_len

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------
!
! The following diagram breaks up a LAM model area into a number
! of subcomponents (the letters will be referred to in the code):
! (Assumes the following model sizes for this example:
!  ROW_LENGTH=ROWS=12
!  RIMWIDTH=3
!  HALO_I=HALO_J=2)
!
!       North
!  aaaaaaaaaaaaaaaa
!  aaaaaaaaaaaaaaaa
!  bbcccddddddeeeff
!  bbcccddddddeeeff
!  bbcccddddddeeeff
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  llmmmnnnnnnooopp
!  llmmmnnnnnnooopp
!  llmmmnnnnnnooopp
!  qqqqqqqqqqqqqqqq
!  qqqqqqqqqqqqqqqq
!       South

!========================================================================
! 1. Set up local constant for later
!========================================================================

      IF (lhook) CALL dr_hook('EG_BALANCE_LBC_VALUES',zhook_in,zhook_handle)

      ext_row_len = row_length + 2*halo_i

! Calculate wet theta in lbc's

      IF( at_extremity(Psouth) .OR. at_extremity(Pnorth) .OR.               &
          at_extremity(Peast)  .OR. at_extremity(Pwest) ) THEN

         DO k = tdims_s%k_start, tdims_s%k_end
            DO i = 1, lenrim
               dry2wet(i,k) = 1.0 + q_lbc(i,k) + qcl_lbc(i,k) + qcf_lbc(i,k)
            END DO
         END DO

         IF( l_mcr_qcf2 ) THEN
            DO k = tdims_s%k_start, tdims_s%k_end
               DO i = 1, lenrim
                  dry2wet(i,k) = dry2wet(i,k) + qcf2_lbc(i,k) 
               END DO
            END DO
         END IF

         IF( l_mcr_qgraup ) THEN
            DO k = tdims_s%k_start, tdims_s%k_end
               DO i = 1, lenrim
                  dry2wet(i,k) = dry2wet(i,k) + qgraup_lbc(i,k) 
               END DO
            END DO
         END IF

         IF( l_mcr_qrain ) THEN
            DO k = tdims_s%k_start, tdims_s%k_end
               DO i = 1, lenrim
                  dry2wet(i,k) = dry2wet(i,k) + qrain_lbc(i,k) 
               END DO
            END DO
         END IF

         DO k = tdims_s%k_start, tdims_s%k_end
            DO i = 1, lenrim
               theta_wet_lbc(i,k) = theta_lbc(i,k)/dry2wet(i,k)
            END DO
         END DO

      END IF

! Initialise rho = constant, this will be overwritten after exner_lbc 
! is computed

      rho_lbc = -1.0
 
! Extend f1 and f2 arrays to have large halos

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          f1_comp_eh(i,j) = f1_comp(i,j)
          f2_comp_eh(i,j) = f2_comp(i,j)
          g_theta_eh(i,j,:) =  g_theta(i,j,:)
        END DO
      END DO 

! DEPENDS ON: swap_bounds
      CALL swap_bounds(f1_comp_eh,row_length,rows,1,halo_i,halo_j,        &
                       fld_type_p,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL swap_bounds(f2_comp_eh,row_length,rows,1,halo_i,halo_j,        &
                       fld_type_p,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL swap_bounds(g_theta_eh,row_length,rows,model_levels+1,         &
                       halo_i,halo_j,fld_type_p,.FALSE.) 

!========================================================================
! 4.  Now apply the boundary conditions to exner, rho and w
!========================================================================

!------------------------------------------------------------------------
! 4.1 Southern LBC region
!------------------------------------------------------------------------

      IF (at_extremity(PSouth)) THEN

!----------------------------------------------------------------
! Define i,j ranges for southern LBC region
!---------------------------------------------------------------- 

        row_start_pt = 1 - halo_i
        row_end_pt   = row_length + halo_i

        first_row    = 1 - halo_j
        last_row     = rimwidth

!----------------------------------------------------------------
! Define length of rows in LBC data for exner and u points  
!----------------------------------------------------------------
        
        lbc_row_len   = ext_row_len
        lbc_row_len_u = ext_row_len
        
        DO j = first_row, last_row
           jj = j + halo_j-1
          DO i = row_start_pt, row_end_pt
            ii = i + halo_i-1

            IF (ns_weights(i,j) /= 0.0) THEN

!----------------------------------------------------------
! Calculate LBC address for given i,j location
!  for exner points, u points and v points
!----------------------------------------------------------

              lbc_address(i,j)   = lbc_start(PSouth)   + jj*lbc_row_len   + ii
              lbc_address_u(i,j) = lbc_start_u(PSouth) + jj*lbc_row_len_u + ii 
              lbc_address_v(i,j) = lbc_start_v(PSouth) + jj*lbc_row_len   + ii

!--------------------------------------------------------
!
! Before horiz. interpolation must calculate
! location of i,j-1 points in LBC data (LBC_address_jm1)
! (required for horiz. interpolation of v).
! Applying BCs as appropriate.
!
!--------------------------------------------------------

              IF (j > first_row) THEN   
                lbc_address_jm1(i,j) = lbc_start_v(PSouth) +              &
                                      (jj-1)*lbc_row_len + ii
              ELSE
                lbc_address_jm1(i,j) = lbc_start_v(PSouth) +              &
                                       jj*lbc_row_len + ii
              END IF

!--------------------------------------------------------
! Calculate LBC_ADDRESS for i-1,j points
! (required for horiz. interpolation of u)      
! Applying BCs as appropriate.
!--------------------------------------------------------

              IF (i > row_start_pt) THEN
                lbc_address_im1(i,j) = lbc_start_u(PSouth) +            &
                                       jj*lbc_row_len_u + ii-1
              END IF 
             
              IF (i == row_start_pt) THEN
                lbc_address_im1(i,j) = lbc_start_u(PSouth) +            &
                                       jj*lbc_row_len_u + ii
              END IF
             
              IF (i == row_end_pt) THEN
                lbc_address_u(i,j) = lbc_start_u(PSouth) +              &
                                     jj*lbc_row_len_u + ii-1
              END IF

            END IF  !ns_weights(i,j) /= 0.0
          END DO
        END DO

        DO j = first_row, last_row
          DO i = row_start_pt, row_end_pt
            IF (ns_weights(i,j) /= 0.0) THEN

!--------------------------------------------------------
! Calculate new balanced values of Exner 
!--------------------------------------------------------

              iloc = lbc_address(i,j)

              work1(i,j,1) = 0.5*(u_lbc(lbc_address_u(i,j),1)         &
                                + u_lbc(lbc_address_im1(i,j),1))      &
                                 *f2_comp_eh(i,j)                     &
                           - 0.5*(v_lbc(lbc_address_v(i,j),1)         &
                                + v_lbc(lbc_address_jm1(i,j),1))      &
                                 *f1_comp_eh(i,j) 
                                   
            END IF
          END DO
        END DO

! NB. Broken into two separate loops to make vectorisable

        DO k = 1, model_levels-1
          kp = k + 1
          DO j = first_row, last_row
            DO i = row_start_pt, row_end_pt
              IF (ns_weights(i,j) /= 0.0) THEN

!--------------------------------------------------------
! Calculate new balanced values of Exner 
!-------------------------------------------------------- 

                iloc = lbc_address(i,j)
                
                work1(i,j,kp) = 0.5*(u_lbc(lbc_address_u(i,j),kp)       &
                                  + u_lbc(lbc_address_im1(i,j),kp))     &
                                   *f2_comp_eh(i,j)                     &
                              - 0.5*(v_lbc(lbc_address_v(i,j),kp)       &
                                  + v_lbc(lbc_address_jm1(i,j),kp))     &
                                   *f1_comp_eh(i,j)
                               
                denom1 = cp * theta_wet_lbc(iloc,k) /                   &
                         (r_rho_levels(i,j,kp) - r_rho_levels(i,j,k))


                IF ( .NOT. l_cartesian ) THEN
                  u_av = intw_rho2w(k,1)*(                              &
                         0.5*(u_lbc(lbc_address_u(i,j),  kp)   +        &
                              u_lbc(lbc_address_im1(i,j),kp))) +        &
                         intw_rho2w(k,2)*(                              &
                         0.5*(u_lbc(lbc_address_u(i,j),  k)    +        &
                              u_lbc(lbc_address_im1(i,j),k)))

                  v_av = intw_rho2w(k,1)*(                              &
                         0.5*(v_lbc(lbc_address_v(i,j),  kp)   +        &
                              v_lbc(lbc_address_jm1(i,j),kp))) +        &
                         intw_rho2w(k,2)*(                              &
                         0.5*(v_lbc(lbc_address_v(i,j),  k)    +        &
                              v_lbc(lbc_address_jm1(i,j),k))) 

                  metric_term = (u_av**2 + v_av**2)/r_theta_levels(i,j,k)
                ELSE
                  metric_term = 0.0 
                END IF   

                dpdz = (metric_term - g_theta_eh(i,j,k)                 &
                     + 0.5*(work1(i,j,kp) + work1(i,j,k)))/denom1

                exner_lbc(iloc,kp) =  exner_lbc(iloc,k) + dpdz

!-------------------------------------------------------- 
! Calculate density to balance the new pressures
! using full non-linear equation of state
!--------------------------------------------------------

                temp = intw_w2rho(k,1)*theta_lbc(iloc,k)                 &
                     + intw_w2rho(k,2)*theta_lbc(iloc,k-1) 
                     
                rho_lbc(iloc,k) = p_zero/(R*temp)                        &
                                 *exner_lbc(iloc,k)**((1.0-kappa)/kappa)

!---------------------------------------------------------
! Set vertical velocities to zero
!---------------------------------------------------------

                W_LBC(iloc,K) = 0.0
           
              END IF
            END DO
          END DO
        END DO

        k = model_levels
        DO j = first_row, last_row
          DO i = row_start_pt, row_end_pt
            IF (ns_weights(i,j) /= 0.0) THEN

!--------------------------------------------------------
! Calculate density to balance the new pressures
! using full non-linear equation of state
!-------------------------------------------------------- 

              iloc = lbc_address(i,j)

              temp = intw_w2rho(k,1)*theta_lbc(iloc,k)                 &
                   + intw_w2rho(k,2)*theta_lbc(iloc,k-1) 

              rho_lbc(iloc,k) = p_zero/(R*temp)                        &
                                 *exner_lbc(iloc,k)**((1.0-kappa)/kappa)

! This shouldnt be needed so just copy the value from below
              exner_lbc(iloc,k+1) = exner_lbc(iloc,k)

            END IF
          END DO
        END DO

      END IF ! IF (at_extremity(PSouth))

!------------------------------------------------------------------------
! 4.2 Northern LBC region
!------------------------------------------------------------------------

      IF (at_extremity(PNorth)) THEN

        row_start_pt = 1 - halo_i
        row_end_pt   = row_length + halo_i
        first_row    = rows - rimwidth + 1
        last_row     = rows + halo_j

        lbc_row_len   = ext_row_len
        lbc_row_len_u = ext_row_len

        IF (at_extremity(PSouth)) THEN
          first_row_of_lbc = rimwidth + 1
        ELSE
          first_row_of_lbc = 1 - halo_j
        END IF

        DO j = first_row, last_row
          weights_j = rows + 1 - j
          jj        = j - (rows - rimwidth) - 1
          DO i = row_start_pt, row_end_pt
            ii =  i + halo_i - 1
            IF (ns_weights(I,weights_J) /= 0.0) THEN

!----------------------------------------------------------
! Calculate LBC address for given i,j location
!  for exner points, u points and v points
!----------------------------------------------------------

              lbc_address(i,j)   = lbc_start(PNorth) + jj*lbc_row_len + ii
              lbc_address_v(i,j) = lbc_start_v(PNorth) + (jj+1)*lbc_row_len + ii
                         
! NB. J+1 above arises because first row of v data in LBCs
!    is one row lower than the first row of exner/u data

              lbc_address_u(i,j) = lbc_start_u(PNorth) + jj*lbc_row_len_u + ii
        
!--------------------------------------------------------
!
! Before horiz. interpolation must calculate
! location of i,j-1 points in LBC data (LBC_address_jm1)
! (required for horiz. interpolation of v)
! Applying BCs as appropriate.
!
!--------------------------------------------------------
!--------------------------------------------------------
! Calculate LBC_ADDRESS for i-1,j points
! (required for horiz. interpolation of u)
! Applying BCs as appropriate.
!--------------------------------------------------------

              IF (i > row_start_pt) THEN
                lbc_address_im1(i,j) = lbc_start_u(PNorth) +                  &
                                       jj*lbc_row_len_u + ii-1
              END IF

              IF (i == row_start_pt) THEN
                lbc_address_im1(i,j) = lbc_start_u(Pnorth) +                  &
                                       jj*lbc_row_len_u + ii
              END IF

              IF (i == row_end_pt) THEN
                lbc_address_u(i,j) = lbc_start_u(PNorth) +                    &
                                     jj*lbc_row_len_u + ii-1
              END IF

!--------------------------------------------------------
! Calculate LBC_address for i,j-1 for 
! (required for interpolation of v)
! Applying BCs as appropriate
!--------------------------------------------------------

              lbc_address_jm1(i,j) = lbc_start_v(PNorth) +                    &
                                     jj*lbc_row_len + ii    

! NB. J not J-1 since first v point in North LBC
!    is one point before first exner point.  

              IF (j == last_row)THEN
                lbc_address_v(i,j) = lbc_start_v(PNorth) +                    &
                                     jj*lbc_row_len + ii
              END IF
            END IF
          END DO
        END DO

! NB. Broken into two separate loops to make vectorisable

        DO j = first_row, last_row
           weights_j = rows + 1 - j 
           DO i = row_start_pt, row_end_pt
             IF (ns_weights(I,weights_J) /= 0.0) THEN

               iloc = LBC_address(i,j)

               work1(i,j,1) = 0.5*(u_lbc(lbc_address_u(i,j),1)        &
                                 + u_lbc(lbc_address_im1(i,j),1))     &
                                  *f2_comp_eh(i,j)                    &
                            - 0.5*(v_lbc(lbc_address_v(i,j),1)        &
                                 + v_lbc(lbc_address_jm1(i,j),1))     &
                                  *f1_comp_eh(i,j)  

             END IF
           END DO
        END DO

        DO k = 1, model_levels-1
          kp = k + 1
          DO J = first_row, last_row
            weights_j = rows + 1 - j
            DO I = row_start_pt, row_end_pt
              IF (ns_weights(I,weights_j) /= 0.0) THEN

!--------------------------------------------------------
! Calculate new balanced values of Exner 
!--------------------------------------------------------

                iloc = LBC_address(i,j)
                
                work1(i,j,kp) = 0.5*(u_lbc(lbc_address_u(i,j),kp)       &
                                  + u_lbc(lbc_address_im1(i,j),kp))     &
                                   *f2_comp_eh(i,j)                     &
                             - 0.5*(v_lbc(lbc_address_v(i,j),kp)        &
                                  + v_lbc(lbc_address_jm1(i,j),kp))     &
                                   *f1_comp_eh(i,j) 

                denom1 = cp * theta_wet_lbc(iloc,k) /                    &
                         (r_rho_levels(i,j,kp) - r_rho_levels(i,j,k))

                IF ( .NOT. l_cartesian ) THEN
                  u_av = intw_rho2w(k,1)*(                              &
                         0.5*(u_lbc(lbc_address_u(i,j),  kp)   +        &
                              u_lbc(lbc_address_im1(i,j),kp))) +        &
                         intw_rho2w(k,2)*(                              &
                         0.5*(u_lbc(lbc_address_u(i,j),  k)    +        &
                              u_lbc(lbc_address_im1(i,j),k)))

                  v_av = intw_rho2w(k,1)*(                              &
                         0.5*(v_lbc(lbc_address_v(i,j),  kp)   +        &
                              v_lbc(lbc_address_jm1(i,j),kp))) +        &
                         intw_rho2w(k,2)*(                              &
                         0.5*(v_lbc(lbc_address_v(i,j),  k)    +        &
                              v_lbc(lbc_address_jm1(i,j),k))) 

                  metric_term = (u_av**2 + v_av**2)/r_theta_levels(i,j,k)
                ELSE
                  metric_term = 0.0 
                END IF 

                dpdz = (metric_term - g_theta_eh(i,j,k)                 &
                     + 0.5*(work1(i,j,kp) + work1(i,j,k)))/denom1

                exner_lbc(iloc,kp) = exner_lbc(iloc,k) + dpdz

!--------------------------------------------------------
! Calculate density to balance the new pressures
! using full non-linear equation of state
!--------------------------------------------------------
                temp = intw_w2rho(k,1)*theta_lbc(iloc,k)                 &
                     + intw_w2rho(k,2)*theta_lbc(iloc,k-1)  

                rho_lbc(iloc,k) = p_zero/(R*temp)                        &
                                 *exner_lbc(iloc,k)**((1.0-kappa)/kappa)
!--------------------------------------------------------
! Set vertical velocities to zero
!--------------------------------------------------------

                W_LBC(iloc,K) = 0.0

              END IF   ! ns_weights(I,weights_j) /= 0.0
            END DO
          END DO
        END DO

        k = model_levels
        DO J = first_row, last_row
          weights_J = ROWS + 1 - J
          DO I = row_start_pt, row_end_pt
            IF (ns_weights(I,weights_J) /= 0.0) THEN

!--------------------------------------------------------
! Calculate density to balance the new pressures
! using full non-linear equation of state
!--------------------------------------------------------

              iloc = LBC_address(i,j) 

              temp = intw_w2rho(k,1)*theta_lbc(iloc,k)                 &
                   + intw_w2rho(k,2)*theta_lbc(iloc,k-1)

              rho_lbc(iloc,k) = p_zero/(R*temp)                        &
                                 *exner_lbc(iloc,k)**((1.0-kappa)/kappa)

! This shouldnt be needed so just copy the value from below
              exner_lbc(iloc,k+1) = exner_lbc(iloc,k) 

            END IF
          END DO
        END DO
      END IF !  IF (AT_EXTREMITY(PNorth))

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! 4.3 Western LBC region
!------------------------------------------------------------------------

      IF (at_extremity(PWest)) THEN

        row_start_pt = 1-halo_i
        row_end_pt   = n_rims_to_do

        IF (at_extremity(PSouth)) THEN
          first_row        = rimwidth + 1
          first_row_of_lbc = rimwidth + 1
        ELSE
          first_row_of_lbc = 1 - halo_j
          first_row        = 1 - halo_j
        END IF

        IF (at_extremity(PNorth)) THEN
          last_row = rows - rimwidth
        ELSE
          last_row = rows + halo_j
        END IF

        lbc_row_len = halo_i + rimwidth

! NB. Loop order k,i,j to improve speed
   
        DO j = first_row, last_row
          gj = datastart(2) + j - 1
          jj = j - first_row_of_lbc
          DO i = row_start_pt, row_end_pt 
            ii = i + halo_i - 1
            IF (ew_weights(i,j) /= 0.0) THEN

!----------------------------------------------------------
! Calculate LBC address for given i,j location
!  for exner points, u points and v points
!----------------------------------------------------------

              lbc_address(i,j)   = lbc_start(PWest)   + jj*lbc_row_len + ii
              lbc_address_u(i,j) = lbc_start_u(PWest) + jj*lbc_row_len + ii
              lbc_address_v(i,j) = lbc_start_v(PWest) + jj*lbc_row_len + ii

!--------------------------------------------------------
!
! Before horiz. interpolation must calculate
! location of i,j-1 points in LBC data
! (LBC_address_jm1(i,j))
! (required for horiz. interpolation of v)
!
!--------------------------------------------------------
!--------------------------------------------------------
! At bottom of West LBC region, 
!    v point south of exner point is in South LBC region.
! At top of West LBC region,
!    v point north of exner point is in North LBC region.
! Since each region requires a different calculation of
! LBC_ADDRESS(i,j) must first determine which LBC region 
! the (i,j) and (i,j-1) points are in
!--------------------------------------------------------
        
              IF (gj == rimwidth + 1) THEN
                lbcarea = 3   ! South LBC
              ELSE IF( gj == global_rows-rimwidth) THEN
                lbcArea = 1   ! North LBC
              ELSE
                lbcarea = 4   ! West LBC
              END IF

!--------------------------------------------------------
! Calculate LBC_address_jm1(i,j) for  
!  appropriate LBC region.
! Applying BCs as appropriate
!--------------------------------------------------------
                
              IF( lbcarea == 4)THEN  ! West LBC
                IF (j /= 1-halo_j) THEN
                  lbc_address_jm1(i,j) = lbc_start_v(PWest) +                 &
                                         (jj-1)*lbc_row_len + ii
                ELSE
                  lbc_address_jm1(i,j) = lbc_start_v(PWest) +                 &
                                         jj*lbc_row_len + ii
                END IF
              END IF

              IF (lbcarea == 3) THEN  ! South LBC 
                lbc_address_jm1(i,j) = lbc_start_v(PSouth) +                  &
                                      (j-1 + halo_j-1)*ext_row_len + ii
              END IF

              IF (lbcarea == 1) THEN  ! North LBC 
                lbc_address_v(i,j) = lbc_start_v(PNorth) +                     &
                                    (j+1 - (rows - rimwidth)-1)*ext_row_len + ii

                lbc_address_jm1(i,j) = lbc_start_v(PWest) +                    &
                                      (jj-1)*lbc_row_len + ii
              END IF
 
!--------------------------------------------------------
! Calculate LBC_ADDRESS for i-1,j points
! (required for horiz. interpolation of u)      
!--------------------------------------------------------

              IF (i > row_start_pt) THEN
                lbc_address_im1(i,j) = lbc_start_u(PWest) +             &
                                       jj* lbc_row_len +  ii-1
              ELSE
                lbc_address_im1(i,j) = lbc_start_u(PWest) +             &
                                       jj*lbc_row_len + ii
              END IF
            END IF
          END DO
        END DO

! NB. Broken into two separate loops to make vectorisable

        DO i = row_start_pt, row_end_pt
          DO j = first_row, last_row
            IF (ew_weights(i,J) /= 0.0) THEN
!--------------------------------------------------------
! Calculate new balanced values of Exner 
!-------------------------------------------------------- 
              iloc = LBC_address(i,j) 
                 
              work1(i,j,1) = 0.5*(u_lbc(lbc_address_u(i,j),1)        &
                                   + u_lbc(lbc_address_im1(i,j),1))     &
                                   *f2_comp_eh(i,j)                     &
                              - 0.5*(v_lbc(lbc_address_v(i,j),1)        &
                                   + v_lbc(lbc_address_jm1(i,j),1))     &
                                   *f1_comp_eh(i,j) 

            END IF
          END DO
        END DO

        DO K = 1, model_levels-1
          kp = k + 1
          DO I = row_start_pt, row_end_pt
            DO J = first_row, last_row
              IF (ew_weights(I,J) /= 0.0) THEN
                 
!--------------------------------------------------------
! Calculate new balanced values of Exner 
!--------------------------------------------------------
                iloc = LBC_address(i,j)
                
                work1(i,j,kp) = 0.5*(u_lbc(lbc_address_u(i,j),kp)       &
                                  + u_lbc(lbc_address_im1(i,j),kp))     &
                                   *f2_comp_eh(i,j)                     &
                             - 0.5*(v_lbc(lbc_address_v(i,j),kp)        &
                                  + v_lbc(lbc_address_jm1(i,j),kp))     &
                                   *f1_comp_eh(i,j)            
                               
                denom1 = cp * theta_wet_lbc(iloc,k) /                    &
                         (r_rho_levels(i,j,kp) - r_rho_levels(i,j,k))
                    
                IF ( .NOT. l_cartesian ) THEN
                  u_av = intw_rho2w(k,1)*(                              &
                         0.5*(u_lbc(lbc_address_u(i,j),  kp)   +        &
                              u_lbc(lbc_address_im1(i,j),kp))) +        &
                         intw_rho2w(k,2)*(                              &
                         0.5*(u_lbc(lbc_address_u(i,j),  k)    +        &
                              u_lbc(lbc_address_im1(i,j),k)))           
                  
                  v_av = intw_rho2w(k,1)*(                              &
                         0.5*(v_lbc(lbc_address_v(i,j),  kp)   +        &
                              v_lbc(lbc_address_jm1(i,j),kp))) +        &
                         intw_rho2w(k,2)*(                              &
                         0.5*(v_lbc(lbc_address_v(i,j),  k)    +        &
                              v_lbc(lbc_address_jm1(i,j),k))) 
                              
                  metric_term = (u_av**2 + v_av**2)/r_theta_levels(i,j,k)
                ELSE
                  metric_term = 0.0 
                END IF 

                dpdz = (metric_term - g_theta_eh(i,j,k)                 &
                     + 0.5*(work1(i,j,kp) + work1(i,j,k)))/denom1

                exner_lbc(iloc,kp) = exner_lbc(iloc,k) + dpdz
              
!--------------------------------------------------------
! Calculate density to balance the new pressures
! using full non-linear equation of state
!--------------------------------------------------------
                temp = intw_w2rho(k,1)*theta_lbc(iloc,k)                 & 
                     + intw_w2rho(k,2)*theta_lbc(iloc,k-1)
                
                rho_lbc(iloc,k) = p_zero/(R*temp)                        &
                                 *exner_lbc(iloc,k)**((1.0-kappa)/kappa)

!--------------------------------------------------------
! Set vertical velocities to zero
!--------------------------------------------------------

                w_lbc(iloc,k) = 0.0

              END IF  !ew_weights(I,J) /= 0.0
        
            END DO
          END DO
        END DO

        k = model_levels
        DO j = first_row, last_row
          DO I = row_start_pt, row_end_pt
            IF (ew_weights(I,J) /= 0.0) THEN
!--------------------------------------------------------
! Calculate density to balance the new pressures
! using full non-linear equation of state
!--------------------------------------------------------
              iloc = LBC_address(i,j)

              temp = intw_w2rho(k,1)*theta_lbc(iloc,k)                 & 
                   + intw_w2rho(k,2)*theta_lbc(iloc,k-1)

              rho_lbc(iloc,k) = p_zero/(R*temp)                        &
                                 *exner_lbc(iloc,k)**((1.0-kappa)/kappa) 

! This shouldnt be needed so just copy the value from below

              exner_lbc(iloc,k+1) = exner_lbc(iloc,k)
            END IF
          END DO
        END DO

      END IF ! IF (at_extremity(PWest))

!------------------------------------------------------------------------
! 4.3 Eastern LBC region
!------------------------------------------------------------------------

      IF (at_extremity(PEast)) THEN

        row_start_pt = row_length - n_rims_to_do + 1
        row_end_pt   = row_length + halo_i

        IF (at_extremity(PSouth)) THEN
          first_row = rimwidth + 1
          first_row_of_lbc = rimwidth + 1
        ELSE
          first_row_of_lbc = 1 - halo_j
          first_row = 1 - halo_j
        END IF

        IF (at_extremity(PNorth)) THEN
          last_row = rows - rimwidth
        ELSE
          last_row = rows + halo_j
        END IF

        lbc_row_len = halo_i + rimwidth
        first_pt_of_lbc = row_length - rimwidth + 1

! NB. Loop order i,j to improve speed
    
        DO j = first_row, last_row
          gj = datastart(2) + j - 1
          jj = j - first_row_of_lbc
          DO i = row_start_pt, row_end_pt
            weights_i = first_pt_of_lbc + rimwidth - i
            ii        = i - first_pt_of_lbc
            IF (ew_weights(weights_I,J) .NE. 0.0) THEN

!----------------------------------------------------------
! Calculate LBC address for given i,j location
!  for exner points, u points and v points
!----------------------------------------------------------

              lbc_address(i,j)   = lbc_start(PEast)   + jj*lbc_row_len + ii
              lbc_address_v(i,j) = lbc_start_v(PEast) + jj*lbc_row_len + ii
              lbc_address_u(i,j) = lbc_start_u(PEast) + jj*lbc_row_len + ii + 1

! NB. I+1 arises because first u point in East LBC
!    is one point before first exner point.  

!--------------------------------------------------------
!
! Before horiz. interpolation must calculate
! location of i,j-1 points in LBC data (LBC_address_jm1)
! (required for horiz. interpolation of v)
!
!--------------------------------------------------------
!--------------------------------------------------------
! At bottom of East LBC region, 
!    v point south of exner point is in South LBC region.
! At top of East LBC region,
!    v point north of exner point is in North LBC region.
! Since each region requires a different calculation of
! LBC_ADDRESS(i,j) must first determine which LBC region 
! the (i,j) and (i,j-1) points are in
!--------------------------------------------------------

              IF (i == row_end_pt) THEN
                lbc_address_u(i,j) = lbc_start_u(PEast) + jj*lbc_row_len + ii
              END IF

              IF (gj == rimwidth+1) THEN
                lbcArea = 3 ! South
              ELSE IF (gj == global_rows - rimwidth) THEN
                lbcarea = 1 ! North
              ELSE 
                lbcarea = 2 ! East
              END IF

!--------------------------------------------------------
! Calculate LBC_address_jm1 appropriate LBC region
! Applying BCs as appropriate
!--------------------------------------------------------

              IF (lbcarea == 2) THEN
                IF (j /= 1-halo_j) THEN

                  lbc_address_jm1(i,j) = lbc_start_v(PEast)+               &
                                        (jj-1)*(halo_i +rimwidth) + ii
                ELSE
                  lbc_address_jm1(i,j) = lbc_start_v(PEast)+               &
                                         jj*(halo_i + rimwidth) + ii
                END IF
              END IF
        
              IF (lbcarea == 3) THEN      
                lbc_address_jm1(i,j) = lbc_start_v(PSouth) +               &
                                      (j-2 + halo_j)*ext_row_len           &
                                      + i + halo_i-1
              END IF
           
              IF (lbcarea == 1) THEN
                lbc_address_v(i,j) = lbc_start_v(PNorth) +                 &
                                     (j+1 - (rows - rimwidth)-1)           &
                                    *ext_row_len + i + halo_i - 1

                lbc_address_jm1(i,j) = lbc_start_v(PEast)+                 &
                                      (jj-1)*lbc_row_len + ii
              END IF

!--------------------------------------------------------
! Calculate LBC_ADDRESS for i-1,j points
! (required for horiz. interpolation of u)
!--------------------------------------------------------

              lbc_address_im1(i,j) = lbc_start_u(PEast)+                   &
                                     jj*lbc_row_len + ii

!NB. I not I-1 since first u point in East LBC
!    is one point before first exner point. 

            END IF
          END DO
        END DO

! NB. Broken into two separate loops to make vectorisable
        DO i = row_start_pt, row_end_pt
          weights_I = first_pt_of_lbc + rimwidth - i
          DO j = first_row, last_row   
            IF (ew_weights(weights_i,j) /= 0.0) THEN

!--------------------------------------------------------
! Calculate new balanced values of Exner 
!--------------------------------------------------------
              iloc = lbc_address(i,j)

              work1(i,j,1) = 0.5*(u_lbc(lbc_address_u(i,j),1)        &
                                  + u_lbc(lbc_address_im1(i,j),1))     &
                                   *f2_comp_eh(i,j)                    &
                             - 0.5*(v_lbc(lbc_address_v(i,j),1)        &
                                  + v_lbc(lbc_address_jm1(i,j),1))     &
                                   *f1_comp_eh(i,j) 

            END IF
          END DO
        END DO

        DO k = 1, model_levels - 1
          kp = k + 1
          DO i = row_start_pt, row_end_pt
            weights_i = first_pt_of_lbc + rimwidth - i
            DO j = first_row, last_row
              IF (ew_weights(weights_i,j) /= 0.0) THEN

!--------------------------------------------------------
! Calculate new balanced values of Exner 
!--------------------------------------------------------

                iloc = lbc_address(i,j)

                work1(i,j,kp) = 0.5*(u_lbc(lbc_address_u(i,j),kp)       &
                                  + u_lbc(lbc_address_im1(i,j),kp))     &
                                   *f2_comp_eh(i,j)                     &
                             - 0.5*(v_lbc(lbc_address_v(i,j),kp)        &
                                  + v_lbc(lbc_address_jm1(i,j),kp))     &
                                   *f1_comp_eh(i,j) 


                denom1 = cp * theta_wet_lbc(iloc,k) /                    &
                        (r_rho_levels(i,j,kp) - r_rho_levels(i,j,k))

                IF ( .NOT. l_cartesian ) THEN
                  u_av = intw_rho2w(k,1)*(                              &
                         0.5*(u_lbc(lbc_address_u(i,j),  kp)   +        &
                              u_lbc(lbc_address_im1(i,j),kp))) +        &
                         intw_rho2w(k,2)*(                              &
                         0.5*(u_lbc(lbc_address_u(i,j),  k)    +        &
                              u_lbc(lbc_address_im1(i,j),k)))           

                  v_av = intw_rho2w(k,1)*(                              &
                         0.5*(v_lbc(lbc_address_v(i,j),  kp)   +        &
                              v_lbc(lbc_address_jm1(i,j),kp))) +        &
                         intw_rho2w(k,2)*(                              &
                         0.5*(v_lbc(lbc_address_v(i,j),  k)    +        &
                              v_lbc(lbc_address_jm1(i,j),k))) 

                  metric_term = (u_av**2 + v_av**2)/r_theta_levels(i,j,k)
                ELSE
                  metric_term = 0.0 
                END IF 

                dpdz = (metric_term - g_theta_eh(i,j,k)                 &
                     + 0.5*(work1(i,j,kp) + work1(i,j,k)))/denom1
                               
                exner_lbc(iloc,kp) = exner_lbc(iloc,k) + dpdz
              

!--------------------------------------------------------
! Calculate density to balance the new pressures
! this is rho_wet * r * r
!--------------------------------------------------------
                temp = intw_w2rho(k,1)*theta_lbc(iloc,k)                 &
                     + intw_w2rho(k,2)*theta_lbc(iloc,k-1)
                     
                rho_lbc(iloc,k) = p_zero/(R*temp)                        &
                                 *exner_lbc(iloc,k)**((1.0-kappa)/kappa)
                   
!--------------------------------------------------------
! Set vertical velocities to zero
!--------------------------------------------------------

                w_lbc(iloc,K) = 0.0

              END IF

            END DO
          END DO
        END DO

        k = model_levels
        DO j = first_row, last_row
          DO i = row_start_pt, row_end_pt
            weights_i = first_pt_of_lbc + rimwidth - i
            IF (ew_weights(weights_i,j) /= 0.0) THEN

!--------------------------------------------------------
! Calculate density to balance the new pressures
! this is rho_wet * r * r
!--------------------------------------------------------

              iloc = LBC_address(i,j)  

              temp = intw_w2rho(k,1)*theta_lbc(iloc,k)                 &
                     + intw_w2rho(k,2)*theta_lbc(iloc,k-1)
                     
              rho_lbc(iloc,k) = p_zero/(R*temp)                        &
                                 *exner_lbc(iloc,k)**((1.0-kappa)/kappa)  

! This shouldnt be needed so just copy the value from below

              exner_lbc(iloc,k+1) = exner_lbc(iloc,k)                
            END IF
          END DO
        END DO
      END IF ! IF (AT_EXTREMITY(PEast))

      IF (lhook) CALL dr_hook('EG_BALANCE_LBC_VALUES',zhook_out,zhook_handle)

      END SUBROUTINE eg_balance_lbc_values
      END MODULE eg_balance_lbc_values_mod
