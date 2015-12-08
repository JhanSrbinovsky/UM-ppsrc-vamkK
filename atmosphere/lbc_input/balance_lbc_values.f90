! *****************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*********************************

!+ Balance Exner in LBCs and set rho using equation of state.
!  Reset w=0

      SUBROUTINE BALANCE_LBC_VALUES(                                    &
        EXNER_LBC, RHO_LBC, THETA_LBC, Q_LBC, W_LBC, W_ADV_LBC,         &
        U_LBC, V_LBC,                                                   &
        R_RHO_LEVELS, R_THETA_LEVELS,                                   &
        ROW_LENGTH, ROWS, model_levels, wet_levels,                     &
        HALO_I, HALO_J, LENRIM,                                         &
        LENRIMU, LENRIMV,                                               &
        LBC_START, LBC_START_U, LBC_START_V,                            &
        RIMWIDTH, N_RIMS_TO_DO, RIMWEIGHTS, AT_EXTREMITY,               &
        DELTA_PHI, DELTA_LAMBDA,                                        &
        BASE_PHI, BASE_LAMBDA,                                          &
        DATASTART, LAT_ROT_NP,                                          &  
        GLOBAL_ROW_LENGTH, GLOBAL_ROWS, L_int_uvw_lbc                   &
        )

USE atmos_constants_mod, ONLY: cp, p_zero, r, c_virtual,                &
                               recip_kappa

USE VECTLIB_MOD, ONLY :                                                 &
      SIN_V,  COS_V

USE earth_constants_mod, ONLY : g, two_omega

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


      USE UM_ParParams
      IMPLICIT NONE

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

      LOGICAL, INTENT(IN) :: L_int_uvw_lbc
                             ! true if interpolated winds are used
                             ! for advection in lateral boundaries

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
       , DATASTART(2)                                                   &  
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
       , THETA_LBC(LENRIM, MODEL_LEVELS)                                & 
       , Q_LBC(LENRIM, WET_LEVELS)                                      & 
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
      INTEGER :: gj         ! j position in global field
      INTEGER :: gi         ! i position in global field
      INTEGER :: weights_I  ! modified I to point to ew_weights array
      INTEGER :: weights_J  ! modified J to point to ns_weights array
      INTEGER :: rim        ! rim number (at corners)
      INTEGER :: rim_I      ! modified I to point to RIMWEIGHTS array
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

      REAL :: ns_weights(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:RIMWIDTH)
      REAL :: ew_weights(1-HALO_I:RIMWIDTH,1-HALO_J:ROWS+HALO_J)
                            ! Arrays which contain the weights
                            ! to apply to the LBCs at each point
      REAL :: U_LBC_X1 !u wind interpolated onto rho points on level k       
      REAL :: U_LBC_X2 !u wind interpolated onto rho points on level k+1     
      REAL :: U_LBC_XZ(1-halo_i:row_length+halo_i,                      &
                     1-halo_j:rows+halo_j, model_levels)               
                       !u wind interpolated onto theta points on level k    
      REAL :: V_LBC_Y1 !v wind interpolated onto rho points on level k       
      REAL :: V_LBC_Y2 !v wind interpolated onto rho points on level k+1     
      REAL :: V_LBC_YZ(1-halo_i:row_length+halo_i,                      &
                   1-halo_j:rows+halo_j, model_levels)
                     !v wind interpolated onto theta points on level k+1 
      REAL :: goverCp                                                     
                                                 
      REAL :: temp          ! virtual potential temperature
      REAL :: weight        ! vertical interpolation weights
      REAL :: latitude_rot  ! latitude on rotated grid
      REAL :: longitude_rot ! longitude on rotated grid                                       
      REAL :: f1_at_theta(1-HALO_I:global_row_length+HALO_I,            &
                        1-HALO_J:global_rows+HALO_J)
      REAL :: f2_at_theta(1-HALO_I:global_row_length+HALO_I,            &
                        1-HALO_J:global_rows+HALO_J)
                        
      REAL :: sin_lat_rot_NP      ! = SIN(lat_rot_NP)
      REAL :: cos_lat_rot_NP      ! = COS(lat_rot_NP)
      REAL :: tmp1(global_row_length+2*(HALO_I)+global_rows+2*(HALO_J))
      REAL :: tmp2(global_row_length+2*(HALO_I)+global_rows+2*(HALO_J))
      REAL :: tmp3(global_row_length+2*(HALO_I)+global_rows+2*(HALO_J))
      INTEGER :: leni
      INTEGER :: lenj
      INTEGER :: total_length

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
!
! Code
!========================================================================
! 1. Set up local constant for later
!========================================================================

      IF (lhook) CALL dr_hook('BALANCE_LBC_VALUES',zhook_in,zhook_handle)

      goverCp = g / Cp

!========================================================================
! 2. Set up the "*_weights" arrays
!       -only want to impose BCs where blending weight is non-zero
!========================================================================
! 2.1 Do the North/South regions

      IF (AT_EXTREMITY(PNorth) .OR. AT_EXTREMITY(PSouth)) THEN
! We will need the ns_weights array

! First we do the halo area (South : qlp, North : abf)

! The NS edge (South: q, North : a)
        DO J=1-HALO_J,0
          DO I=1-HALO_I,ROW_LENGTH+HALO_I
            ns_weights(I,J)=1.0
          ENDDO ! I
        ENDDO ! J

! The Western edge ( South: l, North : b)
        IF (AT_EXTREMITY(PWest)) THEN
        DO J=1,RIMWIDTH
          DO I=1-HALO_I,0
            ns_weights(I,J)=1.0
          ENDDO ! I
        ENDDO ! J
        ENDIF ! IF (AT_EXTREMITY(PWest))

! The Eastern edge (South : p, North : f)
        IF (AT_EXTREMITY(PEast)) THEN
        DO J=1,RIMWIDTH
          DO I=ROW_LENGTH+1,ROW_LENGTH+HALO_I
            ns_weights(I,J)=1.0
          ENDDO ! I
        ENDDO ! J
        ENDIF ! IF (AT_EXTREMITY(PEast))

! And now the boundary areas ( South : mno, North : cde )
! The Western corner (South : m, North : c)
        IF (AT_EXTREMITY(PWest)) THEN
        DO J=1,RIMWIDTH
          DO I=1,RIMWIDTH
            rim=MIN(I,J)
            IF (rim  <=  N_RIMS_TO_DO) THEN
              ns_weights(I,J)=RIMWEIGHTS(rim)
            ELSE ! This rim not required
              ns_weights(I,J)=0.0
            ENDIF !  IF (rim  <=  N_RIMS_TO_DO)
          ENDDO ! I
        ENDDO ! J
        ENDIF ! IF (AT_EXTREMITY(PWest))

! The Eastern corner (South : o, North : e)
        IF (AT_EXTREMITY(PEast)) THEN
        DO J=1,RIMWIDTH
          DO I=ROW_LENGTH-RIMWIDTH+1,ROW_LENGTH
            rim_I=ROW_LENGTH+1-I
            rim=MIN(rim_I,J)
            IF (rim  <=  N_RIMS_TO_DO) THEN
              ns_weights(I,J)=RIMWEIGHTS(rim)
            ELSE ! This rim not required
              ns_weights(I,J)=0.0
            ENDIF !  IF (rim  <=  N_RIMS_TO_DO)
          ENDDO ! I
        ENDDO ! J
        ENDIF ! IF (AT_EXTREMITY(PEast))

! The bit between the corners (South : n, North : d)
        IF (AT_EXTREMITY(PEast)) THEN
          row_end_pt=ROW_LENGTH-RIMWIDTH
        ELSE
          row_end_pt=row_length +halo_i
        ENDIF
        IF (AT_EXTREMITY(Pwest)) THEN
          row_start_pt=RIMWIDTH+1
        ELSE
          row_start_pt=1-halo_i
        ENDIF
        DO J=1,RIMWIDTH
          DO I=row_start_pt,row_end_pt
            IF (J  <=  N_RIMS_TO_DO) THEN
              ns_weights(I,J)=RIMWEIGHTS(J)
            ELSE ! This rim is not required
              ns_weights(I,J)=0.0
            ENDIF ! IF (J  <=  N_RIMS_TO_DO)
          ENDDO ! I
        ENDDO ! J

      ENDIF ! If we're at the North or South edge of the LAM

!-----------------------------------------------------------------
! 2.2 Do the East/West regions

      IF (AT_EXTREMITY(PWest) .OR. AT_EXTREMITY(PEast)) THEN
! We will need the ew_weights array

        IF (AT_EXTREMITY(PSouth)) THEN
          first_row=RIMWIDTH+1
        ELSE ! not at the South
          first_row=1-HALO_J
        ENDIF

        IF (AT_EXTREMITY(PNorth)) THEN
          last_row=ROWS-RIMWIDTH
        ELSE ! not at the North
          last_row=ROWS+HALO_J
        ENDIF

! First the halo area (West : g , East k)

        DO J=first_row,last_row
          DO I=1-HALO_I,0
            ew_weights(I,J)=1.0
          ENDDO ! I
        ENDDO ! J

! And now the boundary area (West : h, East : j)

        DO J=first_row,last_row
          DO I=1,RIMWIDTH
            IF (I  <=  N_RIMS_TO_DO) THEN
              ew_weights(I,J)=RIMWEIGHTS(I)
            ELSE ! this rim not required
              ew_weights(I,J)=0.0
            ENDIF ! IF (I  <=  N_RIMS_TO_DO)
          ENDDO ! I
        ENDDO ! J

      ENDIF ! IF at Western or Eastern edge

!========================================================================
! 3. Calculate coriolis terms on theta points
!========================================================================

      sin_lat_rot_NP = SIN(lat_rot_NP)
      cos_lat_rot_NP = COS(lat_rot_NP)
      leni=global_row_length+2*(HALO_I)
      lenj=global_rows+2*(HALO_J)
      total_length=leni+lenj

      DO GI=1-HALO_I,global_row_length+HALO_I
        tmp1(gi+HALO_I) = (base_lambda + (gi-1) * delta_lambda)
      END DO
      DO GJ=1-HALO_J,global_rows+HALO_J
        tmp1(leni+gj+HALO_J) = (base_phi + (gj-1) * delta_phi)
      END DO
      CALL SIN_V(total_length,tmp1,tmp2)
      CALL COS_V(total_length,tmp1,tmp3)
      DO GJ=1-HALO_J,global_rows+HALO_J
        DO GI=1-HALO_I,global_row_length+HALO_I
          f1_at_theta(gi,gj) = - two_omega *                            &
                           tmp2(gi+HALO_I) * cos_lat_rot_NP
          f2_at_theta(gi,gj) = two_omega *                              &
                         ( tmp3(leni+gj+HALO_J) * sin_lat_rot_NP -      &
                           tmp2(leni+gj+HALO_J) *                       &
                           tmp3(gi+HALO_I) * cos_lat_rot_NP )
        END DO !GI
      END DO !GJ

!========================================================================
! 4.  Now apply the boundary conditions to exner, rho and w
!========================================================================
!------------------------------------------------------------------------
! 4.1 Southern LBC region
!------------------------------------------------------------------------

      IF (AT_EXTREMITY(PSouth)) THEN

        !----------------------------------------------------------------
        ! Define i,j ranges for southern LBC region
        !---------------------------------------------------------------- 

        row_start_pt = 1 - HALO_I
        row_end_pt   = ROW_LENGTH + HALO_I
        
        first_row    = 1 - HALO_J
        last_row     = RIMWIDTH

        !----------------------------------------------------------------
        ! Define length of rows in LBC data for exner and u points  
        !----------------------------------------------------------------
        
        lbc_row_len = ROW_LENGTH + 2 * HALO_I

        IF(AT_EXTREMITY(PEast))THEN
          lbc_row_len_u = ROW_LENGTH + 2 * HALO_I -1
        ELSE
          lbc_row_len_u = ROW_LENGTH + 2 * HALO_I
        ENDIF
        
        DO K = 1, model_levels

          DO J = first_row, last_row
            gj = datastart(2) + j - 1

            DO I = row_start_pt, row_end_pt
              gi = datastart(1) + i - 1

              IF (ns_weights(I,J) .NE. 0.0) THEN

              !----------------------------------------------------------
              ! u and v are required on theta points.
              !
              ! Step A: Horizontal interpolation onto exner points
              ! Step B: Vertical interpolation onto theta points
              ! 
              ! Horizontal grid staggering:
              !
              !       |
              !     v(0,1)                        v(i,j)
              !       |
              !       |
              !   exner(0,1)       u(i-1,j)     exner(i,j)     u(i,j)
              !       |
              !       |
              !     v(0,0)                       v(i,j-1)
              !       |
              !         |
              !   exner(0,0)-------u(0,0)-------exner(1,0)-----u(1,0)
              !
              !                         
              !----------------------------------------------------------
              ! Vertical staggering :
              !
              !  u(i-1,j,k+1)  U_LBC_X2    u(i,j,k+1)        rho levels
              !
              !                U_LBC_XZ                      theta levels
              !
              !  u(i-1,j,k)    U_LBC_X1    u(i,j,k)          rho levels
              !----------------------------------------------------------
              !----------------------------------------------------------
              ! Vertical staggering :
              ! v(i,j-1,k+1)   V_LBC_Y2    v(i,j,k+1)        rho levels
              !
              !                V_LBC_YZ                      theta levels
              ! 
              ! v(i,j-1,k)     V_LBC_Y1    v(i,j,k)          rho levels
              !----------------------------------------------------------
        
              !----------------------------------------------------------
              ! Calculate LBC address for given i,j location
              !  for exner points, u points and v points
              !----------------------------------------------------------

                LBC_address(i,j) = LBC_START(PSouth) +                  &
                            (J + HALO_J-1) * lbc_row_len +              &
                             I + HALO_I-1

                LBC_address_u(i,j) = LBC_START_U(PSouth) +              &
                            (J + HALO_J-1) *lbc_row_len_u+              &
                             I + HALO_I-1

                LBC_address_v(i,j) = LBC_START_V(PSouth) +              &
                            (J + HALO_J-1) *lbc_row_len+                &
                             I + HALO_I-1
        
                !--------------------------------------------------------
                !
                ! Before horiz. interpolation must calculate
                ! location of i,j-1 points in LBC data (LBC_address_jm1)
                ! (required for horiz. interpolation of v).
                ! Applying BCs as appropriate.
                !
                !--------------------------------------------------------

                IF(J > first_row)THEN   
                  LBC_address_jm1(i,j) = LBC_START_V(PSouth) +          &
                               (J-1 + HALO_J-1) *(ROW_LENGTH + 2        &
                               * HALO_I) + I + HALO_I-1      
                ELSE
                  LBC_address_jm1(i,j) = LBC_START_V(PSouth) +          &
                               (J + HALO_J-1) *(ROW_LENGTH + 2          &
                               * HALO_I) + I + HALO_I-1 
                ENDIF ! j > first_row

                !--------------------------------------------------------
                ! Calculate LBC_ADDRESS for i-1,j points
                ! (required for horiz. interpolation of u)      
                ! Applying BCs as appropriate.
                !--------------------------------------------------------

                IF(I > row_start_pt)THEN
                  LBC_address_im1(i,j) = LBC_START_U(PSouth) +          &
                               (J + HALO_J-1) *lbc_row_len_u+           &
                                I-1 + HALO_I-1
                ENDIF 
             
                IF(I == row_start_pt)THEN
                  LBC_address_im1(i,j) = LBC_START_U(PSouth) +          &
                               (J + HALO_J-1) *lbc_row_len_u+           &
                                I + HALO_I-1
                ENDIF !I == row_start_pt
             
                IF(I == row_end_pt)THEN
                  LBC_address_u(i,j) = LBC_START_U(PSouth) +            &
                               (J + HALO_J-1) *lbc_row_len_u+           &
                                I-1 + HALO_I-1
                ENDIF !I == row_end_pt

               !---------------------------------------------------------
               ! Calculate vertical interpolation weights
               !---------------------------------------------------------

                IF(K < model_levels)THEN
                  weight = ( r_theta_levels(I,J,K) -                    &
                             r_rho_levels(I,J,K) ) /                    &
                           (r_rho_levels(I,J,K+1) - r_rho_levels(I,J,K))
                ELSE
                  weight = 0.0
                ENDIF !k < model_levels

                !--------------------------------------------------------
                ! Perform horizontal interpolation of both U and V 
                ! onto rho points for level k and for level k+1
                !--------------------------------------------------------
                
                U_LBC_X1 = 0.5 * U_LBC(LBC_address_im1(i,j),K) +        &
                           0.5 * U_LBC(LBC_address_u(i,j),K)

                IF( K < MODEL_LEVELS)THEN
                  U_LBC_X2 = 0.5 * U_LBC(LBC_address_im1(i,j),K+1) +    &
                             0.5 * U_LBC(LBC_address_u(i,j),K+1)
                ELSE
                  U_LBC_X2 = U_LBC_X1
                ENDIF

                V_LBC_Y1 = 0.5 * V_LBC(LBC_address_jm1(i,j),K) +        &
                           0.5 * V_LBC(LBC_address_v(i,j),K)

                IF( K < MODEL_LEVELS)THEN
                  V_LBC_Y2 = 0.5 * V_LBC(LBC_address_jm1(i,j),K+1) +    & 
                             0.5 * V_LBC(LBC_address_v(i,j),K+1)
                ELSE
                  V_LBC_Y2 = V_LBC_Y1
                ENDIF

                !--------------------------------------------------------
                ! Perform vertical interpolation onto theta points
                !--------------------------------------------------------
                
                U_LBC_XZ(i,j,k)=weight*U_LBC_X2 + (1.-weight)*U_LBC_X1

                V_LBC_YZ(i,j,k)=weight*V_LBC_Y2 + (1.-weight)*V_LBC_Y1

              ENDIF  !ns_weights(I,J) .NE. 0.0
            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

       ! NB. Broken into two separate loops to make vectorisable
        
        DO K = 1, model_levels

          DO J = first_row, last_row
            gj = datastart(2) + j - 1

            DO I = row_start_pt, row_end_pt
              gi = datastart(1) + i - 1

              IF (ns_weights(I,J) .NE. 0.0) THEN

                !--------------------------------------------------------
                ! Calculate new balanced values of Exner 
                !--------------------------------------------------------       

                temp = THETA_LBC(LBC_address(i,j),K) *                  &
                    (1.0 + C_Virtual * Q_LBC(LBC_address(i,j),K))

                IF ( K < model_levels ) THEN

                  EXNER_LBC(LBC_address(i,j),K+1) =                     &
                   EXNER_LBC(LBC_address(i,j),K)                        &
                   + ((r_rho_levels(I,J,K+1) - r_rho_levels(I,J,K)) /   &
                   (CP*temp)) *                                         &
                   (-g + f2_at_theta(gi,gj) * U_LBC_XZ(i,j,k) -         &
                   f1_at_theta(gi,gj) * V_LBC_YZ(i,j,k) +               &
                   ((U_LBC_XZ(i,j,k)**2 + V_LBC_YZ(i,j,k)**2)/          &
                   r_theta_levels(i,j,k)))

                ELSE ! K = model_levels
                  EXNER_LBC(LBC_address(i,j),K+1) =                     &
                   EXNER_LBC(LBC_address(i,j),K) +                      &
                   (2.0*(r_theta_levels(I,J,K)-r_rho_levels(I,J,K))/    & 
                   (CP*temp)) *                                         &
                   (-g +  f2_at_theta(gi,gj) * U_LBC_XZ(i,j,k) -        &
                   f1_at_theta(gi,gj) * V_LBC_YZ(i,j,k) +               &
                   ((U_LBC_XZ(i,j,k)**2 + V_LBC_YZ(i,j,k)**2)/          &
                   r_theta_levels(i,j,k)))                          
                
                ENDIF ! K < model_levels
                
                !--------------------------------------------------------       
                ! Calculate density to balance the new pressures
                ! using full non-linear equation of state
                !--------------------------------------------------------

                IF(K == 1) then
                  RHO_LBC(LBC_address(i,j),K) =                         &
                  (r_rho_levels(I,J,K) * r_rho_levels(I,J,K)) *         &
                  p_zero * EXNER_LBC(LBC_address(i,j),K) ** recip_Kappa/&
                  (R * THETA_LBC(LBC_address(i,j),K) *                  &
                  (1. +  C_Virtual * Q_LBC(LBC_address(i,j),K))         &   
                   * EXNER_LBC(LBC_address(i,j),K) )             
                ELSE ! K > 1
                  weight = ( r_rho_levels(I,J,K) -                      &
                   r_theta_levels(I,J,K-1) ) /                          &
                  (r_theta_levels(I,J,K) - r_theta_levels(I,J,K-1))

                  IF( K-1 > wet_levels ) THEN
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) +   &
                    (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1)
                  ELSEIF( K > wet_levels ) THEN
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) +   &
                     (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1) * &
                     (1.0 +  C_Virtual * Q_LBC(LBC_address(i,j),K-1) )
                  ELSE
                     temp =  weight  * THETA_LBC(LBC_address(i,j),K) *  &
                     (1.0 + C_Virtual * Q_LBC(LBC_address(i,j),K)) +    &
                     (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1) * &
                     (1.0 +  C_Virtual * Q_LBC(LBC_address(i,j),K-1) )
                  ENDIF  ! K-1 > wet_levels
              
                  RHO_LBC(LBC_address(i,j),K) =                         &
                  (r_rho_levels(I,J,K) * r_rho_levels(I,J,K)) *         &
                   p_zero*EXNER_LBC(LBC_address(i,j),K)** recip_Kappa / &
                  (R * temp * EXNER_LBC(LBC_address(i,j),K) )      
                ENDIF !  K == 1
        
               !---------------------------------------------------------
               ! Set vertical velocities to zero
               !---------------------------------------------------------

                W_LBC(LBC_address(i,j),K) = 0.0
                IF( .NOT. L_int_uvw_lbc )                               &
     &             W_ADV_LBC(LBC_address(i,j),K) = 0.0
           
              ENDIF  !ns_weights(I,J) .NE. 0.0
            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

      ENDIF ! IF (AT_EXTREMITY(PSouth))

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! 4.2 Northern LBC region
!------------------------------------------------------------------------

      IF (AT_EXTREMITY(PNorth)) THEN

        row_start_pt = 1 - HALO_I
        row_end_pt   = ROW_LENGTH + HALO_I
        first_row    = ROWS - RIMWIDTH + 1
        last_row     = ROWS + HALO_J

        lbc_row_len = ROW_LENGTH + 2 * HALO_I
       
        IF(AT_EXTREMITY(PEast))THEN
          lbc_row_len_u = ROW_LENGTH + 2 * HALO_I -1
          ! since LAM grid has one less u point than exner points
        
        ELSE
          lbc_row_len_u = ROW_LENGTH + 2 * HALO_I
        ENDIF

        IF (AT_EXTREMITY(PSouth)) THEN
          first_row_of_lbc = RIMWIDTH + 1
        ELSE ! Not at the South
          first_row_of_lbc = 1 - HALO_J
        ENDIF ! IF (AT_EXTREMITY(PSouth))

        DO K =1, model_levels
       
          DO J = first_row, last_row
            weights_J = ROWS + 1 - J
            gj = datastart(2) + j - 1

            DO I = row_start_pt, row_end_pt
              gi = datastart(1) + i - 1

              IF (ns_weights(I,weights_J) .NE. 0.0) THEN

              !----------------------------------------------------------
              ! u and v are required on theta points.
              !
              ! Step A: Horizontal interpolation onto exner points
              ! Step B: Vertical interpolation onto theta points
              !
              ! Horizontal grid staggering:
              !   exner(0,j+1)----u(0,j+1)---exner(1,j+1)---u(1,j+1)
              !       |
              !       |
              !     v(0,j)                    v(i,j)
              !       |
              !       |
              !   exner(0,j)      u(i-1,j)   exner(i,j)     u(i,j)
              !       |
              !       |
              !     v(0,j-1)                  v(i,j-1)
              !       |
              !         |
              !  
              !
              !                         
              !----------------------------------------------------------
              ! Vertical staggering :
              !
              !  u(i-1,j,k+1)  U_LBC_X2    u(i,j,k+1)        rho levels
              !
              !                U_LBC_XZ                      theta levels
              !
              !  u(i-1,j,k)    U_LBC_X1    u(i,j,k)          rho levels
              !----------------------------------------------------------
              !----------------------------------------------------------
              ! Vertical staggering :
              ! v(i,j-1,k+1)   V_LBC_Y2    v(i,j,k+1)        rho levels
              !
              !                V_LBC_YZ                      theta levels
              !
              ! v(i,j-1,k)     V_LBC_Y1    v(i,j,k)          rho levels
              !----------------------------------------------------------

              !----------------------------------------------------------
              ! Calculate LBC address for given i,j location
              !  for exner points, u points and v points
              !----------------------------------------------------------

                LBC_address(i,j) = LBC_START(PNorth) +                  &
                         (J - (ROWS - RIMWIDTH)-1) * lbc_row_len +      &
                         I + HALO_I - 1

                LBC_address_v(i,j) = LBC_START_V(PNorth) +              &
                         (J+1 - (ROWS - RIMWIDTH)-1) * lbc_row_len +    &
                         I + HALO_I - 1
                         
                !NB. J+1 above arises because first row of v data in LBCs
                !    is one row lower than the first row of exner/u data         

                LBC_address_u(i,j) = LBC_START_U(PNorth) +              &
                         (J - (ROWS - RIMWIDTH)-1) *lbc_row_len_u+      &
                         I + HALO_I - 1
        
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

                IF(I > row_start_pt)THEN
                  LBC_address_im1(i,j) = LBC_START_U(PNorth) +          &
                         (J - (ROWS - RIMWIDTH)-1) *lbc_row_len_u+      &
                         I-1 + HALO_I - 1
                ENDIF 
             
                IF(I == row_start_pt)THEN
                  LBC_address_im1(i,j) = LBC_START_U(PNorth) +          &
                         (J - (ROWS - RIMWIDTH)-1) *lbc_row_len_u+      &
                         I + HALO_I - 1
                ENDIF !I == row_start_pt
             
                IF(I == row_end_pt)THEN
                  LBC_address_u(i,j) = LBC_START_U(PNorth) +            &
                         (J - (ROWS - RIMWIDTH)-1) *lbc_row_len_u+      &
                         I-1 + HALO_I - 1
                ENDIF !I == row_end_pt

                !--------------------------------------------------------
                ! Calculate LBC_address for i,j-1 for 
                ! (required for interpolation of v)
                ! Applying BCs as appropriate
                !--------------------------------------------------------

                LBC_address_jm1(i,j) = LBC_START_V(PNorth) +            &
                              (J - (ROWS - RIMWIDTH)-1)*lbc_row_len+    &
                               I + HALO_I - 1     

               !NB. J not J-1 since first v point in North LBC
               !    is one point before first exner point.  

                IF (J == last_row)THEN
                  LBC_address_v(i,j) = LBC_START_V(PNorth) +            &
                       (J - (ROWS - RIMWIDTH)-1) * lbc_row_len +        &
                       I + HALO_I - 1
                ENDIF

               !---------------------------------------------------------
               ! Calculate vertical interpolation weights
               !---------------------------------------------------------

                IF(K < model_levels)THEN
                  weight = ( r_theta_levels(I,J,K) -                    &
                             r_rho_levels(I,J,K) ) /                    &
                           (r_rho_levels(I,J,K+1) - r_rho_levels(I,J,K))
                ELSE
                  weight = 0.0
                ENDIF !k < model_levels

                !--------------------------------------------------------
                ! Perform horizontal interpolation of both U and V 
                ! onto rho points for level k and for level k+1
                !--------------------------------------------------------
                
                U_LBC_X1 = 0.5 * U_LBC(LBC_address_im1(i,j),K) +        &  
                           0.5 * U_LBC(LBC_address_u(i,j),K)

                IF( K < MODEL_LEVELS)THEN
                  U_LBC_X2 = 0.5 * U_LBC(LBC_address_im1(i,j),K+1) +    &
                             0.5 * U_LBC(LBC_address_u(i,j),K+1)
                ELSE
                  U_LBC_X2 = 0.5 * U_LBC(LBC_address_im1(i,j),K) +      &
                             0.5 * U_LBC(LBC_address_u(i,j),K)
                ENDIF
                 
                V_LBC_Y1=0.5*V_LBC(LBC_address_jm1(i,j),K) +            &
                         V_LBC(LBC_address_v(i,j),K)*0.5

                IF( K < MODEL_LEVELS)THEN
                  V_LBC_Y2 = 0.5 * V_LBC(LBC_address_jm1(i,j),K+1) +    &
                             0.5 * V_LBC(LBC_address_v(i,j),K+1)
                ELSE
                  V_LBC_Y2 = 0.5 * V_LBC(LBC_address_jm1(i,j),K) +      &
                             0.5 * V_LBC(LBC_address_v(i,j),K)
                ENDIF

                !--------------------------------------------------------
                ! Perform vertical interpolation onto theta points
                !--------------------------------------------------------

                U_LBC_XZ(i,j,k)= weight*U_LBC_X2 + (1.-weight)*U_LBC_X1

                V_LBC_YZ(i,j,k)= weight*V_LBC_Y2 + (1.-weight)*V_LBC_Y1

              ENDIF   ! ns_weights(I,weights_J) .NE. 0.0
            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

       ! NB. Broken into two separate loops to make vectorisable
       
        DO K =1, model_levels
       
          DO J = first_row, last_row
            weights_J = ROWS + 1 - J
            gj = datastart(2) + j - 1

            DO I = row_start_pt, row_end_pt
              gi = datastart(1) + i - 1

              IF (ns_weights(I,weights_J) .NE. 0.0) THEN

              !----------------------------------------------------------
              ! Calculate new balanced values of Exner 
              !----------------------------------------------------------       

                temp = THETA_LBC(LBC_address(i,j),K) *                  &
                  (1.0 + C_Virtual * Q_LBC(LBC_address(i,j),K))           

                IF ( K < model_levels ) THEN
     
                  EXNER_LBC(LBC_address(i,j),K+1) =                     &
                        EXNER_LBC(LBC_address(i,j),K)                   &
                      + ((r_rho_levels(I,J,K+1) - r_rho_levels(I,J,K)) /&
                        (CP*temp)) *                                    &  
                       (-g + (f2_at_theta(gi,gj) * U_LBC_XZ(i,j,k)) -   &
                       (f1_at_theta(gi,gj) * V_LBC_YZ(i,j,k)) +         &  
                       ((U_LBC_XZ(i,j,k)**2 + V_LBC_YZ(i,j,k)**2)/      &
                       r_theta_levels(i,j,k)))

                ELSE ! K = model_levels
                
                  EXNER_LBC(LBC_address(i,j),K+1) =                     &
                    EXNER_LBC(LBC_address(i,j),K) +                     &
                    (2.0*(r_theta_levels(I,J,K)-r_rho_levels(I,J,K)) /  &
                    (CP*temp)) *                                        &
                    (-g +  (f2_at_theta(gi,gj) * U_LBC_XZ(i,j,k)) -     &
                    (f1_at_theta(gi,gj) * V_LBC_YZ(i,j,k)) +            &
                    ((U_LBC_XZ(i,j,k)**2 + V_LBC_YZ(i,j,k)**2)/         &
                    r_theta_levels(i,j,k)))      

                ENDIF ! K < model_levels        
    
                !--------------------------------------------------------       
                ! Calculate density to balance the new pressures
                ! using full non-linear equation of state
                !--------------------------------------------------------
            
                IF(K == 1) then
                  RHO_LBC(LBC_address(i,j),K) =                         &
                        (r_rho_levels(I,J,K) * r_rho_levels(I,J,K)) *   &
                      p_zero*EXNER_LBC(LBC_address(i,j),K)**recip_Kappa/&
                     (R * THETA_LBC(LBC_address(i,j),K) *               &
                     (1. +  C_Virtual * Q_LBC(LBC_address(i,j),K))      &
                     * EXNER_LBC(LBC_address(i,j),K) )
                ELSE ! K > 1
                  weight = ( r_rho_levels(I,J,K) -                      &
                             r_theta_levels(I,J,K-1) ) /                & 
                       (r_theta_levels(I,J,K) - r_theta_levels(I,J,K-1))

                  IF( K-1 > wet_levels )THEN
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) +   &
                      (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1) 
                  ELSEIF( K > wet_levels )THEN
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) +   &
                     (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1) * &
                     (1.0 +  C_Virtual * Q_LBC(LBC_address(i,j),K-1) )
                  ELSE
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) *   &
                       (1.0 + C_Virtual * Q_LBC(LBC_address(i,j),K)) +  &
                       (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1)*&
                       (1.0 +  C_Virtual * Q_LBC(LBC_address(i,j),K-1) )   
                  ENDIF  ! K-1 > wet_levels
                
                  RHO_LBC(LBC_address(i,j),K) =                         &    
                  (r_rho_levels(I,J,K) * r_rho_levels(I,J,K)) *         &
                  p_zero * EXNER_LBC(LBC_address(i,j),K)**recip_Kappa/  &
                  (R * temp * EXNER_LBC(LBC_address(i,j),K) )
                ENDIF !  K == 1
            
                !--------------------------------------------------------
                ! Set vertical velocities to zero
                !--------------------------------------------------------

                W_LBC(LBC_address(i,j),K) = 0.0
                IF( .NOT. L_int_uvw_lbc )                               &
     &             W_ADV_LBC(LBC_address(i,j),K) = 0.0
           
              ENDIF   ! ns_weights(I,weights_J) .NE. 0.0
            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

      ENDIF !  IF (AT_EXTREMITY(PNorth))

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! 4.3 Western LBC region
!------------------------------------------------------------------------

      IF (AT_EXTREMITY(PWest)) THEN

        row_start_pt = 1-HALO_I
        row_end_pt   = N_RIMS_TO_DO

        IF (AT_EXTREMITY(PSouth)) THEN
          first_row = RIMWIDTH+1
          first_row_of_LBC = RIMWIDTH+1
        ELSE ! Not at the South
          first_row_of_LBC = 1-HALO_J
          first_row = 1-HALO_J
        ENDIF ! IF (AT_EXTREMITY(PSouth))

        IF (AT_EXTREMITY(PNorth)) THEN
          last_row = ROWS - RIMWIDTH
        ELSE ! Not at the North
          last_row = ROWS + HALO_J
        ENDIF ! IF (AT_EXTREMITY(PNorth))

        lbc_row_len = HALO_I + RIMWIDTH

        !NB. Loop order k,i,j to improve speed

        DO K = 1, model_levels
       
          DO I = row_start_pt, row_end_pt
            gi = datastart(1) + i - 1

            DO J = first_row, last_row
              gj = datastart(2) + j - 1

              IF (ew_weights(I,J) .NE. 0.0) THEN

              !----------------------------------------------------------
              ! u and v are required on theta points.
              !
              ! Step A: Horizontal interpolation onto exner points
              ! Step B: Vertical interpolation onto theta points
              !
              ! Horizontal grid staggering:
              !
              !       |
              !     v(0,1)                        v(i,j)
              !       |
              !       |
              !   exner(0,1)       u(i-1,j)     exner(i,j)     u(i,j)
              !       |
              !       |
              !     v(0,0)                       v(i,j-1)
              !       |
              !         |
              !   
              !                         
              !----------------------------------------------------------
              ! Vertical staggering :
              !
              !  u(i-1,j,k+1)  U_LBC_X2    u(i,j,k+1)        rho levels
              !
              !                U_LBC_XZ                      theta levels
              !
              !  u(i-1,j,k)    U_LBC_X1    u(i,j,k)          rho levels
              !----------------------------------------------------------
              !----------------------------------------------------------
              ! Vertical staggering :
              ! v(i,j-1,k+1)   V_LBC_Y2    v(i,j,k+1)        rho levels
              !
              !                V_LBC_YZ                      theta levels
              !
              ! v(i,j-1,k)     V_LBC_Y1    v(i,j,k)          rho levels
              !----------------------------------------------------------

              !----------------------------------------------------------
              ! Calculate LBC address for given i,j location
              !  for exner points, u points and v points
              !----------------------------------------------------------

                LBC_address(i,j) = LBC_START(PWest) +                   &
                         (J - first_row_of_LBC) * lbc_row_len +         &
                         I + HALO_I - 1

                LBC_address_u(i,j) = LBC_START_U(PWest) +               &
                           (J - first_row_of_LBC) * lbc_row_len +       &
                           I + HALO_I - 1

                LBC_address_v(i,j) = LBC_START_V(PWest) +               &
                           (J - first_row_of_LBC) * lbc_row_len +       &
                           I + HALO_I - 1

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
        
                IF(GJ == RIMWIDTH + 1)THEN
                  LBCarea = 3   ! South LBC
                ELSE IF(GJ == GLOBAL_ROWS-RIMWIDTH)THEN
                  LBCarea = 1   ! North LBC
                ELSE
                  LBCarea = 4   ! West LBC
                ENDIF

                !--------------------------------------------------------
                ! Calculate LBC_address_jm1(i,j) for  
                !  appropriate LBC region.
                ! Applying BCs as appropriate
                !--------------------------------------------------------
                
                IF(LBCarea==4)THEN  ! West LBC
                  IF(j.NE.1-HALO_J)THEN

                    LBC_address_jm1(i,j) = LBC_START_V(PWest) +         &
                              (J-1 - first_row_of_LBC) * (HALO_I +      &
                              RIMWIDTH) + I + HALO_I - 1
                  ELSE
                    LBC_address_jm1(i,j) = LBC_START_V(PWest) +         &
                               (J - first_row_of_LBC) * (HALO_I +       &
                               RIMWIDTH) + I + HALO_I - 1
                  ENDIF ! j.NE.1-HALO_J

                ENDIF ! LBCarea==4
        
                IF(LBCarea==3)THEN  ! South LBC 
                  LBC_address_jm1(i,j) = LBC_START_V(PSouth) +          &
                              (J-1 + HALO_J-1) *(ROW_LENGTH + 2         &
                              * HALO_I) + I + HALO_I-1
                ENDIF ! LBCarea==3

                IF(LBCarea==1)THEN  ! North LBC 
                  LBC_address_v(i,j) = LBC_START_V(PNorth) +            &
                               (J+1 - (ROWS - RIMWIDTH)-1)*             &
                               (ROW_LENGTH + 2*HALO_I) +                &
                                I + HALO_I - 1

                  LBC_address_jm1(i,j) = LBC_START_V(PWest) +           &
                               (J-1 - first_row_of_LBC) * (HALO_I +     &
                               RIMWIDTH) + I + HALO_I - 1

                ENDIF ! LBCarea==1
 
                !--------------------------------------------------------
                ! Calculate LBC_ADDRESS for i-1,j points
                ! (required for horiz. interpolation of u)      
                !--------------------------------------------------------

                IF(I > row_start_pt)THEN
                  LBC_address_im1(i,j) = LBC_START_U(PWest) +           &
                              (J - first_row_of_LBC) * lbc_row_len +    &
                              I-1 + HALO_I - 1
                ELSE
                  LBC_address_im1(i,j) = LBC_START_U(PWest) +           &
                              (J - first_row_of_LBC) * lbc_row_len +    &
                              I + HALO_I - 1    
                ENDIF ! row_start_pt

               !---------------------------------------------------------
               ! Calculate vertical interpolation weights
               !---------------------------------------------------------

                IF(K < model_levels)THEN
                  weight = ( r_theta_levels(I,J,K) -                    &
                             r_rho_levels(I,J,K) ) /                    &
                           (r_rho_levels(I,J,K+1) - r_rho_levels(I,J,K))
                ELSE
                  weight = 0.0
                ENDIF !k < model_levels

                !--------------------------------------------------------
                ! Perform horizontal interpolation of both U and V 
                ! onto rho points for level k and for level k+1
                !--------------------------------------------------------
                
                U_LBC_X1 = 0.5 * U_LBC(LBC_address_im1(i,j),K) +        &
                           0.5 * U_LBC(LBC_address_u(i,j),K)
     
                IF( K < MODEL_LEVELS)THEN
                  U_LBC_X2 = 0.5 * U_LBC(LBC_address_im1(i,j),K+1) +    &
                             0.5 * U_LBC(LBC_address_u(i,j),K+1)
                ELSE
                  U_LBC_X2 = U_LBC_X1
                ENDIF
        
                V_LBC_Y1 = 0.5 * V_LBC(LBC_address_jm1(i,j),K) +        &
                           0.5 * V_LBC(LBC_address_v(i,j),K)

                IF( K < MODEL_LEVELS)THEN
                  V_LBC_Y2 = 0.5 * V_LBC(LBC_address_jm1(i,j),K+1) +    &
                             0.5 * V_LBC(LBC_address_v(i,j),K+1)
                ELSE
                  V_LBC_Y2 = V_LBC_Y1
                ENDIF

                !--------------------------------------------------------
                ! Perform vertical interpolation onto theta points
                !--------------------------------------------------------

                U_LBC_XZ(i,j,k)=weight*U_LBC_X2 + (1.-weight)*U_LBC_X1

                V_LBC_YZ(i,j,k)=weight*V_LBC_Y2 + (1.-weight)*V_LBC_Y1

              ENDIF  !ew_weights(I,J) .NE. 0.0
        
            ENDDO ! J
          ENDDO ! I
        ENDDO ! K

        ! NB. Broken into two separate loops to make vectorisable

        DO K = 1, model_levels
       
          DO I = row_start_pt, row_end_pt
            gi = datastart(1) + i - 1
        
            DO J = first_row, last_row
              gj = datastart(2) + j - 1

              IF (ew_weights(I,J) .NE. 0.0) THEN

              !----------------------------------------------------------
              ! Calculate new balanced values of Exner 
              !----------------------------------------------------------

                temp = THETA_LBC(LBC_address(i,j),K) *                  &
                      (1.0 + C_Virtual * Q_LBC(LBC_address(i,j),K))
                                
                IF ( K < model_levels ) THEN

                  EXNER_LBC(LBC_address(i,j),K+1) =                     &
                   EXNER_LBC(LBC_address(i,j),K)                        &
                   + ((r_rho_levels(I,J,K+1) - r_rho_levels(I,J,K)) /   &
                   (CP*temp)) *                                         &
                   (-g + f2_at_theta(gi,gj) * U_LBC_XZ(i,j,k) -         &
                   f1_at_theta(gi,gj) * V_LBC_YZ(i,j,k) +               &
                   ((U_LBC_XZ(i,j,k)**2 + V_LBC_YZ(i,j,k)**2)/          &
                    r_theta_levels(i,j,k)))

                ELSE ! K = model_levels
                  EXNER_LBC(LBC_address(i,j),K+1) =                     &
                    EXNER_LBC(LBC_address(i,j),K) +                     &
                    (2.0*(r_theta_levels(I,J,K)-r_rho_levels(I,J,K)) /  &
                    (CP*temp)) *                                        &
                    (-g +  f2_at_theta(gi,gj) * U_LBC_XZ(i,j,k) -       &
                    f1_at_theta(gi,gj) * V_LBC_YZ(i,j,k) +              &
                    ((U_LBC_XZ(i,j,k)**2 + V_LBC_YZ(i,j,k)**2)/         &
                    r_theta_levels(i,j,k))) 

                ENDIF ! K < model_levels

                !--------------------------------------------------------       
                ! Calculate density to balance the new pressures
                ! using full non-linear equation of state
                !--------------------------------------------------------
          
                IF(K == 1) then
                  RHO_LBC(LBC_address(i,j),K) =                         &
                  (r_rho_levels(I,J,K) * r_rho_levels(I,J,K)) *         &
                   p_zero * EXNER_LBC(LBC_address(i,j),K)**recip_Kappa/ &
                  (R * THETA_LBC(LBC_address(i,j),K) *                  &
                  (1. +  C_Virtual * Q_LBC(LBC_address(i,j),K))         &
                  * EXNER_LBC(LBC_address(i,j),K) )                    
                ELSE ! K > 1
                  weight = ( r_rho_levels(I,J,K) -                      &
                             r_theta_levels(I,J,K-1) ) /                & 
                       (r_theta_levels(I,J,K) - r_theta_levels(I,J,K-1))        

                  IF( K-1 > wet_levels )THEN
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) +   &
                         (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1)
                  ELSEIF( K > wet_levels )THEN
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) +   & 
                    (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1) *  &
                    (1.0 +  C_Virtual * Q_LBC(LBC_address(i,j),K-1) )
                  ELSE 
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) *   &
                       (1.0 + C_Virtual * Q_LBC(LBC_address(i,j),K)) +  &
                       (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1)*&
                       (1.0 +  C_Virtual * Q_LBC(LBC_address(i,j),K-1) )
                  ENDIF ! K-1 > wet_levels
             
                  RHO_LBC(LBC_address(i,j),K) =                         &
                           (r_rho_levels(I,J,K) * r_rho_levels(I,J,K)) *&
                    p_zero * EXNER_LBC(LBC_address(i,j),K)**recip_Kappa/&
                   (R * temp * EXNER_LBC(LBC_address(i,j),K) )             
                ENDIF !  K == 1  

                !--------------------------------------------------------
                ! Set vertical velocities to zero
                !--------------------------------------------------------
           
                W_LBC(LBC_address(i,j),K) = 0.0
                IF( .NOT. L_int_uvw_lbc )                               &
     &             W_ADV_LBC(LBC_address(i,j),K) = 0.0

              ENDIF  !ew_weights(I,J) .NE. 0.0
        
            ENDDO ! J
          ENDDO ! I
        ENDDO ! K
      
      ENDIF ! IF (AT_EXTREMITY(PWest))

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! 4.3 Eastern LBC region
!------------------------------------------------------------------------

      IF (AT_EXTREMITY(PEast)) THEN

        row_start_pt = ROW_LENGTH - N_RIMS_TO_DO + 1
        row_end_pt   = ROW_LENGTH + HALO_I

        IF (AT_EXTREMITY(PSouth)) THEN
          first_row = RIMWIDTH + 1
          first_row_of_LBC = RIMWIDTH + 1
        ELSE ! Not at the South
          first_row_of_LBC = 1 - HALO_J
          first_row = 1 - HALO_J
        ENDIF ! IF (AT_EXTREMITY(PSouth))

        IF (AT_EXTREMITY(PNorth)) THEN
          last_row = ROWS - RIMWIDTH
        ELSE ! Not at the North
          last_row = ROWS + HALO_J
        ENDIF ! IF (AT_EXTREMITY(PNorth))

        lbc_row_len = HALO_I + RIMWIDTH
        first_pt_of_LBC = ROW_LENGTH - RIMWIDTH + 1

        !NB. Loop order k,i,j to improve speed

        DO K = 1, model_levels
        
          DO I = row_start_pt, row_end_pt
            gi = datastart(1) + i - 1
            weights_I = first_pt_of_LBC + RIMWIDTH - I

            DO J = first_row, last_row
              gj = datastart(2) + j - 1

              IF (ew_weights(weights_I,J) .NE. 0.0) THEN

              !----------------------------------------------------------
              ! u and v are required on theta points.
              !
              ! Step A: Horizontal interpolation onto exner points
              ! Step B: Vertical interpolation onto theta points
              !
              ! Horizontal grid staggering:
              !                                               |
              !                                               |
              !                   v(i,j)                   v(i+1,j)
              !                                         |
              !                                               |
              !    u(i-1,j)     exner(i,j)     u(i,j)     exner(i+1,j)
              !                                               |
              !                                               |
              !                   v(i,j-1)                 v(i+1,j-1)
              !                                               |
              !                                                 |
              !            
              !
              !                         
              !----------------------------------------------------------
              ! Vertical staggering :
              !
              !  u(i-1,j,k+1)  U_LBC_X2    u(i,j,k+1)        rho levels
              !
              !                U_LBC_XZ                      theta levels
              !
              !  u(i-1,j,k)    U_LBC_X1    u(i,j,k)          rho levels
              !----------------------------------------------------------
              !----------------------------------------------------------
              ! Vertical staggering :
              ! v(i,j-1,k+1)   V_LBC_Y2    v(i,j,k+1)        rho levels
              !
              !                V_LBC_YZ                      theta levels
              !
              ! v(i,j-1,k)     V_LBC_Y1    v(i,j,k)          rho levels
              !----------------------------------------------------------
              !----------------------------------------------------------
              ! Calculate LBC address for given i,j location
              !  for exner points, u points and v points
              !----------------------------------------------------------

                LBC_address(i,j) = LBC_START(PEast)+                    &
                         (J - first_row_of_LBC) * lbc_row_len +         &
                          I - first_pt_of_LBC

                LBC_address_v(i,j) = LBC_START_V(PEast)+                &
                           (J - first_row_of_LBC) * lbc_row_len +       &
                            I - first_pt_of_LBC

                LBC_address_u(i,j) = LBC_START_U(PEast)+                &  
                            (J - first_row_of_LBC) * lbc_row_len +      &
                            I+1 - first_pt_of_LBC

                !NB. I+1 arises because first u point in East LBC
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

                IF(I == row_end_pt)THEN
                  LBC_address_u(i,j) = LBC_START_U(PEast)+              &
                            (J - first_row_of_LBC) * lbc_row_len +      &
                            I - first_pt_of_LBC
                ENDIF
        
                IF(GJ == RIMWIDTH+1)THEN
                  LBCarea = 3 ! South
                ELSE IF(GJ == GLOBAL_ROWS - RIMWIDTH)THEN
                  LBCarea = 1 ! North
                ELSE 
                  LBCarea = 2 ! East
                ENDIF

                !--------------------------------------------------------
                ! Calculate LBC_address_jm1 appropriate LBC region
                ! Applying BCs as appropriate
                !--------------------------------------------------------
                
                IF(LBCarea==2)THEN
                  IF(J.NE.1-HALO_J)THEN

                    LBC_address_jm1(i,j) = LBC_START_V(PEast)+          &
                               (J-1 - first_row_of_LBC) * (HALO_I +     &
                               RIMWIDTH) + I - first_pt_of_LBC       
                  ELSE
                    LBC_address_jm1(i,j) = LBC_START_V(PEast)+          &
                               (J - first_row_of_LBC) * (HALO_I +       &
                               RIMWIDTH) + I - first_pt_of_LBC
                  ENDIF !J.NE.1-HALO_J
                ENDIF ! LBCarea==2
        
                IF(LBCarea==3)THEN      
                  LBC_address_jm1(i,j) = LBC_START_V(PSouth) +          &
                              (J-1 + HALO_J-1) *(ROW_LENGTH + 2         &
                              * HALO_I) + I + HALO_I-1
                ENDIF ! LBCarea==3
           
                IF(LBCarea==1)THEN
                  LBC_address_v(i,j) = LBC_START_V(PNorth) +            &
                               (J+1 - (ROWS - RIMWIDTH)-1)*             &
                              (ROW_LENGTH + 2*HALO_I) +                 &
                               I + HALO_I - 1

                  LBC_address_jm1(i,j) = LBC_START_V(PEast)+            &
                           (J-1 - first_row_of_LBC) * lbc_row_len +     &
                           I - first_pt_of_LBC
                ENDIF ! LBCarea==1

                !--------------------------------------------------------
                ! Calculate LBC_ADDRESS for i-1,j points
                ! (required for horiz. interpolation of u)      
                !--------------------------------------------------------

                LBC_address_im1(i,j) = LBC_START_U(PEast)+              &
                              (J - first_row_of_LBC) * lbc_row_len +    &
                              I - first_pt_of_LBC

                !NB. I not I-1 since first u point in East LBC
                !    is one point before first exner point.  

                !--------------------------------------------------------
                ! Calculate vertical interpolation weights
                !--------------------------------------------------------

                IF(K < model_levels)THEN
                  weight = ( r_theta_levels(I,J,K) -                    &
                             r_rho_levels(I,J,K) ) /                    &
                           (r_rho_levels(I,J,K+1) - r_rho_levels(I,J,K))
                ELSE
                  weight = 0.0
                ENDIF !k < model_levels

                !--------------------------------------------------------
                ! Perform horizontal interpolation of both U and V 
                ! onto rho points for level k and for level k+1
                !--------------------------------------------------------
                
                U_LBC_X1 = 0.5 * U_LBC(LBC_address_im1(i,j),K) +        &
                           0.5 * U_LBC(LBC_address_u(i,j),K)

                IF( K < MODEL_LEVELS)THEN
                  U_LBC_X2 = 0.5 * U_LBC(LBC_address_im1(i,j),K+1) +    &
                             0.5 * U_LBC(LBC_address_u(i,j),K+1)
                ELSE
                  U_LBC_X2=U_LBC_X1
                ENDIF
                
                V_LBC_Y1 = 0.5 * V_LBC(LBC_address_jm1(i,j),K) +        &
                           0.5 * V_LBC(LBC_address_v(i,j),K)
 
                IF( K < MODEL_LEVELS)THEN
                  V_LBC_Y2 = 0.5 * V_LBC(LBC_address_jm1(i,j),K+1) +    &
                             0.5 * V_LBC(LBC_address_v(i,j),K+1)
                ELSE
                  V_LBC_Y2 = V_LBC_Y1
                ENDIF 

                !--------------------------------------------------------
                ! Perform vertical interpolation onto theta points
                !--------------------------------------------------------

                U_LBC_XZ(i,j,k)= weight*U_LBC_X2 + (1.-weight)*U_LBC_X1

                V_LBC_YZ(i,j,k)= weight*V_LBC_Y2 + (1.-weight)*V_LBC_Y1

              ENDIF  !ew_weights(weights_I,J) .NE. 0.0
       
            ENDDO ! J
          ENDDO ! I
        ENDDO ! K

        ! NB. Broken into two separate loops to make vectorisable

        DO K = 1, model_levels
        
          DO I = row_start_pt, row_end_pt
            gi = datastart(1) + i - 1
            weights_I = first_pt_of_LBC + RIMWIDTH - I

            DO J = first_row, last_row
              gj = datastart(2) + j - 1

              IF (ew_weights(weights_I,J) .NE. 0.0) THEN

              !----------------------------------------------------------
              ! Calculate new balanced values of Exner 
              !----------------------------------------------------------       
                
                temp = THETA_LBC(LBC_address(i,j),K) *                  &
                       (1.0 + C_Virtual * Q_LBC(LBC_address(i,j),K))
                        
                IF ( K < model_levels ) THEN

                  EXNER_LBC(LBC_address(i,j),K+1) =                     &
                                EXNER_LBC(LBC_address(i,j),K)           &    
                     + ((r_rho_levels(I,J,K+1) - r_rho_levels(I,J,K)) / &
                       (CP*temp)) *                                     &
                       (-g + f2_at_theta(gi,gj) * U_LBC_XZ(i,j,k) -     &
                       f1_at_theta(gi,gj) * V_LBC_YZ(i,j,k) +           &
                       ((U_LBC_XZ(i,j,k)**2 + V_LBC_YZ(i,j,k)**2)/      &
                       r_theta_levels(i,j,k)))                   

                ELSE ! K = model_levels
                  EXNER_LBC(LBC_address(i,j),K+1) =                     &
                    EXNER_LBC(LBC_address(i,j),K) +                     &
                    (2.0*(r_theta_levels(I,J,K)-r_rho_levels(I,J,K)) /  &
                    (CP*temp)) *                                        &
                    (-g +  f2_at_theta(gi,gj) * U_LBC_XZ(i,j,k) -       &
                    f1_at_theta(gi,gj) * V_LBC_YZ(i,j,k) +              &
                    ((U_LBC_XZ(i,j,k)**2 + V_LBC_YZ(i,j,k)**2)/         &
                    r_theta_levels(i,j,k)))                

                ENDIF ! K < model_levels

                !--------------------------------------------------------       
                ! Calculate density to balance the new pressures
                ! using full non-linear equation of state
                !--------------------------------------------------------
           
                IF(K == 1) then
                  RHO_LBC(LBC_address(i,j),K) =                         &
                   (r_rho_levels(I,J,K) * r_rho_levels(I,J,K)) *        &
                    p_zero * EXNER_LBC(LBC_address(i,j),K)**recip_Kappa/&
                   (R * THETA_LBC(LBC_address(i,j),K) *                 &
                   (1. +  C_Virtual * Q_LBC(LBC_address(i,j),K))        &
                    * EXNER_LBC(LBC_address(i,j),K) )
                ELSE ! K > 1
                  weight = ( r_rho_levels(I,J,K) -                      &
                                  r_theta_levels(I,J,K-1) ) /           &
                       (r_theta_levels(I,J,K) - r_theta_levels(I,J,K-1)) 

                  IF( K-1 > wet_levels )THEN
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) +   &
                       (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1)      
                  ELSEIF( K > wet_levels )THEN
                     temp =  weight  * THETA_LBC(LBC_address(i,j),K) +  &
                       (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1)*&
                       (1.0 +  C_Virtual * Q_LBC(LBC_address(i,j),K-1) )    
                  ELSE
                    temp =  weight  * THETA_LBC(LBC_address(i,j),K) *   &
                       (1.0 + C_Virtual * Q_LBC(LBC_address(i,j),K)) +  &
                       (1.0 - weight) * THETA_LBC(LBC_address(i,j),K-1)*&
                       (1.0 +  C_Virtual * Q_LBC(LBC_address(i,j),K-1) )
                  ENDIF  ! K-1 > wet_levels
            
                  RHO_LBC(LBC_address(i,j),K) =                         &
                    (r_rho_levels(I,J,K) * r_rho_levels(I,J,K)) *       &
                    p_zero * EXNER_LBC(LBC_address(i,j),K)**recip_Kappa/&
                    (R * temp * EXNER_LBC(LBC_address(i,j),K) )          
                ENDIF !  K == 1

                !--------------------------------------------------------
                ! Set vertical velocities to zero
                !--------------------------------------------------------
          
                W_LBC(LBC_address(i,j),K) = 0.0
                IF( .NOT. L_int_uvw_lbc )                               &
     &             W_ADV_LBC(LBC_address(i,j),K) = 0.0
          
              ENDIF  !ew_weights(weights_I,J) .NE. 0.0
       
            ENDDO ! J
          ENDDO ! I
        ENDDO ! K

      ENDIF ! IF (AT_EXTREMITY(PEast))

       IF (lhook) CALL dr_hook('BALANCE_LBC_VALUES',zhook_out,zhook_handle)
       RETURN
     END SUBROUTINE BALANCE_LBC_VALUES
