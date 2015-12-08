! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_set_initial_humidity

      Subroutine idl_set_initial_humidity(                              &
                            height_domain                               &
      ,                     row_length, rows                            &
      ,                     model_levels, wet_model_levels              &
      ,                     me, halo_i, halo_j                          &
      ,                     qprofile_number, L_dry, q1                  &
      ,                     max_num_profile_data, num_profile_data      &
      ,                     zprofile_data, qprofile_data                &
      ,                     zprofile_orog, idl_interp_option, hf        &
      ,                     r_theta_levels, eta_theta_levels            &
      ,                     q, theta, exner_theta_levels                &
      ,                     q_ref, theta_ref, exner_ref                 &
      ,                     L_code_test)

! Purpose:
!          Sets up initial data for idealised problems.
!          Input from idealised namelist.
!
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + extensions
!   This code is written to UMDP3 programming standards.
!

      USE earth_constants_mod, ONLY: earth_radius
      USE atmos_constants_mod, ONLY: p_zero, recip_kappa

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE qprofile_mod, ONLY: qp_dry, qp_qsat, qp_namelist_rh,          &
                              qp_namelist, qp_dump
      IMPLICIT NONE

! Arguments with Intent (In)

      !  Physical constants
      Real, Intent(In) :: height_domain       ! Height of top of domain

      ! Grid dimensions
      INTEGER, INTENT(IN) :: row_length       ! No. of points on a row
      INTEGER, INTENT(IN) :: rows             ! No. of rows (theta)
      INTEGER, INTENT(IN) :: model_levels     ! No. of model levels
      INTEGER, INTENT(IN) :: wet_model_levels ! No. of wet model levels
      INTEGER, INTENT(IN) :: halo_i           ! Size of halo in x
      INTEGER, INTENT(IN) :: halo_j           ! Size of halo in y

      ! Multi-processor
      INTEGER, INTENT(IN) :: me               ! My processor number

      ! Idealised options
      LOGICAL, INTENT(IN) :: L_code_test      ! User switch
      LOGICAL, INTENT(IN) :: L_dry            ! Dry model switch

      INTEGER, INTENT(IN) :: max_num_profile_data ! max no. profile data
      INTEGER, INTENT(IN) :: num_profile_data ! actual no. data values
      INTEGER, INTENT(IN) :: qprofile_number  ! Humidity profile option
      INTEGER, INTENT(IN) :: idl_interp_option ! Profile interp option

      REAL,    INTENT(IN) :: q1               ! For constant RH
      REAL,    INTENT(IN) :: zprofile_orog    ! Orog hgt of initial prof
      REAL,    INTENT(IN) :: hf               ! Hgt for constant interp

      ! Idealised namelist data profile
      REAL, INTENT(IN) :: zprofile_data(max_num_profile_data) ! Heights
      REAL, INTENT(IN) :: qprofile_data(max_num_profile_data) ! Humidity

      ! Vertical co-ordinate information
      REAL, INTENT(IN) :: r_theta_levels(1-halo_i:row_length+halo_i,    &
                             1-halo_j:rows+halo_j,0:model_levels)
      REAL, INTENT(IN) :: eta_theta_levels(0:model_levels)

! Arguments with Intent (InOut)

      ! 3D prognostic field arrays
      REAL, INTENT(INOUT) ::                                            &
        theta(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,          &
                                         model_levels)                  &
      , exner_theta_levels(1-halo_i:row_length+halo_i,                  &
                           1-halo_j:rows+halo_j, model_levels)          &
      , q(1-halo_i:row_length+halo_i,                                   &
                1-halo_j:rows+halo_j, wet_model_levels)

      ! Vertical reference profile
      REAL, INTENT(INOUT) :: theta_ref(model_levels)
      REAL, INTENT(INOUT) :: exner_ref(model_levels+1)
      REAL, INTENT(INOUT) :: q_ref(wet_model_levels)

! Local variables

      ! Work arrays
      REAL :: p_row(1-halo_i:row_length+halo_i)  ! pressure on a row
      REAL :: t_row(1-halo_i:row_length+halo_i)  ! temperature on a row
      REAL :: qs_row(1-halo_i:row_length+halo_i) ! qsat on a row
      REAL :: p_ref(model_levels)      ! reference pressure profile
      REAL :: t_ref(model_levels+1)    ! reference temperature profile
      REAL :: qs_ref(wet_model_levels) ! reference saturation profile
      REAL :: weight                   ! interpolation weight
      REAL :: z_at_theta               ! height asl of theta level
      REAL :: z_at_orog                ! height of orog
      REAL :: hs                       ! Temporary variable
      REAL :: eta_model                ! eta coord for model levels
      REAL :: eta_profile(num_profile_data) ! eta coord for input prof

      INTEGER :: i, j, k, k2           ! Loop counters
      INTEGER :: length                ! array length

      ! Error reporting
      CHARACTER (LEN=*),  PARAMETER :: RoutineName=                     &
     &                                 'idl_set_initial_humidity'
      CHARACTER (LEN=256)           :: Cmessage
      INTEGER                       :: ErrorStatus

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_SET_INITIAL_HUMIDITY',zhook_in,zhook_handle)
!     IF (me == 0) THEN
!       WRITE (6,*) ' '
!       WRITE (6,*) ' HUMIDITY PROFILE  '
!     END IF

!-------------------------------------------------------------------
!
!                      Dry. Humidity = 0
!
!-------------------------------------------------------------------

      IF (L_dry .or. (qprofile_number == qp_dry)) THEN

        ! set moisture fields to zero.
!       IF (me == 0) THEN
!         Write (6,*) '   Dry model run, q = 0.0'
!       END IF   !(me == 0)

        q(:,:,:) = 0.0
        q_ref(:) = 0.0

!-------------------------------------------------------------------
!
!                  Set constant relative humidity
!
!-------------------------------------------------------------------
      ELSE IF (qprofile_number == qp_qsat) THEN

!       IF (me == 0) THEN
!         WRITE (6,*) '   Setting constant relative humidity',          &
!                     ' ',q1,'% wrt water T>0degC, wrt ice T<0degC'
!       END IF

        ! set q field to constant relative humidity
        length = row_length + halo_i + halo_i
        DO k = 1, wet_model_levels
          DO j = 1-halo_j, rows+halo_j

            ! Calculate temperature and pressure for QSAT
            DO i = 1-halo_i, row_length+halo_i
              t_row(i)=theta(i,j,k) * exner_theta_levels(i,j,k)
              p_row(i)=p_zero*exner_theta_levels(i,j,k)**recip_Kappa
            END DO

! DEPENDS ON: qsat
            CALL qsat( qs_row(1-halo_i), t_row(1-halo_i),               &
                        p_row(1-halo_i), length)

            DO i = 1-halo_i, row_length+halo_i
              q(i,j,k) = q1 * qs_row(i)
            END DO

          END DO
        END DO

        ! Set reference profile
        DO k = 1, wet_model_levels
          t_ref(k) = theta_ref(k) * exner_ref(k)
          p_ref(k) = p_zero*exner_ref(k)**recip_Kappa
        END DO
        length = wet_model_levels
! DEPENDS ON: qsat
        CALL qsat( qs_ref, t_ref, p_ref, length)
        DO k = 1, wet_model_levels
          q_ref(k) = q1 * qs_ref(k)
        END DO

! Use tapered profile
        q = q/q1
        q_ref = q_ref/q1
        DO k = 1, wet_model_levels
          DO j = 1-halo_j, rows+halo_j
            DO i = 1-halo_i, row_length+halo_i
              IF ( r_theta_levels(i,j,k) - Earth_radius < 1000.0 ) THEN
                q(i,j,k) = 0.95*q(i,j,k)
              ELSE IF ( r_theta_levels(i,j,k) - Earth_radius < 20000.0 ) THEN
                q(i,j,k) = (0.95-4.95e-5                                      &
                          *(r_theta_levels(i,j,k) - Earth_radius-1000.0))     &
                          *q(i,j,k)
              ELSE        
                q(i,j,k) = 0.0
              END IF          
            END DO
          END DO
          IF ( r_theta_levels(1,1,k) - Earth_radius < 1000.0 ) THEN
            q_ref(k) = 0.95*q_ref(k)
          ELSE IF ( r_theta_levels(1,1,k) - Earth_radius < 20000.0 ) THEN
            q_ref(k) = (0.95-4.95e-5                                          &
                      *(r_theta_levels(1,1,k)- Earth_radius -1000.0))         &
                      *q_ref(k)
          ELSE        
            q_ref(k) = 0.0
          END IF  
        END DO
! WRITE(6,*) 'Moisture profile',q1
!i=1
!j=1
!DO k=1,wet_model_levels
!  WRITE(6,fmt='(I4,3E16.8)') k,r_theta_levels(i,j,k)- Earth_radius, &
!                             q(i,j,k),q_ref(k)
!END DO
!WRITE(6,*) ' '       

!-------------------------------------------------------------------
!
!            Set field from profile in idealised namelist
!
!-------------------------------------------------------------------
      ELSE IF (qprofile_number == qp_namelist .OR.                      &
               qprofile_number == qp_namelist_rh) THEN

!       IF (me == 0) THEN
!         WRITE (6,*) '   Setting humidity from namelist profile'
!       END IF

        ! Check to make sure the namelist profile data extends
        ! to the top of the model.
        DO j = 1-halo_j, rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            IF (r_theta_levels(i,j,model_levels) - Earth_radius         &
                 >   zprofile_data(num_profile_data)) THEN
              Cmessage =                                                &
                'Idealised namelist vertical profile data'              &
                //'does not extend to the top of the model.'            &
                //'Please modify the namelist data.'
              ErrorStatus = 1

              CALL Ereport( RoutineName, ErrorStatus, Cmessage )
            END IF
          END DO
        END DO

        ! Interpolate q from namelist profile to model levels
        DO k = 1, model_levels
          DO k2 = 1, num_profile_data-1
            DO j = 1-halo_j, rows+halo_j
              DO i = 1-halo_i, row_length+halo_i

                z_at_theta = r_theta_levels(i,j,k) - Earth_radius

                IF (z_at_theta  >   zprofile_data(k2) .AND.             &
                    z_at_theta  <=  zprofile_data(k2+1)) THEN

                  weight = (z_at_theta - zprofile_data(k2))             &
                         /(zprofile_data(k2+1) - zprofile_data(k2))

                  q(i,j,k) = qprofile_data(k2) + weight*                &
                         (qprofile_data(k2+1) - qprofile_data(k2))
                END IF

              END DO
            END DO
          END DO
        END DO


        !---------------------------------------------------------------
        ! Alternative interpolation of initial profile over orography
        !---------------------------------------------------------------
        !
        !  idl_interp_option = 1: constant on height levels
        !                         (default above, no need to modify)
        !  idl_interp_option = 2: hybrid height everywhere up to a
        !                         specified height hf (the same as modl)
        !                         levels are defined). hs=height_domain
        !  idl_interp_option = 3: as option 2 but only when orography is
        !                         less than input profile orography.
        !                         hs = zprofile_orog
        !
        ! Sets up eta coord for each level of the initial profile data
        ! and an eta coordinate for each model column, and interpolates
        ! in eta space if the model level height is less than "hf" and
        ! the model orography height is less than "hs"
        !---------------------------------------------------------------

        IF (idl_interp_option == 2 .OR. idl_interp_option == 3) THEN

          IF (idl_interp_option == 2) hs = height_domain
          IF (idl_interp_option == 3) hs = zprofile_orog

          ! Set up eta coord for each level of the initial profile data
          eta_profile(1) = 0.0
          DO k2 = 2, num_profile_data
            eta_profile(k2) =  (zprofile_data(k2) - zprofile_orog)      &
                              /(hf - zprofile_orog)
          END DO

          ! Interpolate in eta space and overwrite q where appropriate
          DO k = 1, wet_model_levels
            DO k2 = 1, num_profile_data - 1
              DO j = 1-halo_j, rows+halo_j
                DO i = 1-halo_i, row_length+halo_i
                  z_at_theta = r_theta_levels(i,j,k) - Earth_radius
                  z_at_orog  = r_theta_levels(i,j,0) - Earth_radius
                  eta_model = (z_at_theta - z_at_orog)/(hf - z_at_orog)

                  IF ( (z_at_orog <= hs ) .AND. (z_at_theta < hf) .AND. &
                       (eta_model > eta_profile(k2))  .AND.             &
                       (eta_model <= eta_profile(k2+1)) ) THEN

                    weight = (eta_model - eta_profile(k2))              &
                             /(eta_profile(k2+1) - eta_profile(k2))
                    q(i,j,k) = qprofile_data(k2) + weight*              &
                              (qprofile_data(k2+1) - qprofile_data(k2))
                  END IF
                END DO
              END DO
            END DO
          END DO

        END IF ! on idl_interp_option = 2 or 3

        !---------------------------------------------
        ! Set up reference q profile (assumes no orography)
        !---------------------------------------------
        DO k = 1, model_levels
          DO k2 = 1, num_profile_data-1

                z_at_theta = eta_theta_levels(k)*height_domain

                IF (z_at_theta  >   zprofile_data(k2) .AND.             &
                    z_at_theta  <=  zprofile_data(k2+1)) THEN

                  weight = (z_at_theta - zprofile_data(k2))             &
                         /(zprofile_data(k2+1) - zprofile_data(k2))

                  q_ref(k) = qprofile_data(k2) + weight*                &
                            (qprofile_data(k2+1) - qprofile_data(k2))
                END IF

          END DO
        END DO


        !---------------------------------------------
        ! If input data is relative humidity (%/100)
        ! then convert to q (kg/kg)
        !---------------------------------------------
        IF (qprofile_number == qp_namelist_rh) THEN

          ! set q field to constant relative humidity wrt water
!         IF (me == 0) THEN
!           Write (6,*) '   Relative humidity (wrt water) data '
!         END IF

          ! Calculate humidity in kg/kg
          ! q(kg/kg) = RH (held in variable q) * qsat

          length = row_length + halo_i + halo_i
          DO k = 1, wet_model_levels
            DO j = 1-halo_j, rows+halo_j

              ! Calculate temperature and pressure for QSAT
              DO i = 1-halo_i, row_length+halo_i
                t_row(i)=theta(i,j,k) * exner_theta_levels(i,j,k)
                p_row(i)=p_zero*exner_theta_levels(i,j,k)**recip_Kappa
              END DO

! DEPENDS ON: qsat_wat
              CALL qsat_wat( qs_row(1-halo_i), t_row(1-halo_i),         &
                              p_row(1-halo_i), length )

              DO i = 1-halo_i, row_length+halo_i
                q(i,j,k) = q(i,j,k) * qs_row(i)
              END DO

            END DO
          END DO

          ! Set reference profile
          DO k = 1, wet_model_levels
            t_ref(k) = theta_ref(k) * exner_ref(k)
            p_ref(k) = p_zero*exner_ref(k)**recip_Kappa
          END DO
          length = wet_model_levels
! DEPENDS ON: qsat
          CALL qsat( qs_ref, t_ref, p_ref, length)
          DO k = 1, wet_model_levels
            q_ref(k) = q_ref(k) * qs_ref(k)
          END DO

        END IF

!-------------------------------------------------------------------
!
!            Leave q exactly as it is (e.g. from initial dump)
!
!-------------------------------------------------------------------
      ELSE IF (qprofile_number == qp_dump) THEN

!       Write (6,*) '   Start dump 3D humidity field will be used'

!-------------------------------------------------------------------
!
!              qprofile_number not recognised !
!
!-------------------------------------------------------------------
      ELSE

        Cmessage =                                                      &
                'qprofile_number not recognised.'                       &
              //' Please modify the idealised namelist data.'
        ErrorStatus = 1

        CALL Ereport( RoutineName, ErrorStatus, Cmessage )

      END IF

      IF (lhook) CALL dr_hook('IDL_SET_INITIAL_HUMIDITY',zhook_out,zhook_handle)

      END SUBROUTINE idl_set_initial_humidity
