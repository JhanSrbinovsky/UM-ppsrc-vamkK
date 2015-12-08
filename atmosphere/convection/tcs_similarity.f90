! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Module to calculate similarity functions for warm rain tcs.
!
MODULE tcs_similarity


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE
  !
  ! Description:
  ! Module to calculate similarity functions for warm rain tcs.
  !
  ! Method:
  !   Calculate the similarity profiles using standard functions.
  !   Parameters for these functions are determined from a
  !   set predefined in tcs_parameters_warm based on the type of 
  !   convection that has been diagnosed (conv_type).
  !   For specific detail see:
  !   <reference to documentation to go here, once available>
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.1 programming standards.
  !


  ! Working note: This needs to be tidied up a little bit and to make
  ! sure that the <= and < refer exactly to in cloud for rho and theta
  ! levels.  

  ! Working note: Need to check whether these functions are on rho or
  ! theta levels and that they are used consistently in the rest of the
  ! code. 

  ! Working note: Might be nice to get rid of the mask if possible.
  ! (Just need to make sure bounds and initialisations are set
  ! correctly in the rest of the code)



CONTAINS

  SUBROUTINE calc_similarity(eta_theta, eta_rho, conv_type, sim)  

    USE tcs_parameters_warm,       ONLY :              &
       similarity_coefficient, sim_coeff_nonp_sh,      &
       sim_coeff_warm_cg

    USE tcs_class_similarity,      ONLY :              &
       similarity 

    USE tcs_common_warm,           ONLY :              &
       scales

    IMPLICIT NONE
    !-------------------------------------------------------------------
    ! Subroutine Arguments
    !-------------------------------------------------------------------
    REAL, INTENT(in) :: eta_theta(:,:)
                        ! non-dimensional height of cloud levels
    REAL, INTENT(in) :: eta_rho(:,:)  
                        ! non-dimensional height of cloud levels on rho levels
    INTEGER, INTENT(in) :: conv_type(:)
                        ! Indicator of type of convection 
                        !    1=non-precipitating shallow
                        !    2=drizzling shallow
                        !    3=warm congestus

    ! Sim is now declared as intent(inout) so that memory allocation
    ! (and deallocation) can be done in the calling routine.
    TYPE(similarity), INTENT(inout):: sim 

    !-----------------------------------------------------------------
    ! Variables defined locally
    ! Note that sim%n_xx and sim%nlev are dimensions of similarity
    ! functions and are used as shorthand instead of size(sim%...,n)
    !-----------------------------------------------------------------
    REAL :: wexpG(sim%n_xx,sim%nlev)  
    REAL :: wexpK(sim%n_xx,sim%nlev)  
    REAL :: wexpF(sim%n_xx,sim%nlev)  
    REAL :: wexpB(sim%n_xx,sim%nlev)  

    REAL :: arr1(sim%n_xx)   ! An array expression  
    REAL :: arr2(sim%n_xx)   ! An array expression  

    REAL :: zi(sim%n_xx,sim%nlev)     ! Shifted height

    INTEGER :: i,ii,k ! loop counter

    TYPE(similarity_coefficient), POINTER :: SC 
                        ! coefficients for similarity functions

    !
    ! Error reporting variables
    !
    INTEGER :: ErrorStatus

    CHARACTER (Len=80) :: Cmessage

    CHARACTER (Len=*), PARAMETER :: RoutineName='calc_similarity'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('TCS_SIMILARITY:CALC_SIMILARITY',zhook_in,zhook_handle)
    DO i=1,sim%n_xx

      ! Select the appropriate coefficients
      SELECT CASE (conv_type(i))
      CASE (1:2)
        ! shallow convection (NB eventually may want to be split 
        !          into precipitating(1) and non-precipitating(2))
        SC=>sim_coeff_nonp_sh
      CASE (3:4)
        ! congestus convection (NB eventually may want to be split 
        !                       into warm(3) and ice(4))
        SC=>sim_coeff_warm_cg
      CASE default
        ! Should be set to one of the preceding options
        ErrorStatus = -1
        WRITE(cmessage,'(A6,I2)') 'conv_type not recognized: ', conv_type(i)
        
        CALL Ereport(RoutineName, ErrorStatus, Cmessage )
      END SELECT

      ! Working note: a1 and a2 need to be rationalized and the
      ! constants put into tcs_parameters_warm
      arr1(i) = SC%f0(1)*(0.2*scales%wsc_o_mb(i) - 2.)                      &
         + 0.1*scales%wsc_o_mb(i)
      arr2(i) = (SC%f0(1)/SC%f0(2))*(0.2*scales%wsc_o_mb(i) - 2.)

      !-----------------------------------------------------------------------
      !
      !  Functions on rho levels
      !
      !-----------------------------------------------------------------------

      ! Calculate functions
      wexpK(i,:) = 1. - EXP(-SC%K(2)*eta_rho(i,:))
      wexpF(i,:) = 1. - EXP(-SC%Fng(3)*eta_rho(i,:))
      wexpB(i,:) = 1. - EXP(-SC%B(3)*eta_rho(i,:))
      wexpG(i,:) = EXP(-SC%g(5)*eta_rho(i,:))

      DO k=1,sim%nlev
        IF(eta_rho(i,k) < 1.0 - SPACING(eta_rho(i,k)))THEN
          sim%k_func_rho(i,k)   = SC%K(1)*(1. + eta_rho(i,k))*wexpK(i,k)
          sim%fng_func_rho(i,k) = 1. - SC%Fng(1)*(1.+SC%Fng(2)*eta_rho(i,k)) &
             *wexpF(i,k) 
          sim%b_func_rho(i,k)   = SC%b(1)*(1. + eta_rho(i,k)                 &
             - SC%b(2)*eta_rho(i,k)*eta_rho(i,k))                            &
             *wexpB(i,k)
          sim%g_func_rho(i,k)   = SC%g(1)*(eta_rho(i,k)                      &
             - SC%g(2)*eta_rho(i,k)*eta_rho(i,k)                             &
             + (SC%g(3) - SC%g(4)*eta_rho(i,k))                              &
             * wexpG(i,k) - SC%g(3))

          sim%cmask(i,k)=1.
          sim%cmask_rho(i,k)=1.

          IF (eta_rho(i,k) < SC%pc2_d(1))THEN
            sim%pc2_detr(i,k) = SC%pc2_d(2)
          ELSE
            sim%pc2_detr(i,k) = SC%pc2_d(2) +                                &
               (eta_rho(i,k)-SC%pc2_d(1))/(1.-SC%pc2_d(1))*SC%pc2_d(3)
          END IF

        END IF
      END DO

      zi(i,:)= (eta_rho(i,:) - SC%fw(2))/SC%fw(3)
      WHERE(eta_rho(i,:) <= 1.0 + SPACING(eta_rho(i,:)))
        sim%fql_func_rho(i,:) = (1.+eta_rho(i,:)                             &
           - SC%fql(1)*eta_rho(i,:)*eta_rho(i,:))                            &
           *( 1.- SC%fql(2)*EXP(-SC%fql(3)*eta_rho(i,:)))
        sim%gql_func_rho(i,:) = 1. - SC%gql(1)*EXP(-SC%gql(2)*eta_rho(i,:))
        sim%fw_func_rho(i,:)  = SC%fw(1)*EXP(-zi(i,:)*zi(i,:)/2.)            &
           + SC%fw(4) + SC%fw(5)*eta_rho(i,:)                                &
           + SC%fw(6)*eta_rho(i,:)*eta_rho(i,:)
      END WHERE

      !-----------------------------------------------------------------------
      !
      !  Functions on theta levels
      !
      !-----------------------------------------------------------------------



      wexpK(i,:) = 1. - EXP(-SC%K(2)*eta_theta(i,:))
      wexpG(i,:) = EXP(-SC%g(5)*eta_theta(i,:))

      zi(i,:)=(eta_theta(i,:) - SC%fw(2))/SC%fw(3)

      WHERE(eta_theta(i,:) <= 1.0 + SPACING(eta_theta(i,:)))
        sim%g_func(i,:)  = SC%g(1)*(eta_theta(i,:)                           &
           - SC%g(2)*eta_theta(i,:)*eta_theta(i,:)                           &
           + (SC%g(3) - SC%g(4)*eta_theta(i,:))*wexpG(i,:)                   &
           - SC%g(3))

        sim%k_func(i,:)  = SC%k(1)*(1. + eta_theta(i,:))*wexpK(i,:)

        sim%f0_func(i,:) = 1. - arr1(i)*eta_theta(i,:)                       &
           + arr2(i)*(1.-EXP(-SC%f0(2)*eta_theta(i,:)))

        sim%f1_func(i,:) = SC%f1(1)*(eta_theta(i,:)**SC%f1(2))

        sim%ftheta_func(i,:) = SC%fth(1)                                     &
           - SC%fth(2)*EXP(-SC%fth(3)*eta_theta(i,:))
        sim%fw_func(i,:) = SC%fw(1)*EXP(-zi(i,:)*zi(i,:)/2.)                 &
           + SC%fw(4) + SC%fw(5)*eta_theta(i,:)                              &
           + SC%fw(6)*eta_theta(i,:)*eta_theta(i,:)
      END WHERE

    END DO
    IF (lhook) CALL dr_hook('TCS_SIMILARITY:CALC_SIMILARITY',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE calc_similarity

END MODULE tcs_similarity
