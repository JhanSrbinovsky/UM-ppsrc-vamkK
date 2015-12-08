! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module defining the tcs warm rain "interface_input" derived type and 
! subroutines for allocating and deallocating an instance of this class
! and also for assigning values.
!
MODULE tcs_class_interface


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE
  !
  ! Description:
  ! This module defines the tcs warm rain "cloud_input" derived type and 
  ! subroutines for allocating and deallocating an instance of this class.
  !
  !
  ! Method:
  !   Variables of type "interface_input" store arrays of pointers to
  ! arrays sampled on interface levels (e.g.  cloud base, freezing level,
  ! inversion base...)
  ! Each array is defined over convective points.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.2 programming standards.

  TYPE, PUBLIC :: interface_input
     !
     ! These are fields evaluated on interface levels,
     ! e.g. cloudbase (ntml), inversion (ntpar) and levels around these
     !
     INTEGER, POINTER :: n_xx
     ! Number of convectively active points
     INTEGER, POINTER :: ntra           
     ! Number of tracers
     ! These fields will all have shape (n_xx)
     INTEGER, POINTER :: levels(:) ! The levels on which the field is evaluated
     REAL, POINTER :: theta(:)
     REAL, POINTER :: q_mix(:)
     REAL, POINTER :: qse(:)
     REAL, POINTER :: p_theta(:)  ! = p_layer_centres
     REAL, POINTER :: p_rho(:)    ! = p_layer_boundaries
     REAL, POINTER :: exner_theta(:)   ! = exner_layer_centres
     REAL, POINTER :: exner_rho(:)     ! = exner_layer_boundaries
     REAL, POINTER :: z_theta(:)
     REAL, POINTER :: z_rho(:)
     REAL, POINTER :: rho(:)
     !These will have shape (n_xx,ntra)
     REAL, POINTER :: tracer(:,:)
  END TYPE interface_input

CONTAINS

  SUBROUTINE allocate_interface_input(var, n_xx, ntra, initval)
    !
    ! Allocates memory to a interface_input instance "var".   
    ! Inputs "n_xx" and "ntra" define the sizes of the arrays
    ! and if present all fields are initialized to the value "initval"
    !
    IMPLICIT NONE
    TYPE(interface_input) :: var
    INTEGER, OPTIONAL, TARGET :: n_xx, ntra
    REAL, INTENT(in), OPTIONAL :: initval
    LOGICAL :: l_init
    
    CHARACTER(Len=24), PARAMETER ::  RoutineName = 'allocate_interface_input'
    CHARACTER(Len=58) :: Message
    INTEGER :: ErrorStatus ! Return code:
    !   0 = Normal exit
    ! +ve = Fatal Error
    ! -ve = Warning

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('TCS_CLASS_INTERFACE:ALLOCATE_INTERFACE_INPUT',zhook_in,zhook_handle)
    IF (PRESENT(n_xx)) var%n_xx => n_xx
    IF (PRESENT(ntra)) var%ntra => ntra
    IF (.NOT. (ASSOCIATED(var%n_xx) .AND. ASSOCIATED(var%ntra))) THEN
      ErrorStatus=1
      WRITE(Message, '(A54)') &
         ' Error allocating interface_input: Dimensions not defined '
      
      CALL Ereport(RoutineName, ErrorStatus, Message)
    END IF

    ! variables with shape (n_xx)
    ALLOCATE(var%levels(var%n_xx))
    ALLOCATE(var%theta(var%n_xx))
    ALLOCATE(var%q_mix(var%n_xx))
    ALLOCATE(var%qse(var%n_xx))
    ALLOCATE(var%p_theta(var%n_xx))
    ALLOCATE(var%p_rho(var%n_xx))
    ALLOCATE(var%exner_theta(var%n_xx))
    ALLOCATE(var%exner_rho(var%n_xx))
    ALLOCATE(var%z_theta(var%n_xx))
    ALLOCATE(var%z_rho(var%n_xx))
    ALLOCATE(var%rho(var%n_xx))

    ! variables with shape (n_xx, ntra)
    ALLOCATE(var%tracer(var%n_xx,var%ntra))

    ! initialise if required.
    IF (PRESENT(initval))THEN
      l_init=.TRUE. 
    ELSE 
      l_init=.FALSE.
    END IF
    IF (l_init)THEN
      var%theta       =REAL(initval)
      var%q_mix       =REAL(initval)
      var%qse         =REAL(initval)
      var%p_theta     =REAL(initval)
      var%p_rho       =REAL(initval)
      var%exner_theta =REAL(initval)
      var%exner_rho   =REAL(initval)
      var%z_theta     =REAL(initval)
      var%z_rho       =REAL(initval)
      var%rho         =REAL(initval)
      var%tracer      =REAL(initval)
    END IF
    IF (lhook) CALL dr_hook('TCS_CLASS_INTERFACE:ALLOCATE_INTERFACE_INPUT',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE allocate_interface_input

  SUBROUTINE deallocate_interface_input(var)
    !
    ! Deallocates memory assocciated with cloud_input instance "var".   
    !
    
    IMPLICIT NONE
    TYPE(interface_input) :: var

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('TCS_CLASS_INTERFACE:DEALLOCATE_INTERFACE_INPUT',zhook_in,zhook_handle)
    NULLIFY(var%n_xx)
    NULLIFY(var%ntra)
    !    nullify(var%levels)

    DEALLOCATE(var%tracer)

    DEALLOCATE(var%rho)
    DEALLOCATE(var%z_rho)
    DEALLOCATE(var%z_theta)
    DEALLOCATE(var%exner_rho)
    DEALLOCATE(var%exner_theta)
    DEALLOCATE(var%p_rho)
    DEALLOCATE(var%p_theta)
    DEALLOCATE(var%qse)
    DEALLOCATE(var%q_mix)
    DEALLOCATE(var%theta)
    DEALLOCATE(var%levels)
    IF (lhook) CALL dr_hook('TCS_CLASS_INTERFACE:DEALLOCATE_INTERFACE_INPUT',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE deallocate_interface_input

  SUBROUTINE assign_interface_input(var, levels, theta,    &
     q_mix, qse, p_theta, p_rho, exner_theta, exner_rho, &
     z_theta, z_rho, rho, tracer ) 
    !
    ! Assigns values to a interface_input instance "var".   
    !
    IMPLICIT NONE
    TYPE(interface_input), INTENT(inout) :: var
    INTEGER, INTENT(in), TARGET :: levels(:)
    REAL, INTENT(in) :: theta(:,:)
    REAL, INTENT(in) :: q_mix(:,:)
    REAL, INTENT(in) :: qse(:,:)
    REAL, INTENT(in) :: p_theta(:,0:)      ! input starts 0
    REAL, INTENT(in) :: p_rho(:,0:)        ! input starts 0
    REAL, INTENT(in) :: exner_theta(:,0:)  ! input starts 0
    REAL, INTENT(in) :: exner_rho(:,0:)    ! input starts 0
    REAL, INTENT(in) :: z_theta(:,:)
    REAL, INTENT(in) :: z_rho(:,:)
    REAL, INTENT(in) :: rho(:,:)
    REAL, INTENT(in) :: tracer(:,:,:)

    INTEGER :: i

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle


    IF (lhook) CALL dr_hook('TCS_CLASS_INTERFACE:ASSIGN_INTERFACE_INPUT',zhook_in,zhook_handle)
    DO i=1,var%n_xx
      var%levels(i)      = levels(i)
      var%theta(i)       = theta(i,levels(i))
      var%q_mix(i)       = q_mix(i,levels(i))
      var%qse(i)         = qse(i,levels(i))
      var%p_theta(i)     = p_theta(i,levels(i))
      var%p_rho(i)       = p_rho(i,levels(i))
      var%exner_theta(i) = exner_theta(i,levels(i))
      var%exner_rho(i)   = exner_rho(i,levels(i))
      var%z_theta(i)     = z_theta(i,levels(i))
      var%z_rho(i)       = z_rho(i,levels(i))
      var%rho(i)         = rho(i,levels(i))
      ! Working note: need to put on a check that tracer
      ! levels go up to levels(i) - c.f. also assign_cloud_input
      var%tracer(i,:)    = tracer(i,levels(i),:)
    END DO
    IF (lhook) CALL dr_hook('TCS_CLASS_INTERFACE:ASSIGN_INTERFACE_INPUT',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE assign_interface_input

END MODULE tcs_class_interface
