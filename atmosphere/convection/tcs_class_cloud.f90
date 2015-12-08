! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module defining the tcs warm rain cloud_input derived type and 
! subroutines for allocating and deallocating an instance of this class
! and also for assigning values.
!
MODULE tcs_class_cloud


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
  !   Variables of type "cloud_input" store arrays of pointers to
  ! those sections of the existing full model fields, or derivations of 
  ! them, within the cloud layer. Each array is defined over convective
  ! points and in-cloud levels. 
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.2 programming standards.

  TYPE :: cloud_input
     !--------------------------------------------------------------
     ! These are arrays which are sampled on cloud levels
     !--------------------------------------------------------------
     INTEGER, POINTER :: n_xx  
                               ! Number of convectively active points
     INTEGER, POINTER :: nlev  
                               ! Number of levels
     INTEGER, POINTER :: ntra  
                               ! Number of tracers
     ! These will have shape (n_xx,nlev)
     REAL, POINTER :: eta_theta(:,:)
                               ! Non-dimensional height (theta levels)
     REAL, POINTER :: eta_rho(:,:)
                               ! Non-dimensional height (rho levels)
     REAL, POINTER :: theta(:,:)
                               ! theta
     REAL, POINTER :: q_mix(:,:)
                               ! q - mixing ratio
     REAL, POINTER :: qse(:,:)
                               ! saturation mixing ratio
     REAL, POINTER :: exner_theta(:,:)
                               ! exner on theta levels in cloud
     REAL, POINTER :: exner_rho(:,:)
                               ! exner pressure on rho in cloud
     REAL, POINTER :: p_theta(:,:)
                               ! pressure on theta levels in cloud
     REAL, POINTER :: p_rho(:,:)
                               ! pressure on rho levels in cloud
     REAL, POINTER :: z_theta(:,:)
                               ! height of theta levels
     REAL, POINTER :: z_rho(:,:)
                               ! height of rho levels
     REAL, POINTER :: r2rho(:,:)
                               ! r2*rho rho levels (kg/m)
     REAL, POINTER :: r2rho_theta(:,:)
                               ! r2*rho theta levels (kg/m)
     REAL, POINTER :: rho(:,:)
                               ! density on  rho levels
     ! These will have shape (n_xx,nlev,ntra)
     REAL, POINTER :: tracer(:,:,:)  
                               ! tracer in cloud
  END TYPE cloud_input

CONTAINS

  SUBROUTINE allocate_cloud_input(var, n_xx, nlev, ntra, initval)
    !
    ! Allocates memory to a cloud_input instance "var".   
    ! Inputs "n_xx", "nlev" and "ntra" define the sizes of the arrays
    ! and if present all fields are initialized to the value "initval"
    !
    IMPLICIT NONE
    
    TYPE(cloud_input) :: var
    INTEGER, OPTIONAL, TARGET :: n_xx, nlev, ntra
    REAL, INTENT(in), OPTIONAL :: initval
    LOGICAL :: l_init
    
    CHARACTER(Len=20), PARAMETER ::  RoutineName = 'allocate_cloud_input'
    CHARACTER(Len=54) :: Message
    INTEGER :: ErrorStatus ! Return code:
    !   0 = Normal exit
    ! +ve = Fatal Error
    ! -ve = Warning

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('TCS_CLASS_CLOUD:ALLOCATE_CLOUD_INPUT',zhook_in,zhook_handle)

    IF (PRESENT(n_xx)) var%n_xx => n_xx
    IF (PRESENT(nlev)) var%nlev => nlev
    IF (PRESENT(ntra)) var%ntra => ntra
    IF (.NOT. (ASSOCIATED(var%n_xx)          &
               .AND. ASSOCIATED(var%nlev)    &
               .AND. ASSOCIATED(var%ntra))) THEN
      ErrorStatus=1
      WRITE(Message, '(A54)') &
         ' Error allocating cloud_input: Dimensions not defined '
      
      CALL Ereport(RoutineName, ErrorStatus, Message)
    END IF
    
    ALLOCATE(var%eta_theta(var%n_xx,var%nlev))
    ALLOCATE(var%eta_rho(var%n_xx,var%nlev))
    ALLOCATE(var%theta(var%n_xx,var%nlev))
    ALLOCATE(var%q_mix(var%n_xx,var%nlev))
    ALLOCATE(var%qse(var%n_xx,var%nlev))
    ALLOCATE(var%exner_theta(var%n_xx,var%nlev))
    ALLOCATE(var%exner_rho(var%n_xx,var%nlev))
    ALLOCATE(var%p_theta(var%n_xx,var%nlev))
    ALLOCATE(var%p_rho(var%n_xx,var%nlev))
    ALLOCATE(var%z_theta(var%n_xx,var%nlev))
    ALLOCATE(var%z_rho(var%n_xx,var%nlev))
    ALLOCATE(var%r2rho(var%n_xx,var%nlev))
    ALLOCATE(var%r2rho_theta(var%n_xx,var%nlev))
    ALLOCATE(var%rho(var%n_xx,var%nlev))
    ALLOCATE(var%tracer(var%n_xx,var%nlev,var%ntra))
    
    ! initialise if required.
    IF (PRESENT(initval))THEN
       l_init=.TRUE. 
    ELSE 
       l_init=.FALSE.
    END IF
    IF (l_init)THEN
       var%eta_theta   =REAL(initval)
       var%eta_rho     =REAL(initval)
       var%theta       =REAL(initval)
       var%q_mix       =REAL(initval)
       var%qse         =REAL(initval)
       var%exner_theta =REAL(initval)
       var%exner_rho   =REAL(initval)
       var%p_theta     =REAL(initval)
       var%p_rho       =REAL(initval)
       var%z_theta     =REAL(initval)
       var%z_rho       =REAL(initval)
       var%r2rho       =REAL(initval)
       var%r2rho_theta =REAL(initval)
       var%rho         =REAL(initval)
       var%tracer      =REAL(initval)
    END IF
    IF (lhook) CALL dr_hook('TCS_CLASS_CLOUD:ALLOCATE_CLOUD_INPUT',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE allocate_cloud_input

  SUBROUTINE deallocate_cloud_input(var)
    !
    ! Deallocates memory assocciated with cloud_input instance "var".   
    !
    
    IMPLICIT NONE
    TYPE(cloud_input) :: var

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle
    
    IF (lhook) CALL dr_hook('TCS_CLASS_CLOUD:DEALLOCATE_CLOUD_INPUT',zhook_in,zhook_handle)
    NULLIFY(var%n_xx)
    NULLIFY(var%nlev)
    NULLIFY(var%ntra)

    DEALLOCATE(var%tracer)
    DEALLOCATE(var%rho)
    DEALLOCATE(var%r2rho_theta)
    DEALLOCATE(var%r2rho)
    DEALLOCATE(var%z_rho)
    DEALLOCATE(var%z_theta)
    DEALLOCATE(var%p_rho)
    DEALLOCATE(var%p_theta)
    DEALLOCATE(var%exner_rho)
    DEALLOCATE(var%exner_theta)
    DEALLOCATE(var%qse)
    DEALLOCATE(var%q_mix)
    DEALLOCATE(var%theta)
    DEALLOCATE(var%eta_theta)
    DEALLOCATE(var%eta_rho)
    IF (lhook) CALL dr_hook('TCS_CLASS_CLOUD:DEALLOCATE_CLOUD_INPUT',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE deallocate_cloud_input

  SUBROUTINE assign_cloud_input(                                  &
       var, base_levels,  maxlevs , maxtrlevs,                    &
       theta, q_mix, qse, p_theta, p_rho, exner_theta, exner_rho, &
       z_theta, z_rho, r2rho, r2rho_theta, rho, tracer ) 
    !
    ! Assigns values to a cloud_input instance "var".   
    !
  
    USE tcs_common_warm , ONLY : scales, cb
    IMPLICIT NONE
    TYPE(cloud_input), INTENT(inout) :: var
    INTEGER, INTENT(in) :: base_levels(:)  ! Levels at base of region
    ! Working note: maxlevs would normally be 
    ! maxval(levels%ntpar(:) - levels%ntml(:) + 1), but for now
    ! we allow this to be input since it differs in certain 
    ! places in the old code.
    INTEGER, INTENT(in) :: maxlevs         ! Maximum number of levels
                                           ! within region
    INTEGER, INTENT(in) :: maxtrlevs       ! Maximum number of levels
                                           ! for tracers within region
    REAL, INTENT(in) :: theta(:,:)
    REAL, INTENT(in) :: q_mix(:,:)
    REAL, INTENT(in) :: qse(:,:)
    REAL, INTENT(in) :: p_theta(:,0:)      ! input starts 0
    REAL, INTENT(in) :: p_rho(:,0:)        ! input starts 0
    REAL, INTENT(in) :: exner_theta(:,0:)  ! input starts 0
    REAL, INTENT(in) :: exner_rho(:,0:)    ! input starts 0
    REAL, INTENT(in) :: z_theta(:,:)
    REAL, INTENT(in) :: z_rho(:,:)
    REAL, INTENT(in) :: r2rho(:,:)
    REAL, INTENT(in) :: r2rho_theta(:,:)
    REAL, INTENT(in) :: rho(:,:)
    REAL, INTENT(in) :: tracer(:,:,:)
    
    INTEGER :: i, lbase, ltop, maxsize

    CHARACTER(Len=18), PARAMETER ::  RoutineName = 'assign_cloud_input'
    CHARACTER(Len=90) :: Message
    INTEGER :: ErrorStatus ! Return code:
    !   0 = Normal exit
    ! +ve = Fatal Error
    ! -ve = Warning

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('TCS_CLASS_CLOUD:ASSIGN_CLOUD_INPUT',zhook_in,zhook_handle)

    maxsize=SIZE(theta(1,:))
    DO i=1,var%n_xx
       lbase = base_levels(i)
       ltop  = lbase + maxlevs - 1
       IF (ltop >= maxsize)THEN
         ErrorStatus=1
         WRITE(Message, '(A75, 5I3)') &
            'Too many levels in warm convective layer: &
            &i, lbase, ltop, maxlevs, maxsize ='       &
            , i, lbase, ltop, maxlevs, maxsize
         
         CALL Ereport(RoutineName, ErrorStatus, Message)
       END IF
       var%eta_theta(i,1:maxlevs) = (z_theta(i,lbase:ltop) - cb%z_rho(i))   &
                            /scales%zcld(i)
       var%eta_rho(i,1:maxlevs) = (z_rho(i,lbase:ltop) - cb%z_rho(i))       &
                          /scales%zcld(i)
       var%theta(i,1:maxlevs)       = theta(i,lbase:ltop)
       var%q_mix(i,1:maxlevs)       = q_mix(i,lbase:ltop)
       var%qse(i,1:maxlevs)         = qse(i,lbase:ltop)
       var%p_theta(i,1:maxlevs)     = p_theta(i,lbase:ltop)
       var%p_rho(i,1:maxlevs)       = p_rho(i,lbase:ltop)
       var%exner_theta(i,1:maxlevs) = exner_theta(i,lbase:ltop)
       var%exner_rho(i,1:maxlevs)   = exner_rho(i,lbase:ltop)
       var%z_theta(i,1:maxlevs)     = z_theta(i,lbase:ltop)
       var%z_rho(i,1:maxlevs)       = z_rho(i,lbase:ltop)
       var%r2rho(i,1:maxlevs)       = r2rho(i,lbase:ltop)
       var%r2rho_theta(i,1:maxlevs) = r2rho_theta(i,lbase:ltop)
       var%rho(i,1:maxlevs)         = rho(i,lbase:ltop)
       ltop  = lbase + maxtrlevs - 1
       var%tracer(i,1:maxtrlevs,:)    = tracer(i,lbase:ltop,:)
    END DO
    IF (lhook) CALL dr_hook('TCS_CLASS_CLOUD:ASSIGN_CLOUD_INPUT',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE assign_cloud_input

END MODULE tcs_class_cloud
