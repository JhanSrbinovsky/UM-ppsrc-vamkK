! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module defining the tcs warm rain "similarity" derived type and 
! subroutines for allocating and deallocating an instance of this class.
!
MODULE tcs_class_similarity


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE
  !
  ! Description:
  ! This module defines the tcs warm rain "similarity" derived type and 
  ! subroutines for allocating and deallocating an instance of this class.
  !
  !
  ! Method:
  !   Variables of type "similarity" store arrays of in-cloud similarity
  ! functions as described in 
  !   <reference to documentation to go here, once available>
  ! Each array is defined over convective points and in-cloud levels.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.2 programming standards.


  TYPE, PUBLIC :: similarity
     INTEGER, POINTER :: n_xx
     INTEGER, POINTER :: nlev
     !
     ! These will have shape (n_xx,nlev)
     !
     !--------------------------------------------------------------
     ! These are similarity functions
     !--------------------------------------------------------------
     REAL, POINTER :: g_func(:,:)
                                ! g function
     REAL, POINTER :: k_func(:,:)
                                ! k function
     REAL, POINTER :: f0_func(:,:)
                                ! f0 function
     REAL, POINTER :: f1_func(:,:)
                                ! f1 function
     REAL, POINTER :: fw_func(:,:)
                                ! fw function
     REAL, POINTER :: ftheta_func(:,:)
                                ! ftheta function
     REAL, POINTER :: cmask(:,:)
                                ! cloud mask
     REAL, POINTER :: pc2_detr(:,:)
                                ! detrainment rate for pc2
     !
     ! Arrays for functions of height - uv levels
     !
     REAL, POINTER :: fql_func_rho(:,:)
                                ! fql function
     REAL, POINTER :: g_func_rho(:,:)
                                ! g function
     REAL, POINTER :: fw_func_rho(:,:)
                                ! fw function
     REAL, POINTER :: fng_func_rho(:,:)
                                ! fng function
     REAL, POINTER :: k_func_rho(:,:)
                                ! k function
     REAL, POINTER :: b_func_rho(:,:)
                                ! b function
     REAL, POINTER :: gql_func_rho(:,:)
                                ! gql function
     REAL, POINTER :: cmask_rho(:,:)       
                                ! cloud mask
  END TYPE similarity

CONTAINS

    SUBROUTINE allocate_similarity(var, n_xx, nlev, initval)
    !
    ! Allocates memory to a similarity instance "var".   
    ! Inputs "n_xx" and "nlev" define the sizes of the arrays
    ! and if present all fields are initialized to the value "initval"
    !
    IMPLICIT NONE
        
    TYPE(similarity), INTENT(inout) :: var
    INTEGER, INTENT(in), OPTIONAL, TARGET :: n_xx, nlev
    REAL, INTENT(in), OPTIONAL :: initval
    LOGICAL :: l_init
    
    CHARACTER(Len=19), PARAMETER ::  RoutineName = 'allocate_similarity'
    CHARACTER(Len=53) :: Message
    INTEGER :: ErrorStatus ! Return code:
    !   0 = Normal exit
    ! +ve = Fatal Error
    ! -ve = Warning

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle
       
    IF (lhook) CALL dr_hook('TCS_CLASS_SIMILARITY:ALLOCATE_SIMILARITY',zhook_in,zhook_handle)
    IF (PRESENT(n_xx)) var%n_xx => n_xx
    IF (PRESENT(nlev)) var%nlev => nlev
    IF (.NOT. (ASSOCIATED(var%n_xx)          &
               .AND. ASSOCIATED(var%nlev))) THEN
      ErrorStatus=1
      WRITE(Message, '(A53)') &
         ' Error allocating similarity: Dimensions not defined '
      
      CALL Ereport(RoutineName, ErrorStatus, Message)
    END IF
    
    ALLOCATE(var%g_func(var%n_xx,var%nlev))
    ALLOCATE(var%k_func(var%n_xx,var%nlev))
    ALLOCATE(var%f0_func(var%n_xx,var%nlev))
    ALLOCATE(var%f1_func(var%n_xx,var%nlev))
    ALLOCATE(var%fw_func(var%n_xx,var%nlev))
    ALLOCATE(var%ftheta_func(var%n_xx,var%nlev))
    ALLOCATE(var%cmask(var%n_xx,var%nlev))
    ALLOCATE(var%pc2_detr(var%n_xx,var%nlev))
    ALLOCATE(var%fql_func_rho(var%n_xx,var%nlev))
    ALLOCATE(var%g_func_rho(var%n_xx,var%nlev))
    ALLOCATE(var%fw_func_rho(var%n_xx,var%nlev))
    ALLOCATE(var%fng_func_rho(var%n_xx,var%nlev))
    ALLOCATE(var%k_func_rho(var%n_xx,var%nlev))
    ALLOCATE(var%b_func_rho(var%n_xx,var%nlev))
    ALLOCATE(var%gql_func_rho(var%n_xx,var%nlev))
    ALLOCATE(var%cmask_rho(var%n_xx,var%nlev))
    
    ! initialise if required.
    IF (PRESENT(initval))THEN
       l_init=.TRUE. 
    ELSE 
       l_init=.FALSE.
    END IF
    IF (l_init)THEN
       var%g_func       = REAL(initval)
       var%k_func       = REAL(initval)
       var%f0_func      = REAL(initval)
       var%f1_func      = REAL(initval)
       var%fw_func      = REAL(initval)
       var%ftheta_func  = REAL(initval)
       var%cmask        = REAL(initval)
       var%pc2_detr     = REAL(initval)
       var%fql_func_rho = REAL(initval)
       var%g_func_rho   = REAL(initval)
       var%fw_func_rho  = REAL(initval)
       var%fng_func_rho = REAL(initval)
       var%k_func_rho   = REAL(initval)
       var%b_func_rho   = REAL(initval)
       var%gql_func_rho = REAL(initval)
       var%cmask_rho    = REAL(initval)
    END IF
    IF (lhook) CALL dr_hook('TCS_CLASS_SIMILARITY:ALLOCATE_SIMILARITY',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE allocate_similarity

  SUBROUTINE deallocate_similarity(var)
    !
    ! Deallocates memory assocciated with similarity instance "var".   
    !
    IMPLICIT NONE
    TYPE(similarity) :: var

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle
    
    IF (lhook) CALL dr_hook('TCS_CLASS_SIMILARITY:DEALLOCATE_SIMILARITY',zhook_in,zhook_handle)
    NULLIFY(var%n_xx)
    NULLIFY(var%nlev)

    DEALLOCATE(var%cmask_rho)
    DEALLOCATE(var%gql_func_rho)
    DEALLOCATE(var%b_func_rho)
    DEALLOCATE(var%k_func_rho)
    DEALLOCATE(var%fng_func_rho)
    DEALLOCATE(var%fw_func_rho)
    DEALLOCATE(var%g_func_rho)
    DEALLOCATE(var%fql_func_rho)
    DEALLOCATE(var%pc2_detr)
    DEALLOCATE(var%cmask)
    DEALLOCATE(var%ftheta_func)
    DEALLOCATE(var%fw_func)
    DEALLOCATE(var%f1_func)
    DEALLOCATE(var%f0_func)
    DEALLOCATE(var%k_func)
    DEALLOCATE(var%g_func)
    IF (lhook) CALL dr_hook('TCS_CLASS_SIMILARITY:DEALLOCATE_SIMILARITY',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE deallocate_similarity

END MODULE tcs_class_similarity
