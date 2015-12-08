! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module defining the tcs warm rain "scales_conv) derived type and 
! subroutines for allocating and deallocating an instance of this class.
!
MODULE tcs_class_scales


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE
  !
  ! Description:
  ! This module defines the tcs warm rain "scales_conv" derived type and 
  ! subroutines for allocating and deallocating an instance of this class.
  !
  !
  ! Method:
  !   Variables of type "scales_conv" store arrays of pointers to
  ! derived scales (e.g. velocity scales, cloud depth, lcl...)
  ! Each array is defined over convective points.
  ! Note that assignment of these variables is dealt with in the 
  ! tcs_calc_scales
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.2 programming standards.

  TYPE, PUBLIC :: scales_conv
     INTEGER, POINTER :: n_xx
     !
     ! These will have shape (n_xx)
     !
     REAL, POINTER :: wstar_up(:)
                        ! cumulus layer convective velocity scale(m/s)
     REAL, POINTER :: wstar_dn(:)
                        ! subcloud layer velocity scale(m/s)
     REAL, POINTER :: mb(:)
                        ! Cloud base mass flux (m/s)
     REAL, POINTER :: mb_new(:)
                        ! revised mb for incloud calculations (m/s)
     REAL, POINTER :: zcld(:)
                        ! Depth of cloud layer (m)
     REAL, POINTER :: zcld_uv(:)
                        ! Depth of cloud layer (m) CMT cal
     REAL, POINTER :: mb_o_wsc(:)
                        ! mb/wstar_up
     REAL, POINTER :: root_mb_o_wsc(:)
                        ! sqrt of above
     REAL, POINTER :: mb_new_o_wsc(:)
                        ! mb_new/wstar_up
     REAL, POINTER :: root_mb_new_o_wsc(:)
                        ! sqrt of above
     REAL, POINTER :: wstar_up3(:)
                        ! wstar_up**3 * root_mb_new_o_wsc
     REAL, POINTER :: wsc_o_mb(:)
                        ! Convective velocity scale /mb
     REAL, POINTER :: zlcl(:)
                        ! height of the lifting condensation level
     REAL, POINTER :: wup2_cb(:)
                        ! Cloud base updraught velocity squared
                        ! = wup_a1 * wstar_dn(i) * wstar_dn(i)
  END TYPE scales_conv

CONTAINS

  SUBROUTINE allocate_scales(var, n_xx)
    !
    ! Allocates memory to a scales_conv instance "var".   
    ! Inputs "n_xx" defines the sizes of the arrays
    !
    
    IMPLICIT NONE
    
    TYPE(scales_conv) :: var
    INTEGER, OPTIONAL, TARGET :: n_xx
    
    CHARACTER(Len=19), PARAMETER ::  RoutineName = 'allocate_scales'
    CHARACTER(Len=50) :: Message
    INTEGER :: ErrorStatus ! Return code:
                           !   0 = Normal exit
                           ! +ve = Fatal Error
                           ! -ve = Warning

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle
    
    
    IF (lhook) CALL dr_hook('TCS_CLASS_SCALES:ALLOCATE_SCALES',zhook_in,zhook_handle)
    IF (PRESENT(n_xx)) var%n_xx => n_xx
    IF (.NOT. ASSOCIATED(var%n_xx)) THEN
      ErrorStatus=1
      WRITE(Message, '(A50)') &
         ' Error allocating scales: Dimensions not defined '
      
      CALL Ereport(RoutineName, ErrorStatus, Message)
    END IF
    
    ALLOCATE(var%wstar_up(n_xx))
    ALLOCATE(var%wstar_dn(n_xx))
    ALLOCATE(var%mb(n_xx))
    ALLOCATE(var%mb_new(n_xx))
    ALLOCATE(var%zcld(n_xx))
    ALLOCATE(var%zcld_uv(n_xx))
    ALLOCATE(var%mb_o_wsc(n_xx))
    ALLOCATE(var%root_mb_o_wsc(n_xx))
    ALLOCATE(var%mb_new_o_wsc(n_xx))
    ALLOCATE(var%root_mb_new_o_wsc(n_xx))
    ALLOCATE(var%wstar_up3(n_xx))
    ALLOCATE(var%wsc_o_mb(n_xx))
    ALLOCATE(var%zlcl(n_xx))
    ALLOCATE(var%wup2_cb(n_xx))
    IF (lhook) CALL dr_hook('TCS_CLASS_SCALES:ALLOCATE_SCALES',zhook_out,zhook_handle)
    RETURN
 
  END SUBROUTINE allocate_scales

  SUBROUTINE deallocate_scales(var)
    !
    ! Deallocates memory assocciated with scales instance "var".   
    !
    
    IMPLICIT NONE
    TYPE(scales_conv) :: var

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle
    
    IF (lhook) CALL dr_hook('TCS_CLASS_SCALES:DEALLOCATE_SCALES',zhook_in,zhook_handle)
    NULLIFY(var%n_xx)

    DEALLOCATE(var%wup2_cb)
    DEALLOCATE(var%zlcl)
    DEALLOCATE(var%wsc_o_mb)
    DEALLOCATE(var%wstar_up3)
    DEALLOCATE(var%root_mb_new_o_wsc)
    DEALLOCATE(var%mb_new_o_wsc)
    DEALLOCATE(var%root_mb_o_wsc)
    DEALLOCATE(var%mb_o_wsc)
    DEALLOCATE(var%zcld_uv)
    DEALLOCATE(var%zcld)
    DEALLOCATE(var%mb_new)
    DEALLOCATE(var%mb)
    DEALLOCATE(var%wstar_dn)
    DEALLOCATE(var%wstar_up)
    IF (lhook) CALL dr_hook('TCS_CLASS_SCALES:DEALLOCATE_SCALES',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE deallocate_scales

END MODULE tcs_class_scales
