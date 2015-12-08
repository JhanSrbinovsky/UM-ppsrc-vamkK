! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Temporary routine to pad arrays. This is to prevent uninitialised 
! data getting into the SCM output and should only be used until all
! SCM diagnostics are properly initialised
!
MODULE tcs_pad


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  !
  ! Description:
  !   Fills in uninitialised data
  !
  ! Method:
  !   Extends the input array to the size of the output array
  !   and fills extra array entries with some "padding" value.
  !   Subroutines for 1d and 2D arrays are overloaded to a single
  !   interface "pad".
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.1 programming standards.
  !

  INTERFACE pad
     MODULE PROCEDURE pad1
     MODULE PROCEDURE pad2
  END INTERFACE

CONTAINS

  SUBROUTINE pad1(in, padded, padval)
    IMPLICIT NONE
    !------------------------------------------------------------------
    !
    ! Pads array in so that its shape the same as padded 
    ! Rank 1 real version
    !
    !------------------------------------------------------------------
    REAL, INTENT(in)    :: in(:)
    REAL, INTENT(inout) :: padded(:)
    REAL, INTENT(in), OPTIONAL :: padval

    REAL :: val

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('TCS_PAD:PAD1',zhook_in,zhook_handle)
    IF (PRESENT(padval))THEN
      val=padval
    ELSE
      val=0.0
    END IF

    padded=val
    padded(1:SIZE(in,1)) = in(:)
    IF (lhook) CALL dr_hook('TCS_PAD:PAD1',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE pad1

  SUBROUTINE pad2(in, padded, padval)
    IMPLICIT NONE
    !------------------------------------------------------------------
    !
    ! Pads array in so that its shape the same as padded 
    ! Rank 2 real version
    !
    !------------------------------------------------------------------
    REAL, INTENT(in)    :: in(:,:)
    REAL, INTENT(inout) :: padded(:,:)
    REAL, INTENT(in), OPTIONAL :: padval

    REAL :: val

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('TCS_PAD:PAD2',zhook_in,zhook_handle)
    IF (PRESENT(padval))THEN
       val=padval
    ELSE
       val=0.0
    END IF

    padded=val
    padded(1:SIZE(in,1),1:SIZE(in,2)) = in(:,:)
    IF (lhook) CALL dr_hook('TCS_PAD:PAD2',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE pad2

END MODULE tcs_pad
