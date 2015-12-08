! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE stph_readentry2_mod

IMPLICIT NONE

INTERFACE  stph_readentry2
 MODULE PROCEDURE stph_readentry2_int, stph_readentry2_real
END INTERFACE

CONTAINS

! -------------------------------------------------------------------
! Subroutine in which stpharray is passed as a 1D integer array
! -------------------------------------------------------------------
SUBROUTINE stph_readentry2_int( stpharray, sizearray)

! This routine reads in the random seed and stochastic wave coefficients
! from unformatted file written previously at dump steps.

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim

 IMPLICIT NONE

! Arguments
 INTEGER, INTENT(IN)  :: sizearray            ! The size of the array
 INTEGER, INTENT(OUT) :: stpharray(sizearray) ! The random seed

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_READENTRY2_INT',zhook_in,zhook_handle)

! Read array from seed file (unformatted)
 READ(149) stpharray

 IF (lhook) CALL dr_hook('STPH_READENTRY2_INT',zhook_out,zhook_handle)
 RETURN
END SUBROUTINE stph_readentry2_int

! -------------------------------------------------------------------
! Subroutine in which stpharray is passed as a 2D real array
! -------------------------------------------------------------------
SUBROUTINE stph_readentry2_real( stpharray, sizearray)

! This routine reads in the random seed and stochastic wave coefficients
! from unformatted file written previously at dump steps.

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim

 IMPLICIT NONE

! Arguments
 REAL,    INTENT(OUT) :: stpharray(:,:) ! The random seed
 INTEGER, INTENT(IN)  :: sizearray      ! The size of the array

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_READENTRY2_REAL',zhook_in,zhook_handle)

! Read array from seed file (unformatted)
 READ(149) stpharray

 IF (lhook) CALL dr_hook('STPH_READENTRY2_REAL',zhook_out,zhook_handle)
 RETURN
END SUBROUTINE stph_readentry2_real

END MODULE stph_readentry2_mod
