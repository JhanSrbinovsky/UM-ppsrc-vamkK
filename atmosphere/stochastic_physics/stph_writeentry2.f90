! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE stph_writeentry2_mod

IMPLICIT NONE

INTERFACE stph_writeentry2
 MODULE PROCEDURE stph_writeentry2_int, stph_writeentry2_real
END INTERFACE stph_writeentry2

CONTAINS

! -------------------------------------------------------------------
! Subroutine in which stpharray is passed as an 1D integer array
! -------------------------------------------------------------------
SUBROUTINE stph_writeentry2_int( stpharray, sizearray)

! This routine is passed a integer array and its size as arguments
! and writes unformatted to the pre-defined stoch phys
! output file unit=149

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim

 IMPLICIT NONE

! Arguments
 INTEGER, INTENT(IN) :: sizearray            ! The size of the array
 INTEGER, INTENT(IN) :: stpharray(sizearray) ! The random seed

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_WRITEENTRY2_INT',zhook_in,zhook_handle)

! Write the value to the file stream
 WRITE(149) stpharray

 IF (lhook) CALL dr_hook('STPH_WRITEENTRY2_INT',zhook_out,zhook_handle)

END SUBROUTINE stph_writeentry2_int

! -------------------------------------------------------------------
! Subroutine in which stpharray is passed as a 2D real array
! -------------------------------------------------------------------
SUBROUTINE stph_writeentry2_real( stpharray, sizearray)

! This routine is passed a real array and its size as arguments
! and writes unformatted to the pre-defined stoch phys
! output file unit=149

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim

 IMPLICIT NONE

! Arguments
 REAL,    INTENT(IN) :: stpharray(:,:) ! The random seed
 INTEGER, INTENT(IN) :: sizearray      ! The size of the array

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_WRITEENTRY2_REAL',zhook_in,zhook_handle)

! Write the value to the file stream
 WRITE(149) stpharray

 IF (lhook) CALL dr_hook('STPH_WRITEENTRY2_REAL',zhook_out,zhook_handle)

END SUBROUTINE stph_writeentry2_real

END MODULE stph_writeentry2_mod
