! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE FLUX_DIAG----------------------------------------------
!LL
!LL  PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!LL            TO SUM FLUXES INTO THE ATMOSPHERE
!LL            Note global summation now done elsewhere at end of
!LL            energy correction period to save CPU.
!LL
!LL  NOT SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION :
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Energy Correction
      SUBROUTINE FLUX_DIAG (FLUX,AREA,row_length, rows,off_x,off_y,     &
     &                      CONV_FAC,FLUX_SUM, TSTEP)
!
      USE earth_constants_mod, ONLY: earth_radius

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

!----------------------------------------------------------------------
! VECTOR LENGTHS
!----------------------------------------------------------------------
!
      INTEGER row_length, rows                                          &
     &,     off_x,off_y

!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------

      REAL FLUX(row_length, rows) ! IN FLUX TO BE SUMMED
!
      REAL AREA(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y) ! IN AREA OF GRID BOX
!
      REAL TSTEP               ! IN TIMESTEP
!
      REAL CONV_FAC            ! IN CONVERSION FACTOR TO TRANSLATE
                               !    FLUX INTO ENERGY UNITS
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE IN AND OUT
!----------------------------------------------------------------------
!

      REAL FLUX_SUM(row_length,rows)  ! INOUT sum of fluxes into
                                      ! the atmosphere.
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
!
      REAL FACTOR        ! local factor

      INTEGER I,j ! LOOP COUNTER

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!
!----------------------------------------------------------------------
! EXTERNAL SUBROUTINE CALLS  -  NONE
!L---------------------------------------------------------------------
!L Flux_sum holds the summation of fluxes into the atmosphere as a full
!l global field weighted by area. The fluxes are summed over the period
!L A_energysteps. At the end of this period the global sum is formed.
!L Global summation is costly on a MPP machine and is therefore only
!L done at the end of the period instead of every time fluxes are added.
!L Note this appraoch means the array flux_sum must be held in memory
!L adding to the memory overhead.
!L
!L Add additional fluxes to current flux_sum.
!L---------------------------------------------------------------------
!L factor to multiply flux

      IF (lhook) CALL dr_hook('FLUX_DIAG',zhook_in,zhook_handle)
      factor = earth_radius*earth_radius*conv_fac*tstep

      DO j=1,rows
        DO I=1,row_length
          flux_sum(I,j) = flux_sum(i,j)+                                &
     &                           factor*AREA(I,j)*FLUX(I,j)
        END DO
      END DO
!
      IF (lhook) CALL dr_hook('FLUX_DIAG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FLUX_DIAG
