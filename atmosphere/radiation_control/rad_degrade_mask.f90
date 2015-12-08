! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Create a chequerboard array of logicals over the whole domain

!  Purpose: To create a logical array which has a true chequer-board
!           pattern over the whole model DOMAIN (not just the PE).

!  Method:  The datastart variable stores the co-ordinate of the top-
!           left grid box of a PE in terms of the DOMAIN co-ordinates.
!           It can be used in conjunction with the co-ordinates of the
!           grid box in terms of the PE to ensure that the whole domain,
!           and not just a single PE, has a chequer-board pattern of
!           radiation calculations.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

      SUBROUTINE rad_degrade_mask(rad_mask, datastart, ndim_max,        &
                        first_row, last_row, offx, row_length, rows)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: row_length
!            The data length of a row on the PE (control level)
      INTEGER, INTENT(IN) :: rows
!            The number of rows on the PE (control level)
      INTEGER, INTENT(IN) :: first_row
!            first non-polar row of field (misses halo for MPP code)
      INTEGER, INTENT(IN) :: last_row
!            last non-polar row of field (misses halo for MPP code)
      INTEGER, INTENT(IN) :: offx
!            halo size in EW direction
      INTEGER, INTENT(IN) :: ndim_max
!            maximum number of spatial dimensions
      INTEGER, INTENT(IN) :: datastart(ndim_max)
!            position of personal data in global data
!            (set at very high level in control code)

      LOGICAL, INTENT(OUT) :: rad_mask(row_length, rows)
!            A chequer-board array of values such that the whole DOMAIN
!            is a consistent chequerboard


!     Local variables

      INTEGER :: i              ! loop variable
      INTEGER :: j              ! loop variable
      INTEGER :: tempint        ! temporary storage variable

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('RAD_DEGRADE_MASK',zhook_in,zhook_handle)

! Note that all calculations must be done in the data part of the array,
! and that the co-ordinates of such data points are relative to the
! first data point of the whole domain.
! The formula below is a fast method of producing a chequer-board
! pattern over the whole domain.

      DO j=1,last_row
        DO i=1,row_length
          tempint=MOD((datastart(1)+i+datastart(2)+j),2)
          rad_mask(i,j)=(tempint == 0)
        END DO
      END DO

      IF (lhook) CALL dr_hook('RAD_DEGRADE_MASK',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE rad_degrade_mask
