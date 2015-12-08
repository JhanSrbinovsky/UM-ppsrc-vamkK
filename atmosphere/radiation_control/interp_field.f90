! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Interpolation routine for radiation calculations

!  Purpose: To interpolate the fields from grid boxes where radiation
!           calculations were performed to grid boxes where the
!           calculations should have been done, but weren't because
!           of spatial degradation of calculations.

!  Method:  Quite straightforward. FIRST_DATA_INTERP contains the
!           information on where the first grid box, on the first row,
!           needs interpolation, and from the knowledge of a chequer-
!           board pattern we can locate in a simple arithmetic fashion
!           all points which need interpolation. Note that calculations
!           performed inside data part of the PE only (i.e. not halo).

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

      SUBROUTINE interp_field(es_space_interp,field,levels,             &
                              first_data_interp,row_length,num_rows,    &
                              first_row, last_row, offx, offy )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER  ,INTENT(IN) :: levels
!            The no. of levels in the input field
      INTEGER  ,INTENT(IN) :: row_length
!            The length of a row on the PE (control level)
      INTEGER  ,INTENT(IN) :: num_rows
!            The number of rows excluding halo and polar rows
      INTEGER  ,INTENT(IN) :: first_data_interp
!              The first data point, in data co-ords (ie excluding
!              the halo), which needs to be interpolated on a PE
      INTEGER  ,INTENT(IN) :: first_row
!            The first row that needs to be interpolated
      INTEGER  ,INTENT(IN) :: last_row
!            The last row that needs to be interpolated
      INTEGER  ,INTENT(IN) :: offx, offy
!            The halo widths of the input field

      REAL  ,INTENT(IN) :: es_space_interp(4, row_length, num_rows)
!            The coefficients for radiation interpolation

      REAL  ,INTENT(INOUT) :: field(1-offx:row_length+offx,             &
                                        1-offy:num_rows+offy,levels)
!            The field to be interpolated to neighbouring grid boxes

!     Local variables

      INTEGER :: i              ! loop variable
      INTEGER :: j              ! loop variable
      INTEGER :: k              ! loop variable
      INTEGER :: delta          ! temporary storage variable

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('INTERP_FIELD',zhook_in,zhook_handle)
      IF (first_row == 1) THEN
        delta=first_data_interp
      ELSE
        delta=MOD(first_data_interp+1,2)
      END IF

      DO k=1,levels
        DO j=first_row,last_row,2
          DO i=1+delta,row_length,2
            field(i,j,k)=es_space_interp(1,i,j)*field(i,j+1,k)          &
                  + es_space_interp(2,i,j)*field(i+1,j,k)               &
                  + es_space_interp(3,i,j)*field(i,j-1,k)               &
                  + es_space_interp(4,i,j)*field(i-1,j,k)
          END DO
        END DO
        DO j=first_row+1,last_row,2
          DO i=2-delta,row_length,2
            field(i,j,k)=es_space_interp(1,i,j)*field(i,j+1,k)          &
                  + es_space_interp(2,i,j)*field(i+1,j,k)               &
                  + es_space_interp(3,i,j)*field(i,j-1,k)               &
                  + es_space_interp(4,i,j)*field(i-1,j,k)
          END DO
        END DO
      END DO


      IF (lhook) CALL dr_hook('INTERP_FIELD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE interp_field
