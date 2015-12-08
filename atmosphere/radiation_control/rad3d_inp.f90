! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to complete 3D radiation fields.
!
! Subroutine Interface
!
      Subroutine rad3d_inp(                                             &
     &  L_complete_North, L_complete_South, L_complete_deg,             &
     &  row_length, rows, off_x, off_y, first_row, last_row,            &
     &  first_data_interp, ES_space_interp,                             &
     &  levels,                                                         &
     &  RadField                                                        &
     &  )

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      IMPLICIT NONE
!
! Description: At points where full radiation calculations were not
!   done, radiation fields are completed by copying or interpolation.
!   This routine operates on 3D horizontal fields.
!
! Method: Radiation fields are calculated only at the first point on
!   the northern or southern polar rows in a full global configuration:
!   they must be copied to the other points along the row. When spatial
!   degradation is in operation, fields are calculated only at
!   alternate points and must be completed by interpolation.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: radiation control
!
! Declarations:
!
! Global variables:
!
! Arguments of the subroutine
!
      Logical, Intent(IN) :: L_complete_North
!                              Completion required at North pole
      Logical, Intent(IN) :: L_complete_South
!                              Completion required at South pole
      Logical, Intent(IN) :: L_complete_deg
!                              Completion by interpolation is
!                              required at points where radiation
!                              was not calculated because of
!                              spatial degradation
!
      Integer, Intent(IN) :: row_length
!                              Number of points along a row
!                              in the domain
      Integer, Intent(IN) :: rows
!                              Number of rows in the domain
      Integer, Intent(IN) :: levels
!                              Number of vertical levels in the domain
      Integer, Intent(IN) :: off_x
!                              Size of small halo in i-direction
      Integer, Intent(IN) :: off_y
!                              Size of small halo in j-direction
!
      Integer, Intent(IN) :: first_row
!                              First row requiring interpolation
!                              under spatial degradation
      Integer, Intent(IN) :: last_row
!                              Last row requiring interpolation
!                              under spatial degradation
      Integer, Intent(IN) :: first_data_interp
!                              The first data point, in data
!                              coordinates, which needs needs to be
!                              interpolated in the local domain
!
      Real, Intent(IN)    :: ES_space_interp(4, row_length, rows)
!                              Coefficients for spatial interpolation
!                              of radiation quantities
      Real, Intent(INOUT) :: RadField(row_length, rows, levels)
!                              Radiation field to be completed
!
!
!     Local Variables:
!
!
      Integer :: i
!                              Loop variable
      Integer :: j
!                              Loop variable
      Integer :: k
!                              Loop variable
      Real    :: RadField_inp(1-off_x:row_length+off_x,                 &
     &                     1-off_y:rows+off_y, levels)
!                              Field for interpolation including halos

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!- End of header
!
      IF (lhook) CALL dr_hook('RAD3D_INP',zhook_in,zhook_handle)
!
!     Fill in at the north and south poles as required.
      If ( L_complete_North ) Then
        Do k = 1, levels
          Do i = 2, row_length
            RadField(i,rows,k) = RadField(1,rows,k)
          End Do
        End Do
      End If
      If ( L_complete_South ) Then
        Do k = 1, levels
          Do i = 2, row_length
            RadField(i,1,k) = RadField(1,1,k)
          End Do
        End Do
      End If
!
!
!     Fill in fields when spatial degradation is in operation.
      If ( L_complete_deg  ) Then
!
        Do k = 1, levels
          Do j = 1, rows
            Do i = 1, row_length
              RadField_inp(i,j,k) = RadField(i,j,k)
            Enddo
          Enddo
        Enddo
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(RadField_inp, row_length, rows,                &
     &    levels, off_x, off_y, fld_type_p, .FALSE.)

! DEPENDS ON: set_external_halos
        CALL SET_EXTERNAL_HALOS(RadField_inp,row_length, rows,          &
     &                          levels,off_x,off_y,0.0)

! DEPENDS ON: interp_field
        Call Interp_Field(ES_space_interp, RadField_inp,                &
     &    levels, first_data_interp, row_length, rows,                  &
     &    first_row, last_row, off_x, off_y)
!
        Do k = 1, levels
          Do j = 1, rows
            Do i = 1, row_length
              RadField(i,j,k) = RadField_inp(i,j,k)
            Enddo
          Enddo
        Enddo
!
      Endif
!
      IF (lhook) CALL dr_hook('RAD3D_INP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE rad3d_inp
