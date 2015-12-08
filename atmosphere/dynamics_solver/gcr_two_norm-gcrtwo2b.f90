! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_Two_Norm_2B

      SUBROUTINE GCR_Two_Norm_2B(                                       &
     &                        field, row_length, rows,                  &
     &                        model_levels, model_domain, Two_Norm,     &
     &                        offx, offy, at_extremity, n_proc,         &
     &                        gc_proc_col_group, gc_proc_row_group,     &
     &                        number_dif_horiz_points, l_datastart,     &
     &                        i_start, i_stop, j_begin, j_end           &

     &                        )

! Purpose:
!          Calculates the two norm of the input field divided by the
!          number of points in the field.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE global_2d_sums_mod, ONLY : global_2d_sums

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels     ! number of model levels.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

       integer                                                          &
     &  offx, offy                                                      &
     &, gc_proc_col_group                                               &
     &, gc_proc_row_group                                               &
     &, n_proc                                                          &
     &, l_datastart(3)                                                  &
     &, number_dif_horiz_points                                         &
                                     ! domain size set in PE_Helmholtz
     &, i_start, i_stop                                                 &
                                     ! loop bounds set in PE_Helmholtz
     &, j_begin, j_end                                                  &
                                     ! loop bounds set in PE_Helmholtz
     &, istat

      Real                                                              &
     &  field(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
                         !field to find norm of.

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

! Arguments with Intent OUT. ie: variables Output only

      Real                                                              &
     &  Two_Norm

! Local Variables.

      Integer                                                           &
     &  i, j, k

      Real                                                              &
     &  points

! Local arrays for parallel code
      real                                                              &
     &    two_norm_component(row_length,rows)                           &
     &,   two_norm_temp(1)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Calculate Norm.
! ----------------------------------------------------------------------


      IF (lhook) CALL dr_hook('GCR_TWO_NORM_2B',zhook_in,zhook_handle)
      Do j=1,rows
        Do i=1,row_length
          two_norm_component(i,j)=0.0
        End Do
      End Do

      If ( model_domain == mt_Global .and. l_datastart(1) == 1 ) Then
! Global model only calculate norm for one of the polar points.

        If ( at_extremity(PSouth) ) then
          Do k = 1, model_levels
            two_norm_component(1,1) = two_norm_component(1,1) +         &
     &                                field(1,1,k) * field(1,1,k)
          End Do
        End If

        If ( at_extremity(PNorth) ) then
          Do k = 1, model_levels
            two_norm_component(1,rows) = two_norm_component(1,rows) +   &
     &                             field(1,rows,k) * field(1,rows,k)
          End Do
        End If

      End If  !  model_domain == mt_Global .and. l_datastart(1) == 1

      Do k = 1, model_levels
        Do j = j_begin, j_end
          Do i = i_start, i_stop
          two_norm_component(i,j) = two_norm_component(i,j) +           &
     &                              field(i,j,k) * field(i,j,k)
          End Do
        End Do
      End Do ! k = 1, model_levels

      points = model_levels * number_dif_horiz_points

      CALL global_2d_sums(two_norm_component, row_length, rows, 0, 0,   &
                          1, two_norm_temp)

! Normalize norm by dividing by number of points at which it was
! calculated.

      Two_Norm = SQRT ( two_norm_temp(1) / points )

      IF (lhook) CALL dr_hook('GCR_TWO_NORM_2B',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_Two_Norm_2B
