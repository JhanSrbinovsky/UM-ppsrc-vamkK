! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Two_Norm_levels

      Subroutine Two_Norm_levels(                                       &
     &                          field, row_length, rows, levels,        &
     &                          start_level, end_level,                 &
     &                          model_domain, halo_i, halo_j,           &
     &                          off_u, off_v,                           &
     &                          at_extremity, n_procx, n_procy,         &
     &                          global_row_length, global_rows,         &
     &                          gc_proc_col_group, gc_proc_row_group,   &
     &                          l_datastart, L_do_halos,                &
     &                          L_do_rims, rims_to_do, Two_Norm )

! Purpose:
!          Calculates the two norm over a range of levels of the
!          input field.
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
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_do_rims                                                       &
                       ! LAMs: include rims, Global do all polar points
     &, L_do_halos     ! include halos in calculation of norm
!                    NB this means that some points are counted twice

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, levels                                                          &
                         ! number of levels in field.
     &, start_level                                                     &
                         ! start_level for norm calculation
     &, end_level        ! end_level for norm calculation

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, halo_i                                                          &
     &, halo_j                                                          &
     &, off_u                                                           &
                         ! = 1 if field is at u-point
     &, off_v                                                           &
                         ! = 1 if field is at v-point
     &, rims_to_do                                                      &
     &, l_datastart(3)                                                  &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, gc_proc_col_group                                               &
     &, gc_proc_row_group                                               &
     &, n_procx                                                         &
     &, n_procy

      Real                                                              &
     &  field(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

! Arguments with Intent OUT. ie: variables Output only

      Real                                                              &
     &  Two_Norm

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
     &, i_start                                                         &
     &, i_stop                                                          &
     &, j_start                                                         &
     &, j_stop                                                          &
     &, rims2                                                           &
     &, ipoints                                                         &
     &, istat

      Real                                                              &
     &  points           ! number of points used in norm

! Local arrays for parallel code

      real                                                              &
     &  two_norm_component( row_length+2*halo_i, rows+2*halo_j )        &
     &, two_norm_temp( 1 )

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle




!  External Routines:
! ----------------------------------------------------------------------
! Section 1.   Calculate Norm.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('TWO_NORM_LEVELS',zhook_in,zhook_handle)

      rims2 = 2* rims_to_do
      i_start = 1
      i_stop = row_length
      j_start = 1
      j_stop = rows

      if ( L_do_rims ) then

        ipoints = global_row_length * global_rows

        if( L_do_halos ) then
          i_start = 1 - halo_i
          i_stop = row_length + halo_i
          j_start = 1 - halo_j
          j_stop = rows + halo_j
          ipoints = (global_row_length + 2 * n_procx * halo_i ) *       &
     &              (global_rows + 2 * n_procy * halo_j )
        endif ! L_do_halos

        If ( model_domain == mt_global ) Then
          if(at_extremity(PNorth)) j_stop = rows - off_v
          ipoints = global_row_length * (global_rows - off_v)
        Endif ! model_domain == mt_global

        If( model_domain == mt_lam ) then
          if(at_extremity(PNorth)) j_stop = rows - off_v
          if(at_extremity(PEast)) i_stop = row_length - off_u
          ipoints = (global_row_length - off_u) * (global_rows - off_v)
        End If ! model_domain == mt_lam

      else  ! do not do halos or rims

        if( L_do_halos ) then
          write(6,*)' You must do rims if halos are included in L2norms'
        endif ! L_do_halos

        If ( model_domain == mt_global ) Then
          if(at_extremity(PSouth)) j_start = 2
          if(at_extremity(PNorth)) j_stop = rows - 1
          ipoints = global_row_length * (global_rows - 2) + 2
        Endif ! model_domain == mt_global

        If( model_domain == mt_lam ) then
          if(at_extremity(PSouth)) j_start = 1 + rims_to_do
          if(at_extremity(PNorth)) j_stop = rows - rims_to_do
          if(at_extremity(PEast)) i_stop = row_length - rims_to_do
          if(at_extremity(PWest)) i_start = 1 + rims_to_do
          ipoints = (global_row_length - rims2) *                       &
     &              (global_rows  - rims2)

        Elseif (model_domain == mt_cyclic_LAM) then
          if(at_extremity(PSouth)) j_start = 1 + rims_to_do
          if(at_extremity(PNorth)) j_stop = rows - rims_to_do
          ipoints = global_row_length * (global_rows - rims2)

        End If ! model_domain == mt_lam

      endif ! L_do_rims

      points = float(ipoints)

      two_norm_component = 0.0

      If ( model_domain == mt_Global .and. .not. L_do_rims) Then
! Global model only calculate norm for one of the polar points
!  if L_do_rims is .false.

        If(at_extremity(PSouth) .and. (l_datastart(1) == 1)) then
          Do k = start_level, end_level
            two_norm_component(1,1)= two_norm_component(1,1) +          &
     &                               field(1,1,k) * field(1,1,k)
          End Do
        End If
        If(at_extremity(PNorth) .and. (l_datastart(1) == 1)) then
          Do k = start_level, end_level
            two_norm_component(1,rows)= two_norm_component(1,rows) +    &
     &                             field(1,rows,k) * field(1,rows,k)
          End Do
        End If
      End If  !  model_domain == mt_Global .and. .not. L_do_rims

      Do k = start_level, end_level
        Do j = j_start , j_stop
          Do i = i_start, i_stop
            two_norm_component( i + halo_i, j + halo_j ) =              &
     &          two_norm_component( i + halo_i, j + halo_j ) +          &
     &                                 field(i,j,k) * field(i,j,k)
          End Do
        End Do
      End Do  ! k = start_level, end_level

      ! note halos included in sum.
      CALL global_2d_sums(two_norm_component, row_length+2*halo_i,      &
                          rows + 2*halo_j, 0, 0, 1, two_norm_temp)

! Normalize by dividing by number of points

      Two_Norm = SQRT (Two_Norm_temp(1) / points)

      IF (lhook) CALL dr_hook('TWO_NORM_LEVELS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Two_Norm_levels
