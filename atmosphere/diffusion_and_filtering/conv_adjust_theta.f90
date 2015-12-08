! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine conv_adjust_theta

      Subroutine conv_adjust_theta                                      &
     &                     (theta, Exner,                               &
     &                      off_x, off_y,                               &
     &                      rows, row_length, model_levels,             &
     &                      start_level, end_level)

! Purpose:
!          Does a conservative re-arrangement of theta if
!          theta(k) < theta(k-1)
!          Potential energy in the layer is conserved
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, start_level                                                     &
     &, end_level

      Real, Intent(In) ::                                               &
     &  Exner (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels + 1)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real, Intent(InOut) ::                                            &
     &  theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &         model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k               ! Loop indices

      Real                                                              &
     &  exner_ratio

! Local arrays

      Real                                                              &
     &  dthetadz(row_length, rows)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Conservative convective adjustment
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CONV_ADJUST_THETA',zhook_in,zhook_handle)
      Do k = start_level, end_level
        Do j = 1, rows
          Do i = 1, row_length
            dthetadz(i,j) = theta(i,j,k) - theta(i,j,k-1)
             if( dthetadz(i,j) < 0.0 ) then
! If static stability is less than 0 (i.e. theta(k) < theta(k-1))
! then re-arrange and conserve potential energy
              exner_ratio = exner(i,j,k) / exner(i,j,k-1)
              theta(i,j,k-1) = theta(i,j,k-1) + dthetadz(i,j) /         &
     &                            ( 1.0 + 1.0/exner_ratio )
              theta(i,j,k) = theta(i,j,k) - dthetadz(i,j) /             &
     &                            ( 1.0 + exner_ratio )
             endif  ! dthetadz(i,j) < 0.0
           End Do
        End Do
      End Do     ! k = start_level, end_level

      IF (lhook) CALL dr_hook('CONV_ADJUST_THETA',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE conv_adjust_theta
