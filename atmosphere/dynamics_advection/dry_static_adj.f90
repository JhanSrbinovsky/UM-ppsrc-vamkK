! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Dry_static_adj

      Subroutine Dry_static_adj(                                        &
     &                          theta, q,                               &
     &                          rows, row_length, model_levels,         &
     &                          wet_model_levels, off_x, off_y,         &
     &                          diagnostics,                            &
     &                          L_adjust_wet)

! Purpose:
!          Enforces static stability on input theta field.
!          Uses very simple sort algorithm resulting in no mixing of
!          temperature just a simple re-arrangement.
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE atmos_constants_mod, ONLY: recip_epsilon

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_model_levels                                                &
                         ! number of model levels with moist variables.
     &, diagnostics                                                     &
                         ! switch controlling printing of unstable pnts
     &, off_x                                                           &
                     ! Size of small halo in i
     &, off_y        ! Size of small halo in j.

      Logical                                                           &
     &  L_adjust_wet     ! true if code to adjust moist column as well.

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, q (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     wet_model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k, list                                                   &
                          ! Loop indices
     &, count                                                           &
     &, count_old                                                       &
     &, iterations                                                      &
     &, n_search                                                        &
     &, i_search                                                        &
     &, n_match

      Real                                                              &
     &  temp                                                            &
     &, thetav(row_length, rows, model_levels)

! Local arrays

      Integer                                                           &
     &  Unstable(row_length* rows*(model_levels-1))                     &
     &, Unstable_old(row_length* rows*(model_levels-1))                 &
     &, i_prev(row_length* rows)                                        &
     &, j_prev(row_length* rows)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

! ----------------------------------------------------------------------
! Section 1. Check for inertial instability and remove.
!            Iterative sort procedure is used.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('DRY_STATIC_ADJ',zhook_in,zhook_handle)
      iterations = 0
      count = 0

      Do k = 2, model_levels
        Do j = 1, rows
          Do i = 1, row_length

            If (theta(i,j,k)  <   theta(i,j,k-1)) Then
              temp = theta(i,j,k)
              theta(i,j,k) = theta(i,j,k-1)
              theta(i,j,k-1) = temp
              count = count + 1
              unstable(count) = k*1e6+j*1000 + i
            End If
          End Do
        End Do
      End Do

      If (diagnostics  >   0 .and. count  >   0) Then
          print*,' iteration number ',iterations,                       &
     &           ' dry instability at ',count,' points '
        If (diagnostics  >   1) Then
          Do list = 1, count
            k = unstable(list) * 1e-6
            j = (unstable(list) - k * 1e6) *1e-3
            i = unstable(list) - k * 1e6 - j*1e3
            print*,' dry instability at ',i,j,k
          End Do
        End If
      End If

      Do while (count  >   0 .and. iterations  <   20)
        count_old = count
        Do list = 1, count
          unstable_old(list) = unstable(list)
        End Do
        count = 0
        iterations = iterations + 1

        n_search = 0
        n_match = 0
        Do list = 1, count_old
          k = unstable_old(list) * 1e-6
          j = (unstable_old(list) - k * 1e6) *1e-3
          i = unstable_old(list) - k * 1e6 - j*1e3
          Do i_search = 1, n_search
            If (i  ==  i_prev(i_search) .and.                           &
     &          j  ==  j_prev(i_search)) Then
              n_match = 1
            End If
          End Do
          If (n_match  ==  0) Then
            n_search = n_search + 1
            i_prev(n_search) = i
            j_prev(n_search) = j
          End If
          n_match = 0
        End Do

        Do list = 1, n_search
          i = i_prev(list)
          j = j_prev(list)
          Do k = 2, model_levels
            If (theta(i,j,k)  <   theta(i,j,k-1)) Then
              temp = theta(i,j,k)
              theta(i,j,k) = theta(i,j,k-1)
              theta(i,j,k-1) = temp
              count = count + 1
              unstable(count) = k*1e6+j*1000 + i
            End If
          End Do
        End Do

        If (diagnostics  >   1 .and. count  >   0) Then
          print*,' iteration number ',iterations
          Do list = 1, count
            k = unstable(list) * 1e-6
            j = (unstable(list) - k * 1e6) *1e-3
            i = unstable(list) - k * 1e6 - j*1e3
            print*,' dry instability at ',i,j,k
          End Do
        End If

! end do while
      End Do

      If (iterations  ==  20 .and. count  /=  0) Then
        print*,' unable to remove static instability. '
      End If

! ----------------------------------------------------------------------
! Section 2. Check for inertial instability and remove.
!            Iterative sort procedure is used.
! ----------------------------------------------------------------------

      If (L_adjust_wet) Then
        iterations = 0
        count = 0

        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
!ajm              thetav(i,j,k) = theta(i,j,k) * (1. + q(i,j,k) *
!ajm     &                                             recip_epsilon )
!ajm     &                                     / (1. + q(i,j,k))
              thetav(i,j,k) = theta(i,j,k) *                            &
     &        ( 1. +  ( (recip_epsilon -1.)* q(i,j,k) ) )
            End Do
          End Do
        End Do

        Do k = 1 + wet_model_levels, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              thetav(i,j,k) = theta(i,j,k)
            End Do
          End Do
        End Do

        Do k = 2, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              If (thetav(i,j,k)  <   thetav(i,j,k-1)) Then
                temp = thetav(i,j,k)
                thetav(i,j,k) = thetav(i,j,k-1)
                thetav(i,j,k-1) = temp
                count = count + 1
                unstable(count) = k*1e6+j*1000 + i
              End If
            End Do
          End Do
        End Do

        If (diagnostics  >   0 .and. count  >   0) Then
          print*,' iteration number ',iterations,                       &
     &           ' moist instability at ',count,' points '
          If (diagnostics  >   1 ) Then
            Do list = 1, count
              k = unstable(list) * 1e-6
              j = (unstable(list) - k * 1e6) *1e-3
              i = unstable(list) - k * 1e6 - j*1e3
              print*,' moist instability at ',i,j,k
            End Do
          End If
        End If

        Do while (count  >   0 .and. iterations  <   20)
          count_old = count
          Do list = 1, count
            unstable_old(list) = unstable(list)
          End Do
          count = 0
          iterations = iterations + 1

          n_search = 0
          n_match = 0
          Do list = 1, count_old
            k = unstable_old(list) * 1e-6
            j = (unstable_old(list) - k * 1e6) *1e-3
            i = unstable_old(list) - k * 1e6 - j*1e3
            Do i_search = 1, n_search
              If (i  ==  i_prev(i_search) .and.                         &
     &          j  ==  j_prev(i_search)) Then
                n_match = 1
              End If
            End Do
            If (n_match  ==  0) Then
              n_search = n_search + 1
              i_prev(n_search) = i
              j_prev(n_search) = j
            End If
            n_match = 0
          End Do

          Do list = 1, n_search
            i = i_prev(list)
            j = j_prev(list)
            Do k = 2, model_levels
              If (thetav(i,j,k)  <   thetav(i,j,k-1)) Then
                temp = thetav(i,j,k)
                thetav(i,j,k) = thetav(i,j,k-1)
                thetav(i,j,k-1) = temp
                count = count + 1
                unstable(count) = k*1e6+j*1000 + i
              End If
            End Do
          End Do

          If (diagnostics  >   1 .and. count  >   0) Then
            print*,' iteration number ',iterations
            Do list = 1, count
              k = unstable(list) * 1e-6
              j = (unstable(list) - k * 1e6) *1e-3
              i = unstable(list) - k * 1e6 - j*1e3
              print*,' wet instability at ',i,j,k
            End Do
          End If

! end do while
        End Do

        If (iterations  ==  20 .and. count  /=  0) Then
          print*,' unable to remove moist static instability. '
        End If

! now diagnose q from theta and thetav
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              temp = thetav(i,j,k) / theta(i,j,k)
              If (temp  <=  1.0 ) Then
                temp = 1.0
              End If
!ajm              q(i,j,k) = (1. - temp) / (temp - recip_epsilon )
              q(i,j,k) =  (1. - temp) / (1. - recip_epsilon )
            End Do
          End Do
        End Do

      End If

! End of routine
      IF (lhook) CALL dr_hook('DRY_STATIC_ADJ',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Dry_static_adj

