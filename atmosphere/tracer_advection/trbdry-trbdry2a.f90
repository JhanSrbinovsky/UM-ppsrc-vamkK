! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINE TRBDRY-------------------------------------------------
!LL
!LL  Purpose: Special routine to add psuedo source terms to boundary
!LL           data in limited area model.
!LL  Method:  Sets the boundary of aerosol using the
!LL           function BDRYV. This is then copied outward to fill
!LL           the halo. The function BDRYV
!LL           is specific to a model configuration: the current
!LL           version (5.2) is specific to UK MES.
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3,
!LL                        Version 7, dated 11/3/93.
!LL
!LL
!*L  Arguments:---------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Tracer Advection
      Subroutine Trbdry(                                                &
     & row_length, rows, n_rows, model_levels,                          &
     & offx, offy, at_extremity,                                        &
     & pressure,                                                        &
     & u, v,                                                            &
     & murk, timestep                                                   &
     &)

      USE conversions_mod, ONLY: recip_pi_over_180

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                          ! Length of a model row
     &, rows                                                            &
                          ! Number of rows for theta,u fields
     &, n_rows                                                          &
                          ! Number of rows for v fields
     &, model_levels                                                    &
                          ! Number of model levels
     &, offx                                                            &
                          ! Size of "single" halo (EW direction)
     &, offy                                                            &
                          ! Size of "single" halo (NS direction)
     &, timestep          ! Timestep in seconds

      Logical, Intent(In) ::                                            &
     &  at_extremity(4)   ! At an edge?

      Real, Intent(In) ::                                               &
     &  u( 1 - offx : row_length + offx,                                &
                                                  ! U wind
     &     1 - offy : rows + offy,                                      &
     &     model_levels )                                               &
     &, v( 1 - offx : row_length + offx,                                &
                                                  ! V wind
     &     1 - offy : n_rows + offy,                                    &
     &     model_levels )                                               &
     &, pressure( 1 - offx : row_length + offx,                         &
                                                  ! Pressure on
     &            1 - offy : rows + offy,                               &
                                                  ! Rho levels
     &            model_levels)

      Real, Intent(InOut) ::                                            &
     &  murk( 1 - offx : row_length + offx,                             &
                                                ! Aerosol for which the
     &        1 - offy : rows + offy,                                   &
                                                ! boundary will be set
     &        model_levels )


!* Local, including SAVE'd, storage------------------------------------

!  (a) Scalars

      Real ::                                                           &
     & dc, wdir, wspeed, press

      Real ::                                                           &
     & BdryV        ! FUNCTION giving boundary value.

!  (b) Others.
      Integer  ::                                                       &
     & i, j, jj, k  ! Loop counters

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
!L Loop across Northern row.
!-----------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('TRBDRY',zhook_in,zhook_handle)
      If ( at_extremity( PNorth ) ) Then
        Do k = 1, model_levels
          j = rows
          jj = n_rows
          Do i = 1, row_length

            If ( v(i, jj, k) > 0. ) Then       ! Outflow
              murk(i, j, k) = murk(i, j-1, k)
            Else                               ! Inflow
              press = pressure(i, j, k)
              wdir = Atan2( v(i, jj, k), u(i, j, k) ) *                 &
     &                                                Recip_Pi_Over_180
              wspeed = Sqrt( u(i, j, k) * u(i, j, k) +                  &
     &                       v(i, jj, k) * v(i, jj, k) )

! DEPENDS ON: bdryv
              murk(i, j, k) = Bdryv( wdir, wspeed, press )
            End If
          End Do
        End Do
      End If



!-----------------------------------------------------------------------
!L Loop across Southern row.
!-----------------------------------------------------------------------
!
      If ( at_extremity( PSouth ) ) Then
        Do k = 1, model_levels
          j = 1
          Do i = 1, row_length

            If ( v(i, j, k) < 0. ) Then             ! Outflow
              murk(i, j, k) = murk(i, j+1, k)
            Else                                    ! Inflow
              press = pressure(i, j, k)
              wdir = Atan2( v(i, j, k), u(i, j, k) ) * Recip_Pi_Over_180
              wspeed = Sqrt( u(i, j, k) * u(i, j, k) +                  &
     &                       v(i, j, k) * v(i, j, k) )

! DEPENDS ON: bdryv
              murk(i, j, k) = Bdryv( wdir, wspeed, press )
            End If
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
!L Loop across Western column
!-----------------------------------------------------------------------
!
      If ( at_extremity( PWest ) ) Then
        Do k = 1, model_levels
          Do j = 1, rows

            ! jj is used for v indexing
            If ( j > n_rows + offy ) Then
              jj = n_rows +offy
            Else
              jj = j
            End If

            i = 1

            If ( u(i,j,k) < 0. ) Then           ! Outflow
              murk(i, j, k) = murk(i+1, j, k)
            Else                                ! Inflow
              press = pressure(i, j, k)
              wdir = Atan2( v(i, jj, k), u(i, j, k) ) *                 &
     &                                                Recip_Pi_Over_180
              wspeed = Sqrt( u(i, j, k) * u(i, j, k) +                  &
     &                       v(i, jj, k) * v(i, jj, k) )

! DEPENDS ON: bdryv
              murk(i, j, k) = Bdryv( wdir, wspeed, press )
            End If
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
!L Loop across Eastern column
!-----------------------------------------------------------------------
!
      If ( at_extremity( PEast ) ) Then
        Do k = 1, model_levels
          Do j = 1, rows

            ! jj is used for v indexing
            If ( j > n_rows + offy ) Then
              jj = n_rows +offy
            Else
              jj = j
            End If

            i = row_length

            If ( u(i,j,k) > 0. ) Then          ! Outflow
              murk(i, j, k) = murk(i-1, j, k)
            Else                               ! Inflow
              press = pressure(i, j, k)
              wdir = Atan2( v(i, jj, k), u(i, j, k) ) *                 &
     &                                                Recip_Pi_Over_180
              wspeed = Sqrt( u(i, j, k) * u(i, j, k) +                  &
     &                       v(i, jj, k) * v(i, jj, k) )

! DEPENDS ON: bdryv
              murk(i, j, k) = Bdryv( wdir, wspeed, press )
            End If
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! Swap bounds to ensure halos full correctly
!-----------------------------------------------------------------------

! DEPENDS ON: swap_bounds
      Call Swap_Bounds( murk, row_length, rows, model_levels,           &
     &                  offx, offy, fld_type_p, .false. )

!-----------------------------------------------------------------------
!L Now need to fill full extended halos with copies of the
!L calculated data. Start at the North
!-----------------------------------------------------------------------

      If ( at_extremity( PNorth ) ) Then
        Do k = 1, model_levels
          Do j = rows + 1, rows + offy
            Do i = 1 - offx, row_length + offx
              murk(i, j, k) = murk(i, rows, k)
            End Do
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! Fill extended halos for Southern rows
!-----------------------------------------------------------------------

      If ( at_extremity( PSouth ) ) Then
        Do k = 1, model_levels
          Do j = 1 - offy, 0
            Do i = 1 - offx, row_length + offx
              murk(i, j, k) = murk(i, 1, k)
            End Do
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! Fill extended halos for Western columns
!-----------------------------------------------------------------------

      If ( at_extremity( PWest ) ) Then
        Do k = 1, model_levels
          Do j = 1 - offy, rows + offy
            Do i = 1 - offx, 0
              murk(i, j, k) = murk(1, j, k)
            End Do
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! Fill extended halos for Eastern columns
!-----------------------------------------------------------------------

      If ( at_extremity( PEast ) ) Then
        Do k = 1, model_levels
          Do j = 1 - offy, rows + offy
            Do i = row_length + 1, row_length + offx
              murk(i, j, k) = murk (row_length, j, k)
            End Do
          End Do
        End Do
      End If

      IF (lhook) CALL dr_hook('TRBDRY',zhook_out,zhook_handle)
      RETURN
      End Subroutine Trbdry


