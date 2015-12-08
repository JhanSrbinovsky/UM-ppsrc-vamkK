! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

MODULE trophgt1_mod

IMPLICIT NONE
!+ ---------------------------------------------------------------------
!  Description: Limits on the position of the tropopause used in
!               radiative calculations are defined.
!- ----------------------------------------------------------------------
!
!     Set limits on the height of the tropopause. Generous limits are
!     required as the tropopause can be very low in the Antarctic
!     winter. These limits should be reviewed for simulations of
!     climates very different from the present one, such as runs
!     with very high concentrations of CO2.
!
!     These limits are in part chosen to accord with the pressure
!     levels used in earlier configurations of the Unified Model.
!     In setting the upper limit consideration has been given to the
!     paper entitled "The tropical tropopause over the Western
!     Pacifiic: Wave driving, convection and the annual cycle,"
!     (J. Geophys. Res., 1996, Vol. 101 (D16), p. 21223) by
!     G. C. Reid and K. S. Gage.
      Real, parameter :: z_min_trop = 2.0e3 ! Lowest permitted height
      Real, parameter :: z_max_trop = 2.0e4 ! Maximum permitted height
!
! -----------------------------------------------------------------------

END MODULE trophgt1_mod
