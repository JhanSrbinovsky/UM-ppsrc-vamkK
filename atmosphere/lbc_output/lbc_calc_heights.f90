! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate heights at LBC grid points.
!
! Subroutine Interface:

      Subroutine LBC_Calc_Heights (                                     &
     &           lbc_points, model_levels, height_gen_method,           &
     &           first_constant_r_rho_level, z_top_of_model,            &
     &           eta_theta_levels, eta_rho_levels, orog_lbcs,           &
     &           r_theta_levels, r_rho_levels                           &
     & )

      USE earth_constants_mod, ONLY: earth_radius

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE
!
! Description:
!   Calculates a height field for a given set of LBC points
!
! Method:
!   Follows the same algorithm to calculate heights as
!   in subroutine SETCONA. Orography interpolated to the LBC
!   points is used. The height is also derived for the extra rho
!   level above the top theta level.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

      Integer :: lbc_points
      Integer :: model_levels
      Integer :: height_gen_method
      Integer :: first_constant_r_rho_level

      Real    :: z_top_of_model
      Real    :: eta_theta_levels (0:model_levels)
      Real    :: eta_rho_levels (model_levels)
      Real    :: orog_lbcs (lbc_points)

      Real    :: r_theta_levels (lbc_points, 0:model_levels  )
      Real    :: r_rho_levels   (lbc_points,   model_levels+1)

! Local variables

      Integer, Parameter :: height_gen_smooth = 2
      Character (Len=*), Parameter :: RoutineName = 'LBC_Calc_Heights'

      Integer            :: ErrorStatus
      Integer            :: i,k
      Character (Len=80) :: Cmessage

      Real    :: r_ref_theta (0:model_levels)
      Real    :: r_ref_rho   (  model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('LBC_CALC_HEIGHTS',zhook_in,zhook_handle)
      ErrorStatus = 0
      Cmessage = ' '

! ----------------------------
! Set reference height profile
! ----------------------------

      Do k = 1, model_levels
        r_ref_theta(k) = eta_theta_levels(k) * z_top_of_model
        r_ref_rho(k)   = eta_rho_levels(k)   * z_top_of_model
      End Do

! -------------------------------
! Set bottom level, ie: orography
! -------------------------------

      Do i = 1, lbc_points
          r_theta_levels(i,0) = Orog_lbcs(i) + Earth_radius
      End Do

! ----------------------------------------------------------
! Set up heights for levels below first_constant_r_rho_level
! ----------------------------------------------------------

      Select Case ( height_gen_method)

        Case ( height_gen_smooth)

! Smooth quadratic height generation

          Do k = 1, first_constant_r_rho_level-1
            Do i = 1, lbc_points

              r_rho_levels(i,k) = r_ref_rho(k) + Earth_radius +         &
     &        Orog_lbcs(i) * (1.0 -eta_rho_levels(k) /                  &
     &        eta_rho_levels(first_constant_r_rho_level) ) ** 2

              r_theta_levels(i,k) = r_ref_theta(k) + Earth_radius +     &
     &        Orog_lbcs(i) * (1.0 -eta_theta_levels(k) /                &
     &        eta_rho_levels(first_constant_r_rho_level) ) ** 2

            End Do
          End Do

        Case Default

          ErrorStatus = 10
          Write (Cmessage,*) 'Height generation method ',               &
     &    height_gen_method,' not recognised.'

          Call Ereport ( RoutineName, ErrorStatus, Cmessage)

        End Select

! -------------------------------------------------------
! For constant levels set r to be a constant on the level
! -------------------------------------------------------

      Do k = first_constant_r_rho_level, model_levels
        Do i = 1, lbc_points
          r_theta_levels(i,k) = Earth_radius + r_ref_theta(k)
          r_rho_levels(i,k)   = Earth_radius + r_ref_rho(k)
        End Do
      End Do

! ---------------------------------------------------
! Determine height for top rho level (model_levels+1)
! ---------------------------------------------------

      Do i = 1, lbc_points

        r_rho_levels(i,model_levels+1) =                                &
     &  r_theta_levels(i,model_levels) +                                &
     &    ( r_theta_levels(i,model_levels) -                            &
     &      r_rho_levels(i,model_levels) )

      End Do

      IF (lhook) CALL dr_hook('LBC_CALC_HEIGHTS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LBC_Calc_Heights
