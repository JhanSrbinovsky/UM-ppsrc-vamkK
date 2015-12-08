! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Force_equil
      Subroutine IDL_Force_equil(row_length, rows, model_levels         &
     &,                      off_x, off_y, p_zero, theta_star           &
     &,                      cool_rate, const_SST, timestep             &
     &,                      L_physics)


! Purpose: Provides radiative cooling for radiative-convective
!          equilibrium idealised problem
!
! Method:
!          Is described in ;
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, off_x                                                           &
     &, off_y

      Real                                                              &
     &  timestep                                                        &
     &, p_zero

      Real                                                              &
     &  theta_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &        model_levels)                                             &
     &, cool_rate                                                       &
     &, const_SST

      Logical                                                           &
     & L_physics
!
! Local variables
!
      Real                                                              &
     &  cool_rate_ir(model_levels)                                      &
     &, p_refn(model_levels)

      Integer                                                           &
     &  i,j,k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
! define cooling rate by IR radiation
!

      IF (lhook) CALL dr_hook('IDL_FORCE_EQUIL',zhook_in,zhook_handle)
      do k = 1, model_levels
        p_refn(k) = p_zero - (k - 0.5) * (p_zero/model_levels)
        if (p_refn(k)  <   20000.) then
          cool_rate_ir(k) = 0.0
        elseif (p_refn(k)  <   40000.) then
          cool_rate_ir(k) = cool_rate * (p_refn(k) - 20000.)
          cool_rate_ir(k) = cool_rate_ir(k)/20000.
        else
          cool_rate_ir(k) = cool_rate
        endif
!
! divide the cool rate (K/day) by 86400 to get the cool rate in K/second
!
        cool_rate_ir(k) = cool_rate_ir(k)/86400.
      enddo
!
! determine the new value of theta_star for all model levels
      do k = 1, model_levels
        do j = 1 , rows
          do i = 1 , row_length
            theta_star(i,j,k) = theta_star(i,j,k)                       &
     &                        + timestep * cool_rate_ir(k)
          enddo
        enddo
      enddo

!
! if physics not turned on keep the SST constant by keeping theta at
! k=1 constant.
!
      If(.not. L_physics) then
        do j = 1 , rows
          do i = 1 , row_length
            theta_star(i,j,1) = const_SST
          enddo
        enddo
      Endif

      IF (lhook) CALL dr_hook('IDL_FORCE_EQUIL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_Force_equil
