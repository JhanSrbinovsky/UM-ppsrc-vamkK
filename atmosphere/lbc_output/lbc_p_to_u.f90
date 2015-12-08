! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolates u-component from the p-grid to the u-grid.
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output

      Subroutine lbc_p_to_u (u_p,                                       &
                             u,                                         &
                             row_length,                                &
                             rows,                                      &
                             offx,                                      &
                             offy,                                      &
                             L_EG_dump )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit NONE
!
! Description:
!   Interpolates u-component from the p-grid to the u-grid.
!
! Method:
!   Linear Interpolation of u-comp from p-grid to u-grid.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine arguments

      Integer :: row_length    ! No of points
      Integer :: rows          ! No of rows
      Integer :: offx, offy    ! Halo sizes

!     u_p : u-component on p-grid
!     u   : u-component on u-grid

      Real    :: u_p(1-offx:row_length+offx, 1-offy:rows+offy)
      Real    :: u  (1-offx:row_length+offx, 1-offy:rows+offy)
      
      LOGICAL :: L_EG_dump

! Local scalars:

      Integer :: i,j           ! Loop indices

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!- End of header

!     ----------------------------------------------
!     Interpolate u from p-grid to u-grid ; u_p to u
!     ----------------------------------------------

      IF (lhook) CALL dr_hook('LBC_P_TO_U',zhook_in,zhook_handle)
      
      IF (L_EG_dump) THEN
      
        DO j = 1-offy, rows+offy
      
          DO i = 1-offx+1, row_length+offx
            u(i,j) = ( u_p(i-1,j) + u_p(i,j) ) * 0.5
          END DO
        
          ! Copy to first U column
          u(1-offx,j) = u_p(1-offx,j)
        
        END DO
      
      ELSE
      
        DO j = 1-offy, rows+offy

          DO i = 1-offx, row_length+offx-1
            u(i,j) = ( u_p(i,j) +  u_p(i+1,j) ) * 0.5
          END DO

        END DO
      
      END IF

      IF (lhook) CALL dr_hook('LBC_P_TO_U',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE lbc_p_to_u
