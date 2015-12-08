! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolates u-component from the u-grid to the p-grid.
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output

      Subroutine lbc_u_to_p ( u,                                        &
                              u_p,                                      &
                              row_length,                               &
                              rows,                                     &
                              offx,                                     &
                              offy,                                     &
                              L_EG_dump )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit NONE
!
! Description:
!   Interpolates u-component from the u-grid to the p-grid.
!
! Method:
!   Linear Interpolation of u-comp from u-grid to p-grid.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine arguments

      Integer :: row_length    ! No of points
      Integer :: rows          ! No of rows
      Integer :: offx, offy    ! Halo sizes

!     u   : u-component on u-grid
!     u_p : u-component on p-grid

      Real    :: u  (1-offx:row_length+offx, 1-offy:rows+offy)
      Real    :: u_p(1-offx:row_length+offx, 1-offy:rows+offy)
      
      LOGICAL :: L_EG_dump

! Local scalars:

      Integer :: i,j    !  Loop indices

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!- End of header

!     ----------------------------------------------
!     Interpolate u from u-grid to p-grid ; u to u_p
!     ----------------------------------------------

      IF (lhook) CALL dr_hook('LBC_U_TO_P',zhook_in,zhook_handle)
      
      IF (L_EG_dump) THEN
      
        DO j = 1-offy, rows+offy
      
          DO i = 1-offx, row_length+offx-1
            u_p(i,j) = ( u(i,j) + u(i+1,j) ) * 0.5
          END DO
      
          ! Copy last column of U to last column of P
          u_p(row_length+offx,j) = u(row_length+offx,j)
        
        END DO
      
      ELSE
      
        DO j = 1-offy, rows+offy

          DO i = 1-offx+1, row_length+offx
            u_p(i,j) = ( u(i-1,j) +  u(i,j) ) * 0.5
          END DO

          ! Copy for first column of P points.
          u_p(1-offx,j) = u(1-offx,j)

        END DO
      
      END IF

      IF (lhook) CALL dr_hook('LBC_U_TO_P',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE lbc_u_to_p
