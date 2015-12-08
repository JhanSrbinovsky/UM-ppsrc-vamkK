! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolates the v-component from the v-grid to the p-grid.
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output

      Subroutine lbc_v_to_p ( v,                                        &
                              v_p,                                      &
                              row_length,                               &
                              rows,                                     &
                              v_rows,                                   &
                              offx,                                     &
                              offy,                                     &
                              L_EG_dump )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit NONE
!
! Description:
!   Interpolates the v-component from the v-grid to the p-grid.
!
! Method:
!   Linear Interpolation of v-comp from v-grid to p-grid.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine arguments

      Integer :: row_length    ! No of points
      Integer :: rows          ! No of p rows
      Integer :: v_rows        ! No of v rows
      Integer :: offx, offy    ! Halo sizes

!     v   : v-component on v-grid
!     v_p : v-component on p-grid

      Real    :: v  (1-offx:row_length+offx, 1-offy:v_rows+offy)
      Real    :: v_p(1-offx:row_length+offx, 1-offy:  rows+offy)
      
      LOGICAL :: L_EG_dump

! Local scalars:

      Integer :: i,j      ! Loop indices

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header

!     -----------------------------------
!     Interpolate v from v-grid to p-grid
!     -----------------------------------

      IF (lhook) CALL dr_hook('LBC_V_TO_P',zhook_in,zhook_handle)
      
      IF (L_EG_dump) THEN
      
        DO i = 1-offx, row_length+offx
      
          DO j = 1-offy, rows+offy
            v_p(i,j) = ( v(i,j) + v(i,j+1) ) * 0.5
          END DO
        
        END DO
      
      ELSE
      
        DO i = 1-offx, row_length+offx

          DO j = 1-offy+1, v_rows+offy
            v_p(i,j) = ( v(i,j-1) +  v(i,j) ) * 0.5
          END DO

          ! Copy for bottom row of v_p
          v_p(i,1-offy)    = v(i,1-offy)

          ! Copy for top row of v_p
          IF (v_rows < rows) THEN
            v_p(i,rows+offy) = v(i,v_rows+offy)
          END IF

        END DO
      
      END IF

      IF (lhook) CALL dr_hook('LBC_V_TO_P',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE lbc_v_to_p
