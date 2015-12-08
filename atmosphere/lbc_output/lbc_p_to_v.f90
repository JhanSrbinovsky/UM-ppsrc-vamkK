! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolates v-component from the p-grid to the v-grid.
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output

      Subroutine lbc_p_to_v ( v_p,                                      &
                              v,                                        &
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
!   Interpolates v-component from the p-grid to the v-grid.
!
! Method:
!   Linear Interpolation of v-comp from p-grid to v-grid.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine arguments

      Integer :: row_length    !  No of points
      Integer :: rows          !  No of p rows
      Integer :: v_rows        !  No of v rows
      Integer :: offx, offy    !  Halo sizes

!     v_p : v-component on p-grid
!     v   : v-component on v-grid

      Real    :: v_p(1-offx:row_length+offx, 1-offy:  rows+offy)
      Real    :: v  (1-offx:row_length+offx, 1-offy:v_rows+offy)
      
      LOGICAL :: L_EG_dump

! Local scalars:

      Integer :: i,j           !  Loop indices

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!- End of header

!     ----------------------------------------------
!     Interpolate v from p-grid to v-grid ; v_p to v
!     ----------------------------------------------

      IF (lhook) CALL dr_hook('LBC_P_TO_V',zhook_in,zhook_handle)
      
      IF (L_EG_dump) THEN
      
        DO i = 1-offx, row_length+offx
      
          DO j = 1-offy+1, v_rows+offy-1
            v(i,j) = ( v_p(i,j-1) + v_p(i,j) ) * 0.5
          END DO
        
          ! Copy to first and last rows of V
          v(i,1-offy)     = v_p(i,1-offy)
          v(i,v_rows+offy) = v_p(i,rows+offy)
      
        END DO
      
      ELSE
      
        IF (v_rows < rows ) THEN

          DO i = 1-offx, row_length+offx
            DO j = 1-offy, v_rows+offy
              v(i,j) = ( v_p(i,j) +  v_p(i,j+1) ) * 0.5
            END DO
          END DO

        ELSE  !  v_rows == rows

          DO i = 1-offx, row_length+offx
            DO j = 1-offy, v_rows+offy-1
              v(i,j) = ( v_p(i,j) +  v_p(i,j+1) ) * 0.5
            END DO
            v(i,v_rows+offy) = v_p(i,rows+offy)
          END DO

        END IF
      
      END IF

      IF (lhook) CALL dr_hook('LBC_P_TO_V',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE lbc_p_to_v
