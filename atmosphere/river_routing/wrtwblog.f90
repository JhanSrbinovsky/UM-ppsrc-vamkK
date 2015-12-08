! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE wrtwblog_mod

IMPLICIT NONE

CONTAINS



      SUBROUTINE wrtWBlog(iofile, iy, im, id, ih, ndev                  &
     &     , stoall, sto2all, dinall, doutall, drunall, drivall         &
     &     , dt)

      USE PrintStatus_mod
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!     write water balance monitoring
!
!     stoall    : Total river channel storage
!     sto2all   : Total river channel storage at the next time step
!     dinall    : Total inflow to the next grid for the next time step
!     doutall   : Total outflow to the next grid for the next time step
!   [ drunall*dt  : Total inflow to the grid from LSM]
!   [ drivall*dt  : Total inflow to the grid from surrounding grids]
!   [(stoall - sto2all + drunall*dt) / (135.3E12) : Mean runoff to the
!                                                   sea
      INTEGER iofile, iy, im, id, ih, ndev
      REAL rorder, stoall, sto2all, dinall, doutall                     &
     &     , drunall, drivall, dt

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      DATA rorder/1.0E12/
!
      IF (lhook) CALL dr_hook('WRTWBLOG',zhook_in,zhook_handle)
      IF (printstatus >= prstatus_diag) THEN
      WRITE(iofile, '(2A)')' River Routing: ',                          &
           'Error = Boxinflow - (sum of river inflow + box runoff)'

      WRITE (iofile, '(3A)') '         YY/MM/DD/HH ', '  S(t)  S(t+1)'  &
           , '     Din    Dout   Run  Rivflow  Error  Outflow(mm/y) '

      WRITE(iofile, '(7X,I2,A,I2,A,I2,A,I2,8(F7.1))')                      &
           iy, "/", im, "/", id, "/", ih,                               &
           stoall/rorder, sto2all/rorder, dinall/rorder, doutall/rorder &
           , (drunall*dt) / rorder, (drivall*dt) / rorder               &
           , (dinall - (drivall+drunall)*dt), (stoall - sto2all +       &
                drunall*dt) / (135.3E12) * 365.0 * REAL(ndev)

      WRITE (iofile, '(2A,E9.2,A,E9.2,A)')                              &
          "10^12 kg: ",  " Total outflow to sea or inland points is ",  &
     (stoall - sto2all + drunall*dt)* REAL(ndev)/86400, ' Kg/s or ',    &
     (stoall - sto2all + drunall*dt) / (135.3E12)* REAL(ndev), ' mm/day'

      WRITE(iofile, '(A)')' River routing: Water balance monitoring'
      WRITE(iofile, '(2A)')                                             &
      ' (Tot water in - tot water out - diff. water',                   &
      ' storage) summed for all gridboxes'

      WRITE (iofile, '(3A)') '          YY/MM/DD/HH ',                  &
       '  S(t)  S(t+1)', '     Din    Dout   Din-Dout-(S(t+1)- S(t))Kg'

      WRITE(iofile, '(A,I2,A,I2,A,I2,A,I2,5(1X,F7.1))')      &
           "10^12 kg: ", iy, "/", im, "/", id, "/", ih,                 &
           stoall/rorder, sto2all/rorder, dinall/rorder, doutall/rorder &
           , (dinall - doutall -(sto2all - stoall))
      END IF
      IF (lhook) CALL dr_hook('WRTWBLOG',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE wrtWBlog

END MODULE wrtwblog_mod
