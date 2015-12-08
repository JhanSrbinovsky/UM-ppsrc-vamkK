! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE initial0_mod

USE setarea_mod, ONLY: setarea
USE setlen_mod, ONLY: setlen
USE setnext_mod, ONLY: setnext
IMPLICIT NONE

CONTAINS

!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing


      SUBROUTINE initial0(nx, ny, jmax, rmiss, igrcn, iseq              &
     , nseqmax, inextx, inexty, rlen, area, offset_nx, offset_ny)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!     Setting Initial Conditions
!
      INTEGER, INTENT(IN) ::                                            &
       nx                                                               &
        ! number of columns
     , ny                                                               &
        ! number of rows
     , jmax                                                             &
        ! number of rows ! redundant
     , offset_nx                                                        &
        ! column offset from latitude origin
     , offset_ny                                                        &    
        ! row offset from longitude origin  
     , igrcn(nx, ny)                                                    &
        ! river direction 
     , iseq(nx, ny)                                                   
        ! river sequence field 
                     

      INTEGER, INTENT(OUT) ::                                           &
     &  nseqmax                                                         &
       ! max sequence in river subdomain 
     , inextx(nx, ny)                                                   &
       ! next lat downstream points for a river grid point
     , inexty(nx, ny)                                                 
       ! next long downstream point for a river grid point

      REAL, INTENT(OUT) ::                                              &
     &  rlen(nx, ny)                                                      
        ! distance between river grid points

      REAL, INTENT(IN) ::                                               &
     &  rmiss
        ! flag for missing value

      REAL, INTENT(OUT) ::                                              &
     &  area(nx, ny)
        ! river grid point area [m^2]
      INTEGER i, j 
        ! loop counters

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      
      IF (lhook) CALL dr_hook('INITIAL0',zhook_in,zhook_handle)

! get the maximum river sequence in this river subdomain
      nseqmax = 0
      DO i = 1, nx
        DO j= 1, ny
          IF ( iseq(i,j) >  nseqmax) nseqmax = iseq(i,j)
        END DO
      END DO

! set downstream points

      CALL setnext(nx, ny, igrcn, inextx, inexty)


! set the distance between grids [m]

      CALL setlen (nx, ny, igrcn, inextx, inexty, rlen, jmax        &
     , rmiss, offset_nx, offset_ny)

! set area of each river grid point [m^2]

      CALL setarea(nx, ny, area, jmax, offset_ny)


      IF (lhook) CALL dr_hook('INITIAL0',zhook_out,zhook_handle)
      RETURN

    END SUBROUTINE initial0
END MODULE initial0_mod
