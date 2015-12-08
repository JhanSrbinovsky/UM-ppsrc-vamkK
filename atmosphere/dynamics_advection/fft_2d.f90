! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine fft_2d
      SUBROUTINE FFT_2D(global_rows,global_row_length,spectra2          &
     &              ,spectra_im2,INIT,TRIGM,TRIGN,WORK,IFAIL)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Description:
!   Dummy subroutine.  To be replaced by a subroutine which calculates
!   2D FFTs.
!
! Method:
!   Passes the original field back without performing any calculations
!   on it.
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
!+ Dates should have leading zeroes in dd/mm/yy format
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Global variables (#include statements etc):

      Integer                                                           &
     & global_rows                                                      &
     &,global_row_length                                                &
     &,ifail

      Real                                                              &
     & SPECTRA2(global_row_length*global_rows)                          &
     &,SPECTRA_IM2(global_row_length*global_rows)                       &
     &,TRIGN(2*GLOBAL_ROW_LENGTH)                                       &
                                     !  HOLDS TRIGONOMETRIC TERMS
     &,TRIGM(2*GLOBAL_ROWS)                                             & 
                                     !  USED IN FFT'S
     &,WORK(2*GLOBAL_ROW_LENGTH*GLOBAL_ROWS)

      CHARACTER(LEN=1)  INIT

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!---------------------------------------------------------------------
! Section 1.
!---------------------------------------------------------------------

      IF (lhook) CALL dr_hook('FFT_2D',zhook_in,zhook_handle)
      write(6,*)'**************WARNING*****************'
      write(6,*)'******FFT_2D is a dummy routine*******'
      write(6,*)'*****original field is passed back****'
      write(6,*)'**************************************'

!!    END OF ROUTINE FFT_2D
      IF (lhook) CALL dr_hook('FFT_2D',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FFT_2D
