! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE PRINT_INTE--------------------------------------------
!  
!   Purpose: Prints out a integer array to unit 7, formatting as four
!            numbers across a page.
!  
!  
!    Documentation: None
!  
!    -----------------------------------------------------------------
!    Arguments:-------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Small execs
      SUBROUTINE PRINT_INTE(A,N_POINTS,N_FIELD,K,NFTOUT)

      IMPLICIT NONE

      INTEGER                                                           &
     & N_POINTS                                                         &
                      !IN No of values to be printed
     &,N_FIELD                                                          &
                      !IN 1st dimension of array A
     &,K                                                                &
                      !IN Element in 2nd dimension of array A
     &,NFTOUT         !IN Output file unit number

      INTEGER                                                           &
     & A(N_FIELD)     !IN Array to be printed out

!*----------------------------------------------------------------------
!    Local variables:---------------------------------------------------
      INTEGER                                                           &
     & I,J            ! Loop indices

!*----------------------------------------------------------------------

!  1. Print out data modulo 4
        DO I=1,N_POINTS-3,4
          WRITE(NFTOUT,'(1x,4(I5,'':'',I12))')                          &
     &    I,A(I+(K-1)*N_FIELD),                                         &
     &    I+1,A(I+1+(K-1)*N_FIELD),                                     &
     &    I+2,A(I+2+(K-1)*N_FIELD),                                     &
     &    I+3,A(I+3+(K-1)*N_FIELD)
        ENDDO

!  2. Print out remainder of data
        IF(I <= N_POINTS)THEN
          WRITE(NFTOUT,'(1x)',ADVANCE='NO')
          DO J=I,N_POINTS
            WRITE(NFTOUT,'(I5,'':'',I12)',ADVANCE='NO')                 &

     &      J,A(J+(K-1)*N_FIELD)
          ENDDO
          WRITE(NFTOUT,'(/)')
        ENDIF

        RETURN
        END SUBROUTINE PRINT_INTE
