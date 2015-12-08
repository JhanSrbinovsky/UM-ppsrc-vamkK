! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! deck: FIELDADD
!
! contains routines: FieldAdd, etc. ScalarAdd etc.
!
! Purpose: Flux processing routines.
!          Simple arithmetic operations on pp fields
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs


      SUBROUTINE ScalarMult                                             &
     &           (nx, ny, rmdi,                                         &
     &            scalar, in_field,                                     &
     &            out_field,                                            &
     &            icode, cmessage)


!     ScalarMult: subroutine to multiply an array
!     ----------- by a scalar.
!                 Takes account of missing data.

!     Programmer: S J Foreman
!     -----------

!     Method:
!     -------  loop over all elements, multiplying unless either
!              is missing, in which case result is missing.
!              Scalar*array.

      IMPLICIT NONE


!     Input:
!     ------

      INTEGER                                                           &
     &       nx                                                         &
                        ! IN number of columns of array
     &      ,ny         ! IN number of rows in array

       REAL                                                             &
     &       rmdi                                                       &
                        ! IN value of REAL missing data indicator

     &      ,scalar                                                     &
                        ! IN scalar
     &      ,in_field(nx,ny)   ! IN array of input values

!     Output
!     ------

      REAL                                                              &
     &       out_field(nx,ny)  ! OUT results of multiplication

      INTEGER                                                           &
     &       icode      ! OUT Completion code

      CHARACTER(LEN=*)                                                    &
     &       cmessage   ! OUT Error message


!    Local variables
!    ---------------

      INTEGER                                                           &
     &        ix                                                        &
                        ! Loop counter over columns
     &       ,iy        ! Loop counter over rows


! ------------------------------------------------------------------

      icode = 0
      cmessage='ScalarMult: multiplication successful.'

      DO iy = 1, ny     ! iy: loop over rows
        DO ix = 1, nx  ! ix: loop over columns
          IF (    scalar  /=  rmdi                                      &
     &        .AND. in_field(ix,iy)  /=  rmdi )                         &
     &     THEN

             out_field(ix,iy) = scalar * in_field(ix,iy)

           ELSE      ! set missing data if either input is missing

             out_field(ix,iy) = rmdi

           END IF

         END DO  ! ix: loop over columns
      END DO     ! iy: loop over columns

      RETURN
      END SUBROUTINE ScalarMult
!----------------------------------------------------------------------
