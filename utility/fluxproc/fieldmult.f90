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

      SUBROUTINE FieldMult                                              &
     &           (nx, ny, rmdi,                                         &
     &            in_field1, in_field2,                                 &
     &            out_field,                                            &
     &            icode, cmessage)


!     FieldMult: subroutine to multiply two arrays.
!     ---------- Takes account of missing data.

!     Programmer: S J Foreman
!     -----------

!     Method:
!     -------  loop over all elements, adding unless either is
!              missing, in which case result is missing.
!              in_field1*in_field2.

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

     &      ,in_field1(nx,ny)                                           &
                               ! IN first array of input values
     &      ,in_field2(nx,ny)  ! IN second array of input values

!     Output
!     ------

      REAL                                                              &
     &       out_field(nx,ny)  ! OUT results of multipication

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

      icode =0
      cmessage="FieldMult: multiplication successful."

      DO iy = 1, ny     ! iy: loop over rows
        DO ix = 1, nx  ! ix: loop over columns

          IF (        in_field1(ix,iy)  /=  rmdi                        &
     &          .AND. in_field2(ix,iy)  /=  rmdi )                      &

     &      THEN

              out_field(ix,iy) = in_field1(ix,iy) * in_field2(ix,iy)

            ELSE

              out_field(ix,iy) = rmdi

            END IF

         END DO  ! ix: loop over columns
      END DO     ! iy: loop over columns

      RETURN
      END SUBROUTINE FieldMult

!----------------------------------------------------------------------
