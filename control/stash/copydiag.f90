! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Purpose : To copy a single diagnostic field from secondary space to
!             the main data array for stash processing, and to extend
!             the data to a full horizontal field.
!  
!  
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

      SUBROUTINE copydiag(                                              &
     &     diagout,diagin,                                              &
     &     row_length,rows,                                             &
     &     offx_out, offy_out,                                          &
     &     offx_in, offy_in,                                            &
     &     at_extremity,                                                &
     &     im,is,ie,                                                    &
     &     icode,cmessage)

      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook

      IMPLICIT NONE

      INTEGER                                                           &
     &     row_length                                                   &
                                ! Number of points in a row
     &,    rows                                                         &
                                ! Number of rows
     &,    offx_out, offy_out                                           &
                                ! offset dimensions of diagout
     &,    offx_in, offy_in                                             &
                                ! offset dimensions of diagin
     &,    im,is,ie                                                     &
                                ! Model, section, item
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      LOGICAL                                                           &
     &     at_extremity(4)      ! Indicates if this processor is at
                                !  north, south east or west of the
                                !  processor grid
      CHARACTER(LEN=80)  cmessage

! ARGPPX arguments:


      REAL                                                              &
     & diagin(1-offx_in:row_length+offx_in, 1-offy_in:rows+offy_in)     &
                                ! Output field
     &,diagout(1-offx_out:row_length+offx_out,1-offy_out:rows+offy_out)
                                ! Input field


! Parameters
      INTEGER                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      PARAMETER (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

!     Local variables

      INTEGER                                                           &
     &   i,j  ! loop bound
!     &,  gr                     ! grid type of grid
!     &,  fld_type               ! is it a p field or a u field?
!     &,  info                   ! GCOM return code

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Functions called:
!      Integer
!     &   exppxi ,get_fld_type

      IF (lhook) CALL dr_hook('COPYDIAG',zhook_in,zhook_handle)
      icode = 0
      cmessage = ""

! Find out the gridtype of the field
!      gr = exppxi(im,is,ie,ppx_grid_type,
!     &            ICODE,CMESSAGE)

! and use this to find the field type (p field, u field or v field)

!      fld_type = get_fld_type(gr)


!     Copy fields
      DO j = 1, rows
         DO i = 1, row_length
            diagout(i,j) = diagin(i,j)
         END DO
      END DO

      IF (model_domain /= mt_global) THEN

! Marker only for possible extensions of code. Since diagnostics are
! so far defined as having no halos, extra code to populate halos is
! not used at present, so bypassed here on in normal case, when
! offx_out=offy_out=0 .

      IF(offx_out /= 0.OR.offy_out /= 0) THEN ! check no halos

      IF (at_extremity(PNorth)) THEN
         DO i = 1-offx_out, row_length+offx_out
            DO j = rows, rows+offy_out
               diagout(i,j) = diagout(i,rows)
            END DO
         END DO
      END IF

      IF (at_extremity(PSouth)) THEN
         DO i = 1-offx_out, row_length+offx_out
            DO j = 1-offy_out, 1
               diagout(i,j) = diagout(i,1)
            END DO
         END DO
      END IF

!l copy diagnostic information to e and w boundaries

      IF (at_extremity(PEast)) THEN
         DO i = 1-offx_out, 1
            DO j = 1-offy_out, rows+offy_out
               diagout(i,j) = diagout(1,j)
            END DO
         END DO
      END IF


      IF (at_extremity(PWest)) THEN
         DO i = row_length, row_length+offx_out
            DO j = 1-offy_out, rows+offy_out
               diagout(i,j) = diagout(row_length,j)
            END DO
         END DO
      END IF

      ENDIF      ! check no halos
      END IF     ! .not. GLOBAL

      IF (lhook) CALL dr_hook('COPYDIAG',zhook_out,zhook_handle)

      END SUBROUTINE copydiag
