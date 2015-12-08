! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine COPYDIAG_3D -------------------------------------------
!  
!   Purpose : To copy a diagnostic field from secondary space to the
!             main data array for stash processing, and to extend the
!             data to a full horizontal field. Input data of multilevel
!             fields is assumed to be on all model levels. Output data
!             is on the levels required.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

      SUBROUTINE copydiag_3d(                                           &
     &     diagout,diagin,                                              &
     &     row_length,rows,levels,                                      &
     &     offx_out, offy_out,                                          &
     &     offx_in, offy_in,                                            &
     &     at_extremity,                                                &
     &     stlist,len_stlist,stash_levels,                              &
     &     len_stashlevels,                                             &
     &     im,is,ie,                                                    &
     &     icode,cmessage)

      USE um_input_control_mod, ONLY: model_domain
      USE domain_params, ONLY: mt_global

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook

      IMPLICIT NONE

      INTEGER                                                           &
     &       row_length,                                                &
                          ! Number of points in a row
     &       rows,                                                      &
                          ! Number of rows
     &       offx_out, offy_out,                                        &
                                 ! offset dimensions of diagout
     &       offx_in, offy_in,                                          &
                                 ! offset dimensions of diagin
     &       levels,                                                    &
                          ! Number of levels in input data
     &       len_stlist,                                                &
                          !
     &       stlist(len_stlist),                                        &
                                 ! Stash list
     &       len_stashlevels,                                           &
                              !
     &       stash_levels(len_stashlevels,*),                           &
                                              ! Stash levels list.
     &       im,is,ie,                                                  &
                          ! Model, section, item
     &       icode        ! Return code =0 Normal exit
!                                       >1 Error message
      LOGICAL    at_extremity(4)
!! logicals indicating if a processor
                             ! is at the edge of the LPG
! argppx arguments:

      CHARACTER(LEN=80) cmessage

      REAL                                                              &
     &   diagin(1-offx_in:row_length+offx_in, 1-offy_in:rows+offy_in    &
     &          ,levels)                                                &
                                ! Output data
     &,  diagout(1-offx_out:row_length+offx_out,1-offy_out:rows+offy_out&
     &          ,*)             ! Input data

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
     &   i,j,                                                           &
                               !
     &   k,                                                             &
                               !
     &   kout                 !
!     &,  gr                   ! grid type of grid
!     &,  fld_type             ! is it a p field or a u field?
!     &,  info                 ! GCOM return code

! Functions called:
!      Integer
!     &   exppxi,get_fld_type


!      Real
!     &   copy_value_start(levels)  ! value to copy into start of array
!     &,  copy_value_end(levels)    ! value to copy into end of array


      LOGICAL                                                           &
     &   list(levels) !

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('COPYDIAG_3D',zhook_in,zhook_handle)
      icode = 0
      cmessage = ""

! Find out the gridtype of the field
!      gr = exppxi(im,is,ie,ppx_grid_type,
!     &            icode,cmessage)

!      If (icode  >   0) goto 9999

! and use this to find the field type (p field or u field)

!      fld_type=get_fld_type(gr)


! Find the values to copy into the start and end of the arrays

!      Do k=1,levels
!        copy_value_start(k)=diagin(start_point,k)
!        copy_value_end(k)=diagin(end_point,k)
!      Enddo

! If this is the Northern processor row - we must make sure we
! get a consistent value over the polar row - so we take
! the value from PE 0 and use that on all processors
!      If (attop) then
!        Call gcg_rbcast(700,levels,first_comp_pe,gc_proc_row_group,
!     &                  info,copy_value_start)
!      Endif

! If this is the Southern processor row - we must make sure we
! get a consistent value over the polar row - so we take
! the value from PE 0 and use that on all processors
!      If (atbase) then
!        Call gcg_rbcast(701,levels,last_comp_pe,gc_proc_row_group,
!     &                  info,copy_value_end)
!      Endif

! DEPENDS ON: set_levels_list
      CALL set_levels_list(levels,len_stlist,stlist,list,stash_levels,  &
     &      len_stashlevels,icode,cmessage)
      IF(icode >  0) GOTO 9999

!L Move data from DIAGIN to DIAGOUT at levels requested

      kout = 0
      DO k = 1, levels
         IF (list(k)) THEN
            kout = kout + 1

!           Copy fields
            DO j = 1, rows
               DO i = 1, row_length
                  diagout(i,j,kout) = diagin(i,j,k)
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
                     diagout(i,j,kout) = diagout(i,rows,kout)
                  END DO
               END DO
            END IF

            IF (at_extremity(PSouth)) THEN
               DO i = 1-offx_out, row_length+offx_out
                  DO j = 1-offy_out, 1
                     diagout(i,j,kout) = diagout(i,1,kout)
                  END DO
               END DO
            END IF

!l copy diagnostic information to e and w boundaries

            IF (at_extremity(PEast)) THEN
               DO i = 1-offx_out, 1
                  DO j = 1-offy_out, rows+offy_out
                     diagout(i,j,kout) = diagout(1,j,kout)
                  END DO
               END DO
            END IF


            IF (at_extremity(PWest)) THEN
               DO i = row_length, row_length+offx_out
                  DO j = 1-offy_out, rows+offy_out
                     diagout(i,j,kout) = diagout(row_length,j,kout)
                  END DO
               END DO
            END IF

         ENDIF      ! check no halos
         END IF     ! .not. GLOBAL

         END IF
      END DO

 9999 CONTINUE
      IF (lhook) CALL dr_hook('COPYDIAG_3D',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE copydiag_3d
