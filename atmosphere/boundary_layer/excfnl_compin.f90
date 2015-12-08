! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE EXCFNL_COMPIN------------------------------------------
!
!  Purpose: Compute condition for inner interation and
!           compress index array
!
!  Code Description:
!    Language: FORTRAN 77 + common extensions.
!    This code is written to UMDP3 v6 programming standards.
!
!  Documentation: UMDP No.24
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!--------------------------------------------------------------------
SUBROUTINE excfnl_compin (                                              &
! IN fields
 up, wb_ratio, dec_thres, switch,                                       &
! INOUT fields
 c_len_i, ind_todo_i, todo_inner                                        &
 )

  USE atm_fields_bounds_mod, ONLY: pdims
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

! Intent IN:

  INTEGER, INTENT(IN) :: up(pdims%i_end*pdims%j_end)
  REAL,    INTENT(IN) ::                                                &
   wb_ratio(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                    ! WBN_INT/WBP_INT
   dec_thres(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! Local decoupling threshold

  INTEGER, INTENT(IN)           :: switch ! =1 for KSURF, 2 for KTOP

! Intent INOUT:

  INTEGER, INTENT(INOUT)        :: c_len_i
  INTEGER, INTENT(INOUT) :: ind_todo_i(pdims%i_end*pdims%j_end)
  LOGICAL, INTENT(INOUT) :: todo_inner(pdims%i_end*pdims%j_end)

! local variables
  INTEGER                       :: n,m
  INTEGER                       :: l, i1, j1

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('EXCFNL_COMPIN',zhook_in,zhook_handle)

!     Check for active elements
  SELECT CASE (switch)
  CASE(1)
            ! For top of KSURF
    DO n = 1, c_len_i
      l = ind_todo_i(n)
      j1=(l-1)/pdims%i_end+1
      i1=l-(j1-1)*pdims%i_end

      todo_inner(n) = (                                                 &
     (up(l) == 1 .AND. wb_ratio(i1,j1) <  dec_thres(i1,j1)) .OR.        &
                ! keep working up while wb_ratio lt thres
     (up(l) == 0 .AND. wb_ratio(i1,j1) >  dec_thres(i1,j1)) )
                ! keep working down while wb_ratio gt thres

    END DO

  CASE(2)
            ! For base of KTOP
    DO n = 1, c_len_i
      l = ind_todo_i(n)
      j1=(l-1)/pdims%i_end+1
      i1=l-(j1-1)*pdims%i_end

      todo_inner(n) = (                                                 &
     (up(l) == 1 .AND. wb_ratio(i1,j1) >  dec_thres(i1,j1)) .OR.        &
                ! keep working up while wb_ratio gt thres
     (up(l) == 0 .AND. wb_ratio(i1,j1) <  dec_thres(i1,j1)) )
                ! keep working down while wb_ratio lt thres

    END DO

  END SELECT

!     Compress

  m = 0
!CDIR Nodep
  DO n = 1, c_len_i
    IF(todo_inner(n))  THEN
      m = m+1
      todo_inner(m) = todo_inner(n)
      ind_todo_i(m) = ind_todo_i(n)
    END IF
  END DO
  c_len_i = m

  IF (lhook) CALL dr_hook('EXCFNL_COMPIN',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE excfnl_compin
