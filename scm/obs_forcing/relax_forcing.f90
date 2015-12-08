! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
! Routine to relax forcing profiles in the SCM.

SUBROUTINE relax_forcing                                                      &
  ! In
  ( row_length, rows, nlevs, var_n, var_init, var_bg, p, plev, tau, rlx       &

  ! InOut
  , var_inc )


  USE scm_utils, ONLY:                                                        &
    scm_timestep, old_rlx, zhook_in, zhook_out, jprb, lhook, dr_hook          &
  , rlx_none, rlx_init, rlx_bgrd, rlx_inst_init, rlx_inst_bgrd
    
  IMPLICIT NONE

  !
  ! Description:
  !   Returns increments profile (var_inc) which will relax specified SCM
  !   profile (var_n) to background profile (var_bg) or initial profile
  !   (var_init) on pressure levels where p < plev.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Single Column Model
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.3 programming standards.

  !---------------------------------------------------------------------------
  ! Declarations:
  !---------------------------------------------------------------------------

  ! Arguments
  !------------------------------------------------

  INTEGER, INTENT(In) :: &
    row_length           &
  , rows                 &
  , nlevs

  INTEGER, INTENT(In) :: &
    rlx                   ! Option to determine relaxation method

  REAL, INTENT(In) ::    &
    plev                 &! Pressure threshold, above which, apply relaxation
  , tau                   ! Timescale for relaxation (s)

  REAL, INTENT(In) ::                &
    var_n    (row_length,rows,nlevs) &! Profile (start of timestep)
  , var_init (row_length,rows,nlevs) &! Profile (Initial)
  , var_bg   (row_length,rows,nlevs) &! Profile (Background)
  , p        (row_length,rows,nlevs)  ! Pressure


  REAL, INTENT(InOut) ::             &
    var_inc  (row_length,rows,nlevs)  ! Current increment



  ! Local Arrays
  !------------------------------------------------
  INTEGER :: i,j,k

  REAL ::                            &
    rlx_fact (row_length,rows,nlevs) &
  , var_ref  (row_length,rows,nlevs)  ! Reference profile to relax to

  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('RELAX_FORCING',zhook_in,zhook_handle)

  SELECT CASE (rlx)

  CASE (rlx_none) ! no relaxation

  CASE (rlx_init, rlx_bgrd)

    IF (rlx == rlx_init) THEN
      var_ref(:,:,:) = var_init(:,:,:)
    ELSE
      var_ref(:,:,:) = var_bg(:,:,:)
    END IF

    DO k=1, nlevs
      DO j=1, rows
        DO i=1, row_length
          IF (p(i,j,k) < plev) THEN

            IF (old_rlx) THEN
              ! Pre 7.7
              var_inc(i,j,k) = var_inc(i,j,k)                           &
                             + ( var_ref(i,j,k) - var_n(i,j,k) )        &
                             * scm_timestep/tau
            ELSE
              ! Will differ because orignal relaxed wrt to var at
              ! start of timestep. Should really relax wrt to
              ! var at start of timestep + inc up to this point

              var_inc(i,j,k) = var_inc(i,j,k)                           &
                             + (scm_timestep/tau)                       &
                             * (  var_ref(i,j,k)                        &
                                - (var_n(i,j,k) + var_inc(i,j,k)) )

            END IF ! old_rlx

          END IF ! p < plev
        END DO
      END DO
    END DO

  CASE (rlx_inst_init, rlx_inst_bgrd)

    IF (rlx == rlx_inst_init) THEN
      var_ref(:,:,:) = var_init(:,:,:)
    ELSE
      var_ref(:,:,:) = var_bg(:,:,:)
    END IF

    DO k=1, nlevs
      DO j=1, rows
        DO i=1, row_length
          IF (p(i,j,k) < plev) THEN

            IF (old_rlx) THEN
              ! Pre 7.7
              var_inc(i,j,k) = var_inc(i,j,k)                           &
                             + ( var_ref(i,j,k) - var_n(i,j,k) )
            ELSE

              ! Will differ because orignal relaxed wrt to var at
              ! start of timestep. Should really relax wrt to
              ! var at start of timestep + inc up to this point
              var_inc(i,j,k) =  var_ref(i,j,k) - var_n(i,j,k)

            END IF ! old_rlx

          END IF ! p < plev
        END DO
      END DO
    END DO

  END SELECT

  IF (lhook) CALL dr_hook('RELAX_FORCING',zhook_out,zhook_handle)


  RETURN
END SUBROUTINE relax_forcing
