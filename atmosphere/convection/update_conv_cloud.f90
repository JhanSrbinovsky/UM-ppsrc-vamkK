! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------

! Update convective cloud after each substep

SUBROUTINE  update_conv_cloud(                                               &
             rows, row_length, model_levels, wet_levels,                     &
             n_cca_lev, ndecay_steps, n_conv_calls,                          &
             loop_number,                                                    &
             it_ccb, it_cct, it_lcbase, it_lctop, it_ccb0, it_cct0, it_lcbase0,&
             one_over_conv_calls, timestep, decay_time,                      &
             p_layer_boundaries, cld_life_3d,                                &
             it_cca_2d, it_lcca, it_cca, it_ccw, it_cclwp,                   &
             it_cca0_2d, it_cca0,it_ccw0,                                    &
             ! in/out 
             ccb, cct, lcbase,lctop, ccb0, cct0, lcbase0,                    &
             ccb0_local, cct0_local,lcbase0_local,                           &
             lcca,  cca_2d,   cclwp,  cca,  ccw,                             &
             cca0_2d,  cclwp0, cca0, ccw0,                                   &
             it_cclwp0, cclwp0_local, cca0_local, ccw0_local)

USE cv_run_mod,  ONLY:                                                       &
    l_fix_udfactor,  ud_factor,                                              &
    rad_cloud_decay_opt,     cld_life_opt,          cca_min,                 &
    fixed_cld_life, l_dcpl_cld4pc2, l_ccrad

USE cv_param_mod,  ONLY:                                                     &
    cld_life_constant,       cld_life_func_hgt,     rad_decay_off,           &
    rad_decay_full_timestep, rad_decay_conv_substep

USE earth_constants_mod, ONLY: g

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Updates convective cloud diagnostics and prognostics after call to
!  glue each substep of  convection.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection

! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.

! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) ::  &
  row_length            & ! row length
 ,rows                  & ! rows
 ,model_levels          & ! model levels
 ,wet_levels            & ! wet levels
 ,n_cca_lev             & ! Convective cloud levels
 ,ndecay_steps          & ! CCRad variable
 ,n_conv_calls          & ! Number of convection calls
 ,loop_number             ! Counter of substeps

! Subscript 0 on variable name indicates it is the section 0 value
! otherwise section 5 value.

INTEGER, INTENT(IN) ::         &
  it_ccb(row_length,rows)      &! Cnv Cld Base of highest layer on gridpoint
 ,it_cct(row_length,rows)      &! Cnv Cld Top  of highest layer on gridpoint
 ,it_lcbase(row_length,rows)   &! Cnv Cld Base of lowest  layer on gridpoint
 ,it_lctop(row_length,rows)    &! Cnv Cld top of lowest  layer on gridpoint
 ,it_ccb0(row_length,rows)     &! Cnv Cld Base of highest layer on gridpoint
 ,it_cct0(row_length,rows)     &! Cnv Cld Top  of highest layer on gridpoint
 ,it_lcbase0(row_length,rows)   ! Cnv Cld Base of lowest  layer on gridpoint

REAL, INTENT(IN) ::      &
  one_over_conv_calls    & ! 1/n_conv_calls
 ,timestep               & ! Full model timestep (s)
 ,decay_time               ! decay time for conv cloud

REAL, INTENT(IN) ::                                                 &
  p_layer_boundaries(row_length, rows, 0:model_levels)              &
                         ! pressure at layer boundaries. Same as p except at
                         ! bottom level = pstar, and at top = 0.
 ,cld_life_3d(row_length, rows, n_cca_lev) &! CCRad - cloud life time array 
 ,it_cca_2d(row_length,rows)               &! Cnv.Cld Amount (2d) with no anvil
 ,it_cca(row_length,rows,n_cca_lev)        &! Cnv.Cld Amount (0-1)
 ,it_lcca(row_length,rows)                 &! Cnv.Cld Amount (0-1) lowest layer
 ,it_ccw(row_length,rows,wet_levels)       &! Cnv.Cld Water (kg/kg)
 ,it_cca0_2d(row_length,rows)              &! Cnv.Cld Amount (2d) with no anvil
 ,it_cca0(row_length,rows,n_cca_lev)       &! Cnv.Cld Amount (0-1)
 ,it_ccw0(row_length,rows,wet_levels)      &! Cnv.Cld Water (kg/kg) 
 ,it_cclwp(row_length,rows)                 ! Cloud Condensed water path(kg/m^2)

INTEGER, INTENT(INOUT) ::       &
  ccb(row_length,rows)          &! Cnv Cld Base of highest layer on gridpoint
 ,cct(row_length,rows)          &! Cnv Cld Top  of highest layer on gridpoint
 ,lcbase (row_length,rows)      &! Cnv Cld base of lowest  layer on gridpoint
 ,lctop(row_length,rows)        &! Cnv Cld Top  of lowest  layer on  gridpoint
 ,ccb0(row_length,rows)         &! Cnv Cld Base of highest layer on gridpoint
 ,cct0(row_length,rows)         &! Cnv Cld Top  of highest layer on gridpoint
 ,lcbase0 (row_length,rows)     &! Cnv Cld Top  of lowest  layer on gridpoint
 ,ccb0_local(row_length,rows)   &! Holds mean profile until decayed 
 ,cct0_local(row_length,rows)   &! Holds mean profile until decayed 
 ,lcbase0_local(row_length,rows) ! Holds mean profile until decayed 
 
REAL, INTENT(INOUT) ::                   &
  it_cclwp0 (row_length,rows)            &! Cloud Condensed water path(kg/m^2)
 ,lcca(row_length,rows)                  &! Conv. Cloud Amount of low cloud
 ,cca_2d(row_length,rows)                &! Cnv.Cld Amount (2d)
 ,cclwp(row_length,rows)                 &! Cloud Condensed water path(kg/m^2)
                                          ! with no anvil  
 ,cca(row_length,rows,n_cca_lev)         &! Cnv.Cld Amount (0-1)
 ,ccw(row_length,rows,wet_levels)        &! Cnv.Cld Water (kg/kg) 
 ,cca0_2d(row_length,rows)               &! Cnv.Cld Amount (2d)
 ,cclwp0(row_length,rows)                &! Cloud Condensed water path(kg/m^2)
                                          ! with no anvil  
 ,cca0(row_length,rows,n_cca_lev)        &! Cnv.Cld Amount (0-1)
 ,ccw0(row_length,rows,wet_levels)       &! Cnv.Cld Water (kg/kg) 
 ,cca0_local(row_length,rows,n_cca_lev)  &
 ,ccw0_local(row_length,rows,wet_levels) &
 ,cclwp0_local(row_length,rows)


! Local variables
INTEGER  ::         &
  i,j,k                ! loop counters

LOGICAL :: decay_now

 
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook('UPDATE_CONV_CLOUD',zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
!  Application of udfactor
!-----------------------------------------------------------------------------
    IF (l_fix_udfactor) THEN
      ! Use ud_factor of 1 in convection to return CCW and CCLWP
      ! unscaled.  Then the true ud_factor can be applied here to
      ! the whole profile, rather than just to precipitating levels

      ! Updraught factor was used to correct for radiative impacts.
      ! As such it's use should only be reflected in diagnostics
      ! Seen by radiation.

      ! Note: This could have be done simpler by either applying the
      !       factor correctly lower down in the code or only
      !       applying the factor at this level.
      !       It is only done this way because some models may wish
      !       to retain the incorrect version.
      !       If all models are running with either
      !       l_fix_udfactor=.TRUE., or l_ccrad=.TRUE. then this
      !       switch could be removed and the code simplified
      !       by correcting the code lower down.
      !       This switch should NOT be used with l_ccrad

      DO j=1, rows
        DO i=1, row_length
          it_cclwp0(i,j) = ud_factor * it_cclwp0(i,j)
        END DO
      END DO
    END IF       ! l_fix_udfactor

!-----------------------------------------------------------------------------
! Update cloud properties
!-----------------------------------------------------------------------------


    ! NOTE: Section 5 cloud properties were initalised in Atmos_physics2
    !       before physics Loop
    DO j=1, rows
      DO i=1, row_length
        ! Highest convective layer properties
        !------------------------------------
        ! MAX cct across total number of calls to convection
        cct(i,j) = MAX(cct(i,j), it_cct(i,j))

        ! MIN ccb across total number of calls to convection
        ! excluding ccb=0
        IF (ccb(i,j) > 0 .AND. it_ccb(i,j) > 0) THEN
          ccb(i,j) = MIN(ccb(i,j), it_ccb(i,j))
        ELSE
          ccb(i,j) = MAX(ccb(i,j), it_ccb(i,j))
        END IF

        ! Lowest convective layer properties
        !------------------------------------
        ! MAX lctop across total number of calls to convection
        lctop(i,j) = MAX(lctop(i,j), it_lctop(i,j))

        ! MIN lcbase across total number of calls to convection
        ! excluding lcbase=0
        IF (lcbase(i,j) > 0 .AND. it_lcbase(i,j) > 0) THEN
          lcbase(i,j) = MIN(lcbase(i,j), it_lcbase(i,j))
        ELSE
          lcbase(i,j) = MAX(lcbase(i,j), it_lcbase(i,j))
        END IF

        lcca(i,j)   = lcca(i,j)   + one_over_conv_calls*it_lcca(i,j)
        cca_2d(i,j) = cca_2d(i,j) + one_over_conv_calls*it_cca_2d(i,j)
        cclwp(i,j)  = cclwp(i,j)  + one_over_conv_calls*it_cclwp(i,j)

      END DO
    END DO

    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          ccw(i,j,k) = ccw(i,j,k) + one_over_conv_calls*it_ccw(i,j,k)
        END DO
      END DO
    END DO

    DO k=1, n_cca_lev
      DO j=1, rows
        DO i=1, row_length
          cca(i,j,k) = cca(i,j,k) + one_over_conv_calls*it_cca(i,j,k)
        END DO
      END DO
    END DO

    IF (l_ccrad) THEN
      DO k=1, n_cca_lev
        DO j=1, rows
          DO i=1, row_length
            cca(i,j,k) = MIN(cca(i,j,k), 1.0)
          END DO
        END DO
      END DO
    END IF

!=============================================================================
! Radiative Cloud Decay: Section 0 diagnostics
!=============================================================================

! NOTE: CCA0, CCW0 profiles are the resultant profile seen by radiation,
!       they may contain elements of cloud from previous
!       timesteps/conv_timesteps.
!=============================================================================

    SELECT CASE(rad_cloud_decay_opt)
      ! Cloud properties for:
      ! Section 0 = Section 5 + radiative scalings
      !-------------------------------------------
      CASE(rad_decay_off)
        ! Perform same operation as that on seciton 0 diagnostics,
        ! NOTE: Cannot just copy section 5 diagnostics as section 0
        !       diagnostics may not contain all convective cloud types

        ! Radiative decay is off, so these section 0 variables
        ! will have been initialised at the start of atmos_physics2

        DO j=1, rows
          DO i=1, row_length

            ! Highest convective layer properties
            !------------------------------------
            ! MAX cct0 across total number of calls to convection
            cct0(i,j) = MAX(cct0(i,j), it_cct0(i,j))

            ! MIN ccb0 across total number of calls to convection
            ! excluding ccb0=0
            IF (ccb0(i,j) > 0 .AND. it_ccb0(i,j) > 0) THEN
              ccb0(i,j) = MIN(ccb0(i,j), it_ccb0(i,j))
            ELSE
              ccb0(i,j) = MAX(ccb0(i,j), it_ccb0(i,j))
            END IF

            ! Lowest convective layer properties
            !------------------------------------
            ! MIN lcbase0 across total number of calls to convection
            ! excluding lcbase0=0
            IF (lcbase0(i,j) > 0 .AND. it_lcbase0(i,j) > 0) THEN
              lcbase0(i,j) = MIN(lcbase0(i,j), it_lcbase0(i,j))
            ELSE
              lcbase0(i,j) = MAX(lcbase0(i,j), it_lcbase0(i,j))
            END IF

            cclwp0 (i,j) = cclwp0(i,j)                                  &
                         + one_over_conv_calls * it_cclwp0(i,j)
          END DO
        END DO

        DO k=1, wet_levels
          DO j=1, rows
            DO i=1, row_length
              ccw0(i,j,k) = ccw0(i,j,k)                                 &
                          + one_over_conv_calls * it_ccw0(i,j,k)
            END DO
          END DO
        END DO

        DO k=1, n_cca_lev
          DO j=1, rows
            DO i=1, row_length
              cca0(i,j,k) = cca0(i,j,k)                                 &
                          + one_over_conv_calls * it_cca0(i,j,k)
            END DO
          END DO
        END DO

        IF (l_ccrad) THEN
          DO k=1, n_cca_lev
            DO j=1, rows
              DO i=1, row_length
                cca0(i,j,k) = MIN(cca0(i,j,k), 1.0)
              END DO
            END DO
          END DO
        END IF

     !==================================================================
      CASE( rad_decay_conv_substep                                      &
          , rad_decay_full_timestep )
     !==================================================================

        decay_now = .FALSE.

        ! Note: _local arrays were initialised
        !       at start of convection scheme

        IF (ndecay_steps == 1) THEN
          ! The Section 0 cloud properties will be compared against
          ! those output on each substep.
          DO j=1, rows
            DO i=1, row_length
              cct0_local(i,j)   = it_cct0(i,j)
              ccb0_local(i,j)   = it_ccb0(i,j)
              cclwp0_local(i,j) = it_cclwp0(i,j)
            END DO
          END DO

          DO k=1, n_cca_lev
            DO j=1, rows
              DO i=1, row_length
                cca0_local(i,j,k) = it_cca0(i,j,k)
              END DO
            END DO
          END DO

          DO k=1, wet_levels
            DO j=1, rows
              DO i=1, row_length
                ccw0_local(i,j,k) = it_ccw0(i,j,k)
              END DO
            END DO
          END DO

          decay_now = .TRUE.

        ELSE
          ! ndecay_steps /= 1

          ! The Section 0 cloud properties will be compared against
          ! those output on full convection timestep.
          ! ccw/cca are meaned across convection substeps.
          ! ccb, cct are min/max across convection substeps.

          DO j=1, rows
            DO i=1, row_length
              cct0_local(i,j) = MAX(cct0(i,j), it_cct0(i,j))

              IF (ccb0(i,j) > 0 .AND. it_ccb0(i,j) > 0) THEN
                ccb0_local(i,j) = MIN(ccb0(i,j), it_ccb0(i,j))
              ELSE
                ccb0_local(i,j) = MAX(ccb0(i,j), it_ccb0(i,j))
              END IF

              cclwp0_local(i,j)   = cclwp0_local(i,j)                   &
                                  + one_over_conv_calls * it_cclwp0(i,j)
            END DO
          END DO

          DO k=1, n_cca_lev
            DO j=1, rows
              DO i=1, row_length
                cca0_local(i,j,k) = cca0_local(i,j,k)                   &
                                  + one_over_conv_calls * it_cca0(i,j,k)
              END DO
            END DO
          END DO

          DO k=1, wet_levels
            DO j=1, rows
              DO i=1, row_length
                ccw0_local(i,j,k) = ccw0_local(i,j,k)                   &
                                  + one_over_conv_calls * it_ccw0(i,j,k)
              END DO
            END DO
          END DO

          ! Set decay flag if last convection substep
          IF (loop_number == n_conv_calls) THEN
            decay_now = .TRUE.
          END IF

        END IF ! ndecay_steps

        !---------------------------------------------------------------
        ! Decay (Section 0) Convective Cloud amounts for radiation
        !---------------------------------------------------------------
        IF (decay_now) THEN

          ! Re-initialise all parameters sent to radiation as they
          ! should refer to the resultant cloud profiles after cloud
          ! decay. lcbase0, ccb0 and cct0 are all calculated
          ! from the resultant cca0 profile.
          DO j=1, rows
            DO i=1, row_length
              lcbase0 (i,j) = 0
              ccb0    (i,j) = 0
              cct0    (i,j) = 0
            END DO
          END DO


          ! Decay of CCA0 Profile
          !----------------------
          ! Decay current CCA0 profile and compare profile
          ! against mean convection over ndecay_steps.
          ! For a given model level, the higher CCA is retained in
          ! the CCA0 array for radiation. New cloud properties for
          ! radiation are calculated based on the resultant profile.

          DO k=1, n_cca_lev
            DO j=1, rows
              DO i=1, row_length

                ! Determine the higher value of CCA, i.e. from current
                ! call to convection OR decayed CCA0 from previous
                ! timesteps
                cca0(i,j,k) = MAX( cca0(i,j,k)                        &
                                     * (1.0-cld_life_3d(i,j,k))       &
                                 , cca0_local(i,j,k) )

                cca0(i,j,k) = MIN( cca0(i,j,k), 1.0 )

                ! Apply threshold resultant profile for minimum CCA
                IF (cca0(i,j,k) < cca_min) THEN
                  cca0(i,j,k) = 0.0
                END IF

              END DO
            END DO
          END DO


          IF (l_dcpl_cld4pc2) THEN

            ! Identify new ccb0, cct0 from decayed cca profile
            !-------------------------------------------------
            DO k=2, n_cca_lev
              DO j=1, rows
                DO i=1, row_length

                  ! Check for resultant cloud bases for radiation
                  IF (cca0(i,j,k) > 0.0 .AND. cca0(i,j,k-1) == 0.0) THEN
                    ccb0(i,j) = k
                    IF (lcbase0(i,j) == 0) THEN
                      lcbase0(i,j) = k
                    END IF
                  END IF

                  ! Check for resultant cloud tops for radiation
                  IF (cca0(i,j,k) == 0.0 .AND. cca0(i,j,k-1) > 0.0) THEN
                    cct0(i,j) = k-1
                  END IF

                END DO
              END DO
            END DO


            ! Decay of CCW0 Profile on levels which still contain
            ! cca > 0.0
            !-----------------------------------------------------------
            IF (l_ccrad) THEN
              ! Reset cclwp0, as it will based on decay ccw0 pprofile
              DO j=1, rows
                DO i=1, row_length
                  cclwp0(i,j) = 0.0
                END DO
              END DO

              ! Decay of ccw0
              DO k=1, wet_levels
                DO j=1, rows
                  DO i=1, row_length
                    IF (cca0(i,j,k) > 0.0) THEN
                      ccw0(i,j,k) = MAX( ccw0(i,j,k)                    &
                                           * (1.0-cld_life_3d(i,j,k))   &
                                       , ccw0_local(i,j,k) )

                      cclwp0(i,j) = cclwp0(i,j) + ccw0(i,j,k)           &
                                  * (  p_layer_boundaries(i,j,k-1)      &
                                     - p_layer_boundaries(i,j,k) )/g
                    ELSE
                      ccw0(i,j,k) = 0.0
                    END IF
                  END DO
                END DO
              END DO

            ELSE
              ! Original cclwp0 decay used since ccw0 is not held
              ! from one timestep to another unless CCRad is enabled
              DO j=1, rows
                DO i=1, row_length
                  IF (ccb0(i,j) == 0) THEN
                    cclwp0(i,j) = 0.0
                  ELSE
                    cclwp0(i,j) = MAX( cclwp0(i,j)                      &
                                * (1.0-timestep/fixed_cld_life)         &
                                , cclwp0_local(i,j))
                  END IF
                END DO
              END DO
            END IF

          ELSE   ! l_dcpl_cld4pc2
            !----------------------------------------------------------
            ! NOT-RECOMMENEDED FOR USE
            !----------------------------------------------------------
            ! This is to replicate the original code cloud decay code
            ! It is flawed however, and is only here in case the code
            ! needs to be reverted

            IF (rad_cloud_decay_opt == rad_decay_full_timestep) THEN

              DO k=2, n_cca_lev
                DO j=1, rows
                  DO i=1, row_length
                    IF ((ccb0(i,j)   == 0  ) .AND.                     &
                        (cca0(i,j,k) >  0.0)) THEN
                      ccb0(i,j) = k
                      IF (lcbase0(i,j) == 0) THEN
                        lcbase0(i,j) = k
                      END IF
                    END IF

                    IF (cca0(i,j,k) > 0.0) THEN
                      cct0(i,j) = k+1
                    END IF
                  END DO
                END DO
              END DO

              DO j=1, rows
                DO i=1, row_length
                  IF (ccb0(i,j) == 0) THEN
                    cclwp0(i,j) = 0.0
                  ELSE
                    cclwp0(i,j) = MAX( cclwp0(i,j)                      &
                                * (1.0-timestep/fixed_cld_life)         &
                                , cclwp0_local(i,j))
                  END IF
                END DO
              END DO

            ELSE  ! Decay will occurr on each convection substep

              DO k=2, n_cca_lev
                DO j=1, rows
                  DO i=1, row_length

                    IF (cca0(i,j,k) > 0.0 .AND.                         &
                        cca0(i,j,k-1) == 0.0) THEN
                      ccb0(i,j) = k
                      IF (lcbase0(i,j) == 0) THEN
                        lcbase0(i,j) = k
                      END IF
                    END IF

                    IF (cca0(i,j,k) == 0.0 .AND.                        &
                        cca0(i,j,k-1) > 0.0) THEN
                      cct0(i,j) = k-1
                    END IF
                  END DO
                END DO
              END DO

              DO j=1, rows
                DO i=1, row_length
                  cclwp0(i,j) = 0.0
                END DO
              END DO

              DO k=1, wet_levels
                DO j=1, rows
                  DO i=1, row_length
                    IF (cca0(i,j,k) > 0.0) THEN
                      ccw0(i,j,k) = MAX( ccw0(i,j,k)                    &
                                           * (1.0-cld_life_3d(i,j,k))   &
                                       , ccw0_local(i,j,k) )

                      cclwp0(i,j) = cclwp0(i,j) + ccw0(i,j,k)           &
                                  * (  p_layer_boundaries(i,j,k-1)      &
                                     - p_layer_boundaries(i,j,k) )/g
                    ELSE
                      ccw0(i,j,k) = 0.0
                    END IF
                  END DO
                END DO
              END DO
            END IF              ! Rad_cloud_decay_opt

          END IF ! l_dcpl_cld4pc2

        END IF ! decay_now

    END SELECT ! Rad_cloud_decay_opt


!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('UPDATE_CONV_CLOUD',zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
END SUBROUTINE
