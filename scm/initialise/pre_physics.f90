! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up variables before the call to the physics routines.
!
! Subroutine Interface:

SUBROUTINE pre_physics                                                        &
   ! (In)
   ( row_length, rows, model_levels, wet_levels, nfor, ichgf, qcl, qcf        &
   , ch_ug, ch_vg, ilscnt, f_coriolis, lcal360, daycount, stepcount           &
   , a_sw_radstep_diag, a_sw_radstep_prog, l_triffid, npft                    &
   ! (InOut)
   , u, v, ug_scm, vg_scm, npft_trif                                          &
   ! (Out)
   , co2_mmr, l_rad_step, l_rad_step_prog, l_rad_step_diag )

  USE scm_utils, ONLY:                                                        &
    zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_main_force, ONLY:                                                     &
    ug, vg, timestep, obs, geoforce, co2start, co2rate, co2end, ntrad1        &
  , l_geo_centred

  IMPLICIT NONE
!
! Description: This routine sets several variables before the call to
!              the two physics routines.
! Method: It takes the relevant parts from the two physics routines at
!         4.5
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: FORTRAN 90
!
!  ARGUMENTS WITH INTENT In

! Parameters

! Model dimensions
  INTEGER ::               &
    row_length             &
  , rows                   &
  , model_levels           &
  , wet_levels             &
  , npft

  INTEGER ::               &
    nfor                   &! Number of obs. forcing
  , ichgf                  &
  , ilscnt

! Arrays
  REAL ::                                      &
    qcl(row_length,rows,wet_levels)            &
  , qcf(row_length,rows,wet_levels)            &
  , ch_ug(row_length,rows,nfor-1,model_levels) &
  , ch_vg(row_length,rows,nfor-1,model_levels) &
  , f_coriolis(row_length,rows)


! time information
  INTEGER ::               &
    daycount               &
  , stepcount              &
  , a_sw_radstep_diag      &
  , a_sw_radstep_prog


! Logicals

  LOGICAL ::               &
    lcal360                &
  , l_triffid

! ARGUMENTS WITH INTENT In/Out

  REAL ::                                   &
    u(row_length,rows,model_levels)         &
  , v(row_length,rows,model_levels)         &
  , ug_scm(row_length,rows,model_levels)    &
  , vg_scm(row_length,rows,model_levels)

  INTEGER ::               &
    npft_trif

! ARGUMENTS WITH INTENT Out

  REAL ::                  &
    co2_mmr

  LOGICAL ::               &
    l_rad_step             &
  , l_rad_step_prog        &
  , l_rad_step_diag

! local variables.

! loop counters
  INTEGER ::               &
    i, j, k
  REAL ::  utmp      ! Temporary u-wind

  ! Dr Hook
  !=============================================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('PRE_PHYSICS',zhook_in,zhook_handle)

! External routines:

! --------------------------------------------------------------------
  IF (l_triffid) THEN
    npft_trif = npft
  ELSE
    npft_trif = 1
  END IF

!
!---------------------------------------------------------------------
!     If Geostrophic forcing is chosen
!---------------------------------------------------------------------
!
  IF (geoforce) THEN

    IF (ilscnt > 0) THEN
      ug_scm (:,:,:) = ug_scm(:,:,:) + timestep * ch_ug(:,:,ilscnt,:)
      vg_scm (:,:,:) = vg_scm(:,:,:) + timestep * ch_vg(:,:,ilscnt,:)
    END IF

    IF (l_geo_centred) THEN
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            utmp = u(i,j,k)

            u(i,j,k) = ug_scm(i,j,k)                                     &
                + (  (u(i,j,k) - ug_scm(i,j,k))                          &
                   * (1.0 - (0.5*f_coriolis(i,j)*timestep)**2)           &
                   + (v(i,j,k) - vg_scm(i,j,k))                          &
                   * (f_coriolis(i,j)*timestep) )                        &
                / (1.0 + (0.5*f_coriolis(i,j)*timestep)**2)

            v(i,j,k) = vg_scm(i,j,k)                                     &
                + (  (v(i,j,k) - vg_scm(i,j,k))                          &
                   * (1.0 - (0.5*f_coriolis(i,j)*timestep)**2)           &
                   - (utmp-ug_scm(i,j,k))                                &
                   * (f_coriolis(i,j)*timestep) )                        &
                / (1.0 + (0.5*f_coriolis(i,j)*timestep)**2)
          END DO
        END DO
      END DO
    ELSE
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length

!         Store current u-wind to ensure temporally consistent
!         updating of v-wind.
            utmp = u(i,j,k)
            u(i,j,k) = u(i,j,k)                                          &
                      - f_coriolis(i,j) * timestep                       &
                                        * (vg_scm(i,j,k)-v(i,j,k))
            v(i,j,k) = v(i,j,k)                                          &
                      + f_coriolis(i,j) * timestep                       &
                                        * (ug_scm(i,j,k)-utmp)
          END DO
        END DO
      END DO
    END IF            ! l_geo_centred
  END IF            ! geoforce


!---------------------------------------------------------------------
!     Calculate the CO2 mass mixing ratio using the rate of change
!     (per year)
!---------------------------------------------------------------------

  IF (lcal360) THEN
    co2_mmr = co2start + co2rate                                    &
      * ((daycount-1)*86400 + (stepcount-1)*timestep)               &
      / 360*86400
  ELSE
    co2_mmr = co2start + co2rate                                    &
      * ((daycount-1)*86400 + (stepcount-1)*timestep)               &
      / 365*86400
  END IF

  IF (co2_mmr > co2end) THEN
    co2_mmr = co2end
  END IF

!---------------------------------------------------------------------
!     Is this a radiation timestep?
!---------------------------------------------------------------------

  l_rad_step_diag = .FALSE.
  l_rad_step_prog = .FALSE.

  IF (stepcount == ntrad1) THEN
    l_rad_step_diag = .TRUE.
    l_rad_step_prog = .TRUE.
  END IF

  IF (MOD(stepcount-ntrad1, a_sw_radstep_diag) == 0) THEN
    l_rad_step_diag = .TRUE.
  END IF
  IF (MOD(stepcount-ntrad1, a_sw_radstep_prog) == 0) THEN
    l_rad_step_prog = .TRUE.
  END IF

  IF (lhook) CALL dr_hook('PRE_PHYSICS',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE pre_physics

