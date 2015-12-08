! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE stph_skeb1_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE stph_skeb1( rows, row_length, n_rows                         &
 ,                     model_levels                                     &
 ,                     mass, c_dadt                                     &
 ,                     u, v                                             &
 ,                     kdisp                                            &
                     )

! Calculates energy dissipation based on the SKEB1 approach using
! Kinetic Energy of the wind. The dissipation field is added to the
! Smagorinsky and Convection fields, to modulate the SKEB2 stream-
! function random pattern.

 USE dynamics_input_mod, ONLY: l_endgame

! SKEB2 UMUI settings passed in via NameList READ
 USE stochastic_physics_run_mod, ONLY:                                  &
     skeb2_toplev, skeb2_botlev, tot_backscat, br, l_skebprint

! Bounds of arrays
 USE atm_fields_bounds_mod, ONLY:                                       &
     pdims, udims, vdims, udims_s, vdims_s

! Variables related to MPP
 USE proc_info_mod,     ONLY :                                          &
     global_row_length, global_rows

! Routines for global sums
 USE global_2d_sums_mod, ONLY: global_2d_sums

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim

 USE PrintStatus_mod
 USE UM_ParVars
 IMPLICIT NONE

! Variables with Intent (In)
 INTEGER, INTENT(IN) ::                                                 &
       row_length                                                       &
             ! number of points on a row.
     , rows                                                             &
             ! number of rows.
     , n_rows                                                           &
             ! number of v rows.
     , model_levels
             ! number of model levels.

 REAL, INTENT(IN) ::                                                    &
       u     (udims_s%i_start:udims_s%i_end,                            &
              udims_s%j_start:udims_s%j_end,                            &
              udims_s%k_start:udims_s%k_end)                            &
             ! input u-wind
     , v     (vdims_s%i_start:vdims_s%i_end,                            &
              vdims_s%j_start:vdims_s%j_end,                            &
              vdims_s%k_start:vdims_s%k_end)                            &
             ! input v-wind
     , c_dadt(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end)                                &
             ! inverse dx*dy*dt
     , mass  (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)
             ! mass of 3D grid-box

! Variables with Intent (Out)
 REAL, INTENT(OUT) ::                                                   &
       kdisp (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)
             ! (d/dx + d/dy + d/dz) KE Dissipation (m**2.s**-3)

! Local parameters
 INTEGER ::  i
             ! loop over row_length
 INTEGER ::  j
             ! loop over rows
 INTEGER ::  k
             ! loop over model levels
 INTEGER ::  ju
             ! j on u-row
 INTEGER ::  jv
             ! j on v-row
 INTEGER ::  jm1
             ! j minus one
 INTEGER ::  im1
             ! i minus one
 INTEGER ::  istat
             ! return error code
 REAL ::     gltotke
             ! Total global KE (used for calibration of kdisp)
 REAL ::     gltotke_tmp(1)
             ! Temporary array for global sum
 REAL ::     ke_vert(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end)
             ! Mass-weighted vertically integrated KE

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_SKEB2',zhook_in,zhook_handle)

! Initialise kdisp array
!
 DO k = pdims%k_start, skeb2_botlev - 1
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       kdisp(i,j,k) = 0.0
     END DO
   END DO
 END DO
 DO k = skeb2_toplev + 1, pdims%k_end
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       kdisp(i,j,k) = 0.0
     END DO
   END DO
 END DO

! Swap_bounds to get array values in u & v halo
 ! DEPENDS ON: swap_bounds
 CALL swap_bounds( u, row_length, rows,                                 &
                   skeb2_toplev, offx, offy, fld_type_u, .TRUE.  )

 ! DEPENDS ON: swap_bounds
 CALL swap_bounds( v, row_length, n_rows,                               &
                   skeb2_toplev, offx, offy, fld_type_v, .TRUE.  )

! Calculate point KE and globally (volume) integrated KE Flux Rate
! Note: mass = rho * dA * dz
 DO j = pdims%j_start, pdims%j_end
   DO i = pdims%i_start, pdims%i_end
      ke_vert(i,j) = 0.
   END DO
 END DO
 DO k = skeb2_botlev, skeb2_toplev
   DO j = pdims%j_start, pdims%j_end
     jm1 = j - 1
     ju = j
     jv = j
     IF (at_extremity(psouth)) jm1 = MAX(jm1, udims%j_start)
     IF (at_extremity(pnorth)) ju = MIN(ju, udims%j_end - 1)
     IF (at_extremity(pnorth)) jv = MIN(jv, vdims%j_end)

     DO i = pdims%i_start, pdims%i_end
       im1= i - 1
       kdisp(i,j,k) =  ((0.5*( u(im1,ju, k) + u(i, ju, k)))**2  +       &
                        (0.5*( v(i, jm1, k) + v(i, jv, k)))**2)
       ke_vert(i,j) = ke_vert(i,j) + 0.5 * mass(i,j,k) * kdisp(i,j,k)
     END DO !I
   END DO !J
 END DO !K
IF (.NOT. l_endgame) THEN
! Set value at pole = zero (single point with no mass)
 IF (at_extremity(psouth)) THEN
   DO i = pdims%i_start, pdims%i_end
       ke_vert(i,pdims%j_start) = 0.
   END DO
 END IF
 IF (at_extremity(pnorth)) THEN
   DO i = pdims%i_start, pdims%i_end
       ke_vert(i,pdims%j_end) = 0.
   END DO
 END IF
END IF

 DO j = pdims%j_start, pdims%j_end
   DO i = pdims%i_start, pdims%i_end
      ke_vert(i,j) = ke_vert(i,j) * c_dadt(i,j)
   END DO
 END DO

! Calculate global KE flux (sum from all processors)
 CALL global_2d_sums(ke_vert, row_length, rows, 0, 0, 1, gltotke_tmp)

 gltotke = gltotke_tmp(1)/(global_row_length*global_rows)

 IF (l_skebprint) THEN
   WRITE(6,*)
   WRITE(6,'("***** SKEB1 PE-ave KE *****")')
   WRITE(6,'("tot KE Rate (W/m^2)= ",ES22.15)') gltotke
   WRITE(6,'("Calibr. Fac = ",ES22.15)') tot_backscat/(br*gltotke)
 END IF

!
! Calibrate kdisp array to yield "tot_backscat/br" W/m^2
!
 DO k = skeb2_botlev, skeb2_toplev
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       kdisp(i,j,k) = kdisp(i,j,k)*tot_backscat/(br*gltotke)
     END DO
   END DO
 END DO

 IF (lhook) CALL dr_hook('STPH_SKEB1',zhook_out,zhook_handle)
 RETURN

END SUBROUTINE stph_skeb1
END MODULE stph_skeb1_mod
