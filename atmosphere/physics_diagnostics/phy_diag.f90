! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate various diagnostics related to dynamics variables.
!
! Subroutine Interface:
MODULE phy_diag_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE Phy_diag(                                              &
! Primary data: in
     & p_star,p,rho,u,v,w,theta,q,qCL,qCF                               &
     &,p_theta_levels,exner_rho_levels,exner_theta_levels               &
! Grid sizes and definition: in
     &,rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,global_row_length, global_rows                                   &
     &,delta_lambda,delta_phi,sec_theta_latitude                        &
! Control information: in
     &,Model_domain,npmsl_height,l_pmsl_sor                             &
! Pressure levels for output arrays: in
     &,hts_p_press,t_p_press                                            &
     &,RHice_p_press,RHwat_p_press,wbpt_p_press                         &
! Flags to request each diagnostic output field: in
     &,qT_m                                                             &
     &,qhts_theta                                                       &
     &,qhts_p,qt_p,qRHice_p,qwbpt_p                                     &
     &,qqc,qqT                                                          &
     &,qp_MSL                                                           &
     &,qhts_rho                                                         &
     &,qRHwat_p                                                         &
! Diagnostics lengths: in
     &,hts_p_levs,t_p_levs                                              &
     &,RHice_p_levs,RHwat_p_levs,wbpt_p_levs                            &
! Diagnostic arrays: out
     &,T                                                                &
     &,hts_theta                                                        &
     &,hts_p,t_p,RHice_p,wbpt_p                                         &
     &,qc,qT                                                            &
     &,p_MSL                                                            &
     &,hts_rho                                                          &
     &,RHwat_p                                                          &
     &     )

      USE atm_fields_bounds_mod, ONLY:                                  &
          pdims, pdims_s, udims_s, vdims_s, wdims_s,                    &
          tdims_s, wdims_l, qdims_l, pdims_l

      USE level_heights_mod, ONLY:                                      &
                r_theta_levels, r_rho_levels
                
      USE earth_constants_mod, ONLY: g, earth_radius

      USE atmos_constants_mod, ONLY: kappa, p_zero

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE interpor_mod, ONLY: interp_order_linear
      USE thetaw_mod, ONLY: thetaw
      USE calc_pmsl_mod, ONLY: calc_pmsl
      USE t_vert_interp_to_p_mod, ONLY: t_vert_interp_to_p
      USE vert_h_onto_p_mod, ONLY: vert_h_onto_p
      USE vert_interp2_mod, ONLY: vert_interp2
      IMPLICIT NONE
!
! Description:
!   Calculate physics-related diagnostics - held in STASH section 16 -
!   which may include interpolation onto pressure surfaces. Diagnostics
!   currently supported:
!   STASH item
!     4 temperature on model levels
!   201 geopotential height on theta levels
!   202 geopotential height on pressure surfaces
!   203 temperature         on pressure surfaces
!   204 relative humidity wrt ice on pressure surfaces
!   205 wet bulb potential temperature on pressure surfaces
!   206 cloud water content     (qc) on q levels
!   207 total specific humidity (qT) on q levels
!   222 mean sea level pressure
!   255 geopotential height on rho levels
!   256 relative humidity wrt water on pressure surface
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Primary model data is input, and each diagnostic is calculated - if
!   its flag is set - in simple sequential order. Where the extraction
!   of the diagnostic quantity requires further calculation, a lower
!   level diagnostic-specific routine is called.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Physics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!   Documentation: UMDP 80
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
! Grid sizes:
      INTEGER                                                           &
     & rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,global_row_length, global_rows                                   &
! Control information: in
     &,Model_domain                                                     &
                       ! Domain of atmosphere model:
!                                   global,LAM,cyclic LAM,single column

! Diagnostics lengths: IN
     &,hts_p_levs                                                       &
                      ! NO OF LEVS ON WHICH TO INTERP hts_p
     &,t_p_levs                                                         &
                      ! NO OF LEVS ON WHICH TO INTERP t_p
     &,RHice_p_levs                                                     &
                      ! NO OF LEVS ON WHICH TO INTERP RHice_p
     &,RHwat_p_levs                                                     &
                      ! NO OF LEVS ON WHICH TO INTERP RHwat_p
     &,wbpt_p_levs    ! NO OF LEVS ON WHICH TO INTERP wbpt_p

! Grid definition: IN
      REAL                                                              &
     & delta_lambda,delta_phi                                           &
! Orographic height threshold for new pmsl calculation: IN
     &,npmsl_height

      LOGICAL :: l_pmsl_sor

! Flags to request each diagnostic output field: IN
      LOGICAL                                                           &
     & qhts_p                                                           &
                   ! Flag for geopotential heights   on pressure levels
     &,qt_p                                                             &
                   ! Flag for temperature            on pressure levels
     &,qRHice_p                                                         &
                   ! Flag for relatice humidity wrt ice
                   !   on pressure surfaces
     &,qRHwat_p                                                         &
                   ! Flag for relative humidity wrt water
                   !   on pressure surfaces
     &,qwbpt_p                                                          &
                   ! Flag for wet bulb potential T   on pressure levels
     &,qqc                                                              &
                   ! Flag for cloud water content (qc)
     &,qqT                                                              &
                   ! Flag for total specific humidity (qT)
     &,qp_MSL                                                           &
                   ! Flag for mean sea level pressure
     &,qT_m                                                             &
                   ! Flag for temperature on model levels
     &,qhts_theta                                                       &
                   ! Flag for geopotential heights   on theta levels
     &,qhts_rho    ! Flag for geopotential heights   on rho levels
!   Array  arguments with intent(in):
! Primary data: IN
      REAL                                                              &
       p_star(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end)                                &
      ,p     (pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end)                            &
      ,rho   (pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end)                            &
      ,u     (udims_s%i_start:udims_s%i_end,                            &
              udims_s%j_start:udims_s%j_end,                            &
              udims_s%k_start:udims_s%k_end)                            &
      ,v     (vdims_s%i_start:vdims_s%i_end,                            &
              vdims_s%j_start:vdims_s%j_end,                            &
              vdims_s%k_start:vdims_s%k_end)                            &
      ,w     (wdims_s%i_start:wdims_s%i_end,                            &
              wdims_s%j_start:wdims_s%j_end,                            &
              wdims_s%k_start:wdims_s%k_end)                            &
      ,theta (tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end)                            &
      ,q     (qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end)                            &
      ,qCL   (qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end)                            &
      ,qCF   (qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end)                            &
      ,p_theta_levels                                                   &
             (tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end)                            & 
      ,exner_rho_levels                                                 &
             (pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end + 1)                        &
      ,exner_theta_levels                                               &
             (tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end)                            &
! Derived from horizontal grid: IN
     &,sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)      &

! Pressure levels for output arrays: IN
     &,hts_p_press(hts_p_levs)                                          &
                                      ! for heights      on p surfaces
     &,t_p_press  (  t_p_levs)                                          &
                                      ! for temperature  on p surfaces
     &,RHice_p_press ( RHice_p_levs)                                    &
                                      ! for r.h. wrt ice on p surf
     &,RHwat_p_press ( RHwat_p_levs)                                    &
                                      ! for r.h. wrt water on p surf
     &,wbpt_p_press(wbpt_p_levs)      ! for wet bulb T   on p surfaces

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
! Diagnostic arrays: OUT
      REAL                                                              &
     & hts_p(row_length,rows,hts_p_levs)                                &
                                         ! heights           at p levels
     &,t_p  (row_length,rows,t_p_levs)                                  &
                                         ! temperature       at p levels
     &,RHice_p (row_length,rows,RHice_p_levs)                           &
                                              ! rh wrt ice at p levels
     &,RHwat_p (row_length,rows,RHwat_p_levs)                           &
                                              ! rh wrt water at p-levels
     &,wbpt_p(row_length,rows,wbpt_p_levs)                              &
                                          !wet bulb pot temp at p levels
     &,qc(row_length,rows,wet_levels)                                   &
                                   ! cloud water content on q levels
     &,qT(row_length,rows,wet_levels)                                   &
                                   ! total specific humidity on q levels
     &,p_MSL(row_length,rows)                                           &
                                         ! pressure at mean sea level
     &,T(row_length,rows,model_levels)                                  &
                                         ! temperature on model levels
     &,hts_theta(row_length,rows,model_levels)                          &
                                         ! heights on theta levels
     &,hts_rho(row_length,rows,model_levels)
                                         ! heights on rho levels
! Local parameters:
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='Phy_diag')
      LOGICAL, PARAMETER ::                                             &
     & l_potential=.true.              ! Wet bulb potential temperature
                                       ! required from subroutine ThetaW
! Local scalars:
!   ErrorStatus
      INTEGER      ErrorStatus                                          &
                                        ! Error flag (0 = OK)
     &,i,j,k                                                            &
                                        ! loop counters
     &,interp_order    !  order of vertical interpolation

      CHARACTER(LEN=256)                                                     &
     & CMessage         ! Error message if return code >0

      REAL                                                              &
     & dummy                                                            &
                                ! dummy argument - not referenced
     &,pressure_pa                                                      &
                                ! pressure in pascals
     &,pressure_ex              ! exner pressure

! Local dynamic arrays:
      REAL                                                              &
     & p_at_theta                                                       &
     &       (row_length,rows,model_levels)                             &
                                            ! pressure at theta points
     &,exner_at_theta                                                   &
     &       (row_length,rows,model_levels)                             &
                                            ! exner pressure at theta points
     &,rh    (row_length,rows,model_levels)                             &
                                            ! workspace for RH
     &,height(row_length,rows,model_levels)                             &
                                            ! height at rho levels
     &,work1  (row_length,rows)                                         &
                                            ! temporary space
     &,work2  (row_length,rows)                                         &
                                            ! temporary space
     &,work3 (row_length,rows)              ! temporary space

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!- End of header


! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('PHY_DIAG',zhook_in,zhook_handle)

! set Error code to zero
      ErrorStatus = 0

! Set order of vertical interpolation
      interp_order = interp_order_linear

!  Calculate p at theta points. Store in p_at_theta

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            p_at_theta(i,j,k) = p_theta_levels(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

!  Calculate exner at theta points. Store in exner_at_theta

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            exner_at_theta(i,j,k) = exner_theta_levels(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

!   Remove radius of earth from rho levels to create geopotential height
!   of rho level.
!   Required for either 202 (geopotential height on pressure levels) or
!                       255 (geopotential height on rho levels)
      IF (qhts_p .OR. qhts_rho) THEN

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            height(i,j,k) = r_rho_levels(i,j,k) - earth_radius
          END DO ! i
        END DO ! j
      END DO ! k

      ENDIF ! on relevant STASHflags

!   Calculate temperature at theta points
      IF(qt_p .OR. qRHice_p .OR. qRHwat_p .OR. qT_m .OR. qwbpt_p) THEN

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

      ENDIF ! on relevant STASHflags

! ----------------------------------------------------------------------
! STASH item 201: geopotential height      on  theta levels
! ----------------------------------------------------------------------
      IF(qhts_theta) THEN
!   Remove radius of earth from height field
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              hts_theta(i,j,k) = r_theta_levels(i,j,k) - earth_radius
            END DO ! i
          END DO ! j
        END DO ! k
      ENDIF ! on relevant STASHflags

! ----------------------------------------------------------------------
! STASH item 255: geopotential height on  rho levels
! ----------------------------------------------------------------------
      IF(qhts_rho) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              hts_rho(i,j,k) = height(i,j,k)
            END DO ! i
          END DO ! j
        END DO ! k
      ENDIF ! on relevant STASHflags
! ----------------------------------------------------------------------
! STASH item 202: geopotential height      on pressure surfaces
! ----------------------------------------------------------------------
      IF(qhts_p) THEN
        DO k = 1, hts_p_levs
         pressure_pa = hts_p_press(k)*100. ! convert to Pascals
         CALL vert_h_onto_p(                                            &
     &        height(1,1,1), row_length,rows, model_levels              &
     &       ,pressure_pa                                               &
     &       ,p_theta_levels                                            &
     &       ,theta, exner_theta_levels, exner_rho_levels               &
     &       , bl_levels                                                &
     &       ,offx, offy, halo_i, halo_j                                &
     &       ,p, interp_order, hts_p(1,1,k) )
        END DO ! over output STASH pressure levels
      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 203: temperature              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qt_p) THEN
        DO k = 1, t_p_levs
         pressure_pa = t_p_press(k)*100. ! convert to Pascals
         CALL T_vert_interp_to_p(                                       &
     &        T, theta, row_length, rows                                &
     &        ,model_levels, pressure_pa, offx, offy,halo_i,halo_j      &
     &        ,p_theta_levels                                           &
     &        ,bl_levels                                                &
     &        ,exner_theta_levels                                       &
     &        ,T_p(1,1,k) )
        END DO ! over output STASH pressure levels
      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 204, 256 : RH: relative humidity on pressure surfaces
! ----------------------------------------------------------------------

! STASH item 204: Relative humidity wrt ice
      IF(qRHice_p) THEN

        DO k = 1, wet_levels

!  Find humidity saturation at theta points - store in rh
! DEPENDS ON: qsat
          CALL QSAT(rh(1,1,k),T(1,1,k),p_at_theta(1,1,k),               &
     &              theta_field_size)
!  And convert to relative humidity
          DO j = 1, rows
            DO i = 1, row_length
              rh(i,j,k) = q(i,j,k)/rh(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
              IF(rh(i,j,k) <  0.0) THEN
                rh(i,j,k) = 0.
              END IF
            END DO ! i
          END DO ! j

        END DO ! k wet_levels

!  Interpolate
        DO k = 1, RHice_p_levs
          pressure_pa = RHice_p_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          CALL vert_interp2(                                            &
     &          rh, row_length, rows, wet_levels                        &
     &        , pressure_ex                                             &
     &        , 0, 0, 0, 0                                              &
     &        , exner_at_theta, interp_order                            &
     &        , RHice_p(1,1,k) )
        END DO ! k over output STASH pressure levels

      END IF  ! on STASHflag

! STASH item 256: Relative humidity wrt water
      IF(qRHwat_p) THEN

        DO k = 1, wet_levels
!  Find humidity saturation at theta points - store in rh
! DEPENDS ON: qsat_wat
          CALL QSAT_WAT(rh(1,1,k),T(1,1,k),p_at_theta(1,1,k),           &
     &              theta_field_size)
!  And convert to relative humidity
          DO j = 1, rows
            DO i = 1, row_length
              rh(i,j,k) = q(i,j,k)/rh(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
              IF(rh(i,j,k) <  0.0) THEN
                rh(i,j,k) = 0.
              END IF
            END DO ! i
          END DO ! j

        END DO ! k wet_levels

!  Interpolate
        DO k = 1, RHwat_p_levs
          pressure_pa = RHwat_p_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          CALL vert_interp2(                                            &
     &          rh, row_length, rows, wet_levels                        &
     &        , pressure_ex                                             &
     &        , 0, 0, 0, 0                                              &
     &        , exner_at_theta, interp_order                            &
     &        , RHwat_p(1,1,k) )
        END DO ! k over output STASH pressure levels

      END IF  ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 205: wet bulb potential temperature on pressure surfaces
! ----------------------------------------------------------------------
      IF(qwbpt_p) THEN

        DO k = 1, wbpt_p_levs

         pressure_pa = wbpt_p_press(k)*100. ! convert to Pascals

! Interpolate T onto required pressure level (work1)
         CALL T_vert_interp_to_p(                                       &
     &        T, theta, row_length, rows                                &
     &        ,model_levels, pressure_pa, offx, offy, halo_i, halo_j    &
     &        ,p_theta_levels                                           &
     &        ,bl_levels                                                &
     &        ,exner_theta_levels                                       &
     &        ,work1 )

! Interpolate q onto required pressure level (work2)
          pressure_ex = ( pressure_pa /p_zero )**kappa
          CALL vert_interp2(                                            &
     &         q, row_length, rows, wet_levels                          &
     &        ,pressure_ex                                              &
     &        ,halo_i, halo_j, 0, 0                                     &
     &        ,exner_at_theta, interp_order                             &
     &        ,work2 )

! Generate pressure array for required pressure level (work3)
         DO j = 1, rows
           DO i = 1, row_length
             work3(i,j) = pressure_pa
           END DO
         END DO
         CALL Thetaw(                                                   &
     &        theta_field_size,work1,work2,work3,l_potential,           &
                                                               ! in
     &        wbpt_p(1,1,k))                                  ! out

        END DO ! over output STASH pressure levels

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 206: cloud water content (qc) on q levels
! ----------------------------------------------------------------------
      IF(qqc) THEN
        DO k = 1, wet_levels
          DO j = 1, rows
            DO i = 1, row_length
              qc(i,j,k) = qCL(i,j,k) + qCF(i,j,k)
            END DO
          END DO
        END DO
      ENDIF ! on STASHflag
! ----------------------------------------------------------------------
! STASH item 207: total specific humidity (qT) on q levels
! ----------------------------------------------------------------------
      IF(qqT) THEN

        DO k = 1, wet_levels
          DO j = 1, rows
            DO i = 1, row_length
              qT(i,j,k) = q(i,j,k) + qCL(i,j,k) + qCF(i,j,k)
            END DO
          END DO
        END DO

      ENDIF ! on STASHflag
! ----------------------------------------------------------------------
! STASH item 222: mean sea level pressure
! ----------------------------------------------------------------------
      IF(qp_MSL) THEN

         CALL Calc_PMSL(theta, exner_theta_levels, p                    &
     &                 ,row_length, rows, model_levels, bl_levels       &
     &                 ,offx, offy, halo_i, halo_j                      &
     &                 , p_MSL, p_star                                  &
     &                 ,mype, nproc, nproc_x, nproc_y                   &
     &                 ,neighbour, at_extremity                         &
     &                 ,gc_all_proc_group, model_domain                 &
     &                 ,delta_lambda, delta_phi                         &
     &                 ,npmsl_height, l_pmsl_sor, sec_theta_latitude    &
     &                 ,global_row_length, global_rows)

      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! Check error condition
      IF(ErrorStatus >  0) THEN

         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      IF (lhook) CALL dr_hook('PHY_DIAG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Phy_diag
END MODULE phy_diag_mod
