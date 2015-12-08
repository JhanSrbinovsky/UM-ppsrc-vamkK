! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate level of neutral buoyancy, CAPE and CIN for ascent
!
! Subroutine Interface:
SUBROUTINE cv_parcel_neutral_inv(nunstable,model_levels,wet_levels,           &
               nlcl_c,k_plume,                                                &
               z_lcl_c, thv_pert, z_full_c, z_half_c, exner_theta_levels_c,   &
               buoyancy, t_dens_env, dqsatdz, denv_bydz, dpar_bydz,           &
               zh_c, zh2,                                                     &
               k_max,k_neutral,k_inv, kshmin, shmin,                          &
               max_buoy,dt_dens_parc_t,ql_ad_c,delthvu_c,cape_c,cin_c,        &
               dt_dens_parc_t2,dt_dens_parc_tmin,ql_ad2,                      &
               delthvu2)


USE earth_constants_mod, ONLY: g

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
!
! Description:
!   Calculate the level of neutral buoyancy for a parcel ascent.
!   Also calcuate the CAPE and CIN of the ascent
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER, INTENT(IN) :: &
  nunstable            & ! Number of parcel ascents
 ,model_levels         & ! Number of model levels
 ,wet_levels             ! Number of model wet levels

INTEGER, INTENT(IN) :: &
  nlcl_c(nunstable)    & ! Level number of LCL
 ,k_plume(nunstable)     ! start level for surface-driven plume

REAL, INTENT(IN)    ::   &
  z_lcl_c(nunstable)     & ! LCL height rounded to nearest model level
 ,thv_pert(nunstable)      ! threshold thv of parcel  (K)

REAL, INTENT(IN)    ::                           &
  z_full_c(nunstable, model_levels)              & ! Height theta lev (m)
 ,z_half_c(nunstable, model_levels)              & ! Height uv lev    (m)
 ,exner_theta_levels_c(nunstable, model_levels)  & ! Exner on theta lev
 ,buoyancy(nunstable, wet_levels)            & ! undilute parcel buoyancy (K)
 ,t_dens_env(nunstable, wet_levels)          & ! Density potential temperature 
                                               ! of environment. (K)
 ,dqsatdz(nunstable, wet_levels)             & ! dqsat/dz along adiabat
 ,denv_bydz(nunstable, wet_levels)           & ! Gradient of density potential
                                               ! temp in the environment.
 ,dpar_bydz(nunstable, wet_levels)             ! Gradient of density potential
                                               ! temperature of the parcel.


REAL, INTENT(INOUT)    :: &
  zh_c(nunstable)         & ! BL depth compressed
 ,zh2(nunstable)            !  2nd copy

INTEGER, INTENT(OUT) :: &
  k_max(nunstable)      & ! level of max parcel buoyancy
 ,k_neutral(nunstable)  & ! level of neutral parcel buoyancy
 ,k_inv(nunstable)      & ! level from inversion testing
 ,kshmin(nunstable)       ! Position of buoyancy minimum above
                          ! topbl (used for shallow Cu diag)level 

LOGICAL,  INTENT(OUT) :: &
  shmin(nunstable)         ! Flag for finding min in parcel buoyancy below
                           ! 3km (for shallow Cu)
REAL, INTENT(OUT)    ::  &
  max_buoy(nunstable)      ! Maximum buoyancy

REAL, INTENT(OUT)   ::         &
  dt_dens_parc_T(nunstable)    & ! t_dens_parc-t_dens_env at ntpar
 ,dt_dens_parc_Tmin(nunstable) & ! t_dens_parc-t_dens_env at kshmin
 ,cape_c(nunstable)            & ! CAPE from undilute parcel ascent (m2/s2)
 ,cin_c(nunstable)             & ! CIN from undilute parcel ascent (m2/s2) 
 ,delthvu_c(nunstable)         & ! Integral of undilute parcel buoyancy
                                 ! over convective cloud layer (for convection)
 ,ql_ad_c(nunstable)           & ! adiabatic liquid water content at inversion
                                 ! or cloud top (kg/kg)
 ,dt_dens_parc_T2(nunstable)   & ! t_dens_parc-t_dens_env at ntpar
 ,delthvu2(nunstable)          & ! 2nd copy delthuv
 ,ql_ad2(nunstable)              ! ql_ad 2nd copy

!-------------------------------------------------------------------------
! Local variables

INTEGER :: ii, k            !Loop counters 


REAL    :: &
  inc      & ! CAPE in layer
 ,dz         ! layer thickness

REAL    ::                &
  dtv_min(nunstable)      & ! min Tv of parcel in cld layer 1 virtual 
                            ! temperature (K).
 ,dtv_min2(nunstable)       ! 2nd copy min TV of parcel

LOGICAL ::                &
  topbl(nunstable)        & ! Flag set when top of boundary layer is reached.
 ,topprof(nunstable)      & ! Flag set when top of ascent is reached.
 ,above_lcl(nunstable)      ! Flag set when parcel above LCL.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-------------------------------------------------------------------------
! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('cv_parcel_neutral_inv',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialise output arrays

DO ii=1,nunstable
  shmin(ii)   = .FALSE.
  topbl(ii)   = .FALSE.
  topprof(ii) = .FALSE.

  kshmin(ii)     = 1
  k_max(ii)      = 1
  k_neutral(ii)  = 1
  k_inv(ii)      = 1

  dtv_min(ii)  = 0.0
  dtv_min2(ii) = 0.0
  delthvu2(ii) = 0.0
  max_buoy(ii)     = 0.0
  delthvu_c(ii) = 0.0
  ql_ad2(ii)   = 0.0
  ql_ad_c(ii)  = 0.0
  CAPE_c(ii)  = 0.0
  CIN_c(ii)  = 0.0
END DO

  DO  K = 2,wet_levels

    DO ii=1,nunstable


     ! Set flag to true when level below is at least one level above the lcl
     ! and above the lcl transition zone
     ! Code implies ABOVE_LCL at NLCL+3 or greater.

      IF (k-1 >  nlcl_c(ii)+1                                            &
                       .AND. z_full_c(ii,k-1) >  1.1*z_lcl_c(ii)) THEN
        above_lcl(ii)=.TRUE.
      ELSE
        above_lcl(ii)=.FALSE.
      END IF

     !-----------------------------------------------------------------------
     ! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
     !-----------------------------------------------------------------------
     ! Not reached LNB continue testing

      IF ( .NOT.topprof(ii).AND.k >  k_plume(ii) )THEN

        IF (buoyancy(ii,k) >  max_buoy(ii)) THEN
          max_buoy(ii) = buoyancy(ii,k)
          k_max(ii)    = k  
        END IF 

        ! Is parcel still buoyant ?

        IF ( (buoyancy(ii,k)  <=  - thv_pert(ii))                        &
        !                      or reached top of model
                 .OR. (K  >   wet_levels-1)  ) THEN

          k_neutral(ii) = k-1
          topprof(ii) = .TRUE.
          zh_c(ii) = z_half_c(ii,K)

        ! Buoyancy at last buoyant level

          Dt_dens_parc_T(ii) = buoyancy(ii,k-1)

          IF ( delthvu_c(ii)  >   0.0) THEN
          ! compensate for any negative buoyancy of parcel in cloud layer
            delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *                 &
                                      ( z_half_c(ii,K) - z_lcl_c(ii) )
          END IF                                                     
        END IF
      END IF

      !-----------------------------------------------------------------------
      ! Tests applied once found top of parcel ascent.
      ! Aim - to establish if the ascent has an inversion above the top
      !       i.e. the ascent may indicate shallow /congestus convection.
      ! Sets indicator shmin = .true. if conditions met and stops testing.
      !
      ! Conditions are ;
      ! either  denv/dz(k) < dpar/dz(k)
      !   or    par_svl(k-1) -env_svl(k-1) <= 0.0
      !
      !-----------------------------------------------------------------------

      IF ( topbl(ii) .AND. ( denv_bydz(ii,k)  <   dpar_bydz(ii,k)         &
                      .OR. buoyancy(ii,k-1)  <=  0.0 )                    &
                                    .AND. .NOT. shmin(ii) ) THEN
        shmin(ii) = .TRUE.
        dt_dens_parc_tmin(ii) = buoyancy(ii,k-1)
        kshmin(ii) = K-1
      END IF

      !-----------------------------------------------------------------------
      ! Tests applied to find parcel top
      !-----------------------------------------------------------------------

      IF ( .NOT.topbl(ii) .AND. K  >   k_plume(ii) .AND.                  &
          (  ( buoyancy(ii,k) <=  - thv_pert(ii)).OR.                     &

      !           plume non buoyant

          (above_lcl(ii).AND.(denv_bydz(ii,k) >  1.25*dpar_bydz(ii,k)))   &

      !           or environmental virtual temperature gradient
      !           signIficantly larger than parcel gradient
      !           above lIfting condensation level

                 .OR. (K  >   wet_levels-1)                               &
      !                      or reached top of model
               )                                                          &
               ) THEN

        topbl(ii) = .TRUE.
        zh2(ii) = z_half_c(ii,k)
        k_inv(ii) = k-1

        dt_dens_parc_T2(ii) = buoyancy(ii,k-1)
        IF ( delthvu2(ii)  >   0.0) THEN
        ! compensate for any negative buoyancy of parcel in cloud layer
          delthvu2(ii) = delthvu2(ii) - dtv_min2(ii) *             &
                                    ( z_half_c(ii,k) - z_lcl_c(ii) )
        END IF
      END IF          ! test on .not.topbl

      !-----------------------------------------------------------------------
      ! While doing parcel ascent
      ! (a) find minimum buoyancy
      ! (b) integrate CAPE over the ascent
      !-----------------------------------------------------------------------

      IF (k > nlcl_c(ii) .AND. k < wet_levels ) THEN

        dz = z_half_c(ii,k+1) - z_half_c(ii,k)
        inc = g *  buoyancy(ii,k) * dz/t_dens_env(ii,k)

        !---------------------------------------------------------- 
        ! If not reached an inversion or level of neutral buoyancy
        !---------------------------------------------------------- 

        IF (.NOT. topbl(ii)) THEN


        ! adiabatic liquid water content at cloud top or inversion
        !                            = -dqsat/dz * zcld

          ql_ad2(ii) = -1.* dqsatdz(ii,k)                                   &
                          *(z_half_c(ii,K+1) - z_half_c(ii,nlcl_c(ii)+1))

          dtv_min2(ii) = MIN( dtv_min2(ii),                                 &
                             buoyancy(ii,k)/exner_theta_levels_c(ii,k) )

          delthvu2(ii) = delthvu2(ii) + buoyancy(ii,k)*dz                   &
                                             /exner_theta_levels_c(ii,k)
        END IF

        !---------------------------------------------------------- 
        ! If not reached level of neutral buoyancy (top of ascent)
        !---------------------------------------------------------- 

        IF (.NOT. topprof(ii)) THEN
                
        ! Note only calculating CIN and CAPE from ascents reaching
        ! level of neutral buoyancy. This may not always correspond 
        ! to the diagnosed top for the convection scheme.

          IF (inc <  0.0) THEN
            CIN_c(ii)  = CIN_c(ii) + inc
          ELSE      ! CAPE holds only postive part
            CAPE_c(ii) = CAPE_c(ii) + inc
          END IF
          ! adiabatic liquid water content at cloud top
          !                            = -dqsat/dz * zcld
 
          ql_ad_c(ii) = -1.* dqsatdz(ii,k)                        &
                       *(z_half_c(ii,k+1) - z_half_c(ii,nlcl_c(ii)+1))

          dtv_min(ii) = MIN( dtv_min(ii),                         &
                           buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

          delthvu_c(ii) = delthvu_c(ii) + buoyancy(ii,k)* dz      &
                                            /exner_theta_levels_c(ii,k)
 
        END IF    ! test on topprof

      END IF

    END DO   ! ii loop

  END DO     ! level loop


IF (lhook) CALL dr_hook('cv_parcel_neutral_inv',zhook_out,zhook_handle)
RETURN
END SUBROUTINE cv_parcel_neutral_inv
