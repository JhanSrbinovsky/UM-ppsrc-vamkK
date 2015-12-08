! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE init_vert_damp_mod
IMPLICIT NONE
! Subroutine: init_vert_damp
!
! Description: Initialises the height dependent vertical damping 
!              coefficient
!
! Method: The vertical damping term motivated by the results of 
!         Klemp JB, Dudhia J, Hassiotis AD. 2008."An upper gravity-wave
!         absorbing layer for NWP applications." Mon. Wea. Rev. 136, 
!         3987-4004.
!
!         Following Temam and Tribbia (2003) this damping term is 
!         formulated so that it is active in the hydrostatic formulation.
!         Although this means that in the the sponge layer hydrostatic
!         balance will not be strictly enforced.
!
!         Melvin et.al (2009) suggest a piecewise defined profile:
!
!                 /
!                 | 0, \eta < \eta_s
!                 |
!         \mu_w = <
!                 !             \pi  (eta-eta_sa)
!                 | \mu_bar sin(---  ------------)^2, \eta >= \eta_s
!                 |              2    (1-eta_s)
!                 !

!
!        Note that in this implementation eta is defined in terms of
!        xi_3:
!
!              xi_3_at_theta-earth_radius
!        eta = --------------------------,
!               domain_height
!
!        i.e. eta is not necessarily zero at ground level.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: DYNAMICS ADVECTION
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

CONTAINS

SUBROUTINE init_vert_damp(eg_vert_damp_profile, eta_s,                &
                             eg_vert_damp_coeff,                      &
                             height_domain)

USE atm_fields_bounds_mod, ONLY : wdims
USE earth_constants_mod,   ONLY : earth_radius
USE level_heights_mod,     ONLY : xi3_at_theta => r_theta_levels
USE timestep_mod,          ONLY : timestep
USE conversions_mod,       ONLY : pi
USE eg_vert_damp_mod,      ONLY : mu_w

USE proc_info_mod,         ONLY : me

USE parkind1,              ONLY : jprb, jpim
USE yomhook,               ONLY : lhook, dr_hook
USE ereport_mod,           ONLY : ereport
USE PrintStatus_mod
USE horiz_grid_mod,        ONLY : xi2_p
USE missing_data_mod,      ONLY : rmdi
USE um_input_control_mod,  ONLY : model_domain
USE domain_params,         ONLY : mt_lam

IMPLICIT NONE

INTEGER, INTENT(IN) :: eg_vert_damp_profile
REAL,    INTENT(IN) :: eta_s ! "eta" level at which the damping layer starts
                             ! see note about this version of eta above!
REAL,    INTENT(IN) :: eg_vert_damp_coeff
REAL,    INTENT(IN) :: height_domain

!     local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*),PARAMETER :: RoutineName='EG_INIT_VERT_DAMP'    
INTEGER                 :: ErrorStatus                        
CHARACTER(LEN=80)          :: Message   

INTEGER :: i
INTEGER :: j
INTEGER :: k
INTEGER :: ierr

REAL    :: eta_ijk      ! "eta" in gridpoint
                        ! see note about this version of eta above!

REAL    :: recip_one_m_etas   ! 1./(1.-eta_s)


INTEGER, PARAMETER :: zero_vert_damp_prof     = 0 
INTEGER, PARAMETER :: seventyfive_vert_damp_prof = 1 
INTEGER, PARAMETER :: polar_damp_prof         = 2
INTEGER, PARAMETER :: const_damp_prof         = 3
INTEGER, PARAMETER :: polar_damp_prof2        = 4
INTEGER, PARAMETER :: polar_damp_prof2_narrow = 42
INTEGER, PARAMETER :: polar_damp_prof2_extra_narrow = 43
INTEGER, PARAMETER :: polar_damp_prof2_NP     = 41
INTEGER, PARAMETER :: polar_damp_prof3        = 5

REAL :: cos_xi2, cutoff, cutoff_deg, cutoff_rad

ErrorStatus=0 

IF (lhook) CALL dr_hook('INIT_VERT_DAMP',zhook_in,zhook_handle)

ALLOCATE (mu_w(  wdims%i_start:wdims%i_end,                           & 
                 wdims%j_start:wdims%j_end,                           &
                 wdims%k_start:wdims%k_end),STAT=ierr)

IF(ierr/=0) THEN
  WRITE(message,*) 'allocation of mu_w failed'
  CALL Ereport(RoutineName,ierr,message)
END IF
 
SELECT CASE (eg_vert_damp_profile)

! No damping ---------------------------------------------------------
  CASE(zero_vert_damp_prof)

     mu_w = 0

! Standard -----------------------------------------------------------
  CASE(seventyfive_vert_damp_prof)

    recip_one_m_etas = 1.0/(1.0-eta_s)

    DO k=wdims%k_start,wdims%k_end
      DO j=wdims%j_start,wdims%j_end
        DO i=wdims%i_start,wdims%i_end

          eta_ijk = (xi3_at_theta(i,j,k)-earth_radius)                &
                         /height_domain

          IF (eta_ijk>eta_s) THEN
            mu_w(i,j,k) = eg_vert_damp_coeff                          &
                         * SIN( 0.5*pi*(eta_ijk-eta_s)                &
                                *recip_one_m_etas)**2               &
                         * timestep

          ELSE
            mu_w(i,j,k) = 0.0
          END IF

        END DO
      END DO
    END DO

! Polar damping ------------------------------------------------------

  CASE(polar_damp_prof) 

    recip_one_m_etas = 1.0/(1.0-eta_s)

    DO k=wdims%k_start,wdims%k_end
      DO j=wdims%j_start,wdims%j_end
        DO i=wdims%i_start,wdims%i_end

          eta_ijk = (xi3_at_theta(i,j,k)-earth_radius)                &
                         /height_domain


          ! add polar dip

          cos_xi2 = COS(xi2_p(j))
          eta_ijk = cos_xi2*eta_ijk+(1.0-cos_xi2)

          IF (eta_ijk>eta_s) THEN
            mu_w(i,j,k) = eg_vert_damp_coeff                          &
                         * SIN( 0.5*pi*(eta_ijk-eta_s)                &
                                *recip_one_m_etas)**2               &
                         * timestep

          ELSE
            mu_w(i,j,k) = 0.0
          END IF

        END DO
      END DO
    END DO

! Polar damping 2 ----------------------------------------------------

  CASE(polar_damp_prof2) 

    recip_one_m_etas = 1.0/(1.0-eta_s)

    DO k=wdims%k_start,wdims%k_end
      DO j=wdims%j_start,wdims%j_end
        DO i=wdims%i_start,wdims%i_end

          eta_ijk = (xi3_at_theta(i,j,k)-earth_radius)                &
                         /height_domain


          ! add polar dip

          cos_xi2 = COS(xi2_p(j))
          eta_ijk = cos_xi2*eta_ijk+(1.0-cos_xi2)

          IF (eta_ijk>eta_s) THEN
            mu_w(i,j,k) = eg_vert_damp_coeff                          &
                         * (  SIN( 0.5*pi*(eta_ijk-eta_s)             &
                                  *recip_one_m_etas)**2             &
                            + SIN(xi2_p(j))**40.0                     &
                           )                                          &
                         * timestep

          ELSE
            mu_w(i,j,k) = 0.0
          END IF

        END DO
      END DO
    END DO

  CASE(polar_damp_prof2_NP) 

    recip_one_m_etas = 1./(1.-eta_s)

    DO k=wdims%k_start,wdims%k_end
      DO j=wdims%j_start,wdims%j_end
        DO i=wdims%i_start,wdims%i_end

          eta_ijk = (xi3_at_theta(i,j,k)-earth_radius)                &
                         /height_domain


          ! add polar dip (North pole only)

          IF (xi2_p(j)>0.0) THEN
            cos_xi2 = COS(xi2_p(j))
            eta_ijk = cos_xi2*eta_ijk+(1.-cos_xi2)

            IF (eta_ijk>eta_s) THEN
              mu_w(i,j,k) = eg_vert_damp_coeff                        &
                         * (  SIN( 0.5*pi*(eta_ijk-eta_s)             &
                                  *recip_one_m_etas)**2.              &
                            + SIN(xi2_p(j))**40.                      &
                           )                                          &
                         * timestep

            ELSE
              mu_w(i,j,k) = 0.
            END IF
          ELSE
            eta_ijk = eta_ijk

            IF (eta_ijk>eta_s) THEN
              mu_w(i,j,k) = eg_vert_damp_coeff                        &
                         * (  SIN( 0.5*pi*(eta_ijk-eta_s)             &
                                  *recip_one_m_etas)**2.              &
                           )                                          &
                         * timestep

            ELSE
              mu_w(i,j,k) = 0.
            END IF
          END IF

        END DO
      END DO
    END DO

  CASE(polar_damp_prof2_narrow) 

    recip_one_m_etas = 1.0/(1.0-eta_s)

    DO k=wdims%k_start,wdims%k_end
      DO j=wdims%j_start,wdims%j_end
        DO i=wdims%i_start,wdims%i_end

          eta_ijk = (xi3_at_theta(i,j,k)-earth_radius)                &
                         /height_domain


          ! add polar dip

          cos_xi2 = COS(xi2_p(j))
          eta_ijk = cos_xi2*eta_ijk+(1.0-cos_xi2)

          cutoff_deg = 87.0
          cutoff_rad = cutoff_deg/180.*pi

          cutoff = 0.0

          IF ( ABS(xi2_p(j)) > cutoff_rad ) THEN

            cutoff = 1.0-COS(ABS(xi2_p(j)-cutoff_rad)/                &
                     (0.5*pi-cutoff_rad)*0.5*pi)**20

          END IF

          IF (eta_ijk>eta_s) THEN
            mu_w(i,j,k) = eg_vert_damp_coeff                          &
                         * (  SIN( 0.5*pi*(eta_ijk-eta_s)             &
                                  *recip_one_m_etas)**2             &
                            + cutoff                                  &
                           )                                          &
                         * timestep

          ELSE
            mu_w(i,j,k) = 0.0
          END IF

        END DO
      END DO
    END DO

  CASE(polar_damp_prof2_extra_narrow) 

    recip_one_m_etas = 1.0/(1.0-eta_s)

    DO k=wdims%k_start,wdims%k_end
      DO j=wdims%j_start,wdims%j_end
        DO i=wdims%i_start,wdims%i_end

          eta_ijk = (xi3_at_theta(i,j,k)-earth_radius)                &
                         /height_domain


          ! add polar dip

          cos_xi2 = COS(xi2_p(j))**.46
          eta_ijk = cos_xi2*eta_ijk+(1.0-cos_xi2)

          cutoff_deg = 87.0
          cutoff_rad = cutoff_deg/180.*pi

          cutoff = 0.0

          IF ( ABS(xi2_p(j)) > cutoff_rad ) THEN

            cutoff = 1.0-COS(ABS(xi2_p(j)-cutoff_rad)/                &
                     (0.5*pi-cutoff_rad)*0.5*pi)**20

          END IF

          IF (eta_ijk>eta_s) THEN
            mu_w(i,j,k) = eg_vert_damp_coeff                          &
                         * (  SIN( 0.5*pi*(eta_ijk-eta_s)             &
                                  *recip_one_m_etas)**2             &
                            + cutoff                                  &
                           )                                          &
                         * timestep

          ELSE
            mu_w(i,j,k) = 0.0
          END IF

        END DO
      END DO
    END DO

  CASE(polar_damp_prof3) 

    recip_one_m_etas = 1.0/(1.0-eta_s)

    DO k=wdims%k_start,wdims%k_end
      DO j=wdims%j_start,wdims%j_end
        DO i=wdims%i_start,wdims%i_end

          eta_ijk = (xi3_at_theta(i,j,k)-earth_radius)                &
                         /height_domain

          IF (eta_ijk>eta_s) THEN
            mu_w(i,j,k) = eg_vert_damp_coeff                          &
                         * (  SIN( 0.5*pi*(eta_ijk-eta_s)             &
                                  *recip_one_m_etas)**2             &
                           )                                          &
                         * timestep

          ELSE
            mu_w(i,j,k) = 0.0
          END IF

          mu_w(i,j,k) = mu_w(i,j,k) + eg_vert_damp_coeff              &
                            * SIN(xi2_p(j))**40.0                     &
                            * timestep

        END DO
      END DO
    END DO

! Constant -----------------------------------------------------------

  CASE(const_damp_prof) 

    recip_one_m_etas = 1.0/(1.0-eta_s)

    DO k=wdims%k_start,wdims%k_end
      DO j=wdims%j_start,wdims%j_end
        DO i=wdims%i_start,wdims%i_end

           mu_w(i,j,k) = eg_vert_damp_coeff  * timestep

        END DO
      END DO
    END DO

  CASE default

    WRITE(message,*)  'vertical damping profile unknown'

    ErrorStatus = ABS (eg_vert_damp_profile)

    CALL ereport(RoutineName,ErrorStatus,message)

END SELECT


IF ( PrintStatus == PrStatus_Diag .AND. me == 0 ) THEN

  write(6,fmt='(A)') '================================================'
  write(6,fmt='(A)') 'level  | vertical damping coefficient at i=1,j=1'
  write(6,fmt='(A)') '------------------------------------------------'

  DO k=wdims%k_start,wdims%k_end

    write(6,fmt='(I4,A,E15.5)') k,'   |', mu_w(1,1,k)/timestep

  END DO
  write(6,fmt='(A)') '================================================'

END IF

IF (lhook) CALL dr_hook('INIT_VERT_DAMP',zhook_out,zhook_handle)
RETURN

END SUBROUTINE init_vert_damp

END MODULE init_vert_damp_mod
