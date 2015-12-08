! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_generate_grid_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_generate_grid(                                      &
                      first_constant_rho_level                        &
,                     a_ixsts,len_a_ixsts, a_spsts,len_a_spsts        &
!  Grid information
,                     grid_number, grid_flat, height_domain           &
,                     h_o, z_orog_print)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE integrity_mod

USE ereport_mod, ONLY : ereport
USE Field_Types
USE PrintStatus_mod
USE gflat_mod, ONLY: gflat_old, gflat_linear, gflat_linear1, gflat_quadratic
USE eg_swap_bounds_mod
USE atm_fields_bounds_mod

USE level_heights_mod, ONLY : r_theta_levels,r_rho_levels,r_at_u,     &
                              r_at_v,eta_theta_levels,eta_rho_levels, &
                              r_at_u_w,r_at_v_w

USE earth_constants_mod, ONLY : earth_radius
USE proc_info_mod,       ONLY : me, model_domain
USE domain_params,       ONLY : mt_lam

IMPLICIT NONE
!
! Description:
!  
!        Generate r_theta_levels, r_rho_levels given
!        eta_theta_levels, eta_rho_levels and height_domain
!        r_theta_levels(i,j,0) contains orography on input
!        h_o is used to generate reference grid for printing
!        Grid type ( grid_number) and flattening (grid_flat)
!        are required (user options)

! This is a simplified version of the New Dynamics routine.
! Two simple options are allowed: linear or quadratic flatenning.
! The vertical grid spacing (regular, quadratic) is stored in
! eta_theta_levels, eta_rho_levels.

!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


REAL :: height_domain , h_o

INTEGER :: grid_number, grid_flat

INTEGER :: first_constant_rho_level

REAL :: z_orog_print(tdims%k_start:tdims%k_end)

! Stash indexing arrays
INTEGER :: len_a_ixsts
INTEGER :: len_a_spsts
INTEGER :: a_ixsts(len_a_ixsts)     ! stash index array
REAL    :: a_spsts(len_a_spsts)     ! atmos stash array

! local variables

INTEGER :: i, j, k, fcrl, error_code

REAL :: z_ref_theta(pdims%k_start:pdims%k_end)
REAL :: z_ref_rho(pdims%k_start:pdims%k_end)

CHARACTER(LEN=2) :: digit2

TYPE (array_dims) u_w_dims,v_w_dims

! Description: COMDECK containing vertical grid types
!  for use in idealised problems
!
      INTEGER, PARAMETER :: vert_regular=1
      INTEGER, PARAMETER :: vert_quadratic_theta=21
      INTEGER, PARAMETER :: vert_bi_quadratic=22
      INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
      INTEGER, PARAMETER :: vert_schar=3
      INTEGER, PARAMETER :: vert_dwd=4
      INTEGER, PARAMETER :: vert_stretch_plus_regular=5
      INTEGER, PARAMETER :: vert_quad_stretch_thin=6
      INTEGER, PARAMETER :: vert_regular_thin=7
      INTEGER, PARAMETER :: vert_geometric_theta=8
      INTEGER, PARAMETER :: vert_dump=10
      INTEGER, PARAMETER :: vert_idl_um_grid=11


! Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_IDL_GENERATE_GRID',zhook_in,zhook_handle)


! ---------------------------------------------------------------------
! Section 0.  Initialise Data fields.
! ---------------------------------------------------------------------
!   fcrl shorthand for first_constant_rho_level
fcrl = first_constant_rho_level
! ---------------------------------------------------------------------
! Section 1. Set up vertical co-ordinate arrays relative to surface
!            On entry r_theta_levels(level0) contains orographic
!            height relative to Earth's surface
!            Add Earth-radius to  r_theta_levels and r_rho_levels
!            after section 3
! ---------------------------------------------------------------------
! height_domain is either set for problem or changed in vertical grid

!      If (me == 0) Then
!        Write (6,*) ' '
!        Write (6,*) ' VERTICAL GRID'
!      End If
IF(grid_number  ==  vert_dump) THEN
  IF(me  ==  0)THEN
    WRITE(6,fmt='(A,I5)')'grid_number = ',grid_number
    WRITE(6,fmt='(A)')'Vertical grid unchanged - as in input dump '
  END IF   !(me  ==  0)
  IF (lhook) CALL dr_hook('EG_IDL_GENERATE_GRID',zhook_out,zhook_handle)
  RETURN
END IF  !   grid_number  ==  vert_dump


z_orog_print(0) =  h_o
DO k = 1, tdims%k_end
  z_ref_theta(k)  = eta_theta_levels(k) * height_domain
  z_ref_rho(k)    = eta_rho_levels(k)   * height_domain
  z_orog_print(k) = z_ref_theta(k)
END DO


! Set stash level arrays to the new level heights (zero orography)
! This sets blev and bhlev correctly in the PP header.
! At present the C(k) values (brlev and bhrlev) are set to zero
! as there are various different options in this routine.
! To plot on true height above msl, output height field on model levs
! z(i,j,k) = zsea(k) + C(k)*zorog(i,j)

IF (grid_number  /=  vert_dump) THEN
  ! Set theta height at surface to 0.0
  a_spsts(a_ixsts(3)) = 0.0
  a_spsts(a_ixsts(4)) = 0.0
  DO k = 1, tdims%k_end
    ! Set rho levels
    a_spsts(a_ixsts(1) + k-1) = z_ref_rho(k)
    ! Set C(k) to zero at present
    a_spsts(a_ixsts(2) + k-1) = 0.0

    ! Set theta levels
    a_spsts(a_ixsts(3) + k) = z_ref_theta(k)
    ! Set C(k) to zero at present
    a_spsts(a_ixsts(4) + k) = 0.0
  END DO
END IF

IF ( grid_flat  ==  gflat_linear ) THEN

  IF (me  ==  0) THEN
    WRITE(6,fmt='(A)')'** Linear flattening of grid between surface '
    WRITE(6,fmt='(A,I5)')' and first_constant_rho_level grid_flat = ', &
                        grid_flat
  END IF   !(me  ==  0)

  IF ( grid_number  <=  vert_quadratic_theta .OR.                     &
                    grid_number == vert_geometric_theta ) THEN

    DO k = 1, tdims%k_end

      z_orog_print(k) =  z_ref_theta(k) +                             &
              h_o * (1.0 - eta_theta_levels(k))

      DO j = pdims_l%j_start,pdims_l%j_end
        DO i = pdims_l%i_start,pdims_l%i_end

          r_rho_levels(i,j,k) = z_ref_rho(k) +                        &
              r_theta_levels(i,j,0) * (1.0 - eta_rho_levels(k))
          r_theta_levels(i,j,k) = z_ref_theta(k)  +                   &
            r_theta_levels(i,j,0) * (1.0 - eta_theta_levels(k))

        END DO
      END DO

      DO j = udims_l%j_start,udims_l%j_end
        DO i = udims_l%i_start,udims_l%i_end

          r_at_u(i,j,k) = z_ref_rho(k) +                              &
              r_at_u_w(i,j,0) * (1.0 - eta_rho_levels(k))
          r_at_u_w(i,j,k) = z_ref_theta(k) +                          &
              r_at_u_w(i,j,0) * (1.0 - eta_theta_levels(k))

        END DO
      END DO

      DO j = vdims_l%j_start,vdims_l%j_end
        DO i = vdims_l%i_start,vdims_l%i_end

          r_at_v(i,j,k) = z_ref_rho(k) +                              &
              r_at_v_w(i,j,0) * (1.0 - eta_rho_levels(k))
          r_at_v_w(i,j,k) = z_ref_theta(k) +                          &
              r_at_v_w(i,j,0) * (1.0 - eta_theta_levels(k))

        END DO
      END DO

    END DO

  ELSE

    error_code = 333
    WRITE(digit2,'(i2)') grid_number

    CALL ereport("EG_IDL_GENERATE_GRID", error_code,                  &
         "** grid_number "//digit2//" is not supported in ENDGame **")

  END IF

ELSE IF ( grid_flat  ==  gflat_quadratic) THEN

  IF (me == 0) THEN
    WRITE (UNIT=6,FMT='(A47,A36,I2,A1)')                              &
          '   Quadratic flattening of grid over orography ',          &
          '   starting from surface. (grid_flat = ',grid_flat,')'
  END IF

  IF ( grid_number  <=  vert_quadratic_theta .OR.                     &
                    grid_number == vert_geometric_theta ) THEN

! From surface to first_constant_rho_level use quadratic relaxation
! for flattening

    DO k = 1, first_constant_rho_level - 1
      z_orog_print(k) =  z_ref_theta(k) +                             &
              h_o * (1.0 - eta_theta_levels(k) /                      &
              eta_rho_levels(first_constant_rho_level))**2

      DO j = pdims_l%j_start,pdims_l%j_end
        DO i = pdims_l%i_start,pdims_l%i_end

          r_rho_levels(i,j,k) = z_ref_rho(k)  +                       &
               r_theta_levels(i,j,0) * (1.0 - eta_rho_levels(k) /     &
              eta_rho_levels(first_constant_rho_level))**2
          r_theta_levels(i,j,k) = z_ref_theta(k)  +                   &
             r_theta_levels(i,j,0) * (1.0 - eta_theta_levels(k) /     &
              eta_rho_levels(first_constant_rho_level))**2

        END DO
      END DO

      DO j = udims_l%j_start,udims_l%j_end
        DO i = udims_l%i_start,udims_l%i_end

          r_at_u(i,j,k) = z_ref_rho(k)  +                             &
               r_at_u_w(i,j,0) * (1.0 - eta_rho_levels(k) /           &
              eta_rho_levels(first_constant_rho_level))**2
          r_at_u_w(i,j,k) = z_ref_theta(k)  +                         &
               r_at_u_w(i,j,0) * (1.0 - eta_theta_levels(k) /         &
              eta_rho_levels(first_constant_rho_level))**2

        END DO
      END DO

      DO j = vdims_l%j_start,vdims_l%j_end
        DO i = vdims_l%i_start,vdims_l%i_end

          r_at_v(i,j,k) = z_ref_rho(k)  +                             &
               r_at_v_w(i,j,0) * (1.0 - eta_rho_levels(k) /           &
              eta_rho_levels(first_constant_rho_level))**2
          r_at_v_w(i,j,k) = z_ref_theta(k)  +                         &
               r_at_v_w(i,j,0) * (1.0 - eta_theta_levels(k) /         &
              eta_rho_levels(first_constant_rho_level))**2

        END DO
      END DO
    END DO

! For constant levels set r to be a constant on the level
    DO k = first_constant_rho_level, tdims%k_end
      z_orog_print(k) =   z_ref_theta(k)

      DO j = pdims_l%j_start,pdims_l%j_end
        DO i = pdims_l%i_start,pdims_l%i_end
          r_theta_levels(i,j,k) = z_ref_theta(k)
          r_rho_levels(i,j,k)   = z_ref_rho(k)
        END DO
      END DO

      DO j = udims_l%j_start,udims_l%j_end
        DO i = udims_l%i_start,udims_l%i_end
          r_at_u(i,j,k) = z_ref_rho(k)
        END DO
      END DO

      DO j = vdims_l%j_start,vdims_l%j_end
        DO i = vdims_l%i_start,vdims_l%i_end
          r_at_v(i,j,k) = z_ref_rho(k)
        END DO
      END DO
    END DO

  ELSE

    error_code = 444
    WRITE(digit2,'(i2)') grid_number

    CALL ereport("EG_IDL_GENERATE_GRID", error_code,                  &
        "** grid_number "//digit2//" is not supported in ENDGame **")

  END IF

END IF

IF ( grid_number == vert_idl_um_grid ) THEN

  IF (me == 0) THEN
    WRITE (UNIT=6,FMT='(A)')' Idealised UM height gen smoothed grid'
  END IF

! Grid generation from setcona
! For constant levels set r to be a constant on the level

      DO k = first_constant_rho_level, tdims%k_end
        DO j = pdims_l%j_start,pdims_l%j_end
          DO i = pdims_l%i_start,pdims_l%i_end
            r_theta_levels(i,j,k) = z_ref_theta(k)
            r_rho_levels(i,j,k)   = z_ref_rho(k)
          END DO
        END DO

        DO j = udims_l%j_start,udims_l%j_end
          DO i = udims_l%i_start,udims_l%i_end
            r_at_u(i,j,k) = z_ref_rho(k) 
            r_at_u_w(i,j,k) = z_ref_theta(k) 
          END DO
        END DO

        DO j = vdims_l%j_start,vdims_l%j_end
          DO i = vdims_l%i_start,vdims_l%i_end
            r_at_v(i,j,k) = z_ref_rho(k) 
            r_at_v_w(i,j,k) = z_ref_theta(k) 
          END DO
        END DO
      END DO

! A smooth quadratic height generation
      DO k = 1, first_constant_rho_level-1

        DO j = pdims_l%j_start,pdims_l%j_end
          DO i = pdims_l%i_start,pdims_l%i_end

              r_rho_levels(i,j,k) = eta_rho_levels(k) *            &
                height_domain +                                    &
     &          r_theta_levels(i,j,0)* (1.0 - eta_rho_levels(k)    &
     &              /eta_rho_levels(first_constant_rho_level))**2

              r_theta_levels(i,j,k) = eta_theta_levels(k) *        &
     &             height_domain + r_theta_levels(i,j,0) *         &
     &             (1.0 - eta_theta_levels(k) /                    &
     &              eta_rho_levels(first_constant_rho_level))**2

          END DO
        END DO

        DO j = udims_l%j_start,udims_l%j_end
          DO i = udims_l%i_start,udims_l%i_end

              r_at_u(i,j,k) = eta_rho_levels(k) *                  &
                height_domain +                                    &
     &          r_at_u_w(i,j,0)* (1.0 - eta_rho_levels(k)          &
     &              /eta_rho_levels(first_constant_rho_level))**2

              r_at_u_w(i,j,k) = eta_theta_levels(k) *              &
     &             height_domain + r_at_u_w(i,j,0) *               &
     &             (1.0 - eta_theta_levels(k) /                    &
     &              eta_rho_levels(first_constant_rho_level))**2

          END DO
        END DO

        DO j = vdims_l%j_start,vdims_l%j_end
          DO i = vdims_l%i_start,vdims_l%i_end

            r_at_v(i,j,k) = eta_rho_levels(k) *                    &
                height_domain +                                    &
     &          r_at_v_w(i,j,0)* (1.0 - eta_rho_levels(k)          &
     &              /eta_rho_levels(first_constant_rho_level))**2
            r_at_v_w(i,j,k) = eta_theta_levels(k) *                &
     &             height_domain + r_at_v_w(i,j,0) *               &
     &             (1.0 - eta_theta_levels(k) /                    &
     &              eta_rho_levels(first_constant_rho_level))**2

          END DO
        END DO
      END DO

END IF       !  on grid_flat

! ---------------------------------------------------------------------
! Section 2. Add Earth-radius to  r_theta_levels and r_rho_levels
!            calculate r_at_u points and r_at_v points on rho levels
! ---------------------------------------------------------------------


r_theta_levels(:,:,0) = r_theta_levels(:,:,0) + earth_radius
r_at_u_w(:,:,0)   = r_at_u_w(:,:,0) + earth_radius
r_at_v_w(:,:,0)   = r_at_v_w(:,:,0) + earth_radius


DO k = 1, tdims%k_end
  r_theta_levels(:,:,k) = r_theta_levels(:,:,k) + earth_radius
  r_rho_levels(:,:,k)   = r_rho_levels(:,:,k)   + earth_radius
  r_at_u(:,:,k)         = r_at_u(:,:,k)         + earth_radius
  r_at_u_w(:,:,k)       = r_at_u_w(:,:,k)       + earth_radius
  r_at_v(:,:,k)         = r_at_v(:,:,k)         + earth_radius
  r_at_v_w(:,:,k)       = r_at_v_w(:,:,k)       + earth_radius
END DO


IF ( PrintStatus == PrStatus_Diag .AND. me == 0 ) THEN
  write(6,fmt='(A)') '================================================'
  write(6,fmt='(A)') 'level  | r_theta_levls,r_rho_levls,etc @ i=1,j=1'
  write(6,fmt='(A)') '------------------------------------------------'

  write(6,fmt='(I4,A,E15.5)') 0,'   |', r_theta_levels(1,1,0)-earth_radius

  DO k=1,tdims%k_end

    write(6,fmt='(I4,A,6E15.5)') k,'   |',                            &
                              r_theta_levels(1,1,k)-earth_radius,     &
                              r_rho_levels  (1,1,k)-earth_radius,     &
                              r_at_u        (1,1,k)-earth_radius,     &
                              r_at_v        (1,1,k)-earth_radius,     &
                              r_at_u_w      (1,1,k)-earth_radius,     &
                              r_at_v_w      (1,1,k)-earth_radius

  END DO
  write(6,fmt='(A)') '================================================'
END IF

CALL eg_swap_bounds(r_theta_levels,tdims_l,fld_type_p,.FALSE.)
CALL eg_swap_bounds(r_rho_levels,pdims_l,fld_type_p,.FALSE.)
CALL eg_swap_bounds(r_at_u,udims_l,fld_type_u,.FALSE.)
CALL eg_swap_bounds(r_at_v,vdims_l,fld_type_v,.FALSE.)

u_w_dims = udims_l
u_w_dims%k_start = 0
CALL eg_swap_bounds(r_at_u_w,u_w_dims,fld_type_u,.FALSE.)

v_w_dims = vdims_l
v_w_dims%k_start = 0
CALL eg_swap_bounds(r_at_v_w,v_w_dims,fld_type_v,.FALSE.)


IF (integrity_test)                                                   &
  CALL update_hash_m(                                                 &
                  eta_theta_levels,SIZE(eta_theta_levels),'etatl',    &
                  eta_rho_levels,  SIZE(eta_rho_levels),  'etarl')

IF (lhook) CALL dr_hook('EG_IDL_GENERATE_GRID',zhook_out,zhook_handle)

END SUBROUTINE eg_idl_generate_grid
END MODULE eg_idl_generate_grid_mod
