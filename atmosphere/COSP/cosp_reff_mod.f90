! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE cosp_reff_mod
  USE conversions_mod, ONLY: zerodegc
  USE mphys_ice_mod, ONLY: t_agg_min
  USE ereport_mod
  USE gammaf_mod, ONLY: gammaf
  IMPLICIT NONE

! Description:
!   Routine that computes the hydrometeor effective radius of the precipitating
!   hydrometeors.
!
! Method:
!   Uses the analytical solution of the ith moment of the PSD to compute the
!   ratio between the 3rd and 2nd moments, divided by two. This is explained
!   in the COSP user's manual. It used the parameters that describe the PSD
!   in the UM microphysical settings (UMDP26).
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CONTAINS
  SUBROUTINE cosp_reff(flux,X1,X2,X3,X4,A,B,C,D,G, &
                npoints,model_levels,T,rho,mr,Reff)

  IMPLICIT NONE
!----Input arguments
! Precipitation flux is used
  LOGICAL,INTENT(IN) :: flux
! These are the constants that define the PSD in the microphysics 
! scheme (UMDP 26).
! X1,X3: n_ax = X1*exp(-X3*T[degC])
! X2: n_bx
! X4: alpha_x
  REAL,INTENT(IN) :: X1
  REAL,INTENT(IN) :: X2
  REAL,INTENT(IN) :: X3
  REAL,INTENT(IN) :: X4
! These are the constants that define the mass-diameter relationship.
! M_x(D) = a_x*D^(b_x)
! A: a_x
! B: b_x
  REAL,INTENT(IN) :: A
  REAL,INTENT(IN) :: B
! These are the constants that define the terminal fall speed relationship.
! V_x(D) = c_x*D^(d_x)*exp(-h_x)*(rho_0/rho)^g_x
! C: c_x
! D: d_x
! G: g_x
! h_x is always 0
! Abel and Shipway (2007) is not supported
  REAL,INTENT(IN) :: C
  REAL,INTENT(IN) :: D
  REAL,INTENT(IN) :: G
! Dimensions
  INTEGER,INTENT(IN) :: npoints
  INTEGER,INTENT(IN) :: model_levels
! Air temperature [K]
  REAL,INTENT(IN) :: T(npoints,model_levels)
! Air density [kg/m^3]
  REAL,INTENT(IN) :: rho(npoints,model_levels)
! Hydrometeor mixing ratio [kg/kg]
  REAL,INTENT(IN) :: mr(npoints,model_levels)

!----Output arguments
! Effective radius [m]
  REAL,INTENT(OUT) :: Reff(npoints,model_levels)

!----Local variables
  REAL,PARAMETER :: rho_0 = 1.0
  REAL :: gamma_a3,gamma_a4,gamma_ab1,gamma_abd1,frac_exp,gamma_ratio
  REAL :: F_nax(npoints,model_levels)
  CHARACTER(LEN=9) :: routine_name='COSP_REFF'
  INTEGER :: icode,i,k

  icode = 9

  CALL gammaf(3.0+X4,gamma_a3)
  CALL gammaf(4.0+X4,gamma_a4)
  CALL gammaf(1.0+X4+B,gamma_ab1)
  CALL gammaf(1.0+X4+B+D,gamma_abd1)


  Reff = 0.0

  IF (A <= 0.0) CALL Ereport(routine_name, icode," A <= 0.0")

! Compute intercept as function of T, if needed
  IF (X3 /= 0.0) F_nax = EXP(-X3*MAX(T-ZeroDegC,T_AGG_MIN))

! Compute the parameter lambda^-1 of the PSD. stored in variable Reff
  IF (flux) THEN ! precipitation flux. Fall speed needed
      frac_exp = 1.0/(X4+D-X2+4.0)
      IF (C <= 0.0) CALL Ereport(routine_name, icode," C <= 0.0")
      IF (X3 /= 0.0) THEN
            Reff = (mr/(A*C*((rho_0/rho)**G)*X1*F_nax*gamma_abd1))**frac_exp
      ELSE
            Reff = (mr/(A*C*((rho_0/rho)**G)*X1*gamma_abd1))**frac_exp
      END IF
  ELSE ! mixing ratio
      frac_exp = 1.0/(X4+B-X2+1.0)
      IF (X3 /= 0.0) THEN
          Reff = ((rho*mr)/(A*gamma_ab1*X1*F_nax))**frac_exp
      ELSE
          Reff = ((rho*mr)/(A*gamma_ab1*X1))**frac_exp
      END IF
  END IF

! Compute radius and apply sanity check
  gamma_ratio = 0.5*(gamma_a4/gamma_a3)
  DO k = 1, model_levels
    DO i = 1, npoints
      IF (Reff(i,k) > 0.0) Reff(i,k) = gamma_ratio*Reff(i,k)
      IF (Reff(i,k) < 0.0) Reff(i,k) = 0.0
    END DO
  END DO

  RETURN
  END SUBROUTINE cosp_reff
END MODULE cosp_reff_mod
