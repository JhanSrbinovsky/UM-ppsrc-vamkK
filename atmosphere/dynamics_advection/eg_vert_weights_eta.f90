! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
SUBROUTINE eg_vert_weights_eta (                                      &
                        dim_k_in, dim_i_out, dim_j_out, dim_k_out,    &
                        high_order_scheme, monotone_scheme,           &
                        model_domain, l_high, k_int_linear,           &
                        k_out, eta_in, eta_out,                       &
                        coeff_z )

USE um_types,   ONLY: integer32
USE parkind1,   ONLY: jpim, jprb       !DrHook
USE yomhook,    ONLY: lhook, dr_hook   !DrHook
USE highos_mod, ONLY: cubicLagrange, quinticLagrange, hCubic_vLin,    &
                      hQuasiCubic_vQuintic, hCubic_vQuintic,          &
                      hLag3_vHerm3_d2, hLag3_vHerm3_d4

IMPLICIT NONE
!
! Description:Calculates vertical interpolation weights.
!  
!
! Method:
! This subroutine calculates vertical interpolation weights
! appropriate to doing interpolation in eta coordinate.
! Eta coordinates for departure points must be given as inputs.
! Subroutine was developed extracting/adapting code from New Dynamics
! eta_vert_weights
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


INTEGER, INTENT(IN) ::                                                &
  dim_k_in,          &  ! Dimension of Data_in in k direction.
  dim_i_out,         &  ! Dimension of Data_out in i direction.
  dim_j_out,         &  ! Dimension of Data_out in j direction.
  dim_k_out,         &  ! Dimension of Data_out in k direction.
  high_order_scheme, &  ! a code saying which high order scheme to use.    
  monotone_scheme,   &  ! a code saying which monotone scheme to use.
  model_domain          ! holds integer code for model domain

INTEGER, INTENT(IN) :: k_int_linear
                        ! Level below which linear interpolation is used.

LOGICAL, INTENT(IN) :: l_high
                        ! True, if high order interpolation required.

INTEGER (KIND=integer32), INTENT(IN) ::                               &
                               k_out (dim_i_out, dim_j_out, dim_k_out)

REAL, INTENT(IN) ::                                                   &
  eta_in(dim_k_in),                                                   &
                 ! Vertical co-ordinate of input data.
  eta_out(dim_i_out, dim_j_out, dim_k_out)
                 ! Vertical co-ordinate of output data.


REAL, INTENT(OUT) ::                                                  &
  coeff_z(-2:3, dim_i_out, dim_j_out, dim_k_out)

! Local Variables.

! Loop indices
INTEGER :: i, j, k

REAL :: rq3(-1:2,dim_k_in), rq5(-2:3,dim_k_in), rdeta(dim_k_in)

REAL :: e_here_minus2, e_here_minus, e_here, e_here_plus,             &
        e_here_plus2, e_here_plus3, numer_minus, numer,               &
        numer_plus, numer_plus2, numer_minus2, numer_plus3,           &
        weight_eta

! Variables required for hermite cubic
INTEGER(KIND=integer32) :: kD ! Index of grid level immediately below departure
                              ! point

REAL    :: xi ! Non-dimensional coordinate of departure point on interval
              ! [z(kD),z(kD+1)]

REAL    :: deta(dim_k_in-1)    ! Mesh
REAL    :: H1, H2, H3, H4      ! Basis functions for Hermite cubic
REAL    :: rmesh(2:dim_k_in-1) ! ratio of mesh spacings, used for quadratic
                               ! derivative estimate

REAL    :: omega_prime(-2:2,3:dim_k_in-2) ! Derivatives of quartic Lagrange
                                          ! basis functions, evaluated at
                                          ! central point.

! Description: COMDECK containing the allowed
!              monotone scheme options
!
      INTEGER                                                           &
     &     triLinear                                                    &
     &,    mono_quasiCubic

      PARAMETER(                                                        &
     &     triLinear       = 1                                          &
     &,    mono_quasiCubic = 2 )

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_VERT_WEIGHTS_ETA',zhook_in,zhook_handle)


DO k = 1, dim_k_in-1
   deta(k)  = eta_in(k+1) - eta_in(k)
   rdeta(k) = 1.0/deta(k) 
END DO

! ----------------------------------------------------------------------
! Section 1.   Calculate quintic vertical interpolation coefficients
! ----------------------------------------------------------------------

IF ( high_order_scheme == quinticlagrange      .OR.                 &
     high_order_scheme == hquasicubic_vquintic .OR.                 &
     high_order_scheme == hcubic_vquintic          ) THEN

   CALL vquintic_consts()

!$OMP PARALLEL DO PRIVATE(i,j,weight_eta, e_here_minus2,e_here_minus, &
!$OMP&  e_here,e_here_plus, e_here_plus2, e_here_plus3,               &
!$OMP&  numer_minus2,numer_minus,numer,numer_plus,numer_plus2,        &
!$OMP&  numer_plus3)  SCHEDULE(STATIC)
  DO k = 1, dim_k_out
    DO j = 1, dim_j_out
      DO i = 1, dim_i_out

        IF ( k_out(i,j,k) <= k_int_linear .OR.                        &
             k_out(i,j,k) == dim_k_in-1 ) THEN

          weight_eta = (eta_out(i,j,k) - eta_in(k_out(i,j,k)))        &
                      *rdeta(k_out(i,j,k)) 

          coeff_z(-2,i,j,k) = 0.
          coeff_z(-1,i,j,k) = 0.
          coeff_z(0,i,j,k)  = 1.0 - weight_eta
          coeff_z(1,i,j,k)  = weight_eta
          coeff_z(2,i,j,k)  = 0.
          coeff_z(3,i,j,k)  = 0.

        ELSE IF (k_out(i,j,k) == 2 .OR.                               &
                 k_out(i,j,k) == dim_k_in-2) THEN

! use cubic interpolation.
          e_here_minus = eta_in(k_out(i,j,k) - 1)
          e_here       = eta_in(k_out(i,j,k))
          e_here_plus  = eta_in(k_out(i,j,k) + 1)
          e_here_plus2 = eta_in(k_out(i,j,k) + 2)

! Compute coordinate differences in eta.
          numer_minus = eta_out(i,j,k) - e_here_minus
          numer       = eta_out(i,j,k) - e_here
          numer_plus  = eta_out(i,j,k) - e_here_plus
          numer_plus2 = eta_out(i,j,k) - e_here_plus2

          coeff_z(-2,i,j,k)= 0.

          coeff_z(-1,i,j,k)= (numer*numer_plus*numer_plus2)           &
                             *rq5(-1,k_out(i,j,k))

          coeff_z(0,i,j,k) = (numer_minus*numer_plus*numer_plus2 )    &
                            *rq5(0,k_out(i,j,k))

          coeff_z(1,i,j,k) = (numer_minus*numer*numer_plus2 )         &
                            *rq5(1,k_out(i,j,k))

          coeff_z(2,i,j,k) = (numer_minus*numer*numer_plus )          &
                            *rq5(2,k_out(i,j,k))

          coeff_z(3,i,j,k)= 0.

        ELSE
! quintic interpolation

! Re-use the r coordinate variables to hold the corresponding
!     eta coordinate values.

          e_here_minus2 = eta_in(k_out(i,j,k) - 2)
          e_here_minus  = eta_in(k_out(i,j,k) - 1)
          e_here        = eta_in(k_out(i,j,k))
          e_here_plus   = eta_in(k_out(i,j,k) + 1)
          e_here_plus2  = eta_in(k_out(i,j,k) + 2)
          e_here_plus3  = eta_in(k_out(i,j,k) + 3)

!  Compute coordinate differences in eta.
          numer_minus2 = eta_out(i,j,k) - e_here_minus2
          numer_minus  = eta_out(i,j,k) - e_here_minus
          numer        = eta_out(i,j,k) - e_here
          numer_plus   = eta_out(i,j,k) - e_here_plus
          numer_plus2  = eta_out(i,j,k) - e_here_plus2
          numer_plus3  = eta_out(i,j,k) - e_here_plus3

          coeff_z(-2,i,j,k)= ( numer_minus*numer*numer_plus           &
                              *numer_plus2 * numer_plus3 )            &
                            *rq5(-2,k_out(i,j,k))

          coeff_z(-1,i,j,k)= ( numer_minus2*numer*numer_plus          &
                              *numer_plus2 * numer_plus3 )            &
                            *rq5(-1,k_out(i,j,k))

          coeff_z(0,i,j,k) = ( numer_minus2*numer_minus*numer_plus    &
                              *numer_plus2*numer_plus3 )              &
                            *rq5(0,k_out(i,j,k))

          coeff_z(1,i,j,k) = ( numer_minus2*numer_minus*numer         &
                              *numer_plus2*numer_plus3 )              &
                            *rq5(1,k_out(i,j,k))

          coeff_z(2,i,j,k) = ( numer_minus2*numer_minus*numer         &
                              *numer_plus*numer_plus3 )               &
                            *rq5(2,k_out(i,j,k))

          coeff_z(3,i,j,k)= ( numer_minus2*numer_minus*numer          &
                             *numer_plus*numer_plus2 )                &
                            *rq5(3,k_out(i,j,k))
        END IF

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! end quintic in vertical

! ---------------------------------------------------------------------
! Section 2.   Calculate cubic vertical interpolation coefficients
! ---------------------------------------------------------------------

IF ( high_order_scheme == cubiclagrange ) THEN

   CALL vcubic_consts()

!$OMP PARALLEL DO PRIVATE(i,j,weight_eta, e_here_minus2,e_here_minus, &
!$OMP&  e_here,e_here_plus, e_here_plus2,                             &
!$OMP&  numer_minus2,numer_minus,numer,numer_plus,numer_plus2)        &
!$OMP&  SCHEDULE(STATIC)
  DO k = 1, dim_k_out
    DO j = 1, dim_j_out
      DO i = 1, dim_i_out

        IF ( k_out(i,j,k) <= k_int_linear .OR.                        &
             k_out(i,j,k) == dim_k_in-1 ) THEN

! Obtain an eta coordinate for the departure point.
          weight_eta = (eta_out(i,j,k) - eta_in(k_out(i,j,k)))        &
                      *rdeta(k_out(i,j,k))

          coeff_z(-1,i,j,k) = 0.

!  Compute coeffs using eta coord system.

!  Linear interpolation weight already calculated.
          coeff_z(0,i,j,k) = 1.0 - weight_eta
          coeff_z(1,i,j,k) = weight_eta

          coeff_z (2,i,j,k) = 0.

        ELSE ! If not near boundary.
! use cubic interpolation.
! the knots for vertical interpolation are the tabulated eta levels
          e_here_minus = eta_in(k_out(i,j,k) - 1)
          e_here       = eta_in(k_out(i,j,k))
          e_here_plus  = eta_in(k_out(i,j,k) + 1)
          e_here_plus2 = eta_in(k_out(i,j,k) + 2)

! Compute coordinate differences in eta.
          numer_minus = eta_out(i,j,k) - e_here_minus
          numer       = eta_out(i,j,k) - e_here
          numer_plus  = eta_out(i,j,k) - e_here_plus
          numer_plus2 = eta_out(i,j,k) - e_here_plus2

          coeff_z(-1,i,j,k)= ( numer*numer_plus*numer_plus2)           &
                            *rq3(-1,k_out(i,j,k))

          coeff_z(0,i,j,k) = ( numer_minus*numer_plus*numer_plus2 )    &
                             *rq3(0,k_out(i,j,k))

          coeff_z(1,i,j,k) = ( numer_minus*numer*numer_plus2 )         &
                             *rq3(1,k_out(i,j,k))

          coeff_z(2,i,j,k) = ( numer_minus*numer*numer_plus )          &
                             *rq3(2,k_out(i,j,k))

        END IF ! Check on being close to boundary.

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! end cubic in vertical

! ---------------------------------------------------------------------
! Section 3.   Calculate linear vertical interpolation coefficients
! ---------------------------------------------------------------------

IF ( high_order_scheme == hcubic_vlin .OR. ( .NOT. l_high )) THEN

!$OMP PARALLEL DO PRIVATE(i,j,weight_eta)
  DO k=1, dim_k_out
    DO j = 1, dim_j_out
      DO i = 1, dim_i_out

        weight_eta = (eta_out(i,j,k) - eta_in(k_out(i,j,k)))         &
                    *rdeta(k_out(i,j,k))

        coeff_z(0,i,j,k) = 1.0 - weight_eta
        coeff_z(1,i,j,k) = weight_eta

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

! ---------------------------------------------------------------------
! Section 4.   Calculate Hermite cubic with quadratic derivatives
! ---------------------------------------------------------------------

IF ( high_order_scheme == hLag3_vHerm3_d2) THEN

  DO k=2, dim_k_in-1
    rmesh(k)=deta(k)*rdeta(k-1)
  END DO

!$OMP PARALLEL DO PRIVATE(i,j,k,xi,kD,H1,H2,H3,H4)  SCHEDULE(STATIC)
  DO k = 1, dim_k_out
    DO j = 1, dim_j_out
      DO i = 1, dim_i_out

        kD = k_out(i,j,k)
        xi = (eta_out(i,j,k)-eta_in(kD)) / (eta_in(kD+1)-eta_in(kD))

        H1=(1.0+2.0*xi)*(1.0-xi)**2
        H2=xi**2*(3.0-2.0*xi)
        H3=xi*(1.0-xi)**2
        H4=-xi**2*(1.0-xi)

        IF ( kD <= k_int_linear .OR. kD == dim_k_in-1 ) THEN

          coeff_z(-1,i,j,k) = 0.0
          coeff_z(0,i,j,k)  = 1.0-xi
          coeff_z(1,i,j,k)  = xi
          coeff_z(2,i,j,k)  = 0.0

        ELSE
            
          coeff_z(-1,i,j,k)=-H3*rmesh(kD)**2/(1.0+rmesh(kD))

          coeff_z(0,i,j,k) = H1 + H3*(rmesh(kD)-1.0) -                &
                             H4*rmesh(kD+1)/(1.0+rmesh(kD+1))

          coeff_z(1,i,j,k) = H2 + H3/(1.0+rmesh(kD)) +                &
                             H4*(rmesh(kD+1)-1.0)/rmesh(kD+1)

          coeff_z(2,i,j,k) = H4/(rmesh(kD+1)*(1.0+rmesh(kD+1)))

        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

! ---------------------------------------------------------------------
! Section 6.   Calculate Hermite cubic with quartic derivatives
! ---------------------------------------------------------------------

IF ( high_order_scheme == hLag3_vHerm3_d4) THEN

  DO k=2, dim_k_in-1
    rmesh(k)=deta(k)*rdeta(k-1)
  END DO

  CALL basis_derivs_for_hermite()

!$OMP PARALLEL DO PRIVATE(i,j,k,xi,kD,H1,H2,H3,H4)  SCHEDULE(STATIC)
  DO k = 1, dim_k_out
    DO j = 1, dim_j_out
      DO i = 1, dim_i_out

        kD = k_out(i,j,k)
        xi = (eta_out(i,j,k)-eta_in(kD)) / (eta_in(kD+1)-eta_in(kD))

        H1=(1.0+2.0*xi)*(1.0-xi)**2
        H2=xi**2*(3.0-2.0*xi)
        H3=xi*(1.0-xi)**2
        H4=-xi**2*(1.0-xi)

        coeff_z(:,i,j,k) = 0.0

        IF ( kD <= k_int_linear .OR. kD == dim_k_in-1 ) THEN

          coeff_z(-1,i,j,k) = 0.0
          coeff_z(0,i,j,k)  = 1.0-xi
          coeff_z(1,i,j,k)  = xi
          coeff_z(2,i,j,k)  = 0.0

        ELSE IF (kD == 2 .OR. kD == dim_k_in-2) THEN
          
          coeff_z(-1,i,j,k)=-H3*rmesh(kD)**2/(1.0+rmesh(kD))

          coeff_z(0,i,j,k) = H1 + H3*(rmesh(kD)-1.0)                           &
                                - H4*rmesh(kD+1)/(1.0+rmesh(kD+1))

          coeff_z(1,i,j,k) = H2 + H3/(1.0+rmesh(kD))                           &
                                + H4*(rmesh(kD+1)-1.0)/rmesh(kD+1)

          coeff_z(2,i,j,k) = H4/(rmesh(kD+1)*(1.0+rmesh(kD+1)))

        ELSE

          coeff_z(-2,i,j,k) = deta(kD)*H3*omega_prime(-2,kD)

          coeff_z(-1,i,j,k) = deta(kD)*( H3*omega_prime(-1,kD) +               &
                                         H4*omega_prime(-2,kD+1) )

          coeff_z(0,i,j,k) = H1 + deta(kD)*( H3*omega_prime(0,kD) +            &
                                             H4*omega_prime(-1,kD+1) )

          coeff_z(1,i,j,k) = H2 + deta(kD)*( H3*omega_prime(1,kD) +            &
                                             H4*omega_prime(0,kD+1) )

          coeff_z(2,i,j,k) = deta(kD)*( H3*omega_prime(2,kD) +                 &
                                        H4*omega_prime(1,kD+1) )

          coeff_z(3,i,j,k) = deta(kD)*H4*omega_prime(2,kD+1)

        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

! End of routine.

IF (lhook) CALL dr_hook('EG_VERT_WEIGHTS_ETA',zhook_out,zhook_handle)

CONTAINS

   SUBROUTINE basis_derivs_for_hermite()
   IMPLICIT NONE
   INTEGER :: i, j, k

   omega_prime=0.0

   DO k=3,dim_k_in-2
     DO j=-2,2
       IF ( j == 0 ) THEN
         omega_prime(j,k)=0.0
         DO i=1,2
           omega_prime(j,k)=omega_prime(j,k)+1.0/(eta_in(k)-eta_in(k-i))       &
                                            +1.0/(eta_in(k)-eta_in(k+i))
         END DO
       ELSE
         omega_prime(j,k)=1.0/(eta_in(k+j)-eta_in(k))
         DO i=-2,2
           IF (i /= 0 .AND. i /= j) THEN
             omega_prime(j,k)=omega_prime(j,k)*(eta_in(k)-eta_in(k+i)) /       &
                                               (eta_in(k+j)-eta_in(k+i))
           END IF
         END DO
       END IF
     END DO
   END DO

   END SUBROUTINE basis_derivs_for_hermite


   SUBROUTINE vquintic_consts()
   IMPLICIT NONE
   INTEGER  :: k
   REAL     :: p(-2:3), tmp

   rq5 = -1.0e20

   DO k = 3, dim_k_in-3
      p(-2)      = eta_in(k-2)
      p(-1)      = eta_in(k-1)
      p(0)       = eta_in(k)
      p(1)       = eta_in(k+1)
      p(2)       = eta_in(k+2)
      p(3)       = eta_in(k+3)

      tmp        = (p(-2) - p(-1))*(p(-2) - p(0))*(p(-2) - p(1))             &
                  *(p(-2) - p(2))*(p(-2) - p(3))
      rq5(-2,k)  = 1.0/tmp

      tmp        = (p(-1) - p(-2))*(p(-1) - p(0))*(p(-1) - p(1))             &
                  *(p(-1) - p(2))*(p(-1) - p(3))
      rq5(-1,k) = 1.0/tmp

      tmp       = (p(0) - p(-2))*(p(0) - p(-1))*(p(0) - p(1))               &
                 *(p(0) - p(2))*(p(0) - p(3))
      rq5(0,k)  = 1.0/tmp

      tmp       = (p(1) - p(-2))*(p(1) - p(-1))*(p(1) - p(0))               &
                 *(p(1) - p(2))*(p(1) - p(3))
      rq5(1,k)  = 1.0/tmp

      tmp       = (p(2) - p(-2))*(p(2) - p(-1))*(p(2) - p(0))               &
                 *(p(2) - p(1))*(p(2) - p(3))
      rq5(2,k)  = 1.0/tmp

      tmp       = (p(3) - p(-2))*(p(3) - p(-1))*(p(3) - p(0))               &
                 *(p(3) - p(1))*(p(3) - p(2))
      rq5(3,k) = 1.0/tmp
   END DO

   k = 2
   p(-1)     = eta_in(k-1)
   p(0)      = eta_in(k)
   p(1)      = eta_in(k+1)
   p(2)      = eta_in(k+2)

   tmp       = (p(-1) - p(0))*(p(-1) - p(1))*(p(-1) - p(2))
   rq5(-1,k) = 1.0/tmp

   tmp       = (p(0) - p(-1))*(p(0) - p(1))*(p(0) - p(2))
   rq5(0,k)   = 1.0/tmp

   tmp       = (p(1) - p(-1))*(p(1) - p(0))*(p(1) - p(2))
   rq5(1,k)   = 1.0/tmp

   tmp       = (p(2) - p(-1))*(p(2) - p(0))*(p(2) - p(1))
   rq5(2,k)   = 1.0/tmp

   k = dim_k_in-2
   p(-1)     = eta_in(k-1)
   p(0)      = eta_in(k)
   p(1)      = eta_in(k+1)
   p(2)      = eta_in(k+2)

   tmp       = (p(-1) - p(0))*(p(-1) - p(1))*(p(-1) - p(2))
   rq5(-1,k) = 1.0/tmp

   tmp       = (p(0) - p(-1))*(p(0) - p(1))*(p(0) - p(2))
   rq5(0,k)  = 1.0/tmp

   tmp       = (p(1) - p(-1))*(p(1) - p(0))*(p(1) - p(2))
   rq5(1,k)  = 1.0/tmp

   tmp       = (p(2) - p(-1))*(p(2) - p(0))*(p(2) - p(1))
   rq5(2,k)  = 1.0/tmp

   END SUBROUTINE vquintic_consts


   SUBROUTINE vcubic_consts
   IMPLICIT NONE
   INTEGER  :: k
   REAL     :: p(-1:2), tmp

   rq3 = -1.0e20
   
   DO k = 2, dim_k_in-2
      p(-1)      = eta_in(k-1)
      p(0)       = eta_in(k)
      p(1)       = eta_in(k+1)
      p(2)       = eta_in(k+2)

      tmp       = (p(-1) - p(0))*(p(-1) - p(1))*(p(-1) - p(2))
      rq3(-1,k) = 1.0/tmp

      tmp       = (p(0) - p(-1))*(p(0) - p(1))*(p(0) - p(2))
      rq3(0,k)  = 1.0/tmp

      tmp       = (p(1) - p(-1))*(p(1) - p(0))*(p(1) - p(2))
      rq3(1,k)  = 1.0/tmp

      tmp       = (p(2) - p(-1))*(p(2) - p(0))*(p(2) - p(1))
      rq3(2,k)  = 1.0/tmp
   END DO

   END SUBROUTINE vcubic_consts

END SUBROUTINE eg_vert_weights_eta
