! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE eg_precon_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_precon(Px,x,pre_type)


      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE atm_fields_bounds_mod
      USE Field_Types
      USE eg_parameters_mod
      USE tri_sor_mod
      USE eg_swap_bounds_mod
      USE eg_calc_ax_mod

      IMPLICIT NONE

!
! Description: Preconditioners used to accelerate the
!              convergence of the linear solver within
!              Helmholtz equation.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

! Array dimensions

      INTEGER,  INTENT(In)    :: pre_type

! Input vector x and output vector Px (i.e. x multiplied by the
! Pre/Post-conditioner).

      REAL,     INTENT(INOUT) ::                                      &
          Px(pdims_s%i_start:pdims_s%i_end,                           &
             pdims_s%j_start:pdims_s%j_end,                           &
             pdims_s%k_start:pdims_s%k_end),                          &
           x(pdims_s%i_start:pdims_s%i_end,                           &
             pdims_s%j_start:pdims_s%j_end,                           &
             pdims_s%k_start:pdims_s%k_end)


! Local variables

      INTEGER :: i, j, k, n

      REAL tmp(pdims_s%i_start:pdims_s%i_end,                         &
               pdims_s%j_start:pdims_s%j_end,                         &
               pdims_s%k_start:pdims_s%k_end)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('EG_PRECON',zhook_in,zhook_handle)


! Generic number used for pteconditioners.
! Set to 3 currently since lower or greater values aren't worth the 
! computational effort. Note if retained in the final version for the SOR
! routine this will probably stay hardwired for the opmizations (to minimize
! calls to swap_bounds (it's currently 2*Nterms we need to get it down to 1).

      Nterms   = 3

! This code contains all the preconditioners - We could split it
! up at some point to make it easier to read!

      SELECT CASE(pre_type)

         CASE(0)             ! Diagonal preconditioner

            Px = x

         CASE(2)               ! Jacobi Preconditioer
            Px = x 

            DO n = 1, Nterms

               CALL eg_calc_ax(tmp,Px)

               DO k = pdims%k_start, pdims%k_end
                  DO j = pdims%j_start, pdims%j_end
                     DO i = pdims%i_start, pdims%i_end
                        Px(i,j,k) = Px(i,j,k) + x(i,j,k) - tmp(i,j,k)
                     END DO
                  END DO
               END DO


            END DO

         CASE(4)             ! SOR type preconditioner

            CALL tri_sor(Px,x,Nterms)

      END SELECT

      IF (lhook) CALL dr_hook('EG_PRECON',zhook_out,zhook_handle)

      END SUBROUTINE eg_precon
      END MODULE
