! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ PC2 Cloud Scheme: Inhomogeneous forcing
! Subroutine Interface:
SUBROUTINE pc2_inhomogeneous(                                           &
!      Prognostic Fields
 cf, cfl, cff, qcl, qcf,                                                &
!      Forcing quantities for driving the inhomogeneous forcing
 q4_l, q4_f, ls, g_l, g_f)

  USE yomhook,              ONLY: lhook, dr_hook
  USE parkind1,             ONLY: jprb, jpim
  USE atm_fields_bounds_mod,ONLY: qdims

  IMPLICIT NONE

! Purpose:
!   This subroutine calculates the change in liquid cloud fraction,
!   ice cloud fraction and total cloud fraction as a result of
!   inhomogeneous forcing of the gridbox with Q4 increments.

! Method:
!   Uses the method proposed in Annex B of the PC2 cloud scheme project
!   report, partitioned between liquid and ice phases by the method
!   given in Annex D. Condensate forcing is assumed to take place in
!   conjunction with cloudy air. This subroutine does NOT update vapour,
!   temperature and condensate values, which are assumed to have been
!   calculated by the physics scheme which is responsible for the call
!   to this subroutine. NOTE: This subroutine does NOT check that
!   increments calculated are sensible.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation

!  Global Variables:----------------------------------------------------
!      None

!  Subroutine Arguments:------------------------------------------------

  REAL ::                                                               &
                        !, INTENT(IN)
   qcl(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end),                           &
!       Liquid content (kg water per kg air)
   qcf(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end),                           &
!       Ice content (kg water per kg air)
   q4_l(          qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end),                           &
!       Rate of change of liquid content with time from the forcing
!       (kg kg-1 s-1)
   q4_f(          qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end),                           &
!       Rate of change of ice content with time from forcing mechanism
!       (kg kg-1 s-1)
   ls(            qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end),                           &
!       In cloud condensate content of the cloud which is forced into
!       the gridbox (liquid plus ice) (kg kg-1)
   g_l(           qdims%i_start:qdims%i_end,                            &
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end),                           &
!       Volume fraction of injected cloud which contains liquid
   g_f(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end)
!       Volume fraction of injected cloud which contains ice. Note that
!       G_L + G_F need not equal one if there is overlap between liquid
!       and ice

  REAL ::                                                               &
                        !, INTENT(INOUT)
   cf(            qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end),                           &
!       Total cloud fraction (no units)
   cfl(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end),                           &
!       Liquid cloud fraction (no units)
   cff(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  qdims%k_start:qdims%k_end)
!       Liquid cloud fraction (no units)

!  External functions:

!  Local parameters and other physical constants------------------------
  REAL, PARAMETER :: tolerance=1.0e-10 
!       Tolerance to avoid a divide by zero

!  Local scalars--------------------------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
  REAL ::                                                               &
   denom,                                                               &
!      Denominator in calculation
   q4
!      Sum of Q4_L and Q4_F

!  (b) Others.
  INTEGER :: k,i,j ! Loop counters:   K - vertical level index
!                                   I,J - horizontal position index

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

!- End of Header

  IF (lhook) CALL dr_hook('PC2_INHOMOGENEOUS',zhook_in,zhook_handle)

! ==Main Block==--------------------------------------------------------

! Loop round levels to be processed
  DO k = qdims%k_end, 1, -1

    DO j = qdims%j_start, qdims%j_end

      DO i = qdims%i_start, qdims%i_end

! Calculate total Q4. Only perform the calculation if the total Q4 is
! non-zero.

        q4 = q4_l(i,j,k) + q4_f(i,j,k)
        IF (q4  /=  0.0) THEN

! Calculate the change in total cloud fraction

          denom = (ls(i,j,k) - qcl(i,j,k) - qcf(i,j,k) )

          IF ( ABS(denom) > tolerance ) THEN

            denom = 1.0 / denom
            cf(i,j,k)  = cf(i,j,k)  +                                   &
                         ( 1.0        - cf(i,j,k)  ) * q4 * denom
            cfl(i,j,k) = cfl(i,j,k) +                                   &
                         ( g_l(i,j,k) - cfl(i,j,k) ) * q4 * denom
            cff(i,j,k) = cff(i,j,k) +                                   &
                         ( g_f(i,j,k) - cff(i,j,k) ) * q4 * denom

! Otherwise cloud fraction will go to one. In theory, cloud fraction
! can never be reduced by this process.

          ELSE
            cf(i,j,k)  = 1.0
            cfl(i,j,k) = 1.0
            cff(i,j,k) = 1.0
          END IF

        END IF

      END DO !i

    END DO !j

  END DO !k

! End of the subroutine

  IF (lhook) CALL dr_hook('PC2_INHOMOGENEOUS',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_inhomogeneous
