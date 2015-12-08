! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Cloud Scheme: Forcing due to advection (adiabatic cooling)
! Subroutine Interface:
SUBROUTINE pc2_assim(                                                   &
 timestep,                                                              &
                                                  ! Timestep
 l_pc2_homog, l_pc2_cff, l_mixing_ratio,                                &
                                                  ! Control logicals
 t, cf, cfl, cff, q, qcl, qcf, p,                                       &
                                                  ! Prognostic Fields
 deltat, deltaq, deltaqcl, deltaqcf, deltap                             &
                                                  ! Forcing quantities
                                                  ! No diagnostics yet
 )

  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: pdims, qdims, tdims

  IMPLICIT NONE

! Purpose:
!   This subroutine acts to provide the PC2 cloud fraction and
!   condensation increments interfacting to the data assimilation
!   AC and VAR schemes.

! Method:
!   Uses the homogenous forcing routine to calculate condensation
!   and associated cloud fraction changes and performs an estimate
!   of ice cloud fraction changes based on the ice water content
!   forcing. Parts of the code can be switched with logicals to allow
!   future data assimilation schemes only to use the parts that they
!   need.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
! Documentation: New Cloud Scheme Documentation

!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
  LOGICAL ::                                                            &
                        !, INTENT(IN)
   l_pc2_homog, l_pc2_cff,                                              &
!      Control on which sections of code need using
   l_mixing_ratio  
!      Use mixing ratio formulation

  REAL ::                                                               &
                        !, INTENT(IN)
   timestep,                                                            &
!      Model timestep (s)
   deltat(  tdims%i_start:tdims%i_end,                                  & 
            tdims%j_start:tdims%j_end,                                  &  
                        1:tdims%k_end),                                 &
!      Forcing of temperature over the timestep (K)
   deltaq  (qdims%i_start:qdims%i_end,                                  & 
            qdims%j_start:qdims%j_end,                                  &  
                        1:qdims%k_end),                                 &
!      Forcing of vapour (kg kg-1)
   deltaqcl(qdims%i_start:qdims%i_end,                                  & 
            qdims%j_start:qdims%j_end,                                  &  
                        1:qdims%k_end),                                 &
!      Forcing of liquid (kg kg-1)
   deltaqcf(qdims%i_start:qdims%i_end,                                  & 
            qdims%j_start:qdims%j_end,                                  &  
                        1:qdims%k_end),                                 &
!      Forcing of ice    (kg kg-1)
   deltap  (pdims%i_start:pdims%i_end,                                  & 
            pdims%j_start:pdims%j_end,                                  &  
            pdims%k_start:pdims%k_end)  
!      Forcing of pressure    (Pa)

  REAL ::                                                               &
                        !, INTENT(IN)
   qcf(     qdims%i_start:qdims%i_end,                                  & 
            qdims%j_start:qdims%j_end,                                  &  
                        1:qdims%k_end),                                 &
!      Ice content          (kg kg-1)
   p(       pdims%i_start:pdims%i_end,                                  & 
            pdims%j_start:pdims%j_end,                                  &  
            pdims%k_start:pdims%k_end)
!      Pressure at theta levels  (Pa)

  REAL ::                                                               &
                        !, INTENT(INOUT)
   t(       tdims%i_start:tdims%i_end,                                  & 
            tdims%j_start:tdims%j_end,                                  &  
                        1:tdims%k_end),                                 &
!      Dry-bulb temperature (K)
   cf(      qdims%i_start:qdims%i_end,                                  & 
            qdims%j_start:qdims%j_end,                                  &  
                        1:qdims%k_end),                                 &
!      Bulk cloud fraction  (no units)
   cfl(     qdims%i_start:qdims%i_end,                                  & 
            qdims%j_start:qdims%j_end,                                  &  
                        1:qdims%k_end),                                 &
!      Liquid cloud frac    (no units)
   cff(     qdims%i_start:qdims%i_end,                                  & 
            qdims%j_start:qdims%j_end,                                  &  
                        1:qdims%k_end),                                 &
!      Ice cloud fraction   (no units)
   q(       qdims%i_start:qdims%i_end,                                  & 
            qdims%j_start:qdims%j_end,                                  &  
                        1:qdims%k_end),                                 &
!      Vapour content       (kg kg-1)
   qcl(     qdims%i_start:qdims%i_end,                                  & 
            qdims%j_start:qdims%j_end,                                  &  
                        1:qdims%k_end)
!      Liquid content       (kg kg-1)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

!  External subroutine calls: ------------------------------------------
      EXTERNAL pc2_homog_plus_turb
!
!  Local parameters-----------------------------------------------------

  REAL, PARAMETER:: qcf0 = 1.0e-4 
!      Specified in-cloud ice content (kg kg-1)

!  Local scalars--------------------------------------------------------

  INTEGER ::                                                            &
   points,                                                              &
!      Counter for qcf changes
   i,j,k
!      Loop counters

!  Local dynamic arrays-------------------------------------------------
  REAL ::                                                               &
   deltacff_c( (1+qdims%i_end-qdims%i_start)                            &
              *(1+qdims%j_end-qdims%j_start)                            &
              *(  qdims%k_end) ),                                       &
!      Change in ice cloud fraction
   cf_c (      (1+qdims%i_end-qdims%i_start)                            &
              *(1+qdims%j_end-qdims%j_start)                            &
              *(  qdims%k_end) ),                                       &
!      Compressed points
   cfl_c(      (1+qdims%i_end-qdims%i_start)                            &
              *(1+qdims%j_end-qdims%j_start)                            &
              *(  qdims%k_end) ),                                       &
!      Copies of cloud
   cff_c(      (1+qdims%i_end-qdims%i_start)                            &
              *(1+qdims%j_end-qdims%j_start)                            &
              *(  qdims%k_end) ),                                       &
!      Fraction variables
   zeros(      (1+qdims%i_end-qdims%i_start)                            &
              *(1+qdims%j_end-qdims%j_start)                            &
              *(  qdims%k_end) )
!      Array of zeros

  REAL ::                                                               &
   t_copy(      tdims%i_start:tdims%i_end,                              & 
                tdims%j_start:tdims%j_end,                              &  
                            1:tdims%k_end),                             &
! Copy of T
   q_copy(      qdims%i_start:qdims%i_end,                              & 
                qdims%j_start:qdims%j_end,                              &  
                            1:qdims%k_end)    
! Copy of q

!- End of Header

! ==Main Block==--------------------------------------------------------

  IF (lhook) CALL dr_hook('PC2_ASSIM',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! 1. Call homogenous forcing routine if needed
! ----------------------------------------------------------------------
  IF (l_pc2_homog) THEN

! In order that temperature (t) and vapour (q) are not updated by
! the condensation response to assimilation, we take copies of these
! variables to use in the forcing subroutine. The subroutine will
! update T_copy etc., not T. The cloud fractions and qcl *are* updated.

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          t_copy(i,j,k) = t(i,j,k)
          q_copy(i,j,k) = q(i,j,k)
        END DO
      END DO
    END DO

! DEPENDS ON: pc2_homog_plus_turb
    CALL pc2_homog_plus_turb(p, qdims%k_end, timestep,                  &
                             t_copy, cf, cfl, cff, q_copy, qcl,         &
                             deltat, deltaq, deltaqcl, deltap,          &
                             0.0, 0.0,                                  &
                             l_mixing_ratio)

! Now update q and T for the forcings

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          t(i,j,k) = t(i,j,k) + deltat(i,j,k)
          q(i,j,k) = q(i,j,k) + deltaq(i,j,k)
        END DO
      END DO
    END DO

  END IF  ! L_pc2_homog

  IF (l_pc2_cff) THEN

! ----------------------------------------------------------------------
! 2. Calculate changes to ice cloud fraction
! ----------------------------------------------------------------------

! Estimate the in-cloud ice content of the part of the cloud that is
! being added or removed to be qcf_ic = cf * (qcf/cf) + (1-cf) * qcf0
! where qcf0 is a specified in-cloud content.

! Loop over each point

    points=0
    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          IF (deltaqcf(i,j,k) /= 0.0) THEN
            points = points + 1
                ! Store compressed points versions
            cf_c (points) = cf (i,j,k)
            cfl_c(points) = cfl(i,j,k)
            cff_c(points) = cff(i,j,k)
                ! Adjust cff for each point
            deltacff_c(points) = deltaqcf(i,j,k)                        &
                         / ( qcf(i,j,k) + (1.0-cff(i,j,k))*qcf0 )
                ! Check that cff is sensible
            deltacff_c(points) =                                        &
              MAX(  MIN(  deltacff_c(points), (1.0 - cff(i,j,k))  ),    &
                   -cff(i,j,k))
              ! Form zero array for forcing of liquid cloud fraction.
              ! Note that pc2_homog_plus_turb calls pc2_total_cf for the
              ! liquid changes
            zeros(points) = 0.0
          END IF
        END DO
      END DO
    END DO

! ----------------------------------------------------------------------
! 3. Calculate changes to total cloud fraction to update cf_c
! ----------------------------------------------------------------------

! DEPENDS ON: pc2_total_cf
    CALL pc2_total_cf(                                                  &
          points,cfl_c,cff_c,zeros,deltacff_c,cf_c)

! ----------------------------------------------------------------------
! 4. Scatter back changes to ice and total cloud fractions
! ----------------------------------------------------------------------

    points=0
    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          IF (deltaqcf(i,j,k)  /=  0.0) THEN
            points=points+1
            cff(i,j,k) = cff(i,j,k) + deltacff_c(points)
            cf (i,j,k) = cf_c(points)
          END IF
        END DO
      END DO
    END DO

  END IF  ! L_pc2_cff

! ----------------------------------------------------------------------
! 5. Call pc2_checks to make sure values are self consistent
! ----------------------------------------------------------------------
  IF (L_pc2_homog .OR. L_pc2_cff) THEN
! DEPENDS ON: pc2_checks
    CALL pc2_checks(p,                                                  &
        t, cf, cfl, cff, q, qcl, qcf,                                   &
        l_mixing_ratio)

  END IF  ! L_pc2_homog .or. L_pc2_cff

! End of the subroutine

9999 CONTINUE ! Error exit
  IF (lhook) CALL dr_hook('PC2_ASSIM',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_assim
