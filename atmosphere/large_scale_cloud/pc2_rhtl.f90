! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ PC2 Cloud Scheme: Calculation of RH(TL) for later use by initiation
! Subroutine Interface:
SUBROUTINE pc2_rhtl(                                                    &
!      Parallel variables
  halo_i, halo_j, off_x, off_y,                                         &
!      Array dimensions
 wet_levels, row_length, rows,                                          &
!      Prognostic arrays
 theta, exner_theta_levels, q, qcl, p_theta_levels,                     &
!      Output value of RH(TL)
 rhts,                                                                  &
!      Logical control
 l_mixing_ratio                                                         &
 )

  USE water_constants_mod,   ONLY: lc
  USE atmos_constants_mod,   ONLY: cp
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: pdims, qdims, tdims_s, qdims_l, pdims_s

  IMPLICIT NONE

! Purpose:
!   Calculate total relative humidity wrt T_L
!
! Method:
!   Straight calculation using RH_T(TL) = (q+qcl)/qsatwat(TL)
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
!
!  Subroutine Arguments:------------------------------------------------
  INTEGER           , INTENT(IN) ::                                     &
!      Model dimensions
   wet_levels, row_length, rows,                                        &
!      Parallel setup variables
    halo_i,                                                             &
!      Size of halo in i direction.
    halo_j,                                                             &
!      Size of halo in j direction.
    off_x,                                                              &
!      Size of small halo in i.
    off_y      
!      Size of small halo in j.

  REAL              , INTENT(IN) ::                                     &
!      Model prognostics
    theta(             tdims_s%i_start:tdims_s%i_end,                   & 
                       tdims_s%j_start:tdims_s%j_end,                   &  
                       tdims_s%k_start:tdims_s%k_end),                  &
    exner_theta_levels(pdims_s%i_start:pdims_s%i_end,                   & 
                       pdims_s%j_start:pdims_s%j_end,                   &  
                       pdims_s%k_start:pdims_s%k_end),                  &
    q(                 qdims_l%i_start:qdims_l%i_end,                   & 
                       qdims_l%j_start:qdims_l%j_end,                   &  
                       qdims_l%k_start:qdims_l%k_end),                  &
    qcl(               qdims_l%i_start:qdims_l%i_end,                   & 
                       qdims_l%j_start:qdims_l%j_end,                   &  
                       qdims_l%k_start:qdims_l%k_end),                  &
    p_theta_levels(    pdims_s%i_start:pdims_s%i_end,                   & 
                       pdims_s%j_start:pdims_s%j_end,                   &  
                       pdims_s%k_start:pdims_s%k_end)

  LOGICAL           , INTENT(IN) ::                                     &
   l_mixing_ratio     ! Use mixing ratio formulation

  REAL              , INTENT(OUT) ::                                    &
!      RH_T(T_L)
   rhts(row_length,rows,wet_levels)

!  External functions:

!  Local parameters and other physical constants------------------------
  REAL, PARAMETER :: lcrcp=lc/cp
!      Latent heat of condensation divided by heat capacity of air.

!  Local arrays---------------------------------------------------------

  REAL ::                                                               &
   tl(                 qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end),                      &
!      Liquid water temperature (K)
   qsl_tl(             qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end),                      &
!      Qsat wrt liquid water at temp TL
   p_no_halos(         pdims%i_start:pdims%i_end,                       & 
                       pdims%j_start:pdims%j_end)
!      Pressure without halo values (Pa)

!  Local scalars

  INTEGER :: k,i,j ! Loop counters:   K - vertical level index
!                                   I,J - horizontal position index

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! ==Main Block==--------------------------------------------------------

  IF (lhook) CALL dr_hook('PC2_RHTL',zhook_in,zhook_handle)

  DO k = 1, qdims%k_end

    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
!       Calculate liquid temperature TL
        tl(i,j) = theta(i,j,k) * exner_theta_levels(i,j,k)              &
                -lcrcp*qcl(i,j,k)
        p_no_halos(i,j) = p_theta_levels(i,j,k)
      END DO
    END DO

! Calculate qsat(TL) with respect to liquid water
! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix(qsl_tl,tl,p_no_halos,                             &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
!       Calculate RH_T(TL)
        rhts(i,j,k) = ( q(i,j,k) + qcl(i,j,k) ) / qsl_tl(i,j)
      END DO
    END DO

  END DO  ! k

! End of the subroutine

  IF (lhook) CALL dr_hook('PC2_RHTL',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_rhtl
