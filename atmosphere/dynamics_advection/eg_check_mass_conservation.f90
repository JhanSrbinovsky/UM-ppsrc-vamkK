! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.  
! For further details please refer to the file COPYRIGHT.txt  
! which you should have received as part of this distribution.  
! *****************************COPYRIGHT*******************************
!
MODULE eg_check_conserv_mod
IMPLICIT NONE 
!
! Description: 
!
!   Check the error in the mass conservation constraint
!
!-------------------------------------------------------------------------
!   check mass conservation, i.e., check whether the following holds:
! 
!    sum[rho_np1 * av(qs_np1) * volume] = sum[rho_n   * av(qs_n ) * volume ]
!                                       + sum[rho_n   * av(qs_s1) * volume ] 
!                                       + sum[rho_np1 * av(qs_s2) * volume ],
!
! where:
!
! rho_n   == density at level (n)
! qs_n    == tracers at level (n)
! rho_np1 == density at level (n+1)
! qs_np1  == tracers at level (n+1)
! qs_s1   == sources at level (n)  [contribution of physics1, slow-physics]
! qs_s2   == sources at level (n+1)[contribution of physics2, fast-physics]
!
! Note that the argument "qs_n" of this routine contains the sum of qs_n and 
! qs_s1 (i.e., qs_n = qs_n + qs_s1) and qs_s == qs_s2 (physics2 contribution)  
!-------------------------------------------------------------------------
!
! Method: ENDGame formulation version 3.02
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90.
!   This code is written to UMDP3 v8.0 programming standards.
!  

CONTAINS
SUBROUTINE eg_check_mass_conservation(rho_n, rho_np1, qs_n, qs_np1,   &
                                      qs_s, qswitches, number_qs,     &
                                      mype, chara_str )

USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims, tdims_s
USE eg_helmholtz_mod,      ONLY: ec_vol
USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels
USE global_2d_sums_mod,    ONLY: global_2d_sums
USE parkind1,              ONLY: jpim, jprb !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE proc_info_mod,         ONLY: n_proc
USE PrintStatus_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
 
INTEGER, INTENT(IN)           ::  mype, number_qs 

REAL,    INTENT(IN)           ::   qs_n(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end,&
                                        number_qs)

REAL,    INTENT(IN)           ::   qs_s(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end,&
                                        number_qs)

REAL,    INTENT(IN)           :: qs_np1(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end,&
                                        number_qs)
    
REAL,    INTENT(IN)           ::  rho_n(pdims_s%i_start:pdims_s%i_end,&
                                        pdims_s%j_start:pdims_s%j_end,&
                                        pdims_s%k_start:pdims_s%k_end)

REAL,    INTENT(IN)          :: rho_np1(pdims_s%i_start:pdims_s%i_end,&
                                        pdims_s%j_start:pdims_s%j_end,&
                                        pdims_s%k_start:pdims_s%k_end)
     
INTEGER, INTENT(IN)          :: qswitches(number_qs)
CHARACTER(LEN=*)                 :: chara_str
  
! locals 

REAL :: local_sums(tdims%i_start:tdims%i_end,                         &
                   tdims%j_start:tdims%j_end,                         &
                   number_qs,2)

REAL :: global_sums(number_qs,2)
REAL :: err_mass(number_qs) 
REAL :: qsmin(number_qs) 
REAL :: qsmax(number_qs) 
 
REAL :: alfa_za(pdims%k_start:pdims%k_end)
REAL :: beta_za(pdims%k_start:pdims%k_end)

INTEGER  :: i,j,k,kk
REAL     :: t1, t2, t3
INTEGER  :: istat1, istat2
 
IF (lhook) CALL dr_hook('EG_CHECK_MASS_CONSERVATION',zhook_in,        &
                        zhook_handle)

! compute the vertical averaging weights
 
DO k = pdims%k_start, pdims%k_end
  alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) ) /      &
               ( eta_theta_levels(k) - eta_theta_levels(k-1) )
  beta_za(k) = 1.0 - alfa_za(k) 
END DO

! compute the total integrals (mass)
! sums(:,:,:,1) = mass of rho at level n + sources
! sums(:,:,:,2) = mass of rho at level (n+1)

DO kk = 1, number_qs  
  
  local_sums(:,:,kk,1) = 0.0 
  local_sums(:,:,kk,2) = 0.0
   
  IF ( qswitches(kk) == 1 ) THEN   
    DO k = pdims%k_start, pdims%k_end     
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
 
          t1 =  alfa_za(k)*  qs_n(i,j,k,kk) +                         &
                beta_za(k)*  qs_n(i,j,k-1,kk)

          t2 =  alfa_za(k)*  qs_s(i,j,k,kk) +                         &
                beta_za(k)*  qs_s(i,j,k-1,kk)

          t3 =  alfa_za(k)*qs_np1(i,j,k,kk) +                         &
                beta_za(k)*qs_np1(i,j,k-1,kk)

          local_sums(i,j,kk,1)  = local_sums(i,j,kk,1)                &
                                + t1 * rho_n(i,j,k) * ec_vol(i,j,k)   & 
                                + t2 * rho_np1(i,j,k) * ec_vol(i,j,k) 

          local_sums(i,j,kk,2)  = local_sums(i,j,kk,2)                &
                                + t3 * rho_np1(i,j,k) * ec_vol(i,j,k)    
        END DO   
      END DO 
    END DO
  END IF 
END DO
 
CALL global_2d_sums(local_sums, tdims%i_end-tdims%i_start+1,          &
                    tdims%j_end-tdims%j_start+1,                      &
                    0, 0, 2*number_qs, global_sums )

! compute the error in mass conservation and print it if required
! also in here the min/max fields are computed and printed for
! diagnosis/monitoring purposes
 
DO kk = 1, number_qs 

  IF ( qswitches(kk) == 0 ) THEN
    err_mass(kk) =  0.0
  ELSE
    IF ( ABS(global_sums(kk,1)) < 1.0e-30 ) THEN 
      err_mass(kk) =  global_sums(kk,2)-global_sums(kk,1)
    ELSE   
      err_mass(kk) =  (global_sums(kk,2)-global_sums(kk,1))/       &
                      ABS(global_sums(kk,1))
    END IF
  END IF

  IF ( mype == 0 .AND. PrintStatus > PrStatus_Normal ) THEN  
    WRITE(6,fmt='(A,I4,2x,E23.15)')                                  &
  chara_str//' species / mass conservation error  = ',kk, err_mass(kk)
  END IF
END DO
 
DO kk = 1, number_qs 
  IF ( qswitches(kk) == 0 ) THEN 
    qsmin(kk) = 0.0  
    qsmax(kk) = 0.0 
  ELSE                      
    qsmin(kk) = MINVAL(qs_np1(:,:,:,kk))  
    qsmax(kk) = MAXVAL(qs_np1(:,:,:,kk)) 
  END IF  
END DO 
         
CALL gc_rmin(number_qs,n_proc,istat1,qsmin) 
CALL gc_rmax(number_qs,n_proc,istat2,qsmax) 
        
IF ( mype == 0 .AND. PrintStatus > PrStatus_Normal ) THEN         
  DO kk = 1, number_qs 
    WRITE(6,fmt='(A,I4,2(1x,E12.5))') chara_str//' species/min/max =', &
                                      kk, qsmin(kk),qsmax(kk)
  END DO                          
END IF

IF (lhook) CALL dr_hook('EG_CHECK_MASS_CONSERVATION',zhook_out,zhook_handle)
RETURN   
END SUBROUTINE eg_check_mass_conservation
END MODULE eg_check_conserv_mod

