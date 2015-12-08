! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.  
! For further details please refer to the file COPYRIGHT.txt  
! which you should have received as part of this distribution.  
! *****************************COPYRIGHT*******************************
!    
MODULE eg_mass_conserv_mod
IMPLICIT NONE 
 
! Description:  
! ----------------------------------------------------------------------
!  This routine assumes rho_species = rho * q
!
!  therefore, it should work for both:
!
!            (a) rho = rho_dry & q = mixing_ratios, and
!            (b) rho = rho_wet & q = specific humidity
!                              
!
! The quantity we want to conserve is:
! 
! (1) SUM{ volume(i,j,k)*rho_np1(i,j,k)*qbar_np1(i,j,k) } = C (known value)
!
! where volume(i,j,k) is the volume of the box surrounding a rho-point
! or rho(i,j,k) is the centre of volume(i,j,k), and,
!
! (2) qbar(i,j,k) = alfa_za(k) * q(i,j,k) + beta_za(k) * q(i,j,k-1)
!
! Using (1) and (2) we can write the SUM in (1) as:
!
! (3) SUM( psi_np1(i,j,k) * q_np1(i,j,k) )  = C
!
! and
!
! (4) C = SUM(   psi_n(i,j,k) *    qs_n(i,j,k) )  +
!         SUM( psi_np1(i,j,k) *    qs_s(i,j,k) )  
!
!     Note that the source_term qs_s are at (n+1) and refers to Physics2
!     the sources due to physics1 is contained within qs_n
!
! ----------------------------------------------------------------------
!  
! Method: ENDGame formulation version 3.02
! Code Owner: See Unified Model Code Owners HTML page 
! This file belongs in section: Dynamics Advection 
!  
! Code description:
! Language: Fortran 95.  
! This code is written to UMDP3 standards.  

CONTAINS        

SUBROUTINE eg_mass_conservation(rho_n, rho_np1, qs_n, qs_np1, qs_s,   & 
                                qsmin, qswitches, number_qs,          &
                                L_conserv_smooth_lap                  )
                                                                                
     
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims, tdims_s
USE eg_helmholtz_mod,      ONLY: ec_vol
USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels
USE horiz_grid_mod,        ONLY: xi1_p, xi2_p
USE global_2d_sums_mod,    ONLY: global_2d_sums
USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE um_parvars,            ONLY: offx, offy
USE proc_info_mod,         ONLY: n_proc
USE Field_Types

IMPLICIT NONE
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
   
INTEGER, INTENT(IN) :: number_qs             
 
REAL,    INTENT(IN) ::   qs_n(tdims_s%i_start:tdims_s%i_end,          &
                              tdims_s%j_start:tdims_s%j_end,          &
                              tdims_s%k_start:tdims_s%k_end,          &
                              number_qs)

REAL,    INTENT(IN) ::   qs_s(tdims_s%i_start:tdims_s%i_end,          &
                              tdims_s%j_start:tdims_s%j_end,          &
                              tdims_s%k_start:tdims_s%k_end,          &
                              number_qs)

REAL,    INTENT(INOUT) :: qs_np1(tdims_s%i_start:tdims_s%i_end,       &
                                 tdims_s%j_start:tdims_s%j_end,       &
                                 tdims_s%k_start:tdims_s%k_end,       & 
                                 number_qs)

REAL,    INTENT(IN)    ::  rho_n(pdims_s%i_start:pdims_s%i_end,       &
                                 pdims_s%j_start:pdims_s%j_end,       &
                                 pdims_s%k_start:pdims_s%k_end)

REAL,    INTENT(IN)    :: rho_np1(pdims_s%i_start:pdims_s%i_end,      &
                                  pdims_s%j_start:pdims_s%j_end,      &
                                  pdims_s%k_start:pdims_s%k_end)

REAL,    INTENT(IN) :: qsmin (number_qs)
INTEGER, INTENT(IN) :: qswitches (number_qs)
LOGICAL, INTENT(IN) :: L_conserv_smooth_lap
  
! locals
   
REAL :: alfa_za(pdims%k_start:pdims%k_end)
REAL :: beta_za(pdims%k_start:pdims%k_end)

REAL :: psi_n  (tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end)
REAL :: psi_np1(tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end) 
     
REAL :: local_sums(tdims%i_start:tdims%i_end,                         &
                   tdims%j_start:tdims%j_end,                         &
                   number_qs,3)

REAL :: global_sums(number_qs,3)

REAL :: local_sums2(tdims%i_start:tdims%i_end,                        &
                    tdims%j_start:tdims%j_end,                        &
                    number_qs+1)

REAL :: global_sums2(number_qs+1)
         
REAL    :: dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,delsq,gama

REAL    :: dx2(tdims%i_start  :tdims%i_end)
REAL    :: dx3(tdims%i_start-1:tdims%i_end)
REAL    :: dy2(tdims%j_start  :tdims%j_end)
REAL    :: dy3(tdims%j_start-1:tdims%j_end)
REAL    :: dz2(tdims%k_start  :tdims%k_end)
REAL    :: dz3(tdims%k_start  :tdims%k_end)
   
REAL    :: mass_deficit, temp1, temp2
REAL    :: one_over_total_lap

REAL    :: small_tol = 1.0E-30

REAL    :: lap (tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end,number_qs)

REAL    ::  qs_lap(tdims%i_start-2:tdims%i_end+2,                     &
                   tdims%j_start-2:tdims%j_end+2,                     &
                   tdims%k_start:tdims%k_end,                         & 
                   number_qs)

INTEGER :: i,j,k,kk
INTEGER :: filter_data = 1 ! default

IF (lhook) CALL dr_hook('EG_MASS_CONSERVATION',zhook_in,zhook_handle)
 
! ------------------------------------------------------------------
! Section 1 : make qs(n+1) > qsmin if necessary  
! ------------------------------------------------------------------
  
DO kk = 1, number_qs
  qs_np1(:,:,:,kk) = MAX(qs_np1(:,:,:,kk),qsmin(kk))
END DO
     
!------------------------------------------------------------------
! section 2: Compute the vertical averaging weigths for Eq.(2) above
!------------------------------------------------------------------

DO k = pdims%k_start, pdims%k_end

  alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) ) /    &
               ( eta_theta_levels(k) - eta_theta_levels(k-1) )
  beta_za(k) = 1.0 - alfa_za(k)      

END DO

!-------------------------------------------------------------
! section 3: Compute the 3D array psi(i,j,k) for Eqs (3) & (4)
!-------------------------------------------------------------
      
k = tdims%k_start

DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end 
    psi_n(i,j,k) =   rho_n(i,j,k+1) * ec_vol(i,j,k+1) * beta_za(k+1)   
    psi_np1(i,j,k) = rho_np1(i,j,k+1) * ec_vol(i,j,k+1) * beta_za(k+1)      
  END DO
END DO

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)       &
!$OMP& SHARED(tdims,psi_n,rho_n,ec_vol,alfa_za,beta_za,psi_np1,rho_np1)
DO k = tdims%k_start + 1, tdims%k_end - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k)  =  rho_n(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )    &
                    +  rho_n(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
  
      psi_np1(i,j,k) = rho_np1(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )  &
                     + rho_np1(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

k = tdims%k_end
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    psi_n(i,j,k)   =   rho_n(i,j,k) * ec_vol(i,j,k) * alfa_za(k)   
    psi_np1(i,j,k) = rho_np1(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
  END DO
END DO
      
! --------------------------------------------------------------------
!  section 4:  Calculate the laplacian (div) of qs. This div is computed 
!              like Cartesian (and this consistent with interpolation). 
!              These calculations can be simplified for constant mesh  
!              (here a general(variable) mesh is assumed)  
! --------------------------------------------------------------------
  
DO i = tdims%i_start, tdims%i_end  
  dx2(i) = 0.5 * ( xi1_p(i+1) - xi1_p(i-1) )
END DO

DO i = tdims%i_start - 1 , tdims%i_end  
  dx3(i) = xi1_p(i+1) - xi1_p(i)
END DO

DO j = tdims%j_start, tdims%j_end
  dy2(j) = 0.5 * ( xi2_p(j+1) - xi2_p(j-1) )
END DO

DO j = tdims%j_start - 1, tdims%j_end
  dy3(j) = xi2_p(j+1) - xi2_p(j)
END DO

DO k = tdims%k_start+1 , tdims%k_end - 1
  dz2(k) = 0.5 * ( eta_theta_levels(k+1) - eta_theta_levels(k-1) )
END DO

DO k = tdims%k_start , tdims%k_end - 1
  dz3(k) = eta_theta_levels(k+1) - eta_theta_levels(k)  
END DO

lap = 0.0

IF ( filter_data == 1 .AND. L_conserv_smooth_lap )  THEN
  
! Compute the laplacian of filtered qs_np1 (using 1-2-1 filter in 
! all directions). This is equivalent to apply the Laplacian to
! a coarse grid (similar to Multigrid).
! These filtered data are used only in computing "lap" as
! a measure of smoothness -- this is to avoid giving too much
! weight to sub-grid noise

  DO kk = 1, number_qs
!$OMP PARALLEL DO  SCHEDULE(STATIC) DEFAULT(NONE)               &
!$OMP& PRIVATE(i,j,k) SHARED(tdims,qs_lap,qs_np1,kk)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_lap(i,j,k,kk) = qs_np1(i,j,k,kk) 
        END DO
      END DO
    END DO 
!$OMP END PARALLEL DO
  END DO

  i =  1 + tdims%i_end - tdims%i_start
  j =  1 + tdims%j_end - tdims%j_start
  k = (1 + tdims%k_end - tdims%k_start)*number_qs
 
!DEPENDS ON: swap_bounds
  CALL swap_bounds(qs_lap,i,j,k,2,2,fld_type_p,.FALSE.)


  DO kk = 1, number_qs         

    IF ( qswitches(kk) == 1 ) THEN

      k = tdims%k_start
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/ dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/ dx3(i-1)
  
          dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)
  
          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +         &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j)  
           
          lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)
        
        END DO
      END DO
      
      k = tdims%k_start + 1
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 

            dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
            dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/dx3(i-1)

            dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
            dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

            dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
            dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-1,kk))/dz3(k-1)
                              
            delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                     (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                     (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k) 
        
            lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)
        
         END DO
      END DO
      

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP& PRIVATE(i,j,k,dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,delsq)         &
!$OMP& SHARED(tdims,qs_lap,dx3,dx2,dy3,dy2,dz3,dz2,lap,qs_np1,psi_np1,  &
!$OMP& kk)
   
      DO k = tdims%k_start + 2, tdims%k_end - 2
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end 

            dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
            dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/dx3(i-1)

            dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
            dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

            dfdzp  = (qs_lap(i,j,k+2,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
            dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-2,kk))/dz3(k-1)
                              
            delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                     (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                     (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k) 
        
            lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)
         
          END DO
        END DO
      END DO 
      
!$OMP END PARALLEL DO

      
      k = tdims%k_end - 1
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 

            dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
            dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/dx3(i-1)

            dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
            dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

            dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
            dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-1,kk))/dz3(k-1)
                              
            delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                     (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                     (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k) 
        
            lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)
                
         END DO
      END DO
      
      k = tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)

          dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)
  
          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) + &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) 
          
          lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)
        
        END DO
      END DO 

    END IF
  END DO
  
     
ELSE 

    IF ( filter_data == 1  )  THEN

    ! filter qs_np1 (apply 1-2-1 filter in all directions)
    
  DO kk = 1, number_qs
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                      &
!$OMP& PRIVATE(i,j,k) SHARED(tdims,qs_lap,qs_np1,kk)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          qs_lap(i,j,k,kk) = ( qs_np1(i-1,j-1,k,kk)                   &
                                   + qs_np1(i+1,j-1,k,kk)             &
                                   + qs_np1(i-1,j+1,k,kk)             &
                                   + qs_np1(i+1,j+1,k,kk)             &
                                   + 2.0*( qs_np1(i,j-1,k,kk)         &
                                         + qs_np1(i,j+1,k,kk)         &
                                         + qs_np1(i-1,j,k,kk)         &
                                         + qs_np1(i+1,j,k,kk) )       &
                                   + 4.0*qs_np1(i,j,k,kk) )/16.0
        END DO
      END DO
    END DO 
!$OMP END PARALLEL DO
  END DO

   i = 1 + tdims%i_end - tdims%i_start
   j = 1 + tdims%j_end - tdims%j_start
   k = (1 + tdims%k_end - tdims%k_start)*number_qs
  
!DEPENDS ON: swap_bounds
   CALL swap_bounds(qs_lap,i,j,k,2,2,fld_type_p,.FALSE.)

   ELSE !  use unfiltred raw data
  
   DO kk = 1, number_qs
!$OMP PARALLEL DO  SCHEDULE(STATIC) DEFAULT(NONE)                     &
!$OMP& PRIVATE(i,j,k) SHARED(tdims,qs_lap,qs_np1,kk)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start-1, tdims%j_end+1
        DO i = tdims%i_start-1, tdims%i_end+1
          qs_lap(i,j,k,kk) = qs_np1(i,j,k,kk) 
        END DO
      END DO
    END DO 
!$OMP END PARALLEL DO
  END DO

  END IF 
    
  DO kk = 1, number_qs         

    IF ( qswitches(kk) == 1 ) THEN

      k = tdims%k_start
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)
  
          dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)
  
          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +           &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j)  
  
          lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)
         
        END DO
      END DO

!$OMP PARALLEL DO  SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& PRIVATE(i,j,k, dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,              &
!$OMP& delsq) SHARED(tdims,qs_lap,dx3,dx2,dy3,dy2,dz3,                  &
!$OMP& dz2,lap,qs_np1,psi_np1,kk)
   
      DO k = tdims%k_start + 1, tdims%k_end - 1
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end 

            dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
            dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)

            dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
            dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)

            dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
            dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-1,kk))/dz3(k-1)
                              
            delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                     (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                     (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k) 

            lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)
        
          END DO
        END DO
      END DO 
      
!$OMP END PARALLEL DO

      k = tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)

          dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)
  
          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) + &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) 
  
          lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)
           
        END DO
      END DO 

    END IF
  END DO

END IF 

! ----------------------------------------------------------------------
!  section 5:  Calculate the global sums in (3) & (4) and total lap
! ----------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k,kk) DEFAULT(NONE)    &
!$OMP& SHARED(local_sums,tdims,psi_n,psi_np1,qs_n,qs_s,qs_np1,lap,    &
!$OMP& qswitches,number_qs)
DO kk = 1, number_qs

  local_sums(:,:,kk,1) = 0.
  local_sums(:,:,kk,2) = 0.
  local_sums(:,:,kk,3) = 0.

  IF ( qswitches(kk) == 1 ) THEN 
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end 

          local_sums(i,j,kk,1) =  local_sums(i,j,kk,1)                &
                                  +    psi_n(i,j,k) * qs_n(i,j,k,kk)  &
                                  +  psi_np1(i,j,k) * qs_s(i,j,k,kk)

          local_sums(i,j,kk,2) =  local_sums(i,j,kk,2)                &
                                  +  psi_np1(i,j,k) * qs_np1(i,j,k,kk) 

          local_sums(i,j,kk,3) =  local_sums(i,j,kk,3)                &
                                    +  lap(i,j,k,kk)
        END DO
      END DO
    END DO 
  END IF
END DO 
!$OMP END PARALLEL DO


CALL global_2d_sums(local_sums, tdims%i_end-tdims%i_start+1,          &
                                       tdims%j_end-tdims%j_start+1,   &
                                       0, 0, 3*number_qs, global_sums )

! ----------------------------------------------------------------------
!  section 6:  Modify qs to achieve mass conservation
! ---------------------------------------------------------------------- 

DO kk = 1, number_qs 
  IF ( qswitches(kk) == 1 .AND. ABS(global_sums(kk,3)) > small_tol ) THEN 
    mass_deficit = global_sums(kk,1) -  global_sums(kk,2)
    one_over_total_lap = 1.0/global_sums(kk,3)
!$OMP PARALLEL DO  SCHEDULE(STATIC) DEFAULT(NONE)                     & 
!$OMP& PRIVATE(i,j,k, gama) SHARED(tdims,lap,mass_deficit,            &
!$OMP& one_over_total_lap,qs_np1,psi_np1,kk)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end 
          gama = lap(i,j,k,kk) * mass_deficit * one_over_total_lap
          qs_np1(i,j,k,kk) = qs_np1(i,j,k,kk) + gama/psi_np1(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF
END DO 

!=======================================================      
!  section 7: Apply a second correction to force q >= qmin
!      while maintaining mass conservation      
!=======================================================

! clip the field points that are below qsmin

DO kk = 1, number_qs
  qs_np1(:,:,:,kk) = max(qs_np1(:,:,:,kk),qsmin(kk))
END DO
 
! compute new masses at (n+1) after cliping
 
DO kk = 1, number_qs   
  local_sums2(:,:,kk) = 0.0      
  DO k = tdims%k_start, tdims%k_end     
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end     

        local_sums2(i,j,kk) = local_sums2(i,j,kk)                    &
                            + psi_np1(i,j,k) * qs_np1(i,j,k,kk)  
      END DO 
    END DO    
  END DO 
END DO 

! this sum doesn't invlove qs (it's the sum above with q=1)

local_sums2(:,:,number_qs+1) = 0.0      
DO k = tdims%k_start, tdims%k_end     
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end    
      local_sums2(i,j,number_qs+1) = local_sums2(i,j,number_qs+1)    &
                                   + psi_np1(i,j,k)  
    END DO 
  END DO    
END DO  
      
CALL global_2d_sums(local_sums2, tdims%i_end-tdims%i_start+1,         &
                    tdims%j_end-tdims%j_start+1,                      &
                    0, 0, number_qs+1, global_sums2 )

DO kk = 1, number_qs      

  mass_deficit = global_sums2(kk)  - qsmin(kk)*global_sums2(number_qs+1)

  IF ( qswitches(kk) == 1 .AND. ABS(mass_deficit) > small_tol ) THEN
    temp1 = ( global_sums(kk,1)-qsmin(kk)*global_sums2(number_qs+1) )/&
            mass_deficit
    temp2 = qsmin(kk)*(1.0 - temp1)   

    qs_np1(:,:,:,kk) = temp1 * qs_np1(:,:,:,kk) + temp2 

  END IF   
END DO  
      
IF (lhook) CALL dr_hook('EG_MASS_CONSERVATION',zhook_out,zhook_handle)

RETURN

END SUBROUTINE eg_mass_conservation  
END MODULE eg_mass_conserv_mod

