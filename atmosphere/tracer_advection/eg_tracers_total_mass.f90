! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.  
! For further details please refer to the file COPYRIGHT.txt  
! which you should have received as part of this distribution.  
! *****************************COPYRIGHT*******************************
!    
       MODULE eg_tracers_mass_mod
       IMPLICIT NONE 
 
! Description:  
!             This routine computes the total mass of species qs, 
!             massqs = sum[rho*av(q)*volume]. 
!             Note that qs are averaged to rho-points
!
! Method: ENDGame formulation version 3.02
! Code Owner: See Unified Model Code Owners HTML page 
! This file belongs in section: Tracer Advection
!  
! Code description:
! Language: Fortran 95.  
! This code is written to UMDP3 standards.
  
       CONTAINS
       SUBROUTINE eg_tracers_total_mass(rho, qs, massqs, number_qs)   
                   
       USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims, tdims_s
       USE eg_helmholtz_mod,      ONLY: ec_vol
       USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels
       USE horiz_grid_mod,        ONLY: xi1_p, xi2_p
       USE global_2d_sums_mod,    ONLY: global_2d_sums

       IMPLICIT NONE

       INTEGER, INTENT(IN)  :: number_qs                          
       REAL, INTENT(IN) :: rho(pdims_s%i_start:pdims_s%i_end,         &
                               pdims_s%j_start:pdims_s%j_end,         &
                               pdims_s%k_start:pdims_s%k_end)
                             
       REAL, INTENT(IN) :: qs(tdims_s%i_start:tdims_s%i_end,          &
                              tdims_s%j_start:tdims_s%j_end,          &
                              tdims_s%k_start:tdims_s%k_end,          &
                              number_qs)
                                       
       REAL, INTENT(OUT)     :: massqs(number_qs)
       
! locals temp variables

       REAL :: local_sums(pdims%i_start:pdims%i_end,                  &
                          pdims%j_start:pdims%j_end,                  &
                          number_qs)

       REAL :: alfa_za(pdims%k_start:pdims%k_end),                    &
               beta_za(pdims%k_start:pdims%k_end)        
                                                          
       INTEGER  :: i,j,k,kk 
       REAL     :: t1
       
! compute vertical averaging weights
              
       DO k = pdims%k_start, pdims%k_end
         alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) )/&
                      ( eta_theta_levels(k) - eta_theta_levels(k-1) )
         beta_za(k) = 1.0 - alfa_za(k)             
       END DO

! compute total masses of each species q
            
       DO kk = 1, number_qs               
         local_sums(:,:,kk)   = 0.0                                       
         DO k = pdims%k_start, pdims%k_end     
           DO j = pdims%j_start, pdims%j_end  
             DO i = pdims%i_start, pdims%i_end
               t1 =  alfa_za(k)*  qs(i,j,k,kk) +                      &
                     beta_za(k)*  qs(i,j,k-1,kk)    

               local_sums(i,j,kk)  = local_sums(i,j,kk) +             &
                                     t1*rho(i,j,k)*ec_vol(i,j,k) 
             END DO                                                         
           END DO                                                            
         END DO 
       END DO
      
       CALL global_2d_sums(local_sums,pdims%i_end-pdims%i_start+1,    &
                                      pdims%j_end-pdims%j_start+1,    &
                                      0,0,number_qs,massqs         )
                                                                           
       RETURN                                                               
       END SUBROUTINE eg_tracers_total_mass
       END MODULE eg_tracers_mass_mod
