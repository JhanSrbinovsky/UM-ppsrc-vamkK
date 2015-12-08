! *****************************COPYRIGHT******************************* 
! 
! (c) [University of Cambridge] [2008]. All rights reserved. 
! This routine has been licensed to the Met Office for use and 
! distribution under the UKCA collaboration agreement, subject  
! to the terms and conditions set out therein. 
! [Met Office Ref SC138]  
! 
! *****************************COPYRIGHT******************************* 
!                                                                       
! Purpose: Subroutine to call the UKCA lightning emissions routine,     
!          which diagnoses NOx lightning emissions and returns values   
!          to UKCA_EMISSION_CTL.                                        
!          Based on routine provided by Olaf Morgenstern.               
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
!          Called from UKCA_EMISSION_CTL.                               
!                                                                       
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!                                                                       
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
        SUBROUTINE UKCA_LIGHT_CTL(                                     &
              rows, row_length, delta_lambda, delta_phi,               &
              p_levels, conv_base_lev,                                 &
              conv_top_lev, mask, lat, area,                           &
              r_theta_levels, r_rho_levels,                            &
              p_centres, p_boundaries,                                 &
              an)                                                      

        USE earth_constants_mod,   ONLY: g, earth_radius
        USE UKCA_CONSTANTS,        ONLY: c_no2, c_n
        USE yomhook,               ONLY: lhook, dr_hook
        USE parkind1,              ONLY: jprb, jpim
        USE asad_chem_flux_diags
        IMPLICIT NONE                                                   

        INTEGER, INTENT(IN) :: rows                           
! No of latitudes
        INTEGER, INTENT(IN) :: row_length                     
! No of longtitudes
        INTEGER, INTENT(IN) :: p_levels                       
! No of levels
        INTEGER, INTENT(IN) :: conv_base_lev(row_length,rows) 
! Level of convective base
        INTEGER, INTENT(IN) :: conv_top_lev(row_length,rows)  
! Level of convective top
        INTEGER, INTENT(IN) :: mask(row_length, rows)         
! Land/sea mask (1/0)
                                                                        
        REAL, INTENT(IN) :: delta_lambda
        REAL, INTENT(IN) :: delta_phi
        REAL, INTENT(IN) :: lat(row_length, rows)                       
        REAL, INTENT(IN) :: area(row_length, rows)
        REAL, INTENT(IN) :: r_theta_levels(row_length, rows, 0:p_levels)
        REAL, INTENT(IN) :: r_rho_levels(row_length, rows, p_levels)    
        REAL, INTENT(IN) :: p_centres(row_length, rows, p_levels)       
        REAL, INTENT(IN) :: p_boundaries(row_length, rows, 0:p_levels)  
                                                                        
        REAL, INTENT(OUT) :: an(row_length, rows, p_levels)  ! kg(N)/gridbox/s
                                                                        
! Local variables                                                       

        INTEGER :: i                                     ! Loop variable
        INTEGER :: j                                     ! Loop variable
        INTEGER :: k                                     ! Loop variable
        INTEGER :: klc                                   ! Cloud base
        INTEGER :: klt                                   ! Cloud top
        INTEGER :: asurf                                 ! Land/sea mask   

        REAL :: model_half_ht(row_length, rows, 0:p_levels)            
        REAL :: height(row_length, rows, p_levels)                      
        REAL :: cloud_base(row_length, rows)                            
        REAL :: cloud_top(row_length, rows)                             
        REAL :: delp_p(row_length, rows)                                
        REAL :: mass(row_length, rows)                                  
        REAL :: anox(p_levels)                                          
        REAL :: press(p_levels)                                         
        REAL :: hkmb                                                    
        REAL :: hkmt                                                    
        REAL :: asfaera                                                 
        REAL :: adlat                                                   

        INTEGER ::errcode, ierr

!       ...number of flashes per gridcell in flashes/gridcell/s        
        REAL :: total_flash_rate_cell
        
!       ...number of flashes per gridcell cloud to ground in flashes/gridcell/s
        REAL :: cloud2ground_flash_rate_cell
        
!       ...number of flashes per gridcell cloud to cloud in flashes/gridcell/s
        REAL :: cloud2cloud_flash_rate_cell
        
!       ...total L-NOX  column flux dendity (kg(N)/m^2/s)
        REAL :: total_N_cell

!       ...number of flashes /s        
        REAL :: total_flash_rate(row_length, rows)
        
!       ...number of flashes cloud to ground /s        
        REAL :: cloud2ground_flash_rate(row_length, rows) 
        
!       ...number of flashes cloud to cloud /s        
        REAL :: cloud2cloud_flash_rate(row_length, rows)
        
!       ...total N (kg(N)/m^2/s) produced per cell        
        REAL :: total_N(row_length, rows)


        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle
                                                                        
!       ...Initialize variables and arrays
        total_flash_rate = 0.0
        cloud2ground_flash_rate = 0.0
        cloud2cloud_flash_rate = 0.0
        total_N_cell = 0.0
                                                                        
        an(:,:,:) = 0.0
        total_flash_rate(:,:) = 0.0
        cloud2ground_flash_rate(:,:) = 0.0
        cloud2cloud_flash_rate(:,:)  = 0.0
        total_N(:,:)  = 0.0
                                                                        
! Find heights of model half levels                                     
                                                                        
        IF (lhook) CALL dr_hook('UKCA_LIGHT_CTL',zhook_in,zhook_handle)
        DO k = 0,p_levels                                               
          model_half_ht(:,:,k) = r_theta_levels(:,:,k)                 &
                               - r_theta_levels(:,:,0)                  
        END DO                                                           
                                                                        
! Find heights of model levels                                          
                                                                        
        DO k = 1,p_levels                                               
          height(:,:,k) = r_rho_levels(:,:,k) - r_theta_levels(:,:,0)
        END DO                                                           
                                                                        
! Convective cloud base and top height in km                            
                                                                        
        cloud_base = 0.0                                                
        cloud_top  = 0.0                                                
                                                                        
        DO i = 1,rows                                                   
          DO j = 1,row_length                                           
            IF (conv_base_lev(j,i) > 0) THEN                            
              cloud_base(j,i) = height(j,i,conv_base_lev(j,i))/1000.0   
            END IF                                                       
            IF (conv_top_lev(j,i) > 0) THEN                             
              cloud_top(j,i) = height(j,i,conv_top_lev(j,i))/1000.0
            END IF                                                       
          END DO                                                         
        END DO                                                           
                                                                        
        DO i = 1,rows                                                   
          DO j = 1,row_length                                           
            hkmb      = cloud_base(j,i)                                 
            hkmt      = cloud_top(j,i)                                  
            adlat     = lat(j,i)                                        
            asurf     = mask(j,i)                                       
            klc       = conv_base_lev(j,i)                              
            klt       = conv_top_lev(j,i)                               
            asfaera   = area(j,i)             
            anox      = 0.0                                             
            press     = p_centres(j,i,:)                                
            IF (klc > 1 .AND. klt < p_levels) THEN   
!             Call UKCA_LIGHT to calculate NOx emissions (kg N/grid box/s)
! DEPENDS ON: ukca_light
              CALL UKCA_LIGHT(delta_lambda,delta_phi,press,p_levels, &
                              hkmb,hkmt,klc,klt,adlat,               &
                              asfaera,asurf,anox,                    &
                              total_flash_rate_cell,                 &
                              cloud2ground_flash_rate_cell,          &
                              cloud2cloud_flash_rate_cell,           &
                              total_N_cell)                                 

!             ...L-NOx profile here still in kg(N)/gridbox/s
              an(j,i,:)  = anox 

!             ...flash rate in flashes/gridcell/min
              total_flash_rate(j,i) = total_flash_rate_cell

!             ...cld2grd flash rate in flashes/gridcell/min
              cloud2ground_flash_rate(j,i) = cloud2ground_flash_rate_cell
              
!             ...cld2cld flash rate in flashes/gridcell/min              
              cloud2cloud_flash_rate(j,i)  = cloud2cloud_flash_rate_cell
              
!             ...L-NOx column flux density in kg(N)/m^2/s              
              total_N(j,i) = total_N_cell
              
            END IF                                                       
          END DO                                                         
        END DO                                                           

! diagnostic output

!       ...goes to STASH item 50081        
        IF (L_asad_use_chem_diags .AND. L_asad_use_light_diags_tot) &
             CALL asad_lightning_diagnostics( &
             row_length, &
             rows, &
             total_lightning_flashes, &
             total_flash_rate, &
             ierr)
        
!       ...goes to STASH item 50082        
        IF (L_asad_use_chem_diags .AND. L_asad_use_light_diags_c2g) &
             CALL asad_lightning_diagnostics( &
             row_length, &
             rows, &
             lightning_flashes_cloud2ground, &
             cloud2ground_flash_rate, &
             ierr)
        
!       ...goes to STASH item 50083        
        IF (L_asad_use_chem_diags .AND. L_asad_use_light_diags_c2c) &
             CALL asad_lightning_diagnostics( &
             row_length, &
             rows, &
             lightning_flashes_cloud2cloud, &
             cloud2cloud_flash_rate, &
             ierr)
        
!       ...goes to STASH item 50084        
        IF (L_asad_use_chem_diags .AND. L_asad_use_light_diags_N) &
             CALL asad_lightning_diagnostics( &
             row_length, &
             rows, &
             total_N_2D, &
             total_N, &
             ierr)
                                                                         
! Convert NOx emissions to fluxes in mass mixing ratio, i.e., kg(NO2)/kg(air)/s
                                                                        
        DO k = 1,p_levels                                               
          delp_p = p_boundaries(:,:,k-1) - p_boundaries(:,:,k)          

!         Calculate mass of each grid box                               
          mass = area*delp_p/g                

!         convert L-NOx emissions from kg(N)/gridcell/s to kg(NO2)/kg(air)/s
          an(:,:,k) = an(:,:,k)*c_no2/c_n/mass
        END DO                                                           
                                                                       
        IF (lhook) CALL dr_hook('UKCA_LIGHT_CTL',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE UKCA_LIGHT_CTL                                   
