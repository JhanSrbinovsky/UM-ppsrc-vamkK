! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Routine to calculate ozone using cariolle scheme      
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Free Tracers
      
      SUBROUTINE cariolle_o3_psc( ozone_tracer,                      &
                                  o3_prod_loss,   o3_p_l_vmr,         &
                                  o3_vmr,         o3_p_l_temp,        &
                                  o3_temp,        o3_p_l_colo3,       &
                                  o3_colo3,                           &
                                  theta,                              &
                                  p_theta_levels,                     &
                                  off_x,off_y, theta_off_size,        &
                                  rows,row_length,                    &
                                  dt,                                 &
                                  exner_theta_levels,                 &
                                  model_levels)       
                                                            
! Purpose:
!        This routine evaluates the cariolle ozone tracer   
!
! Method:    on input - the actual ozone distribution       
!            on output - the updated ozone distribution   
! p_theta_levels   : pressure on theta levels.                        
! row_length,rows  : number of elements on 1 pe (x and y dimensions)   
!                    (all on one pressure level)                            
! k                : level index, counting from 1 to model_levels            
!                    1        -> surface (997 hPa)                      
!                    model_levels->top(e.g. model_levels=64 -> 0.01hPa)
! dt               : local advection timestep                   
!                    DYN_TIMESTEP/NSWEEPS 
!                    (NSWEEPS=2 if max wind or divergence 
!                    check is switched on and any values greater than
!                    the specified values occur somewhere in the model)
                                                            

        USE atm_fields_bounds_mod, ONLY: tdims_s

        USE yomhook, ONLY: lhook, dr_hook
        USE parkind1, ONLY: jprb, jpim
        USE PrintStatus_mod
        IMPLICIT NONE                                    
                                                     

        INTEGER :: k , rows, row_length        
        INTEGER :: off_x,off_y,theta_off_size                           
        INTEGER :: model_levels                 

        REAL, PARAMETER :: mrconvert =1.6551                         

! The ozone fields

        REAL :: ozone_tracer(tdims_s%i_start:tdims_s%i_end,               &
                             tdims_s%j_start:tdims_s%j_end,               &
                             tdims_s%k_start:tdims_s%k_end)
                             

        REAL :: o3_prod_loss(rows,model_levels)         
        REAL :: o3_p_l_vmr  (rows,model_levels)         
        REAL :: o3_vmr      (rows,model_levels)           
        REAL :: o3_p_l_temp (rows,model_levels)        
        REAL :: o3_temp     (rows,model_levels)         
        REAL :: o3_p_l_colo3(rows,model_levels)        
        REAL :: o3_colo3    (rows,model_levels)

        REAL :: o3_prod_loss_3D(row_length,rows,model_levels)     
        REAL :: o3_p_l_vmr_3D  (row_length,rows,model_levels)       
        REAL :: o3_vmr_3D      (row_length,rows,model_levels)        
        REAL :: o3_p_l_temp_3D (row_length,rows,model_levels)      
        REAL :: o3_temp_3D     (row_length,rows,model_levels)     
        REAL :: o3_p_l_colo3_3D(row_length,rows,model_levels)        
        REAL :: o3_colo3_3D    (row_length,rows,model_levels)


        REAL :: theta(tdims_s%i_start:tdims_s%i_end,                      &
                      tdims_s%j_start:tdims_s%j_end,                      &
                      tdims_s%k_start:tdims_s%k_end)

        REAL :: dt 
        
        
        REAL :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,             &
                               tdims_s%j_start:tdims_s%j_end,             &
                               tdims_s%k_start:tdims_s%k_end)

        REAL :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,         &
                                   tdims_s%j_start:tdims_s%j_end,         &
                                   tdims_s%k_start:tdims_s%k_end)        
                                                                        
                                            
! defining auxiliary variables              
! temp : temperature (derived from theta and pstar)           
! lat  : latitudes for the cariolle parameters           
! dso3 : column integrated o3 differences             
! cosz : cos(solar zenith angle)                      
! dp   : weight for the dso3 integration   
                         
        INTEGER :: i, j  
        REAL    :: temp                                         
        REAL    :: dp, level_dso3, to3     
        REAL    :: dso3(row_length,rows)

!       A Geer: need better way? Maybe not trapezium after all
        REAL :: dtracer_abv(row_length,rows)                   
                                              
        INTEGER :: my_pe                                         

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

       IF (lhook) CALL dr_hook('CARIOLLE_O3_PSC',zhook_in,zhook_handle)

       IF (PrintStatus  >=  PrStatus_Normal) THEN                       
          WRITE(6,*) 'Calculating Cariolle Ozone tracer'
       END IF

! Expand the cariolle parameters to 3-dimensions
! following loops ok due to using tdims fields.

       DO k = 1, model_levels
         DO j = 1, rows
           DO i = 1, row_length
              
              o3_prod_loss_3D(i,j,k)   = o3_prod_loss(j,k)
              o3_p_l_vmr_3D(i,j,k)     = o3_p_l_vmr(j,k)
              o3_vmr_3D(i,j,k)         = o3_vmr(j,k)
              o3_p_l_temp_3D(i,j,k)    = o3_p_l_temp(j,k)
              o3_temp_3D(i,j,k)        = o3_temp(j,k)
              o3_p_l_colo3_3D(i,j,k)   = o3_p_l_colo3(j,k)
              o3_colo3_3D(i,j,k)       = o3_colo3(j,k)

           END DO
         END DO
       END DO


! loop over all field elements                          
       DO k=model_levels,1,-1   ! Loop going down from model levels to 1        
         DO j=1,rows                      
           DO i=1, row_length      
! convert mixing ratio (in kgkg-1) to vmr               
            ozone_tracer(i,j,k) = ozone_tracer(i,j,k)/mrconvert   

! convert theta to temp                                  
            temp=theta(i,j,k) * exner_theta_levels(i,j,k)             
                                                              
! Calculate total column ozone differences above this level          
            IF (k == model_levels) THEN                            
! reset dso3 (vertical column o3 differences above)                
! for k=model_levels (top of model)                             
              dso3(i,j)=0                              
            ELSE                                       
! trapezium rule integration                        
              dp=p_theta_levels(i,j,k)-p_theta_levels(i,j,k+1)        
              to3=(ozone_tracer(i,j,k)-o3_vmr_3D(i,j,k))+dtracer_abv(i,j)
              level_dso3 = (dp*to3/2.)*2.152e25/101300.             
              dso3(i,j)=dso3(i,j)+level_dso3                      
            END IF                                                
                                                   
! A. Geer hack to store tracer values for trapezium rule integration    
! because the layer above has been modified by this scheme already     
! and has been converted back to kgkg-1                                 

            dtracer_abv(i,j) = ozone_tracer(i,j,k)-o3_vmr_3D(i,j,k)     

            IF (k==1) THEN

              ozone_tracer(i,j,k) = ozone_tracer(i,j,k) +              &
                                 (dt*o3_prod_loss_3D(i,j,k)) +         &
                                 (dt*(o3_p_l_temp_3D(i,j,k)*           &
                                     (temp-o3_temp_3D(i,j,k)))) 
         
            ELSE
              ozone_tracer(i,j,k) = ozone_tracer(i,j,k) +              &
                                 (dt*o3_prod_loss_3D(i,j,k)) -         &
                                 (dt*(o3_vmr_3D(i,j,k)*                &
                                          o3_p_l_vmr_3D(i,j,k))) +     &
                                 (dt*(o3_p_l_temp_3D(i,j,k)*           &
                                     (temp-o3_temp_3D(i,j,k)))) +      &
                                 (dt*(o3_p_l_colo3_3D(i,j,k)*dso3(i,j)))
          
            END IF

       ozone_tracer(i,j,k)=ozone_tracer(i,j,k)/(1.-(dt*o3_p_l_vmr_3D(i,j,k)))

! sometimes negative values over Antarctica near the ground               
            IF (ozone_tracer(i,j,k) <= 0) ozone_tracer(i,j,k)=0        
                                                                      
! convert mixing ratio (in ppmv) back to ppm                          
             ozone_tracer(i,j,k) = ozone_tracer(i,j,k)*mrconvert        
           END DO                                                     
         END DO                                                  
       END DO                                                   

       IF (lhook) CALL dr_hook('CARIOLLE_O3_PSC',zhook_out,zhook_handle)
       RETURN
 
END SUBROUTINE cariolle_o3_psc
