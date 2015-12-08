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
!  Purpose: To calculate dry deposition rates for use in UKCA chemistry.
!           Original version from Cambridge TOMCAT model.               
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
!           Called from UKCA_CHEMISTRY_CTL.                             
!                                                                       
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
        SUBROUTINE UKCA_DDEPRT(dayl, tloc, p_fieldda, dzl, bl_levels,   & 
                      z0m, u_s, temp,                                   & 
                      sinlat, i_mon,                                    & 
                      first_point, last_point,                          & 
                      dryrt)                                          

        USE atmos_constants_mod, ONLY: vkman


        USE conversions_mod, ONLY: recip_pi_over_180
        USE ASAD_MOD,        ONLY: depvel
        USE ukca_option_mod, ONLY: jpdd
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE Control_Max_Sizes
        IMPLICIT NONE                                                   

        INTEGER, INTENT(IN) :: p_fieldda             ! No of spatial poi
        INTEGER, INTENT(IN) :: bl_levels             ! No of boundary la
        INTEGER, INTENT(IN) :: i_mon                 ! Month number     
        INTEGER, INTENT(IN) :: first_point           ! first index for s
        INTEGER, INTENT(IN) :: last_point            ! last index for sp
                                                                        
        REAL, INTENT(IN) :: dayl(p_fieldda)          ! Day length       
        REAL, INTENT(IN) :: tloc(p_fieldda)          ! Local time       
        REAL, INTENT(IN) :: sinlat(p_fieldda)        ! Sine(latitude)   
        REAL, INTENT(IN) :: dzl(p_fieldda,bl_levels) ! Height of lowest 
        REAL, INTENT(IN) :: z0m(p_fieldda)           ! Roughness length 
        REAL, INTENT(IN) :: u_s(p_fieldda)           ! Surface friction 
        REAL, INTENT(IN) :: temp(p_fieldda)          ! Surface temperatu
                                                                        
        REAL, INTENT(OUT) :: dryrt(p_fieldda,jpdd)   ! Dry deposition lo
                                                                        
!       Local variables                                                 
                                                                        
        INTEGER :: i                                 ! Loop variable    
        INTEGER :: j                                 ! Loop variable    
        INTEGER :: ns                                ! Loop variable    
        INTEGER :: iday                              ! Integer indicatin
        INTEGER :: isum                              ! Integer indicatin
        INTEGER :: icat                              ! Integer indicatin
                                                                        
        REAL, PARAMETER :: midday = 12.0             ! Midday time      
        REAL, PARAMETER :: temp_summer = 275.0       ! Temp threshold fo
        REAL, PARAMETER :: z0m_min = 1.0e-3          ! Min z0m used to g
        REAL, PARAMETER :: z0m_max = 1.0e-1          ! Max z0m used to g
                                                                        
        REAL :: dawn                                 ! Time of dawn     
        REAL :: dusk                                 ! Time of dusk     
        REAL :: timel                                ! Local time       
        REAL :: vdep1                                ! Depn velocity at 
        REAL :: vdeph                                ! Depn velocity at 
        REAL :: lat(p_fieldda)                       ! Latitudes in degr

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

                                                                        
!       Initialise dryrt to zero                                        
                                                                        
        IF (lhook) CALL dr_hook('UKCA_DDEPRT',zhook_in,zhook_handle)
        dryrt(:,:) = 0.0                                                
                                                                        
!       Calculate latitude in degrees                                   
                                                                        
        DO i = first_point,last_point                                   
          lat(i) = asin(sinlat(i))*Recip_Pi_Over_180                    
        ENDDO                                                           
                                                                        
!       Set up deposition arrays for each species                       
                                                                        
        DO ns = 1,jpdd                     ! Loop over deposited species
          DO i = first_point,last_point    ! Loop over spatial scale    
                                                                        
!           Calculate time of dawn and dusk for particular day.         
                                                                        
            dawn  = midday - (dayl(i)/2.0)                              
            dusk  = midday + (dayl(i)/2.0)                              
            timel = tloc(i)                                             
                                                                        
!           Local night-time/day-time                                   
                                                                        
            IF (timel < dawn .OR. timel > dusk) THEN                    
              iday = 0                                                  
            ELSE                                                        
              iday = 1                                                  
            ENDIF                                                       
                                                                        
!           Find out if grid point is summer/winter depending on tempera
                                                                        
            IF (temp(i) > temp_summer) THEN                             
              isum = 1                                                  
            ELSEIF (temp(i) <= temp_summer) THEN                        
              isum = 0                                                  
            ENDIF                                                       
                                                                        
!           Find out what category (water,forest,grass,desert)          
!           desert not used at present                                  
                                                                        
            IF (z0m(i) < z0m_min) THEN                                  
              icat = 1                      ! water/sea - 0.001         
                                                                        
            ELSEIF (z0m(i) > z0m_max) THEN                              
              icat = 2                      ! forests - 0.1             
                                                                        
            ELSEIF (z0m(i) <= z0m_max .OR. z0m(i) >= z0m_min) THEN      
              icat = 3                      ! all other lands, grass    
                                                                        
            ENDIF                                                       
                                                                        
            IF (lat(i) >= 0.0) THEN         ! Northern hemisphere       
                                                                        
!             Overwrite switch if there is ice cover - code from 2-D mod
                                                                        
              IF((i_mon == 12 .OR. i_mon <= 3).AND.lat(i) >= 45.0)     &   
                icat = 5                                                
              IF((i_mon == 11 .OR. i_mon == 4).AND. lat(i) >= 55.0)    &   
                icat = 5                                                
              IF((i_mon == 10 .OR. i_mon == 5 .OR. i_mon == 6)         &   
              .AND.LAT(I)>=65.0)                                       & 
                icat = 5                                                
              IF((i_mon >= 7.AND. i_mon <= 9) .AND. lat(i) >= 70.0)    &   
                icat = 5                                                
                                                                        
            ELSEIF (lat(i) < 0.0) THEN      ! Southern hemisphere       
                                                                        
              IF ((i_mon == 12.OR. i_mon <= 4) .AND. lat(i) <= -75.0)  &  
                icat = 5                                                
              IF((i_mon >= 5 .AND. i_mon <= 11) .AND. lat(i) <= -70.0) & 
                icat = 5                                                
                                                                        
            ENDIF                                                       
                                                                        
!           Select appropriate dry deposition velocity - 1m values      
!           and convert from cm/s to m/s                                
                                                                        
            IF (iday == 1 .AND. isum == 1) THEN       ! Summer day      
              vdep1 = depvel(1,icat,ns)/100.0                           
                                                                        
            ELSE IF (iday == 0 .AND. isum == 1) THEN  ! Summer night    
              vdep1 = depvel(2,icat,ns)/100.0                           
                                                                        
            ELSE IF (iday == 1 .AND. isum == 0) THEN  ! Winter day      
              vdep1 = depvel(4,icat,ns)/100.0                           
                                                                        
            ELSE IF (iday == 0 .AND. isum == 0) THEN  ! Winter night    
              vdep1 = depvel(5,icat,ns)/100.0                           
                                                                        
            ENDIF                                                       
                                                                        
!           Extrapolate to the middle of lowest model layer             
                                                                        
            vdeph = vdep1/(1.0 + vdep1*alog(dzl(i,1)/2.0)/(vkman*u_s(i)))
            dryrt(i,ns) = vdeph/dzl(i,1)                                
                                                                        
          ENDDO                     ! End of loop over spatial points   
        ENDDO                       ! End of loop over deposited species
                                                                        
        IF (lhook) CALL dr_hook('UKCA_DDEPRT',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE UKCA_DDEPRT                                      
