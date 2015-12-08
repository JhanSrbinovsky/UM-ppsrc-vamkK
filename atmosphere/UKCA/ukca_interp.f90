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
! Purpose: Subroutine to interpolate 2D field onto 3D grid      
!          as zonal mean.             
!          Based on routine from Cambridge TOMCAT model.         
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                          
!          Called from UKCA_STRATF.       
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
!                                                                      
        SUBROUTINE UKCA_INTERP(f2d,pl,lon,lat,lev,p_field,sinlat,f3d)  

        USE conversions_mod, ONLY: pi_over_180

        USE ukca_phot2d
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE Control_Max_Sizes
        IMPLICIT NONE        
                                                                       
        INTEGER, INTENT(IN) :: p_field           ! No of latitudes     
        INTEGER, INTENT(IN) :: lon               ! No of longitudes  
        INTEGER, INTENT(IN) :: lat               ! No of latitudes  
        INTEGER, INTENT(IN) :: lev               ! No of levels   
                                                                       
        REAL, INTENT(IN) :: f2d(nolat,nolev)     ! 2D data field 
        REAL, INTENT(IN) :: pl(lon,lat,lev)      ! 3D pressures    
        REAL, INTENT(IN) :: sinlat(p_field)      ! Sine (latitude)    
                                                                       
        REAL, INTENT(OUT) :: f3d(lon,lat,lev)    ! 3D field output    
                                                 ! from interpolation  
! Local variables    
                                                                       
        INTEGER :: i                             ! Loop variable 
        INTEGER :: j                             ! Loop variable 
        INTEGER :: k                             ! Loop variable  
        INTEGER :: l                             ! Loop variable 
                                                                       
        REAL :: lat2d(nolat)                     ! 2D latitudes 
        REAL :: p2d(nolev)    
        REAL :: press  
        REAL :: lati  
        REAL :: delphi  
        REAL :: wks1(lat,nolev)                  ! Working array   
        REAL :: wks2(nolat)                      ! Working array 
        REAL :: wks3(nolev)                      ! Working array  
        REAL :: wks4(lat,lev)                    ! Working array
        REAL :: UKCA_FLUPJ                       ! Fn to do interpolation

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

                                                                       
!       Set up 2D latitudes. lat=1 pole nord   
!       LAT2D() is the latitude in the centre of the 2D box  
                                                                       
        IF (lhook) CALL dr_hook('UKCA_INTERP',zhook_in,zhook_handle)

        delphi = 180.0/(2.0*nolat)  
                                                                       
        DO i = 1,nolat   
          lat2d(i)=sin((90.0 - (2*i-1)*delphi)*Pi_Over_180) 
        ENDDO 
                                                                        
!       Interpolate latitude (Linear in SIN(LAT))          
                                                                       
        DO j = 1,nolev   
          DO i = 1,nolat 
            wks2(i) = f2d(i,j)    
          ENDDO  
          DO k = 1,lat  
            lati = sinlat((k-1)*lon+1)   
            lati = MAX(lati, lat2d(nolat))  
            lati = MIN(lati, lat2d(    1))  
!DEPENDS ON: ukca_flupj
            wks1(k,j) = UKCA_FLUPJ(lati,lat2d,wks2,nolat)  
            wks1(k,j) = MAX(wks1(k,j), 0.0)  
          ENDDO      
        ENDDO 
             
!       Interpolate vertically     
           
        DO i = 1,lon
          DO k = 1,lat  
            DO j = 1,nolev    
              wks3(j) = wks1(k,j)   
              p2d(j)  = LOG(pr2d(j))    
            ENDDO      
            DO l = 1,lev    
              press = log(pl(i,k,l))   
              press = MAX(press, p2d(nolev))   
              press = MIN(press, p2d(1)) 
!DEPENDS ON: ukca_flupj
              wks4(k,l) = UKCA_FLUPJ(press,p2d,wks3,nolev)  
              wks4(k,l) = MAX(wks4(k,l), 0.0)  
            ENDDO   
          ENDDO   

          DO l = 1,lev     
            DO k = 1,lat    
              f3d(i,k,l) = wks4(k,l)   
            ENDDO      
          ENDDO    
          
        ENDDO    

        IF (lhook) CALL dr_hook('UKCA_INTERP',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE UKCA_INTERP      
