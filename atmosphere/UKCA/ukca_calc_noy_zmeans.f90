! *****************************COPYRIGHT******************************* 
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT******************************* 
!                                                                       
! Purpose: Subroutine to calculate zonal means for NOy species
!                                                                       
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_STRATF
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
      SUBROUTINE UKCA_CALC_NOY_ZMEANS(lon, lat, lev,                &
     &                            glon, glat, first_row,            &
     &                            st, noyzm, nozm, no2zm, no3zm,    &
     &                            n2o5zm, hno3zm, hno4zm)
      USE asad_mod,        ONLY: speci
      USE ukca_option_mod, ONLY: jpspec
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE UM_ParVars
      USE Control_Max_Sizes
      IMPLICIT NONE  

      INTEGER, INTENT(IN) :: lon        ! No of longitudes            
      INTEGER, INTENT(IN) :: lat        ! No of latitudes             
      INTEGER, INTENT(IN) :: lev        ! No of levels      
      INTEGER, INTENT(IN) :: glon       ! No of global longitudes
      INTEGER, INTENT(IN) :: glat       ! No of global latitudes
      INTEGER, INTENT(IN) :: first_row  ! First global latitude       
                                                                              
      REAL, INTENT(IN)    :: st(lon,lat,lev,jpspec)   ! Species

      REAL, INTENT(OUT)   :: noyzm (lat,lev)
      REAL, INTENT(OUT)   :: nozm  (lat,lev)
      REAL, INTENT(OUT)   :: no2zm (lat,lev)
      REAL, INTENT(OUT)   :: no3zm (lat,lev)
      REAL, INTENT(OUT)   :: n2o5zm(lat,lev)
      REAL, INTENT(OUT)   :: hno3zm(lat,lev)
      REAL, INTENT(OUT)   :: hno4zm(lat,lev)
                                                                              
!     Local variables                                                       

      INTEGER :: info           ! Tag used in communication                  
      INTEGER :: i,ij,j,k       ! Loop variables                                                 

      REAL :: fac                         ! Factor used to calculate zmean   
      REAL :: tmp1(glat,lev)              ! Temporary array for zonal 
      REAL :: tmp2(glat,lev)              ! Temporary array for zonal mean   
      REAL :: tmp3(glat,lev)              ! Temporary array for zonal mean   
      REAL :: tmp4(glat,lev)              ! Temporary array for zonal 
      REAL :: tmp5(glat,lev)              ! Temporary array for zonal mean   
      REAL :: tmp6(glat,lev)              ! Temporary array for zonal mean   
      REAL :: tmp7(glat,lev)              ! Temporary array for zonal 

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

                                                                            
!     Initialise zonal means for NOy species 
                                                                        
      IF (lhook) CALL dr_hook('UKCA_CALC_NOY_ZMEANS',zhook_in,zhook_handle)

      noyzm  = 0.0                                                    
      hno3zm = 0.0                                                    
      hno4zm = 0.0                                                    
      nozm   = 0.0                                                    
      no2zm  = 0.0                                                    
      no3zm  = 0.0                                                    
      n2o5zm = 0.0                                                    
                                                                        
      DO k = 1,jpspec                                               
        SELECT CASE (speci(k))                                      
        CASE('N         ','ClONO2    ','HONO      ',         &              
     &   'PAN       ','PPAN      ','MeONO2    ','BrONO2    ')      
          DO j=1,lev
            DO ij=1,lat
              DO i=1,lon
                noyzm(ij,j) = noyzm(ij,j) + st(i,ij,j,k)
              ENDDO
            ENDDO
          ENDDO                              
        CASE('NO        ')                                          
          DO j=1,lev
            DO ij=1,lat
              DO i=1,lon
                noyzm(ij,j) = noyzm(ij,j) + st(i,ij,j,k)
                nozm (ij,j) = nozm (ij,j) + st(i,ij,j,k)
              ENDDO
            ENDDO
          ENDDO                                   
        CASE('NO2       ')                                          
          DO j=1,lev
            DO ij=1,lat
              DO i=1,lon
                noyzm(ij,j) = noyzm(ij,j) + st(i,ij,j,k)
                no2zm(ij,j) = no2zm(ij,j) + st(i,ij,j,k)
              ENDDO
            ENDDO
          ENDDO                                 
        CASE('NO3       ')                                          
          DO j=1,lev
            DO ij=1,lat
              DO i=1,lon
                noyzm(ij,j) = noyzm(ij,j) + st(i,ij,j,k)
                no3zm(ij,j) = no3zm(ij,j) + st(i,ij,j,k)
              ENDDO
            ENDDO
          ENDDO                                
        CASE('N2O5      ')                                          
          DO j=1,lev
            DO ij=1,lat
              DO i=1,lon
                noyzm (ij,j) = noyzm (ij,j) + 2.0*st(i,ij,j,k)
                n2o5zm(ij,j) = n2o5zm(ij,j) +     st(i,ij,j,k)
              ENDDO
            ENDDO
          ENDDO                            
        CASE('HONO2     ')                                          
          DO j=1,lev
            DO ij=1,lat
              DO i=1,lon
                noyzm (ij,j) = noyzm (ij,j) + st(i,ij,j,k)
                hno3zm(ij,j) = hno3zm(ij,j) + st(i,ij,j,k)
              ENDDO
            ENDDO
          ENDDO                                 
        CASE('HO2NO2    ')                                          
          DO j=1,lev
            DO ij=1,lat
              DO i=1,lon
                noyzm (ij,j) = noyzm (ij,j) + st(i,ij,j,k)
                hno4zm(ij,j) = hno4zm(ij,j) + st(i,ij,j,k)
              ENDDO
            ENDDO
          ENDDO                               
        END SELECT                                                  
      ENDDO                           

      fac    = 1.0/real(glon)                                         
      noyzm  = noyzm*fac                                              
      hno3zm = hno3zm*fac                                             
      hno4zm = hno4zm*fac                                             
      nozm   = nozm*fac                                               
      no2zm  = no2zm*fac                                              
      no3zm  = no3zm*fac                                              
      n2o5zm = n2o5zm*fac                                             

!     Need to sum along processor rows                                

      tmp1 = 0.0                                                      
      tmp2 = 0.0                                                      
      tmp3 = 0.0                                                      
      tmp4 = 0.0                                                      
      tmp5 = 0.0                                                      
      tmp6 = 0.0                                                      
      tmp7 = 0.0                                                      

      tmp1(first_row:first_row+lat-1,:) = noyzm                       
      tmp2(first_row:first_row+lat-1,:) = hno3zm                      
      tmp3(first_row:first_row+lat-1,:) = hno4zm                      
      tmp4(first_row:first_row+lat-1,:) = nozm                        
      tmp5(first_row:first_row+lat-1,:) = no2zm                       
      tmp6(first_row:first_row+lat-1,:) = no3zm                       
      tmp7(first_row:first_row+lat-1,:) = n2o5zm                      

      CALL GC_RSUM (glat*lev,nproc,info,tmp1)                         
      CALL GC_RSUM (glat*lev,nproc,info,tmp2)                         
      CALL GC_RSUM (glat*lev,nproc,info,tmp3)                         
      CALL GC_RSUM (glat*lev,nproc,info,tmp4)                         
      CALL GC_RSUM (glat*lev,nproc,info,tmp5)                         
      CALL GC_RSUM (glat*lev,nproc,info,tmp6)                         
      CALL GC_RSUM (glat*lev,nproc,info,tmp7)                         

      noyzm  = tmp1(first_row:first_row+lat-1,:)                      
      hno3zm = tmp2(first_row:first_row+lat-1,:)                      
      hno4zm = tmp3(first_row:first_row+lat-1,:)                      
      nozm   = tmp4(first_row:first_row+lat-1,:)                      
      no2zm  = tmp5(first_row:first_row+lat-1,:)                      
      no3zm  = tmp6(first_row:first_row+lat-1,:)                      
      n2o5zm = tmp7(first_row:first_row+lat-1,:)                      

      IF (lhook) CALL dr_hook('UKCA_CALC_NOY_ZMEANS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALC_NOY_ZMEANS
