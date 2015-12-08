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
! Purpose: Subroutine to calculate NOx lightning emissions for one      
!          vertical column based on cloud height above surface,         
!          latitude & land/ocean surface.                               
!                                                                       
!          Based on light.F from Cambridge TOMCAT model                 
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
! Method:  Cloud height above surface and surface type yield lightning  
!          flash frequency. Latitude yields CG/CC ratio (Price & Rind   
!          1993), assuming that all dH (zero degrees to cloud top) lie  
!          between 5.5-14km. Convert flashes to NOx production for the  
!          1/2 hr dynamical period.                                     
!                                                                       
!          Called from UKCA_LIGHT_CTL.                                  
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
        SUBROUTINE UKCA_LIGHT(delta_lambda,delta_phi,ppress,niv, &
                              hkmb,hkmt,klc,klt,adlat,           &
                              asfaera,asurf,anox,                &
                              total_flash_rate,                  &
                              cloud2ground_flash_rate,           &
                              cloud2cloud_flash_rate,            &
                              total_N)                 

        USE yomhook,         ONLY: lhook, dr_hook
        USE parkind1,        ONLY: jprb, jpim
        USE parcons_mod,     ONLY: rad, deg
        USE ukca_constants,  ONLY: AVC
        IMPLICIT NONE                                                   
                                                                        
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_S,                                                       &
                                 ! relative molecular mass S kg/mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_HO2,                                                     &
                                 ! relative molecular mass HO2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_S    = 3.20E-2,                                    &
                                          ! kg/mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &           RMM_HO2  = 3.30E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------
                                                                        
        INTEGER, INTENT(IN) :: niv ! No of vertical level
        INTEGER, INTENT(IN) :: klc ! Level of cloud base 
        INTEGER, INTENT(IN) :: klt ! Level of cloud top
        INTEGER, INTENT(IN) :: asurf ! Land (1) / sea (0) m
                                                                        
        REAL, INTENT(IN) ::  delta_lambda ! gridbox width (radiants)
        REAL, INTENT(IN) ::  delta_phi ! gridbox height (radiants)
        REAL, INTENT(IN) ::  hkmt ! Height of cloud top 
        REAL, INTENT(IN) ::  hkmb ! Height of cloud base
        REAL, INTENT(IN) ::  adlat ! Latitude
        REAL, INTENT(IN) ::  ppress(niv) ! Pressures at model l
        REAL, INTENT(IN) ::  asfaera ! Surf area * (radius 
                                                                        
!       ...NOx lightning emissions        
        REAL, INTENT(OUT) :: anox(niv) ! kg(N)/gridcell/s

!       ...number of flashes in a gridcell /s        
        REAL, INTENT(OUT) :: total_flash_rate        
        
!       ...number of flashes in a gridcell cloud to ground /s        
        REAL, INTENT(OUT) :: cloud2ground_flash_rate
        
!       ...number of flashes in a gridcell cloud to cloud /s        
        REAL, INTENT(OUT) :: cloud2cloud_flash_rate
        
!       ...lighting N column density in kg(N)/m^2/s
        REAL, INTENT(OUT) :: total_N
                                                                         
! Local variables                                                       
                                                                        
        INTEGER :: jniv ! Loop variable       
        INTEGER :: k ! Loop variable       
                                                                        
!       ...Minimum cloud depth         
        REAL, PARAMETER ::  Min_clouddepth = 5.0
        
!       ...distance (km) per degree at equator        
        REAL, PARAMETER ::  km_per_deg_eq  = 111.11


!       ...NOx production parameters

!       Price et al., J.G.R., 1997 (1 and 2)
!       ...avg energy cld-2-grd flash (J)        
!        REAL, PARAMETER ::  E_cg = 6.7E+09
!       ...avg energy cld-2-cld flash (J) (eq 0.1*E_cg)        
!        REAL, PARAMETER ::  E_cc = 6.7E+08
!       ...avg NO production rate (molecs(NO) J-1)        
!        from Allen and Pickering, JGR, 2002
!        REAL, PARAMETER ::  P_no = 1.0E+17

!       from review of Schumann and Huntrieser, ACP, 2007
!       ...avg energy cld-2-grd flash (J)        
        REAL, PARAMETER ::  E_cg = 3.0E+09
!       ...avg energy cld-2-cld flash (J) 
        REAL, PARAMETER ::  E_cc = 0.9E+09
!       ...avg NO production rate (molecs(NO) J-1)        
        REAL, PARAMETER ::  P_no = 25.0E+16
        
!       ...Molecular mass of N (kg/mol)        
        REAL, PARAMETER ::  M_n  = 14.01*1.0E-03
                                                                        
        REAL :: aflash ! Flash frequency (flashes/min)  
        REAL :: adh                                                     
        REAL :: az ! Cloud-cloud/cloud-ground       
        REAL :: ap ! Cloud-ground flashes / total fl
        REAL :: acgfrac ! Cloud-ground flash frequency (f
        REAL :: accfrac ! Cloud-cloud flash frequency (fl
        REAL :: acgnox                                                  
        REAL :: accnox                                                  
        REAL :: dpcg                                                    
        REAL :: dpcc 
        
        REAL :: gb_area_30N ! area of gridbox at 30N (normalisation factor)
        REAL :: ew_res_deg  ! gridbox width in degrees
        REAL :: ns_res_deg  ! gridbox height in degrees
        REAL :: fr_calib_fac ! model resolution calibration factor

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

                                                                        
!       Initialise variables                                            
                                                                        
        IF (lhook) CALL dr_hook('UKCA_LIGHT',zhook_in,zhook_handle)

        gb_area_30N  = 0.0
        ew_res_deg   = 0.0
        ns_res_deg   = 0.0
        fr_calib_fac = 0.0

        aflash  = 0.0                                                   
        adh     = 0.0                                                   
        az      = 0.0                                                   
        ap      = 0.0                                                   
        acgfrac = 0.0                                                   
        accfrac = 0.0                                                   
        acgnox  = 0.0                                                   
        accnox  = 0.0                                                   
! Initialise to zero and will be value if minimum cloud depth not met.
        anox    = 0.0
        total_flash_rate = 0.0
        cloud2ground_flash_rate = 0.0
        cloud2cloud_flash_rate = 0.0
        total_N = 0.0

!       ...Calculate resolution-dependent calibration factor
!          spatial calibration factor (N96); 
!          cf., Price and Rind, Mon. Weather Rev., 1994 and
!          Allen and Pickering, JGR, 2002, pp. 7-8
        ew_res_deg   = delta_lambda*deg
        ns_res_deg   = delta_phi*deg
        fr_calib_fac = 0.97241*EXP(0.048203*ew_res_deg*ns_res_deg)

!       ...compute gridbox area at 30N (normalisation factor)
        gb_area_30N = (ew_res_deg*km_per_deg_eq*0.87)* & ! gridbox width (km)
                      (ns_res_deg*km_per_deg_eq)*      & ! gridbox height (km)
                      1.0E+06                            ! (km^2 --> m^2)
                                                                        
!       Set minimum cloud depth to 5 km                                 
                                                                        
        IF ((hkmt-hkmb) > Min_clouddepth) THEN                          
                                                                        
!         ...flashes per minute
          IF (asurf == 0) THEN                  ! Ocean                 
            aflash = 6.40e-04*(hkmt**1.73)      ! Ocean flash frequency (1/min)
          ELSE                                  ! Land                  
            aflash = 3.44e-05*(hkmt**4.9)       ! Land flash frequency  (1/min)
          ENDIF                                                         
                                                                        
!         ...Calculate flash rate in flashes/s/gridbox
!            note: c.f., Allen and Pickering, J.G.R., 2002

!         ...calibrate for model resolution          
          aflash = aflash*fr_calib_fac*(asfaera/gb_area_30N)
          
!         ...convert from flases/gridbox/min to flashes/gridbox/s          
          aflash = aflash/60.0
                                                                        
!         ...work out proportion of flashes that are cloud to ground       
          adh = -6.64e-05*(ABS(adlat)**2)-4.73e-03*ABS(adlat)+7.34

          az  = 0.021*(adh**4)-0.648*(adh**3)+7.493*(adh**2) &
              - 36.54*adh+63.09

          ap  = 1./(1.+az)                                              
          acgfrac = aflash*ap
          accfrac = aflash-acgfrac

          ! calculate diagnostics quantities in flashes/gridbox/min
          total_flash_rate         = aflash*60.0
          cloud2ground_flash_rate  = acgfrac*60.0
          cloud2cloud_flash_rate   = accfrac*60.0
                                                                        
!         ...compute NO production in kg(N)/gridbox/s based on flash frequency
!         ...cloud-to-ground NOx
          acgnox = (acgfrac*E_cg*P_no*M_n)/AVC

!         ...cloud-to-cloud NOx
          accnox = (accfrac*E_cc*P_no*M_n)/AVC

!         ...total lightning NOx column density in kg(N)/m^2/s
          total_N = (acgnox + accnox)/asfaera


!         Distribute over the column with each box having same vmr in an
!         multiply by CONVFAC, conversion factor that gives emissions if
!         were 100 flashes s**-1                                        
                                                                        
!         Work out which pressure is closest to 500 hPa                 
                                                                        
          LOOP: DO jniv = niv,1,-1                                      
            IF (ppress(jniv) >= 50000.0) EXIT LOOP                      
          ENDDO LOOP                                                    

!         3 from cloud base to 2 above top                              
!         KLT is the level above cloud top                              
                                                                        
          dpcg = ppress(1) - ppress(jniv)                             
          dpcc = ppress(jniv) - ppress(klt)                             
                                                                        
!         ...construct L-NOx profile in kg(N)/gridcell/s
!         ...first cloud-to-ground L-NOx profiles (kg(N)/gridcell/s)
          IF ((jniv-1) == 1) THEN                                       
            anox(1) = acgnox
          ELSE                                                          
            DO k = 1,jniv-1                                             
              anox(k) = acgnox * ((ppress(k)-ppress(k+1))/dpcg)
            END DO                                                      
          ENDIF                                                         
                                                                        
!         ...then cloud-to-cloud L-NOx profiles (kg(N)/gridcell/s) 
          IF (ppress(jniv) <= ppress(klt)) THEN                         
            anox(klt-1) = anox(klt-1) + accnox                   
          ELSE                                                          
            DO k = jniv,klt-1                                           
              anox(k) = accnox * ((ppress(k)-ppress(k+1))/dpcc)
            END DO                                                      
          ENDIF                                                         
        ENDIF                                                           
                                                                        
        IF (lhook) CALL dr_hook('UKCA_LIGHT',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE UKCA_LIGHT                                       
