! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to contain constants used in UKCA
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      MODULE UKCA_CONSTANTS

USE earth_constants_mod, ONLY: g, earth_radius, two_omega

USE atmos_constants_mod, ONLY: vkman, r

USE conversions_mod,     ONLY: pi, pi_over_180, recip_pi_over_180

USE water_constants_mod, ONLY: rho_water

      IMPLICIT NONE

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
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
! C_RMOL start
!
! Description:
! This contains the molar universal gas constant (8.314 J K-1 MOL-1)
!
!
      REAL,PARAMETER:: RMOL = 8.314 ! molar gas const

! C_RMOL end


! UKCA_ MODE constants

! MODE constants already defined in the UM
! 1) Original values
!      REAL, PARAMETER :: PPI=3.14159265358979323846
!      REAL, PARAMETER :: AVC=6.022E23
!      REAL, PARAMETER :: ZBOLTZ=1.3807E-23
!      REAL, PARAMETER :: VKARMN=0.4
!      REAL, PARAMETER :: RA=287.05
!      REAL, PARAMETER :: GG=9.80665
!      REAL, PARAMETER :: RAD_E=6371229.0
!      REAL, PARAMETER :: MM_DA=AVC*ZBOLTZ/RA
!      REAL, PARAMETER :: RHOSUL=1800.0E0      ! UM is slightly different
!      REAL, PARAMETER :: MMSUL=0.098E0
!      REAL, PARAMETER :: RR=8.314
!      PARAMETER(BCONST=3.44E13)        ! Volume_mode
!      PARAMETER(NU_H2SO4=3.0)          ! Volume_mode
!      PARAMETER(CONVERT=1.0e-18)       ! Volume_mode
!      PARAMETER(RHOW=1000.0,MMW=0.018) ! Volume_mode
!      PARAMETER(RHOSUL=1800.0)`        ! Volume_mode

! 2) UM definitions of MODE constants
      REAL, PARAMETER :: PPI=pi
      REAL, PARAMETER :: AVC=Avogadro
      REAL, PARAMETER :: ZBOLTZ=Boltzmann
      REAL, PARAMETER :: VKARMN=vKman
      REAL, PARAMETER :: RA=R
      REAL, PARAMETER :: RR=rmol
      REAL, PARAMETER :: GG=g
      REAL, PARAMETER :: RAD_E=Earth_Radius
      REAL, PARAMETER :: RHOSUL=rho_so4       ! 1769 or 1800 ?
      REAL, PARAMETER :: RHOW=rho_water       ! Water density (kg/m^3)

! 3) MODE/aerosol chemistry constants not already in UM
      REAL, PARAMETER :: MM_DA=AVC*ZBOLTZ/RA !
      REAL, PARAMETER :: NMOL=1.0E2          !
      REAL, PARAMETER :: TDAYS=0.0           !
      REAL, PARAMETER :: EMS_EPS=1.0E-8      !
      REAL, PARAMETER :: CONC_EPS=1.0E-8     !
      REAL, PARAMETER :: DN_EPS=1.0E-8       !
      REAL, PARAMETER :: BCONST=3.44E13      ! Volume_mode
      REAL, PARAMETER :: NU_H2SO4=3.0        ! Volume_mode
      REAL, PARAMETER :: CONVERT=1.0E-18     ! Volume_mode
      REAL, PARAMETER :: H_plus=1.0E-5       ! cloud/rain H+ concentration

! molecular masses in kg/mol

      LOGICAL, PARAMETER  ::  L_UKCA_DIURNAL_ISOPEMS   = .TRUE.
      REAL, PARAMETER :: m_air        = 0.02897       ! Air
      REAL, PARAMETER :: mmsul        = 0.09808       ! H2SO4
      REAL, PARAMETER :: mmw          = 0.0180154     ! H2O


!  -------------------------------------------------------------------
!  Conversion factor from vmr to mmr for each species.
!             vmr*c_species = mmr
!             c_species = m_species/m_air  (m_air = 28.97)
!
!  Followed by molecular masses for each species (e.g. m_ch4) in g/mol
!  -------------------------------------------------------------------

      REAL, PARAMETER :: C_O3P        = 0.5523
      REAL, PARAMETER :: C_O1D        = 0.5523
      REAL, PARAMETER :: C_O3         = 1.657 
      REAL, PARAMETER :: C_NO         = 1.036                           
      REAL, PARAMETER :: C_NO3        = 2.140 
      REAL, PARAMETER :: C_NO2        = 1.588                         
      REAL, PARAMETER :: C_N2O5       = 3.728 
      REAL, PARAMETER :: C_HO2NO2     = 2.727                     
      REAL, PARAMETER :: C_HONO2      = 2.175 
      REAL, PARAMETER :: C_HNO3       = C_HONO2   ! used in nitrate.F90
      REAL, PARAMETER :: C_OH         = 0.5868                       
      REAL, PARAMETER :: C_HO2        = 1.139 
      REAL, PARAMETER :: C_H2         = 0.06904                        
      REAL, PARAMETER :: C_H2O2       = 1.174 
      REAL, PARAMETER :: C_CH4        = 0.5523                       
      REAL, PARAMETER :: C_C          = 0.4142 
      REAL, PARAMETER :: C_CO         = 0.9665 
      REAL, PARAMETER :: C_CO2        = 1.5188          
      REAL, PARAMETER :: C_HCHO       = 1.036 
      REAL, PARAMETER :: C_MeOO       = 1.622                       
      REAL, PARAMETER :: C_H2O        = 0.6213 
      REAL, PARAMETER :: C_H2OS       = 0.6213 
      REAL, PARAMETER :: C_MeOOH      = 1.657                      
      REAL, PARAMETER :: C_HONO       = 1.622 
      REAL, PARAMETER :: C_O2         = 1.105                         
      REAL, PARAMETER :: C_N2         = 0.9665 
      REAL, PARAMETER :: C_C2H6       = 1.036                        
      REAL, PARAMETER :: C_EtOO       = 2.106 
      REAL, PARAMETER :: C_EtOOH      = 2.140                      
      REAL, PARAMETER :: C_MeCHO      = 1.519
      REAL, PARAMETER :: C_TOTH       = 1.000                        
      REAL, PARAMETER :: C_MeCO3      = 2.589 
      REAL, PARAMETER :: C_PAN        = 4.177                       
      REAL, PARAMETER :: C_C3H8       = 1.519 
      REAL, PARAMETER :: C_PrOO       = 2.589                       
      REAL, PARAMETER :: C_PrOOH      = 2.623 
      REAL, PARAMETER :: C_EtCHO      = 2.002                     
      REAL, PARAMETER :: C_EtCO3      = 3.072 
      REAL, PARAMETER :: C_Me2CO      = 2.002                     
      REAL, PARAMETER :: C_MeCOCH2OO  = 3.072 
      REAL, PARAMETER :: C_MeCOCH2OOH = 3.107            
      REAL, PARAMETER :: C_PPAN       = 4.660 
      REAL, PARAMETER :: C_MeONO2     = 2.658                     
      REAL, PARAMETER :: C_N          = 0.48325 
      REAL, PARAMETER :: C_H          = 0.03452                         
      REAL, PARAMETER :: C_N2O        = 1.5188 
      REAL, PARAMETER :: C_CFCl3      = 4.7480                     
      REAL, PARAMETER :: C_CF2Cl2     = 4.1783 
      REAL, PARAMETER :: C_ClO        = 1.7784                    
      REAL, PARAMETER :: C_HCl        = 1.2604
      REAL, PARAMETER :: C_ClONO2     = 3.3668                    
      REAL, PARAMETER :: C_HOCl       = 1.8129
      REAL, PARAMETER :: C_OClO       = 2.3309                     
      REAL, PARAMETER :: C_BrO        = 3.315
      REAL, PARAMETER :: C_BrONO2     = 4.9034                     
      REAL, PARAMETER :: C_HBr        = 2.7970
      REAL, PARAMETER :: C_HOBr       = 3.3495                      
      REAL, PARAMETER :: C_BrCl       = 3.9884
      REAL, PARAMETER :: C_MeBr       = 3.2805                     
      REAL, PARAMETER :: C_SO2        = 2.2112 
      REAL, PARAMETER :: C_SO3        = 2.7615                      
      REAL, PARAMETER :: C_Me2S       = 2.145
      REAL, PARAMETER :: C_DMS        = 2.145
      REAL, PARAMETER :: C_DMSO       = 2.6965
      REAL, PARAMETER :: C_OCS        = 2.0711                      
      REAL, PARAMETER :: C_COS        = 2.0711                      
      REAL, PARAMETER :: C_H2S        = 1.1766
      REAL, PARAMETER :: C_CS2        = 2.6282
      REAL, PARAMETER :: C_SAD        = 4.1255
      REAL, PARAMETER :: C_MSA        = 3.317
      REAL, PARAMETER :: C_S          = 1.1046                           
      REAL, PARAMETER :: C_H2SO4      = 3.385
      REAL, PARAMETER :: C_CF2ClCFCl2 = 6.4722              
      REAL, PARAMETER :: C_CHF2Cl     = 2.9858
      REAL, PARAMETER :: C_MeCCl3     = 4.6082                 
      REAL, PARAMETER :: C_CCl4       = 5.3158
      REAL, PARAMETER :: C_MeCl       = 1.7432                     
      REAL, PARAMETER :: C_CF2ClBr    = 5.7128
      REAL, PARAMETER :: C_CF3Br      = 5.1432                 
      REAL, PARAMETER :: C_Cl         = 1.2261 
      REAL, PARAMETER :: C_Cl2O2      = 3.5568
      REAL, PARAMETER :: C_Br         = 2.7627
      REAL, PARAMETER :: C_CH2Br2     = 6.0013
      REAL, PARAMETER :: c_mecf2cl    = 3.4673
      REAL, PARAMETER :: c_cf2br2     = 7.2489
      REAL, PARAMETER :: c_cf2brcf2br = 8.9748
      REAL, PARAMETER :: c_cf2clcf3   = 5.3314
      REAL, PARAMETER :: c_cf2clcf2cl = 5.8992
      REAL, PARAMETER :: c_mecfcl2    = 4.0352
      REAL, PARAMETER :: C_C5H8       = 2.3473 
      REAL, PARAMETER :: C_ISO2       = 4.0387       
      REAL, PARAMETER :: C_ISOOH      = 4.0732 
      REAL, PARAMETER :: C_ISON       = 5.3504  ! Might be revised to 5.0052
!                       for RAQ chem where ISON = (NO3)C4H6CHO: C5H7NO4: 145
      REAL, PARAMETER :: C_MACR       = 2.4163 
      REAL, PARAMETER :: C_MACRO2     = 4.1077       
      REAL, PARAMETER :: C_MACROOH    = 4.1422 
      REAL, PARAMETER :: C_MPAN       = 5.0742       
      REAL, PARAMETER :: C_HACET      = 2.5544 
      REAL, PARAMETER :: C_MGLY       = 2.4853       
      REAL, PARAMETER :: C_NALD       = 3.6244 
      REAL, PARAMETER :: C_HCOOH      = 1.5878       
      REAL, PARAMETER :: C_MECO3H     = 2.6234 
      REAL, PARAMETER :: C_MECO2H     = 2.0711      
      REAL, PARAMETER :: C_NH3        = 0.5879
      REAL, PARAMETER :: C_MONOTERP   = 4.7034
      REAL, PARAMETER :: C_SEC_ORG    = 5.1782    ! Molecular weight=150.


!     Extra species for RAQ chemistry
      REAL, PARAMETER :: c_c4h10      = 2.0021
      REAL, PARAMETER :: c_c2h4       = 0.9665
      REAL, PARAMETER :: c_c3h6       = 1.4498
      REAL, PARAMETER :: c_rnc2h4     = 3.6244
      REAL, PARAMETER :: c_rnc3h6     = 4.1077
!     rnc2h4 & rnc3h6 are CH2(NO3)CHO & CH3CH(NO3)CHO
      REAL, PARAMETER :: c_ch3oh      = 1.1046
      REAL, PARAMETER :: c_meoh       = 1.1046
      REAL, PARAMETER :: c_toluene    = 3.1757
      REAL, PARAMETER :: c_oxylene    = 3.6590
      REAL, PARAMETER :: c_memald     = 3.3828
!     MEMALDIAL is CH3-CO-CH=CH-CHO, ring fragmentation
!     product from degradation of toluene and o-xylene
      REAL, PARAMETER :: c_buooh      = 3.1067
      REAL, PARAMETER :: c_mek        = 2.4853
      REAL, PARAMETER :: c_mvk        = 2.4163
      REAL, PARAMETER :: c_mvkooh     = 4.1422
      REAL, PARAMETER :: c_gly        = 2.0021  
      REAL, PARAMETER :: c_rnc5h8     = 5.0052
      REAL, PARAMETER :: c_orgnit     = 5.5230
      REAL, PARAMETER :: c_buoo       = 3.0721
      REAL, PARAMETER :: c_meko2      = 3.5554
      REAL, PARAMETER :: c_hoc2h4o2   = 2.6579
      REAL, PARAMETER :: c_hoc3h6o2   = 3.1412
      REAL, PARAMETER :: c_hoipo2     = 4.0387
      REAL, PARAMETER :: c_tolp1      = 3.7280
      REAL, PARAMETER :: c_oxyl1      = 4.2113
      REAL, PARAMETER :: c_memald1    = 3.9351
      REAL, PARAMETER :: c_homvko2    = 4.1077

!     Extra species for EXTTC chemistry
      REAL, PARAMETER :: C_ISOO       = 4.0387       
      REAL, PARAMETER :: C_MACROO     = 4.1077       
      REAL, PARAMETER :: C_APIN       = 4.9645
      REAL, PARAMETER :: C_PropeOO    = 2.5198
      REAL, PARAMETER :: C_PropeOOH   = 2.5544
      REAL, PARAMETER :: C_OnitU      = 3.5554
      REAL, PARAMETER :: C_MEKOO      = 3.5554
      REAL, PARAMETER :: C_MEKOOH     = 3.5889
      REAL, PARAMETER :: C_EteOO      = 2.0366
      REAL, PARAMETER :: C_ALKA       = 2.0021
      REAL, PARAMETER :: C_ALKAOO     = 3.0721
      REAL, PARAMETER :: C_ALKAOOH    = 3.1067
      REAL, PARAMETER :: C_AROM       = 3.4173
      REAL, PARAMETER :: C_AROMOO     = 4.4874
      REAL, PARAMETER :: C_AROMOOH    = 4.5219
      REAL, PARAMETER :: C_BSVOC1     = 4.9645   ! as APIN
      REAL, PARAMETER :: C_BSVOC2     = 4.9645   ! as APIN
      REAL, PARAMETER :: C_B2ndry     = 4.9645   ! as APIN
      REAL, PARAMETER :: C_ASVOC1     = 3.4173   ! as AROM
      REAL, PARAMETER :: C_ASVOC2     = 3.4173   ! as AROM
      REAL, PARAMETER :: C_A2ndry     = 3.4173   ! as AROM
      REAL, PARAMETER :: C_BSOA       = 5.1778   ! 150.0
      REAL, PARAMETER :: C_ASOA       = 5.1778   ! 150.0
      REAL, PARAMETER :: C_ISOSVOC1   = 2.3473   ! as C5H8
      REAL, PARAMETER :: C_ISOSVOC2   = 2.3473   ! as C5H8
      REAL, PARAMETER :: C_ISOSOA     = 4.4874   ! 130.0

!     molecular masses in g/mol of emitted species, 
!     for budget calculations

      REAL, PARAMETER :: m_ho2     =  33.007
      REAL, PARAMETER :: m_ch4     =  16. 
      REAL, PARAMETER :: m_co      =  28.
      REAL, PARAMETER :: m_hcho    =  30. 
      REAL, PARAMETER :: m_c2h6    =  30.
      REAL, PARAMETER :: m_c3h8    =  44. 
      REAL, PARAMETER :: m_mecho   =  44.
      REAL, PARAMETER :: m_no2     =  46. 
      REAL, PARAMETER :: m_n2o5    = 108.01
      REAL, PARAMETER :: m_me2co   =  58.
      REAL, PARAMETER :: m_isop    =  68. 
      REAL, PARAMETER :: m_no      =  30.
      REAL, PARAMETER :: m_c       =  12.
      REAL, PARAMETER :: m_monoterp =  136.24

!     molecular masses of stratospheric species, for which surface 
!     mmrs are prescribed

      REAL, PARAMETER :: m_hcl        =  36.5
      REAL, PARAMETER :: m_n2o        =  44.
      REAL, PARAMETER :: m_clo        =  51.5
      REAL, PARAMETER :: m_hocl       =  52.5
      REAL, PARAMETER :: m_oclo       =  67.5
      REAL, PARAMETER :: m_clono2     =  97.5
      REAL, PARAMETER :: m_cf2cl2     = 121. 
      REAL, PARAMETER :: m_cfcl3      = 137.5
      REAL, PARAMETER :: m_hbr        =  81. 
      REAL, PARAMETER :: m_mebr       =  95.
      REAL, PARAMETER :: m_bro        =  96. 
      REAL, PARAMETER :: m_hobr       =  97.
      REAL, PARAMETER :: m_brcl       = 115.5
      REAL, PARAMETER :: m_brono2     = 142.
      REAL, PARAMETER :: m_cf2clcfcl2 = 187.5
      REAL, PARAMETER :: m_chf2cl     =  86.5
      REAL, PARAMETER :: m_meccl3     = 133.5
      REAL, PARAMETER :: m_ccl4       = 154.
      REAL, PARAMETER :: m_mecl       =  50.5
      REAL, PARAMETER :: m_cf2clbr    = 165.5
      REAL, PARAMETER :: m_cf3br      = 149.
      REAL, PARAMETER :: m_ch2br2     = 173.835

! sulphur containing, etc.
      REAL, PARAMETER :: m_ocs        =  60. 
      REAL, PARAMETER :: m_cos        =  60. 
      REAL, PARAMETER :: m_h2s        =  34.086
      REAL, PARAMETER :: m_cs2        =  76.14
      REAL, PARAMETER :: m_dms        =  62.1 
      REAL, PARAMETER :: m_dmso       =  78.13
      REAL, PARAMETER :: m_me2s       =  62.1
      REAL, PARAMETER :: m_msa        =  96.1
      REAL, PARAMETER :: m_sec_org    =  150.0
      REAL, PARAMETER :: m_s          =  32.07
      REAL, PARAMETER :: m_so2        =  64.06 
      REAL, PARAMETER :: m_so3        =  80.06
      REAL, PARAMETER :: m_so4        =  96.06
      REAL, PARAMETER :: m_h2so4      =  98.07
      REAL, PARAMETER :: m_nh3        =  17.03
      REAL, PARAMETER :: m_nh42so4    =  132.16

      REAL, PARAMETER :: m_cl         = 35.5 
      REAL, PARAMETER :: m_cl2o2      = 103.
      REAL, PARAMETER :: m_br         = 80.
      REAL, PARAMETER :: m_h2         = 2.016
      REAL, PARAMETER :: m_h2o        = 18.0
      REAL, PARAMETER :: m_mecoch2ooh = 90.0
      REAL, PARAMETER :: m_isooh      = 118.0
      REAL, PARAMETER :: m_mpan       = 147.0
      REAL, PARAMETER :: m_ppan       = 135.0         
      
!     Extra masses for RAQ or other chemistries
      REAL, PARAMETER :: m_c5h8    =  68.0
      REAL, PARAMETER :: m_c4h10   =  58.0
      REAL, PARAMETER :: m_c2h4    =  28.0
      REAL, PARAMETER :: m_c3h6    =  42.0
      REAL, PARAMETER :: m_toluene =  92.0
      REAL, PARAMETER :: m_oxylene = 106.0
      REAL, PARAMETER :: m_ch3oh   =  32.0
      REAL, PARAMETER :: m_meoh    =  32.0
      REAL, PARAMETER :: m_buooh   =  90.0
      REAL, PARAMETER :: m_mvkooh  = 120.0   
      REAL, PARAMETER :: m_orgnit  = 160.0
      REAL, PARAMETER :: m_macrooh = 120.0

      REAL, PARAMETER :: m_hono = 47.0
!     Extra masses for Wesely scheme 
      REAL, PARAMETER :: m_macr    = 70.0
      REAL, PARAMETER :: m_etcho   = 58.0
      REAL, PARAMETER :: m_nald    = 105.0
      REAL, PARAMETER :: m_mgly    = 72.0
      REAL, PARAMETER :: m_hacet  = 74.0


!     Extra masses for EXTTC chemistry
      REAL, PARAMETER :: m_apin     =  136.
      REAL, PARAMETER :: m_mvk      =  70.
      REAL, PARAMETER :: m_mek      =  72.
      REAL, PARAMETER :: m_alka     =  58.        ! as butane
      REAL, PARAMETER :: m_arom     =  99.        ! (toluene + xylene)/2
      REAL, PARAMETER :: m_bsvoc1   = 144.0
      REAL, PARAMETER :: m_bsvoc2   = 144.0
      REAL, PARAMETER :: m_asvoc1   = 99.0
      REAL, PARAMETER :: m_asvoc2   = 99.0
      REAL, PARAMETER :: m_isosvoc1 = 68.0
      REAL, PARAMETER :: m_isosvoc2 = 68.0
      REAL, PARAMETER :: m_onitu    = 102.0
      REAL, PARAMETER :: m_bsoa     = 150.0
      REAL, PARAMETER :: m_asoa     = 150.0
      REAL, PARAMETER :: m_isosoa   = 130.0
      REAL, PARAMETER :: m_alkaooh  = 90.0
      REAL, PARAMETER :: m_aromooh  = 130.0
      REAL, PARAMETER :: m_mekooh   = 104.0
     
!     The mass of organic nitrate is an approximation, 
!     calculated as the average of ORGNIT formed by two 
!     reacs. in UKCA_CHEMCO_RAQ:
!      NO2 + TOLP1 --> ORGNIT (A)
!      NO2 + OXYL1 --> ORGNIT (B)
!      * TOL  = methylbenzene       = C6H5(CH3)
!        OXYL = 1,2-dimethylbenzene = C6H4(CH3)2
!      * TOL  + OH --> TOLP1: C6H4(OH)(CH3)  = methyl phenol 
!        OXYL + OH --> OXYL1: C6H3(OH)(CH3)2 = dimethyl phenol
!      * ORGNIT A: TOLP1 + NO2 ~ C6H3(CH3)(OH)NO2  ~ 
!                  C7H7NO3: methyl nitrophenol   -> 153
!        ORGNIT B: OXYL1 + NO2 ~ C6H2(CH3)2(OH)NO2 ~ 
!                  C8H9NO3: dimethyl nitrophenol -> 167
!  -------------------------------------------------------------------

      END MODULE UKCA_CONSTANTS
