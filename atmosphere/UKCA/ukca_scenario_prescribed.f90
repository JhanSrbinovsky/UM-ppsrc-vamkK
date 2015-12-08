! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved. 
! This routine has been licensed to the Met Office for use and 
! distribution under the UKCA collaboration agreement, subject  
! to the terms and conditions set out therein. 
! [Met Office Ref SC138]  
! *****************************COPYRIGHT*******************************
!
! Description:
!  Specify surface boundary conditions for chemical tracers
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  This subroutine specifies surface MMRs for halogenated source gases   
!  following the A1 scenario of WMO (2006). It also lumps halogen species
!  to account for total chlorine and bromine, where only a restricted
!  number of species is modelled. Trace gas MMrs are taken from the UM/UMUI.
!
!                                                                       
!                                                                       
!  Called from UKCA_MAIN1.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!                     
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
      SUBROUTINE UKCA_SCENARIO_PRESCRIBED(n_lbc_specs, lbc_specs,       &
                            lbc_mmr, perform_lumping)
                                                                        
      USE ukca_trace_gas_mixratio
      USE UKCA_CONSTANTS,    ONLY: c_cf2cl2,c_cf2clcfcl2, c_cf2clcf2cl, &
          c_cf2clcf3, c_cfcl3, c_ccl4, c_meccl3, c_chf2cl, c_mecfcl2,   &
          c_mecf2cl, c_cf2clbr, c_mecl, c_mebr, c_cf2br2, c_cf3br,      &
          c_cf2brcf2br
      USE PrintStatus_mod
      USE yomhook,           ONLY: lhook, dr_hook
      USE parkind1,          ONLY: jprb, jpim
      IMPLICIT NONE                                                     
                                                                        
      INTEGER, INTENT(IN)      :: n_lbc_specs               ! Number of lower BC species
      REAL, INTENT(INOUT)      :: lbc_mmr(n_lbc_specs)      ! Lower BC mass mixing ratios
      CHARACTER(LEN=10), INTENT(IN) :: lbc_specs(n_lbc_specs)    ! LBC species  
      LOGICAL, INTENT(IN)      :: perform_lumping           ! T for lumping of lbc's
                                                                        
! Scenario A1, table 8.5 of WMO(2006), expressed in MMR                 
                                                                        
      LOGICAL, SAVE :: first = .TRUE.                                   
                                                                        
! Cater for 16 species from WMO(2006) and 4 species, listed in order below
                                                                        
      INTEGER,SAVE :: icfcl3=0, icf2cl2=0, icf2clcfcl2=0, icf2clcf2cl=0      
      INTEGER,SAVE :: icf2clcf3=0, iccl4=0, imeccl3=0, ichf2cl=0             
      INTEGER,SAVE :: imecfcl2=0, imecf2cl=0, icf2clbr=0, icf2br2=0          
      INTEGER,SAVE :: icf3br=0, icf2brcf2br=0, imecl=0, imebr=0              
      INTEGER,SAVE :: in2o=0, ich4=0, ih2=0, ico2=0, ich2br2=0               
      INTEGER,SAVE :: icos=0


      INTEGER, PARAMETER :: n_full_spec=18                              
                                                                        
      REAL :: full_lbc(n_full_spec)                                     
                                                                        
! List which species replaces which in case of lumping. 1 = F11, 2 = F12
! 16 = CH3Br                                                            
      INTEGER, SAVE :: replace_Cl(n_full_spec) =                        &
                                        (/0,0,2,2,                      &
                                          2,1,1,2,                      &
                                          1,1,1,0,                      &
                                          0,0,1,0,                      &
                                          0,0 /)

      INTEGER, SAVE :: replace_Br(n_full_spec) =                        &
                                        (/0,0,0,0,                      &
                                          0,0,0,0,                      &
                                          0,0,16,16,                    &
                                          16,16,0,0,                    &
                                          0,0/)                           
                                                                        
      REAL, SAVE :: convfac_Cl(n_full_spec)                             
                                                                        
      REAL, SAVE :: convfac_Br(n_full_spec)                             
                                                                        
      REAL :: CH4              ! Mass mixing ratios                                   
      REAL :: N2O                                   
      REAL :: CO2                                   
      REAL :: CFC11                                        
      REAL :: CFC12                                        
      REAL :: CFC113                                       
      REAL :: CFC114                                       
      REAL :: CFC115                                       
      REAL :: CCl4                                         
      REAL :: MeCCl3                                       
      REAL :: HCFC22                                       
      REAL :: HCFC141b                                     
      REAL :: HCFC142b                                     
      REAL :: H1211                                        
      REAL :: H1202                                        
      REAL :: H1301                                        
      REAL :: H2402                                        
      REAL :: MeBr                                         
      REAL :: MeCl                                         
      REAL :: CH2Br2
      REAL :: H2 ! currently fixed here
      REAL :: COS_mmr  ! used instead of COS because of Fortran intrinsic fn
                                                                        
! local variables                                                       
      REAL :: ftime                                                     
      REAL :: frac                                                      
      INTEGER :: index                                                  
      INTEGER :: i                                                      

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

                                                                        
! Body of subroutine                                                   
      IF (lhook) CALL dr_hook('UKCA_SCENARIO_PRESCRIBED',zhook_in,zhook_handle)

      CH4      = um_ch4_for_ukca
      N2O      = um_n2o_for_ukca
      CO2      = um_co2_for_ukca
      CFC11    = um_cfc11_for_ukca
      CFC12    = um_cfc12_for_ukca
      CFC113   = um_cfc113_for_ukca
      CFC114   = um_cfc114_for_ukca
      CFC115   = um_cfc115_for_ukca
      CCl4     = um_ccl4_for_ukca
      MeCCl3   = um_meccl3_for_ukca
      HCFC22   = um_hcfc22_for_ukca
      HCFC141b = um_hcfc141b_for_ukca
      HCFC142b = um_hcfc142b_for_ukca
      H1211    = um_h1211_for_ukca
      H1202    = um_h1202_for_ukca
      H1301    = um_h1301_for_ukca
      H2402    = um_h2402_for_ukca
      MeBr     = um_mebr_for_ukca
      MeCl     = um_mecl_for_ukca
      CH2Br2   = um_ch2br2_for_ukca
      H2       = um_h2_for_ukca 
      COS_mmr  = um_cos_for_ukca

      IF (first) THEN                                                                  
       ! output values from the UM
         IF (PrintStatus >= PrStatus_Oper) THEN
            write(6,*) 'UM PRESCRIBED CH4      = ',um_ch4_for_ukca
            write(6,*) 'UM PRESCRIBED N2O      = ',um_n2o_for_ukca
            write(6,*) 'UM PRESCRIBED CO2      = ',um_co2_for_ukca
            write(6,*) 'UM PRESCRIBED CFC11    = ',um_cfc11_for_ukca
            write(6,*) 'UM PRESCRIBED CFC12    = ',um_cfc12_for_ukca
            write(6,*) 'UM PRESCRIBED CFC113   = ',um_cfc113_for_ukca
            write(6,*) 'UM PRESCRIBED CFC114   = ',um_cfc114_for_ukca
            write(6,*) 'UM PRESCRIBED CFC115   = ',um_cfc115_for_ukca
            write(6,*) 'UM PRESCRIBED CCl4     = ',um_ccl4_for_ukca
            write(6,*) 'UM PRESCRIBED MeCCl3   = ',um_meccl3_for_ukca
            write(6,*) 'UM PRESCRIBED HCFC22   = ',um_hcfc22_for_ukca
            write(6,*) 'UM PRESCRIBED HCFC141b = ',um_hcfc141b_for_ukca
            write(6,*) 'UM PRESCRIBED HCFC142b = ',um_hcfc142b_for_ukca
            write(6,*) 'UM PRESCRIBED H1211    = ',um_h1211_for_ukca
            write(6,*) 'UM PRESCRIBED H1202    = ',um_h1202_for_ukca
            write(6,*) 'UM PRESCRIBED H1301    = ',um_h1301_for_ukca
            write(6,*) 'UM PRESCRIBED H2402    = ',um_h2402_for_ukca
            write(6,*) 'UM PRESCRIBED MeBr     = ',um_mebr_for_ukca
            write(6,*) 'UM PRESCRIBED MeCl     = ',um_mecl_for_ukca
            write(6,*) 'UM PRESCRIBED CH2Br2   = ',um_ch2br2_for_ukca
            write(6,*) 'UM PRESCRIBED H2       = ',um_h2_for_ukca 
            WRITE(6,'(A,E12.3)') 'UM PRESCRIBED COS      = ',um_cos_for_ukca
         END IF
        first = .FALSE.                                                 
        convfac_Cl = (/                                                 &
          0., 0., 1.5*c_cf2cl2/c_cf2clcfcl2, c_cf2cl2/c_cf2clcf2cl,     &
          0.5*c_cf2cl2/c_cf2clcf3,                                      &
          1.3333*c_cfcl3/c_ccl4, c_cfcl3/c_meccl3,                      &
          0.5*c_cf2cl2/c_chf2cl,                                        &
          0.6667*c_cfcl3/c_mecfcl2, 0.3333*c_cfcl3/c_mecf2cl,           &
          0.3333*c_cfcl3/c_cf2clbr, 0.,                                 &
          0.,0.,0.3333*c_cfcl3/c_mecl,0.,0.,0./)                           
                                                                        
        convfac_Br = (/                                                 &
          0., 0., 0., 0.,                                               &
          0., 0., 0., 0.,                                               &
          0., 0., c_mebr/c_cf2clbr, 2.*c_mebr/c_cf2br2,                 &
          c_mebr/c_cf3br,2.*c_mebr/c_cf2brcf2br, 0., 0., 0., 0./)

      END IF                                                            
                                                                        
      icfcl3 = 0                                                        
      icf2cl2 = 0                                                       
      icf2clcfcl2 = 0                                                   
      icf2clcf2cl = 0                                                   
      icf2clcf3 = 0                                                     
      iccl4 = 0                                                         
      imeccl3 = 0                                                       
      ichf2cl = 0                                                       
      imecf2cl = 0                                                      
      imecfcl2 = 0                                                      
      icf2clbr = 0                                                      
      icf2br2 = 0                                                       
      icf3br = 0                                                        
      icf2brcf2br = 0                                                   
      imebr = 0                                                         
      imecl = 0                                                         
      in2o = 0                                                          
      ico2 = 0                                                          
      ich4 = 0                                                          
      ih2 = 0                                                           
      ich2br2 = 0 
      icos = 0
      DO i=1,n_lbc_specs                                                
        SELECT CASE (lbc_specs(i))                                      
          CASE('CFCl3     ')                                            
            icfcl3=i                                                    
          CASE('CF2Cl2    ')                                            
            icf2cl2=i                                                   
          CASE('CF2ClCFCl2')                                            
            icf2clcfcl2=i                                               
          CASE('CF2ClCF2Cl')                                            
            icf2clcf2cl=i                                               
          CASE('CF2ClCF3  ')                                            
            icf2clcf3=i                                                 
          CASE('CCl4      ')                                            
            iccl4=i                                                     
          CASE('MeCCl3    ')                                            
            imeccl3=i                                                   
          CASE('CHF2Cl    ')                                            
            ichf2cl=i                                                   
          CASE('MeCF2Cl   ')                                            
            imecf2cl=i                                                  
          CASE('MeCFCl2   ')                                            
            imecfcl2=i                                                  
          CASE('CF2ClBr   ')                                            
            icf2clbr=i                                                  
          CASE('CF2Br2    ')                                            
            icf2br2=i                                                   
          CASE('CF3Br     ')                                            
            icf3br=i                                                    
          CASE('CF2BrCF2Br')                                            
            icf2brcf2br=i                                               
          CASE('MeBr      ')                                            
            imebr=i                                                     
          CASE('MeCl      ')                                            
            imecl=i                                                     
          CASE('N2O       ')                                            
            in2o=i                                                      
          CASE('CH4       ')                                            
            ich4=i                                                      
          CASE('H2        ')                                            
            ih2=i                                                       
          CASE('CO2       ')                                            
            ico2=i                                                      
          CASE('CH2Br2    ')                                            
            ich2br2=i                                                   
          CASE('COS       ')                                            
            icos=i                                                   
        END SELECT                                                      
      END DO                                                            
                                                                        
                                                                        
! F11                                                                   
      full_lbc(1) = cfc11
! F12                                                       
      full_lbc(2) = cfc12
! F113                                                      
      full_lbc(3) = cfc113
! F114                                                      
      full_lbc(4) = cfc114
! F115                                                      
      full_lbc(5) = cfc115
! CCl4                                                      
      full_lbc(6) = ccl4
! MeCCl3                                                    
      full_lbc(7) = meccl3
! HCFC-22                                                   
      full_lbc(8) = hcfc22
! HCFC-141b                                                 
      full_lbc(9) = hcfc141b
! HCFC-142b                                                 
      full_lbc(10)= hcfc142b
! H-1211                                                    
      full_lbc(11)= h1211
! H-1202                                                    
      full_lbc(12)= h1202
! H-1301                                                    
      full_lbc(13)= h1301
! H-2402                                                    
      full_lbc(14)= h2402
! MeCl                                                      
      full_lbc(15)= mecl
! MeBr                                                      
      full_lbc(16)= mebr
! CH2Br2                                                                
      full_lbc(17)= ch2br2                                              
! COS
      full_lbc(18)= cos_mmr
! Add 6 pptv of bromine to MeBr to account for the effect of very       
! short-lived bromine gases - !! SHOULD WE STILL DO THIS PRE-INDUSTRIAL?                                             
      IF (ich2br2 == 0) THEN                                            
        full_lbc(16) = full_lbc(16) + 5.0e-12*c_mebr                    
      END IF                                                            
                                                                        
! Perform lumping                                                       
      IF (perform_lumping) THEN                                         
        DO i=1,n_full_spec                                              
          IF (replace_Cl(i) > 0) THEN                                   
            full_lbc(replace_cl(i)) = full_lbc(replace_cl(i)) +         &
                                      convfac_cl(i) * full_lbc(i)       
          END IF                                                        
          IF (replace_Br(i) > 0) THEN                                   
            full_lbc(replace_Br(i)) = full_lbc(replace_Br(i)) +         &
                                      convfac_Br(i) * full_lbc(i)       
          END IF                                                        
        END DO                                                          
      END IF                                                            
                                                                        
! Assign values to LBC array                                            
! Default values to 0                                                   
      lbc_mmr = 0.                                                      
      IF (icfcl3      > 0) lbc_mmr(icfcl3     ) = full_lbc( 1)          
      IF (icf2cl2     > 0) lbc_mmr(icf2cl2    ) = full_lbc( 2)          
      IF (icf2clcfcl2 > 0) lbc_mmr(icf2clcfcl2) = full_lbc( 3)          
      IF (icf2clcf2cl > 0) lbc_mmr(icf2clcf2cl) = full_lbc( 4)          
      IF (icf2clcf3   > 0) lbc_mmr(icf2clcf3  ) = full_lbc( 5)          
      IF (iccl4       > 0) lbc_mmr(iccl4      ) = full_lbc( 6)          
      IF (imeccl3     > 0) lbc_mmr(imeccl3    ) = full_lbc( 7)          
      IF (ichf2cl     > 0) lbc_mmr(ichf2cl    ) = full_lbc( 8)          
      IF (imecfcl2    > 0) lbc_mmr(imecfcl2   ) = full_lbc( 9)          
      IF (imecf2cl    > 0) lbc_mmr(imecf2cl   ) = full_lbc(10)          
      IF (icf2clbr    > 0) lbc_mmr(icf2clbr   ) = full_lbc(11)          
      IF (icf2br2     > 0) lbc_mmr(icf2br2    ) = full_lbc(12)          
      IF (icf3br      > 0) lbc_mmr(icf3br     ) = full_lbc(13)          
      IF (icf2brcf2br > 0) lbc_mmr(icf2brcf2br) = full_lbc(14)          
      IF (imecl       > 0) lbc_mmr(imecl      ) = full_lbc(15)          
      IF (imebr       > 0) lbc_mmr(imebr      ) = full_lbc(16)          
      IF (ich2br2     > 0) lbc_mmr(ich2br2    ) = full_lbc(17)          
      IF (icos        > 0) lbc_mmr(icos       ) = full_lbc(18)          
      IF (in2o        > 0) lbc_mmr(in2o       ) = n2o
      IF (ich4        > 0) lbc_mmr(ich4       ) = ch4
      IF (ih2         > 0) lbc_mmr(ih2        ) = h2
      IF (ico2        > 0) lbc_mmr(ico2       ) = co2
                                                                        
    IF (lhook) CALL dr_hook('UKCA_SCENARIO_PRESCRIBED',zhook_out,zhook_handle)
      RETURN
    END SUBROUTINE UKCA_SCENARIO_PRESCRIBED
