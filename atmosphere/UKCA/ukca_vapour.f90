! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
! *********************************************************************
! To calculate the vapour pressure of H2SO4 before the 
! condensation routine. From Ayers (1980).
!
! Inputs
! NBOX Number of gridboxes
! T    Temperature (K)
! S    Specific humidity (kg/kg)
! SM   grid box mass of air (kg)
! PMID Air pressure at mid point (Pa)
! RP   ?
!
!
! Calculated by subroutine
! TLOG  Natural log of T
! UST   1/T
! BH2O  Partial pressure of water vapour (atm)
! PWLOG Natural log of BH2O
! PATM  Pressure of air (atm)
! WS    Composition of aerosol (wt fraction H2SO4) 
! A,B,XSB,MSB,MUH2SO4 all avriables used to calculate PH2SO4
! PH2SO4 Vapour pressure of H2SO4 (atm)
!
! Outputs
! VPKEL Vapour pressure of H2SO4 (Pa)
! WTS   Catch for WS being less than 40%
! RHOSOL_STRAT Density of particle solution [excl. insoluble cpts] (kg/m3)
! **********************************************************************

      SUBROUTINE UKCA_VAPOUR(NBOX,T,PMID,S,RP,VPKEL,WTS,RHOSOL_STRAT)
               

      USE UKCA_CONSTANTS,      ONLY: mmsul, rr, rhosul
      USE parkind1,            ONLY: jprb, jpim
      USE yomhook,             ONLY: lhook, dr_hook

      IMPLICIT NONE

! Arguments
      INTEGER, INTENT(IN) :: NBOX                   ! Number of gridboxes

      REAL, INTENT(IN)    :: T(NBOX)                ! Temperature (K)
      REAL, INTENT(IN)    :: PMID(NBOX)             ! Air pressure at mid point (Pa)
      REAL, INTENT(IN)    :: S(NBOX)                ! Specific humidity (kg/kg)
      REAL, INTENT(IN)    :: RP(NBOX)               ! ?
      REAL, INTENT(OUT)   :: WTS(NBOX)              ! MAX(41.0, WS*100)
      REAL, INTENT(OUT)   :: VPKEL(NBOX)            ! Vapour pressure of H2SO4 (Pa)
      REAL, INTENT(OUT)   :: RHOSOL_STRAT(NBOX)     ! Density of particle solution (kg/m3)

! Local variables
      INTEGER :: JL

      REAL :: VPH2SO4(NBOX)
      REAL :: BH2O(NBOX), PATM
      REAL :: WS(NBOX)
      REAL :: A,B,XSB(NBOX),MSB(NBOX),PWLOG(NBOX),TLOG(NBOX),UST(NBOX)
      REAL :: MUH2SO4(NBOX)
      REAL :: PH2SO4(NBOX)
      REAL :: PP_ASAD_H2SO4(NBOX)

! select density value from LUT 

      INTEGER :: k
      INTEGER :: ROUND(NBOX)
      REAL    :: T_DIFF(NBOX)
!  Data for density @253K, per Kelvin difference at weight % of H2SO4
      INTEGER, PARAMETER :: ndata = 12 ! no of data points

! percentage of H2SO4 present in aerosol
      INTEGER, PARAMETER :: percent(ndata) =                            &
               (/40,45,50,55,60,65,70,75,80,85,90,95/)

! per Kelvin difference in density
      REAL, PARAMETER :: k_diff(ndata) =                                &
       (/0.80,0.82,0.84,0.86,0.87,0.89,0.92,0.95,0.98,1.02,1.06,1.10/)

! density @ 253K kgm-3
      REAL, PARAMETER :: data253(ndata) =                               &
       (/1333.8,1381.1,1431.0,1483.3,1537.3,1592.4,1647.6,              &
         1701.6,1753.0,1800.1,1840.9,1873.2/)
  
      REAL :: C,D
      REAL, PARAMETER :: XSB_EPS=1.0e-6
      REAL, PARAMETER :: p0=101325.0 

      REAL, PARAMETER :: KS1=-21.661
      REAL, PARAMETER :: KS2= 2724.2
      REAL, PARAMETER :: KS3= 51.81
      REAL, PARAMETER :: KS4=-15732.0    
      REAL, PARAMETER :: KS5=47.004
      REAL, PARAMETER :: KS6=-6969.0
      REAL, PARAMETER :: KS7=-4.6183
      REAL, PARAMETER :: SURFTEN=0.0728
      REAL, PARAMETER :: BMINATM=2.0E-8
      REAL, PARAMETER :: BMAXATM=2.0E-6

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
  
      IF (lhook) CALL dr_hook('UKCA_VAPOUR',zhook_in,zhook_handle)

! Loop around all the boxes
      DO JL=1,NBOX

        PATM=PMID(JL)/p0
 
! calculate vmr and then atm of water vapour

        BH2O(JL)=1.609*S(JL)*PATM
        IF(BH2O(JL) < BMINATM) BH2O(JL)=BMINATM
        IF(BH2O(JL) > BMAXATM) BH2O(JL)=BMAXATM
                             
      
! Calculate H2SO4 vapour pressure from expression for the H2SO4 vp              
!  over H2SO4/H2O solutions from paper of Ayers + al., GRL, 7, 433-436, 1980      
!  The Giaugue expression for chemical potential has been fitted               
!  to a simple expression (see the paper). It is likely to be                  
!  a bit inaccurate for use at low strat conditions, but should be reasonable  
  
!  log(T)and water vapour                                      
        TLOG(JL)=LOG(T(JL)) 
        UST(JL)=1/T(JL) 
        PWLOG(JL)=LOG(BH2O(JL))                               
                                                                          
!  The H2SO4/H2O pure solution concentration                           
                                                                               
!  Mole fraction of H2SO4 in the binary solution                           
        A=KS1 + KS2*UST(JL)                                                   
        B=KS3 + KS4*UST(JL) 

        C=KS5 + KS6*UST(JL) + KS7*TLOG(JL) - PWLOG(JL)
        D=A*A - 4.0*B*C
        IF(D < 0.0) D=0.0
                                                      
!!   XSB(JL)=(-A -SQRT(A*A - 4.0*B*(KS5 + KS6*UST(JL) + KS7*TLOG(JL) - PWLOG(JL))))/(2.0*B)

        XSB(JL)=(-A -SQRT(D))/(2.0*B)
        IF(XSB(JL) < XSB_EPS) XSB(JL)=XSB_EPS

        MSB(JL)=55.51*XSB(JL)/(1.0-XSB(JL))                    
        WS(JL) = MSB(JL)*0.098076/(1.0 + MSB(JL)*0.098076)
                                              
!    wt % of H2SO4                                                              
        WTS(JL)=MAX(41.0, WS(JL)*100.0) 
                                              
        MUH2SO4(JL)=4.184*(1.514E4-286.0*(WTS(JL)-40.0)+1.080*          &
                          (WTS(JL)-40.0)**2-3941.0/(WTS(JL)-40.0)**0.1) 
                                     
!    Equilibrium H2SO4 vapour pressure (atm) using Ayers et al 1980.GRL.7.433-436
!          PH2SO4(JL)=EXP(-10156.0*UST(JL) + 16.2590 - MUH2SO4(JL)/(8.314*T(JL)))

!   Addendum to Ayers et al (1980) by Kulmala and Laaksonen.1990.J.Chem.Phys.93.696-701
!   K&L say Ayers cannot be used outside T range 338-445K. Correction using Ayers
!   at 360K (ln vp=-10156/360+16.259)=-11.9521
!
        PH2SO4(JL)= -11.9521+10156.0*(-(1.0/T(JL))+(1.0/360.0)+         &
               0.38/(905.0-360.0)*(1.0+LOG(360.0/T(JL))-(360.0/T(JL))))
        PH2SO4(JL)= EXP(PH2SO4(JL) - MUH2SO4(JL)/(8.314*T(JL)))


!    and in Pa  
        VPH2SO4(JL)=PH2SO4(JL)*P0

        IF((T(JL) > 240.0).AND.(PMID(JL) > 1.0E4)) VPH2SO4(JL)=0.0 ! set to zero if T>240K
! .. sets vapour pressure to zero if T>240K and pressure greater than 100 hPa
! .. (additional pressure criterion avoids setting it to zero in mid-upper stratosphere. 

! Add the Kelvin effect (still in Pa)
        VPKEL(JL)=VPH2SO4(JL)*EXP((2.0*SURFTEN*MMSUL)/(RHOSUL*RR*T(JL)  &
                    *RP(JL)))

!  New bit of code to calculate the density of an H2SO4-H2O mixture, using Martin
!  et al. 2000. GRL.27. 197-200. Uses weight % of H2SO4 calculated above.
!  1300 kg m-3 is about the minimum RHOSOL_STRAT can be.
        RHOSOL_STRAT(JL)=1300.0
!  need to take WTS and work out nearest integer of 5 to correspond to data
        ROUND(JL)=(NINT(WTS(JL)/5))*5

! find temperature difference - T is allowed to be negative here.
        T_DIFF(JL)=253.0-T(JL)

! read data, find value of n that corresponds to % H2SO4 data number
        DO k=1,ndata
          IF (ROUND(JL)==percent(k))THEN
            RHOSOL_STRAT(JL)=data253(k)+k_diff(k)*T_DIFF(JL)
          END IF ! round(jl)=percent(k)
        END DO ! loop around ndata

      END DO  ! end loop around boxes
                   
                   
      IF (lhook) CALL dr_hook('UKCA_VAPOUR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_VAPOUR
