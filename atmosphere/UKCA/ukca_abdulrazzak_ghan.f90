! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description: Calculates number concentration of aerosol particles
!               which become "activated" into cloud droplets from MODE
!               aerosol results for mass mixing ratio, number and dry radius
!               together with updraught velocities.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford and The Met Office.
!  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
! ----------------------------------------------------------------------

      SUBROUTINE UKCA_ABDULRAZZAK_GHAN(kbdim,   klev,                   &
                                       pesw,    prho,                   &
                                       pn,      pxtm1,                  &
                                       ptm1,    papm1,                  &
                                       pqm1,    prdry,                  &
                                       nwbins,  pwarr,                  &
                                       pwpdf,   pwbin,                  &
                                       psmax,   pwchar,                 &
                                       pcdncactm, pcdncact,             &
                                       pcdncact3)


  

! *aero_activ* calculates the number of activated aerosol 
!              particles from the aerosol size-distribution,
!              composition and ambient supersaturation
!
! Author:
! -------
! Philip Stier, MPI-MET                  2002/2003
! Rosalind West, AOPP, Oxford            2008
!
! Method:
! -------
! The calculation of the activation can be reduced to 3 tasks:
! 
! I)   Calculate the maximum supersaturation
! II)  Calculate the corresponding radius of activation
!      for each mode
! III) Calculate the number of particles that are larger
!      then the radius of activation for each mode.
! 
! III) Calculation of the number of activated particles:
!      See the routine aero_activ_tail below.
!
! References:
! -----------
! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
! Pruppbacher and Klett, Kluewer Ac. Pub., 1997.


      USE UKCA_MODE_SETUP,  ONLY: nmodes,        &
                                  ncp,           &
                                  mm,            &
                                  rhocomp,       &
                                  no_ions,       &
                                  mode,          &
                                  component,     &
                                  sigmag

      USE UKCA_ACTIV_MOD,        ONLY: activclosest
      USE UKCA_CONSTANTS,        ONLY: rmol,     &  ! gas constant
                                       mmw,      &  ! H2O molecular weight kg/mol
                                       m_air        ! Dry air  ------ " ---------
      USE conversions_mod,       ONLY: zerodegc, &  ! in K
                                       pi
      USE earth_constants_mod,   ONLY: g
      USE atmos_constants_mod,   ONLY: cp           ! Specific heat at const p J/(kg.K)
      USE water_constants_mod,   ONLY: tm,        & ! t melt 273.15 K
                                       rho_water, & ! density H2O kg/m^3
                                       lc           ! latent heat of condensation J/kg
      USE PrintStatus_mod
      USE yomhook,               ONLY: lhook, dr_hook
      USE parkind1,              ONLY: jprb, jpim

      IMPLICIT NONE

!--- Arguments:

      INTEGER, INTENT(IN) :: kbdim                    !
      INTEGER, INTENT(IN) :: klev                     !
      INTEGER, INTENT(IN) :: nwbins                   !

      REAL, INTENT(IN)    :: ptm1(kbdim,klev)         ! temperature
      REAL, INTENT(IN)    :: papm1(kbdim,klev)        ! pressure 
      REAL, INTENT(IN)    :: prho(kbdim,klev)         ! air density
      REAL, INTENT(IN)    :: pqm1(kbdim,klev)         ! specific humidity
      REAL, INTENT(IN)    :: pesw(kbdim,klev)         ! saturation water vapour pressure
      REAL, INTENT(IN)    :: prdry(kbdim,klev,nmodes) ! dry radius for each mode    
      REAL, INTENT(IN)    :: pn (kbdim,klev,nmodes)   ! aerosol number concentration
!                                                     ! for each mode [m-3]
      REAL, INTENT(IN)    :: pxtm1(kbdim,klev,nmodes,ncp)! aerosol tracers mass
!                                                        ! mixing ratio
! pdf of updraught velocities:
      REAL,INTENT(IN)     :: pwarr(kbdim,klev,nwbins)  ! lin array of vert vel [m s-1]
      REAL,INTENT(IN)     :: pwpdf(kbdim,klev,nwbins)  ! lin array of pdf of w [m s-1]
      REAL,INTENT(IN)     :: pwbin(kbdim,klev,nwbins)  ! w bin width [m s-1]

      REAL, INTENT(OUT)   :: psmax(kbdim,klev)         ! max supersaturation [fraction]
      REAL, INTENT(OUT)   :: pwchar(kbdim,klev)        ! calculated characteristic
!                                                      ! updraught vel [m s-1] 
      REAL, INTENT(OUT)   :: pcdncactm(kbdim,klev,nmodes) ! expected number of activated
!                                                         !  particles in each mode
      REAL, INTENT(OUT)   :: pcdncact(kbdim,klev)      ! expected number of activated
!                                                      ! particles over all modes
      REAL, INTENT(OUT)   :: pcdncact3(kbdim,klev)     ! <N_d^-1/3> over all modes

!--- Local variables:
      REAL :: zsigmaln(nmodes)   ! ln(geometric std dev)
                
      INTEGER :: jmod            ! loop counters
      INTEGER :: jcp
      INTEGER :: jl
      INTEGER :: jk
      INTEGER :: jw

      REAL :: zmassfrac
      REAL :: zalpha
      REAL :: zeps
      REAL :: zgamma
      REAL :: zf
      REAL :: zg
      REAL :: zxi
      REAL :: zeta
      REAL :: zew
      REAL :: zsum
      REAL :: zgrowth
      REAL :: zdif
      REAL :: zxv
      REAL :: zk
      REAL :: zka
      REAL :: zerf_ratio
      REAL :: zwpwdw                      !
      REAL :: zpwdw                       !

      REAL :: zsmax(kbdim,klev)           ! maximum supersaturation
      REAL :: zsumtop(kbdim,klev)         ! 
      REAL :: zsumbot(kbdim,klev)         !
      REAL :: zw(kbdim,klev)              ! total vertical velocity[m s-1]

      REAL :: zrc(kbdim,klev,nmodes)      ! critical radius of a dry aerosol particle
!                                         ! that becomes activated at the ambient
!                                         ! radius of activation
      REAL :: zfracn(kbdim,klev,nmodes)   ! fraction of activated aerosol numbers
!                                         ! for each mode
      REAL :: zmasssum(kbdim,klev,nmodes) !
      REAL :: za(kbdim,klev,nmodes)       ! curvature parameter A of the Koehler eqn
      REAL :: zb(kbdim,klev,nmodes)       ! hygroscopicity parameter B of Koehler eqn
      REAL :: zsm(kbdim,klev,nmodes)      ! critical supersaturation for activating
!                                         ! particles with the mode No. median radius

      REAL :: zcdnc(kbdim,klev,nwbins)    ! CDNC calculated at each increment of w
      REAL :: zcdncm                      ! CDNC calculated at each increment of w
!                                         ! by mode
      REAL :: zndtopm(kbdim,klev,nmodes)  ! top line integral for calc pcdncactm
      REAL :: zndbotm(kbdim,klev,nmodes)  ! bottom line integral for calc pcdncactm
      REAL :: zndtop(kbdim,klev)          ! top line integral for calc pcdncact
      REAL :: zndbot(kbdim,klev)          ! bottom line integral for calc pcdncact
      REAL :: zndtop3(kbdim,klev)         ! top line integral for calc pcdncact3

! local variables used for finding wchar
      INTEGER :: zclsloc  
      REAL :: zcls 

      REAL, PARAMETER :: zsten = 75.0E-3 ! surface tension of H2O [J m-2] 
                                         !   neglecting salts and temperature
                                         !   (also tried P&K 5.12 - erroneous!)
!--- RW: local variables for diffusivity and conductivity calculation:
 
      REAL, PARAMETER :: p0=101325.0     ! in Pa

!--- RW: osmotic coefficient currently fixed. Is this available for each component?
      REAL, PARAMETER :: zosm  = 1.0
      REAL, PARAMETER :: rv    = rmol/mmw      ! gas constant for water vapour 

!---constants needed which now come from mo_constants:
!      REAL, PARAMETER :: t0=273.15             !in K   =zerodegc
!      REAL, PARAMETER :: rhoh2o= 1000.0        ! =rho_water density of liquid water in kg/m3
!      REAL, PARAMETER :: argas = 8.314409      ! =rmol universal gas constant in J/K/mol  
!      REAL, PARAMETER :: alv   = 2.5008e6     !=lc latent heat for vaporisation in J/kg
!      REAL, PARAMETER :: rv    = 461.51        ! gas constant for water vapour 
!      REAL, PARAMETER :: api   = 3.14159265358979323846 ! =pi
!      REAL, PARAMETER :: cpd   = 1005.46       !=cp specific heat of dry air at constant
!      REAL, PARAMETER :: g     = 9.80665       !=g gravity acceleration in m/s2
!      REAL, PARAMETER :: amw   = 18.0154       !=mmw molecular weight of water vapour
!      REAL, PARAMETER :: amd   = 28.970        !=m_air molecular weight of dry air
! The replacement molecular weights are already in kg/mol units
!      REAL, PARAMETER :: tmelt = 273.15        !=tm melting temperature of ice/snow
!---------------------------------------------------------------
      REAL :: cthomi 

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

  
      IF (lhook) CALL dr_hook('UKCA_ABDULRAZZAK_GHAN',zhook_in,zhook_handle)

      cthomi  = tm-35.0

!--- 0) Initializations:
    
      zmasssum(:,:,:) = 0.0
      zsmax(:,:)      = 0.0
      psmax(:,:)      = 0.0
      zsm(:,:,:)      = 0.0
      za(:,:,:)       = 0.0
      zb(:,:,:)       = 0.0
      zrc(:,:,:)      = 1.0  ! [m] initialized with 1m, only changed if activation occurs
      pcdncact(:,:)   = 0.0
      pcdncact3(:,:)  = 0.0
      pcdncactm(:,:,:)= 0.0
      zcdnc(:,:,:)    = 0.0
      zcdncm          = 0.0
      zndtopm(:,:,:)  = 0.0
      zndbotm(:,:,:)  = 0.0
      zndtop(:,:)     = 0.0
      zndbot(:,:)     = 0.0
      zndtop3(:,:)    = 0.0
  
      zeps=EPSILON(1.0)

!--- calculate ln(sigmag)

      DO jmod=1, nmodes
         zsigmaln(jmod) = LOG(sigmag(jmod))
      END DO

!--- 1) Calculate properties for each aerosol mode:

     !--- 1.1) Calculate the auxiliary parameters A and B of the Koehler equation:

      DO jmod=1, nmodes
         IF(mode(jmod)) THEN
 
      !--- 1.1.0 Initializations:

            zsumtop(:,:)    = 0.0
            zsumbot(:,:)    = 0.0

       !--- 1.1.1) Calculate the mean hygroscopicity parameter B:

       !--- 1) Calculate weighted properties:
       !       (Abdul-Razzak & Ghan, 2000)

            DO jcp=1,ncp
               IF(component(jmod,jcp)) THEN
                  DO jk=1, klev
                     DO jl=1, kbdim
                        zmasssum(jl,jk,jmod) = zmasssum(jl,jk,jmod) +   &
                                                pxtm1(jl,jk,jmod,jcp)
                     END DO
                  END DO
               END IF
            END DO

          
!--- Sum properties over the number of soluble compounds:
!    (nion=0 for insoluble compounds)
! N.B. Molar masses of components in UKCA already given in kg mol-1

            DO jcp=1,ncp
              IF(component(jmod,jcp) .AND.                              &
                 NINT(no_ions(jcp)) > 0) THEN
    
                DO jk=1, klev
                  DO jl=1, kbdim
                    IF (zmasssum(jl,jk,jmod) > zeps) THEN
                         
                        zmassfrac=pxtm1(jl,jk,jmod,jcp)/                &
                                   zmasssum(jl,jk,jmod)
                         
                        zsumtop(jl,jk)=zsumtop(jl,jk)+                  &
                                       pxtm1(jl,jk,jmod,jcp)*           &
                                       no_ions(jcp)*zosm*               &
                                       zmassfrac/mm(jcp)
                        
                        zsumbot(jl,jk)=zsumbot(jl,jk)+                  &
                                       pxtm1(jl,jk,jmod,jcp)/           &
                                       rhocomp(jcp) 
                        
                    END IF     ! zmasssum>0
                  END DO        ! jl
                END DO           ! jk
              END IF              ! no_ions>0 and component
            END DO                 ! loop over cpts

  
         DO jk=1, klev
            DO jl=1, kbdim
               IF (zsumbot(jl,jk)>zeps) THEN
                  !--- 1.1.1) Hygroscopicity parameter B (Eq. 4):
                    
                  zb(jl,jk,jmod)=(mmw*zsumtop(jl,jk))/ &
                                 (rho_water*zsumbot(jl,jk))
              
                  !--- 1.1.2) Calculate the curvature parameter A:

                  za(jl,jk,jmod)=(2.0*zsten*mmw)/ &
                                 (rho_water*rmol*ptm1(jl,jk))

               END IF         !zsumbot > 0
            END DO            !jl
          END DO               !jk
        END IF                  !mode 
      END DO                    !jmod=1, nmodes

!--- 2) Calculate maximum supersaturation at each increment of vertical velocity pdf
    
    DO jw=1, nwbins   

!    --- 2.2) Abbdul-Razzak and Ghan (2000):
!         (Equations numbers from this paper unless otherwise quoted)

      DO jk=1, klev
         DO jl=1, kbdim
            ! get updraught velocity from pdf array
            zw(jl,jk)=pwarr(jl,jk,jw)
             
             !--- Water vapour pressure:

             zew=pqm1(jl,jk)*prho(jl,jk)*rv*ptm1(jl,jk)

             IF(zw(jl,jk)  >zeps   .AND.                                &
                pqm1(jl,jk)>zeps   .AND.                                &
                ptm1(jl,jk)>cthomi       ) THEN  
               
                !--- Abdul-Razzak et al. (1998) (Eq. 11):

                zalpha=(g*mmw*lc)/(cp*rmol*ptm1(jl,jk)**2) -            &
                       (g*m_air)/(rmol*ptm1(jl,jk))

                zgamma=(rmol*ptm1(jl,jk))/(pesw(jl,jk)*mmw) +           &
                       (mmw*lc**2)/(cp* papm1(jl,jk)*m_air*ptm1(jl,jk))

                !--- Diffusivity of water vapour in air (P&K, 13.3) [m2 s-1]:
                !---  RW Code from Steve Ghan
          
                zdif=0.211*(p0/papm1(jl,jk))*                           &
                     ((ptm1(jl,jk)/zerodegc)**1.94)*1.E-4 !in m^2/s
                      
                !--- Thermal conductivity zk (P&K, 13.18) [cal cm-1 s-1 K-1]:

                ! Mole fraction of water:

                zxv=pqm1(jl,jk)*(m_air/mmw)

                zka=(5.69+0.017*(ptm1(jl,jk)-273.15))*1.E-5

                ! Moist air, convert to [J m-1 s-1 K-1]:
                
                zk=zka*4.186*1.E2
               
! --- Abdul-Razzak et al. (1998) (Eq. 16):

                zgrowth=1.0/                                            &
                          ( (rho_water*rmol*ptm1(jl,jk))/               &
                            (pesw(jl,jk)*zdif*mmw) +                    &
                            (lc*rho_water)/(zk*ptm1(jl,jk)) *          & 
                          ((lc*mmw)/(ptm1(jl,jk)*rmol) -1.0) )
           
                !--- Summation for equation (6):

                zsum=0.0

                DO jmod=1, nmodes
                   IF(mode(jmod)) THEN
                      IF (pn(jl,jk,jmod)   >zeps     .AND.              &
                           prdry(jl,jk,jmod)>1.E-9   .AND.              &
                           za(jl,jk,jmod)   >zeps    .AND.              &
                           zb(jl,jk,jmod)   >zeps        ) THEN

                         ! (7):
                         zf=0.5*EXP(2.5*zsigmaln(jmod)**2) 
                         
                         ! (8):
                         zg=1.0+0.25*zsigmaln(jmod)
                         
                         ! (10):
                         zxi=2.0*za(jl,jk,jmod)/3.0 *                   &
                              SQRT(zalpha*zw(jl,jk)/zgrowth)
                         
                         ! (11):
                         zeta=((zalpha*zw(jl,jk)/zgrowth)**1.5) /       &
                              (2.0*pi*rho_water*zgamma*pn(jl,jk,jmod))
                         
                         ! (9):
                         zsm(jl,jk,jmod)=2.0/SQRT(zb(jl,jk,jmod)) *     &
                              (za(jl,jk,jmod)/                          &
                              (3.0*prdry(jl,jk,jmod)))**1.5
                         
                         ! (6): 
                         zsum=zsum + (1.0/zsm(jl,jk,jmod)**2 *          &
                              ( zf*(zxi/zeta)**1.5 +                    &
                              zg*(zsm(jl,jk,jmod)**2/                 &
                              (zeta+3.0*zxi))**0.75) )
                      END IF
                   END IF !mode
                END DO ! jmod

                IF (zsum > zeps) THEN
                   zsmax(jl,jk)=1.0/SQRT(zsum)
                ELSE
                   zsmax(jl,jk)=0.0
                END IF

             ELSE
                zsmax(jl,jk)=0.0
             END IF
             psmax(jl,jk)=zsmax(jl,jk)
          END DO ! jl
       END DO ! jk
   
!--- 3) Calculate activation:

!   ---3.1) Calculate the critical radius (12):
       
       DO jk=1, klev
         DO jl=1, kbdim
          ! Calculate p(w)*dw
           zpwdw=pwpdf(jl,jk,jw)*pwbin(jl,jk,jw)
           ! Calculate w*p(w)*dw
           zwpwdw= pwarr(jl,jk,jw)*zpwdw
           IF (zwpwdw < zeps .AND. Printstatus == prstatus_diag) THEN
             WRITE(6,*) 'RW: zwpwdw(',jl,',',jk,')= ',  zwpwdw
             WRITE(6,*) 'RW: pwpdf(',jl,',',jk,',',jw,')=  ', pwpdf(jl,jk,jw)
             WRITE(6,*) 'RW: pwbin(',jl,',',jk,',',jw,')=  ', pwbin(jl,jk,jw)
             WRITE(6,*) 'RW: pwarr(',jl,',',jk,',',jw,')=  ', pwarr(jl,jk,jw)
          END IF
          DO jmod=1, nmodes
            IF(mode(jmod)) THEN
               IF (psmax(jl,jk)           >zeps  .AND.                  &
                    zsm(jl,jk,jmod)       >zeps  .AND.                  &
                    pn(jl,jk,jmod)        >zeps  .AND.                  &
                    prdry(jl,jk,jmod)     >1.E-9  ) THEN 

                 zrc(jl,jk,jmod)=prdry(jl,jk,jmod)*(zsm(jl,jk,jmod)/    &
                                  psmax(jl,jk))**(2.0/3.0)
                      !--- 3.2) Calculate the total number of activated droplets 
                      !         larger than the critical radii for each mode
                 zerf_ratio=LOG(zrc(jl,jk,jmod)/prdry(jl,jk,jmod))/     &
                             SQRT(2.0)/zsigmaln(jmod)
                 zcdnc(jl,jk,jw)=zcdnc(jl,jk,jw) + 0.5*pn(jl,jk,jmod)*  &
                                  (1-ERF(zerf_ratio))  
                      ! separate out by mode  
                 zcdncm = 0.5*pn(jl,jk,jmod)*(1-ERF(zerf_ratio))
                      !! Calculate the expected value of CDNC 
                      !! for each mode over the whole pdf of w 
                      !!pcdncactm(jl,jk,jmod)=pcdncactm(jl,jk,jmod)+    &
                      !!     zcdncm*pwpdf(jl,jk,jw)*pwbin(jl,jk,jw)
                      ! Calculate top integral for 
                      ! expected value of CDNC for each mode 
                      ! over the whole pdf of w, weighted by w:
                      !zndtopm(jl,jk,jmod)= zndtopm(jl,jk,jmod)+ (zcdncm*zwpwdw)
                 zndtopm(jl,jk,jmod) = zndtopm(jl,jk,jmod)+             &
                                       (zcdncm*zpwdw)
               END IF
                   ! Calculate bottom integral for 
                   ! expected value of CDNC for each mode 
                   ! over the whole pdf of w, weighted by w:
                   !zndbotm(jl,jk,jmod)= zndbotm(jl,jk,jmod)+ zwpwdw
                   zndbotm(jl,jk,jmod)= zndbotm(jl,jk,jmod)+ zpwdw
            END IF !mode
            END DO ! jmod

             !! Calculate the expected value of CDNC
             !! over the pdf of w                
             !!pcdncact(jl,jk)=pcdncact(jl,jk) + &
             !!     zcdnc(jl,jk,jw)*pwpdf(jl,jk,jw)*pwbin(jl,jk,jw)
             ! Calculate top and bottom integrals for 
             ! expected value of CDNC 
             ! over the pdf of w, weighted by w: 
             !zndtop(jl,jk)= zndtop(jl,jk)+ (zcdnc(jl,jk,jw)*zwpwdw)
            zndtop(jl,jk)= zndtop(jl,jk)+ (zcdnc(jl,jk,jw)*zpwdw)
             !zndbot(jl,jk)= zndbot(jl,jk)+ zwpwdw
            zndbot(jl,jk)= zndbot(jl,jk)+ zpwdw
             !! Calculate the expected value of CDNC^-1/3
             !! over the pdf of w
            IF (zcdnc(jl,jk,jw) > zeps) THEN
              zndtop3(jl,jk)= zndtop3(jl,jk)+                            &
                              ((zcdnc(jl,jk,jw)**(-1.00/3.00))*zpwdw)
            ELSE
              zndtop3(jl,jk)= zndtop3(jl,jk)
            END IF
          END DO ! jl
        END DO ! jk
      END DO !jw 

    ! Calculate the normalised expected value of CDNC for each mode
    ! weighted by w
      DO jk=1, klev
        DO jl=1, kbdim
          DO jmod=1, nmodes
             IF(mode(jmod)) THEN
               IF  (zndbotm(jl,jk,jmod) < zeps .AND.                    &
                     PrintStatus == Prstatus_diag) THEN
                   WRITE(6,*) 'RW: zndbotm(',jl,',',jk,',',jmod,')=  ', &
                               zndbotm(jl,jk,jmod)
                ELSE
                   pcdncactm(jl,jk,jmod)= zndtopm(jl,jk,jmod)/          &
                                           zndbotm(jl,jk,jmod)
                END IF
             END IF !mode
          END DO ! jmod
          ! Calculate the normalised expected value of total CDNC and CDNC^-1/3
          ! weighted by w
          IF  (zndbot(jl,jk) > zeps) THEN
             pcdncact(jl,jk)= zndtop(jl,jk)/zndbot(jl,jk)
             pcdncact3(jl,jk)=zndtop3(jl,jk)/zndbot(jl,jk)
          ELSE
            IF (PrintStatus == Prstatus_diag) THEN
             WRITE(6,*) 'RW: zndbot(',jl,',',jk,')=  ',                 &
                        zndbot(jl,jk)
            END IF
          END IF
        END DO ! jl
      END DO ! jk

! Calculate the normalised expected value of total CDNC 
! weighted by w
! pcdncact(:,:)= zndtop(:,:)/zndbot(:,:)

! Now calculate the characteristic updraught velocity:
! i.e. find w* for which E(N_d(w))=N_d(w*)

     DO jk=1, klev
        DO jl=1, kbdim
           CALL ACTIVCLOSEST(zcdnc(jl,jk, :), nwbins, pcdncact(jl,jk),  &
                               zcls, zclsloc)
           pwchar(jl,jk) = pwarr(jl,jk,zclsloc)
        END DO ! jl
     END DO ! jk

     IF (lhook) CALL dr_hook('UKCA_ABDULRAZZAK_GHAN',zhook_out,zhook_handle)

   RETURN
   END SUBROUTINE UKCA_ABDULRAZZAK_GHAN
