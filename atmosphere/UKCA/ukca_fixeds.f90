! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:  Given fixed supersaturation, S, calculates CCN(S) 
!               (number concentration of aerosol particles which would  
!                become "activated" into cloud droplets at that 
!                supersaturation.)
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of Cambridge,
!  University of Leeds, University of Oxford, and The Met Office.
!  See: www.ukca.ac.uk
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
!
!  based on *activ_box* activation box model, but now in 3D!
!
!  Author: Rosalind West, AOPP, Oxford  2008
!  -------
! ----------------------------------------------------------------------

      SUBROUTINE UKCA_FIXEDS(kbdim,                                     &
                             pccn,    pn,      pxtm1,                   &
                             ptm1,    prdry,                            & 
                             psfix,   nsfix              )
  

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
    ! The calculation of the activation at fixed supersaturation
    !    can be reduced to 3 tasks: 
    ! I)   Define the fixed (maximum) supersaturation
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
    ! Pruppacher and Klett, Kluewer Ac. Pub., 1997.


      USE UKCA_MODE_SETUP,  ONLY: nmodes,        &
                                  ncp,           &
                                  mm,            &
                                  rhocomp,       &
                                  no_ions,       &
                                  mode,          &
                                  component,     &
                                  sigmag

      USE UKCA_CONSTANTS,        ONLY: rmol,     &  ! gas constant
                                       mmw          ! H2O molecular weight kg/mol
      USE water_constants_mod,   ONLY: rho_water    ! density H2O kg/m^3

      USE yomhook,               ONLY: lhook, dr_hook
      USE parkind1,              ONLY: jprb, jpim
      IMPLICIT NONE

!--- Arguments:
      INTEGER, INTENT(IN) :: kbdim                ! No. of points
      INTEGER, INTENT(IN) :: nsfix                !

      REAL, INTENT(IN) :: ptm1(kbdim)             ! temperature
      REAL, INTENT(IN) :: prdry(kbdim,nmodes)     ! dry radius for each mode
      REAL, INTENT(IN) :: psfix(nsfix)            ! fixed maximum supersaturation
     
!--- Aerosol tracers mass mixing ratio
      REAL, INTENT(IN) :: pn (kbdim,nmodes)       ! aerosol number concentration
!                                                 ! for each mode [m-3]
      REAL, INTENT(IN) :: pxtm1(kbdim,nmodes,ncp) ! aerosol tracers mass mixing ratio

      REAL, INTENT(OUT):: pccn(kbdim,nsfix)       ! number of activated particles

!--- Local variables:
                
      INTEGER :: jmod                    !
      INTEGER :: jcp                     !
      INTEGER :: jl                      !
      INTEGER :: jsfix                   !

      REAL :: zsigmaln(nmodes)           ! ln(geometric std dev)
      REAL :: zmassfrac                  !
      REAL :: zeps                       !
      REAL :: zsum                       !
      REAL :: zerf_ratio                 !

      REAL :: zsumtop(kbdim)             !
      REAL :: zsumbot(kbdim)             !

      REAL :: zrc(kbdim,nmodes)          ! Critical radius of a dry aerosol
!                                        !  particle that becomes activated at
!                                        !  the ambient radius of activation
      REAL :: zfracn(kbdim,nmodes)       ! Fraction of activated aerosol
!                                        !  numbers for each mode
      REAL :: zmasssum(kbdim,nmodes)
      REAL :: za(kbdim,nmodes)           ! curvature parameter A of the Koehler equation
      REAL :: zb(kbdim,nmodes)           ! hygroscopicity parameter B of the Koehler eqn
            


      REAL, PARAMETER :: zsten = 75.0E-3 ! surface tension of H2O [J m-2] 
!                                        !   neglecting salts and temperature
!                                        !   (also tried P&K 5.12 - erroneous!)
 
!---RW: osmotic coefficient currently fixed. Is this available for each component?

      REAL, PARAMETER :: zosm = 1.0

!---constants needed which used to come from mo_constants:

!      REAL, PARAMETER :: rhoh2o= 1000.0   !=rho_water density of liquid water in kg/m3
!      REAL, PARAMETER :: argas = 8.314409 !=rmol universal gas constant in J/K/mol  
!      REAL, PARAMETER :: rv    = 461.51   ![unused] gas constant for water vapour 
!      REAL, PARAMETER :: api   = 3.14159265358979323846 ![not needed]
!      REAL, PARAMETER :: cpd   = 1005.46   ![not needed]
!      REAL, PARAMETER :: g     = 9.80665   ! [not needed]
!      REAL, PARAMETER :: alv   = 2.5008e6  ! [not needed]
!      REAL, PARAMETER :: amw   = 18.0154  !=mmw molecular weight of water vapour
!      REAL, PARAMETER :: amd   = 28.970   ! [not needed]
!      REAL, PARAMETER :: tmelt = 273.15   ! [not needed]

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_FIXEDS',zhook_in,zhook_handle)

!--- 0) Initializations:
    
      zmasssum(:,:) = 0.0
      za(:,:)       = 0.0
      zb(:,:)       = 0.0
      zrc(:,:)      = 1.0    ! [m] initialized with 1m, only changed if activation occurs
      pccn(:,:)     = 0.0
  
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

            zsumtop(:)    = 0.0
            zsumbot(:)    = 0.0
       
            DO jcp=1,ncp
               IF(component(jmod,jcp)) THEN
                  DO jl=1, kbdim
                     zmasssum(jl,jmod)=zmasssum(jl,jmod)+ &
                          pxtm1(jl,jmod,jcp)
                  END DO !jl
               END IF !cpt=true
            END DO !jcp

          
!--- Sum properties over the number of soluble compounds:
!    (nion=0 for insoluble compounds)
! N.B. Molar masses of components in UKCA already given in kg mol-1

            DO jcp=1,ncp
              IF(component(jmod,jcp) .AND.                              &
                 NINT(no_ions(jcp)) > 0) THEN
    
                DO jl=1, kbdim
                  IF(zmasssum(jl,jmod) > zeps) THEN
                         
                    zmassfrac = pxtm1(jl,jmod,jcp)/zmasssum(jl,jmod)
                        
                    zsumtop(jl) = zsumtop(jl) + pxtm1(jl,jmod,jcp)*     &
                                  no_ions(jcp)*zosm*                    &
                                  zmassfrac/mm(jcp)
                        
                    zsumbot(jl) = zsumbot(jl) + pxtm1(jl,jmod,jcp)/     &
                                  rhocomp(jcp) 
                        
                  END IF  ! zmasssum>0
                END DO     ! jl
              END IF           ! no_ions>0 and component
            END DO              ! loop over cpts

            DO jl=1, kbdim
              IF (zsumbot(jl)>zeps) THEN
               !--- 1.1.1) Hygroscopicity parameter B (Eq. 4):
                     
                zb(jl,jmod) = (mmw*zsumtop(jl))/(rho_water*zsumbot(jl))
                     
               !--- 1.1.2) Calculate the curvature parameter A:
                        
                za(jl,jmod) = (2.0*zsten*mmw)/(rho_water*rmol*ptm1(jl))
                    
              END IF        !zsumbot > 0
            END DO           !jl
          END IF              !mode 
         END DO                !jmod=1, nmodes

      !Rosalind test: try fixing aerosol properties (za and zb)
      
      !za(:,:,:)=0.1148E-08
      !zb(:,:,:)=0.7247
      
!--- 3) Calculate activation:

      DO jmod=1, nmodes
        IF(mode(jmod)) THEN
          DO jsfix=1, nsfix
            DO jl=1, kbdim
              IF (psfix(jsfix)  >zeps  .AND.                            &
                  pn(jl,jmod)   >zeps  .AND.                            &
                  za(jl,jmod)   >zeps  .AND.                            &
                  zb(jl,jmod)   >zeps ) THEN 
               !---3.1) Calculate the critical radius (12):
                zrc(jl,jmod) = (za(jl,jmod)/3.0)*                       &
                               ((4.0/zb(jl,jmod))**(1.0/3.0)/           &
                               (psfix(jsfix)**(2.0/3.0)))
               !--- 3.2) Calculate the total number of CCN larger than each
               !         of the mode critical radii for each value of fixed-S
                zerf_ratio = LOG(zrc(jl,jmod)/prdry(jl,jmod))/          &                
                             SQRT(2.0)/zsigmaln(jmod)
                pccn(jl,jsfix) = pccn(jl,jsfix) + 0.5*pn(jl,jmod)*      &
                                 (1-ERF(zerf_ratio))                              
              END IF
            END DO           ! jl
          END DO              ! jsfix
        END IF                 ! jmod
      END DO                    ! jmod

      IF (lhook) CALL dr_hook('UKCA_FIXEDS',zhook_out,zhook_handle)

      RETURN
      END SUBROUTINE ukca_fixeds
