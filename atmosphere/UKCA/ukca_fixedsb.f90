! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description: Given fixed supersaturation, S, calculates CCN(S) 
!              (number concentration of aerosol particles which would  
!               become "activated" into cloud droplets at that 
!               supersaturation.)
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
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
!  
!
! ----------------------------------------------------------------------

      SUBROUTINE ukca_fixedsb(kbdim,   klev,                            &
                             pccn,    pn,      pxtm1,                   &
                             ptm1,    prdry,                            &
                             psfix,   nsfix  )
  

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
    ! Pruppacher and Klett, Kluewer Ac. Pub., 1997.


      USE UKCA_MODE_SETUP,  ONLY: nmodes,        &
                                  ncp,           &
                                  mm,            &
                                  rhocomp,       &
                                  no_ions,       &
                                  mode,          &
                                  component,     &
                                  sigmag,        &
                                  modesol

      USE UKCA_CONSTANTS,        ONLY: rmol,     &  ! gas constant
                                       mmw          ! H2O molecular weight kg/mol
      USE water_constants_mod,   ONLY: rho_water    ! density H2O kg/m^3

      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

    !--- Arguments:

      INTEGER, INTENT(IN) :: kbdim                   ! No of points
      INTEGER, INTENT(IN) :: klev                    ! No of levels
      INTEGER, INTENT(IN) :: nsfix                   ! 


      REAL, INTENT(IN) :: ptm1(kbdim,klev)           ! temperature
      REAL, INTENT(IN) :: prdry(kbdim,klev,nmodes)   ! dry radius for each mode
      REAL, INTENT(IN) :: psfix(nsfix)               ! fixed maximum supersaturation
     

      REAL, INTENT(IN) :: pn (kbdim,klev,nmodes)     ! aerosol number concentration
!                                                    ! for each mode [m-3]
      REAL, INTENT(IN) :: pxtm1(kbdim,klev,nmodes,ncp)! aerosol tracers mass mix ratio
      REAL, INTENT(OUT):: pccn(kbdim,klev,nsfix)     ! number of activated particles

    !--- Local variables:

      REAL :: zsigmaln(nmodes)   ! ln(geometric std dev)
                
      INTEGER :: jmod,      jcp,         &
                 jl,        jk,          &
                 jsfix        

      REAL ::  zeps
      REAL ::  zerf_ratio

      REAL :: zrc(kbdim,klev,nmodes)       ! critical radius of a dry aerosol particle
!                                          ! that becomes activated at the ambient
!                                          ! radius of activation
      REAL :: zfracn(kbdim,klev,nmodes)    ! fraction of activated aerosol numbers
!                                          ! for each mode
      REAL :: zmasssum(kbdim,klev,nmodes)  !  
      REAL :: za(kbdim,klev,nmodes)        ! curvature parameter A of the Koehler eqn
      REAL :: zb(kbdim,klev,nmodes)        ! hygroscopicity parameter B of Koehler eqn
            


      REAL, PARAMETER :: zsten = 75.0E-3 ! surface tension of H2O [J m-2] 
                                              !   neglecting salts and temperature
                                              !   (also tried P&K 5.12 - erroneous!)
!--- RW: local variables for diffusivity and conductivity calculation:
 
      REAL, PARAMETER :: p0=101325.0        ! in Pa
      REAL, PARAMETER :: t0=273.15          ! in K

!---RW: osmotic coefficient currently fixed. Is this available for each component?

      REAL, PARAMETER :: zosm=1.0

!---constants needed which used to come from mo_constants:

!      REAL, PARAMETER :: rhoh2o= 1000.0        ! density of liquid water in kg/m3
!      REAL, PARAMETER :: argas = 8.314409      ! universal gas constant in J/K/mol  
!      REAL, PARAMETER :: amw   = 18.0154       ! molecular weight of water vapour

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_FIXEDSB',zhook_in,zhook_handle)

!--- 0) Initializations:
    
      zmasssum(:,:,:) = 0.0
      za(:,:,:)       = 0.0
      zb(:,:,:)       = 0.0
      zrc(:,:,:)      = 1.0     ! [m] initialized with 1m, only changed if activation occurs
      pccn(:,:,:)   = 0.0
  
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

            IF(modesol(jmod)==1) THEN
               DO jk=1, klev
                  DO jl=1, kbdim
                        !--- 1.1.1) Hygroscopicity parameter B (Eq. 4):
                        
                     ! Set to entirely soluble sulphate
                     zb(jl,jk,jmod) = (mmw*no_ions(1)*rhocomp(1))/      &
                                      (rho_water*mm(1))
                        !--- 1.1.2) Calculate the curvature parameter A:
                        
                     za(jl,jk,jmod) = (2.0*zsten*mmw)/                  &
                                      (rho_water*rmol*ptm1(jl,jk))
                    
                  END DO           !jl
               END DO              !jk
            END IF               !modesol
         END IF                 !mode 
      END DO                  !jmod=1, nmodes

      !Rosalind test: try fixing aerosol properties (za and zb)
      
      !za(:,:,:)=0.1148E-08
      !zb(:,:,:)=0.7247
      
    !--- 3) Calculate activation:

      DO jmod=1, nmodes
        IF(mode(jmod)) THEN
          DO jsfix=1, nsfix
            DO jk=1, klev
              DO jl=1, kbdim
                IF (psfix(jsfix)      >zeps  .AND. &
                   pn(jl,jk,jmod)   >zeps  .AND. &           
                   za(jl,jk,jmod)   >zeps  .AND. &
                   zb(jl,jk,jmod)   >zeps ) THEN 
                        !---3.1) Calculate the critical radius (12):
                   zrc(jl,jk,jmod) = (za(jl,jk,jmod)/3.0)*              &
                                     ((4.0/zb(jl,jk,jmod))**(1.0/3.0)/  &
                                     (psfix(jsfix)**(2.0/3.0)))
                        !--- 3.2) Calculate the total number of CCN larger than
                        !         each of the mode  critical radii for each value
                        !         of fixed-S
                   zerf_ratio = LOG(zrc(jl,jk,jmod)/prdry(jl,jk,jmod))/ &
                                SQRT(2.0)/zsigmaln(jmod)
                   pccn(jl,jk,jsfix) = pccn(jl,jk,jsfix)+               &
                                       0.5*pn(jl,jk,jmod)*              &
                                       (1-ERF(zerf_ratio))
                END IF
              END DO           ! jl
            END DO             ! jk
          END DO               ! jsfix
        END IF                ! mode(jmod)
      END DO                 ! jmod

      IF (lhook) CALL dr_hook('UKCA_FIXEDSB',zhook_out,zhook_handle)

      RETURN
      END SUBROUTINE UKCA_FIXEDSB
