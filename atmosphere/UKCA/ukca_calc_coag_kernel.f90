! *****************************COPYRIGHT*******************************
! 
! (c) [University of Leeds] [2008]. All rights reserved. 
! This routine has been licensed to the Met Office for use and 
! distribution under the UKCA collaboration agreement, subject  
! to the terms and conditions set out therein. 
! [Met Office Ref SC138]  
! 
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Calculate coagulation kernel for intra-model and inter-modal
!   coagulation between modes.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
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
      SUBROUTINE UKCA_CALC_COAG_KERNEL(NBOX,KII_ARR,KIJ_ARR,            &
       DRYDP,DVOL,WETDP,WVOL,RHOPAR,MFPA,DVISC,T,PMID,                  &
       COAG_ON,ICOAG)
!-------------------------------------------------------------------
!
! Purpose
! -------
! Calculate coagulation kernel for intra-model and inter-modal
! coagulation between modes  --> KII_ARR and KIJ_ARR.
!
! 4 possible methods of calculation chosen via switch ICOAG.
!
! 1) Method as Jacobson, "Fundamentals of Atmospheric Modelling", p. 446
!    for Brownian coagulation between ptls in the transition regime
!
!      K_ij=4*PI*(r_i+r_j)*(D_i+D_j)/(TERM1+TERM2)
!
!    where r_i,r_j are (wet) particle radii
!          D_i,D_j are particle diffusion coefficients
!
!          D_i=ACunn_i*kT/(6*PI*r_i*DVISC)
!
!          ACunn_i is the Cunningham slip correction
!                = 1.0+Kn_i*(1.257+0.4*EXP{-1.1/Kn_i})
!
!          Kn_i is the Knudsen number = MFPA/r_i
!
!          MFPA=2*DVISC/(RHOA*vel_i)
!
!          TERM1=(r_i+r_j)/(r_i+r_j+sqrt{delta_i^2+delta_j^2}
!          TERM2=4*(D_i+D_j)/( sqrt{vel_i^2+vel_j^2}*(r_i+r_j) )
!
!          delta_i,delta_j are mean distances from sphere centre after
!                          travelling distance of particle mfp (mfpp).
!          delta_i=( (2*r_i+mfpp_i)^3 - (4*r_i^2 + mfpp_i^2)^{3/2} )
!
!          mfpp_i= mean free path of particle = 2*D_i/PI/vel_i
!
!          vel_i,vel_j are thermal speed of particles
!                vel_i = sqrt{8kT/PI/mass_i}
!
!          where k is Boltzmann's constant
!                mass_i is the mass of a particle = rho_i*vol_i
!
! 2) Method as in HAM/M7, see Stier et al (2005), Vignati et al (2004)
!
!      K_ij=16*pi*RMID*DMID/(4*DMID/VMID/RMID + RMID/(RMID+DELTAMID))
!
!    where DMID,VMID,DELTAMID are D_i, vel_i and delta_i
!          evaluated for a particle with mean radius RMID=(r_i+r_j)/2
!
!          n.b. I have been informed since then that although this
!          is what is described in the Stier et al (2005) and
!          Vignati et al (2004), that in the latest HAM module,
!          DMID,VMID and DELTAMID are now evaluated more completely as:
!
!          DMID={D(r_i)+D(r_j)}/2.0
!          VMID=sqrt{vel_i^2+vel_j^2}
!          DELTAMID=sqrt{delta_i^2+delta_j^2}/sqrt{2.0}
!
! 3) Method as in UM sulfur scheme as in UMDP 20
!
!      K_ij=(2*k*T/3*DVISC)*(r_i+r_j)/{1/r_i + 1/r_j +
!                                      MFP*ACunn*((1/r_i)^2+(1/r_j)^2)}
!
!    where k=Boltzmann's constant
!          DVISC=dynamic viscosity of air
!          T=air temperature
!          MFP=MFPA=mean free path of air=6.6e-8*p0*T/(p*T0)
!                             [p0=1.01325e5 Pa,T0=293.15K]
!          ACunn=Cunningham slip correction assumed = 1.59 for all modes
!
! 4) Method as in 3) but with MFPP_i,MFPPj and ACunn_i,ACunn_j and MFPA
!                             calculated explicitly as in 1)
!
!      K_ij=(2*k*T/3*DVISC)*(r_i+r_j)/{1/r_i + 1/r_j +
!                         (MFPP_i*ACunn_i/r_i^2)+(MFPP_j*ACunn_j/r_j^2)}
!
! Inputs
! ------
! NBOX      : Number of grid boxes
! DRYDP     : Geometric mean dry diameter for mode (m)
! DVOL      : Dry volume corresponding to DRYDP (m^3)
! WETDP     : Wet diameter corresponding to DRYDP (m)
! WVOL      : Wet volume corresponding to WETDP (m^3)
! RHOPAR    : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
! MFPA      : Mean free path of air (m)
! DVISC     : Dynamic viscosity of air (kg m^-1 s^-1)
! T         : Air temperature (K)
! PMID      : Mid-point pressure (Pa)
! COAG_ON   : Switch for turning coagulation off
! ICOAG     : KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!
! Outputs
! -------
! KII_ARR   : Coag coeff for intra-modal coag (IMODE-IMODE) (cm^3/s)
! KIJ_ARR   : Coag coeff for inter-modal coag (IMODE-JMODE) (cm^3/s)
!
! Local variables
! ---------------
! KII       : Coag coeff for IMODE-IMODE intra-modal coag (cm^3/s)
! KIJ       : Coag coeff for IMODE-JMODE inter-modal coag (cm^3/s)
! RPI       : Geometric mean radius for mode IMODE (m)
! RPJ       : Geometric mean radius for mode JMODE (m)
! VPI       : Volume of particle of radius RPI (m^3)
! VPJ       : Volume of particle of radius RPJ (m^3)
! RHOI      : Density of particles in mode IMODE (kg/m3)
! RHOJ      : Density of particles in mode JMODE (kg/m3)
! MASK1     : Mask over which boxes to calculate kernels
! MASK2     : Mask over which boxes to calculate kernels
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! MODE      : Defines which modes are set
!
!-------------------------------------------------------------------

      USE UKCA_MODE_SETUP
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE
!
! Arguments
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: COAG_ON
      INTEGER, INTENT(IN) :: ICOAG
      REAL, INTENT(IN)    :: DRYDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: DVOL(NBOX,NMODES)
      REAL, INTENT(IN)    :: WETDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: WVOL(NBOX,NMODES)
      REAL, INTENT(IN)    :: RHOPAR(NBOX,NMODES)
      REAL, INTENT(IN)    :: MFPA(NBOX)
      REAL, INTENT(IN)    :: DVISC(NBOX)
      REAL, INTENT(IN)    :: T(NBOX)
      REAL, INTENT(IN)    :: PMID(NBOX)
      REAL, INTENT(OUT)   :: KII_ARR(NBOX,NMODES)
      REAL, INTENT(OUT)   :: KIJ_ARR(NBOX,NMODES,NMODES)
!
! Local variables
      INTEGER :: IMODE
      INTEGER :: JMODE
      REAL    :: KII(NBOX)
      REAL    :: KIJ(NBOX)
      REAL    :: VPI(NBOX)
      REAL    :: RPI(NBOX)
      REAL    :: RHOI(NBOX)
      REAL    :: VPJ(NBOX)
      REAL    :: RPJ(NBOX)
      REAL    :: RHOJ(NBOX)
      LOGICAL :: MASK1(NBOX)
      LOGICAL :: MASK2(NBOX)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
      IF (lhook) CALL dr_hook('UKCA_CALC_COAG_KERNEL',zhook_in,zhook_handle)
      MASK1(:) = .TRUE. ! set to calculate kernels for all boxes
      MASK2(:) = .TRUE. ! set to calculate kernels for all boxes
!
      DO IMODE=1,4 ! KII for intra-coag over soluble modes
       IF(MODE(IMODE)) THEN
!
        VPI(:)=WVOL(:,IMODE)
        RPI(:)=WETDP(:,IMODE)/2.0
        RHOI(:)=RHOPAR(:,IMODE)

! DEPENDS ON: ukca_coag_coff_v
        CALL UKCA_COAG_COFF_V(NBOX,MASK1,RPI,RPI,VPI,VPI,RHOI,RHOI,     &
             MFPA,DVISC,T,PMID,KII,COAG_ON,ICOAG)
!
        KII_ARR(:,IMODE)=KII(:) ! copy to KII_ARR
!
        DO JMODE=(IMODE+1),4 ! KIJ for coag with larger soluble modes

         IF(MODE(JMODE)) THEN
!
          VPJ(:)=WVOL(:,JMODE)
          RPJ(:)=WETDP(:,JMODE)/2.0
          RHOJ(:)=RHOPAR(:,JMODE)

! DEPENDS ON: ukca_coag_coff_v
          CALL UKCA_COAG_COFF_V(NBOX,MASK2,RPI,RPJ,VPI,VPJ,RHOI,RHOJ,   &
               MFPA,DVISC,T,PMID,KIJ,COAG_ON,ICOAG)
!
          KIJ_ARR(:,IMODE,JMODE)=KIJ(:) ! copy to KIJ_ARR
!
         END IF
        END DO ! loop over larger soluble modes
!
        DO JMODE=(IMODE+4),7 ! KIJ for coag with larger insoluble modes

         IF(MODE(JMODE)) THEN
!
          IF(MODESOL(JMODE) == 1) THEN
           RPJ(:)=WETDP(:,JMODE)/2.0
           VPJ(:)=WVOL(:,JMODE)
          ELSE
           RPJ(:)=DRYDP(:,JMODE)/2.0
           VPJ(:)=DVOL(:,JMODE)
          ENDIF
          RHOJ(:)=RHOPAR(:,JMODE)

! DEPENDS ON: ukca_coag_coff_v
          CALL UKCA_COAG_COFF_V(NBOX,MASK2,RPI,RPJ,VPI,VPJ,RHOI,RHOJ,   &
               MFPA,DVISC,T,PMID,KIJ,COAG_ON,ICOAG)
!
          KIJ_ARR(:,IMODE,JMODE)=KIJ(:) ! copy to KIJ_ARR
!
         END IF
        END DO ! loop over larger insoluble modes
!
       END IF
      END DO ! loop over soluble modes
!
      DO IMODE=5,NMODES ! loop over insoluble modes
       IF(MODE(IMODE)) THEN
!
        VPI(:)=WVOL(:,IMODE)
        RPI(:)=WETDP(:,IMODE)/2.0
        RHOI(:)=RHOPAR(:,IMODE)
! DEPENDS ON: ukca_coag_coff_v
        CALL UKCA_COAG_COFF_V(NBOX,MASK1,RPI,RPI,VPI,VPI,RHOI,RHOI,     &
             MFPA,DVISC,T,PMID,KII,COAG_ON,ICOAG)
!
        KII_ARR(:,IMODE)=KII(:) ! copy to KII_ARR
!
        DO JMODE=(IMODE-2),4 ! ins. mode coag with larger soluble modes
         IF(MODE(JMODE)) THEN

          VPJ(:)=WVOL(:,JMODE)
          RPJ(:)=WETDP(:,JMODE)/2.0
          RHOJ(:)=RHOPAR(:,JMODE)
! DEPENDS ON: ukca_coag_coff_v
          CALL UKCA_COAG_COFF_V(NBOX,MASK2,RPI,RPJ,VPI,VPJ,RHOI,RHOJ,   &
               MFPA,DVISC,T,PMID,KIJ,COAG_ON,ICOAG)
!
          KIJ_ARR(:,IMODE,JMODE)=KIJ(:) ! copy to KIJ_ARR
!
         END IF
        END DO
!
       END IF
      END DO

      IF (lhook) CALL dr_hook('UKCA_CALC_COAG_KERNEL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALC_COAG_KERNEL
