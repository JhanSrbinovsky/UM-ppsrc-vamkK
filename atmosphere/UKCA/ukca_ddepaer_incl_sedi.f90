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
!    Calculates aerosol dry deposition and sedimentation.
!    Based on the parameterisation of Zhang et al (2001) which
!    uses the method in the model of Slinn (1982).
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
      SUBROUTINE UKCA_DDEPAER_INCL_SEDI(NBOX,ND,MD,MDT,RHOPAR,ZNOT,     &
       DTC,WETDP,USTR,PMID,PUPPER,PLOWER,T,SURTP,SEAICE,                &
       RHOA,MFPA,DVISC,BUD_AER_MAS,JLABOVE,ILSCAT,SEDI_ON,SM)

!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculates aerosol dry deposition and sedimentation.
! Based on the parameterisation of Zhang et al (2001) which
! uses the method in the model of Slinn (1982).
!
! Sedimentation is done using a simple explicit discretization
! which should be adequate for this process-split method.
!
! Evaluate deposition velocity in lowest level as:
!
! V_dep = V_g + 1/(A_r + S_r)
!
! where V_dep is the deposition velocity
!       V_g   is the gravitational velocity = rho_p*Dp^2*g*CF/(18*DVISC)
!       CF    is the Cunningham slip correction
!       DVISC is the dynamic viscosity
!       A_r   is the aerodynamic resitance
!       S_r   is the surface resistance
!       Dp    is the particle diameter
!       rho_p is the particle density
!       g     is the gravitational acceleration
!
! Evaluate S_r=1/{ 3 * ustar * (EB + EIM + EIN) }
!
! following parameterization by Zhang et al (2001) where
!
! EB,EIM,EIN are collection efficiencies for Brownian diffusion,
! impaction and interception respectively.
!
! EB = Sc^-YR where Sc is the particle Schmidt number = nu/D
!                                where nu = kinematic viscosity of air
!                                      D =  particle diffusion coeff.
!
!  and YR is surface-dependent constant, values as in Table 3 (Zhang01)
!         0.50 over water          (Land use category 13-14)
!         0.56 over forest         (Land use category  1- 5)
!         0.54 over grass and ice  (Land use category  6,12)
!
! EIM = { St/(ALPHA+St) }^2
!
!    where St is the Stokes number = V_g * ustar^2 / DVISC  (z0<1mm)
!                                  = V_g * ustar   / (g*CR) (z0>1mm)
!
!                                    [smooth & rough flow regimes]
!
!      and ALPHA,CR are surface-dependent constant, values as Table 3:
!         ALPHA=100.0, CR=0.0 over water [only divide by CR for veg]
!         ALPHA= 50.0, CR=0.0 over ice   [only divide by CR for veg]
!         ALPHA=  1.0, CR=0.005 over grass
!         ALPHA=  1.2, CR=0.002 over forest
!
! EIN = 0.5*Dp/CR
!
! Evaluates drydep & sedimentation for number & mass using 0th & 3rd
! order moment specific coefficients for modal aerosol as in Appendix 4
! Binkowski & Shankar (1995) JGR, vol 100, no D12, pp. 26,191--26,209.
!
! Note --- in this routine, sedimentation is included at all levels.
!
! Inputs :
! ------
! NBOX      : Number of grid boxes
! ND        : Initial no. concentration of aerosol mode (ptcls/cc)
! MD        : Avg cpt mass of aerosol ptcl in size mode (particle^-1)
! MDT       : Avg tot mass of aerosol ptcl in size mode (particle^-1)
! RHOPAR    : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
! ZNOT      : Roughness length (m)
! DTC       : Chemical timestep (s)
! WETDP     : Wet diameter for ptcl with dry diameter DRYDP (m)
! USTR      : Friction velocity(ms-1)
! PMID      : Centre level pressure (Pa)
! PUPPER    : Pressure at box upper interface (Pa)
! PLOWER    : Pressure at box lower interface (Pa)
! T         : Centre level temperature (K)
! SURTP     : Surface type [0=sea-surf,1=land-surf,2=above-surf]
! SEAICE    : Fraction of horizontal gridbox area containing seaice
! RHOA      : Air density (kg/m3)
! MFPA      : Mean free path of air (m)
! DVISC     : Dynamic viscosity of air (kg m-1 s-1)
! JLABOVE   : Index of box directly above this grid box
! ILSCAT    : Land-use category (1-9 based on UM landsurf types)
! SEDI_ON   : Switch for whether aerosol sedimentation is on/off
! SM        : Grid box mass of air (kg)
!
! Calls subroutine GETROUGH to read roughness length
!
! Outputs
! -------
! Updated particle number density ND (/cm3)
! Updated particle avg mass MD (molecules/particle)
! Updated Avg tot mass of aerosol ptcl in size mode MDT (particle^-1)
! BUD_AER_MAS: Aerosol mass budgets (mlcls/cc/tstep)
!
! Local Variables
! ---------------
! PS_AV_0    : 0th moment avg particle Schmidt Number
! PS_AV_3    : 3rd moment avg particle Schmidt Number
! KVISC      : Kinematic viscosity of air (m2 s-1)
! VGRAV_AV_0 : 0th moment avg grav. settling vel. (m s^-1)
! VGRAV_AV_3 : 3rd moment avg grav. settling vel. (m s^-1)
! VDEP_AV_0  : 0th moment avg deposition velocity (m s^-1)
! VDEP_AV_3  : 3rd moment avg deposition velocity (m s^-1)
! DCOEF_AV_0 : 0th moment avg particle diffusion coefficient(m2 s-1)
! DCOEF_AV_3 : 3rd moment avg particle diffusion coefficient(m2 s-1)
! SN_AV_0    : 0th moment avg Stokes number
! SN_AV_3    : 3rd moment avg Stokes number
! SR_AV_0    : 0th moment avg surface resistance
! SR_AV_3    : 3rd moment avg surface resistance
! EB_AV_0    : 0th moment avg collection eff. for Brownian diffusion
! EB_AV_3    : 3rd moment avg collection eff. for Brownian diffusion
! EIM_AV_0   : 0th moment avg collection eff. for impaction
! EIM_AV_3   : 3rd moment avg collection eff. for impaction
! EIN        : Collection eff. for interception
! AR         : Aerodynamic resistance
! MTOT       : Total aerosol mass conc [all cpts] (molecules/cm3)
! NEWN       : Updated number concentration (/cm3)
! DZ         : Ht difference between box vertical interfaces (m)
! DZMID      : Ht difference between box lower interface & mid-level (m)
! CR,Y,ALPHA: aerosol deposition coefficients
!     [vary with land category & input via DATA statements]
! CR        : Characteristic radius of collectors (m)
! Y         : Parameter for calculating Brownian diffusion
! ALPHA     : Parameter for calculating EIM
! MASK1     : Logical to define regions of domain for roughness categories
! MASK2     : Logical to define regions of domain for JLABOVE
! MASK3     : Logical to define regions with some aerosol in box or 1above
! MASK4     : Logical to define surface regions with aerosol in box/1above
! MASK_SMOO :Logical to define regsions over "smooth"  surface categories
! MASK_VEGE :Logical to define regsions over vegetated surface categories
! MASK_ABOVE_LIM: Logical to define boxes with VGRAV > allowed max value
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! GG        : Gravitational acceleration = 9.80665 ms^-2
! VKARMN    : Von Karman's constant = 0.4
! PPI       : 3.1415927
! RA        : Dry air gas constant = 287.05 Jkg^-1 K^-1
! ZBOLTZ    : Boltzman Constant (kg m2 s-2 K-1 molec-1)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! MODE      : Defines which modes are set
! COMPONENT : Defines which cpts are allowed in each mode
! SIGMAG    : Geometric standard deviation of mode
! MM        : Molar masses of components (kg/mole)
! RHOCOMP   : Densities (dry) of each component (kg/m^3)
! NUM_EPS   : Value of ND_0 below which do not recalculate MD (per cc)
!                                              or carry out process
! CP_SU     : Index of component in which H2SO4 is stored
! CP_BC     : Index of component in which BC is stored
! CP_OC     : Index of component in which 1st OC cpt is stored
! CP_CL     : Index of component in which NaCl is stored
! CP_SO     : Index of component in which 2nd OC cpt is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! Various indices for budget terms in BUD_AER_MAS
!
! References
! ----------
! Slinn, Atmos. En., 1982, 16, 1785-1794
! Zhang et al, Atmos. En., 2001, 35, 549-560
!
!----------------------------------------------------------------------
      USE UKCA_CONSTANTS,   ONLY: GG, VKARMN, PPI, RA, ZBOLTZ
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

! .. Subroutine interface
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: SEDI_ON
      INTEGER, INTENT(IN) :: JLABOVE(NBOX)
      INTEGER, INTENT(IN) :: ILSCAT(NBOX)
      REAL, INTENT(IN)    :: RHOPAR(NBOX,NMODES)
      REAL, INTENT(IN)    :: ZNOT(NBOX)
      REAL, INTENT(IN)    :: DTC
      REAL, INTENT(IN)    :: WETDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: USTR(NBOX)
      REAL, INTENT(IN)    :: PMID(NBOX)
      REAL, INTENT(IN)    :: PUPPER(NBOX)
      REAL, INTENT(IN)    :: PLOWER(NBOX)
      REAL, INTENT(IN)    :: T(NBOX)
      REAL, INTENT(IN)    :: SURTP(NBOX)
      REAL, INTENT(IN)    :: SEAICE(NBOX)
      REAL, INTENT(IN)    :: RHOA(NBOX)
      REAL, INTENT(IN)    :: MFPA(NBOX)
      REAL, INTENT(IN)    :: DVISC(NBOX)
      REAL, INTENT(IN)    :: SM(NBOX)
      REAL, INTENT(INOUT) :: BUD_AER_MAS(NBOX,0:NBUDAER)
      REAL, INTENT(INOUT) :: ND(NBOX,NMODES)
      REAL, INTENT(INOUT) :: MD(NBOX,NMODES,NCP)
      REAL, INTENT(INOUT) :: MDT(NBOX,NMODES)

! .. Local Variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      LOGICAL :: MASK1(NBOX)
      LOGICAL :: MASK2(NBOX)
      LOGICAL :: MASK3(NBOX)
      LOGICAL :: MASK4(NBOX)
      LOGICAL :: MASKSURF(NBOX)
      LOGICAL :: MASK_SMOO(NBOX)
      LOGICAL :: MASK_VEGE(NBOX)
      LOGICAL :: MASK_ABOVE_LIM(NBOX)
      REAL    :: PS_AV_0(NBOX)
      REAL    :: PS_AV_3(NBOX)
      REAL    :: KVISC(NBOX)
      REAL    :: VGRAV_AV_0(NBOX)
      REAL    :: VGRAV_AV_3(NBOX)
      REAL    :: VGRAV_AV_0_UP(NBOX)
      REAL    :: VGRAV_AV_3_UP(NBOX)
      REAL    :: DCOEF_AV_0(NBOX)
      REAL    :: DCOEF_AV_3(NBOX)
      REAL    :: EB_AV_0(NBOX)
      REAL    :: EB_AV_3(NBOX)
      REAL    :: EIM_AV_0(NBOX)
      REAL    :: EIM_AV_3(NBOX)
      REAL    :: EIN(NBOX)
      REAL    :: SN_AV_0(NBOX)
      REAL    :: SN_AV_3(NBOX)
      REAL    :: AR(NBOX)
      REAL    :: SR_AV_0(NBOX)
      REAL    :: SR_AV_3(NBOX)
      REAL    :: VDEP_AV_0(NBOX)
      REAL    :: VDEP_AV_3(NBOX)
      REAL    :: DZMID(NBOX)
      REAL    :: DZ(NBOX)
      REAL    :: DZ_UP(NBOX)
      REAL    :: T_UP(NBOX)
      REAL    :: PLO_UP(NBOX)
      REAL    :: PUP_UP(NBOX)
      REAL    :: SM_UP(NBOX)
      REAL    :: RHOA_UP(NBOX)
      REAL    :: ND0(NBOX,NMODES)
      REAL    :: MD0(NBOX,NMODES,NCP)
      REAL    :: ND0_UP(NBOX,NMODES)
      REAL    :: MD0_UP(NBOX,NMODES,NCP)
      REAL    :: NDNEW(NBOX)
      REAL    :: TERMIN_1(NBOX)
      REAL    :: TERMIN_2(NBOX)
      REAL    :: TERMIN_N(NBOX)
      REAL    :: TERMOUT_N(NBOX)
      REAL    :: TERMIN_M(NBOX,NCP)
      REAL    :: TERMOUT_M(NBOX,NCP)
      REAL    :: DELNSEDI(NBOX)
      REAL    :: DELMSEDI(NBOX,NCP)
      REAL    :: DELMDDEP(NBOX)
      REAL    :: VGRAV_LIM(NBOX)
      REAL    :: VGRAV_LIM_UP(NBOX)
!--------------------------------------------------------------------
!!      REAL, PARAMETER :: YR(9)    = (/0.57,0.56,0.54,0.54,0.54,         &
!!                                      0.56,0.50,0.54,0.54/)
!!      REAL, PARAMETER :: CR(9)    = (/5.0e-3,2.0e-3,2.0e-3,2.0e-3,0.01, &
!!                                      1.5e0 ,0.0e0 ,0.0e0 ,0.0e0/)
!!      REAL, PARAMETER :: ALPHA(9) = (/0.70,1.05,1.20,1.02,1.30,         &
!!                                      1.50,100.0,50.0,50.0/)
      REAL, PARAMETER :: YR(9)    = (/0.56,0.56,0.54,0.54,0.54,         &
                                      0.56,0.50,0.54,0.54/)
      REAL, PARAMETER :: CR(9)    = (/5.0e-3,5.0e-3,2.0e-3,2.0e-3,0.01, &
                                      1.5e0 ,0.0e0 ,0.0e0 ,0.0e0/)
      REAL, PARAMETER :: ALPHA(9) = (/1.00,1.00,1.20,1.02,1.30,         &
                                      1.50,100.0,50.0,50.0/)
!! here changed land-use categories 1/2 (trees) to match values used
!! for forest in 5-LS-category representation (should match then).
!
! Now set to match 9 UM Land-use types (ILSCAT)      YR  ALPHA  CR
!! 1=BL tree  (Zhang cat=avg of 2,4) [Evrgrn,Dec BL] 0.57  0.70 0.005
!! 2=NL tree  (Zhang cat=avg of 1,3) [Evrgrn,Dec NL] 0.56  1.05 0.002
! 1=BL tree  (Zhang cat=avg of 2,4) [Evrgrn,Dec BL] 0.56  1.00 0.005
! 2=NL tree  (Zhang cat=avg of 1,3) [Evrgrn,Dec NL] 0.56  1.00 0.005
! 3=C3 grass (Zhang cat=6)          [grass        ] 0.54  1.20 0.002
! 4=C4 grass (Zhang cat=6)          [grass        ] 0.54  1.20 0.002
! 5=Shrub    (Zhang cat=10)         [shrub i. wood] 0.54  1.30 0.010
! 6=Urban    (Zhang cat=15)         [urban        ] 0.56  1.50 1.500
! 7=Water    (Zhang cat=13/14)      [inl wat/ocean] 0.50 100.0 0.000
! 8=Soil     (Zhang cat=8)          [desert       ] 0.54 50.00 0.000
! 9=Ice      (Zhang cat=12)         [ice cap/glac.] 0.54 50.00 0.000
!
! ILSCAT is now an input to this routine and passed in to UKCA_AERO_STEP
!--------------------------------------------------------------------
!!!      REAL, PARAMETER :: YR(5) = (/  0.5, 0.56, 0.54,0.0,0.54/)
!!!      REAL, PARAMETER :: CR(5) = (/ 0.00,5.E-3,2.E-3,0.0,0.00/)
!!!      REAL, PARAMETER :: ALPHA(5) = (/100.0,  1.0,  1.2,0.0,50.0/)
!
! Previously used code below to set ILSCAT based on roughness length
! as either 1 = water/sea, 2= forest, 3 = all other land
!           4 = desert (not used), 5 = sea ice
! and used values from Zhang et al (2001) for YR (gamma in paper),
!                                             ALPHA
!                                             CR (A in paper)
! BELOW IS OLD TREATMENT                            YR   ALPHA    CR
! 1=water  (Zhang cat=13/14)[inland water/ocean  ] 0.50  100.0   0.000
! 2=forest (Zhang cat=1    )[evergreen needleleaf] 0.56    1.0   0.005
! 3=o land (Zhang cat=6/7  )[grass/crops         ] 0.54    1.2   0.002
! 4=desert (Zhang cat=8    )[desert              ] 0.54   50.0   0.000
! 5=seaice (Zhang cat=12   )[ice cap & glacier   ] 0.54   50.0   0.000

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_DDEPAER_INCL_SEDI',zhook_in,zhook_handle)

! Find out what category (water,forest,grass,desert)
! based on roughness length znot (desert not used at present)
! This should be updated in later version to read land type
! category directly. Desert not represented here.
!
!!!! water/sea - z0<0.001m
!!!      MASK1(:)=(ZNOT(:) < 1.0E-3) ! water/sea
!!!      WHERE(MASK1(:)) ILSCAT(:)=1
!!!
!!!! forests - z0>0.1m
!!!      MASK1(:)=(ZNOT(:) > 1.0E-1) ! forest
!!!      WHERE(MASK1(:)) ILSCAT(:)=2
!!!
!!!! all other lands, grass 0.001<z0<0.1m
!!!      MASK1(:)=((ZNOT(:) >= 1.0E-3).AND.(ZNOT(:) <= 1.0E-1)) ! grass
!!!      WHERE(MASK1(:)) ILSCAT(:)=3
!!!
!!!! If sea ice covers > 50% of sea surface, treat as sea ice
!!!      MASK1(:)=(SEAICE(:) > 0.5) ! seaice
!!!      WHERE(MASK1(:)) ILSCAT(:)=5
!--------------------------------------------------------------------
!
! .. First copy original values of ND,MD to ND0,MD0
      ND0(:,:)=ND(:,:)
      MD0(:,:,:)=MD(:,:,:)

      MASK2(:)=(JLABOVE(:) > 0) ! where using JLABOVE

      DO IMODE=1,NMODES
       WHERE(MASK2(:))
        ND0_UP(:,IMODE)=ND0(JLABOVE(:),IMODE)
       ELSEWHERE
        ND0_UP(:,IMODE)=0.0
       ENDWHERE
       DO ICP=1,NCP
        IF(COMPONENT(IMODE,ICP)) THEN
         WHERE(MASK2(:))
          MD0_UP(:,IMODE,ICP)=MD0(JLABOVE(:),IMODE,ICP)
         ELSEWHERE
          MD0_UP(:,IMODE,ICP)=MD0(:,IMODE,ICP)
         ENDWHERE
        END IF
       END DO
      END DO

      WHERE(MASK2(:))
       T_UP(:)=T(JLABOVE(:))
       PLO_UP(:)=PLOWER(JLABOVE(:))
       PUP_UP(:)=PUPPER(JLABOVE(:))
       SM_UP(:)=SM(JLABOVE(:))
       RHOA_UP(:)=RHOA(JLABOVE(:))
      ELSEWHERE
       T_UP(:)=T(:)
       PLO_UP(:)=PLOWER(:)
       PUP_UP(:)=PUPPER(:)
       SM_UP(:)=SM(:)
       RHOA_UP(:)=RHOA(:)
      ENDWHERE

      DZ(:)=(RA*T(:)/GG)*LOG(PLOWER(:)/PUPPER(:))
      DZMID(:)=(RA*T(:)/GG)*LOG(PLOWER(:)/PMID(:))

      DZ_UP(:)=(RA*T_UP(:)/GG)*LOG(PLO_UP(:)/PUP_UP(:))

      KVISC(:)=DVISC(:)/RHOA(:) ! Calculate kinematic viscosity of air

      MASKSURF(:)=(SURTP(:) < 2.0) ! create mask for boxes at surface.

! .. Calculate aerodynamic resistance (only in gridboxes at surface)
      WHERE(MASKSURF(:)) AR(:)=LOG(DZMID(:)/ZNOT(:))/(VKARMN*USTR(:))

      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN

! DEPENDS ON: ukca_vgrav_av_k
! .. Calculate 0th moment avg. grav. settling velocities
        CALL UKCA_VGRAV_AV_K(NBOX,0,WETDP(:,IMODE),SIGMAG(IMODE),       &
               DVISC(:),MFPA(:),RHOPAR(:,IMODE),VGRAV_AV_0(:))
! DEPENDS ON: ukca_vgrav_av_k
! .. Calculate 3rd moment avg. grav. settling velocities
        CALL UKCA_VGRAV_AV_K(NBOX,3,WETDP(:,IMODE),SIGMAG(IMODE),       &
               DVISC(:),MFPA(:),RHOPAR(:,IMODE),VGRAV_AV_3(:))
!
! .. Store values from box above so that can get flux into box
        WHERE(MASK2(:))
         VGRAV_AV_0_UP(:)=VGRAV_AV_0(JLABOVE(:))
         VGRAV_AV_3_UP(:)=VGRAV_AV_3(JLABOVE(:))
        ELSEWHERE
         VGRAV_AV_0_UP(:)=0.0
         VGRAV_AV_3_UP(:)=0.0
        ENDWHERE

! .. If sedimentation switched off, set all calculated VGRAV to zero
        IF(SEDI_ON == 0) THEN
         VGRAV_AV_0(:)=0.0
         VGRAV_AV_3(:)=0.0
         VGRAV_AV_0_UP(:)=0.0
         VGRAV_AV_3_UP(:)=0.0
        END IF

!       Calculate 0th moment avg particle diffusion coeffs
        CALL UKCA_DCOFF_PAR_AV_K(NBOX,0,WETDP(:,IMODE),SIGMAG(IMODE),   &
                  T(:),DVISC(:),MFPA(:),DCOEF_AV_0(:))

!       Calculate 3rd moment avg particle diffusion coeffs
        CALL UKCA_DCOFF_PAR_AV_K(NBOX,3,WETDP(:,IMODE),SIGMAG(IMODE),   &
                  T(:),DVISC(:),MFPA(:),DCOEF_AV_3(:))

! .. only calculate Schmidt number and collection eff at surface
        WHERE(MASKSURF(:))

!       Calculate 0th and 3rd moment avg. particle Schmidt number
         PS_AV_0(:)=KVISC(:)/DCOEF_AV_0(:)
         PS_AV_3(:)=KVISC(:)/DCOEF_AV_3(:)
!       Calculate particle collection efficiencies
!       -- For Brownian Diffusion
         EB_AV_0(:)=PS_AV_0(:)**(-YR(ILSCAT(:)))
         EB_AV_3(:)=PS_AV_3(:)**(-YR(ILSCAT(:)))
!
        ENDWHERE

! In new version, using 9 UM landsurf types,
! Set smooth surfaces to be water (7), soil (8) or ice (9)
! All other surfaces are vegetated (have CR>0)
        MASK_SMOO=( (ILSCAT(:) >= 7).AND.(ILSCAT(:) <= 9) )
        MASK_VEGE=.NOT.MASK_SMOO(:)

!!! Below is as in previous version using 5 landsurf types
!!        MASK_SMOO=( (ILSCAT(:) == 1).OR.(ILSCAT(:) == 5) )
!!        MASK_VEGE=( (ILSCAT(:) == 2).OR.(ILSCAT(:) == 3) )

! .. only calculate Stokes number at surface
        WHERE(MASK_SMOO(:).AND.MASKSURF(:))
!        Calculate stokes number for smooth surfaces
         SN_AV_0(:)=VGRAV_AV_0(:)*USTR(:)*USTR(:)/DVISC(:)
         SN_AV_3(:)=VGRAV_AV_3(:)*USTR(:)*USTR(:)/DVISC(:)
        ENDWHERE
        WHERE(MASK_VEGE(:).AND.MASKSURF(:))
!        Calculate stokes number for vegetated surfcaes
         SN_AV_0(:)=VGRAV_AV_0(:)*USTR(:)/(GG*CR(ILSCAT(:)))
         SN_AV_3(:)=VGRAV_AV_3(:)*USTR(:)/(GG*CR(ILSCAT(:)))
        ENDWHERE

! .. only calculate impaction collection efficiency at surface
        WHERE(MASKSURF(:))

!       -- For Impaction
         EIM_AV_0(:)=(SN_AV_0(:)/(ALPHA(ILSCAT(:))+SN_AV_0(:)))**2
         EIM_AV_3(:)=(SN_AV_3(:)/(ALPHA(ILSCAT(:))+SN_AV_3(:)))**2

        ENDWHERE

! .. only calculate interception collection eff (smooth) at surface
        WHERE(MASK_SMOO(:).AND.MASKSURF(:))

!       -- For Interception (smooth surfaces)
         EIN(:)=0.0

        ENDWHERE

! .. only calculate interception collection eff (vegetd) at surface
        WHERE(MASK_VEGE(:).AND.MASKSURF(:))

!       -- For Interception (vegetd surfaces)
         EIN(:)=0.5*(WETDP(:,IMODE)*WETDP(:,IMODE)                      &
                              /CR(ILSCAT(:))/CR(ILSCAT(:)))

        ENDWHERE

! .. only calculate surface resistance and deposition vel. at surface
! .. this section also increases VGRAV due to dry dep (at surface)
        WHERE(MASKSURF(:))

!        Calculate surface resistance
         SR_AV_0(:)=1.0/(3.0*USTR(:)*(EB_AV_0(:)+EIM_AV_0(:)+EIN(:)))
         SR_AV_3(:)=1.0/(3.0*USTR(:)*(EB_AV_3(:)+EIM_AV_3(:)+EIN(:)))

!        Calculate deposition velocity
         VDEP_AV_0(:)=VGRAV_AV_0(:)+1.0/(AR(:)+SR_AV_0(:))
         VDEP_AV_3(:)=VGRAV_AV_3(:)+1.0/(AR(:)+SR_AV_3(:))

!        Set gravitational velocity to deposition velocity if in lowest box
         VGRAV_AV_0(:)=VDEP_AV_0(:)
         VGRAV_AV_3(:)=VDEP_AV_3(:)

!        VGRAV_AV_UP never at surface so no need to set to dep. vel.

        ENDWHERE

! .. limit 0&3 grav. settling vel so only falls 0.5 box max [numerical]
        VGRAV_LIM(:)=0.5*DZ(:)/DTC

        MASK_ABOVE_LIM(:)=(VGRAV_AV_0(:) > VGRAV_LIM(:))
        WHERE(MASK_ABOVE_LIM(:)) VGRAV_AV_0(:)=VGRAV_LIM(:)

        MASK_ABOVE_LIM(:)=(VGRAV_AV_3(:) > VGRAV_LIM(:))
        WHERE(MASK_ABOVE_LIM(:)) VGRAV_AV_3(:)=VGRAV_LIM(:)

! .. limit 0&3 grav. settling vel so only falls 0.5 box max [box above]
        VGRAV_LIM(:)=0.5*DZ_UP(:)/DTC

        MASK_ABOVE_LIM(:)=(VGRAV_AV_0_UP(:) > VGRAV_LIM(:))
        WHERE(MASK_ABOVE_LIM(:)) VGRAV_AV_0_UP(:)=VGRAV_LIM(:)

        MASK_ABOVE_LIM(:)=(VGRAV_AV_3_UP(:) > VGRAV_LIM(:))
        WHERE(MASK_ABOVE_LIM(:)) VGRAV_AV_3_UP(:)=VGRAV_LIM(:)

! .. Calculate sedimenting in term (number) using 0th moment VGRAV_AV_0
        WHERE(MASK2(:)) ! where not top model level
         TERMIN_1(:)=ND0_UP(:,IMODE)/DZ_UP(:)
         TERMIN_2(:)=SM_UP(:)*RHOA(:)/SM(:)/RHOA_UP(:)
         TERMIN_N(:)=TERMIN_1(:)*TERMIN_2(:)*DTC*VGRAV_AV_0_UP(:)
        ELSEWHERE
!    If top model level TERMIN_N(:)=0.0
         TERMIN_N(:)=0.0
        ENDWHERE

        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          WHERE(MASK2(:))
           TERMIN_1(:)=ND0_UP(:,IMODE)*MD0_UP(:,IMODE,ICP)/DZ_UP(:)
           TERMIN_M(:,ICP)=TERMIN_1(:)*TERMIN_2(:)*DTC*VGRAV_AV_3_UP(:)
          ELSEWHERE
!    If top model level TERMIN_M(:,ICP)=0.0
           TERMIN_M(:,ICP)=0.0
          ENDWHERE
         END IF
        END DO

! .. Calculate sedimenting out term (number) using 0th moment VGRAV_AV_0
        TERMOUT_N(:)=(ND0(:,IMODE)/DZ(:))*DTC*VGRAV_AV_0(:)

! .. Calculate sedimenting in term (mass  ) using 3rd moment VGRAV_AV_3
! .. Calculate sedimenting out term (mass  ) using 3rd moment VGRAV_AV_3
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          TERMOUT_M(:,ICP)=(ND0(:,IMODE)*MD0(:,IMODE,ICP)/DZ(:))        &
                        *DTC*VGRAV_AV_3(:)
         END IF
        END DO

! .. below calculates net change in number and mass concentration
        DELNSEDI(:)=-(TERMIN_N(:)-TERMOUT_N(:))
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          DELMSEDI(:,ICP)=-(TERMIN_M(:,ICP)-TERMOUT_M(:,ICP))
         END IF
        END DO

! .. below masks boxes with some updates to apply
        MASK3(:)=( (ND0   (:,IMODE) > NUM_EPS(IMODE)) .OR.              &
                   (ND0_UP(:,IMODE) > NUM_EPS(IMODE)) )
        MASK4(:)=MASK3(:).AND.MASKSURF(:) ! some ptcls & also at surface

! .. below sets NDNEW to new value and re-sets MDT in boxes to update
        WHERE(MASK3(:)) ! only do where some particles in/out
         NDNEW(:)=ND0(:,IMODE)-DELNSEDI(:)
         MDT(:,IMODE)=0.0
        ENDWHERE

! .. below updates component masses and MDT for each mode
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          WHERE(MASK3(:)) ! only do where some particles in/out
           MD(:,IMODE,ICP)=                                             &
             (MD0(:,IMODE,ICP)*ND0(:,IMODE)-DELMSEDI(:,ICP))/NDNEW(:)
           MDT(:,IMODE)=MDT(:,IMODE)+MD(:,IMODE,ICP)
          ENDWHERE
         END IF ! IF COMPONENT(ICP)
        END DO ! loop over cpts

! .. below updates number concentration to NDNEW
! ..  (only do where some particles in/out)
        WHERE(MASK3(:)) ND(:,IMODE)=NDNEW(:)

! .. below stores ddep/sedi fluxes to BUD_AER_MAS
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          WHERE(MASK4(:)) DELMDDEP(:)=TERMOUT_M(:,ICP)
          IF(ICP == CP_SU) THEN
           IF((IMODE == 1).AND.(NMASDDEPSUNUCSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSUNUCSOL)=            &
                          BUD_AER_MAS(:,NMASDDEPSUNUCSOL)+DELMDDEP(:)
           IF((IMODE == 2).AND.(NMASDDEPSUAITSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSUAITSOL)=            &
                          BUD_AER_MAS(:,NMASDDEPSUAITSOL)+DELMDDEP(:)
           IF((IMODE == 3).AND.(NMASDDEPSUACCSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSUACCSOL)=            &
                          BUD_AER_MAS(:,NMASDDEPSUACCSOL)+DELMDDEP(:)
           IF((IMODE == 4).AND.(NMASDDEPSUCORSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSUCORSOL)=            &
                          BUD_AER_MAS(:,NMASDDEPSUCORSOL)+DELMDDEP(:)
          END IF
          IF(ICP == CP_BC) THEN
           IF((IMODE == 2).AND.(NMASDDEPBCAITSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPBCAITSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPBCAITSOL)+DELMDDEP(:)
           IF((IMODE == 3).AND.(NMASDDEPBCACCSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPBCACCSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPBCACCSOL)+DELMDDEP(:)
           IF((IMODE == 4).AND.(NMASDDEPBCCORSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPBCCORSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPBCCORSOL)+DELMDDEP(:)
           IF((IMODE == 5).AND.(NMASDDEPBCAITINS > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPBCAITINS)=            &
                       BUD_AER_MAS(:,NMASDDEPBCAITINS)+DELMDDEP(:)
          END IF
          IF(ICP == CP_OC) THEN
           IF((IMODE == 1).AND.(NMASDDEPOCNUCSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPOCNUCSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPOCNUCSOL)+DELMDDEP(:)
           IF((IMODE == 2).AND.(NMASDDEPOCAITSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPOCAITSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPOCAITSOL)+DELMDDEP(:)
           IF((IMODE == 3).AND.(NMASDDEPOCACCSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPOCACCSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPOCACCSOL)+DELMDDEP(:)
           IF((IMODE == 4).AND.(NMASDDEPOCCORSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPOCCORSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPOCCORSOL)+DELMDDEP(:)
           IF((IMODE == 5).AND.(NMASDDEPOCAITINS > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPOCAITINS)=            &
                       BUD_AER_MAS(:,NMASDDEPOCAITINS)+DELMDDEP(:)
          END IF
          IF(ICP == CP_CL) THEN
           IF((IMODE == 3).AND.(NMASDDEPSSACCSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSSACCSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPSSACCSOL)+DELMDDEP(:)
           IF((IMODE == 4).AND.(NMASDDEPSSCORSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSSCORSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPSSCORSOL)+DELMDDEP(:)
          END IF
          IF(ICP == CP_SO) THEN
           IF((IMODE == 1).AND.(NMASDDEPSONUCSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSONUCSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPSONUCSOL)+DELMDDEP(:)
           IF((IMODE == 2).AND.(NMASDDEPSOAITSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSOAITSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPSOAITSOL)+DELMDDEP(:)
           IF((IMODE == 3).AND.(NMASDDEPSOACCSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSOACCSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPSOACCSOL)+DELMDDEP(:)
           IF((IMODE == 4).AND.(NMASDDEPSOCORSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPSOCORSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPSOCORSOL)+DELMDDEP(:)
          END IF
          IF(ICP == CP_DU) THEN
           IF((IMODE == 3).AND.(NMASDDEPDUACCSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPDUACCSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPDUACCSOL)+DELMDDEP(:)
           IF((IMODE == 4).AND.(NMASDDEPDUCORSOL > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPDUCORSOL)=            &
                       BUD_AER_MAS(:,NMASDDEPDUCORSOL)+DELMDDEP(:)
           IF((IMODE == 6).AND.(NMASDDEPDUACCINS > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPDUACCINS)=            &
                       BUD_AER_MAS(:,NMASDDEPDUACCINS)+DELMDDEP(:)
           IF((IMODE == 7).AND.(NMASDDEPDUCORINS > 0))                  &
            WHERE(MASK4(:)) BUD_AER_MAS(:,NMASDDEPDUCORINS)=            &
                       BUD_AER_MAS(:,NMASDDEPDUCORINS)+DELMDDEP(:)
          END IF
         END IF ! if component present in mode
        END DO ! loop over components
       END IF ! IF MODE(IMODE)
      END DO ! loop over modes

      IF (lhook) CALL dr_hook('UKCA_DDEPAER_INCL_SEDI',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_DDEPAER_INCL_SEDI
