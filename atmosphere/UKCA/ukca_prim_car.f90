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
!    Calculates primary carbonaceous aerosol emissions for GLOMAP-mode.
!
!  Method:
!    Use supplied carbonaceous emission arrays and combine these to give
!    modal black carbon and organic carbon mass emission fluxes with an
!    associated number emission.
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
      SUBROUTINE UKCA_PRIM_CAR(row_length, rows, model_levels,          &
                          emc, emcbm, iso2ems, aer_mas_primbc,          &
                          aer_mas_primoc, aer_num_primcar)
!-------------------------------------------------------------
!
!     Calculates primary carbonaceous aerosol emissions
!
!     Assumes Bond et al (2004) emissions (which give mass fluxes
!     of BC and OC for bio-fuel and fossil-fuel emissions on a
!     1x1 degree grid) are for particles consisting of an
!     internal mixture of BC and OC (rather than pure
!     BC particles and pure OC particles).
!
!     EMC  : Carbon Emissions (kgC/m2/s)
!      1: Black carbon   - Bio fuel
!      2: Black carbon   - Fossil fuel
!      3: Organic carbon - Bio fuel
!      4: Organic carbon - Fossil fuel
!
!     EMCBM : Carbon emissions from biomass burning (into heights) (kgC/m2/s)
!      1: Black carbon   - biomass burning
!      2: Organic carbon - biomass burning
!
!     Currently assume all BC and OC insoluble at emission
!
!     Inputs
!     -----
!     EMC        : BC/OC ems rates [bio- & fossil-fuel srcs] (kgC/m2/s)
!     EMCBM      : BC/OC ems rates [biomass burning srcs] (kgC/m2/s)
!     ISO2EMS    : Switch for choice of SO2 emissions
!
!     Outputs
!     -------
!     AER_MAS_PRIMBC  : Updated budget mass fluxes for BC emissions
!     AER_MAS_PRIMOC  : Updated budget mass fluxes for OC emissions
!     AER_NUM_PRIMCAR : Updated budget number fluxes for BC/OC emissions
!
!     Local variables
!     ---------------
!     MODE_DIAM  : Geom. mean diam. of modes into which emit bf/ff ems
!     STDEV      : Geom. std dev. of modes into which emit bf/ff ems
!     LSTDEV     : Natural logarithm of STDEV
!     EMCBC      : BC ems rate into this particular mode (kg/m2/s)
!     EMCOC      : OC ems rate into this particular mode (kg/m2/s)
!     EMCVOL     : Total (BC+OC) emitted particle volume (nm3/m2/s)
!     L          : Loops over 1 (bio-fuel & fire) and 2 (fossil-fuel)
!     FACTOR     : Converts from molecls or partcls per m2/s to kg/m2/s
!     DELN       : Change in particle number concentration in
!                  (equivalent dry air kg/m2/s)
!
!     Inputted by module UKCA_CONSTANTS
!     ---------------------------------
!     PPI        : 3.1415927
!     AVC        : Avogadro's constant (molecules per mole)
!     RA         : Dry air gas constant = 287.05 Jkg^-1 K^-1
!     ZBOLTZ     : Stefan-Boltzmann constant (kg m2 s-2 K-1 molec-1)
!     MM_DA      : Molar mass of dry air
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES     : Number of modes set
!     NCP        : Number of components set
!     MODE       : Logical variable defining which modes are set.
!     COMPONENT  : Logical variable defining which cpt are in which dsts
!     MM         : Molar masses of condensable components (kg/mole)
!     RHOCOMP    : Component mass densities (kg/m3)
!     CP_BC      : index of cpt into which emit BC ptcl mass
!     CP_OC      : index of cpt into which emit OC ptcl mass
!     FRACBCEM   : Fraction of BC ems to go into each mode
!     FRACOCEM   : Fraction of OC ems to go into each mode
!
!
!     References
!     ----------
!     Bond et al (2004), "A technology-based global inventory of black &
!       organic carbon emissions from combustion",
!       JGR Atmos, 109 (D14): Art. No. D14203.
!
!     Stier et al (2005), "The aerosol-climate model ECHAM5-HAM",
!       ACP, 5,  pp. 1125-1156.
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS,      ONLY: ppi, avc, ra, zboltz, mm_da
      USE UKCA_MODE_SETUP,     ONLY: nmodes, mode, cp_bc, cp_oc, mm,    &
                                     fracbcem, fracocem, ncp,           &
                                     component, rhocomp
      USE ereport_mod,         ONLY: ereport
      USE yomhook,             ONLY: lhook, dr_hook
      USE parkind1,            ONLY: jprb, jpim
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: row_length                          
! No of columns
      INTEGER, INTENT(IN) :: rows                                
! No of rows
      INTEGER, INTENT(IN) :: model_levels                        
! No of model levels
      INTEGER, INTENT(IN) :: iso2ems                             
! Switch for SO2 emissions
      REAL, INTENT(IN)    :: emc(row_length,rows,4)              
! Surface carbon emiss kgC/m2/s
      REAL, INTENT(IN)    :: emcbm(row_length,rows,model_levels,2)
                                                                 
! Biomass burn emiss kgC/m2/s
      REAL, INTENT(OUT)   :: aer_mas_primbc(row_length,rows,            &
                                             model_levels,nmodes)
                                                                 
! Mass fluxes for BC emissions
      REAL, INTENT(OUT)   :: aer_mas_primoc(row_length,rows,            &
                                             model_levels,nmodes)
                                                                 
! Mass fluxes for BC emissions
      REAL, INTENT(OUT)   :: aer_num_primcar(row_length,rows,           &
                                             model_levels,nmodes)

! Local variables
      REAL, PARAMETER  :: ocfact = 1.4    ! to convert from C to POM

      INTEGER :: imode                                  
! Loop index for modes
      INTEGER :: icp                                    
! Loop index for components
      INTEGER :: l                                      
! Loop index for sources
      INTEGER :: nsources                               
! Number of sources
      INTEGER :: errcode                                
! Error code
      REAL    :: emcvol(row_length,rows)                
! Total (BC+OC) emitted particle volume (nm3/m2/s)
      REAL    :: emcvol3d(row_length,rows,model_levels)
! Total (BC+OC) emitted particle volume (nm3/m2/s) (3d)
      REAL    :: factor                                 
!  converts from molecules or particles per m2/s to kg/m2/s
      REAL    :: emcbc(row_length,rows)                 
! BC emissn rate into mode (kg/m2/s)
      REAL    :: emcbc3d(row_length,rows,model_levels)  
! 3D ------------ " ---------------
      REAL    :: emcoc(row_length,rows)                 
! OC emissn rate into mode (kg/m2/s)
      REAL    :: emcoc3d(row_length,rows,model_levels)  
! 3D ------------ " ---------------
      REAL    :: deln(row_length,rows)                  
! Change in particle no conc (equivalent dry air kg/m2/s)
      REAL    :: deln3d(row_length,rows,model_levels)   
! Change in particle no conc (equiv-kg/m2/s)
      REAL    :: mode_diam(2)                           
! Geometric mean diameter of modes
      REAL    :: stdev(2)                               
! Geom. std dev. of modes
      REAL    :: lstdev(2)                              
! Log(stdev)
      CHARACTER(LEN=72) :: cmessage                     
! Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_PRIM_CAR',zhook_in,zhook_handle)

      IF(iso2ems >= 1 .AND. iso2ems <= 2) THEN
! If ISO2EMS 1 or 2 then use AEROCOM ACB recommendations for BC/OC size
        mode_diam=(/80.0,30.0/)
        stdev=(/1.8,1.8/)
        nsources=2
      ELSE IF ((iso2ems >= 3).OR.(iso2ems <= 6)) THEN
! If ISO2EMS == 3,4,5,6 use AEROCOM ACB recommendations for BC/OC size
!                       as modified by Stier et al (2005)
        mode_diam=(/150.0,60.0/)
        stdev=(/1.59,1.59/)
        nsources=2
      ELSE IF (iso2ems == 7) THEN
! If ISO2EMS == 7 then use AEROCOM ACB recommendations for BC/OC size
!                          as modified by Stier et al (2005)
! first value is used for wildfires, second value for all anthro emissions
        mode_diam=(/150.0,60.0/)
        stdev=(/1.59,1.59/)
        nsources=2
      ELSE
        cmessage = ' ISO2EMS is out of range'
        errcode = 1
        CALL EREPORT('UKCA_PRIM_CAR',errcode,cmessage)
      END IF


      factor=mm_da/avc   ! converts from molecls or partcls per m2/s to kg/m2/s

      DO imode=1,nmodes
        IF (mode(imode)) THEN

! If some BC/OC is emitted into this mode
          IF ((fracbcem(imode)+fracocem(imode)) > 0.0) THEN

            DO l=1,nsources ! loop over bio-fuel/fires(1) and fossil-fuel(2)
! For ISO2EMS<7
! l=1 accesses EMC(1) and EMC(3) (BC and OC) for bio-fuel & wildfires
! l=2 accesses EMC(2) and EMC(4) (BC and OC) for fossil-fuel
! For ISO2EMS=7
! l=1 accesses EMC(1) and EMC(3) (BC and OC) for wildfires
! l=2 accesses EMC(2) and EMC(4) (BC and OC) for biofuel & fossil-fuel
!-----------------------------------------------------------------------
! .. This section updates due to bio-fuel and fossil-fuel emissions

              lstdev(l) = LOG(stdev(l))

! .. Inject surface emissions (bio & fossil fuel) in kg/m2/s
              emcbc(:,:) = emc(:,:,l)*fracbcem(imode)
              emcoc(:,:) = emc(:,:,l+2)*fracocem(imode)*ocfact
! *ocfact to convert from C to POM

              deln(:,:) = 0.0
! .. Calculate emitted volume in nm3/m2/s
              emcvol(:,:) = 1e27*(emcbc(:,:)/                           &
                              rhocomp(cp_bc)+emcoc(:,:)/rhocomp(cp_oc))

! .. Then divide by mean particle volume and multiply by factor
! ..   to give equivalent (dry air) kg/m2/s 
              deln(:,:) = factor*emcvol(:,:)/((ppi/6.0)*                &
                              ((mode_diam(l))**3)*                    &
                              EXP(4.5*lstdev(l)*lstdev(l)) )


              DO icp=1,ncp
                IF(component(imode,icp)) THEN

! .. Update BC component mass
                  IF (icp == cp_bc) THEN

                    aer_mas_primbc(:,:,1,imode) =                       &
                        aer_mas_primbc(:,:,1,imode) +                   &
                           (factor*avc/mm(icp))*emcbc(:,:)
                  END IF

! .. Update OC component mass
                  IF (icp == cp_oc) THEN

                    aer_mas_primoc(:,:,1,imode) =                       &
                        aer_mas_primoc(:,:,1,imode) +                   &
                           (factor*avc/mm(icp))*emcoc(:,:)
                  END IF

                END IF    ! if component(imode,icp)
              END DO    ! loop over cpts

! .. Sum up each emitted number into mode via budget term (per equiv-kg/m2/s)
              aer_num_primcar(:,:,1,imode) =                            & 
                  aer_num_primcar(:,:,1,imode) + deln(:,:)

!-----------------------------------------------------------------------
! .. This section updates due to biomass burning BC/OC emissions

              IF(l == 1) THEN          ! If l=1, then do biomass burning

! .. Calculate BC emission rate into mode (kg/m2/s)
                emcbc3d(:,:,:) = emcbm(:,:,:,1)*fracbcem(imode)
! .. Calculate OC emission rate into mode (kg/m2/s) as POM
                emcoc3d(:,:,:) = emcbm(:,:,:,2)*fracocem(imode)*ocfact

                deln3d(:,:,:) = 0.0

! .. Calculate emitted particle volume (nm3/m2/s)
                emcvol3d(:,:,:) = 1e27*(emcbc3d(:,:,:)/rhocomp(cp_bc)   &
                                        + emcoc3d(:,:,:)/rhocomp(cp_oc))

! .. Divide by mean particle volume and multiply by factor to give
! ..  equivalent (as dry air) kg/m2/s
                deln3d(:,:,:) = factor*emcvol3d(:,:,:)/((ppi/6.0)*      &
                                  ((mode_diam(l))**3)*                &
                                  EXP(4.5*lstdev(l)*lstdev(l)) )
                  
                DO icp=1,ncp
                  IF(component(imode,icp)) THEN

! .. Update BC component mass
                    IF (icp == cp_bc) THEN

! .. sum up each emitted mass into mode
! ..   emcbc3d is in kg-cpt/m2/s
                      aer_mas_primbc(:,:,:,imode)=                      &
                        aer_mas_primbc(:,:,:,imode)+                    &
                          (factor*avc/mm(icp))*emcbc3d(:,:,:)      ! kg/m2/s

                    END IF      ! if icp=cp_bc

! .. Update OC component mass
                    IF (icp == cp_oc) THEN

! .. sum up each emitted mass into mode
                      aer_mas_primoc(:,:,:,imode)=                      &
                        aer_mas_primoc(:,:,:,imode)+                    &
                          (factor*avc/mm(icp))*emcoc3d(:,:,:)      ! kg/m2/s
! .. emcbc3d is in kg-cpt/m2/s
! .. delmdndcp is in kg-air/m2/s (mass-fraction-->mole-fraction)

                    END IF      ! if icp=cp_oc

                  END IF     ! component
                END DO     ! imode

! .. sum up each emitted number into mode via new budget term equiv-kg/m2/s
                aer_num_primcar(:,:,:,imode) =                          &
                     aer_num_primcar(:,:,:,imode) + deln3d(:,:,:)

              END IF    ! l=1


            END DO    ! l=1,nsources
          END IF   ! BC/OC emitted into this mode
        END IF ! if mode is defined
      END DO ! loop over nmodes

      IF (lhook) CALL dr_hook('UKCA_PRIM_CAR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_PRIM_CAR
