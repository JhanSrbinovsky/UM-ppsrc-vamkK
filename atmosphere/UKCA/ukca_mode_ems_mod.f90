! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
MODULE UKCA_MODE_EMS_MOD

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!  Description:
!   Supplies primary number and mass aerosol emissions for GLOMAP-mode.
!
!  Method:
!   Calls the GLOMAP-mode emission routines in turn using the supplied
!   emission arrays for sulphate, carbonaceous aerosols, and dust, together
!   with the relevant diagnostic arrays for the sea-salt routine.
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

CONTAINS

! Subroutine Interface:
      SUBROUTINE UKCA_MODE_EMS(row_length, rows, model_levels, ndiv,    &
          parfrac, mass, land_fraction, area, rough_length,             &
          seaice_frac, aird, u_scalar_10m,                              &
          emanso2, embiomso2, emvolconso2, emvolexpso2, so2_high_level, &
          emc, emcbm, emcdu,                                            &
          iso2ems, i_mode_ss_scheme,                                    &
          primsu_on,primbcoc_on,primss_on,primdu_on,verbose,            &
          aer_mas_primsu, aer_mas_primbc,                               &
          aer_mas_primoc, aer_mas_primss, aer_mas_primdu,               &
          aer_num_primsu, aer_num_primcar,                              &
          aer_num_primss, aer_num_primdu)

      USE UKCA_MODE_SETUP,  ONLY: nmodes
      USE ereport_mod,      ONLY: ereport
      USE yomhook,          ONLY: lhook, dr_hook
      USE parkind1,         ONLY: jprb, jpim

      IMPLICIT NONE

!  .. Below are input variables
      INTEGER, INTENT(IN) :: row_length             
! No of columns
      INTEGER, INTENT(IN) :: rows                   
! No of rows
      INTEGER, INTENT(IN) :: model_levels           
! No of model levels
      INTEGER, INTENT(IN) :: ndiv                   
! No of dust divisions
      INTEGER, INTENT(IN) :: iso2ems                
!
      INTEGER, INTENT(IN) :: so2_high_level         
! Level for injection of high emissions
      INTEGER, INTENT(IN) :: i_mode_ss_scheme       
!
      INTEGER, INTENT(IN) :: primsu_on              
! 1 for primary sulphate emissions
      INTEGER, INTENT(IN) :: primbcoc_on            
! 1 for primary BC/OC emissions
      INTEGER, INTENT(IN) :: primss_on              
! 1 for primary sea salt emissions
      INTEGER, INTENT(IN) :: primdu_on              
! 1 for primary dust emissions
      INTEGER, INTENT(IN) :: verbose                
! sets level of output

      REAL, INTENT(IN)    :: parfrac                
! fractn of SO2 mass emitted as primary SO4
      REAL, INTENT(IN)    :: mass(row_length,rows,model_levels)           
!
      REAL, INTENT(IN)    :: land_fraction(row_length,rows,model_levels)  
! land fraction
      REAL, INTENT(IN)    :: area(row_length,rows,model_levels)           
! surface area of grid
      REAL, INTENT(IN)    :: rough_length(row_length,rows,model_levels)   
! roughness length
      REAL, INTENT(IN)    :: seaice_frac(row_length,rows,model_levels)    
! sea ice fraction
      REAL, INTENT(IN)    :: aird(row_length,rows,model_levels)           
! No density /cm^3
      REAL, INTENT(IN)    :: u_scalar_10m(row_length,rows,model_levels)   
! scalar 10m wind
      REAL, INTENT(IN)    :: emanso2(row_length,rows,6)                   
! Anthrop. SO2 ems
      REAL, INTENT(IN)    :: embiomso2(row_length,rows,model_levels)      
! SO2 biomass ems
      REAL, INTENT(IN)    :: emvolconso2(row_length,rows,model_levels)    
! Constant Volc ems
      REAL, INTENT(IN)    :: emvolexpso2(row_length,rows,model_levels)    
! Explosive Volc ems
      REAL, INTENT(IN)    :: emc(row_length,rows,4)                       
! BC/OC bf/ff ems
      REAL, INTENT(IN)    :: emcbm(row_length,rows,model_levels,2)        
! BC/OC biomass ems
      REAL, INTENT(IN)    :: emcdu(row_length,rows,ndiv)                  
! dust emission rates (kg m-2 s-1)

! Output emission rates in kg [substance] /m^2/s for mass and
!  equivalent kg/m^2/s for number
      REAL, INTENT(OUT) :: aer_mas_primsu(row_length,rows,model_levels,nmodes)
      REAL, INTENT(OUT) :: aer_mas_primbc(row_length,rows,model_levels,nmodes)
      REAL, INTENT(OUT) :: aer_mas_primoc(row_length,rows,model_levels,nmodes)
      REAL, INTENT(OUT) :: aer_mas_primss(row_length,rows,model_levels,nmodes)
      REAL, INTENT(OUT) :: aer_mas_primdu(row_length,rows,model_levels,nmodes)
      REAL, INTENT(OUT) :: aer_num_primsu(row_length,rows,model_levels,nmodes)
      REAL, INTENT(OUT) :: aer_num_primcar(row_length,rows,model_levels,nmodes)
      REAL, INTENT(OUT) :: aer_num_primss(row_length,rows,model_levels,nmodes)
      REAL, INTENT(OUT) :: aer_num_primdu(row_length,rows,model_levels,nmodes)

! Local variables

      CHARACTER(LEN=72) :: cmessage     ! Error message
      INTEGER :: errcode                ! Error code

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_MODE_EMS',zhook_in,zhook_handle)

!  .. Call UKCA_PRIM_XX routines to return mass and number fluxes
 
 
      aer_mas_primsu(:,:,:,:) = 0.0
      aer_num_primsu(:,:,:,:) = 0.0

      aer_mas_primss(:,:,:,:) = 0.0
      aer_num_primss(:,:,:,:) = 0.0

      aer_mas_primbc(:,:,:,:) = 0.0
      aer_mas_primoc(:,:,:,:) = 0.0
      aer_num_primcar(:,:,:,:) = 0.0

      aer_num_primdu(:,:,:,:) = 0.0
      aer_mas_primdu(:,:,:,:) = 0.0

! .. Calculate primary sulphate emissions (call UKCA_PRIM_SU)
      IF (primsu_on == 1) THEN
 
! DEPENDS ON: ukca_prim_su
        CALL ukca_prim_su(row_length, rows, model_levels,               &
            so2_high_level, parfrac, emanso2, emvolconso2, emvolexpso2, &
            embiomso2, iso2ems, aer_mas_primsu, aer_num_primsu)

      END IF     !  primsu_on=1
      
! .. Calculate primary carbonaceous (BC/OC) aerosol emissions
      IF (primbcoc_ON == 1) THEN
 
! DEPENDS ON: ukca_prim_car
         CALL ukca_prim_car(row_length, rows, model_levels,             &
             emc, emcbm, iso2ems, aer_mas_primbc,                       &
             aer_mas_primoc, aer_num_primcar)

      END IF     !  primbcoc_on=1
 
! .. Calculate primary sea salt aerosol emissions
      IF (primss_on == 1) THEN
 
! DEPENDS ON: ukca_prim_ss
        CALL ukca_prim_ss(row_length, rows, model_levels,               &
            verbose, land_fraction, seaice_frac, mass, u_scalar_10m,    &
            aird, aer_mas_primss, aer_num_primss)

      END IF     ! primss_on=1

! .. Calculate primary dust aerosol emissions
      IF (primdu_on == 1) THEN
 
        cmessage = ' Primary dust emission routine not available'
        errcode = 1
        CALL EREPORT('UKCA_MODE_EMS',errcode,cmessage)
 
      END IF     !  primdu_on=1

      IF (lhook) CALL dr_hook('UKCA_MODE_EMS',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE UKCA_MODE_EMS
END MODULE UKCA_MODE_EMS_MOD
