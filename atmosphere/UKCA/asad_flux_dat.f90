! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module containing the reactions to calculate the fluxes of
!
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!    This code is written to UMDP3 v8 programming standards.
!
! ######################################################################
!
! TP       = RXN for reaction
!            DEP for deposition
!            EMS for emission
!            NET for a tendency of a particular species
!            STE for the Strat-Trop exchange across the 2PVU+380K Tropopause
!                for a particular speices
!            MAS for the mass of air per gridcell
!            PSC for location of type 1/2 PSCs
! STASH#   = First 2 digits are section number, last 3 are item number,
!            5 digits in total, e.g. 34412 etc.
! WHICH    = IF TP==RXN: B = bimolecular         3D
!                        H = heterogeneous       3D
!                        T = termolecular        3D
!                        J = photolysis          3D
!            IF TP==DEP: D = dry deposition      2D
!                        W = wet deposition      3D
!            IF TP==EMS: S = surface emissions   2D
!                        A = aircraft emissions  3D
!                        L = lightning emissions 3D
!                        V = volcanic emissions  3D (only available if have SO2 tracer)
!                        S = 3D SO2   emissions  3D (only available if have SO2 tracer)
!            IF TP==NET  X = not used
!            IF TP==STE  X = not used
!            IF TP==MAS  X = not used 
!            IF TP==PSC  1/2 3D    
! RXN#     = If there are multiple reactions with the same reactants AND products
!            then this number specifies which of those (in order) you wish to budget.
!            Set to 0 otherwise.     
! #SPECIES = number of species. IF TP==RXN this is total of reactants + products
!                               IF TP!=RXN this should be 1 ONLY
! R1       = IF TP!=RXN this is the species being deposited/emitted
!
!
! NOTE: If the same stash number is specified for multiple diagnostics then these
!       will be summed up before outputting into stash for processing. Specifing
!       the same reaction multiple times is allowed, but only works for diagnostics
!       of the same dimensionallity, i.e. 2D *or* 3D
!
! UNITS OF DIAGNOSTICS ARE: moles/gridcell/s
!
! ######################################################################

MODULE asad_flux_dat

      USE UKCA_CHEM_DEFS_MOD,  ONLY: ratb_defs, ratt_defs,              &
                                     rath_defs, ratj_defs, chch_defs,   &
                                     chch_t, ratb_t, rath_t, ratj_t,    &
                                     ratt_t
      USE ASAD_MOD,            ONLY: ndepd, nldepd, ndepw, nldepw,      &
                                     advt, speci
      USE ereport_mod,         ONLY: ereport
      IMPLICIT NONE
  
      PRIVATE
  
! Describes the bimolecular reaction rates
      TYPE asad_flux_defn
      CHARACTER(LEN=3)  :: diag_type         ! which type of flux:RXN,DEP,EMS
      INTEGER           :: stash_number      ! stash number, e.g. 50001 etc.
      CHARACTER(LEN=1)  :: rxn_type          ! which rxn type: B,H,T,J,D,W,S,A,L,V
      LOGICAL           :: tropospheric_mask ! T or F
      INTEGER           :: rxn_location      ! for multiple reacions with the same
                                             ! reactants and products, default=0
      INTEGER           :: num_species       ! number of reactants+products
      CHARACTER(LEN=10) :: reactants(2)      ! list of reactants
      CHARACTER(LEN=10) :: products(4)       ! list of products
      ENDTYPE ASAD_FLUX_DEFN
      PUBLIC ASAD_FLUX_DEFN
  
      PUBLIC :: asad_chemical_fluxes

      CHARACTER(LEN=10) :: blank0 = '          '   ! Defines null product

! Number of chemical fluxes defined below
      INTEGER, PARAMETER :: n_chemical_fluxes = 232

      TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: asad_chemical_fluxes(:)

      TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_prod(20) = &
! Production of Ox
       (/ asad_flux_defn('RXN',50001,'B',.TRUE.,0,4,                    &
       (/'HO2       ','NO        '/),                                   &
       (/'OH        ','NO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50002,'B',.TRUE.,0,5,                       &
       (/'MeOO      ','NO        '/),                                   &
       (/'HO2       ','HCHO      ','NO2       ','          '/)),        &

! NO + RO2 reactions: sum into STASH section 34 item 303
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'EtOO      ','NO        '/),                                   &
       (/'MeCHO     ','HO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'MeCO3     ','NO        '/),                                   &
       (/'MeOO      ','CO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'n-PrOO    ','NO        '/),                                   &
       (/'EtCHO     ','HO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'i-PrOO    ','NO        '/),                                   &
       (/'Me2CO     ','HO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'EtCO3     ','NO        '/),                                   &
       (/'EtOO      ','CO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'MeCOCH2OO ','NO        '/),                                   &
       (/'MeCO3     ','HCHO      ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                       &
       (/'NO        ','ISO2      '/),                                   &
       (/'NO2       ','MACR      ','HCHO      ','HO2       '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                       &
       (/'NO        ','MACRO2    '/),                                   &
       (/'NO2       ','MeCO3     ','HACET     ','CO        '/)),        &

! OH + inorganic acid reactions: sum into STASH section 34 item 304
       asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                       &
       (/'OH        ','HONO2     '/),                                   &
       (/'H2O       ','NO3       ','          ','          '/)),        &
       asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                       &
       (/'OH        ','HONO      '/),                                   &
       (/'H2O       ','NO2       ','          ','          '/)),        &

! OH + organic nitrate reactions: sum into STASH section 34 item 305
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                       &
       (/'OH        ','MeONO2    '/),                                   &
       (/'HCHO      ','NO2       ','H2O       ','          '/)),        &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                       &
       (/'OH        ','NALD      '/),                                   &
       (/'HCHO      ','CO        ','NO2       ','          '/)),        &

! Organic nitrate photolysis: sum into STASH section 34 item 306
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                       &
       (/'MeONO2    ','PHOTON    '/),                                   &
       (/'HO2       ','HCHO      ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                       &
       (/'NALD      ','PHOTON    '/),                                   &
       (/'HCHO      ','CO        ','NO2       ','HO2       '/)),        &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                       &
       (/'ISON      ','PHOTON    '/),                                   &
       (/'NO2       ','MACR      ','HCHO      ','HO2       '/)),        &

! OH + PAN-type reactions: sum into STASH section 34 item 307
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                       &
       (/'OH        ','PAN       '/),                                   &
       (/'HCHO      ','NO2       ','H2O       ','          '/)),        &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                       &
       (/'OH        ','PPAN      '/),                                   &
       (/'MeCHO     ','NO2       ','H2O       ','          '/)),        &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,4,                       &
       (/'OH        ','MPAN      '/),                                   &
       (/'HACET     ','NO2       ','          ','          '/))         &
       /)


      TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_loss01(23) &
! Loss of Ox
       = (/ asad_flux_defn('RXN',50011,'B',.TRUE.,0,4,                  &
       (/'O(1D)     ','H2O       '/),                                   &
       (/'OH        ','OH        ','          ','          '/)),        &

! Minor reactions, should have negligble impact:
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(1D)     ','CH4       '/),                                   &
       (/'OH        ','MeOO      ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(1D)     ','CH4       '/),                                   &
       (/'HCHO      ','H2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,5,                       &
       (/'O(1D)     ','CH4       '/),                                   &
       (/'HCHO      ','HO2       ','HO2       ','          '/)),        &

! Include with above. Each is loss of 2xOx, so include twice to sum
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(3P)     ','O3        '/),                                   &
       (/'O2        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(3P)     ','O3        '/),                                   &
       (/'O2        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(3P)     ','NO2       '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(3P)     ','NO2       '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50013,'B',.TRUE.,0,4,                       &
       (/'HO2       ','O3        '/),                                   &
       (/'OH        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50014,'B',.TRUE.,0,4,                       &
       (/'OH        ','O3        '/),                                   &
       (/'HO2       ','O2        ','          ','          '/)),        &

! O3 + alkene:
! O3 + isoprene reactions
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                       &
       (/'O3        ','C5H8      '/),                                   &
       (/'MACR      ','HCHO      ','MACRO2    ','MeCO3     '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                       &
       (/'O3        ','C5H8      '/),                                   &
       (/'MeOO      ','HCOOH     ','CO        ','H2O2      '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                       &
       (/'O3        ','C5H8      '/),                                   &
       (/'HO2       ','OH        ','          ','          '/)),        &

! O3 + MACR reactions. ratb_defs specifies 2 different rates for each of these,
!  so need to select each one in turn
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'MGLY      ','HCOOH     ','HO2       ','CO        '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,4,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'OH        ','MeCO3     ','          ','          '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,2,6,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'MGLY      ','HCOOH     ','HO2       ','CO        '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,2,4,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'OH        ','MeCO3     ','          ','          '/)),        &

! N2O5 + H20 reaction
       asad_flux_defn('RXN',50016,'B',.TRUE.,0,4,                       &
       (/'N2O5      ','H2O       '/),                                   &
       (/'HONO2     ','HONO2     ','          ','          '/)),        &
! Heterogenious reactions are not included yet
!RXN    50016     H       4       N2O5    H2O    HONO2    HONO2           

! NO3 chemical loss
! sink of 2xOx:
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                       &
       (/'NO3       ','PHOTON    '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                       &
       (/'NO3       ','PHOTON    '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
! these are sinks of 1xOx                                                                          
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
       (/'HO2       ','NO3       '/),                                   &
       (/'OH        ','NO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
       (/'OH        ','NO3       '/),                                   &
       (/'HO2       ','NO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
       (/'MeOO      ','NO3       '/),                                   &
       (/'HO2       ','HCHO      ','NO2       ','          '/))         &
       /)

      TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_loss02(12) &
       = (/ asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                  &
       (/'EtOO      ','NO3       '/),                                   &
       (/'MeCHO     ','HO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
       (/'MeCO3     ','NO3       '/),                                   &
       (/'MeOO      ','CO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
       (/'n-PrOO    ','NO3       '/),                                   &
       (/'EtCHO     ','HO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
       (/'i-PrOO    ','NO3       '/),                                   &
       (/'Me2CO     ','HO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
       (/'EtCO3     ','NO3       '/),                                   &
       (/'EtOO      ','CO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
       (/'MeCOCH2OO ','NO3       '/),                                   &
       (/'MeCO3     ','HCHO      ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
       (/'NO3       ','HCHO      '/),                                   &
       (/'HONO2     ','HO2       ','CO        ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
       (/'NO3       ','MeCHO     '/),                                   &
       (/'HONO2     ','MeCO3     ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
       (/'NO3       ','EtCHO     '/),                                   &
       (/'HONO2     ','EtCO3     ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
       (/'NO3       ','Me2CO     '/),                                   &
       (/'HONO2     ','MeCOCH2OO ','          ','          '/)),        &
! sink of 2xOx                                                               
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                       &
       (/'NO3       ','C5H8      '/),                                   &
       (/'ISON      ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                       &
       (/'NO3       ','C5H8      '/),                                   &
       (/'ISON      ','          ','          ','          '/))         &
       /)

      TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_drydep(11) &
! Dry deposition:
! O3 dry deposition
       = (/asad_flux_defn('DEP',50021,'D',.TRUE.,0,1,                   &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),&
! NOy dry deposition
! sink of 2xOx
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! sink on 3xOx                                         
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'HO2NO2    ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'HONO2     ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'PAN       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'PPAN      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'MPAN      ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

      TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_wetdep(7)= &
! Wet deposition:
! NOy wet deposition
! sink of 2xOx
       (/ asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                    &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! sink on 3xOx                                         
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'HO2NO2    ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'HONO2     ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

      TYPE(asad_flux_defn), PARAMETER :: asad_trop_other_fluxes(16) =   &
! Extra fluxes of interest                                                                            
       (/ asad_flux_defn('RXN',50042,'B',.TRUE.,0,3,                    &
       (/'NO3       ','C5H8      '/),                                   &
       (/'ISON      ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50043,'B',.TRUE.,0,3,                       &
       (/'NO        ','ISO2      '/),                                   &
       (/'ISON      ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50044,'T',.TRUE.,0,4,                       &
       (/'HO2       ','HO2       '/),                                   &
       (/'H2O2      ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,3,                       &
       (/'HO2       ','HO2       '/),                                   &
       (/'H2O2      ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','MeOO      '/),                                   &
       (/'MeOOH     ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','EtOO      '/),                                   &
       (/'EtOOH     ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','MeCO3     '/),                                   &
       (/'MeCO3H    ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                       &
       (/'HO2       ','MeCO3     '/),                                   &
       (/'MeCO2H    ','O3        ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','n-PrOO    '/),                                   &
       (/'n-PrOOH   ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','i-PrOO    '/),                                   &
       (/'i-PrOOH   ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                       &
       (/'HO2       ','EtCO3     '/),                                   &
       (/'O2        ','EtCO3H    ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                       &
       (/'HO2       ','EtCO3     '/),                                   &
       (/'EtCO2H    ','O3        ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','MeCOCH2OO '/),                                   &
       (/'MeCOCH2OOH','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','ISO2      '/),                                   &
       (/'ISOOH     ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','MACRO2    '/),                                   &
       (/'MACROOH   ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50046,'T',.TRUE.,0,4,                       &
       (/'OH        ','NO2       '/),                                   &
       (/'HONO2     ','m         ','          ','          '/))         &
       /)
  
      TYPE(asad_flux_defn), PARAMETER :: asad_general_interest(8) = (/  &
! For CH4 lifetime
       asad_flux_defn('RXN',50041,'B',.TRUE.,0,4,                       &
       (/'OH        ','CH4       '/),                                   &
       (/'H2O       ','MeOO      ','          ','          '/)),        &
! Ozone STE
       asad_flux_defn('STE',50051,'X',.TRUE.,0,1,                       &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Ozone tendency in the troposphere
       asad_flux_defn('NET',50052,'X',.TRUE.,0,1,                       &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Ozone in the troposphere
       asad_flux_defn('OUT',50053,'X',.TRUE.,0,1,                       &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Ozone tendency - whole atmos
       asad_flux_defn('NET',50054,'X',.FALSE.,0,1,                      &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Tropospheric mass
       asad_flux_defn('MAS',50061,'X',.TRUE.,0,1,                       &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Tropospheric mask
       asad_flux_defn('TPM',50062,'X',.TRUE.,0,1,                       &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Total atmospheric mass
       asad_flux_defn('MAS',50063,'X',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /) 


      TYPE(asad_flux_defn), PARAMETER :: asad_trop_co_budget(21) = (/   &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,3,                       &
       (/'OH        ','CO        '/),                                   &
       (/'HO2       ','          ','          ','          '/)),        &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                       &
       (/'NO3       ','HCHO      '/),                                   &
       (/'HONO2     ','HO2       ','CO        ','          '/)),        &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                       &
       (/'OH        ','HCHO      '/),                                   &
       (/'H2O       ','HO2       ','CO        ','          '/)),        &
! MGLY + OH/NO3
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                       &
       (/'OH        ','MGLY      '/),                                   &
       (/'MeCO3     ','CO        ','          ','          '/)),        &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                       &
       (/'NO3       ','MGLY      '/),                                   &
       (/'MeCO3     ','CO        ','HONO2     ','          '/)),        &
! O3 + MACR/ISOP & OTHER FLUXES
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                       &
       (/'O3        ','C5H8      '/),                                   &
       (/'MeOO      ','HCOOH     ','CO        ','H2O2      '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,1,6,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'MGLY      ','HCOOH     ','HO2       ','CO        '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,2,6,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'MGLY      ','HCOOH     ','HO2       ','CO        '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                       &
       (/'NO        ','MACRO2    '/),                                   &
       (/'NO2       ','MeCO3     ','HACET     ','CO        '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                       &
       (/'MACRO2    ','MACRO2    '/),                                   &
       (/'HACET     ','MGLY      ','HCHO      ','CO        '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                       &
       (/'OH        ','NALD      '/),                                   &
       (/'HCHO      ','CO        ','NO2       ','          '/)),        &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                       &
       (/'HCHO      ','PHOTON    '/),                                   &
       (/'HO2       ','HO2       ','CO        ','          '/)),        &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                       &
       (/'HCHO      ','PHOTON    '/),                                   &
       (/'H2        ','CO        ','          ','          '/)),        &
! MGLY photolysis
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                       &
       (/'MGLY      ','PHOTON    '/),                                   &
       (/'MeCO3     ','CO        ','HO2       ','          '/)),        &
! OTHER CO PROD PHOTOLYSIS RXNS
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                       &
       (/'MeCHO     ','PHOTON    '/),                                   &
       (/'MeOO      ','HO2       ','CO        ','          '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,4,                       &
       (/'MeCHO     ','PHOTON    '/),                                   &
       (/'CH4       ','CO        ','          ','          '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                       &
       (/'EtCHO     ','PHOTON    '/),                                   &
       (/'EtOO      ','HO2       ','CO        ','          '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                       &
       (/'MACR      ','PHOTON    '/),                                   &
       (/'MeCO3     ','HCHO      ','CO        ','HO2       '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                       &
       (/'MACROOH   ','PHOTON    '/),                                   &
       (/'HACET     ','CO        ','MGLY      ','HCHO      '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                       &
       (/'NALD      ','PHOTON    '/),                                   &
       (/'HCHO      ','CO        ','NO2       ','HO2       '/)),        &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                       &
       (/'CO        ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

      TYPE(asad_flux_defn), PARAMETER :: asad_lightning_diags(6) = (/   &
! WARNING: Lightning emissions are calculated as NO2, but emitted as NO
       asad_flux_defn('EMS',50081,'L',.FALSE.,0,1,                      &
       (/'NO        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LGT',50082,'T',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LGT',50083,'G',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LGT',50084,'C',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LGT',50085,'N',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LIN',50086,'X',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

       TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                       &
                                            asad_strat_oh_prod(31) = (/ &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,5,                      &
       (/'Cl        ','HOCl      '/),                                   &
       (/'Cl        ','Cl        ','OH        ','          '/)),        & 
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(3P)     ','HBr       '/),                                  &
       (/'OH        ','Br        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(3P)     ','HCl       '/),                                  &
       (/'OH        ','Cl        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(3P)     ','HOCl      '/),                                  &
       (/'OH        ','ClO       ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(1D)     ','HBr       '/),                                  &
       (/'OH        ','Br        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(1D)     ','HCl       '/),                                  &
       (/'OH        ','Cl        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(1D)     ','MeBr      '/),                                  &
       (/'OH        ','Br        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'Cl        ','HO2       '/),                                  &
       (/'ClO       ','OH        ','          ','          '/)),        &   
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'Cl        ','HO2       '/),                                  &
       (/'ClO       ','OH        ','          ','          '/)),        & 
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'NO3       ','HO2       '/),                                  &
       (/'NO2       ','OH        ','          ','          '/)),        &   
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O3        ','HO2       '/),                                  &
       (/'O2        ','OH        ','          ','          '/)),        &   
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(3P)     ','HO2       '/),                                  &
       (/'OH        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(1D)     ','CH4       '/),                                  &
       (/'OH        ','MeOO      ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(1D)     ','H2        '/),                                  &
       (/'OH        ','H         ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(1D)     ','H2O       '/),                                  &
       (/'OH        ','OH        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(1D)     ','H2O       '/),                                  &
       (/'OH        ','OH        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'H         ','HO2       '/),                                  &
       (/'OH        ','OH        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'H         ','HO2       '/),                                  &
       (/'OH        ','OH        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'H         ','O3        '/),                                  &
       (/'OH        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'H         ','NO2       '/),                                  &
       (/'OH        ','NO        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(3P)     ','H2        '/),                                  &
       (/'OH        ','H         ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'O(3P)     ','H2O2      '/),                                  &
       (/'OH        ','HO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,5,                      &
        (/'O(3P)     ','HCHO      '/),                                  &
       (/'OH        ','CO        ','HO2       ','          '/)),        &
       asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
        (/'NO        ','HO2       '/),                                  &
       (/'NO2       ','OH        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
         (/'HOBr      ','PHOTON    '/),                                 &
       (/'OH        ','Br        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
         (/'HOCl      ','PHOTON    '/),                                 &
       (/'OH        ','Cl        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
         (/'HONO2     ','PHOTON    '/),                                 &
       (/'NO2       ','OH        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
         (/'HO2NO2    ','PHOTON    '/),                                 &
       (/'NO3       ','OH        ','          ','          '/)),        &
       asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
         (/'H2O       ','PHOTON    '/),                                 &
       (/'OH        ','H         ','          ','          '/)),        &          
       asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
         (/'H2O2      ','PHOTON    '/),                                 &
       (/'OH        ','OH        ','          ','          '/)),        &  
       asad_flux_defn('RXN',50091,'J',.FALSE.,0,5,                      &
         (/'MeOOH     ','PHOTON    '/),                                 &
       (/'HCHO      ','HO2       ','OH        ','          '/))         &
       /)

! OH loss reactions
     TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                         &
                                            asad_strat_oh_loss(26) = (/ &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'BrO       ','OH        '/),                                   &
       (/'Br        ','HO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'ClO       ','OH        '/),                                   &
       (/'Cl        ','HO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','ClO       '/),                                   &
       (/'HCl       ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','HBr       '/),                                   &
       (/'H2O       ','Br        ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','HCl       '/),                                   &
       (/'H2O       ','Cl        ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','HOCl      '/),                                   &
       (/'ClO       ','H2O       ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','OClO      '/),                                   &
       (/'HOCl      ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','ClONO2    '/),                                   &
       (/'HOCl      ','NO3       ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','MeBr      '/),                                   &
       (/'Br        ','H2O       ','          ','          '/)),        & 
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
        (/'OH        ','O3        '/),                                  &
       (/'HO2       ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'O(3P)     ','OH        '/),                                   &
       (/'O2        ','H         ','          ','          '/)),        & 
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','CH4       '/),                                   &
       (/'H2O       ','MeOO      ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,3,                      &
       (/'OH        ','CO        '/),                                   &
       (/'HO2       ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
        (/'OH        ','NO3       '/),                                  &
       (/'HO2       ','NO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                      &
       (/'OH        ','NO2       '/),                                   &
       (/'HONO2     ','m         ','          ','          '/)),        &
       Asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','HONO2     '/),                                   &
       (/'NO3       ','H2O       ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,5,                      &
       (/'OH        ','HO2NO2    '/),                                   &
       (/'H2O       ','NO2       ','O2        ','          '/)),        &
       asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                      &
        (/'OH        ','OH        '/),                                  &
       (/'H2O2      ','m         ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                      &
       (/'OH        ','OH        '/),                                   &
       (/'H2O2      ','m         ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','OH        '/),                                   &
       (/'H2O       ','O(3P)     ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','OH        '/),                                   &
       (/'H2O       ','O(3P)     ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,3,                      &
       (/'OH        ','HO2       '/),                                   &
       (/'H2O       ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
        (/'OH        ','H2O2      '/),                                  &
       (/'H2O       ','HO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,5,                      &
        (/'OH        ','HCHO      '/),                                  &
       (/'H2O       ','CO        ','HO2       ','          '/)),        &
       Asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','H2        '/),                                   &
       (/'H2O       ','H02       ','          ','          '/)),        &
       asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
       (/'OH        ','MeOOH     '/),                                   &
       (/'MeOO      ','H2O       ','          ','          '/))         &
       /)

! Simple strat ozone budget
       TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                       &
                                          asad_strat_o3_budget(20) = (/ &
       ! production
       asad_flux_defn('RXN',50101,'J',.FALSE.,0,4,                      &
       (/'O2        ','PHOTON    '/),                                   &
       (/'O(3P)     ','O(3P)     ','          ','          '/)),        &
       asad_flux_defn('RXN',50101,'J',.FALSE.,0,4,                      &
       (/'O2        ','PHOTON    '/),                                   &
       (/'O(3P)     ','O(1D)     ','          ','          '/)),        &
       asad_flux_defn('RXN',50102,'B',.FALSE.,0,4,                      &
       (/'HO2       ','NO        '/),                                   &
       (/'OH        ','NO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50103,'B',.FALSE.,0,5,                      &
       (/'MeOO      ','NO        '/),                                   &
       (/'HO2       ','HCHO      ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50104,'B',.FALSE.,0,4,                      &
       (/'OH        ','HONO2     '/),                                   &
       (/'H2O       ','NO3       ','          ','          '/)),        &
       ! loss
       asad_flux_defn('RXN',50111,'J',.FALSE.,0,5,                      &
       (/'Cl2O2     ','PHOTON    '/),                                   &
       (/'Cl        ','Cl        ','O2        ','          '/)),        &
       asad_flux_defn('RXN',50112,'B',.FALSE.,0,4,                      &
       (/'BrO       ','ClO       '/),                                   &
       (/'Br        ','OClO      ','          ','          '/)),        &
       asad_flux_defn('RXN',50113,'B',.FALSE.,0,4,                      &
       (/'HO2       ','O3        '/),                                   &
       (/'OH        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50114,'B',.FALSE.,0,4,                      &
       (/'ClO       ','HO2       '/),                                   &
       (/'HOCl      ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50115,'B',.FALSE.,0,4,                      &
       (/'BrO       ','HO2       '/),                                   &
       (/'HOBr      ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50116,'B',.FALSE.,0,4,                      &
       (/'O(3P)     ','ClO       '/),                                   &
       (/'Cl        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50117,'B',.FALSE.,0,4,                      &
       (/'O(3P)     ','NO2       '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50118,'B',.FALSE.,0,4,                      &
       (/'O(3P)     ','HO2       '/),                                   &
       (/'OH        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50119,'B',.FALSE.,0,4,                      &
       (/'O3        ','H         '/),                                   &
       (/'OH        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50120,'B',.FALSE.,0,4,                      &
       (/'O(3P)     ','O3        '/),                                   &
       (/'O2        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50121,'J',.FALSE.,0,4,                      &
       (/'NO3       ','PHOTON    '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50122,'B',.FALSE.,0,4,                      &
       (/'O(1D)     ','H2O       '/),                                   &
       (/'OH        ','OH        ','          ','          '/)),        &
       asad_flux_defn('RXN',50123,'B',.FALSE.,0,5,                      &
       (/'HO2       ','NO3       '/),                                   &
       (/'OH        ','NO2       ','O2        ','          '/)),        &
       asad_flux_defn('RXN',50124,'B',.FALSE.,0,4,                      &
       (/'OH        ','NO3       '/),                                   &
       (/'HO2       ','NO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50125,'B',.FALSE.,0,5,                      &
       (/'NO3       ','HCHO      '/),                                   &
       (/'HONO2     ','HO2       ','CO        ','          '/))         &
       /)

      TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                        &
                                            asad_strat_o3_misc(15) = (/ &
       asad_flux_defn('DEP',50131,'D',.FALSE.,0,1,                      &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
       (/'HO2NO2    ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
       (/'HONO2     ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
       (/'HO2NO2    ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
       (/'HONO2     ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

! Tropospheric sulphur chemistry
       TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                       &
                                          asad_aerosol_chem(16) = (/    &
       asad_flux_defn('RXN',50140,'B',.FALSE.,0,5,                      &
       (/'DMS       ','OH        '/),                                   &
       (/'SO2       ','MeOO      ','HCHO      ','          '/)),        &
       asad_flux_defn('RXN',50141,'B',.FALSE.,0,5,                      &
       (/'DMS       ','OH        '/),                                   &
       (/'SO2       ','DMSO      ','MeOO      ','          '/)),        &
       asad_flux_defn('RXN',50142,'B',.FALSE.,0,6,                      &
       (/'DMS       ','NO3       '/),                                   &
       (/'SO2       ','HONO2     ','MeOO      ','HCHO      '/)),        &
       asad_flux_defn('RXN',50143,'B',.FALSE.,0,4,                      &
       (/'DMSO      ','OH        '/),                                   &
       (/'SO2       ','MSA       ','          ','          '/)),        &
       asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                      &
       (/'CS2       ','OH        '/),                                   &
       (/'SO2       ','COS       ','          ','          '/)),        &
       asad_flux_defn('RXN',50145,'B',.FALSE.,0,3,                      &
       (/'H2S       ','OH        '/),                                   &
       (/'SO2       ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50146,'B',.FALSE.,0,3,                      &
       (/'COS       ','OH        '/),                                   &
       (/'SO2       ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                      &
       (/'Monoterp  ','OH        '/),                                   &
       (/'Sec_Org   ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                      &
       (/'Monoterp  ','O3        '/),                                   &
       (/'Sec_Org   ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                      &
       (/'Monoterp  ','NO3       '/),                                   &
       (/'Sec_Org   ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50150,'T',.FALSE.,0,4,                      &
       (/'SO2       ','OH        '/),                                   &
       (/'HO2       ','H2SO4     ','          ','          '/)),        &
       asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                      &
       (/'SO2       ','H2O2      '/),                                   &
       (/'NULL0     ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                      &
       (/'SO2       ','O3        '/),                                   &
       (/'NULL1     ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                      &
       (/'SO2       ','O3        '/),                                   &
       (/'NULL2     ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                       &
       (/'SO2       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                       &
       (/'SO2       ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

      PUBLIC :: ASAD_LOAD_DEFAULT_FLUXES
      INTERFACE ASAD_LOAD_DEFAULT_FLUXES
        MODULE PROCEDURE ASAD_LOAD_DEFAULT_FLUXES
      END INTERFACE ASAD_LOAD_DEFAULT_FLUXES

      CONTAINS

! #####################################################################

      SUBROUTINE ASAD_LOAD_DEFAULT_FLUXES
! To load flux information into asad_chemical_fluxes array for use in
!  asad flux diagnostics scheme. Note that both the total number of
!  fluxes and the number in each section (currently: reaction fluxes;
!  photolytic fluxes; dry deposition; and wet deposition) are limited.
!  Reactions etc are either taken from the definition arrays above
!  or are specified in a user-STASHmaster file, and read in using
!  the chemical definition arrays.
!  NOT currently handling the Strat-trop exchange as this should be
!  in section 38.
!  Called from UKCA_MAIN1.

      USE UKCA_D1_DEFS,        ONLY: ukca_sect, mode_diag_sect
      USE yomhook,             ONLY: lhook, dr_hook
      USE parkind1,            ONLY: jprb, jpim
      USE Control_Max_Sizes
      USE version_mod,         ONLY: nproftp, nprofdp, nprofup,         &
                                     ndiagpm, ntimep, NTimSerP,         &
                                     nlevp, npslevp, npslistp,          &
                                     outfile_s, outfile_e, nsectp,      &
                                     nitemp
      USE Submodel_Mod
      USE PrintStatus_mod

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
! COMDECK CSTASH
! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.
!
!
! Declarations:
! Imported global variables:
!    None, but see note above.

! Global parameters:

! Global scalars:
      INTEGER      NDIAG   ! No. of diagnostics
      INTEGER      NTPROF  ! No. of time profiles
      INTEGER      NSERIES ! No. of stash time series
      INTEGER      NDPROF  ! No. of domain profiles
      INTEGER      NUPROF  ! No. of useage profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in VERSION_MOD module
      INTEGER      MODL_B(NDIAGPM)  ! Internal model no.
      INTEGER      ISEC_B(NDIAGPM)  ! Section
      INTEGER      ITEM_B(NDIAGPM)  ! Item
      INTEGER      ITIM_B(NDIAGPM)  ! Time profile number
      INTEGER      IDOM_B(NDIAGPM)  ! Domain profile number
      INTEGER      IUSE_B(NDIAGPM)  ! Useage profile number

!   Time profile information:

      CHARACTER(LEN=8)  TIMPRO(NPROFTP)         ! Name of profile
      INTEGER      ITYP_T(NPROFTP)         ! Type of profile
      INTEGER      INTV_T(NPROFTP)         ! Time Interval
      CHARACTER(LEN=2)  UNT1_T(NPROFTP)         ! Units for time interval
      INTEGER      ISAM_T(NPROFTP)         ! Sampling period
      CHARACTER(LEN=2)  UNT2_T(NPROFTP)         ! Units for sampling period
      INTEGER      IOPT_T(NPROFTP)         ! Output option
      INTEGER      ISTR_T(NPROFTP)         ! Output Start time
      INTEGER      IEND_T(NPROFTP)         ! Output End time
      INTEGER      ISDT_T(6, NPROFTP)      ! Output Start date
      INTEGER      IEDT_T(6, NPROFTP)      ! Output End date
      INTEGER      IFRE_T(NPROFTP)         ! Output frequency
      INTEGER      IOFF_T(NPROFTP)         ! Offset for sampling
      CHARACTER(LEN=2)  UNT3_T(NPROFTP)         ! Units for output times
      INTEGER      ITIM_T(NPROFTP)         ! No. of times in times table
      INTEGER      ISER_T(NTIMEP ,NPROFTP) ! Times table (with units)
      INTEGER      MODL_T(NPROFTP)         ! Indicates internal model
                                           !  for each times table

!   Domain profile information:

      CHARACTER(LEN=8) DOMPRO  (NPROFDP)           ! Name of domain profile
      INTEGER     IOPL_D  (NPROFDP)           ! Levels option
      INTEGER     LEVB_D  (NPROFDP)           ! Bottom level
      INTEGER     LEVT_D  (NPROFDP)           ! Top level
      INTEGER     IOPA_D  (NPROFDP)           ! Area option
      INTEGER     INTH_D  (NPROFDP)           ! North boundary
      INTEGER     ISTH_D  (NPROFDP)           ! South boundary
      INTEGER     IEST_D  (NPROFDP)           ! East boundary
      INTEGER     IWST_D  (NPROFDP)           ! West boundary
      INTEGER     IMSK_D  (NPROFDP)           ! Mask type
      INTEGER     IMN_D   (NPROFDP)           ! Meaning option
      INTEGER     IWT_D   (NPROFDP)           ! Weighting option
      CHARACTER(LEN=1) TS_D    (NPROFDP)           ! Time series profile
      INTEGER     IG_TS
      INTEGER     I1_TS
      INTEGER     I51_TS
      INTEGER     BLIM_TS (NTimSerP)
      INTEGER     TLIM_TS (NTimSerP)
      REAL        BLIMR_TS(NTimSerP)
      REAL        TLIMR_TS(NTimSerP)
      INTEGER     NLIM_TS (NTimSerP)
      INTEGER     SLIM_TS (NTimSerP)
      INTEGER     ELIM_TS (NTimSerP)
      INTEGER     WLIM_TS (NTimSerP)
      INTEGER     ILEV_D  (NPROFDP)           ! Output levels code
      INTEGER     LEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      REAL       RLEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      INTEGER     PLT_D   (NPROFDP)
      INTEGER     PLLEN_D (NPROFDP)
      INTEGER     PLPOS_D (NPROFDP)
      INTEGER     PSLIST_D(NPSLEVP ,NPSLISTP)
      INTEGER     NPSLISTS
      EQUIVALENCE        (RLEVLST_D,LEVLST_D)

! Useage information:

      CHARACTER(LEN=8) USEPRO(NPROFUP)   ! Name of useage profile
      INTEGER     LOCN_U(NPROFUP)   ! Storage location of profile
      INTEGER     IUNT_U(NPROFUP)   ! Unit no.

! Information from ppxref file:

      INTEGER      MODEL_ST       ! Internal model number
      INTEGER      ISPACE         ! Space code
      INTEGER      ITIMA          ! Time availability code
      INTEGER      IGP            ! Grid of data code
      INTEGER      ILEV           ! Level type code
      INTEGER      IBOT           ! First level code
      INTEGER      ITOP           ! Last level code
      INTEGER      IFLAG          ! Level compression flag
      INTEGER      IOPN(6)        ! Sectional option code
      INTEGER      VMSK           ! Integer equiv of bin vers mask
      INTEGER      IPSEUDO        ! Pseudo dimension type
      INTEGER      IPFIRST        ! First pseudo dim code
      INTEGER      IPLAST         ! Last pseudo dim code
      INTEGER      PTR_PROG       ! Section zero point back
      INTEGER      HALO_TYPE      ! Type of halo the field has

! PP output file units
      INTEGER      PPlen2LkUp(OUTFILE_S:OUTFILE_E)
      CHARACTER(LEN=1)  FTOutUnit (OUTFILE_S:OUTFILE_E)

! COMMON blocks:
      COMMON/STCHA/ TIMPRO,UNT1_T,UNT2_T,UNT3_T,DOMPRO,TS_D,            &
     &  USEPRO,FTOutUnit

      COMMON/STSH/                                                      &
     &  NDIAG   ,MODL_B  ,ISEC_B ,ITEM_B  ,ITIM_B  ,IDOM_B  ,IUSE_B,    &
     &  NTPROF  ,ITYP_T  ,INTV_T ,ISAM_T  ,ITIM_T  ,                    &
     &  IOPT_T  ,ISTR_T  ,IEND_T ,IFRE_T  ,IOFF_T, ISER_T  ,MODL_T  ,   &
     &  NDPROF  ,IOPL_D  ,LEVB_D ,ISDT_T  ,IEDT_T  ,                    &
     &  IOPA_D  ,INTH_D  ,ISTH_D ,IEST_D  ,IWST_D  ,                    &
     &  IMSK_D  ,IMN_D   ,IWT_D  ,                                      &
     &  LEVT_D  ,LEVLST_D,                                              &
     &  PLT_D   ,PLLEN_D ,PLPOS_D,PSLIST_D,NPSLISTS,                    &
     &  BLIM_TS ,TLIM_TS ,BLIMR_TS,TLIMR_TS,IG_TS   ,I1_TS   ,          &
     &  NLIM_TS ,SLIM_TS ,ELIM_TS ,WLIM_TS ,I51_TS  ,NSERIES ,          &
     &  NUPROF  ,LOCN_U  ,IUNT_U ,                                      &
     &  MODEL_ST,ISPACE  ,ITIMA  ,IGP     ,                             &
     &  ILEV    ,IBOT    ,ITOP   ,IFLAG   ,IOPN    ,VMSK    ,           &
     &  IPSEUDO ,IPFIRST ,IPLAST ,PTR_PROG, HALO_TYPE,                  &
     & PPlen2LkUp

! CSTASH end
! MODEL Defines model-dependent quantities used by data addressing and
! STASH
!
! submodel_mod must be used before this one
! VERSION_MOD module is required for nsectp, outfile_s and outfile_e
!
      INTEGER, PARAMETER :: AASSETS    = 9
      INTEGER, PARAMETER :: MEAD_TYPES = 4
      INTEGER, PARAMETER :: A_MAX_TRVARS=150 !Max.no.of tracers allowed
      INTEGER, PARAMETER :: A_MAX_UKCAVARS=150 ! Max.no.of UKCA allowed
      INTEGER, PARAMETER :: MAX_AOBS=100

      REAL :: H_A_EWSPACE
      REAL :: H_A_NSSPACE
      REAL :: H_A_FIRSTLAT
      REAL :: H_A_FIRSTLONG
      REAL :: H_A_POLELAT
      REAL :: H_A_POLELONG

      INTEGER :: H_A_GROUP
      INTEGER :: H_OROG_ROUGH
      INTEGER :: A_ASSMGRPS
      INTEGER :: NUM_PVPR

      LOGICAL :: A_RECON
      LOGICAL :: H_OROG_GRAD
      LOGICAL :: ATMODS
      LOGICAL :: CMODS
      LOGICAL :: LMESO

      LOGICAL :: TRACER_A (0:A_MAX_TRVARS)
      LOGICAL :: TR_UKCA_A (0:A_MAX_UKCAVARS)
      LOGICAL :: AASSET   (AASSETS)
      INTEGER :: AASPF    (AASSETS)
      INTEGER :: AASPL    (AASSETS)
      INTEGER :: AOBINC   (MAX_AOBS)
      INTEGER :: AOBGRP   (MAX_AOBS)
      INTEGER :: RUN_TARGET_END( 6)

      COMMON/MODELA/ H_A_EWSPACE,H_A_NSSPACE,H_A_FIRSTLAT,H_A_FIRSTLONG,&
     &  H_A_POLELAT,H_A_POLELONG,A_ASSMGRPS,NUM_PVPR ,A_RECON,H_A_GROUP,&
     &  H_OROG_GRAD,ATMODS,CMODS,LMESO,TRACER_A,TR_UKCA_A,              &
     &  AASSET,AASPF,AASPL

!Total data length for primary fields for each submodel data partition
      INTEGER      LPRIM(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LPRIM
      INTEGER      global_LPRIM(N_SUBMODEL_PARTITION_MAX)
!Total data length for primary fields for each internal model
      INTEGER      LPrimIM(N_INTERNAL_MODEL_MAX)
!Total data length for diagnostic flds for each submodel data partition
! Global (ie. dump on disk) version of LPrimIM
      INTEGER      global_LPrimIM(N_INTERNAL_MODEL_MAX)
      INTEGER      LDUMP(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LDUMP
      INTEGER      global_LDUMP(N_SUBMODEL_PARTITION_MAX)
!Total data length for diagnostic flds for each internal model
      INTEGER      LDumpIM(N_INTERNAL_MODEL_MAX)
! Global (ie. dump on disk) version of LDumpIM
      INTEGER      global_LDumpIM(N_INTERNAL_MODEL_MAX)
!Total data length for secondary flds for each submodel data partition
      INTEGER      LSECD(N_SUBMODEL_PARTITION_MAX)
!Total data length for secondary flds for each internal model
      INTEGER      LSecdIM(N_INTERNAL_MODEL_MAX)
!Total workspace length for each submodel data partition
      INTEGER      LWORK(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each submodel data partition
      INTEGER      NHeadSub(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each internal model
      INTEGER      NHEAD(N_INTERNAL_MODEL_MAX)
!Total length of extra space for each submod. data part.
      INTEGER      LEXTRA(N_SUBMODEL_PARTITION_MAX)
!Data length for dual-time level ocean fields
      INTEGER      LPRIM_O2
      INTEGER      ITEM_MAX_REQ
      INTEGER      ITEM_MAX_ALL

      INTEGER      NRECS_S
      INTEGER      NTIMES_S
      INTEGER      NSERBLK_S
      INTEGER      NSERREC_S
      INTEGER      NLEVL_S
      INTEGER      NMAXLEV_S
      INTEGER      NPSLISTS_S
      INTEGER      NMAXPSL_S
      INTEGER      NHEAD_FILE(OUTFILE_S:OUTFILE_E)
      LOGICAL      LSTUSER

      COMMON/STRET/                                                     &
     &  LPRIM,LDUMP,LSECD,LWORK,NHEAD,LEXTRA,LPRIM_O2,LPrimIM,LDumpIM,  &
     &  LSecdIM,NHeadSub,ITEM_MAX_REQ,ITEM_MAX_ALL,NSERBLK_S,NSERREC_S, &
     &  NLEVL_S,NMAXLEV_S,NPSLISTS_S,NMAXPSL_S,LSTUSER,NRECS_S,NTIMES_S,&
     &  NHEAD_FILE,                                                     &
     &  global_LPRIM,global_LPrimIM,global_LDUMP,global_LDumpIM
      CHARACTER(LEN=1)  H_ATMOS
      CHARACTER(LEN=1)  H_FLOOR
      CHARACTER(LEN=1)  H_STRAT
      CHARACTER(LEN=1)  H_GLOBAL(N_INTERNAL_MODEL_MAX         )
      INTEGER      H_VERS  (N_INTERNAL_MODEL_MAX,0:NSECTP)

      COMMON/CHOICE/ H_ATMOS,H_GLOBAL,H_FLOOR,H_STRAT

      COMMON/HVERS/ H_VERS

! These are set in SETMODL:
      INTEGER MEAN_NUMBER(N_INTERNAL_MODEL_MAX)
      COMMON/MODLMEAN/ MEAN_NUMBER


! Variables read in by namelist and used in SETMODL
      INTEGER      OCAAA 
      REAL         EWSPACEA,NSSPACEA
      REAL         FRSTLATA,FRSTLONA

      LOGICAL      ZonAvOzone
      LOGICAL      ZonAvTppsOzone
      REAL         LATS
      REAL         LONS
      INTEGER      LWBND
      INTEGER      OCALB
      REAL         POLELATA
      REAL         POLELONA
      INTEGER      SWBND
      INTEGER      TCA(A_MAX_TRVARS)
      INTEGER      TCA_LBC(A_MAX_TRVARS)  ! =1 if tracer in lbc file 
      INTEGER      TC_UKCA(A_MAX_UKCAVARS)
      INTEGER      TC_LBC_UKCA(A_MAX_UKCAVARS) ! =1 if tr in lbc file 
      INTEGER      StLevGWdrag
      INTEGER      BotVDiffLev
      INTEGER      TopVDiffLev


      COMMON/STSHCOMM/                                                  &
     &  RUN_TARGET_END,                                                 &
     &  OCAAA,EWSPACEA,POLELATA,FRSTLATA,LATS,                          &
     &  NSSPACEA,POLELONA,FRSTLONA,LONS,                                &
     &  SWBND,LWBND,                                                    &
     &  ZonAvOzone ,ZonAvTppsOzone, AOBINC,  AOBGRP,                    &
     &  StLevGWdrag, BotVDiffLev,TopVDiffLev,                           &
     &  OCALB,TCA,TCA_LBC,TC_UKCA,TC_LBC_UKCA


      CHARACTER(LEN=1) :: LFLOOR
      CHARACTER(LEN=1) :: OROGR
      CHARACTER(LEN=1) :: SWMCR
      CHARACTER(LEN=1) :: MESO

      COMMON/STSHCHAR/                                                  &
     &     LFLOOR,                                                      &
     &  OROGR,   SWMCR, MESO

      NAMELIST/STSHCOMP/                                                &
        RUN_TARGET_END,                                                 &
        OCAAA       ,EWSPACEA    ,POLELATA ,FRSTLATA  ,LATS   ,         &
                     NSSPACEA    ,POLELONA ,FRSTLONA  ,LONS   ,         &
        SWBND       ,LWBND                            ,OROGR  ,         &
        ZonAvOzone  ,SWMCR       ,MESO     ,                            &
        OCALB       ,LFLOOR      ,AOBINC   ,TCA,                        &
        TCA_LBC     ,TC_UKCA     ,TC_LBC_UKCA   ,AOBGRP          
  
! MODEL end

      INTEGER :: n_max_diags     ! max No of diagnostics
      INTEGER :: idiag           ! item number for diagnostic
      INTEGER :: i, j, k         ! counters
      INTEGER :: ierr            ! error code
      CHARACTER(LEN=72) :: cmessage                ! error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_LOAD_DEFAULT_FLUXES',zhook_in,zhook_handle)

      n_max_diags = nitemp

! Use pre-determined reaction list
      ALLOCATE(asad_chemical_fluxes(n_chemical_fluxes))

      asad_chemical_fluxes=(/        &
         asad_trop_ox_budget_prod,   & ! 20
         asad_trop_ox_budget_loss01, & ! 23 43
         asad_trop_ox_budget_loss02, & ! 12 55
         asad_trop_ox_budget_drydep, & ! 11 66
         asad_trop_ox_budget_wetdep, & ! 7  73
         asad_trop_other_fluxes,     & ! 16 89
         asad_general_interest,      & ! 8  97
         asad_trop_co_budget,        & ! 21 118
         asad_lightning_diags,       & ! 6  124
         asad_strat_oh_prod,         & ! 31 155
         asad_strat_oh_loss,         & ! 26 181
         asad_strat_o3_budget,       & ! 20 201
         asad_strat_o3_misc,         & ! 15 216
         asad_aerosol_chem           & ! 16 232
         /)
       
      IF (printstatus >= PrStatus_Normal)  THEN
        WRITE(6,*) 'STASH Requests for UKCA fluxes from ASAD:'
        WRITE(6,*) '========================================='
        WRITE(6,*) ' '
        WRITE(6,*) 'N=STASH=TYPE==R1======R2========P1========P2'//     &
                   '========P3========P4========No.='
        j = 0
        DO i = 1, SIZE(asad_chemical_fluxes)
          IF (asad_chemical_fluxes(i)%stash_number /= IMDI) THEN
            WRITE(6,*) i, asad_chemical_fluxes(i)%stash_number,         &
                        asad_chemical_fluxes(i)%diag_type,' ',          &
                        asad_chemical_fluxes(i)%reactants(1),           &
                        asad_chemical_fluxes(i)%reactants(2),           &
                        asad_chemical_fluxes(i)%products(1),            &
                        asad_chemical_fluxes(i)%products(2),            &
                        asad_chemical_fluxes(i)%products(3),            &
                        asad_chemical_fluxes(i)%products(4),            &
                        asad_chemical_fluxes(i)%num_species
            j = j + 1
          ENDIF
        ENDDO
        WRITE(6,*) 'n_max_diags: ',n_max_diags, ' used: ',j
      ENDIF

      IF (lhook) CALL dr_hook('ASAD_LOAD_DEFAULT_FLUXES',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ASAD_LOAD_DEFAULT_FLUXES

! ######################################################################
      FUNCTION PRCOUNT(i)
! To sum total number of reactants (set to 2) and products

      USE yomhook,             ONLY: lhook, dr_hook
      USE parkind1,            ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i                      ! Reaction index

      INTEGER, PARAMETER  :: n_reactants = 2        ! Default no of reactants
      INTEGER             :: prcount                ! No of reactants + products
      INTEGER             :: n                      ! counter

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_FLUX_DAT:PRCOUNT',zhook_in,zhook_handle)

      n = n_reactants
      IF (asad_chemical_fluxes(i)%products(1) /= blank0) n = n+1
      IF (asad_chemical_fluxes(i)%products(2) /= blank0) n = n+1
      IF (asad_chemical_fluxes(i)%products(3) /= blank0) n = n+1
      IF (asad_chemical_fluxes(i)%products(4) /= blank0) n = n+1

      prcount = n

      IF (lhook) CALL dr_hook('ASAD_FLUX_DAT:PRCOUNT',zhook_out,zhook_handle)
      RETURN
      END FUNCTION PRCOUNT


      END MODULE asad_flux_dat
