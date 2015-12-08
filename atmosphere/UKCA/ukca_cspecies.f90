! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Module to contain tracer and species numbers
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      MODULE UKCA_CSPECIES

      USE ukca_option_mod, ONLY: jpctr, jpspec
      USE yomhook,         ONLY: lhook, dr_hook
      USE parkind1,        ONLY: jprb, jpim
      USE Control_Max_Sizes
      IMPLICIT NONE
      PRIVATE

      INTEGER :: m           ! Loop counter

      INTEGER, SAVE, PUBLIC :: n_ox    ! tracer numbers
      INTEGER, SAVE, PUBLIC :: n_o3
      INTEGER, SAVE, PUBLIC :: n_o3s
      INTEGER, SAVE, PUBLIC :: n_nox
      INTEGER, SAVE, PUBLIC :: n_no
      INTEGER, SAVE, PUBLIC :: n_no2
      INTEGER, SAVE, PUBLIC :: n_no3
      INTEGER, SAVE, PUBLIC :: n_n2o5
      INTEGER, SAVE, PUBLIC :: n_ho2no2
      INTEGER, SAVE, PUBLIC :: n_hono2
      INTEGER, SAVE, PUBLIC :: n_h2o2
      INTEGER, SAVE, PUBLIC :: n_ch4
      INTEGER, SAVE, PUBLIC :: n_sx
      INTEGER, SAVE, PUBLIC :: n_h2
      INTEGER, SAVE, PUBLIC :: n_h2o
      INTEGER, SAVE, PUBLIC :: n_cl
      INTEGER, SAVE, PUBLIC :: n_ox1
      INTEGER, SAVE, PUBLIC :: n_sx1
      INTEGER, SAVE, PUBLIC :: n_dms
      INTEGER, SAVE, PUBLIC :: n_so2
      INTEGER, SAVE, PUBLIC :: n_so3
      INTEGER, SAVE, PUBLIC :: n_h2so4
      INTEGER, SAVE, PUBLIC :: n_sec_org
      INTEGER, SAVE, PUBLIC :: n_cfcl3
      INTEGER, SAVE, PUBLIC :: n_cf2cl2
      INTEGER, SAVE, PUBLIC :: n_bro
      INTEGER, SAVE, PUBLIC :: n_hcl
      INTEGER, SAVE, PUBLIC :: n_o1d
      INTEGER, SAVE, PUBLIC :: n_age ! age of air
      INTEGER, SAVE, PUBLIC :: n_passive ! passive o3


      INTEGER, SAVE, PUBLIC :: nn_o3    ! Species numbers
      INTEGER, SAVE, PUBLIC :: nn_o3s
      INTEGER, SAVE, PUBLIC :: nn_oh
      INTEGER, SAVE, PUBLIC :: nn_ho2
      INTEGER, SAVE, PUBLIC :: nn_h2o2
      INTEGER, SAVE, PUBLIC :: nn_no
      INTEGER, SAVE, PUBLIC :: nn_no2
      INTEGER, SAVE, PUBLIC :: nn_o1d
      INTEGER, SAVE, PUBLIC :: nn_o3p
      INTEGER, SAVE, PUBLIC :: nn_meoo
      INTEGER, SAVE, PUBLIC :: nn_meco3
      INTEGER, SAVE, PUBLIC :: nn_etoo
      INTEGER, SAVE, PUBLIC :: nn_etco3
      INTEGER, SAVE, PUBLIC :: nn_nproo
      INTEGER, SAVE, PUBLIC :: nn_nprooh
      INTEGER, SAVE, PUBLIC :: nn_iProo
      INTEGER, SAVE, PUBLIC :: nn_iProoh
      INTEGER, SAVE, PUBLIC :: nn_mecoch2oo
      INTEGER, SAVE, PUBLIC :: nn_no3
      INTEGER, SAVE, PUBLIC :: nn_n2o5
      INTEGER, SAVE, PUBLIC :: nn_ho2no2
      INTEGER, SAVE, PUBLIC :: nn_hono2
      INTEGER, SAVE, PUBLIC :: nn_ch4
      INTEGER, SAVE, PUBLIC :: nn_so2
      INTEGER, SAVE, PUBLIC :: nn_so3
      INTEGER, SAVE, PUBLIC :: nn_h2so4
      INTEGER, SAVE, PUBLIC :: nn_ohs 
      INTEGER, SAVE, PUBLIC :: nn_ho2s 
      INTEGER, SAVE, PUBLIC :: nn_o1ds 
      INTEGER, SAVE, PUBLIC :: nn_o3ps 
      INTEGER, SAVE, PUBLIC :: nn_n
      INTEGER, SAVE, PUBLIC :: nn_cl
      INTEGER, SAVE, PUBLIC :: nn_clo
      INTEGER, SAVE, PUBLIC :: nn_hcl
      INTEGER, SAVE, PUBLIC :: nn_cl2o2
      INTEGER, SAVE, PUBLIC :: nn_meo
      INTEGER, SAVE, PUBLIC :: nn_bro
      INTEGER, SAVE, PUBLIC :: nn_br
      INTEGER, SAVE, PUBLIC :: nn_buoo        !for RAQ chemistry
      INTEGER, SAVE, PUBLIC :: nn_meko2
      INTEGER, SAVE, PUBLIC :: nn_hoc2h4o2
      INTEGER, SAVE, PUBLIC :: nn_hoc3h6o2
      INTEGER, SAVE, PUBLIC :: nn_oxyl1
      INTEGER, SAVE, PUBLIC :: nn_memald1
      INTEGER, SAVE, PUBLIC :: nn_hoipo2
      INTEGER, SAVE, PUBLIC :: nn_homvko2
      INTEGER, SAVE, PUBLIC :: nn_tolp1
      INTEGER, SAVE, PUBLIC :: nn_cfcl3
      INTEGER, SAVE, PUBLIC :: nn_cf2cl2

! Advected tracers
      REAL, SAVE, ALLOCATABLE, PUBLIC :: c_species(:)
! Non-advected tracers
      REAL, SAVE, ALLOCATABLE, PUBLIC :: c_na_species(:) 

      PUBLIC UKCA_CALC_CSPECIES

      CONTAINS

      SUBROUTINE UKCA_CALC_CSPECIES

      USE UKCA_D1_DEFS,    ONLY: nukca_D1items, UkcaD1Codes, UKCA_sect, &
                                 nm_spec
      USE UKCA_CONSTANTS
      USE ASAD_MOD,        ONLY: advt, nadvt, speci
      USE ereport_mod,     ONLY: ereport
      USE Control_Max_Sizes
      USE PrintStatus_mod
      IMPLICIT NONE

      INTEGER :: errcode                ! Variable passed to ereport
      INTEGER :: i
      CHARACTER(LEN=72) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     Compute list of conversion factors from mmr to vmr.
!     Add new entry if new tracers are introduced, and add value for c_xx
!      in UKCA_CONSTANTS, where xx is the species.

      IF (lhook) CALL dr_hook('UKCA_CSPECIES:UKCA_CALC_CSPECIES',zhook_in,  &
                              zhook_handle)
! advected tracers
      c_species = 0.
      WHERE (advt == 'Ox        ') c_species = c_o3
      WHERE (advt == 'O3        ') c_species = c_o3
      WHERE (advt == 'NO        ') c_species = c_no
      WHERE (advt == 'NO2       ') c_species = c_no2
      WHERE (advt == 'NO3       ') c_species = c_no3
      WHERE (advt == 'NOx       ') c_species = c_no2
      WHERE (advt == 'N2O5      ') c_species = c_n2o5
      WHERE (advt == 'HO2NO2    ') c_species = c_ho2no2
      WHERE (advt == 'HONO2     ') c_species = c_hono2
      WHERE (advt == 'H2O2      ') c_species = c_h2o2
      WHERE (advt == 'CH4       ') c_species = c_ch4
      WHERE (advt == 'CO        ') c_species = c_co
      WHERE (advt == 'HCHO      ') c_species = c_hcho
      WHERE (advt == 'MeOOH     ') c_species = c_meooh
      WHERE (advt == 'HONO      ') c_species = c_hono
      WHERE (advt == 'C2H6      ') c_species = c_c2h6
      WHERE (advt == 'EtOOH     ') c_species = c_etooh
      WHERE (advt == 'MeCHO     ') c_species = c_mecho
      WHERE (advt == 'PAN       ') c_species = c_pan
      WHERE (advt == 'C3H8      ') c_species = c_c3h8
      WHERE (advt == 'i-PrOOH   ') c_species = c_prooh
      WHERE (advt == 'n-PrOOH   ') c_species = c_prooh
      WHERE (advt == 'EtCHO     ') c_species = c_etcho
      WHERE (advt == 'Me2CO     ') c_species = c_me2co
      WHERE (advt == 'MeCOCH2OOH') c_species = c_mecoch2ooh
      WHERE (advt == 'PPAN      ') c_species = c_ppan
      WHERE (advt == 'MeONO2    ') c_species = c_meono2
      WHERE (advt == 'Sx        ') c_species = c_o3
      WHERE (advt == 'O3S       ') c_species = c_o3
      WHERE (advt == 'HOx       ') c_species = c_ho2
      WHERE (advt == 'N2O       ') c_species = c_n2o
      WHERE (advt == 'CFCl3     ') c_species = c_cfcl3
      WHERE (advt == 'CF2Cl2    ') c_species = c_cf2cl2
      WHERE (advt == 'H2O       ') c_species = c_h2o
      WHERE (advt == 'ClONO2    ') c_species = c_clono2
      WHERE (advt == 'Clx       ') c_species = c_clo
      WHERE (advt == 'ClO       ') c_species = c_clo
      WHERE (advt == 'Cl2O2     ') c_species = c_cl2o2
      WHERE (advt == 'Cl        ') c_species = c_cl
      WHERE (advt == 'HCl       ') c_species = c_hcl
      WHERE (advt == 'HOCl      ') c_species = c_hocl
      WHERE (advt == 'OClO      ') c_species = c_oclo
      WHERE (advt == 'Br        ') c_species = c_br
      WHERE (advt == 'BrO       ') c_species = c_bro
      WHERE (advt == 'Brx       ') c_species = c_bro
      WHERE (advt == 'HOBr      ') c_species = c_hobr
      WHERE (advt == 'BrONO2    ') c_species = c_brono2
      WHERE (advt == 'BrCl      ') c_species = c_brcl
      WHERE (advt == 'MeBr      ') c_species = c_mebr
      WHERE (advt == 'HBr       ') c_species = c_hbr
      WHERE (advt == 'CF2ClCFCl2') c_species = c_cf2clcfcl2
      WHERE (advt == 'CHF2Cl    ') c_species = c_chf2cl
      WHERE (advt == 'MeCCl3    ') c_species = c_meccl3
      WHERE (advt == 'CCl4      ') c_species = c_ccl4
      WHERE (advt == 'MeCl      ') c_species = c_mecl
      WHERE (advt == 'CF2ClBr   ') c_species = c_cf2clbr
      WHERE (advt == 'CF3Br     ') c_species = c_cf3br
      WHERE (advt == 'CH2Br2    ') c_species = c_ch2br2
      WHERE (advt == 'C5H8      ') c_species = c_c5h8
      WHERE (advt == 'ISO2      ') c_species = c_iso2
      WHERE (advt == 'ISOOH     ') c_species = c_isooh
      WHERE (advt == 'ISON      ') c_species = c_ison
      WHERE (advt == 'MACR      ') c_species = c_macr
      WHERE (advt == 'MACRO2    ') c_species = c_macro2
      WHERE (advt == 'MACROOH   ') c_species = c_macrooh
      WHERE (advt == 'MPAN      ') c_species = c_mpan
      WHERE (advt == 'HACET     ') c_species = c_hacet
      WHERE (advt == 'MGLY      ') c_species = c_mgly
      WHERE (advt == 'NALD      ') c_species = c_nald
      WHERE (advt == 'HCOOH     ') c_species = c_hcooh
      WHERE (advt == 'MeCO3H    ') c_species = c_meco3h
      WHERE (advt == 'MeCO2H    ') c_species = c_meco2h
      WHERE (advt == 'SO2       ') c_species = c_so2
      WHERE (advt == 'SO3       ') c_species = c_so3
      WHERE (advt == 'DMS       ') c_species = c_dms
      WHERE (advt == 'DMSO      ') c_species = c_dmso
      WHERE (advt == 'Me2S      ') c_species = c_me2s
      WHERE (advt == 'COS       ') c_species = c_cos
      WHERE (advt == 'H2S       ') c_species = c_h2s
      WHERE (advt == 'CS2       ') c_species = c_cs2
      WHERE (advt == 'MSA       ') c_species = c_msa
      WHERE (advt == 'H2SO4     ') c_species = c_h2so4
      WHERE (advt == 'NH3       ') c_species = c_nh3
      WHERE (advt == 'MeOH      ') c_species = c_meoh
      WHERE (advt == 'Monoterp  ') c_species = c_monoterp
      WHERE (advt == 'Sec_Org   ') c_species = c_sec_org
      WHERE (advt == 'O(3P)     ') c_species = c_o3p 
      WHERE (advt == 'OH        ') c_species = c_oh 
      WHERE (advt == 'HO2       ') c_species = c_ho2 
      WHERE (advt == 'MeOO      ') c_species = c_meoo 
      WHERE (advt == 'EtOO      ') c_species = c_etoo 
      WHERE (advt == 'MeCO3     ') c_species = c_meco3 
      WHERE (advt == 'n-PrOO    ') c_species = c_proo 
      WHERE (advt == 'i-PrOO    ') c_species = c_proo 
      WHERE (advt == 'EtCO3     ') c_species = c_etco3 
      WHERE (advt == 'MeCOCH2OO ') c_species = c_mecoch2oo 
      WHERE (advt == 'O(3P)S    ') c_species = c_o3p 
      WHERE (advt == 'O(1D)S    ') c_species = c_o1d 
      WHERE (advt == 'OHS       ') c_species = c_oh 
      WHERE (advt == 'HO2S      ') c_species = c_ho2 
      WHERE (advt == 'N         ') c_species = c_n
      WHERE (advt == 'H         ') c_species = c_h
      WHERE (advt == 'H2        ') c_species = c_h2       !for RAQ chemistry
      WHERE (advt == 'C4H10     ') c_species = c_c4h10
      WHERE (advt == 'MEK       ') c_species = c_mek
      WHERE (advt == 'C2H4      ') c_species = c_c2h4
      WHERE (advt == 'C3H6      ') c_species = c_c3h6
      WHERE (advt == 'oXYLENE   ') c_species = c_oxylene
      WHERE (advt == 's-BuOOH   ') c_species = c_buooh
      WHERE (advt == 'CH3OH     ') c_species = c_ch3oh
      WHERE (advt == 'MeOH      ') c_species = c_ch3oh
      WHERE (advt == 'GLY       ') c_species = c_gly
      WHERE (advt == 'MEMALD    ') c_species = c_memald
      WHERE (advt == 'MVK       ') c_species = c_mvk
      WHERE (advt == 'MVKOOH    ') c_species = c_mvkooh
      WHERE (advt == 'TOLUENE   ') c_species = c_toluene
      WHERE (advt == 'RNC2H4    ') c_species = c_rnc2h4
      WHERE (advt == 'RNC3H6    ') c_species = c_rnc3h6
      WHERE (advt == 'ORGNIT    ') c_species = c_orgnit
      WHERE (advt == 'PASSIVE O3') c_species = 1.0
      WHERE (advt == 'AGE OF AIR') c_species = 1.0

! non-advected tracers 
      c_na_species=0.0 
      WHERE (nadvt == 'OH        ') c_na_species = c_oh 
      WHERE (nadvt == 'HO2       ') c_na_species = c_ho2 
      WHERE (nadvt == 'O(3P)     ') c_na_species = c_o3p 
      WHERE (nadvt == 'O3P       ') c_na_species = c_o3p 
      WHERE (nadvt == 'O(1D)     ') c_na_species = c_o1d 
      WHERE (nadvt == 'O1D       ') c_na_species = c_o1d 
      WHERE (nadvt == 'MeOO      ') c_na_species = c_meoo 
      WHERE (nadvt == 'EtOO      ') c_na_species = c_etoo 
      WHERE (nadvt == 'MeCO3     ') c_na_species = c_meco3 
      WHERE (nadvt == 'n-PrOO    ') c_na_species = c_proo 
      WHERE (nadvt == 'i-PrOO    ') c_na_species = c_proo 
      WHERE (nadvt == 'EtCO3     ') c_na_species = c_etco3 
      WHERE (nadvt == 'MeCOCH2OO ') c_na_species = c_mecoch2oo 
      WHERE (nadvt == 'O(3P)S    ') c_na_species = c_o3p 
      WHERE (nadvt == 'O(1D)S    ') c_na_species = c_o1d 
      WHERE (nadvt == 'OHS       ') c_na_species = c_oh 
      WHERE (nadvt == 'HO2S      ') c_na_species = c_ho2 
      WHERE (nadvt == 's-BuOO    ') c_na_species = c_buoo     !for RAQ chemistry
      WHERE (nadvt == 'MEKO2     ') c_na_species = c_meko2 
      WHERE (nadvt == 'HOC2H4O2  ') c_na_species = c_hoc2h4o2
      WHERE (nadvt == 'HOC3H6O2  ') c_na_species = c_hoc3h6o2
      WHERE (nadvt == 'OXYL1     ') c_na_species = c_oxyl1
      WHERE (nadvt == 'MEMALD1   ') c_na_species = c_memald1
      WHERE (nadvt == 'HOIPO2    ') c_na_species = c_hoipo2 
      WHERE (nadvt == 'HOMVKO2   ') c_na_species = c_homvko2
      WHERE (nadvt == 'TOLP1     ') c_na_species = c_tolp1 
      
!     Initialise tracer numbers

      n_h2o    = 0
      n_cl     = 0
      n_sx     = 0
      n_n2o5   = 0
      n_ox     = 0
      n_o3     = 0
      n_o3s    = 0
      n_nox    = 0
      n_no     = 0
      n_no2    = 0
      n_no3    = 0
      n_ho2no2 = 0
      n_hono2  = 0
      n_h2o2   = 0
      n_ch4    = 0
      n_h2     = 0
      n_dms    = 0
      n_so2    = 0
      n_so3    = 0
      n_h2so4  = 0
      n_sec_org= 0
      n_age    = 0
      n_passive = 0

!     Find tracer numbers

       DO m=1,jpctr
        SELECT CASE (advt(m))
          CASE ('Ox        ')
            n_ox     = m
            n_ox1    = m
          CASE ('O3        ')
            n_o3     = m
          CASE ('O3S       ')
            n_o3s    = m
          CASE ('Sx        ')   ! stratospheric ozone tracer
            n_sx     = m
            n_sx1    = m
          CASE ('NOx       ')
            n_nox    = m
          CASE ('NO        ')
            n_no     = m
          CASE ('NO2       ')
            n_no2    = m
          CASE ('NO3       ')
            n_no3    = m
          CASE ('N2O5      ')
            n_n2o5   = m
          CASE ('HO2NO2    ')
            n_ho2no2 = m
          CASE ('HONO2     ')
            n_hono2  = m
          CASE ('H2O2      ')
            n_h2o2   = m
          CASE ('CH4       ')
            n_ch4    = m
          CASE ('H2O       ')
            n_h2o    = m
          CASE ('Cl        ')
            n_cl     = m
          CASE ('H2        ')
            n_h2     = m
          CASE ('DMS       ')
            n_dms    = m
          CASE ('SO2       ')
            n_so2    = m
          CASE ('SO3       ')
            n_so3    = m
          CASE ('H2SO4     ')
            n_h2so4  = m
          CASE ('Sec_Org   ')
            n_sec_org= m
          CASE ('BrO       ')
            n_bro    = m
          CASE ('HCl       ')
            n_hcl    = m
          CASE ('CFCl3     ')
            n_cfcl3  = m
          CASE ('CF2Cl2    ')
            n_cf2cl2 = m
          CASE ('O(1D)     ')
            n_o1d    = m
        END SELECT
        END DO
        DO i=1,jpctr
          IF ((advt(i) == 'Ox        ') .OR. (advt(i) == 'O3        ')) &
            n_o3 = i
        END DO

      IF (n_o3 == 0) THEN
        cmessage='Ozone not found among chemical tracers.'
        errcode=1
        CALL EREPORT('UKCA_CSPECIES',errcode,cmessage)
      ENDIF

! Age of air & passive o3 don't appear :
      age_passive_loop: DO i=1,Nukca_D1items
        IF ((UkcaD1Codes(I)%item == 149 .AND.                           &
             UkcaD1Codes(I)%section == UKCA_sect) .AND.                 & 
            nm_spec(149) == 'PASSIVE O3') THEN
           n_passive=i
           IF (n_age > 0) EXIT age_passive_loop
        END IF
        IF ((UkcaD1Codes(I)%item == 150 .AND.                           &
             UkcaD1Codes(I)%section == UKCA_sect) .AND.                 & 
            nm_spec(150) == 'AGE OF AIR') THEN
           n_age=i
           IF (n_passive > 0) EXIT age_passive_loop
        END IF
     ENDDO age_passive_loop

!     Initialise species numbers

      nn_o3        = 0
      nn_o3s       = 0
      nn_oh        = 0
      nn_ho2       = 0
      nn_h2o2      = 0
      nn_no        = 0
      nn_no2       = 0
      nn_o1d       = 0
      nn_o3p       = 0
      nn_meoo      = 0
      nn_meco3     = 0
      nn_etoo      = 0
      nn_etco3     = 0
      nn_nproo     = 0
      nn_nprooh    = 0
      nn_iproo     = 0
      nn_iprooh    = 0
      nn_mecoch2oo = 0
      nn_so2       = 0
      nn_so3       = 0
      nn_h2so4     = 0
      nn_ohs       = 0 
      nn_ho2s      = 0 
      nn_o1ds      = 0 
      nn_o3ps      = 0 
      nn_n         = 0
      nn_meo       = 0
      nn_cl        = 0
      nn_clo       = 0
      nn_cl2o2     = 0
      nn_hcl       = 0
      nn_br        = 0
      nn_bro       = 0
      nn_buoo      = 0   ! species for RAQ chemistry
      nn_meko2     = 0
      nn_hoc2h4o2  = 0
      nn_hoc3h6o2  = 0
      nn_oxyl1     = 0
      nn_memald1   = 0
      nn_hoipo2    = 0
      nn_homvko2   = 0
      nn_tolp1     = 0

!     Find species numbers

      DO m=1,jpspec
        SELECT CASE(speci(m))
          CASE ('O3        ')
            nn_o3 = m
          CASE ('O3S       ')
            nn_o3s = m
          CASE ('OH        ')
            nn_oh = m
          CASE ('HO2       ')
            nn_ho2 = m
          CASE ('H2O2      ')
            nn_h2o2 = m
          CASE ('NO        ')
            nn_no = m
          CASE ('NO2       ')
            nn_no2 = m
          CASE ('O(1D)     ')
            nn_o1d = m
          CASE ('O(3P)     ')
            nn_o3p = m
          CASE ('MeOO      ')
            nn_meoo = m
          CASE ('MeCO3     ')
            nn_meco3 = m
          CASE ('EtOO      ')
            nn_etoo = m
          CASE ('EtCO3     ')
            nn_etco3 = m
          CASE ('n-PrOO    ')
            nn_nProo = m
          CASE ('n-PrOOH   ')
            nn_nprooh = m
          CASE ('i-PrOO    ')
            nn_iProo = m
          CASE ('i-PrOOH   ')
            nn_iprooh = m
          CASE ('MeCOCH2OO ')
            nn_mecoch2oo = m
          CASE ('NO3       ')
            nn_no3 = m
          CASE ('N2O5      ')
            nn_n2o5 = m
          CASE ('HO2NO2    ')
            nn_ho2no2 = m
          CASE ('HONO2     ')
            nn_hono2 = m
          CASE ('CH4       ')
            nn_ch4 = m
          CASE ('SO2       ')
            nn_so2 = m
          CASE ('SO3       ')
            nn_so3 = m
          CASE ('H2SO4     ')
            nn_h2so4 = m
          CASE ('O(1D)S    ') 
            nn_o1ds = m 
          CASE ('O(3P)S    ') 
            nn_o3ps = m 
          CASE ('OHS       ') 
            nn_ohs = m 
          CASE ('HO2S      ') 
            nn_ho2s = m 
          CASE ('N         ')
            nn_n = m
          CASE ('MeO       ')
            nn_meo = m
          CASE ('Cl        ')
            nn_cl = m
          CASE ('ClO       ')
            nn_clo = m
          CASE ('Cl2O2     ')
            nn_cl2o2 = m
          CASE ('HCl       ')
            nn_hcl = m
          CASE ('Br        ')
            nn_br = m
          CASE ('BrO       ')
            nn_bro = m
          CASE ('s-BuOO    ')
            nn_buoo = m
          CASE ('MEKO2     ')
            nn_meko2 = m
          CASE ('HOC2H4O2  ')
            nn_hoc2h4o2 = m
          CASE ('HOC3H6O2  ')
            nn_hoc3h6o2 = m
          CASE ('OXYL1     ')
            nn_oxyl1 = m
          CASE ('MEMALD1   ')
            nn_memald1 = m
          CASE ('HOIPO2    ')
            nn_hoipo2 = m
          CASE ('HOMVKO2   ')
            nn_homvko2 = m
          CASE ('TOLP1     ')
            nn_tolp1 = m
          CASE ('CFCl3     ')
            nn_cfcl3 = m
          CASE ('CF2Cl2    ')
            nn_cf2cl2 = m
        END SELECT
      ENDDO
      IF (printstatus > prstatus_oper) THEN
        WRITE(6,'(A)') 'Species indices:'
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'O3 = ',nn_o3,' OH = ',nn_oh
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'HO2 = ',nn_ho2,' NO = ',nn_no
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'NO2 = ',nn_no2,' O(1D) = ',nn_o1d
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'MeOO = ',nn_meoo,' MeCO3 = ',nn_meco3
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'EtOO = ',nn_etoo,' EtCO3 = ',nn_etco3
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'n-Proo = ',nn_nproo,'n-Prooh = ',nn_nprooh
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'i-Proo = ',nn_iproo,'i-Prooh = ',nn_iprooh
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'MeCOCH2OO = ',nn_mecoch2oo
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'SO2: ',nn_so2,' H2O2: ',nn_h2o2
        WRITE(6,'(A,1X,I6,A,1X,I6)') 'SO3: ',nn_so3,' H2SO4: ',nn_h2so4
      END IF

      IF (lhook) CALL dr_hook('UKCA_CSPECIES:UKCA_CALC_CSPECIES',zhook_out, &
                              zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALC_CSPECIES

      END MODULE UKCA_CSPECIES
