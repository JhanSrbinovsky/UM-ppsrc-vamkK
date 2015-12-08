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
!    Module to define setup and indices for gas and aerosol tracers.
!    Many of these are not used in UM but kept for consistency with
!    set-up in TOMCAT.
!    Contains public subroutines:
!      UKCA_INDICES_ORGV1_SOto3
!      UKCA_INDICES_ORGV1_SOto3_COUPL
!      UKCA_INDICES_ORGV1_SOto6
!      UKCA_INDICES_ORGV1_SOto6_COUPL
!      UKCA_INDICES_SV1
!      UKCA_INDICES_SV1_COUPLED
!      UKCA_INDICES_NOCHEM
!      UKCA_INDICES_TRAQU38
!      UKCA_INDICES_TRAQU9
!      UKCA_INDICES_SUSS_4MODE
!      UKCA_INDICES_SUSSBCOC_5MODE
!      UKCA_INDICES_SUSSBCOC_4MODE
!      UKCA_INDICES_SUSSBCOCSO_5MODE
!      UKCA_INDICES_SUSSBCOCSO_4MODE
!      UKCA_INDICES_DUonly_2MODE
!      UKCA_INDICES_DUonly_3MODE (needs to be added at some point)
!      UKCA_INDICES_SUSSBCOCDU_7MODE
!      UKCA_INDICES_SUSSBCOCDU_4MODE
!    which define particular setups.
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
      MODULE UKCA_SETUP_INDICES
! ---------------------------------------------------------------------|
!+ Module to define setup and indices for gas and aerosol tracers.
!+ Many of these are not used in UM but kept for consistency with
!+ set-up in TOMCAT.



      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER :: NTRAER     ! # of aerosol advected tracers
      INTEGER :: NCHEMG     ! # of gas phase chemistry tracers
      INTEGER :: ICHEM      ! 1/0 = gas phase chemistry tracers on/off
      INTEGER :: NADVG      ! # of gas phase advected tracers
      INTEGER :: NOFFOX     ! # of offline oxidant species
      INTEGER :: NTRAG      ! total # of gas phase species
      INTEGER :: BUDGET     ! 1/0 = budget fields on/off
      INTEGER :: NBUDCHEM   ! # of gas chem budget fields
      INTEGER :: IDUSTDEP   ! 1/0 = size-resolved dust dep. fields on/off
      INTEGER :: NDUSTDEP   ! # of size-resolved dust deposition fields
      INTEGER :: NBUDAER    ! # of aerosol budget fields
      INTEGER :: NBUDAERTOT ! total # of aersol budget fields
      INTEGER :: NBUDGET    ! total # of budget terms
      INTEGER :: TRAQU      ! 1/0 = extra diagnostic terms on/off
      INTEGER :: NTRAQU     ! # of extra diagnostic terms
      INTEGER :: GASBUDGET  ! 1/0 = gas phase chemistry fluxes on/off
      INTEGER :: NGASBUDGET ! # of gas phase chemistry fluxes
!
! .. (indices for 40 gas phase tropospheric chemistry species in S0)
      INTEGER :: MOX
      INTEGER :: MNOX
      INTEGER :: MN2O5
      INTEGER :: MHNO4
      INTEGER :: MHNO3
      INTEGER :: MH2O2
      INTEGER :: MCH4
      INTEGER :: MCO
      INTEGER :: MCH2O
      INTEGER :: MMHP
      INTEGER :: MHONO
      INTEGER :: MC2H6
      INTEGER :: METOOH
      INTEGER :: MMECHO
      INTEGER :: MPAN
      INTEGER :: MC3H8
      INTEGER :: MPNOOH
      INTEGER :: MPIOOH
      INTEGER :: METCHO
      INTEGER :: MME2CO
      INTEGER :: MMECOH
      INTEGER :: MPPAN
      INTEGER :: MMENO3
      INTEGER :: MOXS
      INTEGER :: MNOYS
      INTEGER :: MISOP
      INTEGER :: MC2H4
      INTEGER :: MC2H2
      INTEGER :: MISOOH
      INTEGER :: MISON
      INTEGER :: MMACR
      INTEGER :: MMACROOH
      INTEGER :: MMPAN
      INTEGER :: MHACET
      INTEGER :: MMGLY
      INTEGER :: MNALD
      INTEGER :: MHCOOH
      INTEGER :: MMECO3H
      INTEGER :: MMeCO2H
      INTEGER :: MMEOH
!
! .. (indices for advected gas phase tracers in array S0)
! .. (8 sulfur species, H2O2, MONOTER, SEC_ORG, Q3D, PT)
      INTEGER :: MSOTWO     ! for SO2
      INTEGER :: MMeSMe     ! for DMS
      INTEGER :: MH2SO4     ! for H2SO4
      INTEGER :: MDMSO      ! for DMSO
      INTEGER :: MMSA       ! for MSA
      INTEGER :: MCS2       ! for CS2
      INTEGER :: MH2S       ! for H2S
      INTEGER :: MCOS       ! for COS
      INTEGER :: MMONOTER   ! for lumped monoterpene species
      INTEGER :: MSEC_ORG   ! for involatile organic species
      INTEGER :: MH2O2F     ! for H2O2 semi-prognostic
      INTEGER :: MQ3D       ! for water vapour
      INTEGER :: MPT        ! for potential temperature
!
! .. (indices for tropospheric chemistry species in array ST)
      INTEGER :: NO
      INTEGER :: NO1D
      INTEGER :: NO3
      INTEGER :: NNO
      INTEGER :: NNO3
      INTEGER :: NNO2
      INTEGER :: NN2O5
      INTEGER :: NHNO4
      INTEGER :: NHNO3
      INTEGER :: NOH
      INTEGER :: NHO2
      INTEGER :: NH2O2
      INTEGER :: NCH4
      INTEGER :: NCO
      INTEGER :: NCH2O
      INTEGER :: NMEOO
      INTEGER :: NH2O
      INTEGER :: NMHP
      INTEGER :: NHONO
      INTEGER :: NC2H6
      INTEGER :: NETOO
      INTEGER :: NETOOH
      INTEGER :: NMECHO
      INTEGER :: NMECO3
      INTEGER :: NPAN
      INTEGER :: NC3H8
      INTEGER :: NPNOO
      INTEGER :: NPIOO
      INTEGER :: NPNOOH
      INTEGER :: NPIOOH
      INTEGER :: NETCHO
      INTEGER :: NETCO3
      INTEGER :: NME2CO
      INTEGER :: NMECOO
      INTEGER :: NMECOH
      INTEGER :: NPPAN
      INTEGER :: NMENO3
      INTEGER :: NOS
      INTEGER :: NO1DS
      INTEGER :: NO3S
      INTEGER :: NNOXS
      INTEGER :: NHNO3S
      INTEGER :: NNOYS
      INTEGER :: NISOP
      INTEGER :: NC2H4
      INTEGER :: NC2H2
      INTEGER :: NISO2
      INTEGER :: NISOOH
      INTEGER :: NISON
      INTEGER :: NMACR
      INTEGER :: NMACRO2
      INTEGER :: NMACROOH
      INTEGER :: NMPAN
      INTEGER :: NHACET
      INTEGER :: NMGLY
      INTEGER :: NNALD
      INTEGER :: NHCOOH
      INTEGER :: NMECO3H
      INTEGER :: NMECO2H
      INTEGER :: NMEOH
!
! .. (indices for all gas phase tracers in array ST)
      INTEGER :: NSOTWO     ! for SO2
      INTEGER :: NMeSMe     ! for DMS
      INTEGER :: NH2SO4     ! for H2SO4
      INTEGER :: NDMSO      ! for DMSO
      INTEGER :: NMSA       ! for MSA
      INTEGER :: NCS2       ! for CS2
      INTEGER :: NH2S       ! for H2S
      INTEGER :: NCOS       ! for COS
      INTEGER :: NMONOTER   ! for lumped monoterpene species
      INTEGER :: NSEC_ORG   ! for involatile organic species
      INTEGER :: NH2O2F     ! for H2O2 (semi-prognostic)
      INTEGER :: NO3F       ! for O3 (offline)
      INTEGER :: NOHF       ! for OH (offline)
      INTEGER :: NNO3F      ! for NO3 (offline)
      INTEGER :: NQ3D       ! for water vapour
      INTEGER :: NPT        ! for potential temperature
!
! .. (indices for gas phase budget mass fluxes)
      INTEGER :: NDMSEMOC   ! for DMS emissions flux (oceanic sources)
      INTEGER :: NDMSTEND   ! for DMS ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NSO2EMAN   ! for SO2 emissions flux (anthropogenic sources)
      INTEGER :: NSO2EMBM   ! for SO2 emissions flux (biomass burning sources)
      INTEGER :: NSO2EMVL   ! for SO2 emissions flux (volcanic sources)
      INTEGER :: NSO2TEND   ! for SO2 ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NSO2DDEP   ! for SO2 dry deposition flux
      INTEGER :: NSO2WDEP   ! for SO2 wet deposition flux
      INTEGER :: NH2SO4TEND ! for H2SO4 ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NH2SO4DDEP ! for H2SO4 dry deposition flux
      INTEGER :: NCOSEMAN   ! for COS emissions flux (anthropogenic sources)
      INTEGER :: NCOSEMOC   ! for COS emissions flux (oceanic sources)
      INTEGER :: NCOSTEND   ! for COS ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NCS2EMAN   ! for CS2 emissions flux (anthropogenic sources)
      INTEGER :: NCS2EMOC   ! for CS2 emissions flux (oceanic sources)
      INTEGER :: NCS2TEND   ! for CS2 ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NDMSOTEND  ! for DMSO ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NDMSODDEP  ! for DMSO dry deposition flux
      INTEGER :: NMSATEND   ! for MSA ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NMSADDEP   ! for MSA dry deposition flux
      INTEGER :: NTERP_EM   ! for terpene emissions flux (biogenic sources)
      INTEGER :: NTERP_TEND ! for terpene ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NTERP_DDEP ! for terpene dry deposition flux
      INTEGER :: NSORG_TEND ! for sec_org ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NSORG_DDEP ! for sec_org dry deposition flux
      INTEGER :: NSORG_WDEP ! for sec_org wet deposition flux
!
! .. (indices for gas phase chemistry fluxes)
      INTEGER :: IOHDMS1  ! DMS +OH -->   SO2(+   MeOO +HCHO)
      INTEGER :: IOHDMS2  ! DMS +OH -->0.6SO2 +0.4DMSO(+MeOO)
      INTEGER :: INO3DMS  ! DMS +NO3-->   SO2(+  HONO2 +MeOO+HCHO)
      INTEGER :: IDMSOOH1 ! DMSO+OH -->   SO2             
      INTEGER :: IDMSOOH2 ! DMSO+OH -->0.6SO2 +0.4MSA
      INTEGER :: ICS2OH   ! DMS +OH -->   SO2 +COS         
      INTEGER :: IH2SOH   ! DMS +OH -->   SO2
      INTEGER :: ICOSOH   ! DMS +OH -->   SO2
! 
! .. (indices for aerosol budget mass fluxes)
      INTEGER :: NMASPRIMSUAITSOL ! for primary SU ems (all srces) to Aitsol
      INTEGER :: NMASPRIMSUACCSOL ! for primary SU ems (all srces) to accsol
      INTEGER :: NMASPRIMSUCORSOL ! for primary SU ems (all srces) to corsol
      INTEGER :: NMASPRIMSSACCSOL ! for primary SS ems (ocean) to accsol
      INTEGER :: NMASPRIMSSCORSOL ! for primary SS ems (ocean) to corsol
      INTEGER :: NMASPRIMBCAITSOL ! for primary BC ems (all srces) to Aitsol
      INTEGER :: NMASPRIMBCAITINS ! for primary BC ems (all srces) to Aitins
      INTEGER :: NMASPRIMOCAITSOL ! for primary OC ems (all srces) to Aitsol
      INTEGER :: NMASPRIMOCAITINS ! for primary OC ems (all srces) to Aitins
      INTEGER :: NMASDDEPSUNUCSOL ! for drydep of SU from nucsol
      INTEGER :: NMASDDEPSUAITSOL ! for drydep of SU from Aitsol
      INTEGER :: NMASDDEPSUACCSOL ! for drydep of SU from accsol
      INTEGER :: NMASDDEPSUCORSOL ! for drydep of SU from corsol
      INTEGER :: NMASDDEPSSACCSOL ! for drydep of SS from accsol
      INTEGER :: NMASDDEPSSCORSOL ! for drydep of SS from corsol
      INTEGER :: NMASDDEPBCAITSOL ! for drydep of BC from Aitsol
      INTEGER :: NMASDDEPBCACCSOL ! for drydep of BC from accsol
      INTEGER :: NMASDDEPBCCORSOL ! for drydep of BC from corsol
      INTEGER :: NMASDDEPBCAITINS ! for drydep of BC from Aitins
      INTEGER :: NMASDDEPOCNUCSOL ! for drydep of OC from nucsol
      INTEGER :: NMASDDEPOCAITSOL ! for drydep of OC from Aitsol
      INTEGER :: NMASDDEPOCACCSOL ! for drydep of OC from accsol
      INTEGER :: NMASDDEPOCCORSOL ! for drydep of OC from corsol
      INTEGER :: NMASDDEPOCAITINS ! for drydep of OC from Aitins
      INTEGER :: NMASDDEPSONUCSOL ! for drydep of SO from nucsol
      INTEGER :: NMASDDEPSOAITSOL ! for drydep of SO from Aitsol
      INTEGER :: NMASDDEPSOACCSOL ! for drydep of SO from accsol
      INTEGER :: NMASDDEPSOCORSOL ! for drydep of SO from corsol
      INTEGER :: NMASNUSCSUNUCSOL ! for nuscav of SU from nucsol
      INTEGER :: NMASNUSCSUAITSOL ! for nuscav of SU from Aitsol
      INTEGER :: NMASNUSCSUACCSOL ! for nuscav of SU from accsol
      INTEGER :: NMASNUSCSUCORSOL ! for nuscav of SU from corsol
      INTEGER :: NMASNUSCSSACCSOL ! for nuscav of SS from accsol
      INTEGER :: NMASNUSCSSCORSOL ! for nuscav of SS from corsol
      INTEGER :: NMASNUSCBCAITSOL ! for nuscav of BC from Aitsol
      INTEGER :: NMASNUSCBCACCSOL ! for nuscav of BC from accsol
      INTEGER :: NMASNUSCBCCORSOL ! for nuscav of BC from corsol
      INTEGER :: NMASNUSCBCAITINS ! for nuscav of BC from Aitins
      INTEGER :: NMASNUSCOCNUCSOL ! for nuscav of OC from nucsol
      INTEGER :: NMASNUSCOCAITSOL ! for nuscav of OC from Aitsol
      INTEGER :: NMASNUSCOCACCSOL ! for nuscav of OC from accsol
      INTEGER :: NMASNUSCOCCORSOL ! for nuscav of OC from corsol
      INTEGER :: NMASNUSCOCAITINS ! for nuscav of OC from Aitins
      INTEGER :: NMASNUSCSONUCSOL ! for nuscav of SO from nucsol
      INTEGER :: NMASNUSCSOAITSOL ! for nuscav of SO from Aitsol
      INTEGER :: NMASNUSCSOACCSOL ! for nuscav of SO from accsol
      INTEGER :: NMASNUSCSOCORSOL ! for nuscav of SO from corsol
      INTEGER :: NMASIMSCSUNUCSOL ! for imscav of SU from nucsol
      INTEGER :: NMASIMSCSUAITSOL ! for imscav of SU from Aitsol
      INTEGER :: NMASIMSCSUACCSOL ! for imscav of SU from accsol
      INTEGER :: NMASIMSCSUCORSOL ! for imscav of SU from corsol
      INTEGER :: NMASIMSCSSACCSOL ! for imscav of SS from accsol
      INTEGER :: NMASIMSCSSCORSOL ! for imscav of SS from corsol
      INTEGER :: NMASIMSCBCAITSOL ! for imscav of BC from Aitsol
      INTEGER :: NMASIMSCBCACCSOL ! for imscav of BC from accsol
      INTEGER :: NMASIMSCBCCORSOL ! for imscav of BC from corsol
      INTEGER :: NMASIMSCBCAITINS ! for imscav of BC from Aitins
      INTEGER :: NMASIMSCOCNUCSOL ! for imscav of OC from nucsol
      INTEGER :: NMASIMSCOCAITSOL ! for imscav of OC from Aitsol
      INTEGER :: NMASIMSCOCACCSOL ! for imscav of OC from accsol
      INTEGER :: NMASIMSCOCCORSOL ! for imscav of OC from corsol
      INTEGER :: NMASIMSCOCAITINS ! for imscav of OC from Aitins
      INTEGER :: NMASIMSCSONUCSOL ! for imscav of SO from nucsol
      INTEGER :: NMASIMSCSOAITSOL ! for imscav of SO from Aitsol
      INTEGER :: NMASIMSCSOACCSOL ! for imscav of SO from accsol
      INTEGER :: NMASIMSCSOCORSOL ! for imscav of SO from corsol
      INTEGER :: NMASCLPRSUAITSOL1 ! for in-cl. ox. of SO2 by H2O2 to SU in Aitsol
      INTEGER :: NMASCLPRSUACCSOL1 ! for in-cl. ox. of SO2 by H2O2 to SU in accsol
      INTEGER :: NMASCLPRSUCORSOL1 ! for in-cl. ox. of SO2 by H2O2 to SU in corsol
      INTEGER :: NMASCLPRSUAITSOL2 ! for in-cl. ox. of SO2 by O3   to SU in Aitsol
      INTEGER :: NMASCLPRSUACCSOL2 ! for in-cl. ox. of SO2 by O3   to SU in accsol
      INTEGER :: NMASCLPRSUCORSOL2 ! for in-cl. ox. of SO2 by O3   to SU in corsol
      INTEGER :: NMASCONDSUNUCSOL ! for conden. of SU to nucsol
      INTEGER :: NMASCONDSUAITSOL ! for conden. of SU to Aitsol
      INTEGER :: NMASCONDSUACCSOL ! for conden. of SU to accsol
      INTEGER :: NMASCONDSUCORSOL ! for conden. of SU to corsol
      INTEGER :: NMASCONDSUAITINS ! for conden. of SU to Aitins
      INTEGER :: NMASNUCLSUNUCSOL ! for nucln. of SU to nucsol
      INTEGER :: NMASCONDOCNUCSOL ! for conden. of OC to nucsol
      INTEGER :: NMASCONDOCAITSOL ! for conden. of OC to Aitsol
      INTEGER :: NMASCONDOCACCSOL ! for conden. of OC to accsol
      INTEGER :: NMASCONDOCCORSOL ! for conden. of OC to corsol
      INTEGER :: NMASCONDOCAITINS ! for conden. of OC to Aitins
      INTEGER :: NMASCONDSONUCSOL ! for conden. of SO to nucsol
      INTEGER :: NMASCONDSOAITSOL ! for conden. of SO to Aitsol
      INTEGER :: NMASCONDSOACCSOL ! for conden. of SO to accsol
      INTEGER :: NMASCONDSOCORSOL ! for conden. of SO to corsol
      INTEGER :: NMASCONDSOAITINS ! for conden. of SO to Aitins
      INTEGER :: NMASCOAGSUINTR12 ! for inter-modal coag SU nucsol->Aitsol
      INTEGER :: NMASCOAGSUINTR13 ! for inter-modal coag SU nucsol->accsol
      INTEGER :: NMASCOAGSUINTR14 ! for inter-modal coag SU nucsol->corsol
      INTEGER :: NMASCOAGSUINTR15 ! for inter-modal coag SU nucsol->Aitins
      INTEGER :: NMASCOAGOCINTR12 ! for inter-modal coag OC nucsol->Aitsol
      INTEGER :: NMASCOAGOCINTR13 ! for inter-modal coag OC nucsol->accsol
      INTEGER :: NMASCOAGOCINTR14 ! for inter-modal coag OC nucsol->corsol
      INTEGER :: NMASCOAGOCINTR15 ! for inter-modal coag OC nucsol->Aitins
      INTEGER :: NMASCOAGSOINTR12 ! for inter-modal coag SO nucsol->Aitsol
      INTEGER :: NMASCOAGSOINTR13 ! for inter-modal coag SO nucsol->accsol
      INTEGER :: NMASCOAGSOINTR14 ! for inter-modal coag SO nucsol->corsol
      INTEGER :: NMASCOAGSOINTR15 ! for inter-modal coag SO nucsol->Aitins
      INTEGER :: NMASCOAGSUINTR23 ! for inter-modal coag SU Aitsol->accsol
      INTEGER :: NMASCOAGBCINTR23 ! for inter-modal coag BC Aitsol->accsol
      INTEGER :: NMASCOAGOCINTR23 ! for inter-modal coag OC Aitsol->accsol
      INTEGER :: NMASCOAGSOINTR23 ! for inter-modal coag SO Aitsol->accsol
      INTEGER :: NMASCOAGSUINTR24 ! for inter-modal coag SU Aitsol->corsol
      INTEGER :: NMASCOAGBCINTR24 ! for inter-modal coag BC Aitsol->corsol
      INTEGER :: NMASCOAGOCINTR24 ! for inter-modal coag OC Aitsol->corsol
      INTEGER :: NMASCOAGSOINTR24 ! for inter-modal coag SO Aitsol->corsol
      INTEGER :: NMASCOAGSUINTR34 ! for inter-modal coag SU accsol->corsol
      INTEGER :: NMASCOAGBCINTR34 ! for inter-modal coag BC accsol->corsol
      INTEGER :: NMASCOAGOCINTR34 ! for inter-modal coag OC accsol->corsol
      INTEGER :: NMASCOAGSSINTR34 ! for inter-modal coag SS accsol->corsol
      INTEGER :: NMASCOAGSOINTR34 ! for inter-modal coag SO accsol->corsol
      INTEGER :: NMASCOAGBCINTR53 ! for inter-modal coag BC Aitins->accsol
      INTEGER :: NMASCOAGOCINTR53 ! for inter-modal coag OC Aitins->accsol
      INTEGER :: NMASCOAGBCINTR54 ! for inter-modal coag BC Aitins->corsol
      INTEGER :: NMASCOAGOCINTR54 ! for inter-modal coag OC Aitins->corsol
      INTEGER :: NMASAGEDSUINTR52 ! for SU ageing flux Aitins->Aitsol
      INTEGER :: NMASAGEDBCINTR52 ! for BC ageing flux Aitins->Aitsol
      INTEGER :: NMASAGEDOCINTR52 ! for OC ageing flux Aitins->Aitsol
      INTEGER :: NMASAGEDSOINTR52 ! for SO ageing flux Aitins->Aitsol
      INTEGER :: NMASMERGSUINTR12 ! for SU mode-merging flux nucsol->Aitsol
      INTEGER :: NMASMERGOCINTR12 ! for OC mode-merging flux nucsol->Aitsol
      INTEGER :: NMASMERGSOINTR12 ! for SO mode-merging flux nucsol->Aitsol
      INTEGER :: NMASMERGSUINTR23 ! for SU mode-merging flux Aitsol->accsol
      INTEGER :: NMASMERGBCINTR23 ! for BC mode-merging flux Aitsol->accsol
      INTEGER :: NMASMERGOCINTR23 ! for OC mode-merging flux Aitsol->accsol
      INTEGER :: NMASMERGSOINTR23 ! for SO mode-merging flux Aitsol->accsol
      INTEGER :: NMASMERGSUINTR34 ! for SU mode-merging flux accsol->corsol
      INTEGER :: NMASMERGSSINTR34 ! for SS mode-merging flux accsol->corsol
      INTEGER :: NMASMERGBCINTR34 ! for BC mode-merging flux accsol->corsol
      INTEGER :: NMASMERGOCINTR34 ! for OC mode-merging flux accsol->corsol
      INTEGER :: NMASMERGSOINTR34 ! for SO mode-merging flux accsol->corsol
      INTEGER :: NMASPROCSUINTR23 ! for SU cloud-processing Aitsol->accsol
      INTEGER :: NMASPROCBCINTR23 ! for BC cloud-processing Aitsol->accsol
      INTEGER :: NMASPROCOCINTR23 ! for OC cloud-processing Aitsol->accsol
      INTEGER :: NMASPROCSOINTR23 ! for SO cloud-processing Aitsol->accsol
!
! .. below are new ones for dust & modes 6/7 to be integrated
      INTEGER :: NMASPRIMDUACCSOL ! for dust emissions to accsol
      INTEGER :: NMASPRIMDUCORSOL ! for dust emissions to corsol
      INTEGER :: NMASPRIMDUACCINS ! for dust emissions to accins
      INTEGER :: NMASPRIMDUCORINS ! for dust emissions to corins
      INTEGER :: NMASDDEPDUACCSOL ! for dust dry dep from accsol
      INTEGER :: NMASDDEPDUCORSOL ! for dust dry dep from corsol
      INTEGER :: NMASDDEPDUACCINS ! for dust dry dep from accins
      INTEGER :: NMASDDEPDUCORINS ! for dust dry dep from corins
      INTEGER :: NMASNUSCDUACCSOL ! for dust nucscav from accsol
      INTEGER :: NMASNUSCDUCORSOL ! for dust nucscav from corsol
      INTEGER :: NMASNUSCDUACCINS ! for dust nucscav from accins
      INTEGER :: NMASNUSCDUCORINS ! for dust nucscav from corins
      INTEGER :: NMASIMSCDUACCSOL ! for dust impscav from accsol
      INTEGER :: NMASIMSCDUCORSOL ! for dust impscav from corsol
      INTEGER :: NMASIMSCDUACCINS ! for dust impscav from accins
      INTEGER :: NMASIMSCDUCORINS ! for dust impscav from corins
      INTEGER :: NMASCONDSUACCINS ! for conden. of SU to accins
      INTEGER :: NMASCONDSUCORINS ! for conden. of SU to corins
      INTEGER :: NMASCONDOCACCINS ! for conden. of OC to accins
      INTEGER :: NMASCONDOCCORINS ! for conden. of OC to corins
      INTEGER :: NMASCONDSOACCINS ! for conden. of SO to accins
      INTEGER :: NMASCONDSOCORINS ! for conden. of SO to corins
      INTEGER :: NMASCOAGSUINTR16 ! for inter-modal coag SU nucsol->accins
      INTEGER :: NMASCOAGSUINTR17 ! for inter-modal coag SU nucsol->corins
      INTEGER :: NMASCOAGOCINTR16 ! for inter-modal coag OC nucsol->accins
      INTEGER :: NMASCOAGOCINTR17 ! for inter-modal coag OC nucsol->corins
      INTEGER :: NMASCOAGSOINTR16 ! for inter-modal coag SO nucsol->accins
      INTEGER :: NMASCOAGSOINTR17 ! for inter-modal coag SO nucsol->corins
      INTEGER :: NMASCOAGDUINTR34 ! for inter-modal coag DU accsol->corsol
      INTEGER :: NMASCOAGDUINTR64 ! for inter-modal coag DU accins->corsol
      INTEGER :: NMASAGEDSUINTR63 ! for SU ageing flux accins->accsol
      INTEGER :: NMASAGEDDUINTR63 ! for DU ageing flux accins->accsol
      INTEGER :: NMASAGEDOCINTR63 ! for OC ageing flux accins->accsol
      INTEGER :: NMASAGEDSOINTR63 ! for SO ageing flux accins->accsol
      INTEGER :: NMASAGEDSUINTR74 ! for SU ageing flux corins->corsol
      INTEGER :: NMASAGEDDUINTR74 ! for DU ageing flux corins->corsol
      INTEGER :: NMASAGEDOCINTR74 ! for OC ageing flux corins->corsol
      INTEGER :: NMASAGEDSOINTR74 ! for SO ageing flux corins->corsol
      INTEGER :: NMASMERGDUINTR34 ! for DU mode-merging flux accsol->corsol
!
! .. (indices for additional diagnostics)
      INTEGER :: NWTCNT1
      INTEGER :: NWTCNT2
      INTEGER :: NWTCNT3
      INTEGER :: NWTCNT4
      INTEGER :: NLCFRAC
      INTEGER :: NLWC
      INTEGER :: NRAINLS
      INTEGER :: NRAINCS
      INTEGER :: NWINDSP
      INTEGER :: NDRYDP1
      INTEGER :: NDRYDP2
      INTEGER :: NDRYDP3
      INTEGER :: NDRYDP4
      INTEGER :: NDRYDP5
      INTEGER :: NDRYDP6
      INTEGER :: NDRYDP7
! .. note below are partial volumes for SUSSBCOCSO_5MODE
      INTEGER :: NPVOL11
      INTEGER :: NPVOL16
      INTEGER :: NPVOL1W
      INTEGER :: NPVOL21
      INTEGER :: NPVOL22
      INTEGER :: NPVOL23
      INTEGER :: NPVOL26
      INTEGER :: NPVOL2W
      INTEGER :: NPVOL31
      INTEGER :: NPVOL32
      INTEGER :: NPVOL33
      INTEGER :: NPVOL34
      INTEGER :: NPVOL36
      INTEGER :: NPVOL3W
      INTEGER :: NPVOL41
      INTEGER :: NPVOL42
      INTEGER :: NPVOL43
      INTEGER :: NPVOL44
      INTEGER :: NPVOL46
      INTEGER :: NPVOL4W
      INTEGER :: NPVOL52
      INTEGER :: NPVOL53
!
      INTEGER, PARAMETER :: NCHEMGMAX=50 ! max value for NCHEMG
!
! Indices of aerosol components into which condensable gases go to
      INTEGER :: CONDENSABLE_CHOICE(NCHEMGMAX)

! Gas phase species which are condensable (T/F)
      LOGICAL :: CONDENSABLE(NCHEMGMAX)

! Molecular masses of gas phase species (kg/mol)
      REAL :: MM_GAS(NCHEMGMAX)

! Molecular diameter of condensable gas phase species (others = 0)
      REAL :: DIMEN(NCHEMGMAX)

! Interface section

      INTERFACE UKCA_INDICES_NOCHEM
        MODULE PROCEDURE UKCA_INDICES_NOCHEM
      END INTERFACE UKCA_INDICES_NOCHEM

      INTERFACE UKCA_INDICES_ORGV1_SOto3
        MODULE PROCEDURE UKCA_INDICES_ORGV1_SOto3
      END INTERFACE UKCA_INDICES_ORGV1_SOto3

      INTERFACE UKCA_INDICES_ORGV1_SOto6
        MODULE PROCEDURE UKCA_INDICES_ORGV1_SOto6
      END INTERFACE UKCA_INDICES_ORGV1_SOto6

      INTERFACE UKCA_INDICES_SV1
        MODULE PROCEDURE UKCA_INDICES_SV1
      END INTERFACE UKCA_INDICES_SV1

      INTERFACE UKCA_INDICES_TRAQU38
        MODULE PROCEDURE UKCA_INDICES_TRAQU38
      END INTERFACE UKCA_INDICES_TRAQU38

      INTERFACE UKCA_INDICES_TRAQU9
        MODULE PROCEDURE UKCA_INDICES_TRAQU9
      END INTERFACE UKCA_INDICES_TRAQU9

      INTERFACE UKCA_INDICES_SUSSBCOC_5MODE
        MODULE PROCEDURE UKCA_INDICES_SUSSBCOC_5MODE
      END INTERFACE UKCA_INDICES_SUSSBCOC_5MODE

      INTERFACE UKCA_INDICES_SUSSBCOC_4MODE
        MODULE PROCEDURE UKCA_INDICES_SUSSBCOC_4MODE
      END INTERFACE UKCA_INDICES_SUSSBCOC_4MODE

      INTERFACE UKCA_INDICES_SUSSBCOCSO_5MODE
        MODULE PROCEDURE UKCA_INDICES_SUSSBCOCSO_5MODE
      END INTERFACE UKCA_INDICES_SUSSBCOCSO_5MODE

      INTERFACE UKCA_INDICES_SUSSBCOCSO_4MODE
        MODULE PROCEDURE UKCA_INDICES_SUSSBCOCSO_4MODE
      END INTERFACE UKCA_INDICES_SUSSBCOCSO_4MODE

      INTERFACE UKCA_INDICES_SUSS_4MODE
        MODULE PROCEDURE UKCA_INDICES_SUSS_4MODE
      END INTERFACE UKCA_INDICES_SUSS_4MODE

      INTERFACE UKCA_INDICES_DUonly_2MODE
        MODULE PROCEDURE UKCA_INDICES_DUonly_2MODE
      END INTERFACE UKCA_INDICES_DUonly_2MODE

      INTERFACE UKCA_INDICES_ORGV1_SOto3_COUPL
        MODULE PROCEDURE UKCA_INDICES_ORGV1_SOto3_COUPL
      END INTERFACE UKCA_INDICES_ORGV1_SOto3_COUPL

      INTERFACE UKCA_INDICES_ORGV1_SOto6_COUPL
        MODULE PROCEDURE UKCA_INDICES_ORGV1_SOto6_COUPL
      END INTERFACE UKCA_INDICES_ORGV1_SOto6_COUPL

      INTERFACE UKCA_INDICES_SV1_COUPLED
        MODULE PROCEDURE UKCA_INDICES_SV1_COUPLED
      END INTERFACE UKCA_INDICES_SV1_COUPLED

      INTERFACE UKCA_INDICES_SUSSBCOCDU_7MODE
        MODULE PROCEDURE UKCA_INDICES_SUSSBCOCDU_7MODE
      END INTERFACE UKCA_INDICES_SUSSBCOCDU_7MODE

      INTERFACE UKCA_INDICES_SUSSBCOCDU_4MODE
        MODULE PROCEDURE UKCA_INDICES_SUSSBCOCDU_4MODE
      END INTERFACE UKCA_INDICES_SUSSBCOCDU_4MODE

      CONTAINS

! ######################################################################
      SUBROUTINE UKCA_INDICES_NOCHEM

      IMPLICIT NONE

!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:UKCA_INDICES_NOCHEM',  &
                              zhook_in,zhook_handle)

      NCHEMG=0           ! # of gas phase chemistry tracers
      ICHEM=0            ! 1/0 = gas phase chemistry tracers on/off
      NOFFOX=0           ! # of offline oxidant species
      NBUDCHEM=0         ! # of gas chem budget fields
      GASBUDGET=0        ! 1/0 = gas phase chemistry fluxes on/off
      NGASBUDGET=0       ! # of gas phase chemistry fluxes
!
      NADVG=2+NCHEMG     ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX ! total # of gas phase species

! For nochem, NTRAG=2, NADVG= 2
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
!
! .. below are the 40 tropospheric chemistry species
      MOX     = 0 ! not included
      MNOX    = 0 ! not included
      MN2O5   = 0 ! not included
      MHNO4   = 0 ! not included
      MHNO3   = 0 ! not included
      MH2O2   = 0 ! not included
      MCH4    = 0 ! not included
      MCO     = 0 ! not included
      MCH2O   = 0 ! not included
      MMHP    = 0 ! not included
      MHONO   = 0 ! not included
      MC2H6   = 0 ! not included
      METOOH  = 0 ! not included
      MMECHO  = 0 ! not included
      MPAN    = 0 ! not included
      MC3H8   = 0 ! not included
      MPNOOH  = 0 ! not included
      MPIOOH  = 0 ! not included
      METCHO  = 0 ! not included
      MME2CO  = 0 ! not included
      MMECOH  = 0 ! not included
      MPPAN   = 0 ! not included
      MMENO3  = 0 ! not included
      MOXS    = 0 ! not included
      MNOYS   = 0 ! not included
      MISOP   = 0 ! not included
      MC2H4   = 0 ! not included
      MC2H2   = 0 ! not included
      MISOOH  = 0 ! not included
      MISON   = 0 ! not included
      MMACR   = 0 ! not included
      MMACROOH= 0 ! not included
      MMPAN   = 0 ! not included
      MHACET  = 0 ! not included
      MMGLY   = 0 ! not included
      MNALD   = 0 ! not included
      MHCOOH  = 0 ! not included
      MMECO3H = 0 ! not included
      MMECO2H = 0 ! not included
      MMEOH   = 0 ! not included
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
      MSOTWO  = 0 ! not included
      MMeSMe  = 0 ! not included
      MH2SO4  = 0 ! not included
      MDMSO   = 0 ! not included
      MMSA    = 0 ! not included
      MCS2    = 0 ! not included
      MH2S    = 0 ! not included
      MCOS    = 0 ! not included
      MMONOTER= 0 ! not included
      MSEC_ORG= 0 ! not included
      MH2O2F  = 0 ! not included
      MQ3D    = 1
      MPT     = 2
!
! .. molar masses (kg/mol) for gases for nochem 
! .. (48 dummy values)
      MM_GAS=(/0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000/)
!
      CONDENSABLE_CHOICE=(/0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0/)
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY
!
! .. below are the 60 indices for tropospheric chemistry species (ST)
      NO      = 0 ! not included
      NO1D    = 0 ! not included
      NO3     = 0 ! not included
      NNO     = 0 ! not included
      NNO3    = 0 ! not included
      NNO2    = 0 ! not included
      NN2O5   = 0 ! not included
      NHNO4   = 0 ! not included
      NHNO3   = 0 ! not included
      NOH     = 0 ! not included
      NHO2    = 0 ! not included
      NH2O2   = 0 ! not included
      NCH4    = 0 ! not included
      NCO     = 0 ! not included
      NCH2O   = 0 ! not included
      NMEOO   = 0 ! not included
      NH2O    = 0 ! not included
      NMHP    = 0 ! not included
      NHONO   = 0 ! not included
      NC2H6   = 0 ! not included
      NETOO   = 0 ! not included
      NETOOH  = 0 ! not included
      NMECHO  = 0 ! not included
      NMECO3  = 0 ! not included
      NPAN    = 0 ! not included
      NC3H8   = 0 ! not included
      NPNOO   = 0 ! not included
      NPIOO   = 0 ! not included
      NPNOOH  = 0 ! not included
      NPIOOH  = 0 ! not included
      NETCHO  = 0 ! not included
      NETCO3  = 0 ! not included
      NME2CO  = 0 ! not included
      NMECOO  = 0 ! not included
      NMECOH  = 0 ! not included
      NPPAN   = 0 ! not included
      NMENO3  = 0 ! not included
      NOS     = 0 ! not included
      NO1DS   = 0 ! not included
      NO3S    = 0 ! not included
      NNOXS   = 0 ! not included
      NHNO3S  = 0 ! not included
      NNOYS   = 0 ! not included
      NISOP   = 0 ! not included
      NC2H4   = 0 ! not included
      NC2H2   = 0 ! not included
      NISO2   = 0 ! not included
      NISOOH  = 0 ! not included
      NISON   = 0 ! not included
      NMACR   = 0 ! not included
      NMACRO2 = 0 ! not included
      NMACROOH= 0 ! not included
      NMPAN   = 0 ! not included
      NHACET  = 0 ! not included
      NMGLY   = 0 ! not included
      NNALD   = 0 ! not included
      NHCOOH  = 0 ! not included
      NMECO3H = 0 ! not included
      NMECO2H = 0 ! not included
      NMEOH   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
      NSOTWO  = 0 ! not included
      NMeSMe  = 0 ! not included
      NH2SO4  = 0 ! not included
      NDMSO   = 0 ! not included
      NMSA    = 0 ! not included
      NCS2    = 0 ! not included
      NH2S    = 0 ! not included
      NCOS    = 0 ! not included
      NMONOTER= 0 ! not included
      NSEC_ORG= 0 ! not included
      NH2O2F  = 0 ! not included
      NO3F    = 0 ! not included
      NOHF    = 0 ! not included
      NNO3F   = 0 ! not included
      NQ3D    = 1
      NPT     = 2
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for nochem
      NDMSEMOC  = 0 ! not included
      NDMSTEND  = 0 ! not included
      NSO2EMAN  = 0 ! not included
      NSO2EMBM  = 0 ! not included
      NSO2EMVL  = 0 ! not included
      NSO2TEND  = 0 ! not included
      NSO2DDEP  = 0 ! not included
      NSO2WDEP  = 0 ! not included
      NH2SO4TEND= 0 ! not included
      NH2SO4DDEP= 0 ! not included
      NCOSEMAN  = 0 ! not included
      NCOSEMOC  = 0 ! not included
      NCOSTEND  = 0 ! not included
      NCS2EMAN  = 0 ! not included
      NCS2EMOC  = 0 ! not included
      NCS2TEND  = 0 ! not included
      NDMSOTEND = 0 ! not included
      NDMSODDEP = 0 ! not included
      NMSATEND  = 0 ! not included
      NMSADDEP  = 0 ! not included
      NTERP_EM  = 0 ! not included
      NTERP_TEND= 0 ! not included
      NTERP_DDEP= 0 ! not included
      NSORG_TEND= 0 ! not included
      NSORG_DDEP= 0 ! not included
      NSORG_WDEP= 0 ! not included
!
! REACTION INDICES for gas phase chemistry fluxes
      IOHDMS1 =  0
      IOHDMS2 =  0
      INO3DMS =  0
      IDMSOOH1=  0
      IDMSOOH2=  0
      ICS2OH  =  0
      IH2SOH  =  0
      ICOSOH  =  0
 
      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:UKCA_INDICES_NOCHEM',  &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_NOCHEM

! ######################################################################
      SUBROUTINE UKCA_INDICES_ORGV1_SOto3

      IMPLICIT NONE

!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:UKCA_INDICES_ORGV1_SOTO3', &
                               zhook_in,zhook_handle)

      NCHEMG=11          ! # of gas phase chemistry tracers
      ICHEM=1            ! 1/0 = gas phase chemistry tracers on/off
      NOFFOX=3           ! # of offline oxidant species
      NBUDCHEM=26        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
      GASBUDGET=0        ! 1/0 = gas phase chemistry fluxes on/off
      NGASBUDGET=8       ! # of gas phase chemistry fluxes
!
      NADVG=2+NCHEMG     ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX ! total # of gas phase species

! For orgv1, NTRAG=16, NADVG=13
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
!
! .. below are the 40 tropospheric chemistry species
      MOX     = 0 ! not included
      MNOX    = 0 ! not included
      MN2O5   = 0 ! not included
      MHNO4   = 0 ! not included
      MHNO3   = 0 ! not included
      MH2O2   = 0 ! not included
      MCH4    = 0 ! not included
      MCO     = 0 ! not included
      MCH2O   = 0 ! not included
      MMHP    = 0 ! not included
      MHONO   = 0 ! not included
      MC2H6   = 0 ! not included
      METOOH  = 0 ! not included
      MMECHO  = 0 ! not included
      MPAN    = 0 ! not included
      MC3H8   = 0 ! not included
      MPNOOH  = 0 ! not included
      MPIOOH  = 0 ! not included
      METCHO  = 0 ! not included
      MME2CO  = 0 ! not included
      MMECOH  = 0 ! not included
      MPPAN   = 0 ! not included
      MMENO3  = 0 ! not included
      MOXS    = 0 ! not included
      MNOYS   = 0 ! not included
      MISOP   = 0 ! not included
      MC2H4   = 0 ! not included
      MC2H2   = 0 ! not included
      MISOOH  = 0 ! not included
      MISON   = 0 ! not included
      MMACR   = 0 ! not included
      MMACROOH= 0 ! not included
      MMPAN   = 0 ! not included
      MHACET  = 0 ! not included
      MMGLY   = 0 ! not included
      MNALD   = 0 ! not included
      MHCOOH  = 0 ! not included
      MMECO3H = 0 ! not included
      MMECO2H = 0 ! not included
      MMEOH   = 0 ! not included
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
      MSOTWO  = 1
      MMeSMe  = 2
      MH2SO4  = 3
      MDMSO   = 4
      MMSA    = 5
      MCS2    = 6
      MH2S    = 7
      MCOS    = 8
      MMONOTER= 9
      MSEC_ORG=10
      MH2O2F  =11
      MQ3D    =12
      MPT     =13
!
! .. molar masses (kg/mol) for gases for orgv1
! .. (8 S species, terp, sec_org, H2O2F then 37 dummy values)
      MM_GAS=(/0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,         &
               0.136,0.150,0.034,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000/)
!
      CONDENSABLE_CHOICE=(/0,0,1,0,0,0,0,0,0,3,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. SEC_ORG to condense into 3rd aerosol component (CP_OC)

      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,  &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0/)
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY
!
! .. below are the 60 indices for tropospheric chemistry species (ST)
      NO      = 0 ! not included
      NO1D    = 0 ! not included
      NO3     = 0 ! not included
      NNO     = 0 ! not included
      NNO3    = 0 ! not included
      NNO2    = 0 ! not included
      NN2O5   = 0 ! not included
      NHNO4   = 0 ! not included
      NHNO3   = 0 ! not included
      NOH     = 0 ! not included
      NHO2    = 0 ! not included
      NH2O2   = 0 ! not included
      NCH4    = 0 ! not included
      NCO     = 0 ! not included
      NCH2O   = 0 ! not included
      NMEOO   = 0 ! not included
      NH2O    = 0 ! not included
      NMHP    = 0 ! not included
      NHONO   = 0 ! not included
      NC2H6   = 0 ! not included
      NETOO   = 0 ! not included
      NETOOH  = 0 ! not included
      NMECHO  = 0 ! not included
      NMECO3  = 0 ! not included
      NPAN    = 0 ! not included
      NC3H8   = 0 ! not included
      NPNOO   = 0 ! not included
      NPIOO   = 0 ! not included
      NPNOOH  = 0 ! not included
      NPIOOH  = 0 ! not included
      NETCHO  = 0 ! not included
      NETCO3  = 0 ! not included
      NME2CO  = 0 ! not included
      NMECOO  = 0 ! not included
      NMECOH  = 0 ! not included
      NPPAN   = 0 ! not included
      NMENO3  = 0 ! not included
      NOS     = 0 ! not included
      NO1DS   = 0 ! not included
      NO3S    = 0 ! not included
      NNOXS   = 0 ! not included
      NHNO3S  = 0 ! not included
      NNOYS   = 0 ! not included
      NISOP   = 0 ! not included
      NC2H4   = 0 ! not included
      NC2H2   = 0 ! not included
      NISO2   = 0 ! not included
      NISOOH  = 0 ! not included
      NISON   = 0 ! not included
      NMACR   = 0 ! not included
      NMACRO2 = 0 ! not included
      NMACROOH= 0 ! not included
      NMPAN   = 0 ! not included
      NHACET  = 0 ! not included
      NMGLY   = 0 ! not included
      NNALD   = 0 ! not included
      NHCOOH  = 0 ! not included
      NMECO3H = 0 ! not included
      NMECO2H = 0 ! not included
      NMEOH   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
      NSOTWO  = 1
      NMeSMe  = 2
      NH2SO4  = 3
      NDMSO   = 4
      NMSA    = 5
      NCS2    = 6
      NH2S    = 7
      NCOS    = 8
      NMONOTER= 9
      NSEC_ORG=10
      NH2O2F  =11
      NO3F    =12
      NOHF    =13
      NNO3F   =14
      NQ3D    =15
      NPT     =16
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
      NDMSEMOC  = 1
      NDMSTEND  = 2
      NSO2EMAN  = 3
      NSO2EMBM  = 4
      NSO2EMVL  = 5
      NSO2TEND  = 6
      NSO2DDEP  = 7
      NSO2WDEP  = 8
      NH2SO4TEND= 9
      NH2SO4DDEP=10
      NCOSEMAN  =11
      NCOSEMOC  =12
      NCOSTEND  =13
      NCS2EMAN  =14
      NCS2EMOC  =15
      NCS2TEND  =16
      NDMSOTEND =17
      NDMSODDEP =18
      NMSATEND  =19
      NMSADDEP  =20
      NTERP_EM  =21
      NTERP_TEND=22
      NTERP_DDEP=23
      NSORG_TEND=24
      NSORG_DDEP=25
      NSORG_WDEP=26
!
! REACTION INDICES for gas phase chemistry fluxes
      IOHDMS1 =  1
      IOHDMS2 =  2
      INO3DMS =  4
      IDMSOOH1=  0
      IDMSOOH2=  3
      ICS2OH  =  5
      IH2SOH  =  6
      ICOSOH  =  7
! 
!----------------------------------------------------------------

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:UKCA_INDICES_ORGV1_SOTO3', &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_ORGV1_SOto3

! ######################################################################
      SUBROUTINE UKCA_INDICES_ORGV1_SOto3_COUPL

      IMPLICIT NONE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
           'UKCA_INDICES_ORGV1_SOto3_COUPL',zhook_in,zhook_handle)

! Main array lengths and switches

      NCHEMG=50          ! # of gas phase chemistry tracers
      ICHEM=1            ! 1/0 = gas phase chemistry tracers on/off
      NOFFOX=24          ! # of offline oxidant species
      NBUDCHEM=26        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
      GASBUDGET=0        ! 1/0 = gas phase chemistry fluxes on/off
      NGASBUDGET=8       ! # of gas phase chemistry fluxes

      NADVG=2+NCHEMG     ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX ! total # of gas phase species


!   orgv1_soto3_coupled: NTRAG=76, NADVG=52
!
!-----------------------------------------------------------

! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
      MOX     = 1
      MNOX    = 2
      MN2O5   = 3
      MHNO4   = 4
      MHNO3   = 5
      MH2O2   = 6
      MCH4    = 7
      MCO     = 8
      MCH2O   = 9
      MMHP    =10
      MHONO   =11
      MC2H6   =12
      METOOH  =13
      MMECHO  =14
      MPAN    =15
      MC3H8   =16
      MPNOOH  =17
      MPIOOH  =18
      METCHO  =19
      MME2CO  =20
      MMECOH  =21
      MPPAN   =22
      MMENO3  =23
      MOXS    =24
      MNOYS   =25
      MISOP   =26
      MC2H4   =27
      MC2H2   =28
      MISOOH  =29
      MISON   =30
      MMACR   =31
      MMACROOH=32
      MMPAN   =33
      MHACET  =34
      MMGLY   =35
      MNALD   =36
      MHCOOH  =37
      MMECO3H =38
      MMECO2H =39
      MMEOH   =40
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
      MSOTWO  =41
      MMeSMe  =42
      MH2SO4  =43
      MDMSO   =44
      MMSA    =45
      MCS2    =46
      MH2S    =47
      MCOS    =48
      MMONOTER=49
      MSEC_ORG=50
      MH2O2F  = 0 ! not included
      MQ3D    =51
      MPT     =52
!
! .. molar masses (kg/mol) for gases for sv1_coupled
! .. (40 tropospheric gases then 8 sulphur species)
!     Dummy values used for OX =1, NOX  =2, MHP  =10,ETOOH=13, MECHO=14
!     species & no          PAN=15 PNOOH=17,PIOOH=18,ETCHO=19, ME2CO=20
!                           MECOH=21, PPAN=22, MENO3=23, OXS=24
!                           NOYS=25, ISOP=26
      MM_GAS=(/0.048,0.030,0.108,0.079,0.053,0.034,0.016,0.013,         &
               0.030,0.033,0.047,0.030,0.056,0.045,0.088,0.044,         &
               0.045,0.067,0.056,0.087,0.067,0.051,0.034,0.013,         &
               0.018,0.045,0.028,0.026,0.015,0.015,0.015,0.015,         &
               0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,         &
               0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,         &
               0.136,0.150/)

      CONDENSABLE_CHOICE=(/0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,1,0,0,0,0,0,0,3/)
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. SEC_ORG to condense into 3rd aerosol component (CP_OC)
!
      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,      &
              0.0,4.5e-10/)
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY

! .. below are the 60 indices for tropospheric chemistry species (ST)
      NO      = 1
      NO1D    = 2
      NO3     = 3
      NNO     = 4
      NNO3    = 5
      NNO2    = 6
      NN2O5   = 7
      NHNO4   = 8
      NHNO3   = 9
      NOH     =10
      NHO2    =11
      NH2O2   =12
      NCH4    =13
      NCO     =14
      NCH2O   =15
      NMEOO   =16
      NH2O    =17
      NMHP    =18
      NHONO   =19
      NC2H6   =20
      NETOO   =21
      NETOOH  =22
      NMECHO  =23
      NMECO3  =24
      NPAN    =25
      NC3H8   =26
      NPNOO   =27
      NPIOO   =28
      NPNOOH  =29
      NPIOOH  =30
      NETCHO  =31
      NETCO3  =32
      NME2CO  =33
      NMECOO  =34
      NMECOH  =35
      NPPAN   =36
      NMENO3  =37
      NOS     =38
      NO1DS   =39
      NO3S    =40
      NNOXS   =41
      NHNO3S  =42
      NNOYS   =43
      NISOP   =44
      NC2H4   =45
      NC2H2   =46
      NISO2   =47
      NISOOH  =48
      NISON   =49
      NMACR   =50
      NMACRO2 =51
      NMACROOH=52
      NMPAN   =53
      NHACET  =54
      NMGLY   =55
      NNALD   =56
      NHCOOH  =57
      NMECO3H =58
      NMECO2H =59
      NMEOH   =60
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
      NSOTWO  =61
      NMeSMe  =62
      NH2SO4  =63
      NDMSO   =64
      NMSA    =65
      NCS2    =66
      NH2S    =67
      NCOS    =68
      NMONOTER=69
      NSEC_ORG=70
      NH2O2F  =76
      NO3F    =73
      NOHF    =74
      NNO3F   =75
      NQ3D    =71
      NPT     =72
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
      NDMSEMOC  = 1
      NDMSTEND  = 2
      NSO2EMAN  = 3
      NSO2EMBM  = 4
      NSO2EMVL  = 5
      NSO2TEND  = 6
      NSO2DDEP  = 7
      NSO2WDEP  = 8
      NH2SO4TEND= 9
      NH2SO4DDEP=10
      NCOSEMAN  =11
      NCOSEMOC  =12
      NCOSTEND  =13
      NCS2EMAN  =14
      NCS2EMOC  =15
      NCS2TEND  =16
      NDMSOTEND =17
      NDMSODDEP =18
      NMSATEND  =19
      NMSADDEP  =20
      NTERP_EM  =21
      NTERP_TEND=22
      NTERP_DDEP=23
      NSORG_TEND=24
      NSORG_DDEP=25
      NSORG_WDEP=26
!
! REACTION INDICES for gas phase chemistry fluxes
      IOHDMS1 =123
      IOHDMS2 =124
      INO3DMS =127
      IDMSOOH1=125
      IDMSOOH2=126
      ICS2OH  =128
      IH2SOH  =129
      ICOSOH  =130
! 
      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
           'UKCA_INDICES_ORGV1_SOto3_COUPL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_ORGV1_SOto3_COUPL

! ######################################################################
      SUBROUTINE UKCA_INDICES_ORGV1_SOto6

      IMPLICIT NONE

!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:UKCA_INDICES_ORGV1_SOTO6', &
                               zhook_in,zhook_handle)

      NCHEMG=11          ! # of gas phase chemistry tracers
      ICHEM=1            ! 1/0 = gas phase chemistry tracers on/off
      NOFFOX=3           ! # of offline oxidant species
      NBUDCHEM=26        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
      GASBUDGET=0        ! 1/0 = gas phase chemistry fluxes on/off
      NGASBUDGET=8       ! # of gas phase chemistry fluxes
!
      NADVG=2+NCHEMG     ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX ! total # of gas phase species

! For orgv1, NTRAG=16, NADVG=13
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
      MOX     = 0 ! not included
      MNOX    = 0 ! not included
      MN2O5   = 0 ! not included
      MHNO4   = 0 ! not included
      MHNO3   = 0 ! not included
      MH2O2   = 0 ! not included
      MCH4    = 0 ! not included
      MCO     = 0 ! not included
      MCH2O   = 0 ! not included
      MMHP    = 0 ! not included
      MHONO   = 0 ! not included
      MC2H6   = 0 ! not included
      METOOH  = 0 ! not included
      MMECHO  = 0 ! not included
      MPAN    = 0 ! not included
      MC3H8   = 0 ! not included
      MPNOOH  = 0 ! not included
      MPIOOH  = 0 ! not included
      METCHO  = 0 ! not included
      MME2CO  = 0 ! not included
      MMECOH  = 0 ! not included
      MPPAN   = 0 ! not included
      MMENO3  = 0 ! not included
      MOXS    = 0 ! not included
      MNOYS   = 0 ! not included
      MISOP   = 0 ! not included
      MC2H4   = 0 ! not included
      MC2H2   = 0 ! not included
      MISOOH  = 0 ! not included
      MISON   = 0 ! not included
      MMACR   = 0 ! not included
      MMACROOH= 0 ! not included
      MMPAN   = 0 ! not included
      MHACET  = 0 ! not included
      MMGLY   = 0 ! not included
      MNALD   = 0 ! not included
      MHCOOH  = 0 ! not included
      MMECO3H = 0 ! not included
      MMECO2H = 0 ! not included
      MMEOH   = 0 ! not included
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
      MSOTWO  = 1
      MMeSMe  = 2
      MH2SO4  = 3
      MDMSO   = 4
      MMSA    = 5
      MCS2    = 6
      MH2S    = 7
      MCOS    = 8
      MMONOTER= 9
      MSEC_ORG=10
      MH2O2F  =11
      MQ3D    =12
      MPT     =13
!
! .. molar masses (kg/mol) for gases for orgv1
! .. (8 S species, terp, sec_org, H2O2F then 37 dummy values)
      MM_GAS=(/0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,         &
               0.136,0.150,0.034,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000/)
!
      CONDENSABLE_CHOICE=(/0,0,1,0,0,0,0,0,0,6,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. SEC_ORG to condense into 6th aerosol component (CP_SO)

      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,  &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0/)
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY
!
! .. below are the 60 indices for tropospheric chemistry species (ST)
      NO      = 0 ! not included
      NO1D    = 0 ! not included
      NO3     = 0 ! not included
      NNO     = 0 ! not included
      NNO3    = 0 ! not included
      NNO2    = 0 ! not included
      NN2O5   = 0 ! not included
      NHNO4   = 0 ! not included
      NHNO3   = 0 ! not included
      NOH     = 0 ! not included
      NHO2    = 0 ! not included
      NH2O2   = 0 ! not included
      NCH4    = 0 ! not included
      NCO     = 0 ! not included
      NCH2O   = 0 ! not included
      NMEOO   = 0 ! not included
      NH2O    = 0 ! not included
      NMHP    = 0 ! not included
      NHONO   = 0 ! not included
      NC2H6   = 0 ! not included
      NETOO   = 0 ! not included
      NETOOH  = 0 ! not included
      NMECHO  = 0 ! not included
      NMECO3  = 0 ! not included
      NPAN    = 0 ! not included
      NC3H8   = 0 ! not included
      NPNOO   = 0 ! not included
      NPIOO   = 0 ! not included
      NPNOOH  = 0 ! not included
      NPIOOH  = 0 ! not included
      NETCHO  = 0 ! not included
      NETCO3  = 0 ! not included
      NME2CO  = 0 ! not included
      NMECOO  = 0 ! not included
      NMECOH  = 0 ! not included
      NPPAN   = 0 ! not included
      NMENO3  = 0 ! not included
      NOS     = 0 ! not included
      NO1DS   = 0 ! not included
      NO3S    = 0 ! not included
      NNOXS   = 0 ! not included
      NHNO3S  = 0 ! not included
      NNOYS   = 0 ! not included
      NISOP   = 0 ! not included
      NC2H4   = 0 ! not included
      NC2H2   = 0 ! not included
      NISO2   = 0 ! not included
      NISOOH  = 0 ! not included
      NISON   = 0 ! not included
      NMACR   = 0 ! not included
      NMACRO2 = 0 ! not included
      NMACROOH= 0 ! not included
      NMPAN   = 0 ! not included
      NHACET  = 0 ! not included
      NMGLY   = 0 ! not included
      NNALD   = 0 ! not included
      NHCOOH  = 0 ! not included
      NMECO3H = 0 ! not included
      NMECO2H = 0 ! not included
      NMEOH   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
      NSOTWO  = 1
      NMeSMe  = 2
      NH2SO4  = 3
      NDMSO   = 4
      NMSA    = 5
      NCS2    = 6
      NH2S    = 7
      NCOS    = 8
      NMONOTER= 9
      NSEC_ORG=10
      NH2O2F  =11
      NO3F    =12
      NOHF    =13
      NNO3F   =14
      NQ3D    =15
      NPT     =16
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
      NDMSEMOC  = 1
      NDMSTEND  = 2
      NSO2EMAN  = 3
      NSO2EMBM  = 4
      NSO2EMVL  = 5
      NSO2TEND  = 6
      NSO2DDEP  = 7
      NSO2WDEP  = 8
      NH2SO4TEND= 9
      NH2SO4DDEP=10
      NCOSEMAN  =11
      NCOSEMOC  =12
      NCOSTEND  =13
      NCS2EMAN  =14
      NCS2EMOC  =15
      NCS2TEND  =16
      NDMSOTEND =17
      NDMSODDEP =18
      NMSATEND  =19
      NMSADDEP  =20
      NTERP_EM  =21
      NTERP_TEND=22
      NTERP_DDEP=23
      NSORG_TEND=24
      NSORG_DDEP=25
      NSORG_WDEP=26

! REACTION INDICES for gas phase chemistry fluxes
      IOHDMS1 =  1
      IOHDMS2 =  2
      INO3DMS =  4
      IDMSOOH1=  0
      IDMSOOH2=  3
      ICS2OH  =  5
      IH2SOH  =  6
      ICOSOH  =  7
! 
!----------------------------------------------------------------

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:UKCA_INDICES_ORGV1_SOTO6', &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_ORGV1_SOto6

! ######################################################################
      SUBROUTINE UKCA_INDICES_ORGV1_SOto6_COUPL

      IMPLICIT NONE

!-----------------------------------------------------------
! Main array lengths and switches

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_ORGV1_SOTO6_COUPL',zhook_in,zhook_handle)

      NCHEMG=50          ! # of gas phase chemistry tracers
      ICHEM=1            ! 1/0 = gas phase chemistry tracers on/off
      NOFFOX=24          ! # of offline oxidant species
      NBUDCHEM=26        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
      GASBUDGET=0        ! 1/0 = gas phase chemistry fluxes on/off
      NGASBUDGET=8       ! # of gas phase chemistry fluxes
!
      NADVG=2+NCHEMG     ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX ! total # of gas phase species
!
!   orgv1_soto3_coupled: NTRAG=76, NADVG=52
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
      MOX     = 1
      MNOX    = 2
      MN2O5   = 3
      MHNO4   = 4
      MHNO3   = 5
      MH2O2   = 6
      MCH4    = 7
      MCO     = 8
      MCH2O   = 9
      MMHP    =10
      MHONO   =11
      MC2H6   =12
      METOOH  =13
      MMECHO  =14
      MPAN    =15
      MC3H8   =16
      MPNOOH  =17
      MPIOOH  =18
      METCHO  =19
      MME2CO  =20
      MMECOH  =21
      MPPAN   =22
      MMENO3  =23
      MOXS    =24
      MNOYS   =25
      MISOP   =26
      MC2H4   =27
      MC2H2   =28
      MISOOH  =29
      MISON   =30
      MMACR   =31
      MMACROOH=32
      MMPAN   =33
      MHACET  =34
      MMGLY   =35
      MNALD   =36
      MHCOOH  =37
      MMECO3H =38
      MMECO2H =39
      MMEOH   =40
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
      MSOTWO  =41
      MMeSMe  =42
      MH2SO4  =43
      MDMSO   =44
      MMSA    =45
      MCS2    =46
      MH2S    =47
      MCOS    =48
      MMONOTER=49
      MSEC_ORG=50
      MH2O2F  = 0 ! not included
      MQ3D    =51
      MPT     =52
!
! .. molar masses (kg/mol) for gases for sv1_coupled
! .. (40 tropospheric gases then 8 sulphur species)
!     Dummy values used for OX =1, NOX  =2, MHP  =10,ETOOH=13, MECHO=14
!     species & no          PAN=15 PNOOH=17,PIOOH=18,ETCHO=19, ME2CO=20
!                           MECOH=21, PPAN=22, MENO3=23, OXS=24
!                           NOYS=25, ISOP=26
      MM_GAS=(/0.048,0.030,0.108,0.079,0.053,0.034,0.016,0.013,         &
               0.030,0.033,0.047,0.030,0.056,0.045,0.088,0.044,         &
               0.045,0.067,0.056,0.087,0.067,0.051,0.034,0.013,         &
               0.018,0.045,0.028,0.026,0.015,0.015,0.015,0.015,         &
               0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,         &
               0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,         &
               0.136,0.150/)

      CONDENSABLE_CHOICE=(/0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,1,0,0,0,0,0,0,6/)
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. SEC_ORG to condense into 6th aerosol component (CP_SO)
!
      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,      &
              0.0,4.5e-10/)
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY

! .. below are the 60 indices for tropospheric chemistry species (ST)
      NO      = 1
      NO1D    = 2
      NO3     = 3
      NNO     = 4
      NNO3    = 5
      NNO2    = 6
      NN2O5   = 7
      NHNO4   = 8
      NHNO3   = 9
      NOH     =10
      NHO2    =11
      NH2O2   =12
      NCH4    =13
      NCO     =14
      NCH2O   =15
      NMEOO   =16
      NH2O    =17
      NMHP    =18
      NHONO   =19
      NC2H6   =20
      NETOO   =21
      NETOOH  =22
      NMECHO  =23
      NMECO3  =24
      NPAN    =25
      NC3H8   =26
      NPNOO   =27
      NPIOO   =28
      NPNOOH  =29
      NPIOOH  =30
      NETCHO  =31
      NETCO3  =32
      NME2CO  =33
      NMECOO  =34
      NMECOH  =35
      NPPAN   =36
      NMENO3  =37
      NOS     =38
      NO1DS   =39
      NO3S    =40
      NNOXS   =41
      NHNO3S  =42
      NNOYS   =43
      NISOP   =44
      NC2H4   =45
      NC2H2   =46
      NISO2   =47
      NISOOH  =48
      NISON   =49
      NMACR   =50
      NMACRO2 =51
      NMACROOH=52
      NMPAN   =53
      NHACET  =54
      NMGLY   =55
      NNALD   =56
      NHCOOH  =57
      NMECO3H =58
      NMECO2H =59
      NMEOH   =60
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
      NSOTWO  =61
      NMeSMe  =62
      NH2SO4  =63
      NDMSO   =64
      NMSA    =65
      NCS2    =66
      NH2S    =67
      NCOS    =68
      NMONOTER=69
      NSEC_ORG=70
      NH2O2F  =76
      NO3F    =73
      NOHF    =74
      NNO3F   =75
      NQ3D    =71
      NPT     =72
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
      NDMSEMOC  = 1
      NDMSTEND  = 2
      NSO2EMAN  = 3
      NSO2EMBM  = 4
      NSO2EMVL  = 5
      NSO2TEND  = 6
      NSO2DDEP  = 7
      NSO2WDEP  = 8
      NH2SO4TEND= 9
      NH2SO4DDEP=10
      NCOSEMAN  =11
      NCOSEMOC  =12
      NCOSTEND  =13
      NCS2EMAN  =14
      NCS2EMOC  =15
      NCS2TEND  =16
      NDMSOTEND =17
      NDMSODDEP =18
      NMSATEND  =19
      NMSADDEP  =20
      NTERP_EM  =21
      NTERP_TEND=22
      NTERP_DDEP=23
      NSORG_TEND=24
      NSORG_DDEP=25
      NSORG_WDEP=26

! REACTION INDICES for gas phase chemistry fluxes
      IOHDMS1 =123
      IOHDMS2 =124
      INO3DMS =127
      IDMSOOH1=125
      IDMSOOH2=126
      ICS2OH  =128
      IH2SOH  =129
      ICOSOH  =130
 
      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_ORGV1_SOTO6_COUPL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_ORGV1_SOto6_COUPL

! ######################################################################
      SUBROUTINE UKCA_INDICES_SV1
!
      IMPLICIT NONE
!
!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:UKCA_INDICES_SV1',     &
                               zhook_in,zhook_handle)

      NCHEMG=9           ! # of gas phase chemistry tracers
      ICHEM=1            ! 1/0 = gas phase chemistry tracers on/off
      NOFFOX=2           ! # of offline oxidant species
      NBUDCHEM=20        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
      GASBUDGET=0        ! 1/0 = gas phase chemistry fluxes on/off
      NGASBUDGET=8       ! # of gas phase chemistry fluxes
!
      NADVG=2+NCHEMG     ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX ! total # of gas phase species
!
!   sv1: NTRAG=13, NADVG=11
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
      MOX     = 0 ! not included
      MNOX    = 0 ! not included
      MN2O5   = 0 ! not included
      MHNO4   = 0 ! not included
      MHNO3   = 0 ! not included
      MH2O2   = 0 ! not included
      MCH4    = 0 ! not included
      MCO     = 0 ! not included
      MCH2O   = 0 ! not included
      MMHP    = 0 ! not included
      MHONO   = 0 ! not included
      MC2H6   = 0 ! not included
      METOOH  = 0 ! not included
      MMECHO  = 0 ! not included
      MPAN    = 0 ! not included
      MC3H8   = 0 ! not included
      MPNOOH  = 0 ! not included
      MPIOOH  = 0 ! not included
      METCHO  = 0 ! not included
      MME2CO  = 0 ! not included
      MMECOH  = 0 ! not included
      MPPAN   = 0 ! not included
      MMENO3  = 0 ! not included
      MOXS    = 0 ! not included
      MNOYS   = 0 ! not included
      MISOP   = 0 ! not included
      MC2H4   = 0 ! not included
      MC2H2   = 0 ! not included
      MISOOH  = 0 ! not included
      MISON   = 0 ! not included
      MMACR   = 0 ! not included
      MMACROOH= 0 ! not included
      MMPAN   = 0 ! not included
      MHACET  = 0 ! not included
      MMGLY   = 0 ! not included
      MNALD   = 0 ! not included
      MHCOOH  = 0 ! not included
      MMECO3H = 0 ! not included
      MMECO2H = 0 ! not included
      MMEOH   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
      MSOTWO  = 1
      MMeSMe  = 2
      MH2SO4  = 3
      MDMSO   = 4
      MMSA    = 5
      MCS2    = 6
      MH2S    = 7
      MCOS    = 8
      MMONOTER= 0 ! not included
      MSEC_ORG= 0 ! not included
      MH2O2F  = 9
      MQ3D    =10
      MPT     =11
!
! .. molar masses (kg/mol) for gases for sv1
! .. (8 S species, H2O2F then 39 dummy values)
      MM_GAS=(/0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,         &
               0.034,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,         &
               0.000,0.000/) 
!
      CONDENSABLE_CHOICE=(/0,0,1,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
!
      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,      &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0/)
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY

! .. below are the 60 indices for tropospheric chemistry species (ST)
      NO      = 0 ! not included
      NO1D    = 0 ! not included
      NO3     = 0 ! not included
      NNO     = 0 ! not included
      NNO3    = 0 ! not included
      NNO2    = 0 ! not included
      NN2O5   = 0 ! not included
      NHNO4   = 0 ! not included
      NHNO3   = 0 ! not included
      NOH     = 0 ! not included
      NHO2    = 0 ! not included
      NH2O2   = 0 ! not included
      NCH4    = 0 ! not included
      NCO     = 0 ! not included
      NCH2O   = 0 ! not included
      NMEOO   = 0 ! not included
      NH2O    = 0 ! not included
      NMHP    = 0 ! not included
      NHONO   = 0 ! not included
      NC2H6   = 0 ! not included
      NETOO   = 0 ! not included
      NETOOH  = 0 ! not included
      NMECHO  = 0 ! not included
      NMECO3  = 0 ! not included
      NPAN    = 0 ! not included
      NC3H8   = 0 ! not included
      NPNOO   = 0 ! not included
      NPIOO   = 0 ! not included
      NPNOOH  = 0 ! not included
      NPIOOH  = 0 ! not included
      NETCHO  = 0 ! not included
      NETCO3  = 0 ! not included
      NME2CO  = 0 ! not included
      NMECOO  = 0 ! not included
      NMECOH  = 0 ! not included
      NPPAN   = 0 ! not included
      NMENO3  = 0 ! not included
      NOS     = 0 ! not included
      NO1DS   = 0 ! not included
      NO3S    = 0 ! not included
      NNOXS   = 0 ! not included
      NHNO3S  = 0 ! not included
      NNOYS   = 0 ! not included
      NISOP   = 0 ! not included
      NC2H4   = 0 ! not included
      NC2H2   = 0 ! not included
      NISO2   = 0 ! not included
      NISOOH  = 0 ! not included
      NISON   = 0 ! not included
      NMACR   = 0 ! not included
      NMACRO2 = 0 ! not included
      NMACROOH= 0 ! not included
      NMPAN   = 0 ! not included
      NHACET  = 0 ! not included
      NMGLY   = 0 ! not included
      NNALD   = 0 ! not included
      NHCOOH  = 0 ! not included
      NMECO3H = 0 ! not included
      NMECO2H = 0 ! not included
      NMEOH   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
      NSOTWO  = 1
      NMeSMe  = 2
      NH2SO4  = 3
      NDMSO   = 4
      NMSA    = 5
      NCS2    = 6
      NH2S    = 7
      NCOS    = 8
      NMONOTER= 0 ! not included
      NSEC_ORG= 0 ! not included
      NH2O2F  = 9
      NO3F    = 0
      NOHF    =10
      NNO3F   =11
      NQ3D    =12
      NPT     =13
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 20 gas phase budget quantity indices for sv1
      NDMSEMOC  = 1
      NDMSTEND  = 2
      NSO2EMAN  = 3
      NSO2EMBM  = 4
      NSO2EMVL  = 5
      NSO2TEND  = 6
      NSO2DDEP  = 7
      NSO2WDEP  = 8
      NH2SO4TEND= 9
      NH2SO4DDEP=10
      NCOSEMAN  =11
      NCOSEMOC  =12
      NCOSTEND  =13
      NCS2EMAN  =14
      NCS2EMOC  =15
      NCS2TEND  =16
      NDMSOTEND =17
      NDMSODDEP =18
      NMSATEND  =19
      NMSADDEP  =20
      NTERP_EM  = 0 ! not included
      NTERP_TEND= 0 ! not included
      NTERP_DDEP= 0 ! not included
      NSORG_TEND= 0 ! not included
      NSORG_DDEP= 0 ! not included
      NSORG_WDEP= 0 ! not included
!
! REACTION INDICES for gas phase chemistry fluxes
      IOHDMS1 =  1
      IOHDMS2 =  2
      INO3DMS =  4
      IDMSOOH1=  0
      IDMSOOH2=  3
      ICS2OH  =  5
      IH2SOH  =  6
      ICOSOH  =  7
 
      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:UKCA_INDICES_SV1',     &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_SV1

! ######################################################################
      SUBROUTINE UKCA_INDICES_SV1_COUPLED

      IMPLICIT NONE
!
!-----------------------------------------------------------

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_SV1_COUPLED',zhook_in,zhook_handle)

! Main array lengths and switches
!
      NCHEMG=48          ! # of gas phase chemistry tracers
      ICHEM=1            ! 1/0 = gas phase chemistry tracers on/off
      NOFFOX=24          ! # of offline oxidant species
      NBUDCHEM=20        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
      GASBUDGET=0        ! 1/0 = gas phase chemistry fluxes on/off
      NGASBUDGET=8       ! # of gas phase chemistry fluxes
!
      NADVG=2+NCHEMG     ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX ! total # of gas phase species
!
!   sv1_coupled: NTRAG=74, NADVG=50
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
      MOX     = 1
      MNOX    = 2
      MN2O5   = 3
      MHNO4   = 4
      MHNO3   = 5
      MH2O2   = 6
      MCH4    = 7
      MCO     = 8
      MCH2O   = 9
      MMHP    =10
      MHONO   =11
      MC2H6   =12
      METOOH  =13
      MMECHO  =14
      MPAN    =15
      MC3H8   =16
      MPNOOH  =17
      MPIOOH  =18
      METCHO  =19
      MME2CO  =20
      MMECOH  =21
      MPPAN   =22
      MMENO3  =23
      MOXS    =24
      MNOYS   =25
      MISOP   =26
      MC2H4   =27
      MC2H2   =28
      MISOOH  =29
      MISON   =30
      MMACR   =31
      MMACROOH=32
      MMPAN   =33
      MHACET  =34
      MMGLY   =35
      MNALD   =36
      MHCOOH  =37
      MMECO3H =38
      MMECO2H =39
      MMEOH   =40
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
      MSOTWO  =41
      MMeSMe  =42
      MH2SO4  =43
      MDMSO   =44
      MMSA    =45
      MCS2    =46
      MH2S    =47
      MCOS    =48
      MMONOTER= 0 ! not included
      MSEC_ORG= 0 ! not included
      MH2O2F  = 0 ! not included
      MQ3D    =49
      MPT     =50
!
! .. molar masses (kg/mol) for gases for sv1_coupled
! .. (40 tropospheric gases then 8 sulphur species)
!     Dummy values used for OX =1, NOX  =2, MHP  =10,ETOOH=13, MECHO=14
!     species & no          PAN=15 PNOOH=17,PIOOH=18,ETCHO=19, ME2CO=20
!                           MECOH=21, PPAN=22, MENO3=23, OXS=24
!                           NOYS=25, ISOP=26
      MM_GAS=(/0.048,0.030,0.108,0.079,0.053,0.034,0.016,0.013,         &
               0.030,0.033,0.047,0.030,0.056,0.045,0.088,0.044,         &
               0.045,0.067,0.056,0.087,0.067,0.051,0.034,0.013,         &
               0.018,0.045,0.028,0.026,0.015,0.015,0.015,0.015,         &
               0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,         &
               0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,         &
               0.000,0.000/)

      CONDENSABLE_CHOICE=(/0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,0,0,0,0,0,0,                     &
                           0,0,0,0,0,0,1,0,0,0,0,0,0,0/)
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
!
      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
              0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,      &
              0.0,0.0/)
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY

! .. below are the 60 indices for tropospheric chemistry species (ST)
      NO      = 1
      NO1D    = 2
      NO3     = 3
      NNO     = 4
      NNO3    = 5
      NNO2    = 6
      NN2O5   = 7
      NHNO4   = 8
      NHNO3   = 9
      NOH     =10
      NHO2    =11
      NH2O2   =12
      NCH4    =13
      NCO     =14
      NCH2O   =15
      NMEOO   =16
      NH2O    =17
      NMHP    =18
      NHONO   =19
      NC2H6   =20
      NETOO   =21
      NETOOH  =22
      NMECHO  =23
      NMECO3  =24
      NPAN    =25
      NC3H8   =26
      NPNOO   =27
      NPIOO   =28
      NPNOOH  =29
      NPIOOH  =30
      NETCHO  =31
      NETCO3  =32
      NME2CO  =33
      NMECOO  =34
      NMECOH  =35
      NPPAN   =36
      NMENO3  =37
      NOS     =38
      NO1DS   =39
      NO3S    =40
      NNOXS   =41
      NHNO3S  =42
      NNOYS   =43
      NISOP   =44
      NC2H4   =45
      NC2H2   =46
      NISO2   =47
      NISOOH  =48
      NISON   =49
      NMACR   =50
      NMACRO2 =51
      NMACROOH=52
      NMPAN   =53
      NHACET  =54
      NMGLY   =55
      NNALD   =56
      NHCOOH  =57
      NMECO3H =58
      NMECO2H =59
      NMEOH   =60
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
      NSOTWO  =61
      NMeSMe  =62
      NH2SO4  =63
      NDMSO   =64
      NMSA    =65
      NCS2    =66
      NH2S    =67
      NCOS    =68
      NMONOTER= 0 ! not included
      NSEC_ORG= 0 ! not included
      NH2O2F  =74
      NO3F    =71
      NOHF    =72
      NNO3F   =73
      NQ3D    =69
      NPT     =70
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 20 gas phase budget quantity indices for sv1
      NDMSEMOC  = 1
      NDMSTEND  = 2
      NSO2EMAN  = 3
      NSO2EMBM  = 4
      NSO2EMVL  = 5
      NSO2TEND  = 6
      NSO2DDEP  = 7
      NSO2WDEP  = 8
      NH2SO4TEND= 9
      NH2SO4DDEP=10
      NCOSEMAN  =11
      NCOSEMOC  =12
      NCOSTEND  =13
      NCS2EMAN  =14
      NCS2EMOC  =15
      NCS2TEND  =16
      NDMSOTEND =17
      NDMSODDEP =18
      NMSATEND  =19
      NMSADDEP  =20
      NTERP_EM  = 0 ! not included
      NTERP_TEND= 0 ! not included
      NTERP_DDEP= 0 ! not included
      NSORG_TEND= 0 ! not included
      NSORG_DDEP= 0 ! not included
      NSORG_WDEP= 0 ! not included

! REACTION INDICES for gas phase chemistry fluxes
      IOHDMS1 =123
      IOHDMS2 =124
      INO3DMS =127
      IDMSOOH1=125
      IDMSOOH2=126
      ICS2OH  =128
      IH2SOH  =129
      ICOSOH  =130

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_SV1_COUPLED',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_SV1_COUPLED

! ######################################################################
      SUBROUTINE UKCA_INDICES_TRAQU38

      IMPLICIT NONE

!-----------------------------------------------------------

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_TRAQU38',zhook_in,zhook_handle)

! Main array lengths and switches

      TRAQU=1        ! 1/0 = to include extra diagnostics in ST
      NTRAQU=TRAQU*38! # of extra diagnostic fields
      BUDGET=1       ! 1/0 = to include budget terms in ST
      NBUDGET=NBUDCHEM+NBUDAER+NGASBUDGET ! total # of budget terms
!
! EXTRA DIAGNOSTIC INDICES
!
! .. below are 4 water content and 5 cloud, precip & windspd fields (TRAQU)
      NWTCNT1    = 1
      NWTCNT2    = 2
      NWTCNT3    = 3
      NWTCNT4    = 4
      NLCFRAC    = 5
      NLWC       = 6
      NRAINLS    = 7
      NRAINCS    = 8
      NWINDSP    = 9
      NDRYDP1    =10
      NDRYDP2    =11
      NDRYDP3    =12
      NDRYDP4    =13
      NDRYDP5    =14
      NDRYDP6    =15
      NDRYDP7    =16
! .. note below are partial volumes for SUSSBCOCSO_5MODE
      NPVOL11    =17
      NPVOL16    =18
      NPVOL1W    =19
      NPVOL21    =20
      NPVOL22    =21
      NPVOL23    =22
      NPVOL26    =23
      NPVOL2W    =24
      NPVOL31    =25
      NPVOL32    =26
      NPVOL33    =27
      NPVOL34    =28
      NPVOL36    =29
      NPVOL3W    =30
      NPVOL41    =31
      NPVOL42    =32
      NPVOL43    =33
      NPVOL44    =34
      NPVOL46    =35
      NPVOL4W    =36
      NPVOL52    =37
      NPVOL53    =38
!
!----------------------------------------------------------------

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_TRAQU38',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_TRAQU38

! ######################################################################
      SUBROUTINE UKCA_INDICES_TRAQU9

      IMPLICIT NONE

!-----------------------------------------------------------

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_TRAQU9',zhook_in,zhook_handle)

! Main array lengths and switches

      TRAQU=1        ! 1/0 = to include extra diagnostics in ST
      NTRAQU=TRAQU*9 ! # of extra diagnostic fields
      BUDGET=1       ! 1/0 = to include budget terms in ST
      NBUDGET=NBUDCHEM+NBUDAER+NGASBUDGET ! total # of budget terms
!
! EXTRA DIAGNOSTIC INDICES
!
! .. below are 4 water content and 5 cloud, precip & windspd fields (TRAQU)
      NWTCNT1    = 1
      NWTCNT2    = 2
      NWTCNT3    = 3
      NWTCNT4    = 4
      NLCFRAC    = 5
      NLWC       = 6
      NRAINLS    = 7
      NRAINCS    = 8
      NWINDSP    = 9
      NDRYDP1    = 0 ! not included
      NDRYDP2    = 0 ! not included
      NDRYDP3    = 0 ! not included
      NDRYDP4    = 0 ! not included
      NDRYDP5    = 0 ! not included
      NDRYDP6    = 0 ! not included
      NDRYDP7    = 0 ! not included
! .. note below are partial volumes for SUSSBCOCSO_5MODE
      NPVOL11    = 0 ! not included
      NPVOL16    = 0 ! not included
      NPVOL1W    = 0 ! not included
      NPVOL21    = 0 ! not included
      NPVOL22    = 0 ! not included
      NPVOL23    = 0 ! not included
      NPVOL26    = 0 ! not included
      NPVOL2W    = 0 ! not included
      NPVOL31    = 0 ! not included
      NPVOL32    = 0 ! not included
      NPVOL33    = 0 ! not included
      NPVOL34    = 0 ! not included
      NPVOL36    = 0 ! not included
      NPVOL3W    = 0 ! not included
      NPVOL41    = 0 ! not included
      NPVOL42    = 0 ! not included
      NPVOL43    = 0 ! not included
      NPVOL44    = 0 ! not included
      NPVOL46    = 0 ! not included
      NPVOL4W    = 0 ! not included
      NPVOL52    = 0 ! not included
      NPVOL53    = 0 ! not included

!
!----------------------------------------------------------------

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_TRAQU9',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_TRAQU9

! ######################################################################
      SUBROUTINE UKCA_INDICES_SUSSBCOCDU_7MODE

      IMPLICIT NONE
!---------------------------------------------------------------

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_SUSSBCOCDU_7MODE',zhook_in,zhook_handle)

! Main array lengths and switches

      NTRAER=26          ! # of aerosol advected tracers
      NBUDAER=138        ! # of aerosol budget fields
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  26    16    164     9 = 42+164+ 9=  215 (orgv1 ,traqu9 )
! ..  SUSSBCOC  26    16    164    38 = 42+164+38=  244 (orgv1 ,traqu38)
! ..  SUSSBCOC  26    76    164     9 =102+164+ 9=  275 (orgv1c,traqu9 )
! ..  SUSSBCOC  26    76    164    38 =102+164+38=  304 (orgv1c,traqu38)
!                        (138+26)
!
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 26    13  =  39 (orgv1 )
! ..  SUSSBCOC 26    52  =  78 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. 138 aerosol budget indices for SUSSBCOCDU_7MODE [SO in OC]
!
! .. redo these to follow Dominicks approach as categorized under
!
! .. BUD_PRIM, BUD_DDEP, BUD_NUSC, BUD_IMSC which have NBOX,NMODES,NCP
! ..
! .. then also have:
! ..
! .. BUD_CLPR --- NBOX,NMODES,2 (production of sulfate by H2O2, O3)
! .. BUD_COND --- NBOX,NMODES,2 (conden to each modes by H2SO4,SEC_ORG)
! .. BUD_NUCL --- NBOX,2 (BHN and BLN)
! .. BUD_COAG --- NBOX,NMODES,NMODES,NCP
! .. BUD_AGED --- NBOX,3,NCP,2 (ageing of 3 ins modes by H2SO4,SEC_ORG)
! .. BUD_MERG --- NBOX,NMODES,NCP
! .. BUD_PROC --- NBOX,NCP (processing of Aitsol mode to accsol mode)
! ..
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 0 ! BC only emitted to insoluble
      NMASPRIMBCAITINS= 6
      NMASPRIMOCAITSOL= 7
      NMASPRIMOCAITINS= 8
!
      NMASDDEPSUNUCSOL= 9
      NMASDDEPSUAITSOL=10
      NMASDDEPSUACCSOL=11
      NMASDDEPSUCORSOL=12
      NMASDDEPSSACCSOL=13
      NMASDDEPSSCORSOL=14
      NMASDDEPBCAITSOL=15
      NMASDDEPBCACCSOL=16
      NMASDDEPBCCORSOL=17
      NMASDDEPBCAITINS=18
      NMASDDEPOCNUCSOL=19
      NMASDDEPOCAITSOL=20
      NMASDDEPOCACCSOL=21
      NMASDDEPOCCORSOL=22
      NMASDDEPOCAITINS=23
      NMASDDEPSONUCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOAITSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOACCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASNUSCSUNUCSOL=24
      NMASNUSCSUAITSOL=25
      NMASNUSCSUACCSOL=26
      NMASNUSCSUCORSOL=27
      NMASNUSCSSACCSOL=28
      NMASNUSCSSCORSOL=29
      NMASNUSCBCAITSOL=30
      NMASNUSCBCACCSOL=31
      NMASNUSCBCCORSOL=32
      NMASNUSCBCAITINS=33
      NMASNUSCOCNUCSOL=34
      NMASNUSCOCAITSOL=35
      NMASNUSCOCACCSOL=36
      NMASNUSCOCCORSOL=37
      NMASNUSCOCAITINS=38
      NMASNUSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASIMSCSUNUCSOL=39
      NMASIMSCSUAITSOL=40
      NMASIMSCSUACCSOL=41
      NMASIMSCSUCORSOL=42
      NMASIMSCSSACCSOL=43
      NMASIMSCSSCORSOL=44
      NMASIMSCBCAITSOL=45
      NMASIMSCBCACCSOL=46
      NMASIMSCBCCORSOL=47
      NMASIMSCBCAITINS=48
      NMASIMSCOCNUCSOL=49
      NMASIMSCOCAITSOL=50
      NMASIMSCOCACCSOL=51
      NMASIMSCOCCORSOL=52
      NMASIMSCOCAITINS=53
      NMASIMSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASCLPRSUAITSOL1=54
      NMASCLPRSUACCSOL1=55
      NMASCLPRSUCORSOL1=56
      NMASCLPRSUAITSOL2=57
      NMASCLPRSUACCSOL2=58
      NMASCLPRSUCORSOL2=59
!
      NMASCONDSUNUCSOL=60
      NMASCONDSUAITSOL=61
      NMASCONDSUACCSOL=62
      NMASCONDSUCORSOL=63
      NMASCONDSUAITINS=64
      NMASNUCLSUNUCSOL=65
      NMASCONDOCNUCSOL=66
      NMASCONDOCAITSOL=67
      NMASCONDOCACCSOL=68
      NMASCONDOCCORSOL=69
      NMASCONDOCAITINS=70
      NMASCONDSONUCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITSOL= 0 ! SO stored in OC cpt
      NMASCONDSOACCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOCORSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITINS= 0 ! SO stored in OC cpt
!
      NMASCOAGSUINTR12=71
      NMASCOAGSUINTR13=72
      NMASCOAGSUINTR14=73
      NMASCOAGSUINTR15=74
      NMASCOAGOCINTR12=75
      NMASCOAGOCINTR13=76
      NMASCOAGOCINTR14=77
      NMASCOAGOCINTR15=78
      NMASCOAGSOINTR12= 0 ! stored in NMASCOAGOCINTR12
      NMASCOAGSOINTR13= 0 ! stored in NMASCOAGOCINTR13
      NMASCOAGSOINTR14= 0 ! stored in NMASCOAGOCINTR14
      NMASCOAGSOINTR15= 0 ! stored in NMASCOAGOCINTR15
      NMASCOAGSUINTR23=79
      NMASCOAGBCINTR23=80
      NMASCOAGOCINTR23=81
      NMASCOAGSOINTR23= 0 ! stored in NMASCOAGOCINTR23
      NMASCOAGSUINTR24=82
      NMASCOAGBCINTR24=83
      NMASCOAGOCINTR24=84
      NMASCOAGSOINTR24= 0 ! stored in NMASCOAGOCINTR24
      NMASCOAGSUINTR34=85
      NMASCOAGBCINTR34=86
      NMASCOAGOCINTR34=87
      NMASCOAGSSINTR34=88
      NMASCOAGSOINTR34= 0 ! stored in NMASCOAGOCINTR34
!
      NMASCOAGBCINTR53=89
      NMASCOAGOCINTR53=90
      NMASCOAGBCINTR54=91
      NMASCOAGOCINTR54=92
!
      NMASAGEDSUINTR52=93
      NMASAGEDBCINTR52=94
      NMASAGEDOCINTR52=95
      NMASAGEDSOINTR52= 0 ! stored in NMASAGEDOCINTR52
!
      NMASMERGSUINTR12=96
      NMASMERGOCINTR12=97
      NMASMERGSOINTR12= 0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR23=98
      NMASMERGBCINTR23=99
      NMASMERGOCINTR23=100
      NMASMERGSOINTR23=  0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR34=101
      NMASMERGSSINTR34=102
      NMASMERGBCINTR34=103
      NMASMERGOCINTR34=104
      NMASMERGSOINTR34=  0 ! stored in NMASMERGOCINTR34
      NMASPROCSUINTR23=105
      NMASPROCBCINTR23=106
      NMASPROCOCINTR23=107
      NMASPROCSOINTR23=  0 ! stored in NMASPROCOCINTR23
!
! .. below are new ones for dust & modes 6/7 to be integrated
      NMASPRIMDUACCSOL=  0 ! DU only emitted to modes 6/7 here
      NMASPRIMDUCORSOL=  0 ! DU only emitted to modes 6/7 here
      NMASPRIMDUACCINS=108
      NMASPRIMDUCORINS=109
      NMASDDEPDUACCSOL=110
      NMASDDEPDUCORSOL=111
      NMASDDEPDUACCINS=112
      NMASDDEPDUCORINS=113
      NMASNUSCDUACCSOL=114
      NMASNUSCDUCORSOL=115
      NMASNUSCDUACCINS=116
      NMASNUSCDUCORINS=117
      NMASIMSCDUACCSOL=118
      NMASIMSCDUCORSOL=119
      NMASIMSCDUACCINS=120
      NMASIMSCDUCORINS=121
      NMASCONDSUACCINS=122
      NMASCONDSUCORINS=123
      NMASCONDOCACCINS=124
      NMASCONDOCCORINS=125
      NMASCONDSOACCINS=  0 ! secondary organic in OC cpt in this setup
      NMASCONDSOCORINS=  0 ! secondary organic in OC cpt in this setup
      NMASCOAGSUINTR16=126
      NMASCOAGSUINTR17=127
      NMASCOAGOCINTR16=128
      NMASCOAGOCINTR17=129
      NMASCOAGSOINTR16=  0 ! secondary organic in OC cpt in this setup
      NMASCOAGSOINTR17=  0 ! secondary organic in OC cpt in this setup
      NMASCOAGDUINTR34=130
      NMASCOAGDUINTR64=131
      NMASAGEDSUINTR63=132
      NMASAGEDDUINTR63=133
      NMASAGEDOCINTR63=134
      NMASAGEDSOINTR63=  0 ! secondary organic in OC cpt in this setup
      NMASAGEDSUINTR74=135
      NMASAGEDDUINTR74=136
      NMASAGEDOCINTR74=137
      NMASAGEDSOINTR74=  0 ! secondary organic in OC cpt in this setup
      NMASMERGDUINTR34=138
!
!----------------------------------------------------------------

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_SUSSBCOCDU_7MODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_SUSSBCOCDU_7MODE

! ######################################################################
      SUBROUTINE UKCA_INDICES_SUSSBCOCDU_4MODE

      IMPLICIT NONE
!---------------------------------------------------------------

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_SUSSBCOCDU_4MODE',zhook_in,zhook_handle)

! MAIN ARRAY LENGTHS AND SWITCHES
!
      NTRAER=19          ! # of aerosol advected tracers
      NBUDAER=99         ! # of aerosol budget fields
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  19    16    125     9 = 35+125+ 9=  169 (orgv1 ,traqu9 )
! ..  SUSSBCOC  19    16    125    38 = 35+125+38=  198 (orgv1 ,traqu38)
! ..  SUSSBCOC  19    76    125     9 = 95+125+ 9=  229 (orgv1c,traqu9 )
! ..  SUSSBCOC  19    76    125    38 = 95+125+38=  258 (orgv1c,traqu38)
!                        (99+26)
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 19    13  =  32 (orgv1 )
! ..  SUSSBCOC 19    52  =  71 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! ..  99 aerosol budget indices for SUSSBCOCDU_4mode [SO in OC]
!
! .. redo these to follow Dominicks approach as categorized under
!
! .. BUD_PRIM, BUD_DDEP, BUD_NUSC, BUD_IMSC which have NBOX,NMODES,NCP
! ..
! .. then also have:
! ..
! .. BUD_CLPR --- NBOX,NMODES,2 (production of sulfate by H2O2, O3)
! .. BUD_COND --- NBOX,NMODES,2 (conden to each modes by H2SO4,SEC_ORG)
! .. BUD_NUCL --- NBOX,2 (BHN and BLN)
! .. BUD_COAG --- NBOX,NMODES,NMODES,NCP
! .. BUD_AGED --- NBOX,3,NCP,2 (ageing of 3 ins modes by H2SO4,SEC_ORG)
! .. BUD_MERG --- NBOX,NMODES,NCP
! .. BUD_PROC --- NBOX,NCP (processing of Aitsol mode to accsol mode)
! ..
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 6
      NMASPRIMBCAITINS= 0 ! BC only emitted to soluble
      NMASPRIMOCAITSOL= 7
      NMASPRIMOCAITINS= 0 ! OC only emitted to soluble
!
      NMASDDEPSUNUCSOL= 8
      NMASDDEPSUAITSOL= 9
      NMASDDEPSUACCSOL=10
      NMASDDEPSUCORSOL=11
      NMASDDEPSSACCSOL=12
      NMASDDEPSSCORSOL=13
      NMASDDEPBCAITSOL=14
      NMASDDEPBCACCSOL=15
      NMASDDEPBCCORSOL=16
      NMASDDEPBCAITINS= 0 ! BC only present in soluble
      NMASDDEPOCNUCSOL=17
      NMASDDEPOCAITSOL=18
      NMASDDEPOCACCSOL=19
      NMASDDEPOCCORSOL=20
      NMASDDEPOCAITINS= 0 ! OC only present in soluble
      NMASDDEPSONUCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOAITSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOACCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASNUSCSUNUCSOL=21
      NMASNUSCSUAITSOL=22
      NMASNUSCSUACCSOL=23
      NMASNUSCSUCORSOL=24
      NMASNUSCSSACCSOL=25
      NMASNUSCSSCORSOL=26
      NMASNUSCBCAITSOL=27
      NMASNUSCBCACCSOL=28
      NMASNUSCBCCORSOL=29
      NMASNUSCBCAITINS= 0 ! BC only present in soluble
      NMASNUSCOCNUCSOL=30
      NMASNUSCOCAITSOL=31
      NMASNUSCOCACCSOL=32
      NMASNUSCOCCORSOL=33
      NMASNUSCOCAITINS= 0 ! OC only present in soluble
      NMASNUSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASIMSCSUNUCSOL=34
      NMASIMSCSUAITSOL=35
      NMASIMSCSUACCSOL=36
      NMASIMSCSUCORSOL=37
      NMASIMSCSSACCSOL=38
      NMASIMSCSSCORSOL=39
      NMASIMSCBCAITSOL=40
      NMASIMSCBCACCSOL=41
      NMASIMSCBCCORSOL=42
      NMASIMSCBCAITINS= 0 ! BC only present in soluble
      NMASIMSCOCNUCSOL=43
      NMASIMSCOCAITSOL=44
      NMASIMSCOCACCSOL=45
      NMASIMSCOCCORSOL=46
      NMASIMSCOCAITINS= 0 ! OC only present in soluble
      NMASIMSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASCLPRSUAITSOL1=47
      NMASCLPRSUACCSOL1=48
      NMASCLPRSUCORSOL1=49
      NMASCLPRSUAITSOL2=50
      NMASCLPRSUACCSOL2=51
      NMASCLPRSUCORSOL2=52
!
      NMASCONDSUNUCSOL=53
      NMASCONDSUAITSOL=54
      NMASCONDSUACCSOL=55
      NMASCONDSUCORSOL=56
      NMASCONDSUAITINS= 0 ! only soluble modes
      NMASNUCLSUNUCSOL=57
      NMASCONDOCNUCSOL=58
      NMASCONDOCAITSOL=59
      NMASCONDOCACCSOL=60
      NMASCONDOCCORSOL=61
      NMASCONDOCAITINS= 0 ! only soluble modes
      NMASCONDSONUCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITSOL= 0 ! SO stored in OC cpt
      NMASCONDSOACCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOCORSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITINS= 0 ! SO stored in OC cpt
!
      NMASCOAGSUINTR12=62
      NMASCOAGSUINTR13=63
      NMASCOAGSUINTR14=64
      NMASCOAGSUINTR15= 0 ! only soluble modes
      NMASCOAGOCINTR12=65
      NMASCOAGOCINTR13=66
      NMASCOAGOCINTR14=67
      NMASCOAGOCINTR15= 0 ! only soluble modes
      NMASCOAGSOINTR12= 0 ! stored in NMASCOAGOCINTR12
      NMASCOAGSOINTR13= 0 ! stored in NMASCOAGOCINTR13
      NMASCOAGSOINTR14= 0 ! stored in NMASCOAGOCINTR14
      NMASCOAGSOINTR15= 0 ! stored in NMASCOAGOCINTR15
      NMASCOAGSUINTR23=68
      NMASCOAGBCINTR23=69
      NMASCOAGOCINTR23=70
      NMASCOAGSOINTR23= 0 ! stored in NMASCOAGOCINTR23
      NMASCOAGSUINTR24=71
      NMASCOAGBCINTR24=72
      NMASCOAGOCINTR24=73
      NMASCOAGSOINTR24= 0 ! stored in NMASCOAGOCINTR24
      NMASCOAGSUINTR34=74
      NMASCOAGBCINTR34=75
      NMASCOAGOCINTR34=76
      NMASCOAGSSINTR34=77
      NMASCOAGSOINTR34= 0 ! stored in NMASCOAGOCINTR34
!
      NMASCOAGBCINTR53= 0 ! only soluble modes
      NMASCOAGOCINTR53= 0 ! only soluble modes
      NMASCOAGBCINTR54= 0 ! only soluble modes
      NMASCOAGOCINTR54= 0 ! only soluble modes
!
      NMASAGEDSUINTR52= 0 ! only soluble modes
      NMASAGEDBCINTR52= 0 ! only soluble modes
      NMASAGEDOCINTR52= 0 ! only soluble modes
      NMASAGEDSOINTR52= 0 ! only soluble modes
!
      NMASMERGSUINTR12=78
      NMASMERGOCINTR12=79
      NMASMERGSOINTR12= 0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR23=80
      NMASMERGBCINTR23=81
      NMASMERGOCINTR23=82
      NMASMERGSOINTR23= 0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR34=83
      NMASMERGSSINTR34=84
      NMASMERGBCINTR34=85
      NMASMERGOCINTR34=86
      NMASMERGSOINTR34= 0 ! stored in NMASMERGOCINTR34
      NMASPROCSUINTR23=87
      NMASPROCBCINTR23=88
      NMASPROCOCINTR23=89
      NMASPROCSOINTR23= 0 ! stored in NMASPROCOCINTR23
!
! .. below are new ones for dust & modes 6/7 to be integrated
      NMASPRIMDUACCSOL= 90
      NMASPRIMDUCORSOL= 91
      NMASPRIMDUACCINS=  0 ! no acc-ins cor-ins in this setup
      NMASPRIMDUCORINS=  0 ! no acc-ins cor-ins in this setup
      NMASDDEPDUACCSOL= 92
      NMASDDEPDUCORSOL= 93
      NMASDDEPDUACCINS=  0 ! no acc-ins cor-ins in this setup
      NMASDDEPDUCORINS=  0 ! no acc-ins cor-ins in this setup
      NMASNUSCDUACCSOL= 94
      NMASNUSCDUCORSOL= 95
      NMASNUSCDUACCINS=  0 ! no acc-ins cor-ins in this setup
      NMASNUSCDUCORINS=  0 ! no acc-ins cor-ins in this setup
      NMASIMSCDUACCSOL= 96
      NMASIMSCDUCORSOL= 97
      NMASIMSCDUACCINS=  0 ! no acc-ins cor-ins in this setup
      NMASIMSCDUCORINS=  0 ! no acc-ins cor-ins in this setup
      NMASCONDSUACCINS=  0 ! no acc-ins cor-ins in this setup
      NMASCONDSUCORINS=  0 ! no acc-ins cor-ins in this setup
      NMASCONDOCACCINS=  0 ! no acc-ins cor-ins in this setup
      NMASCONDOCCORINS=  0 ! no acc-ins cor-ins in this setup
      NMASCONDSOACCINS=  0 ! no acc-ins cor-ins in this setup
      NMASCONDSOCORINS=  0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR16=  0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR17=  0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR16=  0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR17=  0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR16=  0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR17=  0 ! no acc-ins cor-ins in this setup
      NMASCOAGDUINTR34= 98
      NMASCOAGDUINTR64=  0 ! no acc-ins cor-ins in this setup
      NMASAGEDSUINTR63=  0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR63=  0 ! no acc-ins cor-ins in this setup
      NMASAGEDOCINTR63=  0 ! no acc-ins cor-ins in this setup
      NMASAGEDSOINTR63=  0 ! no acc-ins cor-ins in this setup
      NMASAGEDSUINTR74=  0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR74=  0 ! no acc-ins cor-ins in this setup
      NMASAGEDOCINTR74=  0 ! no acc-ins cor-ins in this setup
      NMASAGEDSOINTR74=  0 ! no acc-ins cor-ins in this setup
      NMASMERGDUINTR34= 99
!
!----------------------------------------------------------------

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
         'UKCA_INDICES_SUSSBCOCDU_4MODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_SUSSBCOCDU_4MODE

! ######################################################################
      SUBROUTINE UKCA_INDICES_SUSSBCOC_5MODE

      IMPLICIT NONE
!---------------------------------------------------------------
!
! Main array lengths and switches
!

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
        'UKCA_INDICES_SUSSBCOC_5MODE',zhook_in,zhook_handle)

      NTRAER=20          ! # of aerosol advected tracers
      NBUDAER=107        ! # of aerosol budget fields
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  20    16    133     9 = 36+133+ 9=  178 (orgv1 ,traqu9 )
! ..  SUSSBCOC  20    16    133    38 = 36+133+38=  207 (orgv1 ,traqu38)
! ..  SUSSBCOC  20    76    133     9 = 96+133+ 9=  238 (orgv1c,traqu9 )
! ..  SUSSBCOC  20    76    133    38 = 96+133+38=  267 (orgv1c,traqu38)
!                         (107+26)
!
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 20    13  =  33 (orgv1 )
! ..  SUSSBCOC 20    52  =  72 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are 107 aerosol budget indices for SUSSBCOC [SO in OC]
!
! .. redo these to follow Dominicks approach as categorized under
!
! .. BUD_PRIM, BUD_DDEP, BUD_NUSC, BUD_IMSC which have NBOX,NMODES,NCP
! ..
! .. then also have:
! ..
! .. BUD_CLPR --- NBOX,NMODES,2 (production of sulfate by H2O2, O3)
! .. BUD_COND --- NBOX,NMODES,2 (conden to each modes by H2SO4,SEC_ORG)
! .. BUD_NUCL --- NBOX,2 (BHN and BLN)
! .. BUD_COAG --- NBOX,NMODES,NMODES,NCP
! .. BUD_AGED --- NBOX,3,NCP,2 (ageing of 3 ins modes by H2SO4,SEC_ORG)
! .. BUD_MERG --- NBOX,NMODES,NCP
! .. BUD_PROC --- NBOX,NCP (processing of Aitsol mode to accsol mode)
! ..
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 0 ! BC only emitted to insoluble
      NMASPRIMBCAITINS= 6
      NMASPRIMOCAITSOL= 7
      NMASPRIMOCAITINS= 8
!
      NMASDDEPSUNUCSOL= 9
      NMASDDEPSUAITSOL=10
      NMASDDEPSUACCSOL=11
      NMASDDEPSUCORSOL=12
      NMASDDEPSSACCSOL=13
      NMASDDEPSSCORSOL=14
      NMASDDEPBCAITSOL=15
      NMASDDEPBCACCSOL=16
      NMASDDEPBCCORSOL=17
      NMASDDEPBCAITINS=18
      NMASDDEPOCNUCSOL=19
      NMASDDEPOCAITSOL=20
      NMASDDEPOCACCSOL=21
      NMASDDEPOCCORSOL=22
      NMASDDEPOCAITINS=23
      NMASDDEPSONUCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOAITSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOACCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASNUSCSUNUCSOL=24
      NMASNUSCSUAITSOL=25
      NMASNUSCSUACCSOL=26
      NMASNUSCSUCORSOL=27
      NMASNUSCSSACCSOL=28
      NMASNUSCSSCORSOL=29
      NMASNUSCBCAITSOL=30
      NMASNUSCBCACCSOL=31
      NMASNUSCBCCORSOL=32
      NMASNUSCBCAITINS=33
      NMASNUSCOCNUCSOL=34
      NMASNUSCOCAITSOL=35
      NMASNUSCOCACCSOL=36
      NMASNUSCOCCORSOL=37
      NMASNUSCOCAITINS=38
      NMASNUSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASIMSCSUNUCSOL=39
      NMASIMSCSUAITSOL=40
      NMASIMSCSUACCSOL=41
      NMASIMSCSUCORSOL=42
      NMASIMSCSSACCSOL=43
      NMASIMSCSSCORSOL=44
      NMASIMSCBCAITSOL=45
      NMASIMSCBCACCSOL=46
      NMASIMSCBCCORSOL=47
      NMASIMSCBCAITINS=48
      NMASIMSCOCNUCSOL=49
      NMASIMSCOCAITSOL=50
      NMASIMSCOCACCSOL=51
      NMASIMSCOCCORSOL=52
      NMASIMSCOCAITINS=53
      NMASIMSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASCLPRSUAITSOL1=54
      NMASCLPRSUACCSOL1=55
      NMASCLPRSUCORSOL1=56
      NMASCLPRSUAITSOL2=57
      NMASCLPRSUACCSOL2=58
      NMASCLPRSUCORSOL2=59
!
      NMASCONDSUNUCSOL=60
      NMASCONDSUAITSOL=61
      NMASCONDSUACCSOL=62
      NMASCONDSUCORSOL=63
      NMASCONDSUAITINS=64
      NMASNUCLSUNUCSOL=65
      NMASCONDOCNUCSOL=66
      NMASCONDOCAITSOL=67
      NMASCONDOCACCSOL=68
      NMASCONDOCCORSOL=69
      NMASCONDOCAITINS=70
      NMASCONDSONUCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITSOL= 0 ! SO stored in OC cpt
      NMASCONDSOACCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOCORSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITINS= 0 ! SO stored in OC cpt
!
      NMASCOAGSUINTR12=71
      NMASCOAGSUINTR13=72
      NMASCOAGSUINTR14=73
      NMASCOAGSUINTR15=74
      NMASCOAGOCINTR12=75
      NMASCOAGOCINTR13=76
      NMASCOAGOCINTR14=77
      NMASCOAGOCINTR15=78
      NMASCOAGSOINTR12= 0 ! stored in NMASCOAGOCINTR12
      NMASCOAGSOINTR13= 0 ! stored in NMASCOAGOCINTR13
      NMASCOAGSOINTR14= 0 ! stored in NMASCOAGOCINTR14
      NMASCOAGSOINTR15= 0 ! stored in NMASCOAGOCINTR15
      NMASCOAGSUINTR23=79
      NMASCOAGBCINTR23=80
      NMASCOAGOCINTR23=81
      NMASCOAGSOINTR23= 0 ! stored in NMASCOAGOCINTR23
      NMASCOAGSUINTR24=82
      NMASCOAGBCINTR24=83
      NMASCOAGOCINTR24=84
      NMASCOAGSOINTR24= 0 ! stored in NMASCOAGOCINTR24
      NMASCOAGSUINTR34=85
      NMASCOAGBCINTR34=86
      NMASCOAGOCINTR34=87
      NMASCOAGSSINTR34=88
      NMASCOAGSOINTR34= 0 ! stored in NMASCOAGOCINTR34
!
      NMASCOAGBCINTR53=89
      NMASCOAGOCINTR53=90
      NMASCOAGBCINTR54=91
      NMASCOAGOCINTR54=92
!
      NMASAGEDSUINTR52=93
      NMASAGEDBCINTR52=94
      NMASAGEDOCINTR52=95
      NMASAGEDSOINTR52= 0 ! stored in NMASAGEDOCINTR52
!
      NMASMERGSUINTR12=96
      NMASMERGOCINTR12=97
      NMASMERGSOINTR12= 0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR23=98
      NMASMERGBCINTR23=99
      NMASMERGOCINTR23=100
      NMASMERGSOINTR23=  0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR34=101
      NMASMERGSSINTR34=102
      NMASMERGBCINTR34=103
      NMASMERGOCINTR34=104
      NMASMERGSOINTR34=  0 ! stored in NMASMERGOCINTR34
      NMASPROCSUINTR23=105
      NMASPROCBCINTR23=106
      NMASPROCOCINTR23=107
      NMASPROCSOINTR23=  0 ! stored in NMASPROCOCINTR23
!
! .. below are new ones for dust & modes 6/7 to be integrated
      NMASPRIMDUACCSOL= 0 ! no DU in this setup
      NMASPRIMDUCORSOL= 0 ! no DU in this setup
      NMASPRIMDUACCINS= 0 ! no DU in this setup
      NMASPRIMDUCORINS= 0 ! no DU in this setup
      NMASDDEPDUACCSOL= 0 ! no DU in this setup
      NMASDDEPDUCORSOL= 0 ! no DU in this setup
      NMASDDEPDUACCINS= 0 ! no DU in this setup
      NMASDDEPDUCORINS= 0 ! no DU in this setup
      NMASNUSCDUACCSOL= 0 ! no DU in this setup
      NMASNUSCDUCORSOL= 0 ! no DU in this setup
      NMASNUSCDUACCINS= 0 ! no DU in this setup
      NMASNUSCDUCORINS= 0 ! no DU in this setup
      NMASIMSCDUACCSOL= 0 ! no DU in this setup
      NMASIMSCDUCORSOL= 0 ! no DU in this setup
      NMASIMSCDUACCINS= 0 ! no DU in this setup
      NMASIMSCDUCORINS= 0 ! no DU in this setup
      NMASCONDSUACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSUCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDOCACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDOCCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSOACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSOCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGDUINTR34= 0 ! no DU in this setup
      NMASCOAGDUINTR64= 0 ! no DU in this setup
      NMASAGEDSUINTR63= 0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR63= 0 ! no DU in this setup
      NMASAGEDOCINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSUINTR74= 0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR74= 0 ! no DU in this setup
      NMASAGEDOCINTR74= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR74= 0 ! no BC/OC/SO in this setup
      NMASMERGDUINTR34= 0 ! no DU in this setup

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
        'UKCA_INDICES_SUSSBCOC_5MODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_SUSSBCOC_5MODE

! ######################################################################
      SUBROUTINE UKCA_INDICES_SUSSBCOC_4MODE

      IMPLICIT NONE
!---------------------------------------------------------------
!
! Main array lengths and switches

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
                'UKCA_INDICES_SUSSBCOC_4MODE',zhook_in,zhook_handle)

      NTRAER=17          ! # of aerosol advected tracers
      NBUDAER=89         ! # of aerosol budget fields
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  17    16    115     9 = 33+115+ 9=  157 (orgv1 ,traqu9 )
! ..  SUSSBCOC  17    16    115    38 = 33+115+38=  186 (orgv1 ,traqu38)
! ..  SUSSBCOC  17    76    115     9 = 63+115+ 9=  217 (orgv1c,traqu9 )
! ..  SUSSBCOC  17    76    115    38 = 63+115+38=  246 (orgv1c,traqu38)
!                        (89+26)
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 17    13  =  30
! ..  SUSSBCOC 17    52  =  69
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! ..  89 aerosol budget indices for SUSSBCOC_4MODE [SO in OC]
!
! .. redo these to follow Dominicks approach as categorized under
!
! .. BUD_PRIM, BUD_DDEP, BUD_NUSC, BUD_IMSC which have NBOX,NMODES,NCP
! ..
! .. then also have:
! ..
! .. BUD_CLPR --- NBOX,NMODES,2 (production of sulfate by H2O2, O3)
! .. BUD_COND --- NBOX,NMODES,2 (conden to each modes by H2SO4,SEC_ORG)
! .. BUD_NUCL --- NBOX,2 (BHN and BLN)
! .. BUD_COAG --- NBOX,NMODES,NMODES,NCP
! .. BUD_AGED --- NBOX,3,NCP,2 (ageing of 3 ins modes by H2SO4,SEC_ORG)
! .. BUD_MERG --- NBOX,NMODES,NCP
! .. BUD_PROC --- NBOX,NCP (processing of Aitsol mode to accsol mode)
! ..
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 6
      NMASPRIMBCAITINS= 0 ! BC only emitted to soluble
      NMASPRIMOCAITSOL= 7
      NMASPRIMOCAITINS= 0 ! OC only emitted to soluble
!
      NMASDDEPSUNUCSOL= 8
      NMASDDEPSUAITSOL= 9
      NMASDDEPSUACCSOL=10
      NMASDDEPSUCORSOL=11
      NMASDDEPSSACCSOL=12
      NMASDDEPSSCORSOL=13
      NMASDDEPBCAITSOL=14
      NMASDDEPBCACCSOL=15
      NMASDDEPBCCORSOL=16
      NMASDDEPBCAITINS= 0 ! BC only present in soluble
      NMASDDEPOCNUCSOL=17
      NMASDDEPOCAITSOL=18
      NMASDDEPOCACCSOL=19
      NMASDDEPOCCORSOL=20
      NMASDDEPOCAITINS= 0 ! OC only present in soluble
      NMASDDEPSONUCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOAITSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOACCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASNUSCSUNUCSOL=21
      NMASNUSCSUAITSOL=22
      NMASNUSCSUACCSOL=23
      NMASNUSCSUCORSOL=24
      NMASNUSCSSACCSOL=25
      NMASNUSCSSCORSOL=26
      NMASNUSCBCAITSOL=27
      NMASNUSCBCACCSOL=28
      NMASNUSCBCCORSOL=29
      NMASNUSCBCAITINS= 0 ! BC only present in soluble
      NMASNUSCOCNUCSOL=30
      NMASNUSCOCAITSOL=31
      NMASNUSCOCACCSOL=32
      NMASNUSCOCCORSOL=33
      NMASNUSCOCAITINS= 0 ! OC only present in soluble
      NMASNUSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASIMSCSUNUCSOL=34
      NMASIMSCSUAITSOL=35
      NMASIMSCSUACCSOL=36
      NMASIMSCSUCORSOL=37
      NMASIMSCSSACCSOL=38
      NMASIMSCSSCORSOL=39
      NMASIMSCBCAITSOL=40
      NMASIMSCBCACCSOL=41
      NMASIMSCBCCORSOL=42
      NMASIMSCBCAITINS= 0 ! BC only present in soluble
      NMASIMSCOCNUCSOL=43
      NMASIMSCOCAITSOL=44
      NMASIMSCOCACCSOL=45
      NMASIMSCOCCORSOL=46
      NMASIMSCOCAITINS= 0 ! OC only present in soluble
      NMASIMSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASCLPRSUAITSOL1=47
      NMASCLPRSUACCSOL1=48
      NMASCLPRSUCORSOL1=49
      NMASCLPRSUAITSOL2=50
      NMASCLPRSUACCSOL2=51
      NMASCLPRSUCORSOL2=52
!
      NMASCONDSUNUCSOL=53
      NMASCONDSUAITSOL=54
      NMASCONDSUACCSOL=55
      NMASCONDSUCORSOL=56
      NMASCONDSUAITINS= 0 ! only soluble modes
      NMASNUCLSUNUCSOL=57
      NMASCONDOCNUCSOL=58
      NMASCONDOCAITSOL=59
      NMASCONDOCACCSOL=60
      NMASCONDOCCORSOL=61
      NMASCONDOCAITINS= 0 ! only soluble modes
      NMASCONDSONUCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITSOL= 0 ! SO stored in OC cpt
      NMASCONDSOACCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOCORSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITINS= 0 ! SO stored in OC cpt
!
      NMASCOAGSUINTR12=62
      NMASCOAGSUINTR13=63
      NMASCOAGSUINTR14=64
      NMASCOAGSUINTR15= 0 ! only soluble modes
      NMASCOAGOCINTR12=65
      NMASCOAGOCINTR13=66
      NMASCOAGOCINTR14=67
      NMASCOAGOCINTR15= 0 ! only soluble modes
      NMASCOAGSOINTR12= 0 ! stored in NMASCOAGOCINTR12
      NMASCOAGSOINTR13= 0 ! stored in NMASCOAGOCINTR13
      NMASCOAGSOINTR14= 0 ! stored in NMASCOAGOCINTR14
      NMASCOAGSOINTR15= 0 ! stored in NMASCOAGOCINTR15
      NMASCOAGSUINTR23=68
      NMASCOAGBCINTR23=69
      NMASCOAGOCINTR23=70
      NMASCOAGSOINTR23= 0 ! stored in NMASCOAGOCINTR23
      NMASCOAGSUINTR24=71
      NMASCOAGBCINTR24=72
      NMASCOAGOCINTR24=73
      NMASCOAGSOINTR24= 0 ! stored in NMASCOAGOCINTR24
      NMASCOAGSUINTR34=74
      NMASCOAGBCINTR34=75
      NMASCOAGOCINTR34=76
      NMASCOAGSSINTR34=77
      NMASCOAGSOINTR34= 0 ! stored in NMASCOAGOCINTR34
!
      NMASCOAGBCINTR53= 0 ! only soluble modes
      NMASCOAGOCINTR53= 0 ! only soluble modes
      NMASCOAGBCINTR54= 0 ! only soluble modes
      NMASCOAGOCINTR54= 0 ! only soluble modes
!
      NMASAGEDSUINTR52= 0 ! only soluble modes
      NMASAGEDBCINTR52= 0 ! only soluble modes
      NMASAGEDOCINTR52= 0 ! only soluble modes
      NMASAGEDSOINTR52= 0 ! only soluble modes
!
      NMASMERGSUINTR12=78
      NMASMERGOCINTR12=79
      NMASMERGSOINTR12= 0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR23=80
      NMASMERGBCINTR23=81
      NMASMERGOCINTR23=82
      NMASMERGSOINTR23= 0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR34=83
      NMASMERGSSINTR34=84
      NMASMERGBCINTR34=85
      NMASMERGOCINTR34=86
      NMASMERGSOINTR34= 0 ! stored in NMASMERGOCINTR34
      NMASPROCSUINTR23=87
      NMASPROCBCINTR23=88
      NMASPROCOCINTR23=89
      NMASPROCSOINTR23= 0 ! stored in NMASPROCOCINTR23
!
! .. below are new ones for dust & modes 6/7 to be integrated
      NMASPRIMDUACCSOL= 0 ! no DU in this setup
      NMASPRIMDUCORSOL= 0 ! no DU in this setup
      NMASPRIMDUACCINS= 0 ! no DU in this setup
      NMASPRIMDUCORINS= 0 ! no DU in this setup
      NMASDDEPDUACCSOL= 0 ! no DU in this setup
      NMASDDEPDUCORSOL= 0 ! no DU in this setup
      NMASDDEPDUACCINS= 0 ! no DU in this setup
      NMASDDEPDUCORINS= 0 ! no DU in this setup
      NMASNUSCDUACCSOL= 0 ! no DU in this setup
      NMASNUSCDUCORSOL= 0 ! no DU in this setup
      NMASNUSCDUACCINS= 0 ! no DU in this setup
      NMASNUSCDUCORINS= 0 ! no DU in this setup
      NMASIMSCDUACCSOL= 0 ! no DU in this setup
      NMASIMSCDUCORSOL= 0 ! no DU in this setup
      NMASIMSCDUACCINS= 0 ! no DU in this setup
      NMASIMSCDUCORINS= 0 ! no DU in this setup
      NMASCONDSUACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSUCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDOCACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDOCCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSOACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSOCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGDUINTR34= 0 ! no DU in this setup
      NMASCOAGDUINTR64= 0 ! no DU in this setup
      NMASAGEDSUINTR63= 0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR63= 0 ! no DU in this setup
      NMASAGEDOCINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSUINTR74= 0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR74= 0 ! no DU in this setup
      NMASAGEDOCINTR74= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR74= 0 ! no BC/OC/SO in this setup
      NMASMERGDUINTR34= 0 ! no DU in this setup

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
                'UKCA_INDICES_SUSSBCOC_4MODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_SUSSBCOC_4MODE

! ######################################################################
      SUBROUTINE UKCA_INDICES_SUSSBCOCSO_5MODE

      IMPLICIT NONE

!---------------------------------------------------------------
!
! Main array lengths and switches

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
                'UKCA_INDICES_SUSSBCOCSO_5MODE',zhook_in,zhook_handle)

      NTRAER=23          ! # of aerosol advected tracers
      NBUDAER=123        ! # of aerosol budget fields

! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  23    16    149     9 = 39+149+ 9=  197 (orgv1 ,traqu9 )
! ..  SUSSBCOC  23    16    149    38 = 39+149+38=  226 (orgv1 ,traqu38)
! ..  SUSSBCOC  23    76    149     9 = 99+149+ 9=  257 (orgv1c,traqu9 )
! ..  SUSSBCOC  23    76    149    38 = 99+149+38=  286 (orgv1c,traqu38)
!                        (123+26)
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 23    13  =  36 (orgv1 )
! ..  SUSSBCOC 23    52  =  75 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are 123 aerosol budget indices for SUSSBCOCSO [SO in SO]
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 0 ! BC only emitted to insoluble
      NMASPRIMBCAITINS= 6
      NMASPRIMOCAITSOL= 7
      NMASPRIMOCAITINS= 8
!
      NMASDDEPSUNUCSOL= 9
      NMASDDEPSUAITSOL=10
      NMASDDEPSUACCSOL=11
      NMASDDEPSUCORSOL=12
      NMASDDEPSSACCSOL=13
      NMASDDEPSSCORSOL=14
      NMASDDEPBCAITSOL=15
      NMASDDEPBCACCSOL=16
      NMASDDEPBCCORSOL=17
      NMASDDEPBCAITINS=18
      NMASDDEPOCNUCSOL= 0  ! stored in NMASDDEPSONUCSOL
      NMASDDEPOCAITSOL=19
      NMASDDEPOCACCSOL=20
      NMASDDEPOCCORSOL=21
      NMASDDEPOCAITINS=22
      NMASDDEPSONUCSOL=23
      NMASDDEPSOAITSOL=24
      NMASDDEPSOACCSOL=25
      NMASDDEPSOCORSOL=26
!
      NMASNUSCSUNUCSOL=27
      NMASNUSCSUAITSOL=28
      NMASNUSCSUACCSOL=29
      NMASNUSCSUCORSOL=30
      NMASNUSCSSACCSOL=31
      NMASNUSCSSCORSOL=32
      NMASNUSCBCAITSOL=33
      NMASNUSCBCACCSOL=34
      NMASNUSCBCCORSOL=35
      NMASNUSCBCAITINS=36
      NMASNUSCOCNUCSOL= 0  ! stored in NMASNUSCSONUCSOL
      NMASNUSCOCAITSOL=37
      NMASNUSCOCACCSOL=38
      NMASNUSCOCCORSOL=39
      NMASNUSCOCAITINS=40
      NMASNUSCSONUCSOL=41
      NMASNUSCSOAITSOL=42
      NMASNUSCSOACCSOL=43
      NMASNUSCSOCORSOL=44
!
      NMASIMSCSUNUCSOL=45
      NMASIMSCSUAITSOL=46
      NMASIMSCSUACCSOL=47
      NMASIMSCSUCORSOL=48
      NMASIMSCSSACCSOL=49
      NMASIMSCSSCORSOL=50
      NMASIMSCBCAITSOL=51
      NMASIMSCBCACCSOL=52
      NMASIMSCBCCORSOL=53
      NMASIMSCBCAITINS=54
      NMASIMSCOCNUCSOL= 0 ! stored in NMASIMSCSONUCSOL
      NMASIMSCOCAITSOL=55
      NMASIMSCOCACCSOL=56
      NMASIMSCOCCORSOL=57
      NMASIMSCOCAITINS=58
      NMASIMSCSONUCSOL=59
      NMASIMSCSOAITSOL=60
      NMASIMSCSOACCSOL=61
      NMASIMSCSOCORSOL=62
!
      NMASCLPRSUAITSOL1=63
      NMASCLPRSUACCSOL1=64
      NMASCLPRSUCORSOL1=65
      NMASCLPRSUAITSOL2=66
      NMASCLPRSUACCSOL2=67
      NMASCLPRSUCORSOL2=68
!
      NMASCONDSUNUCSOL=69
      NMASCONDSUAITSOL=70
      NMASCONDSUACCSOL=71
      NMASCONDSUCORSOL=72
      NMASCONDSUAITINS=73
      NMASNUCLSUNUCSOL=74
      NMASCONDOCNUCSOL= 0 ! stored in NMASCONDSONUCSOL
      NMASCONDOCAITSOL= 0 ! stored in NMASCONDSOAITSOL
      NMASCONDOCACCSOL= 0 ! stored in NMASCONDSOACCSOL
      NMASCONDOCCORSOL= 0 ! stored in NMASCONDSOCORSOL
      NMASCONDOCAITINS= 0 ! stored in NMASCONDSOAITINS
      NMASCONDSONUCSOL=75
      NMASCONDSOAITSOL=76
      NMASCONDSOACCSOL=77
      NMASCONDSOCORSOL=78
      NMASCONDSOAITINS=79
!
      NMASCOAGSUINTR12=80
      NMASCOAGSUINTR13=81
      NMASCOAGSUINTR14=82
      NMASCOAGSUINTR15=83
      NMASCOAGOCINTR12= 0 ! stored in NMASCOAGSOINTR12
      NMASCOAGOCINTR13= 0 ! stored in NMASCOAGSOINTR13
      NMASCOAGOCINTR14= 0 ! stored in NMASCOAGSOINTR14
      NMASCOAGOCINTR15= 0 ! stored in NMASCOAGSOINTR15
      NMASCOAGSOINTR12=84
      NMASCOAGSOINTR13=85
      NMASCOAGSOINTR14=86
      NMASCOAGSOINTR15=87
      NMASCOAGSUINTR23=88
      NMASCOAGBCINTR23=89
      NMASCOAGOCINTR23=90
      NMASCOAGSOINTR23=91
      NMASCOAGSUINTR24=92
      NMASCOAGBCINTR24=93
      NMASCOAGOCINTR24=94
      NMASCOAGSOINTR24=95
      NMASCOAGSUINTR34=96
      NMASCOAGBCINTR34=97
      NMASCOAGOCINTR34=98
      NMASCOAGSSINTR34=99
      NMASCOAGSOINTR34=100
!
      NMASCOAGBCINTR53=101
      NMASCOAGOCINTR53=102
      NMASCOAGBCINTR54=103
      NMASCOAGOCINTR54=104
!
      NMASAGEDSUINTR52=105
      NMASAGEDBCINTR52=106
      NMASAGEDOCINTR52=107
      NMASAGEDSOINTR52=108
!
      NMASMERGSUINTR12=109
      NMASMERGOCINTR12=  0 ! separate SO component
      NMASMERGSOINTR12=110
      NMASMERGSUINTR23=111
      NMASMERGBCINTR23=112
      NMASMERGOCINTR23=113
      NMASMERGSOINTR23=114
      NMASMERGSUINTR34=115
      NMASMERGSSINTR34=116
      NMASMERGBCINTR34=117
      NMASMERGOCINTR34=118
      NMASMERGSOINTR34=119
      NMASPROCSUINTR23=120
      NMASPROCBCINTR23=121
      NMASPROCOCINTR23=122
      NMASPROCSOINTR23=123
!
! .. below are new ones for dust & modes 6/7 to be integrated
      NMASPRIMDUACCSOL= 0 ! no DU in this setup
      NMASPRIMDUCORSOL= 0 ! no DU in this setup
      NMASPRIMDUACCINS= 0 ! no DU in this setup
      NMASPRIMDUCORINS= 0 ! no DU in this setup
      NMASDDEPDUACCSOL= 0 ! no DU in this setup
      NMASDDEPDUCORSOL= 0 ! no DU in this setup
      NMASDDEPDUACCINS= 0 ! no DU in this setup
      NMASDDEPDUCORINS= 0 ! no DU in this setup
      NMASNUSCDUACCSOL= 0 ! no DU in this setup
      NMASNUSCDUCORSOL= 0 ! no DU in this setup
      NMASNUSCDUACCINS= 0 ! no DU in this setup
      NMASNUSCDUCORINS= 0 ! no DU in this setup
      NMASIMSCDUACCSOL= 0 ! no DU in this setup
      NMASIMSCDUCORSOL= 0 ! no DU in this setup
      NMASIMSCDUACCINS= 0 ! no DU in this setup
      NMASIMSCDUCORINS= 0 ! no DU in this setup
      NMASCONDSUACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSUCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDOCACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDOCCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSOACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSOCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGDUINTR34= 0 ! no DU in this setup
      NMASCOAGDUINTR64= 0 ! no DU in this setup
      NMASAGEDSUINTR63= 0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR63= 0 ! no DU in this setup
      NMASAGEDOCINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSUINTR74= 0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR74= 0 ! no DU in this setup
      NMASAGEDOCINTR74= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR74= 0 ! no BC/OC/SO in this setup
      NMASMERGDUINTR34= 0 ! no DU in this setup

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
                'UKCA_INDICES_SUSSBCOCSO_5MODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_SUSSBCOCSO_5MODE

! ######################################################################
      SUBROUTINE UKCA_INDICES_SUSSBCOCSO_4MODE

      IMPLICIT NONE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
                'UKCA_INDICES_SUSSBCOCSO_4MODE',zhook_in,zhook_handle)

! Main array lengths and switches

      NTRAER=20          ! # of aerosol advected tracers
      NBUDAER=104        ! # of aerosol budget fields

! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  20    16    130     9 = 36+130+ 9=  175 (orgv1 ,traqu9 )
! ..  SUSSBCOC  20    16    130    38 = 36+130+38=  204 (orgv1 ,traqu38)
! ..  SUSSBCOC  20    76    130     9 = 36+190+ 9=  235 (orgv1c,traqu9 )
! ..  SUSSBCOC  20    76    130    38 = 36+190+38=  264 (orgv1c,traqu38)
!                         (104+26)
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 20    13  =  33 (orgv1 )
! ..  SUSSBCOC 20    52  =  72 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. 104 aerosol budget indices for SUSSBCOCSO_4MODE [SO in SO]
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 6
      NMASPRIMBCAITINS= 0 ! BC only emitted to soluble
      NMASPRIMOCAITSOL= 7
      NMASPRIMOCAITINS= 0 ! OC only emitted to soluble
!
      NMASDDEPSUNUCSOL= 8
      NMASDDEPSUAITSOL= 9
      NMASDDEPSUACCSOL=10
      NMASDDEPSUCORSOL=11
      NMASDDEPSSACCSOL=12
      NMASDDEPSSCORSOL=13
      NMASDDEPBCAITSOL=14
      NMASDDEPBCACCSOL=15
      NMASDDEPBCCORSOL=16
      NMASDDEPBCAITINS= 0  ! only soluble modes
      NMASDDEPOCNUCSOL= 0  ! stored in NMASDDEPSONUCSOL
      NMASDDEPOCAITSOL=17
      NMASDDEPOCACCSOL=18
      NMASDDEPOCCORSOL=19
      NMASDDEPOCAITINS= 0  ! only soluble modes
      NMASDDEPSONUCSOL=20
      NMASDDEPSOAITSOL=21
      NMASDDEPSOACCSOL=22
      NMASDDEPSOCORSOL=23
!
      NMASNUSCSUNUCSOL=24
      NMASNUSCSUAITSOL=25
      NMASNUSCSUACCSOL=26
      NMASNUSCSUCORSOL=27
      NMASNUSCSSACCSOL=28
      NMASNUSCSSCORSOL=29
      NMASNUSCBCAITSOL=30
      NMASNUSCBCACCSOL=31
      NMASNUSCBCCORSOL=32
      NMASNUSCBCAITINS= 0 ! only soluble modes
      NMASNUSCOCNUCSOL= 0 ! stored in NMASNUSCSONUCSOL
      NMASNUSCOCAITSOL=33
      NMASNUSCOCACCSOL=34
      NMASNUSCOCCORSOL=35
      NMASNUSCOCAITINS= 0 ! only soluble modes
      NMASNUSCSONUCSOL=36
      NMASNUSCSOAITSOL=37
      NMASNUSCSOACCSOL=38
      NMASNUSCSOCORSOL=39
!
      NMASIMSCSUNUCSOL=40
      NMASIMSCSUAITSOL=41
      NMASIMSCSUACCSOL=42
      NMASIMSCSUCORSOL=43
      NMASIMSCSSACCSOL=44
      NMASIMSCSSCORSOL=45
      NMASIMSCBCAITSOL=46
      NMASIMSCBCACCSOL=47
      NMASIMSCBCCORSOL=48
      NMASIMSCBCAITINS= 0 ! only soluble modes
      NMASIMSCOCNUCSOL= 0 ! stored in NMASIMSCSONUCSOL
      NMASIMSCOCAITSOL=49
      NMASIMSCOCACCSOL=50
      NMASIMSCOCCORSOL=51
      NMASIMSCOCAITINS= 0 ! only soluble modes
      NMASIMSCSONUCSOL=52
      NMASIMSCSOAITSOL=53
      NMASIMSCSOACCSOL=54
      NMASIMSCSOCORSOL=55
!
      NMASCLPRSUAITSOL1=56
      NMASCLPRSUACCSOL1=57
      NMASCLPRSUCORSOL1=58
      NMASCLPRSUAITSOL2=59
      NMASCLPRSUACCSOL2=60
      NMASCLPRSUCORSOL2=61
!
      NMASCONDSUNUCSOL=62
      NMASCONDSUAITSOL=63
      NMASCONDSUACCSOL=64
      NMASCONDSUCORSOL=65
      NMASCONDSUAITINS= 0 ! only soluble modes
      NMASNUCLSUNUCSOL=66
      NMASCONDOCNUCSOL= 0 ! stored in NMASCONDSONUCSOL
      NMASCONDOCAITSOL= 0 ! stored in NMASCONDSOAITSOL
      NMASCONDOCACCSOL= 0 ! stored in NMASCONDSOACCSOL
      NMASCONDOCCORSOL= 0 ! stored in NMASCONDSOCORSOL
      NMASCONDOCAITINS= 0 ! stored in NMASCONDSOAITINS
      NMASCONDSONUCSOL=67
      NMASCONDSOAITSOL=68
      NMASCONDSOACCSOL=69
      NMASCONDSOCORSOL=70
      NMASCONDSOAITINS= 0 ! only soluble modes
!
      NMASCOAGSUINTR12=71
      NMASCOAGSUINTR13=72
      NMASCOAGSUINTR14=73
      NMASCOAGSUINTR15= 0 ! only soluble modes
      NMASCOAGOCINTR12= 0 ! stored in NMASCOAGSOINTR12
      NMASCOAGOCINTR13= 0 ! stored in NMASCOAGSOINTR13
      NMASCOAGOCINTR14= 0 ! stored in NMASCOAGSOINTR14
      NMASCOAGOCINTR15= 0 ! stored in NMASCOAGSOINTR15
      NMASCOAGSOINTR12=74
      NMASCOAGSOINTR13=75
      NMASCOAGSOINTR14=76
      NMASCOAGSOINTR15= 0 ! only soluble modes
      NMASCOAGSUINTR23=77
      NMASCOAGBCINTR23=78
      NMASCOAGOCINTR23=79
      NMASCOAGSOINTR23=80
      NMASCOAGSUINTR24=81
      NMASCOAGBCINTR24=82
      NMASCOAGOCINTR24=83
      NMASCOAGSOINTR24=84
      NMASCOAGSUINTR34=85
      NMASCOAGBCINTR34=86
      NMASCOAGOCINTR34=87
      NMASCOAGSSINTR34=88
      NMASCOAGSOINTR34=89
!
      NMASCOAGBCINTR53= 0 ! only soluble modes
      NMASCOAGOCINTR53= 0 ! only soluble modes
      NMASCOAGBCINTR54= 0 ! only soluble modes
      NMASCOAGOCINTR54= 0 ! only soluble modes
!
      NMASAGEDSUINTR52= 0 ! only soluble modes
      NMASAGEDBCINTR52= 0 ! only soluble modes
      NMASAGEDOCINTR52= 0 ! only soluble modes
      NMASAGEDSOINTR52= 0 ! only soluble modes
!
      NMASMERGSUINTR12=90
      NMASMERGOCINTR12= 0 ! separate SO component
      NMASMERGSOINTR12=91
      NMASMERGSUINTR23=92
      NMASMERGBCINTR23=93
      NMASMERGOCINTR23=94
      NMASMERGSOINTR23=95
      NMASMERGSUINTR34=96
      NMASMERGSSINTR34=97
      NMASMERGBCINTR34=98
      NMASMERGOCINTR34=99
      NMASMERGSOINTR34=100
      NMASPROCSUINTR23=101
      NMASPROCBCINTR23=102
      NMASPROCOCINTR23=103
      NMASPROCSOINTR23=104
!
! .. below are new ones for dust & modes 6/7 to be integrated
      NMASPRIMDUACCSOL= 0 ! no DU in this setup
      NMASPRIMDUCORSOL= 0 ! no DU in this setup
      NMASPRIMDUACCINS= 0 ! no DU in this setup
      NMASPRIMDUCORINS= 0 ! no DU in this setup
      NMASDDEPDUACCSOL= 0 ! no DU in this setup
      NMASDDEPDUCORSOL= 0 ! no DU in this setup
      NMASDDEPDUACCINS= 0 ! no DU in this setup
      NMASDDEPDUCORINS= 0 ! no DU in this setup
      NMASNUSCDUACCSOL= 0 ! no DU in this setup
      NMASNUSCDUCORSOL= 0 ! no DU in this setup
      NMASNUSCDUACCINS= 0 ! no DU in this setup
      NMASNUSCDUCORINS= 0 ! no DU in this setup
      NMASIMSCDUACCSOL= 0 ! no DU in this setup
      NMASIMSCDUCORSOL= 0 ! no DU in this setup
      NMASIMSCDUACCINS= 0 ! no DU in this setup
      NMASIMSCDUCORINS= 0 ! no DU in this setup
      NMASCONDSUACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSUCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDOCACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDOCCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSOACCINS= 0 ! no acc-ins cor-ins in this setup
      NMASCONDSOCORINS= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSUINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGOCINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR16= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGSOINTR17= 0 ! no acc-ins cor-ins in this setup
      NMASCOAGDUINTR34= 0 ! no DU in this setup
      NMASCOAGDUINTR64= 0 ! no DU in this setup
      NMASAGEDSUINTR63= 0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR63= 0 ! no DU in this setup
      NMASAGEDOCINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSUINTR74= 0 ! no acc-ins cor-ins in this setup
      NMASAGEDDUINTR74= 0 ! no DU in this setup
      NMASAGEDOCINTR74= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR74= 0 ! no BC/OC/SO in this setup
      NMASMERGDUINTR34= 0 ! no DU in this setup

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
                'UKCA_INDICES_SUSSBCOCSO_4MODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_SUSSBCOCSO_4MODE

! ######################################################################
      SUBROUTINE UKCA_INDICES_SUSS_4MODE

      IMPLICIT NONE


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
            'UKCA_INDICES_SUSS_4MODE',zhook_in,zhook_handle)

! Main array lengths and switches
      NTRAER=10          ! # of aerosol advected tracers
      NBUDAER=46         ! # of aerosol budget fields

!
! For sv1         : NTRAG=13, NADVG=11
! For sv1_coupled : NTRAG=74, NADVG=50
!
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU             = NVTOT
! ..  SUSS      10    13     74     9 = 23 +  74 + 9 = 106(sv1 ,traqu9 )
! ..  SUSS      10    13     74     38= 23 +  74 + 38= 135(sv1 ,traqu38)
! ..  SUSS      10    74     74     9 = 84 +  74 + 9 = 167(sv1c,traqu9 )
! ..  SUSS      10    74     74     38= 84 +  74 + 38= 196(sv1c,traqu38)
!
!            NTRAER+NADVG  NTRA
! ..  SUSS      10   11  =  21 (sv1 )
! ..  SUSS      10   50  =  60 (sv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are 46 aerosol budget variables for SUSS aerosol system
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASPRIMBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASPRIMOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASPRIMOCAITINS= 0 ! no BC/OC/SO in this setup
!
      NMASDDEPSUNUCSOL= 6
      NMASDDEPSUAITSOL= 7
      NMASDDEPSUACCSOL= 8
      NMASDDEPSUCORSOL= 9
      NMASDDEPSSACCSOL=10
      NMASDDEPSSCORSOL=11
      NMASDDEPBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPBCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPBCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASDDEPSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPSOCORSOL= 0 ! no BC/OC/SO in this setup
!
      NMASNUSCSUNUCSOL=12
      NMASNUSCSUAITSOL=13
      NMASNUSCSUACCSOL=14
      NMASNUSCSUCORSOL=15
      NMASNUSCSSACCSOL=16
      NMASNUSCSSCORSOL=17
      NMASNUSCBCAITSOL= 0
      NMASNUSCBCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCBCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASNUSCSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCSOCORSOL= 0 ! no BC/OC/SO in this setup
!
      NMASIMSCSUNUCSOL=18
      NMASIMSCSUAITSOL=19
      NMASIMSCSUACCSOL=20
      NMASIMSCSUCORSOL=21
      NMASIMSCSSACCSOL=22
      NMASIMSCSSCORSOL=23
      NMASIMSCBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCBCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCBCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASIMSCSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCSOCORSOL= 0 ! no BC/OC/SO in this setup
!
      NMASCLPRSUAITSOL1=24
      NMASCLPRSUACCSOL1=25
      NMASCLPRSUCORSOL1=26
      NMASCLPRSUAITSOL2=27
      NMASCLPRSUACCSOL2=28
      NMASCLPRSUCORSOL2=29
!
      NMASCONDSUNUCSOL=30
      NMASCONDSUAITSOL=31
      NMASCONDSUACCSOL=32
      NMASCONDSUCORSOL=33
      NMASCONDSUAITINS= 0 ! only soluble modes in this setup
      NMASNUCLSUNUCSOL=34
      NMASCONDOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASCONDSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOCORSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOAITINS= 0 ! no BC/OC/SO in this setup
!
      NMASCOAGSUINTR12=35
      NMASCOAGSUINTR13=36
      NMASCOAGSUINTR14=37
      NMASCOAGSUINTR15= 0 ! only soluble modes in this setup
      NMASCOAGOCINTR12= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR13= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR14= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR15= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR12= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR13= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR14= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR15= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR23=38
      NMASCOAGBCINTR23= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR23= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR23= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR24=39
      NMASCOAGBCINTR24= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR24= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR24= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR34=40
      NMASCOAGBCINTR34= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR34= 0 ! no BC/OC/SO in this setup
      NMASCOAGSSINTR34=41
      NMASCOAGSOINTR34= 0 ! no BC/OC/SO in this setup
!
      NMASCOAGBCINTR53= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR53= 0 ! no BC/OC/SO in this setup
      NMASCOAGBCINTR54= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR54= 0 ! no BC/OC/SO in this setup
!
      NMASAGEDSUINTR52= 0 ! only soluble modes in this setup
      NMASAGEDBCINTR52= 0 ! no BC/OC/SO in this setup
      NMASAGEDOCINTR52= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR52= 0 ! no BC/OC/SO in this setup
!
      NMASMERGSUINTR12=42
      NMASMERGOCINTR12= 0 ! no BC/OC/SO in this setup
      NMASMERGSOINTR12= 0 ! no BC/OC/SO in this setup
      NMASMERGSUINTR23=43
      NMASMERGBCINTR23= 0 ! no BC/OC/SO in this setup
      NMASMERGOCINTR23= 0 ! no BC/OC/SO in this setup
      NMASMERGSOINTR23= 0 ! no BC/OC/SO in this setup
      NMASMERGSUINTR34=44
      NMASMERGSSINTR34=45
      NMASMERGBCINTR34= 0 ! no BC/OC/SO in this setup
      NMASMERGOCINTR34= 0 ! no BC/OC/SO in this setup
      NMASMERGSOINTR34= 0 ! no BC/OC/SO in this setup
      NMASPROCSUINTR23=46
      NMASPROCBCINTR23= 0 ! no BC/OC/SO in this setup
      NMASPROCOCINTR23= 0 ! no BC/OC/SO in this setup
      NMASPROCSOINTR23= 0 ! no BC/OC/SO in this setup
!
! .. below are new ones for dust & modes 6/7 to be integrated
      NMASPRIMDUACCSOL= 0 ! no DU in this setup
      NMASPRIMDUCORSOL= 0 ! no DU in this setup
      NMASPRIMDUACCINS= 0 ! no DU in this setup
      NMASPRIMDUCORINS= 0 ! no DU in this setup
      NMASDDEPDUACCSOL= 0 ! no DU in this setup
      NMASDDEPDUCORSOL= 0 ! no DU in this setup
      NMASDDEPDUACCINS= 0 ! no DU in this setup
      NMASDDEPDUCORINS= 0 ! no DU in this setup
      NMASNUSCDUACCSOL= 0 ! no DU in this setup
      NMASNUSCDUCORSOL= 0 ! no DU in this setup
      NMASNUSCDUACCINS= 0 ! no DU in this setup
      NMASNUSCDUCORINS= 0 ! no DU in this setup
      NMASIMSCDUACCSOL= 0 ! no DU in this setup
      NMASIMSCDUCORSOL= 0 ! no DU in this setup
      NMASIMSCDUACCINS= 0 ! no DU in this setup
      NMASIMSCDUCORINS= 0 ! no DU in this setup
      NMASCONDSUACCINS= 0 ! only soluble modes in this setup
      NMASCONDSUCORINS= 0 ! only soluble modes in this setup
      NMASCONDOCACCINS= 0 ! only soluble modes in this setup
      NMASCONDOCCORINS= 0 ! only soluble modes in this setup
      NMASCONDSOACCINS= 0 ! only soluble modes in this setup
      NMASCONDSOCORINS= 0 ! only soluble modes in this setup
      NMASCOAGSUINTR16= 0 ! only soluble modes in this setup
      NMASCOAGSUINTR17= 0 ! only soluble modes in this setup
      NMASCOAGOCINTR16= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR17= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR16= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR17= 0 ! no BC/OC/SO in this setup
      NMASCOAGDUINTR34= 0 ! no DU in this setup
      NMASCOAGDUINTR64= 0 ! no DU in this setup
      NMASAGEDSUINTR63= 0 ! only soluble modes in this setup
      NMASAGEDDUINTR63= 0 ! no DU in this setup
      NMASAGEDOCINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR63= 0 ! no BC/OC/SO in this setup
      NMASAGEDSUINTR74= 0 ! only soluble modes in this setup
      NMASAGEDDUINTR74= 0 ! no DU in this setup
      NMASAGEDOCINTR74= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR74= 0 ! no BC/OC/SO in this setup
      NMASMERGDUINTR34= 0 ! no DU in this setup

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
            'UKCA_INDICES_SUSS_4MODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_SUSS_4MODE

! ######################################################################
      SUBROUTINE UKCA_INDICES_DUonly_2MODE

      IMPLICIT NONE


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
               'UKCA_INDICES_DUONLY_2MODE',zhook_in,zhook_handle)

! Main array lengths and switches
      NTRAER=4           ! # of aerosol advected tracers
      NBUDAER=8          ! # of aerosol budget fields
!
! For nochem: NTRAG= 2, NADVG= 2
!
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU             = NVTOT
! ..  DU-only    4     2      8     9 =  6 +   8 + 9 =  23 (nochem)
! ..  DU-only    4     2      8     38=  6 +   8 + 38=  52 (nochem)
!
!            NTRAER+NADVG  NTRA
! ..  DU-only    4    2  =   6
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are  8 aerosol budget variables for DU-only aerosol system
!
      NMASPRIMSUAITSOL= 0 ! no SO4 or SS in this setup
      NMASPRIMSUACCSOL= 0 ! no SO4 or SS in this setup
      NMASPRIMSUCORSOL= 0 ! no SO4 or SS in this setup
      NMASPRIMSSACCSOL= 0 ! no SO4 or SS in this setup
      NMASPRIMSSCORSOL= 0 ! no SO4 or SS in this setup
      NMASPRIMBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASPRIMBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASPRIMOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASPRIMOCAITINS= 0 ! no BC/OC/SO in this setup
!
      NMASDDEPSUNUCSOL= 0 ! no SO4 or SS in this setup
      NMASDDEPSUAITSOL= 0 ! no SO4 or SS in this setup
      NMASDDEPSUACCSOL= 0 ! no SO4 or SS in this setup
      NMASDDEPSUCORSOL= 0 ! no SO4 or SS in this setup
      NMASDDEPSSACCSOL= 0 ! no SO4 or SS in this setup
      NMASDDEPSSCORSOL= 0 ! no SO4 or SS in this setup
      NMASDDEPBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPBCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPBCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASDDEPSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPSOCORSOL= 0 ! no BC/OC/SO in this setup
!
      NMASNUSCSUNUCSOL= 0 ! no SO4 or SS in this setup
      NMASNUSCSUAITSOL= 0 ! no SO4 or SS in this setup
      NMASNUSCSUACCSOL= 0 ! no SO4 or SS in this setup
      NMASNUSCSUCORSOL= 0 ! no SO4 or SS in this setup
      NMASNUSCSSACCSOL= 0 ! no SO4 or SS in this setup
      NMASNUSCSSCORSOL= 0 ! no SO4 or SS in this setup
      NMASNUSCBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCBCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCBCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASNUSCSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCSOCORSOL= 0 ! no BC/OC/SO in this setup
!
      NMASIMSCSUNUCSOL= 0 ! no SO4 or SS in this setup
      NMASIMSCSUAITSOL= 0 ! no SO4 or SS in this setup
      NMASIMSCSUACCSOL= 0 ! no SO4 or SS in this setup
      NMASIMSCSUCORSOL= 0 ! no SO4 or SS in this setup
      NMASIMSCSSACCSOL= 0 ! no SO4 or SS in this setup
      NMASIMSCSSCORSOL= 0 ! no SO4 or SS in this setup
      NMASIMSCBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCBCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCBCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASIMSCSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCSOCORSOL= 0 ! no BC/OC/SO in this setup
!
      NMASCLPRSUAITSOL1=0 ! no SO4 or SS in this setup
      NMASCLPRSUACCSOL1=0 ! no SO4 or SS in this setup
      NMASCLPRSUCORSOL1=0 ! no SO4 or SS in this setup
      NMASCLPRSUAITSOL2=0 ! no SO4 or SS in this setup
      NMASCLPRSUACCSOL2=0 ! no SO4 or SS in this setup
      NMASCLPRSUCORSOL2=0 ! no SO4 or SS in this setup
!
      NMASCONDSUNUCSOL= 0 ! no SO4 or SS in this setup
      NMASCONDSUAITSOL= 0 ! no SO4 or SS in this setup
      NMASCONDSUACCSOL= 0 ! no SO4 or SS in this setup
      NMASCONDSUCORSOL= 0 ! no SO4 or SS in this setup
      NMASCONDSUAITINS= 0 ! only soluble modes in this setup
      NMASNUCLSUNUCSOL= 0 ! no SO4 or SS in this setup
      NMASCONDOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASCONDSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOCORSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOAITINS= 0 ! no BC/OC/SO in this setup
!
      NMASCOAGSUINTR12= 0 ! no SO4 or SS in this setup
      NMASCOAGSUINTR13= 0 ! no SO4 or SS in this setup
      NMASCOAGSUINTR14= 0 ! no SO4 or SS in this setup
      NMASCOAGSUINTR15= 0 ! only soluble modes in this setup
      NMASCOAGOCINTR12= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR13= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR14= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR15= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR12= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR13= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR14= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR15= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR23= 0 ! no SO4 or SS in this setup
      NMASCOAGBCINTR23= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR23= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR23= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR24= 0 ! no SO4 or SS in this setup
      NMASCOAGBCINTR24= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR24= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR24= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR34= 0 ! no SO4 or SS in this setup
      NMASCOAGBCINTR34= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR34= 0 ! no BC/OC/SO in this setup
      NMASCOAGSSINTR34= 0 ! no SO4 or SS in this setup
      NMASCOAGSOINTR34= 0 ! no BC/OC/SO in this setup
!
      NMASCOAGBCINTR53= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR53= 0 ! no BC/OC/SO in this setup
      NMASCOAGBCINTR54= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR54= 0 ! no BC/OC/SO in this setup
!
      NMASAGEDSUINTR52= 0 ! only soluble modes in this setup
      NMASAGEDBCINTR52= 0 ! no BC/OC/SO in this setup
      NMASAGEDOCINTR52= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR52= 0 ! no BC/OC/SO in this setup
!
      NMASMERGSUINTR12= 0 ! no SO4 or SS in this setup
      NMASMERGOCINTR12= 0 ! no BC/OC/SO in this setup
      NMASMERGSOINTR12= 0 ! no BC/OC/SO in this setup
      NMASMERGSUINTR23= 0 ! no SO4 or SS in this setup
      NMASMERGBCINTR23= 0 ! no BC/OC/SO in this setup
      NMASMERGOCINTR23= 0 ! no BC/OC/SO in this setup
      NMASMERGSOINTR23= 0 ! no BC/OC/SO in this setup
      NMASMERGSUINTR34= 0 ! no SO4 or SS in this setup
      NMASMERGSSINTR34= 0 ! no SO4 or SS in this setup
      NMASMERGBCINTR34= 0 ! no BC/OC/SO in this setup
      NMASMERGOCINTR34= 0 ! no BC/OC/SO in this setup
      NMASMERGSOINTR34= 0 ! no BC/OC/SO in this setup
      NMASPROCSUINTR23= 0 ! no SO4 or SS in this setup
      NMASPROCBCINTR23= 0 ! no BC/OC/SO in this setup
      NMASPROCOCINTR23= 0 ! no BC/OC/SO in this setup
      NMASPROCSOINTR23= 0 ! no BC/OC/SO in this setup
!
! .. below are new ones for dust & modes 6/7 to be integrated
      NMASPRIMDUACCSOL= 0 ! DU emitted into insoluble modes
      NMASPRIMDUCORSOL= 0 ! DU emitted into insoluble modes
      NMASPRIMDUACCINS= 1
      NMASPRIMDUCORINS= 2
      NMASDDEPDUACCSOL= 0 ! no aged DU in this setup
      NMASDDEPDUCORSOL= 0 ! no aged DU in this setup
      NMASDDEPDUACCINS= 3
      NMASDDEPDUCORINS= 4
      NMASNUSCDUACCSOL= 0 ! no aged DU in this setup
      NMASNUSCDUCORSOL= 0 ! no aged DU in this setup
      NMASNUSCDUACCINS= 5
      NMASNUSCDUCORINS= 6
      NMASIMSCDUACCSOL= 0 ! no aged DU in this setup
      NMASIMSCDUCORSOL= 0 ! no aged DU in this setup
      NMASIMSCDUACCINS= 7
      NMASIMSCDUCORINS= 8
      NMASCONDSUACCINS= 0 ! no SO4 or SS in this setup
      NMASCONDSUCORINS= 0 ! no SO4 or SS in this setup
      NMASCONDOCACCINS= 0 ! no BC/OC/SO in this setup
      NMASCONDOCCORINS= 0 ! no BC/OC/SO in this setup
      NMASCONDSOACCINS= 0 ! no BC/OC/SO in this setup
      NMASCONDSOCORINS= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR16= 0 ! no SO4 or SS in this setup
      NMASCOAGSUINTR17= 0 ! no SO4 or SS in this setup
      NMASCOAGOCINTR16= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR17= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR16= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR17= 0 ! no BC/OC/SO in this setup
      NMASCOAGDUINTR34= 0 ! no aged DU in this setup
      NMASCOAGDUINTR64= 0 ! no aged DU in this setup
      NMASAGEDSUINTR63= 0 ! no aged DU in this setup
      NMASAGEDDUINTR63= 0 ! no aged DU in this setup
      NMASAGEDOCINTR63= 0 ! no aged DU in this setup
      NMASAGEDSOINTR63= 0 ! no aged DU in this setup
      NMASAGEDSUINTR74= 0 ! no aged DU in this setup
      NMASAGEDDUINTR74= 0 ! no aged DU in this setup
      NMASAGEDOCINTR74= 0 ! no aged DU in this setup
      NMASAGEDSOINTR74= 0 ! no aged DU in this setup
      NMASMERGDUINTR34= 0 ! no aged DU in this setup

      IF (lhook) CALL dr_hook('UKCA_SETUP_INDICES:'//                   &
               'UKCA_INDICES_DUONLY_2MODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INDICES_DUonly_2MODE

      END MODULE UKCA_SETUP_INDICES
