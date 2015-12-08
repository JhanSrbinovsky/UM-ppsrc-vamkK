! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE inijtab_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Subroutine to initialise photolysis tables
!      AJ2A:     O2       +  hv             ->   O     +   O
!      AJ2B:     O2       +  hv             ->   O     +   O(1D)
!      AJ3:      O3       +  hv             ->   O2    +   O(3P)
!      AJ3A:     O3       +  hv(<310nm)     ->   O2    +   O(1D)
!      AJNO:     NO       +  hv             ->   N     +   O
!      AJNO2:    NO2      +  hv             ->   NO    +   O(3P)
!      AJNO31:   NO3      +  hv             ->   NO    +   O2
!      AJNO32:   NO3      +  hv             ->   NO2   +   O
!      AJHNO3:   HNO3     +  hv             ->   NO2   +   OH
!      AJPNA:    HO2NO2   +  hv             ->   NO2   +   HO2 / NO3 + OH
!      AJN2O:    N2O      +  hv             ->   N2    +   O(1D)
!      AJN2O5:   N2O5     +  hv             ->   NO2   +   NO3
!      AJH2O:    H2O      +  hv             ->   OH    +   H
!      AJH2O2:   H2O2     +  hv             ->   OH    +   OH
!      AJCNITA:  ClONO2   +  hv             ->   Cl    +   NO3
!      AJCNITB:  ClONO2   +  hv             ->   ClO   +   NO2
!      AJF11:    CFCl3    +  hv             ->   3Cl
!      AJF12:    CF2Cl2   +  hv             ->   2Cl
!      AJF22:    CHF2Cl   +  hv             ->   Cl
!      AJF113:   CF2CFCl3 +  hv             ->   3Cl
!      AJCH3CL:  CH3Cl    +  hv             ->   CH3   +   Cl
!      AJCCL4:   CCl4     +  hv             ->   4Cl
!      AJHCL:    HCl      +  hv             ->   H     +   Cl
!      AJHOCL:   HOCl     +  hv             ->   OH    +   Cl
!      AJOCLO:   OClO     +  hv             ->   O     +   ClO
!      AJCL2O2:  Cl2O2    +  hv             ->   O2    +   2Cl
!      AJBRO:    BrO      +  hv             ->   Br    +   O
!      AJHOBR:   HOBr     +  hv             ->   Br    +   OH
!      AJBRNO3:  BrNO3    +  hv             ->   Br    +   NO3 / Bro + NO2
!      AJBRCL:   BrCl     +  hv             ->   Br    +   Cl
!      AJC2OA:   CH2O     +  hv             ->   H     +   CHO
!      AJC2OB:   CH2O     +  hv             ->   H2    +   CO
!      AJMHP :   CH3OOH   +  hv             ->   CH3O  +   OH
!      AJCH3BR:  CH3Br    +  hv             ->   CH3   +   Br
!      AJMCFM:   CH3CCl3  +  hv             ->   CHx   +   3Cl
!      AJCH4:    CH4      +  hv             ->   CH3   +   H
!      AJF12B1:  CBrClF2  +  hv             ->   Cxx   +   Cl   + Br
!      AJF13B1:  CBrF3    +  hv             ->   Cxx   +   Br
!      AJCOF2:   COF2     +  hv             ->   Cx    +   2HF
!      AJCOFCL:  COFCl    +  hv             ->   Cx    +   Cl + HF
!      AJCO2:    CO2      +  hv             ->   CO    +   O(3P)
!      AJCOS:    COS      +  hv             ->   CO    +   S
!      AJHONO:   HONO     +  hv             ->   OH    +   NO
!      AJMENA:   MeONO2   +  hv             ->   HO2   +   HCHO + NO2
!      AJCHBr3:  CHBr3    +  hv             ->   HBr   +   2Br
!      AJDBRM :  CH2Br2   +  hv             ->   2Br   +   H2O
!      AJCS2:    CS2      +  hv             ->   COS   +   SO2
!      AJH2SO4:  H2SO4    +  hv             ->   SO3   +   OH
!      AJSO3:    SO3      +  hv             ->   SO2   +   O(3P)


!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CONTAINS
SUBROUTINE inijtab(mype, l_ukca_fastjx)
USE ukca_stdto3_mod, ONLY: stdto3, nin
USE ukca_tbjs_mod, ONLY:     tabj2a   , tabj2b   , tabj3    ,    &
       tabj3a   , tabjno   , tabjno31 , tabjno32 , tabjno2  ,    & 
       tabjn2o5 , tabjhno3 , tabjh2o  , tabjh2o2 , tabjf11  ,    &
       tabjf12  , tabjf22  , tabjf113 , tabjch3cl, tabjccl4 ,    &
       tabjcnita, tabjcnitb, tabjhcl  , tabjhocl , tabjpna  ,    &
       tabjcl2o2, tabjn2o  , tabjbrcl , tabjbrno3, tabjhobr ,    &
       tabjbro  , tabjoclo , tabjc2oa , tabjc2ob , tabjmhp  ,    &
       tabjmcfm , tabjch3br, tabjf12b1, tabjf13b1, tabjcof2 ,    &
       tabjcofcl, tabjch4  , tabjcos  , tabjso2  , tabjco2  ,    &
       tabjhono , tabjmena , tabjchbr3, tabjdbrm , tabjcs2  ,    &
       tabjh2so4, tabjso3  , tabpres  , tabang,    tabt     ,    &
       tabo3,     sao3c
USE ukca_crossec_mod, ONLY: accl4, ach2o, acl2o2, amcfm, acnita, &
    acnitb, aco2 , af11 , af113, af12 , af22 , ah2o , ah2o2,     &
    ahcl , ahno2, ahno3, ahocl, ach3cl, an2o , an2o5, ano  ,     &
    ano2 , ano31, ano32, ao2  , ao2a , ao2b , abrcl, abrno3,     &
    abro , ahobr, aoclo, ac2oa, ac2ob, amhp , ao2sr, ao3  ,      &
    apna , acof2, acofcl, ach3br, af12b1, af13b1, scs  , qeno2,  &
    qeo1d, quanta, ach4 , wavecm, wavenm, aobro, ahono, aocs ,   &
    aso2 , amena, achbr3, adbrm, acs2 , ah2so4, aso3  
USE ukca_parpho_mod, ONLY: jplevp1, jplev, jpchi, jps90, jpchin,               &
                           jpwav, jplo, jphi, jptem, jpo3p, jps90,             &
                           szamax, tmin, tmax, o3min, o3max

! Module procedures used by this subroutine
USE fill_spectra_mod, ONLY:fill_spectra
USE scatcs_mod, ONLY: scatcs
USE setzen_mod, ONLY: setzen
USE settab_mod, ONLY: settab
USE lymana_mod, ONLY: lymana 
USE acsn2o5_mod, ONLY: acsn2o5
USE acsn2o_mod, ONLY: acsn2o
USE acsccl4_mod, ONLY: acsccl4
USE acsf11_mod, ONLY: acsf11
USE acsf12_mod, ONLY: acsf12
USE acsf22_mod, ONLY: acsf22
USE acsmc_mod, ONLY: acsmc
USE quanto12_mod, ONLY: quanto12
USE acso3_mod, ONLY: acso3
USE acshno3_mod, ONLY: acshno3
USE acscnit_mod, ONLY: acscnit
USE acsh2o2_mod, ONLY: acsh2o2
USE acsbrcl_mod, ONLY: acsbrcl
USE acsno2_mod, ONLY: acsno2
USE acscos_mod, ONLY: acscos
USE acsmena_mod, ONLY: acsmena
USE acsdbrm_mod, ONLY: acsdbrm
USE acscs2_mod, ONLY: acscs2
USE acsh2so4_mod, ONLY: acsh2so4
USE acsso3_mod, ONLY: acsso3
USE acsno_mod, ONLY: acsno
USE acssr_mod, ONLY: acssr

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE chsunits_mod, ONLY : nunits

IMPLICIT NONE

!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200



! Subroutine interface
INTEGER, INTENT(IN) :: mype  ! processor number
LOGICAL, INTENT(IN) :: l_ukca_fastjx

! Local variables
! do not read STDTO3 input file. Instead, use values from the module.
LOGICAL, PARAMETER :: read_stdto3 = .FALSE.

REAL, PARAMETER :: gg=9.81,ra=287.06,colfac=2.132e20
!     O2 mixing ratio
REAL, PARAMETER :: r2o2=0.209
!     Ground albedo
REAL, PARAMETER :: albedo=0.30

!     Whether scattering is on or not.
LOGICAL, PARAMETER :: lscat =.TRUE.

INTEGER :: jc
INTEGER :: jl
INTEGER :: jlu
INTEGER :: jll
INTEGER :: jo
INTEGER :: jt
INTEGER :: jw
INTEGER :: l

REAL :: o3above

REAL :: cinf
REAL :: csup
REAL :: dpdif
REAL :: pdif
REAL :: qu

REAL :: factor(jpwav)
REAL :: do2c  (jplev)    ! Average number density of O2 within
!                          each level in molecules per cm^3.
REAL :: do3c  (jplev)    ! Average number density of O3 within
!                          each level in molecules per cm^3.
REAL :: drsc  (jplev)    ! [M] within level
REAL :: do3cu (jplev)    ! Scaled O3 profiles
REAL :: pres  (jplevp1)  ! Pressure (mb) at the edge of each level
REAL :: presc (jplev)    ! Pressure (mb) at the centre of each level
REAL :: o3vmrc(jplev)
REAL :: sao2  (jplevp1)  ! Vertical O2 column above the edge of
!                          each level in molecules per cm^2

REAL :: tabs  (jplev,jpchi,jpwav)
REAL :: tempc (jplev)            ! Temperature (K) at the centre of each level
REAL :: tspo2 (jplev,jpchi)
REAL :: zenmax(jplev)
REAL :: alt   (jplevp1)          ! Altitude at the edge of each level in km
REAL :: altc  (jplev)            ! Altitude at the centre of each level in km
REAL :: dalt  (jplev)            ! Thickness of level (km)

REAL :: rdtemp(nin)
REAL :: rdpres(nin)
REAL :: rdo3v (nin)

INTEGER :: last_wn , ios

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


CHARACTER(LEN=256) :: fstdto3


IF (lhook) CALL dr_hook('INIJTAB',zhook_in,zhook_handle)
IF (read_stdto3) THEN
  !     Read in standard atmosphere temperature and O3
  CALL get_file(ukcasto3_unit, fstdto3, 256, ios)
  OPEN(ukcasto3_unit,file=fstdto3)
  REWIND(ukcasto3_unit)
  DO l=1,nin
    READ(ukcasto3_unit,*) rdpres(l),rdtemp(l),rdo3v(l)
  END DO
  CLOSE(ukcasto3_unit)
ELSE
  ! Read standard atmosphere data from string variable
  DO l=1,nin
    READ(stdto3(l),*) rdpres(l),rdtemp(l),rdo3v(l)
  END DO
END IF

!     Set up j table inter levels equally spaced in log(p)
pres(      1)=rdpres(  1)
pres(jplevp1)=rdpres(nin)
pdif=LOG(rdpres(1)) - LOG(rdpres(nin))
dpdif=pdif/(REAL(jplevp1-1))
DO l=2,jplevp1-1
  pres (l)=EXP(LOG(pres(l-1)) - dpdif)
END DO
!     Pressure at jtable mid levels
DO l=1,jplev
  presc  (l)=EXP(0.5*(LOG(pres(l)) + LOG(pres(l+1))))
  tabpres(l)=presc(l)
END DO

!     Interpolate temperature and O3 to level centres
DO l=1,jplev
  DO jl=1,nin
    IF (rdpres(jl) < presc(l)) EXIT
  END DO
  jlu=jl
  jll=jl-1
  csup=(LOG(rdpres(jll)) - LOG(presc (l  )))/                                  &
    (LOG(rdpres(jll)) - LOG(rdpres(jlu)))
  cinf=1.0 - csup
  tempc (l)=csup*rdtemp(jlu) + cinf*rdtemp(jll)
  o3vmrc(l)=csup*rdo3v (jlu) + cinf*rdo3v (jll)
END DO

!     Work out geopotential height (km). Bottom level is ground.
pres(1)=1013.0
alt (1)=0.0
DO l=2,jplevp1
  alt(l)=alt(l-1) +0.001*ra*tempc(l-1)*LOG(pres(l-1)/pres(l))/gg
END DO
!     Height at centre of level and thickness of level (km).
DO l=1,jplev
  altc(l)=0.5*(alt(l+1) + alt(l))
  dalt(l)=     alt(l+1) - alt(l)
END DO

!     Work out O2 column above interfaces
DO l=1,jplevp1
  sao2(l)=r2o2*100.0*pres(l)*colfac
END DO

!     Average [O2] and [M] within a level
DO l=1,jplev
  do2c(l)=(sao2(l)-sao2(l+1))/(1.0e5*dalt(l))
  drsc(l)= do2c(l)/r2o2
END DO
!     Add O2 above top level to top level
do2c(jplev)=do2c(jplev) + sao2(jplevp1)/(1.0e5*dalt(jplev))

!     Average [O3] within a level
DO l=1,jplev
  do3c(l)=o3vmrc(l)*do2c(l)/r2o2
END DO

!     O3 above top level with scale height of 3.5 km
o3above=0.5*o3vmrc(jplev)*sao2(jplevp1)/r2o2
do3c (jplev)=do3c(jplev) + o3above/(1.0e5*dalt(jplev))

!     Column O3 above centre:
sao3c(jplev)=0.5*(o3above + 1.e5*do3c(jplev)*dalt(jplev))
DO l=jplev-1,1,-1
  sao3c(l)=sao3c(l+1) + 0.5e5*(dalt(l+1)*do3c(l+1)                             &
    +dalt(l  )*do3c(l  ))
END DO
! Call fill_spectra routine here to define spectral information etc.
CALL fill_spectra


! IF FAST-JX is selected, only do first 45 wavelength (up to 177 nm)
IF (l_ukca_fastjx) THEN
  DO jw = 1, jphi
    IF (wavenm(jw) < 177.) last_wn = jw
  END DO
ELSE
  last_wn = jphi
END IF

WRITE (6,'(1X,A,1X,I3,1X,I3,1X,E12.5)') 'last_wn = ',last_wn, jphi,            &
       wavenm(last_wn)

!     Calculate the Rayleigh scattering cross section.
CALL scatcs(jpwav,scs,wavenm)

!     Set zenith angle, max za, temperature and O3 grid.
CALL setzen(tabang,zenmax,altc,tabt,tabo3)

!     Loop over O3 profiles
DO jo=1,jpo3p

  !     Scale O3 profile
  DO jl=1,jplev
    do3cu(jl)=do3c(jl)*tabo3(jo)
  END DO

  CALL settab(alt,altc,do2c,do3cu,tempc,                                       &
    presc,drsc,tabs,tabang,albedo,zenmax,lscat,                                &
    tspo2,tabpres,quanta,wavenm,scs,ao2,ao2sr,ao3,dalt)

  !  Set up the photolysis rate tabulations. 39 rates.

  tabj2a   (:,:,:,jo)=0.0
  tabj2b   (:,:,:,jo)=0.0
  tabj3    (:,:,:,jo)=0.0
  tabj3a   (:,:,:,jo)=0.0
  tabjno2  (:,:,:,jo)=0.0
  tabjno   (:,:,:,jo)=0.0
  tabjno31 (:,:,:,jo)=0.0
  tabjno32 (:,:,:,jo)=0.0
  tabjn2o  (:,:,:,jo)=0.0
  tabjn2o5 (:,:,:,jo)=0.0
  tabjhno3 (:,:,:,jo)=0.0
  tabjpna  (:,:,:,jo)=0.0
  tabjcnita(:,:,:,jo)=0.0
  tabjcnitb(:,:,:,jo)=0.0
  tabjh2o2 (:,:,:,jo)=0.0
  tabjh2o  (:,:,:,jo)=0.0
  tabjhocl (:,:,:,jo)=0.0
  tabjhcl  (:,:,:,jo)=0.0
  tabjcl2o2(:,:,:,jo)=0.0
  tabjch3cl(:,:,:,jo)=0.0
  tabjccl4 (:,:,:,jo)=0.0
  tabjf11  (:,:,:,jo)=0.0
  tabjf12  (:,:,:,jo)=0.0
  tabjf22  (:,:,:,jo)=0.0
  tabjf113 (:,:,:,jo)=0.0
  tabjbro  (:,:,:,jo)=0.0
  tabjhobr (:,:,:,jo)=0.0
  tabjbrno3(:,:,:,jo)=0.0
  tabjbrcl (:,:,:,jo)=0.0
  tabjoclo (:,:,:,jo)=0.0
  tabjc2oa (:,:,:,jo)=0.0
  tabjc2ob (:,:,:,jo)=0.0
  tabjmhp  (:,:,:,jo)=0.0
  tabjmcfm (:,:,:,jo)=0.0
  tabjch3br(:,:,:,jo)=0.0
  tabjf12b1(:,:,:,jo)=0.0
  tabjf13b1(:,:,:,jo)=0.0
  tabjcof2 (:,:,:,jo)=0.0
  tabjcofcl(:,:,:,jo)=0.0
  tabjch4  (:,:,:,jo)=0.0
  tabjco2  (:,:,:,jo)=0.0
  tabjcos  (:,:,:,jo)=0.0
  tabjcos  (:,:,:,jo)=0.0
  tabjhono (:,:,:,jo)=0.0
  tabjmena (:,:,:,jo)=0.0
  tabjchbr3(:,:,:,jo)=0.0
  tabjdbrm (:,:,:,jo)=0.0
  tabjcs2  (:,:,:,jo)=0.0
  tabjh2so4(:,:,:,jo)=0.0
  tabjso3  (:,:,:,jo)=0.0

  !     Lyman alpha photolysis
  CALL lymana(tspo2,tabj2b,tabjh2o,tabjch4,quanta,jo)

  !     Temperature loop.
  DO jt = 1,jptem

    !        Set up the cross-sections for this pressure & temperature.
    !        N2O5
    CALL acsn2o5(tabt(jt),jpwav,wavenm,an2o5)

    !        N2O
    CALL acsn2o (tabt(jt),jpwav,wavenm,an2o)

    !        CCl4
    CALL acsccl4(tabt(jt),jpwav,wavenm,accl4)

    !        F11 (CFCl3)
    CALL acsf11 (tabt(jt),wavenm,af11)

    !        F12 (CF2Cl2)
    CALL acsf12 (tabt(jt),wavenm,af12)

    !        F22 (CHF2Cl)
    CALL acsf22 (tabt(jt),jpwav,wavenm,af22)

    !        CH3Cl
    CALL acsmc  (tabt(jt),jpwav,wavenm,ach3cl)

    !        O[1D] Quantum yield.
    CALL quanto12(tabt(jt),jpwav,wavenm,qeo1d)

    !        O3
    CALL acso3  (tabt(jt),jpwav,ao3)

    !        HNO3
    CALL acshno3(tabt(jt),wavenm,ahno3)

    !        ClONO2
    CALL acscnit(tabt(jt),wavenm,acnita, acnitb)

    !        H2O2
    CALL acsh2o2(tabt(jt),jpwav,wavenm,ah2o2)

    !        BrCl
    CALL acsbrcl(tabt(jt),jpwav,wavenm,abrcl)

    !        NO2
    CALL acsno2 (tabt(jt),wavenm,ano2 )

    !        COS
    CALL acscos (tabt(jt),wavenm,aocs )

    !        MeONO2
    CALL acsmena(tabt(jt),wavenm,amena )

    !        CHBr3
    ! Do not calculate CHBr3. Not needed for CCMVal
    !         CALL ACSCHBR3(TABT(JT),WAVENM,ACHBR3 )

    !        CH2Br2
    CALL acsdbrm(tabt(jt),adbrm )

    !        CS2
    CALL acscs2 (tabt(jt),wavenm,acs2 )
    !        H2SO4
    CALL acsh2so4 (tabt(jt),wavenm,ah2so4 )
    !        SO3
    CALL acsso3 (tabt(jt),wavenm,aso3 )

    !        Level loop.
    DO jl = 1,jplev

      !        Zenith angle loop.
      DO jc = 1 , jpchi

        !          Calculate the NO pressure dependent cross section.
        CALL acsno(tabang(jc),presc(jl),tspo2(jl,1),ano)

        !          Schumann-Runge bands.
        IF ( tabang(jc) <= zenmax(jl) ) THEN
          CALL acssr(tspo2(jl,jc),jpwav,ao2sr)
        ELSE
          DO jw = jplo, jphi
            ao2sr(jw) = 0.0
          END DO
        END IF

        !          Set up enhancement factor array.
        DO jw = jplo , jphi
          factor(jw) = MAX(0.0,tabs(jl,jc,jw))
        END DO

        !          Wavelength loop.
        DO jw = jplo , last_wn

          !           Number of photons into the volume element.
          !           (need to attenuate quanta above top level?)
          qu = quanta(jw)*factor(jw)

          tabj2a (jl,jc,jt,jo)=tabj2a(jl,jc,jt,jo)+qu*(ao2a(jw)+ao2sr(jw))
          tabj2b (jl,jc,jt,jo)=tabj2b(jl,jc,jt,jo)+qu*ao2b(jw)
          tabj3  (jl,jc,jt,jo)=tabj3(jl,jc,jt,jo)+qu*ao3(jw)*(1.0-qeo1d(jw))
          tabj3a (jl,jc,jt,jo)=tabj3a (jl,jc,jt,jo)+qu*ao3 (jw)*qeo1d(jw)
          tabjno2(jl,jc,jt,jo)=tabjno2(jl,jc,jt,jo)+qu*ano2(jw)*qeno2(jw)
          tabjno   (jl,jc,jt,jo)=tabjno   (jl,jc,jt,jo) + qu*ano   (jw)
          tabjno31 (jl,jc,jt,jo)=tabjno31 (jl,jc,jt,jo) + qu*ano31 (jw)
          tabjno32 (jl,jc,jt,jo)=tabjno32 (jl,jc,jt,jo) + qu*ano32 (jw)
          tabjn2o  (jl,jc,jt,jo)=tabjn2o  (jl,jc,jt,jo) + qu*an2o  (jw)
          tabjn2o5 (jl,jc,jt,jo)=tabjn2o5 (jl,jc,jt,jo) + qu*an2o5 (jw)
          tabjhno3 (jl,jc,jt,jo)=tabjhno3 (jl,jc,jt,jo) + qu*ahno3 (jw)
          tabjpna  (jl,jc,jt,jo)=tabjpna  (jl,jc,jt,jo) + qu*apna  (jw)
          tabjcnita(jl,jc,jt,jo)=tabjcnita(jl,jc,jt,jo) + qu*acnita(jw)
          tabjcnitb(jl,jc,jt,jo)=tabjcnitb(jl,jc,jt,jo) + qu*acnitb(jw)
          tabjh2o2 (jl,jc,jt,jo)=tabjh2o2 (jl,jc,jt,jo) + qu*ah2o2 (jw)
          tabjh2o  (jl,jc,jt,jo)=tabjh2o  (jl,jc,jt,jo) + qu*ah2o  (jw)
          tabjhocl (jl,jc,jt,jo)=tabjhocl (jl,jc,jt,jo) + qu*ahocl (jw)
          tabjhcl  (jl,jc,jt,jo)=tabjhcl  (jl,jc,jt,jo) + qu*ahcl  (jw)
          tabjcl2o2(jl,jc,jt,jo)=tabjcl2o2(jl,jc,jt,jo) + qu*acl2o2(jw)
          tabjch3cl(jl,jc,jt,jo)=tabjch3cl(jl,jc,jt,jo) + qu*ach3cl(jw)
          tabjccl4 (jl,jc,jt,jo)=tabjccl4 (jl,jc,jt,jo) + qu*accl4 (jw)
          tabjf11  (jl,jc,jt,jo)=tabjf11  (jl,jc,jt,jo) + qu*af11  (jw)
          tabjf12  (jl,jc,jt,jo)=tabjf12  (jl,jc,jt,jo) + qu*af12  (jw)
          tabjf22  (jl,jc,jt,jo)=tabjf22  (jl,jc,jt,jo) + qu*af22  (jw)
          tabjf113 (jl,jc,jt,jo)=tabjf113 (jl,jc,jt,jo) + qu*af113 (jw)
          tabjbro  (jl,jc,jt,jo)=tabjbro  (jl,jc,jt,jo) + qu*abro  (jw)
          tabjhobr (jl,jc,jt,jo)=tabjhobr (jl,jc,jt,jo) + qu*ahobr (jw)
          tabjbrno3(jl,jc,jt,jo)=tabjbrno3(jl,jc,jt,jo) + qu*abrno3(jw)
          tabjbrcl (jl,jc,jt,jo)=tabjbrcl (jl,jc,jt,jo) + qu*abrcl (jw)
          tabjoclo (jl,jc,jt,jo)=tabjoclo (jl,jc,jt,jo) + qu*aoclo (jw)
          tabjc2oa (jl,jc,jt,jo)=tabjc2oa (jl,jc,jt,jo) + qu*ac2oa (jw)
          tabjc2ob (jl,jc,jt,jo)=tabjc2ob (jl,jc,jt,jo) + qu*ac2ob (jw)
          tabjmhp  (jl,jc,jt,jo)=tabjmhp  (jl,jc,jt,jo) + qu*amhp  (jw)
          tabjch3br(jl,jc,jt,jo)=tabjch3br(jl,jc,jt,jo) + qu*ach3br(jw)
          tabjf12b1(jl,jc,jt,jo)=tabjf12b1(jl,jc,jt,jo) + qu*af12b1(jw)
          tabjf13b1(jl,jc,jt,jo)=tabjf13b1(jl,jc,jt,jo) + qu*af13b1(jw)
          tabjcof2 (jl,jc,jt,jo)=tabjcof2 (jl,jc,jt,jo) + qu*acof2 (jw)
          tabjcofcl(jl,jc,jt,jo)=tabjcofcl(jl,jc,jt,jo) + qu*acofcl(jw)
          tabjmcfm (jl,jc,jt,jo)=tabjmcfm (jl,jc,jt,jo) + qu*amcfm (jw)
          tabjch4  (jl,jc,jt,jo)=tabjch4  (jl,jc,jt,jo) + qu*ach4  (jw)
          tabjco2  (jl,jc,jt,jo)=tabjco2  (jl,jc,jt,jo) + qu*aco2  (jw)
          tabjcos  (jl,jc,jt,jo)=tabjcos  (jl,jc,jt,jo) + qu*aocs  (jw)
          tabjhono (jl,jc,jt,jo)=tabjhono (jl,jc,jt,jo) + qu*ahono (jw)
          tabjmena (jl,jc,jt,jo)=tabjmena (jl,jc,jt,jo) + qu*amena (jw)
          tabjchbr3(jl,jc,jt,jo)=tabjchbr3(jl,jc,jt,jo) + qu*achbr3(jw)
          tabjdbrm (jl,jc,jt,jo)=tabjdbrm (jl,jc,jt,jo) + qu*adbrm (jw)
          tabjcs2  (jl,jc,jt,jo)=tabjcs2  (jl,jc,jt,jo) + qu*acs2  (jw)
          tabjh2so4(jl,jc,jt,jo)=tabjh2so4(jl,jc,jt,jo) + qu*ah2so4(jw)
          tabjso3  (jl,jc,jt,jo)=tabjso3  (jl,jc,jt,jo) + qu*aso3  (jw)

        END DO


        !        End of the zenith angle loop.
      END DO

      !     End of level loop.
    END DO

    !     End of temperature loop.
  END DO

  !     End of O3 profile loop.
END DO
IF (lhook) CALL dr_hook('INIJTAB',zhook_out,zhook_handle)
RETURN



END SUBROUTINE inijtab
END MODULE inijtab_mod
