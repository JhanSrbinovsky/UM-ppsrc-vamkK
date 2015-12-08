! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE calcjs_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Look up photolysis rates from  j tables

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
SUBROUTINE calcjs(ihmin,ihmax,f,ppa,t,o3col,lnt)

USE ukca_dissoc, ONLY: ajhno3, ajpna,  ajh2o2, aj2a,  aj2b, aj3, &
       aj3a, ajcnita,ajcnitb,ajbrno3,ajbrcl, ajoclo, ajcl2o2,    &
       ajhocl, ajno,  ajno2,  ajn2o5, ajno31, ajno32, ajbro,     &
       ajhcl,  ajn2o,  ajhobr, ajf11,  ajf12,  ajh2o,  ajccl4,   &
       ajf113, ajf22,  ajch3cl, ajc2oa, ajc2ob, ajmhp,  ajch3br, &
       ajmcfm,ajch4,  ajf12b1,ajf13b1, ajcof2, ajcofcl, ajco2,   &
       ajcos,  ajhono, ajmena, ajchbr3, ajdbrm, ajcs2, ajh2so4,  &
       ajso3
USE printstatus_mod, ONLY: printstatus, prstatus_diag
USE ukca_constants, ONLY: pi
USE ukca_tbjs_mod,  ONLY:    tabj2a   , tabj2b   , tabj3    ,    &
       tabj3a   , tabjno   , tabjno31 , tabjno32 , tabjno2  ,    & 
       tabjn2o5 , tabjhno3 , tabjh2o  , tabjh2o2 , tabjf11  ,    &
       tabjf12  , tabjf22  , tabjf113 , tabjch3cl, tabjccl4 ,    &
       tabjcnita, tabjcnitb, tabjhcl  , tabjhocl , tabjpna  ,    &
       tabjcl2o2, tabjn2o  , tabjbrcl , tabjbrno3, tabjhobr ,    &
       tabjbro  , tabjoclo , tabjc2oa , tabjc2ob , tabjmhp  ,    &
       tabjmcfm , tabjch3br, tabjf12b1, tabjf13b1, tabjcof2 ,    &
       tabjcofcl, tabjch4  , tabjcos  , tabjso2  , tabjco2  ,    &
       tabjhono , tabjmena , tabjchbr3, tabjdbrm , tabjcs2  ,    &
       tabjh2so4, tabjso3  , tabpres,   tabang,    tabt     ,    &
       tabo3, sao3c
USE ukca_parpho_mod, ONLY: jplevp1, jplev, jpchi, jps90, jpchin,               &
                           jpwav, jplo, jphi, jptem, jpo3p, jps90,             &
                           szamax, tmin, tmax, o3min, o3max
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
INTEGER, INTENT(IN) :: ihmin        ! first point
INTEGER, INTENT(IN) :: ihmax        ! last point
INTEGER, INTENT(IN) :: lnt          ! No of points

REAL, INTENT(IN) :: f(lnt)          ! Cos of zenith angle
REAL, INTENT(IN) :: t(lnt)          ! Temperature
REAL, INTENT(IN) :: o3col(lnt)      ! Ozone column
REAL, INTENT(IN) :: ppa(lnt)        ! Pressure

! Local variables
LOGICAL, SAVE :: init = .FALSE.

INTEGER :: ih, ios

REAL :: csup
REAL :: dlnp
REAL :: tdif
REAL :: o3dif
REAL :: radint
REAL :: cosint

INTEGER :: jx(lnt)
INTEGER :: jxp1(lnt)
INTEGER :: jy(lnt)
INTEGER :: jyp1(lnt)
INTEGER :: jz(lnt)
INTEGER :: jzp1(lnt)
INTEGER :: jo(lnt)
INTEGER :: jop1(lnt)

REAL :: angle(lnt)
REAL :: o3fac(lnt)
REAL :: p(lnt)
REAL :: tu(lnt)
REAL :: zo (lnt)
REAL :: zx (lnt)
REAL :: zy (lnt)
REAL :: zz (lnt)
REAL :: co (lnt)
REAL :: cx (lnt)
REAL :: cy (lnt)
REAL :: cz (lnt)
REAL :: w1 (lnt)
REAL :: w2 (lnt)
REAL :: w3 (lnt)
REAL :: w4 (lnt)
REAL :: w5 (lnt)
REAL :: w6 (lnt)
REAL :: w7 (lnt)
REAL :: w8 (lnt)
REAL :: w9 (lnt)
REAL :: w10(lnt)
REAL :: w11(lnt)
REAL :: w12(lnt)
REAL :: w13(lnt)
REAL :: w14(lnt)
REAL :: w15(lnt)
REAL :: w16(lnt)

LOGICAL, SAVE :: jfirst=.TRUE.

CHARACTER(LEN=256) :: jtable

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! do not read JTABLE file. Instead, call inijtab routine.

IF (lhook) CALL dr_hook('CALCJS',zhook_in,zhook_handle)
IF (init) THEN

  !       Open file for the photolysis tables and and read them.
  CALL get_file(ukcastrd_unit, jtable, 256, ios)
  OPEN(ukcastrd_unit,file=jtable,form='UNFORMATTED',STATUS='OLD')

  REWIND(ukcastrd_unit)

  !     Read the uv/vis radiation field.
  READ(ukcastrd_unit) tabpres ,tabang   ,tabt     ,tabo3 ,sao3c,               &
    tabj2a  ,tabj2b   ,tabj3    ,tabjch4  ,                                    &
    tabj3a  ,tabjccl4 ,tabjcl2o2,tabjcnita,tabjcnitb,                          &
    tabjf11 ,                                                                  &
    tabjf113,tabjf12  ,tabjf22  ,tabjh2o  ,tabjh2o2 ,                          &
    tabjhcl ,tabjhno3 ,tabjhocl ,tabjch3cl,tabjn2o  ,                          &
    tabjn2o5,tabjno   ,tabjno2  ,tabjno31 ,tabjno32 ,                          &
    tabjpna ,tabjbrcl ,tabjbrno3,tabjbro  ,tabjmcfm ,                          &
    tabjhobr,tabjoclo ,tabjc2oa ,tabjc2ob ,tabjmhp  ,                          &
    tabjcof2,tabjch3br,tabjf12b1,tabjf13b1,tabjcofcl,                          &
    tabjco2, tabjcos  ,tabjcs2  ,tabjh2so4,tabjso3

  CLOSE(ukcastrd_unit)
  init = .FALSE.

END IF


IF (jfirst .AND. printstatus == prstatus_diag) THEN

  !       Write out photolysis tables in ASCII format for test purposes.

  WRITE(6,'(A)') 'TABPRES'
  WRITE(6,'(10(E11.3))') tabpres
  WRITE(6,'(A)') 'TABANG'
  WRITE(6,'(10(E11.3))') tabang
  WRITE(6,'(A)') 'TABT'
  WRITE(6,'(10(E11.3))') tabt
  WRITE(6,'(A)') 'TABO3'
  WRITE(6,'(10(E11.3))') tabo3
  WRITE(6,'(A)') 'SAO3C'
  WRITE(6,'(10(E11.3))') sao3c
  WRITE(6,'(A)') 'TABJ2A'
  WRITE(6,'(10(E11.3))') tabj2a
  WRITE(6,'(A)') 'TABJ2B'
  WRITE(6,'(10(E11.3))') tabj2b
  WRITE(6,'(A)') 'TABJ3'
  WRITE(6,'(10(E11.3))') tabj3
  WRITE(6,'(A)') 'TABJCH4'
  WRITE(6,'(10(E11.3))') tabjch4
  WRITE(6,'(A)') 'TABJ3A'
  WRITE(6,'(10(E11.3))') tabj3a
  WRITE(6,'(A)') 'TABCCL4'
  WRITE(6,'(10(E11.3))') tabjccl4
  WRITE(6,'(A)') 'TABCl2O2'
  WRITE(6,'(10(E11.3))') tabjcl2o2
  WRITE(6,'(A)') 'TABCNITA'
  WRITE(6,'(10(E11.3))') tabjcnita
  WRITE(6,'(A)') 'TABCNITB'
  WRITE(6,'(10(E11.3))') tabjcnitb
  WRITE(6,'(A)') 'TABJF11'
  WRITE(6,'(10(E11.3))') tabjf11
  WRITE(6,'(A)') 'TABJF113'
  WRITE(6,'(10(E11.3))') tabjf113
  WRITE(6,'(A)') 'TABJF12'
  WRITE(6,'(10(E11.3))') tabjf12
  WRITE(6,'(A)') 'TABJF22'
  WRITE(6,'(10(E11.3))') tabjf22
  WRITE(6,'(A)') 'TABJH2O'
  WRITE(6,'(10(E11.3))') tabjh2o
  WRITE(6,'(A)') 'TABJH2O2'
  WRITE(6,'(10(E11.3))') tabjh2o2
  WRITE(6,'(A)') 'TABHCL'
  WRITE(6,'(10(E11.3))') tabjhcl
  WRITE(6,'(A)') 'TABHNO3'
  WRITE(6,'(10(E11.3))') tabjhno3
  WRITE(6,'(A)') 'TABJHOCL'
  WRITE(6,'(10(E11.3))') tabjhocl
  WRITE(6,'(A)') 'TABJCH3CL'
  WRITE(6,'(10(E11.3))') tabjch3cl
  WRITE(6,'(A)') 'TABJN2O'
  WRITE(6,'(10(E11.3))') tabjn2o
  WRITE(6,'(A)') 'TABJN2O5'
  WRITE(6,'(10(E11.3))') tabjn2o5
  WRITE(6,'(A)') 'TABJNO'
  WRITE(6,'(10(E11.3))') tabjno
  WRITE(6,'(A)') 'TABJNO2'
  WRITE(6,'(10(E11.3))') tabjno2
  WRITE(6,'(A)') 'TABJNO31'
  WRITE(6,'(10(E11.3))') tabjno31
  WRITE(6,'(A)') 'TABJNO32'
  WRITE(6,'(10(E11.3))') tabjno32
  WRITE(6,'(A)') 'TABJPNA'
  WRITE(6,'(10(E11.3))') tabjpna
  WRITE(6,'(A)') 'TABJBRCL'
  WRITE(6,'(10(E11.3))') tabjbrcl
  WRITE(6,'(A)') 'TABJBRNO3'
  WRITE(6,'(10(E11.3))') tabjbrno3
  WRITE(6,'(A)') 'TABJBRO'
  WRITE(6,'(10(E11.3))') tabjbro
  WRITE(6,'(A)') 'TABJMCFM'
  WRITE(6,'(10(E11.3))') tabjmcfm
  WRITE(6,'(A)') 'TABJHOBR'
  WRITE(6,'(10(E11.3))') tabjhobr
  WRITE(6,'(A)') 'TABJOCLO'
  WRITE(6,'(10(E11.3))') tabjoclo
  WRITE(6,'(A)') 'TABJC2OA'
  WRITE(6,'(10(E11.3))') tabjc2oa
  WRITE(6,'(A)') 'TABJC2OB'
  WRITE(6,'(10(E11.3))') tabjc2ob
  WRITE(6,'(A)') 'TABJMHP'
  WRITE(6,'(10(E11.3))') tabjmhp
  WRITE(6,'(A)') 'TABCOF2'
  WRITE(6,'(10(E11.3))') tabjcof2
  WRITE(6,'(A)') 'TABJCH3BR'
  WRITE(6,'(10(E11.3))') tabjch3br
  WRITE(6,'(A)') 'TABF12B1'
  WRITE(6,'(10(E11.3))') tabjf12b1
  WRITE(6,'(A)') 'TABF13B1'
  WRITE(6,'(10(E11.3))') tabjf13b1
  WRITE(6,'(A)') 'TABJCOFCL'
  WRITE(6,'(10(E11.3))') tabjcofcl
  WRITE(6,'(A)') 'TABJCO2'
  WRITE(6,'(10(E11.3))') tabjco2
  WRITE(6,'(A)') 'TABJCOS'
  WRITE(6,'(10(E11.3))') tabjcos
  WRITE(6,'(A)') 'TABJCS2'
  WRITE(6,'(10(E11.3))') tabjcs2
  WRITE(6,'(A)') 'TABJH2SO4'
  WRITE(6,'(10(E11.3))') tabjh2so4
  WRITE(6,'(A)') 'TABJSO3'
  WRITE(6,'(10(E11.3))') tabjso3

  jfirst=.FALSE.
END IF


!     Spacing of angles in lookup table
cosint=1.0/REAL(jpchi-1-jps90)
radint=(szamax - 90.0)*pi/(180.0*REAL(jps90))

!     Model pressure levels in hPa. SZA in radians. T within range.
DO ih=ihmin,ihmax
  p (ih)=ppa(ih)*0.01
  p (ih)=MIN(p (ih), tabpres(    1))
  p (ih)=MAX(p (ih), tabpres(jplev))
  tu(ih)=t(ih)
  tu(ih)=MIN(tu(ih), tmax)
  tu(ih)=MAX(tu(ih), tmin)
  angle(ih)=ACOS(f(ih))
END DO

!     Find the location of the pressure levels,
!     zenith angle, temperature and O3 in the photolysis look up tables.

!     i)  X. Pressure levels equally spaced in log(p)
!     ii) Y. Zenith angle
dlnp=(LOG(tabpres(1))-LOG(tabpres(jplev)))/REAL(jplev-1)

DO ih = ihmin, ihmax
  jx(ih)=INT((LOG(tabpres(1))-LOG(p(ih)))/dlnp) + 1

  IF (angle(ih) < tabang(jpchin)) THEN
    jy(ih)=INT((1.0 - f(ih))/cosint) + 1
  ELSE
    jy(ih)=jpchin + INT((angle(ih)-0.5*pi)/radint)
  END IF
END DO

!     iii) Z  temperature evenly spaced
IF (jptem == 1) THEN
  DO ih = ihmin, ihmax
    jz(ih)=1
  END DO
ELSE
  tdif=(tmax-tmin)/REAL(jptem-1)
  DO ih = ihmin, ihmax
    jz(ih)=INT((tu(ih)-tmin)/tdif) + 1
  END DO
END IF

!     iv) O3 factor
IF (jpo3p == 1) THEN
  DO ih = ihmin, ihmax
    jo(ih)=1
  END DO
ELSE
  !       Find O3 factor using tabulated O3 column
  o3dif=(o3max-o3min)/REAL(jpo3p-1)
  DO ih=ihmin,ihmax
    csup=(p(ih)            -tabpres(jx(ih)))/                                  &
      (tabpres(jx(ih)+1)-tabpres(jx(ih)))

    !         Ratio O3 to standard atmosphere O3 above model pressure.
    o3fac(ih)=o3col(ih)/                                                       &
      ((1.0-csup)*sao3c(jx(ih))+csup*sao3c(jx(ih)+1))

    !         Check that O3 column ratio is within range
    o3fac(ih)=MIN(o3fac(ih), o3max)
    o3fac(ih)=MAX(o3fac(ih), o3min)

    jo(ih)=INT((o3fac(ih)-o3min)/o3dif) + 1
  END DO
END IF

!     If out of range set pointers accordingly.
DO ih=ihmin,ihmax
  jx(ih) = max0(jx(ih),      1)
  jx(ih) = min0(jx(ih),jplev-1)
  jy(ih) = max0(jy(ih),      1)
  jy(ih) = min0(jy(ih),jpchi-1)
  jz(ih) = max0(jz(ih),      1)
  jz(ih) = min0(jz(ih),jptem-1)
  jo(ih) = max0(jo(ih),      1)
  jo(ih) = min0(jo(ih),jpo3p-1)

  jz(ih) = max0(jz(ih),1)
  jo(ih) = max0(jo(ih),1)

  jxp1(ih) =      jx(ih) + 1
  jyp1(ih) =      jy(ih) + 1
  jzp1(ih) = min0(jz(ih) + 1, jptem)
  jop1(ih) = min0(jo(ih) + 1, jpo3p)
END DO

!     Find points used in interpolation
DO ih=ihmin,ihmax
  !        i) X: pressure
  zx(ih)=(LOG(p      (ih))       - LOG(tabpres(jx(ih))))/                      &
    (LOG(tabpres(jxp1(ih))) - LOG(tabpres(jx(ih))))

  !        ii) Y: zenith angle
  zy(ih)=(angle(ih)         - tabang (jy(ih)))/                                &
    (tabang (jyp1(ih)) - tabang (jy(ih)))

  !        iii) Z: temperature
  IF (jptem == 1) THEN
    zz(ih)=1.0
  ELSE
    zz(ih)=(tu   (ih)       - tabt (jz(ih)))/                                  &
      (tabt (jzp1(ih)) - tabt (jz(ih)))
  END IF

  !        iv) O: O3 profile
  IF (jpo3p == 1) THEN
    zo(ih)=1.0
  ELSE
    zo(ih)=(o3fac(ih)       - tabo3(jo(ih)))/                                  &
      (tabo3(jop1(ih)) - tabo3(jo(ih)))
  END IF

  cx(ih)=1.0 - zx(ih)
  cy(ih)=1.0 - zy(ih)
  cz(ih)=1.0 - zz(ih)
  co(ih)=1.0 - zo(ih)
END DO

!     Calculate weights for interpolation
DO ih=ihmin,ihmax
  w1 (ih)=cx(ih)*cy(ih)*cz(ih)*co(ih)
  w2 (ih)=zx(ih)*cy(ih)*cz(ih)*co(ih)
  w3 (ih)=zx(ih)*zy(ih)*cz(ih)*co(ih)
  w4 (ih)=cx(ih)*zy(ih)*cz(ih)*co(ih)
  w5 (ih)=cx(ih)*cy(ih)*zz(ih)*co(ih)
  w6 (ih)=zx(ih)*cy(ih)*zz(ih)*co(ih)
  w7 (ih)=zx(ih)*zy(ih)*zz(ih)*co(ih)
  w8 (ih)=cx(ih)*zy(ih)*zz(ih)*co(ih)
  w9 (ih)=cx(ih)*cy(ih)*cz(ih)*zo(ih)
  w10(ih)=zx(ih)*cy(ih)*cz(ih)*zo(ih)
  w11(ih)=zx(ih)*zy(ih)*cz(ih)*zo(ih)
  w12(ih)=cx(ih)*zy(ih)*cz(ih)*zo(ih)
  w13(ih)=cx(ih)*cy(ih)*zz(ih)*zo(ih)
  w14(ih)=zx(ih)*cy(ih)*zz(ih)*zo(ih)
  w15(ih)=zx(ih)*zy(ih)*zz(ih)*zo(ih)
  w16(ih)=cx(ih)*zy(ih)*zz(ih)*zo(ih)
END DO

!     Look up photolysis rates. 38 rates.
DO ih=ihmin,ihmax

  !     If angle in range
  IF (angle(ih) <= tabang(jpchi)) THEN
    !1
    aj2a  (ih)=                                                                &
      w1 (ih)*tabj2a   (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabj2a   (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabj2a   (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabj2a   (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabj2a   (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabj2a   (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabj2a   (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabj2a   (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabj2a   (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabj2a   (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabj2a   (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabj2a   (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabj2a   (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabj2a   (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabj2a   (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabj2a   (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !1B
    aj2b  (ih)=                                                                &
      w1 (ih)*tabj2b   (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabj2b   (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabj2b   (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabj2b   (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabj2b   (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabj2b   (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabj2b   (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabj2b   (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabj2b   (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabj2b   (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabj2b   (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabj2b   (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabj2b   (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabj2b   (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabj2b   (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabj2b   (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !2
    aj3   (ih)=                                                                &
      w1 (ih)*tabj3    (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabj3    (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabj3    (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabj3    (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabj3    (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabj3    (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabj3    (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabj3    (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabj3    (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabj3    (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabj3    (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabj3    (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabj3    (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabj3    (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabj3    (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabj3    (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !3
    aj3a  (ih)=                                                                &
      w1 (ih)*tabj3a   (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabj3a   (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabj3a   (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabj3a   (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabj3a   (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabj3a   (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabj3a   (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabj3a   (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabj3a   (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabj3a   (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabj3a   (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabj3a   (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabj3a   (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabj3a   (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabj3a   (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabj3a   (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !4
    ajno  (ih)=                                                                &
      w1 (ih)*tabjno   (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjno   (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjno   (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjno   (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjno   (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjno   (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjno   (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjno   (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjno   (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjno   (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjno   (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjno   (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjno   (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjno   (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjno   (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjno   (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !5
    ajno2 (ih)=                                                                &
      w1 (ih)*tabjno2  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjno2  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjno2  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjno2  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjno2  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjno2  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjno2  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjno2  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjno2  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjno2  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjno2  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjno2  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjno2  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjno2  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjno2  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjno2  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !6
    ajno31(ih)=                                                                &
      w1 (ih)*tabjno31 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjno31 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjno31 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjno31 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjno31 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjno31 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjno31 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjno31 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjno31 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjno31 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjno31 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjno31 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjno31 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjno31 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjno31 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjno31 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !7
    ajno32(ih)=                                                                &
      w1 (ih)*tabjno32 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjno32 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjno32 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjno32 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjno32 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjno32 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjno32 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjno32 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjno32 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjno32 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjno32 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjno32 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjno32 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjno32 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjno32 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjno32 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !8
    ajn2o (ih)=                                                                &
      w1 (ih)*tabjn2o  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjn2o  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjn2o  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjn2o  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjn2o  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjn2o  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjn2o  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjn2o  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjn2o  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjn2o  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjn2o  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjn2o  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjn2o  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjn2o  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjn2o  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjn2o  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !9
    ajn2o5(ih)=                                                                &
      w1 (ih)*tabjn2o5 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjn2o5 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjn2o5 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjn2o5 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjn2o5 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjn2o5 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjn2o5 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjn2o5 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjn2o5 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjn2o5 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjn2o5 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjn2o5 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjn2o5 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjn2o5 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjn2o5 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjn2o5 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !10
    ajhno3(ih)=                                                                &
      w1 (ih)*tabjhno3 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjhno3 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjhno3 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjhno3 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjhno3 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjhno3 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjhno3 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjhno3 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjhno3 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjhno3 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjhno3 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjhno3 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjhno3 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjhno3 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjhno3 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjhno3 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !11
    ajcnita(ih)=                                                               &
      w1 (ih)*tabjcnita(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjcnita(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjcnita(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjcnita(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjcnita(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjcnita(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjcnita(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjcnita(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjcnita(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjcnita(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjcnita(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjcnita(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjcnita(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjcnita(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjcnita(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjcnita(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    ajcnitb(ih)=                                                               &
      w1 (ih)*tabjcnitb(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjcnitb(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjcnitb(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjcnitb(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjcnitb(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjcnitb(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjcnitb(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjcnitb(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjcnitb(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjcnitb(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjcnitb(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjcnitb(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjcnitb(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjcnitb(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjcnitb(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjcnitb(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !12
    ajpna (ih)=                                                                &
      w1 (ih)*tabjpna  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjpna  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjpna  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjpna  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjpna  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjpna  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjpna  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjpna  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjpna  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjpna  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjpna  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjpna  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjpna  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjpna  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjpna  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjpna  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !13
    ajh2o2(ih)=                                                                &
      w1 (ih)*tabjh2o2 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjh2o2 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjh2o2 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjh2o2 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjh2o2 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjh2o2 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjh2o2 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjh2o2 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjh2o2 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjh2o2 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjh2o2 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjh2o2 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjh2o2 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjh2o2 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjh2o2 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjh2o2 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !14
    ajh2o (ih)=                                                                &
      w1 (ih)*tabjh2o  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjh2o  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjh2o  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjh2o  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjh2o  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjh2o  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjh2o  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjh2o  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjh2o  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjh2o  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjh2o  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjh2o  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjh2o  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjh2o  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjh2o  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjh2o  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !15
    ajhocl(ih)=                                                                &
      w1 (ih)*tabjhocl (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjhocl (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjhocl (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjhocl (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjhocl (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjhocl (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjhocl (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjhocl (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjhocl (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjhocl (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjhocl (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjhocl (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjhocl (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjhocl (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjhocl (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjhocl (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !16
    ajhcl (ih)=                                                                &
      w1 (ih)*tabjhcl  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjhcl  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjhcl  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjhcl  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjhcl  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjhcl  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjhcl  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjhcl  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjhcl  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjhcl  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjhcl  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjhcl  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjhcl  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjhcl  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjhcl  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjhcl  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !17
    ajcl2o2(ih)=                                                               &
      w1 (ih)*tabjcl2o2(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjcl2o2(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjcl2o2(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjcl2o2(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjcl2o2(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjcl2o2(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjcl2o2(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjcl2o2(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjcl2o2(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjcl2o2(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjcl2o2(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjcl2o2(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjcl2o2(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjcl2o2(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjcl2o2(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjcl2o2(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !18
    ajbro (ih)=                                                                &
      w1 (ih)*tabjbro  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjbro  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjbro  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjbro  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjbro  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjbro  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjbro  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjbro  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjbro  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjbro  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjbro  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjbro  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjbro  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjbro  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjbro  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjbro  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !19
    ajbrcl(ih)=                                                                &
      w1 (ih)*tabjbrcl (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjbrcl (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjbrcl (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjbrcl (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjbrcl (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjbrcl (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjbrcl (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjbrcl (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjbrcl (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjbrcl (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjbrcl (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjbrcl (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjbrcl (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjbrcl (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjbrcl (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjbrcl (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !20
    ajco2(ih)=                                                                 &
      w1 (ih)*tabjco2  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjco2  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjco2  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjco2  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjco2  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjco2  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjco2  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjco2  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjco2  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjco2  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjco2  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjco2  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjco2  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjco2  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjco2  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjco2  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    ! COS
    ajcos(ih)=                                                                 &
      w1 (ih)*tabjcos  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjcos  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjcos  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjcos  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjcos  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjcos  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjcos  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjcos  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjcos  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjcos  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjcos  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjcos  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjcos  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjcos  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjcos  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjcos  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    ! HONO
    ajhono(ih) =                                                               &
      w1 (ih)*tabjhono (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjhono (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjhono (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjhono (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjhono (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjhono (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjhono (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjhono (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjhono (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjhono (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjhono (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjhono (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjhono (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjhono (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjhono (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjhono (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    ! MeONO2
    ajmena(ih) =                                                               &
      w1 (ih)*tabjmena (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjmena (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjmena (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjmena (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjmena (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjmena (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjmena (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjmena (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjmena (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjmena (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjmena (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjmena (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjmena (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjmena (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjmena (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjmena (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    ! CHBr3
    ajchbr3(ih) =                                                              &
      w1 (ih)*tabjchbr3(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjchbr3(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjchbr3(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjchbr3(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjchbr3(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjchbr3(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjchbr3(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjchbr3(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjchbr3(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjchbr3(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjchbr3(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjchbr3(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjchbr3(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjchbr3(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjchbr3(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjchbr3(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    ! CH2Br2
    ajdbrm (ih) =                                                              &
      w1 (ih)*tabjdbrm (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjdbrm (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjdbrm (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjdbrm (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjdbrm (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjdbrm (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjdbrm (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjdbrm (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjdbrm (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjdbrm (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjdbrm (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjdbrm (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjdbrm (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjdbrm (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjdbrm (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjdbrm (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    ! CS2
    ajcs2 (ih) =                                                               &
      w1 (ih)*tabjcs2 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                    &
      + w2 (ih)*tabjcs2 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                  &
      + w3 (ih)*tabjcs2 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                  &
      + w4 (ih)*tabjcs2 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                  &
      + w5 (ih)*tabjcs2 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                  &
      + w6 (ih)*tabjcs2 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                  &
      + w7 (ih)*tabjcs2 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                  &
      + w8 (ih)*tabjcs2 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                  &
      + w9 (ih)*tabjcs2 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                  &
      + w10(ih)*tabjcs2 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                  &
      + w11(ih)*tabjcs2 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                  &
      + w12(ih)*tabjcs2 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                  &
      + w13(ih)*tabjcs2 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                  &
      + w14(ih)*tabjcs2 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                  &
      + w15(ih)*tabjcs2 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                  &
      + w16(ih)*tabjcs2 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    ! H2SO4
    ajh2so4 (ih) =                                                             &
      w1 (ih)*tabjh2so4 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                  &
      + w2 (ih)*tabjh2so4 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                &
      + w3 (ih)*tabjh2so4 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                &
      + w4 (ih)*tabjh2so4 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                &
      + w5 (ih)*tabjh2so4 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                &
      + w6 (ih)*tabjh2so4 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                &
      + w7 (ih)*tabjh2so4 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                &
      + w8 (ih)*tabjh2so4 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                &
      + w9 (ih)*tabjh2so4 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                &
      + w10(ih)*tabjh2so4 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                &
      + w11(ih)*tabjh2so4 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                &
      + w12(ih)*tabjh2so4 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                &
      + w13(ih)*tabjh2so4 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                &
      + w14(ih)*tabjh2so4 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                &
      + w15(ih)*tabjh2so4 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                &
      + w16(ih)*tabjh2so4 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    ! SO3
    ajso3 (ih) =                                                               &
      w1 (ih)*tabjso3 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                    &
      + w2 (ih)*tabjso3 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                  &
      + w3 (ih)*tabjso3 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                  &
      + w4 (ih)*tabjso3 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                  &
      + w5 (ih)*tabjso3 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                  &
      + w6 (ih)*tabjso3 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                  &
      + w7 (ih)*tabjso3 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                  &
      + w8 (ih)*tabjso3 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                  &
      + w9 (ih)*tabjso3 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                  &
      + w10(ih)*tabjso3 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                  &
      + w11(ih)*tabjso3 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                  &
      + w12(ih)*tabjso3 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                  &
      + w13(ih)*tabjso3 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                  &
      + w14(ih)*tabjso3 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                  &
      + w15(ih)*tabjso3 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                  &
      + w16(ih)*tabjso3 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))

    !     Else darkness
  ELSE
    aj2a   (ih) = 0.0
    aj2b   (ih) = 0.0
    aj3    (ih) = 0.0
    aj3a   (ih) = 0.0
    ajno   (ih) = 0.0
    ajno2  (ih) = 0.0
    ajno31 (ih) = 0.0
    ajno32 (ih) = 0.0
    ajn2o  (ih) = 0.0
    ajn2o5 (ih) = 0.0
    ajhno3 (ih) = 0.0
    ajcnita(ih) = 0.0
    ajcnitb(ih) = 0.0
    ajpna  (ih) = 0.0
    ajh2o2 (ih) = 0.0
    ajh2o  (ih) = 0.0
    ajhocl (ih) = 0.0
    ajhcl  (ih) = 0.0
    ajcl2o2(ih) = 0.0
    ajbro  (ih) = 0.0
    ajbrcl (ih) = 0.0
    ajco2  (ih) = 0.0
    ajcos  (ih) = 0.0
    ajhono (ih) = 0.0
    ajmena (ih) = 0.0
    ajchbr3(ih) = 0.0
    ajdbrm (ih) = 0.0
    ajcs2  (ih) = 0.0
    ajh2so4(ih) = 0.0
    ajso3  (ih) = 0.0
  END IF

END DO

DO ih=ihmin,ihmax

  !     If angle in range
  IF (angle(ih) <= tabang(jpchi)) THEN
    !20
    ajbrno3(ih)=                                                               &
      w1 (ih)*tabjbrno3(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjbrno3(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjbrno3(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjbrno3(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjbrno3(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjbrno3(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjbrno3(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjbrno3(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjbrno3(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjbrno3(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjbrno3(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjbrno3(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjbrno3(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjbrno3(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjbrno3(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjbrno3(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !21
    ajhobr(ih)=                                                                &
      w1 (ih)*tabjhobr (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjhobr (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjhobr (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjhobr (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjhobr (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjhobr (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjhobr (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjhobr (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjhobr (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjhobr (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjhobr (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjhobr (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjhobr (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjhobr (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjhobr (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjhobr (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !22
    ajoclo(ih)=                                                                &
      w1 (ih)*tabjoclo (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjoclo (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjoclo (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjoclo (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjoclo (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjoclo (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjoclo (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjoclo (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjoclo (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjoclo (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjoclo (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjoclo (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjoclo (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjoclo (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjoclo (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjoclo (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !23
    ajc2oa(ih)=                                                                &
      w1 (ih)*tabjc2oa (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjc2oa (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjc2oa (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjc2oa (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjc2oa (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjc2oa (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjc2oa (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjc2oa (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjc2oa (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjc2oa (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjc2oa (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjc2oa (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjc2oa (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjc2oa (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjc2oa (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjc2oa (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !24
    ajc2ob(ih)=                                                                &
      w1 (ih)*tabjc2ob (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjc2ob (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjc2ob (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjc2ob (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjc2ob (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjc2ob (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjc2ob (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjc2ob (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjc2ob (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjc2ob (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjc2ob (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjc2ob (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjc2ob (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjc2ob (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjc2ob (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjc2ob (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !25
    ajmhp (ih)=                                                                &
      w1 (ih)*tabjmhp  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjmhp  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjmhp  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjmhp  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjmhp  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjmhp  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjmhp  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjmhp  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjmhp  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjmhp  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjmhp  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjmhp  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjmhp  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjmhp  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjmhp  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjmhp  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !26
    ajch3cl(ih)=                                                               &
      w1 (ih)*tabjch3cl(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjch3cl(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjch3cl(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjch3cl(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjch3cl(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjch3cl(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjch3cl(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjch3cl(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjch3cl(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjch3cl(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjch3cl(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjch3cl(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjch3cl(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjch3cl(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjch3cl(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjch3cl(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !27
    ajmcfm(ih)=                                                                &
      w1 (ih)*tabjmcfm (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjmcfm (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjmcfm (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjmcfm (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjmcfm (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjmcfm (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjmcfm (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjmcfm (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjmcfm (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjmcfm (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjmcfm (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjmcfm (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjmcfm (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjmcfm (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjmcfm (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjmcfm (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !28
    ajf11 (ih)=                                                                &
      w1 (ih)*tabjf11  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjf11  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjf11  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjf11  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjf11  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjf11  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjf11  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjf11  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjf11  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjf11  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjf11  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjf11  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjf11  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjf11  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjf11  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjf11  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !29
    ajf12 (ih)=                                                                &
      w1 (ih)*tabjf12  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjf12  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjf12  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjf12  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjf12  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjf12  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjf12  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjf12  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjf12  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjf12  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjf12  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjf12  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjf12  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjf12  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjf12  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjf12  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !30
    ajf22 (ih)=                                                                &
      w1 (ih)*tabjf22  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjf22  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjf22  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjf22  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjf22  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjf22  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjf22  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjf22  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjf22  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjf22  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjf22  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjf22  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjf22  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjf22  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjf22  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjf22  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !31
    ajf113(ih)=                                                                &
      w1 (ih)*tabjf113 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjf113 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjf113 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjf113 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjf113 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjf113 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjf113 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjf113 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjf113 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjf113 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjf113 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjf113 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjf113 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjf113 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjf113 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjf113 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !32
    ajccl4(ih)=                                                                &
      w1 (ih)*tabjccl4 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjccl4 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjccl4 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjccl4 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjccl4 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjccl4 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjccl4 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjccl4 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjccl4 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjccl4 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjccl4 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjccl4 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjccl4 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjccl4 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjccl4 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjccl4 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !33
    ajch3br(ih)=                                                               &
      w1 (ih)*tabjch3br(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjch3br(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjch3br(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjch3br(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjch3br(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjch3br(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjch3br(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjch3br(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjch3br(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjch3br(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjch3br(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjch3br(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjch3br(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjch3br(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjch3br(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjch3br(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !34
    ajf12b1(ih)=                                                               &
      w1 (ih)*tabjf12b1(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjf12b1(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjf12b1(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjf12b1(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjf12b1(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjf12b1(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjf12b1(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjf12b1(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjf12b1(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjf12b1(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjf12b1(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjf12b1(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjf12b1(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjf12b1(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjf12b1(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjf12b1(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !35
    ajf13b1(ih)=                                                               &
      w1 (ih)*tabjf13b1(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjf13b1(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjf13b1(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjf13b1(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjf13b1(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjf13b1(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjf13b1(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjf13b1(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjf13b1(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjf13b1(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjf13b1(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjf13b1(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjf13b1(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjf13b1(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjf13b1(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjf13b1(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !36
    ajcof2(ih)=                                                                &
      w1 (ih)*tabjcof2 (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjcof2 (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjcof2 (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjcof2 (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjcof2 (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjcof2 (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjcof2 (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjcof2 (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjcof2 (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjcof2 (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjcof2 (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjcof2 (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjcof2 (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjcof2 (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjcof2 (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjcof2 (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !37
    ajcofcl(ih)=                                                               &
      w1 (ih)*tabjcofcl(jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjcofcl(jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjcofcl(jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjcofcl(jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjcofcl(jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjcofcl(jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjcofcl(jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjcofcl(jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjcofcl(jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjcofcl(jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjcofcl(jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjcofcl(jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjcofcl(jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjcofcl(jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjcofcl(jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjcofcl(jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))
    !38
    ajch4 (ih)=                                                                &
      w1 (ih)*tabjch4  (jx  (ih),jy  (ih),jz  (ih),jo  (ih))                   &
      + w2 (ih)*tabjch4  (jxp1(ih),jy  (ih),jz  (ih),jo  (ih))                 &
      + w3 (ih)*tabjch4  (jxp1(ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w4 (ih)*tabjch4  (jx  (ih),jyp1(ih),jz  (ih),jo  (ih))                 &
      + w5 (ih)*tabjch4  (jx  (ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w6 (ih)*tabjch4  (jxp1(ih),jy  (ih),jzp1(ih),jo  (ih))                 &
      + w7 (ih)*tabjch4  (jxp1(ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w8 (ih)*tabjch4  (jx  (ih),jyp1(ih),jzp1(ih),jo  (ih))                 &
      + w9 (ih)*tabjch4  (jx  (ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w10(ih)*tabjch4  (jxp1(ih),jy  (ih),jz  (ih),jop1(ih))                 &
      + w11(ih)*tabjch4  (jxp1(ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w12(ih)*tabjch4  (jx  (ih),jyp1(ih),jz  (ih),jop1(ih))                 &
      + w13(ih)*tabjch4  (jx  (ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w14(ih)*tabjch4  (jxp1(ih),jy  (ih),jzp1(ih),jop1(ih))                 &
      + w15(ih)*tabjch4  (jxp1(ih),jyp1(ih),jzp1(ih),jop1(ih))                 &
      + w16(ih)*tabjch4  (jx  (ih),jyp1(ih),jzp1(ih),jop1(ih))

    !     Else darkness
  ELSE
    ajbrno3(ih) = 0.0
    ajhobr (ih) = 0.0
    ajoclo (ih) = 0.0
    ajc2oa (ih) = 0.0
    ajc2ob (ih) = 0.0
    ajmhp  (ih) = 0.0
    ajch3cl(ih) = 0.0
    ajmcfm (ih) = 0.0
    ajf11  (ih) = 0.0
    ajf12  (ih) = 0.0
    ajf22  (ih) = 0.0
    ajf113 (ih) = 0.0
    ajccl4 (ih) = 0.0
    ajch3br(ih) = 0.0
    ajf12b1(ih) = 0.0
    ajf13b1(ih) = 0.0
    ajcof2 (ih) = 0.0
    ajcofcl(ih) = 0.0
    ajch4  (ih) = 0.0
  END IF

END DO
IF (lhook) CALL dr_hook('CALCJS',zhook_out,zhook_handle)
RETURN

END SUBROUTINE calcjs
END MODULE calcjs_mod
