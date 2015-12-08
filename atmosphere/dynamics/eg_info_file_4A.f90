! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_info_file_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_info_file (mype,nproc_x,nproc_y,numcycles,innits,              &
                  global_row_length,global_rows,model_levels,                &
                  ih,height_domain,l_shallow,l_rotating,l_rotate_grid,       &
                  gcr_tol_abs,gcr_precon_option,                             &
                  l_inc_solver, l_RK_dps,                                    &
                  l_accel_convergence,l_eliminate_rho,                       &
                  l_init_fnm1,                                               &
                  alpha_u,alpha_v,alpha_w,alpha_theta,                       &
                  alpha_rho,alpha_p,l_rotate_winds,l_baroclinic,l_solid_body,&
                  l_baro_inst,l_baro_perturbed,                              &
                  l_deep_baro_inst, T0_E, T0_P, b_const, k_const,            &
                  L_const_grav,l_expl_horz_drag,                             &
                  l_impl_horz_drag,l_eg_dry_static_adj,l_fix_mass,           &
                  L_conserv_smooth_lap,                                      &
                  eg_vert_damp_coeff,eg_vert_damp_profile,eta_s,grid_np_lon, &
                  grid_np_lat,aa_jet_u0,l_cartesian,surface_type,            &
                  grid_number,tprofile_number,trefer_number,t_surface_in,    &
                  h_o,lambda_fraction,phi_fraction,half_width_x,half_width_y)

USE parkind1,     ONLY: jpim, jprb       !DrHook
USE yomhook,      ONLY: lhook, dr_hook   !DrHook
USE timestep_mod, ONLY: timestep
USE ereport_mod,  ONLY: ereport
USE eg_parameters_mod, ONLY : l_rho_av_zz
USE horiz_grid_mod, ONLY: Nxi1L, Nxi1V, Nxi2L, Nxi2V,                        &
                          delta_xi1_H, delta_xi1_L, delta_xi2_H, delta_xi2_L
USE domain_params
USE PrintStatus_mod
USE proc_info_mod, ONLY: model_domain

USE chsunits_mod, ONLY : nunits

IMPLICIT NONE
!
! Description: Writes summary of salient ENDGame variables in concise
!              format
!  
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER mype,nproc_x,nproc_y,numcycles,innits,global_row_length,      &
        global_rows,model_levels,                                     &
        gcr_precon_option,                                            &
        surface_type,tprofile_number,trefer_number,                   &
        aa_jet_m, aa_jet_n, grid_number,eg_vert_damp_profile,         &
        b_const, k_const


REAL    ih, height_domain, gcr_tol_abs,                               &
        alpha_u,alpha_v,alpha_w,alpha_theta,alpha_rho,                &
        alpha_p,grid_np_lon,grid_np_lat,aa_jet_u0,aa_jet_a,           &
        t_surface_in,h_o,lambda_fraction,phi_fraction,half_width_x,   &
        half_width_y,eg_vert_damp_coeff, eta_s, T0_E, T0_P 

LOGICAL l_shallow,l_rotating,l_rotate_grid,l_accel_convergence,       &
        l_eliminate_rho,l_init_fnm1,                                  &
        l_solid_body,l_baro_inst,                                     &
        l_baro_perturbed,l_baroclinic,l_cartesian,l_rotate_winds,     &
        l_const_grav,l_expl_horz_drag,l_impl_horz_drag,               &
        l_eg_dry_static_adj,l_fix_mass, l_deep_baro_inst,             &
        l_inc_solver, l_RK_dps,                                       &
        L_conserv_smooth_lap 

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


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_INFO_FILE',zhook_in,zhook_handle)
 
        IF( mype == 0 ) THEN

          IF( PrintStatus >= PrStatus_Normal ) THEN
            IF ( l_eliminate_rho ) THEN
              WRITE(6,fmt='(A)') ''
              WRITE(6,fmt='(A)') '*** USING DIAGNOSTIC DENSITY ***'
              WRITE(6,fmt='(A)') 'Mode not compatible with SLICE'
              WRITE(6,"(A,F12.6)") 'alpha_p = ',alpha_p
              WRITE(6,fmt='(A)') ''
            ELSE
              WRITE(6,fmt='(A)') ''
              WRITE(6,fmt='(A)') '*** USING PROGNOSTIC DENSITY ***'
              WRITE(6,fmt='(A)') ''
            END IF

            IF ( l_accel_convergence ) THEN

              CALL ereport('eg_info_file',1,                        &
                   'l_accel_convergence=T not available anymore')

            END IF
          END IF


         OPEN(UNIT=eg_unit,FILE='eg_job.info')
         WRITE(eg_unit,*) '**************************************************'
         WRITE(eg_unit,*) ' Summary of ENDGame parameters'
         WRITE(eg_unit,*)
         WRITE(eg_unit,*) ' Code version (local working copy may differ!)'
         WRITE(eg_unit,*) ' $Revision: 26957 $'
         WRITE(eg_unit,*) ' Note: Rev number based on last eg_atm_step change,'
         WRITE(eg_unit,*) '       not on last branch chage!'
         WRITE(eg_unit,*)
         WRITE(eg_unit,*)
         SELECT CASE(model_domain)
            CASE(mt_global)
               WRITE(eg_unit,*) 'Model Domain...................Global'
             CASE(mt_lam)
               WRITE(eg_unit,*) 'Model Domain...................LAM'   
             CASE(mt_cyclic_lam)
               WRITE(eg_unit,*) 'Model Domain...................cyclic E-W LAM'
             CASE(mt_bi_cyclic_lam)
               WRITE(eg_unit,*) 'Model Domain...................bicyclic LAM'
            CASE DEFAULT
               CALL ereport('eg_info_file',model_domain, 'Invalid model domain')
         END SELECT 
         WRITE(eg_unit,"(A,I3)") ' No. of E-W PEs................',nproc_x
         WRITE(eg_unit,"(A,I3)") ' No. of N-S PEs................',nproc_y
         WRITE(eg_unit,"(A,I5)") ' row_length....................',            &
                                   global_row_length
         WRITE(eg_unit,"(A,I5)") ' rows..........................',global_rows
         WRITE(eg_unit,"(A,I3)") ' model_levels..................',model_levels
         WRITE(eg_unit,"(A,F25.10)") ' height_domain.................',        &
                                   height_domain
         WRITE(eg_unit,"(A,F15.10)") ' timestep......................',timestep
         WRITE(eg_unit,*)
         WRITE(eg_unit,"(A,I3)") ' numcycles.....................',numcycles
         WRITE(eg_unit,"(A,I3)") ' innits........................',innits
         WRITE(eg_unit,*)
         WRITE(eg_unit,*) 'L_shallow.....................',l_shallow
         WRITE(eg_unit,"(A,F15.10)") ' Ih............................',ih
         WRITE(eg_unit,"(A,L1)") ' L_rotating....................',l_rotating
         WRITE(eg_unit,"(A,L1)") ' L_rotate_grid.................',l_rotate_grid
         WRITE(eg_unit,*)
         IF ( l_inc_solver ) WRITE(eg_unit,*) ' Using Incremental Solver '         
         WRITE(eg_unit,"(A,E20.10)") ' solver tol....................',        &
                                      gcr_tol_abs

         SELECT CASE(gcr_precon_option)
            CASE(0)
               WRITE(eg_unit,*) 'pre_type......................Diagonal'
            CASE(1)
               WRITE(eg_unit,*) 'pre_type......................ILU'
               WRITE(eg_unit,*) 'WARNING: This may not work '
               CALL ereport('eg_info_file',gcr_precon_option,                  &
                            'Preconditioner nolonger supported')
            CASE(2)
               WRITE(eg_unit,*) 'pre_type......................Jacobi'
            CASE(3)
               WRITE(eg_unit,*) 'pre_type......................ILU-Neumann'
               WRITE(eg_unit,*) 'WARNING: This may not work '
               CALL ereport('eg_info_file',gcr_precon_option,                  &
                            'Preconditioner nolonger supported')
            CASE(4)
               WRITE(eg_unit,*) 'pre_type......................SOR+tridiag'
            CASE DEFAULT
               CALL ereport('eg_info_file',gcr_precon_option,                  &
                            'Invalid preconditioner')
         END SELECT

         WRITE(eg_unit,*)
         WRITE(eg_unit,"(A,L1)") ' L_accel_convergence...........',            &
                                 l_accel_convergence
         WRITE(eg_unit,"(A,L1)") ' L_eliminate_rho...............',            &
                                 l_eliminate_rho
         WRITE(eg_unit,"(A,L1)") ' L_init_Fnm1...................',l_init_fnm1
         WRITE(eg_unit,*)

         WRITE(eg_unit,*)
         WRITE(eg_unit,"(A,F15.10)") ' ALPHA_U.......................',alpha_u
         WRITE(eg_unit,"(A,F15.10)") ' ALPHA_V.......................',alpha_v
         WRITE(eg_unit,"(A,F15.10)") ' ALPHA_W.......................',alpha_w
         WRITE(eg_unit,"(A,F15.10)") ' ALPHA_THETA...................',        &
                                    alpha_theta
         IF ( l_eliminate_rho ) THEN
           WRITE(eg_unit,"(A,F15.10)") ' ALPHA_P.......................',alpha_p
         ELSE
           WRITE(eg_unit,"(A,F15.10)") ' ALPHA_RHO.....................',      &
                                      alpha_rho
         END IF
         WRITE(eg_unit,*)
         WRITE(eg_unit,"(A,L1)") ' L_rotate_winds................',            &
                                l_rotate_winds
         WRITE(eg_unit,"(A,L1)") ' L_baroclinic..................',l_baroclinic
         WRITE(eg_unit,"(A,L1)") ' L_solid_body..................',l_solid_body
         WRITE(eg_unit,"(A,L1)") ' L_baro_inst...................',l_baro_inst
         WRITE(eg_unit,"(A,L1)") ' L_deep_baro_inst..............',            &
                                l_deep_baro_inst
         IF ( l_deep_baro_inst ) THEN
             WRITE(eg_unit,"(A,F15.10)") '   T0_Pole.....................',T0_P
             WRITE(eg_unit,"(A,F15.10)") '   T0_Equator..................',T0_E   
             WRITE(eg_unit,"(A,I2)")   '   k..........................',k_const
             WRITE(eg_unit,"(A,I2)")   '   b..........................',b_const 
         END IF
         WRITE(eg_unit,"(A,L1)") ' L_baro_perturbed..............',            &
                                l_baro_perturbed
         WRITE(eg_unit,*)
         WRITE(eg_unit,"(A,F15.10)") ' grid_NP_lon...................',        &
                                       grid_np_lon
         WRITE(eg_unit,"(A,F15.10)") ' grid_NP_lat...................',        &
                                       grid_np_lat
         WRITE(eg_unit,"(A,F15.10)") ' AA_jet_u0.....................',aa_jet_u0
         WRITE(eg_unit,"(A,L1)") ' L_Cartesian...................',l_cartesian
         SELECT CASE(surface_type)
           CASE(0)
               WRITE(eg_unit,*) 'surface_type..............surface_flat'
!            CASE(1)
!              WRITE(eg_unit,*) 'surface_type..............surface_ellipse'
!            CASE(2)
!              WRITE(eg_unit,*) 'surface_type..............surface_ridge'
!            CASE(3)
!              WRITE(eg_unit,*) 'surface_type..............surface_plateau'
!            CASE(4)
!              WRITE(eg_unit,*) 'surface_type..............surface_massif'
!            CASE(5)
!              WRITE(eg_unit,*) 'surface_type..............surface_mask'
           CASE(6)
               WRITE(eg_unit,*) 'surface_type..............surface_gauss'
!            CASE(7)
!              WRITE(eg_unit,*) 'surface_type..............surface_ridge_series'
           CASE(8)
               WRITE(eg_unit,*) 'surface_type..............surface_schar_ridge'
           CASE(9)
               WRITE(eg_unit,*) 'surface_type..............surface_baroclinic'
           CASE(10)
               WRITE(eg_unit,*) 'surface_type..............surface_dump'
           CASE DEFAULT
             CALL ereport('eg_info_file',1, 'Invalid surface type')
         END SELECT
         IF(surface_type /= 0 .AND. surface_type /= 10) THEN
           WRITE(eg_unit,"(A,F15.10)") '   h_o.....................',h_o
           WRITE(eg_unit,"(A,F15.10)") '   lambda_fraction.........',       &
                                           lambda_fraction
           WRITE(eg_unit,"(A,F15.10)") '   phi_fraction............',       &
                                           phi_fraction
           WRITE(eg_unit,"(A,F25.10)") '   half_width_x............',       &
                                           half_width_x
           WRITE(eg_unit,"(A,F25.10)") '   half_width_y............',       &
                                           half_width_y
         ENDIF
         WRITE(eg_unit,"(A,I3)") ' grid_number...............',grid_number
         WRITE(eg_unit,"(A,I3)") ' tprofile_number...........',tprofile_number
         WRITE(eg_unit,"(A,I3)") ' trefer_number.............',trefer_number
         WRITE(eg_unit,"(A,F15.10)") ' t_surface_in..............',          &
                                       t_surface_in

         WRITE(eg_unit,*)
         WRITE(eg_unit,"(A,L1)") ' L_rk_dps................. ', l_RK_dps

         WRITE(eg_unit,*)
         IF (l_const_grav) THEN
           WRITE(eg_unit,"(A)") 'CONSTANT GRAVITY SELECTED (l_const_grav)!  ' 
           WRITE(eg_unit,"(A)") 'Whilst this option is required for starting' 
           WRITE(eg_unit,"(A)") 'from ND dumps it is entirely at the users risk'
           WRITE(eg_unit,"(A)") 'It is not supported and known to lead to'
           WRITE(eg_unit,"(A)") 'a spurious 3D divergence tendency'
         ELSE
          WRITE(eg_unit,"(A)") 'gravity is height dependent (note: this is'
          WRITE(eg_unit,"(A)") 'incompatible with running from ND start dumps!)'
         END IF

         WRITE(eg_unit,*)
         WRITE(eg_unit,"(A,L1)") ' l_expl_horz_drag ............. ',           &
                                   l_expl_horz_drag
         WRITE(eg_unit,"(A,L1)") ' l_impl_horz_drag ............. ',           &
                                   l_impl_horz_drag
         WRITE(eg_unit,"(A,L1)") ' l_eg_dry_static_adj .......... ',           &
                                   l_eg_dry_static_adj
         WRITE(eg_unit,"(A,L1)") ' l_fix_mass ................... ',           &
                                   l_fix_mass
         WRITE(eg_unit,"(A,I1)") ' eg_vert_damp_profile ......... ',           &
                                   eg_vert_damp_profile
         WRITE(eg_unit,"(A,F15.10)") ' eta_s ........... ',                    &
                                       eta_s
         WRITE(eg_unit,"(A,E15.10)") ' eg_vert_damp_coeff ........... ',       &
                                       eg_vert_damp_coeff
         WRITE(eg_unit,"(A,I3)") '******************************************'
        CLOSE(eg_unit)
       END IF

IF (lhook) CALL dr_hook('EG_INFO_FILE',zhook_out,zhook_handle)

END SUBROUTINE eg_info_file
END MODULE eg_info_file_mod
