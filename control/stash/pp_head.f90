! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PPHEAD------------------------------------------
!LL
!LL  Creates a 64 word PP header from the the following:-
!LL  1)  PP_XREF (PP cross-reference array record for this sect/item)
!LL  2)  FIXED length header
!LL  3)  INTEGER constants array
!LL  4)  REAL constants array
!LL  5)  Some input arguments
!LL
!LL  Programming standard: U M DOC  Paper 3 vn8.3
!LL
!LLEND-------------------------------------------------------------

!
!*L  INTERFACE and ARGUMENTS:------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

      SUBROUTINE PP_HEAD(                                               &
          im_ident,fixhd,inthd,realhd,                                  &
          len_inthd,len_realhd,ie,is,gr,                                &
          lfullfield,real_level,pseudo_level,                           &
          samples,start,start_or_verif_time,end_or_data_time,pp_len,    &
          extraw,pp_int_head,pp_real_head,n_cols_out,num_words,         &
          len_buf_words,n_rows_out,nrow_in,srow_in,wcol_in,ecol_in,     &
          lbproc_comp,                                                  &
          sample_prd,fcst_prd,comp_accrcy,packing_type,                 &
          st_grid,iwa,zseak_rho,ck_rho,zseak_theta,ck_theta,            &
          model_levels,levindex,rotate,elf,                             &
          icode,cmessage)

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE model_file, ONLY : &
          mf_data_missing
      USE yomhook, ONLY:  &
          lhook,          &
          dr_hook
      USE parkind1, ONLY: &
          jprb,           &
          jpim
      USE IAU_mod, ONLY : &
          L_IAU,          &
          IAU_StartMin,   &
          IAU_EndMin
      USE UM_ParVars
      USE Control_Max_Sizes

      USE c_model_id_mod, ONLY: model_id
      USE lookup_addresses
      USE Submodel_Mod
      USE model_id_mod, ONLY: itab
      USE missing_data_mod, ONLY: imdi, rmdi

      USE cppxref_mod, ONLY:                                            &
          ppx_field_code,ppx_lbvc_code,ppx_lv_code,                     &
          ppx_meto8_fieldcode,ppx_meto8_levelcode,ppx_meto8_surf,       &
          ppx_data_type,ppx_full_level,ppx_half_level
          

      USE nlstcall_mod, ONLY : model_analysis_mins, &
                               expt_id, &
                               job_id, &
                               model_status

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

      INTEGER, PARAMETER :: LEN_FIXHD = 256

      CHARACTER(LEN=80) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE
      CHARACTER(LEN=5)   :: run_id ! To calculate LBEXP in header

      INTEGER, INTENT(IN) ::      &
          start_or_verif_time(7), &! verif time/start time for means etc
          end_or_data_time(7),    &! data time/end time for means etc
          samples,                &! no of samples in period (timeseries)
          im_ident,               &! internal model identifier
          pp_len,                 &! length of the lookup table
          len_inthd,              &! length of the integer constants
          len_realhd,             &! length of the real constants
          fixhd(len_fixhd),       &! array of fixed constants
          inthd(len_inthd),       &! array of integer constants
          st_grid,                &! stash horizontal grid type
          model_levels,           &! no of model levels
          levindex,               &! level index
          n_rows_out,             &! pphoriz_out=n_rows_out*n_cols_out+extra
          n_cols_out,             &! pphoriz_out=n_cols_out*n_rows_out+extra
          nrow_in,srow_in,        &! the most nrthrly/southerly row.
          wcol_in,ecol_in,        &! the most westerly/easterly column
          pseudo_level,           &! output pp pseudo-level
          comp_accrcy,            &! packing accuracy in power of 2
          packing_type,           &! 0 = no packing, 1 = wgdos, 3 = grib
          num_words,              &! number of 64 bit words to hold data
          extraw,                 &! number of extra-data words
          len_buf_words,          &! number of 64 bit words (rounded to 512)
          iwa,                    &! start word address.
          ie,                     &! item number
          is,                     &! section number
          gr,                     &! grid point code
          lbproc_comp(14)          ! subcomponents(0/1) to make up lbproc

      LOGICAL, INTENT(IN) ::      &
          start,                  &! flag to control update for verif/start time
          lfullfield               ! TRUE if output field on full horiz domain

      REAL, INTENT(IN) ::              &
          fcst_prd,                    &! forecast period
          realhd(len_realhd),          &! real header
          real_level,                  &! output pp level(real)
          sample_prd,                  &! sampling period in hours for time mean
          zseak_rho    (model_levels), &! vert coeff zsea on rho levels
          ck_rho       (model_levels), &! vert coeff ck   on rho levels
          zseak_theta(0:model_levels), &! vert coeff zsea on theta levels
          ck_theta   (0:model_levels)   ! vert coeff ck   on theta levels
          
      INTEGER, INTENT(OUT) ::          &
          icode,                       &! return code from the routine
          pp_int_head(pp_len)           ! integer lookup table
      REAL, INTENT(OUT) ::             &
          PP_REAL_HEAD(PP_LEN)          ! Real Lookup table
!
!*---------------------------------------------------------------------
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end
! Contains *CALL VERSION
! --------------------- Comdeck: CHISTORY ----------------------------
!
!  Purpose: COMMON block for history data needed by top level (C0)
!           routines, and passed from run to run.  Mostly set by
!           the User Interface.
!
!           Note that CHISTORY *CALLs ALL individual history comdecks
!
! --------------------------------------------------------------------
!
! ----------------------- Comdeck: IHISTO   ----------------------------
! Description: COMDECK defining Integer History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations


      ! Array containing model data time (Same as MODEL_BASIS_TIME/MODEL
      ! ANALYSIS_HRS depending whether before/after assimilation)
      INTEGER :: model_data_time(6)

      ! Indicator of operational run type
      INTEGER :: run_indic_op

      ! Final target date for the run
      INTEGER :: run_resubmit_target(6)

      ! Last field written/read per FT unit
      INTEGER :: ft_lastfield(20:nunits)

      ! Number of automatically-resubmitted job chunks
      ! Used to name output file
      INTEGER :: run_job_counter

! History Common Block for overall model integers variables.

      COMMON /ihisto/                                                 &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

      NAMELIST /nlihisto/                                             &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

! IHISTO end
! ----------------------- Comdeck: CHISTO   ----------------------------
! Description: COMDECK defining Character History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

  CHARACTER(LEN=8) ::  run_type             ! Type of run
  CHARACTER(LEN=1) ::  ft_active(20:nunits) ! "Y" if file partly written

  LOGICAL :: newrun ! Set to true in NRUN to stop auto-resubmission

  ! History Common Block for overall model character variables.

  COMMON /chisto/                                     &
     run_type,                                        &
     newrun, ft_active

  NAMELIST /nlchisto/                                 &
     run_type,                                        &
     ft_active

! CHISTO end
! ----------------------- Comdeck: IHISTG   ----------------------------
! Description: COMDECK defining Integer History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! This file belongs in section: Top Level
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations
      ! History block copy of A_STEP held in file CTIME
      INTEGER :: h_stepim(n_internal_model_max)

      ! No of means activated
      INTEGER :: mean_offsetim(n_internal_model_max)

      ! Offset between MEAN_REFTIME and model basis time(in model dumps)
      INTEGER :: offset_dumpsim(n_internal_model_max)

      ! No of mean periods chosen
      INTEGER :: mean_numberim(n_internal_model_max)

      ! Indicators used to correct logical units are used for
      ! atmos partial sum dump I/O
      INTEGER :: run_meanctl_indicim(4,n_internal_model_max)

      ! History Common Block for generic model integer variables.

      COMMON /ihistg/                                         &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         mean_numberim, run_meanctl_indicim

      NAMELIST /nlihistg/                                     &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         run_meanctl_indicim

! IHISTG end
! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character variables for
!              managing dump names
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 
!
!   Type declarations
!
! For keeping old restart dump name between creation of new restart dump
! and successful completion of climate means and history file.
CHARACTER(LEN=256) :: save_dumpname_im(n_internal_model_max)
! Name of current restart dump
CHARACTER(LEN=256) :: checkpoint_dump_im(n_internal_model_max)
! Blank name
CHARACTER(LEN=256) :: blank_file_name
!
! History Common Block for generic model characters variables.
!
COMMON /chistg/save_dumpname_im, checkpoint_dump_im, blank_file_name

NAMELIST /nlchistg/checkpoint_dump_im

! CHISTG end
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

!
! Description: This include file contains information needed when
!              generating variable horizontal grid data in the
!              STASH extra data vector. Introduced UM 5.4 - R. Hill
!===================================================================
      LOGICAL :: X_VAR_GRID ! Whether variable grid in E-W direction
      LOGICAL :: Y_VAR_GRID ! and/or in S-N direction

      INTEGER :: VAR_GRID_TYPE ! 0 = none
                               ! 1 = T grid
                               ! 2 = U/V grid

      ! Grid boundaries for T and U,V
      REAL :: X_BOUNDARY(ROW_LENGTH_MAX+1,2)
      REAL :: Y_BOUNDARY(ROWS_MAX+1,2)

      ! Grid Points for T and U,V
      REAL :: X_GRID(ROW_LENGTH_MAX,2)
      REAL :: Y_GRID(ROWS_MAX,2)

      COMMON /OVARGRID/ X_VAR_GRID,Y_VAR_GRID                           &
     & ,X_BOUNDARY,Y_BOUNDARY,X_GRID,Y_GRID,VAR_GRID_TYPE

      ! The following parameters correspond to the extra data
      ! vector descriptors expected, for e.g., in PV-WAVE
      ! plotting routines (e.g. decode_extra.pro). There are
      ! numerous other areas of code where these integer
      ! descriptors must be handled (e.g. FIELDCOS, PPI2H, FTT)
      ! So it is not a trivial matter to introduce new code
      ! descriptors. Furthermore, ieee -32 will destroy these
      ! integers so PP data must always be processed via the
      ! long winded route: QXFIELDCOS -> PUTIBM -> FTT/PPI2H.
      ! (This thoroughly unsatisfactory state of affairs may
      ! be correctable with developments to ieee and convpp).
      INTEGER,PARAMETER :: x_coord_vector=1
                                     ! Indicates that an extra
                                     ! data vector gives LBNPT
                                     ! x-coordinate values
      INTEGER,PARAMETER :: y_coord_vector=2
                                     ! Indicates that an extra
                                     ! data vector gives LBROW
                                     ! Y-coordinate values
      INTEGER,PARAMETER :: x_lbnd_vector=12
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! x-boundary values
      INTEGER,PARAMETER :: x_ubnd_vector=13
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! x-boundary values
      INTEGER,PARAMETER :: y_lbnd_vector=14
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! y-boundary values
      INTEGER,PARAMETER :: y_ubnd_vector=15
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! y-boundary values
! CTIME ----------------------------------------------------
!
!  Purpose: Derived model time/step information including start/end
!           step numbers and frequencies (in steps) of interface field
!           generation, boundary field updating, ancillary field
!           updating; and assimilation start/end times.
!           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!           Also contains current time/date information, current
!           step number (echoed in history file) and steps-per-group.
!
!END -----------------------------------------------------------------

      INTEGER :: I_YEAR               ! Current model time (years)
      INTEGER :: I_MONTH              ! Current model time (months)
      INTEGER :: I_DAY                ! Current model time (days)
      INTEGER :: I_HOUR               ! Current model time (hours)
      INTEGER :: I_MINUTE             ! Current model time (minutes)
      INTEGER :: I_SECOND             ! Current model time (seconds)
      INTEGER :: I_DAY_NUMBER         ! Current model time (day no)
      INTEGER :: PREVIOUS_TIME(7)     ! Model time at previous step
      INTEGER :: IAU_DTResetStep      ! Data time reset step for IAU run

      INTEGER :: BASIS_TIME_DAYS  ! Integral no of days to basis time
      INTEGER :: BASIS_TIME_SECS  ! No of seconds-in-day at basis time

      LOGICAL :: L_C360DY

! UM6.5MODEL_ANALYSIS_HRS changed to REAL - 
!   requires FORECAST_HRS and DATA_MINUS_BASIS_HRS to REAL also 
      REAL    :: FORECAST_HRS     ! Hours since Data Time (ie T+nn)
      REAL    :: DATA_MINUS_BASIS_HRS ! Data time - basis time (hours)

      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
        I_DAY_NUMBER,PREVIOUS_TIME,                                     &
        BASIS_TIME_DAYS,BASIS_TIME_SECS,                                &
        FORECAST_HRS,DATA_MINUS_BASIS_HRS,                              &
        IAU_DTResetStep, L_C360DY

      INTEGER :: STEPim(INTERNAL_ID_MAX)  ! Step no since basis time
      INTEGER :: GROUPim(INTERNAL_ID_MAX) ! Number of steps per group

      ! Finish step number this run
      INTEGER :: TARGET_END_STEPim(INTERNAL_ID_MAX)

      REAL :: SECS_PER_STEPim(INTERNAL_ID_MAX) ! Timestep length in secs

      ! Frequency of interface field generation in steps
      INTEGER :: INTERFACE_STEPSim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Start steps for interface field generation
      INTEGER :: INTERFACE_FSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! End steps for interface field generation
      INTEGER :: INTERFACE_LSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Frequency of  updating boundary fields in steps
      INTEGER :: BOUNDARY_STEPSim(INTERNAL_ID_MAX)

      ! No of steps from boundary data prior to basis time to model
      ! basis time
      INTEGER :: BNDARY_OFFSETim(INTERNAL_ID_MAX)

      ! Lowest frequency for updating of ancillary fields in steps
      INTEGER :: ANCILLARY_STEPSim(INTERNAL_ID_MAX)

      ! Start steps for assimilation
      INTEGER :: ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX)

      ! Number of assimilation steps to analysis
      INTEGER :: ASSIM_STEPSim(INTERNAL_ID_MAX)

      ! Number of assimilation steps after analysis
      INTEGER :: ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)

      COMMON/CTIMEE/                                                    &
     &  STEPim,GROUPim,TARGET_END_STEPim,INTERFACE_STEPSim,             &
     &  INTERFACE_FSTEPim,INTERFACE_LSTEPim,BOUNDARY_STEPSim,           &
     &  BNDARY_OFFSETim,ANCILLARY_STEPSim,ASSIM_FIRSTSTEPim,            &
     &  ASSIM_STEPSim,ASSIM_EXTRASTEPSim,SECS_PER_STEPim

! CTIME end

      EXTERNAL EXPPXI

!    DEFINE LOCAL VARIABLES

      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='PP_head')
      INTEGER                         &
          pp_lbfc,                    &!  m08 level code
          pp_lbtyp,                   &!  m08 field type code
          pp_lblev,                   &!  m08 field level code
          pp_iproj,                   &!  m08 projection number
          pp_lbvc,                    &!  vertical coord type
          ii,                         &!  local counter
          int_level,                  &!  integer value of level
          k,                          &!  local counter
          ia,ib,ic,                   &!  component codes to make up lbtim
          mean_code,                  &!  spatial averaging code derived from gr
          lvcode,                     &!  lv code
          exppxi,                     &!  function to extract ppxref info
          exptcode                     !  integer coded experiment name
      
      INTEGER :: get_um_version_id
      
      ! Local Parameters
      INTEGER, PARAMETER :: unset_exptenc = -5555
      
      ! Local scalars:
      INTEGER,SAVE  ::  expt_enc_answer=unset_exptenc

      LOGICAL                          &
          elf,                         &
          rotate

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!LL   Construct PP header
!
!  Timestamps ----------------------------------------------------------
!
!L
!L Set up time info dependent on start flag.
!L For all but time series start will be TRUE so all time information
!L will be set up from FIXHD in effect, but for time series start
!L will be set up by TEMPORAL and passed in, so that dump headers are
!L set correctly for such fields.
!L Note: end_or_data_time will be updated from current model time in
!L       FIXHD(28-34) for time means/accumulations etc.
!L
      IF (lhook) CALL dr_hook('PP_HEAD',zhook_in,zhook_handle)
      IF (start) THEN    ! start timestep so update start time
        PP_INT_HEAD(LBYR)=start_or_verif_time(1)
        PP_INT_HEAD(LBMON)=start_or_verif_time(2)
        PP_INT_HEAD(LBDAT)=start_or_verif_time(3)
        PP_INT_HEAD(LBHR)=start_or_verif_time(4)
        PP_INT_HEAD(LBMIN)=start_or_verif_time(5)
        PP_INT_HEAD(LBSEC)=start_or_verif_time(6)
!        PP_INT_HEAD(LBDAY)=start_or_verif_time(7)
! If the IAU scheme was used on the previous timestep to add a complete
! increment at the nominal analysis time, reset verification minute to
! zero for fields not in sections 0, 15 or 16. Required so we can obtain
! "analysis" fields from the physics modules.
! UM6.5 - model_analysis_hrs replaced by model_analysis_mins
        IF (L_IAU .AND. im_ident == A_IM) THEN
          IF ( IAU_StartMin == IAU_EndMin              .AND.            &
               IAU_StartMin == model_analysis_mins     .AND.            &
               STEPIM(A_IM) == IAU_DTResetStep + 1     .AND.            &
               IS /= 0                                 .AND.            &
               IS /= 15                                .AND.            &
               IS /= 16 ) THEN
            PP_INT_HEAD(LBMIN) = 0
          END IF
        END IF
      END IF
      PP_INT_HEAD(LBYRD)=end_or_data_time(1)
      PP_INT_HEAD(LBMOND)=end_or_data_time(2)
      PP_INT_HEAD(LBDATD)=end_or_data_time(3)
      PP_INT_HEAD(LBHRD)=end_or_data_time(4)
      PP_INT_HEAD(LBMIND)=end_or_data_time(5)
      PP_INT_HEAD(LBSECD)=end_or_data_time(6)
!      PP_INT_HEAD(LBDAYD)=end_or_data_time(7)
!
!  Secondary time information ------------------------------------------
!
! LBTIM is 100*IA+10*IB+IC - this encodes the time processing type
!
      IA=INT(sample_prd)           ! Sampling period in whole hours
      IF(sample_prd == 0.0) THEN   ! NB: may be a fraction of an hour
        IB=1                       ! Forecast field
      ELSE
        IF (IA == 0) THEN
          IA=1                     ! 0 < sample_prd < 1 counts as 1 hour
        ENDIF
        IB=2                       ! Time mean or accumulation
      ENDIF
      IC=FIXHD(8)                  ! Calendar (1: Gregorian, 2: 360 day)
!
      PP_INT_HEAD(LBTIM)=100*IA+10*IB+IC
      PP_INT_HEAD(LBFT)=FCST_PRD
!
!  Data length ---------------------------------------------------------
!
      PP_INT_HEAD(LBLREC)=NUM_WORDS
!
!  Grid code (determined from dump fixed-length header) ----------------
!
      IF (samples == 0) THEN
!       Field is not a timeseries
        IF(FIXHD(4) <  100) THEN
          IF (VAR_GRID_TYPE == 0) THEN
             PP_INT_HEAD(LBCODE)=1   ! Regular lat/long grid
          ELSE
             ! Set variable  grid code to the same as regular
             ! grid which effectively renders this peice of
             ! information meaningless
             PP_INT_HEAD(LBCODE) = 1   ! Variable grid
          ENDIF
        ELSE
          IF (VAR_GRID_TYPE == 0) THEN
             PP_INT_HEAD(LBCODE)=101   ! Regular lat/long grid
          ELSE
             ! Set variable  grid code to the same as regular
             ! grid which effectively renders this peice of
             ! information meaningless
             PP_INT_HEAD(LBCODE) = 101 ! Variable grid
          ENDIF

        ENDIF
      ELSE
!       Field is a timeseries
        PP_INT_HEAD(LBCODE)=31300
        IF (FIXHD(8) == 1) THEN
!         Calendar --  1: Gregorian
          PP_INT_HEAD(LBCODE)=PP_INT_HEAD(LBCODE)+20
        ELSEIF (FIXHD(8) == 2) THEN
!         Calendar -- 360 day (Model Calendar)
          PP_INT_HEAD(LBCODE)=PP_INT_HEAD(LBCODE)+23
        ELSE
!         Unknown calendar. Fail.
          ICODE=2
      CMESSAGE='PPHEAD: unknown calender type in fixhd(8)'
        ENDIF
      ENDIF
!
!  Hemispheric subregion indicator -------------------------------------
!
      IF (samples >  0 .OR. .NOT.lfullfield) THEN
!  Field is a timeseries/trajectory or subdomain of the full model area
        PP_INT_HEAD(LBHEM)=3
      ELSEIF (FIXHD(4) <  100) THEN
!  Otherwise, use the value for the full model area encoded in the dump
        PP_INT_HEAD(LBHEM)=FIXHD(4)
      ELSE
        PP_INT_HEAD(LBHEM)=FIXHD(4)-100
      ENDIF
!
!  Field dimensions (rows x cols) --------------------------------------
!
      PP_INT_HEAD(LBROW)=N_ROWS_OUT
      PP_INT_HEAD(LBNPT)=N_COLS_OUT
!
!  'Extra data' length (now accomodates timeseries sampling data) ------
!
      PP_INT_HEAD(LBEXT)=extraw
!
!  Packing method indicator (new definition introduced at vn2.8)--------
       IF(PACKING_TYPE == 1)THEN    ! WGDOS packing
         PP_INT_HEAD(LBPACK)=00001
       ELSEIF(PACKING_TYPE == 4)THEN ! Run length encoding
         PP_INT_HEAD(LBPACK)=00004
       ELSEIF(PACKING_TYPE == 3)THEN ! GRIB packing
         PP_INT_HEAD(LBPACK)=00003
       ELSEIF(PACKING_TYPE == 0)THEN ! No packing
         PP_INT_HEAD(LBPACK)=00000
       ELSEIF(PACKING_TYPE == mf_data_missing)THEN ! No idea
         PP_INT_HEAD(LBPACK)=mf_data_missing ! Not set code
       ELSE
         ICODE=1
         CMESSAGE='PPHEAD  Packing type undefined'
         PP_INT_HEAD(LBPACK)=00000
      ENDIF
!
!  PP header release no ------------------------------------------------
!
      PP_INT_HEAD(LBREL)=3   ! for seconds in LBSEC(D), not day_number
!      PP_INT_HEAD(LBREL)=2
!
!  Primary fieldcode (some hardwiring for ELF winds) -------------------
!  Secondary fieldcode not used currently
!
! DEPENDS ON: exppxi
      PP_LBFC=EXPPXI(im_ident, is, ie, ppx_field_code,                  &
                     icode, cmessage)
      IF(ELF.AND..NOT.ROTATE) THEN  ! ELF winds are in x,y direction
        IF(PP_LBFC == 56) PP_LBFC=48
        IF(PP_LBFC == 57) PP_LBFC=49
      ENDIF
      PP_INT_HEAD(LBFC)=PP_LBFC
      PP_INT_HEAD(LBCFC)=0
!
!  Processing code (encodes several things in one field) ---------------
!
      PP_INT_HEAD(LBPROC)=0
      DO II=14,1,-1
        PP_INT_HEAD(LBPROC)=PP_INT_HEAD(LBPROC)*2+LBPROC_COMP(II)
      ENDDO
!
!  Vertical coordinate type --------------------------------------------
!  Vertical coordinate type for reference level not coded
!
! DEPENDS ON: exppxi
      PP_LBVC=EXPPXI(im_ident, is, ie, ppx_lbvc_code,                   &
                     icode, cmessage)
      PP_INT_HEAD(LBVC)=PP_LBVC
      PP_INT_HEAD(LBRVC)=0

! [Note that, although most diagnostics are defined over int_level=
! (1:model_levels), and hence int_level=LevIndex, some variables are
! defined over int_level=(0:model_levels). For this special case,
! LevIndex is (1:model_levels+1), since LevIndex is defined starting at
! 1, and int_level != LevIndex.]
      int_level = real_level+0.00001  ! ensure no rounding problems
!
!  Experiment number coded from expt_id and job_id for non
!  operational set to RUN_INDIC_OP for operational use.
!
      IF (model_status /= 'Operational') THEN
        RUN_ID(1:4)=expt_id
        RUN_ID(5:5)=job_id
        IF (expt_enc_answer == unset_exptenc) THEN
!  Function EXPT_ENC will encode the run_id into a unique integer
! DEPENDS ON: expt_enc
           CALL EXPT_ENC(RUN_ID,EXPTCODE,ICODE,CMESSAGE)
           expt_enc_answer=EXPTCODE

        ELSE
          EXPTCODE=expt_enc_answer
        ENDIF


        PP_INT_HEAD(LBEXP)=EXPTCODE          ! LBEXP
      ELSE
        PP_INT_HEAD(LBEXP)=RUN_INDIC_OP      ! LBEXP (ITAB)
      ENDIF
      
!  new code specifically for Rose. Once migrated to Rose this
!  becomes the default for all model runs and above exptid runid is retired.
      
      IF (itab /= imdi) PP_INT_HEAD(LBEXP)= itab 
!
!  Direct access dataset start address and no of records ---------------
!
      PP_INT_HEAD(LBEGIN)=IWA
      PP_INT_HEAD(LBNREC)=LEN_BUF_WORDS
!
!  Operational fieldsfile projection no, fieldtype + level codes -------
!  These are hardwired according to model resolution
!
      IF(INTHD(6) == 192) THEN
        PP_IPROJ=802
      ELSE IF(INTHD(6) == 288) THEN
        PP_IPROJ=800
      ELSE IF(INTHD(6) == 96) THEN
        PP_IPROJ=870
       ELSE IF(INTHD(6) == 432) THEN
         PP_IPROJ=800
      ELSE
        PP_IPROJ=900
      ENDIF
! DEPENDS ON: exppxi
      PP_LBTYP=EXPPXI(im_ident, is, ie, ppx_meto8_fieldcode,            &
                     icode, cmessage)
! DEPENDS ON: exppxi
      lvcode=EXPPXI(im_ident, is, ie, ppx_lv_code,                      &
                   icode, cmessage)
      IF(real_level == -1.0) THEN
! DEPENDS ON: exppxi
        PP_LBLEV=EXPPXI(im_ident, is, ie, ppx_meto8_levelcode,          &
                     icode, cmessage)   ! levelcode 9999 or 8888
      ELSE
        IF (im_ident  ==  atmos_im) THEN
          IF (lvcode == ppx_half_level .AND. int_level == 0) THEN
!             This is a surface level: reset lblev
            PP_LBLEV=ppx_meto8_surf
          ELSE
            PP_LBLEV=int_level
          ENDIF
        ELSE
          PP_LBLEV=int_level
        ENDIF
      ENDIF
      PP_INT_HEAD(LBPROJ)=PP_IPROJ
      PP_INT_HEAD(LBTYP)=PP_LBTYP
      PP_INT_HEAD(LBLEV)=PP_LBLEV
!
!  Reserved slots for future expansion ---------------------------------
!
      PP_INT_HEAD(LBRSVD1)=0
      PP_INT_HEAD(LBRSVD2)=0
      PP_INT_HEAD(LBRSVD3)=0
      PP_INT_HEAD(LBRSVD4)=0 ! Used externally for Ensemble member number
!
! Generate model version_id
! DEPENDS ON: get_um_version_id
      PP_INT_HEAD(LBSRCE)=get_um_version_id(model_id)
!
! Data type - extract from PPXREF
! DEPENDS ON: exppxi
      PP_INT_HEAD(DATA_TYPE)=EXPPXI(im_ident, is, ie, ppx_data_type,    &
                                    icode, cmessage)
!
!  Address within dump or PP file --------------------------------------
!
      IF (IWA==mf_data_missing) THEN
        ! If we don't know the file location
        ! we also cannot know the data offset
        PP_INT_HEAD(NADDR)=mf_data_missing
      ELSE
        ! The data addr is the addr in the file 
        ! from the data start location, except for fieldsfiles
        ! where is is stated to be the same as LBEGIN.
        PP_INT_HEAD(NADDR)=IWA
      END IF
!
!  LBUSER3 is not currently used (ie set to 0).
!
      PP_INT_HEAD(LBUSER3)=0
!
!  STASH section/item code ---------------------------------------------
!
      PP_INT_HEAD(ITEM_CODE)=IS*1000+IE
!
!  STASH pseudo-level (for fields which have pseudo-levels defined) ----
!
      PP_INT_HEAD(LBPLEV)=pseudo_level
!
!  Spare for user's use ------------------------------------------------
!
      PP_INT_HEAD(LBUSER6)=0
      PP_INT_HEAD(MODEL_CODE) = im_ident
!
!  Reserved for future PP package use ----------------------------------
!
      PP_REAL_HEAD(BRSVD3)=0.0
      PP_REAL_HEAD(BRSVD4)=0.0
      PP_REAL_HEAD(BDATUM)=0.0
      PP_REAL_HEAD(BACC)=COMP_ACCRCY ! packing accuracy stored as real
!
!  Vertical grid description -------------------------------------------
!  Level and reference level
!
      IF(PP_LBVC >= 126.AND.PP_LBVC <= 139) THEN ! Special codes
!                                                  (surf botttom,
!                                                   top all zero)
        PP_REAL_HEAD(BLEV)=0.0
        PP_REAL_HEAD(BHLEV)=0.0
        PP_REAL_HEAD(BRLEV)=0.0
        PP_REAL_HEAD(BHRLEV)=0.0
        PP_REAL_HEAD(BULEV)=0.0
        PP_REAL_HEAD(BHULEV)=0.0
      ELSEIF(PP_LBVC == 9.OR.PP_LBVC == 65) THEN ! Hybrid/ETA levels
! DEPENDS ON: exppxi
        lvcode=EXPPXI(im_ident, is, ie, ppx_lv_code,                    &
          icode, cmessage)

! From vn5.2:
! height of model level k above mean sea level is
!       z(i,j,k) = zsea(k) + C(k)*zorog(i,j)
! bulev,bhulev      zsea,C of upper layer boundary
! blev ,bhlev       zsea,C of level
! brlev,bhrlev      zsea,C of lower level boundary
! The level here can refer to either a theta or rho level, with
! layer boundaries defined by surrounding rho or theta levels.
!
! [Note assumption that the top of the atmosphere model,
!  ie eta_theta_levels(model_levels) = 1.0, is equidistant from a
!  bounding p level above (not represented explicitly) and the p
!  level below (at eta_rho_levels(model_levels) ).]

        IF (lvcode == ppx_half_level) THEN ! theta level (& w)

          IF(int_level >= model_levels) THEN ! top level
            PP_REAL_HEAD(bulev) =  zseak_theta(model_levels) * 2.0      &
                                         - zseak_rho(model_levels)
            PP_REAL_HEAD(bhulev)=     Ck_theta(model_levels) * 2.0      &
                                         -    Ck_rho(model_levels)
          ELSE
            PP_REAL_HEAD(bulev) = zseak_rho(int_level+1)
            PP_REAL_HEAD(bhulev)=    Ck_rho(int_level+1)
          ENDIF                             ! top level

          PP_REAL_HEAD(blev) = zseak_theta(int_level)
          PP_REAL_HEAD(bhlev)=    Ck_theta(int_level)

! Note that the lowest theta level (=1) has a lower layer boundary at
! the surface, as set explicitly in the model interface to physics.
          IF(int_level <= 1) THEN            ! bottom level
             PP_REAL_HEAD(brlev) = 0.     ! zsea at/below surface
             PP_REAL_HEAD(bhrlev)= 1.     ! C    at/below surface
          ELSE
            PP_REAL_HEAD(brlev) = zseak_rho(int_level)
            PP_REAL_HEAD(bhrlev)=    Ck_rho(int_level)
          ENDIF                              ! bottom level

        ELSEIF(lvcode == ppx_full_level) THEN ! rho level (& u,v,p)

          IF(int_level >  model_levels) THEN ! p above top level
            PP_REAL_HEAD(bulev) = zseak_theta(model_levels) * 2.0       &
                                        - zseak_rho(model_levels)
            PP_REAL_HEAD(bhulev)=    Ck_theta(model_levels) * 2.0       &
                                         -   Ck_rho(model_levels)
            PP_REAL_HEAD(blev) = PP_REAL_HEAD(bulev)
           PP_REAL_HEAD(bhlev)= PP_REAL_HEAD(bhulev)
          ELSE
            PP_REAL_HEAD(bulev) = zseak_theta(int_level)
            PP_REAL_HEAD(bhulev)=    Ck_theta(int_level)
            PP_REAL_HEAD(blev)  = zseak_rho(int_level)
            PP_REAL_HEAD(bhlev) =    Ck_rho(int_level)
          ENDIF                              ! p above top level

          IF(int_level <= 0) THEN            ! bottom level
            PP_REAL_HEAD(brlev) = 0.    ! zsea at/below surface
            PP_REAL_HEAD(bhrlev)= 1.    ! C    at/below surface
          ELSE
            PP_REAL_HEAD(brlev) = zseak_theta(int_level-1)
            PP_REAL_HEAD(bhrlev)=    Ck_theta(int_level-1)
          ENDIF                              ! bottom level

        ELSE                ! Illegal lvcode
          ICODE=1
          CMESSAGE=' Inconsistent vertical coordinate codes in '//&
              'STASHmaster for this output field: LevelT indicates '//&
              'that this variable is on model levels, but LBVC used '//&
              'for pp header label does not.'
          WRITE(6,*) RoutineName,CMESSAGE                               &
                  ,' section,item,LevelT,pp_lbvc=',is,ie,lvcode,pp_lbvc

        ENDIF               ! Test on lvcode
      ELSE
        PP_REAL_HEAD(BLEV)=real_level
        PP_REAL_HEAD(BHLEV)=0.0
        PP_REAL_HEAD(BRLEV)=0.0  ! The boundary levels
        PP_REAL_HEAD(BHRLEV)=0.0 ! are not known
        PP_REAL_HEAD(BULEV)=0.0  ! for pressure
        PP_REAL_HEAD(BHULEV)=0.0 ! levels.
      ENDIF
!
!  Horizontal grid description -----------------------------------------
!  Position of pole (from dump fixed-length header)
!  Grid orientation (hardwired 0.0)
!  Origin and spacing of grid (depends on output grid type)
!
      PP_REAL_HEAD(BPLAT)=REALHD(5)
      PP_REAL_HEAD(BPLON)=REALHD(6)
      PP_REAL_HEAD(BGOR)=0.0
      IF (samples >  0) THEN   ! Indicates a timeseries/trajectory
        PP_REAL_HEAD(BZX)=0.0
        PP_REAL_HEAD(BDX)=0.0
        PP_REAL_HEAD(BZY)=0.0
        PP_REAL_HEAD(BDY)=0.0
      ELSE
        IF(st_grid == st_riv_grid)THEN
          PP_REAL_HEAD(BDY) = 180.0/PP_INT_HEAD(LBROW)
          PP_REAL_HEAD(BZY) = REALHD(3) - PP_REAL_HEAD(BDY)*0.5
          PP_REAL_HEAD(BDX) = 360.0/PP_INT_HEAD(LBNPT)
          PP_REAL_HEAD(BZX) = REALHD(4) - PP_REAL_HEAD(BDX)*0.5
        ELSE
        IF (l_vatpoles) THEN
          IF(st_grid == st_uv_grid.OR.st_grid == st_cv_grid .OR.        &
              st_grid == st_zu_grid.OR.st_grid == st_mu_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)    ! V pts
          ELSE
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)/2.0   ! Zeroth Lat BZY
          ENDIF
          
          IF(st_grid == st_uv_grid.OR.st_grid == st_cu_grid .OR.        &
              st_grid == st_zu_grid.OR.st_grid == st_mu_grid) THEN
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)   ! U points
          ELSE
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)/2.0   ! Zeroth Long BZX
          ENDIF
        ELSE
          IF(st_grid == st_uv_grid.OR.st_grid == st_cv_grid .OR.        &
              st_grid == st_zu_grid.OR.st_grid == st_mu_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)/2.0 ! UV pts
          ELSE
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2) ! Zeroth Lat BZY
          ENDIF
          
          IF(st_grid == st_uv_grid.OR.st_grid == st_cu_grid .OR.        &
              st_grid == st_zu_grid.OR.st_grid == st_mu_grid) THEN
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)/2.0 !UV points
          ELSE
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1) ! Zeroth Long BZX
          ENDIF
        ENDIF
        PP_REAL_HEAD(BDX)=REALHD(1) ! Long intvl BDX
        PP_REAL_HEAD(BDY)=REALHD(2) ! Lat intvl BDY
      END IF
!
! Add on offset for fields not starting from the origin (sub-areas)
!
      ! If this is variable horizontal grid information
      ! then the horizontal start points are NOT equally spaced.
      ! The existing code is meaningless in this case so
      ! we must get our actual start point some other way!
        IF (VAR_GRID_TYPE >  0) THEN
          IF (X_VAR_GRID) PP_REAL_HEAD(BZX) =                           &
              X_BOUNDARY(WCOL_IN,VAR_GRID_TYPE)
          IF (Y_VAR_GRID) PP_REAL_HEAD(BZY) =                           &
              Y_BOUNDARY(SROW_IN,VAR_GRID_TYPE)
        ELSE
          IF(st_grid == st_uv_grid.OR.st_grid == st_cv_grid .OR.        &
              st_grid == st_zu_grid.OR.st_grid == st_mu_grid) THEN
            IF(SROW_IN /= 1) THEN
! v area shifted up a row so diagnostic same as pre new dynamics
! when the rows were reversed. ie. N-S
              PP_REAL_HEAD(BZY)=PP_REAL_HEAD(BZY)                       &
                  +(SROW_IN-2)*PP_REAL_HEAD(BDY)
            END IF
          ELSE
            PP_REAL_HEAD(BZY)=PP_REAL_HEAD(BZY)                         &
                +(SROW_IN-1)*PP_REAL_HEAD(BDY)
          END IF
          PP_REAL_HEAD(BZX)=PP_REAL_HEAD(BZX)                           &
              +(WCOL_IN-1)*PP_REAL_HEAD(BDX)
          
        ENDIF
        IF(PP_REAL_HEAD(BZX) >= 360.0) THEN
          PP_REAL_HEAD(BZX)=PP_REAL_HEAD(BZX)-360.0
        ENDIF
        
!
! If horizontal averaging has been applied to the output field,
! set BDX and/or BDY to the full (sub)domain extent which was processed.
! If the input field was intrinsically non-2D (eg. zonal), assume that
! the collapsed dimension(s) covered the full model domain.
!
        mean_code=(GR/block_size)*block_size
        IF (st_grid == st_zt_grid .OR. st_grid == st_zu_grid            &
            .OR. st_grid == st_scalar) THEN
          PP_REAL_HEAD(BDX)=REAL(INTHD(6))*PP_REAL_HEAD(BDX)
        ELSEIF (mean_code == zonal_mean_base .OR.                       &
            mean_code == field_mean_base .OR.                           &
            mean_code == global_mean_base) THEN
          PP_REAL_HEAD(BDX)=ABS(REAL(ECOL_IN-WCOL_IN))*PP_REAL_HEAD(BDX)
        ENDIF
!
        IF (st_grid == st_mt_grid .OR. st_grid == st_mu_grid            &
            .OR. st_grid == st_scalar) THEN
          PP_REAL_HEAD(BDY)=REAL(INTHD(7))*PP_REAL_HEAD(BDY)
        ELSEIF (mean_code == merid_mean_base .OR.                       &
            mean_code == field_mean_base .OR.                           &
            mean_code == global_mean_base) THEN
          PP_REAL_HEAD(BDY)=ABS(REAL(NROW_IN-SROW_IN))*PP_REAL_HEAD(BDY)

        ENDIF
      ENDIF
      IF( (FIXHD(117) > 1 .AND. FIXHD(122) > 1)) THEN 
        PP_REAL_HEAD(BDX) = RMDI
        PP_REAL_HEAD(BDY) = RMDI
        PP_REAL_HEAD(BZX) = RMDI
        PP_REAL_HEAD(BZY) = RMDI
      ENDIF
!
! Missing data indicator (from PARAMETER) ------------------------------
! MKS scaling factor (unity as model uses SI units throughout)
!
      PP_REAL_HEAD(BMDI)=RMDI
      PP_REAL_HEAD(BMKS)=1.0
!

      IF (lhook) CALL dr_hook('PP_HEAD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE PP_HEAD
