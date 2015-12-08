! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! *********************************************************************
!+ SUBROUTINE set_ppindex : Set up addressing of D1 related arrays
!                           such as ppindex 
!
! Subroutine Interface :
!
       SUBROUTINE set_ppindex(                                           &
!   Scalar arguments with intent(in)     
                 len_ppindex, len1_lookup, len2_lookup,                  &
                 tr_levels,no_tracers,                                   &
                 no_tr_ukca,max_progs,l_pc2, l_murk, l_mcr_qcf2,         &
                 l_mcr_qrain, l_mcr_qgraup,                              &
                 l_dust,l_so2,l_so4_aitken,l_so4_accu,l_so4_diss,l_dms,  &
                 l_nh3,l_soot,l_biomass,l_ocff,l_nitrate,                &
                 target_time,                                            &
!   Array arguments with intent(in) :
                 lookup, tracers_active, tr_ukca_active,                 &
!   Scalar arguments with intent(out) :
                 jorog, ju, jv, jw, jrho, jtheta, jq, jqcl, jqcf,        &
                 jexner, ju_adv, jv_adv, jw_adv, jqcf2, jqrain, jqgraup, &
                 jmurk, jcf_bulk, jcf_liquid, jcf_frozen,                &
                 jdust_div1, jdust_div2, jdust_div3, jdust_div4,         &
                 jdust_div5, jdust_div6, jso2, jso4_aitken,              &
                 jso4_accu, jso4_diss, jdms, jnh3, jsoot_new,            &
                 jsoot_agd, jsoot_cld, jbmass_new, jbmass_agd,           &
                 jbmass_cld, jocff_new, jocff_agd, jocff_cld,            &
                 jnitr_acc, jnitr_diss,                                  &
                 calc_length, next_dump,                                 &
!   Array arguments with intent(out) :
                 ppindex, pplengths, lbcdiag_no, d1_pointers,            &
                 halo_size_out,jtracer, jtr_ukca,                        &
                 ppindex_tracer, ppindex_tr_ukca,                        &
                 pplengths_tracer, pplengths_tr_ukca, ndustbin_in,       &
                 ndustbin_out)
!    
! Type for time
      USE makebc_time_mod, ONLY: &
        time
      USE makebc_constants_mod
        
      USE UM_ParVars
      USE Decomp_DB
      USE Submodel_Mod
      IMPLICIT NONE

! Description : Initialise ppindex array and pointers.
!
! Method : This subroutine sets initial value then calls subroutines
!          which go through the lookup table from the input  
!          dump, checks which fields have the correct time and date
!          and sets up ppindex and pplengths. The j pointers 
!          and d1_pointers are then calculated from the sizes of the
!          data. 
!          If there is no data with the correct date 
!          in the dump's lookup table it will set the logical 
!          next_dump to .TRUE. This logical variable controls the
!          loop over times in get_bc which calls this subroutine
!          If there is data in the dump for some but not all
!          required fields then this gives a fatal error
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC
!
! Code Description :
! Language : FORTRAN 90
! This code is written to UMDP3 v8.2 programming standards.
!
! Declarations :
!
! Global Variables :

! Include file : parlbcs.h
!
! Must be called after parvars.h
!
! Description:
!   Contains variables in connection with generating LBCs.
!
! -----------------------------------------------------------
! Stash Codes for LBCs in Section 32 (and section 31 except tracers) 
!
      Integer, Parameter :: lbc_stashcode_orog    = 1 
      Integer, Parameter :: lbc_stashcode_u       = 2 
      Integer, Parameter :: lbc_stashcode_v       = 3 
      Integer, Parameter :: lbc_stashcode_w       = 4 
      Integer, Parameter :: lbc_stashcode_density = 5 
      Integer, Parameter :: lbc_stashcode_theta   = 6 
      Integer, Parameter :: lbc_stashcode_q       = 7 
      Integer, Parameter :: lbc_stashcode_qcl     = 8 
      Integer, Parameter :: lbc_stashcode_qcf     = 9 
      Integer, Parameter :: lbc_stashcode_exner   = 10 
      Integer, Parameter :: lbc_stashcode_u_adv   = 11 
      Integer, Parameter :: lbc_stashcode_v_adv   = 12 
      Integer, Parameter :: lbc_stashcode_w_adv   = 13 
      Integer, Parameter :: lbc_stashcode_qcf2    = 14 
      Integer, Parameter :: lbc_stashcode_qrain   = 15 
      Integer, Parameter :: lbc_stashcode_qgraup  = 16 
      Integer, Parameter :: lbc_stashcode_cf_bulk = 17 
      Integer, Parameter :: lbc_stashcode_cf_liquid = 18 
      Integer, Parameter :: lbc_stashcode_cf_frozen = 19 
      Integer, Parameter :: lbc_stashcode_murk      = 20 
      Integer, Parameter :: lbc_stashcode_free_tracer = 21 
      Integer, Parameter :: lbc_stashcode_ukca_tracer = 22 
      Integer, Parameter :: lbc_stashcode_dust_div1 = 23
      Integer, Parameter :: lbc_stashcode_dust_div2 = 24
      Integer, Parameter :: lbc_stashcode_dust_div3 = 25
      Integer, Parameter :: lbc_stashcode_dust_div4 = 26
      Integer, Parameter :: lbc_stashcode_dust_div5 = 27
      Integer, Parameter :: lbc_stashcode_dust_div6 = 28
      Integer, Parameter :: lbc_stashcode_so2      = 29
      Integer, Parameter :: lbc_stashcode_dms      = 30
      Integer, Parameter :: lbc_stashcode_so4_aitken = 31
      Integer, Parameter :: lbc_stashcode_so4_accu = 32
      Integer, Parameter :: lbc_stashcode_so4_diss = 33
      Integer, Parameter :: lbc_stashcode_nh3      = 35
      Integer, Parameter :: lbc_stashcode_soot_new = 36
      Integer, Parameter :: lbc_stashcode_soot_agd = 37
      Integer, Parameter :: lbc_stashcode_soot_cld = 38
      Integer, Parameter :: lbc_stashcode_bmass_new = 39
      Integer, Parameter :: lbc_stashcode_bmass_agd = 40
      Integer, Parameter :: lbc_stashcode_bmass_cld = 41
      Integer, Parameter :: lbc_stashcode_ocff_new = 42
      Integer, Parameter :: lbc_stashcode_ocff_agd = 43
      Integer, Parameter :: lbc_stashcode_ocff_cld = 44
      Integer, Parameter :: lbc_stashcode_nitr_acc = 45
      Integer, Parameter :: lbc_stashcode_nitr_diss = 46

! -----------------------------------------------------------
!     Data Time for LBC data
      Integer :: LBC_DT_Year
      Integer :: LBC_DT_Month
      Integer :: LBC_DT_Day
      Integer :: LBC_DT_Hour
      Integer :: LBC_DT_Min
      Integer :: LBC_DT_Sec
      Integer :: LBC_DT_DayNo

      COMMON /LBC_DT/ LBC_DT_Year, LBC_DT_Month, LBC_DT_Day,            &
         LBC_DT_Hour, LBC_DT_Min,  LBC_DT_Sec,   LBC_DT_DayNo

! -----------------------------------------------------------

!     Validity Time for LBC data
      Integer :: LBC_VT_Year
      Integer :: LBC_VT_Month
      Integer :: LBC_VT_Day
      Integer :: LBC_VT_Hour
      Integer :: LBC_VT_Min
      Integer :: LBC_VT_Sec
      Integer :: LBC_VT_DayNo

      COMMON /LBC_VT/ LBC_VT_Year, LBC_VT_Month, LBC_VT_Day,            &
         LBC_VT_Hour, LBC_VT_Min,  LBC_VT_Sec,   LBC_VT_DayNo

! -----------------------------------------------------------

      Integer, Parameter :: P_Src_Grid = 2
      Integer, Parameter :: P_LBC_Grid = 4

!     1 : Start Latitude
!     2 : Start Longitude
!     3 : Row Length
!     4 : Rows

      Real :: Src_Grid (Nfld_max, P_Src_Grid)
      Real :: LBC_Grid (Nfld_max, P_LBC_Grid)

      COMMON /LBC_Grids/ Src_Grid, LBC_Grid

! -------------------------------------------------------------

      Integer :: LBC_Global_LenRimA (Nfld_max, Nhalo_max)
      Integer :: LBC_Interp_LenRimA (Nfld_max, Nhalo_max)

      COMMON /LBC_Sizes/ LBC_Global_LenRimA, LBC_Interp_LenRimA

! -------------------------------------------------------------
! CTRACERA start
!  Vn    Date    Modification History
! 6.1  23/06/04  Prognostic tracers now in section 33, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.
! 6.2  13/07/05  Also increase A_MAX_TRVARS to 150. R Barnes.
! 6.2  10/11/05  UKCA tracers put into section 34, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.

      ! First atmospheric tracer (STASH No)
      INTEGER,PARAMETER:: A_TRACER_FIRST = 1
      !First UKCA tracer (STASH No)
      INTEGER,PARAMETER:: A_UKCA_FIRST = 1

      ! Last atmospheric tracer  (STASH No)
      INTEGER,PARAMETER:: A_TRACER_LAST = 150
      !Last UKCA tracer  (STASH No)
      INTEGER,PARAMETER:: A_UKCA_LAST = 150

      ! Maximum number of atmospheric tracers
      INTEGER,PARAMETER:: A_MAX_TRVARS  = 150
      !Maximum number of UKCA tracers
      INTEGER,PARAMETER:: A_MAX_UKCAVARS  = 150

      ! Index to relative position.
      ! A_TR_INDEX(N) gives position in JTRACER for tracer number N.
      ! Set in SET_ATM_POINTERS.
      ! A_TR_INDEX(N) is the position, in the list of tracers
      ! actually present in D1, that tracer number N (in the list
      ! of all tracers selectable from the user interface) occupies,
      ! if it is present.
      ! If tracer number N is absent then A_TR_INDEX(N) is -1.
      ! Similarly for A_UKCA_INDEX.

      INTEGER :: A_TR_INDEX(A_MAX_TRVARS)
      ! A_TR_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: A_TR_StashItem(A_MAX_TRVARS)

      INTEGER :: A_UKCA_INDEX(A_MAX_UKCAVARS)
      ! UKCA_tr_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: UKCA_tr_StashItem(A_MAX_UKCAVARS) 

      ! A_TR_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: A_TR_LBC_StashItem(A_MAX_TRVARS) 
      INTEGER :: A_TR_active_lbc_index(A_MAX_TRVARS) 

      ! UKCA_tr_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: UKCA_tr_LBC_StashItem(A_MAX_UKCAVARS) 
      INTEGER :: UKCA_tr_active_lbc_index(A_MAX_UKCAVARS)

      COMMON/ATRACER/A_TR_INDEX, A_TR_StashItem,                        &
     &               A_TR_LBC_StashItem, A_TR_active_lbc_index,         &
     &               A_UKCA_INDEX, UKCA_tr_StashItem,                   &
     &               UKCA_tr_LBC_StashItem, UKCA_tr_active_lbc_index

! CTRACERA end

! *********************************************************************
! Subroutine arguments
!   Scalar arguments with intent(in) :

      INTEGER, INTENT(IN) :: len_ppindex   ! Dimension of pp_index
      INTEGER, INTENT(IN) :: len1_lookup   ! First dimension of Lookup table
      INTEGER, INTENT(IN) :: len2_lookup   ! Second dimension of Lookup table

      INTEGER, INTENT(IN) :: tr_levels     ! The number of tracer levels
      INTEGER, INTENT(IN) :: no_tracers    ! The number of active tracers
      INTEGER, INTENT(IN) :: no_tr_ukca    ! The number of active tracers

! max no of prognostics
      INTEGER, INTENT(IN) ::   max_progs

      LOGICAL, INTENT(IN) :: l_pc2
      LOGICAL, INTENT(IN) :: l_murk
      LOGICAL, INTENT(IN) :: l_mcr_qcf2
      LOGICAL, INTENT(IN) :: l_mcr_qrain
      LOGICAL, INTENT(IN) :: l_mcr_qgraup
      LOGICAL, INTENT(IN) :: l_dust
      LOGICAL, INTENT(IN) :: l_so2
      LOGICAL, INTENT(IN) :: l_so4_aitken
      LOGICAL, INTENT(IN) :: l_so4_accu
      LOGICAL, INTENT(IN) :: l_so4_diss
      LOGICAL, INTENT(IN) :: l_dms
      LOGICAL, INTENT(IN) :: l_nh3
      LOGICAL, INTENT(IN) :: l_soot
      LOGICAL, INTENT(IN) :: l_biomass
      LOGICAL, INTENT(IN) :: l_ocff
      LOGICAL, INTENT(IN) :: l_nitrate

! The input variable target_time is used to check the 
! date/time of fields in the lookup tables
      TYPE(TIME),INTENT(IN) :: target_time

!   Array arguments with intent(in) :
!
! Lookup table from dump being read
      INTEGER, INTENT(IN) :: lookup(len1_lookup,len2_lookup)  
!
      INTEGER, INTENT(IN) :: tracers_active(max_tracers)
      INTEGER, INTENT(IN) :: tr_ukca_active(max_tr_ukca)

! Number of dust bins in the input dump (in) and output LBCs (out)
      INTEGER, INTENT(IN) :: ndustbin_in, ndustbin_out

!   Scalar arguments with intent(out) :
!   Pointers 
      INTEGER, INTENT(OUT) :: jorog      !pointer to orography
      INTEGER, INTENT(OUT) :: ju         !pointer to u wind component
      INTEGER, INTENT(OUT) :: jv         !pointer to v wind component
      INTEGER, INTENT(OUT) :: jw         !pointer to w wind component (vert.)
      INTEGER, INTENT(OUT) :: jrho       !pointer to density
      INTEGER, INTENT(OUT) :: jtheta     !pointer to theta
      INTEGER, INTENT(OUT) :: jq         !pointer to specific humidity
      INTEGER, INTENT(OUT) :: jqcl       !pointer to qcl liquid water
      INTEGER, INTENT(OUT) :: jqcf       !pointer to qcf frozen water
      INTEGER, INTENT(OUT) :: jexner     !pointer to exner pressure
      INTEGER, INTENT(OUT) :: ju_adv     !pointer to u advection
      INTEGER, INTENT(OUT) :: jv_adv     !pointer to v advection
      INTEGER, INTENT(OUT) :: jw_adv     !pointer to w advection
!   pc2 and murk pointers
      INTEGER, INTENT(OUT) :: jqcf2      !pointer to qcf2 - cloud ice (cry)
      INTEGER, INTENT(OUT) :: jqrain     !pointer to qrain - rain
      INTEGER, INTENT(OUT) :: jqgraup    !pointer to qgraup - graupel horiz.
      INTEGER, INTENT(OUT) :: jmurk      !pointer to murk
      INTEGER, INTENT(OUT) :: jcf_bulk   !pointer to cf_bulk
      INTEGER, INTENT(OUT) :: jcf_liquid !pointer to cf_liquid
      INTEGER, INTENT(OUT) :: jcf_frozen !pointer to cf_frozen
!   Aerosol pointers 
      INTEGER, INTENT(OUT) :: jdust_div1      ! dust mmr, division 1
      INTEGER, INTENT(OUT) :: jdust_div2      ! dust mmr, division 2
      INTEGER, INTENT(OUT) :: jdust_div3      ! dust mmr, division 3
      INTEGER, INTENT(OUT) :: jdust_div4      ! dust mmr, division 4
      INTEGER, INTENT(OUT) :: jdust_div5      ! dust mmr, division 5
      INTEGER, INTENT(OUT) :: jdust_div6      ! dust mmr, division 6
      INTEGER, INTENT(OUT) :: jso2            ! sulphur dioxide gas
      INTEGER, INTENT(OUT) :: jso4_aitken     ! Aitken mode sulphate aer
      INTEGER, INTENT(OUT) :: jso4_accu       ! accumulation mode sulpha
      INTEGER, INTENT(OUT) :: jso4_diss       ! dissloved  sulphate aero
      INTEGER, INTENT(OUT) :: jdms            ! dimethyl sulphide gas
      INTEGER, INTENT(OUT) :: jnh3            ! ammonia gas mmr
      INTEGER, INTENT(OUT) :: jsoot_new       ! fresh soot mmr
      INTEGER, INTENT(OUT) :: jsoot_agd       ! aged soot mmr
      INTEGER, INTENT(OUT) :: jsoot_cld       ! soot in cloud mmr
      INTEGER, INTENT(OUT) :: jbmass_new      ! fresh biomass mmr
      INTEGER, INTENT(OUT) :: jbmass_agd      ! aged biomass mmr
      INTEGER, INTENT(OUT) :: jbmass_cld      ! cloud biomass mmr
      INTEGER, INTENT(OUT) :: jocff_new       ! fresh OCFF mmr
      INTEGER, INTENT(OUT) :: jocff_agd       ! aged OCFF mmr
      INTEGER, INTENT(OUT) :: jocff_cld       ! OCFF in cloud mmr
      INTEGER, INTENT(OUT) :: jnitr_acc       ! Accumulation nitrate aerosol
      INTEGER, INTENT(OUT) :: jnitr_diss      ! Dissolved nitrate aerosol

!   calc_length is the length required to store all prognostics
      INTEGER, INTENT(OUT) ::   calc_length
!
!   logical for get_bc to exit and read next dump      
      LOGICAL, INTENT(OUT) ::   next_dump

!   Array arguments with intent(out) :

! Index of section 0 fields in fields file i.e. ppindex(2,1) 
! gives the first U field in the field file (stash code for
! U is 2)
      INTEGER, INTENT(OUT) :: ppindex(len_ppindex,n_internal_model)     

! Array pplengths allows setppindex to generate extra values required
! by makebc to use readflds to read the dump with um_readdump
! It is the number of 2D fields for each prognostic in the dump (levels)
! It is named num_levels in get_bc
      INTEGER, INTENT(OUT) :: pplengths(len_ppindex,n_internal_model)   

! lbcdiag_no is item number for each prognostic if in section 0 
!
! for section 33/34 tracers however it the number in tracers_active
! or tr_ukca_active
!
! i.e. if tracers_active are 1 3 10 and 11
! then lbcdiag_no(lbcdiag_no_tracerstart) will be 1 
! and lbcdiag_no(lbcdiag_no_tracerstart+1) will be 2
! and lbcdiag_no(lbcdiag_no_tracerstart+3) will be 4
!
! If there are no tracers lbcdiag_no will end after the section zero LBCs
!
! If it is -1 on exit the field is not in the lookup header of the fields
! file
      INTEGER, INTENT(OUT) :: lbcdiag_no(max_progs)
! 
! d1_pointers - the position of each diagnostic in the fields file.
      INTEGER, INTENT(OUT) :: d1_pointers(max_progs)

! halo_size_out - halo size of all prognostics
      INTEGER, INTENT(OUT) :: halo_size_out(max_progs,3)

! Pointer to free and UKCA tracer(s)
      INTEGER, INTENT(OUT) :: jtracer(tr_levels,no_tracers+1)  
      INTEGER, INTENT(OUT) :: jtr_ukca(tr_levels,no_tr_ukca+1)  
! pplengths and ppindex for free and UKCA tracers
      INTEGER, INTENT(OUT) :: ppindex_tracer(no_tracers,1)
      INTEGER, INTENT(OUT) :: ppindex_tr_ukca(no_tr_ukca,1)
      INTEGER, INTENT(OUT) :: pplengths_tracer(no_tracers,1)
      INTEGER, INTENT(OUT) :: pplengths_tr_ukca(no_tr_ukca,1)

! *********************************************************************
!   Local variables

! Set first element of lbcdiag_no to be used for free tracers as parameter
! We always allow enough space for all optional LBCs even if not in use
! Whereas the size of tracer section only large enough for tracers requested
!
      INTEGER, PARAMETER :: lbcdiag_no_tracerstart = &
          num_req_lbcs+num_optional_lbcs_max+1

!
! tracers

      INTEGER :: lbcdiag_no_tracer
      INTEGER :: lbcdiag_no_tracerend


! Set first element of lbcdiag_no to be used for UKCA tracers as parameter
      INTEGER :: lbcdiag_no_trstart_ukca
      INTEGER :: lbcdiag_no_tr_ukca
      INTEGER :: lbcdiag_no_trend_ukca

! Logical for whether data is present for target_time
      LOGICAL :: data_present

!-  End of Header
! *********************************************************************
!
! Set initial values of scalars and arrays
!
! Initialise calc_length to 0 in case not done before 
      calc_length=0

! Initialise d1_pointers to -1
      d1_pointers(:)=-1

! Initialise data_present to .false.
      data_present=.false.
      ju          =-1
      jv          =-1
      jtheta      =-1
      jq          =-1
      jqcf        =-1
      jorog       =-1
      jmurk       =-1
      jw          =-1
      jrho        =-1
      jqcl        =-1
      jexner      =-1
      ju_adv      =-1
      jv_adv      =-1
      jw_adv      =-1
      jcf_bulk    =-1
      jcf_liquid  =-1
      jcf_frozen  =-1
      jqcf2       =-1
      jqrain      =-1
      jqgraup     =-1

! Initialise aerosol pointers
      jdust_div1  =-1
      jdust_div2  =-1
      jdust_div3  =-1
      jdust_div4  =-1
      jdust_div5  =-1
      jdust_div6  =-1
      jso2        =-1
      jso4_aitken =-1
      jso4_accu   =-1
      jso4_diss   =-1
      jdms        =-1
      jnh3        =-1
      jsoot_new   =-1
      jsoot_agd   =-1
      jsoot_cld   =-1
      jbmass_new  =-1
      jbmass_agd  =-1
      jbmass_cld  =-1
      jocff_new   =-1
      jocff_agd   =-1
      jocff_cld   =-1
      jnitr_acc   =-1
      jnitr_diss  =-1

! Initialise tracers
      jtracer(:,:)    =-1
      jtr_ukca(:,:)   =-1
! Give ppindex_tracer and ppindex_tr_ukca so can later use to test for missing data 
      ppindex_tracer  =-1
      ppindex_tr_ukca =-1
      
! Set variable lbcdiag_no_tracer from parameter lbcdiag_no_tracerstart
! to allow dynamic setting of tracer portion of lbcdiag_no

      lbcdiag_no_tracer=lbcdiag_no_tracerstart
      lbcdiag_no_tracerend=lbcdiag_no_tracerstart+no_tracers-1
      
      lbcdiag_no_trstart_ukca= lbcdiag_no_tracerend+1
      lbcdiag_no_tr_ukca=lbcdiag_no_trstart_ukca
      
      lbcdiag_no_trend_ukca=lbcdiag_no_trstart_ukca+no_tr_ukca

! Initialise pplengths such that all elements are set to 0
      pplengths(:,:)=0
      pplengths_tracer(:,:)=0
      pplengths_tr_ukca(:,:)=0

! Initialise lbcdiag_no array as -1
      lbcdiag_no(:)=-1

! Set values of lbcdiag to stash_xxx parameters
! DEPENDS ON: makebc_set_lbcdiag
      CALL makebc_set_lbcdiag(max_progs,ndustbin_in,lbcdiag_no)

! *********************************************************************
! makebc_ppindex_loop loops over all lookup tables in the array lookup
! Get the starting point in each field which is required, and also set
! pplengths to the number each field (levels) and set the jpointer 
!
!
! DEPENDS ON: makebc_ppindex_loop
      CALL makebc_ppindex_loop(max_progs,len_ppindex,len1_lookup,  &
           len2_lookup, no_tracers,no_tr_ukca,lookup,              &
           tracers_active, tr_ukca_active, target_time,            &
           jorog, ju, jv, jw, jrho, jtheta, jq, jqcl, jqcf,        &
           jexner, ju_adv, jv_adv, jw_adv, jqcf2, jqrain, jqgraup, &
           jmurk, jcf_bulk, jcf_liquid, jcf_frozen,                &
           jdust_div1, jdust_div2, jdust_div3, jdust_div4,         &
           jdust_div5, jdust_div6, jso2, jso4_aitken,              &
           jso4_accu, jso4_diss, jdms, jnh3, jsoot_new,            &
           jsoot_agd, jsoot_cld, jbmass_new, jbmass_agd,           &
           jbmass_cld, jocff_new, jocff_agd, jocff_cld,            &
           jnitr_acc, jnitr_diss,lbcdiag_no_tracer,                &
           lbcdiag_no_tr_ukca,data_present, lbcdiag_no,            &
           ppindex, pplengths, ppindex_tracer, ppindex_tr_ukca,    &
           pplengths_tracer, pplengths_tr_ukca )

! *********************************************************************
! Test if data found for this time 
      IF (data_present) THEN
        ! on return to get_bc prevent exit from loop over times
        next_dump=.false.    
!
! Check if any required fields are missing 
! DEPENDS ON: makebc_check_present
        CALL makebc_check_present(                                          &
           jorog, ju, jv, jw, jrho, jtheta, jq, jqcl, jqcf,                 &
           jexner, ju_adv, jv_adv, jw_adv, jqcf2, jqrain, jqgraup,          &
           jmurk, jcf_bulk, jcf_liquid, jcf_frozen,                         &
           jdust_div1, jdust_div2, jdust_div3, jdust_div4,                  &
           jdust_div5, jdust_div6, jso2, jso4_aitken,                       &
           jso4_accu, jso4_diss, jdms, jnh3, jsoot_new,                     &
           jsoot_agd, jsoot_cld, jbmass_new, jbmass_agd,                    &
           jbmass_cld, jocff_new, jocff_agd, jocff_cld,                     &
           jnitr_acc, jnitr_diss, tr_levels, no_tracers, no_tr_ukca,        &
           l_pc2, l_murk, l_mcr_qcf2,l_mcr_qrain, l_mcr_qgraup,             &
           l_dust, l_so2, l_so4_aitken, l_so4_accu, l_so4_diss, l_dms,      &
           l_nh3, l_soot, l_biomass, l_ocff, l_nitrate, ppindex_tracer,     &
           ppindex_tr_ukca, tracers_active, tr_ukca_active, ndustbin_in,    &
           ndustbin_out )
!
! call subroutin makebc_setd1point to loop over the prognostic fields
! and set up d1_pointers, jtracer and jtr_ukca
! 
! DEPENDS ON: makebc_setd1point
        CALL makebc_setd1point( max_progs,len_ppindex,no_tracers,  &
           no_tr_ukca,tr_levels, lbcdiag_no_tracerstart,           &
           lbcdiag_no_tracerend, lbcdiag_no_trstart_ukca,          &
           pplengths,pplengths_tracer,pplengths_tr_ukca,           &
           lbcdiag_no, tracers_active, tr_ukca_active,             &
           jorog, ju, jv, jw, jrho, jtheta, jq, jqcl, jqcf,        &
           jexner, ju_adv, jv_adv, jw_adv, jqcf2, jqrain, jqgraup, &
           jmurk, jcf_bulk, jcf_liquid, jcf_frozen,                &
           jdust_div1, jdust_div2, jdust_div3, jdust_div4,         &
           jdust_div5, jdust_div6, jso2, jso4_aitken,              &
           jso4_accu, jso4_diss, jdms, jnh3, jsoot_new,            &
           jsoot_agd, jsoot_cld, jbmass_new, jbmass_agd,           &
           jbmass_cld, jocff_new, jocff_agd, jocff_cld,            &
           jnitr_acc, jnitr_diss, calc_length, d1_pointers,        &
           jtracer,jtr_ukca,halo_size_out)
    
! *********************************************************************
!
      ELSE ! If no data present
         
        next_dump=.TRUE.   ! on return to get_bc exit from loop over times
          
      END IF  ! if data_present = .true.
      
! *********************************************************************
! Avoid any negative pointers for optional LBCs not required
!
! DEPENDS ON: makebc_avoid_neg_point
      CALL makebc_avoid_neg_point(max_progs,l_pc2, l_murk, l_mcr_qcf2,   &
                 l_mcr_qrain, l_mcr_qgraup,                              &
                 l_dust,l_so2,l_so4_aitken,l_so4_accu,l_so4_diss,l_dms,  &
                 l_nh3,l_soot,l_biomass,l_ocff,l_nitrate,                &
                 tr_levels, no_tracers, no_tr_ukca,                      &
                 lbcdiag_no_tracerstart, lbcdiag_no_trstart_ukca,        &
                 jtracer, jtr_ukca,                                      &
                 jqcf2,  jqrain,  jqgraup, jmurk,  jcf_bulk,             &
                 jcf_liquid,  jcf_frozen,  jdust_div1,  jdust_div2,      &
                 jdust_div3,  jdust_div4,  jdust_div5,  jdust_div6,      &
                 jso2,  jso4_aitken,  jso4_accu,  jso4_diss,  jdms,      &
                 jnh3,  jsoot_new,  jsoot_agd,  jsoot_cld,  jbmass_new,  &
                 jbmass_agd,  jbmass_cld,  jocff_new,  jocff_agd,        &
                 jocff_cld,  jnitr_acc,  jnitr_diss, lbcdiag_no)

      END SUBROUTINE set_ppindex
