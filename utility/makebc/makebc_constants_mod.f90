! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
! A module containing MakeBC specific constants
!
MODULE makebc_constants_mod 

IMPLICIT NONE

! Parameters for the number of LBC - required e.g. u,v and optional e.g. murk
! If the number of required LBCs or optional LBCs changes, change these parameters
!
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: MakeBC
INTEGER, PARAMETER :: num_req_lbcs          = 13 ! number of LBCs always required
INTEGER, PARAMETER :: num_optional_lbcs_max = 30 ! max number of optional LBCs 

! Parameters for the stashnumbers to avoid 'magic numbers'
! will make further modification to the code easier.
! If adding new LBCs need to add the stash numbers here
! See also lbc_index below and the above parameters
!
INTEGER, PARAMETER :: stash_ju         = 2
INTEGER, PARAMETER :: stash_jv         = 3
INTEGER, PARAMETER :: stash_jtheta     = 4
INTEGER, PARAMETER :: stash_jq         = 10
INTEGER, PARAMETER :: stash_jqcf       = 12
INTEGER, PARAMETER :: stash_jorog      = 33
INTEGER, PARAMETER :: stash_jmurk      = 90
INTEGER, PARAMETER :: stash_jw         = 150
INTEGER, PARAMETER :: stash_jrho       = 253
INTEGER, PARAMETER :: stash_jqcl       = 254
INTEGER, PARAMETER :: stash_jexner     = 255
INTEGER, PARAMETER :: stash_ju_adv     = 256
INTEGER, PARAMETER :: stash_jv_adv     = 257
INTEGER, PARAMETER :: stash_jw_adv     = 258
INTEGER, PARAMETER :: stash_jcf_bulk   = 266
INTEGER, PARAMETER :: stash_jcf_liquid = 267
INTEGER, PARAMETER :: stash_jcf_frozen = 268
INTEGER, PARAMETER :: stash_jqcf2      = 271
INTEGER, PARAMETER :: stash_jqrain     = 272
INTEGER, PARAMETER :: stash_jqgraup    = 273
!
! Parameters for aerosol related stashnumbers
INTEGER, PARAMETER :: stash_jdust_div1  = 431
INTEGER, PARAMETER :: stash_jdust_div2  = 432
INTEGER, PARAMETER :: stash_jdust_div3  = 433
INTEGER, PARAMETER :: stash_jdust_div4  = 434
INTEGER, PARAMETER :: stash_jdust_div5  = 435
INTEGER, PARAMETER :: stash_jdust_div6  = 436
INTEGER, PARAMETER :: stash_jso2        = 101
INTEGER, PARAMETER :: stash_jso4_aitken = 103
INTEGER, PARAMETER :: stash_jso4_accu   = 104
INTEGER, PARAMETER :: stash_jso4_diss   = 105
INTEGER, PARAMETER :: stash_jdms        = 102
INTEGER, PARAMETER :: stash_jnh3        = 107
INTEGER, PARAMETER :: stash_jsoot_new   = 108
INTEGER, PARAMETER :: stash_jsoot_agd   = 109
INTEGER, PARAMETER :: stash_jsoot_cld   = 110
INTEGER, PARAMETER :: stash_jbmass_new  = 111
INTEGER, PARAMETER :: stash_jbmass_agd  = 112
INTEGER, PARAMETER :: stash_jbmass_cld  = 113
INTEGER, PARAMETER :: stash_jocff_new   = 114
INTEGER, PARAMETER :: stash_jocff_agd   = 115
INTEGER, PARAMETER :: stash_jocff_cld   = 116
INTEGER, PARAMETER :: stash_jnitr_acc   = 117
INTEGER, PARAMETER :: stash_jnitr_diss  = 118

!
! Indices for lcdiag_no and d1_pointers arrays
! If adding new LBCs need to give it an index here 
!
INTEGER, PARAMETER :: lbcindex_u          = 1
INTEGER, PARAMETER :: lbcindex_v          = 2
INTEGER, PARAMETER :: lbcindex_theta      = 3
INTEGER, PARAMETER :: lbcindex_q          = 4
INTEGER, PARAMETER :: lbcindex_qcf        = 5
INTEGER, PARAMETER :: lbcindex_orog       = 6
INTEGER, PARAMETER :: lbcindex_murk       = 7
INTEGER, PARAMETER :: lbcindex_w          = 8
INTEGER, PARAMETER :: lbcindex_rho        = 9
INTEGER, PARAMETER :: lbcindex_qcl        = 10
INTEGER, PARAMETER :: lbcindex_exner      = 11
INTEGER, PARAMETER :: lbcindex_u_adv      = 12
INTEGER, PARAMETER :: lbcindex_v_adv      = 13
INTEGER, PARAMETER :: lbcindex_w_adv      = 14   
INTEGER, PARAMETER :: lbcindex_cf_bulk    = 15
INTEGER, PARAMETER :: lbcindex_cf_liquid  = 16
INTEGER, PARAMETER :: lbcindex_cf_frozen  = 17
INTEGER, PARAMETER :: lbcindex_qcf2       = 18
INTEGER, PARAMETER :: lbcindex_qrain      = 19
INTEGER, PARAMETER :: lbcindex_qgraup     = 20
INTEGER, PARAMETER :: lbcindex_dust_div1  = 21
INTEGER, PARAMETER :: lbcindex_dust_div2  = 22
INTEGER, PARAMETER :: lbcindex_dust_div3  = 23
INTEGER, PARAMETER :: lbcindex_dust_div4  = 24
INTEGER, PARAMETER :: lbcindex_dust_div5  = 25
INTEGER, PARAMETER :: lbcindex_dust_div6  = 26
INTEGER, PARAMETER :: lbcindex_so2        = 27
INTEGER, PARAMETER :: lbcindex_so4_aitken = 28
INTEGER, PARAMETER :: lbcindex_so4_accu   = 29
INTEGER, PARAMETER :: lbcindex_so4_diss   = 30
INTEGER, PARAMETER :: lbcindex_dms        = 31
INTEGER, PARAMETER :: lbcindex_nh3        = 32
INTEGER, PARAMETER :: lbcindex_soot_new   = 33
INTEGER, PARAMETER :: lbcindex_soot_agd   = 34
INTEGER, PARAMETER :: lbcindex_soot_cld   = 35
INTEGER, PARAMETER :: lbcindex_bmass_new  = 36
INTEGER, PARAMETER :: lbcindex_bmass_agd  = 37
INTEGER, PARAMETER :: lbcindex_bmass_cld  = 38
INTEGER, PARAMETER :: lbcindex_ocff_new   = 39
INTEGER, PARAMETER :: lbcindex_ocff_agd   = 40
INTEGER, PARAMETER :: lbcindex_ocff_cld   = 41
INTEGER, PARAMETER :: lbcindex_nitr_acc   = 42
INTEGER, PARAMETER :: lbcindex_nitr_diss  = 43

! Information on lookup headers from UMDP F3. Set up as parameters
!  1 LBYR Year 
!  2 LBMON Month 
!  3 LBDAT Day of month
!  4 LBHR Hour 
!  5 LBMIN Minute 
! 25 Processing code.  Value of 0 is instantaneous.
! 40 NADDR Start address in DATA.
! 42 ITEM_CODE Stash section number and Item number. Held in the form : 
!    Section Number*1000 + Item Number
! 
INTEGER, PARAMETER :: lo_year   = 1
INTEGER, PARAMETER :: lo_month  = 2
INTEGER, PARAMETER :: lo_day    = 3
INTEGER, PARAMETER :: lo_hour   = 4
INTEGER, PARAMETER :: lo_min    = 5
INTEGER, PARAMETER :: lo_lbproc = 25
INTEGER, PARAMETER :: lo_start  = 40
INTEGER, PARAMETER :: lo_stash  = 42
!
! maximum number of free tracers and UKCA tracers
INTEGER, PARAMETER :: max_tracers=150
INTEGER, PARAMETER :: max_tr_ukca=150
! Section number for free/UKCA tracers
INTEGER, PARAMETER :: sec_freetr =33 
INTEGER, PARAMETER :: sec_ukcatr =34 

END MODULE makebc_constants_mod
