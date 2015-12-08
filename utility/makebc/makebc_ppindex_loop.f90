! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Sets up ppindex and pplengths arrays
!
!+ SUBROUTINE makebc_ppindex_loop

SUBROUTINE makebc_ppindex_loop(max_progs,len_ppindex,len1_lookup,  &
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

! Type for time
USE makebc_time_mod, ONLY: &
  time
USE makebc_constants_mod
USE Submodel_Mod
  
IMPLICIT NONE
!
! *********************************************************************
! Description: sets up ppindex and pplengths arrays from values 
! in the lookup table and sets data_present true if any data found
!
! Method: Loop over all lookup tables 
! Set ppindex to the starting location of field which is required, 
! and also set pplengths to the number of times each 2D field is present 
! (i.e. the number of levels) and set the jpointer to a temporary
! value to indicate that the field has been found.
! The correct value of the jpointers is set in makebc_setd1point
! using the value from the d1pointers array
!
! Do the same for free tracers but set up ppindex_tracer 
! and pplengths_tracer
! For UKCA tracers set up ppindex_tr_ukca and pplengths_tr_ukca
! jtracer and jtr_ukca are not set here. 
! Also set lbcdiag_no for tracer items here
!
! If any data is found for the correct time the logical data_present
! will be set to .true.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.

! *********************************************************************

! Declarations
!
!   Scalar arguments with intent(in) :

! max no of prognostics
INTEGER, INTENT(IN) :: max_progs
INTEGER, INTENT(IN) :: len_ppindex   ! Dimension of pp_index
INTEGER, INTENT(IN) :: len1_lookup   ! First dimension of Lookup table
INTEGER, INTENT(IN) :: len2_lookup   ! Second dimension of Lookup table
INTEGER, INTENT(IN) :: no_tracers    ! The number of active tracers
INTEGER, INTENT(IN) :: no_tr_ukca    ! The number of active tracers
!
!
! Array arguments with intent(in) :
! Lookup table from dump being read
INTEGER, INTENT(IN) :: lookup(len1_lookup,len2_lookup)
!
INTEGER, INTENT(IN) :: tracers_active(max_tracers)
INTEGER, INTENT(IN) :: tr_ukca_active(max_tr_ukca)
!  
! The input variable target_time is used to check the 
! date/time of fields in the lookup tables
TYPE(TIME),INTENT(IN) :: target_time
!
!   Scalar arguments with intent(inout) :
!   Pointers 
INTEGER, INTENT(INOUT) :: jorog      !pointer to orography
INTEGER, INTENT(INOUT) :: ju         !pointer to u wind component
INTEGER, INTENT(INOUT) :: jv         !pointer to v wind component
INTEGER, INTENT(INOUT) :: jw         !pointer to w wind component (vert.)
INTEGER, INTENT(INOUT) :: jrho       !pointer to density
INTEGER, INTENT(INOUT) :: jtheta     !pointer to theta
INTEGER, INTENT(INOUT) :: jq         !pointer to specific humidity
INTEGER, INTENT(INOUT) :: jqcl       !pointer to qcl liquid water
INTEGER, INTENT(INOUT) :: jqcf       !pointer to qcf frozen water
INTEGER, INTENT(INOUT) :: jexner     !pointer to exner pressure
INTEGER, INTENT(INOUT) :: ju_adv     !pointer to u advection
INTEGER, INTENT(INOUT) :: jv_adv     !pointer to v advection
INTEGER, INTENT(INOUT) :: jw_adv     !pointer to w advection
!   pc2 and murk pointers
INTEGER, INTENT(INOUT) :: jqcf2      !pointer to qcf2 - cloud ice (cry)
INTEGER, INTENT(INOUT) :: jqrain     !pointer to qrain - rain
INTEGER, INTENT(INOUT) :: jqgraup    !pointer to qgraup - graupel horiz.
INTEGER, INTENT(INOUT) :: jmurk      !pointer to murk
INTEGER, INTENT(INOUT) :: jcf_bulk   !pointer to cf_bulk
INTEGER, INTENT(INOUT) :: jcf_liquid !pointer to cf_liquid
INTEGER, INTENT(INOUT) :: jcf_frozen !pointer to cf_frozen
!   Aerosol pointers 
INTEGER, INTENT(INOUT) :: jdust_div1      ! dust mmr, division 1
INTEGER, INTENT(INOUT) :: jdust_div2      ! dust mmr, division 2
INTEGER, INTENT(INOUT) :: jdust_div3      ! dust mmr, division 3
INTEGER, INTENT(INOUT) :: jdust_div4      ! dust mmr, division 4
INTEGER, INTENT(INOUT) :: jdust_div5      ! dust mmr, division 5
INTEGER, INTENT(INOUT) :: jdust_div6      ! dust mmr, division 6
INTEGER, INTENT(INOUT) :: jso2            ! sulphur dioxide gas
INTEGER, INTENT(INOUT) :: jso4_aitken     ! Aitken mode sulphate aer
INTEGER, INTENT(INOUT) :: jso4_accu       ! accumulation mode sulpha
INTEGER, INTENT(INOUT) :: jso4_diss       ! dissloved  sulphate aero
INTEGER, INTENT(INOUT) :: jdms            ! dimethyl sulphide gas
INTEGER, INTENT(INOUT) :: jnh3            ! ammonia gas mmr
INTEGER, INTENT(INOUT) :: jsoot_new       ! fresh soot mmr
INTEGER, INTENT(INOUT) :: jsoot_agd       ! aged soot mmr
INTEGER, INTENT(INOUT) :: jsoot_cld       ! soot in cloud mmr
INTEGER, INTENT(INOUT) :: jbmass_new      ! fresh biomass mmr
INTEGER, INTENT(INOUT) :: jbmass_agd      ! aged biomass mmr
INTEGER, INTENT(INOUT) :: jbmass_cld      ! cloud biomass mmr
INTEGER, INTENT(INOUT) :: jocff_new       ! fresh OCFF mmr
INTEGER, INTENT(INOUT) :: jocff_agd       ! aged OCFF mmr
INTEGER, INTENT(INOUT) :: jocff_cld       ! OCFF in cloud mmr
INTEGER, INTENT(INOUT) :: jnitr_acc       ! Accumulation nitrate aerosol
INTEGER, INTENT(INOUT) :: jnitr_diss      ! Dissolved nitrate aerosol
INTEGER, INTENT(INOUT) :: lbcdiag_no_tracer   
INTEGER, INTENT(INOUT) :: lbcdiag_no_tr_ukca  
LOGICAL, INTENT(INOUT) :: data_present
!
!   Array arguments with intent(inout) :
INTEGER, INTENT(INOUT) :: lbcdiag_no(max_progs)

! Index of section 0 fields in fields file i.e. ppindex(2,1) 
! gives the first U field in the field file (stash code for
! U is 2)
INTEGER, INTENT(INOUT) :: ppindex(len_ppindex,n_internal_model)     

! Array pplengths allows setppindex to generate extra values required
! by makebc to use readflds to read the dump with um_readdump
! It is the number of 2D fields for each prognostic in the dump (levels)
! It is named num_levels in get_bc
INTEGER, INTENT(INOUT) :: pplengths(len_ppindex,n_internal_model)   

! pplengths and ppindex for free and UKCA tracers
INTEGER, INTENT(INOUT) :: ppindex_tracer(no_tracers,1)
INTEGER, INTENT(INOUT) :: ppindex_tr_ukca(no_tr_ukca,1)
INTEGER, INTENT(INOUT) :: pplengths_tracer(no_tracers,1)
INTEGER, INTENT(INOUT) :: pplengths_tr_ukca(no_tr_ukca,1)
!
! Local variables
INTEGER :: tracer_no, tracer_no_ukca, stashcode, i
INTEGER :: tracer_lbc_no
INTEGER :: tr_ukca_lbc_no
!
! End of header
! *********************************************************************
!
tracer_lbc_no=0
tr_ukca_lbc_no=0
!

! loop over all items in the lookup headers
DO i=1,len2_lookup

! Test for whether the lookup is for the required time
  IF (lookup(lo_year,i)    == target_time%year    .AND.      &
      lookup(lo_month, i)  == target_time%month .AND.      &
      lookup(lo_day, i)    == target_time%day     .AND.      &
      lookup(lo_hour,i)    == target_time%hour    .AND.      &
      lookup(lo_min, i)    == target_time%min     .AND.      &
      lookup(lo_lbproc,i ) == 0 )          THEN

! If the lookup is for the required time then 
! make sure data_present logical is true   
    data_present=.true.
!
! Find the stash code for this field from the lookup header
    stashcode=lookup(lo_stash,i)
!
! Increment pplengths for the relevant
! field (i.e. count number of fields (levels)
! for each stashcode)
! 
! First deal with all section 0 fields
! The first time any stashcode is found the jpointer
! will equal -1 
! This allow us to find the first place this field appears
! in the input file.Set ppindex equal to this location.
! Also set jpointer equal to lookup(40) 
! to indicate that this stashcode has already been
! found
! The correct value of jpointers are set later on
!
    IF (stashcode == stash_jorog) THEN
      pplengths(stash_jorog,1) = pplengths(stash_jorog,1)+1
      IF (jorog == -1) THEN   
        jorog=lookup(lo_start,i)
        ppindex(stash_jorog,1) = i
      END IF

    ELSE IF (stashcode == stash_ju) THEN
      pplengths(stash_ju,1) = pplengths(stash_ju,1)+1
      IF (ju == -1) THEN
        ju=lookup(lo_start,i)
        ppindex(stash_ju,1) = i
      END IF

    ELSE IF (stashcode == stash_jv) THEN
      pplengths(stash_jv,1) = pplengths(stash_jv,1)+1
      IF (jv == -1) THEN
        jv=lookup(lo_start,i)
        ppindex(stash_jv,1) = i
      END IF

    ELSE IF (stashcode == stash_jw) THEN
      pplengths(stash_jw,1) = pplengths(stash_jw,1)+1
      IF (jw == -1) THEN
        jw=lookup(lo_start,i)
        ppindex(stash_jw,1) = i
      END IF

    ELSE IF (stashcode == stash_jrho) THEN
      pplengths(stash_jrho,1) = pplengths(stash_jrho,1)+1
      IF (jrho == -1) THEN
        jrho=lookup(lo_start,i)
        ppindex(stash_jrho,1) = i
      END IF

    ELSE IF (stashcode == stash_jtheta) THEN
      pplengths(stash_jtheta,1) = pplengths(stash_jtheta,1)+1
      IF (jtheta == -1) THEN
        jtheta=lookup(lo_start,i)
        ppindex(stash_jtheta,1) = i
      END IF

    ELSE IF (stashcode == stash_jq) THEN
      pplengths(stash_jq,1) = pplengths(stash_jq,1)+1
      IF (jq == -1) THEN
        jq=lookup(lo_start,i)
        ppindex(stash_jq,1) = i
      END IF

    ELSE IF (stashcode == stash_jqcl) THEN
      pplengths(stash_jqcl,1) = pplengths(stash_jqcl,1)+1
      if (jqcl == -1) THEN
        jqcl=lookup(lo_start,i)
        ppindex(stash_jqcl,1) = i
      END IF

    ELSE IF (stashcode == stash_jqcf) THEN
      pplengths(stash_jqcf,1) = pplengths(stash_jqcf,1)+1
      IF (jqcf == -1) THEN
        jqcf=lookup(lo_start,i)
        ppindex(stash_jqcf,1) = i
      END IF

    ELSE IF (stashcode == stash_jexner) THEN
      pplengths(stash_jexner,1) = pplengths(stash_jexner,1)+1
      IF (jexner == -1) THEN
        jexner=lookup(lo_start,i)
        ppindex(stash_jexner,1) = i
      END IF

    ELSE IF (stashcode == stash_ju_adv) THEN
      pplengths(stash_ju_adv,1) = pplengths(stash_ju_adv,1)+1
      IF (ju_adv == -1) THEN
        ju_adv=lookup(lo_start,i)
        ppindex(stash_ju_adv,1) = i
      END IF

    ELSE IF (stashcode == stash_jv_adv) THEN
      pplengths(stash_jv_adv,1) = pplengths(stash_jv_adv,1)+1
      IF (jv_adv == -1) THEN
        jv_adv=lookup(lo_start,i)
        ppindex(stash_jv_adv,1) = i
      END IF

    ELSE IF (stashcode == stash_jcf_bulk) THEN
      pplengths(stash_jcf_bulk,1) = pplengths(stash_jcf_bulk,1)+1
      IF (jcf_bulk == -1) THEN
        jcf_bulk=lookup(lo_start,i)
        ppindex(stash_jcf_bulk,1) = i
      END IF

    ELSE IF (stashcode == stash_jcf_liquid) THEN
      pplengths(stash_jcf_liquid,1) = pplengths(stash_jcf_liquid,1)+1
      IF (jcf_liquid == -1) THEN
        jcf_liquid=lookup(lo_start,i)
        ppindex(stash_jcf_liquid,1) = i
      END IF

    ELSE IF (stashcode == stash_jcf_frozen) THEN
      pplengths(stash_jcf_frozen,1) = pplengths(stash_jcf_frozen,1)+1
      IF (jcf_frozen == -1) THEN
        jcf_frozen=lookup(lo_start,i)
        ppindex(stash_jcf_frozen,1) = i
      END IF

    ELSE IF (stashcode == stash_jqcf2) THEN
      pplengths(stash_jqcf2,1) = pplengths(stash_jqcf2,1)+1
      IF (jqcf2 == -1) THEN
        jqcf2=lookup(lo_start,i)
        ppindex(stash_jqcf2,1) = i
      END IF

    ELSE IF (stashcode == stash_jqrain) THEN
      pplengths(stash_jqrain,1) = pplengths(stash_jqrain,1)+1
      IF (jqrain == -1) THEN
        jqrain=lookup(lo_start,i)
        ppindex(stash_jqrain,1) = i
      END IF

    ELSE IF (stashcode == stash_jqgraup) THEN
      pplengths(stash_jqgraup,1) = pplengths(stash_jqgraup,1)+1
      IF (jqgraup == -1) THEN
        jqgraup=lookup(lo_start,i)
        ppindex(stash_jqgraup,1) = i
      END IF

    ELSE IF (stashcode == stash_jmurk) THEN
      pplengths(stash_jmurk,1) = pplengths(stash_jmurk,1)+1
      IF (jmurk == -1) THEN
        jmurk=lookup(lo_start,i)
        ppindex(stash_jmurk,1) = i
      END IF

    ELSE IF (stashcode == stash_jw_adv) THEN
      pplengths(stash_jw_adv,1) = pplengths(stash_jw_adv,1)+1
      IF (jw_adv == -1) THEN
        jw_adv=lookup(lo_start,i)
        ppindex(stash_jw_adv,1) = i
      END IF

    ELSE IF (stashcode == stash_jdust_div1) THEN
      pplengths(stash_jdust_div1,1) = pplengths(stash_jdust_div1,1)+1
      IF (jdust_div1 == -1) THEN
        jdust_div1=lookup(lo_start,i)
        ppindex(stash_jdust_div1,1) = i
      END IF

    ELSE IF (stashcode == stash_jdust_div2) THEN
      pplengths(stash_jdust_div2,1) = pplengths(stash_jdust_div2,1)+1
      IF (jdust_div2 == -1) THEN
        jdust_div2=lookup(lo_start,i)
        ppindex(stash_jdust_div2,1) = i
      END IF

    ELSE IF (stashcode == stash_jdust_div3) THEN
      pplengths(stash_jdust_div3,1) = pplengths(stash_jdust_div3,1)+1
      IF (jdust_div3 == -1) THEN
        jdust_div3=lookup(lo_start,i)
        ppindex(stash_jdust_div3,1) = i
      END IF

    ELSE IF (stashcode == stash_jdust_div4) THEN
      pplengths(stash_jdust_div4,1) = pplengths(stash_jdust_div4,1)+1
      IF (jdust_div4 == -1) THEN
        jdust_div4=lookup(lo_start,i)
        ppindex(stash_jdust_div4,1) = i
      END IF

    ELSE IF (stashcode == stash_jdust_div5) THEN
      pplengths(stash_jdust_div5,1) = pplengths(stash_jdust_div5,1)+1
      IF (jdust_div5 == -1) THEN
        jdust_div5=lookup(lo_start,i)
        ppindex(stash_jdust_div5,1) = i
      END IF

    ELSE IF (stashcode == stash_jdust_div6) THEN
      pplengths(stash_jdust_div6,1) = pplengths(stash_jdust_div6,1)+1
      IF (jdust_div6 == -1) THEN
        jdust_div6=lookup(lo_start,i)
        ppindex(stash_jdust_div6,1) = i
      END IF

    ELSE IF (stashcode == stash_jso2) THEN
      pplengths(stash_jso2,1) = pplengths(stash_jso2,1)+1
      IF (jso2 == -1) THEN
        jso2=lookup(lo_start,i)
        ppindex(stash_jso2,1) = i
      END IF

    ELSE IF (stashcode == stash_jso4_aitken) THEN
      pplengths(stash_jso4_aitken,1) = pplengths(stash_jso4_aitken,1)+1
      IF (jso4_aitken == -1) THEN
        jso4_aitken=lookup(lo_start,i)
        ppindex(stash_jso4_aitken,1) = i
      END IF

    ELSE IF (stashcode == stash_jso4_accu) THEN
      pplengths(stash_jso4_accu,1) = pplengths(stash_jso4_accu,1)+1
      IF (jso4_accu == -1) THEN
        jso4_accu=lookup(lo_start,i)
        ppindex(stash_jso4_accu,1) = i
      END IF

    ELSE IF (stashcode == stash_jso4_diss) THEN
      pplengths(stash_jso4_diss,1) = pplengths(stash_jso4_diss,1)+1
      IF (jso4_diss == -1) THEN
        jso4_diss=lookup(lo_start,i)
        ppindex(stash_jso4_diss,1) = i
      END IF

    ELSE IF (stashcode == stash_jdms) THEN
      pplengths(stash_jdms,1) = pplengths(stash_jdms,1)+1
      IF (jdms == -1) THEN
        jdms=lookup(lo_start,i)
        ppindex(stash_jdms,1) = i
      END IF

    ELSE IF (stashcode == stash_jnh3) THEN
      pplengths(stash_jnh3,1) = pplengths(stash_jnh3,1)+1
      IF (jnh3 == -1) THEN
        jnh3=lookup(lo_start,i)
        ppindex(stash_jnh3,1) = i
      END IF

    ELSE IF (stashcode == stash_jsoot_new) THEN
      pplengths(stash_jsoot_new,1) = pplengths(stash_jsoot_new,1)+1
      IF (jsoot_new == -1) THEN
        jsoot_new=lookup(lo_start,i)
        ppindex(stash_jsoot_new,1) = i
      END IF

    ELSE IF (stashcode == stash_jsoot_agd) THEN
      pplengths(stash_jsoot_agd,1) = pplengths(stash_jsoot_agd,1)+1
      IF (jsoot_agd == -1) THEN
        jsoot_agd=lookup(lo_start,i)
        ppindex(stash_jsoot_agd,1) = i
      END IF

    ELSE IF (stashcode == stash_jsoot_cld) THEN
      pplengths(stash_jsoot_cld,1) = pplengths(stash_jsoot_cld,1)+1
      IF (jsoot_cld == -1) THEN
        jsoot_cld=lookup(lo_start,i)
        ppindex(stash_jsoot_cld,1) = i
      END IF

    ELSE IF (stashcode == stash_jbmass_new) THEN
      pplengths(stash_jbmass_new,1) = pplengths(stash_jbmass_new,1)+1
      IF (jbmass_new == -1) THEN
        jbmass_new=lookup(lo_start,i)
        ppindex(stash_jbmass_new,1) = i
      END IF

    ELSE IF (stashcode == stash_jbmass_agd) THEN
      pplengths(stash_jbmass_agd,1) = pplengths(stash_jbmass_agd,1)+1
      IF (jbmass_agd == -1) THEN
        jbmass_agd=lookup(lo_start,i)
        ppindex(stash_jbmass_agd,1) = i
      END IF

    ELSE IF (stashcode == stash_jbmass_cld) THEN
      pplengths(stash_jbmass_cld,1) = pplengths(stash_jbmass_cld,1)+1
      IF (jbmass_cld == -1) THEN
        jbmass_cld=lookup(lo_start,i)
        ppindex(stash_jbmass_cld,1) = i
      END IF

    ELSE IF (stashcode == stash_jocff_new) THEN
      pplengths(stash_jocff_new,1) = pplengths(stash_jocff_new,1)+1
      IF (jocff_new == -1) THEN
        jocff_new=lookup(lo_start,i)
        ppindex(stash_jocff_new,1) = i
      END IF

    ELSE IF (stashcode == stash_jocff_agd) THEN
      pplengths(stash_jocff_agd,1) = pplengths(stash_jocff_agd,1)+1
      IF (jocff_agd == -1) THEN
        jocff_agd=lookup(lo_start,i)
        ppindex(stash_jocff_agd,1) = i
      END IF

    ELSE IF (stashcode == stash_jocff_cld) THEN
      pplengths(stash_jocff_cld,1) = pplengths(stash_jocff_cld,1)+1
      IF (jocff_cld == -1) THEN
        jocff_cld=lookup(lo_start,i)
        ppindex(stash_jocff_cld,1) = i
      END IF

    ELSE IF (stashcode == stash_jnitr_acc) THEN
      pplengths(stash_jnitr_acc,1) = pplengths(stash_jnitr_acc,1)+1
      IF (jnitr_acc == -1) THEN
        jnitr_acc=lookup(lo_start,i)
        ppindex(stash_jnitr_acc,1) = i
      END IF

    ELSE IF (stashcode == stash_jnitr_diss) THEN
      pplengths(stash_jnitr_diss,1) = pplengths(stash_jnitr_diss,1)+1
      IF (jnitr_diss == -1) THEN
        jnitr_diss=lookup(lo_start,i)
        ppindex(stash_jnitr_diss,1) = i
      END IF

    ! Now deal with free tracers 
    ! First test if this lookup header is a free tracer: 33001 to 33150
    !
    ELSE IF((stashcode >= sec_freetr*1000+1) .AND.              &
       (stashcode <= sec_freetr*1000+max_tracers)) THEN
          !
          ! tracer_no is the item number in S33 for this field 
          !
          tracer_no=stashcode-sec_freetr*1000
          !
          ! check if LBCs for this tracer are required
          IF(tracers_active(tracer_no)==1) THEN
            ! 
            ! Is this tracer different from the previous one?
            ! If so assume it is the first time we have found this tracer
            ! Set ppindex_tracer to this location
            IF(stashcode /= lookup(lo_stash,i-1)) THEN
              !
              ! Increment counter for free tracer LBCs 
              tracer_lbc_no=tracer_lbc_no+1               
              !
              lbcdiag_no(lbcdiag_no_tracer)=tracer_lbc_no 
              
              ! Set ppindex_tracer  
              ppindex_tracer(tracer_lbc_no,1) = i       
              
              lbcdiag_no_tracer=lbcdiag_no_tracer+1
            END IF
            !
            ! Increment pplengths_tracer for this tracer
            pplengths_tracer(tracer_lbc_no,1) =                    &
                      pplengths_tracer(tracer_lbc_no,1)+1
          END IF

    ! Now UKCA tracers 
    ! First test if this lookup header is a UKCA tracer: 34001 to 34150

    ELSE IF((stashcode >= sec_ukcatr*1000+1) .AND.              &
       (stashcode <= sec_ukcatr*1000+max_tr_ukca)) THEN
          !
          ! Get item number in S34 for this lookup from stash code 
          !
          tracer_no_ukca=stashcode-sec_ukcatr*1000
          !
          ! check if LBCs for this tracer are required
          IF(tr_ukca_active(tracer_no_ukca)==1) THEN
            ! 
            ! Is this tracer different from the previous one?
            IF(stashcode /= lookup(lo_stash,i-1)) THEN
              !
              ! Increment counter for UKCA tracer LBCs 
              tr_ukca_lbc_no=tr_ukca_lbc_no+1               
              !
              lbcdiag_no(lbcdiag_no_tr_ukca)=tr_ukca_lbc_no 
              
              ! Set ppindex_tracer  
              ppindex_tr_ukca(tr_ukca_lbc_no,1) = i
              
              lbcdiag_no_tr_ukca=lbcdiag_no_tr_ukca+1
            END IF
            !
            ! Increment pplengths_tracer for this tracer
            pplengths_tr_ukca(tr_ukca_lbc_no,1) =                    &
                      pplengths_tr_ukca(tr_ukca_lbc_no,1)+1
          END IF

    END IF  ! Tests on stashcodes

  END IF ! Test on target_time
END DO ! End of loop over all lookup headers

END SUBROUTINE makebc_ppindex_loop
