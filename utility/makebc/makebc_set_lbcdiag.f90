! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Sets up lbcdiag array
!
SUBROUTINE makebc_set_lbcdiag(max_progs,ndustbin_in,lbcdiag_no)
! *********************************************************************
!
USE makebc_constants_mod

USE domain_params
USE dust_parameters_mod, ONLY: l_dust
USE um_input_control_mod,  ONLY:                                        &
     l_so2,                 l_dms,             l_so4_aitken,            &
     l_so4_accu,            l_so4_diss,        l_nh3,                   &
     l_soot,                l_biomass,         l_ocff,                  &
     l_nitrate

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE cloud_inputs_mod, ONLY: l_pc2
USE murk_inputs_mod,  ONLY: l_murk

IMPLICIT NONE
!
! Description: set values of lbcdiag array to stash_xxx parameters
!              for section 0 items. Note that tracers are treated
!              differently (see makebc_ppindex_loop)
!
! Method : Stash codes and index in lbcdiag_no 
!          come from the module makebc_constants_mod
!          Logicals come from the um_input_control module
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC
!
!
INTEGER, INTENT(IN)    :: max_progs
INTEGER, INTENT(IN)    :: ndustbin_in
INTEGER, INTENT(INOUT) :: lbcdiag_no(max_progs)

            
lbcdiag_no(lbcindex_u    )=stash_ju    
lbcdiag_no(lbcindex_v    )=stash_jv    
lbcdiag_no(lbcindex_theta)=stash_jtheta
lbcdiag_no(lbcindex_q    )=stash_jq    
lbcdiag_no(lbcindex_qcf  )=stash_jqcf  
lbcdiag_no(lbcindex_orog )=stash_jorog 

! Set values for murk in lbcdiag_no only
! if required
IF (l_murk) THEN
  lbcdiag_no(lbcindex_murk)=stash_jmurk
END IF

lbcdiag_no(lbcindex_w    )=stash_jw    
lbcdiag_no(lbcindex_rho  )=stash_jrho  
lbcdiag_no(lbcindex_qcl  )=stash_jqcl  
lbcdiag_no(lbcindex_exner)=stash_jexner
lbcdiag_no(lbcindex_u_adv)=stash_ju_adv
lbcdiag_no(lbcindex_v_adv)=stash_jv_adv
lbcdiag_no(lbcindex_w_adv)=stash_jw_adv

! Set values for pc2 in lbcdiag_no only
! if required
IF (l_pc2) THEN
  lbcdiag_no(lbcindex_cf_bulk  )=stash_jcf_bulk  
  lbcdiag_no(lbcindex_cf_liquid)=stash_jcf_liquid
  lbcdiag_no(lbcindex_cf_frozen)=stash_jcf_frozen
END IF

! Set values for microphysics in lbcdiag_no only
! if required
IF (l_mcr_qcf2) THEN
  lbcdiag_no(lbcindex_qcf2)=stash_jqcf2
END IF

IF (l_mcr_qrain) THEN
  lbcdiag_no(lbcindex_qrain)=stash_jqrain
END IF

IF (l_mcr_qgraup) THEN
  lbcdiag_no(lbcindex_qgraup)=stash_jqgraup
END IF

! Set values for aersols in lbcdiag_no only
! if required 
IF (l_dust) THEN
  lbcdiag_no(lbcindex_dust_div1)=stash_jdust_div1 
  lbcdiag_no(lbcindex_dust_div2)=stash_jdust_div2 
  IF (ndustbin_in == 6) THEN
    lbcdiag_no(lbcindex_dust_div3)=stash_jdust_div3 
    lbcdiag_no(lbcindex_dust_div4)=stash_jdust_div4 
    lbcdiag_no(lbcindex_dust_div5)=stash_jdust_div5 
    lbcdiag_no(lbcindex_dust_div6)=stash_jdust_div6 
  END IF
END IF 

IF (l_so2) THEN
  lbcdiag_no(lbcindex_so2)=stash_jso2
END IF 

IF (l_so4_aitken) THEN
  lbcdiag_no(lbcindex_so4_aitken)=stash_jso4_aitken
END IF

IF (l_so4_accu) THEN
  lbcdiag_no(lbcindex_so4_accu)=stash_jso4_accu
END IF

IF (l_so4_diss) THEN
  lbcdiag_no(lbcindex_so4_diss)=stash_jso4_diss
END IF

IF (l_dms) THEN
  lbcdiag_no(lbcindex_dms)=stash_jdms
END IF

IF (l_nh3) THEN
  lbcdiag_no(lbcindex_nh3 )=stash_jnh3
END IF

IF (l_soot) THEN
  lbcdiag_no(lbcindex_soot_new)=stash_jsoot_new
  lbcdiag_no(lbcindex_soot_agd)=stash_jsoot_agd
  lbcdiag_no(lbcindex_soot_cld)=stash_jsoot_cld
END IF 

IF (l_biomass) THEN
  lbcdiag_no(lbcindex_bmass_new)=stash_jbmass_new
  lbcdiag_no(lbcindex_bmass_agd)=stash_jbmass_agd
  lbcdiag_no(lbcindex_bmass_cld)=stash_jbmass_cld
END IF 

IF (l_ocff) THEN
  lbcdiag_no(lbcindex_ocff_new)=stash_jocff_new
  lbcdiag_no(lbcindex_ocff_agd)=stash_jocff_agd
  lbcdiag_no(lbcindex_ocff_cld)=stash_jocff_cld
END IF 

IF (l_nitrate) THEN
  lbcdiag_no(lbcindex_nitr_acc )=stash_jnitr_acc 
  lbcdiag_no(lbcindex_nitr_diss)=stash_jnitr_diss
END IF 

END SUBROUTINE makebc_set_lbcdiag
