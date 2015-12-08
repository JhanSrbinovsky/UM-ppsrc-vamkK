! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Using the logicals read in from  the name list
! set up aerosol LBC logicals and increments
! 
!
! Subroutine Interface:
SUBROUTINE makebc_lbc_logic(num_optional_lbcs_out,ndustbin_in, ndustbin_out)

  USE PrintStatus_mod
  USE domain_params
  USE dust_parameters_mod, ONLY:                  l_dust,                   &
       l_dust_div1,            l_dust_div2,       l_dust_div3,              &
       l_dust_div4,            l_dust_div5,       l_dust_div6
  USE um_input_control_mod,  ONLY:                                          &
       l_dust_div1_lbc,        l_dust_div1_lbc_out,                         &
       l_dust_div2_lbc,        l_dust_div2_lbc_out,                         &
       l_dust_div3_lbc,        l_dust_div3_lbc_out,                         &
       l_dust_div4_lbc,        l_dust_div4_lbc_out,                         &
       l_dust_div5_lbc,        l_dust_div5_lbc_out,                         &
       l_dust_div6_lbc,        l_dust_div6_lbc_out,                         &
       l_so2,                  l_so2_lbc,                                   &
       l_so2_lbc_out,          l_dms,             l_dms_lbc,                &
       l_dms_lbc_out,          l_so4_aitken,      l_so4_aitken_lbc,         &
       l_so4_aitken_lbc_out,   l_so4_accu,        l_so4_accu_lbc,           &
       l_so4_accu_lbc_out,     l_so4_diss,        l_so4_diss_lbc,           &
       l_so4_diss_lbc_out,     l_nh3,             l_nh3_lbc,                &
       l_nh3_lbc_out,          l_soot_new,        l_soot_new_lbc,           &
       l_soot_new_lbc_out,     l_soot_agd,        l_soot_agd_lbc,           &
       l_soot_agd_lbc_out,     l_soot_cld,        l_soot_cld_lbc,           &
       l_soot_cld_lbc_out,     l_bmass_new,       l_bmass_new_lbc,          &
       l_bmass_new_lbc_out,    l_bmass_agd,       l_bmass_agd_lbc,          &
       l_bmass_agd_lbc_out,    l_bmass_cld,       l_bmass_cld_lbc,          &
       l_bmass_cld_lbc_out,    l_ocff_new,        l_ocff_new_lbc,           &
       l_ocff_new_lbc_out,     l_ocff_agd,        l_ocff_agd_lbc,           &
       l_ocff_agd_lbc_out,     l_ocff_cld,        l_ocff_cld_lbc,           &
       l_ocff_cld_lbc_out,     l_nitr_acc,        l_nitr_acc_lbc,           &
       l_nitr_acc_lbc_out,     l_nitr_diss,       l_nitr_diss_lbc,          &
       l_nitr_diss_lbc_out,    l_soot,            l_biomass,                &
       l_ocff,                 l_nitrate

  IMPLICIT NONE

!
! Description: set up aerosol LBC logicals and increment
! num_optional_lbcs_out for each output lbc field
!
! Method: Unlike other optional LBCs in the UM, aerosol related LBCs
! do not have to be output from the model which means there is more 
! than 1 logical to be set. In addition for some LBCs the aerosols
! have several size bins to be taken care of. All logicals 
! are input and output via um_input_control_mod
!
! If none of the aerosol namelist items are true this routine does
! nothing
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.

! *********************************************************************
! Subroutine arguments

INTEGER, INTENT(INOUT) ::   num_optional_lbcs_out ! Number of optional lbcs required

INTEGER :: ndustbin_in, ndustbin_out  ! Number of dust bins in/out

!
! End of header
! *********************************************************************
!
! *********************************************************************
! 1. Dust
! *********************************************************************
! If dust is on 6 or 2 LBCs are active - one for each of the bins
! We have 2 logicals to set - one to say this tracer is active
! and one to indicate we wish to output it
IF (l_dust) THEN

! Always have at least two dust bins if dust is switched on
  l_dust_div1=.true.
  l_dust_div2=.true.
  IF (ndustbin_in == 2) THEN
    l_dust_div3=.false.
    l_dust_div4=.false.
    l_dust_div5=.false.
    l_dust_div6=.false.
  ELSE 
    l_dust_div3=.true.
    l_dust_div4=.true.
    l_dust_div5=.true.
    l_dust_div6=.true.
  END IF  

 
  l_dust_div1_lbc_out=.true.
  l_dust_div2_lbc_out=.true.

  IF (ndustbin_out == 2) THEN
    l_dust_div3_lbc_out=.false.
    l_dust_div4_lbc_out=.false.
    l_dust_div5_lbc_out=.false.
    l_dust_div6_lbc_out=.false.
    num_optional_lbcs_out=num_optional_lbcs_out+2
  ELSE
    l_dust_div3_lbc_out=.true.
    l_dust_div4_lbc_out=.true.
    l_dust_div5_lbc_out=.true.
    l_dust_div6_lbc_out=.true.
    num_optional_lbcs_out=num_optional_lbcs_out+6
  END IF
  
ELSE ! No dust, no LBCs
  l_dust_div1_lbc_out=.false. 
  l_dust_div2_lbc_out=.false.
  l_dust_div3_lbc_out=.false.
  l_dust_div4_lbc_out=.false.
  l_dust_div5_lbc_out=.false.
  l_dust_div6_lbc_out=.false.
  ! Override namelist settings as l_dust isn't available later on
  ndustbin_out = 0
  ndustbin_in = 0
END IF 

! *********************************************************************
! 2. SO2
! *********************************************************************
! If SO2 is on we have 1 extra LBC 
! 
IF (l_so2) THEN
  l_so2_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+1
END IF

! *********************************************************************
! 3. Sulfate aerosol
! *********************************************************************
! There is no single logical to control the sulfate
! aerosol scheme. Instead we set each logical indpendently
! Increment num_optional_lbcs_out by one for each field

IF (l_so4_aitken) THEN
  l_so4_aitken_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+1
END IF

IF (l_so4_accu) THEN
  l_so4_accu_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+1
END IF

IF (l_so4_diss) THEN
  l_so4_diss_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+1
END IF

! *********************************************************************
! 4. DMS
! *********************************************************************
! If DMS is on we have 1 extra LBC 
!
IF (l_dms) THEN
  l_dms_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+1
END IF

! *********************************************************************
! 5. Ammonia
! *********************************************************************
! If Ammonia is on we have 1 extra LBC 
!
IF (l_nh3) THEN
  l_nh3_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+1
END IF

! *********************************************************************
! 6. Soot aerosol
! *********************************************************************
! If soot is on there are 3 LBCs to be set one for each of the 3 modes
!
IF (l_soot) THEN
  l_soot_new=.true.
  l_soot_agd=.true.
  l_soot_cld=.true.
  l_soot_new_lbc_out=.true.
  l_soot_agd_lbc_out=.true.
  l_soot_cld_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+3
END IF

! *********************************************************************
! 7. Biomass Burning aerosol
! *********************************************************************
! If BB is on there are 3 LBCs to be set one for each of the 3 modes
!
IF (l_biomass) THEN
  l_bmass_new=.true.
  l_bmass_agd=.true.
  l_bmass_cld=.true.
  l_bmass_new_lbc_out=.true.
  l_bmass_agd_lbc_out=.true.
  l_bmass_cld_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+3
END IF

! *********************************************************************
! 8. Organic Carbon Fossil Fuel aerosol
! *********************************************************************
! If OCFF is on there are 3 LBCs to be set one for each of the 3 modes
! 
IF (l_ocff) THEN
  l_ocff_new=.true.
  l_ocff_agd=.true.
  l_ocff_cld=.true.
  l_ocff_new_lbc_out=.true.
  l_ocff_agd_lbc_out=.true.
  l_ocff_cld_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+3
END IF

! *********************************************************************
! 9. Nitrate aerosol
! *********************************************************************
! If soot is on there are 2 LBCs to be set one for each of the 2 modes
IF (l_nitrate) THEN
  l_nitr_acc=.true. 
  l_nitr_diss=.true.
  l_nitr_acc_lbc_out=.true. 
  l_nitr_diss_lbc_out=.true.
  num_optional_lbcs_out=num_optional_lbcs_out+2
END IF

IF (printstatus >= prstatus_diag) THEN
  WRITE(6,*) 'num_optional_lbcs_out=',num_optional_lbcs_out
END IF

END SUBROUTINE makebc_lbc_logic
