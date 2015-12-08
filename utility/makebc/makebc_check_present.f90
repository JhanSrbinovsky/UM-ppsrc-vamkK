! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Interface:
SUBROUTINE makebc_check_present(                                          &
  jorog, ju, jv, jw, jrho, jtheta, jq, jqcl, jqcf,                        &
  jexner, ju_adv, jv_adv, jw_adv, jqcf2, jqrain, jqgraup,                 &
  jmurk, jcf_bulk, jcf_liquid, jcf_frozen,                                &
  jdust_div1, jdust_div2, jdust_div3, jdust_div4,                         &
  jdust_div5, jdust_div6, jso2, jso4_aitken,                              &
  jso4_accu, jso4_diss, jdms, jnh3, jsoot_new,                            &
  jsoot_agd, jsoot_cld, jbmass_new, jbmass_agd,                           &
  jbmass_cld, jocff_new, jocff_agd, jocff_cld,                            &
  jnitr_acc, jnitr_diss, tr_levels, no_tracers, no_tr_ukca,               &
  l_pc2, l_murk, l_mcr_qcf2,l_mcr_qrain, l_mcr_qgraup,                    &
  l_dust, l_so2, l_so4_aitken, l_so4_accu, l_so4_diss, l_dms,             &
  l_nh3, l_soot, l_biomass, l_ocff, l_nitrate, ppindex_tracer,            &
  ppindex_tr_ukca, tracers_active, tr_ukca_active, ndustbin_in,            &
  ndustbin_out )
!
!
! Description: test if any required fields are missing in the input 
!              dump
!
! Method:      If a field has not been found its jpointer will be -1
!              Test if any jpointers for requried fields are -1
!              If any values of -1 for LBCs which are required  
!              this is a fatal error - call Ereport
!
!              For tracers jpointers are not set by makebc_ppindex_loop
!              Instead we test ppindex_tracer and ppindex_tr_ukca
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.

      USE makebc_constants_mod, ONLY: max_tracers,max_tr_ukca
      
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE
! Subroutine arguments
! Scalar arguments with INTENT(IN):
!
! jpointers
INTEGER, INTENT(IN) :: jorog, ju, jv, jw, jrho, jtheta, jq, jqcl,         &
  jqcf, jexner, ju_adv, jv_adv, jw_adv, jqcf2, jqrain, jqgraup,           &
  jmurk, jcf_bulk, jcf_liquid, jcf_frozen,                                &
  jdust_div1, jdust_div2, jdust_div3, jdust_div4,                         &
  jdust_div5, jdust_div6, jso2, jso4_aitken,                              &
  jso4_accu, jso4_diss, jdms, jnh3, jsoot_new,                            &
  jsoot_agd, jsoot_cld, jbmass_new, jbmass_agd,                           &
  jbmass_cld, jocff_new, jocff_agd, jocff_cld,                            &
  jnitr_acc, jnitr_diss                                                   
!
INTEGER, INTENT(IN) :: tr_levels     ! The number of tracer levels
INTEGER, INTENT(IN) :: no_tracers    ! The number of active tracers
INTEGER, INTENT(IN) :: no_tr_ukca    ! The number of active tracers
  
!
! logicals for optional LBCs
LOGICAL, INTENT(IN) :: l_pc2, l_murk, l_mcr_qcf2,                         &
  l_mcr_qrain, l_mcr_qgraup,                                              &
  l_dust,l_so2,l_so4_aitken,l_so4_accu,l_so4_diss,l_dms,                  &
  l_nh3,l_soot,l_biomass,l_ocff,l_nitrate
! 
! Pointer to free and UKCA tracer(s)
INTEGER, INTENT(IN) :: ppindex_tracer(no_tracers,1)
INTEGER, INTENT(IN) :: ppindex_tr_ukca(no_tr_ukca,1)
!
! active tracers
INTEGER, INTENT(IN) :: tracers_active(max_tracers)
INTEGER, INTENT(IN) :: tr_ukca_active(max_tr_ukca)

! Number of dust bins in the input dump (in) and output LBCs (out)
INTEGER, INTENT(IN) :: ndustbin_in, ndustbin_out

!
! Local variables
INTEGER :: errorstatus
CHARACTER(LEN=80) ::  cmessage   !error message
!
! Local parameters 
CHARACTER (LEN=80), PARAMETER :: RoutineName='makebc_check_present'
!
INTEGER  ::  i,k
!
! End of header
! *********************************************************************
! 1. Check if any required LBCs are missing - fatal error 

IF (jorog  == -1 .OR.                                           &
  ju     == -1 .OR.                                             &
  jv     == -1 .OR.                                             &
  jw     == -1 .OR.                                             &
  jrho   == -1 .OR.                                             &
  jtheta == -1 .OR.                                             &
  jq     == -1 .OR.                                             &
  jqcl   == -1 .OR.                                             &
  jqcf   == -1 .OR.                                             &
  jexner == -1 .OR.                                             &
  ju_adv == -1 .OR.                                             &
  jv_adv == -1 .OR.                                             &
  jw_adv == -1                                                  &
!
! *********************************************************************
! 2. Check for pc2, microphysics and murk only if they're required
  .OR. (jqcf2         == -1 .AND. L_mcr_qcf2)                   &
  .OR. (jqrain        == -1 .AND. L_mcr_qrain)                  &
  .OR. (jqgraup       == -1 .AND. L_mcr_qgraup)                 &
  .OR. (jmurk         == -1 .AND. L_murk)                       &
  .OR. (jcf_bulk      == -1 .AND. L_pc2)                        &
  .OR. (jcf_liquid    == -1 .AND. l_pc2)                        &
  .OR. (jcf_frozen    == -1 .AND. l_pc2)                        &
! *********************************************************************
! 3. Check for aerosol related fields if on
  .OR. (jdust_div1    == -1 .AND. l_dust)                       &
  .OR. (jdust_div2    == -1 .AND. l_dust)                       &
  .OR. (jdust_div3    == -1 .AND. ndustbin_in == 6)             &
  .OR. (jdust_div4    == -1 .AND. ndustbin_in == 6)             &
  .OR. (jdust_div5    == -1 .AND. ndustbin_in == 6)             &
  .OR. (jdust_div6    == -1 .AND. ndustbin_in == 6)             &
  .OR. (jso2          == -1 .AND. l_so2)                        &
  .OR. (jso4_aitken   == -1 .AND. l_so4_aitken)                 &
  .OR. (jso4_accu     == -1 .AND. l_so4_accu)                   &
  .OR. (jso4_diss     == -1 .AND. l_so4_diss)                   &
  .OR. (jdms          == -1 .AND. l_dms)                        &
  .OR. (jnh3          == -1 .AND. l_nh3)                        &
  .OR. (jsoot_new     == -1 .AND. l_soot)                       &
  .OR. (jsoot_agd     == -1 .AND. l_soot)                       &
  .OR. (jsoot_cld     == -1 .AND. l_soot)                       &
  .OR. (jbmass_new    == -1 .AND. l_biomass)                    &
  .OR. (jbmass_agd    == -1 .AND. l_biomass)                    &
  .OR. (jbmass_cld    == -1 .AND. l_biomass)                    &
  .OR. (jocff_new     == -1 .AND. l_ocff)                       &
  .OR. (jocff_agd     == -1 .AND. l_ocff)                       &
  .OR. (jocff_cld     == -1 .AND. l_ocff)                       &
  .OR. (jnitr_acc     == -1 .AND. l_nitrate)                    &
  .OR. (jnitr_diss    == -1 .AND. l_nitrate)                    &
) THEN
  errorstatus = 1
  WRITE (6,*) '  Data missing from input file.'
  CMESSAGE  = '  Data missing from input file.'
  WRITE (6,*) ' jorog =    ',jorog
  WRITE (6,*) ' ju=        ',ju
  WRITE (6,*) ' jv=        ',jv
  WRITE (6,*) ' jw=        ',jw
  WRITE (6,*) ' jrho=      ',jrho
  WRITE (6,*) ' jtheta=    ',jtheta
  WRITE (6,*) ' jq=        ',jq
  WRITE (6,*) ' jqcl=      ',jqcl
  WRITE (6,*) ' jqcf=      ',jqcf
  WRITE (6,*) ' jexner=    ',jexner
  WRITE (6,*) ' ju_adv=    ',ju_adv
  WRITE (6,*) ' jv_adv=    ',jv_adv
  WRITE (6,*) ' jw_adv=    ',jw_adv
!
! Note for ease of debugging only write those jpointers which should be set
  IF ( L_mcr_qcf2)   WRITE (6,*) ' jqcf2=     ',jqcf2
  IF ( L_mcr_qrain)  WRITE (6,*) ' jqrain=    ',jqrain
  IF ( L_mcr_qgraup) WRITE (6,*) ' jqgraup=   ',jqgraup
  IF ( L_murk)       WRITE (6,*) ' jmurk=     ',jmurk
  IF ( L_pc2)        WRITE (6,*) ' jcf_bulk=  ',jcf_bulk
  IF ( l_pc2)        WRITE (6,*) ' jcf_liquid=',jcf_liquid
  IF ( l_pc2)        WRITE (6,*) ' jcf_frozen=',jcf_frozen
  IF ( l_dust)       WRITE (6,*) ' jdust_div1=  ',jdust_div1   
  IF ( l_dust)       WRITE (6,*) ' jdust_div2=  ',jdust_div2   
  IF ( ndustbin_in == 6)  WRITE (6,*) ' jdust_div3=  ',jdust_div3   
  IF ( ndustbin_in == 6)  WRITE (6,*) ' jdust_div4=  ',jdust_div4   
  IF ( ndustbin_in == 6)  WRITE (6,*) ' jdust_div5=  ',jdust_div5   
  IF ( ndustbin_in == 6)  WRITE (6,*) ' jdust_div6=  ',jdust_div6   
  IF ( l_so2)        WRITE (6,*) ' jso2=        ',jso2         
  IF ( l_so4_aitken) WRITE (6,*) ' jso4_aitken= ',jso4_aitken  
  IF ( l_so4_accu)   WRITE (6,*) ' jso4_accu=   ',jso4_accu    
  IF ( l_so4_diss)   WRITE (6,*) ' jso4_diss=   ',jso4_diss    
  IF ( l_dms)        WRITE (6,*) ' jdms=        ',jdms         
  IF ( l_nh3)        WRITE (6,*) ' jnh3=        ',jnh3         
  IF ( l_soot)       WRITE (6,*) ' jsoot_new=   ',jsoot_new    
  IF ( l_soot)       WRITE (6,*) ' jsoot_agd=   ',jsoot_agd    
  IF ( l_soot)       WRITE (6,*) ' jsoot_cld=   ',jsoot_cld    
  IF ( l_biomass)    WRITE (6,*) ' jbmass_new=  ',jbmass_new   
  IF ( l_biomass)    WRITE (6,*) ' jbmass_agd=  ',jbmass_agd   
  IF ( l_biomass)    WRITE (6,*) ' jbmass_cld=  ',jbmass_cld   
  IF ( l_ocff)       WRITE (6,*) ' jocff_new=   ',jocff_new    
  IF ( l_ocff)       WRITE (6,*) ' jocff_agd=   ',jocff_agd    
  IF ( l_ocff)       WRITE (6,*) ' jocff_cld=   ',jocff_cld    
  IF ( l_nitrate)    WRITE (6,*) ' jnitr_acc=   ',jnitr_acc    
  IF ( l_nitrate)    WRITE (6,*) ' jnitr_diss=  ',jnitr_diss   


  CALL ereport(routinename,errorstatus,cmessage)
END IF
!
! *********************************************************************
! 4. Check if any required free tracers are missing
!
k=0
DO i=1,max_tracers
  IF (tracers_active(i)==1) THEN        
    k=k+1
    IF(ppindex_tracer(k,1) ==-1) THEN
      errorstatus = 1
      WRITE (6,*) '  Free tracer data missing from input file.'
      CMESSAGE  = '  Free tracer data missing from input file.'
      WRITE (6,*) 'tracers_active=',tracers_active
      WRITE (6,*) 'ppindex_tracer(:,1)=',ppindex_tracer(:,1)

      CALL ereport(routinename,errorstatus,cmessage)
    END IF
  END IF
END DO
!
! *********************************************************************
! 5. Check if any required UKCA tracers are missing
!
k=0
DO i=1,max_tr_ukca
  IF (tr_ukca_active(i)==1) THEN
    k=k+1
    IF (ppindex_tr_ukca(k,1) ==-1) THEN        
      errorstatus = 1
      WRITE (6,*) '  UKCA tracer data missing from input file.'
      CMESSAGE  = '  UKCA tracer data missing from input file.'
      WRITE (6,*) 'tr_ukca_active=',tr_ukca_active
      WRITE (6,*) 'ppindex_tr_ukca(:,1)=',ppindex_tr_ukca(:,1)

      CALL ereport(routinename,errorstatus,cmessage)
    END IF
  END IF
END DO
          
END SUBROUTINE makebc_check_present
