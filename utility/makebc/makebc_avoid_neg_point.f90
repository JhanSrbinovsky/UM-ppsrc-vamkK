! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Avoid negative pointers 
!
! Subroutine Interface:
SUBROUTINE makebc_avoid_neg_point(max_progs,l_pc2, l_murk, l_mcr_qcf2,   &
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

USE makebc_constants_mod
        
USE PrintStatus_mod
IMPLICIT NONE

!
! Description: Avoid negative pointers for optional LBCs which 
!              are unused
!
! Method: Checks which optional LBCs are not in use and set the
!         value of jpointer to 1 and the corrrect element
!         of lbcdiag_no to -1
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.

! Declarations
!
! Global Variables 

! Subroutine arguments
! Scalar arguments with INTENT(IN):
INTEGER, INTENT(IN) :: max_progs
LOGICAL, INTENT(IN) :: l_pc2, l_murk, l_mcr_qcf2,                        &
                 l_mcr_qrain, l_mcr_qgraup,                              &
                 l_dust,l_so2,l_so4_aitken,l_so4_accu,l_so4_diss,l_dms,  &
                 l_nh3,l_soot,l_biomass,l_ocff,l_nitrate
!
INTEGER, INTENT(IN) :: tr_levels     ! The number of tracer levels
INTEGER, INTENT(IN) :: no_tracers    ! The number of active tracers
INTEGER, INTENT(IN) :: no_tr_ukca    ! The number of active tracers
INTEGER, INTENT(IN) :: lbcdiag_no_tracerstart
INTEGER, INTENT(IN) :: lbcdiag_no_trstart_ukca

! Array  arguments with INTENT(INOUT):
! Pointer to free and UKCA tracer(s)
INTEGER, INTENT(INOUT) :: jtracer(tr_levels,no_tracers+1)  
INTEGER, INTENT(INOUT) :: jtr_ukca(tr_levels,no_tr_ukca+1)  
INTEGER, INTENT(INOUT) :: jqcf2,  jqrain,  jqgraup,  jmurk,  jcf_bulk,   &
                 jcf_liquid,  jcf_frozen,  jdust_div1,  jdust_div2,      &
                 jdust_div3,  jdust_div4,  jdust_div5,  jdust_div6,      &
                 jso2,  jso4_aitken,  jso4_accu,  jso4_diss,  jdms,      &
                 jnh3,  jsoot_new,  jsoot_agd,  jsoot_cld,  jbmass_new,  &
                 jbmass_agd,  jbmass_cld,  jocff_new,  jocff_agd,        &
                 jocff_cld,  jnitr_acc,  jnitr_diss
!
! Array  arguments with INTENT(OUT):
INTEGER, INTENT(OUT) :: lbcdiag_no(max_progs)
!
! Local variables
INTEGER :: tracer_loop
INTEGER :: item
!
! End of header
! *********************************************************************
!
! Avoid negative pointers if pc2,microphysics or murk not included

IF (.NOT. l_murk) THEN
  WRITE(6,*)'No murk lbcs'
  jmurk =1
  lbcdiag_no(lbcindex_murk)=-1
END IF
IF (.NOT. l_pc2) THEN
  WRITE(6,*)'No pc2 lbcs'
  jcf_bulk=1
  jcf_liquid=1
  jcf_frozen=1
  lbcdiag_no(lbcindex_cf_bulk  )=-1
  lbcdiag_no(lbcindex_cf_liquid)=-1
  lbcdiag_no(lbcindex_cf_frozen)=-1
END IF
IF (.NOT. l_mcr_qcf2) THEN
  jqcf2  =1
  lbcdiag_no(lbcindex_qcf2)=-1
END IF
IF (.NOT. l_mcr_qrain) THEN
  jqrain =1
  lbcdiag_no(lbcindex_qrain)=-1
END IF
IF (.NOT. l_mcr_qgraup) THEN
  jqgraup=1
  lbcdiag_no(lbcindex_qgraup)=-1
END IF

! Avoid negative pointers for aerosol fields if not included
IF (.NOT. l_dust) THEN
  jdust_div1=1
  jdust_div2=1
  jdust_div3=1
  jdust_div4=1
  jdust_div5=1
  jdust_div6=1
  lbcdiag_no(lbcindex_dust_div1)=-1
  lbcdiag_no(lbcindex_dust_div2)=-1
  lbcdiag_no(lbcindex_dust_div3)=-1
  lbcdiag_no(lbcindex_dust_div4)=-1
  lbcdiag_no(lbcindex_dust_div5)=-1
  lbcdiag_no(lbcindex_dust_div6)=-1
END IF 
IF (.NOT. l_so2) THEN
  jso2 =1
  lbcdiag_no(lbcindex_so2)=-1
END IF
IF (.NOT. l_so4_aitken) THEN
  jso4_aitken=1
  lbcdiag_no(lbcindex_so4_aitken)=-1
END IF
IF (.NOT. l_so4_accu) THEN
  jso4_accu=1
  lbcdiag_no(lbcindex_so4_accu)=-1
END IF
IF (.NOT. l_so4_diss) THEN
  jso4_diss=1
  lbcdiag_no(lbcindex_so4_diss)=-1
END IF
IF (.NOT. l_dms) THEN
  jdms=1
  lbcdiag_no(lbcindex_dms)=-1
END IF
IF (.NOT. l_nh3) THEN
  jnh3=1
  lbcdiag_no(lbcindex_nh3)=-1
END IF
IF (.NOT. l_soot) THEN
  jsoot_new=1
  jsoot_agd=1
  jsoot_cld=1
  lbcdiag_no(lbcindex_soot_new)=-1
  lbcdiag_no(lbcindex_soot_agd)=-1
  lbcdiag_no(lbcindex_soot_cld)=-1
END IF
IF (.NOT. l_biomass) THEN
  jbmass_new=1
  jbmass_agd=1
  jbmass_cld=1
  lbcdiag_no(lbcindex_bmass_new)=-1
  lbcdiag_no(lbcindex_bmass_agd)=-1
  lbcdiag_no(lbcindex_bmass_cld)=-1
END IF
IF (.NOT. l_ocff) THEN
  jocff_new=1
  jocff_agd=1
  jocff_cld=1
  lbcdiag_no(lbcindex_ocff_new)=-1
  lbcdiag_no(lbcindex_ocff_agd)=-1
  lbcdiag_no(lbcindex_ocff_cld)=-1
END IF
IF (.NOT. l_nitrate) THEN
  jnitr_acc=1
  jnitr_diss=1
  lbcdiag_no(lbcindex_nitr_acc)=-1
  lbcdiag_no(lbcindex_nitr_diss)=-1
END IF
! *********************************************************************
!
! If no Tracers, reset JTRACER to prevent negative pointer.
!
IF (no_tracers > 0) then
  DO tracer_loop=lbcdiag_no_tracerstart,                         &
    lbcdiag_no_tracerstart+no_tracers-1
  item=tracer_loop-lbcdiag_no_tracerstart+1
! Avoid negative pointers for all the levels and tracers that are not
! being used
    jtracer(2:tr_levels,item)=1
    IF (jtracer(1,item)==-1) THEN
      IF (printstatus >= prstatus_diag) THEN
        WRITE (6,*) 'avoiding negative pointers for free tracer=',item
        WRITE (6,*) 'settinglbcdiag_no(',tracer_loop,') to -1'
      END IF
      jtracer(1,item)=1
      lbcdiag_no(tracer_loop)=-1
    END IF
  END DO
ELSE
! Avoid negative pointers for whole jtracer array, as tracer lbcs are
! not requested
  jtracer(:,:)=1
  WRITE(6,*) 'No tracer LBCs requested in this run'
END IF
     
! *********************************************************************
!
! If no UKCA Tracers, reset jtr_ukca to prevent negative pointer.
!
IF (no_tr_ukca > 0) then
    DO tracer_loop=lbcdiag_no_trstart_ukca,                         &
                  lbcdiag_no_trstart_ukca+no_tr_ukca-1
! Avoid negative pointers for all the levels and tracers that are not
! being used
      item=tracer_loop-lbcdiag_no_trstart_ukca+1
      jtr_ukca(2:tr_levels,item)=1
      IF (jtr_ukca(1,item)==-1) THEN
        IF (printstatus >= prstatus_diag) THEN
          WRITE (6,*) 'avoiding negative pointers for UKCA tracer=',item
          WRITE (6,*) 'settinglbcdiag_no(',tracer_loop,') to -1'
        END IF
        jtr_ukca(1,item)=1
        lbcdiag_no(tracer_loop)=-1
      END IF
    END DO
ELSE
! Avoid negative pointers for whole jtracer array, as tracer lbcs are
! not requested
  jtr_ukca(:,:)=1
  WRITE(6,*) 'No UKCA tracer LBCs requested in this run'
END IF

RETURN

END SUBROUTINE makebc_avoid_neg_point
