! *****************************COPYRIGHT******************************* 
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************* 
!            
! Purpose: To initialize internal values, addresses and
! other information needed for UKCA
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                          
!          Called from SET_ATM_POINTERS
!                             
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!    
!        
! Code description:            
!   Language: FORTRAN 95       
!   This code is written to UMDP3 programming standards.        
!  
! --------------------------------------------------------------------- 
!
      MODULE ukca_init_mod

      IMPLICIT NONE
      
      CONTAINS

      SUBROUTINE ukca_init

      USE ukca_cdnc_mod,         ONLY: ukca_cdnc_init
      USE ukca_option_mod,       ONLY: l_ukca, l_ukca_aie1,        &
                                       l_ukca_aie2, check_run_ukca
      USE ukca_setup_chem_mod,   ONLY: ukca_setup_chem
      USE parkind1,              ONLY: jprb, jpim
      USE yomhook,               ONLY: lhook, dr_hook
      IMPLICIT NONE        
                                                                       
      INTEGER                       :: errcode=0                 ! Error code: ereport
      CHARACTER(LEN=256)            :: cmessage=' '              ! Error message
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_INIT',zhook_in,zhook_handle)

! Check the logical switches
      CALL check_run_ukca()

! Set internal UKCA values from UKCA namelists
      IF (l_ukca) CALL ukca_setup_chem()

! Set item numbers for CDNC etc
      IF (l_ukca_aie1 .OR. l_ukca_aie2) THEN
        CALL ukca_cdnc_init(errcode,cmessage)
      END IF

      IF (lhook) CALL dr_hook('UKCA_INIT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ukca_init
      END MODULE ukca_init_mod
