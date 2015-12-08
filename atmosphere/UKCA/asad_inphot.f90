! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Dummy routine for initialising photolysis
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CINIT
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INPHOT

        USE yomhook, ONLY: lhook, dr_hook
        USE parkind1, ONLY: jprb, jpim
        USE ereport_mod, ONLY : ereport
        IMPLICIT NONE

        INTEGER                       ::  ICODE
        CHARACTER (Len=72)            ::  CMESSAGE
        CHARACTER (Len=*), Parameter  ::  RoutineName='ASAD_INPHOT'

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


        IF (lhook) CALL dr_hook('ASAD_INPHOT',zhook_in,zhook_handle)
        CMESSAGE = 'Dummy routine called'
        ICODE = -1

        CALL EReport(RoutineName,ICODE,CMESSAGE)

        IF (lhook) CALL dr_hook('ASAD_INPHOT',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_INPHOT
