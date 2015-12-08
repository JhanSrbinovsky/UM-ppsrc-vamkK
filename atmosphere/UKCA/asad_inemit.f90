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
! Purpose: Dummy routine which the user must supply a new version if they
!     need to use and initialise their emissions scheme.
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
        SUBROUTINE ASAD_INEMIT

        USE yomhook, ONLY: lhook, dr_hook
        USE parkind1, ONLY: jprb, jpim
        USE ereport_mod, ONLY : ereport
        IMPLICIT NONE

!       include any commons here for getting any info
!       from the calling model

!       e.g. common /model/ rates(jpnl)

        INTEGER                       ::  ICODE

        CHARACTER (Len=80)            ::  CMESSAGE
        CHARACTER (Len=* ), Parameter ::  RoutineName='ASAD_INEMIT'

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


        IF (lhook) CALL dr_hook('ASAD_INEMIT',zhook_in,zhook_handle)
        CMESSAGE = 'Routine should not be callable'
        ICODE = 1

        CALL EREPORT(RoutineName,ICODE,CMESSAGE)

        IF (lhook) CALL dr_hook('ASAD_INEMIT',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_INEMIT
