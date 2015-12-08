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
! Purpose: Emissions scheme
!
!     Purpose
!     -------
!     Dummy routine which the user must supply a new version if they
!     need to use an emissions scheme. This subroutine is called
!     by ASAD to set the emission rates.
!     Not currently used in UKCA as emissions are added elsewhere.
!
!     Called from asad_cdrive
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
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
        SUBROUTINE ASAD_EMISSN

        USE ASAD_MOD
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE Control_Max_Sizes
        IMPLICIT NONE

!     include any commons here for getting any info
!     from the calling model
!
!    e.g. common /model/ rates(jpnl)
!
!     ---------------------------------------------------------------

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

        IF (lhook) CALL dr_hook('ASAD_EMISSN',zhook_in,zhook_handle)
        RETURN
        IF (lhook) CALL dr_hook('ASAD_EMISSN',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_EMISSN
