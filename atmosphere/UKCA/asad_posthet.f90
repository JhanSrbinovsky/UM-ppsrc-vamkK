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
! Purpose: This is a dummy version of a user-supplied subroutine.
!     The role of this routine is to perform any post-processing
!     or house keeping required by the user's heterogeneous
!     chemistry scheme.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
!     Example
!     -------
!     As an example of what this routine might be used for,
!     in our stratospheric chemistry scheme, we have a
!     heterogeneous chemistry scheme that takes into account
!     the amount of water and nitic acid in the solid phase.
!     Therefore, the routine hetero in our chemistry partitions
!     the gaseous HNO3 into a solid and gaseous form. The amount
!     in the solid phase is stored in a common block. This routine
!     is responsible for adding that solid phase amount of
!     HNO3 back into the gaseous phase at the end of the timestep.
!     ie. routine hetero computes the amount of solid HNO3, subtracts
!     that from the array y(), put it in common, where this routine
!     adds it back to y() before ASAD is exited.
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
        SUBROUTINE asad_posthet


        USE yomhook, ONLY: lhook, dr_hook
        USE parkind1, ONLY: jprb, jpim
        USE ereport_mod, ONLY : ereport
        IMPLICIT NONE

        INTEGER                       ::  ICODE
        CHARACTER (Len=80)            ::  CMESSAGE
        CHARACTER (Len=* ), Parameter ::  RoutineName='ASAD_POSTHET'

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


        IF (lhook) CALL dr_hook('ASAD_POSTHET',zhook_in,zhook_handle)
        CMESSAGE = 'Routine should not be callable'
        ICODE = 1

        CALL EREPORT(RoutineName,ICODE,CMESSAGE)

        IF (lhook) CALL dr_hook('ASAD_POSTHET',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_POSTHET
