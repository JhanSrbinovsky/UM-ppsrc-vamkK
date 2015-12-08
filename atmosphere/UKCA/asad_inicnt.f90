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
! Purpose: Species labelled 'CF': ASAD will treat the species as a constant
!     but will call this routine so that the user may set the
!     values differently at each gridpoint for example. Currently
!     only used for setting the water vapour concentration.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_FYINIT
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Interface
!     On entry, the following will be set:
!              species - character name of species to set.
!                        Will be the same as listed in chch.d file
!              klen    - length of array, y.
!
!     On exit, the following must be set:
!              y       - Array of points to set for the species.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INICNT( species, y, klen )

        USE ASAD_MOD,    ONLY: wp, tnd
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE ereport_mod, ONLY : ereport
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: klen      ! No of spatial points

        CHARACTER (LEN=10), INTENT(IN)  :: species  ! Species char strng

        REAL, INTENT(OUT)   :: y(klen)   ! Species concentration

!       Local variables

        INTEGER :: errcode                ! Variable passed to ereport
        INTEGER :: jl

        CHARACTER (LEN=72) :: cmessage

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1.  Copy water into ASAD array.

        IF (lhook) CALL dr_hook('ASAD_INICNT',zhook_in,zhook_handle)
        IF ( species(1:3) /= 'H2O' ) THEN
           errcode=124
           cmessage= 'Expected species H2O but got '//species

           CALL EREPORT('ASAD_INICNT',errcode,cmessage)
        ENDIF

        DO jl = 1, klen
          y(jl) = wp(jl)*tnd(jl)
        ENDDO

        IF (lhook) CALL dr_hook('ASAD_INICNT',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_INICNT
