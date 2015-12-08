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
! Purpose: Calculates the total number density in gridboxes
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Method
!     ------
!     The total number density at a point is given by: p=nkt
!     p-pressure , n-number density, k-boltzmann's constant,
!     t-temperature.
!
!     local variables
!     ---------------
!     zboltz     boltzmann's constant
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_TOTNUD(n_points)

        USE ASAD_MOD,             ONLY: tnd, p, t, pmintnd, pmin
        USE UKCA_CONSTANTS,       ONLY: zboltz
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        IMPLICIT NONE


        INTEGER, INTENT(IN) :: n_points

!       Local variables

        INTEGER :: jl

        REAL :: zb

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


        IF (lhook) CALL dr_hook('ASAD_TOTNUD',zhook_in,zhook_handle)
        zb = zboltz*1.0e6

!       1. Total number density (1e6 converts numbers to /cm**3).
!          ----- ------ ------- ---- -------- ------- -- --------

        DO jl = 1, n_points
          tnd(jl)     = p(jl) / ( zb * t(jl) )
          pmintnd(jl) = pmin * tnd(jl)
        ENDDO

        IF (lhook) CALL dr_hook('ASAD_TOTNUD',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_TOTNUD
