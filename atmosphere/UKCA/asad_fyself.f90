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
! Purpose: Calculates self-reacting terms for individual species.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_FTOY
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Method
!     ------
!     The reaction table is scanned, and the self-reacting terms are
!     calculated from the appropriate rate coefficients.
!
!     Local variables
!     ---------------
!     isp            Index of self reacting species.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_FYSELF(n_points)

        USE ukca_option_mod, ONLY: jpnr
        USE ASAD_MOD,        ONLY: qa, rk, nstst, nlstst, nspi
        USE parkind1,        ONLY: jprb, jpim
        USE yomhook,         ONLY: lhook, dr_hook
        USE Control_Max_Sizes
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n_points    ! No of spatial points

!       Local variables

        INTEGER :: j                       ! Loop variable
        INTEGER :: jl                      ! Loop variable
        INTEGER :: jr                      ! Loop variable
        INTEGER :: js                      ! Index
        INTEGER :: isp                     ! Index

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1.  Initialisation
!           --------------

        IF (lhook) CALL dr_hook('ASAD_FYSELF',zhook_in,zhook_handle)
        DO j = 1, nstst
          js = nlstst(j)
          DO jl = 1, n_points
            qa(jl,js) = 0.0
          ENDDO
        ENDDO

!       2.  Calculate self-reacting terms
!           --------- ------------- -----

        DO jr = 1, jpnr
          isp = nspi(jr,1)
          IF ( isp == nspi(jr,2) ) THEN
            DO jl = 1, n_points
              qa(jl,isp) = qa(jl,isp) + 2.0 * rk(jl,jr)
            ENDDO
          ENDIF
        ENDDO

        IF (lhook) CALL dr_hook('ASAD_FYSELF',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_FYSELF
