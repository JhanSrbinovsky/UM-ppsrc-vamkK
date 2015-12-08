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
! Purpose: Calculates trimolecular rate coefficients
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
!     See the IUPAC reference material on their website for details on
!     calculation of termolecular rates.
!     http://www.iupac-kinetic.ch.cam.ac.uk/
!
!     Local variables
!     ---------------
!     zo         Low pressure limit to rate*density
!     zi         High pressure limit to rate
!     zr         Ratio of zo/zi
!     iho2       Reaction index for HO2+HO2+M
!     ih2o       Array index for advected tracer H2O
!     in2o5      Reaction index for N2O5+M  
!     ino2no3    Reaction index for NO2+NO3+M
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_TRIMOL(n_points)

        USE ASAD_MOD,        ONLY: rk, at, ntrkx, spt, t300, t, peps,   &
                                   tnd, f, wp, advt
        USE ukca_option_mod, ONLY: jpctr, jptk
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE Control_Max_Sizes
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n_points

!       Local variables

        INTEGER, SAVE :: ih2o              ! Index for h2o in tracer array
        INTEGER       :: iho2              ! Index for ho2+ho2 in rk array
        INTEGER       :: in2o5             ! Reaction index for N2O5+M  
        INTEGER       :: ino2no3           ! Reaction index for NO2+NO3+M
        INTEGER       :: j                 ! Loop variable
        INTEGER       :: jl                ! Loop variable
        INTEGER       :: jtr               ! Loop variable
        INTEGER       :: jr                ! Index

        REAL :: zo                         ! k_0
        REAL :: zi                         ! k_infinity
        REAL :: zfc                        ! F_c
        REAL :: zr                         ! k_0/k_infinity
        REAL, ALLOCATABLE, SAVE :: nf(:)   ! component of broadening factor

        LOGICAL, SAVE :: first = .true.

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1.  Calculate trimolecular rate coefficients
!           --------- ------------ ---- ------------

        IF (lhook) CALL dr_hook('ASAD_TRIMOL',zhook_in,zhook_handle)
        iho2  = 0
        in2o5   = 0  
        ino2no3 = 0

        IF (first) THEN
! Calculate the N factor to be used in broadening factor
          ALLOCATE(nf(jptk+1))
          nf(:) = 1.0
          WHERE (at(:,1) > 1E-3) nf(:) = 0.75 - 1.27*ALOG10(at(:,1))
          
! Check if H2O is an advected tracer
          ih2o = 0
          DO jtr = 1, jpctr
            IF ( advt(jtr)  ==  'H2O       ' ) ih2o = jtr
          END DO
          first = .false.
        END IF

        DO j = 1, jptk
          jr = ntrkx(j)

          IF ( spt(j,1) == 'HO2    '.AND. spt(j,2) == 'HO2    ' )       &
               iho2 = jr

!         N2O5 + M  
!         reset to zero each step  

          in2o5 = 0  
          IF (spt(j,1) == 'N2O5      ') in2o5 = jr  
  
!         NO2 + NO3 + M  
!         reset to zero each step  

          ino2no3 = 0  
          IF ((spt(j,1) == 'NO2       ' .AND. spt(j,2) == 'NO3       ') &   
               .OR.(spt(j,2) == 'NO2       ' .AND.                      &  
               spt(j,1) == 'NO3       ')) ino2no3 = jr  

          DO jl = 1, n_points
            zo = at(j,2) * t300(jl)**at(j,3) *                          &
                           EXP( -at(j,4)/t(jl) ) * tnd(jl)
            zi = at(j,5) * t300(jl)**at(j,6) * EXP( -at(j,7)/t(jl) )
            IF ( zo < peps ) THEN
              rk(jl,jr) = zi
            ELSE IF ( zi < peps ) THEN
              rk(jl,jr) = zo
            ELSE
              IF( at(j,1) <= 1.0 ) THEN
                zfc = at(j,1)
              ELSE IF ( in2o5 /= 0 .OR. ino2no3 /= 0 ) THEN ! dependent 
                                                            ! reactions  
                zfc = 2.5*EXP(-1950.0/t(jl))+0.9*EXP(-t(jl)/at(j,1))
              ELSE              ! temperature dependent Fc
                zfc = EXP( -t(jl)/at(j,1) )
                nf(j) = 0.75 - 1.27*ALOG10(zfc)
              END IF
              zr = zo / zi
              rk(jl,jr) = (zo/(1.0+zr)) *                               &
                           zfc**(1.0/(1.0 + (ALOG10(zr)/nf(j))**2))
            END IF
          END DO
        END DO              ! end of loop over jptk

!       2. Dependent reactions.
!          --------- ----------

! HO2 + HO2 [+ M]
        IF (ih2o /= 0 .AND. iho2 /= 0 ) THEN
!         h2o is an advected tracer
          DO jl = 1, n_points
            rk(jl,iho2) = rk(jl,iho2) *                                &
            ( 1.0 + 1.4E-21*f(jl,ih2o)*EXP(2200./t(jl)) )
          END DO
        ELSE IF (ih2o == 0 .AND. iho2 /= 0) THEN
!         use modelled water concentration
          DO jl = 1, n_points
            rk(jl,iho2) = rk(jl,iho2) *                                &
            ( 1.0 + 1.4E-21*wp(jl)*tnd(jl)*EXP(2200./t(jl)) )
          END DO
        END IF

        IF (lhook) CALL dr_hook('ASAD_TRIMOL',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_TRIMOL
