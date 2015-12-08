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
! Purpose: This routine is intended for use with stiff integrators which
!     require the evaluation of df/dt at any point in time.
!
!     This routine is intended to be passed as an argument to the stiff
!     integrators such as the NAG and SVODE drivers.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_DIFFUN
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Interface
!     ---------
!        tz     - on entry, specifies the time (unchanged).
!        f1     - on entry, contains species at time TZ (unchanged).
!        dfdt   - on exit, must contain time derivative of the gridpt.
!
!     As all the chemistry arrays are in common, we have to copy the
!     input to the common array and copy the tendency to the output
!     argument at the end.
!
!     Method
!     ------
!     IMPORTANT!!! This subroutine assumes that only a single gridpt
!     is being worked on by the integrator calling this subroutine.
!
!     The first step is to copy the passed values of the
!     species back into the ASAD common blocks so that the production
!     and loss can be computed. Since the integrator is only integrating
!     the variables, f, passed between the model and ASAD we only copy t
!     the species that will change during the timestep. Note that we onl
!     copy to the first element in the species array, y.
!
!     When no families are in use, we can use the index array nlf since
!     this stores a list of all the species of type TR and nf = jpctr.
!
!     The next step is to compute the rates of change. The routine diffu
!     is used but we tell it only to compute the first gridpt. On exit
!     from diffun, the fdot array will have been assigned.
!
!     Externals
!     ---------
!       diffun    - to compute rates of change.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_YCN(tz,f1,dfdt)

        USE ASAD_MOD,           ONLY: nf, nlf, y, fdot
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE ereport_mod, ONLY : ereport
        USE Control_Max_Sizes
        USE domain_params
        IMPLICIT NONE


        REAL, INTENT(IN)  :: tz
        REAL, INTENT(IN)  :: f1(jpctr)

        REAL, INTENT(OUT) :: dfdt(jpctr)

!       Local variables

        INTEGER :: j    ! Loop variable
        INTEGER :: js

        LOGICAL, SAVE :: gfirst = .TRUE.

        CHARACTER(LEN=72) :: cmessage  ! Error message

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1. Initialise species array; ONLY COPY TO FIRST ELEMENT!
!          ---------- ------- ------ ---- ---- -- ----- --------

        IF (lhook) CALL dr_hook('ASAD_YCN',zhook_in,zhook_handle)
        IF ( gfirst ) THEN
          IF ( method  >=  10 .and. jpctr  /=  nf ) then
            WRITE (6,*) '** INTERNAL ASAD ERROR: jpctr  /=  nf in ycn',&
            ' There should not be any families in use with the stiff', &
            ' integrators.'
            cmessage = 'jpctr  /=  nf in ycn'

            CALL EREPORT('ASAD_YCN',jpctr,cmessage)
         ENDIF
         gfirst = .false.
        ENDIF

!       n.b. nf = jpctr when no families are in use

        DO j = 1, nf
          js  = nlf(j)
          y(1,js) = f1(j)
        ENDDO

!       2.   Compute rates of change.
!            ------- ----- -- -------

! DEPENDS ON: asad_diffun
        CALL ASAD_DIFFUN( 1 )

        DO j = 1, jpctr
          dfdt(j) = fdot(1,j)
        ENDDO

        IF (lhook) CALL dr_hook('ASAD_YCN',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_YCN
!--------------------------------------------------------------------
