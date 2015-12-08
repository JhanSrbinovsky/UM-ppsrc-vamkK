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
! Purpose: To set up the indexing arrays used in the routine asad_jac.f 
!          to compute the Jacobian matrix elements for the IMPACT time i
!          integration routine
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
!     Interface
!     ---------
!     All variables used and set are in common blocks.
!
!     Method
!     ------
!     The basic procedure is to scan the reactions counting those that
!     will contribute to the Jacobian and storing an entry in a list to
!     that reaction. For species of type 'FM' and 'TR' their
!     contribution will not change during the model run. However, for
!     species that go into and out of a family (type 'FT'), we cannot
!     precompute their contribution, since this is dependent on their
!     lifetime which will vary with time and space. It is more efficient
!     therefore in the Jacobian matrix calculation not to be using
!     species of type 'FT'. Another complication is introduced by the
!     use of chemical families. These can make a positive contribution
!     to the main diagonal of the Jacobian matrix and we take account of
!     these separately (see Carver & Stott, 1997, Annales Geophysicae
!     for more details).
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INIJAC

        USE ASAD_MOD
        USE ukca_option_mod, ONLY: jpctr, jpspec, jpnr
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE ereport_mod, ONLY : ereport
        USE Control_Max_Sizes
        IMPLICIT NONE

!       Local variables

        INTEGER :: ifam(jpmsp)
        INTEGER :: itr(jpmsp)
        INTEGER :: ifneg(2*jpnr,jpctr)
        INTEGER :: ifpos(jppjac,jpctr)
        INTEGER :: ifn(jpctr)
        INTEGER :: ifp(jpctr)
        INTEGER :: ifz(jpctr)
        INTEGER :: ifzer(jpnr,jpctr)
        INTEGER :: ifs(jpctr)
        INTEGER :: ifsss(jpnr,jpctr)
        INTEGER :: ir1                 ! Index
        INTEGER :: ir2                 ! Index
        INTEGER :: ip1                 ! Index
        INTEGER :: ip2                 ! Index
        INTEGER :: ip3                 ! Index
        INTEGER :: is                  ! Index
        INTEGER :: njx                 ! Index
        INTEGER :: j                   ! Loop variable
        INTEGER :: jc                  ! Loop variable
        INTEGER :: je                  ! Loop variable
        INTEGER :: jg                  ! Loop variable
        INTEGER :: jr                  ! Loop variable
        INTEGER :: jtr                 ! Loop variable
        INTEGER :: i
        INTEGER :: irem
        INTEGER :: jpnjcx3

        CHARACTER (LEN=2)  :: itype(jpmsp)
        CHARACTER (LEN=72) :: cmessage     ! Error message

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1. Loop over reactions and set index arrays.
!          ---- ---- --------- --- --- ----- -------

        IF (lhook) CALL dr_hook('ASAD_INIJAC',zhook_in,zhook_handle)

        jpnjcx3 = (jpnr/(3*3))+3*3

        ifn(:) = 0
        ifp(:) = 0
        ifz(:) = 0
        ifs(:) = 0
        ifzer(:,:) = 0
        ifsss(:,:) = 0

!       1.1  Compute contribution from each reaction.

        DO jr = 1, jpnr
          ir1 = nspi(jr,1)
          ir2 = nspi(jr,2)
          ip1 = nspi(jr,3)
          ip2 = nspi(jr,4)
          ip3 = nspi(jr,5)

!         1.1.1  Set indices to determine species type.

          DO j = 1, jpmsp
            ifam(j)  = 0
            itype(j) = '  '
            itr(j)   = 0
          ENDDO
          IF ( ir1 /= 0 ) then
            ifam(1)  = moffam(ir1)
            itype(1) = ctype(ir1)
            itr(1)   = madvtr(ir1)
          ENDIF
          IF ( ir2 /= 0 ) then
            ifam(2)  = moffam(ir2)
            itype(2) = ctype(ir2)
            itr(2)   = madvtr(ir2)
          ENDIF
          IF ( ip1 /= 0 ) then
            ifam(3)  = moffam(ip1)
            itype(3) = ctype(ip1)
            itr(3)   = madvtr(ip1)
          ENDIF
          IF ( ip2 /= 0 ) then
            ifam(4)  = moffam(ip2)
            itype(4) = ctype(ip2)
            itr(4)   = madvtr(ip2)
          ENDIF
          IF ( ip3 /= 0 ) then
            ifam(5)  = moffam(ip3)
            itype(5) = ctype(ip3)
            itr(5)   = madvtr(ip3)
          ENDIF

!         1.2  If first reactant is a family member add the
!              reaction to the negative and positive lists.

          IF ( ifam(1) /= 0 ) then
            is = -nodd(ir1)
            IF ( ifam(2) == ifam(1) .AND. itype(2) == jpfm ) THEN
              is = is - nodd(ir2)
            END IF
            IF ( ifam(3) == ifam(1) .AND. itype(3) == jpfm ) THEN
              is = is + nodd(ip1)
            END IF
            IF ( ifam(4) == ifam(1) .AND. itype(4) == jpfm ) THEN
              is = is + nodd(ip2)
            END IF
            IF ( ifam(5) == ifam(1) .AND. itype(5) == jpfm ) THEN
              is = is + nodd(ip3)
            END IF
            IF ( is < 0 ) THEN
              DO je = 1, abs(is)
                ifn(ifam(1)) = ifn(ifam(1)) + 1
                IF ( ifn(ifam(1)) >  2*jpnr ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifneg dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifneg dimension'

                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifneg(ifn(ifam(1)),ifam(1)) = jr
              END DO
            ELSE IF ( is > 0 ) then
              DO je = 1, is
                ifp(ifam(1)) = ifp(ifam(1)) + 1
                IF ( ifp(ifam(1)) >  jppjac ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifpos dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifpos dimension'

                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifpos(ifp(ifam(1)),ifam(1)) = jr
              END DO
            END IF
          END IF

!         1.3  If second reactant is a family member add to the lists

          IF ( ifam(2) /= 0 ) THEN
            is = -nodd(ir2)
            IF ( ifam(1) == ifam(2) .AND. itype(2) == jpfm ) THEN
              is = is - nodd(ir1)
            END IF
            IF ( ifam(3) == ifam(2) .AND. itype(3) == jpfm ) THEN 
              is = is + nodd(ip1)
            END IF
            IF ( ifam(4) == ifam(2) .AND. itype(4) == jpfm ) THEN
              is = is + nodd(ip2)
            END IF
            IF ( ifam(5) == ifam(2) .AND. itype(5) == jpfm ) THEN
              is = is + nodd(ip3)
            END IF
            IF ( is < 0 ) THEN
              DO je = 1, abs(is)
                ifn(ifam(2)) = ifn(ifam(2)) + 1
                IF ( ifn(ifam(2)) > 2*jpnr ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifneg dimension.'
                  cmessage = 'ARRAY TOO SMALL: Increase ifneg'

                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifneg(ifn(ifam(2)),ifam(2)) = jr
              END DO
            ELSE IF ( is > 0 ) THEN
              DO je = 1, is
                ifp(ifam(2)) = ifp(ifam(2)) + 1
                IF ( ifp(ifam(2)) > jppjac ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifpos dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifpos dimension'

                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifpos(ifp(ifam(2)),ifam(2)) = jr
              END DO
            END IF
          END IF

!         1.4  Calc. contribution from non-family species.

          IF ( itr(1) /= 0 ) then
            is = -1
            IF ( itr(2) == itr(1) .AND. ir2 /= 0 ) THEN
              is = is - 1
            END IF
            IF ( itr(3) == itr(1) .AND. ip1 /= 0 ) THEN
              is = is + 1
            END IF
            IF ( itr(4) == itr(1) .AND. ip2 /= 0 ) THEN 
              is = is + 1
            END IF
            IF ( itr(5) == itr(1) .AND. ip3 /= 0 ) THEN
              is = is + 1
            END IF

            IF ( is < 0 ) then
              DO je = 1, abs(is)
                ifn(itr(1)) = ifn(itr(1)) + 1
                IF ( ifn(itr(1)) > 2*jpnr ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                  ' in inijac.f. Increase ifneg dimension.'
                  cmessage = 'ARRAY TOO SMALL: Increase ifneg'

                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifneg(ifn(itr(1)),itr(1)) = jr
              END DO
            ELSE IF ( is > 0 ) THEN
              DO je = 1, is
                ifp(itr(1)) = ifp(itr(1)) + 1
                IF ( ifp(itr(1)) >  jppjac ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifpos dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifpos dimension'

                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifpos(ifp(itr(1)),itr(1)) = jr
              END DO
            END IF
          END IF

          IF ( itr(2) /= 0 ) then
            is = -1
            IF ( itr(2) == itr(1) .AND. ir1 /= 0 ) THEN 
              is = is - 1
            END IF
            IF ( itr(3) == itr(2) .AND. ip1 /= 0 ) THEN
              is = is + 1
            END IF
            IF ( itr(4) == itr(2) .AND. ip2 /= 0 ) THEN 
              is = is + 1
            END IF
            IF ( itr(5) == itr(2) .AND. ip3 /= 0 ) THEN 
              is = is + 1
            END IF
            
            IF ( is < 0 ) THEN
              DO je = 1, abs(is)
                ifn(itr(2)) = ifn(itr(2)) + 1
                IF ( ifn(itr(2)) >  2*jpnr ) then
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifneg dimension.'
                  cmessage='ARRAY TOO SMALL: Increase igneg dimension'

                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifneg(ifn(itr(2)),itr(2)) = jr
              END DO
            ELSE IF ( is > 0 ) THEN
              DO je = 1, is
                ifp(itr(2)) = ifp(itr(2)) + 1
                IF ( ifp(itr(2)) > jppjac ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                  ' in inijac.f. Increase ifpos dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifpos dimension'

                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifpos(ifp(itr(2)),itr(2)) = jr
              END DO
            END IF
          END IF

! Add bit here for full Jacobian, and save tracer names
        DO j=1,jpmsp
          njcoth(jr,j)=itr(j)
        END DO
        DO j=1,2
          IF(itr(j).ne.0) THEN
            ifz(itr(j)) = ifz(itr(j)) + 1
            ifzer(ifz(itr(j)),itr(j)) = jr
          END IF
        END DO

! Sort out steady state stuff in full Jacobian
        DO j=1,jpmsp
          njcoss(jr,j) = 0
          IF(itype(j) == jpna) THEN
            IF(speci(nspi(jr,j)) == 'O(1D)') THEN
              IF(j <= 2) THEN
                ifs(1) = ifs(1) + 1
                ifsss(ifs(1),1) = jr
              ENDIF
              njcoss(jr,j) = 1
            ELSEIF((o3p_in_ss) .AND.                                    &
               (speci(nspi(jr,j)) == 'O(3P)     ')) THEN
              IF(j <= 2) THEN
                ifs(2) = ifs(2) + 1
                ifsss(ifs(2),2) = jr
              ENDIF
              njcoss(jr,j) = 2
            ELSEIF((n_in_ss) .AND.                                      &
               (speci(nspi(jr,j)) =='N         ')) THEN
              IF(j <= 2) THEN
                ifs(3) = ifs(3) + 1
                ifsss(ifs(3),3) = jr
              ENDIF
              njcoss(jr,j) = 3
            ENDIF
          ENDIF
        ENDDO     ! j=1,jpmsp
      ENDDO      ! jr=1,jpnr



!       2.  Partition lists into 3, 2, and 1 sum loops taking
!           positive and negative contributions into account.
!           -------- -------- ------------- ---- --------

        DO jc = 1, jpctr
          njx = ifn(jc) / 3
          IF ( njx > jpnjcx3 ) THEN
            WRITE (6,*) 'FATAL ASAD ERROR: array njacx3 too ',         &
                     'small in inijac.f. Increase jpnjcx3 to at ',     &
                     'least ', njx
            cmessage = 'ARRAY NJACX3 TOO SMALL: Increase jpnjcx3'

            CALL EREPORT('ASAD_INIJAC',jc,cmessage)
          END IF

! 2.1. Negative 3 term sums.

          i = 1
          njcgrp(jc,3) = njx
          DO jg = 1, njcgrp(jc,3)
            njacx3(1,jg,jc) = ifneg(i,jc)
            njacx3(2,jg,jc) = ifneg(i+1,jc)
            njacx3(3,jg,jc) = ifneg(i+2,jc)
            i               = i + 3
          END DO

! 2.2  Negative 2 and single term sums.

          irem = mod(ifn(jc),3)
          IF ( irem == 2 ) THEN
            njcgrp(jc,2) = 1
            njacx2(1,jc) = ifneg(i,jc)
            njacx2(2,jc) = ifneg(i+1,jc)
            i = i + 2
          ELSE IF ( irem  ==  1 ) THEN
            njcgrp(jc,1) = 1
            njacx1(jc) = ifneg(i,jc)
          END IF
        END DO

! 2.3  Terms which make a positive contribution to the
!      Jacobian. Don't bother to unroll these as there
!      are unlikely to be many of them.

      DO jc = 1, jpctr
        nmpjac(jc) = ifp(jc)
        IF ( ifp(jc) > jppjac ) THEN
          WRITE (6,*) 'FATAL ASAD ERROR: array npjac1 too ',            &
          ' small in injac.f. Increase parameter jppjac to at ',        &
          ' least ',ifp(jc)
          cmessage='ARRAY NPJAC1 TOO SMALL: Increase jppjac'

          CALL EREPORT('ASAD_INIJAC',jc,cmessage)
        END IF
        DO jr = 1, ifp(jc)
          npjac1(jr,jc) = ifpos(jr,jc)
        END DO
      END DO

! 2.4  Terms which make any contribution at all to the
!      Jacobian - for use with full Jacobian only.

      DO jc = 1, jpctr
        nmzjac(jc) = ifz(jc)
        DO jr = 1, ifz(jc)
          nzjac1(jr,jc) = ifzer(jr,jc)
        END DO
      END DO

! 2.5  Terms from steady state species which needed to be
!      considered in the calculation of the full Jacobian.

      DO jc = 1, nstst
        nmsjac(jc) = ifs(jc)
        DO jr = 1, ifs(jc)
          nsjac1(jr,jc) = ifsss(jr,jc)
        END DO
      END DO

      ntro3  = 0
      ntroh  = 0
      ntrho2 = 0
      nspo3  = 0
      nspo1d = 0
      nspo3p = 0
      nspoh  = 0
      nspho2 = 0
      nspno  = 0
      nspn   = 0
      nsph   = 0
      DO jtr = 1, jpctr
        IF (advt (jtr) == 'O3        ') ntro3  = jtr
        IF (advt (jtr) == 'OH        ') ntroh  = jtr
        IF (advt (jtr) == 'HO2       ') ntrho2 = jtr
        IF (advt (jtr) == 'NO        ') ntrno  = jtr
      END DO

      DO jtr=1,jpspec
        IF (speci(jtr) == 'O3        ') nspo3  = jtr
        IF (speci(jtr) == 'O(1D)     ') nspo1d = jtr
        IF (speci(jtr) == 'O(3P)     ') nspo3p = jtr
        IF (speci(jtr) == 'OH        ') nspoh  = jtr
        IF (speci(jtr) == 'HO2       ') nspho2 = jtr
        IF (speci(jtr) == 'NO        ') nspno  = jtr
        IF (speci(jtr) == 'N         ') nspn   = jtr
        IF (speci(jtr) == 'H         ') nsph   = jtr
      END DO

      IF (lhook) CALL dr_hook('ASAD_INIJAC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ASAD_INIJAC
