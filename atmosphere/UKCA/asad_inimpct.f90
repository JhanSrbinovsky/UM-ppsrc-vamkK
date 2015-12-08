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
! Purpose: To be called on the first call to IMPACT. Sets up indexing
!          arrays used by IMPACT routine.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_IMPACT
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Method
!     ------
!     This routine has to compute the lists of species and reactions
!     needed by IMPACT for the correction terms when the Jacobian
!     terms due to photolysis are added in.  It does this by scanning
!     the photolysis ratefile and when it finds a reaction that results
!     in the production of a tracer or family, it stores the relevent
!     information in a list. For a normal tracer (ASAD type 'TR'), we
!     need the reaction no. and the tracer no. For a family, we need
!     the reaction no, tracer no. and species no.
!
!     Arguments
!     ---------
!     Interfaces to ASAD and IMPACT through COMMON variables.
!
!     Local variables
!     ---------------
!     gtr  - Used to indicate this tracer needs correcting
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INIMPCT

        USE ASAD_MOD
        USE ukca_option_mod, ONLY: jpctr, jppj
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE Control_Max_Sizes
        IMPLICIT NONE

!       Local variables

        INTEGER :: j                 ! Loop variable
        INTEGER :: ji                ! Loop variable
        INTEGER :: jii               ! Loop variable
        INTEGER :: jm                ! Loop variable
        INTEGER :: jp                ! Loop variable
        INTEGER :: jt                ! Loop variable
        INTEGER :: ipr
        INTEGER :: ip
        INTEGER :: ientry
        INTEGER :: ntr
        INTEGER :: ientries
        INTEGER :: index
        INTEGER :: inf
        INTEGER :: ipos
        INTEGER :: icount
        INTEGER :: ir
        INTEGER :: ir2
        INTEGER :: irk
        INTEGER :: irk2
        INTEGER :: is
        INTEGER :: is2
        INTEGER :: irtr
        INTEGER :: irtr2
        INTEGER :: inr
        INTEGER :: iprod

        INTEGER :: itr(jpspj-2)
        INTEGER :: ifam(jpspj-2)
        INTEGER :: iwork(0:(jpspj-2)*jppj,2)

        LOGICAL :: gtr(jpctr)

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1.  Loop over photolysis reactions
!           ---- ---- ---------- ---------

        IF (lhook) CALL dr_hook('ASAD_INIMPCT',zhook_in,zhook_handle)
        DO jt = 1, jpctr
          gtr(jt) = .false.
        ENDDO
        DO jt = 0, jpctr
          nltrim(jt,1) = 0
          nltrim(jt,2) = 0
          nltrim(jt,3) = 0
        ENDDO
        DO jt = 1, (jpspj-2)*jppj
          iwork(jt,1) = 0
          iwork(jt,2) = 0
          nlpdv(jt,1) = 0
          nlpdv(jt,2) = 0
        ENDDO
        iwork(0,1) = 0
        iwork(0,2) = 0

        Loop: DO jp = 1, jppj
          ipr = nprkx(jp)
          ir  = nspi(ipr,1)
          IF ( madvtr(ir) == 0 .AND. moffam(ir) == 0 ) CYCLE Loop

          DO jm = 3, jpspj
            ip = nspi(ipr,jm)
            IF ( ip /= 0 ) THEN
              ifam(jm-2) = moffam(ip)
              itr(jm-2)  = madvtr(ip)
            ELSE
              ifam(jm-2) = 0
              itr(jm-2)  = 0
            ENDIF
          ENDDO

!         1.1 Scan the ratefile. Store an entry if; one of the products
!         is tracer/family and the reactant is a tracer/family. For each
!         store the reaction no. and reactant tracer no.
!         For families only, if the product family is the same
!         as the reactant family, we ignore the term since this would
!         appear on the main diagonal.

!         This is complicated a little if the species is of type 'FT'.
!         We need to check if that species is in its family or not
!         during IMPACT itself. But we always include these terms by
!         ensuring the if (itr) statements comes before the test for
!         family membership.

          DO j = 1, jpspj-2
            IF ( itr(j) /= 0 ) THEN
              gtr(itr(j))     = .true.
              iwork(0,1)      = iwork(0,1) + 1
              ientry          = iwork(0,1)
              iwork(ientry,1) = jp
              iwork(ientry,2) = itr(j)
            ELSE IF ( ifam(j) /= 0 .and. ifam(j) /= moffam(ir) ) THEN
              gtr(ifam(j))    = .true.
              iwork(0,1)      = iwork(0,1) + 1
              ientry          = iwork(0,1)
              iwork(ientry,1) = jp
              iwork(ientry,2) = ifam(j)
            ENDIF
          ENDDO
        ENDDO Loop


!       2.  Build lists.
!           ----- ------

        DO j = 1, jpctr
          if ( gtr(j) ) then
            nltrim(0,1) = nltrim(0,1) + 1
            nltrim(nltrim(0,1),1) = j
          endif
        ENDDO
        ntr = nltrim(0,1)

!       2.1  Sort work list; copy entries to the partial
!            derivative list so that the entries are in
!            row, then column, then reaction order when
!            considering the Jacobian matrix. First put
!            them in row and column order.

        ientries = iwork(0,1)
        index    = 1
        DO j = 1, ntr
          jt = nltrim(j,1)
          nltrim(j,3) = index
          DO ji = 1, ientries
            iprod = iwork(ji,2)
            IF ( jt == iprod ) THEN
              nlpdv(index,1) = iwork(ji,1)
              iwork(ji,1)    = 0
              iwork(ji,2)    = 0
              index          = index + 1
            ENDIF
          ENDDO
          nltrim(j,2) = index - nltrim(j,3)
        ENDDO

!       2.2  Now sort each column to group reactions that
!            have the same reactant together. Use iwork as
!            a temporary array while we do the sort.

        DO j = 1, ntr
          jt     = nltrim(j,1)
          inf    = nltrim(j,2)
          ipos   = nltrim(j,3)
          icount = 0
          DO ji = 1, inf
            ir = nlpdv(ipos+ji-1,1)
            IF ( ir /= 0 ) THEN
              irk  = nprkx(ir)
              is   = nspi(irk,1)
              irtr = madvtr(is)
              IF ( irtr == 0 ) irtr = moffam(is)
              icount = icount + 1
              iwork(icount,1) = nlpdv(ipos+ji-1,1)
              nlpdv(ipos+ji-1,1) = 0

!             See if there's any more with the same reactant.

              DO jii = ji+1, inf
                ir2   = nlpdv(ipos+jii-1,1)
                IF ( ir2  /=  0 ) then
                  irk2  = nprkx(ir2)
                  is2   = nspi(irk2,1)
                  irtr2 = madvtr(is2)
                  if ( irtr2 == 0 ) irtr2 = moffam(is2)
                  IF ( irtr2  ==  irtr ) THEN
                    icount = icount + 1
                    iwork(icount,1)     = nlpdv(ipos+jii-1,1)
                    nlpdv(ipos+jii-1,1) = 0
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO

          DO ji = 1, inf
            index = ipos + ji - 1
            nlpdv(index,1) = iwork(ji,1)
            nlpdv(index,2) = 1
            IF ( ji /= inf ) THEN
              DO inr=1,inf-ji
                ir    = iwork(ji+inr,1)
                ir2   = iwork(ji,1)
                irk   = nprkx(ir)
                irk2  = nprkx(ir2)
                is    = nspi(irk,1)
                is2   = nspi(irk2,1)
                irtr  = madvtr(is)
                irtr2 = madvtr(is2)
                IF ( irtr  == 0 ) irtr = moffam(is)
                IF ( irtr2 == 0 ) irtr2 = moffam(is2)
                IF ( irtr /= irtr2 ) THEN
                  EXIT
                ELSE
                  nlpdv(index,2) = nlpdv(index,2) + 1
                ENDIF
              ENDDO     ! inr
            ENDIF
          ENDDO         ! ji

        ENDDO           ! j

        IF (lhook) CALL dr_hook('ASAD_INIMPCT',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_INIMPCT
