! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer. 
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution. 
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Fast-j routine for calculating online photolysis rates
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
        SUBROUTINE FASTJ_CLDSRF(odcol)

        USE     FASTJ_DATA
        USE     FASTJ_MIE,        ONLY: jadsub

        USE yomhook, ONLY: lhook, dr_hook
        USE parkind1, ONLY: jprb, jpim

        IMPLICIT NONE

        INTEGER :: i, j, k, ncfup, ncflo, ix, lx

        REAL    :: odsum
        REAL    :: odmax
        REAL    :: zlbatm
        REAL    :: cfcos
        REAL    :: cftau
        REAL    :: odadd
        REAL    :: odcol(2,kpcx,lpar)
        REAL    :: odtot(lpar)

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

                             
!-----------------------------------------------------------------------
!       This routine sets cloud and surface properties
!-----------------------------------------------------------------------
!       rflect   Surface albedo (Lambertian)
!       odmax    Maximum allowed optical depth, above which they are scaled
!       odcol    Optical depth at each model level
!       odsum    Column optical depth
!       nlbatm   Level of lower photolysis boundary - usually surface
!-----------------------------------------------------------------------

        IF (lhook) CALL dr_hook('FASTJ_CLDSRF',zhook_in,zhook_handle)

!       Default lower photolysis boundary as bottom of level 1
        nlbatm = 1

!       Set surface albedo

        DO ix=1,kpcx
          RFLECT(ix) = MAX(0.0,MIN(1.0,SA(nsl(1,ix),nsl(2,ix))))
        ENDDO

!       Initialise sub-division switch

        DO i = 1,nc
          jadsub(i)=0
        ENDDO

        odcol  = 0.0
        jadsub =   0 

!       Loop over column

        DO ix = 1,kpcx

!       Scale optical depths as appropriate 
!       limit column to 'odmax'

          odmax = 200.0
          odsum =   0.0
          DO i = 1,lwepar
            odcol(1,ix,i) = ODW(nsl(1,ix),nsl(2,ix),i)
            odcol(2,ix,i) = ODI(nsl(1,ix),nsl(2,ix),i)
            odsum    = odsum + odcol(1,ix,i) + odcol(2,ix,i)
          ENDDO
          DO i = lwepar+1,lpar
            odcol(1,ix,i) = 0.0
            odcol(2,ix,i) = 0.0
          ENDDO
          IF (odsum > odmax) THEN
            odsum = odmax/odsum
            DO i = 1,lpar
              odcol(1,ix,i) = odcol(1,ix,i)*odsum
              odcol(2,ix,i) = odcol(2,ix,i)*odsum
            ENDDO
            odsum = odmax
          ENDIF
          DO i = 1,lpar
            odtot(i) = odcol(1,ix,i) + odcol(2,ix,i)
          ENDDO

!        Set sub-division switch if appropriate

          odadd = 0.0
          LOOP_i: DO i = nb-1,1,-1
            k = 2*i
            odadd = odadd + odtot(i)
            IF (odadd > 0.0 .AND. odtot(i) /= 0.0 .AND.               &
     &                            dtausub >   0.0) THEN
              IF (odadd <= dtausub) THEN
                jadsub(k)   = 1
                jadsub(k-1) = 1
              ELSEIF (odadd > dtausub) THEN
                jadsub(k)   = 1
                EXIT LOOP_i 
              ENDIF
            ENDIF
          ENDDO LOOP_i
        ENDDO

        IF (lhook) CALL dr_hook('FASTJ_CLDSRF',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE FASTJ_CLDSRF
