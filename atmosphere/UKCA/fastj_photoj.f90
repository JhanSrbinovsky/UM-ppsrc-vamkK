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
!    new FAST J-Value code, troposphere only (mjprather 6/96)
!    uses special wavelength quadrature spectral data (jv_spec.dat)
!    that includes only 289 nm - 800 nm  (later a single 205 nm add-on)
!    uses special compact Mie code based on Feautrier/Auer/Prather vers.
!
!     zpj      External array providing J-values to main CTM code
!     timej    Offset in hours from start of timestep to time J-values
!              required for - take as half timestep for mid-step Js.
!     solf     Solar distance factor, for scaling; normally given by:
!                      1.0-(0.034*cos(real(iday-186)*2.0*pi/365.))
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
      SUBROUTINE FASTJ_PHOTOJ(zpj,timej)
      USE  FASTJ_DATA
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

      REAL,INTENT(IN)                             :: timej      
      REAL,DIMENSION(:,:,:,:),INTENT(OUT)         :: zpj

!     Local variables

      INTEGER          :: i,j,ix,lx

      REAL             :: solf

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('FASTJ_PHOTOJ',zhook_in,zhook_handle)
      DO j = 1,jppj
        DO i = 1,jjpnl    !! May not need this at all now...
          zj(i,j) = 0.D0
        ENDDO
      ENDDO
      
!     Skip calculation if all in dark

      IF (ALL(SZAFAC(1:kpcx) <= 0.001d0)) THEN 
        IF (lhook) CALL dr_hook('FASTJ_PHOTOJ',zhook_out,zhook_handle)
        RETURN
      END IF

!     Set up profiles on model levels

!DEPENDS ON: fastj_set_prof
      CALL FASTJ_SET_PROF

!DEPENDS ON: fastj_jvalue
      CALL FASTJ_JVALUE

 
!     Include variation in distance to sun
!     solf=1.d0-(0.034d0*cos(dble(iday-186)*2.d0*pi/365.d0))

      solf = 1.d0

!DEPENDS ON: fastj_jratet
      CALL FASTJ_JRATET(solf)


!  "zj" updated in JRATET - pass this back to ASAD as "zpj"
!CDIR nodep

      DO j = 1,jppj
        DO i = 1,jpcl
          DO ix = 1,kpcx
            lx = ix+(i-1)*kpcx
            zpj(nsl(1,ix),nsl(2,ix),i,j)= zj(lx,j)*SZAFAC(ix)
          ENDDO
        ENDDO
      ENDDO

      IF (lhook) CALL dr_hook('FASTJ_PHOTOJ',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_PHOTOJ
