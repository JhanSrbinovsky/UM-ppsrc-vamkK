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
!   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Sol'n of inhomogeneous Rayleigh scattering atmosphere.
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!
!    Takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field!
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
       SUBROUTINE FASTJ_MIESCT
       USE     FASTJ_DATA
       USE     FASTJ_MIE
       USE parkind1, ONLY: jprb, jpim
       USE yomhook, ONLY: lhook, dr_hook
       IMPLICIT NONE

       INTEGER          :: i, id, im, KW, ix, kx

       REAL             :: cmeq1,pm00(2*m_gauss)

       INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
       INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
       REAL(KIND=jprb)               :: zhook_handle

       IF (lhook) CALL dr_hook('FASTJ_MIESCT',zhook_in,zhook_handle)

!      Fix scattering to 4 Gauss pts = 8-stream

!DEPENDS ON: fastj_gaussp
       CALL FASTJ_GAUSSP (n,emu,wt)

!      Solve eqn of R.T. only for first-order M=1
!      ZFLUX = (ZU0*FZ(ND)*ZREFL+FBOT)/(1.0d0+ZREFL)

       DO kw = nwb1,nwb3
         kx        = ((kw-1)/nwb2)+1
         zflux(kw) = (zu0(kx)*fz(kw,nd)*zrefl(kx))/(1.0d0+zrefl(kx))
       ENDDO

      DO i = 1,n
!DEPENDS ON: fastj_legnd0
        CALL FASTJ_LEGND0 (emu(i),pm00,mfit)
        DO im = m,mfit
          pm(i,im) = pm00(im)
        ENDDO
      ENDDO

      CMEQ1 = 0.25D0
      DO ix = 1,kpcx
!DEPENDS ON: fastj_legnd0
        CALL FASTJ_LEGND0 (-zu0(ix),pm00,mfit)
        DO im = m,mfit
          pm0(ix,im) = cmeq1*pm00(im)
        ENDDO
      ENDDO

!DEPENDS ON: fastj_blkslv
      CALL FASTJ_BLKSLV

      DO id = 1,nd,2
        DO kw = nwb1,nwb3
          fj(kw,id) = 4.0d0*fj(kw,id) + fz(kw,id)
        ENDDO
      ENDDO

      IF (lhook) CALL dr_hook('FASTJ_MIESCT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_MIESCT
