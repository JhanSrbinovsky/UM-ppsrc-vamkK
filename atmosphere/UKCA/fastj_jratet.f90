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
!  Calculate and print J-values. Note that the loop in this routine
!  only covers the jpcl levels actually needed by the CTM.
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
      SUBROUTINE FASTJ_JRATET(SOLF)
      USE FASTJ_DATA
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

      INTEGER   :: i, j, k, l, ix, lx, kx
      REAL      :: qo2tot, qo3tot, qo31d, qo33p, qqqt
      REAL      :: solf, tfact(kpcx), zt(kpcx), fastj_flint

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('FASTJ_JRATET',zhook_in,zhook_handle)
      DO I = 1,jpcl
        DO ix = 1,kpcx
          valj(ix,1) = 0.d0
          valj(ix,2) = 0.d0
          valj(ix,3) = 0.d0
          zt(ix)     = t_fastj(nsl(1,ix),nsl(2,ix),i)
        END DO
        DO k = nw1,nw2                       ! Using model 'T's here
          DO ix = 1,kpcx
            kx = (ix-1)*nw2

!     FFF    Actinic flux at each level for each wavelength bin
!     QQQ    Cross sections for species (read in in RD_TJPL)
!     SOLF   Solar distance factor, for scaling; normally given by:
!                      1.0-(0.034*cos(real(iday-186)*2.0*pi/365.))
!                      Assumes aphelion day 186, perihelion day 3.
!     TQQ    Temperatures at which QQQ cross sections supplied

!DEPENDS ON: fastj_flint
            qo2tot = FASTJ_FLINT(zt(ix),tqq(1,1),tqq(2,1),tqq(3,1),     &
     &                          qo2(k,1),qo2(k,2),qo2(k,3))

!DEPENDS ON: fastj_flint
            qo3tot = FASTJ_FLINT(zt(ix),tqq(1,2),tqq(2,2),tqq(3,2),     &
     &                          qo3(k,1),qo3(k,2),qo3(k,3))

!DEPENDS ON: fastj_flint
            qo31d = FASTJ_FLINT(zt(ix),tqq(1,3),tqq(2,3),tqq(3,3),     &
     &                        q1d(k,1),q1d(k,2),q1d(k,3))*qo3tot

            qo33p      = qo3tot - qo31d
            valj(ix,1) = valj(ix,1) + qo2tot*fff(kx+k,i)
            valj(ix,2) = valj(ix,2) + qo33p*fff(kx+k,i)
            valj(ix,3) = valj(ix,3) + qo31d*fff(kx+k,i)
          ENDDO
        ENDDO

!       Calculate remaining J-values with T-dep X-sections 

        DO J = 4,njval
          DO ix = 1,kpcx
            valj(ix,j) = 0.d0
            tfact(ix)  = 0.d0
          ENDDO
          L = jpdep(J)
          IF (tqq(2,j) > tqq(1,j)) THEN
            DO ix = 1,kpcx
              tfact(ix) = MAX(0.e0,MIN(1.e0,                    &
     &                   (zt(ix)-tqq(1,j))/(tqq(2,j)-tqq(1,j))))
            ENDDO
          ENDIF
          IF (L == 0) THEN
            DO k = nw1,nw2
              DO ix=1,kpcx
                kx = (ix-1)*nw2
                qqqt = qqq(k,1,j-3) +                             &
     &                (qqq(k,2,j-3) - qqq(k,1,j-3))*tfact(ix)
                valj(ix,j) = valj(ix,j) + qqqt*fff(kx+k,i)
              ENDDO
            ENDDO
          ELSE
            DO k = nw1,nw2
              DO ix = 1,kpcx
                kx = (ix-1)*NW2
                qqqt = qqq(k,1,j-3) +                             &
     &                (qqq(k,2,j-3) - qqq(k,1,j-3))*tfact(ix)
                valj(ix,j) = valj(ix,j) + qqqt*fff(kx+k,i)*       &
     &                (1.d0+zpdep(k,l)*(pj(ix,i)+pj(ix,i+1))*0.5d0)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        DO j = 1,jppj
          do ix = 1,kpcx
            lx       = ix+(i-1)*kpcx
            zj(lx,j) = VALJ(ix,jind(j))*jfacta(j)*SOLF
          ENDDO
        ENDDO
      ENDDO

      IF (lhook) CALL dr_hook('FASTJ_JRATET',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_JRATET
