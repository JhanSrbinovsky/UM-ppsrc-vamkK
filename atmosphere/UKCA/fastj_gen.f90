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
      SUBROUTINE FASTJ_GEN

      USE FASTJ_DATA
      USE FASTJ_MIE

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      INTEGER  :: id, id0, id1, im
      INTEGER  :: i, j, k, KW, ix, kx  

      REAL     :: sum0, sum1, sum2, sum3
      REAL     :: sum0x(NWB3,nlevs_mie), sum1x(NWB3,2*m_gauss)
      REAL     :: sum2x(2*m_gauss)
      REAL     :: deltau, d1, d2, surfac
      REAL     :: pmpm0(kpcx,m_gauss,2*m_gauss)
      REAL     :: pmpm(m_gauss,m_gauss,2*m_gauss)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



!     Declare S,W,U1,V1,W in here; AA only for ID=1, CC only for ID=ND,
!     Pass A,B,C1,H

      IF (lhook) CALL dr_hook('FASTJ_GEN',zhook_in,zhook_handle)
      DO im = 1,mfit
        DO i = 1,n
          DO ix = 1,kpcx
            pmpm0(ix,i,im) = pm(i,im)*pm0(ix,im)
          ENDDO
        ENDDO
      ENDDO
      DO im = 1,mfit
        DO j = 1,n
         DO i = 1,n
           pmpm(i,j,im) = pm(i,im)*pm(j,im)
         ENDDO
        ENDDO
      ENDDO

      DO i = 1,n
        DO id = 2,nd-1
          DO kw = nwb1,nwb3
            a(kw,i,id)   = emu(i)/(ztau(kw,id+1) - ztau(kw,id-1))
            c1(kw,i,id)  = -a(kw,i,id)
            sum0x(kw,id) = 0.0d0
          ENDDO
        ENDDO
        DO im = m,mfit
          DO kw = nwb1,nwb2
            DO ix = 1,kpcx
              kx              = (ix-1)*nwb2
              sum1x(kx+KW,im) = pmpm0(ix,i,im)
            ENDDO
          ENDDO
        ENDDO
 
!       hard-coded version for m_gauss == 4

        IF (m_gauss == 4)   THEN
          DO id = 2,nd-1,2
            DO kw = nwb1,nwb3
              sum0x(kw,id) = sum0x(kw,id) + pomega(kw,2,id)*sum1x(kw,2)
              sum0x(kw,id) = sum0x(kw,id) + pomega(kw,4,id)*sum1x(kw,4)
              sum0x(kw,id) = sum0x(kw,id) + pomega(kw,6,id)*sum1x(kw,6)
              sum0x(kw,id) = sum0x(kw,id) + pomega(kw,8,id)*sum1x(kw,8)
            ENDDO
          ENDDO
          DO id = 3,nd-1,2
            DO kw = nwb1,nwb3
              sum0x(kw,id) = sum0x(kw,id) + pomega(kw,1,id)*sum1x(kw,1)
              sum0x(kw,id) = sum0x(kw,id) + pomega(kw,3,id)*sum1x(kw,3)
              sum0x(kw,id) = sum0x(kw,id) + pomega(kw,5,id)*sum1x(kw,5)
              sum0x(kw,id) = sum0x(kw,id) + pomega(kw,7,id)*sum1x(kw,7)
            ENDDO
          ENDDO
          
        ELSE
        
          DO id = 2,nd-1,2
            DO im = m+1,mfit,2
              DO kw = nwb1,nwb3
                sum0x(kw,id) = sum0x(kw,id)+pomega(kw,im,id)*sum1x(kw,im)
              ENDDO
            ENDDO
          ENDDO
          do id = 3,nd-1,2
            DO im = m,mfit,2
              DO kw = nwb1,nwb3
                sum0x(kw,id) = sum0x(kw,id)+pomega(kw,im,id)*sum1x(kw,im)
              ENDDO
            ENDDO
          ENDDO

        END IF
        
        DO id = 2,nd-1
          DO kw = nwb1,nwb3
            h(kw,i,id) = sum0x(kw,id)*fz(kw,id)
          ENDDO
        ENDDO

        DO j = 1,i
          DO id = 2,nd-1
            DO kw = nwb1,nwb3
              sum0x(kw,id) = 0.0d0
            ENDDO
          ENDDO
          DO im = m,mfit
            sum2x(im) = pmpm(i,j,im)
          ENDDO

!         hard-coded version for m_gauss == 4

          IF (m_gauss == 4)   THEN
            DO id = 2,nd-1,2
              DO kw = nwb1,nwb3
                sum0x(kw,id) = sum0x(kw,id) + pomega(kw,2,id)*sum2x(2)
                sum0x(kw,id) = sum0x(kw,id) + pomega(kw,4,id)*sum2x(4)
                sum0x(kw,id) = sum0x(kw,id) + pomega(kw,6,id)*sum2x(6)
                sum0x(kw,id) = sum0x(kw,id) + pomega(kw,8,id)*sum2x(8)
              ENDDO
            ENDDO
            DO id = 3,nd-1,2
              do kw = nwb1,nwb3
                sum0x(kw,id) = sum0x(kw,id) + pomega(kw,1,id)*sum2x(1)
                sum0x(kw,id) = sum0x(kw,id) + pomega(kw,3,id)*sum2x(3)
                sum0x(kw,id) = sum0x(kw,id) + pomega(kw,5,id)*sum2x(5)
                sum0x(kw,id) = sum0x(kw,id) + pomega(kw,7,id)*sum2x(7)
              ENDDO
            ENDDO
          ELSE
            DO id = 2,nd-1,2
              DO im = m+1,mfit,2
                DO kw = nwb1,nwb3
                  sum0x(kw,id) = sum0x(kw,id) + pomega(kw,im,id)*sum2x(im)
                ENDDO
              ENDDO
            ENDDO
            DO id = 3,nd-1,2
              DO im = m,mfit,2
                DO kw = nwb1,nwb3
                  sum0x(kw,id) = sum0x(kw,id) + pomega(kw,im,id)*sum2x(im)
                ENDDO
              ENDDO
            ENDDO
          END IF

          DO id = 2,nd-1
            DO kw = nwb1,nwb3
                b(i,j,kw,id) =  - sum0x(kw,id)*wt(j)
                b(j,i,kw,id) =  - sum0x(kw,id)*wt(i)
            ENDDO
          ENDDO
        ENDDO  ! End of J loop
        
        DO id = 2,nd-1
          DO kw = nwb1,nwb3
            b(i,i,kw,id) = b(i,i,kw,id) + 1.0d0
          ENDDO
        ENDDO

      ENDDO   ! End of I loop

!     Calculate generic 2nd-order terms for boundaries
!     for id=1, ID0=1,ID1=2; for id=2, ID0=ND,ID1=ND-1

      DO id = 1,2
        IF (id == 1) THEN
          id0 = 1
          id1 = 2
        ELSE
          id0 = nd
          id1 = nd-1
        ENDIF

        DO i = 1,n
          DO kw = nwb1,nwb3    !!  Lazy?
            ix   =((kw-1)/nwb2)+1
            sum0 = 0.0d0
            sum1 = 0.0d0
            sum2 = 0.0d0
            sum3 = 0.0d0
            DO im = m,mfit,2
              sum0 = sum0 + pomega(kw,im,id0)*pmpm0(ix,i,im)
              sum2 = sum2 + pomega(kw,im,id1)*pmpm0(ix,i,im)
            ENDDO
            DO im = m+1,mfit,2
              sum1 = sum1 + pomega(kw,im,id0)*pmpm0(ix,i,im)
              sum3 = sum3 + pomega(kw,im,id1)*pmpm0(ix,i,im)
            ENDDO
            h(kw,i,id0) = 0.5d0*(sum0*fz(kw,id0) + sum2*fz(kw,id1))
            a(kw,i,id0) = 0.5d0*(sum1*fz(kw,id0) + sum3*fz(kw,id1))
            DO j = 1,i
              sum0 = 0.0d0
              sum1 = 0.0d0
              sum2 = 0.0d0
              sum3 = 0.0d0
              DO im = m,mfit,2
                sum0 = sum0 + pomega(kw,im,id0)*pmpm(i,j,im)
                sum2 = sum2 + pomega(kw,im,id1)*pmpm(i,j,im)
              ENDDO
              DO im = m+1,mfit,2
                sum1 = sum1 + pomega(kw,im,id0)*pmpm(i,j,im)
                sum3 = sum3 + pomega(kw,im,id1)*pmpm(i,j,im)
              ENDDO
              s(kw,i,j)     = - sum2*wt(j)
              s(kw,j,i)     = - sum2*wt(i)
              w(kw,i,j,id)  = - sum1*wt(j)
              w(kw,j,i,id)  = - sum1*wt(i)
              u1(kw,i,j,id) = - sum3*wt(j)
              u1(kw,j,i,id) = - sum3*wt(i)
              sum0          = 0.5d0*(sum0 + sum2)
              b(i,j,kw,id0) = - sum0*wt(j)
              b(j,i,kw,id0) = - sum0*wt(i)
            ENDDO
            s(kw,i,i)     = s(kw,i,i)     + 1.0d0
            w(kw,i,i,id)  = w(kw,i,i,id)  + 1.0d0
            u1(kw,i,i,id) = u1(kw,i,i,id) + 1.0d0
            b(i,i,kw,id0) = b(i,i,kw,id0) + 1.0d0
          ENDDO
        ENDDO

        DO j = 1,n
          sum0 = 1.d0/emu(j)
          DO I = 1,N
            DO kw = nwb1,nwb3
              s(kw,i,j) = s(kw,i,j)*sum0
            ENDDO
          ENDDO
        ENDDO

        DO i = 1,n
          DO kw = nwb1,nwb3
            sum0 = 0.0d0
            DO j = 1,n
              sum0 = sum0 + s(kw,i,j)*a(kw,j,id0)
            ENDDO
            c1(kw,i,id0) = sum0
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,n
            DO kw = nwb1,nwb3
              sum0 = 0.0d0
              sum2 = 0.0d0
              DO k = 1,n
                sum0 = sum0 + s(kw,j,k)*w(kw,k,i,id)
                sum2 = sum2 + s(kw,j,k)*u1(kw,k,i,id)
              ENDDO
              a(kw,j,id0) = sum0
              v1(kw,j)    = sum2
            ENDDO
          ENDDO
          
!         Joining loops gives slightly different answers - care...

          DO J = 1,n
            DO kw = nwb1,nwb3
              w (kw,j,i,id) = a (kw,j,id0)
              u1(kw,j,i,id) = v1(kw,j)
            ENDDO
           ENDDO

         ENDDO

       ENDDO  ! End of id loop

!------Upper boundary, 2nd-order, C-matrix is full (CC)

       DO kw = nwb1,nwb3
         deltau = ztau(kw,2) - ztau(kw,1)
         d2     = 0.25d0*deltau
         DO I = 1,N
           d1 = emu(i)/deltau
           DO J=1,N
             b(i,j,kw,1) = b(i,j,kw,1) + d2*w(kw,i,j,1)
             cc(kw,i,j)  = d2*u1(kw,i,j,1)
           ENDDO
           b(i,i,kw,1) = b(i,i,kw,1) + d1
           cc(kw,i,i)  = cc(kw,i,i) - d1
           h(kw,i,1)   = h(kw,i,1) + 2.0d0*d2*c1(kw,i,1)
           a(kw,i,1)   = 0.0d0
         ENDDO
       ENDDO

!------Lower boundary, 2nd-order, A-matrix is full (AA)

       DO kw = nwb1,nwb3

         ix     =((kw-1)/nwb2)+1
         deltau = ztau(kw,nd) - ztau(kw,nd-1)
         d2     = 0.25d0*deltau
         surfac = 4.0d0*zrefl(ix)/(1.0d0 + zrefl(ix))

         DO i = 1,n
           d1         = emu(i)/deltau
           h(kw,i,nd) = h(kw,i,nd) - 2.0d0*d2*c1(kw,i,nd)
           sum0       = 0.0d0
           DO j = 1,n
             sum0 = sum0 + w(kw,i,j,2)
           ENDDO
           sum0 = d1 + d2*sum0
           sum1 = surfac*sum0
           DO j = 1,n
             b(i,j,kw,nd)=b(i,j,kw,nd)+d2*w(kw,i,j,2)-sum1*emu(j)*wt(j)
           ENDDO
           b(i,i,kw,nd) = b(i,i,kw,nd) + d1
           h(kw,i,nd) = h(kw,i,nd) + sum0*zflux(kw)
           DO J=1,N
             aa(kw,i,j) = - d2*u1(kw,i,j,2)
           ENDDO
           aa(kw,i,i) = aa(kw,i,i) + d1
           c1(kw,i,nd) = 0.0d0
         ENDDO
       ENDDO

       IF (lhook) CALL dr_hook('FASTJ_GEN',zhook_out,zhook_handle)
       RETURN
       END SUBROUTINE FASTJ_GEN
