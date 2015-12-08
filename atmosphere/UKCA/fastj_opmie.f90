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
!  NEW Mie code for J's, only uses 8-term expansion, 4-Gauss pts
!  Currently allow up to NP aerosol phase functions (at all altitudes) to
!  be associated with optical depth AER(1:NC) = aerosol opt.depth @ 1000 nm
!
!  Pick Mie-wavelength with phase function and Qext:
!
!  01 RAYLE = Rayleigh phase
!  02 ISOTR = isotropic
!  03 ABSRB = fully absorbing 'soot', wavelength indep.
!  04 S_Bkg = backgrnd stratos sulfate (n=1.46,log-norm:r=.09um/sigma=.6)
!  05 S_Vol = volcanic stratos sulfate (n=1.46,log-norm:r=.08um/sigma=.8)
!  06 W_H01 = waterhaze (H1/Deirm.)(n=1.335,gamma:r-mode=0.1um /alpha=2)
!  07 W_H04 = waterhaze (H1/Deirm.)(n=1.335,gamma:r-mode=0.4um /alpha=2)
!  08 W_C02 = watercloud (C1/Deirm.)(n=1.335,gamma:r-mode=2.0um /alpha=6)
!  09 W_C04 = watercloud (C1/Deirm.)(n=1.335,gamma:r-mode=4.0um /alpha=6)
!  10 W_C08 = watercloud (C1/Deirm.)(n=1.335,gamma:r-mode=8.0um /alpha=6)
!  11 W_C13 = watercloud (C1/Deirm.)(n=1.335,gamma:r-mode=13.3um /alpha=6)
!  12 W_L06 = watercloud (Lacis) (n=1.335, r-mode=5.5um / alpha=11/3)
!  13 Ice-H = hexagonal ice cloud (Mishchenko)
!  14 Ice-I = irregular ice cloud (Mishchenko)
!
!  Choice of aerosol index MIEDX is made in SET_AER; optical depths are
!  apportioned to the AER array in SET_PROF
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
      SUBROUTINE FASTJ_OPMIE
      USE     FASTJ_DATA
      USE     FASTJ_MIE
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE

!  FUNCTION RAYLAY(WAVE)---RAYLEIGH CROSS-SECTION for wave > 170 nm
!       WSQI = 1.E6/(WAVE*WAVE)
!       REFRM1 = 1.0E-6*(64.328+29498.1/(146.-WSQI)+255.4/(41.-WSQI))
!       RAYLAY = 5.40E-21*(REFRM1*WSQI)**2

!     DTAUX    Local optical depth of each CTM level
!     PIRAY    Contribution of Rayleigh scattering to extinction
!     PIAER    Contribution of Aerosol scattering to extinction
!     TTAU     Optical depth of air vertically above each point 
!              (to top of atm)
!     FTAU     Attenuation of solar beam
!     POMEGA   Scattering phase function
!     FMEAN    Mean actinic flux at desired levels


      CHARACTER(LEN=72) :: cmessage
      INTEGER :: errcode                ! Variable passed to ereport

      INTEGER       :: jndlev(lpar)
      INTEGER       :: jaddlv(nc)
      INTEGER       :: jaddto(nc+1)
      INTEGER       :: KW,km(NWB2),i,j,k,l,ix,j1,kx,jx

      REAL ::  QXMIE(NWB2,MX),SSALB(NWB2,MX)
      REAL ::  xlo2,xlo3,xlray,zk,taudn(NWB3),tauup(NWB3),zk2
      REAL ::  XQO2(NB),XQO3(NB),POMEGAJ(NWB3,2*m_gauss,NC+1)
      REAL ::  DTAUX(NWB3,NB),PIRAY(NWB3,NB),PIAER(NWB3,MX,NB)
      REAL ::  TTAU(NWB3,NC+1),FTAU(NWB3,NC+1)
      REAL ::  ftaulog(NWB3),dttau(NWB3),dpomega(NWB3,2*m_gauss)
      REAL ::  ftaulog2(NWB3),dttau2(NWB3),dpomega2(NWB3,2*m_gauss)
      REAL ::  fastj_flint

      REAL                                      :: XQO3_l,XQO2_l
      REAL,DIMENSION(kpcx,mx)                   :: XLAER
      REAL,DIMENSION(kpcx)                      :: XLRAY_l
      REAL,DIMENSION(nwb3,mfit,mx)              :: paa_v
      REAL,DIMENSION(nwb3,mfit)                 :: paa_1
      REAL,DIMENSION(size(amf,1),nwb2)          :: xltau
      INTEGER,SAVE                              :: nlevs_mie_clear=nlevs_mie

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     Pick nearest Mie wavelength, no interpolation--------------

      IF (lhook) CALL dr_hook('FASTJ_OPMIE',zhook_in,zhook_handle)
      DO kw = nwb1,nwb2
        km(kw) = 1
        IF( wl(kw) > 355.d0 ) km(kw) = 2
        IF( wl(kw) > 500.d0 ) km(kw) = 3
!       IF( wl(kw) > 800.d0 ) km(kw) = 4  !drop the 1000 nm wavelength

!       For Mie code scale extinction at 1000 nm to wavelength WAVEL (QXMIE)
        
        DO i = 1,mx
          qxmie(kw,i) = qaa(km(kw),miedx(i))/qaa(4,miedx(i))
          ssalb(kw,i) = ssa(km(kw),miedx(i))
        ENDDO
      ENDDO

!     Reinitialize arrays

      DO j = 1,nc+1
        DO kw = nwb1,nwb3
          ttau(kw,j) = 0.d0
          ftau(kw,j) = 0.d0
        ENDDO
      ENDDO

!     Set up total optical depth over each CTM level, DTAUX

      j1 = nlbatm
      DO kw = nwb1,nwb2
        DO j = 1,nb
          DO ix = 1,kpcx
            kx=(ix-1)*nwb2+kw

!DEPENDS ON: fastj_flint
            xqo3_l = FASTJ_FLINT(tj(ix,j),tqq(1,2),tqq(2,2),tqq(3,2),    &
     &                             qo3(kw,1),qo3(kw,2),qo3(kw,3))

!DEPENDS ON: fastj_flint
            xqo2_l = FASTJ_FLINT(tj(ix,j),tqq(1,1),tqq(2,1),tqq(3,1),    &
     &                             qo2(kw,1),qo2(kw,2),qo2(kw,3))
            xlo3        = do3(ix,j)*xqo3_l
            xlo2        = dm(ix,j)*xqo2_l*0.20948d0
            xlray_l(ix) = dm(ix,j)*qrayl(kw)
            dtaux(kx,j) = xlo3+xlo2+xlray_l(ix)
            
          END DO

          DO I = 1,MX
            DO ix = 1,kpcx
              kx          = (ix-1)*nwb2+kw
              xlaer(ix,i) = aer(i,ix,j)*qxmie(kw,i)

!             Total optical depth from all elements
              dtaux(kx,j) = dtaux(kx,j)+xlaer(ix,i)

            END DO
          END DO
          DO ix = 1,kpcx
            kx          = (ix-1)*nwb2+kw
            piray(kx,j) = xlray_l(ix)/dtaux(kx,j)
          END DO
          DO I = 1,MX
            DO ix = 1,kpcx
              kx=(ix-1)*NWB2+KW
              
!             Fractional extinction for Rayleigh scattering 
!             and each aerosol type

              piaer(kx,i,j) = ssalb(kw,i)*xlaer(ix,i)/dtaux(kx,j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!     Define the scattering phase fn. with mix of Rayleigh(1) & Mie(MIEDX)
!     No. of quadrature pts fixed at 4 (m_gauss), expansion of phase fn @ 8

      DO i = 1,mfit
        DO kw = nwb1,nwb3
          kx          = MOD(kw-1,nwb2)+1
          paa_1(kw,i) = paa(i,KM(kx),1)
        END DO
      END DO
      DO i = 1,mfit
        DO k = 1,mx
          DO kw = nwb1,nwb3
            kx            = MOD(kw-1,nwb2)+1
            paa_v(kw,i,k) = paa(i,km(kx),miedx(k))
          END DO
        END DO
      END DO
      
      DO j = j1,nb
        DO i = 1,mfit
          DO kw = nwb1,nwb3
            pomegaj(KW,i,j) = piray(kw,j)*paa_1(kw,i)
          END DO
          DO k = 1,mx
            DO kw = nwb1,nwb3
              pomegaj(kw,i,j) = pomegaj(kw,i,j) +                  &
     &                          piaer(kw,k,j)*paa_v(kw,i,k)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!     Calculate attenuated incident beam EXP(-TTAU/U0) and flux on surface

      DO j = j1,nb
        xltau = 0.0D0
        DO kw = nwb1,nwb2
          DO i = 1,nb
            DO ix = 1,kpcx
              kx = (ix-1)*NWB2
              IF (AMF(ix,J,J)> 0.0D0) THEN
                xltau(ix,kw) = xltau(ix,kw) + dtaux(kx+kw,i)*amf(ix,i,j)
              END IF
            END DO
          ENDDO
        END DO

        DO kw = nwb1,nwb2
          DO ix = 1,kpcx
            kx = (ix-1)*NWB2
            IF (amf(ix,j,j) > 0.0D0) THEN
              IF (xltau(ix,kw) > 450.d0) THEN  ! for compilers with 
                                               ! no underflow trapping
                ftau(kx+kw,j) = 0.d0
              ELSE
                ftau(kx+kw,j) = EXP(-xltau(ix,kw))
              ENDIF
            ELSE
              ftau(kx+kw,j) = 0.0d0
            ENDIF
          END DO
        ENDDO
      ENDDO

      DO ix = 1,kpcx
        kx = (ix-1)*NWB2
        IF (u0(ix) > 0.D0) then
          DO kw = nwb1,nwb2
            zflux(kx+kw) = u0(ix)*ftau(kx+kw,j1)*rflect(ix)              &
     &                                           /(1.d0+rflect(ix))
          ENDDO
        ELSE
          DO kw=nwb1,nwb2
            zflux(kx+kw) = 0.d0
          ENDDO
        ENDIF
      ENDDO


!     Take optical properties on CTM layers and convert to a photolysis
!     level grid corresponding to layer centres and boundaries. This is
!     required so that J-values can be calculated for the centre of CTM
!     layers; the index of these layers is kept in the jndlev array.


!     Set lower boundary and levels to calculate J-values at 

      J1 = 2*J1-1
      DO j = 1,lpar
        jndlev(j)=2*j
      ENDDO
      DO j = 1,NC
        jaddlv(j)=0
      ENDDO

!     Calculate column optical depths above each level, TTAU

      do kw = nwb1,nwb3
        ttau(kw,nc+1) = 0.0d0
      ENDDO
      DO j = nc,j1,-1
        i = (j+1)/2
        DO kw = nwb1,nwb3
          ttau(kw,j) = ttau(kw,j+1) + 0.5d0*dtaux(kw,i)
          jaddlv(j)  = MAX(jaddlv(j),INT(0.5d0*dtaux(kw,i)/dtaumax))
        ENDDO

!       Subdivide cloud-top levels if required

        IF (jadsub(j) > 0) THEN
          jadsub(j) = MIN(jaddlv(j)+1,nint(dtausub))*(nint(dsubdiv)-1)
          jaddlv(j) = jaddlv(j)+jadsub(j)
        ENDIF
      ENDDO

!     Calculate attenuated beam, FTAU, level boundaries then level centres

      DO kw = nwb1,nwb3
        ftau(kw,nc+1)=1.0d0
      ENDDO
      DO J = nc-1,j1,-2
        i = (j+1)/2
        DO kw = nwb1,nwb3
          ftau(kw,j) = ftau(kw,i)
        ENDDO
      ENDDO
      DO j = nc,j1,-2
        DO kw = nwb1,nwb3
          ftau(kw,j) = SQRT(ftau(kw,j+1)*ftau(kw,j-1))
        ENDDO
      ENDDO

!     Calculate scattering properties, level centres then level boundaries
!     using an inverse interpolation to give correctly-weighted values

      DO j = NC,J1,-2
        DO i = 1,MFIT
          DO KW = NWB1,NWB3
            pomegaj(KW,i,j) = pomegaj(KW,i,j/2)
          ENDDO
        ENDDO
      ENDDO
      DO j = J1+2,nc,2
        DO KW = NWB1,NWB3
          taudn(KW) = ttau(KW,j-1)-ttau(KW,j)
          tauup(KW) = ttau(KW,j)-ttau(KW,j+1)
        ENDDO
        DO i=1,MFIT
          DO KW=NWB1,NWB3
            pomegaj(KW,i,j) = (pomegaj(KW,i,j-1)*taudn(KW) +             &
     &              pomegaj(KW,i,j+1)*tauup(KW))/(taudn(KW)+tauup(KW))
          ENDDO
        ENDDO
      ENDDO

!     Define lower and upper boundaries

      DO i=1,MFIT
        DO KW=NWB1,NWB3
          pomegaj(KW,i,J1)   = pomegaj(KW,i,J1+1)
          pomegaj(KW,i,nc+1) = pomegaj(KW,i,nc)
        ENDDO
      ENDDO


!     Calculate cumulative total and define levels we want J-values at.
!     Sum upwards for levels, and then downwards for Mie code readjustments.
!
!     jaddlv(i)   Number of new levels to add between (i) and (i+1)
!     jaddto(i)   Total number of new levels to add to and above level (i)
!     jndlev(j)   Level needed for J-value for CTM layer (j)

!     Reinitialize level arrays

      DO j = 1,nc+1
        jaddto(j) = 0
      ENDDO

      jaddto(J1) = jaddlv(J1)
      DO j = J1+1,nc
        jaddto(j) = jaddto(j-1) + jaddlv(j)
      ENDDO
      IF ((jaddto(nc)+nc) > nl) THEN
         errcode=1
         cmessage='opmie_1'
         WRITE(6,*)  jaddto(nc)+nc, 'NL',NL

         CALL EREPORT("FASTJ",errcode,cmessage)
      ENDIF

      DO i = 1,lpar
        jndlev(i) = jndlev(i)+jaddto(jndlev(i)-1)
      ENDDO
      jaddto(nc) = jaddlv(nc)
      DO j = nc-1,J1,-1
        jaddto(j) = jaddto(j+1) + jaddlv(j)
      ENDDO

!---------------------SET UP FOR MIE CODE-------------------------------
!
!  Transpose the ascending TTAU grid to a descending ZTAU grid.
!  Double the resolution - TTAU points become the odd points on the
!  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
!  Odd point added at top of grid for unattenuated beam   (Z='inf')
!  
!        Surface:   TTAU(1)   now use ZTAU(2*NC+1)
!        Top:       TTAU(NC)  now use ZTAU(3)
!        Infinity:            now use ZTAU(1)
!
!  Mie scattering code only used from surface to level NC
!------------------------------------------------------------------------

!     Initialise all Fast-J optical property arrays

      DO k = 1,nlevs_mie_clear
        DO i = 1,MFIT
          DO kw = nwb1,nwb3
            pomega(KW,i,k) = 0.d0
          ENDDO
        ENDDO
      ENDDO
      DO k = 1,nlevs_mie_clear
        DO kw = nwb1,nwb3
          ztau(KW,k) = 0.d0
          fz(KW,k)   = 0.d0
        ENDDO
      ENDDO

!     Ascend through atmosphere transposing grid and adding extra points

      DO j = J1,nc+1
        k = 2*(nc+1-j) + 2*jaddto(j)+1
        DO kw = nwb1,nwb3
          ztau(kw,k)= ttau(kw,j)
          fz(kw,k)  = ftau(kw,j)
          DO i = 1,mfit
            pomega(kw,i,k) = pomegaj(kw,i,j)
          ENDDO
        ENDDO
      ENDDO

!------------------------------------------------------------------------
!  Insert new levels, working downwards from the top of the atmosphere
!  to the surface (down in 'j', up in 'k'). This allows ztau and pomega
!  to be incremented linearly (in a +ve sense), and the flux fz to be
!  attenuated top-down (avoiding problems where lower level fluxes are
!  zero).
!
!    zk        fractional increment in level
!    dttau     change in ttau per increment    (linear, positive)
!    dpomega   change in pomega per increment  (linear)
!    ftaulog   change in ftau per increment    (exponential, normally < 1)
!
!------------------------------------------------------------------------

      DO j=nc,J1,-1
        zk = 0.5d0/(1.d0+dble(jaddlv(j)-jadsub(j)))
        DO kw = nwb1,nwb3
          dttau(KW) = (ttau(KW,j)-ttau(KW,j+1))*zk
          DO i=1,MFIT
            dpomega(KW,i) = (pomegaj(KW,i,j)-pomegaj(KW,i,j+1))*zk
          ENDDO
        ENDDO

!       Filter attenuation factor - set minimum at 1.0d-05

        DO KW = NWB1,NWB3
          IF (ftau(KW,j+1) == 0.d0) then
            ftaulog(KW) = 0.d0
          ELSE
            ftaulog(KW) = ftau(KW,j)/ftau(KW,j+1)
            IF (ftaulog(KW) < 1.d-150) THEN
              ftaulog(KW) = 1.0d-05
            ELSE
              ftaulog(KW) = EXP(log(ftaulog(KW))*zk)
            ENDIF
          ENDIF
        ENDDO
        k = 2*(nc-j+jaddto(j)-jaddlv(j))+1   !  k at level j+1
        l = 0

!       Additional subdivision of first level if required

        IF (jadsub(j) /= 0) THEN
          l   = jadsub(j)/nint(dsubdiv-1)
          zk2 = 1.d0/dsubdiv
          DO kw = nwb1,nwb3
            dttau2(KW)   = dttau(KW)*zk2
            ftaulog2(KW) = ftaulog(KW)**zk2
            DO i = 1,MFIT
              dpomega2(KW,i) = dpomega(KW,i)*zk2
            ENDDO
          ENDDO
          DO jx = 1,2*(jadsub(j)+l)
            DO kw = nwb1,nwb3
              ztau(KW,k+1) = ztau(KW,k) + dttau2(KW)
              fz(KW,k+1)   = fz(KW,k)*ftaulog2(KW)
              DO i=1,MFIT
                pomega(KW,i,k+1) = pomega(KW,i,k) + dpomega2(KW,i)
              ENDDO
            ENDDO
            k = k+1
          ENDDO
        ENDIF
        l = 2*(jaddlv(j)-jadsub(j)-l)+1

!       Add values at all intermediate levels

        DO jx = 1,l
          DO KW = NWB1,NWB3
            ztau(KW,k+1) = ztau(KW,k) + dttau(KW)
            fz(KW,k+1)   = fz(KW,k)*ftaulog(KW)
            DO i = 1,MFIT
              pomega(KW,i,k+1) = pomega(KW,i,k) + dpomega(KW,i)
            ENDDO
          ENDDO
          k = k+1
        ENDDO

      ENDDO

!     Update total number of levels and check doesn't exceed nlevs_mie

      ND = 2*(NC+jaddto(J1)-J1)  + 3
      IF (nd > nlevs_mie) THEN
        errcode=1
        cmessage='opmie_2'
        WRITE(6,*) ND, 'nlevs_mie',nlevs_mie

        CALL EREPORT("FASTJ",errcode,cmessage)
      ENDIF
      nlevs_mie_clear=ND

!     Add boundary/ground layer to ensure no negative J's caused by
!     too large a TTAU-step in the 2nd-order lower b.c.

      DO kw = nwb1,nwb3
        ix            = ((kw-1)/nwb2)+1
        ztau(kw,nd+1) = ztau(kw,nd)*1.000005d0
        ztau(kw,nd+2) = ztau(kw,nd)*1.000010d0
        zk            = MAX(abs(u0(ix)),0.01e0)
        zk            = dexp(-ztau(kw,nd)*5.d-6/zk)
        fz(kw,nd+1)   = fz(kw,nd)*zk
        fz(kw,nd+2)   = fz(kw,nd+1)*zk
        DO i = 1,mfit
          pomega(kw,i,nd+1) = pomega(kw,i,nd)
          pomega(kw,i,nd+2) = pomega(kw,i,nd)
        ENDDO
      ENDDO
      ND = ND+2

      DO ix = 1,kpcx
        zu0(ix)   = u0(ix)
        zrefl(ix) = rflect(ix)
      ENDDO

!DEPENDS ON: fastj_miesct
      CALL FASTJ_MIESCT

!     Accumulate attenuation for selected levels

      l = 2*(NC+jaddto(J1))+3
      DO j = 1,jpcl
        k = l-(2*jndlev(j))
        IF (k <= ND-2) THEN
          DO kw = nwb1,nwb3
            kx        = MOD(kw-1,nwb2)+1
            fff(kw,j) = fff(kw,j) + fl(kx)*fj(kw,k)
          ENDDO
        ENDIF
      ENDDO

      IF (lhook) CALL dr_hook('FASTJ_OPMIE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_OPMIE
