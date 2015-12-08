! *****************************COPYRIGHT*******************************
! (c) [University of California] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

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
!   Fast-jx routine for calculating online photolysis rates
!   Preliminary version of the code used for development.
!
!     zpj      External array providing J-values to main CTM code
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
      SUBROUTINE FASTJX_PHOTOJ(zpj)

      USE  FASTJX_DATA
      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      REAL, INTENT(INOUT)     :: zpj(:,:,:,:)      ! photolysis rates

! Local variables
      INTEGER :: jxtra(lpar*2+3)                   ! Array containing extra levels 

! Temporary variables
      REAL :: wave                                 ! wavelength
      REAL :: ttt                                  ! temperature
      REAL :: xqo2                                 ! O2 cross section
      REAL :: xqo3                                 ! O3 cross section

! Optical depths from o2/o3 and rayleigh scattering
      REAL :: odabs, odray

! Optical density of each layer
      REAL :: dtaux((lpar+1),w_) 

! Solar flux deposited in layer L
      REAL :: flxd((lpar+1),w_) 

! Scattering phase function 
      REAL :: pomegax(sw_phases,(lpar+1),w_)

! Single scattering albedo
      REAL :: ssa(sw_band_aer, (lpar+1))

! 
      REAL :: sleg(sw_phases,  sw_band_aer, (lpar+1))

! mean actinic flux(diff+direct) at model mid-layer
      REAL :: fjact(lpar,w_)

! diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
      REAL :: fjtop(w_)

! diffuse flux onto surface (<0 by definition)
      REAL :: fjbot(w_)

! direct/solar flux onto surface  (<0 by definition)
      REAL :: fsbot(w_)

! diffuse flux across top of model layer L
      REAL :: fjflx(lpar,w_)

! sum of solar flux deposited in atmos
      REAL :: flxd0(w_)

! Mean actinic flux
      REAL :: fff(w_, jpcl)

! Solar distance factor. Initially set to 1.
      REAL :: solf=1.0

! Photolysis rates (before mapping onto ASAD arrays)
      REAL :: valjl(jpcl,njval)

      REAL :: qxmie(sw_band_aer,mx)

      REAL :: rfl(w_)

! Loop variables variables
      INTEGER :: i, j, k, l, ix, lx, kmie

      !-------------------------------------
      ! Shortcuts to speed up code
      REAL :: o3_ix(lpar+1)
      REAL :: dm_ix(lpar+1)
      REAL :: tz_ix(lpar+1)
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      !*************************************
      ! End of Header

      IF (lhook) CALL dr_hook('FASTJX_PHOTOJ',zhook_in,zhook_handle)

      ! Initialise to 0.
      zj(:,:) = 0.

      ! Skip calculation if all in dark
      !IF (ALL(szafac(1:kpcx) <= 0.001E0)) RETURN

      ! Convert GCM arrays to blocked arrays
      DO i = 1, (lpar+1)
        DO ix = 1,kpcx

          dm_block(ix, i) = dm_3d(nsl(1,ix),nsl(2,ix),i)
          o3_block(ix,i) = o3_3d(nsl(1,ix),nsl(2,ix),i)
          tz_block(ix,i) = tz_3d(nsl(1,ix),nsl(2,ix),i)

          odi_block(ix,i) = odi_3d(nsl(1,ix),nsl(2,ix),i)
          odw_block(ix,i) = odw_3d(nsl(1,ix),nsl(2,ix),i)
          ods_block(ix,i) = ods_3d(nsl(1,ix),nsl(2,ix),i)

          pz_block(ix,i) = pz_all(nsl(1,ix),nsl(2,ix),i)
          rz_block(ix,i) = rz_all(nsl(1,ix),nsl(2,ix),i)

          IF(i == (lpar+1)) THEN
            pz_block(ix,(lpar+2)) = pz_all(nsl(1,ix),nsl(2,ix),(lpar+2))
            rz_block(ix,(lpar+2)) = rz_all(nsl(1,ix),nsl(2,ix),(lpar+2))
          ENDIF
        ENDDO
      ENDDO

      !------------------------------------------------
      ! Calculate amf weighting fractions

      DO ix = 1, kpcx

        o3_ix(:) = o3_block(ix,:)
        dm_ix(:) = dm_block(ix,:)
        tz_ix(:) = tz_block(ix,:)

        ! Calculate spherical mass weighting factors
! DEPENDS ON: FASTJX_SPHERE
        CALL FASTJX_SPHERE(u0(ix), rad, rz_block(ix,:), zzht, amf2(:,:,ix), (lpar+1)) 

        ! Calculate the optical properties (opt-depth, single-scat-alb, phase-fn(1:8) 
        ! at the 5 std wavelengths 200-300-400-600-999 nm for cloud+aerosols

        !Initialise to 0.
        pomegax(:,:,:)   = 0.
        ssa(:,:)         = 0.  
        sleg(:,:,:)      = 0.
        od_block(ix,:,:) = 0. 
        od600(ix,:)      = 0.

        DO j=1, sw_band_aer

! Sum water and ice optical depths. Rescale to appropriate wavelengths
! Original calculation is valid for 200-690nm.
! Take this (fairly arbitrarily) to be true at 400nm (closest to middle of band)
! Effect is small as wavelength dependence is small 

          qxmie(j,1) =  (qaa(j,miedx(1))/qaa(4,miedx(1)))
          qxmie(j,2) =  (qaa(j,miedx(2))/qaa(4,miedx(2)))
          qxmie(j,3) =  (qaa(j,miedx(3))/qaa(4,miedx(3)))

          od_block(ix,j,:) = odw_block(ix,:)*qxmie(j,1)  &
                           + odi_block(ix,:)*qxmie(j,2)

          ! Calculate ssa (more elegant way of doing this?)
          ssa( j, :)    = odw_block(ix,:)*qxmie(j,1)*saa(j,miedx(1)) &
                        + odi_block(ix,:)*qxmie(j,2)*saa(j,miedx(2))      

          ! Calculate sleg (more elegant way of doing this?)
          DO k=1, sw_phases
            DO l=1, (lpar+1)
              sleg(k, j,l) = odw_block(ix,l)*qxmie(j,1)*saa(j,miedx(1))   &
                             * paa(k,j,miedx(1))                          &
                             + odi_block(ix,l)*qxmie(j,2)*saa(j,miedx(2)) &
                             * paa(k,j,miedx(2))    
            ENDDO
          ENDDO

        ENDDO

        ! Set 600 od to calculate extra layers
        od600(ix,:) = od_block(ix,4,:)

        ! Calculate where need extra layers due to thick clouds
        ! In principal could add aerosols to this calculation, but don't
! DEPENDS ON: FASTJX_EXTRAL
        CALL FASTJX_EXTRAL(od600(ix,:),(lpar+1),(2*lpar +2),N_,jtaumx,  &
                           atau, atau0, jxtra)

        !---------------------------------------------------------------------------
        ! Add aerosols to optical depth, ssa & sleg
        DO j=1, sw_band_aer
         
          od_block(ix,j,:) = od_block(ix,j,:) + ods_block(ix,:)*qxmie(j,3)

          ssa(j, :)    = ssa(j, :)    + ods_block(ix,:)*qxmie(j,3)      &
                           * saa(j,miedx(3)) 
 
          ! Calculate sleg (more elegant way of doing this?)
          DO k=1, sw_phases
            DO l=1, (lpar+1)
            sleg( k, j,l) = sleg(k, j,l) +  ods_block(ix,l)*qxmie(j,3)  &
                          * saa(j,miedx(3)) * paa(k,j,miedx(3)) 
            END DO
          END DO

          ! Convert from absolute to relative OD
          DO l=1, (lpar+1)

            IF(od_block(ix,j,l) > 0.) THEN
              ssa(j,l) = ssa(j,l)/ od_block(ix,j,l)

              DO k=1, sw_phases
                sleg(k, j,l) = sleg(k, j,l)/ od_block(ix,j,l)
              END DO 
            END IF
          END DO

        END DO

        !----------------------------------------------------------------------------
        ! Loop over wavelength bins
        DO k=1, w_

          wave = wl(k)
          !---Pick nearest Mie wavelength to get scattering properites------------
                              kmie=1   ! use 200 nm prop for <255 nm
          IF( wave > 255.E0 ) kmie=2   ! use 300 nm prop for 255-355 nm
          IF( wave > 355.E0 ) kmie=3   ! use 400 nm prop for 355-500 nm
          IF( wave > 500.E0 ) kmie=4
          IF( wave > 800.E0 ) kmie=5 

          !---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
          !---values at L1_=L_+1 are a pseudo/climatol layer above the top CTM layer (L_)
          DO l = 1,(lpar+1)

! Set temperature value for temp variable
            ttt     = tz_ix(l)

! interpolate o2 and o3 cross sections in temperature using flint function
            xqo2 = FLINT(ttt,tqq(1,1),tqq(2,1),tqq(3,1), qo2(k,1),qo2(k,2),qo2(k,3))
            xqo3 = FLINT(ttt,tqq(1,2),tqq(2,2),tqq(3,2), qo3(k,1),qo3(k,2),qo3(k,3))

            ! Calculate optical depth from sum of ozone and o2
            odabs = xqo3*o3_ix(l) +  xqo2*dm_ix(l)*0.20948E0
            ! Calculate rayleigh scattering optical depth
            odray = dm_ix(l)*qrayl(k) 
 
            ! sum cloud/aerosol, o2/o3 and rayleigh optical depths
            dtaux(l,k) = od_block(ix,kmie,l) + odabs + odray

            ! Loop over phase factors calculating scattering phase fn from clouds
            DO i=1,sw_phases
              pomegax(i,l,k) = sleg(i,kmie,l) *od_block(ix,kmie,l)
            END DO

            ! Add in rayleigh scattering into phase functions
            pomegax(1,l,k) = pomegax(1,L,K) + 1.0E0*odray
            pomegax(3,l,k) = pomegax(3,L,K) + 0.5E0*odray

            ! Loop over phase factors dividing by  optical depth of layers
            DO i=1,sw_phases
              pomegax(i,l,k) = pomegax(i,l,k)/dtaux(l,k)

            END DO

          END DO ! z

        END DO ! wavelengths       

        ! Set surface albedo to be constant across wavelengths        
        rfl(:) = sa_block(ix)

         ! Call Mie scattering routine
! DEPENDS ON: FASTJX_OPMIE
        CALL FASTJX_OPMIE (dtaux, pomegax,u0(ix),rfl, amf2(:,:,ix), jxtra, &
               fjact, fjtop, fjbot, fsbot, fjflx, flxd, flxd0, ix)

        ! Loop over wavelengths and model levels calculating mean actininc flux
        DO K = 1,W_
          DO L = 1,jpcl
            fff(K,L) = solf*fl(k)*fjact(l,k)
          END DO
        END DO  

! DEPENDS ON: FASTJX_JRATET
        CALL FASTJX_JRATET(pz_block(ix,:) ,tz_block(ix,:),fff, valjl)

!---map the J-values from fast-JX onto ASAD ones (use JIND & JFACTA)
        DO L = 1,jpcl
          DO J = 1, jppj
            IF (JIND(J) > 0) then 
              zj(L,J) = VALJL(L,JIND(J))*JFACTA(J)
            ELSE
              zj(L,J) = 0.E0
            END IF

            ! Convert to 4D array expected by ASAD
            zpj(nsl(1,ix),nsl(2,ix),l,j)= zj(l,j)*szafac(ix)       

          END DO
        END DO

      END DO ! xy

      IF (lhook) CALL dr_hook('FASTJX_PHOTOJ',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJX_PHOTOJ
