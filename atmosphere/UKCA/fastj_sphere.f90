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
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!  beam (where tangent height is below altitude J-value desired at).
!
!     GMU     MU, cos(solar zenith angle)
!     RZ      Distance from centre of Earth to each point (cm)
!     RQ      Square of radius ratios
!     TANHT   Tangent height for the current SZA
!     XL      Slant path between points
!     AMF     Air mass factor for slab between level and level above
!     U0      cos(solar zenith angle*pi/180)
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
      SUBROUTINE FASTJ_SPHERE
      USE     FASTJ_DATA
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

      INTEGER                       :: i, j, ii, ix

      LOGICAL,DIMENSION(kpcx)       :: todo

      REAL,DIMENSION(kpcx)          :: gmu, xmu1, xmu2, xl, diff
      REAL,DIMENSION(kpcx)          :: tanht_l
      REAL,DIMENSION(kpcx,NB)       :: rz,rq
      REAL                          :: zbyr

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('FASTJ_SPHERE',zhook_in,zhook_handle)
      diff(:) = 0.d0

      DO i=1,NB
        DO ix=1,kpcx
          rz(ix,i) = rz_3d(nsl(1,ix),nsl(2,ix),i)
        ENDDO
      ENDDO

      gmu(:) = u0(1:kpcx)
      zbyr   = zzht/rad
      DO ii=2,nb
        rq(:,ii-1) = (rz(:,ii-1)/rz(:,ii))**2
      ENDDO
      
      WHERE (GMU(1:kpcx) < 0.0D0) 
        tanht_l(:) = rz(:,nlbatm)/sqrt(1.0d0-gmu(1:kpcx)**2)
      ELSEWHERE
        tanht_l(:) = rz(:,nlbatm)
      ENDWHERE

!     Go up from the surface calculating the slant paths 
!     between each level and the level above, and deriving 
!     the appropriate Air Mass Factor

      amf(:,:,:)=0.D0

      DO j=1,NB

!       Air Mass Factors all zero if below the tangent height

        todo = (rz(1:kpcx,j) >= tanht_l(:))

!       Ascend from layer J calculating AMFs

        xmu1(:)=ABS(gmu(1:kpcx))
        DO i=j,lpar
          WHERE (todo)
            xmu2 (:)          = dsqrt(1.0d0-rq(:,i)*(1.0d0-xmu1**2))
            xl   (:)          = rz(:,i+1)*xmu2-rz(:,i)*xmu1
            amf  (1:kpcx,i,j) = xl/(rz(:,i+1)-rz(:,i))
            xmu1 (:)          = xmu2(:)
          ENDWHERE
        END DO

!       Use function and scale height to provide AMF above top of model

        DO ix = 1,kpcx
          IF (todo(ix))  THEN
            amf(ix,nb,j) = airmas(xmu1(ix),zbyr)
          END IF
        END DO

!       Twilight case - Emergent Beam

        todo = todo .AND. GMU < 0.0D0

        WHERE (todo)
          xmu1(:) = ABS(gmu(1:kpcx))
        ENDWHERE

!       Descend from layer J

        DO ii = J-1,1,-1

          WHERE (todo)
            diff(:)=rz(:,ii+1)*dsqrt(1.0d0-xmu1**2)-rz(:,ii)
          ENDWHERE

          IF (ii == 1) THEN
            diff = MAX(diff,0.e0)   ! filter
          END IF

!         Tangent height below current level - beam passes through twice

          WHERE (diff < 0.0D0 .AND. todo)
            xmu2 (:)           = dsqrt(1.0d0-(1.0d0-xmu1**2)/rq(:,ii))
            xl   (:)           = ABS(rz(:,ii+1)*xmu1-rz(:,ii)*xmu2)
            amf  (1:kpcx,ii,j) = 2.d0*xl/(rz(:,ii+1)-rz(:,ii))
            xmu1 (:)           = xmu2(:)

!         Lowest level intersected by emergent beam

          ELSEWHERE(todo)
            xl  (:)           = rz(:,ii+1)*xmu1*2.0d0
            amf (1:kpcx,ii,j) = xl/(rz(:,ii+1)-rz(:,ii))
            todo              = .false.
          ENDWHERE

        ENDDO

      ENDDO

      IF (lhook) CALL dr_hook('FASTJ_SPHERE',zhook_out,zhook_handle)
      RETURN

      CONTAINS

       FUNCTION AIRMAS(Ux,H)
! Inlined air mass factor function for top of atmosphere
       IMPLICIT NONE
         REAL :: AIRMAS
         REAL :: Ux,H

         INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
         INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
         REAL(KIND=jprb)               :: zhook_handle

       IF (lhook) CALL dr_hook('AIRMAS',zhook_in,zhook_handle)
         AIRMAS = (1.0d0+H)/SQRT(Ux*Ux+2.0d0*H*(1.0d0-                  &
              0.6817d0*EXP(-57.3d0*ABS(Ux)/SQRT(1.0d0+5500.d0*H))/      &
                                               (1.0d0+0.625d0*H)))
       
       IF (lhook) CALL dr_hook('AIRMAS',zhook_out,zhook_handle)
       RETURN
       ENDFUNCTION AIRMAS
      END SUBROUTINE FASTJ_SPHERE
