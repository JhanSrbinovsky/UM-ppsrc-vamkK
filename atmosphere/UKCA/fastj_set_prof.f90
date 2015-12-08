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
!  Routine to set up atmospheric profiles required by Fast-J using a
!  doubled version of the level scheme used in the CTM. First pressure
!  and z* altitude are defined, then O3 and T are taken from the supplied
!  climatology and integrated to the CTM levels (may be overwritten with
!  values directly from the CTM, if desired) and then black carbon and
!  aerosol profiles are constructed.
!
!     pj       Pressure at boundaries of model levels (hPa)
!     z        Altitude of boundaries of model levels (cm)
!     odcol    Optical depth at each model level (1=water,2=ice)
!     masfac   Conversion factor for pressure to column density
!
!     TJ       Temperature profile on model grid
!     DM       Air column for each model level (molecules.cm-2)
!     DO3      Ozone column for each model level (molecules.cm-2)
!     DBC      Mass of Black Carbon at each model level (g.cm-3)  !  .....!
!     PSTD     Approximate pressures of levels for supplied climatology
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
      SUBROUTINE FASTJ_SET_PROF
      USE    FASTJ_DATA
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

      INTEGER          :: i, k, ix

      REAL             :: odcol(2,kpcx,lpar),masfac

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      
!     First calculate pressure at boundaries of CTM levels (lpar+1)

      IF (lhook) CALL dr_hook('FASTJ_SET_PROF',zhook_in,zhook_handle)
      DO i = 1,nb
        DO ix = 1,kpcx
          pj(ix,i) = P(nsl(1,ix),nsl(2,ix),i)
        ENDDO
      ENDDO
      pj(:,nb+1) = 0.0

!     Set up cloud and surface properties

!DEPENDS ON: fastj_cldsrf
      CALL FASTJ_CLDSRF(odcol)

!     Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)

      masfac = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

      DO ix = 1,kpcx

!       Zero if dark

        IF (SZAFAC(ix) == 0.d0) THEN
          DO i = 1,lpar
            odcol(1,ix,i) = 0.d0
            odcol(2,ix,i) = 0.d0
          ENDDO
        ENDIF
 
        DO i = 1,lpar
          tj(ix,i)  = t_fastj(nsl(1,ix),nsl(2,ix),i)
          do3(ix,i) = fj_ozone(nsl(1,ix),nsl(2,ix),i)
        ENDDO
        tj (:,lpar+1) = tj (:,lpar)   ! Above top of model 
        do3(:,lpar+1) = do3(:,lpar)   ! Above top of model 
        
!       Add Aerosol Column - include aerosol types here. Currently 
!       use soot, water and ice; assume black carbon x-section of 
!       10 m2/g, independent of wavelength; assume limiting 
!       temperature for ice of -40 deg C.

        DO i = 1,lpar
          aer(1,ix,i) = odcol(1,ix,i)
          aer(2,ix,i) = odcol(2,ix,i)
          aer(3,ix,i) = sulphur(nsl(1,ix),nsl(2,ix),i)
          DO k = 4,MX
            aer(k,ix,i) = 0.0
          ENDDO
        ENDDO
 
        DO k=1,MX
          AER(k,ix,lpar+1) = 0.0
        ENDDO

      ENDDO

!     Calculate column quantities for Fast-J

      DO i=1,NB
        DO ix=1,kpcx
          dm(ix,i)  = (pj(ix,i)-pj(ix,i+1))*masfac
          do3(ix,i) = do3(ix,i)*dm(ix,i)
        ENDDO
      ENDDO

      IF (lhook) CALL dr_hook('FASTJ_SET_PROF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_SET_PROF
