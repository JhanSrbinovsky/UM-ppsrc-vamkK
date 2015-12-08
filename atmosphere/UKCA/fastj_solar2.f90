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
!   Routine to set up SZA for given lat, lon and time
!
!     timej    Offset in hours from start of timestep to time J-values
!              required for - take as half timestep for mid-step Js.
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
      SUBROUTINE FASTJ_SOLAR2(timej,sinlat,longitude,lcal360)
      
      USE conversions_mod, ONLY: pi_over_180, pi
      
      USE   FASTJ_DATA,    ONLY: iday,tau,szamax,               &
     &                                  sza_2d,SZAFAC_2d
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT  NONE


      LOGICAL,INTENT(IN)             :: lcal360
      
      REAL,INTENT(IN)                :: timej
      REAL,DIMENSION(:,:),INTENT(IN) :: sinlat,longitude 


!     Local variables

      REAL                                          :: sindec
      REAL                                          :: soldek
      REAL                                          :: cosdec
      REAL,DIMENSION(SIZE(sinlat,1),SIZE(sinlat,2)) :: sollat
      REAL,DIMENSION(SIZE(sinlat,1),SIZE(sinlat,2)) :: coslat
      REAL,DIMENSION(SIZE(sinlat,1),SIZE(sinlat,2)) :: cosz
      REAL,DIMENSION(SIZE(sinlat,1),SIZE(sinlat,2)) :: loct

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('FASTJ_SOLAR2',zhook_in,zhook_handle)
      IF (lcal360) THEN
        sindec = 0.3978*SIN(2*pi*(REAL(iday)-80.)/360.)
      ELSE
        sindec = 0.3978*SIN(2*pi*(REAL(iday)-80.)/365.)
      ENDIF

      soldek = ASIN(sindec)
      cosdec = COS(soldek)
      sollat = ASIN(sinlat)
      coslat = COS(sollat)

      loct   = 2*pi*(tau+timej)/24.-pi + longitude
      cosz   = cosdec*coslat*COS(loct) + sindec*sinlat
      sza_2d = ACOS(cosz)/Pi_Over_180

!     Filter for Large-SZA columns

      WHERE (sza_2d <= szamax)
         SZAFAC_2d = 1.d0
      ELSEWHERE
         SZAFAC_2d = 0.d0
         sza_2d    = 89.d0   ! Reduce SZA to avoid airmass problems
      ENDWHERE
 
      IF (lhook) CALL dr_hook('FASTJ_SOLAR2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_SOLAR2
