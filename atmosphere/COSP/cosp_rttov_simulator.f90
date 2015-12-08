! (c) British Crown Copyright 2008-2013, the Met Office.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
!
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!  Description:  Module that prepares the inputs in the correct format,
!                and makes the call to the RTTOV simulator.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP
!

MODULE MOD_COSP_RTTOV_SIMULATOR
  USE cosp_constants_mod
  USE cosp_types_mod
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE COSP_RTTOV_SIMULATOR ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_RTTOV_SIMULATOR(gbx,y)

  ! Arguments
  TYPE(cosp_gridbox),INTENT(IN)  :: gbx ! Gridbox info
  TYPE(cosp_rttov),INTENT(INOUT) :: y   ! RTTOV output

  ! some local variables for profile conversions etc.
  REAL, PARAMETER :: eps    =  0.622
  REAL, PARAMETER :: Mdry   =  28.966
  REAL, PARAMETER :: Mo3    =  47.9983
  REAL, PARAMETER :: Mco2   =  44.0096
  REAL, PARAMETER :: Mch4   =  16.0426
  REAL, PARAMETER :: Mn2o   =  44.0129
  REAL, PARAMETER :: Mco    =  28.0102
  INTEGER, PARAMETER :: MaxLim  =  100

  ! Local variables
  INTEGER :: Npoints
  REAL :: sh(gbx%Npoints, gbx%Nlevels)
  REAL :: pp(gbx%Npoints, gbx%Nlevels)
  REAL :: tt(gbx%Npoints, gbx%Nlevels)
  REAL :: o3(gbx%Npoints, gbx%Nlevels)

  REAL :: co2,ch4,n2o,co
  REAL :: tt_surf(gbx%Npoints) ! 1.5 m T
  REAL :: sh_surf(gbx%Npoints) ! 1.5 m q 
  INTEGER :: nloop,rmod,il
  INTEGER :: istart,istop
  INTEGER :: nprof,nlevels

  Nlevels = gbx%Nlevels
  Npoints = gbx%Npoints
  ! Reverting Levels from TOA to surface
  sh  = gbx%sh(:,Nlevels:1:-1) 
  pp  = gbx%p(:,Nlevels:1:-1) / 100.
  tt  = gbx%t(:,Nlevels:1:-1) 
  o3  = gbx%mr_ozone(:,Nlevels:1:-1)

  ! FIXME: 1.5 m T and q should be added to input
  tt_surf  =  tt(:, Nlevels)
  sh_surf  =  sh(:, Nlevels)

  !Converting Specific Humidity to PPMV
  sh  =  ( sh / ( sh + eps * ( 1. - sh ) ) ) * 1e6

  !Converting Mass mixing ratio of other trace gases to ppmv
  o3   =  ( Mdry / Mo3  ) *     o3  * 1e6
  co2  =  ( Mdry / Mco2 ) * gbx%co2 * 1e6
  ch4  =  ( Mdry / Mch4 ) * gbx%ch4 * 1e6
  n2o  =  ( Mdry / Mn2o ) * gbx%n2o * 1e6
  co   =  ( Mdry / Mco  ) * gbx%co  * 1e6

  !! RTTOV can handle only about 100 profiles at a time
  !! (FIXME: Check this with Roger)
  !! So we are putting a loop of 100
  nloop  =  Npoints / MaxLim
  rmod   =  MOD( Npoints, MaxLim )

  IF ( rmod .ne. 0 ) THEN
     nloop = nloop + 1
  ENDIF

  !! looping over MaxLim number of profiles
  DO il = 1, nloop
     istart  =  (il - 1) * MaxLim + 1
     istop   =  MIN(il * MaxLim, Npoints) 

     IF ( ( il .eq. nloop ) .AND. ( rmod .ne. 0 ) ) THEN
        nprof   =  rmod
     ELSE
        nprof   =  MaxLim
     ENDIF

  ENDDO

END SUBROUTINE COSP_RTTOV_SIMULATOR

END MODULE MOD_COSP_RTTOV_SIMULATOR
