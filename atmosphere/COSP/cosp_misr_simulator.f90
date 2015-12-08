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
!                and makes the call to the MISR simulator.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP
!


MODULE MOD_COSP_MISR_SIMULATOR
  USE cosp_constants_mod
  USE cosp_types_mod
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------- SUBROUTINE COSP_MISR_SIMULATOR -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_MISR_SIMULATOR(gbx,sgx,y)
  IMPLICIT NONE
  ! Arguments
  TYPE(cosp_gridbox),INTENT(IN) :: gbx  ! Gridbox info
  TYPE(cosp_subgrid),INTENT(IN) :: sgx  ! Subgridbox info
  TYPE(cosp_misr),INTENT(INOUT) :: y    ! MISR simulator output

  ! Local variables 
  INTEGER :: Nlevels,Npoints
  REAL :: dtau_s(gbx%Npoints, gbx%Nlevels)
  REAL :: dtau_c(gbx%Npoints, gbx%Nlevels)
  REAL :: at(gbx%Npoints, gbx%Nlevels)
  REAL :: frac_out(gbx%Npoints, gbx%Ncolumns, gbx%Nlevels)
  INTEGER :: sunlit(gbx%Npoints)
  REAL :: zfull(gbx%Npoints, gbx%Nlevels) !  height (in meters) of full model levels (i.e. midpoints)
                                          !  zfull(npoints,1)    is    top level of model
                                          !  zfull(npoints,nlev) is bottom level of model


  Nlevels = gbx%Nlevels
  Npoints = gbx%Npoints
  ! Levels from TOA to surface
  zfull  = gbx%zlev(:,Nlevels:1:-1)
  at     = gbx%T(:,Nlevels:1:-1) 
  dtau_s = gbx%dtau_s(:,Nlevels:1:-1) 
  dtau_c = gbx%dtau_c(:,Nlevels:1:-1) 
  frac_out(1:Npoints,:,1:Nlevels) = sgx%frac_out(1:Npoints,:,Nlevels:1:-1)
  sunlit = INT(gbx%sunlit)

! DEPENDS ON: MISR_simulator
  CALL MISR_simulator(gbx%npoints,gbx%nlevels,gbx%ncolumns,                    &
                     sunlit,zfull,at,dtau_s,dtau_c,frac_out, R_UNDEF,          &
                     y%fq_MISR,y%MISR_dist_model_layertops,y%MISR_meanztop,    &
                     y%MISR_cldarea)

END SUBROUTINE COSP_MISR_SIMULATOR

END MODULE MOD_COSP_MISR_SIMULATOR
