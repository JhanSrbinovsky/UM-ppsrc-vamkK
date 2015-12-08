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
!                and makes the call to the CALIPSO simulator.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP

MODULE MOD_COSP_LIDAR
  USE cosp_constants_mod
  USE cosp_types_mod
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_LIDAR ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_LIDAR(gbx,sgx,sghydro,y)

  ! Arguments
  TYPE(cosp_gridbox),INTENT(IN) :: gbx  ! Gridbox info
  TYPE(cosp_subgrid),INTENT(IN) :: sgx  ! Subgrid info
  TYPE(cosp_sghydro),INTENT(IN) :: sghydro  ! Subgrid info for hydrometeors
  TYPE(cosp_sglidar),INTENT(INOUT) :: y ! Subgrid output

  ! Local variables
  INTEGER :: i
  REAL :: presf(sgx%Npoints, sgx%Nlevels + 1)
  REAL,DIMENSION(sgx%Npoints, sgx%Nlevels) :: lsca,mr_ll,mr_li,mr_cl,mr_ci
  REAL,DIMENSION(sgx%Npoints, sgx%Nlevels) :: beta_tot,tau_tot
  REAL,DIMENSION(sgx%Npoints, PARASOL_NREFL)  :: refle

  presf(:,1:sgx%Nlevels) = gbx%ph
  presf(:,sgx%Nlevels + 1) = 0.0
  lsca = gbx%tca-gbx%cca
  DO i=1,sgx%Ncolumns
      ! Temporary arrays for simulator call
      mr_ll(:,:) = sghydro%mr_hydro(:,i,:,I_LSCLIQ)
      mr_li(:,:) = sghydro%mr_hydro(:,i,:,I_LSCICE)
      mr_cl(:,:) = sghydro%mr_hydro(:,i,:,I_CVCLIQ)
      mr_ci(:,:) = sghydro%mr_hydro(:,i,:,I_CVCICE)
! DEPENDS ON: lidar_simulator
      CALL lidar_simulator(sgx%Npoints, sgx%Nlevels, 4 , PARASOL_NREFL,        &
                 LIDAR_UNDEF, gbx%p, presf, gbx%T, mr_ll, mr_li, mr_cl, mr_ci, &
                 gbx%Reff(:,:,I_LSCLIQ), gbx%Reff(:,:,I_LSCICE),               &
                 gbx%Reff(:,:,I_CVCLIQ), gbx%Reff(:,:,I_CVCICE),               &
                 sgx%frac_out, gbx%lidar_ice_type, y%beta_mol, beta_tot,       &
                 tau_tot, refle )
      y%beta_tot(:,i,:) = beta_tot(:,:)
      y%tau_tot(:,i,:)  = tau_tot(:,:)
      y%refl(:,i,:)     = refle(:,:)
  ENDDO

END SUBROUTINE COSP_LIDAR

END MODULE MOD_COSP_LIDAR
