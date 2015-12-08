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
!                and makes the call to the RADAR simulator.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP


MODULE MOD_COSP_RADAR
  USE cosp_constants_mod
  USE cosp_types_mod
  USE radar_simulator_types
  USE array_lib
  USE format_input
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE radar_simulator(hp,nprof,ngate,undef, &
        hgt_matrix,hm_matrix,re_matrix,Np_matrix, &
        p_matrix,t_matrix,rh_matrix, &
        Ze_non,Ze_ray,g_to_vol,a_to_vol,dBZe, &
        g_to_vol_in,g_to_vol_out)

        USE m_mrgrnk
        USE array_lib
        USE math_lib
        USE optics_lib
        USE radar_simulator_types
        IMPLICIT NONE

        ! ----- INPUTS -----  
        TYPE(class_param) :: hp

        INTEGER, INTENT(IN) :: nprof,ngate

        REAL undef
        REAL*8, DIMENSION(nprof,ngate), INTENT(IN) :: hgt_matrix, p_matrix, &
            t_matrix,rh_matrix
        REAL*8, DIMENSION(hp%nhclass,nprof,ngate), INTENT(IN) :: hm_matrix
        REAL*8, DIMENSION(hp%nhclass,nprof,ngate), INTENT(INOUT) :: re_matrix
        REAL*8, DIMENSION(hp%nhclass,nprof,ngate), INTENT(INOUT) :: Np_matrix

        ! ----- OUTPUTS -----
        REAL*8, DIMENSION(nprof,ngate), INTENT(OUT) :: Ze_non,Ze_ray, &
            g_to_vol,dBZe,a_to_vol
        ! ----- OPTIONAL -----
        REAL*8, OPTIONAL, DIMENSION(nprof,ngate) :: &
            g_to_vol_in,g_to_vol_out
     END SUBROUTINE radar_simulator
  END INTERFACE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_RADAR ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_RADAR(gbx,sgx,sghydro,z)
  IMPLICIT NONE

  ! Arguments
  TYPE(cosp_gridbox),INTENT(INOUT) :: gbx  ! Gridbox info
  TYPE(cosp_subgrid),INTENT(IN) :: sgx  ! Subgrid info
  TYPE(cosp_sghydro),INTENT(IN) :: sghydro  ! Subgrid info for hydrometeors
  TYPE(cosp_sgradar),INTENT(INOUT) :: z ! Output from simulator, subgrid

  ! Local variables 
  INTEGER :: & 
  nsizes            ! num of discrete drop sizes

  REAL*8, DIMENSION(:,:), ALLOCATABLE :: &
  g_to_vol ! integrated atten due to gases, r>v (dB)

  REAL*8, DIMENSION(:,:), ALLOCATABLE :: &
  Ze_non, &         ! radar reflectivity withOUT attenuation (dBZ)
  Ze_ray, &         ! Rayleigh reflectivity (dBZ)
  h_atten_to_vol, &     ! attenuation by hydromets, radar to vol (dB)
  g_atten_to_vol, &     ! gaseous atteunation, radar to vol (dB)
  dBZe, &           ! effective radar reflectivity factor (dBZ)
  hgt_matrix, &         ! height of hydrometeors (km)
  t_matrix, &                   !temperature (k)
  p_matrix, &                   !pressure (hPa)
  rh_matrix                     !relative humidity (%)

  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: &
  hm_matrix, &
! hydrometeor mixing ratio (g kg^-1)
  re_matrix, &
! effective radius (microns).   Optional. 0 ==> use Np_matrix or defaults
  Np_matrix
! total number concentration (kg^-1).   Optional 0==> use defaults 

  INTEGER, PARAMETER :: one = 1
  LOGICAL :: hgt_descending
  INTEGER :: pr,i,j,k,unt,ngate

! ----- main program settings ------

  ! Inputs to Quickbeam
  ALLOCATE(hgt_matrix(gbx%Npoints,gbx%Nlevels),                                &
           p_matrix(gbx%Npoints,gbx%Nlevels),                                  &
           t_matrix(gbx%Npoints,gbx%Nlevels),                                  &
           rh_matrix(gbx%Npoints,gbx%Nlevels))
  ALLOCATE(hm_matrix(gbx%Nhydro,gbx%Npoints,gbx%Nlevels))
  ALLOCATE(re_matrix(gbx%Nhydro,gbx%Npoints,gbx%Nlevels))
  ALLOCATE(Np_matrix(gbx%Nhydro,gbx%Npoints,gbx%Nlevels))

  ! Outputs from Quickbeam
  ALLOCATE(Ze_non(gbx%Npoints,gbx%Nlevels))
  ALLOCATE(Ze_ray(gbx%Npoints,gbx%Nlevels))
  ALLOCATE(h_atten_to_vol(gbx%Npoints,gbx%Nlevels))
  ALLOCATE(g_atten_to_vol(gbx%Npoints,gbx%Nlevels))
  ALLOCATE(dBZe(gbx%Npoints,gbx%Nlevels))

  ! Optional argument. It is computed and returned in the first call to
  ! radar_simulator, and passed as input in the rest
  ALLOCATE(g_to_vol(gbx%Npoints,gbx%Nlevels))

  ! Even if there is no unit conversion, they are needed for type conversion
  p_matrix   = gbx%p/100.0     ! From Pa to hPa
  hgt_matrix = gbx%zlev/1000.0 ! From m to km
  t_matrix   = gbx%T
  rh_matrix  = gbx%q
  re_matrix  = 0.0


  ! set flag denoting position of radar relative to hgt_matrix orientation
  ngate = SIZE(hgt_matrix,2)

  hgt_descending = hgt_matrix(1,1) > hgt_matrix(1,ngate)

  IF ((gbx%surface_radar == 1 .AND. hgt_descending) .OR.                       &
      (gbx%surface_radar == 0 .AND. (.NOT. hgt_descending))) THEN
    gbx%hp%radar_at_layer_one = .FALSE.
  ELSE
    gbx%hp%radar_at_layer_one = .TRUE.
  ENDIF

  ! ----- loop over subcolumns -----
  DO pr=1,sgx%Ncolumns

      !  NOTE:
      !  atmospheric profiles are the same within the same gridbox
      !  only hydrometeor profiles will be different for each subgridbox

         DO i=1,gbx%Nhydro
            hm_matrix(i,:,:) = sghydro%mr_hydro(:,pr,:,i)*1000.0 ! kg/kg to g/kg
            IF (gbx%use_reff) THEN
              re_matrix(i,:,:) = sghydro%Reff(:,pr,:,i)*1.e6 ! m to micron
              Np_matrix(i,:,:) = sghydro%Np(:,pr,:,i)        ! Units [#/kg]
            ENDIF
         ENDDO

      !   ----- call radar simulator -----
      IF (pr == 1) THEN ! Compute gaseous attenuation for all profiles
! DEPENDS ON: radar_simulator
        CALL radar_simulator(gbx%hp,gbx%Npoints,gbx%Nlevels,R_UNDEF,           &
           hgt_matrix,hm_matrix,re_matrix,Np_matrix,                           &
           p_matrix,t_matrix,rh_matrix,                                        &
           Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe,                   &
           g_to_vol_out=g_to_vol)
      ELSE ! Use gaseous atteunuation for pr = 1
! DEPENDS ON: radar_simulator
         CALL radar_simulator(gbx%hp,gbx%Npoints,gbx%Nlevels,R_UNDEF,          &
           hgt_matrix,hm_matrix,re_matrix,Np_matrix,                           &
           p_matrix,t_matrix,rh_matrix,                                        &
           Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe,                   &
           g_to_vol_in=g_to_vol)
      ENDIF

      ! store calculated dBZe values for later output/processing
      z%Ze_tot(:,pr,:)=dBZe(:,:)
  ENDDO !pr

  DEALLOCATE(hgt_matrix,p_matrix,t_matrix,rh_matrix)
  DEALLOCATE(hm_matrix,re_matrix, &
      Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe)
  DEALLOCATE(g_to_vol)
END SUBROUTINE COSP_RADAR

END MODULE MOD_COSP_RADAR
