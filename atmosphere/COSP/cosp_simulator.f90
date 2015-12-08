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


!  Description:  Module that controls the calls to individual simulators in
!                COSP, and converts the outputs to the correct units.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP
!

MODULE MOD_COSP_SIMULATOR
  USE cosp_constants_mod, ONLY: I_RADAR, I_LIDAR, I_ISCCP, I_MISR, I_MODIS, &
                                I_RTTOV, I_STATS, tsim
  USE cosp_types_mod
  USE MOD_COSP_RADAR
  USE MOD_COSP_LIDAR
  USE MOD_COSP_ISCCP_SIMULATOR
  USE MOD_COSP_MODIS_SIMULATOR
  USE MOD_COSP_MISR_SIMULATOR
  USE MOD_COSP_RTTOV_SIMULATOR
  USE MOD_COSP_STATS
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP_SIMULATOR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_SIMULATOR(gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,     &
                          misr,modis,rttov,stradar,stlidar)

  ! Arguments
  TYPE(cosp_gridbox),INTENT(INOUT) :: gbx      ! Grid-box inputs
  TYPE(cosp_subgrid),INTENT(IN) :: sgx      ! Subgrid inputs
  TYPE(cosp_sghydro),INTENT(IN) :: sghydro  ! Subgrid info for hydrometeors
  TYPE(cosp_config),INTENT(IN)  :: cfg      ! Configuration options
  TYPE(cosp_vgrid),INTENT(IN)   :: vgrid    ! Information on vert. grid of stats
  TYPE(cosp_sgradar),INTENT(INOUT) :: sgradar ! Output from radar simulator
  TYPE(cosp_sglidar),INTENT(INOUT) :: sglidar ! Output from lidar simulator
  TYPE(cosp_isccp),INTENT(INOUT)   :: isccp   ! Output from ISCCP simulator
  TYPE(cosp_misr),INTENT(INOUT)    :: misr    ! Output from MISR simulator
  TYPE(cosp_modis),INTENT(INOUT)   :: modis   ! Output from MODIS simulator
  TYPE(cosp_rttov),INTENT(INOUT)    :: rttov    ! Output from RTTOV
  TYPE(cosp_radarstats),INTENT(INOUT) :: stradar ! Summary statistics from radar
  TYPE(cosp_lidarstats),INTENT(INOUT) :: stlidar ! Summary statistics from lidar
  ! Local variables
  INTEGER :: i,j,k,isim
  LOGICAL :: inconsistent
  ! Timing variables
  INTEGER :: t0,t1

  t0 = 0
  t1 = 0

  inconsistent=.FALSE.

  !+++++++++ Radar model ++++++++++
  isim = I_RADAR
  IF (cfg%Lradar_sim) THEN
    CALL SYSTEM_CLOCK(t0)
    CALL cosp_radar(gbx,sgx,sghydro,sgradar)
    CALL SYSTEM_CLOCK(t1)
    tsim(isim) = tsim(isim) + (t1 -t0)
  ENDIF

  !+++++++++ Lidar model ++++++++++
  isim = I_LIDAR
  IF (cfg%Llidar_sim) THEN
    CALL SYSTEM_CLOCK(t0)
    CALL cosp_lidar(gbx,sgx,sghydro,sglidar)
    CALL SYSTEM_CLOCK(t1)
    tsim(isim) = tsim(isim) + (t1 -t0)
  ENDIF

  !+++++++++ ISCCP simulator ++++++++++
  isim = I_ISCCP
  IF (cfg%Lisccp_sim) THEN
    CALL SYSTEM_CLOCK(t0)
    CALL cosp_isccp_simulator(gbx,sgx,isccp)
    CALL SYSTEM_CLOCK(t1)
    tsim(isim) = tsim(isim) + (t1 -t0)
  ENDIF

  !+++++++++ MISR simulator ++++++++++
  isim = I_MISR
  IF (cfg%Lmisr_sim) THEN
    CALL SYSTEM_CLOCK(t0)
    CALL cosp_misr_simulator(gbx,sgx,misr)
    CALL SYSTEM_CLOCK(t1)
    tsim(isim) = tsim(isim) + (t1 -t0)
  ENDIF

  !+++++++++ MODIS simulator ++++++++++
  isim = I_MODIS
  IF (cfg%Lmodis_sim) THEN
    CALL SYSTEM_CLOCK(t0)
    CALL cosp_modis_simulator(gbx,sgx,sghydro,isccp, modis)
    CALL SYSTEM_CLOCK(t1)
    tsim(isim) = tsim(isim) + (t1 -t0)
  ENDIF

  !+++++++++ RTTOV ++++++++++ 
  isim = I_RTTOV
  IF (cfg%Lrttov_sim) THEN
    CALL SYSTEM_CLOCK(t0)
    CALL cosp_rttov_simulator(gbx,rttov)
    CALL SYSTEM_CLOCK(t1)
    tsim(isim) = tsim(isim) + (t1 -t0)
  ENDIF

  !+++++++++++ Summary statistics +++++++++++
  isim = I_STATS
  IF (cfg%Lstats) THEN
    CALL SYSTEM_CLOCK(t0)
    CALL cosp_stats(gbx,sgx,cfg,sgradar,sglidar,vgrid,stradar,stlidar)
    CALL SYSTEM_CLOCK(t1)
    tsim(isim) = tsim(isim) + (t1 -t0)
  ENDIF

  !+++++++++++ Change of units after computation of statistics +++++++++++
  ! This avoids using UDUNITS in CMOR

  ! Cloud fractions from 1 to %
  IF (cfg%Lclcalipso) THEN
    WHERE(stlidar%lidarcld /= R_UNDEF) stlidar%lidarcld = stlidar%lidarcld*100.0
  ENDIF
  IF (cfg%Lcltcalipso.OR.cfg%Lcllcalipso.OR.cfg%Lclmcalipso                    &
     .OR.cfg%Lclhcalipso) THEN
    WHERE(stlidar%cldlayer /= R_UNDEF) stlidar%cldlayer = stlidar%cldlayer*100.0
  ENDIF
  IF (cfg%Lclcalipso2) THEN
    WHERE(stradar%lidar_only_freq_cloud /= R_UNDEF)                            &
      stradar%lidar_only_freq_cloud = stradar%lidar_only_freq_cloud*100.0
  ENDIF

  IF (cfg%Lcltcalipsoliq.OR.cfg%Lcllcalipsoliq.OR.cfg%Lclmcalipsoliq.OR.       &
      cfg%Lclhcalipsoliq.OR.cfg%Lcltcalipsoice.OR.cfg%Lcllcalipsoice.OR.       &
      cfg%Lclmcalipsoice.OR.cfg%Lclhcalipsoice.OR.cfg%Lcltcalipsoun.OR.        &
      cfg%Lcllcalipsoun.OR.cfg%Lclmcalipsoun.OR.cfg%Lclhcalipsoun ) THEN
    WHERE(stlidar%cldlayerphase /= R_UNDEF)                                    &
      stlidar%cldlayerphase = stlidar%cldlayerphase*100.0
  ENDIF
  IF (cfg%Lclcalipsoliq.OR.cfg%Lclcalipsoice.OR.cfg%Lclcalipsoun) THEN
    WHERE(stlidar%lidarcldphase /= R_UNDEF)                                    &
      stlidar%lidarcldphase = stlidar%lidarcldphase*100.0
  ENDIF
  IF (cfg%Lclcalipsotmp.OR.cfg%Lclcalipsotmpliq.OR.cfg%Lclcalipsotmpice.OR.    &
      cfg%Lclcalipsotmpun) THEN
    WHERE(stlidar%lidarcldtmp /= R_UNDEF)                                      &
      stlidar%lidarcldtmp = stlidar%lidarcldtmp*100.0
  ENDIF

  IF (cfg%Lcltisccp) THEN
    WHERE(isccp%totalcldarea /= R_UNDEF)                                       &
      isccp%totalcldarea = isccp%totalcldarea*100.0
  ENDIF
  IF (cfg%Lclisccp) THEN
    WHERE(isccp%fq_isccp /= R_UNDEF)                                           &
      isccp%fq_isccp = isccp%fq_isccp*100.0
  ENDIF

  IF (cfg%LclMISR) THEN
    WHERE(misr%fq_MISR /= R_UNDEF) misr%fq_MISR = misr%fq_MISR*100.0
  ENDIF

  IF (cfg%Lcltlidarradar) THEN
    WHERE(stradar%radar_lidar_tcc /= R_UNDEF)                                  &
      stradar%radar_lidar_tcc = stradar%radar_lidar_tcc*100.0
  ENDIF

  IF (cfg%Lclmodis) THEN
    WHERE(modis%Optical_Thickness_vs_Cloud_Top_Pressure /= R_UNDEF)            &
      modis%Optical_Thickness_vs_Cloud_Top_Pressure =                          &
      modis%Optical_Thickness_vs_Cloud_Top_Pressure*100.0
  ENDIF
  IF (cfg%Lcltmodis) THEN
     WHERE(modis%Cloud_Fraction_Total_Mean /= R_UNDEF)                         &
       modis%Cloud_Fraction_Total_Mean = modis%Cloud_Fraction_Total_Mean*100.0
  ENDIF
  IF (cfg%Lclwmodis) THEN
     WHERE(modis%Cloud_Fraction_Water_Mean /= R_UNDEF)                         &
       modis%Cloud_Fraction_Water_Mean = modis%Cloud_Fraction_Water_Mean*100.0
  ENDIF
  IF (cfg%Lclimodis) THEN
     WHERE(modis%Cloud_Fraction_Ice_Mean /= R_UNDEF)                           &
       modis%Cloud_Fraction_Ice_Mean = modis%Cloud_Fraction_Ice_Mean*100.0
  ENDIF

  IF (cfg%Lclhmodis) THEN
     WHERE(modis%Cloud_Fraction_High_Mean /= R_UNDEF)                          &
       modis%Cloud_Fraction_High_Mean = modis%Cloud_Fraction_High_Mean*100.0
  ENDIF
  IF (cfg%Lclmmodis) THEN
     WHERE(modis%Cloud_Fraction_Mid_Mean /= R_UNDEF)                           &
       modis%Cloud_Fraction_Mid_Mean = modis%Cloud_Fraction_Mid_Mean*100.0
  ENDIF
  IF (cfg%Lcllmodis) THEN
     WHERE(modis%Cloud_Fraction_Low_Mean /= R_UNDEF)                           &
       modis%Cloud_Fraction_Low_Mean = modis%Cloud_Fraction_Low_Mean*100.0
  ENDIF

  ! Change pressure from hPa to Pa.
  IF (cfg%Lboxptopisccp) THEN
    WHERE(isccp%boxptop /= R_UNDEF) isccp%boxptop = isccp%boxptop*100.0
  ENDIF
  IF (cfg%Lpctisccp) THEN
    WHERE(isccp%meanptop /= R_UNDEF) isccp%meanptop = isccp%meanptop*100.0
  ENDIF


END SUBROUTINE COSP_SIMULATOR

END MODULE MOD_COSP_SIMULATOR
