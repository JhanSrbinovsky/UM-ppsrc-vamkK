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

!  Description:  Module that applies some quality control to the inputs,
!                calculates extra inputs if needed, and makes iterative
!                calls to COSP, the routine that calls the individual
!                simulators.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP

MODULE MOD_COSP
  USE cosp_types_mod
  USE mod_cosp_utils, ONLY: cosp_ereport
  USE MOD_COSP_SIMULATOR
  USE MOD_COSP_MODIS_SIMULATOR
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP ---------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,      &
                misr,modis,rttov,stradar,stlidar)

  ! Arguments
  INTEGER,INTENT(IN) :: overlap ! overlap type in SCOPS: 1=max,2=rand,3=max/rand
  INTEGER,INTENT(IN) :: Ncolumns
  TYPE(cosp_config),INTENT(IN) :: cfg   ! Configuration options
  TYPE(cosp_vgrid),INTENT(IN) :: vgrid   ! Information on vertical grid of stats
  TYPE(cosp_gridbox),INTENT(INOUT) :: gbx
  TYPE(cosp_subgrid),INTENT(INOUT) :: sgx   ! Subgrid info
  TYPE(cosp_sgradar),INTENT(INOUT) :: sgradar ! Output from radar simulator
  TYPE(cosp_sglidar),INTENT(INOUT) :: sglidar ! Output from lidar simulator
  TYPE(cosp_isccp),INTENT(INOUT)   :: isccp   ! Output from ISCCP simulator
  TYPE(cosp_misr),INTENT(INOUT)    :: misr    ! Output from MISR simulator
  TYPE(cosp_modis),INTENT(INOUT)   :: modis   ! Output from MODIS simulator
  TYPE(cosp_rttov),INTENT(INOUT)   :: rttov   ! Output from RTTOV
  TYPE(cosp_radarstats),INTENT(INOUT) :: stradar ! Summary statistics from radar
  TYPE(cosp_lidarstats),INTENT(INOUT) :: stlidar ! Summary statistics from lidar

  ! Local variables 
  INTEGER :: Npoints   ! Number of gridpoints
  INTEGER :: Nlevels   ! Number of levels
  INTEGER :: Nhydro    ! Number of hydrometeors
  INTEGER :: Niter     ! Number of calls to cosp_simulator
  INTEGER :: i_first,i_last ! First and last gridbox in each iteration
  INTEGER :: i,Ni
  INTEGER,DIMENSION(2) :: ix,iy
  LOGICAL :: reff_zero
  REAL :: maxp,minp
  INTEGER,DIMENSION(:),ALLOCATABLE :: & ! Dimensions nPoints
                  seed
!  It is recommended that the seed is set to a different value for each model
!  gridbox it is called on, as it is possible that the choice of the same 
!  seed value every time may introduce some statistical bias in the results, 
!  particularly for low values of NCOL.

  ! Types used in one iteration
  TYPE(cosp_gridbox) :: gbx_it
  TYPE(cosp_subgrid) :: sgx_it
  TYPE(cosp_vgrid)   :: vgrid_it
  TYPE(cosp_sgradar) :: sgradar_it
  TYPE(cosp_sglidar) :: sglidar_it
  TYPE(cosp_isccp)   :: isccp_it
  TYPE(cosp_modis)   :: modis_it
  TYPE(cosp_misr)    :: misr_it
  TYPE(cosp_rttov)   :: rttov_it
  TYPE(cosp_radarstats) :: stradar_it
  TYPE(cosp_lidarstats) :: stlidar_it
  ! Error handling
  INTEGER :: icode
  CHARACTER(LEN=200) :: cmessage
  CHARACTER(LEN=*),PARAMETER :: routine_name='COSP'


!++++++++++ Dimensions ++++++++++++
  Npoints  = gbx%Npoints
  Nlevels  = gbx%Nlevels
  Nhydro   = gbx%Nhydro

!++++++++++ Depth of model layers ++++++++++++
  DO i=1,Nlevels-1
    gbx%dlev(:,i) = gbx%zlev_half(:,i+1) - gbx%zlev_half(:,i)
  ENDDO
  gbx%dlev(:,Nlevels) = 2.0*(gbx%zlev(:,Nlevels) - gbx%zlev_half(:,Nlevels))

!++++++++++ Apply sanity checks to inputs ++++++++++
  CALL cosp_check_input('longitude',gbx%longitude,min_val=0.0,max_val=360.0)
  CALL cosp_check_input('latitude',gbx%latitude,min_val=-90.0,max_val=90.0)
  CALL cosp_check_input('dlev',gbx%dlev,min_val=0.0)
  CALL cosp_check_input('p',gbx%p,min_val=0.0)
  CALL cosp_check_input('ph',gbx%ph,min_val=0.0)
  CALL cosp_check_input('T',gbx%T,min_val=0.0)
  CALL cosp_check_input('q',gbx%q,min_val=0.0)
  CALL cosp_check_input('sh',gbx%sh,min_val=0.0)
  CALL cosp_check_input('dtau_s',gbx%dtau_s,min_val=0.0)
  CALL cosp_check_input('dtau_c',gbx%dtau_c,min_val=0.0)
  CALL cosp_check_input('dem_s',gbx%dem_s,min_val=0.0,max_val=1.0)
  CALL cosp_check_input('dem_c',gbx%dem_c,min_val=0.0,max_val=1.0)
  ! Point information (Npoints)
  CALL cosp_check_input('land',gbx%land,min_val=0.0,max_val=1.0)
  CALL cosp_check_input('psfc',gbx%psfc,min_val=0.0)
  CALL cosp_check_input('sunlit',gbx%sunlit,min_val=0.0,max_val=1.0)
  CALL cosp_check_input('skt',gbx%skt,min_val=0.0)
  ! TOTAL and CONV cloud fraction for SCOPS
  CALL cosp_check_input('tca',gbx%tca,min_val=0.0,max_val=1.0)
  CALL cosp_check_input('cca',gbx%cca,min_val=0.0,max_val=1.0)
  ! Precipitation fluxes on model levels
  CALL cosp_check_input('rain_ls',gbx%rain_ls,min_val=0.0)
  CALL cosp_check_input('rain_cv',gbx%rain_cv,min_val=0.0)
  CALL cosp_check_input('snow_ls',gbx%snow_ls,min_val=0.0)
  CALL cosp_check_input('snow_cv',gbx%snow_cv,min_val=0.0)
  CALL cosp_check_input('grpl_ls',gbx%grpl_ls,min_val=0.0)
  ! Hydrometeors concentration and distribution parameters
  CALL cosp_check_input('mr_hydro',gbx%mr_hydro,min_val=0.0)
  ! Effective radius [m]. (Npoints,Nlevels,Nhydro)
  CALL cosp_check_input('Reff',gbx%Reff,min_val=0.0)
  reff_zero=.TRUE.
  IF (ANY(gbx%Reff > 1.e-8)) THEN
     reff_zero=.FALSE.
      ! reff_zero == .false.
      !     and gbx%use_reff == .true.  Reff use in radar and lidar
      !     and reff_zero    == .false. Reff use in lidar and set to 0 for radar
  ENDIF
  IF ((.NOT. gbx%use_reff) .AND. (reff_zero)) THEN
        ! No Reff in radar. Default in lidar
        gbx%Reff = DEFAULT_LIDAR_REFF
        icode = -9
        cmessage = 'Using default Reff in lidar simulations'
        CALL cosp_ereport(routine_name,cmessage,icode)
  ENDIF

  ! Aerosols concentration and distribution parameters
  CALL cosp_check_input('conc_aero',gbx%conc_aero,min_val=0.0)
  ! Checks for CRM mode
  IF (Ncolumns == 1) THEN
     IF (gbx%use_precipitation_fluxes) THEN
        icode = 9
        cmessage = 'Use of precipitation fluxes not supported '//              &
                   'in CRM mode (Ncolumns=1)'
        CALL cosp_ereport(routine_name,cmessage,icode)
     ENDIF
     IF ((MAXVAL(gbx%dtau_c) > 0.0).OR.(MAXVAL(gbx%dem_c) > 0.0)) THEN
        icode = 9
        cmessage = 'dtau_c > 0.0 or dem_c > 0.0. In CRM mode (Ncolumns=1), '// &
                   'the optical depth (emmisivity) of all clouds must be '//   &
                   'passed through dtau_s (dem_s)'
        CALL cosp_ereport(routine_name,cmessage,icode)
     ENDIF
  ENDIF

   ! We base the seed in the decimal part of the surface pressure.
   ALLOCATE(seed(Npoints))
   seed = INT(gbx%psfc) ! This is to avoid division by zero when Npoints = 1
      ! Note: seed value of 0 caused me some problems + I want to
      ! randomize for each CALL to COSP even when Npoints ==1
   minp = MINVAL(gbx%psfc)
   maxp = MAXVAL(gbx%psfc)
   !if (Npoints .gt. 1) seed=int((gbx%psfc-minp)/(maxp-minp)*100000) + 1
   ! Below it's how it was done in the original implementation of the IS.
   ! The one above is better for offline data, when you may have packed data 
   ! that subsamples the decimal fraction of the surface pressure.
   IF (Npoints .gt. 1) seed=(gbx%psfc-INT(gbx%psfc))*1000000


   IF (gbx%Npoints_it >= gbx%Npoints) THEN ! One iteration gbx%Npoints
     CALL cosp_iter(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,      &
                    misr,modis,rttov,stradar,stlidar)
   ELSE ! Several iterations to save memory
        Niter = gbx%Npoints/gbx%Npoints_it ! Integer division
        IF (Niter*gbx%Npoints_it < gbx%Npoints) Niter = Niter + 1
        DO i=1,Niter
            i_first = (i-1)*gbx%Npoints_it + 1
            i_last  = i_first + gbx%Npoints_it - 1
            i_last  = MIN(i_last,gbx%Npoints)
            Ni = i_last - i_first + 1
            IF (i == 1) THEN
                ! Allocate types for all but last iteration
                CALL construct_cosp_gridbox(gbx%time,gbx%time_bnds,            &
                       gbx%radar_freq,gbx%surface_radar,gbx%use_mie_tables,    &
                       gbx%use_gas_abs,gbx%do_ray,gbx%melt_lay,gbx%k2,Ni,      &
                       Nlevels,Ncolumns,N_HYDRO,gbx%Nprmts_max_hydro,          &
                       gbx%Naero,gbx%Nprmts_max_aero,Ni,gbx%lidar_ice_type,    &
                       gbx%isccp_top_height,gbx%isccp_top_height_direction,    &
                       gbx%isccp_overlap,gbx%isccp_emsfc_lw,                   &
                       gbx%use_precipitation_fluxes,gbx%use_reff,              &
                       gbx%plat,gbx%sat,gbx%inst,gbx%nchan,gbx%ZenAng,         &
                       gbx%Ichan(1:gbx%nchan),gbx%surfem(1:gbx%nchan),         &
                       gbx%co2,gbx%ch4,gbx%n2o,gbx%co,gbx_it)
                CALL construct_cosp_vgrid(gbx_it,vgrid%Nlvgrid,vgrid%use_vgrid,&
                       vgrid%csat_vgrid,vgrid_it)
                CALL construct_cosp_subgrid(Ni, Ncolumns, Nlevels, sgx_it)
                CALL construct_cosp_sgradar(cfg,Ni,Ncolumns,Nlevels,N_HYDRO,   &
                       sgradar_it)
                CALL construct_cosp_sglidar(cfg,Ni,Ncolumns,Nlevels,N_HYDRO,   &
                       PARASOL_NREFL,sglidar_it)
                CALL construct_cosp_isccp(cfg,Ni,Ncolumns,Nlevels,isccp_it)
                CALL construct_cosp_modis(cfg, Ni, modis_it)
                CALL construct_cosp_misr(cfg,Ni,misr_it)
                CALL construct_cosp_rttov(cfg,Ni,gbx%nchan,rttov_it)
                CALL construct_cosp_radarstats(cfg,Ni,Ncolumns,vgrid%Nlvgrid,  &
                       N_HYDRO,stradar_it)
                CALL construct_cosp_lidarstats(cfg,Ni,Ncolumns,vgrid%Nlvgrid,  &
                       N_HYDRO,PARASOL_NREFL,stlidar_it)
            ELSEIF (i == Niter) THEN ! last iteration
                CALL free_cosp_gridbox(gbx_it,.TRUE.)
                CALL free_cosp_subgrid(sgx_it)
                CALL free_cosp_vgrid(vgrid_it)
                CALL free_cosp_sgradar(sgradar_it)
                CALL free_cosp_sglidar(sglidar_it)
                CALL free_cosp_isccp(isccp_it)
                CALL free_cosp_modis(modis_it)
                CALL free_cosp_misr(misr_it)
                CALL free_cosp_rttov(rttov_it)
                CALL free_cosp_radarstats(stradar_it)
                CALL free_cosp_lidarstats(stlidar_it)
                ! Allocate types for iterations
                CALL construct_cosp_gridbox(gbx%time,gbx%time_bnds,            &
                       gbx%radar_freq,gbx%surface_radar,gbx%use_mie_tables,    &
                       gbx%use_gas_abs,gbx%do_ray,gbx%melt_lay,gbx%k2,Ni,      &
                       Nlevels,Ncolumns,N_HYDRO,gbx%Nprmts_max_hydro,          &
                       gbx%Naero,gbx%Nprmts_max_aero,Ni,gbx%lidar_ice_type,    &
                       gbx%isccp_top_height,gbx%isccp_top_height_direction,    &
                       gbx%isccp_overlap,gbx%isccp_emsfc_lw,                   &
                       gbx%use_precipitation_fluxes,gbx%use_reff,              &
                       gbx%plat,gbx%sat,gbx%inst,gbx%nchan,gbx%ZenAng,         &
                       gbx%Ichan(1:gbx%nchan),gbx%surfem(1:gbx%nchan),         &
                       gbx%co2,gbx%ch4,gbx%n2o,gbx%co,gbx_it)
                ! --- Copy arrays without Npoints as dimension ---
                gbx_it%dist_prmts_hydro = gbx%dist_prmts_hydro
                gbx_it%dist_type_aero   = gbx_it%dist_type_aero
                CALL construct_cosp_vgrid(gbx_it,vgrid%Nlvgrid,vgrid%use_vgrid,&
                       vgrid%csat_vgrid,vgrid_it)
                CALL construct_cosp_subgrid(Ni, Ncolumns, Nlevels, sgx_it)
                CALL construct_cosp_sgradar(cfg,Ni,Ncolumns,Nlevels,N_HYDRO,   &
                       sgradar_it)
                CALL construct_cosp_sglidar(cfg,Ni,Ncolumns,Nlevels,N_HYDRO,   &
                       PARASOL_NREFL,sglidar_it)
                CALL construct_cosp_isccp(cfg,Ni,Ncolumns,Nlevels,isccp_it)
                CALL construct_cosp_modis(cfg,Ni, modis_it)
                CALL construct_cosp_misr(cfg,Ni,misr_it)
                CALL construct_cosp_rttov(cfg,Ni,gbx%nchan,rttov_it)
                CALL construct_cosp_radarstats(cfg,Ni,Ncolumns,vgrid%Nlvgrid,  &
                       N_HYDRO,stradar_it)
                CALL construct_cosp_lidarstats(cfg,Ni,Ncolumns,vgrid%Nlvgrid,  &
                       N_HYDRO,PARASOL_NREFL,stlidar_it)
            ENDIF
            ! --- Copy sections of arrays with Npoints as dimension ---
            ix=(/i_first,i_last/)
            iy=(/1,Ni/)
            CALL cosp_gridbox_cpsection(ix,iy,gbx,gbx_it)
              ! These serve as initialisation of *_it types
            CALL cosp_subgrid_cpsection(ix,iy,sgx,sgx_it)
            IF (cfg%Lradar_sim)                                                &
                CALL cosp_sgradar_cpsection(ix,iy,sgradar,sgradar_it)
            IF (cfg%Llidar_sim)                                                &
                CALL cosp_sglidar_cpsection(ix,iy,sglidar,sglidar_it)
            IF (cfg%Lisccp_sim) CALL cosp_isccp_cpsection(ix,iy,isccp,isccp_it)
            IF (cfg%Lmodis_sim) CALL cosp_modis_cpsection(ix,iy,modis,modis_it)
            IF (cfg%Lmisr_sim)  CALL cosp_misr_cpsection(ix,iy,misr,misr_it)
            IF (cfg%Lrttov_sim) CALL cosp_rttov_cpsection(ix,iy,rttov,rttov_it)
            IF (cfg%Lradar_sim)                                                &
                CALL cosp_radarstats_cpsection(ix,iy,stradar,stradar_it)
            IF (cfg%Llidar_sim)                                                &
                CALL cosp_lidarstats_cpsection(ix,iy,stlidar,stlidar_it)

            CALL cosp_iter(overlap,seed(ix(1):ix(2)),cfg,vgrid_it,gbx_it,      &
                   sgx_it,sgradar_it,sglidar_it,isccp_it,misr_it,modis_it,     &
                   rttov_it,stradar_it,stlidar_it)
            ! --- Copy results to output structures ---
            ix=(/1,Ni/)
            iy=(/i_first,i_last/)
            CALL cosp_subgrid_cpsection(ix,iy,sgx_it,sgx)
            IF (cfg%Lradar_sim)                                                &
                CALL cosp_sgradar_cpsection(ix,iy,sgradar_it,sgradar)
            IF (cfg%Llidar_sim)                                                &
                CALL cosp_sglidar_cpsection(ix,iy,sglidar_it,sglidar)
            IF (cfg%Lisccp_sim) CALL cosp_isccp_cpsection(ix,iy,isccp_it,isccp)
            IF (cfg%Lmodis_sim) CALL cosp_modis_cpsection(ix,iy,modis_it,modis)
            IF (cfg%Lmisr_sim)  CALL cosp_misr_cpsection(ix,iy,misr_it,misr)
            IF (cfg%Lrttov_sim) CALL cosp_rttov_cpsection(ix,iy,rttov_it,rttov)
            IF (cfg%Lradar_sim)                                                &
                CALL cosp_radarstats_cpsection(ix,iy,stradar_it,stradar)
            IF (cfg%Llidar_sim)                                                &
                CALL cosp_lidarstats_cpsection(ix,iy,stlidar_it,stlidar)
        ENDDO
        ! Deallocate types
        CALL free_cosp_gridbox(gbx_it,.TRUE.)
        CALL free_cosp_subgrid(sgx_it)
        CALL free_cosp_vgrid(vgrid_it)
        CALL free_cosp_sgradar(sgradar_it)
        CALL free_cosp_sglidar(sglidar_it)
        CALL free_cosp_isccp(isccp_it)
        CALL free_cosp_modis(modis_it)
        CALL free_cosp_misr(misr_it)
        CALL free_cosp_rttov(rttov_it)
        CALL free_cosp_radarstats(stradar_it)
        CALL free_cosp_lidarstats(stlidar_it)
   ENDIF
   DEALLOCATE(seed)

END SUBROUTINE COSP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP_ITER ----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_ITER(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,     &
                     misr,modis,rttov,stradar,stlidar)

  ! Arguments
  INTEGER,INTENT(IN) :: overlap ! overlap type in SCOPS: 1=max,2=rand,3=max/rand
  INTEGER,DIMENSION(:),INTENT(IN) :: seed
  TYPE(cosp_config),INTENT(IN) :: cfg   ! Configuration options
  TYPE(cosp_vgrid),INTENT(IN) :: vgrid   ! Information on vertical grid of stats
  TYPE(cosp_gridbox),INTENT(INOUT) :: gbx
  TYPE(cosp_subgrid),INTENT(INOUT) :: sgx   ! Subgrid info
  TYPE(cosp_sgradar),INTENT(INOUT) :: sgradar ! Output from radar simulator
  TYPE(cosp_sglidar),INTENT(INOUT) :: sglidar ! Output from lidar simulator
  TYPE(cosp_isccp),INTENT(INOUT)   :: isccp   ! Output from ISCCP simulator
  TYPE(cosp_misr),INTENT(INOUT)    :: misr    ! Output from MISR simulator
  TYPE(cosp_modis),INTENT(INOUT)   :: modis   ! Output from MODIS simulator
  TYPE(cosp_rttov),INTENT(INOUT)   :: rttov   ! Output from RTTOV
  TYPE(cosp_radarstats),INTENT(INOUT) :: stradar ! Summary statistics from radar
  TYPE(cosp_lidarstats),INTENT(INOUT) :: stlidar ! Summary statistics from lidar

  ! Local variables 
  INTEGER :: Npoints   ! Number of gridpoints
  INTEGER :: Ncolumns  ! Number of subcolumns
  INTEGER :: Nlevels   ! Number of levels
  INTEGER :: Nhydro    ! Number of hydrometeors
  INTEGER :: i,j,k
  INTEGER :: I_HYDRO 
  REAL,DIMENSION(:,:),POINTER :: column_frac_out ! One column of frac_out
  REAL,DIMENSION(:,:),POINTER :: column_prec_out ! One column of prec_frac
  INTEGER :: scops_debug=0    ! non-zero to print out debug info in SCOPS
  REAL,DIMENSION(:, :),ALLOCATABLE :: cca_scops,ls_p_rate,cv_p_rate, &
                     tca_scops ! Cloud cover in each model level 
                               ! (HORIZONTAL gridbox fraction) of total cloud.
                               ! Levels are from TOA to SURFACE. (nPoints, nLev)
  REAL,DIMENSION(:,:),ALLOCATABLE :: frac_ls,prec_ls,frac_cv,prec_cv 
                              ! Cloud/Precipitation fraction in each model level
                              ! Levels are from SURFACE to TOA
  REAL,DIMENSION(:,:),ALLOCATABLE :: rho ! (Npoints, Nlevels). Air density
  TYPE(cosp_sghydro) :: sghydro   ! Subgrid info for hydrometeors


  !++++++++++ Dimensions ++++++++++++
  Npoints  = gbx%Npoints
  Ncolumns = gbx%Ncolumns
  Nlevels  = gbx%Nlevels
  Nhydro   = gbx%Nhydro

  !++++++++++ Climate/NWP mode ++++++++++  
  IF (Ncolumns > 1) THEN
        !++++++++++ Subgrid sampling ++++++++++
        ! Allocate arrays before calling SCOPS
        ALLOCATE(frac_ls(Npoints,Nlevels),frac_cv(Npoints,Nlevels),            &
                 prec_ls(Npoints,Nlevels),prec_cv(Npoints,Nlevels))
        ALLOCATE(tca_scops(Npoints,Nlevels),cca_scops(Npoints,Nlevels), &
                ls_p_rate(Npoints,Nlevels),cv_p_rate(Npoints,Nlevels))
        ! Initialize to zero
        frac_ls=0.0
        prec_ls=0.0
        frac_cv=0.0
        prec_cv=0.0
        ! Cloud fractions for SCOPS from TOA to SFC
        tca_scops = gbx%tca(:,Nlevels:1:-1)
        cca_scops = gbx%cca(:,Nlevels:1:-1)

        ! Call to SCOPS
        ! strat and conv arrays are passed with levels from TOA to SURFACE.
! DEPENDS ON: scops
        CALL scops(Npoints,Nlevels,Ncolumns,seed,tca_scops,cca_scops,overlap,  &
                   sgx%frac_out,scops_debug)

        ! temporarily use prec_ls/cv to transfer information about 
        ! precipitation flux into prec_scops
        IF (gbx%use_precipitation_fluxes) THEN
            ls_p_rate(:,Nlevels:1:-1)=gbx%rain_ls(:,1:Nlevels) +               &
                                      gbx%snow_ls(:,1:Nlevels) +               &
                                      gbx%grpl_ls(:,1:Nlevels)
            cv_p_rate(:,Nlevels:1:-1)=gbx%rain_cv(:,1:Nlevels) +               &
                                      gbx%snow_cv(:,1:Nlevels)
        ELSE
            ls_p_rate(:,Nlevels:1:-1)=gbx%mr_hydro(:,1:Nlevels,I_LSRAIN) + &
                                      gbx%mr_hydro(:,1:Nlevels,I_LSSNOW) + &
                                      gbx%mr_hydro(:,1:Nlevels,I_LSGRPL)
            cv_p_rate(:,Nlevels:1:-1)=gbx%mr_hydro(:,1:Nlevels,I_CVRAIN) + &
                                      gbx%mr_hydro(:,1:Nlevels,I_CVSNOW)
        ENDIF

! DEPENDS ON: prec_scops
        CALL prec_scops(Npoints,Nlevels,Ncolumns,ls_p_rate,cv_p_rate,          &
                        sgx%frac_out,sgx%prec_frac)

        ! Precipitation fraction
        DO j=1,Npoints,1
          DO k=1,Nlevels,1
            DO i=1,Ncolumns,1
                IF (sgx%frac_out (j,i,Nlevels+1-k) == I_LSC)                   &
                    frac_ls(j,k)=frac_ls(j,k)+1.
                IF (sgx%frac_out (j,i,Nlevels+1-k) == I_CVC)                   &
                    frac_cv(j,k)=frac_cv(j,k)+1.
                IF (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 1)                     &
                    prec_ls(j,k)=prec_ls(j,k)+1.
                IF (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 2)                     &
                    prec_cv(j,k)=prec_cv(j,k)+1.
                IF (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 3) THEN
                    prec_cv(j,k)=prec_cv(j,k)+1.
                    prec_ls(j,k)=prec_ls(j,k)+1.
                ENDIF
            ENDDO  !i
            frac_ls(j,k)=frac_ls(j,k)/Ncolumns
            frac_cv(j,k)=frac_cv(j,k)/Ncolumns
            prec_ls(j,k)=prec_ls(j,k)/Ncolumns
            prec_cv(j,k)=prec_cv(j,k)/Ncolumns
          ENDDO  !k
        ENDDO  !j

         ! Levels from SURFACE to TOA.
        IF (Npoints*Ncolumns*Nlevels < 10000) THEN
            sgx%frac_out(1:Npoints,:,1:Nlevels) =                              &
                   sgx%frac_out(1:Npoints,:,Nlevels:1:-1)
            sgx%prec_frac(1:Npoints,:,1:Nlevels) =                             &
                   sgx%prec_frac(1:Npoints,:,Nlevels:1:-1)
        ELSE
            ! This is done within a loop over nPoints to save memory
            DO j=1,Npoints
                sgx%frac_out(j,:,1:Nlevels)  = sgx%frac_out(j,:,Nlevels:1:-1)
                sgx%prec_frac(j,:,1:Nlevels) = sgx%prec_frac(j,:,Nlevels:1:-1)
            ENDDO
        ENDIF

       ! Deallocate arrays that will no longer be used
        DEALLOCATE(tca_scops,cca_scops,ls_p_rate,cv_p_rate)

        ! Populate the subgrid arrays
        CALL construct_cosp_sghydro(Npoints,Ncolumns,Nlevels,Nhydro,sghydro)
        DO k=1,Ncolumns
            !--------- Mixing ratios for clouds and Reff for Clouds and precip -
            column_frac_out => sgx%frac_out(:,k,:)
            WHERE (column_frac_out == I_LSC)     !+++++++++++ LS clouds ++++++++
                sghydro%mr_hydro(:,k,:,I_LSCLIQ) = gbx%mr_hydro(:,:,I_LSCLIQ)
                sghydro%mr_hydro(:,k,:,I_LSCICE) = gbx%mr_hydro(:,:,I_LSCICE)

                sghydro%Reff(:,k,:,I_LSCLIQ)     = gbx%Reff(:,:,I_LSCLIQ)
                sghydro%Reff(:,k,:,I_LSCICE)     = gbx%Reff(:,:,I_LSCICE)

                sghydro%Np(:,k,:,I_LSCLIQ)     = gbx%Np(:,:,I_LSCLIQ)
                sghydro%Np(:,k,:,I_LSCICE)     = gbx%Np(:,:,I_LSCICE)

            ELSEWHERE (column_frac_out == I_CVC) !+++++++++++ CONV clouds ++++++
                sghydro%mr_hydro(:,k,:,I_CVCLIQ) = gbx%mr_hydro(:,:,I_CVCLIQ)
                sghydro%mr_hydro(:,k,:,I_CVCICE) = gbx%mr_hydro(:,:,I_CVCICE)

                sghydro%Reff(:,k,:,I_CVCLIQ)     = gbx%Reff(:,:,I_CVCLIQ)
                sghydro%Reff(:,k,:,I_CVCICE)     = gbx%Reff(:,:,I_CVCICE)

                sghydro%Np(:,k,:,I_CVCLIQ)     = gbx%Np(:,:,I_CVCLIQ)
                sghydro%Np(:,k,:,I_CVCICE)     = gbx%Np(:,:,I_CVCICE)

            END WHERE
            column_prec_out => sgx%prec_frac(:,k,:)
            WHERE ((column_prec_out == 1) .OR. (column_prec_out == 3) )
                !++++ LS precip ++++
                sghydro%Reff(:,k,:,I_LSRAIN) = gbx%Reff(:,:,I_LSRAIN)
                sghydro%Reff(:,k,:,I_LSSNOW) = gbx%Reff(:,:,I_LSSNOW)
                sghydro%Reff(:,k,:,I_LSGRPL) = gbx%Reff(:,:,I_LSGRPL)

                sghydro%Np(:,k,:,I_LSRAIN)     = gbx%Np(:,:,I_LSRAIN)
                sghydro%Np(:,k,:,I_LSSNOW)     = gbx%Np(:,:,I_LSSNOW)
                sghydro%Np(:,k,:,I_LSGRPL)     = gbx%Np(:,:,I_LSGRPL)
            ELSEWHERE ((column_prec_out == 2) .OR. (column_prec_out == 3))
                !++++ CONV precip ++++
                sghydro%Reff(:,k,:,I_CVRAIN) = gbx%Reff(:,:,I_CVRAIN)
                sghydro%Reff(:,k,:,I_CVSNOW) = gbx%Reff(:,:,I_CVSNOW)

                sghydro%Np(:,k,:,I_CVRAIN)     = gbx%Np(:,:,I_CVRAIN)
                sghydro%Np(:,k,:,I_CVSNOW)     = gbx%Np(:,:,I_CVSNOW)
            END WHERE
            !--------- Precip -------
            IF (.NOT. gbx%use_precipitation_fluxes) THEN
                WHERE (column_frac_out == I_LSC)
                    !+++++++++++ LS Precipitation ++++++++
                    sghydro%mr_hydro(:,k,:,I_LSRAIN)= gbx%mr_hydro(:,:,I_LSRAIN)
                    sghydro%mr_hydro(:,k,:,I_LSSNOW)= gbx%mr_hydro(:,:,I_LSSNOW)
                    sghydro%mr_hydro(:,k,:,I_LSGRPL)= gbx%mr_hydro(:,:,I_LSGRPL)
                ELSEWHERE (column_frac_out == I_CVC)
                    !+++++++++++ CONV Precipitation ++++++++
                    sghydro%mr_hydro(:,k,:,I_CVRAIN)= gbx%mr_hydro(:,:,I_CVRAIN)
                    sghydro%mr_hydro(:,k,:,I_CVSNOW)= gbx%mr_hydro(:,:,I_CVSNOW)
                END WHERE
            ENDIF
        ENDDO
        ! convert the mixing ratio and precipitation flux from gridbox mean
        ! to the fraction-based values
        DO k=1,Nlevels
            DO j=1,Npoints
                !--------- Clouds -------
                IF (frac_ls(j,k) .ne. 0.) THEN
                    sghydro%mr_hydro(j,:,k,I_LSCLIQ) =                         &
                           sghydro%mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                    sghydro%mr_hydro(j,:,k,I_LSCICE) =                         &
                           sghydro%mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
                ENDIF
                IF (frac_cv(j,k) .ne. 0.) THEN
                    sghydro%mr_hydro(j,:,k,I_CVCLIQ) =                         &
                           sghydro%mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                    sghydro%mr_hydro(j,:,k,I_CVCICE) =                         &
                           sghydro%mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
                ENDIF
                !--------- Precip -------
                IF (gbx%use_precipitation_fluxes) THEN
                    IF (prec_ls(j,k) .ne. 0.) THEN
                        gbx%rain_ls(j,k) = gbx%rain_ls(j,k)/prec_ls(j,k)
                        gbx%snow_ls(j,k) = gbx%snow_ls(j,k)/prec_ls(j,k)
                        gbx%grpl_ls(j,k) = gbx%grpl_ls(j,k)/prec_ls(j,k)
                    ENDIF
                    IF (prec_cv(j,k) .ne. 0.) THEN
                        gbx%rain_cv(j,k) = gbx%rain_cv(j,k)/prec_cv(j,k)
                        gbx%snow_cv(j,k) = gbx%snow_cv(j,k)/prec_cv(j,k)
                    ENDIF
                ELSE
                    IF (prec_ls(j,k) .ne. 0.) THEN
                        sghydro%mr_hydro(j,:,k,I_LSRAIN) =                     &
                               sghydro%mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                        sghydro%mr_hydro(j,:,k,I_LSSNOW) =                     &
                               sghydro%mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
                        sghydro%mr_hydro(j,:,k,I_LSGRPL) =                     &
                               sghydro%mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
                    ENDIF
                    IF (prec_cv(j,k) .ne. 0.) THEN
                        sghydro%mr_hydro(j,:,k,I_CVRAIN) =                     &
                               sghydro%mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                        sghydro%mr_hydro(j,:,k,I_CVSNOW) =                     &
                               sghydro%mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
                    ENDIF
                ENDIF
            ENDDO !k
        ENDDO !j
        DEALLOCATE(frac_ls,prec_ls,frac_cv,prec_cv)

        IF (gbx%use_precipitation_fluxes) THEN
          ! Density
          ALLOCATE(rho(Npoints,Nlevels))
          I_HYDRO = I_LSRAIN
          CALL cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,       &
                  sgx%prec_frac,1.,n_ax(I_HYDRO),n_bx(I_HYDRO),                &
                  alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO),g_x(I_HYDRO),     &
                  a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                  gamma_3(I_HYDRO),gamma_4(I_HYDRO),gbx%rain_ls,               &
                  sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
          I_HYDRO = I_LSSNOW
          CALL cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,       &
                  sgx%prec_frac,1.,n_ax(I_HYDRO),n_bx(I_HYDRO),                &
                  alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO),g_x(I_HYDRO),     &
                  a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                  gamma_3(I_HYDRO),gamma_4(I_HYDRO),gbx%snow_ls,               &
                  sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
          I_HYDRO = I_CVRAIN
          CALL cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,       &
                  sgx%prec_frac,2.,n_ax(I_HYDRO),n_bx(I_HYDRO),                &
                  alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO),g_x(I_HYDRO),     &
                  a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                  gamma_3(I_HYDRO),gamma_4(I_HYDRO),gbx%rain_cv,               &
                  sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
          I_HYDRO = I_CVSNOW
          CALL cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,       &
                  sgx%prec_frac,2.,n_ax(I_HYDRO),n_bx(I_HYDRO),                &
                  alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO),g_x(I_HYDRO),     &
                  a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                  gamma_3(I_HYDRO),gamma_4(I_HYDRO),gbx%snow_cv,               &
                  sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
          I_HYDRO = I_LSGRPL
          CALL cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,       &
                  sgx%prec_frac,1.,n_ax(I_HYDRO),n_bx(I_HYDRO),                &
                  alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO),g_x(I_HYDRO),     &
                  a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                  gamma_3(I_HYDRO),gamma_4(I_HYDRO),gbx%grpl_ls,               &
                  sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
          IF(ALLOCATED(rho)) DEALLOCATE(rho)
        ENDIF
   !++++++++++ CRM mode ++++++++++
   ELSE
      CALL construct_cosp_sghydro(Npoints,Ncolumns,Nlevels,Nhydro,sghydro)
      sghydro%mr_hydro(:,1,:,:) = gbx%mr_hydro
      sghydro%Reff(:,1,:,:) = gbx%Reff
      sghydro%Np(:,1,:,:) = gbx%Np      ! added by Roj with Quickbeam V3.0

      !--------- Clouds -------
      WHERE ((gbx%dtau_s > 0.0))
        ! Subgrid cloud array. Dimensions (Npoints,Ncolumns,Nlevels)
        sgx%frac_out(:,1,:) = 1
      ENDWHERE
   ENDIF ! Ncolumns > 1

   !++++++++++ Simulator ++++++++++
    CALL cosp_simulator(gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,       &
                        misr,modis,rttov,stradar,stlidar)

    ! Deallocate subgrid arrays
    CALL free_cosp_sghydro(sghydro)
END SUBROUTINE COSP_ITER

END MODULE MOD_COSP
