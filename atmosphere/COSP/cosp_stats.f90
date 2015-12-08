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

MODULE MOD_COSP_STATS
  USE cosp_types_mod
  USE MOD_LLNL_STATS
  USE MOD_LMD_IPSL_STATS
  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_STATS ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_STATS(gbx,sgx,cfg,sgradar,sglidar,vgrid,stradar,stlidar)
   IMPLICIT NONE
   ! Input arguments
   TYPE(cosp_gridbox),INTENT(IN) :: gbx
   TYPE(cosp_subgrid),INTENT(IN) :: sgx
   TYPE(cosp_config),INTENT(IN)  :: cfg
   TYPE(cosp_sgradar),INTENT(IN) :: sgradar
   TYPE(cosp_sglidar),INTENT(IN) :: sglidar
   TYPE(cosp_vgrid),INTENT(IN)   :: vgrid
   ! Output arguments
   TYPE(cosp_radarstats),INTENT(INOUT) :: stradar ! Summary statistics for radar
   TYPE(cosp_lidarstats),INTENT(INOUT) :: stlidar ! Summary statistics for lidar

   ! Local variables
   INTEGER :: Npoints  !# of grid points
   INTEGER :: Nlevels  !# of levels
   INTEGER :: Nhydro   !# of hydrometeors
   INTEGER :: Ncolumns !# of columns
   INTEGER :: Nlr
   LOGICAL :: ok_lidar_cfad = .FALSE.
   REAL,DIMENSION(:,:,:),ALLOCATABLE :: Ze_out,betatot_out,betamol_in,         &
                                        betamol_out,ph_in,ph_out,x3d
   REAL,DIMENSION(:,:),ALLOCATABLE :: ph_c,betamol_c

   Npoints  = gbx%Npoints
   Nlevels  = gbx%Nlevels
   Nhydro   = gbx%Nhydro
   Ncolumns = gbx%Ncolumns
   Nlr      = vgrid%Nlvgrid

   IF (cfg%LcfadLidarsr532) ok_lidar_cfad=.TRUE.

   IF (vgrid%use_vgrid) THEN ! Statistics in a different vertical grid
      ALLOCATE(Ze_out(Npoints,Ncolumns,Nlr),betatot_out(Npoints,Ncolumns,Nlr), &
                betamol_in(Npoints,1,Nlevels),betamol_out(Npoints,1,Nlr),      &
                betamol_c(Npoints,Nlr),ph_in(Npoints,1,Nlevels),               &
                ph_out(Npoints,1,Nlr),ph_c(Npoints,Nlr))
      Ze_out = 0.0
      betatot_out  = 0.0
      betamol_out= 0.0
      betamol_c  = 0.0
      ph_in(:,1,:)  = gbx%ph(:,:)
      ph_out  = 0.0
      ph_c    = 0.0
      !++++++++++++ Radar CFAD ++++++++++++++++
      IF (cfg%Lradar_sim) THEN
        IF (cfg%LcfadDbze94.OR.cfg%Ldbze94gbx)                                 &
            CALL cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,           &
                    gbx%zlev,gbx%zlev_half,sgradar%Ze_tot,Nlr,vgrid%zl,        &
                    vgrid%zu,Ze_out,log_units=.TRUE.)
        IF (cfg%LcfadDbze94)                                                   &
            stradar%cfad_ze = cosp_cfad(Npoints,Ncolumns,Nlr,DBZE_BINS,Ze_out, &
                                 DBZE_MIN,DBZE_MAX,CFAD_ZE_MIN,CFAD_ZE_WIDTH)
        IF (cfg%Ldbze94gbx)                                                    &
            CALL cosp_gridbox_mean(Npoints,Ncolumns,Nlr,gbx%zlev_half(:,1),    &
                    vgrid%zu,.TRUE.,Ze_out,stradar%dbze94gbx,log_units=.TRUE.)
      ENDIF
      !++++++++++++ Lidar-only diagnostics ++++++++++++++++
      IF (cfg%Llidar_sim) THEN
        betamol_in(:,1,:) = sglidar%beta_mol(:,:)
        CALL cosp_change_vertical_grid(Npoints,1,Nlevels,gbx%zlev,             &
                gbx%zlev_half,betamol_in,Nlr,vgrid%zl,vgrid%zu,betamol_out)
        CALL cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,gbx%zlev,      &
                gbx%zlev_half,sglidar%beta_tot,Nlr,vgrid%zl,                   &
                vgrid%zu,betatot_out)
        CALL cosp_change_vertical_grid(Npoints,1,Nlevels,gbx%zlev,             &
                gbx%zlev_half,ph_in,Nlr,vgrid%zl,vgrid%zu,ph_out)
        ph_c(:,:) = ph_out(:,1,:)
        betamol_c(:,:) = betamol_out(:,1,:)
        IF (cfg%Latb532gbx)                                                    &
            CALL cosp_gridbox_mean(Npoints,Ncolumns,Nlr,gbx%zlev_half(:,1),    &
                    vgrid%zu,.TRUE.,betatot_out,stlidar%atb532gbx)
        ! Stats from lidar_stat_summary
        CALL diag_lidar(Npoints,Ncolumns,Nlr,SR_BINS,PARASOL_NREFL,            &
                betatot_out,betamol_c,sglidar%refl,gbx%land,ph_c,LIDAR_UNDEF,  &
                ok_lidar_cfad,stlidar%cfad_sr,stlidar%srbval,LIDAR_NCAT,       &
                stlidar%lidarcld,stlidar%cldlayer,stlidar%parasolrefl)
      ENDIF

      !+++++ Lidar-only cloud amount and lidar&radar total cloud mount ++++++++
      IF (cfg%Lradar_sim.AND.cfg%Llidar_sim)                                   &
          CALL cosp_lidar_only_cloud(Npoints,Ncolumns,Nlr,betatot_out,         &
                  betamol_c,Ze_out,stradar%lidar_only_freq_cloud,              &
                  stradar%radar_lidar_tcc)
      ! Deallocate arrays at coarse resolution
      DEALLOCATE(Ze_out,betatot_out,betamol_in,betamol_out,betamol_c,ph_in,    &
                 ph_out,ph_c)
   ELSE ! Statistics in model levels
      !++++++++++++ Radar-only diagnostics ++++++++++++++++
      IF (cfg%Lradar_sim) THEN
          IF (cfg%LcfadDbze94)                                                 &
              stradar%cfad_ze = cosp_cfad(Npoints,Ncolumns,Nlr,DBZE_BINS,      &
                                   sgradar%Ze_tot,DBZE_MIN,DBZE_MAX,           &
                                   CFAD_ZE_MIN,CFAD_ZE_WIDTH)
          IF (cfg%Ldbze94gbx)                                                  &
              CALL cosp_gridbox_mean(Npoints,Ncolumns,Nlevels,                 &
                      gbx%zlev_half(:,1),vgrid%zu,.FALSE.,sgradar%Ze_tot,      &
                      stradar%dbze94gbx,log_units=.TRUE.)
      ENDIF
      !++++++++++++ Lidar CFAD ++++++++++++++++
      ! Stats from lidar_stat_summary
      IF (cfg%Llidar_sim) THEN
          CALL diag_lidar(Npoints,Ncolumns,Nlr,SR_BINS,PARASOL_NREFL,          &
                  sglidar%beta_tot,sglidar%beta_mol,sglidar%refl,gbx%land,     &
                  gbx%ph,LIDAR_UNDEF,ok_lidar_cfad,stlidar%cfad_sr,            &
                  stlidar%srbval,LIDAR_NCAT,stlidar%lidarcld,stlidar%cldlayer, &
                  stlidar%parasolrefl)
          IF (cfg%Latb532gbx)                                                  &
              CALL cosp_gridbox_mean(Npoints,Ncolumns,Nlevels,                 &
                      gbx%zlev_half(:,1),vgrid%zu,.FALSE.,sglidar%beta_tot,    &
                      stlidar%atb532gbx)

      ENDIF
      !++++++++++++ Lidar-only cloud amount and lidar&radar total cloud mount ++++++++++++++++
      IF (cfg%Lradar_sim.AND.cfg%Llidar_sim)                                   &
          CALL cosp_lidar_only_cloud(Npoints,Ncolumns,Nlr,sglidar%beta_tot,    &
                  sglidar%beta_mol,sgradar%Ze_tot,                             &
                  stradar%lidar_only_freq_cloud,stradar%radar_lidar_tcc)
   ENDIF
   ! Replace undef
   WHERE (stlidar%cfad_sr   == LIDAR_UNDEF) stlidar%cfad_sr   = R_UNDEF
   WHERE (stlidar%lidarcld  == LIDAR_UNDEF) stlidar%lidarcld  = R_UNDEF
   WHERE (stlidar%cldlayer  == LIDAR_UNDEF) stlidar%cldlayer  = R_UNDEF
   WHERE (stlidar%parasolrefl == LIDAR_UNDEF) stlidar%parasolrefl = R_UNDEF

END SUBROUTINE COSP_STATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------- SUBROUTINE COSP_GRIDBOX_MEAN ----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_GRIDBOX_MEAN(Npoints,Ncolumns,Nlevels,zhalf,zu,lstdgrid,x,r,   &
                             log_units)
   IMPLICIT NONE
   ! Input arguments
   INTEGER,INTENT(IN) :: Npoints  !# of grid points
   INTEGER,INTENT(IN) :: Nlevels  !# of grid levels
   INTEGER,INTENT(IN) :: Ncolumns !# of columns
   REAL,DIMENSION(Npoints),INTENT(IN) :: zhalf ! Height at bottom half level [m]
   REAL,DIMENSION(Nlevels),INTENT(IN) :: zu ! Upper boundary of new levels  [m]
   LOGICAL,INTENT(IN) :: lstdgrid ! Using standard grid (40 cloudsat levels)?
   REAL,DIMENSION(Npoints,Ncolumns,Nlevels),INTENT(IN) :: x
                                           !x -  Input variable to be averaged
   LOGICAL,OPTIONAL,INTENT(IN) :: log_units ! Need to convert to linear units
   ! Output
   REAL,DIMENSION(Npoints,Nlevels),INTENT(OUT) :: r
                                           ! Variable averaged over subcolumns

   ! Local variables
   INTEGER :: i,j,k,N
   LOGICAL :: lunits
   REAL :: t

   lunits=.FALSE.
   IF (PRESENT(log_units)) lunits=log_units

   r = 0.0

   DO k=1,Nlevels
     DO i=1,Npoints
       IF (lstdgrid .AND. (zu(k) <= zhalf(i))) THEN
          ! Level of standard grid below model bottom level
          r(i,k) = R_UNDEF
       ELSE
          N = 0 ! Counter of valid points
          IF (lunits) THEN
            DO j=1,Ncolumns
                t = x(i,j,k)
                IF (t /= R_UNDEF) THEN
                  t = 10.0**(t/10.0)
                ELSE
                  t = 0.0
                ENDIF
                r(i,k) = r(i,k) + t
                N = N+1
            ENDDO
            ! Divide and change units back
            IF (N /=0 ) r(i,k) = r(i,k)/N
            IF (r(i,k) <= 0.0) THEN
                r(i,k) = R_UNDEF
            ELSE
                r(i,k) = 10.0*LOG10(r(i,k))
            ENDIF
          ELSE ! lunits=.false
            DO j=1,Ncolumns
                t = x(i,j,k)
                IF (t /= R_UNDEF) THEN
                   r(i,k) = r(i,k) + t
                   N = N+1
                ENDIF
            ENDDO
            ! Compute mean
            IF (N /=0 ) THEN
               r(i,k) = r(i,k)/N
            ELSE
               r(i,k) = R_UNDEF
            ENDIF
          ENDIF ! LUNITS
       ENDIF
     ENDDO
   ENDDO

END SUBROUTINE COSP_GRIDBOX_MEAN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------- SUBROUTINE COSP_CHANGE_VERTICAL_GRID ----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_CHANGE_VERTICAL_GRID(Npoints,Ncolumns,Nlevels,zfull,zhalf,y,   &
                                  Nglevels,newgrid_bot,newgrid_top,r,log_units)
   IMPLICIT NONE
   ! Input arguments
   INTEGER,INTENT(IN) :: Npoints  !# of grid points
   INTEGER,INTENT(IN) :: Nlevels  !# of levels
   INTEGER,INTENT(IN) :: Ncolumns !# of columns
   REAL,DIMENSION(Npoints,Nlevels),INTENT(IN) :: zfull
                  ! Height at model levels [m] (Bottom of model layer)
   REAL,DIMENSION(Npoints,Nlevels),INTENT(IN) :: zhalf
                  ! Height at half model levels [m] (Bottom of model layer)
   REAL,DIMENSION(Npoints,Ncolumns,Nlevels),INTENT(IN) :: y
                  ! Variable to be changed to a different grid
   INTEGER,INTENT(IN) :: Nglevels  !# levels in the new grid
   REAL,DIMENSION(Nglevels),INTENT(IN) :: newgrid_bot
                  ! Lower boundary of new levels  [m]
   REAL,DIMENSION(Nglevels),INTENT(IN) :: newgrid_top
                  ! Upper boundary of new levels  [m]
   LOGICAL,OPTIONAL,INTENT(IN) :: log_units
                  ! log units, need to convert to linear units
   ! Output
   REAL,DIMENSION(Npoints,Ncolumns,Nglevels),INTENT(OUT) :: r
                  ! Variable on new grid

   ! Local variables
   INTEGER :: i,j,k
   LOGICAL :: lunits
   INTEGER :: l
   REAL :: w ! Weight
   REAL :: dbb, dtb, dbt, dtt
             ! Distances between edges of both grids
   INTEGER :: Nw  ! Number of weights
   REAL :: wt  ! Sum of weights
   REAL,DIMENSION(Nlevels) :: oldgrid_bot,oldgrid_top
             ! Lower and upper boundaries of model grid
   REAL :: yp ! Local copy of y at a particular point.
              ! This allows for change of units.

   lunits=.FALSE.
   IF (PRESENT(log_units)) lunits=log_units

   r = 0.0

   DO i=1,Npoints
     ! Calculate tops and bottoms of new and old grids
     oldgrid_bot = zhalf(i,:)
     oldgrid_top(1:Nlevels-1) = oldgrid_bot(2:Nlevels)
     oldgrid_top(Nlevels) = zfull(i,Nlevels) +  zfull(i,Nlevels) -             &
                            zhalf(i,Nlevels) ! Top level symmetric
     l = 0 ! Index of level in the old grid
     ! Loop over levels in the new grid
     DO k = 1,Nglevels
       Nw = 0 ! Number of weigths
       wt = 0.0 ! Sum of weights
       ! Loop over levels in the old grid and accumulate
       ! total for weighted average
       DO
         l = l + 1
         w = 0.0 ! Initialise weight to 0
         ! Distances between edges of both grids
         dbb = oldgrid_bot(l) - newgrid_bot(k)
         dtb = oldgrid_top(l) - newgrid_bot(k)
         dbt = oldgrid_bot(l) - newgrid_top(k)
         dtt = oldgrid_top(l) - newgrid_top(k)
         IF (dbt >= 0.0) EXIT ! Do next level in the new grid
         IF (dtb > 0.0) THEN
           IF (dbb <= 0.0) THEN
             IF (dtt <= 0) THEN
               w = dtb
             ELSE
               w = newgrid_top(k) - newgrid_bot(k)
             ENDIF
           ELSE
             IF (dtt <= 0) THEN
               w = oldgrid_top(l) - oldgrid_bot(l)
             ELSE
               w = -dbt
             ENDIF
           ENDIF
           ! If layers overlap (w/=0), then accumulate
           IF (w /= 0.0) THEN
             Nw = Nw + 1
             wt = wt + w
             DO j=1,Ncolumns
               IF (lunits) THEN
                 IF (y(i,j,l) /= R_UNDEF) THEN
                   yp = 10.0**(y(i,j,l)/10.0)
                 ELSE
                   yp = 0.0
                 ENDIF
               ELSE
                 yp = y(i,j,l)
               ENDIF
               r(i,j,k) = r(i,j,k) + w*yp
             ENDDO
           ENDIF
         ENDIF
       ENDDO
       l = l - 2
       IF (l < 1) l = 0
       ! Calculate average in new grid
       IF (Nw > 0) THEN
         DO j=1,Ncolumns
           r(i,j,k) = r(i,j,k)/wt
         ENDDO
       ENDIF
     ENDDO
   ENDDO

   ! Set points under surface to R_UNDEF, and change to dBZ if necessary
   DO k=1,Nglevels
     DO j=1,Ncolumns
       DO i=1,Npoints
         IF (newgrid_top(k) > zhalf(i,1)) THEN ! Level above model bottom level
           IF (lunits) THEN
             IF (r(i,j,k) <= 0.0) THEN
               r(i,j,k) = R_UNDEF
             ELSE
               r(i,j,k) = 10.0*LOG10(r(i,j,k))
             ENDIF
           ENDIF
         ELSE ! Level below surface
           r(i,j,k) = R_GROUND
         ENDIF
       ENDDO
     ENDDO
   ENDDO

END SUBROUTINE COSP_CHANGE_VERTICAL_GRID

END MODULE MOD_COSP_STATS
