! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE routedbl_mod

USE setcoef_mod, ONLY: setcoef
USE setrval_mod, ONLY: setrval
USE wrtwblog_mod, ONLY: wrtwblog
IMPLICIT NONE

CONTAINS


SUBROUTINE ROUTEDBL(runoff, ru, ratmed, ndev, dt, nx, ny               &
     , igrcn, iseq, sto, jmax, rmiss, dout, drunin, area, offset_nx    &
     , offset_ny)
!
!  Purpose: To route runoff
!  Method:
!  The runoff is routed along predetermined routes using an initial
!  water storage which is updated each timestep. The code and
!  river direction and sequence files were produced by Taikan Oki
!  (March 26, 1998 - see ) as
!  standalone code and was adapted (as little as possible)
!  by Cyndy Bunton to run in the UM GCM. The spinup was removed
!  and initial water storage, river sequence and direction files
!  input via the Atmosphere dump. Gridbox outflow, current water
!  storage and total gridbox runoff are output via stash.

! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE MPL, ONLY : MPL_INTEGER, MPL_MAX, MPL_REAL, MPL_SUM
  USE global_2d_sums_mod, ONLY: global_2d_sums
  USE UM_ParVars
  
  USE cp2_mod, ONLY: cp2
  USE initial0_mod, ONLY: initial0
  USE mmd2kgs_mod, ONLY: mmd2kgs
  USE outflow1_mod, ONLY: outflow1
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::                                                &
        ndev                                                            &
       ! number of river routing timesteps/day
       ,jmax                                                            &
       ! number of rows 
       ,nx                                                              &
       ! number of points in a row
       ,ny                                                              &
       ! number of rows
       ,igrcn(nx,ny)                                                    &
       ! gridbox outflow direction
       ,iseq(NX,NY)                                                     &
       ! river flow sequence
       ,offset_nx                                                       &
       ! offset from global lat 0                      
       ,offset_ny                                                      
       ! offset from global long 0

  REAL, INTENT(IN) ::                                                   &
        runoff(nx,ny)                                                   &
       ! total runoff (Kg/m2/day)
       ,dt                                                              &
       ! river routing timestep (s)
       ,ru                                                              &
       ! effective river flow velocity (m/s)
       ,ratmed                                                          &
       ! River meander ratio
       ,rmiss          
       ! missing data indicator
  
  ! IN/OUT Arguments
  REAL, INTENT(INOUT) ::                                                &
        sto(nx,ny)      ! River water storage (Kg)
  ! OUT arguments
  REAL, INTENT(OUT) ::                                                  &
        DOUT(nx,ny)                                                     &
       ! gridbox outflow (Kg/s)
       ,DRUNIN(nx,ny)   
  ! gridbox inflow (Kg/s)
  
  ! local variablea
  
  
  ! halo size of 2, first halo layer for points which advected into
  ! subdomain, 2nd layer received data from points that advect out
  ! from subdomain (2nd layer not used, just to prevent segfaults) 
  INTEGER, PARAMETER :: halo_x=2, halo_y=2

  ! extended arrays for passed in values  

  REAL                                                                  &
          sto_ext(1-halo_x:nx+halo_x,1-halo_y:ny+halo_y)                &    
       ! River water storage (Kg)
       , drunin_ext(1-halo_x:nx+halo_x,1-halo_y:ny+halo_y)              &
       ! gridbox inflow (Kg/s)
       , dout_ext(1-halo_x:nx+halo_x,1-halo_y:ny+halo_y)                &
       ! gridbox outflow (Kg/s)
       , rc_ext(1-halo_x:nx+halo_x,1-halo_y:ny+halo_y)                  &
       ! transfer coefficient   [1/s]
       , din_ext(1-halo_x:nx+halo_x, 1-halo_y:ny+halo_y)                &
       ! total input to river storage [kg/s]
       , sto2_ext(1-halo_x:nx+halo_x, 1-halo_y:ny+halo_y) 
  ! river channel storage [kg] at end of step
  
  INTEGER                                                               &
       &  igrcn_ext(1-halo_x:nx+halo_x,1-halo_y:ny+halo_y)              &
       ! river direction
       , iseq_ext(1-halo_x:nx+halo_x,1-halo_y:ny+halo_y)                &  
       ! river sequence
       , inextx_ext(1-halo_x:nx+halo_x, 1-halo_y:ny+halo_y)             &
       ! next point downstream (x direction)
       , inexty_ext(1-halo_x:nx+halo_x, 1-halo_y:ny+halo_y)              
  ! next point downstream (y direction)
  
  REAL                                                                  &
         rlen(nx,ny)                                                    &
       ! distance between grids [m]
       , rvel(nx, ny)                                                   &
       ! flow velocity          [m/s]
       , rc(nx, ny)                                                     &
       ! transfer coefficient   [1/s]
       , area(nx, ny)                                                   &
       ! grid area              [m^2]
       , din(nx, ny)                                                    &
       ! total input to river storage [kg/s]
       , sto2(nx, ny)                                                   &                 
       ! river channel storage [kg] at end of step
       , dtotal(nx,ny)                                                  &
       
       , dinput(nx, ny)
  
  REAL stoall, dinall, sto2all, doutall                                 &
       , drunall, drivall, tot_sto2all, tot_drivall                     &
       , tot_doutall, tot_drunall ! totals
  
  INTEGER                                                               &
       & inextx(nx, ny)                                                 &
       ! next point downstream (x direction)
       , inexty(nx, ny) 
  ! next point downstream (x direction)
  
  INTEGER iy, im, idec, ih ! year,month,day,hour
  INTEGER nseqmax        ! maximum points in any river sequence
  INTEGER i, j           ! loop counters
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  INTEGER info, my_comm, ierror

  !** required for inland basin correction
  
  REAL ::                                                               &
       & tot_flux2ocean                                                 &
       ! total global flux to ocean 
       , tot_basin_flux                                                  
  ! total basin flux before advection 

  INTEGER, PARAMETER :: N_LOG_VALS = 6

  REAL ::                                                               &
           &  sum_r(nx, ny, 2)                                          &
           ! temp var for calculating reproducible basin flux global sums
           , sum_all(2)                                                 &
           ! temp var for basin_flux sum (1) and flux2ocean (2) 
           , log_val(N_LOG_VALS), log_val_sum(N_LOG_VALS)
           ! used to calculate sums across decomposition of log variables
  !**
  
  INTEGER, SAVE :: nseqmax_global = 0 
  ! maximum river sequence 
  LOGICAL, SAVE :: first = .TRUE. 
  ! first time entered routine
  
      
  IF (lhook) CALL dr_hook('ROUTEDBL',zhook_in,zhook_handle)
  
  info = 1
  CALL gc_get_communicator(my_comm,info)
  
  iy = 0
  im = 0
  idec = 0

  
  ! initialise halo extended variables to zero, halo initialised to zero also 
  
  inextx_ext = 0
  inexty_ext = 0
  igrcn_ext = 0
  iseq_ext = 0
  
  din_ext = 0.
  sto_ext = 0.
  drunin_ext = 0.
  dout_ext = 0.
  sto2_ext = 0.
  din = 0.
  sum_r = 0.
  sum_all = 0.
  drunin = 0.
  rc_ext = 0.
  
  CALL initial0(nx, ny, jmax, rmiss, igrcn, iseq, nseqmax,          &
       &   inextx, inexty, rlen, area, offset_nx, offset_ny)
  
  ! fix. the velocity here
  
  CALL setrval(nx, ny, rvel, ru, jmax)
  
! set meander coefficient

  CALL setcoef(nx, ny, rlen, rvel, ratmed, rc, jmax)

! convert runoff from mm/day to kg/s

  CALL mmd2kgs(runoff, igrcn, drunin, area, nx, ny, rmiss, jmax)
  
  CALL setrval(nx, ny, dtotal, 0.0, ny)
  
  CALL setrval(nx, ny, dinput, 0.0, ny)
  
  
  !**********************************************************
  ! Inland basin fix
  ! drunin is in kg/s
  
  sum_r = 0.
  DO j=1, ny
    DO i=1, nx
      IF(igrcn(i,j)==10) THEN 
        sum_r(i,j,1) = drunin(i,j) 
      ELSE
        sum_r(i,j,1) = 0.
      END IF
    END DO
  END DO
  
  where(igrcn .eq. 10) drunin=0.
  
  !**********************************************************
   
  ! initialise subdomain of extended fields with original fields
  
  DO j=1, ny
    DO i=1, nx
      sto_ext(i,j) = sto(i,j)
      drunin_ext(i,j) = drunin(i,j) 
      dout_ext(i,j) = dout(i,j) 
      igrcn_ext(i,j) = igrcn(i,j) 
      iseq_ext(i,j) = iseq(i,j)
      
      rc_ext(i,j) =  rc(i,j)
      inextx_ext(i,j) = inextx(i,j)
      inexty_ext(i,j) = inexty(i,j)
      din_ext(i,j) = din(i,j)
      
    END DO
  END DO
  
  ! transform inext halo to account for grid points oriented to advect out
  ! orginal subdomain so they advect inwards 
  ! (i.e. halo's to advect out should reverse direction so they advect
  ! into subdomain)
  
  ! need to retrieve adjacent subdomains row and row_length
  ! for neighbour processes, since river domain may not be 
  ! partitioned equally EW or NS. 

  ! correcting for halo orientation
  DO J = 1, NY
    DO I = 1, NX
      
      IF ((igrcn(i,j) >= 1) .AND. (igrcn(i,j) <= 8)) THEN
      
        ! a zero boundary means prepare bottom halo, of adjacent 
        ! subdomain south, to advect upwards
        IF (inexty(i,j)  ==  0) THEN 
          inexty_ext(i,j) =  g_blsize(2, fld_type_r, neighbour(psouth))
          ! see **
          IF ((inextx(i,j)  /=  0) .AND. (inextx(i,j) /= (NX+1) ))      &
               inextx_ext(i,j) =  inextx(i,J)  
        END IF

        ! a ny+1 boundary means prepare top halo, of adjacent 
        ! subdomain below, to advect downwards
        IF (inexty(i,j)  ==  NY+1) THEN
          inexty_ext(i,j) = 1
          ! see **
          IF ((inextx(i,j)  /=  0) .AND. (inextx(i,j) /= (NX+1) ))      &
               inextx_ext(i,j) =  inextx(i,j)
        END IF
        
        ! a 0 boundary means prepare right halo, of adjacent 
        ! subdomain ea, to advect leftwards
        IF (inextx(i,j)  ==  0) THEN
          inextx_ext(i,j) = g_blsize(1, fld_type_r, neighbour(pwest))
          ! see **
          IF ((inexty(i,j)  /=  0) .AND. (inexty(i,j) /= (NY+1) ))      &
               inexty_ext(i,j) =  inexty(i,j) 
          
        END IF

        ! a nx+1 boundary means prepare left halo, of adjacent 
        ! subdomain right, to advect rightwards
        IF (inextx(i,j)  ==  NX+1) THEN
          inextx_ext(i,j) = 1
          ! see **
          IF ((inexty(i,j)  /=  0) .AND. (inexty(i,j) /= (NY+1) ))      &
               inexty_ext(i,j) =  inexty(i,j)
        END IF
        
      ENDIF
      
    END DO
  END DO

  ! ** don't forget to initialise partner inext if it's not 
  ! advecting across a boundary
  
! DEPEND ON: swap_bounds 
  CALL SWAP_BOUNDS(sto_ext, nx, ny, 1, halo_x, halo_y, fld_type_r,      &
       .FALSE.) 
! DEPEND ON: swap_bounds 
  CALL SWAP_BOUNDS(drunin_ext, nx, ny, 1, halo_x, halo_y,               &
       & fld_type_r, .FALSE.) 
! DEPEND ON: swap_bounds 
  CALL SWAP_BOUNDS(igrcn_ext, nx, ny, 1, halo_x, halo_y,                &
       & fld_type_r, .FALSE.) 
! DEPEND ON: swap_bounds 
  CALL SWAP_BOUNDS(iseq_ext, nx, ny, 1, halo_x, halo_y,                 &
       & fld_type_r, .FALSE.) 
! DEPEND ON: swap_bounds 
  CALL SWAP_BOUNDS(rc_ext, nx, ny, 1, halo_x, halo_y,                   &
       & fld_type_r, .FALSE.) 
! DEPEND ON: swap_bounds 
  CALL SWAP_BOUNDS(inextx_ext, nx, ny, 1, halo_x, halo_y,               &
       & fld_type_r, .FALSE.) 
! DEPEND ON: swap_bounds 
  CALL SWAP_BOUNDS(inexty_ext, nx, ny, 1, halo_x, halo_y,               &
       & fld_type_r, .FALSE.) 
! DEPEND ON: swap_bounds 
  CALL SWAP_BOUNDS(din_ext, nx, ny, 1, halo_x, halo_y,                  &
       & fld_type_r, .FALSE.) 


! reinitialise inext to values prior to halo orientation changes so
! as those changes we're a correction adjacent subdomain

  DO J = 1, NY
    DO I = 1, NX
      inextx_ext(I,J) = inextx(I,J)
      inexty_ext(I,J) = inexty(I,J) 
    ENDDO
  ENDDO  

! find global max river sequence (river sequence field does not change)

  IF (FIRST) THEN
    CALL MPL_ALLREDUCE(nseqmax, nseqmax_global, 1, MPL_INTEGER,         &
         & MPL_MAX, my_comm, ierror)
    FIRST = .FALSE.
  END IF
  
! advect river grid points

  CALL outflow1(sto_ext, drunin_ext, rc_ext, dt, igrcn_ext              &
       , iseq_ext, inextx_ext, inexty_ext, nseqmax_global, nx, ny, jmax &
       , sto2_ext, din_ext, dout_ext, drunall, drivall, stoall, sto2all &
       , doutall, dinall, halo_x, halo_y, nseqmax)
  

! copy back from halo variables

  DO j=1, ny
    DO i=1, nx
      
      dout(i,j) = dout_ext(i,j)
      sto2(i,j) = sto2_ext(i,j)
      din(i,j) = din_ext(i,j) 
      
    END DO
  END DO
  
  !***********************************************************
  !Inland basin fix
  !Add flux from inland basin to river outflow points
  
  DO j=1, ny
    DO i=1, nx
      IF(igrcn(i,j) == 9) THEN 
        sum_r(i,j,2) = dout(i,j) 
      ELSE
        sum_r(i,j,2) = 0.
      END IF
    END DO
  END DO
  
! calculate global sums 

  CALL global_2d_sums(sum_r, nx, ny, 0, 0, 2, sum_all)

  tot_basin_flux = sum_all(1)
  tot_flux2ocean = sum_all(2)

  where (igrcn .eq. 9) dout=dout+(tot_basin_flux/tot_flux2ocean)*dout
  
  !***********************************************************
  
  
  ! Add in the runoff at TRIP seapoints produced by interpolation from
  ! ATMOS runoff to conserve due to mismatch between grids
  ! Note that this will not appear in RINP, ROUT written by
  
  DO i = 1, nx
    DO j = 1, jmax
      IF(iseq(i,j) == rmiss)THEN
        IF(drunin(i,j) /= 0.0)THEN
          dout(i,j) = drunin(i,j)
          din(i,j) = drunin(i,j)
          ! Add in to the totals for water balance check
          drunall = drunall + drunin(i,j)
          dinall = dinall + din(i,j)*dt
          doutall = doutall + dout(i,j)*dt
        ELSE
          dout(i,j) = rmiss
          din(i,j) = rmiss
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  log_val(1) = stoall 
  log_val(2) = sto2all
  log_val(3) = dinall 
  log_val(4) = doutall
  log_val(5) = drunall 
  log_val(6) = drivall 
 
  ! warning: does not bit reproduce on different proc decompositions
  CALL MPL_ALLREDUCE(log_val, log_val_sum, N_LOG_VALS, MPL_REAL,        & 
       & MPL_SUM, my_comm, ierror)

  stoall = log_val_sum(1)
  sto2all = log_val_sum(2)
  dinall = log_val_sum(3) 
  doutall = log_val_sum(4)
  drunall = log_val_sum(5) 
  drivall = log_val_sum(6)

  call wrtWBlog(6, iy, im, idec, 0, ndev, stoall, sto2all, dinall,      &
       doutall, drunall, drivall, dt)

  CALL cp2(sto2, sto, nx, ny, jmax)
  
  IF (lhook) CALL dr_hook('ROUTEDBL',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE ROUTEDBL

END MODULE routedbl_mod
