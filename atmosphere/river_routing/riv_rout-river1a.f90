! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE riv_rout_mod_1A

USE routedbl_mod, ONLY: routedbl
IMPLICIT NONE

CONTAINS


SUBROUTINE RIV_ROUT_1A(global_row_length, global_rows,                 &
     surf_runoffin, sub_runoffin, row_length, rows, invert_atmos, xua, &
     yua, river_row_length, river_rows, frac, rmdi, regrid, ru,        &
     ratmed, dt, cyclic_trip, gridbox_areas, global_trip, invert_trip, & 
     twatstor,                                                         & 
     riverdir, riverseq, riverout_trip, runoff_trip_kgs, riverout_inv, &
     inlandout_trip, inlandout_atm, g_river_row_length, g_river_rows,  &
     global_land_frac,first, cmessage)
  !
  ! Purpose:
  !
  ! Perform the routing of surface and sub-surface runoff in parallel
  !
  ! Method:
  ! This routine regrids the total runoff to the river routing grid
  ! and passes it to the TRIP routines to be routed, then maps the
  ! outflow to seapoints on the ATMOS grid and passes back the
  ! updated water storage. The subroutine DO_MAP_MAX contained in
  ! this deck uses the mapping information obatained by PREAV1.dk
  ! to map the river outflow (Kg/s) produced on the river routing
  ! grid onto sea points in the ATMOS grid. This is then converted
  ! to the usual Kg/m2/s using the ATMOS gridbox areas.
  !

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE UM_ParVars, ONLY: g_datastart_f, fld_type_r
  USE UM_ParCore, ONLY: mype

  USE riv_concerns, ONLY: runoff_send_concern, runoff_recv_concern,     &
       riverout_send_concern, riverout_recv_concern,                    &
       riverout_contribution, outflow_send_concern,                     &
       outflow_recv_concern, outflow_contribution
  USE regrid_alloc_calc, ONLY: calc_alloc_concerns,                     &
       calc_alloc_max_concerns
  USE regrid_utils, ONLY: atmos_grid, riv_grid
  USE regrid_swap, ONLY: swap_atmos_grid, swap_riv_grid,                &
       swap_regrid_data, swap_regrid_data_max
  USE PrintStatus_mod

  USE do_areaver_mod, ONLY: do_areaver
  USE do_map_max_mod, ONLY: do_map_max
  USE pre_areaver_mod, ONLY: pre_areaver
  IMPLICIT NONE
 
  INTEGER, INTENT(IN) ::                                                &
       
         rows                                                           &
       ! atmos rows
       , row_length                                                     &
       ! atmos row length
       , global_rows                                                    &
       ! total rows in global grid across all procs
       , global_row_length                                              &
       ! total row length in global grid across all procs
       , river_row_length                                               &
       ! river row length
       , river_rows                                                     &
       ! river rows
       , g_river_row_length                                             &
       ! global river row length          
       , g_river_rows                                                   
       ! global river rows

  REAL, INTENT(IN) ::                                                   &
       
         surf_runoffin(row_length, rows)                                &
       ! atmos grid subdomain surface runoffin 
       , sub_runoffin(row_length, rows)                                 &
       ! atmos grid subdomain surface runoffin
       , frac(row_length, rows)                                         &
       ! land fractions 
       , riverdir(river_row_length, river_rows)                         &
       ! values of river direction
       , riverseq(river_row_length, river_rows)                         &
       ! values of river sequence 
       , gridbox_areas(row_length, rows)                                & 
       ! atmos gridbox areas
       , global_land_frac(global_row_length, global_rows)              
       ! land seak fraction across all processes 

  REAL, INTENT(INOUT) ::                                                &
       
         riverout_inv(row_length, rows)                                 &
       ! river flow out from each gridbox
       , riverout_trip(river_row_length, river_rows)                    &
       ! gridbox outflow on river grid (Kg/s)
       , runoff_trip_kgs(river_row_length, river_rows)                  &
       ! gridbox runoff on river grid 
       , inlandout_trip(river_row_length, river_rows)                   &
       ! trip outflow from inland basins on trip grid (kg/s)
       , twatstor(river_row_length, river_rows)                         &
       ! initial water storage file 
       , inlandout_atm(row_length, rows)                    
       ! trip outflow from inland basins on atmos grid kg/m2/s

  REAL, INTENT(IN) ::                                                   &
       
         xua(0:global_row_length)                                       &
       ! atmos grid global domain lat coords
       , yua(0:global_rows)                     
       ! atmos grid global domain long coords 

  LOGICAL, INTENT(IN) ::                                                &
       
         invert_atmos                                                   &
       ! invert rows 
       , regrid                                                         &
       ! whether to regrid atmos to river grid 
       , cyclic_trip                                                    &
       ! indicates if trips grid is cyclic 
       , global_trip                                                    &
       ! indicates if atmos grid is cyclic 
       , invert_trip

  LOGICAL, INTENT(IN) :: first 
  ! true if first entry into this routine

  REAL, INTENT(IN) ::                                                   &
       
         rmdi                                                           &
       ! missing value indicator
       , ru                                                             &
       ! effective river flow velocity (m/s)
       , ratmed                                                         &
       ! river meander ratio
       , dt 
       ! river routing time step (s)

  CHARACTER(LEN=80), INTENT(OUT) :: cmessage     
  ! error message 

  !! local variables 

  REAL runoff_trip(river_row_length, river_rows)                        &
                                ! regridded runoff 
       , trip_inlandout(river_row_length, river_rows)                   &
       , riverout_trip_mouth(river_row_length, river_rows)              &
                                ! gridbox outflow and river mouth
       , trip_outflow(river_row_length, river_rows)                     

  REAL, DIMENSION(:), ALLOCATABLE ::                                    &
         weight
  ! for weighting of source grid point 
  INTEGER, DIMENSION(:), ALLOCATABLE ::                                 &
        index_arav
  ! indexes global src grid points

  REAL runoffin(row_length,rows)                                        &
                                ! atmos runoffin with halos for swap bounds 
       , riverout(row_length,rows)                                      &
       , inlandout(row_length, rows)

  INTEGER                                                               & 
       
        count_tr(RIVER_ROW_LENGTH, RIVER_ROWS)                          &
                                ! no. of source boxes in target 
       ,base_tr(river_row_length, river_rows)                           &
       ,trivdir(river_row_length, river_rows)                           &
                                ! river direction file in integers
       ,trivseq(river_row_length, river_rows)                           &
                                ! river sequence file in integers
       , trivdir_inv(river_row_length, river_rows)                      
  LOGICAL amask_all(global_row_length, global_rows)                     &
                                ! set atmos mask to all sea for interp
       
       , trmask_all(river_row_length, river_rows)                       &
                                !  set trip mask to all sea for interp
       , trmask(river_row_length, river_rows)                           &
       , amask_sea(row_length,rows)                                     &  
                                ! ATMOS L/sea mask, F for any sea
       , amask_sea_global(global_row_length, global_rows)               
  ! ATMOS L/sea  global mask, F for any sea                        

  INTEGER lenl
  ! the maximum number of source boxes in weights for all targets

  CHARACTER(LEN=*) routinename
  PARAMETER ( routinename='RIV_ROUT_1A')

  INTEGER                                                               &
         local_rstart_rl, local_rstart_rows                             &
       , local_rstop_rl, local_rl, local_pstop_rl                       &
       , local_rstop_rows, stop_phi, start_lambda, local_rows           &
       , local_pstart_row, local_pstart_rl, local_pstop_row 

  REAL                                                                  &
       
         yut(river_rows+1)                                              &
                                ! TRIP latitude coordinates (limited)
       , xut(0:river_row_length)                                        
  ! TRIP longitude coordinates (limited)
  REAL                                                                  &      
        trip_areas(river_row_length, river_rows)

  INTEGER i, j, k, l  
  ! loop counter variables

  INTEGER jmax                      
  ! max. no. of j values for TRIP

  INTEGER ndev               ! No. of TRIP timesteps/day
  INTEGER error 

  ! debugging variables
  INTEGER ierror, my_comm, info, icode
  REAL total_weight_sum, weight_sum

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('RIV_ROUT_1A',zhook_in,zhook_handle)

  inlandout = 0.0
  inlandout_trip = 0.0
  inlandout_atm = 0.0
  trip_inlandout = 0.
  error = 0
  info = 0

  !! calculate runoff 

  ! Sum the surface and subsurface runoffs
  ! Sum the two types of runoff

  DO j = 1, rows
    DO i = 1, row_length
      IF(surf_runoffin(i,j) /= rmdi                                     &
                  .AND. sub_runoffin(i,j) /= rmdi)THEN
        runoffin(i,j) = (surf_runoffin(i,j)                             &
                                       + sub_runoffin(i,j))*frac(i,j)
      ELSE
        runoffin(i,j) = rmdi
      ENDIF
    ENDDO
  ENDDO

  DO j=1, rows
    DO i=1, row_length
      ! .FALSE. for gridboxes with any sea
      amask_sea(I,J) = (frac(i,j) == 1.0)
    END DO
  END DO


  IF(regrid)THEN

    !! calculate subdomain TRIP ranges with 1 degree res (lat and long) 
    !! TODO: make generic to resolution
    !! latitude is north to south (90 -> -90) 
    !! and longitude is west to east (0 -> 360)

    local_rstart_rl = g_datastart_f(1, fld_type_r, mype)
    local_rstop_rl = local_rstart_rl + RIVER_ROW_LENGTH - 1  

    local_rstart_rows = g_datastart_f(2, fld_type_r, mype)
    local_rstop_rows = local_rstart_rows + RIVER_ROWS - 1

    !! add offset to stop and start point so lat coords reads 90 to -90
    local_rstart_rows = local_rstart_rows - 90
    local_rstop_rows = local_rstop_rows - 90


    !! south to north(1 degree resolution)
    !! account for top of row 
    local_rstart_rows = local_rstart_rows - 1

    j = 1
    DO i = local_rstart_rows, local_rstop_rows
      yut(j)= i
      j = j + 1
    END DO

    !! east to west (1 degree resolution)
    !! account for proc account for last degree 
    !! account for left of row
    local_rstart_rl = local_rstart_rl - 1

    j = 0
    DO i = local_rstart_rl, local_rstop_rl
      xut(j)= i
      j = j + 1 
    END DO

    ! initialise subdomain field        
    DO j=1, global_rows
      DO i=1, global_row_length
        amask_sea_global(i,j) = (global_land_frac(i,j) == 1.0)
      END DO
    END DO

    !! calculate local weights 
    lenl = (global_row_length+river_row_length)*(global_rows+river_rows)

    amask_all = .FALSE.
    trmask_all = .FALSE.

    ALLOCATE(weight(lenl))
    ALLOCATE(index_arav(lenl)) 

    !! calculate area averages 

    CALL PRE_AREAVER(global_row_length,xua,global_rows,yua,             &   
         & .TRUE.,global_row_length,.FALSE.,amask_all,river_row_length, &
         & xut,river_rows,yut,.FALSE.,.TRUE., lenl,count_tr,base_tr,    &
         & index_arav,weight,icode,cmessage)

    ! calculate concerns 
    IF(first) THEN

      ! for forward regridding 
      CALL CALC_ALLOC_CONCERNS(runoff_send_concern, runoff_recv_concern,&
           lenl, index_arav, global_row_length, global_rows, atmos_grid,&
           error, cmessage)

    END IF

    CALL SWAP_REGRID_DATA(runoff_send_concern, runoff_recv_concern,    &
         runoffin, swap_atmos_grid, row_length, rows, icode)

    ! intialise to no data
    DO j=1,river_rows
      DO i=1,river_row_length
        runoff_trip(i,j)=rmdi
      ENDDO
    ENDDO

    i = SIZE(runoff_recv_concern, 1)

    CALL DO_AREAVER(row_length, rows, runoffin,  river_row_length,      &
         river_rows, count_tr, base_tr, river_row_length, .FALSE.,      &  
         trmask_all, index_arav, weight, runoff_trip, global_row_length,& 
         global_rows, atmos_grid, runoff_recv_concern, i)
    
  ELSE

    DO j=1, river_rows
      DO i=1, river_row_length
        runoff_trip(i,j) = runoffin(i,j)
      END DO
    END DO

  ENDIF ! regridding

  DO j=1,river_rows
    DO i=1, river_row_length
      ! Change the runoff from Kg/m2/sec (mm/sec) to mm/day for TRIP routines
      ! Allow Runoff to be at TRIP seapoints for conserve.

      IF(runoff_trip(i,j) >  rmdi)THEN
        runoff_trip(i,j)=runoff_trip(i,j)*86400
      ELSE
        runoff_trip(i,j) = 0.0
      ENDIF

    ENDDO
  ENDDO


  !******************************************************************
  ! Change river direction and sequence arrays to integer as expected by
  ! river routing scheme

  DO j=1, river_rows
    DO i=1, river_row_length
      trivdir(i,j) = INT(riverdir(i,j))
      trivseq(i,j) = INT(riverseq(i,j))
    ENDDO
  ENDDO

  ! Set the number of river routing timestes per day to 1
  ndev = 1

  IF (printstatus >= prstatus_diag) THEN
  WRITE (6, '(A,I6)') 'number of routing timesteps per day = ', ndev
  WRITE (6, '(A,F8.2,A)') 'routing timestep =  ', DT, ' secs.'
  WRITE (6, '(A,F6.4)') 'Meandering Ratio  = ', ratmed
  WRITE (6, '(A,F6.4,A)') 'Effective Velocity ru =', ru, ' (m/s)'
  END IF

  jmax=river_rows
 
  ! need to pass offset values for routine, remember to lat rotate grid
  CALL ROUTEDBL(runoff_trip, ru, ratmed, ndev, dt                       &
       , river_row_length, river_rows, trivdir, trivseq, twatstor       &
       , jmax, rmdi, riverout_trip, runoff_trip_kgs, trip_areas         &
       , local_rstart_rl, (local_rstart_rows+90))

  ! Invert the TRIP OUTFLOW for 'mapping' and remove TRIP_RMDI values

  !AJW fix bug where water is only regridded for ocean pour points thus excluding inland
  !basins. TRIVDIR(I,J) now >= 9 rather than ==9.

  ! Use only outflow at river mouths and sea for the mapping onto ATMOS
  ! grid
  DO j=1, river_rows
    DO i=1, river_row_length
      IF(trivdir(i,j) == rmdi .OR. trivdir(i,j) >= 9 .OR. trivdir(i,j)  &
           == 0)THEN

        riverout_trip_mouth(i,j) = riverout_trip(i,j)
      ELSE
        riverout_trip_mouth(i,j) = 0.0
      ENDIF
    END DO
  END DO

  ! For DO_MAP_MAX, do not allow RMDITRIP values
  DO j=1, river_rows
    DO i=1, river_row_length

      IF(riverout_trip_mouth(i,j) /= rmdi)THEN
        trip_outflow(i,j)=riverout_trip_mouth(i,j)            
      ELSE
        trip_outflow(i,j)=0.0
      ENDIF
    END DO
  END DO


  IF (regrid) THEN

    ! Prepare for mapping and count river mouths on TRIP grid as sea.
    ! Instead of regridding, use count_tr, base_tr and weight to 'map' the
    ! TRIP sea points onto the ATMOS grid where the major part of the TRIP
    ! box lies using do_map_max

    DO j=1, river_rows
      DO i=1, river_row_length
        trmask(i,j)= (riverout_trip(i,j) == rmdi)
      ENDDO
    ENDDO

    riverout = 0
    count_tr = 0
    base_tr = 0

    lenl = (global_row_length+river_row_length)*(global_rows+river_rows)

    CALL PRE_AREAVER(global_row_length,xua,global_rows,yua, .TRUE.,     &   
         global_row_length,.FALSE., amask_sea_global,                   &
         river_row_length,xut,river_rows,yut, .FALSE.,.TRUE.,lenl,      &   
         count_tr,base_tr,index_arav, weight,icode,cmessage)

    IF (first) THEN

      ! for regridding back from trips
      CALL CALC_ALLOC_MAX_CONCERNS(riverout_send_concern,               &
           riverout_recv_concern, riverout_contribution,                &
           base_tr, count_tr, river_row_length, river_rows,             &
           weight, riv_grid, index_arav, lenl, global_row_length,       &
           global_rows, row_length, rows, atmos_grid, trmask, .false.,  &
           g_river_row_length, icode, cmessage)

    END IF

    ! swap src field data 

    CALL SWAP_REGRID_DATA_MAX(trip_outflow, river_row_length,           &
         river_rows, riv_grid, riverout_send_concern,                   &
         riverout_recv_concern, icode, cmessage)


    i = SIZE(riverout_contribution, 1)
    j = SIZE(riverout_recv_concern, 1)

    CALL DO_MAP_MAX(trip_outflow, river_rows, river_row_length,        &
         riverout, row_length, rows, riv_grid, trmask, .false.,        &
         riverout_recv_concern, j, riverout_contribution, i)

    !********************************************************************
    ! Inland basin now 10
    !  calculate  inland basin outflow using grid points where
    !  count_tr=0, trivdir=10, rmdi or 0, and count_tr=0

    DO J=1, river_rows
      DO I=1, river_row_length

        IF(invert_atmos) THEN

          IF((trivdir(i,j) >= 9 .OR. trivdir(i,j) == rmdi.OR.           &
               & trivdir(i,j) == 0) .AND. count_tr(i,j) == 0) THEN

            trip_inlandout(i,j)=trip_outflow(i,j)
            inlandout_trip(i,j)=trip_outflow(i,j)

          ELSE

            trip_inlandout(i,j)=0.0
            inlandout_trip(i,j)=0.0

          ENDIF

        ENDIF

      ENDDO
    ENDDO

    ! call pre_areaver to regrid inland basin outflow
    ! swap src field data 

    weight = 0
    index_arav = 0
    inlandout = 0

    lenl = (global_row_length+river_row_length)*(global_rows+river_rows)

    CALL PRE_AREAVER(global_row_length,xua,global_rows,yua,.TRUE.,      &   
          global_row_length,.FALSE.,amask_all,river_row_length,         &
          xut,river_rows,yut,.FALSE.,.TRUE., lenl,                      &
          count_tr,base_tr,index_arav,weight,icode,cmessage)

    IF (first) THEN

      ! for backwards regridding 
      CALL CALC_ALLOC_MAX_CONCERNS(outflow_send_concern,                &
           outflow_recv_concern, outflow_contribution,                  &
           base_tr, count_tr,                                           &  
           river_row_length, river_rows, weight, riv_grid,              & 
           index_arav, lenl, global_row_length,                         &
           global_rows, row_length, rows, atmos_grid, trmask, .false.,  &
           g_river_row_length, icode, cmessage)

    END IF

    CALL SWAP_REGRID_DATA_MAX(inlandout_trip, river_row_length,         &
         river_rows, riv_grid, outflow_send_concern,                    &
         outflow_recv_concern, icode, cmessage)


    i = SIZE(outflow_contribution, 1)
    j = SIZE(outflow_recv_concern, 1)

    CALL DO_MAP_MAX(inlandout_trip, river_rows, river_row_length,     &
         inlandout, row_length, rows, riv_grid, trmask, .false.,        &
         outflow_recv_concern, j, outflow_contribution, i)

    DEALLOCATE(weight)
    DEALLOCATE(index_arav)

  END IF ! end of regridding 


  !******************************************************************


  ! Change river outflow to mm/s (as expected by the ocean later)

  DO J=1,rows
    DO I=1,row_length
      IF(riverout(i,j) /= rmdi)THEN
        riverout_inv(i,j) =                                       &
             &          riverout(i,j)/gridbox_areas(i,j)
      ELSE
        riverout_inv(i,j) = rmdi
      END IF
    END DO
  END DO


  ! convert inland basin outflow from kg/s to kg/m2/s
  ! as expected by hydrol7a

  DO j=1, rows
    DO i=1, row_length

      IF(inlandout(i,j) /= rmdi)THEN
        inlandout_atm(i,j)= inlandout(i,j)/gridbox_areas(i,j)
      ELSE

        !      set to zero, not rmdi since hydrol7a cannot
        !      use rmdi values

        inlandout_atm(i,j)=0.0
      ENDIF

    END DO
  END DO

  ! Set land points values of inlandout_trip to rmdi post regridding

  DO j=1, river_rows
    DO i=1, river_row_length

      IF(riverout_trip(i,j) == rmdi) inlandout_trip(i,j)=rmdi          

    END DO
  END DO

  ! Set RMDI values in river outflow and runoff on river grid to 0.0 as
  ! points with RMDI are not constant as runoff greater than 0.0 is
  ! 'added in' at some TRIP seapoints due to regridding from
  ! non-coincident grids.
  ! Otherwise postprocessing of fields can give negative values.

  DO j=1,river_rows
    DO i=1, river_row_length
      IF(riverout_trip(i,j) == rmdi) riverout_trip(i,j)=0.0
    END DO
  END DO

  IF (lhook) CALL dr_hook('RIV_ROUT_1A',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE RIV_ROUT_1A

!*********************************************************************

END MODULE riv_rout_mod_1A
