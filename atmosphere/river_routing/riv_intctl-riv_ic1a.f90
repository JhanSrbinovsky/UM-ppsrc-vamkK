! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE riv_intctl_mod_1A

USE riv_rout_mod_1A, ONLY: riv_rout_1A

IMPLICIT NONE

CONTAINS

       SUBROUTINE RIV_INTCTL_1A(                                       &
      XPA, XUA, XVA, YPA, YUA, YVA,                                    &
      G_P_FIELD, G_R_FIELD, N_PROC, ME, RMDI,                          &
      GATHER_PE_TRIP,LAND_POINTS,LAND_INDEX,                           &
      INVERT_ATMOS, ROW_LENGTH, ROWS,                                  &
      GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                                  &
      RIVER_ROW_LENGTH, RIVER_ROWS,                                    &
      GLOBAL_RIVER_ROW_LENGTH, GLOBAL_RIVER_ROWS,                      &
      FLANDG, RIV_STEP, RIV_VEL, RIV_MCOEF,                            &
      TRIVDIR, TRIVSEQ, TWATSTOR, A_BOXAREAS,                          &
      DELTA_PHI, FIRST_ENTRY,                                          &
      R_AREA, SLOPE, FLOWOBS1, R_INEXT, R_JNEXT,                       &
      R_LAND,  SUBSTORE, SURFSTORE, FLOWIN, BFLOWIN,                   &
! IN/OUT accumulated runoff
      TOT_SURF_RUNOFF, TOT_SUB_RUNOFF,                                 &
! OUT
      BOX_OUTFLOW, BOX_INFLOW, RIVEROUT_ATMOS                          &
! Add new arguments for inland basin outflow
! OUT INLAND BASINS
      ,INLANDOUT_ATMOS,INLANDOUT_RIV                                   & 
! Required for soil moisture correction for water conservation 
      ,DSM_LEVELS,ACC_LAKE_EVAP,SMVCST,SMVCWT                          & 
      ,SMCL_DSM,STHU_DSM) 

! Purpose:
! New Control routine for River routing for Global Model.
! Parallel River routing
!

! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!-----------------------------------------------------------------

      USE water_constants_mod, ONLY: rho_water
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE global_2d_sums_mod, ONLY: global_2d_sums
      USE ereport_mod, ONLY: ereport 
      USE UM_ParVars, ONLY: lasize, glsize, fld_type_p,                &
           halo_type_no_halo, gc_all_proc_group
      USE soil_param, ONLY: dzsoil
      IMPLICIT NONE


      INTEGER, INTENT(IN) ::                                           &
       row_length                                                      &
                                 ! NO. OF COLUMNS IN ATMOSPHERE
     , rows                                                            &
                                 ! NO. OF ROWS IN ATMOSPHERE
     , global_row_length                                               &
                                 ! number of points on a row
     , global_rows                                                     &
                                 ! NUMBER OF global rows
     , land_points                                                     &
                                 ! number of landpoints
     , river_row_length                                                &
                                 ! no. of columns in river grid
     , river_rows                                                      &
                                 ! no. of rows in river grid
     , global_river_row_length                                         &
                                 ! global river row length
     , global_river_rows                                               &
                                 ! GLOBAL river rows
     , gather_pe_trip                                                  &
                                 ! pe River routing to be run on
     , n_proc                                                          &
                                 ! Total number of processors
     , me                                                              & 
                                 ! My processor number 
     , dsm_levels               ! no. of deep soil moisture levels 

      INTEGER, INTENT(IN) ::                                           &
       G_P_FIELD                                                       &
                                  ! IN size of global ATMOS field
     , G_R_FIELD                                                       &
                                  ! IN Size of global river field
     , land_index (land_points)  ! IN index of land to global points

      REAL, INTENT(IN) ::                                              &
     & TOT_SURF_RUNOFF(land_points)                                    &
                                   !IN Surf RUNOFF on land pts(KG/M2/S)
     ,TOT_SUB_RUNOFF(land_points)                                      &
                                   ! IN Subsurf.RUNOFF (KG/M2/S)

! Water conservation 
! Remove global lake evaporation from soil moisture store 
     , smvcwt (land_points)                                            & 
!                            ! IN Volumetric wilting point 
                              
     , smvcst (land_points)                                            & 
!                            ! IN Volumetric saturation point 
     ,RMDI                                                             &
                                  ! IN real missing data indicator
     ,XUA(0:global_row_length)                                         &
                                  ! IN Atmosphere UV longitude coords
     ,YUA(global_rows)                                                 &
                                  ! IN Atmosphere latitude coords
     ,XPA(global_row_length+1)                                         &
                                  ! IN Atmosphere longitude coords
     ,YPA(global_rows)                                                 &
                                  ! IN Atmosphere latitude coords
     ,XVA(global_row_length+1)                                         &
                                  ! IN Atmosphere longitude coords
     ,YVA(0:global_rows)                                               &
                                  ! IN Atmosphere latitude coords
     ,A_BOXAREAS(row_length, rows)                                     &
                                  !IN ATMOS gridbox areas
     ,FLANDG(row_length, rows)                           
                                  ! IN Land fraction on global field.    

     ! ********* NOT USED *********!
      REAL, INTENT(IN) ::                                              &
      r_area(global_row_length, global_rows)                           &  
                                  ! IN accumalated areas file
     ,slope(global_row_length, global_rows)                            &
                                  ! IN slopes
     ,r_inext(global_row_length, global_rows)                          &
                                  ! IN X-cordinate of downstream grid pt
     ,r_jnext(global_row_length, global_rows)                          &
                                  ! IN Y-cordinate of downstream grid pt
     ,r_land(global_row_length, global_rows)                           &
                               ! IN land/river depends on value of A_thresh
     ,flowobs1(global_row_length, global_rows)                         &
                               ! IN initialisation for flows
     ,substore(global_row_length, global_rows)                         &
                               ! IN routing subsurface store
     ,surfstore(global_row_length, global_rows)                        &
                               ! IN routing surface store 
     ,flowin(global_row_length, global_rows)                           &
                               ! IN surface lateral inflow
     ,bflowin (global_row_length, global_rows)                         &     
                               ! IN subsurface lateral inflow
     ,delta_phi             
                               ! IN RCM gridsize (radians)

     ! ********* NOT USED *********!

      REAL, INTENT(IN) ::                                              &
      TRIVDIR(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                               ! IN river direction
     ,TRIVSEQ(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                               ! IN river sequence
     ,RIV_VEL                                                          &
                                  ! IN river velocity
     ,RIV_MCOEF                                                        &
                                  ! IN meandering coefficient
     ,RIV_STEP 
                                  ! IN river timestep (secs) 

      REAL, INTENT(INOUT) :: TWATSTOR(RIVER_ROW_LENGTH, RIVER_ROWS)     
      ! water store(Kg)

      LOGICAL, INTENT(IN) ::                                           &
     & INVERT_ATMOS               ! IN True if ATMOS fields are S->N
!                                 ! for regridding runoff from ATMOS.

      REAL, INTENT(OUT) ::                                             &
      RIVEROUT_ATMOS(row_length,rows)                                  &
           ! river flow out from each  gridbox(KG/m2/S)

     ,BOX_OUTFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)                        &
           ! gridbox outflow river grid (Kg/s)

     ,BOX_INFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)                         &
           ! gridbox runoff river grid(Kg/s)

! Declare new variables for inland basin outflow

     ,INLANDOUT_RIV(RIVER_ROW_LENGTH,RIVER_ROWS)                       &
           !TRIP OUTFLOW FROM INLAND BASINS ON TRIP GRID Kg/s

     ,INLANDOUT_ATMOS(ROW_LENGTH,ROWS)               
           ! TRIP OUTFLOW FROM  INLAND BASINS ON atmos GRID Kg/m2/s


! Logical to detect first entry to routing code

       LOGICAL, INTENT(IN) ::  FIRST_ENTRY

! local variables

      INTEGER                                                          &
       I,J,L

      INTEGER                                                          &
       info                                                            &
                                ! Return code from MPP
     , icode                   ! Return code :0 Normal Exit :>0 Error

      CHARACTER(LEN=256)                                                    &
       CMessage         ! Error message if return code >0
      CHARACTER(LEN=*) RoutineName
      Parameter ( RoutineName='RIV_INTCTL_1A')

      LOGICAL                                                          &
       INVERT_TRIP                                                     &
                               ! TRUE WHEN ROW INVERSION IS REQ
     , REGRID                                                          &
                               ! TRUE if TRIP grid different to ATMOS
     , CYCLIC_TRIP                                                     &
                               ! TRUE WHEN THE TRIP MODEL HAS CYCLIC
     , GLOBAL_TRIP            ! TRUE WHEN TRIP GRID SURFACE IS SPHER
       PARAMETER(INVERT_TRIP=.FALSE.,CYCLIC_TRIP=.TRUE.,               &
      GLOBAL_TRIP=.TRUE.,REGRID=.TRUE.)

!      REAL RMDI_TRIP
!      PARAMETER(RMDI_TRIP=-999)

     REAL, INTENT(OUT) ::                                             &
       sthu_dsm(land_points)                                          & 
                             ! Frozen soil moisture content of 
!                            !    bottom layer as a frac of saturation 
     , smcl_dsm(land_points)                                           
                             ! Total soil moisture contents 
!                            !      of bottom layer (kg/m2). 

      REAL, INTENT(IN) ::                                              &
       ACC_LAKE_EVAP(row_length,rows)                                   
!                                    ! Accumulated lake evap over 
!                                    ! river routing timestep (Kg/m2) 

     REAL                                                              &
      SMVCST_G(row_length,rows)                                        & 
!                                    ! SMVCST on global grid      
     ,SMVCWT_G(row_length,rows)                                        & 
!                                    ! SMVCWT on global grid      
     ,SMCL_DSM_G(row_length,rows)                                      & 
!                                    ! SMCL on DSM_LEVELS LAYER 
!                                    ! on global grid      
     ,STHU_DSM_G(row_length,rows)                                      &
!                                    ! STHU on DSM_LEVELS LAYER 
!                                    ! on global grid  
     ,SURF_RUNOFFIN(row_length,rows)                                   &
                                     ! TOTAL RATE OF RUNOFF (KG/M2/S)
     ,SUB_RUNOFFIN(row_length,rows)                                    
!                                    ! TOTAL RATE OF RUNOFF (KG/M2/S) 

      REAL                                                             &
       ACC_LAKE_EVAP_AVG                                               &  
       ! Global daily accumulated total lake evaporation 
      , TOT_LANDAREA                                                    
       ! Total land area with unfrozen soil  moisture > wilting point 

      REAL, DIMENSION(:,:), SAVE, ALLOCATABLE :: land_frac
      ! land mask for global regridding weights calculation

      INTEGER :: my_comm, ierror
      REAL :: acc_lake_evap_avg_g = 0, tot_landarea_g = 0
      REAL :: sum_r(row_length, rows, 2), sum_all(2) 

      LOGICAL, SAVE :: first = .TRUE. 
      ! use this rather than first entry because of continuation runs
      ! and the need to recalculate concerns and gather land fractions
      ! on first fun
      
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      
! gather the TRIP variables to PE0 and call the TRIP river routing
! 1. Gather coupling fields from distributed processors onto a single
!    processor for input to river routing routines.
      
      IF (lhook) CALL dr_hook('RIV_INTCTL_1A',zhook_in,zhook_handle)
      info = 0
      
      sum_r = 0.
      sum_all = 0.

!********************************************************************
      DO J=1,rows
        DO I=1,row_length
          SURF_RUNOFFIN(i,j) = 0.0
          SUB_RUNOFFIN(i,j) = 0.0 
        ENDDO
      ENDDO
      DO J=1,rows
        DO I=1,row_length
          SMVCST_G(i,j) = 0.0 
          SMVCWT_G(i,j) = 0.0 
          SMCL_DSM_G(i,j) = 0.0 
          STHU_DSM_G(i,j) = 0.0
        ENDDO
      ENDDO
      
      ! Copy land points output back to full fields array.
      DO l = 1, land_points
        j=(land_index(l)-1)/row_length + 1
        i=land_index(l) - (j-1)*row_length
        SURF_RUNOFFIN(i,j) = TOT_SURF_RUNOFF(L)
        SUB_RUNOFFIN(i,j) = TOT_SUB_RUNOFF(L)
        SMVCST_G(i,j) = SMVCST(L) 
        SMVCWT_G(i,j) = SMVCWT(L) 
        SMCL_DSM_G(i,j) = SMCL_DSM(L) 
        STHU_DSM_G(i,j) = STHU_DSM(L) 
      ENDDO
      
      
      IF(first) THEN
        ALLOCATE(land_frac(global_row_length, global_rows))

! DEPENDS ON: all_gather_field
        CALL ALL_GATHER_FIELD(flandg, land_frac,                  &
             lasize(1,fld_type_p,halo_type_no_halo),              &
             lasize(2,fld_type_p,halo_type_no_halo),              &
             glsize(1,fld_type_p), glsize(2,fld_type_p),          &
             fld_type_p,halo_type_no_halo,                        &
             GC_ALL_PROC_GROUP,info,cmessage)
        IF(info /= 0) THEN      ! Check return code
          CMESSAGE='ATMPHB2 : ERROR in gather of flandg'
          ICODE=5
          
          Call Ereport(RoutineName,icode,Cmessage)
        END IF
        
      END IF
      
      ! first get number of river grid points on this proc
      
      ! Calculate the water required to be removed from the soil to 
      ! take account of lake evaporation
      
      ! Calculate total land mass where the soil is above the wilting point: 
      tot_landarea = 0.0      
      
      ! Calculate the global total lake evaporation accumulated over a day: 
      
      acc_lake_evap_avg = 0.0
      
      DO j=1, rows
        DO i=1, row_length
          
          IF(sthu_dsm_g(i,j)*smvcst_g(i,j) .GT. smvcwt_g(i,j))            &
                sum_r(i,j,1) =  a_boxareas(i,j)*flandg(i,j)
        
          sum_r(i,j,2)= acc_lake_evap(i,j)*a_boxareas(i,j) * flandg(i,j)
          
        END DO
      END DO
      
      ! calculate global landarea and acc_lake_evap_avg sums 
      CALL global_2d_sums(sum_r, row_length, rows, 0, 0, 2, sum_all)   

      tot_landarea_g = sum_all(1) 
      acc_lake_evap_avg_g = sum_all(2) 
  
      acc_lake_evap_avg = acc_lake_evap_avg_g/tot_landarea_g
      
      ! Apply water correction: 
      
      DO j=1, rows
        DO i=1, row_length
          IF(sthu_dsm_g(i,j)*smvcst_g(i,j) .GT. smvcwt_g(i,j)) THEN
            
              smcl_dsm_g(i,j) = smcl_dsm_g(i,j) - acc_lake_evap_avg
              sthu_dsm_g(i,j) = sthu_dsm_g(i,j) - acc_lake_evap_avg     &
                    / (smvcst_g(i,j)*rho_water*dzsoil(dsm_levels))
            
            IF(sthu_dsm_g(i,j) .LT. 0.0)THEN
              write(6,'(A)')'*warning layer 4 unfr soil moisture frac<0' 
              write(6,'(2(F8.2))')STHU_DSM_g(i,j),STHU_DSM_g(i,j)         & 
                   +ACC_LAKE_EVAP_AVG 
            END IF
          
            IF(smcl_dsm_g(i,j) .LT. 0.0)THEN
              write(6,'(A)')'*warning layer 4 unfrozen soil moisture frac<0' 
              write(6,'(2(F8.2))') smcl_DSM_g(i,j), smcl_DSM_g(i,j)       & 
                   +ACC_LAKE_EVAP_AVG 
            END IF
          END IF
        END DO
      END DO
      
      
      
      CALL RIV_ROUT_1A(global_row_length, global_rows, surf_runoffin,   &
            sub_runoffin, row_length, rows, invert_atmos, xua, yva,     &
            river_row_length, river_rows, flandg, rmdi, regrid, riv_vel,&
            riv_mcoef, riv_step, cyclic_trip, a_boxareas, global_trip,  &
            invert_trip, twatstor, trivdir, trivseq,                    &
            box_outflow, box_inflow, riverout_atmos, inlandout_riv,     &
            inlandout_atmos, global_river_row_length, global_river_rows,&
            land_frac, first, cmessage)
      
      where(flandg == 1.0) riverout_atmos = 0.
      

      ! Copy full fields array back to land points: 
      DO l = 1, land_points 
        j=(land_index(l)-1)/row_length + 1 
        i=land_index(l) - (j-1)*row_length 
        SMCL_DSM(L) = SMCL_DSM_G(i,j) 
        STHU_DSM(L) = STHU_DSM_G(i,j) 
      END DO

      first = .FALSE.
      
      IF (lhook) CALL dr_hook('RIV_INTCTL_1A',zhook_out,zhook_handle)
      RETURN
    END SUBROUTINE RIV_INTCTL_1A
END MODULE riv_intctl_mod_1A
