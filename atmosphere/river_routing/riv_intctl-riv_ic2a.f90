! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE riv_intctl_mod_2A

USE riv_rout_mod_2A, ONLY: riv_rout_2A

IMPLICIT NONE

CONTAINS

       SUBROUTINE RIV_INTCTL_2A(                                        &
       XPA, XUA, XVA, YPA, YUA, YVA,                                    &
       G_P_FIELD, G_R_FIELD, N_PROC, ME, RMDI,                          &
       GATHER_PE_TRIP,LAND_POINTS,LAND_INDEX,                           &
       INVERT_ATMOS, ROW_LENGTH, ROWS,                                  &
       GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                                  &
       RIVER_ROW_LENGTH, RIVER_ROWS,                                    &
       GLOBAL_RIVER_ROW_LENGTH, GLOBAL_RIVER_ROWS,                      &
       FLANDG, RIV_STEP, RIV_VEL, RIV_MCOEF,                            &
       TRIVDIR, TRIVSEQ, TWATSTOR, A_BOXAREAS,                          &
       DELTA_PHI,FIRST,                                                 &
       r_area, slope, flowobs1,r_inext,r_jnext,r_land,                  &
       substore,surfstore,flowin,bflowin,                               &
! IN/OUT accumulated runoff
       TOT_SURF_RUNOFF, TOT_SUB_RUNOFF,                                 &
! OUT
       BOX_OUTFLOW, BOX_INFLOW, RIVEROUT_ATMOS,                         &
! Optional arguments from 1A subroutine needed for interface checking
       inlandout_atmos,inlandout_riv,                                   & 
       dsm_levels,acc_lake_evap,smvcst,smvcwt,                          & 
       smcl_dsm,sthu_dsm) 

! Purpose:
! New Control routine for River routing for Regional Model.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!-----------------------------------------------------------------

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      IMPLICIT NONE

      INTEGER                                                           &
     & row_length                                                       &
                                 ! IN NO. OF COLUMNS IN ATMOSPHERE
     &,rows                                                             &
                                 ! IN NO. OF ROWS IN ATMOSPHERE
     &, global_row_length                                               &
                                 ! number of points on a row
     &, global_rows                                                     &
                                 ! NUMBER OF global rows
     &,land_points                                                      &
                                 ! IN number of landpoints
     &,RIVER_ROW_LENGTH                                                 &
                                 ! IN no. of columns in river grid
     &,RIVER_ROWS                                                       &
                                 ! IN no. of rows in river grid
     &,GLOBAL_RIVER_ROW_LENGTH                                          &
                                 ! IN global river row length
     &,GLOBAL_RIVER_ROWS                                                &
                                 ! IN GLOBAL river rows
     &,gather_pe_trip                                                   &
                                 ! IN pe River routing to be run on
     &, n_proc                                                          &
                                 ! IN Total number of processors
     &, me                       ! IN My processor number

      INTEGER                                                           &
     &  G_P_FIELD                                                       &
                                  ! IN size of global ATMOS field
     &, G_R_FIELD                                                       &
                                  ! IN Size of global river field
     &, land_index (land_points)  ! IN index of land to global points
      REAL                                                              &
     & TOT_SURF_RUNOFF(land_points)                                     &
                                   !IN Surf RUNOFF on land pts(KG/M2/S)
     &,TOT_SUB_RUNOFF(land_points)                                      &
                                   !IN Subsurf.RUNOFF (KG/M2/S)
     &,RMDI                                                             &
                                  ! IN real missing data indicator
     &,XUA(0:row_length)                                                &
                                  ! IN Atmosphere UV longitude coords
     &,YUA(rows)                                                        &
                                  ! IN Atmosphere latitude coords
     &,XPA(row_length+1)                                                &
                                  ! IN Atmosphere longitude coords
     &,YPA(rows)                                                        &
                                  ! IN Atmosphere latitude coords
     &,XVA(row_length+1)                                                &
                                  ! IN Atmosphere longitude coords
     &,YVA(0:rows)                                                      &
                                  ! IN Atmosphere latitude coords
     &,A_BOXAREAS(row_length,rows)                                      &
                                  !IN ATMOS gridbox areas
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
                                  ! IN Land fraction on global field.
     &,DELTA_PHI                  ! RCM gridsize (radians)


      REAL                                                              &
     & TRIVDIR(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                               !IN river direction
     &,TRIVSEQ(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                               !IN river sequence
     &,TWATSTOR(RIVER_ROW_LENGTH, RIVER_ROWS)                           &
                                               !IN/OUT water store(Kg)
     &,RIV_VEL                                                          &
                                  ! IN river velocity
     &,RIV_MCOEF                                                        &
                                  ! IN meandering coefficient
     &,RIV_STEP                   ! IN river timestep (secs)

      LOGICAL                                                           &
     & INVERT_ATMOS               ! IN True if ATMOS fields are S->N
!                                 ! for regridding runoff from ATMOS.

      REAL                                                              &
     & RIVEROUT_ATMOS(row_length,rows)                                  &
                                      ! OUT river flow out from each
!                           ! gridbox(KG/m2/S)
     &,BOX_OUTFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)                        &
                                                 ! OUT gridbox outflow
!                                ! (river gridKg/s)
     &,BOX_INFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)   ! OUT gridbox runoff
!                                ! river grid(Kg/s)


! ancillary variables for grid-to-grid model

       REAL                                                             &
     & R_AREA(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                           &
                       !ACCUMULATED AREAS FILE
     & R_INEXT(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                          &
                       ! X-COORDINATE OF DOWNSTREAM GRID PT
     & R_JNEXT(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                          &
                       ! Y-COORDINATE OF DOWNSTREAM GRID PT
     & SLOPE(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                            &
                       ! SLOPES (NOT USED YET)
     & FLOWOBS1(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                         &
                       ! OPTIONAL INITIALISATION FOR FLOWS
     & R_LAND(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)
                       !LAND/RIVER DEPENDS ON VALUE OF A_THRESH

! PROGNOSTIC VARIABLES FOR GRID-TO-GRID MODEL

       REAL                                                             &
     &  SUBSTORE(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)                         &
                       ! ROUTING SUB_SURFACE STORE (MM)
     & ,SURFSTORE(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)                        &
                       ! ROUTING SURFACE STORE (MM)
     & ,FLOWIN(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)                           &
                       !SURFACE LATERAL INFLOW (MM)
     & ,BFLOWIN(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)
                       ! SUB-SURFACE LATERAL INFLOW (MM)

        LOGICAL FIRST    ! First call to river routing ? (T/F)

! Start of optional arguments from 1A subroutine 
      INTEGER, OPTIONAL :: dsm_levels
      REAL,    OPTIONAL :: inlandout_riv (river_row_length,river_rows)
      REAL,    OPTIONAL :: inlandout_atmos (row_length,rows) 
      REAL,    OPTIONAL :: acc_lake_evap (row_length,rows)
      REAL,    OPTIONAL :: smvcwt (row_length,rows)  
      REAL,    OPTIONAL :: smvcst (row_length,rows) 
      REAL,    OPTIONAL :: sthu_dsm (land_points)
      REAL,    OPTIONAL :: smcl_dsm (land_points)
! End of optional arguments from 1A subroutine 
                              
      INTEGER                                                           &
     & I,J,L
      REAL                                                              &
     & gather_TOT_RUNOFFIN(G_P_FIELD) ! TOTAL RATE OF RUNOFF (KG/M2/S)

      INTEGER                                                           &
     &  info                                                            &
                                ! Return code from MPP
     &, icode                   ! Return code :0 Normal Exit :>0 Error

      CHARACTER(LEN=256)                                                     &
     & CMessage         ! Error message if return code >0

      LOGICAL                                                           &
     & INVERT_TRIP                                                      &
                                ! TRUE WHEN ROW INVERSION IS REQUIRED
     &, REGRID                                                          &
                                ! TRUE if TRIP grid different to ATMOS
     &, CYCLIC_TRIP                                                     &
                                ! TRUE WHEN THE TRIP MODEL HAS CYCLIC
     &, GLOBAL_TRIP             ! TRUE WHEN TRIP GRID SURFACE IS SPHER
       PARAMETER(INVERT_TRIP=.FALSE.,CYCLIC_TRIP=.TRUE.,                &
     & GLOBAL_TRIP=.TRUE.,REGRID=.TRUE.)

!      REAL RMDI_TRIP
!      PARAMETER(RMDI_TRIP=-999)

      integer iarea(row_length,rows)
      integer inext(row_length,rows)
      integer jnext(row_length,rows)
      integer land(row_length,rows)

      REAL                                                              &
     & SURF_RUNOFFIN(row_length,rows)                                   &
                                     !IN TOTAL RATE OF RUNOFF (KG/M2/S)
     &,SUB_RUNOFFIN(row_length,rows) !IN TOTAL RATE OF RUNOFF (KG/M2/S)

! Gathering and Scattering variables:

      REAL gather_riverout_ATMOS(G_P_FIELD)
                                     ! river outflow at seapoints on
!                                    ! the ATMOS grid
      REAL gather_r_area (G_P_FIELD)
                                 ! global field of accumulated area
      REAL gather_r_inext (G_P_FIELD)
                                 ! global field of x-flow directions
      REAL gather_r_jnext (G_P_FIELD)
                                 ! global field of y-flow directions
      REAL gather_slope (G_P_FIELD)
                                 ! global field of slope
      REAL gather_flowobs1 (G_P_FIELD)
                                 ! global field initial flow values
      REAL gather_r_land (G_P_FIELD)
                                 ! global field of land-type
      REAL gather_substore (G_P_FIELD)
                                 ! global field of surface storage
      REAL gather_surfstore (G_P_FIELD)
                                 ! global field sub-surface storage
      REAL gather_flowin (G_P_FIELD)
                                 ! global field of flowin
      REAL gather_bflowin (G_P_FIELD)
                                 ! global field of bflowin

      REAL GATHER_SURF_RUNOFFIN(G_P_FIELD)
                                   ! field for gather runoffin to pe0
      REAL GATHER_SUB_RUNOFFIN(G_P_FIELD)
                                   ! field for gather runoffin to pe0
      REAL gather_flandg(G_P_FIELD)
                                   ! field for gather to land/sea mask

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! gather the TRIP variables to PE0 and call the TRIP river routing
! 1. Gather coupling fields from distributed processors onto a single
!    processor for input to river routing routines.

      IF (lhook) CALL dr_hook('RIV_INTCTL_2A',zhook_in,zhook_handle)
          info = 0

          DO J=1,rows
           DO I=1,row_length
             SURF_RUNOFFIN(i,j) = 0.0
           ENDDO
          ENDDO
          DO J=1,rows
           DO I=1,row_length
             SUB_RUNOFFIN(i,j) = 0.0
           ENDDO
          ENDDO

! Copy land points output back to full fields array.
          DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            SURF_RUNOFFIN(i,j) = TOT_SURF_RUNOFF(L)
            SUB_RUNOFFIN(i,j) = TOT_SUB_RUNOFF(L)
          ENDDO


! DEPENDS ON: gather_field
          CALL GATHER_FIELD(SURF_RUNOFFIN,GATHER_SURF_RUNOFFIN,         &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of SURF_RUNOFFIN'
            ICODE=1
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF


! DEPENDS ON: gather_field
          CALL GATHER_FIELD(SUB_RUNOFFIN,GATHER_SUB_RUNOFFIN,           &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of SUB_RUNOFFIN'
            ICODE=2
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(flandg,gather_flandg,                       &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of tile_frac_m1'
            ICODE=3
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(R_AREA,GATHER_R_AREA,                       &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of R_AREA'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(SLOPE,GATHER_SLOPE,                         &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of SLOPE'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(FLOWOBS1,GATHER_FLOWOBS1,                   &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of R_FLOWOBS1'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(R_INEXT,GATHER_R_INEXT,                     &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of R_INEXT'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(R_JNEXT,GATHER_R_JNEXT,                     &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of R_JNEXT'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(R_LAND,GATHER_R_LAND,                       &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of R_LAND'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(SUBSTORE,GATHER_SUBSTORE,                   &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of SUBSTORE'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(SURFSTORE,GATHER_SURFSTORE,                 &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of SURFSTORE'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(FLOWIN,GATHER_FLOWIN,                       &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of FLOWIN'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(BFLOWIN,GATHER_BFLOWIN,                     &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of BFLOWIN'
            ICODE=5
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF




! Set river routing to run on the last PE
          IF(ME == GATHER_PE_TRIP) THEN

! Convert prognostics from real to integer before calling
! river routing
         DO L = 1, global_rows*global_row_length
           j=(L-1)/global_row_length + 1
           i= L - (j-1)*global_row_length
           iarea(i,j) = nint( gather_r_area(L) )
           land(i,j)  = nint( gather_r_land(L) )
           inext(i,j) = nint( gather_r_inext(L) )
           jnext(i,j) = nint( gather_r_jnext(L) )
         ENDDO

! Call the Grid-to-grid river routing scheme

!------------------------------------------------------------------


! Call River routing
      CALL RIV_ROUT_2A(                                                 &
     &   GATHER_SURF_RUNOFFIN, GATHER_SUB_RUNOFFIN,                     &
     &   GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                                &
     &   DELTA_PHI,FIRST,RIV_STEP,                                      &
! ancillary variables
     &   IAREA, GATHER_SLOPE, GATHER_FLOWOBS1,                          &
     &   INEXT,JNEXT,LAND,                                              &
! prognostic variables
     &   GATHER_SUBSTORE,GATHER_SURFSTORE,                              &
     &   GATHER_FLOWIN,GATHER_BFLOWIN,                                  &
     &   GATHER_RIVEROUT_ATMOS)

! Convert prognostics from INTEGER to REAL before scattering
        DO L = 1, global_rows*global_row_length
           j=(L-1)/global_row_length + 1
           i= L - (j-1)*global_row_length
           gather_r_land(L) = real( land(i,j) )
         ENDDO

! Set all total land pts to 0.0
          WHERE(gather_flandg == 1.0)gather_riverout_ATMOS=0.0

          ENDIF                         ! Single processor
!

! DEPENDS ON: scatter_field
         CALL SCATTER_FIELD(RIVEROUT_ATMOS,gather_riverout_ATMOS,       &
     &    lasize(1,fld_type_p,halo_type_no_halo),                       &
     &    lasize(2,fld_type_p,halo_type_no_halo),                       &
     &    glsize(1,fld_type_p),                                         &
     &    glsize(2,fld_type_p),                                         &
     &    fld_type_p,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of RIVEROUT_ATMOS'
            ICODE=101
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF


! DEPENDS ON: scatter_field
         CALL SCATTER_FIELD(R_LAND,GATHER_R_LAND,                       &
     &    lasize(1,fld_type_p,halo_type_no_halo),                       &
     &    lasize(2,fld_type_p,halo_type_no_halo),                       &
     &    glsize(1,fld_type_p),                                         &
     &    glsize(2,fld_type_p),                                         &
     &    fld_type_p,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of R_LAND'
            ICODE=101
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF


! DEPENDS ON: scatter_field
         CALL SCATTER_FIELD(SUBSTORE,GATHER_SUBSTORE,                   &
     &    lasize(1,fld_type_p,halo_type_no_halo),                       &
     &    lasize(2,fld_type_p,halo_type_no_halo),                       &
     &    glsize(1,fld_type_p),                                         &
     &    glsize(2,fld_type_p),                                         &
     &    fld_type_p,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of SUBSTORE'
            ICODE=101
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: scatter_field
         CALL SCATTER_FIELD(SURFSTORE,GATHER_SURFSTORE,                 &
     &    lasize(1,fld_type_p,halo_type_no_halo),                       &
     &    lasize(2,fld_type_p,halo_type_no_halo),                       &
     &    glsize(1,fld_type_p),                                         &
     &    glsize(2,fld_type_p),                                         &
     &    fld_type_p,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of SURFSTORE'
            ICODE=101
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                       &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: scatter_field
         CALL SCATTER_FIELD(FLOWIN,GATHER_FLOWIN,                       &
     &    lasize(1,fld_type_p,halo_type_no_halo),                       &
     &    lasize(2,fld_type_p,halo_type_no_halo),                       &
     &    glsize(1,fld_type_p),                                         &
     &    glsize(2,fld_type_p),                                         &
     &    fld_type_p,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of FLOWIN'
            ICODE=101
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF

! DEPENDS ON: scatter_field
         CALL SCATTER_FIELD(BFLOWIN,GATHER_BFLOWIN,                     &
     &    lasize(1,fld_type_p,halo_type_no_halo),                       &
     &    lasize(2,fld_type_p,halo_type_no_halo),                       &
     &    glsize(1,fld_type_p),                                         &
     &    glsize(2,fld_type_p),                                         &
     &    fld_type_p,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of BFLOWIN'
            ICODE=101
            write(6,*) CMESSAGE,ICODE

            CALL EREPORT('RIV_INTCTL_2A', ICODE,                        &
     &       CMESSAGE)

          ENDIF



      IF (lhook) CALL dr_hook('RIV_INTCTL_2A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE RIV_INTCTL_2A
END MODULE riv_intctl_mod_2A
