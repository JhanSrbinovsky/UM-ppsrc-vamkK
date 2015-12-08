! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Initialises arrays used in area weighed horizontal interpolation.
      SUBROUTINE h_int_init_aw(icof,idim,max_field_out                  &
      ,                        max_rows_in,max_rows_out                 &
      ,                        row_length_in,row_length_out             &
      ,                        global,FixHd_in,FixHd_out                &
      ,                        RealHd_in,RealHd_out,aw_area_box         &
      ,                        aw_index_targ_lhs,aw_index_targ_top      &
      ,                        bl_index_b_l,bl_index_b_r                &
      ,                        aw_colat_t,aw_long_l                     &
      ,                        weight_t_r,weight_b_r                    &
      ,                        weight_t_l,weight_b_l)
!
! Subroutine Interface:

      USE Ereport_Mod, Only :                                           &
     &    Ereport

      USE Rcf_Grid_Type_Mod, Only :                                     &
     &    Input_Grid,                                                   &
     &    Output_Grid

      USE Rcf_HeadAddress_Mod, ONLY:                                    &
          FH_HorizGrid,     FH_GridStagger,                             &
          RC_LongSpacing,   RC_LatSpacing,                              &
          RC_FirstLat,      RC_FirstLong,                               &
          RC_PoleLat,       RC_PoleLong,                                &
          FH_GridStagger_A, FH_GridStagger_C, FH_GridStagger_Endgame,   &
          FH_RowDepCStart

      USE h_int_co_mod, ONLY: h_int_co
      IMPLICIT NONE
!
! Description:
!   Initialises arrays used in horizontal interpolation.
!   This replaces routine SETWTS1 (A Dickinson) whose function it
!   incorporates.
!
! Method:
!   Sets up gather index and weight arrays for later call to H_INT_BL h
!   Also sets up rotation coefficients for use in ROTATE with rotated
!   grids.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v8 programming standards.
!
! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER  ::  icof             !Second dimension of coefficients ary
      INTEGER  ::  idim             !Second dimension of index arrays
      INTEGER  ::  max_field_out      !No of P pts on target grid
      INTEGER  ::  max_rows_in        !No of P rows on source grid
      INTEGER  ::  max_rows_out       !No of P rows on target grid
      INTEGER  ::  row_length_in    !No of pts per row on source grid
      INTEGER  ::  row_length_out   !No of pts per row on target grid
                                    ! grid
      LOGICAL  ::  global           !T= Global; F= LAM.

!   Array  arguments with intent(in):
      INTEGER  ::  FixHd_in(*)      !Fixed length header for source grid
      INTEGER  ::  FixHd_out(*)     !Fixed length header for target grid
      REAL     ::  RealHd_in(*)     !Real constants from source grid
      REAL     ::  RealHd_out(*)    !Real constants from target grid

!   Array  arguments with intent(Out):

      INTEGER  ::  aw_index_targ_lhs(row_length_out+1,idim)
                                    !Index of source box overlapping
                                    !lhs of target grid-box
      INTEGER  ::  aw_index_targ_top(max_rows_out+1,idim)
                                    !Index of source box overlapping
                                    !top of target grid-box
      INTEGER  ::  bl_index_b_l(max_field_out,idim)
                                    !Gather index for bottom l.h.c of
                                    !source grid box. 1=P-pts; 2=UV-pts
      INTEGER  ::  bl_index_b_r(max_field_out,idim)
                                    !Gather index for bottom r.h.c of
                                    !source grid box. 1=P-pts; 2=UV-pts
      REAL     ::  aw_area_box(idim) !area of grid box in sq units of
                                    !  source grid
      REAL     ::  aw_colat_t(max_rows_out+1,idim)
                                    !Colatitude of top of target grd-box
                                    ! (in units of DELTA_LAT_SRCE)
      REAL     ::  aw_long_l(row_length_out+1,idim)
                                    !Left longitude of target grid-box
                                    ! (in units of DELTA_LONG_SRCE)
      REAL     ::  weight_t_r(max_field_out,idim) ! Weights for bilinear
      REAL     ::  weight_b_r(max_field_out,idim) !\horizontal interpolatn
      REAL     ::  weight_t_l(max_field_out,idim) !/ 1=P-pts; 2=U-pts;
      REAL     ::  weight_b_l(max_field_out,idim) ! 3=V-pts;4=zonal mean


      INTEGER  ::  ErrorStatus          ! Error flag (0 = OK)

! Local parameters:
      Character (Len=*), Parameter :: RoutineName='HINTIAW'

! Local scalars:
      LOGICAL   :: L_endgame_out
      LOGICAL   :: L_endgame_in
      INTEGER   :: i                !\                                 .
      INTEGER   :: ij               ! \ Loop
      INTEGER   :: j                ! /variables
      INTEGER   :: k                !/
      REAL      :: delta_lat_source !\                                 .
      REAL      :: delta_lat_target ! \Grid spacing
      REAL      :: delta_lon_source ! /
      REAL      :: delta_lon_target !/
      REAL      :: npole_lat_source !\                                 .
      REAL      :: npole_lat_target ! \North pole coordinates
      REAL      :: npole_lon_source ! /
      REAL      :: npole_lon_target !/
      REAL      :: start_lat_source !\                                 .
      REAL      :: start_lat_target ! \Coordinates of first data point
      REAL      :: start_lon_source ! /
      REAL      :: start_lon_target !/
      LOGICAL   :: cyclic           !T= Data cyclic
      CHARACTER (LEN=80) :: Cmessage

! Local dynamic arrays:
! i_grid_in/out signifies whether the grid is at the pole:
!   * 1 at pole
!   * 2 not at pole.
! icof=1 is the p grid
! icof=2 is the u grid
! icof=3 is the v grid
      INTEGER   :: i_grid_in(icof)
      INTEGER   :: i_grid_out(icof)
      REAL      :: p_offset_in, p_offset_out
      REAL      :: d_lat_in(icof)
      REAL      :: d_lon_in(icof)
      REAL      :: d_lat_out(icof)
      REAL      :: d_lon_out(icof)
      REAL      :: lambda_in(row_length_in)
                                  ! Latitude coords of source p-grid
      REAL      :: lambda_out(max_field_out)
                                    ! Latitude coords of target grid
      REAL      :: phi_in(max_rows_in)
                                    ! Longitude coords of source p-grid
      REAL      :: phi_out(max_field_out)
                                    ! Longitude coords of target grid
!- End of header

! Set Endgame flags if running Endgame
      L_endgame_in  = FixHd_in(FH_GridStagger) == FH_GridStagger_Endgame
      L_endgame_out = FixHd_out(FH_GridStagger) == FH_GridStagger_Endgame

! Abort if trying variable resolution
      IF (FixHd_in(FH_RowDepCStart) > 0) THEN
        Cmessage = 'Area-weighted interpolation not compatibile with '     &
                    //'variable resolution input grid'
        ErrorStatus = 7
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF
      
      IF (FixHd_out(FH_RowDepCStart) > 0) THEN
        Cmessage = 'Area-weighted interpolation not compatibile with '     &
                    //'variable resolution output grid'
        ErrorStatus = 8
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF

      IF ( L_endgame_in ) THEN
        p_offset_in   = -0.5
      ELSE
        p_offset_in   = -1.0
      END IF

      IF ( L_endgame_out ) THEN
        p_offset_out   = -0.5
      ELSE
        p_offset_out   = -1.0
      END IF
      

! 1: Test if requested horizontal interpolation is sensible

! 1.1: Hemispheric or LAM -> global not allowed
      IF (FixHd_out( FH_HorizGrid ) == 0 .AND.                                &
          FixHd_in( FH_HorizGrid ) >  0) THEN
        WRITE(6,'(A)') ' *ERROR* Trying to interpolate from a hemispheric' // &
                       ' or LAM to a global domain'
        Cmessage = 'Cannot interpolate from hemisphere or LAM to global'
        ErrorStatus = 2
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF

! 1.2: LAM -> hemispheric not allowed
      IF (FixHd_out( FH_HorizGrid ) <  3 .AND.                             &
          FixHd_in( FH_HorizGrid ) >  2) THEN
        WRITE(6,'(A)') ' *ERROR* Trying to interpolate from a limited ' // &
                  'area domain to a global or hemispheric domain'
        Cmessage = 'Cannot interpolate from LAM to global or hemisphere'
        ErrorStatus = 3
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF

! 2: Initialise local constants

! 2.1: Grid spacing
      delta_lat_source=RealHd_in( RC_LatSpacing )
      delta_lat_target=RealHd_out( RC_LatSpacing )
      delta_lon_source=RealHd_in( RC_LongSpacing )
      delta_lon_target=RealHd_out( RC_LongSpacing )

! 2.2: Coordinates of north pole on grid
      npole_lat_source=RealHd_in( RC_PoleLat )
      npole_lat_target=RealHd_out( RC_PoleLat )
      npole_lon_source=RealHd_in( RC_PoleLong )
      npole_lon_target=RealHd_out( RC_PoleLong )

! 2.3: Coordinates of top left hand  p-point on grid
      start_lat_source=RealHd_in( RC_FirstLat )
      start_lat_target=RealHd_out( RC_FirstLat )
      start_lon_source=RealHd_in( RC_FirstLong )
      start_lon_target=RealHd_out( RC_FirstLong )

  ! 2.4: Logical to indicate if input data cyclic
      cyclic=FixHd_in( FH_HorizGrid ) <  3

  ! 2.5: Initialise i_grid_in and IGRID_OUT to 1 ie all points at pole

      i_grid_in (:) = 1
      i_grid_out(:) = 1

! 3: Weights and indices for P points:

  ! 3.1: If source or target grids have different poles
  !      abort with error message

      IF(npole_lat_source /= npole_lat_target.AND.                      &
         npole_lon_source /= npole_lon_target)THEN

        WRITE(6,'(A)') 'Source and target grids have different poles'
        WRITE(6,'(A,A)') 'Reconfigure onto a grid with the same pole as'      &
      ,                  'target grid'
        WRITE(6,'(A)') 'before attempting area weighted interpolation'

        Cmessage = 'Grids have different poles!'
        ErrorStatus = 4
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      ENDIF



! 4: Weights and indices for U and V points:

  ! 4.1: Calculate offsets for source winds
      IF(FixHd_in(fh_gridstagger) == fh_gridstagger_a ) THEN

    ! 4.1.1: Source winds on A grid
        d_lat_in(1)=0.0
        d_lon_in(1)=0.0
        d_lat_in(2)=0.0
        d_lon_in(2)=0.0
        d_lat_in(3)=0.0
        d_lon_in(3)=0.0

      ELSE IF (FixHd_in( fh_gridstagger ) == fh_gridstagger_c ) THEN

    ! 4.1.2: Source winds on C grid (New Dynamics)
        d_lat_in(1)=0.0
        d_lon_in(1)=0.0
        d_lat_in(2)=0.0
        d_lon_in(2)=0.5
        d_lat_in(3)=0.5
        d_lon_in(3)=0.0
        i_grid_in(1)=1
        i_grid_in(2)=1
        i_grid_in(3)=2

      ELSE IF (FixHd_in( fh_gridstagger ) == fh_gridstagger_endgame ) THEN

    ! 4.1.2: Source winds on C grid (Endgame)
    ! Note that the zero point in Endgame is NOT the first p-point
        d_lat_in(1)=0.5
        d_lon_in(1)=0.5
        d_lat_in(2)=0.5
        d_lon_in(2)=0.0
        d_lat_in(3)=0.0
        d_lon_in(3)=0.5
        i_grid_in(1)=2
        i_grid_in(2)=2
        i_grid_in(3)=1

      ELSE

    ! 4.1.3: Source winds on B grid
        d_lat_in(1)=0.0
        d_lon_in(1)=0.0
        d_lat_in(2)=0.5
        d_lon_in(2)=0.5
        d_lat_in(3)=0.5
        d_lon_in(3)=0.5
        i_grid_in(1)=1
        i_grid_in(2)=2
        i_grid_in(3)=2
      END IF

  ! 4.2: Calculate offsets for target winds

      IF (FixHd_out( fh_gridstagger ) == fh_gridstagger_c ) THEN

    ! 4.2.1: Target winds on C grid
        d_lat_out(1)=0.0
        d_lon_out(1)=0.0
        d_lat_out(2)=0.0
        d_lon_out(2)=0.5
        d_lat_out(3)=0.5
        d_lon_out(3)=0.0
        i_grid_out(1)=1
        i_grid_out(2)=1
        i_grid_out(3)=2

      ELSE IF (FixHd_out( fh_gridstagger ) == fh_gridstagger_endgame ) THEN

    ! 4.2.2: Target winds on C grid (Endgame)
    ! Note that the zero point in Endgame is NOT the first p-point
        d_lat_out(1)=0.5
        d_lon_out(1)=0.5
        d_lat_out(2)=0.5
        d_lon_out(2)=0.0
        d_lat_out(3)=0.0
        d_lon_out(3)=0.5
        i_grid_out(1)=2
        i_grid_out(2)=2
        i_grid_out(3)=1

      ELSE

    ! 4.2.3: Target winds on B grid
        d_lat_out(1)=0.0
        d_lon_out(1)=0.0
        d_lat_out(2)=0.5
        d_lon_out(2)=0.5
        d_lat_out(3)=0.5
        d_lon_out(3)=0.5
        i_grid_out(1)=1
        i_grid_out(2)=2
        i_grid_out(3)=2

      ENDIF

  ! 3.2: Calculate area weighted indices for
  !      interpolating from the source grid onto the target grid
  ! P grid

! DEPENDS ON: box_bnd
      CALL box_bnd(aw_index_targ_lhs(1,1),aw_long_l(1,1)                &
      ,            aw_index_targ_top(1,1),aw_colat_t(1,1)               &
      ,            aw_area_box(1)                                       &
      ,            Output_Grid % Glob_p_row_length                      &
      ,            Output_Grid % Glob_p_rows                            &
      ,            Input_Grid % Glob_p_row_length                       &
      ,            Input_Grid % Glob_p_rows                             &
      ,            delta_lon_target,delta_lat_target                    &
      ,            start_lon_target+d_lon_out(1)*delta_lon_target       &
      ,            start_lat_target+d_lat_out(1)*delta_lat_target       &
      ,            delta_lon_source, delta_lat_source                   &
      ,            start_lon_source+d_lon_in(1)*delta_lon_source        &
      ,            start_lat_source+d_lat_in(1)*delta_lat_source        &
      ,            i_grid_out(1),i_grid_in(1),global)

  ! 4.3: Calculate area weighted indices for

  !      interpolating from the source grid onto the target grid

  ! This should be fixed for all C grid stuff now.

  ! U grid
! DEPENDS ON: box_bnd
        CALL box_bnd(aw_index_targ_lhs(1,2),aw_long_l(1,2)              &
      ,              aw_index_targ_top(1,2),aw_colat_t(1,2)             &
      ,              aw_area_box(2)                                     &
      ,              Output_Grid % Glob_u_row_length                    &
      ,              Output_Grid % Glob_u_rows                          &
      ,              Input_Grid % Glob_u_row_length                     &
      ,              Input_Grid % Glob_u_rows                           &
      ,              delta_lon_target,delta_lat_target                  &
      ,              start_lon_target+d_lon_out(2)*delta_lon_target     &
      ,              start_lat_target+d_lat_out(2)*delta_lat_target     &
      ,              delta_lon_source,delta_lat_source                  &
      ,              start_lon_source+d_lon_in(2)*delta_lon_source      &
      ,              start_lat_source+d_lat_in(2)*delta_lat_source      &
      ,              i_grid_out(2),i_grid_in(2),global)

  ! V Grid
! DEPENDS ON: box_bnd
        CALL box_bnd(aw_index_targ_lhs(1,3),aw_long_l(1,3)              &
      ,              aw_index_targ_top(1,3),aw_colat_t(1,3)             &
      ,              aw_area_box(3)                                     &
      ,              Output_Grid % Glob_v_row_length                    &
      ,              Output_Grid % Glob_v_rows                          &
      ,              Input_Grid % Glob_v_row_length                     &
      ,              Input_Grid % Glob_v_rows                           &
      ,              delta_lon_target,delta_lat_target                  &
      ,              start_lon_target+d_lon_out(3)*delta_lon_target     &
      ,              start_lat_target+d_lat_out(3)*delta_lat_target     &
      ,              delta_lon_source,delta_lat_source                  &
      ,              start_lon_source+d_lon_in(3)*delta_lon_source      &
      ,              start_lat_source+d_lat_in(3)*delta_lat_source      &
      ,              i_grid_out(3),i_grid_in(3),global)

! 5: Weights and indices for zonal mean P points:

! 5.2: Calculate area weighted indices for
  !      interpolating from the source grid onto the target grid
! DEPENDS ON: box_bnd
      CALL box_bnd(aw_index_targ_lhs(1,4),aw_long_l(1,4)                &
      ,            aw_index_targ_top(1,4),aw_colat_t(1,4)               &
      ,            aw_area_box(4)                                       &
      ,            1                                                    &
      ,            Output_Grid % Glob_p_rows                            &
      ,            1                                                    &
      ,            Input_Grid % Glob_p_rows                             &
      ,            delta_lon_target, delta_lat_target                   &
      ,            start_lon_target+d_lon_out(1)*delta_lon_target       &
      ,            start_lat_target+d_lat_out(1)*delta_lat_target       &
      ,            delta_lon_source, delta_lat_source                   &
      ,            start_lon_source+d_lon_in(1)*delta_lon_source        &
      ,            start_lat_source+d_lat_in(1)*delta_lat_source        &
      ,            i_grid_out(1),i_grid_in(1),global)


! 6: Weights and indices for Coastal adjustment and Integer fields

  ! 6.1: Lat and lon of target grid
      ij=0
      DO j=1,max_rows_out
        DO i=1,row_length_out
          ij=ij+1
          lambda_out(ij)=start_lon_target+delta_lon_target*(i+p_offset_out)
          phi_out(ij)=start_lat_target+delta_lat_target*(j+p_offset_out)
        END DO
      END DO

  ! 6.2: Lat and lon of source grid
      DO j=1,max_rows_in
        phi_in(j)=start_lat_source+delta_lat_source*(j+p_offset_in)
      END DO
      DO i=1,row_length_in
        lambda_in(i)=start_lon_source+delta_lon_source*(i+p_offset_in)
      END DO

  ! 6.3: Initialise Indices and weights for Bi-linear interpolation

      CALL h_int_co(bl_index_b_l(1,1),bl_index_b_r(1,1)                 &
      ,             weight_t_r(1,1),weight_b_r(1,1)                     &
      ,             weight_t_l(1,1),weight_b_l(1,1)                     &
      ,             lambda_in,phi_in,lambda_out,phi_out                 &
      ,             row_length_in,max_rows_in,max_field_out,cyclic)

      RETURN
      END SUBROUTINE h_int_init_aw
