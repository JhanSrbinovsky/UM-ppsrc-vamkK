! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: This part of the code has been removed from RAD_CTL2
!              in order to make that deck more readable. This part of
!              the code was contained in RAD_CTL2 up to UM 6.1
!
! Purpose: This routine mainly deals with the spatial degradation and
!          the use of segments
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!-----------------------------------------------------------------------
!
  subroutine prelim_lwrad(                                           &
! Parallel variables
       at_extremity, n_proc,                                         &
! Model Dimensions
       row_length,rows,model_levels,                                 &
! Model Switches
       model_domain,mt_global,L_rad_deg,l_extra_top_lw,              &
       L_subsample, L_geostationary,                                 &
! Time stepping Information
       timestep_number,a_lw_radstep,                                 &
! ancillary fields and fields needed to be kept from timestep to
! timestep
       lw_incs,                                                      &
! Satelitte Geometry
       min_view_lon, max_view_lon,                                   &
       min_view_lat, max_view_lat,                                   &
! Number of call
       j_lw,                                                         &
! Other variables
       true_latitude,true_longitude, seconds_since_midnight,         &
       rad_mask,list_lw_points,first_data_interp,                    &
       first_row,last_row,diag_row_list,diag_col_list,               &
       lw_points, OLR, lw_down, lwsea, top_absorption)

! Declare variables


    USE dynamics_grid_mod, ONLY: l_vatpoles

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    Implicit none

    Integer :: N_proc
!     Total number of processors

    Logical :: at_extremity(4)
!     Indicates if this processor is at north, south,
!     east or west of the processor grid.

    Integer, Parameter :: PNorth=1, PSouth=3
!     North processor address in the neighbour array
!     and South processor address in neighbour array


    Integer :: row_length
    Integer :: rows
    Integer :: model_levels

    Integer :: model_domain
    Integer :: mt_global

    Logical :: L_rad_deg
!     Controls the spatial degradation of E-S code

    Logical :: L_subsample
!     Flag to apply spatial subsampling (for satellite footprints)
!     in radiation.
    Logical :: L_geostationary
!     Flag to signal that a geostationary satellite is assumed.

    Logical :: L_extra_top_lw
!     Flag to use an extra top layer in radiative calculations

    Logical :: rad_mask(row_length, rows)
!     A mask which ensures a chequerboard pattern of radiation
!     calculations over the whole domain (not just one PE)

    Real :: min_view_lon
!     Minimum longitude of viewing domain
    Real :: max_view_lon
!     Maximum longitude of viewing domain
    Real :: min_view_lat
!     Minimum latitude of viewing domain
    Real :: max_view_lat
!     Maximum latitude of viewing domain

    Real :: seconds_since_midnight
    Real :: true_longitude(row_length, rows)
    Real :: true_latitude(row_length, rows)


    Real :: LW_incs(row_length, rows, 0:model_levels)
    Real :: OLR(row_length, rows)
    Real :: lw_down(row_length, rows)
    Real :: LWsea(row_length, rows)
    Real :: top_absorption(row_length, rows)

    Integer :: j_lw

    Integer :: timestep_number
    Integer :: a_lw_radstep

    Integer :: list_lw_points(row_length*rows)
    Integer :: lw_points
!     Variables to enable scatter/gather of LW (like SW):
!     creates a LIST of points to do a calculation, and
!     also facilitates segmenting
    Integer :: first_data_interp
!     The first data point, in data co-ords, which needs to be
!     interpolated on a PE (for use in interpolation routine)
    Integer :: first_row,last_row

    Integer :: diag_row_list(row_length*rows)
!     List of row indices of points where diagnostics are
!     calculated
    Integer :: diag_col_list(row_length*rows)
!     List of column indices of points where diagnostics are
!     calculated

! Local Variables

    Logical :: l_viewed

    Integer :: i,j

    Integer :: interp_index
!     Variable which alternates between 0 and 1 every other
!     radiation timestep

! External Routines

    Logical :: in_footprint

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('PRELIM_LWRAD',zhook_in,zhook_handle)
    lw_points = 0
    IF (L_rad_deg) THEN
       interp_index=MOD((timestep_number/A_LW_RADSTEP),2)
!
       If ( ( interp_index == 0 ) .EQV. rad_mask(1,1) ) Then
          first_data_interp = 1
       Else
          first_data_interp = 0
       Endif
!
! Calculations need only be performed at one point on polar
! rows, so special treatment is required in these cases.
!
       first_row = 1
       last_row  = rows
IF (.NOT. l_vatpoles) THEN
       If (model_domain.eq.mt_global) Then
          If (at_extremity(PSouth)) Then
             If (L_subsample) Then
! DEPENDS ON: in_footprint
                l_viewed = in_footprint(j_lw, .FALSE.,               &
                 L_geostationary,                                    &
                 min_view_lon,                                       &
                 max_view_lon,                                       &
                 min_view_lat,                                       &
                 max_view_lat,                                       &
                 true_latitude(1,1),                                 &
                 true_longitude(1,1),                                &
                 seconds_since_midnight)
             Else
                l_viewed=.TRUE.
             Endif
             if (l_viewed) then
               LW_POINTS = LW_POINTS + 1
               LIST_LW_points(LW_POINTS)=1
             endif
             first_row = 2
          Endif
          If (at_extremity(PNorth)) Then
             last_row = rows-1
          Endif
       Endif
END IF ! vatpoles
!
       Do j=first_row, last_row
          Do i=1, row_length
             If ( rad_mask(i,j) .EQV. ( interp_index == 0 ) ) Then
                If (L_subsample) Then
! DEPENDS ON: in_footprint
                   l_viewed = in_footprint(j_lw, .FALSE.,            &
                     L_geostationary,                                &
                     min_view_lon,max_view_lon,                      &
                     min_view_lat,max_view_lat,                      &
                     true_latitude(i,j),                             &
                     true_longitude(i,j),                            &
                     seconds_since_midnight)
                Else
                   l_viewed=.TRUE.
                Endif
                if (l_viewed) Then
                  lw_points                 = lw_points + 1
                  list_lw_points(lw_points) = i+(j-1)*row_length
                endif
             Endif
          Enddo
       Enddo
!
! Include one point on the northern polar row of a global
! domain.
!
IF (.NOT. l_vatpoles) THEN
       If (model_domain.eq.mt_global.and.                            &
                                      at_extremity(PNorth)) Then
          If (L_subsample) Then
! DEPENDS ON: in_footprint
             l_viewed = in_footprint(j_lw, .FALSE.,                  &
            L_geostationary,                                         &
            min_view_lon,max_view_lon,                               &
            min_view_lat,max_view_lat,                               &
            true_latitude(row_length,rows),                          &
            true_longitude(row_length,rows),                         &
            seconds_since_midnight)
          Else
             l_viewed=.TRUE.
          Endif
          if (l_viewed) then
            lw_points                 = lw_points + 1
            LIST_LW_points(lw_points) = (rows-1)*row_length + 1
          Endif
       endif
END IF ! vatpoles
!
! Initialize all output fields.
!
       OLR(:,:)       = 0.0
       lw_down(:,:)   = 0.0
       LWsea(:,:)     = 0.0
       LW_incs(:,:,:) = 0.0
       IF (l_extra_top_lw) Then
          top_absorption(:,:)=0.0
       Endif
!
!
    ELSE  ! SPATIAL DEGRADATION IS SWITCHED OFF
!
! Calculations need only be performed at one point on polar
! rows, so special treatment is required in these cases.
!
       first_row = 1
       last_row  = rows
IF (.NOT. l_vatpoles) THEN
       If (model_domain.eq.mt_global) Then
          If (at_extremity(PSouth)) Then
             If (L_subsample) Then
! DEPENDS ON: in_footprint
                l_viewed = in_footprint(j_lw, .FALSE.,                  &
                    L_geostationary,                                    &
                    min_view_lon,max_view_lon,                          &
                    min_view_lat,max_view_lat,                          &
                    true_latitude(1,1),                                 &
                    true_longitude(1,1),                                &
                    seconds_since_midnight)
             Else
                l_viewed=.TRUE.
             Endif
             if (l_viewed) then
                LW_POINTS = LW_POINTS + 1
                LIST_LW_points(LW_POINTS)=1
             endif
             first_row = 2
          Endif
          If (at_extremity(PNorth)) Then
             last_row = rows-1
          Endif
       Endif
END IF ! vatpoles
!
       Do j=first_row, last_row
          Do i=1, row_length
             If (L_subsample) Then
! DEPENDS ON: in_footprint
                l_viewed = in_footprint(j_lw, .FALSE.,                  &
                    L_geostationary,                                    &
                    min_view_lon,max_view_lon,                          &
                    min_view_lat,max_view_lat,                          &
                    true_latitude(i,j),                                 &
                    true_longitude(i,j),                                &
                    seconds_since_midnight)
             Else
                  l_viewed=.TRUE.
             Endif
             if (l_viewed) then
               lw_points                 = lw_points + 1
               list_lw_points(lw_points) = i+(j-1)*row_length
             endif
          Enddo
       Enddo
!
! Include one point on the northern polar row of a global
! domain.
!
IF (.NOT. l_vatpoles) THEN
       If (model_domain.eq.mt_global.and.                               &
                                      at_extremity(PNorth)) Then
          If (L_subsample) Then
! DEPENDS ON: in_footprint
             l_viewed = in_footprint(j_lw, .FALSE.,                     &
                 L_geostationary,                                       &
                 min_view_lon,max_view_lon,                             &
                 min_view_lat,max_view_lat,                             &
                 true_latitude(row_length,rows),                        &
                 true_longitude(row_length,rows),                       &
                 seconds_since_midnight)
          Else
             l_viewed=.TRUE.
          Endif
          if (l_viewed) then
            lw_points                 = lw_points + 1
            LIST_LW_points(lw_points) = (rows-1)*row_length + 1
          endif
       Endif
END IF ! vatpoles
!
    ENDIF  ! If L_rad_deg block
!
! Infer the row and column indices of the points where
! calculations are required. This is very mildly inefficient,
! but is cleaner than repeating code throughout the preceding
! block.
!
    Do j=1, lw_points
       diag_row_list(j) = (LIST_LW_points(j) - 1)/row_length + 1
       diag_col_list(j) = LIST_LW_points(j) - row_length                &
            * (diag_row_list(j) - 1)
    Enddo
    IF (lhook) CALL dr_hook('PRELIM_LWRAD',zhook_out,zhook_handle)
    RETURN
!
  end subroutine prelim_lwrad
