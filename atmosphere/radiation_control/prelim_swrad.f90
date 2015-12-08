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
  subroutine prelim_swrad(error_code,                                &
! Parallel variables
       at_extremity, n_proc,                                         &
! Model Dimensions
       row_length,rows,model_levels,                                 &
! Model Switches
       model_domain,mt_global,L_rad_deg,                             &
       L_subsample, L_geostationary,                                 &
! Timestepping information
       timestep_number,a_sw_radstep,                                 &
! Satelitte Geometry
       min_view_lon, max_view_lon,                                   &
       min_view_lat, max_view_lat,                                   &
! Number of call
       j_sw,                                                         &
! Other variables
       true_latitude,true_longitude, seconds_since_midnight,         &
       tot_daylight_points,daylight_points,day_fraction,             &
       list_daylight_points, rad_mask, switch,                       &
       first_data_interp, first_data_interp_sw,                      &
       diag_row_list, diag_row_list_sw,                              &
       diag_col_list, diag_col_list_sw )


    USE dynamics_grid_mod, ONLY: l_vatpoles

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    USE ereport_mod, ONLY : ereport
    IMPLICIT NONE

    Integer :: n_proc
!     Total number of processors
    Integer :: error_code

    Logical :: at_extremity(4)
!     Indicates if this processor is at north, south
!     east or west of the processor grid

    Integer, Parameter :: PNorth=1, PSouth=3
!     North processor address in the neighbour array
!     South processor address in the neighbour array

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

    Logical :: rad_mask(row_length, rows)
!     A mask which ensures a chequerboard pattern of radiation
!     calculations over the whole domain (not just one PE)
    Logical :: switch(row_length,rows)

    Integer :: j_sw

    Integer :: A_SW_RADSTEP
    Integer :: timestep_number


    Integer :: daylight_points
    Integer :: tot_daylight_points
!     Total number of daylight points in whole domain
    Integer :: List_daylight_points(row_length*rows)


    Integer :: first_data_interp
!     The first data point, in data co-ords, which needs to be
!     interpolated on a PE (for use in interpolation routine)
    Integer :: first_data_interp_sw
!     Saved value to use with ISCCP diagnostics

    Integer :: diag_row_list(row_length*rows)
!     List of row indices of points where diagnostics are
!     calculated
    Integer :: diag_row_list_sw(row_length*rows)
!     Saved version for ISCCP
    Integer :: diag_col_list(row_length*rows)
!     List of column indices of points where diagnostics are
!     calculated
    Integer :: diag_col_list_sw(row_length*rows)
!     Saved version for ISCCP

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

    Real :: day_fraction(row_length,rows)

! Local Variables

    Logical :: l_viewed

    Integer :: i,j

    Integer :: interp_index
!         Variable which alternates between 0 and 1 every other
!         radiation timestep

! External routines

    Logical :: in_footprint

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('PRELIM_SWRAD',zhook_in,zhook_handle)
    If (L_rad_deg) Then
!
       INTERP_INDEX=MOD((timestep_number/A_SW_RADSTEP),2)
!
       If ( (interp_index == 0) .EQV. rad_mask(1,1)) Then
          first_data_interp = 1
       Else
          first_data_interp = 0
       Endif
!
! Save interpolation parameter for later use with isccp diagnostics
!
       first_data_interp_sw = first_data_interp
!
       Do j=1, rows
          Do i=1, row_length
             switch(i,j) = ( (day_fraction(i,j) > 0.0) .AND.      &
                           (rad_mask(i,j) .EQV.                   &
                           (interp_index == 0) ) )
          Enddo
       Enddo
!
    Else   ! Spatial degradation is switched off
!
       Do j = 1, rows
          Do i = 1, row_length
             switch(i,j) = day_fraction(i,j).GT.0.
          Enddo
       Enddo
!
    Endif   ! If l_rad_deg
!
IF (.NOT. l_vatpoles) THEN
    If (model_domain .eq. mt_global) Then
!
! If at a pole then only calculate first point of polar row if all
! are lit. ie set switch to false for all other points.
!
       If (at_extremity(PNorth) .AND.                             &
          ( day_fraction(1,rows) > 0.0 ) ) Then
          switch(1,rows) = .TRUE.
          Do i = 2, row_length
             switch(i,rows) = .false.
          End Do
       End If
       If (at_extremity(PSouth) .AND.                             &
          ( day_fraction(1,1) > 0.0 ) ) Then
          switch(1,1) = .TRUE.
          Do i = 2, row_length
             switch(i,1) = .false.
          End Do
       End If
    End If
END IF ! vatpoles

    daylight_points = 0
    Do j = 1, rows
       Do i = 1, row_length
          If (L_subsample) Then
! DEPENDS ON: in_footprint
             l_viewed = in_footprint(j_sw, .TRUE.,                &
               L_geostationary,                                   &
               min_view_lon, max_view_lon,                        &
               min_view_lat, max_view_lat,                        &
               true_latitude(i,j),                                &
               true_longitude(i,j),                               &
               seconds_since_midnight)
          Else
                l_viewed = .TRUE.
          Endif
          If (switch(i,j).and.l_viewed) Then
             daylight_points = daylight_points + 1
             list_daylight_points(daylight_points) =              &
                                         (j-1)*row_length + i
!
! The following arrays are 2-D analogues of the above
! used for diagnostic purposes from 5.3.
!
             diag_row_list(daylight_points) = j
             diag_col_list(daylight_points) = i
             diag_row_list_sw(daylight_points) = j
             diag_col_list_sw(daylight_points) = i
          End If
       End Do
    End Do
!
! Add up the total number of daylight points in the whole domain using
! the GC routine ISUM (integer sum).
!
    IF (L_rad_deg) THEN
       tot_daylight_points=daylight_points
       CALL GC_ISUM(1,n_proc,Error_code,tot_daylight_points)
       If ( error_code /= 0 ) Then

           Call Ereport('NI_rad_ctl', Error_code,                 &
                        'Unable to determine total lit points.')
       Endif
    ELSE
       tot_daylight_points=0
    ENDIF   ! If l_rad_deg
    IF (lhook) CALL dr_hook('PRELIM_SWRAD',zhook_out,zhook_handle)
    RETURN
!
  end subroutine prelim_swrad
