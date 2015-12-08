! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine Calc_var_grid -----------------------------------------
!
! Purpose : to derive variable grid parameters at each point
!           lambda's and phi's no-longer have fixed intervals
!
! Method:
!    1. See paper:
!                Setting and testing variable resolution grids
!                inside LAM domains by Terry Davies
!    2. Needs to be called separately for each direction
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level

      SUBROUTINE Calc_var_grid (                                        &
     &                        lamphi_p, lamphi_uv,                      &
     &                        model_domain, L_lambda, L_dprint,         &
     &                        global_pts, datastart, dir_pt,            &
     &                        proc_p_pts, proc_uv_pts, var_pts,         &
     &                        delta_in, origin,                         &
     &                        grid_ratio, var_ratio, position )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE


      INTEGER                                                           &
     &  global_pts                                                      &
                      !IN total number of point in a row or column
     &, proc_p_pts                                                      &
                      !IN  p points in a row or column on this pe
     &, proc_uv_pts                                                     &
                      !IN points in a u-row or v-column on this pe
     &, model_domain                                                    &
     &, datastart(2)                                                    &
                      ! First gridpoints held by this processor
     &, dir_pt                                                          &
                      !IN  datastart(dir_pt) to point to first point
     &, var_pts  
                      ! Size of each variable grid zone

! Input arguments
      REAL                                                              &
     &  delta_in                                                        &
                  ! Original grid spacing in degrees
     &, origin                                                          &
                  ! Lat or lon of first theta point in degrees
     &, grid_ratio                                                      &
                  ! Ratio original LAM gridlength:high res gridlength
     &, var_ratio                                                       &
                  ! Stretching factor for variable grid
     &, position  
                  ! Position of high res zone .5 is middle
!         make position  smaller to move left/down, bigger for right/up

      REAL:: lamphi_p  ( proc_p_pts )
      REAL:: lamphi_uv ( proc_uv_pts )

      Logical                                                           &
     &  L_lambda                                                        &
                  ! True if longitude, false if latitude
     &, L_dprint  ! True if diagnostic prints needed

! local  variables
      CHARACTER(LEN=80)                                                    &
     &       CMESSAGE              ! Error message if ICODE >0
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='Calc_var_grid')

      Real                                                              &
     &  wk_lamphi_p  ( global_pts)                                      &
     &, wk_lamphi_uv ( global_pts)
          
      Integer                                                           &
     &  high_pts                                                        &
     &, reg_pts                                                         &
     &, right_pts                                                       &
     &, left_pts

      REAL                                                              &
     &  delta_high                                                      &
     &, delta_reg                                                       &
     &, delta_var                                                       &
     &, real_pts                                                        &
     &, sum_var                                                         &
     &, denom                                                           &
     &, var_power                                                       &
     &, recip_ratio                                                     &
     &, recip_varm1

! loop bounds/counters
      Integer                                                           &
     &  i, gi                                                           &
     &, start_loop                                                      &
     &, end_loop                                                        &
     &, icode              ! Error message if ICODE >0

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('CALC_VAR_GRID',zhook_in,zhook_handle)
      if (L_dprint) then
        write(*,*) 'global_pts = ' , global_pts
        if ( L_lambda ) then  !longitudes
          if ( origin == 0.0 ) then
            write(*,*) 'First longitude point of grid is at ' , origin, &
     &                 ' degrees  (Greenwich Meridian)'
          elseif ( origin < 180.0 ) then
            write(*,*) 'First longitude point of grid is at ' , origin, &
     &                 ' degrees East'
          else ! origin > 180.0 ) then
            write(*,*) 'First longitude point of grid is at '           &
     &                 , origin - 180.0, ' degrees West'
            write(*,*) 'Value of origin = ', origin
          endif ! origin == 0.0
        else ! latitudes
          if ( origin == 0.0 ) then
            write(*,*) 'First latitude point of grid is at ' , origin,  &
     &                 ' degrees  (Equator)'
          elseif ( origin < 0.0 ) then
            write(*,*) 'First latitude point of grid is at ' , origin,  &
     &                 ' degrees South'
          else ! origin > 0.0 ) then
            write(*,*) 'First latitude point of grid is at ' , origin,  &
     &                 ' degrees North'
          endif ! origin == 0.0
        endif ! L_lambda
        write(*,*) 'delta_in = ' , delta_in
      endif !  L_dprint

! set wk_lamphi_p for all points in LAM domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Section 1:  Check input parameters are sensible
!             and set up constants etc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      grid_ratio = alpha,    var_ratio = r
      if( grid_ratio > 1.0 ) then
        recip_ratio = 1.0 / var_ratio
        recip_varm1 = 1.0 / (var_ratio - 1.0)
        if( var_pts < grid_ratio + 1 ) then
! need to reset width of variable zone as will not have enough points
!    to fill variable domain
          if( var_ratio > 1.0 ) then
            denom = log( 3.0 ) / log( var_ratio )
            var_pts = NINT(denom)
          else !  var_ratio <= 1.0
            icode = 10
            write(*,*) 'var_ratio = ',var_ratio
            Cmessage = 'var_ratio must be set > 1'

            Call Ereport( RoutineName, icode, Cmessage )
          endif !  var_ratio > 1.0
        endif !  var_pts < grid_ratio + 1
        var_power = var_ratio ** var_pts
        denom = 1.0 / ( var_power - 1.0 )
        sum_var = recip_varm1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Section 2: Calculate intervals for new grid
!
!   Old grid    <----------------      N-1              --------------->
!   Old grid    <------------       delta_in            --------------->
! New <-    X1   -><-    Y    -><-    Z     -><-    Y    -><-    X2   ->
! New <-delta_reg-><-delta_var-><-delta_high-><-delta_var-><-delta_reg->
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      real_pts = X1 + X2 = left_pts + right_pts,  global_pts = N
        real_pts = ( (global_pts - 1.0) * (grid_ratio - 1.0) +          &
     &                  var_pts + var_pts ) * denom - 2.0 * sum_var
        reg_pts = NINT(real_pts)
        left_pts = INT(real_pts * position)
        right_pts = reg_pts - left_pts
        sum_var = sum_var * (var_power - 1.0)
!      high_pts = Z,  var_pts = Y
        high_pts = global_pts - left_pts - right_pts -                  &
     &                                       var_pts - var_pts - 1
        delta_high =  delta_in / grid_ratio
        delta_reg = delta_high * var_power
! rescale and recalculate delta_lambdas to make domain lengths match
        delta_high = (global_pts - 1.0) * delta_in /                    &
     &                    ( (left_pts + right_pts) * var_power +        &
     &                                    2.0 * sum_var + high_pts )
        delta_reg = delta_high * var_power
        if (L_dprint) then
          write(*,*) 'high_pts = ', high_pts,' left_pts = ', left_pts,  &
     &             ' right_pts = ', right_pts, ' var_pts = ', var_pts
          write(*,*) 'delta_reg = ' , delta_reg
          write(*,*) 'sum_var = ' , sum_var
          write(*,*) 'delta_high =  ' , delta_high
          write(*,*) 'delta_in / delta_high =  ',                       &
     &                delta_in / delta_high
        endif !  L_dprint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Section 3: Place lambda/phi values for complete domain in wk arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        wk_lamphi_p(1) = origin
        start_loop = 2
        end_loop = left_pts + 1
        do i =  start_loop, end_loop
          wk_lamphi_p(i) = wk_lamphi_p(i-1) + delta_reg
        enddo ! i =  start_loop, end_loop
        start_loop = end_loop + 1
        end_loop = end_loop + var_pts
        delta_var = delta_reg * recip_ratio
        do i =  start_loop, end_loop
          wk_lamphi_p(i) = wk_lamphi_p(i-1) + delta_var
          delta_var = delta_var * recip_ratio
        enddo ! i =  start_loop, end_loop
        start_loop = end_loop + 1
        end_loop = end_loop + high_pts
        do i =  start_loop, end_loop
          wk_lamphi_p(i) = wk_lamphi_p(i-1) + delta_high
        enddo ! i =  start_loop, end_loop
        start_loop = end_loop + 1
        end_loop = end_loop + var_pts
        delta_var = delta_high
        do i =  start_loop, end_loop
          wk_lamphi_p(i) = wk_lamphi_p(i-1) + delta_var
          delta_var = delta_var * var_ratio
        enddo ! i =  start_loop, end_loop
        start_loop = end_loop + 1
        end_loop = end_loop + right_pts
        do i =  start_loop, end_loop
          wk_lamphi_p(i) = wk_lamphi_p(i-1) + delta_reg
        enddo ! i =  start_loop, end_loop
        if (L_dprint) then
          write(*,*) 'Number of lamphi_p points calculated = ', end_loop
        endif  !  L_dprint

      else !  grid_ratio = 1.0
!  regular mesh settings

        delta_reg = delta_in
        wk_lamphi_p(1) = origin
        start_loop = 2
        end_loop = global_pts
        do i =  start_loop, end_loop
          wk_lamphi_p(i) = wk_lamphi_p(i-1) + delta_reg
        enddo ! i =  start_loop, end_loop
        if (L_dprint) then
          write(*,*) '  *****   Regular mesh  ***** '
          write(*,*) ' p points set = ', end_loop
        endif !  L_dprint
!  end of regular mesh settings

      endif ! grid_ratio > 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Section 4: Copy lambda/phi values to local processor
!               and convert to radians
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      start_loop = 1
      end_loop = global_pts - 1
      do i =  start_loop, end_loop
        wk_lamphi_uv(i) = 0.5 * (wk_lamphi_p(i) + wk_lamphi_p(i+1))
      enddo ! i =  start_loop, end_loop
      wk_lamphi_uv(end_loop + 1) = wk_lamphi_uv(end_loop) +  delta_reg

      if (L_dprint) then
        write(*,*) 'wk_lamphi_p =', wk_lamphi_p
        write(*,*) 'wk_lamphi_uv =', wk_lamphi_uv
      endif !  L_dprint

! Do not convert grid-spacing and base values from degrees to radians
! since this will be done later in SETCON
!  Just copy into output array as required
      If ( proc_p_pts == global_pts ) Then

        do i = 1 , proc_p_pts
          lamphi_p(i) =  wk_lamphi_p(i)
        end do

        do i = 1 , proc_uv_pts
          lamphi_uv(i) =  wk_lamphi_uv(i)
        end do

      Else !   proc_p_pts < global_pts

        do i = 1 , proc_p_pts
          gi = datastart(dir_pt) + i - 1
          lamphi_p(i) =  wk_lamphi_p(gi)
        end do

        do i = 1 , proc_uv_pts
          gi = datastart(dir_pt) + i - 1
          lamphi_uv(i) =  wk_lamphi_uv(gi)
        end do

      EndIf !   proc_p_pts == global_pts

      IF (lhook) CALL dr_hook('CALC_VAR_GRID',zhook_out,zhook_handle)
      RETURN
      End SUBROUTINE Calc_var_grid
