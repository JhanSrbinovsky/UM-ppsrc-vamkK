! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine Set_var_look --------------------------------------------
!
! Purpose : to set up look-up table to do variable res searches
!
! Method:
!      1. Set up set of points/intervals at smallest grid length
!      2. Map indices onto variable grid indices and place in
!         look-up arrays
!      3. Needs to be called separately for each direction and
!         variable
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level

      SUBROUTINE Set_var_look(                                          &
     &                        lamphi, rowcol, halo, recip_delta,        &
     &                        look, max_pts, halo_look, L_dprint )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE

! Input arguments

      Logical                                                           &
     &  L_dprint  ! True if diagnostic prints needed
      
      INTEGER                                                           &
     &  rowcol                                                          &
                !IN total number of point in a row or column
     &, halo                                                            &
                !IN halo size of lamphi array
     &, max_pts 
                !IN  max points row or column at highest resolution

      REAL:: lamphi ( 1 - halo : rowcol + halo )

! Output arguments
      REAL                                                              &
     &  recip_delta  ! grid-length for look-up table

      Integer                                                           &
     &  halo_look    ! halo for look-up table

      Integer                                                           &
     &  look( max_pts )

! local  variables
      CHARACTER(LEN=80)                                                    &
     &       CMESSAGE              ! Error message if ICODE >0
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='Set_var_look')

      Real                                                              &
     &  lamphi_fine ( max_pts )

     Integer                                                            &
     &  i, j                                                            &
     &, icode                                                           &
                    ! Error message if ICODE >0
     &, length                                                          &
     &, dim_fine                                                        &
     &, info                                                            &
     &, istat

      REAL                                                              &
     &  delta                                                           &
                   ! grid-length for look-up table
     &, delta_in                                                        &
     &, test                                                            &
     &, tol

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Section 1:  Set arrays and pointers for searching on variable grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (lhook) CALL dr_hook('SET_VAR_LOOK',zhook_in,zhook_handle)
      if (L_dprint) then
        write(*,*)' Set_up look-up table size  ', max_pts
        write(*,*)' Input field has', rowcol                            &
     &           ,' points plus halo size', halo
        write(*,*) 'lamphi'
        write(*,*) lamphi
      endif !  L_dprint

! Initialise values
      tol = 1.0e-8
      Do i = 1, max_pts
        look(i) = 0
      End Do !  i = 1, max_pts

      delta_in = lamphi(2) - lamphi(1)
      if (abs(delta_in - lamphi(rowcol) + lamphi(rowcol-1)) > tol) then
        write(*,*)'              WARNING'
        write(*,*)' Grid length at either end of lamphi array is'       &
     &           ,' different'
      endif ! abs(delta_in - lamphi(rowcol) + lamphi(rowcol-1)) > tol)
      length = rowcol + halo

!  Find minimum delta
      delta = 1000.0  ! big value so it isn't the minimum
      do i = 1, rowcol
        test = lamphi(i+1) - lamphi(i)
        if ( delta > test ) delta = test
      end do

! Find halo_look
!  Align fine-scale array with first point of whole domain
      test = halo * delta_in / delta
      if ( test - halo < tol ) then
        halo_look = halo
      else
        halo_look = test + 1.0
      endif ! test - halo < tol
      if (L_dprint) then
        write(*,*) 'delta = ', delta,'  halo_look = ', halo_look
      endif !  L_dprint

!  Set 1 / delta to save divide in index calculations
      recip_delta = 1.0 / delta

!  Generate fine-scale array corresponding to lamphi
      length = rowcol + halo
      lamphi_fine(1 + halo_look) = lamphi(1)
      i = 1 + halo_look
      Do ! cycle to determine number and values of small boxes
        i = i + 1
        if ( i > max_pts ) then
          icode = 10
          write(*,*)' Look-up table length is too small for this grid'
          Cmessage = 'Increase max_look value in routine CMAXSIZE'

            Call Ereport( RoutineName, icode, Cmessage )
        endif ! i > max_pts
        lamphi_fine(i) = lamphi_fine(i-1) + delta
        if ( lamphi_fine(i) >= lamphi(length) ) Exit
      End Do  ! cycle to determine number and values of small boxes
      if (L_dprint) then
        write(*,*)' Do exit i = ',i,' length = ',length
        write(*,*)' end-point of fine grid > end-point of input grid'
        write(*,*)' lamphi_fine(',i,')=', lamphi_fine(i),'>='           &
     &           ,' lamphi(',length,')=', lamphi(length)
      endif !  L_dprint
      test = lamphi(length) - lamphi_fine(i-1)
      if (test < 0.0 ) then
!       STOP. The difference between the end-point  of input grid
!       and the penultimate point of fine grid should be >= 0.0 
        Cmessage = 'End point of input and fine grids are different'

        Call Ereport( RoutineName, icode, Cmessage )
      else
        if (L_dprint) then
          write(*,*)' Difference between end-point of input grid and '
          write(*,*)' penultimate point of fine grid =', test
        endif !  L_dprint
      endif ! test < 0.0
!  At this point lam_fine spans rowcol+4*halo points of input lamphi
!                    array apart from halo_look values before point 1
      if( test <= tol .and. i-1-halo_look == length  ) then
        dim_fine = i - 1
      else
        dim_fine = i
      endif ! test <= tol.and. i-1-halo_look == length

      istat = 0  !use to set look-up for regular grid
      if( dim_fine - halo_look == length ) then

        if(abs(lamphi_fine(dim_fine) - lamphi(length)) < tol) then
          write(*,*)' grids are the same so ....'
          write(*,*)' Copy lamphi into lamphi_fine'
          Do i = 1, dim_fine
            lamphi_fine(i) = lamphi(i-halo)
          End Do !  i = 1, dim_fine
          
          istat = 1  !use to set look-up for regular grid
          
         endif ! abs(lamphi_fine(dim_fine)-lamphi(length)) < tol

      else  ! still some variable grid values to set
! Need to set halo_look values before point 1
        if (L_dprint) then
          write(*,*)' Grids not the same so set values in halo_look'
        endif !  L_dprint

        i = 1 + halo_look
        Do ! cycle to determine values of small boxes before point 1
          i = i - 1
          lamphi_fine(i) = lamphi_fine(i+1) - delta
          if ( i == 1 ) Exit
        End Do  ! cycle to determine values of small boxes

      endif !  dim_fine - halo_look == length


      Do i = dim_fine + 1, max_pts
        lamphi_fine(i) = 0.0
      End Do !  i = dim_fine + 1, max_pts

      if (L_dprint) then
        write(*,*)' Number of fine grid-points = ', dim_fine
        write(*,*)' Fine grid values '
        write(6,'(7E21.15)') (lamphi_fine(i), i = 1, dim_fine)
      endif !  L_dprint
      
!  Create look-up table
      if (istat == 1) then
!  Copy regular index corresponding to lamphi
      write(*,*)' Setting look for regular grid'
        j = 1-halo_look
        do i = 1, dim_fine
          look(i) = j
          j = j + 1
        enddo ! i = 1, dim_fine

      else  ! istat = 0  variable grid

      if (L_dprint) then
        write(*,*)' Setting look for variable grid'
      endif !  L_dprint

        i = halo_look
        j = 1
        Do     !  cycle loop over small intervals from first point
          i = i + 1
          if ( i  > dim_fine ) Exit
          if ( lamphi_fine(i) >= lamphi(j + 1) ) j = j + 1
          look(i) = j
        End Do   !  cycle loop over small intervals
        i = halo_look + 1
        j = 0
        Do     !  cycle loop over small intervals in halo_look
          i = i - 1
          if ( i  < 1 ) Exit
          if ( lamphi_fine(i) < lamphi(j) ) j = j - 1
          look(i) = j
        End Do   !  cycle loop over small intervals

      endif ! istat == 1

!  assign remaining values of look-up table to bad values
      do i = dim_fine+1, max_pts
        look(i) = -999   ! should never be selected
      enddo !  i = dim_fine+1, max_pts

      if (L_dprint) then
        write(*,*)' look '
        write(6,'(12I5)') (look(i), i = 1, dim_fine)
! DEPENDS ON: um_fort_flush
        call UM_FORT_FLUSH(6,info)
      endif !  L_dprint

      IF (lhook) CALL dr_hook('SET_VAR_LOOK',zhook_out,zhook_handle)
      RETURN
      End SUBROUTINE Set_var_look
