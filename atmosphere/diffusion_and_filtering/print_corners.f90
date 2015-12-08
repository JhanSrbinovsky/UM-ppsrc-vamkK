! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************

! subroutine print_corners

      subroutine print_corners(                                         &
     &                         field, row_length, rows, levels,         &
     &                         level_1, level_2,                        &
     &                         off_x, off_y, off_u, off_v,              &
     &                         l_datastart, n_proc, model_domain )

     USE PR_CORNER_MOD
     USE parkind1, ONLY: jprb, jpim
     USE yomhook, ONLY: lhook, dr_hook

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!     
! Purpose:
!          To print (sub-)domain corner values of a field 
!
! Method:
!          Corner values in PR_CORNER_MOD and are set from
!          the run_diagnostics NAMELIST and are initialised in 
!          init_pr_corner called from  SETCONA
!          Place CALL anywhere in the UM to print arrays
!          e.g. printing theta values from ATM_STEP
!          
!  Whole domain values 
!      dom_w = 1         dom_e = global_row_length
!      dom_s = 1         dom_n = global_rows
      
!   Corner size ix,jy  do not have to be equal
!    Max values  ix = 4 ,    jy = 4         
!     write(6,*)' *** THETA  ***'
!! DEPENDS ON: print_corners
!      call print_corners(                                              &
!     &                    THETA, row_length, rows, model_levels,        &
!     &                    level_1, level_2,
!     &                    offx, offy, 0, 0,                             &
!     &                    datastart, nproc, model_domain )
!
!         off_u  should be set to 1 for calling U fields in a LAM 
!         off_v  should be set to 1 for calling V fields 
!
!  NB Calling routine must contain the argument list values
!     In particular datastart(2), nproc, and model_domain 
!     In some subroutines
!         datastart  may be named l_datastart
!         nproc may be named n_proc 
!         offx, offy may be named off_x, off_y
!
!          The corners are defined by
!                   left(West)    right(East)
!    upper row  =  dom_w,dom_n    dom_e,dom_n
!    lower row  =  dom_w,dom_s    dom_e,dom_s
!
!           Values are printed out as below
!
!                      1                             ix
! upper      left   dom_w,dom_n            ... dom_w+ix-1,dom_n
! upper-jy+1 left   dom_w,dom_n-jy+1       ... dom_w+ix-1,dom_n-jy+1
!
!                      1                             ix
! upper      right  dom_e-ix+1,dom_n       ... dom_e-ix+1,dom_n
! upper-jy+1 right  dom_e-ix+1,dom_n-jy+1  ... dom_e-ix+1,dom_n-jy+1
!
!                      1                             ix
! lower+jy-1 left   dom_w,dom_s+jy-1       ... dom_w+ix-1,dom_s+jy-1
! lower      left   dom_w,dom_s            ... dom_w+ix-1,dom_s
!
!                      1                             ix
! lower+jy-1 right  dom_e-ix+1,dom_s+jy-1  ... dom_e-ix+1,dom_s+jy-1
! lower      right  dom_e-ix+1,dom_s       ... dom_e-ix+1,dom_s

      Implicit None

! Parameters required for dimensioning some of the arguments

! Arguments with Intent IN. ie: Input variables.
      INTEGER, Intent(In) ::                                            &
     &  row_length                                                      &
                           ! in: no of points per local row
     &, rows                                                            &
                           ! in: no of local rows
     &, levels                                                          &
                           ! in: no of levels in input field
     &, level_1                                                         &
                           ! in: level to print
     &, level_2                                                         &
                           ! in: level to print
     &, off_x , off_y                                                   &
                           !  halo sizes of field
     &, off_u                                                           &
                           !  offsets(=1) for LAM u type fields
     &, off_v                                                           &
                           !  offsets(=1) for LAM v type fields
     &, l_datastart(3)                                                  &
     &, n_proc                                                          &
                           ! Total number of processors
     &, model_domain

      Real, Intent(In) ::                                               &
     & field(1-off_x:row_length+off_x,1-off_y:rows+off_y, levels)

!    local variables

! Loop indices/pointers
      Integer                                                           &
     &  i, j, k                                                         &
     &, ii, ji                                                          &
     &, info                                                            &
     &, dom_n, dom_e                                                    &
     &, ie, ies, iee                                                    &
     &, iw, iws, iwe                                                    &
     &, jn, jns, jne                                                    &
     &, js, jss, jse                                                    &
     &, lev_in(2)                                                       &
     &, num_lev                                                         &
     &, off_i                                                           &
                            ! off_set for u if LAM
     &, size

!    local arrays

      Real :: ll(4,4,2)
      Real :: lr(4,4,2)
      Real :: ul(4,4,2)
      Real :: ur(4,4,2)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 1.  Set up pointers
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('PRINT_CORNERS',zhook_in,zhook_handle)

! Subtract offsets if set
        dom_n = dom_no - off_v
        dom_e = dom_eo - off_u

!   Level numbers of field levels to be printed choose 2 values
!   If levels = 1 (i.e. single level field) only 1 level printed

      lev_in(1) = level_1
      lev_in(2) = level_2
      
      num_lev = 2
      if ( lev_in(1) == lev_in(2)) num_lev = 1

      if (model_domain == mt_LAM) then
        off_i = off_u
      else !  model_domain \= mt_LAM
        off_i = 0
      end if !  model_domain == mt_LAM
      
! Set defaults to no-op values
!  (i.e. arrays left empty unless at a corner of desired sub-domain)    
      iws = row_length + 1
      iwe = 0
      ies = row_length + 1
      iee = 0
      jss = rows + 1
      jse = 0
      jns = rows + 1
      jne = 0
      iw = 0
      ie = 0
      jn = 0
      js = 0
            
      if (dom_w > l_datastart(1) - 1 .and.                              &
     &    dom_w < l_datastart(1) + row_length - off_i ) Then
        iws = dom_w - l_datastart(1) + 1
        iwe = iws + ix - 1
        if ( iwe > row_length ) iwe = row_length
      else if (dom_w + ix > l_datastart(1) .and.                        &
     &              dom_w < l_datastart(1)) Then
        iws = 1
        iwe = dom_w + ix - l_datastart(1)
        iw = ix - iwe + iws
      end if !  dom_w > l_datastart(1) + row_length - 1

      if (dom_e - ix + 1 > l_datastart(1) - 1 .and.                     &
     &    dom_e - ix + 1 < l_datastart(1) + row_length - off_i ) Then
        ies = dom_e - ix - l_datastart(1) + 2
        iee = ies + ix - 1
        if ( iee > row_length ) iee = row_length
      else if (dom_e + 1 > l_datastart(1) .and.                         &
     &    dom_e - ix + 1 < l_datastart(1)) Then
        ies = 1
        iee = dom_e + 1 - l_datastart(1)
        ie = ix - iee + ies
      end if !  dom_e - ix + 1 > l_datastart(1) + row_length - 1

      if (dom_s > l_datastart(2) - 1 .and.                              &
     &    dom_s < l_datastart(2) + rows ) Then
        jss = dom_s - l_datastart(2) + 1
        jse = jss + jy - 1
        if ( jse > rows ) jse = rows
      else if (dom_s + jy > l_datastart(2) .and.                        &
     &              dom_s < l_datastart(2)) Then
        jss = 1
        jse = dom_s + jy - l_datastart(2)
        js = jy - jse + jss
      end if !  dom_s > l_datastart(2) + rows - 1

      if (dom_n - jy + 1 > l_datastart(2) - 1 .and.                     &
     &    dom_n - jy + 1 < l_datastart(2) + rows ) Then
        jns = dom_n - jy - l_datastart(2) + 2
        jne = jns + jy - 1
        if ( jne > rows ) jne = rows
      else if (dom_n + 1 > l_datastart(2) .and.                         &
     &    dom_n - jy + 1 < l_datastart(2)) Then
        jns = 1
        jne = dom_n + 1 - l_datastart(2)
        jn = jy - jne + jns
      end if !  dom_n - jy + 1 > l_datastart(2) + rows - 1

! ----------------------------------------------------------------------
! Section 2.  Initialise arrays
! ----------------------------------------------------------------------

      ll = 0.0      ! lower left values
      lr = 0.0      ! lower right values
      ul = 0.0      ! upper left values
      ur = 0.0      ! upper right values

! ----------------------------------------------------------------------
! Section 3.  Fill arrays
! ----------------------------------------------------------------------
! Each array ll, lr, ul, ur is filled by only 1 pe
! Columns are filled in reverse order (North to South)

      do k = 1, num_lev

!   ll is lower left corner
        do j = jss, jse
          ji = jse - j + 1 + js
          do i = iws, iwe
            ii = i - iws + 1 + iw
            ll(ii,ji,k) = field(i,j,lev_in(k))
          end do  !  i = iws, iwe
        end do ! j = jss, jse

!   lr is lower right corner
        do j = jss, jse
          ji = jse - j + 1 + js
          do i = ies, iee
            ii = i - ies + 1 + ie
            lr(ii,ji,k) = field(i,j,lev_in(k))
          end do  !  i = ies, iee
        end do ! j = jss, jse

!   ul is upper left corner
        do j = jns, jne
          ji = jne - j + 1 + jn
          do i = iws, iwe
            ii = i - iws + 1 + iw
            ul(ii,ji,k) = field(i,j,lev_in(k))
          end do  !  i = iws, iwe
        end do ! j = jns, jne

!   ur is upper right corner
        do j = jns, jne
          ji = jne - j + 1 + jn
          do i = ies, iee
            ii = i - ies + 1 + ie
            ur(ii,ji,k) = field(i,j,lev_in(k))
          end do  !  i = ies, iee
        end do ! j = jns, jne

      end do  !  k = 1, num_lev

! Each element of ll, lr, ul, ur is filled by only 1 pe

! Sum over pes to put same value on all processors (hence pe0)
! Length = size = max_x*max_y*num_lev arrays

      size = 4 * 4 * num_lev

      call gc_rsum(size, n_proc, info, ll)
      call gc_rsum(size, n_proc, info, lr)
      call gc_rsum(size, n_proc, info, ul)
      call gc_rsum(size, n_proc, info, ur)

! ----------------------------------------------------------------------
! Section 2.  Print corners
! ----------------------------------------------------------------------

      do k = 1, num_lev
 
        write(6,'(''level '',I3,'' upper left'')') lev_in(k)
        ii = dom_w
        if (ix == 4) then
          write(6,'(''points             '',I4'//                          &
              ',''                     '',I4,''                   '',I4'// &
              ',''                    '',I4)') ii, ii+1, ii+2, ii+3
        else if (ix == 3) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4,''                   '',I4)')&
              ii, ii+1, ii+2
        else if (ix == 2) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4)') ii, ii+1
        else 
          write(6,'(''points             '',I4)')ii
        endif !  ix == 4          
        Do j = 1, jy
          jn = dom_n - j + 1 
          write(6,'(''row '',I4,4E24.16)') jn, (ul(i,j,k), i = 1, ix)
        End Do
        write(6,'(''level '',I3,'' upper right'')') lev_in(k)
        ii = dom_e - ix + 1
        if (ix == 4) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4,''                   '',I4'//&
              ',''                    '',I4)') ii, ii+1, ii+2, ii+3
        else if (ix == 3) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4,''                   '',I4)')&
              ii, ii+1, ii+2
        else if (ix == 2) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4)') ii, ii+1
        else 
          write(6,'(''points             '',I4)')ii
        endif !  ix == 4          
        Do j = 1, jy
          jn = dom_n - j + 1
          write(6,'(''row '',I4,4E24.16)') jn, (ur(i,j,k), i = 1, ix)
        End Do
        write(6,'(''level '',I3,'' lower left'')') lev_in(k)
        ii = dom_w
        if (ix == 4) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4,''                   '',I4'//&
              ',''                    '',I4)') ii, ii+1, ii+2, ii+3
        else if (ix == 3) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4,''                   '',I4)')&
              ii, ii+1, ii+2
        else if (ix == 2) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4)') ii, ii+1
        else 
          write(6,'(''points             '',I4)')ii
        endif !  ix == 4          
        Do j = 1, jy
          jn = dom_s + jy - j 
          write(6,'(''row '',I4,4E24.16)') jn, (ll(i,j,k), i = 1, ix)
        End Do
        write(6,'(''level '',I3,'' lower right'')') lev_in(k)
        ii = dom_e - ix + 1
        if (ix == 4) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4,''                   '',I4'//&
              ',''                    '',I4)') ii, ii+1, ii+2, ii+3
        else if (ix == 3) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4,''                   '',I4)')&
              ii, ii+1, ii+2
        else if (ix == 2) then
          write(6,'(''points             '',I4'//&
              ',''                     '',I4)') ii, ii+1
        else 
          write(6,'(''points             '',I4)')ii
        endif !  ix == 4          
        Do j = 1, jy
          jn = dom_s + jy - j 
          write(6,'(''row '',I4,4E24.16)') jn, (lr(i,j,k), i = 1, ix)
        End Do

      End Do !  k = 1, num_lev

      IF (lhook) CALL dr_hook('PRINT_CORNERS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE print_corners

