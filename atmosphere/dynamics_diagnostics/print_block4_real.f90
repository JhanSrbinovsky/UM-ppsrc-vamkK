! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************

! subroutine print_block4_real

      SUBROUTINE print_block4_real(                                     &
                                   field, row_length, rows, levels,     &
                                   level_1, level_2,                    &
                                   off_x, off_y, off_u, off_v,          &
                                   print_time )

      USE proc_info_mod
      USE pr_block4_mod
      USE timestep_mod

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

! Purpose:
!          To print (sub-)domain corner values of a field

! Method:
!          Corner values in PR_BLOCK4_MOD and are set from
!          the run_diagnostics NAMELIST and are initialised in
!          init_pr_corner CALLed from  SETCONA
!          Place CALL anywhere in the UM to print arrays
!          e.g. printing theta values from ATM_STEP
!          Use tkdiff to compare output from different runs
!
!  Whole domain values
!      dom_w = 1         dom_e = global_row_length
!      dom_s = 1         dom_n = global_rows

!   Corner size ix,jy  do not have to be equal
!    Max values  ix = 4 ,    jy = 4
!     WRITE(6,*)' *** THETA  ***'
!! DEPENDS ON: print_block4_real
!      CALL print_block4_real(
!     &                       THETA,
!     &                       row_length, rows, model_levels,
!     &                       level_1, level_2,
!     &                       offx, offy, off_u, off_v )
!
!         off_u/off_v  should be set to 1 if the right/upper
!         dom_e-1/dom_n-1 wind values values are to be printed.
!         rather than the dom_e/dom_n values
!         e.g for a LAM near the right/upper boundary or for
!             global model v-fields at the North Pole or
!             for comparing values that might have made lbcs
!
!     In some subroutines
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

      IMPLICIT NONE

! Parameters required for dimensioning some of the arguments

! Arguments with Intent IN. ie: Input variables.
      INTEGER, INTENT(IN) ::                                            &
        row_length,                                                     &
                           ! in: no of points per local row
        rows,                                                           &
                           ! in: no of local rows
        levels,                                                         &
                           ! in: no of levels in input field
        level_1,                                                        &
                           ! in: level to print
        level_2,                                                        &
                           ! in: level to print
        off_x , off_y,                                                  &
                           !  halo sizes of field
        off_u,                                                          &
                           !  offsets(=1) for LAM u type fields
        off_v ,                                                         &
                           !  offsets(=1) for LAM v type fields
        print_time         ! timestep_number for print output

      REAL, INTENT(IN) ::                                               &
        field (1-off_x:row_length+off_x,1-off_y:rows+off_y, levels)

!    local variables

! Loop indices/pointers
      INTEGER                                                           &
        i, j, k,                                                        &
        ii, ji,                                                         &
        info,                                                           &
        dom_n, dom_e,                                                   &
        ie, ies, iee,                                                   &
        iw, iws, iwe,                                                   &
        jn, jns, jne,                                                   &
        js, jss, jse,                                                   &
        lev_in(2),                                                      &
        num_lev,                                                        &
        size

!    local arrays

      REAL :: ll(4,4,2)
      REAL :: lr(4,4,2)
      REAL :: ul(4,4,2)
      REAL :: ur(4,4,2)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle   

! ----------------------------------------------------------------------
! Section 1.  Set up pointers
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('PRINT_BLOCK4_REAL',zhook_in,zhook_handle)

      IF ( print_time == timestep_number ) THEN

! Subtract offsets if set
        dom_n = dom_no - off_v
        dom_e = dom_eo - off_u

!   Level numbers of field levels to be printed choose 2 values
!   If levels = 1 (i.e. single level field) only 1 level printed

        lev_in(1) = level_1
        lev_in(2) = level_2

        num_lev = 2
        IF ( lev_in(1) == lev_in(2)) num_lev = 1

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
        jn = jy + 1
        js = jy + 1

        IF (dom_w > datastart(1) - 1 .AND.                              &
            dom_w < datastart(1) + row_length ) THEN
          iws = dom_w - datastart(1) + 1
          iwe = iws + ix - 1
          IF ( iwe > row_length ) iwe = row_length
        ELSE IF (dom_w + ix > datastart(1) .AND.                      &
                      dom_w < datastart(1)) THEN
          iws = 1
          iwe = dom_w + ix - datastart(1)
          iw = ix - iwe + iws - 1
        END IF !  dom_w > datastart(1) + row_length - 1

        IF (dom_e - ix + 1 > datastart(1) - 1 .AND.                   &
            dom_e - ix + 1 < datastart(1) + row_length ) THEN
          ies = dom_e - ix - datastart(1) + 2
          iee = ies + ix - 1
          IF ( iee > row_length ) iee = row_length
        ELSE IF (dom_e + 1 > datastart(1) .AND.                       &
                 dom_e - ix + 1 < datastart(1)) THEN
          ies = 1
          iee = dom_e + 1 - datastart(1)
          ie = ix - iee + ies - 1
        END IF !  dom_e - ix + 1 > datastart(1) - 1

        IF (dom_s > datastart(2) - 1 .AND.                            &
            dom_s < datastart(2) + rows ) THEN
          jss = dom_s - datastart(2) + 1
          jse = jss + jy - 1
          IF ( jse > rows ) jse = rows
        ELSE IF (dom_s + jy > datastart(2) .AND.                      &
                      dom_s < datastart(2)) THEN
          jss = 1
          jse = dom_s + jy - datastart(2)
          js = jse - jss + 2
        END IF !  dom_s > datastart(2) + rows - 1

        IF (dom_n - jy + 1 > datastart(2) - 1 .AND.                   &
            dom_n - jy + 1 < datastart(2) + rows ) THEN
          jns = dom_n - jy - datastart(2) + 2
          jne = jns + jy - 1
          IF ( jne > rows ) jne = rows
        ELSE IF (dom_n + 1 > datastart(2) .AND.                       &
                 dom_n - jy + 1 < datastart(2)) THEN
          jns = 1
          jne = dom_n + 1 - datastart(2)
          jn = jne - jns + 2
        END IF !  dom_n - jy + 1 > datastart(2) + rows - 1

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

        DO k = 1, num_lev

!   ll is lower left corner
          ji = js
          DO j = jss, jse
            ji = ji - 1
            ii = iw
            DO i = iws, iwe
              ii = ii + 1
              ll(ii,ji,k) = field(i,j,lev_in(k))
            END DO  !  i = iws, iwe
          END DO ! j = jss, jse

!   lr is lower right corner
          ji = js
          DO j = jss, jse
            ji = ji - 1
            ii = ie
            DO i = ies, iee
              ii = ii + 1
              lr(ii,ji,k) = field(i,j,lev_in(k))
            END DO  !  i = ies, iee
          END DO ! j = jss, jse

!   ul is upper left corner
          ji = jn
          DO j = jns, jne
            ji = ji - 1
            ii = iw
            DO i = iws, iwe
              ii = ii + 1
              ul(ii,ji,k) = field(i,j,lev_in(k))
            END DO  !  i = iws, iwe
          END DO ! j = jns, jne

!   ur is upper right corner
          ji = jn
          DO j = jns, jne
            ji = ji - 1
            ii = ie
            DO i = ies, iee
              ii = ii + 1
              ur(ii,ji,k) = field(i,j,lev_in(k))
            END DO  !  i = ies, iee
          END DO ! j = jns, jne

        END DO  !  k = 1, num_lev

! Each element of ll, lr, ul, ur is filled by only 1 pe

! Sum over pes to put same value on all processors (hence pe0)
! Length = size = max_x*max_y*num_lev arrays

        size = 4 * 4 * num_lev

        CALL gc_rsum(size, n_proc, info, ll)
        CALL gc_rsum(size, n_proc, info, lr)
        CALL gc_rsum(size, n_proc, info, ul)
        CALL gc_rsum(size, n_proc, info, ur)

! ----------------------------------------------------------------------
! Section 2.  Print corners
! ----------------------------------------------------------------------

        DO k = 1, num_lev

          WRITE(6,'(''level '',I3,'' upper left'')') lev_in(k)
          ii = dom_w
          IF (ix == 4) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4,''                   '', I4,  &
            ''                    '', I4)') ii, ii+1, ii+2, ii+3
          ELSE IF (ix == 3) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4,''                   '', I4)')&
                                                       ii, ii+1, ii+2
          ELSE IF (ix == 2) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4)') ii, ii+1
          ELSE
            WRITE(6,'(''column             '', I4)')ii
          END IF !  ix == 4
          DO j = 1, jy
            jn = dom_n - j + 1
            WRITE(6,'(''row '',I4,4E24.16)') jn, (ul(i,j,k), i = 1, ix)
          END DO
          WRITE(6,'(''level '',I3,'' upper right'')') lev_in(k)
          ii = dom_e - ix + 1
          IF (ix == 4) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4,''                   '', I4,  &
            ''                    '', I4)') ii, ii+1, ii+2, ii+3
          ELSE IF (ix == 3) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4,''                   '', I4)')&
                                          ii, ii+1, ii+2
          ELSE IF (ix == 2) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4)') ii, ii+1
          ELSE
            WRITE(6,'(''column             '', I4)') ii
          END IF !  ix == 4
          DO j = 1, jy
            jn = dom_n - j + 1
            WRITE(6,'(''row '',I4,4E24.16)') jn, (ur(i,j,k), i = 1, ix)
          END DO
          WRITE(6,'(''level '',I3,'' lower left'')') lev_in(k)
          ii = dom_w
          IF (ix == 4) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4,''                   '', I4,  &
            ''                    '', I4)') ii, ii+1, ii+2, ii+3
          ELSE IF (ix == 3) THEN
          WRITE(6,'(''column             '', I4,                        &
          ''                     '', I4,''                   '', I4)')  &
                                         ii, ii+1, ii+2
          ELSE IF (ix == 2) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4)') ii, ii+1
          ELSE
            WRITE(6,'(''column             '', I4)') ii
          END IF !  ix == 4
          DO j = 1, jy
            jn = dom_s + jy - j
            WRITE(6,'(''row '',I4,4E24.16)') jn, (ll(i,j,k), i = 1, ix)
          END DO
          WRITE(6,'(''level '',I3,'' lower right'')') lev_in(k)
          ii = dom_e - ix + 1
          IF (ix == 4) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4,''                   '', I4,  &
            ''                    '', I4)') ii, ii+1, ii+2, ii+3
          ELSE IF (ix == 3) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4,''                   '', I4)')&
                                          ii, ii+1, ii+2
          ELSE IF (ix == 2) THEN
            WRITE(6,'(''column             '', I4,                      &
            ''                     '', I4)') ii, ii+1
          ELSE
            WRITE(6,'(''column             '', I4)')ii
          END IF !  ix == 4
          DO j = 1, jy
            jn = dom_s + jy - j
            WRITE(6,'(''row '',I4,4E24.16)') jn, (lr(i,j,k), i = 1, ix)
          END DO

        END DO !  k = 1, num_lev

      END IF ! print_time == timestep_number

      IF (lhook) CALL dr_hook('PRINT_BLOCK4_REAL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE print_block4_real
