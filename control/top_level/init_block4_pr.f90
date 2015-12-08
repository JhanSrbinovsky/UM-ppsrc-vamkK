! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine init_block4_pr
      SUBROUTINE init_block4_pr(                                        &
                                global_row_length, global_rows,         &
                                blockx_in, blocky_in,                   &
                                dom_w_in, dom_e_in, dom_s_in, dom_n_in)
      USE pr_block4_mod
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
! Purpose:
!          To print 4 blocks of values of a field
!
! Method:
!          Initialise corner values for printing 4 blocks of values of a
!          domain. Input values from Run_Diagnostics namelist
!          Place CALL anywhere in the UM to print arrays
!          e.g. printing theta values from ATM_STEP
!
!          The corners are defined by
!                   left(West)    right(East)
!    upper row  =  dom_w, dom_n    dom_e, dom_n
!    lower row  =  dom_w, dom_s    dom_e, dom_s
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

! Arguments with Intent IN. ie: Input variables.
      INTEGER :: global_row_length
      INTEGER :: global_rows
      INTEGER :: dom_s_in
      INTEGER :: dom_n_in
      INTEGER :: dom_w_in
      INTEGER :: dom_e_in
      INTEGER :: blockx_in
      INTEGER :: blocky_in

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 1.  Set up pointers
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('INIT_BLOCK4_PR',zhook_in,zhook_handle)

!  Whole domain values - edit and re-compile for sub-domains
!      dom_w = 1         dom_eo = global_row_length
!      dom_s = 1         dom_no = global_rows

      IF ( dom_s_in == 0 ) THEN

!        Domain printing corner values need to be set from
!        RUN_Diagnostics NAMELIST
!        Defaulting to 1 for south and west and
!        global_row_length/global_rows for east/north

        dom_w = 1
        dom_eo = global_row_length
        dom_s = 1
        dom_no = global_rows

      ELSE

        dom_w = dom_w_in
        dom_eo = dom_e_in
        dom_s = dom_s_in
        dom_no = dom_n_in

      END IF ! dom_s == 0

!   Corner size ix, jy  do not have to be equal
!    Max values  ix = 4 ,    jy = 4

      IF ( blockx_in == 0 ) THEN

!   Domain printing corner size need to be set from
!   RUN_Diagnostics  NAMELIST.  Defaulting to 4

        ix = 4
        jy = 4

      ELSE

        ix = blockx_in
        jy = blocky_in

      END IF ! blockx_in == 0

      IF (lhook) CALL dr_hook('INIT_BLOCK4_PR',zhook_out,zhook_handle)

      RETURN
      END SUBROUTINE init_block4_pr
