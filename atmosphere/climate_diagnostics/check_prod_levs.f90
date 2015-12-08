! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Checks for shared pressure levels for products
!
! Subroutine Interface:
MODULE check_prod_levs_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE check_prod_levs(                                       &
     &  num_stash_levels,stash_levels,ni,                               &
     &  in1_p_levs,in1_press,                                           &
     &  in2_p_levs,in2_press,                                           &
     &  out_p_levs,out_ind)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Description:
!  Check to see if products required in section 30 share pressure
!  levels with the fields required to produce them. Outputs an array
!  indicating the shared levels
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Climate Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Subroutine arguments
      integer num_stash_levels
      integer out_p_levs,in1_p_levs,in2_p_levs
      integer out_ind(num_stash_levels*2)

      real                                                              &
     &  in1_press(num_stash_levels),                                    &
     &  in2_press(num_stash_levels)
!Stash variables
      integer stash_levels(num_stash_levels+1,*)
      integer ni
! local variables

      integer i,j,k
      real press_levs(num_stash_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('CHECK_PROD_LEVS',zhook_in,zhook_handle)
      out_p_levs=stash_levels(1,ni)
      do k =1,out_p_levs
        press_levs(k)=stash_levels(k+1,ni)/1000.0
        out_ind(k)=0
        out_ind(out_p_levs+k)=0
        do i=1,in1_p_levs
          if (press_levs(k) == in1_press(i)) then
            out_ind(k)=i
          endif
        enddo
        do i=1,in2_p_levs
          if (press_levs(k) == in2_press(i)) then
            out_ind(out_p_levs+k)=i
          endif
        enddo
      enddo
      IF (lhook) CALL dr_hook('CHECK_PROD_LEVS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE check_prod_levs
END MODULE check_prod_levs_mod
