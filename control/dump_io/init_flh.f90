! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE INIT_FLH --------------------------------------
!
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    System component: R30
!
!    Purpose:
!             Initialises the fixed length header to IMDI except
!             for all array dimensions which are set to 1.
!
!    Documentation:
!             Unified Model Documentation Paper No F3
!
!  ------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Dump I/O
SUBROUTINE init_flh (fixhd,len_fixhd)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER                                                           &
 len_fixhd                                                        &
                  ! IN    Length of fixed length header
,fixhd(len_fixhd) ! INOUT Fixed length header

! Local arrays:------------------------------------------------
! None
! -------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
! Local variables:---------------------------------------------
INTEGER j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! -------------------------------------------------------------

! 1.0 Initialise to IMDI
IF (lhook) CALL dr_hook('INIT_FLH',zhook_in,zhook_handle)
DO j = 1,len_fixhd
  fixhd(j) = imdi
END DO

! 2.0 Set all array dimensions to 1
fixhd(101) = 1     !  Integer Constants
fixhd(106) = 1     !  Real Constants
fixhd(111) = 1     !  1st dim - Level dependent constants
fixhd(112) = 1     !  2nd dim - Level dependent constants
fixhd(116) = 1     !  1st dim - Row dependent constants
fixhd(117) = 1     !  2nd dim - Row dependent constants
fixhd(121) = 1     !  1st dim - Column dependent constants
fixhd(122) = 1     !  2nd dim - Column dependent constants
fixhd(126) = 1     !  1st dim - Field of constants
fixhd(127) = 1     !  2nd dim - Field of constants
fixhd(131) = 1     !  Extra constants
fixhd(136) = 1     !  Temp History file
fixhd(141) = 1     !  Compressed field Index 1
fixhd(143) = 1     !  Compressed field Index 2
fixhd(145) = 1     !  Compressed field Index 3
fixhd(151) = 1     !  1st dim - Lookup Table
fixhd(152) = 1     !  2nd dim - Lookup Table
fixhd(161) = 1     !  Data

IF (lhook) CALL dr_hook('INIT_FLH',zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_flh
