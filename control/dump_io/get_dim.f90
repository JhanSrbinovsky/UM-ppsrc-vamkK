! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine GET_DIM
!
!
!    Programming Standard : UM Doc Paper No 3
!
!    Project Task : P3
!
!    Purpose : Get dimensions of data set components
!              from fixed length header.
!
!     ------------------------------------------------------------------
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: Dump I/O
SUBROUTINE get_dim (fixhd,len_fixhd,                              &
                    len_inthd,len_realhd,                         &
                    len1_levdepc,len2_levdepc,                    &
                    len1_rowdepc,len2_rowdepc,                    &
                    len1_coldepc,len2_coldepc,                    &
                    len1_flddepc,len2_flddepc,                    &
                    len_extcnst,len_dumphist,                     &
                    len_cfi1,len_cfi2,len_cfi3,                   &
                    len1_lookup,len2_lookup,                      &
                    len_data)

! ----------------------------------------------------------------------

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER                                                           &
                 !  Dimension of :-
  len_fixhd                                                       &
                 !  Fixed length header
 ,len_inthd                                                       &
                 !  Integer header
 ,len_realhd                                                      &
                 !  Real header
 ,len1_levdepc                                                    &
                 !  Level dependent constants (1st)
 ,len2_levdepc                                                    &
                 !  Level dependent constants (2nd)
 ,len1_rowdepc                                                    &
                 !  Rows  dependent constants (1st)
 ,len2_rowdepc                                                    &
                 !  Rows  dependent constants (2nd)
 ,len1_coldepc                                                    &
                 !  Col   dependent constants (1st)
 ,len2_coldepc                                                    &
                 !  Col   dependent constants (2nd)
 ,len1_flddepc                                                    &
                 !  Field dependent constants (1st)
 ,len2_flddepc                                                    &
                 !  Field dependent constants (2nd)
 ,len_extcnst                                                     &
                 !  Extra constants
 ,len_dumphist                                                    &
                 !  Dump history
 ,len_cfi1                                                        &
                 !  Compressed field index 1
 ,len_cfi2                                                        &
                 !  Compressed field index 2
 ,len_cfi3                                                        &
                 !  Compressed field index 3
 ,len1_lookup                                                     &
                 !  Look up table (1st)
 ,len2_lookup                                                     &
                 !  Look up table (2nd)
 ,len_data       !  Data section

INTEGER  fixhd(len_fixhd)   ! IN  Fixed length header

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('GET_DIM',zhook_in,zhook_handle)
len_inthd    = fixhd(101)
len_realhd   = fixhd(106)
len1_levdepc = fixhd(111)
len2_levdepc = fixhd(112)
len1_rowdepc = fixhd(116)
len2_rowdepc = fixhd(117)
len1_coldepc = fixhd(121)
len2_coldepc = fixhd(122)
len1_flddepc = fixhd(126)
len2_flddepc = fixhd(127)
len_extcnst  = fixhd(131)
len_dumphist = fixhd(136)
len_cfi1     = fixhd(141)
len_cfi2     = fixhd(143)
len_cfi3     = fixhd(145)
len1_lookup  = fixhd(151)
len2_lookup  = fixhd(152)
len_data     = fixhd(161)

IF (lhook) CALL dr_hook('GET_DIM',zhook_out,zhook_handle)
RETURN
END SUBROUTINE get_dim
