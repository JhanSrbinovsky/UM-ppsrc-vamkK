! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module UKMO_GRIB_MOD --------------------------------------
!
! Purpose : Read GRIB 1 & 2 files
!
! Coding Standard : UM documentation paper no. 3
!
! Documentation : None
!
!----------------------------------------------------------------
!
!*L Arguments
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc
MODULE ukmo_grib_mod





USE missing_data_mod, Only: RMDI

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim





IMPLICIT NONE

PRIVATE

PUBLIC :: ukmo_decode

CONTAINS
SUBROUTINE ukmo_decode(fp_data, fp_work, len_fp, num_fp,                     & 
                  vert_coords, len_vert, num_vert,                           &
                  bitmap, len_bitmap, num_bitmap,                            &
                  quasi, len_q, num_q,                                       &
                  width, word_size,                                          &
                  block_0, block_1, block_2, block_3, block_4, block_r,      &
                  mesg, len_mesg, posn, word, off, error,                    & 
                  work_int1, work_int2, work_rel, err_unit, msglvl)

IMPLICIT NONE

! Scalar arguments
INTEGER, INTENT(IN) :: len_fp
INTEGER, INTENT(IN) :: len_vert
INTEGER, INTENT(IN) :: len_bitmap
INTEGER, INTENT(IN) :: len_q
INTEGER, INTENT(IN) :: word_size
INTEGER, INTENT(IN) :: len_mesg
INTEGER, INTENT(IN) :: err_unit
INTEGER, INTENT(IN) :: msglvl

! Array arguments
INTEGER, INTENT(IN) :: mesg(len_mesg)
INTEGER, INTENT(IN) :: posn(4)

! Scalar arguments
INTEGER, INTENT(INOUT) :: word
INTEGER, INTENT(INOUT) :: off

! Array arguments
INTEGER, INTENT(INOUT) :: block_4(2)

! Scalar arguments
INTEGER, INTENT(OUT) :: num_fp
INTEGER, INTENT(OUT) :: num_vert
INTEGER, INTENT(OUT) :: num_bitmap
INTEGER, INTENT(OUT) :: num_q
INTEGER, INTENT(OUT) :: width

! Array arguments
REAL,    INTENT(OUT) :: fp_data(len_fp)
REAL,    INTENT(OUT) :: fp_work(len_fp)
REAL,    INTENT(OUT) :: vert_coords(len_vert)
INTEGER, INTENT(OUT) :: bitmap(len_bitmap)
INTEGER, INTENT(OUT) :: quasi(len_q)
INTEGER, INTENT(OUT) :: block_0(4)
INTEGER, INTENT(OUT) :: block_1(*)
INTEGER, INTENT(OUT) :: block_2(20)
INTEGER, INTENT(OUT) :: block_3(2)
REAL,    INTENT(OUT) :: block_r(20)
INTEGER, INTENT(OUT) :: work_int1(*)
INTEGER, INTENT(OUT) :: work_int2(*)
REAL,    INTENT(OUT) :: work_rel(*)

! Error status (0 good, 2 warning, 3 error)
INTEGER, INTENT(INOUT) :: error
!# DEPENDS ON: decode
CALL decode(fp_data, fp_work, len_fp, num_fp,                          & 
            vert_coords, len_vert, num_vert,                           &
            bitmap, len_bitmap, num_bitmap,                            &
            quasi, len_q, num_q,                                       &
            width, word_size,                                          &
            block_0, block_1, block_2, block_3, block_4, block_r,      &
            mesg, len_mesg, posn, word, off, error,                    & 
            work_int1, work_int2, work_rel, err_unit, msglvl)
IF (lhook) CALL dr_hook('DECODE',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukmo_decode
END MODULE ukmo_grib_mod
