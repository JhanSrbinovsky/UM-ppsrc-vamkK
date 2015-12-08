! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Control variables

!  Description:  

! Fortran module to populate hard wired variables in namelist cntall.

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90

MODULE nlstcall_pop_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE pop_nlstcall(type_letter_1,type_letter_3)

USE chsunits_mod, ONLY : nunits

  IMPLICIT NONE

  INTEGER :: i_arr
  CHARACTER(len=1) :: tl1, tl3

  CHARACTER(len=1) :: type_letter_1(20:NUNITS)
  CHARACTER(len=1) :: type_letter_3(20:NUNITS)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!
  !!! type_letter_1, type_letter_3
  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i_arr=20,nunits,1
    !logic for type_letter_1
    IF (i_arr==90 .OR. i_arr == 140 .OR. i_arr == 141 .OR. i_arr == 142 &
       .OR. i_arr == 143 .OR. i_arr == 144 .OR. i_arr == 145 &
       .OR. i_arr == 146 .OR. i_arr == 174) THEN
      tl1 = 'b'
    ELSE IF (i_arr == 150 .OR. i_arr == 164) THEN
      tl1 = 'c'
    ELSE
      tl1 = 'p'
    END IF
    type_letter_1(i_arr) = tl1
    !end logic for type_letter_1

    !logic for type_letter_3
    IF (i_arr == 60 .OR. i_arr == 140 .OR. i_arr == 150) THEN
      tl3 = 'a'
    ELSE IF (i_arr == 61 .OR. i_arr == 141 .OR. i_arr == 164) THEN
      tl3 = 'b'
    ELSE IF (i_arr == 62 .OR. i_arr == 142) THEN
      tl3 = 'c'
    ELSE IF (i_arr == 63 .OR. i_arr == 143) THEN
      tl3 = 'd'
    ELSE IF (i_arr == 64 .OR. i_arr == 144) THEN
      tl3 = 'e'
    ELSE IF (i_arr == 65 .OR. i_arr == 145) THEN
      tl3 = 'f'
    ELSE IF (i_arr == 66 .OR. i_arr == 146) THEN
      tl3 = 'g'
    ELSE IF (i_arr == 67 .OR. i_arr == 147) THEN
      tl3 = 'h'
    ELSE IF (i_arr == 68) THEN
      tl3 = 'i'
    ELSE IF (i_arr == 69) THEN
      tl3 = 'j'
    ELSE IF (i_arr == 100) THEN
      tl3 = 'l'
    ELSE IF (i_arr == 101) THEN
      tl3 = 'm'
    ELSE IF (i_arr == 102) THEN
      tl3 = 'n'
    ELSE IF (i_arr == 151) THEN
      tl3 = 'k'
    ELSE
      tl3 = ' '
    END IF
    type_letter_3(i_arr) = tl3
    !end logic for type_letter_3
  END DO

END SUBROUTINE pop_nlstcall

END MODULE nlstcall_pop_mod
