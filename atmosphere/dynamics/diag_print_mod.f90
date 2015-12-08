! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Declares diagnostic printing control
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: dynamics

MODULE diag_print_mod

IMPLICIT NONE
! diagnostic printing control
LOGICAL :: L_print_diag
LOGICAL :: L_print_L2helm
LOGICAL :: L_print_allpe
LOGICAL :: L_flush

INTEGER :: norm_start_lev
INTEGER :: norm_end_lev

END MODULE diag_print_mod
