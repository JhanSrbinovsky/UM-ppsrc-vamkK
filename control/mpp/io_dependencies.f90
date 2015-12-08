! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! A module to contain downstream data dependencies
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

MODULE io_dependencies

! Older c layer

! DEPENDS ON : portio2a
! DEPENDS ON : pio_io_timer
! DEPENDS ON : print_from_c


! Newer C99 based layer
   
! DEPENDS ON : pio_data_conv
! DEPENDS ON : portutils
  IMPLICIT NONE

END MODULE io_dependencies
