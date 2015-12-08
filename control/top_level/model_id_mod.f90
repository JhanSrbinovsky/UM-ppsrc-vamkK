! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

! Module for the model configuration identifier
! this can be used to explicitly set LBEXP in UMDP F3
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

MODULE model_id_mod

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

INTEGER   ::    itab = imdi ! model configuration identifier 
                            ! if set then LBEXP will be set to this.

NAMELIST/configid/itab

END MODULE model_id_mod
