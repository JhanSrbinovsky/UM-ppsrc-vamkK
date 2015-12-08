! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE c_skeb2_mod

IMPLICIT NONE
! C_SKEB2 start
! Define, nblock, nfacts for use in skeb2 routines forcing.F90, 
! rpassm.F90, set99.F90 and fourier.F90
 INTEGER, PARAMETER:: nblock=4096
 INTEGER, PARAMETER:: nfacts=10 
! C_SKEB2 end

END MODULE c_skeb2_mod
