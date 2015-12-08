! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
MODULE c_lowcld_mod

IMPLICIT NONE
! C_LOWCLD contains declarations for ceiling and threshold
! constants for lowest cloud layer diagnostics.
      ! Max height asl for 'low' cloud (1000ft)
      REAL,PARAMETER:: STR_CEIL=1000.0

      ! Cloud fraction threshold for low cloud.
      REAL,PARAMETER:: CLOUD_THRESHOLD=0.05
! C_LOWCLD end

END MODULE c_lowcld_mod
