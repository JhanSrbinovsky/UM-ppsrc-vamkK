! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Frames

MODULE frames_mod

IMPLICIT NONE

TYPE frames_lookup_type
  INTEGER, POINTER :: frames_lookup(:,:)
END TYPE frames_lookup_type

! Arrays to hold all headers for every frames output file.
! Allocated in in_intf.
INTEGER, ALLOCATABLE, TARGET :: Frames_FixHds (:,:)
INTEGER, ALLOCATABLE, TARGET :: Frames_IntHds (:,:)
REAL   , ALLOCATABLE, TARGET :: Frames_RealHds(:,:)
REAL   , ALLOCATABLE, TARGET :: Frames_levdepcs(:,:)
TYPE(frames_lookup_type), ALLOCATABLE, TARGET :: Frames_lookups(:)

INTEGER, ALLOCATABLE :: frames_count(:)
INTEGER,SAVE :: Keep_Pack = -1
INTEGER :: frames_total_levels

END MODULE frames_mod
