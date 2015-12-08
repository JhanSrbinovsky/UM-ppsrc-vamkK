! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE c_bm_chm_mod

IMPLICIT NONE
! C_BM_CHM start
! Contains constants required for biomass smoke conversion and
! diffusional scavenging by cloud droplets.
!
!     air parcel lifetime in cloud
      REAL,PARAMETER:: CLOUDTAU = 1.08E4            ! secs (=3 hours)

!     timescale for suspended aerosol to evaporate
      REAL,PARAMETER:: EVAPTAU = 300.0              ! secs  (=5 mins)

!     timescale for accumulation mode particles
      REAL,PARAMETER:: NUCTAU = 30.0                ! secs

!     Cloud liquid water threshold for nucleation scavenging to occur.
      REAL,PARAMETER:: THOLD = 1.0E-8               ! kg/kg

! C_BM_CHM end

END MODULE c_bm_chm_mod
