! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE c_ocff_chm_mod

IMPLICIT NONE
! C_OCFF_CHM start
! Contains constants required for fossil fuel organic carbon
! conversion and diffusional scavenging by cloud droplets.
!
!     air parcel lifetime in cloud
      REAL,PARAMETER:: CLOUD_TAU = 1.08E4           ! secs (=3 hours)

!     timescale for suspended aerosol to evaporate
      REAL,PARAMETER:: EVAP_TAU = 300.0             ! secs  (=5 mins)

!     timescale for accumulation mode particles
      REAL,PARAMETER:: NUC_TAU = 30.0               ! secs

!     Cloud liquid water threshold for nucleation scavenging to occur.
      REAL,PARAMETER:: THOLD = 1.0E-8               ! kg/kg

! C_OCFF_CHM end

END MODULE c_ocff_chm_mod
