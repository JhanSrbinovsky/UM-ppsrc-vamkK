! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE mphys_diags_mod

! Description:
! Holds diagnostic switches and values required by the large-scale
! precipitation scheme that have been integrated over the particle size
! distribution
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

  IMPLICIT NONE

  SAVE

! Switches for aggregate fraction and points used
  LOGICAL :: l_aggfr_diag = .FALSE. ! Aggregate Fraction

  LOGICAL :: l_point_diag = .FALSE. ! Points used

! Microphysical process rate diagnostic logical switches                  
  LOGICAL :: l_psdep_diag = .FALSE.               
                       ! Deposition of vapour to snow agg.
  LOGICAL :: l_psaut_diag = .FALSE.                                      
                       ! Autoconversion of aggregates from cry
  LOGICAL :: l_psacw_diag = .FALSE. 
                       ! Accretion of liq. water by snow agg.
  LOGICAL :: l_psacr_diag = .FALSE.   
                      ! Collection of rain by snow aggregates
  LOGICAL :: l_psaci_diag = .FALSE.   
                       ! Collection of ice crystals by agg.
  LOGICAL :: l_psmlt_diag = .FALSE.   
                       ! Melting of snow aggregates
  LOGICAL :: l_psmltevp_diag = .FALSE.
                       ! Evaporation of melting aggregates
  LOGICAL :: l_praut_diag = .FALSE.   
                       ! Autoconversion of cloud drops to rain
  LOGICAL :: l_pracw_diag = .FALSE.   
                       ! Accretion of liq. water by rain
  LOGICAL :: l_prevp_diag = .FALSE.   
                       ! Evaporation of rain
  LOGICAL :: l_pgaut_diag = .FALSE.   
                       ! Autoconversion of graupel from agg.
  LOGICAL :: l_pgacw_diag = .FALSE.   
                       ! Accretion of liq. water by graupel
  LOGICAL :: l_pgacs_diag = .FALSE.   
                       ! Collection of snow agg. by graupel
  LOGICAL :: l_pgmlt_diag = .FALSE.      
                       ! Melting of graupel
  LOGICAL :: l_pifrw_diag = .FALSE.   
                       ! Homogeneous freezing nucleation
  LOGICAL :: l_piprm_diag = .FALSE.   
                       ! Heterogeneous (primary) nucleation
  LOGICAL :: l_pidep_diag = .FALSE.   
                       ! Deposition of vapour to ice crystals
  LOGICAL :: l_piacw_diag = .FALSE.   
                       ! Accretion of liq. water by ice cry.
  LOGICAL :: l_piacr_diag = .FALSE.   
                       ! Collection of rain by ice crystals
  LOGICAL :: l_pimlt_diag = .FALSE.   
                       ! Melting of ice crystals
  LOGICAL :: l_pimltevp_diag = .FALSE.
                       ! Evaporation of melting ice crystals
  LOGICAL :: l_pifall_diag = .FALSE.  
                       ! Sedimentation of ice crystals
  LOGICAL :: l_psfall_diag = .FALSE.  
                       ! Sedimentation of aggregates
  LOGICAL :: l_prfall_diag = .FALSE.  
                       ! Sedimentation of rain
  LOGICAL :: l_pgfall_diag = .FALSE.  
                       ! Sedimentation of graupel
  LOGICAL :: l_plset_diag = .FALSE.   
                       ! Droplet settling of liquid water
  LOGICAL :: l_plevpset_diag = .FALSE.
                       ! Evaporated settled droplets
  LOGICAL :: l_pifrr_diag = .FALSE.
                       ! Homogeneous freezing of rain

  REAL, ALLOCATABLE ::                                                         &
    psdep(:,:,:),  psaut(:,:,:),   psacw(:,:,:),  psacr(:,:,:), psaci(:,:,:),  &
    psmlt(:,:,:),  psmltevp(:,:,:), praut(:,:,:),  pracw(:,:,:), prevp(:,:,:), &
    frac_agg(:,:,:), pgaut(:,:,:),  pgacw(:,:,:),  pgacs(:,:,:), pgmlt(:,:,:), &
    pifrw(:,:,:),    pifrr(:,:,:),  piprm(:,:,:),  pidep(:,:,:), piacw(:,:,:), &
    piacr(:,:,:), pimlt(:,:,:), pimltevp(:,:,:), pifall(:,:,:), psfall(:,:,:), &
    prfall(:,:,:),  pgfall(:,:,:), plset(:,:,:), plevpset(:,:,:)

  LOGICAL, ALLOCATABLE :: mphys_pts(:,:,:)


END MODULE mphys_diags_mod
