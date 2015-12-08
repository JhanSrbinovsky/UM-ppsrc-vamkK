
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Version 1B of vegetation section: JULES-compatible
!
! Subroutine Interface:
      SUBROUTINE VEG_IC(LAND_PTS,LAND_INDEX,NTILES,CAN_MODEL            &
     &,                 A_STEP,ASTEPS_SINCE_TRIFFID                     &
     &,                 PHENOL_PERIOD,TRIFFID_PERIOD                    &
     &,                 L_PHENOL,L_TRIFFID,L_TRIF_EQ                    &
     &,                 ATIMESTEP,FRAC_DISTURB,SATCON                   &
     &,                 G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC                 &
     &,                 RESP_S_AC,RESP_W_AC                             &
     &,                 CS,FRAC,LAI,HT                                  &
      ,                 CATCH_S,CATCH_T,INFIL_T,Z0_T,z0h_t              &
     &,                 C_VEG,CV,LIT_C,LIT_C_MN,G_LEAF_DAY,G_LEAF_PHEN  &
     &,                 LAI_PHEN,G_LEAF_DR_OUT,NPP_DR_OUT,RESP_W_DR_OUT &
     &,                 RESP_S_DR_OUT                                   &
     &                  )

  USE nstypes
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE
!
! Description:
!   Updates Leaf Area Index for Plant Functional Types (PFTs) and uses
!   this to derive new vegetation parameters for PFTs along with gridbox
!   mean values where appropriate.
!
! Method:
!   Calls PHENOL which models phenology and updates Leaf Area Index
!   (LAI), then passes new LAI into SPARM along with canopy height
!   and fractional cover of Plant Functional Types.  SPARM uses this to
!   derive the vegetation parameters for each PFT, and also derives
!   gridbox means where this is required.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Vegetation
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

      INTEGER                                                           &
     & LAND_PTS                                                         &
                             ! IN Number of land points to be processed.
     &,NTILES                                                           &
                             ! IN Number of land-surface tiles.
     &,CAN_MODEL                                                        &
                                    ! IN Swith for thermal vegetation
     &,A_STEP                                                           &
                             ! IN Atmospheric timestep number.
     &,ASTEPS_SINCE_TRIFFID                                             &
                             ! INOUT Number of atmospheric
!                            !       timesteps since last call
!                            !       to TRIFFID.
     &,PHENOL_PERIOD                                                    &
                             ! IN Phenology period (days).
     &,TRIFFID_PERIOD        ! IN TRIFFID period (days).


      INTEGER                                                           &
     & LAND_INDEX(LAND_PTS)         ! IN I=LAND_INDEX(L) => the Ith
!                                   !    point in P_FIELD is the Lth
!                                   !    land point.
      LOGICAL                                                           &
     & L_PHENOL                                                         &
                                    ! IN .T. for interactive leaf
!                                   !    phenology.
     &,L_TRIFFID                                                        &
                                    ! IN .T. for interactive vegetation.
     &,L_TRIF_EQ                    ! IN .T. for vegetation equilibrium.

      REAL                                                              &
     & ATIMESTEP                                                        &
                                    ! IN Atmospheric timestep (s).
     &,FRAC_DISTURB(LAND_PTS)                                           &
                                    ! IN Fraction of gridbox in which
!                                   !    vegetation is disturbed.
     &,SATCON(LAND_PTS)                                                 &
                                    ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).
     &,G_LEAF_AC(LAND_PTS,NPFT)                                         &
                                    ! INOUT Accumulated leaf turnover
!                                   !       rate.
     &,G_LEAF_PHEN_AC(LAND_PTS,NPFT)                                    &
                                    ! INOUT Accumulated leaf turnover
!                                   !         rate including phenology.
     &,NPP_AC(LAND_PTS,NPFT)                                            &
                                    ! INOUT Accumulated NPP (kg C/m2).
     &,RESP_W_AC(LAND_PTS,NPFT)                                         &
                                    ! INOUT Accumulated wood respiration
!                                   !       (kg C/m2).
     &,RESP_S_AC(LAND_PTS)                                              &
                                    ! INOUT Accumulated soil respiration
!                                   !       (kg C/m2).
     &,CS(LAND_PTS)                                                     &
                                    ! INOUT Soil carbon content
!                                   !       (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                    ! INOUT Fractions of surface types.
     &,LAI(LAND_PTS,NPFT)                                               &
                                    ! INOUT LAI of plant functional
!                                   !       types.
     &,HT(LAND_PTS,NPFT)                                                &
                                    ! INOUT Height of plant functional
!                                   !       types (m).
     &,CATCH_S(LAND_PTS,NTILES)                                         &
                                    ! OUT Canopy snow capacity (kg/m2).
     &,CATCH_T(LAND_PTS,NTILES)                                         &
                                    ! OUT Canopy capacity for tiles
!                                   !     (kg/m2).
     &,INFIL_T(LAND_PTS,NTILES)                                         &
                                    ! OUT Maximum surface infiltration
!                                   !     rate for tiles (kg/m2/s).
     &,G_LEAF_DAY(LAND_PTS,NPFT)                                        &
                                    ! OUT Mean leaf turnover rate for
!                                     !     input to PHENOL (/360days).
     &,G_LEAF_PHEN(LAND_PTS,NPFT)                                       &
                                    ! OUT Mean leaf turnover rate over
!                                     !     phenology period (/360days).
     &,LAI_PHEN(LAND_PTS,NPFT)                                          &
                                    ! OUT LAI of PFTs after phenology.
     &,Z0_T(LAND_PTS,NTILES)                                            &
                                    ! OUT Roughness length for tiles (m)
      ,z0h_t(land_pts,ntiles)                                           &
                                    ! OUT Thermal roughness length for 
                                    !     tiles (m)
     &,C_VEG(LAND_PTS,NPFT)                                             &
                                    ! OUT Total carbon content of
!                                   !     the vegetation (kg C/m2).
     &,CV(LAND_PTS)                                                     &
                                    ! OUT Gridbox mean vegetation
!                                   !     carbon (kg C/m2).
     &,LIT_C(LAND_PTS,NPFT)                                             &
                                    ! OUT Carbon Litter
!                                   !     (kg C/m2/360days).
     &,LIT_C_MN(LAND_PTS)           ! OUT Gridbox mean carbon litter
!                                   !     (kg C/m2/360days).

! External routines called
      EXTERNAL VEG

! DUMMY ARGUMENTS
      REAL                                                              &
     & G_LEAF_DR_OUT(LAND_PTS,NPFT)                                     &
     &,NPP_DR_OUT(LAND_PTS,NPFT)                                        &
     &,RESP_W_DR_OUT(LAND_PTS,NPFT)                                     &
     &,RESP_S_DR_OUT(LAND_PTS)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('VEG_IC',zhook_in,zhook_handle)

! DEPENDS ON: veg
      CALL VEG( LAND_PTS,NTILES,CAN_MODEL                               &
     &,         A_STEP,PHENOL_PERIOD,L_PHENOL                           &
     &,         ATIMESTEP,SATCON                                        &
     &,         G_LEAF_AC,FRAC,LAI,HT                                   &
      ,         CATCH_S,CATCH_T,INFIL_T,Z0_T,z0h_t                      &
     &,         G_LEAF_DAY,G_LEAF_PHEN,LAI_PHEN                         &
     &         )

      IF (lhook) CALL dr_hook('VEG_IC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE VEG_IC
