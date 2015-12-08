! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine HYD_INTCTL------------------------------------------
!
!   Level 3 control routine
!
!   Purpose: Calls HYDROL to calculate and add hydrology increments.
!            Avoids the need for *IF DEF around calls to different
!            versions of HYDROL.
!            Multilayer Version.

!   Programming standard : unified model documentation paper No 3
!
!   System components covered : P25
!
!   System task : P0
!
!   Documentation: Unified Model documentation paper P0
!                  version No 11 dated (26/11/90)
!  END -----------------------------------------------------------------
!   Arguments
! Subroutine Interface:
      SUBROUTINE HYD_INTCTL (                                           &
     &                   LICE_PTS,LICE_INDEX,SOIL_PTS,SOIL_INDEX,       &
     &                   NPNTS,NSHYD,B,CON_RAIN,CON_SNOW,               &
     &                   EXT,HCAP,HCON,LS_RAIN,                         &
     &                   LS_SNOW,SATCON,SATHH,                          &
     &                   SURF_HT_FLUX,                                  &
     &                   TIMESTEP,V_SAT,V_WILT,                         &
     &                   CAN_WCNT_GB,HF_SNOW_MELT,STF_HF_SNOW_MELT,     &
     &                   SMCL,STHF,STHU,TSOIL,LYING_SNOW,               &
     &                   INFIL,SMC,SNOW_MELT,SNOMLT_SUB_HTF,            &
     &                   STF_SUB_SURF_ROFF,SUB_SURF_ROFF,SURF_ROFF,     &
     &                   TOT_TFALL,                                     &

!  Add inland basin outflow to arguments
     &                   INLANDOUT_ATM,L_INLAND,                        &
! Additional variables for MOSES II
     &                   NTILES,TILE_PTS,TILE_INDEX,L_SNOW_ALBEDO,      &
     &                   CAN_MODEL,                                     &
     &                   CAN_CPY,E_CANOPY,INFIL_TILE,                   &
     &                   MELT_TILE,TSTAR_TILE,FRAC,CAN_WCNT,            &
     &                   RGRAIN,SNOW_TILE,                              &
     &                   CATCH_SNOW,SNOW_GRND,                          &

! Additional variables required for large-scale hydrology:
     &                   L_TOP,L_PDM,FEXP,GAMTOT,TI_MEAN,TI_SIG,CS_CH4, &
     &                   DUN_ROFF,FSAT,FWETL,QBASE,QBASE_ZW,ZW,DRAIN,   &
     &                   STHZW,A_FSAT,C_FSAT,A_FWET,C_FWET,FCH4_WETL,   &
     &                   dim_cs1,L_SOIL_SAT_DOWN,l_triffid,             &

     & LTIMER                                                           &
     & )

! Import additional variables for JULES
      USE prognostics, ONLY:                                            &
                                nsnow                                   &
                              , snowdepth

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE
!
! Description:
!     Surface hydrology module which also updates the
!     sub-surface temperatures. Includes soil water phase
!     changes and the effect of soil water and ice on the
!     thermal and hydraulic characteristics of the soil.
!     This version is for use with MOSES II land surface
!     scheme. Calls the following:

!     SFSNOW - to calculate the sub-surface snowmelt
!              and update the lying snow amount    (Essery, 1/95)
!
!     SURF_HYD - to calculate canopy interception and
!                surface runoff         (Allen-Bett, Gregory, 90)
!
!     SOIL_HYD - to update the layer soil moisture contents
!                and calculate the drainage            (Cox 6/95)
!
!     SOIL_HTC - to update the soil layer temperatures and the
!                layer ice contents                    (Cox 6/95)
!
!     ICE_HTC - to update the layer temperatures for land ice
!                                                      (Cox 10/95)
!
!     SOIL_MC - to diagnose the soil moisture in the top metre
!                                                    (Essery 7/97)
!
! Documentation : UM Documentation Paper 25
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Land Surface
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & LICE_PTS                                                         &
                           ! IN Number of land ice points.
     &,NPNTS                                                            &
                           ! IN Number of gridpoints.
     &,NSHYD                                                            &
                           ! IN Number of soil moisture levels.
     &,SOIL_PTS                                                         &
                           ! IN Number of soil points.
     &,DIM_CS1             ! IN soil carbon dimension

      REAL                                                              &
     & TIMESTEP            ! IN Model timestep (s).

      LOGICAL                                                           &
       LTIMER                                                           &
                           ! Logical switch for TIMER diags
      ,L_TRIFFID           ! IN Indicates whether TRIFFID
!                          !    in use.

      LOGICAL                                                           &
     & STF_HF_SNOW_MELT                                                 &
                           ! IN Stash flag for snowmelt heat flux.
     &,STF_SUB_SURF_ROFF   ! IN Stash flag for sub-surface runoff.


!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & LICE_INDEX(NPNTS)                                                &
                           ! IN Array of land ice points.
     &,SOIL_INDEX(NPNTS)   ! IN Array of soil points.

      REAL                                                              &
     & B(NPNTS,NSHYD)                                                   &
                           ! IN Clapp-Hornberger exponent.
     &,CON_RAIN(NPNTS)                                                  &
                           ! IN Convective rain (kg/m2/s).
     &,CON_SNOW(NPNTS)                                                  &
                           ! IN Convective snowfall (kg/m2/s).
     &,EXT(NPNTS,NSHYD)                                                 &
                           ! IN Extraction of water from each soil
!                          !    layer (kg/m2/s).
     &,HCAP(NPNTS,NSHYD)                                                &
                           ! IN Soil heat capacity (J/K/m3).
     &,HCON(NPNTS,0:NSHYD)                                              &
                           ! IN Soil thermal conductivity (W/m/K).
     &,LS_RAIN(NPNTS)                                                   &
                           ! IN Large-scale rain (kg/m2/s).
     &,LS_SNOW(NPNTS)                                                   &
                           ! IN Large-scale snowfall (kg/m2/s).
     &,SATCON(NPNTS,0:NSHYD)                                            &
                           ! IN Saturated hydraulic conductivity
!                          !    (kg/m2/s).
     &,SATHH(NPNTS,NSHYD)                                               &
                           ! IN Saturated soil water pressure (m).
     &,SURF_HT_FLUX(NPNTS)                                              &
                           ! IN Net downward surface heat flux (W/m2)
     &,V_SAT(NPNTS,NSHYD)                                               &
                           ! IN Volumetric soil moisture
!                          !    concentration at saturation
!                          !    (m3 H2O/m3 soil).
     &,V_WILT(NPNTS,NSHYD) ! IN Volumetric soil moisture
!                          !    concentration below which
!                          !    stomata close (m3 H2O/m3 soil).
!
!   Array arguments with intent(INOUT) :
!
      REAL                                                              &
     & SMCL(NPNTS,NSHYD)                                                &
                           ! INOUT Soil moisture content of each
!                          !       layer (kg/m2).
     &,STHF(NPNTS,NSHYD)                                                &
                           ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,STHU(NPNTS,NSHYD)                                                &
                           ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,TSOIL(NPNTS,NSHYD)  ! INOUT Sub-surface temperatures (K).


!   Array arguments with intent(OUT) :
      REAL                                                              &
     & CAN_WCNT_GB(NPNTS)                                               &
                            ! OUT Gridbox canopy water content (kg/m2).
     &,HF_SNOW_MELT(NPNTS)                                              &
                            ! OUT Gridbox snowmelt heat flux (W/m2).
     &,INFIL(NPNTS)                                                     &
                            ! OUT Maximum surface infiltration
!                           !     rate (kg/m2/s).
     &,SMC(NPNTS)                                                       &
                            ! OUT Soil moisture in the top metre (kg/m2)
     &,SNOW_MELT(NPNTS)                                                 &
                            ! OUT Snowmelt (kg/m2/s).
     &,SNOMLT_SUB_HTF(NPNTS)                                            &
                            ! OUT Sub-surface snowmelt heat
!                           !     flux (W/m2).
     &,SUB_SURF_ROFF(NPNTS)                                             &
                            ! OUT Sub-surface runoff (kg/m2/s).
     &,SURF_ROFF(NPNTS)                                                 &
                            ! OUT Surface runoff (kg/m2/s).
     &,TOT_TFALL(NPNTS)                                                 &
                            ! OUT Total throughfall (kg/m2/s).
     &,LYING_SNOW(NPNTS)    ! OUT Gridbox snowmass (kg/m2).


! Additional variables for MOSES II
      INTEGER                                                           &
     & NTILES                                                           &
                           ! IN Number of tiles.
     &,TILE_PTS(NTILES)                                                 &
                           ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTILES)                                         &
!                          ! IN Index of tile points.
     &,CAN_MODEL           ! IN Switch for thermal vegetation canopy.

      LOGICAL                                                           &
     & L_SNOW_ALBEDO       ! IN Flag for prognostic snow albedo

      REAL                                                              &
     & CAN_CPY(NPNTS,NTILES)                                            &
                            !IN Canopy/surface capacity of
!                          !    land tiles (kg/m2).
     &,CATCH_SNOW(NPNTS,NTILES)                                         &
                           ! IN Coniferous canopy snow capacity
!                          !    (kg/m2).
     &,E_CANOPY(NPNTS,NTILES)                                           &
!                          ! IN Canopy evaporation from
!                          !    land tiles (kg/m2/s).
     &,INFIL_TILE(NPNTS,NTILES)                                         &
!                          ! IN Maximum surface infiltration
     &,MELT_TILE(NPNTS,NTILES)                                          &
!                          ! IN Snowmelt on tiles (kg/m2/s).
     &,TSTAR_TILE(NPNTS,NTILES)                                         &
!                          ! IN Surface temperature (K).
     &,FRAC(NPNTS,NTILES)                                               &
                           ! IN Tile fractions.
     &,CAN_WCNT(NPNTS,NTILES)                                           &
!                          ! INOUT Canopy water content for
!                          !       land tiles (kg/m2).
     &,RGRAIN(NPNTS,NTILES)                                             &
                           ! INOUT Snow grain size (microns).
     &,SNOW_TILE(NPNTS,NTILES)                                          &
!                          ! INOUT Snowmass on tiles (kg/m2).
     &,SNOW_GRND(NPNTS)                                                 &
                           ! INOUT Snow below canopy (kg/m2).
!  Declare inland basin outflow variable
     &,INLANDOUT_ATM(npnts)
!                          ! IN INLAND BASIN OUTFLOW kg/m2/s
      LOGICAL                                                           &
     & L_INLAND
                                  ! IN True if re-routing inland
                                  !   basin flow to soil moisture

! Additional variables required for large-scale hydrology:
       LOGICAL                                                          &
     & L_TOP                                                            &
                   ! IN Flag for TOPMODEL-based hydrology.
     &,L_PDM                                                            &
                   ! IN Flag for PDM.
     &,L_SOIL_SAT_DOWN
                   ! IN Direction of super-saturated soil moisture

      REAL                                                              &
     & FEXP(NPNTS)                                                      &
                           ! IN Decay factor in Sat. Conductivity
!                          !    in water table layer.
     &,GAMTOT(NPNTS)                                                    &
                           ! IN Integrated complete Gamma function.
     &,TI_MEAN(NPNTS)                                                   &
                           ! IN Mean topographic index.
     &,TI_SIG(NPNTS)                                                    &
                           ! IN Standard dev. of topographic index.
     &,CS_CH4(NPNTS)                                                    &
                           ! IN Effective soil carbon for
!                          !    CH4 wetlands calc (kg C/m2).
     &,FSAT(NPNTS)                                                      &
                           ! INOUT Surface saturation fraction.
     &,FWETL(NPNTS)                                                     &
                           ! INOUT Wetland fraction.
     &,ZW(NPNTS)                                                        &
                           ! INOUT Water table depth (m).
     &,STHZW(NPNTS)                                                     &
                           ! INOUT soil moist fract. in deep-zw layer.
     &,A_FSAT(NPNTS)                                                    &
                           ! IN Fitting parameter for Fsat in LSH model
     &,C_FSAT(NPNTS)                                                    &
                           ! IN Fitting parameter for Fsat in LSH model
     &,A_FWET(NPNTS)                                                    &
                           ! IN Fitting parameter for Fwet in LSH model
     &,C_FWET(NPNTS)                                                    &
                           ! IN Fitting parameter for Fwet in LSH model

     &,DUN_ROFF(NPNTS)                                                  &
                           ! OUT Dunne part of sfc runoff (kg/m2/s).
     &,QBASE(NPNTS)                                                     &
                           ! OUT Base flow (kg/m2/s).
     &,QBASE_ZW(NPNTS)                                                  &
                           ! OUT Base flow from ZW layer (kg/m2/s).
     &,DRAIN(NPNTS)                                                     & 
                           ! OUT Drainage out of nshyd'th level (kg/m2/s).
     &,FCH4_WETL(NPNTS)    ! OUT Wetland methane flux. (kg C/m2/s).

      Integer                                                           &
     & ErrorStatus            ! error reporting
      Character (Len=*)                                                 &
     & RoutineName
      Parameter (RoutineName='HYD_INTCTL')
      Character (Len=80) Cmessage

!

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('HYD_INTCTL',zhook_in,zhook_handle)

      IF(L_PDM.AND.L_TOP)THEN
        ErrorStatus=10
        Cmessage = 'Cannot have both PDM and LSH switched on'

        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      ENDIF

! DEPENDS ON: hydrol
      CALL HYDROL (                                                     &
! call to JULES2 routine
     &                   LICE_PTS,LICE_INDEX,SOIL_PTS,SOIL_INDEX,NSNOW, &
     &                   NPNTS,NSHYD,B,CAN_CPY,CON_RAIN,                &
     &                   E_CANOPY,EXT,HCAP,HCON,LS_RAIN,                &
     &                   SATCON,SATHH,SNOWDEPTH,                        &
     &                   SURF_HT_FLUX,TIMESTEP,                         &
     &                   V_SAT,V_WILT,                                  &
     &                   CAN_WCNT,STF_HF_SNOW_MELT,                     &
     &                   STF_SUB_SURF_ROFF,SMCL,STHF,STHU,TSOIL,        &
     &                   CAN_WCNT_GB,HF_SNOW_MELT,SMC,                  &
     &                   SNOW_MELT,                                     &
     &                   SUB_SURF_ROFF,SURF_ROFF,TOT_TFALL,             &
! Add inland basin outflow to subroutine call to hydrol
     &                   inlandout_atm,L_INLAND,                        &



! Additional variables for MOSES II
     &                   NTILES,TILE_PTS,TILE_INDEX,                    &
     &                   INFIL_TILE,                                    &
     &                   MELT_TILE,FRAC,                                &
     
! Additional variables required for large-scale hydrology:
     &    L_TOP,L_PDM,FEXP,GAMTOT,TI_MEAN,TI_SIG,CS_CH4,                &
     &    DUN_ROFF,DRAIN,FSAT,FWETL,QBASE,QBASE_ZW,ZW,                  &
     &    STHZW,A_FSAT,C_FSAT,A_FWET,C_FWET,FCH4_WETL,                  &
     &    dim_cs1,L_SOIL_SAT_DOWN,l_triffid,                            &

     & LTIMER                                                           &
     & )

      IF (lhook) CALL dr_hook('HYD_INTCTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE HYD_INTCTL
