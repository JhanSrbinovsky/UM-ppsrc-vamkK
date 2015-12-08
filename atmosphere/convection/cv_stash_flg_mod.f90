! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Set flags for use in convection code

MODULE cv_stash_flg_mod

! Description:
!   Module containing stash flags used in top level convection routine.
!
! Method:
!   Declares all flags and contains a subroutine used to set the flags.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
SAVE

!=====================================================================
! Logical flags used to control output/setting of various convection 
! diagnostics
!=====================================================================

LOGICAL ::           &
  L_apply_diag         ! true if sub-step for calculating diagnostics

LOGICAL ::           &
  flg_up_flx         & ! stash flag for updraught mass flux             
 ,flg_up_flx_half    & ! stash flag for updraught mass flux
 ,flg_dwn_flx        & ! stash flag for downdraught mass flux
 ,flg_entr_up        & ! stash flag for updraught entrainment
 ,flg_entr_dwn       & ! stash flag for downdraught entrainmnent
 ,flg_detr_up        & ! stash flag for updraught detrainment
 ,flg_detr_dwn       & ! stash flag for downdraught detrainmnent
 ,flg_conv_rain_3d   & ! stash flag for 3d conv rainfall
 ,flg_conv_snow_3d     ! stash flag for 3d conv snowfall

! CMT diagnostics

LOGICAL ::           &
  flg_uw_dp          & ! stash flag for deep x-comp stress
 ,flg_vw_dp          & ! stash flag for deep y-comp stress
 ,flg_uw_shall       & ! stash flag for shallow x stress
 ,flg_vw_shall       & ! stash flag for shallow y stress
 ,flg_uw_mid         & ! stash flag for mid-level x stress
 ,flg_vw_mid           ! stash flag for mid-level y stress
! Increment diagnostics:

LOGICAL ::           &
  l_qcl_incr_cinh    & ! liquid cloud condensate qCL (inhom)
 ,l_qcf_incr_cinh    & ! frozen cloud condensate qCF (inhom)
 ,l_cfl_incr_cinh    & ! liquid cloud amount cf_liquid (inhom)
 ,l_cff_incr_cinh    & ! frozen cloud amount cf_frozen (inhom)
 ,l_bcf_incr_cinh    & ! total (bulk) cloud amount bulk_cf (inhom)

 ,l_T_incr_conv      & ! temperature
 ,l_q_incr_conv      & ! humidity
 ,l_qcl_incr_conv    & ! liquid cloud condensate qCL
 ,l_qcf_incr_conv    & ! frozen cloud condensate qCF
 ,l_cfl_incr_conv    & ! liquid cloud amount cf_liquid
 ,l_cff_incr_conv    & ! frozen cloud amount cf_frozen
 ,l_bcf_incr_conv    & ! total (bulk) cloud amount bulk_cf
 ,l_T_conv_only      & ! temperature from convection only
 ,l_q_conv_only      & ! humidity from convection only
 ,l_u_incr_conv      & ! u wind
 ,l_v_incr_conv      & ! v wind
 ,l_up_incr_conv     & ! u wind p-grid
 ,l_vp_incr_conv       ! v wind p_grid

! 5A turbulence based schemes

LOGICAL ::           &
  flg_wqt_flux       & ! stash flag for w'qt'
 ,flg_wql_flux       & ! stash flag for w'ql'
 ,flg_wthetal_flux   & ! stash flag for w'thetal'
 ,flg_wthetav_flux   & ! stash flag for w'thetav'

 ,flg_deep_tops      & ! stash flag for deep tops

 ,flg_mf_deep        & ! stash flag for deep mass flux
 ,flg_mf_congest     & ! stash flag for congestus mass flux
 ,flg_mf_shall       & ! stash flag for shallow mass flux
 ,flg_mf_midlev      & ! stash flag for mid_level mass flux

 ,flg_dt_deep        & ! stash flag for deep dT
 ,flg_dt_congest     & ! stash flag for congestus dT
 ,flg_dt_shall       & ! stash flag for shallow dT
 ,flg_dt_midlev      & ! stash flag for mid_level dT

 ,flg_dq_deep        & ! stash flag for deep dq
 ,flg_dq_congest     & ! stash flag for congestus dq
 ,flg_dq_shall       & ! stash flag for shallow dq
 ,flg_dq_midlev      & ! stash flag for mid_level dq

 ,flg_du_deep        & ! stash flag for deep du
 ,flg_du_congest     & ! stash flag for congestus du
 ,flg_du_shall       & ! stash flag for shallow du
 ,flg_du_midlev      & ! stash flag for mid_level du

 ,flg_dv_deep        & ! stash flag for deep dv
 ,flg_dv_congest     & ! stash flag for congestus dv
 ,flg_dv_shall       & ! stash flag for shallow dv
 ,flg_dv_midlev        ! stash flag for mid_level dv


CONTAINS

!=========================================================================
! Full model version - sets flags according to user requirements
! i.e. stash settings held in array sf
!=========================================================================
SUBROUTINE set_convection_output_flags                                  &

           (nsects,nitems,                                              &
            L_calc_dxek, L_cosp, sf)

USE BL_OPTION_MOD, ONLY: &
  ISrfExCnvGust

IMPLICIT NONE
  
INTEGER, intent(in) ::  &
  nsects                & ! total number of sections
 ,nitems                  ! total numebr of items

LOGICAL, intent(in) ::  &
  l_calc_dxek             ! Switch for calculation of condensate increment

LOGICAL, INTENT(IN) ::  L_cosp ! Flag to request diagnostics needed by COSP

LOGICAL, intent(in) ::  &
  sf(0:nitems,0:nsects)   ! stash flags controlling output

! Start blopt8a

! Description:
!   Permissible settings for BL options.


      INTEGER, PARAMETER :: off = 0  ! Switch disabled
      INTEGER, PARAMETER :: on  = 1  ! Switch enabled

!     Options for non-gradient stress following
      INTEGER, PARAMETER :: BrownGrant97 = 1
      INTEGER, PARAMETER :: BrownGrant97_limited = 2
!       Brown and Grant (1997), version 2 including a limit on its size

!     Options for flux gradient formulation
      INTEGER, PARAMETER :: Locketal2000   = 0
!       Flux gradients as in Lock et al. (2000)
      INTEGER, PARAMETER :: HoltBov1993 = 1
!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
      INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)

!     Options for entrainment enhancement in Sc over Cu
      INTEGER, PARAMETER :: Buoyrev_feedback = 1

!     Options for form drag
      INTEGER, PARAMETER :: No_drag         = 0
      INTEGER, PARAMETER :: Effective_z0    = 1
      INTEGER, PARAMETER :: Explicit_stress = 2

!     Options for marine boundary layers
      INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
      INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory
      INTEGER, PARAMETER :: DynDiag_ZL = 1
      INTEGER, PARAMETER :: DynDiag_ZL_corrn = 2
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used as a dynamic criterion in the
!       diagnosis of BL types: version 2 includes changes to
!       cope with BL_LEVELS >> 3km
      INTEGER, PARAMETER :: DynDiag_ZL_CuOnly = 3
!       As 2 but only applied to points diagnosed with Cumulus 
!       and strictly for sea points (fland<0.01, cf 0.5)
      INTEGER, PARAMETER :: DynDiag_Ribased = 4
!       As 3 but also overrides Cumulus diagnosis if 
!          ZH(Ri) > ZLCL+zhloc_depth_fac*(ZHPAR-ZLCL)
!       Note that here Ri accounts for gradient adjustment by the 
!       non-local scheme.

!     Options for surface exchange
      INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
      INTEGER, PARAMETER :: Limit_ObukhovL = 3
!       Option under the COR_MO_ITER switch for imposing
!       lower limit on L in very stable conditions.
      INTEGER, PARAMETER :: Limit_expl_ustar = 2
!       Option under the COR_UST switch to limit the magnitude of the
!       explicitly calculated ustar
      INTEGER, PARAMETER :: IP_SrfExWithCnv = 1
!       Option to include deep convective gustiness in the surface
!       transfer

!     Options for convective boundary layers
      INTEGER, PARAMETER ::  UM_std     = 0
      INTEGER, PARAMETER ::  neut_cbl   = 1
      INTEGER, PARAMETER ::  LEM_conven = 2
      INTEGER, PARAMETER ::  LEM_std    = 3

!     Options for stable boundary layers
      INTEGER, PARAMETER ::  Long_tails           = 0
      INTEGER, PARAMETER ::  Sharpest             = 1
      INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
      INTEGER, PARAMETER ::  Mes_tails            = 3
      INTEGER, PARAMETER ::  Louis_tails          = 4
      INTEGER, PARAMETER ::  Depth_based          = 5
      INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
      INTEGER, PARAMETER ::  LEM_stability        = 7
      INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8

!     Options for Prandtl number (in local Ri scheme)
      INTEGER, PARAMETER ::  Constant_SBL = 0
      INTEGER, PARAMETER ::  LockMailhot2004 = 1

! End blopt8a

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!---------------------------------------------------------------------------

IF (lhook) CALL dr_hook('SET_CONVECTION_OUTPUT_FLAGS',zhook_in,zhook_handle)

flg_conv_rain_3d = (sf(227,5) .AND. l_apply_diag) .OR. L_cosp
flg_conv_snow_3d = (sf(228,5) .AND. l_apply_diag) .OR. L_cosp

flg_up_flx = sf(250,5) .and. l_apply_diag 

flg_up_flx_half = ( sf(249,5) .OR. sf(246,5) ) .and. l_apply_diag  
flg_dwn_flx  = ( sf(251,5) .OR. (ISrfExCnvGust == IP_SrfExWithCnv) ) &
                 .and. l_apply_diag 
flg_entr_up  = sf(252,5) .and. l_apply_diag 
flg_detr_up  = sf(253,5) .and. l_apply_diag 
flg_entr_dwn = sf(254,5) .and. l_apply_diag 
flg_detr_dwn = sf(255,5) .and. l_apply_diag 

flg_uw_dp    = sf(258,5) .and. l_apply_diag 
flg_vw_dp    = sf(259,5) .and. l_apply_diag 
flg_uw_shall = sf(260,5) .and. l_apply_diag 
flg_vw_shall = sf(261,5) .and. l_apply_diag 
flg_uw_mid   = sf(263,5) .and. l_apply_diag 
flg_vw_mid   = sf(264,5) .and. l_apply_diag 
  
! PC2 increments

l_qcl_incr_cinh = (sf(163,5) .OR. (sf(140,5) .OR. sf(141,5))) &
                   .and. l_apply_diag 
l_qcf_incr_cinh = (sf(164,5) .OR. (sf(142,5) .OR. sf(143,5))) &
                   .and. l_apply_diag 
l_bcf_incr_cinh = sf(172,5) .and. l_apply_diag 
l_cfl_incr_cinh = (sf(173,5) .OR. (sf(146,5) .OR. sf(147,5))) &
                   .and. l_apply_diag 
l_cff_incr_cinh = (sf(174,5) .OR. (sf(148,5) .OR. sf(149,5))) &
                   .and. l_apply_diag 
l_bcf_incr_conv = sf(192,5) .and. l_apply_diag 
l_cfl_incr_conv = (sf(193,5) .OR. (sf(156,5).OR.sf(157,5).OR.sf(158,5)))&
                   .and. l_apply_diag 
l_cff_incr_conv = sf(194,5) .and. l_apply_diag 

! Convection increments

l_T_incr_conv   = (( sf(181,5) .or. sf(187,5) .or. sf(161,5))           &
                     .or. L_calc_dxek) .and. l_apply_diag  
l_q_incr_conv   = (sf(182,5) .or. sf(188,5) .or. sf(162,5)) .and. l_apply_diag 
l_qcl_incr_conv = (sf(183,5) .OR. (sf(150,5).OR.sf(151,5).OR.sf(152,5)))&
                     .and. l_apply_diag 
l_qcf_incr_conv = sf(184,5) .and. l_apply_diag 
l_u_incr_conv   = sf(185,5) .and. l_apply_diag 
l_v_incr_conv   = sf(186,5) .and. l_apply_diag 
l_up_incr_conv  = sf(256,5) .and. l_apply_diag 
l_vp_incr_conv  = sf(257,5) .and. l_apply_diag 

l_T_conv_only   = sf(161,5) .and. l_apply_diag  
l_q_conv_only   = sf(162,5) .and. l_apply_diag 

! 5A convection code 

flg_wqt_flux     = (sf(290,5) .OR. sf(304,5) .OR. sf(306,5) .OR.         &
                 sf(429,5) .OR. sf(430,5) .OR. sf(431,5) .OR. sf(432,5)) &
                   .and. l_apply_diag  
flg_wql_flux     = sf(291,5) .and. l_apply_diag 
flg_wthetal_flux = (sf(292,5) .OR. sf(305,5) .OR. sf(307,5) .OR.         &
                 sf(425,5) .OR. sf(426,5) .OR. sf(427,5) .OR. sf(428,5)) &
                   .and. l_apply_diag 
flg_wthetav_flux = sf(293,5) .and. l_apply_diag 

flg_deep_tops    = sf(319,5) .and. l_apply_diag 

flg_mf_deep      = sf(320,5) .and. l_apply_diag 
flg_mf_congest   = sf(321,5) .and. l_apply_diag 
flg_mf_shall     = (sf(322,5) .OR. sf(417,5) .OR. sf(418,5) .OR.         &
                   sf(419,5) .OR. sf(420,5)) .and. l_apply_diag  
flg_mf_midlev    = sf(323,5) .and. l_apply_diag 
flg_dt_deep      = sf(324,5) .and. l_apply_diag 
flg_dt_congest   = sf(325,5) .and. l_apply_diag 
flg_dt_shall     = (sf(326,5) .OR. sf(409,5) .OR. sf(410,5) .OR.         &
                   sf(411,5) .OR. sf(412,5)) .and. l_apply_diag 
flg_dt_midlev    = sf(327,5) .and. l_apply_diag 
flg_dq_deep      = sf(328,5) .and. l_apply_diag 
flg_dq_congest   = sf(329,5) .and. l_apply_diag 
flg_dq_shall     = (sf(330,5) .OR. sf(413,5) .OR. sf(414,5) .OR.         &
                   sf(415,5) .OR. sf(416,5)) .and. l_apply_diag 
flg_dq_midlev    = sf(331,5) .and. l_apply_diag 
flg_du_deep      = sf(332,5) .and. l_apply_diag 
flg_du_congest   = sf(333,5) .and. l_apply_diag 
flg_du_shall     = sf(334,5) .and. l_apply_diag 
flg_du_midlev    = sf(335,5) .and. l_apply_diag 
flg_dv_deep      = sf(336,5) .and. l_apply_diag 
flg_dv_congest   = sf(337,5) .and. l_apply_diag 
flg_dv_shall     = sf(338,5) .and. l_apply_diag 
flg_dv_midlev    = sf(339,5) .and. l_apply_diag 


IF (lhook) CALL dr_hook('SET_CONVECTION_OUTPUT_FLAGS',zhook_out,zhook_handle)
RETURN
END  SUBROUTINE set_convection_output_flags

END MODULE cv_stash_flg_mod
