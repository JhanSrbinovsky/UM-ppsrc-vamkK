! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the effect of convection upon the large-scale atmosphere
!
! Subroutine Interface:
SUBROUTINE environ (k, npnts, np_full, ntra, sdet_on,                     &
                    l_tracer, l_mom_gk,l_calc_dxek, l_q_interact,         &
                    bwk, bwkp1, bterm, blowst,                            &
                    l_shallow, cumulus, l_mid,                            &
                    timestep,                                             &
                    thek, qek, qclek, qcfek, bcfek, cflek, cffek,         &
                    thekp1,qekp1,qclekp1,qcfekp1,bcfekp1,cflekp1,cffekp1, &
                    uek, vek, uekp1, vekp1, traek, traekp1,               &
                    thpk, qpk, qclpk, qcfpk, flxk, thrk, qrk,             &
                    thpkp1, qpkp1, qclpkp1, qcfpkp1, flxkp1,              &
                    upk, vpk, upkp1, vpkp1, trapk, trapkp1,               &
                    deltak, ekp14, exk, exkp1, delpk, delpkp1,            &
                    delp_uv_k, delp_uv_kp1,                               &
                    thpixs_v, qpixs_v, amdetk, t1_sd, q1_sd,              &
                    dthek, dqek, dqclek, dqcfek, dbcfek, dcflek, dcffek,  &
                    xpk, xpkp1,                                           & 
                    duek, dvek, eflux_u_ud, eflux_v_ud, dtraek,           &
                    dthekp1, dqekp1, dqclekp1, dqcfekp1,                  &
                    dbcfekp1, dcflekp1, dcffekp1,                         &
                    dqek_nonpc2, dthek_nonpc2,                            &
                    dqekp1_nonpc2, dthekp1_nonpc2,                        &
                    duekp1, dvekp1, dtraekp1,                             &
                    max_cfl_c )

USE water_constants_mod, ONLY: lc, lf, tm

USE atmos_constants_mod, ONLY: r, cp

USE cv_run_mod, ONLY:                                                   &
    bl_cnv_mix

USE cloud_inputs_mod, ONLY: i_pc2_conv_coupling,                        &
                            starticeTKelvin, alliceTdegC

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Description : Calculate the effect of convection upon the large-scale 
!                atmosphere.
!
!   Method    : See Unified Model documentation paper 27.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER,  INTENT(IN) :: &
  k                     & ! model layer number    
 ,npnts                 & ! Vector length
 ,np_full               & ! Full vector length
 ,ntra                  & ! Number of tracers
 ,sdet_on                 ! Flag smoothed forced detrainment

LOGICAL, INTENT(IN) ::  &
  l_tracer              & ! Switch for inclusion of tracers
 ,l_mom_gk              & ! Switch for inclusion of momentum transport
                          !  (Gregory-Kershaw scheme)
 ,l_calc_dxek           & ! Switch for calculating condensate inc.
 ,l_q_interact            ! Switch allows overwriting of parcel variables when
                          ! calculating condensate increments (will alter 
                          ! results).

LOGICAL, INTENT(IN) ::  &
  bwk(npnts)            & ! Mask for whether condensate is liquid in layer K
 ,bwkp1(npnts)          & ! and layer K+1
 ,bterm(npnts)          & ! Mask for parcels which terminate in layer k+1
 ,blowst(npnts)         & ! Mask for those points at which stability is low
                          !  enough for convection to occur
 ,l_shallow(npnts)      & ! Switches for type of convection likely to develop
 ,cumulus(npnts)        & ! Cumulus convection flag (deep and shallow)
 ,l_mid(npnts)            ! mid - type 

REAL, INTENT(IN) :: &
  timestep            ! Convection timestep (s) 
                      !         (= model timestep if conv. called once)

! Cloud environment values in layer k
REAL, INTENT(IN) :: &
  thek(npnts)       & ! Potential temperature of cloud environment in layer k(K)
 ,qek(npnts)        & ! Mixing ratio of cloud environment in layer k (kg/kg)
 ,qclek(npnts)      & ! Liquid condensate mixing ratio of cloud environment
                      ! in layer k (kg/kg)
 ,qcfek(npnts)      & ! Frozen condensate mixing ratio of cloud environment
                      ! in layer k (kg/kg)
 ,bcfek(npnts)      & ! Total cloud volume fraction of cloud environment in 
                      ! layer k 
 ,cflek(npnts)      & ! Liquid cloud volume fraction of cloud environment in 
                      ! layer k
 ,cffek(npnts)        ! Frozen cloud volume fraction of cloud environment in 
                      ! layer k

! Cloud environment values in layer k+1
REAL, INTENT(IN) :: &
  thekp1(npnts)     & ! Potential temperature of cloud environment in layer
                      !  k+1 (K)
 ,qekp1(npnts)      & ! Mixing ratio of cloud environment in layer k+1 (kg/kg)
 ,qclekp1(npnts)    & ! Liquid condensate mixing ratio of cloud environment
                      ! in layer k+1 (kg/kg)
 ,qcfekp1(npnts)    & ! Frozen condensate mixing ratio of cloud environment
                      ! in layer k+1 (kg/kg)
 ,bcfekp1(npnts)    & ! Total cloud volume fraction of cloud environment in 
                      ! layer k+1 
 ,cflekp1(npnts)    & ! Liquid cloud volume fraction of cloud environment in 
                      ! layer k+1
 ,cffekp1(npnts)      ! Frozen cloud volume fraction of cloud environment in 
                      ! layer k+1

! Used by CMT
REAL, INTENT(IN)  :: &
  uek(npnts)         & ! Environment U in layer k (m/s) (uv levels)
 ,vek(npnts)         & ! Environment V in layer k (m/s) (uv levels)
 ,uekp1(npnts)       & ! Environment U in layer k+1 (m/s)
 ,vekp1(npnts)         ! Environment v in layer k+1 (m/s)

REAL, INTENT(IN)  ::     &
  traek(np_full,ntra)    & ! Tracer of cloud environment in layer k (kg/kg)
 ,traekp1(np_full,ntra)    ! Tracer of cloud environment in layer k+1 (kg/kg)

! Parcel values in layer k
REAL, INTENT(IN)  :: &
  thpk(npnts)        & ! Parcel potential temperature in layer k (K)
 ,qpk(npnts)         & ! Parcel mixing ratio in layer k (kg/kg)
 ,qclpk(npnts)       & ! Parcel liquid condensate mixing ratio in 
                       ! layer k (kg/kg)
 ,qcfpk(npnts)       & ! Parcel frozen condensate mixing ratio in 
                       ! layer k (kg/kg)
 ,flxk(npnts)        & ! Parcel mass flux in layer k (Pa/s)
 ,thrk(npnts)        & ! Parcel detrainment potential temperature in layer k (K)
 ,qrk(npnts)           ! Parcel detrainment mixing ratio in layer k (kg/kg)

! Parcel values in layer k+1
REAL, INTENT(IN)  :: &
  thpkp1(npnts)      & ! Parcel potential temperature in layer k+1 (K)
 ,qpkp1(npnts)       & ! Parcel mixing ratio in layer k+1 (kg/kg)
 ,qclpkp1(npnts)     & ! Parcel liquid condensate mixing ratio in 
                       ! layer k+! (kg/kg)
 ,qcfpkp1(npnts)     & ! Parcel frozen condensate mixing ratio in 
                       ! layer k+1 (kg/kg)
 ,flxkp1(npnts)        ! Parcel mass flux in layer k+1 (Pa/s)

REAL, INTENT(IN)  :: &
  upk(npnts)         & ! Parcel U in layer k (m/s) (UV levels)
 ,vpk(npnts)         & ! Parcel V in layer k (m/s) (UV levels)
 ,upkp1(npnts)       & ! Parcel U in layer k+1 (m/s) (UV levels)
 ,vpkp1(npnts)         ! Parcel V in layer k+1 (m/s) (UV levels)

REAL, INTENT(IN)  ::     &
  trapk(np_full,ntra)    &  ! Parcel Tracer in layer k (kg/kg)
 ,trapkp1(np_full,ntra)     ! Parcel Tracer in layer k+1 (kg/kg)

REAL, INTENT(IN)  :: &
  deltak(npnts)      & ! Parcel forced detrainment rate in layer k multiplied
                       ! by appropriate layer thickness
 ,ekp14(npnts)       & ! Entrainment rate for level k+1/4 multiplied
                       ! by appropriate layer thickness
 ,exk(npnts)         & ! Exner ratio for mid-point of layer k
 ,exkp1(npnts)       & ! Exner ratio for mid-point of layer k+1
 ,delpk(npnts)       & ! Pressure difference across layer k (Pa)
 ,delpkp1(npnts)     & ! Pressure difference across layer k+1 (Pa)
 ,delp_uv_k(npnts)   & ! Pressure difference across UV layer k (Pa)
 ,delp_uv_kp1(npnts)   ! Pressure difference across UV layer k+1 (Pa)

REAL, INTENT(IN) ::   &
  thpixs_v(npnts)     & ! Parcel excess of theta
 ,qpixs_v(npnts)      & ! Parcel excess of q
 ,amdetk(npnts)       & ! Mixing detrainment at level k multiplied by 
                        ! appropriate layer thickness
 ,t1_sd(npnts)        & ! Standard deviation of turbulent fluctuations of 
                        ! layer 1 temperature (K).
 ,q1_sd(npnts)          ! Standard deviation of turbulent fluctuations of 
                        ! layer 1 humidity (kg/kg).


REAL, INTENT(INOUT) :: &
  dthek(npnts)         & ! IN Increment to model potential temperature in layer
                         !    k due to convection (may be non-zero due to a 
                         !    previous split final detrainment calculation)(K/s)
                         ! OUT Updated increment to model potential temperature 
                         !     in layer k due to convection (K/s)
 ,dqek(npnts)          & ! IN Increment to model mixing ratio in layer k due 
                         !    to convection (may be non-zero due to a previous
                         !    split final detrainment calculation)(kg/kg/s)
                         ! OUT Updated increment to model mixing ratio
                         !     in layer k due to convection (kg/kg/s)
 ,dqclek(npnts)        & ! IN  Increment to model liquid condensate mixing 
                         !     ratio in layer k due to convection (kg/kg/s)
                         !     (may be non-zero due to a previous split final 
                         !      detrainment calculation)
                         ! OUT Updated increment to model liquid condensate
                         !     mixing ratio in layer k  (kg/kg/s)
 ,dqcfek(npnts)        & ! IN  Increment to model frozen condensate mixing 
                         !     ratio in layer k due to convection (kg/kg/s)
                         !     (may be non-zero due to a previous split final 
                         !      detrainment calculation)
                         ! OUT Updated increment to model frozen condensate
                         !     mixing ratio in layer k  (kg/kg/s)
 ,dbcfek(npnts)        & ! IN  Increment to model total cloud volumne fraction
                         !     in layer k due to convection (may be non-zero
                         !     due to a previous split final detrainment 
                         !     calculation) (/s)
                         ! OUT Updated increment to model total cloud volumne
                         !     fraction in layer k due to convection (/s)
 ,dcflek(npnts)        & ! IN  Increment to model liquid cloud volumne fraction
                         !     in layer k due to convection (may be non-zero
                         !     due to a previous split final detrainment 
                         !     calculation) (/s)
                         ! OUT Updated increment to model liquid cloud volumne
                         !     fraction in layer k due to convection (/s)
 ,dcffek(npnts)        & ! IN  Increment to model frozen cloud volumne fraction
                         !     in layer k due to convection (may be non-zero
                         !     due to a previous split final detrainment 
                         !     calculation) (/s)
                         ! OUT Updated increment to model frozen cloud volumne
                         !     fraction in layer k due to convection (/s)
 ,xpk(npnts)          &  ! IN  Parcel cloud water in layer k (kg/kg) 
                         !     Should be redundant if L_q_interact .True.
                         ! OUT sum of qcl and qcf if L_q_interact .True.
 ,xpkp1(npnts)           ! IN Parcel cloud water in layer k+1 (kg/kg) 
                         !    Should be redundant if L_q_interact .True.
                         ! OUT sum of qclp1 and qcfp1 if L_q_interact .True.


REAL, INTENT(INOUT) :: &
  duek(npnts)          & ! IN Increment to model U due to convection m/s
                         ! OUT Updates increment to model U to convection
                         !     (at UV level K)
 ,dvek(npnts)          & ! IN Increment to model V due to convection m/s
                         ! OUT Updates increment to model V to convection
                         !     (at UV level K)
 ,eflux_u_ud(npnts)    & ! IN Eddy flux of momentum at bottom of a layer due to
                         !    Up Draught (UD)
                         ! OUT Eddy flux of momentum at current layer due to UD
 ,eflux_v_ud(npnts)      ! IN Eddy flux of momentum at bottom of a layer due to
                         !    Up Draught (UD)
                         ! OUT Eddy flux of momentum at current layer due to UD   

REAL, INTENT(INOUT) :: &
  dtraek(np_full,ntra)   ! IN Increment to model tracer in layer k due to
                         !    convection (may be non-zero due to a previous
                         !    split final detrainment calculation) (kg/kg/s)
                         ! OUT Updated Increment to model tracer in layer k
                         !      due to convection(kg/kg/s)

REAL, INTENT(INOUT) :: &
  dthekp1(npnts)       & ! Increment to model potential temperature in layer
                         !    k+1 due to convection (K/s)
 ,dqekp1(npnts)        & ! Increment to model mixing ratio in layer k+1 due 
                         !    to convection (kg/kg/s)
 ,dqclekp1(npnts)      & ! Increment to model liquid condensate mixing 
                         !     ratio in layer k+1 due to convection (kg/kg/s)
 ,dqcfekp1(npnts)      & ! Increment to model frozen condensate mixing 
                         !     ratio in layer k+1 due to convection (kg/kg/s)
 ,dbcfekp1(npnts)      & ! Increment to model total cloud volumne fraction
                         !     in layer k+1 due to convection (/s)
 ,dcflekp1(npnts)      & ! Increment to model liquid cloud volumne fraction
                         !     in layer k+1 due to convection  (/s)
 ,dcffekp1(npnts)        ! Increment to model frozen cloud volumne fraction
                         !     in layer k+1 due to convection  (/s)

REAL, INTENT(INOUT) ::  &
  dqekp1_nonpc2(npnts)  & ! Increment to qekp1 if PC2 was not in place
 ,dthekp1_nonpc2(npnts)   ! Increment to thekp1 if PC2 was not in place

REAL, INTENT(INOUT) ::  &
  dtraekp1(np_full,ntra)    ! Increment to model tracer in layer K+1 due to
                            !  convection (kg/kg)
!-----------------------------------------------------------------------
! VARIABLES THAT ARE OUTPUT
!-----------------------------------------------------------------------
REAL, INTENT(OUT) ::    &
  dqek_nonpc2(npnts)    & ! Increment to qek if PC2 was not in place
 ,dthek_nonpc2(npnts)     ! Increment to thek if PC2 was not in place

REAL, INTENT(OUT) ::    &
  duekp1(npnts)         & ! Increment to model U in uv layer k+1 due to
                          !  convection (m/s/s)
 ,dvekp1(npnts)           ! Increment to model V in uv layer k+1 due to
                          !  convection (m/s/s) 

REAL, INTENT(OUT) ::    &
  max_cfl_c(npnts)        ! OUT CFL RATIO

!-----------------------------------------------------------------------
! Variables that are defined locally
!-----------------------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER ::     &
  i,ktra         ! LOOP COUNTERS

LOGICAL ::     &
  l_calc_qcm     ! Switch for mixed phase cloud calc. option

REAL ::              & ! Smoothed forced detrainment variables
  fd_dthek(npnts)    & ! Forced detrainment P.temp inc at level k
 ,tmp_fd_dthek       & ! Forced detrainment P.temp inc across levels k and k+1
 ,fd_dthekp1(npnts)  & ! Forced detrainment P.temp inc at level k+1
 ,fd_dqek(npnts)     & ! Forced detrainment humidity inc in at level k
 ,tmp_fd_dqek        & ! Forced detrainment humidity inc across levels k and k+1
 ,fd_dqekp1(npnts)   & ! Forced detrainment humidity inc in at level k+1
 ,fd_dqclek(npnts)   & ! Forced detrainment qcl inc at level k
 ,tmp_fd_dqclek      & ! Forced detrainment qcl inc across levels k and k+1
 ,fd_dqclekp1(npnts) & ! Forced detrainment qcl inc at level k+1
 ,fd_dqcfek(npnts)   & ! Forced detrainment qcf inc at level k
 ,tmp_fd_dqcfek      & ! Forced detrainment qcf across levels k and k+1
 ,fd_dqcfekp1(npnts) & ! Forced detrainment qcf in at level k+1
 ,tmp_fd_dtraek      & ! Forced detrainment tracer across levels k and k+1
 ,tmp_fd_dthek_nonpc2& ! Forced detrainment non PC2 P.temp inc across levels k and k+1
 ,tmp_fd_dqek_nonpc2   ! Forced detrainment non PC2 humidity inc across levels k and k+1

REAL ::              & 
  el                 & ! Latent heat of condensation or (condensation plus
                       ! fusion) (J/kg)
 ,tempry             & ! Temporary array
 ,temp_com_sub       & ! Temporary array expands to vector.
 ,temp_for_det       & ! Temporary array scalar expands to vector.
 ,thpixs             & ! Parcel excess potential temperature (k)
 ,qpixs              & ! Parcel excess moisture (kg/kg)
 ,flx_u_kp0p5        & ! Flux of zonal momentum in cloud at top of current uv
                       ! layer
 ,flx_v_kp0p5        & ! Flux of meridional momentum in cloud at top of
                       !  current uv layer
 ,tempdcfxk          & ! Temporary workspace scalar expands to vector.
 ,dbydt_of_qclek     & ! Storage space for liquid condensate rate.
 ,dbydt_of_qcfek     & ! Storage space for frozen condensate rate.
 ,dbydt_of_cfek      & ! Storage space for total cloud increment.
 ,cfmek_i            & ! Storage space for mixed phase cloud on K.
 ,cfmekp1_i          & ! Storage space for mixed phase cloud on K+1.
 ,qcmek_i            & ! Storage space for mixed phase condensate on K.
 ,qcmekp1_i          & ! Storage space for mixed phase condensate K+1.
 ,dqcmek_i           & ! Storage space for mixed phase condensate inc K.
 ,deltaxl            & ! 1 If convection parcel is liquid, 0 else
 ,deltaxf            & ! 1 If convection parcel is frozen, 0 else
 ,calc_dxek          & ! 0 If calculating condensate increment, 1 else.
 ,tempdcfxks         & !
 ,tempdcfxkd         & !
 ,dbydt_of_cflek_sub & !
 ,dbydt_of_cffek_sub & !
 ,dbydt_of_cfmek_sub & !
 ,dbydt_of_qclek_det & !
 ,dbydt_of_qcfek_det & ! 
 ,tmp_dcfl           & !
 ,tmp_dcff           & !
 ,tmp_dcfl_d         & !
 ,tmp_dcfl_s         & !
 ,tmp_dcff_d         & !
 ,tmp_dcff_s         & !
 ,tmp_dcfm_d         & !
 ,tmp_dcfm_s         & !
 ,dettmp             & !
 ,fracice            & !
 ,fracliq            & !
 ,fracicep1          & !
 ,fracliqp1            ! 

! Parameters

REAL, PARAMETER ::        &
  tolerance = 1.0e-10,    & ! Threshold to avoid divide by zero errors
  ls0 = 5.0e-5,           & ! Minimum value for ls - l                         
  a_smth = 0.5,           & ! Parameter determining the weighting between
                            ! the increments at k and k-1 used when smoothed
                            ! forced detrainment is selected. 
  wcgrow = 5.e-4            ! Water content of new detrained cloud  

! Several options are available:
INTEGER, PARAMETER :: pc2_conv_original = 1        
! As described in Wilson et al (2008). Condensate increments from 
! detrainment and subsidence advection (combined) are turned into cloud 
! fraction increments using the inhomogeneous forcing method. Note that
! the abrupt change in phase of condensate in the convective plume is
! replicated in the condensate increments due to convection.
INTEGER, PARAMETER :: pc2_conv_maxincldqc = 2
! Unpublished. As above + Protects (in ni_conv_ctl) against 
! generation of inconsistently low cloud fraction implying 
! very high in-cloud condensate amounts.
INTEGER, PARAMETER :: pc2_conv_smooth_liqice = 3
! Unpublished. As above + The phase of the detrained condensate varies 
! smoothly according to the ambient temperature rather than having 
! an abrupt change.
INTEGER, PARAMETER :: pc2_conv_improved = 4      
! Unpublished. As above + Liq and ice cloud fraction increments due to 
! vertical advection by subsidence are calculated separately to 
! detrainment increments. The cloud fraction increments due to 
! detrainment are calculated by assuming the newly created convective 
! clouds have a specified in-cloud condensed water content.
INTEGER, PARAMETER :: pc2_conv_improved_mixed = 5
! As above but mixed-phase cloud fraction is also advected by 
! subsidence and the total cloud fraction updated accordingly.

!---------------------------------------------------------------------

IF (lhook) CALL dr_hook('ENVIRON',zhook_in,zhook_handle)

l_calc_qcm = .TRUE.

IF (l_q_interact) THEN
  calc_dxek = 0.0

!-----------------------------------------------------------------------
!     Xpk is probably undefined in this case and may never be used. At
!     some point it would be best to tidy up and remove it, but not if
!     the .false. option is still being retained. Meanwhile OVERWRITE!
!-----------------------------------------------------------------------

  DO i=1,npnts
!       NB can tell from bwk and bwkp1 which of qcl or qcf is zero.
    xpk(i)   = qclpk(i)   + qcfpk(i)
    xpkp1(i) = qclpkp1(i) + qcfpkp1(i)
  END DO

ELSE
  calc_dxek = 1.0
END IF

DO i=1,npnts

!-----------------------------------------------------------------------
!  Create a vector of latent heats
!-----------------------------------------------------------------------

  IF (bwk(i)) THEN
    el = lc
  ELSE
    el = lc + lf
  END IF

!----------------------------------------------------------------------
! Calculate parcel mass flux divided by the thickness of layer k.
! This value is used in several places in the subroutine.
!----------------------------------------------------------------------

  tempry = flxk(i)/delpk(i)

! L_calc_dxek_if1:
  IF (l_calc_dxek)  THEN

    temp_com_sub = (1.0 + ekp14(i)) * (1.0 - deltak(i)) *         &
                                      (1.0 - amdetk(i))
    temp_for_det = deltak(i) * (1.0 - amdetk(i))

! ----------------------------------------------------------------------
!         DELTAXL should be zero when QCLPK is zero and 1 otherwise,
!         DELTAXF should be zero when QCFPK is zero and 1 otherwise, as
!         currently, parcel has single phase only. It would be wise to
!         think through implications before allowing both = 1 at once!
! ----------------------------------------------------------------------

    IF (bwk(i)) THEN
      deltaxl = 1.
      deltaxf = 0.
    ELSE
      deltaxl = 0.
      deltaxf = 1.
    END IF

  END IF  ! L_calc_dxek_if1

  IF ( blowst(i) .AND. (bl_cnv_mix == 0) ) THEN

! ----------------------------------------------------------------------
!  At the lowest convective layer, the parcel mass flux is a flux from
!  the environment. i.e. the initial mass flux is entrained with excess
!  potential temperature and mixing ratio tpixs, qpixs.
!
!  UM Documentation Paper 27
!  Section (10), equation (39)
! ----------------------------------------------------------------------


    thpixs= thpixs_v(i)
    qpixs = qpixs_v(i)

    dthek(i) = dthek(i) - tempry*thpixs
    dqek(i) = dqek(i) - tempry*qpixs

  END IF

!-----------------------------------------------------------------------
! Calculate the increments due to forced convection and if smoothed
! forced detrainment is selected then split across k and k+1.
!-----------------------------------------------------------------------

! First calculate the forced detrainment for theta and q 
  IF (sdet_on == 1 .OR. sdet_on == 2) THEN
!         Potential temperature increment
    tmp_fd_dthek  = deltak(i) * (1.0-amdetk(i)) *                        &
                    (thrk(i)-thek(i)- (calc_dxek*(el/cp)*xpk(i)/exk(i)))

    fd_dthek(i)   = a_smth*tmp_fd_dthek

    fd_dthekp1(i) = (1.0-a_smth) * exk(i)/exkp1(i) *                     &
                          delpk(i)/delpkp1(i)*tmp_fd_dthek

!         Humidity increment
    tmp_fd_dqek   = deltak(i) * (1.0-amdetk(i)) *                        &
                                        (qrk(i)-qek(i)+(calc_dxek*xpk(i)))

    fd_dqek(i)    = a_smth*tmp_fd_dqek

    fd_dqekp1(i)  = (1.0-a_smth) * delpk(i)/delpkp1(i)*tmp_fd_dqek
  ELSE
!         Potential temperature increment
    fd_dthek(i)   = deltak(i) * (1.0-amdetk(i)) *                        &
                    (thrk(i)-thek(i)- (calc_dxek*(el/cp)*xpk(i)/exk(i)))

    fd_dthekp1(i) = 0.0

!         Q increment
    fd_dqek(i)    = deltak(i) * (1.0-amdetk(i)) *                        &
                                        (qrk(i)-qek(i)+(calc_dxek*xpk(i)))

    fd_dqekp1(i) = 0.0
  END IF

! Now calculate the forced detrainment for qcl and qcf
  IF (l_calc_dxek) THEN
    IF (sdet_on == 2) THEN
!         qcl increment
      tmp_fd_dqclek = temp_for_det * (qclpk(i) - qclek(i))
  
      fd_dqclek(i)  = a_smth*tmp_fd_dqclek
  
      fd_dqclekp1(i)= (1.0-a_smth) * delpk(i)/delpkp1(i)*tmp_fd_dqclek
  
!         qcf increment
      tmp_fd_dqcfek = temp_for_det * (qcfpk(i) - qcfek(i))
  
      fd_dqcfek(i)  = a_smth*tmp_fd_dqcfek
  
      fd_dqcfekp1(i)= (1.0-a_smth) * delpk(i)/delpkp1(i)*tmp_fd_dqcfek
  
  
    ELSE
!         qcl increment
      fd_dqclek(i)   = temp_for_det * (qclpk(i) - qclek(i))
  
      fd_dqclekp1(i) = 0.0
      
!         qcl increment
      fd_dqcfek(i)   = temp_for_det * (qcfpk(i) - qcfek(i))

      fd_dqcfekp1(i) = 0.0

    END IF
  END IF


! ---------------------------------------------------------------------
!  Effect of convection upon potential temperature of layer k
!
!  UM Documentation Paper 27
!  Section(10), equation (38A)
! (Modified to remove re-evaporation of condensate to vapour when the
!  effect of convection upon condensate is formally calculated.)
! --------------------------------------------------------------------
!

  IF (l_q_interact) THEN
    IF (sdet_on == 0 .OR. sdet_on == 1) THEN
!        This is the original old calculation

      dthek_nonpc2(i) = dthek(i) + tempry * (                     &
           (1+ekp14(i)) * (1.0-deltak(i)) *                       &
           (1-amdetk(i)) * (thekp1(i)-thek(i))                    &
                                             ! Compensating subsidence
         + deltak(i) * (1.0-amdetk(i)) *                          &                                                 
           (thrk(i)-thek(i)-((el/cp)*xpk(i)/exk(i)))              &
                                             ! Forced detrainment
         +                                                        &
           amdetk(i) * (thpk(i)-thek(i)-  ((el/cp)*xpk(i)/exk(i))) )  
                                             ! Mixing detrainment                                           
    ELSE IF (sdet_on == 2) THEN
!        This is the old calculation but with smoothed forced detrainment

      tmp_fd_dthek_nonpc2 = deltak(i) * (1.0-amdetk(i)) *         &                                                 
           (thrk(i)-thek(i)-((el/cp)*xpk(i)/exk(i)))              
       
      dthek_nonpc2(i) = dthek(i) + tempry * (                     &

           (1+ekp14(i)) * (1.0-deltak(i)) *                       &
           (1-amdetk(i)) * (thekp1(i)-thek(i))                    &
                                             ! Compensating subsidence
         + a_smth * tmp_fd_dthek_nonpc2                           &
                                             ! Forced detrainment
         +                                                        &
           amdetk(i) * (thpk(i)-thek(i)-  ((el/cp)*xpk(i)/exk(i))) )  
                                             ! Mixing detrainment                                           

      dthekp1_nonpc2(i) = dthekp1_nonpc2(i) +                     &
                          (1.0-a_smth) * exk(i)/exkp1(i) *        &
                          delpk(i)/delpkp1(i)*tmp_fd_dthek
                                             ! Forced det at k+1
    END IF
  END IF

!      This is the modified calculation suitable for PC2 and non-PC2
  dthek(i) = dthek(i) + tempry * (                                 &

           (1+ekp14(i)) * (1.0-deltak(i)) *                        &
           (1-amdetk(i)) * (thekp1(i)-thek(i))                     &
                                             ! Compensating subsidence
         +                                                         &
           fd_dthek(i)                                             &
                                             ! Forced detrainment
         +                                                         &
           amdetk(i) * (thpk(i)-thek(i)-                           &
                                (calc_dxek*(el/cp)*xpk(i)/exk(i))) )
                                             ! Mixing detrainment


  IF (.NOT. l_q_interact) THEN
!        copy dthek to dthek_nonpc2 because it is used later to
!        calculate dCAPEbydt
    dthek_nonpc2(i) = dthek(i)
  END IF

  IF (sdet_on == 1 .or. sdet_on == 2) THEN
    dthekp1(i) = dthekp1(i)+tempry*fd_dthekp1(i)
  END IF


! ---------------------------------------------------------------------
!  Effect of convection upon mixing ratio of layer k
!
!  UM Documentation Paper 27
!  Section (10), equation (38B)
! (Modified to remove re-evaporation of condensate to vapour when the
!  effect of convection upon condensate is formally calculated.)
! --------------------------------------------------------------------
!
  IF (l_q_interact) THEN
    IF (sdet_on == 0 .OR. sdet_on == 1) THEN
!        This is the original old calculation

      dqek_nonpc2(i) = dqek(i) + tempry * (                       &

           (1+ekp14(i)) * (1.0-deltak(i)) *                       &
           (1-amdetk(i)) * (qekp1(i)-qek(i))                      &
                                             ! Compensating subsidence
         +                                                        &
           deltak(i) * (1.0-amdetk(i)) *(qrk(i)-qek(i)+xpk(i))    & 
                                             ! Forced detrainment
         +                                                        &
           amdetk(i) * (qpk(i)-qek(i)+ xpk(i))   )
                                             ! Mixing detrainment
    ELSE IF (sdet_on == 2) THEN
!        This is the old calculation but with smoothed forced detrainment

      tmp_fd_dqek_nonpc2 = deltak(i) * (1.0-amdetk(i)) *          &
          (qrk(i)-qek(i)+xpk(i))
       
      dqek_nonpc2(i) = dqek(i) + tempry * (                       &

           (1+ekp14(i)) * (1.0-deltak(i)) *                       &
           (1-amdetk(i)) * (qekp1(i)-qek(i))                      &
                                             ! Compensating subsidence
         +                                                        &
           a_smth * tmp_fd_dqek_nonpc2                            & 
                                             ! Forced detrainment
         +                                                        &
           amdetk(i) * (qpk(i)-qek(i)+ xpk(i))   )
                                             ! Mixing detrainment
      dqekp1_nonpc2(i) = dqekp1_nonpc2(i) +                       &
         (1.0-a_smth)*delpk(i)/delpkp1(i)*tmp_fd_dqek
                                             ! Forced det at k+1
    END IF
  END IF



!      This is the modified calculation suitable for PC2 and non-PC2
  dqek(i) = dqek(i) + tempry * (                                  &

           (1+ekp14(i)) * (1.0-deltak(i)) *                       &
           (1-amdetk(i)) * (qekp1(i)-qek(i))                      &
                                             ! Compensating subsidence
         +                                                        &
           fd_dqek(i)                                             &
                                             ! Forced detrainment
         +                                                        &
           amdetk(i) * (qpk(i)-qek(i)+ (calc_dxek*xpk(i)))  )
                                             ! Mixing detrainment

  IF (.NOT. l_q_interact) THEN
!        copy dqek to dqek_nonpc2 because it is used later to
!        calculate dCAPEbydt
    dqek_nonpc2(i) = dqek(i)
  END IF


!----------------------------------------------------------------------
! Add the forced detrainment increment to K+1
!----------------------------------------------------------------------
  IF (sdet_on == 1 .or. sdet_on == 2) THEN
    dqekp1(i)  = dqekp1(i)+tempry*fd_dqekp1(i)
  END IF


! ----------------------------------------------------------------------
!  Optional effect of convection upon condensate mixing ratio of layer k
!  (Bushell et al, QJRMS 2003, 129, pp1435-1455)
! ----------------------------------------------------------------------

! L_calc_dxek_if2:
  IF (l_calc_dxek)  THEN

! ----------------------------------------------------------------------
!       Calculate the Q4 rates (liquid, frozen, mixed phase condensate).
! ----------------------------------------------------------------------

!     --Increment to Liquid Condensate Mixing Ratio--
    dbydt_of_qclek = tempry * (                                    &
       temp_com_sub * (qclekp1(i)-qclek(i))                        &
                                             ! Compensating subsidence

      + fd_dqclek(i)                                               &
                                             ! Forced detrainment
!   Forced detrainment of condensate preserves parcel condensate values.

      + amdetk(i) * (qclpk(i) - qclek(i))  )
                                             ! Mixing detrainment

!         --------------------------------------------------------------
!         Initial convecting layer - flux through cloud base is zero.
!         --------------------------------------------------------------

    IF (blowst(i)) THEN
      IF (deltak(i) == 0.0 .AND. amdetk(i) == 0.0) THEN
        ! No detrainment, so set rate to zero:
        dbydt_of_qclek = 0.0
      ELSE
        ! Remove subsidence contribution to rate:
        dbydt_of_qclek = dbydt_of_qclek                               &
                       - tempry * (1+ekp14(i)) * (qclekp1(i)-qclek(i))
      END IF
    END IF 

    fracice= (thek(i)*exk(i)-starticeTKelvin) /                         &
             (alliceTdegC-(starticeTKelvin-tm))
    fracice= MAX(0.0,fracice) 
    fracice = MIN(fracice,1.0) 
    fracliq = 1.0-fracice

    fracicep1= (thekp1(i)*exkp1(i)-starticeTKelvin) /                   &
               (alliceTdegC-(starticeTKelvin-tm))
    fracicep1= MAX(0.0,fracicep1) 
    fracicep1 = MIN(fracicep1,1.0) 
    fracliqp1 = 1.0-fracicep1 

    IF (i_pc2_conv_coupling == pc2_conv_original    .OR.              &
        i_pc2_conv_coupling == pc2_conv_maxincldqc) THEN
      ! Original method
      dqclek(i) = dqclek(i) + dbydt_of_qclek
    ELSE
      ! Convective plume is either liq or ice with an abrupt change. 
      ! Partition the condensate increments detrained from the 
      ! convective plume onto the large-scale more smoothly over
      ! a range of temperature. Adjust for latent heating.
      IF (dbydt_of_qclek > 0.0) THEN 
        dqcfek(i) = dqcfek(i) + fracice*dbydt_of_qclek 
        ! Extra freezing heats 
        dthek(i) =  dthek(i)  + fracice*dbydt_of_qclek*lf/(exk(i)*cp) 
        dqclek(i) = dqclek(i) + fracliq*dbydt_of_qclek
      ELSE
        dqclek(i) = dqclek(i) + dbydt_of_qclek
      END IF
    END IF

!----------------------------------------------------------------------
! Add the forced detrainment increment to K+1
!----------------------------------------------------------------------

    IF (sdet_on == 2) THEN
      dqclekp1(i) = dqclekp1(i)+tempry*fd_dqclekp1(i)
    END IF



!     --Increment rate to Frozen Condensate Mixing Ratio--
    dbydt_of_qcfek = tempry * (                                    &

        temp_com_sub * (qcfekp1(i)-qcfek(i))                       &
                                             ! Compensating subsidence

      + fd_dqcfek(i)                                               &
                                             ! Forced detrainment
!   Forced detrainment of condensate preserves parcel condensate values.

      + amdetk(i) * (qcfpk(i) - qcfek(i))  )
                                             ! Mixing detrainment


!         --------------------------------------------------------------
!         Initial convecting layer - flux through cloud base is zero.
!         --------------------------------------------------------------

    IF (blowst(i)) THEN
      IF (deltak(i) == 0.0 .AND. amdetk(i) == 0.0) THEN
        ! No detrainment, so set rate to zero:
        dbydt_of_qcfek = 0.0
      ELSE
        ! Remove subsidence contribution to rate:
        dbydt_of_qcfek = dbydt_of_qcfek                               &
                       - tempry * (1+ekp14(i)) * (qcfekp1(i)-qcfek(i))
      END IF
    END IF 

    IF (i_pc2_conv_coupling == pc2_conv_original    .OR.              &
        i_pc2_conv_coupling == pc2_conv_maxincldqc) THEN
      ! Original method
      dqcfek(i) = dqcfek(i) + dbydt_of_qcfek
    ELSE
      ! Convective plume is either liq or ice with an abrupt change. 
      ! Partition the condensate increments detrained from the 
      ! convective plume onto the large-scale more smoothly over
      ! a range of temperature. Adjust for latent heating.
      IF (dbydt_of_qcfek > 0.0) THEN 
        dqclek(i) = dqclek(i) + fracliq*dbydt_of_qcfek 
        ! Melting cools 
        dthek(i) = dthek(i)   - fracliq*dbydt_of_qcfek*lf/(exk(i)*cp) 
        dqcfek(i) = dqcfek(i) + fracice*dbydt_of_qcfek
      ELSE
        dqcfek(i) = dqcfek(i) + dbydt_of_qcfek
      END IF
    END IF

!----------------------------------------------------------------------
! Add the forced detrainment increment to K+1
!----------------------------------------------------------------------

    IF (sdet_on == 2) THEN
      dqcfekp1(i) = dqcfekp1(i)+tempry*fd_dqcfekp1(i)
    END IF


!     --Optional: Increment rate to Mixed Condensate Mixing Ratio-------

! L_calc_qcm_if1:
    IF (l_calc_qcm) THEN

      IF (cflek(i)  >   tolerance  .AND.  cffek(i)  >   tolerance) THEN
        cfmek_i   = cflek(i)   + cffek(i)   - bcfek(i)
        qcmek_i   = cfmek_i * ( (qclek(i) / cflek(i)) +            &
                                (qcfek(i) / cffek(i)) )
      ELSE
        cfmek_i   = 0.0
        qcmek_i   = 0.0
      END IF

      IF (cflekp1(i)  >   tolerance .AND. cffekp1(i)  >   tolerance) THEN
        cfmekp1_i = cflekp1(i) + cffekp1(i) - bcfekp1(i)
        qcmekp1_i = cfmekp1_i * ( (qclekp1(i) / cflekp1(i)) +      &
                                  (qcfekp1(i) / cffekp1(i)) )
      ELSE
        cfmekp1_i = 0.0
        qcmekp1_i = 0.0
      END IF

!           Parcel condensate for mixed phase is zero
      dqcmek_i = tempry * (                                        &

        temp_com_sub * (qcmekp1_i - qcmek_i)                       &
                                             ! Compensating subsidence

      - temp_for_det * (qcmek_i)                                   &
                                             ! Forced detrainment
     ! Forced detrainment of condensate preserves parcel condensate values.

      - amdetk(i) * (qcmek_i)  )
                                             ! Mixing detrainment

    END IF  ! L_calc_qcm_if1
  END IF  ! L_calc_dxek_if2

! ----------------------------------------------------------------------
!  Optional effect of convection upon cloud volume fractions of layer k
!  Block separate from condensate for clarity, but actually dependent on
!  calculations of Q4 condensation rates.
! ----------------------------------------------------------------------

! L_calc_dxek_if3:
  IF (l_calc_dxek)  THEN

! ----------------------------------------------------------------------
!       Calculate cloud volume rates (liquid, frozen and total).
! ----------------------------------------------------------------------

    IF (i_pc2_conv_coupling == pc2_conv_original       .OR.             &
        i_pc2_conv_coupling == pc2_conv_maxincldqc     .OR.             &
        i_pc2_conv_coupling == pc2_conv_smooth_liqice) THEN
      !=================================================================
      ! Original PC2 method, as used in Wilson et al (2008a,b). 
      !=================================================================

      ! --Increment to Liquid Cloud Volume Fraction----
      ! Using an explicit form of calculation with enforced limit
      tempdcfxk =      deltaxl * MAX( (qclpk(i)-qclek(i)  ) ,ls0)   &
                + (1.0-deltaxl)*    (          -qclek(i)  )
      !  deltaxl should be zero when qclpk is zero and 1 otherwise.

      IF (ABS(tempdcfxk)  >   tolerance) THEN
        tempdcfxk = dbydt_of_qclek*(deltaxl - cflek(i)) / tempdcfxk
        tempdcfxk = MAX(tempdcfxk, (-1. * cflek(i) / timestep))
        tempdcfxk = MIN(tempdcfxk, ((1. - cflek(i)) / timestep))
      ELSE
        tempdcfxk = 0.0                        ! Zero rate
      END IF

      dcflek(i) = dcflek(i) + tempdcfxk

      ! Initialize total cloud increment
      dbydt_of_cfek = tempdcfxk

      !  --Increment to Frozen Cloud Volume Fraction----
      ! Using an explicit form of calculation with enforced limit
      tempdcfxk =         deltaxf * MAX( (qcfpk(i)-qcfek(i)  ) ,ls0)    &
                   + (1.0-deltaxf)*    (          -qcfek(i)  )
      ! deltaxf should be zero when qcfpk is zero and 1 otherwise.

      IF (ABS(tempdcfxk)  >   tolerance) THEN
        tempdcfxk = dbydt_of_qcfek*(deltaxf - cffek(i)) / tempdcfxk
        tempdcfxk = MAX(tempdcfxk, (-1. * cffek(i) / timestep))
        tempdcfxk = MIN(tempdcfxk, ((1. - cffek(i)) / timestep))
      ELSE
        tempdcfxk = 0.0                        ! Zero rate
      END IF

      dcffek(i) = dcffek(i) + tempdcfxk
      dbydt_of_cfek = dbydt_of_cfek + tempdcfxk

!     --Increment to Mixed Phase Cloud Volume Fraction----
! L_calc_qcm_if2:
      IF (l_calc_qcm) THEN
!           Using an explicit form of calculation with enforced limit
        tempdcfxk = qcmek_i

        IF (ABS(tempdcfxk)   >    tolerance)  THEN
          tempdcfxk = dqcmek_i  * cfmek_i / tempdcfxk
          tempdcfxk = MAX(tempdcfxk, (-1. * cfmek_i / timestep))
          tempdcfxk = MIN(tempdcfxk, ((1. - cfmek_i) / timestep))
        ELSE
          tempdcfxk = 0.0                      ! Zero rate
        END IF

! ----------------------------------------------------------------------
!       Total values.
! ----------------------------------------------------------------------

!     --Increment to Total  Cloud Volume Fraction----
!       METHOD 1: Summation.
        tempdcfxk = dbydt_of_cfek - tempdcfxk
      ELSE  ! L_calc_qcm_if2

!     --Increment to Total  Cloud Volume Fraction----
!       METHOD 2: Direct specification.
!           Using an explicit form of calculation with enforced limit
        tempdcfxk = MAX((qclpk(i)+qcfpk(i)-qclek(i)-qcfek(i)),ls0)
        IF (ABS(tempdcfxk)   >    tolerance)  THEN
          tempdcfxk = (dbydt_of_qclek + dbydt_of_qcfek) *               &
                      (1.0 - bcfek(i)) / tempdcfxk
          tempdcfxk = MAX(tempdcfxk, (-1. * bcfek(i) / timestep))
          tempdcfxk = MIN(tempdcfxk, ((1. - bcfek(i)) / timestep))
        ELSE
          tempdcfxk = 0.0                      ! Zero rate
        END IF

      END IF  ! L_calc_qcm_if2

    ELSE ! (i_pc2_conv_coupling == pc2_conv_improved .OR.
         !  i_pc2_conv_coupling == pc2_conv_improved_mixed )

      ! ================================================================
      ! Improved methods
      ! ================================================================
      !
      ! The change in liquid & ice due to detrainment and subsidence are 
      ! considered separately.
      !
      ! IF i_pc2_conv_coupling == pc2_conv_improved_mixed
      ! changes to mixed-phase cloud fraction due to detrainment and 
      ! subsidence are also calculated separately for use in updating 
      ! the total (bulk) cloud fraction.
      ! ================================================================

      ! Liquid cloud fraction change.
       
      ! First consider detrainment. (Combined liq & ice plume condensate 
      ! as plume is single phase but partitioned using fracliq).
      dbydt_of_qclek_det = tempry * ( temp_for_det + amdetk(i) ) *      &
        fracliq * (qclpk(i) - qclek(i) + qcfpk(i) - qcfek(i))

      IF (dbydt_of_qclek_det > 1.0e-12) THEN
        ! Create some cloud fraction, such that new cloud has
        ! specified in-cloud condenstae amount.
        tempdcfxkd = dbydt_of_qclek_det * (1.0 - cflek(i)) / wcgrow 
      ELSE IF (dbydt_of_qclek_det < -1.0e-12) THEN
        ! Reduce cloud fraction such that in-cloud condensate amount 
        ! stays the same. Unless in-cloud value is v small, in which
        ! reduce cloud fraction so that sensible value is reached. 
        tempdcfxkd = dbydt_of_qclek_det * cflek(i) /                    & 
                       MAX( wcgrow, qclek(i)/MAX(1.e-6,cflek(i)) ) 
      ELSE 
        ! Change in condensate is small. Cloud fraction stays the same.
        tempdcfxkd = 0.0 
      ENDIF 
      ! Update change in liquid cloud fraction with detrainment incr.
      dcflek(i) = dcflek(i) + tempdcfxkd 
      ! Save for later use in mixed-phase calculation
      tmp_dcfl_d= tempdcfxkd 

      ! Subsidence of liquid cloud.
      dbydt_of_cflek_sub = tempry * temp_com_sub * (cflekp1(i)-cflek(i)) 
         
      IF (blowst(i)) dbydt_of_cflek_sub = dbydt_of_cflek_sub -          & 
                       tempry * (1+ekp14(i)) * (cflekp1(i)-cflek(i))

      ! Change in cloud fraction due to sub (considering liq for now)  
      tempdcfxks    = fracliq * dbydt_of_cflek_sub
      ! Detrainment + subsidence increments for liquid cloud fraction
      tempdcfxk     = tempdcfxks + tempdcfxkd
      ! Save for later use in mixed-phase calculation
      tmp_dcfl      = tempdcfxkd 
      ! Update change in liquid cloud fraction with subsidence incr.
      ! (so this now includes det + sub).
      dcflek(i)     = dcflek(i)  + tempdcfxks
      ! Total cloud fraction increment set to sum of detrainment + 
      ! subsidence increments for liquid cloud fraction. 
      dbydt_of_cfek = tempdcfxk 

      ! Ice cloud fraction change.

      ! First consider detrainment. (Combined liq & ice plume condensate 
      ! as plume is single phase but partitioned using fracice).
      dbydt_of_qcfek_det = tempry * ( temp_for_det + amdetk(i) ) *      &
        fracice * (qclpk(i) - qclek(i) + qcfpk(i) - qcfek(i))  

      IF (dbydt_of_qcfek_det > 1.0e-12) THEN 
        tempdcfxkd =  dbydt_of_qcfek_det*(1. - cffek(i))/wcgrow 
      ELSE IF (dbydt_of_qcfek_det < -1.0e-12) THEN 
        tempdcfxkd = dbydt_of_qcfek_det * cffek(i) /                    &
                       MAX( wcgrow, qcfek(i)/MAX(1.e-6,cffek(i)) ) 
      ELSE 
        tempdcfxkd = 0.0 
      ENDIF 
      ! Update change in ice cloud fraction with detrainment incr.
      dcffek(i) = dcffek(i) + tempdcfxkd 
      ! Save for later use in mixed-phase calculation
      tmp_dcff_d= tempdcfxkd 

      ! Subsidence of ice cloud.
      dbydt_of_cffek_sub = tempry * temp_com_sub * (cffekp1(i)-cffek(i)) 

      IF (blowst(i)) dbydt_of_cffek_sub = dbydt_of_cffek_sub -          & 
                       tempry * (1+ekp14(i)) * (cffekp1(i)-cffek(i)) 

      ! Change in cloud fraction due to sub (now considering ice) 
      tempdcfxks    = fracice * dbydt_of_cffek_sub 
      ! Detrainment + subsidence increments for ice cloud fraction
      tempdcfxk     = tempdcfxks + tempdcfxkd 
      ! Save for later use in mixed-phase calculation
      tmp_dcff      = tempdcfxkd 
      ! Update change in ice cloud fraction with subsidence incr.
      ! (so this now includes det + sub).
      dcffek(i)     = dcffek(i) + tempdcfxks 
      ! Update total cloud fraction increment. Now holds sum of 
      ! detrainment & subsidence increments for both liquid and ice 
      ! cloud fraction. 
      dbydt_of_cfek = dbydt_of_cfek + tempdcfxk

      ! Mixed-phase cloud fraction change
      IF (i_pc2_conv_coupling == pc2_conv_improved) THEN
        ! Original code

        IF (l_calc_qcm) THEN
          ! Method 1: Summation.

          ! Using an explicit form of calculation with enforced limit
          tempdcfxk = qcmek_i

          IF (ABS(tempdcfxk)   >    tolerance)  THEN
            tempdcfxk = dqcmek_i  * cfmek_i / tempdcfxk
            tempdcfxk = MAX(tempdcfxk, (-1.0 * cfmek_i / timestep))
            tempdcfxk = MIN(tempdcfxk, ((1.0 - cfmek_i) / timestep))
          ELSE
            tempdcfxk = 0.0 ! Zero rate
          END IF

          ! Increment to be applied to Total Cloud Fraction
          tempdcfxk = dbydt_of_cfek - tempdcfxk

        ELSE  ! l_calc_qcm
          ! Method 2: Direct specification.

          ! Using an explicit form of calculation with enforced limit
          tempdcfxk = MAX((qclpk(i)+qcfpk(i)-qclek(i)-qcfek(i)),ls0)

          IF (ABS(tempdcfxk) > tolerance)  THEN
            tempdcfxk = (dbydt_of_qclek + dbydt_of_qcfek) *             &
                        (1.0 - bcfek(i)) / tempdcfxk
            tempdcfxk = MAX(tempdcfxk, (-1.0 * bcfek(i) / timestep))
            tempdcfxk = MIN(tempdcfxk, ((1.0 - bcfek(i)) / timestep))
          ELSE
            tempdcfxk = 0.0 ! Zero rate
          END IF

        END IF  ! l_calc_qcm

      ELSE ! (i_pc2_conv_coupling == pc2_conv_improved_mixed)
        ! Seperating out detrainment and subsidence.

        IF (l_calc_qcm) THEN
          ! Calculate mixed-phase cloud fraction increments assuming 
          ! convective plume is mixed-phase between 
          ! starticeTKelvin and alliceTdegC.
          IF (fracliq > 0.0 .AND. fracliq < 1.0) THEN
            ! In temperature range where plume is mixed-phase.
            ! Detrained CFL and CFF from convection are assumed to be
            ! completely overlapped.
            tmp_dcfm_d = MIN(tmp_dcfl_d,tmp_dcff_d)
          ELSE 
            ! Plume is all liq or all ice, so detrainment will not change
            ! mixed-phase cloud fraction.
            tmp_dcfm_d = 0.0
          END IF

          ! Subsidence of mixed-phase cloud.
          dbydt_of_cfmek_sub = tempry*temp_com_sub * (cfmekp1_i-cfmek_i)

          IF (blowst(i)) dbydt_of_cfmek_sub = dbydt_of_cfmek_sub -      & 
                           tempry * (1+ekp14(i)) * (cfmekp1_i-cfmek_i) 

          ! Change in cloud fraction due to sub (considering mixed-phase)  
          tmp_dcfm_s = dbydt_of_cfmek_sub

          ! Change in mixed-phase frac is sum of contribs from det & sub.
          tempdcfxk = tmp_dcfm_d + tmp_dcfm_s

          ! Change in bulk (total) cloud fraction is sum of changes from
          ! liquid and ice MINUS the change from mixed-phase.
          tempdcfxk = tmp_dcff + tmp_dcff - tempdcfxk

        ELSE
          ! Do not allow convection to generate mixed-phase cloud frac
          ! or consider the effects of subs on mixed-phase fraction.
          tempdcfxk = 0.0
        END IF ! l_calc_qcm

      END IF

    END IF !i_pc2_conv_coupling

    dbcfek(i) = dbcfek(i) + tempdcfxk
    
  END IF  ! L_calc_dxek_if3

! Calculate courant number using mass flux at half level

  max_cfl_c(i)=tempry*(1+ekp14(i))*(1.0-deltak(i))*(1.0-amdetk(i))

! ----------------------------------------------------------------------
!  Terminal detrainment and subsidence in terminal layer
!
!  UM Documentation Paper 27
!  Section (10), equation (40)
! --------------------------------------------------------------------

  IF ( bterm(i) ) THEN
    IF (bwkp1(i)) THEN
      el = lc
    ELSE
      el = lc + lf
    END IF
    tempry = flxkp1(i)/delpkp1(i)
    IF (l_q_interact) THEN
      ! This is the old calculation
      dthekp1_nonpc2(i) = dthekp1(i)                              &
                            + tempry*((thpkp1(i)-thekp1(i))       &
                            - el*xpkp1(i)/(exkp1(i)*cp))
      dqekp1_nonpc2(i)  = dqekp1(i) + tempry*(qpkp1(i)-qekp1(i) + xpkp1(i))
    END IF
    dthekp1(i) = dthekp1(i) + tempry*((thpkp1(i)-thekp1(i))       &
                    - calc_dxek * el*xpkp1(i)/(exkp1(i)*cp))
    dqekp1(i)  = dqekp1(i) + tempry*(qpkp1(i)-qekp1(i)            &
                                     + calc_dxek * xpkp1(i))
    IF (.NOT. l_q_interact) THEN
      ! For safety, set dthekp1_nonpc2 and dqekp1_nonpc2, even
      ! though they aren't used.
      dthekp1_nonpc2(i) = dthekp1(i)
      dqekp1_nonpc2(i)  = dqekp1(i)
    END IF

    IF (l_calc_dxek)  THEN

      IF (i_pc2_conv_coupling == pc2_conv_original) THEN
      !
      ! Original PC2 method
      !
      dqclekp1(i)  = dqclekp1(i) + tempry*(qclpkp1(i)-qclekp1(i))

      dqcfekp1(i)  = dqcfekp1(i) + tempry*(qcfpkp1(i)-qcfekp1(i))

! Terminal detrainment of liquid

      IF ( (qclpkp1(i)-qclekp1(i))  >   ls0 ) THEN
        ! Go ahead with terminal detraiment
        dcflekp1(i)= dcflekp1(i) + tempry *(deltaxl - cflekp1(i))
      ELSE
        ! Scale down terminal detrainment of liquid cloud fraction
        dcflekp1(i)= dcflekp1(i) + tempry *(deltaxl - cflekp1(i)) &
                   * MAX( (qclpkp1(i)-qclekp1(i)),0.0)/ls0
      END IF

! Terminal detrainment of ice

      IF ( (qcfpkp1(i)-qcfekp1(i))  >   ls0 ) THEN
        ! Go ahead with terminal detraiment
        dcffekp1(i)= dcffekp1(i) + tempry *(deltaxf - cffekp1(i))
      ELSE
        ! Scale down terminal detrainment of frozen cloud fraction
        dcffekp1(i)= dcffekp1(i) + tempry *(deltaxf - cffekp1(i)) &
                   * MAX( (qcfpkp1(i)-qcfekp1(i)),0.0)/ls0
      END IF

! Terminal detrainment of bulk cloud fraction

      IF ( (qcfpkp1(i)-qcfekp1(i)+qclpkp1(i)-qclekp1(i))          &
             >   ls0 ) THEN
        ! Go ahead with terminal detraiment
        dbcfekp1(i)= dbcfekp1(i) + tempry * ( 1.0   - bcfekp1(i))
      ELSE
        ! Scale down terminal detrainment of bulk cloud fraction
        dbcfekp1(i)= dbcfekp1(i) + tempry *( 1.0 - bcfekp1(i))    &
                   * MAX( (qcfpkp1(i)-qcfekp1(i)                  &
                          +qclpkp1(i)-qclekp1(i)),0.0)/ls0
      END IF

      ELSE  ! i_pc2_conv_coupling /= pc2_conv_original
        !
        ! Alternative calculation PC2 cloud fraction changes 
        !
        DetTmp = tempry*(qclpkp1(i)-qclekp1(i)) 
        IF (DetTmp > 0.0) THEN 
          dqcfekp1(i) = dqcfekp1(i) + DetTmp*fracicep1 
          dthekp1(i)  = dthekp1(i)  + DetTmp*fracicep1*lf/(exkp1(i)*cp) 
          DetTmp = DetTmp*fracliqp1 
        END IF 
        dqclekp1(i)  = dqclekp1(I) + DetTmp 
 
        DetTmp = tempry*(qcfpkp1(i)-qcfekp1(i))  
        IF (DetTmp > 0.0) THEN 
          dqclekp1(i) = dqclekp1(i) + DetTmp*fracliqp1 
          dthekp1(i)  = dthekp1(i)  - DetTmp*fracliqp1*lf/(exkp1(i)*cp) 
          DetTmp = DetTmp*fracicep1 
        END IF 
        dqcfekp1(i)  = dqcfekp1(i) + DetTmp 

        DetTmp = tempry*(qclpkp1(i)-qclekp1(i)+qcfpkp1(i)-qcfekp1(i)) 

        IF (fracliqp1*DetTmp > 0.0) THEN 
          dcflekp1(i)= dcflekp1(i) +(1.0 - cflekp1(i))* fracliqp1 *     &
            DetTmp / wcgrow 
        ELSE 
          dcflekp1(i)= dcflekp1(i) + cflekp1(i) * fracliqp1 *           &
            DetTmp / MAX(wcgrow,qclekp1(i)) 
        END IF 

        IF (fracicep1*DetTmp > 0.0) THEN 
          dcffekp1(i)= dcffekp1(i) +(1.0 - cffekp1(i))* fracicep1 *     &
            DetTmp / wcgrow 
        ELSE 
          dcffekp1(i)= dcffekp1(i) + cffekp1(i) * fracicep1 *           &
            DetTmp / MAX(wcgrow,qcfekp1(i)) 
        END IF 

        IF (DetTmp.GT.0.0) THEN 
          dbcfekp1(i)= dbcfekp1(i) +(1.0 - bcfekp1(i)) * DetTmp / wcgrow 
        ELSE 
          dbcfekp1(i)= dbcfekp1(i) + bcfekp1(i) *                       &
            DetTmp /  MAX(wcgrow,(qcfekp1(i)+qclekp1(i))) 
        END IF 

      END IF ! i_pc2_conv_coupling

    END IF   ! endif l_calc_dxek

  END IF    ! endif bterm(i)

END DO

! ---------------------------------------------------------------------
!  Calculate effect of convection upon momentum of layer k and 
!  do terminal detrainment of momentum.
!
!  Rate of change of wind field by convection is estimated using a
!  divergence of vertical eddy momentum flux across the layer.
! --------------------------------------------------------------------
!
! All convective momentum transport calculations for the cumulus convection
! (deep and shallow) done when convection terminates.


IF(l_mom_gk) THEN      ! Gregory-Kershaw CMT

  DO i=1,npnts
!----------------------------------------------------------------------
! Estimate eddy flux at top of current UV layer due to convection
!----------------------------------------------------------------------

    flx_u_kp0p5 = flxk(i) * (1.0-amdetk(i)) * (1.0-deltak(i)) *      &
                            (1.0+ekp14(i)) * (upk(i)-uekp1(i))
    flx_v_kp0p5 = flxk(i) * (1.0-amdetk(i)) * (1.0-deltak(i)) *      &
                            (1.0+ekp14(i)) * (vpk(i)-vekp1(i))

    IF (blowst(i)) THEN
!----------------------------------------------------------------------
! Initial convecting layer - no flux at base of layer
!----------------------------------------------------------------------
      duek(i) = duek(i) - flx_u_kp0p5 / delp_uv_k(i)
      dvek(i) = dvek(i) - flx_v_kp0p5 / delp_uv_k(i)
!----------------------------------------------------------------------
! Store eddy flux at top of current UV layer ready for calculation
! of next layer.
!----------------------------------------------------------------------
      eflux_u_ud(i) = flx_u_kp0p5
      eflux_v_ud(i) = flx_v_kp0p5

    ELSE
!----------------------------------------------------------------------
! Convecting layer - take eddy flux divergence across the layer
!----------------------------------------------------------------------
      duek(i) = duek(i) - ( (flx_u_kp0p5 - eflux_u_ud(i)) /delp_uv_k(i) )  

      dvek(i) = dvek(i) - ( (flx_v_kp0p5 - eflux_v_ud(i)) /delp_uv_k(i) )

!----------------------------------------------------------------------
! Store eddy flux at top of curent UV layer ready for calculation of 
! next layer
!----------------------------------------------------------------------
      eflux_u_ud(i) = flx_u_kp0p5
      eflux_v_ud(i) = flx_v_kp0p5

    END IF

    IF(bterm(i))THEN
!----------------------------------------------------------------------
! Convection terminates - calculate increment due to convection in top 
! layer - no flux out of top layer.
!----------------------------------------------------------------------
      duekp1(i)  = eflux_u_ud(i) / delp_uv_kp1(i)
      dvekp1(i)  = eflux_v_ud(i) / delp_uv_kp1(i)
!----------------------------------------------------------------------
! Zero eddy flux out of top layer.
!----------------------------------------------------------------------
      eflux_u_ud(i) = 0.0
      eflux_v_ud(i) = 0.0

    END IF

  END DO

END IF   ! l_mom_gk

!----------------------------------------------------------------------
!  Effect of convection on tracer content of layer k.
!  (looping over number of tracer variables)
!  and do terminal detrainment of tracer.
!----------------------------------------------------------------------

IF(l_tracer)THEN

  IF (sdet_on == 0 .OR. sdet_on == 1) THEN
    DO ktra = 1,ntra
      DO i = 1,npnts

          tempry = flxk(i)/delpk(i)
          dtraek(i,ktra) = dtraek(i,ktra) + tempry * (                       &

           (1+ekp14(i)) * (1.0-deltak(i)) *                                  &
           (1-amdetk(i)) * (traekp1(i,ktra)-traek(i,ktra))                   &
                                                 ! Compensating subsidence

          + deltak(i) * (1.0-amdetk(i)) *(trapk(i,ktra)-traek(i,ktra))       &
                                                 ! Forced detrainment

          + amdetk(i) * (trapk(i,ktra)-traek(i,ktra)) )
                                                 ! Mixing detrainment

          IF(bterm(i))THEN
            tempry = flxkp1(i)/delpkp1(i)
            dtraekp1(i,ktra) = dtraekp1(i,ktra) +tempry*                     &
                                             (trapkp1(i,ktra)-traekp1(i,ktra))
          END IF

      END DO  ! i
    END DO    ! ktra
  ELSE IF (sdet_on == 2) THEN
    DO ktra = 1,ntra
      DO i = 1,npnts

        tempry = flxk(i)/delpk(i)
        
        tmp_fd_dtraek = deltak(i) * (1.0-amdetk(i)) *(trapk(i,ktra)-traek(i,ktra))
                                                 ! Forced detrainment
        
        dtraek(i,ktra) = dtraek(i,ktra) + tempry * (                         &

          (1+ekp14(i)) * (1.0-deltak(i)) *                                   &
          (1-amdetk(i)) * (traekp1(i,ktra)-traek(i,ktra))                    &
                                                 ! Compensating subsidence

        + a_smth * tmp_fd_dtraek                                              &
                                                 ! Forced detrainment

        + amdetk(i) * (trapk(i,ktra)-traek(i,ktra)) )
                                                 ! Mixing detrainment

        dtraekp1(i,ktra) = dtraekp1(i,ktra) +                                &
          tempry *(1.0-a_smth)*delpk(i)/delpkp1(i)*tmp_fd_dtraek
                                                 ! Forced detrainment at k+1
        IF(bterm(i))THEN
          tempry = flxkp1(i)/delpkp1(i)
          dtraekp1(i,ktra) = dtraekp1(i,ktra) +tempry*                     &
                                           (trapkp1(i,ktra)-traekp1(i,ktra))
        END IF

      END DO  ! i
    END DO    ! ktra

  END IF    ! sdet_on
END IF      ! l_tracer

IF (lhook) CALL dr_hook('ENVIRON',zhook_out,zhook_handle)
RETURN
END SUBROUTINE environ
