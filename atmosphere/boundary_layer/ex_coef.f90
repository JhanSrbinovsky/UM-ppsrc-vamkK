! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE EX_COEF------------------------------------------------
!
!  Purpose: To calculate exchange coefficients for boundary layer
!           subroutine KMKH.
!
!  Programming standard: Unified Model Documentation Paper No 3
!
!  Documentation: UMDP No.24
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer!
!---------------------------------------------------------------------
SUBROUTINE ex_coef (                                                    &
! IN levels/logicals
 bl_levels,k_log_layr,lq_mix_bl,l_subfilter_vert,l_subfilter_horiz,     &
 l_subfilter_blend,rneutml_sq,delta_smag,nSCMDpkgs,L_SCMDiags, BL_diag, &
! IN fields
 sigma_h,flandg,dbdz,dvdzm,ri,rho,z_uv,z_tq,z0m,h_blend,ntpar,ntml_nl,  &
 ntdsc,nbdsc,u_p,v_p,v_s,fb_surf,qw,tl,                                 &
! IN/OUT fields
 cumulus,weight_1dbl,                                                   &
! OUT fields
 lambda_min,zh_local,ntml_local,elm,elh,elh_rho,rhokm,rhokh,fm_3d,fh_3d &
 )

  USE atmos_constants_mod, ONLY: vkman, cp, c_virtual
  USE bl_option_mod, ONLY :  WeightLouisToLong, Variable_RiC, cbl_op,   &
     sg_orog_mixing, RiCrit_sharp, pr_max, l_lambdam2, l_full_lambdas,  &
     local_fa,Prandtl,ishear_bl,L_SBLco,Muw_SBL,Mwt_SBL,sbl_op,         &
     LockMailhot2004, depth_based, LEM_stability, LEM_std, LEM_conven,  & 
     off, on, sharpest, sharp_sea_long_land, sharp_sea_mes_land,        &
     louis_tails, sharp_sea_louis_land, long_tails, mes_tails,          &
     neut_cbl, equilibrium_sbl, lambda_min_nml, lambda_max_nml,         &
     lambda_fac, beta_bl, beta_fa,                                      &
     to_sharp_across_1km, ntml_level_corrn, free_trop_layers
  USE stochastic_physics_run_mod, ONLY:                                 &
     l_rp2, par_mezcla, g0_rp,ricrit_rp,lambda_min_rp   
  USE earth_constants_mod, ONLY: g
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE bl_diags_mod, ONLY: strnewbldiag
  USE atm_fields_bounds_mod, ONLY: pdims,tdims,qdims,pdims_s

  IMPLICIT NONE

! IN arguments
  LOGICAL, INTENT(IN) ::                                                &
   lq_mix_bl,                                                           &
                   ! IN True if mixing ratios used
   l_subfilter_blend,                                                   &
                   ! IN blending BL and Smag coefficients
   l_subfilter_vert,                                                    &
                   ! IN subgrid turbulence scheme in vertical
   l_subfilter_horiz
                   ! IN subgrid turbulence scheme in horizontal

  INTEGER, INTENT(IN) ::                                                &
   bl_levels,                                                           &
                   ! IN maximum number of boundary layer levels
   k_log_layr
                   ! IN num of levs requiring log-profile correction

  INTEGER, INTENT(IN) ::                                                &
   ntml_nl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                   ! IN Number of model layers in the turbulently
                   !    mixed layer as determined from the non-local
                   !    scheme.
   ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                   ! IN Top level of any decoupled Sc
   nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                   ! IN Bottom level of any decoupled Sc layer.
   ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! IN Top level of parcel ascent

  REAL, INTENT(IN) ::                                                   &
   sigma_h(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                   ! IN Standard deviation of subgrid
                   !    orography (m)
   rho(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                   ! IN density on theta levels;
                   !    used in RHOKM so wet density
   dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        2:bl_levels),                                                   &
                   ! IN Buoyancy gradient across lower
                   !    interface of layer.
   u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                   ! IN Westerly wind component horizontally
                   !    interpolated to P-grid. (m/s)
   v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                   ! IN Southerly wind component horizontally
                   !    interpolated to P-grid. (m/s)
   qw(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),   &
                   ! IN Total water content (kg per kg air).
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                   ! IN Liquid/frozen water temperature (K).
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels), &
                   ! IN Z_UV(K) is height of u level k
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                   ! IN Z_TQ(K) is height of T,Q level k
                   !    NOTE: RI(K) is held at Z_TQ(K-1)
   z0m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                   ! IN Roughness length for momentum (m).
   h_blend(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                   ! IN Blending height for effective
                   !    roughness length scheme
   dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         2:bl_levels),                                                  &
                   ! IN Modulus of wind shear.
   ri(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels), &
                   ! IN Local Richardson number.
   v_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                   ! IN Surface friction velocity  (m/s)
   rneutml_sq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              bl_levels),                                               &
                   ! IN Smagorinsky mixing length scale squared 
                   !    (on theta-levels)
   delta_smag(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                   ! IN grid size used in smagorinsky length scale
   fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                   ! IN Surface buoyancy flux over density (m^2/s^3).
   flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   !  IN Land fraction on all tiles.

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN) ::                                                &
    nSCMDpkgs             ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
    L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

! Declaration of new BL diagnostics.
  TYPE (strnewbldiag) :: BL_diag

! INOUT arguments
  LOGICAL, INTENT(INOUT) ::                                             &
   cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! INOUT Flag for boundary layer cumulus.
                   !       Can only be changed if ISHEAR_BL=1

  REAL, INTENT(INOUT) ::                                                &
   weight_1dbl(pdims_s%i_start:pdims_s%i_end,                           &
               pdims_s%j_start:pdims_s%j_end,bl_levels)
                   ! INOUT Weighting applied to 1D BL scheme 
                   !       to blend with Smagorinsky scheme, 
                   !       index k held on theta level (k-1)
! OUT arguments
  REAL, INTENT(OUT) ::                                                  &
   rhokm(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,   &
         bl_levels),                                                    &
                   ! OUT Layer K-1 - to - layer K exchange coefficient
                   !       for momentum, on UV-grid with first and last
                   !       levels set to "missing data"
   rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels),                                                    &
                   ! OUT Layer K-1 - to - layer K exchange coefficient
                   !       for scalars (but currently on th-levels)
                   ! On OUT: still to be multiplied by rho(if lq_mix_bl)
                   !         and, for Ri-based scheme, interpolated to
                   !         rho levels in BDY_EXPL2
   zh_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! OUT Mixing layer height (m).

  INTEGER, INTENT(OUT) ::                                               &
   ntml_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! OUT Number of model layers in the turbulently
                   !     mixed layer as determined from the local
                   !     Richardson number profile.

  REAL, INTENT(OUT) ::                                                  &
   lambda_min,                                                          &
                   ! OUT Min value of length scale LAMBDA.
   fm_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                   ! OUT stability function for momentum transport.
                   !     level 1 value is dummy for use in diagnostics
   fh_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                   ! OUT stability function for heat and moisture.
                   !     level 1 value is dummy for use in diagnostics
   elm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),&
                   ! OUT Mixing length for momentum
   elh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),&
                   ! OUT Mixing length for scalars on theta levels
   elh_rho(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels)
                   ! OUT Mixing length for scalars on rho levels

!-----------------------------------------------------------------------
!    Local and other symbolic constants :-





  REAL :: grcp
  REAL :: eh,em,g0,dh,dm,a_lambda

  PARAMETER (                                                           &
   eh=25.0,                                                             &
                   ! Used in calc of stability function FH.
   em=4.0,                                                              &
                   ! Used in calc of stability function FM.
   a_lambda=2.0,                                                        &
                   ! used in calc of LAMBDA_EFF
   grcp=g/cp                                                            &
                   ! Adiabatic lapse rate.
             )

!  Equilibrium SBL model constants
  REAL ::    RtestMin
  INTEGER :: gn,NGz,kMINh
  PARAMETER (                                                           &
   RtestMin=0.0,                                                        &
                   ! Threshold for Rtest
   gn=19,                                                               &
                   ! Size of "G"-tables (No. HonL values)
   NGz=90,                                                              &
                   ! No. z/h steps in each "G" integration
   kMINh=2                                                              &
                   ! Level of minimum SBL height (>=2)
  )

!  Define local storage.
  REAL ::                                                               &
   RiCrit(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                   ! Critical Richardson number
   func(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                   ! 2D variable for SBL stabiliy function options
   sharp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                   ! 2D variable for SHARP stabiliy function
   invMOsurf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                   ! Inverse of sfce M-O length
                   ! Note: Inverse is used so that neutral conditions
                   !       can be handled (M-O length --> infinity)
   zh_esbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                   ! Ht of equilib SBL (sub-grid
   HLtab(gn),                                                           &
                   ! Lookup tables (Gx calcs
   GHsav(gn,NGz),gmsav(gn,NGz),                                         &
                   ! in equilib SBL scheme)
   THv_TQ(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                   ! Virtual potential temperature on theta levels
   THv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                   ! THv_TQ interpolated to U,V levels
   prandtl_number(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end), &
                   ! = KM/KH
   BL_weight(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             bl_levels),                                                &
                   ! Fractional weight applied to
                   ! BL function, vs free atmos
   turb_length(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               2:bl_levels)
                   ! Turbulent length scale on theta levels, 
                   ! indexed as Ri (m)

  INTEGER ::                                                            &
   ntml_esbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! No. UV-levels inside equilibrium SBL

  LOGICAL ::                                                            &
   topbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! Flag for having reached the
                                  ! top of the turbulently mixed layer.

! Variables for stability function tails
  REAL ::                                                               &
   fm_louis,                                                            &
                   ! FM calculated using Louis
   fm_sharpest,                                                         &
                   ! FM calculated using SHARPEST
   z_scale,                                                             &
                   ! Scale height for interpolation
   g0_orog,                                                             &
                   ! Orog dependent version of G0
   zpr
                   ! z/sigma_h

! Variables for boundary layer depth based formulation
  REAL ::                                                               &
    h_tkeb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                   ! TKE budget based BL depth
    MOsurf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                   ! surface Obukhov length
    diff_min(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

  REAL ::                                                               &
   h_est,                                                               &
   rifb,                                                                &
                    ! Bulk flux Richardson number
   pr_n,                                                                &
                    ! neutral Prandtl number
   r_pr_n,                                                              &
                    ! 1 / neutral Prandtl number
   m_tau,m_buoy,                                                        &
                    ! Indices for implied stress and buoyancy flux profs
   ind,diff

  REAL ::                                                               &
   subb, subc, subg, subh, ric, ricinv, rifac
                    !Constants for LEM stability functions

  REAL ::                                                               &
   f_log,                                                               &
                 ! Temporary in calculation of logarithmic correction
   fh,                                                                  &
                 ! (Value of) stability function for heat & moisture.
   fm,                                                                  &
                 ! (Value of) stability function for momentum transport.
   rtmri,                                                               &
                 ! Temporary in stability function calculation.
   vkz,                                                                 &
                 ! Temporary in calculation of ELH.
   lambdam,                                                             &
                 ! Asymptotic mixing length for turbulent transport
!                  of momentum.
   lambdah,                                                             &
                 ! Asymptotic mixing length for turbulent transport
!                  of heat/moisture.
   lambdah_rho,                                                         &
                 ! Asymptotic mixing length for turbulent transport
!                  of heat/moisture on rho levels
   rlambda_fac,                                                         &
                 ! reciprocal of lambda_fac
   beta,                                                                &
                 ! empirical factor multiplying zh/delta
   zht,                                                                 &
                 ! top of boundary layer mixing
   zfa,                                                                 &
                 ! height to use beta_fa in blendin
   lambda_eff    ! Effective mixing length used with effective
!                  roughness length scheme.

      !Equilibrium SBL model temporary real scalar variables
  REAL ::                                                               &
   Ztop,Zbot,zz,Zprev,dz,hh,                                            &
                                       ! Height variables
   u1,u2,Ubot,Ujmp,                                                     &
                                       ! Velocity variables
   t1,t2,Tbot,Tjmp,                                                     &
                                       ! Temperature variables
   Rib,Rilim,Rtest,                                                     &
                                       ! Ri variables
   UWsfce, WTsfce, invMOsfce,                                           &
                                       ! Surface flux variables
   USequil,WTequil,invMOequil,                                          &
                                       ! Prescribed flux profiles
    km, kh, kuu, ktt, kut, ktu,                                         &
                                       ! Full/partial eddy diffusivities
   pkm,pkh,pkuu,pktt,pkut,pktu,                                         &
                                       ! Scaled arguments of SBLequil
   Zhat,lhat,PHIe,PHIw,RiSBL,                                           &
                                       ! Scaled arguments of SBLequil
   gm,gh,Gs,Gdz,Gz,g1,g2,h1,h2,                                         &
                                       ! Temporaries for Gx calcs
   GHtab1,GHtab2,GMtab1,GMtab2,                                         &
                                       ! Temporaries for Gx calcs
   ce,                                                                  &
                                       ! Constant in dissipation param.
                                       ! (output by SBLequil)
   rpow,cb,cn,                                                          &
                                       ! Constants in mixing length eqn
                                       ! (output by SBLequil)
   tmp,HonL,dbdu,GonTv,ll,Mlamb,                                        &
                                       ! Miscellaneous temporaries
   slope,dz100,Zupper,fG0          ! Miscellaneous temporaries

  INTEGER ::                                                            &
   i,j,                                                                 &
                 ! Loop counter (horizontal field index).
   k, kl,                                                               &
                 ! Loop counters (vertical level index).
   kb, kt,                                                              &
                 ! Base and top level of unstable Ri layers
   mbl           ! Maximum number of model layers below mixed layer top.

      !Equilibrium SBL model temporary integer scalar variables
  INTEGER ::                                                            &
   kZtop,kZbot,GK,kG0,                                                  &
                 ! Temporary loop counters
   iERRSBL   ! SBLequil error status

     !Equilibrium SBL model logical variables
  LOGICAL ::                                                            &
   GcalcDone,                                                           &
                 ! Calculation of Gx values has been performed
   PrevMin,                                                             &
                 ! Previous Utest minimium has been detected
   subcrit,                                                             &
                 ! flag for being in a subcritical ri layer
   subgrid       ! Will perform subgrid SBL depth calculation

     !Switch to enable subgrid SBL depth diagnosis
  LOGICAL ::    sg_enabled
  PARAMETER (sg_enabled=.TRUE.)

     !Equilibrium SBL model SAVED variables
  SAVE HLtab,GHsav,gmsav,GcalcDone

     !Equilibrium SBL model DATA statements
  DATA HLtab /0.0001,0.001,0.002,0.005,0.01,0.02,0.05,                  &
              0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,                   &
              100.0,200.0,500.0/
  DATA GcalcDone /.FALSE./

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('EX_COEF',zhook_in,zhook_handle)
!-----------------------------------------------------------------------
!   IF stochastic physics random parameters is used set the parameter
!   used to vary the stability function to a perturbed value, if not
!   use the standard setting.
!-----------------------------------------------------------------------
  IF (l_rp2) THEN
    g0=g0_rp
  ELSE
    g0=10.0
  END IF

  dh=g0/eh                 ! Used in calc of stability function FH.
  dm=g0/em                 ! Used in calc of stability function FM.

!---------------------------------------------------------------
! Set neutral and default Prandtl number (Pr=KM/KH)
!---------------------------------------------------------------
  pr_n = 1.0
  IF (Prandtl == LockMailhot2004) pr_n = 0.7
  IF (sbl_op  ==  depth_based)    pr_n = 0.7
! Use pr_n=0.7 if any LEM stability function selected
  IF (sbl_op == LEM_stability .OR. cbl_op == LEM_std                    &
                              .OR. cbl_op == LEM_conven ) pr_n = 0.7
  r_pr_n = 1.0 / pr_n
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      prandtl_number(i,j) = pr_n
    END DO
  END DO

! Settings for LEM stability functions
  IF (cbl_op == LEM_conven) THEN
    ! the "conventional" subgrid model, Brown (1999)
    subb = 1.43
    subc = 1.43
  ELSE 
    ! the "standard" LEM subgrid model, Brown (1999)
    subb = 40.0
    subc = 16.0
  END IF
  subg = 1.2
  ric = 0.25
  ricinv = 1./ric

!  Set MBL, "maximum number of boundary layer levels" for the purposes
!  of mixed layer height calculation.

  mbl = bl_levels - 1

!-----------------------------------------------------------------------
! 0. Initialise flag for having reached top of turbulently mixed layer
!-----------------------------------------------------------------------
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      topbl(i,j)     = .FALSE.
    END DO
  END DO
!-----------------------------------------------------------------------
! Set-up a BL weighting function, =1 near the ground (ie in the BL)
!                                 =0 in the free troposphere
! Rate and height at which transition occurs varys depending on choices:
!-----------------------------------------------------------------------
  DO k=1,bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_weight(i,j,k) = 1.0
      END DO
    END DO
  END DO

  IF (local_fa == to_sharp_across_1km) THEN
!---------------------------------------------------------
! Additional code to allow the local Ri scheme to use 
! SHARPEST in the free atmosphere, ie above the BL top,
! regardless of the tail option selected above.
! Set Z_SCALE to 1km to mimic old value of BL_LEVELS, 
!  gives BL_weight~0 by 2km, ~0.95 at 500m
!---------------------------------------------------------
    z_scale = 1000.0
    DO k=2,bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          zpr = z_tq(i,j,k-1)/z_scale
          BL_weight(i,j,k) = 0.5*( 1.0 - TANH(3.*(zpr-1.0) ) )
        END DO
      END DO
    END DO
  END IF

  IF (sg_orog_mixing /= OFF) THEN
!-----------------------------------------------------------------
! Subgrid orographic height dependence for SBL tail (option 1) 
! or orographic dependence of mixing lengths, lambdam,h (opt 2)
! Gives BL_weight~[1,0.95,0.5,0] at ZPR=[0,0.6,1,1.7] 
!----------------------------------------------------------------
    DO k=2,bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF (sigma_h(i,j) > 1.0 ) THEN
            zpr = z_tq(i,j,k-1)/sigma_h(i,j)
            BL_weight(i,j,k) = 0.5*( 1.0 - TANH(4.0*(zpr-1.0) ) )
          END IF
        END DO
      END DO
    END DO
  END IF 
!-----------------------------------------------------------------------
!  Set LAMBDA_MIN 
!-----------------------------------------------------------------------
  IF (l_rp2) THEN 
    lambda_min = lambda_min_rp 
  ELSE 
    lambda_min = lambda_min_nml
  END IF

  DO k=2,bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        turb_length(i,j,k)=lambda_min
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
! Set critical Richardson number
!-----------------------------------------------------------------------
  IF (L_RP2) THEN 
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        RiCrit(i,j) = RiCrit_rp 
      END DO
    END DO
  ELSE 
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        ! Default critical Ri for Long_tails and Louis
        RiCrit(i,j) = 1.0
      END DO
    END DO
  END IF

  IF (Variable_RiC == on) THEN

    SELECT CASE (sbl_op)

!--------------------------------------------
! SHARP TAILS
!--------------------------------------------
    CASE(sharpest)

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          RiCrit(i,j) = RiCrit_sharp
        END DO
      END DO

!--------------------------------------------
! LEM TAILS
!--------------------------------------------
    CASE(LEM_stability)

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          RiCrit(i,j) = ric
        END DO
      END DO

!--------------------------------------------
! SHARP over sea; longer tails over land
!--------------------------------------------
    CASE(sharp_sea_long_land, sharp_sea_mes_land,                       &
         Sharp_sea_Louis_land)

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF (flandg(i,j) < 0.5) THEN
                ! SHARPEST over sea
            RiCrit(i,j) = RiCrit_sharp
          ELSE
                ! Longer tails over land
            IF (L_RP2) THEN
              RiCrit(i,j) = RiCrit_rp
            ELSE
              RiCrit(i,j) = 1.0
            END IF
          END IF
        END DO
      END DO

    END SELECT ! SBL_OP

  END IF
!-----------------------------------------------------------------------
! Initialise 3D stability functions
!-----------------------------------------------------------------------
  IF (l_subfilter_vert .OR. l_subfilter_horiz) THEN
    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          fm_3d(i,j,k) = 0.0
          fh_3d(i,j,k) = 0.0
        END DO
      END DO
    END DO
  END IF
!-----------------------------------------------------------------------
!  1.1 Loop over levels calculating Richardson numbers.
!-----------------------------------------------------------------------

  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

!-----------------------------------------------------------------------
! 1.2 If either a stable layer (Ri>RiCrit) or the maximum BL
!     height has been reached, set boundary layer height (ZH_LOCAL) to
!     the height of the lower boundary of the current layer
!-----------------------------------------------------------------------
        IF ( .NOT.topbl(i,j) .AND.                                      &
             (ri(i,j,k) >  RiCrit(i,j) .OR. k >  mbl) )  THEN
          topbl(i,j) = .TRUE.
          IF (local_fa >= ntml_level_corrn) THEN
            ! Ri(k)>RiC => theta-level(k-1) is supercrit => NTML=k-2
            ntml_local(i,j) = MAX( 1, k-2 )
          ELSE
            ntml_local(i,j) = k-1
          END IF
          zh_local(i,j) = z_uv(i,j,ntml_local(i,j)+1)
        END IF
      END DO  ! Loop over points
    END DO  ! Loop over points
  END DO  ! Loop over levels
! Save original diagnosis
  IF (BL_diag%l_zhlocal) THEN
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%zhlocal(i,j)=zh_local(i,j) 
      END DO  ! Loop over points
    END DO  ! Loop over levels
  END IF
!-----------------------------------------------------------------------
! If NTML_LOCAL is greater than the top of the parcel ascent (NTPAR)
! for a cumulus-capped layer, shear driven mixing is allowed to
! dominate (if ISHEAR_BL=1 selected)
!-----------------------------------------------------------------------
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( ishear_bl  ==  1 .AND.                                       &
           ntml_local(i,j)  >   ntpar(i,j) ) THEN
        cumulus(i,j) = .FALSE.
      END IF
    END DO
  END DO
!-----------------------------------------------------------------------
! In CUMULUS layers the local scheme is capped at the LCL (given in
! this case by NTML_NL).  Save local BL depth as SCM diagnostic.
!-----------------------------------------------------------------------
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( cumulus(i,j) ) THEN
        ntml_local(i,j) = ntml_nl(i,j)
        zh_local(i,j) = z_uv(i,j,ntml_local(i,j)+1)
      END IF
    END DO
  END DO
!-----------------------------------------------------------------------
!! 1.3 Search for sub-critical layers above the PBL and set the 
!      mixing length to scale with these layer depths
!-----------------------------------------------------------------------
  IF (local_fa == free_trop_layers) THEN
    rlambda_fac=1.0/lambda_fac
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        subcrit = .FALSE.
        DO k = 3, bl_levels

          IF ( k > ntml_local(i,j)+1  .AND.                             &
               ! we know Ri(ntml_local(i,j)+2) > RiCrit
               ri(i,j,k) < RiCrit(i,j) .AND. .NOT.subcrit ) THEN
            kb      = k   ! first level of subcritical Ri in layer
            subcrit = .TRUE.
          END IF
          IF (ri(i,j,k) >= RiCrit(i,j) .AND. subcrit ) THEN
            kt      = k-1 ! last level of subcritical ri 
            subcrit = .FALSE.
            !---------------------------------------------------------
            ! turb_length(k) is held, with Ri(k), on th-level(k-1)
            !---------------------------------------------------------
            DO kl=kb,kt
              turb_length(i,j,kl) = z_uv(i,j,kt) - z_uv(i,j,kb-1)
!              turb_length(i,j,kl) = 1.0/(                              &
!                          1.0/(z_tq(i,j,kl-1) - z_uv(i,j,kb-1)) +      &
!                          1.0/(z_uv(i,j,kt)   - z_tq(i,j,kl-1))   )
              turb_length(i,j,kl) = MAX( lambda_min*rlambda_fac,        &
                MIN(turb_length(i,j,kl),lambda_max_nml*rlambda_fac)   )
            END DO
          END IF

        END DO
      END DO
    END DO
  END IF
!-----------------------------------------------------------------------
! When blending with the Smagorinsky scheme set turb_length in the BL
! and  use the DSC layer depth as the length scale within a DSC layer 
! Remember turb_length(k) is held, with Ri(k), on th-level(k-1)
!-----------------------------------------------------------------------
  IF (l_subfilter_blend) THEN
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF ( k-1 <= MAX(ntml_nl(i,j),ntml_local(i,j)) ) THEN
            turb_length(i,j,k) =  MAX( turb_length(i,j,k),              &
                MAX( z_uv(i,j,ntml_nl(i,j)+1), zh_local(i,j) ) )
          END IF
          IF ( k-1 >= nbdsc(i,j) .AND. k-1 <= ntdsc(i,j) ) THEN
            turb_length(i,j,k) = MAX( turb_length(i,j,k),               &
                    ( z_uv(i,j,ntdsc(i,j)+1)-z_uv(i,j,nbdsc(i,j)) ) )
          END IF
        END DO
      END DO
    END DO
  END IF
!-----------------------------------------------------------------------
!! 2.  Richardson Number based local mixing scheme
!-----------------------------------------------------------------------
!! 2.0 Loop round "boundary" levels; calculate the stability-
!!     dependent turbulent mixing coefficients.
!-----------------------------------------------------------------------
! TKE budget based depth diagnosis

! Starting with the definition of the flux Richardson
! number, assuming similarity profiles for
! stress and buoyancy flux, and vertically integrating
! gives an expression for the stable boundary layer
! depth which is based just on surface fluxes and
! the wind speed change across the boundary layer.
! ----------------------------------------------------

  IF (sbl_op  ==  depth_based) THEN

! Index for assumed buoyancy profile
    m_buoy=1.0

! Index for assumed stress profile
    m_tau=1.0

! Effective bulk flux Richardson number
    rifb=0.3

    ind=m_buoy-m_tau+1.0

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

! Set diff_min to a large initial value
        diff_min(i,j)=1000.0

! Surface Obukhov length
        MOsurf(i,j)= -v_s(i,j)*v_s(i,j)*v_s(i,j)                        &
                    /(vkman*fb_surf(i,j))
      END DO
    END DO

    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

              ! The wind speed change from level k to the surface
          u1=SQRT(u_p(i,j,k)*u_p(i,j,k)+v_p(i,j,k)*v_p(i,j,k))

              ! h_est is the estimate of the stable boundary layer
              ! depth using the TKE based formula
          h_est=vkman*MOsurf(i,j)*ind*rifb*u1/v_s(i,j)

              ! Absolute difference between height and estimate
          diff=ABS(z_uv(i,j,k)-h_est)

              ! If h_est is closer than the previous closest value
              ! (diff_min) reset the h_tkeb to h_est

          IF (diff  <   diff_min(i,j)) THEN
            diff_min(i,j)=diff
            h_tkeb(i,j)=h_est
          END IF

        END DO
      END DO
    END DO

  END IF   ! SBL_OP = Depth_based

! ----------------------------------------------------------------
! Main loop over levels
! ----------------------------------------------------------------
  DO k = 2, bl_levels
! ----------------------------------------------------------------
! Load up 2D array FUNC with selected stability function

!  SBL_OP                 Option

!  Long_tails             Long tails
!  Sharpest               SHARPEST function
!  Sharp_sea_long_land    SHARPEST over sea ; Long tails over land
!  Mes_tails              MESOSCALE model: Louis/SHARPEST blend
!  Louis_tails            Louis function
!  Depth_based            Boundary layer depth based formulation
!  Sharp_sea_mes_land     SHARP over sea; Mes over land
!  Sharp_sea_Louis_land   SHARP over sea; Louis over land
! ----------------------------------------------------------------

    SELECT CASE (sbl_op)

!--------------------------------------------
! LONG TAILS
!--------------------------------------------
    CASE(long_tails)

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          func(i,j)=1.0 / ( 1.0 + g0 * ri(i,j,k) )
        END DO
      END DO

!--------------------------------------------
! SHARP TAILS
!--------------------------------------------
    CASE(sharpest)

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF (ri(i,j,k)  <   1.0/g0) THEN
            func(i,j) = 1.0 - 0.5 * g0 * ri(i,j,k)
          ELSE
            func(i,j) = 1.0 / ( 2.0 * g0 * ri(i,j,k) )
          END IF
          func(i,j)=func(i,j)*func(i,j)
        END DO
      END DO

!--------------------------------------------
! LEM TAILS (cut-off at Ric)
!--------------------------------------------
    CASE(LEM_stability)

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF ( ri(i,j,k) >= 0.0 .AND. ri(i,j,k)< ric )  THEN
            rifac = (1.0-ri(i,j,k)*ricinv)**4
            func(i,j) = rifac*(1.0-subg*ri(i,j,k))
          ELSE IF (ri(i,j,k) >= ric) THEN
            func(i,j) = 0.0
          END IF
        END DO
      END DO

!--------------------------------------------
! SHARP over sea; long tails over land
!--------------------------------------------
    CASE(sharp_sea_long_land)

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF (flandg(i,j) < 0.5) THEN
                ! SHARPEST over sea
            IF (ri(i,j,k)  <   1.0/g0) THEN
              func(i,j) = 1.0 - 0.5 * g0 * ri(i,j,k)
            ELSE
              func(i,j) = 1.0 / ( 2.0 * g0 * ri(i,j,k) )
            END IF
            func(i,j)=func(i,j)*func(i,j)
          ELSE
                ! Long tails over land
            func(i,j)= 1.0 / ( 1.0 + g0 * ri(i,j,k) )
          END IF
        END DO
      END DO

!--------------------------------------------
! MESOSCALE MODEL TAILS
!--------------------------------------------
    CASE(mes_tails)

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
              ! Louis function
          fm = 1.0 / ( 1.0 + 0.5 * g0 * ri(i,j,k) )
          fm_louis = fm * fm
              ! code for SHARPEST
          IF (ri(i,j,k)  <   1.0/g0) THEN
            fm = 1.0 - 0.5 * g0 * ri(i,j,k)
          ELSE
            fm = 1.0 / ( 2.0 * g0 * ri(i,j,k) )
          END IF
          fm_sharpest = fm * fm
              ! Linear weighting function giving Louis
              ! at z=0, SHARPEST above Z_SCALE
          z_scale = 200.0
          IF ( z_tq(i,j,k-1)  >=  z_scale ) THEN
            func(i,j) = fm_sharpest
          ELSE
            func(i,j)= fm_louis *( 1.0 - z_tq(i,j,k-1)/z_scale )        &
                       + fm_sharpest * z_tq(i,j,k-1)/z_scale
          END IF
        END DO
      END DO

!--------------------------------------------
! LOUIS TAILS
!--------------------------------------------
    CASE(louis_tails)

          ! LOUIS FUNCTION
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          func(i,j)=1.0 / ( 1.0 + 0.5 * g0 * ri(i,j,k) )
          func(i,j)=func(i,j)*func(i,j)
        END DO
      END DO

!--------------------------------------------
! LONG TAILS FOR USE WITH DEPTH BASED SCHEME
!--------------------------------------------
    CASE(depth_based)
          ! LONG TAILS
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          func(i,j)=1.0 / ( 1.0 + g0 * ri(i,j,k) )
        END DO
      END DO

!--------------------------------------------
! SHARP TAILS OVER SEA; MES TAILS OVER LAND
!--------------------------------------------
    CASE(sharp_sea_mes_land)
          ! SHARP sea; MES land
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
              ! Louis function
          fm = 1.0 / ( 1.0 + 0.5 * g0 * ri(i,j,k) )
          fm_louis = fm * fm
              ! code for SHARPEST
          IF (ri(i,j,k)  <   1.0/g0) THEN
            fm = 1.0 - 0.5 * g0 * ri(i,j,k)
          ELSE
            fm = 1.0 / ( 2.0 * g0 * ri(i,j,k) )
          END IF
          fm_sharpest = fm * fm
              ! Linear weighting function giving Louis at z=0,
              ! SHARPEST above Z_SCALE
          z_scale = 200.0

          IF (flandg(i,j) < 0.5) THEN
                ! SHARP sea
            func(i,j)=fm_sharpest
          ELSE
                ! MES land
            IF ( z_tq(i,j,k-1)  >=  z_scale ) THEN
              func(i,j) = fm_sharpest
            ELSE
              func(i,j)= fm_louis *( 1.0 - z_tq(i,j,k-1)/z_scale )      &
                         + fm_sharpest * z_tq(i,j,k-1)/z_scale
            END IF

          END IF  ! FLANDG(i,j) < 0.5

        END DO ! loop over i
      END DO ! loop over j

!--------------------------------------------
! SHARP TAILS OVER SEA; LOUIS TAILS OVER LAND
!--------------------------------------------
    CASE(Sharp_sea_Louis_land)
          ! SHARP sea; Louis land
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          IF (flandg(i,j) < 0.5) THEN
                ! SHARP sea
            IF (ri(i,j,k)  <   1.0/g0) THEN
              fm = 1.0 - 0.5 * g0 * ri(i,j,k)
            ELSE
              fm = 1.0 / ( 2.0 * g0 * ri(i,j,k) )
            END IF
            func(i,j)=fm * fm
          ELSE
                ! Louis land
            fm_louis = 1.0 / ( 1.0 + 0.5 * g0 * ri(i,j,k) )
            func(i,j)= (1.0 - WeightLouisToLong) *                      &
                       fm_louis * fm_louis +                            &
                       WeightLouisToLong *                              &
                       1.0 / ( 1.0 + g0 * ri(i,j,k) )
          END IF  ! FLANDG(i,j) < 0.5

        END DO ! loop over i
      END DO ! loop over j

    END SELECT ! SBL_OP

!------------------------------------------------------------------
! Additional code to allow the local Ri scheme to use 
! SHARPEST in the free atmosphere, ie above the BL top,
! regardless of the tail option selected above.
!------------------------------------------------------------------
    IF (local_fa == to_sharp_across_1km) THEN

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          !----------------------------
          ! Calculate SHARPEST function
          !----------------------------
          IF (ri(i,j,k)  <   1.0/g0) THEN
            sharp(i,j) = 1.0 - 0.5 * g0 * ri(i,j,k)
          ELSE
            sharp(i,j) = 1.0 / ( 2.0 * g0 * ri(i,j,k) )
          END IF
          sharp(i,j)=sharp(i,j)*sharp(i,j)

          func(i,j) = func(i,j) * BL_weight(i,j,k)                      &
                    + sharp(i,j)*( 1.0 - BL_weight(i,j,k) )

        END DO
      END DO

    END IF
!------------------------------------------------------------------
! Additional code to allow the local Ri scheme stable mixing to 
! depend the size of subgrid orography.  Overwrites values of FUNC 
! as calculated above
!------------------------------------------------------------------
    IF (sg_orog_mixing == 1) THEN

      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
!-------------------------------------------------------
! SBL tail dependent on subgrid orography
!  - use SHARPEST function but with variable coefficient
!    that reduces to sharpest both with height above 
!    orography and as orography gets smaller
!-------------------------------------------------------
          IF ( sigma_h(i,j) > 0.1 ) THEN
            ! Then additional near-surface orographic dependence
            g0_orog = g0 / ( 1.0 +                                      &
                             (sigma_h(i,j)/25.0)*BL_weight(i,j,k) )

            IF (ri(i,j,k) < 1.0/g0_orog) THEN
              func(i,j) = 1.0 - 0.5 * g0_orog * ri(i,j,k)
            ELSE
              func(i,j) = 1.0 / ( 2.0 * g0_orog * ri(i,j,k) )
            END IF
            func(i,j) = func(i,j)*func(i,j)

          END IF
        END DO
      END DO

    END IF
!---------------------------------------------------------------
! Set stable Prandtl number (=KM/KH)
!---------------------------------------------------------------
    IF (Prandtl == LockMailhot2004) THEN
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          prandtl_number(i,j) = MIN( pr_max,                            &
                                     pr_n*(1.0 + 2.0*ri(i,j,k)) )
        END DO
      END DO
    ELSE IF (sbl_op == LEM_stability) THEN
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF ( ri(i,j,k) >= 0.0 .AND. ri(i,j,k) < ric) THEN
            prandtl_number(i,j) = pr_n/(1.0-subg*ri(i,j,k))
          ELSE IF (ri(i,j,k) >= ric) THEN
            prandtl_number(i,j) = pr_n/(1.0-subg*ric)
          ELSE 
            prandtl_number(i,j) = pr_n
          END IF
        END DO
      END DO
    END IF

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
!-----------------------------------------------------------------------
! 2.1 Calculate asymptotic mixing lengths LAMBDAM and LAMBDAH
!     (may be equal or LambdaM=2*LambdaH (operational setting)).
!-----------------------------------------------------------------------
        IF (l_lambdam2) THEN
          IF (l_rp2) THEN
            lambdam = MAX (lambda_min , 2*par_mezcla*zh_local(i,j))
            lambdah = MAX (lambda_min , par_mezcla*zh_local(i,j))
          ELSE
            lambdam = MAX (lambda_min , 0.30*zh_local(i,j))
            lambdah = MAX (lambda_min , 0.15*zh_local(i,j))
          END IF
        ELSE
          IF (l_rp2) THEN
            lambdam = MAX ( lambda_min , par_mezcla*zh_local(i,j) )
          ELSE
            lambdam = MAX ( lambda_min , lambda_fac*zh_local(i,j) )
          END IF
          lambdah = lambdam
        END IF

        IF (sg_orog_mixing == 3) THEN   
!-----------------------------------------------------------------------
! Use orographic mixing length for heat too, and reduce both 
! above sigma_h smoothly
! NOTE: THIS CODE WILL NOT ENHANCE LAMBDAH because it only
! gets used in bdy_expl2 where the calculation is redone as 
! standard - this was a mistake but is now operational in 
! the UKV!
!-----------------------------------------------------------------------
          IF (k >= ntml_local(i,j)+2) THEN
            lambdam = lambda_min
            lambdah = lambda_min
          END IF
          lambdah    = MAX (lambdah,                                    &
                              BL_weight(i,j,k)*a_lambda*h_blend(i,j) )
          lambda_eff = MAX (lambdam,                                    &
                              BL_weight(i,j,k)*a_lambda*h_blend(i,j) )
        ELSE
          lambda_eff = MAX (lambdam, a_lambda*h_blend(i,j) )
!------------------------------------------------------------
! Optionally reduce mixing length above local BL top
!------------------------------------------------------------
          IF (k >= ntml_local(i,j)+2 .AND. .NOT.l_full_lambdas) THEN
            lambdam = lambda_min
            lambdah = lambda_min
            IF (z_tq(i,j,k-1) > a_lambda*h_blend(i,j))                  &
                                              lambda_eff=lambda_min
          END IF
          IF ( k >= ntml_local(i,j)+2 .AND. l_full_lambdas .AND.        &
             local_fa == to_sharp_across_1km) THEN
              ! Weight lambda to lambda_min with height
              ! Assuming only local_fa == to_sharp_across_1km will have 
              ! L_FULL_LAMBDAS. If other LOCAL_FA options are coded here 
              ! then changes must be included in section 5.3 of bdy_expl2

            lambda_eff = lambda_eff * BL_weight(i,j,k)                  &
                       + lambda_min*( 1.0 - BL_weight(i,j,k) )
            lambdah    = lambdah * BL_weight(i,j,k)                     &
                       + lambda_min*( 1.0 - BL_weight(i,j,k) )
          END IF
        END IF

        lambdah_rho  = lambdah

        IF ( local_fa == free_trop_layers ) THEN
          lambda_eff = MAX( lambda_eff, lambda_fac*turb_length(i,j,k) )
          lambdah    = MAX( lambdah,    lambda_fac*turb_length(i,j,k) )
          ! lambdah_rho does not need to be recalculated under 
          ! local_fa option "free_trop_layers" as the full KH profile 
          ! will be interpolated in bdy_expl2
        END IF
!-----------------------------------------------------------------------
! 2.2 Calculate mixing lengths ELH, ELM coincident with RI(K) and so
!     at Z_TQ(K-1)
!-----------------------------------------------------------------------
!  Incorporate log profile corrections to the vertical finite
!  differences into the definitions of ELM and ELH.
!  Note that ELH_RHO is calculated (on rho levels) for direct inclusion 
!  in RHOKH and also (as elh) on theta levels for the unstable 
!  stability functions and inclusion in RHOKH before interpolation 
!  (under local_fa option "free_trop_layers").
!  To save computing logarithms for all K, the values of ELM and ELH
!  are unchanged for K > K_LOG_LAYR.

        IF (k  <=  k_log_layr) THEN
          vkz   = vkman * ( z_uv(i,j,k) - z_uv(i,j,k-1) )
          f_log = LOG( ( z_uv(i,j,k) + z0m(i,j)   ) /                   &
                       ( z_uv(i,j,k-1) + z0m(i,j) ) )
          elm(i,j,k) = vkz / ( f_log + vkz/lambda_eff )
          elh(i,j,k) = vkz / ( f_log + vkz/lambdah )
          vkz   = vkman * ( z_tq(i,j,k) - z_tq(i,j,k-1) )
          f_log = LOG( ( z_tq(i,j,k) + z0m(i,j)   ) /                   &
                       ( z_tq(i,j,k-1) + z0m(i,j) ) )
          elh_rho(i,j,k) = vkz / ( f_log + vkz/lambdah_rho )
        ELSE
          vkz = vkman * ( z_tq(i,j,k-1) + z0m(i,j) )
          elm(i,j,k) = vkz / (1.0 + vkz/lambda_eff )
          elh(i,j,k) = vkz / (1.0 + vkz/lambdah )
          vkz = vkman * ( z_uv(i,j,k) + z0m(i,j) )
          elh_rho(i,j,k) = vkz / (1.0 + vkz/lambdah_rho )
        END IF

        IF (l_subfilter_blend) THEN

          zz = z_tq(i,j,k-1)  ! height of rhokm(k)
          ! Estimate zht = interface between BL and FA, used to set the
          ! blending rate, beta
          zht = MAX( z_uv(i,j,ntml_nl(i,j)+1) , zh_local(i,j) )
          IF (ntdsc(i,j) > 0) zht = MAX( zht, z_uv(i,j,ntdsc(i,j)+1) )
          zfa=zht+1000.0
          IF (zz <= zht) THEN
            beta=beta_bl
          ELSE IF (zz <= zfa) THEN
            beta = beta_bl*(zfa-zz)/(zfa-zht) +                         &
                   beta_fa*(zz-zht)/(zfa-zht)
          ELSE
            beta=beta_fa
          END IF
          ! turb_length is the greater of the local and non-local 
          ! BL depths up to that bl top
          z_scale = MAX( zz, turb_length(i,j,k) )
          ! Need to restrict z_scale to dsc depth within a dsc layer 
          ! (given by turb_length) and to distance from dsc top below the 
          ! dsc layer
          IF ( k-1 <= ntdsc(i,j) )  THEN
            z_scale = MIN( z_scale,                                     &
                MAX( turb_length(i,j,k), z_uv(i,j,ntdsc(i,j)+1)-zz ) )
          END IF

          ! Finally calculate 1D BL weighting factor
          weight_1dbl(i,j,k) =                                          &
                   1.0 - TANH( beta*z_scale/delta_smag(i,j)) *          &
                       MAX( 0.0, (4.0-delta_smag(i,j)/z_scale) )/4.0

          elm(i,j,k) = elm(i,j,k)*weight_1dbl(i,j,k) +                  &
                      SQRT(rneutml_sq(i,j,k-1))*(1.0-weight_1dbl(i,j,k))
          elh(i,j,k) = elh(i,j,k)*weight_1dbl(i,j,k) +                  &
                      SQRT(rneutml_sq(i,j,k-1))*(1.0-weight_1dbl(i,j,k))

        END IF  ! test on l_subfilter_blend

        IF (BL_diag%L_elm3D) BL_diag%elm3D(i,j,k)=elm(i,j,k)

!-----------------------------------------------------------------------
! 2.4 Calculate (values of) stability functions FH, FM.
!-----------------------------------------------------------------------

        IF (ri(i,j,k)  >=  0.0) THEN
            !---------------------------------------------------------
            ! Set FM to the standard requested stability function and
            ! scale FH by the neutral Prandtl number,
            ! rather than keep FH the same and scale FM.
            !---------------------------------------------------------
          fm = func(i,j)
          fh = fm * r_pr_n

            !-----------------------------------------------------------
            ! Then, rescale FM by the full Prandtl number.
            ! The reason for coding in this way is to give FM=1 at
            ! neutral and SHARP+Prandtl gives FH~1/Ri^2 and FM~1/Ri
            !-----------------------------------------------------------
          fm = fh * prandtl_number(i,j)

        ELSE   ! ri < 0
          IF (cbl_op == neut_cbl) THEN
            ! Use neutral stability for unstable mixing
            fm = 1.0 
            fh = r_pr_n
          ELSE IF (cbl_op == LEM_std .OR. cbl_op == LEM_conven) THEN
            fm = SQRT(1.0-subc*ri(i,j,k))
            fh = SQRT(1.0-subb*ri(i,j,k)) * r_pr_n
          ELSE
!           ! UM_std
            rtmri = (elm(i,j,k)/elh(i,j,k)) * SQRT ( -ri(i,j,k) )
            fm = 1.0 - ( g0*ri(i,j,k) / ( 1.0 + dm*rtmri ) )
            fh = ( 1.0 - ( g0*ri(i,j,k) / ( 1.0 + dh*rtmri ) ) ) * r_pr_n
          END IF
        END IF

!-----------------------------------------------------------------------
! 2.5 Calculate exchange coefficients RHO*KM(K), RHO*KH(K)
!     both on TH-level K-1 at this stage (RHOKH will be interpolated
!     onto uv-levels and then be multiplied by ELH in BDY_EXPL2 if 
!     local_fa is not "free_trop_layers")
!-----------------------------------------------------------------------

        IF (l_subfilter_vert .OR. l_subfilter_horiz) THEN
          fm_3d(i,j,k)=fm
          fh_3d(i,j,k)=fh
        END IF
        rhokm(i,j,k) = rho(i,j,k-1) * elm(i,j,k) * elm(i,j,k)           &
                                    * dvdzm(i,j,k) * fm
        IF (lq_mix_bl) THEN
            ! Note "RHO" here is always wet density (RHO_TQ) so
            ! save multiplication of RHOKH to after interpolation
          rhokh(i,j,k) =                elm(i,j,k) * dvdzm(i,j,k) * fh
        ELSE
          rhokh(i,j,k) = rho(i,j,k-1) * elm(i,j,k) * dvdzm(i,j,k) * fh
        END IF
        ! If using the FA mixing length profile it is simplest to 
        ! interpolate the full KH profile, including elh (in bdy_expl2)
        IF (local_fa == free_trop_layers)                               &
                    rhokh(i,j,k) = rhokh(i,j,k) * elh(i,j,k)
! -------------------------------------------
! Boundary layer depth based formulation
! -------------------------------------------

        IF (sbl_op  ==  depth_based .AND.                               &
            fb_surf(i,j)  <=  0.0) THEN
          IF (z_tq(i,j,k-1) < h_tkeb(i,j)) THEN

              ! Formula for diffusion coefficient
              ! see Beare et al 2006, Boundary layer Met.

            km = v_s(i,j) * vkman * z_tq(i,j,k-1) *                     &
                          ( (1.0-z_tq(i,j,k-1)/h_tkeb(i,j))**(1.5) )    &
                          /  (1.0 + 4.7*z_tq(i,j,k-1)/MOsurf(i,j))
            rhokm(i,j,k)= rho(i,j,k-1) * km
            IF (lq_mix_bl) THEN
                ! Note "RHO" here is always wet density (RHO_TQ) so
                ! save multiplication of RHOKH to after interpolation
              rhokh(i,j,k)= km*r_pr_n / elm(i,j,k)
            ELSE
              rhokh(i,j,k) = rho(i,j,k-1)*km*r_pr_n /elm(i,j,k)
            END IF
          ELSE
            rhokm(i,j,k)=0.0
            rhokh(i,j,k)=0.0
          END IF

        END IF   !SBL_OP  ==  Depth_based

      END DO !I
    END DO !j
  END DO ! bl_levels
!-----------------------------------------------------------------------
! 3.  Equilibrium Stable Boundary Layer (SBL) model.
!-----------------------------------------------------------------------
  IF (sbl_op==Equilibrium_SBL) THEN

!-----------------------------------------------------------------------
! 3.1 On first timestep only, calculate "G"-tables.
!     Note that GcalcDone, HLtab, GMsav and GHsav are protected
!     by a SAVE statement, so their values are retained
!     in subsequent calls to this subprogram. Gx values
!     are calculated, for each value of HonL, as contributions
!     to integrals of z/h from 0.1 to 1.0 in NGz steps of size Gdz
!-----------------------------------------------------------------------
    IF (.NOT.GcalcDone) THEN
      Mlamb=1.5*Muw_SBL-Mwt_SBL
      Gdz=0.9/REAL(NGz)
      DO GK = 1, gn
        HonL=HLtab(GK)
        Gz =0.1-0.5*Gdz
        DO k = 1, NGz
          Gz  =Gz+Gdz
          Zhat=Gz*HonL/((1.0-Gz)**(Mlamb))
! DEPENDS ON: sblequil
          CALL sblequil(Zhat,lhat,pkm,pkh,pkuu,pktt,pkut,pktu,          &
                        PHIe,PHIw,RiSBL,ce,rpow,cb,cn,iERRSBL)
          tmp=pkh*lhat*vkman*Zhat
          gh=Gdz*((1.0-Gz)**(2.0*(Mwt_SBL-Muw_SBL)))/tmp
          tmp=pkm*lhat*vkman*Zhat
          gm=Gdz*((1.0-Gz)**(     Mwt_SBL-Muw_SBL ))/tmp
          GHsav(GK,k)=gh
          gmsav(GK,k)=gm
        END DO
      END DO
      GcalcDone=.TRUE.
    END IF

!-----------------------------------------------------------------------
! 3.2 Diagnose depth of equilibrium SBL (ZH_ESBL) on UV-levels (Z_UV).
!     Estimate for ZH_ESBL obtained by sub-grid interpolation.
!     NTML_ESBL is then the first UV-level below ZH_ESBL
!-----------------------------------------------------------------------
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        !Calculate theta_v and interpolate to UV grid
        DO k = 1, bl_levels
          THv_TQ(i,j,k)=(tl(i,j,k)+grcp*z_tq(i,j,k))                    &
                          *(1.0+c_virtual*qw(i,j,k))
        END DO
        k=1
        t1=THv_TQ(i,j,k)
        t2=THv_TQ(i,j,k+1)
        dz=z_tq(i,j,k+1)-z_tq(i,j,k)
        slope=(t2-t1)/dz
        THv(i,j,k+1)=t1+slope*(z_uv(i,j,k+1)-z_tq(i,j,k)) !Lin interp
        THv(i,j,1)=t1-slope*(z_tq(i,j,1)-z_uv(i,j,1))
        DO k = 2, bl_levels-1
          t1=THv_TQ(i,j,k)
          t2=THv_TQ(i,j,k+1)
          dz=z_tq(i,j,k+1)-z_tq(i,j,k)
          slope=(t2-t1)/dz
          THv(i,j,k+1)=t1+slope*(z_uv(i,j,k+1)-z_tq(i,j,k)) !Lin interp
        END DO
        invMOsurf(i,j)=-vkman*fb_surf(i,j)                              &
                     /(v_s(i,j)*v_s(i,j)*v_s(i,j))
        invMOsfce=invMOsurf(i,j)

        topbl(i,j)    =.FALSE.
        zh_esbl(i,j)  =z_uv(i,j,kMINh)+1.0
        ntml_esbl(i,j)=kMINh

        IF (fb_surf(i,j) <= 0.0) THEN !only for stable/neutral cases
          subgrid =.FALSE.
          kZtop   =kMINh-1

          DO WHILE ((.NOT.topbl(i,j)).AND.(kZtop <= bl_levels-1))
            Zprev=z_uv(i,j,kZtop)
            kZtop=kZtop+1
            Ztop =z_uv(i,j,kZtop)
            Zbot =0.1*Ztop !Height of surface layer top
            kZbot=0
            zz   =0.0
            DO WHILE (zz <  Zbot) !find kZbot (UV-level above Zbot)
              kZbot=kZbot+1
              zz=z_uv(i,j,kZbot)
            END DO !find kZbot
            IF (kZbot >  1) THEN
              !Interpolation of U and T to Zbot
              k=kZbot-1
              u1=SQRT(u_p(i,j,k)*u_p(i,j,k)+v_p(i,j,k)*v_p(i,j,k))
              u2=SQRT(u_p(i,j,k+1)*u_p(i,j,k+1)                         &
                     +v_p(i,j,k+1)*v_p(i,j,k+1))
              t1=THv(i,j,k)
              t2=THv(i,j,k+1)
              dz=z_uv(i,j,k+1)-z_uv(i,j,k)
              Ubot=u1+(u2-u1)*(Zbot-z_uv(i,j,k))/dz !Linear interp
              Tbot=t1+(t2-t1)*(Zbot-z_uv(i,j,k))/dz !Linear interp
              kG0=1
            ELSE !kZbot=1
              !Start integration at Z_UV(1) and truncate G-values
              Ubot=SQRT(u_p(i,j,1)*u_p(i,j,1)                           &
                       +v_p(i,j,1)*v_p(i,j,1))
              Tbot=THv(i,j,1)
              tmp=(REAL(NGz)/0.9)*((z_uv(i,j,1)/Ztop)-0.1)
              kG0=MAX(INT(tmp)+1,1)
              fG0=MAX(REAL(kG0)-tmp,0.0)
            END IF

            !Estimate GX values from lookup table
            HonL=Ztop*invMOsfce
            IF (HonL <  HLtab(1)) HonL=HLtab(1)
            GK=2
            DO WHILE ((HonL >  HLtab(GK)).AND.(GK <  gn))
              GK=GK+1
            END DO
            GMtab1=0.0
            GMtab2=0.0
            GHtab1=0.0
            GHtab2=0.0
            IF ((kG0 >  1).AND.(fG0 >  0.0)) THEN
              GMtab1=GMtab1+fG0*gmsav(GK-1,kG0-1)
              GMtab2=GMtab2+fG0*gmsav(GK,  kG0-1)
              GHtab1=GHtab1+fG0*GHsav(GK-1,kG0-1)
              GHtab2=GHtab2+fG0*GHsav(GK,  kG0-1)
            END IF
            DO k = kG0, NGz
              GMtab1=GMtab1+gmsav(GK-1,k)
              GMtab2=GMtab2+gmsav(GK,k)
              GHtab1=GHtab1+GHsav(GK-1,k)
              GHtab2=GHtab2+GHsav(GK,k)
            END DO
            h1=HLtab(GK-1)
            h2=HLtab(GK)
            g1=GMtab1
            g2=GMtab2
            Gs=(g2-g1)/(h2-h1)
            gm=Gs*(HonL-h2)+g2
            g1=GHtab1
            g2=GHtab2
            Gs=(g2-g1)/(h2-h1)
            gh=Gs*(HonL-h2)+g2
            Rilim=gh/(vkman*gm*gm)

            !Calculate Rtest
            u2=SQRT(u_p(i,j,kZtop)*u_p(i,j,kZtop)                       &
                   +v_p(i,j,kZtop)*v_p(i,j,kZtop))
            t2=THv(i,j,kZtop)
            Ujmp=u2-Ubot
            Tjmp=t2-Tbot
            Rtest=RtestMin-0.1e-10
            IF (ABS(Ujmp) >= 0.1e-10) THEN
              Rib=Ztop*(g/THv(i,j,2))*Tjmp/(Ujmp*Ujmp)
              Rtest=Rilim-Rib
            END IF
            IF (Rtest <= RtestMin) THEN
              !Ztop is at or above the actual SBL top
              topbl(i,j)=.TRUE.
              zh_esbl(i,j)=Ztop+1.0
              ntml_esbl(i,j)=kZtop
              !If we are above level kMINh,
              !estimate ZH_ESBL using subgrid scheme (below)
              IF (kZtop >  kMINh) subgrid=.TRUE.
            END IF
          END DO

          IF (subgrid .AND. sg_enabled) THEN !perform subgrid diagnosis
            dz100=(Ztop-Zprev)/100.0
            Rtest=RtestMin+1.0 !Ensures initial entry to WHILE loop
            Ztop=Zprev
            Zupper=zh_esbl(i,j)
            DO WHILE( (Rtest >  RtestMin).AND.                          &
                     (Ztop <= Zupper-1.0-dz100) )
              Zprev=Ztop
              Ztop=Ztop+dz100
              Zbot=0.1*Ztop !Height of surface layer top
              kZbot=0
              zz=0.0
              DO WHILE (zz <  Zbot) !find kZbot (UV-level above Zbot)
                kZbot=kZbot+1
                zz=z_uv(i,j,kZbot)
              END DO !find kZbot
              IF (kZbot >  1) THEN
                !Interpolation of U to Zbot
                k=kZbot-1
                u1=SQRT(u_p(i,j,k)*u_p(i,j,k)+v_p(i,j,k)*v_p(i,j,k))
                u2=SQRT(u_p(i,j,k+1)*u_p(i,j,k+1)                       &
                       +v_p(i,j,k+1)*v_p(i,j,k+1))
                t1=THv(i,j,k)
                t2=THv(i,j,k+1)
                dz=z_uv(i,j,k+1)-z_uv(i,j,k)
                Ubot=u1+(u2-u1)*(Zbot-z_uv(i,j,k))/dz !Linear interp
                Tbot=t1+(t2-t1)*(Zbot-z_uv(i,j,k))/dz !Linear interp
                kG0=1
              ELSE !kZbot=1
                !Start integration at Z_UV(1) and truncate G-values
                Ubot=SQRT(u_p(i,j,1)*u_p(i,j,1)                         &
                         +v_p(i,j,1)*v_p(i,j,1))
                Tbot=THv(i,j,1)
                tmp=(REAL(NGz)/0.9)*((z_uv(i,j,1)/Ztop)-0.1)
                kG0=MAX(INT(tmp)+1,1)
                fG0=MAX(REAL(kG0)-tmp,0.0)
              END IF

              !Estimate GX values from lookup table
              HonL=Ztop*invMOsfce
              IF (HonL <  HLtab(1)) HonL=HLtab(1)
              GK=2
              DO WHILE ((HonL >  HLtab(GK)).AND.(GK <  gn))
                GK=GK+1
              END DO
              GMtab1=0.0
              GMtab2=0.0
              GHtab1=0.0
              GHtab2=0.0
              IF ((kG0 >  1).AND.(fG0 >  0.0)) THEN
                GMtab1=GMtab1+fG0*gmsav(GK-1,kG0-1)
                GMtab2=GMtab2+fG0*gmsav(GK,  kG0-1)
                GHtab1=GHtab1+fG0*GHsav(GK-1,kG0-1)
                GHtab2=GHtab2+fG0*GHsav(GK,  kG0-1)
              END IF
              DO k = kG0, NGz
                GMtab1=GMtab1+gmsav(GK-1,k)
                GMtab2=GMtab2+gmsav(GK,k)
                GHtab1=GHtab1+GHsav(GK-1,k)
                GHtab2=GHtab2+GHsav(GK,k)
              END DO
              h1=HLtab(GK-1)
              h2=HLtab(GK)
              g1=GMtab1
              g2=GMtab2
              Gs=(g2-g1)/(h2-h1)
              gm=Gs*(HonL-h2)+g2
              g1=GHtab1
              g2=GHtab2
              Gs=(g2-g1)/(h2-h1)
              gh=Gs*(HonL-h2)+g2
              Rilim=gh/(vkman*gm*gm)

              !Calculate Rtest
              k=kZtop-1
              u1=SQRT(u_p(i,j,k)*u_p(i,j,k)                             &
                     +v_p(i,j,k)*v_p(i,j,k))
              u2=SQRT(u_p(i,j,k+1)*u_p(i,j,k+1)                         &
                     +v_p(i,j,k+1)*v_p(i,j,k+1))
              t1=THv(i,j,k)
              t2=THv(i,j,k+1)
              dz=z_uv(i,j,k+1)-z_uv(i,j,k)
              slope=(t2-t1)/dz
              t2=t1+slope*(Ztop-z_uv(i,j,k))
              slope=(u2-u1)/dz
              u2=u1+slope*(Ztop-z_uv(i,j,k))
              Ujmp=u2-Ubot
              Tjmp=t2-Tbot
              Rtest=RtestMin-0.1e-10
              IF (ABS(Ujmp) >= 0.1e-10) THEN
                Rib=Ztop*(g/THv(i,j,2))*Tjmp/(Ujmp*Ujmp)
                Rtest=Rilim-Rib
              END IF
              IF (Rtest <= RtestMin) zh_esbl(i,j)=Ztop
            END DO !DO WHILE
            IF (Rtest <= RtestMin) ntml_esbl(i,j)=kZtop-1
          END IF !end subgrid ZH diagnosis

        END IF !only for stable/neutral cases
      END DO !I
    END DO !j

!-----------------------------------------------------------------------
! 3.3 Prescribe local equilibrium fluxes and M-O length;
!     calculate Zhat and invoke equilibrium model for SBL;
!     convert output to diffusivities*density: RHO*Kx(K).
!-----------------------------------------------------------------------
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (fb_surf(i,j) <= 0.0) THEN !only for stable/neutral cases

          DO k = 2, bl_levels!Set Kx back to zero
            rhokm(i,j,k) = 0.0
            rhokh(i,j,k) = 0.0
          END DO
          UWsfce=v_s(i,j)*v_s(i,j)
          WTsfce=fb_surf(i,j) !Actually <w*THv>*(g/THv) at sfce
          invMOsfce=invMOsurf(i,j)
          hh=zh_esbl(i,j)
          DO k = 2, ntml_esbl(i,j)+1
            zz=z_tq(i,j,k-1)
            IF (zz <  hh) THEN !below ZH

              tmp=MAX(1.0-zz/hh,1.0e-3)
              USequil=SQRT(UWsfce*(tmp**Muw_SBL))
              WTequil=WTsfce*(tmp**Mwt_SBL) !Actually <w*THv>*(g/THv)
              invMOequil=invMOsfce/(tmp**(1.5*Muw_SBL-Mwt_SBL))
              Zhat   =zz*invMOequil
! DEPENDS ON: sblequil
              CALL sblequil(Zhat,lhat,pkm,pkh,pkuu,pktt,pkut,pktu,      &
                      PHIe,PHIw,RiSBL,ce,rpow,cb,cn,iERRSBL)
              ll =lhat*vkman*zz
              km =pkm *ll*USequil
              kh =pkh *ll*USequil
                ! Note "RHO" here is always wet density (RHO_TQ)
              rhokm(i,j,k) = rho(i,j,k-1) * km
              IF (lq_mix_bl) THEN
                ! Note "RHO" here is always wet density (RHO_TQ) so
                ! save multiplication of RHOKH to after interpolation
                rhokh(i,j,k) = kh
              ELSE
                rhokh(i,j,k) = rho(i,j,k-1) * kh
              END IF
              dbdu=-1.0
              IF ((dvdzm(i,j,k) >  1.e-10).AND.(dbdz(i,j,k) >  1.e-10)) &
                 dbdu=dbdz(i,j,k)/dvdzm(i,j,k)
              IF (L_SBLco.AND.(dbdu >  0.0)) THEN
                kuu=pkuu*ll*USequil
                ktt=pktt*ll*USequil
                kut=pkut*ll*ll !Actually KUT/(g/THv)
                tmp=WTequil/(USequil*USequil)
                GonTv=g/(  ( tl(i,j,k-1)+ grcp * zz      )              &
                          *( 1.0 + c_virtual*qw(i,j,k-1) ) )
                ktu=-pktu*ll*ll*tmp*tmp*GonTv*GonTv
                                                !Actually KTU*(g/THv)
                !Factor of g/THv in DBDU cancels in equations below,
                !due to factors of 1/(g/THv) in KUT, and g/THv in KTU
                ! Note "RHO" here is always wet density (RHO_TQ)
                rhokm(i,j,k) = rho(i,j,k-1)*(kuu+kut*dbdu)
                IF (lq_mix_bl) THEN
                ! Note "RHO" here is always wet density (RHO_TQ) so
                ! save multiplication of RHOKH to after interpolation
                  rhokh(i,j,k) = (ktu/dbdu+ktt)
                ELSE
                  rhokh(i,j,k) = rho(i,j,k-1)*(ktu/dbdu+ktt)
                END IF
              END IF

            END IF !below ZH
          END DO !K

        END IF !only for stable/neutral cases
      END DO !I
    END DO !j

!-----------------------------------------------------------------------
! 3.4 Set values of ZH_LOCAL and NTML_LOCAL.
!-----------------------------------------------------------------------
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (fb_surf(i,j) <= 0.0) THEN !only for stable/neutral cases
          zh_local(i,j)=zh_esbl(i,j)
          k=ntml_esbl(i,j)
          ntml_local(i,j)=k
          IF (z_tq(i,j,k) >  zh_esbl(i,j)) ntml_local(i,j)=k-1
        END IF !only for stable/neutral cases
      END DO !I
    END DO !j

  END IF !Equilibrium SBL model

!-----------------------------------------------------------------------
! Finish up
!-----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('EX_COEF',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE ex_coef
