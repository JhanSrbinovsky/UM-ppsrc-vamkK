! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************
!  Module mym_option_mod----------------------------------------------
!
!  Purpose: To define options and symbols in the Mellor-Yamada model
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
MODULE mym_option_mod

  USE  missing_data_mod, ONLY: rmdi, imdi
  USE control_max_sizes, ONLY: max_bl_levels
  IMPLICIT NONE

!=======================================================================
! Symbols of switches
!=======================================================================
! Symbols to choose TKE schemes (for BDY_TKE)
!
  INTEGER, PARAMETER ::                                                 &
     NONE = 0,                                                          &
            ! not use TKE schemes
     deardorff = 1,                                                     &
            ! the first order scheme
     mymodel25 = 2,                                                     &
            ! the improved Mellor-Yamada level 2.5 model
     mymodel3  = 3
            ! the improved Mellor-Yamada level 3 model

! Symbols for the option to choose  the mixing length
! in the first order model (for tke_dlen)
  INTEGER, PARAMETER ::                                                 &
     my_length = 1,                                                     &
            ! use the same mixing length in the MY model
     ddf_length = 2,                                                    &
            ! original mixing length suggested by Deardorff(1980)
     non_local_like_length = 3
            ! with correction based on Sun and Chang (1986)

! Symbols for the option to choose the gradient function
! to calculate production terms at the lowest level.
! (for MY_lowest_pd_surf)
  INTEGER, PARAMETER ::                                                 &
     no_pd_surf = 0,                                                    &
            ! No use of surface fluxes
     businger = 1,                                                      &
            ! use the gradient function by Businger
     bh1991 = 2
            ! use the gradient function by Beljaars and Holtslag

!=======================================================================
! MY model options
!=======================================================================
  INTEGER ::                                                            &
     bdy_tke = imdi,                                                    &
            ! switch to choose the TKE schemes, was NONE
     tke_dlen = imdi,                                                   &
            ! was "my_length" 
            ! switch to choose mixing length in the first order model
     my_lowest_pd_surf = no_pd_surf,                                    &
            ! switch related to production terms at the lowest levels
     tke_levels = imdi,                                                 &
            ! was -1
            ! maximum level to predict the prognostic variables in the
            ! TKE scheme.
     shcu_levels = imdi,                                                &
            ! was -1
            ! maximum level to evaluate the non-gradient buoyancy flux
     high_order_scheme_adv_turb = 1,                                    &
            ! choice of the higher order SL advection schemes
            ! for the prognostic variables in the TKE scheme.
     monotone_scheme_adv_turb = 0
            ! choice of the monotone SL advection schemes
            ! for the prognostic variables in the TKE scheme.

  REAL ::                                                               &
     my_z_extra_fact = 0.5,                                             &
            ! IF L_MY_PROD_ADJ == .TRUE., the extra level is assigned at
            ! Z_TQ(:,:,1) * MY_Z_EXTRA_FACT above the surface.
     tke_cm_mx = rmdi,                                                  &
            ! was 0.1                                                   
            ! A proportional coef CM below the top of mixed layer
            ! K = CM * SQRT(E) * L
     tke_cm_fa = rmdi,                                                  &
            ! was 0.1
            ! A proportional coef CM above the top of mixed layer
            ! K = CM * SQRT(E) * L
     my_ini_dbdz_min = rmdi,                                            &
            ! was 1.0e-5
            ! the lower limit for dbdz in the initialization to avoid
            ! to diagnose huge initial values.
     my_z_limit_elb = rmdi,                                             &
            ! was 1.0e10
            ! IF Z_TQ > MY_z_limit_elb, elb is limited less than
            ! vertical grid spacing.
     wb_ng_max = rmdi  ! was 0.05
            ! the maximum limit for the non-gradient buoyancy flux

  REAL :: my_prod_adj_fact(1:max_bl_levels) = rmdi
            ! was 0.225
            ! Factor in production term adjustment related to diffusion.

  LOGICAL ::                                                            &
     l_tke_dlen_blackadar = .FALSE.,                                    &
            ! use the correction to the mixing length by Blackadar
            ! (valid only in the first order model)
     l_my3_improved_closure = .TRUE.,                                   &
            ! use the improved closure constants in the MY model
            ! by Nakanishi and Niino.
            ! If False, use the original value from Mellor-Yamada(1982)
     l_my_condense = .FALSE.,                                           &
            ! was .TRUE.  
            ! If TRUE, buoyancy parameters appearing in TKE production
            ! by buoyancy are evaluated with predicted variants
            ! assuming that fluctuation of heat and moisture can be
            ! described by the bi-normal probability distribution
            ! function. If false, buoyancy parameters calculated in the
            ! large scale clouds scheme (e.g. Smith, PC2...), are
            ! employed  to evaluate the buoyancy flux.
     l_shcu_buoy = .FALSE.,                                             &
            ! was .TRUE.  
            ! If TRUE, non-gradient buoyancy flux based on Lock and
            ! Mailhot(2006) is added.
     l_adv_turb_field = .TRUE.,                                         &
            ! A switch to turn on advection of the prognostic variables
            ! in the TKE scheme (E_TRB, TSQ, QSQ, COV).
     l_high_adv_turb = .TRUE.,                                          &
            ! A switch for higher order schemes of the SL advection of
            ! the prognostic variables in the TKE scheme.
     l_mono_adv_turb = .FALSE.,                                         &
            ! A switch for monotone schemes of the SL advection of
            ! the prognostic variables in the TKE scheme.
     l_conserv_adv_turb = .FALSE.,                                      &
            ! A switch for the conservative SL advection of
            ! the prognostic variables in the TKE scheme.
     l_my_extra_level = .FALSE.,                                        &
            ! If TRUE, the extra level below the lowest level in the
            ! atmosphere is generated for the prognostic variables in
            ! the TKE scheme.
     l_my_prod_adj = .FALSE.,                                           &
            ! was .TRUE.
            ! If TRUE, production terms of covariance brought
            ! by the counter gradient terms are adjusted
            ! so that stability in integration can be secured.
     l_local_above_tkelvs = .FALSE.,                                    &
            ! If TRUE and the parameter "tke_levels" is less
            ! than bl_levels, diffusion coefficients between
            ! tke_levels + 1 and bl_levels are evaluated with the
            ! stability function (i.e. local scheme)
     l_my_initialize = .FALSE.,                                         &
            ! If TRUE, the prognostic variables in the TKE schemes are
            ! initialized assuming balance between production and
            ! dissipation.
            ! If missing values are set to the prognostic variables
            ! through the reconfiguration, the initialization will 
            ! be automatically conducted.
     l_my_ini_zero = .FALSE.,                                           &
            ! In the case that l_my_initialize == .TRUE. or 
            ! missing values are set to the prognostic variables
            ! through the reconfiguration, they are initialized to zeros
            ! if it is true.
     l_print_max_tke = .FALSE.,                                         &
            ! If TRUE, the maximum values of the prognostic variables
            ! are printed.
     l_my_lowest_pd_surf_tqc = .FALSE.
            ! If TRUE, the production terms for the covariaces are
            ! evaluated with surface fluxes. It is valid only when
            ! MY_lowest_pd_surf > 0.

END MODULE mym_option_mod
