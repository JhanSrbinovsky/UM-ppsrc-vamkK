! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates the required EOT climate diagnostics, and fills STASHWORK.
!
! Subroutine Interface:
MODULE eot_diag_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE Eot_diag(                                              &
! Primary data: in
     &  pstar,p,rho,u,v,w,theta,q,qcl,qcf                               &
     &  ,p_theta_levels                                                 &
     &  ,exner_rho_levels,exner_theta_levels                            &
     &  ,energy_corr_now                                                &
     &  ,inc_u, inc_v, inc_w, inc_t                                     &
     &  ,inc_q, inc_qcl, inc_qcf, inc_rho                               &
        ,inc_qrain, inc_qgraup, inc_qcf2                                &
        ,wet_to_dry_n                                                   &
! Grid sizes and definition: in
     &  ,rows,n_rows,row_length,model_levels,wet_levels,bl_levels       &
     &  ,theta_field_size,u_field_size,v_field_size                     &
     &  ,sin_theta_latitude,cos_theta_latitude                          &
     &  ,sin_v_latitude                                                 &
     &  ,sin_theta_longitude,cos_theta_longitude                        &
     &  ,Model_domain                                                   &
     &  ,delta_lambda                                                   &
                                ! grid longitude spacing in radians
     &  ,delta_phi                                                      &
                                ! grid latitude  spacing in radians
     &, global_row_length                                               &
! Grid coordinates: in
     &  ,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole    &
! Pressure levels for output arrays: in
     &  ,u_press,v_press,w_press                                        &
     &  ,T_PRESS,Q_PRESS,RH_PRESS                                       &
     &  ,Z_PRESS,OM_PRESS,HEAVY_PRESS,VSTARBAR_PRESS,WSTARBAR_PRESS     &
     &  ,FY_PRESS,FZ_PRESS,DIVF_PRESS,zvptp_press,zupvp_press           &
! Flags to request each diagnostic output field: in
     &  ,qsf,  qv_m, qw_m,  qt_m, qq_m, qz_m, qke_m, qom_m              &
     &  ,qu_mm,qv_mm,qw_mm ,qt_mm,qq_mm,qz_mm,qke_mm                    &
     &  ,qteot_m,qwbig_m,qwbig2_m,qRH_m,qdry_mass_m                     &
     &  ,qt_inc,qq_inc,qqcl_inc,qqcf_inc                                &
     &  ,qu_inc,qv_inc,qw_inc,qrho_inc                                  &
        ,qqrain_inc, qqgraup_inc, qqqcf2_inc                            &
     &  ,qt2_inc,qq2_inc,qqcl2_inc,qqcf2_inc,qu2_inc,qv2_inc            &
     &  ,qw2_inc,qrho2_inc                                              &
     &  ,qu_p,qv_p,qw_p                                                 &
     &  ,qt_p,qq_p,qrh_p                                                &
     &  ,qz_p,qom_p                                                     &
     &  ,qheavy_p,qtv_p,qtvom_p                                         &
     &  ,qvstarbar_p,qwstarbar_p,qepy_p,qepz_p,qdivep_p                 &
     &  ,qzvptp_p,qzupvp_p                                         &
     &  ,qtot_ke,qtot_ke_w,qtcm                                         &
     &  ,qtcmq,qtcqcl,qtcqcf                                            &
     &  ,qtmf_u,qtmf_v,qtmf_w,qmtorque                                  &
     &  ,qam_m1,qam_m2,qam_m3                                           &
     &  ,qam_w1,qam_w2,qam_w3                                           &
     &  ,qpstar,qpstar_uv                                               &
     &  ,qencorr,qtot_cvt,qtot_gr                                       &
     &  ,qcol_ugz, qcol_vgz, qcol_wgz, qcol_uT,  qcol_vt,  qcol_wt      &
     &  ,qcol_uq,  qcol_vq,  qcol_wq , qcol_uv,  qcol_uw,  qcol_vw      &
     &  ,qcol_uKe, qcol_vKe, qcol_wKe, qcol_u,   qcol_v,   qcol_w       &
     &  ,qm_spd_ew,qm_spd_ns                                            &
     &  ,qtrop_p,qtrop_t,qtrop_z,qtrop_icao                             &
        ,qcol_sat1                                                      &
! Diagnostics lengths: in
     &  ,u_m_levs,v_m_levs,w_m_levs,t_m_levs,q_m_levs,z_m_levs          &
     &  ,ke_m_levs,om_m_levs,u_mm_levs,v_mm_levs,w_mm_levs,t_mm_levs    &
     &  ,q_mm_levs,z_mm_levs,ke_mm_levs                                 &
     &  ,teot_m_levs,wbig_m_levs,wbig2_m_levs,Rh_m_levs,dry_mass_m_levs &
     &  ,t_inc_levs,q_inc_levs,qcl_inc_levs,qcf_inc_levs                &
        ,qrain_inc_levs, qgraup_inc_levs, qqcf2_inc_levs                &
     &  ,u_inc_levs,v_inc_levs                                          &
     &  ,w_inc_levs,rho_inc_levs                                        &
     &  ,t2_inc_levs,q2_inc_levs,qcl2_inc_levs,qcf2_inc_levs            &
     &  ,u2_inc_levs,v2_inc_levs,w2_inc_levs,rho2_inc_levs              &
     &  ,u_p_levs,v_p_levs,w_p_levs,t_p_levs,q_p_levs,rh_p_levs         &
     &  ,z_p_levs,om_p_levs,heavy_p_levs                                &
     &  ,tv_p_levs,tvom_p_levs                                          &
     &  ,vstarbar_p_levs, wstarbar_p_levs                               &
     &  , Fy_p_levs, Fz_p_levs, divF_p_levs                             &
     &  , zvptp_p_levs,zupvp_p_levs                                &
     &  ,prod_p_levs                                                    &
                     ! pressure levels required for each product
     &  ,npress_diags,nmodel_diags                                      &
! Tropopause index bounds
     &  ,min_trop_level,max_trop_level                                  &
! Model levels indicies: in
     &  ,u_m_list,v_m_list,w_m_list,t_m_list,q_m_list,z_m_list,ke_m_list&
     &  ,om_m_list,u_mm_list,v_mm_list,w_mm_list,t_mm_list,q_mm_list    &
     &  ,z_mm_list,ke_mm_list,teot_m_list,wbig_m_list,wbig2_m_list      &
     &  ,RH_m_list,dry_mass_m_list                                      &
     &  ,t_inc_list,q_inc_list,qcl_inc_list,qcf_inc_list                &
        ,qrain_inc_list, qgraup_inc_list, qqcf2_inc_list                &
     &  ,u_inc_list,v_inc_list                                          &
     &  ,w_inc_list,rho_inc_list                                        &
     &  ,t2_inc_list,q2_inc_list,qcl2_inc_list,qcf2_inc_list            &
     &  ,u2_inc_list,v2_inc_list,w2_inc_list,rho2_inc_list              &
     &  ,prod_m_list                                                    &
                     ! index of model levels required for product
! Product levels: in
     &  ,prod_ind                                                       &
! Diagnostic arrays: out
     &  ,si,stashwork                                                   &
     &  ,u_m,v_m,w_m ,t_m,q_m,z_m,ke_m,om_m                             &
     &  ,u_mm,v_mm,w_mm,t_mm,q_mm,z_mm,ke_mm                            &
     &  ,teot_m ,wbig_m,wbig2_m,RH_m,dry_mass_m                         &
     &  ,t_inc,q_inc,qcl_inc,qcf_inc,u_inc,v_inc,w_inc,rho_inc          &
        ,qrain_inc, qgraup_inc, qcft2_inc                               &
     &  ,t2_inc,q2_inc,qcl2_inc,qcf2_inc,u2_inc,v2_inc,w2_inc,rho2_inc  &
     &  ,u_p,v_p,w_p,t_p,q_p,rh_p ,z_p,om_p,heavy_p                     &
     &  ,tv_p,tvom_p                                                    &
     &  ,vstarbar, wstarbar, Fy, Fz, divF, zvptp, zupvp            &
     &  ,tot_ke,tot_ke_w,tcm                                            &
     &  ,tcmq,tcqcl,tcqcf                                               &
     &  ,tmf_u,tmf_v,tmf_w,mtorque                                      &
     &  ,am_m1,am_m2,am_m3                                              &
     &  ,am_w1,am_w2,am_w3                                              &
     &  ,o_pstar,o_pstar_uv                                             &
     &  ,encorr,tot_cvt,tot_gr                                          &
     &  ,col_ugz, col_vgz, col_wgz, col_uT,  col_vt,  col_wt            &
     &  ,col_uq,  col_vq,  col_wq , col_uv,  col_uw,  col_vw            &
     &  ,col_uKe, col_vKe, col_wKe, col_u,   col_v,   col_w             &
     &  ,m_spd_ew,m_spd_ns                                              &
     &  ,trop_p,trop_t,trop_z,trop_icao                                 &
        ,col_sat1                                                       &
     &  )
  
      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE earth_constants_mod, ONLY: g, earth_radius, omega

      USE level_heights_mod, ONLY:                                      &
                      r_theta_levels, r_rho_levels

      USE atm_fields_bounds_mod, ONLY:                                  &
          udims, vdims, wdims, tdims, pdims, qdims,                     &
          udims_s, vdims_s, wdims_s, tdims_s, pdims_s, qdims_l

      USE atmos_constants_mod, ONLY:                                    &
          kappa, p_zero, cp, r, c_virtual

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      USE domain_params
      USE interpor_mod, ONLY: interp_order_linear
      USE vert_eng_massq_mod, ONLY: vert_eng_massq

      USE Field_Types
      USE calc_div_ep_flux_mod, ONLY: calc_div_ep_flux
      USE pc_to_pb_mod, ONLY: pc_to_pb
      USE p_to_t_mod, ONLY: p_to_t
      USE p_to_u_mod, ONLY: p_to_u
      USE p_to_v_mod, ONLY: p_to_v
      USE t_vert_interp_to_p_mod, ONLY: t_vert_interp_to_p
      USE u_to_p_mod, ONLY: u_to_p
      USE uc_to_ub_mod, ONLY: uc_to_ub
      USE v_to_p_mod, ONLY: v_to_p
      USE vert_h_onto_p_mod, ONLY: vert_h_onto_p
      USE vert_interp2_mod, ONLY: vert_interp2
      USE vc_to_vb_mod, ONLY: vc_to_vb

      IMPLICIT NONE
!
! Description:
!   Calculate eot-related diagnostics - held in STASH section 30 -
!   The stash items are split into centuries, as follows:
!
!   000s   Standard Model level fields and products
!   100s   Other Model level fields
!   200s   Standard Pressure level fields and products
!   300s   Other Pressure level fields
!   400s   Surface fields and Column integrals
!
!   Diagnostics on model rho points
!    1  u component of wind
!    2  v component of wind
!    3  w component of wind
!    4  T temperature
!    5  q specific humidity
!    6  z height
!    7  KE = 0.5(u*u+v*v+w*w)
!   Diagnostics on model theta points
!    8  om omega
!   Mass weighted products on model rho levels
!   000+10*(1st code)+(2nd code) ie
!    11 u*u
!    12 u*v
!    17 u*ke
!   mass weighted fields on rho levels
!   where (field * rho *r*r*dr/a*a) a - radius of earth
!   101  u component of wind
!   102  v component of wind
!   103  w component of wind
!   104  T temperature
!   105  q specific humidity
!   106  z height
!   107  KE = 0.5(u*u+v*v+w*w)
!
!   Increments across whole model timestep
!   181  - temperature
!   182  - q
!   183  - qcl
!   184  - qcf
!   185  - U
!   186  - V
!   187  - w
!   188  - density
!   189  - prognostic rain
!   190  - graupel
!   191  - qcf2 - second ice crystal type.
!
!   Diagnostics on pressure Pressure surfaces.
!   STASH item
!   201 u component of wind on pressure surfaces
!   202 v component of wind on pressure surfaces
!   203 w component of wind on pressure surfaces
!   204 t temperature on pressure surfaces
!   205 q specific humidity on pressure surfaces
!   206 rh relative humidity on pressure surfaces
!   207 z geopotential on pressure surfaces
!   208 om omega on pressure surfaces
!
!   Products are also produced, the stashcodes are
!   200+10*(1st code)+(2nd code) ie
!   211 UU
!   212 UV
!   224 VT
!
!
!  Further user and technical documentation can be found under:
!  http://www-hc/~hadsu/
!
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Code 201-208 are processed sepertely, product codes
!   are calculated in a loop over single pressure surface fields
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Climate Diagnostics
!
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables
! Model_domain meaningful names
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

! Subroutine arguments
!   Scalar arguments with intent(in):
! Grid sizes:
      INTEGER                                                           &
     &  rows,n_rows,row_length,model_levels,wet_levels,bl_levels        &
     &  ,theta_field_size,u_field_size,v_field_size                     &
     &  ,Model_domain                                                   &
     &  ,global_row_length ! global row length
! Grid coordinates: in
      REAL                                                              &
     & ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole
! Diagnostics lengths: in
      INTEGER, INTENT(IN) ::                                            &
     & u_m_levs                                                         &
                      ! No of levels for output of U_M
     &,v_m_levs                                                         &
                      ! No of levels for output of V_M
     &,w_m_levs                                                         &
                      ! No of levels for output of W_M
     &,t_m_levs                                                         &
                      ! No of levels for output of T_M
     &,q_m_levs                                                         &
                      ! No of levels for output of Q_M
     &,z_m_levs                                                         &
                      ! No of levels for output of z_M
     &,ke_m_levs                                                        &
                       ! No of levels for output of ke_M
     &,om_m_levs                                                        &
                       ! No of levels for output of om_M
     &,wbig_m_levs                                                      &
                      ! No of levels for output of wbig_m
     &,wbig2_m_levs                                                     &
                       ! No of levels for output of wbig2_m
     &,Rh_m_levs                                                        &
                      ! No of levels for output of RH
     &,dry_mass_m_levs                                                  &
                       ! No of levels for output of dry mass weighting
     &,u_mm_levs                                                        &
                       ! No of levels for output of U_MM
     &,v_mm_levs                                                        &
                       ! No of levels for output of V_MM
     &,w_mm_levs                                                        &
                       ! No of levels for output of W_MM
     &,t_mm_levs                                                        &
                       ! No of levels for output of T_MM
     &,q_mm_levs                                                        &
                       ! No of levels for output of Q_MM
     &,z_mm_levs                                                        &
                       ! No of levels for output of z_MM
     &,ke_mm_levs                                                       &
                        ! No of levels for output of KE_MM
     &,teot_m_levs                                                      &
                         ! No of levels for output of TEOT_M
     &,u_inc_levs                                                       &
                        ! No of levels for output of U_INC
     &,v_inc_levs                                                       &
                        ! No of levels for output of V_INC
     &,w_inc_levs                                                       &
                        ! No of levels for output of W_INC
     &,t_inc_levs                                                       &
                        ! No of levels for output of T_INC
     &,q_inc_levs                                                       &
                        ! No of levels for output of Q_INC
     &,qcl_inc_levs                                                     &
                          ! No of levels for output of Q_INC
     &,qcf_inc_levs                                                     &
                          ! No of levels for output of Q_INC
      ,qrain_inc_levs   & ! No of levels for output of qrain_inc
      ,qgraup_inc_levs  & ! No of levels for output of qgraup_inc
      ,qqcf2_inc_levs   & ! No of levels for output of qcf2_inc
     &,rho_inc_levs                                                     &
                          ! No of levels for output of Q_INC
     &,u2_inc_levs                                                      &
                         ! No of levels for output of U2_INC
     &,v2_inc_levs                                                      &
                         ! No of levels for output of V2_INC
     &,w2_inc_levs                                                      &
                         ! No of levels for output of W2_INC
     &,t2_inc_levs                                                      &
                         ! No of levels for output of T2_INC
     &,q2_inc_levs                                                      &
                         ! No of levels for output of Q2_INC
     &,qcl2_inc_levs                                                    &
                           ! No of levels for output of Qcl2_INC
     &,qcf2_inc_levs                                                    &
                           ! No of levels for output of Qcf2_INC
     &,rho2_inc_levs                                                    &
                           ! No of levels for output of rho2_INC
     &,u_p_levs                                                         &
                      ! No of levels for output of U_P
     &,v_p_levs                                                         &
                      ! No of levels for output of V_P
     &,w_p_levs                                                         &
                      ! No of levels for output of W_P
     &,t_p_levs                                                         &
                      ! No of levels for output of T_P
     &,q_p_levs                                                         &
                      ! No of levels for output of Q_P
     &,rh_p_levs                                                        &
                       ! No of levels for output of RH_P
     &,z_p_levs                                                         &
                      ! No of levels for output of Z_P
     &,om_p_levs                                                        &
                       ! No of levels for output of OM_P
     &,heavy_p_levs                                                     &
                       ! No of levels for output of HEAVY_P
     &,tv_p_levs                                                        &
                        ! No of levels for output of TV_P
     &,tvom_p_levs     ! No of levels for output of TVOM_P             

! No of levels for residual circulation
      INTEGER, INTENT(IN) :: vstarbar_p_levs, wstarbar_p_levs,          &
                             Fy_p_levs, Fz_p_levs, divF_p_levs,         &
                             zvptp_p_levs, zupvp_p_levs


! Flags to request each diagnostic output field: IN
      LOGICAL, INTENT(IN) ::                                            &
     &   qv_m                                                           &
               ! Flag for V wind component on model levels
     &  ,qw_m                                                           &
                 ! Flag for W wind component on model level
     &  ,qt_m                                                           &
                 ! Flag for T on model level
     &  ,qq_m                                                           &
                 ! Flag for Q on model level
     &  ,qz_m                                                           &
                 ! Flag for z on model level
     &  ,qke_m                                                          &
                  ! Flag for ke on model level
     &  ,qom_m                                                          &
                  ! Flag for omega on model level
     &  ,qwbig_m                                                        &
                 ! Flag for wbig on model levels
     &  ,qwbig2_m                                                       &
                  ! Flag for wbig0.1 on model levels
     &  ,qRH_m                                                          &
                 ! Flag for RH on model levels
     &  ,qdry_mass_m                                                    &
                     ! Flag for dry mass weighting on model levels
     &  ,qu_mm                                                          &
                 ! Flag for U wind component on model levels
     &  ,qv_mm                                                          &
                ! Flag for V wind component on model levels
     &  ,qw_mm                                                          &
                  ! Flag for W wind component on model level
     &  ,qt_mm                                                          &
                  ! Flag for T on model level
     &  ,qq_mm                                                          &
                  ! Flag for Q on model level
     &  ,qz_mm                                                          &
                  ! Flag for z on model level
     &  ,qke_mm                                                         &
                   ! Flag for KE on model level
     &  ,qteot_m                                                        &
                    ! Flag for T at eot on model level
     &  ,qu_inc                                                         &
                  ! Flag for U wind component on model levels
     &  ,qv_inc                                                         &
                 ! Flag for V wind component on model levels
     &  ,qw_inc                                                         &
                   ! Flag for W wind component on model level
     &  ,qt_inc                                                         &
                   ! Flag for T on model level
     &  ,qq_inc                                                         &
                   ! Flag for Q on model level
     &  ,qqcl_inc                                                       &
                     ! Flag for QCL on model level
     &  ,qqcf_inc                                                       &
                     ! Flag for QCF on model level
        ,qqrain_inc  & ! Flag for qrain on model level
        ,qqgraup_inc & ! Flag for qgraup on model level
        ,qqqcf2_inc  & ! Flag for qcf2 on model level
     &  ,qrho_inc                                                       &
                     ! Flag for RHO on model level
     &  ,qu2_inc                                                        &
                    ! Flag for U2 wind component on model levels
     &  ,qv2_inc                                                        &
                    ! Flag for V2 wind component on model levels
     &  ,qw2_inc                                                        &
                    ! Flag for W2 wind component on model level
     &  ,qt2_inc                                                        &
                    ! Flag for T2 on model level
     &  ,qq2_inc                                                        &
                    ! Flag for Q2 on model level
     &  ,qqcl2_inc                                                      &
                      ! Flag for Qcl2 on model level
     &  ,qqcf2_inc                                                      &
                      ! Flag for Qcf2 on model level
     &  ,qrho2_inc                                                      &
                      ! Flag for RHO2 on model level
     &  ,qu_p                                                           &
                ! Flag for U wind component on pressure levels
     &  ,qv_p                                                           &
               ! Flag for V wind component on pressure levels
     &  ,qw_p                                                           &
                 ! Flag for W wind component on pressure level
     &  ,qt_p                                                           &
                 ! Flag for T on pressure level
     &  ,qvstarbar_p                                                    &
                 ! Flag for residual mean meridional circulation
     &  ,qwstarbar_p                                                    &
                 ! Flag for residual mean meridional circulation
     &  ,qepy_p                                                         &
                 ! Flag for Eliassen Palm flux (phi component)
     &  ,qepz_p                                                         &
                 ! Flag for Eliassen Palm flux (vertical component)
     &  ,qdivep_p                                                       &
                 ! Flag for divergence of Eliassen Palm flux
     &  ,qzvptp_p                                                       &
                 ! Flag for meridional heat flux
     &  ,qzupvp_p                                                       &
                 ! Flag for meridional momentum flux
     &  ,qq_p                                                           &
                 ! Flag for Q on pressure level
     &  ,qrh_p                                                          &
                 ! Flag for RH on pressure level
     &  ,qz_p                                                           &
                 ! Flag for Z on pressure level
     &  ,qom_p                                                          &
                 ! Flag for Omega on pressure level
     &  ,qheavy_p                                                       &
                   ! Flag for Heavyside fn on pressure level
     &  ,qtv_p                                                          &
                   ! Flag for virtual temp  on pressure level
     &  ,qtvom_p                                                        &
                   ! Flag for virtual temp*om on pressure level
     &  ,qsf(nitems)                                                    &
                     ! Flag for all diags
     &  ,qtot_ke                                                        &
     &  ,qtot_ke_w                                                      &
     &  ,qtot_cvt                                                       &
     &  ,qtot_gr                                                        &
     &  ,qtcm                                                           &
     &  ,qtcmq                                                          &
     &  ,qtcqcl                                                         &
     &  ,qtcqcf                                                         &
     &  ,qtmf_u                                                         &
     &  ,qtmf_v                                                         &
     &  ,qtmf_w                                                         &
     &  ,qmtorque                                                       &
     &  ,qm_spd_ew                                                      &
     &  ,qm_spd_ns                                                      &
     &  ,qam_m1                                                         &
     &  ,qam_m2                                                         &
     &  ,qam_m3                                                         &
     &  ,qam_w1                                                         &
     &  ,qam_w2                                                         &
     &  ,qam_w3                                                         &
     &  ,qpstar                                                         &
     &  ,qpstar_uv                                                      &
     &  ,qencorr                                                        &
     &  ,qtrop_p                                                        &
     &  ,qtrop_t                                                        &
     &  ,qtrop_z                                                        &
     &  ,qtrop_icao                                                     &
     &  ,qcol_ugz, qcol_vgz, qcol_wgz                                   &
     &  ,qcol_uT,  qcol_vt,  qcol_wt                                    &
     &  ,qcol_uq,  qcol_vq,  qcol_wq                                    &
     &  ,qcol_uv,  qcol_uw,  qcol_vw                                    &
     &  ,qcol_uKe, qcol_vKe, qcol_wKe                                   &
     &  ,qcol_u,   qcol_v,   qcol_w                                     &
        ,qcol_sat1

!   Array  arguments with intent(in):
! Primary data: IN
      REAL                                                              &
         u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end, &
           udims_s%k_start:udims_s%k_end)                               &
        ,v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end, &
           vdims_s%k_start:vdims_s%k_end)                               &
        ,w(wdims_s%i_start:wdims_s%i_end,wdims_s%j_start:wdims_s%j_end, &
           wdims_s%k_start:wdims_s%k_end)                               &
        ,theta(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)                           &
        ,q(qdims_l%i_start:qdims_l%i_end,qdims_l%j_start:qdims_l%j_end, &
           qdims_l%k_start:qdims_l%k_end)                               &
        ,qcl(qdims_l%i_start:qdims_l%i_end,                             &
             qdims_l%j_start:qdims_l%j_end,                             &
             qdims_l%k_start:qdims_l%k_end)                             &
        ,qcf(qdims_l%i_start:qdims_l%i_end,                             &
             qdims_l%j_start:qdims_l%j_end,                             &
             qdims_l%k_start:qdims_l%k_end)                             &
        ,rho(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end,                             &
             pdims_s%k_start:pdims_s%k_end)                             &
        ,p(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end, &
           pdims_s%k_start:pdims_s%k_end)                               &
        ,pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,p_theta_levels(tdims_s%i_start:tdims_s%i_end,                  &
                        tdims_s%j_start:tdims_s%j_end,                  &
                        tdims_s%k_start:tdims_s%k_end)                  &
        ,exner_rho_levels(pdims_s%i_start:pdims_s%i_end,                &
                          pdims_s%j_start:pdims_s%j_end,                &
                          pdims_s%k_start:pdims_s%k_end)                &
        ,exner_theta_levels(tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end)              &
        ,energy_corr_now                                                &
      , inc_u(udims_s%i_start:udims_s%i_end,                            &
              udims_s%j_start:udims_s%j_end,                            &
              udims_s%k_start:udims_s%k_end)                            &
      , inc_v(vdims_s%i_start:vdims_s%i_end,                            &
              vdims_s%j_start:vdims_s%j_end,                            &
              vdims_s%k_start:vdims_s%k_end)                            &
      , inc_w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,      &
              wdims%k_start:wdims%k_end)                                &
      , inc_t(tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end)                            &
      , inc_q(qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            &
              1:qdims_l%k_end)                                          &
      , inc_qcl(qdims_l%i_start:qdims_l%i_end,                          &
                qdims_l%j_start:qdims_l%j_end,                          &
                1:qdims_l%k_end)                                        &
      , inc_qcf(qdims_l%i_start:qdims_l%i_end,                          &
                qdims_l%j_start:qdims_l%j_end,                          &
                1:qdims_l%k_end)                                        &
      , inc_qrain(qdims_l%i_start:qdims_l%i_end,                        &
                  qdims_l%j_start:qdims_l%j_end,                        &
                  1:qdims_l%k_end)                                      &
      , inc_qgraup(qdims_l%i_start:qdims_l%i_end,                       &
                   qdims_l%j_start:qdims_l%j_end,                       &
                   1:qdims_l%k_end)                                     &
      , inc_qcf2(qdims_l%i_start:qdims_l%i_end,                         &
                 qdims_l%j_start:qdims_l%j_end,                         &
                 1:qdims_l%k_end)                                       &
      , inc_rho(pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,                          &
                pdims_s%k_start:pdims_s%k_end)                          &
      , wet_to_dry_n(pdims_s%i_start:pdims_s%i_end,                     &
                     pdims_s%j_start:pdims_s%j_end,                     &
                     pdims_s%k_start:pdims_s%k_end)                     &
        ,cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end)              &
                                                 ! cos(latitude)
        ,cos_theta_longitude(tdims%i_start:tdims%i_end,                 &
                             tdims%j_start:tdims%j_end)                 &
        ,sin_theta_longitude(tdims%i_start:tdims%i_end,                 &
                             tdims%j_start:tdims%j_end)                 &
        ,sin_theta_latitude(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end)                  &
        ,sin_v_latitude (vdims%i_start:vdims%i_end,                     &
                         vdims%j_start:vdims%j_end)                     &

! Pressure levels (units mb) for output arrays: IN
     &  ,u_press(u_p_levs)                                              &
                                ! for u wind
     &  ,v_press(v_p_levs)                                              &
                                ! for v wind
     &  ,w_press(w_p_levs)                                              &
                                ! for w wind
     &  ,t_press(t_p_levs)                                              &
                                ! for t
     &  ,q_press(q_p_levs)                                              &
                                ! for q
     &  ,rh_press(rh_p_levs)                                            &
                                ! for rh
     &  ,z_press(z_p_levs)                                              &
                                ! for z
     &  ,om_press(om_p_levs)                                            &
                                ! for om
     &  ,heavy_press(heavy_p_levs)                                      &
                                   ! for heavyside fn
     &  ,vstarbar_press(vstarbar_p_levs)                                &
                                ! Pressure for vstarbar
     &  ,wstarbar_press(wstarbar_p_levs)                                &
                                ! Pressure for wstarbar
     &  ,Fy_press(Fy_p_levs)                                            &
                                ! Pressure for meridional comp EP flux
     &  ,Fz_press(Fz_p_levs)                                            &
                                ! Pressure for vert comp EP flux
     &  ,divF_press(divF_p_levs)                                        &
                                ! Pressure for divergence EP flux
     &  ,zvptp_press(zvptp_p_levs)                                &
                                ! Pressure for meridional heat flux
     &  ,zupvp_press(zupvp_p_levs)                                  &
                                ! Pressure for meridional momentum flux
     &  ,Zlog(u_p_levs)                                            &
                                ! "log-pressure" coordinate
     &  ,delta_lambda                                                   &
                                ! grid longitude spacing in radians
     &  ,delta_phi              ! grid latitude  spacing in radians


      INTEGER, INTENT(in) ::                                            &
     &  npress_diags                                                    &
                                !# of diags on pressure levels
     &  ,nmodel_diags                                                   &
                                 !# of diags on model levels
     & ,prod_IND(NUM_STASH_LEVELS*2,npress_diags,npress_diags)          &
     & ,prod_p_levs(npress_diags,npress_diags)

! Tropopause index bounds.  Should be consistent with
! a statement of the form #include <typcona/typcona.h>.
! and also with Deck O3EXP1.84--87
      INTEGER, INTENT(IN) :: min_trop_level
!         Minimum permitted level of the tropopause
      INTEGER, INTENT(IN) :: max_trop_level
!         Maximum permitted level of the tropopause

      LOGICAL, INTENT(IN) :: u_m_list(pdims%k_start:pdims%k_end)        &
        ,v_m_list(pdims%k_start:pdims%k_end)                            &
        ,w_m_list(pdims%k_start:pdims%k_end)                            &
        ,wbig_m_list(wdims%k_start:wdims%k_end)                         &
        ,wbig2_m_list(wdims%k_start:wdims%k_end)                        &
        ,t_m_list(pdims%k_start:pdims%k_end)                            &
        ,rh_m_list(1:qdims%k_end)                                       &
        ,q_m_list(1:qdims%k_end)                                        &
        ,z_m_list(pdims%k_start:pdims%k_end)                            &
        ,ke_m_list(pdims%k_start:pdims%k_end)                           &
        ,om_m_list(1:tdims%k_end)                                       &
        ,u_mm_list(pdims%k_start:pdims%k_end)                           &
        ,v_mm_list(pdims%k_start:pdims%k_end)                           &
        ,w_mm_list(pdims%k_start:pdims%k_end)                           &
        ,t_mm_list(pdims%k_start:pdims%k_end)                           &
        ,q_mm_list(1:qdims%k_end)                                       &
        ,z_mm_list(pdims%k_start:pdims%k_end)                           &
        ,ke_mm_list(pdims%k_start:pdims%k_end)                          &
        ,teot_m_list(1:tdims%k_end)                                     &
        ,dry_mass_m_list(pdims%k_start:pdims%k_end)                     &
        ,u_inc_list(udims%k_start:udims%k_end)                          &
        ,v_inc_list(vdims%k_start:vdims%k_end)                          &
        ,w_inc_list(wdims%k_start:wdims%k_end)                          &
        ,t_inc_list(1:tdims%k_end)                                      &
        ,q_inc_list(1:qdims%k_end)                                      &
        ,qcl_inc_list(1:qdims%k_end)                                    &
        ,qcf_inc_list(1:qdims%k_end)                                    &
        ,qrain_inc_list(1:qdims%k_end)                                  &
        ,qgraup_inc_list(1:qdims%k_end)                                 &
        ,qqcf2_inc_list(1:qdims%k_end)                                  &
        ,rho_inc_list(pdims%k_start:pdims%k_end)                        &
        ,u2_inc_list(udims%k_start:udims%k_end)                         &
        ,v2_inc_list(vdims%k_start:vdims%k_end)                         &
        ,w2_inc_list(wdims%k_start:wdims%k_end)                         &
        ,t2_inc_list(1:tdims%k_end)                                     &
        ,q2_inc_list(1:qdims%k_end)                                     &
        ,qcl2_inc_list(1:qdims%k_end)                                   &
        ,qcf2_inc_list(1:qdims%k_end)                                   &
        ,rho2_inc_list(pdims%k_start:pdims%k_end)

      LOGICAL prod_m_list(pdims%k_start:pdims%k_end,                    &
                          nmodel_diags,nmodel_diags)


!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):
      INTEGER si(nitems)
      REAL stashwork(*)

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
! Diagnostic arrays: OUT
      REAL, INTENT(OUT), TARGET ::                                      & 
        u_p(udims%i_start:udims%i_end,                                  &
             vdims%j_start:vdims%j_end,u_p_levs)                        &
                                           ! u at selected pressures
        ,v_p(udims%i_start:udims%i_end,                                 &
             vdims%j_start:vdims%j_end,v_p_levs)                        &
                                           ! v at selected pressures
        ,w_p(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,       &
             w_p_levs)                                                  &
                                           ! w at selected pressures
        ,t_p(udims%i_start:udims%i_end,                                 &
             vdims%j_start:vdims%j_end,t_p_levs)                        &
                                           ! t at selected pressures
        ,q_p(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,       &
             q_p_levs)                                                  &
                                           ! q at selected pressures
        ,rh_p(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,      &
              rh_p_levs)                                                &
                                             ! rh at selected pressures
        ,z_p(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,       &
             z_p_levs)                                                  &
                                           ! z at selected pressures
        ,om_p(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,      &
              om_p_levs)
                                             ! om at selected pressures

      REAL, INTENT(OUT) ::                                              &
        u_m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            u_m_levs)                                                   &
                                        ! u at selected model lev
        ,v_m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             v_m_levs)                                                  &
                                         ! v at selected model lev
        ,w_m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             w_m_levs)                                                  &
                                         ! w at selected model lev
        ,t_m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             t_m_levs)                                                  &
                                         ! t at selected model lev
        ,q_m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             q_m_levs)                                                  &
                                         ! q at selected model lev
        ,z_m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             z_m_levs)                                                  &
                                         ! z at selected model lev
        ,ke_m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              ke_m_levs)                                                &
                                          ! ke at selected model lev
        ,om_m(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              om_m_levs)                                                &
                                          ! omega at selected model lev
        ,wbig_m(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,    &
                wbig_m_levs)                                            &
                                             !wbig at selected mod lev
        ,wbig2_m(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,   &
                 wbig2_m_levs)                                          &
                                               !wbig at selected mod lev
        ,RH_m(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,      &
              Rh_m_levs)                                                &
                                           ! RH at selected model lev
        ,dry_mass_m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,&
                    dry_mass_m_levs)                                    &
                                                     ! dry mass mod lev
        ,u_mm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              u_mm_levs)                                                &
                                           ! u at selected model lev
        ,v_mm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              v_mm_levs)                                                &
                                           ! v at selected model lev
        ,w_mm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              w_mm_levs)                                                &
                                           ! w at selected model lev
        ,t_mm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              t_mm_levs)                                                &
                                           ! t at selected model lev
        ,q_mm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              q_mm_levs)                                                &
                                           ! q at selected model lev
        ,z_mm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              z_mm_levs)                                                &
                                           ! z at selected model lev
        ,ke_mm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               ke_mm_levs)                                              &
                                            ! ke at selected model lev
        ,u_inc(udims%i_start:udims%i_end,udims%j_start:udims%j_end,     &
               u_inc_levs)                                              &
                                             ! u at selected model lev
        ,v_inc(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,     &
               v_inc_levs)                                              &
                                               ! v at selected model lev
        ,w_inc(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,     &
               w_inc_levs)                                              &
                                             ! w at selected model lev
        ,t_inc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               t_inc_levs)                                              &
                                             ! t at selected model lev
        ,q_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,     &
               q_inc_levs)                                              &
                                             ! q at selected model lev
        ,qcl_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,   &
                 qcl_inc_levs)                                          &
                                                 ! qcl at  model lev
        ,qcf_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,   &
                 qcf_inc_levs)                                          &
                                                 ! qcf at  model lev
        ,rho_inc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,   &
                 rho_inc_levs)                                          &
                                                 ! rho at  model lev
        ,qrain_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end, &
                   qrain_inc_levs)                                      &
                                             ! qrain at selected model lev
        ,qgraup_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,&
                    qgraup_inc_levs)                                    &
                                             ! qgraup at  model lev
        ,qcft2_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end, &
                   qqcf2_inc_levs)                                      &
                                             ! qcf2 at  model lev
        ,u2_inc(udims%i_start:udims%i_end,udims%j_start:udims%j_end,    &
                u2_inc_levs)                                            &
                                              ! u2 at selected mod lev
        ,v2_inc(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,    &
                v2_inc_levs)                                            &
                                                ! v2 at selected mod lev
        ,w2_inc(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,    &
                w2_inc_levs)                                            &
                                              ! w2 at selected mod lev
        ,t2_inc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                t2_inc_levs)                                            &
                                              ! t2 at selected mod lev
        ,q2_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,    &
                q2_inc_levs)                                            &
                                              ! q2 at selected mod lev
        ,qcl2_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,  &
                  qcl2_inc_levs)                                        &
                                                  ! qcl2 at selected lev
        ,qcf2_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,  &
                  qcf2_inc_levs)                                        &
                                                  ! qcf2 at selected lev
        ,rho2_inc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  rho2_inc_levs)                                        &
                                                  ! rho2 at  model lev
        ,teot_m(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                teot_m_levs)                                            &
                                               ! T at eot/selected
        ,heavy_p(udims%i_start:udims%i_end,                             &
                 vdims%j_start:vdims%j_end,  heavy_p_levs)              &
                                                   ! heavy at pressures
        ,tv_p(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,      &
              tv_p_levs)                                                &
                                             ! tv at pressures
        ,tvom_p(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,    &
                tvom_p_levs)                                            &
                                                 ! tv at pressures
        ,vstarbar(vdims%j_start:vdims%j_end, u_p_levs)                  &
        ,wstarbar(vdims%j_start:vdims%j_end, u_p_levs)                  &
        ,divF(vdims%j_start:vdims%j_end, u_p_levs)                      &
        ,Fy(vdims%j_start:vdims%j_end, u_p_levs)                        &
        ,Fz(vdims%j_start:vdims%j_end, u_p_levs)                        &
        ,zvptp(vdims%j_start:vdims%j_end, u_p_levs)                     &
        ,zupvp(vdims%j_start:vdims%j_end, u_p_levs)                     &
        ,tot_ke(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,tot_ke_w(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
        ,tot_cvt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
        ,tot_gr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,tcm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)       &
        ,tcmq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
        ,tcqcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,tcqcf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,tmf_u(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,tmf_v(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,tmf_w(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,mtorque(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
        ,m_spd_ew(udims%i_start:udims%i_end,udims%j_start:udims%j_end)  &
        ,m_spd_ns(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)  &
        ,am_m1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,am_m2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,am_m3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,am_w1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,am_w2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,am_w3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,o_pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
        ,o_pstar_uv(udims%i_start:udims%i_end,                          &
                    vdims%j_start:vdims%j_end)                          &
        ,encorr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,trop_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,trop_t(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,trop_z(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,trop_icao(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
        ,col_ugz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
        ,col_vgz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
        ,col_wgz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
        ,col_ut(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,col_vt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,col_wt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,col_uq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,col_vq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,col_wq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,col_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,col_uw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,col_vw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,col_uke(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
        ,col_vke(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
        ,col_wke(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
        ,col_u(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,col_v(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,col_w(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
        ,col_sat1(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)


! SCM Dummy variables to keep call to tropin consistent.
       REAL                                                             &
       scm_dummy_1d(1,1)                                                &
      ,scm_dummy_2d(1,1,0:pdims%k_end)

! Local parameters:
      Logical ::                                                        &
     &   qu_m                                                           &
                ! Flag for U wind component on model levels
     &  ,qprod_p(nitems) ! Flag for prods on pressure level

      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='Eot_diag')

! Local scalars:
!   ErrorStatus
      INTEGER      ErrorStatus                                          &
                                        ! Error flag (0 = OK)
        ,i,j,k                                                          &
                                ! loop counters
        ,x,y,z                                                          &
                                ! loop counters
        ,gi,gj                                                          &
        ,interp_order                                                   &
                                !  order of vertical interpolation
        ,iprod,jprod,code30,                                            &
        counter

      CHARACTER(LEN=256)                                                     &
     & CMessage         ! Error message if return code >0

      LOGICAL                                                           &
     & LAM              ! T: limited area model

      REAL                                                              &
     &  dummy                                                           &
                                ! dummy argument - not referenced
     &  ,pressure_pa                                                    &
                                ! pressure in pascals
     &  ,pressure_ex                                                    &
                                ! exner pressure
     &  ,factor                                                         &
                ! scaling factor
     &  ,grid_factor                                                    &
                     ! grid factor
     &  ,cos_long                                                       &
     &  ,sin_long                                                       &
     &  ,sin_lat                                                        &
     &  ,weight1,weight2,weight3,ww1,ww2,temp

! Local dynamic arrays:
      REAL ::                                                           &
         exner_at_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,&
                    udims%k_start:udims%k_end)                          &
                                                    ! exner at u points
        ,exner_at_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,&
                    vdims%k_start:vdims%k_end)                          &
                                                    ! exner at v points
        ,pstar_uv(udims%i_start:udims%i_end,                            &
                  vdims%j_start:vdims%j_end)                            &
                                     ! pstar at uv points 
        ,work_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
                                 ! workspace
        ,u_on_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                pdims%k_start:pdims%k_end)                              &
                                 ! workspace
        ,v_on_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                pdims%k_start:pdims%k_end)                              &
                                 ! workspace
        ,work_om(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 1:tdims%k_end)                                         &
                                              ! workspace for omega
        ,work_rho(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  pdims%k_start:pdims%k_end)                            &
                                                ! workspace
        ,delr_rho(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  pdims%k_start:pdims%k_end)                            &
                                                ! difference
                                ! between rho levels
        ,pstar_halo(pdims_s%i_start:pdims_s%i_end,                      &
                    pdims_s%j_start:pdims_s%j_end)                      &
        ,field_rho(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, &
                   pdims%k_start:pdims%k_end,nmodel_diags)              &
                                            ! Fields on rho grid points
        ,mass_we(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,   &
                 pdims%k_start:pdims%k_end)  ! mass weighting

      INTEGER :: top

      REAL :: temperature(udims%i_start:udims%i_end,                    &
                          vdims%j_start:vdims%j_end, t_p_levs)
                          
                          
      REAL :: work_u(udims%i_start:udims%i_end,                         &
                       udims%j_start:udims%j_end)  
                            ! workspace, used for u & v
             
      REAL :: work_v(vdims%i_start:vdims%i_end,                         &
                       vdims%j_start:vdims%j_end)  
                            ! workspace, used for u & v
             

! Tropopause variables adapted from subroutine O3_to_3D
      INTEGER :: trindx(row_length, rows)
!         Points to the layer boundary below that containing the
!         tropopause

      REAL ::   &
       qsat_work(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,   &
                 1:qdims%k_end)                & ! qsat on model wet levels 
      ,rho_only(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                pdims%k_start:pdims%k_end)     & ! rho only
      ,rho_theta(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 1:tdims%k_end)   ! rho on theta levels

      REAL ::   &
       temp_mass          ! used in col saturation factor calculation

      REAL                                                              &
        cossq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
        ,slcp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
        ,clcp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
        ,spslcp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,spclcp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
        ,r3dr,rocos

      REAL, POINTER ::                                                  &
! prod1 at selected pressures
        field1_p(:,:,:)                                                 &
! prod2 at selected pressures
        ,field2_p(:,:,:)

      REAL                                                              &
     &  mag_vector_np(model_levels)                                     &
                                     ! magnitude of the vector wind NP
     &, dir_vector_np(model_levels)                                     &
                                     ! direction of the vector wind NP
     &, mag_vector_sp(model_levels)                                     &
                                     ! magnitude of the vector wind SP
     &, dir_vector_sp(model_levels)  ! direction of the vector wind SP


! Local dynamic arrays:
      REAL                                                              &
        p_at_theta(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   1:tdims%k_end)                                       &
                                       ! pressure at theta points
        ,exner_at_theta(tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,1:tdims%k_end)        &
                                          ! exner at theta points
        ,rh(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,        &
            1:tdims%k_end)                                              &
                                          ! workspace for RH
        ,T(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           1:tdims%k_end)                                               &
                                         ! T at theta levels
        , true_rho_halo(pdims_s%i_start:pdims_s%i_end,                  &
                        pdims_s%j_start:pdims_s%j_end,                  &
                        pdims_s%k_start:pdims_s%k_end)                  &
        , tmp_data_halo(tdims_s%i_start:tdims_s%i_end,                  &
                        tdims_s%j_start:tdims_s%j_end,                  &
                        1:tdims_s%k_end)


      logical :: qneed_u,qneed_v,qneed_w,qneed_t,qneed_q,qneed_ke       &
     &  , qneed_tropopause,qneed_z
!
! Tropopause variable
      logical QMODEL_HALF_HEIGHT ! Interim local definition
                                 ! Later pass as argument?

      INTEGER N_SUMS
      PARAMETER (N_SUMS=10)   ! number of global summations

      REAL                                                              &
       vert_int_array(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end,n_sums)                 &
                                                 ! array to be summed
      ,flux_array(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)
                                                 ! flux array

      REAL                                                              &
        rho_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                pdims%k_start:pdims%k_end)                              &
                                                    ! rho dry
      , dry_to_wet(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, &
                   pdims%k_start:pdims%k_end)                           &
      , wbig,wbig2                    ! restriction on w

      INTEGER                                                           &
     & n_proc                                                           &
                          ! Total number of processors       (dummy)
     &, n_procx                                                         &
                           ! Number of processors in longitude (dummy)
     &, n_procy            ! Number of processors in latitude  (dummy)

! pointers for vert_int_array
      INTEGER, PARAMETER ::                                             &
     & ip_dry_mass = 1                                                  &
                            ! dry mass
     &,ip_wet_mass = 2                                                  &
                            ! wet mass
     &,ip_cvT      = 3                                                  &
                            ! cvT
     &,ip_gr       = 4                                                  &
                            ! gr
     &,ip_keu      = 5                                                  &
                            ! keu
     &,ip_kev      = 6                                                  &
                            ! kev
     &,ip_kew      = 7                                                  &
                            ! kew
     &,ip_q        = 8                                                  &
                            ! q
     &,ip_qcl      = 9                                                  &
                            ! qcl
     &,ip_qcf      =10                                                  &
                            ! qcf
     &,ip_qu   = 1                                                      &
                        ! qu
     &,ip_qv   = 2                                                      &
                        ! qv
     &,ip_qw   = 3      ! qw

! pointers for rho model level fields and products
      INTEGER, PARAMETER ::                                             &
     &  irho_u = 1                                                      &
     &, irho_v = 2                                                      &
     &, irho_w = 3                                                      &
     &, irho_t = 4                                                      &
     &, irho_q = 5                                                      &
     &, irho_z = 6                                                      &
     &, irho_ke= 7

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header

      IF (lhook) CALL dr_hook('EOT_DIAG',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

! Set Error code to zero
      ErrorStatus = 0

!Set up variable obtained from sf(1,30)
      qprod_p=qsf
      qu_m=qsf(001)

! Factor to scale mass weighted fields by
! this is done because otherwise fields like u*z can take the same value
! as the missing data indicator which prevents meaning.

      factor=1.e-5

! Set order of vertical interpolation
      interp_order = interp_order_linear

! Determine whether limited area model
      IF (model_domain  ==  mt_lam .OR.                                 &
          model_domain  ==  mt_cyclic_lam .OR.                          &
          model_domain  ==  mt_bi_cyclic_lam) THEN
        lam = .TRUE.
      ELSE
        lam = .FALSE.
      END IF

! Set SCM dummy values to zero
       scm_dummy_1d(:,:)   = 0.0
       scm_dummy_2d(:,:,:) = 0.0

!set up logicals for model fields
      qneed_ke=qke_m.OR.qke_mm.OR.qprod_p(17).OR.qprod_p(27).OR.        &
               qprod_p(37).OR.qprod_p(47).OR.qprod_p(57).OR.            &
               qprod_p(67).OR.qprod_p(77).OR.                           &
               qcol_uke.OR.qcol_vke.OR.qcol_wke

      IF (qneed_ke) THEN
        qneed_u=.TRUE.
        qneed_v=.TRUE.
        qneed_w=.TRUE.
      ELSE
        qneed_u=qprod_p(11).OR.qprod_p(12).OR.qprod_p(13).OR.           &
                qprod_p(14).OR.qprod_p(15).OR.qprod_p(16).OR.           &
                qu_m.OR.qu_mm.OR.qcol_u.OR.qcol_ugz.OR.                 &
                qcol_ut.OR.qcol_uq.OR.qcol_uv.OR.qcol_uw

        qneed_v=qprod_p(12).OR.qprod_p(22).OR.qprod_p(23).OR.           &
                qprod_p(24).OR.qprod_p(25).OR.qprod_p(26).OR.           &
                qv_m.OR.qv_mm.OR.qcol_v.OR.qcol_vgz.OR.                 &
                qcol_vt.OR.qcol_vq.OR.qcol_vw.OR.qcol_uv

        qneed_w=qprod_p(13).OR.qprod_p(23).OR.qprod_p(33).OR.           &
                qprod_p(34).OR.qprod_p(35).OR.qprod_p(36).OR.           &
                qw_m.OR.qw_mm .OR.qcol_w.OR.qcol_wgz.OR.                &
                qcol_wt.OR.qcol_wq.OR.qcol_vw.OR.qcol_uw

      END IF  ! test on qneed_ke

      qneed_tropopause = qtrop_p .OR. qtrop_t                           &
        .OR. qtrop_z  .OR. qtrop_icao

      qneed_t=qprod_p(14).OR.qprod_p(24).OR.qprod_p(34).OR.             &
              qprod_p(44).OR.qprod_p(45).OR.qprod_p(46).OR.             &
              qt_m.OR.qt_mm.OR.qcol_wt.OR.                              &
              qcol_ut.OR.qcol_vt.OR. qneed_tropopause .OR.              &
              qcol_sat1 

      qneed_q=qprod_p(15).OR.qprod_p(25).OR.qprod_p(35).OR.             &
              qprod_p(45).OR.qprod_p(55).OR.qprod_p(56).OR.             &
              qq_m.OR.qq_mm .OR.qcol_uq.OR.qcol_vq.OR.qcol_wq

      qneed_z=qprod_p(16).OR.qprod_p(26).OR.qprod_p(36).OR.             &
              qprod_p(46).OR.qprod_p(56).OR.qprod_p(66).OR.             &
              qz_m.OR.qz_mm .OR.qcol_ugz.OR.qcol_vgz.OR.qcol_wgz



!  Calculate exner at theta points
      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            exner_at_theta(i,j,k) = exner_theta_levels(i,j,k)
          END DO                 ! i
        END DO                   ! j
      END DO                     ! k


      IF (qneed_u.OR.qneed_v.OR.qneed_w.OR.                             &
        qneed_t.OR.qneed_q.OR.qdry_mass_m.OR.                           &
        qam_m1.OR.qam_m2.OR.qam_m3.OR.                                  &
        qam_w1.OR.qam_w2.OR.qam_w3) THEN

! DEPENDS ON: swap_bounds
        CALL swap_bounds(u,row_length,rows,model_levels,                &
                         offx,offy,fld_type_u,.TRUE.)

! DEPENDS ON: swap_bounds
        CALL swap_bounds(v,row_length,n_rows,model_levels,              &
                         offx,offy,fld_type_v,.TRUE.)

      CALL u_to_p (u,                                                   &
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity,                      &
                        field_rho(1,1,1,irho_u))

      CALL v_to_p (v,                                                   &
                        vdims_s%i_start,vdims_s%i_end,                  & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity,                      &
                        field_rho(1,1,1,irho_v))

IF (.NOT. l_vatpoles) THEN 
! problem of u & v at poles
        IF (model_domain  ==  1 ) THEN

! DEPENDS ON: polar_vector_wind_n
          CALL Polar_vector_wind_n                                      &
     &      (v,                                                         &
     &      sin_theta_longitude,cos_theta_longitude,                    &
     &      row_length,n_rows, model_levels,                            &
     &      mag_vector_np,dir_vector_np,                                &
     &      mag_vector_sp,dir_vector_sp,                                &
     &      offx, offy, global_row_length,                              &
     &      gc_proc_row_group, at_extremity)

          IF (at_extremity(PSouth) ) THEN
            DO k = pdims%k_start, pdims%k_end
              DO i = pdims%i_start, pdims%i_end
                field_rho(i,1,k,irho_u) = 0.0
              END DO
            END DO
            DO k = pdims%k_start, pdims%k_end
              DO i = pdims%i_start, pdims%i_end
                field_rho(i,1,k,irho_v) = mag_vector_sp(k)
              END DO
            END DO
          END IF

          IF (at_extremity(PNorth) ) THEN
            DO k = pdims%k_start, pdims%k_end
              DO i = pdims%i_start, pdims%i_end
                field_rho(i,rows,k,irho_u) = 0.0
              END DO
            END DO
            DO k = pdims%k_start, pdims%k_end
              DO i = pdims%i_start, pdims%i_end
                field_rho(i,rows,k,irho_v) = mag_vector_np(k)
              END DO
            END DO
          END IF
        END IF
END IF ! vatpoles
!  Calculate delr_rho
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              delr_rho(i,j,k)=r_theta_levels(i,j,k)-                    &
                r_theta_levels(i,j,k-1)
            END DO               ! i
          END DO                 ! j
        END DO                   ! k
!----------------------------------------------------------------------
! convert rho to rho dry in the same way as done in dynamics - Flux_rho
!----------------------------------------------------------------------
        k = 1
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
            temp = 1. - q(i,j,k)- qcl(i,j,k) - qcf(i,j,k)
            rho_dry(i,j,k) = rho(i,j,k) * temp
            dry_to_wet(i,j,k)= 1. / temp
          END DO
        END DO

        DO k = 2, qdims%k_end
          DO j = qdims%j_start, qdims%j_end
            DO i = qdims%i_start, qdims%i_end
              weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
              weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            temp = ( weight2 *                                          &
     &             (1. - q(i,j,k)- qcl(i,j,k) - qcf(i,j,k) ) +          &
     &               weight1 *                                          &
     &             (1. - q(i,j,k-1)- qcl(i,j,k-1) - qcf(i,j,k-1) ) )    &
     &             / weight3
              dry_to_wet(i,j,k)= 1./temp
              rho_dry(i,j,k) = rho(i,j,k) * temp
            END DO
          END DO
        END DO

        DO k = qdims%k_end+1,pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              rho_dry(i,j,k) = rho(i,j,k)
              dry_to_wet(i,j,k)= 1.
            END DO
          END DO
        END DO
!-----------------------------------------------------------------------
! Mass weighting factor for each location ie rho_dry*r*r*dr/(a*a)
! Note currently using rho_dry not rho_wet
! Note rho_dry includes r*r factor
!-----------------------------------------------------------------------

       DO  k=pdims%k_start,pdims%k_end
         DO j=pdims%j_start,pdims%j_end
           DO i=pdims%i_start,pdims%i_end
              mass_we(i,j,k)=rho_dry(i,j,k)*delr_rho(i,j,k)             &
                                      /(earth_radius*earth_radius)
           END DO
         END DO
       END DO
      END IF     ! test on logicals


! ----------------------------------------------------------------------
! 115 dry mass weighting on model rho levels
! ----------------------------------------------------------------------
      IF(qdry_mass_m) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (dry_mass_m_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                dry_mass_m(i,j,counter)=mass_we(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF
!-----------------------------------------------------------------------
!  U and V on pressure levels
!-----------------------------------------------------------------------

      IF (qu_p.OR.qv_p) THEN
! Need call to swapbounds for exner_rho before interpolation
! DEPENDS ON: swap_bounds
        CALL swap_bounds(exner_rho_levels,pdims%i_end,pdims%j_end,      &
                         pdims%k_end,offx,offy,fld_type_p,.FALSE.)
      END IF

      IF(qu_p) THEN

!  Calculate exner at u points. Store in exner_at_u

      CALL p_to_u(exner_rho_levels,                               &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   udims%k_start, udims%k_end,                    &
                   exner_at_u)

      END IF ! on STASHflag

      IF(qv_p) THEN
!  Calculate exner at v points. Store in exner_at_v


      CALL p_to_v(exner_rho_levels,                               &   
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   vdims%k_start, vdims%k_end,                    &
                   exner_at_v) 

      END IF ! on STASHflag

!   Calculate temperature at theta points
      IF(qt_p .OR. qRH_p.OR.qneed_t.OR.qteot_m .OR. qRH_m) THEN

        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
              p_at_theta(i,j,k) = p_theta_levels(i,j,k)
              rh(i,j,k) = 0.
            END DO               ! i
          END DO                 ! j
        END DO                   ! k

      END IF ! on STASHflag

!  height of model rho levels
!   Remove radius of earth from height field
      IF(qz_p.OR.qneed_z) THEN

        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              field_rho(i,j,k,irho_z) = r_rho_levels(i,j,k)             &
     &                                         - earth_radius
            END DO               ! i
          END DO                 ! j
        END DO                   ! k

      END IF ! on relevant STASHflags

! extra points DOne because of LAM DOmains
! Need to come back and check this for both EG and ND
      DO j = 1, rows+1
        DO i = 1,row_length+1
          pstar_halo(i,j) = p(i,j,1) + g * rho(i,j,1)                   &
     &      * (r_rho_levels(i,j,1) -                                    &
     &      r_theta_levels(i,j,0) ) /                                   &
     &      (r_rho_levels(i,j,1) *                                      &
     &      r_rho_levels(i,j,1) )
        END DO
      END DO

! DEPENDS ON: swap_bounds
      CALL Swap_Bounds(pstar_halo,                                      &
        row_length, rows, 1,                                            &
        offx, offy, fld_type_p, .FALSE.)

!CDIR NOUNROLL
      DO j = vdims%j_start, vdims%j_end
        DO i = udims%i_start, udims%i_end
          pstar_uv(i,j)=(pstar_halo(i,j)+                               &
     &      pstar_halo(i+1,j)+                                          &
     &      pstar_halo(i,j+1)+                                          &
     &      pstar_halo(i+1,j+1)) * 0.25
        END DO                   ! i
      END DO                     ! j

! ----------------------------------------------------------------------
! STASH items 001,002 : u,v winds on model levels
! ----------------------------------------------------------------------


      IF(qu_m) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (u_m_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                u_m(i,j,counter)=field_rho(i,j,k,irho_u)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF
      IF(qu_mm) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (u_mm_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                u_mm(i,j,counter)=field_rho(i,j,k,irho_u)*mass_we(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO

      END IF ! on STASHflag

      IF(qv_m) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (v_m_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                v_m(i,j,counter)=field_rho(i,j,k,irho_v)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF
      IF(qv_mm) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (v_mm_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                v_mm(i,j,counter)=field_rho(i,j,k,irho_v)*mass_we(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 003,004,005 : w wind, t and q on model levels
! ----------------------------------------------------------------------
      k=1       !  level 1
      IF (qneed_t) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i=pdims%i_start,pdims%i_end
              field_rho(i,j,k,irho_t)=t(i,j,k)
          END DO
        END DO
      END IF
      IF (qneed_q) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i=pdims%i_start,pdims%i_end
              field_rho(i,j,k,irho_q)=q(i,j,k)
          END DO
        END DO
      END IF
      IF (qneed_w) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i=pdims%i_start,pdims%i_end
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            ww2 = weight2/weight3
            field_rho(i,j,k,irho_w)=ww2 * w(i,j,k)  + ww1 * w(i,j,k-1)
          END DO
        END DO
      END IF

      IF (qneed_w.OR.qneed_t) THEN
        DO k=2,pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i=pdims%i_start,pdims%i_end
              weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
              weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
              ww1 = weight1/weight3
              ww2 = weight2/weight3
              IF (qneed_t) THEN
                field_rho(i,j,k,irho_t)=ww2 * t(i,j,k) + ww1*t(i,j,k-1)
              END IF
              IF (qneed_w) THEN
                field_rho(i,j,k,irho_w)=ww2 * w(i,j,k) + ww1*w(i,j,k-1)
              END IF
            END DO
          END DO
        END DO
      END IF

      IF (qneed_q) THEN
        DO k=2,qdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i=pdims%i_start,pdims%i_end
              weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
              weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
              ww1 = weight1/weight3
              ww2 = weight2/weight3
              field_rho(i,j,k,irho_q)=ww2 * q(i,j,k)+ ww1 * q(i,j,k-1)
            END DO
          END DO
        END DO
      END IF

! ----------------------------------------------------------------------
! 003 w on model rho levels
! ----------------------------------------------------------------------
      IF(qw_m) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (w_m_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                w_m(i,j,counter)=field_rho(i,j,k,irho_w)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF
      IF(qw_mm) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (w_mm_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                w_mm(i,j,counter)=field_rho(i,j,k,irho_w)*mass_we(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF ! on STASHflag
! ----------------------------------------------------------------------
! 004 T on model rho levels
! ----------------------------------------------------------------------
      IF(qt_m) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (t_m_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                t_m(i,j,counter)=field_rho(i,j,k,irho_t)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF
      IF(qt_mm) THEN
            counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (t_mm_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                t_mm(i,j,counter)=field_rho(i,j,k,irho_t)*mass_we(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF ! on STASHflag
! ----------------------------------------------------------------------
! 005 q on model rho levels
! ----------------------------------------------------------------------
      IF(qq_m) THEN
        counter=1
        DO  k=1,qdims%k_end
          IF (q_m_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                q_m(i,j,counter)=field_rho(i,j,k,irho_q)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF
      IF(qq_mm) THEN
        counter=1
        DO  k=1,qdims%k_end
          IF (q_mm_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                q_mm(i,j,counter)=field_rho(i,j,k,irho_q)*mass_we(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF ! on STASHflag
! ----------------------------------------------------------------------
! 006 z on model rho levels
! ----------------------------------------------------------------------
      IF(qz_m) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (z_m_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                z_m(i,j,counter)=field_rho(i,j,k,irho_z)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF
      IF(qz_mm) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (z_mm_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                z_mm(i,j,counter)=field_rho(i,j,k,irho_z)*mass_we(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF ! on STASHflag
! ----------------------------------------------------------------------
! 007 KE on model rho levels
! ----------------------------------------------------------------------
      IF (qneed_ke) THEN
        DO K=pdims%k_start,pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO I=pdims%i_start,pdims%i_end
              field_rho(i,j,k,irho_ke) = 0.5*(                          &
     &        field_rho(i,j,k,irho_u)*field_rho(i,j,k,irho_u)+          &
     &        field_rho(i,j,k,irho_v)*field_rho(i,j,k,irho_v)+          &
     &        field_rho(i,j,k,irho_w)*field_rho(i,j,k,irho_w))
            END DO
          END DO
        END DO
      END IF

      IF(qke_m) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (ke_m_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                ke_m(i,j,counter)=field_rho(i,j,k,irho_ke)
                    END DO
                  END DO
                  counter=counter+1
          END IF
              END DO
      END IF
      IF(qke_mm) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (ke_mm_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
               ke_mm(i,j,counter)=field_rho(i,j,k,irho_ke)*mass_we(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items :products on model levels
! At present all products mass weighted using dry density
! ----------------------------------------------------------------------

      DO iprod=1,nmodel_diags
        DO jprod=iprod,nmodel_diags
          code30=iprod*10+jprod
          IF (qprod_p(code30)) THEN
              counter=1
              DO  k=pdims%k_start,pdims%k_end !loop over levels
                IF (prod_m_list(k,iprod,jprod)) THEN
                  DO j = pdims%j_start, pdims%j_end
                    DO i = pdims%i_start, pdims%i_end
                    stashwork(si(code30)+ (i-pdims%i_start)+            &
                                    (j-pdims%j_start)*row_length+       &
                        (counter-1)*rows*row_length)=                   &
                        field_rho(i,j,k,iprod)*field_rho(i,j,k,jprod)   &
                        *factor*mass_we(i,j,k)
                    END DO
                  END DO
                  counter=counter+1
                END IF
              END DO
            END IF
        END DO
      END DO


! ----------------------------------------------------------------------
! STASH item 111 : T at eot
! ----------------------------------------------------------------------
      IF (qteot_m) THEN
        counter=1
        DO  k=1,tdims%k_end
          IF (teot_m_list(k)) THEN
            DO j=tdims%j_start,tdims%j_end
              DO i=tdims%i_start,tdims%i_end
                teot_m(i,j,counter)=t(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 181 : T increment
! ----------------------------------------------------------------------
      IF (qt_inc) THEN
        counter=1
        DO  k=1,tdims%k_end
          IF (t_inc_list(k)) THEN
            DO j=tdims%j_start,tdims%j_end
              DO i=tdims%i_start,tdims%i_end
                t_inc(i,j,counter)=inc_t(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 182 : Q increment
! ----------------------------------------------------------------------
      IF (qq_inc) THEN
        counter=1
        DO  k=1,qdims%k_end
          IF (q_inc_list(k)) THEN
            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                q_inc(i,j,counter)=inc_q(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 183 : QCL increment
! ----------------------------------------------------------------------
      IF (qqcl_inc) THEN
        counter=1
        DO  k=1,qdims%k_end
          IF (qcl_inc_list(k)) THEN
            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                qcl_inc(i,j,counter)=inc_qcl(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 184 : QCF increment
! ----------------------------------------------------------------------
      IF (qqcf_inc) THEN
        counter=1
        DO  k=1,qdims%k_end
          IF (qcf_inc_list(k)) THEN
            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                qcf_inc(i,j,counter)=inc_qcf(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 185 : U increment
! ----------------------------------------------------------------------
      IF (qu_inc) THEN
        counter=1
        DO  k=udims%k_start,udims%k_end
          IF (u_inc_list(k)) THEN
            DO j=udims%j_start,udims%j_end
              DO i=udims%i_start,udims%i_end
                u_inc(i,j,counter)=inc_u(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 186 : V increment
! ----------------------------------------------------------------------
      IF (qv_inc) THEN
        counter=1
        DO  k=vdims%k_start,vdims%k_end
          IF (v_inc_list(k)) THEN
            DO j=vdims%j_start,vdims%j_end
              DO i=vdims%i_start,vdims%i_end
                v_inc(i,j,counter)=inc_v(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 187 : W increment
! ----------------------------------------------------------------------
      IF (qw_inc) THEN
        counter=1
! Zeroth level not output as w is zero as per boundary conditions
        DO  k=1,wdims%k_end
          IF (w_inc_list(k)) THEN
            DO j=wdims%j_start,wdims%j_end
              DO i=wdims%i_start,wdims%i_end
                w_inc(i,j,counter)=inc_w(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 188 : RHO increment
! ----------------------------------------------------------------------
      IF (qrho_inc) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (rho_inc_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
                rho_inc(i,j,counter)=inc_rho(i,j,k)/                      &
     &            (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 189 : qrain increment
! ----------------------------------------------------------------------
      IF (qqrain_inc) THEN
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(inc_qrain,                                     &
                         qdims%i_end, qdims%j_end, qdims%k_end,         &
                         halo_i, halo_j, fld_type_p, .FALSE.)
        counter=1
        DO  k=1,qdims%k_end
          IF (qrain_inc_list(k)) THEN
            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                qrain_inc(i,j,counter)=inc_qrain(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 190 : qgraup increment
! ----------------------------------------------------------------------
      IF (qqgraup_inc) THEN
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(inc_qgraup,                                    &
                         qdims%i_end, qdims%j_end, qdims%k_end,         &
                         halo_i, halo_j, fld_type_p, .FALSE.)
        counter=1
        DO  k=1,qdims%k_end
          IF (qgraup_inc_list(k)) THEN
            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                qgraup_inc(i,j,counter)=inc_qgraup(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 191 : qcf2 increment
! ----------------------------------------------------------------------
      IF (qqqcf2_inc) THEN
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(inc_qcf2,                                      &
                         qdims%i_end, qdims%j_end, qdims%k_end,         &
                         halo_i, halo_j, fld_type_p, .FALSE.)
        counter=1
        DO  k=1,qdims%k_end
          IF (qqcf2_inc_list(k)) THEN
            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                qcft2_inc(i,j,counter)=inc_qcf2(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 171 : T increment**2
! ----------------------------------------------------------------------
      IF (qt2_inc) THEN
        counter=1
        DO  k=1,tdims%k_end
          IF (t2_inc_list(k)) THEN
            DO j=tdims%j_start,tdims%j_end
              DO i=tdims%i_start,tdims%i_end
                t2_inc(i,j,counter)=inc_t(i,j,k)*inc_t(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 172 : Q increment **2
! ----------------------------------------------------------------------
      IF (qq2_inc) THEN
        counter=1
        DO  k=1,qdims%k_end
          IF (q2_inc_list(k)) THEN
            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                q2_inc(i,j,counter)=inc_q(i,j,k)*inc_q(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 173 : Qcl increment **2
! ----------------------------------------------------------------------
      IF (qqcl2_inc) THEN
        counter=1
        DO  k=1,qdims%k_end
          IF (qcl2_inc_list(k)) THEN
            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                qcl2_inc(i,j,counter)=inc_qcl(i,j,k)*inc_qcl(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 174 : Qcf increment **2
! ----------------------------------------------------------------------
      IF (qqcf2_inc) THEN
        counter=1
        DO  k=1,qdims%k_end
          IF (qcf2_inc_list(k)) THEN
            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                qcf2_inc(i,j,counter)=inc_qcf(i,j,k)*inc_qcf(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 175 : U increment **2
! ----------------------------------------------------------------------
      IF (qu2_inc) THEN
        counter=1
        DO  k=udims%k_start,udims%k_end
          IF (u2_inc_list(k)) THEN
            DO j=udims%j_start,udims%j_end
              DO i=udims%i_start,udims%i_end
                u2_inc(i,j,counter)=inc_u(i,j,k)*inc_u(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 176 : V increment **2
! ----------------------------------------------------------------------
      IF (qv2_inc) THEN
        counter=1
        DO  k=vdims%k_start,vdims%k_end
          IF (v2_inc_list(k)) THEN
            DO j=vdims%j_start,vdims%j_end
              DO i=vdims%i_start,vdims%i_end
                v2_inc(i,j,counter)=inc_v(i,j,k)*inc_v(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 177 : W increment **2
! ----------------------------------------------------------------------
      IF (qw2_inc) THEN
        counter=1
        DO  k=1,wdims%k_end
          IF (w2_inc_list(k)) THEN
            DO j=wdims%j_start,wdims%j_end
              DO i=wdims%i_start,wdims%i_end
                w2_inc(i,j,counter)=inc_w(i,j,k)*inc_w(i,j,k)
              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 178 : RHO increment **2
! ----------------------------------------------------------------------
      IF (qrho2_inc) THEN
        counter=1
        DO  k=pdims%k_start,pdims%k_end
          IF (rho2_inc_list(k)) THEN
            DO j=pdims%j_start,pdims%j_end
              DO i=pdims%i_start,pdims%i_end
               rho2_inc(i,j,counter)=(inc_rho(i,j,k)/                     &
     &            (r_rho_levels(i,j,k)*r_rho_levels(i,j,k)))**2

              END DO
            END DO
            counter=counter+1
          END IF
        END DO
      END IF

! ----------------------------------------------------------------------
! Calculation of RH used by item 206 RH on pressure levels
!                                113 RH on model levels
!                                442 column saturation fraction
! ----------------------------------------------------------------------
      IF(qRH_p .OR. qRH_m  .OR. qcol_sat1) THEN
        DO k = 1, qdims%k_end
!  Find humidity saturation at theta points - store in rh
! DEPENDS ON: qsat
          CALL QSAT(qsat_work(1,1,k),T(1,1,k),p_at_theta(1,1,k),          &
            theta_field_size)
!  And convert to relative humidity
          DO j = qdims%j_start, qdims%j_end
            DO i = qdims%i_start, qdims%i_end

              rh(i,j,k) = q(i,j,k)/qsat_work(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
              IF(rh(i,j,k) <  0.0) THEN
                rh(i,j,k) = 0.
              END IF
            END DO               ! i
          END DO                 ! j
        END DO                   ! k wet_levels
      END IF                     ! on stashflags

! ----------------------------------------------------------------------
! STASH item 008: Omega on model theta levels
!                 Also used by 208 omega on pressure levels
! ----------------------------------------------------------------------
      IF ((qom_p) .OR. (qom_m)) THEN
        DO k = pdims_s%k_start, pdims_s%k_end
          DO j = pdims_s%j_start, pdims_s%j_end
            DO i = pdims_s%i_start, pdims_s%i_end
              true_rho_halo(i,j,k) = rho(i,j,k) /                       &
                (r_rho_levels(i,j,k) *                                  &
                r_rho_levels(i,j,k))
            END DO
          END DO
        END DO

        CALL P_TO_T (                                                   &
          pdims%i_end,pdims%j_end, halo_i, halo_j,                      &
          offx,offy, tdims%k_end-1                                      &
          , r_theta_levels, r_rho_levels                                &
          , true_rho_halo                                               &
          , tmp_data_halo                                               &
          )

! Top theta level set equal to top rho level

        DO j = tdims_s%j_start, tdims_s%j_end
          DO i = tdims_s%i_start, tdims_s%i_end
            tmp_data_halo(i,j,tdims_s%k_end) =                          &
              true_rho_halo(i,j,pdims_s%k_end)
          END DO
        END DO


        DO k = 1, tdims%k_end
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              work_om(i,j,k)= -1*g*tmp_data_halo(i,j,k)*w(i,j,k)
            END DO
          END DO
        END DO

      END IF ! on STASHflag

      IF (qom_m) THEN
! Omega on model theta levels (stash item 008)
        counter = 1
        DO k = 1, tdims%k_end
          IF (om_m_list(k)) THEN
            DO j = tdims%j_start, tdims%j_end
              DO i = tdims%i_start, tdims%i_end
                om_m(i,j,counter) = work_om(i,j,k)
              END DO
            END DO
            counter = counter + 1
          END IF
        END DO
      END IF ! on STASHflag
! ----------------------------------------------------------------------
! Pressure level diagnostics
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! New Dynamics items 201-208 and 301
! ----------------------------------------------------------------------
! STASH items 201 : u wind components on pressure surfaces
! ----------------------------------------------------------------------
      IF(qu_p) THEN

        DO  k=1,u_p_levs

          pressure_pa = u_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          
          CALL vert_interp2 (u, row_length, rows, model_levels          &
                                ,pressure_ex                            &
                                ,offx, offy, 0, 0                       &
                                ,exner_at_u, interp_order               &
                                ,work_u )

                   
! Perform simple horizontal interpolation from 'C' to 'B' grid

          CALL  uC_to_uB(work_u,                                        &
            row_length,rows,n_rows,1,offx,offy,                         &
            u_p(udims%i_start,vdims%j_start,k))

         DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              IF (pstar_uv(i,j) <  u_press(k)*100.) THEN
                u_p(i,j,k)=0.0
              END IF
            END DO
          END DO
        END DO  ! k pressure levels loop
      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 202 : v wind components on pressure surfaces
! ----------------------------------------------------------------------
      IF(qv_p) THEN
      
        DO  k=1,v_p_levs

          pressure_pa = v_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
 
    
          CALL vert_interp2 (v, row_length, n_rows, model_levels        &
                                ,pressure_ex                            &
                                ,offx, offy, 0, 0                       &
                                ,exner_at_v, interp_order               &
                                ,work_v )

! Perform simple horizontal interpolation from 'C' to 'B' grid


          IF ( l_vatpoles ) THEN
            ! In the VATPOLES case, u on the B grid is required by v 
            ! on the B grid to compute polar values, so the pressure 
            ! levels must match.
            
            IF(u_press(k) == v_press(k)) THEN
      
              CALL  vC_to_vB(work_v,                                    &
                rows,row_length,n_rows,1,offx,offy,                     &
                global_row_length,                                      &
                v_p(udims%i_start,vdims%j_start,k)                      &
               ,u_p(udims%i_start,vdims%j_start,k)                      &
              )           

            ELSE
              ErrorStatus = -1        ! Warning
              Cmessage='u and v on B grid pressure levels - '           &
               //'requested pressure levels do not match, diagnostic aborted'

              CALL Ereport(Routinename,ErrorStatus,Cmessage)
              EXIT                    ! jump out of k levels loop            
            END IF

          ELSE

            CALL  vC_to_vB(work_v,                                      &
              rows,row_length,n_rows,1,offx,offy,                       &
              global_row_length,                                        &
              v_p(udims%i_start,vdims%j_start,k)                        &
             )

          END IF  ! vatpoles

          DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              IF (pstar_uv(i,j) <  v_press(k)*100.) THEN
                v_p(i,j,k)=0.0
              END IF
            END DO
          END DO

        END DO  ! k pressure levels loop
      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 203 : w wind components on pressure surfaces
! ----------------------------------------------------------------------
      IF(qw_p) THEN
        DO  k=1,w_p_levs

          pressure_pa = w_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          
          CALL vert_interp2 (w(:,:,1:)                                  &
            ,row_length, rows, model_levels                             &
            ,pressure_ex                                                &
            ,offx, offy, 0,0                                            &
            ,exner_at_theta, interp_order                               &
            ,work_1 )
          
          
! Perform simple horizontal interpolation from 'C' to 'B' grid

          CALL  pC_to_pB(work_1,                                        &
            row_length,rows,n_rows,1,offx,offy,                         &
            w_p(udims%i_start,vdims%j_start,k))
     
          DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              IF (pstar_uv(i,j) <  w_press(k)*100.) THEN
                w_p(i,j,k)=0.0
              END IF
            END DO
          END DO

        END DO  ! k pressure levels loop
      END IF ! on STASHflag


! ----------------------------------------------------------------------
! STASH item 204: temperature              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qt_p) THEN
        DO k = 1, t_p_levs

          pressure_pa = t_press(k)*100. ! convert to Pascals
          
 
          CALL T_vert_interp_to_p(                                      &
            T, theta, row_length, rows                                  &
            ,model_levels, pressure_pa, offx, offy ,halo_i,halo_j       &
            ,p_theta_levels                                             &
            ,bl_levels                                                  &
            ,exner_theta_levels                                         &
            ,work_1 )

          CALL  pC_to_pB(work_1,                                        &
            row_length,rows,n_rows,1,offx,offy,                         &
            t_p(udims%i_start,vdims%j_start,k))


          DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
                temperature(i,j,k) = t_p(i,j,k)
            END DO
          END DO

          DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              IF (pstar_uv(i,j) <  t_press(k)*100.) THEN
                t_p(i,j,k)=0.0
              END IF
            END DO
          END DO

        END DO ! over output STASH pressure levels
      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 205: Q              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qq_p) THEN
        DO k = 1, q_p_levs

          pressure_pa = q_press(k)*100. ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          
          CALL vert_interp2(q(:,:,1:), row_length, rows, model_levels   &
            ,pressure_ex                                                &
            ,halo_i, halo_j, 0,0                                        &
            ,exner_at_theta, interp_order                               &
            ,work_1 )
          

          CALL  pC_to_pB(work_1,                                        &
            row_length,rows,n_rows,1,offx,offy,                         &
            q_p(udims%i_start,vdims%j_start,k))

          DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              IF (pstar_uv(i,j) <  q_press(k)*100.) THEN
                q_p(i,j,k)=0.0
              END IF
            END DO
          END DO

        END DO ! over output STASH pressure levels
      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 206: RH              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qRH_p) THEN

!  Interpolate
        DO k = 1, rh_p_levs

          pressure_pa = RH_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          
          CALL vert_interp2(                                            &
            rh, row_length, rows, model_levels                          &
            ,pressure_ex                                                &
            ,0,0, 0,0                                                   &
            ,exner_at_theta, interp_order                               &
            ,work_1 )

          CALL  pC_to_pB(work_1,                                        &
            row_length,rows,n_rows,1,offx,offy,                         &
            rh_p(udims%i_start,vdims%j_start,k))
     
          DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              IF (pstar_uv(i,j) <  rh_press(k)*100.) THEN
                rh_p(i,j,k)=0.0
              END IF
            END DO
          END DO

        END DO ! k over output STASH pressure levels

      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 207: Z              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qz_p) THEN
        DO k = 1, z_p_levs

         pressure_pa = z_press(k)*100. ! convert to Pascals
         CALL vert_h_onto_p(                                            &
           field_rho(1,1,1,irho_z),row_length, rows, model_levels       &
           ,pressure_pa                                                 &
           ,p_theta_levels                                              &
           ,theta, exner_theta_levels                                   &
           , exner_rho_levels                                           &
           ,bl_levels                                                   &
           ,offx, offy, halo_i, halo_j                                  &
           ,p, interp_order, work_1 )

         CALL  pC_to_pB(work_1,                                         &
           row_length,rows,n_rows,1,offx,offy,                          &
           z_p(udims%i_start,vdims%j_start,k))
     
         DO j = vdims%j_start, vdims%j_end
           DO i = udims%i_start, udims%i_end
             IF (pstar_uv(i,j) <  z_press(k)*100.) THEN
               z_p(i,j,k)=0.0
             END IF
           END DO
         END DO

       END DO                    ! over output STASH pressure levels
      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 208: Omega on pressure surfaces
! ----------------------------------------------------------------------

      IF (qom_p) THEN
! Omega on pressure levels (stash item 208)
      
        DO  k=1,om_p_levs

          pressure_pa = om_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          
          
          CALL vert_interp2 (work_om,                                   &
            row_length, rows, model_levels                              &
            ,pressure_ex                                                &
            ,0,0, 0,0                                                   &
            ,exner_at_theta, interp_order                               &
            ,work_1 )          
           
          
! Perform simple horizontal interpolation from 'C' to 'B' grid

          CALL  pC_to_pB(work_1,                                        &
            row_length,rows,n_rows,1,offx,offy,                         &
            om_p(udims%i_start,vdims%j_start,k))

          DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              IF (pstar_uv(i,j) <  om_press(k)*100.) THEN
                om_p(i,j,k)=0.0
              END IF
            END DO
          END DO
        END DO  ! k pressure levels loop
      END IF ! on STASHflag
! ----------------------------------------------------------------------
! STASH item 301: Heavyside on pressure surfaces
! ----------------------------------------------------------------------
!L The Heavyside function is defined as 1.0 if the pressure level
!L  is above the surface (i.e. pstar) and 0.0 if below. A time mean of
!L  this will give information on the fraction of time a pressure
!L  level is above the land or sea surface.
      IF (qheavy_p) THEN
        DO k=1,heavy_p_levs
          DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              IF (pstar_uv(i,j) <  heavy_press(k)*100.) THEN
                heavy_p(i,j,k)=0.0
              ELSE
                heavy_p(i,j,k)=1.0
              END IF
            END DO
          END DO
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH products items : all pressure surface products
! Loop over the pressure level diags and see if product is required
! If so process
! ----------------------------------------------------------------------
      DO iprod=1,npress_diags ! code of 1st field
        DO jprod=1,npress_diags ! code of 2nd field
          code30=200+iprod*10+jprod ! code of product
          IF(qprod_p(code30)) THEN   ! product required
            IF (iprod  ==  1) THEN
              field1_p=>u_p
            ELSE IF (iprod  ==  2) THEN
              field1_p=>v_p
            ELSE IF (iprod  ==  3) THEN
              field1_p=>w_p
            ELSE IF (iprod  ==  4) THEN
              field1_p=>t_p
            ELSE IF (iprod  ==  5) THEN
              field1_p=>q_p
            ELSE IF (iprod  ==  6) THEN
              field1_p=>rh_p
            ELSE IF (iprod  ==  7) THEN
              field1_p=>z_p
            ELSE IF (iprod  ==  8) THEN
              field1_p=>om_p
            END IF

            IF (jprod  ==  1) THEN
              field2_p=>u_p
            ELSE IF (jprod  ==  2) THEN
              field2_p=>v_p
            ELSE IF (jprod  ==  3) THEN
              field2_p=>w_p
            ELSE IF (jprod  ==  4) THEN
              field2_p=>t_p
            ELSE IF (jprod  ==  5) THEN
              field2_p=>q_p
            ELSE IF (jprod  ==  6) THEN
              field2_p=>rh_p
            ELSE IF (jprod  ==  7) THEN
              field2_p=>z_p
            ELSE IF (jprod  ==  8) THEN
              field2_p=>om_p
            END IF
            
           
            DO  k=1,prod_p_levs(iprod,jprod) !loop over levels
              DO j = vdims%j_start, vdims%j_end
                DO i = udims%i_start, udims%i_end
                  stashwork(si(code30)+                                 &
                    (i-udims%i_start)+(j-vdims%j_start)*row_length+     &
                    (k-1)*n_rows*row_length)=                           &
                    field1_p(i,j,prod_ind(k,iprod,jprod))*              &
                    field2_p(i,j,prod_ind(k+prod_p_levs(iprod,jprod),   &
                    iprod,jprod))
                END DO
              END DO
            END DO
          END IF                 ! on STASHflag
        END DO
      END DO

! ----------------------------------------------------------------------
! STASH item 302 or 282: virtual temp on pressure surfaces
! ----------------------------------------------------------------------
      IF (qt_p .AND. qq_p .AND. qtv_p) THEN
        DO  k=1,prod_p_levs(4,5) !loop over levels shared by t and q
          DO j = vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              tv_p(i,j,k)=t_p(i,j,prod_ind(k,4,5))*                     &
     &          (1+c_virtual*                                           &
     &          q_p(i,j,prod_ind(k+prod_p_levs(4,5),4,5)))
            END DO
          END DO
        END DO
! ----------------------------------------------------------------------
! STASH item 303 or 283: virtual temp* omega on pressure surfaces
! ----------------------------------------------------------------------
        IF (qom_p .AND. qtvom_p) THEN
          DO  k=1,tvom_p_levs
            DO j = vdims%j_start, vdims%j_end
              DO i = udims%i_start, udims%i_end
                tvom_p(i,j,k)=tv_p(i,j,k)*om_p(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
! ----------------------------------------------------------------------
! STASH items 310, 311, 312, 313, 314, 315 and 316
!          or 290, 291, 292, 293, 294, 295 and 296
! Calculate the scaled divergence of the Eliassen-Palm flux from
! u, v, w and temperature on pressure levels.
! Meridional heat flux and momentum fluxes are also calculated.
! ----------------------------------------------------------------------
      IF(qvstarbar_p .OR. qwstarbar_p .OR. qepy_p .OR. qepz_p &
         .OR. qdivep_p .OR. qzvptp_p .OR. qzupvp_p) THEN

       CALL Calc_Div_EP_Flux(                                         & 
          rows, n_rows, row_length, global_row_length, sin_v_latitude,& 
          qvstarbar_p, qwstarbar_p, qepy_p, qepz_p, qdivep_p,         & 
          qzvptp_p, qzupvp_p,                                         & 
          u_press, u_p_levs, v_p_levs, w_p_levs, t_p_levs,            & 
          vstarbar_p_levs, wstarbar_p_levs,  Fy_p_levs, Fz_p_levs,    & 
          divF_p_levs, zvptp_p_levs, zupvp_p_levs,                    & 
          u_p, v_p, w_p, temperature,                                 & 
          vstarbar, wstarbar, Fy, Fz, divF, zvptp, zupvp) 
      END IF

! ----------------------------------------------------------------------
! STASH call routine for all column integral fields
! ----------------------------------------------------------------------

      IF (qtot_ke.OR.qtot_ke_w.OR.qtcm.OR.qtot_cvt.OR.qtot_gr           &
        .OR.qtcmq.OR.qtcqcl.OR.qtcqcf                                   &
        .OR.qtmf_u.OR.qtmf_v.OR.qtmf_w                                  &
        .OR.qencorr ) THEN
        grid_factor=1.0/(earth_radius*earth_radius) 

        CALL vert_eng_massq(                                            &
                            wet_levels,                                 &
           theta ,u,v,w, rho, q,qcl,qcf                                 &
          ,wet_to_dry_n                                                 &
          ,exner_theta_levels                                           &
          ,.TRUE.                                                       &
          ,rho_dry,dry_to_wet                                           &
          ,vert_int_array, flux_array )

      END IF
! ----------------------------------------------------------------------
! STASH item 401: Total KE
! ----------------------------------------------------------------------

      IF (qtot_ke) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tot_ke(i,j)=(vert_int_array(i,j,ip_keu)+                    &
     &        vert_int_array(i,j,ip_kev))*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 402: Total KE with W
! ----------------------------------------------------------------------
      IF (qtot_ke_w) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tot_ke_w(i,j)=(vert_int_array(i,j,ip_keu)+                  &
     &        vert_int_array(i,j,ip_kev)+                               &
     &        vert_int_array(i,j,ip_kew))*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 403: Total column dry mass
! ----------------------------------------------------------------------
      IF (qtcm) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tcm(i,j)=vert_int_array(i,j,ip_dry_mass)*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 404: Total column wet mass
! ----------------------------------------------------------------------
      IF (qtcmq) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tcmq(i,j)=vert_int_array(i,j,ip_wet_mass)*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 405:Total column qcl
! ----------------------------------------------------------------------
      IF (qtcqcl) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tcqcl(i,j)=vert_int_array(i,j,ip_qcl)*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 406: Total column qcf
! ----------------------------------------------------------------------
      IF (qtcqcf) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tcqcf(i,j)=vert_int_array(i,j,ip_qcf)*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 407: Total moisture flux U
! ----------------------------------------------------------------------
      IF (qtmf_u) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tmf_u(i,j)=flux_array(i,j,ip_qu)*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 408: Total moisture flux V
! ----------------------------------------------------------------------
      IF (qtmf_v) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tmf_v(i,j)=flux_array(i,j,ip_qv)*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 409: Total moisture flux W
! ----------------------------------------------------------------------
      IF (qtmf_w) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tmf_w(i,j)=flux_array(i,j,ip_qw)*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 410: Mountain torque at u points on the C-grid
! ----------------------------------------------------------------------

      IF (qmtorque) THEN
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
            mtorque(i,j)= 0.5 * cos_theta_latitude(i,j)                 &
     &                        * ( pstar_halo(i,j) + pstar_halo(i+1,j)  )&
     &                         * (  r_theta_levels(i+1,j,0)             &
     &                            - r_theta_levels( i ,j,0)            )&
     &                         / ( delta_lambda                        )
          END DO
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH item 440 and 441: East-West and North-South components
!                         of the surface pressure drag
! ----------------------------------------------------------------------

      IF (qm_spd_ew) THEN
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
            m_spd_ew(i,j)= 0.5 * ( pstar_halo(i,j) + pstar_halo(i+1,j) )&
     &                         * (  r_theta_levels(i+1,j,0)             &
     &                            - r_theta_levels( i ,j,0)            )&
     &                         / ( delta_lambda * earth_radius         )
          END DO
        END DO
      END IF

      IF (qm_spd_ns) THEN
        DO j = vdims%j_start, vdims%j_end
          DO i = vdims%i_start, vdims%i_end
            m_spd_ns(i,j)= 0.5 * ( pstar_halo(i,j) + pstar_halo(i,j+1) )&
     &                         * (  r_theta_levels(i,j+1,0)             &
     &                            - r_theta_levels(i, j ,0)            )&
     &                         / ( delta_phi    * earth_radius         )
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH items 411-416:Angular momentum
! ----------------------------------------------------------------------
      IF (qam_m1.OR.qam_m2.OR.qam_m3.OR.                                &
        qam_w1.OR.qam_w2.OR.qam_w3) THEN
! constants
        factor=1.e-24

! calculate longitude  & latitude
        DO j=pdims%j_start,pdims%j_end
          gj = datastart(2) + j - 1
          DO i=pdims%i_start,pdims%i_end
            gi = datastart(1) + i - 1
            cossq(i,j)=cos_theta_latitude(i,j)*cos_theta_latitude(i,j)
          END DO
        END DO
        IF (qam_m1.OR.qam_m2.OR.qam_w1.OR.qam_w2) THEN
          DO j=pdims%j_start,pdims%j_end
            DO i=pdims%i_start,pdims%i_end
              spclcp(i,j)=sin_theta_latitude(i,j)*                      &
     &          cos_theta_longitude(i,j)*cos_theta_latitude(i,j)
              spslcp(i,j)=sin_theta_latitude(i,j)*                      &
     &          sin_theta_longitude(i,j)*cos_theta_latitude(i,j)
              clcp(i,j)=cos_theta_longitude(i,j)*cos_theta_latitude(i,j)
              slcp(i,j)=sin_theta_longitude(i,j)*cos_theta_latitude(i,j)
            END DO
          END DO
        END IF
        DO j=pdims%j_start,pdims%j_end
          DO i=pdims%i_start,pdims%i_end
            IF (qam_m1) am_m1(i,j)=0.0
            IF (qam_m2) am_m2(i,j)=0.0
            IF (qam_m3) am_m3(i,j)=0.0
            IF (qam_w1) am_w1(i,j)=0.0
            IF (qam_w2) am_w2(i,j)=0.0
            IF (qam_w3) am_w3(i,j)=0.0
          END DO
        END DO
        DO k=pdims%k_start,pdims%k_end
          DO j=pdims%j_start,pdims%j_end
            DO i=pdims%i_start,pdims%i_end
              r3dr=rho_dry(i,j,k)*r_rho_levels(i,j,k)                   &
     &          *delr_rho(i,j,k)*factor
              rocos=omega*r_rho_levels(i,j,k)*cos_theta_latitude(i,j)

              IF (qam_w1) am_w1(i,j)=am_w1(i,j)                         &
     &          +(field_rho(i,j,k,irho_u)*spclcp(i,j)                   &
     &           -field_rho(i,j,k,irho_v)*slcp(i,j))*r3dr
              IF (qam_w2) am_w2(i,j)=am_w2(i,j)                         &
     &          +(field_rho(i,j,k,irho_u)*spslcp(i,j)                   &
     &           +field_rho(i,j,k,irho_v)*clcp(i,j))*r3dr
              IF (qam_w3) am_w3(i,j)=                                   &
     &          am_w3(i,j)-field_rho(i,j,k,irho_u)*cossq(i,j)*r3dr

              IF (qam_m1) am_m1(i,j)=am_m1(i,j)+rocos*spclcp(i,j)*r3dr
              IF (qam_m2) am_m2(i,j)=am_m2(i,j)+rocos*spslcp(i,j)*r3dr
              IF (qam_m3) am_m3(i,j)=am_m3(i,j)-rocos*cossq(i,j)*r3dr
            END DO
          END DO
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH items 417-418:Pstar
! ----------------------------------------------------------------------
      IF (qpstar) THEN
        DO j=pdims%j_start,pdims%j_end
          DO i=pdims%i_start,pdims%i_end
            o_pstar(i,j)=pstar(i,j)
          END DO
        END DO
      END IF

! pstar at uv-points only in the New Dynamics

      IF (qpstar_uv) THEN
        DO j=vdims%j_start,vdims%j_end
          DO i=udims%i_start,udims%i_end
            o_pstar_uv(i,j)=pstar_uv(i,j)
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 419: Energy Correction
! ----------------------------------------------------------------------
      IF (qencorr) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            encorr(i,j)=(cp-r)*energy_corr_now*                         &
     &        vert_int_array(i,j,ip_dry_mass)*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 420:Total column cvT
! ----------------------------------------------------------------------
      IF (qtot_cvt) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tot_cvt(i,j)=vert_int_array(i,j,ip_cvt)*grid_factor
          END DO
        END DO

      END IF
! ----------------------------------------------------------------------
! STASH item 421:Total column gr
! ----------------------------------------------------------------------
      IF (qtot_gr) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            tot_gr(i,j)=vert_int_array(i,j,ip_gr)*grid_factor
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 422: Column integral ugz dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_ugz) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_ugz(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start , pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_ugz(i,j)= col_ugz(i,j)+field_rho(i,j,k,irho_u)*g*     &
     &                     field_rho(i,j,k,irho_z)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 423: Column integral vgz dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_vgz) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_vgz(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_vgz(i,j)= col_vgz(i,j)+field_rho(i,j,k,irho_v)*g*     &
     &                     field_rho(i,j,k,irho_z)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 424: Column integral wgz dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_wgz) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_wgz(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_wgz(i,j)= col_wgz(i,j)+field_rho(i,j,k,irho_w)*g*     &
     &                     field_rho(i,j,k,irho_z)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 425: Column integral uT dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_ut) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_ut(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_ut(i,j)=col_ut(i,j)+ field_rho(i,j,k,irho_u)*         &
     &                     field_rho(i,j,k,irho_t)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 426: Column integral vt dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_vt) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_vt(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_vt(i,j)=col_vt(i,j)+ field_rho(i,j,k,irho_v)*         &
     &                     field_rho(i,j,k,irho_t)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 427: Column integral wt dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_wt) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_wt(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_wt(i,j)=col_wt(i,j)+ field_rho(i,j,k,irho_w)*         &
     &                     field_rho(i,j,k,irho_t)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 428: Column integral uq dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_uq) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_uq(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_uq(i,j)= col_uq(i,j)+field_rho(i,j,k,irho_u)*         &
     &                     field_rho(i,j,k,irho_q)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 429: Column integral vq dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_vq) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_vq(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_vq(i,j)= col_vq(i,j)+field_rho(i,j,k,irho_v)*         &
     &                     field_rho(i,j,k,irho_q)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 430: Column integral wq dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_wq) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_wq(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_wq(i,j)= col_wq(i,j)+ field_rho(i,j,k,irho_w)*        &
     &                     field_rho(i,j,k,irho_q)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 431: Column integral uv dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_uv) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_uv(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_uv(i,j)= col_uv(i,j)+ field_rho(i,j,k,irho_u)*        &
     &                     field_rho(i,j,k,irho_v)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 432: Column integral uw dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_uw) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_uw(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_uw(i,j)= col_uw(i,j)+ field_rho(i,j,k,irho_u)*        &
     &                     field_rho(i,j,k,irho_w)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 433: Column integral vw dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_vw) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_vw(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_vw(i,j)=col_vw(i,j)+ field_rho(i,j,k,irho_v)*         &
     &                     field_rho(i,j,k,irho_w)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 434: Column integral uke dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_uke) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_uke(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_uke(i,j)=col_uke(i,j)+ field_rho(i,j,k,irho_u)*       &
     &                     field_rho(i,j,k,irho_ke)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 435: Column integral vke dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_vke) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_vke(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_vke(i,j)= col_vke(i,j)+field_rho(i,j,k,irho_v)*       &
     &                     field_rho(i,j,k,irho_ke)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 436: Column integral wke dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_wke) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_wke(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_wke(i,j)=col_wke(i,j)+field_rho(i,j,k,irho_w)*        &
     &                     field_rho(i,j,k,irho_ke)*mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 437: Column integral u dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_u) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_u(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_u(i,j)=col_u(i,j)+field_rho(i,j,k,irho_u)             &
     &                                         *mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 438: Column integral v dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_v) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_v(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_v(i,j)=col_v(i,j)+field_rho(i,j,k,irho_v)             &
     &                                         *mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 439: Column integral w dry mass weighted
! ----------------------------------------------------------------------
      IF (qcol_w) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            col_w(i,j)= 0.0
          END DO
        END DO
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              col_w(i,j)= col_w(i,j)+ field_rho(i,j,k,irho_w)           &
     &                                         *mass_we(i,j,k)
            END DO
          END DO
        END DO
      END IF
! ----------------------------------------------------------------------
!  Stash items 112 : wbig w > wbig on model levels
! ----------------------------------------------------------------------
      IF(qwbig_m) THEN
        wbig=1.0     ! test value
        counter=1
        DO k = 1, wdims%k_end     ! loop over levels
          IF (wbig_m_list(k)) THEN
             DO j = wdims%j_start, wdims%j_end
               DO i = wdims%i_start, wdims%i_end
                 IF (w(i,j,k) >  wbig) THEN
                    wbig_m(i,j,counter) = 1.0
                 ELSE
                    wbig_m(i,j,counter) = 0.0
                 END IF
               END DO
             END DO
             counter=counter+1
          END IF
        END DO
      END IF
! ----------------------------------------------------------------------
!  Stash items 114 : wbig2 w > wbig on model levels
! ----------------------------------------------------------------------
      IF(qwbig2_m) THEN
        wbig2=0.1     ! test value
        counter=1
        DO k = 1, wdims%k_end     ! loop over levels
          IF (wbig2_m_list(k)) THEN
             DO j = wdims%j_start, wdims%j_end
               DO i = wdims%i_start, wdims%i_end
                 IF (w(i,j,k) >  wbig2) THEN
                    wbig2_m(i,j,counter) = 1.0
                 ELSE
                    wbig2_m(i,j,counter) = 0.0
                 END IF
               END DO
             END DO
             counter=counter+1
          END IF
        END DO
      END IF
! ----------------------------------------------------------------------
!  Stash items 113 : RH on model levels - see earlier calculation
!                    Note this is at theta points
! ----------------------------------------------------------------------
      IF(qRH_m) THEN
        counter=1
        DO k = 1, qdims%k_end     ! loop over levels
          IF (RH_m_list(k)) THEN
             DO j = qdims%j_start, qdims%j_end
               DO i = qdims%i_start, qdims%i_end
                   RH_m(i,j,counter) = RH(i,j,k)
               END DO
             END DO
             counter=counter+1
          END IF
        END DO
      END IF
! ----------------------------------------------------------------------
! STASH item 451,452,453,454: Tropopause diagnostics
! ----------------------------------------------------------------------
      IF (qneed_tropopause) THEN
!
! DEPENDS ON: tropin
        CALL tropin (T,exner_rho_levels,                                &
                       exner_theta_levels(:,:,1:tdims%k_end),           &
                       pdims%i_end, pdims%j_end, pdims%k_end,           &
                       offx, offy,                                      &
                       at_extremity,scm_dummy_1d,scm_dummy_2d,          &
                       min_trop_level, max_trop_level, trindx )
!
!           Set the tropopause height at the grid-point.
!           To calculate at the next highest theta level, use
!              r_theta_levels(i, j, trindx(i, j)).
!           To calculate at the model layer boundary, use
!              r_rho_levels(i, j, trindx(i, j)).
!           Note that as TRINDX
!           is counted upward from the surface, while the first rho
!           level is omitted from the physics, no offset to the final
!           index of r_rho_levels is required.
!           (See Deck O3EXP1.301--304)
!           Values of r are measured from the centre of the Earth.
!           Hence subtraction of the surface r value must be performed.
        IF (QTROP_Z) THEN
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              trop_z(i,j) =   r_rho_levels(i, j, trindx(i, j))          &
     &                      - r_theta_levels(i, j, 0)
            END DO
          END DO
        END IF
!
!           Set the tropopause temperature at the grid-point.
!           To calculate at the next highest theta level, set
!              trop_t(i,j)=T(i, j, trindx(i, j))
!           Use linear interpolation in r
!              (consistent with Energy diagnostics).
!           Note that as TRINDX
!           is counted upward from the surface, while the first rho
!           level is omitted from the physics, no offset to the final
!           index of r_rho_levels is required.
!           (See Deck O3EXP1.301--304)
!
        IF (QTROP_T) THEN
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
                k = trindx(i,j)
                weight1=r_theta_levels(i,j,k)-r_rho_levels(i,j,k)
                weight2=r_rho_levels(i,j,k)-r_theta_levels(i,j,k-1)
                weight3=r_theta_levels(i,j,k)-r_theta_levels(i,j,k-1)

                ww1 = weight1/weight3
                ww2 = weight2/weight3

                trop_t(i,j) = ww2 * T(i,j,k) + ww1 * T(i,j,k-1)

            END DO
          END DO
        END IF

!           Set the tropopause pressure at the grid-point.
!           To calculate at the next highest theta level, use
!              p_theta_levels(i, j, trindx(i, j)).
!           To calculate at the model layer boundary, use
!              p(i, j, trindx(i, j)).
!           Note that as TRINDX
!           is counted upward from the surface, while the first rho
!           level is omitted from the physics, no offset to the final
!           index of r_rho_levels is required.
!           (See Deck O3EXP1.301--304)
!
        IF (QTROP_P) THEN
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              trop_p(i,j)=p(i, j, trindx(i, j))
            END DO
          END DO
        END IF

!       Interim initialization of qmodel_half_height
        QMODEL_HALF_HEIGHT = .TRUE.
!
        IF(QMODEL_HALF_HEIGHT) THEN
!-----------------------------------------------------------------------
!      Call tropopause program with quality control
!      Needs CALL to TROP here
!-----------------------------------------------------------------------
!         Calculates ICAO heights from Trop pressure field
!-----------------------------------------------------------------------
          IF(QTROP_ICAO) THEN
            DO j = pdims%j_start, pdims%j_end
              DO i = pdims%i_start, pdims%i_end
                trop_icao(i,j)=13840. ! metres
              END DO
            END DO
! Interim values in the absence of subroutine ICAO_HT
          END IF
!
        ELSE
          WRITE(6,444)
          WRITE(6,444)
 444      FORMAT(' Subroutine TROP not called No MODEL_HALF_HEIGHTS')
        END IF
      END IF

! ---------------------------------------------------------------------
! Column saturation fraction   442 
! ---------------------------------------------------------------------
      IF (qcol_sat1) THEN

! Need to interpolate rho to theta levels, q, qsat and T on theta levels
        DO k = pdims%k_start, pdims%k_end
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              rho_only(i,j,k) = rho(i,j,k)                              &
                              /(r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
            END DO        
          END DO        
        END DO        
        DO k = 1, tdims%k_end - 1
          DO j = tdims%j_start, tdims%j_end
            DO I = tdims%i_start, tdims%i_end
              weight1 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
              weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
              weight3 = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
              ww1 = weight1/weight3
              ww2 = weight2/weight3
              rho_theta(i,j,k) = ww2 * rho_only(i,j,k+1)                &
                               + ww1 * rho_only(i,j,k)
            END DO
          END DO
        END DO
        k = tdims%k_end       !  Top theta level above top rho level
          DO j = tdims%j_start, tdims%j_end
            DO I = tdims%i_start, tdims%i_end
              rho_theta(i,j,k) = rho_only(i,j,k)
            END DO
          END DO

        work_1(:,:)   = 0.0   ! initialise to zero
        col_sat1(:,:) = 0.0

       ! integrate fields
        k = 1                ! special case
          DO j = qdims%j_start, qdims%j_end
            DO I = qdims%i_start, qdims%i_end
              temp_mass = rho_theta(i,j,k)                              &
                       * r_theta_levels(i,j,k)*r_theta_levels(i,j,k)    &
                       * (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,0))
              work_1(i,j) = work_1(i,j) + qsat_work(i,j,k)*temp_mass
              col_sat1(i,j) = col_sat1(i,j) + q(i,j,k)*temp_mass
            END DO
          END DO
        IF (qdims%k_end == tdims%k_end) THEN ! wet_levels == model_levels
          DO k = 2, qdims%k_end - 1
            DO j = qdims%j_start, qdims%j_end
              DO I = qdims%i_start, qdims%i_end
                temp_mass = rho_theta(i,j,k)                            &
                         * r_theta_levels(i,j,k)*r_theta_levels(i,j,k)  &
                         * (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
                work_1(i,j) = work_1(i,j) + qsat_work(i,j,k)*temp_mass
                col_sat1(i,j) = col_sat1(i,j) + q(i,j,k)*temp_mass

              END DO
            END DO
          END DO
          k = qdims%k_end        ! special case 
            DO j = qdims%j_start, qdims%j_end
              DO I = qdims%i_start, qdims%i_end
                temp_mass = rho_theta(i,j,k)                            &
                         * r_theta_levels(i,j,k)*r_theta_levels(i,j,k)  &
                         * (r_theta_levels(i,j,k) - r_rho_levels(i,j,k))
                work_1(i,j) = work_1(i,j) + qsat_work(i,j,k)*temp_mass
                col_sat1(i,j) = col_sat1(i,j) + q(i,j,k)*temp_mass
              END DO
            END DO

        ELSE
          DO k = 2, qdims%k_end
            DO j = qdims%j_start, qdims%j_end
              DO I = qdims%i_start, qdims%i_end
                temp_mass = rho_theta(i,j,k)                            &
                         * r_theta_levels(i,j,k)*r_theta_levels(i,j,k)  &
                         * (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
                work_1(i,j) = work_1(i,j) + qsat_work(i,j,k)*temp_mass
                col_sat1(i,j) = col_sat1(i,j) + q(i,j,k)*temp_mass
              END DO
            END DO
          END DO
        END IF 
        DO j = qdims%j_start, qdims%j_end
          DO I = qdims%i_start, qdims%i_end
             col_sat1(i,j) = col_sat1(i,j) / work_1(i,j)
          END DO
        END DO

      END IF      ! qcol_sat1
! ---------------------------------------------------------------------
! Check error condition
      IF(ErrorStatus >  0) THEN

         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      END IF

      IF (lhook) CALL dr_hook('EOT_DIAG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Eot_diag

END MODULE eot_diag_mod
