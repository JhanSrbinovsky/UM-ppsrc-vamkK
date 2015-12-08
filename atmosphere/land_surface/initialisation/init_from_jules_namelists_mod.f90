! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Initialise model from JULES namelist data
!
MODULE init_from_jules_namelists_mod

IMPLICIT NONE

! Description:
!  init_jules_from_namelists initialises the relevant JULES arrays
!  from data read in via namelists.
!
! Method:
!  Arrays are allocated in the first call, then initialised from namelist data.
!  Should only be called after all the JULES namelists have been read in.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Land
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.

CONTAINS

SUBROUTINE init_from_jules_namelists (land_field, ntiles, sm_levels,   &
                                      nice, nice_use)

USE nstypes, ONLY:                                                     &
  npft,            nnvg

USE switches, ONLY:                                                    &
! jules_triffid:
  l_triffid,       l_phenol,                                           &
! jules_snow_param:
  can_model,                                                           &
! multiple uses:
  l_aggregate

USE ancil_info, ONLY:                                                  &
  nsmax

USE pftparm_io, ONLY:                                                  &
! namelist variables:
  c3_io,           orient_io,        a_wl_io,                          &
  a_ws_io,         albsnc_max_io,    albsnc_min_io,                    &
  albsnf_maxu_io,  albsnf_max_io,    albsnf_maxl_io,                   &
  alpha_io,        alniru_io,        alnir_io,                         &
  alnirl_io,       alparu_io,        alpar_io,                         &
  alparl_io,       b_wl_io,          catch0_io,                        &
  dcatch_dlai_io,  dgl_dm_io,        dgl_dt_io,                        &
  dqcrit_io,       dz0v_dh_io,       eta_sl_io,                        &
  fd_io,           fsmc_of_io,       f0_io,                            &
  g_leaf_0_io,     glmin_io,         infil_f_io,                       &
  kext_io,         kpar_io,          neff_io,                          &
  nl0_io,          nr_nl_io,         ns_nl_io,                         &
  omegau_io,       omega_io,         omegal_io,                        &
  omniru_io,       omnir_io,         omnirl_io,                        &
  r_grow_io,       rootd_ft_io,      sigl_io,                          &
  tleaf_of_io,     tlow_io,          tupp_io,                          &
  emis_pft_io,     z0hm_pft_io,      z0hm_classic_pft_io,              &
  dust_veg_scj_io, fl_o3_ct_io,      dfp_dcuo_io,                      &
  ief_io,          tef_io,           mef_io,                           &
  aef_io

USE pftparm, ONLY:                                                     &
  c3,              orient,           a_wl,                             &
  a_ws,            albsnc_max,       albsnc_min,                       &
  albsnf_maxu,     albsnf_max,       albsnf_maxl,                      &
  alpha,           alniru,           alnir,                            &
  alnirl,          alparu,           alpar,                            &
  alparl,          b_wl,             catch0,                           &
  dcatch_dlai,     dgl_dm,           dgl_dt,                           &
  dqcrit,          dz0v_dh,          eta_sl,                           &
  fd,              fsmc_of,          f0,                               &
  g_leaf_0,        glmin,            infil_f,                          &
  kext,            kpar,             neff,                             &
  nl0,             nr_nl,            ns_nl,                            &
  omegau,          omega,            omegal,                           &
  omniru,          omnir,            omnirl,                           &
  r_grow,          rootd_ft,         sigl,                             &
  tleaf_of,        tlow,             tupp,                             &
  emis_pft,        z0hm_pft,         z0hm_classic_pft,                 &
  dust_veg_scj,    fl_o3_ct,         dfp_dcuo,                         &
  ief,             tef,              mef,                              &
  aef

USE nvegparm_io, ONLY:                                                 &
! namelist variables:
  albsnc_nvg_io,   albsnf_nvgu_io,  albsnf_nvg_io,                     & 
  albsnf_nvgl_io,  catch_nvg_io,    gs_nvg_io,                         &
  infil_nvg_io,    z0_nvg_io,       ch_nvg_io,                         &
  vf_nvg_io,       emis_nvg_io,     z0hm_nvg_io,                       &
  z0hm_classic_nvg_io

USE nvegparm, ONLY:                                                    &
  albsnc_nvg,      albsnf_nvgu,     albsnf_nvg,                        &
  albsnf_nvgl,     catch_nvg,       gs_nvg,                            &
  infil_nvg,       z0_nvg,          ch_nvg,                            &
  vf_nvg,          emis_nvg,        z0hm_nvg ,                         &
  z0hm_classic_nvg

! For use with jules_pftparm and jules_nvegparm
USE c_z0h_z0m, ONLY:                                                   &
  z0h_z0m,         z0h_z0m_classic

USE trif_io, ONLY:                                                     &
! namelist variables:
  crop_io,         g_area_io,        g_grow_io,                        &
  g_root_io,       g_wood_io,        lai_max_io,                       &
  lai_min_io

USE trif, ONLY:                                                        &
  crop,            g_area,           g_grow,                           &
  g_root,          g_wood,           lai_max,                          &
  lai_min

USE snow_param, ONLY:                                                  &
! namelist variables:
  dzsnow_io,       cansnowpft,                                         &
! non-namelist variables:
  dzsnow,          cansnowtile

USE soil_param, ONLY:                                                  &
! namelist variables:
  dzsoil_io,                                                           &
! non-namelist variables:
  dzsoil

USE surf_param, ONLY:                                                  &
! non-namelist variables:
  diff_frac

USE c_elevate, ONLY:                                                   &
! namelist variables:
  surf_hgt_io,                                                         &
! non-namelist variables:
  surf_hgt

USE rad_param, ONLY:                                                   &
! namelist variables:
  dtland,          kland_numerator,                                    &
! non-namelist variables:
  kland,           tcland

! For use with jules_rad_param
USE c_0_dg_c, ONLY:                                                    &
  tm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: land_field   ! Number of land points
INTEGER, INTENT(IN) :: ntiles       ! Number of surface tiles
INTEGER, INTENT(IN) :: sm_levels    ! Number of soil layers
INTEGER, INTENT(IN) :: nice         ! Number of sea ice categories
INTEGER, INTENT(IN) :: nice_use     ! Number of sea ice categories used
                                    ! fully in surface exch and radiation

! Local variables
INTEGER             :: i            ! Looper 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('INIT_FROM_JULES_NAMELISTS',zhook_in,zhook_handle)

! Allocate JULES arrays
! DEPENDS ON: allocate_jules_arrays
CALL allocate_jules_arrays( land_field, ntiles, sm_levels, nice, nice_use )

! Begin initialising arrays from:
! jules_pftparm
c3(:)           = c3_io(1:npft)
orient(:)       = orient_io(1:npft)
a_wl(:)         = a_wl_io(1:npft)
a_ws(:)         = a_ws_io(1:npft)
albsnc_max(:)   = albsnc_max_io(1:npft)
albsnc_min(:)   = albsnc_min_io(1:npft)
albsnf_maxu(:)  = albsnf_maxu_io(1:npft)
albsnf_max(:)   = albsnf_max_io(1:npft)
albsnf_maxl(:)  = albsnf_maxl_io(1:npft)
alpha(:)        = alpha_io(1:npft)
alniru(:)       = alniru_io(1:npft) 
alnir(:)        = alnir_io(1:npft)
alnirl(:)       = alnirl_io(1:npft) 
alparu(:)       = alparu_io(1:npft) 
alpar(:)        = alpar_io(1:npft)
alparl(:)       = alparl_io(1:npft) 
b_wl(:)         = b_wl_io(1:npft)
catch0(:)       = catch0_io(1:npft)
dcatch_dlai(:)  = dcatch_dlai_io(1:npft)
dgl_dm(:)       = dgl_dm_io(1:npft)
dgl_dt(:)       = dgl_dt_io(1:npft)
dqcrit(:)       = dqcrit_io(1:npft)
dz0v_dh(:)      = dz0v_dh_io(1:npft)
eta_sl(:)       = eta_sl_io(1:npft)
fd(:)           = fd_io(1:npft)
fsmc_of(:)      = fsmc_of_io(1:npft)
f0(:)           = f0_io(1:npft)
g_leaf_0(:)     = g_leaf_0_io(1:npft)
glmin(:)        = glmin_io(1:npft)
infil_f(:)      = infil_f_io(1:npft)
kext(:)         = kext_io(1:npft)
kpar(:)         = kpar_io(1:npft)
neff(:)         = neff_io(1:npft)
nl0(:)          = nl0_io(1:npft)
nr_nl(:)        = nr_nl_io(1:npft)
ns_nl(:)        = ns_nl_io(1:npft)
omegau(:)       = omegau_io(1:npft) 
omega(:)        = omega_io(1:npft)
omegal(:)       = omegal_io(1:npft) 
omniru(:)       = omniru_io(1:npft)
omnir(:)        = omnir_io(1:npft)
omnirl(:)       = omnirl_io(1:npft)
r_grow(:)       = r_grow_io(1:npft)
rootd_ft(:)     = rootd_ft_io(1:npft)
sigl(:)         = sigl_io(1:npft)
tleaf_of(:)     = tleaf_of_io(1:npft)
tlow(:)         = tlow_io(1:npft)
tupp(:)         = tupp_io(1:npft)
emis_pft(:)     = emis_pft_io(1:npft)
z0hm_pft(:)     = z0hm_pft_io(1:npft)
z0hm_classic_pft(:)                                                    & 
                = z0hm_classic_pft_io(1:npft)
dust_veg_scj(:) = dust_veg_scj_io(1:npft)
fl_o3_ct(:)     = fl_o3_ct_io(1:npft)
dfp_dcuo(:)     = dfp_dcuo_io(1:npft)
ief(:)          = ief_io(1:npft)
tef(:)          = tef_io(1:npft)
mef(:)          = mef_io(1:npft)
aef(:)          = aef_io(1:npft)

z0h_z0m(1:npft) = z0hm_pft(1:npft) 
z0h_z0m_classic(1:npft)                                                &
                = z0hm_classic_pft(1:npft)       

! jules_nvegparm
albsnc_nvg(:) = albsnc_nvg_io(1:nnvg)
albsnf_nvgu(:)= albsnf_nvgu_io(1:nnvg) 
albsnf_nvg(:) = albsnf_nvg_io(1:nnvg)
albsnf_nvgl(:)= albsnf_nvgl_io(1:nnvg) 
catch_nvg(:)  = catch_nvg_io(1:nnvg)
gs_nvg(:)     = gs_nvg_io(1:nnvg)
infil_nvg(:)  = infil_nvg_io(1:nnvg)
z0_nvg(:)     = z0_nvg_io(1:nnvg)
ch_nvg(:)     = ch_nvg_io(1:nnvg)
vf_nvg(:)     = vf_nvg_io(1:nnvg)
emis_nvg(:)   = emis_nvg_io(1:nnvg)
z0hm_nvg(:)   = z0hm_nvg_io(1:nnvg)
z0hm_classic_nvg(:)                                                    &
              = z0hm_classic_nvg_io(1:nnvg)

z0h_z0m(npft+1:npft+nnvg) = z0hm_nvg(1:nnvg)
z0h_z0m_classic(npft+1:npft+nnvg)                                      &
                = z0hm_classic_nvg(1:nnvg)  


! jules_triffid
! Space only allocated if at least phenology is enabled
IF ( l_triffid .OR. l_phenol ) THEN
  crop(:)    = crop_io(1:npft)
  g_area(:)  = g_area_io(1:npft)
  g_grow(:)  = g_grow_io(1:npft)
  g_root(:)  = g_root_io(1:npft)
  g_wood(:)  = g_wood_io(1:npft)
  lai_max(:) = lai_max_io(1:npft)
  lai_min(:) = lai_min_io(1:npft)
END IF


! jules_snow_param
! dzsnow has only been allocated if nsmax > 0
IF ( nsmax > 0 ) THEN
  dzsnow(:) = dzsnow_io(1:nsmax)
END IF

! Set up canSnowTile
cansnowtile(:) = .FALSE.
IF ( .NOT. l_aggregate .AND. can_model == 4 ) THEN
  cansnowtile(1:npft) = cansnowpft(1:npft)
END IF


! jules_soil_param
dzsoil(:) = dzsoil_io(1:sm_levels)


! jules_surf_param
diff_frac(:) = 0.0


! jules_elevate
! Use same height for all land points.
IF ( l_aggregate ) THEN
  surf_hgt(:,:) = 0.0
ELSE
  DO i=1,ntiles
    surf_hgt(:,i) = surf_hgt_io(i)
  END DO
END IF


! jules_rad_param
! Calculate kland and tcland from namelist values
kland  = kland_numerator / dtland
tcland = tm - dtland


IF (lhook) CALL dr_hook('INIT_FROM_JULES_NAMELISTS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_from_jules_namelists

END MODULE init_from_jules_namelists_mod
