! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ ---------------------------------------------------------------------
!  Sructure containing new boundary layer diagnostics.
!  This permits easier addition of new boundary layer
!  diagnostics without additional passing of arguments
!  though the boundary layer tree.
!  It also does not require the addition
!  of extra subroutine arguments when adding a new diagnostic.
!- ----------------------------------------------------------------------

MODULE bl_diags_mod

IMPLICIT NONE
  TYPE strnewbldiag

! Need to create a flag and a pointer
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer

    LOGICAL :: l_smltop
    LOGICAL :: l_dsctop
    LOGICAL :: l_zhlocal
    LOGICAL :: l_zhpar
    LOGICAL :: l_dscbase
    LOGICAL :: l_cldbase
    LOGICAL :: l_weparm
    LOGICAL :: l_weparm_dsc
    LOGICAL :: l_oblen
    LOGICAL :: l_ustar
    LOGICAL :: l_wbsurf
    LOGICAL :: l_gradrich
    LOGICAL :: l_wstar
    LOGICAL :: l_dbdz
    LOGICAL :: l_dvdzm
    LOGICAL :: l_rhokm
    LOGICAL :: l_rhokh
    LOGICAL :: l_tke
    LOGICAL :: l_ostressx
    LOGICAL :: l_ostressy
    LOGICAL :: l_dtfric
    LOGICAL :: l_rhogamu
    LOGICAL :: l_rhogamv
    LOGICAL :: l_rhogamt
    LOGICAL :: l_rhogamq
    LOGICAL :: l_elm
    LOGICAL :: l_tke_shr_prod
    LOGICAL :: l_tke_boy_prod
    LOGICAL :: l_tke_dissp
    LOGICAL :: l_sm
    LOGICAL :: l_sh
    LOGICAL :: l_wb_ng
    LOGICAL :: l_cf_trb
    LOGICAL :: l_ql_trb
    LOGICAL :: l_sgm_trb
    LOGICAL :: l_elm3D
    LOGICAL :: l_elh3D
    LOGICAL :: l_rhoKmloc
    LOGICAL :: l_rhoKhloc
    LOGICAL :: l_rhoKmsurf
    LOGICAL :: l_rhoKhsurf
    LOGICAL :: l_rhoKmSc
    LOGICAL :: l_rhoKhSc
    LOGICAL :: l_weight1d

    REAL, POINTER :: smltop(:, :)
!                    356      Top of surface mixed layer
    REAL, POINTER :: dsctop(:, :)
!                    357      Top of decoupled stratocu layer
    REAL, POINTER :: zhlocal(:, :)
!                    358      BL depth diagnosed from Ri>RiCrit
    REAL, POINTER :: zhpar(:, :)
!                    359      Height of diagnosis parcel top
    REAL, POINTER :: dscbase(:, :)
!                    360      Height of decoupled layer base
    REAL, POINTER :: cldbase(:, :)
!                    361      Height of stratocumulus cloud base
    REAL, POINTER :: weparm(:, :)
!                    362      Entrainment rate for SML
    REAL, POINTER :: weparm_dsc(:, :)
!                    363      Entrainment rate for DSC
    REAL, POINTER :: oblen(:, :)
!                    464      Surface Obukhov length
    REAL, POINTER :: ustar(:, :)
!                    465      Friction velocity
    REAL, POINTER :: wbsurf(:, :)
!                    467      Surface buoyancy flux
    REAL, POINTER :: gradrich(:, :, :)
!                    468      Gradient Richardson number
    REAL, POINTER :: wstar(:, :)
!                    466      Convective velocity scale
    REAL, POINTER :: dbdz(:, :, :)
!                    469      Vertical buoyancy gradient
    REAL, POINTER :: dvdzm(:, :, :)
!                    470      Modulus of wind shear
    REAL, POINTER :: rhokm(:, :, :)
!                    471      BL Momentum diffusivity
    REAL, POINTER :: rhokh(:, :, :)
!                    472      BL Thermal diffusivity
    REAL, POINTER :: tke(:, :, :)
!                    473      Turbulent kinetic energy
    REAL, POINTER :: ostressx(:, :, :)
!                    474      Orographic stress (x-component)
    REAL, POINTER :: ostressy(:, :, :)
!                    475      Orographic stress (y-component)
    REAL, POINTER :: elm3D(:, :, :)
!                    501      Mixing length for momentum 
    REAL, POINTER :: elh3D(:, :, :)
!                    502      Mixing length for heat and moisture
    REAL, POINTER :: rhoKmloc(:, :, :)
!                    503      Km diffusion coeff from local scheme
    REAL, POINTER :: rhoKhloc(:, :, :)
!                    504      Kh diffusion coeff from local scheme
    REAL, POINTER :: rhoKmsurf(:, :, :)
!                    505      Km diffusion coeff for surface-driven turb
    REAL, POINTER :: rhoKhsurf(:, :, :)
!                    506      Kh diffusion coeff for surface-driven turb
    REAL, POINTER :: rhoKmSc(:, :, :)
!                    507      Km diffusion coeff for  cloud-top-driven turb
    REAL, POINTER :: rhoKhSc(:, :, :)
!                    508      Kh diffusion coeff for  cloud-top-driven turb
    REAL, POINTER :: weight1d(:, :, :)
!                    509      weighting applied to 1D BL scheme in Smag blending
    REAL, POINTER :: dTfric(:, :, :)
!                    188      Heating increment from turbulence dissipation
    REAL, POINTER :: rhogamu(:, :, :)
!                    130      counter gradient term of taux
    REAL, POINTER :: rhogamv(:, :, :)
!                    131      counter gradient term of tauy
    REAL, POINTER :: rhogamt(:, :, :)
!                    132      counter gradient term of ftl
    REAL, POINTER :: rhogamq(:, :, :)
!                    133      counter gradient term of fqw
    REAL, POINTER :: elm(:, :, :)
!                    134      mixing length
    REAL, POINTER :: tke_shr_prod(:, :, :)
!                    135      production rate of TKE by shear
    REAL, POINTER :: tke_boy_prod(:, :, :)
!                    136      production rate of TKE by buoyancy
    REAL, POINTER :: tke_dissp(:, :, :)
!                    137      dissipation rate of TKE
    REAL, POINTER :: sm(:, :, :)
!                    138      non-dimensional diffusion coef. for u, v
    REAL, POINTER :: sh(:, :, :)
!                    139      non-dimensional diffusion coef. for t, q
    REAL, POINTER :: wb_ng(:, :, :)
!                    140      non-gradient buoyancy flux
    REAL, POINTER :: cf_trb(:, :, :)
!                    141      cloud fraction used in the TKE schemes
    REAL, POINTER :: ql_trb(:, :, :)
!                    142      condensed water used in the TKE schemes
    REAL, POINTER :: sgm_trb(:, :, :)
!                    143      standard deviation of the distribution
!                             function in the TKE schemes

  END TYPE strnewbldiag
! ----------------------------------------------------------------------
END MODULE bl_diags_mod
