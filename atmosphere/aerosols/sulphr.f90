! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE sulphr_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE sulphr(                                                 &
  ! Arguments IN
  halo_i, halo_j, off_x, off_y,                                    &
  row_length, rows,                                                &
  model_levels, wet_model_levels,                                  &
  theta_field_size,                                                &
  tstep,                                                           &
  cloudf, cosza2d,                                                 &
  p, t, q, qcl, qcf,                                               &
  oh_conc, h2o2, ho2_conc, o3,                                     &
  n_droplet,                                                       &
  l_sulpc_dms, l_sulpc_newdms,                                     &
  l_sulpc_ozone, l_sulpc_so2_o3_nonbuffered,                       &
  l_sulpc_nh3,                                                     &
  l_sulpc_online_oxidants, l_sulpc_2_way_coupling,                 &
  ! Arguments IN/OUT
  so2, dms, so4_ait, so4_acc, so4_dis,                             &
  nh3, h2o2_mxr,                                                   &
  ! Arguments OUT
  msa, nh3_dep,                                                    &
  f_dms_to_so2, f_dms_to_so4, f_dms_to_msa,                        &
  deltas_dry, deltas_wet, deltas_wet_o3,                           &
  deltas_tot, deltas_dms,                                          &
  deltas_evap, deltas_nucl, deltas_diffuse,                        &
  deltas_merge, deltas_coag, psi )

!---------------------------------------------------------------------
! Purpose: To perform oxidation chemistry of sulphur dioxide to 3 modes
!          of Sulphate aerosol (Aitken, accumulation and dissolved),
!          and dimethyl sulphide to sulphur dioxide and methyl sulphonic
!          acid.
!          There is also exchange between the 3 modes of sulphate
!          aerosol due to nucleation, diffusion and evaporation
!          processes, and Aitken-accumulation mode merging caused by
!          particle growth due to the condensation of newly-formed
!          sulphuric acid onto Aitken mode sulphate particles.
!          Called by Aero_Ctl

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols

! Code Description:
!  Language: Fortran 90.
!  This code is written to UMDP3 v8 programming standards

! Documentation: UM Doc 20

!---------------------------------------------------------------------
USE vectlib_mod, ONLY :                                                 &
  log_v,  rtor_v

USE earth_constants_mod, ONLY: g

USE atmos_constants_mod, ONLY: r



USE water_constants_mod, ONLY: rho_water

USE conversions_mod, ONLY: pi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE


INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle

! Arguments with intent IN:

INTEGER :: row_length           !no. of pts along a row
INTEGER :: rows                 !no. of rows
INTEGER :: model_levels         !no. of model levels
INTEGER :: wet_model_levels     !no. of wet model levels
INTEGER :: halo_i               !EW halo size
INTEGER :: halo_j               !NS halo size
INTEGER :: off_x
INTEGER :: off_y
INTEGER :: theta_field_size

! CONTAINS OTHER REQUIRED PARAMETERS
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_S,                                                       &
                                 ! relative molecular mass S kg/mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_HO2,                                                     &
                                 ! relative molecular mass HO2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_S    = 3.20E-2,                                    &
                                          ! kg/mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &           RMM_HO2  = 3.30E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------
!----------------COMDECK C_SULCOND---------------------------------
! Arrays and variables to be used in interpolation to determine
! dependence of partitioning of dry oxidised SO2 to Aitken and
! accumulation modes:
      REAL                                                              &
     &       REL_HUM_TABLE(101),                                        &
                                 ! Relative humidities between 0 and 1
     &       A_TABLE(101),                                              &
                              ! CdotN_Ait/Mtot_Ait at corresponding RH
     &       B_TABLE(101)     ! CdotN_Acc/Mtot_Acc at corresponding RH
      PARAMETER(REL_HUM_TABLE = (/                                      &
     &        0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,           &
     &        0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15,           &
     &        0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23,           &
     &        0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31,           &
     &        0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,           &
     &        0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47,           &
     &        0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55,           &
     &        0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63,           &
     &        0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71,           &
     &        0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79,           &
     &        0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87,           &
     &        0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95,           &
     &        0.96, 0.97, 0.98, 0.99, 1.00 /) )
      PARAMETER(a_table = (/                                            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.80313E+07,            &
     &   3.87275E+07, 3.94294E+07, 4.01369E+07, 4.08502E+07,            &
     &   4.15693E+07, 4.22940E+07, 4.30243E+07, 4.37602E+07,            &
     &   4.45019E+07, 4.52490E+07, 4.60017E+07, 4.67600E+07,            &
     &   4.75240E+07, 4.82933E+07, 4.90682E+07, 4.98487E+07,            &
     &   5.06346E+07, 5.14260E+07, 5.22229E+07, 5.30252E+07,            &
     &   5.38329E+07, 5.46460E+07, 5.54646E+07, 5.62885E+07,            &
     &   5.71177E+07, 5.79524E+07, 5.87925E+07, 5.96377E+07,            &
     &   6.04882E+07, 6.13442E+07, 6.22053E+07, 6.30717E+07,            &
     &   6.39433E+07, 6.48201E+07, 6.57019E+07, 6.65891E+07,            &
     &   6.74816E+07, 6.83793E+07, 6.92821E+07, 7.01899E+07,            &
     &   7.11028E+07, 7.20209E+07, 7.29444E+07, 7.38725E+07,            &
     &   7.48059E+07, 7.57442E+07, 7.66875E+07, 7.76360E+07,            &
     &   7.85893E+07, 7.95478E+07, 8.13327E+07, 8.33181E+07,            &
     &   8.55388E+07, 8.80381E+07, 9.08706E+07, 9.41064E+07,            &
     &   9.78351E+07, 1.02175E+08, 1.07284E+08, 1.13379E+08,            &
     &   1.20761E+08, 1.29866E+08, 1.41339E+08, 1.56177E+08,            &
     &   1.75991E+08, 2.03530E+08, 2.43830E+08, 3.06957E+08,            &
     &   4.15500E+08  /) )
      PARAMETER(b_table = (/                                            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.59822e+05,            &
     &   8.70240e+05, 8.80669e+05, 8.91104e+05, 9.01551e+05,            &
     &   9.12005e+05, 9.22465e+05, 9.32935e+05, 9.43409e+05,            &
     &   9.53896e+05, 9.64389e+05, 9.74885e+05, 9.85393e+05,            &
     &   9.95906e+05, 1.00643e+06, 1.01695e+06, 1.02748e+06,            &
     &   1.03802e+06, 1.04857e+06, 1.05912e+06, 1.06967e+06,            &
     &   1.08023e+06, 1.09080e+06, 1.10137e+06, 1.11195e+06,            &
     &   1.12253e+06, 1.13312e+06, 1.14371e+06, 1.15431e+06,            &
     &   1.16491e+06, 1.17552e+06, 1.18613e+06, 1.19675e+06,            &
     &   1.20737e+06, 1.21799e+06, 1.22862e+06, 1.23925e+06,            &
     &   1.24989e+06, 1.26053e+06, 1.27117e+06, 1.28182e+06,            &
     &   1.29247e+06, 1.30312e+06, 1.31378e+06, 1.32445e+06,            &
     &   1.33511e+06, 1.34578e+06, 1.35645e+06, 1.36712e+06,            &
     &   1.37780e+06, 1.38848e+06, 1.40824e+06, 1.43001e+06,            &
     &   1.45410e+06, 1.48093e+06, 1.51096e+06, 1.54481e+06,            &
     &   1.58324e+06, 1.62724e+06, 1.67810e+06, 1.73752e+06,            &
     &   1.80784e+06, 1.89228e+06, 1.99547e+06, 2.12424e+06,            &
     &   2.28913e+06, 2.50720e+06, 2.80779e+06, 3.24555e+06,            &
     &   3.93362e+06 /) )
!
      INTEGER                                                           &
     &        nearest_index(row_length, rows, model_levels)
!  Index of element of REL_HUM_TABLE
!  that is closest to actual RH(i,j,k).
!
      REAL                                                              &
     &       dRH,                                                       &
                   ! Interval between values in REL_HUM_TABLE
     &       A_ARRAY(row_length, rows, model_levels),                   &
! A_TABLE(nearest_index) * MMR of Aitken sulphate
     &       B_ARRAY(row_length, rows, model_levels)
! B_TABLE(nearest_index) * MMR of accumulation sulphate
!
! c_sulcond.h contains arrays and variables used in calculation of PSI.

REAL :: tstep                        !chemistry tstep: LE physics tstep

REAL :: cloudf(row_length,rows,wet_model_levels)   !cloud fraction (0-1)
!cloud liquid water
REAL :: qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                 &
            wet_model_levels)
!cloud frozen water
REAL :: qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                 &
            wet_model_levels)
!specific humidity
REAL :: q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,wet_model_levels)
REAL :: cosza2d(row_length,rows)                   !cos zenith angle
!pressure
REAL :: p(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)
REAL :: t(row_length,rows,model_levels)            !temperature
REAL :: oh_conc(row_length,rows,model_levels)      !OH concn
REAL :: ho2_conc(row_length,rows,model_levels)     !HO2 concn
! For prescribed oxidants, this is the H2O2 max limit. 
! For online oxidants, it is the H2O2 mmr from UKCA.
REAL :: h2o2(row_length,rows,model_levels)
REAL :: o3(row_length,rows,model_levels)            !O3 mmr
REAL :: n_droplet(row_length, rows, wet_model_levels) !Drop concentration

LOGICAL :: l_sulpc_dms          !T if DMS chemistry required
LOGICAL :: l_sulpc_newdms       !T if new DMS scheme req'd, else old scheme used
LOGICAL :: l_sulpc_ozone        !T if O3 oxidn required
! l_sulpc_so2_o3_nonbuffered T if SO2+O3 reaction is NOT to be buffered by NH3.
LOGICAL :: l_sulpc_so2_o3_nonbuffered
!
! l_sulpc_nh3 T if NH3 buffering used (always T if L_SULPC_OZONE is T)
LOGICAL :: l_sulpc_nh3
LOGICAL :: l_sulpc_online_oxidants  ! T if oxidants from UKCA are used
! l_sulpc_2_way_coupling T if coupling to UKCA is two way, 
! i.e. if depleted oxidants are passed back to UKCA.
LOGICAL :: l_sulpc_2_way_coupling

! Arguments with intent IN/OUT:

REAL :: so2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                    &
            model_levels)
REAL :: dms(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                    &
              model_levels)
REAL :: so4_ait(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
              model_levels)
REAL :: so4_acc(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
              model_levels)
REAL :: so4_dis(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
              model_levels)
REAL :: nh3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                    &
              model_levels)
REAL :: h2o2_mxr(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
              model_levels)

! Arguments with intent OUT (diagnostics):



REAL :: msa(row_length,rows,model_levels)           !mmr S in MSA
REAL :: f_dms_to_so2(row_length,rows,model_levels)  !frac oxid DMS to SO2
REAL :: f_dms_to_so4(row_length,rows,model_levels)  !frac oxid DMS to SO4
REAL :: f_dms_to_msa(row_length,rows,model_levels)  !frac oxid DMS to MSA
REAL :: nh3_dep(row_length,rows,model_levels)       !NH3 depleted
REAL :: deltas_dry(row_length,rows,model_levels)    !SO2 dry ox per ts
REAL :: deltas_wet(row_length,rows,model_levels)    !SO2 wet ox by H2O2
REAL :: deltas_wet_o3(row_length,rows,model_levels) !SO2 wet ox by O3
REAL :: deltas_tot(row_length,rows,model_levels)    !total SO2 ox per ts
REAL :: deltas_dms(row_length,rows,model_levels)    !DMS dry ox per ts
!SO4_DIS released by evapn of cloud droplets to SO4_ACC per ts
REAL :: deltas_evap(row_length,rows,model_levels)   
!SO4_ACC transfd by nucleation to SO4_DIS per ts
REAL :: deltas_nucl(row_length,rows,model_levels)
!SO4_AIT transfd to SO4_DIS by diffusion per ts
REAL :: deltas_diffuse(row_length,rows,model_levels)
!SO4_AIT transfd by coagulation to SO4_ACC per ts
REAL :: deltas_coag(row_length,rows,model_levels)
!total SO4_AIT transfd per ts
REAL :: deltas_tot_ait(row_length,rows,model_levels)
!fraction of dry ox SO2 becoming SO4_AIT
REAL :: psi(row_length,rows,model_levels)

! Local variables:

INTEGER :: i,j,k               !loop variables

REAL :: dryrate                !rate of dry oxidn SO2 in 1/s
REAL :: wetrate                !rate of wet oxidn SO2 in 1/s
REAL :: dmsrate                !rate of dry oxidn DMS in 1/s

REAL :: clearf(row_length,rows,wet_model_levels)     !clear air fraction
REAL :: qctotal(row_length,rows,wet_model_levels)    !tot condensed water
REAL :: rho_air(row_length,rows,model_levels)        !air density kg/m3
REAL :: air_conc(row_length,rows,model_levels)       !air concn  mcls/cc
REAL :: viscos_air(row_length,rows,model_levels)     !air viscosity (kg/m/s)
REAL :: mean_freep(row_length,rows,model_levels)     !mean free path of air molecules (m)
REAL :: rel_hum(row_length,rows,model_levels)        !relative humidity(0-1)
REAL :: p_array(row_length,rows,model_levels)        !pressure, no halo
REAL :: so2_array(row_length,rows,model_levels)      !SO2, no halo
REAL :: so4_acc_array(row_length,rows,model_levels)  !SO4_ACC,no halo
REAL :: so4_ait_array(row_length,rows,model_levels)  !SO4_AIT,no halo
REAL :: so4_dis_array(row_length,rows,model_levels)  !SO4_DIS,no halo
REAL :: qcl_array(row_length,rows,wet_model_levels)  !QCL, no halo
REAL :: qcf_array(row_length,rows,wet_model_levels)  !QCF, no halo
REAL :: q_array(row_length,rows,wet_model_levels)    !Q, no halo
REAL :: h2o2_mxr_array(row_length,rows,model_levels) !H2O2_MXR,no halo
REAL :: nh3_array(row_length,rows,model_levels)      !NH3, no halo
REAL :: dms_array(row_length,rows,model_levels)      !DMS, no halo

! Removal of oxidants (for 2-way coupling with UKCA):
REAL :: delta_oh_dms(row_length,rows,model_levels)
                                ! Removal of OH in one timestep through reaction with DMS
REAL :: delta_o3_dms(row_length,rows,model_levels)
                                ! Removal of O3 in one timestep through reaction with DMS
REAL :: delta_ho2_dms(row_length,rows,model_levels)
                                ! Removal of HO2 in one timestep through reaction with DMS
REAL :: delta_oh_so2(row_length,rows,model_levels)
                                ! Removal of OH in one timestep through reaction with SO2
REAL :: delta_h2o2_so2(row_length,rows,model_levels)
                                ! Removal of H2O2 in one timestep through reaction with SO2
REAL :: delta_o3_so2(row_length,rows,model_levels)
                                ! Removal of O3 in one timestep through reaction with SO2
REAL :: delta_oh_tot(row_length,rows,model_levels)
REAL :: delta_ho2_tot(row_length,rows,model_levels)
REAL :: delta_o3_tot(row_length,rows,model_levels)
REAL :: delta_h2o2_tot(row_length,rows,model_levels)
!
REAL :: evaptime         ! timescale for cloud droplets to evaporate
REAL :: nuctime          ! timescale for particles to enter a cloud and nucleate.
REAL :: diffuse_tau      ! diffusive lifetime of Aitken mode particles once they enter a cloud
REAL :: rho_cuberoot     ! cube root of density of water.
REAL :: diffuse_tauratio ! CLOUDTAU/DIFFUSE_TAU
REAL :: probdiff_inv     ! inverse of probability of a particle being  diffusionallly scavenged in a cloud.
REAL :: probdiff_fn1     ! PROBDIFF_INV - 0.5
REAL :: probdiff_fn2     ! PROBDIFF_INV*EXP(DIFFUSE_TAURATIO*0.5)
REAL :: probdiff_cloud   ! probability of an Aitken mode particle being in cloud at the start of a step.
REAL :: probdiff_clear   ! probability of an Aitken mode particle being in clear air at the start of a step.
REAL :: lambda_aitken    ! ratio of concentrations of Aitken mode particles in cloudy to clear air.
REAL :: diffrate         ! rate of diffusive capture of Aitken mode particles

REAL :: w_conc           ! H2O concentration in molecules/cc
REAL :: k_so2_oh         ! reaction rate for SO2+OH  cc/mcl/s
REAL :: k_so2oh_lo       ! low_pressure reaction rate limit
REAL :: k_limrat         ! ratio  reaction rate limits (LO/HI)
REAL :: k_ho2_ho2        ! reaction rate for HO2 to H2O2
REAL :: h2o2_con_to_mxr  ! factor = RMM_H2O2/AVOGADRO/AIR DENSITY
REAL :: o2_conc          ! O2 concentration in molecules/cc
REAL :: o3_conc          ! O3 concentration in molecules/cc
REAL :: k1_dms_oh        ! reaction rate for DMS+OH via path 1
REAL :: k2_dms_oh        ! reaction rate for DMS+OH via path 2
REAL :: k3_ch3so2        ! thermal decomposition rate for CH3SO2
REAL :: k6_ch3so3        ! thermal decomposition rate for CH3SO3
REAL :: f_dms_to_ch3so2  ! fraction of oxidised DMS becoming CH3SO2
REAL :: f_dms_to_ch3so3  ! fraction of oxidised DMS becoming CH3SO3
REAL :: f_ch3so2_to_so2  ! fraction of CH3SO2 becoming SO2
REAL :: f_ch3so3_to_so4  ! fraction of CH3SO3 becoming SO4

REAL :: fac2,fac3              !for interp between LO and HI K_SO2OH
!                                 !    reaction rate limits

! dummy variables to assist wet oxidn calculation
REAL :: term1  
REAL :: term2
REAL :: denom
! dummy variables to assist calculation of diffusive capture.
REAL :: term3  
REAL :: term4  
REAL :: tauratio    ! CLOUDTAU/CHEMTAU
REAL :: hftaurat    ! CLOUDTAU/2*CHEMTAU
REAL :: probox_inv  ! 1/probability of oxidn in cloud
REAL :: probox_fn1  ! PROBOX_INV - 0.5
REAL :: probox_fn2  ! PROBOX_INV * EXP(-HFTAURAT)
REAL :: prob_cloud  ! probability SO2 starts in cloud
REAL :: prob_clear  ! probability SO2 starts in clear air
REAL :: lamda       ! ratio SO2 in cloud/clear air
REAL :: h2o2_rep    ! replenished H2O2


!     Extra variables for improved diffusional scavenging.

REAL ::  diffusivity      !Mean diffusion coefficent of Aitken particles
REAL ::  viscosity_air    !Dynamic viscosity of air (kg m-1 s-1).
REAL ::  mean_free_path   !Mean free path of air molecules.
REAL ::  knudsen_weber    !Expression related to the Cunningham slip
!                          flow correction to the diffusion coefficient
REAL ::  log_sigma_ait    !Natural log of SIGMA_AIT (in C_SULCHM).
REAL ::  sq_log_sigma_ait !Square of the previous parameter.
REAL ::  diff_con1        !Term in diffusion coefficent formula.
REAL ::  diff_con2        !Term in diffusion coefficent formula.
REAL ::  pec              !Quantity associated with (not =) Peclet number
REAL ::  work_radius      !Variable linked to average droplet radius.
REAL ::  scavcd           !Scavenging coefficient for in-cloud
                          !advective-diffusive removal of Aitken mode particles

! Extra variables for coagulation code

REAL :: coagrate         !Rate of transfer of SO4_AIT to SO4_ACC
!                         due to coagulation (kg kg-1 s-1)
REAL :: alpha            !Hygroscopic growth parameter in
                         !Fitzgerald's scheme (assume BETA=1)
REAL :: sum_iy           !Sum of COAG_QI_ functions in coag calc
REAL :: sum_iz           !              "
REAL :: yfac             !Coefficient used in coag calcn
REAL :: zfac             !         "

! COAG_QI_ functions
REAL :: coag_qi_1        !COAG_QI for A=1,    B=0
REAL :: coag_qi_2        !   "        A=4/3,  B=-1/3
REAL :: coag_qi_3        !   "        A=2/3,  B=1/3
REAL :: coag_qi_4        !   "        A=2/3,  B=0
REAL :: coag_qi_5        !   "        A=4/3,  B=-2/3
REAL :: coag_qi_6        !   "        A=1/3,  B=1/3
REAL :: coag_qi_7        !   "        A=1,    B=-1/3

REAL,PARAMETER :: acunn = 1.591 ! Cunningham correction factor
REAL :: log_sigma_acc           !nat log of geom std dev of acc mode distn
REAL :: sq_log_sigma_acc        !square of   "

REAL :: mmr_star        ! threshold mmr of SO4_ACC below which PSI=1
! constants required for PSI calcn
REAL :: con_1, con_2
REAL :: scale_factor            ! to prevent negative SO2
REAL :: scale_fact_ait          ! to prevent negative SO4_AIT
REAL :: scale_factor_oh         ! to prevent negative OH_conc
REAL :: scale_factor_o3         ! to prevent negative O3

! Arrays and variables required for mode-merging calculation:
REAL :: ait_production(row_length, rows, model_levels)
REAL :: frac_trans(row_length, rows, model_levels)

! FRAC_TRANS = K_MERGE * AIT_PRODUCTION
REAL,PARAMETER :: k_merge = 3.068e+09 
REAL :: deltas_merge(row_length, rows, model_levels)

! variables for VECTLIB code 
REAL :: tmp(row_length)
REAL :: tmp1(row_length)
REAL :: fca(row_length)
REAL :: parha(row_length)
REAL :: jlog10e
REAL :: jfac1

! External function
REAL, EXTERNAL  :: coag_qi           !Function for coagulation calcn

IF (lhook) CALL dr_hook('SULPHR',zhook_in,zhook_handle)

!--------------------------------------------------------------------
! 0. Set up constants and initialise OUT arrays to 0
!--------------------------------------------------------------------

! The next constant cannot be declared as a PARAMETER because
! it involves non-integer exponents.

rho_cuberoot = rho_water**0.333333

! Parameters used: g, r, acunn, rho_water, sigma_ait, rad_ait, cloudtau,
! chemtau, rmms_air, evaptau, thold
!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(k, j, i, log_sigma_ait,                &
!$OMP& diff_con2, log_sigma_acc, sq_log_sigma_acc, coag_qi_1, diff_con1,       &
!$OMP& coag_qi_2, coag_qi_3, coag_qi_4, coag_qi_5, coag_qi_6, coag_qi_7,       &
!$OMP& sum_iz, sum_iy, tauratio, hftaurat, probox_inv, probox_fn1,             &
!$OMP& probox_fn2, con_1, con_2, prob_clear, prob_cloud, mean_free_path,       &
!$OMP& jlog10e, jfac1, parha, fca, tmp, k_limrat, k_so2oh_lo,                  &
!$OMP& dryrate, fac2, fac3, term1, term2, lamda, denom, wetrate,               &
!$OMP& k1_dms_oh, o2_conc, k2_dms_oh, k3_ch3so2, f_ch3so2_to_so2,              &
!$OMP& k6_ch3so3, dmsrate, f_dms_to_ch3so2, f_ch3so3_to_so4,                   &
!$OMP& evaptime, nuctime, diffuse_tau, diffusivity, pec, knudsen_weber,        &
!$OMP& work_radius, scavcd, diffuse_tauratio, probdiff_inv,                    &
!$OMP& probdiff_fn1, probdiff_fn2, diffrate, zfac, alpha, coagrate,            &
!$OMP& drh, scale_factor, w_conc,viscosity_air, h2o2_con_to_mxr,               &
!$OMP& h2o2_rep, k_ho2_ho2, probdiff_clear, scale_factor_oh, term4,            &
!$OMP& scale_factor_o3,lambda_aitken, probdiff_cloud, k_so2_oh,                &
!$OMP& tmp1, o3_conc, scale_fact_ait, yfac, term3, sq_log_sigma_ait)

!
IF (l_sulpc_online_oxidants) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        IF (oh_conc(i,j,k) .LT. 0.0) oh_conc(i,j,k) = 0.0
        IF (ho2_conc(i,j,k) .LT. 0.0) ho2_conc(i,j,k) = 0.0
        IF (h2o2_mxr(i,j,k) .LT. 0.0) h2o2_mxr(i,j,k) = 0.0
        IF (o3(i,j,k) .LT. 0.0) o3(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!     Extra parameters for improved diffusional scavenging.

log_sigma_ait = LOG(sigma_ait)
sq_log_sigma_ait = log_sigma_ait * log_sigma_ait
diff_con1 = EXP(-2.5*sq_log_sigma_ait)/rad_ait
diff_con2 = EXP(-4.0*sq_log_sigma_ait)/(rad_ait*rad_ait)


!     Extra parameters for coagulation

log_sigma_acc=LOG(sigma)
sq_log_sigma_acc=log_sigma_acc*log_sigma_acc

! Calculate COAG_QI_ functions (not dependent on ALPHA)

! DEPENDS ON: coag_qi
coag_qi_1=coag_qi(1.0, 0.0,                                                    &
  rad_ait,rad_acc,sq_log_sigma_ait,sq_log_sigma_acc)

! DEPENDS ON: coag_qi
coag_qi_2=coag_qi(1.333333, -0.333333,                                         &
  rad_ait,rad_acc,sq_log_sigma_ait,sq_log_sigma_acc)

! DEPENDS ON: coag_qi
coag_qi_3=coag_qi(0.666667, 0.333333,                                          &
  rad_ait,rad_acc,sq_log_sigma_ait,sq_log_sigma_acc)

! DEPENDS ON: coag_qi
coag_qi_4=coag_qi(0.666667, 0.0,                                               &
  rad_ait,rad_acc,sq_log_sigma_ait,sq_log_sigma_acc)

! DEPENDS ON: coag_qi
coag_qi_5=coag_qi(1.333333, -0.666667,                                         &
  rad_ait,rad_acc,sq_log_sigma_ait,sq_log_sigma_acc)

! DEPENDS ON: coag_qi
coag_qi_6=coag_qi(0.333333, 0.333333,                                          &
  rad_ait,rad_acc,sq_log_sigma_ait,sq_log_sigma_acc)

! DEPENDS ON: coag_qi
coag_qi_7=coag_qi(1.0, -0.333333,                                              &
  rad_ait,rad_acc,sq_log_sigma_ait,sq_log_sigma_acc)

! Calculate sums of COAG_QI_ functions

sum_iz = 2.0*coag_qi_1 + coag_qi_2 + coag_qi_3

sum_iy = coag_qi_4 + coag_qi_5 + coag_qi_6 + coag_qi_7

! Calculate parameters which are same for each grid box (for wet oxidn)
! and set up other arrays

tauratio=cloudtau/chemtau
hftaurat=tauratio/2.0
probox_inv=1.0/(1.0-EXP(-tauratio))
probox_fn1=probox_inv-0.5
probox_fn2=probox_inv*EXP(-hftaurat)




!cdir collapse

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      p_array(i,j,k)=p(i,j,k)
      so2_array(i,j,k)=so2(i,j,k)
      so4_acc_array(i,j,k)=so4_acc(i,j,k)
      so4_ait_array(i,j,k)=so4_ait(i,j,k)
      msa(i,j,k)  = 0.0
      f_dms_to_so2(i,j,k)  = 0.0
      f_dms_to_so4(i,j,k)  = 0.0
      f_dms_to_msa(i,j,k)  = 0.0
      nh3_dep(i,j,k)  = 0.0
      deltas_dry(i,j,k)  = 0.0
      deltas_wet(i,j,k)  = 0.0
      deltas_wet_o3(i,j,k) = 0.0
      deltas_tot(i,j,k)  = 0.0
      deltas_dms(i,j,k)  = 0.0
      deltas_evap(i,j,k) = 0.0
      deltas_nucl(i,j,k) = 0.0
      deltas_diffuse(i,j,k) = 0.0
      deltas_coag(i,j,k) = 0.0
      deltas_tot_ait(i,j,k) = 0.0
      IF (l_sulpc_2_way_coupling) THEN
        delta_oh_dms(i,j,k) = 0.0
        delta_o3_dms(i,j,k) = 0.0
        delta_ho2_dms(i,j,k) = 0.0
        delta_oh_so2(i,j,k) = 0.0
        delta_h2o2_so2(i,j,k) = 0.0
        delta_o3_so2(i,j,k) = 0.0
        delta_oh_tot(i,j,k) = 0.0
        delta_o3_tot(i,j,k) = 0.0
      END IF

      rho_air(i,j,k)=p_array(i,j,k)/(r*t(i,j,k))
      air_conc(i,j,k)=                                                         &
        (avogadro*p_array(i,j,k))/(r*t(i,j,k)*rmm_air*1.0e6)

      ! Compute dynamic viscosity of air and mean free path of air molecules
      ! using the formulae on P.417 of Pruppacher and Klett, 2nd edition.

      viscos_air(i,j,k)=(1.718 + (t(i,j,k)-273.15)*0.0049 )*1.0e-5
      mean_freep(i,j,k)=                                                       &
        (mfp_ref*pref_mfp*t(i,j,k))/(tref_mfp*p_array(i,j,k))

    END DO
  END DO
END DO
!$OMP END DO


!---------------------------------------------------------------------
! 1.  This section calculates amount of SO2 dry oxidised using a
!     3_D OH concentration field to control the oxidn rate when
!     cos(zenith angle) is G.T. 10-6  (else rate is zero).
!  Reaction rates are dependent on temperature and pressure
!---------------------------------------------------------------------

con_1 = EXP(4.5*LOG(sigma)*LOG(sigma))
con_2 = 4.0*pi*rho_so4*chi*num_star*(rad_acc)**3

jlog10e=1.0/LOG(10.0)
jfac1=jlog10e/fac1
DO i=1,row_length
  fca(i)=fc
  parha(i)=parh
END DO


!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
!cdir collapse
  DO j=1,rows
    DO i=1,row_length
      tmp(i)=t1/t(i,j,k)
    END DO
    CALL rtor_v(row_length,tmp,parha,tmp)
    DO i=1,row_length
      k_so2oh_lo=air_conc(i,j,k)*k1*tmp(i)
      k_limrat = k_so2oh_lo / k_so2oh_hi
      tmp(i) = k_limrat
      fac2 = k_so2oh_lo /(1.0+ k_limrat)
      tmp1(i) = fac2
    END DO
    CALL log_v(row_length,tmp,tmp)
    DO i=1,row_length
      fac3 = 1.0 + (tmp(i)*jfac1)**2
      tmp(i)=1.0/fac3
      fca(i)=fc
    END DO
    CALL rtor_v(row_length,fca,tmp,fca)
    DO i=1,row_length
      k_so2_oh = tmp1(i) * fca(i)

      !  Only calculate the dry oxidation rate in the daytime.

      IF (cosza2d(i,j) >  1.0e-06) THEN
        dryrate=k_so2_oh * oh_conc(i,j,k)
      ELSE                     ! Zero rate if nightime
        dryrate=0.0
      END IF

      deltas_dry(i,j,k)=dryrate*tstep*so2_array(i,j,k) !AmntSO2dryox

      ! Calculate change in OH_CONC due to this reaction:
      IF (l_sulpc_online_oxidants .AND.                                        &
        l_sulpc_2_way_coupling) THEN
        delta_oh_so2(i,j,k) = deltas_dry(i,j,k) * air_conc(i,j,k)              &
          * rmm_air / rmm_s
      END IF

    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!--------------------------------------------------------------
! 2. This section calculates amount of SO2 wet oxidised using a
!    3_D H2O2 field to control the oxidn rate. As H2O2 is depleted, it
!    is replenished from the HO2 concn field (via a reaction rate), but
!    it may not exceed the H2O2 LIMIT field .
!    The reaction rate is a function of pressure, temp and humidity.
!    It uses D.Roberts' formula for the wet oxidn rate, and R.Colvile's
!    H2O2 and HO2 chemistry.
!    Depletion and replenishment of H2O2 is done at the end of the
!    routine, after scaling to avoid SO2 becoming negative.

!    If H2O2 is limiting, further oxidation by O3 occurs if there is
!    sufficient NH3 available to buffer the reaction; Choularton's
!    suggested procedure is:
!    a) Do oxidn of SO2 by H2O2
!    b) Neutralise sulphate produced with NH3
!    c) If SO2 remains do oxidn by O3 until SO2 or NH3 required to
!       neutralise it is exhausted.

!--------------------------------------------------------------


!$OMP DO SCHEDULE(STATIC)
DO k=1,wet_model_levels
!cdir collapse
  DO j=1,rows
    DO i=1,row_length
      qcl_array(i,j,k)=qcl(i,j,k)
      qcf_array(i,j,k)=qcf(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

IF (l_sulpc_online_oxidants) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
!cdir collapse
    DO j=1,rows
      DO i=1,row_length
        h2o2_mxr_array(i,j,k)=h2o2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

ELSE

!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
!cdir collapse
    DO j=1,rows
      DO i=1,row_length
        h2o2_mxr_array(i,j,k)=h2o2_mxr(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

END IF

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
!cdir collapse
  DO j=1,rows
    DO i=1,row_length
      nh3_array(i,j,k)=nh3(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k=1,wet_model_levels
!cdir collapse
  DO j=1,rows
    DO i=1,row_length

      ! Calculate clear air fraction
      clearf(i,j,k)=1.0-cloudf(i,j,k)

      IF ((qcl_array(i,j,k) >= thold) .AND.                                    &
        &        (cloudf(i,j,k) >  0.0)) THEN

        !   CALCULATE LAMDA
        term1=(clearf(i,j,k)*tauratio)**2
        term1=term1+2.0*tauratio*clearf(i,j,k)*                                &
          &                             (clearf(i,j,k)-cloudf(i,j,k))
        term1=SQRT(1.0+term1)
        term2=2.0*cloudf(i,j,k)-tauratio*clearf(i,j,k)-1.0
        term2=term2+term1
        lamda=term2/(2.0*cloudf(i,j,k))

        !   CALCULATE PROB_CLEAR AND PROB_CLOUD
        denom=clearf(i,j,k)+cloudf(i,j,k)*lamda
        prob_clear=clearf(i,j,k)/denom
        prob_cloud=cloudf(i,j,k)*lamda/denom

        !   CALCULATE EXPECTED LIFETIME OF SO2 (WHICH IS 1/WETRATE)
        term1=probox_fn1*prob_clear
        term2=probox_fn2*prob_cloud
        term2=(term1+term2)*clearf(i,j,k)/cloudf(i,j,k)
        denom=term2*cloudtau+chemtau
        wetrate=1.0/denom

        !!  RC'S H2O2 CHEMISTRY: OXIDATION OF SO2 TO SO4 IS CONTROLLED BY THE
        !                        SMALLER OF THE SO2 AND H2O2 MIX RAT FIELDS

        IF (so2_array(i,j,k) <= (h2o2_mxr_array(i,j,k)*                        &
          relm_s_h2o2)) THEN

          deltas_wet(i,j,k)=(1.0-EXP(-wetrate*tstep))*                         &
            so2_array(i,j,k)

        ELSE

          deltas_wet(i,j,k)=(1.0-EXP(-wetrate*tstep))*                         &
            h2o2_mxr_array(i,j,k)*relm_s_h2o2

          !   H2O2 field is limiting, so oxidise SO2 further with O3, using NH3
          !   as buffer, if sufficient O3 and NH3 available.

          IF (l_sulpc_ozone) THEN
            IF (l_sulpc_so2_o3_nonbuffered) THEN

              ! SO2+O3 NOT buffered by NH3.
              ! When L_SULPC_SO2_O3_NONBUFFERED = T, O3 oxidation proceeds only if
              ! O3 MMR is above a threshold (O3_MIN). The condition on NH3, which
              ! applies when L_SULPC_SO2_O3_NONBUFFERED = F, does not apply when
              ! L_SULPC_SO2_O3_NONBUFFERED = T.

              IF (o3(i,j,k) >= o3_min) THEN

                ! When L_SULPC_SO2_O3_NONBUFFERED = T, O3 oxidation is controlled by SO2
                ! field only.  Note that all SO2 is used in DELTAS_WET_O3 calcn;
                ! adjustment by SCALE_FACTOR at end of routine prevents removing too much.

                deltas_wet_o3(i,j,k)=                                          &
                  (1.0-EXP(-wetrate*tstep))*so2_array(i,j,k)
              END IF               ! End ozone threshold condtion.
            ELSE

              ! SO2+O3 buffered by NH3.

              ! When L_SULPC_SO2_O3_NONBUFFERED = F, O3 oxidation proceeds only if O3
              ! MMR is above a threshold (O3_MIN) and NH3 MMR is large enough that there
              ! is still NH3 left in the gridbox after production of ammonium sulphate
              ! in the aqueous reaction of SO2 with H2O2 (see above).

              IF ((o3(i,j,k) >= o3_min) .AND.                                  &
                (nh3_array(i,j,k) > deltas_wet(i,j,k)/relm_s_2n)) THEN

                ! When L_SULPC_SO2_O3_NONBUFFERED = F, O3 oxidation is controlled by s
                ! smaller of SO2 and NH3 fields. 2 atoms of N required for each S atom
                ! in ammonium sulphate.  Sufficient NH3 must be available to neutralise
                ! all sulphate produced in grid box for O3 reaction to continue (otherwise
                ! PH too low) Note that all SO2 or NH3 is used in DELTAS_WET_O3 calcn;
                ! adjustment by SCALE_FACTOR at end of routine prevents removing too much.

                IF (so2_array(i,j,k)<=(nh3_array(i,j,k)                        &
                  *relm_s_2n)) THEN
                  ! SO2 limiting
                  deltas_wet_o3(i,j,k)=                                        &
                    (1.0-EXP(-wetrate*tstep))*so2_array(i,j,k)

                ELSE
                  ! NH3 limiting
                  deltas_wet_o3(i,j,k)=(1.0-EXP(-wetrate*tstep))*              &
                    nh3_array(i,j,k)*relm_s_2n

                END IF

              END IF    ! End ozone and NH3 thresholds condition.
              !
            END IF   ! End of IF(L_SULPC_SO2_O3_NONBUFFERED)-THEN-ELSE.
            !
            ! Calculate change in O3 MMR due to aqueous reaction with SO2:
            IF (l_sulpc_online_oxidants) THEN
              delta_o3_so2(i,j,k) = deltas_wet_o3(i,j,k)                       &
                * rmm_o3 / rmm_s
            END IF
          END IF     ! End L_SULPC_OZONE test

        END IF           ! End peroxide oxidation block
        !
        ! Calculate change in H2O2 MMR due to aqueous reaction with SO2:
        IF (l_sulpc_online_oxidants) THEN
          delta_h2o2_so2(i,j,k) = deltas_wet(i,j,k)                            &
            * rmm_h2o2 / rmm_s
        END IF
      END IF           ! End cloud present condition


    END DO
  END DO
END DO
!$OMP END DO NOWAIT


!------------------------------------------------------------------
! 3.  This section calculates amount of DMS dry oxidised to SO2
!     and MSA using a 3D OH concentration field to control the oxidn
!     rate when cos zenith angle is G.T. 10-6  (else rate is zero).
!     MSA is not saved, and K_DMS_OH is constant in this version.
!    OR, if O3 available and new DMS chemistry required -
!     This section calculates the amount of DMS dry oxidised to SO2,
!     SO4 and MSA using OH, HO2, O3 and O2 as oxidants. OH and HO2
!     are available in daylight only. Intermediate species have
!     short lifetimes so are not stored between timesteps.
!------------------------------------------------------------------

IF (l_sulpc_dms) THEN


!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
!cdir collapse
    DO j=1,rows
      DO i=1,row_length
        dms_array(i,j,k)=dms(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  IF (l_sulpc_newdms .AND. l_sulpc_ozone) THEN


!$OMP DO SCHEDULE(STATIC)
    DO k=1,model_levels
!cdir collapse
      DO j=1,rows
        DO i=1,row_length

          ! Calculate rate coefficients K1_DMS_OH and K2_DMS_OH for the two
          !    DMS+OH -> CH3SO2 reactions, which are daylight dependent

          IF (cosza2d(i,j) >  1.0e-06) THEN     ! daylight

            k1_dms_oh=1.13e-11*EXP(-254.0/t(i,j,k))

            !  Calculate concentration of O2 as 20.95% of air concentration

            o2_conc=0.2095*air_conc(i,j,k)

            k2_dms_oh=                                                         &
              ( 1.7e-42*o2_conc*EXP(7810.0/t(i,j,k)) ) /                       &
              ( 1.0 + ( 5.5e-31*o2_conc*EXP(7460.0/t(i,j,k)) ) )

            !  Calculate thermal decomposition rate coefficient K3_CH3SO2 for
            !         CH3SO2 -> CH3+SO2

            k3_ch3so2=2.6e11*EXP(-9056.0/t(i,j,k))

            !  Rate coefficient K4_CH3SO2_O3 for CH3SO2+O3 -> CH3SO3+O2 is const
            !         K4_CH3SO2_O3=1.0E-14  in C_SULCHM

            !  Rate coefficient K5_CH3SO3_HO2 for CH3SO3+HO2 -> MSA+O2 is const
            !         K5_CH3SO3_HO2=4.0E-11  in C_SULCHM

            !  Calculate thermal decomposition rate coefficient K6_CH3SO3 for
            !         CH3SO3 -> CH3+SO3

            k6_ch3so3=1.1e17*EXP(-12057.0/t(i,j,k))

            !  Calculate DMS oxidation rate

            dmsrate = (k1_dms_oh + k2_dms_oh) * oh_conc(i,j,k)

            !  Calculate fraction of DELTAS_DMS becoming CH3SO2
            f_dms_to_ch3so2=(k1_dms_oh+0.9*k2_dms_oh)/                         &
              &                    (k1_dms_oh+k2_dms_oh)

            !  Calculate fraction of DELTAS_DMS becoming SO2 (requires O3_CONC)
            o3_conc= o3(i,j,k) * air_conc(i,j,k) * rmm_air/rmm_o3
            f_ch3so2_to_so2=                                                   &
              & k3_ch3so2/(k3_ch3so2+k4_ch3so2_o3*o3_conc)
            f_dms_to_so2(i,j,k)=                                               &
              & f_dms_to_ch3so2 * f_ch3so2_to_so2

            !  Calculate fraction of DELTAS_DMS becoming sulphate
            f_ch3so3_to_so4= k6_ch3so3 /                                       &
              &  (k6_ch3so3+k5_ch3so3_ho2*ho2_conc(i,j,k))
            f_dms_to_so4(i,j,k)=                                               &
              & f_dms_to_ch3so2 * (1.0-f_ch3so2_to_so2)                        &
              &                    * f_ch3so3_to_so4

            !  Calculate fraction of DELTAS_DMS becoming MSA
            f_dms_to_msa(i,j,k)=1.0 - f_dms_to_so2(i,j,k)                      &
              &                  - f_dms_to_so4(i,j,k)

            !  Calculate amount of DMS oxidised in TSTEP (in daylight only)

            !DELTA_DMS
            deltas_dms(i,j,k)=dmsrate*tstep*dms_array(i,j,k)
            !
            ! Calculate changes in oxidants due to reaction with DMS:
            IF (l_sulpc_online_oxidants .AND.                                  &
              & l_sulpc_2_way_coupling) THEN
              delta_oh_dms(i,j,k) =                                            &
                &  deltas_dms(i,j,k) * air_conc(i,j,k)                         &
                &  * rmm_air / rmm_s
              delta_o3_dms(i,j,k) = deltas_dms(i,j,k)                          &
                &  * (1.0 - f_ch3so2_to_so2)  * rmm_o3 /                       &
                & rmm_s
              delta_ho2_dms(i,j,k) = deltas_dms(i,j,k)                         &
                &  * (1.0 - f_ch3so2_to_so2) *                                 &
                & (1.0 - f_ch3so3_to_so4)                                      &
                & * air_conc(i,j,k) * rmm_air / rmm_s
            END IF
          END IF             ! End daylight condition

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT


  ELSE                 ! implement simple DMS scheme

!$OMP DO SCHEDULE(STATIC)
    DO k=1,model_levels
!cdir collapse
      DO j=1,rows
        DO i=1,row_length

          !  Calculate DMS oxidation rate if daytime.
          IF (cosza2d(i,j) >  1.0e-06) THEN
            dmsrate = k_dms_oh * oh_conc(i,j,k)
          ELSE                 ! ZERO RATE IF NIGHTIME
            dmsrate=0.0
          END IF                   ! END COSZA2D IF

          !  CALCULATE FRACTION OF DMS OXIDISED

          !AmtDMSdryox
          deltas_dms(i,j,k)=dmsrate*tstep*dms_array(i,j,k)

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT


  END IF                 ! End L_SULPC_NEWDMS Test

END IF                   ! END L_SULPC_DMS IF

!---------------------------------------------------------------------
! 4. Release of aerosol from evaporating cloud droplets:
!    if no condensed water (liquid + ice) in grid box, release
!    dissolved sulphate as accumulation mode aerosol.
!     If cloud fraction less than 0.95, release some in clear  air.
!--------------------------------------------------------------------


!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      so4_dis_array(i,j,k)=so4_dis(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT




!$OMP DO SCHEDULE(STATIC)
DO k=1,wet_model_levels
!cdir collapse
  DO j=1,rows
    DO i=1,row_length
      qctotal(i,j,k)= qcl_array(i,j,k) + qcf_array(i,j,k)
      IF (qctotal(i,j,k)  <   thold) THEN
        deltas_evap(i,j,k)=so4_dis_array(i,j,k)
      ELSE IF  (cloudf(i,j,k) <  0.95) THEN
        evaptime = evaptau + 0.5*cloudtau
        deltas_evap(i,j,k)=                                                    &
          &   (1.0-EXP(-tstep/evaptime))*so4_dis_array(i,j,k)
      ELSE
        deltas_evap(i,j,k)=0.0
      END IF

    END DO
  END DO
END DO
!$OMP END DO


!     Also release any dissolved aerosol in a non-wet level,
!     where it should not be.

IF (wet_model_levels  <   model_levels)  THEN


!$OMP DO SCHEDULE(STATIC)
  DO k=wet_model_levels+1,model_levels
!cdir collapse
    DO j=1,rows
      DO i=1,row_length
        deltas_evap(i,j,k)= so4_dis_array(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO


END IF

!-------------------------------------------------------------------
! 5. Nucleation of accumulation mode aerosol forming dissolved SO4
!-------------------------------------------------------------------

!    THIS CODE ASSUMES THAT THE PARAMETER NUCTAU, WHICH IS THE
!    TIMESCALE FOR NUCLEATION ONCE A PARTICLE ENTERS A CLOUD, IS
!    VERY SHORT COMPARED WITH CLOUDTAU.

!cdir collapse
!$OMP DO SCHEDULE(STATIC)
DO k=1,wet_model_levels
  DO j=1,rows
    DO i=1,row_length

      IF ((qctotal(i,j,k) >= thold) .AND.                                      &
        &       (cloudf(i,j,k) >  0.0))  THEN
        nuctime=nuctau+((clearf(i,j,k)*cloudtau)/                              &
          & (2.0*cloudf(i,j,k)))
        deltas_nucl(i,j,k)=                                                    &
          &       (1.0-EXP(-tstep/nuctime))*so4_acc_array(i,j,k)
      END IF

    END DO
  END DO
END DO
!$OMP END DO

!-------------------------------------------------------------------
! 6. Diffusional scavenging of Aitken mode SO4 forming dissolved SO4
! Improved code (introduced at vn5.3) allows for:
! (1) The enhancement of diffusional scavenging due to the
! relative motion of cloud droplets and aerosol is modelled.
! (2) The variation of particle diffusivity with ambient
! temperature and pressure is taken into account.
! (3) Averaging over the size distributions of both the cloud
! droplets and the aerosol particles has been done, assuming
! a Khrgian-Mazin distribution for the droplets and a lognormal
! distribution for the aerosol.

!-------------------------------------------------------------------

!    THIS IS A MUCH SLOWER PROCESS THAN NUCLEATION AND THEREFORE WE
!    CANNOT MAKE THE SAME APPROXIMATIONS USED IN THAT CASE.

!cdir collapse

!$OMP DO SCHEDULE(DYNAMIC)
DO k=1,wet_model_levels
  DO j=1,rows
    DO i=1,row_length

      IF ((qctotal(i,j,k) >= thold) .AND.                                      &
        (cloudf(i,j,k) >  0.0)) THEN

        !    FIRST COMPUTE IN-CLOUD TIMESCALE FOR DIFFUSIONAL CAPTURE,
        !    USING TOTAL CONDENSED WATER WITHIN THE CLOUDY PORTION OF THE BOX.
        !    THE DIFFERENCE BETWEEN LIQUID WATER AND ICE IS NEGLECTED HERE.
        !    THIS SHOULD BE IMPROVED ON EVENTUALLY.


        ! Compute mean free path of air molecules.
        !  (See P.417 of Pruppacher and Klett, 2nd edition.)

        mean_free_path=                                                        &
          (mfp_ref*pref_mfp*t(i,j,k))/(tref_mfp*p_array(i,j,k))

        ! Compute the Knudsen and Weber term for a particle of the median
        !  radius (note approximation here: we do not average over the size
        !  distribution). See P.450 of Pruppacher and Klett, 2nd edition.

        knudsen_weber=1.257 + 0.4*                                             &
          EXP(-((1.1*rad_ait)/mean_free_path))

        ! Temporarily use DIFFUSIVITY to store working value.

        diffusivity=diff_con1 + diff_con2*mean_free_path*                      &
          knudsen_weber

        ! Compute dynamic viscosity of air, using an approximate version of
        !  the formula on P.417 of Pruppacher and Klett, 2nd edition.

        viscosity_air=(1.718 + (t(i,j,k)-273.15)*0.0049 )*1.0e-5

        ! Now compute mean diffusion coefficient.

        diffusivity=(boltzmann*t(i,j,k)*diffusivity)/                          &
          (6.0*pi*viscosity_air)

        ! Now compute the term PEC related to (but not equal to) the cube
        !  root of the Peclet Number.

        pec=((4.0*g*rho_water)/                                                &
          (9.0*diffusivity*viscosity_air))**0.333333

        work_radius=qctotal(i,j,k)/                                            &
          (cloudf(i,j,k)*10.0*pi*rho_water*n_droplet(i,j,k))

        work_radius=work_radius**0.333333

        ! We can finally compute the timescale for diffusive
        !  scavenging once inside cloud, DIFFUSE_TAU.

        scavcd=6.0*pi*diffusivity*n_droplet(i,j,k)*work_radius*                &
          (1.0 + pec*work_radius)

        diffuse_tau=1.0/scavcd


        diffuse_tauratio = cloudtau/diffuse_tau
        probdiff_inv = 1.0/( 1.0 - EXP(-diffuse_tauratio) )
        probdiff_fn1 = probdiff_inv - 0.5
        probdiff_fn2 = probdiff_inv*EXP(-(0.5*diffuse_tauratio) )

        !     CALCULATE LAMBDA_AITKEN.

        term3 = (clearf(i,j,k)*diffuse_tauratio)**2
        term3 = term3 + ( 2.0*diffuse_tauratio                                 &
          *clearf(i,j,k)*(clearf(i,j,k)-cloudf(i,j,k)) )
        term3 = SQRT(1.0+term3)
        term4=2.0*cloudf(i,j,k)-diffuse_tauratio*clearf(i,j,k)-1.0
        term4 = term4 + term3
        lambda_aitken = term4/(2.0*cloudf(i,j,k))

        !   CALCULATE PROBDIFF_CLEAR AND PROBDIFF_CLOUD

        denom = clearf(i,j,k)+cloudf(i,j,k)*lambda_aitken
        probdiff_clear = clearf(i,j,k)/denom
        probdiff_cloud = (cloudf(i,j,k)*lambda_aitken)/denom

        !   CALCULATE EXPECTED LIFETIME OF AN AITKEN MODE PARTICLE WITH
        !   RESPECT TO DIFFUSIVE CAPTURE BY CLOUD DROPLETS.

        term3 = probdiff_fn1*probdiff_clear
        term4 = probdiff_fn2*probdiff_cloud
        term4 = (term3+term4)*clearf(i,j,k)/cloudf(i,j,k)
        denom = term4*cloudtau + diffuse_tau
        diffrate = 1.0/denom

        !   NOW COMPUTE THE AMOUNT OF AITKEN MODE SO4 CONVERTED TO DISSOLVED
        !   SO4 BY DIFFUSIVE CAPTURE DURING THE STEP.

        deltas_diffuse(i,j,k)=(1.0-EXP(-diffrate*tstep))*                      &
          so4_ait_array(i,j,k)
      END IF
    END DO
  END DO
END DO
!$OMP END DO NOWAIT


!---------------------------------------------------------------------
! 6.5. Coagulation of Aitken mode with accumulation mode SO4 to form
!      more of the accumulation mode. Allowance is made for hygroscopic
!      growth of particles for values of RH up to 0.97 using
!      Fitzgerald's scheme. For RH > 0.97, coagulation is neglected.
!---------------------------------------------------------------------

! First calculate relative humidity to determine hygroscopic
! growth factor ALPHA for particles; RH is obtained from q and
! qsat using RH=q(1-qsat)/(qsat(1-q))

! Use REL_HUM array to store returned QSAT values in call to QSAT_WAT
! to save space; then calculate fractional REL_HUM


!$OMP DO SCHEDULE(STATIC)
DO k=1,wet_model_levels
  DO j=1,rows
    DO i=1,row_length
      q_array(i,j,k)=q(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO 

!$OMP DO SCHEDULE(DYNAMIC)
DO k=1,wet_model_levels

  ! DEPENDS ON: qsat_wat
  CALL qsat_wat(rel_hum(1,1,k)                                                 &
    ,  t(1,1,k), p_array(1,1,k), theta_field_size)
!cdir collapse
  DO j=1,rows
    DO i=1,row_length

      rel_hum(i,j,k) = q_array(i,j,k)*(1.0-rel_hum(i,j,k)) /                   &
        ( rel_hum(i,j,k)*(1.0-q_array(i,j,k)) )

      IF (rel_hum(i,j,k)  <   0.3) THEN

        ! Assume no hygroscopic growth occurs, so ALPHA=1
        ! Calculate factors multiplying SUM_IZ, SUM_IY (simplified for ALPHA=1)

        zfac = 2.0 * boltzmann * t(i,j,k) *                                    &
          so4_ait_array(i,j,k)*so4_acc_array(i,j,k)
        zfac = zfac * (rho_air(i,j,k)/rho_so4)
        zfac = zfac/( 3.0*viscos_air(i,j,k) )

        yfac = (1.333333*pi)**0.333333
        yfac = yfac*acunn*mean_freep(i,j,k)

        coagrate = zfac*(sum_iz + yfac*sum_iy)

      ELSE IF (rel_hum(i,j,k)  <=  0.97) THEN

        ! Calculate hygroscopic growth factor ALPHA

        ! DEPENDS ON: hygro_fact
        CALL hygro_fact(rel_hum(i,j,k), alpha)

        ! Calculate factors multiplying SUM_IZ, SUM_IY (for ALPHA /= 1)

        zfac = 2.0 * boltzmann * t(i,j,k) *                                    &
          so4_ait_array(i,j,k)*so4_acc_array(i,j,k)
        zfac = zfac * (rho_air(i,j,k)/rho_so4)
        zfac = zfac/( 3.0*viscos_air(i,j,k) )
        zfac = zfac/( alpha**6 )

        yfac = (1.333333*pi)**0.333333
        yfac = yfac*acunn*mean_freep(i,j,k)
        yfac = yfac/alpha

        coagrate = zfac*(sum_iz + yfac*sum_iy)

      ELSE

        ! Neglect coagulation for RH > 0.97

        coagrate = 0.0

      END IF         !End RH condition

      ! Calculate amount of SO4 tranferred from Ait to acc mode so4:

      deltas_coag(i,j,k) = coagrate*tstep

    END DO           !End i loop
  END DO             !End j loop

END DO               !End k loop
!$OMP END DO


!---------------------------------------------------------------------
! 6a.          Calculation of PSI, fraction of dry-oxidied SO2 going
!              into Aitken mode (as opposed to accumulation mode).  The
!              method of calculating PSI has been changed, and the
!              calculation has been moved here from Section 1, as it
!              now depends on relative humidity, which is calculated in
!              Section 6.
!----------------------------------------------------------------------
drh = rel_hum_table(2) - rel_hum_table(1)

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      nearest_index(i,j,k) = (NINT((rel_hum(i,j,k)                             &
        - rel_hum_table(1)) / drh)) + 1
      IF (nearest_index(i,j,k)  >   101) THEN
        nearest_index(i,j,k) = 101
      END IF
      IF (nearest_index(i,j,k)  <   1) THEN
        nearest_index(i,j,k) = 1
      END IF
      a_array(i,j,k) = a_table(nearest_index(i,j,k))                           &
        * so4_ait_array(i,j,k)
      b_array(i,j,k) = b_table(nearest_index(i,j,k))                           &
        * so4_acc_array(i,j,k)
      IF (a_array(i,j,k)  ==  0.0) THEN
        psi(i,j,k) = 1.0
      ELSE
        psi(i,j,k) = a_array(i,j,k) /                                          &
          (a_array(i,j,k) + b_array(i,j,k))
      END IF
    END DO        !End i loop
  END DO           !End j loop
END DO              !End k loop
!$OMP END DO NOWAIT


!-------------------------------------------------------------------
! 7. UPDATE SO2, SO4 modes and DMS, assuming:
!    Dry oxidation SO2 produces SO4 Aitken mode aerosol
!    Wet oxidation SO2 produces dissolved SO4
!    Dry oxidation DMS produces SO2 and MSA,
!                       and SO4 if L_SULPC_NEWDMS=T
!---------------------------------------------------------------------

! CHECK THAT AMOUNTS OF SO2 & DMS OXIDISED NOT GREATER THAN SO2 & DMS
! AVAILABLE AND SCALE INCREMENTS IF NECESSARY


!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length

      deltas_tot(i,j,k)=deltas_dry(i,j,k)+deltas_wet(i,j,k)                    &
        + deltas_wet_o3(i,j,k)

      IF (deltas_tot(i,j,k)  >   so2_array(i,j,k))  THEN
        scale_factor=so2_array(i,j,k)/deltas_tot(i,j,k)
        deltas_dry(i,j,k)=deltas_dry(i,j,k)*scale_factor
        deltas_wet(i,j,k)=deltas_wet(i,j,k)*scale_factor

        IF (l_sulpc_ozone) THEN
          deltas_wet_o3(i,j,k)=deltas_wet_o3(i,j,k)*scale_factor
        END IF

        deltas_tot(i,j,k) = so2_array(i,j,k)
      END IF

      IF ( l_sulpc_dms) THEN
        IF (deltas_dms(i,j,k)  >   dms_array(i,j,k)) THEN
          deltas_dms(i,j,k) = dms_array(i,j,k)
        END IF
      END IF
      ! Calculation of AIT_PRODUCTION, the amount of Aitken mode sulphate
      ! produced in a timestep:
      IF (l_sulpc_dms) THEN
        IF (l_sulpc_newdms .AND. l_sulpc_ozone) THEN
          ait_production(i,j,k) = psi(i,j,k)                                   &
            * (deltas_dry(i,j,k)                                               &
            + deltas_dms(i,j,k) * f_dms_to_so4(i,j,k))
        ELSE
          ait_production(i,j,k) = psi(i,j,k)                                   &
            * deltas_dry(i,j,k)
        END IF
      ELSE
        ait_production(i,j,k) = psi(i,j,k)                                     &
          * deltas_dry(i,j,k)
      END IF

      ! Fraction of Aitken mode transferred to accumulation mode in a
      ! timestep by mode merging:
      frac_trans(i,j,k) = k_merge * ait_production(i,j,k)
      ! Change in MMR of Aitken mode sulphate as a result of this transfer:
      deltas_merge(i,j,k) = frac_trans(i,j,k)                                  &
        * so4_ait_array(i,j,k)

      deltas_tot_ait(i,j,k)=deltas_diffuse(i,j,k)+                             &
        deltas_coag(i,j,k) + deltas_merge(i,j,k)

      IF ( deltas_tot_ait(i,j,k)  >   so4_ait_array(i,j,k) ) THEN
        scale_fact_ait=so4_ait_array(i,j,k)/deltas_tot_ait(i,j,k)
        deltas_diffuse(i,j,k)=deltas_diffuse(i,j,k)*scale_fact_ait
        deltas_coag(i,j,k)=deltas_coag(i,j,k)*scale_fact_ait
        deltas_merge(i,j,k) = deltas_merge(i,j,k) * scale_fact_ait
        deltas_tot_ait(i,j,k)=so4_ait_array(i,j,k)
      END IF

      !  UPDATE SO2, SO4 MODES AND DMS

      IF (l_sulpc_dms) THEN

        IF (l_sulpc_newdms .AND. l_sulpc_ozone) THEN

          so2_array(i,j,k)=                                                    &
            so2_array(i,j,k) +deltas_dms(i,j,k)*f_dms_to_so2(i,j,k)            &
            - deltas_tot(i,j,k)

          so4_ait_array(i,j,k)=                                                &
            so4_ait_array(i,j,k) - deltas_tot_ait(i,j,k)                       &
            + psi(i,j,k)*                                                      &
            (deltas_dry(i,j,k)+deltas_dms(i,j,k)*f_dms_to_so4(i,j,k))

          so4_acc_array(i,j,k)=so4_acc_array(i,j,k)                            &
            + deltas_evap(i,j,k)                                               &
            - deltas_nucl(i,j,k)                                               &
            + (1.0-psi(i,j,k))*                                                &
            (deltas_dry(i,j,k)+deltas_dms(i,j,k)*f_dms_to_so4(i,j,k) )         &
            + deltas_coag(i,j,k)                                               &
            + deltas_merge(i,j,k)

          so4_dis_array(i,j,k)=so4_dis_array(i,j,k)                            &
            + deltas_wet(i,j,k)                                                &
            - deltas_evap(i,j,k)                                               &
            + deltas_nucl(i,j,k)                                               &
            + deltas_diffuse(i,j,k)                                            &
            + deltas_wet_o3(i,j,k)

          dms_array(i,j,k)=dms_array(i,j,k) - deltas_dms(i,j,k)

          msa(i,j,k)=deltas_dms(i,j,k)*f_dms_to_msa(i,j,k)


        ELSE                 ! update for simple DMS scheme

          so2_array(i,j,k)=so2_array(i,j,k)                                    &
            + deltas_dms(i,j,k)*brat_so2                                       &
            - deltas_tot(i,j,k)


          so4_ait_array(i,j,k)=so4_ait_array(i,j,k)                            &
            + psi(i,j,k)*deltas_dry(i,j,k)                                     &
            - deltas_tot_ait(i,j,k)

          so4_acc_array(i,j,k)=so4_acc_array(i,j,k)                            &
            + deltas_evap(i,j,k)                                               &
            - deltas_nucl(i,j,k)                                               &
            + (1.0-psi(i,j,k))*deltas_dry(i,j,k)                               &
            + deltas_coag(i,j,k)                                               &
            + deltas_merge(i,j,k)

          so4_dis_array(i,j,k)=so4_dis_array(i,j,k)                            &
            + deltas_wet(i,j,k)                                                &
            - deltas_evap(i,j,k)                                               &
            + deltas_nucl(i,j,k)                                               &
            + deltas_diffuse(i,j,k)                                            &
            + deltas_wet_o3(i,j,k)

          dms_array(i,j,k) = dms_array(i,j,k) - deltas_dms(i,j,k)

          msa(i,j,k) = deltas_dms(i,j,k)*brat_msa


        END IF            ! END L_SULPC_NEWDMS TEST

      ELSE                   ! no DMS

        so2_array(i,j,k)=so2_array(i,j,k) - deltas_tot(i,j,k)


        so4_ait_array(i,j,k)=so4_ait_array(i,j,k)                              &
          + psi(i,j,k)*deltas_dry(i,j,k)                                       &
          - deltas_tot_ait(i,j,k)

        so4_acc_array(i,j,k)=so4_acc_array(i,j,k)                              &
          + deltas_evap(i,j,k)                                                 &
          - deltas_nucl(i,j,k)                                                 &
          + (1.0-psi(i,j,k))*deltas_dry(i,j,k)                                 &
          + deltas_coag(i,j,k)                                                 &
          + deltas_merge(i,j,k)

        so4_dis_array(i,j,k)=so4_dis_array(i,j,k)                              &
          + deltas_wet(i,j,k)                                                  &
          - deltas_evap(i,j,k)                                                 &
          + deltas_nucl(i,j,k)                                                 &
          + deltas_diffuse(i,j,k)                                              &
          + deltas_wet_o3(i,j,k)


      END IF                 ! End L_SULPC_DMS condition

    END DO
  END DO
END DO
!$OMP END DO

IF (l_sulpc_dms) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        dms(i,j,k)=dms_array(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

END IF

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length

      so2(i,j,k)=so2_array(i,j,k)
      so4_acc(i,j,k)=so4_acc_array(i,j,k)
      so4_ait(i,j,k)=so4_ait_array(i,j,k)
      so4_dis(i,j,k)=so4_dis_array(i,j,k)

    END DO
  END DO
END DO
!$OMP END DO NOWAIT


!--------------------------------------------------------------------
! 8.  DEPLETE  NH3.
!     IF ONLINE OXIDANTS ARE USED, DEPLETE OH, HO2, O3 AND H2O2.
!     OTHERWISE, DEPLETE H2O2, AND REPLENISH FROM HO2.
!--------------------------------------------------------------------

IF (l_sulpc_online_oxidants .AND. l_sulpc_2_way_coupling) THEN
!cdir collapse

!$OMP DO SCHEDULE(STATIC)
  DO k=1,wet_model_levels
    DO j=1,rows
      DO i=1,row_length
        delta_oh_tot(i,j,k) = delta_oh_dms(i,j,k)                              &
          + delta_oh_so2(i,j,k)
        IF (delta_oh_tot(i,j,k) .GT. oh_conc(i,j,k)) THEN
          scale_factor_oh = oh_conc(i,j,k) / delta_oh_tot(i,j,k)

          delta_oh_dms(i,j,k) = delta_oh_dms(i,j,k)                            &
            * scale_factor_oh
          delta_oh_so2(i,j,k) = delta_oh_so2(i,j,k)                            &
            * scale_factor_oh
        END IF
        delta_o3_tot(i,j,k) = delta_o3_dms(i,j,k)                              &
          + delta_o3_so2(i,j,k)

        IF (delta_o3_tot(i,j,k) .GT. o3(i,j,k)) THEN
          scale_factor_o3 =  o3(i,j,k) / delta_o3_tot(i,j,k)
          delta_o3_dms(i,j,k) = delta_o3_dms(i,j,k)                            &
            * scale_factor_o3
          delta_o3_so2(i,j,k) = delta_o3_so2(i,j,k)                            &
            * scale_factor_o3
          IF (delta_o3_tot(i,j,k) .EQ. 0.0) THEN
          END IF
        END IF

        IF (delta_ho2_dms(i,j,k) .GT. ho2_conc(i,j,k)) THEN
          delta_ho2_dms(i,j,k) = ho2_conc(i,j,k)
        END IF

        IF (delta_h2o2_so2(i,j,k) .GT. h2o2_mxr_array(i,j,k))                  &
          THEN
          delta_h2o2_so2(i,j,k) = h2o2_mxr_array(i,j,k)
        END IF

        oh_conc(i,j,k) = oh_conc(i,j,k) - delta_oh_dms(i,j,k)                  &
          - delta_oh_so2(i,j,k)
        ho2_conc(i,j,k) = ho2_conc(i,j,k) -  delta_ho2_dms(i,j,k)
        h2o2_mxr_array(i,j,k) = h2o2_mxr_array(i,j,k)                          &
          - delta_h2o2_so2(i,j,k)
        o3(i,j,k) = o3(i,j,k)  - delta_o3_dms(i,j,k)                           &
          - delta_o3_so2(i,j,k)
        IF (oh_conc(i,j,k) .LT. 0.0) oh_conc(i,j,k) = 0.0
        IF (ho2_conc(i,j,k) .LT. 0.0) ho2_conc(i,j,k) = 0.0
        IF (h2o2_mxr_array(i,j,k) .LT. 0.0)                                    &
          h2o2_mxr_array(i,j,k) = 0.0
        IF (o3(i,j,k) .LT. 0.0) o3(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO


ELSE


  ! Case of prescribed oxidants, or online oxidants without 2-way coupling

!cdir collapse

!$OMP DO SCHEDULE(STATIC)
  DO k=1,wet_model_levels
    DO j=1,rows
      DO i=1,row_length

        !  Deplete H2O2

        IF (deltas_wet(i,j,k)  >   0.0)  THEN
          h2o2_mxr_array(i,j,k)=                                               &
            h2o2_mxr_array(i,j,k)-deltas_wet(i,j,k)/relm_s_h2o2
        END IF

        !   Replenish H2O2_MXR field from HO2 field in dry part of grid box
        !   (whether or not oxidn has occurred) in daylight only.

        IF ( (h2o2_mxr_array(i,j,k)  <   h2o2(i,j,k)) .AND.                    &
          (cosza2d(i,j)  >   1.0e-06) )  THEN

          !   Calculate replenishment in concn, then convert to mix ratio
          !   Reaction rate K_HO2_HO2 is a fn of T, P, AND Q

          w_conc=q_array(i,j,k)*air_conc(i,j,k)*rmm_air/rmm_w !MCL/CC

          h2o2_con_to_mxr=rmm_h2o2/(rmm_air*air_conc(i,j,k))   !CC/MCL

          k_ho2_ho2=( k2*EXP(t2/t(i,j,k)) +                                    &
            air_conc(i,j,k)*k3*EXP(t3/t(i,j,k)) )                              &
            * (1.0 + w_conc*k4*EXP(t4/t(i,j,k)) )

          h2o2_rep = ho2_conc(i,j,k)*ho2_conc(i,j,k)*k_ho2_ho2                 &
            * tstep*clearf(i,j,k)

          h2o2_rep = h2o2_con_to_mxr * h2o2_rep

          !  Increment H2O2_MXR up to H2O2 value

          h2o2_mxr_array(i,j,k) = h2o2_mxr_array(i,j,k) + h2o2_rep

          IF ( h2o2_mxr_array(i,j,k)  >   h2o2(i,j,k) ) THEN
            h2o2_mxr_array(i,j,k) = h2o2(i,j,k)
          END IF

        END IF                     !End replenishment condition
      END DO
    END DO
  END DO
!$OMP END DO

END IF ! L_SULPC_ONLINE_OXIDANTS condition

! Deplete NH3 if ozone oxidation included and wet oxidn has occurred

IF (l_sulpc_ozone) THEN

!cdir collapse

!$OMP DO SCHEDULE(STATIC)
  DO k=1,wet_model_levels
    DO j=1,rows
      DO i=1,row_length

        IF ( (deltas_wet(i,j,k)  >   0.0) .OR.                                 &
          (deltas_wet_o3(i,j,k)  >   0.0) )  THEN

          nh3_dep(i,j,k) = (deltas_wet(i,j,k)+deltas_wet_o3(i,j,k))            &
            / relm_s_2n

          IF ( nh3_dep(i,j,k)  >   nh3_array(i,j,k) ) THEN
            nh3_dep(i,j,k) = nh3_array(i,j,k)
          END IF

          nh3_array(i,j,k) = nh3_array(i,j,k) - nh3_dep(i,j,k)

        END IF       ! end wet oxidn test

      END DO
    END DO
  END DO
!$OMP END DO


END IF       ! End L_SULPC_OZONE Test


!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      h2o2_mxr(i,j,k)=h2o2_mxr_array(i,j,k)
      nh3(i,j,k)=nh3_array(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

IF (l_sulpc_online_oxidants) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        h2o2(i,j,k)  = h2o2_mxr_array(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

END IF

!-----------------------------------------------------------------------


!$OMP END PARALLEL
IF (lhook) CALL dr_hook('SULPHR',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sulphr
END MODULE sulphr_mod
