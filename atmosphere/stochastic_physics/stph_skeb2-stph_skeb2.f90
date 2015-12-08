! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:  Stochastic Kinetic Energy Backscatter V2 (SKEB2)
!     This code provides a function routine for computing horizontally
!     non-divergent wind increments due to a supposed upscale energy
!     transfer from the mesoscale. The forcing is computed in the
!     spherical domain by expanding a streamfunction forcing function
!     in spherical harmonics. Each harmonic mode has a time variation
!     described by a first-order Markov process with in-built memory
!     corresponding to a typical mesoscale eddy turnover time (tau
!     ~ 0.5 -> 1 day).
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: Stochastic Physics
MODULE stph_skeb2_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE stph_skeb2(                                                  &
! in
           row_length, rows, n_rows, model_levels,                      &
           delta_phi, delta_lambda,                                     &
           rho_r2, u, v,                                                &
! in/out
           r_u, r_v,                                                    &
! STASH array sf
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
           stashwork35, first_atmstep_call)

 USE dynamics_input_mod, ONLY: l_endgame

! SKEB2 UMUI settings passed in via NameList READ
 USE stochastic_physics_run_mod, ONLY:                                  &
     n1, n2,skeb2_toplev, skeb2_botlev, br                              &
  ,  tot_backscat, tau, alphac, l_skeb2_psicdisp                        &
  ,  l_skeb2_psisdisp, l_skeb2_skeb1disp                                &
  ,  skeb2_cdisp, type_cape, type_mflx, l_skeb2_velpot                  &
  ,  skeb2_up_flux, skeb2_dwn_flux, skeb2_cape                          &
  ,  sdispfac, cdispfac, kdispfac, nsmooth                              &
  ,  offx_stph, offy_stph, mask_pdamp, mask_smooth                      &
  ,  l_skebsmooth_adv, l_skebprint, l_stphseed_read, l_stphseed_write   

! Sructure containing stochastic physics diagnostics
 USE stph_diag_mod, ONLY :                                              &
     strstphdiag
! Sructure containing field details for multi-variable SWAP_BOUNDS
 USE swapable_field_mod, ONLY:                                          &
     swapable_field_pointer_type
! Type needed for gather field of FFT run on multiple procs
 USE mpl, ONLY:                                                         &
     mpl_real

! Model level-height modules
 USE level_heights_mod,     ONLY:                                       &
     r_rho_levels, r_theta_levels, eta_theta_levels

! Bounds of arrays
 USE atm_fields_bounds_mod, ONLY:                                       &
     pdims, pdims_s,                                                    &
     udims, udims_s,                                                    &
     vdims, vdims_s,                                                    &
     tdims, wdims_s, stphdims_l

! Model grid trigonometry
 USE trignometric_mod, ONLY:                                            &
     cos_theta_latitude, sin_theta_latitude,                            &
     cos_v_latitude, sin_v_latitude

! Variables related to MPP
 USE proc_info_mod,     ONLY :                                          &
     global_row_length, global_rows
! Potential future use
!         mype=>me,                                                     &
!         nproc=>n_proc,                                                &
!         at_extremity, model_domain

 USE timestep_mod,      ONLY :                                          &
     timestep, timestep_number

! Routines for global sums
 USE global_2d_sums_mod, ONLY:                                          &
     global_2d_sums

! * tau:          Decorrelation time (~5.5 hrs in this case)
! * tot_backscat: Global-mean rate of energy backscatter in m**2.s**(-3)
! * br:           Backscatter ratio of dissipation Energy backscattered
! * alphac:       Updraught proportion of gridbox (0.2%)

 USE earth_constants_mod, ONLY: g, earth_radius

 USE conversions_mod, ONLY: pi

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim
 USE ereport_mod, ONLY : ereport
 USE PrintStatus_mod
 USE UM_ParVars
 USE backscatter_spectrum_mod, ONLY: backscatter_spectrum
 USE diagnostics_stph_mod, ONLY: diagnostics_stph
 USE skeb_forcing_mod,     ONLY: skeb_forcing
 USE skeb_smagorinsky_mod, ONLY: skeb_smagorinsky
 USE stph_closeinput_mod,  ONLY: stph_closeinput
 USE stph_closeoutput_mod, ONLY: stph_closeoutput
 USE stph_openinput_mod,   ONLY: stph_openinput
 USE stph_openoutput_mod,  ONLY: stph_openoutput
 USE stph_readentry2_mod,  ONLY: stph_readentry2
 USE stph_skeb1_mod,       ONLY: stph_skeb1
 USE stph_writeentry2_mod, ONLY: stph_writeentry2
 USE update_dpsidt_mod,    ONLY: update_dpsidt
 USE um_input_control_mod, ONLY: model_domain
                                 
 USE domain_params, ONLY: mt_global                            
 USE Submodel_Mod

 USE nlstcall_mod, ONLY : ldump

 USE chsunits_mod, ONLY : nunits

 IMPLICIT NONE

! This include contains: logical l_timer,
! This include contains STASH type declarations matching argsts.h
! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
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

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end
! This include contains: l_stphseed_read, l_stphseed_write
!LL  Comdeck: CCONTROL -------------------------------------------------
!LL
!LL  Purpose: COMMON block for top level switches and 2nd level switches
!LL           needed by the top level (C0) and 2nd level routines, but
!LL           not held in the history COMMON block.
!LL
!LLEND ---------------------------------------------------------------

!#include "cntlall.h"
! cntlgen.h was replaced by control/top_level/nlstgen_mod.F90
! #include "cntlgen.h"


 INTEGER, INTENT(IN) ::                                                 &
     row_length                                                         &
              ! local number of points on a row
 ,   rows                                                               &
              ! local number of rows for u
 ,   n_rows                                                             &
              ! local number of rows for v
 ,   model_levels
              ! model levels

 REAL, INTENT (IN) ::                                                   &
     delta_lambda                                                       &
              ! EW (x) grid spacing in radians
 ,   delta_phi                                                          &
              ! NS (y) grid spacing in radians
 ,   rho_r2(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end,                              &
            pdims_s%k_start:pdims_s%k_end)                              &
              ! density * square of radius of earth
 ,   u     (udims_s%i_start:udims_s%i_end,                              &
            udims_s%j_start:udims_s%j_end,                              &
            udims_s%k_start:udims_s%k_end)                              &
              ! main prog. u field
 ,   v     (vdims_s%i_start:vdims_s%i_end,                              &
            vdims_s%j_start:vdims_s%j_end,                              &
            vdims_s%k_start:vdims_s%k_end)
              ! main prog. v field

 REAL, INTENT (INOUT) ::                                                &
     r_u   (udims_s%i_start:udims_s%i_end,                              &
            udims_s%j_start:udims_s%j_end,                              &
            udims_s%k_start:udims_s%k_end)                              &
              ! main physics u increment (SKEB2 added to this)
 ,   r_v   (vdims_s%i_start:vdims_s%i_end,                              &
            vdims_s%j_start:vdims_s%j_end,                              &
            vdims_s%k_start:vdims_s%k_end)                              &
              ! main physics v increment (SKEB2 added to this)
 ,   stashwork35(*)
              ! Array containing requested sec35 STASH diagnostics

 REAL   ::                                                              &
     skeb2_urot(udims%i_start:udims%i_end,                              &
                udims%j_start:udims%j_end,                              &
                udims%k_start:udims%k_end)                              &
              ! rotational u increment from SKEB2
 ,   skeb2_vrot(vdims%i_start:vdims%i_end,                              &
                vdims%j_start:vdims%j_end,                              &
                vdims%k_start:vdims%k_end)                              &
              ! rotational v increment from SKEB2 on V-points
 ,   skeb2_udiv(udims%i_start:udims%i_end,                              &
                udims%j_start:udims%j_end,                              &
                udims%k_start:udims%k_end)                              &
              ! divergent u increment from SKEB2
 ,   skeb2_vdiv(vdims%i_start:vdims%i_end,                              &
                vdims%j_start:vdims%j_end,                              &
                vdims%k_start:vdims%k_end)
              ! divergent v increment from SKEB2 on V-points


! ----------------------------------------------------------------
!     NLIM is the spectral truncation. its possible values are
!     constrained by the fft routine which requires it to have
!     no prime factors > 19 and the total number of primes
!     factors (including repetitions) must not exceed 20.

!     NLIM, NLAT are given values based on  model dimensions
! ----------------------------------------------------------------

 INTEGER, SAVE ::                                                       &
     nlim                                                               &
              ! Spectral equivalent global row_length
 ,   nlat                                                               &
              ! Spectral equivalent number of latitude pairs * 2
 ,   level2km
              ! model level at ~2km (=10 for L38)

 REAL, SAVE ::                                                          &
     max_tot_backscat                                                   &
              ! Maximum allowed backscatter
 ,   crit_tot_backscat                                                  &
              ! Critical backscatter amount requiring damping
 ,   logscale                                                           &
              ! Multiplication factor for LOG10 arguments (see levfac)
 ,   levfac
              ! Storage array for LOG10(model_level), used to reduce

! Allocatable work arrays
 REAL, ALLOCATABLE ::                                                   &
     my_dpsidt    (:,:)                                                 &
              ! Local PE d(psi)/d(t) -- in grid-space
 ,   my_dpsidtc   (:,:)                                                 &
              ! Local PE Fourier Space COS coeffs
 ,   my_dpsidts   (:,:)                                                 &
              ! Local PE Fourier Space SIN coeffs
 ,   my_dpsidtr   (:,:)                                                 &
              ! Local PE Fourier Space Radius
 ,   my_phi       (:,:)                                                 &
              ! Local PE Fourier Space degree  
 ,   my_phishft   (:,:)                                                 &
              ! Local PE Fourier Space rotation
 ,   g_rand_nums  (:,:,:)                                               &
              ! Global version of rand_nums, used to generate global
              ! random field which is SCATTERED to each PE to ensure
              ! reproducibility across different PE configurations.
 ,   work_p_halo  (:,:,:)                                               &
              ! Dummy variable to swap bounds and hold data on p-grid
 ,   work_z_halo  (:,:,:)                                               &
              ! Dummy variable to swap bounds and hold data on vort-grid
 ,   vert_int_work(:,:)                                                 &
              ! vertically integrated flux W.m^-2 (working array)
 ,   dpsidt       (:)                                                   &
              ! d(psi)/d(t) used for Markov process integration
 ,   sum_temp     (:)                                                   &
              ! Temporary array to hold row sums for rvecsumr
 ,   psif0        (:,:)                                                 &
              ! streamfunction forcing field with halo
 ,   sm_tot_disp0 (:,:,:)
              ! holding variable for smoothing of modulating field

! Allocatable variables (needing to be saved between timesteps)
 REAL, ALLOCATABLE, SAVE ::                                             &
     dpsidtc (:,:)                                                      &
              ! 2D version of d(psi)/d(t) COS coeffs in Fourier
 ,   dpsidts (:,:)                                                      &
              ! 2D version of d(psi)/d(t) SIN coeffs in Fourier
 ,   dx_theta(:,:)                                                      &
              ! Delta-X on u points
 ,   dy_theta(:,:)                                                      &
              ! Delta-Y on u points
 ,   dx_v    (:,:)                                                      &
              ! Delta-X on v points
 ,   dy_v    (:,:)                                                      &
              ! Delta-Y on v points
 ,   darea   (:,:)                                                      &
              ! Horizontal area of gridboxes
 ,   c_dadt  (:,:)                                                      &
              ! Horizontal area of gridboxes * timestep (inverted)
 ,   gspect  (:)                                                        &
              ! Wave-number dependent noise amplitude
 ,   psi     (:,:)
              ! Global single-level version of SF pattern psif

 INTEGER, ALLOCATABLE ::                                                &
     iranseed(:)
              ! Random seed size used for positioning read statement
              ! of dpsidtc/s from the seed file (unit=149)

! Local arrays and scalars for dissipation calculations
 REAL  ::                                                               &
     psif     (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! streamfunction forcing field (E=tot_backscat)
 ,   m_psif   (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! Modulated stream_function forcing field
 ,   rand_nums(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! 3 dimensional random field
 ,   mu                                                                 &
              ! mu is sin(latitude)
 ,   alpha                                                              &
              ! autoregressive process parameter related to tau
 ,   delr_rho (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! thickness between r_rho_levels
 ,   mass     (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! mass (kg) of 3D grid boxes
 ,   mass_2d  (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end)                               &
              ! mass (kg) of column
 ,   kr                                                                 &
              ! Ratio in z direction
 ,   smallp = 1.e-08                                                    &
              ! Small real number
 ,   local_crit_tot_backscat                                            &
              ! crit_tot_backscat at local level
 ,   gltotke                                                            &
              ! Used for global sum of quantity for printing
 ,   gltotke_tmp(1)
              ! Array in global sum

! Loop counters
 INTEGER ::                                                             &
     i                                                                  &
              ! loop index over x direction
 ,   ilat                                                               &
              ! loop index over latitude
 ,   j                                                                  &
              ! loop index over y direction
 ,   k                                                                  &
              ! loop index over z direction
 ,   m                                                                  &
              ! loop index over EW wavespace
 ,   n                                                                  &
              ! loop index over NS wavespace
 ,   ii                                                                 &
              ! loop index ii inside i-loop
 ,   jj                                                                 &
              ! loop index jj inside j-loop
 ,   i1                                                                 &
              ! loop index i1 inside i-loop
 ,   i2                                                                 &
              ! loop index i2 inside j-loop
 ,   j1                                                                 &
              ! loop index j1 inside j-loop
 ,   j2                                                                 &
              ! loop index j2 inside j-loop
 ,   ip1                                                                &
              ! i plus one
 ,   im1                                                                &
              ! i minus one
 ,   jp1                                                                &
              ! j plus one
 ,   jm1                                                                &
              ! j minus one
 ,   ju                                                                 &
              ! j on u-grid (limit to rows if at north pole)
 ,   jv                                                                 &
              ! j on v-grid (limit to n_rows if at north pole)
 ,   ismooth                                                            &
              ! loop index over smoothing iterations
 ,   icode = 0                                                          &
              ! Return code for error reporting
 ,   info  = 0                                                          &
              ! return code from GC stuff
 ,   ifield                                                             &
              ! number of fields used in swap_bounds_mv
 ,   firsttimestep                                                      &
              ! Set to 1 if first call (otherwise 2)
 ,   sndcount                                                           &
              ! Size of real array bcast to all proc's using gc_rbcast
 ,   tam
              ! scalar holding size of the random seed

 REAL ::                                                                &
     rho      (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! density on each model level
 ,   sm_tot_disp(pdims%i_start:pdims%i_end,                             &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! smoothed total instantaneous dissipation rate
 ,   diff_flux(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! difference of the up-dwn mass flux
 ,   sdisp    (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! dissipitation in SMAGORINSKY bit
 ,   cdisp    (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! dissipitation from Convection
 ,   kdisp    (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)                               &
              ! dissipitation in SKEB1-type (KE-based) calculation
 ,   cape_kk  (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end)
              ! local modified cape field

 LOGICAL ::                                                             &
     l_latwgt = .FALSE.                                                 &
              ! Set to TRUE will apply sin(lat)/cos(lat) weights to
              ! SKEB2 rot/div wind increments (default=.FALSE.)
 ,   first_atmstep_call
              ! Used to set firsttimestep == 1 or 2
              ! Is true for first step of: NRUN and each CRUN

 CHARACTER(LEN=256) :: cmessage
              ! OUT error message

! Parameters
 INTEGER, PARAMETER :: zero = 0 ! used for identifying zero'th PE
 CHARACTER(LEN=*)       :: routinename
 PARAMETER( routinename='stph_skeb2' )

! Local targets
 REAL, TARGET ::                                                        &
     skeb2_urot0(udims_s%i_start:udims_s%i_end,                         &
                 udims_s%j_start:udims_s%j_end,                         &
                 udims_s%k_start:udims_s%k_end)                         &
              ! rotational u increment on U-points
 ,   skeb2_vrot0(vdims_s%i_start:vdims_s%i_end,                         &
                 vdims_s%j_start:vdims_s%j_end,                         &
                 vdims_s%k_start:vdims_s%k_end)                         &
              !  rotational v increment on V-points (before smoothing)
 ,   skeb2_udiv0(udims_s%i_start:udims_s%i_end,                         &
                 udims_s%j_start:udims_s%j_end,                         &
                 udims_s%k_start:udims_s%k_end)                         &
              ! divergent u increment on U-points
 ,   skeb2_vdiv0(vdims_s%i_start:vdims_s%i_end,                         &
                 vdims_s%j_start:vdims_s%j_end,                         &
                 vdims_s%k_start:vdims_s%k_end)
              ! divergent v increment on V-points (before smoothing)

! Local variables for Latitude decomposition in FFT CALL (looping
! over latitude)
 INTEGER :: fft_rows(0:nproc_max), my_rows
 INTEGER :: rowcounts(0:nproc_max)
 INTEGER :: displs(0:nproc_max), my_displs

! Declaration of Stochastic Physics diagnostics.
 TYPE (strstphdiag)                 :: stph_diag
! Declaration of fields used in swap_bounds_mv
 TYPE (swapable_field_pointer_type) :: fields_to_swap(8)

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

! NOTE on Arakawa C-grid
!      Smagorinsky and Convection estimates of Energy Dissipation
!      are stored on P-points. Final streamfunction forcing (S) used
!      to calculate winds needs to be interpolated to the vorticty
!      (Z) points for rot winds to avoid the U|V grids being swapped.
!      Rot winds are from PSI on Z-points: (u;v) = (-dS/dy; dS/dx)
!      Div winds are from PSI on P-points: (u;v) = ( dS/dx; dS/dy)
!
!  LON=  0   1/2dX  dX   3/2dX
!        |     :     |     :   
!        |     :     |     :   
!      --V-----Z-----V-----Z-- 3/2dY
!        |     :     |     :   
!        |     :     |     :   
!      ==P=====U=====P=====U== dY
!        |     :     |     :   
!        |     :     |     :   
!      --V-----Z-----V-----Z-- 1/2dY
!        |     :     |     :   
!        |     :     |     :   
!      ==P=====U=====P=====U== SOUTH POLE
!
! ------------------------------------------------------------------
!      END OF VARIABLE DECLARATIONS - START OF THE CODE
! ------------------------------------------------------------------

 IF (lhook) CALL dr_hook('STPH_SKEB2',zhook_in,zhook_handle)

! Error trap if this is a LAM configuration (SKEB2 should not be called)
 IF (model_domain /= mt_global) THEN
 WRITE(6,'("**ERROR**: SKEB2 not available in a Limited Area Model")')
 WRITE(6,'("  Section 35: check UMUI settings")')
 icode   = 1
 WRITE (cmessage,*) 'STPH_SKEB2: Not configured to run in LAMS. ',      &
   &                'Turn scheme off in UMUI.'

 
 CALL ereport(routinename, icode, cmessage)
 END IF  ! .NOT. GLOBAL

! Check bounds of SKEB2_TOPLEV and SKEB2_BOTLEV
! SKEB2 cannot extend to top and bottom model levels because of
! calculations involving levels (k-1) and (k+1) in the code
! and physical realism of the thickness of the backscattered layer
 IF  (skeb2_toplev >= model_levels) THEN
   WRITE(6,'("**ERROR**: SKEB2 TOP LEVEL = OR > MODEL TOP")')
   WRITE(6,'("  Section 35: check UMUI settings")')
   icode   = 1
   WRITE (cmessage,*) 'STPH_SKEB2: skeb2_toplev (',skeb2_toplev,') is ',&
     &                '>= model_levels (',model_levels,')'

   CALL ereport(routinename, icode, cmessage)
 END IF
 IF  (skeb2_botlev <= 1) THEN
   WRITE(6,'("ERROR: SKEB2 BOTTOM LEVEL = OR < MODEL BASE")')
   WRITE(6,'("  Section 35: check UMUI settings")')
   icode   = 1
   WRITE (cmessage,*) 'STPH_SKEB2: skeb2_botlev (',skeb2_botlev,') is ',&
     &                '<= 1'

   CALL ereport(routinename, icode, cmessage)
 END IF
 IF  (skeb2_toplev < skeb2_botlev + 3) THEN
   WRITE(6,'("ERROR: SKEB2 VERTICAL EXTENT TOO THIN")')
   WRITE(6,'("       SKEB2 TOP LEVEL < BOTTOM LEVEL+4")')
   WRITE(6,'("  Section 35: check UMUI settings")')
   icode   = 1
   WRITE (cmessage,*) 'STPH_SKEB2: Need at least 3 levels between ',    &
     &                'skeb2_toplev (',skeb2_toplev,') and ',           &
     &                'skeb2_botlev (',skeb2_botlev,')'

   CALL ereport(routinename, icode, cmessage)
 END IF


! Settings required only at the first time-step
! firsttimestep = 1 at first call (LOGICAL first_atmstep_call = T)
! firsttimestep = 2 otherwise (LOGICAL first_atmstep_call = F)
 firsttimestep = 2
 IF  (first_atmstep_call) firsttimestep = 1
 IF  (firsttimestep == 1) THEN

! Initialize variables from UM data for the sph.harm calculations
! For the UM, not being an spectral model, NLIM:

   IF  (mod(global_row_length,2) == 0) THEN
     nlim=global_row_length/2
   ELSE
     nlim=(global_row_length+1)/2
   END IF

! nlat should be equal to global_rows (and even)
! Uses SCATTER_FIELD for psi => psif
   IF  (mod(global_rows,2) == 0) THEN
     nlat = global_rows
     IF  (.NOT.ALLOCATED(psi)) THEN
       ALLOCATE (psi(1:2*nlim,nlat))
       DO j = 1, nlat
         DO i = 1, 2*nlim
           psi(i,j) = 0.0
         END DO
       END DO
     END IF
   ELSE
     nlat = global_rows-1
     IF  (.NOT.ALLOCATED(psi)) THEN
       ALLOCATE (psi(1:2*nlim,nlat+1))
       DO j = 1, nlat+1
         DO i = 1, 2*nlim
           psi(i,j) = 0.0
         END DO
       END DO
     END IF
   END IF

! The values of these variables are saved so I only need
! to calculate/initialise them in the first time-step
   IF  (.NOT.ALLOCATED(gspect)) THEN
     ALLOCATE (gspect(nlim))
   END IF
   IF  (.NOT.ALLOCATED(dpsidtc)) THEN
      ALLOCATE (dpsidtc(0:n2,n1:n2))
   END IF
   IF  (.NOT.ALLOCATED(dpsidts)) THEN
      ALLOCATE (dpsidts(0:n2,n1:n2))
   END IF
   IF  (.NOT.ALLOCATED(dx_theta)) THEN
      ALLOCATE (dx_theta(udims%i_start:udims%i_end,                     &
                         udims%j_start:udims%j_end))
   END IF
   IF  (.NOT.ALLOCATED(dy_theta)) THEN
      ALLOCATE (dy_theta(udims%i_start:udims%i_end,                     &
                         udims%j_start:udims%j_end))
   END IF
   IF  (.NOT.ALLOCATED(dx_v)) THEN
      ALLOCATE (dx_v    (vdims%i_start:vdims%i_end,                     &
                         vdims%j_start:vdims%j_end))
   END IF
   IF  (.NOT.ALLOCATED(dy_v)) THEN
      ALLOCATE (dy_v    (vdims%i_start:vdims%i_end,                     &
                         vdims%j_start:vdims%j_end))
   END IF
   IF  (.NOT.ALLOCATED(darea)) THEN
      ALLOCATE (darea   (pdims%i_start:pdims%i_end,                     &
                         pdims%j_start:pdims%j_end))
   END IF
   IF  (.NOT.ALLOCATED(c_dadt)) THEN
      ALLOCATE (c_dadt  (pdims%i_start:pdims%i_end,                     &
                         pdims%j_start:pdims%j_end))
   END IF

!  These arrays have sections that may not be properly initialised
   DO n = n1, n2
     DO m = 0, n2
       dpsidtc(m, n) = 0.0
       dpsidts(m, n) = 0.0
     END DO
   END DO

!  Calculate and store grid increment values
   IF (l_latwgt) THEN
     ! Weight wind incr by sin/cos latitude to decrease impact of
     ! rot-wind in tropics and div-wind in high latitudes
     ! Note: dx = dlam * a * cos(lat), so "cos(lat)" term cancels out
     DO j = udims%j_start, udims%j_end
       DO i = udims%i_start, udims%i_end
         dx_theta(i,j) = delta_lambda * earth_radius
         dy_theta(i,j) = (delta_phi * earth_radius)/                    &
                          sin_theta_latitude(1,j)
       END DO
     END DO
     DO j = vdims%j_start, vdims%j_end
       DO i = vdims%i_start, vdims%i_end
         dx_v(i,j) = (delta_lambda * earth_radius * cos_v_latitude(i,j))&
                     / sin_v_latitude(i,j)
         dy_v(i,j) = (delta_phi * earth_radius)/ cos_v_latitude(i,j)
       END DO
     END DO
   ELSE
     ! No wind incr weighting (this is the default)
     DO j = udims%j_start, udims%j_end
       DO i = udims%i_start, udims%i_end
         dx_theta(i,j) = delta_lambda * earth_radius
         dy_theta(i,j) = delta_phi * earth_radius
       END DO
     END DO
     DO j = vdims%j_start, vdims%j_end
       DO i = vdims%i_start, vdims%i_end
         dx_v(i,j) = delta_lambda * earth_radius
         dy_v(i,j) = delta_phi * earth_radius
       END DO
     END DO
   END IF   ! l_latwgt

!  Limit ratio between dX & dY to 1:3 to eliminate generating large
!  wind increments near the poles
!  WHERE(dx_theta < 0.333*dy_theta) dx_theta = 0.333*dy_theta
!  WHERE(dx_v < 0.333*dy_v) dx_v = 0.333*dy_v

!  Horizontal area of gridbox (weighted by latitude)
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       darea(i,j) = delta_lambda * earth_radius * delta_phi *           &
                    earth_radius * cos_theta_latitude(i,j)
     END DO
   END DO

!  Set minimum grid-box area (to avoid divide-by-zero at poles)
   WHERE(darea < 1.e3) darea = 1.e3

!  Invert array so that it can be used as a multiplication factor
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       c_dadt(i,j) = 1./(darea(i,j)*timestep)
     END DO
   END DO

   ! Keep a record of the settings for this run in the PE output files
   IF (printstatus  >=  prstatus_normal) THEN
     WRITE(6,*) 'L_SKEB2_PSISDISP = ', l_skeb2_psisdisp
     WRITE(6,*) 'L_SKEB2_PSICDISP = ', l_skeb2_psicdisp
     WRITE(6,*) 'L_SKEB2_SKEB1DISP = ', l_skeb2_skeb1disp
     WRITE(6,'(A,I3)') 'SKEB2_CDISP type = ', skeb2_cdisp
     WRITE(6,'(A,L3)') 'L_SKEB2_VELPOT = ', l_skeb2_velpot
     IF (l_skeb2_psicdisp) THEN
       IF (skeb2_cdisp.NE.type_cape .AND. skeb2_cdisp.NE.type_mflx) THEN
         WRITE(6,'("**warning**: skeb2 convective dissipation is on, '//&
             'but no valid scheme is selected")')
         WRITE(6,'("  Section 35: check Rose settings")')
         icode    = -1
         WRITE (cmessage,*) 'STPH_SKEB2: no convective dissipation ',   &
                            'scheme is selected - will return zeroes'
         CALL ereport(routinename, icode, cmessage)
       END IF
     END IF
   END IF

  ! Find 1st model theta level above 2km
  ! Model levels are constant in time
   DO k = wdims_s%k_start, wdims_s%k_end
     level2km = k
     IF (eta_theta_levels(k) * (r_theta_levels(1,1,wdims_s%k_end) -     &
         Earth_radius) > 2000.) EXIT
   END DO

  ! Convert range [mod_lev=1; mod_lev=level2km] => [1; 10]
   logscale = 10./level2km
   IF (printstatus  ==  prstatus_diag) THEN
     WRITE(6,'(" ---------------------------------------- ")')
     WRITE(6,'(" Vertical ramp decrease of wind increment ")')
     WRITE(6,'(" Level of 2km = ",I4)') level2km
     WRITE(6,'(" Logscale = ",ES12.4)') logscale
   END IF

 END IF ! firsttimestep ==1

! Allocate work variables
 IF  (.NOT.ALLOCATED(dpsidt)) THEN
   ALLOCATE (dpsidt(0:2*nlim+1))
 END IF
 IF  (.NOT.ALLOCATED(vert_int_work)) THEN
   ALLOCATE (vert_int_work(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end))
 END IF
 IF  (.NOT.ALLOCATED(work_p_halo)) THEN
   ALLOCATE (work_p_halo(  pdims_s%i_start:pdims_s%i_end,               &
                           pdims_s%j_start:pdims_s%j_end,               &
                           pdims_s%k_start:pdims_s%k_end))
 END IF
 IF  (.NOT.ALLOCATED(work_z_halo)) THEN
   ALLOCATE (work_z_halo(  vdims_s%i_start:vdims_s%i_end,               &
                           vdims_s%j_start:vdims_s%j_end,               &
                           vdims_s%k_start:vdims_s%k_end))
 END IF

! Prepare STASH diagnostics where requested
! Allocate STASH diagnostic arrays using Structure TYPE stph_diag
 stph_diag%l_skeb2_u              = sf(1,35)
 stph_diag%l_skeb2_v              = sf(2,35)
 stph_diag%l_skeb2_u_incr         = sf(3,35)
 stph_diag%l_skeb2_v_incr         = sf(4,35)
 stph_diag%l_skeb2_u_rot          = sf(5,35)
 stph_diag%l_skeb2_v_rot          = sf(6,35)
 stph_diag%l_skeb2_u_div          = sf(7,35)
 stph_diag%l_skeb2_v_div          = sf(8,35)
 stph_diag%l_skeb2_disp_smag      = sf(9,35)
 stph_diag%l_skeb2_disp_conv      = sf(10,35)
 stph_diag%l_skeb2_disp_skeb1     = sf(11,35)
 stph_diag%l_skeb2_smodfield      = sf(12,35)
 stph_diag%l_skeb2_streamfunction = sf(13,35)
 stph_diag%l_skeb2_random_pattern = sf(14,35)
 stph_diag%l_skeb2_ke_psif        = sf(15,35)
 stph_diag%l_skeb2_ke_sdisp       = sf(16,35)
 stph_diag%l_skeb2_ke_cdisp       = sf(17,35)
 stph_diag%l_skeb2_ke_kdisp       = sf(18,35)
 stph_diag%l_skeb2_ke_m_psif      = sf(19,35)
 stph_diag%l_skeb2_ke_prewindincr = sf(20,35)
 stph_diag%l_skeb2_ke_windincr    = sf(21,35)
 stph_diag%l_skeb2_ke_postwindincr= sf(22,35)
    
 IF (sf(0,35)) THEN

! u after skeb2
   IF (stph_diag%L_skeb2_u) then
     ALLOCATE(stph_diag%skeb2_u(udims%i_start:udims%i_end,              &
                                udims%j_start:udims%j_end,              &
                                udims%k_start:udims%k_end))
   ELSE
     ! Code to allocate unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_u(1,1,1))
     stph_diag%skeb2_u(:,:,:) = 0.0
   END IF

! v after skeb2
   IF (stph_diag%L_skeb2_v) then
     ALLOCATE(stph_diag%skeb2_v(vdims%i_start:vdims%i_end,              &
                                vdims%j_start:vdims%j_end,              &
                                vdims%k_start:vdims%k_end))
   ELSE
     ! Code to allocate unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_v(1,1,1))
     stph_diag%skeb2_v(:,:,:) = 0.0
   END IF

! u increment diagnostic (store increment before SKEB2)
   IF (stph_diag%l_skeb2_u_incr) THEN
     ALLOCATE(stph_diag%skeb2_u_incr(udims%i_start:udims%i_end,         &
                                     udims%j_start:udims%j_end,         &
                                     udims%k_start:udims%k_end))
     DO k = udims%k_start, udims%k_end
       DO j = udims%j_start, udims%j_end
         DO i = udims%i_start, udims%i_end
           stph_diag%skeb2_u_incr(i,j,k) = r_u(i,j,k)
         END DO
       END DO
     END DO
   ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
     ALLOCATE(stph_diag%skeb2_u_incr(1,1,1))
     stph_diag%skeb2_u_incr(:,:,:) = 0.0
   END IF

! v increment diagnostic (store increment before SKEB2)
   IF (stph_diag%L_skeb2_v_incr) THEN
     ALLOCATE(stph_diag%skeb2_v_incr(vdims%i_start:vdims%i_end,         &
                                     vdims%j_start:vdims%j_end,         &
                                     vdims%k_start:vdims%k_end))
     DO k = vdims%k_start, vdims%k_end
       DO j = vdims%j_start, vdims%j_end
         DO i = vdims%i_start, vdims%i_end
           stph_diag%skeb2_v_incr(i,j,k) = r_v(i,j,k)
         END DO
       END DO
     END Do
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_v_incr(1,1,1))
     stph_diag%skeb2_v_incr(:,:,:) = 0.0
   END IF

! rotational u increments from SKEB2 
   IF (stph_diag%L_skeb2_u_rot) THEN
     ALLOCATE(stph_diag%skeb2_u_rot(udims%i_start:udims%i_end,          &
                                    udims%j_start:udims%j_end,          &
                                    udims%k_start:udims%k_end))
     stph_diag%skeb2_u_rot(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_u_rot(1,1,1))
     stph_diag%skeb2_u_rot(:,:,:) = 0.0
   END IF

! rotational v increments from SKEB2
   IF (stph_diag%L_skeb2_v_rot) THEN
     ALLOCATE(stph_diag%skeb2_v_rot(vdims%i_start:vdims%i_end,          &
                                    vdims%j_start:vdims%j_end,          &
                                    vdims%k_start:vdims%k_end))
     stph_diag%skeb2_v_rot(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_v_rot(1,1,1))
     stph_diag%skeb2_v_rot(:,:,:) = 0.0
   END IF

! divergent u increments from SKEB2 
   IF (stph_diag%L_skeb2_u_div) THEN
     ALLOCATE(stph_diag%skeb2_u_div(udims%i_start:udims%i_end,          &
                                    udims%j_start:udims%j_end,          &
                                    udims%k_start:udims%k_end))
     stph_diag%skeb2_u_div(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_u_div(1,1,1))
     stph_diag%skeb2_u_div(:,:,:) = 0.0
   END IF

! divergent v increments from SKEB2 
   IF (stph_diag%L_skeb2_v_div) THEN
     ALLOCATE(stph_diag%skeb2_v_div(vdims%i_start:vdims%i_end,          &
                                    vdims%j_start:vdims%j_end,          &
                                    vdims%k_start:vdims%k_end))
     stph_diag%skeb2_v_div(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_v_div(1,1,1))
     stph_diag%skeb2_v_div(:,:,:) = 0.0
   END IF

! SKEB2: dissipation field from smagorinsky
   IF (stph_diag%L_skeb2_disp_smag) THEN
     ALLOCATE(stph_diag%skeb2_disp_smag(pdims%i_start:pdims%i_end,      &
                                        pdims%j_start:pdims%j_end,      &
                                        pdims%k_start:pdims%k_end))
     stph_diag%skeb2_disp_smag(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_disp_smag(1,1,1))
     stph_diag%skeb2_disp_smag(:,:,:) = 0.0
   END IF

! SKEB2: dissipation field from convection
   IF (stph_diag%L_skeb2_disp_conv) THEN
     ALLOCATE(stph_diag%skeb2_disp_conv(pdims%i_start:pdims%i_end,      &
                                        pdims%j_start:pdims%j_end,      &
                                        pdims%k_start:pdims%k_end))
     stph_diag%skeb2_disp_conv(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_disp_conv(1,1,1))
     stph_diag%skeb2_disp_conv(:,:,:) = 0.0
   END IF

! SKEB2: dissipation field from SKEB1 KE dissipation
   IF (stph_diag%L_skeb2_disp_skeb1) THEN
     ALLOCATE(stph_diag%skeb2_disp_skeb1(pdims%i_start:pdims%i_end,     &
                                         pdims%j_start:pdims%j_end,     &
                                         pdims%k_start:pdims%k_end))
     stph_diag%skeb2_disp_skeb1(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_disp_skeb1(1,1,1))
     stph_diag%skeb2_disp_skeb1(:,:,:) = 0.0
   END IF

! SKEB2: Smoothed Modulating field
   IF (stph_diag%L_skeb2_smodfield) THEN
     ALLOCATE(stph_diag%skeb2_smodfield(pdims%i_start:pdims%i_end,      &
                                        pdims%j_start:pdims%j_end,      &
                                        pdims%k_start:pdims%k_end))
     stph_diag%skeb2_smodfield(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_smodfield(1,1,1))
     stph_diag%skeb2_smodfield(:,:,:) = 0.0
   END IF

! SKEB2: final stream function forcing field
   IF (stph_diag%L_skeb2_streamfunction) THEN
     ALLOCATE(stph_diag%skeb2_streamfunction(pdims%i_start:pdims%i_end, &
                                             pdims%j_start:pdims%j_end, &
                                             pdims%k_start:pdims%k_end))
     stph_diag%skeb2_streamfunction(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_streamfunction(1,1,1))
     stph_diag%skeb2_streamfunction(:,:,:) = 0.0
   END IF

! SKEB2: initial random pattern
   IF (stph_diag%L_skeb2_random_pattern) THEN
     ALLOCATE(stph_diag%skeb2_random_pattern(pdims%i_start:pdims%i_end, &
                                             pdims%j_start:pdims%j_end, &
                                             pdims%k_start:pdims%k_end))
     stph_diag%skeb2_random_pattern(:,:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_random_pattern(1,1,1))
     stph_diag%skeb2_random_pattern(:,:,:) = 0.0
   END IF

! SKEB2: Vert Integ. KE of initial SF forcing
   IF (stph_diag%L_skeb2_ke_psif) THEN
     ALLOCATE(stph_diag%skeb2_ke_psif(pdims%i_start:pdims%i_end,        &
                                      pdims%j_start:pdims%j_end))
     stph_diag%skeb2_KE_psif(:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_ke_psif(1,1))
     stph_diag%skeb2_KE_psif(:,:) = 0.0
   END IF

! SKEB2: Vert Integ. KE of numerical diss
   IF (stph_diag%L_skeb2_ke_sdisp) THEN
     ALLOCATE(stph_diag%skeb2_ke_sdisp(pdims%i_start:pdims%i_end,       &
                                       pdims%j_start:pdims%j_end))
     stph_diag%skeb2_KE_sdisp(:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_ke_sdisp(1,1))
     stph_diag%skeb2_KE_sdisp(:,:) = 0.0
   END IF

! SKEB2: Vert Integ. KE of convection diss
   IF (stph_diag%L_skeb2_ke_cdisp) THEN
     ALLOCATE(stph_diag%skeb2_ke_cdisp(pdims%i_start:pdims%i_end,       &
                                       pdims%j_start:pdims%j_end))
     stph_diag%skeb2_KE_cdisp(:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_ke_cdisp(1,1))
     stph_diag%skeb2_KE_cdisp(:,:) = 0.0
   END IF

! SKEB2: Vert Integ. KE of 2nd convection diss
   IF (stph_diag%L_skeb2_ke_kdisp) THEN
     ALLOCATE(stph_diag%skeb2_ke_kdisp(pdims%i_start:pdims%i_end,       &
                                       pdims%j_start:pdims%j_end))
     stph_diag%skeb2_KE_kdisp(:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_ke_kdisp(1,1))
     stph_diag%skeb2_KE_kdisp(:,:) = 0.0
   END IF

! SKEB2: Vert Integ. KE of modulated SF forcing
   IF (stph_diag%L_skeb2_ke_m_psif) THEN
     ALLOCATE(stph_diag%skeb2_ke_m_psif(pdims%i_start:pdims%i_end,      &
                                        pdims%j_start:pdims%j_end))
     stph_diag%skeb2_KE_m_psif(:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_ke_m_psif(1,1))
     stph_diag%skeb2_ke_m_psif(:,:) = 0.0
   END IF

! SKEB2: Vert Integ. KE of total wind incr before SKEB2
   IF (stph_diag%L_skeb2_ke_prewindincr) THEN
     ALLOCATE(stph_diag%skeb2_ke_prewindincr(pdims%i_start:pdims%i_end, &
                                             pdims%j_start:pdims%j_end))
     stph_diag%skeb2_ke_prewindincr(:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_ke_prewindincr(1,1))
     stph_diag%skeb2_ke_prewindincr(:,:) = 0.0
   END IF

! SKEB2: Vert Integ. KE of wind incr from SKEB2
   IF (stph_diag%L_skeb2_ke_windincr) THEN
     ALLOCATE(stph_diag%skeb2_ke_windincr(pdims%i_start:pdims%i_end,    &
                                          pdims%j_start:pdims%j_end))
     stph_diag%skeb2_ke_windincr(:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_ke_windincr(1,1))
     stph_diag%skeb2_ke_windincr(:,:) = 0.0
   END IF

! SKEB2: Vert Integ. KE of total wind incr after SKEB2
   IF (stph_diag%L_skeb2_ke_postwindincr) THEN
     ALLOCATE(stph_diag%skeb2_ke_postwindincr(                          &
                                      pdims%i_start:pdims%i_end,        &
                                      pdims%j_start:pdims%j_end))
     stph_diag%skeb2_ke_postwindincr(:,:) = 0.0
   ELSE
     ! Code to ALLOCATE unity arrays when not
     ! used (for portability)
     ALLOCATE(stph_diag%skeb2_ke_postwindincr(1,1))
     stph_diag%skeb2_ke_postwindincr(:,:) = 0.0
   END IF

 END IF ! SF(0,35)

! -----------------------------------------------------------------
!  This bit is to calculate the spherical harm. forcing
!  compute the spectral weighting function g(n) - which contains
!  the power law  "chi= n**(-1.54)" that Glenn saw in the CRM -
!  new recommended value (used in SPBS at ECMWF) is -1.54
!  **************
!  A vertical variation of the pattern is done by changing the wave
!  coefficients with height at an angular speed dependent on wave
!  number. The coefficients evolve in time according to a 1st-order
!  Markov process expression.
! -----------------------------------------------------------------


! Perform work on PE=0 to avoid extra calls to random_number which
! affects bit reproducibility by changing the random seed according
! to the number of PEs
 IF  (firsttimestep == 1) THEN

  ! In this version of SKEB2, alpha is not wavenumber dependent
   alpha= 1.- EXP(-timestep/tau)

  ! Set maximum backscatter amount
   max_tot_backscat = 1.e3 * tot_backscat

  ! Set backscatter level requiring damping of forcing pattern
   crit_tot_backscat = 1.e4 * tot_backscat

   CALL backscatter_spectrum( gspect, nlim, n1, n2, tau, timestep,      &
                              tot_backscat  ,alpha)

   ! CRUN restart (must read in dpsidtc/s from seed file - if available)
   ! Fallback is to remain with initialised zero's
   IF (timestep_number > 1 .AND. l_stphseed_read) THEN
     IF (mype == 0) THEN
       CALL stph_openinput(.true.)
   
       ! Get random seed size and reset random seed to this value
       CALL random_seed(SIZE=tam)
       IF (.NOT.ALLOCATED(iranseed)) ALLOCATE(iranseed(tam))
       CALL stph_readentry2(iranseed, tam)
       CALL random_seed(PUT=iranseed(1:tam))
       DEALLOCATE(iranseed)

       CALL stph_readentry2(dpsidtc, SIZE(dpsidtc))

       CALL stph_readentry2(dpsidts, SIZE(dpsidts))

       CALL stph_closeinput()

     END IF

     ! Scatter restart coeffs to all processors
     sndcount = (n2+1)*(n2-n1+1)
     CALL gc_rbcast(3247, sndcount, zero, nproc, icode, dpsidtc)
     CALL gc_rbcast(3248, sndcount, zero, nproc, icode, dpsidts)

   END IF       ! CRUN with timestep_number > 1
 END IF         ! firsttimestep

 CALL update_dpsidt( alpha, dpsidtc, dpsidts, n1, n2, nlim, gspect,     &
                     zero, icode, info, firsttimestep)

! Compute decomposition of latitude loop
! Array "displs" keeps track of the real latitudes offset on each PE
 displs(0) = 0
 DO i = 0, nproc-1
   fft_rows(i) = nlat/nproc
   IF (i >= nproc - MOD(nlat,nproc)) Then
     fft_rows(i) = fft_rows(i) + 1
   END IF
   rowcounts(i) = fft_rows(i)
   IF (i > 0) THEN
     displs(i) = SUM(rowcounts(0:i-1))
   END IF
 END DO

! Number of rows on this PE
 my_rows = fft_rows(mype)
! Number of rows before this on lower PEs
 my_displs = displs(mype)

! Allocate local data arrays
! Real space, Fourier Sin, Cos, Radius Coeffs
 ALLOCATE (my_dpsidt(1:2*nlim, my_rows))
 ALLOCATE (my_dpsidtc(0:n2,n1:n2))
 ALLOCATE (my_dpsidts(0:n2,n1:n2))
 ALLOCATE (my_dpsidtr(0:n2,n1:n2))
 ALLOCATE (my_phi(0:n2,n1:n2))
 ALLOCATE (my_phishft(0:n2,n1:n2))

!  These arrays have sections that may not be properly initialised
 DO n = n1, n2
   DO m = n+1, n2
     my_dpsidtc(m, n) = 0.0
     my_dpsidts(m, n) = 0.0
   END DO
 END DO
 DO n = n1, n2-1
   DO m = 0, n2
     my_dpsidtr(m, n) = 0.0
     my_phishft(m, n) = 0.0
     my_phi(m, n) = 0.0
   END DO
 END DO

 DO n = n1, n2
   DO m = 0, n
     my_dpsidtr(m, n) = SQRT(dpsidtc(m, n)**2 + dpsidts(m, n)**2)
!    Determine angle from SIN & COS wave components (single step)
     my_phi(m, n) = ATAN2(dpsidts(m, n), dpsidtc(m, n))
!    Max shift ranges from 0 <-> pi  for wavenumbers n1 <-> n2
     my_phishft(m, n) = (n2 - MAX(m, n)) * pi/ (n2 - n1)
   END DO
 END DO

 DO k = skeb2_botlev, skeb2_toplev
   ! Adjust coefficients in the vertical
   ! Level 1 = no change -> 12km Level = max change (=pi)
   !  cycles around above that level
   kr = eta_theta_levels(k) * (r_theta_levels(1,1,tdims%k_end) -        &
       Earth_radius)/ 12000.
   DO n = n1, n2
     DO m = 0, n

       my_dpsidtc(m, n) = my_dpsidtr(m, n) * COS(my_phi(m, n) +         &
                          kr * my_phishft(m, n))
       my_dpsidts(m, n) = my_dpsidtr(m, n) * SIN(my_phi(m, n) +         &
                          kr * my_phishft(m, n))

     END DO
   END DO

   ii=1
   DO ilat = 1, my_rows
     mu= 2.0*(ilat-nlat/2+my_displs)/nlat

     ! Calculates dpsidt in grid space
     CALL skeb_forcing( my_dpsidtc, my_dpsidts, n1, n2, nlim, dpsidt,   &
                        mu, nlat, firsttimestep, ii)

     ! Copy dpsidt to 2-D array
     my_dpsidt(:,ilat)= dpsidt(0:2*nlim-1)
   END DO

   ! Gather field split by latitude and then distribute as normal
   CALL gc_get_communicator(icode, info)

   CALL mpl_gatherv( my_dpsidt(1,1), my_rows*2*nlim, mpl_real,          &
                     psi, rowcounts*2*nlim, displs*2*nlim, mpl_real,    &
                     zero, icode, info)

   ! Scatter global field psi to local array psif on each processor
   cmessage=''
! DEPENDS ON: scatter_field
   CALL scatter_field( psif(pdims%i_start, pdims%j_start, k), psi(1,1), &
                       row_length, rows, global_row_length, global_rows,&
                       fld_type_p, halo_type_no_halo, zero,             &
                       gc_all_proc_group, icode, cmessage)

 END DO

 ! Deallocate work arrays
 IF (ALLOCATED(my_dpsidt)) DEALLOCATE(my_dpsidt)
 IF (ALLOCATED(my_dpsidtc)) DEALLOCATE(my_dpsidtc)
 IF (ALLOCATED(my_dpsidts)) DEALLOCATE(my_dpsidts)
 IF (ALLOCATED(my_dpsidtr)) DEALLOCATE(my_dpsidtr)
 IF (ALLOCATED(my_phi)) DEALLOCATE(my_phi)
 IF (ALLOCATED(my_phishft)) DEALLOCATE(my_phishft)

 IF (l_endgame) THEN
! Reset r_u and r_v in EG configuration
   DO k = skeb2_botlev, skeb2_toplev
     DO j = udims%j_start, udims%j_end
       DO i = udims%i_start, udims%i_end
         r_u(i,j,k)     = 0.0
       END DO
     END DO
     DO j = vdims%j_start, vdims%j_end
       DO i = vdims%i_start, vdims%i_end
         r_v(i,j,k)     = 0.0
       END DO
     END DO
   END DO
 END IF

! Set outside part of arrays to zero for StdOut and STASH consistency
 DO k = pdims%k_start, skeb2_botlev - 1
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       psif(i,j,k)       = 0.0
       m_psif(i,j,k)     = 0.0
     END DO
   END DO
 END DO
 DO k = udims%k_start, skeb2_botlev - 1
   DO j = udims%j_start, udims%j_end
     DO i = udims%i_start, udims%i_end
       skeb2_urot(i,j,k) = 0.0
       skeb2_udiv(i,j,k) = 0.0
     END DO
   END DO
 END DO
 DO k = vdims%k_start, skeb2_botlev - 1
   DO j = vdims%j_start, vdims%j_end
     DO i = vdims%i_start, vdims%i_end
       skeb2_vrot(i,j,k) = 0.0
       skeb2_vdiv(i,j,k) = 0.0
     END DO
   END DO
 END DO
 DO k = skeb2_toplev + 1, pdims%k_end
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       psif(i,j,k)       = 0.0
       m_psif(i,j,k)     = 0.0
     END DO
   END DO
 END DO
 DO k = skeb2_toplev + 1, udims%k_end
   DO j = udims%j_start, udims%j_end
     DO i = udims%i_start, udims%i_end
       skeb2_urot(i,j,k) = 0.0
       skeb2_udiv(i,j,k) = 0.0
     END DO
   END DO
 END DO
 DO k = skeb2_toplev + 1, vdims%k_end
   DO j = vdims%j_start, vdims%j_end
     DO i = vdims%i_start, vdims%i_end
       skeb2_vrot(i,j,k) = 0.0
       skeb2_vdiv(i,j,k) = 0.0
     END DO
   END DO
 END DO

 IF  (printstatus  ==  prstatus_diag) THEN
   WRITE(6,*)' ---------------------------------------- '
   WRITE(6,*) 'max(dpsidtc)= ',MAXVAL(dpsidtc)
   WRITE(6,*) 'min(dpsidtc)= ',MINVAL(dpsidtc)
   WRITE(6,*) 'max(dpsidts)= ',MAXVAL(dpsidts)
   WRITE(6,*) 'min(dpsidts)= ',MINVAL(dpsidts)
   WRITE(6,*) 'max(psif)= ',MAXVAL(psif)
   WRITE(6,*) 'min(psif)= ',MINVAL(psif)
   WRITE(6,*)' ---------------------------------------- '
 END IF


! Calculate Density, mass and dZ for each level (k denotes top of level)
! These values are used in vertical integrals of Convective Dissipation
! and Kinetic Energy diagnostics at the end of this subroutine
 DO j = pdims%j_start, pdims%j_end
   DO i = pdims%i_start, pdims%i_end
     mass_2d(i,j) = 0.
   END DO
 END DO
 DO k = skeb2_botlev, skeb2_toplev
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       delr_rho(i,j,k) = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
       rho(i,j,k) = rho_r2(i,j,k)/(r_rho_levels(i,j,k) *                &
                    r_rho_levels(i,j,k))
       mass(i,j,k) = darea(i,j) * rho(i,j,k) * delr_rho(i,j,k)
       mass_2d(i,j) = mass_2d(i,j) + mass(i,j,k)
     END DO
   END DO
 END DO
! Zero mass at poles reset to unity (to avoid divide-by-zero)
 WHERE(mass_2d < 1.e3) mass_2d = 1.

! Set outside part of arrays to zero
 DO k = pdims%k_start, skeb2_botlev - 1
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       rho(i,j,k)       = 0.0
       delr_rho(i,j,k)     = 0.0
     END DO
   END DO
 END DO
 DO k = skeb2_toplev + 1, pdims%k_end
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       rho(i,j,k)       = 0.0
       delr_rho(i,j,k)     = 0.0
     END DO
   END DO
 END DO
 IF (printstatus  ==  prstatus_diag) THEN
   WRITE(6,'("***** CAPE")')
   WRITE(6,'("***** UNITS: m^2")')
   WRITE(6,'("max(skeb2_cape)=",ES22.15,12X,                            &
           & "min(skeb2_cape)=",ES22.15)')                              &
              MAXVAL(skeb2_cape),MINVAL(skeb2_cape)
   WRITE(6,'("maxloc(skeb2_cape)=",3I7,10X,                             &
           & "minloc(skeb2_cape)=",3I7)')                               &
              MAXLOC(skeb2_cape),MINLOC(skeb2_cape)
   WRITE(6,'("mean(skeb2_cape)= ",ES22.15)')                            &
              SUM(skeb2_cape)/SIZE(skeb2_cape)
   WRITE(6,'("***********************************")')
   WRITE(6,'("**SKEB2** up_flux")')
   WRITE(6,'("***** UNITS: m^2")')
   WRITE(6,'("max(skeb2_up_flux)=",ES22.15,12X,                         &
           & "min(skeb2_up_flux)=",ES22.15)')                           &
         MAXVAL(skeb2_up_flux),MINVAL(skeb2_up_flux)
   WRITE(6,'("maxloc(skeb2_up_flux)=",3I7,10X,                          &
           & "minloc(skeb2_up_flux)=",3I7)')                            &
         MAXLOC(skeb2_up_flux),MINLOC(skeb2_up_flux)
   WRITE(6,'("mean(skeb2_up_flux)= ",ES22.15)')                         &
              SUM(skeb2_up_flux)/SIZE(skeb2_up_flux)
   WRITE(6,'("***********************************")')
   WRITE(6,'("**SKEB2** dwn_flux")')
   WRITE(6,'("***** UNITS: m^2")')
   WRITE(6,'("max(skeb2_dwn_flux)=",ES22.15,12X,                        &
           & "min(skeb2_dwn_flux)=",ES22.15)')                          &
       MAXVAL(skeb2_dwn_flux),MINVAL(skeb2_dwn_flux)
   WRITE(6,'("maxloc(skeb2_dwn_flux)=",3I7,10X,                         &
           & "minloc(skeb2_dwn_flux)=",3I7)')                           &
       MAXLOC(skeb2_dwn_flux),MINLOC(skeb2_dwn_flux)
   WRITE(6,'("mean(skeb2_dwn_flux)= ",ES22.15)')                        &
            SUM(skeb2_dwn_flux)/SIZE(skeb2_dwn_flux)
   WRITE(6,'("***********************************")')
   WRITE(6,*)
   WRITE(6,'("***** darea")')
   WRITE(6,'("***** UNITS: m^2")')
   WRITE(6,'("max(darea)=",ES22.15,12X,"min(darea)=",ES22.15)')         &
              MAXVAL(darea),MINVAL(darea)
   WRITE(6,'("maxloc(darea)=",2I7,17X,"minloc(darea)=",2I7)')           &
              MAXLOC(darea),MINLOC(darea)
   WRITE(6,'("mean(darea)= ",ES22.15)') SUM(darea)/SIZE(darea)
   WRITE(6,'("***********************************")')
   WRITE(6,*)
   WRITE(6,'("***** rho")')
   WRITE(6,'("***** UNITS: kg/m^3")')
   WRITE(6,'("max(rho)=",ES22.15,12X,"min(rho)=",ES22.15)')             &
              MAXVAL(rho),MINVAL(rho)
   WRITE(6,'("maxloc(rho)=",3I7,10X,"minloc(rho)=",3I7)')               &
              MAXLOC(rho),MINLOC(rho)
   WRITE(6,'("mean(rho)= ",ES22.15)') SUM(rho)/SIZE(rho)
   WRITE(6,'("***********************************")')
   WRITE(6,*)
   WRITE(6,'("***** dZ")')
   WRITE(6,'("***** UNITS: m")')
   WRITE(6,'("max(delr_rho)=",ES22.15,12X,"min(delr_rho)=",ES22.15)')   &
              MAXVAL(delr_rho),MINVAL(delr_rho)
   WRITE(6,'("maxloc(delr_rho)=",3I7,10X,"minloc(delr_rho)=",3I7)')     &
              MAXLOC(delr_rho),MINLOC(delr_rho)
   WRITE(6,'("mean(delr_rho)=",ES22.15)') SUM(delr_rho)/SIZE(delr_rho)
   WRITE(6,'("***********************************")')
   WRITE(6,*)
   WRITE(6,'("***** 2D Mass")')
   WRITE(6,'("***** UNITS: kg")')
   WRITE(6,'("max(mass_2d)=",ES22.15,12X,"min(mass_2d)=",ES22.15)')     &
              MAXVAL(mass_2d),MINVAL(mass_2d)
   WRITE(6,'("maxloc(mass_2d)=",2I7,17X,"minloc(mass_2d)=",2I7)')       &
              MAXLOC(mass_2d),MINLOC(mass_2d)
   WRITE(6,'("mean(mass_2d)= ",ES22.15)') SUM(mass_2d)/SIZE(mass_2d)
   WRITE(6,'("***********************************")')
   WRITE(6,*)
 END IF

! Initialise all dissipation fields. They are conditionally filled based
! on logicals l_skeb2_psisdisp, l_skeb2_psicdisp and l_skeb2_skeb1disp
 DO k = pdims%k_start, pdims%k_end
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       sdisp(i,j,k) = 0.0
       cdisp(i,j,k) = 0.0
       kdisp(i,j,k) = 0.0
     END DO
   END DO
 END DO


! --------------------------------------------------------------------
! Calculate energy dissipated from advection (using Smagorinsky
!  diffusion)
! --------------------------------------------------------------------

 IF (l_skeb2_psisdisp) THEN
   CALL skeb_smagorinsky(   rows, row_length, n_rows                    &
 ,                          model_levels                                &
 ,                          u, v                                        &
 ,                          delta_lambda, delta_phi                     &
 ,                          sdisp                                       &
 ,                          firsttimestep                               &
                          )

   IF  (printstatus  >  prstatus_normal) THEN
     WRITE(6,'("***** 3D Energy Dissipation Rate per volume in Smag.")')
     WRITE(6,'("***** UNITS: m^2.s^-3 or kg.m^2.s^-2/s/kg")')
     WRITE(6,'("max(sdisp)=",ES22.15,12X,"min(sdisp)=",ES22.15)')       &
                MAXVAL(sdisp),MINVAL(sdisp)
     WRITE(6,'("maxloc(sdisp)=",3I7,10X,"minloc(sdisp)=",3I7)')         &
                MAXLOC(sdisp),MINLOC(sdisp)
     WRITE(6,'("mean(sdisp)= ",ES22.15)') SUM(sdisp)/SIZE(sdisp)
     WRITE(6,'("***********************************")')
     WRITE(6,*)
   END IF

 END IF     ! l_skeb2_psisdisp


! --------------------------------------------------------------------
!  Calculate energy dissipated by convection (CAPE*mass_flux)
!  to re-scale the amount of energy feedback which is
!  10E-4 (see variable: tot_backscat) from the sph.harm. field
! --------------------------------------------------------------------

 IF (l_skeb2_psicdisp) THEN
   IF (skeb2_cdisp .EQ. type_cape) THEN
     ! ---------------------------------------------------------------
     ! Calculate modulating field (based on CAPE and with a vertical
     ! profile derived from dM/dz)
     ! ---------------------------------------------------------------

     ! Calculate diff flux
     DO k = skeb2_botlev - 1, skeb2_toplev + 1
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
           diff_flux(i,j,k) = skeb2_up_flux(i,j,k) -                    &
                              skeb2_dwn_flux(i,j,k)
         END DO
       END DO
     END DO
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         cape_kk(i,j) = skeb2_cape(i,j)
       END DO
     END DO
     ! Remove negative/noise values in CAPE array
     WHERE (cape_kk < 0.01) cape_kk = 0.0

     ! ---------------------------------------------------------------
     ! Calculate sm_tot_disp =  vertical profile from centred derivative
     !                     of dIF f_flux
     !                     Horizontal (xy) field from cape
     ! --------------------------------------------------------------------

     DO k = skeb2_botlev, skeb2_toplev
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

     !   Need to divide by (g*rho*dz) to convert Pa.s^-1 to s^-1
     !   Then: multiply by CAPE (J.kg^-1) => J.kg^-1.s^-1  or  m^2.s^-3
     !   Remember: delr_rho is a model-level thickness and diff_flux
     !             is differenced over 2 levels, => delr_rho@(k)+(k+1)
           cdisp(i,j,k)=(ABS(diff_flux(i,j,k+1)-diff_flux(i,j,k-1))     &
                        /((delr_rho(i,j,k)+delr_rho(i,j,k+1))           &
                        *g*rho(i,j,k)))*cape_kk(i,j)
         END DO
       END DO
     END DO

     IF  (printstatus  >  prstatus_normal) THEN
       WRITE(6,'("*** 3D Energy Diss Rate per vol in Conv M-FLX*CAPE")')
       WRITE(6,'("*** UNITS: m^2.s^-3 or kg.m^2.s^-2/s/kg")')
       WRITE(6,'("max(cdisp)=",ES22.15,12X,"min(cdisp)=",ES22.15)')     &
                  MAXVAL(cdisp),MINVAL(cdisp)
       WRITE(6,'("maxloc(cdisp)=",3I7,10X,"minloc(cdisp)=",3I7)')       &
                  MAXLOC(cdisp),MINLOC(cdisp)
       WRITE(6,'("mean(cdisp)= ",ES22.15)') SUM(cdisp)/SIZE(cdisp)
       WRITE(6,'("***********************************")')
       WRITE(6,*)
     END IF

!!!   END IF  ! (l_skeb2_cdisp_cape)

   ELSE IF (skeb2_cdisp .EQ. type_mflx) THEN
     DO k = skeb2_botlev,skeb2_toplev
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

     !   KE= 0.5 * m * v^2
     !     = 1/g*rho * d(upflux)/dz * (upflux/g*rho*alphac)**2
           cdisp(i,j,k)=(MAX(skeb2_up_flux(i,j,k+1) -                   &
                             skeb2_up_flux(i,j,k-1),0.)/                &
                     (g*rho(i,j,k)*(delr_rho(i,j,k)+delr_rho(i,j,k+1))))&
                        *(MAX(skeb2_up_flux(i,j,k),0.)/                 &
                        (g*rho(i,j,k)*alphac))**2
         END DO
       END DO
     END DO

     IF  (printstatus  >  prstatus_normal) THEN
       WRITE(6,'("*** 3D Energy Diss Rate per vol in Conv M-FLX=>w")')
       WRITE(6,'("*** UNITS: m^2.s^-3 or kg.m^2.s^-2/s/kg")')
       WRITE(6,'("max(cdisp)=",ES22.15,12X,"min(cdisp)=",ES22.15)')     &
               MAXVAL(cdisp),MINVAL(cdisp)
       WRITE(6,'("maxloc(cdisp)=",3I7,10X,"minloc(cdisp)=",3I7)')       &
               MAXLOC(cdisp),MINLOC(cdisp)
       WRITE(6,'("mean(cdisp)= ",ES22.15)') SUM(cdisp)/SIZE(cdisp)
       WRITE(6,'("***********************************")')
       WRITE(6,*)
     END IF

   END IF  ! (skeb2_cdisp)
 END IF  ! (l_skeb2_psicdisp)


! -----------------------------------------------------------
! Calculate assumed energy dissipation (using SKEB1-type KE)
! This is similar to Smagorinsky turbulence, but focuses on
!   the Kinetic Energy field.
! -----------------------------------------------------------

 IF (l_skeb2_skeb1disp) THEN
     CALL stph_skeb1( rows, row_length, n_rows                          &
 ,                          model_levels                                &
 ,                          mass, c_dadt                                &
 ,                          u, v                                        &
 ,                          kdisp                                       &
                          )

   IF  (printstatus  >  prstatus_normal) THEN
     WRITE(6,'("***** 3D Energy Dissipation Rate per volume in Smag.")')
     WRITE(6,'("***** UNITS: m^2.s^-3 or kg.m^2.s^-2/s/kg")')
     WRITE(6,'("max(kdisp)=",ES22.15,12X,"min(kdisp)=",ES22.15)')       &
                MAXVAL(kdisp),MINVAL(kdisp)
     WRITE(6,'("maxloc(kdisp)=",3I7,10X,"minloc(kdisp)=",3I7)')         &
                MAXLOC(kdisp),MINLOC(kdisp)
     WRITE(6,'("mean(kdisp)= ",ES22.15)') SUM(kdisp)/SIZE(kdisp)
     WRITE(6,'("***********************************")')
     WRITE(6,*)
   END IF

 END IF  ! (l_skeb2_skeb1disp)


! Calculate total dissipation field
 DO k = skeb2_botlev, skeb2_toplev
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       sm_tot_disp(i,j,k) = cdispfac*cdisp(i,j,k) +                     &
                   sdispfac*sdisp(i,j,k) + kdispfac*kdisp(i,j,k)
     END DO
   END DO
 END DO

 IF  (printstatus  ==  prstatus_diag) THEN
    WRITE(6,'("***** Total Dissipation Field *************")')
    WRITE(6,'("Convection Diss. Factor (cdispfac) = ",F6.3)') cdispfac
    WRITE(6,'(" Numerical Diss. Factor (sdispfac) = ",F6.3)') sdispfac
    WRITE(6,'("  SKEB1-KE Diss. Factor (kdispfac) = ",F6.3)') kdispfac
 END IF

! Iterative cycling spreads isolated convective dissipation elements
! more evenly (number of iterations controlled by nsmooth)
 IF (.NOT. l_skebsmooth_adv) THEN

! Limiter on energy backscattered
   WHERE(sm_tot_disp > max_tot_backscat)                                &
                       sm_tot_disp = max_tot_backscat

   IF  (printstatus  ==  prstatus_diag) THEN
  !  Fill perimeter of array with zero for diagnostic print consistency
     DO k = pdims%k_start, skeb2_botlev - 1
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
           sm_tot_disp(i,j,k) = 0.
         END DO
       END DO
     END DO
     DO k = skeb2_toplev + 1, pdims%k_end
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
           sm_tot_disp(i,j,k) = 0.
         END DO
       END DO
     END DO
     WRITE(6,'("max(sm_tot_disp)=",ES22.15,12x,"min(sm_tot_disp)=",     &
                &ES22.15)') MAXVAL(sm_tot_disp),MINVAL(sm_tot_disp)
     WRITE(6,'("maxloc(sm_tot_disp)=",3I7,10X,"minloc(sm_tot_disp)=",   &
                &3I7)') MAXLOC(sm_tot_disp),MINLOC(sm_tot_disp)
     WRITE(6,'("mean(sm_tot_disp)= ",ES22.15)')                         &
                 SUM(sm_tot_disp)/SIZE(sm_tot_disp)
     WRITE(6,'("***********************************")')
   END IF

  !  Allocate smoothing buffer array
   IF (.NOT. ALLOCATED(sm_tot_disp0)) THEN
     ALLOCATE(sm_tot_disp0(pdims_s%i_start:pdims_s%i_end,               &
                           pdims_s%j_start:pdims_s%j_end,               &
                           pdims_s%k_start:pdims_s%k_end))
   END IF

   DO ismooth = 1, nsmooth

     IF (printstatus  ==  prstatus_diag) THEN
       WRITE(6,*) 'starting ismooth = ', ismooth ,' of ', nsmooth
     END IF

  !  Fill buffer array with latest data and exchange halos
     DO k = skeb2_botlev,skeb2_toplev
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
           sm_tot_disp0(i,j,k)=sm_tot_disp(i,j,k)
         END DO
       END DO
       IF (.NOT. l_endgame) THEN
         ! Set dissipation to zero at pole for ND grid
         IF (at_extremity(psouth)) THEN
           DO i = pdims%i_start, pdims%i_end
             sm_tot_disp0(i,pdims%j_start,k) = 0.
           END DO
         END IF
         IF (at_extremity(pnorth)) THEN
           DO i = pdims%i_start, pdims%i_end
             sm_tot_disp0(i,pdims%j_end,k) = 0.
           END DO
         END IF
       END IF
     END DO

  ! DEPENDS ON: swap_bounds
     CALL swap_bounds( sm_tot_disp0, row_length, rows,                    &
                       skeb2_toplev, offx, offy, fld_type_p,              &
                       .FALSE.  )


  !  2D 1-2-1 spatial filter
     DO k = skeb2_botlev, skeb2_toplev
       DO j = pdims%j_start, pdims%j_end
         jm1 = j - 1
         jp1 = j + 1
         IF (at_extremity(psouth)) jm1 = MAX(jm1, pdims%j_start)
         IF (at_extremity(pnorth)) jp1 = MIN(jp1, pdims%j_end)
         DO i = pdims%i_start, pdims%i_end
           im1 = i - 1
           ip1 = i + 1
           sm_tot_disp(i,j,k) = 0.0625*(sm_tot_disp0(im1,jp1,k) +         &
                                        sm_tot_disp0(ip1,jp1,k) +         &
                                        sm_tot_disp0(im1,jm1,k) +         &
                                        sm_tot_disp0(ip1,jm1,k) +         &
                                        2*( sm_tot_disp0(i,jp1,k) +       &
                                            sm_tot_disp0(im1,j,k) +       &
                                            sm_tot_disp0(ip1,j,k) +       &
                                            sm_tot_disp0(i,jm1,k) ) +     &
                                          4*sm_tot_disp0(i,j,k) )
         END DO  ! i
       END DO  ! j
     END DO  ! k

     IF (printstatus  ==  prstatus_diag) THEN
       WRITE(6,'("***** Smoothed Dissipation Field *************")')
       WRITE(6,'("max(sm_tot_disp)=",ES22.15,12X,"min(sm_tot_disp)=",     &
                &ES22.15)') MAXVAL(sm_tot_disp),MINVAL(sm_tot_disp)
       WRITE(6,'("maxloc(sm_tot_disp)=",3I7,10X,"minloc(sm_tot_disp)=",   &
                &3I7)') MAXLOC(sm_tot_disp),MINLOC(sm_tot_disp)
       WRITE(6,'("mean(sm_tot_disp)= ",ES22.15)')                         &
                 SUM(sm_tot_disp)/SIZE(sm_tot_disp)
       WRITE(6,'("***********************************")')
     END IF

   END DO   ! ismooth iteration

   IF (.NOT. l_endgame) THEN
     ! Set dissipation to zero at pole for ND grid
     DO k = skeb2_botlev, skeb2_toplev
       IF (at_extremity(psouth)) THEN
         DO i = pdims%i_start, pdims%i_end
           sm_tot_disp(i,pdims%j_start,k) = 0.
         END DO
       END IF
       IF (at_extremity(pnorth)) THEN
         DO i = pdims%i_start, pdims%i_end
           sm_tot_disp(i,pdims%j_end,k) = 0.
         END DO
       END IF
     END DO
   END IF

   ! DEALLOCATE smoothing buffer array
   DEALLOCATE(sm_tot_disp0)

 ELSE

! Advanced smoothing option: using pre-defined smoothing array
! Removes backscatter in larger region of instability

  !  Allocate smoothing buffer array (using STPH-defined halo)
   IF (.NOT. ALLOCATED(sm_tot_disp0)) THEN
     ALLOCATE(sm_tot_disp0(stphdims_l%i_start:stphdims_l%i_end,         &
                           stphdims_l%j_start:stphdims_l%j_end,         &
                           stphdims_l%k_start:stphdims_l%k_end))
   END IF

  !  Allocate pattern smoothing array mask
   IF (.NOT. ALLOCATED(psif0)) THEN
     ALLOCATE(psif0(stphdims_l%i_start:stphdims_l%i_end,                &
                    stphdims_l%j_start:stphdims_l%j_end))
   END IF

  ! Initialise mask to zero
  ! This is a 2-D array and will be set to one if instability is
  ! detected at any level in the column
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       psif0(i,j) = 0.0
     END DO  ! i
   END DO  ! j

  !  Fill buffer array with latest data and exchange halos
  ! When excessive backscatter diagnosed, do not backscatter in this
  ! column by setting forcing pattern to zero in the halo region
  ! [offx_stph:offy_stph]
  ! ... useful diagnostic to find instabilities leading to GPStorms
   DO k = skeb2_botlev,skeb2_toplev
     levfac = MAX(0.1, LOG10(logscale*k))   ! Max factor=10
     local_crit_tot_backscat = crit_tot_backscat/levfac
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         ! Case of severe instability
         IF (sm_tot_disp(i,j,k) > local_crit_tot_backscat) THEN
           IF (l_skebprint) THEN
             WRITE(6,'("WARNING: SKEB Backscatter = ",ES12.5,'//&
                 '" at [",I4,";",I4,";",I4,"] EXCEEDS limit of ",ES12.5)')   &
                        sm_tot_disp(i,j,k),i,j,k,local_crit_tot_backscat
           END IF
           sm_tot_disp(i,j,k) = 0.0
           psif0(i,j) = 1.0
         END IF
         sm_tot_disp(i,j,k) = MIN(max_tot_backscat, sm_tot_disp(i,j,k))
         sm_tot_disp0(i,j,k) = sm_tot_disp(i,j,k)
         ! Reset to zero for iterative smoothing at next step
         sm_tot_disp(i,j,k) = 0.0
       END DO
     END DO
   END DO

   ! DEPENDS ON: swap_bounds
   CALL swap_bounds( psif0, row_length, rows, 1,                        &
                     offx_stph, offy_stph, fld_type_p, .FALSE.  )

   ! Reduce dissipation in region around instability
   ! Area of influence can cross PE boundaries, so central point of
   ! smoothing can lie in the halo region.
   DO j = stphdims_l%j_start, stphdims_l%j_end
     DO i = stphdims_l%i_start, stphdims_l%i_end

       IF (psif0(i,j) > 0.5) THEN
         ! Apply at all SKEB levels
         i1 = MAX(pdims%i_start, i-offx_stph)
         i2 = MIN(pdims%i_end,   i+offx_stph)
         j1 = MAX(pdims%j_start, j-offy_stph)
         j2 = MIN(pdims%j_end,   j+offy_stph)
         DO k = skeb2_botlev, skeb2_toplev
           DO jj = j1, j2
             DO ii = i1, i2
               sm_tot_disp0(ii,jj,k) = mask_pdamp(ii-i,jj-j) *          &
                                           sm_tot_disp0(ii,jj,k)
             END DO      ! ii
           END DO      ! jj
         END DO      ! k
       END IF      ! Unstable point (can also be in halo)

     END DO      ! i
   END DO      ! j

   ! DEPENDS ON: swap_bounds
   CALL swap_bounds( sm_tot_disp0, row_length, rows,                    &
                     skeb2_toplev, offx_stph, offy_stph, fld_type_p,    &
                     .FALSE.  )

  !  2D 1-2-1 spatial filter (applied nsmooth times)
   DO k = skeb2_botlev, skeb2_toplev
  !  Apply smoother
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         DO jj = -offy_stph, offy_stph
           DO ii = -offx_stph, offx_stph
             sm_tot_disp(i,j,k) = sm_tot_disp(i,j,k) +                  &
                mask_smooth(ii,jj) * sm_tot_disp0(ii+i,jj+j,k)
           END DO  ! ii
         END DO  ! jj
       END DO  ! i
     END DO  ! j
     IF (.NOT. l_endgame) THEN
       ! Set dissipation to zero at pole for ND grid
       IF (at_extremity(psouth)) THEN
         DO i = pdims%i_start, pdims%i_end
           sm_tot_disp(i,pdims%j_start,k) = 0.
         END DO
       END IF
       IF (at_extremity(pnorth)) THEN
         DO i = pdims%i_start, pdims%i_end
           sm_tot_disp(i,pdims%j_end,k) = 0.
         END DO
       END IF
     END IF
   END DO  ! SKEB levels

  !  Fill perimeter of array with zero for diagnostic print consistency
   DO k = pdims%k_start, skeb2_botlev - 1
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         sm_tot_disp(i,j,k) = 0.
       END DO
     END DO
   END DO
   DO k = skeb2_toplev + 1, pdims%k_end
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         sm_tot_disp(i,j,k) = 0.
       END DO
     END DO
   END DO
   IF (printstatus  ==  prstatus_diag) THEN
     WRITE(6,'("***** Smoothed Dissipation Field *************")')
     WRITE(6,'("max(sm_tot_disp)=",ES22.15,12X,"min(sm_tot_disp)=",'//  &
         'ES22.15)') MAXVAL(sm_tot_disp),MINVAL(sm_tot_disp)
     WRITE(6,'("maxloc(sm_tot_disp)=",3I7,10X,"minloc(sm_tot_disp)=",'//&
         '3I7)') MAXLOC(sm_tot_disp),MINLOC(sm_tot_disp)
     WRITE(6,'("mean(sm_tot_disp)= ",ES22.15)')                         &
               SUM(sm_tot_disp)/SIZE(sm_tot_disp)
     WRITE(6,'("***********************************")')
   END IF

   ! DEALLOCATE smoothing buffer and pattern mask arrays
   DEALLOCATE(sm_tot_disp0)
   DEALLOCATE(psif0)

 END IF   ! Smoothing option: l_skebsmooth_adv


! --------------------------------------------------------------------
!              Calculate u,v increments
! --------------------------------------------------------------------


! Modulate streamfunction forcing field with smoothed E. Diss Rate
!
 DO k = skeb2_botlev, skeb2_toplev
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       m_psif(i,j,k) = psif(i,j,k) *                                    &
                       SQRT(br*sm_tot_disp(i,j,k)/tot_backscat)
     END DO  ! i
   END DO    ! j
 END DO      ! k


! SWAP BOUNDS: used to communicate information from neighbouring
!              processors to make the smoothing correctly

 DO k = skeb2_botlev, skeb2_toplev
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       work_p_halo(i,j,k) = m_psif(i,j,k)
     END DO  ! i
   END DO    ! j
 END DO      ! k

! DEPENDS ON: swap_bounds
 CALL swap_bounds( work_p_halo, row_length, rows,                       &
                   skeb2_toplev, offx, offy, fld_type_p,                &
                   .FALSE.  )

 DO k = skeb2_botlev, skeb2_toplev
   DO j = vdims%j_start, vdims%j_end
     jp1 = j + 1
     DO i = vdims%i_start, vdims%i_end
       ip1 = i + 1
       work_z_halo(i,j,k) = 0.25*(                                      &
                        work_p_halo(i,j,k) + work_p_halo(ip1,j,k) +     &
                        work_p_halo(i,jp1,k) + work_p_halo(ip1,jp1,k))
     END DO  ! i
   END DO    ! j
   ! Set dissipation to zero on Pole row
   IF (l_endgame) THEN
     IF (at_extremity(psouth)) THEN
       DO i = vdims%i_start, vdims%i_end
         work_z_halo(i,vdims%j_start,k) = 0.
       END DO  ! i
     END IF
     IF (at_extremity(pnorth)) THEN
       DO i = vdims%i_start, vdims%i_end
         work_z_halo(i,vdims%j_end,k) = 0.
       END DO  ! i
     END IF
   END IF
 END DO      ! k

! DEPENDS ON: swap_bounds
 CALL swap_bounds( work_z_halo, row_length, n_rows,                     &
                   skeb2_toplev, offx, offy, fld_type_v,                &
                   .FALSE.  )

 DO k = skeb2_botlev, skeb2_toplev
   ! Wind Incr on U-grid
   DO j = udims%j_start, udims%j_end
     jm1 = j - 1
     jv = j
     IF (.NOT. l_endgame) THEN
       ! Avoid using work_z_halo values beyond limits of V-grid
       IF (at_extremity(psouth)) jm1 = MAX( jm1, vdims%j_start)
       IF (at_extremity(pnorth)) jv = MIN( jv, vdims%j_end)
     END IF
     DO i = udims%i_start, udims%i_end
       ! -------------------------------------------------------------
       ! These are the ROTATIONAL incr :: from m_psif on vort points
       ! -------------------------------------------------------------
       skeb2_urot0(i,j,k)=(work_z_halo(i,jm1,k)-work_z_halo(i,jv,k))    &
                           /dy_theta(i,j)
     END DO    ! I
     IF (l_skeb2_velpot) THEN
       DO i = udims%i_start, udims%i_end
         ip1 = i + 1
         ! --------------------------------------------------
         ! These are the DIVERGENT increments
         !  Careful: +ve Vel Pot = conv in both hemispheres
         !           +ve StreamF = clockwise in both H
         !           Need (div)convergence in (anti)cyclone
         !           # sin(lat) < 0 in S.H. => (-1)*div-comp
         !           # using single column from sin_theta_latitude
         ! In future releases the Velocity Potential field  
         !   will be calculated independently               
         ! --------------------------------------------------
         skeb2_udiv0(i,j,k)=((work_p_halo(i,j,k)-work_p_halo(ip1,j,k))* &
                  SIGN(1.,sin_theta_latitude(1,j)))/ dx_theta(i,j)

       END DO    ! I
     END IF
   END DO     ! J
   ! Wind Incr on V-grid
   DO j = vdims%j_start, vdims%j_end
     DO i = vdims%i_start, vdims%i_end
       im1 = i - 1
       skeb2_vrot0(i,j,k)=(work_z_halo(i,j,k)-work_z_halo(im1,j,k))     &
                           /dx_v(i,j)

     END DO    ! I
     IF (l_skeb2_velpot) THEN
       jp1 = j + 1
       jv = j
       IF (l_endgame) THEN
         ! Avoid using work_p_halo values beyond limits of P-grid at pole
         IF (at_extremity(pnorth)) jp1 = MIN( jp1, pdims%j_end)
         IF (at_extremity(psouth)) jv = MAX( jv, pdims%j_start)
       END IF
       DO i = vdims%i_start, vdims%i_end
         skeb2_vdiv0(i,j,k)=((work_p_halo(i,jv,k)-work_p_halo(i,jp1,k))*&
                  SIGN(1.,sin_v_latitude(i,jv)))/ dy_v(i,jv)
       END DO    ! I
     END IF
   END DO     ! J
 END DO      ! K

! Deallocate work arrays
 IF (ALLOCATED(work_p_halo)) DEALLOCATE(work_p_halo)
 IF (ALLOCATED(work_z_halo)) DEALLOCATE(work_z_halo)

! Use Multi-variable swap-bounds to exploit MPP
 ifield = 0
 ifield = ifield + 1
 fields_to_swap(ifield) % field       => skeb2_urot0(:,:,:)
 fields_to_swap(ifield) % field_type  =  fld_type_u
 fields_to_swap(ifield) % levels      =  model_levels
 fields_to_swap(ifield) % rows        =  rows
 fields_to_swap(ifield) % vector      =  .TRUE.

 ifield = ifield + 1
 fields_to_swap(ifield) % field       => skeb2_vrot0(:,:,:)
 fields_to_swap(ifield) % field_type  =  fld_type_v
 fields_to_swap(ifield) % levels      =  model_levels
 fields_to_swap(ifield) % rows        =  n_rows
 fields_to_swap(ifield) % vector      =  .TRUE.

 IF (l_skeb2_velpot) THEN
   ifield = ifield + 1
   fields_to_swap(ifield) % field       => skeb2_udiv0(:,:,:)
   fields_to_swap(ifield) % field_type  =  fld_type_u
   fields_to_swap(ifield) % levels      =  model_levels
   fields_to_swap(ifield) % rows        =  rows
   fields_to_swap(ifield) % vector      =  .TRUE.

   ifield = ifield + 1
   fields_to_swap(ifield) % field       => skeb2_vdiv0(:,:,:)
   fields_to_swap(ifield) % field_type  =  fld_type_v
   fields_to_swap(ifield) % levels      =  model_levels
   fields_to_swap(ifield) % rows        =  n_rows
   fields_to_swap(ifield) % vector      =  .TRUE.
 END IF

! DEPENDS ON: swap_bounds_mv
 CALL swap_bounds_mv( fields_to_swap, ifield, row_length, offx, offy )


! --------------------------------------------------------------------
!  Smoothing of the wind fields to reduce noise (1-2-1 filter)
! --------------------------------------------------------------------

 DO k = skeb2_botlev, skeb2_toplev
   ! u-component
   DO j = udims%j_start, udims%j_end
     jm1= j - 1
     jp1= j + 1
     IF (at_extremity(psouth)) jm1 = MAX(jm1, udims%j_start)
     IF (at_extremity(pnorth)) jp1 = MIN(jp1, udims%j_end)
     DO i = udims%i_start, udims%i_end
       im1= i - 1
       ip1= i + 1
       ! -----------------------
       ! Rotational part
       ! -----------------------
       skeb2_urot(i,j,k)=0.0625*(skeb2_urot0(im1,jp1,k)+                &
                                   skeb2_urot0(ip1,jp1,k)+              &
                                   skeb2_urot0(im1,jm1,k)+              &
                                   skeb2_urot0(ip1,jm1,k)+              &
                                   2*( skeb2_urot0(i,jp1,k)+            &
                                       skeb2_urot0(im1,j,k)+            &
                                       skeb2_urot0(ip1,j,k)+            &
                                       skeb2_urot0(i,jm1,k) )+          &
                                   4*skeb2_urot0(i,j,k) )*timestep
     END DO !I
     IF (l_skeb2_velpot) THEN
       DO i = udims%i_start, udims%i_end
         im1= i - 1
         ip1= i + 1
         ! -----------------------
         ! Divergent part
         ! -----------------------
         skeb2_udiv(i,j,k)=0.0625*(skeb2_udiv0(im1,jp1,k)+              &
                                       skeb2_udiv0(ip1,jp1,k)+          &
                                       skeb2_udiv0(im1,jm1,k)+          &
                                       skeb2_udiv0(ip1,jm1,k)+          &
                                       2*( skeb2_udiv0(i,jp1,k)+        &
                                           skeb2_udiv0(im1,j,k)+        &
                                           skeb2_udiv0(ip1,j,k)+        &
                                           skeb2_udiv0(i,jm1,k) )+      &
                                       4*skeb2_udiv0(i,j,k) )*timestep
       END DO !I
     ELSE
       ! This array is not initialised otherwise
       DO i = udims%i_start, udims%i_end
         skeb2_udiv(i,j,k)=0.0
       END DO !I
     END IF
   END DO !J
   ! v-component !
   DO j = vdims%j_start, vdims%j_end
     jm1= j - 1
     jp1= j + 1
     IF (at_extremity(psouth)) jm1 = MAX(jm1, vdims%j_start)
     IF (at_extremity(pnorth)) jp1 = MIN(jp1, vdims%j_end)
     DO i = vdims%i_start, vdims%i_end
       im1= i - 1
       ip1= i + 1
       ! -----------------------
       ! Rotational part
       ! -----------------------
       skeb2_vrot(i,j,k)=0.0625*(skeb2_vrot0(im1,jp1,k)+                &
                                   skeb2_vrot0(ip1,jp1,k)+              &
                                   skeb2_vrot0(im1,jm1,k)+              &
                                   skeb2_vrot0(ip1,jm1,k)+              &
                                   2*( skeb2_vrot0(i,jp1,k)+            &
                                       skeb2_vrot0(im1,j,k)+            &
                                       skeb2_vrot0(ip1,j,k)+            &
                                       skeb2_vrot0(i,jm1,k) )+          &
                                  4*skeb2_vrot0(i,j,k) )*timestep
     END DO !I
     IF (l_skeb2_velpot) THEN
       DO i = vdims%i_start, vdims%i_end
         im1= i - 1
         ip1= i + 1
         ! -----------------------
         ! Divergent part
         ! -----------------------
         skeb2_vdiv(i,j,k)=0.0625*(skeb2_vdiv0(im1,jp1,k)+              &
                                     skeb2_vdiv0(ip1,jp1,k)+            &
                                     skeb2_vdiv0(im1,jm1,k)+            &
                                     skeb2_vdiv0(ip1,jm1,k)+            &
                                     2*( skeb2_vdiv0(i,jp1,k)+          &
                                         skeb2_vdiv0(im1,j,k)+          &
                                         skeb2_vdiv0(ip1,j,k)+          &
                                         skeb2_vdiv0(i,jm1,k) )+        &
                                     4*skeb2_vdiv0(i,j,k) )*timestep
       END DO
     ELSE
       ! This array is not initialised otherwise
       DO i = vdims%i_start, vdims%i_end
         skeb2_vdiv(i,j,k)=0.0
       END DO
     END IF
   END DO
   IF (l_endgame) THEN
     ! Set value of v-wind incr at pole = zero
     IF (at_extremity(psouth)) THEN
       DO i = vdims%i_start, vdims%i_end
         skeb2_vrot(i,vdims%j_start,k)=0.
         skeb2_vdiv(i,vdims%j_start,k)=0.
       END DO
     END IF
     IF (at_extremity(pnorth)) THEN
       DO i = vdims%i_start, vdims%i_end
         skeb2_vrot(i,vdims%j_end,k)=0.
         skeb2_vdiv(i,vdims%j_end,k)=0.
       END DO
     END IF
   ELSE
     ! Set value of u-wind incr at pole = zero
     IF (at_extremity(psouth)) THEN
       DO i = udims%i_start, udims%i_end
         skeb2_urot(i,udims%j_start,k)=0.
         skeb2_udiv(i,udims%j_start,k)=0.
       END DO
     END IF
     IF (at_extremity(pnorth)) THEN
       DO i = udims%i_start, udims%i_end
         skeb2_urot(i,udims%j_end,k)=0.
         skeb2_udiv(i,udims%j_end,k)=0.
       END DO
     END IF
   END IF
 END DO !k


! --------------------------------------------------------------------
!  Reduce wind increments to 0 at level=1 from level~2km using LOG10
!  This avoids creating large values near the top of the boundary lyr
! --------------------------------------------------------------------
 DO k = skeb2_botlev, level2km
   levfac = MAX(0.,LOG10(logscale*k))   ! log decrease to zero
   IF (printstatus  >  prstatus_oper) THEN
     WRITE(6,'(" Factor at level (",I4,") = ",ES12.4)') k, levfac
   END IF
   DO j = udims%j_start, udims%j_end
     DO i = udims%i_start, udims%i_end
       skeb2_urot(i,j,k) = skeb2_urot(i,j,k) * levfac
       skeb2_udiv(i,j,k) = skeb2_udiv(i,j,k) * levfac
     END DO
   END DO
   DO j = vdims%j_start, vdims%j_end
     DO i = vdims%i_start, vdims%i_end
       skeb2_vrot(i,j,k) = skeb2_vrot(i,j,k) * levfac
       skeb2_vdiv(i,j,k) = skeb2_vdiv(i,j,k) * levfac
     END DO
   END DO
 END DO
  
! Ramp off backscatter (linear) at top 3 levels of SKEB2 range
 IF (l_skebsmooth_adv) THEN
   DO k = skeb2_toplev-2, skeb2_toplev
     levfac = 0.25 * (skeb2_toplev - k + 1)
     DO j = udims%j_start, udims%j_end
       DO i = udims%i_start, udims%i_end
         skeb2_urot(i,j,k) = skeb2_urot(i,j,k) * levfac
         skeb2_udiv(i,j,k) = skeb2_udiv(i,j,k) * levfac
       END DO
     END DO
     DO j = vdims%j_start, vdims%j_end
       DO i = vdims%i_start, vdims%i_end
         skeb2_vrot(i,j,k) = skeb2_vrot(i,j,k) * levfac
         skeb2_vdiv(i,j,k) = skeb2_vdiv(i,j,k) * levfac
       END DO
     END DO
   END DO
 END IF   ! Advanced smoothing option


 IF (printstatus  ==  prstatus_diag) THEN
   WRITE(6,*)' ---------------------------------------- '
   WRITE(6,'("max(R_u)=",ES22.15,12X,"min(R_u)=",ES22.15)')             &
              MAXVAL(r_u),MINVAL(r_u)
   WRITE(6,'("maxloc(R_u)=",3I7,10X,"minloc(R_u)=",3I7)')               &
              MAXLOC(r_u),MINLOC(r_u)
   WRITE(6,'("mean(R_u)= ",ES22.15)') SUM(ABS(r_u))/SIZE(r_u)

   WRITE(6,'("max(R_v)=",ES22.15,12X,"min(R_v)=",ES22.15)')             &
              MAXVAL(r_v),MINVAL(r_v)
   WRITE(6,'("maxloc(R_v)=",3I7,10X,"minloc(R_v)=",3I7)')               &
              MAXLOC(r_v),MINLOC(r_v)
   WRITE(6,'("mean(R_v)= ",ES22.15)') SUM(ABS(r_v))/SIZE(r_v)

   WRITE(6,'("max(urot)=",ES22.15,12X,"min(urot)=",ES22.15)')           &
              MAXVAL(skeb2_urot),MINVAL(skeb2_urot)
   WRITE(6,'("maxloc(urot)=",3I7,10X,"minloc(urot)=",3I7)')             &
              MAXLOC(skeb2_urot),MINLOC(skeb2_urot)
   WRITE(6,'("mean(urot)= ",ES22.15)') SUM(ABS(skeb2_urot))             &
              /SIZE(skeb2_urot)

   WRITE(6,'("max(vrot)=",ES22.15,12X,"min(vrot)=",ES22.15)')           &
              MAXVAL(skeb2_vrot),MINVAL(skeb2_vrot)
   WRITE(6,'("maxloc(vrot)=",3I7,10X,"minloc(vrot)=",3I7)')             &
              MAXLOC(skeb2_vrot),MINLOC(skeb2_vrot)
   WRITE(6,'("mean(vrot)= ",ES22.15)') SUM(ABS(skeb2_vrot))             &
              /SIZE(skeb2_vrot)

   WRITE(6,'("max(udiv)=",ES22.15,12X,"min(udiv)=",ES22.15)')           &
              MAXVAL(skeb2_udiv),MINVAL(skeb2_udiv)
   WRITE(6,'("maxloc(udiv)=",3I7,10X,"minloc(udiv)=",3I7)')             &
              MAXLOC(skeb2_udiv),MINLOC(skeb2_udiv)
   WRITE(6,'("mean(udiv)= ",ES22.15)') SUM(ABS(skeb2_udiv))             &
              /SIZE(skeb2_udiv)

   WRITE(6,'("max(vdiv)=",ES22.15,12X,"min(vdiv)=",ES22.15)')           &
              MAXVAL(skeb2_vdiv),MINVAL(skeb2_vdiv)
   WRITE(6,'("maxloc(vdiv)=",3I7,10X,"minloc(vdiv)=",3I7)')             &
              MAXLOC(skeb2_vdiv),MINLOC(skeb2_vdiv)
   WRITE(6,'("mean(vdiv)= ",ES22.15)') SUM(ABS(skeb2_vdiv))             &
              /SIZE(skeb2_vdiv)

   WRITE(6,'("max(sdisp)=",ES22.15,12X,"min(sdisp)=",ES22.15)')         &
              MAXVAL(sdisp),MINVAL(sdisp)
   WRITE(6,'("maxloc(sdisp)=",3I7,10X,"minloc(sdisp)=",3I7)')           &
              MAXLOC(sdisp),MINLOC(sdisp)
   WRITE(6,'("mean(sdisp)= ",ES22.15)') SUM(sdisp)/SIZE(sdisp)

   WRITE(6,'("max(cdisp)=",ES22.15,12X,"min(cdisp)=",ES22.15)')         &
              MAXVAL(cdisp),MINVAL(cdisp)
   WRITE(6,'("maxloc(cdisp)=",3I7,10X,"minloc(cdisp)=",3I7)')           &
              MAXLOC(cdisp),MINLOC(cdisp)
   WRITE(6,'("mean(cdisp)= ",ES22.15)') SUM(cdisp)/SIZE(cdisp)

   WRITE(6,'("max(kdisp)=",ES22.15,12X,"min(kdisp)=",ES22.15)')         &
              MAXVAL(kdisp),MINVAL(kdisp)
   WRITE(6,'("maxloc(kdisp)=",3I7,10X,"minloc(kdisp)=",3I7)')           &
              MAXLOC(kdisp),MINLOC(kdisp)
   WRITE(6,'("mean(kdisp)= ",ES22.15)') SUM(kdisp)/SIZE(kdisp)

 END IF


! SKEB2: Vert Integ. KE of total wind incr before SKEB2

 IF (stph_diag%l_skeb2_ke_prewindincr) THEN
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       vert_int_work(i,j) = 0.
     END DO
   END DO

! DEPENDS ON: swap_bounds
   CALL swap_bounds( r_u, row_length, rows,                             &
                     skeb2_toplev, offx, offy, fld_type_u, .TRUE.  )

! DEPENDS ON: swap_bounds
   CALL swap_bounds( r_v, row_length, n_rows,                           &
                     skeb2_toplev, offx, offy, fld_type_v, .TRUE.  )

   ! Interpolate squared velocity increments to pi-points
   DO k = skeb2_botlev, skeb2_toplev
     DO j = pdims%j_start, pdims%j_end
       jm1= j - 1
       ju = j
       jv = j
       DO i = pdims%i_start, pdims%i_end
         im1 = i - 1
         vert_int_work(i,j) = vert_int_work(i,j) + (0.5 * mass(i,j,k) * &
                    (0.5*(r_u(im1,ju,k)**2 + r_u(i,ju,k)**2) +          &
                     0.5*(r_v( i,jm1,k)**2 + r_v(i,jv,k)**2)))
       END DO !i
     END DO !j
   END DO !k
   IF (.NOT. l_endgame) THEN
     ! Set value at pole = zero
     IF (at_extremity(psouth)) THEN
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,pdims%j_start) = 0.
       END DO
     END IF
     IF (at_extremity(pnorth)) THEN
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,pdims%j_end) = 0.
       END DO
     END IF
   END IF

   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       stph_diag%skeb2_ke_prewindincr(i,j) = vert_int_work(i,j) *      &
                                              c_dadt(i,j)
     END DO
   END DO
 END IF


! --------------------------------------------------------------------

! Pass increments from SKEB2 to R_u, R_v

! --------------------------------------------------------------------
 DO k = skeb2_botlev, skeb2_toplev
   DO j = udims%j_start, udims%j_end
     DO i = udims%i_start, udims%i_end
       r_u(i,j,k) = r_u(i,j,k) + skeb2_urot(i,j,k) + skeb2_udiv(i,j,k)
     END DO ! i
   END DO ! j
   DO j = vdims%j_start, vdims%j_end
     DO i = vdims%i_start, vdims%i_end
       r_v(i,j,k) = r_v(i,j,k) + skeb2_vrot(i,j,k) + skeb2_vdiv(i,j,k)
     END DO ! i
   END DO ! j
 END DO ! k


! --------------------------------------------------------------------

! Output Stash Diagnostics

! --------------------------------------------------------------------

     
! u after skeb2
 IF (stph_diag%L_skeb2_u) then
   DO k = udims%k_start, udims%k_end
     DO j = udims%j_start, udims%j_end
       DO i = udims%i_start, udims%i_end
             stph_diag%skeb2_u(i,j,k) = u(i,j,k) + r_u(i,j,k)
       END DO
     END DO
   END DO
 END IF


! v after skeb2
 IF (stph_diag%L_skeb2_v) then
   DO k = vdims%k_start, vdims%k_end
     DO j = vdims%j_start, vdims%j_end
       DO i = vdims%i_start, vdims%i_end
             stph_diag%skeb2_v(i,j,k) = v(i,j,k) + r_v(i,j,k)
       END DO
     END DO
   END DO
 END IF

    
! u increment diagnostic
 IF (stph_diag%L_skeb2_u_incr) then
   DO k = udims%k_start, udims%k_end
     DO j = udims%j_start, udims%j_end
       DO i = udims%i_start, udims%i_end
         stph_diag%skeb2_u_incr(i,j,k) = r_u(i,j,k) -                   &
                     stph_diag%skeb2_u_incr(i,j,k)
       END DO
     END DO 
   END DO
 END IF


! v increment diagnostic
 IF (stph_diag%L_skeb2_v_incr) then
   DO k = vdims%k_start, vdims%k_end
     DO j = vdims%j_start, vdims%j_end
       DO i = vdims%i_start, vdims%i_end
         stph_diag%skeb2_v_incr(i,j,k) = R_v(i,j,k) -                   &
                     stph_diag%skeb2_v_incr(i,j,k)
       END DO 
     END DO 
   END DO
 END IF
                                                 

! rotational u increments from SKEB2
 IF (stph_diag%l_skeb2_u_rot) THEN
   DO k = udims%k_start, udims%k_end
     DO j = udims%j_start, udims%j_end
       DO i = udims%i_start, udims%i_end
         stph_diag%skeb2_u_rot(i,j,k) = skeb2_urot(i,j,k)
       END DO
     END DO
   END DO
 END IF


! rotational v increments from SKEB2
 IF (stph_diag%l_skeb2_v_rot) THEN
   DO k = vdims%k_start, vdims%k_end
     DO j = vdims%j_start, vdims%j_end
       DO i = vdims%i_start, vdims%i_end
         stph_diag%skeb2_v_rot(i,j,k) = skeb2_vrot(i,j,k)
       END DO
     END DO
   END DO
 END IF


! divergent u increments from SKEB2
 IF (stph_diag%l_skeb2_u_div) THEN
   DO k = udims%k_start, udims%k_end
     DO j = udims%j_start, udims%j_end
       DO i = udims%i_start, udims%i_end
         stph_diag%skeb2_u_div(i,j,k) = skeb2_udiv(i,j,k)
       END DO
     END DO
   END DO
 END IF


! divergent v increments from SKEB2
 IF (stph_diag%l_skeb2_v_div) THEN
   DO k = vdims%k_start, vdims%k_end
     DO j = vdims%j_start, vdims%j_end
       DO i = vdims%i_start, vdims%i_end
         stph_diag%skeb2_v_div(i,j,k) = skeb2_vdiv(i,j,k)
       END DO
     END DO
   END DO
 END IF


! dissipation from smagorinsky
 IF (stph_diag%l_skeb2_disp_smag) THEN
   DO k = pdims%k_start, pdims%k_end
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         stph_diag%skeb2_disp_smag(i,j,k) = sdisp(i,j,k)
       END DO
     END DO
   END DO
 END IF


! dissipation from convection
 IF (stph_diag%l_skeb2_disp_conv) THEN
   IF (skeb2_cdisp .EQ. type_cape) THEN
     WRITE(6,*)
     WRITE(6,'("**INFO**: SKEB2 CAPE Convective Dissipation")')
   ELSE IF (skeb2_cdisp .EQ. type_mflx) THEN
     WRITE(6,*)
     WRITE(6,'("**INFO**: SKEB2 MFLX Convective Dissipation")')
   END IF
   DO k = pdims%k_start, pdims%k_end
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         stph_diag%skeb2_disp_conv(i,j,k) = cdisp(i,j,k)
       END DO
     END DO
   END DO
 END IF


! Dissipation from SKEB1
 IF (stph_diag%l_skeb2_disp_skeb1) THEN
   DO k = pdims%k_start, pdims%k_end
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         stph_diag%skeb2_disp_skeb1(i,j,k) = kdisp(i,j,k)
       END DO
     END DO
   END DO
 END IF


! Smoothed dissipation field
 IF (stph_diag%l_skeb2_smodfield) THEN
   DO k = pdims%k_start, pdims%k_end
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         stph_diag%skeb2_smodfield(i,j,k) = sm_tot_disp(i,j,k)
       END DO
     END DO
   END DO
 END IF


! Streamfunction (modulated stream function forcing field)
 IF (stph_diag%l_skeb2_streamfunction) THEN
   DO k = pdims%k_start, pdims%k_end
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         stph_diag%skeb2_streamfunction(i,j,k) = m_psif(i,j,k)
       END DO
     END DO
   END DO
 END IF


! Streamfunction Forcing Field (initial from FFT)
 IF (stph_diag%l_skeb2_random_pattern) THEN
   DO k = pdims%k_start, pdims%k_end
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         stph_diag%skeb2_random_pattern(i,j,k) = psif(i,j,k)
       END DO
     END DO
   END DO
 END IF


! SKEB2: Mass-weighted Vert Integ. of modulated SF forcing field
 IF (stph_diag%l_skeb2_ke_psif) THEN
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       vert_int_work(i,j) = 0.
     END DO
   END DO
   DO k = skeb2_botlev, skeb2_toplev
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,j)=vert_int_work(i,j) +                        &
                 (m_psif(i,j,k) * mass(i,j,k))
       END DO
     END DO
   END DO
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       stph_diag%skeb2_ke_psif(i,j) = vert_int_work(i,j) * c_dadt(i,j)
     END DO
   END DO
 END IF


! SKEB2: Mass-weighted Vert Integ. of numerical diss
 IF (stph_diag%l_skeb2_ke_sdisp) THEN
   IF (.NOT.l_skeb2_psisdisp) THEN
     WRITE(6,*)
     WRITE(6,'("**WARNING**: SKEB2 Numerical Dissipation off")')
     WRITE(6,'("  Requested STASH diagnostic (0,35,16) not available")')
   ELSE
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,j) = 0.
       END DO
     END DO
     DO k = skeb2_botlev, skeb2_toplev
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
           vert_int_work(i,j)=vert_int_work(i,j) +                      &
                  (sdisp(i,j,k) * mass(i,j,k))
         END DO
       END DO
     END DO
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
         stph_diag%skeb2_ke_sdisp(i,j) = vert_int_work(i,j)
       END DO
     END DO

! Print global mean value of Numerical Energy Dissipation
     IF (l_skebprint) THEN

       CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,    &
                           gltotke_tmp)

       gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
       WRITE(6,*)
       WRITE(6,'("***** SKEB2 Global-total sdisp *****")')
       WRITE(6,'("tot(W/m^2)= ",ES22.15)') gltotke
     END IF    ! l_skebprint
   END IF
 END IF


! SKEB2: Mass-weighted Vert Integ. of convection diss
 IF (stph_diag%l_skeb2_ke_cdisp) THEN
   IF (skeb2_cdisp .NE. type_cape) THEN
     WRITE(6,*)
     WRITE(6,'("**WARNING**: SKEB2 CAPE Convective Dissipation off")')
     WRITE(6,'("  Requested STASH diagnostic (0,35,17) not available")')
   ELSE
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,j) = 0.
       END DO
     END DO
     DO k = skeb2_botlev, skeb2_toplev
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
           vert_int_work(i,j)=vert_int_work(i,j) +                      &
                    (cdisp(i,j,k) * mass(i,j,k))
         END DO
       END DO
     END DO
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
         stph_diag%skeb2_ke_cdisp(i,j) = vert_int_work(i,j)
       END DO
     END DO

! Print global mean value of Numerical Energy Dissipation
     IF (l_skebprint) THEN
       CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,    &
                           gltotke_tmp)

       gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
       WRITE(6,*)
       WRITE(6,'("***** SKEB2 Global-total cdisp *****")')
       WRITE(6,'("tot(W/m^2)= ",ES22.15)') gltotke
     END IF    ! l_skebprint
   END IF
 END IF


! SKEB2: Mass-weighted Vert Integ. of SKEB1 KE dissipation
 IF (stph_diag%l_skeb2_ke_kdisp) THEN
   IF (.NOT.l_skeb2_skeb1disp) THEN
     WRITE(6,*)
     WRITE(6,'("**WARNING**: SKEB1 Dissipation l_skeb2_skeb1disp=F")')
     WRITE(6,'("  Requested STASH diagnostic (0,35,18) not available")')
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         stph_diag%skeb2_ke_kdisp(i,j) = 0.
       END DO
     END DO
   ELSE
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,j) = 0.
       END DO
     END DO
     DO k = skeb2_botlev, skeb2_toplev
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
           vert_int_work(i,j)=vert_int_work(i,j) +                      &
                  (kdisp(i,j,k) * mass(i,j,k))
         END DO
       END DO
     END DO
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
         stph_diag%skeb2_ke_kdisp(i,j) = vert_int_work(i,j)
       END DO
     END DO

! Print global mean value of Numerical Energy Dissipation
     IF (l_skebprint) THEN
       CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,    &
                           gltotke_tmp)

       gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
       WRITE(6,*)
       WRITE(6,'("***** SKEB2 Global-total kdisp *****")')
       WRITE(6,'("tot(W/m^2)= ",ES22.15)') gltotke
     END IF    ! l_skebprint
   END IF
 END IF


! SKEB2: Mass-weighted Vert Integ. of smoothed dissipation field
! Values for this diagnostic are printed to PE output for monitoring
! purposes
 DO j = pdims%j_start, pdims%j_end
   DO i = pdims%i_start, pdims%i_end
     vert_int_work(i,j) = 0.
   END DO
 END DO
 DO k = skeb2_botlev, skeb2_toplev
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       vert_int_work(i,j)=vert_int_work(i,j) +                        &
                (sm_tot_disp(i,j,k) * mass(i,j,k))
     END DO
   END DO
 END DO
 DO j = pdims%j_start, pdims%j_end
   DO i = pdims%i_start, pdims%i_end
     vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
   END DO
 END DO

! Print global mean value of Numerical Energy Dissipation
 IF (l_skebprint) THEN
   CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,    &
                       gltotke_tmp)

   gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
   WRITE(6,*)
   WRITE(6,'("***** SKEB2 Global-total energy dissipation *****")')
   WRITE(6,'("tot(W/m^2)= ",ES22.15)') gltotke
 END IF    ! l_skebprint
 IF (stph_diag%l_skeb2_ke_m_psif) THEN
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       stph_diag%skeb2_ke_m_psif(i,j) = vert_int_work(i,j)
     END DO
   END DO
 END IF


! SKEB2: Mass-weighted Vert Integ. KE of total wind incr after SKEB2
 IF (stph_diag%l_skeb2_ke_postwindincr) THEN
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       vert_int_work(i,j) = 0.
     END DO
   END DO

! DEPENDS ON: swap_bounds
    CALL swap_bounds( r_u, row_length, rows,                            &
                      skeb2_toplev, offx, offy, fld_type_u, .TRUE.  )

! DEPENDS ON: swap_bounds
    CALL swap_bounds( r_v, row_length, n_rows,                          &
                      skeb2_toplev, offx, offy, fld_type_v, .TRUE.  )

   ! Interpolate squared windspeed (i.e. scalar field) to P-gridpoints
   !  (to avoid over-smoothing of a vector field)
   DO k = skeb2_botlev, skeb2_toplev
     DO j = pdims%j_start, pdims%j_end
       jm1 = j - 1
       ju = j
       jv = j
       DO i = pdims%i_start, pdims%i_end
         im1 = i - 1
         vert_int_work(i,j) = vert_int_work(i,j) +                      &
                    (0.5 * mass(i,j,k) *                                &
                    (0.5*(r_u(im1,ju,k)**2 + r_u(i,ju,k)**2) +          &
                     0.5*(r_v( i,jm1,k)**2 + r_v(i,jv,k)**2)))
       END DO !I
     END DO !J
   END DO !K
   IF (.NOT. l_endgame) THEN
     ! Set value at pole = zero
     IF (at_extremity(psouth)) THEN
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,pdims%j_start) = 0.
       END DO
     END IF
     IF (at_extremity(pnorth)) THEN
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,pdims%j_end) = 0.
       END DO
     END IF
   END IF

   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       stph_diag%skeb2_ke_postwindincr(i,j) = vert_int_work(i,j) *      &
                                              c_dadt(i,j)
     END DO
   END DO
 END IF


! SKEB2: Vert Integ. KE of wind incr from SKEB2
 IF (stph_diag%l_skeb2_ke_windincr) THEN
   DO k = udims%k_start, udims%k_end
     DO j = udims%j_start, udims%j_end
       DO i = udims%i_start, udims%i_end
         skeb2_urot0(i,j,k) = skeb2_urot(i,j,k)
         skeb2_udiv0(i,j,k) = skeb2_udiv(i,j,k)
       END DO
     END DO
     DO j = vdims%j_start, vdims%j_end
       DO i = vdims%i_start, vdims%i_end
         skeb2_vrot0(i,j,k) = skeb2_vrot(i,j,k)
         skeb2_vdiv0(i,j,k) = skeb2_vdiv(i,j,k)
       END DO
     END DO
   END DO

! Use Multi-variable swap-bounds to exploit MPP
   ifield = 0
   ifield = ifield + 1
   fields_to_swap(ifield) % field       => skeb2_urot0(:,:,:)
   fields_to_swap(ifield) % field_type  =  fld_type_u
   fields_to_swap(ifield) % levels      =  model_levels
   fields_to_swap(ifield) % rows        =  rows
   fields_to_swap(ifield) % vector      =  .TRUE.

   ifield = ifield + 1
   fields_to_swap(ifield) % field       => skeb2_udiv0(:,:,:)
   fields_to_swap(ifield) % field_type  =  fld_type_u
   fields_to_swap(ifield) % levels      =  model_levels
   fields_to_swap(ifield) % rows        =  rows
   fields_to_swap(ifield) % vector      =  .TRUE.

   ifield = ifield + 1
   fields_to_swap(ifield) % field       => skeb2_vrot0(:,:,:)
   fields_to_swap(ifield) % field_type  =  fld_type_v
   fields_to_swap(ifield) % levels      =  model_levels
   fields_to_swap(ifield) % rows        =  n_rows
   fields_to_swap(ifield) % vector      =  .TRUE.

   ifield = ifield + 1
   fields_to_swap(ifield) % field       => skeb2_vdiv0(:,:,:)
   fields_to_swap(ifield) % field_type  =  fld_type_v
   fields_to_swap(ifield) % levels      =  model_levels
   fields_to_swap(ifield) % rows        =  n_rows
   fields_to_swap(ifield) % vector      =  .TRUE.

! DEPENDS ON: swap_bounds_mv
   CALL swap_bounds_mv( fields_to_swap, ifield, row_length,             &
                           offx, offy )

! KE = 0.5*M*V**2, but because the u and v components are on different
!      grids, we convert it all to the P-grid (see diagram at top)
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       vert_int_work(i,j) = 0.
     END DO
   END DO
   DO k = skeb2_botlev, skeb2_toplev
     DO j = pdims%j_start, pdims%j_end
       jm1 = j - 1
       ju = j
       jv = j
       DO i = pdims%i_start, pdims%i_end
         im1= i - 1
         vert_int_work(i,j) = vert_int_work(i,j) +                      &
              (0.5 * mass(i,j,k) *                                      &
              (0.5*((skeb2_urot0(im1,ju,k)+ skeb2_udiv0(im1,ju,k))**2 + &
                    (skeb2_urot0(i,ju,k)  + skeb2_udiv0(i,ju,k))**2)  + &
               0.5*((skeb2_vrot0(i,jm1,k) + skeb2_vdiv0(i,jm1,k))**2  + &
                    (skeb2_vrot0(i,jv,k)  + skeb2_vdiv0(i,jv,k))**2)))
       END DO !I
     END DO !J
   END DO !K
   IF (.NOT. l_endgame) THEN
     ! Set value at pole = zero (single point with no mass)
     IF (at_extremity(psouth)) THEN
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,pdims%j_start) = 0.
       END DO
     END IF
     IF (at_extremity(pnorth)) THEN
       DO i = pdims%i_start, pdims%i_end
         vert_int_work(i,pdims%j_end) = 0.
       END DO
     END IF
   END IF

   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
       stph_diag%skeb2_ke_windincr(i,j) = vert_int_work(i,j)
     END DO
   END DO

! Print global mean value of Numerical Energy Dissipation
   IF (l_skebprint) THEN
     CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,    &
                         gltotke_tmp)

     gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
     WRITE(6,*)
     WRITE(6,'("***** SKEB2 Global-total incr KE *****")')
     WRITE(6,'("tot(W/m^2)= ",ES22.15)') gltotke
   END IF
 END IF

! Deallocate work arrays
 IF (ALLOCATED(vert_int_work)) DEALLOCATE(vert_int_work)

 IF (sf(0,35)) THEN
   CALL  diagnostics_stph( row_length, rows, model_levels,              &
                           n_rows, at_extremity, stph_diag,             &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                           stashwork35)
   ! ------------------------
   ! Tidy allocatable arrays 
   ! ------------------------
   DEALLOCATE(stph_diag%skeb2_u)
   DEALLOCATE(stph_diag%skeb2_v)
   DEALLOCATE(stph_diag%skeb2_u_incr)
   DEALLOCATE(stph_diag%skeb2_v_incr)
   DEALLOCATE(stph_diag%skeb2_u_rot)
   DEALLOCATE(stph_diag%skeb2_v_rot)                          
   DEALLOCATE(stph_diag%skeb2_u_div)  
   DEALLOCATE(stph_diag%skeb2_v_div)
   DEALLOCATE(stph_diag%skeb2_disp_smag)  
   DEALLOCATE(stph_diag%skeb2_disp_conv)
   DEALLOCATE(stph_diag%skeb2_disp_skeb1)
   DEALLOCATE(stph_diag%skeb2_smodfield)
   DEALLOCATE(stph_diag%skeb2_streamfunction)
   DEALLOCATE(stph_diag%skeb2_random_pattern)
   DEALLOCATE(stph_diag%skeb2_ke_psif)
   DEALLOCATE(stph_diag%skeb2_ke_sdisp)
   DEALLOCATE(stph_diag%skeb2_ke_cdisp)
   DEALLOCATE(stph_diag%skeb2_ke_kdisp)
   DEALLOCATE(stph_diag%skeb2_ke_m_psif)
   DEALLOCATE(stph_diag%skeb2_ke_prewindincr)
   DEALLOCATE(stph_diag%skeb2_ke_windincr)
   DEALLOCATE(stph_diag%skeb2_ke_postwindincr)
 END IF
 IF (ALLOCATED(skeb2_up_flux)) DEALLOCATE(skeb2_up_flux)
 IF (ALLOCATED(skeb2_dwn_flux)) DEALLOCATE(skeb2_dwn_flux)
 IF (ALLOCATED(skeb2_cape)) DEALLOCATE(skeb2_cape)

 ! Append dpsidtc/s to seed file (if active) at dump times.
 ! Generally used for CRUNs.
 IF (ldump) THEN
   IF (l_stphseed_read.OR.l_stphseed_write) THEN
     IF (mype == 0) THEN
       ! Get latest seed value
       ! All calls to random for this timestep are complete at this point
       CALL random_seed(SIZE=tam)
       IF (.NOT.ALLOCATED(iranseed)) ALLOCATE(iranseed(tam))
       CALL random_seed(GET=iranseed(1:tam))

       CALL stph_openoutput("rewind    ", .true.)

       CALL stph_writeentry2(iranseed, tam)
       DEALLOCATE(iranseed)

       CALL stph_writeentry2(dpsidtc, SIZE(dpsidtc))

       CALL stph_writeentry2(dpsidts, SIZE(dpsidts))

       CALL stph_closeoutput()
     END IF
   END IF
 END IF

 IF (lhook) CALL dr_hook('STPH_SKEB2',zhook_out,zhook_handle)
 RETURN

! --------------------------------------------------------------------

END SUBROUTINE stph_skeb2
END MODULE stph_skeb2_mod
