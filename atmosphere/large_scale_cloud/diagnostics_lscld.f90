! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Diagnostic output area for Large Scale Cloud (Section 9) Diagnostics.
! Subroutine Interface:
      Subroutine diagnostics_lscld(                                     &
                             row_length, rows, model_levels,            &
                             rhc_row_length, rhc_rows,                  &
                             wet_model_levels, boundary_layer_levels,   &
                             cloud_levels,                              &
                             n_rows, global_row_length, global_rows,    &
                             halo_i, halo_j, off_x, off_y, me,          &
                             n_proc, n_procx, n_procy,                  &
                             g_rows, g_row_length,                      &
                             at_extremity, p_theta_levels,              &
                             p,                                         &
                             T_incr_diagnostic, q_incr_diagnostic,      &
                             qcl_incr_diagnostic,                       &
                             T, q, qcl, qcf,                            &
                             area_cloud_fraction, bulk_cloud_fraction,  &
                             cfl, cff, rho,                             &
                             p_star, rhcpt, combined_cloud, L_murk,     &
                             Aerosol, RHcrit, l_mixing_ratio,           &
                             delta_lambda, delta_phi,                   &
                             FV_cos_theta_latitude,                     &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                             STASHwork,                                 &
                             l_combi_cld                                &
                             )

USE atm_fields_bounds_mod, ONLY:                                        &
    tdims, tdims_s, qdims

USE conversions_mod, ONLY: &
    zerodegc,              &
    ft2m

USE level_heights_mod, ONLY:                                            &
    r_theta_levels, r_rho_levels

Use ac_diagnostics_mod, Only :  &
    cf_lsc, qcl_lsc

USE earth_constants_mod, ONLY: earth_radius
USE pc2_constants_mod,   ONLY: cloud_pc2_tol_2,condensate_limit
USE cloud_inputs_mod,    ONLY: l_filter_cloud, tau_thresh
USE rad_pcf,             ONLY: ip_cloud_mix_max, ip_cloud_mix_random

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE Control_Max_Sizes
USE c_lowcld_mod, ONLY: str_ceil, cloud_threshold
USE thetaw_mod, ONLY: thetaw
USE gammaf_mod, ONLY: gammaf,  gammaf_lower_incomp
USE Submodel_Mod

IMPLICIT NONE
!
! Purpose:
!          Calculates diagnostics associated with large scale cloud
! (section 9) for output through STASH system. This routine is
! positioned at the end of moist processes calculation, ie the end
! of the moist timestep, and is therefore the appropriate location
! for output of moisture related diagnostics.
!
! Method:
!  Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the calling
! routine. After minor calculations, each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing.
!
!  Diagnostics currently available: (in order calculated)
!
! STASH item (all section 9 )
!   4 Temperature on model levels after large scale cloud
! 181 T   increment over boundary layer + large scale cloud routines
! 182 q   increment over boundary layer + large scale cloud routines
! 183 qcl increment over boundary layer + large scale cloud routines
! 201 bulk cloud fraction
! 202 very low cloud amount
! 203 low      cloud amount
! 204 medium   cloud amount
! 205 high     cloud amount
! 206 qcl after large scale cloud. Note this is only needed to provide
!     assimilation increments for latent heat nudging. Otherwise should
!     use prognostic (0,254).
! 208 cloud base for cover >  0.1 octa kft
! 209 cloud base for cover >  1.5 octa kft
! 210 cloud base for cover >  2.5 octa kft
! 211 cloud base for cover >  3.5 octa kft
! 212 cloud base for cover >  4.5 octa kft
! 213 cloud base for cover >  5.5 octa kft
! 214 cloud base for cover >  6.5 octa kft
! 215 cloud base for cover >  7.9 octa kft
! 216 total cloud ; random overlap
! 217 total cloud ; max/random overlap
! 218 cloud fraction below 1000 ft asl
! 219 low cloud base ft asl
! 220 low cloud top ft asl
! 221 wet bulb freezing level
! 222 wet bulb temperature
! 223 total cloud top height kft
! 229 relative humidity on model levels. Note: percentage units.
! 226 layer cloud frequency
! 228 critical relative humidity on model levels
! 230 Visibility on model levels. Note that if murk aerosol is included
!     then vis depends on aerosol as well as humidity. Visibility is
!     available over all model levels but aerosol is only actively
!     mixed over boundary levels.
! 231 combined cloud amount on model levels
! 232 total cloud ; random overlap
!     (cloud above ceilometer range filtered out)
! 233 total cloud ; max/random overlap
!     (cloud above ceilometer range filtered out)
! 234 combined cloud amount on model levels
!     (cloud above ceilometer range filtered out)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: large scale cloud
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
!  Global Variables:----------------------------------------------------
!     None.
!
! Arguments with Intent IN. ie: Input variables.

      LOGICAL ,INTENT(IN) ::                                            &
        at_extremity(4),                                                &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
        l_combi_cld      ! Indicates whether any of the combined cloud
                         ! diagnostics are requested.

      INTEGER ,INTENT(IN) ::                                            &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, rhc_row_length                                                  &
                         ! row_length if diagnostic RHcrit ON, else 1.
     &, rhc_rows                                                        &
                         ! rows       if diagnostic RHcrit ON, else 1.
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, cloud_levels                                                    &
                         ! number of cloud model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, boundary_layer_levels
                         ! number of boundary layer levels


      INTEGER ,INTENT(IN) ::                                            &
     &  global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows                                                     &
                            !IN. NUMBER OF global rows
     &, me                                                              &
                            !IN. Processor number
     &, halo_i                                                          &
                            !IN. size of large halo in x direction
     &, halo_j                                                          &
                            !IN. size of large halo in y direction
     &, off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y                                                           &
                            !IN. size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       

      LOGICAL ,INTENT(IN) ::                                            &
     &  L_murk                                                          &
                            ! murk aerosol included if T
     &, L_mixing_ratio      ! Use mixing ratio formulation

! Primary Arrays used in all models
      Real ,INTENT(IN) ::                                               &
        p_theta_levels(row_length, rows, model_levels),                 &
        p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
            model_levels),                                              &
        T(row_length, rows, model_levels),                              &
        q(row_length,rows, wet_model_levels),                           &
        qcl(row_length, rows, wet_model_levels),                        &
        qcf(row_length, rows, wet_model_levels),                        &
        bulk_cloud_fraction(row_length, rows, wet_model_levels),        &
        cfl(qdims%i_start:qdims%i_end,                                  &
            qdims%j_start:qdims%j_end,                                  &
                        1:qdims%k_end),                                 &
        cff(qdims%i_start:qdims%i_end,                                  &
            qdims%j_start:qdims%j_end,                                  &
                        1:qdims%k_end),                                 &
        rho(tdims%i_start:tdims%i_end,                                  &
            tdims%j_start:tdims%j_end,                                  &
                        1:tdims%k_end-1),                               &
        p_star(row_length, rows),                                       &
        RHCPT(rhc_row_length, rhc_rows, wet_model_levels),              &

        T_incr_diagnostic(row_length, rows, model_levels),              &
        q_incr_diagnostic(row_length,rows, wet_model_levels),           &
        qcl_incr_diagnostic(row_length, rows, wet_model_levels),        &
        rhcrit(wet_model_levels),                                       &
                                    ! critical RH for cloud formation
        Aerosol(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
                model_levels),                                          &
                                    ! murk aerosol
        FV_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end)

      REAL ,INTENT(INOUT) ::                                            &
        area_cloud_fraction(row_length, rows, wet_model_levels),        &
                                    ! Area cloud fraction
        combined_cloud(row_length, rows, wet_model_levels)
                                    ! Combined cloud fraction using area
                                    ! cloud fraction and convective.

! Grid-size information
      REAL, INTENT(IN) :: delta_lambda, delta_phi

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
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA.
!
      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
        h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
        HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end
!         cconsts defines LOW_BOT_LEVEL to HIGH_TOP_LEVEL cloud splits.
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
!
! arguments with intent in/out. ie: input variables changed on output.
!
! Diagnostics info
      Real                                                              &
     & STASHwork(*)     ! STASH workspace
!
!  Local parameters and other physical constants------------------------
!
      Integer sect ! STASH section for diagnostics
      Parameter ( sect = 9 )
!      ( LS Cloud diagnostics are output under STASH Section 09 )
!
      CHARACTER(LEN=*) RoutineName
      Parameter  ( RoutineName='diagnostics_lscld')
!
!  Local scalars--------------------------------------------------------
!
      Integer                                                           &
     &  i, j, k, l, ji                                                  &
     &, kinvert                                                         &
                                ! vertical index for inverted arrays.
     &, item                                                            &
                                ! STASH item of individual diagnostics.
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      REAL, PARAMETER:: vis_probability=0.5 ! median visibility diag

      CHARACTER(LEN=80)  cmessage
!
      Integer                                                           &
     &  im_index                                                        &
                        ! internal model index
!                  ( LS Cloud is part of the Atmosphere internal model )
     &, nclds      ! Number of radiation cloud levels ( <=  wet levels)
!  Local scalars required to calculate Diagnostics 09208 - 09215
      INTEGER, PARAMETER :: noktas=8
      REAL,    PARAMETER :: roktas=8.0
      REAL    :: cloud_base_height         ! Height of cloud base (m)
      REAL    :: combined_cloud_in_oktas   ! Combined cloud in oktas

!  Local dynamic arrays required to calculate Diagnostics 09208 - 09215
      INTEGER ::                                                        &
     & cloud_base_level(row_length, rows, noktas)
      REAL ::                                                           &
     & cloud_base(row_length, rows, noktas)
      REAL :: oktas(noktas)
      DATA oktas/0.1,1.5,2.5,3.5,4.5,5.5,6.5,7.9/

!  Local scalars required to calculate Diagnostics 09218, 09219, 09220
      REAL ::                                                           &
     & pu,pl                                                            &
                                       ! Upper and lower half
                                       !       level pressure
     &,pt                                                               &
                                       ! Pressure thickness accumulator
     &,ft                                                               &
                                       ! Cloud fraction accumulator
     &,dp                                                               &
                                       ! Later pressure thickness
     &,h_asl                                                            &
                                       ! Layer base height asl
     &,h_asln                                                           &
                                       ! Layer top height asl
     &,fr                              ! Layer fraction below ceiling

      REAL, PARAMETER ::                                                &
     & str_ceilm = str_ceil * ft2m     ! STR_CEIL in metres

!  Local scalars required to calculate Diagnostic 09221 and 09222
      LOGICAL, PARAMETER ::                                             &
     & l_potential=.false.             ! Wet bulb temperature required
                                       ! from subroutine ThetaW
      REAL ::                                                           &
     & frac

!  Local scalars required to calculate Diagnostic 09223
      INTEGER ::                                                        &
     & cld_top_lev                     ! Top layer in cloud
      REAL, PARAMETER::                                                 &
     & thresh = 0.0627,                                                 &
     & m_to_kft = (1./ft2m)*0.001

!  Local scalars required to calculate Diagnostics 09232 - 09234
      REAL, PARAMETER :: ceilometer_range = 6000.0
                                       ! Max range of ceilometers used
                                       ! in filtered cloud diagnostics
!
! For aircraft icing algorithms (Icing Index: ii)
      REAL, PARAMETER    :: iitmin       = zerodegc-20.0 ! Min temp for icing
      REAL, PARAMETER    :: iitmax       = zerodegc      ! Max temp for icing
      REAL, PARAMETER    :: qcl_thresh   = 1.0e-4        ! Critical LWC
      REAL               :: mean_qcl    ! In-cloud mean qcl (in kg/m3)
      REAL               :: variance_qcl! Variance of qcl
      REAL               :: tmp         ! Temporary storage variable
      REAL               :: mu          ! For calculating log-normal PDF
      REAL               :: sigma       ! For calculating log-normal PDF
      REAL               :: focat       ! Fraction of Cloud Above Threshold 
      REAL               :: erf_in(1)   ! Input to Error function routine
      REAL               :: erf_out(1)  ! Output from Error function routine
      REAL               :: kay         ! For calculating log-normal PDF
      REAL               :: thet        ! For calculating log-normal PDF
      REAL               :: BigGamma    ! Output from complete gamma fn calc
      REAL               :: LittleGamma ! Output from lower incomplete 
                                        ! gamma fn calc
      REAL               :: horiz_scale ! Size of the grid-box for finding
                                        ! FSD of LWC.
      REAL               :: phic        ! Used in calculating FSD
      REAL               :: fsd         ! Fractional Standard Deviation of LWC
      INTEGER            :: my_index, ref_item, v_item
!
! For filtering of optically thin ice cloud
      REAL               :: ice_fraction! Fraction cloud which is ice (not liq)
      REAL               :: tau         ! Optical depth
      REAL               :: filtered_ice_cloud 
                                        ! Cloud fraction after filter 
!
!  Local dynamic arrays-------------------------------------------------
       Real                                                             &
        work2d_1(row_length, rows),                                     &
                                        ! Single-level work array (cloud)
        work2d_2(row_length, rows),                                     &
        work2d_3(row_length, rows),                                     &
        work3d_1(row_length, rows, wet_model_levels)
                                        ! Full work array (cloud/moist).

       REAL, ALLOCATABLE :: original_combined_cloud(:,:,:)
! Copy of original combined cloud amount (before any filtering to remove 
! optically thin ice cloud) for checking that new filtered version
! does not lead to an increase in total cloud amount.
      REAL, ALLOCATABLE :: comb_ceilometer_cloud(:,:,:)
                                       ! Mixed CCA and CF per gridbox
                                       ! without cloud above ceilometer_range

!        filtered_area_cloud_fraction(row_length, rows, wet_model_levels)
! Area cloud amount filtered to remove optically thin ice cloud

       REAL, ALLOCATABLE :: icing3d(:,:,:)

       INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
       INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
       REAL(KIND=jprb)               :: zhook_handle
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
      IF (lhook) CALL dr_hook('DIAGNOSTICS_LSCLD',zhook_in,zhook_handle)

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

      IF (l_filter_cloud) THEN

        IF (l_combi_cld) THEN
          ! At least one of the diagnostics which uses combined cloud fraction
          ! has being requested.
          ALLOCATE ( original_combined_cloud( (qdims%i_end-qdims%i_start)+1, &
                                              (qdims%j_end-qdims%j_start)+1, &
                                               qdims%k_end ) )

          DO k = 1, qdims%k_end
            DO j = qdims%j_start,qdims%j_end
              DO i = qdims%i_start,qdims%i_end
                original_combined_cloud(i,j,k) = combined_cloud(i,j,k)
              END DO
            END DO
          END DO

          ! Filter the combined cloud amount by removing sub-visual
          ! (optically thin) ice cloud before calculating other diagnostics
          ! such as total cloud amount.
          DO k = 1, qdims%k_end
            DO j = qdims%j_start,qdims%j_end
              DO i = qdims%i_start,qdims%i_end

                kinvert = wet_model_levels+1-k

                IF ( combined_cloud(i,j,kinvert)>cloud_pc2_tol_2 )THEN
                  ! Extinction (m-1) taken from Heymsfield et al 2003.
                  sigma=(((qcf(i,j,k)*rho(i,j,k)*1000.0) / &
                        combined_cloud(i,j,kinvert))**0.9)*0.02

                  ! Optical depth=integral of sigma over depth of model level.
                  tau = sigma*(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))

                  IF (tau < tau_thresh) THEN
                    filtered_ice_cloud = 0.0
                  ELSE
                    filtered_ice_cloud = combined_cloud(i,j,kinvert)
                  END IF

                  ice_fraction = cff(i,j,k)/combined_cloud(i,j,kinvert)

                  combined_cloud(i,j,kinvert)=(ice_fraction*filtered_ice_cloud)&
                           +((1.0-ice_fraction)*combined_cloud(i,j,kinvert))
                ELSE
                  ! If there is practically no cloud
                  ! set filtered combined cloud to zero.
                  combined_cloud(i,j,kinvert) = 0.0
                END IF

              END DO
            END DO
          END DO

        END IF ! l_combi_cld

        IF ( sf(202,9) .OR. sf(203,9) .OR. sf(204,9) .OR. sf(205,9) .OR. &
             sf(218,9) .OR. sf(219,9) .OR. sf(220,9) .OR. sf(226,9) ) THEN
          ! At least one of the diagnostics which uses area cloud fraction
          ! has being requested.

          ! Filter the area cloud amount by removing sub-visual
          ! (optically thin) ice cloud before calculating other diagnostics
          ! such as very low, low, medium and high cloud amount.
          DO k = 1, qdims%k_end
            DO j = qdims%j_start,qdims%j_end
              DO i = qdims%i_start,qdims%i_end

                IF ( area_cloud_fraction(i,j,k)>cloud_pc2_tol_2 ) THEN

                  ! Extinction (m-1) taken from Heymsfield et al 2003.
                  sigma=(((qcf(i,j,k)*rho(i,j,k)*1000.0) / &
                        area_cloud_fraction(i,j,k))**0.9)*0.02

                  ! Optical depth=integral of sigma over depth of model level.
                  tau = sigma*(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))

                  IF (tau < tau_thresh) THEN
                    filtered_ice_cloud = 0.0
                  ELSE
                    filtered_ice_cloud = area_cloud_fraction(i,j,k)
                  END IF

                  ice_fraction = cff(i,j,k)/area_cloud_fraction(i,j,k)

                  area_cloud_fraction(i,j,k) = &
                           (ice_fraction*filtered_ice_cloud) &
                           +((1.0-ice_fraction)*area_cloud_fraction(i,j,k))
                ELSE
                  ! If there is practically no cloud
                  ! set filtered area cloud to zero.
                  area_cloud_fraction(i,j,k) = 0.0
                END IF

              END DO
            END DO
          END DO

        END IF ! sf(202,9) .OR. sf(203,9) .OR. sf(204,9) .OR. sf(205,9) etc.

      END IF ! l_filter_cloud
!
! ----------------------------------------------------------------------
! DIAG.09004 Copy T from Main Cloud to stashwork
! ----------------------------------------------------------------------
      item = 4
! Diag09004_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),T,         &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
         If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(T Main Cloud)"
         End if
      End if  ! Diag09004_if1
!
! ----------------------------------------------------------------------
! DIAG.09181 Copy T INC: bdy layer + ls cld to stashwork
! ----------------------------------------------------------------------
      item = 181
! Diag09181_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       T_incr_diagnostic,                                         &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(T INC: bl + ls cld)"
        End if
      End if  ! Diag09181_if1
!
! ----------------------------------------------------------------------
! DIAG.09182 Copy q INC: bdy layer + ls cld to stashwork
! ----------------------------------------------------------------------
      item = 182
! Diag09182_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       q_incr_diagnostic,                                         &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Q INC: bl + ls cld)"
        End if
      End if  ! Diag09182_if1
!
! ----------------------------------------------------------------------
! DIAG.09183 Copy qCL INC: bdy layer + ls cld to stashwork
! ----------------------------------------------------------------------
      item = 183
! Diag09183_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d ( stashwork(si(item,sect,im_index)),           &
     &       qcl_incr_diagnostic,                                       &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
        cmessage="cld_ctl  : error in copydiag_3d(Qcl INC: bl + ls cld)"
        End if
      End if  ! Diag09183_if1
!
! ----------------------------------------------------------------------
! DIAG.09201 Copy Bulk cloud fraction to stashwork
! ----------------------------------------------------------------------
! Needed for Assimilation, otherwise duplicates 00266 output (preferred)
      item = 201
! Diag09201_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then

      If (.not.allocated(cf_lsc)) then
        Allocate ( cf_lsc(row_length*rows,wet_model_levels) )
      End If
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        bulk_cloud_fraction,                                      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
! Copy bulk_cloud_fraction into cf_lsc for ac_diagnostics_mod

      do k = 1,wet_model_levels
        do j = 1,rows
          do i = 1,row_length
            ji = (j-1)*row_length+i
            cf_lsc(ji,k) = bulk_cloud_fraction(i,j,k)
          end do
        end do
      end do

        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Bulk CF Main Cloud)"
        End if
      End if  ! Diag09201_if1
!
! ----------------------------------------------------------------------
! DIAG.09202 Copy Very Low Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 202
! Diag09202_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09202_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i, j) = area_cloud_fraction(i, j, 1)
          End Do
        End Do  ! Diag09203_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09202_do2:
        Do k = 1, MAX((LOW_BOT_LEVEL - 1), 1)
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
                   MAX(work2d_1(i,j),area_cloud_fraction(i,j,k))
            End Do
          End Do
        End Do  ! Diag09202_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(LOW Cloud Amount)"
         End if
      End if  ! Diag09203_if1
!
! ----------------------------------------------------------------------
! DIAG.09203 Copy Low Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 203
! Diag09203_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09203_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i,j) = area_cloud_fraction(i,j,low_bot_level)
          End Do
        End Do  ! Diag09203_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09203_do2:
        Do k = LOW_BOT_LEVEL + 1, LOW_TOP_LEVEL
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
                   MAX(work2d_1(i,j),area_cloud_fraction(i,j,k))
            End Do
          End Do
        End Do  ! Diag09203_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(LOW Cloud Amount)"
         End if
      End if  ! Diag09203_if1
!
! ----------------------------------------------------------------------
! DIAG.09204 Copy Medium Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 204
! Diag09204_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09204_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i,j) = area_cloud_fraction(i,j,med_bot_level)
          End Do
        End Do  ! Diag09204_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09204_do2:
        Do k = MED_BOT_LEVEL + 1, MED_TOP_LEVEL
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
                   MAX(work2d_1(i,j),area_cloud_fraction(i,j,k))
            End Do
          End Do
        End Do  ! Diag09204_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(MEDIUM Cloud Amount)"
        End if
      End if  ! Diag09204_if1
!
! ----------------------------------------------------------------------
! DIAG.09205 Copy High Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 205
! Diag09205_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09205_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i,j) = area_cloud_fraction(i,j,high_bot_level)
          End Do
        End Do  ! Diag09205_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09205_do2:
        Do k = HIGH_BOT_LEVEL + 1, HIGH_TOP_LEVEL
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
                   MAX(work2d_1(i,j),area_cloud_fraction(i,j,k))
            End Do
          End Do
        End Do  ! Diag09205_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(HIGH Cloud Amount)"
        End if
      End if  ! Diag09205_if1
!
! ----------------------------------------------------------------------
! DIAG.09206 Copy cloud liquid condensate to stashwork
! ----------------------------------------------------------------------
! Needed for Assimilation, otherwise duplicates 00254 output (preferred)
      item = 206
! Diag09206_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
      If (.not.allocated(qcl_lsc)) Then
        Allocate ( qcl_lsc(row_length*rows,wet_model_levels) )
      End If

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(stashwork(si(item,sect,im_index)),qcl,        &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

! Copy qcl to qcl_lsc for ac_diagnostics_mod

      do k = 1,wet_model_levels
        do j = 1,rows
          do i = 1,row_length
            ji = (j-1)*row_length+i
            qcl_lsc(ji,k) = qcl(i,j,k)
          end do
        end do
      end do

         If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Qcl Main Cloud)"
         End if
      End if  ! Diag09206_if1
!
! ----------------------------------------------------------------------
! DIAG.09207 Copy cloud frozen condensate to stashwork
! ----------------------------------------------------------------------
      If (.false.) Then
! Prevent access to diagnostic after repositioning of routine until a
! final decision is made on if/how to ingest it. Timetable UM5.1 -> 5.2
      item = 207
! Diag09207_if1:
!     If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(stashwork(si(item,sect,im_index)),qcf,        &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage="cld_ctl  : error in copydiag_3d(cloud ice)"
         End if
      End if  ! Diag09207_if1
!
! ----------------------------------------------------------------------
! Find cloud base for pre defined cloud cover threshholds
! for Diagnostics 09208 - 09215.
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. ( sf(208,sect) .OR.                       &
     &                          sf(209,sect) .OR.                       &
     &                          sf(210,sect) .OR.                       &
     &                          sf(211,sect) .OR.                       &
     &                          sf(212,sect) .OR.                       &
     &                          sf(213,sect) .OR.                       &
     &                          sf(214,sect) .OR.                       &
     &                          sf(215,sect))) THEN
! Initialise output arrays
        cloud_base(:,:,:)=RMDI
        cloud_base_level(:,:,:)=IMDI

! Set cloud base to model levels
        DO l=1,noktas                 ! Loop over threshholds
          DO k=1,wet_model_levels     ! Loop over levels
            DO j=1,rows               ! Loop over rows
              DO i=1,row_length       ! Loop over points
                kinvert = wet_model_levels + 1 - k
  ! Convert to oktas
                combined_cloud_in_oktas=combined_cloud(i,j,kinvert)&
     &                                  *roktas
  ! Calculate cloud_base_level
                IF(combined_cloud_in_oktas >  oktas(l))THEN
                  IF(cloud_base_level(i,j,l) <  0)THEN
                    cloud_base_level(i,j,l) = k
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO

! Convert level numbers to heights (M converted to Kft)
        DO l=1,noktas                 ! Loop over threshholds
          DO j=1,rows                 ! Loop over rows
            DO i=1,row_length         ! Loop over points
              IF(cloud_base_level(i,j,l) > 0)THEN
                cloud_base_height =                                     &
     &                 r_rho_levels(i,j,cloud_base_level(i,j,l))        &
     &                 - earth_radius
                cloud_base(i,j,l)=cloud_base_height*m_to_kft
              END IF
            END DO
          END DO
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09208 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 0.1 okta
! ----------------------------------------------------------------------
      item = 208
! Diag09208_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,1),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 0.1 okta)"
        END IF
      END IF  ! Diag09208_if1
!
! ----------------------------------------------------------------------
! DIAG.09209 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 1.5 okta
! ----------------------------------------------------------------------
      item = 209
! Diag09209_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,2),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 1.5 okta)"
        END IF
      END IF  ! Diag09209_if1
!
! ----------------------------------------------------------------------
! DIAG.09210 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 2.5 okta
! ----------------------------------------------------------------------
      item = 210
! Diag09210_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,3),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 2.5 okta)"
        END IF
      END IF  ! Diag09210_if1
!
! ----------------------------------------------------------------------
! DIAG.09211 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 3.5 okta
! ----------------------------------------------------------------------
      item = 211
! Diag09211_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,4),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 3.5 okta)"
        END IF
      END IF  ! Diag09211_if1
!
! ----------------------------------------------------------------------
! DIAG.09212 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 4.5 okta
! ----------------------------------------------------------------------
      item = 212
! Diag09212_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,5),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 4.5 okta)"
        END IF
      END IF  ! Diag09212_if1
!
! ----------------------------------------------------------------------
! DIAG.09213 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 5.5 okta
! ----------------------------------------------------------------------
      item = 213
! Diag09213_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,6),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 5.5 okta)"
        END IF
      END IF  ! Diag09213_if1
!
! ----------------------------------------------------------------------
! DIAG.09214 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 6.5 okta
! ----------------------------------------------------------------------
      item = 214
! Diag09214_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,7),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 6.5 okta)"
        END IF
      END IF  ! Diag09214_if1
!
! ----------------------------------------------------------------------
! DIAG.09215 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 7.9 okta
! ----------------------------------------------------------------------
      item = 215
! Diag09215_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,8),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 7.9 okta)"
        END IF
      END IF  ! Diag09215_if1
!
! ***** Code matches definition in COMMON RADIATION Section 70. ****
      nclds = MIN(cloud_levels, wet_model_levels)
!
! ----------------------------------------------------------------------
! DIAG.09216 Copy Total Cloud Amount RANDOM Overlap to stashwork
! ----------------------------------------------------------------------
      item = 216
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
! Diag09216_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: r2_calc_total_cloud_cover
        CALL r2_calc_total_cloud_cover(                                   &
               row_length*rows, wet_model_levels, nclds                   &
             , ip_cloud_mix_random,combined_cloud(1,1,1),work2d_1         &
             , row_length*rows, wet_model_levels                          &
              )
!
! DEPENDS ON: copydiag
        CALL copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
              row_length,rows,0,0,0,0, at_extremity,                    &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         If (icode  >   0) then
            cmessage="cld_ctl  : error in copydiag(total cloud Random)"
         Endif
      Endif  ! Diag09216_if1
!
! ----------------------------------------------------------------------
! DIAG.09217 Copy Total Cloud Amount MAX/RANDOM Overlap to stashwork
! ----------------------------------------------------------------------
      item = 217
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
! Diag09217_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
!
! DEPENDS ON: r2_calc_total_cloud_cover
        CALL r2_calc_total_cloud_cover(                                 &
               row_length*rows, wet_model_levels, nclds,                &
               ip_cloud_mix_max, combined_cloud(1,1,1), work2d_1,       &
               row_length*rows, wet_model_levels                        &
              )

        IF (l_filter_cloud) THEN
!
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL r2_calc_total_cloud_cover(                               &
               row_length*rows, wet_model_levels, nclds,                &
               ip_cloud_mix_max,original_combined_cloud(1,1,1),work2d_2,&
               row_length*rows, wet_model_levels                        &
              )
          ! If doing filtering, calculate total cloud amount with filtered
          ! cloud field, then re-calculate it with unfiltered data.
          ! Then ensure new total cloud amount is not any bigger than
          ! original. This could happen due to formulation of max/random
          ! overlap if a chunk gets taken out of a continuous cloud deck
          ! and the previously max overlapped layer are now randomly overlapped.
          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end
              work2d_1(i,j)=MIN(work2d_1(i,j),work2d_2(i,j))
            END DO
          END DO

        END IF
!
! DEPENDS ON: copydiag
        CALL copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
              row_length,rows,0,0,0,0, at_extremity,                    &
              atmos_im,sect,item,                                       &
              icode,cmessage)

        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(total cloud Max/Rand)"
        Endif
      Endif  ! Diag09217_if1
! ----------------------------------------------------------------------
! Find Cloud Fraction in air &lt; STR_CEIL, Low Cloud Base and Top
! for Diagnostics 09218 and/or 09219 and/or 09220
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. ( sf(218,sect) .OR.                       &
     &                          sf(219,sect) .OR.                       &
     &                          sf(220,sect))) THEN
        DO j=1,rows               ! Loop over rows
          DO i=1,row_length       ! Loop over points
            pt=0.0
            ft=0.0
            work2d_1(i,j)=RMDI    ! Initialise work arrays
            work2d_2(i,j)=RMDI
            work2d_3(i,j)=RMDI
            h_asln=r_theta_levels(i,j,0) - earth_radius
            pu=p_star(i,j)
            DO k=1,min(wet_model_levels,model_levels-1)
              h_asl=h_asln
              h_asln=r_rho_levels(i,j,k+1) - earth_radius

! Check if have not already found low cloud base
              IF(work2d_2(i,j) == RMDI) THEN
! IF not, and cloud cover in this layer &gt; threshhold
                IF(area_cloud_fraction(i,j,k) >= cloud_threshold) THEN
! THEN call the bottom of this layer the base of low cloud
                  work2d_2(i,j) = h_asl / ft2m
                END IF
              END IF

! Check if already found low cloud base but not top
              IF((work2d_2(i,j) /= RMDI).AND.(work2d_3(i,j) == RMDI))   &
     &          THEN
! IF not, and cloud cover in this layer &lt; threshhold
                IF(area_cloud_fraction(i,j,k) <= cloud_threshold) THEN
! THEN call the bottom of this layer the top of low cloud
                  work2d_3(i,j) = h_asl / ft2m
                END IF
              END IF

! IF bottom of layer is below low cloud ceiling (1000ft)
               IF(h_asl  <   str_ceilm) THEN
! Calculate top and bottom layer pressures
                 pl = pu
                 pu = p(i,j,k+1)
! And accumulate presure thickness and pressure weighted cloud amount
                 dp = pu - pl
! IF whole layer below ceiling, simply accumulate whole layer.
                 IF(h_asln  <   str_ceilm) THEN
                   pt = pt + dp
                   ft = ft + dp * area_cloud_fraction(i,j,k)
                 ELSE
! Otherwise height interpolate
                  fr = (str_ceilm - h_asl ) / (h_asln - h_asl )
                  pt = pt + dp * fr
                  ft = ft + dp * area_cloud_fraction(i,j,k) * fr
! And set results
                  work2d_1(i,j) = ft / pt

                 END IF
               END IF
             END DO               ! k over max
           END DO                 ! END loop over points
         END DO                   ! END loop over rows
       END IF
!
! ----------------------------------------------------------------------
! DIAG.09218 Copy Cloud fraction below 1000ft ASL to stashwork.
! ----------------------------------------------------------------------
      item = 218
! Diag09218_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_1,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     & "cld_ctl  : error in copydiag(Cloud fraction below 1000ft ASL)"
        END IF
      END IF  ! Diag09218_if1
!
! ----------------------------------------------------------------------
! DIAG.09219 Copy Low Cloud Base to stashwork.
! ----------------------------------------------------------------------
      item = 219
! Diag09219_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_2,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Low cloud base)"
        END IF
      END IF  ! Diag09219_if1
!
! ----------------------------------------------------------------------
! DIAG.09220 Copy Low Cloud Top to stashwork.
! ----------------------------------------------------------------------
      item = 220
! Diag09220_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_3,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Low cloud top)"
        END IF
      END IF  ! Diag09220_if1
!
! ----------------------------------------------------------------------
! Calculate wet bulb temperature for Diagnostics 09221 and/or 09222
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. ( sf(221,sect) .OR.                       &
     &                          sf(222,sect))) THEN

        DO k = 1, wet_model_levels

! Initialize 2d work arrays
          DO j = 1, rows
            DO i = 1, row_length
              work2d_1(i,j) = t(i,j,k)
              work2d_2(i,j) = q(i,j,k)
              work2d_3(i,j) = p_theta_levels(i,j,k)
            END DO
          END DO

! Calculate tw
          Call ThetaW(row_length*rows,                                  &
     &                work2d_1, work2d_2, work2d_3, l_potential,        &
     &                work3d_1(1,1,k))
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09222 Copy Wet Bulb Temperature to stashwork.
! ----------------------------------------------------------------------
      item = 222
! Diag09222_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d(stashwork(si(item,sect,im_index)), work3d_1,   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag_3d(Wet Bulb Temperature)"
        END IF
      END IF  ! Diag09222_if1
!
! ----------------------------------------------------------------------
! Calculate wet bulb freezing level for Diagnostic 09221
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND.  sf(221,sect)) THEN
        DO j = 1, rows                 ! Loop over rows
          DO i = 1, row_length         ! Loop over points
            DO k = 1, wet_model_levels ! Loop over all wet-levels
              IF (work3d_1(i,j,k)  /=  rmdi) THEN
                IF (work3d_1(i,j,k)  <=  zerodegc) THEN
                  IF ( k  ==  1) THEN
                    work2d_1(i,j)=r_theta_levels(i,j,k) -               &
     &                            earth_radius
                  ELSE
                    frac = (zerodegc - work3d_1(i,j,k-1))/              &
     &                     (work3d_1(i,j,k)-work3d_1(i,j,k-1))
                    work2d_1(i,j)=r_theta_levels(i,j,k)*frac +          &
     &                            r_theta_levels(i,j,k-1)*(1.0-frac)    &
     &                            - earth_radius
                  END IF
                  EXIT
                END IF
              ELSE
                work2d_1(i,j) = rmdi
              END IF
            END DO
          END DO
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09221 Copy Wet Bulb Freezing Level to stashwork.
! ----------------------------------------------------------------------
      item = 221
! Diag09221_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_1,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Wet Bulb Freezing Level)"
        END IF
      END IF  ! Diag09221_if1

!
! ----------------------------------------------------------------------
! Calculate Total Cloud Top Height for Diagnostic 09223
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. sf(223,sect)) THEN
        DO j = 1, rows                   ! Loop over rows
          DO i = 1, row_length           ! Loop over points
            cld_top_lev = IMDI
            DO k = 1, wet_model_levels   ! Loop over all wet-levels
              IF (combined_cloud(i,j,k) > thresh) THEN
                cld_top_lev = wet_model_levels + 1 - k
              END IF
              IF (cld_top_lev  >   0) EXIT
            END DO
            IF (cld_top_lev  >   0) THEN
              work2d_1(i,j) = r_theta_levels(i,j,cld_top_lev)           &
     &                        - earth_radius
! Convert to Kiloft
              work2d_1(i,j) = work2d_1(i,j) * m_to_kft
            ELSE
              work2d_1(i,j) = RMDI
            END IF
          END DO
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09223 Copy Total Cloud Top Height to stashwork.
! ----------------------------------------------------------------------
      item = 223
! Diag09223_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_1,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Total cloud top height)"
        END IF
      END IF  ! Diag09223_if1
!
! ----------------------------------------------------------------------
! DIAG.09229 Copy Relative Humidity to stashwork
! ----------------------------------------------------------------------
      item = 229
! Diag09229_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Calculate Qsat then work out Relative Humidity.
! Diag09229_do1:
        Do k = 1, wet_model_levels
! DEPENDS ON: qsat_mix
          Call qsat_mix(work2d_1,T(1,1,K),p_theta_levels(1,1,k),        &
     &              row_length*rows,l_mixing_ratio)
!
          Do j = 1, rows
            Do i = 1, row_length
              ! relative humidity in per cent
              work3d_1(i, j, k) = q(i, j, k) / work2d_1(i, j) *100.0

!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
              If ( work3d_1(i, j, k) < 0.0) Then
                  work3d_1(i, j, k) = 0.0
              Endif

            End Do
          End Do
!
        End Do  ! Diag09229_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(RH Main Cloud)"
        End if
      End if  ! Diag09229_if1
!
! ----------------------------------------------------------------------
! DIAG.09226 Copy Layer Cloud Frequency to stashwork.
! ----------------------------------------------------------------------
      item = 226
! Diag09226_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
!
! Freq_k_do1:
        DO k = 1, wet_model_levels
          DO j = 1, rows
            DO i = 1, row_length
              IF (area_cloud_fraction(i, j, k)  <=  0.) THEN 
                work3d_1(i, j, k) = 0.
              ELSE
                work3d_1(i, j, k) = 1.
              END IF
            END DO
          END DO
        END DO  ! Freq_k_do1
!
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),work3d_1,   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage="cld_ctl  : error in copydiag_3d(Layer Cloud Freq)"
        END IF
      END IF  ! Diag09226_if1
!
! ----------------------------------------------------------------------
! DIAG.09228 Copy Critical Relative Humidity to stashwork.
! ----------------------------------------------------------------------
      item = 228
! Diag09228_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!     STASH will only permit this diagnostic to be chosen if 3D RHcrit
!     diagnosis is active. So RHCPT should be (row_length,rows,wet_ml).
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),RHCPT,      &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Critical RH)"
        End if
      End if  ! Diag09228_if1
!
! ----------------------------------------------------------------------
! DIAG.09230 Calculate visibility on model levels and copy to stashwork
! ----------------------------------------------------------------------
      item = 230
! Note that the following calculates the visibility for all model levels
! and then extracts just those levels that have been requested. If the
! diagnostic is normally called for a small subset of levels it will be
! more efficient to introduce a call to set_levels_list and calculate
! only those levels that are needed.
! Diag09230_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then

        do k=1,wet_model_levels

          ! visibility is calculated using prescribed critical RH
          ! (not diagnostic RH)
! DEPENDS ON: visbty
          Call Visbty(                                                  &
     &     p_theta_levels(1,1,k), T(1,1,k), q(1,1,k)                    &
                                                           !INPUT
     &     ,qcl(1,1,k), qcf(1,1,k)                                      &
                                                           !INPUT
     &     ,Aerosol(1:row_length,1:rows,k)                              &
                                                           !INPUT
     &     ,vis_probability, RHcrit(k), L_murk                          &
                                                           !INPUT
     &     ,row_length * rows                                           &
                                                           !INPUT
     &     ,work3d_1(1,1,k))                               !OUTPUT

        End do ! k   wet_model_levels

! DEPENDS ON: copydiag_3d
        Call copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)


      End if  ! Diag09230_if1

! ----------------------------------------------------------------------
! DIAG.09231 Copy Combined Cloud fraction to stashwork
! ----------------------------------------------------------------------
      item = 231
! Diag09231_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Diag09231_do1:
        Do k = 1, wet_model_levels
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
!     Combined_cloud is calculated this way but re-inverted for STASH.
!
          kinvert = wet_model_levels+1-k
          Do j = 1, rows
            Do i = 1, row_length
              work3d_1(i, j, k) = combined_cloud(i,j,kinvert)
            End Do
          End Do
        End Do  ! Diag09231_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)), work3d_1,  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Combined Cld On Lv)"
        End if
      End if  ! Diag09231_if1
! ----------------------------------------------------------------------
! DIAGS.09232-00234 Filter out cloud above ceilometer range
! N.B. If l_filter_cloud=.TRUE. then these diagnostics will operate on the
! cloud fields that have ALSO been filtered to remove optically thin clouds.
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND.                                           &
           ( sf(232,sect) .OR. sf(233,sect) .OR. sf(234,sect) )) THEN

        ALLOCATE ( comb_ceilometer_cloud( (qdims%i_end-qdims%i_start)+1,&
                                          (qdims%j_end-qdims%j_start)+1,&
                                           qdims%k_end ) )

! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
!     Combined_cloud is calculated this way but re-inverted for STASH.
!
        DO k = 1, wet_model_levels
          kinvert = wet_model_levels+1-k
          DO j = 1, rows
            DO i = 1, row_length

! Set cloud on model levels above ceilometer range to zero
              IF ( (r_theta_levels(i,j,kinvert)-r_theta_levels(i,j,0))  &
                   > ceilometer_range ) THEN
                comb_ceilometer_cloud(i, j, k) = 0.0
              ELSE
                comb_ceilometer_cloud(i, j, k) = combined_cloud(i,j,k)
              END IF
            END DO
          END DO
        END DO

! ----------------------------------------------------------------------
! DIAG.09232 Copy Ceilometer Filtered Cloud Amount RANDOM Overlap 
!                                                    to stashwork
! ----------------------------------------------------------------------
        item = 232
        IF (icode <= 0 .AND. sf(item,sect)) THEN
!
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                               &
               row_length*rows, wet_model_levels, nclds,                &
               IP_CLOUD_MIX_RANDOM, comb_ceilometer_cloud(1,1,1),       &
               work2d_1,row_length*rows, wet_model_levels)
!
! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)),work2d_1,     &
                row_length,rows,0,0,0,0, at_extremity,                  &
                atmos_im,sect,item,                                     &
                icode,cmessage)

          IF (icode  >   0) THEN
            cmessage =                                                  &
                 "cld_ctl  : error in copydiag(total ceil cld Random)"
          END IF
          
        END IF

! ----------------------------------------------------------------------
! DIAG.09233 Copy Ceilometer Filtered Cloud Amount MAX/RANDOM Overlap 
!                                                        to stashwork
! ----------------------------------------------------------------------
        item = 233
        IF (icode  <=  0  .AND.  sf(item,sect)) THEN
!
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                               &
               row_length*rows, wet_model_levels, nclds,                &
               IP_CLOUD_MIX_MAX, comb_ceilometer_cloud(1,1,1),work2d_1, &
               row_length*rows, wet_model_levels)
!
! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)),work2d_1,     &
               row_length,rows,0,0,0,0, at_extremity,                   &
               atmos_im,sect,item,                                      &
               icode,cmessage)

          IF (icode  >   0) THEN
            cmessage =                                                  &
               "cld_ctl  : error in copydiag(total ceil cloud Max/Rand)"
        END IF
      END IF ! Diag09217_if1
! ----------------------------------------------------------------------
! DIAG.09234 Copy Ceilometer Filtered Cloud fraction to stashwork
! ----------------------------------------------------------------------
        item = 234
        IF (icode <= 0 .AND. sf(item,sect)) THEN

          DO k = 1, wet_model_levels
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
! comb_ceilometer_cloud is calculated this way but re-inverted for STASH.
!
            kinvert = wet_model_levels+1-k
            DO j = 1, rows
              DO i = 1, row_length
                work3d_1(i, j, k) = comb_ceilometer_cloud(i,j,kinvert)
              END DO
            END DO
          END DO

          CALL copydiag_3d (stashwork(si(item,sect,im_index)), work3d_1,&
              row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
!
          IF (icode > 0) THEN
            cmessage="cld_ctl  : error in copydiag_3d(Ceilometer Cld On Lv)"
          END IF
        END IF

        DEALLOCATE (comb_ceilometer_cloud)

      END IF

      IF (l_filter_cloud .AND. l_combi_cld) THEN
        DEALLOCATE (original_combined_cloud)
      END IF
!
!=======================================================================
! START AIRCRAFT ICING ALGORITHMS
!=======================================================================
!
! All icing products are conditional upon being In Right Temperature
! Range (IRTR). Additional tests are done on top of this.
! It is assumed that if you request FOO of an icing index above a 
! threshold, that you have also request that icing index!
! ----------------------------------------------------------------------
! DIAG.09310: Relative Humidity (IRTR), 0 otherwise.
! II(RH)
! ----------------------------------------------------------------------

      IF (icode  <=  0  .AND. ( sf(310,sect).OR.sf(311,sect).OR.           &
                                sf(312,sect).OR.sf(313,sect).OR.           &
                                sf(314,sect).OR.sf(315,sect).OR.           &
                                sf(316,sect).OR.sf(317,sect).OR.           &
                                sf(318,sect).OR.sf(319,sect).OR.           &
                                sf(391,sect)                     ) ) THEN

        ALLOCATE ( icing3d( (qdims%i_end-qdims%i_start)+1,                 &
                            (qdims%j_end-qdims%j_start)+1,                 &
                             qdims%k_end ) )
        icing3d = 0.0

        DO k = 1, qdims%k_end
! DEPENDS ON: qsat_mix
          CALL qsat_mix(work2d_1,t(1,1,k),p_theta_levels(1,1,k),           &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               l_mixing_ratio)

          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end
              IF ((t(i,j,k) < iitmax) .AND. (t(i,j,k) > iitmin)) THEN
                icing3d(i,j,k) = q(i,j,k) / work2d_1(i, j)           
              END IF
            END DO
          END DO

        END DO

        item = 310
        v_item=391

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work3d_1 = icing3d

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,    &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
              qdims%k_end,0,0,0,0, at_extremity,                          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,         &
              stash_levels,num_stash_levels+1,                            &
              atmos_im,sect,item,                                         &
              icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag_3d(IcingIndex)"
          END IF

        END IF

        ref_item=item

        DO my_index=1,9

        item = ref_item+my_index

          IF (icode  <=  0  .AND.  sf(item,sect)) THEN
            work3d_1 = 0.0

            DO k = 1, qdims%k_end
              DO j = qdims%j_start,qdims%j_end
                DO i = qdims%i_start,qdims%i_end
                  IF (icing3d(i,j,k) >= REAL(my_index)*0.1) THEN
                    work3d_1(i,j,k) = 1.0
                  END IF
                END DO
              END DO
            END DO

! DEPENDS ON: copydiag_3d
            CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,  &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
              qdims%k_end,0,0,0,0, at_extremity,                          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,         &
              stash_levels,num_stash_levels+1,                            &
              atmos_im,sect,item,                                         &
              icode,cmessage)

            IF (icode >  0) THEN
              cmessage="cld_ctl  : error in copydiag_3d(Icing loop)"
            END IF

          END IF

        END DO ! my_index=1,9 

        item = v_item

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work2d_1=0.0

! Calculate a vertically-"integrated" icing index using maximum/random overlap
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                                  &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end, nclds,                                         &
               ip_cloud_mix_max, icing3d(1,1,1), work2d_1,                 &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end)

! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)), work2d_1,       &
               (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
               0,0,0,0, at_extremity,                                      &
               atmos_im,sect,item,                                         &
               icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag(Max/Rand Icing)"
          END IF

        END IF

        DEALLOCATE ( icing3d )

      END IF

! ----------------------------------------------------------------------
! DIAG.09320: Relative Humidity if bulk cloud present (IRTR), 
!             0 otherwise.
! II(RHCldPres)
! ----------------------------------------------------------------------

      IF (icode  <=  0  .AND. ( sf(320,sect).OR.sf(321,sect).OR.           &
                                sf(322,sect).OR.sf(323,sect).OR.           &
                                sf(324,sect).OR.sf(325,sect).OR.           &
                                sf(326,sect).OR.sf(327,sect).OR.           &
                                sf(328,sect).OR.sf(329,sect).OR.           &
                                sf(392,sect)                     ) ) THEN

        ALLOCATE ( icing3d( (qdims%i_end-qdims%i_start)+1,                 &
                            (qdims%j_end-qdims%j_start)+1,                 &
                             qdims%k_end ) )
        icing3d = 0.0

        DO k = 1, qdims%k_end
! DEPENDS ON: qsat_mix
          CALL qsat_mix(work2d_1,t(1,1,k),p_theta_levels(1,1,k),           &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               l_mixing_ratio)

          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end
              IF ((t(i,j,k) < iitmax) .AND. (t(i,j,k) > iitmin)) THEN
                IF (bulk_cloud_fraction(i,j,k) > 0.0) THEN
                   icing3d(i,j,k) = q(i,j,k) / work2d_1(i, j)           
                END IF
              END IF
            END DO
          END DO
        END DO

        item = 320
        v_item=392

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work3d_1 = icing3d

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,    &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
              qdims%k_end,0,0,0,0, at_extremity,                          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,         &
              stash_levels,num_stash_levels+1,                            &
              atmos_im,sect,item,                                         &
              icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag_3d(IcingIndex)"
          END IF

        END IF

        ref_item=item

        DO my_index=1,9

        item = ref_item+my_index

          IF (icode  <=  0  .AND.  sf(item,sect)) THEN
            work3d_1 = 0.0

            DO k = 1, qdims%k_end
              DO j = qdims%j_start,qdims%j_end
                DO i = qdims%i_start,qdims%i_end
                  IF (icing3d(i,j,k) >= REAL(my_index)*0.1) THEN
                    work3d_1(i,j,k) = 1.0
                  END IF
                END DO
              END DO
            END DO

! DEPENDS ON: copydiag_3d
            CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,  &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
              qdims%k_end,0,0,0,0, at_extremity,                          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,         &
              stash_levels,num_stash_levels+1,                            &
              atmos_im,sect,item,                                         &
              icode,cmessage)

            IF (icode >  0) THEN
              cmessage="cld_ctl  : error in copydiag_3d(Icing loop)"
            END IF

          END IF

        END DO ! my_index=1,9 

        item = v_item

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work2d_1=0.0

! Calculate a vertically-"integrated" icing index using maximum/random overlap
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                                  &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end, nclds,                                         &
               ip_cloud_mix_max, icing3d(1,1,1), work2d_1,                 &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end)

! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)), work2d_1,       &
               (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
               0,0,0,0, at_extremity,                                      &
               atmos_im,sect,item,                                         &
               icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag(Max/Rand Icing)"
          END IF

        END IF

        DEALLOCATE ( icing3d )

      END IF

! ----------------------------------------------------------------------
! DIAG.09330: Bulk cloud fraction (IRTR)
! II(CF)
! ----------------------------------------------------------------------

      IF (icode  <=  0  .AND. ( sf(330,sect).OR.sf(331,sect).OR.          &
                                sf(332,sect).OR.sf(333,sect).OR.          &
                                sf(334,sect).OR.sf(335,sect).OR.          &
                                sf(336,sect).OR.sf(337,sect).OR.          &
                                sf(338,sect).OR.sf(339,sect).OR.          &
                                sf(393,sect)                     ) ) THEN

        ALLOCATE ( icing3d( (qdims%i_end-qdims%i_start)+1,                &
                            (qdims%j_end-qdims%j_start)+1,                &
                             qdims%k_end ) )
        icing3d = 0.0

        DO k = 1, qdims%k_end
          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end
              IF ((t(i,j,k) < iitmax) .AND. (t(i,j,k) > iitmin)) THEN
                icing3d(i,j,k) = bulk_cloud_fraction(i,j,k)
              END IF
            END DO
          END DO
        END DO

        item = 330
        v_item=393

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work3d_1 = icing3d

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,    &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
              qdims%k_end,0,0,0,0, at_extremity,                          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,         &
              stash_levels,num_stash_levels+1,                            &
              atmos_im,sect,item,                                         &
              icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag_3d(IcingIndex)"
          END IF

        END IF

        ref_item=item

        DO my_index=1,9

        item = ref_item+my_index

          IF (icode  <=  0  .AND.  sf(item,sect)) THEN
            work3d_1 = 0.0

            DO k = 1, qdims%k_end
              DO j = qdims%j_start,qdims%j_end
                DO i = qdims%i_start,qdims%i_end
                  IF (icing3d(i,j,k) >= REAL(my_index)*0.1) THEN
                    work3d_1(i,j,k) = 1.0
                  END IF
                END DO
              END DO
            END DO

! DEPENDS ON: copydiag_3d
            CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,  &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
              qdims%k_end,0,0,0,0, at_extremity,                          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,         &
              stash_levels,num_stash_levels+1,                            &
              atmos_im,sect,item,                                         &
              icode,cmessage)

            IF (icode >  0) THEN
              cmessage="cld_ctl  : error in copydiag_3d(Icing loop)"
            END IF

          END IF

        END DO ! my_index=1,9 

        item = v_item

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work2d_1=0.0

! Calculate a vertically-"integrated" icing index using maximum/random overlap
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                                  &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end, nclds,                                         &
               ip_cloud_mix_max, icing3d(1,1,1), work2d_1,                 &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end)

! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)), work2d_1,       &
               (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
               0,0,0,0, at_extremity,                                      &
               atmos_im,sect,item,                                         &
               icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag(Max/Rand Icing)"
          END IF

        END IF

        DEALLOCATE ( icing3d )

      END IF

! ----------------------------------------------------------------------
! DIAG.09340: Liquid cloud fraction (IRTR)
! II(LCF)
! ----------------------------------------------------------------------

      IF (icode  <=  0  .AND. ( sf(340,sect).OR.sf(341,sect).OR.           &
                                sf(342,sect).OR.sf(343,sect).OR.           &
                                sf(344,sect).OR.sf(345,sect).OR.           &
                                sf(346,sect).OR.sf(347,sect).OR.           &
                                sf(348,sect).OR.sf(349,sect).OR.           &
                                sf(394,sect)                     ) ) THEN

        ALLOCATE ( icing3d( (qdims%i_end-qdims%i_start)+1,                 &
                            (qdims%j_end-qdims%j_start)+1,                 &
                             qdims%k_end ) )
        icing3d = 0.0

        DO k = 1, qdims%k_end
          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end
              IF ((t(i,j,k) < iitmax) .AND. (t(i,j,k) > iitmin)) THEN
                icing3d(i,j,k) = MIN(MAX(cfl(i,j,k),0.0),1.0)
              END IF
            END DO
          END DO
        END DO

        item = 340
        v_item=394

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work3d_1 = icing3d

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,     &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag_3d(IcingIndex)"
          END IF

        END IF

        ref_item=item

        DO my_index=1,9

        item = ref_item+my_index

          IF (icode  <=  0  .AND.  sf(item,sect)) THEN
            work3d_1 = 0.0

            DO k = 1, qdims%k_end
              DO j = qdims%j_start,qdims%j_end
                DO i = qdims%i_start,qdims%i_end
                  IF (icing3d(i,j,k) >= REAL(my_index)*0.1) THEN
                    work3d_1(i,j,k) = 1.0
                  END IF
                END DO
              END DO
            END DO

! DEPENDS ON: copydiag_3d
            CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,   &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

            IF (icode >  0) THEN
              cmessage="cld_ctl  : error in copydiag_3d(Icing loop)"
            END IF

          END IF

        END DO ! my_index=1,9 

        item = v_item

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work2d_1=0.0

! Calculate a vertically-"integrated" icing index using maximum/random overlap
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                                  &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end, nclds,                                         &
               ip_cloud_mix_max, icing3d(1,1,1), work2d_1,                 &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end)

! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)), work2d_1,       &
               (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
               0,0,0,0, at_extremity,                                      &
               atmos_im,sect,item,                                         &
               icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag(Max/Rand Icing)"
          END IF

        END IF

        DEALLOCATE ( icing3d )

      END IF

! ----------------------------------------------------------------------
! DIAG.09350: Grid-box-mean Liquid Water Content (IRTR) kg m-3 not kg kg-1
! normalised by a threshold value such that index=0.5 if LWC=threshold and
! capped so cannot be greater than 1.0
! II(LWC)
! ----------------------------------------------------------------------

      IF (icode  <=  0  .AND. ( sf(350,sect).OR.sf(351,sect).OR.           &
                                sf(352,sect).OR.sf(353,sect).OR.           &
                                sf(354,sect).OR.sf(355,sect).OR.           &
                                sf(356,sect).OR.sf(357,sect).OR.           &
                                sf(358,sect).OR.sf(359,sect).OR.           &
                                sf(395,sect)                     ) ) THEN

        ALLOCATE ( icing3d( (qdims%i_end-qdims%i_start)+1,                 &
                            (qdims%j_end-qdims%j_start)+1,                 &
                             qdims%k_end ) )
        icing3d = 0.0

        DO k = 1, qdims%k_end-1
          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end
              IF ((t(i,j,k) < iitmax) .AND. (t(i,j,k) > iitmin)) THEN

                icing3d(i,j,k) = ( qcl(i,j,k) * rho(i,j,k) ) /             &
                                 ( 2.0 * qcl_thresh )

                icing3d(i,j,k) = MIN(icing3d(i,j,k),1.0)

              END IF
            END DO
          END DO
        END DO

        item = 350
        v_item=395

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work3d_1 = icing3d

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,     &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag_3d(IcingIndex)"
          END IF

        END IF

        ref_item=item

        DO my_index=1,9

        item = ref_item+my_index

          IF (icode  <=  0  .AND.  sf(item,sect)) THEN
            work3d_1 = 0.0

            DO k = 1, qdims%k_end
              DO j = qdims%j_start,qdims%j_end
                DO i = qdims%i_start,qdims%i_end
                  IF (icing3d(i,j,k) >= REAL(my_index)*0.1) THEN
                    work3d_1(i,j,k) = 1.0
                  END IF
                END DO
              END DO
            END DO

! DEPENDS ON: copydiag_3d
            CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,   &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

            IF (icode >  0) THEN
              cmessage="cld_ctl  : error in copydiag_3d(Icing loop)"
            END IF

          END IF

        END DO ! my_index=1,9 

        item = v_item

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work2d_1=0.0

! Calculate a vertically-"integrated" icing index using maximum/random overlap
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                                  &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end, nclds,                                         &
               ip_cloud_mix_max, icing3d(1,1,1), work2d_1,                 &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end)

! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)), work2d_1,       &
               (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
               0,0,0,0, at_extremity,                                      &
               atmos_im,sect,item,                                         &
               icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag(Max/Rand Icing)"
          END IF

        END IF

        DEALLOCATE ( icing3d )

      END IF

! ----------------------------------------------------------------------
! DIAG.09360: In-cloud-mean Liquid Water Content (IRTR) kg m-3 not kg kg-1
! normalised by a threshold value such that index=0.5 if LWC=threshold and
! capped so cannot be greater than 1.0
! II(icLWC)
! ----------------------------------------------------------------------

      IF (icode  <=  0  .AND. ( sf(360,sect).OR.sf(361,sect).OR.        &
                                sf(362,sect).OR.sf(363,sect).OR.        &
                                sf(364,sect).OR.sf(365,sect).OR.        &
                                sf(366,sect).OR.sf(367,sect).OR.        &
                                sf(368,sect).OR.sf(369,sect).OR.        &
                                sf(396,sect)                     ) ) THEN

        ALLOCATE ( icing3d( (qdims%i_end-qdims%i_start)+1,              &
                            (qdims%j_end-qdims%j_start)+1,              &
                             qdims%k_end ) )
        icing3d = 0.0

        DO k = 1, qdims%k_end-1
          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end
              IF ((t(i,j,k) < iitmax) .AND. (t(i,j,k) > iitmin)) THEN
                IF (cfl(i,j,k) > 0.01) THEN

                  icing3d(i,j,k) = ( qcl(i,j,k) * rho(i,j,k) ) /        &
                                   ( cfl(i,j,k) * 2.0 * qcl_thresh )

                  icing3d(i,j,k) = MIN(icing3d(i,j,k),1.0)

                END IF
              END IF
            END DO
          END DO
        END DO

        item = 360
        v_item=396

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work3d_1 = icing3d

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,     &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag_3d(IcingIndex)"
          END IF

        END IF

        ref_item=item

        DO my_index=1,9

        item = ref_item+my_index

          IF (icode  <=  0  .AND.  sf(item,sect)) THEN
            work3d_1 = 0.0

            DO k = 1, qdims%k_end
              DO j = qdims%j_start,qdims%j_end
                DO i = qdims%i_start,qdims%i_end
                  IF (icing3d(i,j,k) >= REAL(my_index)*0.1) THEN
                    work3d_1(i,j,k) = 1.0
                  END IF
                END DO
              END DO
            END DO

! DEPENDS ON: copydiag_3d
            CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,   &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

            IF (icode >  0) THEN
              cmessage="cld_ctl  : error in copydiag_3d(Icing loop)"
            END IF

          END IF

        END DO ! my_index=1,9 

        item = v_item

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work2d_1=0.0

! Calculate a vertically-"integrated" icing index using maximum/random overlap
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                                  &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end, nclds,                                         &
               ip_cloud_mix_max, icing3d(1,1,1), work2d_1,                 &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end)

! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)), work2d_1,       &
             (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,  &
             0,0,0,0, at_extremity,                                        &
             atmos_im,sect,item,                                           &
             icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag(Max/Rand Icing)"
          END IF

        END IF

        DEALLOCATE ( icing3d )

      END IF

! ----------------------------------------------------------------------
! DIAG.09370: Fraction Of Sky above Threshold (FrOST) 
!             liquid water content, using a log-normal distribution.
! II(FROST-LN)
! ----------------------------------------------------------------------

      IF (icode  <=  0  .AND. ( sf(370,sect).OR.sf(371,sect).OR.           &
                                sf(372,sect).OR.sf(373,sect).OR.           &
                                sf(374,sect).OR.sf(375,sect).OR.           &
                                sf(376,sect).OR.sf(377,sect).OR.           &
                                sf(378,sect).OR.sf(379,sect).OR.           &
                                sf(397,sect)                     ) ) THEN

        ALLOCATE ( icing3d( (qdims%i_end-qdims%i_start)+1,                 &
                            (qdims%j_end-qdims%j_start)+1,                 &
                             qdims%k_end ) )
        icing3d = 0.0

        DO k = 1, qdims%k_end-1
          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end

              IF ((t(i,j,k) < iitmax) .AND. (t(i,j,k) > iitmin)) THEN

                IF ( (cfl(i,j,k) > 0.01) .AND. (qcl(i,j,k) > 1.0e-10) ) THEN
                  ! Calculate mean in-cloud LWC
                  mean_qcl = ( qcl(i,j,k) * rho(i,j,k) ) / cfl(i,j,k)  

                  ! Calculate a grid-box size in km
                  horiz_scale = SQRT ( r_theta_levels(i,j,k) * delta_lambda   &
                                     * r_theta_levels(i,j,k) * delta_phi      &
                                     * fv_cos_theta_latitude(i,j) )           &
                                     * 1.0e-3
                  ! Calculate fractional standard deviation (fsd) from 
                  ! parametrization given in Boutle et al (2013).
                  phic = ( (horiz_scale*cfl(i,j,k)) ** 0.3333 )               &
                       * ((((0.06*horiz_scale*cfl(i,j,k))**1.5)+1.0)**(-0.17))

                  IF (cfl(i,j,k) < 1.0) THEN
                    fsd = (0.45-(0.25*cfl(i,j,k)))*phic
                  ELSE
                    fsd = 0.11 * phic
                  END IF

                  ! Rescale by sqrt(2.0) to account for 1d to 2d
                  fsd=fsd*1.414

                  ! Find variance of LWC PDF from the fractional 
                  ! standard deviation.
                  variance_qcl = (mean_qcl*fsd)**2

                  ! Assume in-cloud LWC is given by a log-normal PDF.
                  tmp   = (variance_qcl/EXP(2.0*LOG(mean_qcl)))+1.0
                  mu    = LOG(mean_qcl)-(0.5*LOG(tmp))
                  sigma = SQRT((2.0*LOG(mean_qcl))-(2.0*mu))
                  ! Calculate Fraction Of Cloud Above Threshold (FOCAT).
                  ! This is given by the integral of the PDF from 
                  ! qcl_thresh to infinity. Or equivalently by 1-CDF 
                  ! where the cumulative distribution function 
                  ! CDF is the integral from 0 to qcl_thresh.

                  ! For a log-normal PDF, cumulative distribution 
                  ! function CDF is related to ERF.
                  erf_in(1) = (LOG(qcl_thresh)-mu)/SQRT(2.0*sigma*sigma)

! DEPENDS ON: mym_errfunc
                  CALL mym_errfunc(1, erf_in(1), erf_out(1))

                  ! FOCAT = 1 - CDF = 1 - (0.5+0.5*erf_out)
                  focat = 0.5-(0.5*erf_out(1))
                  focat=MIN(MAX(focat,0.0),1.0)

                  ! Fraction Of Sky above Threshold is then found from 
                  ! product of Fraction Of Cloud Above Threshold (FOCAT)
                  ! and Fraction of sky which is cloudy (CFL).
                  icing3d(i,j,k) = focat*cfl(i,j,k)

                END IF
              END IF

            END DO
          END DO
        END DO

        item = 370
        v_item=397

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work3d_1 = icing3d

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,     &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag_3d(IcingIndex)"
          END IF

        END IF

        ref_item=item

        DO my_index=1,9

        item = ref_item+my_index

          IF (icode  <=  0  .AND.  sf(item,sect)) THEN
            work3d_1 = 0.0

            DO k = 1, qdims%k_end
              DO j = qdims%j_start,qdims%j_end
                DO i = qdims%i_start,qdims%i_end
                  IF (icing3d(i,j,k) >= REAL(my_index)*0.1) THEN
                    work3d_1(i,j,k) = 1.0
                  END IF
                END DO
              END DO
            END DO

! DEPENDS ON: copydiag_3d
            CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,   &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

            IF (icode >  0) THEN
              cmessage="cld_ctl  : error in copydiag_3d(Icing loop)"
            END IF

          END IF

        END DO ! my_index=1,9 

        item = v_item

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work2d_1=0.0

! Calculate a vertically-"integrated" icing index using maximum/random overlap
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                                  &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end, nclds,                                         &
               ip_cloud_mix_max, icing3d(1,1,1), work2d_1,                 &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end)

! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)), work2d_1,       &
               (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
               0,0,0,0, at_extremity,                                      &
               atmos_im,sect,item,                                         &
               icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag(Max/Rand Icing)"
          END IF

        END IF

        DEALLOCATE ( icing3d )

      END IF

! ----------------------------------------------------------------------
! DIAG.09380: Fraction Of Sky above Threshold (FrOST) 
!             liquid water content, using a gamma distribution.
! II(FROST-Gamma)
! ----------------------------------------------------------------------

      IF (icode  <=  0  .AND. ( sf(380,sect).OR.sf(381,sect).OR.        &
                                sf(382,sect).OR.sf(383,sect).OR.        &
                                sf(384,sect).OR.sf(385,sect).OR.        &
                                sf(386,sect).OR.sf(387,sect).OR.        &
                                sf(388,sect).OR.sf(389,sect).OR.        &
                                sf(398,sect)                     ) ) THEN

        ALLOCATE ( icing3d( (qdims%i_end-qdims%i_start)+1,              &
                            (qdims%j_end-qdims%j_start)+1,              &
                             qdims%k_end ) )
        icing3d = 0.0

        DO k = 1, qdims%k_end-1
          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end

              IF ((t(i,j,k) < iitmax) .AND. (t(i,j,k) > iitmin)) THEN

                IF ( (cfl(i,j,k) > 0.01) .AND. (qcl(i,j,k) > 1.0e-10) ) THEN
                  ! Calculate mean in-cloud LWC
                  mean_qcl = ( qcl(i,j,k) * rho(i,j,k) ) / cfl(i,j,k)

                  ! Calculate a grid-box size in km
                  horiz_scale = SQRT ( r_theta_levels(i,j,k) * delta_lambda   &
                                     * r_theta_levels(i,j,k) * delta_phi      &
                                     * fv_cos_theta_latitude(i,j) )           &
                                     * 1.0e-3
                  ! Calculate fractional standard deviation (fsd) from 
                  ! parametrization given in Boutle et al (2013).
                  phic = ( (horiz_scale*cfl(i,j,k)) ** 0.3333 )               &
                       * ((((0.06*horiz_scale*cfl(i,j,k))**1.5)+1.0)**(-0.17))

                  IF (cfl(i,j,k) < 1.0) THEN
                    fsd = (0.45-(0.25*cfl(i,j,k)))*phic
                  ELSE
                    fsd = 0.11 * phic
                  END IF

                  ! Rescale by sqrt(2.0) to account for 1d to 2d
                  fsd=fsd*1.414

                  ! Find variance of LWC PDF from the fractional 
                  ! standard deviation.
                  variance_qcl = (mean_qcl*fsd)**2
                  ! Assume in-cloud LWC is given by a gamma PDF.
                  thet = variance_qcl/mean_qcl
                  kay  = mean_qcl/thet

                  ! Calculate Fraction Of Cloud Above Threshold (FOCAT).
                  ! This is given by the integral of the PDF from 
                  ! qcl_thresh to infinity. Or equivalently by 1-CDF 
                  ! where the cumulative distribution function 
                  ! CDF is the integral from 0 to qcl_thresh.

                  ! For a gamma PDF, cumulative distribution 
                  ! function CDF is related to the complete and 
                  ! incomplete gamma function (Big Gamma and 
                  ! LittleGamma respectively).
                  CALL gammaf(kay, BigGamma)

                  CALL gammaf_lower_incomp(kay,qcl_thresh/thet,LittleGamma,20)
                  focat = 1.0 - ( LittleGamma/BigGamma )
                  focat=MIN(MAX(focat,0.0),1.0)

                  ! Fraction Of Sky above Threshold is then found from 
                  ! product of Fraction Of Cloud Above Threshold (FOCAT) 
                  ! and Fraction of sky which is cloudy (CFL).
                  icing3d(i,j,k)=focat*cfl(i,j,k)

                END IF
              END IF

            END DO
          END DO
        END DO

        item = 380
        v_item=398

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work3d_1 = icing3d

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,     &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag_3d(IcingIndex)"
          END IF

        END IF

        ref_item=item

        DO my_index=1,9

        item = ref_item+my_index

          IF (icode  <=  0  .AND.  sf(item,sect)) THEN
            work3d_1 = 0.0

            DO k = 1, qdims%k_end
              DO j = qdims%j_start,qdims%j_end
                DO i = qdims%i_start,qdims%i_end
                  IF (icing3d(i,j,k) >= REAL(my_index)*0.1) THEN
                    work3d_1(i,j,k) = 1.0
                  END IF
                END DO
              END DO
            END DO

! DEPENDS ON: copydiag_3d
            CALL copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,   &
              (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1, &
              qdims%k_end,0,0,0,0, at_extremity,                           &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                             &
              atmos_im,sect,item,                                          &
              icode,cmessage)

            IF (icode >  0) THEN
              cmessage="cld_ctl  : error in copydiag_3d(Icing loop)"
            END IF

          END IF

        END DO ! my_index=1,9 

        item = v_item

        IF (icode  <=  0  .AND.  sf(item,sect)) THEN

          work2d_1=0.0

! Calculate a vertically-"integrated" icing index using maximum/random overlap
! DEPENDS ON: r2_calc_total_cloud_cover
          CALL R2_calc_total_cloud_cover(                                  &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end, nclds,                                         &
               ip_cloud_mix_max, icing3d(1,1,1), work2d_1,                 &
               (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),&
               qdims%k_end)

! DEPENDS ON: copydiag
          CALL copydiag(stashwork(si(item,sect,im_index)), work2d_1,       &
               (qdims%i_end-qdims%i_start)+1,(qdims%j_end-qdims%j_start)+1,&
               0,0,0,0, at_extremity,                                      &
               atmos_im,sect,item,                                         &
               icode,cmessage)

          IF (icode >  0) THEN
            cmessage="cld_ctl  : error in copydiag(Max/Rand Icing)"
          END IF

        END IF

        DEALLOCATE ( icing3d )

      END IF
!
!=======================================================================
! END AIRCRAFT ICING ALGORITHMS
!=======================================================================
!
! ----------------------------------------------------------------------
      If (icode  /=  0) Then

        Call Ereport(RoutineName,icode,Cmessage)
      End if
!
      IF (lhook) CALL dr_hook('DIAGNOSTICS_LSCLD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE diagnostics_lscld
