! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! ALLOCATE and initialise arrays for convective diagnostics
!
SUBROUTINE cv_alloc_diag_array( row_length, rows, l_calc_dxek, l_cosp,        &
                               l_dust, l_sulpc_so2, l_sulpc_nh3, l_soot,      &
                               l_biomass, l_ocff, l_nitrate,                  &
! Stash info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                               ntml,ntpar,                                    &
                               theta_inc, q_inc, qcl_inc, qcf_inc,            &
                               cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,     &
                               conv_rain, conv_snow)


! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                          &
   udims, vdims, wdims, tdims, pdims, qdims

USE cv_run_mod,  ONLY:                                                    &
   n_conv_calls, l_mom

USE cv_diagnostic_array_mod , ONLY:                                     &
        precip_deep, precip_shall ,precip_mid ,precip_cong ,cape_out    &
       ,deep_ind ,shallow_ind ,congestus_ind ,congestus_ind2 ,mid_ind   &
       ,ntml_diag ,ntpar_diag ,freeze_diag ,kterm_diag                  &
       ,wstar_up_diag ,wstar_dn_diag ,mb1_diag ,mb2_diag ,wqt_cb        &
       ,wthetal_cb ,wqt_inv ,wthetal_inv ,sh_top ,sh_base ,cg_top       &
       ,cg_base ,cg_term ,cape_ts_diag ,ind_cape_reduced_diag           &
       ,conscav_dust ,conscav_so4ait ,conscav_so4acc ,conscav_so4dis    &
       ,conscav_agedsoot ,conscav_agedbmass ,conscav_agedocff           &
       ,conscav_nitracc ,conscav_nitrdiss ,conwash_so2 ,conwash_nh3     &
       ,uw_dp, vw_dp ,uw_shall ,vw_shall ,uw_mid ,vw_mid ,up_flux       &
       ,dwn_flux ,entrain_up ,detrain_up ,entrain_dwn ,detrain_dwn      &
       ,up_flux_half ,T_incr_diag_conv ,q_incr_diag_conv                &
       ,qcl_incr_diag_conv ,qcf_incr_diag_conv                          & 
       ,cf_liquid_incr_diag_conv ,cf_frozen_incr_diag_conv              &
       ,bulk_cf_incr_diag_conv,u_incr_diag_conv ,v_incr_diag_conv       &
       ,theta_diag, q_diag                                              &
       ,mf_deep ,mf_congest ,mf_shall ,mf_midlev                        &
       ,dt_deep ,dt_congest ,dt_shall ,dt_midlev                        &
       ,dq_deep ,dq_congest ,dq_shall ,dq_midlev                        &
       ,du_deep ,du_congest ,du_shall ,du_midlev                        &
       ,dv_deep ,dv_congest ,dv_shall ,dv_midlev                        &
       ,wqt_flux_sh ,wql_flux_sh ,wthetal_flux_sh ,wthetav_flux_sh      &
       ,dubydt_pout ,dvbydt_pout ,conv_rain_3d ,conv_snow_3d            &
       ,t_incr_conv_only,q_incr_conv_only                               &
       ,qcl_incr_inhom_diag ,qcf_incr_inhom_diag                        &
       ,bulk_cf_incr_inhom_diag ,cf_liquid_incr_inhom_diag              &
       ,cf_frozen_incr_inhom_diag, deep_cfl_limited, mid_cfl_limited    &
       ,deep_tops

USE cv_stash_flg_mod, ONLY:                                             &
    l_apply_diag, l_qcl_incr_cinh, l_qcf_incr_cinh, l_cfl_incr_cinh     &
   ,l_cff_incr_cinh, l_bcf_incr_cinh

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Submodel_Mod

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   ALLOCATE and initialise  arrays required by convection for a full model 
!   timestep on the first convection substep of a model timestep.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::  &
  row_length            & ! Local number of points on a row
 ,rows                    ! Local number of rows in a theta field

LOGICAL, INTENT(IN)  ::  &
  l_calc_dxek            & ! Switch for calculation of condensate increment
 ,l_dust                 & ! Switch for mineral dust
 ,l_sulpc_so2            & ! Switch for sulphur cycle
 ,l_sulpc_nh3            & ! Switch for NH3 in S-cycle
 ,l_soot                 & ! Switch for soot cycle
 ,l_biomass              & ! Switch for biomass aerosol scheme
 ,l_ocff                 & ! Switch for fossil-fuel organic carbon scheme
 ,l_nitrate                ! Switch for ammonium nitrate aerosol

! Switch for COSP
LOGICAL, INTENT(IN)  ::  l_cosp

INTEGER, INTENT(IN) ::      &
  ntml(row_length, rows)    & ! Top level of surface mixed layer
 ,ntpar(row_length, rows)     ! Top level of initial parcel ascent

REAL, INTENT(IN)  ::                              &
  theta_inc    (tdims%i_start:tdims%i_end,        &
                tdims%j_start:tdims%j_end,        &
                            1:tdims%k_end)        &
 ,q_inc        (qdims%i_start:qdims%i_end,        &
                qdims%j_start:qdims%j_end,        &
                            1:qdims%k_end)        &
 ,qcl_inc      (qdims%i_start:qdims%i_end,        &
                qdims%j_start:qdims%j_end,        &
                            1:qdims%k_end)        &
 ,qcf_inc      (qdims%i_start:qdims%i_end,        &
                qdims%j_start:qdims%j_end,        &
                            1:qdims%k_end)        &
 ,cf_liquid_inc(qdims%i_start:qdims%i_end,        &
                qdims%j_start:qdims%j_end,        &
                            1:qdims%k_end)        &
 ,cf_frozen_inc(qdims%i_start:qdims%i_end,        &
                qdims%j_start:qdims%j_end,        &
                            1:qdims%k_end)        &
 ,bulk_cf_inc  (qdims%i_start:qdims%i_end,        &
                qdims%j_start:qdims%j_end,        &
                            1:qdims%k_end)

REAL,INTENT(OUT) ::                       &
  conv_rain(row_length,rows)              & ! convective rainfall
, conv_snow(row_length,rows)                ! convective snowfall

! required for stash variables in  "argsts.h"
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

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k                  ! loop counters

REAL ::                  &
  one_over_conv_calls      ! 1/n_conv_calls

! Required by Dr hook 
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------


IF (lhook) CALL dr_hook('CV_ALLOC_DIAG_ARRAY',zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Allocate 2d arrays for diagnostics 

 ALLOCATE(precip_deep(row_length,rows) )
 ALLOCATE(precip_shall(row_length,rows) )
 ALLOCATE(precip_mid(row_length,rows) ) 
 ALLOCATE(precip_cong(row_length,rows) )
 ALLOCATE(cape_out(row_length,rows) )
 ALLOCATE(deep_ind(row_length,rows) )
 ALLOCATE(shallow_ind(row_length,rows) )
 ALLOCATE(congestus_ind(row_length,rows) )
 ALLOCATE(congestus_ind2(row_length,rows) )
 ALLOCATE(mid_ind(row_length,rows) )
 ALLOCATE(ntml_diag(row_length,rows) ) 
 ALLOCATE(ntpar_diag(row_length,rows) )
 ALLOCATE(freeze_diag(row_length,rows) )
 ALLOCATE(kterm_diag(row_length,rows) )
 ALLOCATE(wstar_up_diag(row_length,rows) )
 ALLOCATE(wstar_dn_diag(row_length,rows) )
 ALLOCATE(mb1_diag(row_length,rows) )
 ALLOCATE(mb2_diag(row_length,rows) )
 ALLOCATE(wqt_cb(row_length, rows) )
 ALLOCATE(wthetal_cb(row_length, rows) ) 
 ALLOCATE(wqt_inv(row_length, rows) )
 ALLOCATE(wthetal_inv(row_length, rows) )
 ALLOCATE(sh_top(row_length, rows) )
 ALLOCATE(sh_base(row_length, rows) )
 ALLOCATE(cg_top(row_length, rows) )
 ALLOCATE(cg_base(row_length, rows) )
 ALLOCATE(cg_term(row_length, rows) )
 ALLOCATE(cape_ts_diag(row_length, rows) )
 ALLOCATE(ind_cape_reduced_diag(row_length, rows) )
 ALLOCATE(deep_cfl_limited(row_length, rows) )
 ALLOCATE(mid_cfl_limited(row_length, rows) )


! Initialise diagnostics

  one_over_conv_calls = 1.0/(n_conv_calls*1.0)

  DO j = 1, rows
    DO i = 1, row_length
      ntml_diag(i,j)  = 0.0
      ntpar_diag(i,j) = 0.0
      conv_rain(i,j) = 0.0
      conv_snow(i,j) = 0.0
      precip_deep(i,j)  = 0.0
      precip_shall(i,j) = 0.0
      precip_mid(i,j)   = 0.0
      precip_cong(i,j)  = 0.0
      cape_out(i,j) = 0.0
      kterm_diag(i,j) = 0.0
      freeze_diag(i,j) = 0.0
      deep_ind(i,j)       = 0.0
      shallow_ind(i,j)    = 0.0
      congestus_ind(i,j)  = 0.0
      congestus_ind2(i,j) = 0.0
      mid_ind(i,j)        = 0.0
      ind_cape_reduced_diag(i,j) = 0.0
      cape_ts_diag(i,j) = 0.0
      deep_cfl_limited(i,j) = 0.0
      mid_cfl_limited(i,j)  = 0.0

! 5A & 6A diagnostics
      wstar_up_diag(i,j) = 0.0
      wstar_dn_diag(i,j) = 0.0
      mb1_diag(i,j) = 0.0
      mb2_diag(i,j) = 0.0
      cg_term(i,j) = 0.0
      cg_top(i,j)  = 0.0
      cg_base(i,j) = 0.0
      wqt_cb(i,j)     = 0.0
      wthetal_cb(i,j) = 0.0
      wqt_inv(i,j)     = 0.0
      wthetal_inv(i,j) = 0.0
      sh_top(i,j)  = 0.0
      sh_base(i,j) = 0.0
    END DO
  END DO

! Convective aerosol relative diagnostics - allocation dependent on whether
! active or not.

IF (l_dust) THEN
  ALLOCATE(conscav_dust(row_length, rows, 6))
ELSE 
  ALLOCATE(conscav_dust(1,1,1))
END IF

IF (l_soot) THEN
  ALLOCATE(conscav_agedsoot(row_length, rows) )
ELSE
  ALLOCATE(conscav_agedsoot(1,1) )
END IF

IF (l_biomass) THEN
  ALLOCATE(conscav_agedbmass(row_length, rows) )
ELSE
  ALLOCATE(conscav_agedbmass(1,1) )
END IF
IF (l_ocff) THEN
  ALLOCATE(conscav_agedocff(row_length, rows) )
ELSE
  ALLOCATE(conscav_agedocff(1,1) )
END IF

IF (l_nitrate) THEN
  ALLOCATE(conscav_nitracc(row_length, rows) )
  ALLOCATE(conscav_nitrdiss(row_length, rows) )
ELSE
  ALLOCATE(conscav_nitracc(1,1) )
  ALLOCATE(conscav_nitrdiss(1,1) )
END IF                                    
IF (l_sulpc_so2) THEN
  ALLOCATE(conwash_so2(row_length, rows) )
  ALLOCATE(conscav_so4ait(row_length, rows) )
  ALLOCATE(conscav_so4acc(row_length, rows) )
  ALLOCATE(conscav_so4dis(row_length, rows) )
ELSE
  ALLOCATE(conwash_so2(1,1) )
  ALLOCATE(conscav_so4ait(1,1) )
  ALLOCATE(conscav_so4acc(1,1) )
  ALLOCATE(conscav_so4dis(1,1) )
END IF
IF (l_sulpc_nh3 .AND. l_sulpc_so2 ) THEN
  ALLOCATE(conwash_nh3(row_length, rows) )
ELSE
  ALLOCATE(conwash_nh3(1,1) )
END IF

! End of timestep values of theta and q after convection

 ALLOCATE ( theta_diag(tdims%i_start:tdims%i_end,        &
                       tdims%j_start:tdims%j_end,        &
                                   1:tdims%k_end) )
 ALLOCATE ( q_diag(qdims%i_start:qdims%i_end,        &
                   qdims%j_start:qdims%j_end,        &
                               1:qdims%k_end) )



!---------------------------------------------------------------------
! Intercept values of physics increments before being updated by
! convection for optional output of convective increments
!---------------------------------------------------------------------

! Convection increments

  IF (((sf(181,5).OR.sf(187,5).OR.sf(161,5)) .AND. L_apply_diag)     &
                 .OR. L_calc_dxek) THEN
    ALLOCATE ( T_incr_diag_conv(tdims%i_start:tdims%i_end,    &
                                tdims%j_start:tdims%j_end,    &
                                            1:tdims%k_end))
    DO k =             1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          T_incr_diag_conv(i,j,k) = theta_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

  ELSE
    ALLOCATE ( T_incr_diag_conv(1,1,1) )
  END IF

  IF (((sf(182,5).OR.sf(188,5).OR.sf(162,5)) .AND. L_apply_diag)      &
                 .OR. L_calc_dxek) THEN
    ALLOCATE (q_incr_diag_conv(qdims%i_start:qdims%i_end,    &
                               qdims%j_start:qdims%j_end,    &
                                           1:qdims%k_end))
    ! Hold input q increment
    DO k =             1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          q_incr_diag_conv(i,j,k) = q_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( q_incr_diag_conv(1,1,1) )
  END IF

  IF ( ( sf(183,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) THEN
    ALLOCATE (qcl_incr_diag_conv(qdims%i_start:qdims%i_end,    &
                                 qdims%j_start:qdims%j_end,    &
                                             1:qdims%k_end))
    ! Hold input qcl increment
    DO k =             1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          qcl_incr_diag_conv(i,j,k) = qcl_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE (qcl_incr_diag_conv(1,1,1))
  END IF                   ! on STASHflag

  IF ( ( sf(184,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) THEN
    ALLOCATE (qcf_incr_diag_conv(qdims%i_start:qdims%i_end,    &
                                 qdims%j_start:qdims%j_end,    &
                                             1:qdims%k_end))

    ! Hold input qcf increment
    DO k =             1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          qcf_incr_diag_conv(i,j,k) = qcf_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE (qcf_incr_diag_conv(1,1,1))
  END IF                   ! on STASHflag

  IF ( ( sf(193,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) THEN
    ALLOCATE(cf_liquid_incr_diag_conv(qdims%i_start:qdims%i_end,    &
                                      qdims%j_start:qdims%j_end,    &
                                                  1:qdims%k_end))

    ! Hold input cf_liquid increment
    DO k =             1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          cf_liquid_incr_diag_conv(i,j,k) = cf_liquid_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE(cf_liquid_incr_diag_conv(1,1,1))
  END IF                   ! on STASHflag

  IF ( ( sf(194,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) THEN
    ALLOCATE(cf_frozen_incr_diag_conv(qdims%i_start:qdims%i_end,    &
                                      qdims%j_start:qdims%j_end,    &
                                                  1:qdims%k_end))

    ! Hold input cf_frozen increment
    DO k =             1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          cf_frozen_incr_diag_conv(i,j,k) = cf_frozen_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE(cf_frozen_incr_diag_conv(1,1,1))
  END IF                   ! on STASHflag

  IF ( ( sf(195,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) THEN
    ALLOCATE(bulk_cf_incr_diag_conv(qdims%i_start:qdims%i_end,    &
                                    qdims%j_start:qdims%j_end,    &
                                                1:qdims%k_end))

    ! Hold input bulk_cf increment
    DO k =             1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          bulk_cf_incr_diag_conv(i,j,k) = bulk_cf_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE(bulk_cf_incr_diag_conv(1,1,1))
  END IF                   ! on STASHflag

  IF ( sf(185,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( u_incr_diag_conv(udims%i_start:udims%i_end,      &
                                udims%j_start:udims%j_end,      &
                                udims%k_start:udims%k_end) )
    DO k = udims%k_start, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          u_incr_diag_conv(i,j,k) = 0.0
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( u_incr_diag_conv(1,1,1) )
  END IF

  IF ( sf(186,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( v_incr_diag_conv(vdims%i_start:vdims%i_end,      & 
                                vdims%j_start:vdims%j_end,      & 
                                vdims%k_start:vdims%k_end) )
    DO k = vdims%k_start, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          v_incr_diag_conv(i,j,k) = 0.0
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( v_incr_diag_conv(1,1,1) )
  END IF

  IF (sf(161,5) .AND. L_apply_diag) THEN
    ALLOCATE ( T_incr_conv_only(tdims%i_start:tdims%i_end,    &
                                tdims%j_start:tdims%j_end,    &
                                            1:tdims%k_end))
    DO k =             1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          T_incr_conv_only(i,j,k) = 0.0
        END DO ! i
      END DO ! j
    END DO ! k

  ELSE
    ALLOCATE ( T_incr_conv_only(1,1,1) )
  END IF

  IF (sf(162,5) .AND. L_apply_diag) THEN
    ALLOCATE (q_incr_conv_only(qdims%i_start:qdims%i_end,    &
                               qdims%j_start:qdims%j_end,    &
                                           1:qdims%k_end))
    ! Hold input q increment
    DO k =             1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          q_incr_conv_only(i,j,k) = 0.0
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( q_incr_conv_only(1,1,1) )
  END IF

! Other model level diagnostics :
! setup allocatable arrays to limit memory over head if diagnostics not
! required.

  IF ( ( sf(227,5) .AND. L_apply_diag ) .OR. L_cosp ) THEN 
    ALLOCATE ( conv_rain_3d(qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end) )
    DO k=1, qdims%k_end
      DO j=qdims%j_start,qdims%j_end
        DO i=qdims%i_start,qdims%i_end
          conv_rain_3d(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( conv_rain_3d(1,1,1) )
  END IF
  IF ( ( sf(228,5) .AND. L_apply_diag ) .OR. L_cosp ) THEN
    ALLOCATE ( conv_snow_3d(qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end) )
    DO k=1, qdims%k_end
      DO j=qdims%j_start,qdims%j_end
        DO i=qdims%i_start,qdims%i_end
          conv_snow_3d(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( conv_snow_3d(1,1,1) )
  END IF

  IF ( (sf(290,5) .OR. sf(429,5) .OR. sf(430,5) .OR.          &
        sf(431,5) .OR. sf(432,5) ) .AND. L_apply_diag ) THEN
    ALLOCATE ( wqt_flux_sh(qdims%i_start:qdims%i_end,         &
                           qdims%j_start:qdims%j_end,         &
                                       1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
          wqt_flux_sh(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( wqt_flux_sh(1,1,1) )
  END IF

  IF ( sf(291,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( wql_flux_sh(qdims%i_start:qdims%i_end,        &
                           qdims%j_start:qdims%j_end,        &
                                       1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
          wql_flux_sh(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( wql_flux_sh(1,1,1) )
  END IF

  IF ( (sf(292,5) .OR. sf(425,5) .OR. sf(426,5) .OR.          &
        sf(427,5) .OR. sf(428,5) ) .AND. L_apply_diag ) THEN
    ALLOCATE ( wthetal_flux_sh(qdims%i_start:qdims%i_end,     &
                               qdims%j_start:qdims%j_end,     &
                                           1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
          wthetal_flux_sh(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( wthetal_flux_sh(1,1,1) )
  END IF

  IF ( sf(293,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( wthetav_flux_sh(qdims%i_start:qdims%i_end,        &
                               qdims%j_start:qdims%j_end,        &
                                           1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
          wthetav_flux_sh(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( wthetav_flux_sh(1,1,1) )
  END IF

 ! Require these arrays whether diagnostics required or not as used to
 ! update r_u and r_v
  IF (l_mom) THEN
    ALLOCATE ( dubydt_pout(tdims%i_start:tdims%i_end,        &
                           tdims%j_start:tdims%j_end,        &
                                       1:tdims%k_end) )
    ALLOCATE ( dvbydt_pout(tdims%i_start:tdims%i_end,        &
                           tdims%j_start:tdims%j_end,        &
                                       1:tdims%k_end) )
    DO k=1,tdims%k_end
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
           dubydt_pout(i,j,k) = 0.0
           dvbydt_pout(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE     ! No CMT so no wind increments
    ALLOCATE ( dubydt_pout(1,1,1) )
    ALLOCATE ( dvbydt_pout(1,1,1) )
  END IF   ! l_mom

  IF ( sf(320,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( mf_deep(qdims%i_start:qdims%i_end,        &
                       qdims%j_start:qdims%j_end,        &
                                   1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           mf_deep(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( mf_deep(1,1,1) )
  END IF
  IF ( sf(321,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( mf_congest(qdims%i_start:qdims%i_end,        &
                          qdims%j_start:qdims%j_end,        &
                                      1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           mf_congest(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( mf_congest(1,1,1) )
  END IF

  IF ( ( sf(322,5) .OR. sf(417,5) .OR. sf(418,5) .OR.         &
         sf(419,5) .OR. sf(420,5) ) .AND. L_apply_diag ) THEN
    ALLOCATE ( mf_shall(qdims%i_start:qdims%i_end,            &
                        qdims%j_start:qdims%j_end,            &
                                    1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           mf_shall(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( mf_shall(1,1,1) )
  END IF

  IF ( sf(323,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( mf_midlev(qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           mf_midlev(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( mf_midlev(1,1,1) )
  END IF
  IF ( sf(324,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( dt_deep(qdims%i_start:qdims%i_end,        &
                       qdims%j_start:qdims%j_end,        &
                                   1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dt_deep(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dt_deep(1,1,1) )
  END IF
  IF ( sf(325,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( dt_congest(qdims%i_start:qdims%i_end,        &
                          qdims%j_start:qdims%j_end,        &
                                      1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dt_congest(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dt_congest(1,1,1) )
  END IF

  IF ( ( sf(326,5) .OR. sf(409,5) .OR. sf(410,5) .OR.         &
             sf(411,5) .OR. sf(412,5) ) .AND. L_apply_diag ) THEN
    ALLOCATE ( dt_shall(qdims%i_start:qdims%i_end,            &
                        qdims%j_start:qdims%j_end,            &
                                    1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dt_shall(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dt_shall(1,1,1) )
  END IF

  IF ( sf(327,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( dt_midlev(qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dt_midlev(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dt_midlev(1,1,1) )
  END IF
  IF ( sf(328,5).AND. L_apply_diag ) THEN
    ALLOCATE ( dq_deep(qdims%i_start:qdims%i_end,        &
                       qdims%j_start:qdims%j_end,        &
                                   1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dq_deep(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dq_deep(1,1,1) )
  END IF
  IF ( sf(329,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( dq_congest(qdims%i_start:qdims%i_end,        &
                          qdims%j_start:qdims%j_end,        &
                                      1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dq_congest(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dq_congest(1,1,1) )
  END IF

  IF ( ( sf(330,5) .OR. sf(413,5) .OR. sf(414,5) .OR.         &
         sf(415,5) .OR. sf(416,5) ) .AND. L_apply_diag ) THEN
    ALLOCATE ( dq_shall(qdims%i_start:qdims%i_end,            &
                        qdims%j_start:qdims%j_end,            &
                                    1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dq_shall(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dq_shall(1,1,1) )
  END IF

  IF ( sf(331,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( dq_midlev(qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dq_midlev(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dq_midlev(1,1,1) )
  END IF
  IF ( sf(332,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( du_deep(qdims%i_start:qdims%i_end,        &
                       qdims%j_start:qdims%j_end,        &
                                   1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           du_deep(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( du_deep(1,1,1) )
  END IF
  IF ( sf(333,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( du_congest(qdims%i_start:qdims%i_end,        &
                          qdims%j_start:qdims%j_end,        &
                                      1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           du_congest(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( du_congest(1,1,1) )
  END IF
  IF ( sf(334,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( du_shall(qdims%i_start:qdims%i_end,        &
                        qdims%j_start:qdims%j_end,        &
                                    1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           du_shall(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( du_shall(1,1,1) )
  END IF
  IF ( sf(335,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( du_midlev(qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           du_midlev(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( du_midlev(1,1,1) )
  END IF
  IF ( sf(336,5).AND. L_apply_diag ) THEN
    ALLOCATE ( dv_deep(qdims%i_start:qdims%i_end,        &
                       qdims%j_start:qdims%j_end,        &
                                   1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dv_deep(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dv_deep(1,1,1) )
  END IF
  IF ( sf(337,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( dv_congest(qdims%i_start:qdims%i_end,        &
                          qdims%j_start:qdims%j_end,        &
                                      1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dv_congest(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dv_congest(1,1,1) )
  END IF
  IF ( sf(338,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( dv_shall(qdims%i_start:qdims%i_end,        &
                        qdims%j_start:qdims%j_end,        &
                                    1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dv_shall(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dv_shall(1,1,1) )
  END IF
  IF ( sf(339,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( dv_midlev(qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start,qdims%i_end
           dv_midlev(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( dv_midlev(1,1,1) )
  END IF

  IF ( (sf(249,5) .OR. sf(246,5) ) .AND. L_apply_diag ) THEN
    ALLOCATE ( up_flux_half(pdims%i_start:pdims%i_end,        &
                            pdims%j_start:pdims%j_end,        &
                            pdims%k_start:pdims%k_end) )
    up_flux_half(:,:,:) = 0.0
  ELSE
    ALLOCATE ( up_flux_half(1,1,1) )
  END IF

! Note always require mass flux up and down as not just required in 
! atmos_physics2 for output diagnostics

  ALLOCATE ( up_flux(qdims%i_start:qdims%i_end,           &
                     qdims%j_start:qdims%j_end,           &
                                 1:qdims%k_end) )
  ALLOCATE ( dwn_flux(qdims%i_start:qdims%i_end,          &
                      qdims%j_start:qdims%j_end,          &
                                  1:qdims%k_end) )

  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start,qdims%i_end
        up_flux(i,j,k)=0.0
        dwn_flux(i,j,k) = 0.0
      END DO
    END DO
  END DO

  IF ( sf(252,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( entrain_up(qdims%i_start:qdims%i_end,        &
                          qdims%j_start:qdims%j_end,        &
                                      1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j=qdims%j_start,qdims%j_end
        DO i= qdims%i_start,qdims%i_end
          entrain_up(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( entrain_up(1,1,1) )
  END IF

  IF ( sf(253,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( detrain_up(qdims%i_start:qdims%i_end,        &
                          qdims%j_start:qdims%j_end,        &
                                      1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j=qdims%j_start,qdims%j_end
        DO i= qdims%i_start,qdims%i_end
          detrain_up(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( detrain_up(1,1,1) )
  END IF

  IF ( sf(254,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( entrain_dwn(qdims%i_start:qdims%i_end,        &
                           qdims%j_start:qdims%j_end,        &
                                       1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j=qdims%j_start,qdims%j_end
        DO i= qdims%i_start,qdims%i_end
          entrain_dwn(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( entrain_dwn(1,1,1) )
  END IF

  IF ( sf(255,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( detrain_dwn(qdims%i_start:qdims%i_end,        &
                           qdims%j_start:qdims%j_end,        &
                                       1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j=qdims%j_start,qdims%j_end
        DO i= qdims%i_start,qdims%i_end
          detrain_dwn(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( detrain_dwn(1,1,1) )
  END IF

  IF ( sf(258,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( uw_dp(tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                                 1:tdims%k_end) )
    DO k=1,tdims%k_end
      DO j=tdims%j_start,tdims%j_end
        DO i= tdims%i_start,tdims%i_end
          uw_dp(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( uw_dp(1,1,1) )
  END IF


  IF ( sf(259,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( vw_dp(tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                                 1:tdims%k_end) )
    DO k=1,tdims%k_end
      DO j=tdims%j_start,tdims%j_end
        DO i= tdims%i_start,tdims%i_end
          vw_dp(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( vw_dp(1,1,1) )
  END IF

  IF ( sf(260,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( uw_shall(tdims%i_start:tdims%i_end,        &
                        tdims%j_start:tdims%j_end,        &
                                    1:tdims%k_end) )
    DO k=1,tdims%k_end
      DO j=tdims%j_start,tdims%j_end
        DO i= tdims%i_start,tdims%i_end
          uw_shall(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( uw_shall(1,1,1) )
  END IF


  IF ( sf(261,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( vw_shall(tdims%i_start:tdims%i_end,        &
                        tdims%j_start:tdims%j_end,        &
                                    1:tdims%k_end) )
    DO k=1,tdims%k_end
      DO j=tdims%j_start,tdims%j_end
        DO i= tdims%i_start,tdims%i_end
          vw_shall(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( vw_shall(1,1,1) )
  END IF

  IF ( sf(263,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( uw_mid(tdims%i_start:tdims%i_end,        &
                      tdims%j_start:tdims%j_end,        &
                                  1:tdims%k_end) )
    DO k=1,tdims%k_end
      DO j=tdims%j_start,tdims%j_end
        DO i= tdims%i_start,tdims%i_end
          uw_mid(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( uw_mid(1,1,1) )
  END IF


  IF ( sf(264,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( vw_mid(tdims%i_start:tdims%i_end,        &
                      tdims%j_start:tdims%j_end,        &
                                  1:tdims%k_end) )
    DO k=1,tdims%k_end
      DO j=tdims%j_start,tdims%j_end
        DO i= tdims%i_start,tdims%i_end
          vw_mid(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( vw_mid(1,1,1) )
  END IF

  IF ( sf(319,5) .AND. L_apply_diag ) THEN
    ALLOCATE ( deep_tops(qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end) )
    DO k=1,qdims%k_end
      DO j=qdims%j_start,qdims%j_end
        DO i= qdims%i_start,qdims%i_end
          deep_tops(i,j,k) = 0.0
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( deep_tops(1,1,1) )
  END IF



! inhomogenous diagnostics

  IF(l_qcl_incr_cinh) THEN

    ALLOCATE (qcl_incr_inhom_diag(qdims%i_start:qdims%i_end,        &
                                  qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end) )

    ! Hold input qcl increment (almost certainly zero)
    DO k=1,qdims%k_end
      DO j=1,rows
        DO i=1,row_length
          qcl_incr_inhom_diag(i,j,k) = qcl_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE (qcl_incr_inhom_diag(1,1,1))
  END IF                   ! on STASHflag

  IF(l_qcf_incr_cinh) THEN

    ALLOCATE (qcf_incr_inhom_diag(qdims%i_start:qdims%i_end,        &
                                  qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end) )
    ! Hold input qcf increment (almost certainly zero)
    DO k=1,qdims%k_end
      DO j=1,rows
        DO i=1,row_length
          qcf_incr_inhom_diag(i,j,k) = qcf_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE (qcf_incr_inhom_diag(1,1,1))

  END IF                   ! on STASHflag

  IF(l_bcf_incr_cinh) THEN

    ALLOCATE (bulk_cf_incr_inhom_diag(qdims%i_start:qdims%i_end,        &
                                      qdims%j_start:qdims%j_end,        &
                                                  1:qdims%k_end) )

    ! Hold input bulk_cf increment (almost certainly zero)
    DO k=1,qdims%k_end
      DO j=1,rows
        DO i=1,row_length
          bulk_cf_incr_inhom_diag(i,j,k) = bulk_cf_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE (bulk_cf_incr_inhom_diag(1,1,1))

  END IF                   ! on STASHflag

  IF(l_cfl_incr_cinh .OR. l_calc_dxek) THEN

    ALLOCATE (cf_liquid_incr_inhom_diag(qdims%i_start:qdims%i_end,        &
                                        qdims%j_start:qdims%j_end,        &
                                                    1:qdims%k_end) )
    ! Hold input bulk_cf increment (almost certainly zero)
    DO k=1,qdims%k_end
      DO j=1,rows
        DO i=1,row_length
          cf_liquid_incr_inhom_diag(i,j,k) = cf_liquid_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE (cf_liquid_incr_inhom_diag(1,1,1))

  END IF                   ! on STASHflag

  IF(l_cff_incr_cinh) THEN

    ALLOCATE (cf_frozen_incr_inhom_diag(qdims%i_start:qdims%i_end,        &
                                        qdims%j_start:qdims%j_end,        &
                                                    1:qdims%k_end) )

    ! Hold input bulk_cf increment (almost certainly zero)
    DO k=1,qdims%k_end
      DO j=1,rows
        DO i=1,row_length
          cf_frozen_incr_inhom_diag(i,j,k) = cf_frozen_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE (cf_frozen_incr_inhom_diag(1,1,1))

  END IF                   ! on STASHflag


!------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('CV_ALLOC_DIAG_ARRAY',zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE cv_alloc_diag_array
