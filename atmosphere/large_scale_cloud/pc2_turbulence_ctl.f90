! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the effect of changing the width of the PDF used in PC2.
!
SUBROUTINE pc2_turbulence_ctl (                                         &
! Primary fields passed in/out - these will be unchanged if calling from
! AP1 and incremented if calling from microphys_ctl
 T, q, qcl, cf, cfl, cff,                                               &
 p_theta_levels,                                                        &
! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
 STASHwork4,                                                            &
! SCM diagnostics switches (dummy in full UM)
 nSCMDpkgs, L_SCMDiags,                                                 &
! increments passed in/out - these are unchanged if calling
! from microphys_ctl but updated if calling from AP1
 T_inc, q_inc, qcl_inc, cf_inc, cfl_inc)
!
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE cloud_inputs_mod, ONLY: dbsdtbs_turb_0, l_micro_eros,             &
                              i_pc2_erosion_method
  USE pc2_constants_mod, ONLY: dbsdtbs_turb_1, pc2eros_hybrid_allfaces
  USE atm_fields_bounds_mod, ONLY: tdims, qdims, pdims
  USE um_input_control_mod, ONLY: l_mr_physics1
  USE timestep_mod, ONLY: timestep
  USE proc_info_mod, ONLY: at_extremity
  USE ereport_mod, ONLY : ereport
  USE Submodel_Mod
!
  IMPLICIT NONE
!
! Description: Condensation and cloud fraction changes in the PC2
!              framework as a result of changing the width of the
!              moisture PDF without changing its shape.
!
! Method:      Uses the equations outlined in the PC2 cloud scheme
!              documentation UMDP029A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! Primary fields passed in
  REAL                                                                  &
   T(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tdims%k_end),  &
   q(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,qdims%k_end),  &
   qcl(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,             &
         qdims%k_end),                                                  &
   cf(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,              &
        qdims%k_end),                                                   &
   cfl(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,             &
         qdims%k_end),                                                  &
   cff(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,             &
         qdims%k_end),                                                  &
   p_theta_levels(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  pdims%k_end)
!
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
!
! Diagnostics info
  REAL                                                                  &
   STASHwork4(*)     ! STASH workspace
!
! Additional variables for SCM diagnostics
  INTEGER                                                               &
   nSCMDpkgs              ! No of SCM diagnostics packages
!
  LOGICAL                                                               &
   L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages
!
  REAL                                                                  &
   T_inc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tdims%k_end),                                                  &
   q_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,           &
         qdims%k_end),                                                  &
   qcl_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,         &
           qdims%k_end),                                                &
   cf_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,          &
          qdims%k_end),                                                 &
   cfl_inc(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,         &
           qdims%k_end)
!
! Local variables
!
  CHARACTER(LEN=*), Parameter ::  RoutineName = 'pc2_turbulence_ctl'
  CHARACTER(LEN=80) cmessage
!
  INTEGER                                                               &
   i,j,k,                             & ! Loop counters
   icode,item,im_index
!
  INTEGER, PARAMETER :: sect = 4        ! STASH section for microphys
!
  REAL                                                                  &
   T_work(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tdims%k_end),                                                 &
   q_work(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,          &
          qdims%k_end),                                                 &
   qcl_work(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,        &
            qdims%k_end),                                               &
   cf_work(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,         &
           qdims%k_end),                                                &
   cfl_work(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,        &
            qdims%k_end),                                               &
   zeros(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,qdims%k_end)
!
  REAL,ALLOCATABLE::                                                    &
   cf_above(:,:,:), cf_below(:,:,:)
!
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
!
!- End of header
!
  IF (lhook) CALL dr_hook('PC2_TURBULENCE_CTL',zhook_in,zhook_handle)
!
! Define zeros array
  zeros(:,:,:) = 0.0
!
  IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces) THEN
    ALLOCATE(cf_above(qdims%i_start:qdims%i_end,                        &
                      qdims%j_start:qdims%j_end,qdims%k_end))
    ALLOCATE(cf_below(qdims%i_start:qdims%i_end,                        &
                      qdims%j_start:qdims%j_end,qdims%k_end))
    k = 1
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        cf_above(i,j,k) = cf(i,j,k+1)
        cf_below(i,j,k) = cf(i,j,k)
      END DO
    END DO
    DO k = 2, qdims%k_end-1
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          cf_above(i,j,k) = cf(i,j,k+1)
          cf_below(i,j,k) = cf(i,j,k-1)
        END DO
      END DO
    END DO
    k = qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        cf_above(i,j,k) = cf(i,j,k)
        cf_below(i,j,k) = cf(i,j,k-1)
      END DO
    END DO
  ELSE
    ALLOCATE(cf_above(1,1,1))
    ALLOCATE(cf_below(1,1,1))
  END IF
!
! Call homogenous forcing routine - this does not add the increments
! and outputs them separately
! DEPENDS ON: pc2_hom_conv
  CALL PC2_HOM_CONV(p_theta_levels,qdims%k_end,timestep,                &
                    T, q, qcl, cf, cfl, cff,                            &
                    zeros, zeros, zeros, zeros, zeros,                  &
                    cf_above, cf_below,                                 &
                    t_work, q_work, qcl_work, cf_work, cfl_work,        &
                    dbsdtbs_turb_0, dbsdtbs_turb_1, l_mr_physics1)
!
  DEALLOCATE(cf_below)
  DEALLOCATE(cf_above)
!
  IF (l_micro_eros) THEN
! if calling from microphys ctl, update the main fields
    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          q(i,j,k)   = q(i,j,k)   + q_work(i,j,k)
          qcl(i,j,k) = qcl(i,j,k) + qcl_work(i,j,k)
          cf(i,j,k)  = cf(i,j,k)  + cf_work(i,j,k)
          cfl(i,j,k) = cfl(i,j,k) + cfl_work(i,j,k)
          t(i,j,k)   = t(i,j,k)   + t_work(i,j,k)
        END DO
      END DO
    END DO 
  ELSE
! if calling from AP1,  update the increment fields
    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          q_inc(i,j,k)   = q_work(i,j,k)   + q_inc(i,j,k)
          qcl_inc(i,j,k) = qcl_work(i,j,k) + qcl_inc(i,j,k)
          cf_inc(i,j,k)  = cf_work(i,j,k)  + cf_inc(i,j,k)
          cfl_inc(i,j,k) = cfl_work(i,j,k) + cfl_inc(i,j,k)
          T_inc(i,j,k)   = T_work(i,j,k)   + T_inc(i,j,k)
        END DO
      END DO
    END DO
  END IF
!
! ----------------------------------------------------------------------
! Output Diagnostics
! ----------------------------------------------------------------------
!
!
  icode = 0 ! Initialise error status
  im_index = internal_model_index(atmos_im)
! ----------------------------------------------------------------------
! DIAG.04281 Copy hom_conv qcl incr to stashwork
! ----------------------------------------------------------------------
  item = 281
  IF (sf(item,sect)) THEN
!
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),               &
        t_work,                                                         &
        tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)
!
    IF (icode /=  0) THEN
      CALL ereport(RoutineName,icode,cmessage)
      icode = 0
    END IF
  END IF
!
! ----------------------------------------------------------------------
! DIAG.04282 Copy hom_conv qcl incr to stashwork
! ----------------------------------------------------------------------
  item = 282
  IF (sf(item,sect)) THEN
!
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),               &
        q_work,                                                         &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)
!
    IF (icode /=  0) THEN
      CALL ereport(RoutineName,icode,cmessage)
      icode = 0
    END IF
  END IF
!
! ----------------------------------------------------------------------
! DIAG.04283 Copy hom_conv qcl incr to stashwork
! ----------------------------------------------------------------------
  item = 283
  IF (sf(item,sect)) THEN
!
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),               &
        qcl_work,                                                       &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)
!
    IF (icode /=  0) THEN
      CALL ereport(RoutineName,icode,cmessage)
      icode = 0
    END IF
  END IF
!
! ----------------------------------------------------------------------
! DIAG.04292 Copy hom_conv qcl incr to stashwork
! ----------------------------------------------------------------------
  item = 292
  IF (sf(item,sect)) THEN
!
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),               &
        cf_work,                                                        &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)
!
    IF (icode /=  0) THEN
      CALL ereport(RoutineName,icode,cmessage)
      icode = 0
    END IF
  END IF
!
! ----------------------------------------------------------------------
! DIAG.04293 Copy hom_conv cfl incr to stashwork
! ----------------------------------------------------------------------
  item = 293
  IF (sf(item,sect)) THEN
!
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),               &
       cfl_work,                                                        &
       qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)
!
    IF (icode /=  0) THEN
      CALL ereport(RoutineName,icode,cmessage)
      icode = 0
    END IF
  END IF
!
!
  IF (lhook) CALL dr_hook('PC2_TURBULENCE_CTL',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_turbulence_ctl
