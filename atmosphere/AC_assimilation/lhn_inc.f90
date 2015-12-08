! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE LHN_INC ------------------------------------------------
!LL
!LL  Purpose : Latent Heat nudging of model heating profiles
!LL
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL  Documentation :  FR  WP 171
!LL
!LLEND------------------------------------------------------------------
!L*  ARGUMENTS
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE lhn_inc_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE LHN_INC(CONV_HEAT,LS_HEAT,LS_RAIN,LS_SNOW,CONV_RAIN,   &
     &                   CONV_SNOW,P_FIELD,Q_LEVELS,PR_INC,TIMESTEP,    &
     &                   IROWS,ICOLS,                                   &
     &                   RANGE,phi_p,lambda_p,L_regular,                &
     &                   THINCS,LENMG,ICODE,CMESSAGE)

      USE earth_constants_mod, ONLY: earth_radius

      USE conversions_mod, ONLY: pi_over_180
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      USE lhn_search_mod, ONLY: lhn_search
      USE rfcsl_mod, ONLY: rfcsl
      USE rfcslr_mod, ONLY: rfcslr
      IMPLICIT NONE

!-----AC common blocks
! MPPAC start
      ! dimensions for obs allocated in subroutine AC
      ! dynamically allocated in AC and RDOBS

      ! dimension inxdim allocated in subroutine HORINF
      INTEGER,PARAMETER:: inxdim    = 15000

      ! common for Statistics Calcs in DIAGO ; Prints in RDOBS,GETOBS
      REAL :: R_STAT(MODEL_LEVELS_MAX,0:8)
      REAL :: S_STAT(MODEL_LEVELS_MAX,0:8)
      INTEGER :: COUNTA(NOBTYPMX)
      INTEGER :: COUNTB(NOBTYPMX)
      INTEGER :: COUNTC(NOBTYPMX)

      COMMON /mpp_ac/ R_STAT,S_STAT,COUNTA,COUNTB,COUNTC

      ! common to pass longitudes and latitudes for edges of local
      ! box from setcona to RDOBS and HINTCF

      REAL :: LONG_E
      REAL :: LONG_W
      REAL :: LAT_N,LAT_S
      REAL :: LONG_W_MODEL
      REAL :: LONG_E_MODEL
      COMMON/latlonmax/                                                 &
     &  LONG_E,LONG_W,LAT_N,LAT_S,LONG_W_MODEL,LONG_E_MODEL
! MPPAC end
!-----------------------------------------------------------------------
!*LCOMDECK COMACP
! Description:
!   Declares variables for the common block COMACP. This controls the
!  execution of the assimilation code.
!-----------------------------------------------------------------------

! Declare parameters:
      INTEGER      MODEACP
        PARAMETER (MODEACP = 36)
      INTEGER      NANALTYP
        PARAMETER (NANALTYP = 30)
      INTEGER      NRADARS
        PARAMETER (NRADARS  = 15)

! Declare variables:
      INTEGER NACT,                        NPROG
      INTEGER AC_OBS_TYPES(NOBTYPMX),      LACT(NOBTYPMX)
      INTEGER GROUP_INDEX(NOBTYPMX),       TYPE_INDEX(NOBTYPMX)
      INTEGER GROUP_FIRST(NOBTYPMX),       GROUP_LAST(NOBTYPMX)
      INTEGER OBS_UNITNO,                  OBS_FORMAT
      INTEGER NO_OBS_FILES,                DIAG_RDOBS
      INTEGER IUNITNO,                     MGEOWT
      INTEGER N_GROUPS,                    GROUP_NO(NOBTYPMX)
      INTEGER MHCORFN,                     MACDIAG(MODEACP)
      INTEGER MWTFN,                       MDATADFN
      INTEGER NPASS_RF,                    NSLABS_SCFACT(MODEL_LEVELS_MAX)
      INTEGER NO_SCFACT(NOBTYPMX),         IOMITOBS(NANALTYP)
      INTEGER MASTER_AC_TYPES(NOBTYPMX),   DEF_AC_ORDER(NOBTYPMX)
      INTEGER DEF_NO_ITERATIONS(NOBTYPMX), DEF_INTERVAL_ITER(NOBTYPMX)
      INTEGER DEF_NO_ANAL_LEVS(NOBTYPMX),  DEF_NO_WT_LEVS(NOBTYPMX)
      INTEGER DEF_MODE_HANAL(NOBTYPMX),    LENACT(NOBTYPMX)
      INTEGER DEF_OBTHIN(NOBTYPMX),        MVINT205
      INTEGER MRAMPFN,                     MGLOSSFN
      INTEGER      LHN_RANGE
      INTEGER      NPASS_RF_LHN

      INTEGER WB_LonOffset,                WB_LonPts
      INTEGER WB_LatOffset,                WB_LatPts

      REAL    OBTIME_NOM
      REAL    VERT_FILT
      REAL    GEOWT_H(ROWS_MAX -1)
      REAL    TROPLAT
      REAL    GEOWT_V(MODEL_LEVELS_MAX)
      REAL    VERT_COR_SCALE(MODEL_LEVELS_MAX, 4)
      REAL    VERT_CUTOFF_SL
      REAL    VERT_CUTOFF_BW
      REAL    VERT_CUTOFF_BH
      REAL    NON_DIV_COR
      REAL NON_DIV_COR_10M
      REAL    SPEED_LIMIT305
      REAL    TROPINT
      REAL    TIMEF_START
      REAL    TIMEF_OBTIME
      REAL    TIMEF_END
      REAL    CSCFACT_H(ROWS_MAX)
      REAL    CSCFACT_V(MODEL_LEVELS_MAX)
      REAL    DEF_TIMEB(NOBTYPMX)
      REAL    DEF_TIMEA(NOBTYPMX)
      REAL    DEF_TGETOBB(NOBTYPMX)
      REAL    DEF_TGETOBA(NOBTYPMX)
      REAL    DEF_CSCALE_START(NOBTYPMX)
      REAL    DEF_CSCALE_OBTIME(NOBTYPMX)
      REAL    DEF_CSCALE_END(NOBTYPMX)
      REAL    DEF_RADINF(NOBTYPMX)
      REAL    WB_LAT_CC(ROWS_MAX)
      REAL    WB_VERT_V(MODEL_LEVELS_MAX)
      REAL    WB_LAND_FACTOR
      REAL         RADAR_LAT(NRADARS)
      REAL         RADAR_LON(NRADARS)
      REAL         RADAR_RANGE_MAX
      REAL         EPSILON_LHN
      REAL         RELAX_CF_LHN
      REAL         F1_506 , F2_506 , F3_506
      REAL         ALPHA_LHN
      REAL         LHN_LIMIT
      REAL         FI_SCALE_LHN

      REAL    DEF_NUDGE_NH(NOBTYPMX)
      REAL    DEF_NUDGE_TR(NOBTYPMX)
      REAL    DEF_NUDGE_SH(NOBTYPMX)
      REAL    DEF_NUDGE_LAM(NOBTYPMX)

      REAL    DEF_FI_VAR_FACTOR(NOBTYPMX)
      REAL    FI_SCALE
      REAL    FI_SCALE_FACTOR(MODEL_LEVELS_MAX)
      REAL    DF_SCALE
      REAL    DF_SCALE_LEV(MODEL_LEVELS_MAX)
      REAL    DF_COEFF(MODEL_LEVELS_MAX)
      REAL    THRESH_DL
      REAL    THRESH_LM
      REAL    THRESH_MH
      REAL    THRESH_RMSF
      REAL    RADAR_RANGE
      REAL    NORTHLAT, SOUTHLAT, WESTLON, EASTLON
      REAL    VERT_COR_AERO

      LOGICAL LGEO
      LOGICAL LHYDR
      LOGICAL LHYDROL
      LOGICAL LSYN
      LOGICAL LTIMER_AC
      LOGICAL LAC_UARS
      LOGICAL LAC_MES
      LOGICAL LWBAL_SF,     LWBAL_UA
      LOGICAL WB_THETA_UA, WB_LAND_SCALE, WB_THETA_SF
      LOGICAL LRADAR (NRADARS)
      LOGICAL L_LATLON_PRVER
      LOGICAL L_MOPS_EQUALS_RH
      LOGICAL LCHECK_GRID
      LOGICAL      L_506_OBERR
      LOGICAL      L_LHN , L_LHN_SCALE
      LOGICAL      L_LHN_SEARCH , LHN_DIAG
      LOGICAL      L_VERIF_RANGE
      LOGICAL      L_LHN_LIMIT
      LOGICAL      L_LHN_FACT
      LOGICAL      L_LHN_FILT
      LOGICAL L_OBS_CHECK
      LOGICAL  REMOVE_NEG_LH
      LOGICAL USE_CONV_IN_MOPS

      COMMON /COMACP/ NACT,N_GROUPS,NPROG,                              &
     &  AC_OBS_TYPES,     LACT,              GROUP_NO,                  &
     &  LENACT,           LWBAL_SF,          LWBAL_UA,                  &
     &  LTIMER_AC,        LGEO,              LHYDR,                     &
     &  MGEOWT,           LSYN,              LAC_UARS,                  &
     &  OBS_UNITNO,       OBS_FORMAT,        NO_OBS_FILES,              &
     &  L_OBS_CHECK,                                                    &
     &  DIAG_RDOBS,       IUNITNO,           MVINT205,                  &
     &  MHCORFN,          MACDIAG,                                      &
     &  DEF_AC_ORDER,     DEF_NO_ITERATIONS, DEF_INTERVAL_ITER,         &
     &  MWTFN,            MDATADFN,          NSLABS_SCFACT,             &
     &  NO_SCFACT,        NPASS_RF,          MRAMPFN,                   &
     &  IOMITOBS,         TROPINT,           SPEED_LIMIT305,            &
     &  GEOWT_H,          GEOWT_V,           MGLOSSFN,                  &
     &  NON_DIV_COR,      TROPLAT,           VERT_FILT,                 &
     &  NON_DIV_COR_10M,                                                &
     &  VERT_COR_SCALE,                                                 &
     &  VERT_CUTOFF_SL,   VERT_CUTOFF_BW,    VERT_CUTOFF_BH,            &
     &  TIMEF_START,      TIMEF_OBTIME,      TIMEF_END,                 &
     &  CSCFACT_H,        CSCFACT_V,                                    &
     &  MASTER_AC_TYPES,                                                &
     &  DEF_NO_ANAL_LEVS, DEF_NO_WT_LEVS,    DEF_MODE_HANAL,            &
     &  DEF_TIMEB,        DEF_TIMEA,         DEF_TGETOBB,               &
     &  DEF_TGETOBA,      OBTIME_NOM,        DEF_OBTHIN,                &
     &  DEF_RADINF,       DEF_CSCALE_START,  DEF_CSCALE_OBTIME,         &
     &  DEF_CSCALE_END,                                                 &
     &  DEF_NUDGE_NH,     DEF_NUDGE_TR,      DEF_NUDGE_SH,              &
     &  DEF_NUDGE_LAM,                                                  &
     &  WB_LonOffset,     WB_LonPts,         WB_LatOffset,              &
     &  WB_LatPts,                                                      &
     &  GROUP_INDEX,      GROUP_FIRST,       GROUP_LAST,                &
     &  TYPE_INDEX,                                                     &
     &  FI_SCALE,         FI_SCALE_FACTOR,   DEF_FI_VAR_FACTOR,         &
     &  DF_SCALE,         DF_SCALE_LEV,      DF_COEFF,                  &
     &  LAC_MES,                                                        &
     &  THRESH_DL,        THRESH_LM,         THRESH_MH,                 &
     &  THRESH_RMSF,                                                    &
     &  RADAR_RANGE,      LRADAR,            LHYDROL,                   &
     &  L_LATLON_PRVER,   NORTHLAT,          SOUTHLAT,                  &
     &  WESTLON,          EASTLON,           L_MOPS_EQUALS_RH,          &
     &  LHN_RANGE ,  L_LHN , L_LHN_SCALE ,                              &
     &  L_LHN_SEARCH , LHN_DIAG , REMOVE_NEG_LH,                        &
     &  RADAR_LAT , RADAR_LON , RADAR_RANGE_MAX ,                       &
     &  EPSILON_LHN , ALPHA_LHN , RELAX_CF_LHN , LHN_LIMIT ,            &
     &  F1_506 , F2_506 , F3_506 ,                                      &
     &  L_506_OBERR , L_VERIF_RANGE , L_LHN_LIMIT , L_LHN_FACT ,        &
     &  L_LHN_FILT , FI_SCALE_LHN , NPASS_RF_LHN ,                      &
     &  VERT_COR_AERO,    LCHECK_GRID,                                  &
     &  WB_LAT_CC,        WB_VERT_V,         WB_LAND_FACTOR,            &
     & WB_THETA_UA,      WB_LAND_SCALE,   WB_THETA_SF, USE_CONV_IN_MOPS

! COMACP end
! COMACDG

      INTEGER,PARAMETER :: NLDACP=5
      INTEGER,PARAMETER :: NDACOP=200
      INTEGER,PARAMETER :: NDACP=100
      INTEGER,PARAMETER :: NDACVP=11

      LOGICAL :: LDIAGAC
      LOGICAL :: LLDAC(NLDACP)
      LOGICAL :: LLDAG0
      LOGICAL :: LLDAG(NLDACP)
      LOGICAL :: LLBAND
      LOGICAL :: LTITLE
      LOGICAL :: LRMS
      LOGICAL :: LNORMF
      LOGICAL :: LTEMP
      LOGICAL :: LVERIF

      INTEGER :: MDIAGTPS(NOBTYPMX)
      INTEGER :: MDIAG
      INTEGER :: NDGPRT
      INTEGER :: NDACO
      INTEGER :: MDACO(NDACOP)
      INTEGER :: NDAC
      INTEGER :: NDACPRT
      INTEGER :: MDAC(NDACP)
      INTEGER :: NDACV
      INTEGER :: MODACO

      REAL :: DLATN
      REAL :: DLATS
      REAL :: DLONGW
      REAL :: DLONGE
      REAL :: DAGLAT
      REAL :: DAGLON
      REAL :: DACV(NDACP,NDACVP)

      CHARACTER(LEN=10) :: CACV(NDACVP)
      CHARACTER(LEN=12) :: CACT1
      CHARACTER(LEN=32) :: CACT2
      CHARACTER(LEN=30) :: CACT3
      CHARACTER(LEN=55) :: CACT4

      COMMON /COMACDG/ LDIAGAC,LLDAC,LLDAG0,LLDAG,                      &
     &  LRMS,LNORMF,LTEMP,LVERIF,                                       &
     &  NDACO,MDACO,NDAC,NDACPRT,MDAC,NDACV,MODACO,                     &
     &  MDIAGTPS,MDIAG,LLBAND,LTITLE,NDGPRT,                            &
     &  DLATN,DLATS,DLONGW,DLONGE,DAGLAT,DAGLON,                        &
     &  DACV,CACV,CACT1,CACT2,CACT3,CACT4

! COMACDG end
! COMAG start

      INTEGER :: DEF_AGRES_ROWS(NOBTYPMX)
      INTEGER :: DEF_AGRES_PTS(NOBTYPMX)
      INTEGER :: NROWSAG
      INTEGER :: NPTSAGMX
      INTEGER :: NPTSAGMN
      INTEGER :: NPTSAG(ROWS_MAX)
      INTEGER :: NPTS0AG(ROWS_MAX+1)
      INTEGER :: MIN_AGPTS

      REAL :: STAGROW1
      REAL :: STAGPT1
      REAL :: ROW1MG
      REAL :: ROW1MGTH
      REAL :: ROW1MGUV
      REAL :: ROW1AG
      REAL :: DLATAG
      REAL :: DLATMG
      REAL :: PT1MGTH
      REAL :: PT1MGUV
      REAL :: PT1AG
      REAL :: DLONGMG
      REAL :: AGLATDEC
      REAL :: AGROWLEN
      REAL :: DLONGAG(ROWS_MAX)
      REAL :: COSROWAG(ROWS_MAX)

      LOGICAL :: LAGNP
      LOGICAL :: LAGSP

      COMMON /COMAG/ DEF_AGRES_ROWS, DEF_AGRES_PTS,                     &
      ! Analysis grid variables
     &  NROWSAG,                                                        &
     &  NPTSAG,   NPTS0AG,  NPTSAGMX,  NPTSAGMN,                        &
     &  MIN_AGPTS,                                                      &
     &  LAGNP,    LAGSP,                                                &
     &  STAGROW1, STAGPT1,                                              &
     &  ROW1AG,   PT1AG,                                                &
     &  DLATAG,   DLONGAG,                                              &
     &  AGLATDEC, AGROWLEN, COSROWAG,                                   &
      ! Model grid variables
     &  ROW1MG,   ROW1MGTH, ROW1MGUV,                                   &
     &  PT1MGTH,  PT1MGUV,                                              &
     &  DLATMG,   DLONGMG

! COMAG end
!--------------------------------------------------------------------
!LCOMDECK COMMG
!L-------------
      REAL          DLAT,DLONG,XLATN,XLONGW
      REAL          ELFPLAT,ELFPLON
      COMMON/COMMG/ DLAT,DLONG,XLATN,XLONGW,ELFPLAT,ELFPLON
!--------------------------------------------------------------------

!-----DECLARE VARIABLES
      INTEGER       P_FIELD,Q_LEVELS,LENMG
      INTEGER       IROWS,ICOLS,RANGE
      REAL          PR_INC(LENMG),TIMESTEP,                             &
     &              LS_SNOW(P_FIELD),LS_RAIN(P_FIELD),                  &
     &              CONV_SNOW(P_FIELD),CONV_RAIN(P_FIELD),              &
     &              CONV_HEAT(P_FIELD,Q_LEVELS),                        &
     &              LS_HEAT(P_FIELD,Q_LEVELS),                          &
     &              THINCS(P_FIELD,Q_LEVELS)

      REAL lambda_p(1-halo_i:icols+halo_i)
      Real phi_p(1-halo_i:icols+halo_i, 1-halo_j:irows+halo_j)
      Logical L_regular
!*
      INTEGER       ICODE
      CHARACTER(LEN=256) CMESSAGE

!-INTENT=IN---------------------------------------------------
!     P_FIELD                    - No. of points
!     Q_LEVELS                   - No. of wet levels
!     IROWS                      - No. of rows in model grid
!     ICOLS                      - No. of pts in a row
!     RANGE                      - Search range in grid points
!     LENMG                      - Length of model grid
!     PR_INC(LENMG)              - Rainfall incrs on model grid
!     TIMESTEP                   - Timestep in seconds
!     LS_SNOW(P_FIELD)          \                                      .
!     LS_RAIN(P_FIELD)           \  Large scale and convective
!     CONV_SNOW(P_FIELD)         /    rain and snow rates
!     CONV_RAIN(P_FIELD)        /           (diagnostic)
!     CONV_HEAT(P_FIELD,Q_LEVELS)- L.H. incrs to theta due to conv'n
!     LS_HEAT(P_FIELD,Q_LEVELS)  - L.H. incrs to theta due to dynamics
!-INTENT=INOUT-----------------------------------------------
!
!-INTENT=OUT-------------------------------------------------
!     THINCS                     - Calculated increments to theta
!     ICODE                      - Non-zero for failure
!     CMESSAGE                   - Reason for failure
!*-----------------------------------------------------------
!*L External subroutine calls
      EXTERNAL TIMER  

! Local arrays and variables
      INTEGER       SEARCH(4*RANGE*(RANGE+1),2)
                                             ! Search Template
      INTEGER       JPTS , JLEVS             ! Loop counters over model
                                             !       points and levels
      INTEGER       NEAR(P_FIELD)            ! Nearest point found
                                             !    by LHN_SEARCH
      INTEGER       NO_INCRS                 ! Diagnostic - no. of incrs
      INTEGER       NO_SEARCH                !   ""    - no. of searches
      INTEGER       RADIUS(5)                !   "" - breakdown of
                                             !        search results
      INTEGER       M_GRID                   ! Grid type for filter

      REAL          TOT_PR(P_FIELD)          ! Total model precip rate
      REAL          ANAL_PR(P_FIELD)         ! Obs ppn on model grid
      REAL          TOT_LH(P_FIELD,Q_LEVELS) ! Total LH profile
      REAL          LIMIT                    ! Timestep limit of incrs
      REAL          Z                        ! used to set COS_LAT
      REAL          COS_LAT(ROWS_MAX)        ! used in filter routine
      REAL          FI_LHN                   ! filter scale in radians
      LOGICAL       L_FIRST_SEARCH           ! True if first search of
                                             !       the timestep

      REAL    pr_inc_global(Max2DFieldSize)
      REAL    tot_pr_global(Max2DFieldSize)
      REAL    anal_pr_global(Max2DFieldSize)
      INTEGER near_global(Max2DFieldSize)
      REAL    tot_lh_global(glsize(1,fld_type_p)*glsize(2,fld_type_p),  &
     &                      Q_LEVELS/nproc+1)
      REAL    thincs_global(glsize(1,fld_type_p)*glsize(2,fld_type_p),  &
     &                      Q_LEVELS/nproc+1)

      integer row_length_g, p_rows_g, p_field_g, k, istat, iproc
      INTEGER                                                           &
     &  MAP(Q_LEVELS)                                                   &
                            ! processor number for level
     &, N_LEVS_ON_PROC(0:nproc-1)                                       &
                                  ! number of levels on each processor
     &, lev_on_gatpe(Q_LEVELS)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Variables required for variable resolution grids
! GRID_LATS_G are the latitudes for the entire domain. 
! Since the total number of rows is unknown here it is 
! made allocatable (although they could be dimensioned using glsize
! as tot_lh_global and thincs_global
      REAL, ALLOCATABLE :: GRID_LATS_G(:)
      REAL PHI_P_LOCAL(ICOLS,IROWS)
      REAL LAMBDA_P_LOCAL(ICOLS,IROWS)
      REAL, ALLOCATABLE :: PHI_P_GLOBAL(:,:)
      REAL, ALLOCATABLE :: LAMBDA_P_GLOBAL(:,:)
      REAL, ALLOCATABLE :: DELTA_LAT(:)
      REAL, ALLOCATABLE :: DELTA_LON(:)
      INTEGER I,J

!
!  set up global dimensions
!
      IF (lhook) CALL dr_hook('LHN_INC',zhook_in,zhook_handle)
      row_length_g=glsize(1,fld_type_p)
      p_rows_g=glsize(2,fld_type_p)
      p_field_g= row_length_g*p_rows_g


      ALLOCATE (GRID_LATS_G(p_rows_g))

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('LHN_INC ',3)
      IF(LDIAGAC.AND.mype == 0) THEN
        WRITE (6,*) "***** Starting LHN_INC"
      END IF
!
!C 1.0   set parameters and variables
!
      IF (L_LHN_SEARCH .AND. RANGE == 0) THEN
        if(mype == 0)then
          WRITE (6,*) "WARNING : LHN_RANGE=0,"
          WRITE (6,*) "   therefore, setting L_LHN_SEARCH to FALSE"
        endif
        L_LHN_SEARCH = .FALSE.
      ENDIF

! Grid latitudes. If regular grid calculate otherwise obtain from phi_p
      IF(L_regular) THEN

        DO JPTS=1,p_rows_g
          GRID_LATS_G(JPTS)=ROW1MGTH+REAL(JPTS-1)*DLATMG
        ENDDO
 
      ELSE

        ALLOCATE (PHI_P_GLOBAL(row_length_g,p_rows_g))
        ALLOCATE (LAMBDA_P_GLOBAL(row_length_g,p_rows_g))

! Copy information from PHI_P/LAMBDA_P to PHI_P_LOCAL/LAMBDA_P_LOCAL
!  excluding halos
         DO I=1,IROWS
           DO J=1,ICOLS
             PHI_P_LOCAL(J,I)=PHI_P(J,I)
             LAMBDA_P_LOCAL(J,I)=LAMBDA_P(J)
           ENDDO
         ENDDO

! Gather up all PHI_P values into global array and distribute to all PEs
! DEPENDS ON: all_gather_field
        CALL ALL_GATHER_FIELD(PHI_P_LOCAL,   PHI_P_GLOBAL,              &
     &                    ICOLS,         IROWS,                         &
     &                    row_length_g,  p_rows_g,                      &
     &                    fld_type_p,    halo_type_no_halo,             &
     &                    gc_all_proc_group,                            &
     &                    ICODE,         CMESSAGE)

! DEPENDS ON: all_gather_field
        CALL ALL_GATHER_FIELD(LAMBDA_P_LOCAL,   LAMBDA_P_GLOBAL,        &
     &                    ICOLS,         IROWS,                         &
     &                    row_length_g,  p_rows_g,                      &
     &                    fld_type_p,    halo_type_no_halo,             &
     &                    gc_all_proc_group,                            &
     &                    ICODE,         CMESSAGE)
      
! Extract one column of phi values into GRID_LATS_G
! inverting so colats run in correct direction
        DO JPTS=1,P_ROWS_G
          GRID_LATS_G(JPTS)=PHI_P_GLOBAL(1,P_ROWS_G-JPTS+1)/          &
     &                                    PI_OVER_180
           GRID_LATS_G(JPTS)=90.0-GRID_LATS_G(JPTS)
           GRID_LATS_G(JPTS)=GRID_LATS_G(JPTS)*PI_OVER_180
         ENDDO

      ENDIF

! Populate DELTA_LAT and DELTA_LON arrays

      ALLOCATE (DELTA_LAT(p_rows_g))
      ALLOCATE (DELTA_LON(row_length_g))

      IF(L_regular) THEN

! These are not used at present but will be needed if RFCSL is used
! instead of RFCSLR but bit comparibilty will be lost.
!        DELTA_LAT(:)=DLATMG
!        DELTA_LON(:)=DLONGMG

      ELSE

        DO JPTS=1,p_rows_g-1
          DELTA_LAT(JPTS)=PHI_P_GLOBAL(1,JPTS+1)-PHI_P_GLOBAL(1,JPTS)
        ENDDO
        DELTA_LAT(p_rows_g)=DELTA_LAT(p_rows_g-1)
        DO JPTS=1,row_length_g-1
          DELTA_LON(JPTS)=LAMBDA_P_GLOBAL(JPTS+1,1)-                    &
     &                   LAMBDA_P_GLOBAL(JPTS,1)
        ENDDO
        DELTA_LON(row_length_g)=DELTA_LON(row_length_g-1)

      ENDIF  

      IF (L_LHN_LIMIT) LIMIT = LHN_LIMIT * TIMESTEP / 86400.0

      L_FIRST_SEARCH = .TRUE.
      FI_LHN    = FI_SCALE_LHN / Earth_Radius
      NO_INCRS  = 0
      NO_SEARCH = 0
      RADIUS(1) = 0
      RADIUS(2) = 0
      RADIUS(3) = 0
      RADIUS(4) = 0
      RADIUS(5) = 0

!
!C 2     initial calculations
!
!C 2.1   precipitation variables

!
      DO JPTS = 1 , P_FIELD
        TOT_PR(JPTS) = LS_RAIN(JPTS) + LS_SNOW(JPTS)+                   &
     &                          CONV_RAIN(JPTS) + CONV_SNOW(JPTS)
        IF (TOT_PR(JPTS)  <   0.0) TOT_PR(JPTS)=0.0
        ANAL_PR(JPTS) = TOT_PR(JPTS) + PR_INC(JPTS)
        IF (ANAL_PR(JPTS)  <   0.0) THEN
          ANAL_PR(JPTS) = 0.0
          PR_INC(JPTS)   = ANAL_PR(JPTS) - TOT_PR(JPTS)
        ENDIF
      ENDDO     ! JPTS      

!
!C 2.2   heating profiles, (FR Working Paper 171, eq2)
!
      DO JLEVS = 1 , Q_LEVELS
        DO JPTS = 1 , P_FIELD
          TOT_LH(JPTS,JLEVS) = ( CONV_HEAT(JPTS,JLEVS) +                &
     &                     LS_HEAT(JPTS,JLEVS) ) * TIMESTEP
      IF ( REMOVE_NEG_LH .AND. (TOT_LH(JPTS,JLEVS)                      &
     &       <   0.0) ) THEN
          TOT_LH(JPTS,JLEVS) = 0.0
      ENDIF
        ENDDO   ! JPTS
      ENDDO     ! JLEVS

!CCCC
!
!C 3      Calculate theta increments
!
!C 3.1    Set up array NEAR to decide where to get profile for scaling
!         NEAR=-1 implies no scaling necessary, 0 implies scale here
!         greater than 0 means scale a nearby profile, or not at all
!
!  PR_INC, TOT_PR and ANAL_PR are needed on all PEs,
!  so gather them to all pes.

! DEPENDS ON: all_gather_field
      CALL ALL_GATHER_FIELD(PR_INC,       PR_INC_global,                &
     &                  ICOLS,        IROWS,                            &
     &                  row_length_g, p_rows_g,                         &
     &                  fld_type_p,   halo_type_no_halo,                &
     &                  gc_all_proc_group,                              &
     &                  ICODE,        CMESSAGE)


! DEPENDS ON: all_gather_field
      CALL ALL_GATHER_FIELD(TOT_PR,        TOT_PR_global,               &
     &                  ICOLS,         IROWS,                           &
     &                  row_length_g,  p_rows_g,                        &
     &                  fld_type_p,    halo_type_no_halo,               &
     &                  gc_all_proc_group,                              &
     &                  ICODE,         CMESSAGE)


! DEPENDS ON: all_gather_field
      CALL ALL_GATHER_FIELD(ANAL_PR,       ANAL_PR_global,              &
     &                  ICOLS,         IROWS,                           &
     &                  row_length_g,  p_rows_g,                        &
     &                  fld_type_p,    halo_type_no_halo,               &
     &                  gc_all_proc_group,                              &
     &                  ICODE,         CMESSAGE)



! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('LHN_SRCH',3)

!  Calculate NEAR_global on all PEs      
      DO JPTS = 1 , p_field_g
        NEAR_global(JPTS) = -1
        IF ( PR_INC_global(JPTS)  /=  0.0 ) THEN
          NEAR_global(JPTS) = 0
          NO_INCRS   = NO_INCRS + 1
          IF ( TOT_PR_global(JPTS)  <                                   &
     &          (EPSILON_LHN * ANAL_PR_global(JPTS)) ) THEN
            NEAR_global(JPTS) = JPTS
            NO_SEARCH = NO_SEARCH + 1
            IF (L_LHN_SEARCH) THEN
!  search for suitable profile
              CALL LHN_SEARCH(JPTS, NEAR_global(JPTS),                  &
     &           RANGE, SEARCH, p_rows_g, row_length_g,                 &
     &           L_FIRST_SEARCH, TOT_PR_global,                         &
     &           ANAL_PR_global(JPTS), p_field_g, RADIUS,jpts)
            ENDIF
          ENDIF
        ENDIF
      ENDDO    ! JPTS

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('LHN_SRCH',4)
!
!C 3.2   calculate the increments, and scale by the relaxation coeff
!

! Calculate the mapping of which processor each level will go to
      DO K=0,NPROC-1
        N_LEVS_ON_PROC(K)=0
      ENDDO

      DO K=1, Q_LEVELS
       ! assumes first PE is PE 0
        MAP(K)                 = MOD((K-1),NPROC)
        N_LEVS_ON_PROC(MAP(K)) = N_LEVS_ON_PROC(MAP(K))+1
        LEV_ON_GATPE(K)        = N_LEVS_ON_PROC(MAP(K))
      ENDDO

! Distribute TOT_LH_GLOBAL over the processors
      DO K=1, Q_LEVELS
! DEPENDS ON: gather_field
        CALL gather_field( tot_lh(:,k),                                 &
                           tot_lh_global(:,lev_on_gatpe(k)),            &
                           icols, irows,                                &
                           row_length_g, p_rows_g,                      &
                           fld_type_p, halo_type_no_halo,               &
                           map(k), gc_all_proc_group,                   &
                           icode, cmessage )
      END DO


      DO JLEVS = 1 , N_LEVS_ON_PROC(mype)

        DO JPTS = 1 , p_field_g
          IF (NEAR_global(JPTS)  <   0) THEN
!  No scaling necessary
            THINCS_global(JPTS,JLEVS) = 0.0
          ELSEIF (NEAR_global(JPTS)  ==  0) THEN
!  Scale the profile at the point itself
!  (FR WP 171, eq 5)
            THINCS_global(JPTS,JLEVS) =                                 &
     &               RELAX_CF_LHN * TOT_LH_global(JPTS,JLEVS) *         &
     &               PR_INC_global(JPTS) / TOT_PR_global(JPTS)
          ELSEIF (NEAR_global(JPTS)  >   0 .AND.                        &
     &                               NEAR_global(JPTS)  /=  JPTS) THEN
!  Scale a nearby profile
!  (FR WP 171, eq 7)
            THINCS_global(JPTS,JLEVS) = RELAX_CF_LHN *                  &
     &         (ANAL_PR_global(JPTS)/TOT_PR_global(NEAR_global(JPTS)) * &
     &         TOT_LH_global(NEAR_global(JPTS),JLEVS) -                 &
     &         TOT_LH_global(JPTS,JLEVS))
          ELSEIF (NEAR_global(JPTS)  ==  JPTS .AND. L_LHN_SCALE) THEN
!  Scale by 1/EPSILON_LHN if no suitable profile available
!  (FR WP 171, eq 8)
            THINCS_global(JPTS,JLEVS) = RELAX_CF_LHN *                  &
     &         (( 1.0/EPSILON_LHN) - 1.0) * TOT_LH_global(JPTS,JLEVS)
          ELSE
!  If none of the above apply, ignore the ob
            THINCS_global(JPTS,JLEVS) = 0.0
          ENDIF
        ENDDO ! JPTS

      ENDDO ! JLEVS

      IF (.NOT. L_LHN_FILT) THEN
! Copy calculated THINCS values back to individual PEs
        DO K=1, Q_LEVELS
! DEPENDS ON: gather_field
          CALL scatter_field(thincs(:,k),                               &
                             thincs_global(:,lev_on_gatpe(k)),          &
                             icols, irows,                              &
                             row_length_g, p_rows_g,                    &
                             fld_type_p, halo_type_no_halo,             &
                             map(k), gc_all_proc_group,                 &
                             icode, cmessage )
        END DO

!
!  Impose limit by factor of alpha
!
        IF (L_LHN_FACT) THEN
          DO JLEVS = 1 , Q_LEVELS
            DO JPTS = 1 , P_FIELD
              IF (PR_INC(JPTS)  <   0.0) THEN
                IF ((PR_INC(JPTS)/TOT_PR(JPTS)) <  (ALPHA_LHN-1.0)) THEN
                  THINCS(JPTS,JLEVS) = (ALPHA_LHN - 1.0) *              &
     &                                 TOT_LH(JPTS,JLEVS) * RELAX_CF_LHN
                ENDIF
              ENDIF
            ENDDO   ! JPTS
          ENDDO     ! JLEVS
        ENDIF
!
!C 3.3   Recursive filtering of increments
!
      ELSE    ! L_LHN_FILT is true

        IF (L_LHN_FACT) THEN 
          DO JLEVS = 1 , n_levs_on_proc(mype) 
            DO JPTS = 1 , P_FIELD_g 
              IF (PR_INC_global(JPTS)  <   0.0) THEN 
                IF ((PR_INC_global(JPTS)/TOT_PR_global(JPTS)) <         &
                                            (ALPHA_LHN-1.0)) THEN 

                  THINCS_global(JPTS,JLEVS) = (ALPHA_LHN - 1.0) *       & 
                               TOT_LH_global(JPTS,JLEVS) * RELAX_CF_LHN 
                END IF 
              END IF 
            END DO   ! JPTS 
          END DO     ! JLEVS 
        END IF 

!       Initialise COS_LAT
        IF(L_regular) THEN
          Z = ROW1MGTH
          DO JPTS = 1, p_rows_g
            COS_LAT(JPTS) = SIN(Z)
            Z = Z + DLATMG
          ENDDO  ! JPTS
        ELSE
          DO JPTS = 1, p_rows_g
            COS_LAT(JPTS)=SIN(GRID_LATS_G(JPTS))
          ENDDO  ! JPTS
        ENDIF


        M_GRID = 0

      IF(L_regular) THEN
        DO JLEVS = 1 , N_LEVS_ON_PROC(mype)
        CALL RFCSLR(THINCS_global(1,JLEVS),p_rows_g,row_length_g,M_GRID,&
     &             0.0,COS_LAT,DLATMG,DLONGMG,FI_LHN,NPASS_RF_LHN)
        ENDDO

      ELSE
                
        DO JLEVS = 1 , N_LEVS_ON_PROC(mype)
        CALL RFCSL(THINCS_global(1,JLEVS),p_rows_g,row_length_g,M_GRID, &
     &             0.0,COS_LAT,DELTA_LAT,DELTA_LON,FI_LHN,NPASS_RF_LHN)
        ENDDO

      ENDIF

! Copy calculated THINCS values back to individual PEs
      DO K=1, Q_LEVELS
! DEPENDS ON: gather_field
        CALL scatter_field(thincs(:,k),                                 &
                           thincs_global(:,lev_on_gatpe(k)),            &
                           icols, irows,                                &
                           row_length_g, p_rows_g,                      &
                           fld_type_p, halo_type_no_halo,               &
                           map(k), gc_all_proc_group,                   &
                           icode, cmessage )
      END DO

      ENDIF   ! L_LHN_FILT


!
!C Impose absolute limit
!
          IF (L_LHN_LIMIT) THEN
            DO JLEVS = 1 , Q_LEVELS
              DO JPTS = 1 , P_FIELD
                IF (THINCS(JPTS,JLEVS)  <   -1.0 * LIMIT)               &
     &                 THINCS(JPTS,JLEVS) = -1.0 * LIMIT
                IF (THINCS(JPTS,JLEVS)  >   LIMIT)                      &
     &                 THINCS(JPTS,JLEVS) = LIMIT
              ENDDO  ! JPTS
            ENDDO    ! JLEVS
          ENDIF

!
!C  4  Diagnostics
!
      IF (LDIAGAC .AND. LHN_DIAG) THEN
       if(mype == 0)then
        WRITE (6,*) ' '
        WRITE (6,*) "Latent Heat Nudging Scheme, LHN_INC, diagnostics"
        WRITE (6,*) "    Parameters set   : EPSILON_LHN = ",EPSILON_LHN
        WRITE (6,*) "                     : LHN_RANGE   = ",RANGE
        IF (L_LHN_FACT) THEN
          WRITE (6,*) "    Limit increments to scale down rain rate",   &
     &                " by at most factor of 1/ ALPHA"
          WRITE (6,*) "       ALPHA=",ALPHA_LHN
        ENDIF
        IF (L_LHN_FILT) THEN
          WRITE (6,*) "    Filtering of increments performed"
          WRITE (6,*) "       filter scale = ",FI_SCALE_LHN             &
     &                                                  / 1000.0," Km"
        ENDIF
        IF (L_LHN_LIMIT) THEN
          WRITE (6,*) "    Limiting of increments set to ",LIMIT,       &
     &                " degrees per timestep = ",LHN_LIMIT," degrees",  &
     &                                                     " per day"
        ENDIF
        WRITE (6,*) ' '
        WRITE (6,*) "Number of increments required was ",NO_INCRS
        IF (L_LHN_SEARCH) THEN
          RADIUS(1) = RADIUS(2) + RADIUS(3) + RADIUS(4)
          WRITE (6,*) "LHN_SEARCH was called for ",NO_SEARCH," of them"
          WRITE (6,*) "The search was successful ",RADIUS(1)," times"
          WRITE (6,*) "It found :"
          WRITE (6,*) "     ",RADIUS(2)," pts at 1 point away"
          WRITE (6,*) "     ",RADIUS(3)," pts at 2 points away"
          WRITE (6,*) "     ",RADIUS(4)," pts at 3 or more points away"
          WRITE (6,*) "     ",RADIUS(5)," total pts searched"
        ELSE
          WRITE (6,*) NO_SEARCH," points required the search, but"
          WRITE (6,*) " the search algorithm was disabled"
        ENDIF
        WRITE (6,*) "Scaling at points failing the search was set to: ",&
     &                        L_LHN_SCALE
        WRITE (6,*) ' '
       endif
      ENDIF

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('LHN_INC ',4)

      IF(ALLOCATED(GRID_LATS_G)) DEALLOCATE(GRID_LATS_G)
      IF(ALLOCATED(PHI_P_GLOBAL)) DEALLOCATE (PHI_P_GLOBAL)
      IF(ALLOCATED(LAMBDA_P_GLOBAL)) DEALLOCATE (LAMBDA_P_GLOBAL)
      IF(ALLOCATED(DELTA_LAT)) DEALLOCATE (DELTA_LAT)
      IF(ALLOCATED(DELTA_LON)) DEALLOCATE (DELTA_LON)

      IF (lhook) CALL dr_hook('LHN_INC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LHN_INC
END MODULE lhn_inc_mod
