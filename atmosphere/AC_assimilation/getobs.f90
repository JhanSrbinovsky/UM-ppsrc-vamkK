! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE GETOBS,GETOB2,GETOB3----------------------------------
!LL
!LL  Purpose : Extract appropriate observations from the COMOBS common
!LL            block and the OBS array (passed from RDOBS via argument
!LL            list)
!LL            GETOBS called from AC gets list of obs in time window
!LL            GETOB2 called from AC2 gets lat, long, time and
!LL                   model obs type for an obs type.
!LL            GETOB3 called from VERTANL gets data specific to a type
!LL                   (eg data values and assoc error ratios)
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: AC Assimilation
MODULE getobs_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE GETOBS (KACT,ASSM_TIME,INOBS,OBS_NO,LENOBT,            &
     &                   OBS,OBS_FLAG,                                  &
     &                   TIMESTEP_NO,DG_ONLY,                           &
     &                   INOBSDIM,ICODE,CMESSAGE)
!L
!L    --------------------------------------------------------
!L  CALLED BY SUBROUTINE AC
!L  THIS OBTAINS ALL OBSERVATIONS TOGETHER WITH THEIR POSITION
!L  OF TYPE LACT(KACT) WHICH ARE WITHIN THE INSERTION PERIOD RELATIVE
!L  TO THE ASSIMILATION TIME
!L  (for scatwinds and sat120 a subset may be chosen)
!L

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx, obs_info
      IMPLICIT NONE


!     UM AC Comdecks
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
      INTEGER Iproc
!-----------------------------------------------------------------------
      INTEGER KACT      ! IN Index for this obs type (pos in LACT etc)
      INTEGER INOBS           ! IN  Total no of obs
      INTEGER INOBSDIM        ! IN  Total no of obs (for dimensioning)
      INTEGER OBS_NO(INOBSDIM)   ! IN  Pointers to obs to be assimilated
      INTEGER OBS_FLAG(INOBSDIM) ! IN  Observation flags
      INTEGER TIMESTEP_NO     ! IN  Timestep number
      INTEGER LENOBT          ! OUT No of obs to be assimilated

      REAL    OBS(INOBSDIM,*) ! IN  Observation data for this type
      REAL    ASSM_TIME       ! IN  Assm time relative to start

      LOGICAL DG_ONLY         ! IN  Diagnostic only switch
      INTEGER ICODE           ! IN  Return code
      CHARACTER(LEN=256) CMESSAGE  ! IN  Error message
!-----------------------------------------------------------------------
!     Local variables
      INTEGER INB       !  No of obs before insertion period
      INTEGER INA       !  No of obs after insertion period
      INTEGER INF       !  No of obs that have been flagged
      INTEGER INT       !  No of obs thinned out (skipped)
      INTEGER JOB       !  Loop counter over no of obs
      INTEGER ISTART    !  First obs
      INTEGER ISKIP     !  Skip ISKIP obs in loop over obs
      INTEGER IP_TIME   !  Pointer to times in observation data
      INTEGER OBTHIN    !  Thinning factor for this ob type
      INTEGER ISTAT     !  for calls to GCOM routines
      REAL    TDIFF     !  Time difference between obs and assm time
      REAL    TGETOBB   !  Insertion period before obs time
      REAL    TGETOBA   !  Insertion period after obs time

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
      EXTERNAL TIMER
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GETOBS',zhook_in,zhook_handle)
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOBS  ',3)

      INB=0
      INA=0
      INF=0
      INT=0

! SET UP INDEX VECTOR FOR GATHERING OBS
!  ALLOWING FOR OPTION TO THIN OBS :likely use GLOSS/SAT120/ERS-1 DATA
!  BUT NOT ON DIAGNOSTIC ONLY ITERATION
!     Get OBTHIN for this observation type
      OBTHIN = DEF_OBTHIN( TYPE_INDEX(KACT) )
      IF (DG_ONLY) THEN
        ISTART = 1
        ISKIP  = 1
      ELSE
        IF (OBTHIN >  1 .AND. INOBS >  OBTHIN) THEN
          ISTART = 1 + MOD(TIMESTEP_NO,OBTHIN)
          ISKIP  = OBTHIN
        ELSE
          ISTART = 1
          ISKIP  = 1
        ENDIF
      ENDIF

!     Get TGETOBB and TGETOBA for this observation type
      TGETOBB = DEF_TGETOBB( TYPE_INDEX(KACT) )
      TGETOBA = DEF_TGETOBA( TYPE_INDEX(KACT) )

!     Pointer to observation time for this obs type
      IP_TIME = 3

      IF(INOBS >  0)THEN


      DO JOB = ISTART,INOBS,ISKIP

!     Get time difference between obs time and current assm time
      TDIFF = OBS(JOB,IP_TIME)-ASSM_TIME

!     Use obs if : time diff <  tgetobb-0.5
!         and if : time diff > -(tgetoba-0.5)
!         and if : it has not been flagged
!   -0.5 minutes helps make results same on different machines

      IF (TDIFF >= TGETOBB-0.5) THEN    !  Obs yet to be assimilated
        INB=INB+1
      ELSEIF (TDIFF <= -TGETOBA+0.5) THEN !  Obs has been assimilated
        INA=INA+1
      ELSEIF (OBS_FLAG(JOB) /= 0) THEN  !  Obs has been flagged
        INF=INF+1
      ELSE                              !  Assimilate this observation
        LENOBT=LENOBT+1
        OBS_NO(LENOBT) = OBS_INFO % OBS_NO_ST(KACT)+JOB
      END IF
!
      ENDDO

      ENDIF  !INOBS >  0

      INT=INOBS-(LENOBT+INB+INA+INF)
!
      IF (LDIAGAC) THEN
        CountA(1)=LENOBT
        CountA(2)=INB
        CountA(3)=INA
        CountA(4)=INF
        CountA(5)=INT
        If(mype == 0)Then
          Do JOB=1,5
           CountC(JOB)=CountA(JOB)
          EndDo
        Endif


         Do iproc=1,nproc-1
          IF(mype == 0) THEN
! PE0 receives
            CALL GC_IRECV(IPROC,5,IPROC,ISTAT,COUNTB,COUNTA)
            Do JOB=1,5
              CountC(JOB)=CountC(JOB)+CountB(JOB)
             EndDo
           ELSEIF(mype == IPROC) THEN
! other PEs send
             CALL GC_ISEND(IPROC,5,0,ISTAT,COUNTB,COUNTA)
           ENDIF
         EndDo


         PRINT '(A,I4,A,I6,A,I6,A,I6,A,I6,A,I6,A)',                     &
     &   ' TYPE',LACT(KACT),'  LENOBT',COUNTC(1),'  OMITTED WERE :',    &
     &   COUNTC(2),' (BEFORE) ',COUNTC(3),' (AFTER) ',                  &
     &   COUNTC(4),' (FLAGGED) ',COUNTC(5),' (THINNED) '

      ENDIF
!
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOBS  ',4)
      IF (lhook) CALL dr_hook('GETOBS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GETOBS
!
END MODULE getobs_mod
