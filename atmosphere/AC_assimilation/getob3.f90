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
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: AC Assimilation
MODULE getob3_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE GETOB3 (KACT,OBS,OBDATA,                               &
     &                   OBS_LAT,OBS_LONG,                              &
     &                   LENOBT,INOBS,NDV,OBS_NO,ICODE,CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Atmos_Max_Sizes
      USE UM_ParParams
      USE comobs_mod, ONLY: nobtypmx, obs_info
      IMPLICIT NONE

      
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
!--------------------------------------------------------------------
!LCOMDECK COMMG
!L-------------
      REAL          DLAT,DLONG,XLATN,XLONGW
      REAL          ELFPLAT,ELFPLON
      COMMON/COMMG/ DLAT,DLONG,XLATN,XLONGW,ELFPLAT,ELFPLON
!--------------------------------------------------------------------
!-----------------------------------------------------------------------
      INTEGER KACT                   ! IN  Index to obs type
      INTEGER LENOBT                 ! IN  No of obs to be assimilated
      INTEGER INOBS                  ! IN  No of obs for this type
      INTEGER NDV                    ! IN  No of data values
!                                    !     excluding header section
      REAL OBS (INOBS,*)             ! IN  OBS array from RDOBS
      REAL OBDATA (LENOBT,NDV)       ! OUT ob values and error ratios
      REAL OBS_LAT (LENOBT)          ! IN  latitudes
      REAL OBS_LONG(LENOBT)          ! IN  longitudes
      INTEGER OBS_NO(LENOBT)         ! IN  pointers to obs.
      INTEGER ICODE                  ! OUT error code and message
      CHARACTER(LEN=256) CMESSAGE
!-----------------------------------------------------------------------
!L local work arrays
      REAL WLT    (LENOBT)
      REAL WLN    (LENOBT)
      REAL WLN2   (LENOBT)
      REAL COEFF1 (LENOBT)
      REAL COEFF2 (LENOBT)
      REAL UWK    (LENOBT)
      REAL VWK    (LENOBT)
!-----------------------------------------------------------------------
      REAL ZMINEP
      PARAMETER(ZMINEP=1.E-10)
!     ZMINEP IS THE MINIMUM ALLOWED OBSERVATIONAL ERROR TO AVOID /0.
      CHARACTER(LEN=10) LABEL
      INTEGER NEROPT,KTYPE,NLEV,JDV,JOB,KDV

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
      EXTERNAL TIMER
      EXTERNAL EQTOLL,W_COEFF,W_EQTOLL
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GETOB3',zhook_in,zhook_handle)

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOB3  ',3)

      KTYPE  = LACT (KACT)                !  AC Observation Type
      NLEV   = OBS_INFO % NOBLEV(KACT)    !  No of levels for this type

!     Get observation data and error ratios
      DO JDV=1,NDV
        DO JOB=1,LENOBT
          OBDATA(JOB,JDV) = OBS(OBS_NO(JOB)-OBS_INFO % OBS_NO_ST(KACT), &
                                OBS_INFO % NDVHDR+JDV)
        ENDDO
      ENDDO

      IF (LDIAGAC) THEN


        IF(NDAC >  0)THEN

!      DO LAT/LON TRANSFORMATIONS AND GET COEFF1,COEFF2 OUTSIDE JDV LOOP

         DO JDV=1,NDV
         IF (KTYPE == 101) THEN
                 IF (JDV == 1) LABEL = 'P* DATA'
                 IF (JDV == 2) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 201) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'T  DATA',JDV
                 ELSE
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ENDIF
         ELSEIF (KTYPE == 202 .OR. KTYPE == 204) THEN
               IF (JDV == 1) LABEL = 'T DATA'
               IF (JDV == 2) LABEL = 'EO/EP'
         ELSE IF(KTYPE == 203) THEN
               IF (JDV == 1) LABEL = 'OBS LEV'
               IF (JDV == 2) LABEL = 'T DATA'
               IF (JDV == 3) LABEL = 'EO/EP'
         ELSE IF(KTYPE == 205.OR.KTYPE == 206.OR.                       &
     &           KTYPE == 207.OR.KTYPE == 209.OR.                       &
     &           KTYPE == 211) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'T  DATA',JDV
                 ELSE
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ENDIF
         ELSEIF (KTYPE == 208) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'T  DATA',JDV
                 ELSEIF (JDV >  NLEV.AND.JDV <= 2*NLEV) THEN
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'T F/G  ',KDV
                 ELSEIF (JDV >  2*NLEV.AND.JDV <= 3*NLEV) THEN
                    KDV=JDV-NLEV*2
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ELSE
                    LABEL='QUALIND'
                 ENDIF
         ELSEIF (KTYPE == 301 .OR. KTYPE == 311) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'U  DATA',JDV
                 ELSEIF (JDV <= NLEV*2)THEN
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'V  DATA',KDV
                 ELSE
                    KDV=JDV-NLEV*2
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ENDIF
         ELSEIF (KTYPE == 302 .OR. KTYPE == 304 .OR.                    &
     &           KTYPE == 305 .OR. KTYPE == 306 ) THEN
                 IF (JDV == 1) LABEL = 'U DATA'
                 IF (JDV == 2) LABEL = 'V DATA'
                 IF (JDV == 3) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 303) THEN
                 IF (JDV == 1) LABEL = 'OBS LEV'
                 IF (JDV == 2) LABEL = 'U DATA'
                 IF (JDV == 3) LABEL = 'V DATA'
                 IF (JDV == 4) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 401 .OR. KTYPE == 405 .OR. KTYPE == 406) THEN
                 IF (JDV <= NLEV) THEN
                    WRITE(LABEL,21)'RH DATA',JDV
                 ELSE
                    KDV=JDV-NLEV
                    WRITE(LABEL,21)'EO/EP  ',KDV
                 ENDIF
         ELSEIF (KTYPE == 402 .OR. KTYPE == 404) THEN
                 IF (JDV == 1) LABEL = 'RH DATA'
                 IF (JDV == 2) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 403) THEN
                 IF (JDV == 1) LABEL = 'OBS LEV'
                 IF (JDV == 2) LABEL = 'RH DATA'
                 IF (JDV == 3) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 407) THEN
                 WRITE(LABEL,21)'CTT CNT',JDV
         ELSEIF (KTYPE == 506) THEN
                 IF (JDV == 1) LABEL = 'PR DATA'
                 IF (JDV == 2) LABEL = 'EO/EP'
         ELSEIF (KTYPE == 901) THEN
                 IF (JDV == 1) LABEL = 'VIS DATA'
                 IF (JDV == 2) LABEL = 'EO/EP'
         ELSE
         ICODE=1
         CMESSAGE=' GETOB3 : ILLEGAL OB TYPE'
         GOTO 999
         ENDIF
!
   21 FORMAT(A7,I3)
!
       ENDDO   !  Loop over JDV
       END IF
      END IF


      NEROPT = OBS_INFO % NERLEV1(KACT) - OBS_INFO % NDVHDR

      DO JDV=NEROPT,NEROPT+NLEV-1

!       ENSURE ALL OB ERRORS ARE >= ZMINEP TO AVOID
!       DIVIDING BY ZERO IN VAN## Routines - is this necessary?

        DO JOB=1,LENOBT
          IF (OBDATA(JOB,JDV) <  ZMINEP .AND. OBDATA(JOB,JDV) /=        &
              OBS_INFO % MISSD) OBDATA(JOB,JDV) = ZMINEP
        ENDDO

      ENDDO

999   CONTINUE
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOB3  ',4)
      IF (lhook) CALL dr_hook('GETOB3',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GETOB3
END MODULE getob3_mod
