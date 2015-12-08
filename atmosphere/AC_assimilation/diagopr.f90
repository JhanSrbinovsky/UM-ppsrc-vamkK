! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DIAGOPR ------------------------------------------------
!LL
!LL  Purpose : Provide Statistics on Precipitation Rates
!LL
!LL            Calculate and print out the following :
!LL            - Contingency tables of Observation vs Forecast data
!LL            - Various skill scores
!LL
!LLEND------------------------------------------------------------------

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE diagopr_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE DIAGOPR  (lenobt_dim,                                  &
     &                     INTPR,OBDATA,LENOBT_local,LENOB,NDV,         &
     &                     LMISSD,NOBIPT,NPTOBT,NO_ANAL_LEVS)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE

      integer lenobt_dim
      INTEGER LENOBT_local,LENOB,NDV,NOBIPT,NPTOBT,NO_ANAL_LEVS
      REAL    INTPR  (lenobt_dim),OBDATA (lenobt_dim,NDV)
      LOGICAL LMISSD (LENOB+1,NO_ANAL_LEVS)

!   INTENT=INOUT-----------------------------------------------
!   INTPR        - interpolated precipitation rates
!   OBDATA       - observation data
!   LENOBT       - number of observations of this type
!   NDV          - number of data values
!   NOBIPT       - pointer to data in OBDATA
!   LMISSD       - logical array to point to obs with missing data
!   LENOB        - number of obs in group of ob types
!   NO_ANAL_LEVS - no of analysis levels
!   NPTOBT       - pointer to present ob type within group
!
!   INTENT=OUT-------------------------------------------------
!* ------------------------------------------------------------


!-----AC common blocks
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


!*L   External subroutine calls
      EXTERNAL TIMER

!     ----------------------------------------------------
!     Local arrays and variables
      integer LENOBT
      integer iproc

      INTEGER OBCAT,FCAT,JOB,JCAT,JCAT1,                                &
     &        MISSING,OBS_USED,                                         &
     &        JBDY,COUNTER,COUNTERF                                     &
     &        ,COUNTERF2
      REAL    AVG_FCRATE,AVG_OBRATE,                                    &
     &        OBS_PCNT,FC_PCNT,CRCT_PCNT,                               &
     &        RMS,RMSF,CC,                                              &
     &        RMSF2,                                                    &
     &        HR(3),FAR(3),BIAS_RATE,BIAS_AREA,                         &
     &        HKS(3),TS(3),                                             &
     &        TAB_RATE(5,5),TAB_YN(2,2)
      INTEGER                                                           &
     &        MISSING_rem(0:NPROC_MAX),                                   &
     &        COUNTER_rem(0:NPROC_MAX),COUNTERF_rem(0:NPROC_MAX)            &
     &        ,COUNTERF2_rem(0:NPROC_MAX)                                 &
     &        ,lenobt_rem(0:NPROC_MAX)
      INTEGER ISTAT  ! for calls to GCOM
      REAL    AVG_FCRATE_rem(0:NPROC_MAX),AVG_OBRATE_rem(0:NPROC_MAX),      &
     &        RMS_rem(NPROC_MAX),RMSF_rem(0:NPROC_MAX),                     &
     &        RMSF2_rem(0:NPROC_MAX),                                     &
     &        TAB_RATE_rem(5,5,0:NPROC_MAX)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
!   TAB_RATE   - precipitation magnitude contingency table
!   TAB_YN     - 2x2 'yes/no' table compiled from TAB_RATE
!   OBCAT,FCAT - pointers within TAB_RATE, TAB_PHASE
!   OBS_USED   - total number of observations used
!   MISSING    - counter for no of obs with missing data
!   AVG_OBRATE - mean observed precipitation rate
!   AVG_FCRATE - mean forecast precipitation rate at obs points
!   OBS_PCNT   - percentage of observations which give precipitation
!   FC_PCNT    - percentage of forecast values which give precipitation
!   CRCT_PCNT  - percentage of correct forecasts
!   RMS        - Root mean square error of forecast rate v obs
!   RMSF       - Root mean square factor
!   RMSF2      - Root mean square factor, for either value > threshold
!   CC         - Tetrachoric Correlation Coefficient
!   HR         - Hit Rate               }  arrays
!   FAR        - False Alarm Rate       }  with values based
!   HKS        - Hanssen & Kuiper Score }  on light, moderate
!   TS         - Threat Score           }  and heavy thresholds
!   BIAS_RATE  - Bias Ratio for precipitation rates
!   BIAS_AREA  - Bias Ratio for number of points with rain
!   JOB        - loop counter over obs
!   JCAT    }  - loop counters
!   JCAT1   }  -               over categories of obs/fc precip
!   JBDY       - loop counter for precip rate category boundaries
!   COUNTER    - counter for obs contributing to rms rate error
!   COUNTERF   - counter for obs contributing to rms factor score
!   COUNTERF2  - counter for obs contributing to RMSF2 score
!

      IF (lhook) CALL dr_hook('DIAGOPR',zhook_in,zhook_handle)

      if(mype == 0)then
      PRINT '('' '')'
      PRINT '(''  PRECIPITATION verification routine DIAGOPR'')'
      PRINT '('' '')'
      PRINT '(''  Called from VANRAIN'')'
      IF (L_VERIF_RANGE) THEN
      PRINT  '(''Verification range is '',F12.6,'' km from nearest radar'')',&
          RADAR_RANGE
      ELSE
      PRINT  '(''Verification range is '',F12.6,'' km from nearest radar'')',&
          RADAR_RANGE_MAX
      ENDIF
      PRINT '('' '')'
      endif

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('DIAGOPR ',3)

!L*** 1.0  Initialize variables
      DO JCAT=1,5
        DO JCAT1=1,5
          TAB_RATE(JCAT,JCAT1) = 0
        END DO
      END DO
      LENOBT=0
      MISSING=0
      AVG_FCRATE=0.0
      AVG_OBRATE=0.0
      COUNTER=0
      COUNTERF=0
      COUNTERF2=0
      RMSF2=0.0
      RMS=0.0
      RMSF=0.0
      CC =9999999.9

!
!L***2.0   Calculate average forecast and observed rates,
!L         rms error in fc rate and rms factor.
!          Ignore points where ob is missing data.
!
      DO JOB=1,LENOBT_local
         IF (LMISSD(NPTOBT+JOB,1)) THEN
            MISSING = MISSING + 1
         ELSE
            AVG_FCRATE = AVG_FCRATE + INTPR(JOB)
            AVG_OBRATE = AVG_OBRATE + OBDATA(JOB,NOBIPT)
            IF ( INTPR(JOB)  >   THRESH_DL .AND.                        &
     &           OBDATA(JOB,NOBIPT)  >   THRESH_DL ) THEN
               COUNTER = COUNTER + 1
              RMS = RMS + ( INTPR(JOB) - OBDATA(JOB,NOBIPT) ) ** 2
            ENDIF
            IF ( INTPR(JOB)  >   THRESH_RMSF .AND.                      &
     &           OBDATA(JOB,NOBIPT)  >   THRESH_RMSF ) THEN
               COUNTERF = COUNTERF + 1
               RMSF=RMSF+(LOG(ABS( INTPR(JOB)/OBDATA(JOB,NOBIPT) )))**2
            ENDIF
      IF ( INTPR(JOB)  >= THRESH_RMSF .OR.                              &
     &                 OBDATA(JOB,NOBIPT)  >= THRESH_RMSF ) THEN
        COUNTERF2 = COUNTERF2 + 1
        RMSF2=RMSF2 + (LOG( MAX(THRESH_RMSF/2.0,INTPR(JOB))             &
     &            / MAX(THRESH_RMSF/2.0,OBDATA(JOB,NOBIPT)) ))**2

      ENDIF
         ENDIF
         END DO ! Job

!
!L*** 3.0   calculate 4x4 table for rain magnitudes
!
      DO JOB = 1,LENOBT_local
         IF ( .NOT. LMISSD(NPTOBT+JOB,1) ) THEN
            FCAT = 1
            IF ( INTPR(JOB)  >   THRESH_DL ) FCAT=2
            IF ( INTPR(JOB)  >   THRESH_LM ) FCAT=3
            IF ( INTPR(JOB)  >   THRESH_MH ) FCAT=4
            OBCAT = 1
            IF ( OBDATA(JOB,NOBIPT)  >   THRESH_DL ) OBCAT=2
            IF ( OBDATA(JOB,NOBIPT)  >   THRESH_LM ) OBCAT=3
            IF ( OBDATA(JOB,NOBIPT)  >   THRESH_MH ) OBCAT=4
            TAB_RATE(OBCAT,FCAT) = TAB_RATE(OBCAT,FCAT) + 1
         ENDIF
       END DO ! Job

      IF(mype == 0) THEN
        lenobt_rem(mype)=lenobt_local
        missing_rem(mype)=missing
        counter_rem(mype)=counter
        counterf_rem(mype)=counterf
        counterf2_rem(mype)=counterf2
        rms_rem(mype)=rms
        rmsf_rem(mype)=rmsf
        rmsf2_rem(mype)=rmsf2
        AVG_FCRATE_rem(mype)=AVG_FCRATE
        AVG_OBRATE_rem(mype)=AVG_OBRATE
        JBDY=0
        DO JCAT=1,5
          DO JCAT1=1,5
            JBDY=JBDY+1
            tab_rate_rem(1,1,JBDY)=tab_rate(JCAT1,JCAT)
          ENDDO
        ENDDO
      ENDIF

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_IRECV(IPROC,1,IPROC,ISTAT,LENOBT_REM(IPROC),          &
     &                  LENOBT_LOCAL)
         ELSEIF(mype == IPROC) THEN
           CALL GC_ISEND(IPROC,1,0,ISTAT,LENOBT_REM(IPROC),             &
     &                  LENOBT_LOCAL)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_IRECV(IPROC,1,IPROC,ISTAT,MISSING_REM(IPROC),         &
     &                  MISSING)
         ELSEIF(mype == IPROC) THEN
           CALL GC_ISEND(IPROC,1,0,ISTAT,MISSING_REM(IPROC),            &
     &                  MISSING)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_IRECV(IPROC,1,IPROC,ISTAT,COUNTER_REM(IPROC),         &
     &                  COUNTER)
         ELSEIF(mype == IPROC) THEN
           CALL GC_ISEND(IPROC,1,0,ISTAT,COUNTER_REM(IPROC),            &
     &                  COUNTER)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_IRECV(IPROC,1,IPROC,ISTAT,COUNTERF_REM(IPROC),        &
     &                  COUNTERF)
         ELSEIF(mype == IPROC) THEN
           CALL GC_ISEND(IPROC,1,0,ISTAT,COUNTERF_REM(IPROC),           &
     &                  COUNTERF)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_IRECV(IPROC,1,IPROC,ISTAT,COUNTERF2_REM(IPROC),       &
     &                  COUNTERF2)
         ELSEIF(mype == IPROC) THEN
           CALL GC_ISEND(IPROC,1,0,ISTAT,COUNTERF2_REM(IPROC),          &
     &                  COUNTERF2)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_RRECV(IPROC,1,IPROC,ISTAT,RMS_REM(IPROC),RMS)
         ELSEIF(mype == IPROC) THEN
           CALL GC_RSEND(IPROC,1,0,ISTAT,RMS_REM(IPROC),RMS)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_RRECV(IPROC,1,IPROC,ISTAT,RMSF_REM(IPROC),RMSF)
         ELSEIF(mype == IPROC) THEN
           CALL GC_RSEND(IPROC,1,0,ISTAT,RMSF_REM(IPROC),RMSF)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_RRECV(IPROC,1,IPROC,ISTAT,RMSF2_REM(IPROC),RMSF2)
         ELSEIF(mype == IPROC) THEN
           CALL GC_RSEND(IPROC,1,0,ISTAT,RMSF2_REM(IPROC),RMSF2)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_RRECV(IPROC,1,IPROC,ISTAT,AVG_FCRATE_REM(IPROC)       &
     &                 ,AVG_FCRATE)
         ELSEIF(mype == IPROC) THEN
           CALL GC_RSEND(IPROC,1,0,ISTAT,AVG_FCRATE_REM(IPROC)          &
     &                  ,AVG_FCRATE)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_RRECV(IPROC,1,IPROC,ISTAT,AVG_OBRATE_REM(IPROC)       &
     &                 ,AVG_OBRATE)
         ELSEIF(mype == IPROC) THEN
           CALL GC_RSEND(IPROC,1,0,ISTAT,AVG_OBRATE_REM(IPROC)          &
     &                  ,AVG_OBRATE)
         ENDIF
      ENDDO

      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN
          CALL GC_RRECV(IPROC,25,IPROC,ISTAT,TAB_RATE_REM(1,1,IPROC)    &
     &                 ,TAB_RATE)
         ELSEIF(mype == IPROC) THEN
           CALL GC_RSEND(IPROC,25,0,ISTAT,TAB_RATE_REM                  &
     &                  ,TAB_RATE)
         ENDIF
      ENDDO

      if(mype == 0)then

      lenobt=lenobt_local

      do iproc=1,nproc-1

      lenobt=lenobt+lenobt_rem(iproc)
      missing=missing+missing_rem(iproc)
      counter=counter+counter_rem(iproc)
      counterf=counterf+counterf_rem(iproc)
      counterf2=counterf2+counterf2_rem(iproc)
      rms=rms+rms_rem(iproc)
      rmsf=rmsf+rmsf_rem(iproc)
      rmsf2=rmsf2+rmsf2_rem(iproc)
      AVG_FCRATE=AVG_FCRATE+AVG_FCRATE_rem(iproc)
      AVG_OBRATE=AVG_OBRATE+AVG_OBRATE_rem(iproc)

      do jcat=1,4
      do jcat1=1,4
      tab_rate(jcat,jcat1)= tab_rate(jcat,jcat1)                        &
     &                +tab_rate_rem(jcat,jcat1,iproc)
      enddo
      enddo

      enddo !iproc=0,nproc-1

      OBS_USED = LENOBT - MISSING
      IF (OBS_USED == 0) THEN
         PRINT '(''  *** all obs are missing data ***'')'
         GOTO 9999
      ENDIF

! calculate average rates, rms rate error and rms factor
! convert from model units of  kg m-2 s-1  to  mm/hr
!
      AVG_FCRATE = AVG_FCRATE / OBS_USED * 3600.0
      AVG_OBRATE = AVG_OBRATE / OBS_USED * 3600.0

      IF (COUNTER  ==  0) THEN
         RMS = 99999999.
      ELSE
         RMS = SQRT(RMS/COUNTER)
      ENDIF
      IF (COUNTERF  ==  0) THEN
         RMSF= 99999999.
      ELSE
         RMSF = EXP( SQRT(RMSF/COUNTERF) )
      ENDIF
      IF (COUNTERF2  ==  0) THEN
         RMSF2 = 99999999.
      ELSE
         RMSF2 = EXP( SQRT(RMSF2/COUNTERF2) )
      ENDIF


!
!   Convert tables to percentages, and calculate row totals
!
      DO JCAT=1,4
        DO JCAT1=1,4
          TAB_RATE(JCAT,JCAT1) = TAB_RATE(JCAT,JCAT1)*100.0/OBS_USED
          TAB_RATE(5,JCAT1) = TAB_RATE(5,JCAT1) + TAB_RATE(JCAT,JCAT1)
          TAB_RATE(JCAT,5) = TAB_RATE(JCAT,5) + TAB_RATE(JCAT,JCAT1)
        END DO
      END DO

!
!L*** 4.0  Condense 4x4 TAB_RATE into various 2x2 tables, TAB_YN
!
      DO JBDY=1,3
!
!  Calculate TAB_YN for each boundary
!
         IF (JBDY  ==  1) THEN
!   dry / light  boundary
           TAB_YN(1,1) = TAB_RATE(1,1)
           TAB_YN(2,1) = TAB_RATE(2,1) + TAB_RATE(3,1) + TAB_RATE(4,1)
           TAB_YN(1,2) = TAB_RATE(1,2) + TAB_RATE(1,3) + TAB_RATE(1,4)
           TAB_YN(2,2) = TAB_RATE(2,2) + TAB_RATE(3,2) + TAB_RATE(4,2)  &
     &                 + TAB_RATE(2,3) + TAB_RATE(3,3) + TAB_RATE(4,3)  &
     &                 + TAB_RATE(2,4) + TAB_RATE(3,4) + TAB_RATE(4,4)
         ELSEIF (JBDY  ==  2) THEN
!   light / moderate  boundary
           TAB_YN(1,1) = TAB_RATE(1,1) + TAB_RATE(2,1)                  &
     &                 + TAB_RATE(1,2) + TAB_RATE(2,2)
           TAB_YN(2,1) = TAB_RATE(3,1) + TAB_RATE(3,2)                  &
     &                 + TAB_RATE(4,1) + TAB_RATE(4,2)
           TAB_YN(1,2) = TAB_RATE(1,3) + TAB_RATE(2,3)                  &
     &                 + TAB_RATE(1,4) + TAB_RATE(2,4)
           TAB_YN(2,2) = TAB_RATE(3,3) + TAB_RATE(3,4)                  &
     &                 + TAB_RATE(4,3) + TAB_RATE(4,4)
         ELSE
!   moderate / heavy  boundary
           TAB_YN(2,2) = TAB_RATE(4,4)
           TAB_YN(2,1) = TAB_RATE(4,1) + TAB_RATE(4,2) + TAB_RATE(4,3)
           TAB_YN(1,2) = TAB_RATE(1,4) + TAB_RATE(2,4) + TAB_RATE(3,4)
           TAB_YN(1,1) = TAB_RATE(1,1) + TAB_RATE(1,2) + TAB_RATE(1,3)  &
     &                 + TAB_RATE(2,1) + TAB_RATE(2,2) + TAB_RATE(2,3)  &
     &                 + TAB_RATE(3,1) + TAB_RATE(3,2) + TAB_RATE(3,3)
         ENDIF
!
!L*** 4.1  Now calculate statistics for each set of TAB_YN
!

!  Hit Rate : percentage of observed precipitation correctly forecast
         IF ( (TAB_YN(2,1)+TAB_YN(2,2))  >   0) THEN
                HR(JBDY) = TAB_YN(2,2) * 100.0                          &
     &                  / ( TAB_YN(2,1) + TAB_YN(2,2) )
         ELSE
            HR(JBDY) = 999999.9
         ENDIF

!  False Alarm Rate : percentage of forecast precipitation observed
!                     to be dry
         IF ( (TAB_YN(1,2)+TAB_YN(2,2))  >   0) THEN
                FAR(JBDY)= TAB_YN(1,2) * 100.0                          &
     &                  / ( TAB_YN(1,2) + TAB_YN(2,2) )
         ELSE
            FAR(JBDY)= 9999999.9
         ENDIF

!  Hanssen & Kuiper : Hit Rate - percentage of dry observations
!                              incorrectly forecast
         IF ( (TAB_YN(1,1)+TAB_YN(1,2))  >   0) THEN
              HKS(JBDY) = HR(JBDY) - TAB_YN(1,2) * 100.0                &
     &                           / ( TAB_YN(1,1) + TAB_YN(1,2) )
         ELSE
            HKS(JBDY) = 999999.9
         ENDIF

!  Threat Score : Correctly forecast precipitation
!                    / ( total observed precipitation
!                      + total forecast precipitation
!                      - correct forecast precipitation )
!  threat score is also known as CSI - critical success index
         IF ( (TAB_YN(1,2)+TAB_YN(2,1)+TAB_YN(2,2))  >   0.) THEN
                TS(JBDY) = TAB_YN(2,2) * 100.0                          &
     &                  / ( TAB_YN(1,2) + TAB_YN(2,1) + TAB_YN(2,2) )
         ELSE
            TS(JBDY) = 999999.9
         ENDIF

!   Tetrachoric correlation coefficient
        IF ( (JBDY  ==  1) .AND. ( (TAB_YN(2,1)+TAB_YN(2,2))  >   0.0)  &
     &       .AND.  ( (TAB_YN(1,2)+TAB_YN(2,2))  >   0.0 )              &
     &       .AND.  ( (TAB_YN(1,1)+TAB_YN(1,2))  >   0.0 )              &
     &       .AND.  ( (TAB_YN(1,1)+TAB_YN(2,1))  >   0.0 ) ) THEN
          CC = ( TAB_YN(1,1)*TAB_YN(2,2) - TAB_YN(1,2)*TAB_YN(2,1) ) /  &
     &       SQRT( (TAB_YN(2,1)+TAB_YN(2,2))*(TAB_YN(1,2)+TAB_YN(2,2))  &
     &       * (TAB_YN(1,1)+TAB_YN(1,2)) * (TAB_YN(1,1)+TAB_YN(2,1)) )
        ENDIF
      END DO ! jbdy

!
!L*** 5.0 Calculate overall scores
!
!  percentage of all obs which give precipitation
      OBS_PCNT = 100.0 - ( TAB_RATE(1,1) + TAB_RATE(1,2) +              &
     &                        TAB_RATE(1,3) + TAB_RATE(1,4) )

!  percentage of all forecasts which give precipitation
      FC_PCNT = 100.0 - ( TAB_RATE(1,1) + TAB_RATE(2,1) +               &
     &                        TAB_RATE(3,1) + TAB_RATE(4,1) )

!  percentage of forecasts in correct category
      CRCT_PCNT = TAB_RATE(1,1) + TAB_RATE(2,2) +                       &
     &                   TAB_RATE(3,3) + TAB_RATE(4,4)

!  Bias Ratio : forecast precipitation / observed precipitation
!
!  BIAS_RATE : by mean rate
      IF (AVG_OBRATE  >   0.0) THEN
         BIAS_RATE = AVG_FCRATE / AVG_OBRATE
      ELSE
         BIAS_RATE = 99999.9
      ENDIF
!  BIAS_AREA : by area
      IF (OBS_PCNT  >   0.0) THEN
         BIAS_AREA = FC_PCNT / OBS_PCNT
      ELSE
         BIAS_AREA = 99999.9
      ENDIF

!L*** 6.0 Print results

      PRINT '('' '')'
      PRINT '(''  * table entries are % of obs'')'
      PRINT '('' '')'
      PRINT '(''  Precipitation Rate Table'')'
      PRINT '(''  ------------------------'')'
      PRINT '('' '')'
      PRINT '(21X,''OBS'')'
      PRINT '('' '')'
      PRINT '(16X,'' Dry    Light  Moderate  Heavy   TOTAL'')'
      PRINT '(13X,''------------------------------------------'')'
      PRINT '(''  F.C.    Dry| '',F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)', &
          TAB_RATE(1,1),TAB_RATE(2,1),TAB_RATE(3,1),                      &
          TAB_RATE(4,1),TAB_RATE(5,1)
      PRINT '(''        Light| '',F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)', &
          TAB_RATE(1,2),TAB_RATE(2,2),TAB_RATE(3,2),                       &
          TAB_RATE(4,2),TAB_RATE(5,2)
      PRINT '(''     Moderate| '',F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)', &
          TAB_RATE(1,3),TAB_RATE(2,3),TAB_RATE(3,3),                       &
          TAB_RATE(4,3),TAB_RATE(5,3)
      PRINT '(''        Heavy| '',F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)',&
          TAB_RATE(1,4),TAB_RATE(2,4),TAB_RATE(3,4),                       &
          TAB_RATE(4,4),TAB_RATE(5,4)
      PRINT '(''       TOTAL | '',F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)', &
          TAB_RATE(1,5),TAB_RATE(2,5),TAB_RATE(3,5),                       &
          TAB_RATE(4,5),100.0
      PRINT '('' '')'
      PRINT '(''  STATISTICS'')'
      PRINT '(''  ----------'')'
      PRINT '('' '')'
      PRINT '(''    threshold used   :      d/l     l/m    m/h'')'
      PRINT '(''  Hit Rate           HR= '',3F7.1)', HR
      PRINT '(''  False Alarm Rate  FAR= '',3F7.1)', FAR
      PRINT '(''  Hanssen & Kuiper  HKS= '',3F7.1)', HKS
      PRINT '(''  Threat Score       TS= '',3F7.1)', TS
      PRINT '(''    (above scores are %)'')'
      PRINT '('' '')'
      PRINT '(''  % of obs with precip   '',F5.1)', OBS_PCNT
      PRINT '(''  % fc with precip       '',F5.1)', FC_PCNT
      PRINT '(''  Mean Forecast Rate     '',F6.3,''  mm/hour'')',       &
     &                    AVG_FCRATE
      PRINT '(''  Mean Observation Rate  '',F6.3,''  mm/hour'')',       &
     &                    AVG_OBRATE
      PRINT '(''  Bias Ratio (rate)      = '',F5.3)', BIAS_RATE
      PRINT '(''             (area)      = '',F5.3)', BIAS_AREA
      PRINT '(''  % of fc in correct cat = '',F5.1)', CRCT_PCNT
      PRINT '(''  RMS rate error         = '',F7.3,''  mm/hr'')',       &
     &                     RMS*3600
      PRINT '(''  RMS factor (intersection)  RMSF= '',F7.3)', RMSF
      PRINT '(''    (calculated from threshold '',F6.3,'' mm/hr)'')',   &
     &                     THRESH_RMSF*3600
      PRINT '(''     No. of obs passing threshold = '',I6)', COUNTERF
      PRINT '(''  RMS factor (union)         RMSF2= '',F7.3)', RMSF2
      PRINT '(''     No. of obs for RMSF2    = '',I6)', COUNTERF2
      PRINT '(''  Correlation Coeff  CC= '',F7.3)', CC
      PRINT '(''    (based on rain/no rain)'')'
      PRINT '('' '')'
      PRINT '('' '')'
      PRINT '(''  Number of Obs used     '',I6)', OBS_USED
      IF (L_LATLON_PRVER) THEN
       PRINT '(''  Verification in area :'')'
       PRINT '(''    lat '',F6.3,'' deg'',''  to  '',F6.3,'' deg'')',   &
     &                    SOUTHLAT,NORTHLAT
       PRINT '(''    lon '',F6.3,'' deg'',''  to  '',F6.3,'' deg'')',   &
     &                    WESTLON,EASTLON
      ENDIF
      PRINT '('' '')'
      PRINT '(''  Threshold values used :'')'
      PRINT '(''     Dry / Light         '',F6.3,''  mm/hour'')',       &
     &                    THRESH_DL*3600
      PRINT '(''     Light / Moderate    '',F6.3,''  mm/hour'')',       &
     &                    THRESH_LM*3600
      PRINT '(''     Moderate / Heavy    '',F6.3,''  mm/hour'')',       &
     &                    THRESH_MH*3600
9999   CONTINUE

      endif ! mype == 0


! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('DIAGOPR ',4)
      IF (lhook) CALL dr_hook('DIAGOPR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DIAGOPR
END MODULE diagopr_mod
