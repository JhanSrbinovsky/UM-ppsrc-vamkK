! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  4 Subroutines in deck : RDOBS, RDOBS2, RDOBS3 and DAYS -----
!LL
!LL  Purpose : Read from ACOBS Files,reformat and place OBS header
!LL            details in COMOBS. The bulk of the required OBS data
!LL            is put into dynamic work array OBS for transmission via
!LL            argument list to GETOBS. OBS is written out to a cache
!LL            file for subsequent reading at later timesteps.
!LL            Thus reread of ACOBS files only required intermittently
!LL            (The routine DAYS does a dd/mm/yy to dayno)
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: AC Assimilation
MODULE rdobs2_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE RDOBS2(NFILES,TIMESTEP_NO,OBS,OBS_FLAG,TNDV,           &
                        TNOBS,TNDVMAX,NOBSMAX,P_LEVELS,Q_LEVELS,TIMEREL,&
                        lambda_p,phi_p,L_regular,P_ROWS,ROW_LENGTH,     &
                     ICODE,CMESSAGE)
!L----------------------------------------------------------------------
!L   INTENT IN:
!L      NFILES   : NO OF AC OBSERVATION FILES TO BE READ
!L      TIMESTEP_NO : TIMESTEP NUMBER
!L      TIMEREL     : relative time for this timestep
!L      TNDVMAX  : MAX SIZE OF OBS ARRAY
!L      NOBSMAX  : MAX NO OF OBS (FOR DIMENSIONING)
!L      P_LEVELS : NUMBER OF MODEL LEVELS
!L      Q_LEVELS : NUMBER OF WET MODEL LEVELS
!L   INTENT OUT:
!L      TNDV     : ACTUAL SIZE OF OBS ARRAY
!L      TNOBS    : ACTUAL NO OF OBS
!L      OBS      : OBS array
!L      OBS_FLAG : do not use flags
!L      ICODE/CMESSAGE: for error processing
!L
!L ONLY APPLICABLE TO LAMS NOW
!L MOPS DATA ONLY - NO WIND OBS
!L----------------------------------------------------------------------
!L
!
      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global

      USE conversions_mod, ONLY: pi_over_180,pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE IO
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx, ndatavmx, obs_info, used_files
      USE days_mod, ONLY: days
      USE rdobs3_mod, ONLY: rdobs3
      USE setdac_mod, ONLY: setdac
      USE lltoeq_mod, ONLY: lltoeq
      IMPLICIT NONE

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
!--------------------------------------------------------------------
!LCOMDECK COMMG
!L-------------
      REAL          DLAT,DLONG,XLATN,XLONGW
      REAL          ELFPLAT,ELFPLON
      COMMON/COMMG/ DLAT,DLONG,XLATN,XLONGW,ELFPLAT,ELFPLON
!--------------------------------------------------------------------

      integer common_length

!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
      INTEGER TNDV,TNOBS,MAX_NDV
      INTEGER P_LEVELS,Q_LEVELS
      INTEGER TNDVMAX,NOBSMAX
      INTEGER NFILES
      INTEGER TIMESTEP_NO
      INTEGER ICODE
      CHARACTER(LEN=256) CMESSAGE

      INTEGER P_ROWS,ROW_LENGTH
      Real lambda_p(1-halo_i:row_length+halo_i)
      Real phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)
      LOGICAL L_regular

!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER KTYPE
      INTEGER IO_STAT
      INTEGER JF,JOBT,JOBT2,JDV,J,JOB,JTYP,JLEV,JACT,JMOT
      INTEGER ITOT0,ITOT1,ITOTAL
      INTEGER IRDAY,IFDAY
      INTEGER IPT,IPC,IFILE,I,IWORKSP,IOBT,IJF
      INTEGER IPT_THIS_FILE
      INTEGER MAX_NDATAV
      INTEGER INDVMAX
      REAL                                                              &
       ZZLATMX,ZZLATMN,ZZLONGMX,ZZLONGMN,                               &
       TIMEREL,TIMEADJ,TWSTART,TWEND,TGETOBB,TGETOBA
      COMMON/ZZ/ZZLATMX,ZZLATMN,ZZLONGMX,ZZLONGMN
      INTEGER iproc,istat,imsg
      REAL W_LIMIT,E_LIMIT,N_LIMIT,S_LIMIT,DELTA_PHI,DELTA_LAMBDA
!-----------------------------------------------------------------------
      REAL DATALEVS (P_LEVELS+1,NOBTYPMX)
      INTEGER IREF     (NOBTYPMX,NDATAVMX,NFILES)
      INTEGER INOBS    (NOBTYPMX,NFILES)
      INTEGER IOBSTYP  (NOBTYPMX,NFILES)
      INTEGER INDATAV  (NOBTYPMX,NFILES)
      INTEGER INOBLEV  (NOBTYPMX,NFILES)
      INTEGER IOBLVTYP (NOBTYPMX,NFILES)
      INTEGER INDVHDR  (NFILES)
      INTEGER INOBTYP  (NFILES)
      INTEGER IMAXNLEV (NFILES)
      REAL PLEVELS     (P_LEVELS+1,NOBTYPMX,NFILES)
      INTEGER KOBSTYP  (NOBTYPMX)
      LOGICAL LEMPTY   (NFILES)
!-----------------------------------------------------------------------
      INTEGER OBS_FILE_YY,OBS_FILE_MM,OBS_FILE_DD
      INTEGER OBS_FILE_HH,OBS_FILE_MIN,OBS_FILE_SEC
      INTEGER IP_LAT,IP_LONG,IP_TIME,IP_TYPE,IP_MOT,IP_U,IP_V
      INTEGER IP_MOT1,IP_MOT2,IP_TYPE1,IP_TYPE2
      INTEGER IDUMMY
!-----------------------------------------------------------------------
!*L----------------- COMDECK DUMP_LEN --------------------------------

      INTEGER LEN_FIXHD
      INTEGER LEN_INTHD
      INTEGER LEN_REALHD
      INTEGER LEN1_LEVDEPC, LEN2_LEVDEPC
      INTEGER LEN1_ROWDEPC, LEN2_ROWDEPC
      INTEGER LEN1_COLDEPC, LEN2_COLDEPC
      INTEGER LEN1_FLDDEPC, LEN2_FLDDEPC
      INTEGER LEN_EXTCNST
      INTEGER LEN_DUMPHIST
      INTEGER LEN_CFI1, LEN_CFI2, LEN_CFI3
      INTEGER LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS
!*----------------------------------------------------------------------
      PARAMETER (LEN_FIXHD=256)
      INTEGER FIXHD(LEN_FIXHD)
      INTEGER LEN_DATA
!-----------------------------------------------------------------------
!     DYNAMIC ALLOCATION
      INTEGER OBS_FLAG(NOBSMAX)
      INTEGER IWORK(NOBSMAX)
      REAL OBS(TNDVMAX)
      REAL WORK(TNDVMAX+2048)
      REAL U_WRK(NOBSMAX),V_WRK(NOBSMAX)
      REAL WRKLAT(NOBSMAX),WRKLON(NOBSMAX)
      REAL COEFF1(NOBSMAX),COEFF2(NOBSMAX)
      INTEGER ISPARE(7)
      INTEGER ENVVAR

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Calculate limits of domain in degrees
      IF(L_regular) THEN

        N_LIMIT=lat_n
        S_LIMIT=lat_s-0.5*dlat
        W_LIMIT=long_w_model
        E_LIMIT=long_e_model+0.5*dlong

      ELSE

        DELTA_PHI=(PHI_P(1,2)-PHI_P(1,1))/PI_OVER_180
        DELTA_LAMBDA=(LAMBDA_P(ROW_LENGTH)-LAMBDA_P(ROW_LENGTH-1)) &
                 /PI_OVER_180

        N_LIMIT=PHI_P(1,P_ROWS)/PI_OVER_180
        S_LIMIT=PHI_P(1,1)/PI_OVER_180-0.5*DELTA_PHI
        W_LIMIT=LAMBDA_P(1)/PI_OVER_180
        E_LIMIT=LAMBDA_P(ROW_LENGTH)/PI_OVER_180+0.5*DELTA_LAMBDA

      ENDIF
   
!-----------------------------------------------------------------------
!L           SECTION 1: COPY INPUT FILES TO WORK AREA, ETC.
!-----------------------------------------------------------------------

!       Read in the AC Observation files and merge the observations.

      IF (lhook) CALL dr_hook('RDOBS2',zhook_in,zhook_handle)
        IF (DIAG_RDOBS >= 1.AND.mype == 0) THEN
          PRINT *, ' '
          PRINT *, 'READ IN AC OBS FILES - TIMESTEP : ',TIMESTEP_NO
        ENDIF

!       Set up time for next read of observation files
        OBS_INFO % TIMENEXT = TIMEREL + OBS_INFO % TIMEINT

!       Initilaise to zero to allow for empty files.
        DO JF=1,NFILES
          LEMPTY  (JF) = .FALSE.
          INDVHDR (JF) = 0
          INOBTYP (JF) = 0
          IMAXNLEV(JF) = 0
          DO JOBT=1,NOBTYPMX
            IOBSTYP  (JOBT,JF) = 0
            INDATAV  (JOBT,JF) = 0
            INOBS    (JOBT,JF) = 0
            INOBLEV  (JOBT,JF) = 0
            IOBLVTYP (JOBT,JF) = 0
            DO JLEV=1,P_LEVELS+1
              PLEVELS(JLEV,JOBT,JF) = 0.0
            ENDDO
          ENDDO
        ENDDO

!       Initialise IREF
        DO JF = 1,NFILES
          DO JDV = 1,NDATAVMX
            DO JOBT = 1,NOBTYPMX
              IREF(JOBT,JDV,JF) = 0
            ENDDO
          ENDDO
        ENDDO

!       IPT IS THE AMOUNT OF WORK SPACE USED IN THE ARRAY WORK SO FAR
        IPT=0

!       JF LOOP: ONE CYCLE FOR EACH INPUT FILE.
        DO JF=1,OBS_INFO % NUM_USED_FILES

          IFILE = OBS_UNITNO

!       Get unit number of observation file

          ENVVAR = 1
          CALL FILE_OPEN(IFILE,USED_FILES(JF),OBS_INFO % FILENAME_LEN(JF),&
                         0,ENVVAR,ICODE)


!       READ FILE CONTENTS TO BUFFER.
!       IF FILE IS EMPTY - PROCEED TO READ NEXT FILE

          OBS_INFO % NOBTYP = 0
          TNDV   = 0
          TNOBS  = 0
          OBS_INFO % NDVHDR = 0
          DO JOBT=1,NOBTYPMX
            OBS_INFO % OBSTYP  (JOBT) = 0
            OBS_INFO % NDATAV  (JOBT) = 0
            OBS_INFO % NOBLEV  (JOBT) = 0
            OBS_INFO % OBLEVTYP(JOBT) = 0
            OBS_INFO % NOBS    (JOBT) = 0
          ENDDO

!-----------------------------------------------------------------------
          IF (OBS_FORMAT == 2 .OR. OBS_FORMAT == 3) THEN

!         Go to start of obs file
! read in fixed length header

! DEPENDS ON: read_flh
            CALL READ_FLH (IFILE,FIXHD,LEN_FIXHD,ICODE,CMESSAGE)
            IF (ICODE >  0) GOTO 9999

!         Get dimensions of data set components in Fixed length header
! DEPENDS ON: get_dim
            CALL GET_DIM (FIXHD,                                          &
! DUMP_AR2 array dimensions
     &  LEN_FIXHD, LEN_INTHD, LEN_REALHD,                               &
     &  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,         &
     &  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC,LEN2_FLDDEPC,          &
     &  LEN_EXTCNST,  LEN_DUMPHIST,                                     &
     &  LEN_CFI1, LEN_CFI2, LEN_CFI3,                                   &
     &  LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS,                               &
! DUMP_AR2 end
                        LEN_DATA)

            OBS_FILE_YY  = FIXHD(21)   ! )
            OBS_FILE_MM  = FIXHD(22)   ! ) Date/time
            OBS_FILE_DD  = FIXHD(23)   ! ) of
            OBS_FILE_HH  = FIXHD(24)   ! ) Observation
            OBS_FILE_MIN = FIXHD(25)   ! ) File
            OBS_FILE_SEC = FIXHD(26)   ! )

            MAX_NDV      = FIXHD(162)  !  Max total no of data values

!         If MAX_NDV = 0, set to 1. (Obs file with no observations)
!         Prevents dynamic allocation with 0 in READACOBS
            IF (MAX_NDV == 0) MAX_NDV = 1

!
!  check dimension of WORK
!
            common_length = 2048+IPT+ OBS_INFO % PER_FILE_TNDVMAX(JF)
! 2048 allows space for headers when no or very few obs; was 3000.
            IF (common_length  >   TNDVMAX+2048) THEN
              WRITE(6,*)'RDOBS2 - common_length,IPT,PER_FILE_TNDVMAX,'  &
                  ,'TNDVMAX+2048 ',common_length,IPT,                   &
                  OBS_INFO % PER_FILE_TNDVMAX,  TNDVMAX+2048
!
! important failure messages to PEn,OUT and OUT2 output streams
              WRITE(0,*)' RDOBS2 : Insufficient space in WORK array'
              WRITE(0,*) ' Recode or rerun with fewer obs'
              WRITE(0,*)'TNDVMAX',TNDVMAX,' should be >= ',common_length-2048
              WRITE(6,*)' RDOBS2 : Insufficient space in WORK array'
              WRITE(6,*) ' Recode or rerun with fewer obs'
              WRITE(6,*)'TNDVMAX',TNDVMAX,' should be >= ',common_length-2048
              WRITE(7,*)' RDOBS2 : Insufficient space in WORK array'
              WRITE(7,*) ' Recode or rerun with fewer obs'
              ICODE = 1
              CMESSAGE =' RDOBS2 : Insufficient space in WORK array'
              GO TO 9999
            ELSE
              IF (DIAG_RDOBS >= 1) THEN
! diagnostic message directed to every PEn output stream
                IPC=(common_length*100)/(TNDVMAX+2048)
                WRITE(6,*)'RDOBS2:% of WORK used so far for reading obs=',IPC
              ENDIF
            END IF
          CALL RDOBS3(IFILE,NOBTYPMX,OBS_INFO % NOBTYP,OBS_INFO % OBSTYP,  &
                      OBS_INFO % NOBS, OBS_INFO % NDATAV,                  &
                      OBS_INFO % NOBLEV,OBS_INFO % OBLEVTYP,DATALEVS,      &
                      WORK(IPT+1),                                         &
! DUMP_AR2 array dimensions
     &  LEN_FIXHD, LEN_INTHD, LEN_REALHD,                               &
     &  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,         &
     &  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC,LEN2_FLDDEPC,          &
     &  LEN_EXTCNST,  LEN_DUMPHIST,                                     &
     &  LEN_CFI1, LEN_CFI2, LEN_CFI3,                                   &
     &  LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS,                               &
! DUMP_AR2 end
                      LEN_DATA,MAX_NDV,TNOBS,OBS_INFO % MISSD,P_LEVELS,    &
                      ICODE,CMESSAGE                                       &
                               ,IPT                                        &
                                   )

          IF (ICODE >  0) GOTO 9999
          
          OBS_INFO % NDVHDR = 5                 ! No of data values in header
          OBS_INFO % MAXNLEV1 = LEN1_LEVDEPC-2  ! Max no of levels + 1
          
        ELSE
          ICODE=1
          CMESSAGE='RDOBS2: ILLEGAL OBS_FORMAT'
          GOTO 9999
        ENDIF

!       Convert any old obs type numbers to new obs type numbers.
!       For version 2.6 onwards.
        IF (OBS_INFO % NOBTYP >  0) THEN
          DO JOBT=1,OBS_INFO % NOBTYP
            IF (OBS_INFO % OBSTYP(JOBT) == 501) THEN
              OBS_INFO % OBSTYP(JOBT) = 302
             if(mype == 0)PRINT *, 'Type 501 in Obs file changed to 302'
            ENDIF
            IF (OBS_INFO % OBSTYP(JOBT) == 502) THEN
              OBS_INFO % OBSTYP(JOBT) = 305
             if(mype == 0)PRINT *, 'Type 502 in Obs file changed to 305'
            ENDIF
          ENDDO
        ENDIF

!       PRINT CONTENTS OF FILE HEADER.
        IF (DIAG_RDOBS >= 1.AND.mype == 0) THEN
          PRINT *, ' '
          PRINT '(A,I8)', ' AC OBS FILE - UNIT NO :',IFILE
          PRINT '(A,I8)', ' NO OF OBS TYPES       :',OBS_INFO % NOBTYP
          PRINT '(A,I8)', ' TOTAL NO OF OBS       :',TNOBS
          PRINT *, ' '
          PRINT '(A,I3.2,A,I2.2,A,I4)',                                 &
            ' DATE :',OBS_FILE_DD,'/',OBS_FILE_MM,'/',OBS_FILE_YY
          PRINT '(A,I3.2,I2.2,A)',                                      &
            ' TIME :',OBS_FILE_HH,OBS_FILE_MIN,'Z'
          PRINT *, ' '

          IF (OBS_INFO % NOBTYP >  0) THEN

            PRINT '(A,T30,18I7/(T30,18I7))',                            &
                 ' AC OBS TYPES      :',(OBS_INFO % OBSTYP(I),          &
                  I=1,OBS_INFO % NOBTYP)
            PRINT '(A,T30,18I7/(T30,18I7))',                            &
                  ' NO OF LEVELS      :',(OBS_INFO % NOBLEV(I),         &
                  I=1,OBS_INFO % NOBTYP)
            PRINT '(A,T30,18I7/(T30,18I7))',                            &
                  ' NO OF DATA VALUES :',(OBS_INFO % NDATAV(I),         &
                  I=1,OBS_INFO % NOBTYP)
            PRINT '(A,T30,18I7/(T30,18I7))',                            &
                  ' OBS LEVEL TYPE    :',(OBS_INFO % OBLEVTYP(I),       &
                  I=1,OBS_INFO % NOBTYP)
            PRINT '(A,T30,18I7/(T30,18I7))',                            &
                 ' NO OF OBS         :',                                &
                   (OBS_INFO % NOBS(I),I=1,OBS_INFO % NOBTYP)

          ENDIF

        ENDIF

        IF (TIMESTEP_NO == 1 .AND. TNOBS == 0) THEN
          IF (OBS_INFO % OBSTYP(1) == 406) THEN
            WRITE (7,'(A,I3.2,I2.2,A,I3.2,A,I2.2,A,I4,A)')              &
              ' AC Observation File (MOPS) - ',                         &
              OBS_FILE_HH,OBS_FILE_MIN,'Z',                             &
              OBS_FILE_DD,'/',OBS_FILE_MM,'/',OBS_FILE_YY,              &
              ' - No MOPS ACOBS data'
          ELSE
          WRITE (7,'(A,I3.2,I2.2,A,I3.2,A,I2.2,A,I4,A)')                &
              ' AC Observation File - ',                                &
              OBS_FILE_HH,OBS_FILE_MIN,'Z',                             &
              OBS_FILE_DD,'/',OBS_FILE_MM,'/',OBS_FILE_YY,              &
              ' - No standard ACOBS data'
          ENDIF
        ENDIF

!       Check that no of observation types in file (NOBTYP) does not
!       exceed maximum allowed. (NOBTYPMX)
        IF (OBS_INFO % NOBTYP >  NOBTYPMX) THEN
          ICODE = 1
          CMESSAGE = ' RDOBS2 : TOO MANY OBSERVATION TYPES IN OBS FILE'
          if(mype == 0)then
          PRINT *,' RDOBS2 : TOO MANY OBSERVATION TYPES IN OBS FILE'
          PRINT '(A,I5)',' NO OF OBS TYPES = ',OBS_INFO % NOBTYP
          PRINT '(A,I5)',' MAXIMUM ALLOWED = ',NOBTYPMX
          endif
          GO TO 9999
        ENDIF

!       Check that no of data values for each obs type (NDATAV) does
!       not exceed maximum allowed. (NDATAVMX)
        IF (OBS_INFO % NOBTYP >  0) THEN
          DO JOBT=1,OBS_INFO % NOBTYP
          IF (OBS_INFO % NDATAV(JOBT) >  NDATAVMX) THEN
            ICODE = 1
            CMESSAGE = ' RDOBS2 : TOO MANY DATA VALUES IN OBS FILE'
            if(mype == 0)then
            PRINT *, ' RDOBS2 : Too many Data values in Obs file'
            PRINT *, ' Observation Types =  ',OBS_INFO % OBSTYP(JOBT)
            PRINT *, ' No of Data Values =  ',OBS_INFO % NDATAV(JOBT)
            PRINT *, ' Maximum allowed   =  ',NDATAVMX
            endif
            GO TO 9999
          ENDIF
          ENDDO
        ENDIF

!       Store information for this obs file
        INOBTYP(JF) = OBS_INFO % NOBTYP
        INDVHDR(JF) = OBS_INFO % NDVHDR
        IMAXNLEV(JF) = OBS_INFO % MAXNLEV1
        DO JOBT=1,OBS_INFO % NOBTYP
          IOBSTYP(JOBT,JF)  = OBS_INFO % OBSTYP(JOBT)
          INDATAV(JOBT,JF)  = OBS_INFO % NDATAV(JOBT)
          INOBLEV(JOBT,JF)  = OBS_INFO % NOBLEV(JOBT)
          IOBLVTYP(JOBT,JF) = OBS_INFO % OBLEVTYP(JOBT)
          DO JLEV=1, OBS_INFO % MAXNLEV1
            PLEVELS(JLEV,JOBT,JF) = DATALEVS(JLEV,JOBT)
          ENDDO
        ENDDO

!       EVALUATE POINTERS TO EACH SUB-VECTOR.
!       IPT=0 BEFORE READING FIRST FILE

        IPT_THIS_FILE = IPT

          DO JOBT=1,OBS_INFO % NOBTYP
            DO JDV=1,OBS_INFO % NDATAV(JOBT)
              IREF(JOBT,JDV,JF)=IPT
              IPT = IPT+ OBS_INFO % NOBS(JOBT)
            ENDDO
          ENDDO

!       CONVERT REFERENCE DATE & FILE DATE TO ELAPSED DAYS

        CALL DAYS (OBS_INFO % OBS_REF_DD, OBS_INFO % OBS_REF_MM, &
                   OBS_INFO % OBS_REF_YY, IRDAY)
        CALL DAYS (OBS_FILE_DD,OBS_FILE_MM,OBS_FILE_YY,IFDAY)

!       FIND RELATIVE TIME ADJUSTMENT.

        TIMEADJ = 1440*(IFDAY-IRDAY) +                                  &
         60*(OBS_FILE_HH-OBS_INFO % OBS_REF_HH) + OBS_FILE_MIN - OBS_INFO % OBS_REF_MIN

        IF (OBS_INFO % NOBTYP >  0) THEN

!         LOOP OVER OBSERVATION TYPES

          DO JOBT=1,OBS_INFO % NOBTYP

            TGETOBB = 0.0
            TGETOBA = 0.0
            DO JOBT2=1,NOBTYPMX
              IF (OBS_INFO % OBSTYP(JOBT) == MASTER_AC_TYPES(JOBT2)) THEN
                TGETOBB = DEF_TGETOBB(JOBT2)
                TGETOBA = DEF_TGETOBA(JOBT2)
              ENDIF
            ENDDO
            IF (TGETOBB == 0.0 .OR. TGETOBA == 0.0) THEN
              ICODE = 1
              CMESSAGE = 'RDOBS2 : TGETOBB and/or TGETOBA = 0.0 ?'
              GO TO 9999
            ENDIF

!         Set up time window to get observations required for assim.
          TWSTART  = TIMEREL - TGETOBA
          TWEND    = TIMEREL + OBS_INFO % TIMEINT + TGETOBB

        ZZLATMX = N_LIMIT
        ZZLATMN = S_LIMIT
        if(at_extremity(PNorth))ZZLATMX= OBS_INFO % OBS_LAT_N
        if(at_extremity(PSouth))ZZLATMN= OBS_INFO % OBS_LAT_S
         ZZLONGMX= E_LIMIT
         ZZLONGMN= W_LIMIT
        if(at_extremity(PEast))ZZLONGMX= OBS_INFO % OBS_LONG_E
        if(at_extremity(PWest))ZZLONGMN= OBS_INFO % OBS_LONG_W
        if(ZZLONGMN >= 360.0)ZZLONGMN=ZZLONGMN-360.0
        if(ZZLONGMX >= 360.0)ZZLONGMX=ZZLONGMX-360.0

        IP_LAT  = 1
        IP_LONG = 2
        IP_TIME = 3

        IF (OBS_INFO % NOBS(JOBT) >  0) THEN
!
!  Ensure boundaries are numerically the same.

! process PEs not along eastern boundary
        IF(.NOT.AT_EXTREMITY(PEAST))THEN

          DO IPROC=1,NPROC-1

            IMSG=JOBT*100+IPROC  ! message tag

            IF(mype == IPROC) THEN
              IF(.NOT.AT_EXTREMITY(PWEST)) THEN
! this PE receives if not along western boundary
                CALL GC_RRECV(IMSG,1,IPROC-1,ISTAT,ZZLONGMN,ZZLONGMX)
              ENDIF
            ELSEIF(mype == IPROC-1) THEN
! PEs not along eastern boundary send to adjacent PE
              CALL GC_RSEND(IMSG,1,IPROC,ISTAT,ZZLONGMN,ZZLONGMX)
            ENDIF

          ENDDO

        ELSEIF(.NOT.AT_EXTREMITY(PWEST)) THEN
! PEs along eastern boundary can receive only

          DO IPROC=1,NPROC-1
            IF(mype == IPROC) THEN        ! this PE receives
              IMSG=JOBT*100+IPROC
              CALL GC_RRECV(IMSG,1,IPROC-1,ISTAT,ZZLONGMN,ZZLONGMX)
            ENDIF
          ENDDO

        ENDIF

! Now along southern/northern boundaries

! process PEs not along southern boundary
        IF(.NOT.AT_EXTREMITY(PSOUTH))THEN

          DO IPROC=NPROC-1,0,-1

            IMSG=JOBT*1000+IPROC  ! message tag

            IF(mype == IPROC-nproc_x) THEN
              IF(.NOT.AT_EXTREMITY(PNORTH)) THEN
! this PE receives if not along northern boundary
                CALL GC_RRECV(IMSG,1,IPROC,ISTAT,ZZLATMX,ZZLATMN)
              ENDIF
            ELSEIF(mype == IPROC) THEN
! PEs not along southern boundary send to adjacent PE
              CALL GC_RSEND(IMSG,1,IPROC-nproc_x,ISTAT,ZZLATMX,ZZLATMN)
            ENDIF

          ENDDO

        ELSE
! PEs along southern boundary can receive only

          DO IPROC=NPROC-1,0,-1
            IF(mype == IPROC-nproc_x) THEN        ! this PE receives
              IMSG=JOBT*1000+IPROC
              CALL GC_RRECV(IMSG,1,IPROC,ISTAT,ZZLATMX,ZZLATMN)
            ENDIF
          ENDDO

        ENDIF


        IF (model_domain /= mt_global) THEN
!       ROTATE REAL LAT/LON OF OBS TO ELF CO-ORDINATES
!       WRITE TRANSFORMED LAT,LON BACK TO SAME AREA OF WORK ARRAY
!       OUTPUT LONGITUDES FROM LLTOEQ ARE IN RANGE 0-->360 DEGREES.
        CALL LLTOEQ                                                     &
        (WORK(IREF(JOBT,IP_LAT,JF)+1),WORK(IREF(JOBT,IP_LONG,JF)+1),    &
         WORK(IREF(JOBT,IP_LAT,JF)+1),WORK(IREF(JOBT,IP_LONG,JF)+1),    &
         ELFPLAT,ELFPLON,OBS_INFO % NOBS(JOBT) )
        END IF

        DO JOB=1,OBS_INFO % NOBS(JOBT)

!       Reset Observation time so times are relative to start of assm.

        WORK(IREF(JOBT,IP_TIME,JF)+JOB) =                               &
        WORK(IREF(JOBT,IP_TIME,JF)+JOB)+TIMEADJ

!       Test if observation in area and time window.
!   -0.5 on timewindow helps make results same on different machines

         IF ( ZZLONGMN  <   ZZLONGMX ) THEN

          IF ( WORK(IREF(JOBT,IP_LAT, JF)+JOB)  <   ZZLATMX .AND.       &
               WORK(IREF(JOBT,IP_LAT, JF)+JOB)  >=  ZZLATMN .AND.       &
               WORK(IREF(JOBT,IP_LONG,JF)+JOB)  <   ZZLONGMX  .AND.     &
               WORK(IREF(JOBT,IP_LONG,JF)+JOB)  >=  ZZLONGMN  .AND.     &
               WORK(IREF(JOBT,IP_TIME,JF)+JOB)  <   TWEND-0.5 .AND.     &
               WORK(IREF(JOBT,IP_TIME,JF)+JOB)  >   TWSTART+0.5  ) THEN

!           Count no of obs in area and time window and
!           record those observations.

            INOBS(JOBT,JF)        = INOBS(JOBT,JF)+1
            IWORK(INOBS(JOBT,JF)) = JOB

          ENDIF

        ELSE

          IF ( WORK(IREF(JOBT,IP_LAT, JF)+JOB)  <   ZZLATMX  .AND.      &
               WORK(IREF(JOBT,IP_LAT, JF)+JOB)  >=  ZZLATMN  .AND.      &
             ( WORK(IREF(JOBT,IP_LONG,JF)+JOB)  <   ZZLONGMX .OR.       &
               WORK(IREF(JOBT,IP_LONG,JF)+JOB)  >=  ZZLONGMN ) .AND.    &
               WORK(IREF(JOBT,IP_TIME,JF)+JOB)  <   TWEND-0.5 .AND.     &
               WORK(IREF(JOBT,IP_TIME,JF)+JOB)  >   TWSTART+0.5  ) THEN

!           Count no of obs in area and time window and
!           record those observations.

            INOBS(JOBT,JF)        = INOBS(JOBT,JF)+1
            IWORK(INOBS(JOBT,JF)) = JOB

          ENDIF
        ENDIF

        END DO

!       Compress out observations not required.
        DO JDV=1,OBS_INFO % NDATAV(JOBT)
        DO JOB=1,INOBS(JOBT,JF)
        WORK(IREF(JOBT,JDV,JF)+JOB)=WORK(IREF(JOBT,JDV,JF)+IWORK(JOB))
        ENDDO
        ENDDO

        ENDIF

        END DO ! jobt

        ENDIF

!       Print no of observations in time window for file JF.
        DO JOBT=1,OBS_INFO % NOBTYP
        CountA(JOBT)=INOBS(JOBT,JF)
        ENDDO

        If(mype == 0)Then
          Do JOBT=1,OBS_INFO % NOBTYP
            CountC(JOBT)=0
          EndDo
        Endif

        Do iproc=0,nproc-1
          IF(mype == 0) THEN
            CALL GC_IRECV(IPROC,OBS_INFO % NOBTYP,IPROC,ISTAT,COUNTB,COUNTA)
            Do JOB=1,OBS_INFO % NOBTYP
              CountC(JOB)=CountC(JOB)+CountB(JOB)
            EndDo
          ELSEIF(mype == IPROC) THEN
            CALL GC_ISEND(IPROC,OBS_INFO % NOBTYP,0,ISTAT,COUNTB,COUNTA)
          ENDIF
!          IF (DIAG_RDOBS >  1.AND.mype == 0) THEN
            PRINT '(A,I3,T30,18I7/(T30,18I7))',                         &
                  ' NO OF OBS IN T.W: pe=',iproc,                       &
                  (CountB(JOBT),JOBT=1,OBS_INFO % NOBTYP)
!          ENDIF
        EndDo
!        If(mype == 0) THEN
          PRINT *, ' '
          PRINT '(A,T30,18I7/(T30,18I7))',                              &
                ' NO OF OBS IN T.W: total:',                            &
                (CountC(JOBT),JOBT=1,OBS_INFO % NOBTYP)
!        ENDIF

!       Compress out unused work space
!       ==============================

!       Following commented line may be brought back in future

        IPT = IPT_THIS_FILE
        DO JOBT=1,OBS_INFO % NOBTYP
          DO JDV=1,OBS_INFO % NDATAV(JOBT)
            IF (INOBS(JOBT,JF) >  0) THEN
              DO JOB=1,INOBS(JOBT,JF)
              WORK(IPT+JOB) = WORK(IREF(JOBT,JDV,JF)+JOB)
              ENDDO
            ENDIF
            IREF(JOBT,JDV,JF)=IPT
            IPT = IPT+INOBS(JOBT,JF)
          ENDDO
        ENDDO
        IF (DIAG_RDOBS >= 1) THEN
          IPC = (IPT*100)/TNDVMAX
! diagnostic message directed to every PEn output stream
          WRITE(6,*)'RDOBS2:% of OBS array required=',IPC
        ENDIF

        CALL FILE_CLOSE(IFILE,USED_FILES(JF),OBS_INFO % FILENAME_LEN(JF),&
                        ENVVAR,0,ICODE)

        END DO ! JF=1,NUM_USED_FILES

!-----------------------------------------------------------------------

!         DETERMINE LIST OF AC OBS TYPES FOR MERGED FILES
!         NOBTYP IS NOW THE NO OF AC OBS TYPES IN THE MERGED LIST

!         INITIAL LIST IS LIST FOR FIRST FILE
          OBS_INFO % NOBTYP=INOBTYP(1)
          IF (OBS_INFO % NOBTYP >  0) THEN
            DO JOBT=1,OBS_INFO % NOBTYP
              KOBSTYP(JOBT) = IOBSTYP(JOBT,1)
            ENDDO
          ENDIF

        IF (NFILES >  1) THEN

!         GET FULL LIST FROM OTHER FILES
          DO JF=2,NFILES
          IF (.NOT.LEMPTY(JF)) THEN
          DO JTYP=1,INOBTYP(JF)
          DO JOBT=1,OBS_INFO % NOBTYP
            IF (IOBSTYP(JTYP,JF) == KOBSTYP(JOBT)) GO TO 2020
          END DO
          OBS_INFO % NOBTYP = OBS_INFO % NOBTYP+1
          IF (OBS_INFO % NOBTYP <= NOBTYPMX) THEN
            KOBSTYP(OBS_INFO % NOBTYP) = IOBSTYP(JTYP,JF)
          ELSE
            ICODE = 1
            CMESSAGE =                                                  &
            ' RDOBS2 : MAX NO OF AC OBS TYPES REACHED FOR MERGED FILES'
            if(mype == 0)PRINT *,                                       &
            ' RDOBS : MAX NO OF AC OBS TYPES REACHED FOR MERGED FILES'
            GO TO 9999
          ENDIF
          
 2020     CONTINUE
          END DO ! jtyp
          ENDIF
          ENDDO

        ENDIF

!       NOBTYP is now set to the no of obs types to be assimilated
!       which is the same as NACT.
        OBS_INFO % NOBTYP = NACT

        DO JOBT=1,OBS_INFO % NOBTYP
          OBS_INFO % OBSTYP(JOBT)   = -1
          OBS_INFO % NOBS  (JOBT)   =  0
          OBS_INFO % NDATAV(JOBT)   = -1
          OBS_INFO % NOBLEV(JOBT)   = -1
          OBS_INFO % OBLEVTYP(JOBT) = -1
          DO JLEV=1, OBS_INFO % MAXNLEV1
            DATALEVS(JLEV,JOBT) = OBS_INFO % MISSD
          ENDDO
        ENDDO

        IF (OBS_INFO % NOBTYP >  0) THEN

!         OBSTYP is now the list of obs types to be assimilated
!         which is the same as LACT.
          DO JACT=1,NACT
            OBS_INFO % OBSTYP(JACT) = LACT(JACT)
          ENDDO

          DO JOBT=1,OBS_INFO % NOBTYP
            DO JF=1,NFILES
              IF (.NOT.LEMPTY(JF)) THEN
              DO JTYP=1,INOBTYP(JF)
                IF (IOBSTYP(JTYP,JF) == OBS_INFO % OBSTYP(JOBT)) THEN

                  OBS_INFO % NDATAV(JOBT)   = INDATAV(JTYP,JF)
                  OBS_INFO % NOBLEV(JOBT)   = INOBLEV(JTYP,JF)
                  OBS_INFO % OBLEVTYP(JOBT) = IOBLVTYP(JTYP,JF)

                  DO JLEV=1, OBS_INFO % MAXNLEV1
                    DATALEVS(JLEV,JOBT) = PLEVELS(JLEV,JTYP,JF)
                  ENDDO
                  GO TO 2110

                ENDIF
              ENDDO
              ENDIF
            ENDDO
 2110       CONTINUE
          ENDDO

        ENDIF

!       Check that data is consistent for obs types being assimilated.

        IF (NFILES >  1 .AND. OBS_INFO % NOBTYP >  0) THEN

          DO JF=1,NFILES

!         ChecK NDVHDR

          IF (.NOT.LEMPTY(JF) .AND. INDVHDR(JF) /= OBS_INFO % NDVHDR) THEN
            ICODE = 1
            CMESSAGE =' RDOBS2 : MIS-MATCH IN AC OBS FILES - NDVHDR'
            if(mype == 0)then
            PRINT *,' '
            PRINT *,' RDOBS2 : DIFFERENT VALUES OF NDVHDR ?'
            PRINT '(A,5I5)',' NDVHDR =',(INDVHDR(IJF),IJF=1,NFILES)
            endif
            GO TO 9999
          ENDIF

!         Check MAXNLEV1

          IF (.NOT.LEMPTY(JF) .AND. IMAXNLEV(JF) /= OBS_INFO % MAXNLEV1) THEN
            ICODE = 1
            CMESSAGE =' RDOBS2 : MIS-MATCH IN AC OBS FILES - MAXNLEV1'
            if(mype == 0)then
            PRINT *,' '
            PRINT *,' RDOBS2 : DIFFERENT VALUES OF MAXNLEV1'
            PRINT '(A,5I5)',' MAXNLEV1 =',(IMAXNLEV(IJF),IJF=1,NFILES)
            endif
            GO TO 9999
          ENDIF

          ENDDO

          DO JOBT=1,OBS_INFO % NOBTYP
            DO JF=1,NFILES
            IF (INOBTYP(JF) >  0) THEN
              DO JTYP=1,INOBTYP(JF)

                IF (IOBSTYP(JTYP,JF) == OBS_INFO % OBSTYP(JOBT)) THEN

!               Check no of data values (NDATAV)

                IF (INDATAV(JTYP,JF) /= OBS_INFO % NDATAV(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =' RDOBS2 : MIS-MATCH IN OBS FILES - NDATAV'
                  if(mype == 0)then
                  PRINT *, ' RDOBS2 : Different No of Data Values.'
                  PRINT *, '        : See Obs Type ',OBS_INFO % OBSTYP(JOBT)
                  endif
                  GO TO 9999
                ENDIF

!               Check no of observation levels (NOBLEV)

                IF (INOBLEV(JTYP,JF) /= OBS_INFO % NOBLEV(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =' RDOBS2 : MIS-MATCH IN OBS FILES - NOBLEV'
                  if(mype == 0)then
                  PRINT *, ' RDOBS2 : Different No of Obs levels.'
                  PRINT *, '        : See Obs Type ',OBS_INFO % OBSTYP(JOBT)
                  endif
                  GO TO 9999
                ENDIF

!               Check Observation level type (OBLEVTYP)

                IF (IOBLVTYP(JTYP,JF) /= OBS_INFO % OBLEVTYP(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
                  ' RDOBS2 : MIS-MATCH IN OBS FILES - OBLEVTYP'
                  if(mype == 0)then
                  PRINT *,                                              &
                  ' RDOBS2 : Different Observation Level Type'
                  PRINT *, '        : See Obs Type ',OBS_INFO % OBSTYP(JOBT)
                  endif
                  GO TO 9999
                ENDIF

!               Check pressure levels of observations (DATALEVS)

                DO JLEV=1, OBS_INFO % MAXNLEV1
                 IF (PLEVELS(JLEV,JTYP,JF) /= DATALEVS(JLEV,JOBT)) THEN
                   ICODE = 1
                   CMESSAGE =                                           &
                   ' RDOBS2 : MIS-MATCH IN OBS FILES - DATALEVS'
                  if(mype == 0)then
                   PRINT *, ' RDOBS2 : Different levels in Obs files'
                   PRINT *, '        : See Obs Type ',OBS_INFO % OBSTYP(JOBT)
                   PRINT '(A,I5,A,2F10.1)', ' LEVEL',JLEV,              &
                   ' ; Pressures =',PLEVELS(JLEV,JTYP,JF),              &
                   DATALEVS(JLEV,JOBT)
                   endif
                   GO TO 9999
                 ENDIF
                ENDDO

                ENDIF

              ENDDO  !  End of loop over JTYP
            ENDIF
            ENDDO  !  End of loop over JF
          ENDDO  !  End of loop over JOBT

        ENDIF

!-----------------------------------------------------------------------
!L           SECTION 2: MERGE INPUT DATA & SET UP HEADER.
!-----------------------------------------------------------------------

!       Merge observation data and put into output buffer OBS
        TNDV=0
        DO JOBT=1,OBS_INFO % NOBTYP
          DO JDV=1,OBS_INFO % NDATAV(JOBT)
            DO JF=1,NFILES
              DO JTYP=1,INOBTYP(JF)
                IF (IOBSTYP(JTYP,JF) == OBS_INFO % OBSTYP(JOBT) .AND.   &
                    INOBS(JTYP,JF)   >  0) THEN

                  IF (TNDV+INOBS(JTYP,JF) <= TNDVMAX) THEN

                    DO JOB=1,INOBS(JTYP,JF)
                      OBS(TNDV+JOB) = WORK(IREF(JTYP,JDV,JF)+JOB)
                    ENDDO
                    TNDV = TNDV+INOBS(JTYP,JF)

                  ELSE

                    ICODE = 1
                    CMESSAGE =                                          &
                    ' RDOBS2 : Insufficient space in OBS array'
! important failure messages to PEn,OUT and OUT2 output streams
                    WRITE (0,*)                                         &
                    ' RDOBS2 : Insufficient space in OBS array'
                    WRITE (0,*) ' Recode or rerun with fewer obs'
                    WRITE (6,*)                                         &
                    ' RDOBS2 : Insufficient space in OBS array'
                    WRITE (6,*) ' Recode or rerun with fewer obs'
                    WRITE (7,*)                                         &
                    ' RDOBS2 : Insufficient space in OBS array'
                    WRITE (7,*) ' Recode or rerun with fewer obs'
                    GO TO 9999

                  ENDIF
                ENDIF
              ENDDO   !  Loop over JTYP
            ENDDO   !  Loop over JF
          ENDDO   !  Loop over JDV
        ENDDO   !  Loop over JOBT
      WRITE(6,*)'temp check OBS ',TNDV,TNDVMAX,(TNDV*TNDVMAX)*100.0,'%'

!       Count total no of obs for each obs type
        DO JOBT=1,OBS_INFO % NOBTYP
          OBS_INFO % NOBS(JOBT)=0
          DO JF=1,NFILES
            DO JTYP=1,INOBTYP(JF)
            IF (IOBSTYP(JTYP,JF) == OBS_INFO % OBSTYP(JOBT)) THEN
              OBS_INFO % NOBS(JOBT) = OBS_INFO % NOBS(JOBT)+INOBS(JTYP,JF)
            ENDIF
            ENDDO
          ENDDO
        ENDDO

!       Get total no of obs (TNOBS)
        TNOBS = 0
        DO JOBT=1,OBS_INFO % NOBTYP
          TNOBS = TNOBS + OBS_INFO % NOBS(JOBT)
        ENDDO

!       Set up pointers to start of data for each obs type
        OBS_INFO % MDISPOBT(1)=0
        DO JOBT=2,OBS_INFO % NOBTYP
          OBS_INFO % MDISPOBT(JOBT) = OBS_INFO % MDISPOBT(JOBT-1) +       &
                                      OBS_INFO % NOBS(JOBT-1) *           &
                                      OBS_INFO % NDATAV(JOBT-1)
        ENDDO

!       Offset to first obs for each type
        OBS_INFO % OBS_NO_ST(1)=0
        DO JOBT=2,OBS_INFO % NOBTYP
          OBS_INFO % OBS_NO_ST(JOBT) = OBS_INFO % OBS_NO_ST(JOBT-1) +   &
                                       OBS_INFO % NOBS(JOBT-1)
        ENDDO


!       Check no of observations against maximum allowed (NOBSMAX)
        IF (TNOBS >  NOBSMAX) THEN
          ICODE = 1
          CMESSAGE = ' RDOBS2 : Insufficient space in OBS_FLAG array'
! check in every pe?
          WRITE(0,*)' RDOBS2 : Insufficient space in OBS_FLAG array'
          WRITE(0,*)' Recode or rerun with fewer obs'
          WRITE(0,'(A,I8)')' NO OF OBS   = ',TNOBS
          WRITE(6,*)' RDOBS2 : Insufficient space in OBS_FLAG array'
          WRITE(6,*)' Recode or rerun with fewer obs'
          WRITE(6,'(A,I8)')' NO OF OBS   = ',TNOBS
          WRITE(7,*)' RDOBS2 : Insufficient space in OBS_FLAG array'
          WRITE(7,*)' Recode or rerun with fewer obs'
          WRITE(7,'(A,I8)')' NO OF OBS   = ',TNOBS
          WRITE(0,'(A,I8)')' MAX ALLOWED = ',NOBSMAX
          WRITE(6,'(A,I8)')' MAX ALLOWED = ',NOBSMAX
          WRITE(7,'(A,I8)')' MAX ALLOWED = ',NOBSMAX
          GO TO 9999
         ELSE
          IF (DIAG_RDOBS >= 1) THEN
           IPC = (TNOBS*100)/NOBSMAX
! diagnostic message directed to every PEn output stream
          WRITE(6,*)'RDOBS2:% of OBS_FLAG array required=',IPC
         ENDIF
        ENDIF

!       Print summary of output file.
        IF (DIAG_RDOBS >  0) THEN
        DO JOBT=1,OBS_INFO % NOBTYP
        CountA(JOBT)=OBS_INFO % NOBS(JOBT)
        ENDDO

        If(mype == 0)Then
          PRINT *, ' '
          PRINT '(A,I3.2,A,I2.2,A,I4)',                                 &
          ' REF DATE :',OBS_INFO % OBS_REF_DD,'/',                      &
                        OBS_INFO % OBS_REF_MM,'/',                      &
                        OBS_INFO % OBS_REF_YY
          PRINT '(A,I3.2,I2.2,A)',                                      &
          ' REF TIME :',OBS_INFO % OBS_REF_HH, OBS_INFO % OBS_REF_MIN,'Z'
          PRINT *, ' '
          PRINT '(A,F8.1,A)', ' REL TIME      :',TIMEREL,'M'
          PRINT '(A,F8.1,A)', ' REL T.W START :',TWSTART,'M'
          PRINT '(A,F8.1,A)', ' REL T.W END   :',TWEND,'M'
          PRINT *, ' '
          PRINT '(A,T30,18I6/(T30,18I6))',                              &
                ' AC OBS TYPES :',(OBS_INFO % OBSTYP(JOBT),             &
                                   JOBT=1,OBS_INFO % NOBTYP)
          PRINT *, ' '
          Do JOBT=1,OBS_INFO % NOBTYP
            CountC(JOBT)=0
          EndDo
        ENDIF

        Do iproc=0,nproc-1
          IF(mype == 0) THEN
            CALL GC_IRECV(IPROC,OBS_INFO % NOBTYP,IPROC,ISTAT,COUNTB,COUNTA)
            Do JOB=1,5
              CountC(JOB)=CountC(JOB)+CountB(JOB)
            EndDo
          ELSEIF(mype == IPROC) THEN
            CALL GC_ISEND(IPROC,OBS_INFO % NOBTYP,0,ISTAT,COUNTB,COUNTA)
          ENDIF

          IF (DIAG_RDOBS >  1) THEN
          PRINT '(A,I3,T30,18I6/(T30,18I6))',                           &
              ' NO OF OBS    : pe=',iproc,(CountB(JOBT),JOBT=1,         &
               OBS_INFO % NOBTYP)
          ENDIF
        EndDo
        if(mype == 0) PRINT '(A,T30,18I6/(T30,18I6))',                  &
               ' NO OF OBS    : total:',(CountC(JOBT),JOBT=1,           &
                OBS_INFO % NOBTYP)
!        Get total no of obs (all pe's)
        If(OBS_INFO % NOBTYP >= 2)Then
          Do JOBT=2,OBS_INFO % NOBTYP
            CountC(1) = CountC(1) + CountC(JOBT)
          Enddo
        Endif  !NOBTYP >= 2
        IF(mype == 0) THEN
          PRINT *, ' '
          PRINT '(A,I8,A)', ' TOTAL NO OF OBS IN T.W :',CountC(1)
!        Any observations to assimilate ?
         If (CountC(1) == 0) Then
          PRINT *, ' '
          PRINT *, 'Timestep ',TIMESTEP_NO
          PRINT *, 'There are no observations to be assimilated.'
          WRITE (7,*) 'Timestep ',TIMESTEP_NO
          WRITE (7,*) 'There are no observations to be assimilated.'
         Endif !TNOBS == 0
        ENDIF

        ENDIF

        IF (OBS_INFO % NOBTYP >  0) THEN

          DO JOBT=1,OBS_INFO % NOBTYP

            IF (OBS_INFO % NOBS(JOBT) >  0) THEN

!             Convert Obs Latitudes (degrees) to Co-latitudes (radians)
              IP_LAT = OBS_INFO % MDISPOBT(JOBT)
              DO JOB=1,OBS_INFO % NOBS(JOBT)
                OBS(IP_LAT+JOB) = (90.0-OBS(IP_LAT+JOB))*PI_OVER_180
              ENDDO

!             Convert Obs Longitudes to Radians in range 0 to 2*PI.
              IP_LONG = OBS_INFO % MDISPOBT(JOBT) + OBS_INFO % NOBS(JOBT)
              DO JOB=1,OBS_INFO % NOBS(JOBT)
                IF (OBS(IP_LONG+JOB)  <   0.0)                          &
                OBS(IP_LONG+JOB) = OBS(IP_LONG+JOB)+360.0
                OBS(IP_LONG+JOB) = OBS(IP_LONG+JOB)*PI_OVER_180
              ENDDO

            ENDIF

          ENDDO

!L        SET UP OBLEVELS AND OBLAYERB (LEVELS IN PASCALS)
!         ------------------------------------------------
          DO JOBT=1,OBS_INFO % NOBTYP

          IF (OBS_INFO % OBLEVTYP(JOBT) == 3 .OR.                       &
              OBS_INFO % OBLEVTYP(JOBT) == 4) THEN

          IF (OBS_INFO % OBLEVTYP(JOBT) == 3) THEN

            DO JLEV=1,OBS_INFO % NOBLEV(JOBT)
              OBS_INFO % OBLEVELS(JLEV,JOBT) = DATALEVS(JLEV,JOBT)
            ENDDO

            DO JLEV=2,OBS_INFO % NOBLEV(JOBT)
              OBS_INFO % OBLAYERB(JLEV,JOBT) =                          &
              SQRT ( OBS_INFO % OBLEVELS(JLEV-1,JOBT) *                 &
                     OBS_INFO % OBLEVELS(JLEV,JOBT) )
            ENDDO

            OBS_INFO % OBLAYERB(1,JOBT) =  OBS_INFO % OBLEVELS(1,JOBT) *&
             OBS_INFO % OBLEVELS(1,JOBT) / OBS_INFO % OBLAYERB(2,JOBT)

            OBS_INFO % OBLAYERB(OBS_INFO % NOBLEV(JOBT)+1,JOBT) =       &
                   OBS_INFO % OBLEVELS(OBS_INFO % NOBLEV(JOBT),JOBT) *  &
                   OBS_INFO % OBLEVELS(OBS_INFO % NOBLEV(JOBT),JOBT) /  &
                   OBS_INFO % OBLAYERB(OBS_INFO % NOBLEV(JOBT),JOBT)

          ENDIF

          IF (OBS_INFO % OBLEVTYP(JOBT) == 4) THEN

            DO JLEV=1,OBS_INFO % NOBLEV(JOBT)+1
              OBS_INFO % OBLAYERB(JLEV,JOBT) = DATALEVS(JLEV,JOBT)
            ENDDO

            DO JLEV=1,OBS_INFO % NOBLEV(JOBT)
              OBS_INFO % OBLEVELS(JLEV,JOBT) =                          &
              SQRT ( OBS_INFO % OBLAYERB(JLEV,JOBT) *                   &
                     OBS_INFO % OBLAYERB(JLEV+1,JOBT) )
            ENDDO

          ENDIF

          if(mype == 0)then
          PRINT *, ' '
          PRINT *, ' Observation Type ',OBS_INFO % OBSTYP(JOBT)
          PRINT *, ' '
          PRINT '(A,/,(1X,5F8.1))', '  Levels (mb) =',                  &
                 (OBS_INFO % OBLEVELS(JLEV,JOBT)*0.01,                  &
                 JLEV=1,OBS_INFO % NOBLEV(JOBT))
          PRINT '(A,/,(1X,5F8.1))', '  Layer boundaries (mb) =',        &
                 (OBS_INFO % OBLAYERB(JLEV,JOBT)*0.01,                  &
                  JLEV=1,OBS_INFO % NOBLEV(JOBT)+1)

          endif
          ENDIF

          ENDDO

      ENDIF
!     ====================================================
!     SET UP ARRAY NERLEV1 WHICH POINTS THE FIRST DATA VALUE
!     CORRESPONDING TO THE FIRST LEVEL OF ERROR RATIO FOR THE OBS TYPE

      IF (OBS_INFO % NOBTYP >  0) THEN
        DO JOBT=1,OBS_INFO % NOBTYP
          OBS_INFO % NERLEV1(JOBT) = OBS_INFO % NDATAV(JOBT) -          &
                                     OBS_INFO % NOBLEV(JOBT)+1
        ENDDO
      ENDIF
!     ====================================================
      IF (NACT >  0) THEN

        if(mype == 0)then
        PRINT   10, (LACT(J),J=1,NACT)
 10     FORMAT('0TYPES TO BE PROCESSED : LACT    =',15I5)
        PRINT   12, (OBS_INFO % NOBLEV(J),J=1,NACT)
 12     FORMAT('          NO OF LEVELS : NOBLEV  =',15I5)
        PRINT   13, (OBS_INFO % NDATAV(J),J=1,NACT)
 13     FORMAT('  NO OF DATA VARIABLES : NDATAV  =',15I5)
        PRINT   14, (OBS_INFO % NERLEV1(J),J=1,NACT)
 14     FORMAT('       FIRST ERROR LEV : NERLEV1 =',15I5)
        endif

      ENDIF
!     ====================================================
!     CALL SETDAC TO SET UP FOR ANY DIAGNOSTICS REQUIRED
      CALL SETDAC (OBS,TNDV)
!     ====================================================

!     RE-USE IWORK IN THIS SECTION

      IF (TNOBS >  0) THEN

        DO JOB=1,TNOBS
          OBS_FLAG(JOB)=0
        ENDDO

!       Convert Analysis Types into INTEGER.

        IP_TYPE = 4
        DO JOBT=1,OBS_INFO % NOBTYP
          IF (OBS_INFO % NOBS(JOBT) >  0) THEN
            IP_MOT = OBS_INFO % MDISPOBT(JOBT) +                        &
                     (IP_TYPE-1) * OBS_INFO % NOBS(JOBT)
            DO JOB=1,OBS_INFO % NOBS(JOBT)
              IWORK(OBS_INFO % OBS_NO_ST(JOBT)+JOB) = NINT( OBS(IP_MOT+JOB) )
            ENDDO
          ENDIF
        ENDDO

        ITOT0=0
        ITOT1=0
        DO JMOT=1,NANALTYP
          IF (IOMITOBS(JMOT) >  0) THEN
            DO JOB=1,TNOBS
              IF (IWORK(JOB) == IOMITOBS(JMOT)) OBS_FLAG(JOB)=1
            ENDDO
            IF (DIAG_RDOBS >= 1) THEN
              ITOT1=0
              DO JOB=1,TNOBS
                IF (OBS_FLAG(JOB) == 1) ITOT1=ITOT1+1
              ENDDO
              ITOTAL=ITOT1-ITOT0
! in every pe?
              WRITE(6,'(A,I5,A,I6)')' ANALYSIS TYPE ',IOMITOBS(JMOT),   &
              ' : NO OF OBSERVATIONS FLAGGED (DO NOT USE)  -',ITOTAL
              ITOT0=ITOT1
            ENDIF
          ENDIF
        ENDDO


!       OPTIONAL PRINT OUT.
        IF (DIAG_RDOBS >= 1 .AND. ITOT1 /= 0) THEN
! in every pe?
          WRITE(6, '(A,I6)')                                            &
          ' TOTAL NO OF OBSERVATIONS FLAGGED (DO NOT USE)  -',ITOT1

        ENDIF

      ENDIF
!-----------------------------------------------------------------------


 9999 CONTINUE
      IF (lhook) CALL dr_hook('RDOBS2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE RDOBS2
END MODULE rdobs2_mod
