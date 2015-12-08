! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINES ACP_NAMEL,ACDIAG_NAMEL  -------------------------------
!LL
!LL  Purpose : Read in AC Namelist (&ACP) and process.
!LL
!LL ACDIAG_NAMEL:
!LL Set defaults for ACDIAG namelist variables.
!LL                  Read in namelist and process.
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical system components covered:
!LL
!LL  Project Task : P3
!LL
!LL  External documentation:
!LL
!LLEND -------------------------------------------------------------
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE acdiag_namel_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE ACDIAG_NAMEL (ICODE,CMESSAGE)

      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global
      USE conversions_mod, ONLY: pi_over_180

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      USE lltoeq_mod, ONLY: lltoeq
      IMPLICIT NONE

! Global parameters:

! AC Comdecks:
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

!*
! Exported variables (INTENT=OUT)
      INTEGER       ICODE              ! Non zero for failure

      CHARACTER(LEN=256) CMESSAGE           ! Reason for failure

      REAL LAT(4),LON(4)
      REAL :: DAGLAT_temp(1)
      REAL :: DAGLON_temp(1) 
      REAL XTEMP
      REAL ZLAT,ZLONG
      INTEGER J

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Namelist acdiag
      NAMELIST /ACDIAG/ LDIAGAC,LLDAC,MODACO,MDACO,                     &
     &          DLATN,DLATS,DLONGW,DLONGE,NDACPRT,                      &
     &          MDIAGTPS,LLBAND,MDIAG,                                  &
     &          DAGLAT,DAGLON,LLDAG0,LLDAG,                             &
     &          LNORMF,LRMS,LTEMP,LVERIF
! End of header.

      IF (lhook) CALL dr_hook('ACDIAG_NAMEL',zhook_in,zhook_handle)

!L 1.1 Set namelist defaults for variables defined by UI

      LDIAGAC=.FALSE.

!L 1.2 Read in Namelist set by UI

! rewind required to ensure first ACDIAG is read
      REWIND 5
      READ (5, ACDIAG, END=120, ERR=120)

! QA Fortran recommends use of IOSTAT but it doesn't work with namelist
      GOTO 121
120   CONTINUE
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL : Error reading 1st ACDIAG namelist'
        if(mype == 0) PRINT *, CMESSAGE
        GOTO 999
121   CONTINUE

!L 1.3 Set namelist defaults for variables defined by user
      LLDAC(1) = .FALSE.
      LLDAC(2) = .TRUE.
      LLDAC(3) = .FALSE.
      LLDAC(4) = .FALSE.
      LLDAC(5) = .FALSE.

      LLDAG0   = .FALSE.
      LLDAG(1) = .FALSE.
      LLDAG(2) = .FALSE.
      LLDAG(3) = .FALSE.
      LLDAG(4) = .FALSE.
      LLDAG(5) = .FALSE.

      LNORMF = .TRUE.
      LRMS   = .TRUE.
      LTEMP  = .TRUE.
      LVERIF = .FALSE.

      DO J = 1, NDACOP
        MDACO(J) = 0

      END DO
      MODACO = 1

      DO J = 1, NOBTYPMX
        MDIAGTPS(J) = 0

      END DO

      LLBAND  = .FALSE.
      MDIAG   = 0
      NDACPRT = 10

      DLATN   = 55.0
      DLATS   = 50.0
      DLONGW  =-10.0
      DLONGE  =- 5.0

      DAGLAT  = 57.0
      DAGLON  = 340.0

      NDGPRT  = 6

!L 1.4 Read in Namelist set by user
      READ (5, ACDIAG, END=140, ERR=140)

! QA Fortran recommends use of IOSTAT but it doesn't work with namelist
! Jump over NAMELIST read error trap.
      GOTO 141

! Namelist read error trap:
140   CONTINUE
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL: Error reading 2nd ACDIAG namelist'
        if(mype == 0) PRINT *, CMESSAGE
        GOTO 999

141   CONTINUE

!L 2. Print out Namelist for this run
      if(mype == 0)Then
       PRINT *, ' '
       PRINT *, ' Parameters in Namelist ACDIAG for this run'
       WRITE(6,ACDIAG)
      endif

!L 3. Process the namelist
! Process ACDIAG namelist parameters and set up variables in
! COMACDG for diagnostic control in the assimilation.

! Check MODACO
      IF (MODACO <  1 .OR. MODACO >  3) THEN
        LLDAC(1) = .FALSE.
        if(mype == 0)PRINT *, ' Invalid value for MODACO - ',MODACO,    &
     &                                           ' ; LLDAC(1) set to F'

      END IF

! If 'obs-model' statistics required for verification purposes
! then reset various parameters in ACP and ACDIAG namelist.
      IF (LVERIF) THEN
        LGEO    = .FALSE.   ! Do not calculate geostrophic increments
        LHYDR   = .FALSE.   ! Do not calculate hydrostatic increments
        LHYDROL = .FALSE.   ! Do not calculate hydrology   increments
      L_LHN = .FALSE.       ! Do not calculate LHN increments
        LNORMF  = .FALSE.   ! Set normalisation factors to 1.0
        LWBAL_SF= .FALSE.   ! Do not calculate P* & theta from windincs
        LWBAL_UA= .FALSE.   !   "        "        "         "     "
        LRMS    = .TRUE.    ! Print RMS values.
        LTEMP   = .TRUE.    ! Use Temp Increments

        DO J = 1, NOBTYPMX
          DEF_TGETOBB(J) = 180.0   ! Time Window of Obs - Before
          DEF_TGETOBA(J) = 181.0   ! Time Window of Obs - After
          DEF_OBTHIN(J) = 1        ! Dont reduce data vols by thinning

        END DO
      END IF

      IF (model_domain /= mt_global) THEN
      LAT(1) = DLATN
      LON(1) = DLONGW
      LAT(2) = LAT(1)
      LON(2) = DLONGE
      LAT(3) = DLATS
      LON(3) = LON(2)
      LAT(4) = LAT(3)
      LON(4) = LON(1)

! Transform corners of diagnostic area to elf co-ordinates
      CALL LLTOEQ (LAT,LON,LAT,LON,ELFPLAT,ELFPLON,4)

! Adjust area to be rectangular in elf lat/lon
      DLATS  = 0.5* ( LAT(3) + LAT(4) )
      DLONGW = 0.5* ( LON(1) + LON(4) )
      DLATN  = 0.5* ( LAT(1) + LAT(2) )
      DLONGE = 0.5* ( LON(2) + LON(3) )

! Make sure dlatn >  dlats in new elf co-ords
      IF (DLATS  >   DLATN) THEN
        XTEMP = DLATS
        DLATS = DLATN
        DLATN = XTEMP

      END IF

! The search for obs permits dlongw >  dlonge by assuming
! that an area crossing the meridian is intended, so no test
! as for latitudes.

! Transform diagnostic point to elf co-ordinates
! Move DAGLAT/LON into temp array as LLTOEQ expects arrays
      DAGLAT_temp(1) = DAGLAT 
      DAGLON_temp(1) = DAGLON 
      CALL LLTOEQ (DAGLAT_temp,DAGLON_temp,DAGLAT_temp,DAGLON_temp,ELFPLAT, &
           ELFPLON,1)
      DAGLAT = DAGLAT_temp(1) 
      DAGLON = DAGLON_temp(1) 
      END IF  ! .NOT. GLOBAL

! Convert Latitudes into Co-Latitudes (0-180) (NP-SP)
      DLATS = 90.0 - DLATS
      DLATN = 90.0 - DLATN

! Convert Co-Latitudes to Radians
      DLATS = DLATS * PI_OVER_180
      DLATN = DLATN * PI_OVER_180

! Check Longitudes within Range (0-360 DEGREES)
      IF (DLONGW <  0.0) DLONGW = DLONGW+360.0
      IF (DLONGE <  0.0) DLONGE = DLONGE+360.0

! Convert Longitudes to Radians
      DLONGW = DLONGW * PI_OVER_180
      DLONGE = DLONGE * PI_OVER_180

999   CONTINUE
      IF (lhook) CALL dr_hook('ACDIAG_NAMEL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ACDIAG_NAMEL
END MODULE acdiag_namel_mod
