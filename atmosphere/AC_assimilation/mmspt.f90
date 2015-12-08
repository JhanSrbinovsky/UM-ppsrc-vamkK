! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE MMSPT   MMSPTW   ---------------------------------------
!LL
!LL   2 Subroutines in deck : MMSPT and MMSPTW
!LL   MMSPTW is same as MMSPT for Limited Area Winds.
!LL
!LL  Purpose : Provide weighted Mean and S.D. plus extremes
!LL            on Model grids. Used for Weights and Increments
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: AC Assimilation
MODULE mmspt_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE MMSPT (PVALS,KLEV,KGRID,PNTLAB,                        &
     &                  ROW_LENGTH,P_ROWS,                              &
     &                  phi_p,L_regular)
!
!L    CALCULATE AREA-WEIGHTED MEAN & MEAN-SQUARE
!L    OF A FIELD ON THE MODEL GRID
!L
!L    KLEV : LEVEL OF FIELD (USED IN PRINT OUT ONLY)
!L
!L    IF KGRID=0 DATA ON P*/THETA MODEL GRID
!L            =1 DATA ON WIND MODEL GRID
!L
!L    16 CHARACTER PNTLAB IS PRINTED TO IDENTIFY MEAN & MEAN SQ.
!Lx
      USE conversions_mod, ONLY: pi_over_180
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE
!
      EXTERNAL TIMER
!
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
! only for LTIMER_AC
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
!-----------------------------------------------------------------------
      INTEGER KLEV,KGRID                !IN LEVEL AND GRID IDENTIFER
      INTEGER ROW_LENGTH,P_ROWS         !IN MODEL DIMENSIONS
      REAL    PVALS(ROW_LENGTH,*)       !IN MODEL FIELD
      CHARACTER(LEN=16) PNTLAB               !IN CHARACTER IDENTIFER OF FIELD

      Logical L_regular

      Real phi_p(1-halo_i:ROW_LENGTH+halo_i, 1-halo_j:P_ROWS+halo_j)

!     LOCAL VARIABLES
      REAL PVALS_G(Max2DFieldSize)
      INTEGER row_length_global
      INTEGER p_rows_global
      REAL WT,SUMWT,SUM,ZROW1,ZLAT,ZDLAT,ZM,ZMS,ZMAX,ZMIN,ZRMS
      INTEGER JPTF,JPTL,JROWF,JROWL,JROW,JPT,NPTS
      INTEGER IMAXPT,IMAXRO,IMINPT,IMINRO

      INTEGER ICODE
      CHARACTER (LEN=80) :: CMESSAGE

! Variables required for variable resolution grids
! Since the total number of rows is unknown here it is 
! made allocatable (although glsize could be used)
      REAL PHI_P_LOCAL(ROW_LENGTH,P_ROWS)
      REAL, ALLOCATABLE :: PHI_P_GLOBAL(:,:)
      INTEGER I,J
      REAL CPU_START,CPU_TEST,GET_CPU_TIME

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('MMSPT',zhook_in,zhook_handle)

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('MMSPT   ',3)
      row_length_global=glsize(1,fld_type_p)
      p_rows_global=glsize(2,fld_type_p)

! KGRID must be 0
      IF(KGRID.NE.0) THEN
        WRITE(6,*) 'INVALID KGRID IN MMSPT ',KGRID
        RETURN
      ENDIF

!     JPTF  = FIRST POINT IN ROW
!     JPTL  = LAST  POINT IN ROW
!     JROWF = FIRST ROW
!     JROWL = LAST  ROW

!     OUTSIDE TWO BOUNDARY POINTS OF ELF GRID NOT USED
      JPTF  = 3
      JPTL  = row_length_global-2
      JROWF = 3
      JROWL = p_rows_global-2

! Gather P_VALS onto a global field PVALS_G
! DEPENDS ON: gather_field
      Call Gather_Field( PVALS, PVALS_G, row_length, p_rows,            &
     &                   row_length_global, p_rows_global,              &
     &                   fld_type_p, halo_type_no_halo,                 &
     &                   0, gc_all_proc_group,                          &
     &                   icode, cmessage)

! If variable grid collect together phi coordinates of entire grid
      IF(.NOT.L_regular) THEN

        ALLOCATE (PHI_P_GLOBAL(row_length_global,p_rows_global))

! Copy information from PHI_P to PHI_P_LOCAL
!  excluding halos
         DO I=1,P_ROWS
           DO J=1,ROW_LENGTH
             PHI_P_LOCAL(J,I)=PHI_P(J,I)
           ENDDO
         ENDDO

! Gather up all PHI_P values into global array
! DEPENDS ON: gather_field
        CALL GATHER_FIELD(PHI_P_LOCAL,   PHI_P_GLOBAL,                  &
     &                    row_length,    p_rows,                        &
     &                    row_length_global,  p_rows_global,            &
     &                    fld_type_p,    halo_type_no_halo,             &
     &                    0,             gc_all_proc_group,             &
     &                    ICODE,         CMESSAGE)


      ENDIF


      if(mype == 0)then

      ZROW1 = XLATN - DLAT*(JROWL-JROWF+1)

      ZLAT   = ZROW1*PI_OVER_180
      ZDLAT  = DLAT *PI_OVER_180

!L SET ACCUMULATORS
      ZMAX=pvals_g(JPTF+(JROWF-1)*row_length_global)
      ZMIN=pvals_g(JPTF+(JROWF-1)*row_length_global)

      ZM    = 0.0
      ZMS   = 0.0
      ZRMS  = 0.0
      SUMWT = 0.0
      WT    = 0.0
      IMAXPT= JPTF
      IMAXRO= JROWF
      IMINPT= JPTF
      IMINRO= JROWF

      DO JROW = JROWF,JROWL
        IF(L_regular) THEN
          WT = COS(ZLAT)
        ELSE
         WT = COS(phi_p_GLOBAL(1,JROW+1)) ! note offset by 1
        ENDIF

        SUMWT = SUMWT+WT

!     CALCULATE MEAN

      SUM = 0.0
      DO JPT=JPTF,JPTL
        SUM = SUM + pvals_g(JPT+(JROW-1)*row_length_global)
      ENDDO
      ZM = ZM + SUM*WT

!     CALCULATE MEAN SQUARE

      SUM = 0.0
      DO JPT=JPTF,JPTL
        SUM = SUM + pvals_g(JPT+(JROW-1)*row_length_global)*            &
     &              pvals_g(JPT+(JROW-1)*row_length_global)
      ENDDO
      ZMS = ZMS + SUM*WT

!     CALCULATE MAX

      DO JPT=JPTF,JPTL
        IF(pvals_g(JPT+(JROW-1)*row_length_global) >  ZMAX)THEN
          ZMAX=pvals_g(JPT+(JROW-1)*row_length_global)
          IMAXPT=JPT
          IMAXRO=JROW
        ENDIF
      ENDDO

!     CALCULATE MIN

      DO JPT=JPTF,JPTL
        IF(pvals_g(JPT+(JROW-1)*row_length_global) <  ZMIN)THEN
          ZMIN=pvals_g(JPT+(JROW-1)*row_length_global)
        IMINPT=JPT
        IMINRO=JROW
        ENDIF
      ENDDO

      ZLAT = ZLAT + ZDLAT
      END DO ! jrow

!     EVALUATE STATS AND WRITE OUT

      NPTS = JPTL-JPTF+1
      IF(SUMWT /= 0.) WT   = 1.0/(SUMWT*NPTS)
      ZM   = ZM  * WT
      ZMS  = ZMS * WT
      IF(ZMS >  0.)   ZRMS = SQRT(ZMS)

      WRITE(6,62)                                                       &
     &  PNTLAB,KLEV,ZM,ZRMS,ZMAX,IMAXRO,IMAXPT,ZMIN,IMINRO,IMINPT

62    FORMAT(1X,A16,2X,I4,' MEAN=',G12.5,' RMS=',G12.5,                 &
     & ' MAX=',G12.5,' AT (',I4,',',I4,')',                             &
     & ' MIN=',G12.5,' AT (',I4,',',I4,')')

      endif      ! mype == 0

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('MMSPT   ',4)
      IF (lhook) CALL dr_hook('MMSPT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE MMSPT
END MODULE mmspt_mod
