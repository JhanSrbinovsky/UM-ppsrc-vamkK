! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose : Calculate no of observations and no of observation
!            values in AC observation files. These values are
!            used to dimension any observation data arrays in the
!            assimilation code.
!
!  Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!
!  Project Task : P3
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE check_obs_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE CHECK_OBS (OBS_ROW_LENGTH,OBS_P_ROWS,OBS_P_LEVELS,     &
     &                      OBS_Q_LEVELS,OBS_AK,OBS_BK,OBS_LONG_RES,    &
     &                      OBS_LAT_RES,OBS_START_LAT,OBS_START_LONG,   &
     &                      OBS_LAT_PSEUDO_POLE,OBS_LONG_PSEUDO_POLE,   &
     &                      P_LEVELS,Q_LEVELS,P_ROWS,ROW_LENGTH,AK,BK,  &
     &                      REALHD1,REALHD2,REALHD3,REALHD4,            &
     &                      REALHD5,REALHD6,ICODE,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParVars
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
      INTEGER OBS_ROW_LENGTH,OBS_P_ROWS,OBS_P_LEVELS,OBS_Q_LEVELS
      REAL OBS_LONG_RES,OBS_LAT_RES,OBS_START_LAT
      REAL OBS_START_LONG,OBS_LAT_PSEUDO_POLE,OBS_LONG_PSEUDO_POLE
      INTEGER P_LEVELS      !  No of model levels
      INTEGER Q_LEVELS      !  No of model wet levels
      INTEGER P_ROWS        !  No of model rows
      INTEGER ROW_LENGTH    !  No of points on row
      REAL    OBS_AK(P_LEVELS),OBS_BK(P_LEVELS)
      REAL    AK(P_LEVELS)  !  Vertical grid
      REAL    BK(P_LEVELS)
      REAL REALHD1,REALHD2  !  Horizontal grid
      REAL REALHD3,REALHD4
      REAL REALHD5,REALHD6
      INTEGER ICODE           !  Return code
      CHARACTER(LEN=256) CMESSAGE  !  Error message if ICODE > 0
      INTEGER ROW_LENGTH_GLOBAL,P_ROWS_GLOBAL
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER JLEV

      REAL P1,P2
      LOGICAL LNER

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      ! Implied function
      LNER(P1,P2) = ((ABS(P1-P2))  >   (1.E-6*ABS(P1+P2)))

!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CHECK_OBS',zhook_in,zhook_handle)

      ROW_LENGTH_GLOBAL=glsize(1,fld_type_p)
      P_ROWS_GLOBAL=glsize(2,fld_type_p)
      PRINT *, ' '
      PRINT *, 'NUM_OBS : Checking consistency of level and grid'
      PRINT *, '          between ACOBS file and MODEL'
      PRINT *, ' '
      ICODE = 0
      IF (OBS_ROW_LENGTH /= ROW_LENGTH_GLOBAL) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'row_length_acobs      = ',obs_row_length
        PRINT *, 'row_length_global      = ',row_length_global
        ICODE = 1
      ENDIF
      IF (OBS_P_ROWS /= P_ROWS_GLOBAL) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'p_rows_acobs          = ',obs_p_rows
        PRINT *, 'p_rows_global          = ',p_rows_global
        ICODE = 1
      ENDIF

      IF (OBS_P_LEVELS /= P_LEVELS) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'p_levels_acobs        = ',obs_p_levels
        PRINT *, 'p_levels_model        = ',p_levels
        ICODE = 1
      ENDIF
      IF (OBS_Q_LEVELS /= Q_LEVELS) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'q_levels_acobs        = ',obs_q_levels
        PRINT *, 'q_levels_model        = ',q_levels
        ICODE = 1
      ENDIF

      DO JLEV=1,P_LEVELS
        IF ((OBS_AK(JLEV) /= AK(JLEV)) .OR.                             &
     &     (OBS_BK(JLEV) /= BK(JLEV))) THEN
          PRINT *, 'Inconsistency found:'
          PRINT *, 'Level ',JLEV
          PRINT *, 'ak_acobs            = ',obs_ak(jlev)
          PRINT *, 'ak_model            = ',ak(jlev)
          PRINT *, 'bk_acobs            = ',obs_bk(jlev)
          PRINT *, 'bk_model            = ',bk(jlev)
          IF (LNER(OBS_AK(JLEV),AK(JLEV)) .OR.                          &
     &        LNER(OBS_BK(JLEV),BK(JLEV))) THEN
            ICODE=1
          ELSE
            PRINT *, 'But inconsistency is within 32 bit accuracy'
          ENDIF
        ENDIF
      ENDDO

      IF (OBS_LONG_RES /= REALHD1) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'long_res_acobs        = ',obs_long_res
        PRINT *, 'long_res_model        = ',realhd1
       IF (LNER(OBS_LONG_RES,REALHD1)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_LAT_RES /= REALHD2) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'lat_res_acobs         = ',obs_lat_res
        PRINT *, 'lat_res_model         = ',realhd2
       IF (LNER(OBS_LAT_RES,REALHD2)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_START_LAT /= REALHD3) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'start_lat_acobs       = ',obs_start_lat
        PRINT *, 'start_lat_model       = ',realhd3
       IF (LNER(OBS_START_LAT,REALHD3)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_START_LONG /= REALHD4) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'start_long_acobs      = ',obs_start_long
        PRINT *, 'start_long_model      = ',realhd4
       IF (LNER(OBS_START_LONG,REALHD4)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_LAT_PSEUDO_POLE /= REALHD5) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'lat_pseudo_pole_acobs = ',obs_lat_pseudo_pole
        PRINT *, 'lat_pseudo_pole_model = ',realhd5
      IF (LNER(OBS_LAT_PSEUDO_POLE,REALHD5)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF
      IF (OBS_LONG_PSEUDO_POLE /= REALHD6) THEN
        PRINT *, 'Inconsistency found:'
        PRINT *, 'long_pseudo_pole_acobs= ',obs_long_pseudo_pole
        PRINT *, 'long_pseudo_pole_model= ',realhd6
      IF (LNER(OBS_LONG_PSEUDO_POLE,REALHD6)) THEN
          ICODE=1
       ELSE
          PRINT *, 'But inconsistency is within 32 bit accuracy'
       ENDIF
      ENDIF

      IF (ICODE == 1) THEN
        CMESSAGE = 'NUM_OBS : Failure in NUM_OBS - '//&
            'mismatch between model and acobs level/grid information'
      ELSE
        PRINT *, 'NUM_OBS : ACOBS and MODEL are consistent'
      ENDIF

      IF (lhook) CALL dr_hook('CHECK_OBS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE CHECK_OBS
END MODULE check_obs_mod
