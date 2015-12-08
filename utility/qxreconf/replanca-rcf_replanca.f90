! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read in Ancillary Fields
!
!  Subroutine replanca_rcf_replanca  - Read in Ancillary Fields
!
! Description:
!   Read in Ancillary Fields
!
! Method:
!   This routine reads in the required ancillary fields. Any time
!   interpolation is done as required. Both 360 and 365 day calendars
!   are catered for. The ANCILmaster files are used in the
!   reconfiguration.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!

  SUBROUTINE replanca_rcf_replanca(I_YEAR,I_MONTH,I_DAY,I_HOUR,    &
                      I_MINUTE,I_SECOND,                           &
                      P_FIELD,P_ROWS,U_FIELD,V_FIELD,RR_FIELD,     &
                      LAND_FIELD,                                  &
                      D1,LAND,                                     &
                      ICE_FRACTION,TSTAR,FLAND,                    &
                      TSTAR_LAND_CTILE,TSTAR_SEA_CTILE,            &
                      TSTAR_SICE_CTILE,                            &
                      TSTAR_ANOM,                                  &
                      NS_SPACE,FIRST_LAT,                          &
                      LEN1_LOOKUP,LEN_FIXHD,LEN_INTHD,             &
                      LEN_D1,FIXHD,INTHD,                          &
                      LOOKUP,RLOOKUP,LOOKUP_START,                 &
                      NLOOKUPS,                                    &
                      ICODE,CMESSAGE)

  USE ancil_mod, ONLY :                                             &
      Anc_File,       ancFiles,                                     &
      Anc_Record,     ancRecs,                                      &
      AncF_UnitNo,    ancil_add,                                    &
      StashAncil,     Levels,                                       &
      Lookup_Step,    NLookup

  USE Locate_Anc_mod, ONLY :                                        &
      Locate_Anc_Field,                                             &
      Locate_Anc_File

  USE decomp_params, ONLY :                                         &
      Decomp_rcf_output

  USE um_input_control_mod, ONLY :                                  &
      lcal360

  USE ancilcta_namelist_mod, ONLY :                                 &
      lamipii,                                                      &
      l_sstanom

  USE UM_ParVars, ONLY :                                            &
      mype

  USE PrintStatus_mod, ONLY :                                       &
      PrintStatus, PrStatus_Diag

  USE Rcf_CntlAtm_Mod, ONLY :                                       &
      ltleads 
      
  USE ukca_option_mod, ONLY :                                       &
      l_ukca
     
  USE switches, ONLY: l_ctile

  USE IO
     
  USE water_constants_mod, ONLY: tfs,tm
     
      USE lookup_addresses

      Use Rcf_Exppx_Mod, Only :      &
          Rcf_Exppx
      
      Use Rcf_Ppx_Info_Mod, Only :    &
          STM_Record_Type

      Use Rcf_Grid_Type_Mod, Only : &
          Output_Grid

      Use Rcf_Level_Code_Mod, Only : &
          Rcf_Level_Code


      USE t_int_mod, ONLY: t_int
      USE t_int_c_mod, ONLY: t_int_c
      USE science_fixes_mod, ONLY: l_error_ancil_struct
      IMPLICIT NONE

! Arguments
  Integer :: I_YEAR          ! Curent Model Time
  Integer :: I_MONTH         !   "      "     "
  Integer :: I_DAY           !   "      "     "
  Integer :: I_HOUR          !   "      "     "
  Integer :: I_MINUTE        !   "      "     "
  Integer :: I_SECOND        !   "      "     "

  Integer :: P_FIELD         ! Size of horizontal fields
  Integer :: P_ROWS          !
  Integer :: U_FIELD         !   "  "      "         "
  Integer :: V_FIELD         !   "  "      "         "
  Integer :: RR_FIELD        !   "  "      "         "
  Integer :: Land_Field      !   "  "      "         "

  Integer :: NLOOKUPS        ! Number of lookup tables
  Integer :: Len_FixHd       ! Length of Fixed Header
  Integer :: Len_IntHd       ! Length of Integer Header
  Integer :: Len1_Lookup     ! First dimension of Lookup
  Integer :: LEN_D1          ! Size of primary data array

  Integer FIXHD(LEN_FIXHD,ancFiles)    ! Anc fixed headers
  Integer INTHD(LEN_INTHD,ancFiles)    ! Anc Integer Headers
  Integer LOOKUP(LEN1_LOOKUP,NLOOKUPS) ! Anc Lookup Tables
  Integer LOOKUP_START(ancFiles)       ! Start of lookup tables for
                                       ! anc files opened.
  Integer ICODE                        ! Return Code
  Integer IPOS_111                     ! Land fraction
!                                      ! ancil position


  Real NS_Space              ! NS latitude spacing
  Real First_Lat             ! Latitude of first gridpoint
  Real D1(Len_D1)            ! INOUT Data array to hold fields
                             !       except TStar and Ice Fraction
  Real Ice_Fraction(P_Field) ! INOUT Ice Fraction, updated if
                             !       requested
  Real TStar (P_Field)       ! INOUT T Star, updated if requested
  Real TStar_Land_Ctile (P_Field)
!                            ! INOUT T*_land, updated if requested
  Real TStar_Sea_Ctile (P_Field)
!                            ! INOUT T*_sea, updated if requested
  Real TStar_Sice_Ctile (P_Field)
!                            ! INOUT T*_sice, updated if requested
  Real TStar_Anom(P_Field)   ! INOUT SST Anomaly, formed in Recon.
                             !       Added in model run if requested

  Real RLookup(Len1_Lookup,NLookups) ! REAL copy of Lookup Table

  Logical Land (P_Field)     ! Land mask

  Character (Len=80) :: CMessage

! Local Variables

! Buffers to hold values of ancillary data for time interpolation.
! Field of ancillary data held prior to selective updating.
  Real, Dimension (:), Allocatable :: ANCIL1
  Real, Dimension (:), Allocatable :: ANCIL2
  Real, Dimension (:), Allocatable :: ANCIL_DATA

  Real SNOW_CHANGE(P_FIELD)   ! Fractional time of change of
                              ! snow cover
  Real ICE_EXTENT(P_FIELD,2)  ! Fractional time of change
                              ! of ice cover
  Real PRES_VALUE(P_FIELD)    ! Prescribed value of data when
                              ! controlling field is zero.
  Real NO_ICE_EXTENT(P_FIELD) ! Indicator for no sea ice
                              ! =0 if ice cover

  Real TStar_Land (P_Field)  ! T*_land, updated if requested
  Real TStar_Sea (P_Field)   ! T*_sea, updated if requested
  Real TStar_Sice (P_Field)  ! T*_sice, updated if requested
  Real TStar_Ssi (P_Field)   ! T*_ssi, updated if requested
  Real Flandg (P_Field)      ! Land fraction in gridbox
  Real Fland (Land_Field)    ! Land fraction in gridbox
!                            ! (on land-only points).
  Integer I,J,L              ! Loop indices

  Integer I1,I2,I3,II       ! Array indices
  Integer ID                !
  Integer IM                !
  Integer IY                !
  Integer FIELD             ! Current Ancil Ref Number.
  Integer FILE              ! Anc File number
  Integer IREC              ! Anc record number

  Integer INTERVAL          ! Interval between data times
  Integer STEP              ! Number of data times skipped.
  Integer MONTHS            ! )Used in calculation of position
  Integer HOURS             ! )of data required
  Integer PERIOD            ! Period of periodic data
  Integer START_MONTH       !
  Integer LEVEL             ! Loop index for levels
  Integer ANCIL_REF_DAYS    ! Ancil.reference time in whole days
  Integer ANCIL_REF_SECS    ! Ancil.reference time in extra seconds
  Integer DAY,SEC           ! Times relative to reference time
  Integer LEN
  Integer ROW_LENGTH
  Integer I_YEAR1          ! Copy of Curent Model Time year
  Integer I_MONTH1         !   "      "     "          month
  Integer I_DAY1           !   "      "     "          day
  Integer I_HOUR1          !   "      "     "          hour

  Integer len_anc_env_var  !  Length of env var
  Integer anc_file_open    !  Stores which ancillary file is open.
  Integer ipos_27          !  Pos of anc record for Sea Ice Frac
  Integer ipos_28          !  Pos of anc record for SST
  Integer ipos_30          !  Pos of anc record for surf u-current
  Integer ipos_31          !  Pos of anc record for surf v-current

  Integer FILEANCIL(ancRecs) ! Anc File No for anc record

  Integer FIELD_SIZE       !  Stores field size allowing use of
                           !  P_FIELD or RR_FIELD in calls.

  Logical LINTERPOLATE      ! Indicates whether time
                            ! interpolation needed.
  Logical LT_INT_C          ! Indicates use of controlled time
                            ! interpolation
  Logical LMISMATCH         ! Used in header checks
  Logical SINGLE_TIME       ! Indicates that only one time is
                            ! available in data set
  Logical PERIODIC          ! Data set is periodic
  Logical REGULAR           ! Interval between data times in
                            ! dataset is regular in model timesteps.
  Logical LICE_DEPTH        ! T : Field is Ice Depth
  Logical LICE_FRACTION     ! T : Field is Ice Fraction
  Logical LSNOW_DEPTH       ! T : Field is Snow depth

  Logical UPDATE (ancRecs)  ! T : Anc Field to be updated

  Logical Sea (P_Field)      ! Sea mask

  Real ZERO         !
  Real TIME1        !  Times if data used in time interpolation
  Real TIME2        !
  Real TIME         !  Target time for time interpolation
  Real LAT_P        !  Latitude of point


! STASHmaster entry
      Type (STM_Record_Type), Pointer :: STM_Rec
      INTEGER :: bot_level


! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!   0. Set up local arrays for ANCIL data to required size.
  FIELD_SIZE=MAX(P_FIELD,RR_FIELD)

  IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
    WRITE(6,'(A,I12)') "Rcf_replanca_rcf_replanca: chosen FIELD_SIZE=",&
          FIELD_SIZE
  END IF

  ALLOCATE ( ANCIL1(FIELD_SIZE) )
  ALLOCATE ( ANCIL2(FIELD_SIZE) )
  ALLOCATE ( ANCIL_DATA(FIELD_SIZE) )

!L  1.  Initialisation for atmosphere

  ICODE=0
  CMESSAGE=' '
  anc_file_open = 0

! Initialise ANCIL1/2. Includes Halos for MPP runs.
  ANCIL1(:)=0.0
  ANCIL2(:)=0.0

! Read in fractional land field first:
! use ANCIL_DATA as temporary storage

  IF (L_CTILE) THEN

    CALL locate_anc_field (111, ipos_111)

    IF(anc_record(ipos_111) % anc_field_read == 1)THEN !read ancil

      CALL locate_anc_file(anc_record(ipos_111) % anc_file_number, FILE)
      len_anc_env_var = LEN_TRIM ( anc_file(file) % anc_env_var )

      CALL file_open (AncF_UnitNo,                                 &
                       anc_file(file) % anc_env_var,               &
                       len_anc_env_var,0,0,icode)

! DEPENDS ON: rcf_readflds
      CALL Rcf_ReadFlds                                            &
                  (AncF_UnitNo,1,NLOOKUP(ipos_111),                &
                   LOOKUP(1,LOOKUP_START(FILE)),                   &
                   LEN1_LOOKUP,ANCIL_DATA,P_FIELD,FIXHD(1,FILE),   &
                   ICODE,CMESSAGE)
  
      CALL file_close (AncF_UnitNo,                                &
                       anc_file(file) % anc_env_var,               &
                       len_anc_env_var,0,0,icode)
 
      WRITE(6,'(A)') 'READ IN LAND FRACTION'

      DO I=1,P_FIELD
        D1(ANCIL_ADD(IPOS_111)+I-1)=ANCIL_DATA(I)
      END DO

      DO I=1,P_FIELD
        FLANDG(I)=0.0
        IF (LAND(I)) FLANDG(I)=ANCIL_DATA(I)
      END DO

    ELSE             ! Land frac already read in from input dump

      L=0
      DO I=1,P_FIELD
        FLANDG(I)=0.0
        IF (LAND(I)) THEN
          L=L+1
          FLANDG(I)=FLAND(L)
        END IF
      END DO
    END IF

    DO I=1,P_FIELD
! If land or sea fraction is less than machine tolerance print warning
      IF (LAND(I).AND.FLANDG(I) <  EPSILON(1.0)) THEN
        WRITE(6,'(A)')'*****************WARNING********************'
        WRITE(6,'(A)')'LAND FRACTION IS LESS THAN MACHINE TOLERANCE'
      END IF
      IF (.NOT.LAND(I).AND.1.0-FLANDG(I) <  EPSILON(1.0)) THEN
        WRITE(6,'(A)')'*****************WARNING********************'
        WRITE(6,'(A)')'SEA FRACTION IS LESS THAN MACHINE TOLERANCE'
      END IF

      IF (FLANDG(I) <= 0.0.AND.LAND(I)) THEN
        WRITE(6,'(A)')'*ERROR* a) LAND FRAC & LAND MASK INCONSISTENT'
        ICODE = 800
        CMESSAGE='replanca_rcf_replanca:ERR:LAND FRAC & MASK ARE ' &
        //'INCONSISTENT'
      END IF
      IF (FLANDG(I) >  0.0.AND..NOT.LAND(I)) THEN
        WRITE(6,'(A)')'*ERROR* b) LAND FRAC & LAND MASK INCONSISTENT'
        ICODE = 801
        CMESSAGE='replanca_rcf_replanca:ERR:LAND FRAC & MASK ARE ' &
        //'INCONSISTENT'
      END IF

    END DO

  ELSE                     ! Not coastal tiling:
    DO I=1,P_FIELD
      IF (LAND(I)) THEN
        FLANDG(I)=1.0
      ELSE
        FLANDG(I)=0.0
      END IF
    END DO
  END IF                   ! End of coastal tiling loop


  DO I=1,P_FIELD
    IF (FLANDG(I) <  1.0) THEN
      SEA(I)=.TRUE.
    ELSE
      SEA(I)=.FALSE.
    END IF
  END DO


! Set up surface temperatures:
  IF (L_CTILE) THEN
    DO I=1,P_FIELD
      TSTAR_LAND(I)=TSTAR_LAND_CTILE(I)
      TSTAR_SEA(I)=TSTAR_SEA_CTILE(I)
      TSTAR_SICE(I)=TSTAR_SICE_CTILE(I)
      IF (ICE_FRACTION(I) <= 0.0) THEN
        TSTAR_SSI(I)=TSTAR_SEA(I)
      ELSE
        TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                &
               +(1.0-ICE_FRACTION(I))*TSTAR_SEA(I)
      END IF
    END DO
  ELSE
    DO I=1,P_FIELD
      TSTAR_LAND(I)=TSTAR(I)
      TSTAR_SSI(I)=TSTAR(I)
    END DO
  END IF

!L  1.1 Set logical UPDATE for each ancillary field independently

! FILEANCIL not required. Remove later.

  DO i=1,ancRecs
    CALL locate_anc_file( anc_record(i) % anc_file_number, FILEANCIL(i) )
    UPDATE(i)    = anc_record(i) % anc_field_read > 0
  END DO

! Initialise for valid time interpolation in reconfiguration mode.
  ANCIL_REF_DAYS = 0
  ANCIL_REF_SECS = 0

!L 1.2 Allow for dependencies between fields
! Sea surface temperature must be updated when sea ice is updated

  CALL locate_anc_field (27, ipos_27)
  CALL locate_anc_field (28, ipos_28)
  CALL locate_anc_field (30, ipos_30)
  CALL locate_anc_field (31, ipos_31)

  UPDATE(ipos_28) = UPDATE(ipos_27) .OR. UPDATE(ipos_28)

! Both surface current components must be updated together

  UPDATE(ipos_30) = UPDATE(ipos_30) .OR. UPDATE(ipos_31)
  UPDATE(ipos_31) = UPDATE(ipos_30)

!L Select method of time interpolation for SST. The interpolation
!L allows for sea ice if ice data is available at the same times
!L as the temperature data. Otherwise linear interpolation is used.

  LT_INT_C=.TRUE.

  IF (UPDATE(ipos_28)) THEN

    IF (FIXHD(10,FILEANCIL(ipos_27)) == 0) LT_INT_C=.FALSE.
    IF (LT_INT_C) THEN
      DO I=21,41
        IF (FIXHD(I,FILEANCIL(ipos_27)) /=                          &
            FIXHD(I,FILEANCIL(ipos_28))) THEN
          LT_INT_C=.FALSE.
          WRITE(6,'(3A)')' WARNING: controlled time interpolation ',&
          'for SST not available: Mismatch in SST and SEA-ICE ',    &
          'ancillary data times in FIXED HEADER'
          WRITE(6,'(A,I4,A,I6)')' position=',I,                     &
                      ' SEA-ICE=',FIXHD(I,FILEANCIL(ipos_27))
          WRITE(6,'(A,I4,A,I6)')' position=',I,                     &
                      ' SST    =',FIXHD(I,FILEANCIL(ipos_28))
        END IF
      END DO
    END IF
  END IF

!L Loop over ancillary fields(atmosphere)

  DO irec=1,ancRecs

    field = anc_record(irec) % ancil_ref_number

    LICE_DEPTH=field == 29  ! required for LAMIPII

    IF (UPDATE(irec)) THEN  ! (1st level IF)

      CALL locate_anc_file(anc_record(irec) % anc_file_number, FILE)

      IF (file  /=  anc_file_open) THEN ! File is not open
            
! Only close file if anc_file_open is greater than 0. (e.g. first open).
        IF (anc_file_open > 0) THEN
          CALL file_close (AncF_UnitNo,                             &
                           anc_file(anc_file_open) % anc_env_var,   &
                           len_anc_env_var,0,0,icode)
        END IF

        len_anc_env_var = LEN_TRIM ( anc_file(file) % anc_env_var )


        CALL file_open (AncF_UnitNo,                                &
                        anc_file(file) % anc_env_var,               &
                        len_anc_env_var,0,0,icode)

        IF (icode /= 0) THEN
          WRITE (6,'(A)') ' problem with file_open.'
          WRITE (6,'(A,I5)') ' icode ',icode
          WRITE (6,'(2A)') ' cmessage ',cmessage
          GO TO 9999
        END IF

!       Record ancillary file opened
        anc_file_open = file

      END IF ! If file not already opened

      IF (LICE_DEPTH.AND.LAMIPII) THEN

! Uses ice fraction set earlier in field loop.
! WARNING this will fail if the order of ancillary fields is ever
! changed so that ice-depth preceeds ice fraction
! Note : For complete sea ice cover
!        Arctic ice depth    = 2m
!        Antarctic ice depth = 1m
! For ice concentrations less than 1. ice depth is 1 or 2 times conc.
! This results in similar values to those from runs using ancillary
! files containing ice depths set to 1 or 2m.

        ROW_LENGTH=P_FIELD/P_ROWS
        DO I=1,P_ROWS
! work out latitude in radians
          LAT_P=FIRST_LAT-NS_SPACE*(I-1)
          DO J=1,ROW_LENGTH
            II=J+(I-1)*ROW_LENGTH
            ANCIL_DATA(II)=0.0
            IF (ICE_FRACTION(II) >  0.0) THEN
              IF (LAT_P >  0.0) THEN   ! Arctic ice depth
                ANCIL_DATA(II)=2.*ICE_FRACTION(II)
              ELSE                     ! Antarctic ice depth
                ANCIL_DATA(II)=1.*ICE_FRACTION(II)
              END IF
            END IF
          END DO
        END DO
!L      Sea ice thickness
!L      Update over all sea points (all sea ice points are the only
!L      ones strictly required, but this cannot be determined easily)

        DO I=1,P_FIELD
          IF (SEA(I)) THEN
            D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
          END IF
        END DO
      ELSE
!     Update required for field

        IF ( mype == 0 ) THEN
          WRITE(6,'(A,I6,2A)')'replanca_rcf_replanca: UPDATE REQUIRED'&
          //' FOR FIELD',FIELD,' : ', anc_record(irec) % anc_name
        END IF
        IF ( FIXHD(10,FILE) < 0 .OR. FIXHD(10,FILE) > 2 ) THEN
          WRITE (6,'(A,I5)') ' file ',file
          WRITE (6,'(A,I5)') ' FIXHD(10,file) ',FIXHD(10,file)
          ICODE = 700 + file
          CMESSAGE = 'replanca_rcf_replanca: Error in fixed '// &
          'header(10) of ancillary file'
          GO TO 9999
        END IF

!L    Check whether more than one data time available in data set

        SINGLE_TIME=FIXHD(10,FILE) == 0

!L    Set default values for time interpolation

        LINTERPOLATE=.TRUE.

        IF (SINGLE_TIME) THEN
          LINTERPOLATE=.FALSE.
        END IF

        IF (FIELD >  9 .AND. FIELD <  19) THEN
          LINTERPOLATE=.FALSE.
        END IF

!L 2.1 Find position of input record

!L    Default settings of search parameters if only one time present

        IF (SINGLE_TIME) THEN
          STEP=0
        ELSE

          PERIODIC=FIXHD(10,FILE) == 2
          REGULAR=.TRUE.

          IF (.NOT. LCAL360) THEN
            REGULAR=FIXHD(35,FILE) == 0.AND.FIXHD(36,FILE) == 0
! i.e. data at intervals of days/hours & non-periodic
            IF (PERIODIC) REGULAR=REGULAR.AND.FIXHD(37,FILE) == 0
! i.e. data at intervals of hours & periodic
          END IF

!         Error checking on time information.

          IF ( FIXHD(35,FILE) < 0 .OR.                            &
               FIXHD(36,FILE) < 0 .OR. FIXHD(36,FILE) > 12 .OR.   &
       REGULAR .AND. (FIXHD(37,FILE) < 0 .OR. FIXHD(37,FILE) > 31 &
          .OR. FIXHD(38,FILE) < 0 .OR. FIXHD(38,FILE) > 24) ) THEN
!           FIXHD(39-40) are not used by replanca_rcf_replanca.
!           FIXHD(35-37) have already been used if not CAL360.
            ICODE = 700 + FIELD
            CMESSAGE = 'replanca_rcf_replanca: Error in validity '&
            //'time interval given in ancillary file'
            RETURN
          END IF

          IF ( FIXHD(21,FILE) < 0 .AND. .NOT. PERIODIC              &
         .OR. .NOT. ( REGULAR .AND. PERIODIC ) .AND.                &
!     If it is REGULAR & PERIODIC more detailed check is applied below
             ( FIXHD(22,FILE) < 0 .OR. FIXHD(22,FILE) > 12 .OR.     &
               FIXHD(23,FILE) < 0 .OR. FIXHD(23,FILE) > 31 .OR.     &
               FIXHD(24,FILE) < 0 .OR. FIXHD(24,FILE) > 24 .OR.     &
               FIXHD(25,FILE) < 0 .OR. FIXHD(25,FILE) > 60 .OR.     &
               FIXHD(26,FILE) < 0 .OR. FIXHD(26,FILE) > 60 ) ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'replanca_rcf_replanca: Error in first '// &
            'validity time given in ancillary file'
            RETURN
          END IF

          IF (.NOT.PERIODIC) THEN

!L        If data taken from full time series of input data.

! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR                &
                         ,I_MINUTE,I_SECOND                          &
                         ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC      &
                         ,LCAL360)

            IF (REGULAR) THEN

!L 2.1.1  Standard cases:360 day calender;
!L 2.1.1  or Gregorian calendar with
!L        interval between data times in days or hours
!L        updating interval may be regular in model timesteps,
!L        or (LGREG_MONTHLY=T) irregular in model timesteps,

              HOURS=SEC/3600+DAY*24
!L FInd time(in hours) of first ancillary data on file
! DEPENDS ON: time2sec
              CALL TIME2SEC(FIXHD(21,FILE),FIXHD(22,FILE),          &
                     FIXHD(23,FILE),FIXHD(24,FILE),                 &
                     FIXHD(25,FILE),FIXHD(26,FILE),                 &
                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,         &
                     LCAL360)
              HOURS=HOURS-SEC/3600-DAY*24

              IF (HOURS <  0) THEN
                ICODE=400+FIELD
                CMESSAGE='replanca_rcf_replanca: Current time '&
                //'precedes start time of data'
                RETURN
              END IF

!L FInd interval(in hours) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+      &
                       FIXHD(37,FILE)*24+FIXHD(38,FILE)

! Do not interpolate in time if data time exactly matches model time

              IF (MOD(HOURS,INTERVAL) == 0) THEN
                LINTERPOLATE=.FALSE.
              END IF

              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE

!L 2.1.2 Gregorian calender;ancillary data interval is in months or
!L       years,which is irregular in model timesteps.
!L original code is inaccurate for this section - corrected code under
!L LAMIPII makes use of dates in lookup headers
!L For a real calendar year the mid-point of each month is different
!L in terms of its hour and day. The old inaccurate method assumes
!L the hour and day are taken from the fixhd values. These are only
!L usually correct for the first month on the ancillary file.

!L FInd interval(in months) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*12+FIXHD(36,FILE)
              MONTHS=I_YEAR*12+I_MONTH
              START_MONTH=FIXHD(21,FILE)*12+FIXHD(22,FILE)
              MONTHS=MONTHS-START_MONTH
!  Check for time within month
              IF (LAMIPII) THEN   ! corrected code uses pp header
                STEP=MONTHS/INTERVAL
                I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*STEP
                I1=I2+LOOKUP_START(FILE)-1
! Check against day and hour of actual lookup header not first field
                IF ((I_DAY*24+I_HOUR) <                             &
                    (LOOKUP(3,I1)*24+LOOKUP(4,I1))) THEN
                  MONTHS=MONTHS-1
                END IF
              ELSE           ! old less accurate code uses FIXHD
                IF ((I_DAY*24+I_HOUR) <                             &
                    (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                  MONTHS=MONTHS-1
                END IF
              END IF ! LAMIPII

              IF (MONTHS <  0) THEN
                ICODE=400+FIELD
                CMESSAGE='replanca_rcf_replanca: Current time '&
                //'precedes start time of data'
                RETURN
              END IF


              STEP=MONTHS/INTERVAL

              IF (LAMIPII) THEN       ! corrected code
                TIME=REAL(SEC)/3600+REAL(DAY*24)
! correct calculation of dates uses lookup table dates not fixhd date
                I2=NLOOKUP(Irec)+LOOKUP_STEP(irec)*STEP
                I1=I2+LOOKUP_START(FILE)-1
                I_YEAR1=lookup(1,i1)
                I_MONTH1=lookup(2,i1)
                I_DAY1=lookup(3,i1)
                I_HOUR1=lookup(4,i1)
! DEPENDS ON: time2sec
                CALL TIME2SEC(I_YEAR1,I_MONTH1,I_DAY1,I_HOUR1,      &
                    FIXHD(25,FILE),FIXHD(26,FILE),                  &
                    ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,          &
                    LCAL360)
                TIME1=REAL(SEC)/3600+REAL(DAY*24)
! I1+LOOKUP_STEP(irec) correct pointer to next field as multi-level fields 
! are possible.
                IF (i1+lookup_step(irec) < nlookups) THEN
              ! Second check next lookup is the same field.
                  IF (lookup(item_code,i1) ==                       &
                      lookup(item_code,i1+lookup_step(irec))) THEN
                    i_year1=lookup(1,i1+lookup_step(irec))
                    i_month1=lookup(2,i1+lookup_step(irec))
                    i_day1=lookup(3,i1+lookup_step(irec))
                    i_hour1=lookup(4,i1+lookup_step(irec))
                  ELSE
                    icode = 500
                    WRITE(cmessage,'(A,I6)') 'REPLANCA: error finding'&
                    //' next lookup entry for field:', field
                    GO TO 9999
                  END IF
                ELSE
                  icode = 500
                  WRITE(cmessage,'(A,I6)') 'REPLANCA: error due to '&
                  //'lookup out of range for field:', field
                  GO TO 9999
                END IF

! DEPENDS ON: time2sec
                CALL TIME2SEC(I_YEAR1,I_MONTH1,I_DAY1,I_HOUR1,      &
                    FIXHD(25,FILE),FIXHD(26,FILE),                  &
                    ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,          &
                    LCAL360)
                TIME2=REAL(SEC)/3600+REAL(DAY*24)

              ELSE   ! LAMIPII test - old inaccurate code using FIXHD
! NB INTERVAL may be > 1 month
                MONTHS=STEP*INTERVAL
! Calculate data times for time interpolation
                TIME=REAL(SEC)/3600+REAL(DAY*24)
                IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
                IY=FIXHD(21,FILE)+(MONTHS+FIXHD(22,FILE)-1)/12
! DEPENDS ON: time2sec
                CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),  &
                    FIXHD(25,FILE),FIXHD(26,FILE),                  &
                    ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,          &
                    LCAL360)
                TIME1=REAL(SEC)/3600+REAL(DAY*24)
                IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
                IY=FIXHD(21,FILE)+(MONTHS+INTERVAL+FIXHD(22,FILE)-1)/12
! DEPENDS ON: time2sec
                CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),  &
                    FIXHD(25,FILE),FIXHD(26,FILE),                  &
                    ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,          &
                    LCAL360)
                TIME2=REAL(SEC)/3600+REAL(DAY*24)
              END IF     ! end LAMIPII test

! Do not interpolate in time if data time exactly matches model time

              IF (TIME == TIME1) THEN
                LINTERPOLATE=.FALSE.
              END IF

            END IF ! End of REGULAR/not REGULAR

          ELSE  ! PERIODIC data

!L 2.2   If data is taken from ancillary periodic data.

! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR,             &
                          I_MINUTE,I_SECOND,                       &
                          ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,   &
                          LCAL360)

            IF (REGULAR) THEN
!L 2.2.1 Standard cases:1) 360 day calender, with allowed periods of
!L       1 day, 1 month or 1 year;
!L
!L       2) Gregorian calender with update in hours,and period of
!L       data 1 day.
!L
!L       For both updating interval and number of
!L       data times to be skipped in data set calculated in hours.

              HOURS=SEC/3600+DAY*24
              INTERVAL=FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+      &
                       FIXHD(37,FILE)*24+FIXHD(38,FILE)

              PERIOD=INTHD(3,FILE)*INTERVAL

!L   Do not allow non-standard periods
              IF (LCAL360) THEN
                IF (PERIOD/=8640.AND.PERIOD/=720.AND.PERIOD/=24) THEN
                  ICODE=600+FIELD
                  CMESSAGE='replanca_rcf_replanca: Non-standard '&
                  //'period for periodic data'
                  RETURN
                END IF
              ELSE
                IF (PERIOD /= 24) THEN
                  ICODE=600+FIELD
                  CMESSAGE='replanca_rcf_replanca: Non-standard '&
                  //'period for periodic data'
                  RETURN
                END IF
              END IF
              IF (PERIOD == 24) THEN
! Ancillary data interval in hour(s), period is 1 day

                IY=I_YEAR
                IM=I_MONTH
                ID=I_DAY
                IF (I_HOUR < FIXHD(24,FILE)) HOURS=HOURS+24

              ELSE IF (PERIOD == 720) THEN
! Ancillary data interval in day(s) or hours , period is 1 month

                IY=I_YEAR
                IM=I_MONTH
                ID=FIXHD(23,FILE)
                IF ((I_DAY*24+I_HOUR) <                            &
                    (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                  HOURS=HOURS+720
                END IF

              ELSE IF (PERIOD == 8640) THEN
! Ancillary data interval in month(s)or days or hours, period is 1 year

                IY=I_YEAR
                IM=FIXHD(22,FILE)
                ID=FIXHD(23,FILE)
                IF ((I_MONTH*720+I_DAY*24+I_HOUR)<(FIXHD(22,FILE)*720 &
                              +FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                  HOURS=HOURS+8640
                END IF 

              END IF

! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,ID,FIXHD(24,FILE),                &
                       FIXHD(25,FILE),FIXHD(26,FILE),               &
                       ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,       &
                       LCAL360)
              HOURS=HOURS-SEC/3600-DAY*24

! Do not interpolate in time if data time exactly matches model time

              IF (MOD(HOURS,INTERVAL) == 0) THEN
                LINTERPOLATE=.FALSE.
              END IF
              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE  ! non regular case

!L 2.2.2 Gregorian calender,and data interval is in months,
!L       period is 1 year
!L       Updating interval and number of data times to be skipped
!L       calculated in months.

              TIME=REAL(SEC)/3600+REAL(DAY*24)
              INTERVAL=FIXHD(36,FILE)+FIXHD(35,FILE)*12
              PERIOD=INTHD(3,FILE)*INTERVAL
              IF (PERIOD /= 12) THEN
                ICODE=600+FIELD
                CMESSAGE='replanca_rcf_replanca: Non-standard '// &
                'period for periodic data'
                RETURN
              END IF
!  Difference between date now (month) & first date ancil file (month)
              MONTHS=I_MONTH-FIXHD(22,FILE)

              IF (LAMIPII) THEN ! correct code to use lookup header dates
! Correctly use day and hour from lookup header not fixhd which
! contains values for first field on ancillary file only.
                step=months/INTERVAL
                I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*step
                I1=I2+LOOKUP_START(FILE)-1
!  Check for time within month - using ppheader information
                IF ((I_DAY*24+I_HOUR)<(lookup(3,i1)*24+lookup(4,i1))) THEN
                  MONTHS=MONTHS-1
                END IF
                IF (MONTHS <  0) THEN
                  MONTHS=MONTHS+12
                END IF
! recalculate STEP
                STEP=MONTHS/INTERVAL
! NB INTERVAL qmay be > 1 month
                MONTHS=STEP*INTERVAL
                IY=I_YEAR
                IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
                IF (IM >  I_MONTH) IY=IY-1
                I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*STEP
                I1=I2+LOOKUP_START(FILE)-1
! DEPENDS ON: time2sec
                CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),      &
                    FIXHD(25,FILE),FIXHD(26,FILE),                  &
                    ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
                TIME1=REAL(SEC)/3600+REAL(DAY*24)
!  Calculate  TIME2 for second ancillary data time
!  set IY correctly for time interpolation calculations
                IY=I_YEAR
                IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
                IF (IM <  I_MONTH) IY=IY+1
                I1=(IM-1)/INTERVAL
                I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*I1
                I1=I2+LOOKUP_START(FILE)-1
! DEPENDS ON: time2sec
                CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),      &
                    FIXHD(25,FILE),FIXHD(26,FILE),                  &
                    ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
                TIME2=REAL(SEC)/3600+REAL(DAY*24)

              ELSE   ! original code inaccurate use of FIXHD dates
!  Check for time within month
                IF ((I_DAY*24+I_HOUR) <                             &
                    (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                  MONTHS=MONTHS-1
                END IF
                IF (MONTHS <  0) THEN
                  MONTHS=MONTHS+12
                END IF

                STEP=MONTHS/INTERVAL
! NB INTERVAL may be > 1 month
                MONTHS=STEP*INTERVAL
!  Calculate TIME1 for first ancillary data time
!  set IY correctly for time interpolation calculations
                IY=I_YEAR
                IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
                IF (IM >  I_MONTH) IY=IY-1
! DEPENDS ON: time2sec
                CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),  &
                    FIXHD(25,FILE),FIXHD(26,FILE),                  &
                    ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,          &
                    LCAL360)
                TIME1=REAL(SEC)/3600+REAL(DAY*24)
!  Calculate  TIME2 for second ancillary data time
!  set IY correctly for time interpolation calculations
                IY=I_YEAR
                IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
                IF (IM <  I_MONTH) IY=IY+1
! DEPENDS ON: time2sec
                CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),  &
                    FIXHD(25,FILE),FIXHD(26,FILE),                  &
                    ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,          &
                    LCAL360)
                TIME2=REAL(SEC)/3600+REAL(DAY*24)
              END IF  ! end LAMIPII test

! Do not interpolate in time if data time exactly matches model time

              IF (TIME == TIME1) THEN
                LINTERPOLATE=.FALSE.
              END IF

            END IF  ! regular/non-regular

          END IF  ! non-periodic/periodic

        END IF ! singletime/non-singletime

!L 2.3   Check STASH Code

        I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*STEP

        I1=LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)-1)

        LMISMATCH=.FALSE.

        IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
          WRITE(6,'(A,I12)') ' Information used in checking ancillary '&
          //'data set: position of lookup table in dataset:',I2
          WRITE(6,'(A,I12)') ' Position of first lookup table referring'&
          //' to data type ',NLOOKUP(irec)
          WRITE(6,'(2(A,I8))')' Interval between lookup tables referring'&
          //' to data type ', LOOKUP_STEP(irec),' Number of steps', STEP
          WRITE(6,'(2(A,I8))') ' STASH code in dataset ',I1,          &
                              '  STASH code requested ',STASHANCIL(irec)
          WRITE(6,'(A,I12)') '''Start'' position of lookup tables for '&
          //'dataset in overall lookup array ' ,LOOKUP_START(FILE)
        END IF

        IF (I1 /= STASHANCIL(irec)) THEN
          WRITE(6,'(A,3I8)')' I1,STASHANCIL(irec)',I1,STASHANCIL(irec),irec
          LMISMATCH=.TRUE.
        END IF

!L Error exit if checks fail

        IF (LMISMATCH) THEN
          ICODE=200+FIELD
          CMESSAGE='replanca_rcf_replanca: PP HEADERS ON ANCILLARY '//&
          'FILE DO NOT MATCH'
          RETURN
        END IF

        IF (LINTERPOLATE.AND..NOT.SINGLE_TIME) THEN
!L Check time interpolation factors
          IF (TIME < TIME1 .OR. TIME > TIME2) THEN
            WRITE(6,'(A)')' Information used in interpolation/replacement:'
            WRITE(6,'(A,I8)')' Time of first data=', TIME1
            WRITE(6,'(A,I8)')' Validity Time for update=', TIME
            WRITE(6,'(A,I8)')' Time of second data=', TIME2

            ICODE=500+FIELD
            CMESSAGE='replanca_rcf_replanca: TIME INTERPOLATION ERROR'
            RETURN
          END IF
        END IF

!L 3   Loop over levels of ancillary data for field I
!L Reset pointer for dataset

!L Includes loop over X and Y components of surface currents

        LICE_FRACTION=FIELD == 27
        LSNOW_DEPTH=FIELD == 9
        LICE_DEPTH=FIELD == 29

! Find whether we need to handle the extra level due to ancillary not supporting
! theta level 0.
         STM_Rec => Rcf_Exppx( 1, stashancil(irec)/1000, &
                               MOD(stashancil(irec),1000) )
         
         IF ( STM_Rec % lb_code > 0) THEN
           CALL Rcf_Level_Code( STM_Rec % lb_code, bot_level, Output_Grid )
         ELSE
           bot_level = 1
         END IF

! Force bottom level to be 1 (if not 0)
         IF (bot_level /= 0) THEN
           bot_level = 1
         END IF

         ! Only want to loop upto levels in ancillary.
        DO LEVEL=1,LEVELS(IREC)

!L Do not go through loop for ice edge or snow edge

          IF (.NOT.((LICE_FRACTION.OR.LSNOW_DEPTH).AND.LEVEL == 2)) THEN

!       Check to see if field is one of the River Routing ones.
            IF ( FIELD  ==  124 .OR.                                    &
                FIELD  ==  125 .OR.                                     &
                FIELD  ==  126 ) THEN
              
              IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
                WRITE (6,'(A)') 'Rcf_replanca_rcf_replanca: '&
                    //'Resetting Ancil field size to RR_FIELD'
              END IF

              FIELD_SIZE=RR_FIELD

            ELSE ! If not River Routing, then assume field size is PFIELD

              FIELD_SIZE=P_FIELD

            END IF

!L 3.1 Read data for single level of ancillary field.
            
            IF (.NOT.LICE_FRACTION) THEN
! AMIPII case ice depth field not read from ancillary file
              IF (.NOT.(LICE_DEPTH.AND.LAMIPII)) THEN

! DEPENDS ON: rcf_readflds
                CALL Rcf_ReadFlds                                   &
                    (AncF_UnitNo,1,I2,LOOKUP(1,LOOKUP_START(FILE)), &
                    LEN1_LOOKUP,ANCIL1,FIELD_SIZE,FIXHD(1,FILE),    &
                    ICODE,CMESSAGE)

              END IF

              IF (ICODE > 0) THEN
                ICODE=FIELD+100
                CMESSAGE='replanca_rcf_replanca :I/O ERROR '
                GO TO 9999
              END IF
              
            ELSE

!L If ice-fraction,read fractional time field as well
!L       UNLESS IT IS A SINGLE TIME FIELD
!L If snow-depth,read fractional time field as well only if time
!L interpolation required.

              IF (.NOT.SINGLE_TIME.AND..NOT.LAMIPII) THEN
                IF (LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)) == 38) THEN

! DEPENDS ON: rcf_readflds
                  CALL Rcf_ReadFlds                                   &
                      (AncF_UnitNo,2,I2,LOOKUP(1,LOOKUP_START(FILE)), &
                      LEN1_LOOKUP,ICE_EXTENT,FIELD_SIZE,FIXHD(1,FILE),&
                      ICODE,CMESSAGE)

                  IF (ICODE > 0) THEN
                    ICODE=FIELD+100
                    CMESSAGE='replanca_rcf_replanca :I/O ERROR '
                    GO TO 9999
                  END IF

                ELSE
                  ICODE=FIELD+100
                  CMESSAGE='replanca_rcf_replanca :ICE CHANGE DATA MISSING'
                  GO TO 9999
                END IF
              ELSE    ! single time or LAMIPII - ie no time change field

! DEPENDS ON: rcf_readflds
                CALL Rcf_ReadFlds                                   &
                    (AncF_UnitNo,1,I2,LOOKUP(1,LOOKUP_START(FILE)), &
                    LEN1_LOOKUP,ICE_EXTENT,FIELD_SIZE,FIXHD(1,FILE),&
                    ICODE,CMESSAGE)

                IF (ICODE >  0) THEN
                  ICODE=FIELD+100
                  CMESSAGE='replanca_rcf_replanca :I/O ERROR '
                  GO TO 9999
                END IF
              END IF
            END IF

            IF (LSNOW_DEPTH.AND.LINTERPOLATE) THEN
              IF (LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)) == 27) THEN

! DEPENDS ON: rcf_readflds
                CALL Rcf_ReadFlds                                      &
                    (AncF_UnitNo,1,I2+1,                               &
                    LOOKUP(1,LOOKUP_START(FILE)),                      &
                    LEN1_LOOKUP,SNOW_CHANGE,FIELD_SIZE,FIXHD(1,FILE),  &
                    ICODE,CMESSAGE)

                IF (ICODE >  0) THEN
                  ICODE=FIELD+100
                  CMESSAGE='replanca_rcf_replanca :I/O ERROR '
                  GO TO 9999
                END IF

              ELSE
                ICODE=FIELD+100
                CMESSAGE='replanca_rcf_replanca :SNOW CHANGE DATA MISSING'
                GO TO 9999
              END IF
            END IF

!L If sea surface temperature or other ice fields, read ice fraction
!L and fractional time field if not already pressent and if required
!L by time interpolation.  

            IF (FIELD == 29.OR.(FIELD == 28.AND.LT_INT_C) ) THEN
              
              IF (.NOT.UPDATE(ipos_27)) THEN
                I3 = NLOOKUP(ipos_27) + LOOKUP_STEP(ipos_27)*STEP     &
                    + LOOKUP_START(FILEANCIL(ipos_27))

                IF ( LOOKUP(ITEM_CODE,I3)  ==  38 ) THEN

! DEPENDS ON: rcf_readflds
                  CALL Rcf_ReadFlds                                 &
                      (AncF_UnitNo,2,                               &
                      NLOOKUP(ipos_27)+LOOKUP_STEP(ipos_27)*STEP,   &
                      LOOKUP(1,LOOKUP_START(FILEANCIL(ipos_27))),   &
                      LEN1_LOOKUP,ICE_EXTENT,                       &
                      FIELD_SIZE,FIXHD(1,FILEANCIL(ipos_27)),       &
                      ICODE,CMESSAGE)

                  IF (ICODE /= 0) THEN
                    ICODE=FIELD+100
                    CMESSAGE='replanca_rcf_replanca :I/O ERROR '
                    GO TO 9999
                  END IF
                  IF ( RLOOKUP(BMDI,I3-1)  /=  RMDI ) THEN
                    ICODE = 700 + FIELD
                    CMESSAGE = 'replanca_rcf_replanca: nonstandard lookup'&
                        //' RMDI in ancil file sea-ice chge times'
                    GO TO 9999
                  END IF

                ELSE
                  ICODE=FIELD+100
                  CMESSAGE='replanca_rcf_replanca :ICE FIELD DATA MISSING'
                  GO TO 9999
                END IF
              END IF
            END IF

!L 3.3 If time interpolation required, read second record

            IF (LINTERPOLATE) THEN

              I1=I2+ LOOKUP_STEP(irec)
              IF (I1 <= FIXHD(152,FILE)) THEN

                IF ( lookup(item_code,lookup_start(file)+i1-1) /= &
                     lookup(item_code,lookup_start(file)+i2-1) ) THEN
                  IF (l_error_ancil_struct) THEN
                    icode=field+100
                    cmessage='replanca_rcf_replanca: start and end fields' // &
                             ' are different.'
                    GO TO 9999
                  END IF
                END IF

! AMIP II and ice depth don't read in ice depth field
                IF (.NOT.(LAMIPII.AND.LICE_DEPTH)) THEN

! DEPENDS ON: rcf_readflds
                  CALL Rcf_ReadFlds                                 &
                      (AncF_UnitNo,1,I1,                            &
                      LOOKUP(1,LOOKUP_START(FILE)),                 &
                      LEN1_LOOKUP,ANCIL2,FIELD_SIZE,FIXHD(1,FILE),  &
                      ICODE,CMESSAGE)
                  
                END IF

                IF (ICODE /= 0) THEN
                  ICODE=FIELD+100
                  CMESSAGE='replanca_rcf_replanca :I/O ERROR '
                  GO TO 9999
                END IF

              ELSE !end of data on file

!L  If end of data has been reached go back to the start.If data is
!L  periodic.
!L  Otherwise cancel time interpolation

                IF (PERIODIC) THEN

                  I1 = NLOOKUP(irec) + LEVEL - 1
                  IF ( lookup(item_code,lookup_start(file)+i1-1) /= &
                       lookup(item_code,lookup_start(file)+i2-1) ) THEN
                    IF (l_error_ancil_struct) THEN
                      icode=field+100
                      cmessage='replanca_rcf_replanca: start and end fields'// &
                               ' are different.'
                      GO TO 9999
                    END IF
                  END IF

! DEPENDS ON: rcf_readflds
                  CALL Rcf_ReadFlds                               &
                      (AncF_UnitNo,1,I1,                          &
                      LOOKUP(1,LOOKUP_START(FILE)),               &
                      LEN1_LOOKUP,ANCIL2,FIELD_SIZE,FIXHD(1,FILE),&
                      ICODE,CMESSAGE)

                  IF (ICODE /= 0) THEN
                    ICODE=FIELD+100
                    CMESSAGE='replanca_rcf_replanca :I/O ERROR '
                    GO TO 9999
                  END IF
                ELSE
                  WRITE(6,'(A,A)') 'REPLANCA: Reached end of ancillary ',   &
                                   'switched off time interpolation.'
                  LINTERPOLATE=.FALSE.
                END IF
              END IF ! End of position on file test
              
              ICODE=0
            END IF ! End LINTERPOLATE

!L 3.4 Perform time interpolation
            
            IF (LINTERPOLATE) THEN
              
              ZERO=0.0

!L Select appropriate time interpolation for each field
!  Snowdepth: set equal to zero if no snow cover

              IF (LSNOW_DEPTH) THEN
                DO I=1,P_FIELD
                  PRES_VALUE(I)=ZERO
                END DO

! For the call to T_INT_C, need to know BMDI is OK for SNOW_CHANGE
!  which was read in from position I2+1.
                IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I2) /= RMDI ) THEN
                  ICODE = 700 + FIELD
                  CMESSAGE = 'replanca_rcf_replanca: nonstandard lookup '&
                      //'RMDI in ancil file snow chge times'
                  GO TO 9999
                END IF

                CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,      &
                    TIME,FIELD_SIZE,SNOW_CHANGE,ANCIL1,PRES_VALUE)

! Ice fraction: ice depth set equal to zero if no ice

              ELSE IF (FIELD == 27.OR.FIELD == 29) THEN
                IF (FIELD == 27) THEN
! For the call to T_INT_C, need to know BMDI is OK for ICE_EXTENT(1,2)
!  which was read in from position I1+1
                  IF (.NOT.LAMIPII) THEN
                    IF (RLOOKUP(BMDI,LOOKUP_START(FILE)+I1) /= RMDI) THEN
                      ICODE = 700 + FIELD
                      CMESSAGE = 'replanca_rcf_replanca: nonstandard lookup'&
                          //' RMDI in ancil file sea-ice chge times'
                      GO TO 9999
                    END IF
                  END IF

                  IF (LAMIPII) THEN
! linear uncontrolled time interpolation
                    CALL T_INT (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA, &
                        TIME,FIELD_SIZE)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.
                    
                    DO I=1,FIELD_SIZE
                      IF (ANCIL_DATA(I) <  0.3) ANCIL_DATA(I)=0.0
                      IF (ANCIL_DATA(I) >  1.0) ANCIL_DATA(I)=1.0
                    END DO

                  ELSE       ! non AMIPII option
                    DO I=1,FIELD_SIZE
                      PRES_VALUE(I)=0
                    END DO

                    CALL T_INT_C (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA, &
                        TIME,FIELD_SIZE,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)

                  END IF     ! end AMIPII test
                
                ELSE IF (FIELD == 29) THEN

                  DO I=1,FIELD_SIZE
                    PRES_VALUE(I)=0
                  END DO

                  CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,       &
                      TIME,FIELD_SIZE,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)

                END IF

! Sea surface temperature, set equal to TFS if ice present

              ELSE IF (FIELD == 28.AND.LT_INT_C) THEN
                IF (LAMIPII) THEN

                  CALL T_INT (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,         &
                      TIME,FIELD_SIZE)
! remove any T below TFS
                  DO I=1,FIELD_SIZE
                    IF (ANCIL_DATA(i) <  TFS)  ANCIL_DATA(I)=TFS
                  END DO

                ELSE     ! non AMIPII option

                  DO I=1,FIELD_SIZE
                    PRES_VALUE(I)=TFS

! Set no_ice_extent indicator for controlled SST interpolation
                    IF (ICE_EXTENT(I,1) == 0) THEN
                      NO_ICE_EXTENT(I)=1.0
                    ELSE
                      NO_ICE_EXTENT(I)=0.0
                    END IF
                  END DO

                  CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,       &
                      TIME,FIELD_SIZE,ICE_EXTENT(1,2),NO_ICE_EXTENT,PRES_VALUE)

                END IF   ! end AMIPII test
! Otherwise linear interpolation in time, unless missing data indicator
! present at either time.

              ELSE

! Time interpolation checks the data against the standard missing data
!   indicator - check that the field is labelled as using the same one.
!  (It is to have the right I1 here that I3 is used above.)
                IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1) /= RMDI .OR.    &
                    RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1) /= RMDI ) THEN
                  WRITE (6,'(A,2F12.1)') 'LOOKUPS:',                       &
                      RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1),                &
                      RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1)
                  ICODE = 700 + FIELD
                  CMESSAGE = 'replanca_rcf_replanca: Nonstandard MDI in '&
                      //'lookup of ancil file'
                  GO TO 9999
                END IF

                LEN=FIELD_SIZE
!L  Ozone, test for zonal mean or full field
                IF (FIELD == 7) THEN
                  IF (LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
                    LEN=P_ROWS
                  END IF
!   Tropopause-based ozone, test for zonal mean or full field.
!   Currently the same test as for conventional ozone.
                ELSE IF (FIELD == 110) THEN
                  IF (LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
                    LEN=P_ROWS
                  END IF
!   Cariolle ozone, test for zonal mean or full field. 
!   Currently same test as for conventional ozone.
                ELSE IF (FIELD >= 178.AND.FIELD <= 185) THEN
                  IF (LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
                    LEN=P_ROWS
                  END IF
                END IF

                CALL T_INT(ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,            &
                    TIME,LEN)

              END IF ! End Lsnow_depth
            
! If no interpolation, copy data into final array

            ELSE ! no interpolation
              IF (LICE_FRACTION) THEN
                IF (LAMIPII) THEN
                  DO I=1,FIELD_SIZE
                    
                    ANCIL_DATA(I)=ICE_EXTENT(I,1)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

                    IF (ANCIL_DATA(I) <  0.3) ANCIL_DATA(I)=0.0
                    IF (ANCIL_DATA(I) >  1.0) ANCIL_DATA(I)=1.0
                    
                  END DO
                ELSE           ! non AMIP II option
                  DO I=1,FIELD_SIZE
                    ANCIL_DATA(I)=ICE_EXTENT(I,1)
                  END DO
                END IF           ! end of AMIPII test
              ELSE IF (LAMIPII.AND.FIELD == 28) THEN
                DO I=1,FIELD_SIZE
                  ANCIL_DATA(I)=ANCIL1(I)
                  IF (ANCIL_DATA(I) <  TFS) ANCIL_DATA(I)=TFS
                END DO
              ELSE
                DO I=1,FIELD_SIZE
                  ANCIL_DATA(I)=ANCIL1(I)
                END DO
                
              END IF
            END IF !End interpolate/no interpolate

!L 3.5 Updating action for each field at each level
!L     Fields replaced except that Sea Surface Temperature may be
!L     incremented. Take appropriate action for each field.

            IF(FIELD <= 2.OR.FIELD == 7.OR.FIELD == 39.OR.FIELD == 40      &
                .OR.FIELD == 41.OR.FIELD == 42.OR.FIELD == 43              &
                .OR.FIELD == 44.OR.FIELD == 45                             &
                                         ! multi-level murk
                .OR.(FIELD >= 48 .AND. FIELD <= 67 .AND. L_UKCA)           &
                                         ! single-level user ancillaries
                .OR.(FIELD >= 68.AND.FIELD <= 70)                          &
                                         !NH3,soot aerosol emissions
                .OR.(FIELD >= 72.AND.FIELD <= 77)                          &
                                         !Sulphur cycle
                .OR.FIELD == 78                                            &
                                         !CO2 EMISSIONS
                .OR.FIELD == 82                                            &
                                         !HADCM2 sulphate aerosol
                .OR.(FIELD >= 90.AND.FIELD <= 109)                         &
                                         !multi-level user ancillaries
                .OR.(FIELD >= 112.AND.FIELD <= 120)                        &
                                         !mineral dust fields
                .OR.(FIELD >= 121.AND.FIELD <= 122)                        &
                                         !Biomass emissions
                .OR.FIELD == 123                                           &
                                         !Seawater DMS concentration
                .OR.(FIELD >= 124.AND.FIELD <=126)                         &
                                         !River routing
                .OR.(FIELD >= 157.AND.FIELD <= 177)                        &
                                         !Aerosol climatologies
                .OR.(FIELD >= 178.AND.FIELD <= 185)                        &
                                         !Cariolle ozone ancillaries
                .OR.(FIELD >= 186.AND.FIELD <= 187)                        &
                                         !OCFF emissions
                .OR.FIELD == 188                                           &
                                         !Unfiltered orography
                ) THEN

!L 3.5.0 Updates at all points

              LEN=FIELD_SIZE
!L  Ozone, test for zonal mean or full field
              IF (FIELD == 7) THEN
                IF (LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
                  LEN=P_ROWS
                END IF
!   Tropopause-based ozone, test for zonal mean or full field.
!   Currently the same test as for conventional ozone.
              ELSE IF (FIELD == 110) THEN
                IF (LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
                  LEN=P_ROWS
                END IF
!   Cariolle ozone, test for zonal mean or full field. 
!   Currently same test as for conventional ozone.
              ELSE IF (FIELD >= 178.AND.FIELD <= 185) THEN
                IF (LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
                  LEN=P_ROWS
                END IF
              END IF

              DO I=1,LEN
                D1(ANCIL_ADD(irec)+I-1+(LEVEL-bot_level)*LEN)=ANCIL_DATA(I)
              END DO

          ! Copy theta level 1 to theta level 0
              IF (level == 1 .AND. bot_level == 0) THEN
                Do I=1,LEN
                  D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
                End Do
              END IF

!L 3.5.1 Updates over all land points

            ELSE IF((field > 2.AND.field < 7)                               &
                .OR.(field > 7.AND.field < 13)                              &
                .OR.(field > 13.AND.field < 27)                             &
                .OR.(field == 32).OR. (field >= 35.AND.field <= 36)         &
                .OR.(field >= 48 .AND. field <= 67 .AND. .NOT. l_ukca)      &
                                      ! single level user ancillaries
                .OR.(field >= 46 .AND. field <= 47)                         &
                                      ! Orographic roughness
                .OR.(field >= 155.AND. field <= 156)                        &
                                      ! Orographic X & Y gradients
                .OR.(field >= 79 .AND. field <= 81)                         &
                                      ! MOSES-I
                .OR.(field >= 83 .AND. field <= 89)                         &
                                      ! MOSES-II
                .OR.(field >= 33 .AND. field <= 34)                         &
                                      ! LSH Topographic index fields
                .OR.(field == 111)                                          &
                                      !COASTAL TILING
                .OR.(field == 193)                                          &
                                      !Snow depth on tiles.
                .OR.(field >= 194.AND.field <= 196)                         &
                                      ! Land surface albedos obs/clim 
               ) THEN
!  Set default value of Z0 over sea
              IF (FIELD == 26) THEN
                DO I=1,P_FIELD
                  IF (SEA(I)) THEN
                    D1(ANCIL_ADD(irec)+I-1+(LEVEL-1)*P_FIELD)=10.E-4
                  END IF
                END DO
              END IF

              DO I=1,P_FIELD
                IF (LAND(I)) THEN
                  D1(ANCIL_ADD(irec)+I-1+(LEVEL-1)*P_FIELD)=ANCIL_DATA(I)
                END IF
              END DO

!L  Reset TSTAR to TM if snow cover present

              IF (LSNOW_DEPTH) THEN
                DO I=1,P_FIELD
                  IF (LAND(I).AND.ANCIL_DATA(I) >  0.0) THEN
                    IF (TSTAR_LAND(I) >  TM) TSTAR_LAND(I)=TM
                  END IF
                END DO
              END IF

! Iceberg calving for the OASIS coupler:
            ELSE IF (field == 13) THEN
              DO i=1,u_field
                d1(ancil_add(irec)+i-1)=ancil_data(i)
              END DO

!L 3.5.2 Ice fraction

            ELSE IF (FIELD == 27) THEN
              
              DO I=1,P_FIELD
                ICE_FRACTION(I)=0.0
                IF (SEA(I)) THEN
                  ICE_FRACTION(I)=ANCIL_DATA(I)
                END IF
              END DO

!L Reduce TSTAR to TFS where ice fraction greater than zero
! Required at present because radiation and boundary layer codes
! assume T* is TFS and ignore any value set in TSTAR.

              IF (.NOT.LTLEADS) THEN
                DO I=1,P_FIELD
                  IF (ICE_FRACTION(I) >  0.0) THEN
                    TSTAR_SSI(I)=MIN(TSTAR_SSI(I),TFS)
                  END IF
                END DO
              END IF

!L 3.5.3 Sea surface temperatures for atmosphere, allow fields to be
!L       incremented rather than replaced

            ELSE IF (FIELD == 28) THEN

              IF (L_SSTANOM) THEN
                DO I=1,P_FIELD
                  TSTAR_ANOM(I)=0.0
                END DO
              END IF

              DO I=1,P_FIELD
!           Calculate SST anomalies over all open sea points,
!           but ignore anomalies over grid points with sea-ice during
!           the forecast stage when updating SSTs
                IF (.NOT. LAND(I)) THEN
                  IF (L_SSTANOM) THEN
                    TSTAR_ANOM(I)=TSTAR_SSI(I)-ANCIL_DATA(I)
                  ELSE
                    TSTAR_SEA(I)=ANCIL_DATA(I)
                    IF (ICE_FRACTION(I) <= 0.0) TSTAR_SSI(I)=TSTAR_SEA(I)
                  END IF
                END IF
              END DO

!L 3.5.4 Sea ice thickness
!L       Update over all sea points (all sea ice points are the only
!L       ones strictly required, but this cannot be determined easily)

            ELSE IF (FIELD == 29) THEN

              DO I=1,P_FIELD
                IF (SEA(I)) THEN
                  D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
                END IF
              END DO

!L 3.5.5 Surface currents

            ELSE IF (FIELD == 30) THEN
              DO I=1,U_FIELD
                D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
              END DO

            ELSE IF (FIELD == 31) THEN
              DO I=1,V_FIELD
                D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
              END DO

            ELSE

              WRITE(6,'(A,I6,A)')' replanca_rcf_replanca: ERROR - FIELD ', &
                  FIELD,' omitted from update block'

            END IF !End tests on FIELD numbers

!L End loop over levels

            I2=I2+1

          END IF
        END DO

!L End loop over ancillary fields (atmosphere)
      END IF ! LAMIPII and ice depth test

    END IF    ! End UPDATE(irec) test     level 1 IF

  END DO

  IF(L_CTILE)THEN
    DO I=1,P_FIELD
      IF(SEA(I).AND.ICE_FRACTION(I) >  0.0)THEN

        IF(LTLEADS.OR.LAMIPII)THEN
          TSTAR_SICE(I)=AMIN1(TSTAR_SICE(I),TFS)
          TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                &
            +(1.-ICE_FRACTION(I))*TSTAR_SEA(I)
        ELSE
          TSTAR_SEA(I)=TFS
          TSTAR_SICE(I)=(TSTAR_SSI(I)                               &
            -(1.-ICE_FRACTION(I))*TSTAR_SEA(I))/ICE_FRACTION(I)
        END IF

      END IF
!
      TSTAR(I)=FLANDG(I)*TSTAR_LAND(I)                              &
        +(1.-FLANDG(I))*TSTAR_SSI(I)
    END DO
  ELSE
    DO I=1,P_FIELD
      IF(LAND(I))THEN
        TSTAR(I)=TSTAR_LAND(I)
      ELSE
        TSTAR(I)=TSTAR_SSI(I)
      END IF
    END DO
  END IF

! Set up surface temperatures:
  IF(L_CTILE)THEN
    DO I=1,P_FIELD
      TSTAR_LAND_CTILE(I)=TSTAR_LAND(I)
      TSTAR_SEA_CTILE(I)=TSTAR_SEA(I)
      ! Ensure consistency with equivalent code in 
      ! replanca_rcf_replanca-rpanca1a.F90. Also helps to avoid 
      ! crazy values of TSTAR_SICE in reconfigured dumps.
      ! TSTAR_SICE_CTILE(I)=TSTAR_SICE(I)
    END DO
  END IF

! Deallocate the temporary storage arrays used for ancil data.
  DEALLOCATE ( ANCIL1 )
  DEALLOCATE ( ANCIL2 )
  DEALLOCATE ( ANCIL_DATA )

! Only close file if anc_file_open is greater than 0.
  IF (anc_file_open > 0) THEN
    CALL file_close (AncF_UnitNo,                             &
                     anc_file(anc_file_open) % anc_env_var,   &
                     len_anc_env_var, ioNameInEnv, ioNoDelete, icode)
    anc_file_open = 0
  END IF

9999  CONTINUE
  RETURN
  END SUBROUTINE replanca_rcf_replanca
