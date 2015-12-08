! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Top level routine for reading, handling and writeing GRIB data.

MODULE Rcf_Grib_Control_Mod

IMPLICIT NONE

! SUBROUTINE Rcf_Grib_Control
!
! Description: This is the top level routine for reading in GRIB data
!              and generating a UM style intermediary dump containing
!              that data
!
! Method: read in the GRIB records one at a time storing the information
!         in a set of dynamically allocated lists. Using this info, set
!         up the corresponding UM style headers.
!         re-read records from the GRIB file, this time performing
!         simple transformations (if necessary) on the data before
!         writing it out to the intermediary UM dump.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
CONTAINS
  SUBROUTINE Rcf_Grib_Control()

! uses variables and routines from other modules

    USE UM_ParVars, ONLY :                                                     &
        mype

    USE Rcf_GRIB_Block_Params_Mod

    USE Rcf_GRIB_FldSort_Mod, ONLY  :                                          &
        Rcf_GRIB_FldSort

    USE Rcf_StashCodes_Mod

    USE Rcf_HeadAddress_Mod

    USE Rcf_GRIB_Lookups_Mod               ! Contains cross ref table for
                                       ! Stash => ECMWF parameter Id's
                                       ! and appropriate params

    USE Rcf_Grib_Debug_Tools_Mod, ONLY :                                       &
        Grib_Debug_Print_Basics,                                               &
        Grib_Debug_Print_Blocks,                                               &
        Grib_Debug_ListCounts

    USE Rcf_Grib_Assign_Mod, ONLY :                                            &
        Rcf_Grib_Assign

    USE Rcf_Grib_Check_Mod, ONLY :                                             &
        Rcf_Grib_Check

    USE Rcf_Grib_SetHdr_Mod, ONLY :                                            &
        Rcf_Grib_SetHdr

    USE Rcf_Grib_Dest_List_Mod, ONLY :                                         &
        Rcf_Grib_Dest_List

    USE Rcf_UMhead_Mod, ONLY :                                                 &
        LenFixHd,                                                              &
        um_header_type            ! Derived containing UM header info

    USE EReport_Mod, ONLY :                                                    &
        EReport

    USE Rcf_WriteUMhdr_Mod, ONLY :                                             &
        Rcf_WriteUMhdr

    USE io

    USE Rcf_FreeUMhdr_Mod, ONLY :                                              &
        Rcf_FreeUMhdr

    USE Rcf_FortranIO_Mod, ONLY :                                              &
        Rcf_Get_Unit,                                                          &
        RCF_Free_Unit

    USE Rcf_Grib_Read_Data_Mod, ONLY :                                         &
        Rcf_Grib_Read_Data

    USE Rcf_Grib_Spcl_Ctl_Mod, ONLY :                                          &
        Rcf_Grib_Spcl_Ctl

    USE Rcf_Grib_Spcl_Hdr_Mod, ONLY :                                          &
        Rcf_Grib_Spcl_Hdr

    USE Rcf_Grid_Type_Mod, ONLY :                                              &
        Grid_type,                                                             &
        Output_Grid

    USE PrintStatus_mod, ONLY :                                                &
        PrintStatus,                                                           &
        PrStatus_Min,           &       ! =1 Minimum output
        PrStatus_Normal,        &       ! =2 Short informative output
        PrStatus_Oper,          &       ! =3 Full informative output
        PrStatus_Diag                   ! =4 Extra Diagnostic output

    USE lookup_addresses
    IMPLICIT NONE

! Comdecks
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
! contains LBLREC (amongst others)

! Local variables

    TYPE (Um_Header_type)            :: Hdr_Dmy
    TYPE (Um_Header_type)            :: Hdr_Itm
! local var to preserve Output_Grid values during 1st stage
    TYPE (Grid_type)                 :: Storage_Grid

    CHARACTER (LEN=*), PARAMETER     :: RoutineName='Rcf_GRIB_Control'
    CHARACTER (LEN=80)               :: Cmessage(2)   ! used for EReport
    CHARACTER (LEN=20)               :: cFormat
    INTEGER                          :: ErrorStatus   ! used for EReport
    INTEGER                          :: errstat
    REAL                             :: rerrstat
    INTEGER                          :: i
    INTEGER                          :: cnter,Criteria
    INTEGER                          :: dummy1
    INTEGER                          :: dummy2
    INTEGER                          :: disk_address

!-----------------------------------------------------------------------
!  Variables to do specifically with the GRIB record
!-----------------------------------------------------------------------

    TYPE (Grib_Record),POINTER       :: Current       !\ Pointer to current
                                                  !/ grib record

!An array of pointer pairs to form the head and tail of the field lists
    TYPE (List_Marker)               :: Lists(0:grib_max_fields)

    INTEGER, PARAMETER               :: len_max_c = 256
    CHARACTER (LEN=len_max_c)        :: c_Tmp_Buff(1)
    INTEGER                          :: act_io_len
    INTEGER                          :: pos_in_file
    INTEGER                          :: readable_io_len

    LOGICAL                          :: Order
    LOGICAL, PARAMETER               :: Ascending  = .TRUE.
    LOGICAL, PARAMETER               :: Descending = .FALSE.

!=======================================================================
!  Initialise variables
!=======================================================================

    Storage_grid = Output_Grid

    NULLIFY(Current)

! Allocate a dummy header (used to hold values other routines would copy
! from the input dump) and the header for the Intermediate dump
    ALLOCATE (Hdr_Dmy % FixHd(LenFixHd))
    ALLOCATE (Hdr_Itm % FixHd(LenFixHd))

! Nullify the components of Hdr_Dmy that are not being set anywhere
! to avoid later problems with deallocation.
    NULLIFY (Hdr_Dmy % IntC)
    NULLIFY (Hdr_Dmy % CompFldI1)
    NULLIFY (Hdr_Dmy % CompFldI2)
    NULLIFY (Hdr_Dmy % CompFldI3)
    NULLIFY (Hdr_Dmy % Lookup)
    NULLIFY (Hdr_Dmy % RealC)
    NULLIFY (Hdr_Dmy % LevDepC)
    NULLIFY (Hdr_Dmy % RowDepC)
    NULLIFY (Hdr_Dmy % ColDepC)
    NULLIFY (Hdr_Dmy % FldsOfC)
    NULLIFY (Hdr_Dmy % ExtraC)
    NULLIFY (Hdr_Dmy % HistFile)

! Make sure the list pointers are NULL and the counts are zero
    DO i = 0,grib_max_fields
      NULLIFY( Lists(i) % Begin )
      NULLIFY( Lists(i) % End )
      Lists(i) % LstCount = 0
    END DO

    pos_in_file = 3   ! Offset to allow for small rewind in each loop

!=======================================================================
!  Open GRIB File
!=======================================================================

    CALL Rcf_Get_Unit( Hdr_Dmy % UnitNum )


    CALL File_Open( Hdr_Dmy % UnitNum, 'AINITIAL', 8, 0, 0, errstat)

    IF ( errstat /= 0 ) THEN
      Cmessage(1) = 'Failed to Open GRIB file'
      ErrorStatus = 10
      CALL EReport( RoutineName, ErrorStatus, Cmessage(1) )
    ELSE
      IF ( PrintStatus >= PrStatus_Diag  ) THEN
        IF ( mype == 0 ) THEN
          WRITE (6,'(A)') "Opened Grib Data File"
        END IF
      END IF
    END IF

!=======================================================================
!  Loop, reading blocks from file, until no more records found
!=======================================================================

    Reading: DO                      ! Endless loop until exit supplied

      pos_in_file = pos_in_file - 3  ! Small rewind to allow Index to seek
                                 ! for GRIB properly

!=======================================================================
!  Find the begining of a record by seeking 'GRIB'
!=======================================================================

  ! Set the position within the file to 'pos_in_file'.
      CALL SetPos8(Hdr_Dmy % UnitNum, pos_in_file,errstat)
      IF ( errstat /= 0 ) THEN
        WRITE(Cmessage(1),'(A,I2)')   'Failed trying to SetPos to ',           &
            pos_in_file
        ErrorStatus = 20
        CALL EReport( RoutineName, ErrorStatus, Cmessage(1) )
      END IF

  ! Clear the buffer and read in section of file to search for 'GRIB'
      c_Tmp_Buff(1) = " "
  ! We need to read in data (even if its less than we expected) so lets find
  ! smallest value of possible number of words and the max amount of data we
  ! would like to read.
      readable_io_len = MIN(readableWords(Hdr_Dmy % UnitNum, wordLength=1),    &
          len_max_c)
      IF (readable_io_len > 0) THEN
        CALL Buffin (Hdr_Dmy % UnitNum, c_Tmp_Buff(:), 1, readable_io_len,     &
            act_io_len, rerrstat)
      ELSE
    ! The amount of readable data seems to be 0 or less.  Set values to be
    ! consistent on all processors.
        act_io_len = 0
        rerrstat = -1.0
      END IF

      IF ( rerrstat > 0 ) THEN
        Cmessage(1)    = 'Failed trying to read short block from GRIB file'
        ErrorStatus = 30
        CALL EReport( RoutineName, ErrorStatus, Cmessage(1) )
      END IF

      i = INDEX (c_Tmp_Buff(1), "GRIB")   ! Serch the block for text 'GRIB'

!=======================================================================
!        -> _Didn't_ Find 'GRIB' - try again
!=======================================================================

      IF (i == 0) THEN               ! Index didn't find GRIB in the buffer

        IF (PrintStatus >= PrStatus_Diag) THEN
          IF ( mype == 0 ) THEN
            WRITE (6,*) "Buffer block read contained no 'GRIB' indicator"
          END IF
        END IF

    ! If data read is smaller than data requested then must have read
    ! until the end of the file. => Exit Do loop
        IF (len_max_c > readable_io_len) THEN
          IF (PrintStatus >= PrStatus_Diag) THEN
            IF ( mype == 0 ) THEN
              WRITE (6,'(2(A,I3),A)') 'len_max_c =', len_max_c,                &
                  ' readable_io_len =', readable_io_len,                       &
                  'Exiting from Read loop'
            END IF
          END IF

          EXIT Reading             ! This forces exit from reading the loop
        ELSE
      ! Force (if some what slow) progress through the file
          pos_in_file = pos_in_file + act_io_len
        END IF

!=======================================================================
!        -> _Found_ 'GRIB' - read and decode data
!=======================================================================

      ELSE                                   ! GRIB WAS found

!=======================================================================
!  Allocate 'GRIB Record' to store header Info
!=======================================================================

        ALLOCATE(Current)

        NULLIFY(Current % VertCoords)

        Current % Block_1(:) = 0
        Current % Block_2(:) = 0
        Current % Block_3(:) = 0
        Current % Block_4(:) = 0
        Current % Block_R(:) = 0.000

        Current % Start_pos  = 0
        Current % StashCode  = -1
        Current % Data_Type  = Grb_Data_Real
        Current % Num_Fp     = 0
        Current % Num_Vert   = 0
        Current % Num_Bitmap = 0
        Current % Num_Quasi  = 0

!=======================================================================
!  Read binary form of record into CRecord work array
!=======================================================================

        pos_in_file = pos_in_file + i - 1      ! set the pos in the file
        Current % Start_pos = pos_in_file      ! record start pos in header

        CALL Rcf_Grib_Read_Data(Hdr_Dmy % UnitNum,Current,FpData,              &
            LenArrayMax,pos_in_file)

!=======================================================================
!  Test record and assign to correct list
!=======================================================================

        CALL Rcf_Grib_Assign(Current, Lists)

   ! Set position in file to end of record read
        pos_in_file = pos_in_file + Current % Block_0(p_Mes_Len)

!=======================================================================
!  End Loop "until no more records found"
!=======================================================================

      END IF                               ! Test to find 'GRIB' in buffer

    END DO Reading

!=======================================================================
!  Sort the 'Lists' into sensible orders
!=======================================================================

    cFormat = "(3A,I2,A)"
    DO i = 1, grib_max_fields

      IF (ASSOCIATED(Lists(i) % Begin) .AND.                                   &
          Lists(i) % LstCount > 1) THEN    ! Checks list has members

    ! For each case, set the order desired, the criteria (currently
    ! any value in Block 1) on which to perform the ordering, and
    ! the output message text.
        SELECT CASE(i)

      ! Pressure level fields ordered away from surface.
        CASE (grib_U_field,                                                    &
            grib_V_field,                                                      &
            grib_W_field,                                                      &
            grib_Temp_field,                                                   &
            grib_Q_field,                                                      &
            grib_QCL_field,                                                    &
            grib_QCF_field,                                                    &
            grib_CC_field,                                                     &
            grib_ozone_field,                                                  &
            grib_NOX_field,                                                    &
            grib_CH4_field,                                                    &
            grib_CO_field,                                                     &
            grib_HCHO_field,                                                   &
            grib_GO3_field,                                                    &
            grib_NO2_field,                                                    &
            grib_NO_field,                                                     &
            grib_SO2_field,                                                    &
            grib_HNO3_field,                                                   &
            grib_PAN_field,                                                    &
            grib_C2H6_field,                                                   &
            grib_C3H8_field,                                                   &
            grib_OMFRSH_field,                                                 &
            grib_OMAGD_field,                                                  &
            grib_BCFRSH_field,                                                 &
            grib_BCAGD_field,                                                  &
            grib_UMDUST1_field,                                                &
            grib_UMDUST2_field,                                                &
            grib_UMDUST3_field,                                                &
            grib_UMDUST4_field,                                                &
            grib_UMDUST5_field,                                                &
            grib_UMDUST6_field,                                                &
            grib_UMSO4AITK_field,                                              &
            grib_UMSO4ACCU_field)

          Order = Descending
          Criteria = p_Lvl_Desc_1     ! 1st Level description parameter
          WRITE(cMessage(1),cFormat) "Re-ordering ",                           &
              TRIM(Lists(i) % Begin % Desc),                                   &
              ". ", Lists(i) % LstCount,  " entries"

        ! Soil levels ordered down from surface
        CASE (grib_Soil_Temp_field,                                            &
            grib_Soil_Moist_field)

          Order = Ascending
          Criteria = p_Lvl_Desc_1     ! 1st Level description parameter
          WRITE(cMessage(1),cFormat) "Re-ordering ",                           &
              TRIM(Lists(i) % Begin % Desc),                                   &
              ". ", Lists(i) % LstCount,  " entries"

        CASE Default
          Order = Descending
          Criteria = p_Lvl_Desc_1
          WRITE(cMessage(1),'(A,I2,A,I2,A)')                                   &
              "List ", i, ", ", Lists(i) % LstCount,                           &
              " entries re-ordered using defaults. No Specific rule exists"
        END SELECT

        CALL rcf_Grib_FldSort(Lists(i), Criteria, Order)

      ELSE
        WRITE(cMessage(1),'(A,I2,A,I2,A)') "List ", i ,                        &
            " Not re-ordered because it contained ", Lists(i) % LstCount,      &
            " entry"
      END IF

      IF (PrintStatus >= PrStatus_Diag) THEN
        IF ( mype == 0 ) THEN
          WRITE (6,'(A)') cMessage(1)
        END IF
      END IF

    END DO ! I over no. of fields

!=======================================================================
! Loop through lists showing count
!=======================================================================

    IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
      WRITE (6,*)
      CALL Grib_Debug_ListCounts(Lists)
      WRITE (6,*)
    END IF

!=======================================================================
! Loop through lists displaying basic info
!=======================================================================

    IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
      WRITE (6,*)
      CALL Grib_Debug_Print_Basics(Lists)
      WRITE (6,*)
    END IF

!=======================================================================
! Perform Basic Data Integrity Checks
!=======================================================================

    CALL Rcf_Grib_Check(Lists)

!=======================================================================
! Setting up the headers.
!=======================================================================

    CALL Rcf_Grib_SetHdr(Lists,Output_Grid,Hdr_Dmy,Hdr_Itm)

!=======================================================================
!  Open the Output Dump for the Intermediary Dump
!=======================================================================

    CALL Rcf_Get_Unit( Hdr_Itm % UnitNum )

    IF ( PrintStatus >= PrStatus_Diag  ) THEN
      IF ( mype == 0 ) THEN
        WRITE (6,*) "Opening Intermediary Dump"
      END IF
    END IF

    CALL File_Open( Hdr_Itm % UnitNum, 'RECONTMP', 8, 1, 0, errstat)

    IF ( errstat /= 0 ) THEN
      Cmessage(1)    = 'Failed to Open Intermediary Dump file'
      ErrorStatus = 99
      CALL EReport( RoutineName, ErrorStatus, Cmessage(1) )
    END IF

!=======================================================================
!  Set up addressing info used when writing Data
!=======================================================================

! DEPENDS ON: set_dumpfile_address
    CALL Set_Dumpfile_Address( Hdr_Itm % FixHd,                                &
        LenFixHd,                                                              &
        Hdr_Itm % Lookup,                                                      &
        Hdr_Itm % Len1Lookup,                                                  &
        Hdr_Itm % Len2Lookup,                                                  &
        dummy1, dummy2, disk_address)

!=======================================================================
!  Write the Field data
!=======================================================================

    cnter = 0

!Loop across all lists
    DO i = 1, grib_max_fields

      IF (ASSOCIATED(Lists(i) % Begin) ) THEN

        Current => Lists(i) % Begin
        DO While (ASSOCIATED(Current))

          cnter = cnter + 1

!=======================================================================
!  Read binary form of record into CRecord work array
!=======================================================================

      CALL Rcf_Grib_Read_Data(Hdr_Dmy % UnitNum,                 &
                              Current,FpData,                    &
                              lenarraymax,                       &
                              Current % Start_pos)

      ! Call to routine which calls for 'special' handling of given data

          CALL Rcf_Grib_Spcl_Ctl(FpData,LgData,Lists,Current,                  &
              Hdr_Dmy,Hdr_Itm,i,cnter)

      ! Set .FALSE. for multi_pe indicates data is held on a single PE
      ! Since the GRIB reconfigruation is creating a dump which will be read
      ! back in later by the reconfiguration as an input dump we should make
      ! sure we write this dump using the input STASHmaster.  Its more for
      ! consistency than anything else so set stash_in to .TRUE.. This was
      ! needed for Endgame where GRIB reconfiguration assumes New Dynamics but
      ! the STASHmaster's were different - this was fixed when Endgame used the
      ! same STASHmaster as ND with the new level bottom code.
          IF (Current % Data_Type == Grb_Data_Real ) THEN
! DEPENDS ON: rcf_writflds
            CALL Rcf_WritFlds( Hdr_Itm % UnitNum, 1,cnter,                     &
                Hdr_Itm % Lookup, Hdr_Itm % Len1Lookup, FpData,                &
                Hdr_Itm % Lookup( lblrec , cnter ),                            &
                Hdr_Itm % FixHd,                                               &
                ErrorStatus, Cmessage, .FALSE., .TRUE. )

          ELSE IF (Current % Data_Type == Grb_Data_Log ) THEN
! DEPENDS ON: rcf_writflds
            CALL Rcf_WritFlds( Hdr_Itm % UnitNum, 1,cnter,                     &
                Hdr_Itm % Lookup, Hdr_Itm % Len1Lookup, LgData,                &
                Hdr_Itm % Lookup( lblrec , cnter ),                            &
                Hdr_Itm % FixHd,                                               &
                ErrorStatus, Cmessage, .FALSE., .TRUE. )

          ELSE
            Cmessage(1) = 'Failed to write to temp dump :Unknown Data Type'
            ErrorStatus = 99
            CALL EReport( RoutineName, ErrorStatus, Cmessage(1) )
          END IF

          Current => Current % Next

        END DO  ! members of a list

      END IF ! Associated(Lists(I) % Begin)

    END DO  ! (I) loop over all lists

!=======================================================================
!  Write the Header on the Intermediary Dump
!=======================================================================

    CALL Rcf_WriteUMhdr( Hdr_Itm )

!=======================================================================
!  Closing the Grib File
!=======================================================================


    CALL File_Close ( Hdr_Dmy % UnitNum, 'AINITIAL', 8, 0, 0, errstat)
    CALL Rcf_Free_Unit ( Hdr_Dmy % UnitNum )


    CALL File_Close ( Hdr_Itm % UnitNum, 'RECONTMP', 8, 0, 0, errstat)
    CALL Rcf_Free_Unit ( Hdr_Itm % UnitNum )

!=======================================================================
!  Quick tidy up behind ourselves.
!=======================================================================

! Free up space taken by the headers
    CALL Rcf_FreeUMhdr ( Hdr_Dmy )
    CALL Rcf_FreeUMhdr ( Hdr_Itm )

! Deallocate the lists
    DO i = 0, grib_max_fields
      CALL Rcf_Grib_Dest_List (Lists(i))
    END DO

! return Output_Grid values to those read in by namelist
    Output_Grid = Storage_Grid

    RETURN

  END SUBROUTINE Rcf_Grib_Control
END MODULE Rcf_Grib_Control_Mod
