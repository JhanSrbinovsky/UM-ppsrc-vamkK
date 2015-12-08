! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read in Headers from Ancillary Files
!
!  Subroutine inancila_rcf_inancila  - Read in Headers from Ancillary Files
!
! Description:
!   Read in the headers & lookuptables from ancillary files.
!
! Method:
!    For ancillary files that are required, the files are opened
!    and the headers and look-up tables read in.
!
! --------------------------------------------------------------------
! This RCF_ specific routine is required until the UM version of 
! inancila_rcf_inancila is updated to use ANCILmaster files.
! --------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

      SUBROUTINE inancila_rcf_inancila(LEN_FIXHD,LEN_INTHD,LEN_REALHD,               &
     &                    LEN1_LEVDEPC,LEN2_LEVDEPC,                    &
     &                    FIXHD,INTHD,REALHD,LOOKUP,                    &
     &                    A_REALHD,A_LEVDEPC,                           &
     &                    NLOOKUPS,                                     &
     &                    LOOKUP_START,LEN1_LOOKUP,ROW_LENGTH,          &
     &                    loc_row_length,                               &
     &                    P_ROWS,loc_p_rows,                            &
     &                    u_rows,                                       &
     &                    r_row_length,r_rows,                          &
     &                    loc_r_row_length,loc_r_rows,                  &
     &                    p_levels,                                     &
     &                    TR_LEVELS,ST_LEVELS,SM_LEVELS,                &
     &                    OZONE_LEVELS, tpps_ozone_levels,              &
     &                    Ancil_Add,                                    &
     &                    ICODE,CMESSAGE)

      Use Rcf_Model_Mod, Only :                                         &
     &    ZonAvOzone, ZonAvTppsOzone

      Use Rcf_Stashcodes_Mod, Only :                                    &
          stashcode_prog_sec,                                           &
          stashcode_soil_moist

      Use Rcf_Data_Source_Mod, Only :                                   &
          Ancillary_File

      Use Rcf_Items_Mod, Only :                                         &
          Source_Array,                                                 &
          Sctn_Array,                                                   & 
          Item_Array,                                                   &
          Num_Items

      Use Ancil_mod, Only :                                             &
     &    anc_record,      anc_file,        ancRecs,                    &
     &    AncF_UnitNo,     ancFiles,                                    &
     &    levels,          stashancil,                                  &
     &    nlookup,         lookup_step

      Use Locate_Anc_mod, Only :                                        &
     &    Locate_Anc_Field,                                             &
     &    Locate_Anc_File

      USE UM_ParVars, Only :                                            &
     &    mype

      USE PrintStatus_mod, Only :                                       &
     &    PrintStatus,                                                  &
     &    PrStatus_Diag

      USE nlsizes_namelist_mod, ONLY :                                  &
          ntiles

      USE um_input_control_mod, ONLY :                                  &
     &    lcal360

      USE ancilcta_namelist_mod, ONLY :                                 &
          l_sstanom

      Use Ereport_Mod, Only :                                           &
     &    Ereport

      Use Rcf_HeadAddress_Mod, Only :                                   &
     &    SoilDepths,                                                   &
     &    FH_Dataset,                                                   &
     &    FH_Dataset_Ancil

      Use Rcf_FortranIO_Mod, Only :                                     &
     &    Max_Filename_Len

      USE nstypes 
      USE clmchfcg_scenario_mod, ONLY: nsulpat
      USE IO
      USE lookup_addresses

      
      Use Rcf_Exppx_Mod, Only :      &
          Rcf_Exppx
      
      Use Rcf_Ppx_Info_Mod, Only :    &
          STM_Record_Type

      Use Rcf_Grid_Type_Mod, Only : &
          Output_Grid

      Use Rcf_Level_Code_Mod, Only : &
          Rcf_Level_Code

      USE science_fixes_mod, ONLY : &
          l_error_ancil_struct

      IMPLICIT NONE

! Arguments
      Integer :: Len_FixHd     ! Length of fixed header   ) in
      Integer :: Len_IntHd     ! Length of Integer header ) ancillary
      Integer :: Len_RealHd    ! Length of Real header    ) files
      Integer :: Len1_LevDepC  ! ) First and second dimension of model
      Integer :: Len2_LevDepC  ! ) level dependent constants array
      Integer :: Len1_Lookup   ! First dimension for lookup tables.

      Integer :: NLookups       ! Total no of lookup entries in
                                ! ancillary files to be read in.

!     Arrays for ancillary file headers
      Integer, dimension(Len_FixHd,ancFiles)  :: FixHd  ! Fixed H
      Integer, dimension(Len_IntHd,ancFiles)  :: IntHd  ! Integer H
      Real   , dimension(Len_RealHd,ancFiles) :: RealHd ! Real H
      Integer, dimension(Len1_Lookup,NLookups):: Lookup ! Lookup Table

!     Arrays for dump headers
      Real,  dimension(Len_RealHd) :: A_RealHd   ! Real Header
      Real,  dimension(Len1_LevDepC,Len2_LevDepC) :: A_LevDepC
                                ! Level dependent constants array

      Integer :: Row_Length     ! Global row length
      Integer :: Loc_Row_Length ! Local  row length

      Integer :: P_Rows         ! Global no of p rows
      Integer :: Loc_P_Rows     ! Local  no of p rows
      Integer :: U_Rows         ! Global no of u rows
      Integer :: R_Rows         ! Global no of r rows
      Integer :: R_Row_Length   ! Global length of r rows
      Integer :: loc_R_Rows     ! Local no of r rows
      Integer :: loc_R_Row_Length ! Local length of r rows
      Integer :: P_Levels       ! No of model levels
      Integer :: TR_Levels      ! No of Tracer levels
      Integer :: ST_Levels      ! No of Soil Temperature levels
      Integer :: SM_Levels      ! No of Soil Moisture levels
      Integer :: Ozone_Levels   ! No of Ozone levels
      Integer ::tpps_ozone_levels
!                  No of ozone levels in tropopause-based ozone dataset

      Integer, dimension(ancFiles) :: Lookup_Start
                                ! Pointer to first lookup entry for
                                ! each ancillary file
      Integer, dimension(ancRecs) :: Ancil_Add
                                ! Addresses in work space for anc data

      Integer  :: Icode         ! Return Code

      Character (Len=80) :: CMessage  !  Error Message if ICode > 0

! Local variables
      Integer :: I,J,J1,K      ! Loop indices
      Integer :: Irec          ! Loop over record
      Integer :: Section
      Integer :: Item
      Integer :: Lookups
      Integer :: Start_Block
      Integer :: Ozone_Row_Length
      Integer :: tpps_ozone_row_length !! row length for tpps_ozone
      Integer, Parameter :: Dummy = 1
      Integer :: anc_file_no
      Integer :: len_anc_env_var
      Integer :: ipos_27
      Integer :: ipos_9

      Real,    dimension ( (P_levels+1) * 4 ) :: LevDepC
      REAL :: ColDepC(ROW_LENGTH+1),RowDepC(P_ROWS+1)

      Logical :: l_vert_mismatch
      Logical :: Check_Fail    ! flag for more complex checks.
      Logical :: l_soilm_ancil = .FALSE.  ! soil moisture from ancilary

      Character (Len=*), Parameter :: RoutineName = 'inancila_rcf_inancila'
      Character (Len=Max_Filename_Len) ::  AncFileName  ! Anc file name

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
! STASHmaster entry
      Type (STM_Record_Type), Pointer :: STM_Rec
      INTEGER :: bot_level
      INTEGER :: um_levels

      REAL P1,P2
      LOGICAL LNER
      LNER(P1,P2) = ((ABS(P1-P2))  >   (1.E-6*ABS(P1+P2)))

! CL Internal Structure

      ICODE=0
      CMESSAGE=' '
! C
! CL  1.  Initialisation for atmosphere model

      Do I=1,ancRecs
        STASHANCIL(I) = 1000 * ( anc_record(I) % section_number ) +  &
                        anc_record(I) % item_number        
      End Do

!     Set no of levels for ancillary fields.

      DO irec=1,ancrecs

        i = anc_record(irec) % ancil_ref_number

        SELECT CASE (I)
          CASE (7)                         ! Ozone
            levels (irec) = ozone_levels
          CASE (10)                        ! Deep Soil Temps
            levels (irec) = st_levels
          CASE (36)                        ! Multi-layer hydrology
            levels (irec) = sm_levels
          CASE (41, 42, 43)                ! Mulit-layer aerosols
            levels (irec) = tr_levels
          CASE (44, 45, 72 : 76, 90 : 109, 157 : 177, 178 : 185)
            levels (irec) = p_levels       
          CASE (82)
            levels (irec) = nsulpat
          CASE (83)
            levels (irec) = ntype
          CASE (84, 85)
            levels (irec) = npft
          CASE (110)                       ! TppsOzone
            levels(irec) = tpps_ozone_levels           
          CASE (193)
            levels(irec) = ntiles
          CASE DEFAULT                     ! Single level
            levels (irec) = 1
        END SELECT

      END DO

! CL 1.4 Read headers

      LOOKUPS=0

      Do I=1,ancFiles

! C  Initialise LOOKUP_START (=0 implies file I not required)
        LOOKUP_START(I)=0

! Open Ancillary file if required and read in headers & lookup table

      If ( anc_file(I) % anc_file_open  ==  1 ) then

      If (mype == 0) then

       Write(6,*) ' '
       Write(6,*) 'Ancillary File ',anc_file(I) % anc_file_number,' : ',&
     &             anc_file(I) % anc_file_title

       If (PrintStatus >= PrStatus_Diag) Then  ! get file name & print
        Call Fort_Get_Env( anc_file(i) % anc_env_var, 8, AncFileName,   &
     &                     Max_Filename_Len, icode)
        Write (6,*) 'File Name : ',AncFileName(1:len_trim(AncFileName))

       End If

      End If

! Read headers for physical files required

!      Get len of env var
       len_anc_env_var = len_trim ( anc_file(i) % anc_env_var )

!      Open the ancillary file

       call File_Open (AncF_UnitNo, anc_file(i) % anc_env_var,          &
     &                 len_anc_env_var,0,0,icode)

       If (icode /= 0) then
         write (Cmessage, '(A, I4)') 'Problem opening Ancillary file ', &
                                     anc_file(I) % anc_file_number
         Call Ereport ( RoutineName, icode, Cmessage )
       End If

! 1.4.1 Buffer in fixed length header record


        call Setpos (AncF_UnitNo, 0, icode)

        If (icode /= 0) then
          write (Cmessage, '(A, I4)') 'Problem in SETPOS for Ancillary file ', &
                                      anc_file(I) % anc_file_number
          Call Ereport ( RoutineName, icode, Cmessage )
        End If! Check icode

!       Read in fixed header to get array dimensions
! DEPENDS ON: read_flh
        CALL READ_FLH(AncF_UnitNo,FIXHD(1,I),LEN_FIXHD,ICODE,CMESSAGE)
        If (ICODE >  0) THEN
          WRITE (6,*) ' Error in reading fixed header for file ',       &
                      anc_file(I) % anc_file_number
          Go To 9999   !  Return
        End If

! C       Check for negative dimensions
        If (FIXHD(101,I) <= 0) FIXHD(101,I)=1
        If (FIXHD(106,I) <= 0) FIXHD(106,I)=1
        If (FIXHD(111,I) <= 0) FIXHD(111,I)=1
        If (FIXHD(112,I) <= 0) FIXHD(112,I)=1
        If (FIXHD(151,I) <= 0) FIXHD(151,I)=1
        If (FIXHD(152,I) <= 0) FIXHD(152,I)=1
        If (FIXHD(161,I) <= 0) FIXHD(161,I)=1

! Check for valid Ancil file.
        If (FIXHD(FH_Dataset,I) /= FH_Dataset_Ancil) Then
          Write(CMessage,'(A, I4)') 'Invalid fixed header for ancillary file ',&
                                    anc_file(I) % anc_file_number
          
          If (mype == 0) Then
            Call Fort_Get_Env( anc_file(i) % anc_env_var, 8, AncFileName,   &
     &                         Max_Filename_Len, icode)
            Write(6,*) CMessage
            Write(6,*) 'Filename : ', AncFileName(1:Len_Trim(AncFileName))
          End If
          
          icode = 10
          Call Ereport ( RoutineName, icode, Cmessage )
        End If

! C Set start position of boundary fields for file
        LOOKUP_START(I)=LOOKUPS+1

        If (LOOKUPS+FIXHD(152,I) >  NLOOKUPS) THEN
          WRITE (6,*) 'No room in LOOKUP table for Ancillary File ',    &
                      anc_file(I) % anc_file_number
          CMESSAGE='inancila_rcf_inancila: Insufficient space for LOOKUP headers'
          ICODE=14
          Go To 9999   !  Return
        End If


        CALL SETPOS(AncF_UnitNo, 0, ICODE)

        If (icode /= 0) then
          write (Cmessage, '(A, I4, A)')                         &                
     &      'Problem in SETPOS for Ancillary file ',             &
     &      anc_file(I) % anc_file_number,' before READHEAD.'
          Call Ereport ( RoutineName, icode, Cmessage )
        End If! Check icode


! DEPENDS ON: readhead
      If (fixhd(116,I) > 1)THEN
        CALL READHEAD(AncF_UnitNo,                                      &
     &                FIXHD(1,I),LEN_FIXHD,                             &
     &                INTHD(1,I),FIXHD(101,I),                          &
     &                REALHD(1,I),FIXHD(106,I),                         &
     &                LEVDEPC,FIXHD(111,I),FIXHD(112,I),                &
     &                ROWDEPC,FIXHD(116,I),FIXHD(117,I),                &
     &                COLDEPC,FIXHD(121,I),FIXHD(122,I),                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                LOOKUP(1,LOOKUPS+1),FIXHD(151,I),FIXHD(152,I),    &
     &                FIXHD(161,I),                                     &
     &                START_BLOCK,ICODE,CMESSAGE)
      Else
        CALL READHEAD(AncF_UnitNo,                                      &
     &                FIXHD(1,I),LEN_FIXHD,                             &
     &                INTHD(1,I),FIXHD(101,I),                          &
     &                REALHD(1,I),FIXHD(106,I),                         &
     &                LEVDEPC,FIXHD(111,I),FIXHD(112,I),                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                LOOKUP(1,LOOKUPS+1),FIXHD(151,I),FIXHD(152,I),    &
     &                FIXHD(161,I),                                     &
     &                START_BLOCK,ICODE,CMESSAGE)
       End If
       
        If (ICODE >  0) THEN
           WRITE(6,*) 'ERROR in READHEAD for Ancillary File ',          &
                      anc_file(I) % anc_file_number
           Go To 9999   !   Return
        End If

!      Close the ancillary file

       call File_Close (AncF_UnitNo, anc_file(i) % anc_env_var,         &
     &                 len_anc_env_var,0,0,icode)

       If (icode /= 0) then
         write (Cmessage, '(A, I4)') 'Problem closing Ancillary file ', &
                                     anc_file(I) % anc_file_number
         Call Ereport ( RoutineName, icode, Cmessage )
       End If

!   Check calendar indicator is correct if this is set in the ancillary
       IF (fixhd(8,i) /=  imdi) THEN
         IF ((     lcal360 .AND. fixhd(8,i) /= 2) .OR.                  &
              (.NOT.lcal360 .AND. fixhd(8,i) /= 1) ) THEN
           icode = 100 + anc_file(i) % anc_file_number
           cmessage = 'inancila_rcf_inancila : '                        &
                // 'Wrong calendar set in Ancillary File'
           IF (mype == 0) THEN
             WRITE (6,'(a)')                                            &
                  '  ******** Error in inancila_rcf_inancila ********'
             WRITE (6,'(a,i4)')                                         &
                  '  Wrong calendar setting in Ancillary File ',        &
                  anc_file(i) % anc_file_number
             IF (lcal360) THEN
               WRITE (6,'(a)')                                          &
                    '  Model run is set up for 360 day calendar.'
               WRITE (6,'(a)')                                          &
                    '  Ancillary File is for 365 day calendar.'
             ELSE
               WRITE (6,'(a)')                                          &
                    '  Model run is set up for 365 day calendar.'
               WRITE (6,'(a)')                                          &
                    '  Ancillary File is for 360 day calendar.'
             END IF
             WRITE (6,'(a)') '  Rerun with correct ancillary file.'
           END IF
           GO TO 9999   !  Return
         END IF
       ELSE
         IF (mype == 0) THEN
           WRITE (6,'(a)') ' Unspecified calendar type in ancillary file'
         END IF
       END IF
        


! CL 1.4.2 Buffer in integer constants

           If (FIXHD(100,I) >  0) THEN

! Check validity of integer constants array

            If (INTHD(6,I) /= ROW_LENGTH) Then

              ! set fail to true unless reset by special cases.
              Check_Fail=.True.

              Select Case (anc_file(I) % anc_file_number)

              Case (1)         ! Ozone file.
                ! Ozone can have rows of length 1 too. (Zonal Ozone)
                If ( INTHD(6,I) == 1 ) Then
                  Check_Fail=.False.
                End If
              
              Case (30,31)     ! River Routing files
                ! River routing has it's own grid and it's own row len.
                If (INTHD(6,I) == r_row_length) Then
                  Check_Fail=.False.
                End If

              Case (46)         ! Ozone tracer file (Cariolle)
                ! Ozone can have rows of length 1 too. (Zonal Ozone)
                If ( INTHD(6,I) == 1 ) Then
                  Check_Fail=.False.
                End If
             
              End Select

              If (Check_Fail) then
                ICODE=4
                CMESSAGE='inancila_rcf_inancila:integer header error - row length'
                WRITE(6,*) ' INTHD(6) : ',INTHD(6,I),' ?'
                Goto 9999
              End If
            End If
              
            If(INTHD(7,I) /= P_ROWS) Then
              Check_Fail=.True.
              Select Case (anc_file(I) % anc_file_number)

              Case (9)     ! L.S.M
                If (INTHD(6,I) == u_rows) Then
                  Check_Fail=.False.
                End If

              Case (30,31)     ! River Routing files
                ! River Routing has own grid with own no. of rows.
                If (INTHD(7,I) == r_rows) Then
                  Check_Fail=.False.
                End If
              
              End Select

              If (Check_Fail) then
                ICODE=5
                CMESSAGE='inancila_rcf_inancila:integer header error - no. of rows'
                WRITE(6,*) ' INTHD(7) : ',INTHD(7,I),' ?'
                Goto 9999
              End If
            End If

          End If

! CL 1.4.3 Buffer in real constants

          If(FIXHD(105,I) > 0 .AND.                                     &
     &         FIXHD(117,I) < 1 .AND. FIXHD(122,I) < 1) THEN 

! C Check validity of real header and print out information

! Only check if ancillary and model is on C grid (with P at poles)
            IF (fixhd(9,I) == 3 .AND. output_grid % grid_stagger == 3) THEN
              Do J=1,6
                If (LNER(REALHD(J,I),A_REALHD(J))) Then
                If(anc_file(I) % anc_file_number /= 1                      &
                   .OR.(J /= 1.AND.J /= 4))THEN
                  write (6,*) ' Inconsistency in Real Headers.'
                  write (6,*) ' Real Header Values.'
                  Do K=1,6
                  WRITE(6,*) K,' Anc File ',REALHD(K,I),                   &
     &                         ' Model Dump ',A_REALHD(K)
                  End Do
                  ICODE=8
                  CMESSAGE='inancila_rcf_inancila: REAL header Error.'
                  Goto 9999
                End If
                End If
              End Do
            END IF

         End If

! CL 1.4.4 Buffer in level dependent constants if required
! C        Not retained in model after initial check

         If(FIXHD(110,I) >  0) THEN

! CL Only files 1 (Ozone), 25 (tpps_ozone) and
! CL 3 (Soil temperature)should contain multi-level data.
! CL File 2 (Soil moisture,snow depth,fractional snow time
! CL and soil moisture in layers) may possibly also have multi level
! CL data.
! CL FILES 13,14,16 (aerosols, murkiness, user ancil.) may also have
! CL  multi level data.

           If ( anc_file(I) % anc_file_number == 1   .or.               &
                                  !  Ozone File
     &          anc_file(I) % anc_file_number == 14  .or.               &
                                  !  Murkiness File
     &          anc_file(I) % anc_file_number == 16  .or.               &
                                  !  User Ancillary file
     &          anc_file(I) % anc_file_number == 46) Then     
                                  !  cariolle ozone file

! Check that ancillary file is set up for correct vertical levels

            If (fixhd(111,I)-1 /= p_levels) Then
              icode=110
              write (CMESSAGE, '(A, I4, A, I4)')                          &
     &        'Ancillary File set up for wrong no of model levels. Anc ', &
     &        fixhd(111,I)-1, ' Model ',p_levels
              Call Ereport ( RoutineName, icode, Cmessage )
            End If

            l_vert_mismatch = .false.

! Check eta_theta and eta_rho

            Do j=1,p_levels+1
              If (LNER( LEVDEPC(J), A_LEVDEPC(J,1) )) Then
                l_vert_mismatch = .true.
                exit
              End If
            End Do

            Do j=1,p_levels
              If (LNER( LEVDEPC(FIXHD(111,I)+J), A_LEVDEPC(J,2) )) Then
                l_vert_mismatch = .true.
                exit
              End If
            End Do

! Abort if there is a mis-match

            If (l_vert_mismatch) then
              write (6,*) 'Mismatch in vertical levels between model ', &
     &                    'and Ancillary File.'
              write (6,*) 'Anc File : ',anc_file(I) % anc_file_title
              write (6,*) 'Eta_Theta - Model'
              write (6,'(5F10.7)') (a_levdepc(k,1),k=1,p_levels+1)
              write (6,*) 'Eta_Theta - Anc File'
              write (6,'(5F10.7)') (levdepc(k),k=1,p_levels+1)
              write (6,*) 'Eta_Rho   - Model'
              write (6,'(5F10.7)') (a_levdepc(k,2),k=1,p_levels)
              write (6,*) 'Eta_Rho   - Anc File'
              write (6,'(5F10.7)') (levdepc(p_levels+1+k),k=1,p_levels)
                   ICODE=11
              Write (CMESSAGE, '(2A)') 'Mismatch in LEVDEPC ',          &
     &        'between model and Ancillary File.'
              Call Ereport ( RoutineName, icode, Cmessage )
            End If

           !! tropopause-based ozone
           Else If (anc_file(I) % anc_file_number  ==  25) then 
             !! no checks to run

           !   Soil Moisture File
           Else If (anc_file(I) % anc_file_number == 2) Then    
! Check Soil Moisture levels

          If (PrintStatus >= PrStatus_Diag .and. mype == 0 )then
            write (6,*)
            write (6,*) 'SoilDepths = ',SoilDepths
            write (6,*) 'SM_Levels  = ',sm_levels
            do j=1,sm_levels
            write (6,*) 'model ',A_LEVDEPC(J,SoilDepths),               &
     &                   ' anc ',LEVDEPC(fixhd(111,I)*3+J)
            End Do
          End If

! The smow ancillary file can be used for snow but not soil moisture
! so only check soil moisture if it is required.
            DO k=1,num_items
              IF (     (  sctn_array(k)==stashcode_prog_sec  )          &
                  .AND.(  item_array(k)==stashcode_soil_moist)          &
                  .AND.(source_array(k)==Ancillary_File      ))then
                l_soilm_ancil=.TRUE.
              END IF
            END DO

                Do J=1,SM_LEVELS
                  If (LNER(LEVDEPC(fixhd(111,I)*3+J),                   &
     &                     A_LEVDEPC(J,SoilDepths))                     &
     &                     .AND.l_soilm_ancil) THEN
                    ICODE=12
                    CMESSAGE='inancila_rcf_inancila: error in LEVDEPC.'
                   Goto 9999
                  End If
                End Do

          !   Deep Soil Temperature File
          Else If (anc_file(I) % anc_file_number == 3) Then    
          If (PrintStatus >= PrStatus_Diag .and. mype == 0)then
            write (6,*)
            write (6,*) 'SoilDepths = ',SoilDepths
            write (6,*) 'st_levels  = ',st_levels
            do j=1,st_levels
            write (6,*) 'model ',A_LEVDEPC(J,SoilDepths),               &
     &                   ' anc ',LEVDEPC(fixhd(111,I)*3+J)
            End Do
          End If

               Do J=1,ST_LEVELS
                 If (LNER(LEVDEPC(fixhd(111,I)*3+J),                    &
     &                    A_LEVDEPC(J,SoilDepths)))THEN
                   ICODE=12
                   CMESSAGE='inancila_rcf_inancila: error in LEVDEPC.'
                   Goto 9999
                 End If
               End Do

! CL If aerosol file, check against model levels

           !  Aerosol Tracer File
           Else If (anc_file(I) % anc_file_number == 13) THEN    
             Do J=1,TR_LEVELS
               Do J1=1,4
                 If(LNER(LEVDEPC(J+(J1-1)*FIXHD(111,I)),A_LEVDEPC       &
     &                   (J,J1))) THEN
                   WRITE(6,*)'Error in level dependent constants'
                   WRITE(6,*)'Level=',J,' Position=',J1
                   WRITE(6,*)'Value in model =',A_LEVDEPC(J,J1)
                   WRITE(6,*)'Value in ancillary data =',LEVDEPC(J+     &
     &                             (J1-1)*FIXHD(111,I))
                   ICODE=16
                   CMESSAGE='inancila_rcf_inancila: error in LEVDEPC.'
                   Goto 9999
                 End If
               End Do
             End Do

           End If  !  If I

         End If  !  If Fixhd(110,I) > 0

         !  Ozone file or tpps_ozone
         If (anc_file(I) % anc_file_number == 1                         &
             .or. anc_file(I) % anc_file_number == 25 ) THEN   
            WRITE (6,*) ' '
            Ozone_row_length = lookup (lbnpt,lookups+1)
            If (Ozone_row_length == 1) Then
              If (mype == 0) then
                WRITE (6,*) 'OZONE file contains zonal mean data.'
              End If
            Else If (Ozone_row_length == ROW_LENGTH) Then
              If (mype == 0) then
                WRITE (6,*) 'OZONE file contains full fields.'
              End If
            End If

! Check that correct ozone file has been provided.

            If (ZonAvOzone .or. ZonAvTppsOzone) THEN
              If (Ozone_row_length /= 1)THEN
                WRITE (6,*)                                             &
     &          'Error : Zonal Ozone Data is expected for ',P_ROWS,     &
     &          ' rows'
                ICODE = 51
                write (Cmessage, '(A)') 'Wrong Ozone data provided.'
                Call Ereport ( RoutineName, icode, Cmessage )
              End If
            Else
              If (Ozone_Row_Length /= Row_Length)THEN
                WRITE (6,*) 'Error : Ozone Data is expected for ',      &
     &          ROW_LENGTH,' points x ',P_ROWS,' rows.'
                ICODE = 52
                write (Cmessage, '(A)') 'Wrong Ozone data provided.'
                Call Ereport ( RoutineName, icode, Cmessage )
              End If
            End If

          End If  ! If I == 1

!        Add on no of lookup entries
         LOOKUPS=LOOKUPS+FIXHD(152,I)

       End If    !  If Anc file to be opened and fields read in

      End Do    ! Loop over ancillary files (I)

! CL 1.5 Set positions in main data blocks

      ITEM=1
      Do I=1,ancrecs
        If(anc_record (I) % anc_field_read  ==  1)THEN
          Ancil_add(I)=ITEM

! Find whether we need to handle the extra level due to ancillary not supporting
! theta level 0.
          STM_Rec => Rcf_Exppx( 1, stashancil(i)/1000, &
                                MOD(stashancil(i),1000) )
          
          IF ( STM_Rec % lb_code > 0) THEN
            CALL Rcf_Level_Code( STM_Rec % lb_code, bot_level, Output_Grid )
          ELSE
            bot_level = 1
          END IF

! Force bottom level to be 1 (if not 0)
          IF (bot_level /= 0) THEN
            bot_level = 1
          END IF
! um_levels is the number of levels in the UM (not in the ancillary)
          um_levels = levels(i)-bot_level+1

          Select Case ( STASHANCIL(I) )
            Case (60)        ! Ozone
              If( ZonAvOzone )then
                ITEM = ITEM + loc_p_rows * um_levels
              Else
                ITEM = ITEM + loc_row_length * loc_p_rows * um_levels
              End If

            Case (341)       ! tropopause-based ozone
              If( ZonAvTppsOzone )then
                ITEM = ITEM + loc_p_rows * um_levels
              Else
                ITEM = ITEM + loc_row_length * loc_p_rows * um_levels
              End If

            Case (151,152,153) ! River Routing
              ITEM = ITEM + loc_R_Row_Length * loc_R_Rows * um_levels

            Case default ! All 'normal' ancil fields.
              ITEM = ITEM + loc_row_length * loc_p_rows * um_levels

          End Select
        End If
      End Do


! CL 1.6 Set positions of data
      Do I=1,ancRecs

        call locate_anc_file(anc_record (i) % anc_file_number, anc_file_no)

        NLOOKUP(I) =0
        LOOKUP_STEP(I)=0

! C If LOOKUP_START=0 for anc_file_no, no fields required.
        If (LOOKUP_START(anc_file_no) >  0) THEN

          If (  PrintStatus > PrStatus_Diag .and. mype == 0 ) Then
            write (6,*) ' lookup_start > 0 for Anc Ref No ',            &
     &      anc_record(i) % ancil_ref_number
          End If

          Do J=LOOKUP_START(anc_file_no),LOOKUPS

            If (LOOKUP(ITEM_CODE,J) == STASHANCIL(I)) THEN
              NLOOKUP(I)=J-LOOKUP_START(anc_file_no)+1
              Exit
            End If

          End Do

! C Find second occurence of data to set LOOKUP_STEP

          LOOKUP_STEP(I)=0

          If(J <  LOOKUPS) THEN

            Do J1=J+LEVELS(I),LOOKUPS
              If (LOOKUP(ITEM_CODE,J1) == STASHANCIL(I)) THEN
                LOOKUP_STEP(I)=                                         &
     &          J1-NLOOKUP(I)-LOOKUP_START(anc_file_no)+1
                Exit
              End If
            End Do

          End If
          ! Check Ancillary files are consistent since we assume:
          ! for month in month_list:
          !   for STASH in STASH_list:
          !     for level in level_psuedolevel_list:
          ! First check if its periodic since we can only check these at
          ! this point and whether it has multiple fields defined.
          IF (fixhd(10,anc_file_no) == 2 .AND. lookup_step(i) > 0) THEN
            ! Loop over the fields which should only differ by time.
            DO j = lookup_start(anc_file_no) + nlookup(i)-1,               & 
                   lookup_start(anc_file_no) + fixhd(152,anc_file_no) - 1, &
                   lookup_step(i)
              ! If the STASH is different we have a problem.
              IF (lookup(item_code,j) /= stashancil(i) ) THEN
                IF (l_error_ancil_struct) THEN
                  icode = 53
                ELSE
                  icode = -53
                END IF
                WRITE (cmessage,'(A,I4)')                                  &
                  'Incorrect structure for ancillary file ',               &
                  anc_record (i) % anc_file_number
                CALL ereport ( routinename, icode, cmessage )
                EXIT
              END IF
            END DO
          END IF
        END IF

      End Do

! CL SET LEVELS=2 FOR ICE FRACTION AND SNOW DEPTH, TO INDICATE PRESCENCE
! CL fractional time fields

      call locate_anc_field (9 , ipos_9 )
      call locate_anc_field (27 , ipos_27 )
      
      LEVELS(ipos_9 )=2      !   Snow Depth
      LEVELS(ipos_27)=2      !   Ice Fraction

      If (  PrintStatus >= PrStatus_Diag .and. mype == 0 ) Then

      write (6,*) ' '
      write (6,*) ' Summary from inancila_rcf_inancila '
      write (6,*) ' LOOKUP_START '
      write (6,*)   LOOKUP_START
      write (6,*) ' Ancil Ref Numbers '
      write (6,'(10I5)') anc_record(:) % ancil_ref_number
      write (6,*) ' Anc_field_read '
      write (6,'(10I5)') anc_record(:) % anc_field_read
      write (6,*) ' LOOKUP_STEP '
      write (6,'(10I5)') LOOKUP_STEP
      write (6,*) ' NLOOKUP '
      write (6,'(10I5)') NLOOKUP
      write (6,*) ' LEVELS '
      write (6,'(10I5)') LEVELS
      write (6,*) ' ANCIL_ADD '
      write (6,'(10I7)') ANCIL_ADD
      write (6,*) ' '

      End If

 9999 CONTINUE
      RETURN
      END SUBROUTINE inancila_rcf_inancila
