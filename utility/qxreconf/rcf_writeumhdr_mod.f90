! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
!   Simple wrapper routine for WritHead
!
! Method:
!   Wrapper routine for WritHead - setting the addressing correctly
!   and converting from the GRID data-type.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

MODULE Rcf_WriteUMhdr_mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Rcf_WriteUMhdr( UMhdr )

    USE UM_Types
    USE writhead_mod
    USE Rcf_UMhead_Mod, ONLY :                                                 &
        UM_header_type,                                                        &
        LenFixHd
    USE Ereport_Mod, ONLY :                                                    &
        Ereport
    USE UM_ParVars, ONLY :                                                     &
        mype
    USE io
    IMPLICIT NONE

    ! Arguments
    TYPE (UM_header_type), INTENT( INOUT ) :: UMhdr

    ! Local constant
    CHARACTER (LEN=*), PARAMETER :: RoutineName = 'Rcf_WriteUMhdr'

    ! Local variables
    CHARACTER (LEN=80)           :: Cmessage
    INTEGER                      :: ErrorStatus
    INTEGER   ::  WordAddress        ! Position on file, used in SETPOS
    INTEGER   ::  number_of_data_words_in_memory
    INTEGER   ::  number_of_data_words_on_disk
    INTEGER   ::  disk_address       ! Current rounded disk address
                                     ! and final data length
    INTEGER   :: start_block_dummy   ! Writhead calculates start of data but not
                                     ! on sector boundary - use dummy instead.

    !----------------------------------------------------------------

    WordAddress = 0
    ErrorStatus = 0

    CALL SetPos (UMhdr % UnitNum, WordAddress, ErrorStatus)

    IF ( ErrorStatus /= 0) THEN
      Cmessage = 'SETPOS failure, repositioning to start of dump!'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

    ! DEPENDS ON: set_dumpfile_address
    CALL Set_Dumpfile_Address( UMhdr % FixHd, LenFixHd, UMhdr % Lookup,        &
        UMhdr % Len1Lookup,  UMhdr % Len2Lookup,                               &
        number_of_data_words_in_memory,                                        &
        number_of_data_words_on_disk, disk_address)

    CALL WritHead(  UMhdr % UnitNum,                                           &
        UMhdr % FixHd,     LenFixHd,                                           &
        UMhdr % IntC,      UMhdr % LenIntC,                                    &
        UMhdr % RealC,     UMhdr % LenRealC,                                   &
        UMhdr % LevDepC,   UMhdr % Len1LevDepC,                                &
        UMhdr % Len2LevDepC,                                                   &
        UMhdr % RowDepC,   UMhdr % Len1RowDepC,                                &
        UMhdr % Len2RowDepC,                                                   &
        UMhdr % ColDepC,   UMhdr % Len1ColDepC,                                &
        UMhdr % Len2ColDepC,                                                   &
        UMhdr % FldsOfC,   UMhdr % Len1FldsOfC,                                &
        UMhdr % Len2FldsOfC,                                                   &
        UMhdr % ExtraC,    UMhdr % LenExtraC,                                  &
        UMhdr % HistFile,  UMhdr % LenHistFile,                                &
        UMhdr % CompFldI1, UMhdr % LenCompFldI1,                               &
        UMhdr % CompFldI2, UMhdr % LenCompFldI2,                               &
        UMhdr % CompFldI3, UMhdr % LenCompFldI3,                               &
        UMhdr % Lookup,    UMhdr % Len1Lookup,                                 &
        UMhdr % Len2Lookup,                                                    &
        UMhdr % LenData,                                                       &
        umFortranIntegerSize()*8,                                              &
        start_block_dummy,                                                     &
        ErrorStatus,        CMessage )

    IF (ErrorStatus /= 0) THEN
      Cmessage='Error in Rcf_WritHead!'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

    RETURN

  END SUBROUTINE Rcf_WriteUMhdr

END MODULE Rcf_WriteUMhdr_Mod
