! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

! Description:
!  This routine writes a number of fields out to dump.
!
! Method:
!  Setpos is used to set the position in the dump to write too.
!  The addressing is calculated and a number of variables
!  (grid code and lbc_levels) required by lower level routines are set.
!  Data is written by rcf_write_multi


!-------------------------------------------------------------------
! Note that this routine does not have a module as it is called
! will different type arguents for D1
!-------------------------------------------------------------------

SUBROUTINE Rcf_WritFlds( nftout,   number_of_fields,                           &
    position, lookup,len1_lookup,                                              &
    d1,   len_buf,  fixhd,                                                     &
    errorstatus,  cmessage, multi_pe, stash_in)

  USE UM_ParVars, ONLY :                                                       &
      mype,                                                                    &
      blsizeu,                                                                 &
      blsizep,             blsizev

  USE PrintStatus_mod, ONLY :                                                  &
      PrintStatus,                                                             &
      PrStatus_Diag

  USE Rcf_Exppx_Mod, ONLY :                                                    &
      Rcf_Exppx

  USE Rcf_Ppx_Info_Mod, ONLY :                                                 &
      STM_record_type


  USE Rcf_Write_Multi_Mod, ONLY :                                              &
      Rcf_Write_Multi

  USE Ereport_Mod, ONLY :                                                      &
      Ereport

  USE Rcf_Level_Code_Mod, ONLY :                                               &
      Rcf_Level_Code

  USE Rcf_Grid_Type_Mod, ONLY :                                                &
      Output_Grid

  USE Rcf_Global_To_Local_Mod, ONLY :                                          &
      Rcf_Get_Fld_Type

  USE io
  USE lookup_addresses

  USE cppxref_mod, ONLY :                                                      &
       ppx_atm_ozone,                                                          &
       ppx_atm_tzonal
  IMPLICIT NONE

! Arguments
  INTEGER, INTENT(IN)  :: nftout            ! Unit number for I/O
  INTEGER, INTENT(IN)  :: number_of_fields ! No of fields to be
                                          ! written
  INTEGER, INTENT(IN)  :: len_buf           ! Length of I/O buffer
  INTEGER, INTENT(IN)  :: position          ! Field number from
                                          ! which to begin I/O
  INTEGER, INTENT(IN)  :: fixhd(*)          ! fixed length header
  INTEGER, INTENT(IN)  :: Len1_Lookup       ! 1st dim of lookup
  INTEGER, INTENT(IN)  :: Lookup(Len1_Lookup,*) ! lookup table
  REAL,    INTENT(IN)  :: d1(*)             ! Data to write


  INTEGER, INTENT(OUT)           :: ErrorStatus
  CHARACTER(LEN=80), INTENT(OUT) :: Cmessage

  Logical, Intent(In)            :: Multi_PE
  Logical, Intent(In)            :: stash_in ! Use input STASHmaster

! Local variables
  INTEGER                :: k            ! index
  INTEGER                :: D1_Off       ! Offset in D1 array
  INTEGER                :: len_io
  INTEGER                :: Field_Start  ! word address to begin I/O
  INTEGER                :: Data_Write_Size ! data to write
  INTEGER                :: Fld_Type
  INTEGER                :: local_len    ! local size of field
  INTEGER                :: grid_type
  INTEGER                :: LBC_levels
  INTEGER                :: field_model
  INTEGER                :: field_sect
  INTEGER                :: field_item
  CHARACTER (LEN=*), PARAMETER :: RoutineName='Writeflds'
  
  REAL                   :: A_IO
  TYPE (STM_RECORD_TYPE) :: Stash_Record

  errorstatus = 0
  cmessage    = ' '

!L 2. Buffer out NUMBER_OF_FIELDS fields of data:
  D1_Off = 0
  DO k = position, position + number_of_fields-1

! Location on disk from which to begin I/O
    Field_Start=lookup(lbegin,k)

    Data_Write_Size = Lookup( lbnrec, k)

! Position file pointer

    CALL setpos( nftout, Field_Start, ErrorStatus)


! Get some information about this field
    field_item  = MOD(lookup(42,k),1000)
    field_sect  = (lookup(42,k)-field_item)/1000
    field_model = lookup(45,k)
  
    stash_record = rcf_exppx( field_model, field_sect, field_item, &
                              stash_in_arg = stash_in )
    grid_type = stash_record % grid_type

! For atmosphere zonal ozone fields - set to zonal grid type
    IF ( grid_type == ppx_atm_ozone .AND. lookup(lbnpt,k) == 1) THEN
      Stash_Record % grid_type = ppx_atm_tzonal
    END IF

    CALL Rcf_Write_Multi( nftout, d1( D1_Off + 1), Data_Write_Size,            &
        len_io, local_len, a_io,                                               &
        lookup(1,k), fixhd(12),                                                &
        Stash_Record, Multi_PE )

! Reset the STASHmaster grid-type to what it really should be
    Stash_Record % grid_type = grid_type

! Check for I/O errors
    IF (a_io /= -1.0 .OR. len_io /= Data_Write_Size ) THEN
      WRITE(6,'('' *ERROR* Writing field no'',I5)')k
      IF (fixhd(5) < 6 .OR. fixhd(5) > 10) THEN ! Not AC/Cx/Cov/ObSt
! DEPENDS ON: pr_look
        CALL pr_look(lookup,lookup,len1_lookup,k)
      END IF

! DEPENDS ON: ioerror
      CALL ioerror('rcf_writfields: buffer out of real data', &
          a_io, len_io,Data_Write_Size)
      ErrorStatus = NINT(a_io)+1
      cmessage    = 'Rcf_WRITFLDS:I/O error'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

! Data summary used to be here - removed for time being

! Increment offset by size of level written out
    D1_Off = D1_Off + Local_Len

  END DO

  RETURN
END SUBROUTINE Rcf_WritFlds


