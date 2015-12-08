! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers a field from many processors to one processor

Module Rcf_Gather_Field_Mod

!  Subroutine Rcf_Gather_field_real  - gather real fields
!  Subroutine Rcf_Gather_field_log   - gathers logical fields
!
! Description:
!  Takes a model field that has been decomposed over a group of
!  processors, and gathers the data together so that one processor
!  contains the entire global field.
!
! Method:
!  A send and receive map is constructed which instructs the GCOM
!  permute operation to do a gather from all processors in the
!  group to the GATHER_PE
!--------------------------------------------------------------------
!   ************************************************************
!      This subroutine has an overloaded interface so as to
!      cope cleanly with a number of different data types.
!   ***********************************************************
!--------------------------------------------------------------------
!
! Derived from UM4.5 code.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


  IMPLICIT NONE

Interface Rcf_Gather_Field
  Module Procedure Rcf_Gather_Field_Log, Rcf_Gather_Field_Real
End Interface

Contains

!--------------------------------------------------------------------
!   Thus subroutine copes with the REAL case
!--------------------------------------------------------------------
Subroutine Rcf_Gather_Field_Real(LOCAL_FIELD,    GLOBAL_FIELD, &
                                 LOCAL_ROW_LEN,  LOCAL_ROWS,   &
                                 GLOBAL_ROW_LEN, GLOBAL_ROWS,  &
                                 GATHER_PE,      PROC_GROUP )

Use Ereport_Mod, Only : &
    Ereport

USE UM_ParVars   ! Use a lot of this
USE gcom_mod

IMPLICIT NONE

! Subroutine Arguments:
Integer, Intent(In)    :: Local_Row_Len  ! row length in local field
Integer, Intent(In)    :: Global_Row_Len ! row length in global field
Integer, Intent(In)    :: Local_Rows     ! number of rows locally
Integer, Intent(In)    :: Global_Rows    ! number of rows globally
Integer, Intent(In)    :: Gather_PE      ! PE on which to gather data
Integer, Intent(In)    :: Proc_Group     ! group of PEs involved

Real, Intent(In)       :: Local_Field( Local_Row_Len * Local_Rows )
Real, Intent(Out)      :: Global_Field( Global_Row_Len * Global_Rows)

! Local variables
Integer                :: Info        ! return code from Gcom
Integer                :: Fld_Type    ! P, U or V
Integer                :: Flag        ! dummy for Gcom
Integer                :: iproc       ! processor number
Integer                :: ErrorStatus
Character (Len=80)     :: Cmessage
Character (Len=*), Parameter :: RoutineName = 'Rcf_Gather_Field_Real'

Integer, Save          :: n_mess_to_rec
Integer, Save          :: send_map(7,1)
Integer, ALLOCATABLE, SAVE :: receive_map(:,:)
Integer, Save          :: Old_Global_Row_Len       = -1234  ! old values
Integer, Save          :: Old_Global_Rows          = -1234  ! from
Integer, Save          :: Old_Proc_Group           = -1234  ! previous
Integer, Save          :: Old_Gather_PE            = -1234  ! calls to
Integer, Save          :: Old_Decomp               = -1234  ! routine

!-----------------------------------------------------------------
! 0.0 Can we use the same send/receive map that we calculated
!     last time round?
!-----------------------------------------------------------------

IF (.NOT.ALLOCATED(receive_map)) &
    ALLOCATE(receive_map(7,NPROC_MAX))

IF ((GLOBAL_ROW_LEN .NE. old_GLOBAL_ROW_LEN) .OR. &
    (GLOBAL_ROWS    .NE. old_GLOBAL_ROWS   ) .OR. &
    (PROC_GROUP     .NE. old_PROC_GROUP    ) .OR. &
    (GATHER_PE      .NE. old_GATHER_PE     ) .OR. &
    (current_decomp_type .NE. old_DECOMP   )) THEN
!       Different arguments from the last call so we need
!       to calculate a new send/receive map

! 1.0 Find the type of field (P or U) being done

  IF (GLOBAL_ROWS    .EQ. glsizep(2) .AND. &
      GLOBAL_ROW_LEN .EQ. glsizep(1) ) THEN
    ! This should cover P field and U (global) field for C grid
    fld_type=fld_type_p
  ELSEIF (GLOBAL_ROWS    .EQ. glsizev(2) .AND. &
          GLOBAL_ROW_LEN .EQ. glsizev(1) ) THEN
    fld_type=fld_type_v
  ELSEIF (GLOBAL_ROWS    .EQ. glsizeu(2) .AND. &
          GLOBAL_ROW_LEN .EQ. glsizeu(1) ) THEN
    fld_type=fld_type_u
  ELSEIF (GLOBAL_ROWS    .EQ. glsizer(2) .AND. &
          GLOBAL_ROW_LEN .EQ. glsizer(1) ) THEN
    fld_type=fld_type_r
  ELSE
    Cmessage = 'Unrecognised field type'
    ErrorStatus = 10
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  ENDIF


! 2.0 Set up send map

  send_map(S_DESTINATION_PE,1) = GATHER_PE !  processor to send to

  send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,1) = Offy*LOCAL_ROW_LEN+1+Offx
!       first data to send

  IF (fld_type .EQ. fld_type_p) THEN
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)= blsizep(2)
    send_map(S_STRIDE_IN_SEND_ARRAY,1) = blsizep(1) + 2 * Offx
    send_map(S_ELEMENT_LENGTH,1) = blsizep(1)
    send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
           datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  ELSE IF (fld_type .EQ. fld_type_u) THEN
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)= blsizeu(2)
    send_map(S_STRIDE_IN_SEND_ARRAY,1) = blsizeu(1) + 2 * Offx
    send_map(S_ELEMENT_LENGTH,1) = blsizeu(1)
    send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
           datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  ELSE IF (fld_type .EQ. fld_type_v) THEN
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)= blsizev(2)
    send_map(S_STRIDE_IN_SEND_ARRAY,1) = blsizev(1) + 2 * Offx
    send_map(S_ELEMENT_LENGTH,1) = blsizev(1)
    send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
           datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  ELSE IF (fld_type .EQ. fld_type_r) THEN
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)= blsizer(2)
    send_map(S_STRIDE_IN_SEND_ARRAY,1) = blsizer(1) + 2 * Offx
    send_map(S_ELEMENT_LENGTH,1) = blsizer(1)
    send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
           datastartr(1)+(datastartr(2)-1)*GLOBAL_ROW_LEN
  ENDIF

  send_map(S_STRIDE_IN_RECV_ARRAY,1) = GLOBAL_ROW_LEN
!       stride between rows in global data


! 3.0 Set up the receive map (for PE GATHER_PE only)

! Assume here that this group consists of all processors
! We'll get some new GCG functionality soon to improve this

  n_mess_to_rec=0

  IF (mype .EQ. GATHER_PE) THEN
    DO iproc=0,nproc-1
      receive_map(R_SOURCE_PE,iproc+1) = iproc

      IF (fld_type .EQ. fld_type_p) THEN
        receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizep(1)
        receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                g_blsizep(2,iproc)
        receive_map(R_ELEMENT_LENGTH,iproc+1) = g_blsizep(1,iproc)
        receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) = g_blsizep(1,iproc)
        receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                Offy*g_blsizep(1,iproc)+Offx+1

      ELSEIF (fld_type .EQ. fld_type_u) THEN
        receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizeu(1)
        receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                g_blsizeu(2,iproc)
        receive_map(R_ELEMENT_LENGTH,iproc+1) = g_blsizeu(1,iproc)
        receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) = &
                g_blsizeu(1,iproc)
        receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                Offy*g_blsizeu(1,iproc)+Offx+1

      ELSEIF (fld_type .EQ. fld_type_v) THEN
      receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
              g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizev(1)
        receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                g_blsizev(2,iproc)
        receive_map(R_ELEMENT_LENGTH,iproc+1) = g_blsizev(1,iproc)
        receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) = &
                g_blsizev(1,iproc)
        receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                Offy*g_blsizev(1,iproc)+Offx+1

      ELSEIF (fld_type .EQ. fld_type_r) THEN
      receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                g_datastartr(1,iproc) + (g_datastartr(2,iproc) -1 ) &
                * glsizer(1)
        receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                g_blsizer(2,iproc)
        receive_map(R_ELEMENT_LENGTH,iproc+1) = g_blsizer(1,iproc)
        receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) = g_blsizer(1,iproc)
        receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                Offy*g_blsizer(1,iproc)+Offx+1
      ENDIF

! Should the following go in the above list?
      receive_map(R_STRIDE_IN_RECV_ARRAY,iproc+1) = GLOBAL_ROW_LEN


    ENDDO
    n_mess_to_rec=nproc
  ENDIF

  old_GLOBAL_ROW_LEN=GLOBAL_ROW_LEN
  old_GLOBAL_ROWS=GLOBAL_ROWS
  old_PROC_GROUP=PROC_GROUP
  old_GATHER_PE=GATHER_PE
  old_DECOMP=current_decomp_type

ENDIF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

flag=0  ! This is currently ignored at GCG v1.1
info=GC_NONE

CALL GCG_RALLTOALLE(LOCAL_FIELD,send_map,1, &
                    LOCAL_ROW_LEN*LOCAL_ROWS, &
                    GLOBAL_FIELD,receive_map,n_mess_to_rec, &
                    GLOBAL_ROW_LEN*GLOBAL_ROWS, &
                    PROC_GROUP,flag,info)

Return
End Subroutine Rcf_Gather_Field_Real

!--------------------------------------------------------------------
!   Thus subroutine copes with the LOGICAL case
!--------------------------------------------------------------------
Subroutine Rcf_Gather_Field_Log(LOCAL_FIELD,    GLOBAL_FIELD, &
                                LOCAL_ROW_LEN,  LOCAL_ROWS,   &
                                GLOBAL_ROW_LEN, GLOBAL_ROWS,  &
                                GATHER_PE,      PROC_GROUP )

Use Ereport_Mod, Only : &
    Ereport

USE UM_ParVars   ! Use a lot of this
USE gcom_mod

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
Integer, Intent(In)    :: Local_Row_Len  ! row length in local field
Integer, Intent(In)    :: Global_Row_Len ! row length in global field
Integer, Intent(In)    :: Local_Rows     ! number of rows locally
Integer, Intent(In)    :: Global_Rows    ! number of rows globally
Integer, Intent(In)    :: Gather_PE      ! PE on which to gather data
Integer, Intent(In)    :: Proc_Group     ! group of PEs involved

Logical, Intent(In)    :: Local_Field( Local_Row_Len * Local_Rows )
Logical, Intent(Out)   :: Global_Field( Global_Row_Len * Global_Rows)

! Local variables
Integer                :: Info        ! return code from Gcom
Integer                :: Fld_Type    ! P, U or V
Integer                :: Flag        ! dummy for Gcom
Integer                :: iproc       ! processor number
Integer                :: ErrorStatus
Character (Len=80)     :: Cmessage
Character (Len=*), Parameter :: RoutineName = 'Rcf_Gather_Field_Log'

Integer, Save          :: n_mess_to_rec
Integer, Save          :: send_map(7,1)
Integer, ALLOCATABLE, SAVE :: receive_map(:,:)
Integer, Save          :: Old_Global_Row_Len       = -1234  ! old values
Integer, Save          :: Old_Global_Rows          = -1234  ! from
Integer, Save          :: Old_Proc_Group           = -1234  ! previous
Integer, Save          :: Old_Gather_PE            = -1234  ! calls to
Integer, Save          :: Old_Decomp               = -1234  ! routine

!-----------------------------------------------------------------
! 0.0 Can we use the same send/receive map that we calculated
!     last time round?
!-----------------------------------------------------------------

IF (.NOT.ALLOCATED(receive_map)) &
    ALLOCATE(receive_map(7,NPROC_MAX))

IF ((GLOBAL_ROW_LEN .NE. old_GLOBAL_ROW_LEN) .OR. &
    (GLOBAL_ROWS    .NE. old_GLOBAL_ROWS   ) .OR. &
    (PROC_GROUP     .NE. old_PROC_GROUP    ) .OR. &
    (GATHER_PE      .NE. old_GATHER_PE     ) .OR. &
    (current_decomp_type .NE. old_DECOMP   )) THEN
!       Different arguments from the last call so we need
!       to calculate a new send/receive map

! 1.0 Find the type of field (P or U) being done

  IF (GLOBAL_ROWS    .EQ. glsizep(2) .AND. &
      GLOBAL_ROW_LEN .EQ. glsizep(1) ) THEN
    ! This should cover P field and U (global) field for C grid
    fld_type=fld_type_p
  ELSEIF (GLOBAL_ROWS    .EQ. glsizev(2) .AND. &
          GLOBAL_ROW_LEN .EQ. glsizev(1) ) THEN
    fld_type=fld_type_v
  ELSEIF (GLOBAL_ROWS    .EQ. glsizeu(2) .AND. &
          GLOBAL_ROW_LEN .EQ. glsizeu(1) ) THEN
    fld_type=fld_type_u
  ELSEIF (GLOBAL_ROWS    .EQ. glsizer(2) .AND. &
          GLOBAL_ROW_LEN .EQ. glsizer(1) ) THEN
    fld_type=fld_type_r
  ELSE
    Cmessage = 'Unrecognised field type'
    ErrorStatus = 10
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  ENDIF


! 2.0 Set up send map

  send_map(S_DESTINATION_PE,1) = GATHER_PE !  processor to send to

  send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,1) = Offy*LOCAL_ROW_LEN+1+Offx
!       first data to send

  IF (fld_type .EQ. fld_type_p) THEN
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)= blsizep(2)
    send_map(S_STRIDE_IN_SEND_ARRAY,1) = blsizep(1) + 2 * Offx
    send_map(S_ELEMENT_LENGTH,1) = blsizep(1)
    send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
           datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  ELSE IF (fld_type .EQ. fld_type_u) THEN
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)= blsizeu(2)
    send_map(S_STRIDE_IN_SEND_ARRAY,1) = blsizeu(1) + 2 * Offx
    send_map(S_ELEMENT_LENGTH,1) = blsizeu(1)
    send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
           datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  ELSE IF (fld_type .EQ. fld_type_v) THEN
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)= blsizev(2)
    send_map(S_STRIDE_IN_SEND_ARRAY,1) = blsizev(1) + 2 * Offx
    send_map(S_ELEMENT_LENGTH,1) = blsizev(1)
    send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
           datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  ELSE IF (fld_type .EQ. fld_type_r) THEN
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)= blsizer(2)
    send_map(S_STRIDE_IN_SEND_ARRAY,1) = blsizer(1) + 2 * Offx
    send_map(S_ELEMENT_LENGTH,1) = blsizer(1)
    send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
           datastartr(1)+(datastartr(2)-1)*GLOBAL_ROW_LEN
  ENDIF


  send_map(S_STRIDE_IN_RECV_ARRAY,1) = GLOBAL_ROW_LEN
!       stride between rows in global data


! 3.0 Set up the receive map (for PE GATHER_PE only)

! Assume here that this group consists of all processors
! We'll get some new GCG functionality soon to improve this

  n_mess_to_rec=0

  IF (mype .EQ. GATHER_PE) THEN
    DO iproc=0,nproc-1
      receive_map(R_SOURCE_PE,iproc+1) = iproc

      IF (fld_type .EQ. fld_type_p) THEN
        receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizep(1)
        receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                g_blsizep(2,iproc)
        receive_map(R_ELEMENT_LENGTH,iproc+1) = g_blsizep(1,iproc)
        receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) = g_blsizep(1,iproc)
        receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                Offy*g_blsizep(1,iproc)+Offx+1

      ELSEIF (fld_type .EQ. fld_type_u) THEN
        receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizeu(1)
        receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                g_blsizeu(2,iproc)
        receive_map(R_ELEMENT_LENGTH,iproc+1) = g_blsizeu(1,iproc)
        receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) = &
                g_blsizeu(1,iproc)
        receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                Offy*g_blsizeu(1,iproc)+Offx+1

      ELSEIF (fld_type .EQ. fld_type_v) THEN
      receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
              g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizev(1)
        receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                g_blsizev(2,iproc)
        receive_map(R_ELEMENT_LENGTH,iproc+1) = g_blsizev(1,iproc)
        receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) = &
                g_blsizev(1,iproc)
        receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                Offy*g_blsizev(1,iproc)+Offx+1

      ELSEIF (fld_type .EQ. fld_type_r) THEN
      receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                g_datastartr(1,iproc) + (g_datastartr(2,iproc) -1 ) &
                * glsizer(1)
        receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                g_blsizer(2,iproc)
        receive_map(R_ELEMENT_LENGTH,iproc+1) = g_blsizer(1,iproc)
        receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) = g_blsizer(1,iproc)
        receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                Offy*g_blsizer(1,iproc)+Offx+1
      ENDIF

! Should the following go in the above list?
      receive_map(R_STRIDE_IN_RECV_ARRAY,iproc+1) = GLOBAL_ROW_LEN


    ENDDO
    n_mess_to_rec=nproc
  ENDIF

  old_GLOBAL_ROW_LEN=GLOBAL_ROW_LEN
  old_GLOBAL_ROWS=GLOBAL_ROWS
  old_PROC_GROUP=PROC_GROUP
  old_GATHER_PE=GATHER_PE
  old_DECOMP=current_decomp_type

ENDIF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

flag=0  ! This is currently ignored at GCG v1.1
info=GC_NONE

CALL GCG_RALLTOALLE(LOCAL_FIELD,send_map,1, &
                    LOCAL_ROW_LEN*LOCAL_ROWS, &
                    GLOBAL_FIELD,receive_map,n_mess_to_rec, &
                    GLOBAL_ROW_LEN*GLOBAL_ROWS, &
                    PROC_GROUP,flag,info)

Return
End Subroutine Rcf_Gather_Field_Log

End Module Rcf_Gather_Field_Mod


