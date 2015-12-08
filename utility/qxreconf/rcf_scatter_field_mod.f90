! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Scatters a field from one processor to many processors

Module Rcf_Scatter_Field_Mod

!  Subroutine Rcf_Scatter_Field_Real - for real data
!  Subroutine Rcf_Scatter_field_Log  - for logical data
!
! Description:
!  Takes a model field which is stored entirely on one processor
!  and distributes it over a group of processors.
!
! Method:
!  A send and receive map is constructed which instructs the GCOM
!  permute operation to do a scatter to all processors in the
!  group from the SCATTER_PE
!
!  Derived from UM4.5 code.
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

!--------------------------------------------------------------------
!   ************************************************************
!      This subroutine has an overloaded interface so as to
!      cope cleanly with a number of different data types.
!   ***********************************************************
!--------------------------------------------------------------------
  USE ereport_mod, ONLY : ereport
IMPLICIT NONE

Interface Rcf_Scatter_Field
  Module Procedure Rcf_Scatter_Field_Log, Rcf_Scatter_Field_Real
End Interface

Contains

!-----------------------------------------------------------------
! This deals with real data
!-----------------------------------------------------------------

Subroutine Rcf_Scatter_Field_Real(LOCAL_FIELD,    GLOBAL_FIELD, &
                                  LOCAL_ROW_LEN,  LOCAL_ROWS,   &
                                  GLOBAL_ROW_LEN, GLOBAL_ROWS,  &
                                  SCATTER_PE,     PROC_GROUP    )

USE UM_ParVars   !

USE gcom_mod

Use Ereport_Mod, Only : &
    Ereport

  USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
Integer, Intent(In)   :: Local_Row_Len   ! len of rows in local field
Integer, Intent(In)   :: Local_Rows      ! num of rows in local field
Integer, Intent(In)   :: Global_Row_Len  ! len of rows in global field
Integer, Intent(In)   :: Global_Rows     ! num of rows in global field
Integer, Intent(In)   :: Scatter_PE      ! processor to scatter from
Integer, Intent(In)   :: Proc_Group      ! group ID of involved PEs

Real, Intent(In)      :: Global_Field( global_row_len * global_rows )
Real, Intent(Out)     :: Local_Field ( local_row_len * local_rows )

! Local variables
Integer               :: Info            ! return code from GCOM
Integer               :: fld_type
Integer               :: iproc
Integer               :: flag
Integer               :: ErrorStatus
Character (Len=80)    :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Scatter_Field_Real'

Integer, Save         :: n_mess_to_send
Integer, ALLOCATABLE, SAVE :: send_map(:,:)
Integer, Save         :: receive_map(7,1)
Integer, Save         :: Old_Global_Row_Len   = -1234   ! old values
Integer, Save         :: Old_Global_Rows      = -1234   ! from previous
Integer, Save         :: Old_Proc_Group       = -1234   ! calls to this
Integer, Save         :: Old_Scatter_PE       = -1234   ! routine
Integer, Save         :: Old_Decomp           = -1234   ! ...

IF (.NOT.ALLOCATED(send_map)) &
    ALLOCATE(send_map(7,NPROC_MAX))

!-------------------------------------------------------

! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

IF ((GLOBAL_ROW_LEN .NE. old_GLOBAL_ROW_LEN) .OR. &
    (GLOBAL_ROWS    .NE. old_GLOBAL_ROWS   ) .OR. &
    (PROC_GROUP     .NE. old_PROC_GROUP    ) .OR. &
    (SCATTER_PE     .NE. old_SCATTER_PE    ) .OR. &
    (current_decomp_type .NE. old_DECOMP  )) THEN
!       Different arguments from the last call so we need
!       to calculate a new send/receive map

! 1.0 Find the type of field (P or U) being done

  IF (GLOBAL_ROWS        .EQ. glsizep(2) .AND. &
      GLOBAL_ROW_LEN     .EQ. glsizep(1) ) THEN
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
    write(6,'(A,I4,I4)')'Rows/Row length provided: ',global_rows,global_row_len
    write(6,'(A,I4,I4)')'P - fields:', glsizep(2), glsizep(1)
    write(6,'(A,I4,I4)')'V - fields:', glsizev(2), glsizev(1)
    write(6,'(A,I4,I4)')'U - fields:', glsizeu(2), glsizeu(1)
    write(6,'(A,I4,I4)')'R - fields:', glsizer(2), glsizer(1)
    Cmessage = 'Unable to determine field type from supplied rows/row_length'
    ErrorStatus = 10
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  ENDIF

! 2.0 Set up the send map (for PE SCATTER_PE only)

! Assume here that this group consists of all processors
! We'll get some new GCG functionality soon to improve this

  n_mess_to_send=0

  IF (mype .EQ. SCATTER_PE) THEN
    DO iproc=0,nproc-1
      send_map(S_DESTINATION_PE,iproc+1) = iproc

      IF (fld_type .EQ. fld_type_p) THEN
        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*GLOBAL_ROW_LEN
        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                 g_blsizep(2,iproc)
        send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizep(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizep(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = &
                 g_blsizep(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizep(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizep(1,iproc)

      ELSEIF (fld_type .EQ. fld_type_u) THEN
        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*GLOBAL_ROW_LEN
        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                 g_blsizeu(2,iproc)
        send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizeu(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizeu(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = &
                 g_blsizeu(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizeu(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizeu(1,iproc)

      ELSEIF (fld_type .EQ. fld_type_v) THEN
        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*GLOBAL_ROW_LEN
        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                 g_blsizev(2,iproc)
        send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizev(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizev(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizev(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizev(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizev(1,iproc)

      ELSEIF (fld_type .EQ. fld_type_r) THEN
        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                 g_datastartr(1,iproc)+                  &
                (g_datastartr(2,iproc)-1)*GLOBAL_ROW_LEN
        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                 g_blsizer(2,iproc)
        send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizer(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizer(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = &
                 g_blsizer(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizer(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizer(1,iproc)

      ENDIF

      send_map(S_STRIDE_IN_SEND_ARRAY,iproc+1) = GLOBAL_ROW_LEN
! not sure about above line....

    ENDDO
    n_mess_to_send=nproc
  ENDIF

! 3.0 Set up the receive map

  receive_map(R_SOURCE_PE,1) = SCATTER_PE
  receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
              Offy*LOCAL_ROW_LEN+1+Offx

  IF (fld_type .EQ. fld_type_p) THEN
    receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) =  blsizep(2)
    receive_map(R_STRIDE_IN_RECV_ARRAY,1) = blsizep(1) + 2 * Offx
    receive_map(R_ELEMENT_LENGTH,1) = blsizep(1)
    receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,1) = &
              datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  ELSE IF (fld_type .EQ. fld_type_u) THEN
    receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) = blsizeu(2)
    receive_map(R_STRIDE_IN_RECV_ARRAY,1) = blsizeu(1) + 2 * Offx
    receive_map(R_ELEMENT_LENGTH,1) = blsizeu(1)
    receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,1) = &
              datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  ELSE IF (fld_type .EQ. fld_type_v) THEN
    receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) = blsizev(2)
    receive_map(R_STRIDE_IN_RECV_ARRAY,1) = blsizev(1) + 2 * Offx
    receive_map(R_ELEMENT_LENGTH,1) = blsizev(1)
    receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,1) = &
              datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  ELSE IF (fld_type .EQ. fld_type_r) THEN
    receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) = blsizer(2)
    receive_map(R_STRIDE_IN_RECV_ARRAY,1) = blsizer(1) + 2 * Offx
    receive_map(R_ELEMENT_LENGTH,1) = blsizer(1)
    receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,1) = &
              datastartr(1)+(datastartr(2)-1)*GLOBAL_ROW_LEN
  ENDIF

  receive_map(R_STRIDE_IN_SEND_ARRAY,1) = GLOBAL_ROW_LEN


  old_GLOBAL_ROW_LEN=GLOBAL_ROW_LEN
  old_GLOBAL_ROWS=GLOBAL_ROWS
  old_PROC_GROUP=PROC_GROUP
  old_SCATTER_PE=SCATTER_PE
  old_DECOMP=current_decomp_type

ENDIF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

flag=GC_NONE  ! This is currently ignored at GCG v1.1
info=GC_NONE

CALL GCG_RALLTOALLE(GLOBAL_FIELD, send_map,n_mess_to_send, &
                    GLOBAL_ROW_LEN*GLOBAL_ROWS,            &
                    LOCAL_FIELD,  receive_map,  1,         &
                    LOCAL_ROW_LEN*LOCAL_ROWS,              &
                    PROC_GROUP, flag, info)


RETURN
End Subroutine Rcf_Scatter_Field_Real


!------------------------------------------------------------------
! This deals with logical fields
!------------------------------------------------------------------
Subroutine Rcf_Scatter_Field_Log(LOCAL_FIELD,    GLOBAL_FIELD, &
                                 LOCAL_ROW_LEN,  LOCAL_ROWS,   &
                                 GLOBAL_ROW_LEN, GLOBAL_ROWS,  &
                                 SCATTER_PE,     PROC_GROUP    )

USE UM_ParVars   !
USE gcom_mod

Use Ereport_Mod, Only : &
    Ereport

IMPLICIT NONE

! Subroutine Arguments:
Integer, Intent(In)   :: Local_Row_Len   ! len of rows in local field
Integer, Intent(In)   :: Local_Rows      ! num of rows in local field
Integer, Intent(In)   :: Global_Row_Len  ! len of rows in global field
Integer, Intent(In)   :: Global_Rows     ! num of rows in global field
Integer, Intent(In)   :: Scatter_PE      ! processor to scatter from
Integer, Intent(In)   :: Proc_Group      ! group ID of involved PEs

Logical, Intent(In)      :: Global_Field( global_row_len * global_rows )
Logical, Intent(Out)     :: Local_Field ( local_row_len * local_rows )

! Local variables
Integer               :: Info            ! return code from GCOM
Integer               :: fld_type
Integer               :: iproc
Integer               :: flag
Integer               :: ErrorStatus
Character (Len=80)    :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Scatter_Field_Log'

Integer, Save         :: n_mess_to_send
Integer, ALLOCATABLE, SAVE :: send_map(:,:)
Integer, Save         :: receive_map(7,1)
Integer, Save         :: Old_Global_Row_Len   = -1234   ! old values
Integer, Save         :: Old_Global_Rows      = -1234   ! from previous
Integer, Save         :: Old_Proc_Group       = -1234   ! calls to this
Integer, Save         :: Old_Scatter_PE       = -1234   ! routine
Integer, Save         :: Old_Decomp           = -1234   ! ...

IF (.NOT.ALLOCATED(send_map)) &
    ALLOCATE(send_map(7,NPROC_MAX))

!-------------------------------------------------------

! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

IF ((GLOBAL_ROW_LEN .NE. old_GLOBAL_ROW_LEN) .OR. &
    (GLOBAL_ROWS    .NE. old_GLOBAL_ROWS   ) .OR. &
    (PROC_GROUP     .NE. old_PROC_GROUP    ) .OR. &
    (SCATTER_PE     .NE. old_SCATTER_PE    ) .OR. &
    (current_decomp_type .NE. old_DECOMP  )) THEN
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

! 2.0 Set up the send map (for PE SCATTER_PE only)

! Assume here that this group consists of all processors
! We'll get some new GCG functionality soon to improve this

  n_mess_to_send=0

  IF (mype .EQ. SCATTER_PE) THEN
    DO iproc=0,nproc-1
      send_map(S_DESTINATION_PE,iproc+1) = iproc

      IF (fld_type .EQ. fld_type_p) THEN
        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*GLOBAL_ROW_LEN
        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                 g_blsizep(2,iproc)
        send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizep(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizep(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = &
                 g_blsizep(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizep(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizep(1,iproc)

      ELSEIF (fld_type .EQ. fld_type_u) THEN
        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*GLOBAL_ROW_LEN
        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                 g_blsizeu(2,iproc)
        send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizeu(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizeu(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = &
                 g_blsizeu(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizeu(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizeu(1,iproc)

      ELSEIF (fld_type .EQ. fld_type_v) THEN
        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*GLOBAL_ROW_LEN
        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                 g_blsizev(2,iproc)
        send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizev(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizev(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizev(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizev(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizev(1,iproc)

      ELSEIF (fld_type .EQ. fld_type_r) THEN
        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
                 g_datastartr(1,iproc)+                  &
                (g_datastartr(2,iproc)-1)*GLOBAL_ROW_LEN
        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = &
                 g_blsizer(2,iproc)
        send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizer(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizer(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizer(1,iproc)
        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = &
                 Offy*g_blsizer(1,iproc)+Offx+1
        send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = g_blsizer(1,iproc)

      ENDIF

      send_map(S_STRIDE_IN_SEND_ARRAY,iproc+1) = GLOBAL_ROW_LEN
! not sure about above line....

    ENDDO
    n_mess_to_send=nproc
  ENDIF

! 3.0 Set up the receive map

  receive_map(R_SOURCE_PE,1) = SCATTER_PE
  receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,1) = &
              Offy*LOCAL_ROW_LEN+1+Offx

  IF (fld_type .EQ. fld_type_p) THEN
    receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) =  blsizep(2)
    receive_map(R_STRIDE_IN_RECV_ARRAY,1) = blsizep(1) + 2 * Offx
    receive_map(R_ELEMENT_LENGTH,1) = blsizep(1)
  ELSE IF (fld_type .EQ. fld_type_u) THEN
    receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) = blsizeu(2)
    receive_map(R_STRIDE_IN_RECV_ARRAY,1) = blsizeu(1) + 2 * Offx
    receive_map(R_ELEMENT_LENGTH,1) = blsizeu(1)
  ELSE
    receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) = blsizev(2)
    receive_map(R_STRIDE_IN_RECV_ARRAY,1) = blsizev(1) + 2 * Offx
    receive_map(R_ELEMENT_LENGTH,1) = blsizev(1)
  ENDIF

  receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,1) = &
              datastart(1)+(datastart(2)-1)*GLOBAL_ROW_LEN
  receive_map(R_STRIDE_IN_SEND_ARRAY,1) = GLOBAL_ROW_LEN


  old_GLOBAL_ROW_LEN=GLOBAL_ROW_LEN
  old_GLOBAL_ROWS=GLOBAL_ROWS
  old_PROC_GROUP=PROC_GROUP
  old_SCATTER_PE=SCATTER_PE
  old_DECOMP=current_decomp_type

ENDIF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

flag=GC_NONE  ! This is currently ignored at GCG v1.1
info=GC_NONE

CALL GCG_RALLTOALLE(GLOBAL_FIELD, send_map,n_mess_to_send, &
                    GLOBAL_ROW_LEN*GLOBAL_ROWS,            &
                    LOCAL_FIELD,  receive_map,  1,         &
                    LOCAL_ROW_LEN*LOCAL_ROWS,              &
                    PROC_GROUP, flag, info)


RETURN
End Subroutine Rcf_Scatter_Field_Log

End Module Rcf_Scatter_Field_Mod


