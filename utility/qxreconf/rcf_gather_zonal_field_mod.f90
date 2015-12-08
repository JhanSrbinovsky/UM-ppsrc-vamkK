! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! + Gathers a zonal field from many processors to one processor

Module Rcf_Gather_Zonal_Field_Mod

!  Subroutine Rcf_GAther_Zonal_Field -  Gathers a field onto 1 pe
!
! Description:
! Takes a zonal field decomposed on many processors and gathers
! it onto a single specified PE.
!
! Method:
!  Calculates send and receive maps to use with GCG_RALLETOALLE
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Gather_Zonal_Field ( LOCAL_FIELD, GLOBAL_FIELD,  &
                                    LOCAL_SIZE,  GLOBAL_SIZE,   &
                                    LEVELS,      GRID_TYPE,     &
                                    GATHER_PE)

USE UM_ParVars    ! Most of this used
USE gcom_mod
USE cppxref_mod, ONLY : ppx_atm_tzonal, ppx_atm_ozone

IMPLICIT NONE

! Arguments
Integer, Intent(In)   :: Local_Size
Integer, Intent(In)   :: Global_Size
Integer, Intent(In)   :: Levels
Integer, Intent(In)   :: Grid_Type
Integer, Intent(In)   :: Gather_PE
Real, Intent(In)      :: Local_Field ( Local_Size, Levels )
Real, Intent(Out)     :: Global_Field( Global_Size, Levels )

! Local variables
Integer               :: k                 ! looper
Integer               :: flag              ! Gcom argument
Integer               :: fld_type          ! P or U field
Integer               :: info              ! Gcom return code
Integer               :: iproc             ! loop counter
Integer               :: send_map(7,1)
INTEGER, ALLOCATABLE,SAVE  :: receive_map(:,:)
Integer               :: n_mess_to_send
Integer               :: n_mess_to_receive

!====================================================================

IF (.NOT.ALLOCATED(receive_map)) &
    ALLOCATE(receive_map(7,NPROC_MAX))

send_map(:,:)    = 0
receive_map(:,:) = 0

! Note no V zonal grids....
IF ((grid_type .EQ. ppx_atm_tzonal) .OR. &
    (grid_type .EQ. ppx_atm_ozone)) THEN
  fld_type=fld_type_p
ELSE
  fld_type=fld_type_u
ENDIF

!--------------------------------------------------------------------

n_mess_to_receive=0

IF (mype .EQ. GATHER_PE) THEN
  DO iproc=0,nproc-1
    IF (g_gridpos(1,iproc) .EQ. 0) THEN
!           Only one processor per LPG row needs to send the data
!           as it will be the same for each processor along the
!           row.
      receive_map(R_SOURCE_PE,n_mess_to_receive+1) = iproc
      receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY, n_mess_to_receive+1) = &
                  g_datastart(2,iproc)
      receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM, n_mess_to_receive+1) = 1
      receive_map(R_STRIDE_IN_RECV_ARRAY, n_mess_to_receive+1) = 0
      IF (fld_type .EQ. fld_type_p) THEN
        receive_map(R_ELEMENT_LENGTH,n_mess_to_receive+1) = &
                    g_blsizep(2,iproc)
      ELSE
        receive_map(R_ELEMENT_LENGTH,n_mess_to_receive+1) = &
                    g_blsizeu(2,iproc)
      ENDIF
      receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY, &
                  n_mess_to_receive+1) = Offy+1
      receive_map(R_STRIDE_IN_SEND_ARRAY, n_mess_to_receive+1) = 0
      n_mess_to_receive=n_mess_to_receive+1
    ENDIF
  ENDDO
ENDIF

n_mess_to_send=0
  IF (atwest) THEN ! only processors at the left of the LPG will
                   ! send anything
    send_map(S_DESTINATION_PE,1) = GATHER_PE
    send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,1) = Offy+1
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1) = 1
    send_map(S_STRIDE_IN_SEND_ARRAY,1) = 0
    IF (fld_type .EQ. fld_type_p) THEN
      send_map(S_ELEMENT_LENGTH,1) = blsizep(2)
    ELSE
     send_map(S_ELEMENT_LENGTH,1) = blsizeu(2)
    ENDIF
    send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = datastart(2)
    send_map(S_STRIDE_IN_RECV_ARRAY,1) = 0

    n_mess_to_send=1
  ENDIF


DO k=1,LEVELS

  info=GC_NONE
  flag=GC_NONE

  IF (fld_type .EQ. fld_type_p) THEN
    CALL GCG_RALLTOALLE(                                 &
        LOCAL_FIELD(1,k),send_map,n_mess_to_send,        &
        lasize(2,fld_type_p,halo_type_single),          &
        GLOBAL_FIELD(1,k),receive_map,n_mess_to_receive, &
        glsizep(2),GC_ALL_PROC_GROUP,flag,info)
  ELSE
    CALL GCG_RALLTOALLE(                                 &
        LOCAL_FIELD(1,k),send_map,n_mess_to_send,        &
        lasize(2,fld_type_p,halo_type_single),          &
        GLOBAL_FIELD(1,k),receive_map,n_mess_to_receive, &
        glsizep(2)-1,GC_ALL_PROC_GROUP,flag,info)
  ENDIF

ENDDO

Return

END Subroutine Rcf_Gather_Zonal_Field
End Module Rcf_Gather_Zonal_Field_Mod


