! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Provides client interfaces for async stash operations

MODULE IOS_stash

  USE IOS_common
  USE IOS_stash_common
  USE IOS_model_geometry
  USE IOS_geometry_utils
  USE IOS_Client_Queue
  USE IOS_Client_Coupler, ONLY :                                               &
      target_ios_rank
  USE yomhook, ONLY :                                                          &
      lhook, dr_hook
  USE parkind1, ONLY :                                                         &
      jprb, jpim

  IMPLICIT NONE

  INTEGER                :: handle
  INTEGER                :: IOS_block_size
  INTEGER                :: IOS_block_start
  LOGICAL                :: disable_subdomaining

  INTEGER, PARAMETER     :: async_tag_base            = 100000
  INTEGER, PARAMETER     :: IOS_async_max_handle      = 499
! params for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1
  CHARACTER (LEN=*), PARAMETER, PRIVATE  :: RoutineName = 'IOS_Stash_Client'
  CHARACTER (LEN=200), PRIVATE           :: IOS_st_message

CONTAINS

!
! Initialise client-side variables
!
  SUBROUTINE IOS_stash_client_init(stash_active,dump_active,                   &
      disable_subdomaining_in)
    USE mpl
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: stash_active
    LOGICAL, INTENT(IN) :: dump_active
    LOGICAL, INTENT(IN) :: disable_subdomaining_in
    INTEGER             :: i
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH:CLIENT_INIT',zhook_in,zhook_handle)
    CALL useAsyncStash(stash_active)
    CALL useAsyncDump(dump_active)
    disable_subdomaining=disable_subdomaining_in
    ! Never deallocated....

    ALLOCATE(slot(IOS_AsyncNumSlots))
    ALLOCATE(astash_mpi_requests(IOS_AsyncNumSlots))


    DO i=1,IOS_AsyncNumSlots
      CALL resetBuffer(i)
    END DO
    handle=0

    IF (lhook) CALL dr_hook('IOS_STASH:CLIENT_INIT',zhook_out,zhook_handle)

  END SUBROUTINE IOS_stash_client_init

!Just wait for all dispatch slots to complete.
  SUBROUTINE IOS_stash_client_fini()

    USE mpl, ONLY : mpl_status_size
    IMPLICIT NONE
    INTEGER              :: indx
    LOGICAL              :: msg
    INTEGER              :: ierror
    INTEGER              :: mplStat(mpl_status_size)
    REAL, EXTERNAL       :: get_wallclock_time
    REAL                 :: t1
    REAL                 :: t2
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH:CLIENT_FINI',                           &
        zhook_in,zhook_handle)

    msg=.FALSE.
    IF (IOS_Verbosity>=IOS_prstatus_diag)                                      &
        msg=.TRUE.

! DEPENDS ON : get_wallclock_time
    t1=Get_Wallclock_time()

    IF (IOS_Verbosity>=IOS_prstatus_oper)                                      &
        WRITE(6,'(A,A)')                                                       &
        'IOS: Info: IOS_Stash: ',                                              &
        'Waiting for all stash asynchronous requests to complete: '

    DO indx=1,IOS_AsyncNumSlots
      IF (slot(indx)%state==ios_queue_slot_unused) THEN
        IF (msg)                                                               &
            WRITE(6,'(A,I4,A,I3,A,I4,A,I3,A,I3)')                              &
            'Unused slot: ',indx,                                              &
            ' unit=',slot(indx)%unit,                                          &
            ' handle=',slot(indx)%handle,                                      &
            ' state=',slot(indx)%state,                                        &
            ' type=',slot(indx)%bufferType
      ELSE IF (slot(indx)%state==ios_queue_slot_dispatched) THEN
        IF (msg)                                                               &
            WRITE(6,'(A,I4,A,I3,A,I4,A,I3,A,I3)')                              &
            'Waiting for dispatched slot: ',indx,                              &
            ' unit=',slot(indx)%unit,                                          &
            ' handle=',slot(indx)%handle,                                      &
            ' state=',slot(indx)%state,                                        &
            ' type=',slot(indx)%bufferType

        CALL mpl_waitall(1,astash_mpi_requests(indx),mplStat,ierror)
      ELSE
        IF (msg)                                                               &
            WRITE(6,'(A,I4,A,I3,A,I4,A,I3,A,I3)')                              &
            'Not testing ',indx,                                               &
            ' unit=',slot(indx)%unit,                                          &
            ' handle=',slot(indx)%handle,                                      &
            ' state=',slot(indx)%state,                                        &
            ' type=',slot(indx)%bufferType
      END IF

    END DO

! DEPENDS ON : get_wallclock_time
    t2=Get_Wallclock_time()

    IF (IOS_Verbosity>=IOS_prstatus_oper)                                      &
        WRITE(6,'(A,A,F8.3)')                                                  &
        'IOS: Info: IOS_Stash: ',                                              &
        'Time waiting for all requests to complete: ',                         &
        t2-t1

    IF (lhook) CALL dr_hook('IOS_STASH:CLIENT_FINI',                           &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_stash_client_fini

  !
  !- Prepares a fresh buffer ready for packing.
  !- Prepared a fresh control structure for describing the data
  !
  SUBROUTINE IOS_stash_next_buffer( unit , s , buftype)

    USE mpl, ONLY : mpl_status_size
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit
    INTEGER, INTENT(OUT) :: s
    INTEGER, OPTIONAL    :: buftype
    INTEGER              :: ierror
    INTEGER              :: indx
    INTEGER              :: IOS_pe_for_unit
    LOGICAL              :: flag
    LOGICAL              :: msg
    LOGICAL              :: need_timer
    REAL(KIND=jprb)      :: zhook_handle
    REAL                 :: t1
    REAL                 :: t2
    INTEGER              :: mplStat(mpl_status_size)
    REAL, EXTERNAL       :: get_wallclock_time

    IF (lhook) CALL dr_hook('IOS_STASH:NEXT_BUFFER',zhook_in,zhook_handle)

    msg=.FALSE.
    IF (IOS_Verbosity>=IOS_prstatus_debug) msg=.TRUE.

    IF (msg)                                                                   &
        WRITE(6,'(A,A,I3)')'IOS: Info: ios_stash: ',                           &
        'Searching for free buffer for unit ',unit
! First check for outstanding buffers for this unit

    DO indx=1,IOS_AsyncNumSlots
      IF (slot(indx)%unit==unit) THEN
        IF (slot(indx)%state /= ios_queue_slot_dispatched) THEN
          WRITE(6,'(A,I3)')'Next buffer request for unit=',unit
          WRITE(6,'(A,I4,A)')                                                  &
              'But an undispatched slot ',indx,' for this unit'
          WRITE(6,'(A)')'exists!'
          CALL IOS_Ereport( RoutineName, 99 ,                                  &
              "Async dispatch slot screw up :-(" )
        END IF
      END IF
    END DO


! Are all slots dispatched?
! If not we have to wait for someone....
    s=-1
    indx=1

    ! IS there anything we can use or wait for?
    flag=.FALSE.
    DO WHILE (s==-1.AND.indx<=IOS_AsyncNumSlots)
      IF (slot(indx)%state==ios_queue_slot_dispatched)flag=.TRUE.
      IF (slot(indx)%state==ios_queue_slot_unused)flag=.TRUE.
      indx=indx+1
    END DO

    IF (flag) THEN

      flag=.FALSE.
      indx=1
      need_timer=.FALSE.

      DO WHILE (.NOT.flag)
        IF (slot(indx)%state==ios_queue_slot_unused) THEN
          IF (msg)                                                             &
              WRITE(6,'(A,I4,A,I3,A,I4,A,I3,A,I3)')                            &
              'unused slot: ',indx,                                            &
              ' unit=',slot(indx)%unit,                                        &
              ' handle=',slot(indx)%handle,                                    &
              ' state=',slot(indx)%state,                                      &
              ' type=',slot(indx)%bufferType
          flag=.TRUE.
          s=indx
        ELSE IF (slot(indx)%state==ios_queue_slot_dispatched) THEN
          IF (msg)                                                             &
              WRITE(6,'(A,I4,A,I3,A,I4,A,I3,A,I3)')                            &
              'testing dispatched slot: ',indx,                                &
              ' unit=',slot(indx)%unit,                                        &
              ' handle=',slot(indx)%handle,                                    &
              ' state=',slot(indx)%state,                                      &
              ' type=',slot(indx)%bufferType

          CALL mpl_testall(1,astash_mpi_requests(indx),flag,mplStat,ierror)
          IF (flag)s=indx
        ELSE
          IF (msg)                                                             &
              WRITE(6,'(A,I4,A,I3,A,I4,A,I3,A,I3)')                            &
              'not testing ',indx,                                             &
              ' unit=',slot(indx)%unit,                                        &
              ' handle=',slot(indx)%handle,                                    &
              ' state=',slot(indx)%state,                                      &
              ' type=',slot(indx)%bufferType
        END IF
        indx=indx+1
        IF (indx>IOS_AsyncNumSlots) THEN
          msg=.FALSE. ! shut of messages lest we fill the filesystem
          indx=1
          IF (.NOT.need_timer) THEN
! DEPENDS ON : get_wallclock_time
            t1=Get_Wallclock_time()
            need_timer=.TRUE.
          END IF
        END IF

        IF (need_timer) THEN
! DEPENDS ON : get_wallclock_time
          t2=Get_Wallclock_time()
          IF (t2-t1 > IOS_Timeout) THEN
            WRITE(6,'(A)')'Timeout value tripped waiting for stash'
            WRITE(IOS_St_message,'(A)')                                        &
                'Timeout value tripped waiting for stash'
            CALL IOS_Ereport( RoutineName, 99 , IOS_st_message )

          END IF
        END IF
      END DO

      IF (need_timer) THEN
! DEPENDS ON : get_wallclock_time
        t2=Get_Wallclock_time()
        IF (IOS_Verbosity>=IOS_prstatus_oper)                                  &
            WRITE(6,'(A,A,F10.3,A,F8.3)')                                      &
            'IOS: Info: IOS_Stash: ',                                          &
            'Stall getting asynchronous stash slot at ',                       &
            t1-IOS_Start_time,' of ',t2-t1
      END IF

    ELSE
      IF (IOS_Verbosity>=IOS_prstatus_normal)                                  &
          WRITE(6,'(A)')                                                       &
          'IOS_Stash: no slots in dispatched or unused state'
    END IF

    IF (s>0) THEN
      IF (IOS_Verbosity>=IOS_prstatus_debug)                                   &
          WRITE(6,'(A,A,I3,A)')'IOS: Info: ios_stash: ',                       &
          'Slot ',s,' is now free.'

      handle=handle+1
      IF (handle > IOS_async_max_handle)handle=1

      CALL resetBuffer(s)
      slot(s)%unit   =unit
      slot(s)%handle =handle
      slot(s)%state  =ios_queue_slot_initialized
      IF (PRESENT(bufType))THEN
        slot(s)%bufferType=bufType
      ELSE
        slot(s)%bufferType=IOS_Action_StashWritePPData
      END IF

    ELSE

      IF (IOS_Verbosity>=IOS_prstatus_normal) THEN
        WRITE(6,'(A,A)')'IOS: Info: IOS_stash: ',                              &
            'No slot could be freed, something needs dispatching'
        WRITE(6,'(A)')'IOS: Info: IOS_stash:   or more slots allocating'
      END IF

    END IF

    IF (lhook) CALL dr_hook('IOS_STASH:NEXT_BUFFER',zhook_out,zhook_handle)

  END SUBROUTINE IOS_stash_next_buffer

  SUBROUTINE resetBuffer(i)
    USE mpl, ONLY : mpl_request_null
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    slot(i)%buffer_offset       =0
    slot(i)%control_offset      =0
    slot(i)%num_records_in_pack =0
    slot(i)%state               =ios_queue_slot_unused
    slot(i)%tag                 =-1
    slot(i)%unit                =-1
    slot(i)%handle              =-1
    slot(i)%time_create         =0
    slot(i)%time_update         =0
    slot(i)%time_complete       =0
    slot(i)%bufferType          =-1
    astash_mpi_requests(i)      =mpl_request_null
  END SUBROUTINE resetBuffer

  FUNCTION IOS_getCurrentSlot(unit) RESULT(s)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER :: s
    INTEGER :: indx

    s=-1
    indx=1
    DO WHILE (s==-1.AND.indx<=IOS_AsyncNumSlots)
      IF (slot(indx)%unit==unit) THEN
        IF (                                                                   &
            slot(indx)%state==ios_queue_slot_initialized .OR.                  &
            slot(indx)%state==ios_queue_slot_partfilled                        &
            ) s=indx
      END IF
      indx=indx+1
    END DO

  END FUNCTION IOS_getCurrentSlot

  FUNCTION IOS_getLevsInPack(s) RESULT(levs)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: s
    INTEGER             :: levs

    levs=slot(s)%num_records_in_pack

  END FUNCTION IOS_getLevsInPack

  FUNCTION IOS_getSlotUnit(s) RESULT(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: s
    INTEGER             :: unit

    unit=slot(s)%unit

  END FUNCTION IOS_getSlotUnit

  FUNCTION IOS_getSlotState(s) RESULT(state)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: s
    INTEGER             :: state

    state=slot(s)%state

  END FUNCTION IOS_getSlotState

  SUBROUTINE IOS_Stash_Discard_Buffer(s)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: s
    CALL ResetBuffer(s)
  END SUBROUTINE IOS_Stash_Discard_Buffer

!
! Pack a stash dispatch buffer
!
! Note in this version we send whole fields, this is to simplify things and
! get the method robust before we start adding complexity. We munge the data
! such that we have a whole tile that we can logically forward to the server
! even if there is no relevant data on the tile.
!
  SUBROUTINE IOS_stash_pack_pp_data (                                          &
      s,                                                                       &
      array,                                                                   &
      arrayLen,                                                                &
      arg_fld_type,                                                            &
      arg_halo_type,                                                           &
      arg_preprocess_flag,                                                     &
      arg_packing_flag,                                                        &
      arg_subdomain_flag,                                                      &
      arg_pack_type,                                                           &
      arg_comp_accry,                                                          &
      arg_dmi,                                                                 &
      arg_disk_block,                                                          &
      arg_seekpos,                                                             &
      arg_S_boundary,                                                          &
      arg_N_boundary,                                                          &
      arg_W_boundary,                                                          &
      arg_E_boundary,                                                          &
      arg_section,                                                             &
      arg_code,                                                                &
      arg_level                                                                &
      )

    USE UM_Parvars, ONLY :                                                     &
        halosize,                                                              &
        blsize,                                                                &
        glsize,                                                                &
        lasize

    IMPLICIT NONE

    INTEGER,INTENT(IN)              :: s!slot
    REAL   ,INTENT(IN)              :: array(:)
    INTEGER,INTENT(IN)              :: arrayLen
    INTEGER,          INTENT(IN)    :: arg_fld_type
    INTEGER,          INTENT(IN)    :: arg_halo_type
    INTEGER,          INTENT(IN)    :: arg_preprocess_flag
    INTEGER,          INTENT(IN)    :: arg_subdomain_flag
    INTEGER,          INTENT(IN)    :: arg_packing_flag
    INTEGER,          INTENT(IN)    :: arg_S_boundary
    INTEGER,          INTENT(IN)    :: arg_N_boundary
    INTEGER,          INTENT(IN)    :: arg_W_boundary
    INTEGER,          INTENT(IN)    :: arg_E_boundary
    INTEGER,          INTENT(IN)    :: arg_pack_type
    INTEGER,          INTENT(IN)    :: arg_comp_accry
    REAL,             INTENT(IN)    :: arg_dmi
    INTEGER,          INTENT(IN)    :: arg_disk_block
    INTEGER,          INTENT(IN)    :: arg_seekpos
    INTEGER,          INTENT(IN)    :: arg_section
    INTEGER,          INTENT(IN)    :: arg_code
    INTEGER,          INTENT(IN)    :: arg_level

    INTEGER(KIND=integer64), ALLOCATABLE :: control(:)
    INTEGER                         :: halo_EW
    INTEGER                         :: halo_NS

! In this version we will subdomain on the server, but ultimately we do
! not want to send unneeded data across the net, so the subdomain array
! will be polulated with the whole of the local field, but omitting halos.
! local_len will be computed as the size of the data transmitted.
    REAL, ALLOCATABLE               :: subdomain(:)
    INTEGER                         :: local_len
    INTEGER                         :: subdomain_type
    TYPE(box)                       :: domainBox
    TYPE(box)                       :: localBox
    TYPE(box)                       :: localBox_W
    TYPE(box)                       :: localBox_E

! General local indices/bounds for loops
    INTEGER                         :: col
    INTEGER                         :: src_row
    INTEGER                         :: src_point_row_start
    INTEGER                         :: target_pnt
    INTEGER                         :: i

    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH:PACK_PP_DATA',zhook_in,zhook_handle)

    ALLOCATE(control(1:IOS_async_control_sz_max))
    control(:)=IOS_unset_int
    control(loc_record_type)=IOS_stash_distributed_field

    control(loc_fld_type)         =arg_fld_type
    control(loc_preprocess_flag)  =arg_preprocess_flag
    control(loc_seek_target)      =arg_seekpos

    ! Check the array size against the claim
    IF (arrayLen /= SIZE(array)) THEN !someone lied
      WRITE(IOS_st_message,'(A,I14,A,I14,A)')                                  &
          'mismatch in size between array length(',                            &
          SIZE(array),') and len argument(',arrayLen,')'
      WRITE(6,'(A)')IOS_st_message
      WRITE(6,'(A,A)')                                                         &
          'Please pass in the correctly sized rank-1 array ',                  &
          '(and not the sloppy f77 thing you just did)'
      CALL IOS_Ereport(RoutineName, 99, IOS_st_message)
    END IF

    ! Check for required arguments for a preprocessed field
    IF (arg_preprocess_flag==IOS_stash_preprocess) THEN
      control(loc_packing_flag)      = arg_packing_flag
      control(loc_subdomain_flag)    = arg_subdomain_flag
      control(loc_landmaskcompress)  = IOS_packing_type_nolandmask
    END IF

    ! Compute the halos and thus the relevant data size/address
    halo_EW=halosize(1,arg_halo_type)
    halo_NS=halosize(2,arg_halo_type)
    local_len=                                                                 &
        blsize(1,arg_fld_type)*blsize(2,arg_fld_type)

    IF (arg_subdomain_flag==IOS_partial_field) THEN
      domainBox%n = arg_N_boundary
      domainBox%s = arg_S_boundary
      domainBox%e = arg_E_boundary
      domainBox%w = arg_W_boundary
    ELSE
      domainBox%n = glsize(2,arg_fld_type)
      domainBox%s = 1
      domainBox%e = glsize(1,arg_fld_type)
      domainBox%w = 1
    END IF

    CALL IOS_pack4(control(loc_boundary),                                      &
        domainBox%n,                                                           &
        domainBox%s,                                                           &
        domainBox%e,                                                           &
        domainBox%w)

    CALL IOS_pack4(control(loc_id),                                            &
        arg_section,                                                           &
        arg_code,                                                              &
        arg_level,                                                             &
        0)

    ALLOCATE(subdomain(blsize(1,arg_fld_type)*blsize(2,arg_fld_type)))

    ! If subdomaining is off ensure that we initialise the grid to dmi.
    IF (disable_subdomaining) THEN
      subdomain(:)=arg_dmi
    END IF

    IF (arg_subdomain_flag==IOS_partial_field) THEN
      subdomain_type=classify_subdomain(model_rank,                            &
          arg_fld_type,domainBox)
    ELSE
      subdomain_type=complete_intersection
    END IF

    IF (arg_subdomain_flag==IOS_full_field) THEN
      ! The input buffer is the whole tile, so strip any halos and copy it
      !to the output buffer
      target_pnt=1
      DO src_row=halo_NS+1,blsize(2,arg_fld_type)+halo_NS
        src_point_row_start=(src_row-1)*lasize(1,arg_fld_type,arg_halo_type)   &
            +1+halo_EW
        DO col=1,blsize(1,arg_fld_type)
          subdomain(target_pnt)=array(src_point_row_start+col-1)
          target_pnt=target_pnt+1
        END DO
      END DO
      IF(target_pnt-1 > SIZE(subdomain)) THEN
        CALL IOS_Ereport('pack_pp',98,                                         &
            'memory overwrite in copying to halo-free buffer')
      END IF
    ELSE IF (subdomain_type==complete_intersection) THEN
      ! The input buffer is the whole tile, but halo free so
      ! copy it to the output buffer
      target_pnt=1
      DO src_row=1,blsize(2,arg_fld_type)
        src_point_row_start=(src_row-1)*blsize(1,arg_fld_type)+1
        DO col=1,blsize(1,arg_fld_type)
          subdomain(target_pnt)=array(src_point_row_start+col-1)
          target_pnt=target_pnt+1
        END DO
      END DO
      IF(target_pnt-1 > SIZE(subdomain)) THEN
        CALL IOS_Ereport('pack_pp',99,                                         &
            'memory overwrite in copying to halo-free buffer')
      END IF
    ELSE ! partial field that doesn"t overlap completely with this proc
      IF (subdomain_type==1) THEN
        ! there is one encroachment of the subdomain onto my tile

        ! get the local coords in my domain to copy into
        CALL getLocalSubdomainBounds                                           &
            (model_rank, arg_fld_type,domainBox,localBox)
        ! Debug
        ! subdomain=10.0
        target_pnt=1
        DO src_row=localBox%s,localBox%n
          DO col=localBox%w,localBox%e
            subdomain(col+(src_row-1)*blsize(1,arg_fld_type))=                 &
                array(target_pnt)
            target_pnt=target_pnt+1
          END DO
        END DO
      ELSE IF (subdomain_type==2) THEN
        ! there are 2 bits of the requested tile on my processor,
        ! this happens when a wrap-around subdomain starts on my processor
        ! wraps all the way around the global model and ends on my processor
        ! or the west boundary is more easterly than the west boundary and
        ! both boundaries belong to me.

        target_pnt=1
        ! Get the westerly part of the subdomain
        ! (eastmost in the local processor domain)
        CALL getLocalSubdomainBoundsCyclic                                     &
            (model_rank,arg_fld_type,                                          &
            domainBox,localBox_W,WesterlyZone)
        ! Get the easterly part of the subdomain
        ! (westmost in the local processor domain)
        CALL getLocalSubdomainBoundsCyclic                                     &
            (model_rank,arg_fld_type,                                          &
            domainBox,localBox_E,EasterlyZone)
        ! The north south localBox boundaries should be the same
        ! in both cases.

        DO src_row=localBox_W%s,localBox_W%n
          DO col=localBox_W%w,localBox_W%e
            subdomain(col+(src_row-1)*blsize(1,arg_fld_type))=                 &
                array(target_pnt)
            target_pnt=target_pnt+1
          END DO
          DO col=localBox_E%w,localBox_E%e
            subdomain(col+(src_row-1)*blsize(1,arg_fld_type))=                 &
                array(target_pnt)
            target_pnt=target_pnt+1
          END DO
        END DO
        ! Debug
        ! ELSE IF (subdomain_type==no_intersection) THEN
        ! subdomain=20.0
      END IF
    END IF

    IF (arg_preprocess_flag==IOS_stash_preprocess) THEN
      IF (arg_packing_flag==IOS_packing) THEN
        control(loc_pack_type)   =  IOS_packing_type_wgdos
        control(loc_comp_accry)   = arg_comp_accry
      END IF
      CALL UM_memcpy_f(control(loc_dmi) , arg_dmi, 1)

      control(loc_disk_block)=arg_disk_block
    END IF ! We need to preprocess

    ! Now that we have a a homogenised tile with no halos, we pass it on.
    CALL IOS_pack_async_buffer                                                 &
        (s,subdomain,local_len,control,IOS_async_control_sz_max,subdomain_type)

    DEALLOCATE(subdomain)
    DEALLOCATE(control)

    IF (lhook) CALL dr_hook('IOS_STASH:PACK_PP_DATA',zhook_out,zhook_handle)

  END SUBROUTINE IOS_stash_pack_pp_data

!
! Receive a homogenised tile, together with metadata for the tile, and
! add it into our aggregated dispatch buffer.
!

  SUBROUTINE IOS_pack_async_buffer( s , array, arrayLen, cont, clen, domType)
    ! array: data to be packed into send buffers, the control
    ! record will be updated
    IMPLICIT NONE
    INTEGER,INTENT(IN)     :: s        !slot
    REAL   ,INTENT(IN)     :: array(:) !the tile
    INTEGER(KIND=integer64),INTENT(IN)     :: cont(:)  !metadata
    INTEGER,INTENT(IN)     :: arrayLen !data length
    INTEGER,INTENT(IN)     :: clen     !metadata length
    INTEGER,INTENT(IN)     :: domType  !domain classification

    INTEGER                :: i
    INTEGER                :: j
    INTEGER                :: err
    LOGICAL                :: isRepeat
    INTEGER, SAVE          :: lastRecord
    INTEGER                :: lastRecordLen
    REAL(KIND=jprb)        :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH:PACK_BUFFER',zhook_in,zhook_handle)

    slot(s)%num_records_in_pack=slot(s)%num_records_in_pack+1

    ! Determine if this request is identical to the previous one....
    IF (slot(s)%num_records_in_pack==1) THEN
      isRepeat=.FALSE.
    ELSE
      lastRecordLen=-1*slot(s)%control(lastRecord+4)
      IF (clen==lastRecordLen) THEN
        isRepeat=.TRUE.
        DO i=1,clen
          IF(slot(s)%control(lastRecord+i+5) /= cont(i)) THEN
            isRepeat=.FALSE.
          END IF
        END DO
      ELSE
        isRepeat=.FALSE.
      END IF
    END IF

    IF (isRepeat) THEN

! Check we have space...
      IF (slot(s)%control_offset+ios_stash_control_auto_len+1 >                &
          SIZE(slot(s)%control)) THEN
        WRITE(IOS_st_message,'(A,I4,A,I6,A,I10)')                              &
            'Adding ',ios_stash_control_auto_len+1,                            &
            'elements to ',slot(s)%control_offset,                             &
            ' exceeds control size of ',                                       &
            SIZE(slot(s)%control)
        CALL IOS_Ereport( RoutineName, 98 , IOS_st_message )
      END IF

      slot(s)%control(slot(s)%control_offset+loc_record_start)=                &
          IOS_stash_record_start
      slot(s)%control(slot(s)%control_offset+loc_record_len_control)=          &
          -1! record control len

      slot(s)%control(slot(s)%control_offset+                                  &
          ios_stash_control_auto_len+1)=                                       &
          ios_repeat_record

      slot(s)%control_offset=slot(s)%control_offset+                           &
          ios_stash_control_auto_len+1 ! +1 for the repeat record.

    ELSE

      lastRecord=slot(s)%control_offset+1

! Check we have space...
      IF (slot(s)%control_offset+clen+ios_stash_control_auto_len >             &
          SIZE(slot(s)%control)) THEN
        WRITE(IOS_st_message,'(A,I4,A,I6,A,I10)')                              &
            'Adding ',clen+ios_stash_control_auto_len,                         &
            'elements to ',slot(s)%control_offset,                             &
            ' exceeds control size of ',                                       &
            SIZE(slot(s)%control)
        CALL IOS_Ereport( RoutineName, 99 , IOS_st_message )
      END IF

! Fill out the mandatory bit of the metadata (metametadata)
      slot(s)%control(slot(s)%control_offset+loc_record_start)=                &
          IOS_stash_record_start
      slot(s)%control(slot(s)%control_offset+loc_record_len_control)=          &
          -clen ! record control length

      slot(s)%control_offset=                                                  &
          slot(s)%control_offset+ios_stash_control_auto_len

      slot(s)%control(slot(s)%control_offset+1:                                &
          slot(s)%control_offset+clen) = cont(1:clen)

      slot(s)%control_offset=slot(s)%control_offset+clen
    END IF

    ! Check we have space for data....
    IF (slot(s)%buffer_offset+arrayLen>SIZE(slot(s)%buffer)) THEN

      WRITE(6,'(A)'     )'ATTEMPTED BUFFER OVERRUN!'
      WRITE(6,'(A,I6)'  )'adding field #',slot(s)%num_records_in_pack+1
      WRITE(6,'(A,I5,A)')'to buffer of slot ',s,' OUCH!'
      WRITE(6,'(A,I7)'  )'each field expected to have max size of ',           &
          MaxFieldDomain
      WRITE(IOS_st_message,'(A,I8,A,I10,A,I10)')                               &
          'Adding ',arrayLen, 'elements to ',                                  &
          slot(s)%buffer_offset,                                               &
          ' exceeds buffer size of ',SIZE(slot(s)%buffer)
      CALL IOS_Ereport( RoutineName, 99 , IOS_st_message )
    END IF


    IF (DomType /= no_intersection .OR. IOS_AsyncSendNull) THEN
      slot(s)%buffer(slot(s)%buffer_offset+1:slot(s)%buffer_offset+arrayLen) = &
          array(1:arrayLen)
      slot(s)%buffer_offset=slot(s)%buffer_offset+arrayLen
    END IF

    ! Fix up the status of the slot
    slot(s)%state=ios_queue_slot_partfilled

    ! Everything is nicely stored away in the slot,
    ! so we can go back to running the model.
    IF (lhook) CALL dr_hook('IOS_STASH:PACK_BUFFER',zhook_out,zhook_handle)
  END SUBROUTINE IOS_pack_async_buffer

!
! Send a control message block to IOS
!
  SUBROUTINE IOS_stash_write_pp_data( s )
    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: s!slot
    INTEGER                    :: tag
    INTEGER                    :: ierr
    INTEGER                    :: targetRank
    INTEGER                    :: qHandle
    TYPE(IOS_metadata_type),                                                   &
        POINTER                :: metadata
    INTEGER(KIND=ios_tukind),                                                  &
        POINTER                :: sendBuffer(:)

    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH:WRITE_PP_DATA',zhook_in,zhook_handle)

    IF (slot(s)%bufferType /= IOS_Action_StashWritePPData)THEN
      CALL IOS_Ereport("stash_write_pp_data", 99 ,                             &
          "write PP data called with wrong buffer type" )
    END IF

    IF (model_rank==0) THEN
      DO targetRank=1,IOS_tasks_per_server
        qHandle = ios_init_md(                                                 &
            slot(s)%unit,-1,                                                   &
            IOS_Action_StashWritePPData,                                       &
            dataSize=slot(s)%control_offset,rank=targetRank)
        metadata => IOS_Attach_Metadata(qHandle)
        metadata % handle    = slot(s)%handle
        sendBuffer => IOS_attach_SendBuffer(qHandle)
        CALL um_memcpy64(sendBuffer,                                           &
            slot(s)%control,                                                   &
            slot(s)%control_offset)
        CALL IOS_Send(qHandle)
      END DO
    END IF

    IF (lhook) CALL dr_hook('IOS_STASH:WRITE_PP_DATA',zhook_out,zhook_handle)
  END SUBROUTINE IOS_stash_write_pp_data


  SUBROUTINE IOS_stash_dispatch( s )
    USE mpl, ONLY :                                                            &
        mpl_real,                                                              &
        mpl_request_null
    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: s ! slot
    INTEGER                    :: tag
    INTEGER                    :: ierr
    INTEGER                    :: targetPE
    INTEGER                    :: targetServer
    INTEGER                    :: targetRank
    TYPE(IOS_metadata_type)    :: metadata
    REAL(KIND=jprb)            :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH:DISPATCH',zhook_in,zhook_handle)

    targetPE = io_server_for_unit(slot(s)%unit)
    targetServer=ioServerNo(targetPe)
    targetRank=io_servers(targetServer,target_ios_rank+1)

    tag=slot(s)%handle

    IF (IOS_Verbosity>=IOS_prstatus_diag)                                      &
        WRITE(6,'(A,A,I3,A,I3,A,I4,A,I3,A,I4)')                                &
        'IOS: Info: IOS_stash: ',                                              &
        'Dispatch: unit=',slot(s)%unit,                                        &
        ' slt=',s,' hndle=',slot(s)%handle,                                    &
        ' flds=',slot(s)%num_records_in_pack,                                  &
        ' srvr=',targetRank

    IF (slot(s)%buffer_offset > 0 ) THEN
      CALL MPL_Isend(slot(s)%buffer,                                           &
          slot(s)%buffer_offset,                                               &
          MPL_Real,                                                            &
          targetRank,                                                          &
          tag,                                                                 &
          IOS_async_comm,                                                      &
          astash_mpi_requests(s),                                              &
          ierr )
    ELSE
      astash_mpi_requests(s)=mpl_request_null
    END IF

    slot(s)%state=ios_queue_slot_dispatched

    IF (lhook) CALL dr_hook('IOS_STASH:DISPATCH',zhook_out,zhook_handle)

  END SUBROUTINE IOS_stash_dispatch

END MODULE IOS_stash
