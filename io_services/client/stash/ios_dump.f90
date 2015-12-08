! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Provides client interfaces for async stash operations

MODULE ios_dump
  USE ios_stash
  USE cppxref_mod, ONLY: ppx_atm_tzonal,ppx_atm_uzonal
  USE um_types
  IMPLICIT NONE
  INTEGER                :: IOS_dump_pack_flag
  INTEGER                :: IOS_dump_pack_type
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1
  CHARACTER (LEN=*), PARAMETER, PRIVATE  :: RoutineName = 'IOS_Dump_Client'

CONTAINS

!
! Pack a dump dispatch buffer
!
  SUBROUTINE IOS_dump_pack_data (                                              &
      s,                                                                       &
      array,                                                                   &
      arg_grid_type,                                                           &
      arg_halo_type,                                                           &
      arg_preprocess_flag,                                                     &
      arg_subdomain_flag,                                                      &
      arg_seekpos,                                                             &
      arg_S_boundary,                                                          &
      arg_N_boundary,                                                          &
      arg_W_boundary,                                                          &
      arg_E_boundary,                                                          &
      arg_landmask_compress                                                    &
      )

    USE UM_Parvars, ONLY :                                                     &
        halosize,                                                              &
        blsize,                                                                &
        glsize,                                                                &
        lasize

    IMPLICIT NONE

    INTEGER,           INTENT(IN)    :: s!slot
    REAL   ,           INTENT(IN)    :: array(:)
    INTEGER, OPTIONAL, INTENT(IN)    :: arg_grid_type
    INTEGER, OPTIONAL, INTENT(IN)    :: arg_halo_type
    INTEGER, OPTIONAL, INTENT(IN)    :: arg_preprocess_flag
    INTEGER, OPTIONAL, INTENT(IN)    :: arg_subdomain_flag
    INTEGER, OPTIONAL, INTENT(IN)    :: arg_S_boundary
    INTEGER, OPTIONAL, INTENT(IN)    :: arg_N_boundary
    INTEGER, OPTIONAL, INTENT(IN)    :: arg_W_boundary
    INTEGER, OPTIONAL, INTENT(IN)    :: arg_E_boundary
    INTEGER, OPTIONAL, INTENT(IN)    :: arg_seekpos
    LOGICAL, OPTIONAL, INTENT(IN)    :: arg_landmask_compress

    INTEGER(KIND=integer64), ALLOCATABLE             :: control(:)
    INTEGER                          :: halo_EW
    INTEGER                          :: halo_NS
    INTEGER                          :: fld_type
    INTEGER                          :: get_fld_type

! In this version we will subdomain on the server, but ultimately we do
! not want to send unneeded data across the net, so the subdomain array
! will be populated with the whole of the local field, but omitting halos.
! local_len will be computed as the size of the data transmitted.
    REAL, ALLOCATABLE                :: subdomain(:)
    INTEGER                          :: local_len
    INTEGER                          :: lngth
    INTEGER                          :: subdomain_type
    TYPE(box)                        :: domainBox
    TYPE(box)                        :: localBox
    TYPE(box)                        :: localBox_W
    TYPE(box)                        :: localBox_E

! General local indices/bounds for loops
    INTEGER                          :: col
    INTEGER                          :: src_row
    INTEGER                          :: src_point_row_start
    INTEGER                          :: target_pnt
    INTEGER                          :: i

    REAL(KIND=jprb)                  :: zhook_handle
    CHARACTER (LEN=200)              :: IOS_dmp_message

    IF (lhook) CALL dr_hook('IOS_DUMP:PACK_DUMP_DATA',                         &
        zhook_in,zhook_handle)

    ALLOCATE(control(1:IOS_async_control_sz_max))

    control(:)=IOS_unset_int
    control(loc_record_type)=IOS_stash_distributed_field
    IF  (.NOT.PRESENT(arg_grid_type) .OR.                                      &
        .NOT.PRESENT(arg_halo_type) .OR.                                       &
        .NOT.PRESENT(arg_preprocess_flag) ) THEN
      WRITE(IOS_dmp_message,'(A)')                                             &
          'Missing attributes for a distributed field'
      CALL IOS_Ereport( RoutineName, 99 , IOS_dmp_message )
    END IF

    fld_type = get_fld_type(arg_grid_type)
    control(loc_fld_type)         =fld_type
    control(loc_preprocess_flag)  =arg_preprocess_flag
    ! Check for required arguments for a preprocessed field
    IF (arg_preprocess_flag==IOS_stash_preprocess) THEN
      IF  (.NOT.PRESENT(arg_subdomain_flag)) THEN
        WRITE(IOS_dmp_message,'(A,A)')                                         &
            'Missing attributes for preprocessed ',                            &
            'field packing/subdomain'
        CALL IOS_Ereport( RoutineName, 99 , IOS_dmp_message )
      END IF
      control(loc_subdomain_flag) = arg_subdomain_flag
    END IF
    control(loc_packing_flag)     = IOS_dump_pack_flag
    control(loc_pack_type)        = IOS_dump_pack_type
    control(loc_landmaskcompress) = IOS_packing_type_nolandmask

    IF (PRESENT(arg_landmask_compress)) THEN
      IF (arg_landmask_compress)                                               &
          control(loc_landmaskcompress) = IOS_packing_type_landmask
    END IF

    IF (IOS_dump_pack_type /= IOS_No_Packing) THEN
      IF (IOS_dump_pack_flag == IOS_No_Packing) THEN
        WRITE(IOS_dmp_message,'(A)')                                           &
            'Packing type is set, but packing flag is not set.'
        CALL IOS_Ereport( RoutineName, 99 , IOS_dmp_message )
      END IF
    END IF

    IF (control(loc_subdomain_flag)==IOS_partial_field) THEN
      IF (.NOT.PRESENT(arg_N_boundary) .OR.                                    &
          .NOT.PRESENT(arg_S_boundary) .OR.                                    &
          .NOT.PRESENT(arg_E_boundary) .OR.                                    &
          .NOT.PRESENT(arg_W_boundary)) THEN
        WRITE(IOS_dmp_message,'(A)')                                           &
            'Missing attributes for subdomain: N/S/E/W'
        CALL IOS_Ereport( RoutineName, 99 , IOS_dmp_message )
      END IF

      domainBox%n = arg_N_boundary
      domainBox%s = arg_S_boundary
      domainBox%e = arg_E_boundary
      domainBox%w = arg_W_boundary
    ELSE
      domainBox%n = glsize(2,fld_type)
      domainBox%s = 1
      domainBox%e = glsize(1,fld_type)
      domainBox%w = 1
    END IF

    ! Correction for zonal fields
    IF ((arg_grid_type  ==  ppx_atm_tzonal) .OR.                               &
        (arg_grid_type  ==  ppx_atm_uzonal)) THEN
      domainBox%e               = 1
      domainBox%w               = 1
      control(loc_subdomain_flag)=ios_partial_field
    END IF

    CALL IOS_pack4(control(loc_boundary),                                      &
        domainBox%n,                                                           &
        domainBox%s,                                                           &
        domainBox%e,                                                           &
        domainBox%w)

    ! Compute the halos and thus the relevant data size/address
    halo_EW=halosize(1,arg_halo_type)
    halo_NS=halosize(2,arg_halo_type)
    local_len=                                                                 &
        blsize(1,fld_type)*blsize(2,fld_type)

    lngth=(blsize(1,fld_type)+2*halo_EW)*                                      &
        (blsize(2,fld_type)+2*halo_NS)

    ! Check the array size against the claim
    IF (lngth > SIZE(array).AND.                                               &
        control(loc_subdomain_flag)==IOS_full_field) THEN
      !someone lied
      WRITE(IOS_dmp_message,'(A,I14,A,I14)')                                   &
          'mismatch in size between array length(',                            &
          SIZE(array),') and data size (',lngth,')'
      CALL IOS_Ereport(RoutineName, 99, IOS_dmp_message)
    END IF

    ALLOCATE(subdomain(blsize(1,fld_type)*blsize(2,fld_type)))

    IF (control(loc_subdomain_flag)==IOS_partial_field) THEN
      subdomain_type=classify_subdomain(model_rank,                            &
          fld_type,domainBox)
    ELSE
      subdomain_type=complete_intersection
    END IF
    subdomain=0.0

    IF (control(loc_subdomain_flag)==IOS_full_field.OR.                        &
        subdomain_type==complete_intersection) THEN
      target_pnt=1
      DO src_row=halo_NS+1,blsize(2,fld_type)+halo_NS
        src_point_row_start=(src_row-1)*                                       &
            lasize(1,fld_type,arg_halo_type)                                   &
            +1+halo_EW
        DO col=1,blsize(1,fld_type)
          subdomain(target_pnt)=array(src_point_row_start+col-1)
          target_pnt=target_pnt+1
        END DO
      END DO
      IF(target_pnt-1 > SIZE(subdomain))                                       &
          CALL IOS_Ereport('pack_pp',99,                                       &
          'memory overwrite in copying to halo-free buffer')
    ELSE ! partial field that doesn"t overlap completely with this proc
      IF (subdomain_type==no_intersection) THEN
        ! zero the field for sanity
        subdomain=20.0
      ELSE IF (subdomain_type==1) THEN
        ! get the local coords in my domain to copy into
        CALL getLocalSubdomainBounds                                           &
            (model_rank,fld_type,domainBox,localBox)
        target_pnt=1
        DO src_row=localBox%s,localBox%n
          DO col=localBox%w,localBox%e
            subdomain(col+(src_row-1)*blsize(1,fld_type))=                     &
                array(target_pnt)
            target_pnt=target_pnt+1
          END DO
        END DO
      ELSE IF (subdomain_type==2) THEN

        target_pnt=1
        ! Get the westerly part of the subdomain
        ! (eastmost in the local processor domain)
        CALL getLocalSubdomainBoundsCyclic                                     &
            (model_rank,fld_type,                                              &
            domainBox,localBox_W,WesterlyZone)
        CALL getLocalSubdomainBoundsCyclic                                     &
            (model_rank,fld_type,                                              &
            domainBox,localBox_E,EasterlyZone)
        ! The north south localBox boundaries should be the same

        DO src_row=localBox_W%s,localBox_W%n
          DO col=localBox_W%w,localBox_W%e
            subdomain(col+(src_row-1)*blsize(1,fld_type))=                     &
                array(target_pnt)
            target_pnt=target_pnt+1
          END DO
          DO col=localBox_E%w,localBox_E%e
            subdomain(col+(src_row-1)*blsize(1,fld_type))=                     &
                array(target_pnt)
            target_pnt=target_pnt+1
          END DO
        END DO
      END IF
    END IF

    control(loc_dmi)=0 ! thats a real, but it does not matter
    control(loc_disk_block)=IOS_block_size
    control(loc_disk_block_start)=IOS_block_start
    CALL IOS_pack_async_buffer                                                 &
        (s,subdomain,local_len,control,                                        &
        IOS_async_control_sz_max,subdomain_type)

    DEALLOCATE(subdomain)
    DEALLOCATE(control)

    IF (lhook) CALL dr_hook('IOS_DUMP:PACK_DUMP_DATA',                         &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_dump_pack_data


  SUBROUTINE  IOS_dumpSetFieldParams(flag,packType,blockSize,location)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: flag
    INTEGER, INTENT(IN) :: packType
    INTEGER, INTENT(IN) :: blockSize
    INTEGER, INTENT(IN) :: location
    IOS_dump_pack_flag = flag
    IOS_dump_pack_type = packType
    IOS_block_size     = blockSize
    IOS_block_start    = location
  END SUBROUTINE IOS_dumpSetFieldParams

!
! Send a control message block to IOS
!
  SUBROUTINE IOS_dump_write_data( s )
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
    REAL(KIND=jprb)            :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_DUMP:DUMP_WRITE_DATA',                        &
        zhook_in,zhook_handle)

    IF (slot(s)%bufferType /= IOS_Action_StashWriteDumpData)THEN
      CALL IOS_Ereport("dump_write_data", 99 ,                                 &
          "write PP data called with wrong buffer type" )
    END IF

    IF (model_rank==0) THEN
      DO targetRank=1,IOS_tasks_per_server

        qHandle = ios_init_md(                                                 &
            slot(s)%unit,-1,                                                   &
            IOS_Action_StashWriteDumpData,                                     &
            dataSize=slot(s)%control_offset,                                   &
            rank=targetRank)

        metadata => IOS_Attach_Metadata(qHandle)
        metadata % handle    = slot(s)%handle
        sendBuffer => IOS_attach_SendBuffer(qHandle)
        CALL um_memcpy64(sendBuffer,                                           &
            slot(s)%control,                                                   &
            slot(s)%control_offset)
        CALL IOS_Send(qHandle)
      END DO
    END IF

    IF (lhook) CALL dr_hook('IOS_DUMP:DUMP_WRITE_DATA',                        &
        zhook_out,zhook_handle)
  END SUBROUTINE IOS_dump_write_data

END MODULE ios_dump
