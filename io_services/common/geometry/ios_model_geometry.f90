! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

!
! Provides data about the decomposition and dimensions of a model
!

MODULE IOS_Model_Geometry

  USE IOS_Stash_Common
  USE yomhook, ONLY: lhook, dr_hook  
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

! params/vars  for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1
    
  INTEGER, PARAMETER   :: ndims=2
  INTEGER, PARAMETER   :: xdir=1
  INTEGER, PARAMETER   :: ydir=2
  INTEGER, POINTER     :: atm_global_points(:,:)=>NULL()   
  INTEGER, POINTER     :: processor_map(:,:)=>NULL()! x,y for each atm pe
  INTEGER, POINTER     :: offset_map(:,:,:)=>NULL() ! start point for each atm pe
  INTEGER, POINTER     :: size_map(:,:,:)=>NULL()   ! local_size for each atm pe
  INTEGER              :: nfld_types
  INTEGER              :: atm_numprocs
  INTEGER              :: ios_model_atm_nprocx
  INTEGER              :: ios_model_atm_nprocy
  INTEGER              :: MaxFieldDomain ! Maximum field on any proc
  LOGICAL, POINTER     :: land_mask(:)=>NULL()
  LOGICAL, POINTER     :: local_land_mask(:)=>NULL()
  CHARACTER (LEN=132),&
      PRIVATE          :: geom_message

  PRIVATE allocate_geometry
  PRIVATE compute_offsets
  PRIVATE setMaxFieldDomain

CONTAINS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign storage for model geometry information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE allocate_geometry()    
    IMPLICIT NONE
    INTEGER  :: atm_procs

    atm_procs=ios_model_atm_nprocx*ios_model_atm_nprocy

    ALLOCATE (atm_global_points(ndims,nfld_types              ))
    ALLOCATE (    processor_map(ndims           ,0:atm_procs-1))
    ALLOCATE (         size_map(ndims,nfld_types,0:atm_procs-1))
    ALLOCATE (       offset_map(ndims,nfld_types,0:atm_procs-1))
  END SUBROUTINE allocate_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set my geometry data both locally and transmit to all IO
! servers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ios_dump_init(mask, numpoints)
    USE IOS_Client_Queue
    USE MPL
    IMPLICIT NONE
    LOGICAL, INTENT(IN)       :: mask(:)
    INTEGER, INTENT(IN)       :: numpoints
    TYPE(IOS_metadata_type), &
        POINTER               :: metadata
    INTEGER                   :: server
    INTEGER                   :: sub_task
    INTEGER                   :: error
    INTEGER                   :: qHandle
    REAL(KIND=jprb)           :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER               :: sb(:)

    IF (lhook) CALL dr_hook              &
        ('IOS_MODEL_GEOMETRY:DUMP_INIT', &
        zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      DO server=1,IOS_Server_Groups
         DO sub_task=1,IOS_tasks_per_server
            qHandle=IOS_Init_MD(-1*io_servers(server,sub_task), &
                -1,IOS_Action_DumpInitModel,                    &
                dataSize=numpoints)
            sb => IOS_Attach_SendBuffer(qHandle)
            CALL um_memcpy64(sb,mask,numpoints)
            CALL IOS_Send(qHandle)
            IF (IOS_Verbosity >= IOS_PrStatus_Diag) THEN
              WRITE(6,'(A,I4,A,I4)')'Sending landmask to server ', &
                  server,', rank ',io_servers(server,sub_task)
            END IF
         END DO
      END DO
    END IF

    IF (lhook) CALL dr_hook              &
        ('IOS_MODEL_GEOMETRY:DUMP_INIT', &
        zhook_out,zhook_handle)

  END SUBROUTINE ios_dump_init

  SUBROUTINE ios_DumpInitModel(mask, numpoints)
    IMPLICIT NONE
    INTEGER, INTENT(IN)       :: numpoints
    REAL, POINTER             :: mask(:) ! Note this is really logical data
    REAL(KIND=jprb)           :: zhook_handle
    
    IF (lhook) CALL dr_hook &
        ('IOS_MODEL_GEOMETRY:DUMPINITMODEL', &
        zhook_in,zhook_handle)

    IF (ASSOCIATED(land_mask)) THEN
      CALL IOS_Ereport('DumpInitModel',99,&
          'Land Mask already allocated! Called twice?') 
    END IF
  
    ALLOCATE (land_mask(numpoints))!Never deallocate
    land_mask=TRANSFER(mask,land_mask(1),numpoints)
    WRITE(6,'(A,I8)')'Allocated Global Land Mask on IOS Points=', &
        numpoints
    IF (lhook) CALL dr_hook &
        ('IOS_MODEL_GEOMETRY:DUMPINITMODEL', &
        zhook_out,zhook_handle)

  END SUBROUTINE ios_DumpInitModel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set my geometry data both locally and transmit to all IO
! servers, having flattened all needed data to a 1D array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_Client_Geometry_Init(deco)
    USE IOS_Client_Queue    
    USE IOS_Decompose
    USE IOS_Client_Coupler
! This is the only interfacing to the UM data structures.... hopefully. !
    USE decomp_db
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: deco
    TYPE(IOS_metadata_type), &
        POINTER             :: metadata
    INTEGER, ALLOCATABLE    :: arguments(:)! data to write to ios
    INTEGER                 :: ierror
    INTEGER                 :: numArguments
    INTEGER                 :: ierr
    INTEGER                 :: server
    INTEGER                 :: sub_task
    INTEGER                 :: qHandle
    INTEGER                 :: cpu
    INTEGER                 :: fld_type
    INTEGER                 :: dim
    INTEGER                 :: ii
    INTEGER                 :: ibuf
    INTEGER                 :: tag
    INTEGER(KIND=IOS_TUKind), &
        POINTER             :: sb(:)
    REAL(KIND=jprb)         :: zhook_handle

    IF (lhook) CALL dr_hook &
        ('IOS_MODEL_GEOMETRY:GEOMETRY_INIT', &
        zhook_in,zhook_handle)

    nfld_types=nfld_max

    IF (.NOT.L_IOS_Active()) THEN
      CALL IOS_Ereport('ios_model_geometry:geometry_init',99, &
          'IOS is not active') 
    END IF
    numArguments=3+2*(model_procs+1)*nfld_types
    ALLOCATE(arguments(numArguments))

    arguments(1:2)=decompDB(deco)%gridsize(1:2)
    arguments(3)=nfld_types

    ii=4

! Provide global field sizes
    DO fld_type=1,nfld_types
      arguments(ii:ii+1)=decompDB(deco)%glSize(1:2,fld_type)
      ii=ii+2
    END DO

! Provide local (i.e. per cpu) grid sizes
    DO fld_type=1,nfld_types
      DO cpu=0,model_procs-1
        arguments(ii:ii+1)= &
            decompDB(deco)%g_blsize(1:2,fld_type,cpu)
        ii=ii+2
      END DO
    END DO

    IF (model_rank == 0 ) THEN
      DO server=1,IOS_Server_Groups
         DO sub_task=1,IOS_tasks_per_server
            qHandle=IOS_Init_MD(-1*io_servers(server,sub_task), &
                -1,IOS_Action_StashInitModel,                   &
                dataSize=numArguments)
            sb => IOS_Attach_SendBuffer(qHandle)
            CALL um_memcpy64(sb,arguments,numArguments)
            CALL IOS_Send(qHandle)
         END DO
      END DO
    END IF

! Now initialise the geometry locally
! This data is all on the client, but we'd like IOS clients and servers 
! to use the same code for computation of quantities for seperataion
! reasons
    ios_model_atm_nprocy=decompDB(deco)%gridsize(2)
    ios_model_atm_nprocx=decompDB(deco)%gridsize(1)
    atm_numprocs=ios_model_atm_nprocx*ios_model_atm_nprocy

    CALL allocate_geometry()
    
! Populate our data structures from Decomp_DB
    atm_global_points(1:ndims,1:nfld_types)= &
        decompDB(deco)%glsize(1:ndims,1:nfld_types)
    size_map(1:ndims,1:nfld_types,0:atm_numprocs-1)= &
        decompDB(deco)%g_blsize(1:ndims,1:nfld_types,0:atm_numprocs-1)

! Compute quantities to quickly lookup out location in the grid
    CALL compute_offsets()

! Also compute the local maximum
    CALL setMaxFieldDomain()


! Allocate sendbuffers.
    DO ibuf=1,IOS_AsyncNumSlots
      ALLOCATE(slot(ibuf)%buffer (MaxFieldDomain*IOS_AsyncMaxFieldsInPack))
      ALLOCATE(slot(ibuf)%control(                                &
          ios_max_cntrlwords_per_field * &
          IOS_AsyncMaxFieldsInPack))
    END DO


    WRITE(6,'(A,I9,A)')                                      &
        'IOS_Geometry allocated ',IOS_AsyncNumSlots*         &
        MaxFieldDomain*IOS_AsyncMaxFieldsInPack,' words on', &
        ' each client for send buffers'


    target_ios_rank=locate_in_range        &
        (1,ios_model_atm_nprocy,           &
        model_rank/ios_model_atm_nprocx+1, &
        IOS_tasks_per_server)

    IF (IOS_Verbosity >= IOS_PrStatus_Oper)              &
        WRITE(6,'(A,I2)')                                &
        'IOS: Info: Asynchronous message target rank: ', &
        target_ios_rank

! This should be pretty when printed

    IF (IOS_Verbosity >= IOS_PrStatus_Diag) THEN
      WRITE(6,*)
      WRITE(6,'(A)')'^ North      IOS Geometry initisialisation'
      WRITE(6,'(A)')'|            -----------------------------'
      WRITE(6,'(A)')'|'
      WRITE(6,'(A)')'+-----> East'
      WRITE(6,'(A)')

      DO fld_type=1,nfld_types
        WRITE(6,*)
        WRITE(6,'(A,I2)')'Field Type ',fld_type
        WRITE(6,'(A,I5)')'  CPU =',model_rank 
        WRITE(6,'(A,I5,A,I5)')'  Grid=',        &
            atm_global_points(1,fld_type),'x',  &
            atm_global_points(2,fld_type)
        WRITE(6,'(A,I5,A,I5)')' LSize=',        &
            size_map(1,fld_type,model_rank),'x',      &
            size_map(2,fld_type,model_rank)     
        WRITE(6,'(A,i5,A,i5)')                      &
            '        ',offset_map(1,fld_type,model_rank), &
            '       ',offset_map(1,fld_type,model_rank)+  &
            size_map(1,fld_type,model_rank)-1
        WRITE(6,'(A)')'         |          |         '
        WRITE(6,'(A,I5)')'---------+----------+---------', &
            offset_map(2,fld_type,model_rank)+size_map(2,fld_type,model_rank)-1
        WRITE(6,'(A)')'         |          |         '
        WRITE(6,'(A)')'         |          |         '
        WRITE(6,'(A,I5,A)')'         |    ',model_rank,' |         '
        WRITE(6,'(A)')'         |          |         '
        WRITE(6,'(A)')'         |          |         '
        WRITE(6,'(A,I5)')'---------+----------+---------',&
            offset_map(2,fld_type,model_rank)
        WRITE(6,'(A)')'         |          |         '
        WRITE(6,*)
      END DO

    END IF

    DEALLOCATE(arguments)
    IF (lhook) CALL dr_hook                  &
        ('IOS_MODEL_GEOMETRY:GEOMETRY_INIT', &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_Client_Geometry_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the first global coordinate in each subdomain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE compute_offsets()
    IMPLICIT NONE
    INTEGER :: iy
    INTEGER :: ix
    INTEGER :: ifld
    INTEGER :: p
    INTEGER :: p2

    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook                    &
        ('IOS_MODEL_GEOMETRY:COMPUTE_OFFSETS', &
        zhook_in,zhook_handle)

    !now create the offset map X 
    !1st column is all ones (no offset)
    DO iy=1,ios_model_atm_nprocy
       p=processor_from_XY(1,iy)
       DO ifld=1,nfld_types
          offset_map(1,ifld,p)=1
       END DO
    END DO
    
    !then just add the last offset to the previous data size
    DO ix=2,ios_model_atm_nprocx
       DO iy=1,ios_model_atm_nprocy
          p=processor_from_XY(ix,iy)
          p2=processor_from_XY(ix-1,iy)
          DO ifld=1,nfld_types
             offset_map(1,ifld,p)= &
                 offset_map(1,ifld,p2)+size_map(1,ifld,p2)
          END DO
       END DO
    END DO
    
    !now create the offset map Y
    !1st row is all ones (no offset)
    DO ix=1,ios_model_atm_nprocx
       p=processor_from_XY(ix,1)
       DO ifld=1,nfld_types
          offset_map(2,ifld,p)=1
       END DO
    END DO

    !then just add the last offset to the previous data size
    DO iy=2,ios_model_atm_nprocy
       DO ix=1,ios_model_atm_nprocx
          p=processor_from_XY(ix,iy)
          p2=processor_from_XY(ix,iy-1)
          DO ifld=1,nfld_types
             offset_map(2,ifld,p)= &
                 offset_map(2,ifld,p2)+size_map(2,ifld,p2)
          END DO
       END DO
    END DO

    IF (lhook) CALL dr_hook &
        ('IOS_MODEL_GEOMETRY:COMPUTE_OFFSETS', &
        zhook_out,zhook_handle)
    
  END SUBROUTINE compute_offsets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return an atm cpu number (0..N-1) from X,Y address (1,N) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER FUNCTION processor_from_XY(x,y)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: x,y
    processor_from_XY=(y-1)*ios_model_atm_nprocx+(x-1)
  END FUNCTION processor_from_XY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return an X,Y address (1..N) from an atm cpu number (0..N-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE XY_from_processor(x,y,proc)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: proc
    INTEGER,INTENT(OUT):: x,y

    INTEGER            :: fproc
    fproc=proc+1
    y=proc/ios_model_atm_nprocx+1
    x=fproc-ios_model_atm_nprocx*(y-1)
  END SUBROUTINE XY_from_processor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialise geometry on IO servers. Only when each atm cpu has 
! sent data do we calculate our offset tables etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_Server_Geometry_Init(setupData)
    IMPLICIT NONE
    INTEGER :: setupData(:)
    INTEGER :: client
    INTEGER :: fld_type
    INTEGER :: ii
    INTEGER :: x,y           ! Sanity checks
    INTEGER :: newproc,proc  ! Sanity checks
    REAL(KIND=jprb)           :: zhook_handle
    LOGICAL :: first_call=.TRUE.
    IF (lhook) CALL dr_hook &
        ('IOS_MODEL_GEOMETRY:SET_SERVER', &
        zhook_in,zhook_handle)

    IF (SIZE(setupData)<3)                                  &
        CALL IOS_Ereport('ios_geometry_set_server',99, &
        'Data buffer unexpectedly small') 
    
    ios_model_atm_nprocx = setupData(1)
    ios_model_atm_nprocy = setupData(2)
    nfld_types           = setupData(3)

    IF (SIZE(setupData) /= nfld_types*2*(atm_numprocs+1)+3) THEN
      WRITE(geom_message,*)'Data buffer wrong size: Expected ',    &
          nfld_types*2*(atm_numprocs+1)+3,' items, but received ', &
          SIZE(setupData)
      CALL IOS_Ereport('ios_geometry_set_server',99, &
          geom_message) 
    END IF

    IF (first_call) THEN
      CALL allocate_geometry()
      first_call=.FALSE.
    ELSE
        CALL IOS_Ereport('ios_geometry_set_server',99, &
        'Function was wrongly called twice.')       
    END IF

    ii=4
    DO fld_type=1,nfld_types
      atm_global_points(1:2,fld_type) = setupData(ii:ii+1)
      ii=ii+2
    END DO

    DO fld_type=1,nfld_types
      DO client=0,atm_numprocs-1
        size_map(1:2,fld_type,client) = setupData(ii:ii+1)
        ii=ii+2
      END DO
    END DO

    ! Sanity check the lookup tables
    DO proc=0,atm_numprocs-1
      CALL XY_from_processor(x,y,proc) 
      newproc=processor_from_XY(x,y) 
      IF (newproc /= proc)THEN
        WRITE(geom_message,*)'table error: reverse lookup gave cpu '&
            ,newproc
        CALL IOS_Ereport( 'geometry_set_server', 99, geom_message )
      END IF
    END DO
    
    CALL compute_offsets()
    
    CALL setMaxFieldDomain()
    
    IF (IOS_Verbosity >= IOS_PrStatus_Diag) THEN
      
      WRITE(6,*) 
      WRITE(6,'(A,A)') &
          'Processor geometry; ', &
          'g=global, l=local, off=offset in global'
      WRITE(6,*) 
      WRITE(6,'(A3,1x,A3,1x,A4,1x,A4,1x,A4,1x,A4,1x,A4,1x,A4)') &
          'CPU','fld','gEW','gNS','lEW','lNS','offE','offN'
      WRITE(6,'(A)')'---------------------------------------'
      DO fld_type=1,nfld_types
        DO proc=0,atm_numprocs-1
          WRITE(6,'(I3,1x,I3,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4)') &
              proc,fld_type,                 &
              atm_global_points(1,fld_type), &
              atm_global_points(2,fld_type), &            
              size_map(1,fld_type,proc),     &
              size_map(2,fld_type,proc),     &
              offset_map(1,fld_type,proc),   &
              offset_map(2,fld_type,proc)
        END DO
      END DO
      WRITE(6,'(A)')'---------------------------------------'
      WRITE(6,*) 
    END IF ! Printing out geometry
    
    
    IF (lhook) CALL dr_hook &
        ('IOS_MODEL_GEOMETRY:SET_SERVER', &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_Server_Geometry_Init

  SUBROUTINE setMaxFieldDomain()
    IMPLICIT NONE
    INTEGER :: ifld,mfld
    INTEGER :: icpu,mcpu
    REAL(KIND=jprb)           :: zhook_handle

    IF (lhook) CALL dr_hook &
        ('IOS_MODEL_GEOMETRY:SETMAXFIELDDOMAIN', &
        zhook_in,zhook_handle)

    MaxFieldDomain=0
    
    DO ifld=1,nfld_types
      DO icpu=0,ios_model_atm_nprocx*ios_model_atm_nprocy-1
         IF ( size_map(1,ifld,icpu)* &
              size_map(2,ifld,icpu) > MaxFieldDomain) THEN
            MaxFieldDomain=size_map(1,ifld,icpu)*size_map(2,ifld,icpu)
            mfld=ifld
            mcpu=icpu
         END IF
      END DO
    END DO

    WRITE(6,'(A,I7,A)')'IOS_Model_Geometry: Info: Max Field Domain=', &
        MaxFieldDomain,' words'
    WRITE(6,'(A,I6,A,I2)')'IOS_Model_Geometry: Info: for CPU=',mcpu,  &
         ' Field type=',mfld
    
    IF (lhook) CALL dr_hook &
        ('IOS_MODEL_GEOMETRY:SETMAXFIELDDOMAIN', &
        zhook_out,zhook_handle)

  END SUBROUTINE setMaxFieldDomain

  INTEGER FUNCTION getMaxFieldDomain()
    IMPLICIT NONE
    getMaxFieldDomain=MaxFieldDomain
  END FUNCTION getMaxFieldDomain

END MODULE IOS_MODEL_GEOMETRY
