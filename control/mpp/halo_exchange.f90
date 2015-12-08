! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

!
! This module provides methods for updating boundaries (halos)
! on fields, either by extending the field or exchanging data
! with neighbouring processors depending on the model.
!
! Communications take place in private communicators which musch be 
! initialised with a call to halo_exchange_init() before other methods
! can be used. 
!
! Logic and communicator setup is non-trivial for global models.

MODULE Halo_Exchange
  
  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE mpl, ONLY :      &
      mpl_real,        &
      mpl_status_size

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE UM_ParVars
  USE domain_params
  IMPLICIT NONE

  INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_out = 1    

  TYPE domain
    INTEGER :: x ! aka row_length
    INTEGER :: y ! aka rows
    INTEGER :: hx ! aka halo_x
    INTEGER :: hy ! aka halo_y      
    INTEGER :: levels
  END TYPE domain

! Communicators to use for comms
  INTEGER, PRIVATE :: global_comm
  INTEGER, PRIVATE :: ew_comm
  INTEGER, PRIVATE :: ns_comm


! grid coordinates in the model
  INTEGER          :: my_proc_row
  INTEGER          :: my_proc_col

! Rank in the row/col communicators
  INTEGER          :: mype_ew
  INTEGER          :: mype_ns

! Processor index !!! in the row/column communicator !!eleventy1!!
! or nodomain if not sensible
  INTEGER          :: he_neighbour(4)
! Tags
  INTEGER, PARAMETER :: base_tag_west_to_east=1
  INTEGER, PARAMETER :: base_tag_east_to_west=2
  INTEGER, PARAMETER :: base_tag_north_to_south=3
  INTEGER, PARAMETER :: base_tag_south_to_north=4

  LOGICAL       :: HE_Initialised=.FALSE.

! Usual convention
!  increasing x=i=more easterly = fast moving
!  increasing y=j=more northerly = slow moving

CONTAINS

! Purpose:
! Initialise the module, must be called before other routines
! Set up communicators for halo exchange ops
  SUBROUTINE Halo_Exchange_Init()
    IMPLICIT NONE
    INTEGER             :: fld
    INTEGER             :: he_rows_max
    INTEGER             :: he_row_len_max
    INTEGER             :: he_halo_x_max
    INTEGER             :: he_halo_y_max
    INTEGER             :: i,ierror,row_start
    INTEGER             :: ewsize,nssize,ns_comm_size
    INTEGER             :: dummy_comm
    INTEGER             :: ns_index_delta
    INTEGER             :: key,colour,rr,cr
    INTEGER             :: proc_id_north
    INTEGER             :: proc_id_south
    INTEGER             :: proc_id_east
    INTEGER             :: proc_id_west
    INTEGER, PARAMETER  :: model_levels=100 ! temporary
    LOGICAL             :: isOverpole
    REAL(KIND=jprb)     :: zhook_handle



    IF (lhook) CALL dr_hook('HALO_EXCHANGE:HALO_EXCHANGE_INIT', &
        zhook_in,zhook_handle)
    
    CALL gc_get_communicator(global_comm,ierror)
    row_start      = (mype/nproc_x) * nproc_x
    my_proc_col    = mype - row_start
    my_proc_row    = (mype/nproc_x)
    
    ! Define a row communicator in natural order, moving from west (0)
    ! to east.
    colour         =my_proc_row
    key            =my_proc_col
    CALL MPL_Comm_Split(global_comm, colour, key, ew_comm, ierror)


    ! Define a column communicator in natural order, moving from south (0), in
    ! an initially northerly direction.

    ! For the global model a column is more complex. A column is a complete
    ! south->north->south circumnavigation. Columns start at the south pole
    ! in the west hemisphere traverse the north pole and end at the south 
    ! pole in the east hemisphere, thus the set of 'columns' forms a set of 
    ! nproc_x/2 rings.
    !
    ! LAM and Torus2D global models form more conventional columns without/with
    ! wrapping boundary conditions respectively. 

    colour         =my_proc_col
    key            =my_proc_row
    isOverpole     =.FALSE.  
    IF (nproc_x>1) THEN  ! it is an even number by definition
      IF (sb_model_domain==mt_global.AND.bound(2) /= bc_cyclic)THEN
        IF (my_proc_col>=nproc_x/2) THEN 
          ! Easter hemisphere processors want to join the communicator 
          ! that loops over the north pole from the 
          ! western hemisphere
          colour=my_proc_col-nproc_x/2
          ! and I want to index such that the +1
          ! direction is south (the continued direction of north
          ! as you move over the pole)
          key   =2*nproc_y-my_proc_row-1
          isOverPole=.TRUE.
        END IF
      END IF
    END IF
    
    CALL MPL_Comm_Split(global_comm, colour, key, ns_comm, ierror)    
    CALL MPL_Comm_Size(ns_comm ,ns_comm_size  ,ierror)
    CALL MPL_Comm_Rank(ew_comm ,mype_ew  ,ierror)
    CALL MPL_Comm_Rank(ns_comm ,mype_ns  ,ierror)
    
    ! Set up convenience indexing of processor neighbours    
    proc_id_north  =nodomain
    proc_id_south  =nodomain
    proc_id_east   =nodomain
    proc_id_west   =nodomain
    
    ! EW neighbours are trivial due to the natural ordering of
    ! the cummunicator
    IF (neighbour(peast) /= nodomain) THEN      
      proc_id_east=my_proc_col+1
      IF (proc_id_east==nproc_x) proc_id_east=0
    END IF
    
    IF (neighbour(pwest) /= nodomain) THEN
      proc_id_west=my_proc_col-1
      IF (proc_id_west==-1) proc_id_west=nproc_x-1
    END IF
    
    ! For NS the direction of northerly flips at the pole
    ! so we change our 'increment' to be negative.
    ns_index_delta=1                 ! normally north is +1
    IF (isOverPole)ns_index_delta=-1 ! except in the e. hemisphere
    !                                ! in global models
    
    ! If there is only me, then my neighbour must be me, also.
    IF (nproc_y==1 .AND. nproc_x==1)ns_index_delta=0  
    
    IF (neighbour(pnorth) /= nodomain) THEN
       proc_id_north=key+ns_index_delta      
      
       ! If we are at maximal processor in the south->north 
       ! direction, then we are either 
       ! - a global model, and hence at the southerly pole 
       !   in east hemisphere, so traversing the south pole
       !   takes us back to rank 0 (southerly in the western 
       !   hemisphere)
       ! - a cyclic model on the northerly edge, and north
       !   wraps around to the south edge, and hence back 
       !   to rank 0
       IF (proc_id_north==ns_comm_size) proc_id_north=0
      
    END IF
    
    IF (neighbour(psouth) /= nodomain) THEN
      proc_id_south=key-ns_index_delta
            
      ! adjust s. hem southerly row south to max in communicator
      ! wrapping under the pole for global etc etc
      IF (proc_id_south==-1)           proc_id_south=ns_comm_size-1
      IF (proc_id_south==ns_comm_size) proc_id_south=0
    END IF

    ! Store these in the module.
    he_neighbour(pEast)=proc_id_east
    he_neighbour(pWest)=proc_id_west
    he_neighbour(pNorth)=proc_id_north
    he_neighbour(pSouth)=proc_id_south

    he_rows_max=1
    he_row_len_max=1
    he_halo_x_max=1
    he_halo_y_max=1
    DO i=1,Nfld_max
      IF (blsize(1,i)>he_row_len_max)  he_row_len_max=blsize(1,i)
      IF (blsize(2,i)>he_rows_max)     he_rows_max=blsize(2,i)
    END DO

    DO i=1,Nhalo_max
      IF (halosize(1,i)>he_halo_x_max)he_halo_x_max=halosize(1,i)
      IF (halosize(2,i)>he_halo_y_max)he_halo_y_max=halosize(2,i)
    END DO

    ewsize=model_levels* &
        ((he_rows_max*he_halo_y_max)    +(2*he_halo_x_max*he_halo_y_max))
    nssize=model_levels* &
        ((he_row_len_max*he_halo_y_max) +(2*he_halo_x_max*he_halo_y_max))

    HE_Initialised=.TRUE.

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:HALO_EXCHANGE_INIT', &
        zhook_out,zhook_handle)
    
  END SUBROUTINE Halo_Exchange_Init

! Purpose:
! Update halos of the supplied field
!
! This version is used if the module is in use, otherwise 
! the f77 interface will be used.
  SUBROUTINE swap_bounds(                &
      field, row_length, rows, levels,   & ! field
      halo_x, halo_y,                    & ! halos
      field_type, l_vector               & ! supporting information
      )
    
    IMPLICIT NONE  
      
    INTEGER, INTENT(IN) ::  row_length ! number of points on a row
    !                                  ! (not including halos)
    INTEGER, INTENT(IN) ::  rows       ! number of rows in a theta field
    !                                  ! (not including halos)
    INTEGER, INTENT(IN) ::  levels     ! number of model levels
    INTEGER, INTENT(IN) ::  halo_x     ! size of halo in "i" direction
    INTEGER, INTENT(IN) ::  halo_y     ! size of halo in "j" direction
    INTEGER, INTENT(IN) ::  field_type ! Defines the grid interpolation type
    !                                  ! of the input FIELD (u,v or w)
    LOGICAL, INTENT(IN) :: l_vector    ! TRUE:  Data is a horizontal vector
    !                                  !     component
    !                                  !     FALSE: Data is a scalar
    
    REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
        1-halo_y:rows+halo_y,levels)   ! Field to have its halos updated
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS', &
        zhook_in,zhook_handle)
  
    !------------------------------------------------------------------
    ! 2.0 Exchanges in EW then NS directions
    IF (halo_x > 0)                        &
        
        CALL swap_bounds_EW(               &
        field, row_length, rows, levels,   & ! field
        halo_x, halo_y                     & ! halos
        )
  
    IF (halo_y > 0)                        &
        CALL swap_bounds_NS(               &
        field, row_length, rows, levels,   & ! field
        halo_x, halo_y,                    & ! halos
        field_type, l_vector               & ! supporting information
        )

    IF (sb_model_domain == mt_lam) THEN
!DEPENDS ON: fill_external_halos
      CALL fill_external_halos(field, row_length, rows, levels, &
                               halo_x, halo_y)
    END IF
    
    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS', &
        zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE swap_bounds


! Purpose:
!   This subroutine takes care of all boundary swapping in EW and
!   extending of arrays at the global boundaries in EW.
  SUBROUTINE swap_bounds_ew(             &
      field, row_length, rows, levels,   & ! field
      halo_x, halo_y                     & ! halos
      )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: row_length ! number of points on a row
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: levels     ! number of model levels
    INTEGER, INTENT(IN) :: halo_x     ! size of halo in "i" direction
    INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction  
    REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
        1-halo_y:rows+halo_y,levels)  !Field to have its halos updated

    ! Comdecks
    
    ! Local variables
    
    ! The send and recieve buffers, replacing common blocks
    REAL ::                                                       &
        send_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y))),   &
        receive_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y)))

    INTEGER  ::                  &
        full_rows,               & ! number of rows including halo rows
        ew_halo_size,            & ! size of EW halo region
        index_2_start              ! start address of index2 in buffers


    ! Variables used by MPL/GCOM
    INTEGER  :: ierror
    INTEGER  :: istat(mpl_status_size,4)
    INTEGER  :: ireq(4)
    INTEGER  :: nreq
    INTEGER  :: nr
    

    INTEGER ::  i,j,k              ! Spatial loop counters

    ! Statement functions for addressing the buffer arrays

    INTEGER :: ew_address

    ! Variables used in statment functions

    INTEGER :: row,point,halo,level,indy

    REAL(KIND=jprb)               :: zhook_handle

    ew_address(row,halo,level,indy)=                                    &
        (indy-1)*(levels*ew_halo_size) +                                &
        (level-1)*ew_halo_size +                                        &
        (halo-1)*full_rows +                                            &
        row+halo_y

    !------------------------------------------------------------------
    ! 0.0 Check if there is anything to do

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_EW', &
        zhook_in,zhook_handle)

    !------------------------------------------------------------------
    ! 0.1 Select an optimised kernel if available
    
     IF  (nproc_x == 1) THEN
       CALL swap_bounds_ew_serial(          &
           field, row_length, rows, levels, & ! field
           halo_x,halo_y                    & ! halos
           )
     ELSE IF (halo_x == 1) THEN
       CALL swap_bounds_ew_h1(              &
           field, row_length, rows, levels, & ! field
           halo_y                           & ! halos
           ) 
     ELSE
      
      !------------------------------------------------------------------
      ! 1.0 Initialise variables
      
      full_rows       = rows + 2*halo_y      
      ew_halo_size=full_rows*halo_x

      !------------------------------------------------------------------
      ! 2.0 East-West communications

      ! index_2_start points to the start of the second index
      ! within the buffer arrays. The first index contains
      ! data sent from the Western side, the second index
      ! contains data sent from the Eastern side

      nreq = 0
      index_2_start=ew_address(1-halo_y,1,1,2)

      !---------------------------------------
      ! 2.1.1 Post receives


      IF (he_neighbour(pwest)  /=  nodomain) THEN
        ! Receive West
        nreq=nreq+1
        CALL mpl_irecv(receive_buffer(index_2_start),      &
            ew_halo_size*levels,                           &
            mpl_real,he_neighbour(pwest),502,              &
            ew_comm,ireq(nreq),ierror)
      END IF

      IF (he_neighbour(peast)  /=  nodomain) THEN
        ! Receive East
        nreq=nreq+1
        CALL mpl_irecv(receive_buffer,ew_halo_size*levels, &
            mpl_real,he_neighbour(peast),501,              &
            ew_comm,ireq(nreq),ierror)
      END IF

      !---------------------------------------
      ! 2.1.2 Copy the data into buf_send
      
      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          DO i=1,halo_x

            ! Copy stuff from the Western side of the grid
            send_buffer(ew_address(j,i,k,1))=field(i,j,k)

            ! Copy stuff from the Eastern side of the grid
            send_buffer(ew_address(j,i,k,2))=                         &
                field(row_length-halo_x+i,j,k)

          END DO ! I
        END DO ! J
      END DO ! K

      !---------------------------------------
      ! 2.1.3 Send and the data and wait for completion

      IF (he_neighbour(pwest)  /=  nodomain) THEN
        ! Send West
        
        nreq=nreq+1
        CALL mpl_isend(send_buffer,ew_halo_size*levels,  &
            mpl_real,he_neighbour(pwest),501,            &
            ew_comm,ireq(nreq),ierror)
        
      END IF
      
      IF (he_neighbour(peast)  /=  nodomain) THEN
        ! Send East
        nreq=nreq+1
        CALL mpl_isend(send_buffer(index_2_start),       &
            ew_halo_size*levels,                         &
            mpl_real,he_neighbour(peast),502,            &
            ew_comm,ireq(nreq),ierror)        
      END IF

      CALL mpl_waitall(nreq,ireq,istat,ierror)

      !---------------------------------------
      ! 2.1.4 Fill the halos with data

      IF (he_neighbour(peast)  /=  nodomain) THEN

        ! unpack data from receive_buffer into FIELD

        DO k=1,levels
          DO j=1-halo_y,rows+halo_y
            DO i=1,halo_x
              field(row_length+i,j,k)=                                &
                  receive_buffer(ew_address(j,i,k,1))
            END DO
          END DO
        END DO

      ELSE IF (sb_model_domain  ==  mt_global) THEN
        ! No neighbour to my East (ie. at edge of the domain
        ! and it's not a cyclic boundary condition)
        ! Just copy data from last column
        ! NOTE: This block of code should never be executed. It's being left
        !       in so that it can be used in the future by changing the logic
        !       in the ELSE statement above.

        DO k=1,levels
          DO j=1-halo_y,rows+halo_y
            DO i=1,halo_x
              field(row_length+i,j,k)=field(row_length,j,k)
            END DO
          END DO
        END DO

      END IF ! IF (he_neighbour(PEast)  /=  NoDomain)

      IF (he_neighbour(pwest)  /=  nodomain) THEN

        ! unpack data from receive_buffer into FIELD

        DO k=1,levels
          DO j=1-halo_y,rows+halo_y
            DO i=1,halo_x
              field(i-halo_x,j,k)=                                    &
                  receive_buffer(ew_address(j,i,k,2))
            END DO
          END DO
        END DO

      ELSE IF (sb_model_domain  ==  mt_global) THEN
        ! No he_neighbour to my West (ie. at edge of the domain
        ! and it's not a cyclic boundary condition)
        ! Just copy data from first column
        ! NOTE: This block of code should never be executed. It's being left
        !       in so that it can be used in the future by changing the logic
        !       in the ELSE statement above.

        DO k=1,levels
          DO j=1-halo_y,rows+halo_y
            DO i=1,halo_x
              field(i-halo_x,j,k)=field(1,j,k)
            END DO
          END DO
        END DO

      END IF ! IF (he_neighbour(PWest)  /=  NoDomain)        

    END IF ! Optimised Kernel Switch
    
    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_EW', &
        zhook_out,zhook_handle)
    RETURN
    
  END SUBROUTINE swap_bounds_ew
  
! Purpose:
!   This subroutine takes care of all boundary swapping in EW and
!   extending of arrays at the global boundaries in EW.
!   special case only for halo_x=1
  SUBROUTINE swap_bounds_ew_h1(             &
      field, row_length, rows, levels,   & ! field
      halo_y                             & ! halos
      )
    
    IMPLICIT NONE
    INTEGER, PARAMETER  :: halo_x=1   ! size of halo in "i" direction

    INTEGER, INTENT(IN) :: row_length ! number of points on a row
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: levels     ! number of model levels
    INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction  
    REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
        1-halo_y:rows+halo_y,levels)  !Field to have its halos updated

    ! Comdecks
    
    ! Local variables
    
    ! The send and recieve buffers, replacing common blocks
    REAL ::                                                       &
        send_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y))),   &
        receive_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y)))

    INTEGER  ::                  &
        full_row_length,         & ! length of row including halos
        full_rows,               & ! number of rows including halo rows
        ew_halo_size,            & ! size of EW halo region
        ns_halo_size,            & ! size of NS halo region
        west_halo_source,        & ! first column of data for west halo
        east_halo_source,        & ! first column of data for east halo
        south_halo_source,       &
        north_halo_source,       &
        half_full_row_length,    & ! half of the total EW dimension
        half_row_length,         & ! half of the data (no halo) EW dim
        index_2_start,           & ! start address of index2 in buffers
        north_off,               & ! Offsets to use when copying data
        south_off                  ! to send around poles


    ! Variables used by MPL/GCOM
    INTEGER  :: ierror
    INTEGER  :: istat(mpl_status_size,4)
    INTEGER  :: ireq(4)
    INTEGER  :: nreq
    INTEGER  :: nr

    INTEGER  :: i,j,k              ! Spatial loop counters
       
    ! Statement functions for addressing the buffer arrays
    
    INTEGER :: ew_address
    
    ! Variables used in statment functions
    
    INTEGER :: row,point,halo,level,indy
    
    REAL(KIND=jprb)               :: zhook_handle
    
    ew_address(row,halo,level,indy)=                                    &
        (indy-1)*(levels*ew_halo_size) +                                &
        (level-1)*ew_halo_size +                                        &
        (halo-1)*full_rows +                                            &
        row+halo_y
    
    !------------------------------------------------------------------
    ! 0.0 Check if there is anything to do
    
    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_EW_H1', &
        zhook_in,zhook_handle)

    !------------------------------------------------------------------
    ! 1.0 Initialise variables
      
    full_rows       = rows + 2*halo_y    
    ew_halo_size    =full_rows
    
    !------------------------------------------------------------------
    ! 2.0 East-West communications
    
    ! index_2_start points to the start of the second index
    ! within the buffer arrays. The first index contains
    ! data sent from the Western side, the second index
    ! contains data sent from the Eastern side
    
    nreq = 0
    index_2_start=ew_address(1-halo_y,1,1,2)
    
    !---------------------------------------
    ! 2.1.1 Receive the data
    
    
    IF (he_neighbour(pwest)  /=  nodomain) THEN      
      ! Receive West
      nreq=nreq+1
      CALL mpl_irecv(receive_buffer(index_2_start),      &
          ew_halo_size*levels,                           &
          mpl_real,he_neighbour(pwest),502,              &
          ew_comm,ireq(nreq),ierror)
    END IF

    IF (he_neighbour(peast)  /=  nodomain) THEN
      ! Receive East
      nreq=nreq+1
      CALL mpl_irecv(receive_buffer,ew_halo_size*levels, &
          mpl_real,he_neighbour(peast),501,              &
          ew_comm,ireq(nreq),ierror)
    END IF

    !---------------------------------------
    ! 2.1.2 Copy the data into buf_send
    
    DO k=1,levels
      DO j=1-halo_y,rows+halo_y
        ! Copy stuff from the Western side of the grid
        send_buffer(ew_address(j,1,k,1))=field(1,j,k)
        ! Copy stuff from the Eastern side of the grid
        send_buffer(ew_address(j,1,k,2))=field(row_length,j,k)
      END DO ! J
    END DO ! K

    !---------------------------------------
    ! 2.1.3 Send data and wait    
    
    IF (he_neighbour(pwest)  /=  nodomain) THEN      
      ! Send West
      nreq=nreq+1
      CALL mpl_isend(send_buffer,ew_halo_size*levels,     &
          mpl_real,he_neighbour(pwest),501,               &
          ew_comm,ireq(nreq),ierror)      
    END IF

    IF (he_neighbour(peast)  /=  nodomain) THEN
      ! Send East
      nreq=nreq+1
      CALL mpl_isend(send_buffer(index_2_start),          &
          ew_halo_size*levels,                            &
          mpl_real,he_neighbour(peast),502,               &
          ew_comm,ireq(nreq),ierror)      
    END IF

    CALL mpl_waitall(nreq,ireq,istat,ierror)
    
    !---------------------------------------
    ! 2.1.4 Fill the halos with data
    
    IF (he_neighbour(peast)  /=  nodomain) THEN
      
      ! unpack data from receive_buffer into FIELD
      
      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          field(row_length+1,j,k)=                                &
              receive_buffer(ew_address(j,1,k,1))
        END DO
      END DO

    ELSE IF (sb_model_domain  ==  mt_global) THEN
      ! No neighbour to my East (ie. at edge of the domain
      ! and it's not a cyclic boundary condition)
      ! Just copy data from last column
      ! NOTE: This block of code should never be executed. It's being left
      !       in so that it can be used in the future by changing the logic
      !       in the ELSE statement above.

      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          field(row_length+1,j,k)=field(row_length,j,k)
        END DO
      END DO
      
    END IF ! IF (he_neighbour(PEast)  /=  NoDomain)

    IF (he_neighbour(pwest)  /=  nodomain) THEN
      
      ! unpack data from receive_buffer into FIELD
      
      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          field(0,j,k)=receive_buffer(ew_address(j,1,k,2))
        END DO
      END DO
      
    ELSE IF (sb_model_domain  ==  mt_global) THEN
      ! No neighbour to my West (ie. at edge of the domain
      ! and it's not a cyclic boundary condition)
      ! Just copy data from first column
      ! NOTE: This block of code should never be executed. It's being left
      !       in so that it can be used in the future by changing the logic
      !       in the ELSE statement above.
      
      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          field(0,j,k)=field(1,j,k)
        END DO
      END DO
      
    END IF ! IF (he_neighbour(PWest)  /=  NoDomain)

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_EW_H1', &
        zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE swap_bounds_ew_h1

! Purpose:
!   This subroutine takes care of all boundary swapping in EW and
!   extending of arrays at the global boundaries in EW.
!   special case only for nproc_x=1
  SUBROUTINE swap_bounds_ew_serial(      &
      field, row_length, rows, levels,   & ! field
      halo_x, halo_y                     & ! halos
      )
    
    IMPLICIT NONE
        
    INTEGER, INTENT(IN) :: row_length ! number of points on a row
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: levels     ! number of model levels
    INTEGER, INTENT(IN) :: halo_x     ! size of halo in "i" direction
    INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction  
    REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
        1-halo_y:rows+halo_y,levels)  !Field to have its halos updated

    ! Comdecks
    
    ! Local variables
    
    ! The send and recieve buffers, replacing common blocks
    REAL ::                                                      &
        send_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y))),   &
        receive_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y)))

    INTEGER  ::                  &
        west_halo_source,        & ! first column of data for west halo
        east_halo_source           ! first column of data for east halo

    INTEGER :: i,j,k               ! Spatial loop counters

    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_EW_SERIAL', &
        zhook_in,zhook_handle)

    !------------------------------------------------------------------
    ! 2.0 East-West communications on 1 processor
    
    IF(bound(1)  ==  bc_cyclic) THEN

      west_halo_source=row_length-halo_x+1 ! copy from opposite end
      east_halo_source=1                   ! of each row

!CDIR NOVECTOR
      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          DO i=1,halo_x

            ! Fill Western halo
            field(i-halo_x,j,k)=field(west_halo_source+i-1,j,k)

            ! Fill Eastern halo
            field(row_length+i,j,k)=field(east_halo_source+i-1,j,k)

          END DO ! J
        END DO ! K
      END DO ! I
    END IF    !if cyclic

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_EW_SERIAL', &
        zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE swap_bounds_ew_serial


! Purpose:
!   This subroutine takes care of boundary swapping from East
!   extending of arrays at the global boundaries easterly.
  SUBROUTINE swap_bounds_e(              &
      field, row_length, rows, levels,   & ! field
      halo_x, halo_y                     & ! halos
      )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: row_length ! number of points on a row
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: levels     ! number of model levels
    INTEGER, INTENT(IN) :: halo_x     ! size of halo in "i" direction
    INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction  
    REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
        1-halo_y:rows+halo_y,levels)  !Field to have its halos updated

    ! Comdecks
    
    ! Local variables
    
    ! The send and recieve buffers, replacing common blocks
    REAL ::                                                       &
        send_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y))),  &
        receive_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y)))

    INTEGER  ::                  &
        full_rows,               & ! number of rows including halo rows
        ew_halo_size               ! size of EW halo region


    ! Variables used by MPL/GCOM
    INTEGER  :: ierror
    INTEGER  :: istat(mpl_status_size,4)
    INTEGER  :: ireq(4)
    INTEGER  :: nreq
    INTEGER  :: nr    
    INTEGER ::  i,j,k              ! Spatial loop counters

    ! Statement functions for addressing the buffer arrays

    INTEGER :: ew_address

    ! Variables used in statment functions

    INTEGER :: row,point,halo,level,indy

    REAL(KIND=jprb)               :: zhook_handle

    ew_address(row,halo,level,indy)=                                    &
        (indy-1)*(levels*ew_halo_size) +                                &
        (level-1)*ew_halo_size +                                        &
        (halo-1)*full_rows +                                            &
        row+halo_y

    !------------------------------------------------------------------
    ! 0.0 Check if there is anything to do

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_E', &
        zhook_in,zhook_handle)


    !------------------------------------------------------------------
    ! 1.0 Initialise variables

    full_rows       = rows + 2*halo_y

    ew_halo_size=full_rows*halo_x

    !------------------------------------------------------------------
    ! 2.0 East-West communications
    
    nreq = 0

    !---------------------------------------
    ! 2.1.1 post receives

    IF (he_neighbour(peast)  /=  nodomain) THEN
      ! Receive East
      nreq=nreq+1
      CALL mpl_irecv(receive_buffer,ew_halo_size*levels,          &
          mpl_real,he_neighbour(peast),501,               &
          ew_comm,ireq(nreq),ierror)
    END IF
    !---------------------------------------
    ! 2.1.2 Copy the data into buf_send
    
    DO k=1,levels
      DO j=1-halo_y,rows+halo_y
        DO i=1,halo_x

          ! Copy stuff from the Western side of the grid
          send_buffer(ew_address(j,i,k,1))=field(i,j,k)

        END DO ! I
      END DO ! J
    END DO ! K

    !---------------------------------------
    ! 2.1.3 Send the data and wait for completion

    IF (he_neighbour(pwest)  /=  nodomain) THEN

      ! Send West
      nreq=nreq+1
      CALL mpl_isend(send_buffer,ew_halo_size*levels,             &
          mpl_real,he_neighbour(pwest),501,               &
          ew_comm,ireq(nreq),ierror)

    END IF

    CALL mpl_waitall(nreq,ireq,istat,ierror)

    !---------------------------------------
    ! 2.1.4 Fill the halos with data

    IF (he_neighbour(peast)  /=  nodomain) THEN

      ! unpack data from receive_buffer into FIELD

      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          DO i=1,halo_x
            field(row_length+i,j,k)=                                &
                receive_buffer(ew_address(j,i,k,1))
          END DO
        END DO
      END DO

    ELSE IF (sb_model_domain  ==  mt_global) THEN
      ! No neighbour to my East (ie. at edge of the domain
      ! and it's not a cyclic boundary condition)
      ! Just copy data from last column
      ! NOTE: This block of code should never be executed. It's being left
      !       in so that it can be used in the future by changing the logic
      !       in the ELSE statement above.

      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          DO i=1,halo_x
            field(row_length+i,j,k)=field(row_length,j,k)
          END DO
        END DO
      END DO

    END IF ! IF (he_neighbour(PEast)  /=  NoDomain)

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_E', &
        zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE swap_bounds_e
  
! Purpose:
!   This subroutine takes care of boundary swapping from West
!   extending of arrays at the global boundaries westerly.
  SUBROUTINE swap_bounds_w(              &
      field, row_length, rows, levels,   & ! field
      halo_x, halo_y                     & ! halos
      )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: row_length ! number of points on a row
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: levels     ! number of model levels
    INTEGER, INTENT(IN) :: halo_x     ! size of halo in "i" direction
    INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction  
    REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
        1-halo_y:rows+halo_y,levels)  !Field to have its halos updated

    ! Comdecks
    
    ! Local variables
    
    ! The send and recieve buffers, replacing common blocks
    REAL ::                                                        &
        send_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y))),   &
        receive_buffer(2*levels*((rows*halo_x)+(2*halo_x*halo_y)))

    INTEGER  ::                  &
        full_rows,               & ! number of rows including halo rows
        ew_halo_size,            & ! size of EW halo region
        index_2_start              ! start address of index2 in buffers

    ! Variables used by MPL/GCOM
    INTEGER  :: ierror
    INTEGER  :: istat(mpl_status_size,4)
    INTEGER  :: ireq(4)
    INTEGER  :: nreq
    INTEGER  :: nr
    

    INTEGER ::i,j,k                ! Spatial loop counters


    ! Statement functions for addressing the buffer arrays

    INTEGER :: ew_address

    ! Variables used in statment functions

    INTEGER :: row,point,halo,level,indy

    REAL(KIND=jprb)               :: zhook_handle

    ew_address(row,halo,level,indy)=                                    &
        (indy-1)*(levels*ew_halo_size) +                                &
        (level-1)*ew_halo_size +                                        &
        (halo-1)*full_rows +                                            &
        row+halo_y

    !------------------------------------------------------------------
    ! 0.0 Check if there is anything to do

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_W', &
        zhook_in,zhook_handle)

    
    !------------------------------------------------------------------
    ! 1.0 Initialise variables

    full_rows       = rows + 2*halo_y

    ew_halo_size=full_rows*halo_x

    !------------------------------------------------------------------
    ! 2.0 East-West communications
    
    !---------------------------------------
    ! 2.1.2 Copy the data into buf_send
    
    DO k=1,levels
      DO j=1-halo_y,rows+halo_y
        DO i=1,halo_x

          ! Copy stuff from the Eastern side of the grid
          send_buffer(ew_address(j,i,k,2))=                         &
              field(row_length-halo_x+i,j,k)

        END DO ! I
      END DO ! J
    END DO ! K

    !---------------------------------------
    ! 2.1.3 Send and receive the data


    ! index_2_start points to the start of the second index
    ! within the buffer arrays. The first index contains
    ! data sent from the Western side, the second index
    ! contains data sent from the Eastern side

    nreq = 0
    index_2_start=ew_address(1-halo_y,1,1,2)

    IF (he_neighbour(pwest)  /=  nodomain) THEN

      ! Receive West
      nreq=nreq+1
      CALL mpl_irecv(receive_buffer(index_2_start),    &
          ew_halo_size*levels,                         &
          mpl_real,he_neighbour(pwest),502,            &
          ew_comm,ireq(nreq),ierror)

    END IF

    IF (he_neighbour(peast)  /=  nodomain) THEN

      ! Send East
      nreq=nreq+1
      CALL mpl_isend(send_buffer(index_2_start),       &
          ew_halo_size*levels,                         &
          mpl_real,he_neighbour(peast),502,            &
          ew_comm,ireq(nreq),ierror)

    END IF

    CALL mpl_waitall(nreq,ireq,istat,ierror)


    !---------------------------------------
    ! 2.1.4 Fill the halos with data

    IF (he_neighbour(pwest)  /=  nodomain) THEN

      ! unpack data from receive_buffer into FIELD

      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          DO i=1,halo_x
            field(i-halo_x,j,k)=                                    &
                receive_buffer(ew_address(j,i,k,2))
          END DO
        END DO
      END DO

    ELSE IF (sb_model_domain  ==  mt_global) THEN
      ! No neighbour to my West (ie. at edge of the domain
      ! and it's not a cyclic boundary condition)
      ! Just copy data from first column
      ! NOTE: This block of code should never be executed. It's being left
      !       in so that it can be used in the future by changing the logic
      !       in the ELSE statement above.

      DO k=1,levels
        DO j=1-halo_y,rows+halo_y
          DO i=1,halo_x
            field(i-halo_x,j,k)=field(1,j,k)
          END DO
        END DO
      END DO

    END IF ! IF (he_neighbour(PWest)  /=  NoDomain)

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_W', &
        zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE swap_bounds_w
  
! Purpose:
!   This subroutine takes care of boundary swapping in NS directions
!   Data is swapped across the poles for any non-zero halo size in 
!   the y direction if it is a global model or extending of arrays 
!   at the global boundaries NS.
  SUBROUTINE swap_bounds_ns(             &
      field, row_length, rows, levels,   & ! field
      halo_x, halo_y,                    & ! halos
      field_type, l_vector               & ! supporting information
      )
    
    IMPLICIT NONE

! Implementation
!   The logic flow is non-trivial!
!   The across pole differencing must be  handled carefully!
!   The basic idea is to copy the data to be transferred into a
!   send_buffer array, which is sent to the receiving processor
!   where it arrives in receive_buffer. The receiving processor
!   then copies this data into the appropriate halo region.

! Arguments:
    
    INTEGER, INTENT(IN) :: row_length ! number of points on a row
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
    !                                 !     (not including halos)
    INTEGER, INTENT(IN) :: levels     ! number of model levels
    INTEGER, INTENT(IN) :: halo_x     ! size of halo in "i" direction
    INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction  
    INTEGER, INTENT(IN) :: field_type ! Defines the grid interpolation type
    !                                 !     of the input FIELD (u,v or w)  
    LOGICAL, INTENT(IN) :: L_vector   ! 
    REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
        1-halo_y:rows+halo_y,levels)  !Field to have its halos updated

    
    ! Local variables
    
    ! The send and recieve buffers
    REAL :: send_buffer(2*levels*(row_length*halo_y+(2*halo_x*halo_y)))
    REAL :: receive_buffer(2*levels*(row_length*halo_y+(2*halo_x*halo_y)))
    
    INTEGER  ::                  &
        full_row_length,         & ! length of row including halos
        full_rows,               & ! number of rows including halo rows
        ns_halo_size,            & ! size of NS halo region
        west_halo_source,        & ! first column of data for west halo
        east_halo_source,        & ! first column of data for east halo
        south_halo_source,       &
        north_halo_source,       &
        half_full_row_length,    & ! half of the total EW dimension
        half_row_length,         & ! half of the data (no halo) EW dim
        index_2_start,           & ! start address of index2 in buffers
        north_off,               & ! Offsets to use when copying data
        south_off                  ! to send around poles


    ! Variables used by MPL/GCOM
    INTEGER  :: ierror
    INTEGER  :: istat(mpl_status_size,4)
    INTEGER  :: ireq(4)
    INTEGER  :: nreq
    INTEGER  :: nr
    

    INTEGER  :: i,j,k              ! Spatial loop counters


    LOGICAL :: change_sign         ! .TRUE. if sign change across pole

    ! Statement functions for addressing the buffer arrays

    INTEGER :: ns_address, ew_address

    ! Variables used in statment functions

    INTEGER :: row,point,halo,level,indy

    REAL(KIND=jprb)               :: zhook_handle

    ns_address(point,halo,level,indy)=                                  &
        (indy-1)*(levels*ns_halo_size) +                                &
        (level-1)*ns_halo_size +                                        &
        (halo-1)*full_row_length +                                      &
        point+halo_x

    !------------------------------------------------------------------
    ! 0.0 Check if there is anything to do
    
    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_NS', &
        zhook_in,zhook_handle)

    !------------------------------------------------------------------
    ! 1.0 Initialise variables

    full_row_length = row_length + 2*halo_x
    full_rows       = rows + 2*halo_y

    half_full_row_length=full_row_length/2
    half_row_length=row_length/2

    ns_halo_size=full_row_length*halo_y

    !------------------------------------------------------------------
    ! 3.0 North-South communications

    !  section of code for cyclic N-S boundary conditions
    !  copied from above
    IF(bound(2)  ==  bc_cyclic) THEN

      IF (nproc_y  ==  1) THEN ! only 1 processor north-south
        south_halo_source=rows-halo_y+1       ! copy from opposite end
        north_halo_source=1                   ! of each column
        ! 2.1 Simple case of only one processor

        DO k=1,levels
          DO j=1,halo_y
            DO i=1-halo_x ,row_length+halo_x

              ! Fill southern halo
              field(i,1-halo_y+j-1,k)=field(i,south_halo_source+j-1,k)

              ! Fill northern halo
              field(i,rows+j,k)=field(i,north_halo_source+j-1,k)

            END DO ! I
          END DO ! J
        END DO ! K

        !---------------------------------------
        ! 2.1 Now the more common case of having
        !     a number of processors in the
        !     North-South direction


      ELSE ! If there is more than 1 processor north_south

        !---------------------------------------
        ! 2.1.2 Copy the data into buf_send

        DO k=1,levels
          DO j=1,halo_y
            DO i=1-halo_x  , row_length +halo_x

              ! Copy stuff from the southern side of the grid
              send_buffer(ns_address(i,j,k,1))=field(i,1+j-1,k)

              ! Copy stuff from the northern side of the grid
              send_buffer(ns_address(i,j,k,2))=                         &
                  field(i,rows-halo_y+j,k)

            END DO ! I
          END DO ! J
        END DO ! K

        !---------------------------------------
        ! 2.1.3 Send and receive the data

        !---------------------------------------
        ! 2.1.3.1 Special case of 2 processors
        !         north-south - both sides are
        !         sent to the same processor
        !         (if cyclic BC)

        IF (nproc_y  ==  2 ) THEN

          CALL MPL_Isend(send_buffer,2*ns_halo_size*levels,MPL_Real,    &
              he_neighbour(pnorth),41,ns_comm,ireq(1),ierror)
          CALL MPL_Irecv(receive_buffer,2*ns_halo_size*levels,MPL_Real, &
              he_neighbour(psouth),41,ns_comm,ireq(2),ierror)
          nreq=2

          !---------------------------------------
          ! 2.1.3.2 More general case when there
          !         are more than 2 processors
          !         in the NS direction. Each
          !         halo can be sent seperately
          !         as there is no danger of them
          !         both being sent to the same
          !         processor.

        ELSE ! more than 2 processors North-South

          ! index_2_start points to the start of the second index
          ! within the buffer arrays. The first index contains
          ! data sent from the southern side, the second index
          ! contains data sent from the northern side

          index_2_start=ns_address(1-halo_x,1,1,2)
          nreq=0

          IF (neighbour(psouth)  /=  nodomain) THEN
            nreq=nreq+1
            ! Send south
            CALL MPL_Isend(send_buffer,ns_halo_size*levels,MPL_Real, &
                he_neighbour(psouth),31,ns_comm,ireq(nreq),ierror)

          END IF

          IF (neighbour(pnorth)  /=  nodomain) THEN
            nreq=nreq+1
            ! Send north
            CALL MPL_Isend(send_buffer(index_2_start), &
                ns_halo_size*levels,MPL_Real, &
                he_neighbour(pnorth),32,ns_comm,ireq(nreq),ierror)
          END IF

          IF (neighbour(pnorth)  /=  nodomain) THEN

            nreq=nreq+1
            ! Receive from north
            CALL MPL_Irecv(receive_buffer,    &
                ns_halo_size*levels,MPL_Real, &
                he_neighbour(pnorth),31,ns_comm,ireq(nreq),ierror)

          END IF

          IF (neighbour(psouth)  /=  nodomain) THEN

            nreq=nreq+1
            ! Receive from south
            CALL MPL_Irecv(receive_buffer(index_2_start),    &
                ns_halo_size*levels,MPL_Real,                &
                he_neighbour(psouth),32,ns_comm,ireq(nreq),ierror)

          END IF

        END IF ! test on numbers of processors north-south

        CALL MPL_WaitAll(nreq, ireq, istat, ierror )

        !---------------------------------------
        ! 2.1.4 Fill the halos with data

        IF (neighbour(pnorth)  /=  nodomain) THEN

          ! unpack data from receive_buffer into FIELD

          DO k=1,levels
            DO j=1,halo_y
              DO i=1-halo_x,row_length+halo_x
                field(i,j+rows,k)=                               &
                    receive_buffer(ns_address(i,j,k,1))
              END DO
            END DO
          END DO

        END IF ! IF (neighbour(Pnorth)  /=  NoDomain)

        IF (neighbour(psouth)  /=  nodomain) THEN

          ! unpack data from receive_buffer into FIELD

          DO k=1,levels
            DO j=1,halo_y
              DO i=1-halo_x,row_length+halo_x
                field(i,j-halo_y,k)=                             &
                    receive_buffer(ns_address(i,j,k,2))
              END DO
            END DO
          END DO

        END IF ! IF (neighbour(Psouth)  /=  NoDomain)

      END IF ! IF (nproc_y  ==  1)

    ELSE !  .NOT. bc_cyclic in NS
      ! Set up some variables

      ! Set up the offsets. When copying data that is to be passed over
      ! the pole, on wind (u or v) grid, then copy data one row away
      ! from the pole

      north_off=0
      south_off=0
     IF (l_vatpoles) THEN
      IF ((sb_model_domain  ==  mt_global) .AND. &
          ( field_type  ==  fld_type_v)          &
          ) THEN
        IF (at_extremity(pnorth)) north_off=1
        IF (at_extremity(psouth)) south_off=1
      END IF
     ELSE
      IF ((sb_model_domain  ==  mt_global) .AND. &
          ((field_type  ==  fld_type_p) .OR.     &
          ( field_type  ==  fld_type_u))         &
          ) THEN
        IF (at_extremity(pnorth)) north_off=1
        IF (at_extremity(psouth)) south_off=1
      END IF
     END IF ! vatpoles

      ! Set up the sign factor. If L_VECTOR is true and data has been passed
      ! over the poles in a global model, then the variables must change
      ! sign

      IF (.NOT. ((l_vector) .AND.                &
          (sb_model_domain  ==  mt_global))) THEN
        change_sign=.FALSE.
      ELSE
        change_sign=.TRUE.
      END IF

      !---------------------------------------
      ! 3.1 Copy data into the send_buffer
      !     But not if:
      !       - Not a global model and the data is at the North/South edge
      !       - A global model at the North/South edge but only 1 processor
      !         in the East-West direction.

      IF (.NOT. (at_extremity(psouth) .AND.      &
          ((nproc_x  ==  1)  .OR.                &
          (sb_model_domain  /=  mt_global))))    &
          THEN

        DO k=1,levels
          DO j=1,halo_y
            DO i=1-halo_x,row_length+halo_x

              send_buffer(ns_address(i,j,k,1))=  &
                  field(i,j+south_off,k)

            END DO
          END DO
        END DO
      END IF

      IF (.NOT. (at_extremity(pnorth) .AND.      &
          ((nproc_x  ==  1) .OR.                 &
          (sb_model_domain  /=  mt_global))))    &
          THEN

        DO k=1,levels
          DO j=1,halo_y
            DO i=1-halo_x,row_length+halo_x

              send_buffer(ns_address(i,j,k,2))=  &
                  field(i,rows-halo_y-north_off+j,k)

            END DO
          END DO
        END DO
      END IF

      !---------------------------------------
      ! 3.2 Send and receive the data

      !---------------------------------------
      ! 3.2.1 The special case where nproc_y=1
      !       Both buffers are sent to the
      !       same processor as one message

      IF ((nproc_y  ==  1) .AND. (sb_model_domain  ==  mt_global) .AND. &
          (nproc_x  >   1)) THEN

        CALL MPL_Isend(send_buffer,         &
            2*ns_halo_size*levels,MPL_Real, &
            he_neighbour(psouth),32,ns_comm,ireq(1),ierror)

        CALL MPL_IRecv(receive_buffer,      &
            2*ns_halo_size*levels,MPL_Real, &
            he_neighbour(psouth),32,ns_comm,ireq(2),ierror)

        CALL MPL_WaitAll(2, ireq, istat, ierror )

      END IF

      !---------------------------------------
      ! 3.2.2 The more general case, where
      !       each buffer is sent to a
      !       different processor

      index_2_start=ns_address(1-halo_x,1,1,2)

      nreq=0
      IF (at_extremity(psouth)) THEN

        IF ((neighbour(psouth)  /=  nodomain) .AND.                &
            (neighbour(psouth)  /=  mype)) THEN

          nreq=nreq+1
          CALL mpl_irecv(receive_buffer(index_2_start),            &
              ns_halo_size*levels,                                 &
              mpl_real,he_neighbour(psouth),11,ns_comm,ireq(nreq), &
              ierror)

        END IF

      ELSE ! not at the South

        nreq=nreq+1
        CALL mpl_irecv(receive_buffer(index_2_start),ns_halo_size*levels, &
            mpl_real,he_neighbour(psouth),14,ns_comm,ireq(nreq),          &
            ierror)

      END IF

      IF (at_extremity(pnorth)) THEN

        IF ((neighbour(pnorth)  /=  nodomain) .AND.                  &
            (neighbour(pnorth)  /=  mype)) THEN

          nreq=nreq+1
          CALL mpl_irecv(receive_buffer,ns_halo_size*levels,         &
              mpl_real,he_neighbour(pnorth),13,ns_comm,ireq(nreq),   &
              ierror)

        END IF

      ELSE

        nreq=nreq+1
        CALL mpl_irecv(receive_buffer,ns_halo_size*levels,           &
            mpl_real,he_neighbour(pnorth),12,ns_comm,ireq(nreq),     &
            ierror)

      END IF

      IF (at_extremity(psouth)) THEN

        IF ((neighbour(psouth)  /=  nodomain) .AND.                  &
            (neighbour(psouth)  /=  mype)) THEN

          nreq=nreq+1
          CALL mpl_isend(send_buffer,ns_halo_size*levels,mpl_real,   &
              he_neighbour(psouth),11,ns_comm,ireq(nreq),ierror)

        END IF

      ELSE ! not at the South

        nreq=nreq+1
        CALL mpl_isend(send_buffer,ns_halo_size*levels,mpl_real,     &
            he_neighbour(psouth),12,ns_comm,ireq(nreq),ierror)

      END IF

      IF (at_extremity(pnorth)) THEN

        IF ((neighbour(pnorth)  /=  nodomain) .AND.                      &
            (neighbour(pnorth)  /=  mype)) THEN

          nreq=nreq+1
          CALL mpl_isend(send_buffer(index_2_start),ns_halo_size*levels, &
              mpl_real,he_neighbour(pnorth),13,ns_comm,ireq(nreq),       &
              ierror)

        END IF

      ELSE ! not at the North

        nreq=nreq+1
        CALL mpl_isend(send_buffer(index_2_start),ns_halo_size*levels,   &
            mpl_real,he_neighbour(pnorth),14,ns_comm,ireq(nreq),         &
            ierror)

      END IF

      CALL mpl_waitall( nreq, ireq, istat, ierror )

      !---------------------------------------
      ! 3.3 Fill the halos with data

      !---------------------------------------
      ! 3.3.1 Southern halo

      IF (at_extremity(psouth)) THEN

        IF (neighbour(psouth)  ==  nodomain) THEN
          IF (sb_model_domain  /=  mt_lam) THEN

            ! Just copy adjacent rows into halo area
            ! NOTE: This block of code should never be executed. 
            !       It's being left in so that it can be used in the 
            !       future by changing the logic in the IF statement above.

            DO k=1,levels
              DO j=1,halo_y
                DO i=1-halo_x,row_length+halo_x
                  field(i,j-halo_y,k)=field(i,1,k)
                END DO
              END DO
            END DO
          END IF ! IF (sb_Model_domain  /=  mt_lam)

        ELSE IF (neighbour(psouth)  ==  mype) THEN
          ! Local across pole difference

          IF (change_sign) THEN
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,half_full_row_length
                  field(half_row_length+i,1-j,k)=                     &
                      -field(i,j+south_off,k)
                  field(i-halo_x,1-j,k)=                              &
                      -field(half_row_length+i-halo_x,j+south_off,k)
                END DO
              END DO
            END DO
          ELSE ! don't change sign
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,half_full_row_length
                  field(half_row_length+i,1-j,k)=                     &
                      field(i,j+south_off,k)
                  field(i-halo_x,1-j,k)=                              &
                      field(half_row_length+i-halo_x,j+south_off,k)
                END DO
              END DO
            END DO
          END IF ! IF (change_sign)

        ELSE ! data is receive_buffer(index 2)

          IF (change_sign) THEN
            DO k=1,levels
              DO j=1,halo_y
                DO i=1-halo_x,row_length+halo_x
                  field(i,1-j,k)=                                       &
                      -receive_buffer(ns_address(i,j,k,2))
                END DO
              END DO
            END DO
          ELSE ! don't change sign
            DO k=1,levels
              DO j=1,halo_y
                DO i=1-halo_x,row_length+halo_x
                  field(i,1-j,k)=                                       &
                      receive_buffer(ns_address(i,j,k,2))
                END DO
              END DO
            END DO
          END IF !  IF (change_sign)

        END IF ! What type of South extremity

      ELSE ! IF (at_extremity(PSouth)

        ! not at a South extremity

        DO k=1,levels
          DO j=1,halo_y
            DO i=1-halo_x,row_length+halo_x
              field(i,j-halo_y,k)=                                      &
                  receive_buffer(ns_address(i,j,k,2))
            END DO
          END DO
        END DO

      END IF ! IF (at_extremity(PSouth)


      !---------------------------------------
      ! 3.3.2 Northern halo

      IF (at_extremity(pnorth)) THEN

        IF (neighbour(pnorth)  ==  nodomain) THEN
          IF (sb_model_domain  /=  mt_lam) THEN
            ! Just copy adjacent rows into halo area
            ! NOTE: This block of code should never be executed. It's being left
            !       in so that it can be used in the future by changing the logic
            !       in the IF statement above.

            DO k=1,levels
              DO j=1,halo_y
                DO i=1-halo_x,row_length+halo_x
                  field(i,rows+j,k)=field(i,rows,k)
                END DO
              END DO
            END DO
          END IF ! IF (sb_Model_domain  /=  mt_lam)

        ELSE IF (neighbour(pnorth)  ==  mype) THEN
          ! Local across pole difference

          IF (change_sign) THEN
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,half_full_row_length

                  field(half_row_length+i,rows+j,k)=                    &
                      -field(i,rows-j+1-north_off,k)
                  field(i-halo_x,rows+j,k)=                             &
                      -field(half_row_length+i-halo_x,                  &
                      rows-j+1-north_off,k)

                END DO
              END DO
            END DO
          ELSE ! don't change sign
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,half_full_row_length

                  field(half_row_length+i,rows+j,k)=                    &
                      field(i,rows-j+1-north_off,k)
                  field(i-halo_x,rows+j,k)=                             &
                      field(half_row_length+i-halo_x,                   &
                      rows-j+1-north_off,k)

                END DO
              END DO
            END DO
          END IF ! IF (change_sign)

        ELSE ! data is receive_buffer(index 1)

          IF (change_sign) THEN
            DO k=1,levels
              DO j=1,halo_y
                DO i=1-halo_x,row_length+halo_x
                  field(i,rows+j,k)=                                    &
                      -receive_buffer(ns_address(i,halo_y-j+1,k,1))
                END DO
              END DO
            END DO
          ELSE ! don't change sign
            DO k=1,levels
              DO j=1,halo_y
                DO i=1-halo_x,row_length+halo_x
                  field(i,rows+j,k)=                                    &
                      receive_buffer(ns_address(i,halo_y-j+1,k,1))
                END DO
              END DO
            END DO
          END IF ! IF (change_sign)

        END IF ! What type of North extremity

      ELSE ! IF (at_extremity(PNorth)

        ! not at a North extremity

        DO k=1,levels
          DO j=1,halo_y
            DO i=1-halo_x,row_length+halo_x
              field(i,rows+j,k)=                                        &
                  receive_buffer(ns_address(i,j,k,1))
            END DO
          END DO
        END DO

      END IF ! IF (at_extremity(PNorth)
    END IF

    IF (lhook) CALL dr_hook('HALO_EXCHANGE:SWAP_BOUNDS_NS', &
        zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE swap_bounds_ns

END MODULE Halo_Exchange
