! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Direct MPI version of swap_bounds to deal with multiple variables
! at once

SUBROUTINE swap_bounds_2d_mv(  &
   input_fields,               &        ! Fields to be swapped
   n_multi,                    &        ! The number of Fields to swap
   row_length,                 &        ! field size
   halo_x, halo_y)                   ! halos

USE dynamics_grid_mod, ONLY: l_vatpoles

USE mpl, ONLY :             &
         mpl_real,          &
         mpl_status_size
USE swapable_field_mod, ONLY : &
    swapable_field_pointer_type

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParVars
USE domain_params
IMPLICIT NONE

!  This code performs the same process as swap_bounds but
!  can swap multiple variables at once and hence utilise bandwidth
!  better. This is the MPL version.

!  Note that it can only deal with multiple variables with the same
!  row_length and same halo size.

! Purpose:
!   This subroutine takes care of all boundary swapping and
!   extending of arrays at the global boundaries. Data is swapped
!   across the poles for any non-zero halo size in the y direction
!   if it is a global model.

! Implementation
!   The logic flow is non-trivial!
!   The across pole differencing in particular must be handled carefully!
!   The basic idea is to copy the data to be transferred into a
!   send_buffer array, which is sent to the receiving processor
!   where it arrives in receive_buffer. The receiving processor
!   then copies this data into the appropriate halo region.
!   The East/West halos are done first, then the North/South
!   halos.

!   Note that due to the fact that pointers to the data are used,
!   addressing is (1:row_length+2*halo_x, 1:rows+2*halo_y) rather
!   than (1-halo_x:row_length+halo_x, 1-halo_y:rows+halo_y) as used
!   elsewhere. This is unavoidable as pointers change the addressing
!   mode.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Comdecks

! Arguments:

INTEGER, INTENT(IN)  :: n_multi      ! number of fields to be swapped
INTEGER, INTENT(IN)  :: row_length   ! number of points on row (no halos)
INTEGER, INTENT(IN)  :: halo_x       ! sixze of "i" halo
INTEGER, INTENT(IN)  :: halo_y       ! sixze of "j" halo

! Fields to swap
TYPE(swapable_field_pointer_type), TARGET :: input_fields(n_multi)

! Local scalar variables
INTEGER :: rows               ! number of rows in field (no halos)
INTEGER :: field_type    ! The grid type of the field (u,v,p)
INTEGER :: i,j,k         ! Spatial loop counters
INTEGER :: info          ! GCOM return code
INTEGER :: length
INTEGER :: i_field
INTEGER :: ierror        ! MPI return code
INTEGER :: nreq_r_ew     ! number of MPI recv requests ew
INTEGER :: nreq_s_ew     ! number of MPI send requests ew
INTEGER :: nreq_r_ns     ! number of MPI recv requests ns
INTEGER :: nreq_s_ns     ! number of MPI send requests ns
INTEGER :: full_row_length       ! length of row including halos
INTEGER :: max_full_rows         ! max no of rows including halo rows
INTEGER :: ew_halo_size          ! size of EW halo region
INTEGER :: ns_halo_size          ! size of NS halo region
INTEGER :: west_halo_source      ! first column of data for west halo
INTEGER :: east_halo_source      ! first column of data for east halo
INTEGER :: south_halo_source
INTEGER :: north_halo_source
INTEGER :: half_full_row_length  ! half of the total EW dimension
INTEGER :: half_row_length       ! half of the data (no halo) EW dim
INTEGER :: index_2_start         ! start address of index2 in buffers
INTEGER :: max_rows              ! maximum rows used over all fields
INTEGER :: buffer_size           ! size of send/recv buffers
INTEGER :: my_comm               ! Communicator
INTEGER :: errorstatus           ! error code

LOGICAL :: l_vector      ! TRUE if a vector field


! Local arrays
INTEGER :: istat(mpl_status_size,4)  ! MPI status
INTEGER :: ireq_r_ew(4)              ! MPI requests
INTEGER :: ireq_s_ew(4)              ! MPI requests
INTEGER :: ireq_r_ns(4)              ! MPI requests
INTEGER :: ireq_s_ns(4)              ! MPI requests
INTEGER :: north_off(n_multi)        ! Offsets to use when copying data
INTEGER :: south_off(n_multi)        ! to send around poles

CHARACTER (LEN=80) :: cmessage   ! Error message
CHARACTER (LEN=*), PARAMETER  :: routinename='swap_bounds_2d_mv'  ! routine name

REAL, POINTER :: field(:, :)

! Send and receive buffers to be allocated dynamically
REAL, ALLOCATABLE :: send_buffer(:)
REAL, ALLOCATABLE :: receive_buffer(:)

LOGICAL :: change_sign(n_multi)     ! .TRUE. if sign change across pole


! Statement functions for addressing the buffer arrays
INTEGER :: ew_address
INTEGER :: ns_address

! Variables used in statement functions

INTEGER :: row
INTEGER :: point
INTEGER :: halo
INTEGER :: index
INTEGER :: j_field

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

ew_address(row,halo,index,j_field)=                                       &
   n_multi*(index-1)* ew_halo_size  +                                     &
   (j_field-1)* ew_halo_size  +                                           &
   (halo-1)*max_full_rows +                                               &
   row

ns_address(point,halo,index,j_field)=                                     &
   n_multi*(index-1)* ns_halo_size  +                                     &
   (j_field-1)* ns_halo_size  +                                           &
   (halo-1)*full_row_length +                                             &
   point


IF (lhook) CALL dr_hook('SWAP_BOUNDS_2D_MV',zhook_in,zhook_handle)

!------------------------------------------------------------------
! 0.0 Check if there is anything to do
IF (((halo_x == 0) .AND. (halo_y == 0)) .OR. n_multi == 0) GO TO 9999

!------------------------------------------------------------------
! 1.0 Initialise variables

! Maximum rows and levels
max_rows   = 0
DO i_field=1, n_multi
  max_rows  = MAX(input_fields(i_field) % rows,   max_rows)

  IF (input_fields(i_field) % levels /= 1) THEN

    CALL ereport(routinename, errorstatus, cmessage)
  END IF
END DO

full_row_length = row_length + 2*halo_x
max_full_rows   = max_rows + 2*halo_y

half_full_row_length = full_row_length/2
half_row_length      = row_length/2

ew_halo_size         = max_full_rows   * halo_x
ns_halo_size         = full_row_length * halo_y


! Allocate buffers for communication
buffer_size =  n_multi * 2*                                               &
               (MAX(max_rows*halo_x, row_length*halo_y) +                 &
               (2 * halo_x * halo_y))

ALLOCATE( send_buffer(buffer_size) )
ALLOCATE( receive_buffer(buffer_size) )

! Get communicator we will be using from GCOM
CALL gc_get_communicator(my_comm,ierror)

!------------------------------------------------------------------
! 2.0 East-West communications
!---------------------------------------
! 2.1 Simple case of only one processor
!     in the East-West direction

IF (halo_x > 0) THEN

nreq_s_ew = 0
nreq_s_ns = 0
IF (nproc_x == 1) THEN ! only 1 processor East-West

  IF (bound(1) == bc_cyclic) THEN        ! cyclic boundary conditions
    west_halo_source=row_length+1        ! copy from opposite end
    east_halo_source=halo_x+1            ! of each row

    DO i_field=1, n_multi
      field       => input_fields(i_field) % field_2d
      rows        = input_fields(i_field) % rows
!CDIR NOVECTOR
      DO i=1,halo_x
        DO j=1,rows + (2 * halo_y)

            ! Fill Western halo
            field(i,j) = field(west_halo_source+i-1,j)

            ! Fill Eastern halo
            field(row_length+halo_x+i,j) =                              &
                            field(east_halo_source+i-1,j)

        END DO ! J
      END DO ! I

    END DO ! loop over fields

  END IF   !  bound(1) == BC_CYCLIC

!---------------------------------------
! 2.1 Now the more common case of having
!     a number of processors in the
!     East-West direction

ELSE ! If there is more than 1 processor East-West

!---------------------------------------
! 2.1.1 Copy the data into send_buffer

  DO i_field=1, n_multi
    field       => input_fields(i_field) % field_2d
    rows        = input_fields(i_field) % rows

    DO j=1,rows + (2*halo_y)
      DO i=1,halo_x

        ! Copy stuff from the Western side of the grid
        send_buffer(ew_address(j,i,1,i_field)) =                        &
                    field(halo_x+i,j)

        ! Copy stuff from the Eastern side of the grid
        send_buffer(ew_address(j,i,2,i_field)) =                        &
                    field(row_length+i,j)

      END DO ! I
    END DO ! J
  END DO ! loop over fields

!---------------------------------------
! 2.1.2 Send and receive the data

!---------------------------------------
! 2.1.2.1 Special case of 2 processors
!         East-West - both sides are
!         sent to the same processor
!         (if cyclic BC)

  IF ((nproc_x == 2) .AND. (bound(1) == bc_cyclic)) THEN

    length = 2 * n_multi * ew_halo_size

    CALL gc_rsend(1,length,neighbour(peast),info,                         &
                   receive_buffer, send_buffer)

    CALL gc_rrecv(1,length,neighbour(pwest),info,                         &
                   receive_buffer, send_buffer)

!---------------------------------------
! 2.1.2.2 More general case when there
!         are more than 2 processors
!         in the EW direction. Each
!         halo can be sent seperately
!         as there is no danger of them
!         both being sent to the same
!         processor.

  ELSE ! more than 2 processors East-West

    ! index_2_start points to the start of the second index
    ! within the buffer arrays. The first index contains
    ! data sent from the Western side, the second index
    ! contains data sent from the Eastern side

    index_2_start=ew_address(1,1,2,1)

    length = n_multi * ew_halo_size

    nreq_r_ew=0
    IF (neighbour(peast) /= nodomain) THEN

      ! Receive from East
      nreq_r_ew=nreq_r_ew+1
      CALL mpl_irecv(receive_buffer,                                      &
                     length, mpl_real, neighbour(peast), 2,               &
                     my_comm, ireq_r_ew(nreq_r_ew), ierror)

    END IF

    IF (neighbour(pwest) /= nodomain) THEN

      ! Receive from West
      nreq_r_ew=nreq_r_ew+1
      CALL mpl_irecv(receive_buffer(index_2_start),                       &
                     length, mpl_real, neighbour(pwest), 3,               &
                     my_comm, ireq_r_ew(nreq_r_ew), ierror)

    END IF

    nreq_s_ew=0
    IF (neighbour(peast) /= nodomain) THEN

      ! Send East
      nreq_s_ew=nreq_s_ew+1
      CALL mpl_isend(send_buffer(index_2_start),                          &
                     length, mpl_real, neighbour(peast), 3,               &
                     my_comm, ireq_s_ew(nreq_s_ew), ierror)
    END IF

    IF (neighbour(pwest) /= nodomain) THEN

      ! Send West
      nreq_s_ew=nreq_s_ew+1
      CALL mpl_isend(send_buffer,                                         &
                     length, mpl_real, neighbour(pwest), 2,               &
                     my_comm, ireq_s_ew(nreq_s_ew), ierror)

    END IF

    CALL mpl_waitall( nreq_r_ew, ireq_r_ew, istat, ierror )

  END IF ! test on numbers of processors East-West

  CALL mpl_waitall( nreq_s_ew, ireq_s_ew, istat, ierror )

!---------------------------------------
! 2.1.2 Fill the halos with data

  IF (neighbour(peast) /= nodomain) THEN

    ! unpack data from receive_buffer into field

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      rows        = input_fields(i_field) % rows

      DO j=1,rows + (2*halo_y)
        DO i=1,halo_x
          field(row_length+halo_x+i,j)=                                 &
             receive_buffer(ew_address(j,i,1,i_field))
        END DO
      END DO
    END DO ! loop over fields

  ELSE IF (sb_model_domain == mt_global) THEN
       ! No neighbour to my East (ie. at edge of the domain
       ! and it's not a cyclic boundary condition)
       ! Just copy data from last column
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the ELSE statement above.

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      rows        = input_fields(i_field) % rows

      DO j=1,rows+ (2*halo_y)
        DO i=1,halo_x
          field(row_length+halo_x+i,j)=field(row_length+halo_x,j)
        END DO
      END DO
    END DO ! loop over fields

  END IF ! IF (neighbour(PEast) /= NoDomain)

  IF (neighbour(pwest) /= nodomain) THEN

    ! unpack data from receive_buffer into field

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      rows        = input_fields(i_field) % rows

      DO j=1,rows+ (2*halo_y)
        DO i=1,halo_x
          field(i,j) = receive_buffer(ew_address(j,i,2,i_field))
        END DO
      END DO
    END DO ! loop over fields

  ELSE IF (sb_model_domain == mt_global) THEN
       ! No neighbour to my West (ie. at edge of the domain
       ! and it's not a cyclic boundary condition)
       ! Just copy data from first column
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the ELSE statement above.

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      rows        = input_fields(i_field) % rows

      DO j=1,rows+ (2*halo_y)
        DO i=1,halo_x
          field(i,j) = field(1+halo_x,j)
        END DO
      END DO
    END DO ! loop over fields

  END IF ! IF (neighbour(PWest) /= NoDomain)

END IF ! IF (nproc_x == 1)

END IF ! halo_x > 0

!------------------------------------------------------------------
! 3.0 North-South communications
!  section of code for cyclic N-S boundary conditions
IF (halo_y > 0) THEN

IF(bound(2) == bc_cyclic) THEN

  IF (nproc_y == 1) THEN ! only 1 processor north-south
  ! 2.1 Simple case of only one processor

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      rows        = input_fields(i_field) % rows

      south_halo_source=rows+1              ! copy from opposite end
      north_halo_source=1                   ! of each column

      DO j=1,halo_y
        DO i=1 ,row_length+ (2*halo_x)

          ! Fill southern halo
          field(i,j)=field(i,south_halo_source+j-1)

          ! Fill northern halo
          field(i,rows+halo_y+j)=field(i,north_halo_source+j-1)

        END DO ! I
      END DO ! J
    END DO ! loop over fields

  !---------------------------------------
  ! 2.1 Now the more common case of having
  !     a number of processors in the
  !     North-South direction

  ELSE ! If there is more than 1 processor north_south

  !---------------------------------------
  ! 2.1.1 Copy the data into buf_send

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      rows        = input_fields(i_field) % rows

      DO j=1,halo_y
        DO i=1, row_length + (2*halo_x)

          ! Copy stuff from the southern side of the grid
          send_buffer(ns_address(i,j,1,i_field)) = field(i,halo_y+j)

          ! Copy stuff from the northern side of the grid
          send_buffer(ns_address(i,j,2,i_field)) = field(i,rows+j)

        END DO ! I
      END DO ! J
    END DO ! loop over fields

  !---------------------------------------
  ! 2.1.2 Send and receive the data

  !---------------------------------------
  ! 2.1.2.1 Special case of 2 processors
  !         north-south - both sides are
  !         sent to the same processor
  !         as cyclic BC

   IF ( nproc_y == 2 ) THEN

     length=2*n_multi*ns_halo_size

      CALL gc_rsend(1,length,neighbour(pnorth),info,                      &
                     receive_buffer, send_buffer)

      CALL gc_rrecv(1,length,neighbour(psouth),info,                      &
                     receive_buffer, send_buffer)

  !---------------------------------------
  ! 2.1.2.2 More general case when there
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

      index_2_start=ns_address(1,1,2,1)

      length = n_multi * ns_halo_size

      IF (neighbour(psouth) /= nodomain) THEN

        ! Send south
        CALL gc_rsend(2,length,neighbour(psouth),info,                    &
                       receive_buffer, send_buffer)
      END IF

      IF (neighbour(pnorth) /= nodomain) THEN

        ! Send north
        CALL gc_rsend(3,length,neighbour(pnorth),info,                    &
                       receive_buffer(index_2_start),                     &
                       send_buffer(index_2_start))
      END IF

      IF (neighbour(pnorth) /= nodomain) THEN

        ! Receive from north
        CALL gc_rrecv(2,length,neighbour(pnorth),info,                    &
                       receive_buffer, send_buffer)
      END IF

      IF (neighbour(psouth) /= nodomain) THEN

        ! Receive from south
        CALL gc_rrecv(3,length,neighbour(psouth),info,                    &
                       receive_buffer(index_2_start),                     &
                       send_buffer(index_2_start))

      END IF

    END IF ! test on numbers of processors north-south

  !---------------------------------------
  ! 2.1.2 Fill the halos with data

    IF (neighbour(pnorth) /= nodomain) THEN

      ! unpack data from receive_buffer into field

      DO i_field=1, n_multi
        field      => input_fields(i_field) % field_2d
        rows        = input_fields(i_field) % rows

        DO j=1,halo_y
          DO i=1,row_length+ (2*halo_x)
            field(i,j+halo_y+rows)=                                   &
            receive_buffer(ns_address(i,j,1,i_field))

          END DO
        END DO
      END DO ! loop over fields

    END IF ! IF (neighbour(Pnorth) /= NoDomain)

    IF (neighbour(psouth) /= nodomain) THEN

      ! unpack data from receive_buffer into field

      DO i_field=1, n_multi
        field      => input_fields(i_field) % field_2d
        DO j=1,halo_y
          DO i=1,row_length+ (2*halo_x)
            field(i,j) = receive_buffer(ns_address(i,j,2,i_field))
          END DO
        END DO
      END DO ! loop over fields


    END IF ! IF (neighbour(Psouth) /= NoDomain)

  END IF ! IF (nproc_y == 1)

ELSE                 !!! bc_cyclic in NS

  ! Set up some variables

  ! Set up the offsets. When copying data that is to be passed over
  ! the pole, on wind (u or v) grid, then copy data one row away
  ! from the pole

  DO i_field=1, n_multi
    north_off(i_field)=0
    south_off(i_field)=0
    field_type = input_fields(i_field) % field_type
    IF (l_vatpoles) THEN
      IF ((sb_model_domain == mt_global) .AND.                              &
           (field_type  ==  fld_type_v)) THEN
        IF (at_extremity(pnorth)) north_off(i_field)=1
        IF (at_extremity(psouth)) south_off(i_field)=1
      END IF
    ELSE
      IF ((sb_model_domain == mt_global) .AND.                              &
              ((field_type == fld_type_p) .OR.                              &
               (field_type == fld_type_u))) THEN
        IF (at_extremity(pnorth)) north_off(i_field)=1
        IF (at_extremity(psouth)) south_off(i_field)=1
      END IF
    END IF ! vatpoles

  ! Set up the sign factor. If l_vector is true and data has been passed
  ! over the poles in a global model, then the variables must change
  ! sign

    l_vector = input_fields(i_field) % vector
    IF (.NOT. ((l_vector) .AND. (sb_model_domain == mt_global))) THEN
      change_sign(i_field)=.FALSE.
    ELSE
      change_sign(i_field)=.TRUE.
    END IF
  END DO

  !---------------------------------------
  ! 3.1 Copy data into the send_buffer
  !     But not if:
  !       - Not a global model and the data is at the North/South edge
  !       - A global model at the North/South edge but only 1 processor
  !         in the East-West direction.

  IF (.NOT. (at_extremity(psouth) .AND.                                   &
              ((nproc_x == 1)  .OR.                                       &
              (sb_model_domain /= mt_global))))                           &
     THEN

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      DO j=1,halo_y
        DO i=1,row_length+ (2*halo_x)

          send_buffer(ns_address(i,j,1,i_field)) =                        &
             field(i,j+halo_y+south_off(i_field))

        END DO
      END DO
    END DO ! loop over fields
  END IF

  IF (.NOT. (at_extremity(pnorth) .AND.                                   &
              ((nproc_x == 1) .OR.                                        &
              (sb_model_domain /= mt_global))))                           &
     THEN

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      rows        = input_fields(i_field) % rows
      DO j=1,halo_y
        DO i=1,row_length+ (2*halo_x)

          send_buffer(ns_address(i,j,2,i_field)) =                        &
             field(i,rows-north_off(i_field)+j)

        END DO
      END DO
    END DO ! loop over fields
  END IF

  !---------------------------------------
  ! 3.2 Send and receive the data

  !---------------------------------------
  ! 3.2.1 The special case where nproc_y=1
  !       Both buffers are sent to the
  !       same processor as one message

  IF ((nproc_y == 1) .AND. (sb_model_domain == mt_global) .AND.           &
       (nproc_x > 1)) THEN

    length = 2 * n_multi * ns_halo_size

    CALL gc_rsend(10,length,neighbour(psouth),info,                       &
                     receive_buffer, send_buffer)

    CALL gc_rrecv(10,length,neighbour(psouth),info,                       &
                   receive_buffer, send_buffer)
  END IF

  !---------------------------------------
  ! 3.2.2 The more general case, where
  !       each buffer is sent to a
  !       different processor

  index_2_start=ns_address(1,1,2,1)

  nreq_r_ns=0

  length = n_multi * ns_halo_size

  IF (at_extremity(psouth)) THEN

    IF ((neighbour(psouth) /= nodomain) .AND.                             &
         (neighbour(psouth) /= mype)) THEN

       nreq_r_ns=nreq_r_ns+1
       CALL mpl_irecv(receive_buffer(index_2_start),                      &
         length, mpl_real, neighbour(psouth), 11, my_comm,                &
         ireq_r_ns(nreq_r_ns), ierror)

    END IF

  ELSE ! not at the South

     nreq_r_ns=nreq_r_ns+1
     CALL mpl_irecv(receive_buffer(index_2_start),                        &
       length, mpl_real, neighbour(psouth), 14, my_comm,                  &
       ireq_r_ns(nreq_r_ns), ierror)

  END IF

  IF (at_extremity(pnorth)) THEN

    IF ((neighbour(pnorth) /= nodomain) .AND.                             &
         (neighbour(pnorth) /= mype)) THEN

       nreq_r_ns=nreq_r_ns+1
       CALL mpl_irecv(receive_buffer,                                     &
         length, mpl_real, neighbour(pnorth), 13, my_comm,                &
         ireq_r_ns(nreq_r_ns), ierror)

    END IF

  ELSE

       nreq_r_ns=nreq_r_ns+1
       CALL mpl_irecv(receive_buffer,                                     &
         length, mpl_real, neighbour(pnorth), 12, my_comm,                &
         ireq_r_ns(nreq_r_ns), ierror)

  END IF

  nreq_s_ns=0
  IF (at_extremity(psouth)) THEN

    IF ((neighbour(psouth) /= nodomain) .AND.                             &
         (neighbour(psouth) /= mype)) THEN

       nreq_s_ns=nreq_s_ns+1
       CALL mpl_isend(send_buffer,                                        &
         length, mpl_real, neighbour(psouth), 11, my_comm,                &
         ireq_s_ns(nreq_s_ns),ierror)

    END IF

  ELSE ! not at the South

       nreq_s_ns=nreq_s_ns+1
       CALL mpl_isend(send_buffer,                                        &
         length, mpl_real, neighbour(psouth), 12, my_comm,                &
         ireq_s_ns(nreq_s_ns),ierror)

  END IF

  IF (at_extremity(pnorth)) THEN

    IF ((neighbour(pnorth) /= nodomain) .AND.                             &
         (neighbour(pnorth) /= mype)) THEN

     nreq_s_ns=nreq_s_ns+1
     CALL mpl_isend(send_buffer(index_2_start),                           &
       length, mpl_real, neighbour(pnorth), 13, my_comm,                  &
       ireq_s_ns(nreq_s_ns),ierror)

    END IF

  ELSE ! not at the North

     nreq_s_ns=nreq_s_ns+1
     CALL mpl_isend(send_buffer(index_2_start),                           &
       length, mpl_real, neighbour(pnorth), 14, my_comm,                  &
       ireq_s_ns(nreq_s_ns),ierror)

  END IF

  CALL mpl_waitall( nreq_r_ns, ireq_r_ns, istat, ierror )
  !---------------------------------------
  ! 3.3 Fill the halos with data

  !---------------------------------------
  ! 3.3.1 Southern halo

  IF (at_extremity(psouth)) THEN

    IF (neighbour(psouth) == nodomain) THEN
      IF (sb_model_domain /= mt_lam) THEN

      ! Just copy adjacent rows into halo area
  ! NOTE: This block of code should never be executed. It's being left
  !       in so that it can be used in the future by changing the logic
  !       in the IF statement above.

        DO i_field=1, n_multi
          field      => input_fields(i_field) % field_2d
          DO j=1,halo_y
            DO i=1,row_length+(2*halo_x)
              field(i,j) = field(i,halo_y+1)
            END DO
          END DO
        END DO ! loop over fields
      END IF ! IF (sb_Model_domain /= mt_lam)

    ELSE IF (neighbour(psouth) == mype) THEN
    ! Local across pole difference

      DO i_field=1, n_multi
        field      => input_fields(i_field) % field_2d
        IF (change_sign(i_field)) THEN
          DO j=1,halo_y
            DO i=1,half_full_row_length
              field(half_row_length+halo_x+i,1-j+halo_y)=               &
                 -field(i+halo_x,j+halo_y+south_off(i_field))
              field(i,1-j+halo_y)=                                      &
                 -field(half_row_length+i,j+halo_y+                     &
                 south_off(i_field))
            END DO
          END DO
        ELSE ! don't change sign
          DO j=1,halo_y
            DO i=1,half_full_row_length
              field(half_row_length+halo_x+i,1-j+halo_y)=               &
                 field(i+halo_x,j+halo_y+south_off(i_field))
              field(i,1-j+halo_y)=                                      &
                 field(half_row_length+i,j+halo_y+                      &
                 south_off(i_field))
            END DO
          END DO
        END IF ! IF (change_sign(i_field))
      END DO ! loop over fields

    ELSE ! data is receive_buffer(index 2)

      DO i_field=1, n_multi
        field      => input_fields(i_field) % field_2d
        IF (change_sign(i_field)) THEN
          DO j=1,halo_y
            DO i=1,row_length+(2*halo_x)
              field(i,1-j+halo_y)=                                      &
              -receive_buffer(ns_address(i,j,2,i_field))
            END DO
          END DO
        ELSE ! don't change sign
          DO j=1,halo_y
            DO i=1,row_length+(2*halo_x)
              field(i,1-j+halo_y)=                                      &
               receive_buffer(ns_address(i,j,2,i_field))
            END DO
          END DO
        END IF !  IF (change_sign(i_field))
      END DO ! loop over fields

    END IF ! What type of South extremity

  ELSE ! IF (at_extremity(PSouth)

    ! not at a South extremity

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      DO j=1,halo_y
        DO i=1,row_length+(2*halo_x)
          field(i,j)=                                                   &
             receive_buffer(ns_address(i,j,2,i_field))
        END DO
      END DO
    END DO ! loop over fields

  END IF ! IF (at_extremity(PSouth)


  !---------------------------------------
  ! 3.3.2 Northern halo

  IF (at_extremity(pnorth)) THEN

    IF (neighbour(pnorth) == nodomain) THEN
      IF (sb_model_domain /= mt_lam) THEN
      ! Just copy adjacent rows into halo area
  ! NOTE: This block of code should never be executed. It's being left
  !       in so that it can be used in the future by changing the logic
  !       in the IF statement above.

        DO i_field=1, n_multi
          field      => input_fields(i_field) % field_2d
          rows        = input_fields(i_field) % rows
          DO j=1,halo_y
          DO i=1,row_length+(2*halo_x)
            field(i,rows+halo_y+j) = field(i,rows+halo_y)
          END DO
          END DO
        END DO ! loop over fields
      END IF

    ELSE IF (neighbour(pnorth) == mype) THEN
    ! Local across pole difference

      DO i_field=1, n_multi
        field      => input_fields(i_field) % field_2d
        rows        = input_fields(i_field) % rows
        IF (change_sign(i_field)) THEN
          DO j=1,halo_y
            DO i=1,half_full_row_length

              field(half_row_length+halo_x+i,rows+halo_y+j)=            &
                 -field(i+halo_x,rows+halo_y-j+1-north_off(i_field))
              field(i,rows+j+halo_y)=                                   &
                 -field(half_row_length+i,                              &
                      rows-j+halo_y+1-north_off(i_field))

            END DO
          END DO
        ELSE ! don't change sign
          DO j=1,halo_y
            DO i=1,half_full_row_length

              field(half_row_length+halo_x+i,rows+halo_y+j)=            &
                 field(i+halo_x,rows+halo_y-j+1-north_off(i_field))
              field(i,rows+j+halo_y)=                                   &
                 field(half_row_length+i,                               &
                     rows-j+1+halo_y-north_off(i_field))

              END DO
            END DO
        END IF ! IF (change_sign(i_field))
      END DO ! loop over fields

    ELSE ! data is receive_buffer(index 1)

      DO i_field=1, n_multi
        field      => input_fields(i_field) % field_2d
        rows        = input_fields(i_field) % rows
        IF (change_sign(i_field)) THEN
          DO j=1,halo_y
            DO i=1,row_length+(2*halo_x)
              field(i,rows+halo_y+j)=                                   &
                 -receive_buffer(                                       &
                 ns_address(i,halo_y-j+1,1,i_field))
            END DO
          END DO
        ELSE ! don't change sign
          DO j=1,halo_y
            DO i=1,row_length+(2*halo_x)
              field(i,rows+halo_y+j)=                                   &
                 receive_buffer(                                        &
                 ns_address(i,halo_y-j+1,1,i_field))
            END DO
          END DO
        END IF ! IF (change_sign(i_field))
      END DO ! loop over fields

    END IF ! What type of North extremity

  ELSE ! IF (at_extremity(PNorth)

    ! not at a North extremity

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field_2d
      rows        = input_fields(i_field) % rows
      DO j=1,halo_y
        DO i=1,row_length+ (2*halo_x)
          field(i,rows+halo_y+j)=                                       &
             receive_buffer(ns_address(i,j,1,i_field))
        END DO
      END DO
    END DO ! loop over fields
  END IF ! IF (at_extremity(PNorth)
END IF ! BC_CYCLIC  for NS

CALL mpl_waitall( nreq_s_ns, ireq_s_ns, istat, ierror )

END IF ! halo_y > 0

9999 CONTINUE

IF (ALLOCATED(send_buffer))    DEALLOCATE(send_buffer)
IF (ALLOCATED(receive_buffer)) DEALLOCATE(receive_buffer)

! Fill external halos for LAMs
IF (sb_model_domain == mt_lam) THEN
!DEPENDS ON: fill_external_halos
  DO i_field = 1, n_multi
    field  => input_fields(i_field) % field_2d
    rows   =  input_fields(i_field) % rows
      CALL fill_external_halos(field, row_length, rows, 1, &
                               halo_x, halo_y)
   END DO
END IF

IF (lhook) CALL dr_hook('SWAP_BOUNDS_2D_MV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE swap_bounds_2d_mv
