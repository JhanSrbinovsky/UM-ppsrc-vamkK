! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine Bi_Linear_H.

SUBROUTINE bi_linear_h(                                                        &
  data_in, lambda_out, phi_out,                                                &
  dim_i_in, dim_j_in, dim_k_in,                                                &
  dim_i_out, dim_j_out, dim_k_out,                                             &
  row_length, rows,                                                            &
  i_out_in, j_out_in,                                                          &
  weight_lambda_in, weight_phi_in,                                             &
  model_domain,                                                                &
  me, n_procx, n_procy, nproc,                                                 &
  halo_i, halo_j, datastart,                                                   &
  g_row_length, g_i_pe, at_extremity,                                          &
  g_rows,       g_j_pe, l_2dcomm,                                              &
  size_2dcomm, group_2dcomm,                                                   &
  pole_handling, proc_all_group,                                               &
  l_sl_halo_reprod, l_regular,                                                 &
  data_out)

! Purpose:
!          Performs bi-linear horizontal interpolation of a field
!          defined on a w grid to another grid. Input data can be on a
!          sphere or be a rectangular box. Requested output points must
!          lie inside the region defined for the input Data.

! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.


USE mpl, ONLY :                                                                &
  mpl_integer,                                                                 &
  mpl_real,                                                                    &
  mpl_status_size


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParParams
USE um_types, ONLY: integer32

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(in) :: dim_i_in  ! Dimension of Data_in in i direction.
INTEGER, INTENT(in) :: dim_j_in  ! Dimension of Data_in in j direction.
INTEGER, INTENT(in) :: dim_k_in  ! Dimension of Data_in in k direction.
INTEGER, INTENT(in) :: dim_i_out ! Dimension of Data_out in i direction.
INTEGER, INTENT(in) :: dim_j_out ! Dimension of Data_out in j direction.
INTEGER, INTENT(in) :: dim_k_out ! Dimension of Data_out in j direction.
INTEGER, INTENT(in) :: me        ! My processor number
INTEGER, INTENT(in) :: n_procx   ! Number of processors E-W
INTEGER, INTENT(in) :: nproc     ! Number of processors in total
INTEGER, INTENT(in) :: n_procy   ! Number of processors N-S
INTEGER, INTENT(in) :: halo_i    ! Size of halo in i direction.
INTEGER, INTENT(in) :: halo_j    ! Size of halo in j direction.
INTEGER, INTENT(in) :: datastart(3)  ! First gridpoints held by this processor.
INTEGER, INTENT(in) :: g_row_length  ! global number of points on a row
INTEGER, INTENT(in) :: g_rows        ! global number of rows
INTEGER, INTENT(in) :: g_i_pe(1-halo_i:g_row_length+halo_i)  ! processor on my
! processor-row holding a given value in i direction
INTEGER, INTENT(in) :: g_j_pe(1-halo_j:g_rows+halo_j)        ! processor on my
! processor-column holding a given value in j direction
INTEGER, INTENT(in) :: pole_handling  ! How to treat the poles:
!   0 - no calculations at the poles
!   1 - poles in one point
!   2 - poles in all points
INTEGER, INTENT(in) :: proc_all_group ! Group id for all processors
INTEGER, INTENT(in) :: row_length    ! points in a row on this pe
INTEGER, INTENT(in) :: rows          ! rows on this pe
INTEGER, INTENT(in) :: model_domain  ! integer code for model domain

LOGICAL, INTENT(in) :: l_regular     ! false if variable resolution
LOGICAL, INTENT(in) :: l_sl_halo_reprod  ! if true then sl code bit repoducible
! with any sensible halo size
LOGICAL l_2dcomm
INTEGER size_2dcomm
INTEGER group_2dcomm

REAL, INTENT(in) ::                                                            &
  data_in (1-halo_i:dim_i_in+halo_i, 1-halo_j:dim_j_in+halo_j,                 &
  dim_k_in)                 ! data to be interpolated

REAL, INTENT(in) :: lambda_out (dim_i_out, dim_j_out,dim_k_in)
! Lambda co-ordinate of output data on input.
REAL, INTENT(in) ::  phi_out (dim_i_out, dim_j_out, dim_k_in)
! Phi co-ordinate of output data on input.

REAL, INTENT(in) :: weight_lambda_in (dim_i_out, dim_j_out, dim_k_out)
REAL, INTENT(in) :: weight_phi_in (dim_i_out, dim_j_out, dim_k_out)

INTEGER, INTENT(in) :: i_out_in (dim_i_out, dim_j_out, dim_k_out)
INTEGER, INTENT(in) :: j_out_in (dim_i_out, dim_j_out, dim_k_out)

LOGICAL, INTENT(in) :: at_extremity(4)  ! Indicates if this processor is at north,
! south, east or west of the processor grid

! Arguments with Intent OUT. ie: Output variables.

REAL, INTENT(out) ::  data_out (dim_i_out, dim_j_out, dim_k_out)
! data interpolated to desired locations.

! Local Variables.

! scalars

INTEGER :: i, j, k, kk      ! Loop indices
INTEGER :: j0, j1
INTEGER :: index
INTEGER :: errorstatus

INTEGER :: count

! arrays

INTEGER (KIND=integer32) :: i_out (dim_i_out, dim_j_out, dim_k_out)
INTEGER (KIND=integer32) :: j_out (dim_i_out, dim_j_out, dim_k_out)

REAL :: weight_lambda (dim_i_out, dim_j_out, dim_k_out)
REAL :: weight_phi (dim_i_out, dim_j_out, dim_k_out)

! Varibles applied in the "compute-on-demand" strategy

INTEGER :: ibase,ime

INTEGER                                                                        &
  irecv, my_imin, my_imax, dim_e_out                                           &
  , nsend, nrecv, info, len, itmp, sender                                      &
  , my_iminp, my_imaxp,my_jmin, my_jmax

INTEGER (KIND=integer32) :: i_store(dim_k_out*dim_i_out*dim_j_out,0:size_2dcomm)
INTEGER (KIND=integer32) :: j_store(dim_k_out*dim_i_out*dim_j_out,0:size_2dcomm)
INTEGER (KIND=integer32) :: k_store(dim_k_out*dim_i_out*dim_j_out,0:size_2dcomm)
INTEGER (KIND=integer32) :: i_out_e(dim_k_out*dim_i_out*dim_j_out)
INTEGER (KIND=integer32) :: j_out_e(dim_k_out*dim_i_out*dim_j_out)

INTEGER :: n_sendto(0:size_2dcomm)
INTEGER :: n_recvfrom(0:size_2dcomm)
INTEGER :: sp_send(0:size_2dcomm)
INTEGER :: sp_levels(0:size_2dcomm,dim_k_out)
INTEGER :: np_send(0:size_2dcomm)
INTEGER :: np_levels(0:size_2dcomm,dim_k_out)
REAL :: weight_lambda_e(dim_k_out*dim_i_out*dim_j_out)
REAL :: weight_phi_e(dim_k_out*dim_i_out*dim_j_out)
REAL :: data_out_e(dim_k_out*dim_i_out*dim_j_out)
REAL :: recv_data(dim_k_out*dim_i_out*dim_j_out,0:size_2dcomm)
REAL :: bcast_data(4*dim_k_out)


! Stuff for improved comms
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
INTEGER :: my_comm
INTEGER :: num_reqs
INTEGER :: request(0:size_2dcomm)
INTEGER :: send_stat(mpl_status_size,0:size_2dcomm)
INTEGER :: recv_stat(mpl_status_size)
INTEGER, SAVE :: mpl_send_type = imdi
INTEGER :: oldtypes(0:1)
INTEGER :: blockcounts(0:1)
INTEGER :: offsets(0:1)
INTEGER :: extent

TYPE sendrecv_type
  SEQUENCE
  INTEGER  :: i_out
  INTEGER  :: j_out
  REAL     :: weight_lambda
  REAL     :: weight_phi
END TYPE sendrecv_type

!!! may need to increase this size by g_rows, nproc-1
!     Type (sendrecv_type) :: send_array(dim_k_out*g_row_length*g_rows,0:nproc-1)
!     Type (sendrecv_type) :: recv_array(dim_k_out*g_row_length*g_rows,0:nproc-1)
TYPE (sendrecv_type) :: send_array(dim_k_out*g_row_length,0:size_2dcomm)
TYPE (sendrecv_type) :: recv_array(dim_k_out*g_row_length,0:size_2dcomm)

INTEGER :: dsm1x
INTEGER :: dsm1y

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle


! Functions: None

IF (lhook) CALL dr_hook('BI_LINEAR_H',zhook_in,zhook_handle)

IF (.NOT. l_2dcomm) THEN

  ! ----------------------------------------------------------------------
  ! Section 1.   Extend input data and r arrays to bigger area to
  !              allow interpolation to be done without having to re-do
  !              any end points.
  ! ----------------------------------------------------------------------

  ! No extension required in the parallel version - taken care of by
  ! swap_bounds.

  ! ----------------------------------------------------------------------
  ! Section 2.   For each output point find i,j so that the point on the
  !              output grid lies between i and i+1, j and j+1
  ! ----------------------------------------------------------------------

  ! i_out and j_out must be copied because they change in this subroutine
  DO k = 1, dim_k_out
    DO j = 1, dim_j_out
      DO i = 1, dim_i_out
        i_out(i,j,k) = i_out_in(i,j,k)
        j_out(i,j,k) = j_out_in(i,j,k)
      END DO
    END DO
  END DO

  IF ( n_procx > 1 ) THEN
    IF ( model_domain == mt_global ) THEN
      ! Send the points outside my region to the appropriate processor for
      ! interpolation. Only performed if the domain is decomposed in the
      ! i direction and not performed for LAM versions of the model.

      ! The first and last point I can interpolate in, based on available
      ! data on this processor

      my_imin = datastart(1) - halo_i + 2
      my_imax = datastart(1) + dim_i_out - 1 + halo_i - 2

      ! values for use in polar row to ensure pole is only calculated on one
      ! processor.
      my_iminp = datastart(1)
      my_imaxp = datastart(1)+dim_i_out-1
      IF (at_extremity(pwest) ) my_iminp = my_imin
      IF (at_extremity(peast) ) my_imaxp = my_imax

      ! The base processor on this row, and my address relative to that
      ! processor

      ibase = (me/n_procx) * n_procx
      ime = me - ibase

      DO i = 0, n_procx-1
        n_sendto(i) = 0
      END DO

      j0 = 1
      j1 = dim_j_out

      ! If the pole values not are going to be used, we we can save a
      ! large amount of communication by evaluating local dummys at the poles

      IF (pole_handling  ==  0) THEN

        IF (at_extremity(psouth)) THEN
          j0 = 2
          DO k = 1, dim_k_out
            DO i = 1, dim_i_out
              i_out(i,1,k) = i
            END DO
          END DO
        END IF

        IF (at_extremity(pnorth)) THEN
          j1 = dim_j_out-1
          DO k = 1, dim_k_out
            DO i = 1, dim_i_out
              i_out(i,dim_j_out,k) = i
            END DO
          END DO
        END IF

      END IF

      ! If all values along a polar row is evaluated in one point, we can
      ! save a large amount of communication by evaluating only one point
      ! correctly and then broadcast this value.

      IF (pole_handling  ==  1) THEN

        DO i = 0, n_procx-1
          sp_send(i) = 0
        END DO
        IF (at_extremity(psouth)) THEN
          j0 = 2
          DO k = 1, dim_k_out
            IF (i_out(1,1,k)  >=  my_iminp .AND.                               &
              i_out(1,1,k)  <=  my_imaxp) THEN
              i_out(1,1,k) = i_out(1,1,k) - datastart(1) + 1
              sp_send(ime) = sp_send(ime) + 1
              sp_levels(ime,sp_send(ime)) = k
              DO i = 2, dim_i_out
                i_out(i,1,k) = i ! i_out(i,1,k) - datastart(1) + 1
              END DO
            ELSE
              sender = g_i_pe(i_out(1,1,k))
              sp_send(sender) = sp_send(sender) + 1
              sp_levels(sender,sp_send(sender)) = k
              DO i = 1, dim_i_out
                i_out(i,1,k) = i
              END DO
            END IF
          END DO
        END IF

        DO i = 0, n_procx-1
          np_send(i) = 0
        END DO
        IF (at_extremity(pnorth)) THEN
          j1 = dim_j_out-1
          DO k = 1, dim_k_out
            IF (i_out(1,dim_j_out,k)  >=  my_iminp .AND.                       &
              i_out(1,dim_j_out,k)  <=  my_imaxp)  THEN
              np_send(ime) = np_send(ime) + 1
              np_levels(ime,np_send(ime)) = k
              i_out(1,dim_j_out,k) =                                           &
                i_out(1,dim_j_out,k) - datastart(1) + 1
              DO i = 2, dim_i_out
                i_out(i,dim_j_out,k) = i
              END DO
            ELSE
              sender = g_i_pe(i_out(1,dim_j_out,k))
              np_send(sender) = np_send(sender) + 1
              np_levels(sender,np_send(sender)) = k
              DO i = 1, dim_i_out
                i_out(i,dim_j_out,k) = i
              END DO
            END IF
          END DO
        END IF

      END IF

      IF ( l_sl_halo_reprod) THEN

        ! On the global boundaries, use i_out < 1 or i_out > g_row_length
        ! if that makes local computation possible. Not required when
        ! L_sl_halo_reprod is false is other logic ensures this is done.

        ! This code unsafe if applied at poles, where it isn't required.

        IF (at_extremity(pwest)) THEN
          DO k = 1, dim_k_out
            DO j = j0, j1
              DO i = 1, halo_i
                IF (i_out(i,j,k)  >   g_row_length-halo_i+2)                   &
                  i_out(i,j,k) = i_out(i,j,k) - g_row_length
              END DO
            END DO
          END DO
        END IF
        IF (at_extremity(peast)) THEN
          DO k = 1, dim_k_out
            DO j = j0, j1
              DO i = dim_i_out-halo_i+1, dim_i_out
                IF (i_out(i,j,k)  <   halo_i-1)                                &
                  i_out(i,j,k) = i_out(i,j,k) + g_row_length
              END DO
            END DO
          END DO
        END IF

      END IF ! on L_sl_halo_reprod

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, j, k, itmp,     &
!$OMP& irecv)
      DO k = 1, dim_k_out
        DO j = j0, j1
          DO i = 1, dim_i_out
            IF (i_out(i,j,k)  >=  my_imin .AND.                                &
              i_out(i,j,k)  <=  my_imax) THEN
              ! Process locally, so find the local destination
              i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
            ELSE
              ! Send to a remote processor, given by the array g_i_pe
              irecv = g_i_pe(i_out(i,j,k))

! prevent case of another thread modifying n_sendto proc point count
!$OMP CRITICAL 
              n_sendto(irecv) = n_sendto(irecv) + 1
              itmp = n_sendto(irecv)
!$OMP END CRITICAL 
              send_array(itmp,irecv) % i_out = i_out(i,j,k)
              send_array(itmp,irecv) % j_out = j_out(i,j,k)
              send_array(itmp,irecv) % weight_lambda = weight_lambda_in(i,j,k)
              send_array(itmp,irecv) % weight_phi = weight_phi_in(i,j,k)

              i_store(itmp,irecv) = i
              j_store(itmp,irecv) = j
              k_store(itmp,irecv) = k
!!! fix to stop out of memory address calls
              i_out(i,j,k)=i
              j_out(i,j,k)=j
            END IF
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      ! Counts can be distributed via an alltoall with the row communicator
      CALL mpl_alltoall(n_sendto,       1,    mpl_integer,                     &
        n_recvfrom,     1,    mpl_integer,                                     &
        group_2dcomm  , info)

      ! Get types setup if not done
      IF (mpl_send_type == imdi) THEN
        offsets    (0) = 0
        oldtypes   (0) = mpl_integer
        blockcounts(0) = 2

        CALL mpl_type_extent(mpl_integer, extent, info)

        offsets    (1) = 2 * extent
        oldtypes   (1) = mpl_real
        blockcounts(1) = 2

        CALL mpl_type_struct(2, blockcounts, offsets, oldtypes,                &
          mpl_send_type, info)
        CALL mpl_type_commit(mpl_send_type, info)
      END IF

      ! Send/Recv data in one hit using isend and recv
      num_reqs = 0
      DO i = 0,n_procx-1
        IF (n_sendto(i)  >   0) THEN
          CALL mpl_isend(send_array(1,i), n_sendto(i), mpl_send_type,          &
            i, 10, group_2dcomm  , request(num_reqs), info )
          num_reqs = num_reqs + 1
        END IF
      END DO

      DO i = 0,n_procx-1
        IF (n_recvfrom(i)  >   0) THEN
          CALL mpl_recv(recv_array(1,i), n_recvfrom(i), mpl_send_type,         &
            i, 10, group_2dcomm  , recv_stat, info )
        END IF
      END DO

      IF (num_reqs > 0) THEN
        CALL mpl_waitall(num_reqs, request, send_stat, info)
      END IF


      !      Else If ( L_regular ) then   ! for LAMs only
    ELSE      ! for LAMs only

      DO k = 1, dim_k_out
      DO j = 1, dim_j_out
      DO i = 1, dim_i_out
        i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
      END DO  ! i = 1, dim_i_out * dim_j_out * dim_k_out
      END DO
      END DO

    END IF ! model_domain == mt_Global

  END IF ! n_procx > 1

  ! ----------------------------------------------------------------------
  ! Section 3.   Perform required Interpolation.
  ! ----------------------------------------------------------------------

  ! DEPENDS ON: bi_linear
  CALL bi_linear (dim_i_out, dim_j_out, dim_k_out,                             &
    dim_i_in, dim_j_in, dim_k_in,                                              &
    halo_i, halo_j, data_in,                                                   &
    i_out, j_out, weight_lambda_in, weight_phi_in,                             &
    data_out)

  IF ( n_procx > 1 .AND. model_domain == mt_global ) THEN

    dim_e_out = 0
    DO i = 0, n_procx-1
      IF (n_recvfrom(i)  >   0) THEN
        DO j = 1, n_recvfrom(i)
          dim_e_out = dim_e_out + 1
          i_out_e(dim_e_out) = recv_array(j,i) % i_out -datastart(1)+1
          j_out_e(dim_e_out) = recv_array(j,i) % j_out
          weight_lambda_e(dim_e_out) = recv_array(j,i) % weight_lambda
          weight_phi_e(dim_e_out) = recv_array(j,i) % weight_phi
        END DO
      END IF
    END DO

    IF (dim_e_out  >   dim_k_out*dim_i_out*dim_j_out) THEN
      errorstatus = 10
      
      CALL ereport("Bi_linear_h", ErrorStatus,                                 &
        "over-writing due to dim_e_out size" )
    END IF

    IF (dim_e_out  >   0) THEN
      ! DEPENDS ON: bi_linear
      CALL bi_linear (dim_e_out, 1, 1,                                         &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        halo_i, halo_j, data_in,                                               &
        i_out_e, j_out_e,                                                      &
        weight_lambda_e, weight_phi_e,                                         &
        data_out_e)
    END IF


    nsend = 1
    DO i = 0, n_procx-1
      IF (n_recvfrom(i)  >   0) THEN
        len = n_recvfrom(i)
        CALL gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,                   &
          recv_data(1,ime), data_out_e(nsend))
        nsend = nsend + n_recvfrom(i)
      END IF
    END DO

    DO i = 0, n_procx-1
      IF (n_sendto(i)  >   0) THEN
        len = n_sendto(i)
        CALL gc_rrecv(40*(ibase+i+1)+me, len, ibase+i, info,                   &
          recv_data(1,i), data_out_e)
        DO j = 1, n_sendto(i)
          data_out(i_store(j,i), j_store(j,i),                                 &
            k_store(j,i)) = recv_data(j,i)
        END DO
      END IF
    END DO

    ! Now distribute the pole values.

    IF (pole_handling  ==  1) THEN

      IF (at_extremity(psouth)) THEN
        DO j = 0, n_procx-1
          IF (sp_send(j)  >   0) THEN
            DO kk = 1, sp_send(j)
              k = sp_levels(j,kk)
              bcast_data(kk) = data_out(1,1,k)
            END DO
            len = sp_send(j)
            CALL gcg_rbcast(201, len, ibase+j,                                 &
              group_2dcomm  , info, bcast_data)
            DO kk = 1, sp_send(j)
              k = sp_levels(j,kk)
              DO i = 1, dim_i_out
                data_out(i,1,k) = bcast_data(kk)
              END DO
            END DO
          END IF
        END DO
      END IF

      IF (at_extremity(pnorth)) THEN
        DO j = 0, n_procx-1
          IF (np_send(j)  >   0) THEN
            DO kk = 1, np_send(j)
              k = np_levels(j,kk)
              bcast_data(kk) = data_out(1,dim_j_out,k)
            END DO
            len = np_send(j)
            CALL gcg_rbcast(201, len, ibase+j,                                 &
              !    &                 proc_row_group, info, bcast_data)
              group_2dcomm  , info, bcast_data)
            DO kk = 1, np_send(j)
              k = np_levels(j,kk)
              DO i = 1, dim_i_out
                data_out(i,dim_j_out,k) = bcast_data(kk)
              END DO
            END DO
          END IF
        END DO
      END IF

    END IF

  END IF !  n_procx > 1 .and. model_domain == mt_Global

ELSE     ! L_2dcomm = TRUE

  ! ----------------------------------------------------------------------
  ! Section 1.   Extend input data and r arrays to bigger area to
  !              allow interpolation to be done without having to re-do
  !              any end points.
  ! ----------------------------------------------------------------------

  ! No extension required in the parallel version - taken care of by
  ! swap_bounds.

  ! ----------------------------------------------------------------------
  ! Section 2.   For each output point find i,j so that the point on the
  !              output grid lies between i and i+1, j and j+1
  ! ----------------------------------------------------------------------

  ! i_out and j_out must be copied because they change in this subroutine
  ! i_out is global , j_out is local
  dsm1y= datastart(2)-1
  dsm1x= datastart(1)-1
  DO k = 1, dim_k_out
    DO j = 1, dim_j_out
      DO i = 1, dim_i_out
        i_out(i,j,k) = i_out_in(i,j,k)
        j_out(i,j,k) = j_out_in(i,j,k) + dsm1y
        weight_lambda(i,j,k) = weight_lambda_in(i,j,k)
        weight_phi(i,j,k)    = weight_phi_in(i,j,k)
      END DO
    END DO
  END DO
  ! i_out is now global , j_out is now global

  !     if(nproc=1)then
  !       ! no action needed, point on processor
  !     elseif (n_procy >1 .and. n_procx==1) then
  !       ! only check on j_out
  !       ! do procedure
  !     elseif (n_procx >1 .and. n_procy==1) then
  !       ! only check on i_out
  !       ! do procedure
  !       ! equivalent of current code
  !     else   ! n_procx >1 .and. n_procy >1
  !       ! check i_out and then j_out
  !       ! do procedure
  !     endif

  ! procedure is - 1) find acceptable limits on this processor
  !                2) sort out what happens at pole (pole_handling) (global only)
  !                3) L_sl_halo_reprod section
  !                4a) if point on processor reset i_out to be local
  !                4b) if point off processor set up communication arrays
  !                5) send and receive communication arrays

  !       call bi_linear

  !     if(nproc /= 1)then
  !      ! setup local i_out_e etc arrays from communication arrays
  !      ! call bi_linear again
  !      ! send data back to processors
  !      ! sort out what happens at pole (pole_handling) (global only)
  !     end if

  ! For LAMs we need to check what we want to do if the point is at an edge
  ! and outside the acceptable range.

  IF ( n_procy>1 .OR. model_domain==mt_lam) THEN
    ! The first and last point I can interpolate in, based on available
    ! data on this processor
    my_jmin = datastart(2) - halo_j + 2
    my_jmax = datastart(2) + dim_j_out - 1 + halo_j - 2
  END IF

  IF ( n_procx>1 .OR. model_domain==mt_lam) THEN
    ! The first and last point I can interpolate in, based on available
    ! data on this processor
    my_imin = datastart(1) - halo_i + 2
    my_imax = datastart(1) + dim_i_out - 1 + halo_i - 2
  END IF

  IF (nproc > 1) THEN
    IF ( model_domain == mt_global .AND. pole_handling==1) THEN
      ! values for use in polar row to ensure pole is only calculated on one
      ! processor.
      my_iminp = datastart(1)
      my_imaxp = datastart(1)+dim_i_out-1
      IF (at_extremity(pwest) ) my_iminp = my_imin
      IF (at_extremity(peast) ) my_imaxp = my_imax
    END IF

    IF (model_domain == mt_lam) THEN
!!! check if need to reset weight_lambda/phi
      IF (at_extremity(pwest)) THEN
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              IF (i_out(i,j,k) < my_imin) THEN
                i_out(i,j,k)=my_imin
                weight_lambda(i,j,k)=0.0
              END IF
            END DO
          END DO
        END DO
      END IF
      IF (at_extremity(peast)) THEN
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              IF (i_out(i,j,k) > my_imax) THEN
                i_out(i,j,k)=my_imax
                weight_lambda(i,j,k)=1.0
              END IF
            END DO
          END DO
        END DO
      END IF
      IF (at_extremity(pnorth)) THEN
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              IF (j_out(i,j,k) > my_jmax) THEN
                j_out(i,j,k)=my_jmax
                weight_phi(i,j,k)=1.0
              END IF
            END DO
          END DO
        END DO
      END IF
      IF (at_extremity(psouth)) THEN
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              IF (j_out(i,j,k) < my_jmin) THEN
                j_out(i,j,k)=my_jmin
                weight_phi(i,j,k)=0.0
              END IF
            END DO
          END DO
        END DO
      END IF
    END IF  ! model_domain=mt_lam

    DO i = 0, nproc-1
      n_sendto(i) = 0
    END DO

    j0 = 1
    j1 = dim_j_out

    ! If the pole values not are going to be used, we we can save a
    ! large amount of communication by evaluating local dummys at the poles

    IF (model_domain==mt_global) THEN
      IF (pole_handling  ==  0) THEN
        IF (at_extremity(psouth)) THEN
          DO k = 1, dim_k_out
            DO i = 1, dim_i_out
              i_out(i,1,k) = datastart(1)+i-1
              j_out(i,1,k) = datastart(2)
            END DO
          END DO
        END IF
        IF (at_extremity(pnorth)) THEN
          DO k = 1, dim_k_out
            DO i = 1, dim_i_out
              i_out(i,dim_j_out,k) = datastart(1)+i-1
              j_out(i,dim_j_out,k) = datastart(2)+dim_j_out-1
            END DO
          END DO
        END IF
      END IF ! pole_handling=0

      ! If all values along a polar row is evaluated in one point, we can
      ! save a large amount of communication by evaluating only one point
      ! correctly and then broadcast this value.


      IF (pole_handling  ==  1) THEN
        IF (at_extremity(psouth)) THEN
          DO k = 1, dim_k_out
            DO i = 2, dim_i_out
              i_out(i,1,k) = i + datastart(1) - 1
              j_out(i,1,k) =  datastart(2)
            END DO
          END DO
        END IF

        IF (at_extremity(pnorth)) THEN
          DO k = 1, dim_k_out
            DO i = 2, dim_i_out
              i_out(i,dim_j_out,k) = i+datastart(1) - 1
              j_out(i,dim_j_out,k) = datastart(2)+dim_j_out-1
            END DO
          END DO
        END IF

      END IF      ! pole_handling=1
    END IF ! model_domain == mt_global
  END IF ! nproc>1

  IF (nproc > 1) THEN
    IF (model_domain==mt_global) THEN
      IF ( l_sl_halo_reprod) THEN

        ! On the global boundaries, use i_out < 1 or i_out > g_row_length
        ! if that makes local computation possible. Not required when
        ! L_sl_halo_reprod is false is other logic ensures this is done.

        ! This code unsafe if applied at poles, where it isn't required.
!!! need to check the comment above

        IF (at_extremity(pwest)) THEN
          DO k = 1, dim_k_out
            DO j = j0, j1
              DO i = 1, halo_i
                IF (i_out(i,j,k)  >   g_row_length-halo_i+2)                   &
                  i_out(i,j,k) = i_out(i,j,k) - g_row_length
              END DO
            END DO
          END DO
        END IF
        IF (at_extremity(peast)) THEN
          DO k = 1, dim_k_out
            DO j = j0, j1
              DO i = dim_i_out-halo_i+1, dim_i_out
                IF (i_out(i,j,k)  <   halo_i-1)                                &
                  i_out(i,j,k) = i_out(i,j,k) + g_row_length
              END DO
            END DO
          END DO
        END IF

      END IF ! on L_sl_halo_reprod
    END IF ! model_domain == mt_global
  END IF ! nproc>1

  IF (n_procx>1 .AND. n_procy==1) THEN
    DO k = 1, dim_k_out
      DO j = j0, j1
        DO i = 1, dim_i_out
          IF (i_out(i,j,k)  >=  my_imin .AND.                                  &
            i_out(i,j,k)  <=  my_imax) THEN
            i_out(i,j,k) = i_out(i,j,k) - dsm1x
            j_out(i,j,k) = j_out(i,j,k) - dsm1y
          ELSE
            irecv = g_i_pe(i_out(i,j,k))
            n_sendto(irecv) = n_sendto(irecv) + 1
            itmp = n_sendto(irecv)
            send_array(itmp,irecv) % i_out = i_out(i,j,k)
            send_array(itmp,irecv) % j_out = j_out(i,j,k)
            send_array(itmp,irecv) % weight_lambda = weight_lambda(i,j,k)
            send_array(itmp,irecv) % weight_phi = weight_phi(i,j,k)
            i_store(itmp,irecv) = i
            j_store(itmp,irecv) = j
            k_store(itmp,irecv) = k
            i_out(i,j,k)=i
            j_out(i,j,k)=j
          END IF
        END DO
      END DO
    END DO
  END IF ! n_procx>1 .and. n_procy==1

  IF (n_procy>1 .AND. n_procx==1) THEN
    DO k = 1, dim_k_out
      DO j = j0, j1
        DO i = 1, dim_i_out
          IF (j_out(i,j,k)  >=  my_jmin .AND.                                  &
            j_out(i,j,k)  <=  my_jmax) THEN
            i_out(i,j,k) = i_out(i,j,k) - dsm1x
            j_out(i,j,k) = j_out(i,j,k) - dsm1y
          ELSE
            irecv = g_j_pe(j_out(i,j,k))
            n_sendto(irecv) = n_sendto(irecv) + 1
            itmp = n_sendto(irecv)
            send_array(itmp,irecv) % i_out = i_out(i,j,k)
            send_array(itmp,irecv) % j_out = j_out(i,j,k)
            send_array(itmp,irecv) % weight_lambda = weight_lambda(i,j,k)
            send_array(itmp,irecv) % weight_phi = weight_phi(i,j,k)
            i_store(itmp,irecv) = i
            j_store(itmp,irecv) = j
            k_store(itmp,irecv) = k
            i_out(i,j,k)=i
            j_out(i,j,k)=j
          END IF
        END DO
      END DO
    END DO
  END IF  ! n_procy>1 .and. n_procx==1

  IF (n_procy>1 .AND. n_procx>1) THEN
    DO k = 1, dim_k_out
      DO j = j0, j1
        DO i = 1, dim_i_out
          IF ((j_out(i,j,k)  >=  my_jmin .AND.                                 &
            j_out(i,j,k)  <=  my_jmax) .AND.                                   &
            (i_out(i,j,k)  >=  my_imin .AND.                                   &
            i_out(i,j,k)  <=  my_imax) ) THEN
            i_out(i,j,k) = i_out(i,j,k) - dsm1x
            j_out(i,j,k) = j_out(i,j,k) - dsm1y
          ELSE
            !! need to check this
            irecv = g_i_pe(i_out(i,j,k)) +                                     &
              n_procx* g_j_pe(j_out(i,j,k))
            n_sendto(irecv) = n_sendto(irecv) + 1
            itmp = n_sendto(irecv)
            send_array(itmp,irecv) % i_out = i_out(i,j,k)
            send_array(itmp,irecv) % j_out = j_out(i,j,k)
            send_array(itmp,irecv) % weight_lambda = weight_lambda(i,j,k)
            send_array(itmp,irecv) % weight_phi = weight_phi(i,j,k)
            i_store(itmp,irecv) = i
            j_store(itmp,irecv) = j
            k_store(itmp,irecv) = k
            i_out(i,j,k)=i
            j_out(i,j,k)=j
          END IF
        END DO
      END DO
    END DO
  END IF ! n_procy>1 .and. n_procx>1

  IF ( nproc > 1 ) THEN
    !! original mpl section
    ! Counts can be distributed via an alltoall with the row communicator
!!! needs to be changed so does all processors not just row
    CALL mpl_alltoall(n_sendto,       1,    mpl_integer,                       &
      n_recvfrom,     1,    mpl_integer,                                       &
      proc_all_group, info)

    ! Get types setup if not done
    IF (mpl_send_type == imdi) THEN
      offsets    (0) = 0
      oldtypes   (0) = mpl_integer
      blockcounts(0) = 2

      CALL mpl_type_extent(mpl_integer, extent, info)

      offsets    (1) = 2 * extent
      oldtypes   (1) = mpl_real
      blockcounts(1) = 2

      CALL mpl_type_struct(2, blockcounts, offsets, oldtypes,                  &
        mpl_send_type, info)
      CALL mpl_type_commit(mpl_send_type, info)
    END IF

    ! Send/Recv data in one hit using isend and recv
    num_reqs = 0
    !       Do i = 0,n_procx-1
    DO i = 0,nproc-1
      IF (n_sendto(i)  >   0) THEN
        CALL mpl_isend(send_array(1,i), n_sendto(i), mpl_send_type,            &
          i, 10, proc_all_group, request(num_reqs), info )
        num_reqs = num_reqs + 1
      END IF
    END DO

    !       Do i = 0,n_procx-1
    DO i = 0,nproc-1
      IF (n_recvfrom(i)  >   0) THEN
        CALL mpl_recv(recv_array(1,i), n_recvfrom(i), mpl_send_type,           &
          i, 10, proc_all_group, recv_stat, info )
      END IF
    END DO

    IF (num_reqs > 0) THEN
      CALL mpl_waitall(num_reqs, request, send_stat, info)
    END IF
    !! end original mpl section
  END IF ! nproc >1

  ! ----------------------------------------------------------------------
  ! Section 3.   Perform required Interpolation.
  ! ----------------------------------------------------------------------

  ! DEPENDS ON: bi_linear
  CALL bi_linear (dim_i_out, dim_j_out, dim_k_out,                             &
    dim_i_in, dim_j_in, dim_k_in,                                              &
    halo_i, halo_j, data_in,                                                   &
    i_out, j_out, weight_lambda, weight_phi,                                   &
    data_out)

  IF ( nproc > 1 ) THEN

    dim_e_out = 0
    DO i = 0, nproc-1
      IF (n_recvfrom(i)  >   0) THEN
        DO j = 1, n_recvfrom(i)
          dim_e_out = dim_e_out + 1
          i_out_e(dim_e_out) = recv_array(j,i) % i_out -dsm1x
          j_out_e(dim_e_out) = recv_array(j,i) % j_out -dsm1y
          weight_lambda_e(dim_e_out) = recv_array(j,i) % weight_lambda
          weight_phi_e(dim_e_out) = recv_array(j,i) % weight_phi
        END DO
      END IF
    END DO

    IF (dim_e_out  >   dim_k_out*dim_i_out*dim_j_out) THEN
      errorstatus = 10
      
      CALL ereport("Bi_linear_h", ErrorStatus,                                 &
        "over-writing due to dim_e_out size" )
    END IF

    IF (dim_e_out  >   0) THEN
      ! DEPENDS ON: bi_linear
      CALL bi_linear (dim_e_out, 1, 1,                                         &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        halo_i, halo_j, data_in,                                               &
        i_out_e, j_out_e,                                                      &
        weight_lambda_e, weight_phi_e,                                         &
        data_out_e)
    END IF

    nsend = 1
    DO i = 0, nproc-1
      IF (n_recvfrom(i)  >   0) THEN
        len = n_recvfrom(i)
        CALL gc_rsend(40*(me+1)+i, len, i, info,                               &
          recv_data(1,me), data_out_e(nsend))
        nsend = nsend + n_recvfrom(i)
      END IF
    END DO
    DO i = 0, nproc-1
      IF (n_sendto(i)  >   0) THEN
        len = n_sendto(i)
        CALL gc_rrecv(40*(i+1)+me, len, i, info,                               &
          recv_data(1,i), data_out_e)
        DO j = 1, n_sendto(i)
          data_out(i_store(j,i), j_store(j,i),                                 &
            k_store(j,i)) = recv_data(j,i)
        END DO
      END IF
    END DO

    IF (model_domain== mt_global) THEN
      !! need to check this
      IF (pole_handling  ==  1) THEN

        IF (at_extremity(psouth)) THEN
          DO k = 1, dim_k_out
            DO i = 2, dim_i_out
              data_out(i,1,k) = data_out(1,1,k)
            END DO
          END DO
        END IF

        CALL gc_gsync(nproc,info)

        IF (at_extremity(pnorth)) THEN
          DO k = 1, dim_k_out
            DO i = 2, dim_i_out
              data_out(i,dim_j_out,k) = data_out(1,dim_j_out,k)
            END DO
          END DO
        END IF

      END IF  ! pole_handling==1
    END IF    !   model_domain == mt_Global
  END IF ! nproc >1

END IF    ! l_2dcomm

! End of routine.
IF (lhook) CALL dr_hook('BI_LINEAR_H',zhook_out,zhook_handle)
RETURN
END SUBROUTINE bi_linear_h

