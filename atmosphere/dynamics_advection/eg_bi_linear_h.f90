! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
SUBROUTINE eg_bi_linear_h(                                            &
                       data_in, lambda_out, phi_out, pnt_type,        &
                       dim_i_in, dim_j_in, dim_k_in,                  &
                       dim_i_out, dim_j_out, dim_k_out,               &
                       model_domain,                                  &
                       me, n_procx,                                   &
                       halo_i, halo_j, datastart,                     &
                       g_row_length, g_i_pe, at_extremity,            &
                       proc_row_group, data_out)

USE mpl, ONLY :                                                       &
  mpl_integer,                                                        &
  mpl_real,                                                           &
  mpl_status_size

USE locate_hdps_mod

USE missing_data_mod, ONLY:                                           &
    imdi

USE um_types, ONLY: integer32

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod, ONLY : ereport
USE UM_ParParams
IMPLICIT NONE
!
! Description:
!          Performs bi-linear horizontal interpolation of a field
!          defined on a u,v,w grid to another grid. Input data can be on a
!          sphere or be a rectangular box. Requested output points
!          lie inside the region defined for the input Data.
!  
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! Arguments with Intent IN. ie: Input variables.

INTEGER                                                                &
  dim_i_in                                                             &
              ! Dimension of Data_in in i direction.
, dim_j_in                                                             &
              ! Dimension of Data_in in j direction.
, dim_k_in                                                             &
              ! Dimension of Data_in in k direction.
, dim_i_out                                                            &
              ! Dimension of Data_out in i direction.
, dim_j_out                                                            &
              ! Dimension of Data_out in j direction.
, dim_k_out                                                            &
              ! Dimension of Data_out in j direction.
, pnt_type                                                             &
              ! grid point type
, me                                                                   &
              ! My processor number
, n_procx                                                              &
              ! Number of processors in longitude
, halo_i                                                               &
              ! Size of halo in i direction.
, halo_j                                                               &
              ! Size of halo in j direction.
, datastart(3)                                                         &
               ! First gridpoints held by this processor.
, g_row_length                                                         &
               ! global number of points on a row
, g_i_pe(1-halo_i:g_row_length+halo_i)                                 &
                                       ! processor on my
             ! processor-row holding a given value in i direction
, proc_row_group ! Group id for processors on the same row

INTEGER                                                                &
  model_domain     ! holds integer code for model domain

REAL                                                                   &
  data_in (1-halo_i:dim_i_in+halo_i, 1-halo_j:dim_j_in+halo_j,         &
           dim_k_in)                 ! data to be interpolated

REAL                                                                   &
  lambda_out (dim_i_out, dim_j_out,dim_k_in)                           &
                                              ! Lambda
                                     ! co-ordinate of
                                     ! output data on
                                     ! input.
, phi_out (dim_i_out, dim_j_out, dim_k_in)     ! Phi Co-ordinate
                                     ! of output data
                                     ! on input.

LOGICAL                                                                &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid


! Arguments with Intent OUT. ie: Output variables.

REAL                                                                   &
  data_out (dim_i_out, dim_j_out, dim_k_out) ! data interpolated
                                             ! to desired locatns.

! Local Variables.

! scalars

INTEGER                                                                &
  i, j, k, kk                                                          &
                       ! Loop indices
, j0, j1                                                               &
, INDEX

REAL                                                                   &
  recip_delta_lambda                                                   &
, recip_delta_phi                                                      &
, temp                                                                 &
, wk

! arrays

INTEGER (KIND=integer32) ::                                            &
  i_out (dim_i_out, dim_j_out, dim_k_out)                              &
, j_out (dim_i_out, dim_j_out, dim_k_out)

REAL                                                                   &
  weight_lambda (dim_i_out, dim_j_out, dim_k_out)                      &
, weight_phi (dim_i_out, dim_j_out, dim_k_out)

! Varibles applied in the "compute-on-demand" strategy

INTEGER                                                                &
  ime, ibase, irecv, my_imin, my_imax, dim_e_out                       &
, nsend, nrecv, info, comm_len, itmp, sender                           &
, my_iminp, my_imaxp

INTEGER, PARAMETER :: on_demand_size = 80

INTEGER (KIND=integer32) ::                                            &
  i_store(dim_k_out*dim_i_out*on_demand_size,0:n_procx-1)              &
, j_store(dim_k_out*dim_i_out*on_demand_size,0:n_procx-1)              &
, k_store(dim_k_out*dim_i_out*on_demand_size,0:n_procx-1)              &
, i_out_e(dim_k_out*dim_i_out*on_demand_size)                          &
, j_out_e(dim_k_out*dim_i_out*on_demand_size)

INTEGER                                                                &
  n_sendto(0:n_procx-1), n_recvfrom(0:n_procx-1)                       &
, sp_send(0:n_procx-1), sp_levels(0:n_procx-1,dim_k_out)               &
, np_send(0:n_procx-1), np_levels(0:n_procx-1,dim_k_out)

REAL                                                                   &
  weight_lambda_e(dim_k_out*dim_i_out*on_demand_size)                  &
, weight_phi_e(dim_k_out*dim_i_out*on_demand_size)                     &
, data_out_e(dim_k_out*dim_i_out*on_demand_size)                       &
, recv_data(dim_k_out*dim_i_out*on_demand_size,0:n_procx-1)            &
, bcast_data(4*dim_k_out)

INTEGER dsm1

! Stuff for improved comms 
INTEGER :: nsend_msg
INTEGER :: nrecv_msg
INTEGER :: extent 
INTEGER :: num_reqs 
INTEGER :: request(0:n_procx-1) 
INTEGER :: rec_request(0:n_procx-1) 
INTEGER :: send_stat(mpl_status_size,0:n_procx-1) 
INTEGER :: recv_stats(mpl_status_size,0:n_procx-1) 
INTEGER :: recv_stat(mpl_status_size) 
INTEGER, SAVE :: mpl_send_type = imdi 
INTEGER :: oldtypes(0:1) 
INTEGER :: blockcounts(0:1) 
INTEGER :: offsets(0:1) 

TYPE sendrecv_type 
  SEQUENCE 
  INTEGER  :: i_out 
  INTEGER  :: j_out 
  REAL     :: weight_lambda 
  REAL     :: weight_phi 
END TYPE 

TYPE (sendrecv_type) :: send_array(dim_k_out*g_row_length,0:n_procx) 
TYPE (sendrecv_type) :: recv_array(dim_k_out*g_row_length,0:n_procx) 

! Functions: None


! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook('EG_BI_LINEAR_H',zhook_in,zhook_handle)



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

! w points.

dsm1 = datastart(2)-1

CALL locate_hdps(i_out, j_out, weight_lambda, weight_phi,                   &
                 lambda_out, phi_out, dsm1,                                 &
                 pnt_type, dim_i_out, dim_j_out, dim_k_out)


IF ( n_procx > 1 ) THEN
  IF ( model_domain == mt_global ) THEN
! Send the points outside my region to the appropriate processor for
! interpolation. Only performed if the domain is decomposed in the
! i direction and not performed for LAM versions of the model.

! The first and last point I can interpolate in, based on available
! data on this processor

    my_imin = datastart(1) - halo_i + 2
    my_imax = datastart(1) + dim_i_out - 1 + halo_i - 2

!   values for use in polar row to ensure pole is only calculated on one
!   processor.
    my_iminp = datastart(1)
    my_imaxp = datastart(1)+dim_i_out-1
    IF(at_extremity(pwest) ) my_iminp = my_imin
    IF(at_extremity(peast) ) my_imaxp = my_imax

!   The base processor on this row, and my address relative to that
!   processor

    ibase = (me/n_procx) * n_procx
    ime = me - ibase

    DO i = 0, n_procx-1
      n_sendto(i) = 0
    END DO

    j0 = 1
    j1 = dim_j_out

  DO k = 1, dim_k_out
    DO j = j0, j1
      DO i = 1, dim_i_out
        IF (i_out(i,j,k)  >=  my_imin .AND.                           &
            i_out(i,j,k)  <=  my_imax) THEN
! Process locally, so find the local destination
          i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
        ELSE
! Send to a remote processor, given by the array g_i_pe
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
!!! fix to stop out of memory address calls
          i_out(i,j,k)=i
          j_out(i,j,k)=j
        END IF
      END DO
    END DO
  END DO

  ! Counts can be distributed via an alltoall with the row communicator
  CALL mpl_alltoall(n_sendto,       1,    mpl_integer,               &
                    n_recvfrom,     1,    mpl_integer,               &
                    proc_row_group, info)

  ! Get types setup if not done
  IF (mpl_send_type == imdi) THEN
    offsets    (0) = 0
    oldtypes   (0) = mpl_integer
    blockcounts(0) = 2

    CALL mpl_type_extent(mpl_integer, extent, info)

    offsets    (1) = 2 * extent
    oldtypes   (1) = mpl_real
    blockcounts(1) = 2

    CALL mpl_type_struct(2, blockcounts, offsets, oldtypes,          &
                         mpl_send_type, info)
    CALL mpl_type_commit(mpl_send_type, info)
  END IF

  ! Send/Recv data in one hit using isend and recv
  num_reqs = 0
  DO i = 0, n_procx-1
    IF (n_sendto(i)  >   0) THEN
      CALL mpl_isend(send_array(1,i), n_sendto(i), mpl_send_type,    &
                     i, 10, proc_row_group, request(num_reqs), info )
      num_reqs = num_reqs + 1
    END IF
  END DO

  DO i = 0, n_procx-1
    IF (n_recvfrom(i)  >   0) THEN
      CALL mpl_recv(recv_array(1,i), n_recvfrom(i), mpl_send_type, &
                    i, 10, proc_row_group, recv_stat, info )
    END IF
  END DO

  IF (num_reqs > 0) THEN
    CALL mpl_waitall(num_reqs, request, send_stat, info)
  END IF

ELSE

  DO k = 1, dim_k_out
    DO j = 1, dim_j_out
       DO i = 1, dim_i_out
        i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
      END DO
    END DO
  END DO

END IF ! model_domain == mt_Global
END IF ! n_procx > 1


! ----------------------------------------------------------------------
! Section 3.   Perform required Interpolation.
! ----------------------------------------------------------------------

! DEPENDS ON: bi_linear
CALL bi_linear (dim_i_out, dim_j_out, dim_k_out,                      &
                dim_i_in, dim_j_in, dim_k_in,                         &
                halo_i, halo_j, data_in,                              &
                i_out, j_out, weight_lambda, weight_phi,              &
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

  IF (dim_e_out  >   dim_k_out*dim_i_out*on_demand_size) THEN

     Call ereport("eg_bi_linear_h", dim_e_out,                        &
                  "over-writing due to dim_e_out size" )
  END IF

  IF (dim_e_out  >   0) THEN
! DEPENDS ON: bi_linear
    CALL bi_linear (dim_e_out, 1, 1,                                  &
                    dim_i_in, dim_j_in, dim_k_in,                     &
                    halo_i, halo_j, data_in,                          &
                    i_out_e, j_out_e,                                 &
                    weight_lambda_e, weight_phi_e,                    &
                    data_out_e)
  END IF


  nsend = 1
  nsend_msg = 0
  DO i = 0, n_procx-1
    IF (n_recvfrom(i)  >   0) THEN
      comm_len = n_recvfrom(i)
      CALL mpl_isend(data_out_e(nsend), comm_len, mpl_real, i, 40,    &
                     proc_row_group, request(nsend_msg), info)

      nsend_msg = nsend_msg + 1
      nsend = nsend + n_recvfrom(i)
    END IF
  END DO

  nrecv_msg = 0
  DO i = 0, n_procx-1
    IF (n_sendto(i)  >   0) THEN
      comm_len = n_sendto(i)
      CALL mpl_irecv(recv_data(1,i), comm_len, mpl_real, i, 40,       &
                     proc_row_group, rec_request(nrecv_msg), info)

      nrecv_msg = nrecv_msg + 1
    END IF
  END DO

! Potential to move this wait to end of routine.
  IF (nsend_msg > 0) THEN
    CALL mpl_waitall(nsend_msg, request, send_stat, info)
  END IF
 
! This wait could possibly utilise a waitany so copies can happen
! before all data is received
  IF (nrecv_msg > 0) THEN
    CALL mpl_waitall(nrecv_msg, rec_request, recv_stats, info)
    DO i = 0, n_procx-1
      IF (n_sendto(i)  >   0) THEN
        DO j = 1, n_sendto(i)
          data_out(i_store(j,i), j_store(j,i), k_store(j,i)) = recv_data(j,i)
        END DO
      END IF
    END DO
  END IF

END IF !  n_procx > 1 .and. model_domain == mt_Global

! End of routine.
IF (lhook) CALL dr_hook('EG_BI_LINEAR_H',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_bi_linear_h
