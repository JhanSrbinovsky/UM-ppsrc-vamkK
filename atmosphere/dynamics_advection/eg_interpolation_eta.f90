! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_interpolation_eta_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_interpolation_eta(                                      &
                          eta_in,  pnt_type,                          &
                          number_of_inputs,                           &
                          dim_i_in, dim_j_in, dim_k_in,               &
                          dim_j_in_w,                                 &
                          dim_i_out, dim_j_out, dim_k_out,            &
                          high_order_scheme, monotone_scheme,         &
                          model_domain, l_high, l_mono,               &
                          eta_out, lambda_out, phi_out,               &
                          me, n_proc, n_procx, n_procy,               &
                          halo_i, halo_j, g_row_length,               &
                          datastart, at_extremity, g_i_pe,            &
                          proc_row_group, proc_col_group,             &
                          halo_data_out_i, halo_data_out_j,           &
                          off_x, off_y,                               &
                          error_code,                                 &
                          data_in, data_out,                          &
                          delta_r_in,                                 &
                          integrity_key, integrity_input,             &
                          k_int_linear_in)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE integrity_mod

USE ereport_mod, ONLY : ereport
USE PrintStatus_mod
USE UM_ParParams

USE locate_hdps_mod

USE interp_grid_const_mod

USE mpl, ONLY:     &
    mpl_integer,   &
    mpl_real,      &
    mpl_status_size

USE missing_data_mod, ONLY: &
    imdi, rmdi

USE um_types,   ONLY: integer32
USE highos_mod, ONLY: cubicLagrange, quinticLagrange,                 &
                      ECMWF_quasiCubic, ECMWF_mono_quasiCubic,        &
                      hCubic_vLin, hCubic_vQuintic,                   &
                      hLag3_vHerm3_d2, hLag3_vHerm3_d4

IMPLICIT NONE
!
! Description:
!  
!          Performs interpolation of a field or fields defined on one
!          grid to another grid.
!          Requested output points must
!          lie inside the region defined for the input Data.
!          The number of fields interpolated is either 1,2 or 3 and is
!          controlled by the switch number_of_inputs. If only one field
!          is to be interpolated the fields Data_in2/3 and Data_out2/3
!          should be set to dummy arguments as they are not used.
!
! Method:
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

INTEGER                                                               &
  dim_i_in                                                            &
              ! Dimension of Data_in in i direction.
, dim_j_in                                                            &
              ! Dimension of Data_in in j direction.
, dim_j_in_w                                                          &
              ! Dimension of Data_in in j direction.
, dim_k_in                                                            &
              ! Dimension of Data_in in k direction.
, dim_i_out                                                           &
              ! Dimension of Data_out in i direction.
, dim_j_out                                                           &
              ! Dimension of Data_out in j direction.
, dim_k_out                                                           &
              ! Dimension of Data_out in k direction.
, me                                                                  &
              ! My processor number
, n_proc                                                              &
              ! Total number of processors
, n_procx                                                             &
              ! Number of processors in longitude
, n_procy                                                             &
              ! Number of processors in latitude
, halo_i                                                              &
              ! Size of halo in i direction.
, halo_j                                                              &
              ! Size of halo in j direction.
, off_x                                                               &
, off_y                                                               &
, halo_data_out_i                                                     &
                  ! size of data out halo in i direction
, halo_data_out_j                                                     &
                  ! size of data out halo in j direction
, proc_row_group                                                      &
                 ! Group id for processors on the same row
, proc_col_group                                                      &
                 ! Group id for processors on the same column
, g_row_length                                                        &
               ! global number of points on a row
, datastart(3)                                                        &
               ! First gridpoints held by this processor.
, g_i_pe(1-halo_i:g_row_length+halo_i) ! processor on my procr-row
                       ! holding a given value in i direction

INTEGER                                                               &
  pnt_type                                                            &
              ! Defines via an integer code the nature of the
              ! interpolation in terms of which grid the input
              ! Data is on. The codes are given in
              ! terms of a primary variable that would be held
              ! at that point and are u=1, v=2, w=3.
, high_order_scheme                                                   &
                     ! a code saying which high order scheme to
                     ! use.
, monotone_scheme                                                     &
                  ! a code saying which monotone scheme to use.
, number_of_inputs                                                    
                   !number of fields to interpolate.

INTEGER, OPTIONAL :: k_int_linear_in ! Level up to which linear 
                                     ! interpolation is used

LOGICAL                                                               &
  l_high                                                              &
                 ! True, if high order interpolation required.
, l_mono
                 ! True, if interpolation required to be monotone.
INTEGER                                                               &
  model_domain     ! holds integer code for model domain

REAL                                                                  &
  eta_in(dim_k_in) ! eta coordinate levels.

REAL                                                                  &
  data_in  (1-halo_i:dim_i_in+halo_i,                                 &
            1-halo_j:dim_j_in+halo_j, dim_k_in,                       &
            number_of_inputs )                                       
                                                ! data to be
                                                ! interpolated

REAL, OPTIONAL ::                                                     &
      delta_r_in (dim_i_in, dim_j_in, dim_k_in)
                                                ! Vertical
                                                ! layer thickness
                                                ! of input data.

REAL                                                                  &
  lambda_out (dim_i_out, dim_j_out, dim_k_out)                        &
                                                ! Lambda
                                                ! co-ordinate of
                                                ! output data on
                                                ! input.
, phi_out (dim_i_out, dim_j_out, dim_k_out)                           &
                                                ! Phi Co-ordinate
                                                ! of output data
                                                ! on input.
, eta_out (dim_i_out, dim_j_out, dim_k_out)     ! Vertical
                                                ! co-ordinate
                                                ! of output data.

LOGICAL                                                               &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid


! Arguments with Intent OUT. ie: Output variables.

REAL                                                                  &
     ! data interpolated to desired locations.
  data_out (1-halo_data_out_i:dim_i_out+halo_data_out_i,              &
            1-halo_data_out_j:dim_j_out+halo_data_out_j,              &
            dim_k_out, number_of_inputs)

INTEGER                                                               &
  error_code     ! Non-zero on exit if error detected.

! Local Variables.

! scalars

INTEGER   i, j, k, n, id, jd
INTEGER :: k_int_linear ! Level up to which linear interpolation is used
                        ! The value of this may be changed, optionally,
                        ! by k_int_linear_in

REAL      rdel, temp

! arrays

INTEGER (KIND=integer32) ::                                           &
  i_out (dim_i_out, dim_j_out, dim_k_out)                             &
, j_out (dim_i_out, dim_j_out, dim_k_out)                             &
, k_out (dim_i_out, dim_j_out, dim_k_out)


REAL                                                                  &
  ext_data(1-halo_i:dim_i_in+halo_i+1,1-halo_j:dim_j_in+halo_j,       &
           -1:dim_k_in+2,number_of_inputs)                            &
, data_out_intp (dim_i_out, dim_j_out, dim_k_out,                     &
                 number_of_inputs)                                    &
, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                     &
, weight_phi (dim_i_out, dim_j_out, dim_k_out)                        &
, coeff_z(-2:3, dim_i_out, dim_j_out, dim_k_out)


! Varibles applied in the "compute-on-demand" strategy

INTEGER                                                               &
  ime,                                                                &
  ibase,                                                              &
  irecv,                                                              &
  my_imin,                                                            &
  my_imax,                                                            &
  dim_e_out,                                                          &
  h_factor,                                                           &
  nsend,                                                              &
  nrecv,                                                              &
  info,                                                               &
  comm_len,                                                           &
  itmp,                                                               &
  j0,                                                                 &
  j1,                                                                 &
  kk,                                                                 &
  dsm1


INTEGER (KIND=integer32) ::                                           &
  i_store(dim_k_out*g_row_length,0:n_procx-1)                         &
, j_store(dim_k_out*g_row_length,0:n_procx-1)                         &
, k_store(dim_k_out*g_row_length,0:n_procx-1)                         &
, i_out_e(dim_k_out*g_row_length)                                     &
, j_out_e(dim_k_out*g_row_length)                                     &
, k_out_e(dim_k_out*g_row_length)

INTEGER                                                               &
  n_sendto(0:n_procx-1), n_recvfrom(0:n_procx-1)                      &
, sp_send(0:n_procx-1), sp_levels(0:n_procx-1,dim_k_out)              &
, np_send(0:n_procx-1), np_levels(0:n_procx-1,dim_k_out)

REAL                                                                  &
  weight_lambda_e(dim_k_out*g_row_length)                             &
, weight_phi_e(dim_k_out*g_row_length)                                &
, coeff_z_e(-2:3, dim_k_out*g_row_length)                             &
, eta_out_e(dim_k_out*g_row_length)                                   &
, data_out_intp_e(dim_k_out*g_row_length*number_of_inputs)            &
, send_data(number_of_inputs*dim_k_out*g_row_length,0:n_procx-1)      &
, recv_data(number_of_inputs*dim_k_out*g_row_length,0:n_procx-1)      &
, bcast_data(4*number_of_inputs*dim_k_out)

REAL :: x_off, y_off

! Varibles used for index vectors

INTEGER                                                               &
  ic,                                                                 &
  icnt,                                                               &
  inx(dim_i_out)

! pointer for strecthed grid data
  REAL, POINTER :: s_xi1(:), t_xi1(:), s_xi2(:), t_xi2(:),            &
                   q_xi1(:,:), q_xi2(:,:)

! Stuff for improved comms 
INTEGER :: nsend_msg
INTEGER :: nrecv_msg
INTEGER :: extent 
INTEGER :: num_reqs 
INTEGER :: request(0:n_procx-1) 
INTEGER :: rec_request(0:n_procx-1) 
INTEGER :: send_stat(mpl_status_size,0:n_procx-1) 
INTEGER :: recv_stat(mpl_status_size) 
INTEGER :: recv_stats(mpl_status_size,0:n_procx-1) 
INTEGER, SAVE :: mpl_send_type = imdi 
INTEGER :: oldtypes(0:1) 
INTEGER :: blockcounts(0:1) 
INTEGER :: offsets(0:1) 

TYPE sendrecv_type 
  SEQUENCE 
  INTEGER  :: i_out 
  INTEGER  :: j_out 
  INTEGER  :: k_out
  REAL     :: weight_lambda 
  REAL     :: weight_phi 
  REAL     :: eta_out 
END TYPE 

TYPE (sendrecv_type) :: send_array(dim_k_out*g_row_length,0:n_procx) 
TYPE (sendrecv_type) :: recv_array(dim_k_out*g_row_length,0:n_procx) 

INTEGER d_imin,d_imax, lev_ext
INTEGER k_found
INTEGER g_rows              !!!! Currently a dummy value = global_rows !!!!

CHARACTER (len=5), OPTIONAL :: integrity_key
INTEGER          , OPTIONAL :: integrity_input

LOGICAL                     :: l_cubic_interp

! Description: COMDECK containing the allowed
!              monotone scheme options
!
      INTEGER                                                           &
     &     triLinear                                                    &
     &,    mono_quasiCubic

      PARAMETER(                                                        &
     &     triLinear       = 1                                          &
     &,    mono_quasiCubic = 2 )

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_INTERPOLATION_ETA',zhook_in,zhook_handle)

IF(integrity_test.AND.PRESENT(integrity_key)) THEN

  IF(integrity_input.eq.1)                                            &
         CALL check_hash_m(data_in,SIZE(data_in),integrity_key)

END IF

! Check if linear interpolation is to be used on more than one layer
! near the lower boundary.
IF (PRESENT(k_int_linear_in)) THEN
  k_int_linear=k_int_linear_in
ELSE
  k_int_linear=1
END IF

! ---------------------------------------------------------------------
! Section 2.   Call appropriate routine to extend data.
!              Minimum amount of extending done to cope with highest
!              order interpolation scheme requested. This can leave
!              unset values in Ext_Data and Ext_r_in, but
!              these values will not be used.

!              Parallel version: No extension in r required.
! ---------------------------------------------------------------------

IF( high_order_scheme  ==  ecmwf_quasicubic      .OR.               &
    high_order_scheme  ==  ecmwf_mono_quasicubic .OR.               &
    monotone_scheme    ==  mono_quasicubic           ) THEN

   CALL ereport("eg_interpolation_eta", 1,                          &
                "Interpolation options not available within EG" )
END IF

l_cubic_interp = .FALSE.
IF( high_order_scheme == cubiclagrange     .OR.                      &
    high_order_scheme == hcubic_vquintic   .OR.                      &
    high_order_scheme == hcubic_vlin       .OR.                      &
    high_order_scheme == hLag3_vHerm3_d2   .OR.                      &
    high_order_scheme == hLag3_vHerm3_d4        ) THEN

   l_cubic_interp = .TRUE.
END IF

lev_ext = 0
IF ( high_order_scheme  ==  quinticlagrange .OR.                     &
     high_order_scheme  ==  hcubic_vquintic .OR.                     &
     high_order_scheme  ==  hLag3_vHerm3_d4      ) THEN

  lev_ext = 2

ELSE IF ( high_order_scheme  ==  cubiclagrange .OR.                  &
          high_order_scheme  ==  hLag3_vHerm3_d2    ) THEN

  lev_ext = 1

END IF

! Copy core data

  DO n = 1, number_of_inputs
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE) SHARED(dim_k_in,halo_j, &
!$OMP& halo_i,dim_i_in,dim_j_in,Ext_Data, Data_in,n) SCHEDULE(STATIC)
     DO k = 1, dim_k_in
        DO j = 1-halo_j, dim_j_in+halo_j
          DO i = 1-halo_i, dim_i_in+halo_i
            Ext_Data(i,j,k,n) = Data_in(i,j,k,n)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END DO

! Extend bottom

DO n = 1, number_of_inputs
   DO k = 1-lev_ext, 0
      DO j = 1-halo_j, dim_j_in+halo_j
         DO i = 1-halo_i, dim_i_in+halo_i
            Ext_Data(i,j,k,n) = rmdi
          END DO
      END DO
   END DO

! Extend Top

   DO k = 1+dim_k_in, dim_k_in + lev_ext
      DO j = 1-halo_j, dim_j_in+halo_j
         DO i = 1-halo_i, dim_i_in+halo_i
            Ext_Data (i,j,k,n) = rmdi
         END DO
      END DO
   END DO
END DO


! ---------------------------------------------------------------------
! Section 3.  For each output point find i,j,k so that the point on the
!             output grid lies between i and i+1, j and j+1, k and k+1
! ---------------------------------------------------------------------

SELECT CASE(pnt_type)
   CASE(fld_type_u)
     s_xi1 => sig_xi1_u
     t_xi1 => tau_xi1_u
     s_xi2 => sig_xi2_p
     t_xi2 => tau_xi2_p
     q_xi1 => q1_u
     q_xi2 => q2_p
   CASE(fld_type_v)
     s_xi1 => sig_xi1_p
     t_xi1 => tau_xi1_p
     s_xi2 => sig_xi2_v
     t_xi2 => tau_xi2_v
     q_xi1 => q1_p
     q_xi2 => q2_v
   CASE DEFAULT
     s_xi1 => sig_xi1_p
     t_xi1 => tau_xi1_p
     s_xi2 => sig_xi2_p
     t_xi2 => tau_xi2_p
     q_xi1 => q1_p
     q_xi2 => q2_p
END SELECT

dsm1 = datastart(2)-1

! Find i and j point.

CALL locate_hdps(i_out, j_out, weight_lambda, weight_phi,                   &
                 lambda_out, phi_out, dsm1,                                 &
                 pnt_type, dim_i_out, dim_j_out, dim_k_out)

! New search algorithm using the continuous nature of the search
! space to optimise

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(static)                        &
!$OMP&            PRIVATE(k,j,i,kk, k_found) SHARED(dim_k_out, dim_j_out,&
!$OMP&                    dim_i_out, dim_k_in, k_out, eta_out, eta_in)
DO k = 1, dim_k_out
   DO j = 1, dim_j_out
      DO i = 1, dim_i_out

! Start searching from the input data point
         k_found = k
         IF (eta_out(i,j,k) <= eta_in(k_found)) THEN

!search down
            DO kk = k_found, 1, -1
               k_found = kk
               IF (eta_out(i,j,k) > eta_in(kk)) EXIT
            END DO

         ELSE

! search up
            IF (k_found == dim_k_in) THEN
               k_found = k_found-1
            ELSE

               DO kk = k_found+1, dim_k_in-1
                  IF (eta_out(i,j,k) < eta_in(kk)) THEN
                     k_found = kk-1
                     EXIT
                  END IF
               END DO
            END IF
         END IF

        k_out(i,j,k) = k_found

      END DO
   END DO
END DO
!$OMP END PARALLEL DO


IF ( n_procx > 1 ) THEN
  IF ( model_domain == mt_global ) THEN

! Send the points outside my region to the appropriate processor for
! interpolation. Only performed if the domain is decomposed in the
! i direction.

! The first and last point I can interpolate in, based on available
! data on this processor (minus/plus one to avoid use of ge/le)

    h_factor = 2
    IF (high_order_scheme  ==  quinticlagrange ) THEN
      h_factor = 3
    END IF
    my_imin = datastart(1) - halo_i + h_factor - 1
    my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor + 1

! The base processor on this row, and my address relative to that
! processor

    ibase = (me/n_procx) * n_procx
    ime = me - ibase

    DO i = 0, n_procx-1
      n_sendto(i) = 0
    END DO

    j0 = 1
    j1 = dim_j_out

! And now decide where a point should be evaluated

    d_imin=my_imin-datastart(1) + 1
    d_imax=my_imax-datastart(1) + 1

    DO k = 1, dim_k_out
      DO j = j0, j1
        DO i = 1, dim_i_out
            i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
        END DO
      END DO
    END DO

    itmp=0
    DO k = 1, dim_k_out
      DO j = j0, j1
        DO i = 1, dim_i_out
          IF (i_out(i,j,k)  <=  d_imin .OR.                           &
              i_out(i,j,k)  >=  d_imax) THEN
            itmp=itmp+1
          END IF
        END DO
      END DO
    END DO

    DO k = 1, dim_k_out
      DO j = j0, j1
        icnt=0
        DO i = 1, dim_i_out
          IF (i_out(i,j,k)  <=  d_imin .OR.                           &
              i_out(i,j,k)  >=  d_imax) THEN
              icnt = icnt+1
              inx(icnt) = i
          END IF
        END DO

! Process locally, so find the local destination
        DO ic=1, icnt
           i = inx(ic)
           i_out(i,j,k) = i_out(i,j,k) + datastart(1) - 1

! Send to a remote processor, given by the array g_i_pe

!   CODE TO STOP BIT NON-REPRODUCIBILITY
           IF(i_out(i,j,k) >  g_row_length+halo_i-h_factor)THEN
              i_out(i,j,k) = i_out(i,j,k)-g_row_length
           END IF
           IF(i_out(i,j,k) <  1-halo_i+h_factor)THEN
              i_out(i,j,k) = i_out(i,j,k)+g_row_length
           END IF
!   END CODE TO STOP BIT NON-REPRODUCIBILITY

           irecv = g_i_pe(i_out(i,j,k))
           n_sendto(irecv) = n_sendto(irecv) + 1
           itmp = n_sendto(irecv)
           send_array(itmp, irecv) % i_out         = i_out(i,j,k)
           send_array(itmp, irecv) % j_out         = j_out(i,j,k)
           send_array(itmp, irecv) % k_out         = k_out(i,j,k)
           send_array(itmp, irecv) % weight_lambda = weight_lambda(i,j,k)
           send_array(itmp, irecv) % weight_phi    = weight_phi(i,j,k)
           send_array(itmp, irecv) % eta_out       = eta_out(i,j,k)
           i_store(itmp,irecv) = i
           j_store(itmp,irecv) = j
           k_store(itmp,irecv) = k
           i_out(i,j,k) = i

        END DO
      END DO
    END DO

! Send the points to be evaluated outside my region
! Counts can be distributed via an alltoall with the row communicator 
    CALL mpl_alltoall(n_sendto,       1,    mpl_integer,              & 
                      n_recvfrom,     1,    mpl_integer,              & 
                      proc_row_group, info) 

! Get types setup if not done 
    IF (mpl_send_type == imdi) THEN 
      offsets    (0) = 0 
      oldtypes   (0) = mpl_integer 
      blockcounts(0) = 3 

      CALL mpl_type_extent(mpl_integer, extent, info) 

      offsets    (1) = 3 * extent 
      oldtypes   (1) = mpl_real 
      blockcounts(1) = 3 

      CALL mpl_type_struct(2, blockcounts, offsets, oldtypes,         & 
                           mpl_send_type, info) 
      CALL mpl_type_commit(mpl_send_type, info) 
    END IF 

! Send/Recv data in one hit using isend and recv 
    num_reqs = 0 
    DO i = 0,n_procx-1 
      IF (n_sendto(i)  >   0) THEN 
        CALL mpl_isend(send_array(1,i), n_sendto(i), mpl_send_type, & 
                       i, 10, proc_row_group, request(num_reqs), info ) 
        num_reqs = num_reqs + 1 
      END IF 
    END DO 

    DO i = 0,n_procx-1 
      IF (n_recvfrom(i)  >   0) THEN 
        CALL mpl_recv(recv_array(1,i), n_recvfrom(i), mpl_send_type, & 
                      i, 10, proc_row_group, recv_stat, info ) 
      END IF 
    END DO 

    IF (num_reqs > 0) THEN 
      CALL mpl_waitall(num_reqs, request, send_stat, info) 
    END IF 

    dim_e_out = 0
    DO i = 0, n_procx-1
      IF (n_recvfrom(i)  >   0) THEN
        DO j = 1, n_recvfrom(i)
          dim_e_out = dim_e_out + 1
          i_out_e(dim_e_out)         = recv_array(j,i) % i_out       &
                                     - datastart(1) + 1
          j_out_e(dim_e_out)         = recv_array(j,i) % j_out
          k_out_e(dim_e_out)         = recv_array(j,i) % k_out
          weight_lambda_e(dim_e_out) = recv_array(j,i) % weight_lambda
          weight_phi_e(dim_e_out)    = recv_array(j,i) % weight_phi
          eta_out_e(dim_e_out)       = recv_array(j,i) % eta_out
        END DO
      END IF
    END DO

    IF (dim_e_out  >   dim_k_out*g_row_length) THEN
      CALL ereport("eg_interpolation_eta",dim_e_out,                  &
                  "warning over-writing due to dim_e_out size" )
    END IF

    IF (dim_e_out  >   0) THEN

! DEPENDS ON: eg_vert_weights_eta
       CALL eg_vert_weights_eta (                                     &
                      dim_k_in, dim_e_out, 1, 1,                      &
                      high_order_scheme, monotone_scheme,             &
                      model_domain, l_high, k_int_linear,             &
                      k_out_e, eta_in, eta_out_e,                     &
                      coeff_z_e)
    END IF

  ELSE  !  model_domain is NOT GLOBAL

      DO k = 1, dim_k_out
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
          END DO
        END DO
      END DO

  END IF ! model_domain == mt_Global
END IF ! n_procx > 1

! CODE TO STOP BIT NON-REPRODUCIBILITY
IF(model_domain == mt_global .AND. n_procx == 1)THEN
  h_factor = 2
  IF (high_order_scheme  ==  quinticlagrange ) THEN
    h_factor = 3
  END IF
  my_imin = datastart(1) - halo_i + h_factor - 1
  my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor + 1
  DO k = 1, dim_k_out
    DO j = 1,dim_j_out
      DO i = 1, dim_i_out
        IF(i_out(i,j,k) >= my_imax)THEN
          i_out(i,j,k)=i_out(i,j,k) - g_row_length
        END IF
        IF(i_out(i,j,k) <= my_imin )THEN
          i_out(i,j,k)=i_out(i,j,k)+g_row_length
        END IF
      END DO
    END DO
  END DO
END IF   !model_domain == 1 .and. n_procx == 1
!  END CODE TO STOP BIT NON-REPRODUCIBILITY

! calculate level vertical interpolation weights.

! DEPENDS ON: eg_vert_weights_eta
CALL eg_vert_weights_eta(                                             &
                      dim_k_in, dim_i_out, dim_j_out, dim_k_out,      &
                      high_order_scheme, monotone_scheme,             &
                      model_domain, l_high, k_int_linear,             &
                      k_out, eta_in, eta_out,                         &
                      coeff_z )

! ----------------------------------------------------------------------
! Section 4.   Perform required Interpolations.
! ----------------------------------------------------------------------

! Call high order scheme if required.

IF (l_cubic_interp) THEN
! DEPENDS ON: eg_cubic_lagrange
   CALL eg_cubic_lagrange(                                              &
        ext_data,                                                       &
        dim_i_in, dim_j_in, dim_k_in,                                   &
        dim_i_out, dim_j_out, dim_k_out,                                &
        halo_i, halo_j, number_of_inputs,                               &
        weight_lambda, weight_phi,                                      &
        s_xi1, s_xi2, t_xi1, t_xi2, q_xi1, q_xi2,                       &
        i_out, j_out, k_out,                                            &
        coeff_z, lev_ext,                                               &
        data_out_intp)

ELSE IF (high_order_scheme  ==  quinticlagrange) THEN

! DEPENDS ON: eg_quintic_lagrange
     CALL eg_quintic_lagrange(                                        &
                               ext_data,                              &
                               dim_i_in, dim_j_in, dim_k_in,          &
                               dim_i_out, dim_j_out, dim_k_out,       &
                               halo_i, halo_j, number_of_inputs,      &
                               weight_lambda, weight_phi,             &
                               i_out, j_out, k_out,                   &
                               coeff_z,                               &
                               data_out_intp )
END IF

IF( l_mono ) THEN
  IF (l_high ) THEN
! DEPENDS ON: mono_enforce
    CALL mono_enforce(                                                &
                      ext_data, number_of_inputs,                     &
                      dim_i_in, dim_j_in, dim_k_in,                   &
                      dim_i_out, dim_j_out, dim_k_out,                &
                      halo_i, halo_j,                                 &
                      i_out, j_out, k_out,                            &
                      data_out_intp)
  ELSE
! DEPENDS ON: eg_tri_linear
    CALL eg_tri_linear(                                               &
                         ext_data,                                    &
                         dim_i_in, dim_j_in, dim_k_in,                &
                         dim_i_out, dim_j_out, dim_k_out,             &
                         halo_i, halo_j, number_of_inputs,            &
                         weight_lambda, weight_phi,                   &
                         i_out, j_out, k_out,                         &
                         coeff_z,data_out_intp )
  END IF
END IF

! Repeat the interpolation procedure for the "compute-on-demand" points

IF (n_procx  >   1 .AND. model_domain  ==  mt_global) THEN

! Call high order scheme if required.

  IF (dim_e_out  >   0) THEN

  IF (l_cubic_interp) THEN
! DEPENDS ON: eg_cubic_lagrange
          CALL eg_cubic_lagrange(                                     &
            ext_data,                                                 &
            dim_i_in, dim_j_in, dim_k_in,                             &
            dim_e_out, 1, 1,                                          &
            halo_i, halo_j, number_of_inputs,                         &
            weight_lambda_e,weight_phi_e,                             &
            s_xi1, s_xi2, t_xi1, t_xi2, q_xi1, q_xi2,                 &
            i_out_e, j_out_e, k_out_e,                                &
            coeff_z_e, lev_ext,                                       &
            data_out_intp_e )

  ELSE IF (high_order_scheme  ==  quinticlagrange) THEN

! DEPENDS ON: eg_quintic_lagrange
      CALL eg_quintic_lagrange(                                       &
                             ext_data,                                &
                             dim_i_in, dim_j_in, dim_k_in,            &
                             dim_e_out, 1, 1,                         &
                             halo_i, halo_j, number_of_inputs,        &
                             weight_lambda_e,weight_phi_e,            &
                             i_out_e, j_out_e, k_out_e,               &
                             coeff_z_e,                               &
                             data_out_intp_e)

  END IF

  IF( l_mono ) THEN
    IF ( l_high ) THEN

! DEPENDS ON: mono_enforce
      CALL mono_enforce(                                              &
                      ext_data, number_of_inputs,                     &
                      dim_i_in, dim_j_in, dim_k_in,                   &
                      dim_e_out, 1, 1,                                &
                      halo_i, halo_j,                                 &
                      i_out_e, j_out_e, k_out_e,                      &
                      data_out_intp_e)
    ELSE

! DEPENDS ON: eg_tri_linear
      CALL eg_tri_linear(                                             &
                         ext_data,                                    &
                         dim_i_in, dim_j_in, dim_k_in,                &
                         dim_e_out, 1, 1,                             &
                         halo_i, halo_j, number_of_inputs,            &
                         weight_lambda_e,weight_phi_e,                &
                         i_out_e, j_out_e, k_out_e,                   &
                         coeff_z_e, data_out_intp_e)
    END IF
  END IF

END IF ! (dim_e_out  >   0)

! Time to return the points I have computed on demand.

nsend = 0
nsend_msg = 0

DO i = 0, n_procx-1
   IF (n_recvfrom(i)  >   0) THEN
        comm_len = n_recvfrom(i)*number_of_inputs
        DO k = 1, number_of_inputs
          DO j = 1, n_recvfrom(i)
            send_data((k-1)*n_recvfrom(i)+j,i)                        &
                 = data_out_intp_e(nsend+j+(k-1)*dim_e_out)
          END DO
        END DO

      CALL mpl_isend(send_data(1,i), comm_len, mpl_real, i, 40,       &
                     proc_row_group, request(nsend_msg), info)
       
      nsend_msg = nsend_msg + 1
      nsend = nsend + n_recvfrom(i)

    END IF
  END DO

  nrecv_msg = 0
  DO i = 0, n_procx-1
    IF (n_sendto(i)  >   0) THEN
      comm_len = n_sendto(i)*number_of_inputs

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

          DO n = 1, number_of_inputs
            DO j = 1, n_sendto(i)
              data_out_intp(i_store(j,i), j_store(j,i),                 &
                   k_store(j,i), n)                                     &
                   = recv_data((n-1)*n_sendto(i)+j,i)
            END DO
          END DO
      END IF
    END DO
  END IF

END IF ! (n_procx  >   1 .and. model_domain  ==  mt_Global)

! ----------------------------------------------------------------------
! Section 5.  Put interpolated field into Data_out.
! ----------------------------------------------------------------------

    DO n = 1, number_of_inputs
      DO k = 1, dim_k_out
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            data_out (i,j,k,n) = data_out_intp (i,j,k,n)
          END DO
        END DO
      END DO
    END DO 

IF (lhook) CALL dr_hook('EG_INTERPOLATION_ETA',zhook_out,zhook_handle)

END SUBROUTINE eg_interpolation_eta
END MODULE eg_interpolation_eta_mod
