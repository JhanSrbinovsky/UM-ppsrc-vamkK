! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine Calc_non_int_sl_theta


SUBROUTINE calc_non_int_sl_theta(                                              &
  w, theta, theta_lbcs,                                                        &
  r_theta_levels, r_rho_levels,                                                &
  check_bottom_levels,                                                         &
  interp_vertical_search_tol,                                                  &
  row_length, rows, model_levels,                                              &
  rimwidth, rimweights,                                                        &
  lenrim,lbc_size, lbc_start,                                                  &
  n_ext_fields,                                                                &
  lambda_rm, lambda_rp, phi_rm, phi_rp,                                        &
  recip_lambda_m, recip_lambda_0,                                              &
  recip_lambda_p, recip_lambda_p2,                                             &
  recip_phi_m, recip_phi_0,                                                    &
  recip_phi_p, recip_phi_p2,                                                   &
  i_out_in, j_out_in,                                                          &
  weight_lambda_in, weight_phi_in,                                             &
  model_domain, timestep, alpha_2,                                             &
  high_order_scheme_theta,                                                     &
  monotone_scheme_theta,                                                       &
  l_high_theta, l_mono_theta,                                                  &
  l_sl_halo_reprod, l_lbc_old,                                                 &
  l_regular, l_new_tdisc, cycleno,                                             &
  r_out, me, n_proc, n_procx, n_procy,                                         &
  off_x, off_y, halo_i, halo_j,                                                &
  datastart,  at_extremity,                                                    &
  g_row_length, g_rows,                                                        &
  proc_row_group, proc_col_group, proc_all_group,                              &
  g_i_pe, g_j_pe, l_2dcomm,                                                    &
  size_2dcomm, group_2dcomm, max_comm_size,                                    &
  theta_star, theta_np1,                                                       &
  error_code)

! Purpose:
!          Performs non-interpolating in the vertical semi-Lagrangian
!          advection of theta.

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
  mpl_real,                                                                    &
  mpl_integer,                                                                 &
  mpl_status_size

USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParParams
USE highos_mod,  ONLY: cubicLagrange, quinticLagrange, ECMWF_quasiCubic,       &
                       ECMWF_mono_quasiCubic, hCubic_vLin,                     &
                       hQuasiCubic_vQuintic, hCubic_vQuintic

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(in) :: me        ! My processor number
INTEGER, INTENT(in) :: n_procx   ! Number of processors E-W
INTEGER, INTENT(in) :: n_proc     ! Number of processors in total
INTEGER, INTENT(in) :: n_procy   ! Number of processors N-S
INTEGER, INTENT(in) :: halo_i    ! Size of halo in i direction.
INTEGER, INTENT(in) :: halo_j    ! Size of halo in j direction.
INTEGER, INTENT(in) :: off_x     ! Size of small halo in i direction.
INTEGER, INTENT(in) :: off_y     ! Size of small halo in j direction.
INTEGER, INTENT(in) :: row_length    ! points in a row on this pe
INTEGER, INTENT(in) :: rows          ! rows on this pe
INTEGER, INTENT(in) :: datastart(3)  ! First gridpoints held by this processor.
INTEGER, INTENT(in) :: g_row_length  ! global number of points on a row
INTEGER, INTENT(in) :: g_rows        ! global number of rows
INTEGER, INTENT(in) :: g_i_pe(1-halo_i:g_row_length+halo_i)  ! processor on my
! processor-row holding a given value in i direction
INTEGER, INTENT(in) :: g_j_pe(1-halo_j:g_rows+halo_j)        ! processor on my
! processor-column holding a given value in j direction
INTEGER, INTENT(in) :: proc_all_group ! Group id for all processors
INTEGER, INTENT(in) :: proc_row_group ! Group id for processors on thesame row
INTEGER, INTENT(in) :: proc_col_group ! Group id for processors on thesame column

LOGICAL :: l_2dcomm      ! T = doing 2D comms on demand
INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand
INTEGER :: max_comm_size ! error check size for comms on demand

INTEGER                                                                        &
  model_levels                                                                 &
                                ! Dimension of Data_in in k direction.
  , rimwidth                                                                   &
                                ! Width of boundaries in LBCs
  , rimweights(rimwidth)                                                       &
                                ! Weights to apply to the LBCs
  , lenrim                                                                     &
                                ! Size of single level of LBC data
  , lbc_size(4)                                                                &
                                ! Size of each side of LBC data
  , lbc_start(4)                                                               &
                                ! Start of each side in LBC data
  , cycleno

LOGICAL                                                                        &
  l_sl_halo_reprod                                                             &
                                ! if true then sl code bit repoducible with
                                ! any sensible halo size
  , l_regular                                                                  &
                                ! false if variable resolution
  , l_lbc_old                                                                  &
                                !  false for new lbc treatment
  , l_new_tdisc

INTEGER                                                                        &
  interp_vertical_search_tol                                                   &
                                !number of levels either side of
                                ! default level to search.
  , check_bottom_levels                                                        &
                                ! used in interpolation code, and is
                                ! the number of levels to check to see
                                ! if the departure point lies inside the
                                ! orography.
  , n_ext_fields    ! number of ext_data fields required, 1 for
! global model, 2 for LAM.

INTEGER                                                                        &
  model_domain     ! holds integer code for model domain

INTEGER                                                                        &
  high_order_scheme_theta                                                      &
                                ! a code saying which high order
                                ! scheme to use for theta.
  , monotone_scheme_theta ! a code saying which monotone
! scheme to use for theta.

LOGICAL                                                                        &
  l_high_theta                                                                 &
                                ! True, if high order interpolation required
                                !       for theta.
  , l_mono_theta   ! True, if interpolation required to be monotone
!       for theta.

REAL                                                                           &
  timestep                                                                     &
  , alpha_2

!VarRes horizontal co-ordinate information
REAL                                                                           &
  lambda_rm   ( 1-halo_i : row_length+halo_i )                                 &
  , lambda_rp   ( 1-halo_i : row_length+halo_i )                               &
  , phi_rm      ( 1-halo_i : row_length + halo_i                               &
  ,               1-halo_j : rows + halo_j )                                   &
  , phi_rp      ( 1-halo_i : row_length + halo_i                               &
  ,               1-halo_j : rows + halo_j )                                   &
  , recip_lambda_m (1-halo_i : row_length+halo_i)                              &
  , recip_lambda_0 (1-halo_i : row_length+halo_i)                              &
  , recip_lambda_p (1-halo_i : row_length+halo_i)                              &
  , recip_lambda_p2(1-halo_i : row_length+halo_i)                              &
  , recip_phi_m  ( 1-halo_i : row_length + halo_i                              &
  ,                1-halo_j : rows + halo_j )                                  &
  , recip_phi_0  ( 1-halo_i : row_length + halo_i                              &
  ,                1-halo_j : rows + halo_j )                                  &
  , recip_phi_p  ( 1-halo_i : row_length + halo_i                              &
  ,                1-halo_j : rows + halo_j )                                  &
  , recip_phi_p2 ( 1-halo_i : row_length + halo_i                              &
  ,                1-halo_j : rows + halo_j )

REAL                                                                           &
  r_rho_levels (1-halo_i:row_length+halo_i,                                    &
  1-halo_j:rows+halo_j, model_levels)                                          &
  , r_theta_levels (1-halo_i:row_length+halo_i,                                &
  1-halo_j:rows+halo_j, 0:model_levels)

REAL                                                                           &
  weight_lambda (row_length, rows, model_levels)                               &
  , weight_phi (row_length, rows, model_levels)
REAL                                                                           &
  weight_lambda_in (row_length, rows, model_levels)                            &
  , weight_phi_in (row_length, rows, model_levels)

INTEGER                                                                        &
  i_out_in (row_length, rows, model_levels)                                    &
  , j_out_in (row_length, rows, model_levels)

REAL                                                                           &
  theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                         &
  model_levels)                                                                &
  , theta_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                   &
  model_levels)                                                                &
  , w (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                           &
  0:model_levels)

REAL                                                                           &
  theta_lbcs(lenrim, model_levels)

LOGICAL                                                                        &
  at_extremity(4)  ! Indicates if this processor is at north,
! south, east or west of the processor grid

! Arguments with Intent IN/OUT. ie: Input and Output variables.
REAL                                                                           &
  r_out (row_length, rows, model_levels)        ! Vertical
! co-ordinate
! of output data.

! Arguments with Intent OUT. ie: Output variables.

REAL                                                                           &
  theta_star(1-off_x:row_length+off_x,                                         &
  1-off_y:rows+off_y, model_levels)

INTEGER                                                                        &
  error_code     ! Non-zero on exit if error detected.

! Local Variables.

! Description: COMDECK containing the allowed
!              monotone scheme options
!
      INTEGER                                                           &
     &     triLinear                                                    &
     &,    mono_quasiCubic

      PARAMETER(                                                        &
     &     triLinear       = 1                                          &
     &,    mono_quasiCubic = 2 )

INTEGER ixl(rows*row_length),jxl(rows*row_length),ijlc
INTEGER ixu(rows*row_length),jxu(rows*row_length),ijuc
INTEGER ij
! scalars

INTEGER                                                                        &
  i, j, k                                                                      &
                                ! Loop indices
  , index                                                                      &
  , lower_limit, upper_limit                                                   &
  , itype                                                                      &
  , number_of_inputs

INTEGER :: errorstatus

REAL                                                                           &
  r_below                                                                      &
  , r_above                                                                    &
  , r_here                                                                     &
  , r_belowi(row_length,rows)                                                  &
  , r_abovei(row_length,rows)                                                  &
  , r_herei(row_length,rows)                                                   &
  , r_mid_above                                                                &
  , r_mid_below                                                                &
  , r_mid_abovei(row_length,rows)                                              &
  , r_mid_belowi(row_length,rows)                                              &
  , max_mono                                                                   &
  , min_mono

LOGICAL                                                                        &
  l_vector                                                                     &
  , l_conserv

! arrays

LOGICAL                                                                        &
  l_continue(row_length, rows, model_levels)

INTEGER                                                                        &
  i_out (row_length, rows, model_levels)                                       &
  , j_out (row_length, rows, model_levels)

INTEGER                                                                        &
  depart_level(row_length, rows, model_levels)                                 &
                                ! closest model lev
                                ! to departure point
  , level_below(row_length, rows, model_levels) ! model level just
! below departure point

REAL                                                                           &
  w_star(row_length, rows, model_levels) ! vertical velocity
! required to move parcel
! from closest model level
! to departure point to arrival point model level.

REAL                                                                           &
  work (row_length, rows, model_levels)                                        &
  , work2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                   &
  model_levels)                                                                &
  , w_minus_wstar (row_length, rows, model_levels)                             &
  , ext_data (1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,                 &
  -1:model_levels+2,n_ext_fields )                                             &
  , r_d_s (row_length, rows)                                                   &
  , a_coeff_a (row_length, rows)                                               &
  , coeff_a                                                                    &
  , coeff_b                                                                    &
  , coeff_c                                                                    &
  , coeff_d                                                                    &
  , coeff_ai(row_length,rows)                                                  &
  , coeff_bi(row_length,rows)                                                  &
  , coeff_ci(row_length,rows)                                                  &
  , coeff_di(row_length,rows)                                                  &
  , theta_d_max(row_length, rows, model_levels)                                &
  , theta_d_min(row_length, rows, model_levels)

LOGICAL                                                                        &
  l_do_halos                                                                   &
                                ! update the halos?
  , l_do_boundaries   ! update the boundaries?

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
  INTEGER  :: k
  REAL     :: weight_lambda
  REAL     :: weight_phi
  REAL     :: r_out
  REAL     :: r_theta_levels
  REAL     :: w
END TYPE sendrecv_type

TYPE (sendrecv_type) :: send_array(max_comm_size,0:size_2dcomm)
TYPE (sendrecv_type) :: recv_array(max_comm_size,0:size_2dcomm)

! Variables applied in the "compute-on-demand" strategy

INTEGER                                                                        &
  ime, ibase, irecv, my_imin, my_imax, dim_e_out, h_factor                     &
  , nsend, nrecv, info, len, itmp, j0, j1                                      &
  , my_iminp, my_imaxp, my_jmin, my_jmax
LOGICAL l_continue_e
REAL r_d_s_e(max_comm_size)

INTEGER :: sp_send(0:n_procx-1)
INTEGER :: sp_levels(0:n_procx-1,model_levels)
INTEGER :: np_send(0:n_procx-1)
INTEGER :: np_levels(0:n_procx-1,model_levels)

INTEGER kk, sender
REAL ctmp1(2,model_levels)

INTEGER                                                                        &
  n_sendto(0:size_2dcomm), n_recvfrom(0:size_2dcomm)                           &
  , i_store(max_comm_size               ,0:size_2dcomm)                        &
  , j_store(max_comm_size               ,0:size_2dcomm)                        &
  , k_store(max_comm_size               ,0:size_2dcomm)                        &
  , i_out_e(max_comm_size               )                                      &
  , j_out_e(max_comm_size               )                                      &
  , depart_level_e(max_comm_size               )                               &
  , level_below_e(max_comm_size               )                                &
  , k_e(max_comm_size               )                                          &
  , isend_arr(3,max_comm_size               ,0:size_2dcomm)                    &
  , irecv_arr(3,max_comm_size               ,0:size_2dcomm)
REAL                                                                           &
  rsend_arr(5,max_comm_size               ,0:size_2dcomm)                      &
  , rrecv_arr(5,max_comm_size               ,0:size_2dcomm)                    &
  , weight_lambda_e(max_comm_size               )                              &
  , weight_phi_e(max_comm_size               )                                 &
  , coeff_a_e(max_comm_size               )                                    &
  , coeff_b_e(max_comm_size               )                                    &
  , coeff_c_e(max_comm_size               )                                    &
  , coeff_d_e(max_comm_size               )                                    &
  , r_out_e(max_comm_size               )                                      &
  , w_e(max_comm_size               )                                          &
  , w_star_e(max_comm_size               )                                     &
  , theta_d_max_e(max_comm_size               )                                &
  , theta_d_min_e(max_comm_size               )                                &
  , theta_star_e(max_comm_size               )                                 &
  , r_theta_levels_e(max_comm_size               )                             &
  , send_data(max_comm_size               ,0:size_2dcomm)                      &
  , recv_data(max_comm_size               ,0:size_2dcomm)                      &
  , rsend_arr2(2,max_comm_size               )                                 &
  , rrecv_arr2(2,max_comm_size               ,0:size_2dcomm)

REAL :: wrk

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle


! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Set some control variables.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CALC_NON_INT_SL_THETA',zhook_in,zhook_handle)
itype = 3 ! data at theta points
l_vector = .FALSE. ! field to be interpolated is not a horizontal
!                          vector component.
l_conserv = .FALSE. ! conservation not possible with this scheme
number_of_inputs = 1 ! only one field to interpolate

! Execute rest of routine only if error code is still zero.

IF (error_code  ==  0 ) THEN

  ! ----------------------------------------------------------------------
  ! Section 3.   For each output point find i,j so that the point on the
  !              output grid lies between i and i+1, j and j+1.
  ! ----------------------------------------------------------------------
  IF (.NOT.l_2dcomm) THEN

    ! i_out and j_out must be copied because they change in this subroutine
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i)                                    &
!$OMP& SHARED(rows,row_length,model_levels,i_out,j_out,i_out_in,j_out_in)      &
!$OMP& SHARED(weight_lambda,weight_phi,weight_lambda_in,weight_phi_in)
!$OMP DO SCHEDULE(DYNAMIC)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          i_out(i,j,k) = i_out_in(i,j,k)
          j_out(i,j,k) = j_out_in(i,j,k)
          weight_lambda(i,j,k) = weight_lambda_in(i,j,k)
          weight_phi(i,j,k)    = weight_phi_in(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (n_procx > 1 ) THEN
      IF (model_domain == mt_global ) THEN
        ! Send the points outside my region to the appropriate processor for
        ! interpolation. Only performed if the domain is decomposed in the
        ! i direction.

        ! The first and last point I can interpolate in, based on available
        ! data on this processor

        h_factor = 2
        IF ( high_order_scheme_theta  ==  quinticlagrange .AND.                &
          l_high_theta )                                                       &
          h_factor = 3
        my_imin = datastart(1) - halo_i + h_factor
        my_imax = datastart(1) + row_length - 1 + halo_i - h_factor

        ! values for use in polar row to ensure pole is only calculated on one
        ! processor.
        my_iminp = datastart(1)
        my_imaxp = datastart(1)+row_length-1
        IF (at_extremity(pwest) ) my_iminp = my_imin
        IF (at_extremity(peast) ) my_imaxp = my_imax

        ! The base processor on this row, and my address relative to that
        ! processor

        ibase = (me/n_procx) * n_procx
        ime = me - ibase

        DO i = 0, n_procx-1
          n_sendto(i) = 0
        END DO

        DO i = 0, n_procx-1
          sp_send(i) = 0
          sp_levels(i,:) = 0
        END DO
        IF (at_extremity(psouth)) THEN
          DO k = 1, model_levels
            IF (i_out(1,1,k)  >=  my_iminp .AND.                               &
              i_out(1,1,k)  <=  my_imaxp)  THEN
              i_out(1,1,k) = i_out(1,1,k) - datastart(1) + 1
              sp_send(ime) = sp_send(ime) + 1
              sp_levels(ime,sp_send(ime)) = k
              DO i = 2, row_length
                i_out(i,1,k) = i ! i_out(i,1,k) - datastart(1) + 1
              END DO
            ELSE
              sender = g_i_pe(i_out(1,1,k))
              sp_send(sender) = sp_send(sender) + 1
              sp_levels(sender,sp_send(sender)) = k
              DO i = 1, row_length
                i_out(i,1,k) = i
              END DO
            END IF
          END DO
        END IF
        DO i = 0, n_procx-1
          np_send(i) = 0
          np_levels(i,:) = 0
        END DO
        IF (at_extremity(pnorth)) THEN
          DO k = 1, model_levels
            IF (i_out(1,rows,k)  >=  my_iminp .AND.                            &
              i_out(1,rows,k)  <=  my_imaxp)  THEN
              np_send(ime) = np_send(ime) + 1
              np_levels(ime,np_send(ime)) = k
              i_out(1,rows,k) = i_out(1,rows,k)-datastart(1)+1
              DO i = 2, row_length
                i_out(i,rows,k) = i ! i_out(i,rows,k)-datastart(1)+1
              END DO
            ELSE
              sender = g_i_pe(i_out(1,rows,k))
              np_send(sender) = np_send(sender) + 1
              np_levels(sender,np_send(sender)) = k
              DO i = 1, row_length
                i_out(i,rows,k) = i
              END DO
            END IF
          END DO
        END IF

        j0 = 1
        j1 = rows
        IF (at_extremity(psouth)) j0 = 2
        IF (at_extremity(pnorth)) j1 = rows-1

        IF ( l_sl_halo_reprod) THEN

          ! On the global boundaries, use i_out < 1 or i_out > g_row_length
          ! if that makes local computation possible. Not required when
          ! L_sl_halo_reprod is false is other logic ensures this is done.

          ! This code unsafe if applied at poles, where it isn't required.

          IF (at_extremity(pwest)) THEN
            DO k = 1, model_levels
              DO j = j0, j1
                DO i = 1, halo_i
                  IF (i_out(i,j,k)  >   g_row_length-halo_i+h_factor)          &
                    i_out(i,j,k) = i_out(i,j,k) - g_row_length
                END DO
              END DO
            END DO
          END IF
          IF (at_extremity(peast)) THEN
            DO k = 1, model_levels
              DO j = j0, j1
                DO i = row_length-halo_i+1, row_length
                  IF (i_out(i,j,k)  <   halo_i-h_factor+1) THEN
                    i_out(i,j,k) = i_out(i,j,k) + g_row_length
                  END IF
                END DO
              END DO
            END DO
          END IF

        END IF ! on L_sl_halo_reprod

        DO k = 1, model_levels
          DO j = j0, j1
            DO i = 1, row_length
              IF (i_out(i,j,k)  >=  my_imin .AND.                              &
                i_out(i,j,k)  <=  my_imax) THEN
                ! Process locally, so find the local destination
                i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
              ELSE
                !     CODE TO STOP BIT NON-REPRODUCIBILITY
                IF (i_out(i,j,k) > g_row_length+halo_i-h_factor) THEN
                  i_out(i,j,k)=i_out(i,j,k)-g_row_length
                END IF
                IF (i_out(i,j,k) < 1-halo_i+h_factor) THEN
                  i_out(i,j,k)=i_out(i,j,k)+g_row_length
                END IF
                !     END CODE TO STOP BIT NON-REPRODUCIBILITY
                ! Send to a remote processor, given by the array g_i_pe
                irecv = g_i_pe(i_out(i,j,k))
                n_sendto(irecv) = n_sendto(irecv) + 1
                itmp = n_sendto(irecv)
                send_array(itmp, irecv) % i_out           = i_out(i,j,k)
                send_array(itmp, irecv) % j_out           = j_out(i,j,k)
                send_array(itmp, irecv) % k               = k
                send_array(itmp, irecv) % weight_lambda   =                    &
                  weight_lambda(i,j,k)
                send_array(itmp, irecv) % weight_phi      = weight_phi(i,j,k)
                send_array(itmp, irecv) % r_out           = r_out(i,j,k)
                send_array(itmp, irecv) % r_theta_levels  =                    &
                  r_theta_levels(i,j,k)
                send_array(itmp, irecv) % w               = w(i,j,k)

                i_store(itmp,irecv) = i
                j_store(itmp,irecv) = j
                k_store(itmp,irecv) = k
                i_out(i,j,k) = i
              END IF
            END DO
          END DO
        END DO

        ! Counts can be distributed via an alltoall with the row communicator
        CALL mpl_alltoall(n_sendto,       1,    mpl_integer,                   &
          n_recvfrom,     1,    mpl_integer,                                   &
          proc_row_group, info)

        ! Get types setup if not done
        IF (mpl_send_type == imdi) THEN
          offsets    (0) = 0
          oldtypes   (0) = mpl_integer
          blockcounts(0) = 3

          CALL mpl_type_extent(mpl_integer, extent, info)

          offsets    (1) = 3 * extent
          oldtypes   (1) = mpl_real
          blockcounts(1) = 5

          CALL mpl_type_struct(2, blockcounts, offsets, oldtypes,              &
            mpl_send_type, info)
          CALL mpl_type_commit(mpl_send_type, info)
        END IF

        ! Send/Recv data in one hit using isend and recv
        num_reqs = 0
        DO i = 0,n_procx-1
          IF (n_sendto(i)  >   0) THEN
            CALL mpl_isend(send_array(1,i), n_sendto(i), mpl_send_type,        &
              i, 10, proc_row_group, request(num_reqs), info )
            num_reqs = num_reqs + 1
          END IF
        END DO

        DO i = 0,n_procx-1
          IF (n_recvfrom(i)  >   0) THEN
            CALL mpl_recv(recv_array(1,i), n_recvfrom(i), mpl_send_type,       &
              i, 10, proc_row_group, recv_stat, info )
          END IF
        END DO

        IF (num_reqs > 0) THEN
          CALL mpl_waitall(num_reqs, request, send_stat, info)
        END IF


      ELSE   ! model is a type of LAM

        DO k = 1, model_levels
        DO j = 1, rows
        DO i = 1, row_length

          i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
        END DO  ! i = row_length * rows * model_levels
        END DO  ! j
        END DO  ! k

      END IF  !  model_domain == mt_Global

    END IF  ! n_procx > 1

    !     CODE TO STOP BIT NON-REPRODUCIBILITY
    IF (model_domain == mt_global .AND. n_procx == 1) THEN
      h_factor = 2
      IF (high_order_scheme_theta == quinticlagrange .AND.                     &
        l_high_theta) THEN
        h_factor = 3
      END IF
      my_imin = datastart(1) - halo_i + h_factor - 1
      my_imax =datastart(1) + row_length - 1 + halo_i - h_factor + 1
      DO k = 1, model_levels
        DO j = 1,rows
          DO i = 1, row_length
            IF (i_out(i,j,k) >= my_imax) THEN
              i_out(i,j,k)=i_out(i,j,k) - g_row_length
            END IF
            IF (i_out(i,j,k) <= my_imin ) THEN
              i_out(i,j,k)=i_out(i,j,k)+g_row_length
            END IF
          END DO
        END DO
      END DO
    END IF   !model_domain == 1 .and. n_procx == 1
    !     END CODE TO STOP BIT NON-REPRODUCIBILITY

  ELSE !l_2dcomm

    ! i_out and j_out must be copied because they change in this subroutine
    ! i_out is global , j_out is local
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i)                                    &
!$OMP& SHARED(rows,row_length,model_levels,i_out,j_out,i_out_in,j_out_in)      &
!$OMP& SHARED(weight_lambda,weight_phi,weight_lambda_in,weight_phi_in)         &
!$OMP& SHARED(datastart)
!$OMP DO SCHEDULE(DYNAMIC)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          i_out(i,j,k) = i_out_in(i,j,k)
          j_out(i,j,k) = j_out_in(i,j,k) + datastart(2) -1
          weight_lambda(i,j,k) = weight_lambda_in(i,j,k)
          weight_phi(i,j,k)    = weight_phi_in(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL
    ! i_out is now global , j_out is now global

    ! For LAMs we need to check what we want to do if the point is at an edge
    ! and outside the acceptable range.
    IF ( n_procy>1 .OR. model_domain==mt_lam) THEN
      ! The first and last point I can interpolate in, based on available
      ! data on this processor
      h_factor = 2
      IF (high_order_scheme_theta  ==  quinticlagrange .AND.                   &
        l_high_theta)    h_factor = 3
      !ajm         End If
      my_jmin = datastart(2) - halo_j + h_factor
      my_jmax = datastart(2) + rows - 1 + halo_j - h_factor
    END IF
    IF ( n_procx>1 .OR. model_domain==mt_lam) THEN
      ! The first and last point I can interpolate in, based on available
      ! data on this processor (minus/plus one to avoid use of ge/le)
      h_factor = 2
      IF (high_order_scheme_theta  ==  quinticlagrange .AND.                   &
        l_high_theta)   h_factor = 3
      my_imin = datastart(1) - halo_i + h_factor
      my_imax = datastart(1) + row_length - 1 + halo_i - h_factor
    END IF
    IF (n_proc >1) THEN
      IF ( model_domain == mt_global ) THEN
        ! values for use in polar row to ensure pole is only calculated on one
        ! processor,
        my_iminp = datastart(1)
        my_imaxp = datastart(1)+row_length-1
        IF (at_extremity(pwest) ) my_iminp = my_imin
        IF (at_extremity(peast) ) my_imaxp = my_imax
      END IF

      IF (model_domain == mt_lam) THEN
!!! check if need to reset weight_lambda/phi
        IF (at_extremity(pwest)) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                IF (i_out(i,j,k) < my_imin) THEN
                  i_out(i,j,k)=my_imin
                  weight_lambda(i,j,k)=0.0
                END IF
              END DO
            END DO
          END DO
        END IF
        IF (at_extremity(peast)) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                IF (i_out(i,j,k) > my_imax) THEN
                  i_out(i,j,k)=my_imax
                  weight_lambda(i,j,k)=1.0
                END IF
              END DO
            END DO
          END DO
        END IF
        IF (at_extremity(pnorth)) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                IF (j_out(i,j,k) > my_jmax) THEN
                  j_out(i,j,k)=my_jmax
                  weight_phi(i,j,k)=1.0
                END IF
              END DO
            END DO
          END DO
        END IF
        IF (at_extremity(psouth)) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                IF (j_out(i,j,k) < my_jmin) THEN
                  j_out(i,j,k)=my_jmin
                  weight_phi(i,j,k)=0.0
                END IF
              END DO
            END DO
          END DO
        END IF
      END IF  ! model_domain=mt_lam

      ! The base processor on this row, and my address relative to that
      ! processor

      ibase = (me/n_procx) * n_procx
      ime = me - ibase

      DO i = 0, n_proc-1
        n_sendto(i) = 0
      END DO

      IF ( model_domain == mt_global ) THEN
        IF (at_extremity(psouth)) THEN
          DO k = 1, model_levels
            DO i = 2, row_length
              i_out(i,1,k) = i + datastart(1) - 1
              j_out(i,1,k) = datastart(2)
            END DO
          END DO
        END IF
        IF (at_extremity(pnorth)) THEN
          DO k = 1, model_levels
            DO i = 2, row_length
              i_out(i,rows,k) = i + datastart(1) - 1
              j_out(i,rows,k) = datastart(2) + rows - 1
            END DO
          END DO
        END IF

        IF ( l_sl_halo_reprod) THEN

          ! On the global boundaries, use i_out < 1 or i_out > g_row_length
          ! if that makes local computation possible. Not required when
          ! L_sl_halo_reprod is false is other logic ensures this is done.

          ! This code unsafe if applied at poles, where it isn't required.

          IF (at_extremity(pwest)) THEN
            DO k = 1, model_levels
              DO j = j0, j1
                DO i = 1, halo_i
                  IF (i_out(i,j,k)  >   g_row_length-halo_i+h_factor)          &
                    i_out(i,j,k) = i_out(i,j,k) - g_row_length
                END DO
              END DO
            END DO
          END IF
          IF (at_extremity(peast)) THEN
            DO k = 1, model_levels
              DO j = j0, j1
                DO i = row_length-halo_i+1, row_length
                  IF (i_out(i,j,k)  <   halo_i-h_factor+1) THEN
                    i_out(i,j,k) = i_out(i,j,k) + g_row_length
                  END IF
                END DO
              END DO
            END DO
          END IF

        END IF ! on L_sl_halo_reprod
      END IF ! model_domain==mt_global
    END IF ! n_proc>1

    ! And now decide where a point should be evaluated
    IF (n_procx>1 .AND. n_procy==1) THEN
      DO k = 1, model_levels
        DO j = 1,rows
          DO i = 1, row_length
            IF (i_out(i,j,k)  >=  my_imin .AND.                                &
              i_out(i,j,k)  <=  my_imax) THEN
              ! Process locally, so find the local destination
              i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
              j_out(i,j,k) = j_out(i,j,k) - datastart(2) + 1
            ELSE
              ! Send to a remote processor, given by the array g_i_pe
              !     CODE TO STOP BIT NON-REPRODUCIBILITY
              IF (i_out(i,j,k) >  g_row_length+halo_i-h_factor) THEN
                i_out(i,j,k)=i_out(i,j,k)-g_row_length
              END IF
              IF (i_out(i,j,k) <  1-halo_i+h_factor) THEN
                i_out(i,j,k)=i_out(i,j,k)+g_row_length
              END IF
              !     END CODE TO STOP BIT NON-REPRODUCIBILITY
              irecv = g_i_pe(i_out(i,j,k))
              n_sendto(irecv) = n_sendto(irecv) + 1
              itmp = n_sendto(irecv)
              send_array(itmp, irecv) % i_out           = i_out(i,j,k)
              send_array(itmp, irecv) % j_out           = j_out(i,j,k)
              send_array(itmp, irecv) % k               = k
              send_array(itmp, irecv) % weight_lambda   =                      &
                weight_lambda(i,j,k)
              send_array(itmp, irecv) % weight_phi      = weight_phi(i,j,k)
              send_array(itmp, irecv) % r_out           = r_out(i,j,k)
              send_array(itmp, irecv) % r_theta_levels  =                      &
                r_theta_levels(i,j,k)
              send_array(itmp, irecv) % w               = w(i,j,k)
              i_store(itmp,irecv) = i
              j_store(itmp,irecv) = j
              k_store(itmp,irecv) = k
              i_out(i,j,k) = i
              j_out(i,j,k) = j
            END IF
          END DO
        END DO
      END DO
    END IF ! n_procx>1 .and. n_procy==1

    IF (n_procy>1 .AND. n_procx==1) THEN
      DO k = 1, model_levels
        DO j = 1,rows
          DO i = 1, row_length
            IF (j_out(i,j,k)  >=  my_jmin .AND.                                &
              j_out(i,j,k)  <=  my_jmax) THEN
              ! Process locally, so find the local destination
              i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
              j_out(i,j,k) = j_out(i,j,k) - datastart(2) + 1
            ELSE
              ! Send to a remote processor, given by the array g_i_pe
              irecv = g_j_pe(j_out(i,j,k))
              n_sendto(irecv) = n_sendto(irecv) + 1
              itmp = n_sendto(irecv)
              send_array(itmp, irecv) % i_out           = i_out(i,j,k)
              send_array(itmp, irecv) % j_out           = j_out(i,j,k)
              send_array(itmp, irecv) % k               = k
              send_array(itmp, irecv) % weight_lambda   =                      &
                weight_lambda(i,j,k)
              send_array(itmp, irecv) % weight_phi      = weight_phi(i,j,k)
              send_array(itmp, irecv) % r_out           = r_out(i,j,k)
              send_array(itmp, irecv) % r_theta_levels  =                      &
                r_theta_levels(i,j,k)
              send_array(itmp, irecv) % w               = w(i,j,k)
              i_store(itmp,irecv) = i
              j_store(itmp,irecv) = j
              k_store(itmp,irecv) = k
              i_out(i,j,k) = i
              j_out(i,j,k) = j
            END IF
          END DO
        END DO
      END DO
    END IF ! n_procy>1 .and. n_procx==1

    IF (n_procx>1 .AND. n_procy>1) THEN
      DO k = 1, model_levels
        DO j = 1,rows
          DO i = 1, row_length
            IF (i_out(i,j,k)  >=  my_imin .AND.                                &
              i_out(i,j,k)  <=  my_imax .AND.                                  &
              j_out(i,j,k)  >=  my_jmin .AND.                                  &
              j_out(i,j,k)  <=  my_jmax) THEN
              ! Process locally, so find the local destination
              i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
              j_out(i,j,k) = j_out(i,j,k) - datastart(2) + 1
            ELSE
              ! Send to a remote processor, given by the array g_i_pe
              !     CODE TO STOP BIT NON-REPRODUCIBILITY
              IF (i_out(i,j,k) >  g_row_length+halo_i-h_factor) THEN
                i_out(i,j,k)=i_out(i,j,k)-g_row_length
              END IF
              IF (i_out(i,j,k) <  1-halo_i+h_factor) THEN
                i_out(i,j,k)=i_out(i,j,k)+g_row_length
              END IF
              !     END CODE TO STOP BIT NON-REPRODUCIBILITY
              irecv = g_i_pe(i_out(i,j,k)) +                                   &
                n_procx* g_j_pe(j_out(i,j,k))
              n_sendto(irecv) = n_sendto(irecv) + 1
              itmp = n_sendto(irecv)
              send_array(itmp, irecv) % i_out           = i_out(i,j,k)
              send_array(itmp, irecv) % j_out           = j_out(i,j,k)
              send_array(itmp, irecv) % k               = k
              send_array(itmp, irecv) % weight_lambda   =                      &
                weight_lambda(i,j,k)
              send_array(itmp, irecv) % weight_phi      = weight_phi(i,j,k)
              send_array(itmp, irecv) % r_out           = r_out(i,j,k)
              send_array(itmp, irecv) % r_theta_levels  =                      &
                r_theta_levels(i,j,k)
              send_array(itmp, irecv) % w               = w(i,j,k)
              i_store(itmp,irecv) = i
              j_store(itmp,irecv) = j
              k_store(itmp,irecv) = k
              i_out(i,j,k) = i
              j_out(i,j,k) = j
            END IF
          END DO
        END DO
      END DO
    END IF ! n_procx>1 .and. n_procy>1

    IF (n_proc > 1) THEN

      ! Counts can be distributed via an alltoall with the row communicator
      CALL mpl_alltoall(n_sendto,       1,    mpl_integer,                     &
        n_recvfrom,     1,    mpl_integer,                                     &
        proc_all_group, info)

      ! Get types setup if not done
      IF (mpl_send_type == imdi) THEN
        offsets    (0) = 0
        oldtypes   (0) = mpl_integer
        blockcounts(0) = 3

        CALL mpl_type_extent(mpl_integer, extent, info)

        offsets    (1) = 3 * extent
        oldtypes   (1) = mpl_real
        blockcounts(1) = 5

        CALL mpl_type_struct(2, blockcounts, offsets, oldtypes,                &
          mpl_send_type, info)
        CALL mpl_type_commit(mpl_send_type, info)
      END IF


      ! Send/Recv data in one hit using isend and recv
      num_reqs = 0
      DO i = 0,n_proc-1
        IF (n_sendto(i)  >   0) THEN
          CALL mpl_isend(send_array(1,i), n_sendto(i), mpl_send_type,          &
            i, 10, proc_all_group, request(num_reqs), info )
          num_reqs = num_reqs + 1
        END IF
      END DO

      DO i = 0,n_proc-1
        IF (n_recvfrom(i)  >   0) THEN
          CALL mpl_recv(recv_array(1,i), n_recvfrom(i), mpl_send_type,         &
            i, 10, proc_all_group, recv_stat, info )
        END IF
      END DO

      IF (num_reqs > 0) THEN
        CALL mpl_waitall(num_reqs, request, send_stat, info)
      END IF

    END IF ! n_proc >1

    !     CODE TO STOP BIT NON-REPRODUCIBILITY
    IF (model_domain == mt_global .AND. n_proc == 1) THEN
      h_factor = 2
      IF (high_order_scheme_theta == quinticlagrange .AND.                     &
        l_high_theta) THEN
        h_factor = 3
      END IF
      my_imin = datastart(1) - halo_i + h_factor - 1
      my_imax =datastart(1) + row_length - 1 + halo_i - h_factor + 1
      DO k = 1, model_levels
        DO j = 1,rows
          DO i = 1, row_length
            IF (i_out(i,j,k) >= my_imax) THEN
              i_out(i,j,k)=i_out(i,j,k) - g_row_length
            END IF
            IF (i_out(i,j,k) <= my_imin ) THEN
              i_out(i,j,k)=i_out(i,j,k)+g_row_length
            END IF
          END DO
        END DO
      END DO
    END IF   !model_domain == 1 .and. n_procx == 1
    !     END CODE TO STOP BIT NON-REPRODUCIBILITY

  END IF ! l_2dcomm

  ! ----------------------------------------------------------------------
  ! Section 4.   Find closest model level to departure point and
  !              calculate w_star.
  !              Limit trajectories so that they do not go below bottom
  !              data level.
  !              Find max/min values of theta that surround departure
  !              point. (only if L_theta_mono)
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  ! Section 4.1  W-star and trajectory limits for points on processor
  ! ----------------------------------------------------------------------

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED)                             &
!$OMP& PRIVATE(i, j, k, coeff_a, coeff_b, coeff_c, coeff_d, r_d_s, lower_limit, &
!$OMP& upper_limit, r_here, r_mid_above, r_above, ijlc, ijuc, r_below, coeff_ai,&
!$OMP& coeff_ci, coeff_di, coeff_bi, r_belowi, r_herei, r_abovei, r_mid_abovei, &
!$OMP& r_mid_belowi, r_mid_below, ixu, jxu, ixl, jxl, index)
  DO k = 1, model_levels

    ! calculate horizontal interpolation weights
    DO j = 1, rows
      DO i = 1, row_length
        depart_level(i,j,k) = k
        w_star(i,j,k) = 0.
        l_continue(i,j,k) = .TRUE.
      END DO
    END DO

    IF ( k  <=  check_bottom_levels ) THEN
      ! Perform check for below bottom data surface.

      DO j = 1, rows
        DO i = 1, row_length

          coeff_a = (1.-weight_lambda(i,j,k)) *                                &
            (1.-weight_phi(i,j,k))
          coeff_b = weight_lambda(i,j,k) *                                     &
            (1.-weight_phi(i,j,k))
          coeff_c = (1.-weight_lambda(i,j,k)) *                                &
            weight_phi(i,j,k)
          coeff_d = weight_lambda(i,j,k)*weight_phi(i,j,k)

          r_d_s(i,j) = coeff_a *                                               &
            r_theta_levels (i_out(i,j,k),j_out(i,j,k),1)                       &
            + coeff_b                                                          &
            * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k),1)                   &
            + coeff_c                                                          &
            * r_theta_levels (i_out(i,j,k),j_out(i,j,k)+1,1)                   &
            + coeff_d *                                                        &
            r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)+1,1)

          IF (r_out(i,j,k) <   r_d_s(i,j) ) THEN
            ! move trajectory up to lowest level of data.
            r_out(i,j,k) = r_d_s(i,j)
            ! Set logical switch to say it was done
            l_continue(i,j,k) = .FALSE.
            ! set w_star = w to turn off vertical adjustment term
            w_star(i,j,k) = w(i,j,k)
            ! depart_level = 1 to get bottom most value.
            depart_level(i,j,k) = 1
          END IF

        END DO
      END DO

    END IF

    ! Find k point.
    ! use search over restricted levels.
    ! Find level which is just below r_out value, min possible is
    ! max(1,k- interp_vertical_search_tol), max possible is
    ! min(model_levels, k+interp_vertical_search_tol)-1

    IF ( k  <=  check_bottom_levels ) THEN
      lower_limit = MAX(1,k - check_bottom_levels)
      upper_limit = MIN(model_levels,                                          &
        k + check_bottom_levels)
    ELSE
      lower_limit = MAX(1,k - interp_vertical_search_tol)
      upper_limit = MIN(model_levels,                                          &
        k + interp_vertical_search_tol)
    END IF

    IF (k  ==  1) THEN
      ! level 1 Only performs level check and upward search
      DO j = 1, rows
        DO i = 1, row_length
          coeff_a = (1.-weight_lambda(i,j,k)) *                                &
            (1.-weight_phi(i,j,k))
          coeff_b = weight_lambda(i,j,k) *                                     &
            (1.-weight_phi(i,j,k))
          coeff_c = (1.-weight_lambda(i,j,k)) *                                &
            weight_phi(i,j,k)
          coeff_d = weight_lambda(i,j,k)*weight_phi(i,j,k)

          IF (l_continue(i,j,k)) THEN
            index = k
            r_above =                                                          &
              coeff_a * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)  ,index+1)  &
              + coeff_b * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)  ,index+1)&
              + coeff_c * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)+1,index+1)&
              + coeff_d * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,index+1)

            r_here = r_d_s(i,j)
            r_mid_above= (r_above + r_here )/2.

            IF (r_out(i,j,k)  <   r_mid_above ) THEN
              ! No need to check r_below as bottom adjust has done this.
              ! set w_star

              w_star(i,j,k) = (r_theta_levels(i,j,k) - r_d_s(i,j))             &
                / timestep

            ELSE
              ! upward search
              DO WHILE (index  <   upper_limit .AND.                           &
                l_continue(i,j,k) )
                index = index + 1

                r_here = r_above
                r_mid_below= r_mid_above
                r_above =                                                      &
                  coeff_a * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)  ,index+1)&
                  + coeff_b * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)  ,index+1)&
                  + coeff_c * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)+1,index+1)&
                  + coeff_d * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,index+1)

                r_mid_above= (r_here + r_above)/2.
                IF (r_out(i,j,k)  >=  r_mid_below .AND.                        &
                  r_out(i,j,k)  <   r_mid_above ) THEN

                  depart_level(i,j,k) = index
                  w_star(i,j,k) = (r_theta_levels(i,j,k) - r_here)             &
                    / timestep

                  l_continue(i,j,k) = .FALSE.
                END IF
              END DO

            END IF ! on level search

          END IF ! on bottom adjustment
        END DO
      END DO

    ELSE IF (k  ==  model_levels) THEN
      ! AT upper limit data stays at top level so no work required.
      ! w_star and depart_level are already set
      DO j = 1, rows
        DO i = 1, row_length
          l_continue(i,j,k) = .FALSE.
        END DO
      END DO

    ELSE
      ! Interior levels
      ! Calculate level k value
      ijlc=0
      ijuc=0
      DO j = 1, rows
        DO i = 1, row_length
          coeff_a = (1.-weight_lambda(i,j,k)) *                                &
            (1.-weight_phi(i,j,k))
          coeff_b = weight_lambda(i,j,k) *                                     &
            (1.-weight_phi(i,j,k))
          coeff_c = (1.-weight_lambda(i,j,k)) *                                &
            weight_phi(i,j,k)
          coeff_d = weight_lambda(i,j,k)*weight_phi(i,j,k)

          IF (l_continue(i,j,k) ) THEN
            r_above =                                                          &
              coeff_a * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)  ,k+1)     &
              + coeff_b * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)  ,k+1)   &
              + coeff_c * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)+1,k+1)   &
              + coeff_d * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)+1,k+1)
            r_here =                                                           &
              coeff_a * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)  ,k)       &
              + coeff_b * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)  ,k)     &
              + coeff_c * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)+1,k)     &
              + coeff_d * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)+1,k)
            r_below =                                                          &
              coeff_a * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)  ,k-1)     &
              + coeff_b * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)  ,k-1)   &
              + coeff_c * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)+1,k-1)   &
              + coeff_d * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)+1,k-1)

            coeff_ai(i,j) = coeff_a
            coeff_bi(i,j) = coeff_b
            coeff_ci(i,j) = coeff_c
            coeff_di(i,j) = coeff_d
            r_belowi(i,j) = r_below
            r_herei(i,j) = r_here
            r_abovei(i,j) = r_above

            r_mid_abovei(i,j) = (r_above+r_here)/2.
            r_mid_belowi(i,j) = (r_below+r_here)/2.
          END IF
        END DO
      END DO

      DO j = 1, rows
        DO i = 1, row_length
          IF (l_continue(i,j,k) ) THEN
            r_here = r_herei(i,j)
            r_mid_above = r_mid_abovei(i,j)
            r_mid_below = r_mid_belowi(i,j)

            IF (r_out(i,j,k)  >   r_mid_above ) THEN
              ijlc=ijlc+1
              ixl(ijlc)=i
              jxl(ijlc)=j

            ELSE IF (r_out(i,j,k)  <   r_mid_below ) THEN
              ijuc=ijuc+1
              ixu(ijuc)=i
              jxu(ijuc)=j


            ELSE
              ! set w_star values as level index is correct
              w_star(i,j,k) = (r_theta_levels(i,j,k) -                         &
                r_here                                                         &
                ) / timestep
            END IF

          END IF ! end if on bottom adjustment
        END DO
      END DO

      DO ij=1,ijlc
        i=ixl(ij)
        j=jxl(ij)
        coeff_a = coeff_ai(i,j)
        coeff_b = coeff_bi(i,j)
        coeff_c = coeff_ci(i,j)
        coeff_d = coeff_di(i,j)
        r_below = r_belowi(i,j)
        r_here = r_herei(i,j)
        r_above = r_abovei(i,j)
        r_mid_above = r_mid_abovei(i,j)
        r_mid_below = r_mid_belowi(i,j)

        index = k
        ! upward search
        DO WHILE (index  <   upper_limit .AND.                                 &
          l_continue(i,j,k) )
          index = index + 1

          r_below     = r_here
          r_here      = r_above
          r_mid_below = r_mid_above
          r_above =                                                            &
            coeff_a * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)  ,index+1)    &
            + coeff_b * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)  ,index+1)  &
            + coeff_c * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)+1,index+1)  &
            + coeff_d * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,index+1)
          r_mid_above=(r_above+r_here)/2.
          IF (r_out(i,j,k)  >=  r_mid_below .AND.                              &
            r_out(i,j,k)  <   r_mid_above ) THEN

            depart_level(i,j,k) = index
            w_star(i,j,k) = (r_theta_levels(i,j,k) - r_here)                   &
              / timestep

            l_continue(i,j,k) = .FALSE.
          END IF
        END DO
      END DO

      DO ij=1,ijuc
        i=ixu(ij)
        j=jxu(ij)
        coeff_a = coeff_ai(i,j)
        coeff_b = coeff_bi(i,j)
        coeff_c = coeff_ci(i,j)
        coeff_d = coeff_di(i,j)
        r_below = r_belowi(i,j)
        r_here = r_herei(i,j)
        r_above = r_abovei(i,j)
        r_mid_above = r_mid_abovei(i,j)
        r_mid_below = r_mid_belowi(i,j)

        index = k

        DO WHILE (index  >   lower_limit .AND.                                 &
          l_continue(i,j,k))
          index = index - 1
          r_above = r_here
          r_here= r_below
          r_mid_above=r_mid_below
          r_below =                                                            &
            coeff_a * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)  ,index-1)    &
            + coeff_b * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)  ,index-1)  &
            + coeff_c * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)+1,index-1)  &
            + coeff_d * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,index-1)
          r_mid_below=(r_below+r_here)/2.

          IF (r_out(i,j,k)  >=  r_mid_below .AND.                              &
            r_out(i,j,k)  <   r_mid_above ) THEN

            depart_level(i,j,k) = index
            w_star(i,j,k) = (r_theta_levels(i,j,k) - r_here)                   &
              / timestep

            l_continue(i,j,k) = .FALSE.
          END IF
        END DO
      END DO

    END IF ! on which model level

    ! ----------------------------------------------------------------------
    ! Section 4.2  max/min values for points on processor
    ! ----------------------------------------------------------------------

    ! if monotonocity required then find level just below departure point
    ! depart_level must be either just below or just above desired level
    ! so work out which it is.
    IF (l_mono_theta ) THEN
      DO j = 1, rows
        DO i = 1, row_length
          coeff_a = (1.-weight_lambda(i,j,k)) *                                &
            (1.-weight_phi(i,j,k))
          coeff_b =      weight_lambda(i,j,k) *                                &
            (1.-weight_phi(i,j,k))
          coeff_c = (1.-weight_lambda(i,j,k)) *                                &
            weight_phi(i,j,k)
          coeff_d = weight_lambda(i,j,k)*weight_phi(i,j,k)

          r_below = coeff_a *                                                  &
            r_theta_levels(i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k))      &
            + coeff_b *                                                        &
            r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k))   &
            + coeff_c *                                                        &
            r_theta_levels (i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k))   &
            + coeff_d *                                                        &
            r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k))
          IF (r_below  <   r_out(i,j,k) .AND.                                  &
            depart_level(i,j,k)  /=  1) THEN
            level_below(i,j,k) = depart_level(i,j,k)
          ELSE
            level_below(i,j,k) = depart_level(i,j,k) -1
          END IF
        END DO
      END DO
    END IF  !  L_mono_theta

  END DO ! k = 1, model_levels
  
!$OMP END PARALLEL DO 

  IF (l_mono_theta ) THEN
    ! set up theta in ext_data
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          ext_data (i,j,k,1) = theta(i,j,k)
        END DO
      END DO
    END DO

    ! Ensure data above the top level is initialised for later max/min
    DO j = 1,rows
      DO i = 1,row_length
        ext_data (i,j,model_levels+1,1) = theta(i,j,model_levels)
      END DO
    END DO

    ! Ensure data below the bottom level is initialised for later max/min
    DO j = 1,rows
      DO i = 1,row_length
        ext_data (i,j,0,1) = theta(i,j,1)
      END DO
    END DO
    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(                                                          &
      ext_data(1-halo_i,1-halo_j,0,1),                                         &
      row_length, rows, model_levels+2,                                        &
      halo_i, halo_j, fld_type_p, l_vector)

    IF (model_domain == mt_lam) THEN

      IF (l_lbc_old) THEN
        ! Set data on edge processor edge haloes from lateral boundary data
        l_do_halos=.TRUE.
        l_do_boundaries=.FALSE.

        ! DEPENDS ON: set_lateral_boundaries
        CALL set_lateral_boundaries(                                           &
          row_length,rows,halo_i,halo_j,                                       &
          model_levels,fld_type_p,                                             &
          ext_data(1-halo_i,1-halo_j,1,1),                                     &
          lenrim,lbc_size,lbc_start,halo_i,halo_j,theta_lbcs,                  &
          rimwidth,rimwidth,rimweights,at_extremity,                           &
          l_do_boundaries,l_do_halos)

        ! Ensure data below the bottom level is initialised for later max/min
        ! DEPENDS ON: set_lateral_boundaries
        CALL set_lateral_boundaries(                                           &
          row_length,rows,halo_i,halo_j,                                       &
          1,fld_type_p,                                                        &
          ext_data(1-halo_i,1-halo_j,0,1),                                     &
          lenrim,lbc_size,lbc_start,halo_i,halo_j,theta_lbcs,                  &
          rimwidth,rimwidth,rimweights,at_extremity,                           &
          l_do_boundaries,l_do_halos)

        ! Ensure data above the top level is initialised for later max/min
        ! DEPENDS ON: set_lateral_boundaries
        CALL set_lateral_boundaries(                                           &
          row_length,rows,halo_i,halo_j,                                       &
          1,fld_type_p,                                                        &
          ext_data(1-halo_i,1-halo_j,model_levels+1,1),                        &
          lenrim,lbc_size,lbc_start,halo_i,halo_j,                             &
          theta_lbcs(1,model_levels),                                          &
          rimwidth,rimwidth,rimweights,at_extremity,                           &
          l_do_boundaries,l_do_halos)

      ELSE ! new lbcs just fill external halos

        ! DEPENDS ON: fill_external_halos
        CALL fill_external_halos(                                              &
          ext_data(1-halo_i,1-halo_j,0,1),                                     &
          row_length, rows, model_levels+2,                                    &
          halo_i, halo_j)

      END IF !  L_lbc_old

    END IF ! model_domain == mt_LAM

    ! calculate max/min theta values at departure points
    DO k = 1, model_levels - 1

      DO j = 1, rows
        DO i = 1, row_length
          theta_d_max(i,j,k) = MAX (                                           &
            ext_data(i_out(i,j,k),j_out(i,j,k),level_below(i,j,k),1),          &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k),level_below(i,j,k),1),        &
            ext_data(i_out(i,j,k),j_out(i,j,k)+1,level_below(i,j,k),1),        &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k)+1,level_below(i,j,k),1),      &
            ext_data(i_out(i,j,k),j_out(i,j,k),level_below(i,j,k)+1,1),        &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k),level_below(i,j,k)+1,1),      &
            ext_data(i_out(i,j,k),j_out(i,j,k)+1,level_below(i,j,k)+1,1),      &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k)+1,level_below(i,j,k)+1,1) )
          theta_d_min(i,j,k) = MIN (                                           &
            ext_data(i_out(i,j,k),j_out(i,j,k),level_below(i,j,k),1),          &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k),level_below(i,j,k),1),        &
            ext_data(i_out(i,j,k),j_out(i,j,k)+1,level_below(i,j,k),1),        &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k)+1,level_below(i,j,k),1),      &
            ext_data(i_out(i,j,k),j_out(i,j,k),level_below(i,j,k)+1,1),        &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k),level_below(i,j,k)+1,1),      &
            ext_data(i_out(i,j,k),j_out(i,j,k)+1,level_below(i,j,k)+1,1),      &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k)+1,level_below(i,j,k)+1,1) )
        END DO
      END DO

    END DO  !  k = 1, model_levels - 1

  END IF !  L_mono_theta

  ! ----------------------------------------------------------------------
  ! Section 4.3  W-star and trajectory limits for points requested from
  !              other processors
  ! ----------------------------------------------------------------------

  ! receive extra points and process them

  IF ((n_proc > 1  .AND. l_2dcomm)  .OR.                                       &
    (n_procx >1 .AND. model_domain == mt_global .AND.                          &
    .NOT.l_2dcomm) ) THEN

    IF (.NOT. l_2dcomm) THEN
      dim_e_out = 0
      DO i = 0, n_procx-1
        IF (n_recvfrom(i)  >   0) THEN
          DO j = 1, n_recvfrom(i)
            dim_e_out = dim_e_out + 1
            i_out_e(dim_e_out) =                                               &
              recv_array(j,i) % i_out - datastart(1) + 1
            j_out_e(dim_e_out) = recv_array(j,i) % j_out
            k_e(dim_e_out) = recv_array(j,i) % k
            weight_lambda_e(dim_e_out) = recv_array(j,i) % weight_lambda
            weight_phi_e(dim_e_out) = recv_array(j,i) % weight_phi
            r_out_e(dim_e_out) = recv_array(j,i) % r_out
            r_theta_levels_e(dim_e_out) = recv_array(j,i) % r_theta_levels
            w_e(dim_e_out) = recv_array(j,i) % w

          END DO
        END IF
      END DO

    ELSE  !l_2dcomm
      dim_e_out = 0
      DO i = 0, n_proc-1
        IF (n_recvfrom(i)  >   0) THEN
          DO j = 1, n_recvfrom(i)
            dim_e_out = dim_e_out + 1
            i_out_e(dim_e_out) =                                               &
              recv_array(j,i) % i_out - datastart(1) + 1
            j_out_e(dim_e_out) = recv_array(j,i) % j_out                       &
              - datastart(2) + 1
            k_e(dim_e_out) = recv_array(j,i) % k
            weight_lambda_e(dim_e_out) = recv_array(j,i) % weight_lambda
            weight_phi_e(dim_e_out) = recv_array(j,i) % weight_phi
            r_out_e(dim_e_out) = recv_array(j,i) % r_out
            r_theta_levels_e(dim_e_out) = recv_array(j,i) % r_theta_levels
            w_e(dim_e_out) = recv_array(j,i) % w
          END DO
        END IF
      END DO

    END IF ! L_2dcomm

    IF (dim_e_out  >  max_comm_size            ) THEN
      errorstatus = 10
      
      CALL ereport("calc_non_int_sl_theta", ErrorStatus,                       &
        "over-writing due to dim_e_out size" )
    END IF

    DO i = 1, dim_e_out

      ! calculate horizontal interpolation weights
      coeff_a_e(i) = (1.-weight_lambda_e(i)) *                                 &
        (1.-weight_phi_e(i))
      coeff_b_e(i) = weight_lambda_e(i) *                                      &
        (1.-weight_phi_e(i))
      coeff_c_e(i) = (1.-weight_lambda_e(i)) *                                 &
        weight_phi_e(i)
      coeff_d_e(i) = weight_lambda_e(i)*weight_phi_e(i)
      depart_level_e(i) = k_e(i)
      w_star_e(i) = 0.
      l_continue_e = .TRUE.

      IF ( k_e(i)  <=  check_bottom_levels ) THEN
        ! Perform check for below bottom data surface.

        r_d_s_e(i) =                                                           &
          coeff_a_e(i) * r_theta_levels(i_out_e(i)  ,j_out_e(i)  ,1)           &
          + coeff_b_e(i) * r_theta_levels(i_out_e(i)+1,j_out_e(i)  ,1)         &
          + coeff_c_e(i) * r_theta_levels(i_out_e(i)  ,j_out_e(i)+1,1)         &
          + coeff_d_e(i) * r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,1)

        IF (r_out_e(i) <   r_d_s_e(i) ) THEN
          ! move trajectory up to lowest level of data.
          r_out_e(i) = r_d_s_e(i)
          ! Set logical switch to say it was done
          l_continue_e = .FALSE.
          ! set w_star = w to turn off vertical adjustment term
          w_star_e(i) = w_e(i)
          ! depart_level = 1 to get bottom most value.
          depart_level_e(i) = 1
        END IF

      END IF

      ! Find k point.
      ! use search over restricted levels.
      ! Find level which is just below r_out value, min possible is
      ! max(1,k- interp_vertical_search_tol), max possible is
      ! min(model_levels, k+interp_vertical_search_tol)-1

      IF ( k_e(i)  <=  check_bottom_levels ) THEN
        lower_limit = MAX(1,k_e(i) - check_bottom_levels)
        upper_limit = MIN(model_levels,                                        &
          k_e(i) + check_bottom_levels)
      ELSE
        lower_limit = MAX(1,k_e(i)-interp_vertical_search_tol)
        upper_limit = MIN(model_levels,                                        &
          k_e(i) + interp_vertical_search_tol)
      END IF

      IF (k_e(i)  ==  1) THEN
        ! level 1 Only performs level check and upward search
        IF (l_continue_e) THEN
          index = k_e(i)
          r_above =                                                            &
            coeff_a_e(i) *r_theta_levels(i_out_e(i)  ,j_out_e(i)  ,index+1)    &
            + coeff_b_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)  ,index+1)  &
            + coeff_c_e(i) *r_theta_levels(i_out_e(i)  ,j_out_e(i)+1,index+1)  &
            + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index+1)
          r_here=r_d_s_e(i)
          r_mid_above= (r_above +r_here)/2.

          IF (r_out_e(i)  <   r_mid_above ) THEN
            ! No need to check r_below as bottom adjust has done this.
            ! set w_star
            w_star_e(i) = (r_theta_levels_e(i) - r_d_s_e(i) )                  &
              / timestep

          ELSE
            ! upward search
            DO WHILE (index  <   upper_limit .AND.                             &
              l_continue_e )
              index = index + 1

              r_here=r_above
              r_mid_below = r_mid_above
              r_above =                                                        &
                coeff_a_e(i) *r_theta_levels(i_out_e(i)  ,j_out_e(i)  ,index+1) &
                +coeff_b_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)  ,index+1) &
                +coeff_c_e(i) *r_theta_levels(i_out_e(i)  ,j_out_e(i)+1,index+1) &
                +coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index+1)
              r_mid_above=(r_here+r_above)/2.

              IF (r_out_e(i)  >=  r_mid_below .AND.                            &
                r_out_e(i)  <   r_mid_above ) THEN

                depart_level_e(i) = index
                w_star_e(i) = (r_theta_levels_e(i) - r_here )                  &
                  / timestep

                l_continue_e = .FALSE.
              END IF
            END DO

          END IF ! on level search

        END IF ! on bottom adjustment

      ELSE IF (k_e(i)  ==  model_levels) THEN
        ! AT upper limit data stays at top level so no work required.
        ! w_star and depart_level are already set
        l_continue_e = .FALSE.
        depart_level_e(i) = k_e(i)
        w_star_e(i) = 0.0

      ELSE
        ! Interior levels
        ! Calculate level k value
        IF (l_continue_e ) THEN
          index = k_e(i)
          r_above =                                                            &
            coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index+1)    &
            + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index+1)  &
            + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index+1)  &
            + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index+1)
          r_here =                                                             &
            coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index)      &
            + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index)    &
            + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index)    &
            + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index)
          r_below =                                                            &
            coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index-1)    &
            + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index-1)  &
            + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index-1)  &
            + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index-1)
          r_mid_above=(r_here+r_above)/2.
          r_mid_below=(r_here+r_below)/2.

          IF (r_out_e(i)  >   r_mid_above ) THEN
            ! upward search
            DO WHILE (index  <   upper_limit .AND.                             &
              l_continue_e )
              index = index + 1

              r_below = r_here
              r_here  = r_above
              r_mid_below = r_mid_above
              r_above =                                                        &
                coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index+1)&
                + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index+1)&
                + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index+1)&
                + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index+1)
              r_mid_above=(r_here+r_above)/2.

              IF (r_out_e(i)  >=  r_mid_below .AND.                            &
                r_out_e(i)  <   r_mid_above ) THEN

                depart_level_e(i) = index
                w_star_e(i) = (r_theta_levels_e(i) - r_here )                  &
                  / timestep

                l_continue_e = .FALSE.
              END IF
            END DO

          ELSE IF (r_out_e(i)  <   r_mid_below ) THEN

            DO WHILE (index  >   lower_limit .AND.                             &
              l_continue_e)
              index = index - 1
              r_above = r_here
              r_here  = r_below
              r_mid_above = r_mid_below
              r_below =                                                        &
                coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index-1)&
                + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index-1)&
                + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index-1)&
                + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index-1)
              r_mid_below=(r_here+r_below)/2.

              IF (r_out_e(i)  >=  r_mid_below .AND.                            &
                r_out_e(i)  <   r_mid_above ) THEN

                depart_level_e(i) = index
                w_star_e(i) = (r_theta_levels_e(i) - r_here)                   &
                  / timestep

                l_continue_e = .FALSE.
              END IF
            END DO

          ELSE
            ! set w_star values as level index is correct
            w_star_e(i) = (r_theta_levels_e(i) - r_here)                       &
              / timestep
          END IF

        END IF ! end if on bottom adjustment

      END IF ! on which model level

    END DO ! end loop over communication on demand points

    ! ----------------------------------------------------------------------
    ! Section 4.4  Max/min values for points requested by other processors
    ! ----------------------------------------------------------------------

    ! if monotonocity required then find level just below departure point
    ! depart_level must be either just below or just above desired level
    ! so work out which it is.
    IF (l_mono_theta ) THEN
      DO i = 1, dim_e_out
        r_below = coeff_a_e(i) *                                               &
          r_theta_levels(i_out_e(i),j_out_e(i),depart_level_e(i))              &
          + coeff_b_e(i) *                                                     &
          r_theta_levels (i_out_e(i)+1,j_out_e(i),depart_level_e(i))           &
          + coeff_c_e(i) *                                                     &
          r_theta_levels (i_out_e(i),j_out_e(i)+1,depart_level_e(i))           &
          + coeff_d_e(i) *                                                     &
          r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,depart_level_e(i))
        IF (r_below  <   r_out_e(i) .AND.                                      &
          depart_level_e(i)  /=  1) THEN
          level_below_e(i) = depart_level_e(i)
        ELSE
          level_below_e(i) = depart_level_e(i)  - 1
        END IF

        ! calculate max/min theta values at departure points
        theta_d_max_e(i) = MAX (                                               &
          ext_data(i_out_e(i),j_out_e(i),level_below_e(i),1),                  &
          ext_data(i_out_e(i)+1,j_out_e(i),level_below_e(i),1),                &
          ext_data(i_out_e(i),j_out_e(i)+1,level_below_e(i),1),                &
          ext_data(i_out_e(i)+1,j_out_e(i)+1,level_below_e(i),1),              &
          ext_data(i_out_e(i),j_out_e(i),level_below_e(i)+1,1),                &
          ext_data(i_out_e(i)+1,j_out_e(i),level_below_e(i)+1,1),              &
          ext_data(i_out_e(i),j_out_e(i)+1,level_below_e(i)+1,1),              &
          ext_data(i_out_e(i)+1,j_out_e(i)+1,level_below_e(i)+1,1) )
        theta_d_min_e(i) = MIN (                                               &
          ext_data(i_out_e(i),j_out_e(i),level_below_e(i),1),                  &
          ext_data(i_out_e(i)+1,j_out_e(i),level_below_e(i),1),                &
          ext_data(i_out_e(i),j_out_e(i)+1,level_below_e(i),1),                &
          ext_data(i_out_e(i)+1,j_out_e(i)+1,level_below_e(i),1),              &
          ext_data(i_out_e(i),j_out_e(i),level_below_e(i)+1,1),                &
          ext_data(i_out_e(i)+1,j_out_e(i),level_below_e(i)+1,1),              &
          ext_data(i_out_e(i),j_out_e(i)+1,level_below_e(i)+1,1),              &
          ext_data(i_out_e(i)+1,j_out_e(i)+1,level_below_e(i)+1,1) )

      END DO ! end loop over communication on demand points
    END IF

    ! ----------------------------------------------------------------------
    ! Section 4.5  send data requested back to other processors
    ! ----------------------------------------------------------------------

    IF (.NOT.l_2dcomm) THEN

      ! send answers back to processor that asked for them
      ! 1. w_star
      nsend = 1
      DO i = 0, n_procx-1
        IF (n_recvfrom(i)  >   0) THEN
          len = n_recvfrom(i)
          CALL gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,                 &
            recv_data(1,ime), w_star_e(nsend))
          nsend = nsend + n_recvfrom(i)
        END IF
      END DO

      CALL gcg_ssync(proc_row_group, info)

      DO i = 0, n_procx-1
        IF (n_sendto(i)  >   0) THEN
          len = n_sendto(i)
          CALL gc_rrecv(40*(ibase+i+1)+me, len, ibase+i, info,                 &
            recv_data(1,i), w_star_e)
          DO j = 1, n_sendto(i)
            w_star(i_store(j,i), j_store(j,i),                                 &
              k_store(j,i)) = recv_data(j,i)
          END DO
        END IF
      END DO

      ! Note that we only need to send levels 1-model_levels of W for polar
      ! distribution - hence the addressing of W being the first point of the
      ! level 1
      ! DEPENDS ON: cnislt_distribute_poles
      CALL cnislt_distribute_poles(                                            &
        sp_send, sp_levels, np_send, np_levels,                                &
        n_procx, ibase, proc_row_group,                                        &
        rows, row_length, model_levels,                                        &
        off_x, off_y, 0, 0,                                                    &
        ime,  at_extremity, .TRUE.,                                            &
        w(1-off_x, 1-off_y, 1) , w_star)

      ! 2. max/min theta
      IF (l_mono_theta) THEN
        DO i = 1, dim_e_out
          rsend_arr2(1,i) = theta_d_max_e(i)
          rsend_arr2(2,i) = theta_d_min_e(i)
        END DO

        nsend = 1
        DO i = 0, n_procx-1
          IF (n_recvfrom(i)  >   0) THEN
            len = 2 * n_recvfrom(i)
            CALL gc_rsend(50*(me+1)+ibase+i, len, ibase+i, info,               &
              rrecv_arr2(1,1,ime), rsend_arr2(1,nsend) )
            nsend = nsend + n_recvfrom(i)
          END IF
        END DO

        CALL gcg_ssync(proc_row_group, info)

        DO i = 0, n_procx-1
          IF (n_sendto(i)  >   0) THEN
            len = n_sendto(i) * 2
            CALL gc_rrecv(50*(ibase+i+1)+me, len, ibase+i, info,               &
              rrecv_arr2(1,1,i), rsend_arr2 )
            DO j = 1, n_sendto(i)
              theta_d_max(i_store(j,i), j_store(j,i),                          &
                k_store(j,i)) = rrecv_arr2(1,j,i)
              theta_d_min(i_store(j,i), j_store(j,i),                          &
                k_store(j,i)) = rrecv_arr2(2,j,i)
            END DO
          END IF
        END DO

        ! DEPENDS ON: cnislt_distribute_poles
        CALL cnislt_distribute_poles(                                          &
          sp_send, sp_levels, np_send, np_levels,                              &
          n_procx, ibase, proc_row_group,                                      &
          rows, row_length, model_levels,                                      &
          0, 0, 0, 0,                                                          &
          ime,  at_extremity, .FALSE.,                                         &
          theta_d_max, theta_d_min)

      END IF   ! on (L_mono_theta)


    ELSE !l_2dcomm
      ! send answers back to processor that asked for them
      ! 1. w_star
      nsend = 1
      DO i = 0, n_proc-1
        IF (n_recvfrom(i)  >   0) THEN
          len = n_recvfrom(i)
          CALL gc_rsend(40*(me+1)+i, len, i, info,                             &
            recv_data(1,me), w_star_e(nsend))
          nsend = nsend + n_recvfrom(i)
        END IF
      END DO

      CALL gcg_ssync(proc_all_group, info)

      DO i = 0, n_proc-1
        IF (n_sendto(i)  >   0) THEN
          len = n_sendto(i)
          CALL gc_rrecv(40*(i+1)+me, len, i, info,                             &
            recv_data(1,i), w_star_e)
          DO j = 1, n_sendto(i)
            w_star(i_store(j,i), j_store(j,i),                                 &
              k_store(j,i)) = recv_data(j,i)
          END DO
        END IF
      END DO

      ! Note that we only need to send levels 1-model_levels of W for polar
      ! distribution - hence the addressing of W being the first point of the
      ! level 1
      ! Now distribute the pole values.

      IF (model_domain == mt_global) THEN
        IF (at_extremity(psouth)) THEN
          DO k = 1, model_levels
            DO i = 2, row_length
              w_star(i,1,k) = w_star(1,1,k)
            END DO
          END DO
        END IF
        IF (at_extremity(pnorth)) THEN
          DO k = 1, model_levels
            DO i = 2, row_length
              w_star(i,rows,k) = w_star(1,rows,k)
            END DO
          END DO
        END IF
      END IF ! (model_domain  ==  mt_Global)

      ! 2. max/min theta
      IF (l_mono_theta) THEN
        DO i = 1, dim_e_out
          rsend_arr2(1,i) = theta_d_max_e(i)
          rsend_arr2(2,i) = theta_d_min_e(i)
        END DO

        nsend = 1
        DO i = 0, n_proc-1
          IF (n_recvfrom(i)  >   0) THEN
            len = 2 * n_recvfrom(i)
            CALL gc_rsend(50*(me+1)+i, len, i, info,                           &
              rrecv_arr2(1,1,me), rsend_arr2(1,nsend) )
            nsend = nsend + n_recvfrom(i)
          END IF
        END DO

        CALL gcg_ssync(proc_all_group, info)

        DO i = 0, n_proc-1
          IF (n_sendto(i)  >   0) THEN
            len = n_sendto(i) * 2
            CALL gc_rrecv(50*(i+1)+me, len, i, info,                           &
              rrecv_arr2(1,1,i), rsend_arr2 )
            DO j = 1, n_sendto(i)
              theta_d_max(i_store(j,i), j_store(j,i),                          &
                k_store(j,i)) = rrecv_arr2(1,j,i)
              theta_d_min(i_store(j,i), j_store(j,i),                          &
                k_store(j,i)) = rrecv_arr2(2,j,i)
            END DO
          END IF
        END DO


        IF (model_domain == mt_global) THEN
          IF (at_extremity(psouth)) THEN
            DO k = 1, model_levels
              DO i = 2, row_length
                theta_d_max(i,1   ,k) = theta_d_max(1,1   ,k)
                theta_d_min(i,1   ,k) = theta_d_min(1,1   ,k)
              END DO
            END DO
          END IF
          IF (at_extremity(pnorth)) THEN
            DO k = 1, model_levels
              DO i = 2, row_length
                theta_d_max(i,rows,k) = theta_d_max(1,rows,k)
                theta_d_min(i,rows,k) = theta_d_min(1,rows,k)
              END DO
            END DO
          END IF
        END IF ! (model_domain  ==  mt_Global)
      END IF   ! on (L_mono_theta)

    END IF  ! l_2dcomm

  END IF  ! on nproc>1 and l_2dcomm    or
  ! n_procx>1 .and. model_domain==mt_Global .and.not.l_2dcomm

  ! ----------------------------------------------------------------------
  ! Section 5.    Calculate Eulerian vertical advection term.
  ! ----------------------------------------------------------------------

  k = 1
  DO j = 1, rows
    DO i = 1, row_length
      ! CFL=(w-w*) * dt /dz
      ! CFL limit is 1
      w_minus_wstar(i,j,k) = timestep *                                        &
        ( w(i,j,k) - w_star(i,j,k) )                                           &
        / (r_theta_levels(i,j,k+1) -                                           &
        r_theta_levels(i,j,k))
      IF (w_minus_wstar(i,j,k)  >   1.0) THEN
        w_minus_wstar(i,j,k) = 1.0
      ELSE IF (w_minus_wstar(i,j,k)  <   -1.0) THEN
        w_minus_wstar(i,j,k) = -1.0
      END IF
      work(i,j,k) = w_minus_wstar(i,j,k) * (theta(i,j,k+1) -                   &
        theta(i,j,k))
    END DO
  END DO
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i,wrk)            &
!$OMP& SHARED(model_levels,rows,row_length,timestep,w,w_star,                  &
!$OMP& r_theta_levels,theta,work,w_minus_wstar)
  DO k = 2, model_levels - 1
    DO j = 1, rows
      DO i = 1, row_length
        ! CFL=(w-w*) * dt / dz
        ! In this case it is 2 dz so CFL limit is 0.5 not 1
        wrk = timestep * ( w(i,j,k) - w_star(i,j,k) )                          &
          / (r_theta_levels(i,j,k+1) -                                         &
          r_theta_levels(i,j,k-1))
        wrk = MAX(-0.5,MIN(0.5,wrk))
        work(i,j,k) = wrk * (theta(i,j,k+1) - theta(i,j,k-1))
        w_minus_wstar(i,j,k) = wrk
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! ----------------------------------------------------------------------
  ! Section 6.    Perform Horizontal Interpolation
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  ! Section 6.1   Extend data arrays.
  ! ----------------------------------------------------------------------

  IF ( model_domain  ==  mt_global .OR.                                        &
    model_domain  ==  mt_cyclic_lam .OR.                                       &
    model_domain  ==  mt_bi_cyclic_lam) THEN

    DO k = 1, model_levels - 1
      DO j = 1, rows
        DO i = 1, row_length
          ext_data (i,j,k,1) = theta(i,j,k) -                                  &
            (1.-alpha_2) * work(i,j,k)
        END DO
      END DO
    END DO
    k = model_levels
    DO j = 1, rows
      DO i = 1, row_length
        ext_data (i,j,k,1) = theta(i,j,k)
      END DO
    END DO

    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(                                                          &
      ext_data(1-halo_i,1-halo_j,1,1),                                         &
      row_length, rows, model_levels,                                          &
      halo_i, halo_j, fld_type_p, l_vector)

  ELSE IF (model_domain  ==  mt_lam ) THEN

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          ext_data (i,j,k,1) = theta(i,j,k)
        END DO
      END DO
    END DO

    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(                                                          &
      ext_data(1-halo_i,1-halo_j,1,1),                                         &
      row_length, rows, model_levels,                                          &
      halo_i, halo_j, fld_type_p, l_vector)

    ! Set data on edge processor edge haloes from lateral boundary data

    l_do_halos = .TRUE.
    l_do_boundaries = .FALSE.

    ! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(                                               &
      row_length,rows, halo_i, halo_j,                                         &
      model_levels, fld_type_p,                                                &
      ext_data(1-halo_i,1-halo_j,1,1),                                         &
      lenrim, lbc_size, lbc_start,                                             &
      halo_i, halo_j, theta_lbcs,                                              &
      rimwidth, rimwidth, rimweights,                                          &
      at_extremity,                                                            &
      l_do_boundaries, l_do_halos)

    DO k = 1, model_levels - 1

      DO j = 1, rows
        DO i = 1, row_length
          ext_data (i,j,k,2) = work(i,j,k)
        END DO
      END DO
    END DO
    DO j = 1, rows
      DO i = 1, row_length
        ext_data (i,j,model_levels,2) = 0.
      END DO
    END DO

    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(                                                          &
      ext_data(1-halo_i,1-halo_j,1,2),                                         &
      row_length, rows, model_levels,                                          &
      halo_i, halo_j, fld_type_p, l_vector)

    ! Set data on edge processor edge haloes to zero

    l_do_halos = .TRUE.
    l_do_boundaries = .FALSE.

    ! DEPENDS ON: zero_lateral_boundaries
    CALL zero_lateral_boundaries(                                              &
      row_length, rows,                                                        &
      halo_i, halo_j,                                                          &
      model_levels, fld_type_p,                                                &
      ext_data(1-halo_i,1-halo_j,1,2),                                         &
      rimwidth, at_extremity,                                                  &
      l_do_boundaries, l_do_halos)

  END IF  !   model_domain

  ! ----------------------------------------------------------------------
  ! Section 6.2   Call interpolation code.
  ! ----------------------------------------------------------------------

  ! Call high order scheme if required.

  IF (l_high_theta ) THEN

    IF (high_order_scheme_theta  ==  cubiclagrange  .OR.                       &
      high_order_scheme_theta  ==  hcubic_vlin .OR.                            &
      high_order_scheme_theta  ==  hcubic_vquintic) THEN

      ! DEPENDS ON: cubic_lagrange_niv
      CALL cubic_lagrange_niv (                                                &
        ext_data,                                                              &
        row_length, rows, model_levels,                                        &
        row_length, rows, model_levels,                                        &
        halo_i, halo_j,                                                        &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, depart_level,                                            &
        row_length, rows,                                                      &
        lambda_rm, lambda_rp, phi_rm, phi_rp,                                  &
        recip_lambda_m, recip_lambda_0,                                        &
        recip_lambda_p, recip_lambda_p2,                                       &
        recip_phi_m, recip_phi_0,                                              &
        recip_phi_p, recip_phi_p2,                                             &
        l_regular,                                                             &
        model_domain,                                                          &
        at_extremity, n_procx, n_procy,                                        &
        g_row_length, g_rows,                                                  &
        proc_col_group, proc_row_group,                                        &
        datastart,                                                             &
        off_x, off_y,                                                          &
        theta_star)

      IF (model_domain  ==  mt_lam) THEN

        ! DEPENDS ON: cubic_lagrange_niv
        CALL cubic_lagrange_niv (                                              &
          ext_data(1-halo_i,1-halo_j,-1,2),                                    &
          row_length, rows, model_levels,                                      &
          row_length, rows, model_levels,                                      &
          halo_i, halo_j,                                                      &
          weight_lambda, weight_phi,                                           &
          i_out, j_out, depart_level,                                          &
          row_length, rows,                                                    &
          lambda_rm, lambda_rp, phi_rm, phi_rp,                                &
          recip_lambda_m, recip_lambda_0,                                      &
          recip_lambda_p, recip_lambda_p2,                                     &
          recip_phi_m, recip_phi_0,                                            &
          recip_phi_p, recip_phi_p2,                                           &
          l_regular,                                                           &
          model_domain,                                                        &
          at_extremity, n_procx, n_procy,                                      &
          g_row_length, g_rows,                                                &
          proc_col_group, proc_row_group,                                      &
          datastart,                                                           &
          halo_i, halo_j,                                                      &
          work2)

      END IF

    ELSE IF (high_order_scheme_theta  ==  quinticlagrange) THEN

      ! DEPENDS ON: quintic_lagrange_niv
      CALL quintic_lagrange_niv (                                              &
        ext_data,                                                              &
        row_length, rows, model_levels,                                        &
        row_length, rows, model_levels,                                        &
        halo_i, halo_j,                                                        &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, depart_level,                                            &
        off_x, off_y,                                                          &
        theta_star)

      IF (model_domain  ==  mt_lam) THEN

        ! DEPENDS ON: quintic_lagrange_niv
        CALL quintic_lagrange_niv (                                            &
          ext_data(1-halo_i,1-halo_j,-1,2),                                    &
          row_length, rows, model_levels,                                      &
          row_length, rows, model_levels,                                      &
          halo_i, halo_j,                                                      &
          weight_lambda, weight_phi,                                           &
          i_out, j_out, depart_level,                                          &
          halo_i, halo_j,                                                      &
          work2)

      END IF

    ELSE IF (high_order_scheme_theta  ==  ecmwf_quasicubic .OR.                &
      high_order_scheme_theta  ==  hquasicubic_vquintic)                       &
      THEN

      ! DEPENDS ON: ecmwf_quasi_cubic_niv
      CALL ecmwf_quasi_cubic_niv (                                             &
        ext_data,                                                              &
        row_length, rows, model_levels,                                        &
        row_length, rows, model_levels,                                        &
        halo_i, halo_j,                                                        &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, depart_level,                                            &
        off_x, off_y,                                                          &
        theta_star)

      IF (model_domain  ==  mt_lam) THEN

        ! DEPENDS ON: ecmwf_quasi_cubic_niv
        CALL ecmwf_quasi_cubic_niv (                                           &
          ext_data(1-halo_i,1-halo_j,-1,2),                                    &
          row_length, rows, model_levels,                                      &
          row_length, rows, model_levels,                                      &
          halo_i, halo_j,                                                      &
          weight_lambda, weight_phi,                                           &
          i_out, j_out, depart_level,                                          &
          halo_i, halo_j,                                                      &
          work2)

      END IF

    END IF

  END IF  !  L_high_theta

  IF (l_mono_theta .AND. l_high_theta) THEN
    ! Perform montonicity enforcement if high order scheme
    ! used and monotonicity required.

    DO k = 1, model_levels
      DO j = 1,rows
        DO i = 1, row_length

          ! Find max and min monotone values for the point concerned.

          max_mono = MAX (                                                     &
            ext_data(i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k),1),         &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k),1),       &
            ext_data(i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k),1),       &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k),1) )

          min_mono = MIN (                                                     &
            ext_data(i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k),1),         &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k),1),       &
            ext_data(i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k),1),       &
            ext_data(i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k),1) )

          IF (theta_star(i,j,k)  >   max_mono )                                &
            theta_star(i,j,k) = max_mono
          IF (theta_star(i,j,k)  <   min_mono )                                &
            theta_star(i,j,k) = min_mono

        END DO
      END DO
    END DO

  ELSE IF (l_mono_theta ) THEN
    ! Call monotone scheme if required.

    IF (monotone_scheme_theta  ==  trilinear ) THEN

      ! DEPENDS ON: bi_linear_niv
      CALL bi_linear_niv (                                                     &
        ext_data,                                                              &
        row_length, rows, model_levels,                                        &
        row_length, rows, model_levels,                                        &
        halo_i, halo_j,                                                        &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, depart_level,                                            &
        off_x, off_y,                                                          &
        theta_star)

      IF (model_domain  ==  mt_lam) THEN

        ! DEPENDS ON: bi_linear_niv
        CALL bi_linear_niv (                                                   &
          ext_data(1-halo_i,1-halo_j,-1,2),                                    &
          row_length, rows, model_levels,                                      &
          row_length, rows, model_levels,                                      &
          halo_i, halo_j,                                                      &
          weight_lambda, weight_phi,                                           &
          i_out, j_out, depart_level,                                          &
          halo_i, halo_j,                                                      &
          work2)

      END IF

    ELSE IF (monotone_scheme_theta  ==  mono_quasicubic) THEN

      ! DEPENDS ON: ecmwf_mono_quasi_cubic_niv
      CALL ecmwf_mono_quasi_cubic_niv (                                        &
        ext_data,                                                              &
        row_length, rows, model_levels,                                        &
        row_length, rows, model_levels,                                        &
        halo_i, halo_j,                                                        &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, depart_level,                                            &
        off_x, off_y,                                                          &
        theta_star)

      IF (model_domain  ==  mt_lam) THEN

        ! DEPENDS ON: ecmwf_mono_quasi_cubic_niv
        CALL ecmwf_mono_quasi_cubic_niv (                                      &
          ext_data(1-halo_i,1-halo_j,-1,2),                                    &
          row_length, rows, model_levels,                                      &
          row_length, rows, model_levels,                                      &
          halo_i, halo_j,                                                      &
          weight_lambda, weight_phi,                                           &
          i_out, j_out, depart_level,                                          &
          halo_i, halo_j,                                                      &
          work2)

      END IF  !  model_domain  ==  mt_LAM

    END IF  !  monotone_scheme_theta

  END IF    !  L_mono_theta .and. L_high_theta

  ! And now do the extra points

  IF (n_proc > 1 .AND. l_2dcomm  .OR.                                          &
    n_procx > 1 .AND. model_domain == mt_global .AND.                          &
    .NOT.l_2dcomm) THEN

    IF (dim_e_out  >   0) THEN

      IF (l_high_theta ) THEN

        IF (high_order_scheme_theta  ==  cubiclagrange  .OR.                   &
          high_order_scheme_theta  ==  hcubic_vlin .OR.                        &
          high_order_scheme_theta  ==  hcubic_vquintic) THEN

          ! DEPENDS ON: cubic_lagrange_niv
          CALL cubic_lagrange_niv (                                            &
            ext_data,                                                          &
            row_length, rows, model_levels,                                    &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j,                                                    &
            weight_lambda_e, weight_phi_e,                                     &
            i_out_e, j_out_e, depart_level_e,                                  &
            row_length, rows,                                                  &
            lambda_rm, lambda_rp, phi_rm, phi_rp,                              &
            recip_lambda_m, recip_lambda_0,                                    &
            recip_lambda_p, recip_lambda_p2,                                   &
            recip_phi_m, recip_phi_0,                                          &
            recip_phi_p, recip_phi_p2,                                         &
            l_regular,                                                         &
            model_domain,                                                      &
            at_extremity, n_procx, n_procy,                                    &
            g_row_length, g_rows,                                              &
            proc_col_group, proc_row_group,                                    &
            datastart,                                                         &
            0, 0,                                                              &
            theta_star_e)


        ELSE IF (high_order_scheme_theta  ==  quinticlagrange) THEN

          ! DEPENDS ON: quintic_lagrange_niv
          CALL quintic_lagrange_niv (                                          &
            ext_data,                                                          &
            row_length, rows, model_levels,                                    &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j,                                                    &
            weight_lambda_e, weight_phi_e,                                     &
            i_out_e, j_out_e, depart_level_e,                                  &
            0, 0,                                                              &
            theta_star_e)


        ELSE IF (high_order_scheme_theta  ==  ecmwf_quasicubic .OR.            &
          high_order_scheme_theta  ==  hquasicubic_vquintic)                   &
          THEN

          ! DEPENDS ON: ecmwf_quasi_cubic_niv
          CALL ecmwf_quasi_cubic_niv (                                         &
            ext_data,                                                          &
            row_length, rows, model_levels,                                    &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j,                                                    &
            weight_lambda_e, weight_phi_e,                                     &
            i_out_e, j_out_e, depart_level_e,                                  &
            0, 0,                                                              &
            theta_star_e)


        END IF

      END IF

      IF (l_mono_theta .AND. l_high_theta) THEN
        ! Perform montonicity enforcement if high order scheme
        ! used and monotonicity required.

        DO i = 1, dim_e_out

          ! Find max and min monotone values for the point concerned.

          max_mono = MAX (                                                     &
            ext_data(i_out_e(i),j_out_e(i),depart_level_e(i),1),               &
            ext_data(i_out_e(i)+1,j_out_e(i),depart_level_e(i),1),             &
            ext_data(i_out_e(i),j_out_e(i)+1,depart_level_e(i),1),             &
            ext_data(i_out_e(i)+1,j_out_e(i)+1,depart_level_e(i),1) )

          min_mono = MIN (                                                     &
            ext_data(i_out_e(i),j_out_e(i),depart_level_e(i),1),               &
            ext_data(i_out_e(i)+1,j_out_e(i),depart_level_e(i),1),             &
            ext_data(i_out_e(i),j_out_e(i)+1,depart_level_e(i),1),             &
            ext_data(i_out_e(i)+1,j_out_e(i)+1,depart_level_e(i),1) )

          IF (theta_star_e(i)  >   max_mono )                                  &
            theta_star_e(i) = max_mono
          IF (theta_star_e(i)  <   min_mono )                                  &
            theta_star_e(i) = min_mono

        END DO

      ELSE IF (l_mono_theta ) THEN
        ! Call monotone scheme if required.

        IF (monotone_scheme_theta  ==  trilinear) THEN

          ! DEPENDS ON: bi_linear_niv
          CALL bi_linear_niv (                                                 &
            ext_data,                                                          &
            row_length, rows, model_levels,                                    &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j,                                                    &
            weight_lambda_e, weight_phi_e,                                     &
            i_out_e, j_out_e, depart_level_e,                                  &
            0, 0,                                                              &
            theta_star_e)


        ELSE IF (monotone_scheme_theta  ==  mono_quasicubic) THEN

          ! DEPENDS ON: ecmwf_mono_quasi_cubic_niv
          CALL ecmwf_mono_quasi_cubic_niv (                                    &
            ext_data,                                                          &
            row_length, rows, model_levels,                                    &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j,                                                    &
            weight_lambda_e, weight_phi_e,                                     &
            i_out_e, j_out_e, depart_level_e,                                  &
            0, 0,                                                              &
            theta_star_e)


        END IF

      END IF

    END IF ! (dim_e_out  >   0)

    IF ( .NOT. l_2dcomm) THEN

      ! NB: Message passing does not cope with model_domain eq 2 as it fails
      !     to return work2_e to the originating processor.

      nsend = 0
      DO i = 0, n_procx-1
        IF (n_recvfrom(i)  >   0) THEN
          len = n_recvfrom(i)
          DO j = 1, len
            send_data(j,i) = theta_star_e(nsend+j)
            send_data(len+j,i) = r_out_e(nsend+j)
          END DO
          CALL gc_rsend(60*(me+1)+ibase+i, 2*len, ibase+i, info,               &
            recv_data(1,ime), send_data(1,i))
          nsend = nsend + n_recvfrom(i)
        END IF
      END DO

      CALL gcg_ssync(proc_row_group,info)

      DO i = 0, n_procx-1
        IF (n_sendto(i)  >   0) THEN
          len = n_sendto(i)
          CALL gc_rrecv(60*(ibase+i+1)+me, 2*len, ibase+i, info,               &
            recv_data(1,i), send_data(1,ime))
          DO j = 1, len
            theta_star(i_store(j,i), j_store(j,i),                             &
              k_store(j,i)) = recv_data(j,i)
            r_out(i_store(j,i), j_store(j,i),                                  &
              k_store(j,i)) = recv_data(len+j,i)
          END DO
        END IF
      END DO

      ! DEPENDS ON: cnislt_distribute_poles
      CALL cnislt_distribute_poles(                                            &
        sp_send, sp_levels, np_send, np_levels,                                &
        n_procx, ibase, proc_row_group,                                        &
        rows, row_length, model_levels,                                        &
        off_x, off_y, 0, 0,                                                    &
        ime,  at_extremity, .FALSE.,                                           &
        theta_star, r_out)


    ELSE  ! L_2dcomm=TRUE
      ! NB: Message passing does not cope with model_domain eq 2 as it fails
      !     to return work2_e to the originating processor.

      nsend = 0
      DO i = 0, n_proc-1
        IF (n_recvfrom(i)  >   0) THEN
          len = n_recvfrom(i)
          DO j = 1, len
            send_data(j,i) = theta_star_e(nsend+j)
            send_data(len+j,i) = r_out_e(nsend+j)
          END DO
          CALL gc_rsend(60*(me+1)+i, 2*len, i, info,                           &
            recv_data(1,me), send_data(1,i))
          nsend = nsend + n_recvfrom(i)
        END IF
      END DO

      CALL gcg_ssync(proc_all_group,info)

      DO i = 0, n_proc-1
        IF (n_sendto(i)  >   0) THEN
          len = n_sendto(i)
          CALL gc_rrecv(60*(i+1)+me, 2*len, i, info,                           &
            recv_data(1,i), send_data(1,me))
          DO j = 1, len
            theta_star(i_store(j,i), j_store(j,i),                             &
              k_store(j,i)) = recv_data(j,i)
            r_out(i_store(j,i), j_store(j,i),                                  &
              k_store(j,i)) = recv_data(len+j,i)
          END DO
        END IF
      END DO

      IF (model_domain == mt_global) THEN
        IF (at_extremity(psouth)) THEN
          DO k = 1, model_levels
            DO i = 2, row_length
              theta_star(i,1   ,k) = theta_star(1,1   ,k)
              r_out(i,1   ,k) = r_out(1,1   ,k)
            END DO
          END DO
        END IF
        IF (at_extremity(pnorth)) THEN
          DO k = 1, model_levels
            DO i = 2, row_length
              theta_star(i,rows,k) = theta_star (1,rows,k)
              r_out(i,rows,k) = r_out(1,rows,k)
            END DO
          END DO
        END IF
      END IF ! (model_domain  ==  mt_Global)

    END IF  ! .not. L_2dcomm

  END IF  ! nproc>1 and l_2dcomm or nprocx>1 and not.l_2dcomm

  ! ----------------------------------------------------------------------
  ! Section 7.    Add on Eulerian vertical advection term.
  ! ----------------------------------------------------------------------

  DO j=1,rows
    DO i = 1, row_length
      a_coeff_a(i,j) = alpha_2
    END DO
  END DO

  IF (model_domain  ==  mt_global .OR.                                         &
    model_domain  ==  mt_cyclic_lam .OR.                                       &
    model_domain  ==  mt_bi_cyclic_lam) THEN

    IF (model_domain  ==  mt_cyclic_lam) THEN
      IF (at_extremity(psouth)) THEN
        DO i = 1, row_length
          a_coeff_a(i,1) = 1.0
        END DO
      ELSE IF (at_extremity(pnorth)) THEN
        DO i = 1, row_length
          a_coeff_a(i,rows) = 1.0
        END DO
      END IF
    END IF

    IF ( ( .NOT. l_new_tdisc ) .OR. cycleno == 1 ) THEN

      ! calculate first approximation to theta_star

!$OMP  PARALLEL DEFAULT(NONE) SHARED(theta_star, work2, a_coeff_a,             &
!$OMP& w_minus_wstar, work, model_levels, rows, row_length) PRIVATE(i, j, k)
 
!$OMP DO SCHEDULE(STATIC)      
      DO k = 1, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            work2(i,j,k) = theta_star(i,j,k)                                   &
              - a_coeff_a(i,j) * work(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO

      ! ----------------------------------------------------------------
      !  Theta on top surface is advected only in the horizontal.
      ! ----------------------------------------------------------------
      k = model_levels
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          work2(i,j,k) = theta_star(i,j,k)
        END DO
      END DO
!$OMP END DO

      ! calculate second approximation to theta_star
      k = 1
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          theta_star(i,j,k) = work2(i,j,k)                                     &
            - a_coeff_a(i,j) * w_minus_wstar(i,j,k)                            &
            * (work2(i,j,k+1) - work2(i,j,k))                                  &
            + a_coeff_a(i,j) * work(i,j,k)
        END DO
      END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      DO k = 2, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            theta_star(i,j,k) = work2(i,j,k)                                   &
              - a_coeff_a(i,j) * w_minus_wstar(i,j,k)                          &
              * (work2(i,j,k+1) - work2(i,j,k-1))                              &
              + a_coeff_a(i,j) * work(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO 

!$OMP END PARALLEL 
    ELSE ! CycleNo >1 .and. L_new_tdisc=T

      ! calculate second approximation to theta_star
      ! Compute theta* = theta* - a2*dt*[w^n-w*]d2r(theta^(1))
      ! where theta^(1) the theta^n+1 estimate from the 1st sweep.

      k = 1
      DO j = 1, rows
        DO i = 1, row_length
          theta_star(i,j,k) = theta_star(i,j,k)                                &
            - a_coeff_a(i,j) * w_minus_wstar(i,j,k)                            &
            * (theta_np1(i,j,k+1) - theta_np1(i,j,k))
        END DO
      END DO

    
      DO k = 2, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            theta_star(i,j,k) = theta_star(i,j,k)                              &
              - a_coeff_a(i,j) * w_minus_wstar(i,j,k)                          &
              * (theta_np1(i,j,k+1) - theta_np1(i,j,k-1))
          END DO
        END DO
      END DO
    END IF  ! .NOT. L_new_tdisc  .OR.  CycleNo == 1

  ELSE IF (model_domain  ==  mt_lam ) THEN

    IF (l_lbc_old) THEN

      ! set up array for alpha_2 to take into account LAM boundaries
      IF (at_extremity(pnorth)) THEN
        DO j = rows-halo_j+1, rows
          DO i = 1, row_length
            a_coeff_a(i,j) = 1.0
          END DO
        END DO
      END IF
      IF (at_extremity(psouth)) THEN
        DO j = 1, halo_j
          DO i = 1, row_length
            a_coeff_a(i,j) = 1.0
          END DO
        END DO
      END IF
      IF (at_extremity(pwest)) THEN
        DO j = 1, rows
          DO i = 1, halo_i
            a_coeff_a(i,j) = 1.0
          END DO
        END DO
      END IF
      IF (at_extremity(peast)) THEN
        DO j = 1, rows
          DO i = row_length-halo_i+1, row_length
            a_coeff_a(i,j) = 1.0
          END DO
        END DO
      END IF

    END IF !  L_lbc_old

    IF ( ( .NOT. l_new_tdisc ) .OR. cycleno == 1 ) THEN

      ! calculate first approximation to theta_star
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)                &
!$OMP& SHARED(model_levels,rows,row_length,work2,theta_star,a_coeff_a,         &
!$OMP& work)
      DO k = 1, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            work2(i,j,k) = theta_star(i,j,k)                                   &
              - a_coeff_a(i,j) * work(i,j,k)                                   &
              - (1.0-a_coeff_a(i,j)) * work2(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      ! ----------------------------------------------------------------
      !  Theta on top surface is advected only in the horizontal.
      ! ----------------------------------------------------------------
      k = model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work2(i,j,k) = theta_star(i,j,k)
        END DO
      END DO

      ! calculate second approximation to theta_star
      k = 1
      DO j = 1, rows
        DO i = 1, row_length
          theta_star(i,j,k) = work2(i,j,k)                                     &
            - a_coeff_a(i,j) * w_minus_wstar(i,j,k)                            &
            * (work2(i,j,k+1) - work2(i,j,k))                                  &
            + a_coeff_a(i,j) * work(i,j,k)
        END DO
      END DO

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)                &
!$OMP& SHARED(model_levels,rows,row_length,theta_star,work2,                   &
!$OMP& a_coeff_a,w_minus_wstar,work)
      DO k = 2, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            theta_star(i,j,k) = work2(i,j,k)                                   &
              - a_coeff_a(i,j) * w_minus_wstar(i,j,k)                          &
              * (work2(i,j,k+1) - work2(i,j,k-1))                              &
              + a_coeff_a(i,j) * work(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

    ELSE  ! L_new_tdisc .and. CycleNo > 1

      ! calculate work2=theta_dl - (1-a2)*dt*[(w^n-w^*)*d2r(theta)]_dl

      DO k = 1, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            work2(i,j,k) = theta_star(i,j,k)                                   &
              - (1.0-a_coeff_a(i,j)) * work2(i,j,k)
          END DO
        END DO
      END DO

      ! calculate second approximation to theta_star:
      ! theta* = theta_dl - (1-a2)*dt*[(w^n-w^*)*d2r(theta)]_dl
      !                   - a2*dt*(w^(1)_w*)*d2r(theta^(1))

      k = 1
      DO j = 1, rows
        DO i = 1, row_length
          theta_star(i,j,k) = work2(i,j,k)                                     &
            - a_coeff_a(i,j) * w_minus_wstar(i,j,k)                            &
            * (theta_np1(i,j,k+1) - theta_np1(i,j,k))
        END DO
      END DO

      DO k = 2, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            theta_star(i,j,k) = work2(i,j,k)                                   &
              - a_coeff_a(i,j) * w_minus_wstar(i,j,k)                          &
              * (theta_np1(i,j,k+1) - theta_np1(i,j,k-1))
          END DO
        END DO
      END DO

    END IF !  .NOT. L_new_tdisc .OR. CycleNo == 1

  END IF ! on model_domain


  ! end conditional on zero error code
END IF

! End of routine.
IF (lhook) CALL dr_hook('CALC_NON_INT_SL_THETA',zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_non_int_sl_theta

