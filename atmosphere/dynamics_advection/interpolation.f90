! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine Interpolation.


SUBROUTINE interpolation(                                                      &
  data_in1, data_in2, data_in3,                                                &
  eta_in,                                                                      &
  r_in, delta_r_in, r_in_w, itype,                                             &
  number_of_inputs,                                                            &
  check_bottom_levels,                                                         &
  interp_vertical_search_tol,                                                  &
  first_flat_level_in,                                                         &
  dim_i_in, dim_j_in, dim_k_in,                                                &
  dim_j_in_w,                                                                  &
  dim_i_out, dim_j_out, dim_k_out,                                             &
  delta_lambda, delta_phi,                                                     &
  base_lambda, base_phi,                                                       &
  glambda_u, phi_v,                                                            &
  gdlambda_p, dphi_p, gdlambda_u, dphi_v,                                      &
  grecip_dlamu, recip_dphiv,                                                   &
  lambda_rm, lambda_rp,                                                        &
  phi_rm, phi_rp,                                                              &
  recip_lambda_m, recip_lambda_0,                                              &
  recip_lambda_p, recip_lambda_p2,                                             &
  recip_phi_m, recip_phi_0,                                                    &
  recip_phi_p, recip_phi_p2,                                                   &
  i_out_w_in, j_out_w_in,                                                      &
  weight_lambda_w_in, weight_phi_w_in,                                         &
  high_order_scheme, monotone_scheme,                                          &
  cos_latitude, l_regular,                                                     &
  l_vector, model_domain, l_high, l_mono,                                      &
  l_conserv,                                                                   &
  r_out, lambda_out, phi_out,                                                  &
  me, n_proc, n_procx, n_procy,                                                &
  halo_i, halo_j, g_row_length, g_rows,                                        &
  row_length, rows, n_rows, rows_g,                                            &
  datastart, at_extremity, g_i_pe,                                             &
  g_j_pe, l_2dcomm, size_2dcomm,                                               &
  group_2dcomm, size_int, proc_all_group,                                      &
  proc_row_group, proc_col_group,                                              &
  pole_handling,                                                               &
  halo_data_out_i, halo_data_out_j,                                            &
  l_sl_halo_reprod, off_x, off_y,                                              &
  data_out1, data_out2, data_out3,                                             &
  error_code)

! Purpose:
!          Performs interpolation of a field or fields defined on one
!          grid to another grid. The current grids that can be handled
!          can be found in the documentation as can the current options
!          for interpolation schemes. The desired interpolation can be
!          monotone and conservative if desired. Input data can be on a
!          sphere or be a rectangular box. Requested output points must
!          lie inside the region defined for the input Data.
!          The number of fields interpolated is either 1,2 or 3 and is
!          controlled by the switch number_of_inputs. If only one field
!          is to be interpolated the fields Data_in2/3 and Data_out2/3
!          should be set to dummy arguments as they are not used.

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


USE dynamics_grid_mod, ONLY: l_vatpoles

USE mpl, ONLY :                                                                &
  mpl_integer,                                                                 &
  mpl_real,                                                                    &
  mpl_status_size


USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParParams
USE um_types,    ONLY: integer32
USE highos_mod,  ONLY: cubicLagrange, quinticLagrange, ECMWF_quasiCubic,       &
                       ECMWF_mono_quasiCubic, hCubic_vLin,                     &
                       hQuasiCubic_vQuintic, hCubic_vQuintic
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(in) :: dim_i_in  ! Dimension of Data_in in i direction.
INTEGER, INTENT(in) :: dim_j_in  ! Dimension of Data_in in j direction.
INTEGER, INTENT(in) :: dim_j_in_w  ! Dimension of Data_in in j direction.
INTEGER, INTENT(in) :: dim_k_in  ! Dimension of Data_in in k direction.
INTEGER, INTENT(in) :: dim_i_out ! Dimension of Data_out in i direction.
INTEGER, INTENT(in) :: dim_j_out ! Dimension of Data_out in j direction.
INTEGER, INTENT(in) :: dim_k_out ! Dimension of Data_out in j direction.
INTEGER, INTENT(in) :: me        ! My processor number
INTEGER, INTENT(in) :: n_procx   ! Number of processors E-W
INTEGER, INTENT(in) :: n_proc     ! Number of processors in total
INTEGER, INTENT(in) :: n_procy   ! Number of processors N-S
INTEGER, INTENT(in) :: halo_i    ! Size of halo in i direction.
INTEGER, INTENT(in) :: halo_j    ! Size of halo in j direction.
INTEGER, INTENT(in) :: off_x     ! Size of small halo in i direction.
INTEGER, INTENT(in) :: off_y     ! Size of small halo in j direction.
INTEGER, INTENT(in) :: pole_handling ! How to treat the poles:
!   0 - no calculations at the poles
!   1 - poles in one point
!   2 - poles in all points
INTEGER, INTENT(in) :: row_length    ! points in a row on this pe
INTEGER, INTENT(in) :: rows          ! rows on this pe
INTEGER, INTENT(in) :: n_rows        ! rows for v-field dynamic arrays
INTEGER, INTENT(in) :: datastart(3)  ! First gridpoints held by this processor.
INTEGER, INTENT(in) :: g_row_length  ! global number of points on a row
INTEGER, INTENT(in) :: g_rows        ! global number of rows
INTEGER, INTENT(in) :: g_i_pe(1-halo_i:g_row_length+halo_i)  ! processor on my
! processor-row holding a given value in i direction
INTEGER, INTENT(in) :: g_j_pe(1-halo_j:g_rows+halo_j)        ! processor on my
! processor-column holding a given value in j direction
INTEGER, INTENT(in) :: proc_all_group ! Group id for all processors
INTEGER, INTENT(in) :: proc_row_group ! Group id for processors on the same row
INTEGER, INTENT(in) :: proc_col_group ! Group id for processors on the same column
INTEGER, INTENT(in) :: rows_g        ! rows for input field dynamic arrays

INTEGER, INTENT(in) :: halo_data_out_i  ! size of data out halo in i direction
INTEGER, INTENT(in) :: halo_data_out_j  ! size of data out halo in j direction

LOGICAL, INTENT(in) :: l_sl_halo_reprod  ! if true then sl code bit repoducible
! with any sensible halo size
LOGICAL, INTENT(in) :: l_2dcomm
INTEGER, INTENT(in) :: size_2dcomm
INTEGER, INTENT(in) :: size_int
INTEGER, INTENT(in) :: group_2dcomm

INTEGER                                                                        &
  itype                                                                        &
                                ! Defines via an integer code the nature of the
                                ! interpolation in terms of which grid the input
                                ! Data is on. The codes are given in
                                ! terms of a primary variable that would be held
                                ! at that point and are u=1, v=2, w=3.
  , high_order_scheme                                                          &
                                ! a code saying which high order scheme to
                                ! use.
  , monotone_scheme                                                            &
                                ! a code saying which monotone scheme to use.
  , number_of_inputs                                                           &
                                !number of fields to interpolate.
  , interp_vertical_search_tol                                                 &
                                !number of levels either side of
                                ! default level to search.
  , check_bottom_levels                                                        &
                                ! used in interpolation code, and is
                                ! the number of levels to check to see
                                ! if the departure point lies inside the
                                ! orography.
  , first_flat_level_in

LOGICAL                                                                        &
  l_vector                                                                     &
                                ! True, if data is a horizontal vector component,
                                ! False, then data is a scalar.
  , l_high                                                                     &
                                ! True, if high order interpolation required.
  , l_mono                                                                     &
                                ! True, if interpolation required to be monotone.
  , l_conserv                                                                  &
                                ! True, if interpolation to be monotone and
                                !       conservative.
  , l_regular      ! False if variable resolution

INTEGER                                                                        &
  model_domain     ! holds integer code for model domain

REAL                                                                           &
  delta_lambda                                                              &
                                ! holds spacing between points in the i
                                ! direction for the input data field.
  , delta_phi                                                               &
                                ! holds spacing between points in the j
                                ! direction for the input data field.
  , base_lambda                                                                &
  , base_phi

!  VarRes horizontal co-ordinate spacing etc
REAL                                                                           &
  glambda_u(1-halo_i : g_row_length + halo_i)                                  &
  , phi_v( 1-halo_i : row_length + halo_i                                      &
  ,        1-halo_j : n_rows + halo_j )                                        &
  , gdlambda_p( 1-halo_i : g_row_length + halo_i)                              &
  , gdlambda_u( 1-halo_i : g_row_length + halo_i)                              &
  , dphi_p    ( 1-halo_i : row_length + halo_i                                 &
  ,             1-halo_j : rows + halo_j )                                     &
  , dphi_v    ( 1-halo_i : row_length + halo_i                                 &
  ,             1-halo_j : n_rows + halo_j )                                   &
  , grecip_dlamu(1-halo_i : g_row_length + halo_i)                             &
  , recip_dphiv( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : n_rows + halo_j )                                  &
  , lambda_rm (1-halo_i : row_length + halo_i)                                 &
  , lambda_rp (1-halo_i : row_length + halo_i)                                 &
  , phi_rm    ( 1-halo_i : row_length + halo_i                                 &
  ,             1-halo_j : rows_g + halo_j )                                   &
  , phi_rp    ( 1-halo_i : row_length + halo_i                                 &
  ,             1-halo_j : rows_g + halo_j )                                   &
  , recip_lambda_m(1-halo_i : row_length + halo_i)                             &
  , recip_lambda_0(1-halo_i : row_length + halo_i)                             &
  , recip_lambda_p(1-halo_i : row_length + halo_i)                             &
  , recip_lambda_p2(1-halo_i : row_length + halo_i)                            &
  , recip_phi_m( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : rows_g + halo_j )                                  &
  , recip_phi_0( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : rows_g + halo_j )                                  &
  , recip_phi_p( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : rows_g + halo_j )                                  &
  , recip_phi_p2( 1-halo_i : row_length + halo_i                               &
  ,               1-halo_j : rows_g + halo_j )

REAL                                                                           &
  eta_in(dim_k_in) ! eta coordinate levels.

REAL                                                                           &
  data_in1 (1-halo_i:dim_i_in+halo_i,                                          &
  1-halo_j:dim_j_in+halo_j, dim_k_in)                                          &
                                ! data to be
                                ! interpolated
  , data_in2 (1-halo_i:dim_i_in+halo_i,                                        &
  1-halo_j:dim_j_in+halo_j, dim_k_in)                                          &
                                ! optional second
                                ! field of data to
                                ! be interpolated
  , data_in3 (1-halo_i:dim_i_in+halo_i,                                        &
  1-halo_j:dim_j_in+halo_j, dim_k_in)                                          &
                                ! optional third
                                ! field of data to
                                ! be interpolated
  , r_in (1-halo_i:dim_i_in+halo_i,                                            &
  1-halo_j:dim_j_in+halo_j, dim_k_in)                                          &
                                ! Vertical
                                ! co-ordinate
                                ! of input data.
  , r_in_w(1-halo_i:dim_i_in+halo_i,                                           &
  1-halo_j:dim_j_in_w+halo_j, dim_k_in)                                        &
                                ! Vertical
                                ! co-ordinate
                                ! of input data
                                ! on w grid
  , delta_r_in (dim_i_in, dim_j_in, dim_k_in)                                  &
                                ! Vertical
                                ! layer thickness
                                ! of input data.
  , cos_latitude (1-off_x:dim_i_in+off_x,                                      &
  1-off_y:dim_j_in+off_y)
REAL                                                                           &
  lambda_out (dim_i_out, dim_j_out, dim_k_out)                                 &
                                ! Lambda
                                ! co-ordinate of
                                ! output data on
                                ! input.
  , phi_out (dim_i_out, dim_j_out, dim_k_out)                                  &
                                ! Phi Co-ordinate
                                ! of output data
                                ! on input.
  , r_out (dim_i_out, dim_j_out, dim_k_out)       ! Vertical
! co-ordinate
! of output data.

LOGICAL                                                                        &
  at_extremity(4)  ! Indicates if this processor is at north,
! south, east or west of the processor grid

REAL                                                                           &
  weight_lambda_w_in (dim_i_out, dim_j_out, dim_k_out)                         &
  , weight_phi_w_in (dim_i_out, dim_j_out, dim_k_out)
REAL                                                                           &
  weight_lambda_w (dim_i_out, dim_j_out, dim_k_out)                            &
  , weight_phi_w (dim_i_out, dim_j_out, dim_k_out)

INTEGER                                                                        &
  i_out_w_in (dim_i_out, dim_j_out, dim_k_out)                                 &
  , j_out_w_in (dim_i_out, dim_j_out, dim_k_out)

! Arguments with Intent OUT. ie: Output variables.

REAL                                                                           &
                                ! data interpolated to desired locations.
  data_out1(1-halo_data_out_i:dim_i_out+halo_data_out_i,                       &
  1-halo_data_out_j:dim_j_out+halo_data_out_j,                                 &
  dim_k_out)                                                                   &
  , data_out2(1-halo_data_out_i:dim_i_out+halo_data_out_i,                     &
  1-halo_data_out_j:dim_j_out+halo_data_out_j,                                 &
  dim_k_out)                                                                   &
  , data_out3(1-halo_data_out_i:dim_i_out+halo_data_out_i,                     &
  1-halo_data_out_j:dim_j_out+halo_data_out_j,                                 &
  dim_k_out)

INTEGER                                                                        &
  error_code     ! Non-zero on exit if error detected.

! Local Variables.

! scalars

LOGICAL                                                                        &
  conserv_fail

INTEGER                                                                        &
  i, j, k, n                                                                   &
                                ! Loop indices
  , dsm1x                                                                      &
  , dsm1y                                                                      &
  , index                                                                      &
  , lambda_start
! pointer for start of lambda_p/lambda_u on this pe

REAL                                                                           &
  recip_delta_lambda                                                        &
  , recip_delta_phi                                                         &
  , temp

! arrays

INTEGER (KIND=integer32) ::                                                    &
    i_out (dim_i_out, dim_j_out, dim_k_out)                                    &
  , j_out (dim_i_out, dim_j_out, dim_k_out)                                    &
  , k_out (dim_i_out, dim_j_out, dim_k_out)                                    &
  , i_out_w (dim_i_out, dim_j_out, dim_k_out)                                  &
  , j_out_w (dim_i_out, dim_j_out, dim_k_out)


REAL                                                                           &
  ext_data(1-halo_i:dim_i_in+halo_i+1,1-halo_j:dim_j_in+halo_j,                &
  -1:dim_k_in+2,number_of_inputs)                                              &
  , data_out_high (dim_i_out, dim_j_out, dim_k_out,                            &
  number_of_inputs)                                                            &
  , data_out_mono (dim_i_out, dim_j_out, dim_k_out,                            &
  number_of_inputs)                                                            &
  , weight_lambda (dim_i_out, dim_j_out, dim_k_out)                            &
  , weight_phi (dim_i_out, dim_j_out, dim_k_out)                               &
  , coeff_z(dim_i_out, dim_j_out, dim_k_out, -2:3)                             &
  , coeff_z_lin(dim_i_out, dim_j_out, dim_k_out, 0:1)

! Varibles applied in the "compute-on-demand" strategy

INTEGER                                                                        &
  ime, ibase, irecv, my_imin, my_imax, dim_e_out, h_factor                     &
  , nsend, nrecv, info, len, itmp, j0, j1, kk, sender                          &
  , my_iminp, my_imaxp, my_jmin, my_jmax


INTEGER (KIND=integer32) :: i_store   (size_int,0:size_2dcomm)
INTEGER (KIND=integer32) :: j_store   (size_int,0:size_2dcomm)
INTEGER (KIND=integer32) :: k_store   (size_int,0:size_2dcomm)
INTEGER (KIND=integer32) :: i_out_e   (size_int)
INTEGER (KIND=integer32) :: j_out_e   (size_int)
INTEGER (KIND=integer32) :: k_out_e   (size_int)
INTEGER (KIND=integer32) :: i_out_w_e (size_int)
INTEGER (KIND=integer32) :: j_out_w_e (size_int)

INTEGER :: n_sendto  (0:size_2dcomm)
INTEGER :: n_recvfrom(0:size_2dcomm)
INTEGER :: k_level_e (size_int)
INTEGER :: isend_arr (5,size_int,0:size_2dcomm)
INTEGER :: irecv_arr (5,size_int,0:size_2dcomm)
INTEGER :: sp_send   (0:size_2dcomm)
INTEGER :: sp_levels (0:size_2dcomm,dim_k_out)
INTEGER :: np_send   (0:size_2dcomm)
INTEGER :: np_levels (0:size_2dcomm,dim_k_out)
REAL :: rsend_arr        (5,size_int,0:size_2dcomm)
REAL :: rrecv_arr        (5,size_int,0:size_2dcomm)
REAL :: weight_lambda_e  (size_int)
REAL :: weight_phi_e     (size_int)
REAL :: weight_lambda_w_e(size_int)
REAL :: weight_phi_w_e   (size_int)
REAL :: data_out_high_e  (size_int*number_of_inputs)
REAL :: data_out_mono_e  (size_int*number_of_inputs)
REAL :: send_data        (2*number_of_inputs*size_int,0:size_2dcomm)
REAL :: recv_data        (2*number_of_inputs*size_int,0:size_2dcomm)
REAL :: coeff_z_e        (size_int, -2:3)
REAL :: coeff_z_lin_e    (size_int, 0:1)
REAL :: r_out_e          (size_int)

INTEGER :: errorstatus

! Variables used for index vectors

INTEGER                                                                        &
  ic,icnt,inx(dim_i_out)

INTEGER d_imin,d_imax, d_jmin, d_jmax
INTEGER                                                                        &
  lon_dem(3), gon_dem(3)
COMMON / ondem / lon_dem, gon_dem

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
  INTEGER  :: i_out_w
  INTEGER  :: j_out_w
  INTEGER  :: k
  REAL     :: weight_lambda
  REAL     :: weight_phi
  REAL     :: r_out
  REAL     :: weight_lambda_w
  REAL     :: weight_phi_w
END TYPE sendrecv_type

TYPE (sendrecv_type) :: send_array(size_int,0:size_2dcomm)
TYPE (sendrecv_type) :: recv_array(size_int,0:size_2dcomm)

! Statement function for optimising floor
! Can only be used when the result is cast into an integer
REAL :: real_temp_var
REAL :: my_floor
my_floor(real_temp_var) = real_temp_var - 0.5 + SIGN(0.5,real_temp_var)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Description: COMDECK containing the allowed
!              monotone scheme options
!
      INTEGER                                                           &
     &     triLinear                                                    &
     &,    mono_quasiCubic

      PARAMETER(                                                        &
     &     triLinear       = 1                                          &
     &,    mono_quasiCubic = 2 )

! Functions: None

! ----------------------------------------------------------------------
!  Section 1.   Error trap Un-supported options.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('INTERPOLATION',zhook_in,zhook_handle)

dsm1x = datastart(1)-1
dsm1y = datastart(2)-1

error_code = 0

! Execute rest of routine only if error code is still zero.

IF (error_code  ==  0 ) THEN

  ! ----------------------------------------------------------------------
  ! Section 2.   Call appropriate routine to extend data.
  !              Minimum amount of extending done to cope with highest
  !              order interpolation scheme requested. This can leave
  !              unset values in Ext_Data and Ext_r_in, but
  !              these values will not be used.

  !              Parallel version: No extension in r required.
  ! ----------------------------------------------------------------------

  IF ( (high_order_scheme  ==  quinticlagrange .OR.                            &
    high_order_scheme  ==  hquasicubic_vquintic .OR.                           &
    high_order_scheme  ==  hcubic_vquintic) .AND. l_high )                     &
    THEN

    ! DEPENDS ON: extend_data_quintic
    CALL extend_data_quintic (data_in1, data_in2, data_in3,                    &
      dim_i_in, dim_j_in, dim_k_in,                                            &
      number_of_inputs, halo_i, halo_j,                                        &
      ext_data)

  ELSE IF ( ( ( high_order_scheme  ==  cubiclagrange                           &
    .OR. high_order_scheme  ==  ecmwf_quasicubic                               &
    .OR. high_order_scheme  ==  ecmwf_mono_quasicubic                          &
    .OR. high_order_scheme  ==  hcubic_vlin )                                  &
    .AND. l_high )                                                             &
    .OR. ( monotone_scheme  ==  mono_quasicubic                                &
    .AND. l_mono ) ) THEN

    ! DEPENDS ON: extend_data_cubic
    CALL extend_data_cubic  (data_in1, data_in2, data_in3,                     &
      dim_i_in, dim_j_in, dim_k_in,                                            &
      number_of_inputs, halo_i, halo_j,                                        &
      ext_data)

  ELSE

    ! DEPENDS ON: extend_data_linear
    CALL extend_data_linear (data_in1, data_in2, data_in3,                     &
      dim_i_in, dim_j_in, dim_k_in,                                            &
      number_of_inputs, halo_i, halo_j,                                        &
      ext_data)

  END IF

  ! ----------------------------------------------------------------------
  ! Section 3.   For each output point find i,j,k so that the point on the
  !              output grid lies between i and i+1, j and j+1, k and k+1
  ! ----------------------------------------------------------------------

  IF ( l_regular ) THEN
    recip_delta_lambda = 1./delta_lambda
    recip_delta_phi = 1./delta_phi
  END IF ! L_regular

  IF (.NOT. l_2dcomm) THEN

! i_out_w and j_out_w must be copied because they change in this subroutine

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(i_out_w, j_out_w, &
!$OMP& j_out_w_in, i_out_w_in, weight_lambda_w, weight_phi_w,              &
!$OMP& weight_lambda_w_in, weight_phi_w_in, dim_k_out, dim_j_out,          &
!$OMP& dim_i_out) PRIVATE(k, j, i)
    DO k = 1, dim_k_out
      DO j = 1, dim_j_out
        DO i = 1, dim_i_out
          i_out_w(i,j,k) = i_out_w_in(i,j,k)
          j_out_w(i,j,k) = j_out_w_in(i,j,k)
          weight_lambda_w(i,j,k) = weight_lambda_w_in(i,j,k)
          weight_phi_w(i,j,k)    = weight_phi_w_in(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 

    ! Find i and j point.

    IF ( l_regular ) THEN
      IF ( itype == 1) THEN            ! u points.

       IF (l_vatpoles) THEN
!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,temp)           &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,lambda_out,i_out,                  &
!$OMP&     recip_delta_lambda,weight_lambda,j_out,j_out_w,              &
!$OMP&     weight_phi, weight_phi_w)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              temp = lambda_out(i,j,k) * recip_delta_lambda
              i_out(i,j,k) = my_floor(1. + temp)
              weight_lambda(i,j,k) = temp - (i_out(i,j,k) - 1.0)
            END DO
          END DO
          !         End Do        ! loop removed to improve speedup

!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out      ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              j_out(i,j,k)   = j_out_w(i,j,k)
              weight_phi(i,j,k) = weight_phi_w(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
       ELSE
!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,temp)           &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,lambda_out,i_out,                  &
!$OMP&     recip_delta_lambda,weight_lambda,j_out,j_out_w,              &
!$OMP&     weight_phi, weight_phi_w)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              temp = lambda_out(i,j,k) * recip_delta_lambda
              i_out(i,j,k) = my_floor(.5 + temp)
              weight_lambda(i,j,k) = temp - ( i_out(i,j,k) - 0.5 )
            END DO
          END DO
          !         End Do        ! loop removed to improve speedup

!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out      ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              j_out(i,j,k)   = j_out_w(i,j,k)
              weight_phi(i,j,k) = weight_phi_w(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
       END IF ! vatpoles

      ELSE IF ( itype == 2) THEN       ! v points.

       IF (l_vatpoles) THEN
!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,temp)    &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out, i_out_w, i_out,            &
!$OMP&    j_out, weight_lambda_w,phi_out, base_phi,                     &
!$OMP&    weight_lambda, weight_phi, recip_delta_phi,dsm1y)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              i_out(i,j,k)   = i_out_w(i,j,k)
              weight_lambda(i,j,k) = weight_lambda_w(i,j,k)
            END DO
          END DO
          !         End Do   ! loop removed to improve speedup
!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out     ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              temp = (phi_out(i,j,k) - base_phi) * recip_delta_phi
              j_out(i,j,k) = my_floor(1.0 + temp)
              weight_phi(i,j,k) = temp - (j_out(i,j,k) - 1.0)
              j_out(i,j,k) = j_out(i,j,k) - dsm1y
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
       ELSE
!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,temp)    &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out, i_out_w, i_out,            &
!$OMP&    j_out, weight_lambda_w,phi_out, base_phi,                     &
!$OMP&    weight_lambda, weight_phi, recip_delta_phi,dsm1y)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              i_out(i,j,k)   = i_out_w(i,j,k)
              weight_lambda(i,j,k) = weight_lambda_w(i,j,k)
            END DO
          END DO
          !         End Do   ! loop removed to improve speedup
!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out     ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              temp = (phi_out(i,j,k) - base_phi) * recip_delta_phi
              j_out(i,j,k) = my_floor(.5 + temp)
              weight_phi(i,j,k) = temp - ( j_out(i,j,k) - 0.5)
              j_out(i,j,k) = j_out(i,j,k) - dsm1y
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
       END IF ! vatpoles

      ELSE                            ! w points.

!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)         &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,i_out,i_out_w,              &
!$OMP&     weight_lambda_w, weight_lambda,j_out,j_out_w,                &
!$OMP&     weight_phi,weight_phi_w)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              i_out(i,j,k)   = i_out_w(i,j,k)
              weight_lambda(i,j,k) = weight_lambda_w(i,j,k)
            END DO
          END DO
          !         End Do     ! loop removed to improve speedup
!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out    ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              j_out(i,j,k)   = j_out_w(i,j,k)
              weight_phi(i,j,k) = weight_phi_w(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO

      END IF  ! itype == 1

    ELSE  ! variable resolution

      IF ( itype == 1) THEN    ! u points

        DO k = 1, dim_k_out
        DO j = 1, dim_j_out
        DO i = 1, dim_i_out
          IF ( lambda_out(i,j,k) + base_lambda <                               &
            glambda_u(i_out_w(i,j,k)) ) THEN
            i_out(i,j,k) = i_out_w(i,j,k) - 1
          ELSE   ! depart point is between w(i) and u(i)
            i_out(i,j,k) = i_out_w(i,j,k)
          END IF ! depart point is before u(i)
          weight_lambda(i,j,k) = ( lambda_out(i,j,k) + base_lambda -           &
            glambda_u(i_out(i,j,k)) ) *                                        &
            grecip_dlamu(i_out(i,j,k))
          j_out(i,j,k) = j_out_w(i,j,k)
          weight_phi(i,j,k) = weight_phi_w(i,j,k)
        END DO  ! i loop
        END DO  ! j loop
        END DO  ! k loop

      ELSE IF ( itype == 2 ) THEN    ! v points

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)                &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,i_out_w,i_out,                     &
!$OMP&     weight_phi,j_out_w,     weight_lambda_w,phi_out,phi_v,              &
!$OMP&     weight_lambda,j_out,recip_dphiv)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              IF ( phi_out(i,j,k) <                                            &
                phi_v(i, j_out_w(i,j,k)) ) THEN
                j_out(i,j,k) = j_out_w(i,j,k) - 1
              ELSE         ! depart point is between w(j) and v(j)
                j_out(i,j,k) = j_out_w(i,j,k)
              END IF        ! depart point is before v(j)
              weight_phi(i,j,k) = ( phi_out(i,j,k) -                           &
                phi_v(i, j_out(i,j,k)) ) *                                     &
                recip_dphiv(i, j_out(i,j,k))
              i_out(i,j,k) = i_out_w(i,j,k)
              weight_lambda(i,j,k) =  weight_lambda_w(i,j,k)
            END DO
          END DO
        END DO  !  k = 1, dim_k_out
!$OMP END PARALLEL DO

      ELSE   !  itype == 3   w - points
        DO k = 1, dim_k_out
        DO j = 1, dim_j_out
        DO i = 1, dim_i_out
          i_out(i,j,k) = i_out_w(i,j,k)
          j_out(i,j,k) = j_out_w(i,j,k)
          weight_lambda(i,j,k) =  weight_lambda_w(i,j,k)
          weight_phi(i,j,k) = weight_phi_w(i,j,k)
        END DO  ! i loop
        END DO  ! j loop
        END DO  ! k loop
      END IF   !  itype == 1

    END IF ! L_regular


    IF ( n_procx > 1 ) THEN

      IF ( model_domain == mt_global ) THEN

        ! Send the points outside my region to the appropriate processor for
        ! interpolation. Only performed if the domain is decomposed in the
        ! i direction.

        ! The first and last point I can interpolate in, based on available
        ! data on this processor (minus/plus one to avoid use of ge/le)

        h_factor = 2
        IF (high_order_scheme  ==  quinticlagrange .AND. l_high) THEN
          h_factor = 3
        END IF
        my_imin = datastart(1) - halo_i + h_factor - 1
        my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor             &
          + 1

        ! values for use in polar row to ensure pole is only calculated on one
        ! processor, (minus/plus one to avoid use of ge/le)
        my_iminp = datastart(1) - 1
        my_imaxp = datastart(1)+dim_i_out
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
                i_out_w(i,1,k) = i
              END DO
            END DO

          END IF

          IF (at_extremity(pnorth)) THEN
            j1 = dim_j_out-1

            DO k = 1, dim_k_out
              DO i = 1, dim_i_out
                i_out(i,dim_j_out,k) = i
                i_out_w(i,dim_j_out,k) = i
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
              IF (i_out(1,1,k)  >   my_iminp .AND.                             &
                i_out(1,1,k)  <   my_imaxp) THEN
                i_out(1,1,k) = i_out(1,1,k) - datastart(1) + 1
                i_out_w(1,1,k) = i_out_w(1,1,k) - datastart(1) + 1
                sp_send(ime) = sp_send(ime) + 1
                sp_levels(ime,sp_send(ime)) = k
                DO i = 2, dim_i_out
                  i_out(i,1,k) = i ! i_out(i,1,k) - datastart(1) + 1
                  i_out_w(i,1,k) = i
                END DO
              ELSE
                sender = g_i_pe(i_out(1,1,k))
                sp_send(sender) = sp_send(sender) + 1
                sp_levels(sender,sp_send(sender)) = k
                DO i = 1, dim_i_out
                  i_out(i,1,k) = i
                  i_out_w(i,1,k) = i
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
              IF (i_out(1,dim_j_out,k)  >   my_iminp .AND.                     &
                i_out(1,dim_j_out,k)  <   my_imaxp)  THEN

                np_send(ime) = np_send(ime) + 1
                np_levels(ime,np_send(ime)) = k
                i_out(1,dim_j_out,k) =                                         &
                  i_out(1,dim_j_out,k) - datastart(1) + 1
                i_out_w(1,dim_j_out,k) =                                       &
                  i_out_w(1,dim_j_out,k) - datastart(1) + 1
                DO i = 2, dim_i_out
                  i_out(i,dim_j_out,k) = i
                  i_out_w(i,dim_j_out,k) = i
                END DO
              ELSE
                sender = g_i_pe(i_out(1,dim_j_out,k))
                np_send(sender) = np_send(sender) + 1
                np_levels(sender,np_send(sender)) = k
                DO i = 1, dim_i_out
                  i_out(i,dim_j_out,k) = i
                  i_out_w(i,dim_j_out,k) = i
                END DO
              END IF
            END DO
          END IF

        END IF ! on pole_handling

        IF ( l_sl_halo_reprod) THEN

          ! On the global boundaries, use i_out < 1 or i_out > g_row_length
          ! if that makes local computation possible. Not required when
          ! L_sl_halo_reprod is false is other logic ensures this is done.

          ! This code unsafe if applied at poles, where it isn't required.

          IF (at_extremity(pwest)) THEN
            DO k = 1, dim_k_out
              DO j = j0, j1
                DO i = 1, halo_i
                  IF (i_out(i,j,k)  >   g_row_length-halo_i+h_factor)          &
                    THEN
                    i_out(i,j,k) = i_out(i,j,k) - g_row_length
                    i_out_w(i,j,k) = i_out_w(i,j,k) - g_row_length
                  END IF
                END DO
              END DO
            END DO
          END IF
          IF (at_extremity(peast)) THEN
            DO k = 1, dim_k_out
              DO j = j0, j1
                DO i = dim_i_out-halo_i+1, dim_i_out
                  IF (i_out(i,j,k)  <   halo_i-h_factor+1)                     &
                    THEN
                    i_out(i,j,k) = i_out(i,j,k) + g_row_length
                    i_out_w(i,j,k) = i_out_w(i,j,k) + g_row_length
                  END IF
                END DO
              END DO
            END DO
          END IF

        END IF ! on L_sl_halo_reprod

        ! And now decide where a point should be evaluated

        d_imin=my_imin-datastart(1) + 1
        d_imax=my_imax-datastart(1) + 1

        !               If (i_out(i,j,k)  >   my_imin .and.
        !    &              i_out(i,j,k)  <   my_imax) then
       
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(k, j, i, inx, icnt,&
!$OMP& ic, irecv, itmp)
        DO k = 1, dim_k_out
          DO j = j0, j1
            icnt=0
            !OMP index x locations and how many of them exist (+reduction) 
            DO i = 1, dim_i_out 
              i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
              i_out_w(i,j,k) = i_out_w(i,j,k) - datastart(1) + 1
              IF (i_out(i,j,k)   <=  d_imin .OR.   &
                i_out(i,j,k)   >=  d_imax) THEN
                
                icnt = icnt+1
                inx(icnt) = i
              END IF
            END DO

            !OMP use index locations to determine 
            ! Process locally, so find the local destination
            DO ic=1, icnt
              i = inx(ic)
              i_out(i,j,k) = i_out(i,j,k) + datastart(1) - 1
              i_out_w(i,j,k) = i_out_w(i,j,k) + datastart(1) - 1
              ! Send to a remote processor, given by the array g_i_pe
              !     CODE TO STOP BIT NON-REPRODUCIBILITY
              IF (i_out(i,j,k) >  g_row_length+halo_i-h_factor) THEN
                i_out(i,j,k)=i_out(i,j,k)-g_row_length
                i_out_w(i,j,k)=i_out_w(i,j,k)-g_row_length
              END IF
              IF (i_out(i,j,k) <  1-halo_i+h_factor) THEN
                i_out(i,j,k)=i_out(i,j,k)+g_row_length
                i_out_w(i,j,k)=i_out_w(i,j,k)+g_row_length
              END IF
              !     END CODE TO STOP BIT NON-REPRODUCIBILITY
              irecv = g_i_pe(i_out(i,j,k)) ! find proc to send to (no lock)

! each thread updates the sendto array indicating the irecv proc has 
! an additional point to process  

!$OMP CRITICAL 
! note number of point 'irecv' proc received 
              n_sendto(irecv) = n_sendto(irecv) + 1 
              itmp = n_sendto(irecv) ! get point 
!$OMP END CRITICAL 

              send_array(itmp, irecv) % i_out           = i_out(i,j,k) 
              send_array(itmp, irecv) % j_out           = j_out(i,j,k)
              send_array(itmp, irecv) % i_out_w         = i_out_w(i,j,k)
              send_array(itmp, irecv) % j_out_w         = j_out_w(i,j,k)
              send_array(itmp, irecv) % k               = k
              send_array(itmp, irecv) % weight_lambda   =                      &
                weight_lambda(i,j,k)
              send_array(itmp, irecv) % weight_phi      = weight_phi(i,j,k)
              send_array(itmp, irecv) % r_out           = r_out(i,j,k)
              send_array(itmp, irecv) % weight_lambda_w =                      &
                weight_lambda_w(i,j,k)
              send_array(itmp, irecv) % weight_phi_w    =                      &
                weight_phi_w(i,j,k)
              i_store(itmp,irecv) = i
              j_store(itmp,irecv) = j
              k_store(itmp,irecv) = k
              i_out(i,j,k) = i
              i_out_w(i,j,k) = i
              !                End If
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO 

        ! Send the points to be evaluated outside my region

        ! Counts can be distributed via an alltoall with the row communicator
        CALL mpl_alltoall(n_sendto,       1,    mpl_integer,                   &
          n_recvfrom,     1,    mpl_integer,                                   &
          proc_row_group, info)

        ! Get types setup if not done
        IF (mpl_send_type == imdi) THEN
          offsets    (0) = 0
          oldtypes   (0) = mpl_integer
          blockcounts(0) = 5

          CALL mpl_type_extent(mpl_integer, extent, info)

          offsets    (1) = 5 * extent
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


        dim_e_out = 0
        DO i = 0, n_procx-1
          IF (n_recvfrom(i)  >   0) THEN
            DO j = 1, n_recvfrom(i)
              dim_e_out = dim_e_out + 1
              i_out_e(dim_e_out) = recv_array(j,i) % i_out - datastart(1)+1
              j_out_e(dim_e_out) = recv_array(j,i) % j_out
              i_out_w_e(dim_e_out) = recv_array(j,i) % i_out_w                 &
                - datastart(1)+1
              j_out_w_e(dim_e_out) = recv_array(j,i) % j_out_w
              k_level_e(dim_e_out) = recv_array(j,i) % k
              weight_lambda_e(dim_e_out) = recv_array(j,i) % weight_lambda
              weight_phi_e(dim_e_out) = recv_array(j,i) % weight_phi
              r_out_e(dim_e_out) = recv_array(j,i) % r_out
              weight_lambda_w_e(dim_e_out) = recv_array(j,i) % weight_lambda_w
              weight_phi_w_e(dim_e_out) = recv_array(j,i) % weight_phi_w
            END DO
          END IF
        END DO


        IF (dim_e_out  >   size_int)               THEN
          errorstatus = 10
          
          CALL ereport("Interpolation", ErrorStatus,                           &
            "over-writing due to dim_e_out size" )
        END IF

        IF (dim_e_out  >   0) THEN
          ! calculate level just below departure point in vertical and
          ! vertical interpolation weights.
          ! DEPENDS ON: eta_vert_weights_e
          CALL eta_vert_weights_e(eta_in, r_in_w,                              &
            check_bottom_levels,                                               &
            interp_vertical_search_tol,                                        &
            first_flat_level_in,                                               &
            dim_i_in, dim_j_in_w, dim_k_in,                                    &
            dim_k_out, dim_e_out,                                              &
            high_order_scheme, monotone_scheme,                                &
            model_domain, l_high, l_mono,                                      &
            l_conserv,                                                         &
            halo_i, halo_j,                                                    &
            r_out_e, k_level_e,                                                &
            coeff_z_e, coeff_z_lin_e,                                          &
            k_out_e, i_out_w_e, j_out_w_e,                                     &
            weight_lambda_w_e, weight_phi_w_e )

        END IF

      ELSE  !  model_domain is NOT GLOBAL
        DO k = 1, dim_k_out
        DO j = 1, dim_j_out
        DO i = 1, dim_i_out
          i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
          i_out_w(i,j,k) = i_out_w(i,j,k) - datastart(1) + 1
        END DO !  i loop
        END DO !  j loop
        END DO !  k loop

      END IF ! model_domain == mt_Global
    END IF ! n_procx > 1

    


    !     CODE TO STOP BIT NON-REPRODUCIBILITY
    IF (model_domain == mt_global .AND. n_procx == 1) THEN
      h_factor = 2
      IF (high_order_scheme  ==  quinticlagrange .AND. l_high) THEN
        h_factor = 3
      END IF
      my_imin = datastart(1) - halo_i + h_factor - 1
      my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor + 1
      DO k = 1, dim_k_out
        DO j = 1,dim_j_out
          DO i = 1, dim_i_out
            IF (i_out(i,j,k) >= my_imax) THEN
              i_out(i,j,k)=i_out(i,j,k) - g_row_length
              i_out_w(i,j,k)=i_out_w(i,j,k) - g_row_length
            END IF
            IF (i_out(i,j,k) <= my_imin ) THEN
              i_out(i,j,k)=i_out(i,j,k)+g_row_length
              i_out_w(i,j,k)=i_out_w(i,j,k)+g_row_length
            END IF
          END DO
        END DO
      END DO
    END IF   !model_domain == 1 .and. n_procx == 1
    !     END CODE TO STOP BIT NON-REPRODUCIBILITY


  ELSE ! l_2dcomm

    ! i_out_w and j_out_w must be copied because they change in this subroutine
    ! i_out is global , j_out is local
    DO k = 1, dim_k_out
      DO j = 1, dim_j_out
        DO i = 1, dim_i_out
          i_out_w(i,j,k) = i_out_w_in(i,j,k)
          j_out_w(i,j,k) = j_out_w_in(i,j,k) + dsm1y
          weight_lambda_w(i,j,k) = weight_lambda_w_in(i,j,k)
          weight_phi_w(i,j,k)    = weight_phi_w_in(i,j,k)
        END DO
      END DO
    END DO
    ! i_out is now global , j_out is now global

    ! Find i and j point.

    

    IF ( l_regular ) THEN
      IF ( itype == 1) THEN            ! u points.

       IF (l_vatpoles) THEN
!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,temp)           &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,lambda_out,i_out,                  &
!$OMP&     recip_delta_lambda,weight_lambda,j_out,j_out_w,                  &
!$OMP&     weight_phi,weight_phi_w)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              temp = lambda_out(i,j,k) * recip_delta_lambda
              i_out(i,j,k) = my_floor(1.0 + temp)
              weight_lambda(i,j,k) = temp - (i_out(i,j,k) - 1.0)
            END DO
          END DO
          !         End Do        ! loop removed to improve speedup

!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out      ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              j_out(i,j,k) = j_out_w(i,j,k)
              weight_phi(i,j,k) = weight_phi_w(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
       ELSE
!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,temp)           &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,lambda_out,i_out,                  &
!$OMP&     recip_delta_lambda,weight_lambda,j_out,j_out_w,                  &
!$OMP&     weight_phi,weight_phi_w)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              temp = lambda_out(i,j,k) * recip_delta_lambda
              i_out(i,j,k) = my_floor(.5 + temp)
              weight_lambda(i,j,k) = temp - ( i_out(i,j,k) - 0.5 )
            END DO
          END DO
          !         End Do        ! loop removed to improve speedup

!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out      ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              j_out(i,j,k) = j_out_w(i,j,k)
              weight_phi(i,j,k) = weight_phi_w(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
       END IF ! vatpoles

      ELSE IF ( itype == 2) THEN       ! v points.

       IF (l_vatpoles) THEN
!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,temp)           &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,i_out_w,i_out,                     &
!$OMP&     weight_phi,             weight_lambda_w,phi_out,base_phi,           &
!$OMP&     weight_lambda,j_out,  recip_delta_phi,dsm1y)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              i_out(i,j,k) = i_out_w(i,j,k)
              weight_lambda(i,j,k) = weight_lambda_w(i,j,k)
            END DO
          END DO
          !         End Do   ! loop removed to improve speedup
!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out     ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              temp = (phi_out(i,j,k) - base_phi) * recip_delta_phi
              j_out(i,j,k) = my_floor(1.0 + temp)
              weight_phi(i,j,k) = temp - (j_out(i,j,k) - 1.0)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
       ELSE
!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,temp)           &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,i_out_w,i_out,                     &
!$OMP&     weight_phi,             weight_lambda_w,phi_out,base_phi,           &
!$OMP&     weight_lambda,j_out,  recip_delta_phi,dsm1y)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              i_out(i,j,k) = i_out_w(i,j,k)
              weight_lambda(i,j,k) = weight_lambda_w(i,j,k)
            END DO
          END DO
          !         End Do   ! loop removed to improve speedup
!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out     ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              temp = (phi_out(i,j,k) - base_phi) * recip_delta_phi
              j_out(i,j,k) = my_floor(.5 + temp)
              weight_phi(i,j,k) = temp - ( j_out(i,j,k) - 0.5)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
       END IF ! vatpoles

      ELSE                            ! w points.

!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)                &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,i_out,i_out_w,                     &
!$OMP&     weight_lambda_w      ,weight_lambda,j_out,j_out_w,                  &
!$OMP&     weight_phi,weight_phi_w)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              i_out(i,j,k) = i_out_w(i,j,k)
              weight_lambda(i,j,k) = weight_lambda_w(i,j,k)
            END DO
          END DO
          !         End Do     ! loop removed to improve speedup
!cdir nounroll
!cdir collapse
          !         Do k = 1, dim_k_out    ! loop removed to improve speedup
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              j_out(i,j,k) = j_out_w(i,j,k)
              weight_phi(i,j,k) = weight_phi_w(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO

      END IF  ! itype == 1

    ELSE  ! variable resolution

      IF ( itype == 1) THEN    ! u points

        DO k = 1, dim_k_out
        DO j = 1, dim_j_out
        DO i = 1, dim_i_out





          IF ( lambda_out(i,j,k) + base_lambda <                               &
            glambda_u(i_out_w(i,j,k)) ) THEN
            i_out(i,j,k) = i_out_w(i,j,k) - 1
          ELSE   ! depart point is between w(i) and u(i)
            i_out(i,j,k) = i_out_w(i,j,k)
          END IF ! depart point is before u(i)
          weight_lambda(i,j,k) = ( lambda_out(i,j,k) + base_lambda -           &
            glambda_u(i_out(i,j,k)) ) *                                        &
            grecip_dlamu(i_out(i,j,k))
          j_out(i,j,k) = j_out_w(i,j,k)
          weight_phi(i,j,k) = weight_phi_w(i,j,k)
        END DO  !  i loop

        END DO  !  j loop
        END DO  !  k loop


      ELSE IF ( itype == 2 ) THEN    ! v points

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)                &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,i_out_w,i_out,                     &
!$OMP&     weight_phi,j_out_w,     weight_lambda_w,phi_out,phi_v,              &
!$OMP&     weight_lambda,j_out,recip_dphiv)
        DO k = 1, dim_k_out
          DO j = 1, dim_j_out
            DO i = 1, dim_i_out
              IF ( phi_out(i,j,k) <                                            &
                phi_v(i, j_out_w(i,j,k)) ) THEN
                j_out(i,j,k) = j_out_w(i,j,k) - 1
              ELSE         ! depart point is between w(j) and v(j)
                j_out(i,j,k) = j_out_w(i,j,k)
              END IF        ! depart point is before v(j)
              weight_phi(i,j,k) = ( phi_out(i,j,k) -                           &
                phi_v(i, j_out(i,j,k)) ) *                                     &
                recip_dphiv(i, j_out(i,j,k))
              i_out(i,j,k) = i_out_w(i,j,k)
              weight_lambda(i,j,k) =  weight_lambda_w(i,j,k)
            END DO
          END DO
        END DO  !  k = 1, dim_k_out
!$OMP END PARALLEL DO

      ELSE   !  itype == 3   w - points

        DO k = 1, dim_k_out
        DO j = 1, dim_j_out
        DO i = 1, dim_i_out





          i_out(i,j,k) = i_out_w(i,j,k)
          j_out(i,j,k) = j_out_w(i,j,k)
          weight_lambda(i,j,k) =  weight_lambda_w(i,j,k)
          weight_phi(i,j,k) = weight_phi_w(i,j,k)
        END DO  ! i loop

        END DO  ! j loop
        END DO  ! k loop

      END IF   !  itype == 1

    END IF ! L_regular

    

    ! For LAMs we need to check what we want to do if the point is at an edge
    ! and outside the acceptable range.

    IF ( n_procy>1 .OR. model_domain==mt_lam) THEN
      ! The first and last point I can interpolate in, based on available
      ! data on this processor (minus/plus one to avoid use of ge/le)
      h_factor = 2
      IF (high_order_scheme  ==  quinticlagrange .AND. l_high) h_factor = 3
      my_jmin = datastart(2) - halo_j + h_factor - 1
      my_jmax = datastart(2) + dim_j_out - 1 + halo_j - h_factor + 1
    END IF

    IF ( n_procx>1 .OR. model_domain==mt_lam) THEN
      ! The first and last point I can interpolate in, based on available
      ! data on this processor (minus/plus one to avoid use of ge/le)
      h_factor = 2
      IF (high_order_scheme  ==  quinticlagrange .AND. l_high) h_factor = 3
      my_imin = datastart(1) - halo_i + h_factor - 1
      my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor + 1
    END IF

    IF (n_proc > 1) THEN
      IF ( model_domain == mt_global .AND. pole_handling==1) THEN
        ! values for use in polar row to ensure pole is only calculated on one
        ! processor.

        ! values for use in polar row to ensure pole is only calculated on one
        ! processor, (minus/plus one to avoid use of ge/le)
        my_iminp = datastart(1) - 1
        my_imaxp = datastart(1)+dim_i_out
        IF (at_extremity(pwest) ) my_iminp = my_imin
        IF (at_extremity(peast) ) my_imaxp = my_imax
      END IF
    END IF ! n_proc > 1

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
    ! The base processor on this row, and my address relative to that
    ! processor
    IF (n_proc > 1) THEN

      ibase = (me/n_procx) * n_procx
      ime = me - ibase

      DO i = 0, n_proc-1
        n_sendto(i) = 0
      END DO

      j0 = 1
      j1 = dim_j_out

      ! If the pole values not are going to be used, we we can save a
      ! large amount of communication by evaluating local dummys at the poles
      IF ( model_domain == mt_global ) THEN
        IF (pole_handling  ==  0) THEN
          IF (at_extremity(psouth)) THEN
            DO k = 1, dim_k_out
              DO i = 1, dim_i_out
                i_out(i,1,k) = datastart(1)+i-1
                i_out_w(i,1,k) = datastart(1)+i-1
                !ajmpole          j_out(i,1,k) = datastart(2)
                !ajmpole          j_out_w(i,1,k) = datastart(2)
              END DO
            END DO
          END IF
          IF (at_extremity(pnorth)) THEN
            DO k = 1, dim_k_out
              DO i = 1, dim_i_out
                i_out(i,dim_j_out,k) = datastart(1)+i-1
                i_out_w(i,dim_j_out,k) = datastart(1)+i-1
                !ajmpole          j_out(i,dim_j_out,k) = datastart(2)+dim_j_out-1
                !ajmpole          j_out_w(i,dim_j_out,k) = datastart(2)+dim_j_out-1
              END DO
            END DO
          END IF
        END IF !pole_handling=0

        ! If all values along a polar row is evaluated in one point, we can
        ! save a large amount of communication by evaluating only one point
        ! correctly and then broadcast this value.

        IF (pole_handling  ==  1) THEN
          IF (at_extremity(psouth)) THEN
            DO k = 1, dim_k_out
              DO i = 2, dim_i_out
                i_out(i,1,k) = i + datastart(1) - 1
                i_out_w(i,1,k) = i + datastart(1) - 1
                j_out(i,1,k) =  datastart(2)
                j_out_w(i,1,k) =  datastart(2)
              END DO
            END DO
          END IF

          IF (at_extremity(pnorth)) THEN
            DO k = 1, dim_k_out
              DO i = 2, dim_i_out
                i_out(i,dim_j_out,k) = i+datastart(1) - 1
                i_out_w(i,dim_j_out,k) = i+datastart(1) - 1
                j_out(i,dim_j_out,k) = datastart(2)+dim_j_out-1
                j_out_w(i,dim_j_out,k) = datastart(2)+dim_j_out-1
              END DO
            END DO
          END IF

        END IF ! pole_handling=1
      END IF  !model_domain == mt_global
    END IF  ! n_proc > 1

    IF (n_proc > 1) THEN
      IF (model_domain==mt_global) THEN
        IF ( l_sl_halo_reprod) THEN

          ! On the global boundaries, use i_out < 1 or i_out > g_row_length
          ! if that makes local computation possible. Not required when
          ! L_sl_halo_reprod is false is other logic ensures this is done.

          ! This code unsafe if applied at poles, where it isn't required.

          IF (at_extremity(pwest)) THEN
            DO k = 1, dim_k_out
              DO j = j0, j1
                DO i = 1, halo_i
                  IF (i_out(i,j,k)  >   g_row_length-halo_i+h_factor)          &
                    THEN
                    i_out(i,j,k) = i_out(i,j,k) - g_row_length
                    i_out_w(i,j,k) = i_out_w(i,j,k) - g_row_length
                  END IF
                END DO
              END DO
            END DO
          END IF
          IF (at_extremity(peast)) THEN
            DO k = 1, dim_k_out
              DO j = j0, j1
                DO i = dim_i_out-halo_i+1, dim_i_out
                  IF (i_out(i,j,k)  <   halo_i-h_factor+1)                     &
                    THEN
                    i_out(i,j,k) = i_out(i,j,k) + g_row_length
                    i_out_w(i,j,k) = i_out_w(i,j,k) + g_row_length
                  END IF
                END DO
              END DO
            END DO
          END IF

        END IF ! on L_sl_halo_reprod
      END IF ! model_domain==mt_global
    END IF  ! n_proc > 1


    ! And now decide where a point should be evaluated
    IF (n_procx>1 .AND. n_procy==1) THEN
      d_imin=my_imin-dsm1x
      d_imax=my_imax-dsm1x

      !               If (i_out(i,j,k)  >   my_imin .and.
      !    &              i_out(i,j,k)  <   my_imax) then

      DO k = 1, dim_k_out
        DO j = j0, j1
          DO i = 1, dim_i_out
            i_out(i,j,k) = i_out(i,j,k) - dsm1x
            i_out_w(i,j,k) = i_out_w(i,j,k) - dsm1x
            j_out(i,j,k) = j_out(i,j,k) - dsm1y
            j_out_w(i,j,k) = j_out_w(i,j,k) - dsm1y
          END DO
        END DO
      END DO

      DO k = 1, dim_k_out
        DO j = j0, j1
          icnt=0
          DO i = 1, dim_i_out
            IF (i_out(i,j,k)  <=  d_imin .OR.                                  &
              i_out(i,j,k)  >=  d_imax) THEN
              icnt = icnt+1
              inx(icnt) = i
            END IF
          END DO

          ! Process locally, so find the local destination
          DO ic=1, icnt
            i = inx(ic)
            i_out(i,j,k) = i_out(i,j,k) + dsm1x
            i_out_w(i,j,k) = i_out_w(i,j,k) + dsm1x
            j_out(i,j,k) = j_out(i,j,k) + dsm1y
            j_out_w(i,j,k) = j_out_w(i,j,k) + dsm1y
            ! Send to a remote processor, given by the array g_i_pe
            !     CODE TO STOP BIT NON-REPRODUCIBILITY
            IF (i_out(i,j,k) >  g_row_length+halo_i-h_factor) THEN
              i_out(i,j,k)=i_out(i,j,k)-g_row_length
              i_out_w(i,j,k)=i_out_w(i,j,k)-g_row_length
            END IF
            IF (i_out(i,j,k) <  1-halo_i+h_factor) THEN
              i_out(i,j,k)=i_out(i,j,k)+g_row_length
              i_out_w(i,j,k)=i_out_w(i,j,k)+g_row_length
            END IF
            !     END CODE TO STOP BIT NON-REPRODUCIBILITY
            irecv = g_i_pe(i_out(i,j,k))
            n_sendto(irecv) = n_sendto(irecv) + 1
            itmp = n_sendto(irecv)
            send_array(itmp, irecv) % i_out           = i_out(i,j,k)
            send_array(itmp, irecv) % j_out           = j_out(i,j,k)
            send_array(itmp, irecv) % i_out_w         = i_out_w(i,j,k)
            send_array(itmp, irecv) % j_out_w         = j_out_w(i,j,k)
            send_array(itmp, irecv) % k               = k
            send_array(itmp, irecv) % weight_lambda   =                        &
              weight_lambda(i,j,k)
            send_array(itmp, irecv) % weight_phi      = weight_phi(i,j,k)
            send_array(itmp, irecv) % r_out           = r_out(i,j,k)
            send_array(itmp, irecv) % weight_lambda_w =                        &
              weight_lambda_w(i,j,k)
            send_array(itmp, irecv) % weight_phi_w    =                        &
              weight_phi_w(i,j,k)
            i_store(itmp,irecv) = i
            j_store(itmp,irecv) = j
            k_store(itmp,irecv) = k
            i_out(i,j,k) = i
            i_out_w(i,j,k) = i
            j_out(i,j,k) = j
            j_out_w(i,j,k) = j
            !                End If
          END DO
        END DO
      END DO
    END IF ! n_procx>1 .and. n_procy==1

!!!!  2 if Blocks below here to be done
    IF (n_procy>1 .AND. n_procx==1) THEN
      d_jmin=my_jmin-datastart(2) + 1
      d_jmax=my_jmax-datastart(2) + 1

      !               If (i_out(i,j,k)  >   my_imin .and.
      !    &              i_out(i,j,k)  <   my_imax) then

      DO k = 1, dim_k_out
        DO j = j0, j1
          DO i = 1, dim_i_out
            i_out(i,j,k) = i_out(i,j,k) - dsm1x
            i_out_w(i,j,k) = i_out_w(i,j,k) - dsm1x
            j_out(i,j,k) = j_out(i,j,k) - dsm1y
            j_out_w(i,j,k) = j_out_w(i,j,k) - dsm1y
          END DO
        END DO
      END DO

      DO k = 1, dim_k_out
        DO j = j0, j1
          icnt=0
          DO i = 1, dim_i_out
            IF (j_out(i,j,k)  <=  d_jmin .OR.                                  &
              j_out(i,j,k)  >=  d_jmax) THEN
              icnt = icnt+1
              inx(icnt) = i
            END IF
          END DO

          ! Process locally, so find the local destination
          DO ic=1, icnt
            i = inx(ic)
            i_out(i,j,k) = i_out(i,j,k) + dsm1x
            i_out_w(i,j,k) = i_out_w(i,j,k) + dsm1x
            j_out(i,j,k) = j_out(i,j,k) + dsm1y
            j_out_w(i,j,k) = j_out_w(i,j,k) + dsm1y
            ! Send to a remote processor, given by the array g_i_pe
            !     CODE TO STOP BIT NON-REPRODUCIBILITY
            IF (j_out(i,j,k) >  g_rows+halo_j-h_factor) THEN
              j_out(i,j,k)=j_out(i,j,k)-g_rows
              j_out_w(i,j,k)=j_out_w(i,j,k)-g_rows
            END IF
            IF (j_out(i,j,k) <  1-halo_j+h_factor) THEN
              j_out(i,j,k)=j_out(i,j,k)+g_rows
              j_out_w(i,j,k)=j_out_w(i,j,k)+g_rows
            END IF
            !     END CODE TO STOP BIT NON-REPRODUCIBILITY
            irecv = g_j_pe(j_out(i,j,k))
            n_sendto(irecv) = n_sendto(irecv) + 1
            itmp = n_sendto(irecv)
            send_array(itmp, irecv) % i_out           = i_out(i,j,k)
            send_array(itmp, irecv) % j_out           = j_out(i,j,k)
            send_array(itmp, irecv) % i_out_w         = i_out_w(i,j,k)
            send_array(itmp, irecv) % j_out_w         = j_out_w(i,j,k)
            send_array(itmp, irecv) % k               = k
            send_array(itmp, irecv) % weight_lambda   =                        &
              weight_lambda(i,j,k)
            send_array(itmp, irecv) % weight_phi      = weight_phi(i,j,k)
            send_array(itmp, irecv) % r_out           = r_out(i,j,k)
            send_array(itmp, irecv) % weight_lambda_w =                        &
              weight_lambda_w(i,j,k)
            send_array(itmp, irecv) % weight_phi_w    =                        &
              weight_phi_w(i,j,k)
            i_store(itmp,irecv) = i
            j_store(itmp,irecv) = j
            k_store(itmp,irecv) = k
            i_out(i,j,k) = i
            i_out_w(i,j,k) = i
            j_out(i,j,k) = j
            j_out_w(i,j,k) = j
            !                End If
          END DO
        END DO
      END DO
    END IF ! n_procx=1 .and. n_procy> 1


    IF (n_procx>1 .AND. n_procy>1) THEN
      d_imin=my_imin-dsm1x
      d_imax=my_imax-dsm1x
      d_jmin=my_jmin-dsm1y
      d_jmax=my_jmax-dsm1y

      !               If (i_out(i,j,k)  >   my_imin .and.
      !    &              i_out(i,j,k)  <   my_imax) then

      DO k = 1, dim_k_out
        DO j = j0, j1
          DO i = 1, dim_i_out
            i_out(i,j,k) = i_out(i,j,k) - dsm1x
            i_out_w(i,j,k) = i_out_w(i,j,k) - dsm1x
            j_out(i,j,k) = j_out(i,j,k) - dsm1y
            j_out_w(i,j,k) = j_out_w(i,j,k) - dsm1y
          END DO
        END DO
      END DO

      ! count points outside this process for a single 
      ! 1-d line of atmos
      DO k = 1, dim_k_out
        DO j = j0, j1
          icnt=0
          DO i = 1, dim_i_out
            IF (i_out(i,j,k)  <=  d_imin .OR.                                  &
              i_out(i,j,k)  >=  d_imax .OR.                                    &
              j_out(i,j,k)  <=  d_jmin .OR.                                    &
              j_out(i,j,k)  >=  d_jmax ) THEN
              icnt = icnt+1
              inx(icnt) = i
            END IF
          END DO

          ! Process locally, so find the local destination
          DO ic=1, icnt
            i = inx(ic)
            
            ! transfrom back to global coords
            i_out(i,j,k) = i_out(i,j,k) + dsm1x
            i_out_w(i,j,k) = i_out_w(i,j,k) + dsm1x
            j_out(i,j,k) = j_out(i,j,k) + dsm1y
            j_out_w(i,j,k) = j_out_w(i,j,k) + dsm1y

            ! Send to a remote processor, given by the array g_i_pe
            !     CODE TO STOP BIT NON-REPRODUCIBILITY
            IF (i_out(i,j,k) >  g_row_length+halo_i-h_factor) THEN
              i_out(i,j,k)=i_out(i,j,k)-g_row_length
              i_out_w(i,j,k)=i_out_w(i,j,k)-g_row_length
            END IF
            IF (i_out(i,j,k) <  1-halo_i+h_factor) THEN
              i_out(i,j,k)=i_out(i,j,k)+g_row_length
              i_out_w(i,j,k)=i_out_w(i,j,k)+g_row_length
            END IF
            !     END CODE TO STOP BIT NON-REPRODUCIBILITY

            irecv = g_i_pe(i_out(i,j,k)) +                                     &
              n_procx* g_j_pe(j_out(i,j,k)) ! calc process to recv

            ! store list of location in send array
            n_sendto(irecv) = n_sendto(irecv) + 1 
            ! store position of process to recv in send array
            itmp = n_sendto(irecv) 

            ! data to be sent to process
            send_array(itmp, irecv) % i_out           = i_out(i,j,k) 
            send_array(itmp, irecv) % j_out           = j_out(i,j,k)
            send_array(itmp, irecv) % i_out_w         = i_out_w(i,j,k)
            send_array(itmp, irecv) % j_out_w         = j_out_w(i,j,k)
            send_array(itmp, irecv) % k               = k
            send_array(itmp, irecv) % weight_lambda   =                        &
              weight_lambda(i,j,k)
            send_array(itmp, irecv) % weight_phi      = weight_phi(i,j,k)
            send_array(itmp, irecv) % r_out           = r_out(i,j,k)
            send_array(itmp, irecv) % weight_lambda_w =                        &
              weight_lambda_w(i,j,k)
            send_array(itmp, irecv) % weight_phi_w    =                        &
              weight_phi_w(i,j,k)
            i_store(itmp,irecv) = i
            j_store(itmp,irecv) = j
            k_store(itmp,irecv) = k
            
            ! setting advected positions back to normal 
            i_out(i,j,k) = i 
            i_out_w(i,j,k) = i
            j_out(i,j,k) = j
            j_out_w(i,j,k) = j
            !                End If
          END DO
        END DO
      END DO
    END IF ! n_procx>1 .and. n_procy==1
!!!! end of temporary code sections to be revised

    IF (n_proc >1) THEN
      ! Send the points to be evaluated outside my region

      ! Counts can be distributed via an alltoall with the row communicator
      CALL mpl_alltoall(n_sendto,       1,    mpl_integer,                     &
        n_recvfrom,     1,    mpl_integer,                                     &
        proc_all_group, info)

      ! Get types setup if not done
      IF (mpl_send_type == imdi) THEN
        offsets    (0) = 0
        oldtypes   (0) = mpl_integer
        blockcounts(0) = 5

        CALL mpl_type_extent(mpl_integer, extent, info)

        offsets    (1) = 5 * extent
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

      dim_e_out = 0
      DO i = 0, n_proc-1
        IF (n_recvfrom(i)  >   0) THEN
          DO j = 1, n_recvfrom(i)
            dim_e_out = dim_e_out + 1
            i_out_e(dim_e_out) = recv_array(j,i) % i_out - dsm1x
            j_out_e(dim_e_out) = recv_array(j,i) % j_out -dsm1y
            i_out_w_e(dim_e_out) = recv_array(j,i) % i_out_w - dsm1x
            j_out_w_e(dim_e_out) = recv_array(j,i) % j_out_w - dsm1y
            k_level_e(dim_e_out) = recv_array(j,i) % k
            weight_lambda_e(dim_e_out) = recv_array(j,i) % weight_lambda
            weight_phi_e(dim_e_out) = recv_array(j,i) % weight_phi
            r_out_e(dim_e_out) = recv_array(j,i) % r_out
            weight_lambda_w_e(dim_e_out) = recv_array(j,i) % weight_lambda_w
            weight_phi_w_e(dim_e_out) = recv_array(j,i) % weight_phi_w
          END DO
        END IF
      END DO

      IF (dim_e_out  >   size_int) THEN
        errorstatus = 10
        
        CALL ereport("Interpolation", ErrorStatus,                             &
          "over-writing due to dim_e_out size" )
      END IF

      IF (dim_e_out  >   0) THEN
        ! calculate level just below departure point in vertical and
        ! vertical interpolation weights.
        ! DEPENDS ON: eta_vert_weights_e
        CALL eta_vert_weights_e(eta_in, r_in_w,                                &
          check_bottom_levels,                                                 &
          interp_vertical_search_tol,                                          &
          first_flat_level_in,                                                 &
          dim_i_in, dim_j_in_w, dim_k_in,                                      &
          dim_k_out, dim_e_out,                                                &
          high_order_scheme, monotone_scheme,                                  &
          model_domain, l_high, l_mono,                                        &
          l_conserv,                                                           &
          halo_i, halo_j,                                                      &
          r_out_e, k_level_e,                                                  &
          coeff_z_e, coeff_z_lin_e,                                            &
          k_out_e, i_out_w_e, j_out_w_e,                                       &
          weight_lambda_w_e, weight_phi_w_e )

      END IF

    END IF ! n_proc > 1

  END IF ! l_2dcomm

!!! DO WE NEED THIS SECTION ANYMORE
  !     CODE TO STOP BIT NON-REPRODUCIBILITY
  IF (model_domain == mt_global .AND. n_proc == 1 .AND. l_2dcomm .OR.          &
    model_domain == mt_global .AND. n_procx == 1 .AND. .NOT.l_2dcomm           &
    ) THEN
    h_factor = 2
    IF (high_order_scheme  ==  quinticlagrange .AND. l_high) THEN
      h_factor = 3
    END IF
    my_imin = datastart(1) - halo_i + h_factor - 1
    my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor + 1
    DO k = 1, dim_k_out
      DO j = 1,dim_j_out
        DO i = 1, dim_i_out
          IF (i_out(i,j,k) >= my_imax) THEN
            i_out(i,j,k)=i_out(i,j,k) - g_row_length
            i_out_w(i,j,k)=i_out_w(i,j,k) - g_row_length
          END IF
          IF (i_out(i,j,k) <= my_imin ) THEN
            i_out(i,j,k)=i_out(i,j,k)+g_row_length
            i_out_w(i,j,k)=i_out_w(i,j,k)+g_row_length
          END IF
        END DO
      END DO
    END DO
  END IF   !model_domain == 1 .and. n_procx == 1

  
  !     END CODE TO STOP BIT NON-REPRODUCIBILITY



  ! calculate level just below departure point in vertical and
  ! vertical interpolation weights.
  ! DEPENDS ON: eta_vert_weights
  CALL eta_vert_weights    (eta_in, r_in_w,                                    &
    check_bottom_levels,                                                       &
    interp_vertical_search_tol,                                                &
    first_flat_level_in,                                                       &
    dim_i_in, dim_j_in_w, dim_k_in,                                            &
    dim_i_out, dim_j_out, dim_k_out,                                           &
    high_order_scheme, monotone_scheme,                                        &
    model_domain, l_high, l_mono,                                              &
    l_conserv,                                                                 &
    halo_i, halo_j,                                                            &
    r_out,                                                                     &
    coeff_z, coeff_z_lin,                                                      &
    k_out, i_out_w, j_out_w,                                                   &
    weight_lambda_w, weight_phi_w )


  ! ----------------------------------------------------------------------
  ! Section 4.   Perform required Interpolations.
  ! ----------------------------------------------------------------------


  ! Call high order scheme if required.

  IF (l_high ) THEN

    IF (high_order_scheme  ==  cubiclagrange) THEN

      ! DEPENDS ON: cubic_lagrange
      CALL cubic_lagrange (                                                    &
        ext_data,                                                              &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        halo_i, halo_j, number_of_inputs,                                      &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, k_out,                                                   &
        row_length, rows_g,                                                    &
        lambda_rm, lambda_rp, phi_rm, phi_rp,                                  &
        recip_lambda_m, recip_lambda_0,                                        &
        recip_lambda_p, recip_lambda_p2,                                       &
        recip_phi_m, recip_phi_0,                                              &
        recip_phi_p, recip_phi_p2,                                             &
        coeff_z, l_regular,                                                    &
        model_domain,                                                          &
        at_extremity, n_procx, n_procy,                                        &
        g_row_length, g_rows,                                                  &
        proc_col_group, proc_row_group,                                        &
        datastart,                                                             &
        data_out_high)

    ELSE IF (high_order_scheme  ==  quinticlagrange) THEN

      ! DEPENDS ON: quintic_lagrange
      CALL quintic_lagrange (                                                  &
        ext_data,                                                              &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        halo_i, halo_j, number_of_inputs,                                      &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, k_out,                                                   &
        coeff_z,                                                               &
        data_out_high )

    ELSE IF (high_order_scheme  ==  ecmwf_quasicubic) THEN

      ! DEPENDS ON: ecmwf_quasi_cubic
      CALL ecmwf_quasi_cubic (                                                 &
        ext_data,                                                              &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        halo_i, halo_j, number_of_inputs,                                      &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, k_out,                                                   &
        coeff_z,                                                               &
        data_out_high )

    ELSE IF (high_order_scheme  ==  ecmwf_mono_quasicubic) THEN

      ! DEPENDS ON: ecmwf_mono_quasi_cubic
      CALL ecmwf_mono_quasi_cubic (                                            &
        ext_data,                                                              &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        halo_i, halo_j, number_of_inputs,                                      &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, k_out,                                                   &
        coeff_z,                                                               &
        data_out_high)

    ELSE IF (high_order_scheme  ==  hcubic_vlin) THEN

      ! DEPENDS ON: h_cubic_v_linear
      CALL h_cubic_v_linear (                                                  &
        ext_data,                                                              &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        halo_i, halo_j, number_of_inputs,                                      &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, k_out,                                                   &
        coeff_z_lin,                                                           &
        data_out_high)

    ELSE IF (high_order_scheme  ==  hquasicubic_vquintic) THEN

      ! DEPENDS ON: h_quasi_cubic_v_quintic
      CALL h_quasi_cubic_v_quintic (                                           &
        ext_data,                                                              &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        halo_i, halo_j, number_of_inputs,                                      &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, k_out,                                                   &
        coeff_z,                                                               &
        data_out_high)

    ELSE IF (high_order_scheme  ==  hcubic_vquintic) THEN

      ! DEPENDS ON: h_cubic_v_quintic
      CALL h_cubic_v_quintic (                                                 &
        ext_data,                                                              &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        halo_i, halo_j, number_of_inputs,                                      &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, k_out,                                                   &
        coeff_z,                                                               &
        data_out_high)



    END IF

  END IF

  IF (((l_high .AND. l_mono) .OR. l_conserv)                                   &
    .AND. high_order_scheme  /=  ecmwf_mono_quasicubic) THEN

    ! DEPENDS ON: mono_enforce
    CALL mono_enforce(                                                         &
      ext_data, number_of_inputs,                                              &
      dim_i_in, dim_j_in, dim_k_in,                                            &
      dim_i_out, dim_j_out, dim_k_out,                                         &
      halo_i, halo_j,                                                          &
      i_out, j_out, k_out,                                                     &
      data_out_high)

  END IF

  ! Call monotone scheme if required.

  IF ( (l_mono .AND. .NOT. l_high ) .OR. l_conserv ) THEN

    IF (monotone_scheme  ==  trilinear) THEN

      ! DEPENDS ON: tri_linear
      CALL tri_linear (                                                        &
        ext_data,                                                              &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        halo_i, halo_j, number_of_inputs,                                      &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, k_out,                                                   &
        coeff_z_lin,                                                           &
        data_out_mono )

    ELSE IF (monotone_scheme  ==  mono_quasicubic) THEN

      ! DEPENDS ON: ecmwf_mono_quasi_cubic
      CALL ecmwf_mono_quasi_cubic (                                            &
        ext_data,                                                              &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        halo_i, halo_j, number_of_inputs,                                      &
        weight_lambda, weight_phi,                                             &
        i_out, j_out, k_out,                                                   &
        coeff_z,                                                               &
        data_out_mono)

    END IF

  END IF
  ! Repeat the interpolation procedure for the "compute-on-demand" points

  IF (n_proc  >   1  .AND. l_2dcomm   .OR.                                     &
    n_procx > 1 .AND. model_domain == mt_global .AND. .NOT.l_2dcomm)           &
    THEN

    ! Call high order scheme if required.

    IF (dim_e_out  >   0) THEN

      IF (l_high ) THEN

        IF (high_order_scheme  ==  cubiclagrange) THEN

          ! DEPENDS ON: cubic_lagrange
          CALL cubic_lagrange (                                                &
            ext_data,                                                          &
            dim_i_in, dim_j_in, dim_k_in,                                      &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j, number_of_inputs,                                  &
            weight_lambda_e,weight_phi_e,                                      &
            i_out_e, j_out_e, k_out_e,                                         &
            row_length, rows_g,                                                &
            lambda_rm, lambda_rp, phi_rm, phi_rp,                              &
            recip_lambda_m, recip_lambda_0,                                    &
            recip_lambda_p, recip_lambda_p2,                                   &
            recip_phi_m, recip_phi_0,                                          &
            recip_phi_p, recip_phi_p2,                                         &
            coeff_z_e, l_regular,                                              &
            model_domain,                                                      &
            at_extremity, n_procx, n_procy,                                    &
            g_row_length, g_rows,                                              &
            proc_col_group, proc_row_group,                                    &
            datastart,                                                         &
            data_out_high_e )

        ELSE IF (high_order_scheme  ==  quinticlagrange) THEN

          ! DEPENDS ON: quintic_lagrange
          CALL quintic_lagrange (                                              &
            ext_data,                                                          &
            dim_i_in, dim_j_in, dim_k_in,                                      &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j, number_of_inputs,                                  &
            weight_lambda_e,weight_phi_e,                                      &
            i_out_e, j_out_e, k_out_e,                                         &
            coeff_z_e,                                                         &
            data_out_high_e)

        ELSE IF (high_order_scheme  ==  ecmwf_quasicubic) THEN

          ! DEPENDS ON: ecmwf_quasi_cubic
          CALL ecmwf_quasi_cubic (                                             &
            ext_data,                                                          &
            dim_i_in, dim_j_in, dim_k_in,                                      &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j, number_of_inputs,                                  &
            weight_lambda_e,weight_phi_e,                                      &
            i_out_e, j_out_e, k_out_e,                                         &
            coeff_z_e,                                                         &
            data_out_high_e)

        ELSE IF (high_order_scheme  ==  ecmwf_mono_quasicubic) THEN

          ! DEPENDS ON: ecmwf_mono_quasi_cubic
          CALL ecmwf_mono_quasi_cubic (                                        &
            ext_data,                                                          &
            dim_i_in, dim_j_in, dim_k_in,                                      &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j, number_of_inputs,                                  &
            weight_lambda_e,weight_phi_e,                                      &
            i_out_e, j_out_e, k_out_e,                                         &
            coeff_z_e,                                                         &
            data_out_high_e)

        ELSE IF (high_order_scheme  ==  hcubic_vlin) THEN


          ! DEPENDS ON: h_cubic_v_linear
          CALL h_cubic_v_linear (                                              &
            ext_data,                                                          &
            dim_i_in, dim_j_in, dim_k_in,                                      &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j, number_of_inputs,                                  &
            weight_lambda_e,weight_phi_e,                                      &
            i_out_e, j_out_e, k_out_e,                                         &
            coeff_z_lin_e,                                                     &
            data_out_high_e)

        ELSE IF (high_order_scheme  ==  hquasicubic_vquintic) THEN

          ! DEPENDS ON: h_quasi_cubic_v_quintic
          CALL h_quasi_cubic_v_quintic (                                       &
            ext_data,                                                          &
            dim_i_in, dim_j_in, dim_k_in,                                      &
            dim_e_out, 1,1,                                                    &
            halo_i, halo_j, number_of_inputs,                                  &
            weight_lambda_e, weight_phi_e,                                     &
            i_out_e, j_out_e, k_out_e,                                         &
            coeff_z_e,                                                         &
            data_out_high_e)

        ELSE IF (high_order_scheme  ==  hcubic_vquintic) THEN

          ! DEPENDS ON: h_cubic_v_quintic
          CALL h_cubic_v_quintic (                                             &
            ext_data,                                                          &
            dim_i_in, dim_j_in, dim_k_in,                                      &
            dim_e_out, 1,1,                                                    &
            halo_i, halo_j, number_of_inputs,                                  &
            weight_lambda_e, weight_phi_e,                                     &
            i_out_e, j_out_e, k_out_e,                                         &
            coeff_z_e,                                                         &
            data_out_high_e)



        END IF

      END IF

      IF (((l_high .AND. l_mono) .OR. l_conserv)                               &
        .AND. high_order_scheme  /=  ecmwf_mono_quasicubic) THEN

        ! DEPENDS ON: mono_enforce
        CALL mono_enforce(                                                     &
          ext_data, number_of_inputs,                                          &
          dim_i_in, dim_j_in, dim_k_in,                                        &
          dim_e_out, 1, 1,                                                     &
          halo_i, halo_j,                                                      &
          i_out_e, j_out_e, k_out_e,                                           &
          data_out_high_e)

      END IF

      ! Call monotone scheme if required.

      IF ( (l_mono .AND. .NOT. l_high ) .OR. l_conserv ) THEN

        IF (monotone_scheme  ==  trilinear) THEN

          ! DEPENDS ON: tri_linear
          CALL tri_linear (                                                    &
            ext_data,                                                          &
            dim_i_in, dim_j_in, dim_k_in,                                      &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j, number_of_inputs,                                  &
            weight_lambda_e,weight_phi_e,                                      &
            i_out_e, j_out_e, k_out_e,                                         &
            coeff_z_lin_e,                                                     &
            data_out_mono_e)

        ELSE IF (monotone_scheme  ==  mono_quasicubic) THEN

          ! DEPENDS ON: ecmwf_mono_quasi_cubic
          CALL ecmwf_mono_quasi_cubic (                                        &
            ext_data,                                                          &
            dim_i_in, dim_j_in, dim_k_in,                                      &
            dim_e_out, 1, 1,                                                   &
            halo_i, halo_j, number_of_inputs,                                  &
            weight_lambda_e,weight_phi_e,                                      &
            i_out_e, j_out_e, k_out_e,                                         &
            coeff_z_e,                                                         &
            data_out_mono_e)

        END IF

      END IF

    END IF ! (dim_e_out  >   0)

    ! Time to return the points I have computed on demand.
  
    nsend = 0

    IF (l_2dcomm) THEN

      DO i = 0, n_proc-1
        IF (n_recvfrom(i)  >   0) THEN
          IF (l_high .AND. l_mono .AND. l_conserv) THEN
            len = 2*n_recvfrom(i)*number_of_inputs
            DO k = 1, number_of_inputs
              DO j = 1, n_recvfrom(i)
                send_data(2*(k-1)*n_recvfrom(i)+2*(j-1)+1,i)                   &
                  = data_out_mono_e(nsend+j+(k-1)*dim_e_out)
                send_data(2*(k-1)*n_recvfrom(i)+2*(j-1)+2,i)                   &
                  = data_out_high_e(nsend+j+(k-1)*dim_e_out)
              END DO
            END DO
            CALL gc_rsend(40*(me+1)+i, len, i, info,                           &
              recv_data(1,me), send_data(1,i))
            nsend = nsend + n_recvfrom(i)
          ELSE IF (l_high) THEN
            len = n_recvfrom(i)*number_of_inputs
            DO k = 1, number_of_inputs
              DO j = 1, n_recvfrom(i)
                send_data((k-1)*n_recvfrom(i)+j,i)                             &
                  = data_out_high_e(nsend+j+(k-1)*dim_e_out)
              END DO
            END DO
            CALL gc_rsend(40*(me+1)+i, len, i, info,                           &
              recv_data(1,me), send_data(1,i))
            nsend = nsend + n_recvfrom(i)
          ELSE
            len = n_recvfrom(i)*number_of_inputs
            DO k = 1, number_of_inputs
              DO j = 1, n_recvfrom(i)
                send_data((k-1)*n_recvfrom(i)+j,i)                             &
                  = data_out_mono_e(nsend+j+(k-1)*dim_e_out)
              END DO
            END DO
            CALL gc_rsend(40*(me+1)+i, len, i, info,                           &
              recv_data(1,me), send_data(1,i))
            nsend = nsend + n_recvfrom(i)
          END IF
        END IF
      END DO

      CALL gcg_ssync(proc_row_group, info)

      DO i = 0, n_proc-1
        IF (n_sendto(i)  >   0) THEN
          IF (l_high .AND. l_mono .AND. l_conserv) THEN
            len = 2*n_sendto(i)*number_of_inputs
          ELSE
            len = n_sendto(i)*number_of_inputs
          END IF

          CALL gc_rrecv(40*(i+1)+me, len, i, info,                             &
            recv_data(1,i), send_data(1,me))


          IF (l_high .AND. l_mono .AND. l_conserv) THEN
            DO n = 1, number_of_inputs
              DO j = 1, n_sendto(i)
                data_out_mono(i_store(j,i), j_store(j,i),                      &
                  k_store(j,i), n)                                             &
                  = recv_data(2*(n-1)*n_sendto(i)+2*(j-1)+1,i)
                data_out_high(i_store(j,i), j_store(j,i),                      &
                  k_store(j,i), n)                                             &
                  = recv_data(2*(n-1)*n_sendto(i)+2*(j-1)+2,i)
              END DO
            END DO
          ELSE IF ( l_high ) THEN
            DO n = 1, number_of_inputs
              DO j = 1, n_sendto(i)
                data_out_high(i_store(j,i), j_store(j,i),                      &
                  k_store(j,i), n)                                             &
                  = recv_data((n-1)*n_sendto(i)+j,i)
              END DO
            END DO
          ELSE
            DO n = 1, number_of_inputs
              DO j = 1, n_sendto(i)
                data_out_mono(i_store(j,i), j_store(j,i),                      &
                  k_store(j,i), n)                                             &
                  = recv_data((n-1)*n_sendto(i)+j,i)
              END DO
            END DO
          END IF
        END IF
      END DO

    ELSE   ! .not. L_2dcomm

      DO i = 0, n_procx-1
        IF (n_recvfrom(i)  >   0) THEN
          IF (l_high .AND. l_mono .AND. l_conserv) THEN
            len = 2*n_recvfrom(i)*number_of_inputs
            DO k = 1, number_of_inputs
              DO j = 1, n_recvfrom(i)
                send_data(2*(k-1)*n_recvfrom(i)+2*(j-1)+1,i)                   &
                  = data_out_mono_e(nsend+j+(k-1)*dim_e_out)
                send_data(2*(k-1)*n_recvfrom(i)+2*(j-1)+2,i)                   &
                  = data_out_high_e(nsend+j+(k-1)*dim_e_out)
              END DO
            END DO
            CALL gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,               &
              recv_data(1,ime), send_data(1,i))
            nsend = nsend + n_recvfrom(i)
          ELSE IF (l_high) THEN
            len = n_recvfrom(i)*number_of_inputs
            DO k = 1, number_of_inputs
              DO j = 1, n_recvfrom(i)
                send_data((k-1)*n_recvfrom(i)+j,i)                             &
                  = data_out_high_e(nsend+j+(k-1)*dim_e_out)
              END DO
            END DO
            CALL gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,               &
              recv_data(1,ime), send_data(1,i))
            nsend = nsend + n_recvfrom(i)
          ELSE
            len = n_recvfrom(i)*number_of_inputs
            DO k = 1, number_of_inputs
              DO j = 1, n_recvfrom(i)
                send_data((k-1)*n_recvfrom(i)+j,i)                             &
                  = data_out_mono_e(nsend+j+(k-1)*dim_e_out)
              END DO
            END DO
            CALL gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,               &
              recv_data(1,ime), send_data(1,i))
            nsend = nsend + n_recvfrom(i)
          END IF
        END IF
      END DO

      CALL gcg_ssync(proc_row_group, info)

      DO i = 0, n_procx-1
        IF (n_sendto(i)  >   0) THEN
          IF (l_high .AND. l_mono .AND. l_conserv) THEN
            len = 2*n_sendto(i)*number_of_inputs
          ELSE
            len = n_sendto(i)*number_of_inputs
          END IF

          CALL gc_rrecv(40*(ibase+i+1)+me, len, ibase+i, info,                 &
            recv_data(1,i), send_data(1,ime))

          IF (l_high .AND. l_mono .AND. l_conserv) THEN
            DO n = 1, number_of_inputs
              DO j = 1, n_sendto(i)
                data_out_mono(i_store(j,i), j_store(j,i),                      &
                  k_store(j,i), n)                                             &
                  = recv_data(2*(n-1)*n_sendto(i)+2*(j-1)+1,i)
                data_out_high(i_store(j,i), j_store(j,i),                      &
                  k_store(j,i), n)                                             &
                  = recv_data(2*(n-1)*n_sendto(i)+2*(j-1)+2,i)
              END DO
            END DO
          ELSE IF ( l_high ) THEN
            DO n = 1, number_of_inputs
              DO j = 1, n_sendto(i)
                data_out_high(i_store(j,i), j_store(j,i),                      &
                  k_store(j,i), n)                                             &
                  = recv_data((n-1)*n_sendto(i)+j,i)
              END DO
            END DO
          ELSE
            DO n = 1, number_of_inputs
              DO j = 1, n_sendto(i)
                data_out_mono(i_store(j,i), j_store(j,i),                      &
                  k_store(j,i), n)                                             &
                  = recv_data((n-1)*n_sendto(i)+j,i)
              END DO
            END DO
          END IF
        END IF
      END DO

    END IF  ! l_2dcomm



    ! Now distribute the pole values.

    IF (l_2dcomm) THEN
      IF (model_domain == mt_global) THEN
        IF (pole_handling  ==  1) THEN
          IF (at_extremity(psouth)) THEN
            DO n = 1,number_of_inputs
              DO k = 1, dim_k_out
                DO i = 2, dim_i_out
                  data_out_high(i,1,k,n) = data_out_high(1,1,k,n)
                  data_out_mono(i,1,k,n) = data_out_mono(1,1,k,n)
                END DO
              END DO
            END DO
          END IF
          IF (at_extremity(pnorth)) THEN
            DO n = 1,number_of_inputs
              DO k = 1, dim_k_out
                DO i = 2, dim_i_out
                  data_out_high(i,dim_j_out,k,n) = data_out_high(1,dim_j_out,k,n)
                  data_out_mono(i,dim_j_out,k,n) = data_out_mono(1,dim_j_out,k,n)
                END DO
              END DO
            END DO
          END IF

        END IF ! (pole_handling  ==  1)
      END IF ! (model_domain  ==  mt_Global)

    ELSE  ! .not. L-2dcomm
      IF (pole_handling  ==  1) THEN

        ! DEPENDS ON: interpolation_distribute_poles
        CALL interpolation_distribute_poles(                                   &
          sp_send, sp_levels, np_send, np_levels,                              &
          number_of_inputs, n_procx, ibase, proc_row_group,                    &
          ime, dim_i_out, dim_j_out, dim_k_out,                                &
          at_extremity, l_high, l_mono, l_conserv,                             &
          data_out_mono, data_out_high)

      END IF ! (pole_handling  ==  1)
    END IF  !l_2dcomm

  END IF ! n_proc>1 .and.l_2dcomm .or.
  ! n_procx>1 .and. model_domain==mt_Global .and. .not.l_2dcomm

  ! ----------------------------------------------------------------------
  ! Section 5.   Perform conservation if desired.
  !              Put interpolated field into Data_out.
  ! ----------------------------------------------------------------------

  IF ( l_high .AND. l_mono .AND. l_conserv ) THEN

    lambda_start = datastart(1) - halo_i

    IF ( itype == 1 ) THEN  ! for monotone u

      ! DEPENDS ON: mono_conserv
      CALL mono_conserv (                                                      &
        data_out_high(1,1,1,1),                                                &
        data_out_mono(1,1,1,1),                                                &
        ext_data(1-halo_i,1-halo_j,-1,1), r_in,                                &
        delta_r_in,                                                            &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        delta_lambda, delta_phi,                                         &
        row_length, n_rows, 0, 1,                                              &
        gdlambda_p(lambda_start), dphi_v,                                      &
        l_regular,                                                             &
        cos_latitude, off_x, off_y,                                            &
        n_proc, halo_i, halo_j,                                                &
        me, proc_row_group, proc_col_group,                                    &
        halo_data_out_i, halo_data_out_j,                                      &
        data_out1, conserv_fail)

      IF ( conserv_fail .AND. me == 0 ) THEN
        PRINT*,' non-conservation for field 1.'
      END IF

      IF (number_of_inputs >= 2) THEN


        ! DEPENDS ON: mono_conserv
        CALL mono_conserv (                                                    &
          data_out_high(1,1,1,2),                                              &
          data_out_mono(1,1,1,2),                                              &
          ext_data(1-halo_i,1-halo_j,-1,2), r_in,                              &
          delta_r_in,                                                          &
          dim_i_in, dim_j_in, dim_k_in,                                        &
          dim_i_out, dim_j_out, dim_k_out,                                     &
          delta_lambda, delta_phi,                                       &
          row_length, n_rows, 0, 1,                                            &
          gdlambda_p(lambda_start), dphi_v,                                    &
          l_regular,                                                           &
          cos_latitude, off_x, off_y,                                          &
          n_proc, halo_i, halo_j,                                              &
          me, proc_row_group, proc_col_group,                                  &
          halo_data_out_i, halo_data_out_j,                                    &
          data_out2, conserv_fail)

        IF ( conserv_fail .AND. me == 0 ) THEN
          PRINT*,' non-conservation for field 2.'
        END IF

      END IF  !  number_of_inputs >= 2

      IF (number_of_inputs >= 3) THEN


        ! DEPENDS ON: mono_conserv
        CALL mono_conserv (                                                    &
          data_out_high(1,1,1,3),                                              &
          data_out_mono(1,1,1,3),                                              &
          ext_data(1-halo_i,1-halo_j,-1,3), r_in,                              &
          delta_r_in,                                                          &
          dim_i_in, dim_j_in, dim_k_in,                                        &
          dim_i_out, dim_j_out, dim_k_out,                                     &
          delta_lambda, delta_phi,                                       &
          row_length, n_rows, 0, 1,                                            &
          gdlambda_p(lambda_start), dphi_v,                                    &
          l_regular,                                                           &
          cos_latitude, off_x, off_y,                                          &
          n_proc, halo_i, halo_j,                                              &
          me, proc_row_group, proc_col_group,                                  &
          halo_data_out_i, halo_data_out_j,                                    &
          data_out3, conserv_fail)

        IF ( conserv_fail .AND. me == 0 ) THEN
          PRINT*,' non-conservation for field 3.'
        END IF

      END IF  !  number_of_inputs >= 3

    ELSE IF ( itype == 2 ) THEN  ! for monotone v


      ! DEPENDS ON: mono_conserv
      CALL mono_conserv (                                                      &
        data_out_high(1,1,1,1),                                                &
        data_out_mono(1,1,1,1),                                                &
        ext_data(1-halo_i,1-halo_j,-1,1), r_in,                                &
        delta_r_in,                                                            &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        delta_lambda, delta_phi,                                         &
        row_length, rows, 1, 0,                                                &
        gdlambda_u(lambda_start), dphi_p,                                      &
        l_regular,                                                             &
        cos_latitude, off_x, off_y,                                            &
        n_proc, halo_i, halo_j,                                                &
        me, proc_row_group, proc_col_group,                                    &
        halo_data_out_i, halo_data_out_j,                                      &
        data_out1, conserv_fail)

      IF ( conserv_fail .AND. me == 0 ) THEN
        PRINT*,' non-conservation for field 1.'
      END IF

      IF (number_of_inputs >= 2) THEN

        ! DEPENDS ON: mono_conserv
        CALL mono_conserv (                                                    &
          data_out_high(1,1,1,2),                                              &
          data_out_mono(1,1,1,2),                                              &
          ext_data(1-halo_i,1-halo_j,-1,2), r_in,                              &
          delta_r_in,                                                          &
          dim_i_in, dim_j_in, dim_k_in,                                        &
          dim_i_out, dim_j_out, dim_k_out,                                     &
          delta_lambda, delta_phi,                                       &
          row_length, rows, 1, 0,                                              &
          gdlambda_u(lambda_start), dphi_p,                                    &
          l_regular,                                                           &
          cos_latitude, off_x, off_y,                                          &
          n_proc, halo_i, halo_j,                                              &
          me, proc_row_group, proc_col_group,                                  &
          halo_data_out_i, halo_data_out_j,                                    &
          data_out2, conserv_fail)

        IF ( conserv_fail .AND. me == 0 ) THEN
          PRINT*,' non-conservation for field 2.'
        END IF

      END IF  !  number_of_inputs >= 2

      IF (number_of_inputs >= 3) THEN

        ! DEPENDS ON: mono_conserv
        CALL mono_conserv (                                                    &
          data_out_high(1,1,1,3),                                              &
          data_out_mono(1,1,1,3),                                              &
          ext_data(1-halo_i,1-halo_j,-1,3), r_in,                              &
          delta_r_in,                                                          &
          dim_i_in, dim_j_in, dim_k_in,                                        &
          dim_i_out, dim_j_out, dim_k_out,                                     &
          delta_lambda, delta_phi,                                       &
          row_length, rows, 1, 0,                                              &
          gdlambda_u(lambda_start), dphi_p,                                    &
          l_regular,                                                           &
          cos_latitude, off_x, off_y,                                          &
          n_proc, halo_i, halo_j,                                              &
          me, proc_row_group, proc_col_group,                                  &
          halo_data_out_i, halo_data_out_j,                                    &
          data_out3, conserv_fail)

        IF ( conserv_fail .AND. me == 0 ) THEN
          PRINT*,' non-conservation for field 3.'
        END IF

      END IF  !  number_of_inputs >= 3

    ELSE                       ! for monotone p-fields

      ! DEPENDS ON: mono_conserv
      CALL mono_conserv (                                                      &
        data_out_high(1,1,1,1),                                                &
        data_out_mono(1,1,1,1),                                                &
        ext_data(1-halo_i,1-halo_j,-1,1), r_in,                                &
        delta_r_in,                                                            &
        dim_i_in, dim_j_in, dim_k_in,                                          &
        dim_i_out, dim_j_out, dim_k_out,                                       &
        delta_lambda, delta_phi,                                         &
        row_length, n_rows, 1, 1,                                              &
        gdlambda_u(lambda_start), dphi_v,                                      &
        l_regular,                                                             &
        cos_latitude, off_x, off_y,                                            &
        n_proc, halo_i, halo_j,                                                &
        me, proc_row_group, proc_col_group,                                    &
        halo_data_out_i, halo_data_out_j,                                      &
        data_out1, conserv_fail)

      IF ( conserv_fail .AND. me  ==  0 ) THEN
        PRINT*,' non-conservation for field 1.'
      END IF

      IF (number_of_inputs  >=  2) THEN

        ! DEPENDS ON: mono_conserv
        CALL mono_conserv (                                                    &
          data_out_high(1,1,1,2),                                              &
          data_out_mono(1,1,1,2),                                              &
          ext_data(1-halo_i,1-halo_j,-1,2), r_in,                              &
          delta_r_in,                                                          &
          dim_i_in, dim_j_in, dim_k_in,                                        &
          dim_i_out, dim_j_out, dim_k_out,                                     &
          delta_lambda, delta_phi,                                       &
          row_length, n_rows, 1, 1,                                            &
          gdlambda_u(lambda_start), dphi_v,                                    &
          l_regular,                                                           &
          cos_latitude, off_x, off_y,                                          &
          n_proc, halo_i, halo_j,                                              &
          me, proc_row_group, proc_col_group,                                  &
          halo_data_out_i, halo_data_out_j,                                    &
          data_out2, conserv_fail)

        IF ( conserv_fail .AND. me  ==  0 ) THEN
          PRINT*,' non-conservation for field 2.'
        END IF

      END IF  !  number_of_inputs  >=  2

      IF (number_of_inputs  >=  3) THEN

        ! DEPENDS ON: mono_conserv
        CALL mono_conserv (                                                    &
          data_out_high(1,1,1,3),                                              &
          data_out_mono(1,1,1,3),                                              &
          ext_data(1-halo_i,1-halo_j,-1,3), r_in,                              &
          delta_r_in,                                                          &
          dim_i_in, dim_j_in, dim_k_in,                                        &
          dim_i_out, dim_j_out, dim_k_out,                                     &
          delta_lambda, delta_phi,                                       &
          row_length, n_rows, 1, 1,                                            &
          gdlambda_u(lambda_start), dphi_v,                                    &
          l_regular,                                                           &
          cos_latitude, off_x, off_y,                                          &
          n_proc, halo_i, halo_j,                                              &
          me, proc_row_group, proc_col_group,                                  &
          halo_data_out_i, halo_data_out_j,                                    &
          data_out3, conserv_fail)

        IF ( conserv_fail .AND. me  ==  0 ) THEN
          PRINT*,' non-conservation for field 3.'
        END IF

      END IF !  number_of_inputs  >=  3

    END IF ! itype == 1

  ELSE IF ( l_high ) THEN

    DO k = 1, dim_k_out
      DO j = 1, dim_j_out
        DO i = 1, dim_i_out
          data_out1 (i,j,k) = data_out_high (i,j,k,1)
        END DO
      END DO
    END DO

    IF (number_of_inputs  >=  2) THEN
      DO k = 1, dim_k_out
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            data_out2 (i,j,k) = data_out_high (i,j,k,2)
          END DO
        END DO
      END DO
    END IF

    IF (number_of_inputs  >=  3) THEN
      DO k = 1, dim_k_out
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            data_out3 (i,j,k) = data_out_high (i,j,k,3)
          END DO
        END DO
      END DO
    END IF

  ELSE

    DO k = 1, dim_k_out
      DO j = 1, dim_j_out
        DO i = 1, dim_i_out
          data_out1 (i,j,k) = data_out_mono (i,j,k,1)
        END DO
      END DO
    END DO

    IF (number_of_inputs  >=  2) THEN
      DO k = 1, dim_k_out
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            data_out2 (i,j,k) = data_out_mono (i,j,k,2)
          END DO
        END DO
      END DO
    END IF

    IF (number_of_inputs  >=  3) THEN
      DO k = 1, dim_k_out
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            data_out3 (i,j,k) = data_out_mono (i,j,k,3)
          END DO
        END DO
      END DO
    END IF

  END IF

END IF          ! End if statement on Error Code

IF (lhook) CALL dr_hook('INTERPOLATION',zhook_out,zhook_handle)

RETURN
END SUBROUTINE interpolation

