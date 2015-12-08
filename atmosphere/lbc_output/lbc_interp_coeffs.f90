! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates the Horizontal Interpolation Coefficients for LBCs.
!
! Subroutine Interface:

      SUBROUTINE LBC_Interp_Coeffs (                                    &
           lbc_size                                                     &
      ,    src_row_len                                                  &
      ,    src_rows                                                     &
      ,    lambda_source                                                &
      ,    phi_source                                                   &
      ,    src_pole_lat                                                 &
      ,    src_pole_long                                                &
      ,    src_cyclic                                                   &
      ,    src_rotated                                                  &
      ,    lbc_row_len                                                  &
      ,    lbc_rows                                                     &
      ,    same_rotation                                                &
      ,    l_eg_grid                                                    &
      ,    lambda_in                                                    &
      ,    phi_in                                                       &
      ,    lbc_delta_lat                                                &
      ,    lbc_delta_long                                               &
      ,    lbc_first_lat                                                &
      ,    lbc_first_long                                               &
      ,    lbc_pole_lat                                                 &
      ,    lbc_pole_long                                                &
      ,    rimwidth                                                     &
      ,    lbc_halo_x                                                   &
      ,    lbc_halo_y                                                   &
      ,    lbc_halo_x_coords                                            &
      ,    lbc_halo_y_coords                                            &
      ,    lbc_index_bl                                                 &
      ,    lbc_index_br                                                 &
      ,    lbc_weights_tr                                               &
      ,    lbc_weights_br                                               &
      ,    lbc_weights_bl                                               &
      ,    lbc_weights_tl                                               &
      ,    lbc_coeff1                                                   &
      ,    lbc_coeff2                                                   &
      ,    i_uv                                                         &
       )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      USE eqtoll_mod, ONLY: eqtoll
      USE h_int_co_mod, ONLY: h_int_co
      USE lltoeq_mod, ONLY: lltoeq
      USE w_coeff_mod, ONLY: w_coeff
      USE lam_inclusion_mod, ONLY: lam_inclusion
      IMPLICIT NONE
!
! Description:
!   Calculates the Horizontal Interpolation Coefficients for LBCs.
!
! Method:
!   1. Calculate lat/longs of lbc points on rotated grid.
!   2. Call EQTOLL to get corresponding true lat/longs.
!   3. Call W_COEFF to calculate coefficients to rotate the winds.
!   4. Set up lat/longs for the source grid.
!   5. Call H_INT_CO to calculate the horizontal interpolation
!      coefficients.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables :

!*L------------------ COMDECK CINTFA ----------------------------------
!L CMAXSIZE should be called first.
!
!   Contains Variables, Headers and Index blocks for control of
!   generation of boundary information for the limited area model.
!
!   Interfaces to all other models are handled by STASH, and there is
!   no explicit coding written for them in the model.
!
! Interface variables initialised through INTFCNSTA
! Namelist read in the interface control routine INTF_CTL.

      INTEGER                                                           &
        intf_row_length                                                 &
                         ! Interface field row length
       ,intf_p_rows                                                     &
                         ! Interface field no of rows
       ,intf_p_levels                                                   &
                         ! Interface field no of levels
       ,intf_q_levels                                                   &
                         ! Interface field no of wet levels
       ,intf_tr_levels                                                  &
                         ! Interface field no of tracer levels
       ,intfwidtha                                                      &
                         ! Width of interface zone (atmosphere)
       ,intf_exthalo_ns                                                 &
                         ! Extended Halo in NS direction
       ,intf_exthalo_ew                                                 &
                         ! Extended Halo in EW direction
       ,a_intf_start_hr                                                 &
                         ! ) Start and End time in
       ,a_intf_freq_hr                                                  &
                         ! ) hours, Frequency in h,m,s for which
       ,a_intf_freq_mn                                                  &
                         ! ) atmosphere interface data
       ,a_intf_freq_sc                                                  &
                         ! ) is to be generated.
       ,a_intf_end_hr                                                   &
                         ! )
       ,intf_pack                                                       &
                         ! Packing Indicator for boundary data
       ,lbc_stream_a                                                    &
                         ! Output streams in UMUI
       ,lbc_unit_no_a                                                   &
                         ! Unit Nos for Atmos Boundary Dataset
       ,lbc_first_r_rho                                                 &
                         ! First rho level at which height is constant
       ,intf_v_int_order(max_n_intf_a)

      REAL                                                              &
        intf_ewspace                                                    &
                         ! E-W grid spacing (degrees)
       ,intf_nsspace                                                    &
                         ! N-S grid spacing (degrees)
       ,intf_firstlat                                                   &
                         ! Latitude of first row (degrees)
       ,intf_firstlong                                                  &
                         ! Longitude of first row (degrees)
       ,intf_polelat                                                    &
                         ! Real latitude of coordinate pole (degrees)
       ,intf_polelong                                                   &
                         ! Real longitude of coordinate pole (degrees)
       ,lbc_z_top_model                                                 &
                         ! Height of top of model
       ,lbc_q_min                                                       &
                         ! Minimum value for q
!
! VarRes grid spacing
      , lambda_intf_p(max_intf_lbcrow_length, max_n_intf_a)             &
      , lambda_intf_u(max_intf_lbcrow_length, max_n_intf_a)             &    
      , phi_intf_p(max_intf_lbcrows, max_n_intf_a)                      &
      , phi_intf_v(max_intf_lbcrows, max_n_intf_a)

      LOGICAL                                                           &
        intf_vert_interp                                                &
                         ! Switch to request vertical interpolation
       ,lnewbnd          ! True for initialising new boundary data file

! Switch for variable resolution LBC output
      LOGICAL  intf_l_var_lbc(max_n_intf_a)

! Switch to not rotate if input and output grids have same poles.
      LOGICAL intf_avoid_rot(MAX_N_INTF_A)

! Switch to output LBC for Endgame
      LOGICAL intf_l_eg_out(MAX_N_INTF_A)

! Files for VERTLEVS namelist     
      CHARACTER(LEN=256) :: intf_vertlevs

! Files for HorzGrid namelist  
      CHARACTER(LEN=256) :: intf_HorzGrid(max_n_intf_a)
!*----------------------------------------------------------------------
      COMMON /INTFCTL_ATMOS/                                            &
        intf_ewspace(max_n_intf_a)    ,intf_nsspace(max_n_intf_a),      &
        intf_firstlat(max_n_intf_a)   ,intf_firstlong(max_n_intf_a),    &
        intf_polelat(max_n_intf_a)    ,intf_polelong(max_n_intf_a),     &
        intf_row_length(max_n_intf_a) ,intf_p_rows(max_n_intf_a),       &
        intf_p_levels(max_n_intf_a)   ,intf_q_levels(max_n_intf_a),     &
        intf_tr_levels(max_n_intf_a)  ,intfwidtha(max_n_intf_a),        &
        intf_exthalo_ns(max_n_intf_a) ,intf_exthalo_ew(max_n_intf_a),   &
        a_intf_start_hr(max_n_intf_a) ,a_intf_freq_hr(max_n_intf_a),    &
        a_intf_freq_mn(max_n_intf_a)  ,a_intf_freq_sc(max_n_intf_a),    &
        a_intf_end_hr(max_n_intf_a)   ,                                 & 
        lnewbnd(max_n_intf_a)         ,intf_vert_interp(max_n_intf_a),  &
        intf_pack(max_n_intf_a)       ,lbc_stream_a(max_n_intf_a),      &
        lbc_unit_no_a(max_n_intf_a)   ,lbc_first_r_rho(max_n_intf_a),   &
        lbc_z_top_model(max_n_intf_a) ,                                 &
        intf_vertlevs(max_n_intf_a)   ,lbc_q_min,                       &
        intf_l_var_lbc                ,intf_horzgrid,                   &
        lambda_intf_p                 ,lambda_intf_u,                   &
        phi_intf_p                    ,phi_intf_v,                      &
        intf_avoid_rot                ,intf_v_int_order,                &
        intf_l_eg_out
!---------------------------------------------------------------------

! Subroutine arguments

      Integer :: lbc_size

!     Source Grid
      Integer :: src_row_len
      Integer :: src_rows
      Real    :: src_pole_lat
      Real    :: src_pole_long
      Logical :: src_cyclic
      Logical :: src_rotated

! Arrays containing the source grid
      REAL :: lambda_source(src_row_len)
      REAL :: phi_source(src_rows)

!     LBC Grid
      Integer :: lbc_row_len
      Integer :: lbc_rows
      Integer :: lbc_halo_x
      Integer :: lbc_halo_y
      INTEGER :: lbc_halo_x_coords
      INTEGER :: lbc_halo_y_coords
      Real    :: lbc_delta_lat
      Real    :: lbc_delta_long
      Real    :: lbc_first_lat
      Real    :: lbc_first_long
      Real    :: lbc_pole_lat
      Real    :: lbc_pole_long
      
      Logical same_rotation
! Input VarRes grid info in degrees 
      REAL    :: Lambda_in(1-lbc_halo_x_coords:lbc_row_len+lbc_halo_x_coords)  
      REAL    :: Phi_in   (1-lbc_halo_y_coords:lbc_rows+lbc_halo_y_coords)
      
      Integer :: RimWidth
      Integer, dimension (lbc_size) :: lbc_index_bl
      Integer, dimension (lbc_size) :: lbc_index_br
      Real,    dimension (lbc_size) :: lbc_weights_tr
      Real,    dimension (lbc_size) :: lbc_weights_br
      Real,    dimension (lbc_size) :: lbc_weights_tl
      Real,    dimension (lbc_size) :: lbc_weights_bl
      Real,    dimension (lbc_size) :: lbc_coeff1
      Real,    dimension (lbc_size) :: lbc_coeff2

      Integer :: i_uv   !  If 1 => u, 2 => v, Otherwise 0.
      
      LOGICAL :: l_eg_grid

! Local parameters:

      Character (Len=*), Parameter :: RoutineName = 'LBC_Interp_Coeffs'

! Local scalars:

      Integer :: ipt         ! LBC point number
      Integer :: row,pt      ! Loop indices for row and point
      Integer :: iside       ! Loop index for LBC sides
      Integer :: lbc_len     ! Computed no of lbc points
      Integer :: ErrorStatus ! Error Code
      ! Offsets due to u/v LBC grid.
      INTEGER ::  v_extra_nt   &
                , v_extra_nb   &
                , v_extra_et   &
                , v_extra_eb   &
                , v_extra_st   &
                , v_extra_sb   &
                , v_extra_wt   &
                , v_extra_wb   &
                , u_extra_nl   &
                , u_extra_nr   &
                , u_extra_el   &
                , u_extra_er   &
                , u_extra_sl   &
                , u_extra_sr   &
                , u_extra_wl   &
                , u_extra_wr

      Character (Len=80) :: CMessage

! Local dynamic arrays:

      Real, dimension (:), allocatable :: lambda_lbc
      Real, dimension (:), allocatable :: phi_lbc
      Real, dimension (:), allocatable :: lambda_targ
      Real, dimension (:), allocatable :: phi_targ

      REAL :: phi_targ_min, phi_targ_max

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header

      IF (lhook) CALL dr_hook('LBC_INTERP_COEFFS',zhook_in,zhook_handle)
      ErrorStatus = 0
      CMessage = ' '

! If u/v field we need to take into extra point in LBC
! Names denote location: n, e, s, w - north, east, south, west
!                        t, b, l, r - top, bottom, left, right
! E.g. v_extra_nt: how many extra rows for interpolation to V
! field are required at the top of the north section of the LBC?
      IF (i_uv > 0) THEN
        IF (l_eg_grid) THEN
          v_extra_nt = 1
          v_extra_nb = 0
          v_extra_et = 1
          v_extra_eb = 1
          v_extra_st = 0
          v_extra_sb = 1
          v_extra_wt = 1
          v_extra_wb = 1
          u_extra_nl = 1
          u_extra_nr = 0
          u_extra_el = 1
          u_extra_er = 0
          u_extra_sl = 1
          u_extra_sr = 0
          u_extra_wl = 1
          u_extra_wr = 0
        ELSE
          v_extra_nt = 0
          v_extra_nb = 1
          v_extra_et = 0
          v_extra_eb = 0
          v_extra_st = 1
          v_extra_sb = 0
          v_extra_wt = 0
          v_extra_wb = 0
          u_extra_nl = 0
          u_extra_nr = 0
          u_extra_el = 1
          u_extra_er = 0
          u_extra_sl = 0
          u_extra_sr = 0
          u_extra_wl = 0
          u_extra_wr = 1
        END IF
      ELSE
          v_extra_nt = 0
          v_extra_nb = 0
          v_extra_et = 0
          v_extra_eb = 0
          v_extra_st = 0
          v_extra_sb = 0
          v_extra_wt = 0
          v_extra_wb = 0
          u_extra_nl = 0
          u_extra_nr = 0
          u_extra_el = 0
          u_extra_er = 0
          u_extra_sl = 0
          u_extra_sr = 0
          u_extra_wl = 0
          u_extra_wr = 0
      END IF

! --------------------
! Allocate work arrays
! --------------------

      ALLOCATE ( lambda_lbc  (lbc_size) )
      ALLOCATE ( phi_lbc     (lbc_size) )

! ---------------------------
! Get lat/Longs of LBC points
! ---------------------------

      ipt=0  
   
      DO ISide = 1, 4 
 
        IF (ISide == PSouth) THEN   !  Southern Boundary (SP)  
  
          DO Row = 1 - lbc_halo_y - v_extra_sb, RimWidth + v_extra_st  
            DO Pt = 1 - lbc_halo_x - u_extra_sl, lbc_row_len + lbc_halo_x + u_extra_sr  
   
              ipt = ipt + 1  
              Lambda_lbc(ipt) = Lambda_in(pt) 
              Phi_lbc(ipt)    = Phi_in(row)  

            END DO 
          END DO

        END IF  !  South  
  
        IF (ISide == PEast) THEN   !  Eastern Boundary   
  
          DO Row = RimWidth + 1 - v_extra_eb, lbc_rows - RimWidth + v_extra_et  
            DO Pt = lbc_row_len - u_extra_el - RimWidth + 1,                &
                  lbc_row_len + lbc_halo_x + u_extra_er  
  
              ipt = ipt + 1  
              Lambda_lbc(ipt) = Lambda_in(pt) 
              Phi_lbc(ipt)    = Phi_in(row)      

            END DO          
          END DO            
   
        END IF  ! East   
 
        IF (ISide == PNorth) THEN  !  Northern Boundary (NP)  
          DO Row = lbc_rows - RimWidth + 1 - v_extra_nb,                    &
                 lbc_rows + lbc_halo_y + v_extra_nt  
            DO Pt = 1-lbc_halo_x - u_extra_nl, lbc_row_len+lbc_halo_x+u_extra_nr  
  
            ipt = ipt + 1 
              Lambda_lbc(ipt) = Lambda_in(pt) 
              Phi_lbc(ipt)    = phi_in(row)             

            END DO             
          END DO
 
        END IF  ! North 
  
        IF (ISide == PWest) THEN  !  Western Boundary  
  
          DO Row = RimWidth + 1 - v_extra_wb, lbc_rows - RimWidth + v_extra_wt
            DO Pt = 1 - lbc_halo_x - u_extra_wl, RimWidth + u_extra_wr 
    
              ipt = ipt + 1  
              Lambda_lbc(ipt) = Lambda_in(pt) 
              Phi_lbc(ipt)    = Phi_in(row) 

            END DO
          END DO
  
        END IF  ! West 
  
      END DO !  ISide 

      lbc_len = ipt

! -----------------------------------------------------
! Check no of lbc points computed ; must match lbc_size
! -----------------------------------------------------

      IF (lbc_len /= lbc_size) THEN
        WRITE (6,*) ' Mismatch in number of LBC points.'
        WRITE (6,*) ' Expected no of lbc points ',lbc_size
        WRITE (6,*) ' Computed no of lbc points ',lbc_len
        WRITE (CMessage,*) ' Mismatch in number of LBC points.'
        ErrorStatus = 10

        CALL Ereport ( RoutineName, ErrorStatus, CMessage)
      END IF

! -------------------------------------
! Get true lat and longs for lbc points
! -------------------------------------

      IF ( .NOT. same_rotation ) THEN

        ALLOCATE ( lambda_targ (lbc_size) )
        ALLOCATE ( phi_targ    (lbc_size) )

        CALL EqToLL (                                                   &
                     phi_lbc,  lambda_lbc,                              &
                     phi_targ, lambda_targ,                             &
                     lbc_pole_lat, lbc_pole_long, lbc_len )

! ----------------------------------------------
! Calculate the coefficients to rotate the winds
! ----------------------------------------------

     ! This is only needed when v field is calculated.
        IF (i_uv == 2) THEN

          CALL W_Coeff (                                                &
                        lbc_coeff1, lbc_coeff2,                         &
                        lambda_targ, lambda_lbc,                        &
                        lbc_pole_lat, lbc_pole_long, lbc_len )

        END IF  !  i_uv == 2

! -----------------------------------------------------------
! For rotated model grids, get lat/longs w.r.t the model grid
! -----------------------------------------------------------

        IF (src_rotated) THEN

          CALL LLToEq (                                                 &
             phi_targ                                                   &
      ,      lambda_targ                                                &
      ,      phi_targ                                                   &
      ,      lambda_targ                                                &
      ,      src_pole_lat                                               &
      ,      src_pole_long                                              &
      ,      lbc_len                                                    &
       )
        END IF
      END IF ! same_rotation

! Now we check LBC is inside the source grid.
! Only do the checking if src_cyclic is false, i.e. we are NOT a global model
      IF (.NOT. src_cyclic) THEN
        IF (.NOT. same_rotation) THEN
          ! As not same rotation use rotated phi coords.
          phi_targ_min = MINVAL(phi_targ)
          phi_targ_max = MAXVAL(phi_targ)
        ELSE
          ! If same rotation use the phi_lbc directly
          phi_targ_min = MINVAL(phi_lbc)
          phi_targ_max = MAXVAL(phi_lbc)
        END IF

        CALL lam_inclusion(src_rows, src_row_len, 0, 0,           &
          lambda_source, phi_source, src_pole_lat, src_pole_long, &
          lbc_rows, lbc_row_len, lbc_halo_x, lbc_halo_y,          &
          lambda_in, phi_in, phi_targ_min, phi_targ_max,          &
          lbc_pole_lat, lbc_pole_long,                            &
          same_rotation, src_rotated)

      END IF ! .NOT. src_cyclic

! --------------------------------------------------
! Calculate the horizontal interpolation cofficients
! --------------------------------------------------

      IF ( same_rotation ) THEN
        Call H_Int_Co (                                                 &
         lbc_index_bl,                                                  &
         lbc_index_br,                                                  &
         lbc_weights_tr                                                 &
      ,  lbc_weights_br                                                 &
      ,  lbc_weights_tl                                                 &
      ,  lbc_weights_bl                                                 &
      ,  Lambda_Source                                                  &
      ,  Phi_Source                                                     &
      ,  Lambda_lbc                                                     &
      ,  Phi_lbc                                                        &
      ,  src_row_len                                                    &
      ,  src_rows                                                       &
      ,  lbc_len                                                        &
      ,  src_cyclic                                                     &
       )
      ELSE
        Call H_Int_Co (                                                 &
         lbc_index_bl,                                                  &
         lbc_index_br,                                                  &
         lbc_weights_tr                                                 &
      ,  lbc_weights_br                                                 &
      ,  lbc_weights_tl                                                 &
      ,  lbc_weights_bl                                                 &
      ,  Lambda_Source                                                  &
      ,  Phi_Source                                                     &
      ,  Lambda_Targ                                                    &
      ,  Phi_Targ                                                       &
      ,  src_row_len                                                    &
      ,  src_rows                                                       &
      ,  lbc_len                                                        &
      ,  src_cyclic                                                     &
       )
       DEALLOCATE (Lambda_Targ)
       DEALLOCATE (Phi_Targ)
        
      END IF ! same rotation
 

      deallocate (phi_lbc)
      deallocate (lambda_lbc)

      IF (lhook) CALL dr_hook('LBC_INTERP_COEFFS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LBC_Interp_Coeffs
