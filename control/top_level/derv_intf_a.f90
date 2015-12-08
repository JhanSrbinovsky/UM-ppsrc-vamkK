! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine DERV_INTF_A : Calculates Interface array dimensions.
!
! Subroutine Interface :
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level
      SUBROUTINE derv_intf_a (max_intf_model_levels,         &
                              max_lbcrow_length,max_lbcrows, &
                              n_intf_a)

      USE check_iostat_mod
      USE ereport_mod, ONLY : ereport
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Control_Max_Sizes
      USE printstatus_mod
      IMPLICIT NONE
!
! Description : Calculate array dimensions for boundary data output.
!
! Method : Reads in INTFCNSTA namelist to get grid dimensions of
!          interface area. Calculates array dimensions for boundary
!          data. Also sets dimensions to 1 if no interface areas
!          required to prevent zero dynamic allocation.
!
! Code Description :
! Language : FORTRAN 77 + common extensions
!
! Declarations :

!     Arguments
      Integer MAX_INTF_MODEL_LEVELS ! OUT  Max no of lbc levels
      Integer MAX_LBCROW_LENGTH     ! OUT  Max no of lbc row_length
      Integer MAX_LBCROWS           ! OUT  Max no of lbc rows
      Integer N_INTF_A          ! IN   No of interface areas

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

!  Namelist for atmos interface constants
!+ COMDECK CNAMINFA
!
!    Description:
!       This COMDECK contains the INTFCNSTA namelist which
!       contains all the variables required to define the grids
!       of the interface areas for Atmosphere Boundary data.
!
!       All variables are set up in the UMUI.
!       All variables are declared in comdeck CINTFA

      NAMELIST/INTFCNSTA/                                               &
     &         INTF_EWSPACE,INTF_NSSPACE,INTF_FIRSTLAT,INTF_FIRSTLONG,  &
     &         INTF_POLELAT,INTF_POLELONG,                              &
     &         INTF_ROW_LENGTH,INTF_P_ROWS,INTF_P_LEVELS,INTF_Q_LEVELS, &
     &         INTF_TR_LEVELS,                                          &
     &         INTFWIDTHA, Intf_ExtHalo_NS, Intf_ExtHalo_EW,            &
     &         Intf_Pack, LBC_Q_MIN,                                    &
     &         A_INTF_FREQ_HR,A_INTF_FREQ_MN,A_INTF_FREQ_SC,            &
     &         A_INTF_START_HR,A_INTF_END_HR                            &
     &        ,LBC_Stream_A, Intf_VertLevs                              &
     &        ,INTF_L_VAR_LBC, INTF_HorzGrid                            &
     &        ,intf_avoid_rot, intf_v_int_order, intf_l_eg_out
 
      ! -----------------------------
      ! OLDVERT Namelist for 4.5 LBCs
      ! -----------------------------

      Logical Vert_Interp
      Integer Meth_Lev_Calc
      Integer Max_sig_hlev, Min_prs_hlev
      Real Etah (model_levels_max+1)

      Namelist /OLDVERT/ Vert_Interp, Meth_Lev_Calc,                    &
                         Max_sig_hlev, Min_prs_hlev, etah

!- End of comdeck CNAMINFA

!     Local variables
      INTEGER JINTF  !  loop index
      INTEGER    ::   errorstatus
      INTEGER    ::   icode
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      CHARACTER(LEN=256)            :: cmessage
      CHARACTER(LEN=*), PARAMETER   ::  RoutineName = 'derv_intf_a'

!     Read in INTFCNSTA namelist to get output grids for
!     generating boundary data.

      IF (lhook) CALL dr_hook('DERV_INTF_A',zhook_in,zhook_handle)
      REWIND 5
      READ (UNIT=5, NML=INTFCNSTA, IOSTAT=ErrorStatus) 
      CALL check_iostat(errorstatus, "namelist INTFCNSTA")
      REWIND 5

! Find the number of atmosphere interface areas by counting the non-zero
! stream numbers of the LBC streams
      N_INTF_A=0
      DO JINTF=1,MAX_N_INTF_A
        IF(LBC_STREAM_A(JINTF) /= 0) N_INTF_A = N_INTF_A + 1
      END DO

      IF (N_INTF_A >  0) THEN

!       Boundary data to be generated in this run.

        MAX_INTF_MODEL_LEVELS = 0
        MAX_LBCROW_LENGTH     = 0
        MAX_LBCROWS           = 0
        DO JINTF=1,N_INTF_A

          MAX_INTF_MODEL_LEVELS =                                       &
              MAX ( MAX_INTF_MODEL_LEVELS , INTF_P_LEVELS(JINTF)+1 )
          MAX_LBCROW_LENGTH =                                           &
              MAX ( MAX_LBCROW_LENGTH , INTF_ROW_LENGTH(JINTF)+1 )
          MAX_LBCROWS = MAX ( MAX_LBCROWS , INTF_P_ROWS(JINTF)+1 )  
        END DO

!       Check >= 1 to avoid zero dynamic allocation
        MAX_INTF_MODEL_LEVELS = MAX ( MAX_INTF_MODEL_LEVELS , 1)
        MAX_LBCROW_LENGTH     = MAX ( MAX_LBCROW_LENGTH , 1)
        MAX_LBCROWS           = MAX ( MAX_LBCROWS , 1)
        
        IF (PrintStatus >= PrStatus_Diag) THEN
          WRITE (6,'(A,I9)') ' n_intf_a              ',n_intf_a
          WRITE (6,'(A,I9)') ' max_intf_model_levels ',max_intf_model_levels
          WRITE (6,'(A,I9)') ' max_lbcrow_length     ',max_lbcrow_length
          WRITE (6,'(A,I9)') ' max_lbcrows           ',max_lbcrows
        END IF

        ! We are generating LBCs - this is not supported within a UM
        ! run so report a warning.
        icode = -1
        WRITE(cmessage,'(A)')                                          &
          'Generating LBCs from within a UM run is not supported.  '// &
          'Please use makebc.'
        CALL ereport(routinename,icode,cmessage)

      ELSE

!       No boundary conditions to be generated.
!       Initialise to prevent zero length dynamic allocation.

        IF (PrintStatus >= PrStatus_Diag) THEN
          WRITE (6,'(A,I3)') ' old n_intf_a ',n_intf_a
        END IF

        N_INTF_A = 1
        MAX_INTF_MODEL_LEVELS = 1
        MAX_LBCROW_LENGTH = 1
        MAX_LBCROWS = 1

        IF (PrintStatus >= PrStatus_Diag) THEN
          WRITE (6,'(A,I3)') ' n_intf_a ',n_intf_a
          WRITE (6,'(A,I4)') ' max_intf_model_levels ',max_intf_model_levels
          WRITE (6,'(A,I5)') ' max_lbcrow_length ',max_lbcrow_length
          WRITE (6,'(A,I5)') ' max_lbcrows ',max_lbcrows
        END IF

      END IF

      IF (lhook) CALL dr_hook('DERV_INTF_A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DERV_INTF_A
