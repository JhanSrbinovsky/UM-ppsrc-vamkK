! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Performs transforms after a LBC field has been interpolated
!
! Subroutine Interface:

      Subroutine LBC_Post_Interp_Transf (                               &
                 lbc_size                                               &
      ,          lbc_first_level                                        &
      ,          lbc_last_level                                         &
      ,          lbc_sectioncode                                        &
      ,          lbc_stashcode                                          &
      ,          lbc_data                                               &
      ,          L_extra_surface_level                                  &
      ,          L_vi)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      IMPLICIT NONE
!
! Description:
!   Perform tranformations/processing on a LBC field after it has
!   been interpolated.
!
! Method:
!   Choice of transform/processing is based on stashcode. Data passes
!   through unchanged if no processing done.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

! Arguments
      Integer  ::  lbc_size
      Integer  ::  lbc_first_level
      Integer  ::  lbc_last_level
      Integer  ::  lbc_stashcode
      Integer  ::  lbc_sectioncode

      Real  ::  lbc_data (lbc_size, lbc_first_level:lbc_last_level)
      
      LOGICAL :: L_extra_surface_level
      LOGICAL :: L_vi

!Local Parameters

      Integer,           Parameter :: Sect32   = 32

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
! Include file : parlbcs.h
!
! Must be called after parvars.h
!
! Description:
!   Contains variables in connection with generating LBCs.
!
! -----------------------------------------------------------
! Stash Codes for LBCs in Section 32 (and section 31 except tracers) 
!
      Integer, Parameter :: lbc_stashcode_orog    = 1 
      Integer, Parameter :: lbc_stashcode_u       = 2 
      Integer, Parameter :: lbc_stashcode_v       = 3 
      Integer, Parameter :: lbc_stashcode_w       = 4 
      Integer, Parameter :: lbc_stashcode_density = 5 
      Integer, Parameter :: lbc_stashcode_theta   = 6 
      Integer, Parameter :: lbc_stashcode_q       = 7 
      Integer, Parameter :: lbc_stashcode_qcl     = 8 
      Integer, Parameter :: lbc_stashcode_qcf     = 9 
      Integer, Parameter :: lbc_stashcode_exner   = 10 
      Integer, Parameter :: lbc_stashcode_u_adv   = 11 
      Integer, Parameter :: lbc_stashcode_v_adv   = 12 
      Integer, Parameter :: lbc_stashcode_w_adv   = 13 
      Integer, Parameter :: lbc_stashcode_qcf2    = 14 
      Integer, Parameter :: lbc_stashcode_qrain   = 15 
      Integer, Parameter :: lbc_stashcode_qgraup  = 16 
      Integer, Parameter :: lbc_stashcode_cf_bulk = 17 
      Integer, Parameter :: lbc_stashcode_cf_liquid = 18 
      Integer, Parameter :: lbc_stashcode_cf_frozen = 19 
      Integer, Parameter :: lbc_stashcode_murk      = 20 
      Integer, Parameter :: lbc_stashcode_free_tracer = 21 
      Integer, Parameter :: lbc_stashcode_ukca_tracer = 22 
      Integer, Parameter :: lbc_stashcode_dust_div1 = 23
      Integer, Parameter :: lbc_stashcode_dust_div2 = 24
      Integer, Parameter :: lbc_stashcode_dust_div3 = 25
      Integer, Parameter :: lbc_stashcode_dust_div4 = 26
      Integer, Parameter :: lbc_stashcode_dust_div5 = 27
      Integer, Parameter :: lbc_stashcode_dust_div6 = 28
      Integer, Parameter :: lbc_stashcode_so2      = 29
      Integer, Parameter :: lbc_stashcode_dms      = 30
      Integer, Parameter :: lbc_stashcode_so4_aitken = 31
      Integer, Parameter :: lbc_stashcode_so4_accu = 32
      Integer, Parameter :: lbc_stashcode_so4_diss = 33
      Integer, Parameter :: lbc_stashcode_nh3      = 35
      Integer, Parameter :: lbc_stashcode_soot_new = 36
      Integer, Parameter :: lbc_stashcode_soot_agd = 37
      Integer, Parameter :: lbc_stashcode_soot_cld = 38
      Integer, Parameter :: lbc_stashcode_bmass_new = 39
      Integer, Parameter :: lbc_stashcode_bmass_agd = 40
      Integer, Parameter :: lbc_stashcode_bmass_cld = 41
      Integer, Parameter :: lbc_stashcode_ocff_new = 42
      Integer, Parameter :: lbc_stashcode_ocff_agd = 43
      Integer, Parameter :: lbc_stashcode_ocff_cld = 44
      Integer, Parameter :: lbc_stashcode_nitr_acc = 45
      Integer, Parameter :: lbc_stashcode_nitr_diss = 46

! -----------------------------------------------------------
!     Data Time for LBC data
      Integer :: LBC_DT_Year
      Integer :: LBC_DT_Month
      Integer :: LBC_DT_Day
      Integer :: LBC_DT_Hour
      Integer :: LBC_DT_Min
      Integer :: LBC_DT_Sec
      Integer :: LBC_DT_DayNo

      COMMON /LBC_DT/ LBC_DT_Year, LBC_DT_Month, LBC_DT_Day,            &
         LBC_DT_Hour, LBC_DT_Min,  LBC_DT_Sec,   LBC_DT_DayNo

! -----------------------------------------------------------

!     Validity Time for LBC data
      Integer :: LBC_VT_Year
      Integer :: LBC_VT_Month
      Integer :: LBC_VT_Day
      Integer :: LBC_VT_Hour
      Integer :: LBC_VT_Min
      Integer :: LBC_VT_Sec
      Integer :: LBC_VT_DayNo

      COMMON /LBC_VT/ LBC_VT_Year, LBC_VT_Month, LBC_VT_Day,            &
         LBC_VT_Hour, LBC_VT_Min,  LBC_VT_Sec,   LBC_VT_DayNo

! -----------------------------------------------------------

      Integer, Parameter :: P_Src_Grid = 2
      Integer, Parameter :: P_LBC_Grid = 4

!     1 : Start Latitude
!     2 : Start Longitude
!     3 : Row Length
!     4 : Rows

      Real :: Src_Grid (Nfld_max, P_Src_Grid)
      Real :: LBC_Grid (Nfld_max, P_LBC_Grid)

      COMMON /LBC_Grids/ Src_Grid, LBC_Grid

! -------------------------------------------------------------

      Integer :: LBC_Global_LenRimA (Nfld_max, Nhalo_max)
      Integer :: LBC_Interp_LenRimA (Nfld_max, Nhalo_max)

      COMMON /LBC_Sizes/ LBC_Global_LenRimA, LBC_Interp_LenRimA

! -------------------------------------------------------------

! Local
      Integer :: Level
      Integer :: i
      Integer :: ErrorStatus      ! Return code

      Character (Len=256)          :: CMESSAGE          ! Error message
      Character (Len=*), Parameter :: RoutineName =                     &
     &                                  'lbc_post_interp_transf'

      Logical :: negative_tr_value            ! true if tracer negative

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! -------------------------
! Do transforms as required
! -------------------------

      IF (lhook) CALL dr_hook('LBC_POST_INTERP_TRANSF',zhook_in,zhook_handle)
      If (lbc_sectioncode == Sect32) Then
      
        IF (L_extra_surface_level) THEN
      
          IF (L_vi) THEN
          
          ! Vertically interpolated field already has all levels in the right
          ! place. We copy the first level to the surface (surface
          ! data is extrapolated from first two levels in lbc_vert_interp).
          
            DO i = 1, lbc_size
               lbc_data(i,lbc_first_level) = lbc_data(i,lbc_first_level+1)
            END DO
          
          ELSE
          
          ! If there has been no vertical interpolation (i.e. the ND dump
          ! has the same level set as the EG LAM except for the extra surface
          ! level) then we need to restructure the data correctly by shifting
          ! everything up one level, leaving the surface as a copy of the
          ! first level.
      
            DO level = lbc_last_level, lbc_first_level + 1, -1
              DO i = 1, lbc_size
                 lbc_data(i, level) = lbc_data(i, level - 1)
              END DO
            END DO
            
          END IF ! L_vi
        
        END IF ! L_extra_surface_level

        Select Case ( lbc_stashcode )

          Case (lbc_stashcode_w, lbc_stashcode_w_adv)

            lbc_data(:,lbc_first_level) = 0.0
            lbc_data(:,lbc_last_level)  = 0.0

          Case (lbc_stashcode_q)

            Do level = lbc_first_level, lbc_last_level
              Do i = 1, lbc_size
                If ( lbc_data(i,level) < lbc_q_min ) Then
                  lbc_data(i,level) = lbc_q_min
                End If
              End Do
            End Do

          Case (lbc_stashcode_qcf, lbc_stashcode_qcl,                   &
     &          lbc_stashcode_qcf2, lbc_stashcode_qrain,                &
     &          lbc_stashcode_qgraup)

            Do level = lbc_first_level, lbc_last_level
              Do i = 1, lbc_size
                If ( lbc_data(i,level) < 0.0 ) Then
                  lbc_data(i,level) = 0.0
                End If
              End Do
            End Do

          Case (lbc_stashcode_cf_bulk, lbc_stashcode_cf_liquid,         &
     &          lbc_stashcode_cf_frozen)

            Do level = lbc_first_level, lbc_last_level
              Do i = 1, lbc_size
                If (lbc_data(i,level) < 0.0 ) Then
                  lbc_data(i,level) = 0.0
                End If
                If (lbc_data(i,level) > 1.0 ) Then
                  lbc_data(i,level) = 1.0
                End If
              End Do
            End Do

          Case (lbc_stashcode_murk)

          ! Reset murk aerosol if less than 0.1 after interpolation

            Do level = lbc_first_level, lbc_last_level
              Do i = 1, lbc_size
                If ( lbc_data(i,level) < 0.1 ) Then
                  lbc_data(i,level) = 0.1
                End If
              End Do
            End Do

          ! Ensure tracer values are not negative
          ! after interpolation
          Case (lbc_stashcode_free_tracer)

            Do level = lbc_first_level, lbc_last_level
              Do i = 1, lbc_size
                If ( lbc_data(i,level) < 0.0 ) Then
                  lbc_data(i,level) = 0.0
                End If
              End Do
            End Do

        Case (lbc_stashcode_dust_div1,lbc_stashcode_dust_div2,          &
     &        lbc_stashcode_dust_div3,lbc_stashcode_dust_div4,          &
     &        lbc_stashcode_dust_div5,lbc_stashcode_dust_div6,          &
     &        lbc_stashcode_so2,lbc_stashcode_dms,                      &
     &        lbc_stashcode_so4_aitken,lbc_stashcode_so4_accu,          &
     &        lbc_stashcode_so4_diss,                                   &
     &        lbc_stashcode_nh3,lbc_stashcode_soot_new,                 &
     &        lbc_stashcode_soot_agd,lbc_stashcode_soot_cld,            &
     &        lbc_stashcode_bmass_new,lbc_stashcode_bmass_agd,          &
     &        lbc_stashcode_bmass_cld,lbc_stashcode_ocff_new,           &
     &        lbc_stashcode_ocff_agd,lbc_stashcode_ocff_cld,            &
     &        lbc_stashcode_nitr_acc,lbc_stashcode_nitr_diss)

          ErrorStatus = 0
          negative_tr_value=.FALSE.

        ! Reset tracer if less than 0.0 after interpolation

          Do level = lbc_first_level, lbc_last_level
            Do i = 1, lbc_size
              If ( lbc_data(i,level) < 0.0 ) Then
                negative_tr_value=.TRUE.
                lbc_data(i,level) = 0.0
              End If
            End Do
          End Do
          If ( negative_tr_value ) Then
            WRITE (Cmessage,*) 'lbc_post_interp_transf: negative        &
     &          tracer concentration ', lbc_stashcode
                Errorstatus = -1
        
                Call Ereport(RoutineName,Errorstatus,Cmessage)
          End If


          ! Ensure UKCA tracer values are not negative
          ! after interpolation
          Case (lbc_stashcode_ukca_tracer)

            Do level = lbc_first_level, lbc_last_level
              Do i = 1, lbc_size
                If ( lbc_data(i,level) < 0.0 ) Then
                  lbc_data(i,level) = 0.0
                End If
              End Do
            End Do

        End Select

      End if

      IF (lhook) CALL dr_hook('LBC_POST_INTERP_TRANSF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LBC_Post_Interp_Transf
