! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Sets up lookup-table for boundary file
!
! Subroutine Interface:

      SUBROUTINE LBC_SetUp_LookUp (                                     &
     &           lbc_lookUp,                                            &
     &           len1_lbc_lookup,                                       &
     &           n_lbc_vars,                                            &
     &           lbc_item_codes,                                        &
     &           lbc_stash_codes,                                       &
     &           jintf,                                                 &
     &           ntime,                                                 &
     &           lbc_fixhd,                                             &
     &           lbc_rim_size,                                          &
     &           lbc_rim_type,                                          &
     &           lbc_halo_type,                                         &
     &           lbc_fld_type,                                          &
     &           intf_halosize,                                         &
     &           lbc_levels                                             &
     & )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE io_configuration_mod, ONLY : io_field_padding
      USE UM_ParVars
      USE Control_Max_Sizes
      USE lbc_mod
      USE c_model_id_mod, ONLY: model_id
      USE lookup_addresses
      USE Submodel_Mod

      USE cppxref_mod, ONLY: ppx_dump_packing

      IMPLICIT NONE
!
! Description:
!   Set up the lookup headers for lbc variables for one area.
!   Called for each validity time of lbc data.
!
! Method:
!   Set up lookup-table as defined in UMDP F3.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output
!
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):

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

! Subroutine arguments

      Integer  :: len1_lbc_lookup
      Integer  :: n_lbc_vars
      Integer  :: lbc_lookup(len1_lbc_Lookup, n_lbc_vars)
      Integer  :: jintf
      Integer  :: ntime
      Integer  :: lbc_fixhd(256)
      Integer  :: lbc_item_codes(n_lbc_vars)
      Integer  :: lbc_stash_codes(n_lbc_vars)      
      Integer  :: lbc_rim_size  (n_lbc_vars)
      Integer  :: lbc_rim_type  (n_lbc_vars)
      Integer  :: lbc_halo_type (n_lbc_vars)
      Integer  :: lbc_fld_type  (n_lbc_vars)
      Integer  :: lbc_levels    (n_lbc_vars)
      Integer  :: intf_halosize (2,Nhalo_max)

! Local parameters:

!      Integer,           Parameter :: Sect31   = 31
      Integer,           Parameter :: Sect32   = 32
      Character (Len=*), Parameter :: RoutineName= 'LBC_SetUp_LookUp'

!     Local variables

      Integer :: Start_Address   !  w.r.t start of data
      Integer :: Disk_Address    !  w.r.t start of file
      Integer :: Len_data        !  actual length of data
      Integer :: Disk_length     !  length of data including rounding
                                 !  to sector boundary
      Integer :: N1,N2,N3,N4,N5  !  components of packing indicator
      Integer :: R1R0            !  lbc rimwidth
      Integer :: Y1Y0            !  lbc halo - NS
      Integer :: X1X0            !  lbc halo - EW
      Integer :: Var, Var1       !  loop indexes
      Integer :: lbc_item        !  Item Code

      Real    :: Lat0, Lon0      ! Zeroth Lat & Long

      Integer :: ErrorStatus          !  error code
      Character (Len=80) :: Cmessage  !  error message

! Local dynamic arrays:

!     The lookup table is not stored for previous LBC data times so
!     retain disk and start addresses to avoid unnecessary i/o.


! Function & Subroutine calls:
      Integer Exppxi, Get_um_version_id

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header


!     Initialise disk address and start address

!     - If LBCs already exist in the LBC file then LBC_LOOKUP on entry
!     contains the lookup headers for the last batch of LBCs generated.

!     - At the start of CRUNs LBC_LOOKUP will contain the lookup header
!     corresponding to the last LBC variable before the CRUN restart
!     position. (Lookup Header read in in IN_INTF)

      IF (lhook) CALL dr_hook('LBC_SETUP_LOOKUP',zhook_in,zhook_handle)
      If (ntime == 1) Then
        disk_address  = lbc_fixhd(160)-1
        start_address = 1
      Else
        disk_address  = LBC_Lookup(lbegin,n_lbc_vars) +                 &
     &                  LBC_Lookup(lbnrec,n_lbc_vars)
        start_address = LBC_Lookup(naddr ,n_lbc_vars) +                 &
     &                  LBC_Lookup(lblrec,n_lbc_vars)
      End If

      If (ntime == 1) Then  !  Include Orography
        Var1 = 1
      Else                  !  Exclude Orography
        Var1 = 2
      Endif

      LBC_Lookup(:,:) = 0

      Do Var = Var1, n_lbc_vars

        lbc_item = lbc_item_codes(var)

!       Determine length of data for this lbc variable

        len_data =                                                      &
     &  lbc_global_lenrima(lbc_fld_type(var),lbc_halo_type(var) ) *     &
     &  lbc_levels(var)

!       Validity Time
        LBC_Lookup (lbyr ,var) = LBC_VT_Year
        LBC_Lookup (lbmon,var) = LBC_VT_Month
        LBC_Lookup (lbdat,var) = LBC_VT_Day
        LBC_Lookup (lbhr ,var) = LBC_VT_Hour
        LBC_Lookup (lbmin,var) = LBC_VT_Min
        LBC_Lookup (lbsec,var) = LBC_VT_Sec

!       Data Time
        LBC_Lookup (lbyrd ,var) = LBC_DT_Year
        LBC_Lookup (lbmond,var) = LBC_DT_Month
        LBC_Lookup (lbdatd,var) = LBC_DT_Day
        LBC_Lookup (lbhrd ,var) = LBC_DT_Hour
        LBC_Lookup (lbmind,var) = LBC_DT_Min
        LBC_Lookup (lbsecd,var) = LBC_DT_Sec

!       Length of data
        LBC_Lookup (lblrec,var) = LEN_DATA

!       Grid Type - Rotated regular lat/long
        LBC_Lookup (lbcode,var) = 100+1

!       Hemisphere Indicator : Not used for LBC files
!       Use to store no of lbc levels : 100+(no of levels)
        LBC_Lookup (lbhem,var) = 100+lbc_levels(var)

        IF (intf_l_eg_out(jintf)) THEN
        
        IF (lbc_fld_type(var) == fld_type_v) THEN
          LBC_LOOKUP(lbrow,var) = intf_p_rows(jintf)+1
        ELSE
          LBC_LOOKUP(lbrow,var) = intf_p_rows(jintf)
        END IF
        
        LBC_LOOKUP(lbnpt,var) = intf_row_length(jintf)
        
        ELSE ! EndGame Lookup

!       No of rows - not including haloes
        If (lbc_fld_type(var) == fld_type_v) Then  ! v-comp
          LBC_Lookup (lbrow,var) = intf_p_rows(jintf)-1
        Else
          LBC_Lookup (lbrow,var) = intf_p_rows(jintf)
        End If

!       No of points - not including haloes
        If (lbc_fld_type(var) == fld_type_u) Then  ! u-comp
          LBC_Lookup (lbnpt,var) = intf_row_length(jintf)-1
        Else
          LBC_Lookup (lbnpt,var) = intf_row_length(jintf)
        End If

        END IF

!       Packing Indicator
        Select Case ( intf_pack(jintf) )

          Case (0)   !  No packing
            N1 = 0

          Case (1)   !  32-bit packing
            N1 = 2

          Case (2)   !  Pack according to Stashmaster Record
! DEPENDS ON: exppxi
            N1 = Exppxi(                                                &
     &           atmos_im,Sect32,lbc_item,ppx_dump_packing,             &
     &           ErrorStatus, cmessage)

        End Select

        N2 = 0   !  Data not compressed
        N3 = 0   !  Compression definition
        N4 = 0   !  Number format
        N5 = 0   !  Not used
        LBC_Lookup (lbpack,var) =                                       &
     &  N5*10000 + N4*1000 +N3*100 + N2*10 + N1

        LBC_Lookup (lbrel,var) = 3

!       Disk Address
        LBC_Lookup (lbegin,var) = disk_address

!       Disk Length

!       Get length of data
        if (mod(lbc_lookup(lbpack, var), 10) == 2) then
          disk_length= (LBC_Lookup(lblrec,var)+1)/2
        else
          disk_length = LBC_Lookup(lblrec,var)
        endif

!       round up to sector boundary
        disk_length=((disk_length+io_field_padding-1)/                    &
     &               io_field_padding)*io_field_padding

!       Disk length : includes rounding-up to sector boundary
        LBC_Lookup (lbnrec,var) = disk_length

!       Level code : 7777 for multi-level field.
        LBC_Lookup (lblev,var) = 7777

! DEPENDS ON: get_um_version_id
        LBC_Lookup (lbsrce,var) = get_um_version_id(model_id)

!       Data Type
        LBC_Lookup (data_type,var) = 1

!       Start Address
        LBC_Lookup (naddr,var) = Start_Address

!       LBUSER(3) stores rimwidth and haloes as R1R0Y1Y0X1X0
!       R1R0 : Rimwidth, Y1Y0 : halo_y, X1X0 : halo_x

        R1R0 = lbc_rim_size(lbc_rim_type(var))
        Y1Y0 = Intf_HaloSize(2,lbc_halo_type(var))
        X1X0 = Intf_HaloSize(1,lbc_halo_type(var))
        LBC_Lookup (lbuser3,var) = R1R0 * 10000 + Y1Y0 * 100 + X1X0

!       Item Code
        LBC_Lookup (item_code,var) = lbc_stash_codes(var)

!       Model Code
        LBC_Lookup (model_code,var) = atmos_im

!       REAL section of lookup.

!       Lat & Long of rotated pole
        LBC_Lookup (bplat,var) = Transfer (intf_polelat (jintf),1)
        LBC_Lookup (bplon,var) = Transfer (intf_polelong(jintf),1)

! Check whether fixed or variable resolution.
        IF (intf_ewspace(jintf) < 0) THEN
          LBC_Lookup (bzy,var) = TRANSFER(rmdi, 1)
          LBC_Lookup (bdy,var) = TRANSFER(rmdi, 1)
          LBC_Lookup (bzx,var) = TRANSFER(rmdi, 1)
          LBC_Lookup (bdx,var) = TRANSFER(rmdi, 1)
        ELSE
!       Zeroth latitude
          Lat0 = LBC_Grid(lbc_fld_type(var), 1) - intf_nsspace(jintf)
          LBC_Lookup (bzy,var) = TRANSFER (Lat0, 1)
  
!       Latitude Interval
          LBC_Lookup (bdy,var) = TRANSFER (intf_nsspace(jintf),1)

!       Zeroth longitude
          Lon0 = LBC_Grid(lbc_fld_type(var), 2) - intf_ewspace(jintf)
          LBC_Lookup (bzx,var) = TRANSFER (Lon0, 1)

!       Longitude Interval
          LBC_Lookup (bdx,var) = TRANSFER (intf_ewspace(jintf), 1)
        END IF


!       Update disk and start address
        disk_address  = disk_address  + LBC_Lookup(lbnrec,var)
        start_address = start_address + LBC_Lookup(lblrec,var)

!       Assume always real data (could be checked using data_type (set to 1
!       above which is real data anyway)
        LBC_Lookup (bmdi,var) = Transfer (rmdi,1)

      End Do
      
        Do Var = Var1, n_lbc_vars
        End Do
      IF (lhook) CALL dr_hook('LBC_SETUP_LOOKUP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LBC_SetUp_LookUp
