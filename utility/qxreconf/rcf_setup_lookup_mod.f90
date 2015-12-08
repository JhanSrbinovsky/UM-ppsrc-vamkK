! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up the output dump lookup tables

MODULE Rcf_Setup_Lookup_Mod
  IMPLICIT NONE

!  Subroutine Rcf_Setup_Lookup - sets up lookups for output dump
!
! Description:
!   The lookup headers are filled in - with calculated addressing
!   and other namelist derived variables
!
! Method:
!    UMDP F3 defines the lookups.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

  SUBROUTINE Rcf_Setup_Lookup( Hdr_In, Hdr_Out )

    USE Rcf_UMhead_Mod, ONLY :                                                 &
        um_header_type
!!!    LenFixHd

!!!Use Rcf_headers_Mod, ONLY :&
!!!    FixHd

    USE Rcf_Exppx_mod, ONLY :                                                  &
        Rcf_Exppx

    USE Rcf_Ppx_Info_mod, ONLY :                                               &
        STM_Record_Type

    USE Submodel_Mod, ONLY :                                                   &
        Internal_model_list,                                                   &
        N_Internal_model

    USE Rcf_NRecon_Mod, ONLY :                                                 &
        Recondat_Node,                                                         &
        RecondatList,                                                          &
        DumpProgLevs

    USE Ereport_Mod, ONLY :                                                    &
        Ereport

    USE Rcf_Model_Mod, ONLY :                                                  &
        ZonAvOzone

    USE Rcf_HeadAddress_Mod, ONLY :                                            &
        FH_HorizGrid,         RC_LatSpacing,        RC_LongSpacing,            &
        RC_PoleLat,           RC_PoleLong,          FH_VertCoord,              &
        FH_VertCoord_CP,      RC_FirstLat,          LDC_ZseaTheta,             &
        LDC_CkTheta,          LDC_ZseaRho,          ldc_ckrho,                 &
        FH_GridStagger,       FH_GridStagger_Endgame,                          &
        FH_VTYear,            FH_VTSecond,                                     &
        FH_DTYear,            FH_DTSecond

    USE Rcf_generate_heights_mod, ONLY :                                       &
        height_gen_original,                                                   &
        height_gen_smooth

    USE Rcf_Readnl_horizont_Mod, ONLY :                                        &
        Iproj

    USE Rcf_Grid_Type_Mod, ONLY :                                              &
        Output_Grid

    USE Rcf_Items_Mod, ONLY :                                                  &
        num_items,          area_array,             item_array

    USE Rcf_Recon_Mod, ONLY :                                                  &
        Dump_Pack,                                                             &
        Var_Recon

    USE um_input_control_mod, ONLY :                                           &
        lcal360

    USE lbc_mod, ONLY :                                                        &
        Rimwidtha

    USE RimTypes

    USE Rcf_Level_Code_Mod, ONLY :                                             &
        Rcf_Level_Code

    USE Rcf_Stashcodes_Mod, ONLY :                                             &
        stashcode_u

    USE c_model_id_mod, ONLY: model_id
    USE lookup_addresses

    USE version_mod, ONLY : NSectP
    USE cppxref_mod, ONLY :                                                    &
        ppx_theta_level,      ppx_rho_level,       ppx_atm_ozone,              &
        ppx_atm_river,        ppx_atm_lbc_u,       ppx_atm_lbc_v,              &
        ppx_atm_lbc_theta,    ppx_atm_cuall,       ppx_atm_cvall

    IMPLICIT NONE

! Arguments
    TYPE (Um_Header_type), INTENT(IN) :: Hdr_In
    TYPE (Um_Header_type), Target     :: Hdr_Out

! Comdecks
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

! Local vars.
    INTEGER                :: int_val = 1      ! used for transfers
    INTEGER, POINTER       :: Lookup(:,:)
    TYPE (STM_Record_Type) :: STM_Record
    INTEGER                :: i, j, jj, k ! Loopers
    INTEGER                :: icount     ! Counter for calc. Lookup(40)
    INTEGER                :: sec_item
    INTEGER                :: item
    INTEGER                :: section
    INTEGER                :: model
    INTEGER                :: k_out
    INTEGER                :: whole      ! Whole packing indicator -
                                     ! Lookup(21) - see UMDP F3
    INTEGER                :: n1, n2, n3 ! Packing and compression parts
                                     ! of above
    INTEGER                :: length
    INTEGER                :: start_address
    INTEGER                :: n_levels
    INTEGER                :: n_plevels
    INTEGER                :: lev          ! level
    INTEGER                :: bot_lev      ! bottom level for field
    INTEGER                :: lblev_val    ! value for lblev
    INTEGER                :: area
    INTEGER                :: ppxref_grid_type
    INTEGER                :: ErrorStatus
    INTEGER                :: Area_Expand( Hdr_Out % Len2Lookup )
    CHARACTER (LEN=*), PARAMETER :: RoutineName='Rcf_Setup_Lookup'
    CHARACTER (LEN=80)     :: Cmessage
    REAL                   :: depth
    REAL                   :: level( Hdr_Out % Len1LevDepC )
    REAL                   :: zsea  !} height defining constants
    REAL                   :: Ck    !}
    REAL                   :: riv_bzx
    REAL                   :: riv_bzy
    REAL                   :: grid_offset ! To deal with ND/ENDGAME grids

! Function & Subroutine calls:
    INTEGER            :: get_um_version_id

!------------------------------------------------------------------
! First set all elements of lookup to 0 - except dates which
! come from the input lookup (there must be 1).  If reconfiguring
! for VAR then need to make sure date and time is from an
! instantaneous field (problem when time accumulated fields are at
! the beginning) therefore find where u velocity is.
!------------------------------------------------------------------
    Lookup => Hdr_Out % Lookup
    k = 1
    IF (Var_Recon) THEN
      DO i = 1, Hdr_In % Len2Lookup
        ! Lets assume we have stashcode_u available.
        IF (Hdr_In % Lookup(item_code, i) == stashcode_u) THEN
          k = i
          EXIT
        END IF
      END DO
    END IF

    ! Ideally we should copy the fixed header values across but some jobs seem
    ! to have different data times in the fixed headers than in the lookups.
    ! This should be investigated to see if the date and times are set
    ! correctly before making the change here.

    DO i = lbyr, lbsecd
      Lookup( i, : ) = Hdr_In % Lookup( i, k )
    END DO
    Lookup( lbsec, : )   = Hdr_In % FixHd( FH_DTSecond )
    Lookup( lbsecd, : )  = Hdr_In % FixHd( FH_VTSecond )
    Lookup( 13 : 45, : ) = 0
    Lookup( 46 : 64, : ) = TRANSFER( 0.0, int_val)

    k_out=1

! 5: Loop through prognostic items and initialise Lookup
    DO j=1,n_internal_model
      DO jj=0,NSectP
        recondat_node => RecondatList(j,jj)
        DO While (ASSOCIATED(recondat_node % recondat_info))
      ! 5.1: Extract addressing and number of levels from linked list
          sec_item      = recondat_node % recondat_info % sec_item
          n_levels      = recondat_node % recondat_info % rlevs
          length        = recondat_node % recondat_info % LEN
          start_address = recondat_node % recondat_info % raddress
          n_plevels     = recondat_node % recondat_info % rplevs

          area = 1
          DO k = 1, num_items
            IF ( sec_item == item_array(k) ) THEN
              area = area_array(k)
            END IF
          END DO

          icount=0
          DO  k=k_out,k_out+n_levels*n_plevels-1

            area_expand(k) = area

        ! 5.3.1: Initialise Lookup
            Lookup( item_code, k ) = sec_item
            Lookup( model_code,k )=internal_model_list(j)

            item    = MOD(sec_item,1000)
            section = (Lookup(item_code,k) - item)/1000
            model   = internal_model_list(j)
            STM_record = Rcf_Exppx( model, section, item )

        ! 5.3.2: Set addressing information
            Lookup( lblrec,k)=length/(n_levels*n_plevels)

            Lookup( naddr, k)=start_address+icount
            icount=icount+Lookup( lblrec,k)

        !-----------------------------------------------------------
        !        Calculate levels from level dependent constants
        !        Set levels for multi-level fields only
        !        Levels not set for single level fields
        !-----------------------------------------------------------
            IF (n_levels > 1) THEN

          ! Set the level
              CALL Rcf_Level_Code( STM_record % lb_code, bot_lev,              &
                  Output_Grid)

          ! Set level number. For ozone we only use top of the
          ! atmosphere levels. Otherwise we count from the surface.
              lev = MOD(k - k_out, N_levels) + bot_lev

              IF (lev == 0) THEN
                lblev_val = 9999                      ! surface
              ELSE
                lblev_val = lev
              END IF

              IF (STM_record % grid_type == ppx_atm_ozone) THEN
                lev = lev + Output_Grid % model_levels -                       &
                    Output_Grid % ozone_levels
              END IF

              lookup( lblev,k ) = lblev_val

              IF ( Output_Grid % height_gen_method == height_gen_smooth) THEN

                IF (STM_record % lbvc_code == 9 .OR.                           &
                    STM_record % lbvc_code ==65 ) THEN  ! Hybrid/Eta levels

                  IF (STM_record % lv_code == ppx_theta_level) THEN

                    IF (lev == Output_Grid % model_levels ) THEN ! Top level

                  ! Level is current theta level
                      zsea = Hdr_Out % LevDepC( lev+1, LDC_ZseaTheta )
                      Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkTheta )
                      Lookup( blev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhlev, k ) = TRANSFER( Ck, int_val )

                  ! Lower boundary is rho level below
                      zsea = Hdr_Out % LevDepC( lev, LDC_ZseaRho )
                      Ck   = Hdr_Out % LevDepC( lev, ldc_ckrho )
                      Lookup( brlev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhrlev, k ) = TRANSFER( Ck, int_val )

                  ! Upper boundary is rho level above - calculated
                      zsea  = 2.0 * Hdr_Out % LevDepC( lev+1, LDC_ZseaTheta)   &
                          - Hdr_Out % LevDepC( lev, LDC_ZseaRho )
                      Ck    = 0.0    ! above 1st constant rho level
                      Lookup( bulev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhulev, k ) = TRANSFER( Ck, int_val )

                    ELSE IF ( lev > Output_Grid % model_levels ) THEN

                      Lookup( bulev, lev )  = TRANSFER( 0., int_val )
                      Lookup( blev, lev )   = TRANSFER( 0., int_val )
                      Lookup( brlev, lev )  = TRANSFER( 0., int_val )
                      Lookup( bhulev, lev ) = TRANSFER( 0., int_val )
                      Lookup( bhlev, lev )  = TRANSFER( 0., int_val )
                      Lookup( bhrlev, lev ) = TRANSFER( 0., int_val )

                    ELSE

                  ! Level is current theta level
                      zsea = Hdr_Out % LevDepC( lev+1, LDC_ZseaTheta )
                      Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkTheta )
                      Lookup( blev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhlev, k ) = TRANSFER( Ck, int_val )

                  ! Upper boundary is rho level above
                      zsea = Hdr_Out % LevDepC( lev+1, LDC_ZseaRho )
                      Ck   = Hdr_Out % LevDepC( lev+1, ldc_ckrho )
                      Lookup( bulev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhulev, k ) = TRANSFER( Ck, int_val )

                  ! Set the lower boundary
                      IF (lev <= 1) THEN
                    ! Lowest theta level has a lower boundary of the
                    ! physical surface, as defined in the physics.
                        zsea = 0.0
                        Ck   = 1.0
                      ELSE
                    ! Lower boundary is rho level below
                        zsea = Hdr_Out % LevDepC( lev, LDC_ZseaRho )
                        Ck   = Hdr_Out % LevDepC( lev, ldc_ckrho )
                      END IF
                      Lookup( brlev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhrlev, k ) = TRANSFER( Ck, int_val )
                    END IF

                  ELSE IF (STM_record % lv_code == ppx_rho_level ) THEN
                    IF ( lev == Output_Grid % model_levels + 1) THEN

                  ! Level is top rho level - calculate
                      zsea = 2.0 * Hdr_Out % LevDepC( lev, LDC_ZseaTheta ) -   &
                          Hdr_Out % LevDepC( lev - 1, LDC_ZseaRho )
                      Ck   = 0.0      ! Above 1st const rho by definition
                      Lookup( blev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhlev, k ) = TRANSFER( Ck, int_val )

                  ! Upper level boundary same as level
                      Lookup( bulev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhulev, k ) = TRANSFER( Ck, int_val )

                  ! Lower boundary is top theta level
                      zsea = Hdr_Out % LevDepC( lev, LDC_ZseaTheta )
                      Ck   = Hdr_Out % LevDepC( lev, LDC_CkTheta )
                      Lookup( brlev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhrlev, k ) = TRANSFER( Ck, int_val )

                    ELSE

                  ! Level above is Theta level
                      zsea = Hdr_Out % LevDepC( lev+1, LDC_ZSeaTheta )
                      Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkTheta )
                      Lookup( bulev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhulev, k ) = TRANSFER( Ck, int_val )

                  ! Level is Rho level
                      zsea = Hdr_Out % LevDepC( lev, ldc_zsearho )
                      Ck   = Hdr_Out % LevDepC( lev, ldc_ckrho )
                      Lookup( blev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhlev, k ) = TRANSFER( Ck, int_val )

                  ! Level below is Theta level
                      IF (lev <=0) THEN  ! orography
                        zsea = 0.0
                        Ck   = 1.0
                      ELSE
                        zsea = Hdr_Out % LevDepC( lev, LDC_ZSeaTheta )
                        Ck   = Hdr_Out % LevDepC( lev, LDC_CkTheta )
                      END IF
                      Lookup( brlev, k )  = TRANSFER( zsea, int_val )
                      Lookup( bhrlev, k ) = TRANSFER( Ck, int_val )

                    END IF    ! number of levels
                  END IF      ! rho levels
                END IF        ! Hybrid Levels

              ELSE IF (Output_Grid % height_gen_method ==                      &
                  height_gen_original) THEN

            ! Alternative lbvc code=65 enabled in case adopted in future
                IF (STM_record % lbvc_code == 9 .OR.                           &
                    STM_record % lbvc_code ==65 )  THEN   ! Hybrid/Eta levels
                  Lookup( bhlev, k )  = TRANSFER( rmdi, int_val )
                  Lookup( bhulev, k ) = TRANSFER( rmdi, int_val )
                  Lookup( bhrlev, k ) = TRANSFER( rmdi, int_val )

                  IF (STM_record % lv_code == ppx_theta_level) THEN

                    IF ( lev == Output_Grid % model_levels ) THEN ! Top level
                      Lookup( bulev, k ) = TRANSFER(                           &
                          2.0 * Output_Grid % eta_theta_levels( lev ) -        &
                          Output_Grid % eta_rho_levels( lev ), int_val )
                      Lookup( blev, k ) = TRANSFER(                            &
                          Output_Grid % eta_theta_levels( Lev), int_val )
                      Lookup( brlev, k ) = TRANSFER(                           &
                          Output_Grid % eta_rho_levels  ( Lev), int_val )
                    ELSE IF ( lev > Output_Grid % model_levels ) THEN
                      Lookup( bulev, k ) = TRANSFER( 0., int_val )
                      Lookup( blev, k  ) = TRANSFER( 0., int_val )
                      Lookup( brlev,k )  = TRANSFER( 0., int_val )
                    ELSE
                      Lookup( bulev, k ) = TRANSFER(                           &
                          Output_Grid % eta_rho_levels( Lev + 1), int_val )
                      Lookup( blev, k ) = TRANSFER(                            &
                          Output_Grid % eta_theta_levels( Lev), int_val )

                      IF (lev <= 0) THEN
                        Lookup( brlev, k ) = TRANSFER( 0., int_val )
                      ELSE
                        Lookup( brlev, k ) = TRANSFER(                         &
                            Output_Grid % eta_rho_levels  ( Lev), int_val )
                      END IF
                    END IF

                  ELSE IF ( STM_record % lv_code == ppx_rho_level ) THEN

                    IF ( lev == Output_Grid % model_levels + 1) THEN
                      Lookup( bulev, k ) = TRANSFER(                           &
                          2.0 * Output_Grid % eta_theta_levels(lev - 1) -      &
                          Output_Grid % eta_rho_levels( lev - 1), int_val )
                      Lookup( blev, k  ) = Lookup( bulev, k )
                      Lookup( brlev,k )  = TRANSFER(                           &
                          Output_Grid % eta_theta_levels( lev - 1), int_val)
                    ELSE
                      Lookup( bulev, k ) = TRANSFER(                           &
                          Output_Grid % eta_theta_levels( Lev), int_val )
                      Lookup( blev, k ) = TRANSFER(                            &
                          Output_Grid % eta_rho_levels  ( Lev), int_val )

                      IF ( lev <= 0 ) THEN   ! bottom level
                        Lookup( brlev, k ) = TRANSFER( 1.0, int_val )
                      ELSE
                        Lookup( brlev, k ) = TRANSFER(                         &
                            Output_Grid % eta_theta_levels( Lev - 1), int_val)
                      END IF
                    END IF
                  ELSE
                    WRITE (6,*) 'ERROR:- lv_code=',STM_record % lv_code
                    WRITE (6,*) 'model, section, item = ', model, section, item
                    ErrorStatus = 10
                    Cmessage = 'lv_code not right from STASHmaster'
                    CALL Ereport(  RoutineName, ErrorStatus, Cmessage )

                  END IF

                END IF   ! Hybrid levels
              END IF   ! Height Gen methods
            END IF    ! N_LEVELS > 1

        !--------------------------------------------------------
        ! Set pseudo level number if variable is on pseudo levels
        !--------------------------------------------------------
            IF (STM_Record % pt_code > 0) THEN
              Lookup(lbplev,k) = k - k_out + 1
            END IF

          END DO ! K=K_OUT...
          k_out=k_out+(n_levels*n_plevels)
          recondat_node => recondat_node % next
        END DO ! While associated
      END DO ! All Sections
    END DO ! All internal/sub models
!-------------------------------------------------------------------
! Initialise LOOKUP fields from PPXREF
!-------------------------------------------------------------------

    DO k=1, Hdr_Out % Len2Lookup
      item=MOD(Lookup(item_code,k),1000)
      section=(Lookup(item_code,k)-item)/1000
      model=Lookup(model_code,k)
      STM_Record = Rcf_Exppx( model, section, item)

      IF ( Hdr_Out % fixhd( FH_HorizGrid ) < 100) THEN
        Lookup( lbcode,k)= 1
        Lookup( lbhem, k)= Hdr_Out % FixHd( FH_HorizGrid )
      ELSE
        Lookup( lbcode,k)= 101 !100 added for non-standard polar axis
        Lookup( lbhem, k)= Hdr_Out % FixHd( FH_HorizGrid ) - 100
      END IF

  ! LBCs have an lbhem value of 99
      IF (STM_Record % grid_type == ppx_atm_lbc_theta .OR.                     &
          STM_Record % grid_type == ppx_atm_lbc_u      .OR.                    &
          STM_Record % grid_type == ppx_atm_lbc_v) THEN
        Lookup( lbhem, k)= 99
      END IF

      Lookup( lbext,k )  = 0 ! No extra data
      Lookup( lbrel,k )  = 3 ! Header release number currently 3
      Lookup( lbfc,k )   = STM_Record % field_code
      Lookup( lbvc,k )   = STM_Record % lbvc_code
      Lookup( lbegin,k ) = 0
      Lookup( lbnrec,k ) = 0
      Lookup( lbproj,k ) = iproj
      Lookup( lbtyp,k )  = STM_Record % cf_fieldcode

      IF (Lookup( lblev,k )  ==  0) THEN
        Lookup( lblev,k ) = STM_Record % cf_levelcode
      END IF

  ! DEPENDS ON: get_um_version_id
      Lookup( lbsrce,k )=get_um_version_id(model_id)

      IF (Lookup( data_type,k )  ==  0) THEN
        Lookup( data_type,k ) = STM_Record % data_type
      END IF

      IF (Lookup( lbpack,k )  ==  0) THEN
        Lookup( lbpack,k ) = STM_Record % dump_packing
        IF (dump_pack==2 .OR. dump_pack==3 ) THEN
      ! Do not pack data ; Override packing indicator from PPXREF
          n1 = 0   !   No packing
          Lookup( lbpack,k ) = (Lookup( lbpack,k )/10)*10 + n1
        END IF
      END IF

      Lookup( bmdi,k ) = TRANSFER( rmdi, int_val )
      Lookup( bmks,k ) = TRANSFER( 1.0,  int_val )
    END DO

!-------------------------------------------------------------------
! Change LOOKUP to allow for change in horizontal dimensions
!-------------------------------------------------------------------

    DO k=1, Hdr_Out % Len2Lookup

      item=MOD(Lookup( item_code,k ),1000)
      section=(Lookup( item_code,k )-item)/1000
      model=Lookup( model_code,k )
      STM_Record = Rcf_Exppx( model, section, item)

      IF (AREA_Expand(k) == 1) THEN

    ! Get N2 and N3 from whole value of LBPACK
        whole=Lookup( lbpack,k )
        n2=MOD(INT(whole/10),10)
        n3=MOD(INT(whole/100),10)

    ! Find grid type to calculate grid info
        ppxref_grid_type = STM_Record % grid_type

        IF (n2 == 2 .AND. n3 == 1) THEN
          Lookup( lbrow,k ) = 0
          Lookup( lbnpt,k ) = 0
        ELSE IF (n2 == 1 .AND. n3 == 1) THEN
          Lookup( lbrow,k ) = 0
          Lookup( lbnpt,k ) = 0
        ELSE
      ! Only deal with C-grid
          SELECT CASE( ppxref_grid_type )
          CASE (ppx_atm_cuall)
            Lookup( lbrow,k ) = Output_Grid % Glob_u_rows
            Lookup( lbnpt,k ) = Output_Grid % Glob_u_row_length

          CASE (ppx_atm_cvall)
            Lookup( lbrow,k ) = Output_Grid % Glob_v_rows
            Lookup( lbnpt,k ) = Output_Grid % Glob_v_row_length

          CASE (ppx_atm_lbc_theta) ! LBC theta points
            Lookup( lbrow,k ) = rimwidtha(Rima_type_norm)
            Lookup( lbnpt,k ) = Output_Grid % Glob_p_row_length

          CASE (ppx_atm_lbc_u) ! LBC U points
            Lookup( lbrow,k ) = rimwidtha(Rima_type_norm)
            Lookup( lbnpt,k ) = Output_Grid % Glob_u_row_length - 1
                            ! Note this is actual rather than
                            ! oversized and only exists for LAM

          CASE (ppx_atm_lbc_v) ! LBC V points
            Lookup( lbrow,k ) = rimwidtha(Rima_type_norm)
            Lookup( lbnpt,k ) = Output_Grid % Glob_v_row_length

          CASE (ppx_atm_river)
            Lookup( lbrow,k ) = Output_Grid % Glob_r_rows
            Lookup( lbnpt,k ) = Output_Grid % Glob_r_row_length

          CASE (ppx_atm_ozone)
            Lookup( lbrow,k ) = Output_Grid % Glob_p_rows
            IF (ZonAvOzone) THEN
              Lookup( lbnpt,k ) = 1
            ELSE
              Lookup( lbnpt,k ) = Output_Grid % Glob_p_row_length
            END IF
          CASE Default
            Lookup( lbrow,k ) = Output_Grid % Glob_p_rows
            Lookup( lbnpt,k ) = Output_Grid % Glob_p_row_length
          END SELECT
        END IF

        Lookup( bplat,k ) = TRANSFER( Hdr_Out % RealC( RC_PoleLat  ),          &
            int_val)
        Lookup( bplon,k ) = TRANSFER( Hdr_Out % RealC( RC_PoleLong ),          &
            int_val)
        Lookup( bgor, k ) = TRANSFER( 0.                            ,          &
            int_val)

        Lookup( bdy,k ) = TRANSFER( Hdr_Out % RealC(2), int_val )
        IF (Lookup( lbnpt,k ) ==  1) THEN
          Lookup( bdx,k ) = TRANSFER( 360., int_val )
        ELSE
          Lookup( bdx,k ) = TRANSFER( Hdr_Out % RealC( RC_LongSpacing),        &
              int_val )
        END IF

! First set grid offset
        IF (Hdr_Out % Fixhd(FH_GridStagger) == FH_Gridstagger_Endgame) THEN
          grid_offset = 0.5
        ELSE
          grid_offset = 0.0
        END IF

! Assume P grid
        Lookup( bzy,k ) = TRANSFER( Hdr_Out % RealC(3) -                       &
            (1.-grid_offset)*Hdr_Out % RealC(2), int_val )

        Lookup( bzx,k ) = TRANSFER( Hdr_Out % RealC(4) -                       &
            (1.-grid_offset)*Hdr_Out % RealC(1), int_val)

! Deal with special grids.
        IF (ppxref_grid_type == ppx_atm_cuall) THEN

          Lookup( bzy,k ) = TRANSFER( Hdr_Out % RealC(3) -                     &
              (1.0-grid_offset)*Hdr_Out % RealC(2),int_val )
          Lookup( bzx,k ) = TRANSFER( Hdr_Out % RealC(4) -                     &
              (0.5+grid_offset)*Hdr_Out % RealC(1),int_val )
        ELSE IF (ppxref_grid_type == ppx_atm_cvall) THEN

          Lookup( bzy,k ) = TRANSFER( Hdr_Out % RealC(3) -                     &
              (0.5+grid_offset)*Hdr_Out % RealC(2),int_val )
          Lookup( bzx,k ) = TRANSFER( Hdr_Out % RealC(4) -                     &
              (1.0-grid_offset)*Hdr_Out % RealC(1),int_val )
        ELSE IF (ppxref_grid_type == ppx_atm_river) THEN
          Lookup( bdy,k ) = TRANSFER( 180.0/Output_Grid % Glob_r_rows, int_val)
          riv_bzy = Hdr_Out % RealC(3)                                         &
              - 0.5*( 180.0/Output_Grid % Glob_r_rows )
          Lookup( bzy,k ) = TRANSFER( riv_bzy,int_val)
          Lookup( bdx,k ) =                                                    &
              TRANSFER( 360.0/Output_Grid % Glob_r_row_length,                 &
              int_val)
          riv_bzx = Hdr_Out % RealC(4)                                         &
              - 0.5*(360.0/Output_Grid % Glob_r_row_length)
          Lookup( bzx,k ) = TRANSFER( riv_bzx, int_val )
        END IF
      END IF

  ! Set lookup 13 from LCAL360.
      IF (lcal360) THEN
        Lookup( lbtim,k ) = 2
      ELSE
        Lookup( lbtim,k ) = 1
      END IF

    END DO

! Set BZY and BZX to RMDI in variable resolution header
    IF (Hdr_Out %RealC(1)< 0) THEN
      DO k=1, Hdr_Out % Len2Lookup
        Lookup( bzy,k ) = TRANSFER( rmdi, int_val )
        Lookup( bzx,k ) = TRANSFER( rmdi, int_val )
      END DO
    END IF

    RETURN
  END SUBROUTINE Rcf_Setup_Lookup
END MODULE Rcf_Setup_Lookup_Mod

