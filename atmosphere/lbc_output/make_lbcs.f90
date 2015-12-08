! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generate the LBC data from the model prognostics
!
! Subroutine Interface:

      SUBROUTINE make_lbcs(                                             &
! Prognostic
        SOURCE_FIELD,                                                   &
        LOCAL_ROW_LENGTH,                                               &
        LOCAL_ROWS,                                                     &
        GLOBAL_ROW_LENGTH,                                              &
        GLOBAL_ROWS,                                                    &
        SOURCE_HALO_X,                                                  &
        SOURCE_HALO_Y,                                                  &
        SOURCE_LEVELS,                                                  &
        SOURCE_FLD_TYPE,                                                &
        SOURCE_HALO_TYPE,                                               &
        lambda_source,                                                  &
        phi_source,                                                     &
        source_pole_lat,                                                &
        source_pole_long,                                               &
        source_cyclic,                                                  &
        source_rotated,                                                 &
! Orography Field
        OROG_FIELD,                                                     &
        OROG_LOCAL_ROW_LENGTH,                                          &
        OROG_LOCAL_ROWS,                                                &
        OROG_GLOBAL_ROW_LENGTH,                                         &
        OROG_GLOBAL_ROWS,                                               &
        OROG_FLD_TYPE,                                                  &
        OROG_HALO_TYPE,                                                 &
        OROG_HALO_X,                                                    &
        OROG_HALO_Y,                                                    &
        lambda_source_orog,                                             &
        phi_source_orog,                                                &
! LBC field
        lbc_row_len,                                                    &
        lbc_rows,                                                       &
        lbc_levels,                                                     &
        lbc_delta_lat,                                                  &
        lbc_delta_long,                                                 &
        lbc_first_lat,                                                  &
        lbc_first_long,                                                 &
        lbc_pole_lat,                                                   &
        lbc_pole_long,                                                  &
        lbc_halo_x,                                                     &
        lbc_halo_y,                                                     &
        lbc_halo_x_coords,                                              &
        lbc_halo_y_coords,                                              &
        rimwidth,                                                       &
        LBC,                                                            &
        coeff1,                                                         &
        coeff2,                                                         &
        Lbc_size,                                                       &
        level_type,                                                     &
! VarRes grid info in degrees
        same_rotation,                                                  &
        l_eg_grid,                                                      &
        L_extra_surface_level_dump,                                     &
        Lambda_in,                                                      &
        Phi_in,                                                         &
! hi - indexes and weights
        l_calc_lbc_wts,                                                 &
        lbc_index_bl  , lbc_index_br  ,                                 &
        lbc_weights_tr, lbc_weights_br,                                 &
        lbc_weights_bl, lbc_weights_tl,                                 &
! vi
        max_seg_size,                                                   &
        N_SEGS,                                                         &
        MAX_LEVELS_PER_PE,                                              &
        GATHER_PE                                                       &
      , L_VI                                                            &
! vi - src
      ,       src_model_levels                                          &
      ,       src_ht_gen_method                                         &
      ,       src_first_r_rho                                           &
      ,       src_z_top_model                                           &
      ,       src_eta_theta                                             &
      ,       src_eta_rho                                               &
      ,       src_first_level                                           &
      ,       src_last_level                                            &
! vi - lbc
      ,       v_interp_order                                            &
      ,       lbc_model_levels                                          &
      ,       lbc_ht_gen_method                                         &
      ,       lbc_first_r_rho                                           &
      ,       lbc_z_top_model                                           &
      ,       lbc_eta_theta                                             &
      ,       lbc_eta_rho                                               &
      ,       lbc_first_level                                           &
      ,       lbc_last_level                                            &
      ,       i_uv                                                      &
       )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE gcom_mod
      USE h_int_bl_mod, ONLY: h_int_bl



      IMPLICIT NONE

!
! Description:
!   Generates the Atmosphere LBCs from the model prognostics. Performs
!   horizontal and vertical interpolation as required.
!
! Method:
!   1. Horizontal Interpolation. (HI)
!      HI from source grid to lbc grid is done by level on a PE.
!      a) Determine no of levels to be interpolated on each PE.
!      b) Gather global source data on PEs from local data.
!      c) Determine the horizontal interpolation coefficients.
!      d) DO the HI through suboutine H_INT_BL.
!      e) If no vertical interpolation, gather all lbc data on PE 0.
!
!   2. Vertical interpolation (VI) (to follow)
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

! Subroutine arguments

      INTEGER                                                           &
        LOCAL_ROW_LENGTH                                                &
                      ! IN  : East-West size of source field
      , LOCAL_ROWS                                                      &
                      ! IN  : North-South size of source field
      , GLOBAL_ROW_LENGTH                                               &
                      ! IN  : East-West Size of full gathered level
                      !     : glsize(1,fld_type)
      , GLOBAL_ROWS                                                     &
                      ! IN  : North-South Size of full gathered level
                      !     : glsize(2,fld_type)
      , SOURCE_HALO_X                                                   &
                      ! IN  : Size of East-West halo for source field
      , SOURCE_HALO_Y                                                   &
                      ! IN  : Size of North-South halo for source field

      , SOURCE_LEVELS                                                   &
                      ! IN  : Number of levels for source field
      , SOURCE_FLD_TYPE                                                 &
                        ! IN  : Field type of source field
      , SOURCE_HALO_TYPE ! IN  : Halo type of source field

      Logical                                                           &
        SOURCE_CYCLIC                                                   &
                      ! IN   : T if Source Grid is cyclic
      , SOURCE_ROTATED! IN   : T if Source Grid is rotated

       Logical                                                          &
        L_CALC_LBC_WTS ! IN  : T if Interp coeffs to be calculated

      INTEGER                                                           &
        OROG_LOCAL_ROW_LENGTH                                           &
                      ! IN  : East-West size of source field
      , OROG_LOCAL_ROWS                                                 &
                      ! IN  : North-South size of source field
      , OROG_GLOBAL_ROW_LENGTH                                          &
                      ! IN  : East-West Size of full gathered level
                      !     : glsize(1,fld_type)
      , OROG_GLOBAL_ROWS                                                &
                      ! IN  : North-South Size of full gathered level
                      !     : glsize(2,fld_type)
      , OROG_HALO_X                                                     &
                    ! IN  : Size of East-West halo for source field
      , OROG_HALO_Y                                                     &
                    ! IN  : Size of North-South halo for source field
      , OROG_HALO_TYPE                                                  &
      , OROG_FLD_TYPE

      REAL ::                                                           &
        source_pole_lat                                                 &
      , source_pole_long                                                &
      , lbc_delta_lat                                                   &
      , lbc_delta_long                                                  &
      , lbc_first_lat                                                   &
      , lbc_first_long                                                  &
      , lbc_pole_lat                                                    &
      , lbc_pole_long

! If target grid on same pole as input grid?
      Logical same_rotation
      LOGICAL :: l_eg_grid
      LOGICAL :: L_extra_surface_level_dump
! Input grid info in degrees      
      REAL                                                              &
        Lambda_in( 1-lbc_halo_x_coords: lbc_row_len+lbc_halo_x_coords ) & 
      , Phi_in( 1-lbc_halo_y_coords: lbc_rows+lbc_halo_y_coords )
      
      REAL :: lambda_source(global_row_length)
      REAL :: phi_source(global_rows)
      REAL :: lambda_source_orog(orog_global_row_length)
      REAL :: phi_source_orog(orog_global_rows)


      Integer                                                           &
        rimwidth                                                        &
      , lbc_halo_x                                                      &
      , lbc_halo_y                                                      &
      , lbc_halo_x_coords                                               &
      , lbc_halo_y_coords                                               &
      , LBC_SIZE                                                        &
                      ! IN  : Size of single level of LBC array
      , LBC_LEVELS                                                      &
                      ! IN  : Number of levels of LBC data
      , level_type                                                      &
      , max_seg_size                                                    &
                      ! IN  : Size of data on each processor for vert.
                      !       int. step. Something like:
                      ! MAX(minimum_segment_size,((LBC_SIZE-1)/nproc)+1)
      , N_SEGS                                                          &
                      ! IN  : Number of segments,
                      ! given by ((LBC_SIZE-1)/LBC_SEG_SIZE)+1
      , MAX_LEVELS_PER_PE                                               &
                      ! IN  : Maximum no of levels per PE - basiCALLy:
                      !       MAX_LEVELS_PER_PE=(SOURCE_LEVELS-1)/NPROC
      , GATHER_PE                                                       &
                      ! IN  : Processor that will contain the final
                      !       3D LBC field
                      !       (probably usually PE 0)
      , i_uv          ! IN  : 1 => u field, 2 => v field, 0 => other

      REAL                                                              &
       Source_Field(1-SOURCE_HALO_X:LOCAL_ROW_LENGTH+SOURCE_HALO_X,     &
                    1-SOURCE_HALO_Y:LOCAL_ROWS+SOURCE_HALO_Y,           &
                    SOURCE_LEVELS)

      REAL                                                              &
       Orog_Field(1-OROG_HALO_X:OROG_LOCAL_ROW_LENGTH+OROG_HALO_X,      &
                  1-OROG_HALO_Y:OROG_LOCAL_ROWS+OROG_HALO_Y )

      REAL                                                              &
        Coeff1 (lbc_size)                                               &
                            ! \ Coefficients for rotating winds
      , Coeff2 (lbc_size)   ! / on lbc points


      REAL lbc (lbc_size, lbc_levels)  !  OUT : LBC data

! vi
      LOGICAL                                                           &
       l_vi           ! IN  : Perform vertical interpolation?

! vi - src
      INTEGER   :: src_model_levels
      INTEGER   :: src_first_level
      INTEGER   :: src_last_level
      INTEGER   :: src_ht_gen_method
      INTEGER   :: src_first_r_rho
      REAL      :: src_z_top_model
      REAL      :: src_eta_theta (0:src_model_levels)
      REAL      :: src_eta_rho   (  src_model_levels)

! vi - lbc
      INTEGER   :: v_interp_order
      INTEGER   :: lbc_model_levels
      INTEGER   :: lbc_first_level
      INTEGER   :: lbc_last_level
      INTEGER   :: lbc_ht_gen_method
      INTEGER   :: lbc_first_r_rho
      REAL      :: lbc_z_top_model
      REAL      :: lbc_eta_theta (0:lbc_model_levels)
      REAL      :: lbc_eta_rho   (  lbc_model_levels)

! Local parameters:

      Character (Len=*), Parameter :: RoutineName= 'Make_LBCs'

! Local scalars:

! vi
      INTEGER                                                           &
        send_map(7,N_SEGS*MAX_LEVELS_PER_PE)                            &
      , recv_map(7,N_SEGS*SOURCE_LEVELS)

      INTEGER                                                           &
        map(SOURCE_LEVELS)                                              &
                           ! mapping of processors holding which level
      , LEVEL_ON_GAT_PE(SOURCE_LEVELS)                                  &
                                       ! what is the level number on the
                           !            gathering PE of this level
      , LEVEL_COUNT(0:nproc-1)  ! Count of levels on each PE

      INTEGER :: i, j   ! loop counters

      INTEGER                                                           &
        iproc                                                           &
                    ! loop counter over processors
      , iseg                                                            &
                    ! loop counter over segments
      , seg_start                                                       &
                    ! first point of segment
      , seg_end                                                         &
                    ! last point of segment
      , seg_size_pe                                                     &
                    ! local segment size
      , seg_size(n_segs)                                                &
      , n_send                                                          &
                    ! number of items of data to send
      , n_recv                                                          &
                    ! number of items of data to receive
      , k                                                               &
                    ! loop counter over levels
      , offset                                                          &
                    ! may need to skip extra surface level
                    ! in EG fields if going EG -> ND.
      , info                                                            &
                    ! return code
      , flag                                                            &
                    ! GCOM input code
      , ErrorStatus ! Return Code

      Character (Len=80)  :: CMessage

! Local dynamic arrays:

! Whole levels of the source field gathered onto processors

      Real, dimension (:,:,:), allocatable :: Gathered_Source_Field
      Real, dimension (  :,:), allocatable :: Gathered_Orog_Field

! LBC data after horizontal interpolation

      Real, dimension (  :,:), allocatable :: lbc_data
      Real, dimension (    :), allocatable :: lbc_orog

! Indexes for bottom left & right points in bi-linear interpolation




      Integer, dimension (lbc_size) :: lbc_index_bl
      Integer, dimension (lbc_size) :: lbc_index_br

! Weights for 4 points in bi-linear interpolation

      Real,    dimension (lbc_size) :: lbc_weights_tr
      Real,    dimension (lbc_size) :: lbc_weights_br
      Real,    dimension (lbc_size) :: lbc_weights_tl
      Real,    dimension (lbc_size) :: lbc_weights_bl

! LBC data before/after vertical interpolation

      Real,    dimension (:,:), allocatable :: lbc_vi_data_in
      Real,    dimension (:,:), allocatable :: lbc_vi_data_out
      Real,    dimension (  :), allocatable :: lbc_vi_orog

! Temp

      Integer kk,ipt  ! temp for write statement.
      integer halo_x, halo_y, lbc_row_len, lbc_rows
      integer row,pt,level
      logical lprint

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header
!--------------------------------------------------------------------

! 1.0 Set up arrays to map which levels are on which processors

      IF (lhook) CALL dr_hook('MAKE_LBCS',zhook_in,zhook_handle)
      Do iproc=0,nproc-1
        level_count(iproc)=0
      Enddo

      Do k=1,source_levels
        map(k)              = MOD(k-1,nproc)
        level_count(map(k)) = level_count(map(k))+1
        level_on_gat_pe(k)  = level_count(map(k))
      Enddo

! 2.0 Gather levels to different PEs





      allocate ( Gathered_Source_Field ( global_row_length,             &
                                         global_rows,                   &
                                         max_levels_per_pe ) )

! DEPENDS ON: gather_field_ml
      CALL Gather_Field_ML (                                            &
        Source_Field, Gathered_Source_Field,                            &
        Local_Row_Length + 2*Source_Halo_x,                             &
        Local_Rows + 2*Source_Halo_y,                                   &
        Source_Levels,                                                  &
        Global_Row_Length, Global_Rows, Max_Levels_per_pe,              &
        Map, Source_Fld_Type, Source_Halo_Type)



      If ( l_vi ) Then





! Set up a global field of orography on PE 0

        ALLOCATE ( Gathered_Orog_Field ( orog_global_row_length,        &
                                         orog_global_rows ) )

! DEPENDS ON: gather_field
        CALL Gather_Field (Orog_Field,                                  &
                           Gathered_Orog_Field,                         &
                           Orog_Local_Row_Length, Orog_Local_Rows,      &
                           Orog_Global_Row_Length, Orog_Global_Rows,    &
                           Orog_Fld_Type, Orog_Halo_Type,               &
                           0, gc_all_proc_group, info, cmessage )

! 3.0 Perform horizontal interpolation.
!     This processor has LEVEL_COUNT(mype) levels to do

! 3.1 Allocate workspace for the interpolation coefficients

!      work space ALLOCATEd in GEN_INTF_A

! 3.2 Allocate work space to horizontally interpolated orography

        ALLOCATE ( lbc_orog (lbc_size) )

        IF (mype == 0) THEN

! 3.3 Calculate the interpolation coefficients for orography field

! NB. LBC_Interp_Coeffs must be calculated for the orography before
!     the prognostic so that on exit COEFF1 and COEFF2 correspond to
!     the prognostic.

          IF (l_calc_lbc_wts ) THEN

! DEPENDS ON: lbc_interp_coeffs
            CALL LBC_interp_coeffs (                                    &
               lbc_size                                                 &
      ,        orog_global_row_length                                   &
      ,        orog_global_rows                                         &
      ,        lambda_source_orog                                       &
      ,        phi_source_orog                                          &
      ,        source_pole_lat                                          &
      ,        source_pole_long                                         &
      ,        source_cyclic                                            &
      ,        source_rotated                                           &
      ,        lbc_row_len                                              &
      ,        lbc_rows                                                 &
      ,        same_rotation                                            &
      ,        l_eg_grid                                                &
      ,        Lambda_in                                                &
      ,        Phi_in                                                   &
      ,        lbc_delta_lat                                            &
      ,        lbc_delta_long                                           &
      ,        lbc_first_lat                                            &
      ,        lbc_first_long                                           &
      ,        lbc_pole_lat                                             &
      ,        lbc_pole_long                                            &
      ,        rimwidth                                                 &
      ,        lbc_halo_x                                               &
      ,        lbc_halo_y                                               &
      ,        lbc_halo_x_coords                                        &
      ,        lbc_halo_y_coords                                        &
      ,        lbc_index_bl                                             &
      ,        lbc_index_br                                             &
      ,        lbc_weights_tr                                           &
      ,        lbc_weights_br                                           &
      ,        lbc_weights_bl                                           &
      ,        lbc_weights_tl                                           &
      ,        coeff1                                                   &
      ,        coeff2                                                   &
      ,        i_uv                                                     &
       )

          End If  !  l_calc_lbc_wts

! 3.4 Do the bi-linear horizontal interpolation for orography

          CALL h_int_bl (                                               &
               orog_global_rows                                         &
      ,        orog_global_row_length                                   &
      ,        lbc_size                                                 &
      ,        lbc_index_bl                                             &
      ,        lbc_index_br                                             &
      ,        gathered_orog_field                                      &
      ,        lbc_weights_bl                                           &
      ,        lbc_weights_br                                           &
      ,        lbc_weights_tl                                           &
      ,        lbc_weights_tr                                           &
      ,        lbc_orog                                                 &
       )

        End If  !  If mype=0
! Need to deallocate field outside mype block since all processors
! allocate it. Check on allocated status as it is not allocated on 
! the MAKEBC preprocessor define.
        IF (allocated ( gathered_orog_field)) &
            deallocate ( gathered_orog_field)

      END IF
      
! 3.2 Calculate the interpolation coefficients for the prognostic

      IF ( l_calc_lbc_wts) THEN

! DEPENDS ON: lbc_interp_coeffs
      CALL LBC_interp_coeffs (                                          &
           lbc_size                                                     &
      ,    global_row_length                                            &
      ,    global_rows                                                  &
      ,    lambda_source                                                &
      ,    phi_source                                                   &
      ,    source_pole_lat                                              &
      ,    source_pole_long                                             &
      ,    source_cyclic                                                &
      ,    source_rotated                                               &
      ,    lbc_row_len                                                  &
      ,    lbc_rows                                                     &
      ,    same_rotation                                                &
      ,    l_eg_grid                                                    &
      ,    Lambda_in                                                    &
      ,    Phi_in                                                       &
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
      ,    coeff1                                                       &
      ,    coeff2                                                       &
      ,    i_uv                                                         &
       )

      END IF   !  l_calc_lbc_wts

! 3.3 Allocate work space for the horizontally interpolated data

      ALLOCATE ( lbc_data (lbc_size, max_levels_per_pe) )

! 3.4 DO the bi-linear horizontal interpolation

      Do level = 1, level_count(mype)  !  Levels to do on this pe.
        CALL h_int_bl (                                                 &
             global_rows                                                &
                                   !  glsize(2)
      ,      global_row_length                                          &
                                   !  glsize(1)
      ,      lbc_size                                                   &
      ,      lbc_index_bl                                               &
      ,      lbc_index_br                                               &
      ,      gathered_source_field (1,1,level)                          &
      ,      lbc_weights_bl                                             &
      ,      lbc_weights_br                                             &
      ,      lbc_weights_tl                                             &
      ,      lbc_weights_tr                                             &
      ,      lbc_data (1,level)                                         &
       )

      END DO  !  Loop over levels

! 3.5 DEALLOCATE space reserved for interpolation coefficients

!      work space deallocated in GEN_INTF_A

!     Deallocate space reserved for gathered source data

      deallocate ( gathered_source_field )

! We should now have an array LBC_DATA(:,:) where we have filled in
! 1:LEVEL_COUNT(mype) levels of data (this may be zero if there are
! more processors than levels !)

! 4.0 Perform Vertical Interpolation.

      IF (L_VI) THEN

!       Determine segment size for each segment.
!       Max_Seg_Size is maximum size of any segment.

        DO iseg = 1, n_segs
          seg_start=((iseg-1)*Max_Seg_Size)+1
          seg_end=MIN(LBC_SIZE,seg_start+Max_Seg_Size-1)
          seg_size(iseg) = seg_end - seg_start + 1
        END DO

!       Get segment size for this pe.
        IF (mype < n_segs) THEN
          seg_size_pe = seg_size(mype+1)
        Else
          seg_size_pe = 1
        END IF

        ! Move data so that each processor has a segment of data, on
        ! all levels

        n_send=0    ! number of messages I send
        n_recv=0    ! number of messages I receive
        send_map (:,:) = 0
        recv_map (:,:) = 0

        DO iseg=1,N_SEGS

          iproc=iseg-1  ! This is the processor doing vertical
!                       ! interpolation on segment iseg

          seg_start=((iseg-1)*Max_Seg_Size)+1
          seg_end=MIN(LBC_SIZE,seg_start+Max_Seg_Size-1)

          DO k=1,SOURCE_LEVELS

            IF (MAP(k) == mype) THEN  ! I've got data to send
              n_send=n_send+1
              send_map(S_DESTINATION_PE,n_send)=iproc
              send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)=            &
                seg_start+((LEVEL_ON_GAT_PE(k)-1)*LBC_SIZE)
              send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=1
              send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=1
              send_map(S_ELEMENT_LENGTH,n_send)=seg_end-seg_start+1
              send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=            &
                1+(k-1)*seg_size(iseg)
              send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=1
            END IF

            IF (iproc == mype) THEN  ! I am receiving data
              n_recv=n_recv+1
              recv_map(R_SOURCE_PE,n_recv)=MAP(k)
              recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=            &
                1+(k-1)*seg_size(iseg)

              recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=1
              recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=1
              recv_map(R_ELEMENT_LENGTH,n_recv)=seg_end-seg_start+1
              recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=            &
                seg_start+((LEVEL_ON_GAT_PE(k)-1)*LBC_SIZE)
              recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=1
            END IF

          END DO ! k
        END DO ! iseg

        ALLOCATE ( lbc_vi_data_in (seg_size_pe, source_levels) )

        flag=GC_NONE
        info=GC_NONE

        CALL GCG_RALLTOALLE(LBC_DATA,send_map,n_send,                   &
                            LBC_SIZE*MAX_LEVELS_PER_PE,                 &
                            LBC_VI_DATA_IN,recv_map,n_recv,             &
                            SEG_SIZE_PE*SOURCE_LEVELS,                  &
                            gc_all_proc_group,flag,info)

        IF (info /= GC_NONE) THEN
          ErrorStatus = 10
          WRITE (cmessage,*) 'GCG_RALLTOALLE failed'

          CALL Ereport ( RoutineName, ErrorStatus, cmessage)
        END IF

        ! Now the prognostic is in LBC_VI_DATA_IN(:,SOURCE_LEVELS)

! We now need to scatter the orography to the same number of segments

        n_send=0    ! number of messages I send
        n_recv=0    ! number of messages I receive
        send_map (:,:) = 0
        recv_map (:,:) = 0

        DO iseg=1,N_SEGS

          iproc=iseg-1  ! This is the processor doing vertical
!                       ! interpolation on segment iseg

          seg_start=((iseg-1)*Max_Seg_Size)+1
          seg_end=MIN(LBC_SIZE,seg_start+Max_Seg_Size-1)

            IF ( mype == 0 ) THEN  ! Only PE 0 has orog to send out
              n_send=n_send+1
              send_map(S_DESTINATION_PE,n_send)=iproc
              send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)= seg_start
              send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=1
              send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=1
              send_map(S_ELEMENT_LENGTH,n_send)=seg_end-seg_start+1
              send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=1
              send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=1
            END IF

            IF (iproc == mype) THEN  ! This PE to receive data
              n_recv=n_recv+1
              recv_map(R_SOURCE_PE,n_recv)=0
              recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=1
              recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=1
              recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=1
              recv_map(R_ELEMENT_LENGTH,n_recv)=seg_end-seg_start+1
              recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=seg_start
              recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=1
            END IF

        END DO ! iseg

        ALLOCATE ( lbc_vi_orog (seg_size_pe) )

        flag=GC_NONE
        info=GC_NONE

        CALL GCG_RALLTOALLE(LBC_OROG,send_map,n_send,                   &
                            LBC_SIZE,                                   &
                            LBC_VI_OROG,recv_map,n_recv,                &
                            seg_size_pe,                                &
                            gc_all_proc_group,flag,info)

        IF (info /= GC_NONE) THEN
          ErrorStatus = 10
          WRITE (cmessage,*) 'GCG_RALLTOALLE failed'

          CALL Ereport ( RoutineName, ErrorStatus, cmessage)
        END IF

        ! Now the orography is in LBC_VI_OROG(:)

        DEALLOCATE (lbc_orog)
        ALLOCATE ( lbc_vi_data_out (seg_size_pe, lbc_levels) )

        IF (mype < N_SEGS) THEN

          ! I have to do some vertical interpolation...

! DEPENDS ON: lbc_vert_interp
          CALL LBC_Vert_Interp (                                        &
               LBC_VI_Data_in,                                          &
               LBC_VI_Orog,                                             &
               LBC_VI_Data_out,                                         &
               seg_size_pe,                                             &
               level_type                                               &
! src
      ,        src_model_levels                                         &
      ,        source_levels                                            &
      ,        src_first_level                                          &
      ,        src_last_level                                           &
      ,        src_ht_gen_method                                        &
      ,        src_first_r_rho                                          &
      ,        src_z_top_model                                          &
      ,        src_eta_theta                                            &
      ,        src_eta_rho                                              &
! lbc
      ,        v_interp_order                                           &
      ,        lbc_model_levels                                         &
      ,        lbc_levels                                               &
      ,        lbc_first_level                                          &
      ,        lbc_last_level                                           &
      ,        lbc_ht_gen_method                                        &
      ,        lbc_first_r_rho                                          &
      ,        lbc_z_top_model                                          &
      ,        lbc_eta_theta                                            &
      ,        lbc_eta_rho                                              &
       )

        END IF

        ! Now the data is in LBC_VI_DATA_OUT(:,LBC_LEVELS)
        ! This needs to be moved to LBC on processor GATHER_PE

        n_send=0
        n_recv=0
        send_map (:,:) = 0
        recv_map (:,:) = 0

        DO iseg=1,N_SEGS

          iproc=iseg-1  ! This is the processor doing vertical
!                       ! interpolation on segment iseg

          seg_start=((iseg-1)*Max_Seg_Size)+1
          seg_end=MIN(LBC_SIZE,seg_start+Max_Seg_Size-1)

          IF (mype == iproc) THEN ! I've got data to send
            n_send=n_send+1
            send_map(S_DESTINATION_PE,n_send)=GATHER_PE
            send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)=1
            send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=LBC_LEVELS
            send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=SEG_SIZE(iseg)
            send_map(S_ELEMENT_LENGTH,n_send)=seg_end-seg_start+1
            send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=seg_start
            send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=LBC_SIZE
          END IF

          IF (mype == GATHER_PE) THEN ! I've got data to receive
            n_recv=n_recv+1
            recv_map(R_SOURCE_PE,n_recv)=iproc
            recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=seg_start
            recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=LBC_LEVELS
            recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=LBC_SIZE
            recv_map(R_ELEMENT_LENGTH,n_recv)=seg_end-seg_start+1
            recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=1
            recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=SEG_SIZE(iseg)
          END IF

        END DO ! iseg

        flag=GC_NONE
        info=GC_NONE

        CALL GCG_RALLTOALLE(LBC_VI_DATA_OUT,send_map,n_send,            &
                            SEG_SIZE_PE*LBC_LEVELS,                     &
                            LBC,recv_map,n_recv,                        &
                            LBC_SIZE*LBC_LEVELS,                        &
                            gc_all_proc_group,flag,info)

        IF (info /= GC_NONE) THEN
          ErrorStatus = 20
          WRITE (cmessage,*) 'GCG_RALLTOALLE failed'

          CALL Ereport ( RoutineName, ErrorStatus, cmessage)
        END IF

! DEALLOCATE arrays
        DEALLOCATE ( lbc_vi_data_in  )
        DEALLOCATE ( lbc_vi_data_out )
        DEALLOCATE ( lbc_vi_orog     )

      Else ! .NOT. L_VI - ie. no vertical interpolation

        ! Need to move the data from LBC_DATA with different levels
        ! on different PEs to LBC with all the data on GATHER_PE
        
        ! ENDGAME:
        !     Need to ignore extra surface level if no vertical 
        !     interpolation. The VI is switched off ONLY if the
        !     source and LAM have the same level set. Therefore
        !     source_levels will always equal lbc_levels (except
        !     for extra surface level for certain ENDGAME fields.)
        
        IF (L_extra_surface_level_dump) THEN
           offset = 2
        ELSE
           offset = 1
        END IF

        n_send=0
        n_recv=0
        send_map (:,:) = 0
        recv_map (:,:) = 0

        DO k=offset,source_levels

          IF (Map(k) == mype) THEN ! I need to send data
            n_send=n_send+1
            send_map(S_DESTINATION_PE,n_send)=GATHER_PE
            send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)=              &
              1+((LEVEL_ON_GAT_PE(k)-1)*LBC_SIZE)
            send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=1
            send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=1
            send_map(S_ELEMENT_LENGTH,n_send)=LBC_SIZE
            send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=              &
              1+(k-offset)*LBC_SIZE
            send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=1
          END IF

          IF (mype == GATHER_PE) THEN ! I need to receive data
            n_recv=n_recv+1
            recv_map(R_SOURCE_PE,n_recv)=MAP(k)
            recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=              &
              1+(k-offset)*LBC_SIZE
            recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=1
            recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=1
            recv_map(R_ELEMENT_LENGTH,n_recv)=LBC_SIZE
            recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=              &
              1+((LEVEL_ON_GAT_PE(k)-1)*LBC_SIZE)
            recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=1
          Endif

        Enddo ! k

        flag=GC_NONE
        info=GC_NONE

        Call GCG_RALLTOALLE(LBC_DATA,send_map,n_send,                   &
                            LBC_SIZE*MAX_LEVELS_PER_PE,                 &
                            LBC,recv_map,n_recv,                        &
                            LBC_SIZE*LBC_LEVELS,                        &
                            gc_all_proc_group,flag,info)

        If (info /= GC_NONE) THEN
          ErrorStatus = 30
          Write (cmessage,*) 'GCG_RALLTOALLE failed'

          Call Ereport ( RoutineName, ErrorStatus, cmessage)
        Endif

      Endif ! IF (L_VI)

      deallocate ( lbc_data )

      ! The array LBC on processor GATHER_PE should now contain the full
      ! 3D LBC array

      IF (lhook) CALL dr_hook('MAKE_LBCS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE MAKE_LBCS
