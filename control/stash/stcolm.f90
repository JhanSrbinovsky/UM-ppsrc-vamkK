! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine: STCOLM ---------------------------------------------------
!
!  Purpose: Calculate weighted column mean within a region specified
!           by a lower left hand and upper right hand corner.
!           (STASH service routine).
!
!
!  Tested under compiler:   cft77
!  Tested under OS version: UNICOS 5.1
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D712
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!  ---------------------------------------------------------------
!
!  Interface and arguments: ------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: STASH
      SUBROUTINE STCOLM(fieldin,vx,vy,vz,fld_type,halo_type,            &
     &                  lwrap,lmasswt,                                  &
     &                  xstart,ystart,xend,yend,                        &
     &                  fieldout,index_lev,zsize,                       &
     &                  pstar_weight,                                   &
     &                  area_weight,mask,                               &
     &                  level_code,mask_code,weight_code,rmdi,          &
     &                  icode,cmessage)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE sterr_mod, ONLY: st_upper_less_lower, st_not_supported,       &
                           st_no_data,st_nd, st_bad_array_param,        &
                           st_bad_address, st_unknown,                  &
                           st_bad_wraparound, st_illegal_weight,        &
                           unknown_weight, unknown_mask,                &
                           unknown_processing, nonsense
      IMPLICIT NONE
!
      INTEGER                                                           &
     &    vx,vy,vz,                                                     &
                                                ! IN  input field size
     &    fld_type,                                                     &
                                                ! IN  field type(u/v/p)
     &    halo_type,                                                    &
                                                ! IN  halo type
     &    xstart,ystart,                                                &
                                                ! IN  lower LH corner
     &    xend,yend,                                                    &
                                                ! IN  upper RH corner
     &    zsize,                                                        &
                                      ! IN no of horiz levels to process
     &    index_lev(zsize),                                             &
                                      ! IN offset for each horiz level
     &    level_code,                                                   &
                                                ! IN  input level code
     &    mask_code,                                                    &
                                                ! IN  masking code
     &    weight_code,                                                  &
                                                ! IN  weighting code
     &    icode                                 ! OUT error return code
      CHARACTER(LEN=*)                                                     &
     &    cmessage                              ! OUT error return msg
      LOGICAL                                                           &
     &    lwrap,                                                        &
                                                ! IN  TRUE if wraparound
     &    lmasswt,                                                      &
                                                ! IN  TRUE if masswts OK
     &    mask(vx+1,vy)                         ! IN  mask array
      REAL                                                              &
     &    fieldin(vx,vy,vz),                                            &
                                                ! IN  input field
     &    fieldout(xstart:xend,ystart:yend),                            &
                                                ! OUT output field
     &    pstar_weight(vx+1,vy,zsize),                                  &
                                                ! IN  mass weight factor
     &    area_weight(vx+1,vy),                                         &
                                                ! IN  area weight factor
     &    rmdi                                  ! IN  missing data indic
!*----------------------------------------------------------------------
!
! External subroutines called
!
!
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end
!
! Local variables
!
        INTEGER i,j,k,ii,kk       ! ARRAY INDICES FOR VARIABLE

        REAL SUMCTOP(xstart:xend,ystart:yend)
        REAL SUMCBOT(xstart:xend,ystart:yend)

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

!L----------------------------------------------------------------------
!L 0. Initialise sums
!L
      IF (lhook) CALL dr_hook('STCOLM',zhook_in,zhook_handle)
      DO i=xstart,xend
        DO j=ystart,yend
          SUMCTOP(i,j)=0.0
          SUMCBOT(i,j)=0.0
        END DO
      END DO
!L----------------------------------------------------------------------
!L 1. Form column sums
!L
!L 1.1 NULL weighting or area weighting
!L
      IF (weight_code == stash_weight_null_code.or.                     &
     &    weight_code == stash_weight_area_code) THEN
        DO kk=1,zsize
          k=index_lev(kk)
          DO i=xstart,xend

            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
            ELSE
              ii=i
            ENDIF

            DO j=ystart,yend
              SUMCBOT(i,j)=SUMCBOT(i,j)+1.0
              SUMCTOP(i,j)=SUMCTOP(i,j)+fieldin(ii,j,k)
            END DO
          END DO
        END DO
!L
!L 1.2 mass weighting
!L
      ELSEIF (weight_code == stash_weight_mass_code) THEN
        IF (.NOT.lmasswt) THEN
! Mass-weighting on level types with no mass-weight defined is not
! supported - should be prevented by UI
          cmessage='STCOLM  : mass-weights not defined for this diag'
          icode=st_illegal_weight
          goto 999
        ELSE
          DO kk=1,zsize
            k=index_lev(kk)
            DO i=xstart,xend
            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
              ELSE
                ii=i
              END IF
              DO j=ystart,yend
                SUMCBOT(i,j)=SUMCBOT(i,j) +  pstar_weight(ii,j,kk)
                SUMCTOP(i,j)=SUMCTOP(i,j) +  fieldin(ii,j,k)*           &
     &                                       pstar_weight(ii,j,kk)
              END DO
            END DO
          END DO
        END IF
      ELSE
        cmessage='STCOLM  : Invalid weighting code detected'
        icode=unknown_weight
        goto 999
      ENDIF
!L----------------------------------------------------------------------
!L 2. Perform masking (set missing data at masked points) - compute mean
!L
      DO i=xstart,xend
            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
        ELSE
          ii=i
        ENDIF
        DO j=ystart,yend
          IF (mask(ii,j)) THEN
            fieldout(i,j)=SUMCTOP(i,j)/SUMCBOT(i,j)
          ELSE
            fieldout(i,j)=rmdi
          ENDIF
        ENDDO
      ENDDO
!L
  999 CONTINUE
      IF (lhook) CALL dr_hook('STCOLM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE STCOLM
