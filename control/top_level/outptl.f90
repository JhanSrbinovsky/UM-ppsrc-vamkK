! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calc stash list output lens; reset boundary spec for full area output.

! Subroutine Interface:

      SUBROUTINE OUTPTL(                                                &
                        NRECS,ErrorStatus,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE Decomp_DB

      USE cppxref_mod, ONLY:                                            &
          ppx_grid_type, ppx_halo_type, ppx_pt_code, ppx_lv_code
      USE version_mod, ONLY:                                            &
          nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e

      USE Submodel_Mod
      USE stextend_mod, ONLY: NRECS_TS, LIST_S, LEVLST_S, LENPLST

      IMPLICIT NONE

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 90 
!    Written to UMDP3 programming standards version 8.3.
!
! Global variables:
! COMDECK CSTASH
! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.
!
!
! Declarations:
! Imported global variables:
!    None, but see note above.

! Global parameters:

! Global scalars:
      INTEGER      NDIAG   ! No. of diagnostics
      INTEGER      NTPROF  ! No. of time profiles
      INTEGER      NSERIES ! No. of stash time series
      INTEGER      NDPROF  ! No. of domain profiles
      INTEGER      NUPROF  ! No. of useage profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in VERSION_MOD module
      INTEGER      MODL_B(NDIAGPM)  ! Internal model no.
      INTEGER      ISEC_B(NDIAGPM)  ! Section
      INTEGER      ITEM_B(NDIAGPM)  ! Item
      INTEGER      ITIM_B(NDIAGPM)  ! Time profile number
      INTEGER      IDOM_B(NDIAGPM)  ! Domain profile number
      INTEGER      IUSE_B(NDIAGPM)  ! Useage profile number

!   Time profile information:

      CHARACTER(LEN=8)  TIMPRO(NPROFTP)         ! Name of profile
      INTEGER      ITYP_T(NPROFTP)         ! Type of profile
      INTEGER      INTV_T(NPROFTP)         ! Time Interval
      CHARACTER(LEN=2)  UNT1_T(NPROFTP)         ! Units for time interval
      INTEGER      ISAM_T(NPROFTP)         ! Sampling period
      CHARACTER(LEN=2)  UNT2_T(NPROFTP)         ! Units for sampling period
      INTEGER      IOPT_T(NPROFTP)         ! Output option
      INTEGER      ISTR_T(NPROFTP)         ! Output Start time
      INTEGER      IEND_T(NPROFTP)         ! Output End time
      INTEGER      ISDT_T(6, NPROFTP)      ! Output Start date
      INTEGER      IEDT_T(6, NPROFTP)      ! Output End date
      INTEGER      IFRE_T(NPROFTP)         ! Output frequency
      INTEGER      IOFF_T(NPROFTP)         ! Offset for sampling
      CHARACTER(LEN=2)  UNT3_T(NPROFTP)         ! Units for output times
      INTEGER      ITIM_T(NPROFTP)         ! No. of times in times table
      INTEGER      ISER_T(NTIMEP ,NPROFTP) ! Times table (with units)
      INTEGER      MODL_T(NPROFTP)         ! Indicates internal model
                                           !  for each times table

!   Domain profile information:

      CHARACTER(LEN=8) DOMPRO  (NPROFDP)           ! Name of domain profile
      INTEGER     IOPL_D  (NPROFDP)           ! Levels option
      INTEGER     LEVB_D  (NPROFDP)           ! Bottom level
      INTEGER     LEVT_D  (NPROFDP)           ! Top level
      INTEGER     IOPA_D  (NPROFDP)           ! Area option
      INTEGER     INTH_D  (NPROFDP)           ! North boundary
      INTEGER     ISTH_D  (NPROFDP)           ! South boundary
      INTEGER     IEST_D  (NPROFDP)           ! East boundary
      INTEGER     IWST_D  (NPROFDP)           ! West boundary
      INTEGER     IMSK_D  (NPROFDP)           ! Mask type
      INTEGER     IMN_D   (NPROFDP)           ! Meaning option
      INTEGER     IWT_D   (NPROFDP)           ! Weighting option
      CHARACTER(LEN=1) TS_D    (NPROFDP)           ! Time series profile
      INTEGER     IG_TS
      INTEGER     I1_TS
      INTEGER     I51_TS
      INTEGER     BLIM_TS (NTimSerP)
      INTEGER     TLIM_TS (NTimSerP)
      REAL        BLIMR_TS(NTimSerP)
      REAL        TLIMR_TS(NTimSerP)
      INTEGER     NLIM_TS (NTimSerP)
      INTEGER     SLIM_TS (NTimSerP)
      INTEGER     ELIM_TS (NTimSerP)
      INTEGER     WLIM_TS (NTimSerP)
      INTEGER     ILEV_D  (NPROFDP)           ! Output levels code
      INTEGER     LEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      REAL       RLEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      INTEGER     PLT_D   (NPROFDP)
      INTEGER     PLLEN_D (NPROFDP)
      INTEGER     PLPOS_D (NPROFDP)
      INTEGER     PSLIST_D(NPSLEVP ,NPSLISTP)
      INTEGER     NPSLISTS
      EQUIVALENCE        (RLEVLST_D,LEVLST_D)

! Useage information:

      CHARACTER(LEN=8) USEPRO(NPROFUP)   ! Name of useage profile
      INTEGER     LOCN_U(NPROFUP)   ! Storage location of profile
      INTEGER     IUNT_U(NPROFUP)   ! Unit no.

! Information from ppxref file:

      INTEGER      MODEL_ST       ! Internal model number
      INTEGER      ISPACE         ! Space code
      INTEGER      ITIMA          ! Time availability code
      INTEGER      IGP            ! Grid of data code
      INTEGER      ILEV           ! Level type code
      INTEGER      IBOT           ! First level code
      INTEGER      ITOP           ! Last level code
      INTEGER      IFLAG          ! Level compression flag
      INTEGER      IOPN(6)        ! Sectional option code
      INTEGER      VMSK           ! Integer equiv of bin vers mask
      INTEGER      IPSEUDO        ! Pseudo dimension type
      INTEGER      IPFIRST        ! First pseudo dim code
      INTEGER      IPLAST         ! Last pseudo dim code
      INTEGER      PTR_PROG       ! Section zero point back
      INTEGER      HALO_TYPE      ! Type of halo the field has

! PP output file units
      INTEGER      PPlen2LkUp(OUTFILE_S:OUTFILE_E)
      CHARACTER(LEN=1)  FTOutUnit (OUTFILE_S:OUTFILE_E)

! COMMON blocks:
      COMMON/STCHA/ TIMPRO,UNT1_T,UNT2_T,UNT3_T,DOMPRO,TS_D,            &
     &  USEPRO,FTOutUnit

      COMMON/STSH/                                                      &
     &  NDIAG   ,MODL_B  ,ISEC_B ,ITEM_B  ,ITIM_B  ,IDOM_B  ,IUSE_B,    &
     &  NTPROF  ,ITYP_T  ,INTV_T ,ISAM_T  ,ITIM_T  ,                    &
     &  IOPT_T  ,ISTR_T  ,IEND_T ,IFRE_T  ,IOFF_T, ISER_T  ,MODL_T  ,   &
     &  NDPROF  ,IOPL_D  ,LEVB_D ,ISDT_T  ,IEDT_T  ,                    &
     &  IOPA_D  ,INTH_D  ,ISTH_D ,IEST_D  ,IWST_D  ,                    &
     &  IMSK_D  ,IMN_D   ,IWT_D  ,                                      &
     &  LEVT_D  ,LEVLST_D,                                              &
     &  PLT_D   ,PLLEN_D ,PLPOS_D,PSLIST_D,NPSLISTS,                    &
     &  BLIM_TS ,TLIM_TS ,BLIMR_TS,TLIMR_TS,IG_TS   ,I1_TS   ,          &
     &  NLIM_TS ,SLIM_TS ,ELIM_TS ,WLIM_TS ,I51_TS  ,NSERIES ,          &
     &  NUPROF  ,LOCN_U  ,IUNT_U ,                                      &
     &  MODEL_ST,ISPACE  ,ITIMA  ,IGP     ,                             &
     &  ILEV    ,IBOT    ,ITOP   ,IFLAG   ,IOPN    ,VMSK    ,           &
     &  IPSEUDO ,IPFIRST ,IPLAST ,PTR_PROG, HALO_TYPE,                  &
     & PPlen2LkUp

! CSTASH end
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

! Subroutine arguments
!   Array arguments with intent(in):
      INTEGER NRECS

!   Array arguments with intent(out):
      CHARACTER(LEN=80) CMESSAGE

! ErrorStatus:
      INTEGER ErrorStatus

! Local variables
      INTEGER output_length
      INTEGER IE
      INTEGER IN
      INTEGER IP_DIM
      INTEGER IREC
      INTEGER IS
      INTEGER MODL
      INTEGER ISEC
      INTEGER ITEM
      INTEGER IT_DIM
      INTEGER IW
      INTEGER IX_DIM
      INTEGER IY_DIM
      INTEGER IZ_DIM
      INTEGER                                                           &
! local versions of the global subdomain boundaries
        local_north,local_east,local_south,local_west                   &
      , local_IN,local_IE,local_IS,local_IW                             &
! global versions of the X and Y horizontal dimensions, and
! total output size
      , global_IX_DIM,global_IY_DIM,global_output_length                &
! variables indicating the decomposition type at various stages
      , orig_decomp,decomp_type

! Function and subroutine calls:
      INTEGER  EXPPXI

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      EXTERNAL EXPPXI

!- End of Header --------------------------------------------------

      IF (lhook) CALL dr_hook('OUTPTL',zhook_in,zhook_handle)
      orig_decomp=current_decomp_type

! Loop over STASH records
      DO IREC=1,NRECS

! Obtain model, section, item for this record
        MODL = LIST_S(st_model_code  ,IREC)
        ISEC = LIST_S(st_sect_no_code,IREC)
        ITEM = LIST_S(st_item_code   ,IREC)

! Set the correct decomposition type for this model
      IF (MODL  ==  ATMOS_IM) THEN
        decomp_type=decomp_standard_atmos
      ELSE
!       Shouldn't get to this
        decomp_type=decomp_unset
        WRITE(6,'(A)') 'OUTPTL : Error'
        WRITE(6,'(A,I4,A)') 'Unsupported Model ',MODL,' for MPP code'
        CMESSAGE='Unsupported Model for MPP code'
        ErrorStatus=-1
        GO TO 9999
      END IF

      IF (current_decomp_type  /=  decomp_type) THEN
        CALL CHANGE_DECOMPOSITION(decomp_type,ErrorStatus)
        IF (ErrorStatus  /=  0) THEN
          WRITE(6,'(A)') 'OUTPUTL : Error'
          WRITE(6,'(A,A,I4)') 'Call to CHANGE_DECOMPOSITION failed ',  &
                     'with decomposition type ',decomp_type
          CMESSAGE='Unsupported decomposition for MPP code'
          GO TO 9999
        END IF
      END IF


! Extract level code, grid type code from ppx lookup array
! DEPENDS ON: exppxi
        ILEV    = EXPPXI(MODL,ISEC,ITEM,ppx_lv_code,                    &
                                      ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IGP     = EXPPXI(MODL,ISEC,ITEM,ppx_grid_type,                  &
                                      ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPSEUDO = EXPPXI(MODL ,ISEC ,ITEM,ppx_pt_code      ,            &
                                      ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        HALO_TYPE = EXPPXI(MODL,ISEC,ITEM,ppx_halo_type,                &
                                      ErrorStatus,CMESSAGE)
        IF (LIST_S(st_proc_no_code,IREC) == 0) THEN
! Dummy record - output length zero
          LIST_S(st_output_length,IREC)=0
        ELSE IF(LIST_S(st_input_code,IREC) <  0.and.                    &
               LIST_S(st_proc_no_code,IREC) /= 8) THEN
! Child record - get output length from parent
          LIST_S(st_output_length,IREC)=                                &
          LIST_S(st_output_length,-LIST_S(st_input_code,IREC))

        ELSE
! Neither dummy nor child - calculate output length
!   T dimension (equals 1 except for the time series case)
          IF((LIST_S(st_proc_no_code,IREC) == 1).OR.                    &
             (LIST_S(st_proc_no_code,IREC) == 2).OR.                    &
             (LIST_S(st_proc_no_code,IREC) == 3).OR.                    &
             (LIST_S(st_proc_no_code,IREC) == 5).OR.                    &
             (LIST_S(st_proc_no_code,IREC) == 6))     THEN
            IT_DIM=1
          ELSE IF (LIST_S(st_proc_no_code,IREC) == 4.or.                &
                  LIST_S(st_proc_no_code,IREC) == 8) THEN
! Time series case
            IT_DIM=                                                     &
            LIST_S(st_period_code,IREC)/LIST_S(st_freq_code,IREC)

          ELSE
            WRITE(6,'(A,A,I4,A,I6)')'OUTPTL: ',  &
            'ERROR UNEXPECTED PROCESSING CODE',  &
            LIST_S(st_proc_no_code,IREC),'   FOR RECORD ',IREC
          END IF

          IF (LIST_S(st_series_ptr,IREC) == 0) THEN
! Set up local versions of the boundaries of the subdomain

! DEPENDS ON: global_to_local_subdomain
            CALL GLOBAL_TO_LOCAL_SUBDOMAIN( .TRUE., .TRUE.,             &
                                     IGP,HALO_TYPE,mype,                &
                                     LIST_S(st_south_code,IREC),        &
                                     LIST_S(st_east_code,IREC),         &
                                     LIST_S(st_north_code,IREC),        &
                                     LIST_S(st_west_code,IREC),         &
                                     local_south,local_east,            &
                                     local_north,local_west)

! Not a time series profile
!   X dimension
            IF (LIST_S(st_gridpoint_code,IREC) <  0) THEN
              WRITE(6,'(A,A,I4,A,I6)')'OUTPTL: ',  &
              'ERROR UNEXPECTED GRIDPOINT CODE',  &
              LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <  20) THEN
              IX_DIM=local_east-local_west+1
              global_IX_DIM=LIST_S(st_east_code,IREC)-                  &
                            LIST_S(st_west_code,IREC)+1
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <  30) THEN
              IX_DIM=1
              global_IX_DIM=1
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <  40) THEN
              IX_DIM=local_east-local_west+1
              global_IX_DIM=LIST_S(st_east_code,IREC)-                  &
                            LIST_S(st_west_code,IREC)+1
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <= 43) THEN
              IX_DIM=1
              global_IX_DIM=1
            ELSE
              WRITE(6,'(A,A,I4,A,I6)')'OUTPTL: ',  &
              'ERROR UNEXPECTED GRIDPOINT CODE',  &
              LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            END IF  ! X dim

            IF (IX_DIM <  1) THEN
! Area cut by global model
! DEPENDS ON: lltorc
              CALL LLTORC(IGP,90,-90,0,360,IN,IS,IW,IE)
! DEPENDS ON: global_to_local_subdomain
              CALL GLOBAL_TO_LOCAL_SUBDOMAIN(                           &
                .TRUE. , .TRUE. , IGP ,halo_type, mype ,                &
                IN,IE,IS,IW,                                            &
                local_IN,local_IE,local_IS,local_IW)
              IX_DIM=IX_DIM+local_IE-2*halosize(1,halo_type)
! Subtract two halos, because we don't want wrap around to include
! the halo at the end, and the beginning of field

            END IF
            IF (global_IX_DIM <  1) THEN
! DEPENDS ON: lltorc
              CALL LLTORC(IGP,90,-90,0,360,IN,IS,IW,IE)
              global_IX_DIM=global_IX_DIM+IE
            END IF

!   Y dimension
            IF (LIST_S(st_gridpoint_code,IREC) <  0) THEN
              WRITE(6,'(A,A,I4,A,I6)')'OUTPTL: ',  &
              'ERROR UNEXPECTED GRIDPOINT CODE',  &
              LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <  30) THEN
! Atmos grid - first lat is southern most
              IY_DIM=local_north-local_south+1
              global_IY_DIM=LIST_S(st_north_code,IREC)-                 &
                            LIST_S(st_south_code,IREC)+1
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <= 40) THEN
              IY_DIM=1
              global_IY_DIM=1
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <= 43) THEN
              IY_DIM=1
              global_IY_DIM=1
            ELSE
              WRITE(6,'(A,A,I4,A,I6)')'OUTPTL: ',  &
              'ERROR UNEXPECTED GRIDPOINT CODE',  &
              LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            END IF  ! Y dim

!   Z dimension
            IF (LIST_S(st_gridpoint_code,IREC) <  0) THEN
              WRITE(6,'(A,A,I4,A,I6)')'OUTPTL: ',  &
              'ERROR UNEXPECTED GRIDPOINT CODE',  &
              LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <  10) THEN
              IF (ILEV == 5) THEN
                IZ_DIM=1
              ELSE IF (LIST_S(st_output_bottom,IREC) <  0) THEN
                IZ_DIM=LEVLST_S(1,-LIST_S(st_output_bottom,IREC))
              ELSE
                IZ_DIM=LIST_S(st_output_top,IREC)-                      &
                       LIST_S(st_output_bottom,IREC)+1
              END IF
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <  20) THEN
              IZ_DIM=1
            ELSE IF (LIST_S(st_gridpoint_code,IREC) <= 43) THEN
              IF (ILEV == 5) THEN
                IZ_DIM=1
              ELSE IF (LIST_S(st_output_bottom,IREC) <  0) THEN
                IZ_DIM=LEVLST_S(1,-LIST_S(st_output_bottom,IREC))
              ELSE
                IZ_DIM=LIST_S(st_output_top,IREC)-                      &
                       LIST_S(st_output_bottom,IREC)+1
              END IF
            ELSE
              WRITE(6,'(A,A,I4,A,I6)')'OUTPTL: ',  &
              'ERROR UNEXPECTED GRIDPOINT CODE',  &
              LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            END IF  ! Z dim

!   P dimension - pseudo levels
            IF (IPSEUDO >  0) THEN
              IP_DIM=LENPLST(LIST_S(st_pseudo_out,IREC))
            ELSE
              IP_DIM=1
            END IF

! Output length - total number of points
            output_length = IT_DIM*IX_DIM*IY_DIM*IZ_DIM*IP_DIM

            global_output_length =                                      &
              IT_DIM*global_IX_DIM*global_IY_DIM*IZ_DIM*IP_DIM
            LIST_S(st_output_length,IREC) = output_length
            LIST_S(st_dump_output_length,IREC) =                        &
              global_output_length
            LIST_S(st_dump_level_output_length,IREC) =                  &
              global_IX_DIM*global_IY_DIM  ! size of horizontal field
          ELSE    ! Time series profile
            LIST_S(st_output_length,IREC)=                              &
             NRECS_TS(LIST_S(st_series_ptr,IREC))*IT_DIM+               &
            (NRECS_TS(LIST_S(st_series_ptr,IREC))+1)*6

            LIST_S(st_dump_output_length,IREC)=                         &
              LIST_S(st_output_length,IREC)
          END IF
        END IF    ! Neither dummy nor child
      END DO      ! Loop over STASH records
      IF ((orig_decomp  /=  current_decomp_type) .AND.                  &
          (orig_decomp  /=  decomp_unset)) THEN
        CALL CHANGE_DECOMPOSITION(orig_decomp,ErrorStatus)
      END IF

9999  CONTINUE
      IF (lhook) CALL dr_hook('OUTPTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE OUTPTL
