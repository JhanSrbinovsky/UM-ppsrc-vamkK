! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Deletes duplicate diags & times; checks for overlap levs & times.
!
! Subroutine Interface:

      SUBROUTINE DUPLIC(NRECS,NTIMES,NLEVELS)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE version_mod, ONLY:                                            &
          nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e
      USE Submodel_Mod
      USE stextend_mod, ONLY: LLISTTY, INDX_S, LIST_S, ITIM_S,  &
                              RLEVLST_S, LEVLST_S

      IMPLICIT NONE
!
! Description:
!   Deletes duplicate diagnostic entries from STASH list; deletes
!   duplicate STASH times;
!   checks for overlap of levels and times, to some extent.
!   Called by STPROC.
!   Input : NRECS   No. of STASH list records
!           NTIMES  No. of STASH times
!           NLEVELS No. of STASH levels
!           LIST_S  STASH list array with prelim. pointer system
!           ITIM_S  STASH times array
!   Output: NRECS   Reduced no. of STASH list records
!           NTIMES  Reduced no. of STASH times
!           NLEVEL  Reduced no. of STASH levels
!           ITIM_S  Reduced STASH times array
!           LIST_S  Reduced STASH list with prelim. pointers,
!                           consistent with STASH times.
!
! Method:
!
!   (a) STASH times tables in ITIM_S.
!       The times at which STASH processing occurs for a diagnostic
!   IREC may be specified by the entries (start_time,end_time,period)
!   in LIST_S.
!       Alternatively, if LIST_S(st_freq_code,IREC) has value '-n',
!   then STASH processing times for this diagnostic are given by a
!   'times table' in ITIM_S.
!   (In such a case, the above 3 entries in LIST_S are ignored).
!   The times table is given by column 'n' of ITIM_S, i.e.,
!   ITIM_S(time,n).
!   In this routine, the logical array entry LTUSE(n) is set to
!   .TRUE. if col 'n' of ITIM_S contains a times table. Any column
!   of ITIM_S which does not contain a times table is filled with
!   -1's. The cols which contain times tables are then shuffled along,
!   so that they occupy the first NTIMES cols of ITIM_S. The pointers
!   in LIST_S(st_freq_code,IREC) are altered accordingly.
!
!   (b) STASH levels lists in (R)LEVLST_S.
!       The levels on which STASH processing occurs for a diagnostic
!   IREC is specified by the entries (output_bottom, output_top) in
!   LIST_S.
!     If LIST_S(bot)=m, then output is on a range of model levels,
!   with level m as the bottom level, and LIST_S(top) points to the
!   top output model level.
!     If LIST_S(bot)=-n, then there is a levels list in col 'n' of
!   LEVLST_S, and LIST_S(top) contains a code value indicating the
!   type of levels (model, pressures, heights or theta). Each levels
!   list also has a corresponding entry in LLISTTY, indicating whether
!   the list is real or integer.
!     In this routine, the cols of LEVLST_S which contain levels lists
!   are shuffled along so that they occupy the first NLEVELS cols of
!   LEVLST_S. The pointers in LIST_S(output_bottom,IREC), and the
!   entries in LLISTTY, are altered accordingly.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
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

! Subroutine Arguments:
!
!   Scalar arguments with intent(InOut):

      INTEGER NRECS      ! No. of STASH list records
      INTEGER NTIMES     ! No. of STASH times
      INTEGER NLEVELS    ! No. of STASH levels

! Local scalars:

      LOGICAL LTRPT
      LOGICAL TESTFLAG
      INTEGER I
      INTEGER I1
      INTEGER I2
      INTEGER IEND
      INTEGER IITM
      INTEGER IL
      INTEGER IREC
      INTEGER ISEC
      INTEGER ISTR
      INTEGER IT
      INTEGER ITAGS1
      INTEGER ITAGS2
      INTEGER ITAGU1
      INTEGER ITAGU2
      INTEGER ITEND1
      INTEGER ITEND2
      INTEGER MODL
      INTEGER NLEVSW
      INTEGER NRECSW
      INTEGER NTIMESW

! Local arrays:

      LOGICAL LTUSE(2*NPROFTP+2)  ! LTUSE(n) set to .T. if column n
                                  ! in ITIM_S contains a STASH times
                                  ! table.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of Header -----------------------------------------------------


! Initialise LTUSE array

      IF (lhook) CALL dr_hook('DUPLIC',zhook_in,zhook_handle)
      DO I=1,NTIMES
        LTUSE(I)=.FALSE.
      END DO


! Blank out unused STASH times

      DO IREC=1,NRECS
        IF (LIST_S(st_freq_code,IREC) <  0) THEN   ! STASH times table
          LTUSE(-LIST_S(st_freq_code,IREC))=.TRUE. !  exists for IREC
        END IF
      END DO

      DO I=1,NTIMES
        IF(.NOT.LTUSE(I)) ITIM_S(1,I)=-1  ! Fill unused columns in
      END DO                              ! ITIM_S with -1 in each row.


! Delete blank STASH times

      NTIMESW=1

      DO IT=1,NTIMES

! If col 'IT' contains a times table, find
! corresponding record IREC in LIST_S, and replace entry '-IT'
! by '-NTIMESW'. In each case, NTIMESW <= IT.

        IF(ITIM_S(1,IT) /= -1) THEN
          DO IREC=1,NRECS

            IF (LIST_S(st_freq_code,IREC) == -IT) THEN
                LIST_S(st_freq_code,IREC)=-NTIMESW
            END IF

          END DO

          IF (IT /= NTIMESW) THEN
 ! Move times table in col 'IT' to col 'NTIMESW'. Hence array
 ! ITIM_S is compressed.
            DO I=1,NTIMEP
              ITIM_S(I,NTIMESW)=ITIM_S(I,IT)
            END DO
          END IF

          NTIMESW=NTIMESW+1
        END IF
      END DO

      NTIMES=NTIMESW-1  ! No. of STASH-times tables remaining, so far


! Delete blank STASH levels

      NLEVSW=1

      DO IL=1,NLEVELS

! If col 'IL' of LEVLST_S contains a levs list, then find corresponding
!  record IREC in LIST_S and replace entry '-IL' by '-NLEVSW'. In each
! case, NLEVSW <= IL.

        IF(LEVLST_S(1,IL) /= 0) THEN
          DO IREC=1,NRECS
            IF (LIST_S(st_output_bottom,IREC) == -IL) THEN
                LIST_S(st_output_bottom,IREC)=-NLEVSW
            END IF
          END DO
          IF(IL /= NLEVSW) THEN
 ! Move levels list in col 'IL' to col 'NLEVSW'. Hence array
 ! LEVLST_S is compressed.
            DO I=1,NLEVP_S
              LEVLST_S(I,NLEVSW)=LEVLST_S(I,IL)
              RLEVLST_S(I,NLEVSW)=RLEVLST_S(I,IL)
            END DO
            LLISTTY(NLEVSW)=LLISTTY(IL) ! Move corresponding entry in
          END IF                        ! LLISTTY

          NLEVSW=NLEVSW+1
        END IF
      END DO

      NLEVELS=NLEVSW-1


! Check for duplication/overlap of STASH levels

      NRECSW=NRECS

      DO MODL  = 1,N_INTERNAL_MODEL_MAX
      DO ISEC  = 0,NSECTP
      DO IITM  = 1,NITEMP

        IF(INDX_S(2,MODL,ISEC,IITM) >= 2) THEN  !More than one STASH rec
                                                !  for (model,sec,item)
          ISTR=     INDX_S(1,MODL,ISEC,IITM)    !1st record with m,s,i
          IEND=ISTR+INDX_S(2,MODL,ISEC,IITM)-1  !Last record with m,s,i

          DO I1=ISTR,IEND-1

           ITAGS1=LIST_S(st_macrotag,I1)/1000              ! System tag
           ITAGU1=LIST_S(st_macrotag,I1)-1000*ITAGS1       ! User tag

           IF (LIST_S(st_model_code,I1) <= N_INTERNAL_MODEL_MAX) THEN
!              Not flagged redundant
            DO I2=I1+1,IEND

              ITAGS2=LIST_S(st_macrotag,I2)/1000           ! System tag
              ITAGU2=LIST_S(st_macrotag,I2)-1000*ITAGS2    ! User tag

              IF((LIST_S(st_proc_no_code,I1) ==                         &
     &            LIST_S(st_proc_no_code,I2)).AND.                      &
     &           (LIST_S(st_freq_code,I1) ==                            &
     &            LIST_S(st_freq_code,I2)).AND.                         &
     &           (LIST_S(st_period_code,I1) ==                          &
     &            LIST_S(st_period_code,I2)).AND.                       &
     &           (LIST_S(st_gridpoint_code,I1) ==                       &
     &            LIST_S(st_gridpoint_code,I2)).AND.                    &
     &           (LIST_S(st_weight_code,I1) ==                          &
     &            LIST_S(st_weight_code,I2)).AND.                       &
     &           (LIST_S(st_north_code,I1) ==                           &
     &            LIST_S(st_north_code,I2)).AND.                        &
     &           (LIST_S(st_south_code,I1) ==                           &
     &            LIST_S(st_south_code,I2)).AND.                        &
     &           (LIST_S(st_west_code,I1) ==                            &
     &            LIST_S(st_west_code,I2)).AND.                         &
     &           (LIST_S(st_east_code,I1) ==                            &
     &            LIST_S(st_east_code,I2)).AND.                         &
     &           (LIST_S(st_input_code,I1) ==                           &
     &            LIST_S(st_input_code,I2)).AND.                        &
     &           (LIST_S(st_output_code,I1) ==                          &
     &            LIST_S(st_output_code,I2)).AND.                       &
     &           (LIST_S(st_series_ptr,I1) ==                           &
     &            LIST_S(st_series_ptr,I2)).AND.                        &
     &           (LIST_S(st_pseudo_out,I1) ==                           &
     &            LIST_S(st_pseudo_out,I2)).AND.                        &
     &        ((ITAGS1 == ITAGS2).OR.(ITAGS1 == 0).OR.                  &
     &        (ITAGS2 == 0)).AND.                                       &
     &        ((ITAGU1 == ITAGU2).OR.(ITAGU1 == 0).OR.                  &
     &        (ITAGU2 == 0)).AND.                                       &
     &      (LIST_S(st_model_code,I2) <= N_INTERNAL_MODEL_MAX)) THEN
!            Not flagged redundant

! If they are the same in all but time and level

                ITEND1=LIST_S(st_end_time_code,I1)
                ITEND2=LIST_S(st_end_time_code,I2)

                IF(ITEND1 == -1) ITEND1=                                &
     &             LIST_S(st_start_time_code,I2)+1     ! Force overlap

                IF(ITEND2 == -1) ITEND2=                                &
     &             LIST_S(st_start_time_code,I1)+1     ! Force overlap

! Where period_code is zero we have to prevent second part from
! being evaluated so break this OR statement into two parts:
                TESTFLAG=.FALSE.
                IF((LIST_S(st_period_code,I1) == 0).OR.                 &
     &            (LIST_S(st_period_code,I1) == -1)) THEN
                  TESTFLAG=.TRUE.
                ELSEIF((MOD(LIST_S(st_start_time_code,I2)-              &
     &               LIST_S(st_start_time_code,I1),                     &
     &               LIST_S(st_period_code,I1)) == 0)) THEN
                  TESTFLAG=.TRUE.
                ENDIF

                IF((.NOT.((LIST_S(st_start_time_code,I1)                &
     &                                            >  ITEND2).OR.        &
     &            (ITEND1 <  LIST_S(st_start_time_code,I2))).OR.        &
     &              (LIST_S(st_output_code,I1) >  0)).AND.              &
     &          (MOD(LIST_S(st_start_time_code,I2)-                     &
     &               LIST_S(st_start_time_code,I1),                     &
     &               LIST_S(st_freq_code,I1)) == 0).AND.                &
     &               TESTFLAG.AND.                                      &
     &              (LIST_S(st_output_bottom,I2) ==                     &
     &               LIST_S(st_output_bottom,I1)).AND.                  &
     &              (LIST_S(st_output_top,I2) ==                        &
     &               LIST_S(st_output_top,I1))) THEN

! (Times overlap or in dump) and overlay in freq & period
!                            and levels the same

                  IF(ITAGU1 == 0) THEN
                    LIST_S(st_macrotag,I1)=ITAGU2
                    ITAGU1=ITAGU2
                  ELSE
                    LIST_S(st_macrotag,I1)=ITAGU1
                  END IF

                  IF(ITAGS1 == 0) THEN
                    LIST_S(st_macrotag,I1)=                             &
     &              ITAGS2*1000+LIST_S(st_macrotag,I1)
                    ITAGS1=ITAGS2
                  ELSE
                    LIST_S(st_macrotag,I1)=                             &
     &              ITAGS1*1000+LIST_S(st_macrotag,I1)
                  END IF

                  LIST_S(st_start_time_code,I1)=                        &
     &            MIN(LIST_S(st_start_time_code,I1),                    &
     &                LIST_S(st_start_time_code,I2))

                  IF((LIST_S(st_end_time_code,I1) == -1).OR.            &
     &               (LIST_S(st_end_time_code,I2) == -1)) THEN
                      LIST_S(st_end_time_code,I1)=-1
                  ELSE
                    LIST_S(st_end_time_code,I1)=                        &
     &              MAX(LIST_S(st_end_time_code,I1),                    &
     &                  LIST_S(st_end_time_code,I2))
                  END IF

                  LIST_S(st_model_code,I2)=N_INTERNAL_MODEL_MAX+1
! Sets model id to be greater than no of models,
! so that this diag is put at the end of any sorted list.

                  NRECSW=NRECSW-1

                  DO I=ISTR,IEND                      !Change pointers
                    IF(LIST_S(st_input_code,I )   ==                    &
     &                -LIST_S(NELEMP+1     ,I2)) THEN
                       LIST_S(st_input_code,I ) =                       &
     &                -LIST_S(NELEMP+1     ,I1)
                    END IF
                  END DO
                END IF

              END IF   ! I1,I2 comparison
            END DO     ! I2
           END IF      ! I1 Not flagged redundant
          END DO       ! I1

        END IF   ! More than one STASH record for m,s,i
      END DO     ! Items
      END DO     ! Sections
      END DO     ! Models

! Remove unwanted records (i.e., those flagged redundant)

! DEPENDS ON: order
      CALL ORDER(NRECS)
      NRECS=NRECSW
! DEPENDS ON: sindx
      CALL SINDX(NRECS)
!
      IF (lhook) CALL dr_hook('DUPLIC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DUPLIC
