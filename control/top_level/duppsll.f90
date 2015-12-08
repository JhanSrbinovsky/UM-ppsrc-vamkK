! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Determine whether a pseudo lev list is a duplicate of another one
! Subroutine Interface:

      SUBROUTINE DUPPSLL(LDUPLL,NDUPLL)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE version_mod, ONLY:                                            &
          nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e
      USE Submodel_Mod
      USE stextend_mod, ONLY: LENPLST

      IMPLICIT NONE
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

! Subroutine arguments:

!   Scalar arguments with intent(out):

      LOGICAL LDUPLL
      LOGICAL LLOCAL
      INTEGER NDUPLL

! Local scalars:

      INTEGER I
      INTEGER J

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of Header -------------------------------------------------


      IF (lhook) CALL dr_hook('DUPPSLL',zhook_in,zhook_handle)
      LDUPLL=.FALSE.
      NDUPLL=0
      DO I=1,NPSLISTS-1
        IF(LENPLST(I) == LENPLST(NPSLISTS)) THEN
          LLOCAL=.TRUE.
          DO J=1,LENPLST(NPSLISTS)
            IF(PSLIST_D(J,I) /= PSLIST_D(J,NPSLISTS)) THEN
              LLOCAL=.FALSE.
              EXIT
            END IF
          END DO
          
          IF(LLOCAL) THEN
            LDUPLL=.TRUE.
            NDUPLL=I
            GOTO 9999
          END IF
        END IF
      END DO

 9999 CONTINUE
      IF (lhook) CALL dr_hook('DUPPSLL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DUPPSLL
