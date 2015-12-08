! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: get_lookup
!
! Purpose: Flux processing routine.
!          Return lookup table for one file
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine get_lookup (UnitIn, icode,                             &
! DUMP_AR2 array dimensions
     &  LEN_FIXHD, LEN_INTHD, LEN_REALHD,                               &
     &  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,         &
     &  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC,LEN2_FLDDEPC,          &
     &  LEN_EXTCNST,  LEN_DUMPHIST,                                     &
     &  LEN_CFI1, LEN_CFI2, LEN_CFI3,                                   &
     &  LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS,                               &
! DUMP_AR2 end
     & Len_data, LOOKUP)

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      implicit none

! declaration of argument list
      integer UnitIn  ! IN unit number of file to read
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected
! all arguments in DUMP_LEN are IN
!*L----------------- COMDECK DUMP_LEN --------------------------------

      INTEGER LEN_FIXHD
      INTEGER LEN_INTHD
      INTEGER LEN_REALHD
      INTEGER LEN1_LEVDEPC, LEN2_LEVDEPC
      INTEGER LEN1_ROWDEPC, LEN2_ROWDEPC
      INTEGER LEN1_COLDEPC, LEN2_COLDEPC
      INTEGER LEN1_FLDDEPC, LEN2_FLDDEPC
      INTEGER LEN_EXTCNST
      INTEGER LEN_DUMPHIST
      INTEGER LEN_CFI1, LEN_CFI2, LEN_CFI3
      INTEGER LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS
!*----------------------------------------------------------------------

      integer Len_data ! IN length of data in file

! LOOKUP intent is OUT and is declared by DUMP_DIM

! declaration of local arrays
! DUMP_DIM  declares local arrays for ancillary file header
!*L----------------- COMDECK DUMP_DIM ----------------------------------

      INTEGER FIXHD(LEN_FIXHD)
      INTEGER INTHD(LEN_INTHD)
      REAL    REALHD(LEN_REALHD)
      REAL    LEVDEPC(LEN1_LEVDEPC,LEN2_LEVDEPC)
      REAL    ROWDEPC(LEN1_ROWDEPC,LEN2_ROWDEPC)
      REAL    COLDEPC(LEN1_COLDEPC,LEN2_COLDEPC)
      REAL    FLDDEPC(LEN1_FLDDEPC,LEN2_FLDDEPC)
      REAL    EXTCNST(LEN_EXTCNST)
      REAL    DUMPHIST(LEN_DUMPHIST)
      INTEGER CFI1(LEN_CFI1), CFI2(LEN_CFI2), CFI3(LEN_CFI3)
      INTEGER LOOKUP(LEN1_LOOKUP_OBS,LEN2_LOOKUP_OBS)
!*----------------------------------------------------------------------

! declaration of local scalars
      integer START_BLOCK       !  start of data block
      CHARACTER(LEN=256) CMESSAGE    !  error message

      external READHEAD

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'get_lookup'  ! subroutine name for error messages

! 1. Read dump / ancillary file header

! DEPENDS ON: readhead
      call READHEAD(UnitIn,                                             &
!*L----------------- COMDECK DUMP_AR1--------------------------------

!      Array addresses and dimensions
     &  FIXHD,    LEN_FIXHD,                                            &
     &  INTHD,    LEN_INTHD,                                            &
     &  REALHD,   LEN_REALHD,                                           &
     &  LEVDEPC,  LEN1_LEVDEPC,LEN2_LEVDEPC,                            &
     &  ROWDEPC,  LEN1_ROWDEPC,LEN2_ROWDEPC,                            &
     &  COLDEPC,  LEN1_COLDEPC,LEN2_COLDEPC,                            &
     &  FLDDEPC,  LEN1_FLDDEPC,LEN2_FLDDEPC,                            &
     &  EXTCNST,  LEN_EXTCNST,                                          &
     &  DUMPHIST, LEN_DUMPHIST,                                         &
     &  CFI1,     LEN_CFI1,                                             &
     &  CFI2,     LEN_CFI2,                                             &
     &  CFI3,     LEN_CFI3,                                             &
     &  LOOKUP,   LEN1_LOOKUP_OBS,LEN2_LOOKUP_OBS,                      &
!*----------------------------------------------------------------------
     & Len_data,                                                        &
     & START_BLOCK,ICODE,CMESSAGE)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &  ' step 1. failed to read ancillary file header; cmessage is ',  &
     &  cmessage
        icode = 32
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE get_lookup
!----------------------------------------------------------------------
