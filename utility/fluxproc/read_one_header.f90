! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: read_one_header
!
! Purpose: Flux processing routine.
!          Opens and reads fixed header and lookup table of one file
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_one_header ( InUnit, ICode,                       &
     & Len_FixHd_P, Len1_Lookup_P, Len2_Lookup_P,                       &
     & Len2_Lookup_Actual, FixHd,                                       &
     & Lookup)
      use io
      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      USE chsunits_mod, ONLY : nunits
      implicit none

! declaration of argument list
      integer InUnit  ! IN     input unit number
      integer ICode   ! IN/OUT error code ; > 0 => fatal error detected

! dimensions used to declare arrays to be read in
      integer Len_FixHd_P   ! IN length of fixed header
      integer Len1_Lookup_P ! IN length of first dimension of Lookup
      integer Len2_Lookup_P ! IN max length of 2nd dimension of Lookup

! fixed header and lookup tables: intent OUT
      integer Len2_Lookup_Actual    ! OUT actual 2nd dimension of Lookup
      integer FixHd(Len_FixHd_P)                     ! OUT fixed header
      integer Lookup(Len1_Lookup_P, Len2_Lookup_P)   ! OUT lookup tables

!  declaration of globals
!----------------------------------------------------------------------
! comdeck: CENVIRON
! Purpose: defines environment variables for units
!          which have to be opened by open_file
!----------------------------------------------------------------------
! declaration of parameters

! declarations of common blocks
      common / LEnviron /    LEnv
      common / CEnviron /    CEnv

! declarations of variables
      integer LEnv(NUnits)      ! lengths of environment variable names
      CHARACTER(LEN=15) CEnv(NUnits) ! names of environment variables
!----------------------------------------------------------------------

! declaration of local arrays

! lengths of headers in fields file (local arrays)
! (this declares LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS)
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

! declaration of local scalars
      integer Len_data  ! length of data in file
      CHARACTER(LEN=256) CMessage ! error messages

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_one_header'  ! subroutine name for error messages
      CMessage = ' '

! 1 open  file
      call file_open(InUnit, CEnv(InUnit), LEnv(InUnit), 0, 0, icode)

      if (icode  >   0) then
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 1. failed to open file with environment name ',          &
     &  CEnv(InUnit)
        icode = 27
        go to 9999
      end if

! 2. Read fixed header
! DEPENDS ON: read_flh
      CALL READ_FLH(InUnit,FIXHD,LEN_FIXHD_P,ICODE,CMESSAGE)

      if ( icode  >   0) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2. unable to read fixed header; cmessage is ',      &
     &       cmessage
        icode = 28
        go to 9999
      end if

! 3. get dimensions from lookup table, check them and set actual
!    2nd dimension of lookup table

! 3.0  get dimensions from lookup table
      LEN_FIXHD = LEN_FIXHD_P
! DEPENDS ON: get_dim
      CALL GET_DIM(FIXHD,                                               &
! DUMP_AR2 array dimensions
     &  LEN_FIXHD, LEN_INTHD, LEN_REALHD,                               &
     &  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,         &
     &  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC,LEN2_FLDDEPC,          &
     &  LEN_EXTCNST,  LEN_DUMPHIST,                                     &
     &  LEN_CFI1, LEN_CFI2, LEN_CFI3,                                   &
     &  LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS,                               &
! DUMP_AR2 end
     & Len_data)

! 3.1 Set to zero any dimensions of headers which are less than zero
!     (readhead etc. fail if this is not done)

      if ( LEN1_COLDEPC  <   0) LEN1_COLDEPC = 0
      if ( LEN2_COLDEPC  <   0) LEN2_COLDEPC = 0
      if ( LEN1_FLDDEPC  <   0) LEN1_FLDDEPC = 0
      if ( LEN2_FLDDEPC  <   0) LEN2_FLDDEPC = 0
      if ( LEN_EXTCNST  <   0)  LEN_EXTCNST  = 0
      if ( LEN_DUMPHIST  <   0) LEN_DUMPHIST = 0
      if ( LEN_CFI1  <   0)     LEN_CFI1 = 0
      if ( LEN_CFI2  <   0)     LEN_CFI2 = 0
      if ( LEN_CFI3  <   0)     LEN_CFI3 = 0

! 3.2 check lookup 2nd dimensions are not too large
      if ( Len2_Lookup_P  <   Len2_Lookup_Obs ) then
        write(UnErr,*)CErr,CSub,                                        &
     &  ' step 3.2 lookup table is not big enough; Len2_Lookup_P = ',   &
     &  Len2_Lookup_P,'; Len2_Lookup = ', Len2_Lookup_Obs
        icode = 29
        go to 9999
      end if

! 3.3 check first dimensions of lookup tables match
      if ( Len1_Lookup_P  /=  Len1_Lookup_Obs ) then
        write(UnErr,*)CErr,CSub,                                        &
     &  ' step 3.3 lookup first dimensions do not match ;',             &
     &  ' Len1_Lookup_P = ', Len1_Lookup_P,'; Len1_Lookup = ',          &
     &  Len1_Lookup_Obs
        icode = 30
        go to 9999
      end if

! 3.4 set actual 2nd dimension of lookup table
      Len2_Lookup_Actual = Len2_Lookup_Obs

! 4. set position to start of file to re-read the header
      call setpos(InUnit, 0, icode)

      if ( icode  >   0) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 4. setpos failed; icode = ', icode
        icode = 31
        go to 9999
      end if

! 5. get the lookup header

! DEPENDS ON: get_lookup
      call get_lookup(InUnit, icode,                                    &
! DUMP_AR2 array dimensions
     &  LEN_FIXHD, LEN_INTHD, LEN_REALHD,                               &
     &  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,         &
     &  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC,LEN2_FLDDEPC,          &
     &  LEN_EXTCNST,  LEN_DUMPHIST,                                     &
     &  LEN_CFI1, LEN_CFI2, LEN_CFI3,                                   &
     &  LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS,                               &
! DUMP_AR2 end
     & Len_data, LOOKUP)

      if ( icode  >   0) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 5. failed read lookup table '
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE read_one_header
!----------------------------------------------------------------------
