! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: read_lsm_anc
!
! Purpose: Flux processing routine.
!          Reads in land sea mask from an ancillary file and sets
!          grid coordinates from row and column dependent constants
!          if they are present or from lookup table if not
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_lsm_anc(InUnit, Len_FixHd_P, Len1_Lookup_P,       &
     &           FixHd, Lookup, ncols, nrows, lsm, lambda, phi,         &
     &           icode)
      USE io
      USE Decomp_DB, ONLY : decompose
      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      USE lookup_addresses

      IMPLICIT NONE

! parameters used in argument list
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

! declaration of argument list
      integer InUnit        ! IN unit number of file
      integer Len_FixHd_P   ! IN length of fixed header
      integer Len1_Lookup_P ! IN length of first dimension of Lookup
      integer FixHd(Len_FixHd_P)    ! IN fixed header
      integer Lookup(Len1_Lookup_P) ! IN lookup table from file
      integer ncols       ! IN   number of columns in grid
      integer nrows       ! IN   number of rows in grid
      integer lsm(ncols,nrows) ! OUT land / sea mask
      real lambda(ncols)  ! OUT coords of longitudes
      real phi(nrows)     ! OUT coords of latitudes
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of local arrays
      real rowdepc(nrows)   ! row dependent constants
      real coldepc(ncols)   ! column dependent constants
      logical ll_lsm(ncols,nrows) ! land / sea mask T = land points
      real flt_lsm(ncols,nrows)   ! land / sea mask 1.0 = land point

! declaration of local scalars
      real DPhi    ! latitude interval
      real Phi0    ! Zeroth latitude
      real DLambda ! Zeroth longitude
      real Lambda0 ! Longitude interval

      integer jrow, icol  ! loop indices for rows and columns
      integer fld_no   ! field number in file
      integer Len_data  ! length of data in file
      integer gl_row_len  ! # of columns in global grid
      integer gl_n_rows   ! # of rows in global grid

      CHARACTER(LEN=256) CMessage ! error messages

      external copy_to_real
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_lsm_anc'  ! subroutine name for error messages


! 1. get dimensions from lookup table

      call setpos(InUnit, 0, icode)   ! reset to start of file
      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &   ' Step 1. Unable to reset ancillary file'
        icode = 20
        goto 9999
      endif

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

! 2. if row and column dependent constants have non-zero length
!    extract them from the lsm file

! N.B. this code is schematic only at this stage and would need
!      proper checking if it was to be used !!!!
!      I assume that rowdepc and coldepc are set up as implied in
!      section 4. below.

      if ( LEN1_COLDEPC  >   1 .and. LEN1_ROWDEPC  >   1) then

! 2.1 check that the lengths match nrows and ncols

        if (      nrows  /=  LEN1_ROWDEPC                               &
     &       .or. ncols  /=  LEN1_COLDEPC) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.1. dimensions do not match  '
          go to 9999
        end if

! 2.2 get the row and column dependent constants (which are row
!     and column spacings

        call setpos(InUnit, 0, icode)
        if ( icode  /=  0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &     ' Step 2.2. Unable to reset ancillary file'
          icode = 20
          goto 9999
        endif

! DEPENDS ON: get_rows_cols
        call get_rows_cols(InUnit, icode,                               &
! DUMP_AR2 array dimensions
     &  LEN_FIXHD, LEN_INTHD, LEN_REALHD,                               &
     &  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,         &
     &  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC,LEN2_FLDDEPC,          &
     &  LEN_EXTCNST,  LEN_DUMPHIST,                                     &
     &  LEN_CFI1, LEN_CFI2, LEN_CFI3,                                   &
     &  LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS,                               &
! DUMP_AR2 end
     &  Len_data, rowdepc, coldepc)

        if ( icode  /=  0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &     ' Step 2.2. Failed to retreive grid spacings'
          goto 9999
        endif

! 3. else take the row and column spacings from the lookup table
      else   !  LEN1_COLDEPC etc.

! DEPENDS ON: copy_to_real
        call copy_to_real ( Lookup(BDY), DPhi )
        do jrow = 1, nrows
          rowdepc(jrow) = DPhi
        end do

! DEPENDS ON: copy_to_real
        call copy_to_real ( Lookup(BDX), DLambda )
        do icol = 1, ncols
          coldepc(icol) = DLambda
        end do

      end if ! LEN1_COLDEPC etc.


! 4. copy to row and column coordinates
! DEPENDS ON: copy_to_real
      call copy_to_real ( Lookup(BZY), Phi0 )
      Phi(1) = Phi0 + rowdepc(1)
      do jrow = 2, nrows
        Phi(jrow) = Phi(jrow-1) + rowdepc(jrow)
      end do

! DEPENDS ON: copy_to_real
      call copy_to_real ( Lookup(BZX), Lambda0 )
      Lambda(1) = Lambda0 + coldepc(1)
      do icol = 2, ncols
        Lambda(icol) = Lambda(icol-1) + coldepc(icol)
      end do

! 5. read in the land sea mask itself; assumed to be the first field
!    in the file

      fld_no = 1

! Reset STASH code to avoid failure when reading field
      lookup(ITEM_CODE)=101

! 5.1 if land sea mask data is of type logical read it into a
!     temporary array and convert T = land => 1 and F = sea => 0

      gl_row_len=lookup(19)
      gl_n_rows=lookup(18)

      CALL decompose(gl_row_len, gl_n_rows,0,0,1)

      if ( Lookup(data_type)  ==  3) then  !  logical

! DEPENDS ON: readflds
        call readflds (InUnit , 1, fld_no, Lookup,                      &
     &  Len1_Lookup_P, ll_lsm, ncols*nrows, FIXHD,                      &
     &  icode, cmessage)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 5.1 unable to read logicals land sea mask: ',          &
     &    ' cmessage is ', cmessage
          icode = 23
          go to 9999
        end if

        do jrow = 1, nrows
          do icol = 1, ncols
            if ( ll_lsm(icol, jrow) ) then
              lsm(icol, jrow) = 1
            else
              lsm(icol, jrow) = 0
            end if
          end do  ! icol
        end do  ! jrow

! 5.2 else if land sea mask data is of type integer read it

      else if( Lookup(data_type)  ==  2) then  !  integer

! DEPENDS ON: readflds
        call readflds (InUnit , 1, fld_no, Lookup,                      &
     &  Len1_Lookup_P, lsm, ncols*nrows, FIXHD,                         &
     &  icode, cmessage)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 5.2 unable to read integer land sea mask:',            &
     &    ' cmessage is ',   cmessage
          icode = 24
          go to 9999
        end if

! 5.3 else if land sea mask data is of type real, read it and
!     convert: 1.0 = land => 1  0.0 = sea => 0

      else if( Lookup(data_type)  ==  1) then  !  real
! DEPENDS ON: readflds
        call readflds (InUnit , 1, fld_no, Lookup,                      &
     &  Len1_Lookup_P, flt_lsm, ncols*nrows, FIXHD,                     &
     &  icode, cmessage)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 5.2 unable to read real land sea mask:',               &
     &    ' cmessage is ',   cmessage
          icode = 25
          go to 9999
        end if

        do jrow = 1, nrows
          do icol = 1, ncols
            if ( flt_lsm(icol, jrow)  /=  0.0 ) then
              lsm(icol, jrow) = 1
            else
              lsm(icol, jrow) = 0
            end if
          end do  ! icol
        end do  ! jrow

! 5.4 else there is an error in data type of land sea mask

      else
        icode = 26
        write(UnErr,*)CErr,CSub,                                        &
     &  ' step 5.3 land sea mask is of data type:', Lookup(data_type)   &
     &  , '. Change this to indicator for integer or logical data'

      end if   ! Lookup(data_type)

9999  continue
      return
      END SUBROUTINE read_lsm_anc
!----------------------------------------------------------------------
