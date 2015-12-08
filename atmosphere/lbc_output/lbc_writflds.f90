! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Writes out LBC data
!
! Subroutine Interface:

      Subroutine lbc_writflds (                                         &
     &           nftout,lbc_data,len_lbc_data,lookup,fixhd,ltimer )

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE IO
      USE ereport_mod, ONLY : ereport
      USE PrintStatus_mod
      USE UM_ParVars
      USE lookup_addresses

      IMPLICIT NONE

!
! Description:
!   Writes out LBC data to a boundary file. Data is packed first
!   if required.
!
! Method:
!   1. Pack data to 32 bits if required
!   2. Write data out using BUFFOUT.
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
! Subroutine arguments

      Integer :: nftout       ! Unit no for boundary file
      Integer :: len_lbc_data ! Length of data buffer
      Integer :: lookup(64)   ! Lookup header for data field
      Integer :: fixhd(256)   ! Fixed header

      Real    :: lbc_data(len_lbc_data) !  LBC data

      Logical :: ltimer       ! controls call to timer
! Local parameters:

      Character (Len=*), Parameter :: RoutineName= 'LBC_Writflds'

! Local scalars:

      Integer   :: disk_address  ! disk address
      Integer   :: len_data_disk ! length of data to be written out
      Integer   :: len_data      ! exact length of data
      Integer   :: pack_data     ! packing indicator
      Integer   :: len_io        ! length of data actually written
      Real      :: a_io          ! return code from buffout
      Integer   :: ErrorStatus   ! error flag
      Character (Len=80)  :: CMessage   !  error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Function & Subroutine calls:

!- End of header

! -------------------------------------------
! Get required information from lookup header
! -------------------------------------------
      IF (lhook) CALL dr_hook('LBC_WRITFLDS',zhook_in,zhook_handle)
      disk_address  = lookup (lbegin)
      len_data_disk = lookup (lbnrec)
      len_data      = lookup (lblrec)
      pack_data     = lookup (lbpack)

!     len_data_disk includes rounding up to sector boundaries

!     Add check that len_lbc_data >= len_data_disk ?

! -------------
! Pack the data
! -------------

      If (Mod (pack_data,10) == 2) Then

        If (mype == 0) Then
! DEPENDS ON: pack21
          Call Pack21 (len_data,lbc_data,lbc_data)
        Endif

        If (PrintStatus >= PrStatus_Diag) Then
          write (6,*) ' Packing LBC_data for Item Code ',               &
     &    lookup(item_code)
        Endif

      Endif

! ------------------
! Write out the data
! ------------------

      If (PrintStatus >= PrStatus_Diag) Then
        write (6,*) ' Writing LBC_data for Item Code ',                 &
     &  lookup(item_code),' length of data ',len_data_disk
      Endif


      call setpos  (nftout, disk_address, ErrorStatus)


      If (Mod (pack_data,10) == 2) Then
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*len_data_disk 32 bit words using BUFFOUT32_F77 (because we supply
! a 64 bit array)

! DEPENDS ON : buffout32_f77
        call buffout32_f77 &
            (nftout, lbc_data, 2*len_data_disk, len_io, a_io)
! And then halve len_io to satisfy tests against len_data_disk
        len_io = len_io/2
      Else
! For non-packed data

        call buffout (nftout, lbc_data, len_data_disk, len_io, a_io)
      Endif





      If (a_io /= -1.0) Then

        Write (6,*) ' Return Code from BUFFOUT    ',a_io
        Write (6,*) ' Length of data transferred  ',len_io
        Write (6,*) ' Expected transferred length ',len_data_disk

        ErrorStatus = 10
        write (cmessage,*) ' Problem in BUFFOUT for LBCs : Item ',      &
     &  lookup(item_code),' being written out.'

        call ereport ( RoutineName, ErrorStatus, CMessage )

      Endif

      IF (lhook) CALL dr_hook('LBC_WRITFLDS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE lbc_writflds
