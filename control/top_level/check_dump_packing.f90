! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Checks the packing codes in a model dump.
!
! Subroutine Interface:

      Subroutine Check_Dump_Packing (                                   &
     &           FixHd, Len_FixHd,                                      &
     &           Lookup, Len1_Lookup, Len2_Lookup,                      &
     &           Dump_Pack, IM_Ident )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE PrintStatus_mod
      USE lookup_addresses
      USE cppxref_mod, ONLY: ppx_dump_packing
      USE Submodel_Mod

      IMPLICIT NONE
!
! Description:
!   Checks consistency of packing codes in a model dump.
!
! Method:
!   1. Packing in a dump is controlled through DUMP_PACKim
!   2. The packing codes in the LookUp table (Word 21) is checked for
!      consistency with DUMP_PACKim.
!   3. If inconsistent, the packing code in LookUp table is updated
!      to match DUMP_PACKim.
!   4. Only 32 bit packing or no packing is catered for.
!      (ie. N1=0 or N1=2 , see UMDP F3)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Subroutine arguments

      Integer :: Len_FixHd
      Integer :: Len1_Lookup
      Integer :: Len2_Lookup
      Integer :: FixHd (Len_FixHd)
      Integer :: LookUp(Len1_Lookup,Len2_Lookup)
      Integer :: Dump_Pack
      Integer :: IM_Ident

! Comdecks/common blocks/parameters

! Diagnostic message switch

! Local variables

      Integer, Parameter :: No_packing = 0
      Integer, Parameter :: Pack_32bit = 2

      Integer            :: i        ! Loop index
      Integer            :: item     ! Stash Item No
      Integer            :: sect     ! Stash Section No
      Integer            :: number_of_data_words_in_memory
      Integer            :: number_of_data_words_on_disk
      Integer            :: disk_address
      Integer            :: pack_code       ! Packing code in lookup
      Integer            :: ErrorStatus

      Logical            :: prognostic      ! T : Item is prognostic
      Logical            :: diagnostic      ! T : Item is diagnostic
      Logical            :: packing_changed ! T : packing codes changed
      Logical            :: real_data       ! T : real data

      Character (Len=80) :: CMessage
      Character (Len=*), Parameter :: RoutineName='Check_Dump_Packing'

! Function
      Integer :: exppxi

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!- End of header

      IF (lhook) CALL dr_hook('CHECK_DUMP_PACKING',zhook_in,zhook_handle)
      ErrorStatus = 0
      CMessage = ' '

      packing_changed = .false.

      Do i = 1, Len2_Lookup

!       Extract details from LookUp Table

        pack_code = Mod (Lookup(lbpack,i), 10)
        sect      = Lookup(item_code,i) / 1000
        item      = Mod ( Lookup(item_code,i) ,1000 )
        real_data = (Lookup(data_type,i) == 1)

        prognostic = (sect == 0)
        diagnostic = (sect >  0)

        Select Case (pack_code)  !  Packing code in LookUp Table

          Case (No_Packing)

            If ( Real_data      .and.                                   &
     &         ( Dump_Pack == 1 .or.                                    &
                                                  ! Prog & Diag packed
     &          (Dump_Pack == 2 .and. Diagnostic))                      &
                                                  ! Diag packed
     &         ) Then

!             Set packing indicator according to StashMaster record

              Lookup(lbpack,i) =                                        &
! DEPENDS ON: exppxi
     &               exppxi (im_ident, sect, item, ppx_dump_packing,    &
     &               errorstatus, cmessage)

              packing_changed = .true.

            End If

          Case (Pack_32bit)

            If (  Real_data      .and.                                  &
     &         ( (Dump_Pack == 2 .and. Prognostic) .or.                 &
                                                        ! Prog unpacked
     &            Dump_Pack == 3 )                                      &
                                                 ! Prog & Diag unpacked
     &         ) Then

!             Set packing indicator to no packing

              Lookup(lbpack,i) = ( Lookup(lbpack,i)/10 ) * 10

              packing_changed = .true.

            End If

          Case Default

            ErrorStatus = 10
            Write(CMessage,*)' Unexpected Packing code in dump.',       &
     &                       ' Field No ',i,' Packing Code ',Pack_Code

            Call EReport (RoutineName, ErrorStatus, CMessage)

        End Select

      End Do

!     -------------------------------------------------------
!     If any of the packing indicators have been changed then
!     reset the disk addresses and lengths
!     -------------------------------------------------------

      If (packing_changed) Then

! DEPENDS ON: set_dumpfile_address
        Call Set_DumpFile_Address (FixHd, Len_FixHd,                    &
     &                             Lookup, Len1_Lookup, Len2_Lookup,    &
     &                             number_of_data_words_in_memory,      &
     &                             number_of_data_words_on_disk,        &
     &                             disk_address)

        If (PrintStatus >= PrStatus_Normal) Then
          ErrorStatus = -20   !  Warning
          Write (CMessage,*) ' Packing codes in dump inconsistent',     &
     &    ' with DUMP_PACKim. Packing codes updated.'

          Call EReport (RoutineName, ErrorStatus, CMessage)
        End If

      End If

      IF (lhook) CALL dr_hook('CHECK_DUMP_PACKING',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Check_Dump_Packing
