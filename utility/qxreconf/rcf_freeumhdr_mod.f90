! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Clears space from UM_header_type

Module Rcf_FreeUMhdr_Mod

!  Subroutine Rcf_FreeUMhdr - frees up used arrays etc.
!
! Description:
!   Allocated pointers to hold header components are released and
!   sizes set to MDI as required.
!
! Method:
!   Deallocate Pointers, Set sizes to IMDI, Nullify Pointers.
!
! Inspired by VAR code.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None. R.Sharp
!   5.4   12/06/02   GRIB Prep. Protected Deallocate statements with
!                    an If allocated check. Rod.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains

Subroutine Rcf_FreeUMhdr ( UMhdr )

Use Rcf_UMhead_Mod, Only : &
    UM_header_type

Implicit None

! Arguments
Type (UM_header_type), Intent( InOut ) ::  UMhdr

! Comdecks
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!-----------------------------------------------------------------
! 1. Set constants to IMDI
!-----------------------------------------------------------------
UMhdr % LenIntC      = IMDI
UMhdr % LenRealC     = IMDI
UMhdr % Len1LevDepC  = IMDI
UMhdr % Len2LevDepC  = IMDI
UMhdr % Len1RowDepC  = IMDI
UMhdr % Len2RowDepC  = IMDI
UMhdr % Len1ColDepC  = IMDI
UMhdr % Len2ColDepC  = IMDI
UMhdr % Len1FldsOfC  = IMDI
UMhdr % Len2FldsOfC  = IMDI
UMhdr % LenExtraC    = IMDI
UMhdr % LenHistFile  = IMDI
UMhdr % LenCompFldI1 = IMDI
UMhdr % LenCompFldI2 = IMDI
UMhdr % LenCompFldI3 = IMDI
UMhdr % Len1Lookup   = IMDI
UMhdr % Len2Lookup   = IMDI
UMhdr % StartData    = IMDI
UMhdr % LenData      = IMDI
UMhdr % NumFlds      = IMDI
UMhdr % MaxFldSize   = IMDI
UMhdr % UnitNum      = IMDI

!--------------------------------------------------------------
! 2. Deallocate pointers
!--------------------------------------------------------------

If (Associated( UMhdr % FixHd     )) Deallocate( UMhdr % FixHd     )
If (Associated( UMhdr % IntC      )) Deallocate( UMhdr % IntC      )
If (Associated( UMhdr % CompFldI1 )) Deallocate( UMhdr % CompFldI1 )
If (Associated( UMhdr % CompFldI2 )) Deallocate( UMhdr % CompFldI2 )
If (Associated( UMhdr % CompFldI3 )) Deallocate( UMhdr % CompFldI3 )
If (Associated( UMhdr % Lookup    )) Deallocate( UMhdr % Lookup    )
If (Associated( UMhdr % RealC     )) Deallocate( UMhdr % RealC     )
If (Associated( UMhdr % LevDepC   )) Deallocate( UMhdr % LevDepC   )
If (Associated( UMhdr % RowDepC   )) Deallocate( UMhdr % RowDepC   )
If (Associated( UMhdr % ColDepC   )) Deallocate( UMhdr % ColDepC   )
If (Associated( UMhdr % FldsOfC   )) Deallocate( UMhdr % FldsOfC   )
If (Associated( UMhdr % ExtraC    )) Deallocate( UMhdr % ExtraC    )
If (Associated( UMhdr % HistFile  )) Deallocate( UMhdr % HistFile  )

!---------------------------------------------------------------
! 3. Nullify pointers (for safety...)
!---------------------------------------------------------------

Nullify( UMhdr % Fixhd     )
Nullify( UMhdr % IntC      )
Nullify( UMhdr % CompFldI1 )
Nullify( UMhdr % CompFldI2 )
Nullify( UMhdr % CompFldI3 )
Nullify( UMhdr % Lookup    )
Nullify( UMhdr % RealC     )
Nullify( UMhdr % LevDepC   )
Nullify( UMhdr % RowDepC   )
Nullify( UMhdr % ColDepC   )
Nullify( UMhdr % FldsofC   )
Nullify( UMhdr % ExtraC    )
Nullify( UMhdr % HistFile  )

Return

End Subroutine Rcf_FreeUMhdr

End Module Rcf_FreeUMhdr_Mod
