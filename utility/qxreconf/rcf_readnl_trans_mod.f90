! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads the TRANS namelists

Module Rcf_ReadNL_Trans_Mod
IMPLICIT NONE

!  Subroutine Rcf_Readnl_Trans - reads the TRANS namelists
!
! Description:
!   Reads the TRANS namelist for controling transplanting fields.
!
! Method:
!   Reads into temporary arrays and then copies into the correct
!   dynamically allocated arrays.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_ReadNL_Trans( nft )

Use Rcf_Trans_Mod         ! All of it

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Ereport_Mod, Only : &
    Ereport

USE UM_ParVars, Only : &
    mype

IMPLICIT NONE

! Arguments
Integer, Intent(In)       :: nft      ! unit number

! Local variables
Integer, Parameter        :: max_trans = 100 ! maximum number of items
Integer                   :: iostatus
Integer                   :: itemc
Integer                   :: sctnc
Integer                   :: lev1
Integer                   :: lev2
Integer                   :: row1
Integer                   :: row2
Integer                   :: col1
Integer                   :: col2

! temporary arrays
Integer                   :: itemc_temp( max_trans )
Integer                   :: sctnc_temp( max_trans )
Integer                   :: lev1_temp( max_trans )
Integer                   :: lev2_temp( max_trans )
Integer                   :: row1_temp( max_trans )
Integer                   :: row2_temp( max_trans )
Integer                   :: col1_temp( max_trans )
Integer                   :: col2_temp( max_trans )

Integer                   :: ErrorStatus
Character (Len=*), Parameter :: RoutineName = 'Rcf_ReadNL_Trans'
Character (Len=80)        :: Cmessage

Namelist /Trans/ sctnc, itemc, lev1, lev2, col1, col2, row1, row2

!-----------------------------------------------------------------
! Not all platforms will necessarily require this rewind
!-----------------------------------------------------------------
Rewind( Unit=nft )

!-----------------------------------------------------------------
! Set num_trans to 0, loop for reading in namelist
!-----------------------------------------------------------------
num_trans = 0
iostatus  = 0

Do While ( iostatus == 0 )
  Read( Unit=nft, Nml=trans, Iostat=iostatus )

  If (iostatus == 0 ) Then
    num_trans = num_trans + 1

    If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
      Write ( Unit = 6, Nml = trans )
    End If

    If ( num_trans > max_trans ) Then
      Cmessage = 'Maximum number of transplanted'//&
          ' fields exceeded - please use a mod to increase'
      Write (6,*) Cmessage
      Write (6,*) 'Number     of TRANS namelists read in : ',num_trans
      write (6,*) 'Maximum no of TRANS namelists allowed : ',max_trans
      ErrorStatus = 10
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    If (lev1 == 0 .OR. lev2 == 0 .OR. row1 == 0 .OR. row2 == 0 .OR. &
        col1 == 0 .OR. col2 == 0 ) Then
      Write (6,*) 'Namelist TRANS has some zero values'
      Write (6,*) 'lev1 = ', lev1, ' lev2 = ', lev2
      Write (6,*) 'row1 = ', row1, ' row2 = ', row2
      Write (6,*) 'col1 = ', col1, ' col2 = ', col2

      Cmessage = 'Zero values in TRANS namelist'
      ErrorStatus = 20
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    sctnc_temp( num_trans ) = sctnc
    itemc_temp( num_trans ) = itemc
    lev1_temp( num_trans ) = lev1
    lev2_temp( num_trans ) = lev2
    col1_temp( num_trans ) = col1
    col2_temp( num_trans ) = col2
    row1_temp( num_trans ) = row1
    row2_temp( num_trans ) = row2
  End If

End Do

If ( PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
  Write (6,*) 'Num_Trans = ', num_trans
End if

!-----------------------------------------------------------------
! Allocate proper space for trans list and copy data
!-----------------------------------------------------------------
Allocate( sctnc_array( num_trans ) )
Allocate( itemc_array( num_trans ) )
Allocate( lev1_array( num_trans ) )
Allocate( lev2_array( num_trans ) )
Allocate( col1_array( num_trans ) )
Allocate( col2_array( num_trans ) )
Allocate( row1_array( num_trans ) )
Allocate( row2_array( num_trans ) )

sctnc_array( 1 : num_trans ) = sctnc_temp( 1: num_trans )
itemc_array( 1 : num_trans ) = itemc_temp( 1: num_trans )
lev1_array( 1 : num_trans ) = lev1_temp( 1: num_trans )
lev2_array( 1 : num_trans ) = lev2_temp( 1: num_trans )
col1_array( 1 : num_trans ) = col1_temp( 1: num_trans )
col2_array( 1 : num_trans ) = col2_temp( 1: num_trans )
row1_array( 1 : num_trans ) = row1_temp( 1: num_trans )
row2_array( 1 : num_trans ) = row2_temp( 1: num_trans )

Return
End Subroutine Rcf_ReadNL_Trans
End Module Rcf_ReadNL_Trans_Mod
