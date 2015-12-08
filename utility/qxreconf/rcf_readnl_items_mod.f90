! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads the ITEMS namelists

MODULE Rcf_readnl_items_Mod
IMPLICIT NONE

!  Subroutine Rcf_Readnl_Items - Reads the ITEMS namelists
!
! Description:
!   Reads the ITEMS namelists to control the source of required fields.
!
! Method:
!   Read initially into fixed size array and transferred to
!   to correctly sized dynamically allocated array later.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

Subroutine rcf_readnl_items (nft)

Use Rcf_Items_Mod, Only : &
    Source_Array,         &
    Sctn_Array,           & 
    Item_Array,           &
    Area_Array,           &
    Upas_Array,           &
    Upaa_Array,           &
    Uprc_Array,           &
    Upaf_Array,           &
    Num_Items

USE PrintStatus_mod, Only : &
    PrintStatus,          &
    PrStatus_Oper

Use Ereport_Mod, Only : &
    Ereport

USE UM_ParVars, Only : &
    mype

USE rcf_recon_mod, ONLY : &
    l_interp_input_only

IMPLICIT NONE

! Arguments
Integer      ::  nft

! Local variables
! increase number of STASH items
      integer, Parameter      :: max_items = 10000 ! max no of items

! temporary arrays
Integer                   :: sctn_temp ( max_items )
Integer                   :: item_temp ( max_items )
Integer                   :: srce_temp ( max_items )
Integer                   :: area_temp ( max_items )
Integer                   :: upas_temp ( max_items )
Integer                   :: upaa_temp ( max_items )
Real                      :: uprc_temp ( max_items )
Character (Len=80)        :: upaf_temp ( max_items )

! Namelist ITEMS
Integer  :: Section  !  Section (Stash) Code
Integer  :: Item     !  Item (Stash) Code
Integer  :: Source   !  Source of data
Integer  :: Domain   !  Areal coverage of data
Integer  :: User_Prog_Ancil_SctnC
Integer  :: User_Prog_Ancil_ItemC
Real     :: User_Prog_RConst
Character (Len=80) :: User_Prog_Ancil_File

Namelist /ITEMS/                                                  &
   Section, Item, Source, Domain, User_Prog_Ancil_SctnC,          &
   User_Prog_Ancil_ItemC, User_Prog_RConst, User_Prog_Ancil_File

Integer  ::      Errorstatus
Integer  ::      iostatus
Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Readnl_Items'

num_items = 0
iostatus  = 0

Do While ( iostatus == 0 .AND. .NOT. l_interp_input_only)

! Initialise namelist
  Section = 0
  Item    = 0
  Source  = 0
  Domain  = 0
  User_Prog_Ancil_SctnC = 0
  User_Prog_Ancil_ItemC = 0
  User_Prog_Ancil_File  = ' '
  User_Prog_RConst      = 0.0

  Read( Unit=nft, Nml=items, Iostat=iostatus )

  If (iostatus == 0 ) Then
    num_items = num_items + 1

    If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
      Write ( Unit = 6, Nml = items )
    End If

    If ( num_items > max_items ) Then
      Cmessage = 'Maximum number of ITEMS namelists '//&
          'exceeded - please use a mod to increase.'
      Write (6,*) Cmessage
      Write (6,*) 'Number     of ITEMS namelists read in : ',num_items
      write (6,*) 'Maximum no of ITEMS namelists allowed : ',max_items
      ErrorStatus = 10
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    sctn_temp ( num_items ) = section
    item_temp ( num_items ) = item
    srce_temp ( num_items ) = source
    area_temp ( num_items ) = domain
    upas_temp ( num_items ) = user_prog_ancil_sctnc
    upaa_temp ( num_items ) = user_prog_ancil_itemc
    upaf_temp ( num_items ) = user_prog_ancil_file
    uprc_temp ( num_items ) = user_prog_rconst

  End If

End Do

If ( PrintStatus >= PrStatus_Oper ) Then
  write (6,*) ' No of ITEMS namelists read in : ',num_items
End If

if ( Num_Items > 0 ) then

!-----------------------------------------------------------------
! Allocate proper space for items list and copy data
!-----------------------------------------------------------------
Allocate( sctn_array  ( num_items ) )
Allocate( item_array  ( num_items ) )
Allocate( source_array( num_items ) )
Allocate( area_array  ( num_items ) )
Allocate( uprc_array  ( num_items ) )
Allocate( upaf_array  ( num_items ) )
Allocate( upaa_array  ( num_items ) )
Allocate( upas_array  ( num_items ) )

sctn_array  ( 1 : num_items ) = sctn_temp( 1: num_items )
item_array  ( 1 : num_items ) = item_temp( 1: num_items )
source_array( 1 : num_items ) = srce_temp( 1: num_items )
area_array  ( 1 : num_items ) = area_temp( 1: num_items )
uprc_array  ( 1 : num_items ) = uprc_temp( 1: num_items )
upaf_array  ( 1 : num_items ) = upaf_temp( 1: num_items )
upaa_array  ( 1 : num_items ) = upaa_temp( 1: num_items )
upas_array  ( 1 : num_items ) = upas_temp( 1: num_items )

endif  !  if num_items > 0

return

END SUBROUTINE  Rcf_readnl_items

END MODULE Rcf_readnl_items_Mod
