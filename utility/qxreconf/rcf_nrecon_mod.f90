! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Output dump addressing arrays

Module Rcf_NRecon_Mod

! Description:
!   Arrays used for output dump addressing.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None. R.Sharp
!   6.1   22/06/04   Extend Recondat for Tracer Prognostics. R Barnes
!   6.2   10/11/05   Extend Recondat for UKCA Prognostics. R Barnes
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Use version_mod, Only       : &
    Nitemp,                   &
    NSectP

Use Submodel_Mod, Only  :    &
    N_Internal_Model_Max

Implicit None

Type :: Type_Recondat
  Integer                    :: sec_item
  Integer                    :: rlevs
  Integer                    :: len
  Integer                    :: raddress
  Integer                    :: rplevs
End Type Type_Recondat

Type :: Type_Recondat_Node
  Type ( Type_Recondat ), Pointer    :: recondat_info => Null()
  Type ( Type_Recondat_Node ), Pointer :: next => Null()
End Type Type_Recondat_Node

Type ( Type_Recondat_Node ), Pointer :: recondat_node => Null()

Type ( Type_Recondat_Node ), Save, Target &
  :: RecondatList( N_Internal_Model_Max, 0 : NSectP )
Integer                     :: PrimDataLen ( N_Internal_Model_Max )
Integer                     :: DumpProgLevs( N_Internal_Model_Max )

End Module Rcf_NRecon_Mod
