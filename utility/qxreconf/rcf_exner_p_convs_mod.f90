! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Converts between P and exner (on rho levels)

Module Rcf_Exner_P_Convs_Mod

!  Subroutine Rcf_Conv_Exner_P - converts exner to P
!  Subroutine Rcf_Conv_P_Exner - converts P to exner
!
! Description:
!   Performs conversions between exner and P on rho levels.
!
! Method:
!   Data is stored in the *original* field data pointer!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE
Contains

!******************************************************************
! This routine converts Exner to P
!******************************************************************
Subroutine Rcf_Conv_Exner_P( exner_field )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_prog_sec

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Submodel_Mod, Only :  &
    Atmos_IM
    
USE atmos_constants_mod, ONLY: pref, recip_kappa    

IMPLICIT NONE

! Arguments
Type( field_type), Intent( InOut )    :: exner_field


!  Local variables
Integer                        :: i
Integer                        :: k
Integer                        :: ErrorStatus
Character (Len=*), Parameter   :: RoutineName = 'Rcf_Conv_Exner_P'
Character (Len=80)             :: Cmessage

!-----------------------------------------------------------------
! Make sure field actually is exner...
!-----------------------------------------------------------------
If ( exner_field % stashmaster % section /= stashcode_prog_sec .OR.  &
     exner_field % stashmaster % item    /= stashcode_exner ) Then
  ErrorStatus = 10
  Cmessage = 'Stashcode does not match expected (exner) data'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!------------------------------------------------------------------
! Convert exner to P
!------------------------------------------------------------------
Do k = 1, exner_field % levels
  Do i = 1, exner_field % level_size
    exner_field % Data(i,k) =                                       &
                  (exner_field % Data(i,k) ** recip_kappa)* pref
  End Do
End Do

!------------------------------------------------------------------
! Set field stashmaster to P
!------------------------------------------------------------------
exner_field % stashmaster => Rcf_Exppx( Atmos_IM, stashcode_prog_sec, &
                                        stashcode_p )

Return
End Subroutine Rcf_Conv_Exner_P


!*******************************************************************
! Routine to convert P to exner
!*******************************************************************

Subroutine Rcf_Conv_P_Exner( exner_field )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_prog_sec

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Submodel_Mod, Only : &
    Atmos_IM

USE atmos_constants_mod, ONLY: pref, kappa

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Arguments
Type( field_type), Intent( InOut )    :: exner_field

!  Local variables
Integer                        :: i
Integer                        :: k
Integer                        :: ErrorStatus
Real, Parameter                :: recip_p_zero=1./pref
Character (Len=*), Parameter   :: RoutineName = 'Rcf_Conv_P_Exner'
Character (Len=80)             :: Cmessage

!-----------------------------------------------------------------
! Make sure field actually is P...
!-----------------------------------------------------------------
If ( exner_field % stashmaster % section /= stashcode_prog_sec .OR.   &
     exner_field % stashmaster % item    /= stashcode_p ) Then
  ErrorStatus = 10
  Cmessage = 'Stashcode does not match expected (pressure) data'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!------------------------------------------------------------------
! Convert P to Exner
!------------------------------------------------------------------
Do k = 1, exner_field % levels
  Do i = 1, exner_field % level_size
    exner_field % Data(i,k) =                                       &
                  (exner_field % Data(i,k) * recip_p_zero) ** kappa
  End Do
End Do

!------------------------------------------------------------------
! Set field stashmaster to exner
!------------------------------------------------------------------
exner_field % stashmaster => Rcf_Exppx(Atmos_IM, stashcode_prog_sec,  &
                                       stashcode_exner)

Return
End Subroutine Rcf_Conv_P_Exner

End Module Rcf_Exner_P_Convs_Mod
