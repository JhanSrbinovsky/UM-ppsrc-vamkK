! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Top level program for reconfiguration
! Description:
!  Top level program for reconfiguration
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Program Reconfigure

Use Rcf_Initialise_Mod, Only : &
    Rcf_Initialise

Use Rcf_Finalise_Mod, Only : &
    Rcf_Finalise

Use Rcf_Read_Namelists_Mod, Only : &
    Rcf_Read_Namelists

Use Rcf_Control_Mod, Only : &
    Rcf_Control

Use Rcf_UMhead_Mod, Only : &
    um_header_type

USE UM_ParVars, Only : &
    mype,               &
    nproc

Use Ereport_mod, Only :&
    Ereport, Ereport_finalise

Use IO, Only : ioShutdown,ioInit

USE UM_Config, ONLY : &
    appInit, &
    exe_RCF

Implicit None

Type (um_header_type)        :: hdr_in       ! header from input dump
Type (um_header_type)        :: hdr_out      ! header from output dump
Integer                      :: ErrorStatus
Integer                      :: info         ! GCOM dummy
Character (Len=*), Parameter :: RoutineName = 'Reconfigure'
Character (Len=80)           :: Cmessage

Integer, Parameter :: gc_alltoall_version = 2    ! From GCOM 
Integer, Parameter :: gc_alltoall_multi   = 2    ! From GCOM 

!------------------------------------------------------------------
! Initialise GCOM
!------------------------------------------------------------------
CALL gc_init( ' ', mype, nproc )
If ( nproc < 0 ) Then
  ErrorStatus = 20
  Cmessage = 'Parallel Initialisation Failed!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If
CALL appInit(exe_RCF)

! Set GCOM to use the alternative version of RALLTOALLE 
! throughout the run 
Call Gc_Setopt(gc_alltoall_version, gc_alltoall_multi, Errorstatus) 

Call ioInit()

!-------------------------------------------------------------------
! Perform initialisation
!-------------------------------------------------------------------
Call Rcf_Initialise( hdr_in, hdr_out )

!-------------------------------------------------------------------
! Do the real work
!-------------------------------------------------------------------
Call Rcf_Control( hdr_in, hdr_out )

!------------------------------------------------------------------
! Tidy Up
!------------------------------------------------------------------
Call Rcf_Finalise( hdr_in, hdr_out )

!------------------------------------------------------------------
! Report on IO if configured
!------------------------------------------------------------------

Call ioShutdown()

CALL Ereport_finalise( )

Write (6,*) ' End of rcf program reached. PE ',mype
!------------------------------------------------------------------
! Kill off gcom
!------------------------------------------------------------------
Call gc_exit()

End Program Reconfigure
