! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE readwritd_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE READWRITD()

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Purpose: Time-step control of WRITEDATA. Was originally read from
!          NAMLST file, but usage of this feature is being depreciated,
!          so it is now hard-coded to be switched off in this routine.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc
!
!
! This Comdeck declares and stores the logical and integer
! variables used in the time-step control for writing general
! data.
!
!
! Switch which activates output of arrays for Peer utility
      LOGICAL L_PEER

! Switches which activate dump writing
      LOGICAL                                                           &
     &  L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                  &
     &  ,L_WRIT_INIT

! Timesteps for dump writing
      INTEGER                                                           &
     &  T_WRITD1_START                                                  &
                              ! First timestep
     &  ,T_WRITD1_END                                                   &
                              ! Last timestep
     &  ,T_WRITD1_INT         ! Timestep interval between dumps

      INTEGER                                                           &
     &  PEER_VN                  ! Version of PEER utility

      NAMELIST/NLSTWRITEDATA/                                           &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT

      COMMON/WRITEDATA/                                                 &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT
!
!
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
! --------------------------------------------------------------------
!
!
!  Read control variables for WRITE_DATA
!
      IF (lhook) CALL dr_hook('READWRITD',zhook_in,zhook_handle)

! Information for time-step control of WRITEDATA no longer controlled 
! from NAMLST file read from READCNTL. Hard-code to "off" settings here: 
      L_WRIT_INIT=.FALSE. 
      T_WRITD1_START=0 
      T_WRITD1_END=0 
      T_WRITD1_INT=0 

      IF (lhook) CALL dr_hook('READWRITD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE READWRITD

END MODULE readwritd_mod
