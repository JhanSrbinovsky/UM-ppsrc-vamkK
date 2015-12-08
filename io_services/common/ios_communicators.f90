! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Module defining communicators used by the IOS.

MODULE IOS_communicators
  IMPLICIT NONE

! Global: all processes, atmos and all all IOS roles 
  INTEGER :: global_rank
  INTEGER :: global_comm
  INTEGER :: global_procs

! model: the IO server ot atmosphere team I am in
  INTEGER :: model_rank
  INTEGER :: model_comm
  INTEGER :: model_procs
! IO: the set of all IO related processes
  INTEGER :: io_rank
  INTEGER :: io_comm
  INTEGER :: io_procs

! Leaders : The set of rank 0 in models
  INTEGER :: leader_rank
  INTEGER :: leader_comm
  INTEGER :: leader_procs
  
! Bcast : the leader of an IOS and all atmos
  INTEGER :: bcast_rank
  INTEGER :: bcast_procs

  INTEGER          :: IOS_BCast_Server_Comm
  INTEGER, POINTER :: IOS_BCast_Comm(:)=>NULL()


END MODULE IOS_communicators

