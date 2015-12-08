! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Description:
!   Provides data to client routines relevent to describe coupling to
!   the parallel server
!
! Method:
!   In the current version we only decompose by row-blocks, and thus this
!   coupling from the client side ammounts to the id of the owning
!   server component.

MODULE IOS_Client_Coupler
  IMPLICIT NONE
  INTEGER              :: target_ios_rank
END MODULE IOS_Client_Coupler
