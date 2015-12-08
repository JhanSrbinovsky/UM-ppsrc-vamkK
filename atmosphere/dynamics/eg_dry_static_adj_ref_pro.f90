! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_dry_static_adj_ref_pro_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_dry_static_adj_ref_pro(theta)

USE parkind1,              ONLY : jpim, jprb       !DrHook
USE yomhook,               ONLY : lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod

IMPLICIT NONE
!
! Description:
!          Enforces static stability on input theta field.
!          Uses insertion sort algorithm resulting in no mixing of
!          temperature just a simple re-arrangement.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

REAL :: theta (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end) 

! Local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                                                               &
  i                                                                   &
, j                                                                   &
, k                                                                   &
                    ! Loop indices
, cnt    

REAL :: temp

! No External Routines

IF (lhook) CALL dr_hook('EG_DRY_STATIC_ADJ_REF_PRO',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 1. Check for inertial instability and remove.
!            Insertion sort algorithm used.    
!            No print output as algorithm doesn't make note of points modified.
!            No error checking as algorithm guarantees a sorted field.
! ----------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,cnt,temp)  &
!$OMP&  SHARED(tdims,theta)
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
        DO k=1,tdims%k_end 
          temp=theta(i,j,k)        !  save to make a hole at index k
          cnt=k
         ! keep moving the hole to next smaller index until theta(i,j,cnt-1)
         ! is <= theta(i,j,k)
          DO WHILE( theta(i,j,cnt-1) > temp)
            theta(i,j,cnt)=theta(i,j,cnt-1)  ! move hole to next smaller index
            cnt=cnt-1
            IF (cnt == 0) EXIT
          END DO
          theta(i,j,cnt)=temp    ! put item in the hole
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

! End of routine
IF (lhook) CALL dr_hook('DRY_STATIC_ADJ',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_dry_static_adj_ref_pro
END MODULE eg_dry_static_adj_ref_pro_mod
