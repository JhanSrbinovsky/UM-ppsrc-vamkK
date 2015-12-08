! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reformates a field from rho points to uv position on 'C' grid
!
! Subroutine Interface:
MODULE pc_to_pb_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE pc_to_pb(pc,                                           &
     & row_length,rows,n_rows,levels,offx,offy,                         &
     & pb)


      USE atm_fields_bounds_mod, ONLY:  &
                udims, udims_s, vdims, vdims_s


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      IMPLICIT NONE
!
! Description:
!   pC_to_pB performs a simple horizontal interpolation of a field pC
!   at the p position for Arakawa 'C' grid staggering onto the uv
!   position of the 'B' grid.
! Method:
!   1. Copy field (without halos) to an array with halos.
!   2. Populate halos using swapbounds.
!   3. Make a simple average of adjacent rows and copy to output field
!     (without halos).
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Climate Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & row_length,rows                                                  &
                        ! horizontal dimensions
     &,n_rows                                                           &
                        ! rows for last (N) row of pes
     &,levels                                                           &
                        ! vertical levels
     &,offx,offy        ! fixed halo sizes

!   Array  arguments with intent(in):
      REAL                                                              &
       pc(row_length,rows,levels) ! Field on p points on C grid

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL                                           &
       pb(udims%i_start:udims%i_end,          &  
          vdims%j_start:vdims%j_end,          &  
          levels) !                 ! Field on uv points on B grid

! Local parameters:
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='pc_to_pb')

! Local scalars:
      INTEGER                                                           &
     & i,j,k     !  loop indices

! Local dynamic arrays:
! pc_halo dims are valid for both EG and ND grids.
      REAL                                                              &
     & pc_halo(1-offx:row_length+offx,1-offy:rows+offy,levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('PC_TO_PB',zhook_in,zhook_handle)

!- End of header

!   1. Copy field (without halos) to an array with halos
      DO k=1,levels
        DO j=1,rows
          DO i=1,row_length
            pc_halo(i,j,k)=pC(i,j,k)
          END DO                 ! i
        END DO                   ! j
      END DO                     ! k

!   2. Populate halos using swapbounds.

! DEPENDS ON: swap_bounds
      CALL Swap_bounds(pc_halo,row_length,rows,levels,                  &
                       offx,offy,fld_type_p,.FALSE.)

!   3. Make a simple average of adjacent rows and copy to output field
!     (without halos).

      DO k=1,levels
!CDIR NOUNROLL
        DO j=vdims%j_start,vdims%j_end
          DO i=udims%i_start,udims%i_end
            pb(i,j,k)=(pc_halo(i,j,k)+pc_halo(i,j+1,k)+                 &
              pc_halo(i+1,j,k)+pc_halo(i+1,j+1,k)) * 0.25
          END DO                 ! i
        END DO                   ! j
      END DO                     ! k

      IF (lhook) CALL dr_hook('PC_TO_PB',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE pc_to_pb
END MODULE pc_to_pb_mod
