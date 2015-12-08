! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: set_lsmu
!
! Purpose: Flux processing routine.
!          Sets u velocity atmosphere land / sea mask
!
! In the atmosphere code C_D | dv/dz | is linearly interpolated from
! tracer points to velocity points. Velocity points next to a land
! point  are thus "contaminated" with these values from land points
! which could be much larger than those at sea points. So only velocity
! points surrounded by tracer sea points are treated as land points
! in code to form fluxes for the operational FOAM model.
!
! The routines for rotating wind components assume that they are both
! defined at the same points - as they are on a B-grid.
!
! We have chosen to average the C-grid velocities to the B-grid vely
! points. This is preferable to averaging to the tracer grid in that
! averaging of the land/sea mask and near the pole of the atmosphere
! grid is less severe.
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine set_lsmu ( ncols, nrowst, nrowsuv, LCyclic,            &
     &                      lambda_t, phi_t, lsmt, ILandPt,             &
     &                      lsmuv, lambda_uv, phi_uv )

      implicit none

! declaration of arguments

      integer ncols                   ! IN # of columns (east-west)
      integer nrowst                  ! IN # of rows on tracer grid
      integer nrowsuv                 ! IN no of rows on grid with
                                      !      coincident velocity points
      logical LCyclic                 ! IN T => atmosphere grid cyclic
      real lambda_t(ncols)            ! IN atmos tracer longitudes
      real phi_t(nrowst)              ! IN atmos tracer latitudes
      integer lsmt(ncols, nrowst)     ! IN land/sea mask for tracers
      integer ILandPt                 ! IN value of land points

      integer lsmuv(ncols, nrowsuv)  ! OUT land/sea mask for velocities
      real lambda_uv(ncols)          ! OUT grid coordinates (east-west)
      real phi_uv(nrowsuv)           ! OUT grid coords (north-south)

! declaration of local arrays
      integer lsm_full ( ncols+1, nrowst) ! includes "wrap points"

! declaration of local scalars
      integer jrow, icol  ! loop indices for rows and columns

! declaration of externals
      external copy_to_real

!----------------------------------------------------------------------

! 1. Build an extended tracer land / sea mask

      do jrow = 1, nrowst
        do icol = 1, ncols
          lsm_full(icol,jrow) = lsmt(icol,jrow)
        end do
      end do

      if ( LCyclic ) then
        do jrow = 1, nrowst
          lsm_full(ncols+1,jrow) = lsm_full(1,jrow)
        end do
      else
! points lie outside limited model's area
        do jrow = 1, nrowst
          lsm_full(ncols+1,jrow) = ILandPt
        end do
      end if  ! LCyclic

! 2. Convert tracer land/sea mask to velocity grid land/sea mask

      do jrow = 1, nrowsuv
        do icol = 1, ncols
          lsmuv(icol,jrow) = max ( lsm_full(icol+1,jrow+1),             &
     &                             lsm_full(icol  ,jrow+1),             &
     &                             lsm_full(icol+1,jrow  ),             &
     &                             lsm_full(icol  ,jrow  )  )
        end do
      end do

! 3. Set  latitudes (Phi_uv)
      do jrow = 1, nrowsuv
        phi_uv(jrow) = 0.5 * ( phi_t(jrow) + phi_t(jrow+1) )
      end do

! 4. Set longitudes (Lambda_uv)

      do icol = 1, ncols-1
        lambda_uv(icol) =  0.5 * ( lambda_t(icol) + lambda_t(icol+1) )
      end do

      lambda_uv(ncols) = 2 * lambda_uv(ncols-1) - lambda_uv(ncols-2)

      return
      END SUBROUTINE set_lsmu
!----------------------------------------------------------------------
