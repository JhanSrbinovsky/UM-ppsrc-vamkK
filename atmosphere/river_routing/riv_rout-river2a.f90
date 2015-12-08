! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE riv_rout_mod_2A

IMPLICIT NONE

CONTAINS



      SUBROUTINE RIV_ROUT_2A(                                           &
     &   SURF_RUNOFF,SUB_SURF_RUNOFF,ICOLS,JROWS,                       &
     &   DELTA_PHI,FIRST,G2G_TIMESTEP,                                  &
! ancillary variables
     &   IAREA, SLOPE, FLOWOBS1,INEXT,JNEXT,LAND,                       &
! prognostic variables
     &   SUBSTORE,SURFSTORE,FLOWIN,BFLOWIN,RIVFLOW)


!
! Purpose:
!
! Perform the routing of surface and sub-surface runoff for Regional
!  model.
!
! Method:
! This subroutine routes surface and subsurface runoff
! for the RCM. The routing procedure is performed on the
! RCM grid, so no regridding is needed. The routing
! procedure is currently based on the Kinematic Wave
! routing model used at CEH Wallingford.
!
!-----------------------------------------------------------------

      USE earth_constants_mod, ONLY: earth_radius

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE c_grid2grid_mod, ONLY: cland, criver, cbland, cbriver,        &
                                 runoff_factor, retl, retr, slfac, a_thresh
      IMPLICIT NONE

! IN Arguments

! External variables
!=====================

      INTEGER                                                           &
     & icols                                                            &
                  ! IN NO. OF COLUMNS IN ATMOSPHERE (icols)
     &,jrows      ! IN NO. OF ROWS IN ATMOSPHERE (TP GRID) : jrows
      REAL                                                              &
     & SURF_RUNOFF(icols,jrows)                                         &
        ! IN TOTAL RATE OF surface RUNOFF (KG/M2/S=mm/s)
     &,SUB_SURF_RUNOFF(icols,jrows)                                     &
        ! IN TOTAL RATE OF  sub-surf RUNOFF (KG/M2/S=mm/s)
     &,DELTA_PHI                                                        &
        ! RCM gridsize (radians)
     &,G2G_TIMESTEP   ! IN river timestep (secs)

      REAL                                                              &
     & RIVFLOW(icols,jrows)
! OUT river flow out from each gridbox (KG/m2/S = mm/s)

      logical FIRST  ! first call to river routing ?(T/F)
! ancillary variables

       integer                                                          &
     &  iarea(icols,jrows)                                              &
                            !accumulated areas file
     &, inext(icols,jrows)                                              &
                            !x-coordinate of downstream grid point
     &, jnext(icols,jrows)                                              &
                            !y-coordinate of downstream grid point
     &, land(icols,jrows)   !land/river depends on value a_thresh ?

       real                                                             &
     &  slope(icols,jrows)                                              &
                                !slopes (not used yet)
     &  ,flowobs1(icols,jrows)  !optional initialisation for flows


! prognostic variables

       real                                                             &
     &  substore(icols,jrows)                                           &
                                !routing sub_surface store (mm)
     & ,surfstore(icols,jrows)                                          &
                                !routing surface store (mm)
     & ,flowin(icols,jrows)                                             &
                                !surface lateral inflow (mm)
     & ,bflowin(icols,jrows)    !sub-surface lateral inflow (mm)

! internal variables
!====================
       integer                                                          &
     & landtype                                                         &
                        !local for land type
     &,in,jn                                                            &
                        !local co-ords of downstream point
     &,i,j              !co-ordinate counters in do loops


       real                                                             &
     &  landtheta,rivertheta                                            &
           !parameters- surface wave speeds
     &, sublandtheta, subrivertheta                                     &
           !parameters- sub-surface wave speeds
!    &, retl,retr
!          !parameters - return flow for land and river
     &, substore_n(icols,jrows)                                         &
           !internal array for surface store at next timestep(mm)
     &, surfstore_n(icols,jrows)                                        &
           !internal array for subsurface store at next timestep(mm)
     &, flowin_n(icols,jrows)                                           &
           !internal, surface lateral inflow at next timestep(mm)
     &, bflowin_n(icols,jrows)                                          &
           !internal, sub-surface lateral inflow at next timestep(mm)
     &, returnflow                                                      &
           !internal variable for returnflow
!    &, slfac
!          !currently unused slope variables
     &, real_area(icols,jrows)                                          &
           !internal, real value of accumulated area (iarea)
     &, baseflow(icols,jrows)                                           &
           !sub-surface flow - currently internal
     &, dt                                                              &
           !internal value of g2g model t/s
     &, surf_roff(icols,jrows)                                          &
           !internal value of surf_runoff
     &, sub_surf_roff(icols,jrows)                                      &
           !internal value of sub_surf_runoff
     &, dx
           !RCM grid spacing in metres

       INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
       INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
       REAL(KIND=jprb)               :: zhook_handle



! include file to setup grid-to-grid model parameters

       IF (lhook) CALL dr_hook('RIV_ROUT_2A',zhook_in,zhook_handle)
        dt = g2g_timestep            ! grid-to-grid model timestep(s)
        dx = earth_radius*delta_phi  ! horizontal gridsize (m)

! Set up routing parameters

         landtheta = cland*dt /dx
         rivertheta = criver*dt /dx
         sublandtheta = cbland*dt /dx
         subrivertheta = cbriver*dt /dx

         if (landtheta >  1. .or. rivertheta >  1.or.                   &
     &      sublandtheta >  1. .or. subrivertheta >  1.) then
             write(6,*)'Unstable theta: reduce wavespeed in routing'
             write(6,*) 'The theta values will be set to zero'
            landtheta = 0.
            rivertheta = 0.
            sublandtheta = 0.
            subrivertheta = 0.
         end if

! Initialise variables at first timestep
!-----------------------------------------

        if (FIRST) then

! From the cumulative catchment areas dataset, determine which
! squares are land or river
!------------------------------------------------------------------
! land(i,j) is set at the first entry, and the values must be
! remembered for subsequent timesteps
!-------------------------------------------------------------------
        do  j = 1, jrows
           do i = 1, icols
             if (iarea(i,j) <  0) then
                land(i,j) = 0          !Sea
             else if (iarea(i,j) >  a_thresh) then
                land(i,j) = 1          !river
             else
                land(i,j) = 2          !land
             end if
             real_area(i,j) = real(iarea(i,j)) +1.
             ! include current point, so add 1
           end do
        end do


! initialisation of variables

        do  j = 1, jrows
           do i = 1, icols

! this is where we will eventually initialise the routing stores
! using flow observations; here they are just set to zero....
              surfstore(i,j) = 0.
              substore(i,j) = 0.

! initialise surface and sub-surface inflows
              flowin(i,j) = 0.
              bflowin(i,j) = 0.
              rivflow(i,j) = 0.

           end do
        end do

        end if     ! FIRST
!----------------------------------------------------------------


! For each timestep...
!-----------------------------

! convert runoff from mm/sec to mm (over the model timestep)
! assume dx in metres and dt in seconds

       do j = 1,jrows
          do i = 1,icols
            if (surf_runoff(i,j) >= 0.and.sub_surf_runoff(i,j)          &
     &               >= 0.) then
               surf_roff(i,j)    =runoff_factor*surf_runoff(i,j)*dt
               sub_surf_roff(i,j)=runoff_factor*sub_surf_runoff(i,j)*dt
            else
               surf_roff(i,j) = 0.       ! weed out dodgy data
               sub_surf_roff(i,j) = 0.
             end if
          end do
       end do


! Initialise accumulated inflows and stores for the next timestep

        do  j = 1, jrows
           do i = 1, icols
              flowin_n(i,j) = 0.
              bflowin_n(i,j) = 0.
              surfstore_n(i,j) = 0.
              substore_n(i,j)= 0.
           end do
        end do

! Route runoff


      do j = 1,jrows
         do i = 1,icols

                in = inext(i,j)
                jn = jnext(i,j)
!               sl = slope(i,j)       ! not used at the moment
                landtype = land(i,j)

! route runoff using simple kw model

                if (landtype == 2) then  !land
!-----------------------------------------------------------------
! land surface
                   surfstore_n(i,j) = (1.-landtheta )*surfstore(i,j)+   &
     &                flowin(i,j) +  surf_roff(i,j)
                   flowin_n(in,jn) =flowin_n(in,jn) +                   &
     &                   landtheta*surfstore_n(i,j)

! land subsurface
                   substore_n(i,j) = (1.-sublandtheta )*substore(i,j)   &
     &                + bflowin(i,j) + sub_surf_roff(i,j)
                   bflowin_n(in,jn) = bflowin_n(in,jn)                  &
     &                + sublandtheta*substore_n(i,j)

! return flow
                  returnflow = amax1(abs(substore_n(i,j)*retl),0.)
                  substore_n(i,j) = substore_n(i,j) - returnflow
                  surfstore_n(i,j) = surfstore_n(i,j) +returnflow


                else if (landtype == 1) then  !river
!-------------------------------------------------------------------

! river subsurface
                  substore_n(i,j) = (1.-subrivertheta )*substore(i,j)+  &
     &                bflowin(i,j) + sub_surf_roff(i,j)
                  bflowin_n(in,jn) = bflowin_n(in,jn)                   &
     &                 + subrivertheta*substore_n(i,j)


! river surface
                  surfstore_n(i,j) = (1.-rivertheta )*surfstore(i,j)+   &
     &                flowin(i,j) + surf_roff(i,j)
                  flowin_n(in,jn) = flowin_n(in,jn)                     &
     &                + rivertheta*surfstore_n(i,j)

! return flow
                  returnflow = amax1(abs(substore_n(i,j)*retr),0.)
                  substore_n(i,j) = substore_n(i,j) - returnflow
                  surfstore_n(i,j) = surfstore_n(i,j) +returnflow

          end if
!------------------------------------------------------------------
            end do
         end do   ! end of routing loop

! housekeeping for next timestep

        do  j = 1, jrows
          do i = 1, icols

! keep inflows for next timestep

              flowin(i,j) = flowin_n(i,j)
              bflowin(i,j) = bflowin_n(i,j)

! calculate flow in rivers (mm/s = kg/m2/s)

              if (land(i,j) == 1) then
                  rivflow(i,j) = rivertheta/dt*surfstore(i,j)
                  baseflow(i,j) = subrivertheta/dt*substore(i,j)

                  in = inext(i,j)
                  jn = jnext(i,j)

! calculate flow into the sea if the next point is sea
!(assume it's equal to river flow in the adjacent land pt)

                  if (land(in,jn) == 0)                                 &
     &               rivflow(in,jn) = rivflow(i,j)

              end if

! keep routing stores for next timestep
              surfstore(i,j) = surfstore_n(i,j)
              substore(i,j) = substore_n(i,j)
           end do
        end do



        IF (lhook) CALL dr_hook('RIV_ROUT_2A',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE RIV_ROUT_2A

END MODULE riv_rout_mod_2A
