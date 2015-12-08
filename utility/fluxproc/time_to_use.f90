! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: time_to_use
!
! Purpose: Flux processing routine.
!          Determines whether to use data time or validity time when
!          searching in the lookup table
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine time_to_use ( itemvalue, l_climate_field, l_data_time)

      USE cfdcodes_mod, ONLY: stcwindstressu, stcwindstressv,           &
                              stcwindmixeng, stcsw, stcsw1, stclongwave,&
                              stcsensibleheat, stcsublim, stctopmelt,   &
                              stcbotmelt, stcevaporation, stcdrain,     &
                              stcconvrain, stcdsnow, stcconvsnow,       &
                              stcsst, stcsss, stchice, stcaice, stcssp, &
                              stcwindspeedu, stcwindspeedv, outstctaux, &
                              outstctauy, outstcwme, outstcsol,         &
                              outstchtn, outstcple, outstcsno,          &
                              outstcsub, outstctop, outstcbot,          &
                              outstcsst, outstcsss, outstchice,         &
                              outstcssp, outstcwspx, outstcwspy, fftaux,&
                              fftauy, ffwme, ffsol, ffhtn, ffple, ffsno,&
                              ffsub, fftop, ffbot, ffsst, ffsss, ffhice,&
                              ffssp, ffwspx, ffwspy, pptaux, pptauy,    &
                              ppwme, ppsol, pphtn, ppple, ppsno, ppsub, &
                              pptop, ppbot, ppsst, ppsss, pphice, ppssp,&
                              ppwspx, ppwspy
      implicit none

! declaration of argument list
      integer itemvalue       ! IN    (stash) item being looked for
      logical l_climate_field ! IN   T => trying to read climate field
                              !      F => trying to read NWP field
      logical l_data_time     ! OUT  T => use data time
                              !      F => use validity time
! declarations of parameters
! no globals, local arrays or scalars

!----------------------------------------------------------------------

! 1. Set l_data_time

      if (   itemvalue  ==  StCAICE .or.                                &
     &       itemvalue  ==  StCSST  .or.                                &
     &       itemvalue  ==  StCHICE .or.                                &
     &       itemvalue  ==  StCSSS  .or.                                &
     &       itemvalue  ==  StCSSP  .or.                                &
     &       itemvalue  ==  StCWindSpeedU  .or.                         &
     &       itemvalue  ==  StCWindSpeedV  .or.                         &
     &       l_climate_field   ) then
        l_data_time = .false.

      else
        l_data_time = .true.

      end if

      return
      END SUBROUTINE time_to_use
!----------------------------------------------------------------------
