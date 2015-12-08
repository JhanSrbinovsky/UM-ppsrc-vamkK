! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Init_turb_diff
!

      SUBROUTINE Init_turb_diff(                                        &
                                 model_levels, bl_levels,               &
                                 start_hor_lev, end_hor_lev,            &
                                 start_vert_lev, end_vert_lev,          &
                                 global_u_filter, global_v_filter,      &
                                 diff_factor_in, mix_factor_in,         &
                                 L_horiz, L_vert, L_blend )

! Purpose:
!          Initialises turbulent diffusion parameters
!          and declares arrays. Called from setcona
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

USE proc_info_mod, ONLY: model_domain
USE turb_diff_mod, ONLY: l_subfilter_horiz, l_subfilter_vert,         &
     l_subfilter_blend, turb_startlev_horiz,turb_endlev_horiz,        &
     turb_startlev_vert,turb_endlev_vert, diff_factor, mix_factor
USE UM_ParParams
USE ereport_mod, ONLY : ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE bl_option_mod, ONLY: Keep_Ri_FA, except_disc_inv, Kprof_cu,       &
                         klcl_entr

      IMPLICIT NONE

! Input variables
      LOGICAL  :: L_horiz
      LOGICAL  :: L_vert
      LOGICAL  :: L_blend

      INTEGER  :: model_levels
      INTEGER  :: bl_levels
      INTEGER  :: start_hor_lev
      INTEGER  :: end_hor_lev
      INTEGER  :: start_vert_lev
      INTEGER  :: end_vert_lev
      INTEGER  :: global_u_filter
      INTEGER  :: global_v_filter

      REAL  :: diff_factor_in
      REAL  :: mix_factor_in

! Loop indices.
      INTEGER  ::  i, j
      INTEGER  ::  errorstatus

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('INIT_TURB_DIFF',zhook_in,zhook_handle)

! Code to prepare for Smagorinsky turbulence scheme
      L_subfilter_horiz = L_horiz
      L_subfilter_vert  = L_vert
      L_subfilter_blend = L_blend

      turb_startlev_horiz = start_hor_lev
      turb_endlev_horiz   = end_hor_lev
      turb_startlev_vert  = start_vert_lev
      turb_endlev_vert    = end_vert_lev

      diff_factor = diff_factor_in
      mix_factor  = mix_factor_in

      WRITE(6,*) ' '
      IF (bl_levels /= model_levels - 1) THEN
          WRITE(6,*)'BL_LEVELS = ',bl_levels,                           &
                    'MODEL_LEVELS = ',model_levels
          errorstatus=123
        CALL ereport("INIT_TURB_DIFF", errorstatus,                     &
              "The number of boundary layer levels must equal " //      &
              "the number of model levels minus one when the " //       &
              "3D subgrid turbulence scheme is selected.")
      ELSE ! bl_levels = model_levels - 1
        WRITE(6,*)'3D subgrid turbulence scheme is active '
        WRITE(6,*)' bl_levels = ', bl_levels
      END IF ! bl_levels /= model_levels - 1

      IF (diff_factor < 0.0 .OR. diff_factor > 1.0) THEN
        WRITE(6,*)' diff_factor = ', diff_factor
        errorstatus=123
       CALL ereport("INIT_TURB_DIFF", errorstatus,                      &
              "diff_factor must have a value greater than " //          &
              "zero and less than or equal to 1.0 so that the "//       &
              "numerical stability is maintained")
      ELSE    ! 0.0 <= diff_factor <= 1.0 
        WRITE(6,*)' diff_factor = ', diff_factor
      END IF !  diff_factor < 0.0 .OR. diff_factor > 1.0

      IF (mix_factor <= 0.0 .OR. mix_factor > 1.0) THEN
        WRITE(6,*)' mix_factor =',mix_factor
        errorstatus=123
        CALL ereport("INIT_TURB_DIFF", errorstatus,                     &
              "mix_factor should have a value greater or equal " //     &
              "to zero and less than or equal to 1.0 ")
      ELSE    ! 0.0 <= mix_factor <= 1.0 
        WRITE(6,*)' mix_factor =', mix_factor
      END IF !  mix_factor <= 0.0 .OR. mix_factor > 1.0

      IF (l_subfilter_horiz) THEN

        IF (turb_startlev_horiz > turb_endlev_horiz) THEN
            errorstatus=123
         CALL ereport("INIT_TURB_DIFF", errorstatus,                    &
                "The start level for the horizontal turbulence  " //    &
                "scheme is greater than the end level! ")
        END IF ! turb_startlev_horiz > turb_endlev_horiz

! The levels over which the turbulence scheme acts must be
! between 2 and model_levels-1

        IF (turb_startlev_horiz < 2) turb_startlev_horiz = 2
        IF (turb_endlev_horiz > model_levels - 1)                       &
            turb_endlev_horiz = model_levels - 1

        WRITE(6,*)' Horizontal subgrid turbulence'
        WRITE(6,*)'  turb_startlev_horiz = ', turb_startlev_horiz
        WRITE(6,*)'  turb_endlev_horiz = ', turb_endlev_horiz
        IF (model_domain == mt_global) THEN
          WRITE(6,'(A)') ' Horizontal subgrid turbulence will not' //   &
                         ' be applied in longitudinal'
          WRITE(6,'(A)') ' direction where polar filter is applied,'
          WRITE(6,'(A, I3, A, I3, A)') ' i.e. for first and last',      &
                            global_u_filter+1, ' rows for u/theta and', &
                                       global_v_filter, ' rows for v.'
        END IF !  model_domain == mt_global

      ELSE  !  l_subfilter_horiz
 
        WRITE(6,*)' Horizontal subgrid turbulence is not active '
 
      END IF !  l_subfilter_horiz

      IF (l_subfilter_vert) THEN

        IF (turb_startlev_vert > turb_endlev_vert) THEN
            errorstatus=123
          CALL ereport("INIT_TURB_DIFF", errorstatus,                   &
                 "The start level for the vertical turbulence  " //     &
                 "scheme is greater than the end level! ")
        END IF !  turb_startlev_vert > turb_endlev_vert

! The levels over which the turbulence scheme acts must be
! between 2 and model_levels-1

        IF ( turb_startlev_vert < 2) turb_startlev_vert = 2
        IF ( turb_endlev_vert > model_levels - 1)                       &
             turb_endlev_vert = model_levels - 1

        WRITE(6,*)' Vertical subgrid turbulence '
        WRITE(6,*)'  turb_startlev_vert = ', turb_startlev_vert
        WRITE(6,*)'  turb_endlev_vert = ', turb_endlev_vert

      ELSE  !  l_subfilter_vert
 
        WRITE(6,*)' Vertical subgrid turbulence is not active '
 
      END IF  !  l_subfilter_vert

      IF (l_subfilter_blend) THEN

        ! Set BL options appropriately for blending
        Keep_Ri_FA = except_disc_inv
        Kprof_cu   = klcl_entr

      END IF

IF (lhook) CALL dr_hook('INIT_TURB_DIFF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE init_turb_diff
