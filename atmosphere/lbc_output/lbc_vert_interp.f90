! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Vertical Interpolation of LBC data
!
! Subroutine Interface:

      SUBROUTINE lbc_vert_interp (                                      &
                 lbc_vi_data_in,                                        &
                 lbc_vi_orog,                                           &
                 lbc_vi_data_out,                                       &
                 lbc_seg_size,                                          &
                 level_type,                                            &
! src
                 src_model_levels,                                      &
                 src_levels,                                            &
                 src_first_level,                                       &
                 src_last_level,                                        &
                 src_ht_gen_method,                                     &
                 src_first_r_rho,                                       &
                 src_z_top_model,                                       &
                 src_eta_theta,                                         &
                 src_eta_rho,                                           &
! lbc
                 interp_order,                                          &
                 lbc_model_levels,                                      &
                 lbc_levels,                                            &
                 lbc_first_level,                                       &
                 lbc_last_level,                                        &
                 lbc_ht_gen_method,                                     &
                 lbc_first_r_rho,                                       &
                 lbc_z_top_model,                                       &
                 lbc_eta_theta,                                         &
                 lbc_eta_rho                                            &
       )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE vert_interp_mod, ONLY: vert_interp
      IMPLICIT NONE

!
! Description:
!   Performs Vertical interpolation of Lateral Boundary Conditions (LBC)
!   from model to LBC levels.
!
! Method:
!   1. Heights of rho & theta levels are calculated for both model and
!      LBC levels at the LBC points.
!   2. VERT_INTERP is called to interpolate the LBCs. Linear
!      interpolation used for all variables.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

      INTEGER  ::  lbc_seg_size
      INTEGER  ::  level_type     !  Theta or Rho levels

! source vertical grid
      INTEGER  :: src_model_levels
      INTEGER  :: src_levels
      INTEGER  :: src_first_level
      INTEGER  :: src_last_level
      INTEGER  :: src_ht_gen_method
      INTEGER  :: src_first_r_rho
      REAL     :: src_z_top_model
      REAL     :: src_eta_theta (0:src_model_levels)
      REAL     :: src_eta_rho   (  src_model_levels)
! lbc vertical grid
      INTEGER  :: interp_order
      INTEGER  :: lbc_model_levels
      INTEGER  :: lbc_levels
      INTEGER  :: lbc_first_level
      INTEGER  :: lbc_last_level
      INTEGER  :: lbc_ht_gen_method
      INTEGER  :: lbc_first_r_rho
      REAL     :: lbc_z_top_model
      REAL     :: lbc_eta_theta (0:lbc_model_levels)
      REAL     :: lbc_eta_rho   (  lbc_model_levels)

! data
      REAL  ::  lbc_vi_data_in  (lbc_seg_size,                          &
                                 src_first_level:src_last_level)
      REAL  ::  lbc_vi_data_out (lbc_seg_size,                          &
                                 lbc_first_level:lbc_last_level)
      REAL  ::  lbc_vi_orog     (lbc_seg_size)

! Local variables
      INTEGER, PARAMETER :: rho_levels   = 1
      INTEGER, PARAMETER :: theta_levels = 2

      CHARACTER (LEN=*), PARAMETER ::  routinename = 'LBC_Vert_Interp'

      INTEGER :: level
      INTEGER :: errorstatus
      CHARACTER (LEN=80) :: cmessage

! Height fields to be calculated
      REAL, ALLOCATABLE :: src_r_theta_levels(:,:)
      REAL, ALLOCATABLE :: src_r_rho_levels  (:,:)
      REAL, ALLOCATABLE :: lbc_r_theta_levels(:,:)
      REAL, ALLOCATABLE :: lbc_r_rho_levels  (:,:)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------
! Check validity of height algorithm
! ----------------------------------

      IF (lhook) CALL dr_hook('LBC_VERT_INTERP',zhook_in,zhook_handle)

!     Height algorithm for model and lbcs must be the same.
      If (src_ht_gen_method /= lbc_ht_gen_method) Then
          Write (CMessage,*) 'Mismatch in height algorithm for ',       &
     &    'model ',src_ht_gen_method,' and LBCs ',lbc_ht_gen_method
          ErrorStatus = 10

          Call Ereport (RoutineName, ErrorStatus, CMessage)
      End If

! ------------------------------------
! Allocate space for the height fields
! ------------------------------------

      allocate (src_r_theta_levels(lbc_seg_size,0:src_model_levels)  )
      allocate (src_r_rho_levels  (lbc_seg_size,  src_model_levels+1))

      allocate (lbc_r_theta_levels(lbc_seg_size,0:lbc_model_levels)  )
      allocate (lbc_r_rho_levels  (lbc_seg_size,  lbc_model_levels+1))

! ---------------------------------------
! Calculate heights for the source levels
! ---------------------------------------

! DEPENDS ON: lbc_calc_heights
      Call lbc_calc_heights ( lbc_seg_size,                             &
     &                        src_model_levels, src_ht_gen_method,      &
     &                        src_first_r_rho, src_z_top_model,         &
     &                        src_eta_theta, src_eta_rho,               &
     &                        lbc_vi_orog,                              &
     &                        src_r_theta_levels, src_r_rho_levels )

! ------------------------------------
! Calculate heights for the lbc levels
! ------------------------------------

! DEPENDS ON: lbc_calc_heights
      Call lbc_calc_heights ( lbc_seg_size,                             &
     &                        lbc_model_levels, lbc_ht_gen_method,      &
     &                        lbc_first_r_rho, lbc_z_top_model,         &
     &                        lbc_eta_theta, lbc_eta_rho,               &
     &                        lbc_vi_orog,                              &
     &                        lbc_r_theta_levels, lbc_r_rho_levels )

! -----------------------------
! Do the Vertical Interpolation
! -----------------------------

      Select Case (Level_Type)

        Case (Rho_Levels)

          Do level = lbc_first_level, lbc_last_level

!           write (6,*) ' VI on Rho levels for level ',level

            Call Vert_Interp (lbc_vi_data_in, lbc_seg_size,             &
     &                        src_levels, lbc_r_rho_levels(1:,level),   &
     &                        src_r_rho_levels(1:,src_first_level:),    &
     &                        interp_order,                             &
     &                        lbc_vi_data_out(:,level) )

          End Do

        Case (Theta_Levels)

          Do level = lbc_first_level, lbc_last_level

!           write (6,*) ' VI on Theta levels for level ',level

            Call Vert_Interp (lbc_vi_data_in, lbc_seg_size,             &
     &                        src_levels, lbc_r_theta_levels(1:,level), &
     &                        src_r_theta_levels(1:,src_first_level:),  &
     &                        interp_order,                             &
     &                        lbc_vi_data_out(:,level) )

          End Do

        Case Default

          Write (CMessage,*) 'LBC Vertical Interpolation not ',         &
     &                       'catered for Level Type ',Level_Type
          ErrorStatus = 20

          Call Ereport (RoutineName, ErrorStatus, CMessage)

      End Select

! --------------------
! Deallocate workspace
! --------------------

      deallocate ( src_r_theta_levels )
      deallocate ( src_r_rho_levels   )
      deallocate ( lbc_r_theta_levels )
      deallocate ( lbc_r_rho_levels   )

      IF (lhook) CALL dr_hook('LBC_VERT_INTERP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE lbc_vert_interp
