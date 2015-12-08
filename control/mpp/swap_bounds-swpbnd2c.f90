
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose:
!   This subroutine takes care of all boundary swapping and
!   extending of arrays at the global boundaries. It is a wrapper, 
!   used when 'USE halo_exchange' is absent from the calling routine.
!   Core functionality now in halo_exchange.F90
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: MPP

SUBROUTINE swap_bounds(                &
    field, row_length, rows, levels,   & ! field
    halo_x, halo_y,                    & ! halos
    field_type, l_vector               & ! supporting information
    )

  USE mpl, ONLY :                      &
      mpl_real,                        &
      mpl_status_size
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE Halo_Exchange, ONLY : &
      Swap_bounds_EW,       &
      Swap_Bounds_NS,       &
      he_initialised
  USE ereport_mod,   ONLY : ereport
  USE um_parvars,    ONLY : sb_model_domain
  USE domain_params, ONLY : mt_lam

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) ::  row_length ! number of points on a row
                                     ! (not including halos)
  INTEGER, INTENT(IN) ::  rows       ! number of rows in a theta field
                                     ! (not including halos)
  INTEGER, INTENT(IN) ::  levels     ! number of model levels
  INTEGER, INTENT(IN) ::  halo_x     ! size of halo in "i" direction
  INTEGER, INTENT(IN) ::  halo_y     ! size of halo in "j" direction
  INTEGER, INTENT(IN) ::  field_type ! Defines the grid interpolation type
                                     ! of the input FIELD (u,v or w)
  LOGICAL, INTENT(IN) :: l_vector    ! TRUE:  Data is a horizontal vector
                                     !        component
                                     ! FALSE: Data is a scalar

  REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
      1-halo_y:rows+halo_y,levels)   ! Field to have its halos updated

  INTEGER :: errorstatus
  
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1    
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('SWAP_BOUNDS',zhook_in,zhook_handle)

  errorstatus = 0

  IF (halo_x == 0 .AND. halo_y == 0 .AND. rows == 1 &
                                    .AND. row_length == 1) THEN
      !Choose not to output warning as the SCM output becomes massive
      !errorstatus = -1
      ! CALL ereport('SWAP_BOUNDS',errorstatus,       &
      !            ' not initialised as this is a SCM.')
      HE_Initialised = .TRUE.
  END IF    


  IF (.NOT. HE_Initialised) THEN
      errorstatus = 99
      CALL ereport('SWAP_BOUNDS',errorstatus,' not initialised')
  END IF    


  !------------------------------------------------------------------
  ! 2.0 Exchanges in EW then NS directions
  IF (halo_x > 0)                        &
      
      CALL swap_bounds_EW(               &
      field, row_length, rows, levels,   & ! field
      halo_x, halo_y                     & ! halos
      )
  
  IF (halo_y > 0)                        &
      CALL swap_bounds_NS(               &
      field, row_length, rows, levels,   & ! field
      halo_x, halo_y,                    & ! halos
      field_type, l_vector               & ! supporting information
      )
  
  IF (sb_model_domain == mt_lam) THEN
!DEPENDS ON: fill_external_halos
    CALL fill_external_halos(field, row_length, rows, levels, &
                             halo_x, halo_y)
  END IF

  IF (lhook) CALL dr_hook('SWAP_BOUNDS',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE swap_bounds
