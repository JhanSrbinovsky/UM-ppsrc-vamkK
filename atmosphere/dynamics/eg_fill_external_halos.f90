! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_fill_external_halos_mod

IMPLICIT NONE

CONTAINS

! Purpose:
!   wrapper around fill_external_halos to eliminate grid variables and
!   do some basic error checking
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Dynamics

SUBROUTINE eg_fill_external_halos(field, dims, l_test)

  USE parkind1, ONLY: jpim, jprb       !DrHook
  USE yomhook,  ONLY: lhook, dr_hook   !DrHook

  USE atm_fields_bounds_mod, ONLY: array_dims
  USE integrity_mod,         ONLY: integrity_test
  USE ereport_mod,           ONLY: ereport
  USE proc_info_mod,         ONLY: me,n_proc





  IMPLICIT NONE

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  
  TYPE (array_dims) , INTENT(IN) :: dims

  REAL, INTENT(INOUT) :: field(dims%i_start:dims%i_end, &
                               dims%j_start:dims%j_end, &
                               dims%k_start:dims%k_end)
                                     ! Field to have its halos updated


  LOGICAL, OPTIONAL :: l_test ! fill halo with NaNs to test if swapbound
                              ! was required

  INTEGER ::  row_length ! number of points on a row
                                     ! (not including halos)
  INTEGER ::  rows       ! number of rows in a theta field
                                     ! (not including halos)
  INTEGER ::  levels     ! number of model levels
  INTEGER ::  halo_x     ! size of halo in "i" direction
  INTEGER ::  halo_y     ! size of halo in "j" direction

  INTEGER ::  ierr,istat

  CHARACTER (len=65) :: hash_in,hash_out
  CHARACTER (len=80) :: cmessage

  REAL,    PARAMETER :: inf = O'0777600000000000000000'
  REAL,    PARAMETER :: nan = O'0777610000000000000000'
  INTEGER, PARAMETER :: bigint = O'17777777777'

  IF (lhook) CALL dr_hook('EG_FILL_EXTERNAL_HALOS',zhook_in,zhook_handle)

  halo_x = dims%halo_i
  halo_y = dims%halo_j

  row_length = dims%i_end-dims%i_start+1-2*halo_x
  rows       = dims%j_end-dims%j_start+1-2*halo_y

  levels     = dims%k_end-dims%k_start+1

  ! DEPENDS ON: eg_hash

  IF ( integrity_test ) CALL eg_hash(field,8*SIZE(field),hash_in)

! DEPENDS ON: swap_bounds
  CALL fill_external_halos(            &
    field, row_length, rows, levels,   & ! field
    halo_x, halo_y                     & ! halos
    )

  IF ( integrity_test ) THEN
    CALL eg_hash(field,8*SIZE(field),hash_out)

    ierr = 0

    IF (hash_in == hash_out) ierr = 1

    CALL GC_ISUM(1, n_proc, istat, ierr)

    IF (ierr.eq.n_proc) THEN 
     cmessage = 'potentially unnecessary fill external halo detected '
     IF (me == 0) THEN
       WRITE(0,fmt='(2A,L1,2I8)') cmessage,hash_out,hash_in.eq.hash_out,  &
                                 ierr,n_proc
       ierr = -1
       CALL ereport("eg_fill_external_halos",ierr,cmessage)






      END IF
    END IF
  END IF

  IF (PRESENT(l_test)) THEN

    IF (l_test) THEN

      field(dims%i_start:dims%i_start+dims%halo_i-1,                  &
            dims%j_start:dims%j_end,                                  &
            dims%k_start:dims%k_end)                 = nan

      field(dims%i_end-dims%halo_i+1:dims%i_end,                      &
            dims%j_start:dims%j_end,                                  &
            dims%k_start:dims%k_end)                 = nan

      field(dims%i_start:dims%i_end,                                  &
            dims%j_start:dims%j_start+dims%halo_j-1,                  &
            dims%k_start:dims%k_end)                 = nan

      field(dims%i_start:dims%i_end,                                  &
            dims%j_end-dims%halo_j+1:dims%j_end,                      &
            dims%k_start:dims%k_end)                 = nan

    END IF
  END IF

  IF (lhook) CALL dr_hook('EG_FILL_EXTERNAL_HALOS',zhook_out,zhook_handle)

END SUBROUTINE eg_fill_external_halos

END MODULE eg_fill_external_halos_mod
