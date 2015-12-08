! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details pleaaase refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Purpose:
!   Dummy routine to resolve calls for ISCCP routine

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!   Language: FORTRAN 90.
! -----------------------------ATTENTION--------------------------------
! The full isccp_cloud_types subroutine is available under the GNU
! Lesser General Public License (http://www.gnu.org) and is not
! currently distributed with the UM.
! If you require this file for use with the UM please contact the code
! owner at: mark.webb@metoffice.gov.uk
! ----------------------------------------------------------------------

SUBROUTINE isccp_cloud_types()


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

INTEGER                      ::  icode
CHARACTER (LEN=80)           ::  cmessage
CHARACTER (LEN=*), PARAMETER ::  RoutineName='isccp_cloud_types'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('ISCCP_CLOUD_TYPES',zhook_in,zhook_handle)
cmessage = 'Routine should not be callable'
icode = 1

CALL ereport(RoutineName,icode,cmessage)

IF (lhook) CALL dr_hook('ISCCP_CLOUD_TYPES',zhook_out,zhook_handle)
RETURN
END SUBROUTINE isccp_cloud_types
