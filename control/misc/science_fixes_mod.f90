! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc.
!
! This module declares 'short-term' temporary logicals used to protect
! science bug fixes that lead to significant alterations in science results.
! It is expected that these logicals will be short lived as the preference
! should be for all configurations to use the corrected code. But
! to maintain short term reprocucibility of results across UM versions
! the fixes are protected by logicals until the fixes become the default
! in all model configurations and the logical retired.
!
! All logicals below should have a review period attached to them for future
! retirement of both the logical and the broken code.
!
! ! ticket #xxxx
! LOGICAL :: fix_me = .FALSE.   ! review again MM YYYY
! 
! then add logical to namelist /temp_fixes/
! and add code to subroutine warn_temp_fixes to report when fix is not used
! ie when the logical = .FALSE.
!
! -----------------------------------------------------------
! -----------------------------------------------------------

MODULE science_fixes_mod

IMPLICIT NONE

! ticket #5493
LOGICAL :: l_use_old_mass_fix = .FALSE. ! Review Again May 2014
! ticket #4800
LOGICAL :: L_p2t_weight_fix = .FALSE.   ! Review Again June 2014  
! ticket #4982
LOGICAL :: l_rm_neg_par     = .FALSE.   ! Review Again June 2014
! ticket #4612
LOGICAL :: l_error_ancil_struct = .FALSE. ! Review Again June 2013
! ticket #5168
LOGICAL :: l_exclude_fix_theta_source = .FALSE.! Review Again May 2014
! ticket #5414
LOGICAL :: l_riverinland_fix = .FALSE.  ! Review Again June 2014
! ticket #5551
LOGICAL :: l_emis_ssi_full = .FALSE.   ! Review again June 2014

NAMELIST/temp_fixes/l_use_old_mass_fix,                       &
        L_p2t_weight_fix, l_rm_neg_par, l_error_ancil_struct, &
        l_exclude_fix_theta_source, l_riverinland_fix,        &
        l_emis_ssi_full
        
! -----------------------------------------------------------        
! -----------------------------------------------------------
        
CONTAINS

SUBROUTINE warn_temp_fixes()

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY : ereport

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER :: ErrorStatus            ! Return code : 0 Normal Exit : >0 Error
CHARACTER(LEN=256) :: CMessage    ! Error message if Errorstatus /=0

IF (lhook) CALL dr_hook('warn_temp_fixes',zhook_in,zhook_handle)

ErrorStatus = 0


! -----------------------------------------------------------        
! -----------------------------------------------------------
! define whether the fix is appropriate to RCF, UM RUN or both
! and warn the user if the fix is not used.


 
IF (.NOT. L_p2t_weight_fix) THEN
  ErrorStatus = -100 
  CMessage    = 'Model run excludes ticket #4800 as'// &
                ' L_p2t_weight_fix=.FALSE.' 
                          
  CALL ereport('warn_temp_fixes', ErrorStatus, CMessage) 
END IF
IF (.NOT. l_rm_neg_par) THEN
  ErrorStatus = -100 
  CMessage    = 'Model run excludes ticket #4982 as'// &
                ' l_rm_neg_par=.FALSE.' 
                          
  CALL ereport('warn_temp_fixes', ErrorStatus, CMessage) 
END IF
IF (.NOT. l_error_ancil_struct) THEN
  ErrorStatus = -100 
  CMessage    = 'Model run excludes ticket #4612 as'// &
                ' L_error_ancil_struct=.FALSE.' 
                          
  CALL ereport('warn_temp_fixes', ErrorStatus, CMessage) 
END IF
IF (l_exclude_fix_theta_source) THEN
  ErrorStatus = -100 
  CMessage    = 'Model run excludes ticket #5168 as'// &
                ' l_exclude_fix_theta_source=.TRUE.'// &
                ' your results will be WRONG'
                          
  CALL ereport('warn_temp_fixes', ErrorStatus, CMessage) 
END IF
IF (.NOT. l_riverinland_fix) THEN
  ErrorStatus = -100 
  CMessage    = 'Model run excludes ticket #5414 as'// &
                ' l_riverinland_fix=.FALSE.'
  
  CALL ereport('warn_temp_fixes', ErrorStatus, CMessage) 
END IF
IF (.NOT.l_emis_ssi_full) THEN
  ErrorStatus = -100 
  CMessage    = 'Model run excludes ticket #5551 as'// &
                ' l_emis_ssi_full=.FALSE.'// &
                ' your results will be WRONG'// &
                ' if the emissivity of sea-ice is not 1.'
                          
  CALL ereport('warn_temp_fixes', ErrorStatus, CMessage) 
END IF

IF (l_use_old_mass_fix) THEN
  ErrorStatus = -100
  CMessage    = 'Model run excludes ticket #5493 as'// &
                ' l_use_old_mass_fix=.TRUE.'

  CALL ereport('warn_temp_fixes', ErrorStatus, CMessage)
END IF




! -----------------------------------------------------------        
! -----------------------------------------------------------

        
IF (lhook) CALL dr_hook('warn_temp_fixes',zhook_out,zhook_handle)
RETURN        
END SUBROUTINE warn_temp_fixes

END MODULE science_fixes_mod

