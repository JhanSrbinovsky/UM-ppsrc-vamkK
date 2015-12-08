! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Append sub-step number to diagnostic's short name

FUNCTION add_substep_to_sname (SCMop, I)

! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE

! Description:
!   Diagnostics may be defined on multiple sub-steps, in which
!   case the different entries for the different sub-steps need to
!   have the sub-step number appended to their short name before
!   being output so that they can be distinguished.  This routine
!   returns the short name of diagnostic entry with, if it's
!   defined within a sub-stepped part of the code, its sub-step
!   number appended to the end.

! Method:

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
      ! Language: Fortran 90

  ! InOut: The derived-type structure containing
  !        all the diagnostic information
  TYPE(SCMop_type) :: SCMop

  ! In: The entry in SCMop we want the (possibly altered) short name of.
  INTEGER :: I

  ! The output string should be somewhat longer than a normal sname
  ! to account for the appended substep number.
  CHARACTER (len=lsname+10) :: Add_Substep_To_Sname
  CHARACTER (len=8) :: c_substep

  ! Add the sub-step number onto the diagnostic entry's short name
  ! if its sub-step number is greater than one, or it is equal to
  ! one but it is expected there will be multiple sub-steps
  IF (SCMop%substep(I) > 1 .OR.                                         &
     (SCMop%substep(I) == 1 .AND. SCMop%num_substeps > 1)) THEN
    WRITE(c_substep,'(I8)') SCMop%substep(I)
    WRITE(Add_Substep_To_Sname,'(A)')                                   &
      TRIM(SCMop%sname(I))//'_'//TRIM(ADJUSTL(c_substep))
  ELSE

    ! This diagnostic is defined on only one substep, or is outside
    ! of the sub-stepped part(s) of the code. We'll leave its name
    ! as it is then.
    Add_Substep_To_Sname = SCMop%sname(I)
  END IF

  RETURN
END FUNCTION add_substep_to_sname

