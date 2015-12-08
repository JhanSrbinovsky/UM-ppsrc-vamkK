! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     subroutine U_TO_p  for calculating variables held at u points
!     at p points,  on both the new dynamics and v-at-poles/endgame grids
!
!     This routine does interior points of array not halos,
!     but requires halo information to be set.
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: Grids
MODULE u_to_p_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE u_to_p(array_on_u_points,                              &
                  ini_start,ini_end,                              &
                  inj_start,inj_end,                              &
                  outi_start,outi_end,                            &
                  outj_start,outj_end,                            &
                  levels,model_domain,                            &
                  at_extremity, array_on_p_points)

  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE domain_params
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN)  ::   ini_start,ini_end
INTEGER, INTENT(IN)  ::   inj_start,inj_end
INTEGER, INTENT(IN)  ::   outi_start,outi_end
INTEGER, INTENT(IN)  ::   outj_start,outj_end
INTEGER, INTENT(IN)  ::   levels
INTEGER, INTENT(IN)  ::   model_domain
        

LOGICAL, INTENT(IN)  ::   at_extremity(4)  
                   ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid


REAL, INTENT(IN)  :: array_on_u_points(ini_start:ini_end,         &
                                       inj_start:inj_end,         &
                                       levels)

REAL, INTENT(OUT) :: array_on_p_points(outi_start:outi_end,       &
                                       outj_start:outj_end,       &
                                       levels) 

! local variables

INTEGER            :: i,j,k
INTEGER            :: j0,j1  ! local variables

! Parameters
INTEGER, PARAMETER ::  pnorth =1  
                          ! North processor address in the neighbor array
INTEGER, PARAMETER ::  peast  =2  
                          ! East processor address in the neighbor array
INTEGER, PARAMETER ::  psouth =3
                          ! South processor address in the neighbor array
INTEGER, PARAMETER ::  pwest  =4  
                          ! West processor address in the neighbor array

INTEGER, PARAMETER ::  nodomain =-1     
                           ! Value in neighbor array if the domain has
                           ! no neighbor in this direction. Otherwise
                           ! the value will be the tid of the neighbor

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('U_TO_P',zhook_in,zhook_handle)

j0 = outj_start
j1 = outj_end

IF (.NOT. l_vatpoles) THEN
  IF (model_domain  ==  mt_global) THEN
  ! Do not do poles as this will be done by polar vector wind.
    IF (at_extremity(psouth) ) THEN
      j0 = outj_start+1
    END IF
    IF (at_extremity(pnorth) ) THEN
      j1 = outj_end-1
    END IF
  END IF
END IF ! vatpoles

DO k=1,levels
  DO j= j0, j1

    DO i= outi_start,outi_end

      array_on_p_points(i,j,k)= 0.5 *                             &
    ( array_on_u_points(i,j,k) + array_on_u_points(i-1,j,k) )




    END DO

  END DO
END DO

IF (lhook) CALL dr_hook('U_TO_P',zhook_out,zhook_handle)
RETURN
END SUBROUTINE u_to_p
END MODULE u_to_p_mod
