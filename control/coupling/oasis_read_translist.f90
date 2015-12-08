! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE OASIS_READ_TRANSLIST()

!
! Description:
! Read translist to help us set up a full list of
! fields involved in coupling.
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Coupling
! 
!==================================================================

USE oasis_atm_data_mod, ONLY: TRANSCOUNT,    &
                              transient_in,  &  
                              transient_out, &
                              c_in,          & 
                              n_in,          &
                              c_out,         & 
                              n_out,         &
                              i_number,      & 
                              f_id,          &
                              my_component,  &
                              transfld,      &
                              lb,            & 
                              max_transients 
                              
                                
USE Field_Types, ONLY : fld_type_u, fld_type_v

  IMPLICIT NONE

  INTEGER :: iostatus

  OPEN(LB,file="translist",status="OLD",Iostat=iostatus)

  ! Obtain the maximum index used by our coupling transients
  READ( Unit=LB, Nml=TRANSCOUNT, Iostat=iostatus )

  ALLOCATE (transient_in(max_transients))
  ALLOCATE (transient_out(max_transients))

  ! Initialise all our transient details with 
  ! default conditions. 
  transient_out(:)%indx = -999
  transient_in(:)%indx = -999

  transient_out(:)%grid = ""
  transient_in(:)%grid = ""

  transient_out(:)%name = ""
  transient_in(:)%name = ""

  transient_out(:)%field_id = -999
  transient_in(:)%field_id = -999
       
  ! Set default processing conditions:
  ! The ENDGAME case does not need to do polar meaning and nor do LAMs
  transient_in(:)%polar_mean = .FALSE.

 
  ! Ice processing applies regardless of grid type
  transient_in(:)%ice_min = .FALSE.

  ! Read all the transient fields from our namelist and save the 
  ! types and names in the appropriate in or out arrays. 
  DO WHILE ( iostatus == 0 )

    READ( Unit=LB, Nml=TRANSFLD, Iostat=iostatus )

    IF (c_out(1:3) == my_component) THEN
       transient_out(i_number)%grid=c_out(4:4)
       transient_out(i_number)%name=n_out
       transient_out(i_number)%indx=i_number
       transient_out(i_number)%field_id=f_id
    END IF

    IF (c_in(1:3) == my_component) THEN
       transient_in(i_number)%grid=c_in(4:4)
       transient_in(i_number)%name=n_in
       transient_in(i_number)%indx=i_number
       transient_in(i_number)%field_id=f_id
    END IF

  END DO

  CLOSE(LB)

  END SUBROUTINE OASIS_READ_TRANSLIST


