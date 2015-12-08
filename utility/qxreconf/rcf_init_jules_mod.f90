! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ initialises JULES variables from namelists

MODULE rcf_init_jules_mod
IMPLICIT NONE

! Description:
! Initialises JULES variables using the namelist information
! already read in.

! Method:
!  Allocates and maps variables as required.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE rcf_init_snow_param ( )

USE switches, ONLY:                                                    &
  can_model,                                                           &
  l_aggregate

USE nstypes, ONLY : &
    npft

USE nlsizes_namelist_mod, ONLY : &
    ntiles

USE ancil_info, ONLY:                                                  &
  nsmax

USE snow_param, ONLY:                                                  &
! namelist variables:
  dzsnow_io,       cansnowpft,                                         &
! non-namelist variables:
  dzsnow,          cansnowtile

IMPLICIT NONE

! jules_snow_param
! dzsnow is only allocated if nsmax > 0
IF ( nsmax > 0 ) THEN
  ALLOCATE(dzsnow( nsmax ))
  dzsnow(:) = dzsnow_io(1:nsmax)
END IF

! Set up cansnowtile
ALLOCATE(cansnowtile( ntiles ))
cansnowtile(:) = .FALSE.
IF ( .NOT. l_aggregate .AND. can_model == 4 ) THEN
  cansnowtile(1:npft) = cansnowpft(1:npft)
END IF

RETURN
END SUBROUTINE rcf_init_snow_param

END MODULE rcf_init_jules_mod
