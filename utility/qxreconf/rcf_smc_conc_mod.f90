! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_smc_conc_mod

!  Subroutine Rcf_smc_conc - calculate the soil moisture conc

! Description:
! This module calculate soil moisture conc fields,
! and then interpolate smc_conc into the new domain.
! This smc_conc will be converted into the soil moisture
! in the post process routine.
!

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CONTAINS

SUBROUTINE rcf_smc_conc(smc_field)

USE printstatus_mod, ONLY : &
    printstatus,            &
    prstatus_normal

USE um_parvars, ONLY : &
    mype

USE rcf_grid_type_mod, ONLY : &
    input_grid

USE rcf_field_type_mod, ONLY : &
    field_type

USE missing_data_mod, ONLY : &
    rmdi

IMPLICIT NONE

! Arguments
TYPE( field_type), INTENT( INOUT )   :: smc_field

! Local variables
INTEGER                              :: j
INTEGER                              :: k

IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE (6,'(a)') ' Converting SMCL to soil-moisture concentration'
END IF

! Calculate smc_conc and use the new data for smc:
DO j = 1, smc_field % levels
    DO k=1, smc_field % level_size

      IF (ABS(smc_field % DATA(k,j) - rmdi) < EPSILON(rmdi)) THEN
        smc_field % DATA(k,j) = rmdi
      ELSE
        smc_field % DATA(k,j) = smc_field % DATA(k,j)/ &
        ( input_grid % soil_depths(j) )
      END IF

    END DO
END DO

RETURN
END SUBROUTINE rcf_smc_conc
END MODULE rcf_smc_conc_mod
