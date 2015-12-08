! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Convert snow amount from m to mm

MODULE rcf_grib_spcl_snow_mod

! SUBROUTINE Rcf_Grib_Spcl_Snow

! Description:
!   This routine multiplies the values by 1000 to convert
!   snow amount from m to mm.
!   At present there is a double check to ensure the data
!   came from ECMWF as not all input data will be volumetric.
!
! Method:
!   Multiply by 1000.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.4 programming standards.

CONTAINS
SUBROUTINE rcf_grib_spcl_snow(current,fielddata)

USE rcf_grib_block_params_mod, ONLY :                                          &
  grib_record,                                                                 &
  lenarraymax,                                                                 &
  p_orig_cntr

USE rcf_grib_lookups_mod, ONLY :                                               &
  grborigecmwf

USE printstatus_mod, ONLY :                                                    &
  printstatus,                                                                 &
  prstatus_normal,                                                             &
  prstatus_diag

USE um_parvars, ONLY :                                                         &
  mype

USE ereport_mod, ONLY :                                                        &
  ereport

IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(InOut):>
TYPE (grib_record),POINTER          :: current
REAL, INTENT(inout)                 :: fielddata(lenarraymax)

! Local constants
CHARACTER (len=*), PARAMETER     :: routinename='Rcf_GRIB_Spcl_Snow'

! Local variables

CHARACTER (len=80)               :: cmessage(2)   ! used for EReport
INTEGER                          :: errorstatus   ! used for EReport

!=======================================================================
!  Routine Code Start :
!=======================================================================

! double check it's an ECMWF field
IF (current % block_1(p_orig_cntr) == grborigecmwf) THEN

  IF ( printstatus >= prstatus_diag  ) THEN
    IF ( mype == 0 ) THEN
      WRITE (6,'(A)') "Converting snow amount from m of water to kg/m2"
    END IF
  END IF

  ! Multiply values by 1000
  fielddata(:) = 1000 * fielddata(:)

ELSE
  WRITE (cmessage(1),*) 'Data did not come from ECMWF'
  errorstatus = 10
  CALL ereport( routinename, errorstatus, cmessage(1) )
END IF

RETURN

END SUBROUTINE rcf_grib_spcl_snow
END MODULE rcf_grib_spcl_snow_mod
