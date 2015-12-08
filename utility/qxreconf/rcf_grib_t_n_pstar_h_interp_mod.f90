! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Module used to hold data to allow recall later in execution tree.

Module Rcf_GRIB_T_n_Pstar_H_Interp_Mod

! Description:
!   Stores the horizontally interpolated versions of T and Pstar used
!   to generate heights when reconfiguring from ECMWF pressure levels.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5   30/08/02   Original code.  R.Sharp
!   6.2   20/06/05   Added ECMWF model level definition
!                    variables. Paul Earnshaw (frpe)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Use Rcf_Field_Type_Mod, Only :  &
  field_type

Implicit None

Type (field_type), Target, Save :: GRIB_T
Type (field_type), Target, Save :: GRIB_Pstar

Real, Target, Allocatable, Save :: GRIB_Levels(:)

! Store ecmwf model level definitions
Real, Target, Allocatable, Save :: ak(:)
Real, Target, Allocatable, Save :: bk(:)

End Module Rcf_GRIB_T_n_Pstar_H_Interp_Mod
