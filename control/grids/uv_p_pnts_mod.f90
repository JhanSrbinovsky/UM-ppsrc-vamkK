! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids

MODULE uv_p_pnts_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE uv_p_pnts(u_single,v_single,cos_theta_longitude,sin_theta_longitude, &
  model_domain,global_row_length,proc_row_group,u_p,v_p)

! Purpose:
! Take a single level field of u and v and get u and v on the p grid.
! For use with diagnostic 3.230, set_seasalt and dms_flux
! and 3.436 Wind gust

! Code Description:
!   Language:Fortran 90.
!   This code is written to UM programming standards version 8.3.

USE dynamics_grid_mod, ONLY: l_vatpoles

USE atm_fields_bounds_mod, ONLY:   udims, vdims, udims_s, vdims_s,     &
                                   pdims, pdims_s, tdims
USE swapable_field_mod, ONLY :     swapable_field_pointer_type
USE proc_info_mod, ONLY:           at_extremity
USE um_parparams, ONLY:            pnorth, peast, psouth, pwest
USE Field_Types, ONLY:             fld_type_u, fld_type_v

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE domain_params
USE u_to_p_mod, ONLY: u_to_p
USE v_to_p_mod, ONLY: v_to_p

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

REAL,INTENT(IN)    :: u_single                                         &
  (udims%i_start:udims%i_end,udims%j_start:udims%j_end)
REAL,INTENT(IN)    :: v_single                                         &
  (vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
REAL,INTENT(IN)    :: cos_theta_longitude                              &
  (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL,INTENT(IN)    :: sin_theta_longitude                              &
  (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

INTEGER,INTENT(IN) :: model_domain
INTEGER,INTENT(IN) :: global_row_length ! number of points on a global row
INTEGER,INTENT(IN) :: proc_row_group
REAL,INTENT(OUT)   :: u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL,INTENT(OUT)   :: v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Local alloctable arrays
REAL, ALLOCATABLE, TARGET :: u_halo(:,:)
REAL, ALLOCATABLE, TARGET :: v_halo(:,:)

! The following are declared as arrays rather than scalars to maintain
! consistency of argument type, since they are passed to a general 
! purpose subroutine in which they are array arguments
REAL :: mag_vector_np(1)
REAL :: mag_vector_sp(1)
REAL :: dir_vector_np(1)
REAL :: dir_vector_sp(1)

INTEGER :: i,j,i_field,y,z     ! counters

! Temporarily create these as other subroutines use them
! can be removed when all vatpoles changes are made
INTEGER :: offx,offy,row_length,n_rows

TYPE(swapable_field_pointer_type) :: fields_to_swap(2)   ! for multivar
                                                         ! swap_bounds
INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UV_P_PNTS',zhook_in,zhook_handle)

! Temporarily create these as other subroutines use them
! can be removed when all vatpoles changes are made
offx=pdims_s%i_end-pdims%i_end
offy=pdims_s%j_end-pdims%j_end
row_length=pdims%i_end
n_rows=vdims%j_end-vdims%j_start+1

! Set allocated arrays to zero to make sure any data left in the memory
! allocated to them is purged
u_p(:,:)   =0.0
v_p(:,:)   =0.0

! Allocate arrays to hold u and v on u/v grid with halos
ALLOCATE(u_halo(udims_s%i_start:udims_s%i_end,                                 &
  udims_s%j_start:udims_s%j_end))
ALLOCATE(v_halo(vdims_s%i_start:vdims_s%i_end,                                 &
  vdims_s%j_start:vdims_s%j_end))


! Set allocated arrays to zero to make sure any data left in the memory
! allocated to them is purged
u_halo(:,:)=0.0
v_halo(:,:)=0.0


! Set non halo points to input values of u and v
DO y=udims%j_start,udims%j_end
  DO z=udims%i_start,udims%i_end
    u_halo(z,y)=u_single(z,y)
  END DO
END DO
DO y=vdims%j_start,vdims%j_end
  DO z=vdims%i_start,vdims%i_end
    v_halo(z,y)=v_single(z,y)
  END DO
END DO

! Update halos for u_halo and v_halo
i_field = 0
i_field = i_field + 1
fields_to_swap(i_field) % field_2d    => u_halo(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_u
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  udims%j_end
fields_to_swap(i_field) % vector      =  .TRUE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d    => v_halo(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_v
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  n_rows
fields_to_swap(i_field) % vector      =  .TRUE.

! DEPENDS ON: swap_bounds_2d_mv
CALL swap_bounds_2d_mv( fields_to_swap, i_field,                               &
  pdims%i_end, offx, offy)

! interpolate u and v to p grid.
CALL u_to_p(u_halo,                                                            &
  udims_s%i_start,udims_s%i_end,                                               &
  udims_s%j_start,udims_s%j_end,                                               &
  pdims%i_start,pdims%i_end,                                                   &
  pdims%j_start,pdims%j_end,                                                   &
  1,                                                                           &
  model_domain,at_extremity,u_p)


CALL v_to_p(v_halo,                                                            &
  vdims_s%i_start,vdims_s%i_end,                                               &
  vdims_s%j_start,vdims_s%j_end,                                               &
  pdims%i_start,pdims%i_end,                                                   &
  pdims%j_start,pdims%j_end,                                                   &
  1,                                                                           &
  model_domain,at_extremity,v_p)

! Global models on the New dynamics grid need to have special code
! to handle u and v at the poles
IF (.NOT. l_vatpoles) THEN
IF (model_domain  ==  mt_global) THEN
  ! Overwrite values of U_P, V_P at the poles with the magnitude of
  ! the vector wind.
  ! DEPENDS ON: polar_vector_wind_n
  CALL polar_vector_wind_n(                                                    &
    v_halo,                                                                    &
    sin_theta_longitude,                                                       &
    cos_theta_longitude,vdims%i_end,                                           &
    n_rows, 1 , mag_vector_np,                                                 &
    dir_vector_np, mag_vector_sp,                                              &
    dir_vector_sp,                                                             &
    offx, offy, global_row_length,                                             &
    proc_row_group, at_extremity)

  IF (at_extremity(psouth) ) THEN
    DO i=pdims%i_start,pdims%i_end
      v_p(i,1) = mag_vector_sp(1)
      u_p(i,1) = 0.0
    END DO
  END IF

  IF (at_extremity(pnorth) ) THEN
    DO i=pdims%i_start,pdims%i_end
      v_p(i,pdims%j_end) = mag_vector_np(1)
      u_p(i,pdims%j_end) = 0.0
    END DO
  END IF

END IF
END IF ! vatpoles

DEALLOCATE(u_halo)
DEALLOCATE(v_halo)

IF (lhook) CALL dr_hook('UV_P_PNTS',zhook_out,zhook_handle)

END SUBROUTINE uv_p_pnts

END MODULE uv_p_pnts_mod
