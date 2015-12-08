! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_lsh_land_ice_chk_mod

!  Subroutine rcf_lsh_land_ice_chk

! Description:
!   Ensures that land-ice are consistent with LSH fields.

! Method:
! Some LSH fields are sensitive to land-ice values so lets make sure we are
! consistent.

! Current Code Owner: T. Green

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE rcf_lsh_land_ice_chk( fields,      &
                                 field_count, &
                                 decomp, hdr )

USE rcf_field_type_mod, ONLY : &
    field_type

USE rcf_umhead_mod, ONLY : &
    um_header_type

USE rcf_locate_mod, ONLY : &
    rcf_locate

USE rcf_read_field_mod, ONLY : &
    rcf_read_field

USE rcf_write_field_mod, ONLY : &
    rcf_write_field

USE rcf_alloc_field_mod, ONLY : &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE rcf_stashcodes_mod, ONLY : &
    stashcode_frac_surf_type,  &
    stashcode_zw,              &
    stashcode_fsat,            &
    stashcode_fwetl,           &
    stashcode_prog_sec

USE um_parvars, ONLY :         &
    nproc,                     &
    mype

USE printstatus_mod, ONLY :    &
    printstatus,               &
    prstatus_normal

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields(:)
TYPE( um_header_type ), INTENT(IN) :: hdr
INTEGER, INTENT(IN)                :: field_count
INTEGER, INTENT(IN)                :: decomp

! Local variables
INTEGER, PARAMETER                 :: land_ice_lvl = 9
INTEGER, PARAMETER                 :: stashlist_size = 3
INTEGER, PARAMETER                 :: stashlist(stashlist_size) = &
                                      (/ stashcode_zw,            &
                                         stashcode_fsat,          &
                                         stashcode_fwetl /)

INTEGER                            :: i, j
INTEGER                            :: count
INTEGER                            :: pos
INTEGER                            :: istat
REAL                               :: value_on_land_ice

TYPE( field_type ), POINTER        :: frac_surf
TYPE( field_type ), POINTER        :: field

! Include files
! C_TOPOG start
!
! Topographic index increment:
      REAL,PARAMETER :: DTI = 0.2
! Maximum topographic index considered:
      REAL,PARAMETER :: TI_MAX = 10.0
! Maximum allowed water table depth (m):
      REAL,PARAMETER :: ZW_MAX = 6.0
! Standard deviation of LOG(Ksat(0)):
      REAL,PARAMETER :: SIGMA_LOGK = 0.0
! Parameter to remove very high water tables
! from the calculated wetland fraction:
      REAL,PARAMETER :: TI_WETL = 1.5
!
! C_TOPOG end

!-------------------------------------------------------------------
! Locate required parameters in output dump
!-------------------------------------------------------------------

CALL rcf_locate( stashcode_prog_sec, stashcode_frac_surf_type,         &
                 fields, field_count, pos, .TRUE. )

! If we have surface fractions available then lets check land ice level.
IF (pos /= 0 ) THEN
  frac_surf => fields( pos )
  CALL rcf_alloc_field( frac_surf )
  CALL rcf_read_field( frac_surf, hdr, decomp )

  ! Loop over all stashcodes which we want to check.
  DO j = 1, stashlist_size
    count = 0
    CALL rcf_locate( stashcode_prog_sec, stashlist(j),                 &
                     fields, field_count, pos, .TRUE.)
    IF (pos /= 0) THEN
      field => fields(pos)
      CALL rcf_alloc_field(field)
      CALL rcf_read_field(field, hdr, decomp)

      IF (stashlist(j) == stashcode_zw) THEN
! Use value from c_topog header.
        value_on_land_ice = zw_max
      ELSE
! Currently other stashcodes checks set to zero.
        value_on_land_ice = 0.0
      END IF

      DO i = 1, field % level_size
        ! Land-ice fraction should be 0 or 1 so lets just check for positive.
        IF ( frac_surf % data ( i, land_ice_lvl ) > 0.0 .AND. &
             field % data ( i, 1 ) /= value_on_land_ice ) THEN
          ! Set new value in field over land ice.
          field % data (i, 1) = value_on_land_ice
          count = count + 1
        END IF
      END DO

      CALL gc_isum (1, nproc, istat, count)
      IF (printstatus >= prstatus_normal .AND. mype == 0) THEN
        WRITE (6,*) 'LSH land-ice check: No of values reset ', &
                     count, ' for stashcode ', stashlist(j)
      END IF
!---------------------------------------------------------------------
! Write out changed field
!---------------------------------------------------------------------
      IF (count > 0) THEN
        CALL rcf_write_field( field, hdr, decomp )
      END IF
      CALL rcf_dealloc_field( field )
    END IF
  END DO
END IF

RETURN
END SUBROUTINE rcf_lsh_land_ice_chk
END MODULE rcf_lsh_land_ice_chk_mod
