! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE Ancil_Mod

!
! Description:
! Define structure to hold Ancillary FIELD information read in
! from ANCILmaster records

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
      IMPLICIT NONE
!
TYPE ANC_record_type
  SEQUENCE
  INTEGER                  :: ancil_ref_number
  INTEGER                  :: model_number
  INTEGER                  :: section_number
  INTEGER                  :: item_number
  INTEGER                  :: anc_file_number
  INTEGER                  :: anc_field_read

  CHARACTER (LEN=36)       :: anc_name
  CHARACTER (LEN=1)        :: anc_flag
END TYPE ANC_record_type

INTEGER, PARAMETER         :: AncRec_IntDataLen  = 5
INTEGER, PARAMETER         :: AncRec_CharDataLen = 36

TYPE (ANC_record_type), ALLOCATABLE, SAVE  :: anc_record(:)

INTEGER, SAVE              :: ancRecs     = 0     ! No. of ancil records
INTEGER, SAVE              :: max_ancRecs = 0 ! Max no of ancil records


! Define structure to hold Ancillary FILE information read in
! from ANCILmaster records

TYPE ANC_file_type
  SEQUENCE
  INTEGER                  :: anc_file_number
  INTEGER                  :: model_number
  INTEGER                  :: anc_file_open
  CHARACTER (LEN=8)        :: anc_env_var
  CHARACTER (LEN=41)       :: anc_file_title
  CHARACTER (LEN=1)        :: anc_flag
END TYPE ANC_file_type

INTEGER, PARAMETER         :: AncFile_IntDataLen  = 2
INTEGER, PARAMETER         :: AncFile_CharDataLen = 49

TYPE (ANC_file_type), ALLOCATABLE, SAVE  :: anc_file(:)

INTEGER , SAVE :: ancFiles     = 0     ! No of anc files
INTEGER , SAVE :: max_ancFiles = 0     ! Max no of anc files
INTEGER , SAVE :: AncF_UnitNo          ! Unit No for ancillary files


! Namelist for user ANCILmaster file

INTEGER, SAVE              :: n_uancil
CHARACTER (LEN=80), SAVE   :: uancfils(20)

! For ROSE - disabling this code as the name list isn't suitable for
! ROSE in its present form. Have commented out the code rather than deleting it,
! in case there is a need to reactivate it in the future - although the relation
! between name list and user interface would have to be re-designed.

! Namelist /uancnum/ &
! n_uancil, uancfils


! Allocatable work arrays for ancillary processing

INTEGER, ALLOCATABLE :: nlookup(:)
INTEGER, ALLOCATABLE :: lookup_step(:)
INTEGER, ALLOCATABLE :: stashancil(:)
! INTEGER, ALLOCATABLE :: fieldcode(:)
INTEGER, ALLOCATABLE :: levels(:)
INTEGER, ALLOCATABLE :: ancil_add(:)

END MODULE Ancil_Mod

