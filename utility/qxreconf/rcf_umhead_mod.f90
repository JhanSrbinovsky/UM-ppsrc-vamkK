! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Data module defining um_header_type

MODULE Rcf_UMhead_Mod

! Description:
!  Contains parameters used in UM dumps,
!  a structure to hold a UM header
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

!   Global Constants
INTEGER, PARAMETER :: LenFixHd = 256   ! Length of Fixed Length Header.

!   Global Type Definitions:
TYPE UM_header_type
  SEQUENCE

  INTEGER ::      LenIntC         ! Length of Integer Constants array.
  INTEGER ::      LenRealC        ! Length of Real Constants array.
  INTEGER ::      Len1LevDepC     ! 1st dimension \  of Level Dependent
  INTEGER ::      Len2LevDepC     ! 2nd dimension /  Constants array.
  INTEGER ::      Len1RowDepC     ! 1st dimension \  of Row Dependent
  INTEGER ::      Len2RowDepC     ! 2nd dimension /  Constants array.
  INTEGER ::      Len1ColDepC     ! 1st dimension \  of Column Dependent
  INTEGER ::      Len2ColDepC     ! 2nd dimension /  Constants array.
  INTEGER ::      Len1FldsOfC     ! 1st dimension \  of Fields of
  INTEGER ::      Len2FldsOfC     ! 2nd dimension /  Constants array.
  INTEGER ::      LenExtraC       ! Length of Extra Constants array.
  INTEGER ::      LenHistFile     ! Length of History File.
  INTEGER ::      LenCompFldI1    ! Length of Compressed Field Index 1.
  INTEGER ::      LenCompFldI2    ! Length of Compressed Field Index 2.
  INTEGER ::      LenCompFldI3    ! Length of Compressed Field Index 3.
  INTEGER ::      Len1Lookup      ! 1st dimension of Lookup table.
  INTEGER ::      Len2Lookup      ! 2nd dimension of Lookup table.
  INTEGER ::      LenData         ! Length of Data array.
  INTEGER ::      StartData       ! Position of start of Data array.
  INTEGER ::      NumFlds         ! Number of data fields.
  INTEGER ::      MaxFldSize      ! Maximum size of field.
  INTEGER ::      UnitNum         ! Unit number associated with UM dump.

  INTEGER, POINTER :: FixHd    (:) => NULL() ! Fixed length header.
  INTEGER, POINTER :: IntC     (:) => NULL() ! Integer Constants array.
  INTEGER, POINTER :: CompFldI1(:) => NULL() ! Compressed Field Index array 1.
  INTEGER, POINTER :: CompFldI2(:) => NULL() ! Compressed Field Index array 2.
  INTEGER, POINTER :: CompFldI3(:) => NULL() ! Compressed Field Index array 3.
  INTEGER, POINTER :: Lookup(:,:)  => NULL() ! Lookup table.

  REAL, POINTER    :: RealC  (:)   => NULL() ! Real Constants array.
  REAL, POINTER    :: LevDepC(:,:) => NULL() ! Level Dependent Constants array.

!Use 2-dimentional Row and Column dependent consts.
  REAL, POINTER    :: RowDepC(:,:) => NULL() ! Row Dependent Const. array.
  REAL, POINTER    :: ColDepC(:,:) => NULL() ! Column Dependent Const. array.

  REAL, POINTER    :: FldsOfC (:)  => NULL() ! Field Dependent Constants array.
  REAL, POINTER    :: ExtraC  (:)  => NULL() ! Extra Constants array.
  REAL, POINTER    :: HistFile(:)  => NULL() ! History File.

END TYPE UM_header_type

END MODULE Rcf_UMhead_Mod
