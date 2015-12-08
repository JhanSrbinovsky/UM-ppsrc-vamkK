! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module providing UKCA with type definitions describing chemistry
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE UKCA_CHEM_DEFS_MOD
IMPLICIT NONE

! Describes the chemical species used in the model
TYPE CHCH_T
  INTEGER           :: idum      ! dummy integer
  CHARACTER(LEN=10) :: speci     ! species name
  INTEGER           :: nodd      ! no of odd atoms
  CHARACTER(LEN=10) :: ctype     ! Species type
  CHARACTER(LEN=10) :: family    ! Family
  INTEGER           :: switch1   ! 1 if dry deposits
  INTEGER           :: switch2   ! 1 if wet deposits
  INTEGER           :: switch3   ! > 0 if an emitted species
ENDTYPE CHCH_T
PUBLIC CHCH_T

! Describes the bimolecular reaction rates
TYPE RATB_T
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL              :: K0        ! rate coeff param K0
  REAL              :: alpha     ! rate coeff param alpha
  REAL              :: beta      ! rate coeff param beta
  REAL              :: pyield1   ! product1 fractional yield
  REAL              :: pyield2   ! product2 fractional yield
  REAL              :: pyield3   ! product3 fractional yield
  REAL              :: pyield4   ! product4 fractional yield
ENDTYPE RATB_T
PUBLIC RATB_T

! Describes heterogenous reactions
TYPE RATH_T
  CHARACTER(LEN=10) :: react1    ! reactant1 name
  CHARACTER(LEN=10) :: react2    ! reactant2 name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL              :: pyield1   ! product yield
  REAL              :: pyield2   ! product yield
  REAL              :: pyield3   ! product yield
  REAL              :: pyield4   ! product yield
ENDTYPE RATH_T
PUBLIC RATH_T

! describes photolytic reactions
TYPE RATJ_T
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL              :: pyield1   ! product yield
  REAL              :: pyield2   ! product yield
  REAL              :: pyield3   ! product yield
  REAL              :: pyield4   ! product yield
  REAL              :: jfacta    ! quantum yield
  CHARACTER(LEN=10) :: fname     ! file name/label
ENDTYPE RATJ_T
PUBLIC RATJ_T

! Describes termolecular reactions
TYPE RATT_T
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  REAL              :: F         ! rate coeff param F
  REAL              :: K1        ! rate coeff param K1
  REAL              :: alpha1    ! rate coeff param alpha1
  REAL              :: beta1     ! rate coeff param beta1
  REAL              :: K2        ! rate coeff param K2
  REAL              :: alpha2    ! rate coeff param alpha2
  REAL              :: beta2     ! rate coeff param beta2
  REAL              :: pyield1   ! product1 fractional yield
  REAL              :: pyield2   ! product2 fractional yield
ENDTYPE RATT_T
PUBLIC RATT_T

! These are set, depending on chemistry selection, in routine
! chem1_init (contained in ukca_chem1_dat)
TYPE(CHCH_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: chch_defs
TYPE(RATB_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: ratb_defs
TYPE(RATH_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: rath_defs
TYPE(RATJ_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: ratj_defs
TYPE(RATT_T), DIMENSION(:),     ALLOCATABLE, PUBLIC, SAVE :: ratt_defs
REAL,         DIMENSION(:,:,:), ALLOCATABLE, PUBLIC, SAVE :: depvel_defs
REAL,         DIMENSION(:,:),   ALLOCATABLE, PUBLIC, SAVE :: henry_defs

END MODULE
