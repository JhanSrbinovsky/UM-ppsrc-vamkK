! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Stores magic numbers used in energy correction
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Energy Correction
!

MODULE eng_mass_param_mod

      INTEGER, PARAMETER :: ip_dry_mass = 1
      INTEGER, PARAMETER :: ip_wet_mass = 2
      INTEGER, PARAMETER :: ip_cvT      = 3
      INTEGER, PARAMETER :: ip_gr       = 4
      INTEGER, PARAMETER :: ip_keu      = 5
      INTEGER, PARAMETER :: ip_kev      = 6
      INTEGER, PARAMETER :: ip_kew      = 7
      INTEGER, PARAMETER :: ip_q        = 8
      INTEGER, PARAMETER :: ip_qcl      = 9
      INTEGER, PARAMETER :: ip_qcf      = 10
      INTEGER, PARAMETER :: ip_qu       = 1
      INTEGER, PARAMETER :: ip_qv       = 2
      INTEGER, PARAMETER :: ip_qw       = 3

      INTEGER, PARAMETER :: n_sums = 10
      INTEGER, PARAMETER :: n_flux = 3

END MODULE
