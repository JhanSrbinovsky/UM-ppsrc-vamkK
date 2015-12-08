! *****************************COPYRIGHT******************************* 
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************* 
!            
! Purpose: To hold magic numbers associated with UKCA chemistry 
! schemes
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                          
!          Called from SET_ATM_POINTERS
!                             
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!    
!        
! Code description:            
!   Language: FORTRAN 95       
!   This code is written to UMDP3 programming standards.        
!  
! --------------------------------------------------------------------- 
!
MODULE ukca_chem_schemes_mod

IMPLICIT NONE

! Solver types
INTEGER, PARAMETER :: int_method_impact = 1
INTEGER, PARAMETER :: int_method_nr = 3
INTEGER, PARAMETER :: int_method_be = 5
INTEGER, PARAMETER :: int_method_be_explicit = 10

! Magic numbers for each chemistry scheme 
INTEGER, PARAMETER :: i_ukca_chem_off        = 0
INTEGER, PARAMETER :: i_ukca_chem_ageofair   = 1
!
! Backward Euler schemes
INTEGER, PARAMETER :: i_ukca_chem_trop       = 11
INTEGER, PARAMETER :: i_ukca_chem_raq        = 13
!
! Newton Raphson schemes
INTEGER, PARAMETER :: i_ukca_chem_tropisop   = 50
INTEGER, PARAMETER :: i_ukca_chem_strattrop  = 51
INTEGER, PARAMETER :: i_ukca_chem_strat      = 52
INTEGER, PARAMETER :: i_ukca_chem_std_trop   = 53

END MODULE ukca_chem_schemes_mod


