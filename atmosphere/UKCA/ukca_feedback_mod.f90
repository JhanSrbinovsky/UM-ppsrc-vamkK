! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!          To define positions of individual greenhouse 
!          gases in greenhouse gas array, grgas_field
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90 
!   This code is written to UMDP3 v8.5 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_feedback_mod
IMPLICIT NONE

INTEGER, PARAMETER :: p_o3   = 1
INTEGER, PARAMETER :: p_ch4  = 2
INTEGER, PARAMETER :: p_n2o  = 3
INTEGER, PARAMETER :: p_f11  = 4
INTEGER, PARAMETER :: p_f12  = 5
INTEGER, PARAMETER :: p_f113 = 6
INTEGER, PARAMETER :: p_f22  = 7
INTEGER, PARAMETER :: p_h2os = 8

END MODULE ukca_feedback_mod
