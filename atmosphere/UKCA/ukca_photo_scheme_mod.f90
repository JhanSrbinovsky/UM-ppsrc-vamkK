! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module to hold magic numbers relating to photolysis
!          scheme options
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8.1 programming standards.
!
! ---------------------------------------------------------------------
MODULE ukca_photo_scheme_mod
IMPLICIT NONE
INTEGER,PARAMETER :: i_ukca_nophot  = 0 ! Photolysis off. Default value
INTEGER,PARAMETER :: i_ukca_phot2d  = 1 ! offline 2D photolysis scheme
INTEGER,PARAMETER :: i_ukca_fastj   = 2 ! original version of Fast-J
INTEGER,PARAMETER :: i_ukca_fastjx  = 3 ! Fast-JX - updated and improved 
   
END MODULE ukca_photo_scheme_mod
