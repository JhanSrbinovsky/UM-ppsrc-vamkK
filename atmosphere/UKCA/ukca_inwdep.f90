! *****************************COPYRIGHT******************************* 
! 
! (c) [University of Cambridge] [2008]. All rights reserved. 
! This routine has been licensed to the Met Office for use and 
! distribution under the UKCA collaboration agreement, subject  
! to the terms and conditions set out therein. 
! [Met Office Ref SC138]  
! 
! *****************************COPYRIGHT******************************* 
!                                                                       
! Purpose: Subroutine to read in coefficients for calculating           
!          effective Henry's Law coefficients. Original version         
!          taken from the Cambridge TOMCAT model.                       
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
!          Called from ASAD routine ASAD_cinit.                         
!                                                                       
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!                                                                       
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
        SUBROUTINE UKCA_INWDEP                                          

        USE ASAD_MOD,             ONLY: kd298, k298, ddhr, dhr
        USE UKCA_CHEM_DEFS_MOD,   ONLY: henry_defs
        USE ukca_option_mod,      ONLY: jpdw
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE ereport_mod, ONLY : ereport
        USE Control_Max_Sizes
        IMPLICIT NONE                                                   
                                                                        
        INTEGER :: errcode                ! Variable passed to ereport
        INTEGER            :: ns       ! Loop variable                      

        CHARACTER (LEN=72) :: cmessage ! String for error message

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

        IF (lhook) CALL dr_hook('UKCA_INWDEP',zhook_in,zhook_handle)

!       Use module to define Henry constants

        IF (SIZE(henry_defs) /= jpdw*6) THEN
          cmessage='jpdw and henry_defs are inconsistent'
          errcode=1
          CALL EREPORT('UKCA_INWDEP',errcode,cmessage)
        ENDIF

        DO ns=1,jpdw
          k298(ns)    = henry_defs(1,ns)
          dhr(ns)     = henry_defs(2,ns)
          kd298(ns,1) = henry_defs(3,ns)
          ddhr(ns,1)  = henry_defs(4,ns)
          kd298(ns,2) = henry_defs(5,ns)
          ddhr(ns,2)  = henry_defs(6,ns)
        ENDDO

        IF (lhook) CALL dr_hook('UKCA_INWDEP',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE UKCA_INWDEP  
