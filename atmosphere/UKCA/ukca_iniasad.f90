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
! Purpose: Subroutine to initialise ASAD and fill ldepd and ldepw arrays
!          Adapted from original version written by Olaf Morgenstern.   
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_MAIN1.                                      
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
        SUBROUTINE UKCA_INIASAD(npoints)                                

        USE UKCA_CHEM_DEFS_MOD,   ONLY: chch_t, chch_defs
        USE ukca_option_mod,      ONLY: jpctr, jpspec, jpnr, jpbk, jptk,  &
                                        jppj, jphk, jpdd, jpdw
        USE ASAD_MOD,             ONLY: ldepd, ldepw, ih_o3, ih_h2o2,     &
                                        ih_hno3, ih_so2, asad_mod_init
        USE parkind1,             ONLY: jprb, jpim
        USE yomhook,              ONLY: lhook, dr_hook
        USE ereport_mod,          ONLY: ereport
        USE PrintStatus_mod
        USE Control_Max_Sizes
        IMPLICIT NONE                                                   
                                                                        
        INTEGER, INTENT(IN) :: npoints   ! no of spatial points         

!       Local variables
                                                                        
        INTEGER :: errcode              ! Variable passed to ereport
        INTEGER :: k                    ! Loop variable                 
        INTEGER :: iw                   ! Loop variable                 
        INTEGER :: jerr(jpspec)         ! For error analysis
                                                                        
        CHARACTER (LEN=72) :: cmessage  ! Error message                 

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

        IF (lhook) CALL dr_hook('UKCA_INIASAD',zhook_in,zhook_handle)

        IF (printstatus >= prstatus_oper) THEN 
          write(6,*) 'ASAD initialised from namelist using:'
          write(6,*) 'jpctr: ',jpctr, ' jpspec: ',jpspec, ' jpnr: ',jpnr
          write(6,*) 'jpbk: ', jpbk,  ' jptk: ',  jptk,   ' jppj: ',jppj
          write(6,*) 'jphk: ', jphk,  ' jpdd: ',  jpdd,   ' jpdw: ',jpdw
        END IF

        CALL ASAD_MOD_INIT

! Set up dry and wet deposition logicals using module switches
        ldepd(:) = .false.
        ldepw(:) = .false.
        DO k=1,jpspec
          ldepd(k) = (chch_defs(k)%switch1 == 1)
          ldepw(k) = (chch_defs(k)%switch2 == 1)
        END DO

! Indentify index of henry_defs array for SO2 etc to use in asad_hetero
        iw = 0
        ih_o3 = 0
        ih_h2o2 = 0
        ih_hno3 = 0
        ih_so2 = 0
        DO k=1,jpspec
          IF (ldepw(k)) iw = iw + 1
          IF (chch_defs(k)%speci == 'O3        ') ih_o3 = iw
          IF (chch_defs(k)%speci == 'H2O2      ') ih_h2o2 = iw
          IF (chch_defs(k)%speci == 'HONO2     ' .OR.                   &
              chch_defs(k)%speci == 'HNO3      ') ih_hno3 = iw
          IF (chch_defs(k)%speci == 'SO2       ') ih_so2 = iw
        END DO

! Check if module sizes are compatible with UMUI values

        IF (size(chch_defs) /= jpspec) THEN
          cmessage='size of chch_defs inconsistent with jpspec'
          write(6,*) cmessage, size(chch_defs), jpspec
          errcode = 1
          CALL EREPORT('UKCA_INIASAD',errcode,cmessage)
        ENDIF

        jerr(:)=0
        WHERE (ldepd) jerr=1
        IF (SUM(jerr) /= jpdd) THEN
          cmessage='chch_defs%switch1 values inconsistent with jpdd'
          write(6,*) cmessage, SUM(jerr), jpdd
          errcode = 1
          CALL EREPORT('UKCA_INIASAD',errcode,cmessage)
        ENDIF

        jerr(:)=0
        WHERE (ldepw) jerr=1
        IF (SUM(jerr) /= jpdw) THEN
          cmessage='chch_defs%switch2 values inconsistent with jpdw'
          write(6,*) cmessage, SUM(jerr), jpdw
          errcode = 1
          CALL EREPORT('UKCA_INIASAD',errcode,cmessage)
        ENDIF

! DEPENDS ON: asad_cinit
        CALL ASAD_CINIT(npoints)
                                                                        
      IF (lhook) CALL dr_hook('UKCA_INIASAD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INIASAD                                    
