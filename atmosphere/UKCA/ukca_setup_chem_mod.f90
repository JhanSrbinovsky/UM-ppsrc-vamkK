! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  A routine to set up values such as the number of reactions used 
!  by each chemistry scheme
!
! Method: 
!  i_ukca_chem, l_ukca_chem_aero and l_ukca_trophet are read 
!  in from namelists
!  This routine uses these values to set the number or reactants, and the 
!  number of each type of reaction, plus the number of dry and wet 
!  deposited species
!  Also includes error checking for schemes which have not yet
!  been set up
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
!   Language:  FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_setup_chem_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE ukca_setup_chem

USE ukca_option_mod, ONLY:  i_ukca_chem, l_ukca_chem_aero,              &
   l_ukca_achem, l_ukca_aerchem, l_ukca_raq, l_ukca_trop,               &
   l_ukca_tropisop, l_ukca_std_trop, l_ukca_strattrop, l_ukca_stratcfc, &
   l_ukca_strat, l_ukca_ageair, l_ukca_chem, l_ukca_het_psc,            &
   l_ukca_mode, l_ukca_advh2o, l_ukca_trophet, ukca_int_method,         &
   jpctr, jpspec, jpbk, jptk, jphk, jppj, jpdd, jpdw, jpnr
USE ukca_chem_schemes_mod, ONLY: int_method_be_explicit, int_method_nr, &
   i_ukca_chem_off, i_ukca_chem_ageofair,                               & 
   i_ukca_chem_trop, i_ukca_chem_raq,                                   &
   i_ukca_chem_tropisop, i_ukca_chem_strattrop, i_ukca_chem_strat,      &
   i_ukca_chem_std_trop
USE ukca_d1_defs,    ONLY: ukca_item_sulpc
USE ereport_mod,     ONLY: ereport
USE PrintStatus_mod, ONLY: PrintStatus, PrStatus_Diag
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER (LEN=70) :: cmessage        ! Error message
INTEGER            :: errcode         ! Variable passed to ereport

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('ukca_setup_chem',zhook_in,zhook_handle)


SELECT CASE (i_ukca_chem)

  CASE (i_ukca_chem_off)
  ! Chemistry off completely
    l_ukca_chem     = .FALSE.
    ukca_int_method = 0
    jpctr           = 0
    jpspec          = 0
    jpbk            = 0
    jptk            = 0
    jppj            = 0
    jphk            = 0
    jpdd            = 0
    jpdw            = 0

  CASE (i_ukca_chem_ageofair)
  ! Age of air
    l_ukca_chem     = .FALSE.
    l_ukca_ageair   = .TRUE.
    ukca_int_method = 0
    jpctr           = 0
    jpspec          = 0
    jpbk            = 0
    jptk            = 0
    jppj            = 0
    jphk            = 0
    jpdd            = 0
    jpdw            = 0

  CASE (i_ukca_chem_trop)
    
    l_ukca_chem     = .TRUE.
    l_ukca_trop     = .TRUE.
    ukca_int_method = int_method_BE_explicit
    IF (l_ukca_chem_aero) THEN
      ! Tropospheric chemistry with Aerosols (BE)
      l_ukca_aerchem  = .TRUE.
      jpctr           = 33
      jpspec          = 53
      jpbk            = 88
      jptk            = 14
      jppj            = 20
      jphk            = 0
      jpdd            = 26
      jpdw            = 19
   
    ELSE
      ! Tropopspheric chemistry (BE) without aerosols
      jpctr            = 26
      jpspec           = 46
      jpbk             = 88
      jptk             = 14
      jppj             = 20
      jphk             = 0
      jpdd             = 22
      jpdw             = 15
    END IF

  CASE (i_ukca_chem_raq)
    ! Regional Air Quality (BE)
    l_ukca_chem     = .TRUE.
    l_ukca_raq      = .TRUE.
    ukca_int_method = int_method_BE_explicit
    jpctr           = 40
    jpspec          = 58
    jpbk            = 113
    jptk            = 12
    jppj            = 23
    jphk            = 0
    jpdd            = 16
    jpdw            = 19
    IF (l_ukca_chem_aero) THEN
      ! this is not available yet - abort
         WRITE(6,'(A)') 'Aerosol chemistry not available ' &
        // 'for RAQ yet'
         cmessage='Unsupported option combination'
         errcode=4
         CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF

  CASE (i_ukca_chem_tropisop)
    l_ukca_chem     = .TRUE.
    l_ukca_tropisop = .TRUE.
    ukca_int_method = int_method_nr
    jpctr           = 49
    jpspec          = 56
    jpbk            = 117
    jptk            = 15
    jppj            = 35
    jphk            = 0
    jpdd            = 33
    jpdw            = 24
    ! If aerosol chemistry on, add additional reactions
    IF (l_ukca_chem_aero) THEN
       l_ukca_achem     = .TRUE.
       jpctr            = jpctr + 11
       jpspec           = jpspec + 11
       jpbk             = jpbk + 10
       jptk             = jptk + 1
       jphk             = jphk + 5
       jpdd             = jpdd + 7
       jpdw             = jpdw + 4
       ! incremement further if Tropospheric het chemistry on
       IF (l_ukca_trophet) THEN
         jphk             = jphk + 2
       END IF
    ELSE 
       ! Can't have trophet reactions without aerosol chemistry
       IF (l_ukca_trophet) THEN
         WRITE(6,'(A)') 'trop het chem requires aerosol chemistry on'
         cmessage='Unsupported option choice'
         errcode=ABS(i_ukca_chem)
         CALL ereport('ukca_setup_chem',errcode,cmessage)
       END IF
    END IF

  CASE (i_ukca_chem_strattrop)
    l_ukca_chem      = .TRUE.
    l_ukca_strattrop = .TRUE.
    l_ukca_advh2o    = .TRUE.
    ukca_int_method = int_method_nr
    jpctr           = 71
    jpspec          = 75
    jpbk            = 198
    jptk            = 24
    jppj            = 56
    jphk            = 0
    jpdd            = 36
    jpdw            = 29
    ! If aerosol chemistry on, add additional reactions
    IF (l_ukca_chem_aero) THEN
       l_ukca_achem     =.TRUE.
       jpctr            = jpctr + 12
       jpspec           = jpspec + 12
       jpbk             = jpbk + 15
       jptk             = jptk + 1
       jppj             = jppj + 4
       jphk             = jphk + 3
       jpdd             = jpdd + 5
       jpdw             = jpdw + 5
       ! incremement further if Het PSC on 
       IF (l_ukca_het_psc) THEN
         jphk             = jphk + 5
       END IF
       ! incremement again if Tropospheric het chemistry on
       IF (l_ukca_trophet) THEN
         jphk             = jphk + 2
       END IF
    ELSE 
       ! Can't have trophet reactions without aerosol chemistry
       IF (l_ukca_trophet) THEN
         WRITE(6,'(A)') 'trop het chem requires aerosol chemistry on'
         cmessage='Unsupported option choice'
         errcode=ABS(i_ukca_chem)
         CALL ereport('ukca_setup_chem',errcode,cmessage)
       END IF
    END IF

  CASE (i_ukca_chem_strat)
    l_ukca_chem     = .TRUE.
    l_ukca_strat    = .TRUE.
    l_ukca_advh2o    = .TRUE.
    ukca_int_method = int_method_nr
    jpctr           = 37
    jpspec          = 41
    jpbk            = 113
    jptk            = 17
    jppj            = 34
    jphk            = 0
    jpdd            = 15
    jpdw            = 15
    ! If aerosol chemistry on, add additional reactions
    IF (l_ukca_chem_aero) THEN
       l_ukca_achem     =.TRUE.
       jpctr            = jpctr + 8
       jpspec           = jpspec + 8
       jpbk             = jpbk + 12
       jptk             = jptk + 1
       jppj             = jppj + 4
       jphk             = jphk + 3
       jpdd             = jpdd + 1
       jpdw             = jpdw + 1
       IF (l_ukca_het_psc) THEN
         jphk             = jphk + 5
       END IF
    END IF
    IF (l_ukca_trophet) THEN
     WRITE(6,'(A)') 'Strat chem does not support trop het chem'
     cmessage='Unsupported option choice'
     errcode=ABS(i_ukca_chem)
     CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF

  CASE DEFAULT
! If we get here, we don't know how to deal with this chemistry scheme,
! stop
     WRITE(6,'(A,1X,I6)') 'Unknown chemistry scheme: ',i_ukca_chem
     cmessage='Unknown chemistry scheme'
     errcode=ABS(i_ukca_chem)
     CALL ereport('ukca_setup_chem',errcode,cmessage)

END SELECT

! calculate jpnr from the sum of the reactions
jpnr = jpbk + jptk + jppj + jphk

! Set item numbers for UKCA oxidants used in the CLASSIC sulphur cycle:
! Settings for B-E solvers with OH and HO2 as non-transported species
!       001 : O3 MASS MIXING RATIO AFTER TSTEP
!       007 : HONO2 MASS MIXING RATIO AFTER TSTEP
!       008 : H2O2 MASS MIXING RATIO AFTER TSTEP
!       153 : OH MASS MIXING RATIO  AFTER TSTEP
!       154 : HO2 MASS MIXING RATIO AFTER TSTEP

! Settings for N-R solvers with OH and HO2 as transported species
!        81 : OH MASS MIXING RATIO  AFTER TSTEP
!        82 : HO2 MASS MIXING RATIO AFTER TSTEP

      IF (ukca_int_method == int_method_nr) THEN
        ukca_item_sulpc(:) = (/1,7,8,81,82/)                       ! N-R solver
      ELSE
        ukca_item_sulpc(:) = (/1,7,8,153,154/)
      END IF

! Error checking 
IF (ukca_int_method  == int_method_BE_explicit) THEN
    ! If aerosol chemistry or tropospheric het chem on abort - not supported
    IF (l_ukca_chem_aero .OR. l_ukca_trophet) THEN
     WRITE(6,'(A)') 'BE schemes do not support add on aerosol ' &
        // 'chemistry or tropospheric heterogeneous chemistry'
     cmessage='Unsupported option combination'
     errcode=ABS(i_ukca_chem)
     CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF
END IF

! Tropospheric heterogeneous chemistry requires GLOMAP-mode
IF (l_ukca_trophet .AND. .NOT. l_ukca_mode) THEN
     WRITE(6,'(A)') 'Need Glomap MODE on for ' &
        // 'tropospheric heterogeneous chemistry'
     cmessage='Unsupported option combination'
     errcode=1
     CALL ereport('ukca_setup_chem',errcode,cmessage)
END IF

! Check chemistry supports GLOMAP-mode
IF (l_ukca_mode) THEN

   ! Newton Raphson 
   ! Require that l_ukca_achem true
   IF (ukca_int_method  == int_method_nr) THEN
       IF (.NOT. l_ukca_achem) THEN
         WRITE(6,'(A)') 'Need aerosol chemistry on for ' &
        // 'Glomap MODE'
         cmessage='Unsupported option combination'
         errcode=2
         CALL ereport('ukca_setup_chem',errcode,cmessage)
       END IF

   ! Backward Euler - only Trop + Aerosols supports 
   ! GLOMAP-mode
   ELSE
     IF (.NOT. l_ukca_aerchem) THEN
         WRITE(6,'(A)') 'Need aerosol chemistry on for ' &
        // 'Glomap MODE'
         cmessage='Unsupported option combination'
         errcode=3
         CALL ereport('ukca_setup_chem',errcode,cmessage)
     END IF
   END IF  
 
END IF


! Print the values we have set
IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(6,'(A)')       'Logicals and parameters set by ukca_setup_chem:'
  WRITE(6,'(A,1X,L7)') 'l_ukca_chem      = ', l_ukca_chem
  WRITE(6,'(A,1X,L7)') 'l_ukca_ageair    = ', l_ukca_ageair
  WRITE(6,'(A,1X,L7)') 'l_ukca_trop      = ', l_ukca_trop
  WRITE(6,'(A,1X,L7)') 'l_ukca_achem   = ',   l_ukca_achem
  WRITE(6,'(A,1X,L7)') 'l_ukca_aerchem   = ', l_ukca_aerchem
  WRITE(6,'(A,1X,L7)') 'l_ukca_tropisop  = ', l_ukca_tropisop
  WRITE(6,'(A,1X,L7)') 'l_ukca_strattrop = ', l_ukca_strattrop 
  WRITE(6,'(A,1X,L7)') 'l_ukca_strat     = ', l_ukca_strat
  WRITE(6,'(A,1X,I5)') 'ukca_int_method  = ', ukca_int_method
  WRITE(6,'(A,1X,I5)') 'jpctr            = ', jpctr
  WRITE(6,'(A,1X,I5)') 'jpspec           = ', jpspec
  WRITE(6,'(A,1X,I5)') 'jpbk             = ', jpbk
  WRITE(6,'(A,1X,I5)') 'jptk             = ', jptk
  WRITE(6,'(A,1X,I5)') 'jppj             = ', jppj
  WRITE(6,'(A,1X,I5)') 'jphk             = ', jphk 
  WRITE(6,'(A,1X,I5)') 'jpnr             = ', jpnr
  WRITE(6,'(A,1X,I5)') 'jpdd             = ', jpdd
  WRITE(6,'(A,1X,I5)') 'jpdw             = ', jpdw
END IF

IF (lhook) CALL dr_hook('ukca_setup_chem',zhook_out,zhook_handle)

END SUBROUTINE ukca_setup_chem

END MODULE ukca_setup_chem_mod
