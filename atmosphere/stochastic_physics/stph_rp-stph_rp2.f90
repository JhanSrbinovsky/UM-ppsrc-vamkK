! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Stochastic Physics (sect35) Random Parameters Ver. 2
MODULE stph_rp2_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE stph_rp2(model_levels,                                       &
                    rhcrit, rhcrit_max, rhcrit_min,                     &
                    gwd_frc, gwd_frc_max, gwd_frc_min,                  &
                    kay_gwave, kay_gwave_max, kay_gwave_min,            &
                    m_ci, m_ci_max, m_ci_min,                           &
                    Charnock)
     
 USE cv_run_mod, ONLY:                                                  &
     cape_timescale, i_convection_vn,i_convection_vn_4a     
     

 USE g_wave_input_mod, ONLY : i_gwd_vn,             &
                              i_gwd_vn_4a

 USE stochastic_physics_run_mod, ONLY:                                  &
     ran_max, rhcrit_ref_level, ran_count,                              &
     cape_timescale_min, cape_timescale_max,                            &
     entcoef_min, entcoef_max, g0_rp, g0_rp_max, g0_rp_min,             &
     par_mezcla, par_mezcla_max, par_mezcla_min,                        &
     charnock_max, charnock_min,                                        &
     lambda_min_rp_max, lambda_min_rp, lambda_min_rp_min,               &
     ricrit_rp, ricrit_rp_max, ricrit_rp_min,                           &
     a_ent_1_rp, a_ent_1_rp_max, a_ent_1_rp_min,                        &
     g1_rp, g1_rp_max, g1_rp_min
     
 USE entcoef_mod, ONLY: entcoef    

 USE yomhook, ONLY: lhook, dr_hook           
 USE parkind1, ONLY: jprb, jpim

 USE PrintStatus_mod
 USE UM_ParVars

 USE domain_params
 USE stph_rp_pert_mod, ONLY: stph_rp_pert
 IMPLICIT NONE

!
! Description:  Introduce stochastic perturbation in some physics params
!               to take into account model error
!
! Method:       The value of a given PARAMETER is chosen randomly
!               between a maximum and minimum values. PARAMETER's values
!               are temporally correlated (1st order markov process)
!
! Code Owner:   See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
!
! Code Description:
!   Language:   FORTRAN 90
!   This code is written to UMDP3 version 8.3 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER, INTENT(IN) :: ExampleVariable  !Description of variable
!
! Global variables (#include statements etc):

!-------------------------------------------------------------
! Variable definition
!-------------------------------------------------------------
!
! IN variables
!
 INTEGER, INTENT(IN) :: model_levels
!
! Default, Maximum and minimum values for the STPH_RP scheme
! Large Scale Precipitation
! Note: rhcrit scalars refer to reference level=rhcrit_ref_level
!
 REAL, INTENT(IN) :: rhcrit_max        ! Max value critical rh
 REAL, INTENT(IN) :: rhcrit_min        ! Min value critical rh
 REAL, SAVE       :: rhcrit_0          ! Def value critical rh
 REAL, SAVE       :: rhcrit_p          ! Perturbed critical rh
 REAL, INTENT(IN) :: m_ci_max          ! Max value of multiplication
                                       ! factor for CI  
 REAL, INTENT(IN) :: m_ci_min          ! Min value of multiplication
                                       ! factor for CI                             
 REAL, SAVE       :: m_ci_0            ! Def value of multiplication
                                       ! factor for CI                             
!
! Gravity Wave drag
!
 REAL, INTENT(IN) :: gwd_frc_max       ! Max value critical Froude number
 REAL, INTENT(IN) :: gwd_frc_min       ! Min value critical Froude number
 REAL, SAVE       :: gwd_frc_0         ! Def value critical Froude number
 REAL, INTENT(IN) :: kay_gwave_max     ! Max value gravity wave parameter
 REAL, INTENT(IN) :: kay_gwave_min     ! Min value gravity wave parameter
 REAL, SAVE       :: kay_gwave_0       ! Def value gravity wave parameter
!
! Convection
!
 REAL, SAVE       :: entcoef_0         ! Def Entrainment rate coef
 REAL, SAVE       :: cape_timescale_0  ! Def cape closure timescale
!
! Boundary Layer
!
 REAL, SAVE       :: par_mezcla_0      ! Def value neutral mixing length
                                       ! parameter
 REAL, SAVE       :: g0_rp_0           ! Def value of stability functions
 REAL, SAVE       :: lambda_min_rp_0   ! Def value of min mixing length
 REAL, SAVE       :: ricrit_rp_0       ! Def value of critical Ri 
 REAL, SAVE       :: a_ent_1_rp_0      ! Def value of entrainment A1 
 REAL, SAVE       :: g1_rp_0           ! Def value of velocity scale
 REAL, SAVE       :: charnock_0        ! Def value Charnock param.                             

!                                  
! IN/OUT variables
!
 REAL, INTENT(INOUT) :: kay_gwave      ! Surface stress constant GWD
 REAL, INTENT(INOUT) :: rhcrit(model_levels)
                                       ! Crit relative humidity for
                                       ! layer cloud formation
 REAL, INTENT(INOUT) :: gwd_frc        ! Critical Froude Number 
 REAL, INTENT(INOUT) :: m_ci           ! Fall-speed in 3C/D microphysics
 REAL, INTENT(INOUT) :: Charnock       ! Charnock Parameter


! LOCAL VARIABLES
!
! Local variables to assign values to rhcrit, which is level
! dependant
 INTEGER, ALLOCATABLE,SAVE :: interval(:)
 REAL,    ALLOCATABLE,SAVE :: drhcrit(:)
                                       ! used to store differences
                                       ! in rhcrit between model levels

! Other local variables
 INTEGER, SAVE   :: max_interval
 LOGICAL, SAVE   :: l_firstcall = .TRUE.
                                       ! 1st call to RP?
 LOGICAL         :: local_switch       ! Logical to calculate the first
                                       ! new rhcrit value.
 INTEGER         :: i,j                ! loop variables

! Variables associated with random number generator
 REAL            :: rp_rand(ran_max)   ! random number for each parameter
 REAL            :: rp_rand0(ran_max)  ! random number for initial values
 INTEGER         :: istat              ! To record call in the synch.

! Define mean values
 REAL, PARAMETER :: entcoefmean=3.0    ! entrainment rate coef mean

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_RP2',zhook_in,zhook_handle)

 IF (l_firstcall) THEN
   ! ---------------------------------------
   ! Allocate interval the first time STPH_RP is called only
   ! This will be used when perturbing rhcrit.
   ! ---------------------------------------
   IF (.NOT.ALLOCATED(interval)) THEN
     ALLOCATE (interval(model_levels))
     interval(:) = 0
   END IF
   IF (.NOT.ALLOCATED(drhcrit)) THEN
     ALLOCATE (drhcrit(model_levels))
     drhcrit(:) = 0
   END IF
   j=1
   DO i = rhcrit_ref_level, model_levels
     IF(rhcrit(i) <  rhcrit(i-1)) THEN          
       drhcrit(i) = rhcrit(i-1) - rhcrit(i)  
       interval(i) = j
       j = j + 1
     END IF
   ENDDO
   max_interval = j-1
!
   DO i = rhcrit_ref_level + 1, model_levels
     IF(interval(i) <  interval(i-1)) THEN
          interval(i)=interval(i-1)
     END IF
     IF(drhcrit(i) <  drhcrit(i-1)) THEN
          drhcrit(i)=drhcrit(i-1)
     END IF
   ENDDO
   rhcrit_0 = rhcrit(rhcrit_ref_level)       ! critical RH mean value 

! Setup random number for initial perturbed values
   IF(mype == 0) THEN
     CALL random_number(rp_rand0)
   END IF
   CALL gc_rbcast(3145, ran_max, 0, nproc, istat, rp_rand0)
 END IF
!----------------------------------------------------------------------
! Each time random parameters is called broadcast the random number
! to all processors
!----------------------------------------------------------------------
 IF(mype == 0) THEN
   CALL random_number(rp_rand)
 END IF
 CALL gc_rbcast(3145, ran_max, 0, nproc, istat, rp_rand)
 IF (printstatus  >  prstatus_normal) THEN
   WRITE(6,*)
   WRITE(6,'("*** STPH_RP2 NEW PARAM VALUES ***")')
 END IF
! Start using 1st random number from array rp_rand
 ran_count = 1

!
! Convection
!
IF (i_convection_vn == i_convection_vn_4a) THEN
! Perturb entrainment rate coefficient (Convection)
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, entcoef_0,          &
                    entcoef_max, entcoef_min, entcoef)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("ENTCOEF .......... ",2ES16.8)') entcoef, entcoef_0

! Perturb cape_timescale closure (Convection)
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, cape_timescale_0,   &
           cape_timescale_max, cape_timescale_min, cape_timescale)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("CAPE_CLOSURE...... ",2ES16.8)') cape_timescale,           &
                                              cape_timescale_0
END IF

!
! Large-scale Precip
!
! Perturb m_ci  (Ice fall speed) - 3C/D Microphysics LSPCON3C
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, m_ci_0,             &
                    m_ci_max, m_ci_min, m_ci)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("M_CI ............. ",2ES16.8)') m_ci, m_ci_0

!
! Boundary Layer
!

! Perturb par_mezcla (neutral mixing length) - EXCOEF (PBL)
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, par_mezcla_0,       &
                    par_mezcla_max, par_mezcla_min, par_mezcla)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("PAR_MEZCLA ....... ",2ES16.8)') par_mezcla, par_mezcla_0

! Perturb g0_rp (stability function parameter) - EXCOEF (PBL)
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, g0_rp_0,            &
                    g0_rp_max, g0_rp_min, g0_rp)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("G0_RP ............ ",2ES16.8)') g0_rp, g0_rp_0

! Perturb lambda_min_rp (minimum mixing length) - (PBL) 
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, lambda_min_rp_0,    &
                    lambda_min_rp_max, lambda_min_rp_min, lambda_min_rp)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("LAMBDA_MIN_RP..... ",2ES16.8)') lambda_min_rp,            &
                                              lambda_min_rp_0

! Perturb ricrit_rp (critical Ri) - (PBL) 
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, ricrit_rp_0,        &
                    ricrit_rp_max, ricrit_rp_min, ricrit_rp)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("RICRIT_RP ........ ",2ES16.8)') ricrit_rp, ricrit_rp_0

! Perturb a_ent_1_rp (entrainment parameter A1) - (PBL) 
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, a_ent_1_rp_0,       &
                    a_ent_1_rp_max, a_ent_1_rp_min, a_ent_1_rp)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("A_ENT_1_RP ....... ",2ES16.8)') a_ent_1_rp, a_ent_1_rp_0

! Perturb g1_rp (velocity scale parameter g1) - (PBL) 
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, g1_rp_0,            &
                    g1_rp_max, g1_rp_min, g1_rp)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("G1_RP ............ ",2ES16.8)') g1_rp, g1_rp_0

! Perturb charnock (Charnock parameter) - (PBL)
 CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, charnock_0,         &
                    charnock_max, charnock_min, charnock)
 IF (printstatus  >  prstatus_normal)                                   &
   WRITE(6,'("CHARNOCK ......... ",2ES16.8)') charnock, charnock_0


!
! Gravity-wave drag
!
 IF ( i_gwd_vn == i_gwd_vn_4a ) THEN
!   Perturb gwd_frc (critical Froude number) - (GWD)
    CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, gwd_frc_0,          &
                       gwd_frc_max, gwd_frc_min, gwd_frc)
    IF (printstatus  >  prstatus_normal)                                   &
      WRITE(6,'("GWD_FRC .......... ",2ES16.8)') gwd_frc, gwd_frc_0

!   Perturb kay_gwave (surface stress constant) - (GWD)
    CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, kay_gwave_0,        &
                       kay_gwave_max, kay_gwave_min, kay_gwave)
    IF (printstatus  >  prstatus_normal)                                   &
      WRITE(6,'("KAY_GWAVE ........ ",2ES16.8)') kay_gwave, kay_gwave_0
 END IF ! i_gwd_vn

! Perturb rhcrit (critical relative humidity) - (LSP)
 IF (l_firstcall) rhcrit_p = rp_rand0(ran_count) *                      &
                             (rhcrit_max-rhcrit_min) + rhcrit_min
!
 local_switch = .TRUE. 
 DO j=1,max_interval
   DO i=1,model_levels
     IF (j == 1) THEN
       IF (interval(i) == j) THEN
         IF (local_switch) THEN 
           CALL stph_rp_pert( .FALSE., rp_rand, rp_rand, rhcrit_0,      &
                              rhcrit_max, rhcrit_min, rhcrit_p)
           rhcrit(i) = rhcrit_p
           local_switch = .FALSE.
         ELSE
           rhcrit(i) = rhcrit(i-1)
         END IF 
       END IF
     ELSE
       IF (interval(i) == j) THEN
         rhcrit(i) = rhcrit_p - (drhcrit(i)*(j-1))
       END IF
     END IF
   ENDDO
 ENDDO
 IF (printstatus  >  prstatus_normal) THEN
   WRITE(6,'("RHCRIT(4,7,10) ... ",4ES16.8)') rhcrit(4), rhcrit(7),     &
                                          rhcrit(10), rhcrit_0
   WRITE(6,'("RHCRIT(all) ... ",8ES16.8)')                              &
                                 (rhcrit(i), i = 1, model_levels)
 END IF

! Terminate firstcall condition
 l_firstcall=.FALSE.
 
 IF (lhook) CALL dr_hook('STPH_RP2',zhook_out,zhook_handle)
 RETURN
END SUBROUTINE stph_rp2
END MODULE stph_rp2_mod
