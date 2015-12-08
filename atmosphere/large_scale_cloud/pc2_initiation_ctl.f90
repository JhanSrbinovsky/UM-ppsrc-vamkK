! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initiate cloud and liquid water within the PC2 Cloud Scheme.

SUBROUTINE pc2_initiation_ctl (                                         &
! Dimensions of Rh crit array
  rhc_row_length, rhc_rows,                                             &

! Model switches
  ltimer, l_mixing_ratio,                                               &
  l_acf_cusack, l_cld_area,                                             &
! in time stepping information.
  timestep,                                                             &

! SCM diagnostics switches
  nSCMDpkgs,L_SCMDiags,                                                 &

! Primary fields passed IN/OUT
  t,q,qcl,qcf,cf,cfl,cff,rhts,tlts,qtts,ptts,cf_area,                   &

! Primary fields passed IN
  p,pstar,p_theta_levels,ccb,cumulus,rhcrit,                            &

! Output increments
  t_work,q_work,qcl_work,qcf_work,cf_work,cfl_work,cff_work             &
 )

  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: pdims, qdims, tdims, pdims_s

  IMPLICIT NONE

! Description:
!   Initiate a small amount of cloud fraction and liquid
!   water content for the PC2 Cloud Scheme. Check that moisture
!   variables are consistent with each other.
!
! Method:
!   See the PC2 documentation.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Code Description:
!   Language: fortran 77 + common extensions
!   This code is written to UMDP3 v6 programming standards.

! Declarations:
! Arguments with intent in. ie: input variables.

! Model dimensions
  INTEGER ::         &
    rhc_row_length,  & ! Dimensions of RHcrit variable
    rhc_rows           ! Dimensions of RHcrit variable
!
  LOGICAL ::                            &
    ltimer,                             & ! Output timing information.
    cumulus(qdims%i_start:qdims%i_end,  & ! 
            qdims%j_start:qdims%j_end), & ! Convection is occurring.
    l_cld_area,                         & ! Switch for area cloud 
                                          ! fraction(ACF)parametrization.
    l_acf_cusack,                       & ! ... to select Cusack ACF.
    l_mixing_ratio                        ! Use mixing ratio formulation

! time information for current timestep
  REAL ::                                                               &
    timestep

! Primary fields passed in/out
  REAL ::                                                               &
    t(                 tdims%i_start:tdims%i_end,                       & 
                       tdims%j_start:tdims%j_end,                       &  
                                   1:tdims%k_end),                      &
    q(                 qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    qcl(               qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    qcf(               qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    cf(                qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    cfl(               qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    cff(               qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    cf_area(           qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end)

! Primary fields passed in
  REAL ::                                                               &
    p(                 pdims_s%i_start:pdims_s%i_end,                   & 
                       pdims_s%j_start:pdims_s%j_end,                   &  
                       pdims_s%k_start:pdims_s%k_end+1),                &
    pstar(             pdims%i_start  :pdims%i_end,                     & 
                       pdims%j_start  :pdims%j_end),                    &
    p_theta_levels(    pdims%i_start:pdims%i_end,                       & 
                       pdims%j_start:pdims%j_end,                       &  
                       pdims%k_start:pdims%k_end),                      &
    rhcrit(            rhc_row_length,                                  &
                       rhc_rows,                                        &
                                   1:qdims%k_end),                      &
    rhts(              qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    tlts(              qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
!       TL at start of timestep
    qtts(              qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
!       qT at start of timestep
    ptts(              pdims%i_start:pdims%i_end,                       & 
                       pdims%j_start:pdims%j_end,                       &  
                       pdims%k_start:pdims%k_end)
!       Pressure at theta levels at start of timestep

  INTEGER ::                                                            &
    ccb(               qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end)

  INTEGER ::                                                            &
    nSCMDpkgs             ! No of SCM diagnostics packages

  LOGICAL ::                                                            &
    L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

! Local variables


! Output increment diagnostics
  REAL ::                                                               &
    t_work(            tdims%i_start:tdims%i_end,                       & 
                       tdims%j_start:tdims%j_end,                       &  
                                   1:tdims%k_end),                      &
    q_work(            qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    qcl_work(          qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    qcf_work(          qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    cf_work(           qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    cfl_work(          qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    cff_work(          qdims%i_start:qdims%i_end,                       & 
                       qdims%j_start:qdims%j_end,                       &  
                                   1:qdims%k_end),                      &
    p_layer_boundaries(pdims%i_start:pdims%i_end,                       & 
                       pdims%j_start:pdims%j_end,                       &  
                                   0:pdims%k_end),                      &
!       Pressure at layer boundaries. Same as p except at
!       bottom level = pstar, and at top = 0.
    p_layer_centres(   pdims%i_start:pdims%i_end,                       & 
                       pdims%j_start:pdims%j_end,                       &  
                                   0:pdims%k_end)
!       Pressure at layer centres. Same as p_theta_levels
!       except bottom level = pstar, and at top = 0.

  INTEGER ::                                                            &
    i,j,k,                                                              &
!       Loop counters
    large_levels,                                                       &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((wet_model_levels - 2)*levels_per_level) + 2
    levels_per_level
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops

  CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'pc2_initiation_ctl'

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

! External Functions:

!- End of header


  IF (lhook) CALL dr_hook('PC2_INITIATION_CTL',zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end

! Work fields are set to starting fields so diagnostics
! can be calculated

        q_work(i,j,k)   = q(i,j,k)
        qcl_work(i,j,k) = qcl(i,j,k)
        qcf_work(i,j,k) = qcf(i,j,k)
        cf_work(i,j,k)  = cf(i,j,k)
        cfl_work(i,j,k) = cfl(i,j,k)
        cff_work(i,j,k) = cff(i,j,k)

      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t_work(i,j,k)   = t(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL 

! Call checking routine

! DEPENDS ON: pc2_checks
  CALL pc2_checks(p_theta_levels,                                       &
      t, cf, cfl, cff, q, qcl, qcf,                                     &
      l_mixing_ratio)


! Call initiation routine

! using area cloud parameterisation
  IF (l_cld_area) THEN
! use Cusack interpolation option
    IF (l_acf_cusack) THEN

! set p at layer boundaries.
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          p_layer_boundaries(i,j,0) = pstar(i,j)
          p_layer_centres(i,j,0) = pstar(i,j)
        END DO
      END DO
      DO k = pdims%k_start, pdims%k_end - 1
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            p_layer_boundaries(i,j,k) = p(i,j,k+1)
            p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
          END DO
        END DO
      END DO
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          p_layer_boundaries(i,j,pdims%k_end) = 0.0
          p_layer_centres(i,j,pdims%k_end) =                           &
                           p_theta_levels(i,j,pdims%k_end)
        END DO
      END DO

! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
      levels_per_level = 3
      large_levels = ((pdims%k_end - 2)*levels_per_level) + 2

! DEPENDS ON: pc2_arcld
      CALL pc2_arcld(p_layer_centres,p_layer_boundaries,                &
        ccb,cumulus,rhcrit,                                             &
        rhc_row_length,rhc_rows,                                        &
        large_levels,levels_per_level,cf_area,                          &
        t,cf,cfl,cff,q,qcl,qcf,rhts,tlts,qtts,ptts,l_mixing_ratio)

    END IF !L_ACF_Cusack

  ELSE !l_cld_area

! DEPENDS ON: pc2_initiate
    CALL pc2_initiate(p_theta_levels,ccb,cumulus,rhcrit,                &
      qdims%k_end, rhc_row_length,rhc_rows,                             &
      t,cf,cfl,cff,q,qcl,rhts,l_mixing_ratio)

  END IF !l_cld_area


! Call second checking routine

! DEPENDS ON: pc2_checks2
  CALL pc2_checks2(p_theta_levels,rhcrit,                               &
      rhc_row_length,rhc_rows,                                          &
      t, cf, cfl, cff, q, qcl, l_mixing_ratio)

! Call first checking routine again

! DEPENDS ON: pc2_checks
  CALL pc2_checks(p_theta_levels,                                       &
      t, cf, cfl, cff, q, qcl, qcf,                                     &
      l_mixing_ratio)

! using area cloud parameterisation
  IF (l_cld_area) THEN
! use Cusack interpolation option
    IF (l_acf_cusack) THEN
! DEPENDS ON: pc2_hom_arcld
      CALL pc2_hom_arcld(p_layer_centres,p_layer_boundaries,            &
         large_levels,levels_per_level,                                 &
         cf_area,t,cf,cfl,cff,q,qcl,qcf,                                &
         l_mixing_ratio)
    END IF  ! L_ACF_Cusack
  END IF    ! L_cld_area

! Update work array to hold net increment from the above routines

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        q_work(i,j,k)   = q(i,j,k)   - q_work(i,j,k)
        qcl_work(i,j,k) = qcl(i,j,k) - qcl_work(i,j,k)
        qcf_work(i,j,k) = qcf(i,j,k) - qcf_work(i,j,k)
        cf_work(i,j,k)  = cf(i,j,k)  - cf_work(i,j,k)
        cfl_work(i,j,k) = cfl(i,j,k) - cfl_work(i,j,k)
        cff_work(i,j,k) = cff(i,j,k) - cff_work(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t_work(i,j,k)   = t(i,j,k)   - t_work(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL 


! End of routine initiation_ctl

  IF (lhook) CALL dr_hook('PC2_INITIATION_CTL',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_initiation_ctl
