! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Convert a proportion of fresh soot to aged soot
!-----------------------------------------------------------------------
!+ Perform diffusional scavenging of aged soot to cloud soot
MODULE sootdiffscav_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE sootdiffscav(                                                       &
  ! Arguments IN
  rows, row_length, off_x, off_y, halo_i, halo_j,                              &
  model_levels, wet_model_levels, timestep,                                    &
  cloudf, qcl, qcf, p, t,                                                      &
  n_droplet,                                                                   &
  ! Arguments INOUT
  soot_agd, soot_cld,                                                          &
  ! Arguments OUT
  delta_sootdiffscav                                                           &
  )

! Purpose:
!   To perform diffusional scavenging of aged soot to form
!   the third mode of soot, soot in cloud water.
!
!   Called by Aero_ctl
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
!
! Code Description:
!  Language: Fortran 90.
!  This code is written to UMDP3 v8 programming standards
!
! Documentation: UMDP20

USE earth_constants_mod, ONLY: g 
USE atmos_constants_mod, ONLY: r 

USE water_constants_mod, ONLY: rho_water 

USE conversions_mod, ONLY: pi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE c_st_chm_mod, ONLY: cloudtau, evaptau, nuctau, thold
IMPLICIT NONE

!  includes parameters for rate of soot nucleation/evaporation

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Arguments with intent IN:

INTEGER :: rows                 !no. of rows
INTEGER :: row_length           !no. of pts along a row
INTEGER :: off_x                !size of small halo in i
INTEGER :: off_y                !size of small halo in j
INTEGER :: halo_i               !EW halo size
INTEGER :: halo_j               !NS halo size
INTEGER :: model_levels         !no. of model levels
INTEGER :: wet_model_levels     !no. of wet model levels

REAL :: cloudf(row_length,rows,wet_model_levels)   !Decimal cloud fraction
REAL :: qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                  &
            wet_model_levels)                      !Cloud liquid water
REAL :: qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                  &
            wet_model_levels)                      !Cloud frozen water
REAL :: p(1-off_x:row_length+off_x,1-off_y:rows+off_y,                         &
          model_levels)                            !pressure
REAL :: t(row_length,rows,model_levels)            !temperature
REAL :: timestep                                   !Atmos model timestep
REAL :: n_droplet(row_length,rows,wet_model_levels)!droplet concentration

! Arguments with intent IN:
REAL :: soot_agd(1-off_x:row_length+off_x,1-off_y:rows+off_y,                  &
                 model_levels)  !mmr of aged soot
REAL :: soot_cld(1-off_x:row_length+off_x,1-off_y:rows+off_y,                  &
  model_levels)    !mmr of soot-in-cloud

! Arguments with intent OUT:

! cloud soot increment due to diffusional scavenging
REAL :: delta_sootdiffscav(row_length,rows,model_levels)

! Local variables:

INTEGER  ::  i,j,k  ! loop counters

REAL                                                                           &
  delta_nuc(row_length,rows,model_levels),                                     &
  !              Increment to cloud soot
  delta_evap(row_length,rows,model_levels),                                    &
  !              Increment to aged soot
  qctotal(row_length,rows,wet_model_levels),                                   &
  !              Total cloud water
  clear_frac
!              Clear fraction of grid box (1.0 - cloudf)

!    Variables for improved diffusional scavenging.

REAL :: diffusivity             ! Mean diffusion coefficent of aged soot.
REAL :: viscosity_air           ! Dynamic viscosity of air (kg m-1 s-1).
REAL :: mean_free_path          ! Mean free path of air molecules.
REAL :: knudsen_weber           ! Expression related to the Cunningham slip
                                ! flow correction to the diffusion coefficient.
REAL :: log_sigma_age           ! Natural log of the sigma_age parameter.
REAL :: sq_log_sigma_age        ! Square of the sigma_age parameter.
REAL :: diff_con1               ! Term in diffusion coefficent formula.
REAL :: diff_con2               ! Term in diffusion coefficent formula.
REAL :: pec                     ! Quantity associated with (not equal to) Peclet number.
REAL :: work_radius             ! Variable linked to average droplet radius.
REAL :: scavcd                  ! Scavenging coefficient for in-cloud advective-diffusive 
                                ! removal of aged soot particles.

REAL, PARAMETER :: boltzmann = 1.3804e-23 ! Boltzmanns constant. J K-1
REAL, PARAMETER :: mfp_ref = 6.6e-8       ! Reference value of mean free path. m
REAL, PARAMETER :: tref_mfp = 293.15      ! Reference temperature for mean free path. K
REAL, PARAMETER :: pref_mfp = 1.01325e5   ! Reference pressure for mean free path. Pa
REAL, PARAMETER :: sigma_age = 2.0        ! Geometric standard deviation of the aged soot
                                          ! mode distribution.
REAL, PARAMETER :: rad_age = 4.0e-8       ! Median radius of aged soot distribution, (m).


REAL :: evaptime         ! timescale for cloud droplets to evaporate
REAL :: nuctime          ! timescale for particles to enter a cloud and nucleate.
REAL :: diffuse_tau      ! diffusive lifetime of soot particles once they enter a cloud
REAL :: rho_cuberoot     ! cube root of density of water.
REAL :: diffuse_tauratio ! Cloudtau/Diffuse_tau
REAL :: probdiff_inv     ! inverse of probability of a particle being diffusionallly scavenged in a cloud.
REAL :: probdiff_fn1     ! Probdiff_inv - 0.5
REAL :: probdiff_fn2     ! Probdiff_inv*EXP(diffuse_tauratio*0.5)
REAL :: probdiff_cloud   ! probability of a soot particle being in cloud at the start of a step.
REAL :: probdiff_clear   ! probability of a soot particle being in clear air at the start of a step.
REAL :: lambda_soot      ! ratio of concentrations of soot  particles in cloudy to clear air.
REAL :: diffrate         ! rate of diffusive capture of soot particles

REAL :: term3            ! Local workspace
REAL :: term4            ! Local workspace
REAL :: denom            ! Local workspace

IF (lhook) CALL dr_hook('SOOTDIFFSCAV',zhook_in,zhook_handle)

!     Extra parameters for improved diffusional scavenging.
log_sigma_age = LOG(sigma_age)
sq_log_sigma_age = log_sigma_age * log_sigma_age
diff_con1 = EXP(-2.5*sq_log_sigma_age)/rad_age
diff_con2 = EXP(-4.0*sq_log_sigma_age)/(rad_age*rad_age)

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,evaptime,clear_frac,      &
!$OMP& mean_free_path,knudsen_weber, diffusivity,viscosity_air, pec,    &
!$OMP& work_radius, scavcd, diffuse_tau, diffuse_tauratio, probdiff_inv,&
!$OMP& probdiff_fn1, probdiff_fn2, term3, term4, denom, diffrate,       &
!$OMP& lambda_soot,probdiff_clear,probdiff_cloud)

!-----------------------------------------------------------------------
! 1. Initialise increments to zero
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      delta_nuc(i,j,k)=0.0
      delta_evap(i,j,k)=0.0
    END DO
  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 2. Release soot from evaporating cloud droplets in partly cloudy grid
!    boxes. Also release any soot-in-cloud in cloud-free grid boxes.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(DYNAMIC)
DO k=1,wet_model_levels
  DO j=1,rows
    DO i=1,row_length
      qctotal(i,j,k)=qcl(i,j,k) + qcf(i,j,k)
      IF (qctotal(i,j,k)  <   thold) THEN
        delta_evap(i,j,k) = soot_cld(i,j,k)
        !                      evaporate all the cloud soot in this grid box
      ELSE IF (cloudf(i,j,k)  <   0.95) THEN
        evaptime=evaptau + 0.5*cloudtau
        delta_evap(i,j,k) = (1.0 - EXP(-timestep/evaptime)) *                  &
          soot_cld(i,j,k)
      ELSE
        delta_evap(i,j,k) = 0.0
      END IF
    END DO
  END DO
END DO
!$OMP END DO

!  Also evaporate any soot-in-cloud on non-wet model levels

IF (wet_model_levels  <   model_levels) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=wet_model_levels+1,model_levels
    DO j=1,rows
      DO i=1,row_length
        delta_evap(i,j,k) = soot_cld(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
END IF

!-----------------------------------------------------------------------
! 3. In-cloud scavenging of aged soot particles by cloud droplets.
!    It is assumed that the soot particles do not themselves nucleate
!    droplets. The parametrisation is the same as that used for Aitken
!    mode sulphate aerosol.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(DYNAMIC) 
DO k=1,wet_model_levels
  DO j=1,rows
    DO i=1,row_length

      clear_frac = 1.0 - cloudf(i,j,k)

      IF ((qctotal(i,j,k)  >=  thold) .AND.                                    &
        (cloudf(i,j,k)  >  0.0)) THEN

        !     First compute in-cloud timescale for diffusional capture,
        !     using total condensed water within the cloudy portion of the
        !     grid box. The difference between liquid water and ice is
        !     neglected here. This should be improved on eventually.

        !     Compute mean free path of air molecules.
        !     (See P.417 of Pruppacher and Klett, 2nd edition.)

        mean_free_path = (mfp_ref*pref_mfp*t(i,j,k))/                          &
          (tref_mfp*p(i,j,k))

        !     Compute the Knudsen and Weber term for a particle of the median
        !     radius (note approximation here: we do not average over the size
        !     distribution). See P.450 of Pruppacher and Klett, 2nd edition.

        knudsen_weber = 1.257 + 0.4*                                           &
          EXP(- ((1.1*rad_age)/mean_free_path) )

        !     Temporarily use DIFFUSIVITY to store working value.

        diffusivity = diff_con1 +                                              &
          diff_con2*mean_free_path*knudsen_weber

        !     Compute dynamic viscosity of air, using an approximate version of
        !     the formula on P.417 of Pruppacher and Klett, 2nd edition.

        viscosity_air = ( 1.718 + (t(i,j,k) - 273.15)*0.0049 )*                &
          1.0e-5

        !     Now compute mean diffusion coefficient.

        diffusivity = (boltzmann*t(i,j,k)*diffusivity)/                        &
          (6.0*pi*viscosity_air)

        !     Now compute the term PEC related to (but not equal to) the cube
        !     root of the Peclet Number.

        pec= ((4.0*g*rho_water)/(9.0*diffusivity*viscosity_air))               &
          **0.333333

        work_radius = qctotal(i,j,k)/                                          &
          (cloudf(i,j,k)*10.0*pi*rho_water*n_droplet(i,j,k))
        work_radius = work_radius**0.333333

        !     We can finally compute the timescale for diffusive
        !     scavenging once inside cloud, DIFFUSE_TAU.

        scavcd = 6.0*pi*diffusivity*n_droplet(i,j,k)*work_radius*              &
          (1.0 + pec*work_radius)
        diffuse_tau = 1.0/scavcd

        diffuse_tauratio = cloudtau/diffuse_tau
        probdiff_inv = 1.0/( 1.0 - EXP(-diffuse_tauratio) )
        probdiff_fn1 = probdiff_inv - 0.5
        probdiff_fn2 = probdiff_inv*EXP(-(0.5*diffuse_tauratio) )
        !   CALCULATE LAMBDA_SOOT.

        term3 = (clear_frac*diffuse_tauratio)**2
        term3 = term3 + ( 2.0*diffuse_tauratio                                 &
          *clear_frac*(clear_frac-cloudf(i,j,k)))
        term3 = SQRT(1.0+term3)
        term4=2.0*cloudf(i,j,k)-                                               &
          diffuse_tauratio*clear_frac-1.0
        term4 = term4 + term3
        lambda_soot = term4/(2.0*cloudf(i,j,k))

        !   CALCULATE PROBDIFF_CLEAR AND PROBDIFF_CLOUD

        denom = clear_frac+cloudf(i,j,k)*lambda_soot
        probdiff_clear = clear_frac/denom
        probdiff_cloud = (cloudf(i,j,k)*lambda_soot)/denom

        !   CALCULATE EXPECTED LIFETIME OF A SOOT PARTICLE WITH
        !   RESPECT TO DIFFUSIVE CAPTURE BY CLOUD DROPLETS.

        term3 = probdiff_fn1*probdiff_clear
        term4 = probdiff_fn2*probdiff_cloud
        term4 = (term3+term4)*clear_frac/cloudf(i,j,k)
        denom = term4*cloudtau + diffuse_tau
        diffrate = 1.0/denom

        !   Now compute the amount of aged soot converted to cloud soot
        !   in this timestep.
        delta_nuc(i,j,k)=(1.0-EXP(-diffrate*timestep))*                        &
          soot_agd(i,j,k)

      END IF ! Test on qctotal and thold
    END DO
  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 4. Calculate total increment for output
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      delta_sootdiffscav(i,j,k) =                                              &
        delta_nuc(i,j,k) - delta_evap(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO

!$OMP END PARALLEL 

IF (lhook) CALL dr_hook('SOOTDIFFSCAV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sootdiffscav
END MODULE sootdiffscav_mod
