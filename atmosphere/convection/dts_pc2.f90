! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates increments to qcl and qcf if PC2 on
!

SUBROUTINE dts_pc2(n_dp,nlev,dts_ntpar,                                       &
                   l_calc_dxek,                                               &
                   timestep,z_theta,zlcl,h_ad,ztop,rho,rho_theta,dr_across_rh,&
                   dr_across_th,temperature,qse,up_flux,q,                    &
                   qcl_env,qcf_env,cf_frozen,cf_liquid,qcl_plume,qcf_plume,   &
                   dcffbydt,dcflbydt,dqcfbydt,dqclbydt,                       &
                   depint,freezint,condint,qclint,qcfint)

! Modules none used


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Calculates increments to qcl and qcf if PC2 being used.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.!
! ------------------------------------------------------------------------------

! Subroutine arguments
 
INTEGER, INTENT(IN) :: &
  n_dp                 &  ! number of points
 ,nlev                    ! Number of levels

INTEGER, INTENT(IN) :: &
  dts_ntpar(n_dp)         ! level at which convection terminates

LOGICAL, INTENT(IN) :: & 
  l_calc_dxek             ! ture if PC2

REAL, INTENT(IN) ::        &
  timestep                 & ! convection timestep (s)
 ,z_theta(n_dp,nlev)       & ! height of theta levels (m)
 ,rho_theta(n_dp,nlev)     & ! density on theta levels
 ,dr_across_rh(n_dp,nlev)  & ! rho layer thickess (m)
 ,dr_across_th(n_dp,nlev)  & ! theta layer thickess (m)
 ,rho(n_dp,nlev)           & ! density on rho levels
 ,zlcl(n_dp)               & ! height of lcl
 ,h_ad(n_dp)               & ! height of maximum buoyancy excess (m)
 ,ztop(n_dp)              & ! height of cloud top (m)
 ,up_flux(n_dp,nlev)       & ! mass flux
 ,q(n_dp,nlev)             & ! environmental water vapour
 ,qcl_env(n_dp,nlev)       & ! Large scale cloud from PC2
 ,qcf_env(n_dp,nlev)       & ! Large scale ice from PC2
 ,qse(n_dp,nlev)           & ! qsat relative to environment
 ,temperature(n_dp,nlev)   & ! temperature (K)         
 ,qcl_plume(n_dp,nlev)     & ! qcl of plume (kg/kg)
 ,qcf_plume(n_dp,nlev)     & ! qcf of plume (kg/kg)      
 ,depint(n_dp)             & ! column integral of deposition rate 
 ,freezint(n_dp)           & ! column integral of freezing rate
 ,condint(n_dp)              ! column integral of condensation rate

            
! Other PC2 prognostics  (qcl and qcf need to be updated if PC2)

REAL, INTENT(IN) ::        &
  cf_frozen(n_dp,nlev)     &  ! Frozen water cloud volume ( )
 ,cf_liquid(n_dp,nlev)        ! Liq water cloud volume

! process fields to give required pc2 outputs
REAL, INTENT(OUT) ::       &
  dcffbydt(n_dp,nlev)      & ! Increments to ice cloud fr  (/s)
 ,dcflbydt(n_dp,nlev)      & ! Increments to liq fr (/s)
 ,dqcfbydt(n_dp,nlev)      & ! Increments to ice condensate (kg/kg/s)
 ,dqclbydt(n_dp,nlev)      & ! Increments to liq condensate (kg/kg/s)
 ,qcfint(n_dp)             & ! Integral of dqcfbydt (takes no account of r**2)
 ,qclint(n_dp)               ! Integral of dqclbydt (takes no account of r**2)    

! Local variables

INTEGER ::        & 
  i_dp,k,kk,i     &  ! loop counters
 ,inew               ! switch

REAL ::                &
  dmdz(n_dp,nlev)      & ! mass flux vertical gradient on theta levels
 ,dmdz_rho(n_dp,nlev)  & ! mass flux gradient on rho levels
 ,entr(n_dp,nlev)      & ! Entrainment rate 
 ,detr(n_dp,nlev)      & ! Detrainment rate
 ,rh(n_dp,nlev)        & ! relative humidity
 ,tmpint               & ! height below which 1g/kg for qcl_plume
 ,temp,factor          &
 ,ent_rate             &
 , det_rate

REAL :: temp1(n_dp)    ! work array for checking qcl and qcf

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
    
!---------------------------------------------------------------------------
! variables for setting up the plume water content
      
  IF (lhook) CALL dr_hook('DTS_PC2',zhook_in,zhook_handle)

  dqclbydt(:,:) = 0.0
  dcflbydt(:,:) = 0.0
  dqcfbydt(:,:) = 0.0
  dcffbydt(:,:) = 0.0
  qcfint(:) = 0.0
  qclint(:) = 0.0

!--------------------------------------------------------------------
! First calculate the mass flux gradient:

  dmdz(:,:) = 0.0
  dmdz_rho(:,:) = 0.0

  DO k=2,nlev
    DO i_dp=1,n_dp
    ! calculate gradient and put on rho levels
      dmdz_rho(i_dp,k) = (up_flux(i_dp,k)-up_flux(i_dp,k-1))/          &
                                                 dr_across_rh(i_dp,k)
    END DO
  END DO

! Added level 1 value 
  k=1
    DO i_dp=1,n_dp
      dmdz_rho(i_dp,k) = up_flux(i_dp,k)/dr_across_rh(i_dp,k)
    END DO


! now interpolate onto theta levels
!        do k=2,nlev-1
! Convection may go from level 1 therefore need a dmdz value? 
  DO k=1,nlev-1
    DO i_dp=1,n_dp
      dmdz(i_dp,k) = (dmdz_rho(i_dp,k)*rho(i_dp,k)               &
                      +dmdz_rho(i_dp,k+1)*rho(i_dp,k+1))/2.      &
                       /rho_theta(i_dp,k)

    END DO
  END DO
  DO k=1,nlev
    DO i_dp=1,n_dp

      rh(i_dp,k) = q(i_dp,k)/qse(i_dp,k)
    END DO
  END DO
!----------------------------------------------------------------------------
  inew = 0

  IF (inew == 1) THEN  ! experimental code option
!----------------------------------------------------------------------------
!  Assume mixing detrainment below  h_ad
!    d = const   small 
! 
!  Above this increase the detrainment to account for forced detrainment?
!    
!----------------------------------------------------------------------------

    det_rate = 1.e-4
    DO k=1,nlev
      DO i_dp=1,n_dp
        IF (up_flux(i_dp,k) > 0.0) THEN  
          ! Not really a measure of detrainment
          IF (z_theta(i_dp,k) <= h_ad(i_dp) ) THEN  
            detr(i_dp,k) = det_rate*up_flux(i_dp,k)
          ELSE  ! plume air becoming cloud detrained as buoyancy falls
            detr(i_dp,k) = up_flux(i_dp,k)*det_rate                      &
                                  *(1.0+2.0*(z_theta(i_dp,k)-h_ad(i_dp)) &
                              /(ztop(i_dp) - h_ad(i_dp))  )

          ENDIF
        ELSE  ! outside convective plume
          detr(i_dp,k) = 0.0
        ENDIF      

      END DO
    END DO

  ELSE    ! Original code - does not work very well 

! Set up a very  simple entrainment rate profile
! nb may want to be more sophisticated in setting this up

    ent_rate = 0.5     ! original
 
    DO k=1,nlev
      DO i_dp=1,n_dp
        entr(i_dp,k) = ent_rate*up_flux(i_dp,k)/z_theta(i_dp,k)
      END DO
    END DO

! Now deduce a detrainment profile based on dmdz = E-D, D = E-dmdz
    DO k=1,nlev
      DO i_dp=1,n_dp
        detr(i_dp,k) = entr(i_dp,k)-dmdz(i_dp,k)
        ! remove any negative detrainment caused by wiggles in mass flux
        IF(detr(i_dp,k) < 0.0) THEN
          detr(i_dp,k) = 0.0
        END IF
      END DO
    END DO

  END IF  ! test on inew


! now put together the liquid cloud fraction and detrained water
! content:
  IF(l_calc_dxek) THEN
    DO k=1,nlev-1
      DO i_dp=1,n_dp
!--------------------------------------------------------------
! Liquid
! Note problems can occur if qcl_plume < qcl_env

! Should it detrain all plume in top layer ? May avoid problems

!       IF (k == dts_ntpar(i_dp)) THEN
!         dqclbydt(i_dp,k) = detr(i_dp,k)*qcl_plume(i_dp,k)       &
!                                                      /rho_theta(i_dp,k)
!       ELSE 
          dqclbydt(i_dp,k) = detr(i_dp,k)*(qcl_plume(i_dp,k) -&
                                      qcl_env(i_dp,k))/rho_theta(i_dp,k)
!       END IF
        dcflbydt(i_dp,k) = detr(i_dp,k)*(1.-cf_liquid(i_dp,k))&
                                            /rho_theta(i_dp,k)
!--------------------------------------------------------------
! Ice
! deduce the same way as liquid:
! Should it detrain all plume in top layer ? May avoid problems
!       IF (k == dts_ntpar(i_dp)) THEN
!         dqcfbydt(i_dp,k) = detr(i_dp,k)*qcf_plume(i_dp,k)         &
!                                                        /rho_theta(i_dp,k)
!       ELSE 
          dqcfbydt(i_dp,k) = detr(i_dp,k)*(qcf_plume(i_dp,k) -      &
                                      qcf_env(i_dp,k))/rho_theta(i_dp,k)
!       END IF
        dcffbydt(i_dp,k) = detr(i_dp,k)*(1.-cf_frozen(i_dp,k))&
                                              /rho_theta(i_dp,k)
!--------------------------------------------------------------
        qcfint(i_dp) = qcfint(i_dp) + dqcfbydt(i_dp,k)              &
                         *rho_theta(i_dp,k)*dr_across_th(i_dp,k)

        qclint(i_dp) = qclint(i_dp) + dqclbydt(i_dp,k)              &
                         *rho_theta(i_dp,k)*dr_across_th(i_dp,k)
      END DO
    END DO
  END IF     ! test on l_calc_dxekl

! LOOPS WRONG WAY AROUND BELOW NEED to sort out?

!---------
! Now check pc2 conversion is less than the deposition and freezing
  DO i_dp=1,n_dp
    IF(depint(i_dp)+freezint(i_dp) < qcfint(i_dp)) THEN 
               
      ! If not, renormalise
      DO k=1,nlev
        dqcfbydt(i_dp,k) = dqcfbydt(i_dp,k)*(depint(i_dp)+freezint(i_dp))  &
                             /qcfint(i_dp)
      END DO
      ! reset pc2 integral as well
      qcfint(i_dp) = depint(i_dp)+freezint(i_dp)
    END IF
  END DO
!--------
! Do same for liquid water:
  DO i_dp=1,n_dp
    tmpint = condint(i_dp)
    IF(tmpint < qclint(i_dp)) THEN 
               
      ! renormalise
      DO k=1,nlev
        dqclbydt(i_dp,k) = dqclbydt(i_dp,k)*(tmpint/qclint(i_dp))
      END DO
      ! reset pc2 integral as well
      qclint(i_dp) = tmpint
    END IF
  END DO


!----------------------------------------------------------------------------
! Checking for negatives after increments applied to env values
! If the environment values are negative after the increments then the 
! model may well blowup.
! Code present to help debug if the model fails
!----------------------------------------------------------------------------
  DO k = 1,nlev
    DO i = 1,n_dp
      IF (dqclbydt(i,k) /= 0.0) THEN
        temp1(i)=qcl_env(i,k) + dqclbydt(i,k) * timestep
        IF (temp1(i)  <   0.0) THEN

         ! Only warning print statement not doing any thing
          write(6,*) ' negative qcl deep',i,k,dts_ntpar(i), &
                       temp1(i),dqclbydt(i,k)
          write(6,*) 'dqcl',(dqclbydt(i,kk),kk=1,nlev)
          write(6,*) 'qcl',(qcl_env(i,kk),kk=1,nlev)

        END IF

      END IF
      IF (dqcfbydt(i,k) /= 0.0) THEN
        temp1(i)=qcf_env(i,k) + dqcfbydt(i,k) * timestep
        IF (temp1(i)  <   0.0) THEN

         ! Only warning print statement not doing any thing
          write(6,*) ' negative qcf deep',i,k,dts_ntpar(i), &
                           temp1(i),dqcfbydt(i,k)
          write(6,*) 'dqcf',(dqcfbydt(i,kk),kk=1,nlev)
          write(6,*) 'qcf',(qcf_env(i,kk),kk=1,nlev)

        END IF
      END IF
    END DO ! n_dp loop
  END DO  ! nlev
  IF (lhook) CALL dr_hook('DTS_PC2',zhook_out,zhook_handle)
  RETURN

!----------------------------------------------------------------------------

END SUBROUTINE dts_pc2
