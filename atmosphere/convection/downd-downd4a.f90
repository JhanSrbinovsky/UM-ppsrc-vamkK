! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calls downdraught calculation or checks for change of phase
!
SUBROUTINE downd_4a5a (npnts,np_full,k,kct,ntra,nddon_a, iccb,    &
                  l_tracer,bwater_k,                              &
                  timestep,the_k,the_km1,qe_k,qe_km1,qse_km1,     &
                  p_km1,delpk,                                    &
                  delpkm1,exk,exkm1,deltd,delqd,amdetk,ekm14,     &
                  ekm34,flx_ud_k,cca, trae_k,trae_km1,deltrad,    &
                  bdd_start,bddwt_k,bddwt_km1,bdd_on,             &
                  thdd_k,qdd_k,dthbydt_k,dthbydt_km1,dqbydt_k,    &
                  dqbydt_km1,rain,snow,precip_k,rain_env,snow_env,&
                  rain_dd,snow_dd,flx_dd_k,lr_ud_ref,             &
                  tradd_k,dtrabydt_k,dtrabydt_km1)

USE earth_constants_mod, ONLY: g

USE cv_run_mod, ONLY:                                             &
    i_convection_vn,                                              &
    i_convection_vn_4a,                                           &
    i_convection_vn_5a,                                           &
    dd_opt
USE cv_param_mod, ONLY:                                           &
    ddptef

USE ereport_mod, ONLY : Ereport

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
! ------------------------------------------------------------------------------
! Description:
! Calls downdraught calculation.
! Change of phase calculation where no downdraught
!
!  See UM Documentation paper No. 27
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,np_full              & ! Full vector length
 ,k                    & ! Present model layer
 ,kct                  & ! Convective cloud top layer 
 ,ntra                 & ! Number of tracers
 ,nddon_a                ! Number of points at which downdraught does occur

INTEGER, INTENT(IN) :: &
  iccb(npnts)            ! Cloud base model level

LOGICAL, INTENT(IN) :: &
  l_tracer               ! Switch for tracers

LOGICAL, INTENT(IN) :: &
  bwater_k(npnts)        ! Mask for points at which condensate is water
                         ! in layer k

REAL, INTENT(IN) ::    &
  timestep               ! timestep

REAL, INTENT(IN) :: &
  the_k(npnts)      & ! potential temperature of environment in layer k (K)
 ,the_km1(npnts)    & ! Potential temperature of environment in layer k-1 (K)
 ,qe_k(npnts)       & ! environment mixing ratio of layer k (kg/kg)
 ,qe_km1(npnts)     & ! environment mixing ratio of layer k-1 (kg/kg)
 ,qse_km1(npnts)    & ! environment qsat mixing ratio of layer k-1 (kg/kg)
 ,p_km1(npnts)      & ! Pressure of layer k-1 (Pa)
 ,delpk(npnts)      & ! Pressure difference across layer K  (Pa)
 ,delpkm1(npnts)    & ! Pressure difference across layer K-1 (Pa)
 ,exk(npnts)        & ! exner ratio for layer K
 ,exkm1(npnts)      & ! exner ratio for layer K-1
 ,deltd(npnts)      & ! Cooling necessary to achieve saturation (K)
 ,delqd(npnts)      & ! moistening necessary to achieve saturation (kg/kg)
 ,amdetk(npnts)     & ! Mixing detrainment at level k multiplied by 
                      ! appropriate layer thickness
 ,ekm14(npnts)      & ! exner ratio at layer k-1/4
 ,ekm34(npnts)      & ! exner ratio at layer k-3/4
 ,flx_ud_k(npnts)   & ! updraught mass flux at layer K
 ,cca(npnts)          ! convective cloud amount

REAL, INTENT(IN) ::    &
  trae_k(npnts,ntra)   & ! tracer of environment in layer k (kg/kg)
 ,trae_km1(npnts,ntra) & ! tracer of environment in layer k-1 (kg/kg)
 ,deltrad(npnts,ntra)    ! Depletion of environment tracer due to 
                         ! downdraught formation (kg/kg)


LOGICAL, INTENT(INOUT) :: &
  bdd_start(npnts)        & !IN Mask for those points where downdraught may
                            !   form in layer k
                            !OUT Mask for those points where downdraught may
                            !   form in layer k-1
 ,bddwt_k(npnts)          & !   Mask for those points in downdraught where
                            !   precipitation is liquid in layer k
 ,bddwt_km1(npnts)        & !   Mask for those points in downdraught where
                            !   precipitation is liquid in layer k-1
 ,bdd_on(npnts)             !IN Mask for those points where DD has continued
                            !   from previous layer 
                            !OUT Mask for those points where downdraught  
                            !    continues to layer k-1

REAL, INTENT(INOUT) :: & 
  thdd_k(npnts)        & ! IN Model potential temperature of downdraught 
                         !   in layer k (K)
                         !OUT Updated potential temperature of downdraught 
                         !   in layer k (K)
 ,qdd_k(npnts)         & ! IN mixing ratio of downdraught in layer k (kg/kg)
                         !OUT Updated mixing ratio of downdraught in layer k
 ,dthbydt_k(npnts)     & ! IN Increment to model potential temperature of 
                         !    layer k (K/s)
                         !OUT Updated increment to model potential temperature
                         !    of layer k (K/s)
 ,dthbydt_km1(npnts)   & ! IN Increment to model potential temperature of 
                         !    layer k-1 (K/s)
                         !OUT Updated increment to model potential temperature
                         !    of layer-1 (K/s)
 ,dqbydt_k(npnts)      & ! IN Increment to model mixing ratio of layer k
                         !    (kg/kg/s)
                         !OUT Updated increment to  mixing ratio of layer k
                         !    (kg/kg/s)
 ,dqbydt_km1(npnts)    & ! IN Increment to model mixing ratio of layer k-1
                         !    (kg/kg/s)
                         !OUT Updated increment to  mixing ratio of layer k-1
                         !    (kg/kg/s)
 ,rain_env(npnts)      & ! Amount of rainfall passing through environment 
                         ! (kg/m**2/s)
 ,snow_env(npnts)      & ! Amount of snowfall passing through environment 
                         ! (kg/m**2/s)
 ,rain_dd(npnts)       & ! Amount of rainfall passing through
                         ! downdraught (KG/M**2/S)
 ,snow_dd(npnts)       & ! Amount of snowfall passing through
                         ! downdraught (KG/M**2/S)
 ,rain(npnts)          & ! IN Initialised rainfall (kg/m**2/s)
                         ! OUT Surface rainfall (kg/m**2/s)
 ,snow(npnts)          & ! IN Initialised snowfall (kg/m**2/s)
                         ! OUT Surface snowfall (kg/m**2/s)
 ,precip_k(npnts)      & ! Precipitation added when descending from layer
                         ! k to k-1 (kg/m**2/s)
 ,flx_dd_k(npnts)      & ! Downdraught initial mass flux (Pa/s)
 ,lr_ud_ref(npnts)       ! precipitation mixing ratio at lowest 
                         ! precipitating level of UD
 
REAL, INTENT(INOUT) ::      &
  tradd_k(npnts,ntra)       & ! Model tracer of downdraught in layer k (kg/kg)
 ,dtrabydt_k(np_full,ntra)  & ! IN Increment to model tracer of layer k
                              !    (kg/kg/s)
                              ! OUT Updated increment to model tracer of layer k
                              !    (kg/kg/s)
 ,dtrabydt_km1(np_full,ntra)  ! IN Increment to model tracer of layer k-1
                              !    (kg/kg/s)
                              ! OUT Updated increment to model tracer of 
                              !    layer k-1 (kg/kg/s)

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::        & 
  i,ktra          & ! loop counter
 ,nddon             ! Number of points at which downdraught does occur
                 
INTEGER ::        &
  index1(nddon_a)   !  Index for compressing points

LOGICAL ::          &
  bwork(nddon_a,5)  & !  Work space for 'bit' masks
 ,b_dd_end(npnts)     !  mask for points where downdraught has ended


REAL ::                &
  work(nddon_a,38)      ! Work space  for compression

REAL ::                &
  qse_km1_c(nddon_a)      ! qsat k-1 for compression

REAL ::                    &
  tradd_k_c(nddon_a,ntra)  & ! Tracer content in downdraught at layer k
                             ! compressed (kg/kg)
 ,trae_k_c(nddon_a,ntra)   & ! Tracer content of environment layer k
                             ! compressed (kg/kg) 
 ,trae_km1_c(nddon_a,ntra) & ! Tracer content of environment layer k-1
                             ! compressed (kg/kg) 
 ,dtra_k_c(nddon_a,ntra)   & ! Increment to model tracer in layer k
                             ! compressed (kg/kg) 
 ,dtra_km1_c(nddon_a,ntra) & ! Increment to model tracer in layer k-1
                             ! compressed (kg/kg) 
 ,deltrad_c(nddon_a,ntra)    ! Depletion of environmen  tracer due to 
                             ! downdraught formation compressed

REAL ::                   &
  ppn_mix_dd_c(nddon_a)     ! Precip mixing ratio in DD

REAL ::             &
  factor            & ! Proportion of rainfall going into downdraught from UD
 ,factor_env        & ! Proportion of rainfall going into DD from falling ppn
 ,ppn_dd_ref          ! Reference DD ppn mass

INTEGER ::          &
  errorstatus         ! Return code : 0 Normal Exit : >0 Error

CHARACTER(LEN=256)  &
  cmessage            ! Error message if return code >0

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (Len=*), PARAMETER ::    &
  RoutineName='DOWND_4a5a'

!-----------------------------------------------------------------------
! Start of main loop
!   Update precipitation and calculate mask for where precipitation
!   is liquid.
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('DOWND_4A5A',zhook_in,zhook_handle)

DO i=1,npnts
  b_dd_end(i) = .FALSE.
END DO

IF (k == kct+1) THEN
  DO i=1,npnts
    rain_dd(i) = 0.0
    rain_env(i) = 0.0
    snow_dd(i) = 0.0
    snow_env(i) = 0.0
  END DO
END IF

!----------------------------------------------------------------------
! Injection of precipitation from UD at level k
!----------------------------------------------------------------------

IF (dd_opt == 1) THEN
  !----------------------------------------------------------------------
  ! Corrected DD
  ! Transfer all the precipitation created in the UD into the environment
  !----------------------------------------------------------------------
  DO i=1,npnts

    IF (bwater_k(i)) THEN
      rain_env(i) = rain_env(i) + precip_k(i)
    ELSE
      snow_env(i) = snow_env(i) + precip_k(i)
    END IF
    precip_k(i) = 0.0
  END DO

ELSE
  !----------------------------------------------------------------------
  ! Original DD
  !----------------------------------------------------------------------
  DO i=1,npnts
    factor = 0.0
    IF (bdd_on(i) .AND. flx_ud_k(i) >  0.0) THEN
      factor = g * flx_dd_k(i)/flx_ud_k(i)
      factor = AMIN1(factor, 1.0)
    END IF

    IF (bwater_k(i)) THEN
      rain_dd(i)  = rain_dd(i) + precip_k(i)*factor
      rain_env(i) = rain_env(i) + precip_k(i)*(1.0-factor)
    ELSE
      snow_dd(i)  = snow_dd(i) + precip_k(i)*factor
      snow_env(i) = snow_env(i) + precip_k(i)*(1.0-factor)
    END IF

  END DO
END IF   ! dd_opt

!----------------------------------------------------------------------
! Interaction of downdraught with reserve of precipitation outside
! downdraught

! Based upon continuuity of precipitation mixing ratio within
! downdraught - either after injection of rain from UD in level
! k or with PPN mixing ratio in lowest preciptating layer

! IF downdraught increases in mass then water injected
! IF downdraught decreases in mass then water is removed

!----------------------------------------------------------------------


IF (dd_opt == 1) THEN
!----------------------------------------------------------------------
! Corrected DD
!----------------------------------------------------------------------
  DO i=1,npnts

    IF (bdd_on(i)) THEN

      IF (flx_ud_k(i) > 0.0) THEN
        factor_env = ddptef*flx_dd_k(i)/flx_ud_k(i)*  &
        delpk(i)/5000.0
      ELSE
        factor_env = 1.0*delpk(i)/5000.0
      END IF
      factor_env = AMIN1(factor_env,1.0)


      IF (factor_env >  0.0) THEN
        rain_dd(i)  = rain_dd(i) + rain_env(i)*factor_env
        rain_env(i) = rain_env(i) * (1.0-factor_env)
        snow_dd(i)  = snow_dd(i) + snow_env(i)*factor_env
        snow_env(i) = snow_env(i) * (1.0-factor_env)
      ELSE
        rain_env(i) = rain_env(i) - rain_dd(i)*factor_env
        rain_dd(i)  = rain_dd(i) * (1.0+factor_env)
        snow_env(i) = snow_env(i) - snow_dd(i)*factor_env
        snow_dd(i)  = snow_dd(i) * (1.0+factor_env)
      END IF

    END IF

  END DO   ! i

ELSE
!----------------------------------------------------------------------
! Original DD
!----------------------------------------------------------------------
  DO i=1,npnts

    IF (bdd_on(i)) THEN

      factor_env = 0.0
      IF (precip_k(i) >  0.0 .AND. flx_ud_k(i) >  0.0) THEN

      !---------------------------------------------------------------------
      ! Calculate new reference PPN mixing ratio
      ! DD PPN mixing ratio in layer km1 based on continuity
      ! with that in layer K
      !---------------------------------------------------------------------

        lr_ud_ref(i) = g * precip_k(i)/flx_ud_k(i)
        ppn_dd_ref   = rain_dd(i)+snow_dd(i)
      ELSE

      !---------------------------------------------------------------------
      ! DD PPN mixing ratio in layer km1 based on continuity
      ! with that in last preciptating UD layer
      !---------------------------------------------------------------------

        ppn_dd_ref = lr_ud_ref(i) * flx_dd_k(i)
      END IF

      !--------------------------------------------------------------------
      ! Inject ppn into DD from ppn falling outside of the DD
      !--------------------------------------------------------------------

      IF ((rain_env(i) + snow_env(i))  >   0.0) THEN
      !-------Already inside IF ( BDD_ON(I)) block----------------------------
        factor_env = ( (ppn_dd_ref * (1.0+ekm14(i))*                   &
                         (1.0+ekm34(i))*(1.0-amdetk(i))) -             &
                              (rain_dd(i)+snow_dd(i)) ) /              &
                              (rain_env(i)+snow_env(i))
        factor_env = AMIN1(factor_env,1.0)
        factor_env = AMAX1(factor_env,-1.0)
      END IF

      IF (factor_env >  0.0) THEN
        rain_dd(i)  = rain_dd(i) + rain_env(i)*factor_env
        rain_env(i) = rain_env(i) * (1.0-factor_env)
        snow_dd(i)  = snow_dd(i) + snow_env(i)*factor_env
        snow_env(i) = snow_env(i) * (1.0-factor_env)
      ELSE
        rain_env(i) = rain_env(i) - rain_dd(i)*factor_env
        rain_dd(i)  = rain_dd(i) * (1.0+factor_env)
        snow_env(i) = snow_env(i) - snow_dd(i)*factor_env
        snow_dd(i)  = snow_dd(i) * (1.0+factor_env)
      END IF

    END IF

    !--------------------------------------------------------------------
    ! Zero precipitation rate in layer k
    !--------------------------------------------------------------------

    precip_k(i) = 0.0

  END DO

END IF   !(DD_OPT)


!-----------------------------------------------------------------------
! Compress out on basis of bit vector BDDON - those points with a
! downdraught
!-----------------------------------------------------------------------

nddon=0

DO i=1,npnts
  IF (bdd_on(i)) THEN
     nddon = nddon+1
     index1(nddon) = i
  END IF
END DO

IF (nddon  /=  0) THEN
  DO i=1,nddon
    work(i,1) = thdd_k(index1(i))
    work(i,2) = qdd_k(index1(i))
    work(i,3) = the_k(index1(i))
    work(i,4) = the_km1(index1(i))
    work(i,5) = qe_k(index1(i))
    work(i,6) = qe_km1(index1(i))
    work(i,7) = dthbydt_k(index1(i))
    work(i,8) = dthbydt_km1(index1(i))
    work(i,9) = dqbydt_k(index1(i))
    work(i,10) = dqbydt_km1(index1(i))
    work(i,11) = flx_dd_k(index1(i))
    work(i,12) = p_km1(index1(i))
    work(i,13) = delpk(index1(i))
    work(i,14) = delpkm1(index1(i))
    work(i,15) = exk(index1(i))
    work(i,16) = exkm1(index1(i))
    work(i,17) = deltd(index1(i))
    work(i,18) = delqd(index1(i))
    work(i,19) = amdetk(index1(i))
    work(i,20) = ekm14(index1(i))
    work(i,21) = ekm34(index1(i))
    work(i,22) = rain_dd(index1(i))
    work(i,23) = snow_dd(index1(i))
    work(i,24) = cca(index1(i))
    qse_km1_c(i) = qse_km1(index1(i))

    bwork(i,1) = bdd_start(index1(i))
    bwork(i,2) = bddwt_k(index1(i))
    bwork(i,3) = bddwt_km1(index1(i))
    bwork(i,4) = bdd_on(index1(i))
    bwork(i,5) = b_dd_end(index1(i))
  END DO

  DO i=1,nddon
    ppn_mix_dd_c(i) = g*(rain_dd(index1(i)) +                               &
                         snow_dd(index1(i)))/flx_dd_k(index1(i))
  END DO

  IF (l_tracer) THEN

    DO ktra=1,ntra
      DO i=1,nddon
        tradd_k_c(i,ktra) = tradd_k(index1(i),ktra)
        trae_k_c(i,ktra) = trae_k(index1(i),ktra)
        trae_km1_c(i,ktra) = trae_km1(index1(i),ktra)
        dtra_k_c(i,ktra)  = dtrabydt_k(index1(i),ktra)
        dtra_km1_c(i,ktra) = dtrabydt_km1(index1(i),ktra)
        deltrad_c(i,ktra) = deltrad(index1(i),ktra)
      END DO
    END DO
  
  END IF


  !-----------------------------------------------------------------------
  ! Start downdraught calculation
  !-----------------------------------------------------------------------

  SELECT CASE ( i_convection_vn )
     CASE ( i_convection_vn_4a )

! DEPENDS ON: ddraught_4a
        CALL ddraught_4a (nddon,nddon_a,k,kct,ntra,                       &
                          l_tracer,                                       &
                          work(1,1),work(1,2),                            &
                          work(1,3),work(1,4),work(1,5),work(1,6),        &
                          qse_km1_c,                                      &
                          work(1,7),work(1,8),work(1,9),work(1,10),       &
                          work(1,11),work(1,12),work(1,13),work(1,14),    &
                          work(1,15),work(1,16),work(1,17),work(1,18),    &
                          work(1,19),work(1,20),work(1,21),work(1,22),    &
                          work(1,23),bwork(1,1),bwork(1,2),bwork(1,3),    &
                          bwork(1,4),bwork(1,5),                          &
                          deltrad_c, work(1,24),ppn_mix_dd_c,             &
                          tradd_k_c,trae_k_c,trae_km1_c,                  &
                          dtra_k_c,dtra_km1_c)

     CASE ( i_convection_vn_5a )

! DEPENDS ON: ddraught_5a
        CALL ddraught_5a (nddon,nddon_a,k,kct,ntra,                       &
                          l_tracer,                                       &
                          work(1,1),work(1,2),                            &
                          work(1,3),work(1,4),work(1,5),work(1,6),        &
                          qse_km1_c,                                      &
                          work(1,7),work(1,8),work(1,9),work(1,10),       &
                          work(1,11),work(1,12),work(1,13),work(1,14),    &
                          work(1,15),work(1,16),work(1,17),work(1,18),    &
                          work(1,19),work(1,20),work(1,21),work(1,22),    &
                          work(1,23),bwork(1,1),bwork(1,2),bwork(1,3),    &
                          bwork(1,4),bwork(1,5),                          &
                          deltrad_c, work(1,24),ppn_mix_dd_c,             &
                          tradd_k_c,trae_k_c,trae_km1_c,                  &
                          dtra_k_c,dtra_km1_c)

     CASE DEFAULT ! i_convection_vn

        errorstatus = 10
        WRITE (cmessage,'(A)') 'Convection scheme version value not recognised'
        WRITE (cmessage,'(A,I6)') '   i_convection_vn = ',i_convection_vn
        CALL Ereport ( RoutineName, errorstatus, cmessage)

  END SELECT ! i_convection_vn

  !-----------------------------------------------------------------------
  ! Expand requried vectors back to full fields
  !-----------------------------------------------------------------------

  DO i=1,nddon
    thdd_k(index1(i)) = work(i,1)
    qdd_k(index1(i)) = work(i,2)
    dthbydt_k(index1(i)) = work(i,7)
    dthbydt_km1(index1(i)) = work(i,8)
    dqbydt_k(index1(i)) = work(i,9)
    dqbydt_km1(index1(i)) = work(i,10)
    flx_dd_k(index1(i)) = work(i,11)
    rain_dd(index1(i)) = work(i,22)
    snow_dd(index1(i)) = work(i,23)

    bdd_start(index1(i)) = bwork(i,1)
    bddwt_k(index1(i)) = bwork(i,2)
    bddwt_km1(index1(i)) = bwork(i,3)
    bdd_on(index1(i)) = bwork(i,4)
    b_dd_end(index1(i)) = bwork(i,5)
  END DO

  IF (l_tracer) THEN

    DO ktra=1,ntra
      DO i=1,nddon
        tradd_k(index1(i),ktra) = tradd_k_c(i,ktra)
        dtrabydt_k(index1(i),ktra) = dtra_k_c(i,ktra)
        dtrabydt_km1(index1(i),ktra) = dtra_km1_c(i,ktra)
      END DO
   END DO

  END IF

END IF   ! nddon > 0

!-----------------------------------------------------------------------
! Reset precipitation falling through environment if downdraught
! did not form
!-----------------------------------------------------------------------

DO i=1,npnts
  IF (.NOT.bdd_on(i).AND..NOT.b_dd_end(i)) THEN
    rain_env(i) = rain_env(i)+rain_dd(i)
    snow_env(i) = snow_env(i)+snow_dd(i)
    rain_dd(i) = 0.0
    snow_dd(i) = 0.0
  END IF
END DO

!-----------------------------------------------------------------------
! Carry out change of phase calculation for precipitation falling
! through environment
!-----------------------------------------------------------------------

! DEPENDS ON: chg_phse
CALL chg_phse (npnts,k,rain_env,snow_env,dthbydt_km1,             &
                  exk,exkm1,delpkm1,the_k,the_km1,timestep,cca)

!-----------------------------------------------------------------------
! Evaporate rain falling through environment if layer k below
! cloud base
!-----------------------------------------------------------------------

SELECT CASE ( i_convection_vn )
   CASE ( i_convection_vn_4a )

! DEPENDS ON: pevp_bcb_4a
      CALL pevp_bcb_4a (npnts,k-1,iccb,the_km1,p_km1,qe_km1,qse_km1,delpkm1,  &
                        rain_env,snow_env,dthbydt_km1,dqbydt_km1,             &
                        exkm1,timestep,cca)

   CASE ( i_convection_vn_5a )

! DEPENDS ON: pevp_bcb
      CALL pevp_bcb (npnts,k-1,iccb,the_km1,p_km1,qe_km1,qse_km1,delpkm1,  &
                     rain_env,snow_env,dthbydt_km1,dqbydt_km1,             &
                     exkm1,timestep,cca)

   CASE DEFAULT ! i_convection_vn

      errorstatus = 10
      WRITE (cmessage,'(A)') 'Convection scheme version value not recognised'
      WRITE (cmessage,'(A,I6)') '   i_convection_vn = ',i_convection_vn
      CALL Ereport ( RoutineName, errorstatus, cmessage)

END SELECT ! i_convection_vn

!-----------------------------------------------------------------------
! Reset precipitation falling through environment if downdraught
! terminates
!-----------------------------------------------------------------------

DO i=1,npnts
  IF (b_dd_end(i)) THEN
    rain_env(i) = rain_env(i)+rain_dd(i)
    snow_env(i) = snow_env(i)+snow_dd(i)
    rain_dd(i) = 0.0
    snow_dd(i) = 0.0
  END IF
END DO

!-----------------------------------------------------------------------
! Update rain and snow
!-----------------------------------------------------------------------

IF (k == 2) THEN
  DO i=1,npnts
    rain(i) = rain(i)+rain_dd(i)+rain_env(i)
    snow(i) = snow(i)+snow_dd(i)+snow_env(i)
  END DO
END IF

IF (lhook) CALL dr_hook('DOWND_4A5A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE downd_4a5a
