! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the effect of convection upon the large-scale atmosphere
!
MODULE environ_6a_mod

IMPLICIT NONE

! Description:
!   Calculate the effect of convection upon the large-scale atmosphere
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.

CONTAINS

! Subroutine Interface:
SUBROUTINE environ_6a (k, npnts, np_full, ntra,                              &
                    delpk, delpkp1, delpkp12, delp_uv_k, delp_uv_kp1,     &
                    exk, exkp1,                                           &
                    thek, thekp1, qek, qekp1,                             &
                    qclek, qclekp1, qcfek, qcfekp1,                       &
                    thpk, thpkp1, qpk, qpkp1,                             &
                    qclpk, qclpkp1, qcfpk, qcfpkp1,                       &
                    thrk, qrk,                                            &
                    ekp14, amdetk, deltak, flxk,                          &
                    uek, uekp1, vek, vekp1,                               &
                    upk, upkp1, vpk, vpkp1,                               &
                    traek, traekp1, trapk, trapkp1,                       &
                    l_q_interact, l_mom_gk, l_tracer,                     &
                    blowst, bterm,                                        &
                    ! In/out
                    dthek, dqek, dqclek, dqcfek,                          &
                    eflux_u_ud, eflux_v_ud,                               &
                    duek, dvek,                                           &
                    dtraek,                                               &
                    ! Out
                    dthekp1, dqekp1, dqclekp1, dqcfekp1,                  &
                    duekp1, dvekp1,                                       &
                    dtraekp1)
                    

USE water_constants_mod, ONLY: lc, lf
USE cv_derived_constants_mod, ONLY: ls, lcrcp, lfrcp, lsrcp
USE atmos_constants_mod, ONLY: r, cp

USE cv_run_mod, ONLY: bl_cnv_mix

USE cloud_inputs_mod, ONLY: i_pc2_conv_coupling,                        &
                            starticeTKelvin, alliceTdegC

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER,INTENT(IN) :: k             ! present model layer
INTEGER,INTENT(IN) :: npnts         ! Number of points
INTEGER,INTENT(IN) :: np_full       ! Full vector length
INTEGER,INTENT(IN) :: ntra          ! Number of tracer variables

REAL,INTENT(IN) :: delpk(npnts)     ! pressure difference across layer k (Pa)
REAL,INTENT(IN) :: delpkp1(npnts)   ! pressure difference across layer k+1 (Pa)
REAL,INTENT(IN) :: delpkp12(npnts)  ! pressure diff. across layer k+1/2 (Pa)
REAL,INTENT(IN) :: delp_uv_k(npnts) ! pressure difference across UV 
                                    ! layer k (Pa)
REAL,INTENT(IN) :: delp_uv_kp1(npnts) ! pressure difference across UV
                                    ! layer k+1 (Pa)
REAL,INTENT(IN) :: exk(npnts)       ! Exner ratio at mid-point of layer k
REAL,INTENT(IN) :: exkp1(npnts)     ! Exner ratio at mid-point of layer k+1

REAL,INTENT(IN) :: thek(npnts)      ! Env. pot. temperature in layer k (K)
REAL,INTENT(IN) :: thekp1(npnts)    ! Env. pot. temperature in layer k+1 (K)
REAL,INTENT(IN) :: qek(npnts)       ! Env. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qekp1(npnts)     ! Env. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qclek(npnts)     ! Env. qcl in layer k (kg/kg)
REAL,INTENT(IN) :: qclekp1(npnts)   ! Env. qcl in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qcfek(npnts)     ! Env. qcf in layer k (kg/kg)
REAL,INTENT(IN) :: qcfekp1(npnts)   ! Env. qcf in layer k+1 (kg/kg)
REAL,INTENT(IN) :: thpk(npnts)      ! Par. pot. temperature in layer k (K)
REAL,INTENT(IN) :: thpkp1(npnts)    ! Par. pot. temperature in layer k+1 (K)
REAL,INTENT(IN) :: qpk(npnts)       ! Par. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qpkp1(npnts)     ! Par. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qclpk(npnts)     ! Par. qcl in layer k (kg/kg)
REAL,INTENT(IN) :: qclpkp1(npnts)   ! Par. qcl in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qcfpk(npnts)     ! Par. qcf in layer k (kg/kg)
REAL,INTENT(IN) :: qcfpkp1(npnts)   ! Par. qcf in layer k+1 (kg/kg)
REAL,INTENT(IN) :: thrk(npnts)      ! pot. temperature of forced detrained
                                    ! parcel in layer k (K)
REAL,INTENT(IN) :: qrk(npnts)       ! Specific humidity of forced detrained
                                    ! parcel in layer k (kg/kg)
REAL,INTENT(IN) :: ekp14(npnts)     ! Entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: amdetk(npnts)    ! Mixing detrainment coefficient at level k 
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: deltak(npnts)    ! Parcel forced detrainment rate in 
                                    ! layer k multiplied by layer thickness
REAL,INTENT(IN) :: flxk(npnts)      ! Parcel massflux in layer k (Pa/s) 
REAL,INTENT(IN) :: uek(npnts)       ! Env. U in layer k (m/s)
REAL,INTENT(IN) :: uekp1(npnts)     ! Env. U in layer k+1 (m/s)
REAL,INTENT(IN) :: vek(npnts)       ! Env. V in layer k (m/s)
REAL,INTENT(IN) :: vekp1(npnts)     ! Env. V in layer k+1 (m/s)
REAL,INTENT(IN) :: upk(npnts)       ! Par. U in layer k (m/s)
REAL,INTENT(IN) :: upkp1(npnts)     ! Par. U in layer k+1 (m/s)
REAL,INTENT(IN) :: vpk(npnts)       ! Par. V in layer k (m/s)
REAL,INTENT(IN) :: vpkp1(npnts)     ! Par. V in layer k+1 (m/s)
REAL,INTENT(IN) :: traek(np_full,ntra)    ! Env. tracer content 
                                          ! in layer k (kg/kg)
REAL,INTENT(IN) :: traekp1(np_full,ntra)  ! Env. tracer content 
                                          ! in layer k+1 (kg/kg)
REAL,INTENT(IN) :: trapk(np_full,ntra)    ! Par. tracer content 
                                          ! in layer k (kg/kg)
REAL,INTENT(IN) :: trapkp1(np_full,ntra)  ! Par. tracer content 
                                          ! in layer k+1 (kg/kg)

LOGICAL,INTENT(IN) :: l_q_interact  ! True if PC2 is switched on
LOGICAL,INTENT(IN) :: l_mom_gk      ! Switch for inclusion of Gregory-Kershaw
                                    ! CMT
LOGICAL,INTENT(IN) :: l_tracer      ! Switch for tracers

LOGICAL,INTENT(IN) :: blowst(npnts) ! mask for those points at which stability
                                    ! is low enough for convection to occur
LOGICAL,INTENT(IN) :: bterm(npnts)  ! Mask for parcels which terminate 
                                    ! in layer k+1
!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------

! Convection increments to model fields at level k
! IN before processes at level k. This may be non-zero because of smoothed
! forced detrainment and the initial perturbation
! OUT after processes at level k
REAL,INTENT(INOUT) :: dthek(npnts)  ! Increment to p. temperature
                                    ! in layer k (K/s)
REAL,INTENT(INOUT) :: dqek(npnts)   ! Increment to spec. humidity 
                                    ! in layer k (kg/kg/s)
REAL,INTENT(INOUT) :: dqclek(npnts) ! Increment to qcl in layer k (kg/kg/s)
REAL,INTENT(INOUT) :: dqcfek(npnts) ! Increment to qcf in layer k (kg/kg/s)

REAL,INTENT(INOUT) :: eflux_u_ud(npnts) ! Eddy flux of momentum to UD
                                    ! IN at bottom of layer
                                    ! OUT at top of layer
REAL,INTENT(INOUT) :: eflux_v_ud(npnts) ! Eddy flux of momentum to UD
                                    ! IN at bottom of layer
                                    ! OUT at top of layer
! Convection increments to winds and tracers at level k
! IN before processes at level k. This may be non-zero because of smoothed
! forced detrainment and the initial perturbation
! OUT after processes at level k
REAL,INTENT(INOUT) :: duek(npnts)   ! Increment to U in layer k (m/s**2)
REAL,INTENT(INOUT) :: dvek(npnts)   ! Increment to V in layer k (m/s**2)
REAL,INTENT(INOUT) :: dtraek(np_full,ntra)  ! Increment to tracer in layer k
                                            ! (kg/kg/s) 

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
! Convection increments to model fields at level k+1
REAL,INTENT(OUT) :: dthekp1(npnts)  ! Increment to p. temperature
                                      ! in layer k+1 (K/s)
REAL,INTENT(OUT) :: dqekp1(npnts)   ! Increment to spec. humidity 
                                      ! in layer k+1 (kg/kg/s)
REAL,INTENT(OUT) :: dqclekp1(npnts) ! Increment to qcl in layer k+1 (kg/kg/s)
REAL,INTENT(OUT) :: dqcfekp1(npnts) ! Increment to qcf in layer k+1 (kg/kg/s)
REAL,INTENT(OUT) :: duekp1(npnts)   ! Increment to U in layer k+1 (m/s**2)
REAL,INTENT(OUT) :: dvekp1(npnts)   ! Increment to V in layer k+1 (m/s**2)
REAL,INTENT(OUT) :: dtraekp1(np_full,ntra) ! Increment to tracer in layer k+1
                                           ! (kg/kg/s)

!-----------------------------------------------------------------------
! Variables that are defined locally
!-----------------------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER :: i              ! loop counter 
INTEGER :: ktra           ! Loop counter for tracers

! Smoothed forced detrainment variables
REAL :: tmp_fd_dthek      ! F. detrainment P.temp inc across levels k and k+1
REAL :: tmp_fd_dqek       ! F. detrainment humidity inc across levels k and k+1
REAL :: tmp_fd_dqclek     ! F. detrainment qcl inc across levels k and k+1
REAL :: tmp_fd_dqcfek     ! F. detrainment qcf across levels k and k+1
REAL :: tmp_fd_dtraek     ! F. detrainment tracer across levels k and k+1

REAL :: el                ! Latent heat of condensation or (condensation plus
                          ! fusion) (J/kg)
REAL :: flxbydpk(npnts)   ! mass flux divided by layer thickness (1/s)
REAL :: flx_u_kp0p5       ! Flux of zonal momentum in cloud at top of current
                          ! uv layer
REAL :: flx_v_kp0p5       ! Flux of meridional momentum in cloud at top of
                          ! current uv layer
REAL :: tmp_dqclek        ! Storage space for liquid condensate rate.
REAL :: tmp_dqcfek        ! Storage space for frozen condensate rate.
REAL :: dth_detcond(npnts)! Used to calculate the contribution to the theta
                          ! increment from the evaporation of detrained 
                          ! 0.0 if PC2 is on.
REAL :: dq_detcond(npnts) ! Used to calculate the contribution to the humidity
                          ! increment from the evaporation of detrained 
                          ! 0.0 if PC2 is on.
REAL :: frac_icek         ! Fraction of ice water at level k
REAL :: frac_liqk         ! Fraction of liquid water at level k

! Parameters

REAL,PARAMETER :: a_smth    = 0.5     ! Parameter determining the weighting
                                      ! between the forced detrainment 
                                      ! increments at k and k-1 

! Several options are available:
INTEGER,PARAMETER :: pc2_conv_original = 1        
! As described in Wilson et al (2008). Condensate increments from 
! detrainment and subsidence advection (combined) are turned into cloud 
! fraction increments using the inhomogeneous forcing method. Note that
! the abrupt change in phase of condensate in the convective plume is
! replicated in the condensate increments due to convection.
INTEGER,PARAMETER :: pc2_conv_maxincldqc = 2
! Unpublished. As above + Protects (in ni_conv_ctl) against 
! generation of inconsistently low cloud fraction implying 
! very high in-cloud condensate amounts.
INTEGER,PARAMETER :: pc2_conv_smooth_liqice = 3
! Unpublished. As above + The phase of the detrained condensate varies 
! smoothly according to the ambient temperature rather than having 
! an abrupt change.

!---------------------------------------------------------------------

IF (lhook) CALL dr_hook('ENVIRON_6A',zhook_in,zhook_handle)

IF (l_q_interact) THEN
! PC2 is on and therefore the detrained condensate is not evaporated
! and does not affect the theta and q increments
  DO i=1,npnts
    dth_detcond(i)  = 0.0
    dq_detcond(i)   = 0.0
  END DO
ELSE
! PC2 is off and therefore the detrained condensate is evaporated
! and will affect the theta and q increments
  DO i=1,npnts
    dth_detcond(i)  = (lcrcp*qclpk(i) + lsrcp*qcfpk(i))/exk(i)
    dq_detcond(i)   =        qclpk(i) +       qcfpk(i)
  END DO  
END IF

DO i=1,npnts

!----------------------------------------------------------------------
! Calculate parcel mass flux divided by the thickness of layer k.
! This value is used in several places in the subroutine.
!----------------------------------------------------------------------
  flxbydpk(i)   = flxk(i)/delpk(i)

!-----------------------------------------------------------------------
!         Potential temperature increment
!-----------------------------------------------------------------------
!                Smoothed forced detrainment term
  tmp_fd_dthek  = deltak(i) * (1.0-amdetk(i))                          &
                * ( thrk(i) - thek(i) - dth_detcond(i) )

!           Theta increment at level k
  dthek(i)      = dthek(i) + flxbydpk(i)                               &
!                Compensating subsidence
                * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))         &
                * ( thekp1(i)-thek(i) )                                &
!                Smoothed forced detrainment  
                + a_smth*tmp_fd_dthek                                  &
!                Mixing detrainment
                + amdetk(i)                                            &
                * ( thpk(i) - thek(i) - dth_detcond(i) ) )
                
!           Theta increment at level k+1 due 
!           to smoothed  forced detrainment
  dthekp1(i)    = flxbydpk(i) * (1.0-a_smth)                           &
                * exk(i)/exkp1(i)                                      &
                * delpk(i)/delpkp1(i) * tmp_fd_dthek

!-----------------------------------------------------------------------
!         Humidity increment
!-----------------------------------------------------------------------
!                Smoothed forced detrainment term
  tmp_fd_dqek   = deltak(i) * (1.0-amdetk(i))                          &
                * ( qrk(i) - qek(i) + dq_detcond(i) )

!           q increment at level k
  dqek(i)       = dqek(i) + flxbydpk(i)                                &
!                Compensating subsidence
                * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))         &
                * ( qekp1(i)-qek(i) )                                  &
!                Smoothed forced detrainment  
                + a_smth*tmp_fd_dqek                                   &
!                Mixing detrainment
                + amdetk(i)                                            &
                * ( qpk(i) - qek(i) + dq_detcond(i) ) )

!           q increment at level k+1 due 
!           to smoothed  forced detrainment
  dqekp1(i)     = flxbydpk(i) * (1.0-a_smth)                           &
                * delpk(i)/delpkp1(i)*tmp_fd_dqek
END DO  !npnts

IF (l_q_interact) THEN
  DO i=1,npnts
!-----------------------------------------------------------------------
!         PC2 is on. Set the qcl and qcf increments appropriately
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!         Liquid condensate increment
!         A temporary variable is used because PC2 modifies the
!         increment
!-----------------------------------------------------------------------
!                Smoothed forced detrainment term
    tmp_fd_dqclek = deltak(i) * (1.0-amdetk(i))                        &
                  * ( qclpk(i)-qclek(i) )

!           qcl increment at level k
!           previous contribution at level k not included here
!           but is added later because of PC2 messing about with
!           this increment.
    tmp_dqclek    = flxbydpk(i)                                        &
!                  Compensating subsidence
                  * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))       &
                  * ( qclekp1(i)-qclek(i) )                            &
!                  Smoothed forced detrainment  
                  + a_smth*tmp_fd_dqclek                               &
!                  Mixing detrainment
                  + amdetk(i)                                          &
                  * ( qclpk(i) - qclek(i) ) )
                
!           qcl increment at level k+1 due 
!           to smoothed  forced detrainment
    dqclekp1(i)   = flxbydpk(i) * (1.0-a_smth)                         &
                  * delpk(i)/delpkp1(i)*tmp_fd_dqclek


!-----------------------------------------------------------------------
!         Frozen condensate increment
!         A temporary variable is used because PC2 modifies the
!         increment
!-----------------------------------------------------------------------
!                Smoothed forced detrainment term
    tmp_fd_dqcfek = deltak(i) * (1.0-amdetk(i))                        &
                  * ( qcfpk(i)-qcfek(i) )

!           qcf increment at level k
!           previous contribution at level k not included here
!           but is added later because of PC2 messing about with
!           this increment.
    tmp_dqcfek    = flxbydpk(i)                                        &
!                  Compensating subsidence
                  * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))       &
                  * ( qcfekp1(i)-qcfek(i) )                            &
!                  Smoothed forced detrainment  
                  + a_smth*tmp_fd_dqcfek                               &
!                  Mixing detrainment
                  + amdetk(i)                                          &
                  * ( qcfpk(i) - qcfek(i) ) )


!           qcf increment at level k+1 due 
!           to smoothed  forced detrainment
    dqcfekp1(i)   = flxbydpk(i) * (1.0-a_smth)                         &
                  * delpk(i)/delpkp1(i)*tmp_fd_dqcfek

! ----------------------------------------------------------------------
!       Adjust temperature increment and condensate increments to
!       take account of a smoother transistion between water and ice
! ----------------------------------------------------------------------

    IF (i_pc2_conv_coupling == pc2_conv_original    .OR.               &
        i_pc2_conv_coupling == pc2_conv_maxincldqc) THEN
      ! Original method
      dqclek(i) = dqclek(i) + tmp_dqclek
      dqcfek(i) = dqcfek(i) + tmp_dqcfek
    ELSE
      ! Convective plume is either liq or ice with an abrupt change. 
      ! Partition the condensate increments detrained from the 
      ! convective plume onto the large-scale more smoothly over
      ! a range of temperature. Adjust for latent heating.
      
      frac_icek   = (thek(i)*exk(i)-starticeTKelvin)/alliceTdegC 
      frac_icek   = MAX(0.0,frac_icek) 
      frac_icek   = MIN(frac_icek,1.0) 
      frac_liqk   = 1.0-frac_icek

      IF (tmp_dqclek > 0.0) THEN 
        dqcfek(i) = dqcfek(i) + frac_icek*tmp_dqclek 
        ! Extra freezing heats 
        dthek(i)  = dthek(i)  + frac_icek*tmp_dqclek*lfrcp/exk(i)
        dqclek(i) = dqclek(i) + frac_liqk*tmp_dqclek
      ELSE
        dqclek(i) = dqclek(i) + tmp_dqclek
      END IF

      IF (tmp_dqcfek > 0.0) THEN 
        dqclek(i) = dqclek(i) + frac_liqk*tmp_dqcfek 
        ! Melting cools 
        dthek(i)  = dthek(i)  - frac_liqk*tmp_dqcfek*lfrcp/exk(i)
        dqcfek(i) = dqcfek(i) + frac_icek*tmp_dqcfek
      ELSE
        dqcfek(i) = dqcfek(i) + tmp_dqcfek
      END IF

    END IF

  END DO  !npnts
ELSE
!-----------------------------------------------------------------------
!         PC2 is off. Set the qcl and qcf increments to zero.
!-----------------------------------------------------------------------
  DO i=1,npnts
    dqclek(i)     = 0.0
    dqclekp1(i)   = 0.0

    dqcfek(i)     = 0.0
    dqcfekp1(i)   = 0.0
  END DO
END IF  !l_q_interact

! ---------------------------------------------------------------------
!  Calculate effect of convection upon momentum of layer k and 
!  do terminal detrainment of momentum.
!
!  Rate of change of wind field by convection is estimated using a
!  divergence of vertical eddy momentum flux across the layer.
! --------------------------------------------------------------------
!
! All convective momentum transport calculations for the cumulus convection
! (deep and shallow) done when convection terminates.


IF(l_mom_gk) THEN      ! Gregory-Kershaw CMT

  DO i=1,npnts
!----------------------------------------------------------------------
! Estimate eddy flux at top of current UV layer due to convection
!----------------------------------------------------------------------

    flx_u_kp0p5 = flxk(i) * (1.0-amdetk(i)) * (1.0-deltak(i)) *      &
                            (1.0+ekp14(i)) * (upk(i)-uekp1(i))
    flx_v_kp0p5 = flxk(i) * (1.0-amdetk(i)) * (1.0-deltak(i)) *      &
                            (1.0+ekp14(i)) * (vpk(i)-vekp1(i))

    IF (blowst(i)) THEN
!----------------------------------------------------------------------
! Initial convecting layer - no flux at base of layer
!----------------------------------------------------------------------
      duek(i) = duek(i) - flx_u_kp0p5 / delp_uv_k(i)
      dvek(i) = dvek(i) - flx_v_kp0p5 / delp_uv_k(i)
!----------------------------------------------------------------------
! Store eddy flux at top of current UV layer ready for calculation
! of next layer.
!----------------------------------------------------------------------
      eflux_u_ud(i) = flx_u_kp0p5
      eflux_v_ud(i) = flx_v_kp0p5

    ELSE
!----------------------------------------------------------------------
! Convecting layer - take eddy flux divergence across the layer
!----------------------------------------------------------------------
      duek(i) = duek(i) - ( (flx_u_kp0p5 - eflux_u_ud(i)) /delp_uv_k(i) )  

      dvek(i) = dvek(i) - ( (flx_v_kp0p5 - eflux_v_ud(i)) /delp_uv_k(i) )

!----------------------------------------------------------------------
! Store eddy flux at top of curent UV layer ready for calculation of 
! next layer
!----------------------------------------------------------------------
      eflux_u_ud(i) = flx_u_kp0p5
      eflux_v_ud(i) = flx_v_kp0p5

    END IF

    IF(bterm(i))THEN
!----------------------------------------------------------------------
! Convection terminates - calculate increment due to convection in top 
! layer - no flux out of top layer.
!----------------------------------------------------------------------
      duekp1(i)  = eflux_u_ud(i) / delp_uv_kp1(i)
      dvekp1(i)  = eflux_v_ud(i) / delp_uv_kp1(i)
!----------------------------------------------------------------------
! Zero eddy flux out of top layer.
!----------------------------------------------------------------------
      eflux_u_ud(i) = 0.0
      eflux_v_ud(i) = 0.0
    ELSE
      duekp1(i)  = 0.0
      dvekp1(i)  = 0.0
    END IF

  END DO

END IF   ! l_mom_gk

!----------------------------------------------------------------------
!  Effect of convection on tracer content of layer k.
!  (looping over number of tracer variables)
!  and do terminal detrainment of tracer.
!----------------------------------------------------------------------

IF(l_tracer)THEN
  DO ktra = 1,ntra
    DO i = 1,npnts

!                   Smoothed forced detrainment term
      tmp_fd_dtraek   =  deltak(i) * (1.0-amdetk(i))                   &
                      * ( trapk(i,ktra)-traek(i,ktra) )

!                   tracer increment at level k
      dtraek(i,ktra)  = dtraek(i,ktra) + flxbydpk(i)                   &
!                   Compensating subsidence
                      * ( (1-amdetk(i))*(1.0-deltak(i))*(1+ekp14(i))   &
                      * (traekp1(i,ktra)-traek(i,ktra))                &
!                   Smoothed forced detrainment  
                      + a_smth*tmp_fd_dtraek                           &
!                   Mixing detrainment
                      + amdetk(i)                                      &
                      * ( trapk(i,ktra)-traek(i,ktra) ) )

!                   tracer increment at level k+1 due 
!                   to smoothed  forced detrainment
      dtraekp1(i,ktra)= flxbydpk(i) *(1.0-a_smth)                      &
                      * delpk(i)/delpkp1(i)*tmp_fd_dtraek

    END DO  ! i
  END DO    ! ktra
END IF      ! l_tracer

IF (lhook) CALL dr_hook('ENVIRON_6A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE environ_6a
END MODULE environ_6a_mod
