! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates layer dependent constants for layer k, downdraught code
!
SUBROUTINE layer_dd_4a5a(npnts,k,kct,the_k,the_km1,flx_strt,           &
                    p_layer_centres,                              &
                    p_layer_boundaries,                           &
                    exner_layer_centres,                          &
                    exner_km12,                                   &
                    exner_kp12,exner_km32,pstar,pk,pkm1,delpk,    &
                    delpkm1,exk,exkm1,amdetk,ekm14,ekm34,kmin,    &
                    bddi, recip_pstar )

USE cv_run_mod, ONLY:                                             &
    dd_opt

USE cv_param_mod, ONLY:                                             &
    ae2, ddcoef1, ddcoef2, det_lyr

USE water_constants_mod, ONLY: tm

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Calculates layer dependent constants for layer k, downdraught code
!
!  See UM Documentation paper No 27
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,k                    & ! present model layer
 ,kct                    ! convective cloud top layer       


LOGICAL, INTENT(IN) ::  &
  bddi(npnts)             ! Mask for points where downdraught may initiate

REAL, INTENT(IN) ::                      &
  p_layer_centres(npnts,0:kct+2)         & ! Pressure (Pa)
 ,p_layer_boundaries(npnts,0:kct+1)      & ! Pressure at half level above
                                           ! p_layer_centres (Pa)
 ,exner_layer_centres(npnts,0:kct+1)       ! exner pressure

REAL, INTENT(IN) ::  &
  pstar(npnts)       & ! Surface pressure (Pa)
 ,exner_km12(npnts)  & ! Exner function at layer k-1/2
 ,exner_kp12(npnts)  & ! Exner function at layer k+1/2
 ,exner_km32(npnts)  & ! Exner function at layer k-3/2
 ,flx_strt(npnts)    & ! Updraught mas flux at level where downdraught starts
                       !  (Pa/s)
 ,the_k(npnts)       & ! Pontential temperature of environment 
                       ! in layer k (K)
 ,the_km1(npnts)     & ! Pontential temperature of environment 
                       ! in layer k-1 (K)
 ,recip_pstar(npnts)   ! 1/pstar (Pa)


!---------------------------------------------------------------------
! Variables which are output
!---------------------------------------------------------------------

INTEGER , INTENT(OUT) :: &
  kmin(npnts)              ! freezing level

REAL, INTENT(OUT) ::   &  
  pk(npnts)            & ! Pressure of layer k (Pa)
 ,pkm1(npnts)          & ! Pressure of layer k-1 (Pa)
 ,delpk(npnts)         & ! Pressure difference across layer k   (Pa)
 ,delpkm1(npnts)       & ! Pressure difference across layer k-1 (Pa)
 ,amdetk(npnts)        & ! Mixing detrainment at level k multiplied by 
                         ! appropriate layer thickness
 ,ekm14(npnts)         & ! exner ratio at layer k-1/4
 ,ekm34(npnts)         & ! exner ratio at layer k-3/4
 ,exk(npnts)           & ! Exner ratio at for layer k
 ,exkm1(npnts)           ! Exner ratio at for layer k-1

!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER ::        & 
  i                 ! loop counter 

REAL  ::     &
  ttk        & ! Temperature store at layer k
 ,ttkm1      & ! Temperature store at layer k-1
 ,thkm12     & ! Potential temperature store at layer k-1/2
 ,ttkm12     & ! Temperature store at layer k-1/2
 ,incr_fac   & ! Increment factor for entrainment rates at freezing level
 ,pu         &
 ,pl         &
 ,ddcoef2a     ! coefficient used in calculation of downdraught
               ! entrainment rates

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('LAYER_DD_4A5A',zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Set kmin to initial value
!  Calculate PK, DELPK and EXNER function - If k = KCT then
!  values for previous pass through routine at (k-1)+1 are taken.
!----------------------------------------------------------------------

IF (k == kct+1) THEN
  DO i=1,npnts
    kmin(i) = kct+2
    pk(i) = p_layer_centres(i,k)
    delpk(i) =  -(p_layer_boundaries(i,k) - p_layer_boundaries(i,k-1))
    exk(i) = exner_layer_centres(i,k)
  END DO
ELSE
  DO i=1,npnts
    pk(i) = pkm1(i)
    delpk(i) = delpkm1(i)
    exk(i) = exkm1(i)
  END DO
END IF

! ---------------------------------------------------------------------
!  Calculate PKM1, DELPKM1
!  Calculate EXNER functions at mid-layer k and k-1, and
!  difference of exner function across layer k
! ---------------------------------------------------------------------

DO i=1,npnts
  pkm1(i) = p_layer_centres(i,k-1)
  delpkm1(i) =  -(p_layer_boundaries(i,k-1)                       &
                - p_layer_boundaries(i,k-2))
  exkm1(i) = exner_layer_centres(i,k-1)
END DO

!---------------------------------------------------------------------
! Set DDCOEF2A depending upon which revision of the DD code is used.
!---------------------------------------------------------------------

IF (dd_opt == 1) THEN
  ddcoef2a=2.0
ELSE
  ddcoef2a=ddcoef2
END IF
!
! ---------------------------------------------------------------------
!  Calculate freezing level : Check if freezing level in this layer
! ---------------------------------------------------------------------
!
DO i=1,npnts
  IF (kmin(i) == kct+2) THEN
    ttk = the_k(i)*exk(i)
    ttkm1 = the_km1(i)*exkm1(i)
    thkm12 = (the_km1(i)+the_k(i))*0.5
    ttkm12 = thkm12*exner_km12(i)
    IF (ttkm12  >=  tm .AND. ttk  <   tm) THEN
      kmin(i) = k
    ELSE IF (ttkm1  >=  tm .AND. ttkm12  <   tm) THEN
      kmin(i) = k-1
    END IF
  END IF


! ---------------------------------------------------------------------
!  Calculate entrainment coefficients multiplied by
!  appropriate layer thickness
!
!  Calculate mixing detrainment coefficient multiplied by
!  appropriate layer thickness
!
!  UM DOCUMENTATION PAPER 27
!  Section (2C), Equation(14)
! ---------------------------------------------------------------------

  IF (pk(i) <  pstar(i)-det_lyr) THEN
    ekm14(i) = ae2 * (p_layer_boundaries(i,k-1)-pk(i)) * recip_pstar(i)

    ekm34(i) = ae2 * (pkm1(i)-p_layer_boundaries(i,k-1)) * recip_pstar(i)

    amdetk(i) = (ekm14(i)+ekm34(i)) * (1.0-1.0/ae2)
  ELSE
    ekm14(i) = 0.0
    ekm34(i) = 0.0
    amdetk(i) = delpk(i) /(pstar(i)- p_layer_boundaries(i,k))
  END IF

  IF (bddi(i)) THEN

    IF (k == kmin(i) .AND. pk(i) <  pstar(i)-det_lyr) THEN
      incr_fac = flx_strt(i)*ddcoef1*recip_pstar(i)
      IF (incr_fac >  6.0) incr_fac=6.0
      ekm14(i) = ekm14(i)*incr_fac
      ekm34(i) = ekm34(i)*incr_fac
    ELSE
      ekm14(i) = ekm14(i)*ddcoef2a
      ekm34(i) = ekm34(i)*ddcoef2a
      IF ((dd_opt == 1) .OR. (kmin(i) /= kct+2 .AND. k <  kmin(i)            & 
               .AND. pk(i) < pstar(i)-det_lyr)) THEN
        amdetk(i) = amdetk(i)*ddcoef2a
      END IF
    END IF

  END IF     ! bddi
END DO       ! i

IF (lhook) CALL dr_hook('LAYER_DD_4A5A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE layer_dd_4a5a
