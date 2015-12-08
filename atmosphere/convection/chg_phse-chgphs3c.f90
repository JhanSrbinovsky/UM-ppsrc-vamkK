! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Change phase for points where no downdraught occurring.
!
SUBROUTINE chg_phse (npnts,k,rain,snow,dthbydt_km1,               &
                     exk,exkm1,delpkm1,the_k,the_km1,timestep,    &
                     cca)

USE earth_constants_mod, ONLY: g

USE cv_run_mod, ONLY:                                             &
    l_phase_lim
USE cv_param_mod, ONLY:                                           &
    cldarea
USE water_constants_mod, ONLY: lf, tm
USE atmos_constants_mod, ONLY: cp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Change phase for points where no downdraught occurring.
! Updates poential temperature of layer k as precipitation changes phase
! in situ.
! Add latent heating where precipittaion crosses a melting of frezing level.
!
! See UM Documentation paper 27.
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
 ,k                      ! Model layer

REAL, INTENT(IN) ::    &
  timestep               ! Timestep

REAL, INTENT(IN) :: &
  exk(npnts)        & ! Exner ratio for layer k
 ,exkm1(npnts)      & ! Exner ratio for layer k-1
 ,delpkm1(npnts)    & ! Pressure difference across layer k-1 (Pa)
 ,the_k(npnts)      & ! potential temperature of environment in layer k (K)
 ,the_km1(npnts)    & ! potential temperature of environment in layer k-1 (K)
 ,cca(npnts)          ! Convective cloud amount


!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------
REAL, INTENT(INOUT) :: &
  rain(npnts)          & ! IN Amount of falling rain (kg/m**2/s)
                         ! OUT Updated amount of falling rain (kg/m**2/s)
 ,snow(npnts)          & ! IN Amount of falling Snowfall (kg/m**2/s)
                         ! OUT Updated amount of falling snow (kg/m**2/s)
 ,dthbydt_km1(npnts)     ! IN Increment to model potential temperature in
                         !    layer k-1
                         ! OUT Updated increment to model potential 
                         !     temperature in layer k-1 due to change of phase

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::        & 
  i                 ! loop counter

LOGICAL ::       &
  bppnwt_k       & ! Mask where precipitation is liquid in layer k
 ,bppnwt_km1       ! Mask where precipitation is liquid in layer k-1


REAL ::           &
  factor          & ! Used in the calculation of change of phase of falling
                    ! precipitation.
 ,wpc             & ! Amount of precipitation which can change phase
 ,ca              & ! Local cloud area
 ,the_km1_new       ! THE_KM1 updated with all increments prior to this 
                    ! subroutine

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CHG_PHSE',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!   Add latent heating where precip crosses a melting or freezing level
!
!   UM Documentation paper 27
!   Section (11), equation (42)
! ----------------------------------------------------------------------
!

IF (l_phase_lim) THEN

  DO i=1,npnts
    the_km1_new = the_km1(i) + timestep*dthbydt_km1(i)
    bppnwt_k = the_k(i) >  tm/exk(i)
    bppnwt_km1 = the_km1_new  >   tm/exkm1(i)
    factor = lf*g/(exkm1(i)*cp*delpkm1(i))
    ca = cldarea * cca(i)

    ! Freeze

    IF (.NOT.bppnwt_km1.AND.(bppnwt_k.OR.rain(i) >  0.0)) THEN
      wpc = MIN( rain(i),                                         &
               ca * (tm/exkm1(i) - the_km1_new) /                 &
                    (timestep * factor) )
      dthbydt_km1(i) = dthbydt_km1(i) + wpc * factor
      snow(i) = snow(i) + wpc
      rain(i) = rain(i) - wpc
    END IF

    ! Melt

    IF (bppnwt_km1.AND.(.NOT.bppnwt_k.OR.snow(i) >  0.0)) THEN
      wpc = MIN( snow(i),                                         &
               ca * (the_km1_new - tm/exkm1(i)) /                 &
                    (timestep * factor) )
      dthbydt_km1(i) = dthbydt_km1(i) - wpc * factor
      rain(i) = rain(i) + wpc
      snow(i) = snow(i) - wpc
    END IF
  END DO

ELSE

  DO i=1,npnts
    bppnwt_k = the_k(i)*exk(i) >  tm
    bppnwt_km1 = the_km1(i)*exkm1(i) >  tm
    factor = lf*g/(exkm1(i)*cp*delpkm1(i))
    ! Freeze
    IF (.NOT.bppnwt_km1.AND.(bppnwt_k.OR.rain(i) >  0.0)) THEN
      dthbydt_km1(i) = dthbydt_km1(i)+rain(i)*factor
      snow(i) = snow(i)+rain(i)
      rain(i) = 0.0
    END IF
    ! Melt
    IF (bppnwt_km1.AND.(.NOT.bppnwt_k.OR.snow(i) >  0.0)) THEN
      dthbydt_km1(i) = dthbydt_km1(i)-snow(i)*factor
      rain(i) = rain(i)+snow(i)
      snow(i) = 0.0
    END IF
  END DO

END IF     ! IF L_PHASE_LIM

IF (lhook) CALL dr_hook('CHG_PHSE',zhook_out,zhook_handle)

RETURN
END SUBROUTINE chg_phse

