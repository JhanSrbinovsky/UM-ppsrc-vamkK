! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! forced detrainment calculation
!
SUBROUTINE detrain_4a5a (npnts,  bwkp1,bgmk,                   &
                         thek,qek,thpk,qpk,qsek,dqsk,          &
                         thekp1,qekp1,qsekp1,dqskp1,           &
                         ekp14,ekp34,pk,pkp1,exk,exkp1,xsbmin, &
                         ! In/out
                         bgmkp1, thpkp1,qpkp1,xsqkp1,          &
                         ! Out
                         deltak,thrk,qrk )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!    forced detrainment calculation
!
!  Subroutine thp_det calculates the potential temperature of the parcel
!  in layer k+1 after forced detrainment.
!
!  Subroutine thetar calculates the potential temperature of the air in
!   layer k undergoing forced detrainment.
!
!  Subroutine det_rate calculates the forced detrainment rate of the
!   ensemble in layer k.
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
  npnts                  ! Number of points

LOGICAL, INTENT(IN) :: &
  bwkp1(npnts)           ! Mask for whether condensate is liquid in layer k+1

REAL,INTENT(IN) :: &
  thek(npnts)      & ! potential temperature of cloud environment in layer k (K)
 ,qek(npnts)       & ! mixing ratio of cloud environment in layer k (kg/kg)
 ,thpk(npnts)      & ! parcel potential temperature in layer k (K)
 ,qpk(npnts)       & ! parcel mixing ratio in layer k (kg/kg)
 ,qsek(npnts)      & ! saturation mixing ratio of cloud environment 
                     ! in layer k (kg/kg)
 ,dqsk(npnts)      & ! gradient of saturation mixing ratio with potential 
                     ! temperature for the cloud environment of layer k
                     ! (kg/kg/K)
 ,thekp1(npnts)    & ! potential temperature of cloud environment in layer k+1
 ,qekp1(npnts)     & ! mixing ratio of cloud environment in layer k+1 (kg/kg)
 ,qsekp1(npnts)    & ! saturation mixing ratio of cloud environment
                     ! in layer k+1 (kg/kg)
 ,dqskp1(npnts)    & ! gradient of saturation mixing ratio with potential  
                     ! temperature for the cloud environment in layer k+1
                     ! (kg/kg/K)
 ,ekp14(npnts)     & ! Entrainment coefficient at level k+1/4 multiplied by 
                     ! appropriate layer thickness
 ,ekp34(npnts)     & ! Entrainment coefficient at level k+3/4 multiplied by 
                     ! appropriate layer thickness
 ,exkp1(npnts)     & ! Exner ratio at mid-point of layer k+1
 ,exk(npnts)       & ! Exner ratio at mid-point of layer k
 ,pkp1(npnts)      & ! pressure at mid-point of layer k+1   (Pa)
 ,pk(npnts)        & ! pressure at mid-point of layer k (Pa)
 ,xsbmin(npnts)      ! Threshold buoyancy for forced detrainment
                     ! Function of delta P


!-----------------------------------------------------------------------
! Variables which are input and output
!-----------------------------------------------------------------------
LOGICAL, INTENT(INOUT) :: &
  bgmkp1(npnts)        & !IN Mask for parcels which are saturated in layer k+1
                         ! after entrainment and latent heating
                         !OUT Mask for parcels which are saturated in layer k+1
                         ! after forced detrainment
 ,bgmk(npnts)            !IN Mask for parcels which are saturated in layer k
                         !OUT Modified by thetar call.

REAL,INTENT(INOUT) ::  &

  thpkp1(npnts)        & ! IN parcel potential temperature in layer k+1 (K)
                         ! after entrainment and latent heating
                         ! OUT parcel potential temperature in layer k+1 (K)
                         ! after forced detrainment
 ,qpkp1(npnts)         & ! IN parcel mixing ratio in layer k+1 (kg/kg)
                         ! after entrainment and latent heating
                         ! OUT parcel mixing ratio in layer k+1 (K)
                         ! after forced detrainment
 ,xsqkp1(npnts)          ! IN Excess water in parcel after lifting layer k to
                         !  k+1 after entrainment and latent heating (kg/kg)
                         ! OUT Excess water in parcel after lifting layer k to
                         !  k+1 after  forced detrainment 

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------

REAL, INTENT(OUT) :: &
  thrk(npnts)        & ! parcel detrainment potential temperature in layer k (K)
 ,qrk(npnts)         & ! parcel detrainment mixing ratio in layer k (kg/kg)
 ,deltak(npnts)        ! parcel forced detrainment rate layer k multiplied
                       ! by appropriate layer thickness


!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::        & 
  i               & ! loop counter
 ,nredo             ! Number of points for which forced detrainment calculation
                    ! must be, as the processes either causes the parcel to
                    ! become saturated or sub-saturated.

LOGICAL ::       &
  bdetk(npnts)   & ! Mask for parcels which are undergoing forced detrainment
                   ! in their ascent from layer k to k+1 
 ,bgkp1w(npnts)  & ! Mask for saturation in layer k+1

 ,brecal(npnts)    ! Mask for those points at which the detrainment
                   !  calculation needs repeating
REAL ::          &
  epss             ! (1+EKP14)*(1+EKP34)

REAL ::           &
  xsqr(npnts)     & ! Excess parcel water vapour during detrainment (kg/kg)
 ,thpkp1w(npnts)  & ! temporary storage for parcel potential temperature (K)
 ,qpkp1w(npnts)   & ! temporary storage for parcel mixing ratio (kg/kg)
 ,xsqk1w(npnts)   & ! temporary storage for parcel excess water vapour
 ,tt(npnts)         ! temporary store for temperature for the calculation of
                    ! saturated mixing ratio (K)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('DETRAIN_4A5A',zhook_in,zhook_handle)

DO i=1,npnts

  !----------------------------------------------------------------------
  ! At start of routine forced detrainment done at all points so 
  ! set array bdetk equal to .TRUE.
  ! Set forced detrainment rate equal to zero.
  !----------------------------------------------------------------------

  bdetk(i) = .TRUE.
  deltak(i) = 0.0

  !-----------------------------------------------------------------------
  !   Save the current values of qpkp1, xsqkp1 and  bgmkp1
  !-----------------------------------------------------------------------

  thpkp1w(i) = thpkp1(i)
  qpkp1w(i) = qpkp1(i)
  xsqk1w(i) = xsqkp1(i)
  bgkp1w(i) = bgmkp1(i)

  !-----------------------------------------------------------------------
  !   Add the excess water vapour into the detraining parcels
  !-----------------------------------------------------------------------

  qpkp1(i) = qpkp1(i) + xsqkp1(i)
END DO

! ----------------------------------------------------------------------
!   Calculate the ensemble average potential temperature in layer k+1
!   at the points where detrainment is taking place
!
!   SUBROUTINE THP_DET
!
!   UM Documentation paper P27
!   Section (6), equation (28)
! ----------------------------------------------------------------------

! DEPENDS ON: thp_det_4a5a
CALL thp_det_4a5a (npnts,bgmkp1,bdetk,thekp1,qpkp1,qekp1,qsekp1,dqskp1, &
                   xsbmin,thpkp1)

! ---------------------------------------------------------------------
!   Check to see if sufficient excess water vapour in the initial dry
!   ascent to allow parcel to be saturated in layer k+1 after
!   forced detrainment
!
!   UM Documentation paper P27
!   Section (6), equation (29)
!
!   NOTE : Only allow parcel to be saturated in layer k+1 if
!          saturated initially.  It is possible for small
!          supersaturations to if subroutine LATENT_H causes
!          parcel to be come unsaturated.  In this case treat
!          the parcel as unsaturated in layer k+1
! ---------------------------------------------------------------------

!-----------------------------------------------------------------------
!   Calculate the excess water vapour in layer k+1 and Recalculate
!   BGMKP1 and QPKP1.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Convert potential temperature to temperature and calculate
! pressure of layer k for calculation of saturated mixing ratio
!-----------------------------------------------------------------------

DO i = 1,npnts
  tt(i) = thpkp1(i)*exkp1(i)
END DO

! DEPENDS ON: qsat
CALL qsat (xsqkp1,tt,pkp1,npnts)

DO i=1,npnts
  xsqkp1(i) = qpkp1(i) - xsqkp1(i)

  brecal(i) = bgmkp1(i)

  !----------------------------------------------------------------------
  ! Only allow parcel to be saturated in initial BGMKP1 = .TRUE.
  ! (stored in BRECAL at this point)
  !----------------------------------------------------------------------

  IF ( bgmk(i) .OR.( (xsqkp1(i)  >   0.) .AND. brecal(i) ) ) THEN
    bgmkp1(i) = .TRUE.
  ELSE
    bgmkp1(i) = .FALSE.
    xsqkp1(i) = 0.0
  END IF

  qpkp1(i) = qpkp1(i) - xsqkp1(i)

  ! ----------------------------------------------------------------------
  ! Recalculate the ensemble average potential temperature at points
  ! where the ensemblehas become unsaturated
  !
  !   UM Documentation paper P27
  !   Section (6), equation (28)
  ! ----------------------------------------------------------------------

  brecal(i) = bdetk(i) .AND. brecal(i) .AND. .NOT.bgmkp1(i)
END DO

! DEPENDS ON: thp_det_4a5a
CALL thp_det_4a5a (npnts,bgmkp1,brecal,thekp1,qpkp1,qekp1,qsekp1,dqskp1, &
                   xsbmin,thpkp1)

! ----------------------------------------------------------------------
! Because of the removal of latent heating, the new parcel potential
! temperature may be lower than its value before the detrainment
! calculation. In this case abandon the detrainment calculation.
! ----------------------------------------------------------------------

DO i=1,npnts
  bdetk(i) = thpkp1(i)  >   thpkp1w(i)
END DO

! ----------------------------------------------------------------------
!   Calculate the potential temperature and mxing ratio of detraining
!   air and the excess water vapour condensed from detraining air.
!
!   UM Documentation paper P27
!   Section (6), equation (26)
! ----------------------------------------------------------------------

! DEPENDS ON: thetar_4a5a
CALL thetar_4a5a (npnts, bwkp1,thek,qek,qpk,qsek,dqsk,exk,pk,     &
                  bgmk,thrk,qrk,xsqr)

! ----------------------------------------------------------------------
!   Calculate the detrainment rate, DELTAK.
!
!   UM Documentation paper P27
!   Section (6), equation (31)
! ----------------------------------------------------------------------

! DEPENDS ON: det_rate_4a5a
CALL det_rate_4a5a (npnts,bwkp1,bdetk,thrk,xsqr,thpk,thek,thekp1,           &
                    xsqkp1,thpkp1,ekp14,ekp34,exk,exkp1,deltak)

nredo = 0

! ----------------------------------------------------------------------
! Add water vapour which was removed from detraining air into XSQKP1
!
!   UM Documentation paper P27
!   Section 86), equation (11C)
! ----------------------------------------------------------------------
!
DO i=1,npnts

  epss = (1.0+ekp14(i))*(1.0+ekp34(i))

  IF (bdetk(i)) THEN           
    xsqkp1(i) = xsqkp1(i) + (deltak(i)*xsqr(i)/(epss*(1.0-deltak(i))))  
  END IF

  ! ----------------------------------------------------------------------
  !   If the excess water vapour in layer k+1 is less than zero
  !   i.e. the parcel has become unsaturated through the forced
  !   detrainment process than abandon the calculation.
  ! ----------------------------------------------------------------------

  brecal(i) = bgmkp1(i)
 
  bgmkp1(i) = xsqkp1(i)  >   0.0

  brecal(i) = bdetk(i) .AND. brecal(i) .AND. .NOT.bgmkp1(i)

  IF (brecal(i)) THEN
    qpkp1(i)  = qpkp1(i) + xsqkp1(i)                              &
                   - (deltak(i)*xsqr(i)/(epss*(1.0-deltak(i))))
    xsqkp1(i) = 0.0
  END IF

  !----------------------------------------------------------------------
  ! Count points at which detrainment calculation needs repeating
  !----------------------------------------------------------------------

  IF (brecal(i)) nredo = nredo + 1
END DO

! ---------------------------------------------------------------------
!   Repeat calculation of parcel potential temperature, detrainment
!   rate and excess parcel water if the parcel becomes unsaturated
!   in layer k+1 after forced detrainment.
! ---------------------------------------------------------------------

IF (nredo  >   0) THEN

!----------------------------------------------------------------------
!  Calculate new parcel potential temperature in layer k+1
!  after forced detrainment
!----------------------------------------------------------------------

! DEPENDS ON: thp_det_4a5a
  CALL thp_det_4a5a (npnts,bgmkp1,brecal,thekp1,qpkp1,qekp1,qsekp1,dqskp1, &
                     xsbmin,thpkp1)

!----------------------------------------------------------------------
!  Check if forced detrainment still possible and reset recalculation
!  mask to FALSE if it is not
!----------------------------------------------------------------------

  DO i=1,npnts
    IF (brecal(i)) THEN
      bdetk(i) = thpkp1(i)  >   thpkp1w(i)
      brecal(i) = bdetk(i)
    END IF
  END DO

!----------------------------------------------------------------------
!  Recalculate forced detrainment rate
!----------------------------------------------------------------------

! DEPENDS ON: det_rate_4a5a
  CALL det_rate_4a5a (npnts,bwkp1,brecal,thrk,xsqr,thpk,thek,thekp1, &
                      xsqkp1,thpkp1,ekp14,ekp34,exk,exkp1,deltak)

!----------------------------------------------------------------------
!  Recalculate excess water vapour in layer k+1
!  after forced detrainment
!----------------------------------------------------------------------

  DO i=1,npnts
    IF (brecal(i)) THEN
      epss = (1.0+ekp14(i))*(1.0+ekp34(i))
      xsqkp1(i) = xsqkp1(i) + (deltak(i)*xsqr(i)/(epss*(1.0-deltak(i)))) 
    END IF
  END DO

END IF

! ----------------------------------------------------------------------
!   Make sure that the detrainment rate is between 0 and 1
!
!   If <0 then no detrainment occurs and original values are
!   restored.
!
!   If >1 then sset to 1 and  THRK = THPK, QRK = QPK and values
!   in layer k+1 arerestored.  Although these are not used
!   in any thermodynamic calculation they are used to specify
!   cloud top in subroutine CONRAD.
! ----------------------------------------------------------------------

DO i=1,npnts

  IF (bdetk(i)) THEN

    IF (deltak(i) <= 0.0) THEN
      bdetk(i) = .FALSE.
      thpkp1(i) = thpkp1w(i)
      qpkp1 (i) = qpkp1w(i)
      xsqkp1(i) = xsqk1w(i)
      bgmkp1(i) = bgkp1w(i)
      deltak(i) = 0.0
    ELSE IF (deltak(i) >  1.0) THEN
      deltak(i) = 1.0
      thrk(i) = thpk(i)
      qrk(i) = qpk(i)
      thpkp1(i) = thpkp1w(i)
      qpkp1 (i) = qpkp1w(i)
      xsqkp1(i) = xsqk1w(i)
      bgmkp1(i) = bgkp1w(i)
    END IF

  ELSE

    thpkp1(i) = thpkp1w(i)
    qpkp1 (i) = qpkp1w(i)
    xsqkp1(i) = xsqk1w(i)
    bgmkp1(i) = bgkp1w(i)
    deltak(i) = 0.0

  END IF    ! test on bdetk
END DO

IF (lhook) CALL dr_hook('DETRAIN_4A5A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE detrain_4a5a
