! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Forced detrainment calculation

MODULE detrain_6a_mod

IMPLICIT NONE

!
! Description:
!   Forced detrainment calculation
!
!   Subroutine thp_det calculates the potential temperature of the parcel
!   in layer k+1 after forced detrainment.
!
!   Subroutine thetar calculates the potential temperature of the air in
!   layer k undergoing forced detrainment.
!
!   Subroutine det_rate calculates the forced detrainment rate of the
!   ensemble in layer k.
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

SUBROUTINE detrain_6a (npnts, pk, pkp1, exk, exkp1,                 &
                       thek, thekp1, qek, qekp1,                    &
                       thpk, qpk,                                   &
                       watldek, watldpk, watldekp1, watldpkp1,      & 
                       Qlkp1, Qfkp1, Frezkp1,                       &
                       ekp14, ekp34, xsbmin,                        &
                       bwk, bwkp1, bgmk, bgmkp1,                    &
                       thpkp1, qpkp1,                               &
                       deltak, thrk, qrk)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE thetar_6a_mod
USE thp_det_6a_mod
USE det_rate_6a_mod

IMPLICIT NONE

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

INTEGER,INTENT(IN) :: npnts         ! Number of points

REAL,INTENT(IN) :: exkp1(npnts)     ! Exner ratio at mid-point of layer k+1
REAL,INTENT(IN) :: exk(npnts)       ! Exner ratio at mid-point of layer k
REAL,INTENT(IN) :: pkp1(npnts)      ! pressure at mid-point of layer k+1 (Pa)
REAL,INTENT(IN) :: pk(npnts)        ! pressure at mid-point of layer k (Pa)
REAL,INTENT(IN) :: thek(npnts)      ! Env. potential temperature in layer k (K)
REAL,INTENT(IN) :: thekp1(npnts)    ! Env. potential temperature in layer k+1
REAL,INTENT(IN) :: qek(npnts)       ! Env. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qekp1(npnts)     ! Env. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: thpk(npnts)      ! Par. potential temperature in layer k (K)
REAL,INTENT(IN) :: qpk(npnts)       ! Par. mixing ratio in layer k (kg/kg)
REAL,INTENT(IN) :: watldek(npnts)   ! Env. water loading in layer k (kg/kg)
REAL,INTENT(IN) :: watldpk(npnts)   ! Par. water loading in layer k (kg/kg)
REAL,INTENT(IN) :: watldekp1(npnts) ! Env. water loading in layer k+1 (kg/kg)
REAL,INTENT(IN) :: watldpkp1(npnts) ! Par. water loading in layer k+1 (kg/kg)
REAL,INTENT(IN) :: Qlkp1(npnts)     ! Amount of condensation to liquid water 
                                    ! in the parcel (kg/kg)
REAL,INTENT(IN) :: Qfkp1(npnts)     ! Amount of deposition to ice water
                                    ! in the parcel (kg/kg)
REAL,INTENT(IN) :: Frezkp1(npnts)   ! Amount of freezing from liquid 
                                    ! to ice water in the parcel (kg/kg)

REAL,INTENT(IN) :: ekp14(npnts)     ! Entrainment coefficient at level k+1/4 
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: ekp34(npnts)     ! Entrainment coefficient at level k+3/4 
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: xsbmin(npnts)    ! Threshold buoyancy for forced 
                                    ! detrainment (K)

LOGICAL,INTENT(IN) :: bwk(npnts)    ! Mask for whether condensate is 
                                    ! liquid in layer k
LOGICAL,INTENT(IN) :: bwkp1(npnts)  ! Mask for whether condensate is 
                                    ! liquid in layer k+1
LOGICAL,INTENT(IN) :: bgmk(npnts)   ! Mask for parcels whcih are saturated
                                    ! in layer k
LOGICAL,INTENT(IN) :: bgmkp1(npnts) ! Mask for parcels which are saturated
                                    ! in layer k+1

!-----------------------------------------------------------------------
! Variables which are input and output
!-----------------------------------------------------------------------
REAL,INTENT(INOUT) :: thpkp1(npnts) ! Par. pot. temperature in layer k+1 (K)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment
REAL,INTENT(INOUT) :: qpkp1(npnts)  ! Par. spec. humidity in layer k+1 (kg/kg)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------

REAL,INTENT(OUT) :: thrk(npnts)     ! pot. temperature of forced detrained
                                    ! parcel in layer k (K)
REAL,INTENT(OUT) :: qrk(npnts)      ! Specific humidity of forced detrained
                                    ! parcel in layer k (kg/kg)
REAL,INTENT(OUT) :: deltak(npnts)   ! Parcel forced detrainment rate in 
                                    ! layer k multiplied by layer thickness

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER :: i  ! loop counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('DETRAIN_6A',zhook_in,zhook_handle)

! Initialise deltak
DO i=1,npnts
  deltak(i) = 0.0
  thrk(i)   = 0.0
  qrk(i)    = 0.0
END DO

! ----------------------------------------------------------------------
!   Calculate the potential temperature and humidity of the parcel 
!   undergoing forced detrainment
! ----------------------------------------------------------------------

CALL thetar_6a (npnts, exk, pk, thek, qek, qpk, watldek, watldpk,         &
                bwk, bgmk, thrk, qrk)

! ----------------------------------------------------------------------
!   Calculate the potential temperature and humidity in layer k+1
!   at the points where detrainment is taking place
! ----------------------------------------------------------------------

CALL thp_det_6a (npnts, exkp1, pkp1, thekp1, qekp1, watldekp1, watldpkp1, &
                 xsbmin, bwkp1, bgmkp1, qpkp1, thpkp1)

!----------------------------------------------------------------------
!  calculate forced detrainment rate
!----------------------------------------------------------------------

CALL det_rate_6a (npnts, qek, qekp1, qpk, qpkp1, qrk, Qlkp1, Qfkp1,       &
                  ekp14, ekp34, deltak)

IF (lhook) CALL dr_hook('DETRAIN_6A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE detrain_6a
END MODULE detrain_6a_mod
