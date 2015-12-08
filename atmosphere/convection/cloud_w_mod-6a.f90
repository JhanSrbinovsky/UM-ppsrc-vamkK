! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Convective cloud microphysics routine
!
MODULE cloud_w_6a_mod

IMPLICIT NONE

! Description:
!   Calculates preciptation produced in lifting parcel from layer k to k+1
!
!   Calls CON_RAD to calculate parameters for radiation calculation
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

SUBROUTINE cloud_w_6a (k, npnts, start_lev,                                &
                       flxkp1, qclpk, qcfpk,                               &
                       thekp1, qekp1, qsekp1,                              &
                       ekp14, ekp34, delpkp1,                              &
                       blowst, bwkp1, bland, bterm, l_q_interact,          &
                       lcbase, lctop,                                      &
                       qclpkp1, qcfpkp1,                                   &
                       tcw, depth, cclwp, lcca,                            &
                       iccb, icct, prekp1, cca, ccwkp1)

USE atmos_constants_mod, ONLY: cp, c_virtual

USE cv_run_mod, Only:                                                   &
    mparwtr, qlmin, fac_qsat, ccw_for_precip_opt

USE earth_constants_mod, ONLY: g

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE con_rad_6a_mod

IMPLICIT NONE

! Subroutine arguments 

!----------------------------------------------------------------------
! Variables that are input
!----------------------------------------------------------------------
INTEGER,INTENT(IN) :: k                   ! Present model layer number
INTEGER,INTENT(IN) :: npnts               ! Vector length

INTEGER,INTENT(IN) :: start_lev(npnts)    ! Level at which convection initiated

REAL,INTENT(IN) :: flxkp1(npnts)    ! Parcel mass flux in layer k+1 (Pa/s)
REAL,INTENT(IN) :: qclpk(npnts)     ! Par. qcl in layer k (kg/kg)
REAL,INTENT(IN) :: qcfpk(npnts)     ! Par. qcf in layer k (kg/kg)
REAL,INTENT(IN) :: thekp1(npnts)    ! Env. pot. temperature in layer k+1 (K)
REAL,INTENT(IN) :: qekp1(npnts)     ! Env. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qsekp1(npnts)    ! Env. saturated specific humidity in 
                                    ! in layer k+1 (kg/kg)
REAL,INTENT(IN) :: ekp14(npnts)     ! Entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: ekp34(npnts)     ! Entrainment coefficient at level k+3/4 
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: delpkp1(npnts)   ! pressure difference across layer k+1 (Pa)

LOGICAL,INTENT(IN) :: blowst(npnts) ! mask for those points at which stability
                                    ! is low enough for convection to occur
LOGICAL,INTENT(IN) :: bwkp1(npnts)  ! mask for parcels which have liquid 
                                    ! condensate in layer k+1
LOGICAL,INTENT(IN) :: bland(npnts)  ! Land/sea mask
LOGICAL,INTENT(IN) :: bterm(npnts)  ! Mask for parcels which terminate 
                                    ! in layer k+1
LOGICAL,INTENT(IN) :: l_q_interact  ! True if PC2 is switched on

!----------------------------------------------------------------------
! Variables that are input and output
!----------------------------------------------------------------------
INTEGER,INTENT(INOUT) :: lcbase(npnts)! Lowest conv. cloud base level
INTEGER,INTENT(INOUT) :: lctop(npnts) ! Lowest conv. cloud top level

REAL,INTENT(INOUT) :: qclpkp1(npnts)! Parcel liquid condensate mixing ratio
                                    ! in layer k+1 (kg/kg)
                                    ! IN:  before precipitation
                                    ! OUT: after precipitation
REAL,INTENT(INOUT) :: qcfpkp1(npnts)! Parcel liquid condensate mixing ratio
                                    ! in layer k+1 (kg/kg)
                                    ! IN:  before precipitation
                                    ! OUT: after precipitation
REAL,INTENT(INOUT) :: tcw(npnts)    ! Total condensed water (kg/m**2/s)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: depth(npnts)  ! Depth of convective cloud (m)
                                    ! IN summed to layer k 
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: cclwp(npnts)  ! Condensed water path (kg/m**2)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: lcca(npnts)   ! Lowest conv. cloud amount (%)

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
INTEGER,INTENT(OUT) :: iccb(npnts)  ! convective cloud base_level
INTEGER,INTENT(OUT) :: icct(npnts)  ! convective cloud top level

REAL,INTENT(OUT) :: prekp1(npnts)   ! precipitation from parcel as it rises 
                                    ! from layer k to k+1 (kg/m**2/s)
REAL,INTENT(OUT) :: cca(npnts)      ! convective cloud amount (%)
REAL,INTENT(OUT) :: ccwkp1(npnts)   ! Total condensate in level k+1 (kg/kg)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

INTEGER :: i                ! Loop counter

REAL :: mparmult            ! Factor used to multiply mparwtr, value between
                            ! 1. and ~3. being 1.0 for deeper clouds.
REAL :: tmp_ccwkp1(npnts)   ! Total condensate in level k+1 before precipitation
                            ! (kg/kg)
REAL :: ccwk(npnts)         ! Total condensate in level k (kg/kg)
                    
REAL :: xmin(npnts)         ! Maxmimum amount of cloud water retained by the 
                            ! parcel before the parcel precipitates (kg/kg)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CLOUD_W_6A',zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! Store convective cloud liquid water before Precipitation
!-------------------------------------------------------------------------------
!
DO i=1,npnts
  tmp_ccwkp1(i) = qclpkp1(i) + qcfpkp1(i)  ! Total condensate in level k+1 
                                           ! before precipitation
  ccwk(i)       = qclpk(i)   + qcfpk(i)    ! Total condensate in level k
  prekp1(i)     = 0.0                      ! initialise precipitation to zero
END DO

!-------------------------------------------------------------------------------
! Calculate Precipitation from layer k+1 and adjust cloud water
!-------------------------------------------------------------------------------

SELECT CASE (ccw_for_precip_opt)

  CASE (4)            ! xmin profile based on qsat with user defined 
                      ! minimum, maximum and qsat scaling.
    DO i=1,npnts
      
      xmin(i) = MAX(MIN(mparwtr, fac_qsat*qsekp1(i)), qlmin)

    END DO

  CASE (3)            ! Manoj's function for congestus
                                          
    DO i=1,npnts

      mparmult = 1.5 + 0.5*TANH((3000.0 -depth(i))/800.0)
      xmin(i)  = MIN(mparwtr*mparmult, fac_qsat*qsekp1(i)) 

      ! If a land point and liquid water in the layer 
      ! increase the minimum cloud water for Precipitation
      ! The reasons for this are an attempt to take some account of 
      ! more aerosols over land leading to more small cloud drops and
      ! therefore more cloud water before Precipitation.  

      IF (bwkp1(i) .AND. bland(i)) xmin(i) =xmin(i)*2.0

      ! limit max value to 0.003 kg/kg
      xmin(i) = MIN(xmin(i),0.003)          

      IF (l_q_interact) THEN   ! PC2
      ! Limit xmin(i)
        xmin(i) = MAX(xmin(i), qlmin)
      END IF

    END DO

  CASE (2)            ! Manoj's function dependent on depth of cloud

    DO i=1,npnts

      mparmult = 2.0 + 1.0*TANH((1500.0 -depth(i))/1000.0)
      xmin(i)  = MIN(mparwtr*mparmult, fac_qsat*qsekp1(i)) 

      ! If a land point and liquid water in the layer 
      ! increase the minimum cloud water for Precipitation
      ! The reasons for this are an attempt to take some account of 
      ! more aerosols over land leading to more small cloud drops and
      ! therefore more cloud water before Precipitation.  

      ! If a land point and liquid water in the layer 
      IF (bwkp1(i) .AND. bland(i)) xmin(i) = xmin(i)*2.0

      IF (l_q_interact) THEN
      ! Limit xmin
        xmin(i) = MAX(xmin(i), qlmin)
      END IF

    END DO

  CASE (1)            ! No test on a critical depth for Precipitation
                      ! Option in use for HadGEM2
    DO i=1,npnts
      xmin(i) = MIN(mparwtr, fac_qsat*qsekp1(i)) 

      ! If a land point and liquid water in the layer 
      IF (bwkp1(i) .AND. bland(i)) xmin(i) = xmin(i)*2.0

      IF (l_q_interact) THEN
      ! Limit xmin
        xmin(i) = MAX(xmin(i), qlmin)
      END IF

    END DO

END SELECT      ! test on value of ccw_for_precip_opt 

DO i=1,npnts

  ! Precipitate if cloud water in the layer > xmin

  IF ( bwkp1(i) ) THEN  !If condensate is liquid
      IF ( qclpkp1(i) > xmin(i) ) THEN
        prekp1(i)  = (qclpkp1(i) - xmin(i)) * flxkp1(i) /g
        qclpkp1(i) = xmin(i)
      END IF
  ELSE !Condensate is frozen
      IF ( qcfpkp1(i) > xmin(i) ) THEN
        prekp1(i)  = (qcfpkp1(i) - xmin(i)) * flxkp1(i) /g
        qcfpkp1(i) = xmin(i)
      END IF
  END IF
  
  !Update the convective cloud water. Do not permit negative CCW.
  ccwkp1(i) = MAX(qclpkp1(i) + qcfpkp1(i), 0.0)

END DO

!-------------------------------------------------------------------------------
! Calculate convective cloud base, convective cloud top, total condensed
! water/ice and convective cloud amount
!
! SUBROUTINE CON_RAD - now called after rainout has occurred
! UM Documentation paper 27
! Section (9)
!-------------------------------------------------------------------------------

CALL con_rad_6a (k, npnts, start_lev,                        &
                 ccwk, ccwkp1, tmp_ccwkp1, flxkp1, delpkp1,  &
                 l_q_interact, bterm,                        &
                 tcw, cclwp, lcca, lcbase, lctop,            &
                 iccb, icct, cca)

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('CLOUD_W_6A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE cloud_w_6a
END MODULE cloud_w_6a_mod
