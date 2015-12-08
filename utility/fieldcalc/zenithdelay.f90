! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the Zenith Total Delay
!
subroutine zenithdelay(NumLevs, num_wetlevs,                 &
                       pressure_on_theta,                    &
                       exner_on_rho,                         &
                       rho_heights,                          &
                       q_on_theta,                           &
                       zenith_total_delay,                   &
                       errorstatus                           )



! Description: Calculates the zenith total delay at the surface
!
! Method: Based on the algorithm used in the OPS code.
!
! The following steps occur:
!
! 1. Calculate signal delay above model top
!
! 2. Calculate zenith total delay (ZTD)
!
! 3. Apply height correction to observation
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!

USE atmos_constants_mod, ONLY:  &
  r, cp, pref, kappa, c_virtual, repsilon

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type

USE FldCodes_Mod, ONLY:   &
  ST_ztd,     &
  PP_ztd,     &
  MO8_ztd,    &
  LV_special, &
  VC_Surface

USE Err_Mod, ONLY:        &
  StatusOK

USE earth_constants_mod, ONLY: g

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN)              :: NumLevs
INTEGER, INTENT(IN)              :: num_wetlevs
TYPE(PP_Field_type), INTENT(IN)  :: pressure_on_theta(numlevs)
TYPE(PP_Field_type), INTENT(IN)  :: exner_on_rho(numlevs+1)
TYPE(PP_Field_type), INTENT(IN)  :: rho_heights(numlevs+1)
TYPE(PP_Field_type), INTENT(IN)  :: q_on_theta(NumLevs)
TYPE(PP_Field_type), INTENT(OUT) :: zenith_total_delay
INTEGER, INTENT(INOUT)           :: ErrorStatus


REAL              :: exner_on_theta
REAL, ALLOCATABLE :: p_on_toprho(:,:)
REAL              :: T
REAL              :: Nwet
REAL              :: Ndry

REAL :: temp1,temp2

INTEGER :: row_len
INTEGER :: n_rows
! Local general purpose e.g. loop counters
INTEGER :: i
INTEGER :: j
INTEGER :: k

CHARACTER(LEN=*), PARAMETER :: RoutineName = "ZenithDelay"


! Included constants/Parameters from UM constants files

! Parameters set up for the header
REAL, PARAMETER :: nalpha=77.6
REAL, PARAMETER :: nbeta=3.73E5
REAL, PARAMETER :: govercp=g/cp

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! If Zenith delay is already allocated, deallocate before reallocating
IF ( ASSOCIATED( zenith_total_delay % RData ) )then
  DEALLOCATE( zenith_total_delay % RData )
END IF

! Derive size of grid from header of pressure_on_theta
row_len = pressure_on_theta(1) % Hdr % NumCols
n_rows  = pressure_on_theta(1) % Hdr % NumRows

! Allocate zenith_total_delay
ALLOCATE(zenith_total_delay % RData(row_len,n_rows))

! Allocate locally derived items for correct grid size
ALLOCATE(p_on_toprho(row_len,n_rows))

! Zero reals.
exner_on_theta           =  0.0
T                        =  0.0
Nwet                     =  0.0
Ndry                     =  0.0

! Set up a header for Zenith Delay
! Copy header from an input argument on the same grid as ZTD
zenith_total_delay % Hdr = pressure_on_theta(1) % Hdr

! Overwrite those items that need to be different.
zenith_total_delay % Hdr % STCode   = ST_ztd
zenith_total_delay % Hdr % PPCode   = PP_ztd
zenith_total_delay % Hdr % MO8Level = LV_Special
zenith_total_delay % Hdr % MO8Type  = MO8_ztd
zenith_total_delay % Hdr % LBVC     = VC_Surface

p_on_toprho(:,:) = (exner_on_rho(NumLevs+1) % RData(:,:) ** (1/kappa)) * Pref

! Calculate the correction for delay above top of model
zenith_total_delay % RData(:,:) = 1E-6 * Nalpha /                            &
                                  100.0 * R * (p_on_toprho(:,:) / g)

IF (NumLevs > Num_wetlevs) THEN
  DO k=NumLevs,Num_wetlevs+1,-1
    DO j=1, n_rows
      DO i=1, row_len
        exner_on_theta = (pressure_on_theta(k) % RData(i,j) / Pref) ** kappa
        temp1 = govercp*(rho_heights(k+1) % RData(i,j)-                      &
          rho_heights(k) % RData(i,j))
        temp2 = exner_on_rho(k) % RData(i,j) - exner_on_rho(k+1) % RData(i,j)

        T =  exner_on_theta*temp1/temp2

        ! Calculate the dry refractivity for this level
        Ndry = 0.01*Nalpha*pressure_on_theta(k) % RData(i,j) / T

        ! Calculate zenith delay for this level
        zenith_total_delay % Rdata(i,j) = zenith_total_delay % rdata(i,j) +  &
                        ( 1E-6 * (Ndry) *                                    &
                          ( rho_heights(k+1) % RData(i,j) -                  &
                            rho_heights(k)   % RData(i,j) ) )

      END DO
    END DO
  END DO
ELSE ! Numlevs > num_wetlevs
  DO k=Num_wetlevs,1,-1
    DO j=1, n_rows
      DO i=1, row_len
        Nwet  = 0.0
        exner_on_theta = (pressure_on_theta(k) % RData(i,j) / Pref) ** kappa
        temp1 = govercp*(rho_heights(k+1) % RData(i,j) -                     &
          rho_heights(k) % RData(i,j))
        temp2 = exner_on_rho(k) % RData(i,j) - exner_on_rho(k+1) % RData(i,j)

        T =  exner_on_theta*temp1/temp2
        ! If a wet level calculate the wet refractivity
        T = T / (1.0 + C_virtual * q_on_theta(k) % RData(i,j))
        temp1 = 0.01 * NBeta * pressure_on_theta(k) % RData(i,j) *           &
                q_on_theta(k) % RData(i,j)
        temp2 = ( T**2 *                                                     &
                  ( repsilon + (1-repsilon) * q_on_theta(k) % RData(i,j) ) )
        Nwet  = temp1/temp2

        ! Calculate the dry refractivity for this level
        Ndry = 0.01*Nalpha*pressure_on_theta(k) % RData(i,j) / T

        ! Calculate zenith delay for this level
        zenith_total_delay % Rdata(i,j) = zenith_total_delay % rdata(i,j) +    &
                        ( 1E-6 * (Ndry+Nwet) *                                 &
                          ( rho_heights(k+1) % RData(i,j) -                    &
                            rho_heights(k)   % RData(i,j) ) )

      END DO
    END DO
  END DO
END IF

DEALLOCATE(p_on_toprho)

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

RETURN
END SUBROUTINE zenithdelay

