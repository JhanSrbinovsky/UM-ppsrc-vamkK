! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE Eq_pot_temp( NumLevs,      &  ! in
                        TFields,      &  ! in
                        PFields,      &  ! in
                        EqPotTemp,    &  ! inout
                        ErrorStatus )    ! inout

! Description:
!  Routine to calculate equivalent potential temperature (theta)
!  on model levels given pressure and temperature model level fields.
!
! Method:
!   1) Allocate set of output fields (EqPotTemp) to hold
!      equivalent potential temperature on each model level of interest
!   2) Loop over all gridpoints and model levels. Calculate:
!
!         i) the potential temperature at this point
!              theta = temp * (p0/pressure)**0.286
! For later vn:
!              theta = temp * (Pref/pressure)**Kappa
!        ii) the saturated vapour pressure (es) for each gridpoint using
!            the Clausius-Clapeyron equation
!               es = eso * EXP{(LC/RV)*(( 1/T0 ) - ( 1/Temp) )}
!                where eso=611 Pa at T0=273 K
!                (see e.g. I. James (1994): "Introduction to Circulating
!                 atmospheres", Cambridge University Press)
!       iii) the saturated mixing ratio (ws) for each gridpoint using
!                ws = 0.622 * ( es / ( pressure-es ) )
! For later vn:
!                ws = (R/RV) * ( es / ( pressure-es ) )
!        iv) the equivalent potential temperature (theta_e) using
!            theta_e = theta * EXP( (L*ws)/(CP*Temp) )
! For later vn:
!            theta_e = theta * EXP( (LC*ws)/(CP*Temp) )
!
!   3) The values of equivalent potential temperature are then put into the
!      fields EqPotTemp and returned to the main program.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_mod, ONLY:          &
  PP_Header_type,          &
  PP_Field_type
USE Err_Mod, ONLY:         &
  StatusOK,                &
  StatusWarning
USE FldCodes_mod, ONLY:    &
  ST_MnCATPt, MO8_MnCATPt, &
  PP_MnCATPt, VC_MnCATPt
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs                          ! No. levels
TYPE(PP_Field_type), INTENT(IN)   :: TFields(NumLevs)   ! temperature (K),
TYPE(PP_Field_type), INTENT(IN)   :: PFields(NumLevs)   ! pressure (Pa) &
TYPE(PP_Field_type), INTENT(INOUT):: EqPotTemp(NumLevs) ! equivalent
                                                        ! potential temp
                                                        ! on theta levels
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Eq_pot_temp"
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! For later vn USE atmos_constants_mod
! For later vn USE water_constants_mod
! For later vn #include "rv.h"

INTEGER         :: i,j,k         ! Loop counters

! To be replaced in later vn
REAL, PARAMETER :: Cp  = 1004.0   ! Specific heat at constant pressure in J
REAL, PARAMETER :: L   = 2.5E6    ! Latent heat of condensation in Jkg-1
REAL, PARAMETER :: Rv  = 461.5    ! Gas constant for water vapour in Jkg-1
REAL, PARAMETER :: eso = 611.0    ! Sat. vapour pressure at reference tempe
REAL, PARAMETER :: T0  = 273.0    ! Reference temperature in K
REAL, PARAMETER :: p0  = 100000.0 ! Standard pressure in Pa

! For later vn REAL, PARAMETER :: eso = 611.0  ! Sat. vapour pressure in Pa at
! For later vn                                 !  reference temp T0=273.0
! For later vn REAL, PARAMETER :: T0  = 273.15 ! Reference temperature in K

REAL            :: x             ! exponiential part of saturated vapour
                                 !   pressure equation
REAL            :: es            ! saturated vapour pressure in Pa
REAL            :: ws            ! saturated mixing ratio in
REAL            :: theta         ! Potential temperature in K
REAL            :: theta_e       ! Equivalent potential temperature in K

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF


!Check that temperature and pressure fields are on the same grid
DO k = 1, NumLevs
  IF ( (TFields(k) % Hdr % NumCols   /=                             &
        PFields(k) % Hdr % NumCols) .OR.                            &
       (TFields(k) % Hdr % NumRows   /=                             &
        PFields(k) % Hdr % NumRows) .OR.                            &
       (TFields(k) % Hdr % ZerothLon /=                             &
        PFields(k) % Hdr % ZerothLon) .OR.                          &
       (TFields(k) % Hdr % ZerothLat /=                             &
        PFields(k) % Hdr % ZerothLat) ) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus,                         &
                  "Model level fields supplied are on different grids" )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF
END DO


!Set up output field

!Copy header from temperature fields
!(doesn't matter that fieldcode is correct as this field won't be output)
EqPotTemp(1:NumLevs) % Hdr           = TFields(1:NumLevs) % Hdr
EqPotTemp(1:NumLevs) % Hdr % BMDI    = RMDI

DO k = 1, NumLevs

  IF ( ASSOCIATED( EqPotTemp(k) % RData ) ) THEN
    DEALLOCATE( EqPotTemp(k) % RData )
    NULLIFY( EqPotTemp(k) % RData )
  END IF

  IF ( .NOT. ASSOCIATED( EqPotTemp(k) % RData ) ) THEN
    ALLOCATE( EqPotTemp(k) % RData(EqPotTemp(1) % Hdr % NumCols,     &
                                   EqPotTemp(1) % Hdr % NumRows) )
  END IF

END DO


! Calculate the potential temperature, saturated vapour pressure,
! saturated mixing ratio then equivalent potential temperature for each
! gridpoint & level where possible. Set equivalent potential temperature
! to missing data indictor if the terrain makes calculation impossible.

DO k = 1, NumLevs
  DO j = 1, TFields(1) % Hdr % NumRows
    DO i = 1, TFields(1) % Hdr % NumCols
      IF ((TFields(k) % RData(i,j) /= TFields(1) % Hdr % BMDI) .AND.  &
          (PFields(k) % RData(i,j) /= PFields(1) % Hdr % BMDI)) THEN
        !Calculation is possible

        !Calculate the potential temperature at this point
        theta = TFields(k) % RData(i,j) * &
                ( ( p0/PFields(k) % RData(i,j) )**0.286 )
! For later vn     theta = TFields(k) % RData(i,j) * &
! For later vn             ( ( Pref/PFields(k) % RData(i,j) )**Kappa )

        !Calculate the saturated vapour pressure for this point
        x = (L/RV) * ( ( 1/T0 ) - ( 1/TFields(k) % RData(i,j) ) )
! For later vn      x = (LC/RV) * ( ( 1/T0 ) - ( 1/TFields(k) % RData(i,j) ) )
        es = eso * EXP(x)

        !Calculate the saturated mixing ratio at this point
        ws = 0.622 * ( es / ( PFields(k) % RData(i,j)-es ) )
! For later vn      ws = (R/RV) * ( es / ( PFields(k) % RData(i,j)-es ) )

        !Calculate the equivalent potential temperature at this point
        !and put into the output field element

        theta_e = theta * EXP( (L*ws)/(CP*TFields(k) % RData(i,j)) )
        EqPotTemp(k) % RData(i,j) = theta_e

     ELSE
        !Terrain makes calculation impossible
        EqPotTemp(k) % RData(i,j) = EqPotTemp(k) % Hdr % BMDI

     END IF

    END DO
  END DO
END DO

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Eq_pot_temp
