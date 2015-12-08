! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate contrail prediction diagnostic

SUBROUTINE Contrail( NumLevs,      &  ! in
                     PStar,        &  ! in
                     PFields,      &  ! in
                     TFields,      &  ! in
                     ZFields,      &  ! in
                     ContrailL,    &  ! inout
                     ContrailU,    &  ! inout
                     ErrorStatus )    ! inout

! Description:
!   This subroutine uses pressure, temperature and height to
!   calculate an upper and lower limit within which contrails are
!   expected to form.
!
! Method:
!   The aim is to find the intersection of the temperature profile and
!   the Mintra -14 deg 'environment curve' at each point.  At the
!   the intersection: Pressure (p_m) = (( Temperature)/A0)**(1/A1)
!   Also,             Temperature    = t(k)*(p_m/p(k))**(R*lapse/g)
!   So  p_m = ( (t(k)/A0)*p(k)**(-R*lapse/g) )**(1/(A1-R*lapse/g))
!   The routine loops through the input temperature fields to find the
!   2 intervals where the curve is crossed.  Interpolation is done
!   within these intervals to find the 2 intersections with the curve,
!   which are the upper and lower contrail limits.
!   Some information on Mintra lines can be found in the Forecasters
!   Reference Book.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3

USE atmos_constants_mod, ONLY: r

USE IO_mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_mod, ONLY:   &
  ST_ContrL, MO8_ContrL,  &
  ST_ContrU, MO8_ContrU,  &
  PP_Contr,  LV_Special,  &
  VC_Lower,  VC_Upper

USE earth_constants_mod, ONLY : g

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
TYPE(PP_Field_type), INTENT(IN) :: PStar            ! Surface pressure
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs) ! Pressure, temp
TYPE(PP_Field_type), INTENT(IN) :: TFields(NumLevs) !   and height
TYPE(PP_Field_type), INTENT(IN) :: ZFields(NumLevs) !   on theta levels

TYPE(PP_Field_type), INTENT(INOUT) :: ContrailL  ! Lower CT ICAO hght
TYPE(PP_Field_type), INTENT(INOUT) :: ContrailU  ! Upper CT ICAO hght
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Contrail"

REAL, PARAMETER :: G_over_R = G / R
REAL, PARAMETER :: A0 = 137.14816 ! Constants for approximating the
REAL, PARAMETER :: A1 = 0.046822  ! Mintra - 14 deg "environment curve"

! Local Variables:
INTEGER :: i, j, k, lev     ! Loop counters
INTEGER :: LevLwr (PStar%Hdr%NumCols,PStar%Hdr%NumRows)
INTEGER :: LevUpr (PStar%Hdr%NumCols,PStar%Hdr%NumRows)
REAL :: lapse               ! Lapse rate in layer containing
                            ! Intersection
REAL :: temp, exponent      ! temporary scalars
REAL :: mintra_k  (PStar%Hdr%NumCols,PStar%Hdr%NumRows) ! Mintra line
REAL :: mintra_km1(PStar%Hdr%NumCols,PStar%Hdr%NumRows) !   - 14 degrees
REAL :: p_m       (PStar%Hdr%NumCols,PStar%Hdr%NumRows) ! P at
                                                        ! Intersection
TYPE(PP_Field_type) :: PLwr ! Max pressure for contrail formation
TYPE(PP_Field_type) :: PUpr ! Min pressure for contrail formation

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Set up headers & allocate memory
PLwr % Hdr = PStar % Hdr
PLwr % Hdr % LBVC = VC_Lower
PLwr % Hdr % MO8Level = LV_Special
PUpr % Hdr = PLwr % Hdr
PUpr % Hdr % LBVC = VC_Upper
ALLOCATE( PLwr % RData(PLwr % Hdr % NumCols,  &
                       PLwr % Hdr % NumRows) )
ALLOCATE( PUpr % RData(PUpr % Hdr % NumCols,  &
                       PUpr % Hdr % NumRows) )
PLwr % RData(:,:) = 0.0
PUpr % RData(:,:) = 0.0

! Check lowest level for contrails
LevLwr(:,:) = 0
LevUpr(:,:) = 0
Mintra_k = A0*(PFields(1) % RData **A1) - TFields(1) % RData
WHERE( Mintra_k > 0 )
  LevLwr = 1              ! Contrails forming at surface
  PLwr % RData = PStar % RData
END WHERE

! Find the 2 intervals containing the Mintra -14deg curve
DO k = 2,NumLevs
  Mintra_km1 = Mintra_k
  Mintra_k = A0*(PFields(k) % RData **A1) - TFields(k) % RData

  WHERE( (LevLwr > 0)                  .AND. &
         (PFields(k) % RData > 5000.0) .AND. &
         (Mintra_k * Mintra_km1 <= 0.0) )
    LevUpr = k            ! Top of contrail layer (bottom already found)
  END WHERE
  WHERE( (LevLwr == 0)                 .AND. &
         (PFields(k) % RData > 5000.0) .AND. &
         (Mintra_k * Mintra_km1 <= 0.0) )
    LevLwr = k            ! Bottom of contrail layer
  END WHERE
END DO

! Interpolate within the intervals to find the curve exactly
DO lev = 1,2              ! Lwr then Upr
  DO j = 1,PLwr % Hdr % NumRows
    DO i = 1,PLwr % Hdr % NumCols

      k = LevLwr(i,j)
      IF ( k > 1 ) THEN
        lapse = (TFields(k-1) % RData(i,j) - TFields(k) % RData(i,j)) &
                 / (ZFields(k) % RData(i,j) - ZFields(k-1) % RData(i,j))
        temp  = (PFields(k) % RData(i,j)**lapse) &
                              * ((A0/TFields(k) % RData(i,j))**G_over_R)
        exponent = 1./(lapse - A1*G_over_R)
        p_m(i,j) = temp**exponent

      END IF

    END DO
  END DO
  IF ( lev == 1 ) THEN
    ! First intersection is lower limit
    WHERE( LevLwr > 1 )
      PLwr % RData = p_m
    END WHERE
  ELSE
    ! Second intersection is upper limit
    WHERE( LevUpr > 1 )
      PUpr % RData = p_m
    END WHERE
  END IF
  LevLwr = LevUpr
END DO

! Convert pressure to ICAO height
! DEPENDS ON: icaoheight
CALL ICAOHeight( PLwr, ContrailL, ErrorStatus )
! DEPENDS ON: icaoheight
CALL ICAOHeight( PUpr, ContrailU, ErrorStatus )
ContrailL % Hdr % PPCode  =  PP_Contr
ContrailL % Hdr % MO8Type = MO8_ContrL
ContrailL % Hdr % STCode  =  ST_ContrL
ContrailU % Hdr % PPCode  =  PP_Contr
ContrailU % Hdr % MO8Type = MO8_ContrU
ContrailU % Hdr % STCode  =  ST_ContrU
DEALLOCATE( PLwr % RData )
DEALLOCATE( PUpr % RData )

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Contrail

