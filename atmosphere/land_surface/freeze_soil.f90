! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE FREEZE_SOIL -----------------------------------------

!
! Subroutine Interface:
      SUBROUTINE FREEZE_SOIL (NPNTS,NSHYD,B,DZ                          &
     &,                       SATHH,SMCL,TSOIL,V_SAT,STHU,STHF)

      USE water_constants_mod, ONLY: rho_water
      USE conversions_mod, ONLY: zerodegc

      USE ereport_mod, ONLY : Ereport

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Description:
!     Calculates the unfrozen and frozen water within a soil layer
!     as a fraction of saturation.                          (Cox, 6/95)
!
! Documentation : UM Documentation Paper 25
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: land surface
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!


! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NSHYD                ! IN Number of soil layers.


!   Array arguments with intent(IN) :

      REAL                                                              &
     & B(NPNTS)                                                         &
                            ! IN Clapp-Hornberger exponent.
     &,DZ(NSHYD)                                                        &
                            ! IN Thicknesses of the soil layers (m).
     &,SATHH(NPNTS)                                                     &
                            ! IN Saturated soil water pressure (m).
     &,SMCL(NPNTS,NSHYD)                                                &
                            ! IN Soil moisture content of
                            !    layers (kg/m2).
     &,TSOIL(NPNTS,NSHYD)                                               &
                            ! IN Sub-surface temperatures (K).
     &,V_SAT(NPNTS)         ! IN Volumetric soil moisture
                            !    concentration at saturation
                            !    (m3 H2O/m3 soil).

!   Array arguments with intent(OUT) :
      REAL                                                              &
     & STHF(NPNTS,NSHYD)                                                &
                            ! OUT Frozen soil moisture content of
                            !     the layers as a fraction of
                            !     saturation.
     &,STHU(NPNTS,NSHYD)    ! OUT Unfrozen soil moisture content of
                            !     the layers as a fraction of
                            !     saturation.

! Local scalars:
      INTEGER                                                           &
     & I,J,N                ! WORK Loop counters.

! Local arrays:
      REAL                                                              &
     & SMCLF(NPNTS,NSHYD)                                               &
                            ! WORK Frozen moisture content of the
                            !      soil layers (kg/m2).
     &,SMCLU(NPNTS,NSHYD)                                               &
                            ! WORK Unfrozen moisture content of the
                            !      soil layers (kg/m2).
     &,SMCLSAT(NPNTS,NSHYD)                                             &
                            ! WORK The saturation moisture content of
                            !      the layers (kg/m2).
     &,TMAX(NPNTS)                                                      &
                            ! WORK Temperature above which all water is
                            !      unfrozen (Celsius)
     &,TSL(NPNTS,NSHYD)     ! WORK Soil layer temperatures (Celsius).

      REAL WORK1(NPNTS,NSHYD)  ! Stores (SMCL(I,N)/SMCLSAT(I,N))**(B(I))
      
      LOGICAL L_DZ_ZERO     ! Warn if DZ contains zero values
      
! Error reporting
      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'freeze_soil'
      CHARACTER(LEN=256)       :: message
      INTEGER                  :: errorstatus

      REAL tiny_0, small_value

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
       
! C_PERMA start

      ! Specific heat capacity of water vapour (J/kg/K)
      REAL,PARAMETER:: HCAPV=1850.0

      ! Specific heat capacity of water (J/kg/K)
      REAL,PARAMETER:: HCAPW=4180.0

      ! Specific heat capacity of ice (J/kg/K)
      REAL,PARAMETER:: HCAPI=2100.0

      ! Density of ice (kg/m3)
      REAL,PARAMETER:: RHO_ICE=917

      ! Rate of change of ice potential with temperature
      ! RHO_ICE*LF/ZERODEGC*1/(RHO_WATER*G) (m/K)
      REAL,PARAMETER:: DPSIDT=114.3

! C_PERMA end

       IF (lhook) CALL dr_hook('FREEZE_SOIL',zhook_in,zhook_handle)
      tiny_0=TINY(0.0)
      small_value=EPSILON(0.0)

      WRITE(6,*)                                                        &
     &'Freeze_Soil: Calculating unfrozen and frozen soil water fracs.'           

      L_DZ_ZERO=.false.
            
      ! Warn if soil layer thickness less than tiny_0 .
      DO N=1,NSHYD
         IF ( DZ(N)<tiny_0 ) THEN
             L_DZ_ZERO=.true.
         ENDIF
      ENDDO !NSHYD
      IF (L_DZ_ZERO) THEN
         errorstatus = 42
         message = "Soil layer thickness is too small."
         CALL ereport(RoutineName, errorstatus, message)
      ENDIF
      
      DO N=1,NSHYD

        DO I=1,NPNTS
!-----------------------------------------------------------------------
! Calculate TMAX, the temperature above which all soil water is
! unfrozen
!-----------------------------------------------------------------------
          SMCLSAT(I,N)=RHO_WATER*DZ(N)*V_SAT(I)
          TSL(I,N)=TSOIL(I,N)-ZERODEGC
          IF ((V_SAT(I) > tiny_0) .and. (SMCL(I,N) > small_value)) THEN
                 WORK1(I,N)=(SMCL(I,N)/SMCLSAT(I,N))**(B(I))
                 IF ( WORK1(I,N) > small_value ) THEN
                      TMAX(I)=-SATHH(I)/(DPSIDT*WORK1(I,N))
                 ELSE
                      TMAX(I)=-273.15
                 ENDIF  
          ELSE
                 TMAX(I)=-273.15
          ENDIF

!--------------------------------------------------------------------
! Diagnose unfrozen and frozen water contents
!--------------------------------------------------------------------
          IF (TSL(I,N) >= TMAX(I)) THEN
            SMCLU(I,N)=SMCL(I,N)
            SMCLF(I,N)=0.0
          ELSE
!-----------------------------------------------------------------
! For ice points (V_SAT=0) set SMCLU=0.0 and SMCLF=0.0
!-----------------------------------------------------------------
            IF (V_SAT(I) == 0.0) THEN
              SMCLU(I,N)=0.0
              SMCLF(I,N)=0.0
            ELSE
              SMCLU(I,N)=SMCLSAT(I,N)                                   &
     &                    *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/B(I))
              SMCLF(I,N)=SMCL(I,N)-SMCLU(I,N)
            ENDIF
          ENDIF
          IF (SMCLSAT(I,N) >  0.0) THEN
            STHF(I,N)=SMCLF(I,N)/SMCLSAT(I,N)
            STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)
          ELSE
            STHF(I,N)=0.0
            STHU(I,N)=0.0
          ENDIF
        ENDDO

      ENDDO

      IF (lhook) CALL dr_hook('FREEZE_SOIL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FREEZE_SOIL
