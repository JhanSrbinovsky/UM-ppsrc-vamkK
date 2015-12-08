! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Calculates aerodynamic and quasi-laminar resistances.
!  Returns Resa, Rb
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_CHEMISTRY_CTL.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!                     
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
      SUBROUTINE UKCA_AEROD(row_length, rows, nwater,                  &
        t0, p0, hflx, u_s, canht, gsf,                                 &
        zbl, z0tile, resa, rb, so4_vd)

      USE atmos_constants_mod, ONLY: vkman
      USE ASAD_MOD,            ONLY: ndepd, speci, nldepd
      USE ukca_option_mod,     ONLY: jpdd
      USE UKCA_CONSTANTS
      USE nstypes,             ONLY: ntype, npft
      USE parkind1,            ONLY: jprb, jpim
      USE yomhook,             ONLY: lhook, dr_hook
      USE Control_Max_Sizes
      USE PrintStatus_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: row_length, rows
      INTEGER, INTENT(IN) :: nwater
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: t0
                                   ! Surface temperature (K)
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: p0       
                                   ! Surface pressure (Pa)
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: hflx    
                                   ! Surface heat flux (W m-2)
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: u_s   
                                   ! Surface friction velocity (m s-1)
      REAL, DIMENSION(row_length,rows,npft), INTENT(IN) :: canht
                                   ! Canopy height of vegetation (m)
      REAL, DIMENSION(row_length,rows,ntype), INTENT(IN) :: gsf
                                   ! Surface tile fractions
      REAL, DIMENSION(row_length,rows),       INTENT(IN) :: zbl
                                   ! Boundary layer depth
      REAL, DIMENSION(row_length,rows,ntype), INTENT(INOUT) :: z0tile
                                   ! Roughness length on tiles (m)
      REAL, DIMENSION(row_length,rows,ntype), INTENT(OUT) :: resa
                                   ! Aerodynamic resistance (s m-1)
      REAL, DIMENSION(row_length,rows,jpdd), INTENT(OUT) :: rb
                                   ! Quasi-laminar resistance (s m-1)
      REAL, DIMENSION(row_length,rows), INTENT(OUT) :: so4_vd
                                   ! Sulphate aerosol dep velocity (m s-1)

      LOGICAL, SAVE :: first = .true.

      INTEGER :: i, j, k, n     ! Loop counts

      REAL :: cp = 1004.67         ! Heat capacity of air under const p (J K-1 kg-1)
      REAL :: pr = 0.72            ! Prandtl number of air
      REAL :: vair296 = 1830.0e-8  ! Dynamic viscosity of air at 296 K (kg m-1 s-1).
      REAL :: twoth = 2.0 / 3.0    ! 2/3
      REAL :: zref = 50.0          ! Reference height for dry deposition (m).
      REAL :: mair = 0.02897       ! RMM(air), kg mol-1
      REAL :: d_h2o = 2.08e-5      ! Diffusion coefficent of water in air (m2 s-1)

      REAL :: b0                   ! Temporary store
      REAL :: b1                   ! Temporary store
      REAL :: sc_pr                ! Ratio of Schmidt and Prandtl numbers

      REAL, SAVE :: reus                ! Holds -cp/(vkman*g)
      REAL, ALLOCATABLE, SAVE :: d0(:)  ! Diffusion coefficients

      REAL, DIMENSION(row_length,rows) :: l         ! Monin-Obukhov length (m)
      REAL, DIMENSION(row_length,rows) :: rho_air   ! Density of air at surface (kg m-3)
      REAL, DIMENSION(row_length,rows) :: kva       ! Kinematic velocity of air (m2 s-1)
      REAL, DIMENSION(row_length,rows) :: ustar     ! Friction velocity [checked] (m s-1)

      REAL, DIMENSION(row_length,rows,npft) :: d    ! Zero-plane displacement (m)

      REAL, DIMENSION(row_length,rows,ntype) :: z   ! Reference height (m)

      REAL, DIMENSION(row_length,rows,ntype) :: psi ! Businger function (dimensionless)
      
      REAL, PARAMETER :: r_null = 1.0e50            ! Null resist to depos (1/r_null~0)

      REAL :: bl_l                                  ! Boundary layer height divided by
                                                    ! Monin-Obukhov length

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     Assign diffusion coefficients, units m2 s-1. Set to -1 
!     unless species dry deposits. If no value found in literature, 
!     D0 calculated using: D(X) = D(H2O) * SQRT[RMM(H2O)/RMM(X)], where 
!     X is the species in question and D(H2O) = 2.08 x 10^-5 m2 s-1 
!     (Marrero & Mason, J Phys Chem Ref Dat, 1972).

!     The values of d0 will be used to flag those species 
!     that dry deposit

      IF (lhook) CALL dr_hook('UKCA_AEROD',zhook_in,zhook_handle)
      IF (first) THEN
        ALLOCATE(d0(jpdd))
        d0(:) = -1.0
        DO j = 1, ndepd
          SELECT CASE (speci(nldepd(j)))
            CASE ('O3        ','NO2       ','O3S       ','NO3       ')
              d0(j) = 1.4e-5
            ! NO = 6*NO2 (following Gian. 1998)
            CASE ('NO        ')
              d0(j) = 8.4e-5
            CASE ('HNO3      ','HONO2     ','ISON      ','B2ndry    ', &
                  'A2ndry    ','N2O5      ','HO2NO2    ','HNO4      ')
              d0(j) = 1.2e-5
            CASE ('H2O2      ','HOOH      ')
              d0(j) = 1.46e-5
            CASE ('CH3OOH    ','MeOOH     ','HCOOH     ')
              d0(j) = 1.27e-5
            CASE ('C2H5OOH   ','EtOOH     ','MeCO3H    ','MeCO2H    ', &
                  'HACET     ')
              d0(j) = 1.12e-5
            CASE ('n_C3H7OOH ','i_C3H7OOH ','n-PrOOH   ','i-PrOOH   ', &
                  'PropeOOH  ')
              d0(j) = 1.01e-5
            CASE ('MeCOCH2OOH')
              d0(j) = d_h2o * SQRT(m_h2o / m_mecoch2ooh)
            CASE ('ISOOH     ')
              d0(j) = d_h2o * SQRT(m_h2o / m_isooh)
            CASE ('HONO      ')
              d0(j) = d_h2o * SQRT(m_h2o / m_hono)
            CASE ('MACROOH   ')
              d0(j) = d_h2o * SQRT(m_h2o / m_macrooh)
            CASE ('MEKOOH    ')
              d0(j) = d_h2o * SQRT(m_h2o / m_mekooh)
            CASE ('ALKAOOH   ')
              d0(j) = d_h2o * SQRT(m_h2o / m_alkaooh)
            CASE ('AROMOOH   ')
              d0(j) = d_h2o * SQRT(m_h2o / m_aromooh)
            CASE ('BSVOC1    ')
              d0(j) = d_h2o * SQRT(m_h2o / m_bsvoc1)
            CASE ('BSVOC2    ')
              d0(j) = d_h2o * SQRT(m_h2o / m_bsvoc2)
            CASE ('ASVOC1    ')
              d0(j) = d_h2o * SQRT(m_h2o / m_asvoc1)
            CASE ('ASVOC2    ')
              d0(j) = d_h2o * SQRT(m_h2o / m_asvoc2)
            CASE ('ISOSVOC1  ')
              d0(j) = d_h2o * SQRT(m_h2o / m_isosvoc1)
            CASE ('ISOSVOC2  ')
              d0(j) = d_h2o * SQRT(m_h2o / m_isosvoc2)
            CASE ('PAN       ')
              d0(j) = 0.31e-5
            CASE ('PPAN      ')
              d0(j) = d_h2o * SQRT(m_h2o / m_ppan)
            CASE ('MPAN      ')
              d0(j) = d_h2o * SQRT(m_h2o / m_mpan)
            CASE ('ONITU     ')
              d0(j) = d_h2o * SQRT(m_h2o / m_onitu)
            CASE ('CO        ')
              d0(j) = 1.86e-5
            CASE ('CH4       ')
              d0(j) = 5.74e-5
            CASE ('NH3       ')
              d0(j) = 2.08e-5
            CASE ('H2        ')
              d0(j) = 6.7e-5
            CASE ('SO2       ')
              d0(j) = 1.2e-5
            CASE ('DMSO      ')
              d0(j) = d_h2o * SQRT(m_h2o / m_dmso)
            CASE ('SEC_ORG   ','Sec_Org   ')
              d0(j) = d_h2o * SQRT(m_h2o / m_sec_org)
            CASE ('H2SO4     ')
              d0(j) = d_h2o * SQRT(m_h2o / m_h2so4)
            CASE ('MSA       ')
              d0(j) = d_h2o * SQRT(m_h2o / m_msa)
            CASE ('MVKOOH    ')
              d0(j) = d_h2o * SQRT(m_h2o / m_mvkooh)
            CASE ('s-BuOOH   ')
              d0(j) = d_h2o * SQRT(m_h2o / m_buooh)
            CASE ('ORGNIT    ')
              d0(j) = 1.0e-5
            CASE ('BSOA      ')
              d0(j) = d_h2o * SQRT(m_h2o / m_bsoa)
            CASE ('ASOA      ')
              d0(j) = d_h2o * SQRT(m_h2o / m_asoa)
            CASE ('ISOSOA    ')
              d0(j) = d_h2o * SQRT(m_h2o / m_isosoa)
            ! Species added for sonsitency with 2D scheme
            ! Use simple approximation a la ExTC
            CASE ('HCHO')
              d0(j) = d_h2o * SQRT(m_h2o / m_hcho) 
            CASE ('MeCHO') 
              d0(j) = d_h2o * SQRT(m_h2o / m_mecho)
            CASE ('EtCHO')
              d0(j) = d_h2o * SQRT(m_h2o / m_etcho)
            CASE ('MACR')
              d0(j) = d_h2o * SQRT(m_h2o / m_macr)
            CASE ('NALD') 
              d0(j) = d_h2o * SQRT(m_h2o / m_nald)
            CASE ('MGLY')
              d0(j) = d_h2o * SQRT(m_h2o / m_mgly)
            CASE ('MeOH')
              d0(j) = d_h2o * SQRT(m_h2o / m_meoh)
            CASE ('Monoterp')
              d0(j) = d_h2o * SQRT(m_h2o / m_monoterp)
          END SELECT
        END DO
!
        DO j = 1, ndepd
          IF (d0(j) < 0.0 .AND. printstatus > Prstatus_oper) THEN
            WRITE(6,'(a31,a10,a20)') 'Warning: No dry deposition for ', &
              speci(nldepd(j)), ' will be calculated.'
          END IF
        END DO

!       Set up other constants
        reus = -cp / (vkman * g)

        first = .FALSE.

      END IF

!     Calculate air density, rho-air. Set ustar to minimum value if
!     undefined

      DO k = 1, rows
        DO i = 1, row_length
          rho_air(i,k) = mair * p0(i,k) / (rmol * t0(i,k))
          ustar(i,k) = u_s(i,k)
          IF (ustar(i,k) <= 0.0) ustar(i,k) = 1.0e-2
        END DO
      END DO

!     Set zero-plane displacement to 0.7 of canopy height 
!     (Smith et al., 2000).

      DO n = 1, npft
        DO k = 1, rows
          DO i = 1, row_length
            IF (canht(i,k,n) > 0.0) THEN
              d(i,k,n) = canht(i,k,n) * 0.70
            ELSE
              d(i,k,n) = 0.0
            END IF
          END DO
        END DO
      END DO

!     Ensure undefined parts of Z0 are set to 0.001 to avoid any
!     floating-point exceptions. Initialise z to reference height.

      DO n = 1, ntype
        DO k = 1, rows
          DO i = 1, row_length
            IF (z0tile(i,k,n) <= 0.0) z0tile(i,k,n) = 0.001
            z(i,k,n) = zref
          END DO
        END DO
      END DO

!     Set Z to height above surface minus zero plane displacement
!     for vegetated tiles.

      DO n = 1, npft
        DO k = 1, rows
          DO i = 1, row_length
            z(i,k,n) = zref - d(i,k,n)
          END DO
        END DO
      END DO

!     Calculate kinematic viscosity of air, 
!     KVA; = dynamic viscosity / density
!     Formula from Kaye and Laby, 1952, p.43.

      DO k = 1, rows
        DO i = 1, row_length
          kva(i,k) = (vair296 - 4.83e-8 * (296.0-t0(i,k))) /           &
                      rho_air(i,k)

!         Calculate roughness length over oceans using Charnock 
!         formula in form given by Hummelshoj et al., 1992.

          IF (gsf(i,k,nwater) > 0.0) THEN
            z0tile(i,k,nwater) = (kva(i,k) / (9.1 * ustar(i,k))) +     &
              0.016 * ustar(i,k) * ustar(i,k) / g
          END IF
        END DO
      END DO

!     Calculate Monin-Obukhov length L

      DO k = 1, rows
        DO i = 1, row_length
          IF (hflx(i,k) /= 0.0) THEN
!           Stable or unstable b.l.
            l(i,k) = reus * t0(i,k) * rho_air(i,k) * (ustar(i,k)**3) / &
              hflx(i,k)
          ELSE
!           Neutral b.l.
            l(i,k) = 10000.0
          END IF
        END DO
      END DO

!     Calculate Businger functions (PSI(ZETA))

      DO k = 1, rows
        DO i = 1, row_length
          IF (l(i,k) > 0.0) THEN
!           Stable b.l.
!CDIR EXPAND=9
            DO n = 1, ntype
              psi(i,k,n) = -5.0 * (z(i,k,n) - z0tile(i,k,n)) / l(i,k)
            END DO
          ELSE
!           Unstable b.l.
!CDIR EXPAND=9
            DO n = 1, ntype
              b1 = (1.0 - 16.0 * (z(i,k,n) / l(i,k)))**0.5
              b0 = (1.0 - 16.0 * (z0tile(i,k,n) / l(i,k)))**0.5
              psi(i,k,n) = 2.0 * LOG((1.0+b1) / (1.0+b0))
            END DO
          END IF
        END DO
      END DO

!   1. Resa.  Calculate aerodynamic resistance (Resa) for all species

      DO n = 1, ntype
        DO k = 1, rows
          DO i = 1, row_length
            resa(i,k,n) = (LOG(z(i,k,n) / z0tile(i,k,n)) - psi(i,k,n)) /  &
              (vkman * ustar(i,k))
          END DO
        END DO
      END DO

!   2. Rb

!    Calculate Schmidt number divided by Prandtl number for each 
!    species. Schmidt number = Kinematic viscosity of air divided 
!    by molecular diffusivity. Calculate quasi-laminar boundary 
!    layer resistance (Hicks et al., 1987)

      DO j = 1, ndepd
      
        ! First set rb for aerosols to 1.0 as done in STOCHEM code
        IF ( (speci(nldepd(j)) == 'ORGNIT    ') .OR.                    &
             (speci(nldepd(j)) == 'ASOA      ') .OR.                    &
             (speci(nldepd(j)) == 'BSOA      ') .OR.                    &
             (speci(nldepd(j)) == 'ISOSOA    ')                         &
           ) THEN
          DO k = 1, rows
            DO i = 1, row_length
              rb(i,k,j) =  1.0
            END DO
          END DO
  
        ELSE IF (d0(j) > 0.0) THEN
          DO k = 1, rows
            DO i = 1, row_length
              sc_pr = kva(i,k) / (pr * d0(j))
              rb(i,k,j) = (sc_pr**twoth) / (vkman * ustar(i,k))
            END DO
          END DO
        ELSE
          DO k = 1, rows
            DO i = 1, row_length
              rb(i,k,j) = 0.0
            END DO
          END DO
        END IF
      END DO

!    3. SO4_VD: Needed for ORGNIT or other species treated as aerosols
! 
!     Calculate surface deposition velocity term for sulphate particles
!     using simple parameterisation of Wesely in the form given by
!     Zhang et al. (2001), Atmos Environ, 35, 549-560.
!     These values can be used for different aerosol types  
!     and ORGNIT is considered to be an aerosol.
!    
      so4_vd = 1.0 / r_null
      
      DO j = 1, ndepd
        IF ( (speci(nldepd(j)) == 'ORGNIT    ') .OR.                    &
             (speci(nldepd(j)) == 'ASOA      ') .OR.                    &
             (speci(nldepd(j)) == 'BSOA      ') .OR.                    &
             (speci(nldepd(j)) == 'ISOSOA    ')                         &
           ) THEN

          DO k = 1, rows
            DO i = 1, row_length
              so4_vd(i,k) = 0.0
              bl_l = zbl(i,k) / l(i,k)
              IF (l(i,k) >= 0.0) THEN
                so4_vd(i,k) = 0.002 * ustar(i,k)
              ENDIF
              IF (bl_l < -30.0) THEN
                so4_vd(i,k) = 0.0009 * ustar(i,k) * ((-bl_l)**twoth)
              ENDIF
              IF (bl_l >= -30.0 .AND. l(i,k) < 0.0) THEN
               so4_vd(i,k) = 0.002*ustar(i,k)*(1.0 +                    &
                                (-300.0/l(i,k))**twoth)
              ENDIF
            END DO
          END DO
  
          EXIT
        END IF
      END DO
  
      IF (lhook) CALL dr_hook('UKCA_AEROD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_AEROD
