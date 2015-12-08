! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Convert aerosol climatology surface areas into aerosol total mass 
!    mixing ratios [kg/kg]. Two options to retrieve MMRs are
!    currently implemented but switches are not placed 
!    in namelist yet.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
      SUBROUTINE UKCA_STRAT_AERO_CLIM_MMR(row_length,rows,model_levels, &
                 r_rho_levels, rho_r2, so4_sa, so4_reff, so4_sa_mmr) 
      
      USE UKCA_CONSTANTS,  ONLY: m_s, m_h2so4
      USE yomhook,         ONLY: lhook, dr_hook
      USE parkind1,        ONLY: jprb, jpim
      IMPLICIT NONE
      

      INTEGER, INTENT(IN) :: row_length                      ! Model dimensions
      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: model_levels

! local vertical co-ordinate information
      REAL, INTENT(IN)    :: r_rho_levels(row_length,rows,model_levels)
! Density*R*R
      REAL, INTENT(IN)    :: rho_r2(row_length,rows, model_levels)
! Sulphate climatology surface area density
      REAL, INTENT(IN)    :: so4_sa(row_length, rows, model_levels)
! Sulphate climatology effective radius
      REAL, INTENT(IN)    :: so4_reff(row_length, rows, model_levels)      

! Sulphate climatology mass mixing ratio
      REAL, INTENT(OUT)   :: so4_sa_mmr(row_length, rows, model_levels)

! Local variables
      REAL              :: vdclim(row_length, rows, model_levels) ! volume density
      REAL              :: rho_air(row_length,rows, model_levels) ! air density
      
! Declare conversion factor(s) SAD -> MMR (R. Hommel 03/2009):
      REAL, PARAMETER   :: revert  = 1.0E+8  ! cm2 cm-3 to micron2 cm-3
      REAL, PARAMETER   :: fac     = 1.0E-12 ! micron3 cm-3 to SI
      REAL, PARAMETER   :: rho_aero= 1700.  ! sulphuric acid density [kg m-3]
      REAL, PARAMETER   :: weight  = 0.75   ! sulphuric acid weight fraction [0-1]


      INTEGER                       :: i, j, k        ! loop counters
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_STRAT_AERO_CLIM_MMR',zhook_in,zhook_handle)

! Calculate layer thickness for rho layers and calculate rho
      DO k=1,model_levels
        DO j = 1, rows
          DO i=1,row_length
            rho_air(i,j,k)  = rho_r2(i,j,k)/ (r_rho_levels(i,j,k)*      &
                              r_rho_levels(i,j,k))
          END DO
        END DO
      END DO

      
!-- SAD to volume density (VD) conversion after 
!   Grainger et al., JGR, 100, 16507, 1995, Eq 2 and
!   subsequent VD to MMR conversion, MMR in total sulphur [kg(S)/kg] 
!   by R. Hommel 03/2009
!
! Grainger et al, JGR, 1995, Eq 2: SAD = 8.406 VD^(0.751)
! Massie   et al, JGR, 1998, Eq 1: SAD = 7.26  VD^(0.86)

! --> VD  = (SAD/8.406)^(1/0.751) ... Grainger
!
!     MMR = VD*fac*rho_aero*w%*M_S / (M_H2SO4*rho_air), with
!
!     fac = e * 1.e-12 ... convert micron3/cm3 to SI units
!
!     rho_aero = 1700.   sulphuric acid density [kg m-3]
!     w%       = 0.75    sulphuric acid weight %   [0-1]
!     M_S      = 32.064  molar weight of S     [g mol-1]
!     M_H2SO4  = 98.0734 molar weight of H2SO4 [g mol-1]
!
!      REAL, PARAMETER :: sad2mmr_fac= 1.e-12*1700.*0.75*32.064/98.0734
!      so4_sa_mmr = so4_sa*sad2mmr


!-- another option is to estimate an effective radius and
!   convert SADs to VDs via R_eff = 3*VD/SAD
!   --> VD = R_eff*SAD/3
!   Note that estimating an so4_Reff which is "typical" is not trivial
!   since the parameter ranges from 0.09 to ~0.24 micron in 
!   stratospheric background and is enhanced after VEI>4 
!   eruptions ranging from 0.4 to 1.0 micron or even higher 
!   [e.g. Russel et al, JGR,101, 18745, 1996], not taking the negative 
!   bias of SAGEII SAD's into account [SPARC ASAP 2006; 
!   Thomason et al., ACP, 8, 983, 2008]

      vdclim(:,:,:) = so4_reff(:,:,:)*so4_sa(:,:,:)*revert / 3.
      so4_sa_mmr(:,:,:) = vdclim(:,:,:)*fac*rho_aero*weight*m_s /       &
                           (m_h2so4*rho_air(:,:,:))
      
      IF (lhook) CALL dr_hook('UKCA_STRAT_AERO_CLIM_MMR',zhook_out,     &
                               zhook_handle)
      RETURN
      END SUBROUTINE UKCA_STRAT_AERO_CLIM_MMR
