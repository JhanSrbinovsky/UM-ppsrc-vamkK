! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module containing subroutine to calculate the pressure of
!          the 2.0pvu surface and the pressure of the 380K surface.
!          Routine combines them to calculate the pressure of the
!          tropopause and set L_troposphere to .true. for those 
!          gridboxes in the troposphere.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      MODULE UKCA_TROPOPAUSE

      USE UKCA_CONSTANTS
      USE PrintStatus_mod 
      USE ereport_mod, ONLY: ereport 
      USE parkind1,    ONLY: jprb, jpim
      USE yomhook,     ONLY: lhook, dr_hook
      IMPLICIT NONE
      SAVE
      PRIVATE

      REAL, PARAMETER :: tpv = 2.0e-6        ! tropopause PV (pvu)
      REAL, PARAMETER :: tpt = 380.0         ! tropopause theta (K)

      INTEGER,ALLOCATABLE,DIMENSION(:,:),PUBLIC  ::                    &
                                         tropopause_level 

!     Pressures of the theta, pv, and combined tropopause

      REAL,ALLOCATABLE,DIMENSION(:,:),PUBLIC  :: p_tropopause
      REAL,ALLOCATABLE,DIMENSION(:,:),PUBLIC  :: theta_trop
      REAL,ALLOCATABLE,DIMENSION(:,:),PUBLIC  :: pv_trop

!     Logical set to true for gridpoints within the troposphere

      LOGICAL, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: L_troposphere

      CHARACTER(LEN=72) :: cmessage           ! Error message
      INTEGER           :: ierr               !   "   code

      PUBLIC UKCA_CALC_TROPOPAUSE

      CONTAINS

      SUBROUTINE UKCA_CALC_TROPOPAUSE(                                 &
        row_length, rows, model_levels,                                &
        sin_theta_latitude, theta, pv, pr_boundaries, pr_levels,       &
        p_tropopause, tropopause_level)

!      Description:
!       Subroutine to calculate p_tropopause. This is a weighted
!       average of 2 tropopause definitions. In the extratropics
!       (>= 28 deg), p_tropopause is the pressure (in Pa) of the
!       2.0 PVU surface. In the tropics (<= 13.0 deg), it is the
!       pressure of the 380K isentropic surface. Between the two
!       latitudes, the weighting function is
!       W = A * sech (lat) + B * lat**2 + C and then
!       p_tropopause = W *(PV tropopause) + (1-W) *(380K tropopause)
!
!       This function is from Hoerling, Monthly Weather Review,
!       Vol 121 162-172. "A Global Analysis of STE during Northern
!       Winter."
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: row_length
      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: model_levels

      INTEGER, INTENT(OUT), DIMENSION(row_length,rows) ::              &
                            tropopause_level
      REAL, PARAMETER :: lat_etropics = 28.0    ! define extratropics
      REAL, PARAMETER :: lat_tropics  = 13.0    ! define tropics
      REAL, PARAMETER :: tpv          = 2.0e-6  ! tropopause PV
      REAL, PARAMETER :: tpt          = 380.0   ! tropopause theta
      REAL, PARAMETER :: fixed_pres   = 40000.0 ! default tropopause

      REAL, INTENT(IN):: sin_theta_latitude(row_length,rows) ! sine(latitude)
      REAL, INTENT(IN):: theta(row_length,rows,model_levels) ! theta
      REAL, INTENT(IN):: pv(row_length,rows,model_levels)    ! pv on model levs
      REAL, INTENT(IN):: pr_boundaries(row_length,rows,0:model_levels) 
                                ! pressure at layer boundaries
      REAL, INTENT(IN):: pr_levels(row_length,rows,1:model_levels) 
                                ! pressure at theta levels

      REAL, INTENT(OUT), DIMENSION(row_length,rows) :: p_tropopause

!     Local variables

      INTEGER                          :: max_trop_level 
      INTEGER                          :: i,j,l           ! Loop counters
      INTEGER                          :: jll,jlu         ! Level indices

      REAL                             :: thalph

      REAL, DIMENSION(rows)            :: wt              ! weighting fn
      REAL, DIMENSION(model_levels)    :: wth             ! theta profile
      REAL, DIMENSION(model_levels)    :: wpl             ! pressure profile
      REAL, DIMENSION(model_levels)    :: wpv             ! PV profile
      REAL, DIMENSION(row_length,rows) :: theta_latitude  ! lat in radians
      REAL, DIMENSION(row_length,rows) :: latitude        ! lat in degrees

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     Calculate weighting function from Hoerling, 1993 paper
!     W = a sech (lat) + B (lat)**2 + C

      IF (lhook) CALL dr_hook('UKCA_TROPOPAUSE:UKCA_CALC_TROPOPAUSE',zhook_in,zhook_handle)
      theta_latitude = ASIN(sin_theta_latitude)         ! latitude in rads
      latitude       = theta_latitude*Recip_Pi_Over_180 ! latitude in deg

      DO i = 1,rows
        IF (ABS(latitude(1,i)) >= lat_etropics) THEN
          wt(i) = 1.0                                   ! extratropcs
        ELSE IF (ABS(latitude(1,i)) <= lat_tropics) THEN
          wt(i) = 0.0                                   ! tropics
        ELSE
          wt(i) = 5560.74*2                                            &
                /(EXP(latitude(1,i))+EXP(-1.0*latitude(1,i)))          &
                + 1.67e-3*(latitude(1,i))**2 - 0.307    ! sub-tropics
        END IF
      END DO

!     Calculate theta and pv tropopauses

      max_trop_level = model_levels-2
      DO j = 1,rows
        DO i = 1,row_length
          DO l=1,model_levels
            wth(l) = theta(i,j,l)         ! theta profile
            wpl(l) = pr_levels(i,j,l)     ! pressure profile
            wpv(l) = pv(i,j,l)            ! PV profile
          ENDDO

!         Find theta levels which straddle tpt

          DO l = max_trop_level,2,-1
            IF (wth(l) <= tpt .AND. wth(l+1) >= tpt) THEN
              jll = l
              jlu = l+1
              EXIT
            ENDIF
          END DO

!         Calculate pressure of theta tropopause

          thalph = (tpt-wth(jll))/(wth(jlu)-wth(jll))
          theta_trop(i,j) = (1.0-thalph)*log(wpl(jll))                 &
                          + thalph*log(wpl(jlu))
          IF (wpl(jll) < 0.0 .OR. wpl(jlu) < 0.0 .OR. thalph < 0.0     &
            .OR. thalph > 1.0) THEN
            IF (PrintStatus >= PrStatus_Diag) THEN
              DO l=1,model_levels
                WRITE(6,*) i,j,l,theta(i,j,l),                         &
                           pr_levels(i,j,l), pv(i,j,l)
              END DO
              WRITE(6,*) jll,jlu,thalph,theta_trop(i,j)
            END IF
          END IF
          theta_trop(i,j) = EXP(theta_trop(i,j))

!         Find theta levels which straddle tpv

          DO l = max_trop_level,2,-1
            IF (ABS(wpv(l)) <= tpv .AND. ABS(wpv(l+1)) >= tpv) THEN
              jll = l
              jlu = l+1
              EXIT
            ENDIF
          END DO

!         Calculate pressure of pv tropopause

          thalph = (tpv-ABS(wpv(jll)))/(ABS(wpv(jlu))-ABS(wpv(jll)))
          pv_trop(i,j) = (1.0-thalph)*log(wpl(jll))                    &
                       + thalph*log(wpl(jlu))
          IF (wpl(jll) < 0.0 .OR. wpl(jlu) < 0.0 .OR. thalph < 0.0     & 
            .OR. thalph > 1.0) THEN
            IF (PrintStatus >= PrStatus_Diag) THEN
              WRITE(6,'(A20)') 'UKCA_TROPOPAUSE:'
              DO l=1,model_levels
                WRITE(6,'(3I6,3E12.3)') i,j,l,theta(i,j,l),            &
                           pr_levels(i,j,l), pv(i,j,l)
              END DO
            END IF
            cmessage = ' Difficulty diagnosing pv tropopause, '//      &
                       'Reverting to default tropopause pressure'
            ierr = -1
            CALL ereport('UKCA_TROPOPAUSE::UKCA_CALC_TROPOPAUSE',      &
                          ierr,cmessage)
           
            pv_trop(i,j) = log(fixed_pres) 
          END IF

          pv_trop(i,j) = EXP(pv_trop(i,j))

        END DO
      END DO

!     Calculate combined tropopause based on weighting function

      L_troposphere = .false.
      DO j = 1,rows
        DO i = 1,row_length
          p_tropopause(i,j) = (wt(j)*pv_trop(i,j))                     &
                            + ((1.0-wt(j))*theta_trop(i,j))

!         Check for whole gridboxes below tropopause and find
!         the model level which contains the tropopause

          LOOP_trop: DO l=1,model_levels
            IF (pr_boundaries(i,j,l) >= p_tropopause(i,j)) THEN
              L_troposphere(i,j,l)  = .true.
            ELSE
              tropopause_level(i,j) = l
              EXIT LOOP_trop
            ENDIF
          END DO LOOP_trop

        END DO
      END DO

      IF (lhook) CALL dr_hook('UKCA_TROPOPAUSE:UKCA_CALC_TROPOPAUSE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALC_TROPOPAUSE

      END MODULE UKCA_TROPOPAUSE
