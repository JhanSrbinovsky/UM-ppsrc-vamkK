! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Calculates the fraction of each species in each
!   dissolved state.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE UKCA_FDISS(n_points, qcl_min, t, p, qcl, fdiss)


      USE ASAD_MOD,            ONLY: ddhr, dhr, kd298, k298, jpeq
      USE UKCA_CONSTANTS,      ONLY: avogadro, rmol, H_plus, m_air
      USE ukca_option_mod,     ONLY: jpdw
      USE parkind1,            ONLY: jprb, jpim
      USE yomhook,             ONLY: lhook, dr_hook
      USE Control_Max_Sizes
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n_points                 ! No of points

      REAL, INTENT(IN) :: qcl_min                     ! do calcs when qcl > qcl_min
      REAL, INTENT(IN) :: t(n_points)                 ! Temperature (K)
      REAL, INTENT(IN) :: p(n_points)                 ! Pressure (Pa)
      REAL, INTENT(IN) :: qcl(n_points)               ! Cloud liquid water (kg/kg)

      REAL, INTENT(OUT) :: fdiss(n_points, jpdw, jpeq+1)

! Local variables
      INTEGER :: ns                   ! Loop variable

      REAL, PARAMETER :: p0      = 1.0135E5   ! P0 in Pa
      REAL, PARAMETER :: inv_298 = 1.0/298.0  ! 1/298

      REAL :: tmp1(n_points)                  ! Temporary variable
      REAL :: tmp2(n_points)                  ! Temporary variable
      REAL :: Henry(n_points)                 ! Henry's Law consts
      REAL :: Henry_eff(n_points)             ! Effective  ---"---
      REAL :: qcla(n_points)                  ! qcl in kg/m^3
      LOGICAL :: todo(n_points)               ! defines qcl > qcl_min

! Aqueous equilibrium constants:
      REAL :: eq_con(n_points,jpeq)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_FDISS',zhook_in,zhook_handle)


      fdiss(:,:,:) = 0.0
      todo(:) = .false.
      WHERE (qcl(:) > qcl_min)
        todo(:) = .true.
        qcla(:) = qcl(:)*M_air*p(:)/(Rmol*t(:))
        tmp1(:) = (1.0/t(:)) - inv_298
        tmp2(:) = Rmol*qcla(:)*t(:)
      ENDWHERE

      DO ns=1,jpdw
        IF (kd298(ns,2) > 1E-12) THEN      ! 2nd dissociation constant
          WHERE (todo(:))
            Henry(:) = k298(ns) * EXP(dhr(ns)*tmp1(:))
            eq_con(:,1) = kd298(ns,1) * EXP(ddhr(ns,1)*tmp1(:))
            eq_con(:,2) = kd298(ns,2) * EXP(ddhr(ns,2)*tmp1(:))
            Henry_eff(:) = Henry(:)*(1.0+eq_con(:,1)/H_plus)
            fdiss(:,ns,1) = 1.0/(1.0 + (P0/(tmp2(:)*Henry(:))))
            fdiss(:,ns,2) = 1.0/(1.0 + (P0/(tmp2(:)*Henry_eff(:))))
            fdiss(:,ns,2) = fdiss(:,ns,2) - fdiss(:,ns,1)
            Henry_eff(:) = Henry(:) * (1.0 + eq_con(:,1)/               &
                     H_plus + eq_con(:,1)*eq_con(:,2)/H_plus**2)
            fdiss(:,ns,3) = 1.0/(1.0 + (P0/(tmp2(:)*Henry_eff(:))))
            fdiss(:,ns,3) = fdiss(:,ns,3)- fdiss(:,ns,2) - fdiss(:,ns,1)
          ENDWHERE
        ELSE IF (kd298(ns,1) > 1E-12) THEN  ! one dissociation
          WHERE (todo(:))
            Henry(:) = k298(ns) * EXP(dhr(ns)*tmp1(:))
            eq_con(:,1) = kd298(ns,1) * EXP(ddhr(ns,1)*tmp1(:))
            Henry_eff(:)=Henry(:)*(1.0+eq_con(:,1)/H_plus)
            tmp2(:) = Rmol*qcla(:)*t(:)
            fdiss(:,ns,1) = 1.0/(1.0 + (p0/(tmp2(:)*Henry(:))))
            fdiss(:,ns,2) = 1.0/(1.0 + (p0/(tmp2(:)*Henry_eff(:))))
            fdiss(:,ns,2) = fdiss(:,ns,2) - fdiss(:,ns,1)
          ENDWHERE
        ELSE                                  ! No dissociation
          IF (k298(ns) > 1E-12) THEN          ! avoid glitch with BrONO2
            WHERE (todo(:))
              Henry(:) = k298(ns) * EXP(dhr(ns)*tmp1(:))
              fdiss(:,ns,1) = 1.0/(1.0 + (p0/(tmp2(:)*Henry(:))))
            ENDWHERE
          END IF
        END IF
      END DO

      IF (lhook) CALL dr_hook('UKCA_FDISS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_FDISS
