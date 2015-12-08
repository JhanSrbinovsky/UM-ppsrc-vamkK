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
      SUBROUTINE UKCA_FRACDISS(row_length, rows, model_levels,          &
                               wet_levels, temp, p_theta_levels, rh,    &
                               qcl, frac_diss,kp_nh)


      USE ASAD_MOD,            ONLY: ddhr, dhr, kd298, k298, jpeq
      USE UKCA_CONSTANTS,      ONLY: avogadro, rmol, H_plus, m_air
      USE ukca_option_mod,     ONLY: jpdw
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE Control_Max_Sizes
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: row_length     ! No of points per row
      INTEGER, INTENT(IN) :: rows           ! No of rows
      INTEGER, INTENT(IN) :: model_levels   ! No of vertical levels
      INTEGER, INTENT(IN) :: wet_levels     ! No of wet levels

      REAL, INTENT(IN) :: temp(row_length,rows,model_levels)           ! Temperature
      REAL, INTENT(IN) :: p_theta_levels(row_length,rows,model_levels) ! Pressure (Pa)
      REAL, INTENT(IN) :: rh(row_length,rows,wet_levels)  ! Relative Humidity
      REAL, INTENT(IN) :: qcl(row_length,rows,wet_levels) ! Cloud liquid water (kg/kg)

      REAL, INTENT(OUT) :: frac_diss(row_length,rows,model_levels,      &
                                     jpdw,jpeq+1)
      REAL, INTENT(OUT) :: kp_nh(row_length,rows,model_levels)

! Local variables
      INTEGER :: i                    ! Loop variable
      INTEGER :: k                    ! Loop variable
      INTEGER :: ns                   ! Loop variable

      REAL, PARAMETER :: qcl_min=1.0E-12  ! do calcs when qcl > qcl_min

      REAL :: inv_298                        ! 1/298
      REAL :: tmp1(row_length,rows,wet_levels)        ! Temporary variable
      REAL :: tmp2(row_length,rows,wet_levels)        ! Temporary variable
      REAL :: Henry(row_length,rows,wet_levels)       ! Henry's Law consts
      REAL :: Henry_eff(row_length,rows,wet_levels)   ! Effective  ---"---
      REAL :: twet(row_length,rows,wet_levels)        ! Temperature on wet levels
      REAL :: qcla(row_length,rows,wet_levels)        ! qcl in kg/m^3
      LOGICAL :: todo(row_length,rows,wet_levels)     ! defines qcl > qcl_min

! Aqueous equilibrium constants:
      REAL :: eq_con(row_length,rows,wet_levels,jpeq)
      REAL :: frac_aq(row_length,rows,wet_levels)     ! Aqueous fraction

! For ammonium nitrate dissociation
      REAL :: rhd(row_length,rows,wet_levels)        ! rhd deliquesance RH
      REAL :: kp(row_length,rows)                    ! kp
      REAL :: kp1(row_length,rows)                   ! kp1
      REAL :: p1(row_length,rows)                    ! p1
      REAL :: p2(row_length,rows)                    ! p2
      REAL :: p3(row_length,rows)                    ! p3

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_FRACDISS',zhook_in,zhook_handle)

      inv_298 = 1.0/298.0

      frac_diss(:,:,:,:,:)=0.0
      todo(:,:,:) = .false.
      DO k=1,wet_levels
        WHERE (qcl(:,:,k) > qcl_min)
          todo(:,:,k) = .true.
          twet(:,:,k) = temp(:,:,k)
          qcla(:,:,k)=qcl(:,:,k)*M_air*p_theta_levels(:,:,k)/           &
                      (Rmol*temp(:,:,k))
        ENDWHERE
      ENDDO

      DO ns=1,jpdw
        IF (kd298(ns,2) > 1E-12) THEN      ! 2nd dissociation constant
          WHERE (todo(:,:,:))
            tmp1(:,:,:) = (1.0/twet(:,:,:)) - inv_298
            Henry(:,:,:) = k298(ns) * EXP(dhr(ns)*tmp1(:,:,:))
            eq_con(:,:,:,1) = kd298(ns,1) * EXP(ddhr(ns,1)*tmp1(:,:,:))
            eq_con(:,:,:,2) = kd298(ns,2) * EXP(ddhr(ns,2)*tmp1(:,:,:))
            Henry_eff(:,:,:)=Henry(:,:,:)*(1.0+eq_con(:,:,:,1)/H_plus)
            tmp2(:,:,:) = Rmol*qcla(:,:,:)*twet(:,:,:)
            frac_diss(:,:,:,ns,1) = 1.0/(1.0 + (1.013e5/                &
                                    (tmp2(:,:,:)*Henry(:,:,:))))
            frac_diss(:,:,:,ns,2) = 1.0/(1.0 + (1.013e5/                &
                                   (tmp2(:,:,:)*Henry_eff(:,:,:))))
            frac_diss(:,:,:,ns,2) = frac_diss(:,:,:,ns,2)-              &
                                    frac_diss(:,:,:,ns,1)
            Henry_eff(:,:,:) = Henry(:,:,:) * (1.0 + eq_con(:,:,:,1)/   &
                     H_plus + eq_con(:,:,:,1)*eq_con(:,:,:,2)/H_plus**2)
            frac_diss(:,:,:,ns,3) = 1.0/(1.0 + (1.013e5/                &
                                   (tmp2(:,:,:)*Henry_eff(:,:,:))))
            frac_diss(:,:,:,ns,3)=frac_diss(:,:,:,ns,3)-                &
                     frac_diss(:,:,:,ns,2) - frac_diss(:,:,:,ns,1)
          ENDWHERE
        ELSE IF (kd298(ns,1) > 1E-12) THEN  ! one dissociation
          WHERE (todo(:,:,:))
            tmp1(:,:,:) = (1.0/twet(:,:,:)) - inv_298
            Henry(:,:,:) = k298(ns) * EXP(dhr(ns)*tmp1(:,:,:))
            eq_con(:,:,:,1) = kd298(ns,1) * EXP(ddhr(ns,1)*tmp1(:,:,:))
            Henry_eff(:,:,:)=Henry(:,:,:)*(1.0+eq_con(:,:,:,1)/H_plus)
            tmp2(:,:,:) = Rmol*qcla(:,:,:)*twet(:,:,:)
            frac_diss(:,:,:,ns,1) = 1.0/(1.0 + (1.013e5/                &
                                    (tmp2(:,:,:)*Henry(:,:,:))))
            frac_diss(:,:,:,ns,2) = 1.0/(1.0 + (1.013e5/                &
                                    (tmp2(:,:,:)*Henry_eff(:,:,:))))
            frac_diss(:,:,:,ns,2) = frac_diss(:,:,:,ns,2)-              &
                                    frac_diss(:,:,:,ns,1)
          ENDWHERE
        ELSE                                  ! No dissociation
          WHERE (todo(:,:,:))
            tmp1(:,:,:) = (1.0/twet(:,:,:)) - inv_298
            Henry(:,:,:) = k298(ns) * EXP(dhr(ns)*tmp1(:,:,:))
            tmp2(:,:,:)=qcla(:,:,:)*Henry(:,:,:)*Rmol*twet(:,:,:)
            frac_diss(:,:,:,ns,1) = 1.0/(1.0 + (1.013e5/tmp2(:,:,:)))
          ENDWHERE
        ENDIF
      ENDDO

! Calculate NH4NO3 dissociation constant
      DO k=1,wet_levels
        rhd(:,:,k) = exp((618.3/temp(:,:,k))-2.551)
      ENDDO

      kp_nh(:,:,:)=0.0
      todo(:,:,:)=.false.
      todo = (rh < rhd)
      DO k=1,wet_levels
          WHERE (todo(:,:,k))
            kp(:,:)=(temp(:,:,k)**(-6.025))*EXP(118.87-24084./          &
                     temp(:,:,k))
          ELSEWHERE
            p1(:,:)=exp(-135.94+(8763/temp(:,:,k)))*temp(:,:,k)**19.12
            p2(:,:)=exp(-122.65+(9969/temp(:,:,k)))*temp(:,:,k)**16.22
            p3(:,:)=exp(-182.61+(13875/temp(:,:,k)))*temp(:,:,k)**24.46
            kp1(:,:)=(temp(:,:,k)**(-6.025))*EXP(118.87-24084./         &
                      temp(:,:,k))
            kp(:,:)=(p1(:,:)-p2(:,:)*(1.0-rh(:,:,k))+p3(:,:)*           &
                    (1.0-rh(:,:,k))**2)*((1.0-rh(:,:,k))**1.75)*kp1(:,:)
          ENDWHERE
          kp_nh(:,:,k)=kp(:,:)*1.0E-8/((rmol*1.0E6*temp(:,:,k)/         &
                       avogadro)**2)
      ENDDO

      IF (lhook) CALL dr_hook('UKCA_FRACDISS',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE UKCA_FRACDISS
