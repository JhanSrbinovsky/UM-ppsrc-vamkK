! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
  SUBROUTINE OASIS_POINT_TRANSLIST(nice,nice_use)
!
! Set up pointers to transients involved in coupling.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Coupling
!
!===========================================================

  USE oasis_atm_data_mod

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nice     ! Dims for multi-cat/layer ice fields
  INTEGER, INTENT(IN) :: nice_use ! Dims for multi-cat/layer ice fields

  ! Local loop indices
  INTEGER :: tc
  INTEGER :: n
 
  ! Now set up pointers to make processing of incoming 
  ! and outgoing arrays easier later on. Array field_id must 
  ! strictly match namcouple field ID number and be unique.

  DO tc = 1, max_transients

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check potential outgoing fields
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (transient_out(tc)%field_id == vind_heatflux) THEN
      l_heatflux = .TRUE.
      transient_out(tc)%field => heatflux(:,:)
    END IF

    IF (transient_out(tc)%field_id == vind_bluerad) THEN
      l_bluerad = .TRUE.
      transient_out(tc)%field => blue2d(:,:)
    END IF

    IF (transient_out(tc)%field_id == vind_solar) THEN
      l_solar = .TRUE.
      transient_out(tc)%field => solar2d(:,:)
    END IF

    IF (transient_out(tc)%field_id == vind_runoff) THEN
      l_runoff = .TRUE.
      transient_out(tc)%field => riverout(:,:)
    END IF
 
    IF (transient_out(tc)%field_id == vind_w10) THEN
      l_w10 = .TRUE.
      transient_out(tc)%field => w10(:,:)
    END IF

    IF (transient_out(tc)%field_id == vind_train)       THEN
      l_train = .TRUE.
      transient_out(tc)%field => totalrain(:,:)
    END IF

    IF (transient_out(tc)%field_id == vind_tsnow)       THEN
      l_tsnow = .TRUE.
      transient_out(tc)%field => totalsnow(:,:)
    END IF

    IF (transient_out(tc)%field_id == vind_evap2d)       THEN
      l_evap2d = .TRUE.
      transient_out(tc)%field => evap2d(:,:)
    END IF

    DO n = 1 , nice
      IF (transient_out(tc)%field_id == vind_topmeltn(n)) THEN
        l_topmeltn(n) = .TRUE.
        transient_out(tc)%field => topmeltn(:,:,n)
      END IF

      IF (transient_out(tc)%field_id == vind_fcondtopn(n)) THEN
        l_fcondtopn(n) = .TRUE.
        transient_out(tc)%field => fcondtopn(:,:,n)
      END IF
    END DO  ! Over n

    DO n = 1 , nice_use
      IF (transient_out(tc)%field_id == vind_sublim(n)) THEN
        l_sublim(n) = .TRUE.
        transient_out(tc)%field => sublim(:,:,n)
      END IF
    END DO  ! Over n

    IF (transient_out(tc)%field_id == vind_taux) THEN
      l_taux = .TRUE.
      transient_out(tc)%field => taux(:,:)
    END IF

    IF (transient_out(tc)%field_id == vind_tauy) THEN
      l_tauy = .TRUE.
      transient_out(tc)%field => tauy(:,:)
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check potential incoming fields
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (transient_in(tc)%field_id == vind_ocn_sst) THEN
      l_ocn_sst = .TRUE.
      transient_in(tc)%field => ocn_sst(:,:)
    END IF

    DO n = 1 , nice
      IF (transient_in(tc)%field_id == vind_ocn_freezen(n)) THEN
        l_ocn_freezen(n) = .TRUE.
        transient_in(tc)%field => ocn_freezen(:,:,n)
      END IF

      IF (transient_in(tc)%field_id == vind_ocn_snowthickn(n)) THEN
        l_ocn_snowthickn(n) = .TRUE.
        transient_in(tc)%field => ocn_snowthickn(:,:,n)
        transient_in(tc)%ice_min = .TRUE.
      END IF
 
      IF (transient_in(tc)%field_id == vind_ocn_hicen(n)) THEN
        l_ocn_hicen(n) = .TRUE.
        transient_in(tc)%field => ocn_hicen(:,:,n)
        transient_in(tc)%ice_min = .TRUE.
      END IF

      IF (transient_in(tc)%field_id == vind_ocn_icetn(n)) THEN
        l_ocn_icetn(n) = .TRUE.
        transient_in(tc)%field => ocn_icetn(:,:,n)
        transient_in(tc)%ice_min = .TRUE.
      END IF

      IF (transient_in(tc)%field_id == vind_ocn_icekn(n)) THEN
        l_ocn_icekn(n) = .TRUE.
        transient_in(tc)%field => ocn_icekn(:,:,n)
        transient_in(tc)%ice_min = .TRUE.
      END IF
    END DO 

    IF (transient_in(tc)%field_id == vind_ocn_u) THEN    
      l_ocn_u = .TRUE.
      transient_in(tc)%field => ocn_u(:,:)
      transient_in(tc)%polar_mean = .FALSE.
    END IF

    IF (transient_in(tc)%field_id == vind_ocn_v) THEN
      l_ocn_v = .TRUE.
      transient_in(tc)%field => ocn_v(:,:)
      transient_in(tc)%polar_mean = .FALSE.
    END IF

    IF (transient_in(tc)%field_id == vind_ocn_freeze) THEN
      l_ocn_freeze = .TRUE.
      transient_in(tc)%field => ocn_freeze(:,:)
    END IF

  END DO ! Over tc

  END SUBROUTINE OASIS_POINT_TRANSLIST
