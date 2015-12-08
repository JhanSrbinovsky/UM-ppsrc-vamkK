! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Module to calculate pc2 increments
!
MODULE tcs_pc2


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  !
  ! Description:
  ! Module to calculate pc2 increments
  !
  ! Method:
  !   Currenly uses a fairly ad hoc treatment based on in-cloud
  !   liquid water flux and detrainment rates defined in 
  !   tcs_parameters_warm
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.1 programming standards.
  !

CONTAINS

  SUBROUTINE calc_pc2( n_xx, max_cldlev, timestep, ntml, ntpar,       &
     mf_h_cld, ql_up_cld, wql_cld, qcl, qcf, cf_liquid, cf_frozen,    &
     cld_in, sim, dqclbydt, dqcfbydt, dbcfbydt, dcflbydt, dcffbydt)

    !-------------------------------------------------------------------
    !
    ! Description:
    !   Calculate the increments to PC2 variables.
    !
    !-------------------------------------------------------------------

    USE tcs_parameters_warm,    ONLY : Aepsilon
    USE tcs_class_similarity,   ONLY : similarity
    USE tcs_class_cloud,        ONLY : cloud_input
    USE tcs_common_warm,        ONLY : scales

    IMPLICIT NONE
    !----------------------------------------------------------------
    ! Subroutine Arguments 
    !----------------------------------------------------------------
    INTEGER, INTENT(in) ::                                            &
       n_xx                                                           &
                          ! No. of congestus convection points
       , max_cldlev                                                     
                          ! Maximum number of convective cloud levels

    REAL, INTENT(in) ::                                               &
       timestep
                                ! Model timestep (s)

    INTEGER, INTENT(in) ::                                            &
       ntml(:)                                                        &
                                ! Top level of surface mixed layer defined
                                ! relative to theta,q grid
       ,ntpar(:)                   
    ! Top level of surface mixed layer defined
    ! relative to theta,q grid                     

    REAL, INTENT(in) ::                                               &
       mf_h_cld(:,:)                                                  &
                                ! mass flux in cld (rho levels)
       ,ql_up_cld(:,:)                                                &
                                ! liquid water content in the updraught
       ,wql_cld(:,:)
                                ! liquid water flux 

    REAL, INTENT(inout) ::                                            &
       qcl(:,:)                                                       &
                                ! Liq condensate mix ratio (kg/kg)
       , qcf(:,:)                                                     &
                                ! Ice condensate mix ratio (kg/kg)
       , cf_liquid(:,:)                                               &
                                ! Frozen water cloud volume             
       , cf_frozen(:,:)                                         
    ! Liq water cloud volume

    TYPE(cloud_input), INTENT(in) :: cld_in ! fields on cloud levels

    TYPE(similarity), INTENT(in) :: sim !Similarity functions

    REAL, INTENT(inout) ::                                            &
       dqclbydt(:,:)                                                  &
       ,dqcfbydt(:,:)                                                 &
       ,dbcfbydt(:,:)                                                 &
       ,dcflbydt(:,:)                                                 & 
       ,dcffbydt(:,:) 

    !-------------------------------------------------------------------
    ! Variables defined locally
    !-------------------------------------------------------------------

    REAL, ALLOCATABLE :: dmf(:,:)   ! mass flux divergence
    REAL, ALLOCATABLE :: entr(:,:)  ! entrainment rate
    REAL, ALLOCATABLE :: detr(:,:)  ! detrainment rate

    !-------------------------
    ! Loop counters
    !-------------------------
    INTEGER :: i,k,lbase

    !=================
    ! Select method 
    !=================
    INTEGER, PARAMETER :: ipc2_method=2

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !-------------------------------------------------------------------
    ! 1.0 Allocate and initialise arrays
    !-------------------------------------------------------------------

    IF (lhook) CALL dr_hook('TCS_PC2:CALC_PC2',zhook_in,zhook_handle)

    ALLOCATE(dmf(n_xx,max_cldlev))
    ALLOCATE(detr(n_xx,max_cldlev))
    ALLOCATE(entr(n_xx,max_cldlev))

    dqclbydt(:,:) = 0.0
    dqcfbydt(:,:) = 0.0
    dbcfbydt(:,:) = 0.0
    dcflbydt(:,:) = 0.0
    dcffbydt(:,:) = 0.0

    IF (ipc2_method==1)THEN

      !-----------------------------------------------------------------
      ! 2.0 calculate mass flux divergence
      !-----------------------------------------------------------------

      DO k=1,max_cldlev-1
        DO i = 1,n_xx
          dmf(i,k) = (mf_h_cld(i,k+1) - mf_h_cld(i,k))/                     &
             (cld_in%z_rho(i,k+1)-cld_in%z_rho(i,k))
        END DO
      END DO

      k=max_cldlev
      DO i = 1,n_xx
        dmf(i,k) = mf_h_cld(i,k)/cld_in%z_rho(i,k)
      END DO


      !-----------------------------------------------------------------
      ! 3.0 calculate entrainment and detrainment rates
      !-----------------------------------------------------------------

      DO k=1,max_cldlev-1
        DO i = 1,n_xx
          entr(i,k) = Aepsilon*scales%wstar_up(i)/scales%mb(i)              &
             /(scales%zcld(i) * cld_in%eta_theta(i,k))
          detr(i,k) = entr(i,k)-dmf(i,k)
        END DO
      END DO

      !-----------------------------------------------------------------
      ! 3.1 Put sensible limits on detr
      !-----------------------------------------------------------------

      DO k=1,max_cldlev-1
        DO i = 1,n_xx
          IF (k > ntpar(i) - ntml(i) + 1 &
             .OR. k < ntml(i) + 1        &
             .OR. detr(i,k) < 0.) detr(i,k) = 0.0
          IF (detr(i,k) > .5) detr(i,k) = .5
        END DO
      END DO


      !-----------------------------------------------------------------
      ! 4.0 calculate detrained water
      !-----------------------------------------------------------------

      ! working note: for now we just calculate detrained liquid - if the 
      ! scheme is being run in a cold environment this may not be 
      ! appropriate.
      ! working note: we don't detrain anything in the lowest cloud level
      DO i = 1,n_xx
        DO k=2,ntpar(i)-ntml(i)
          lbase = ntml(i)+k
          dqclbydt(i,lbase) = detr(i,k)*(ql_up_cld(i,k-1)-qcl(i,lbase))      &
             /cld_in%rho(i,k)
          dcflbydt(i,lbase) = detr(i,k)*(1.-cf_liquid(i,lbase))/cld_in%rho(i,k)
        END DO
        ! Enhanced detrainment at cloud top
        k=ntpar(i)-ntml(i)+1
        lbase = ntml(i)+k
        dqclbydt(i,lbase) = detr(i,k)*(ql_up_cld(i,k-1))/cld_in%rho(i,k)
        dcflbydt(i,lbase) = detr(i,k)*(1.-cf_liquid(i,k))/cld_in%rho(i,k)
      END DO


    ELSEIF (ipc2_method==2)THEN
      DO i = 1,n_xx
        DO k=1,ntpar(i)-ntml(i)+1
          lbase = ntml(i)+k
          detr(i,k) = sim%pc2_detr(i,k)
          dqclbydt(i,lbase) = detr(i,k)*wql_cld(i,k)
        END DO
      END DO
    END IF

    !-------------------------------------------------------------------
    ! 5.0 Put sensible limits on increments and work out change in cloud 
    !     fraction
    !-------------------------------------------------------------------

    DO i = 1,n_xx
      DO k=1,ntpar(i)-ntml(i)
        lbase = ntml(i)+k
        dqclbydt(i,lbase) = MAX(dqclbydt(i,lbase), -qcl(i,lbase)/timestep)
        ! dqclbydt(i,lbase) = min(dqclbydt(i,lbase), ql_up_cld(i,k)/timestep)
        dcflbydt(i,lbase) = dqclbydt(i,lbase)*               &
           (1.-cf_liquid(i,k))/(ql_up_cld(i,k)-qcl(i,lbase))
        dcflbydt(i,lbase) = MAX(MIN(dcflbydt(i,lbase),1.),0.)
      END DO
      k=ntpar(i)-ntml(i)+1
      lbase = ntml(i)+k
      dqclbydt(i,lbase) = MAX(dqclbydt(i,lbase), -qcl(i,lbase)/timestep)
      ! dqclbydt(i,lbase) = min(dqclbydt(i,lbase), ql_up_cld(i,k-1)/timestep)
      dcflbydt(i,lbase) = dqclbydt(i,lbase)*               &
         (1.-cf_liquid(i,k))/(ql_up_cld(i,k)-qcl(i,lbase) + 1e-10)
      dcflbydt(i,lbase) = MAX(MIN(dcflbydt(i,lbase),1.),0.)
    END DO

    !-------------------------------------------------------------------
    ! 6.0 Deallocation
    !-------------------------------------------------------------------

    DEALLOCATE(entr)
    DEALLOCATE(detr)
    DEALLOCATE(dmf)
    IF (lhook) CALL dr_hook('TCS_PC2:CALC_PC2',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE calc_pc2
    
END MODULE tcs_pc2
