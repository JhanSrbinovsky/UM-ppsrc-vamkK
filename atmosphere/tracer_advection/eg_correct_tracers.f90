! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.  
! For further details please refer to the file COPYRIGHT.txt  
! which you should have received as part of this distribution.  
! *****************************COPYRIGHT*******************************
!    
      MODULE eg_conserv_tracers_mod   
      IMPLICIT NONE 
 
! Description:  
!
!            This routine enforces the mass conservation of tracers.
!            It should work with either mixing ratios and dry density 
!            or with specific humdities and wet density. 
!            If L_conserve_tracers=.FALSE., then the routine make 
!            sure that the tracers are strictly positive only and
!            doesn't impose any mass conservation constraint.
!  
! Method: ENDGame formulation version 3.02
! Code Owner: See Unified Model Code Owners HTML page 
! This file belongs in section: Tracer Advection
!  
! Code description:
! Language: Fortran 95.  
! This code is written to UMDP3 standards.  

      CONTAINS
      SUBROUTINE eg_correct_tracers(                                  &        
                               mype, super_array_size,                &
                               super_tracer_phys1, super_tracer_phys2,&
                               rho_n, rho_np1,                        &
                               CO2, L_CO2_interactive,                &
                               murk, L_murk_advect,                   &
                               soot_new, soot_agd, soot_cld, L_soot,  &
                               bmass_new, bmass_agd, bmass_cld,       &
                               L_biomass,                             &
                               ocff_new, ocff_agd, ocff_cld, l_ocff,  &
                               DUST_DIV1,DUST_DIV2,DUST_DIV3,         &
                               DUST_DIV4,DUST_DIV5,DUST_DIV6,         &
                               L_DUST,                                &
                               so2, so4_aitken, so4_accu,             &
                               so4_diss, nh3, dms,                    &
                               L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms, &
                               nitr_acc, nitr_diss, L_nitrate,        &
                               L_USE_CARIOLLE, OZONE_TRACER,          &
                               tr_ukca, tracer_ukca,                  &
                               L_conserve_tracers,                    &
                               L_conserv_smooth_lap                   )

      USE eg_group_tracers_mod
      USE eg_ungroup_tracers_mod
      USE eg_check_conserv_mod
      USE eg_mass_conserv_mod
      USE dust_parameters_mod, ONLY: L_twobin_dust
      USE atm_fields_bounds_mod
      USE PrintStatus_mod 
      USE Field_Types
      USE parkind1, ONLY: jpim, jprb       !DrHook
      USE yomhook,  ONLY: lhook, dr_hook   !DrHook
   
      IMPLICIT NONE
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER, INTENT(IN)  :: mype 
      INTEGER, INTENT(IN)  :: super_array_size
      INTEGER, INTENT(IN)  :: tr_ukca  
      LOGICAL, INTENT(IN)  :: L_CO2_interactive  
      LOGICAL, INTENT(IN)  :: L_murk_advect, L_Soot
      LOGICAL, INTENT(IN)  :: L_biomass, L_ocff, L_DUST, L_sulpc_so2
      LOGICAL, INTENT(IN)  :: L_sulpc_nh3, l_sulpc_dms, L_USE_CARIOLLE
      LOGICAL, INTENT(IN)  :: L_nitrate, L_conserve_tracers
      LOGICAL, INTENT(IN)  :: L_conserv_smooth_lap
 
      REAL, INTENT(INOUT) ::  CO2      (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  MURK     (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  SOOT_NEW (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  SOOT_AGD (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  SOOT_CLD (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  SO2      (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  SO4_AITKEN                              & 
                                       (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  SO4_ACCU (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  SO4_DISS (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  NH3      (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  DMS      (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  DUST_DIV1(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  DUST_DIV2(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  DUST_DIV3(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  DUST_DIV4(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  DUST_DIV5(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  DUST_DIV6(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  BMASS_NEW(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  BMASS_AGD(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  BMASS_CLD(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  OCFF_NEW (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  OCFF_AGD (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  OCFF_CLD (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  NITR_ACC (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  NITR_DISS(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::  OZONE_TRACER                            & 
                                       (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) :: tracer_ukca                              &
                                       (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end,&
                                        tr_ukca)
 
      REAL, INTENT(IN) :: super_tracer_phys1                          &
                     (tdims_l%i_start:tdims_l%i_end,                  &
                      tdims_l%j_start:tdims_l%j_end,                  &
                      tdims_l%k_start:tdims_l%k_end,                  &
                      super_array_size)
       
      REAL, INTENT(IN) :: super_tracer_phys2                          &
                     (tdims%i_start:tdims%i_end,                      &
                      tdims%j_start:tdims%j_end,                      &
                      tdims%k_start:tdims%k_end,                      &
                      super_array_size)   
       
      REAL, INTENT(IN)   :: rho_n                                     &
                     (pdims_s%i_start:pdims_s%i_end,                  &
                      pdims_s%j_start:pdims_s%j_end,                  &
                      pdims_s%k_start:pdims_s%k_end)

      REAL, INTENT(IN)   :: rho_np1                                   &
                     (pdims_s%i_start:pdims_s%i_end,                  &
                      pdims_s%j_start:pdims_s%j_end,                  &
                      pdims_s%k_start:pdims_s%k_end)     
     
! Local Variables

      REAL ::    super_array (tdims_s%i_start:tdims_s%i_end,          &
                              tdims_s%j_start:tdims_s%j_end,          &
                              tdims_s%k_start:tdims_s%k_end,          &
                              super_array_size) 
    
      INTEGER :: tracers_switches(super_array_size) 
      INTEGER :: ist,iend, number_qs, total_number_tracers 
      INTEGER :: i, mpierr

      REAL,    ALLOCATABLE :: qs_n(:,:,:,:)
      REAL,    ALLOCATABLE :: qs_np1(:,:,:,:)
      REAL,    ALLOCATABLE :: qs_s(:,:,:,:) 
      REAL,    ALLOCATABLE :: qsmin(:)
      INTEGER, ALLOCATABLE :: qswitches(:)
 
      IF (lhook) CALL dr_hook('EG_CORRECT_TRACERS',zhook_in,          &
                                                   zhook_handle)

! pack/group all the tracers sources into a superarray
                
      CALL eg_group_tracers(                                          &
                          super_array_size,                           &
                          total_number_tracers,                       &
                          super_array,                                &
                          L_CO2_interactive, CO2,                     &
                          L_Murk_advect, murk,                        &
                          L_Soot, soot_new, soot_agd, soot_cld,       &
                          L_SULPC_SO2, SO2, SO4_aitken, so4_accu,     &
                                       so4_diss,                      &
                          L_sulpc_nh3, nh3,                           &
                          L_sulpc_dms, dms,                           &
                          L_dust, DUST_DIV1, DUST_DIV2, DUST_DIV3,    &
                                  DUST_DIV4, DUST_DIV5, DUST_DIV6,    &
                          L_biomass, bmass_new, bmass_agd, bmass_cld, &
                          L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                          L_nitrate, nitr_acc, nitr_diss,             &
                          L_USE_CARIOLLE, OZONE_TRACER,               &
                          tr_ukca, tracer_ukca,                       &
                          L_twobin_dust, tracers_switches             )
        
        
      IF (L_conserve_tracers ) THEN
 
        number_qs =  total_number_tracers  
             
        ALLOCATE(                                                     &
               qs_n(tdims_s%i_start:tdims_s%i_end,                    &
                    tdims_s%j_start:tdims_s%j_end,                    &
                    tdims_s%k_start:tdims_s%k_end,number_qs))

        ALLOCATE(                                                     &
             qs_np1(tdims_s%i_start:tdims_s%i_end,                    &
                    tdims_s%j_start:tdims_s%j_end,                    &
                    tdims_s%k_start:tdims_s%k_end,number_qs))

        ALLOCATE(                                                     &
               qs_s(tdims_s%i_start:tdims_s%i_end,                    &
                    tdims_s%j_start:tdims_s%j_end,                    &
                    tdims_s%k_start:tdims_s%k_end,number_qs))

         ALLOCATE(                                                    &
             qsmin(number_qs))

        ALLOCATE(                                                     &
              qswitches(number_qs))
     
  
        qswitches = tracers_switches(1:number_qs)
        qsmin     = 0.0           
        qs_np1    = super_array(:,:,:,1:number_qs)

        ist  = 1
        iend = total_number_tracers 

! put the sources into one superarray conistent with the field
! superarray
   
        qs_n(:,:,:,ist:iend) = super_tracer_phys1(                    &
                             tdims_s%i_start:tdims_s%i_end,           &
                             tdims_s%j_start:tdims_s%j_end,           &
                             tdims_s%k_start:tdims_s%k_end, ist:iend )
        
        qs_s( tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                              &
              tdims%k_start:tdims%k_end, ist:iend ) =                 &
        super_tracer_phys2(                                           &      
                 tdims%i_start:tdims%i_end,                           &
                 tdims%j_start:tdims%j_end,                           &
                 tdims%k_start:tdims%k_end, ist:iend      )
 
! Enforce the boundary condition filed(surface)=field(level 1)
! 
! These are the assumption used in computing the total mass.
! In some array they are forced earlier on and in some arrays
! they are not --- and this is why they are forced here 
! so this routine doesn't miss-diagnose the total mass of tracers
  
          qs_n(:,:,tdims_s%k_start,:) =   qs_n(:,:,tdims_s%k_start+1,:)
        qs_np1(:,:,tdims_s%k_start,:) = qs_np1(:,:,tdims_s%k_start+1,:)
          qs_s(:,:,tdims_s%k_start,:) =   qs_s(:,:,tdims_s%k_start+1,:)
 
! maybe this swap_bounds is not necessary, but it's here for safety after
! the boundary conditions are explicitly imposed -- because some tracers 
! don't have anything at level 0 higher-up and this swap-bounds to make 
! sure all the points have sensible data so mass is not missdiagnosed.    

        CALL Swap_Bounds(qs_np1,                                          &
             tdims%i_end-tdims%i_start+1,tdims%j_end-tdims%j_start+1,     &
             number_qs*(tdims%k_end-tdims%k_start+1),                     &
             tdims%i_start-tdims_s%i_start,tdims%j_start-tdims_s%j_start, &
             fld_type_p, .FALSE.                                          ) 
    
! print mass conservation error before correction if needed 

        IF (printstatus > prstatus_normal) THEN
          IF( mype == 0) THEN
            WRITE(6,fmt='(A)')                                        &
             'Error in mass conservation for tracers before correction'
          END IF                 
          CALL eg_check_mass_conservation(rho_n, rho_np1, qs_n,       &
                                          qs_np1, qs_s,               & 
                       qswitches, number_qs, mype,'eg_correct_tracers') 
        END IF        

! enfore mass conservation

        CALL eg_mass_conservation(rho_n, rho_np1, qs_n, qs_np1,       &
                                  qs_s, qsmin,qswitches, number_qs,   &
                                  L_conserv_smooth_lap                ) 

! print mass conservation error after correction if needed 

           IF (printstatus > prstatus_normal) THEN
              IF( mype == 0) THEN
                WRITE(6,fmt='(A)') &
                'Error in mass conservation for tracers after correction'
              END IF                 
              CALL eg_check_mass_conservation( rho_n, rho_np1,      &
                                               qs_n, qs_np1, qs_s,  & 
                                        qswitches, number_qs, mype, &
                                               'eg_correct_tracers' ) 
           END IF 

           super_array(:,:,:,1:number_qs) = qs_np1
      
           ! may be there is no need for this swap-bounds if the tracers
           ! are swap-bounded later on -- but it's safer to leave it
            
           CALL Swap_Bounds(super_array,                                    &
               tdims%i_end-tdims%i_start+1,tdims%j_end-tdims%j_start+1,     &
               number_qs*(tdims%k_end-tdims%k_start+1),                     &
               tdims%i_start-tdims_s%i_start,tdims%j_start-tdims_s%j_start, &
               fld_type_p, .FALSE.                                          ) 
      
        DEALLOCATE( qswitches )
        DEALLOCATE( qsmin     )
        DEALLOCATE( qs_s      )
        DEALLOCATE( qs_np1    )
        DEALLOCATE( qs_n      )
      
      ELSE     

! if L_conserve_tracers==.false. then simply make sure that the tracers
! are strictly positive

        super_array = max(super_array, 0.0 )
        
      END IF ! L_conserve_tracers

! unpack/ungroup the superarray. Copy the superarray fields into their 
! appropriate arrays (i.e.,reverse the operation in "eg_group_tracers")

      CALL eg_ungroup_tracers(                                        &
                          super_array_size,                           & 
                          super_array,                                &
                          L_CO2_interactive, CO2,                     &
                          L_Murk_advect, murk,                        &
                          L_Soot, soot_new, soot_agd, soot_cld,       &
                          L_SULPC_SO2, SO2, SO4_aitken, so4_accu,     &
                                       so4_diss,                      &
                          L_sulpc_nh3, nh3,                           &
                          L_sulpc_dms, dms,                           &
                          L_dust, DUST_DIV1, DUST_DIV2, DUST_DIV3,    &
                                  DUST_DIV4, DUST_DIV5, DUST_DIV6,    &
                          L_biomass, bmass_new, bmass_agd, bmass_cld, &
                          L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                          L_nitrate, nitr_acc, nitr_diss,             &
                          L_USE_CARIOLLE, OZONE_TRACER,               &
                          tr_ukca, tracer_ukca,                       &
                          L_twobin_dust                               )

      IF (lhook) CALL dr_hook('EG_CORRECT_TRACERS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eg_correct_tracers

      END MODULE eg_conserv_tracers_mod
