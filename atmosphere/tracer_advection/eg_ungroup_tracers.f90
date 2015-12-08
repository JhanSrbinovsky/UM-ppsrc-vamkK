! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.  
! For further details please refer to the file COPYRIGHT.txt  
! which you should have received as part of this distribution.  
! *****************************COPYRIGHT*******************************
!  
!  
      MODULE eg_ungroup_tracers_mod     
      IMPLICIT NONE 
 
! Description: 
!      This routine ungroup/unpack a superarray into individual
!      tracers array -- this is the reverse operation of the
!      the routine "eg_group_tracers".
!  
! Method: ENDGame formulation version 3.02
! Code Owner: See Unified Model Code Owners HTML page 
! This file belongs in section: Tracer Advection 
!  
! Code description:
! Language: Fortran 95.  
! This code is written to UMDP3 standards.  

      CONTAINS
      SUBROUTINE eg_ungroup_tracers(                                  &
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
      
      USE atm_fields_bounds_mod
      USE parkind1, ONLY: jpim, jprb       !DrHook
      USE yomhook,  ONLY: lhook, dr_hook   !DrHook
       
      IMPLICIT NONE
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      
      INTEGER, INTENT(IN) :: super_array_size
      INTEGER, INTENT(IN) :: tr_ukca   
       
      REAL, INTENT(OUT)   ::  CO2      (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  MURK     (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  SOOT_NEW (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(OUT)   ::  SOOT_AGD (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  SOOT_CLD (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  SO2      (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  SO4_AITKEN                              & 
                                       (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  SO4_ACCU (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  SO4_DISS (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  NH3      (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  DMS      (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  DUST_DIV1(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  DUST_DIV2(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  DUST_DIV3(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  DUST_DIV4(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  DUST_DIV5(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  DUST_DIV6(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  BMASS_NEW(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  BMASS_AGD(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  BMASS_CLD(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  OCFF_NEW (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  OCFF_AGD (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  OCFF_CLD (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  NITR_ACC (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  NITR_DISS(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  OZONE_TRACER                            & 
                                       (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(OUT)   ::  tracer_ukca                             &
                                       (tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end,&
                                        tr_ukca )
      
      REAL, INTENT(IN)  ::  super_array                               &
                             (tdims_s%i_start:tdims_s%i_end,          &
                              tdims_s%j_start:tdims_s%j_end,          &
                              tdims_s%k_start:tdims_s%k_end,          &
                              super_array_size)  

      LOGICAL, INTENT(IN) ::                                           & 
                 L_CO2_interactive,L_Murk_advect,L_soot,L_SULPC_SO2,   & 
                 L_sulpc_nh3,L_sulpc_dms,L_biomass,L_dust, L_ocff,     &
                 L_USE_CARIOLLE, L_nitrate, L_twobin_dust
                 
      INTEGER ::  i,j,k,l,tr_count,array_size_count 

      IF (lhook) CALL dr_hook('EG_UNGROUP_TRACERS',zhook_in,zhook_handle)
      
      array_size_count=0

! ----------------------------------------------------------------------
! Section 1  carbon cycle.
! ----------------------------------------------------------------------

      IF (L_CO2_interactive) THEN
        array_size_count = array_size_count + 1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = tdims_s%j_start, tdims_s%j_end
            DO i = tdims_s%i_start, tdims_s%i_end
              CO2(i,j,k) = super_array(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF   ! L_CO2_interactive

! ----------------------------------------------------------------------
! Section 2  Soot cycle.
! ----------------------------------------------------------------------
      IF (L_soot) THEN
        array_size_count = array_size_count + 1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = tdims_s%j_start, tdims_s%j_end
            DO i = tdims_s%i_start, tdims_s%i_end
              soot_new(i,j,k) = super_array(i,j,k,array_size_count)
              soot_agd(i,j,k) = super_array(i,j,k,array_size_count+1)
              soot_cld(i,j,k) = super_array(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF    ! L_soot

! ----------------------------------------------------------------------
! Section 3  Biomass aerosol.
! ----------------------------------------------------------------------
      IF (l_Biomass) THEN
        array_size_count = array_size_count +1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = tdims_s%j_start, tdims_s%j_end
            DO i = tdims_s%i_start, tdims_s%i_end
              bmass_new(i,j,k) = super_array(i,j,k,array_size_count) 
              bmass_agd(i,j,k) = super_array(i,j,k,array_size_count+1)
              bmass_cld(i,j,k) = super_array(i,j,k,array_size_count+2) 
            END DO
          END DO
        END DO
        array_size_count=array_size_count + 2
      endif  ! l_biomass

! ----------------------------------------------------------------------
! Section 4  sulphur cycle.
! ----------------------------------------------------------------------
      IF (L_SULPC_SO2) THEN 
        array_size_count = array_size_count + 1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = tdims_s%j_start, tdims_s%j_end
            DO i = tdims_s%i_start, tdims_s%i_end
            so4_aitken(i,j,k) = super_array(i,j,k,array_size_count)  
              so4_accu(i,j,k) = super_array(i,j,k,array_size_count+1)
              so4_diss(i,j,k) = super_array(i,j,k,array_size_count+2)  
                   so2(i,j,k) = super_array(i,j,k,array_size_count+3) 
            END DO
          END DO
        END DO
        array_size_count=array_size_count +3

        IF(L_sulpc_nh3) THEN
          array_size_count = array_size_count + 1
          DO k = tdims_s%k_start, tdims_s%k_end
            DO j = tdims_s%j_start, tdims_s%j_end
              DO i = tdims_s%i_start, tdims_s%i_end
                nh3(i,j,k) = super_array(i,j,k,array_size_count)  
              END DO
            END DO
          END DO
        END IF  ! L_sulpc_nh3

        IF(L_sulpc_dms) THEN
          array_size_count = array_size_count + 1
          DO k = tdims_s%k_start, tdims_s%k_end
            DO j = tdims_s%j_start, tdims_s%j_end
              DO i = tdims_s%i_start, tdims_s%i_end
                dms(i,j,k) = super_array(i,j,k,array_size_count)  
              END DO
            END DO
          END DO
        END IF  ! L_sulpc_dms
      END IF  ! L_sulpc_SO2

! ----------------------------------------------------------------------
! Section 5  Mineral dust.
! ----------------------------------------------------------------------
      IF (L_dust) THEN
        array_size_count=array_size_count + 1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = tdims_s%j_start, tdims_s%j_end
            DO i = tdims_s%i_start, tdims_s%i_end
              dust_div1(i,j,k) = super_array(i,j,k,array_size_count  ) 
              dust_div2(i,j,k) = super_array(i,j,k,array_size_count+1) 
             END DO
          END DO
        END DO
        array_size_count = array_size_count + 1 
        IF(.NOT.L_twobin_dust) THEN
            array_size_count = array_size_count + 1
            DO k = tdims_s%k_start, tdims_s%k_end
              DO j = tdims_s%j_start, tdims_s%j_end
                 DO i = tdims_s%i_start, tdims_s%i_end
                  dust_div3(i,j,k) = super_array(i,j,k,array_size_count  ) 
                  dust_div4(i,j,k) = super_array(i,j,k,array_size_count+1)  
                  dust_div5(i,j,k) = super_array(i,j,k,array_size_count+2)  
                  dust_div6(i,j,k) = super_array(i,j,k,array_size_count+3) 
                 END DO
               END DO
             END DO
             array_size_count=array_size_count + 3
         END IF
      END IF    ! L_dust

! ----------------------------------------------------------------------
! Section 6  Fossil-fuel organic carbon aerosol
! ----------------------------------------------------------------------
      IF (L_ocff) THEN
        array_size_count=array_size_count + 1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = tdims_s%j_start, tdims_s%j_end
            DO i = tdims_s%i_start, tdims_s%i_end
              ocff_new(i,j,k) = super_array(i,j,k,array_size_count)  
              ocff_agd(i,j,k) = super_array(i,j,k,array_size_count+1)  
              ocff_cld(i,j,k) = super_array(i,j,k,array_size_count+2) 
            END DO
          END DO
        END DO
        array_size_count=array_size_count + 2
      END IF    ! L_ocff

! ----------------------------------------------------------------------
! Section 7  Cariolle ozone tracer.
! ----------------------------------------------------------------------
      IF (L_USE_CARIOLLE) THEN
        array_size_count=array_size_count + 1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = tdims_s%j_start, tdims_s%j_end
            DO i = tdims_s%i_start, tdims_s%i_end
              OZONE_TRACER(i,j,k) = super_array(i,j,k,array_size_count)  
            END DO
          END DO
        END DO
      END IF    ! L_USE_CARIOLLE

! ----------------------------------------------------------------------  
! Section 8  Ammonium nitrate aerosol. 
! ----------------------------------------------------------------------  
      IF (L_nitrate) THEN  
        array_size_count=array_size_count +1  
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = tdims_s%j_start, tdims_s%j_end
            DO i = tdims_s%i_start, tdims_s%i_end  
              nitr_acc(i,j,k) = super_array(i,j,k,array_size_count)  
              nitr_diss(i,j,k) = super_array(i,j,k,array_size_count+1) 
            END DO  
          END DO  
        END DO  
        array_size_count=array_size_count +1 
      END IF    ! L_nitrate 

! ----------------------------------------------------------------------
! Section 9  Murk cycle  
! ----------------------------------------------------------------------
      IF (L_Murk_advect) THEN
        array_size_count = array_size_count + 1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = tdims_s%j_start, tdims_s%j_end
            DO i = tdims_s%i_start, tdims_s%i_end
              murk(i,j,k) = super_array(i,j,k,array_size_count)  
            END DO
          END DO
        END DO 
      END IF  ! L_Murk_advect

! add further stracers if needed 
! ---------------------------------------------------------------------- 
! Section 10  UKCA tracers 
! ---------------------------------------------------------------------- 
      IF ( tr_ukca > 0 ) THEN 
        DO tr_count=1,tr_ukca 
          array_size_count=array_size_count +1 
          DO k = tdims_s%k_start,tdims_s%k_end 
            DO j = tdims_s%j_start,tdims_s%j_end 
              DO i = tdims_s%i_start,tdims_s%i_end 
                tracer_ukca(i,j,k,tr_count) =                            & 
                       super_array(i,j,k,array_size_count) 
              END DO 
            END DO 
          END DO 
        END DO 
 
      END IF  ! tr_ukca > 0 
 
      IF (lhook) CALL dr_hook('EG_UNGROUP_TRACERS',zhook_out,zhook_handle)
         
      RETURN
      END SUBROUTINE eg_ungroup_tracers

      END MODULE eg_ungroup_tracers_mod

