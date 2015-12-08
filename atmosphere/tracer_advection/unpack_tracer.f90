! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
      MODULE unpack_tracer_mod
      IMPLICIT NONE
      CONTAINS

      SUBROUTINE unpack_tracer(                                       &
                               super_array_size,                      &
                               super_array, supertrdims ,             &
                               co2, l_co2_interactive,                &
                               murk, l_murk_advect,                   &
                               soot_new, soot_agd, soot_cld, l_soot,  &
                               bmass_new, bmass_agd, bmass_cld,       &
                               l_biomass,                             &
                               ocff_new, ocff_agd, ocff_cld, l_ocff,  &
                               dust_div1,dust_div2,dust_div3,         &
                               dust_div4,dust_div5,dust_div6,         &
                               l_dust,                                &
                               so2, so4_aitken, so4_accu,             &
                               so4_diss, nh3, dms,                    &
                               l_sulpc_so2, l_sulpc_nh3, l_sulpc_dms, &
                               nitr_acc, nitr_diss, l_nitrate,        &
                               tracers, tr_levels, tr_vars,           &
                               tracers_ukca,tr_ukca,                  &
                               l_use_cariolle, ozone_tracer,          &
                               error_code)

! Purpose:
!          Unpacks tracer super array
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Tracer Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE dust_parameters_mod, ONLY: l_twobin_dust
      USE yomhook,             ONLY: lhook, dr_hook
      USE parkind1,            ONLY: jprb, jpim

      USE UM_ParParams
      USE atm_fields_bounds_mod
      USE Field_Types

      IMPLICIT NONE

      TYPE (array_dims)  supertrdims 

      INTEGER  tr_levels, tr_vars,  tr_ukca
! Arguments with Intent IN. ie: Input variables.

      LOGICAL, Intent(IN) ::                                          &
        l_co2_interactive,                                            &
        l_murk_advect,                                                &
        l_soot,                                                       &
        l_biomass,                                                    &
        l_ocff,                                                       &
        l_dust,                                                       &
        l_sulpc_so2, l_sulpc_nh3, l_sulpc_dms,                        &
        l_use_cariolle,                                               &
        l_nitrate
     
      REAL, INTENT(INOUT) :: co2                                      &
                        (tdims_s%i_start:tdims_s%i_end,               &
                         tdims_s%j_start:tdims_s%j_end,               &
                         tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: murk                                    &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: soot_new                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: soot_agd                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: soot_cld                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: so2                                     &  
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: so4_aitken                              &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: so4_accu                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: so4_diss                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT)  :: nh3                                     &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dms                                     &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT)  :: dust_div1                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div2                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div3                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div4                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div5                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div6                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: bmass_new                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: bmass_agd                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: bmass_cld                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT)  :: ocff_new                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: ocff_agd                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: ocff_cld                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: nitr_acc                                &
                       (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT)  :: nitr_diss                               &
                       (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

       REAL, INTENT(INOUT)  ::                                        &
                tracers(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tr_levels,tr_vars),                           &
           tracers_ukca(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end,tr_ukca),       &
! Add cariolle specific parameters for ozone tracer     
           OZONE_TRACER(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)                                          


      INTEGER                                                         &
        error_code     ! non-zero on exit if error detected.


! scalars

      INTEGER                                                         &
        i, j, k, l ,count                 

      INTEGER super_array_size, array_size_count

      REAL :: super_array (supertrdims%i_start:supertrdims%i_end,     &
                           supertrdims%j_start:supertrdims%j_end,     &
                           supertrdims%k_start:supertrdims%k_end,     &
                           super_array_size) 
   
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
!  Section 0.    Initialise array_size_count
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SL_TRACER2',zhook_in,zhook_handle)
      array_size_count=0

! ----------------------------------------------------------------------
! Section 2.1  carbon cycle.
! ----------------------------------------------------------------------
      array_size_count=0
      IF(l_CO2_interactive)THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = tdims%j_start,tdims%j_end
            DO i = tdims%i_start,tdims%i_end
              co2(i,j,k) = super_array(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF  ! l_CO2_interactive

! ----------------------------------------------------------------------
! Section 2.2  Soot cycle.
! ----------------------------------------------------------------------
      IF (l_Soot) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = tdims%j_start,tdims%j_end
            DO i = tdims%i_start,tdims%i_end
              soot_new(i,j,k) = super_array(i,j,k,array_size_count)
              soot_agd(i,j,k) = super_array(i,j,k,array_size_count+1)
              soot_cld(i,j,k) = super_array(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF  ! l_soot

! ----------------------------------------------------------------------
! Section 2.3  Biomass aerosol.
! ----------------------------------------------------------------------
      IF (l_Biomass) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = tdims%j_start,tdims%j_end
            DO i = tdims%i_start,tdims%i_end
              bmass_new(i,j,k) = super_array(i,j,k,array_size_count)
              bmass_agd(i,j,k) = super_array(i,j,k,array_size_count+1)
              bmass_cld(i,j,k) = super_array(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF  ! l_biomass

! ----------------------------------------------------------------------
! Section 2.4  sulphur cycle.
! ----------------------------------------------------------------------
      IF (l_Sulpc_so2) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = tdims%j_start,tdims%j_end
            DO i = tdims%i_start,tdims%i_end
              so4_aitken(i,j,k) = super_array(i,j,k,array_size_count)
              so4_accu(i,j,k) = super_array(i,j,k,array_size_count+1)
              so4_diss(i,j,k) = super_array(i,j,k,array_size_count+2)
              so2(i,j,k) = super_array(i,j,k,array_size_count+3)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +3

        IF(L_sulpc_nh3)THEN
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = tdims%j_start,tdims%j_end
              DO i = tdims%i_start,tdims%i_end
                nh3(i,j,k) = super_array(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        endif

        IF(L_sulpc_dms)THEN
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = tdims%j_start,tdims%j_end
              DO i = tdims%i_start,tdims%i_end
                dms(i,j,k) = super_array(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        endif

      END IF  ! l_sulpc_so2

! ----------------------------------------------------------------------
! Section 2.5  mineral dust.
! ----------------------------------------------------------------------
      IF (L_DUST) THEN
        array_size_count=array_size_count +1
        DO K = TDIMS%K_START, TDIMS%K_END
          DO j = tdims%j_start,tdims%j_end
            DO i = tdims%i_start,tdims%i_end
              DUST_DIV1(i,j,k) = super_array(i,j,k,array_size_count)
              DUST_DIV2(i,j,k) = super_array(i,j,k,array_size_count+1)
              IF(.NOT.l_twobin_dust) THEN
                DUST_DIV3(i,j,k) = super_array(i,j,k,array_size_count+2)
                DUST_DIV4(i,j,k) = super_array(i,j,k,array_size_count+3)
                DUST_DIV5(i,j,k) = super_array(i,j,k,array_size_count+4)
                DUST_DIV6(i,j,k) = super_array(i,j,k,array_size_count+5)
              END IF 
            END DO
          END DO
        END DO
        IF (l_twobin_dust) THEN
          array_size_count=array_size_count + 1
        ELSE
          array_size_count=array_size_count + 5
        END IF
      ENDIF  ! L_DUST

! ----------------------------------------------------------------------
! New addition  Fossil-fuel organic carbon aerosol
! ----------------------------------------------------------------------
      IF (L_OCFF) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = tdims%j_start,tdims%j_end
            DO i = tdims%i_start,tdims%i_end
              ocff_new(i,j,k) = super_array(i,j,k,array_size_count)
              ocff_agd(i,j,k) = super_array(i,j,k,array_size_count+1)
              ocff_cld(i,j,k) = super_array(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF  ! l_ocff

! ----------------------------------------------------------------------
! Section 2.5.1  Cariolle ozone tracer.
! ----------------------------------------------------------------------
      IF (L_USE_CARIOLLE) THEN
        array_size_count=array_size_count +1
        DO K = TDIMS%K_START, TDIMS%K_END
          DO j = tdims%j_start,tdims%j_end
            DO i = tdims%i_start,tdims%i_end
       OZONE_TRACER(i,j,k) = super_array(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      ENDIF  ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
! Section 2.5.2  Ammonium nitrate aerosol
! ----------------------------------------------------------------------
      IF (L_nitrate) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = tdims%j_start,tdims%j_end
            DO i = tdims%i_start,tdims%i_end
              nitr_acc(i,j,k) = super_array(i,j,k,array_size_count)
              nitr_diss(i,j,k)= super_array(i,j,k,array_size_count+1)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +1
      END IF  ! l_ocff

! ----------------------------------------------------------------------
!  ANY NEW NAMED TRACER SPECIES SHOULD BE ADDED HERE ^^^
!  and in the same location/order in sl_tracer1 and tr_set_phys
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 2.6.a  free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF ( tr_vars>0 ) THEN
        DO count=1,tr_vars
          array_size_count=array_size_count +1
          DO K = 1, tr_levels
            DO j = tdims%j_start,tdims%j_end
              DO i = tdims%i_start,tdims%i_end
              tracers(i,j,k,count) = super_array(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        END DO
      END IF  ! tr_vars > 0

! ----------------------------------------------------------------------
! Section 2.7  UKCA tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF ( tr_ukca > 0 ) THEN
        DO count=1,tr_ukca
          array_size_count=array_size_count +1
          DO k = tdims%k_start,tdims%k_end
            DO j = tdims%j_start,tdims%j_end
              DO i = tdims%i_start,tdims%i_end
              tracers_ukca(i,j,k,count) =                             &
                 super_array(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        END DO

      END IF  ! tr_ukca > 0

! ----------------------------------------------------------------------
! Section 2.99  Murk cycle.  This must be the last Full level field in
!                           the super_array
! ----------------------------------------------------------------------
      IF (L_Murk_advect) then
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = tdims%j_start,tdims%j_end
            DO i = tdims%i_start,tdims%i_end
              murk(i,j,k) = super_array(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF  ! L_Murk_advect

      IF (lhook) CALL dr_hook('SL_TRACER2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE unpack_tracer
      END MODULE unpack_tracer_mod
