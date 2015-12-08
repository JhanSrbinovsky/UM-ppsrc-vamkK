! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE tr_set_phys_mod
IMPLICIT NONE
CONTAINS

      SUBROUTINE tr_set_phys_4A(                                      &
                          super_array_size, super_tracer_phys,        &
                          l_co2_interactive, co2,                     &
                          l_murk_advect, murk,                        &
                          l_soot, soot_new, soot_agd, soot_cld,       &
                          l_sulpc_so2, so2, so4_aitken, so4_accu,     &
                                       so4_diss,                      &
                          l_sulpc_nh3, nh3,                           &
                          l_sulpc_dms, dms,                           &
                          l_dust, dust_div1, dust_div2, dust_div3,    &
                                  dust_div4, dust_div5, dust_div6,    &
                          l_biomass, bmass_new, bmass_agd, bmass_cld, &
                          l_ocff, ocff_new, ocff_agd, ocff_cld,       &
                          l_nitrate, nitr_acc, nitr_diss,             &
                          l_use_cariolle, ozone_tracer,               &
                          tracer_phys1, tracer,                       &
                          ukca_tracer_phys1, tracer_ukca,             &
                          row_length, rows,                           &
                          model_levels, tr_levels, tr_vars, tr_ukca,  &
                          offx, offy, model_domain,                   &
                          l_init, supertrdims                         &
                            )

! Purpose: Interface routine to initialise tracer fields
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
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      USE atm_fields_bounds_mod

      IMPLICIT NONE


      TYPE (array_dims)  supertrdims 

! Arguments with INTENT IN. ie: Input variables.

! Arguments with INTENT IN. ie: Input variables.

      INTEGER                                                         &
              ! model dimensions
        row_length                                                    &
                         ! number of points on a row
      , rows                                                          &
                         ! number of rows in a theta field
      , model_levels                                                  &
                         ! number of model levels
      , tr_levels                                                     &
                         ! number of tracer levels
      , tr_vars                                                       &
                         ! number of tracers
      , tr_ukca                                                       &
                         ! number of ukca tracers
      , super_array_size                                              &
      , offx                                                          &
      , offy

      INTEGER, INTENT(IN) :: model_domain

      REAL, INTENT(IN) :: co2                                         &
                        (tdims_s%i_start:tdims_s%i_end,               &
                         tdims_s%j_start:tdims_s%j_end,               &
                         tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: murk                                       &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: soot_new                                   &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: soot_agd                                   &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: soot_cld                                   &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: so2                                        &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: so4_aitken                                 &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: so4_accu                                   &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: so4_diss                                   &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: nh3                                        &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dms                                        &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: dust_div1                                  &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div2                                  &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div3                                  &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div4                                  &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div5                                  &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div6                                  &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: bmass_new                                  &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: bmass_agd                                  &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: bmass_cld                                  &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: ocff_new                                   &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: ocff_agd                                   &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: ocff_cld                                   &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: nitr_acc                                   &
                       (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: nitr_diss                                  &
                       (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::                                             &
        tracer      (1-offx:row_length+offx, 1-offy:rows+offy,        &
                     tr_levels,1:tr_vars),                            &
        tracer_ukca(tdims_s%i_start:tdims_s%i_end,                    &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end,1:tr_ukca)

! Add cariolle specific parameters for ozone tracer    
      REAL, INTENT(IN) ::                                             &
        ozone_tracer(tdims_s%i_start:tdims_s%i_end,                   &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end) 


      REAL, INTENT(OUT) :: tracer_phys1                               &
                        (trdims_ltl%i_start:trdims_ltl%i_end,         &
                         trdims_ltl%j_start:trdims_ltl%j_end,         &
                         trdims_ltl%k_start:trdims_ltl%k_end, tr_vars)

      REAL, INTENT(OUT) :: ukca_tracer_phys1                          &
                         (trdims_ltl%i_start:trdims_ltl%i_end,        &
                          trdims_ltl%j_start:trdims_ltl%j_end,        &
                          trdims_ltl%k_start:trdims_ltl%k_end,        &
                          tr_ukca)

      REAL, INTENT(OUT) :: super_tracer_phys                          &
                          (supertrdims%i_start:supertrdims%i_end,     &
                           supertrdims%j_start:supertrdims%j_end,     &
                           supertrdims%k_start:supertrdims%k_end,     &
                           super_array_size) 

      LOGICAL                                                         &
        l_co2_interactive,                                            &
        l_murk_advect,                                                &
        l_soot,                                                       &
        l_sulpc_so2,                                                  &
        l_sulpc_nh3,                                                  &
        l_sulpc_dms,                                                  &
        l_biomass,                                                    &
        l_dust,                                                       &
        l_ocff,                                                       &
        l_use_cariolle,                                               &
        l_nitrate

      LOGICAL :: L_init   ! flag for setting halo values 

!      local variables
      INTEGER                                                         &
        i,j,k,l,counter        !loop variables

      INTEGER :: array_size_count

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('TR_SET_PHYS_4A',zhook_in,zhook_handle)
      array_size_count=0

! ----------------------------------------------------------------------
! Section 1.1  carbon cycle.
! ----------------------------------------------------------------------
      IF (L_CO2_interactive) THEN
        array_size_count=array_size_count +1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) = CO2(i,j,k)
            END DO
          END DO
        END DO
      END IF   ! L_CO2_interactive

! ----------------------------------------------------------------------
! Section 1.2  Soot cycle.
! ----------------------------------------------------------------------
      IF (l_soot) THEN
        array_size_count=array_size_count +1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count)   =  soot_new(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  soot_agd(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+2) =  soot_cld(i,j,k)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF    ! L_soot

! ----------------------------------------------------------------------
! Section 1.3  Biomass aerosol.
! ----------------------------------------------------------------------
      IF (l_biomass) THEN
        array_size_count=array_size_count +1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count)   =  bmass_new(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  bmass_agd(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+2) =  bmass_cld(i,j,k)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF  ! l_biomass

! ----------------------------------------------------------------------
! Section 1.4  sulphur cycle.
! ----------------------------------------------------------------------
      IF (l_sulpc_so2) THEN
        array_size_count=array_size_count +1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count)   =  so4_aitken(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  so4_accu(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+2) =  so4_diss(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+3) =  so2(i,j,k)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +3

        IF(l_sulpc_nh3) THEN
          array_size_count=array_size_count +1
          DO k = tdims_s%k_start, tdims_s%k_end
            DO j = 1, rows
              DO i = 1, row_length
                super_tracer_phys(i,j,k,array_size_count) =  nh3(i,j,k)
              END DO
            END DO
          END DO
        END IF  ! L_sulpc_nh3

        IF(l_sulpc_dms) THEN
          array_size_count=array_size_count +1
          DO k = tdims_s%k_start, tdims_s%k_end
            DO j = 1, rows
              DO i = 1, row_length
                super_tracer_phys(i,j,k,array_size_count) =  dms(i,j,k)
              END DO
            END DO
          END DO
        END IF  ! L_sulpc_dms
      END IF  ! L_sulpc_SO2

! ----------------------------------------------------------------------
! Section 1.5  Mineral dust.
! ----------------------------------------------------------------------
      IF (l_dust) THEN
        array_size_count=array_size_count +1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  dust_div1(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  dust_div2(i,j,k)
              IF(.NOT.l_twobin_dust) THEN
                super_tracer_phys(i,j,k,array_size_count+2) =  dust_div3(i,j,k)
                super_tracer_phys(i,j,k,array_size_count+3) =  dust_div4(i,j,k)
                super_tracer_phys(i,j,k,array_size_count+4) =  dust_div5(i,j,k)
                super_tracer_phys(i,j,k,array_size_count+5) =  dust_div6(i,j,k)
              END IF
            END DO
          END DO
        END DO

        IF(l_twobin_dust) THEN
          array_size_count=array_size_count + 1
        ELSE
          array_size_count=array_size_count + 5
        END IF
      END IF    ! L_dust

! ----------------------------------------------------------------------
! Section 1.6  Fossil-fuel organic carbon aerosol
! ----------------------------------------------------------------------
      IF (l_ocff) THEN
        array_size_count=array_size_count +1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  ocff_new(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  ocff_agd(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+2) =  ocff_cld(i,j,k)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF    ! L_ocff

! ----------------------------------------------------------------------
! Section 1.7  Cariolle ozone tracer.
! ----------------------------------------------------------------------
      IF (l_use_cariolle) THEN
        array_size_count=array_size_count +1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  OZONE_TRACER(i,j,k)
            END DO
          END DO
        END DO
      END IF    ! L_USE_CARIOLLE

! ----------------------------------------------------------------------  
! Section 1.8  Ammonium nitrate aerosol. 
! ----------------------------------------------------------------------  
      IF (l_nitrate) THEN  
        array_size_count=array_size_count +1  
        DO k = tdims_s%k_start, tdims_s%k_end  
          DO j = 1, rows  
            DO i = 1, row_length  
              super_tracer_phys(i,j,k,array_size_count) = nitr_acc(i,j,k)  
              super_tracer_phys(i,j,k,array_size_count+1) = nitr_diss(i,j,k)  
            END DO  
          END DO  
        END DO  
        array_size_count=array_size_count +1 
      END IF    ! L_nitrate 

! ----------------------------------------------------------------------
!  ANY NEW NAMED TRACER SPECIES SHOULD BE ADDED HERE
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Section 1.33  Free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF (model_levels==tr_levels.AND. tr_vars>0) THEN
        DO counter=1,tr_vars
          array_size_count=array_size_count +1
          DO K = 1, MODEL_LEVELS
            DO J = 1, ROWS
              DO I = 1, ROW_LENGTH
                super_tracer_phys(i,j,k,array_size_count) =          &
                                                  tracer(i,j,k,counter)
              END DO
            END DO
          END DO
        END DO

!       IF(L_init) THEN
!       ! Call SWAPBOUNDS and set external halos for free tracers
! DEPENDS ON: swap_bounds
!         Call Swap_Bounds(                                          &
! &     super_tracer_phys(1-halo_i_in,1-halo_j_in,1,              &
! &                       super_array_size-tr_vars+1),            &
! &       row_length, rows,                                       &
! &       model_levels*tr_vars,                                   &
! &       halo_i_in, halo_j_in, fld_type_p,  .FALSE. )
!         IF (model_domain == mt_lam) THEN
! DEPENDS ON: set_external_halos
!           Call SET_EXTERNAL_HALOS(                                 &
!    &      super_tracer_phys(1-halo_i_in,1-halo_j_in,1,             &
!!   &                        super_array_size-tr_vars+1),           &
!    &      row_length, rows,                                        &
!    &      tr_vars*model_levels, halo_i_in, halo_j_in, 0.0)
!         END IF
!       END IF ! L_INIT
      END IF  ! tr_vars>0

! ----------------------------------------------------------------------
! Section 1.34  UKCA tracers  
! ----------------------------------------------------------------------
      IF ( tr_ukca > 0) THEN
        DO counter=1,tr_ukca
          array_size_count=array_size_count +1
          DO k = tdims_s%k_start,tdims_s%k_end
            DO j = 1, rows
              DO i = 1, row_length
                super_tracer_phys(i,j,k,array_size_count) =          &
                     tracer_ukca(i,j,k,counter)
              END DO
            END DO
          END DO
        END DO
      END IF     ! tr_ukca>0

! ----------------------------------------------------------------------
! Section 1.99  Murk cycle.  This must be the last Full level field in
!                            the super_array
! ----------------------------------------------------------------------
      IF (l_murk_advect) THEN
        array_size_count=array_size_count +1
        DO k = tdims_s%k_start, tdims_s%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  murk(i,j,k)
            END DO
          END DO
        END DO
      END IF  ! L_Murk_advect

       super_tracer_phys(:,:,0,:)=super_tracer_phys(:,:,1,:)

      IF(l_init)THEN
! ----------------------------------------------------------------------
! Call SWAPBOUNDS and set external halos for all arrays if required
! ----------------------------------------------------------------------

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                             &
             super_tracer_phys,                                       &
              row_length, rows,                                       &
         (supertrdims%k_end-supertrdims%k_start+1)*array_size_count,  &
          supertrdims%halo_i,supertrdims%halo_j, fld_type_p,  .FALSE.)

        IF (model_domain == mt_lam) THEN
! DEPENDS ON: set_external_halos
          CALL set_external_halos(                                    &
              super_tracer_phys,                                      &
              row_length, rows,                                       &
              (supertrdims%k_end-supertrdims%k_start+1)*array_size_count,  &
              supertrdims%halo_i, supertrdims%halo_j, 0.0)
        END IF

      END IF  !L_INIT
        IF (l_murk_advect) THEN
          IF (model_domain == mt_lam) THEN
! DEPENDS ON: fill_external_halos
            CALL fill_external_halos(                                &
             super_tracer_phys(supertrdims%i_start,                  &
                               supertrdims%j_start,                  &
                               supertrdims%k_start,                  &
                               array_size_count)                     &
                         , row_length, rows,                         &
                           supertrdims%k_end-supertrdims%k_start+1,  &
                           supertrdims%halo_i, supertrdims%halo_j)
          END IF
        END IF  ! L_Murk_advect
      
!     Endif  !L_init
! ----------------------------------------------------------------------
! Section 1.6.b  free tracers  (model_levels/=tr_levels)
! ----------------------------------------------------------------------
      IF (tr_vars > 0.AND. tr_levels /= model_levels) THEN
        DO l = 1 , tr_vars
          DO k = trdims_ltl%k_start,trdims_ltl%k_end
            DO j = 1, rows
              DO i = 1, row_length
                tracer_phys1(i,j,k,l) = TRACER(i,j,k,l)
              END DO
            END DO
          END DO
        END DO

        tracer_phys1(:,:,0,:)=tracer_phys1(:,:,1,:)
        IF(l_init)THEN
! Call SWAPBOUNDS and set external halos for free tracers if necessary
! DEPENDS ON: swap_bounds
          CALL swap_bounds(                                           &
               tracer_phys1,                                          &
               row_length, rows,                                      &
              (trdims_ltl%k_end-trdims_ltl%k_start+1)*tr_vars,        &
               supertrdims%halo_i, supertrdims%halo_j , fld_type_p,   &
               .FALSE. )

          IF (model_domain == mt_lam) THEN
! DEPENDS ON: set_external_halos
           CALL set_external_halos(                                   &
              tracer_phys1,                                           &
              row_length, rows,                                       &
              (trdims_ltl%k_end-trdims_ltl%k_start+1)*tr_vars,        &
              supertrdims%halo_i, supertrdims%halo_j, 0.0)
          END IF
        END IF !L_init
      END IF  ! tr_vars > 0


! ----------------------------------------------------------------------
! Section 1.6.c  UKCA tracers  (model_levels/=tr_levels)
! ----------------------------------------------------------------------
      IF (tr_ukca > 0.AND. tr_levels /= model_levels) THEN
        DO l = 1 , tr_ukca
          DO k = 1, tr_levels
            DO j = 1, rows
              DO i = 1, row_length
                ukca_tracer_phys1(i,j,k,l) = TRACER_UKCA(i,j,k,l)
              END DO
            END DO
          END DO
        END DO

        IF(l_init)THEN
! Call SWAPBOUNDS and set external halos for UKCA tracers if necessary
! DEPENDS ON: swap_bounds
          CALL swap_bounds(                                          &
               ukca_tracer_phys1,                                    &
               row_length, rows, tr_levels*tr_ukca,                  &
               supertrdims%halo_i, supertrdims%halo_j, fld_type_p,   &
                .FALSE. )
          IF (model_domain == mt_lam) THEN
! DEPENDS ON: set_external_halos
           CALL set_external_halos(ukca_tracer_phys1, row_length,    &
              rows, tr_ukca*tr_levels, supertrdims%halo_i,           &
              supertrdims%halo_j, 0.0)
          END IF
        END IF !L_init
      END IF  ! tr_ukca > 0

      IF (lhook) CALL dr_hook('TR_SET_PHYS_4A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE TR_Set_Phys_4A
      END MODULE tr_set_phys_mod
