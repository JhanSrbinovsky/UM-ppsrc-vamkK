! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

       SUBROUTINE TR_Reset(                                            &
                            super_array_size, super_tracer_phys,        &
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
                            tracer_phys2, tracer,                       &
                            ukca_tracer_phys2, tracer_ukca,             &
                            row_length, rows,                           &
                            model_levels, tr_levels, tr_vars, tr_ukca,  &
                            offx, offy                                  &
                            )

! Purpose: Interface routine to calculate increase in tracer fields
!          during the call to atmos_physics2
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Tracer Advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      USE dust_parameters_mod, ONLY: l_twobin_dust
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE atm_fields_bounds_mod

      Implicit None

! Arguments with Intent IN. ie: Input variables.

! Arguments with Intent IN. ie: Input variables.

      INTEGER, Intent(IN) :: row_length  ! number of points on a row
      INTEGER, Intent(IN) :: rows        ! number of rows in a theta field
      INTEGER, Intent(IN) :: model_levels   ! number of model levels         
      INTEGER, Intent(IN) :: tr_levels   ! number of tracer levels         
      INTEGER, Intent(IN) :: tr_vars   ! number of tracers         
      INTEGER, Intent(IN) :: tr_ukca   ! number of ukca tracers         
      INTEGER, Intent(IN) :: super_array_size  
      INTEGER, Intent(IN) :: offx      ! halo size in x direction
      INTEGER, Intent(IN) :: offy      ! halo size in y direction

      REAL, INTENT(IN) :: CO2       (tdims_s%i_start:tdims_s%i_end,                  &
                         tdims_s%j_start:tdims_s%j_end,                  &
                         tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: MURK       (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: SOOT_NEW   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: SOOT_AGD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: SOOT_CLD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: SO2        (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: SO4_AITKEN (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: SO4_ACCU   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: SO4_DISS   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: NH3        (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: DMS        (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: DUST_DIV1  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: DUST_DIV2  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: DUST_DIV3  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: DUST_DIV4  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: DUST_DIV5  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: DUST_DIV6  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: BMASS_NEW  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: BMASS_AGD  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: BMASS_CLD  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: OCFF_NEW   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: OCFF_AGD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: OCFF_CLD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: NITR_ACC (tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: NITR_DISS(tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::                                               &
        tracer      (1-offx:row_length+offx, 1-offy:rows+offy,          &
                     tr_levels,1:tr_vars),                              &
        tracer_ukca (tdims_s%i_start:tdims_s%i_end,                     &
                       tdims_s%j_start:tdims_s%j_end,                   &
                       tdims_s%k_start:tdims_s%k_end,1:tr_ukca)    
 
! Add cariolle specific parameters for ozone tracer     
      REAL, INTENT(IN) ::                                               &
        OZONE_TRACER(tdims_s%i_start:tdims_s%i_end,                     &
                           tdims_s%j_start:tdims_s%j_end,               &
                           tdims_s%k_start:tdims_s%k_end)                                          



      REAL, Intent(INOUT) ::                                       &
        super_tracer_phys(tdims%i_start:tdims%i_end,               &
                           tdims%j_start:tdims%j_end,              &
                           tdims%k_start:tdims%k_end,              &
                           super_array_size)                       &
     &, tracer_phys2(trdims_xstl%i_start:trdims_xstl%i_end,        &
                     trdims_xstl%j_start:trdims_xstl%j_end,        &
                     trdims_xstl%k_start:trdims_xstl%k_end,        &
                     tr_vars)                                      &
     &, ukca_tracer_phys2(trdims_xstl%i_start:trdims_xstl%i_end,   &
                          trdims_xstl%j_start:trdims_xstl%j_end,   &
                          trdims_xstl%k_start:trdims_xstl%k_end,   &
                          tr_ukca)

      LOGICAL, Intent(IN):: L_CO2_interactive 
      LOGICAL, Intent(IN):: L_Murk_advect     
      LOGICAL, Intent(IN):: L_soot
      LOGICAL, Intent(IN):: L_SULPC_SO2
      LOGICAL, Intent(IN):: L_SULPC_nh3
      LOGICAL, Intent(IN):: L_SULPC_dms
      LOGICAL, Intent(IN):: L_biomass  
      LOGICAL, Intent(IN):: L_dust
      LOGICAL, Intent(IN):: L_ocff
      LOGICAL, Intent(IN):: L_USE_CARIOLLE
      LOGICAL, Intent(IN):: L_nitrate

!      local variables
      INTEGER                                                           &
        i,j,k,l,count        !loop variables

      INTEGER :: array_size_count

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('TR_RESET',zhook_in,zhook_handle)
      array_size_count=0

! ----------------------------------------------------------------------
! Section 1.1  carbon cycle.
! ----------------------------------------------------------------------
      If (L_CO2_interactive) Then
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) = CO2(i,j,k)    &
                   -super_tracer_phys(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF   ! L_CO2_interactive

! ----------------------------------------------------------------------
! Section 1.2  Soot cycle.
! ----------------------------------------------------------------------
      If (L_soot) Then
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
         super_tracer_phys(i,j,k,array_size_count) =  soot_new(i,j,k)   &
                  -super_tracer_phys(i,j,k,array_size_count)
         super_tracer_phys(i,j,k,array_size_count+1) =  soot_agd(i,j,k) &
                  -super_tracer_phys(i,j,k,array_size_count+1)
         super_tracer_phys(i,j,k,array_size_count+2) =  soot_cld(i,j,k) &
                  -super_tracer_phys(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF    ! L_soot

! ----------------------------------------------------------------------
! Section 1.3  Biomass aerosol.
! ----------------------------------------------------------------------
      If (l_Biomass) Then
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  bmass_new(i,j,k)   &
                  -super_tracer_phys(i,j,k,array_size_count)
        super_tracer_phys(i,j,k,array_size_count+1) =  bmass_agd(i,j,k) &
                  -super_tracer_phys(i,j,k,array_size_count+1)
        super_tracer_phys(i,j,k,array_size_count+2) =  bmass_cld(i,j,k) &
                  -super_tracer_phys(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF  ! l_biomass

! ----------------------------------------------------------------------
! Section 1.4  sulphur cycle.
! ----------------------------------------------------------------------
      If (L_SULPC_SO2) Then
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  so4_aitken(i,j,k)  &
                   -super_tracer_phys(i,j,k,array_size_count)
        super_tracer_phys(i,j,k,array_size_count+1) =  so4_accu(i,j,k)  &
                   -super_tracer_phys(i,j,k,array_size_count+1)
        super_tracer_phys(i,j,k,array_size_count+2) =  so4_diss(i,j,k)  &
                   -super_tracer_phys(i,j,k,array_size_count+2)
        super_tracer_phys(i,j,k,array_size_count+3) =  so2(i,j,k)       &
                   -super_tracer_phys(i,j,k,array_size_count+3)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +3

        If (L_sulpc_nh3) Then
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
          super_tracer_phys(i,j,k,array_size_count) =  nh3(i,j,k)       &
                  -super_tracer_phys(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        END IF  ! L_sulpc_nh3

        If (L_sulpc_dms) Then
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  dms(i,j,k)         &
                  -super_tracer_phys(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        END IF  ! L_sulpc_dms
      END IF  ! L_sulpc_SO2

! ----------------------------------------------------------------------
! Section 1.5  mineral dust.
! ----------------------------------------------------------------------
      If (L_dust) Then
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  dust_div1(i,j,k)   &
                  -super_tracer_phys(i,j,k,array_size_count)
        super_tracer_phys(i,j,k,array_size_count+1) =  dust_div2(i,j,k) &
                  -super_tracer_phys(i,j,k,array_size_count+1)
               IF(.NOT.l_twobin_dust) THEN
        super_tracer_phys(i,j,k,array_size_count+2) =  dust_div3(i,j,k) &
                  -super_tracer_phys(i,j,k,array_size_count+2)
        super_tracer_phys(i,j,k,array_size_count+3) =  dust_div4(i,j,k) &
                  -super_tracer_phys(i,j,k,array_size_count+3)
        super_tracer_phys(i,j,k,array_size_count+4) =  dust_div5(i,j,k) &
                  -super_tracer_phys(i,j,k,array_size_count+4)
        super_tracer_phys(i,j,k,array_size_count+5) =  dust_div6(i,j,k) &
                  -super_tracer_phys(i,j,k,array_size_count+5)
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
! New addition. Fossil-fuel organic carbon aerosol.
! ----------------------------------------------------------------------
      If (l_ocff) Then
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  ocff_new(i,j,k)    &
                  -super_tracer_phys(i,j,k,array_size_count)
        super_tracer_phys(i,j,k,array_size_count+1) =  ocff_agd(i,j,k)  &
                  -super_tracer_phys(i,j,k,array_size_count+1)
        super_tracer_phys(i,j,k,array_size_count+2) =  ocff_cld(i,j,k)  &
                  -super_tracer_phys(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF  ! l_ocff

! ----------------------------------------------------------------------
! Section 1.5.1  cariolle ozone tracer.
! ----------------------------------------------------------------------
      If (L_USE_CARIOLLE) Then
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  OZONE_TRACER(i,j,k)&
                 -super_tracer_phys(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF    ! L_USE_CARIOLLE

! ---------------------------------------------------------------------- 
! New addition. Ammonium nitrate aerosol.  
! ----------------------------------------------------------------------  
      If (l_nitrate) Then  
        array_size_count=array_size_count +1  
         DO k = tdims%k_start, tdims%k_end  
          DO j = 1, rows  
            DO i = 1, row_length  
        super_tracer_phys(i,j,k,array_size_count) =  nitr_acc(i,j,k)    &  
                  -super_tracer_phys(i,j,k,array_size_count)  
        super_tracer_phys(i,j,k,array_size_count+1) = nitr_diss(i,j,k)  &  
                  -super_tracer_phys(i,j,k,array_size_count+1)  
            END DO  
          END DO  
        END DO  
        array_size_count=array_size_count +1 
      END IF  ! l_nitrate 

! ----------------------------------------------------------------------
!  ANY NEW NAMED TRACER SPECIES SHOULD BE ADDED HERE
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 1.6.a  free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      If (model_levels == tr_levels .and. tr_vars > 0) Then
        DO count=1,tr_vars
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) = tracer(i,j,k,count) &
                  -super_tracer_phys(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        END DO
      END IF  ! tr_vars > 0

! ----------------------------------------------------------------------
! Section 1.7  ukca tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF ( tr_ukca > 0) THEN
        DO count=1,tr_ukca
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
                super_tracer_phys(i,j,k,array_size_count) =             &
                   tracer_ukca(i,j,k,count)                             &
                  -super_tracer_phys(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        END DO
      END IF  ! tr_ukca > 0

! ----------------------------------------------------------------------
! Section 1.99  Murk cycle.  This must be the last Full level field in
!                           the super_array
! ----------------------------------------------------------------------
      If (L_Murk_advect) Then
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  murk(i,j,k)        &
                  -super_tracer_phys(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF  ! L_Murk_advect

! ----------------------------------------------------------------------
! Section 1.7.b  free tracers  (model_levels/=tr_levels)
! ----------------------------------------------------------------------
      If (tr_vars > 0 .and. tr_levels /= model_levels) Then
        DO l = 1 , tr_vars
          DO k = 1, tr_levels
            DO j = 1, rows
              DO i = 1, row_length
                tracer_phys2(i,j,k,l) = tracer(i,j,k,l)                 &
                  -tracer_phys2(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
      END IF  ! tr_vars > 0

! ----------------------------------------------------------------------
! Section 1.7.c  UKCA tracers  (model_levels/=tr_levels)
! ----------------------------------------------------------------------
      If (tr_ukca > 0 .and. tr_levels /= model_levels) Then
        DO l = 1 , tr_ukca
          DO k = 1, tr_levels
            DO j = 1, rows
              DO i = 1, row_length
                ukca_tracer_phys2(i,j,k,l) = tracer_ukca(i,j,k,l)       &
                  -ukca_tracer_phys2(i,j,k,l)
              END DO
            END DO
          END DO
        END DO
      END IF  ! tr_ukca > 0


      IF (lhook) CALL dr_hook('TR_RESET',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE TR_Reset

