! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE TR_ReSet_mod
IMPLICIT NONE
CONTAINS

       SUBROUTINE TR_Reset_4A(                                          &
                            super_array_size, super_tracer_phys,        &
                            L_CO2_interactive, CO2,                     &
                            L_murk_advect, murk,                        &
                            L_Soot, soot_new, soot_agd, soot_cld,       &
                            L_sulpc_SO2, SO2, SO4_aitken, so4_accu,     &
                                         so4_diss,                      &
                            L_sulpc_nh3, nh3,                           &
                            L_sulpc_dms, dms,                           &
                            L_dust, dust_div1, dust_div2, dust_div3,    &
                                    dust_div4, dust_div5, dust_div6,    &
                            L_biomass, bmass_new, bmass_agd, bmass_cld, &
                            L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                            L_nitrate, nitr_acc, nitr_diss,             &
                            l_use_cariolle, ozone_tracer,               &
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

      IMPLICIT NONE

! Arguments with INTENT IN. ie: Input variables.

! Arguments with INTENT IN. ie: Input variables.

      INTEGER, INTENT(IN) :: row_length  ! number of points on a row
      INTEGER, INTENT(IN) :: rows        ! number of rows in a theta field
      INTEGER, INTENT(IN) :: model_levels   ! number of model levels         
      INTEGER, INTENT(IN) :: tr_levels   ! number of tracer levels         
      INTEGER, INTENT(IN) :: tr_vars   ! number of tracers         
      INTEGER, INTENT(IN) :: tr_ukca   ! number of ukca tracers         
      INTEGER, INTENT(IN) :: super_array_size  
      INTEGER, INTENT(IN) :: offx      ! halo size in x direction
      INTEGER, INTENT(IN) :: offy      ! halo size in y direction

      REAL, INTENT(IN) :: CO2       (tdims_s%i_start:tdims_s%i_end,      &
                         tdims_s%j_start:tdims_s%j_end,                  &
                         tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: murk       (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: soot_new   (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: soot_agd   (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: soot_cld   (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: SO2        (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: SO4_aitken (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: so4_accu   (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: so4_diss   (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: NH3        (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: DMS        (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: dust_div1  (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div2  (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div3  (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div4  (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div5  (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: dust_div6  (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: bmass_new  (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: bmass_agd  (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: bmass_cld  (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: ocff_new   (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: ocff_agd   (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: ocff_cld   (tdims_s%i_start:tdims_s%i_end,    &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN)  :: nitr_acc (tdims_s%i_start:tdims_s%i_end,      &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN)  :: nitr_diss(tdims_s%i_start:tdims_s%i_end,      &
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
        ozone_tracer(tdims_s%i_start:tdims_s%i_end,                     &
                           tdims_s%j_start:tdims_s%j_end,               &
                           tdims_s%k_start:tdims_s%k_end)                      



      REAL, INTENT(INOUT) ::                                       &
        super_tracer_phys(tdims%i_start:tdims%i_end,               &
                           tdims%j_start:tdims%j_end,              &
                           tdims%k_start:tdims%k_end,              &
                           super_array_size)                       &
      , tracer_phys2(trdims_xstl%i_start:trdims_xstl%i_end,        &
                     trdims_xstl%j_start:trdims_xstl%j_end,        &
                     trdims_xstl%k_start:trdims_xstl%k_end,        &
                     tr_vars)                                      &
      , ukca_tracer_phys2(trdims_xstl%i_start:trdims_xstl%i_end,   &
                          trdims_xstl%j_start:trdims_xstl%j_end,   &
                          trdims_xstl%k_start:trdims_xstl%k_end,   &
                          tr_ukca)

      LOGICAL, INTENT(IN):: L_CO2_interactive 
      LOGICAL, INTENT(IN):: L_murk_advect     
      LOGICAL, INTENT(IN):: L_soot
      LOGICAL, INTENT(IN):: L_sulpc_SO2
      LOGICAL, INTENT(IN):: L_SULPC_nh3
      LOGICAL, INTENT(IN):: L_SULPC_dms
      LOGICAL, INTENT(IN):: L_biomass  
      LOGICAL, INTENT(IN):: L_dust
      LOGICAL, INTENT(IN):: L_ocff
      LOGICAL, INTENT(IN):: l_use_cariolle
      LOGICAL, INTENT(IN):: L_nitrate

!      local variables
      INTEGER                                                           &
        i,j,k,l,counter        !loop variables

      INTEGER :: array_size_count

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('TR_RESET_4A',zhook_in,zhook_handle)
      array_size_count=0

! ----------------------------------------------------------------------
! Section 1.1  carbon cycle.
! ----------------------------------------------------------------------
      IF (L_CO2_interactive) THEN
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
      IF (L_soot) THEN
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
      IF (l_Biomass) THEN
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
      IF (L_sulpc_SO2) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  SO4_aitken(i,j,k)  &
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

        IF (L_sulpc_nh3) THEN
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

        IF (L_sulpc_dms) THEN
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
      IF (L_dust) THEN
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
      IF (l_ocff) THEN
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
      IF (l_use_cariolle) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  ozone_tracer(i,j,k)&
                 -super_tracer_phys(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF    ! l_use_cariolle

! ---------------------------------------------------------------------- 
! New addition. Ammonium nitrate aerosol.  
! ----------------------------------------------------------------------  
      IF (l_nitrate) THEN  
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
      IF (model_levels == tr_levels .AND. tr_vars > 0) THEN
        DO counter=1,tr_vars
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) = tracer(i,j,k,counter) &
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
        DO counter=1,tr_ukca
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
                super_tracer_phys(i,j,k,array_size_count) =             &
                   tracer_ukca(i,j,k,counter)                             &
                  -super_tracer_phys(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        END DO
      END IF  ! tr_ukca > 0

! ----------------------------------------------------------------------
! Section 1.99  murk cycle.  This must be the last Full level field in
!                           the super_array
! ----------------------------------------------------------------------
      IF (L_murk_advect) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  murk(i,j,k)        &
                  -super_tracer_phys(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF  ! L_murk_advect

 super_tracer_phys(:,:,0,1:array_size_count) =                          &
                        super_tracer_phys(:,:,1,1:array_size_count) 

! ----------------------------------------------------------------------
! Section 1.7.b  free tracers  (model_levels/=tr_levels)
! ----------------------------------------------------------------------
      IF (tr_vars > 0 .AND. tr_levels /= model_levels) THEN
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
      IF (tr_ukca > 0 .AND. tr_levels /= model_levels) THEN
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


      IF (lhook) CALL dr_hook('TR_RESET_4A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE TR_Reset_4A
END MODULE TR_ReSet_mod

