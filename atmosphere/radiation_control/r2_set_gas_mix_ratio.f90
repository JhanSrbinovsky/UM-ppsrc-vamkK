! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to set the mixing ratios of gases.

! Purpose:
!   The full array of mass mixing ratios of gases is filled.

! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
      SUBROUTINE r2_set_gas_mix_ratio(ierr                              &
         , n_profile, nlevs, n_layer, nwet, nozone                      &
         , i_gather, l_extra_top                                        &
         , n_absorb, type_absorb                                        &
         , l_n2o, l_ch4, l_cfc11, l_cfc12, l_o2                         &
         , l_cfc113, l_cfc114, l_hcfc22, l_hfc125, l_hfc134a            &
         , h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio                   &
         , c11_mix_ratio, c12_mix_ratio, o2_mix_ratio                   &
         , c113_mix_ratio, c114_mix_ratio, hcfc22_mix_ratio             &
         , hfc125_mix_ratio, hfc134a_mix_ratio                          &
         , gas_mix_ratio                                                &
         , co2_dim1, co2_dim2, co2_3d, l_co2_3d                         &
         , npd_field, npd_profile, npd_layer, npd_species, first_layer  &
! chemical greenhouse gas fields
         , ngrgas, grgas_field                                          &
         )

      USE rad_pcf
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE ereport_mod, ONLY : ereport
      USE gasid3a_mod, ONLY: npd_gases, ip_h2o, ip_co2, ip_o3, ip_n2o,   &
                             ip_co, ip_ch4, ip_o2, ip_no, ip_so2, ip_no2,&
                             ip_nh3, ip_hno3, ip_n2, ip_cfc11, ip_cfc12, &
                             ip_cfc113, ip_hcfc22, ip_hfc125, ip_hfc134a,&
                             ip_cfc114
      USE ukca_feedback_mod, ONLY: p_ch4, p_n2o, p_f11, p_f12, p_f113,   &
                                   p_f22 

      IMPLICIT NONE

!     DUMMY ARGUMENTS.

      INTEGER                                                           &
                !, INTENT(OUT)
           ierr
!             ERROR FLAG

!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
           npd_field                                                    &
!             SIZE OF ARRAY FROM UM
         , npd_profile                                                  &
!             SIZE OF ARRAY
         , npd_layer                                                    &
!             SIZE OF ARRAY
         , npd_species                                                  &
!             SIZE OF ARRAY
         , first_layer
!             First layer for some variables
!             0 for flux_calc(3A/C), 1 for radiance_calc(3Z)

!     SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
           n_profile                                                    &
!             NUMBER OF PROFILES
         , nlevs                                                        &
!             Number of layers in the main model
         , n_layer                                                      &
!             Number of radiative layers
         , nwet                                                         &
!             NUMBER OF WET LEVELS
         , nozone
!             NUMBER OF OZONE LEVELS

!     GATHERING ARRAY:
      INTEGER                                                           &
                !, INTENT(IN)
           i_gather(npd_field)
!             LIST OF POINTS TO BE GATHERED
      LOGICAL                                                           &
                !, INTENT(IN)
           l_extra_top
!             Flag to use an extra top layer in radiative calculations

!     TYPES OF GASES:
      INTEGER                                                           &
                !, INTENT(IN)
           n_absorb                                                     &
!             NUMBER OF ABSORBERS
         , type_absorb(npd_species)
!             TYPES OF ABSORBERS

!     FLAGS FOR MINOR GASES:
      LOGICAL                                                           &
                !,INTENT(IN)
           l_n2o                                                        &
!             FLAG FOR NITROUS OXIDE
         , l_ch4                                                        &
!             FLAG FOR METHANE
         , l_cfc11                                                      &
!             FLAG FOR CFC11
         , l_cfc12                                                      &
!             FLAG FOR CFC12
         , l_o2                                                         &
!             FLAG FOR O2
         , l_cfc113                                                     &
!             FLAG FOR CFC113
         , l_cfc114                                                     &
!             FLAG FOR CFC114
         , l_hcfc22                                                     &
!             FLAG FOR HCFC22
         , l_hfc125                                                     &
!             FLAG FOR HFC125
         , l_hfc134a
!             FLAG FOR HFC134A

!     MIXING RATIOS SUPPLIED:
      INTEGER  co2_dim1, co2_dim2   ! dimensions of CO2_3D field
      LOGICAL  l_co2_3d    !  controls use of 3D co2 field
      REAL                                                              &
                !, INTENT(IN)
           h2o(npd_field, nwet)                                         &
!             MASS MIXING RATIO OF WATER VAPOUR
         , co2                                                          &
!             MASS MIXING RATIO OF CARBON DIOXIDE
         , co2_3d(co2_dim1, co2_dim2)                                   &
!             3D MASS MIXING RATIO OF CO2 (full field)
         , o3(npd_field, nozone)                                        &
!             MASS MIXING RATIO OF OZONE
         , n2o_mix_ratio                                                &
!             MASS MIXING RATIO OF NITROUS OXIDE
         , ch4_mix_ratio                                                &
!             MASS MIXING RATIO OF METHANE
         , c11_mix_ratio                                                &
!             MASS MIXING RATIO OF CFC11
         , c12_mix_ratio                                                &
!             MASS MIXING RATIO OF CFC12
         , o2_mix_ratio                                                 &
!             MASS MIXING RATIO OF O2
         , c113_mix_ratio                                               &
!             MASS MIXING RATIO OF CFC113
         , c114_mix_ratio                                               &
!             MASS MIXING RATIO OF CFC114
         , hcfc22_mix_ratio                                             &
!             MASS MIXING RATIO OF HCFC22
         , hfc125_mix_ratio                                             &
!             MASS MIXING RATIO OF HFC125
         , hfc134a_mix_ratio
!             MASS MIXING RATIO OF HFC134A

!  ngrgas is either non-zero, if called from LWRAD, or 0
!  if called from swrad. In the latter case, the supplied
!  fields are ignored.
      INTEGER, INTENT(IN) :: ngrgas
      REAL, INTENT(IN) :: grgas_field(npd_field, nlevs, ngrgas)

!     ARRAY OF ASSIGNED MXING RATIOS:
      REAL                                                              &
                !, INTENT(OUT)
           gas_mix_ratio(npd_profile, first_layer:npd_layer, npd_species)
!             MIXING RATIOS

!     LOCAL VARIABLES.

!     POINTERS TO GASES:
      INTEGER                                                           &
           iump_h2o                                                     &
!             POINTER TO WATER VAPOUR
         , iump_co2                                                     &
!             POINTER TO CARBON DIOXIDE
         , iump_o3                                                      &
!             POINTER TO OZONE
         , iump_n2o                                                     &
!             POINTER TO NITOUS OXIDE
         , iump_ch4                                                     &
!             POINTER TO METHANE
         , iump_cfc11                                                   &
!             POINTER TO CFC11
         , iump_cfc12                                                   &
!             POINTER TO CFC12
         , iump_o2                                                      &
!             POINTER TO O2
         , iump_cfc113                                                  &
!             POINTER TO CFC113
         , iump_cfc114                                                  &
!             POINTER TO CFC114
         , iump_hcfc22                                                  &
!             POINTER TO HCFC22
         , iump_hfc125                                                  &
!             POINTER TO HFC125
         , iump_hfc134a
!             POINTER TO HFC134A
      INTEGER                                                           &
           i                                                            &
!             LOOP VARIABLE
         , l                                                            &
!             LOOP VARIABLE
         , lg                                                           &
!             CORRESPONDING UNGATHERED INDEX
         , i_top_copy
!             Topmost layer where properties are set by copying the
!             input fields.


      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'r2_set_gas_mix_ratio'
      CHARACTER (LEN=240) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('R2_SET_GAS_MIX_RATIO',zhook_in,zhook_handle)

!     MATCH THE INDEXING NUMBERS OF GASEOUS SPECIES IN THE SPECTRAL
!     FILE WITH ACTUAL TYPES OF GASES KNOWN TO THE UM.

!     SET ALL POINTERS TO 0 INITIALLY TO FLAG MISSING GASES.
      iump_h2o=0
      iump_co2=0
      iump_o3=0
      iump_n2o=0
      iump_ch4=0
      iump_cfc11=0
      iump_cfc12=0
      iump_o2=0
      iump_cfc113=0
      iump_cfc114=0
      iump_hcfc22=0
      iump_hfc125=0
      iump_hfc134a=0


      DO i=1, n_absorb

         IF (type_absorb(i) == ip_h2o) THEN
            iump_h2o=i
         ELSE IF (type_absorb(i) == ip_co2) THEN
            iump_co2=i
         ELSE IF (type_absorb(i) == ip_o3) THEN
            iump_o3=i
         ELSE IF (type_absorb(i) == ip_n2o) THEN
            iump_n2o=i
         ELSE IF (type_absorb(i) == ip_ch4) THEN
            iump_ch4=i
         ELSE IF (type_absorb(i) == ip_cfc11) THEN
            iump_cfc11=i
         ELSE IF (type_absorb(i) == ip_cfc12) THEN
            iump_cfc12=i
         ELSE IF (type_absorb(i) == ip_o2) THEN
            iump_o2=i
         ELSE IF (type_absorb(i) == ip_cfc113) THEN
            iump_cfc113=i
         ELSE IF (type_absorb(i) == ip_cfc114) THEN
            iump_cfc114=i
         ELSE IF (type_absorb(i) == ip_hcfc22) THEN
            iump_hcfc22=i
         ELSE IF (type_absorb(i) == ip_hfc125) THEN
            iump_hfc125=i
         ELSE IF (type_absorb(i) == ip_hfc134a) THEN
            iump_hfc134a=i
         END IF

      END DO


      IF (l_extra_top) THEN
!       The second radiative layer will be the first to have properties
!       set by copying input fields.
        i_top_copy=2
      ELSE
!       The first radiative layer will be the first to have properties
!       set by copying input fields.
        i_top_copy=1
      END IF


!     ASSIGN MIXING RATIOS OF THE GASES TO THE MAIN ARRAYS.

!     WATER VAPOUR:

      IF (iump_h2o >  0) THEN
!        No water exists above the wet levels.
         DO i=1, n_layer-nwet
            DO l=1, n_profile
               gas_mix_ratio(l, i, iump_h2o)=0.0e+00
            END DO
         END DO
         DO i=n_layer-nwet+1, n_layer
            DO l=1, n_profile
               lg=i_gather(l)
               gas_mix_ratio(l, i, iump_h2o)                            &
                 =MAX(h2o(lg, n_layer-i+1), 0.0e+00)
            END DO
         END DO
      END IF

!     CARBON DIOXIDE:

      IF (iump_co2 >  0) THEN
         DO i=1, n_layer
           IF (l_co2_3d) THEN
             DO l=1, n_profile
               lg=i_gather(l)
               gas_mix_ratio(l, i, iump_co2)=co2_3d(lg, n_layer-i+1)
             END DO
           ELSE
             DO l=1, n_profile
               gas_mix_ratio(l, i, iump_co2)=co2
             END DO
           END IF
         END DO
      END IF

!     OZONE:

      IF (iump_o3 >  0) THEN
!        The input field of ozone is supplied on NOZONE levels.
!        These values apply to the upper layers used by the full UM.
!        If NOZONE is smaller than NLEVS, the mixing ratio on the
!        bottom level supplied is copied to lower levels. If an
!        extra top level is used its mixing ratio is set by copying
!        the value for the top non-radiative level.
         IF (l_extra_top) THEN
           DO l=1, n_profile
             lg=i_gather(l)
             gas_mix_ratio(l, 1, iump_o3)=o3(lg, nozone)
           END DO
         END IF
         DO i=i_top_copy, nozone+i_top_copy-1
            DO l=1, n_profile
               lg=i_gather(l)
               gas_mix_ratio(l, i, iump_o3)=o3(lg, nozone+i_top_copy-i)
            END DO
         END DO
         DO i=nozone+i_top_copy, n_layer
            DO l=1, n_profile
               lg=i_gather(l)
               gas_mix_ratio(l, i, iump_o3)=o3(lg, 1)
            END DO
         END DO
         gas_mix_ratio(1:n_profile, 1:n_layer, iump_o3) =               &
            MAX(gas_mix_ratio(1:n_profile, 1:n_layer, iump_o3), 0.0)
      END IF



!     OTHER TRACE GASES:

!     THESE GASES ARE NOT ALWAYS INCLUDED IN THE CALCULATION.
!     TESTING IS THEREFORE MORE INTRICATE.

      IF (iump_n2o >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_n2o) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  IF (ngrgas >= p_n2o) THEN
                    lg=i_gather(l)
                    gas_mix_ratio(l, i, iump_n2o)=                      &
                          grgas_field(lg, n_layer-i+1,p_n2o)
                  ELSE
                    gas_mix_ratio(l, i, iump_n2o)=n2o_mix_ratio
                  END IF
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_n2o)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_n2o) THEN
            cmessage =                                                  &
              '*** ERROR: NITROUS OXIDE IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF

      IF (iump_ch4 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_ch4) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  IF (ngrgas >= p_ch4) THEN
                    lg=i_gather(l)
                    gas_mix_ratio(l, i, iump_ch4)=                      &
                          grgas_field(lg, n_layer-i+1,p_ch4)
                  ELSE
                    gas_mix_ratio(l, i, iump_ch4)=ch4_mix_ratio
                  END IF
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_ch4)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_ch4) THEN
            cmessage =                                                  &
              '*** ERROR: METHANE IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF

      IF (iump_cfc11 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_cfc11) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  IF (ngrgas >= p_f11) THEN
                    lg=i_gather(l)
                    gas_mix_ratio(l, i, iump_cfc11)=                    &
                          grgas_field(lg, n_layer-i+1,p_f11)
                  ELSE
                    gas_mix_ratio(l, i, iump_cfc11)=c11_mix_ratio
                  END IF
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc11)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_cfc11) THEN
            cmessage =                                                  &
              '*** ERROR: CFC11 IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF

      IF (iump_cfc12 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_cfc12) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  IF (ngrgas >= p_f12) THEN
                    lg=i_gather(l)
                    gas_mix_ratio(l, i, iump_cfc12)=                    &
                          grgas_field(lg, n_layer-i+1,p_f12)
                  ELSE
                    gas_mix_ratio(l, i, iump_cfc12)=c12_mix_ratio
                  END IF
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc12)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_cfc12) THEN
            cmessage =                                                  &
              '*** ERROR: CFC12 IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF

      IF (iump_o2 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_o2) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_o2)=o2_mix_ratio
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_o2)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_o2) THEN
            cmessage =                                                  &
              '*** ERROR: O2 IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF

      IF (iump_cfc113 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_cfc113) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  IF (ngrgas >= p_f113) THEN
                    lg=i_gather(l)
                    gas_mix_ratio(l, i, iump_cfc113)=                   &
                          grgas_field(lg, n_layer-i+1,p_f113)
                  ELSE
                    gas_mix_ratio(l, i, iump_cfc113)=c113_mix_ratio
                  END IF
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc113)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_cfc113) THEN
            cmessage =                                                  &
              '*** ERROR: CFC113 IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF

      IF (iump_cfc114 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_cfc114) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc114)=c114_mix_ratio
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc114)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_cfc114) THEN
            cmessage =                                                  &
              '*** ERROR: CFC114 IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF

      IF (iump_hcfc22 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_hcfc22) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  IF (ngrgas >= p_f22) THEN
                    lg=i_gather(l)
                    gas_mix_ratio(l, i, iump_hcfc22)=                   &
                          grgas_field(lg, n_layer-i+1,p_f22)
                  ELSE
                    gas_mix_ratio(l, i, iump_hcfc22)=hcfc22_mix_ratio
                  END IF
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_hcfc22)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_hcfc22) THEN
            cmessage =                                                  &
              '*** ERROR: HCFC22 IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF

      IF (iump_hfc125 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_hfc125) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_hfc125)=hfc125_mix_ratio
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_hfc125)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_hfc125) THEN
            cmessage =                                                  &
              '*** ERROR: HFC125 IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF

      IF (iump_hfc134a >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (l_hfc134a) THEN
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_hfc134a)=hfc134a_mix_ratio
               END DO
            END DO
         ELSE
            DO i=1, n_layer
               DO l=1, n_profile
                  gas_mix_ratio(l, i, iump_hfc134a)=0.0e+00
               END DO
            END DO
         END IF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (l_hfc134a) THEN
            cmessage =                                                  &
              '*** ERROR: HFC134A IS NOT IN THE SPECTRAL FILE.'
            ierr=i_err_fatal
            GO TO 9999
         END IF
      END IF


 9999 CONTINUE
! Check error condition
      IF (ierr /= i_normal) THEN

        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook('R2_SET_GAS_MIX_RATIO',zhook_out,zhook_handle)
      END SUBROUTINE r2_set_gas_mix_ratio
