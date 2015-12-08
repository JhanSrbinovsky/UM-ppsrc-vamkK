! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Returns information about a given grid type

      SUBROUTINE GT_DECODE(                                             &
     &                      GRID_TYPE,                                  &
     &                      MODEL_TYPE,CONTENT,COVERAGE,DOMAIN,CYCLIC)

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE ereport_mod, ONLY : ereport
      USE grdtypes_mod, ONLY: gt_unset, gt_atmos, gt_ocean, gt_wave,    &
                              gt_thetamass, gt_velocity, gt_u_c, gt_v_c,&
                              gt_hybrid, gt_river, gt_allpts, gt_land,  &
                              gt_sea, gt_full, gt_zonal, gt_meridional, &
                              gt_ozone, gt_scalar, gt_compressed,       &
                              gt_lbc, gt_nocyclic, gt_optcyclic, gt_cyclic
      IMPLICIT NONE

! Given an input GRIDTYPE, this routine will return a value
! for each of MODEL_TYPE, CONTENT, COVERAGE, DOMAIN and CYCLIC
! from the gt_* variables defined in the GRDTYPES comdeck


! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Arguments:

      INTEGER                                                           &
     &  GRID_TYPE                                                       &
                   ! IN  : Grid type to investigate
     &, MODEL_TYPE                                                      &
                   ! OUT : What model type grid is for
     &, CONTENT                                                         &
                   ! OUT : What type of data the grid is for
     &, COVERAGE                                                        &
                   ! OUT : What type of points are on the grid
     &, DOMAIN                                                          &
                   ! OUT : What subset of points are on the grid
     &, CYCLIC     ! OUT : If the grid includes extra cyclic wrap
                   !       around points at the start and end of
                   !       each row

! Comdecks

! Local variables

      INTEGER max_grid_types
      PARAMETER (max_grid_types=65)

      INTEGER GRID_DATA(5,max_grid_types)
!             Array holding all the descriptions of the different
!             grid types.5 descriptors for each grid type.

      INTEGER                                                           &
     &  ICODE      ! Error code

      CHARACTER(LEN=80)                                                      &
     &  CMESSAGE   ! Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Definition of the field types
      DATA                                                              &
     &  GRID_DATA(1:5,1)                                                &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_full,gt_nocyclic/          &
     &, GRID_DATA(1:5,2)                                                &
     &   /gt_atmos,gt_thetamass,gt_land,gt_full,gt_nocyclic/            &
     &, GRID_DATA(1:5,3)                                                &
     &   /gt_atmos,gt_thetamass,gt_sea,gt_full,gt_nocyclic/             &
     &, GRID_DATA(1:5,4)                                                &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_zonal,gt_nocyclic/         &
     &, GRID_DATA(1:5,5)                                                &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_meridional,gt_nocyclic/    &
     &, GRID_DATA(1:5,6)  /5*gt_unset/                                  &
     &, GRID_DATA(1:5,7)  /5*gt_unset/                                  &
     &, GRID_DATA(1:5,8)  /5*gt_unset/                                  &
     &, GRID_DATA(1:5,9)  /5*gt_unset/                                  &
     &, GRID_DATA(1:5,10) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,11)                                               &
     &   /gt_atmos,gt_velocity,gt_allpts,gt_full,gt_nocyclic/           &
     &, GRID_DATA(1:5,12)                                               &
     &   /gt_atmos,gt_velocity,gt_land,gt_full,gt_nocyclic/             &
     &, GRID_DATA(1:5,13)                                               &
     &   /gt_atmos,gt_velocity,gt_sea,gt_full,gt_nocyclic/              &
     &, GRID_DATA(1:5,14)                                               &
     &   /gt_atmos,gt_velocity,gt_allpts,gt_zonal,gt_nocyclic/          &
     &, GRID_DATA(1:5,15)                                               &
     &   /gt_atmos,gt_velocity,gt_allpts,gt_meridional,gt_nocyclic/     &
     &, GRID_DATA(1:5,16) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,17)                                               &
     &   /gt_atmos,gt_hybrid,gt_allpts,gt_scalar,gt_nocyclic/           &
     &, GRID_DATA(1:5,18)                                               &
     &   /gt_atmos,gt_U_C,gt_allpts,gt_full,gt_nocyclic/                &
     &, GRID_DATA(1:5,19)                                               &
     &   /gt_atmos,gt_V_C,gt_allpts,gt_full,gt_nocyclic/                &
     &, GRID_DATA(1:5,20) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,21)                                               &
     &   /gt_atmos,gt_thetamass,gt_land,gt_compressed,gt_nocyclic/      &
     &, GRID_DATA(1:5,22)                                               &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_ozone,gt_nocyclic/         &
     &, GRID_DATA(1:5,23)                                               &
     &   /gt_atmos,gt_river,gt_allpts,gt_full,gt_nocyclic/              &
     &, GRID_DATA(1:5,24) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,25)                                               &
     &   /gt_atmos,gt_hybrid,gt_allpts,gt_LBC,gt_nocyclic/              &
     &, GRID_DATA(1:5,26)                                               &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_LBC,gt_nocyclic/           &
     &, GRID_DATA(1:5,27)                                               &
     &   /gt_atmos,gt_U_C,gt_allpts,gt_LBC,gt_nocyclic/                 &
     &, GRID_DATA(1:5,28)                                               &
     &   /gt_atmos,gt_V_C,gt_allpts,gt_LBC,gt_nocyclic/                 &
     &, GRID_DATA(1:5,29)                                               &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_LBC,gt_nocyclic/           &
     &, GRID_DATA(1:5,30) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,31)                                               &
     &   /gt_ocean,gt_thetamass,gt_sea,gt_compressed,gt_nocyclic/       &
     &, GRID_DATA(1:5,32)                                               &
     &   /gt_ocean,gt_velocity,gt_sea,gt_compressed,gt_nocyclic/        &
     &, GRID_DATA(1:5,33) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,34) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,35) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,36)                                               &
     &   /gt_ocean,gt_thetamass,gt_allpts,gt_full,gt_optcyclic/         &
     &, GRID_DATA(1:5,37)                                               &
     &   /gt_ocean,gt_velocity,gt_allpts,gt_full,gt_optcyclic/          &
     &, GRID_DATA(1:5,38)                                               &
     &   /gt_ocean,gt_U_C,gt_allpts,gt_full,gt_optcyclic/               &
     &, GRID_DATA(1:5,39)                                               &
     &   /gt_ocean,gt_V_C,gt_allpts,gt_full,gt_optcyclic/               &
     &, GRID_DATA(1:5,40) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,41)                                               &
     &   /gt_ocean,gt_thetamass,gt_allpts,gt_full,gt_cyclic/            &
     &, GRID_DATA(1:5,42)                                               &
     &   /gt_ocean,gt_velocity,gt_allpts,gt_full,gt_cyclic/             &
     &, GRID_DATA(1:5,43)                                               &
     &   /gt_ocean,gt_thetamass,gt_allpts,gt_zonal,gt_nocyclic/         &
     &, GRID_DATA(1:5,44)                                               &
     &   /gt_ocean,gt_velocity,gt_allpts,gt_zonal,gt_nocyclic/          &
     &, GRID_DATA(1:5,45)                                               &
     &   /gt_ocean,gt_thetamass,gt_allpts,gt_meridional,gt_nocyclic/    &
     &, GRID_DATA(1:5,46)                                               &
     &   /gt_ocean,gt_velocity,gt_allpts,gt_meridional,gt_nocyclic/     &
     &, GRID_DATA(1:5,47)                                               &
     &   /gt_ocean,gt_hybrid,gt_allpts,gt_scalar,gt_nocyclic/           &
     &, GRID_DATA(1:5,48) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,49) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,50) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,51)                                               &
     &   /gt_ocean,gt_hybrid,gt_allpts,gt_LBC,gt_nocyclic/              &
     &, GRID_DATA(1:5,52) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,53) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,54) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,55) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,56) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,57) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,58) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,59) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,60)                                               &
     &   /gt_wave,gt_thetamass,gt_allpts,gt_full,gt_nocyclic/           &
     &, GRID_DATA(1:5,61) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,62)                                               &
     &   /gt_wave,gt_thetamass,gt_sea,gt_compressed,gt_nocyclic/        &
     &, GRID_DATA(1:5,63) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,64) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,65)                                               &
     &   /gt_wave,gt_hybrid,gt_allpts,gt_LBC,gt_nocyclic/

! And the code

      IF (lhook) CALL dr_hook('GT_DECODE',zhook_in,zhook_handle)
      IF ((GRID_TYPE  <   1) .OR.                                       &
     &    (GRID_TYPE  >   max_grid_types)) THEN

        ICODE=1
        WRITE(CMESSAGE,*) 'Invalid GRID_TYPE ',GRID_TYPE,               &
     &                    ' passed to GT_DECODE'

        CALL EREPORT('GT_DECODE',ICODE,CMESSAGE)

      ENDIF

      MODEL_TYPE = GRID_DATA(1,GRID_TYPE)
      CONTENT    = GRID_DATA(2,GRID_TYPE)
      COVERAGE   = GRID_DATA(3,GRID_TYPE)
      DOMAIN     = GRID_DATA(4,GRID_TYPE)
      CYCLIC     = GRID_DATA(5,GRID_TYPE)

      IF (MODEL_TYPE  ==  gt_unset) THEN

        ICODE=2
        WRITE(CMESSAGE,*) 'GRID_TYPE ',GRID_TYPE,                       &
     &                    ' undefined in GT_DECODE'
        CALL EREPORT('GT_DECODE',ICODE,CMESSAGE)

      ENDIF

      IF (lhook) CALL dr_hook('GT_DECODE',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE GT_DECODE
