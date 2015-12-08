! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initialises variables of multi-layer snow scheme.

MODULE rcf_init_ml_snow_mod
IMPLICIT NONE

! Description:
!     This subroutine initialises all multi-layer snow scheme variables,
!     as a preconditioning step for the main model run. It combines the
!     functions of total_snow and layersnow in JULES standalone.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE rcf_init_ml_snow( fields_out, field_count_out, hdr_out, &
                             ml_snow_field, field_stashcode )

USE ereport_mod, ONLY : &
    ereport

USE rcf_locate_mod, ONLY : &
    rcf_locate

USE rcf_grid_type_mod, ONLY : &
    grid_type

USE rcf_alloc_field_mod, ONLY : &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE rcf_stashcodes_mod, ONLY :    &
    stashcode_prog_sec          , &
    stashcode_soil_temp         , &
    stashcode_snow_tile         , &
    stashcode_snow_grnd         , &
    stashcode_snowdep_grd_tile  , &
    stashcode_snowpack_bk_dens  , &
    stashcode_nsnow_layrs_tiles , &
    stashcode_snow_laythk_tiles , &
    stashcode_snow_ice_tile     , &
    stashcode_snow_liq_tile     , &
    stashcode_snow_t_tile       , &
    stashcode_snow_laydns_tiles

USE rcf_field_type_mod, ONLY : &
    field_type

USE rcf_umhead_mod, ONLY : &
    um_header_type

USE printstatus_mod, ONLY : &
    printstatus,            &
    prstatus_normal

USE um_parvars, ONLY : &
    mype

USE decomp_params, ONLY : &
    decomp_rcf_output

USE rcf_read_field_mod, ONLY : &
    rcf_read_field

USE rcf_write_field_mod, ONLY : &
    rcf_write_field

USE rcf_lsm_mod, ONLY : &
    local_land_out

USE water_constants_mod, ONLY:  &
    tm

USE nlsizes_namelist_mod, ONLY : &
    ntiles

USE ancil_info, ONLY:  &
    nsmax

USE snow_param, ONLY : &
    rho_snow_const,    &
    rho_snow_fresh,    &
    cansnowtile,       &
    dzsnow

USE switches, ONLY: &
  can_model

USE rcf_recon_mod, ONLY : &
    l_canopy_snow_throughfall

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( field_type ), INTENT(INOUT) :: ml_snow_field
TYPE( um_header_type), INTENT(IN) :: hdr_out
INTEGER, INTENT(IN)               :: field_count_out
INTEGER, INTENT(IN)               :: field_stashcode

! Internal variables
TYPE( field_type ), POINTER  :: soil_temp
TYPE( field_type ), POINTER  :: snow_tile
TYPE( field_type ), POINTER  :: snow_grnd

! These fields will be saved
REAL, ALLOCATABLE, SAVE  :: snowdepth  (:,:)
REAL, ALLOCATABLE, SAVE  :: bulk_dens  (:,:)
REAL, ALLOCATABLE, SAVE  :: nsnow      (:,:)
REAL, ALLOCATABLE, SAVE  :: ds         (:,:)
REAL, ALLOCATABLE, SAVE  :: snow_ice   (:,:)
REAL, ALLOCATABLE, SAVE  :: snow_liq   (:,:)
REAL, ALLOCATABLE, SAVE  :: tsnow      (:,:)
REAL, ALLOCATABLE, SAVE  :: layer_dens (:,:)

REAL         :: total_snow_ground( local_land_out, ntiles )
REAL         :: snow_ground_store( local_land_out, ntiles )

REAL         ::  remains     ! Remaining depth of snow
                             ! for other layers.
INTEGER      ::  errorstatus
INTEGER      ::  pos         ! position in array
INTEGER      ::  i,j,n       ! loop counters
INTEGER      ::  pseudo      ! position in pseudo-level

CHARACTER (LEN=*), PARAMETER      :: RoutineName='Rcf_Init_ML_Snow'
CHARACTER (LEN=80)                :: Cmessage

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
    IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
      WRITE (6,'(a,i4)') ' Initialising multi-layer snow variable: Stashcode' &
                         ,field_stashcode
    END IF

!----------------------------------------------------------------------
! Only calculate when no saved fields are allocated.
!----------------------------------------------------------------------
IF (.NOT.(     ALLOCATED(snowdepth ) &
           .OR.ALLOCATED(bulk_dens ) &
           .OR.ALLOCATED(nsnow     ) &
           .OR.ALLOCATED(ds        ) &
           .OR.ALLOCATED(snow_ice  ) &
           .OR.ALLOCATED(snow_liq  ) &
           .OR.ALLOCATED(tsnow     ) &
           .OR.ALLOCATED(layer_dens)  ) ) THEN
  
!----------------------------------------------------------------------
! Find required fields in output dump and read them in
!----------------------------------------------------------------------
! soil temperatures
CALL rcf_locate( stashcode_prog_sec, stashcode_soil_temp,                   &
           fields_out, field_count_out, pos)
soil_temp => fields_out(pos)
CALL rcf_alloc_field( soil_temp )
CALL rcf_read_field( soil_temp, hdr_out, decomp_rcf_output )

! snow on tile
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_tile,                   &
           fields_out, field_count_out, pos)
snow_tile => fields_out(pos)
CALL rcf_alloc_field( snow_tile )
CALL rcf_read_field( snow_tile, hdr_out, decomp_rcf_output )

IF (can_model == 4) THEN
! snow under canopy
  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_grnd,                 &
             fields_out, field_count_out, pos)
  snow_grnd => fields_out(pos)
  CALL rcf_alloc_field( snow_grnd )
  CALL rcf_read_field( snow_grnd, hdr_out, decomp_rcf_output )
END IF

!-------------------------------------------------------------------------------
! Do sanity checks because these fields are used interchangeably.
!-------------------------------------------------------------------------------
IF ( snow_tile % level_size /= local_land_out ) THEN
          errorstatus = 43
          WRITE (cmessage, '(A, I8, A, I8)')                &
                'Field sizes do not agree: level_size: ',   &
                snow_tile % level_size,                     &
                ' local_land_out: ',                            &
                local_land_out
          CALL ereport( routinename, errorstatus, cmessage )
END IF
IF ( snow_tile % levels /= ntiles ) THEN
          errorstatus = 43
          WRITE (cmessage, '(A, I3, A, I3)')                &
                'Field sizes do not agree: levels: ',       &
                snow_tile % levels,                         &
                ' ntiles: ',                                &
                ntiles
          CALL ereport( routinename, errorstatus, cmessage )
END IF

!-------------------------------------------------------------------------------
! Allocate fields to save.
!-------------------------------------------------------------------------------
ALLOCATE(snowdepth (local_land_out, ntiles) )
ALLOCATE(bulk_dens (local_land_out, ntiles) )
ALLOCATE(nsnow     (local_land_out, ntiles) )
IF ( nsmax > 0 ) THEN
  ALLOCATE(ds        (local_land_out, ntiles * nsmax) )
  ALLOCATE(snow_ice  (local_land_out, ntiles * nsmax) )
  ALLOCATE(snow_liq  (local_land_out, ntiles * nsmax) )
  ALLOCATE(tsnow     (local_land_out, ntiles * nsmax) )
  ALLOCATE(layer_dens(local_land_out, ntiles * nsmax) )
END IF

!-------------------------------------------------------------------------------
! Record the snow amount on the ground,
! putting all snow onto the ground and zeroing canopy snow if required.
!-------------------------------------------------------------------------------
 IF (can_model == 4) THEN
   snow_ground_store = snow_grnd % data
 ELSE
   snow_ground_store = 0.0
 END IF

 IF (l_canopy_snow_throughfall) THEN

  IF (mype == 0 .AND. printstatus >= prstatus_normal) THEN
    WRITE (6,'(a)') ' Moving snow from canopy to ground'
  END IF

  DO n=1,ntiles
    IF ( cansnowtile(n).AND.(can_model == 4) ) THEN
      DO i=1,local_land_out
        total_snow_ground(i,n) = MAX(  snow_tile % data(i,n), 0.0 ) &
                               + MAX( snow_ground_store(i,n), 0.0 )
      END DO
    ELSE
      DO i=1,local_land_out
        total_snow_ground(i,n) = MAX(  snow_tile % data(i,n), 0.0 )
      END DO
    END IF
  END DO
! Initialise stores to zero.
  snow_ground_store = 0.0
  snow_tile % data  = 0.0

  DO n=1,ntiles
    IF ( canSnowTile(n).AND.(can_model == 4) ) THEN
      DO i=1,local_land_out
        snow_ground_store(i,n) = total_snow_ground(i,n)
      END DO
    ELSE
      DO i=1,local_land_out
         snow_tile % data(i,n) = total_snow_ground(i,n)
      END DO
    END IF
  END DO

 ELSE ! no throughfall

  DO n=1,ntiles
    IF ( cansnowtile(n).AND.(can_model == 4) ) THEN
      DO i=1,local_land_out
        total_snow_ground(i,n) = MAX( snow_ground_store(i,n), 0.0 )
      END DO
    ELSE
      DO i=1,local_land_out
        total_snow_ground(i,n) = MAX(  snow_tile % data(i,n), 0.0 )
      END DO
    END IF
  END DO

 END IF ! throughfall test

!-------------------------------------------------------------------------------
! Set snowpack bulk density and calculate snow depth from mass per unit area.
!-------------------------------------------------------------------------------
! Initialise with values that will be retained where there is
! no tile fraction - using "sensible" values for when these are printed.
!
  snowdepth = 0.0

  IF ( nsmax < 1 ) THEN
    bulk_dens = rho_snow_const
  ELSE
    bulk_dens = rho_snow_fresh
  END IF

  DO n=1,ntiles
    DO i=1,local_land_out
!     Use the constant (snowpack) density for nsmax=0 and if there is an
!     existing pack. If nsmax>0 and there is no pack, initialise the density
!     to the fresh snow value so that this value is used when/if a snowpack
!     next develops.
      IF ( nsmax < 1 .OR.  &
           ( total_snow_ground(i,n) > EPSILON(total_snow_ground) ) ) THEN
        bulk_dens(i,n) = rho_snow_const
      ELSE
        bulk_dens(i,n) = rho_snow_fresh
      END IF
      snowdepth(i,n) = MAX(total_snow_ground(i,n),0.0) &
                       / bulk_dens(i,n)
    END DO
  END DO

!----------------------------------------------------------------------
! Calculate snow layer thicknesses.
! Set layer frozen and liquid content, density, temperature.
!----------------------------------------------------------------------
IF ( nsmax > 0 ) THEN

! Initialise layer numbers
  nsnow = 0.0

! Initialise layer depths
! (this value will persist at locations where tile frac=0).
  ds = 0.0

  DO j = 1,ntiles
    DO i = 1,local_land_out

!   Only divide into layers if depth is >= a threshold.
      IF ( snowdepth(i,j) >= dzsnow(1) ) THEN
        remains = snowdepth(i,j)

        DO n=1,nsmax

! pseudo-level index for NTILES*NSMAX dimension
          pseudo = j + (n-1) * ntiles

          ds(i,pseudo) = dzsnow(n)
          remains = remains - dzsnow(n)
          IF ( remains <= dzsnow(n) .OR. n == nsmax ) THEN
            ds(i,pseudo) = ds(i,pseudo) + remains
            EXIT
          END IF
        END DO

! Set number of layers.
        nsnow(i,j) = REAL(n)
      END IF    !  >dzSnow(1)

    END DO
  END DO

  snow_ice   = 0.0
  snow_liq   = 0.0
  layer_dens = 0.0
  tsnow      = tm

! layer ice-mass and density
  DO j = 1,ntiles
    DO i = 1,local_land_out
      DO n=1,NINT(nsnow(i,j))

! pseudo-level index for NTILES*NSMAX dimension
        pseudo = j + (n-1) * ntiles

        IF (snowdepth(i,j) > 0.0) THEN
          snow_ice(i,pseudo) = &
           total_snow_ground(i,j) * ds(i,pseudo) / snowdepth(i,j)
          layer_dens(i,pseudo) = bulk_dens(i,j)
        END IF

      END DO
    END DO
  END DO

! layer temperatures
  DO i=1,local_land_out
    tsnow(i,:) = MIN( soil_temp % data(i,1), tm )
  END DO

END IF ! nsmax

!----------------------------------------------------------------------
! Write out the calculated snow fields.
!----------------------------------------------------------------------
IF (l_canopy_snow_throughfall) THEN
  CALL rcf_write_field( snow_tile , hdr_out, decomp_rcf_output )
  IF (can_model == 4) THEN
    snow_grnd % data = snow_ground_store
    CALL rcf_write_field( snow_grnd , hdr_out, decomp_rcf_output )
  END IF
END IF

!----------------------------------------------------------------------
! Clear up dynamic memory used.
!----------------------------------------------------------------------
CALL rcf_dealloc_field( soil_temp )
CALL rcf_dealloc_field( snow_tile )
IF (can_model == 4) THEN
  CALL rcf_dealloc_field( snow_grnd )
END IF

END IF ! None of saved fields allocated

!----------------------------------------------------------------------
! Initialise output field and deallocate saved field.
!----------------------------------------------------------------------
SELECT CASE( field_stashcode )

! snow depth
  CASE( stashcode_snowdep_grd_tile  )
    ml_snow_field % data = snowdepth
    DEALLOCATE(snowdepth )

! snowpack bulk density
  CASE( stashcode_snowpack_bk_dens  )
    ml_snow_field % data = bulk_dens
    DEALLOCATE(bulk_dens )

! number of snow layers on each tile
  CASE( stashcode_nsnow_layrs_tiles )
    ml_snow_field % data = nsnow
    DEALLOCATE(nsnow     )

! snow layer thicknesses on each tile
  CASE( stashcode_snow_laythk_tiles )
    ml_snow_field % data = ds
    DEALLOCATE(ds        )

! solid-phase snow amount in each layer on each tile
  CASE( stashcode_snow_ice_tile     )
    ml_snow_field % data = snow_ice
    DEALLOCATE(snow_ice  )

! liquid-phase snow amount in each layer on each tile
  CASE( stashcode_snow_liq_tile     )
    ml_snow_field % data = snow_liq
    DEALLOCATE(snow_liq  )

! temperature of snow layers on each tile
  CASE( stashcode_snow_t_tile       )
    ml_snow_field % data = tsnow
    DEALLOCATE(tsnow     )

! density of snow layers on each tile
  CASE( stashcode_snow_laydns_tiles )
    ml_snow_field % data = layer_dens
    DEALLOCATE(layer_dens)

END SELECT


RETURN
END SUBROUTINE rcf_init_ml_snow
END MODULE rcf_init_ml_snow_mod
