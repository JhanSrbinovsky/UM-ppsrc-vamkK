! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Set up a UM header describing the GRIB data

MODULE Rcf_Grib_Sethdr_Mod

! SUBROUTINE Rcf_Grib_SetHdr

! Description: Set up the UM style header using the GRIB header
!              information.

! Method: Copies information from the GRIB headers into the various
!         derived types used for UM headers and then calls the normal
!         header set up routines.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CONTAINS
SUBROUTINE Rcf_Grib_Sethdr(lists,output_grid,hdr_dmy,hdr_itm)

USE Rcf_Grib_Block_Params_Mod
! Has derived types of List_Marker and Grib_record
! plus also containd most parameter definitions starting with p_

USE Rcf_Grib_Lookups_Mod, ONLY :                                               &
  grib_max_fields,                                                             &
  grib_u_field,                                                                &
  grib_soil_temp_field,                                                        &
  grib_soil_moist_field

USE Rcf_Grid_Type_Mod, ONLY :                                                  &
  grid_type

USE Rcf_Umhead_Mod, ONLY :                                                     &
  lenfixhd,                                                                    &
  um_header_type            ! Derived containing UM header info

USE Rcf_Headaddress_Mod

USE Rcf_Headers_Mod, ONLY :                                                    &
  fixhd

USE PrintStatus_Mod, ONLY : &
  printstatus,              &
  prstatus_min,             & ! =1 Minimum output
  prstatus_normal,          & ! =2 Short informative output
  prstatus_oper,            & ! =3 Full informative output
  prstatus_diag               ! =4 Extra Diagnostic output

USE UM_ParVars, ONLY :                                                    &
  mype

USE Submodel_Mod, ONLY :                                                       &
  atmos_im

USE Ereport_Mod, ONLY :                                                        &
  Ereport

USE Rcf_Allochdr_Mod, ONLY :                                                   &
  Rcf_Allochdr

USE Rcf_Setup_Fixhd_Mod, ONLY :                                                &
  Rcf_Setup_Fixhd

USE Rcf_Generate_Heights_Mod, ONLY :                                           &
  height_gen_ecmwf_press,                                                      &
  height_gen_ecmwf_hybrd

USE lookup_addresses

IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(in):>
TYPE (list_marker)               :: lists(0:grib_max_fields)

!< Scalar arguments with intent(InOut):>
TYPE (um_header_type)                    :: hdr_dmy
TYPE (um_header_type)                    :: hdr_itm
TYPE (grid_type)                         :: output_grid

! Comdecks
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Local variables

TYPE (grib_record), POINTER      :: current
TYPE (grib_record), POINTER      :: first_multi  ! first multi level
! list for hdr info

INTEGER                          :: namelst_temp(lenfixhd)

CHARACTER (len=*), PARAMETER     :: routinename='Rcf_GRIB_SetHdr'
CHARACTER (len=80)               :: cmessage(2)   ! used for EReport
INTEGER                          :: errorstatus   ! used for EReport
INTEGER                          :: i,count
INTEGER                          :: sign

INTEGER                          :: fields_tot
INTEGER                          :: int_val  ! Used for mold in transfer
INTEGER                          :: first_multi_list
INTEGER                          :: first_multi_type

!=======================================================================
! Preparation for setting up the fixed header.
!=======================================================================

hdr_dmy % fixhd (:)     = IMDI

! Count the number of records in all the lists _except_ the misc list
! Also, set pointer to the first entry in the first pressure field
! list encountered.

fields_tot       = 0
first_multi_list = 0
NULLIFY (first_multi)

DO i = 1, grib_max_fields
  IF (ASSOCIATED(lists(i) % begin) ) THEN
    fields_tot = fields_tot + lists(i) % lstcount

    IF ( (.NOT. ASSOCIATED(first_multi)) .AND.  (                              &
      (lists(i) % begin % block_1(p_lvl_type) == tb3_pressure)  .OR.           &
      (lists(i) % begin % block_1(p_lvl_type) == tb3_hybrid  )                 &
      )) THEN
      first_multi      => lists(i) % begin
      first_multi_list = i
      first_multi_type = lists(i) % begin % block_1(p_lvl_type)
    END IF
  END IF
END DO

! Double check that 'First_Multi' was assigned
IF (.NOT. ASSOCIATED(first_multi)) THEN
  cmessage(1) = "No multi-level fields found in GRIB data"
  errorstatus = 10
  CALL Ereport( routinename, errorstatus, cmessage(1))
END IF

! Set Values that would normally have been taken from the start dump
output_grid % glob_p_row_length = first_multi % block_2(p_pnts_prll)
output_grid % glob_p_rows       = first_multi % block_2(p_pnts_merid)
output_grid % model_levels      = lists(first_multi_list) % lstcount

! Check for Soil Temp Levels before using their count.
IF (ASSOCIATED(lists(grib_soil_temp_field) % begin)) THEN
  output_grid % st_levels     = lists(grib_soil_temp_field) % lstcount
END IF

! Check for Soil Moisture Levels before using their count.
IF (ASSOCIATED(lists(grib_soil_moist_field) % begin)) THEN
  output_grid % sm_levels     = lists(grib_soil_moist_field) % lstcount
END IF

hdr_dmy % fixhd (fh_runid)      = 0
hdr_dmy % fixhd (fh_exptno)     = IMDI
hdr_dmy % fixhd (fh_submodel)   = atmos_im

hdr_itm % lenintc               = p_len_intc
hdr_itm % lenrealc              = p_len_realc
hdr_itm % len1levdepc           = lists(grib_u_field) % lstcount + 1
hdr_itm % len2levdepc           = p_len2_levdepc
hdr_itm % len1rowdepc           = p_len1_rowdepc
hdr_itm % len2rowdepc           = p_len2_rowdepc
hdr_itm % len1coldepc           = p_len1_coldepc
hdr_itm % len2coldepc           = p_len2_coldepc
hdr_itm % len1fldsofc           = p_len1_fldsofc
hdr_itm % len2fldsofc           = p_len2_fldsofc
hdr_itm % lenextrac             = p_len_extrac
hdr_itm % lenhistfile           = p_len_histfile
hdr_itm % lencompfldi1          = p_len_compfldi1
hdr_itm % lencompfldi2          = p_len_compfldi2
hdr_itm % lencompfldi3          = p_len_compfldi3
hdr_itm % len1lookup            = p_len1_lookup
hdr_itm % len2lookup            = fields_tot
hdr_itm % lendata               = fields_tot *                                 &
  first_multi % block_2(p_pnts_prll) *                                         &
  first_multi % block_2(p_pnts_merid)

! Set up the Initial data Time and Validity Time
hdr_dmy % fixhd (fh_dtYear)   = first_multi % block_1 (p_ref_year)             &
  + ((first_multi % block_1 (p_ref_cent) -1) * 100)
hdr_dmy % fixhd (fh_dtMonth:fh_dtMinute) =                                     &
  first_multi % block_1 (p_ref_month:p_ref_min)
hdr_dmy % fixhd (fh_dtSecond) = 0
hdr_dmy % fixhd (fh_dtDayNo)  = 0

hdr_dmy % fixhd (fh_vtYear)   = first_multi % block_1 (p_ref_year)             &
  + ((first_multi % block_1 (p_ref_cent) -1) * 100)
hdr_dmy % fixhd (fh_vtMonth:fh_vtMinute) =                                     &
  first_multi % block_1 (p_ref_month:p_ref_min)
hdr_dmy % fixhd (fh_vtSecond) = 0
hdr_dmy % fixhd (fh_vtDayNo)  = 0

IF (printstatus >= prstatus_diag) THEN
  IF ( mype == 0 ) THEN
    WRITE (6,'(A,I4,2("/",I2),/,A,2(I2,":"),I2,/,A,I3)')                       &
      "Date from GRIB data appears to be (yyyy/mm/dd) : ",                     &
      hdr_dmy % fixhd (fh_dtYear:fh_dtDay),                                    &
      "And the Time is : ",hdr_dmy % fixhd (fh_dtHour:fh_dtSecond),            &
      "Dayno. is :",hdr_dmy % fixhd (fh_dtDayNo)
  END IF
END IF

! Set the 'global' flag if appropriate
! Check that the grid represents 360 by 180 degrees
!
! p_Pnts_Prll is the number of points along a parallel 
! p_IncrPrll is the increment along a parallel in millidegrees
! The product of these two should be close to 360000
! millidegrees (360 degrees) if it is a global grid
!
! p_Pnts_Merid is the number of points along a meridian 
! p_IncrMerid is the increment along a meridian in millidegrees
! The product of these two should be close to 180000
! millidegrees (180 degrees) if it is a global grid
! 
IF ( ( (first_multi % block_2(p_pnts_prll) *                                   &
  first_multi % block_2(p_incrprll )  ) >= 359999 ) .AND.                      &
  ( (first_multi % block_2(p_pnts_merid) *                                     &
  first_multi % block_2(p_incrmerid)   ) >= 179999 ) )  THEN

  output_grid % global = .TRUE.

ELSE
  ! Areas of the GRIB handling code _assume_ the incoming grid is a
  ! global grid and thus an abort is placed here
  ! Provide information on lat/lon grid to help debugging
  cmessage(1) = "Does not appear to be Global grid"
  WRITE (6,*) cmessage(1)
  WRITE (6,*) 'No of lons, longitude increment, product=',                   &
    first_multi % block_2(p_pnts_prll),first_multi % block_2(p_incrprll ),   &
    first_multi % block_2(p_pnts_prll)*first_multi % block_2(p_incrprll )
  WRITE (6,*) 'No of lats, latitude increment, product=',                    &
    first_multi % block_2(p_pnts_merid),first_multi % block_2(p_incrmerid ), &
    first_multi % block_2(p_pnts_merid)*first_multi % block_2(p_incrmerid )
  errorstatus = 20
  CALL ereport( routinename, errorstatus, cmessage(1))

END IF  ! test to see if using global grid

! Call to Set up fixed header.

! First - store Recon Namelist and send dummy blank one to prevent
!         overwriting of data
namelst_temp(:) = fixhd(:)
fixhd(:) = IMDI

CALL Rcf_Setup_Fixhd(hdr_dmy, hdr_itm)

! Reset values describing the grid which are 'hardwired' in above call
hdr_itm % fixhd ( fh_vertcoord )   = fh_vertcoord_pressure
hdr_itm % fixhd ( fh_gridstagger ) = fh_gridstagger_a

! Not forgetting to restore the Recon Namelist
fixhd(:) = namelst_temp(:)

! Call to Allocate memory for rest of header info
CALL Rcf_Allochdr (hdr_itm)

!=======================================================================
! Initialising Integer Constants
!=======================================================================
! Emulating what was done at 4.5

hdr_itm % intc (:) = IMDI

hdr_itm % intc (ic_xlen)          = output_grid % glob_p_row_length
hdr_itm % intc (ic_ylen)          = output_grid % glob_p_rows
hdr_itm % intc (ic_plevels)       = output_grid % model_levels
hdr_itm % intc (ic_wetlevels)     = output_grid % model_levels
hdr_itm % intc (ic_soiltlevels)   = output_grid % st_levels
hdr_itm % intc (ic_soilmoistlevs) = output_grid % sm_levels
! The following requires valid values since there are used to calculate level
! information within rcf_level_code.
hdr_itm % intc (ic_tracerlevs)    = output_grid % model_levels

! Why set number of boundary layer levels like this?
! Hardwired to pressure levels!
hdr_itm % intc (ic_blevels)     = 0
current => lists(first_multi_list) % begin  ! just to keep following
DO WHILE (ASSOCIATED(current))
  IF (current % block_1 (p_lvl_desc_1) > 850 ) THEN
    hdr_itm % intc (ic_blevels) = hdr_itm % intc (ic_blevels) +1
  END IF
  current => current % next
END DO

hdr_itm % intc (ic_mdi)         = IMDI

! Set height_generation_method for Pressure or Hybrid
SELECT CASE (first_multi_type)

CASE (tb3_pressure)
  hdr_itm % intc (ic_heightmethod)  = height_gen_ecmwf_press

CASE (tb3_hybrid  )
  hdr_itm % intc (ic_heightmethod)  = height_gen_ecmwf_hybrd

CASE default
  cmessage(1)    =                                                             &
    'Multi-level field type not set for height generation method'
  errorstatus = 55
  CALL Ereport( routinename, errorstatus, cmessage(1) )

END SELECT

!=======================================================================
! Initialising Real Constants
!=======================================================================

hdr_itm % realc (:) = rmdi

! First latitude and longitude points in Real Consts Header
hdr_itm % realc (rc_firstlat)    = first_multi % block_2(p_latgrdpnt1)         &
  / 1000.000
hdr_itm % realc (rc_firstlong)   = first_multi % block_2(p_longrdpnt1)         &
  / 1000.000

! Latitude and longitude spacing values for Real Const (normally +ve)
hdr_itm % realc (rc_latspacing) = first_multi % block_2(p_incrmerid)           &
  / 1000.000
hdr_itm % realc (rc_longspacing) = first_multi % block_2(p_incrprll)           &
  / 1000.000

! GRIB data holds the Lat and Long of the Southern Pole _iff_ the grid
! is rotated - otherwise both values are zero.
IF ( (first_multi % block_2(p_s_pole_lat) == 0.00 ) .AND.                      &
  (first_multi % block_2(p_s_pole_lon) == 0.00 ) ) THEN
  hdr_itm % realc (rc_polelat)     = 90.0
  hdr_itm % realc (rc_polelong)    = 00.0

  IF (printstatus >= prstatus_diag) THEN
    IF (mype == 0 ) THEN
      WRITE (6,*) "Using default value for North Pole lat and Long"
    END IF
  END IF
ELSE
  IF (printstatus >= prstatus_diag) THEN
    IF (mype == 0 ) THEN
      WRITE (6,*) "Using calculated North Pole lat and Long"
    END IF
  END IF

  hdr_itm % realc (rc_polelat) =                                               &
    90.0 - (first_multi % block_2(p_s_pole_lat) / 1000 )
  hdr_itm % realc (rc_polelong) =                                              &
    (first_multi % block_2(p_s_pole_lon) / 1000 ) + 180.00
  IF ( hdr_itm % realc (rc_polelong)>= 360.00 ) THEN
    hdr_itm % realc (rc_polelong)= hdr_itm % realc (rc_polelong)- 360.00
  END IF

END IF

!=======================================================================
!  Initialising Level Dependant Constants
!=======================================================================

hdr_itm % levdepc (:,:) = rmdi
current => lists(first_multi_list) % begin
! Set I = 1 so it can be used in the While loop in Case statment
! to index the LevDepC array.
i = 1
SELECT CASE (first_multi_type)

CASE (tb3_pressure)
  ! *100 to go from HPa to Pa
  !cdir novector
  DO WHILE (ASSOCIATED(current))
    hdr_itm % levdepc (i,ldc_pressure) =                                       &
      current % block_1 ( p_lvl_desc_1 ) * 100
    current => current % next
    i = i + 1
  END DO

CASE (tb3_hybrid  )
  ! Model levels keep model level index from ECMWF GRIB file
  !cdir novector
  DO WHILE (ASSOCIATED(current))
    hdr_itm % levdepc (i,ldc_mlindex) =                                        &
      current % block_1 ( p_lvl_desc_1 )
    current => current % next
    i = i + 1
  END DO

CASE default
  cmessage(1)    =                                                             &
    'Multi-level field type not set to correct value'
  errorstatus = 56
  CALL Ereport( routinename, errorstatus, cmessage(1) )

END SELECT

! Set soil level dependent constants
! Note: rcf_grib_check sets soil index in p_Lvl_Desc_1 and
!       soil depths in p_Lvl_Desc_2
IF (ASSOCIATED(lists(grib_soil_temp_field) % begin)) THEN
  current => lists(grib_soil_temp_field) % begin
  !cdir novector
  DO i = 1 , hdr_itm % intc (ic_soiltlevels)
    ! Level depth stored in cm, convert to m
    hdr_itm % levdepc (i,soildepths) =                                         &
      current % block_1 ( p_lvl_desc_2 ) / 100.
    current => current % next
  END DO
ELSE IF (ASSOCIATED(lists(grib_soil_moist_field) % begin)) THEN
  current => lists(grib_soil_moist_field) % begin
  !cdir novector
  DO i = 1 , hdr_itm % intc (ic_soilmoistlevs)
    ! Level depth stored in cm, convert to m
    hdr_itm % levdepc (i,soildepths) =                                         &
      current % block_1 ( p_lvl_desc_2 ) / 100.
    current => current % next
  END DO
END IF

!=======================================================================
!  Initialising Row Dependant Constants
!=======================================================================

hdr_itm % rowdepc( :,: ) = rmdi

!=======================================================================
!  Initialising Col Dependant Constants
!=======================================================================

hdr_itm % coldepc( :,: ) = rmdi

!=======================================================================
!  Initialising Addresses and Lengths in Lookup
!=======================================================================

!Clear the decks
hdr_itm % lookup(  1 : 45, : ) = 0
hdr_itm % lookup( 46 : 64, : ) = TRANSFER( 0.0, int_val)

!Set the data time - twice
hdr_itm % lookup( lbyr   , : ) = hdr_dmy % fixhd (fh_dtYear)
hdr_itm % lookup( lbmon  , : ) = hdr_dmy % fixhd (fh_dtMonth)
hdr_itm % lookup( lbdat  , : ) = hdr_dmy % fixhd (fh_dtDay)
hdr_itm % lookup( lbhr   , : ) = hdr_dmy % fixhd (fh_dtHour)
hdr_itm % lookup( lbmin  , : ) = hdr_dmy % fixhd (fh_dtMinute)
hdr_itm % lookup( lbyrd  , : ) = hdr_dmy % fixhd (fh_dtYear)
hdr_itm % lookup( lbmond , : ) = hdr_dmy % fixhd (fh_dtMonth)
hdr_itm % lookup( lbdatd , : ) = hdr_dmy % fixhd (fh_dtDay)
hdr_itm % lookup( lbhrd  , : ) = hdr_dmy % fixhd (fh_dtHour)
hdr_itm % lookup( lbmind , : ) = hdr_dmy % fixhd (fh_dtMinute)

!Set Field dimensions
hdr_itm % lookup( lbnpt  , : ) = output_grid % glob_p_row_length
hdr_itm % lookup( lbrow  , : ) = output_grid % glob_p_rows
hdr_itm % lookup( lblrec , : ) = output_grid % glob_p_rows *                   &
  output_grid % glob_p_row_length

!Set grid type
hdr_itm % lookup( lbcode , : )    = 1

!Set Data type to 'real' accross the board.
hdr_itm % lookup( data_type , : ) = 1

!Set Packing Type to 'no packing' accross the board.
hdr_itm % lookup( lbpack    , : ) = 0

!Set Header Release Number
hdr_itm % lookup( lbrel     , : ) = 2

!Set Internal Model No.
hdr_itm % lookup( model_code , : ) = 1

!Words 46 to 64 are reals stored as integers

! Real Lat and Long of Pseudo North Pole
hdr_itm % lookup( bplat , : ) = TRANSFER (                                     &
  hdr_itm % realc (rc_polelat) , int_val )

hdr_itm % lookup( bplon , : ) = TRANSFER (                                     &
  hdr_itm % realc (rc_polelong), int_val )

! Lat and Long spacing
! calculate sign of latitude spacing (from 'last - first' grid point)
sign = ( first_multi % block_2(p_latextrmpt) -                                 &
  first_multi % block_2(p_latgrdpnt1) )
sign = sign / ABS(sign)

hdr_itm % lookup( bdy   , : ) = TRANSFER (                                     &
  sign * hdr_itm % realc (rc_latspacing), int_val )

hdr_itm % lookup( bdx   , : ) = TRANSFER (                                     &
  hdr_itm % realc (rc_longspacing), int_val )

! Zeroth Lat and Long - (Zeroth not First, Hence subtraction)
hdr_itm % lookup( bzy   , : ) = TRANSFER(                                      &
  hdr_itm % realc (rc_firstlat) -                                              &
  ( sign * hdr_itm % realc (rc_latspacing) )                                   &
  , int_val )

hdr_itm % lookup( bzx   , : ) = TRANSFER (                                     &
  hdr_itm % realc (rc_firstlong) - hdr_itm % realc (rc_longspacing)            &
  , int_val )

! Real missing data indicator
hdr_itm % lookup( bmdi  , : ) = TRANSFER ( rmdi, int_val )

! Set the Lookup which holds the stash code

count = 0

! Loop across all lists
DO i = 1, grib_max_fields

  IF (ASSOCIATED(lists(i) % begin) ) THEN

    current => lists(i) % begin
    DO WHILE (ASSOCIATED(current))

      count = count + 1

      hdr_itm % lookup( item_code , count ) = current % stashcode
      hdr_itm % lookup( blev      , count ) =                                  &
        TRANSFER (REAL(current % block_1 (p_lvl_desc_1)),int_val)

      current => current % next

    END DO  ! members of a list

  END IF ! Associated(Lists(I) % Begin)

END DO ! I over all lists


RETURN

END SUBROUTINE Rcf_Grib_Sethdr
END MODULE Rcf_Grib_Sethdr_Mod
