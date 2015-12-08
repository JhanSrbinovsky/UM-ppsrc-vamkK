! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!   Subroutine REPLANCA ---------------------------------------------
!
!   Purpose:  Updates ancillary fields as requested in FIELDCODE array.
!     Tests whether update is required for each field, allowing for
!     dependencies between fields. Uses LOOKUP array to find data for
!     appropriate time, reads a record and checks for current data
!     type. Reads second record if time interpolation required. Updates
!     the field. Under DEF RECON, the interface to the routine is
!     modified for use in the reconfiguration rather than the model.
!     Under DEF CAL360 the 360 day rather than the Gregorian calender
!     is used.
!     -------------------------------------------------------------
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: Ancillaries

 SUBROUTINE replanca(i_year,i_month,i_day,i_hour,                              &
                     i_minute,i_second,i_day_number,                           &
                     ancil_reftime,offset_steps,                               &
                     p_field,p_rows,u_field,v_field,d1,land,                   &
                     a_step,land_field,steps_per_hr,                           &
                     fland_ctile,                                              &
                     tstar_land_ctile,tstar_sea_ctile,                         &
                     tstar_sice_ctile,                                         &
                     ice_fraction,tstar,tstar_anom,                            &
                     sm_levels,dz_soil,smc_updated,                            &
                     ns_space,first_lat,                                       &
                     len1_lookup,len_fixhd,len_inthd,                          &
                     len_realhd,len_d1,fixhd,inthd,realhd,                     &
                     lookup,rlookup,ftnancil,lookup_start,                     &
                     ndatasets,nlookups,                                       &
                     icode,cmessage,lcal360)        ! Intent Out

USE water_constants_mod, ONLY: tfs, tm
USE mask_compression, ONLY : compress_to_mask

USE ancilcta_namelist_mod, ONLY:                                               &
  l_sstanom, lamipii

USE switches, ONLY: l_ctile
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE domain_params
USE lookup_addresses
USE t_int_mod,   ONLY: t_int
USE t_int_c_mod, ONLY: t_int_c
USE um_input_control_mod,  ONLY: ltleads
USE ukca_option_mod, ONLY: l_ukca
USE science_fixes_mod, ONLY: l_error_ancil_struct
USE Submodel_Mod

IMPLICIT NONE

LOGICAL lcal360

INTEGER                                                                        &
       i_year,                                                                 &
                          ! Current Model Time
       i_month,                                                                &
                          !   "      "     "
       i_day,                                                                  &
                          !   "      "     "
       i_hour,                                                                 &
                          !   "      "     "
       i_minute,                                                               &
                          !   "      "     "
       i_second,                                                               &
                          !   "      "     "
       i_day_number,                                                           &
       ancil_reftime(6),                                                       &
                          ! Reference time for ancillary updating
       offset_steps,                                                           &
                          ! Offset in timesteps of ref. from basis


       a_step,land_field,steps_per_hr,                                         &


       p_field,                                                                &
                          ! Size of horizontal fields
       p_rows,                                                                 &
                          !
       u_field,                                                                &
                          !   "  "      "         "
       v_field,                                                                &
                          !   "  "      "         "
       ndatasets,                                                              &
                          ! Number of ancillary datasets
       nlookups,                                                               &
                          ! Number of lookup tables
       len_d1             ! Size of primary data array

INTEGER                                                                        &
       len1_lookup,                                                            &
                          ! First dimension of lookup table
       len_fixhd,                                                              &
                          ! Length of headers in data sets
       len_inthd,                                                              &
                          !
       len_realhd,                                                             &
                          !
       fixhd(len_fixhd,ndatasets),                                             &
                          ! Data set headers
       inthd(len_inthd,ndatasets),                                             &
                          !
       lookup(len1_lookup,nlookups),                                           &
                          ! Data set lookup tables
       ftnancil(ndatasets),                                                    &
                          ! FTN numbers of data sets
       lookup_start(ndatasets),                                                &
                          ! Start of lookup tables
                          ! referring to each data set.
       sm_levels          ! number of soil levels

REAL                                                                           &
       d1(len_d1),                                                             &
                          ! INOUT  Primary data array used to hold
                          !        all fields except TSTAR and
                          !        ICE_FRACTION
       ice_fraction(p_field),                                                  &
                          ! INOUT  Ice frac of sea part of grid
                          !        box, updated if requested
       fland_ctile(land_field),                                                &
                          ! IN     Fractional land on land pts.
       fland_g(p_field),                                                       &
                          ! WORK   Frac land over all points.
       tstar(p_field),                                                         &
                          ! INOUT  TSTAR:updated if requested
       tstar_land_ctile(p_field),                                              &
                          ! INOUT  as above, but for land.
       tstar_sea_ctile(p_field),                                               &
                          ! INOUT  as above, but for open sea.
       tstar_sice_ctile(p_field),                                              &
                          ! INOUT  as above, but for sea-ice.
       tstar_anom(p_field),                                                    &
                          ! INOUT  SST anomaly,formed in recon;
                          !        added if requested in model run
       realhd(len_realhd,ndatasets),                                           &
       rlookup(len1_lookup,nlookups),                                          &
       ns_space,                                                               &
                          ! NS latitude spacing
       first_lat,                                                              &
                          ! Latitude of first gridpoint
       dz_soil(sm_levels) ! OUT soil thicknesses

LOGICAL                                                                        &
       land(p_field),                                                          &
                          ! WORK LAND mask
       sea(p_field),                                                           &
                          ! WORK SEA mask
       ltstar_sice,                                                            &
                          ! IN TRUE if TSTAR_SICE has been read in
                          !         from input dump.
                          !         If FALSE set to TSTAR_SEA.
       smc_updated        ! OUT T if smc updated

INTEGER                                                                        &
       icode,                                                                  &
                          ! Return code
       iounit             !OUT I/O unit passed out in RECON mode

CHARACTER(LEN=80)                                                              &
       cmessage           ! Error message
! ----------------------------------------------------------------
!*L--------------------COMDECK  CANCILA ---------------------------
!
! Purpose : Contains index blocks for control of update of
!           ancillary fields.
!
! -------------------------------------------------------------------
!
! CANCMAXA Store maximum total no of atmosphere/slab ancillary fields.
! -------------------------------------------------------------------
! Type Declarations

      INTEGER NANCIL_FIELDS  ! No of Atmosphere & Slab Ancillary fields
PARAMETER (NANCIL_FIELDS = 196)
! CANCMAXA end
! Type Declarations

      INTEGER                                                           &
     &  FILEANCIL,                                                      &
                         ! File number associated with ancillary fields
     &  NLOOKUP,                                                        &
                         ! Position of ancillary field in lookup tables.
     &  LOOKUP_STEP,                                                    &
                         ! Interval between PP Headers refering to
!                        ! to the same ancillary fields at diferent time
     &  LEVELS,                                                         &
                         ! Number of levels of data in each ancillary
!                        ! field.
     &  STASHANCIL,                                                     &
                         ! Stash codes for ancillary files
     &  D1_ANCILADD      ! Address of ancillary field in main data block


      COMMON/IXANCILA/ FILEANCIL(NANCIL_FIELDS),                        &
     &           NLOOKUP(NANCIL_FIELDS),                                &
     &           LOOKUP_STEP(NANCIL_FIELDS),                            &
     &           LEVELS(NANCIL_FIELDS),                                 &
     &           STASHANCIL(NANCIL_FIELDS),                             &
     &           D1_ANCILADD(NANCIL_FIELDS)

!*L---------- Control data calculated from NAMELIST-------------------
      LOGICAL                                                           &
     &         UPDATE

      INTEGER  FIELDCODE,                                               &
     &         STEPS
!*----------------------------------------------------------------------
      COMMON/CTANCILA/                                                  &
     &         FIELDCODE(2,NANCIL_FIELDS),                              &
     &         STEPS(NANCIL_FIELDS),UPDATE(NANCIL_FIELDS)
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

!     Local real arrays

REAL                                                                           &
       ancil1(p_field),                                                        &
                          ! Buffers to hold values of ancillary
                          ! data for time interpolation.
       ancil2(p_field),                                                        &
                          !
       ancil_data(p_field),                                                    &
                          ! Field of ancillary data held prior
                          ! to selective updating.
       snow_change(p_field),                                                   &
                          ! Fractional time of change of
                          ! snow cover
       ice_extent(p_field,2),                                                  &
                          ! Fractional time of change
                          ! of ice cover
       pres_value(p_field),                                                    &
                          ! Prescribed value of data when
                          ! controlling field is zero.
       no_ice_extent(p_field),                                                 &
                          ! Indicator for no sea ice
                          ! =0 if ice cover
       tstar_land(p_field),                                                    &
                          !Temporary store for land surface temp.
       tstar_sea(p_field),                                                     &
                          !as above, but for open sea.
       tstar_sice(p_field),                                                    &
                          !as above, but for sea-ice.
       tstar_ssi(p_field)
                          !as above, but for sea mean.

!     Local variables

INTEGER                                                                        &
       i,                                                                      &
                          !
       i1,                                                                     &
                          !
       i2,                                                                     &
                          !
       i3,                                                                     &
       id,                                                                     &
                          !
       im,                                                                     &
                          !
       iy,                                                                     &
                          !
       k,                                                                      &
       l,                                                                      &
                          ! Land index
       field,                                                                  &
                          ! Current field number.
       file               !

INTEGER                                                                        &
       interval,                                                               &
                          ! Interval between data times
       step,                                                                   &
                          ! Number of data times skipped.
       months,                                                                 &
                          ! Used in calculation of position
                          ! of data required.
       hours,                                                                  &
                          !
       period,                                                                 &
                          ! Period of periodic data
       start_month,                                                            &
                          !
       level,                                                                  &
                          !
       nftin,                                                                  &
                          ! Current FTN number for ancillary field
       ancil_ref_days,                                                         &
                          ! Ancil.reference time in whole days
       ancil_ref_secs,                                                         &
                          ! Ancil.reference time in extra seconds
       day,sec,                                                                &
                          ! Times relative to reference time
       day1,sec1,                                                              &
                          ! Times relative to reference time
       incr_sec,                                                               &
                          ! Increment in sec
       LEN,                                                                    &
       iend,                                                                   &
       ii,row_length,j
INTEGER                                                                        &
       i_year1,                                                                &
                          ! Copy of Current Model Time year
       i_month1,                                                               &
                          !   "      "     "          month
       i_day1,                                                                 &
                          !   "      "     "          day
       i_hour1            !   "      "     "          hour

INTEGER                                                                        &
       update_months      ! update frequency (months) if Gregorian
LOGICAL                                                                        &
       lgreg_monthly      ! True for Gregorian monthly updating

! *IF -DEF,CAL360

INTEGER                                                                        &
       i_year_basis,                                                           &
                          ! Basis Model Time
       i_month_basis,                                                          &
                          !   "     "     "
       i_day_basis,                                                            &
                          !   "     "     "
       i_hour_basis,                                                           &
                          !   "     "     "
       i_minute_basis,                                                         &
                          !   "     "     "
       i_second_basis,                                                         &
                          !   "     "     "
       i_day_number_basis

! *ENDIF


INTEGER                                                                        &
       i_year_ref,                                                             &
                          ! Reference Time
       i_month_ref,                                                            &
                          !    "       "
       i_day_ref,                                                              &
                          !    "       "
       i_hour_ref,                                                             &
                          !    "       "
       i_minute_ref,                                                           &
                          !    "       "
       i_second_ref
                          !    "       "


LOGICAL                                                                        &
       linterpolate,                                                           &
                          ! Indicates whether time
                          ! interpolation needed.
       lt_int_c,                                                               &
                          ! Indicates use of controlled time
                          ! interpolation
       lmismatch,                                                              &
                          ! Used in header checks
       lice_fraction,                                                          &
                          !
       lsnow_depth,                                                            &
                          !
       single_time,                                                            &
                          ! Indicates that only one time is
                          ! available in data set
       periodic,                                                               &
                          ! Data set is periodic
       regular,                                                                &
                          ! Interval between data times in
                          ! dataset is regular in model timesteps.
       lice_depth

REAL                                                                           &
       zero,                                                                   &
                          !
       time1,                                                                  &
                          ! Times if data used in time interpolation
       time2,                                                                  &
                          !
       time,                                                                   &
                          !Target time for time interpolation
       lat_p
                          ! latitude of point

INTEGER :: lb_code
INTEGER :: bot_level
INTEGER :: exppxi
EXTERNAL exppxi

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!     Internal structure

!   List of Atmosphere & Slab Ancillary fields.
! CANCLSTA start
!
! Purpose : Cross-Reference List of Ancillary Fields
!           Atmosphere. Ancillary Reference numbers,
!           Stash Codes, Model Codes and Logical file Numbers
!
! -------------------------------------------------------------------
!   Column A  : Ancillary Reference Number
!          B  : Internal Model Number
!          C  : Stash Item code
!          D  : Logical file number
!          E  : Field Description
!
!   A  B    C   D  E
!   ----------------------------------------------------------------
!   1  1   30   9  Land Sea Mask
!   2  1   33  10  Orography
!   3  1   34  10  Orographic Variance
!   4  1   35  10  Orographic gradient XX
!   5  1   36  10  Orographic gradient XY
!   6  1   37  10  Orographic gradient YY
!   7  1   60   1  Ozone
!   8              Not used
!   9  1   23   2  Snow Depth
!  10  1   20   3  Deep Soil Temperature
!  11  1   40   4  Vol SMC at Wilting
!  12  1   41   4  Vol SMC at Critical Point
!  13  1  190  11  Oasis Iceberg Calving
!  14  1   43   4  Vol SMC at Saturation
!  15  1   44   4  Saturated Soil Conductivity
!  16              Not used
!  17  1   46   4  Thermal Capacity
!  18  1   47   4  Thermal Conductivity
!  19              Not used
!  20              Not used
!  21              Not used
!  22              Not used
!  23              Not used
!  24              Not used
!  25              Not used
!  26  1   26   5  Surface Roughness (used for sea points on first timestep)
!  27  1   31   7  Sea-ice Fraction
!  28  1   24   6  Sea-surface Temperature
!  29  1   32   7  Sea-ice Thickness
!  30  1   28   8  Surface Currents : u-component
!  31  1   29   8  Surface Currents : v-component
!  32  1   93   9  Runoff coastal outflow point
!  33  1  274  98  Mean topographic index
!  34  1  275  99  Standard deviation of topographic index
!  35  1   48   4  Saturated soil water suction
!  36  1    9   2  Soil moisture in layers
!  37              Not used (SLAB) 
!  38              Not used (SLAB) 
!  39  1   58  12  Sulphur dioxide emission
!  40  1   59  12  Dimethyl sulphur emission
!  41  1   88  13  Sulphate aerosol mass mixing ratio
!  42  1   87  13  Sulphuric acid aerosol mass mixing ratio
!  43  1   85  13  Soot aerosol mass mixing ratio
!  44  1   57  14  Multi-level murk source term emission
!  45  1   90  14  Multi-level murk concentration
!  46  1   17  10  Silhouette area of orography (orog. roughness)
!  47  1   18  10  Peak to trough height (for orog. roughness scheme)
!  48  1  301  15  User ancillary field 1
!  49  1  302  15  User ancillary field 2
!  50  1  303  15  User ancillary field 3
!  51  1  304  15  User ancillary field 4
!  52  1  305  15  User ancillary field 5
!  53  1  306  15  User ancillary field 6
!  54  1  307  15  User ancillary field 7
!  55  1  308  15  User ancillary field 8
!  56  1  309  15  User ancillary field 9
!  57  1  310  15  User ancillary field 10
!  58  1  311  15  User ancillary field 11
!  59  1  312  15  User ancillary field 12
!  60  1  313  15  User ancillary field 13
!  61  1  314  15  User ancillary field 14
!  62  1  315  15  User ancillary field 15
!  63  1  316  15  User ancillary field 16
!  64  1  317  15  User ancillary field 17
!  65  1  318  15  User ancillary field 18
!  66  1  319  15  User ancillary field 19
!  67  1  320  15  User ancillary field 20
!  68  1  127  12  NH3 (ammonia) aerosol emission
!  69  1  128  23  Surface fresh soot aerosol emission
!  70  1  129  23  Elevated fresh soot aerosol emission
!  71              Not used
!  72  1  121  17  Natural Sulphur dioxide emissions
!  73  1  122  18  OH concentrations
!  74  1  123  18  HO2 concentrations
!  75  1  124  18  H2O2 concentrations
!  76  1  125  18  Ozone (CHEM) concentrations
!  77  1  126  12  Sulphur dioxide high level emission
!  78  1  251  24  Surface CO2 emissions
!  79  1  207   4  Clapp-Hornberger parameter
!  80              Not used
!  81              Not used
!  82  1  160  19  Aerosol data for radiative forcing of climate change
!  83  1  216  20  Initial fractions of surface types
!  84  1  217  21  Initial leaf area index of plant functional types
!  85  1  218  21  Initial canopy height of plant functional types
!  86  1  213  21  Initial gridbox mean canopy conductance
!  87  1  219  22  Fraction of vegetation subject to disturbance
!  88  1  220   4  Snow free albedo of bare soil
!  89  1  223   4  Soil carbon content
!  90  1  321  16  User ancillary multi 1
!  91  1  322  16  User ancillary multi 2
!  92  1  323  16  User ancillary multi 3
!  93  1  324  16  User ancillary multi 4
!  94  1  325  16  User ancillary multi 5
!  95  1  326  16  User ancillary multi 6
!  96  1  327  16  User ancillary multi 7
!  97  1  328  16  User ancillary multi 8
!  98  1  329  16  User ancillary multi 9
!  99  1  330  16  User ancillary multi 10
! 100  1  331  16  User ancillary multi 11
! 101  1  332  16  User ancillary multi 12
! 102  1  333  16  User ancillary multi 13
! 103  1  334  16  User ancillary multi 14
! 104  1  335  16  User ancillary multi 15
! 105  1  336  16  User ancillary multi 16
! 106  1  337  16  User ancillary multi 17
! 107  1  338  16  User ancillary multi 18
! 108  1  339  16  User ancillary multi 19
! 109  1  340  16  User ancillary multi 20
! 110  1  341  25  Tropopause-based Ozone
! 111  1  505  26  Land Fraction
! 112  1  418  27  Dust parent soil clay fraction
! 113  1  419  27  Dust parent soil silt fraction
! 114  1  420  27  Dust parent soil sand fraction
! 115  1  421  27  Dust soil mass fraction div 1
! 116  1  422  27  Dust soil mass fraction div 2
! 117  1  423  27  Dust soil mass fraction div 3
! 118  1  424  27  Dust soil mass fraction div 4
! 119  1  425  27  Dust soil mass fraction div 5
! 120  1  426  27  Dust soil mass fraction div 6
! 121  1  130  28  Biomass surface emissions
! 122  1  131  28  Biomass elevated emissions
! 123  1  132  29  DMS concentration in seawater
! 124  1  153  30  River Water Storage
! 125  1  151  31  Riv Channel Sequence
! 126  1  152  31  Riv Channel Direction
!
! 127-154 not used in UM
!
! 155  1    5  10  Orographic X gradient
! 156  1    6  10  Orographic Y gradient 
! 157  1  351  38  Climatalogical biogenic aerosol mmr
! 158  1  352  39  Clim Biomass-burning (fresh) mmr
! 159  1  353  39  Clim Biomass-burning (aged) mmr
! 160  1  354  39  Clim Biomass-burning (in-cloud) mmr
! 161  1  355  40  Clim Black Carbon (fresh) mmr
! 162  1  356  40  Clim Black Carbon (aged) mmr
! 163  1  357  41  Clim Sea-salt (film mode) npm3
! 164  1  358  41  Clim Sea-salt (jet mode) npm3
! 165  1  359  42  Clim Sulphate (accumulation mode) mmr
! 166  1  360  42  Clim Sulphate (Aitken mode) mmr
! 167  1  361  42  Clim Sulphate (dissolved) mmr
! 168  1  362  43  Clim Dust size division 1 mmr
! 169  1  363  43  Clim Dust size division 2 mmr
! 170  1  364  43  Clim Dust size division 3 mmr
! 171  1  365  43  Clim Dust size division 4 mmr
! 172  1  366  43  Clim Dust size division 5 mmr
! 173  1  367  43  Clim Dust size division 6 mmr
! 174  1  368  44  Clim Organic Carbon from Fossil Fuels (fresh) mmr
! 175  1  369  44  Clim Organic Carbon from Fossil Fuels (aged) mmr
! 176  1  370  44  Clim Organic Carbon from Fossil Fuels (in-cloud) mmr
! 177  1  371  45  Clim Delta Aerosol mmr
! 178  1  480  46  Prognostic Ozone Tracer Cariolle Scheme
! 179  1  481  46  Cariolle Ozone Production - Loss (P-L)
! 180  1  482  46  Cariolle Ozone P-L wrt Ozone Mix Ratio
! 181  1  483  46  Cariolle Ozone Volume Mixing Ratio
! 182  1  484  46  Cariolle Ozone P-L wrt Temperature
! 183  1  485  46  Cariolle Ozone Clim Temp
! 184  1  486  46  Cariolle Ozone P-L wrt Ozone above a point
! 185  1  487  46  Cariolle Ozone Column above a point
! 186  1  134  47  Surface fresh fossil-fuel OC aerosol emissions
! 187  1  135  47  Elevated fresh fossil-fuel OC aerosol emissions
! 188  1    7  10  Unfiltered orography
! 194  1  243   5  Obs/clim snow-free surface sw albedo  
! 195  1  244   5  Obs/clim snow-free surface vis albedo  
! 196  1  245   5  Obs/clim snow-free surface nir albedo
!  ------------------------------------------------------------------
! CANCLSTA end


!   1.  Initialisation for atmosphere
IF (lhook) CALL dr_hook('REPLANCA',zhook_in,zhook_handle)
icode=0
iounit=0
smc_updated=.FALSE.
update_months=0
incr_sec = 0
! These are printed out and possibly not always required.
hours    = imdi
interval = imdi
period   = imdi

!     Set up surface temperatures:

IF(l_ctile)THEN
  DO i=1,p_field
    tstar_land(i)=tstar_land_ctile(i)
    tstar_sea(i)=tstar_sea_ctile(i)
    tstar_sice(i)=tstar_sice_ctile(i)
    IF(ice_fraction(i) <= 0.0)THEN
      tstar_ssi(i)=tstar_sea(i)
    ELSE
      tstar_ssi(i)=ice_fraction(i)*tstar_sice(i)                               &
        +(1.0-ice_fraction(i))*tstar_sea(i)
    END IF
  END DO
ELSE
  DO i=1,p_field
    tstar_land(i)=tstar(i)
    tstar_ssi(i)=tstar(i)
  END DO
END IF


!     Initialise ANCIL1/2. Includes Halos for MPP runs.
DO i=1,p_field
  ancil1(i)=0.0
  ancil2(i)=0.0
END DO
!   1.1 Set logical UPDATE for each ancillary field independently

DO field=1,nancil_fields


  update(field)=.FALSE.
  IF(steps(field) /= 0) THEN
!         UPDATE(FIELD)=MOD(A_STEP,STEPS(FIELD)) == 0
  update(field)=(MOD(a_step+offset_steps,steps(field)) == 0                    &
                 .OR.a_step == 0)                                              &
                  .AND.fieldcode(1,field) >  0                                 &
                   .AND.d1_anciladd(field) >  1
  END IF

!   1.05 Copy ancillary updating reference time to local variables
  i_year_ref   = ancil_reftime(1)
  i_month_ref  = ancil_reftime(2)
  i_day_ref    = ancil_reftime(3)
  i_hour_ref   = ancil_reftime(4)
  i_minute_ref = ancil_reftime(5)
  i_second_ref = ancil_reftime(6)
!        and convert to reference days & secs
! DEPENDS ON: time2sec
  CALL time2sec(i_year_ref,i_month_ref,i_day_ref,                              &
                i_hour_ref,i_minute_ref,i_second_ref,                          &
                0,0,ancil_ref_days,ancil_ref_secs,lcal360)

  IF (.NOT. lcal360) THEN

!   1.11 Set logical UPDATE for Gregorian calender updates at monthly
!        or yearly intervals. NB STEPS value set to 1 day in INANCILA
    IF(fieldcode(1,field) == 1.OR.fieldcode(1,field) == 2) THEN
      months=i_month+i_year*12-(i_month_ref+i_year_ref*12)
      update_months= fieldcode(2,field)*                                       &
     ((3-fieldcode(1,field))/2 *12+ 1-(3-fieldcode(1,field))/2)
     update(field)=MOD(months,update_months) == 0.AND.i_day == 1
    END IF
  END IF !  (.NOT.LCAL360)

END DO

!  1.2 Allow for dependencies between fields
! Sea surface temperature must be updated when sea ice is updated

update(28)=update(27).OR.update(28)

! Both surface current components must be updated together

update(30)=update(30).OR.update(31)
update(31)=update(30)

!  Select method of time interpolation for SST. The interpolation
!  allows for sea ice if ice data is available at the same times
!  as the temperature data. Otherwise linear interpolation is used.

lt_int_c=.TRUE.

IF(update(28)) THEN
IF(fixhd(10,fileancil(27)) == 0) lt_int_c=.FALSE.
  IF(lt_int_c) THEN
  DO i=21,41
    IF(fixhd(i,fileancil(27)) /= fixhd(i,                                      &
      fileancil(28))) THEN
      lt_int_c=.FALSE.
      WRITE(6,'(A,A,A)')' WARNING:controlled time interpolation for SST', &
          ' not available: Mismatch in SST and SEA-ICE ancillary data',&
          ' times in FIXED HEADER'
      WRITE(6,*)' position=',i,' SEA-ICE=',fixhd(i,fileancil(27))
      WRITE(6,*)' position=',i,' SST=',fixhd(i,fileancil(28))
    END IF
  END DO
  END IF
END IF


! Read in fractional land field
! Set up global fractional land field
IF(l_ctile)THEN
  l=0
  DO i=1,p_field
    fland_g(i)=0.0
    IF(land(i))THEN
      l=l+1
      fland_g(i)=fland_ctile(l)
   END IF
  END DO
ELSE
  DO i=1,p_field
    IF(land(i))THEN
      fland_g(i)=1.0
    ELSE
      fland_g(i)=0.0
    END IF
  END DO
END IF


DO i=1,p_field
  sea(i)=.FALSE.
  IF(fland_g(i) <  1.0)sea(i)=.TRUE.
END DO

!  Loop over ancillary fields(atmosphere)

DO field=1,nancil_fields

  lice_depth=field == 29  ! required for LAMIPII

IF (update(field)) THEN  ! (1st level IF)
  FILE=fileancil(field)
  nftin=ftnancil(FILE)

  IF(lice_depth.AND.lamipii) THEN

! Uses ice fraction set earlier in field loop.
! WARNING this will fail if the order of ancillary fields is ever
! changed so that ice-depth preceeds ice fraction
! Note : For complete sea ice cover
!        Arctic ice depth    = 2m
!        Antarctic ice depth = 1m
! For ice concentrations less than 1. ice depth is 1 or 2 times conc.
! This results in similar values to those from runs using ancillary
! files containing ice depths set to 1 or 2m.

    row_length=p_field/p_rows
    DO i=1,p_rows
! work out latitude in radians
      lat_p=first_lat+ns_space*(i+datastart(2)-offy-1)
      DO j=1,row_length
        ii=j+(i-1)*row_length
        ancil_data(ii)=0.0
        IF (ice_fraction(ii) >  0.0) THEN
          IF (lat_p >  0.0) THEN   ! Arctic ice depth
            ancil_data(ii)=2.*ice_fraction(ii)
          ELSE                     ! Antarctic ice depth
            ancil_data(ii)=1.*ice_fraction(ii)
          END IF
        END IF
      END DO
    END DO
!      Sea ice thickness
!        Update over all sea points (all sea ice points are the only
!        ones strictly required, but this cannot be determined easily)

    DO i=1,p_field
      IF(sea(i)) THEN
        d1(d1_anciladd(field)+i-1)=ancil_data(i)
      END IF
    END DO
  ELSE
!     Update required for field

  WRITE(6,*)'REPLANCA: UPDATE REQUIRED FOR FIELD',field

    IF ( fixhd(10,FILE)  <   0 .OR. fixhd(10,FILE)  >   2 ) THEN
      icode = 700 + field
      cmessage = 'REPLANCA: error in fixed header(10) of ancillary file'
      GOTO 9999
    END IF

!     Check whether more than one data time available in data set

  single_time=fixhd(10,FILE) == 0

!     Set default values for time interpolation

  linterpolate=.TRUE.
  IF(single_time) THEN
    linterpolate=.FALSE.
  END IF

  IF (field >  9 .AND. field <  19) THEN
    linterpolate=.FALSE.
  END IF

!  2.1 Find position of input record

!     Default settings of search parameters if only one time present

  IF(single_time) THEN
    step=0
  ELSE


    lgreg_monthly=.FALSE.

IF (.NOT. lcal360) THEN
    IF(fieldcode(1,field) == 1.OR.fieldcode(1,field) == 2) THEN
      lgreg_monthly=.TRUE.
      update_months= fieldcode(2,field)*                                       &
      ((3-fieldcode(1,field))/2 *12+ 1-(3-fieldcode(1,field))/2)
    END IF
END IF



    periodic=fixhd(10,FILE) == 2
    regular=.TRUE.


IF (.NOT. lcal360) THEN
    regular=fixhd(35,FILE) == 0.AND.fixhd(36,FILE) == 0
! i.e. data at intervals of days/hours & non-periodic
    IF(periodic) regular=regular.AND.fixhd(37,FILE) == 0
! i.e. data at intervals of hours & periodic
END IF


!         Error checking on time information.

    IF ( fixhd(35,FILE)  <   0 .OR.                                            &
         fixhd(36,FILE)  <   0 .OR. fixhd(36,FILE) >  12 .OR.                  &
         regular .AND. ( fixhd(37,FILE) <  0 .OR. fixhd(37,FILE) > 31          &
        .OR. fixhd(38,FILE)  <   0 .OR. fixhd(38,FILE) >  24 ) ) THEN
!           FIXHD(39-40) are not used by REPLANCA.
!           FIXHD(35-37) have already been used if not CAL360.
      icode = 700 + field
      WRITE(cmessage,*) 'REPLANCA: Error in validity time interval',           &
                        ' given in ancil file (FIXHD(35-38))'
      GOTO 9999
    END IF

    IF ( fixhd(21,FILE)  <   0 .AND. .NOT. periodic                            &
  .OR. .NOT. ( regular .AND. periodic ) .AND.                                  &
!    !  If it is REGULAR & PERIODIC more detailed check is applied below
     ( fixhd(22,FILE)  <   0 .OR. fixhd(22,FILE)  >   12 .OR.                  &
       fixhd(23,FILE)  <   0 .OR. fixhd(23,FILE)  >   31 .OR.                  &
       fixhd(24,FILE)  <   0 .OR. fixhd(24,FILE)  >   24 .OR.                  &
       fixhd(25,FILE)  <   0 .OR. fixhd(25,FILE)  >   60 .OR.                  &
       fixhd(26,FILE)  <   0 .OR. fixhd(26,FILE)  >   60 ) ) THEN
      icode = 700 + field
      WRITE(cmessage,*) 'REPLANCA: Error in first validity time given',        &
                        ' in ancillary file (FIXHD(21-26))'
      GOTO 9999
    END IF

    period=-1
    IF(.NOT.periodic) THEN

!             If data taken from full time series of input data.

! DEPENDS ON: time2sec
      CALL time2sec(i_year,i_month,i_day,i_hour                                &
                    ,i_minute,i_second                                         &
                    ,ancil_ref_days,ancil_ref_secs,day,sec                     &
                    ,lcal360)


!  Adjust time to middle of updating interval

      IF(.NOT.lgreg_monthly) THEN
        sec=sec+steps(field)*1800/steps_per_hr

!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
        IF (a_step == 0) THEN
          day1 = day
          sec1 = sec
       incr_sec=-3600*MOD(offset_steps,steps(field))/steps_per_hr
! DEPENDS ON: time_df
          CALL time_df(day1,sec1,0,incr_sec,day,sec)
        END IF

      ELSE
        im=MOD(i_month+update_months-1,12) + 1
        iy=i_year+(i_month+update_months-1)/12
! DEPENDS ON: time2sec
        CALL time2sec(iy,im,i_day,i_hour                                       &
                    ,i_minute,i_second                                         &
                    ,ancil_ref_days,ancil_ref_secs,day1,sec1                   &
                    ,lcal360)
        IF (MOD(day+day1,2) == 0) THEN
          day=(day+day1)/2
          sec=(sec+sec1)/2
        ELSE
          day=(day+day1-1)/2
          sec=(sec+sec1+86400)/2
        END IF
!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
        IF (a_step == 0) THEN
          day1 = day
          sec1 = sec
       incr_sec=-3600*MOD(offset_steps,steps(field))/steps_per_hr
! DEPENDS ON: time_df
          CALL time_df(day1,sec1,0,incr_sec,day,sec)
        END IF
      END IF


      IF(regular) THEN
!  2.1.1  Standard cases:360 day calender;
!  2.1.1  or Gregorian calendar with
!         interval between data times in days or hours
!         updating interval may be regular in model timesteps,
!         or (LGREG_MONTHLY=T) irregular in model timesteps,

        hours=sec/3600+day*24
!  FInd time(in hours) of first ancillary data on file
! DEPENDS ON: time2sec
        CALL time2sec(fixhd(21,FILE),fixhd(22,FILE),                           &
                   fixhd(23,FILE),fixhd(24,FILE),                              &
                   fixhd(25,FILE),fixhd(26,FILE),                              &
                   ancil_ref_days,ancil_ref_secs,day,sec,                      &
                   lcal360)
        hours=hours-sec/3600-day*24

        IF(hours <  0) THEN
          icode=400+field
          cmessage='REPLANCA: Current time precedes start time of data'
          GOTO 9999
        END IF

!  FInd interval(in hours) between ancillary data on file
        interval=fixhd(35,FILE)*8640+fixhd(36,FILE)*720+                       &
                fixhd(37,FILE)*24+fixhd(38,FILE)

! Do not interpolate in time if data time exactly matches model time

        IF(MOD(hours,interval) == 0) THEN
          linterpolate=.FALSE.
        END IF

        step=hours/interval
        time=REAL(hours)
        time1=step*interval
        time2=(step+1)*interval

      ELSE

!  2.1.2 Gregorian calender;ancillary data interval is in months or
!        years,which is irregular in model timesteps.
!  original code is inaccurate for this section - corrected code under
!  LAMIPII makes use of dates in lookup headers
!  For a real calendar year the mid-point of each month is different
!  in terms of its hour and day. The old inaccurate method assumes
!  the hour and day are taken from the fixhd values. These are only
!  usually correct for the first month on the ancillary file.


!  Adjust YMD time to middle of updating interval

        i_year1=i_year
        i_month1=i_month
        i_day1=i_day
        i_hour1=i_hour
! DEPENDS ON: sec2time
        CALL sec2time(day,sec,ancil_ref_days,ancil_ref_secs,                   &
                     i_year,i_month,i_day,                                     &
                     i_hour,i_minute,i_second,i_day_number,                    &
                     lcal360)


!  FInd interval(in months) between ancillary data on file
        interval=fixhd(35,FILE)*12+fixhd(36,FILE)
        months=i_year*12+i_month
        start_month=fixhd(21,FILE)*12+fixhd(22,FILE)
        months=months-start_month
!  Check for time within month
     IF (lamipii) THEN   ! corrected code uses pp header
        step=months/interval
        i2=nlookup(field)+lookup_step(field)*step
        i1=i2+lookup_start(FILE)-1
! Check against day and hour of actual lookup header not first field
        IF((i_day*24+i_hour) <                                                 &
           (lookup(3,i1)*24+lookup(4,i1))) THEN
          months=months-1
        END IF
     ELSE              ! old less accurate code uses FIXHD
        IF((i_day*24+i_hour) <                                                 &
           (fixhd(23,FILE)*24+fixhd(24,FILE))) THEN
          months=months-1
        END IF
     END IF ! LAMIPII

        IF(months <  0) THEN
          icode=400+field
          cmessage='REPLANCA: Current time precedes start time of data'
          GOTO 9999
        END IF


!  Adjust YMD time back to start of updating interval

        i_year=i_year1
        i_month=i_month1
        i_day=i_day1
        i_hour=i_hour1



        step=months/interval

     IF (lamipii) THEN       ! corrected code
        time=REAL(sec)/3600+REAL(day*24)
! correct calculation of dates uses lookup table dates not fixhd date
        i2=nlookup(field)+lookup_step(field)*step
        i1=i2+lookup_start(FILE)-1
        i_year1=lookup(1,i1)
        i_month1=lookup(2,i1)
        i_day1=lookup(3,i1)
        i_hour1=lookup(4,i1)
! DEPENDS ON: time2sec
        CALL time2sec(i_year1,i_month1,i_day1,i_hour1,                         &
              fixhd(25,FILE),fixhd(26,FILE),                                   &
              ancil_ref_days,ancil_ref_secs,day,sec,                           &
              lcal360)
        time1=REAL(sec)/3600+REAL(day*24)
! I1+LOOKUP_STEP(irec) correct pointer to next field as multi-level fields &
! are possible.
        ! First check next lookup is within lookups
        IF (i1+lookup_step(field) < nlookups) THEN
          ! Second check next lookup is the same field.
          IF (lookup(item_code,i1) ==                                          &
              lookup(item_code,i1+lookup_step(field))) THEN
            i_year1=lookup(1,i1+lookup_step(field))
            i_month1=lookup(2,i1+lookup_step(field))
            i_day1=lookup(3,i1+lookup_step(field))
            i_hour1=lookup(4,i1+lookup_step(field))
          ELSE
            icode = 500
            WRITE(cmessage,'(A,I5)')                                           &
              'REPLANCA: error finding next lookup entry for field:', field
            GOTO 9999
          END IF
        ELSE
          icode = 500
          WRITE(cmessage,'(A,I5)')                                             &
            'REPLANCA: error due to lookup out of range for field:', field
          GOTO 9999
        END IF



! DEPENDS ON: time2sec
        CALL time2sec(i_year1,i_month1,i_day1,i_hour1,                         &
              fixhd(25,FILE),fixhd(26,FILE),                                   &
              ancil_ref_days,ancil_ref_secs,day,sec,                           &
              lcal360)
        time2=REAL(sec)/3600+REAL(day*24)

     ELSE   ! LAMIPII test - old inaccurate code using FIXHD
! NB INTERVAL may be > 1 month
        months=step*interval
! Calculate data times for time interpolation
        time=REAL(sec)/3600+REAL(day*24)
        im=MOD(fixhd(22,FILE)+months-1,12)+1
        iy=fixhd(21,FILE)+(months+fixhd(22,FILE)-1)/12
! DEPENDS ON: time2sec
        CALL time2sec(iy,im,fixhd(23,FILE),fixhd(24,FILE),                     &
              fixhd(25,FILE),fixhd(26,FILE),                                   &
              ancil_ref_days,ancil_ref_secs,day,sec,                           &
              lcal360)
        time1=REAL(sec)/3600+REAL(day*24)
        im=MOD(fixhd(22,FILE)+months+interval-1,12)+1
        iy=fixhd(21,FILE)+(months+interval+fixhd(22,FILE)-1)/12
! DEPENDS ON: time2sec
        CALL time2sec(iy,im,fixhd(23,FILE),fixhd(24,FILE),                     &
              fixhd(25,FILE),fixhd(26,FILE),                                   &
              ancil_ref_days,ancil_ref_secs,day,sec,                           &
              lcal360)
        time2=REAL(sec)/3600+REAL(day*24)
     END IF     ! end LAMIPII test

! Do not interpolate in time if data time exactly matches model time

        IF(time == time1) THEN
          linterpolate=.FALSE.
        END IF

      END IF ! End of REGULAR/not REGULAR

    ELSE  ! PERIODIC data

!  2.2   If data is taken from ancillary periodic data.

! DEPENDS ON: time2sec
      CALL time2sec(i_year,i_month,i_day,i_hour,                               &
                     i_minute,i_second,                                        &
                     ancil_ref_days,ancil_ref_secs,day,sec,                    &
                     lcal360)


!  Adjust time to middle of updating interval

      IF(.NOT.lgreg_monthly) THEN
        sec=sec+steps(field)*1800/steps_per_hr

!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
        IF (a_step == 0) THEN
          day1 = day
          sec1 = sec
       incr_sec=-3600*MOD(offset_steps,steps(field))/steps_per_hr
! DEPENDS ON: time_df
          CALL time_df(day1,sec1,0,incr_sec,day,sec)
        END IF

      ELSE
        im=MOD(i_month+update_months-1,12) + 1
        iy=i_year+(i_month+update_months-1)/12
! DEPENDS ON: time2sec
        CALL time2sec(iy,im,i_day,i_hour                                       &
                    ,i_minute,i_second                                         &
                    ,ancil_ref_days,ancil_ref_secs,day1,sec1                   &
                    ,lcal360)
        IF (MOD(day+day1,2) == 0) THEN
          day=(day+day1)/2
          sec=(sec+sec1)/2
        ELSE
          day=(day+day1-1)/2
          sec=(sec+sec1+86400)/2
        END IF
!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
        IF (a_step == 0) THEN
          day1 = day
          sec1 = sec
       incr_sec=-3600*MOD(offset_steps,steps(field))/steps_per_hr
! DEPENDS ON: time_df
          CALL time_df(day1,sec1,0,incr_sec,day,sec)
        END IF
      END IF


!  Adjust YMD time to middle of updating interval

      i_year1=i_year
      i_month1=i_month
      i_day1=i_day
      i_hour1=i_hour
! DEPENDS ON: sec2time
      CALL sec2time(day,sec,ancil_ref_days,ancil_ref_secs,                     &
                     i_year,i_month,i_day,                                     &
                     i_hour,i_minute,i_second,i_day_number,                    &
                     lcal360)



      IF (regular) THEN
!  2.2.1 Standard cases:1) 360 day calender, with allowed periods of
!        1 day, 1 month or 1 year;
!
!        2) Gregorian calender with update in hours,and period of
!        data 1 day.
!
!        For both updating interval and number of
!        data times to be skipped in data set calculated in hours.

        hours=sec/3600+day*24
        interval=fixhd(35,FILE)*8640+fixhd(36,FILE)*720+                       &
                fixhd(37,FILE)*24+fixhd(38,FILE)

        period=inthd(3,FILE)*interval

!    Do not allow non-standard periods
IF (lcal360) THEN
        IF(period /= 8640.AND.period /= 720.AND.period /= 24)THEN
          icode=600+field
          cmessage='REPLANCA: Non-standard period for periodic data'
          GOTO 9999
        END IF
ELSE
        IF(period /= 24)THEN
          icode=600+field
          cmessage='REPLANCA: Non-standard period for periodic data'
          GOTO 9999
        END IF
END IF
        IF(period == 24)THEN
! Ancillary data interval in hour(s), period is 1 day

          iy=i_year
          im=i_month
          id=i_day
          IF(i_hour <  fixhd(24,FILE)) hours=hours+24

        ELSE IF(period == 720)THEN
! Ancillary data interval in day(s) or hours , period is 1 month

          iy=i_year
          im=i_month
          id=fixhd(23,FILE)
          IF((i_day*24+i_hour) <                                               &
             (fixhd(23,FILE)*24+fixhd(24,FILE)))                               &
           hours=hours+720

        ELSE IF(period == 8640)THEN
! Ancillary data interval in month(s)or days or hours, period is 1 year

          iy=i_year
          im=fixhd(22,FILE)
          id=fixhd(23,FILE)
          IF((i_month*720+i_day*24+i_hour) <                                   &
          (fixhd(22,FILE)*720+fixhd(23,FILE)*24+fixhd(24,FILE)))               &
           hours=hours+8640

        END IF

! DEPENDS ON: time2sec
        CALL time2sec(iy,im,id,fixhd(24,FILE),                                 &
                     fixhd(25,FILE),fixhd(26,FILE),                            &
                     ancil_ref_days,ancil_ref_secs,day,sec,                    &
                     lcal360)
        hours=hours-sec/3600-day*24

! Do not interpolate in time if data time exactly matches model time

        IF(MOD(hours,interval) == 0) THEN
          linterpolate=.FALSE.
        END IF
        step=hours/interval
        time=REAL(hours)
        time1=step*interval
        time2=(step+1)*interval

      ELSE  ! non regular case

!  2.2.2 Gregorian calender,and data interval is in months,
!        period is 1 year
!        Updating interval and number of data times to be skipped
!        calculated in months.

        time=REAL(sec)/3600+REAL(day*24)
        interval=fixhd(36,FILE)+fixhd(35,FILE)*12
        period=inthd(3,FILE)*interval
        IF(period /= 12)THEN
          icode=600+field
          cmessage='REPLANCA: Non-standard period for periodic data'
          GOTO 9999
        END IF
!  Difference between date now (month) & first date ancil file (month)
        months=i_month-fixhd(22,FILE)


     IF (lamipii) THEN ! correct code to use lookup header dates
! Correctly use day and hour from lookup header not fixhd which
! contains values for first field on ancillary file only.
       step=months/interval
       i2=nlookup(field)+lookup_step(field)*step
       i1=i2+lookup_start(FILE)-1
!  Check for time within month - using ppheader information
    IF((i_day*24+i_hour) <  (lookup(3,i1)*24+lookup(4,i1))) THEN
          months=months-1
    END IF
       IF(months <  0) THEN
          months=months+12
       END IF
! recalculate STEP
        step=months/interval
! NB INTERVAL may be > 1 month
        months=step*interval
        iy=i_year
        im=MOD(fixhd(22,FILE)+months-1,12)+1
        IF(im >  i_month) iy=iy-1
        i2=nlookup(field)+lookup_step(field)*step
        i1=i2+lookup_start(FILE)-1
! DEPENDS ON: time2sec
        CALL time2sec(iy,im,lookup(3,i1),lookup(4,i1),                         &
              fixhd(25,FILE),fixhd(26,FILE),                                   &
              ancil_ref_days,ancil_ref_secs,day,sec,lcal360)
        time1=REAL(sec)/3600+REAL(day*24)
!  Calculate  TIME2 for second ancillary data time
!  set IY correctly for time interpolation calculations
        iy=i_year
        im=MOD(fixhd(22,FILE)+months+interval-1,12)+1
        IF(im <  i_month) iy=iy+1
        i1=(im-1)/interval
        i2=nlookup(field)+lookup_step(field)*i1
        i1=i2+lookup_start(FILE)-1
! DEPENDS ON: time2sec
        CALL time2sec(iy,im,lookup(3,i1),lookup(4,i1),                         &
              fixhd(25,FILE),fixhd(26,FILE),                                   &
              ancil_ref_days,ancil_ref_secs,day,sec,lcal360)
        time2=REAL(sec)/3600+REAL(day*24)

     ELSE   ! original code inaccurate use of FIXHD dates
!  Check for time within month
        IF((i_day*24+i_hour) <                                                 &
           (fixhd(23,FILE)*24+fixhd(24,FILE))) THEN
          months=months-1
        END IF
        IF(months <  0) THEN
          months=months+12
        END IF

        step=months/interval
! NB INTERVAL may be > 1 month
        months=step*interval
!  Calculate TIME1 for first ancillary data time
!  set IY correctly for time interpolation calculations
        iy=i_year
        im=MOD(fixhd(22,FILE)+months-1,12)+1
        IF(im >  i_month) iy=iy-1
! DEPENDS ON: time2sec
        CALL time2sec(iy,im,fixhd(23,FILE),fixhd(24,FILE),                     &
              fixhd(25,FILE),fixhd(26,FILE),                                   &
              ancil_ref_days,ancil_ref_secs,day,sec,                           &
              lcal360)
        time1=REAL(sec)/3600+REAL(day*24)
!  Calculate  TIME2 for second ancillary data time
!  set IY correctly for time interpolation calculations
        iy=i_year
        im=MOD(fixhd(22,FILE)+months+interval-1,12)+1
        IF(im <  i_month) iy=iy+1
! DEPENDS ON: time2sec
        CALL time2sec(iy,im,fixhd(23,FILE),fixhd(24,FILE),                     &
              fixhd(25,FILE),fixhd(26,FILE),                                   &
              ancil_ref_days,ancil_ref_secs,day,sec,                           &
              lcal360)
        time2=REAL(sec)/3600+REAL(day*24)
     END IF  ! end LAMIPII test

! Do not interpolate in time if data time exactly matches model time

        IF(time == time1) THEN
          linterpolate=.FALSE.
        END IF

      END IF  ! regular/non-regular


!  Adjust YMD time back to start of updating interval

      i_year=i_year1
      i_month=i_month1
      i_day=i_day1
      i_hour=i_hour1


    END IF  ! non-periodic/periodic

  IF (linterpolate) THEN
  WRITE(6,*)' REPLANCA - time interpolation for field ',field
  WRITE(6,*)' time,time1,time2 ',time,time1,time2
  WRITE(6,*)' hours,int,period ',hours,interval,period
  END IF

  END IF ! singletime/non-singletime

!  2.3   Check STASH Code

  i2=nlookup(field)+lookup_step(field)*step

  i1=lookup(item_code,i2+lookup_start(FILE)-1)

  lmismatch=.FALSE.
  WRITE(6,*)' Information used ',                                              &
      'in checking ancillary data set:',                                       &
      ' position of lookup table in dataset:',i2
  WRITE(6,*)' Position of first ',                                             &
      'lookup table referring to ',                                            &
      'data type ',nlookup(field)
  WRITE(6,*)' Interval between lookup ',                                       &
      'tables referring to data ',                                             &
      'type ', LOOKUP_STEP(FIELD),' Number of steps', STEP
  WRITE(6,*)' STASH code in dataset ',i1,                                      &
      '  STASH code requested ',stashancil(field)
  WRITE(6,*)'''start'' position of lookup tables for dataset ',                &
      'in overall lookup array ' ,lookup_start(FILE)

  IF(i1 /= stashancil(field)) THEN
  WRITE(6,*)i1,stashancil(field),field
    lmismatch=.TRUE.
  END IF

!  Error exit if checks fail

  IF(lmismatch) THEN
    icode=200+field
    cmessage='REPLANCA: PP HEADERS ON ANCILLARY FILE DO NOT MATCH'
    GOTO 9999
  END IF

  IF(linterpolate.AND..NOT.single_time) THEN
!  Check time interpolation factors
    IF(time <  time1.OR.time >  time2) THEN
     WRITE(6,*)' Information used in interpolation/replacement:'
     WRITE(6,*)' Time of first data=', time1
     WRITE(6,*)' Validity Time for update=', time
     WRITE(6,*)' Time of second data=', time2

     icode=500+field
     cmessage='REPLANCA: TIME INTERPOLATION ERROR'
     GOTO 9999
    END IF
  END IF

!  3   Loop over levels of ancillary data for field I
!  Reset pointer for dataset


!  Includes loop over X and Y components of surface currents

   lice_fraction=field == 27
   lsnow_depth=field == 9
   lice_depth=field == 29

! Find whether we are on theta level 0
   bot_level = -1
   lb_code = -1
   icode   = 0
! DEPENDS ON: exppxi
   lb_code = EXPPXI(1,stashancil(field)/1000,MOD(stashancil(field),1000),9, &
                    icode,CMESSAGE)
   
   IF (lb_code > 0) THEN
! DEPENDS ON: levcod
     CALL LEVCOD(lb_code,bot_level,icode,CMESSAGE)
   END IF
 
   IF (icode /= 0) THEN
     RETURN
   END IF

   IF (bot_level /= 0) THEN
     bot_level = 1
   END IF
   
  DO level=1,levels(field)

!  Do not go through loop for ice edge or snow edge

  IF((lice_fraction.OR.lsnow_depth).AND.level == 2) THEN
    CYCLE
  END IF

!  3.1 Read data for single level of ancillary field.

  IF(.NOT.lice_fraction) THEN
! AMIPII case ice depth field not read from ancillary file
   IF(.NOT.(lice_depth.AND.lamipii)) THEN
! DEPENDS ON: readflds
    CALL readflds(nftin,1,i2,lookup(1,lookup_start(FILE)),                     &
                  len1_lookup,ancil1,p_field,fixhd(1,FILE),                    &
                  icode,cmessage)

   END IF
    IF(icode /= 0)THEN
      icode=field+100
      iounit=nftin
      cmessage='REPLANCA :I/O ERROR '
      GOTO 9999
    END IF

  ELSE

!  If ice-fraction,read fractional time field as well
!        UNLESS IT IS A SINGLE TIME FIELD
!  If snow-depth,read fractional time field as well only if time
!  interpolation required.

IF(.NOT.single_time.AND..NOT.lamipii) THEN
   IF(lookup(item_code,i2+lookup_start(FILE)) == 38) THEN
! DEPENDS ON: readflds
    CALL readflds(nftin,2,i2,lookup(1,lookup_start(FILE)),                     &
                  len1_lookup,ice_extent,p_field,fixhd(1,FILE),                &
                  icode,cmessage)
    IF(icode /= 0)THEN
      icode=field+100
      iounit=nftin
      cmessage='REPLANCA :I/O ERROR '
      GOTO 9999
    END IF

   ELSE
     icode=field+100
     iounit=nftin
     cmessage='REPLANCA :ICE CHANGE DATA MISSING'
     GOTO 9999
   END IF
  ELSE    ! single time or LAMIPII - ie no time change field
! DEPENDS ON: readflds
    CALL readflds(nftin,1,i2,lookup(1,lookup_start(FILE)),                     &
                  len1_lookup,ice_extent,p_field,fixhd(1,FILE),                &
                  icode,cmessage)
    IF(icode /= 0)THEN
      icode=field+100
      iounit=nftin
      cmessage='REPLANCA: I/O ERROR'
      GOTO 9999
    END IF
  END IF
END IF

  IF(lsnow_depth.AND.linterpolate) THEN
IF(lookup(item_code,i2+lookup_start(FILE)) == 27) THEN

! DEPENDS ON: readflds
     CALL readflds(nftin,1,i2+1,lookup(1,lookup_start(FILE)),                  &
                   len1_lookup,snow_change,p_field,fixhd(1,FILE),              &
                   icode,cmessage)
    IF(icode /= 0)THEN
       icode=field+100
       iounit=nftin
       cmessage='REPLANCA :I/O ERROR '
       GOTO 9999
     END IF

   ELSE
     icode=field+100
     iounit=nftin
     cmessage='REPLANCA :SNOW CHANGE DATA MISSING'
     GOTO 9999
   END IF
  END IF

!  If sea surface temperature or other ice fields, read ice fraction
!  and fractional time field if not already pressent and if required
!  by time interpolation.  Similar if SLAB ref SST or ice depth needed.

  IF(field == 29.OR.(field == 28.AND.lt_int_c).OR.                             &
     field == 38)                                                              &
    THEN

   IF(.NOT.update(27)) THEN
    i3 = nlookup(27) + lookup_step(27)*step + lookup_start(                    &
       fileancil(27))
    IF ( lookup(item_code,i3)  ==  38 ) THEN

! DEPENDS ON: readflds
      CALL readflds(ftnancil(fileancil(27)),2,                                 &
                    nlookup(27)+lookup_step(27)*step,                          &
                    lookup(1,lookup_start(fileancil(27))),                     &
                    len1_lookup,ice_extent,                                    &
                    p_field,fixhd(1,fileancil(27)),                            &
                    icode,cmessage)
    IF(icode /= 0)THEN
        icode=field+100
        iounit=nftin
        cmessage='REPLANCA :I/O ERROR '
        GOTO 9999
      END IF
    IF ( rlookup(bmdi,i3-1)  /=  rmdi ) THEN
      icode = 700 + field
      WRITE(cmessage,*) 'REPLANCA: RMDI in lookup of ancillary file of times', &
                        ' of sea-ice chge not standard'
      GOTO 9999
    END IF


    ELSE
      icode=field+100
      iounit=nftin
      cmessage='REPLANCA :ICE FIELD DATA MISSING'
      GOTO 9999
    END IF
   END IF
  END IF

!  3.3 If time interpolation required, read second record

  IF(linterpolate) THEN

    i1=i2+ lookup_step(field)
    IF(i1 <= fixhd(152,FILE)) THEN

    ! If the two fields are different we have a problem with the ancillary.
      IF ( lookup(item_code,lookup_start(file)+i1-1) /= &
           lookup(item_code,lookup_start(file)+i2-1) ) THEN
        IF(l_error_ancil_struct) THEN
          icode=field+100
          cmessage='REPLANCA: start and end fields are different.'
          GO TO 9999
        END IF
      END IF

! AMIP II and ice depth don't read in ice depth field
      IF (.NOT.(lamipii.AND.lice_depth)) THEN

! DEPENDS ON: readflds
        CALL readflds(nftin,1,i1,lookup(1,lookup_start(FILE)),                 &
                      len1_lookup,ancil2,p_field,fixhd(1,FILE),                &
                    icode,cmessage)
      END IF
      IF(icode /= 0)THEN
        icode=field+300
        iounit=nftin
        cmessage='REPLANCA :I/O ERROR '
        GOTO 9999
      END IF

    ELSE !end of data on file

!   If end of data has been reached go back to the start.If data is
!   periodic.
!   Otherwise cancel time interpolation

      IF(periodic) THEN

        i1 = nlookup(field) + level - 1
      ! If the two fields are different we have a problem with the ancillary.
        IF ( lookup(item_code,lookup_start(file)+i1-1) /= &
             lookup(item_code,lookup_start(file)+i2-1) ) THEN
          IF(l_error_ancil_struct) THEN
            icode=field+100
            cmessage='REPLANCA: start and end fields are different.'
            GO TO 9999
          END IF
        END IF
! DEPENDS ON: readflds
        CALL readflds(nftin,1,i1,lookup(1,lookup_start(FILE)),                 &
                      len1_lookup,ancil2,p_field,fixhd(1,FILE),                &
                      icode,cmessage)
        IF(icode /= 0)THEN
          icode=field+300
          iounit=nftin
          cmessage='REPLANCA :I/O ERROR '
          GOTO 9999
        END IF
      ELSE
        ! We switch off time interpolation if we have reached 
        ! the end of the file.
        WRITE(6,'(A)')                                                 &
          'REPLANCA: Reached end of ancillary, switched off time interpolation.'
        linterpolate=.FALSE.
      END IF
    END IF! End of position on file test

    icode=0
  END IF ! End LINTERPOLATE

!  3.4 Perform time interpolation

  IF(linterpolate) THEN

    zero=0.0

!  Select appropriate time interpolation for each field
!  Snowdepth: set equal to zero if no snow cover

    IF(lsnow_depth) THEN
      DO i=1,p_field
        pres_value(i)=zero
      END DO

! For the call to T_INT_C, need to know BMDI is OK for SNOW_CHANGE
!  which was read in from position I2+1.
    IF ( rlookup(bmdi,lookup_start(FILE)+i2)  /=  rmdi ) THEN
      icode = 700 + field
      WRITE(cmessage,*) 'REPLANCA: RMDI in lookup of ancil file of times',     &
                        ' of snow change non-standard '
      GOTO 9999
    END IF

      CALL t_int_c (ancil1,time1,ancil2,time2,ancil_data,                      &
           time,p_field,snow_change,ancil1,pres_value)

! Ice fraction: ice depth set equal to zero if no ice

    ELSE IF(field == 27.OR.field == 29.OR.field == 38) THEN
      IF(field == 27) THEN
! For the call to T_INT_C, need to know BMDI is OK for ICE_EXTENT(1,2)
!  which was read in from position I1+1
    IF(.NOT.lamipii) THEN
    IF ( rlookup(bmdi,lookup_start(FILE)+i1)  /=  rmdi ) THEN
      icode = 700 + field
      WRITE(cmessage,*) 'REPLANCA: RMDI in lookup of ancil file of times',     &
                        ' of sea-ice chge non-standard'
      GOTO 9999
    END IF
    END IF

       IF (lamipii) THEN
! linear uncontrolled time interpolation
        CALL t_int (ice_extent,time1,ancil2,time2,ancil_data,                  &
             time,p_field)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

        DO i=1,p_field
          IF (ancil_data(i) <  0.3) ancil_data(i)=0.0
          IF (ancil_data(i) >  1.0) ancil_data(i)=1.0
        END DO

       ELSE       ! non AMIPII option
        DO i=1,p_field
          pres_value(i)=0
        END DO

        CALL t_int_c (ice_extent,time1,ancil2,time2,ancil_data,                &
             time,p_field,ice_extent(1,2),ice_extent,pres_value)

       END IF     ! end AMIPII test

      ELSE IF (field == 29.OR.field == 38) THEN

        DO i=1,p_field
          pres_value(i)=0
        END DO

        CALL t_int_c (ancil1,time1,ancil2,time2,ancil_data,                    &
             time,p_field,ice_extent(1,2),ice_extent,pres_value)


      END IF


! Sea surface temperature, set equal to TFS if ice present

    ELSE IF (field == 28.AND.lt_int_c) THEN
     IF (lamipii) THEN

      CALL t_int (ancil1,time1,ancil2,time2,ancil_data,                        &
              time,p_field)
! remove any T below TFS
      DO i=1,p_field
        IF (ancil_data(i) <  tfs)  ancil_data(i)=tfs
      END DO

     ELSE     ! non AMIPII option

      IF(.NOT.ltleads)THEN
      DO i=1,p_field
          pres_value(i)=tfs

! Set no_ice_extent indicator for controlled SST interpolation
          IF(ice_extent(i,1) == 0) THEN
            no_ice_extent(i)=1.0
          ELSE
            no_ice_extent(i)=0.0
          END IF
      END DO

      CALL t_int_c (ancil1,time1,ancil2,time2,ancil_data,                      &
           time,p_field,ice_extent(1,2),no_ice_extent,pres_value)
      ELSE
      CALL t_int (ancil1,time1,ancil2,time2,ancil_data,                        &
              time,p_field)
      END IF

     END IF   ! end AMIPII test
! Otherwise linear interpolation in time, unless missing data indicator
! present at either time.

    ELSE

! Time interpolation checks the data against the standard missing data
!   indicator - check that the field is labelled as using the same one.
!  (It is to have the right I1 here that I3 is used above.)
    IF ( rlookup(bmdi,lookup_start(FILE)+i1-1)  /=  rmdi .OR.                  &
         rlookup(bmdi,lookup_start(FILE)+i2-1)  /=  rmdi ) THEN
      WRITE (6, *) 'LOOKUPS:',                                                 &
         rlookup(bmdi,lookup_start(FILE)+i1-1),                                &
         rlookup(bmdi,lookup_start(FILE)+i2-1)
      icode = 700 + field
      cmessage = 'REPLANCA: MDI in lookup of ancil file is non-standard'
      GOTO 9999
    END IF

    LEN=p_field
!   Ozone, test for zonal mean or full field
    IF(field == 7) THEN
      IF(lookup(lbnpt,lookup_start(FILE)+i2-1) == 1) THEN
        LEN=p_rows
      END IF
!   Cariolle ozone, test for zonal mean or full field.
!   Currently same test as for conventional ozone.
    ELSE IF (field >= 178.AND.field <= 185) THEN
      IF(lookup(lbnpt,lookup_start(FILE)+i2-1) == 1) THEN
        LEN=p_rows
      END IF
    END IF

      CALL t_int(ancil1,time1,ancil2,time2,ancil_data,                         &
                 time,LEN)

    END IF ! End Lsnow_depth

! If no interpolation, copy data into final array

  ELSE ! no interpolation
   IF(lice_fraction) THEN
    IF (lamipii) THEN
    DO i=1,p_field

    ancil_data(i)=ice_extent(i,1)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

       IF (ancil_data(i) <  0.3) ancil_data(i)=0.0
       IF (ancil_data(i) >  1.0) ancil_data(i)=1.0

    END DO
    ELSE           ! non AMIP II option
      DO i=1,p_field
       ancil_data(i)=ice_extent(i,1)
      END DO
    END IF           ! end of AMIPII test
   ELSE IF (lamipii.AND.field == 28) THEN
    DO i=1,p_field
      ancil_data(i)=ancil1(i)
      IF (ancil_data(i) <  tfs) ancil_data(i)=tfs
    END DO
   ELSE
    DO i=1,p_field
      ancil_data(i)=ancil1(i)
    END DO
   END IF
  END IF !End interpolate/no interpolate

!  3.5 Updating action for each field at each level
!      Fields replaced except that Sea Surface Temperature may be
!      incremented. Take apropriate action for each field.

  IF(field <= 2.OR.field == 7.OR.field == 39.OR.field == 40                    &
  .OR.field == 41.OR.field == 42.OR.field == 43                                &
  .OR.field == 44.OR.field == 45                                               &
                                    ! multi-level murk
  .OR.(field >= 48 .AND. field <= 67 .AND. l_ukca)                             &
                                    ! single-level user ancillaries
  .OR.(field >= 68.AND.field <= 70)                                            &
                                    !NH3,soot aerosol emissions
  .OR.(field >= 72.AND.field <= 77)                                            &
                                    !Sulphur cycle
  .OR.field == 78                                                              &
                                    !CO2 EMISSIONS
  .OR.field == 82                                                              &
                                    !HADCM2 sulphate aerosol
  .OR.(field >= 90.AND.field <= 109)                                           &
                                     !multi-level user ancillaries
  .OR.(field >= 112.AND.field <= 120)                                          &
                                      !mineral dust fields
  .OR.(field >= 121.AND.field <= 122)                                          &
                                      !Biomass emissions
  .OR.field == 123                                                             &
                                     !Seawater DMS concentration
  .OR.(field >= 157.AND.field <= 177)                                          &
                                     !Aerosol climatologies
  .OR.(field >=178 .AND.field<=185)                                            &
                                     !Cariolle ozone ancillaries
  .OR.(field >= 186.AND.field <= 187)                                          &
                                     !OCFF emissions
  .OR.field == 188                                                             &
                                     !Unfiltered orography
  )THEN

!  3.5.0 Updates at all points

    LEN=p_field
!   Ozone, test for zonal mean or full field
    IF(field == 7) THEN
      IF(lookup(lbnpt,lookup_start(FILE)+i2-1) == 1) THEN
        LEN=p_rows
      END IF

!   Cariolle ozone, test for zonal mean or full field.
!   Currently same test as for conventional ozone.
    ELSE IF (field >= 178.AND.field <= 185) THEN
      IF(lookup(lbnpt,lookup_start(FILE)+i2-1) == 1) THEN
        LEN=p_rows
      END IF
    END IF

    DO i=1,LEN
      d1(d1_anciladd(field)+i-1+(level-bot_level)*LEN)=ancil_data(i)
    END DO
    
    ! Copy level 1 to level 0
    IF (bot_level == 0 .AND. level == 1) THEN
      DO i=1,LEN
        d1(d1_anciladd(field)+i-1)=ancil_data(i)
      END DO
    END IF

!  3.5.1 Updates over all land points

  ELSE IF((field > 2.AND.field < 7)                                            &
   .OR.(field > 7.AND.field < 13)                                              &
   .OR.(field > 13.AND.field < 27)                                             &
   .OR.(field == 32).OR.(field >= 35.AND.field <= 36)                          &
   .OR.(field >= 48 .AND. field <= 67 .AND. .NOT. l_ukca)                      &
                                           ! single level user ancillaries
   .OR.(field >= 46.AND.field <= 47)                                           &
                                           !Orographic roughness
   .OR.(field >= 155.AND.field <= 156)                                         &
                                           ! Orographic X & Y gradients
   .OR.(field == 79)                                                           &
                                           !Clapp-Hornberger soil parameter
   .OR.(field >= 83.AND.field <= 89)                                           &
                                           !MOSES-II 
   .OR.(field >= 33.AND.field <= 34)                                           &
                                           ! LSH Topographic index fields
   .OR.(field == 193)                                                          &
                                           ! Snow depth on tiles.
   .OR.(field >= 194.AND.field <= 196)                                         &
                                           ! Land surface albedos obs/clim 
  ) THEN


!  If not reconfiguration, set snowdepth values at all land points
!  Reset TSTAR to TM if snow cover present

    IF(lsnow_depth) THEN
      DO i=1,p_field
        IF(land(i)) THEN
          d1(d1_anciladd(field)+i-1)=ancil_data(i)
          IF(tstar_land(i) >  tm.AND.ancil_data(i) >  0.0) THEN
            tstar_land(i)=tm
          END IF
        END IF
      END DO

!  Set all other fields , which are stored at land points only

    ELSE
!  If field is a single time soil moisture, only update if
!  model time matches data time and then deactivate to prevent
!  any further updating

      IF(field == 36.AND.fixhd(10,FILE) == 0) THEN
        IF(lookup(lbyr,lookup_start(FILE)) == i_year.AND.                      &
           lookup(lbmon,lookup_start(FILE)) == i_month.AND.                    &
           lookup(lbdat,lookup_start(FILE)) == i_day.AND.                      &
           lookup(lbhr,lookup_start(FILE)) == i_hour) THEN
          WRITE(6,*) 'Updating soil moisture at ',                             &
                       i_year,i_month,i_day,i_hour,                            &
          ' for level ',rlookup(blev,lookup_start(FILE)+level-1)
          dz_soil(level)=rlookup(blev,lookup_start(FILE)+level-1)

          CALL compress_to_mask(ancil_data,d1(d1_anciladd(field)+ &
              (level-1)*land_field),land,p_field,i)
! Switch off to prevent further attempts to update
          fieldcode(1,field)=0
          steps(field)=0
! Set flag to indicate that soil moisture has been updated
          smc_updated=.TRUE.
        ELSE
          WRITE(6,*) 'Update of soil moisture skipped'
        END IF

      ELSE
! other fields
        CALL compress_to_mask(ancil_data,d1(d1_anciladd(field)+  &
            (level-1)*land_field),land,p_field,i)
      END IF

    END IF



! Iceberg calving for the OASIS coupler:
  ELSE IF (field == 13) THEN
    DO i=1,u_field
      d1(d1_anciladd(field)+i-1)=ancil_data(i)
    END DO

!  3.5.2 Ice fraction
  ELSE IF(field == 27) THEN
    DO i=1,p_field
      ice_fraction(i)=0.
      IF (sea(i)) THEN
        ice_fraction(i)=ancil_data(i)
      END IF
    END DO

! Reduce TSTAR to TFS where ice fraction greater than zero
! Required at present because radiation and boundary layer codes
! assume T* is TFS and ignore any value set in TSTAR.

    IF(.NOT.ltleads)THEN
      DO i=1,p_field
        IF(ice_fraction(i) >  0.0) THEN
          tstar_ssi(i)=MIN(tstar_ssi(i),tfs)
        END IF
      END DO
    END IF

!  3.5.3 Sea surface temperatures for atmosphere, adds anomaly field
!        (if defined) to updated ancillary where no sea-ice present,
!        otherwise updates field directly from ancillary

  ELSE IF (field == 28) THEN


    DO i=1,p_field
      IF (l_ctile.OR.ice_fraction(i) <= 0.0) THEN
        IF (sea(i)) THEN
          IF (l_sstanom) THEN
            tstar_sea(i)=ancil_data(i)+tstar_anom(i)
          ELSE
            tstar_sea(i)=ancil_data(i)
          END IF
          IF (ice_fraction(i) <= 0.0) tstar_ssi(i)=tstar_sea(i)
        END IF
      END IF
    END DO

!  3.5.3.1 Reference SSTs for SLAB model

  ELSE IF (field == 37) THEN

    DO i=1,p_field
      IF (sea(i)) THEN
        d1(d1_anciladd(field)+i-1)=ancil_data(i)
      ELSE
        d1(d1_anciladd(field)+i-1)=rmdi

      END IF
    END DO

!  3.5.4 Sea ice thickness/Reference seaice thickness for SLAB
!        Update over all sea points (all sea ice points are the only
!        ones strictly required, but this cannot be determined easily)

  ELSE IF (field == 29.OR.field == 38) THEN

    DO i=1,p_field
      IF(sea(i)) THEN
        d1(d1_anciladd(field)+i-1)=ancil_data(i)
      END IF
    END DO

!  3.5.5 Surface currents

  ELSE IF (field == 30) THEN
    DO i=1,u_field
      d1(d1_anciladd(field)+i-1)=ancil_data(i)
    END DO

  ELSE IF (field == 31) THEN
    DO i=1,v_field
      d1(d1_anciladd(field)+i-1)=ancil_data(i)
    END DO

  ELSE IF (field >= 124.AND.field <= 126) THEN   ! Riv Routing
    icode=750
    cmessage='REPLANCA: ERROR Trying to use Riv Route Ancils'
!         There is no code yet to support in model updating of the
!         River Routing fields.
    GOTO 9999

  ELSE

  WRITE(6,*)' REPLANCA: ERROR - FIELD ',field,                                 &
  ' omitted from update block'

  END IF !End tests on FIELD numbers

!  End loop over levels

i2=i2+1

  END DO

!  End loop over ancillary fields (atmosphere)
 END IF ! LAMIPII and ice depth test

END IF    ! End UPDATE(field) test     level 1 IF


END DO


IF(l_ctile)THEN
  DO i=1,p_field
    IF(sea(i).AND.ice_fraction(i) >  0.0)THEN
      IF(ltleads.OR.lamipii)THEN

        tstar_ssi(i)=ice_fraction(i)*tstar_sice(i)                             &
          +(1.-ice_fraction(i))*tstar_sea(i)

      ELSE

        tstar_sea(i)=tfs
        tstar_sice(i)=(tstar_ssi(i)                                            &
          -(1.-ice_fraction(i))*tstar_sea(i))/ice_fraction(i)

      END IF
    END IF

    tstar(i)=fland_g(i)*tstar_land(i)                                          &
      +(1.-fland_g(i))*tstar_ssi(i)
  END DO
ELSE
  DO i=1,p_field
    IF(land(i))THEN
      tstar(i)=tstar_land(i)
    ELSE
      tstar(i)=tstar_ssi(i)
    END IF
  END DO
END IF

!     Set up surface temperatures:
IF(l_ctile)THEN
  DO i=1,p_field
    tstar_land_ctile(i)=tstar_land(i)
    tstar_sea_ctile(i)=tstar_sea(i)
    ! The use of TSTAR_SICE appears to cause problems in 
    ! some configurations (e.g. seasonal). Possibly because
    ! of the use of inconsistent ancillary fields.
    ! Hence this is commented out but retained for reference.
    ! [See also equivalent change in replanca-rcf_replanca.F90]
    !tstar_sice_ctile(i)=tstar_sice(i)
  END DO
END IF

9999 CONTINUE
    
IF (lhook) CALL dr_hook('REPLANCA',zhook_out,zhook_handle)
RETURN

END SUBROUTINE replanca
