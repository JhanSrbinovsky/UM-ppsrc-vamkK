! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Establish what information is contained within the GRIB file

MODULE Rcf_Grib_Spcl_Ctl_Mod

! SUBROUTINE Rcf_Grib_Spcl_Ctl
!
! Description: This routine holds the 'logic' which results in calls
!              to further routines which perform 'minor' transformations
!              on the data read from the GRIB file.
!              e.g. to transform geopotential to Orography by dividing
!                   by the  gravitational acceleration g.
!
! Method: The logic is split into three blocks -
!         Block 1:
!            Simple transformations which affect only a single field,
!            not requiring any other fields to exist. Thus the 'switch'
!            _will_ be the STASH code.
!            e.g. The transformation of geopotential to Orography
!
!         Block 2:
!            Transformations which may apply to _all_ fields based on a
!            given criteria rather than the STASH code.
!            e.g. The reversing of N->S grids to get S->N grids
!
!         Block 3:
!            More complex transformations which will rely on the
!            existance of multiple fields within the dump.
!            e.g. converting T to TH requires Pstar
!            **NOTE** Any fields 'created' here (i.e. not _replacing_ a
!            particular field read in from the GRIB data) Should also
!            have a corresponding entry in Rcf_Grib_Spcl_Hdr to allocate
!            header space for the field.
!            -Note- This is done last specifically so that any
!            calculations requiring previously read fields have data
!            in correct order (i.e. already through Block 2).
!
!         REMEMBER - all 3 blocks are passed through _every_ time this
!                    routine is called. (i.e. for every field read)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!

CONTAINS
  SUBROUTINE Rcf_Grib_Spcl_Ctl(FieldData,L_FldData,Lists,Current,              &
      Hdr_Dmy,Hdr_Itm,List_no,Lkps_Posn)

    USE Rcf_GRIB_Block_Params_Mod, ONLY :                                      &
        List_Marker,     Grib_Record,                                          &
        p_LatGrdPnt1,    p_LatExtrmPt,                                         &
        p_Param_ID,      LenArrayMax,                                          &
        Grb_Data_Log,    Grb_Data_Int,                                         &
        EID_Surf_Press,  EID_Log_Surf_Press

    USE rcf_GRIB_lookups_Mod, ONLY :                                           &
        grib_max_fields,                                                       &
        grborigecmwf

    USE Rcf_HeadAddress_Mod, ONLY :                                            &
        ic_ylen,IC_Xlen

    USE Rcf_UMhead_Mod, ONLY :                                                 &
        um_header_type            ! Derived containing UM header info

    USE PrintStatus_mod, ONLY :                                                &
        PrintStatus,                                                           &
        PrStatus_Diag

    USE Rcf_StashCodes_Mod, ONLY :                                             &
        stashcode_orog,          stashcode_exner,                              &
        stashcode_lsm,                                                         &
        stashcode_pstar,         stashcode_soil_temp,                          &
        stashcode_soil_moist,    stashcode_mean_snow

    USE UM_ParVars, ONLY :                                                     &
        mype

    USE Rcf_Grib_Spcl_Orog_Mod, ONLY :                                         &
        Rcf_Grib_Spcl_Orog

    USE Rcf_Grib_Spcl_LSM_Mod, ONLY :                                          &
        Rcf_Grib_Spcl_LSM

    USE Rcf_Grib_Spcl_Exner_Mod, ONLY :                                        &
        Rcf_Grib_Spcl_Exner

    USE Rcf_Grib_Spcl_LPstar_Mod, ONLY :                                       &
        Rcf_Grib_Spcl_LPstar

    USE Rcf_Grib_Spcl_SoilM_Mod, ONLY :                                        &
        Rcf_Grib_Spcl_SoilM

    USE Rcf_Grib_Spcl_Snow_Mod, ONLY :                                         &
        Rcf_Grib_Spcl_Snow

    USE EReport_Mod, ONLY :                                                    &
        EReport

    USE Rcf_Reverse_Field_Mod, ONLY :                                          &
        Rcf_Reverse_Field

    USE earth_constants_mod, ONLY: g

    USE mask_compression, ONLY: compress_to_mask

    USE Lookup_addresses

    IMPLICIT NONE

! Declarations:
!   These are of the form:-
!     INTEGER, INTENT(IN) :: ExampleVariable    !Description of variable
!
! Global variables (#include statements etc):

! Subroutine arguments

!< Scalar arguments with intent(In):>
    INTEGER, INTENT(IN)                 :: Lkps_Posn
    INTEGER, INTENT(IN)                 :: List_no

!< Array  arguments with intent(In):>
    TYPE (List_Marker), INTENT(IN)      :: Lists(0:grib_max_fields)
    TYPE (Um_Header_type),INTENT(IN)    :: Hdr_Dmy

!< Scalar arguments with intent(InOut):>

!< Array  arguments with intent(InOut):>
    REAL, INTENT(INOUT)                 :: FieldData(LenArrayMax)
    LOGICAL, INTENT(INOUT)              :: L_FldData(LenArrayMax)
    TYPE (Grib_Record),POINTER          :: Current
    TYPE (Um_Header_type),INTENT(INOUT) :: Hdr_Itm

!< Scalar arguments with intent(out):>

!< Array  arguments with intent(out):>

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
! contains LBLREC (amongst others)

! Local constants
    CHARACTER (LEN=*), PARAMETER     :: RoutineName='Rcf_GRIB_Spcl_Ctl'

! Local variables
    CHARACTER (LEN=80)               :: Cmessage(2)   ! used for EReport
    INTEGER                          :: ErrorStatus   ! used for EReport

    INTEGER                          :: i
    INTEGER                          :: Criteria

    REAL                             :: LandPacked(LenArrayMax)
    INTEGER                          :: LandPoints

! Local variables for saved land sea mask
! Note: although L_GotLSM is initialised as false, it acts as a saved
! variable from then on.
    LOGICAL,SAVE  :: LandSeaMask(LenArrayMax)
    LOGICAL,SAVE  :: l_gotlsm=.FALSE.

!=======================================================================
!  Routine Code Start :
!=======================================================================

!=======================================================================
!  Block 1 : Simple transformations
!=======================================================================

    SELECT CASE (Current % Stashcode)

    CASE (stashcode_orog)
      IF ( mype == 0 ) THEN
        WRITE(6,'(A)') "Special treatment for orography"
      END IF
      CALL Rcf_Grib_Spcl_Orog(Current,FieldData)

    CASE (stashcode_lsm)
      IF ( mype == 0 ) THEN
        WRITE(6,'(A)') "Special treatment for landsea mask"
      END IF
      CALL Rcf_Grib_Spcl_LSM(Current,FieldData,L_FldData,Hdr_Itm,              &
          Lkps_Posn)

    CASE (stashcode_exner)
      IF ( mype == 0 ) THEN
        WRITE(6,'(A)') "Special treatment for exner pressure"
      END IF
      CALL Rcf_Grib_Spcl_Exner(Current,FieldData)

    CASE (stashcode_pstar)
      IF (Current % Block_1(p_Param_ID) == EID_Log_Surf_Press) THEN
        IF ( mype == 0 ) THEN
          WRITE(6,'(A)') "Special treatment for log surface pressure"
        END IF
        CALL Rcf_Grib_Spcl_LPstar(Current,FieldData)
      END IF

    CASE (stashcode_soil_moist)
      IF ( mype == 0 ) THEN
        WRITE(6,'(A)') "Special treatment for soil moisture"
      END IF
      CALL Rcf_Grib_Spcl_SoilM(Current,FieldData)

    CASE (stashcode_mean_snow)
      IF ( mype == 0 ) THEN
        WRITE(6,'(A)') "Special treatment for snow amount"
      END IF
      CALL Rcf_Grib_Spcl_Snow(Current,FieldData)

    END SELECT

!=======================================================================
!  Block 2 : Global transformations
!=======================================================================

! Check data is stored South to North not North to South
    IF ( Current % Block_2(p_LatGrdPnt1) >                                     &
        Current % Block_2(p_LatExtrmPt)) THEN
      IF ( PrintStatus >= PrStatus_Diag   ) THEN
        IF ( mype == 0 ) THEN
          WRITE (6,'(2A)') "About to reverse latitudes of ", Current % Desc
        END IF
      END IF

      SELECT CASE(Current % Data_Type)

      CASE (Grb_Data_Log)
        CALL Rcf_Reverse_Field(L_FldData,     &     ! Data to reverse
            Hdr_Itm % IntC (ic_xlen), & ! Row Length
            Hdr_Itm % IntC (ic_ylen), & ! No. of rows
            1,             &     ! No. of levels in Data
            Lkps_Posn,     &     ! Position in lookups
            Hdr_Itm)             ! UM Hdr for updating

      CASE (Grb_Data_Int)
    ! At present there is no code to flip an Integer data field
        Cmessage(1) = 'Tried to Flip Integer Data Type without any code'
        ErrorStatus = 10
        CALL EReport( RoutineName, ErrorStatus, Cmessage(1) )


      CASE Default
        CALL Rcf_Reverse_Field(FieldData,     &     ! Data to reverse
            Hdr_Itm % IntC (ic_xlen), & ! Row Length
            Hdr_Itm % IntC (ic_ylen), & ! No. of rows
            1,             &     ! No. of levels in Data
            Lkps_Posn,     &     ! Position in lookups
            Hdr_Itm)             ! UM Hdr for updating

      END SELECT
    END IF

!=======================================================================
!  Block 3 : Complex transformations
!=======================================================================

! If data needs to be land packed, then do so.
! Must be done after any other field manipulations.
! - Will be dependent on stash code
    SELECT CASE (Current % Stashcode)

    CASE (stashcode_lsm)
    ! Store land/sea mask in case needed for land packing below
      IF ( mype == 0 ) THEN
        WRITE(6,'(A)') "Store Land/Sea Mask for landpacking"
      END IF
      LandSeaMask=L_FldData
      l_gotlsm=.TRUE.

    CASE (stashcode_soil_temp,                                                 &
        stashcode_soil_moist)

      IF (.NOT. l_gotlsm) THEN
        Cmessage(1) = 'LSM not available for land packing'
        ErrorStatus = 10
        CALL Ereport( RoutineName, ErrorStatus, Cmessage(1) )
      END IF

      IF ( mype == 0 ) THEN
        WRITE(6,'(A,I5)') "Performing land-packing for stashcode ",            &
            Current % Stashcode
      END IF

    ! Use main UM routine compress_to_mask to do the land packing.
      CALL compress_to_mask(FieldData,    & ! Data to landpack
          LandPacked,   & ! Landpacked data
          LandSeaMask,  & ! land/sea mask
          LenArrayMax,  & ! size of input/output arrays
          LandPoints)     ! number of landpoints

    ! Put the land packed values at the head of the FieldData
    ! array so that the reconfiguration sees the correct values.
      FieldData(1:LandPoints)=LandPacked(1:LandPoints)

    END SELECT

!=======================================================================
!  Cleaning Up afterwards
!=======================================================================

    RETURN

  END SUBROUTINE Rcf_Grib_Spcl_Ctl
END MODULE Rcf_Grib_Spcl_Ctl_Mod
