! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Convert GRIB LSM to Logical UM one
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!

MODULE Rcf_Grib_Spcl_LSM_Mod

! SUBROUTINE Rcf_Grib_Spcl_LSM
!
! Description: Converts LSM from 'real' type to 'logical'
!
! Method: flip flop the logical based on data > 0.9999 for .True.

CONTAINS
  SUBROUTINE Rcf_Grib_Spcl_LSM(Current,FpData,LgData,UM_Hdr,Marker)

    USE Rcf_GRIB_Block_Params_Mod, ONLY :                                      &
        Grib_Record,                                                           &
        LenArrayMax,                                                           &
        Grb_Data_Real,                                                         &
        Grb_Data_Log

    USE PrintStatus_mod, ONLY :                                                &
        PrintStatus,                                                           &
        PrStatus_Diag

    USE Rcf_UMhead_Mod, ONLY :                                                 &
        um_header_type            ! Derived containing UM header info

    USE Rcf_HeadAddress_Mod, ONLY :                                            &
        IC_NumLandPoints

    USE EReport_Mod, ONLY :                                                    &
        EReport

    USE UM_ParVars, ONLY :                                                     &
        mype

    USE lookup_addresses

    IMPLICIT NONE

! Declarations:
! Global variables (#include statements etc):
! contains LBLREC (amongst others)

! Subroutine arguments

!< Scalar arguments with intent(In):>
    INTEGER                             :: Marker  ! position within lookups

!< Array  arguments with intent(InOut):>
    TYPE (Grib_Record),POINTER          :: Current
    REAL, INTENT(INOUT)                 :: FpData(LenArrayMax)
    LOGICAL, INTENT(INOUT)              :: LgData(LenArrayMax)
    TYPE (Um_Header_type)               :: UM_Hdr

! Local constants
    CHARACTER (LEN=*), PARAMETER     :: RoutineName='Rcf_GRIB_Spcl_LSM'

! Local variables

    CHARACTER (LEN=80)               :: Cmessage(2)   ! used for EReport
    INTEGER                          :: ErrorStatus   ! used for EReport

    INTEGER                          :: cnter,i       ! used to count no.
                                                  ! of land points

!=======================================================================
!  Routine Code Start :
!=======================================================================

    cnter = 0

! Double Check the datatype is still 'real'
    IF (Current % Data_Type == Grb_Data_Real ) THEN

      IF ( PrintStatus >= PrStatus_Diag  ) THEN
        IF ( mype == 0 ) THEN
          WRITE (6,'(A)') 'Converting Real LSM to Logical LSM'
        END IF
      END IF

      WHERE (FpData >0.99999)
      LgData = .TRUE.  ! It's a land point
      ElseWhere
      LgData = .FALSE. ! It's not
      END WHERE

      Current % Data_Type = Grb_Data_Log

    ELSE   ! It wasn't real data I recieved.
      Cmessage(1) = 'Couldnt convert data type recieved to logical'
      ErrorStatus = 10
      CALL EReport( RoutineName, ErrorStatus, Cmessage(1) )

    END IF

    DO i = 1, UM_Hdr % Lookup (lblrec,Marker)
      IF (LgData(i)) THEN
        cnter =cnter + 1
      END IF
    END DO

! record the no. of land points in the header
    UM_Hdr % IntC( IC_NumLandPoints ) = cnter

    RETURN

  END SUBROUTINE Rcf_Grib_Spcl_LSM
END MODULE Rcf_Grib_Spcl_LSM_Mod
