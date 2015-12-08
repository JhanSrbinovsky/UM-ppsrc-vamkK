! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

!+ reverse a N-S grid to give a S-N one

MODULE Rcf_Reverse_Field_Mod

! Description:
!              Reverse the position of "n_rows" rows for each level of
! the 3D data array passed into the routine such that the first row
! swaps place with the row at position n_rows. The 2nd row swaps place
! with th row at position n_rows-1 and so on. The lookup headers for
! the field are altered to be consistent with the new row ordering.
! Based on the routine - PF_REVERSE at UM4.5
!
! Method:
!
! -Note- This routine has been overloaded to handle both real and
! logical fields. Please ensure you duplicate 'fixes' performed in one
! routine in the other routine where necessary

  INTERFACE Rcf_Reverse_Field
    MODULE PROCEDURE Rcf_Reverse_Field_Log, Rcf_Reverse_Field_Real
  END INTERFACE

CONTAINS

  SUBROUTINE Rcf_Reverse_Field_Real(RData,row_length,n_rows,levels,            &
      field_no,UM_Hdr)

    USE Rcf_UMhead_Mod, ONLY :                                                 &
        UM_Header_Type            ! Derived containing UM header info

    USE Rcf_HeadAddress_Mod, ONLY  :                                           &
        RC_FirstLat,                                                           &
        RC_LatSpacing

    USE EReport_Mod, ONLY :                                                    &
        EReport

    USE lookup_addresses
    IMPLICIT NONE
! Subroutine arguments

!< Scalar arguments with intent(in):>
    INTEGER, INTENT(IN)              :: row_length  ! length of each row
    INTEGER, INTENT(IN)              :: n_rows      ! no. of rows per level
    INTEGER, INTENT(IN)              :: levels      ! no. of levels
    INTEGER, INTENT(IN)              :: field_no    ! pos in lookup of
                                                ! first field

!< Array  arguments with intent(InOut):>
    REAL,    INTENT(INOUT)              :: RData(*) ! the data to be reversd
    TYPE (UM_Header_Type),INTENT(INOUT) :: UM_Hdr   ! UM header info

! Comdecks
! defines BDY and BZY

! Local variables

    INTEGER                          :: Int_Val
    INTEGER                          :: top_row
    INTEGER                          :: bottom_row
    INTEGER                          :: point
    INTEGER                          :: lev_shift
    INTEGER                          :: b_row_shift
    INTEGER                          :: t_row_shift
    INTEGER                          :: point_loc
    INTEGER                          :: k
    REAL                             :: temp
    REAL                             :: Zeroth_Lat
    REAL                             :: Lat_Spacing
    REAL                             :: Real_Val

    INTEGER                          :: ErrorStatus
    CHARACTER(LEN=80)                :: cMessage
    CHARACTER (LEN=*), PARAMETER     :: RoutineName='Rcf_Reverse_Field_Real'

!=======================================================================
! Subroutine Code.
!=======================================================================
!loop over levels
    DO k=1,levels

      lev_shift = (k-1) * (row_length * n_rows)

      DO top_row=1,n_rows/2

        bottom_row = n_rows + 1 - top_row

        b_row_shift = lev_shift + ((bottom_row - 1) * row_length)
        t_row_shift = lev_shift + ((top_row    - 1) * row_length)

        DO point=1,row_length

      ! Swap over corresponding point at top and bottom of data.
          temp                      = RData(b_row_shift + point)
          RData(b_row_shift + point) = RData(t_row_shift + point)
          RData(t_row_shift + point) = temp

        END DO ! point

      END DO ! top_row

  ! Now update the Latitude info held in Lookups
      Lat_Spacing = TRANSFER( UM_Hdr % Lookup( bdy, field_no ) , Real_Val)
      Zeroth_Lat  = TRANSFER( UM_Hdr % Lookup( bzy, field_no ) , Real_Val)

  ! Transpose Zeroth latitude.
      Zeroth_Lat  = Zeroth_Lat + ( (n_rows + 1) * Lat_Spacing )

  ! Sign of the latitude interval is changed.
      Lat_Spacing = - Lat_Spacing


  ! Put new values back into the Lookups
      UM_Hdr % Lookup( bdy, field_no ) = TRANSFER(Lat_Spacing , Int_Val)
      UM_Hdr % Lookup( bzy, field_no ) = TRANSFER(Zeroth_Lat  , Int_Val)

  ! double check the values in the Real Headers
      IF (Um_Hdr % RealC (RC_FirstLat) /= ( Zeroth_Lat + Lat_Spacing)) THEN
        Um_Hdr % RealC (RC_FirstLat) = ( Zeroth_Lat + Lat_Spacing)
      END IF
      IF (Um_Hdr % RealC (RC_LatSpacing) /= ABS(Lat_Spacing) ) THEN
    !Um_Hdr % RealC (RC_LatSpacing) = Abs(Lat_Spacing)
        WRITE(cMessage,'(A)') "Latitude Spacing appears to have changed"
        ErrorStatus = 20
        CALL EReport( RoutineName, ErrorStatus, Cmessage)

      END IF


  ! ***************************************************
  ! At UM4.5 something was done here to re-align u or v
  ! onto the correct spacings.
  ! ***************************************************

    END DO ! k over levels


    RETURN

  END SUBROUTINE Rcf_Reverse_Field_Real


  SUBROUTINE Rcf_Reverse_Field_Log(LData,row_length,n_rows,levels,             &
      field_no,UM_Hdr)

    USE Rcf_UMhead_Mod, ONLY :                                                 &
        UM_Header_Type            ! Derived containing UM header info

    USE Rcf_HeadAddress_Mod, ONLY  :                                           &
        RC_FirstLat,                                                           &
        RC_LatSpacing

    USE EReport_Mod, ONLY :                                                    &
        EReport

    USE ereport_mod, ONLY : ereport
    USE lookup_addresses
    IMPLICIT NONE
! Subroutine arguments

!< Scalar arguments with intent(in):>
    INTEGER, INTENT(IN)              :: row_length  ! length of each row
    INTEGER, INTENT(IN)              :: n_rows      ! no. of rows per level
    INTEGER, INTENT(IN)              :: levels      ! no. of levels
    INTEGER, INTENT(IN)              :: field_no    ! pos in lookup of
                                                ! first field

!< Array  arguments with intent(InOut):>
    LOGICAL, INTENT(INOUT)              :: LData(*) ! the data to be reversd
    TYPE (UM_Header_Type),INTENT(INOUT) :: UM_Hdr   ! UM header info

! Comdecks
! defines BDY and BZY

! Local variables

    INTEGER                          :: Int_Val
    INTEGER                          :: top_row
    INTEGER                          :: bottom_row
    INTEGER                          :: point
    INTEGER                          :: lev_shift
    INTEGER                          :: b_row_shift
    INTEGER                          :: t_row_shift
    INTEGER                          :: point_loc
    INTEGER                          :: k
    REAL                             :: Zeroth_Lat
    REAL                             :: Lat_Spacing
    REAL                             :: Real_Val
    LOGICAL                          :: temp

    INTEGER                          :: ErrorStatus
    CHARACTER(LEN=80)                :: cMessage
    CHARACTER (LEN=*), PARAMETER     :: RoutineName='Rcf_Reverse_Field_Log'

!=======================================================================
! Subroutine Code.
!=======================================================================
!loop over levels
    DO k=1,levels

      lev_shift = (k-1) * (row_length * n_rows)

      DO top_row=1,n_rows/2

        bottom_row = n_rows + 1 - top_row

        b_row_shift = lev_shift + ((bottom_row - 1) * row_length)
        t_row_shift = lev_shift + ((top_row    - 1) * row_length)

        DO point=1,row_length

      ! Swap over corresponding point at top and bottom of data.
          temp                      = LData(b_row_shift + point)
          LData(b_row_shift + point) = LData(t_row_shift + point)
          LData(t_row_shift + point) = temp

        END DO ! point

      END DO ! top_row

  ! Now update the Latitude info held in Lookups
      Lat_Spacing = TRANSFER( UM_Hdr % Lookup( bdy, field_no ) , Real_Val)
      Zeroth_Lat  = TRANSFER( UM_Hdr % Lookup( bzy, field_no ) , Real_Val)

  ! Transpose Zeroth latitude.
      Zeroth_Lat  = Zeroth_Lat + ( (n_rows + 1) * Lat_Spacing )

  ! Sign of the latitude interval is changed.
      Lat_Spacing = - Lat_Spacing


  ! Put new values back into the Lookups
      UM_Hdr % Lookup( bdy, field_no ) = TRANSFER(Lat_Spacing , Int_Val)
      UM_Hdr % Lookup( bzy, field_no ) = TRANSFER(Zeroth_Lat  , Int_Val)

  ! double check the values in the Real Headers
      IF (Um_Hdr % RealC (RC_FirstLat) /= ( Zeroth_Lat + Lat_Spacing)) THEN
        Um_Hdr % RealC (RC_FirstLat) = ( Zeroth_Lat + Lat_Spacing)
      END IF
      IF (Um_Hdr % RealC (RC_LatSpacing) /= ABS(Lat_Spacing) ) THEN
    !Um_Hdr % RealC (RC_LatSpacing) = Abs(Lat_Spacing)
        WRITE(cMessage,'(A)') "Latitude Spacing appears to have changed"
        ErrorStatus = 20
        CALL EReport( RoutineName, ErrorStatus, Cmessage)

      END IF


  ! ***************************************************
  ! At UM4.5 something was done here to re-align u or v
  ! onto the correct spacings.
  ! ***************************************************

    END DO ! k over levels

    RETURN

  END SUBROUTINE Rcf_Reverse_Field_Log

END MODULE Rcf_Reverse_Field_Mod
