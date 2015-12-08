! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Packing Tools

MODULE Mask_Compression

IMPLICIT NONE
   
CONTAINS

!    Subroutine Compress_To_Mask --------------------------------------
!    Purpose:  Selects land points values from input horizontal field 
!              and writes then as contiguous values to array
!              COMPRESSEDFIELD.
!    ------------------------------------------------------------------

  SUBROUTINE Compress_To_Mask &
      (field,compressedField,mask,points,maskPoints)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)  :: points          !Total no of points (uncompressed)
    INTEGER, INTENT(OUT) :: maskPoints      !No of land points
    LOGICAL, INTENT(IN)  :: mask(points)    ! Land-sea mask
    REAL, INTENT(OUT)    :: compressedField(points) !Data (compressed)
    REAL, INTENT(IN)     :: field(points)   !Data (uncompressed)
    
    ! Local varables
    INTEGER              :: index(points)   !Gather index
    INTEGER              :: i 
    
    !  Compute gather index for land points    
    maskPoints = 0
    DO i=1,points
      IF (mask(i)) THEN
        maskPoints=maskPoints + 1
        index(maskPoints) = i
      END IF
    END DO
    
    !  Gather land points from input array field
    DO i=1,maskPoints
      compressedField(i)=field(index(i))
    END DO
    
    RETURN
  END SUBROUTINE Compress_To_Mask
  
  !    Subroutine Expand_From_Mask --------------------------------------
  !    Purpose:  Selects successive values from input array
  !              COMPRESSEDFIELD and writes them to masked points
  !              in horizontal field array FIELD.
  !    ------------------------------------------------------------------

  SUBROUTINE Expand_From_Mask &
      (field,compressedField,mask,points,maskPoints)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)  :: points          ! Total no of points (expanded)
    INTEGER, INTENT(OUT) :: maskPoints      ! No of points in mask
    LOGICAL, INTENT(IN)  :: mask(points)    ! mask
    REAL, INTENT(IN)     :: compressedField(points) ! Data on land points only
    REAL, INTENT(OUT)    :: field(points)   ! expanded data
    
    ! Local varables
    INTEGER              :: index(points)   !Scatter index for land points
    INTEGER              :: i 
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

    !  Initialise sea points to MDI
    !  Calculate scatter index for land points
    maskPoints = 0
    DO i=1,points
      field(i)=rmdi
      IF (mask(i)) THEN
        maskPoints=maskPoints + 1
        index(maskPoints) = i
      END IF
    END DO

    !  Scatter points to array field
!DIR$ IVDEP
    DO i=1,maskPoints
      field(index(i))=compressedField(i)
    END DO

    RETURN
  END SUBROUTINE Expand_From_Mask
END MODULE Mask_Compression
