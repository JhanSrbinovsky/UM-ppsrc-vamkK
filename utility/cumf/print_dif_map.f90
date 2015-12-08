! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Program  MAIN_COMPARE and Subroutine COMPARE
!LL
!LL  Purpose: Compares two UM atmosphere, ocean, or ancillary files.
!LL           MAIN_COMPARE reads in fixed length and integer
!LL           headers of UM files to be compared, extracts dimensions
!LL           of each file and then passes these values to
!LL           subroutine COMPARE.
!LL
!LL            COMPARE subroutine:
!LL          Compares two UM atmosphere, ocean, or ancillary files.
!LL          COMPARE reads in headers and data fields from files on
!LL          NFTIN1 and NFTIN2, comparing values.
!LL          UNIT 6: If an exact compare is found the message 'OK'
!LL          is written out, otherwise
!LL          i)  if header, all differring values are printed
!LL          ii) if field, 1st 10 differring values are printed plus
!LL              the maximum difference between the fields.
!LL          iii) if field only present in one file, a warning message
!LL               is displayed
!LL          UNIT 7: Number of differences displayed for each header.
!LL                  Number of fields with differences is also
!LL                  displayed along with the number of differences
!LL                  for each field which has differences
!LL
!LL  Programming standard:
!LL
!LL  Logical components covered:
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  -----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs

      SUBROUTINE PRINT_DIF_MAP(DIFF,ROWS,COLS,KEY,                      &
     &  GRID_TYPE1,GRID_TYPE2)
! Declarations
!LL Writes out a map of the differences between two fields - with one
!LL character per point. This allows points in two fields which are
!LL different to be quickly identified.
!LL Writes to UNIT10 - opened in COMPARE - filename must be supplied
!LL by UNIT10 environment variable via the cumf script

        USE mask_compression, ONLY: expand_from_mask
      USE UM_ParVars
      IMPLICIT NONE

      INTEGER                                                           &
     &  ROWS                                                            &
              ! IN : number of rows in field
     &, COLS                                                            &
              ! IN : number of cols in field
     &, GRID_TYPE1                                                      &
                      !IN grid type for file 1
     &, GRID_TYPE2    !IN grid type for file 2

      CHARACTER(LEN=1)                                                       &
     &  DIFF(ROWS*COLS)  ! IN : difference map field to be output

      CHARACTER(LEN=*)                                                     &
     &  KEY  ! IN : key to difference map

! Local variables
      INTEGER X,Y,Z
      integer i ,j                                                      &
     &, LAND_POINTS                                                     &
                      !no of land points
     &, POINTS        !total no of land & sea points to be processed

! Constants from comdecks:----------------------
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
!====================== COMDECK ATM_LSM ========================
! Description:
!   This comdeck contains a COMMON block which contains the
!   atmosphere land sea mask - both the full field, and the
!   local subdomain on this processor.
!   This data is required for various compression/decompression
!   algorithms.
!
!   Requires AMAXSIZE comdeck to be called first for Max2DFieldSize
!

      LOGICAL                                                           &
!  Full-grid land-sea mask:
     &  atmos_landmask(Max2DFieldSize)                                  &
! Local subdomain area land-sea mask:
     &, atmos_landmask_local(Max2DFieldSize)

      INTEGER atmos_number_of_landpts ! total number of land points

      COMMON /Atmos_LSM_Common/                                         &
     &  atmos_landmask                                                  &
     &, atmos_landmask_local                                            &
     &, atmos_number_of_landpts

! End of comdeck ATM_LSM

      REAL WORK_DIFF(ROWS*COLS),NWORK_DIFF(ROWS*COLS)
      CHARACTER(LEN=1) numb(10), blank

      data numb/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
      data blank /' '/

      WRITE(10,'(/a/)') KEY


      POINTS=ROWS*COLS

!L 1. Converting symbols to real numbers (input to expand_from_mask)
!  -----------------------------------------------------------------
          IF((GRID_TYPE1 == 21).AND.(GRID_TYPE2 == 21)) THEN
            do I=1,rows*cols
              if (diff(I) == "#") work_diff(I)=5.0
              if (diff(I) == "X") work_diff(I)=4.0
              if (diff(I) == "O") work_diff(I)=3.0
              if (diff(I) == "o") work_diff(I)=2.0
              if (diff(I) == ":") work_diff(I)=1.0
              if (diff(I) == ".") work_diff(I)=0.0
            enddo

! Pass work_diff to expand_from_mask
            CALL expand_from_mask(NWORK_DIFF,WORK_DIFF,ATMOS_LANDMASK,  &
     &                            POINTS,LAND_POINTS)

!L 2. Converting real numbers to symbols (to be output)
!  ----------------------------------------------------
            do I=1,points
              if (nwork_diff(I) == 5.0) diff(I)="#"
              if (nwork_diff(I) == 4.0) diff(I)="X"
              if (nwork_diff(I) == 3.0) diff(I)="O"
              if (nwork_diff(I) == 2.0) diff(I)="o"
              if (nwork_diff(I) == 1.0) diff(I)=":"
              if (nwork_diff(I) == 0.0) diff(I)="."
              if (nwork_diff(I) == rmdi) diff(I)="~"
            enddo
          ENDIF

      write(10,'(6x,120a1)') ((blank, j=1,9),                           &
     & numb(mod((i+10)/10, 10)+1), i=1, cols, 10)
!
      write(10,'(6x,120a1)') (numb(mod(i, 10)+1), i=1,cols)

      do y=1,rows
        z=(y-1)*cols
        if(cols == 120) then
          write(10,123) y,(diff(x+z),x=1,cols)
 123      format(1x,i3,'->',120a1)
        else
          write(10,124) y,(diff(x+z),x=1,cols)
 124      format(1x,i3,'->',120a1/(6x,120a1))
        endif
      enddo

      RETURN
      END SUBROUTINE PRINT_DIF_MAP
