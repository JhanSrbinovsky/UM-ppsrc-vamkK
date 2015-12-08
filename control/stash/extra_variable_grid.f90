! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH
      SUBROUTINE EXTRA_VARIABLE_GRID(buf3,LEN_IN                        &
     &            ,SROW_OUT,N_ROWS_OUT,WCOL_OUT,N_COLS_OUT)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      IMPLICIT NONE
!
! Description: This routine controls the addition of
!              variable horizontal grid information in the
!              "extra data" part of a PP field.
!
!==========================================================
! Description: This include file contains information needed when
!              generating variable horizontal grid data in the
!              STASH extra data vector. Introduced UM 5.4 - R. Hill
!===================================================================
      LOGICAL :: X_VAR_GRID ! Whether variable grid in E-W direction
      LOGICAL :: Y_VAR_GRID ! and/or in S-N direction

      INTEGER :: VAR_GRID_TYPE ! 0 = none
                               ! 1 = T grid
                               ! 2 = U/V grid

      ! Grid boundaries for T and U,V
      REAL :: X_BOUNDARY(ROW_LENGTH_MAX+1,2)
      REAL :: Y_BOUNDARY(ROWS_MAX+1,2)

      ! Grid Points for T and U,V
      REAL :: X_GRID(ROW_LENGTH_MAX,2)
      REAL :: Y_GRID(ROWS_MAX,2)

      COMMON /OVARGRID/ X_VAR_GRID,Y_VAR_GRID                           &
     & ,X_BOUNDARY,Y_BOUNDARY,X_GRID,Y_GRID,VAR_GRID_TYPE

      ! The following parameters correspond to the extra data
      ! vector descriptors expected, for e.g., in PV-WAVE
      ! plotting routines (e.g. decode_extra.pro). There are
      ! numerous other areas of code where these integer
      ! descriptors must be handled (e.g. FIELDCOS, PPI2H, FTT)
      ! So it is not a trivial matter to introduce new code
      ! descriptors. Furthermore, ieee -32 will destroy these
      ! integers so PP data must always be processed via the
      ! long winded route: QXFIELDCOS -> PUTIBM -> FTT/PPI2H.
      ! (This thoroughly unsatisfactory state of affairs may
      ! be correctable with developments to ieee and convpp).
      INTEGER,PARAMETER :: x_coord_vector=1
                                     ! Indicates that an extra
                                     ! data vector gives LBNPT
                                     ! x-coordinate values
      INTEGER,PARAMETER :: y_coord_vector=2
                                     ! Indicates that an extra
                                     ! data vector gives LBROW
                                     ! Y-coordinate values
      INTEGER,PARAMETER :: x_lbnd_vector=12
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! x-boundary values
      INTEGER,PARAMETER :: x_ubnd_vector=13
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! x-boundary values
      INTEGER,PARAMETER :: y_lbnd_vector=14
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! y-boundary values
      INTEGER,PARAMETER :: y_ubnd_vector=15
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! y-boundary values

      INTEGER :: LEN_IN              ! IN size of buf3
      INTEGER :: SROW_OUT            ! IN 1st southern row
      INTEGER :: N_ROWS_OUT          ! IN No. of rows to output
      INTEGER :: WCOL_OUT            ! IN 1st western col
      INTEGER :: N_COLS_OUT          ! IN No. of cols to output

      REAL :: BUF3(LEN_IN)           ! IN/OUT data for output

      INTEGER :: EXTRA_START         ! Local start point for
                                     ! extra data
      INTEGER :: VECTOR_SIZE         ! Local size of vector

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      ! The start point for any extra data is at the first point
      ! after the incoming data.
      IF (lhook) CALL dr_hook('EXTRA_VARIABLE_GRID',zhook_in,zhook_handle)
      EXTRA_START = 1

      ! If grid is variable in the X direction, add appropriate
      ! extra data.
      IF (X_VAR_GRID) THEN

         ! Coordinates first
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_COLS_OUT             &
     &              ,X_GRID(WCOL_OUT,VAR_GRID_TYPE),VECTOR_SIZE         &
     &              ,x_coord_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE

         ! Lower boundaries
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_COLS_OUT             &
     &              ,X_BOUNDARY(WCOL_OUT,VAR_GRID_TYPE),VECTOR_SIZE     &
     &              ,x_lbnd_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE

         ! Upper boundaries
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_COLS_OUT             &
     &              ,X_BOUNDARY(WCOL_OUT+1,VAR_GRID_TYPE),VECTOR_SIZE   &
     &              ,x_ubnd_vector)

         EXTRA_START = EXTRA_START + VECTOR_SIZE
      ENDIF

      IF (Y_VAR_GRID) THEN
         ! Coordinates first
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_ROWS_OUT             &
     &              ,Y_GRID(SROW_OUT,VAR_GRID_TYPE),VECTOR_SIZE         &
     &              ,y_coord_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE

         ! Lower boundaries
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_ROWS_OUT             &
     &              ,Y_BOUNDARY(SROW_OUT,VAR_GRID_TYPE),VECTOR_SIZE     &
     &              ,y_lbnd_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE

         ! Upper boundaries
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_ROWS_OUT             &
     &              ,Y_BOUNDARY(SROW_OUT+1,VAR_GRID_TYPE),VECTOR_SIZE   &
     &              ,y_ubnd_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE
      ENDIF

      IF (lhook) CALL dr_hook('EXTRA_VARIABLE_GRID',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE EXTRA_VARIABLE_GRID

