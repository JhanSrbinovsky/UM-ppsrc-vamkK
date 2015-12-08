! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Fills the LBC region of a LAM model with zeros

      SUBROUTINE ZERO_LATERAL_BOUNDARIES(                               &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,LEVELS,FLD_TYPE,FIELD,            &
     &  RIMWIDTH,AT_EXTREMITY,                                          &
     &  L_ZERO_BOUNDARIES,L_ZERO_HALOS)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParParams
      IMPLICIT NONE

!  Fills the LBC region of a LAM model with zeros. The region can
!  be split into two areas:
!  1) Boundary - true edge of the model data
!  2) Halo - extended "halo" edge of the model data
!  Both of these areas can be zeroed independently
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input
!
! Arguments:

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                          ! IN : number of points in row
     &, ROWS                                                            &
                          ! IN : number of rows
     &, HALO_I                                                          &
                          ! IN : size of halo in EW direction
     &, HALO_J                                                          &
                          ! IN : size of halo in NS direction
     &, LEVELS                                                          &
                          ! IN : number of vertical levels
     &, FLD_TYPE                                                        &
                          ! IN : type of input field (P,U or V)
     &, RIMWIDTH          ! IN : size (width) of the boundary area

      LOGICAL                                                           &
     &  AT_EXTREMITY(4)                                                 &
                          ! IN : Indicates if this processor is at
                          !      the edge (North,East,South,West)
                          !      of the processor grid
     &, L_ZERO_BOUNDARIES                                               &
                          ! IN : Do we want to zero the
                          !      boundary region
     &, L_ZERO_HALOS      ! IN : Do we want to zero the halo
                          !      region


      REAL                                                              &
     &  FIELD(1-HALO_I:ROW_LENGTH+HALO_I,                               &
                                           ! IN/OUT : field to
     &        1-HALO_J:ROWS+HALO_J,                                     &
                                           !          update
     &        LEVELS)

! Comdecks

! Local variables
      INTEGER                                                           &
     &  row_start_pt                                                    &
                          ! first point along row to update
     &, row_end_pt                                                      &
                          ! last point along row to update
     &, first_row                                                       &
                          ! first row to update
     &, last_row                                                        &
                          ! last row to update
     &, I,J,K             ! loop counters (along row, row)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Code

!--------------------------------------------------------------------
! Northern region

      IF (lhook) CALL dr_hook('ZERO_LATERAL_BOUNDARIES',zhook_in,zhook_handle)

      IF (AT_EXTREMITY(PNorth)) THEN  ! This processor has a northern
                                      ! LBC region
        IF (L_ZERO_HALOS) THEN        ! Zero out the halo area of the
                                      ! LBC region
          row_start_pt=1-HALO_I
          row_end_pt=ROW_LENGTH+HALO_I
          first_row=ROWS+1
          last_row=ROWS+HALO_J

          DO K=1,LEVELS
            DO J=first_row,last_row
              DO I=row_start_pt,row_end_pt
                FIELD(I,J,K)=0.0
              ENDDO
            ENDDO
          ENDDO

        ENDIF !  IF (L_ZERO_HALOS)

        IF (L_ZERO_BOUNDARIES) THEN   ! Zero out the boundary area
                                      ! of the LBC region
      If (AT_EXTREMITY(PWest)) Then
        row_start_pt = 1
      Else
        row_start_pt = 1-halo_i
      End If

      If (AT_EXTREMITY(PEast)) Then
        If (FLD_TYPE  ==  fld_type_u) Then
          row_end_pt=ROW_LENGTH-1
        Else
          row_end_pt=ROW_LENGTH
        End If
      Else
        If (FLD_TYPE  ==  fld_type_u) Then
          row_end_pt=ROW_LENGTH+halo_i-1
        Else
          row_end_pt=ROW_LENGTH+halo_i
        End If
      End If

          first_row=ROWS-RIMWIDTH+1
          last_row=ROWS

          DO K=1,LEVELS
            DO J=first_row,last_row
              DO I=row_start_pt,row_end_pt
                FIELD(I,J,K)=0.0
              ENDDO
            ENDDO
          ENDDO

        ENDIF ! IF (L_ZERO_BOUNDARIES)

      ENDIF ! IF (AT_EXTREMITY(PNorth))

!--------------------------------------------------------------------
! Eastern region

      IF (AT_EXTREMITY(PEast)) THEN   ! This processor has a eastern
                                      ! LBC region
        IF (L_ZERO_HALOS) THEN        ! Zero out the halo area of the
                                      ! LBC region

          IF (FLD_TYPE  ==  fld_type_u) THEN
            row_start_pt=ROW_LENGTH
          ELSE
            row_start_pt=ROW_LENGTH+1
          ENDIF
          row_end_pt=ROW_LENGTH+HALO_I

          IF (AT_EXTREMITY(PSouth)) THEN
            first_row=1
          ELSE
            first_row=1-HALO_J
          ENDIF
          IF (AT_EXTREMITY(PNorth)) THEN
            last_row=ROWS
          ELSE
            last_row=ROWS+HALO_J
          ENDIF

          DO K=1,LEVELS
            DO J=first_row,last_row
              DO I=row_start_pt,row_end_pt
                FIELD(I,J,K)=0.0
              ENDDO
            ENDDO
          ENDDO

        ENDIF !  IF (L_ZERO_HALOS)

        IF (L_ZERO_BOUNDARIES) THEN   ! Zero out the boundary area
                                      ! of the LBC region
          IF (FLD_TYPE  ==  fld_type_u) THEN
            row_start_pt=ROW_LENGTH-RIMWIDTH
            row_end_pt=ROW_LENGTH-1
          ELSE
            row_start_pt=ROW_LENGTH-RIMWIDTH+1
            row_end_pt=ROW_LENGTH
          ENDIF
      If (AT_EXTREMITY(PSouth)) Then
        first_row = 1
      Else
        first_row = 1-halo_j
      End If

      If (AT_EXTREMITY(PNorth)) Then
        last_row = rows
      Else
        last_row = rows+halo_j
      End IF


          DO K=1,LEVELS
            DO J=first_row,last_row
              DO I=row_start_pt,row_end_pt
                FIELD(I,J,K)=0.0
              ENDDO
            ENDDO
          ENDDO

        ENDIF ! IF (L_ZERO_BOUNDARIES)

      ENDIF ! IF (AT_EXTREMITY(PEast))

!--------------------------------------------------------------------
! Southern region

      IF (AT_EXTREMITY(PSouth)) THEN   ! This processor has a eastern
                                      ! LBC region
        IF (L_ZERO_HALOS) THEN        ! Zero out the halo area of the
                                      ! LBC region
          row_start_pt=1-HALO_I
          row_end_pt=ROW_LENGTH+HALO_I
          first_row=1-HALO_J
          last_row=0

          DO K=1,LEVELS
            DO J=first_row,last_row
              DO I=row_start_pt,row_end_pt
                FIELD(I,J,K)=0.0
              ENDDO
            ENDDO
          ENDDO

        ENDIF !  IF (L_ZERO_HALOS)

        IF (L_ZERO_BOUNDARIES) THEN   ! Zero out the boundary area
                                      ! of the LBC region
      If (AT_EXTREMITY(PWest)) Then
        row_start_pt = 1
      Else
        row_start_pt = 1-halo_i
      End If

      If (AT_EXTREMITY(PEast)) Then
        If (FLD_TYPE  ==  fld_type_u) Then
          row_end_pt=ROW_LENGTH-1
        Else
          row_end_pt=ROW_LENGTH
        End If
      Else
        If (FLD_TYPE  ==  fld_type_u) Then
          row_end_pt=ROW_LENGTH+halo_i-1
        Else
          row_end_pt=ROW_LENGTH+halo_i
        End If
      End If

          first_row=1
          last_row=RIMWIDTH

          DO K=1,LEVELS
            DO J=first_row,last_row
              DO I=row_start_pt,row_end_pt
                FIELD(I,J,K)=0.0
              ENDDO
            ENDDO
          ENDDO

        ENDIF ! IF (L_ZERO_BOUNDARIES)

      ENDIF ! IF (AT_EXTREMITY(PNorth))

!--------------------------------------------------------------------
! Western region

      IF (AT_EXTREMITY(PWest)) THEN   ! This processor has a western
                                      ! LBC region
        IF (L_ZERO_HALOS) THEN        ! Zero out the halo area of the
                                      ! LBC region
          row_start_pt=1-HALO_I
          row_end_pt=0

          IF (AT_EXTREMITY(PSouth)) THEN
            first_row=1
          ELSE
            first_row=1-HALO_J
          ENDIF
          IF (AT_EXTREMITY(PNorth)) THEN
            last_row=ROWS
          ELSE
            last_row=ROWS+HALO_J
          ENDIF

          DO K=1,LEVELS
            DO J=first_row,last_row
              DO I=row_start_pt,row_end_pt
                FIELD(I,J,K)=0.0
              ENDDO
            ENDDO
          ENDDO

        ENDIF !  IF (L_ZERO_HALOS)

        IF (L_ZERO_BOUNDARIES) THEN   ! Zero out the boundary area
                                      ! of the LBC region
          row_start_pt=1
          row_end_pt=RIMWIDTH
      If (AT_EXTREMITY(PSouth)) Then
        first_row = 1
      Else
        first_row = 1-halo_j
      End If

      If (AT_EXTREMITY(PNorth)) Then
        last_row = rows
      Else
        last_row = rows+halo_j
      End IF


          DO K=1,LEVELS
            DO J=first_row,last_row
              DO I=row_start_pt,row_end_pt
                FIELD(I,J,K)=0.0
              ENDDO
            ENDDO
          ENDDO

        ENDIF ! IF (L_ZERO_BOUNDARIES)

      ENDIF ! IF (AT_EXTREMITY(PEast))

      IF (lhook) CALL dr_hook('ZERO_LATERAL_BOUNDARIES',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ZERO_LATERAL_BOUNDARIES
