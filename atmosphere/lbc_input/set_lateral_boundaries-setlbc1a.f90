! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Fills the LBC region of a LAM model with data from LBCs

      SUBROUTINE SET_LATERAL_BOUNDARIES(                                &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,LEVELS,FLD_TYPE,FIELD,            &
     &  LENRIM,LBC_SIZE,LBC_START,LBC_HALO_I,LBC_HALO_J,LBC,            &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_DO_BOUNDARIES,L_DO_HALOS)

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParParams
      IMPLICIT NONE

! Fills the LBC region of a LAM model FIELD with data from
! the LBC record.
! The N_RIMS_TO_DO variable allows only a certain depth of
! rim width to be done.
! The RIMWEIGHTS array contains a weighting factor that
! controls what proportion of FIELD is used and what
! proportion of LBC. In the halo area only LBC data is used.
!
! If RIMWIDTH is zero then the routine assumes there are no
! LBCs available and exits without updating any fields
!
! The logicals L_DO_BOUNDARIES and L_DO_HALOS allow each
! of the two components of the LBC region to be updated
! independently.
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
                          ! IN : size of FIELD halo in EW direction
     &, HALO_J                                                          &
                          ! IN : size of FIELD halo in NS direction
     &, LEVELS                                                          &
                          ! IN : number of vertical levels
     &, FLD_TYPE                                                        &
                          ! IN : type of input field (P,U or V)
     &, LENRIM                                                          &
                          ! IN : size of one level of LBC data
     &, LBC_SIZE(4)                                                     &
                          ! IN : size of each side (NESW) of LBC data
     &, LBC_START(4)                                                    &
                          ! IN : offset of each side of LBC data
     &, LBC_HALO_I                                                      &
                          ! IN : size of LBC halo in EW direction
     &, LBC_HALO_J                                                      &
                          ! IN : size of LBC halo in NS direction
     &, RIMWIDTH                                                        &
                          ! IN : size (width) of the boundary area
     &, N_RIMS_TO_DO      ! IN : number of rims to do (counting from
                          !      the outside in)

      LOGICAL                                                           &
     &  AT_EXTREMITY(4)                                                 &
                          ! IN : Indicates if this processor is at
                          !      the edge (North,East,South,West)
                          !      of the processor grid
     &, L_DO_BOUNDARIES                                                 &
                          ! IN : .TRUE. if the boundary region is
                          !      to be updated
     &, L_DO_HALOS        ! IN : .TRUE. if the halo region is to
                          !      updated

      REAL                                                              &
     &  RIMWEIGHTS(RIMWIDTH)                                            &
                          ! IN : Weight to apply to each successive
                          !      rim
     &, LBC(LENRIM,LEVELS)                                              &
                          ! IN : LBC to apply to FIELD
     &, FIELD(1-HALO_I:ROW_LENGTH+HALO_I,                               &
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
     &, lbc_row_len                                                     &
                          ! Length of a row of LBC data
     &, first_pt_of_lbc                                                 &
                          ! First point on row contained in LBC data
     &, first_row_of_lbc                                                &
                          ! First row contained in LBC data
     &, lbc_address                                                     &
                          ! address in the LBC array
     &, weights_I                                                       &
                          ! modified I to point to ew_weights array
     &, weights_J                                                       &
                          ! modified J to point to ns_weights array
     &, rim                                                             &
                          ! rim number (at corners)
     &, I,J,K                                                           &
                          ! loop counters (along row, row, level)
     &, rim_I             ! modified I to point to RIMWEIGHTS array

      REAL                                                              &
     &  ns_weights(1-HALO_I:ROW_LENGTH+HALO_I,                          &
     &             1-HALO_J:RIMWIDTH)                                   &
     &, ew_weights(1-HALO_I:RIMWIDTH,                                   &
     &             1-HALO_J:ROWS+HALO_J)
                          ! Arrays which contain the weights
                          ! to apply to the LBCs at each point

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!---------------------------------------------------------------
!
! The following diagram breaks up a LAM model area into a number
! of subcomponents (the letters will be referred to in the code):
! (Assumes the following model sizes for this example:
!  ROW_LENGTH=ROWS=12
!  RIMWIDTH=3
!  HALO_I=HALO_J=2)
!
!       North
!  aaaaaaaaaaaaaaaa
!  aaaaaaaaaaaaaaaa
!  bbcccddddddeeeff
!  bbcccddddddeeeff
!  bbcccddddddeeeff
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  llmmmnnnnnnooopp
!  llmmmnnnnnnooopp
!  llmmmnnnnnnooopp
!  qqqqqqqqqqqqqqqq
!  qqqqqqqqqqqqqqqq
!       South

! Code

      IF (lhook) CALL dr_hook('SET_LATERAL_BOUNDARIES',zhook_in,zhook_handle)

!=================================================================
! 1.0 Check RIMWIDTH. If RIMWIDTH is zero then there are no LBCs
!     and we can return without doing anything

      IF (RIMWIDTH  ==  0) GOTO 9999

!=================================================================
! 2. Set up the "*_weights" arrays

!-----------------------------------------------------------------
! 2.1 Do the North/South regions

      IF (AT_EXTREMITY(PNorth) .OR. AT_EXTREMITY(PSouth)) THEN
! We will need the ns_weights array

! First we do the halo area (South : qlp, North : abf)

! The NS edge (South: q, North : a)
        DO J=1-HALO_J,0
          DO I=1-HALO_I,ROW_LENGTH+HALO_I
            IF (L_DO_HALOS) THEN
              ns_weights(I,J)=1.0
            ELSE ! don't update halos
              ns_weights(I,J)=0.0
            ENDIF ! IF (L_DO_HALOS
          ENDDO ! I
        ENDDO ! J

! The Western edge ( South: l, North : b)
        IF (AT_EXTREMITY(PWest)) THEN
        DO J=1,RIMWIDTH
          DO I=1-HALO_I,0
              IF (L_DO_HALOS) THEN
                ns_weights(I,J)=1.0
              ELSE ! don't update halos
                ns_weights(I,J)=0.0
              ENDIF ! IF (L_DO_HALOS)
          ENDDO ! I
        ENDDO ! J
        ENDIF ! IF (AT_EXTREMITY(PWest))

! The Eastern edge (South : p, North : f)
        IF (AT_EXTREMITY(PEast)) THEN
        DO J=1,RIMWIDTH
          DO I=ROW_LENGTH+1,ROW_LENGTH+HALO_I
              IF (L_DO_HALOS) THEN
                ns_weights(I,J)=1.0
              ELSE ! don't update halos
                ns_weights(I,J)=0.0
              ENDIF ! IF (L_DO_HALOS)
          ENDDO ! I
        ENDDO ! J
        ENDIF ! IF (AT_EXTREMITY(PEast))

! And now the boundary areas ( South : mno, North : cde )
! The Western corner (South : m, North : c)
        IF (AT_EXTREMITY(PWest)) THEN
        DO J=1,RIMWIDTH
          DO I=1,RIMWIDTH
            IF (L_DO_BOUNDARIES) THEN
                rim=MIN(I,J)
                IF (rim  <=  N_RIMS_TO_DO) THEN
                  ns_weights(I,J)=RIMWEIGHTS(rim)
                ELSE ! This rim not required
                  ns_weights(I,J)=0.0
                ENDIF !  IF (rim  <=  N_RIMS_TO_DO)
            ELSE ! don't update boundaries
              ns_weights(I,J)=0.0
            ENDIF ! IF (L_DO_BOUNDARIES)
          ENDDO ! I
        ENDDO ! J
        ENDIF ! IF (AT_EXTREMITY(PWest))

! The Eastern corner (South : o, North : e)
        IF (AT_EXTREMITY(PEast)) THEN
        DO J=1,RIMWIDTH
          DO I=ROW_LENGTH-RIMWIDTH+1,ROW_LENGTH
            rim_I=ROW_LENGTH+1-I
            IF (L_DO_BOUNDARIES) THEN
                rim=MIN(rim_I,J)
                IF (rim  <=  N_RIMS_TO_DO) THEN
                  ns_weights(I,J)=RIMWEIGHTS(rim)
                ELSE ! This rim not required
                  ns_weights(I,J)=0.0
                ENDIF !  IF (rim  <=  N_RIMS_TO_DO)
            ELSE ! don't update boundaries
              ns_weights(I,J)=0.0
            ENDIF ! IF (L_DO_BOUNDARIES)
          ENDDO ! I
        ENDDO ! J
        ENDIF ! IF (AT_EXTREMITY(PEast))

! The bit between the corners (South : n, North : d)
        IF (AT_EXTREMITY(PEast)) THEN
          row_end_pt=ROW_LENGTH-RIMWIDTH
        else
          row_end_pt=row_length +halo_i
        endif
        IF (AT_EXTREMITY(Pwest)) THEN
          row_start_pt=RIMWIDTH+1
        else
          row_start_pt=1-halo_i
        endif
        DO J=1,RIMWIDTH
          DO I=row_start_pt,row_end_pt
            IF (L_DO_BOUNDARIES) THEN
              IF (J  <=  N_RIMS_TO_DO) THEN
                ns_weights(I,J)=RIMWEIGHTS(J)
              ELSE ! This rim is not required
                ns_weights(I,J)=0.0
              ENDIF ! IF (J  <=  N_RIMS_TO_DO)
           ELSE ! don't update boundaries
              ns_weights(I,J)=0.0
            ENDIF ! IF (L_DO_BOUNDARIES)
          ENDDO ! I
        ENDDO ! J

! Finally, for New dynamics we need to check if this is U-grid data and 
! if so, if this processor is at the far right (ie. East) of the grid. In
! which case we need to shift all the end of row data left by
! one because the LAM U fields have one less point of data on
! each row. This does not apply to Endgame.
        IF (.NOT. l_vatpoles) THEN
        IF ((FLD_TYPE  ==  fld_type_u) .AND.                            &
     &      (AT_EXTREMITY(PEast))) THEN

          DO J=1-HALO_J,RIMWIDTH
            DO I=ROW_LENGTH-RIMWIDTH,ROW_LENGTH+HALO_I-1
              ns_weights(I,J)=ns_weights(I+1,J)
            ENDDO ! I
          ENDDO ! J

        ENDIF ! If U-grid and at east edge of LAM
        END IF  ! vatpoles

      ENDIF ! If we're at the North or South edge of the LAM

!-----------------------------------------------------------------
! 2.2 Do the East/West regions

      IF (AT_EXTREMITY(PWest) .OR. AT_EXTREMITY(PEast)) THEN
! We will need the ew_weights array

        IF (AT_EXTREMITY(PSouth)) THEN
          first_row=RIMWIDTH+1
        ELSE ! not at the South
          first_row=1-HALO_J
        ENDIF

        IF (AT_EXTREMITY(PNorth)) THEN
          last_row=ROWS-RIMWIDTH
        ELSE ! not at the North
          last_row=ROWS+HALO_J
        ENDIF

! First the halo area (West : g , East k)

        DO J=first_row,last_row
          DO I=1-HALO_I,0
            IF (L_DO_HALOS) THEN
              ew_weights(I,J)=1.0
            ELSE ! don't update halos
              ew_weights(I,J)=0.0
            ENDIF ! IF (L_DO_HALOS)
          ENDDO ! I
        ENDDO ! J

! And now the boundary area (West : h, East : j)

        DO J=first_row,last_row
          DO I=1,RIMWIDTH
            IF (L_DO_BOUNDARIES) THEN
              IF (I  <=  N_RIMS_TO_DO) THEN
                ew_weights(I,J)=RIMWEIGHTS(I)
              ELSE ! this rim not required
                ew_weights(I,J)=0.0
              ENDIF ! IF (I  <=  N_RIMS_TO_DO)
            ELSE ! don't update boundaries
              ew_weights(I,J)=0.0
            ENDIF ! IF (L_DO_BOUNDARIES)
          ENDDO ! I
        ENDDO ! J

      ENDIF ! IF at Western or Eastern edge

!=================================================================
! 3 Now apply the boundary conditions

!-----------------------------------------------------------------
! 3.1 Southern region

      IF (AT_EXTREMITY(PSouth)) THEN

        IF (L_DO_HALOS) THEN
          first_row=1-HALO_J
          row_start_pt=1-HALO_I
          row_end_pt=ROW_LENGTH+HALO_I
        ELSE
          first_row=1
        If (at_extremity(PWest)) Then
          row_start_pt = 1
        Else
          row_start_pt=1-halo_i
        End If
        If (at_extremity(PEast)) Then
          row_end_pt = row_length
        Else
          row_end_pt = row_length+halo_i
        End If
        ENDIF

        IF (L_DO_BOUNDARIES) THEN
          last_row=RIMWIDTH
        ELSE
          last_row=0
        ENDIF

        IF ( l_vatpoles ) THEN
        lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I
        ELSE
        IF ((FLD_TYPE  ==  fld_type_u) .AND.                            &
     &    (AT_EXTREMITY(PEast))) THEN
          lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I - 1
          row_end_pt=row_end_pt-1
        ELSE
          lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I
        ENDIF
        END IF  ! vatpoles

        DO K=1,LEVELS
          DO J=first_row,last_row
            DO I=row_start_pt,row_end_pt
              IF (ns_weights(I,J)  /=  0.0) THEN

                LBC_address=LBC_START(PSouth) +                         &
     &                      (J+LBC_HALO_J-1)*lbc_row_len +              &
     &                      I+LBC_HALO_I-1

                IF (ns_weights(I,J)  ==  1.0) THEN
                  FIELD(I,J,K)= LBC(LBC_address,K)
                ELSE
                  FIELD(I,J,K)=FIELD(I,J,K)*(1.0-ns_weights(I,J)) +     &
     &                         LBC(LBC_address,K)*ns_weights(I,J)
                ENDIF
              ENDIF
            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

        IF((.not. L_DO_BOUNDARIES) .and. L_DO_HALOS) Then
          first_row = 1
          last_row = rimwidth

!      boundary area l
          If (AT_EXTREMITY(PWest)) Then
            row_start_pt = 1-halo_i
            row_end_pt = 0
            DO K=1,LEVELS
              DO J=first_row,last_row
                DO I=row_start_pt,row_end_pt
                  IF (ns_weights(I,J)  /=  0.0) THEN
                    LBC_address=LBC_START(PSouth) +                     &
     &                      (J+LBC_HALO_J-1)*lbc_row_len +              &
     &                      I+LBC_HALO_I-1
                    FIELD(I,J,K)=LBC(LBC_address,K)
                  END IF
                END DO
              END DO
            END DO
          ENDIF   !  (AT_EXTREMITY(PWest))

!      boundary area p
          IF (AT_EXTREMITY(PEast)) then
            row_start_pt = row_length+1
            row_end_pt = row_length+halo_i
            IF (.NOT. l_vatpoles) THEN
            IF (fld_type  ==  fld_type_u) THEN
              row_start_pt = row_start_pt-1
              row_end_pt = row_end_pt-1
            END IF
            END IF  ! vatpoles
            DO K=1,LEVELS
              DO J=first_row,last_row
                DO I=row_start_pt,row_end_pt
                  IF (ns_weights(I,J)  /=  0.0) THEN
                    LBC_address=LBC_START(PSouth) +                     &
     &                      (J+LBC_HALO_J-1)*lbc_row_len +              &
     &                      I+LBC_HALO_I-1
                    FIELD(I,J,K)=LBC(LBC_address,K)
                  END IF
                END DO
              END DO
            END DO
          END IF    ! (AT_EXTREMITY(PEast))

        END IF      ! ((.not. L_DO_BOUNDARIES) .and. L_DO_HALOS)

      ENDIF ! IF (AT_EXTREMITY(PSouth))

!-----------------------------------------------------------------
! 3.2 Northern region

      IF (AT_EXTREMITY(PNorth)) THEN

        IF (L_DO_HALOS) THEN
          last_row=ROWS+HALO_J
          row_start_pt=1-HALO_I
          row_end_pt=ROW_LENGTH+HALO_I
        ELSE
          last_row=ROWS
        If (at_extremity(PWest)) Then
          row_start_pt = 1
        Else
          row_start_pt=1-halo_i
        End If
        If (at_extremity(PEast)) Then
          row_end_pt = row_length
        Else
          row_end_pt = row_length+halo_i
        End If

        ENDIF

        IF (L_DO_BOUNDARIES) THEN
          first_row=ROWS-RIMWIDTH+1
        ELSE
          first_row=ROWS+1
        ENDIF

        IF ( l_vatpoles ) THEN
        lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I
        ELSE
        IF ((FLD_TYPE  ==  fld_type_u) .AND.                            &
     &    (AT_EXTREMITY(PEast))) THEN
          lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I - 1
          row_end_pt=row_end_pt-1
        ELSE
          lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I
        ENDIF
        END IF  ! vatpoles

        DO K=1,LEVELS
          DO J=first_row,last_row
            weights_J=ROWS+1-J
            DO I=row_start_pt,row_end_pt
              IF (ns_weights(I,weights_J)  /=  0.0) THEN

                LBC_address=LBC_START(PNorth) +                         &
     &                      (J-(ROWS-RIMWIDTH)-1)*lbc_row_len +         &
     &                      I+LBC_HALO_I-1

                IF (ns_weights(I,weights_J)  ==  1.0) THEN
                  FIELD(I,J,K)=LBC(LBC_address,K)
                ELSE
                  FIELD(I,J,K)=                                         &
     &              FIELD(I,J,K)*(1.0-ns_weights(I,weights_J)) +        &
     &              LBC(LBC_address,K)*ns_weights(I,weights_J)
                ENDIF
              ENDIF
            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

        IF((.not. L_DO_BOUNDARIES) .and. L_DO_HALOS) Then
          first_row = rows-rimwidth+1
          last_row = rows

!      boundary area b
          If (AT_EXTREMITY(PWest)) Then
            row_start_pt = 1-halo_i
            row_end_pt = 0
            DO K=1,LEVELS
              DO J=first_row,last_row
                weights_J=ROWS+1-J
                DO I=row_start_pt,row_end_pt
                  IF (ns_weights(I,weights_J)  /=  0.0) THEN
                    LBC_address=LBC_START(PNorth) +                     &
     &                      (J-(ROWS-RIMWIDTH)-1)*lbc_row_len +         &
     &                      I+LBC_HALO_I-1
                    FIELD(I,J,K)=LBC(LBC_address,K)
                  END IF
                END DO
              END DO
            END DO
          ENDIF   !  (AT_EXTREMITY(PWest))

!      boundary area f
          IF (AT_EXTREMITY(PEast)) then
            row_start_pt = row_length+1
            row_end_pt = row_length+halo_i
            IF (.NOT. l_vatpoles) THEN
            If (FLD_TYPE  ==  fld_type_u) then
              row_start_pt = row_start_pt-1
              row_end_pt = row_end_pt-1
            End IF
            END IF  ! vatpoles
            DO K=1,LEVELS
              DO J=first_row,last_row
                weights_J=ROWS+1-J
                DO I=row_start_pt,row_end_pt
                  IF (ns_weights(I,weights_J)  /=  0.0) THEN
                    LBC_address=LBC_START(PNorth) +                     &
     &                      (J-(ROWS-RIMWIDTH)-1)*lbc_row_len +         &
     &                      I+LBC_HALO_I-1
                   FIELD(I,J,K)=LBC(LBC_address,K)
                  END IF
                END DO
              END DO
            END DO
          END IF    ! (AT_EXTREMITY(PEast))

        END IF      ! ((.not. L_DO_BOUNDARIES) .and. L_DO_HALOS)

      ENDIF !  IF (AT_EXTREMITY(PNorth))

!-----------------------------------------------------------------
! 3.3 Western region

      IF (AT_EXTREMITY(PWest)) THEN

        IF (L_DO_HALOS) THEN
          row_start_pt=1-HALO_I
        ELSE
          row_start_pt=1
        ENDIF

        IF (L_DO_BOUNDARIES) THEN
          row_end_pt=N_RIMS_TO_DO
        ELSE
          row_end_pt=0
        ENDIF

        IF (AT_EXTREMITY(PSouth)) THEN
          first_row=RIMWIDTH+1
          first_row_of_LBC=RIMWIDTH+1
        ELSE ! Not at the South
          first_row_of_LBC=1-LBC_HALO_J
            first_row=1-HALO_J
        ENDIF ! IF (AT_EXTREMITY(PSouth))

        IF (AT_EXTREMITY(PNorth)) THEN
          last_row=ROWS-RIMWIDTH
        ELSE ! Not at the North
            last_row=ROWS+HALO_J
        ENDIF ! IF (AT_EXTREMITY(PNorth))

        lbc_row_len=LBC_HALO_I+RIMWIDTH

        DO K=1,LEVELS
          DO J=first_row,last_row
            DO I=row_start_pt,row_end_pt
              IF (ew_weights(I,J)  /=  0.0) THEN

                LBC_address=LBC_START(PWest)+                           &
     &                      (J-first_row_of_LBC)*lbc_row_len +          &
     &                      I+LBC_HALO_I-1

                IF (ew_weights(I,J)  ==  1.0) THEN
                  FIELD(I,J,K)=LBC(LBC_address,K)
                ELSE
                  FIELD(I,J,K)=FIELD(I,J,K)*(1.0 - ew_weights(I,J)) +   &
     &                         LBC(LBC_address,K)*ew_weights(I,J)
                ENDIF
              ENDIF
            ENDDO ! I
          ENDDO ! J
        ENDDO ! 0K

      ENDIF ! IF (AT_EXTREMITY(PWest))

!-----------------------------------------------------------------
! 3.3 Eastern region

      IF (AT_EXTREMITY(PEast)) THEN

        IF (L_DO_HALOS) THEN
          row_end_pt=ROW_LENGTH+HALO_I
        ELSE
          row_end_pt=ROW_LENGTH
        ENDIF

        IF (L_DO_BOUNDARIES) THEN
          row_start_pt=ROW_LENGTH-N_RIMS_TO_DO+1
        ELSE
          row_start_pt=ROW_LENGTH+1
        ENDIF

        IF (AT_EXTREMITY(PSouth)) THEN
          first_row=RIMWIDTH+1
          first_row_of_LBC=RIMWIDTH+1
        ELSE ! Not at the South
          first_row_of_LBC=1-LBC_HALO_J
            first_row=1-HALO_J
        ENDIF ! IF (AT_EXTREMITY(PSouth))

        IF (AT_EXTREMITY(PNorth)) THEN
          last_row=ROWS-RIMWIDTH
        ELSE ! Not at the North
            last_row=ROWS+HALO_J
        ENDIF ! IF (AT_EXTREMITY(PNorth))

        lbc_row_len=LBC_HALO_I+RIMWIDTH

        IF ( l_vatpoles ) THEN
        first_pt_of_LBC=ROW_LENGTH-RIMWIDTH+1
        ELSE
        IF (FLD_TYPE  ==  fld_type_u) THEN
          row_start_pt=row_start_pt-1
          row_end_pt=row_end_pt-1
          first_pt_of_LBC=ROW_LENGTH-RIMWIDTH
        ELSE
          first_pt_of_LBC=ROW_LENGTH-RIMWIDTH+1
        END IF
        END IF  ! vatpoles

        DO K=1,LEVELS
          DO J=first_row,last_row
            DO I=row_start_pt,row_end_pt
              weights_I=first_pt_of_LBC+RIMWIDTH-I
              IF (ew_weights(weights_I,J)  /=  0.0) THEN

                LBC_address=LBC_START(PEast)+                           &
     &                      (J-first_row_of_LBC)*lbc_row_len +          &
     &                      I-first_pt_of_LBC

                IF (ew_weights(weights_I,J)  ==  1.0) THEN
                  FIELD(I,J,K)=LBC(LBC_address,K)
                ELSE
                  FIELD(I,J,K)=                                         &
     &              FIELD(I,J,K)*(1.0 - ew_weights(weights_I,J)) +      &
     &              LBC(LBC_address,K)*ew_weights(weights_I,J)
                ENDIF
              ENDIF
            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

      ENDIF ! IF (AT_EXTREMITY(PEast))

!=================================================================
! 4 Tidy up
!  Endgame: 
!   If we're at the Eastern edge of a V-grid or P-grid field, the last point
!   on each row has not been updated. V-grid fields will be filled with zeros,
!   and P-grid fields with the value of their western neighbour.
!  New Dynamics:
!   If we're at the Eastern edge of a U grid field, the last point
!   on each row has not been updated, so we fill it with zeros here

      IF ( l_vatpoles ) THEN
! Zero last v for Endgame
      IF ( AT_EXTREMITY(PEast) .AND.                                    &
           (fld_type  ==  fld_type_v)) THEN
        DO k = 1, levels
          DO j = 1-halo_j,rows+halo_j
            FIELD(row_length+halo_i,J,K)=0.0
          END DO ! J
        END DO ! K
      END IF ! If at East and a V field (Endgame)

! Copy adjacent value for P field in Endgame
      IF ( AT_EXTREMITY(PEast) .AND.                                    &
           (fld_type  ==  fld_type_p)) THEN
        DO k = 1, levels
          DO j = 1-halo_j,rows+halo_j
            FIELD(row_length+halo_i,J,K)=FIELD(row_length+halo_i-1,j,K)
          END DO ! J
        END DO ! K
      END IF ! If at East and a P field (Endgame)


      ELSE
! Zero last U for New Dynamics
      IF ( AT_EXTREMITY(PEast) .AND.                                    &
           (fld_type  ==  fld_type_u)) THEN
        DO k = 1, levels
          DO j = 1-halo_j,rows+halo_j
            FIELD(row_length+halo_i,J,K)=0.0
          END DO ! J
        END DO ! K
      END IF ! If at East and a U field (ND)
      END IF  ! vatpoles

 9999 CONTINUE

      IF (lhook) CALL dr_hook('SET_LATERAL_BOUNDARIES',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SET_LATERAL_BOUNDARIES
