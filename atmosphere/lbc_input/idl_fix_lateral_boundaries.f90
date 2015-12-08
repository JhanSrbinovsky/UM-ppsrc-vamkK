! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Sets the LBCs with data from the LBC region of a LAM model

      SUBROUTINE IDL_FIX_LATERAL_BOUNDARIES(                            &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,LEVELS,FLD_TYPE,FIELD,            &
     &  LENRIM,LBC_SIZE,LBC_START,LBC_HALO_I,LBC_HALO_J,LBC,            &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_DO_BOUNDARIES,L_DO_HALOS)

! Purpose:
! Fills the LBC record with data from the LBC region of
! a LAM model FIELD. The FIELD has halo_i, halo_j and may be
! holding a field of different type for work purposes
!
! N_RIMS_TO_DO should be set to RIMWIDTH
!
! If RIMWIDTH is zero then the routine assumes there are no
! LBCs available and exits without updating any fields
!
! The logicals L_DO_BOUNDARIES and L_DO_HALOS should both be .true.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParParams
      IMPLICIT NONE

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
                                           ! INPUT: field
     &        1-HALO_J:ROWS+HALO_J,                                     &
                                           !  to put in lbc
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
     &, I,J,K             ! loop counters (along row, row, level)

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
!       South
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
!       North

! Code

      IF (lhook) CALL dr_hook('IDL_FIX_LATERAL_BOUNDARIES',zhook_in,zhook_handle)

!=================================================================
! 1.0 Check RIMWIDTH. If RIMWIDTH is zero then there are no LBCs
!     and we can return without doing anything

      IF (RIMWIDTH  ==  0) THEN 
        IF (lhook) CALL dr_hook('IDL_FIX_LATERAL_BOUNDARIES',zhook_out,zhook_handle)
        RETURN
      END IF

!=================================================================
! 2 Apply the boundary conditions

!-----------------------------------------------------------------
! 2.1 Southern region
!-----------------------------------------------------------------

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

        IF ((FLD_TYPE  ==  fld_type_u) .AND.                            &
     &    (AT_EXTREMITY(PEast))) THEN
          lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I - 1
          row_end_pt=row_end_pt-1
        ELSE
          lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I
        ENDIF

        DO K=1,LEVELS
          DO J=first_row,last_row
            DO I=row_start_pt,row_end_pt

              LBC_address=LBC_START(PSouth) +                           &
     &                      (J+LBC_HALO_J-1)*lbc_row_len +              &
     &                      I+LBC_HALO_I-1

              LBC(LBC_address,K) = FIELD(I,J,K)

            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

      ENDIF ! IF (AT_EXTREMITY(PSouth))

!-----------------------------------------------------------------
! 2.2 Northern region
!-----------------------------------------------------------------

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

        IF ((FLD_TYPE  ==  fld_type_u) .AND.                            &
     &    (AT_EXTREMITY(PEast))) THEN
          lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I - 1
          row_end_pt=row_end_pt-1
        ELSE
          lbc_row_len=ROW_LENGTH + 2*LBC_HALO_I
        ENDIF

        DO K=1,LEVELS
          DO J=first_row,last_row
            DO I=row_start_pt,row_end_pt

              LBC_address=LBC_START(PNorth) +                           &
     &                      (J-(ROWS-RIMWIDTH)-1)*lbc_row_len +         &
     &                      I+LBC_HALO_I-1

              LBC(LBC_address,K) = FIELD(I,J,K)

            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

      ENDIF !  IF (AT_EXTREMITY(PNorth))

!-----------------------------------------------------------------
! 2.3 Western region
!-----------------------------------------------------------------

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

              LBC_address=LBC_START(PWest)+                             &
     &                      (J-first_row_of_LBC)*lbc_row_len +          &
     &                      I+LBC_HALO_I-1

              LBC(LBC_address,K) = FIELD(I,J,K)

            ENDDO ! I
          ENDDO ! J
        ENDDO ! 0K

      ENDIF ! IF (AT_EXTREMITY(PWest))

!-----------------------------------------------------------------
! 2.4 Eastern region
!-----------------------------------------------------------------

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

        IF (FLD_TYPE  ==  fld_type_u) THEN
          row_start_pt=row_start_pt-1
          row_end_pt=row_end_pt-1
          first_pt_of_LBC=ROW_LENGTH-RIMWIDTH
        ELSE
          first_pt_of_LBC=ROW_LENGTH-RIMWIDTH+1
        ENDIF

        DO K=1,LEVELS
          DO J=first_row,last_row
            DO I=row_start_pt,row_end_pt

              LBC_address=LBC_START(PEast)+                             &
     &                      (J-first_row_of_LBC)*lbc_row_len +          &
     &                      I-first_pt_of_LBC

              LBC(LBC_address,K) = FIELD(I,J,K)

            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

      ENDIF ! IF (AT_EXTREMITY(PEast))

      IF (lhook) CALL dr_hook('IDL_FIX_LATERAL_BOUNDARIES',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_FIX_LATERAL_BOUNDARIES
