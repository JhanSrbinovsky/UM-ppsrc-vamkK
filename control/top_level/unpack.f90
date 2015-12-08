! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routines: UNPACK --------------------------------------------------
!
! Purpose:  This subroutine uses the index arrays
!           INDEX_COMP(N_SEG), INDEX_EXP(N_SEG) and
!           INDEX_TO_ROWS(MAX_ROW,MAX_LEVEL) to unpack data
!           FROM T_COMP(N_COMP) into T_EXP(IMAX,JMAX,KMAX)
!
! Programming standard : UMDP PAPER 3
!
! -----------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

      SUBROUTINE UNPACK(J1,J2,K1,K2,MAX_ROW,MAX_LEVEL,                  &
     &IMAX,JMAX,KMAX,INDEX_COMP,INDEX_EXP,N_SEG,INDEX_TO_ROWS,          &
     &N_COMP,T_COMP,T_EXP,REAL_MDI,CYCLIC)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!Argument variables
      INTEGER                                                           &
     & J1                                                               &
           ! (IN) FIRST ROW OF DATA TO BE UNPACKED
     &,J2                                                               &
           ! (IN) LAST ROW OF DATA TO BE UNPACKED
     &,K1                                                               &
           ! (IN) FIRST LEVEL OF DATA TO BE UNPACKED
     &,K2                                                               &
           ! (IN) LAST LEVEL OF DATA TO BE UNPACKED
     &,IMAX                                                             &
             ! (IN) NUMBER OF POINTS EAST-WEST
     &,JMAX                                                             &
             ! (IN) NUMBER OF POINTS NORTH-SOUTH
     &,KMAX                                                             &
             ! (IN) NUMBER OF POINTS IN VERTICAL
     &,N_SEG                                                            &
             ! (IN) TOTAL NUMBER OF SEA SEGMENTS
     &,N_COMP                                                           &
                ! (IN) DIMENSION OF COMPRESSED ARRAY
     &,MAX_ROW                                                          &
                ! (IN) UPPER LIMIT ON JMAX
     &,MAX_LEVEL                                                        &
                  ! (IN) UPPER LIMIT ON KMAX
     &,INDEX_COMP(N_SEG)                                                &
                          ! (IN) CONTAINS POSITIONS IN COMPRESSED ARRAY
!                         ! OF START OF EACH SEA SEGMENT
     &,INDEX_EXP(N_SEG)                                                 &
                         ! (IN) CONTAINS POSITIONS IN 1-DIMENSIONAL
!                        ! EXPANDED ARRAY
!                        ! T_EXP_1_D(IMAX*MAX_ROW*MAX_LEVEL)
!                        ! OF START OF EACH SEA SEGMENT
     &,INDEX_TO_ROWS(MAX_ROW,MAX_LEVEL)                                 &
                                         ! (IN) CONTAINS NUMBER OF FIRST
!                                        ! SEA SEGMENT IN EACH ROW
!                                        ! AT EACH LEVEL IF THERE IS A
!                                        ! SEA SEGMENT IN THE ROW
!                                        ! CONTAINS NUMBER OF NEXT
!                                        ! SEA SEGMENT OTHERWISE
     &,I_DATA  ! NUMBER OF DISTINCT DATA POINTS EAST-WEST

      REAL                                                              &
     & REAL_MDI                                                         &
                 ! (IN) REAL MISSING DATA INDICATOR
     &,T_COMP(N_COMP)                                                   &
                       ! (IN) COMPRESSED ARRAY
     &,T_EXP(IMAX,JMAX,KMAX)  ! (OUT) 3-DIMENSIONAL EXPANDED ARRAY

      LOGICAL                                                           &
     & CYCLIC  ! (IN) INDICATES WHETHER T_EXP(IMAX,JMAX,KMAX) INCLUDES
!              ! DATA FOR CYCLIC WRAP-AROUND POINTS

!Local variables
      INTEGER                                                           &
     & J60                                                              &
                    ! LOCAL LOOP INDEX FOR 60
     &,JSPAN                                                            &
                    ! LOCAL NUMBER OF ITERATIONS OF ROWS
     &,KSPAN                                                            &
                    ! LOCAL NUMBER OF ITERATIONS OF LEVELS
     &,I                                                                &
                    ! LOCAL LOOP INDEX FOR COLUMNS
     &,J                                                                &
                    ! LOCAL LOOP INDEX FOR ROWS
     &,K                                                                &
                    ! LOCAL LOOP INDEX FOR LEVELS
     &,JN,KN                                                            &
                    ! LOCAL VARIABLES
     &,NUM_SEG                                                          &
                    ! NUMBER OF SEA SEGMENTS IN PRESENT ROW
     &,LEN_SEG                                                          &
                    ! LENGTH OF PRESENT SEA SEGMENT
     &,COUNT                                                            &
                    ! LOCAL COUNTER FOR POINTS IN A SEA SEGMENT
     &,SEG                                                              &
                    ! LOCAL LOOP INDEX FOR SEA SEGMENTS IN PRESENT ROW
     &,SEG_POS                                                          &
                    ! LOCAL COUNTER FOR SEA SEGMENTS
     &,X_POS                                                            &
                    ! LOCAL COUNTER FOR POINTS IN A ROW
     &,IPOINT_EXP                                                       &
                    ! LOCAL POINTER TO EXPANDED ARRAY
     &,IPOINT_COMP  ! LOCAL POINTER TO COMPRESSED ARRAY

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------
!1. Fill the expanded array with real missing data indicators

      IF (lhook) CALL dr_hook('UNPACK',zhook_in,zhook_handle)
      DO K=1,KMAX
         DO J=1,JMAX
            DO I=1,IMAX
               T_EXP(I,J,K)=REAL_MDI
            END DO
         END DO
      END DO
!---------------------------------------------------------------------
!2. Set the wrap-around parameters

      IF (CYCLIC) THEN
         I_DATA=IMAX-2
      ELSE
         I_DATA=IMAX
      END IF
!---------------------------------------------------------------------
!3. Unpack levels and rows
!3.1 Set the number of iterations over rows and levels
      JSPAN = J2 - J1 + 1
      KSPAN = K2 - K1 + 1

      DO J60=0,JSPAN*KSPAN-1

!3.3 Set up 2-d loop indices

         K = J60/JSPAN
         J = J60 - K*JSPAN + J1
         K = K + K1

!3.4 Define the next row and level

         IF (J == MAX_ROW) THEN
            JN=1
            KN=K+1
         ELSE
           JN=J+1
           KN=K
         END IF

!3.5 Calculate the number of segments in the present row

         IF (KN >  MAX_LEVEL) THEN
            NUM_SEG=N_SEG-INDEX_TO_ROWS(J,K)+1
         ELSE
            NUM_SEG=INDEX_TO_ROWS(JN,KN)-INDEX_TO_ROWS(J,K)
         END IF

         DO SEG=1,NUM_SEG
            SEG_POS=INDEX_TO_ROWS(J,K)+SEG-1

!3.6 Calculate the length of the present sea segment

            IF (SEG_POS <  N_SEG) THEN
               LEN_SEG=INDEX_COMP(SEG_POS+1)-INDEX_COMP(SEG_POS)
            ELSE
               LEN_SEG=N_COMP-INDEX_COMP(SEG_POS)+1
            END IF

!3.7 Calculate t_exp(i,j,k) for each point in the segment

            DO COUNT=1,LEN_SEG
               IPOINT_EXP=INDEX_EXP(SEG_POS)+COUNT-1
               IPOINT_COMP=INDEX_COMP(SEG_POS)+COUNT-1
               X_POS=IPOINT_EXP-(K-1)*I_DATA*MAX_ROW-(J-1)*I_DATA
               T_EXP(X_POS,J-J1+1,K-K1+1)=T_COMP(IPOINT_COMP)
            END DO
         END DO
      END DO
!---------------------------------------------------------------------
!4. Calculate t_exp(i,j,k) for wrap-around points if necessary

      IF (CYCLIC) THEN
         DO K=1,KMAX
            DO J=1,JMAX
               T_EXP(IMAX-1,J,K)=T_EXP(1,J,K)
               T_EXP(IMAX,J,K)=T_EXP(2,J,K)
            END DO
         END DO
      END IF

      IF (lhook) CALL dr_hook('UNPACK',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UNPACK
