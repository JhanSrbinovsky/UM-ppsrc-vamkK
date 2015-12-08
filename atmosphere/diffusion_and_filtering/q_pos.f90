! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!----------------------------------------------------------------------
! Subroutine Q_POS - Does the work.
!----------------------------------------------------------------------
      subroutine q_pos                                                  &
     &    (q, row_length, rows, q_levels,                               &
     &     l_q_pos_local, qlimit_in, model_domain)
!
!   PURPOSE:   REMOVES VALUES OF q BELOW qlimit_in.  THERE ARE TWO 
!              ALTERNATIVE METHODS:- 
!              METHOD 2 IS USED IF L_q_POS_LOCAL = .TRUE.
!              METHOD 1 IS THE ORIGINAL SCHEME AND IS USED IF METHOD
!              2 FAILS OR IF L_q_POS_LOCAL = .FALSE.
!              METHOD 2 IS DESIGNED TO ELIMINATE THE
!              VERY SLOW CLIMATE DRIFT IN q IN UPPER MODEL LEVELS.
!              IT IS UNNECESSARILY COMPLICATED (AND EXPENSIVE) FOR
!              FORECAST USE.
!
!        METHOD 1:  RESCALES q by SUBTRACTING q_limit_in. ACCUMULATES 
!              TOTALS OF NEGATIVE AND POSITIVE VALUES ON A LEVEL, 
!              ZEROING ALL NEGATIVE POINTS AND PROPORTIONALLY REMOVING 
!              THE SUM OF THE NEGATIVES FROM ALL POSITIVE POINTS.
!              THEN RESCALES q BY ADDING q_limit_in.
!
!        METHOD 2:  THE FOLLOWING METHOD IS APPLIED AT EACH LAYER:
!            RESCALE q by SUBTRACTING q_limit_in.
!            STEP 1 - LOOP ROUND ALL POINTS. IF q IS NEGATIVE,
!              (a):
!              SUM q ROUND NEIGHBOURING POINTS WITHIN ONE ROW AND
!              ONE COLUMN OF POINT IN CURRENT LAYER.  IF THIS VALUE
!              IS POSITIVE AND LARGE ENOUGH,THE CENTRE NEGATIVE VALUE
!              IS SET EQUAL TO ZERO AND THE NEIGHBOURS ARE SCALED
!              TO CONSERVE q. IF (a) DOES NOT WORK
!              THEN (b):
!              EXTEND SUMMATION WITH A LARGER SEARCH RADIUS AND REPEAT
!              PROCESS. (SEARCH RADIUS = 4 FOUND TO BE SUFFICIENT)
!            STEP 2 - IF GLOBAL MODEL, CHECK FOR NEGATIVE q AT POLES
!              AND CORRECT USING ALL VALUES IN NEIGHBOURING ROW.
!              RESCALE q BY ADDING q_limit_in.
!            STEP 3 - IF ANY NEGATIVE VALUES REMAIN (VERY UNCOMMON)!
!              PERFORM GLOBAL CORRECTION AS IN METHOD 1.
!
!              IN BOTH METHODS IF THERE ARE INSUFFICIENT
!              POSITIVE VALUES WITHIN THE WHOLE LAYER AN ERROR
!              MESSAGE IS PRODUCED.  PROGRAM CONTINUES IF L_NEG_Q
!              IS TRUE IN WHICH CASE q IS NOT CONSERVED.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering


      USE domain_params
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      
      implicit none

! Arguments with Intent IN. ie: Input variables.

      integer                                                           &
     &  row_length                                                      &
                            ! number of points per row
     &, rows                                                            &
     &, q_levels                                                        &
                            ! number of moist levels on this processor.
     &, model_domain

! Arguments with Intent IN/OUT.

      real                                                              &
     & q(row_length*rows,q_levels)    ! q field.

! Local Variables.

      integer                                                           &
     &  i,j,k                                                           &
                           ! loop variables
     &, no_neg                                                          &
                           ! count of number negative in step 3
     &, no_neg2                                                         &
                           ! count of number negative in step 3
     &, P_FIELD,II,MAX_SEARCH                                           &
     &, JE,JW,JJ,NPT,NN,IN,IS,POINTER,ROW                               &
     &, point(row_length*rows),l

      real                                                              &
     &  sum_positive                                                    &
                           ! global sum of positive q
     &, sum_negative                                                    &
                           ! global sum of negative q
     &, SUM_NEIGHBOURS                                                  &
     &, factor                                                          &
                           ! factor used to rescale neighbours of q
     &, TEMP1,TEMP2                                                     &
     &, points             ! total number of points real(P_field)

      logical                                                           &
     &  FAIL                                                            &
     &, l_q_pos_local                                                   &
                           ! logical : true to do Method 2,
                           !         : false to do Method 1
     &, L_SP,L_NP                                                       &
     &, l_reset                                                         &
     &, L_local_fail

      real                                                              &
     &  qlimit_in                                                       &
                        !    lowest allowed value of q
     &, qlimit          !    lowest allowed value of q

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('Q_POS',zhook_in,zhook_handle)
      P_FIELD = ROWS * ROW_LENGTH
      points=real(P_FIELD)
      MAX_SEARCH=4

!---------------------------------------------------------------------
!    BLOCK FOR GLOBAL MODELS (has extra code for poles)
!---------------------------------------------------------------------
      if(model_domain == mt_global)then

      DO K= 1,Q_LEVELS
        FAIL=.FALSE.
        l_reset=.false.

        IF (L_q_POS_LOCAL) THEN
          no_neg=0
!---------------------------------------------------------------------
!    METHOD 2   ! local search
!---------------------------------------------------------------------

!   Loop round non-polar points
          DO I = ROW_LENGTH+1, P_FIELD-ROW_LENGTH
            IF(q(I,K)  <   qlimit_in) THEN
              no_neg = no_neg+1
              point(no_neg)=i
            endif
          end do
! rescale q if there are negative values in the field.
          if((no_neg >  0 .or. q(1,k) <  qlimit_in .or.                 &
     &      q(p_field,k) <  qlimit_in) .and. qlimit_in /=  0.0 )then
            q(:,k)=q(:,k)-qlimit_in
            l_reset =.true.
          endif

!---------------------------------------------------------------------
!    SECTION 2. THIS DOES STEP 1. WHERE q IS NEGATIVE LOOP ROUND
!              8 NEIGHBOURS TO PERFORM LOCAL CORRECTION.  THIS OCCURS
!              AT APPROXIMATELY 2% OF THE POINTS.  THE FOLLOWING CODE
!              IN SECTIONS 2 & 3 IS RATHER LONG-WINDED, BUT CPU IS
!              SIGNIFICANTLY REDUCED OVER A SIMPLER APPROACH.
!---------------------------------------------------------------------
! loop over negative points
          do l=1,no_neg
            i=point(l)
            qlimit=0.
!   If q is negative find 8 neighbours ...
            if(q(i,k) < 0.0)then
            ROW=(I-1)/ROW_LENGTH +1
            NPT=I-(ROW-1)*ROW_LENGTH
            JE=1
            JW=-1
            IF(NPT == ROW_LENGTH) JE=1-ROW_LENGTH
            IF(NPT == 1) JW=ROW_LENGTH-1
!   ... and sum neighbouring 8 values
            SUM_NEIGHBOURS = q(I+JE,K) + q(I+JW,K)
            IF(ROW >  2) THEN
              II=I-ROW_LENGTH
              SUM_NEIGHBOURS = SUM_NEIGHBOURS +                         &
     &                         q(II,K) + q(II+JE,K) + q(II+JW,K)
            END IF
            IF(ROW <  ROWS-1) THEN
              II=I+ROW_LENGTH
              SUM_NEIGHBOURS = SUM_NEIGHBOURS +                         &
     &                         q(II,K) + q(II+JE,K) + q(II+JW,K)
            END IF
! AJM pole section
            IF(ROW == 2) THEN
              II=I-ROW_LENGTH
              SUM_NEIGHBOURS = SUM_NEIGHBOURS +  q(II,K)
            END IF
            IF(ROW == ROWS-1) THEN
              II=I+ROW_LENGTH
              SUM_NEIGHBOURS = SUM_NEIGHBOURS +  q(II,K)
            END IF
!
!   If possible set q=0 at point with negative value and adjust
!   neighbouring values by FACTOR to conserve q
            IF(SUM_NEIGHBOURS +q(I,k)  >=  0.0) THEN
              FACTOR = 1.0 + q(I,K) / SUM_NEIGHBOURS
              IF(FACTOR >= 0.0) THEN
                q(I,K) = qlimit
                q(I+JE,K) = q(I+JE,K) * FACTOR
                q(I+JW,K) = q(I+JW,K) * FACTOR
                IF(ROW  >   2) THEN
                  II=I-ROW_LENGTH
                  q(II,K) = q(II,K)*FACTOR
                  q(II+JE,K) = q(II+JE,K) * FACTOR
                  q(II+JW,K) = q(II+JW,K) * FACTOR
                END IF
                IF(ROW <  ROWS-1) THEN
                  II=I+ROW_LENGTH
                  q(II,K) = q(II,K)*FACTOR
                  q(II+JE,K) = q(II+JE,K)*FACTOR
                  q(II+JW,K) = q(II+JW,K)*FACTOR
                END IF
! AJM pole section
                IF(ROW == 2) THEN
                  II=I-ROW_LENGTH
                  q(II,K) = q(II,K)*FACTOR
                  do JJ=1,row_length
                    q(JJ,K)=q(II,K)
                  end do
                END IF
                IF(ROW == ROWS-1) THEN
                  II=I+ROW_LENGTH
                  q(II,K) = q(II,K)*FACTOR
                  do jj=P_FIELD-row_length+1,P_field
                    q(JJ,K)=q(II,K)
                  end do
                END IF
              END IF    !factor
            END IF      !sum_neighbours >= 0.0
            end if   ! point still negative
          end do        ! loop over negative points (l)
! find points which are still negative
          no_neg2=0
          do i=1,no_neg
            if (q(point(i),k) <  0.0)then
              no_neg2=no_neg2+1
              point(no_neg2)=point(i)
            endif
          end do

!
!---------------------------------------------------------------------
!    SECTION 3. THIS DOES STEP 2. WHERE q IS NEGATIVE AND STEP 1 
!            CORRECTIONS HAVE FAILED, ATTEMPT LOCAL CORRECTION BY 
!            EXTENDING SEARCH TO 24, 48 OR 80 SURROUNDING POINTS (SEARCH 
!            RADIUS EQUALS 2, 3 OR 4), OR FURTHER. THIS OCCURS LESS THAN
!            0.1% OF THE TIME.  ALTHOUGH SMALL IN NUMBER, LOCAL
!            CORRECTIONS AT THESE POINTS ARE STILL REQUIRED TO AVOID
!            SLOW CLIMATE DRIFT.
!---------------------------------------------------------------------
          do l=1,no_neg2
            i=point(l)
            if(q(i,k) < 0.0)then
            DO NN=2,MAX_SEARCH
              ROW=(I-1)/ROW_LENGTH +1
              NPT=I-(ROW-1)*ROW_LENGTH
              L_SP=.false.
              L_NP=.false.
              JW=NPT-NN
              JE=NPT+NN
              IS=ROW-NN
              IN=ROW+NN
              IF(IS <  1) IS=1
              IF(IN >= ROWS+1) IN=ROWS
              if(IS == 1) then
                L_SP=.true.
                IS=2
              end if
              if(IN == rows) then
                L_NP=.true.
                IN=rows-1
              end if
!   ... and sum q at neighbouring points
!  the loops mean we include the centre point so set start value to
!  minus centre point
              SUM_NEIGHBOURS=-q(I,K)
              DO II=IS,IN
                DO JJ=JW,JE
                  POINTER = (II-1)*ROW_LENGTH + JJ
                  IF(JJ >  ROW_LENGTH) POINTER = POINTER-ROW_LENGTH
                  IF(JJ <  1) POINTER = POINTER+ROW_LENGTH
                  SUM_NEIGHBOURS = SUM_NEIGHBOURS + q(POINTER,K)
                END DO
              END DO
              if(l_SP)then
                SUM_NEIGHBOURS = SUM_NEIGHBOURS + q(npt,k)
              endif
              if(l_NP)then
                SUM_NEIGHBOURS = SUM_NEIGHBOURS +                       &
     &                           q(p_field-row_length+npt,k)
              endif

!
!   If possible set q=0 at point with negative value and adjust
!   neighbouring values by FACTOR to conserve q
              IF(SUM_NEIGHBOURS +q(i,k)  >=  0.0) THEN
                FACTOR = 1.0 + q(I,K) / SUM_NEIGHBOURS
                IF(FACTOR >= 0.0) THEN
                  DO II=IS,IN
                    DO JJ=JW,JE
                      POINTER = (II-1)*ROW_LENGTH + JJ
                      IF(JJ >  ROW_LENGTH)POINTER=POINTER-ROW_LENGTH
                      IF(JJ <  1) POINTER = POINTER+ROW_LENGTH
                      q(POINTER,K) = q(POINTER,K)*FACTOR
                    END DO
                  END DO
                  if(l_SP)then
                    q(npt,k)=q(npt,k)*FACTOR
                    do jj=1,row_length
                      q(jj,k)=q(npt,k)
                    end do
                  endif
                  if(l_NP)then
                    q(p_field-row_length+npt,k)=                        &
     &                   q(p_field-row_length+npt,k)*factor
                    do jj=p_field-row_length+1,p_field
                      q(jj,k)=q(p_field-row_length+npt,k)
                    end do
                  endif
                  q(I,K) = qlimit
                  GOTO 300    ! Successful correction performed
                END IF      ! factor
              END IF        !sum_neighbours
            END DO     ! End loop on search radius (NN)
!
            FAIL=.TRUE.
  300       CONTINUE
            end if   ! point still negative
          END DO         ! End loop on points (l)

!
!  Loop round polar points
!
!---------------------------------------------------------------------
!    SECTION 4.  THIS DOES STEP 2 FOR THE POLE POINTS WORKING ONE ROW 
!                AT A TIME.
!---------------------------------------------------------------------
!  Check for negative q at south pole
          IF(q(1,K)  <   0.0) THEN
            SUM_NEIGHBOURS = 0.0
            DO I = 1, ROW_LENGTH
              SUM_NEIGHBOURS = SUM_NEIGHBOURS + q(ROW_LENGTH+I,K)
            END DO
            IF(SUM_NEIGHBOURS+q(1,k) >= 0.0) THEN
              FACTOR = 1.0 + q(1,K) / SUM_NEIGHBOURS
              IF(FACTOR >= 0.0) THEN
                DO I = 1, ROW_LENGTH
                  q(I,K) = 0.0
                  q(ROW_LENGTH+I,K) = q(ROW_LENGTH+I,K)*FACTOR
                END DO
                GOTO 400    ! Correction performed
              END IF
            END IF
            FAIL=.TRUE.
          END IF
  400     CONTINUE
!  Check for negative q at north pole
          IF(q(P_FIELD-1,K) <  0.0) THEN
            SUM_NEIGHBOURS = 0.0
            DO I = 1, ROW_LENGTH
              SUM_NEIGHBOURS = SUM_NEIGHBOURS +                         &
     &                           q(P_FIELD+1-ROW_LENGTH-I,K)
            END DO
            IF(SUM_NEIGHBOURS+q(P_FIELD-1,K) >                          &
     &                 0.0) THEN
              FACTOR = 1.0 + q(P_FIELD-1,K)  /SUM_NEIGHBOURS
              IF(FACTOR >= 0.0) THEN
                DO I=1,ROW_LENGTH
                  q(P_FIELD-ROW_LENGTH+I,K) = 0.0
                  q(P_FIELD+1-ROW_LENGTH-I,K) =                         &
     &                       q(P_FIELD+1-ROW_LENGTH-I,K) * FACTOR
                END DO
                GOTO 500    ! Correction performed
              END IF
            END IF
            FAIL=.TRUE.
          END IF
  500     CONTINUE
!  RESCALE q by adding q_limit_in if necessary
          if(L_reset)q(:,k)=q(:,k)+qlimit_in
        END IF          !  on L_q_POS_LOCAL  = true

!---------------------------------------------------------------------
!   SECTION 5.  METHOD 1. (ALSO METHOD 2 STEP 3, IF METHOD 2 HAS FAILED)
!              REMOVE REMAINING NEGATIVE q BY SUMMING
!              NEGATIVE AND POSITIVE VALUES SEPARATELY ON THE
!              LEVEL AND PERFORMING A GLOBAL CORRECTION.
!---------------------------------------------------------------------
        L_LOCAL_FAIL=FAIL
        L_reset=.false.
        IF(FAIL .OR. .NOT. L_q_POS_LOCAL) THEN
          IF(qlimit_in /= 0.0)THEN
            q(:,k)=q(:,k)-qlimit_in
            l_reset =.true.
          END IF
          NO_NEG = 0            !zero count of -ve points
          SUM_POSITIVE = 0.0
          SUM_NEGATIVE = 0.0
!         include one point off each polar row in loop
          DO I=row_length,P_field-row_length+1
            IF(q(I,K)  <  0.0 )THEN
              NO_NEG = NO_NEG+1           !count number of -ve points
              SUM_NEGATIVE = SUM_NEGATIVE + q(I,K)
              point(no_neg)=i
            ELSE
              SUM_POSITIVE = SUM_POSITIVE + q(I,K)
            END IF
          END DO

! If a negative value is found re-scale all positive points
!
          IF (NO_NEG  >  0) THEN
            FAIL=.TRUE.
            if(sum_positive+sum_negative  >= 0.0                        &
     &          .and. SUM_POSITIVE /= 0.0) then
              factor =1. +(sum_negative / sum_positive ) 
            else
              factor=-1.0
            endif
            IF(FACTOR <  0)THEN
      WRITE(6,*)'WARNING q_POS : UNABLE TO RESET VALUES CONSERVATIVELY'
      WRITE(6,*)'WARNING q_POS :',NO_NEG,' points were less than ',     &
     &        qlimit_in, 'and have been reset to ',qlimit_in
      WRITE(6,*)'WARNING q_POS : All other points unchanged'
              FACTOR = 1.
            ENDIF
          ELSE
            FAIL=.FALSE.
          END IF
        END IF   ! FAIL .OR. .NOT. L_q_POS_LOCAL
!
!  To rescale (i.e. conserve?)  after adding back lost moisture
        IF (FAIL) THEN
          DO I= row_length,P_field-row_length+1
            q(I,K) = q(I,K) * FACTOR
          END DO
          do i=1,no_neg
            q(point(i),k)=0.0
          end do
        ENDIF
!
!---------------------------------------------------------------------
!    SECTION 6.     IN GLOBAL MODEL, IF METHOD 1 USED
!                   SET ALL POLAR VALUES TO BE THE SAME
!---------------------------------------------------------------------

        IF (L_LOCAL_FAIL .or. .NOT.L_q_POS_LOCAL) THEN
!  Set polar points to equal values
          DO I=1,ROW_LENGTH-1
            q(I,K) = q(ROW_LENGTH,K)
            q(P_FIELD+1-I,K) = q(P_FIELD+1-ROW_LENGTH,K)
          END DO
        END IF

        IF (l_reset)q(:,k)=q(:,k)+qlimit_in

! END LOOP OVER LEVELS
      END DO

!---------------------------------------------------------------------
!    BLOCK FOR LIMITED AREA MODELS 
!---------------------------------------------------------------------
      elseif(model_domain == mt_lam)then

      DO K= 1,Q_LEVELS
        FAIL=.FALSE.
        l_reset=.false.

        IF (L_q_POS_LOCAL) THEN
          no_neg=0
!---------------------------------------------------------------------
!    METHOD 2  ! local search
!---------------------------------------------------------------------

! find points which need resetting
          DO I = 1, P_FIELD
            IF(q(I,K)  <   qlimit_in) THEN
              no_neg = no_neg+1
              point(no_neg)=i
            endif
          end do
! reset field so target is zero not qlimit
          if(no_neg >  0 .and. qlimit_in /= 0.0)then
            q(:,k)=q(:,k)-qlimit_in
            l_reset =.true.
          endif
!
!---------------------------------------------------------------------
!    SECTION 2.  THIS DOES STEP 1. WHERE q IS NEGATIVE LOOP ROUND
!              8 NEIGHBOURS TO PERFORM LOCAL CORRECTION.  THIS OCCURS
!              AT APPROXIMATELY 2% OF THE POINTS.  THE FOLLOWING CODE
!              IN SECTIONS 2 & 3 IS RATHER LONG-WINDED, BUT CPU IS
!              SIGNIFICANTLY REDUCED OVER A SIMPLER APPROACH.
!---------------------------------------------------------------------
! loop over negative points
          do l=1,no_neg
            i=point(l)
            qlimit=0.
!   If q is negative find 8 neighbours ...
            if(q(i,k) < 0.0)then
            ROW=(I-1)/ROW_LENGTH +1
            NPT=I-(ROW-1)*ROW_LENGTH
            JE=NPT+1
            JW=NPT-1
            IS=ROW-1
            IN=ROW+1
            IF(JE >  ROW_LENGTH) JE=ROW_LENGTH
            IF(JW <  1) JW=1
            IF(IS <  1) IS=1
            IF(IN >  rows) IN= rows
!   ... and sum neighbouring 8 values
            sum_neighbours =0.0
            do ii=is,in
              do jj=jw,je
                j=(ii-1)*row_length +jj
                SUM_NEIGHBOURS = SUM_NEIGHBOURS +q(J,K)
              end do
            end do
!
!   If possible set q=0 at point with negative value and adjust
!   neighbouring values by FACTOR to conserve q
            IF(SUM_NEIGHBOURS  >=  0.0) THEN
              FACTOR = 1.0 + q(I,K) / (SUM_NEIGHBOURS-q(i,k))
              do ii=is,in
                do jj=jw,je
                  j=(ii-1)*row_length +jj
                  q(j,K) = q(j,K)*FACTOR
                end do
              end do
              q(I,K) = qlimit
            END IF
            end if   ! point still negative
          end do  ! loop over negative points
! find points which are still negative
          no_neg2=0
          do i=1,no_neg
            if (q(point(i),k) <  0.0)then
              no_neg2=no_neg2+1
              point(no_neg2)=point(i)
            endif
          end do
!
!---------------------------------------------------------------------
!    SECTION 3. THIS DOES STEP 2. WHERE q IS NEGATIVE AND STEP 1 
!            CORRECTIONS HAVE FAILED, ATTEMPT LOCAL CORRECTION BY 
!            EXTENDING SEARCH TO 24, 48 OR 80 SURROUNDING POINTS (SEARCH 
!            RADIUS EQUALS 2, 3 OR 4), OR FURTHER. THIS OCCURS LESS THAN
!            0.1% OF THE TIME.  ALTHOUGH SMALL IN NUMBER, LOCAL
!            CORRECTIONS AT THESE POINTS ARE STILL REQUIRED TO AVOID
!            SLOW CLIMATE DRIFT.
!---------------------------------------------------------------------
          do l=1,no_neg2
            i=point(l)
            if(q(i,k) < 0.0)then
            DO NN=2,MAX_SEARCH
              ROW=(I-1)/ROW_LENGTH +1
              NPT=I-(ROW-1)*ROW_LENGTH
              JW=NPT-NN
              JE=NPT+NN
              IS=ROW-NN
              IN=ROW+NN
              IF(JE >  ROW_LENGTH) JE=ROW_LENGTH
              IF(JW <  1) JW=1
              IF(IS <  1) IS=1
              IF(IN >  rows) IN= rows
!   ... and sum q at neighbouring points
! the loops mean we include the centre point. Unlike the global code,
! where we set start value to minus centre point, we set the start value 
! to zero and adjust sum_neighbours where required.
              sum_neighbours =0.0
              do ii=is,in
                do jj=jw,je
                  j=(ii-1)*row_length +jj
                  SUM_NEIGHBOURS = SUM_NEIGHBOURS +q(J,K)
                end do
              end do
!
!   If possible set q=0 at point with negative value and adjust
!   neighbouring values by FACTOR to conserve q.  
              IF(SUM_NEIGHBOURS  >=  0.0) THEN
                FACTOR = 1.0 + q(I,K) / (SUM_NEIGHBOURS-q(i,k))
                do ii=is,in
                  do jj=jw,je
                    j=(ii-1)*row_length +jj
                    q(j,K) = q(j,K)*FACTOR
                  end do
                end do
                q(I,K) = qlimit
                GOTO 301    ! Successful correction performed
              END IF
            END DO     ! End loop on search radius (NN)
            FAIL=.TRUE.
  301       CONTINUE
            end if   ! point still negative
          END DO         ! End loop on points (I)
! rescale q by addingf q_limit_in if necessary
          if(L_reset)q(:,k)=q(:,k)+qlimit_in
        END IF          !  on L_q_POS_LOCAL  = true

!---------------------------------------------------------------------
!    SECTION 5.  METHOD 1.(ALSO METHOD 2 STEP 3, IF METHOD 2 HAS FAILED)
!              REMOVE REMAINING NEGATIVE q BY SUMMING
!              NEGATIVE AND POSITIVE VALUES SEPARATELY ON THE
!              LEVEL AND PERFORMING A GLOBAL CORRECTION.
!---------------------------------------------------------------------
        L_reset=.false.
        IF(FAIL .OR. .NOT. L_q_POS_LOCAL) THEN
          IF(qlimit_in /= 0.0)THEN 
            q(:,k)=q(:,k)-qlimit_in
            l_reset =.true.
          ENDIF
          NO_NEG = 0            !zero count of -ve points
          SUM_POSITIVE = 0.0
          SUM_NEGATIVE = 0.0
          DO I=1, p_field
            IF(q(I,K)  <  0.0)THEN
              NO_NEG = NO_NEG+1           !count number of -ve points
              SUM_NEGATIVE = SUM_NEGATIVE + q(I,K)
              point(no_neg)=i
            ELSE
              SUM_POSITIVE = SUM_POSITIVE + q(I,K)
            END IF
          END DO

! If a negative value is found re-scale all positive points
!
          IF (NO_NEG  >  0) THEN
            FAIL=.TRUE.
            if(sum_positive+sum_negative  >=  0.0                       &
     &         .and. SUM_POSITIVE /= 0.0) then
              factor =1. +(sum_negative / sum_positive)
            else
              factor=-1.0
            endif
            IF(FACTOR <  0)THEN
      WRITE(6,*)'WARNING q_POS : UNABLE TO RESET VALUES CONSERVATIVELY'
      WRITE(6,*)'WARNING q_POS :',NO_NEG,' points were less than ',     &
     &        qlimit_in, 'and have been reset to ',qlimit_in
      WRITE(6,*)'WARNING q_POS : All other points unchanged'
              FACTOR = 1.
            ENDIF
          ELSE
            FAIL=.FALSE.
          END IF
        END IF
!
!  To rescale (i.e. conserve?)  after adding back lost moisture
        IF (FAIL) THEN
          DO I= 1,P_FIELD
            q(I,K) = q(I,K) * FACTOR
          END DO
          do i=1,no_neg
            q(point(i),k)=0.0
          end do
        ENDIF

        IF(l_reset)q(:,k)=q(:,k)+qlimit_in

! END LOOP OVER LEVELS
      END DO

      else  ! not mt_global or mt_lam
        write(6,*)'WARNING : QPOS NOT CODED FOR THIS DOMAIN.'
        write(6,*)'ORIGINAL FIELD RETURNED'
      end if

      IF (lhook) CALL dr_hook('Q_POS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE q_pos

