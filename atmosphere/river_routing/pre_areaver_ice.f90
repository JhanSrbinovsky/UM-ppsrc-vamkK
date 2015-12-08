! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: River Routing
MODULE pre_areaver_ice_mod

IMPLICIT NONE

CONTAINS


      SUBROUTINE PRE_AREAVER_ICE(GAPS_LAMBDA_SRCE,LAMBDA_SRCE           &
     &,GAPS_PHI_SRCE,PHI_SRCE,CYCLIC_SRCE,LROW_SRCE,WANT,MASK_SRCE      &
     &,GAPS_LAMBDA_TARG,LAMBDA_TARG,GAPS_PHI_TARG,PHI_TARG              &
     &,CYCLIC_TARG,SPHERICAL,FRAC_SRCE                                  &
     &,MAXL,COUNT_TARG,BASE_TARG,INDEX_SRCE,WEIGHT,ICODE,CMESSAGE)
!
!  Subroutine PRE_AREAVER_ICE
!
!  ****THIS IS A SLIGHTLY MODIFIED VERSION OF PRE_AREAVER - ANY
!  CHANGES MADE TO PRE_AREAVER SHOULD BE MIRRORED IN THIS ROUTINE****
!
!  Calculate weights for area-averaging data on the source grid to
!  data on the target grid taking into account the ice or leads area.
!  Call includes extra argument FRAC_SRCE
!
!  For details of this program, see PRE_AREAVER.
!
      USE conversions_mod, ONLY: pi_over_180
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
      
      INTEGER                                                           &
     & GAPS_LAMBDA_SRCE                                                 &
                               !IN number of lambda gaps in source grid
     &,GAPS_PHI_SRCE                                                    &
                               !IN number of phi gaps in source grid
     &,GAPS_LAMBDA_TARG                                                 &
                               !IN number of lambda gaps in target grid
     &,GAPS_PHI_TARG                                                    &
                               !IN number of phi gaps in target grid
     &,LROW_SRCE                                                        &
                               !IN first dimension of MASK_SRCE
     &,MAXL                                                             &
                               !INOUT maximum entries in output lists
     &,COUNT_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                       &
!                              !OUT no. of source boxes in target box
     &,BASE_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                        &
!                              !OUT first index in list for target box
     &,INDEX_SRCE(MAXL)                                                 &
                               !OUT list of source box indices
     &,ICODE                   !OUT return code
      LOGICAL                                                           &
     & CYCLIC_SRCE                                                      &
                               !IN source grid is cyclic
     &,CYCLIC_TARG                                                      &
                               !IN target grid is cyclic
     &,WANT                                                             &
                               !IN indicator of wanted points in mask
     &,MASK_SRCE(LROW_SRCE,*)                                           &
                               !IN land/sea mask for source grid
     &,SPHERICAL               !IN calculate weights for a sphere
      REAL                                                              &
     & LAMBDA_SRCE(*)                                                   &
                               !IN source lambda line coordinates
     &,PHI_SRCE(*)                                                      &
                               !IN source phi line coordinates
     &,LAMBDA_TARG(*)                                                   &
                               !IN target lambda line coordinates
     &,PHI_TARG(*)                                                      &
                               !IN target phi line coordinates
     &,FRAC_SRCE(LROW_SRCE,*)                                           &
                               !IN source ice or leads fraction
     &,WEIGHT(MAXL)            !OUT list of weights for source boxes
      CHARACTER(LEN=80) CMESSAGE  !OUT error message
!
      INTEGER                                                           &
     & LINES_LAMBDA_SRCE                                                &
                               ! number of source lambda lines
     &,LINES_PHI_SRCE                                                   &
                               ! number of source phi lines
     &,LINES_LAMBDA_TARG                                                &
                               ! number of target lambda lines
     &,LINES_PHI_TARG                                                   &
                               ! number of target phi lines
     &,COUNT_LAMBDA(GAPS_LAMBDA_TARG)                                   &
!                              ! number of source lambda gaps per target
     &,COUNT_PHI(GAPS_PHI_TARG)                                         &
!                              ! number of source phi gaps per target
     &,INDEX_LAMBDA(GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)                  &
!                              ! source lambda gap indices
     &,INDEX_PHI(GAPS_PHI_SRCE+GAPS_PHI_TARG)                           &
!                              ! source phi gap indices
     &,IX1,IY1,IX2,IY2                                                  &
                               ! working SRCE/TARG LAMBDA/PHI indices
     &,IX1N,IX1W                                                        &
                               ! working indices
     &,IXL(GAPS_LAMBDA_TARG+1)                                          &
                               ! source line on the left of target line
     &,IX2N                                                             &
                               ! target gap to the right of IX2
     &,IYT(GAPS_PHI_TARG+1)                                             &
                               ! source line above target line
     &,IXP,IYP,IP                                                       &
                               ! pointers into lists
     &,IX,IY,I                 ! loop indices
      REAL                                                              &
     & BASLAM                                                           &
                               ! minimum lambda for TEMP coordinates
     &,BTARG                                                            &
                               ! width of target gap
     &,BLO,BHI                                                          &
                               ! limits of gap
     &,TEMP_SRCE(GAPS_LAMBDA_SRCE+1)                                    &
!                              ! adjusted version of LAMBDA_SRCE
     &,TEMP_TARG(GAPS_LAMBDA_TARG+1)                                    &
!                              ! adjusted version of LAMBDA_TARG
     &,FRAC_LAMBDA(GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)                   &
!                              ! fractions of target lambda gaps
     &,FRAC_PHI(GAPS_PHI_SRCE+GAPS_PHI_TARG)                            &
!                              ! fractions of target phi gaps
     &,SUM                     ! sum of weights

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!   1  Set up the lambda coordinates to make them easier to handle.
!
!   1.1  Produce in TEMP_SRCE a monotonically increasing set of angles
!   equivalent to LAMBDA_SRCE i.e. equal under modulo 360.
!
      IF (lhook) CALL dr_hook('PRE_AREAVER_ICE',zhook_in,zhook_handle)
      IF (CYCLIC_SRCE) THEN
        LINES_LAMBDA_SRCE=GAPS_LAMBDA_SRCE
      ELSE
        LINES_LAMBDA_SRCE=GAPS_LAMBDA_SRCE+1
      ENDIF
      BASLAM=LAMBDA_SRCE(1)
      DO IX1=1,LINES_LAMBDA_SRCE
        IF (LAMBDA_SRCE(IX1) <  BASLAM) THEN
          TEMP_SRCE(IX1)=LAMBDA_SRCE(IX1)+360.
        ELSE
          TEMP_SRCE(IX1)=LAMBDA_SRCE(IX1)
        ENDIF
      END DO ! Labelled DO for the sake of fpp
!
!   1.2  Produce in TEMP_TARG a set of angles equivalent to
!   LAMBDA_TARG i.e. equal under modulo 360, but all in the range
!   BASLAM to BASLAM+360, where BASLAM=min(TEMP_LAMBDA).
!
      IF (CYCLIC_TARG) THEN
        LINES_LAMBDA_TARG=GAPS_LAMBDA_TARG
      ELSE
        LINES_LAMBDA_TARG=GAPS_LAMBDA_TARG+1
      ENDIF
      DO IX2=1,LINES_LAMBDA_TARG
        TEMP_TARG(IX2)=MOD(LAMBDA_TARG(IX2)-BASLAM,360.)
        IF (TEMP_TARG(IX2) <  0.) TEMP_TARG(IX2)=TEMP_TARG(IX2)+360.
        TEMP_TARG(IX2)=TEMP_TARG(IX2)+BASLAM
      END DO ! Labelled DO for the sake of fpp
!
!   2  For each target lambda line, find the index of the next source
!   lambda line to the left.
!
      DO IX2=1,LINES_LAMBDA_TARG
        DO IX1=1,LINES_LAMBDA_SRCE
          IF (TEMP_TARG(IX2) >= TEMP_SRCE(IX1)) IXL(IX2)=IX1
        ENDDO
      ENDDO
!
!   3  Find which source lambda gaps cover each target gap and the
!   fractions they contribute.
!
!     At this point IXL(target_line) gives the index of the next source
!     lambda line to the left of the target lambda line, wrapping round
!     if the source grid is cyclic. This is also the index of the source
!     gap in which the target line falls. Similarly, the index of the
!     target line is also that of the target gap of which it is the
!     left-hand limit. Therefore also IXL(target_gap+1, wrapping round
!     if reqd.), is the index of the source gap which contains the
!     right-hand limit of the target gap. For each target gap, we loop
!     over all source gaps and find the fraction covered by each. Record
!     the fraction and the source index in cumulative lists. If the
!     source grid is not cyclic, parts of the target gap lying outside
!     the source grid are neglected.
!
      IXP=0
      DO IX2=1,GAPS_LAMBDA_TARG
        IX=0
        IX2N=MOD(IX2,LINES_LAMBDA_TARG)+1
        BTARG=TEMP_TARG(IX2N)-TEMP_TARG(IX2)
        IF (BTARG <  0.) THEN
          BTARG=BTARG+360.
          IX1N=IXL(IX2N)+LINES_LAMBDA_SRCE
        ELSE
          IX1N=IXL(IX2N)
        ENDIF
        DO IX1W=IXL(IX2),IX1N
          IX1=MOD(IX1W-1,LINES_LAMBDA_SRCE)+1
          IF (CYCLIC_SRCE.OR.IX1 /= LINES_LAMBDA_SRCE) THEN
            IF (IX1W == IXL(IX2)) THEN
              BLO=TEMP_TARG(IX2)
            ELSE
              BLO=TEMP_SRCE(IX1)
            ENDIF
            IF (IX1W == IX1N) THEN
              BHI=TEMP_TARG(IX2N)
            ELSE
              BHI=TEMP_SRCE(MOD(IX1,LINES_LAMBDA_SRCE)+1)
            ENDIF
            IF (BHI <  BLO) BHI=BHI+360.
            IF (ABS(BHI-BLO) >  1E-7*ABS(BTARG)) THEN
              IX=IX+1
              INDEX_LAMBDA(IXP+IX)=IX1
              FRAC_LAMBDA(IXP+IX)=(BHI-BLO)/BTARG
            ENDIF
          ENDIF
        ENDDO
        COUNT_LAMBDA(IX2)=IX
        IXP=IXP+COUNT_LAMBDA(IX2)
      ENDDO
!
!   4  For each target phi line, find the index of the next source phi
!   line above. Comments as for section 2, without wrap-round. Note
!   that this assumes that the atmosphere gridbox ordering is from
!   south to north; a north to south atmosphere would require LE
!   instead of GE in the following IF test.
!
      LINES_PHI_SRCE=GAPS_PHI_SRCE+1
      LINES_PHI_TARG=GAPS_PHI_TARG+1
      DO IY2=1,LINES_PHI_TARG
        IYT(IY2)=0
        DO IY1=1,LINES_PHI_SRCE
          IF (PHI_TARG(IY2) >= PHI_SRCE(IY1)) IYT(IY2)=IY1
        ENDDO
      ENDDO
!
!   5  Find which source phi gaps cover each target gap and the
!   fractions they contribute. Comments as for section 3, without
!   wrap-round.
!
      IYP=0
      DO IY2=1,GAPS_PHI_TARG
        IY=0
        IF (SPHERICAL) THEN
!     Contain angle between +-90. There is no real area outside
!     these limits on a sphere.
          BTARG=SIN(MAX(MIN(PHI_TARG(IY2+1),90.),-90.)*PI_OVER_180)     &
     &    -SIN(MAX(MIN(PHI_TARG(IY2),90.),-90.)*PI_OVER_180)
        ELSE
           BTARG=PHI_TARG(IY2+1)-PHI_TARG(IY2)
        ENDIF
        DO IY1=IYT(IY2),IYT(IY2+1)
          IF (.NOT.(IY1 == 0.OR.IY1 == LINES_PHI_SRCE)) THEN
            IF (IY1 == IYT(IY2)) THEN
              BLO=PHI_TARG(IY2)
            ELSE
              BLO=PHI_SRCE(IY1)
            ENDIF
            IF (IY1 == IYT(IY2+1)) THEN
              BHI=PHI_TARG(IY2+1)
            ELSE
              BHI=PHI_SRCE(IY1+1)
            ENDIF
            IF (SPHERICAL) THEN
              BLO=MAX(MIN(BLO,90.),-90.)
              BHI=MAX(MIN(BHI,90.),-90.)
            ENDIF
            IF (ABS(BHI-BLO) >  1E-7*ABS(BTARG)) THEN
              IY=IY+1
              INDEX_PHI(IYP+IY)=IY1
!             Both numerator and denominator in the following are -ve.
              IF (SPHERICAL) THEN
                FRAC_PHI(IYP+IY)                                        &
     &          =(SIN(BHI*PI_OVER_180)-SIN(BLO*PI_OVER_180))/BTARG
              ELSE
                FRAC_PHI(IYP+IY)=(BHI-BLO)/BTARG
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        COUNT_PHI(IY2)=IY
        IYP=IYP+COUNT_PHI(IY2)
      ENDDO
!
!   6  For each target box, loop over contributing source boxes and
!   calculate the weights for each one, ignoring land boxes.
!
!     After the first pass for each target box, go back and normalise
!     the weights to compensate for land source boxes and any outside
!     the source grid. Record the source box index and the weight in
!     cumulative lists.
!
      IP=0
      IYP=0
      DO IY2=1,GAPS_PHI_TARG
        IXP=0
        DO IX2=1,GAPS_LAMBDA_TARG
          I=0
          SUM=0
          DO IY=IYP+1,IYP+COUNT_PHI(IY2)
            DO IX=IXP+1,IXP+COUNT_LAMBDA(IX2)
              IF (MASK_SRCE(INDEX_LAMBDA(IX),INDEX_PHI(IY))             &
     &        .EQV.WANT) THEN
                I=I+1
                INDEX_SRCE(IP+I)=INDEX_LAMBDA(IX)                       &
     &          +(INDEX_PHI(IY)-1)*GAPS_LAMBDA_SRCE
                WEIGHT(IP+I)=FRAC_LAMBDA(IX)*FRAC_PHI(IY)*              &
     &                    FRAC_SRCE(INDEX_LAMBDA(IX),INDEX_PHI(IY))
                SUM=SUM+WEIGHT(IP+I)
              ENDIF
            ENDDO
          ENDDO
          COUNT_TARG(IX2,IY2)=I
          BASE_TARG(IX2,IY2)=IP
          DO I=1,COUNT_TARG(IX2,IY2)
            IF(SUM /=  0.) THEN
              WEIGHT(IP+I)=WEIGHT(IP+I)/SUM
            ELSE
              WEIGHT(IP+I)=0.
            ENDIF
          ENDDO
          IP=IP+COUNT_TARG(IX2,IY2)
          IXP=IXP+COUNT_LAMBDA(IX2)
        ENDDO
        IYP=IYP+COUNT_PHI(IY2)
      ENDDO
      MAXL=IP
!
      ICODE=0
      CMESSAGE=' '
      IF (lhook) CALL dr_hook('PRE_AREAVER_ICE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE PRE_AREAVER_ICE
END MODULE pre_areaver_ice_mod
