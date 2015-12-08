! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE con_rad_4a5a
!
!  PURPOSE : Calculates convective cloud top, base and amount
!
!  Suitable for single column model use
!
!  PROGRAMMING STANDARDS : Unified Model Documentation Paper NO. 4
!  VERSION NO. 1
!
!  LOGICAL COMPONENT NUMBER: P27
!
!  SYSTEM TASK : P27
!
!  DOCUMENTATION : Unified Model Documentation Paper P27
!                  Section (9)
!
!-----------------------------------------------------------------
!
!
!+ Convective Cloud Top, Base and Amount Calculation Scheme.
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection

      SUBROUTINE con_rad_4a5a (                                        &
                K,XPK,XPKP1,FLXKP1,BTERM,CCA,ICCB,ICCT,START_LEV,      &
                TCW,CCW,CCLWP,DELPKP1,LCCA,LCBASE,LCTOP,NPNTS,         &
                L_Q_INTERACT)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      
      USE earth_constants_mod, ONLY : g 
      
      IMPLICIT NONE



!
!----------------------------------------------------------------------
! VECTOR LENGTH AND LOOP VARIABLES
!----------------------------------------------------------------------
!
      INTEGER NPNTS        ! IN VECTOR LENGTH
!
      INTEGER K            ! IN PRESENT MODEL LAYER
!
      INTEGER I            ! LOOP COUNTER
!
! ----------------------------------------------------------------------
! Arguments with Intent IN.  ie:  Input variables.
! ----------------------------------------------------------------------
!
      REAL XPK(NPNTS)      ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
!
      REAL XPKP1(NPNTS)    ! IN PARCEL CLOUD WATER IN LAYER K+1 (KG/KG)
!
      LOGICAL BTERM(NPNTS) ! IN MASK FOR POINTS WHERE CONVECTION
                           !    IS ENDING
!
      REAL FLXKP1(NPNTS)   ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
!
      REAL DELPKP1(NPNTS)  ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
!
      REAL CCW(NPNTS)      ! IN PARCEL CLOUD WATER AS CALCULATED BEFORE
                           !    PRECIPITATION. LAYER K+1 (KG/KG)
!
!
      INTEGER START_LEV(NPNTS)  ! IN LEVEL AT WHICH CONVECTION INITIATED
!
      LOGICAL L_Q_INTERACT ! IN Switch (PC2) allows overwriting of
                           ! parcel variables
!
! ----------------------------------------------------------------------
! Arguments with Intent IN/OUT. ie: input variables changed on output.
! ----------------------------------------------------------------------
!
      REAL TCW(NPNTS)      ! INOUT
                           ! IN  TOTAL CONDENSED WATER SUMMED TO
                           !     LAYER K (KG/M**2/S)
                           ! OUT TOTAL CONDENSED WATER SUMMED TO
                           !     LAYER K+1 OR IF CONVECTION HAS
                           !     TERMINATED ZEROED (KG/M**2/S)
!
      REAL CCLWP(NPNTS)    ! INOUT
                           ! IN  TOTAL CLOUD LIQUID WATER PATH
                           !     SUMMED TO LAYER K  (KG/M**2)
                           ! OUT TOTAL CLOUD LIQUID WATER PATH
                           !     SUMMED TO LAYER K+1 (KG/M**2)
      REAL LCCA(NPNTS)      ! INOUT LOWEST CONV.CLOUD AMOUNT (%)
!
      INTEGER LCBASE(NPNTS) ! INOUT LOWEST CONV.CLOUD BASE LEVEL
!
      INTEGER LCTOP(NPNTS)  ! INOUT LOWEST CONV.CLOUD TOP LEVEL
!
! ----------------------------------------------------------------------
! Arguments with Intent OUT. ie: Output variables.
! ----------------------------------------------------------------------
!
      REAL CCA(NPNTS)      ! OUT CONVECTIVE CLOUD AMOUNT (%)
!
      INTEGER ICCB(NPNTS)   ! OUT CONVECTIVE CLOUD BASE LEVEL
!
      INTEGER ICCT(NPNTS)   ! OUT CONVECTIVE CLOUD TOP LEVEL
!


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!  External subroutine calls: ------------------------------------------
!
!     EXTERNAL None.
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
!
! ----------------------------------------------------------------------
!  CALCULATE Cloud Base and Lowest Cloud Base
!
!  WHEN CLOUD BASE SET ZERO TOTAL CONDENSED WATER
! ----------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('CON_RAD_4A5A',zhook_in,zhook_handle)
      DO I=1, NPNTS
        IF ( XPK(I)  <=  0.0 .AND. CCW(I)  >   0.0 ) THEN
!       Assuming initial parcel condensate is zero
          ICCB(I)=K+1
          CCLWP(I)=0.0
!
          IF ( LCBASE(I)  ==  0 ) THEN
!         Lowest cloud base
            LCBASE(I) = K+1
          END IF  ! If_lcbase_1
!
        ELSE IF ( L_Q_INTERACT .AND. XPK(I)  >   0.0 .AND.              &
     &           K  ==  START_LEV(I) ) THEN

          ! Do not know why this is done,
          ! but it lowers the ccb so that it is no longer consistent
          ! with ccw. Needs looking into.

          ! Non-zero initial parcel condensate (initialized to environment)
          ICCB(I)  = K
          CCLWP(I) = 0.0
!
          IF ( LCBASE(I)  ==  0 ) THEN
!         Lowest cloud base
            LCBASE(I) = K
          END IF  ! If_lcbase_2

        END IF
!
! ----------------------------------------------------------------------
!  CALCULATE Cloud Top and Lowest Cloud Top
! ----------------------------------------------------------------------
!
        IF (BTERM(I) .AND.                                              &
     &      ((CCW(I) >  0.0).OR.(XPK(I) >  0.0)) ) ICCT(I) = K+1

        IF (BTERM(I) .AND.  LCTOP(I) == 0 .AND.                         &
     &      ((CCW(I) >  0.0).OR.(XPK(I) >  0.0)) ) THEN
          LCTOP(I) = K+1
        END IF
!
        IF ( FLXKP1(I)  >   0.0) THEN
!
!---------------------------------------------------------------------
! SUM TOTAL CONDENSED WATER PER SECOND - ASSUMES THAT THE INITIAL
! CONVECTIVE LAYER IS UNSATURATED
!---------------------------------------------------------------------
!
          TCW(I) = TCW(I) + FLXKP1(I) * CCW(I) / G
!
!---------------------------------------------------------------------
! SUM CONV CONDENSED WATER PATH - ASSUMES THAT THE INITIAL
! CONVECTIVE LAYER IS UNSATURATED
!---------------------------------------------------------------------
!
          CCLWP(I) = CCLWP(I) + XPKP1(I) * DELPKP1(I) / G

        END IF
!
!---------------------------------------------------------------------
! CALCULATE CONVECTIVE CLOUD AMOUNT IF CONVECTION TERMINATES IN
! LAYER K AND TOTAL CONDENSED WATER PATH OVER A TIME STEP
!
! UM DOCUMENTATION PAPER P27
! SECTION (9), EQUATION (37)
!---------------------------------------------------------------------
!
        IF( BTERM(I) .AND. TCW(I) >  0.0 ) THEN
!
          IF ( TCW(I)  <   2.002E-6 ) TCW(I) = 2.002E-6
!
          CCA(I) = 0.7873 + 0.06 * LOG(TCW(I))
          IF (CCA(I)  >   1.0) CCA(I) = 1.0
!
          IF (LCCA(I) <= 0.0) THEN
            LCCA(I) = 0.7873 + 0.06 * LOG(TCW(I))
            IF (LCCA(I)  >   1.0) LCCA(I) = 1.0
          END IF
!
          TCW(I) = 0.0
!
        END IF
      END DO ! I loop over NPNTS
!
      IF (lhook) CALL dr_hook('CON_RAD_4A5A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE con_rad_4a5a
