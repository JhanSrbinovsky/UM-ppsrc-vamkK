! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_FIT_FSAT-------------------------------------------
!
!   Purpose: To speed up the large scale hydrology code (LTOP=TRUE)
!            dramatically. This is done by fitting exponential
!            functions to the incomplete gamma function for each grid
!            box and the complete range of possible "water table"
!            (top_crit) cases - see documentation.
!            Estimates the fitted parameters for Fsat=function(ZW)
!            and  Fwet=function(ZW) for each land grid point.
!            (Calculating the incomplete gamma function for each grid
!            box at each time step was very time consuming).
!                                                             !
! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration

      SUBROUTINE CALC_FIT_FSAT(SOIL_PTS,SOIL_INDEX,NPNTS                &
     &  ,FEXP,TI_MEAN,TI_SIG,GAMTOT,ZDEPTH                              &
     &  ,A_FSAT,C_FSAT,A_FWET,C_FWET)
     
    Use PrintStatus_Mod, Only : &
        PrintStatus,                &
        PrStatus_Diag

     Use Ereport_Mod, Only :                                            &
       Ereport 

     IMPLICIT NONE

! C_TOPOG start
!
! Topographic index increment:
      REAL,PARAMETER :: DTI = 0.2
! Maximum topographic index considered:
      REAL,PARAMETER :: TI_MAX = 10.0
! Maximum allowed water table depth (m):
      REAL,PARAMETER :: ZW_MAX = 6.0
! Standard deviation of LOG(Ksat(0)):
      REAL,PARAMETER :: SIGMA_LOGK = 0.0
! Parameter to remove very high water tables
! from the calculated wetland fraction:
      REAL,PARAMETER :: TI_WETL = 1.5
!
! C_TOPOG end

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS                                                            &
                        ! IN No. of land points.
     &,SOIL_PTS         ! IN No. of land soil points.

      REAL                                                              &
     & ZDEPTH           ! IN Standard Soil model DEPTH.

!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)! IN Array of soil points.

      REAL                                                              &
     & TI_MEAN(NPNTS)                                                   &
                        ! IN Gridbox mean topographic index.
     &,TI_SIG(NPNTS)                                                    &
                        ! IN Std. deviation in topographic index.
     &,FEXP(NPNTS)                                                      &
                        ! IN Exp. decay in deep layer.
     &,GAMTOT(NPNTS)    ! IN Total gamma function.

!   Array arguments with intent(OUT) :
      REAL                                                              &
     &   A_FSAT(NPNTS)                                                  &
                        ! OUT Fitting parameter for Fsat.
     &  ,C_FSAT(NPNTS)                                                  &
                        ! OUT Fitting parameter for Fsat.
     &  ,A_FWET(NPNTS)                                                  &
                        ! OUT Fitting parameter for Fwet.
     &  ,C_FWET(NPNTS)  ! OUT Fitting parameter for Fwet.

! Local scalars:
      INTEGER NZW  ! Number of ZW values used in fit.
      PARAMETER(NZW=200)  ! Maximum value for a significant improvement
!                         ! in the fit.

      INTEGER                                                           &
     & I,J,IZ                                                           &
                   ! Loop counters.
     &,IFITA                                                            &
                   ! Loop counters.
     &,NFITA       ! Number of loops for fitting.
      PARAMETER(NFITA=20)

      REAL DZW     ! WORK ZW increment ; defined by ZW_MAX and NZW.

      REAL                                                              &
     & RMS                                                              &
                   ! WORK RMS errors for given fsat fit values.
     &,RMSW                                                             &
                   ! WORK RMS errors for given fwet fit values.
     &,RMSOLD                                                           &
                   ! WORK RMS errors for given fsat fit values.
!                  !      for best fit so far.
     &,RMSWOLD                                                          &
                   ! WORK RMS errors for given fwet fit values
!                  !      for best fit so far.
     &,CFITMIN                                                          &
                   ! WORK Minimum possible value for Cfit.
     &,CFITMAX                                                          &
                   ! WORK Maximum possible value for Cfit.
     &,CFIT                                                             &
                   ! WORK CFit value for given loop.
     &,THR_ERR     ! WORK Error threshold value


      PARAMETER(CFITMIN=0.0,CFITMAX=3.00)
      PARAMETER(THR_ERR=5.0E-3)

! Local arrays:
      REAL                                                              &
     & FSAT_CALC(NZW)                                                   &
                              ! WORK Surface saturation fraction.
     &,FSAT_FIT(NZW)                                                    &
                              ! WORK Fitted surface saturation fraction.
     &,FWET_CALC(NZW)                                                   &
                              ! WORK Wetland fraction.
     &,FWET_FIT(NZW)                                                    &
                              ! WORK Fitted wetland fraction.
     &,DUMZW(NZW)                                                       &
                              ! WORK Dummy water table depth (m).
     &,DUMFSAT                                                          &
                              ! WORK Dummy surface saturation fraction.
     &,DUMFWETL               ! WORK Dummy wetland fraction.

      REAL                                                              &
     & TOP_CRIT(NZW)                                                    &
                              ! WORK LOG(QBASE_MAX/QBASE) -see document.
     &,TOP_CRIT1Z                                                       &
                              ! WORK Aas above but for an individual zw.
     &,TOP_MIN                                                          &
                              ! WORK value for when zw=zw_max.
     &,WUTOT                  ! WORK Dummy (set to 1.0).

      REAL temp1
      INTEGER ERRORSTATUS
      CHARACTER (LEN=80)  CMESSAGE
      CHARACTER (LEN=13)  ROUTINENAME

      ROUTINENAME='CALC_FIT_FSAT'

      ! Define the water table depths to be used in the fitting process:
      DZW=1.0/REAL(NZW)*ZW_MAX
      DO IZ=1,NZW
        DUMZW(IZ)=REAL(IZ-1)*DZW
      ENDDO
      WUTOT=1.0
      DUMFSAT=0.0
      DUMFWETL=0.0

! Calculate TOP_CRIT for the water tables depths:
      DO J=1,SOIL_PTS
         I=SOIL_INDEX(J)

         IF(TI_MEAN(I) >  0.0 .AND. TI_SIG(I) >  0.0 )THEN
            TOP_MIN=1.0/FEXP(I)*EXP(-FEXP(I)*(ZW_MAX-ZDEPTH))

            temp1 = ZDEPTH+1.0/FEXP(I)-TOP_MIN

            DO IZ=1,NZW

               IF(DUMZW(IZ) <= ZDEPTH)TOP_CRIT1Z=                    &
                    -LOG(1-(DUMZW(IZ)/temp1))
               IF(DUMZW(IZ) >  ZDEPTH)TOP_CRIT1Z=                    &
     &               LOG(temp1/                                      &
     &                   (1./FEXP(I)*EXP(-FEXP(I)*(DUMZW(IZ)-ZDEPTH))&
     &                    -TOP_MIN))

! Calculate FSAT and FWET for one ZW at one soil pnt:
! DEPENDS ON: calc_fsat
               CALL CALC_FSAT(.FALSE.,1,1,1,TI_MEAN(I),TI_SIG(I),&
                              WUTOT,TOP_CRIT1Z,GAMTOT(I),DUMFSAT,DUMFWETL)

               FSAT_CALC(IZ)=DUMFSAT
               FWET_CALC(IZ)=DUMFWETL
               TOP_CRIT(IZ)=TOP_CRIT1Z
               IF(IZ == 1)THEN   ! Values at zw=0m
                  A_FSAT(I)=FSAT_CALC(IZ)
                  A_FWET(I)=FWET_CALC(IZ)
               ENDIF
            ENDDO

            RMSOLD=1.0E10
            RMSWOLD=1.0E10

            DO IFITA=1,NFITA
               CFIT=CFITMAX*(IFITA)/FLOAT(NFITA)

! This isnt really root mean squared - just a measure of fitness.
               RMS=0.0
               RMSW=0.0
!TOP_CRIT=TI_MAX when zw=zw_max
               DO IZ=1,NZW
                  FSAT_FIT(IZ)=A_FSAT(I)*EXP(-CFIT*TOP_CRIT(IZ))
                  FWET_FIT(IZ)=A_FWET(I)*EXP(-CFIT*TOP_CRIT(IZ))
                  RMS=RMS+(FSAT_CALC(IZ)-FSAT_FIT(IZ))**2
                  RMSW=RMSW+(FWET_CALC(IZ)-FWET_FIT(IZ))**2
               ENDDO            !ZW
               RMS=RMS/FLOAT(NZW)
               RMSW=RMSW/FLOAT(NZW)

               IF(RMS < RMSOLD)THEN
                  RMSOLD=RMS
                  C_FSAT(I)=CFIT
               ENDIF
               IF(RMSW <  RMSWOLD)THEN
                  RMSWOLD=RMSW
                  C_FWET(I)=CFIT
               ENDIF
            ENDDO

            IF(RMSOLD >= THR_ERR**2)THEN
             IF(C_FSAT(I) <= CFITMIN.OR.C_FSAT(I) >= CFITMAX)THEN
               
               DO IZ=1,NZW
                  FSAT_FIT(IZ)=A_FSAT(I)*EXP(-C_FSAT(I)*TOP_CRIT(IZ))
               ENDDO               !ZW

               Write(6,'(a,i4,3f20.10)')                                &
                 'ERROR CFIT FSAT',I,C_FSAT(I),CFITMIN,CFITMAX
               Write(6,'(a,3f20.10)')'fsat_calc='                       &
     &           ,FSAT_CALC(1),FSAT_CALC(3),FSAT_CALC(5)
               Write(6,'(a,3f20.10)')'fsat_fit='                        &
     &           ,FSAT_FIT(1),FSAT_FIT(3),FSAT_FIT(5)
               Write(6,'(a,f20.10)')'RMS=',SQRT(RMSOLD)
               ErrorStatus = 35
               Write(CMessage, '(A)') 'Error in CFIT FSAT in LSH model setup'

                 Call Ereport ( RoutineName, ErrorStatus, CMessage)
             ENDIF
            ENDIF

            IF(RMSWOLD >= THR_ERR**2)THEN
             IF(C_FWET(I) <= CFITMIN.OR.C_FWET(I) >= CFITMAX)THEN
 
               DO IZ=1,NZW
                  FSAT_FIT(IZ)=A_FSAT(I)*EXP(-C_FSAT(I)*TOP_CRIT(IZ))
                  FWET_FIT(IZ)=A_FWET(I)*EXP(-C_FWET(I)*TOP_CRIT(IZ))
               ENDDO               !ZW

               Write(6,'(a,i4,3f20.10)')                                &
                 'ERROR CFIT FWET',I,C_FWET(I),CFITMIN,CFITMAX
               Write(6,'(a,3f20.10)')'fwet_calc='                       &
     &          ,FWET_CALC(1),FWET_CALC(3),FWET_CALC(5)
               Write(6,'(a,3f20.10)')'fwet_fit='                        &
     &          ,FWET_FIT(1),FWET_FIT(3),FWET_FIT(5)
               Write(6,'(a,f20.10)')'RMSW=',SQRT(RMSWOLD)
               Write(6,'(a,3f20.10)')'(fsat_calc=)'                     &
     &          ,FSAT_CALC(1),FSAT_CALC(3),FSAT_CALC(5)
               Write(6,'(a,3f20.10)')'(fsat_fit=)'                      &
     &          ,FSAT_FIT(1),FSAT_FIT(3),FSAT_FIT(5)
               Write(6,'(a,f20.10)')'(RMS=)',SQRT(RMSOLD)
               ErrorStatus = 40
               Write(CMessage, '(A)') 'Error in CFIT FWET in LSH model setup'

               Call Ereport ( RoutineName, ErrorStatus, CMessage)
             ENDIF
            ENDIF
            IF(     (RMSOLD >= THR_ERR**2 .OR. RMSWOLD >= THR_ERR**2)   &
               .AND.(PrintStatus >= PrStatus_Diag)) THEN
              Write(6,'(a,2f20.10)')                                    &
                'Warning LSH RMS Error in fit:'                         &
               ,SQRT(RMSOLD),SQRT(RMSWOLD)
            END IF
         ENDIF

      ENDDO                     ! NPNTS

      END SUBROUTINE CALC_FIT_FSAT
