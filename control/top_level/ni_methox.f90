! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! + Calculate chemical changes to water vapour due to methane oxidation
!   and hydrogen photolysis, using the method used at ECMWF
!   (Untch et al (ECMWF Newsletter No 87, winter 1998/99, pp 2-8) and
!   Simmons (pers. comm.)).
!
! Subroutine NI_methox. Called by atmos_physics1. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

      Subroutine NI_methox(                                             &
! Parallel variables
     &  halo_i, halo_j                                                  &

! model dimensions.
     &, row_length, rows                                                &
     &, model_levels, wet_model_levels                                  &

! model levels
     &, eta_theta_levels                                                &

! in time stepping information.
     &, timestep                                                        &

! in/out
     &,  q_n,q_inc                                                      &

! error information
     &, Error_code  )

      USE atmos_constants_mod, ONLY: pref
      USE conversions_mod, ONLY: pi

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE domain_params
      USE level_heights_mod, ONLY: z_top_of_model => z_top_theta

      IMPLICIT NONE

! purpose: Interface to Atmospheric Physics methane oxidation  code.
!
! method: Follows (Untch et al (ECMWF Newsletter No 87, winter
!   1998/99, pp 2-8) and Simmons (pers. comm.)). The model methane
!    mixing ratio is implicit, and derived from the assumption that
!    2 [CH4] + [H2O] = 3.75  ppmm throughout the stratosphere.
!   The methane oxidation and hydrogen photolysis rate coefficients
!   vary only with pressure, and are calculated
!   within the code assuming a surface pressure of 1000 hPa.
!
! code description:
!   language: fortran 77 + cray extensions
!
! arguments with intent in. ie: input variables.

! Parallel setup variables
      Integer                                                           &
     &  halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j     ! Size of halo in j direction.

! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, wet_model_levels
! model levels
      Real                                                              &
     &  eta_theta_levels(model_levels)
          
! The above looks incorrect as here dimensioned (1:model_levels)
! Calling routine it is dimensioned (0:model_levels)  
! Thus indexing is shifted by 1. To be investigated.    

! Model switches

! model parameters
      Real                                                              &
     &  timestep


! arguments with intent in/out. ie: input variables changed on output.
      Real                                                              &
     &  q_n(row_length, rows, wet_model_levels)                         &
     &, q_inc(row_length, rows, wet_model_levels)

! arguments with intent out. ie: output variables.

      Integer                                                           &
     &  Error_code

! local variables.

! loop counters
      Integer                                                           &
     &  i, j, k

! physical constants
      Real                                                              &
     &  max2MpW         ! maximum of 2CH4+H2O used to imply
                        ! methane amount
      PARAMETER(max2MpW=3.75E-06)
! methane oxidation (ak1) / hydrogen photolysis (ak1) coefficients.
      Real                                                              &
     &  ak1(wet_model_levels),ak2(wet_model_levels),                    &
     &  press(wet_model_levels)

!  calculation of ak1, ak2
      REAL                                                              &
     & ALP1,ALP2             ! used in calculation of methane oxidation
                             ! coefficients ak1 and ak2, respectively.
       REAL                                                             &
     & COEFF1_K1, COEFF2_K1                                             &
                             ! empirical coeffs used in calculation
     &,COEFF3_K1, COEFF4_K1                                             &
                             ! of ak1
     &,COEFF1_K2, COEFF2_K2                                             &
                             ! empirical coeffs used in calculation
     &,COEFF3_K2, COEFF4_K2                                             &
                             ! of ak2
     &,COEFF5_K2, COEFF6_K2                                             &
     &,P1_K1,P2_K1                                                      &
                             ! pressures for ak1 calculation
     &,P1_K2,P2_K2                                                      &
                             ! pressures for ak2 calculation
     &,SECS_IN_DAY                                                      &
                             ! No of seconds in day
     &,TAU1,TAU2                                                        &
                             ! Chemical timescales
     &,SCALE_HT          ! Scale height (in m)

       INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
       INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
       REAL(KIND=jprb)               :: zhook_handle
       PARAMETER(COEFF1_K1=19.0 ,COEFF2_K1=10.0 ,COEFF3_K1=20.0         &
     &    ,COEFF4_K1=100.0 ,COEFF1_K2=0.3333 ,COEFF2_K2=0.01            &
     &    ,COEFF3_K2=3.0 ,COEFF4_K2=100.0 ,COEFF5_K2=0.005              &
     &    ,COEFF6_K2=0.01 ,P1_K1=50.0 ,P2_K1=10000.0                    &
     &    ,P1_K2=0.1 ,P2_K2=20.0 ,SECS_IN_DAY=86400.00                  &
     &    ,SCALE_HT=7000.0)

      IF (lhook) CALL dr_hook('NI_METHOX',zhook_in,zhook_handle)
      If (error_code  ==  0 ) Then

! calculate ak1 and ak2 coefficients for methane oxidation.
! ----------------------------------------------------------
! Follows the approach of Simmons (pers. comm). Analytical
! functions are calculated that approximately match chemical
! lifetime plots shown in Brasseur and Solomon (1986) (the rate
! coefficient is the reciprocal of the chemical lifetime).
! It is assumed that the rate coefficients only vary with pressure.
          ALP1=COEFF1_K1*ALOG(COEFF2_K1)/(ALOG(COEFF3_K1))**4
          ALP2=ALOG(COEFF1_K2+COEFF2_K2)
!   calculate pressure
          DO J=1,WET_MODEL_LEVELS
             PRESS(J)=PREF*EXP(-(z_top_of_model*ETA_THETA_LEVELS(J)/SCALE_HT))
          END DO
!   calculate k1 and k2 coefficients
          DO J=1,WET_MODEL_LEVELS
             IF (PRESS(J) <= P1_K1) THEN
               AK1(J)=1/(SECS_IN_DAY*COEFF4_K1)
             ENDIF
             IF ((PRESS(J) >  P1_K1).AND.                               &
     &         (PRESS(J) <  P2_K1)) THEN
               TAU1=COEFF4_K1*(1+ALP1*                                  &
     &         ((ALOG(PRESS(J)/P1_K1))**4/                              &
     &         ALOG(P2_K1/PRESS(J))))
               AK1(j)=1/(SECS_IN_DAY*TAU1)
             ENDIF
             IF (PRESS(J) >= P2_K1) THEN
               AK1(j)=0.0
             ENDIF
             IF (PRESS(J) <= P1_K2) THEN
               AK2(J)=1/(SECS_IN_DAY*COEFF3_K2)
             ENDIF
             IF ((PRESS(J) >  P1_K2).AND.                               &
     &         (PRESS(J) <  P2_K2)) THEN
               TAU2=1/(EXP(ALP2-0.5*(ALOG(COEFF4_K2)+ALP2)*             &
     &         (1+COS(PI*ALOG(PRESS(J)/P2_K2)/                          &
     &         ALOG(COEFF5_K2))))-COEFF6_K2)
               AK2(J)=1/(SECS_IN_DAY*TAU2)
             ENDIF
             IF (PRESS(J) >= P2_K2) THEN
               AK2(j)=0.0
             ENDIF
          ENDDO

! ----------------------------------------------------------------------
! Add on q increments
! ----------------------------------------------------------------------

! add on increments to  q for next methox call
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                 If (q_n(i,j,k) >  0.0.and.q_n(i,j,k) <  max2MpW) Then
                   q_inc(i,j,k) = q_inc(i,j,k) + (ak1(k)*(max2MpW-      &
     &             q_n(i,j,k))-ak2(k)*q_n(i,j,k))*timestep
                 End If
              End Do
            End Do
          End Do

      End If ! on error code equal to zero

      IF (lhook) CALL dr_hook('NI_METHOX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_methox
