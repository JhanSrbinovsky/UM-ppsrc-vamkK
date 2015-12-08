! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Perform a chemical integration using the Backward Euler
!  solver for the UKCA model with RAQ chemistry (based on STOCHEM).
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method: Chemistry integration (Backward Euler) for RAQ chemistry
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
! 
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      SUBROUTINE UKCA_DERIV_RAQ(nr_therm, n_be_calls,                  &
                          n_pnts,                                      &
                          rc, dw, dd, dj,                              &
                          h2o, m, o2,                                  &
                          dts, y)

      USE ukca_option_mod, ONLY: nit, jppj, jpspec
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Control_Max_Sizes
      USE VECTLIB_MOD, ONLY :                                          &
      ONEOVER_V
      IMPLICIT NONE
      
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_S,                                                       &
                                 ! relative molecular mass S kg/mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_HO2,                                                     &
                                 ! relative molecular mass HO2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_S    = 3.20E-2,                                    &
                                          ! kg/mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &           RMM_HO2  = 3.30E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------
      
      INTEGER, INTENT(IN) :: nr_therm   ! No. thermal reactions
      INTEGER, INTENT(IN) :: n_be_calls ! No. chemical steps
      INTEGER, INTENT(IN) :: n_pnts     ! Actual no. calculations

      REAL, INTENT(IN)    :: dts                 ! timestep
      REAL, INTENT(IN)    :: rc(n_pnts,nr_therm) ! rxn rate coeffs
      REAL, INTENT(IN)    :: dj(n_pnts,jppj)     ! photol rates
      REAL, INTENT(IN)    :: dw(n_pnts,jpspec)   ! wet dep rates
      REAL, INTENT(IN)    :: dd(n_pnts,jpspec)   ! dry dep rates
      REAL, INTENT(IN)    :: h2o(n_pnts)         ! h2o concn
      REAL, INTENT(IN)    :: m(n_pnts)           ! air density
      REAL, INTENT(IN)    :: o2(n_pnts)          ! o2 density

      REAL, INTENT(INOUT) :: y(n_pnts,jpspec)    ! species concn

!     Local variables 

      INTEGER :: i, j, n                       ! loop variables
            
      REAL :: p(n_pnts)                        ! production rate
      REAL :: r1(n_pnts)                       ! working array
      REAL :: r2(n_pnts)                       ! working array
      REAL :: l(n_pnts)                        ! loss rate
      REAL :: l1(n_pnts)                       ! working array
      REAL :: l2(n_pnts)                       ! working array
      REAL :: l3(n_pnts)                       ! working array
      REAL :: yp(n_pnts,jpspec)                ! previous species concn
      REAL :: tmp1(n_pnts)                     ! work array for vector reciprocal
      REAL :: tmp2(n_pnts)                     ! work array for vector reciprocal  


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_DERIV_RAQ',zhook_in,zhook_handle)

!     Initialise flux diagnostics to zero
      tmp1 = 0.0
      tmp2 = 0.0

      DO n = 1, n_be_calls           ! loop over chemical timesteps
        
        DO j = 1,jpspec
          yp(:,j) = y(:,j)
        END DO

!       Iteration start

        DO i=1,nit                   ! loop over iterations
!
! This section written automatically by MECH9GEN
! from the file newmech9.txt
! with 178 equations and 59 species (including h2o).
!
! Note: Loss terms of O3 [y(4)] have been copied manually for O3S [y(28)].
! After that (DD(:,4)  ) has to be changed to (DD(:,28)  ) there. 
!
! Note: After the automatic generation, the code has been re-written to
! split additions for chemical production and loss terms to six per line 
! maximum. 
! This significantly increases the speed of the additions on the IBM
!
!
!          O(3P)        y( 1)

      p = 0.0                                                           &
      +(DJ(:,3)  *y(:,7 ))        +(DJ(:,15) *y(:,6 ))                  &
      +(RC(:,7)  *y(:,2 ))        +(DJ(:,1)  *y(:,4 ))                  
      l = 0.0                                                           &
      +(RC(:,1)  )        +(RC(:,5)  *y(:,5 ))+(RC(:,16) *y(:,7 )) 
         
      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)
 
      y(:, 1) = (yp(:, 1)+dts*p)*tmp1

!
!          O(1D)        y( 2)
      p = 0.0                                                           &
      +(DJ(:,2)  *y(:,4 ))                                              
      l = 0.0                                                           &
      +(RC(:,7)  )        +(RC(:,8)  *h2o)     

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:, 2) = (yp(:, 2)+dts*p)*tmp1
!
!          OH           y( 3)
      p = 0.0                                                           &
      +(DJ(:,23) *y(:,35))                                              &
      +(DJ(:,21) *y(:,52))        +(DJ(:,22) *y(:,31))                  &
      +(DJ(:,19) *y(:,21))        +(DJ(:,20) *y(:,26))                  &
      +(DJ(:,5)  *y(:,10))
      p = p                                                             & 
      +(DJ(:,17) *y(:,17))                                              &
      +(DJ(:,4)  *y(:,11))        +(DJ(:,4)  *y(:,11))                  &
      +(RC(:,137)*y(:,35)*y(:,3 ))+(RC(:,182)*y(:,15)*y(:,55))          &
      +(RC(:,129)*y(:,4 )*y(:,34)*0.36)
      p = p                                                             &
      +(RC(:,135)*y(:,31)*y(:,3 ))                                      &
      +(RC(:,123)*y(:,4 )*y(:,50)*0.28)+(RC(:,128)*y(:,4 )*y(:,29)*0.27)&
      +(RC(:,102)*y(:,3 )*y(:,26))+(RC(:,104)*y(:,3 )*y(:,52))          &
      +(RC(:,45) *y(:,3 )*y(:,17))
      p = p                                                             &
      +(RC(:,100)*y(:,3 )*y(:,21))                                      &
      +(RC(:,17) *y(:,5 )*y(:,15))+(RC(:,34) *y(:,6 )*y(:,15))          &
      +(RC(:,8)  *y(:,2 )*h2o*2.00)    +(RC(:,14) *y(:,15)*y(:,4 ))     
      l = 0.0                                                           &
      +(RC(:,177)*y(:,39))                                              &
      +(RC(:,174)*y(:,54))+(RC(:,175)*y(:,54))+(RC(:,176)*y(:,58))      &
      +(RC(:,160)*y(:,32))+(RC(:,170)*y(:,58))
      l = l                                                             &
      +(RC(:,172)*y(:,56))                                              &
      +(RC(:,143)*y(:,34))+(RC(:,154)*y(:,45))+(RC(:,156)*y(:,47))      &
      +(RC(:,138)*y(:,33))+(RC(:,139)*y(:,57))
      l = l                                                             &
      +(RC(:,141)*y(:,29))                                              &
      +(RC(:,125)*y(:,50))+(RC(:,135)*y(:,31))+(RC(:,137)*y(:,35))      &
      +(RC(:,102)*y(:,26))+(RC(:,104)*y(:,52))
      l = l                                                             &
      +(RC(:,109)*y(:,49))                                              &
      +(RC(:,94) *y(:,27))+(RC(:,98) *y(:,23))+(RC(:,100)*y(:,21))      &
      +(RC(:,81) *y(:,51))+(RC(:,86) *y(:,53))
      l = l                                                             &
      +(RC(:,92) *y(:,25))                                              &
      +(RC(:,70) *y(:,13))+(RC(:,71) *y(:,19))+(RC(:,75) *y(:,22))      &
      +(RC(:,59) *y(:,12))+(RC(:,63) *y(:,41))
      l = l                                                             &
      +(RC(:,66) *y(:,14))                                              &
      +(RC(:,44) *y(:,21))+(RC(:,44) *y(:,26))+(RC(:,44) *y(:,52))      &
      +(RC(:,35) *y(:,10))+(RC(:,44) *y(:,17))
      l = l                                                             &
      +(RC(:,45) *y(:,17))                                              &
      +(RC(:,30) *y(:,15))+(RC(:,31) *y(:,11))
      l = l                                                             &
      +(RC(:,33) *y(:,43))                                              &
      +(RC(:,13) *y(:,4 ))+(RC(:,21) *y(:,7 ))+(RC(:,24) *y(:,9 ))  

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:, 3) = (yp(:, 3)+dts*p)*tmp1
!
!          O3           y( 4)
      p = 0.0                                                           &
      +(RC(:,1)  *y(:,1 ))        +(RC(:,131)*y(:,20)*y(:,15)*0.3 )     
      l = 0.0                                                           &
      +(DJ(:,1)  )        +(DJ(:,2)  )        +(DD(:,4)  )              &
      +(RC(:,124)*y(:,50))+(RC(:,128)*y(:,29))+(RC(:,129)*y(:,34))
      l = l                                                             &
      +(RC(:,14) *y(:,15))+(RC(:,112)*y(:,49))+(RC(:,123)*y(:,50))      &
      +(RC(:,11) *y(:,5 ))+(RC(:,12) *y(:,7 ))+(RC(:,13) *y(:,3 )) 

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:, 4) = (yp(:, 4)+dts*p)*tmp1
!
!          NO           y( 5)
      p = 0.0                                                           &
      +(DJ(:,3)  *y(:,7 ))        +(DJ(:,14) *y(:,6 ))                  &
      +(RC(:,16) *y(:,7 )*y(:,1 ))+(RC(:,19) *y(:,7 )*y(:,6 ))          
      l = 0.0                                                           &
      +(RC(:,173)*y(:,44))                                              &
      +(RC(:,126)*y(:,40))+(RC(:,142)*y(:,46))+(RC(:,144)*y(:,48))      &
      +(RC(:,95) *y(:,36))+(RC(:,105)*y(:,37))
      l = l                                                             &
      +(RC(:,110)*y(:,38))                                              &
      +(RC(:,79) *y(:,20))+(RC(:,83) *y(:,24))+(RC(:,93) *y(:,30))      &
      +(RC(:,17) *y(:,15))
      l = l                                                             &
      +(RC(:,60) *y(:,16))+(RC(:,72) *y(:,18))                          &
      +(RC(:,5)  *y(:,1 ))+(RC(:,11) *y(:,4 ))+(RC(:,15) *y(:,6 ))   
   
      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:, 5) = (yp(:, 5)+dts*p)*tmp1
!
!          NO3          y( 6)
      p = 0.0                                                           &
      +(DJ(:,16) *y(:,8 ))                                              &
      +(RC(:,35) *y(:,3 )*y(:,10))+(RC(:,98) *y(:,3 )*y(:,23))          &
      +(RC(:,12) *y(:,7 )*y(:,4 ))+(RC(:,29) *y(:,8 ))                  
      l = 0.0                                                           &
      +(DJ(:,14) )        +(DJ(:,15) )        +(DW(:,6)  )              &
      +(RC(:,158)*y(:,22))+(RC(:,159)*y(:,29))+(RC(:,178)*y(:,39))
      l = l                                                             &
      +(RC(:,152)*y(:,51))+(RC(:,153)*y(:,49))+(RC(:,155)*y(:,50))      &
      +(RC(:,34) *y(:,15))+(RC(:,67) *y(:,14))+(RC(:,151)*y(:,19))
      l = l                                                             &
      +(RC(:,27) *y(:,6 ))+(RC(:,27) *y(:,6 ))+(RC(:,32) *y(:,15))      &
      +(RC(:,15) *y(:,5 ))+(RC(:,19) *y(:,7 ))+(RC(:,20) *y(:,7 ))   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:, 6) = (yp(:, 6)+dts*p)*tmp1
!
!          NO2          y( 7)
      p = 0.0                                                           &
      +(DJ(:,18) *y(:,23))                                              &
      +(DJ(:,15) *y(:,6 ))        +(DJ(:,16) *y(:,8 ))                  &
      +(DJ(:,5)  *y(:,10))        +(DJ(:,11) *y(:,9 ))                  &
      +(RC(:,177)*y(:,3 )*y(:,39))
      p = p                                                             &
      +(RC(:,178)*y(:,6 )*y(:,39)*2.00)                                 &
      +(RC(:,160)*y(:,32)*y(:,3 ))+(RC(:,173)*y(:,44)*y(:,5 ))          &
      +(RC(:,154)*y(:,45)*y(:,3 ))+(RC(:,156)*y(:,47)*y(:,3 ))          &
      +(RC(:,142)*y(:,46)*y(:,5 ))
      p = p                                                             &
      +(RC(:,144)*y(:,48)*y(:,5 ))                                      &
      +(RC(:,110)*y(:,38)*y(:,5 ))+(RC(:,126)*y(:,40)*y(:,5 ))          &
      +(RC(:,95) *y(:,36)*y(:,5 ))+(RC(:,105)*y(:,5 )*y(:,37))          &
      +(RC(:,83) *y(:,24)*y(:,5 ))
      p = p                                                             &
      +(RC(:,93) *y(:,30)*y(:,5 ))                                      &
      +(RC(:,78) *y(:,23))        +(RC(:,79) *y(:,20)*y(:,5 ))          &
      +(RC(:,60) *y(:,5 )*y(:,16))+(RC(:,72) *y(:,18)*y(:,5 ))          &
      +(RC(:,29) *y(:,8 ))        
      p = p                                                             &
      +(RC(:,34) *y(:,6 )*y(:,15))                                      &
      +(RC(:,24) *y(:,9 )*y(:,3 ))+(RC(:,27) *y(:,6 )*y(:,6 )*2.00)     &
      +(RC(:,19) *y(:,7 )*y(:,6 ))+(RC(:,23) *y(:,9 ))
      p = p                                                             &
      +(RC(:,15) *y(:,5 )*y(:,6 )*2.00)+(RC(:,17) *y(:,5 )*y(:,15))     &
      +(RC(:,5)  *y(:,1 )*y(:,5 ))+(RC(:,11) *y(:,5 )*y(:,4 ))          
      l = 0.0                                                           &
      +(DJ(:,3)  )        +(DD(:,7)  )                                  &
      +(RC(:,77) *y(:,20))+(RC(:,171)*y(:,42))+(RC(:,183)*y(:,55))
      l = l                                                             &
      +(RC(:,20) *y(:,6 ))+(RC(:,21) *y(:,3 ))+(RC(:,22) *y(:,15))      &
      +(RC(:,12) *y(:,4 ))+(RC(:,16) *y(:,1 ))+(RC(:,19) *y(:,6 ))
      
      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:, 7) = (yp(:, 7)+dts*p)*tmp1
!
!          N2O5         y( 8)
      p = 0.0                                                           &
      +(RC(:,20) *y(:,7 )*y(:,6 ))                                      
      l = 0.0                                                           &
      +(RC(:,29) )        +(DJ(:,16) )        +(DW(:,8)  )     

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:, 8) = (yp(:, 8)+dts*p)*tmp1
!
!          HO2NO2       y( 9)
      p = 0.0                                                           &
      +(RC(:,22) *y(:,7 )*y(:,15))                                      
      l = 0.0                                                           &
      +(DW(:,9)  )                                                      &
      +(RC(:,23) )        +(RC(:,24) *y(:,3 ))+(DJ(:,11) )   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:, 9) = (yp(:, 9)+dts*p)*tmp1
!
!          HONO2        y(10)
      p = 0.0                                                           &
      +(RC(:,152)*y(:,6 )*y(:,51))+(RC(:,158)*y(:,6 )*y(:,22))          &
      +(RC(:,67) *y(:,6 )*y(:,14))+(RC(:,151)*y(:,6 )*y(:,19))          &
      +(RC(:,21) *y(:,7 )*y(:,3 ))+(RC(:,32) *y(:,6 )*y(:,15))          
      l = 0.0                                                           &
      +(DD(:,10) )                                                      &
      +(RC(:,35) *y(:,3 ))+(DJ(:,5)  )        +(DW(:,10) ) 

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,10) = (yp(:,10)+dts*p)*tmp1
!
!          H2O2         y(11)
      p = 0.0                                                           &
      +(RC(:,36) *y(:,15)*y(:,15))                                      
      l = 0.0                                                           &
      +(DD(:,11) )                                                      &
      +(RC(:,31) *y(:,3 ))+(DJ(:,4)  )        +(DW(:,11) )   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,11) = (yp(:,11)+dts*p)*tmp1
!
!          CH4          y(12)
      p = 0.0                                                           &
      +(RC(:,123)*y(:,4 )*y(:,50)*0.30)                                 
      l = 0.0                                                           &
      +(RC(:,59) *y(:,3 ))+(DD(:,12) )   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,12) = (yp(:,12)+dts*p)*tmp1
!
!          CO           y(13)
      p = 0.0                                                           &
      +(DJ(:,13) *y(:,57))                                              &
      +(DJ(:,12) *y(:,33))        +(DJ(:,13) *y(:,57))                  &
      +(DJ(:,7)  *y(:,14))        +(DJ(:,8)  *y(:,22))                  &
      +(RC(:,139)*y(:,57)*y(:,3 )*2.0 )
      p = p                                                             &
      +(DJ(:,6)  *y(:,14))                                              &
      +(RC(:,129)*y(:,4 )*y(:,34)*0.76)+(RC(:,138)*y(:,33)*y(:,3 ))     &
      +(RC(:,124)*y(:,4 )*y(:,50)*0.58)+(RC(:,128)*y(:,4 )*y(:,29)*0.78)
      p = p                                                             &
      +(RC(:,112)*y(:,4 )*y(:,49)*0.31)+(RC(:,123)*y(:,4 )*y(:,50)*0.40)&
      +(RC(:,66) *y(:,3 )*y(:,14))+(RC(:,67) *y(:,6 )*y(:,14))          
      l = 0.0                                                           &
      +(RC(:,70) *y(:,3 ))+(DD(:,13) )            

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,13) = (yp(:,13)+dts*p)*tmp1
!
!          HCHO         y(14)
      p = 0.0                                                           &
      +(DJ(:,22) *y(:,31))        +(DJ(:,23) *y(:,35))                  &
      +(RC(:,181)*y(:,44)*y(:,16))+(DJ(:,17) *y(:,17))                  &
      +(RC(:,144)*y(:,48)*y(:,5 ))+(RC(:,154)*y(:,45)*y(:,3 ))
      p = p                                                             &
      +(RC(:,137)*y(:,35)*y(:,3 ))+(RC(:,142)*y(:,46)*y(:,5 ))          &
      +(RC(:,133)*y(:,48)*y(:,16)*1.0 )+(RC(:,135)*y(:,31)*y(:,3 ))     &
      +(RC(:,129)*y(:,4 )*y(:,34)*0.24)+(RC(:,132)*y(:,46)*y(:,16)*1.0 )
      p = p                                                             &
      +(RC(:,127)*y(:,16)*y(:,40)*2.00)+(RC(:,128)*y(:,4 )*y(:,29)*0.22)&
      +(RC(:,123)*y(:,4 )*y(:,50))+(RC(:,126)*y(:,40)*y(:,5 ))          &
      +(RC(:,112)*y(:,4 )*y(:,49))+(RC(:,112)*y(:,4 )*y(:,49)*0.47)
      p = p                                                             &
      +(RC(:,110)*y(:,38)*y(:,5 )*2.00)+(RC(:,111)*y(:,16)*y(:,38)*3.00)&
      +(RC(:,98) *y(:,3 )*y(:,23))+(RC(:,106)*y(:,16)*y(:,37))          &
      +(RC(:,96) *y(:,36)*y(:,16)*2.00)+(RC(:,97) *y(:,30)*y(:,16))
      p = p                                                             &
      +(RC(:,84) *y(:,24)*y(:,16))+(RC(:,95) *y(:,36)*y(:,5 ))          &
      +(RC(:,74) *y(:,16)*y(:,20)*2.00)+(RC(:,80) *y(:,16)*y(:,20))     &
      +(RC(:,63) *y(:,41)*y(:,3 ))
      p = p                                                             &
      +(RC(:,73) *y(:,18)*y(:,16))                                      &
      +(RC(:,61) *y(:,16)*y(:,16)*2.00)+(RC(:,62) *y(:,16)*y(:,16))     &
      +(RC(:,45) *y(:,3 )*y(:,17))+(RC(:,60) *y(:,5 )*y(:,16))          
      l = 0.0                                                           &
      +(DJ(:,7)  )        +(DW(:,14) )                                  &
      +(RC(:,66) *y(:,3 ))+(RC(:,67) *y(:,6 ))+(DJ(:,6)  )  

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,14) = (yp(:,14)+dts*p)*tmp1
!
!          HO2          y(15)
      p = 0.0                                                           &
      +(DJ(:,22) *y(:,31))        +(DJ(:,23) *y(:,35))                  &
      +(DJ(:,20) *y(:,26))        +(DJ(:,21) *y(:,52))                  &
      +(DJ(:,17) *y(:,17))        +(DJ(:,19) *y(:,21))
      p = p                                                             &
      +(DJ(:,13) *y(:,57))        +(DJ(:,13) *y(:,57))                  &
      +(DJ(:,11) *y(:,9 ))        +(DJ(:,12) *y(:,33))                  &
      +(DJ(:,6)  *y(:,14))        +(DJ(:,8)  *y(:,22))
      p = p                                                             &
      +(RC(:,181)*y(:,44)*y(:,16)*2.00)+(DJ(:,6)  *y(:,14))             &
      +(RC(:,173)*y(:,44)*y(:,5 ))+(RC(:,174)*y(:,3 )*y(:,54))          &
      +(RC(:,159)*y(:,6 )*y(:,29))+(RC(:,170)*y(:,3 )*y(:,58))
      p = p                                                             &
      +(RC(:,153)*y(:,6 )*y(:,49))+(RC(:,155)*y(:,6 )*y(:,50))          &
      +(RC(:,142)*y(:,46)*y(:,5 ))+(RC(:,144)*y(:,48)*y(:,5 ))          &
      +(RC(:,133)*y(:,48)*y(:,16)*2.0 )+(RC(:,139)*y(:,57)*y(:,3 )*1.0 )
      p = p                                                             &
      +(RC(:,129)*y(:,4 )*y(:,34)*0.36)+(RC(:,132)*y(:,46)*y(:,16)*2.0 )&
      +(RC(:,127)*y(:,16)*y(:,40)*2.00)+(RC(:,128)*y(:,4 )*y(:,29)*0.27)&
      +(RC(:,124)*y(:,4 )*y(:,50)*0.18)+(RC(:,126)*y(:,40)*y(:,5 ))
      p = p                                                             &
      +(RC(:,112)*y(:,4 )*y(:,49)*0.20)+(RC(:,123)*y(:,4 )*y(:,50)*0.30)&
      +(RC(:,110)*y(:,38)*y(:,5 ))+(RC(:,111)*y(:,16)*y(:,38)*2.00)     &
      +(RC(:,97) *y(:,30)*y(:,16)*2.00)+(RC(:,106)*y(:,16)*y(:,37))
      p = p                                                             &
      +(RC(:,93) *y(:,30)*y(:,5 ))+(RC(:,96) *y(:,36)*y(:,16))          &
      +(RC(:,84) *y(:,24)*y(:,16)*2.00)+(RC(:,90) *y(:,18)*y(:,18)*2.00)&
      +(RC(:,80) *y(:,16)*y(:,20))+(RC(:,83) *y(:,24)*y(:,5 ))
      p = p                                                             &
      +(RC(:,72) *y(:,18)*y(:,5 ))+(RC(:,73) *y(:,18)*y(:,16)*2.00)     &
      +(RC(:,67) *y(:,6 )*y(:,14))+(RC(:,70) *y(:,3 )*y(:,13))          &
      +(RC(:,63) *y(:,41)*y(:,3 ))+(RC(:,66) *y(:,3 )*y(:,14))
      p = p                                                             &
      +(RC(:,60) *y(:,5 )*y(:,16))+(RC(:,61) *y(:,16)*y(:,16)*2.00)     &
      +(RC(:,31) *y(:,3 )*y(:,11))+(RC(:,33) *y(:,3 )*y(:,43))          &
      +(RC(:,13) *y(:,3 )*y(:,4 ))+(RC(:,23) *y(:,9 ))                  
      l = 0.0                                                           &
      +(RC(:,180)*y(:,42))+(RC(:,182)*y(:,55))+(DW(:,16) )              &
      +(RC(:,131)*y(:,20))+(RC(:,134)*y(:,46))+(RC(:,136)*y(:,48))
      l = l                                                             &
      +(RC(:,99) *y(:,18))+(RC(:,101)*y(:,30))+(RC(:,103)*y(:,24))      &
      +(RC(:,36) *y(:,15))+(RC(:,36) *y(:,15))+(RC(:,65) *y(:,16))
      l = l                                                             &
      +(RC(:,30) *y(:,3 ))+(RC(:,32) *y(:,6 ))+(RC(:,34) *y(:,6 ))      &
      +(RC(:,14) *y(:,4 ))+(RC(:,17) *y(:,5 ))+(RC(:,22) *y(:,7 ))

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,15) = (yp(:,15)+dts*p)*tmp1
!
!          MeOO         y(16)
      p = 0.0                                                           &
      +(DJ(:,10) *y(:,27))                                              &
      +(RC(:,131)*y(:,20)*y(:,15)*0.8 )+(DJ(:,8)  *y(:,22))             &
      +(RC(:,91) *y(:,20)*y(:,20)*2.00)+(RC(:,123)*y(:,4 )*y(:,50)*0.58)
      p =p                                                              &
      +(RC(:,79) *y(:,20)*y(:,5 ))+(RC(:,80) *y(:,16)*y(:,20))          &
      +(RC(:,44) *y(:,3 )*y(:,17))+(RC(:,59) *y(:,3 )*y(:,12))          
      l = 0.0                                                           &
      +(DW(:,15) )                                                      &
      +(RC(:,132)*y(:,46))+(RC(:,133)*y(:,48))+(RC(:,181)*y(:,44))
      l = l                                                             &
      +(RC(:,106)*y(:,37))+(RC(:,111)*y(:,38))+(RC(:,127)*y(:,40))      &
      +(RC(:,84) *y(:,24))+(RC(:,96) *y(:,36))+(RC(:,97) *y(:,30))
      l = l                                                             &
      +(RC(:,73) *y(:,18))+(RC(:,74) *y(:,20))+(RC(:,80) *y(:,20))      &
      +(RC(:,62) *y(:,16))
      l = l                                                             &
      +(RC(:,62) *y(:,16))+(RC(:,65) *y(:,15))                          &
      +(RC(:,60) *y(:,5 ))+(RC(:,61) *y(:,16))+(RC(:,61) *y(:,16))      

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,16) = (yp(:,16)+dts*p)*tmp1
!
!          MeOOH        y(17)
      p = 0.0                                                           &
      +(RC(:,65) *y(:,16)*y(:,15))                                      
      l = 0.0                                                           &
      +(DW(:,17) )        +(DD(:,17) )                                  &
      +(RC(:,44) *y(:,3 ))+(RC(:,45) *y(:,3 ))+(DJ(:,17) )     

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,17) = (yp(:,17)+dts*p)*tmp1
!
!          EtOO         y(18)
      p = 0.0                                                           &
      +(RC(:,151)*y(:,6 )*y(:,19))+(DJ(:,9)  *y(:,53))                  &
      +(RC(:,44) *y(:,3 )*y(:,21))+(RC(:,71) *y(:,3 )*y(:,19))          
      l = 0.0                                                           &
      +(RC(:,90) *y(:,18))+(RC(:,99) *y(:,15))                          &
      +(RC(:,72) *y(:,5 ))+(RC(:,73) *y(:,16))+(RC(:,90) *y(:,18))   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,18) = (yp(:,18)+dts*p)*tmp1
!
!          C2H6         y(19)
      p = 0.0                                                           
      l = 0.0                                                           &
      +(RC(:,71) *y(:,3 ))+(RC(:,151)*y(:,6 ))  

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,19) = (yp(:,19)+dts*p)*tmp1
!
!          MeCO3        y(20)
      p = 0.0                                                           &
      +(DJ(:,18) *y(:,23))                                              &
      +(DJ(:,10) *y(:,27))        +(DJ(:,12) *y(:,33))                  &
      +(RC(:,158)*y(:,6 )*y(:,22))+(DJ(:,9)  *y(:,53))                  &
      +(RC(:,131)*y(:,20)*y(:,15)*0.2 )
      p = p                                                             &
      +(RC(:,138)*y(:,33)*y(:,3 ))                                      &
      +(RC(:,105)*y(:,5 )*y(:,37))+(RC(:,106)*y(:,16)*y(:,37))
      p = p                                                             &
      +(RC(:,95) *y(:,36)*y(:,5 ))+(RC(:,96) *y(:,36)*y(:,16))          &
      +(RC(:,75) *y(:,3 )*y(:,22))+(RC(:,78) *y(:,23))                  
      l = 0.0                                                           &
      +(RC(:,131)*y(:,15))                                              &
      +(RC(:,80) *y(:,16))+(RC(:,91) *y(:,20))+(RC(:,91) *y(:,20))      &
      +(RC(:,74) *y(:,16))+(RC(:,77) *y(:,7 ))+(RC(:,79) *y(:,5 ))      

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,20) = (yp(:,20)+dts*p)*tmp1
!
!          EtOOH        y(21)
      p = 0.0                                                           &
      +(RC(:,99) *y(:,15)*y(:,18))                                      
      l = 0.0                                                           &
      +(DW(:,21) )        +(DD(:,21) )                                  &
      +(RC(:,44) *y(:,3 ))+(RC(:,100)*y(:,3 ))+(DJ(:,19) )   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,21) = (yp(:,21)+dts*p)*tmp1
!
!          MeCHO        y(22)
      p = 0.0                                                           &
      +(DJ(:,19) *y(:,21))                                              &
      +(RC(:,127)*y(:,16)*y(:,40))+(RC(:,156)*y(:,47)*y(:,3 ))          &
      +(RC(:,124)*y(:,4 )*y(:,50))+(RC(:,126)*y(:,40)*y(:,5 ))
      p = p                                                             &
      +(RC(:,105)*y(:,5 )*y(:,37))+(RC(:,106)*y(:,16)*y(:,37))          &
      +(RC(:,90) *y(:,18)*y(:,18)*2.00)+(RC(:,100)*y(:,3 )*y(:,21))     &
      +(RC(:,72) *y(:,18)*y(:,5 ))+(RC(:,73) *y(:,18)*y(:,16))          
      l = 0.0                                                           &
      +(RC(:,75) *y(:,3 ))+(RC(:,158)*y(:,6 ))+(DJ(:,8)  )    

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,22) = (yp(:,22)+dts*p)*tmp1
!
!          PAN          y(23)
      p = 0.0                                                           &
      +(RC(:,77) *y(:,20)*y(:,7 ))                                      
      l = 0.0                                                           &
      +(DD(:,23) )                                                      &
      +(RC(:,78) )        +(RC(:,98) *y(:,3 ))+(DJ(:,18) )    

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,23) = (yp(:,23)+dts*p)*tmp1
!
!          s-BuOO       y(24)
      p = 0.0                                                           &
      +(RC(:,152)*y(:,6 )*y(:,51))                                      &
      +(RC(:,44) *y(:,3 )*y(:,52))+(RC(:,81) *y(:,3 )*y(:,51))          
      l = 0.0                                                           &
      +(RC(:,83) *y(:,5 ))+(RC(:,84) *y(:,16))+(RC(:,103)*y(:,15))      

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,24) = (yp(:,24)+dts*p)*tmp1
!
!          C3H8         y(25)
      p = 0.0                                                           
      l = 0.0                                                           &
      +(RC(:,92) *y(:,3 ))                                              

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,25) = (yp(:,25)+dts*p)*tmp1
!
!          i-PrOOH      y(26)
      p = 0.0                                                           &
      +(RC(:,101)*y(:,15)*y(:,30))                                      
      l = 0.0                                                           &
      +(DW(:,26) )        +(DD(:,26) )                                  &
      +(RC(:,44) *y(:,3 ))+(RC(:,102)*y(:,3 ))+(DJ(:,20) ) 

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,26) = (yp(:,26)+dts*p)*tmp1
!
!          Me2CO        y(27)
      p = 0.0                                                           &
      +(RC(:,102)*y(:,3 )*y(:,26))+(DJ(:,20) *y(:,26))                  &
      +(RC(:,93) *y(:,30)*y(:,5 ))+(RC(:,97) *y(:,30)*y(:,16))          
      l = 0.0                                                           &
      +(RC(:,94) *y(:,3 ))+(DJ(:,10) )    

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,27) = (yp(:,27)+dts*p)*tmp1
!
!          O3S          y(28)
      p = 0.0                                                           
      l = 0.0                                                           &
      +(DJ(:,1)  )        +(DJ(:,2)  )        +(DD(:,28)  )             &
      +(RC(:,124)*y(:,50))+(RC(:,128)*y(:,29))+(RC(:,129)*y(:,34))
      l = l                                                             &
      +(RC(:,14) *y(:,15))+(RC(:,112)*y(:,49))+(RC(:,123)*y(:,50))      &
      +(RC(:,11) *y(:,5 ))+(RC(:,12) *y(:,7 ))+(RC(:,13) *y(:,3 )) 

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,28) = (yp(:,28)+dts*p)*tmp1
!
!          C5H8         y(29)
      p = 0.0                                                           
      l = 0.0                                                           &
      +(RC(:,128)*y(:,4 ))+(RC(:,141)*y(:,3 ))+(RC(:,159)*y(:,6 ))    

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,29) = (yp(:,29)+dts*p)*tmp1
!
!          i-PrOO       y(30)
      p = 0.0                                                           &
      +(RC(:,44) *y(:,3 )*y(:,26))+(RC(:,92) *y(:,3 )*y(:,25))          
      l = 0.0                                                           &
      +(RC(:,93) *y(:,5 ))+(RC(:,97) *y(:,16))+(RC(:,101)*y(:,15))

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,30) = (yp(:,30)+dts*p)*tmp1
!
!          ISOOH        y(31)
      p = 0.0                                                           &
      +(RC(:,134)*y(:,46)*y(:,15))                                      
      l = 0.0                                                           &
      +(DD(:,31) )                                                      &
      +(RC(:,135)*y(:,3 ))+(DJ(:,22) )        +(DW(:,31) )       

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,31) = (yp(:,31)+dts*p)*tmp1
!
!          ISON         y(32)
      p = 0.0                                                           &
      +(RC(:,159)*y(:,6 )*y(:,29))                                      
      l = 0.0                                                           &
      +(RC(:,160)*y(:,3 ))+(DW(:,32) )        

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,32) = (yp(:,32)+dts*p)*tmp1
!
!          MGLY         y(33)
      p = 0.0                                                           &
      +(DJ(:,23) *y(:,35))                                              &
      +(RC(:,180)*y(:,42)*y(:,15))+(RC(:,181)*y(:,44)*y(:,16))          &
      +(RC(:,170)*y(:,3 )*y(:,58)*0.80)+(RC(:,173)*y(:,44)*y(:,5 ))
      p = p                                                             &
      +(RC(:,137)*y(:,35)*y(:,3 ))+(RC(:,144)*y(:,48)*y(:,5 ))          &
      +(RC(:,129)*y(:,4 )*y(:,34))+(RC(:,133)*y(:,48)*y(:,16)*1.0 )     
      l = 0.0                                                           &
      +(RC(:,138)*y(:,3 ))+(DJ(:,12) )        +(DW(:,33) )   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,33) = (yp(:,33)+dts*p)*tmp1
!
!          MVK          y(34)
      p = 0.0                                                           &
      +(RC(:,160)*y(:,32)*y(:,3 ))+(DJ(:,22) *y(:,31))                  &
      +(RC(:,135)*y(:,31)*y(:,3 ))+(RC(:,142)*y(:,46)*y(:,5 ))          &
      +(RC(:,128)*y(:,4 )*y(:,29))+(RC(:,132)*y(:,46)*y(:,16)*1.0 )     
      l = 0.0                                                           &
      +(RC(:,129)*y(:,4 ))+(RC(:,143)*y(:,3 ))  

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,34) = (yp(:,34)+dts*p)*tmp1
!
!          MVKOOH       y(35)
      p = 0.0                                                           &
      +(RC(:,136)*y(:,48)*y(:,15))                                      
      l = 0.0                                                           &
      +(DD(:,35) )                                                      &
      +(RC(:,137)*y(:,3 ))+(DJ(:,23) )        +(DW(:,35) )    

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,35) = (yp(:,35)+dts*p)*tmp1
!
!          ACETO2       y(36)
      p = 0.0                                                           &
      +(RC(:,94) *y(:,3 )*y(:,27))                                      
      l = 0.0                                                           &
      +(RC(:,95) *y(:,5 ))+(RC(:,96) *y(:,16)) 

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,36) = (yp(:,36)+dts*p)*tmp1
!
!          MEKO2        y(37)
      p = 0.0                                                           &
      +(RC(:,86) *y(:,3 )*y(:,53))                                      
      l = 0.0                                                           &
      +(RC(:,105)*y(:,5 ))+(RC(:,106)*y(:,16))    

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,37) = (yp(:,37)+dts*p)*tmp1
!
!          HOC2H4O2     y(38)
      p = 0.0                                                           &
      +(RC(:,109)*y(:,3 )*y(:,49))                                      
      l = 0.0                                                           &
      +(RC(:,110)*y(:,5 ))+(RC(:,111)*y(:,16))    

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,38) = (yp(:,38)+dts*p)*tmp1
!
!          ORGNIT       y(39)
      p = 0.0                                                           &
      +(RC(:,171)*y(:,7 )*y(:,42))+(RC(:,183)*y(:,7 )*y(:,55))          
      l = 0.0                                                           &
      +(DD(:,39) )                                                      &
      +(RC(:,177)*y(:,3 ))+(RC(:,178)*y(:,6 ))+(DW(:,39) )              

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,39) = (yp(:,39)+dts*p)*tmp1
!
!          HOC3H6O2     y(40)
      p = 0.0                                                           &
      +(RC(:,125)*y(:,3 )*y(:,50))                                      
      l = 0.0                                                           &
      +(RC(:,126)*y(:,5 ))+(RC(:,127)*y(:,16))

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,40) = (yp(:,40)+dts*p)*tmp1
!
!          CH3OH        y(41)
      p = 0.0                                                           &
      +(RC(:,62) *y(:,16)*y(:,16))+(RC(:,123)*y(:,4 )*y(:,50)*0.12)     
      l = 0.0                                                           &
      +(RC(:,63) *y(:,3 ))+(DW(:,41) )

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,41) = (yp(:,41)+dts*p)*tmp1
!
!          OXYL1        y(42)
      p = 0.0                                                           &
      +(RC(:,176)*y(:,3 )*y(:,58))                                      
      l = 0.0                                                           &
      +(RC(:,171)*y(:,7 ))+(RC(:,180)*y(:,15))   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,42) = (yp(:,42)+dts*p)*tmp1
!
!          H2           y(43)
      p = 0.0                                                           &
      +(DJ(:,7)  *y(:,14))                                              &
      +(RC(:,112)*y(:,4 )*y(:,49)*0.13)+(RC(:,124)*y(:,4 )*y(:,50)*0.24)
      l = 0.0                                                           &
      +(RC(:,33) *y(:,3 ))+(DD(:,43) )    

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,43) = (yp(:,43)+dts*p)*tmp1
!
!          MEMALD1      y(44)
      p = 0.0                                                           &
      +(RC(:,172)*y(:,3 )*y(:,56))                                      
      l = 0.0                                                           &
      +(RC(:,173)*y(:,5 ))+(RC(:,181)*y(:,16))    

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,44) = (yp(:,44)+dts*p)*tmp1
!
!          RNC2H4       y(45)
      p = 0.0                                                           &
      +(RC(:,153)*y(:,6 )*y(:,49))                                      
      l = 0.0                                                           &
      +(RC(:,154)*y(:,3 ))  

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,45) = (yp(:,45)+dts*p)*tmp1
!
!          HOIPO2       y(46)
      p = 0.0                                                           &
      +(RC(:,141)*y(:,3 )*y(:,29))                                      
      l = 0.0                                                           &
      +(RC(:,132)*y(:,16))+(RC(:,134)*y(:,15))+(RC(:,142)*y(:,5 ))     

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,46) = (yp(:,46)+dts*p)*tmp1
!
!          RNC3H6       y(47)
      p = 0.0                                                           &
      +(RC(:,155)*y(:,6 )*y(:,50))                                      
      l = 0.0                                                           &
      +(RC(:,156)*y(:,3 ))               

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,47) = (yp(:,47)+dts*p)*tmp1
!
!          HOMVKO2      y(48)
      p = 0.0                                                           &
      +(RC(:,143)*y(:,3 )*y(:,34))                                      
      l = 0.0                                                           &
      +(RC(:,133)*y(:,16))+(RC(:,136)*y(:,15))+(RC(:,144)*y(:,5 ))     

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,48) = (yp(:,48)+dts*p)*tmp1
!
!          C2H4         y(49)
      p = 0.0                                                           
      l = 0.0                                                           &
      +(RC(:,109)*y(:,3 ))+(RC(:,112)*y(:,4 ))+(RC(:,153)*y(:,6 ))   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,49) = (yp(:,49)+dts*p)*tmp1
!
!          C3H6         y(50)
      p = 0.0                                                           
      l = 0.0                                                           &
      +(RC(:,155)*y(:,6 ))                                              &
      +(RC(:,123)*y(:,4 ))+(RC(:,124)*y(:,4 ))+(RC(:,125)*y(:,3 )) 

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,50) = (yp(:,50)+dts*p)*tmp1
!
!          NC4H10       y(51)
      p = 0.0                                                           
      l = 0.0                                                           &
      +(RC(:,81) *y(:,3 ))+(RC(:,152)*y(:,6 ))       

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,51) = (yp(:,51)+dts*p)*tmp1
!
!          s-BuOOH      y(52)
      p = 0.0                                                           &
      +(RC(:,103)*y(:,15)*y(:,24))                                      
      l = 0.0                                                           &
      +(DW(:,52) )        +(DD(:,52) )                                  &
      +(RC(:,44) *y(:,3 ))+(RC(:,104)*y(:,3 ))+(DJ(:,21) )  

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,52) = (yp(:,52)+dts*p)*tmp1
!
!          MEK          y(53)
      p = 0.0                                                           &
      +(RC(:,104)*y(:,3 )*y(:,52))+(DJ(:,21) *y(:,52))                  &
      +(RC(:,83) *y(:,24)*y(:,5 ))+(RC(:,84) *y(:,24)*y(:,16))          
      l = 0.0                                                           &
      +(RC(:,86) *y(:,3 ))+(DJ(:,9)  )            

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,53) = (yp(:,53)+dts*p)*tmp1
!
!          TOLUEN       y(54)
      p = 0.0                                                           
      l = 0.0                                                           &
      +(RC(:,174)*y(:,3 ))+(RC(:,175)*y(:,3 )) 

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,54) = (yp(:,54)+dts*p)*tmp1
!
!          TOLP1        y(55)
      p = 0.0                                                           &
      +(RC(:,175)*y(:,3 )*y(:,54))                                      
      l = 0.0                                                           &
      +(RC(:,182)*y(:,15))+(RC(:,183)*y(:,7 )) 

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,55) = (yp(:,55)+dts*p)*tmp1
!
!          MEMALD       y(56)
      p = 0.0                                                           &
      +(RC(:,180)*y(:,42)*y(:,15))+(RC(:,182)*y(:,15)*y(:,55))          &
      +(RC(:,177)*y(:,3 )*y(:,39))+(RC(:,178)*y(:,6 )*y(:,39))          &
      +(RC(:,170)*y(:,3 )*y(:,58)*0.80)+(RC(:,174)*y(:,3 )*y(:,54))     
      l = 0.0                                                           &
      +(RC(:,172)*y(:,3 ))               

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,56) = (yp(:,56)+dts*p)*tmp1
!
!          GLY          y(57)
      p = 0.0                                                           &
      +(RC(:,181)*y(:,44)*y(:,16))+(RC(:,182)*y(:,15)*y(:,55))          &
      +(RC(:,177)*y(:,3 )*y(:,39))+(RC(:,178)*y(:,6 )*y(:,39))          &
      +(RC(:,173)*y(:,44)*y(:,5 ))+(RC(:,174)*y(:,3 )*y(:,54))          
      l = 0.0                                                           &
      +(RC(:,139)*y(:,3 ))+(DJ(:,13) )        +(DW(:,57) )

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,57) = (yp(:,57)+dts*p)*tmp1
!
!          OXYL         y(58)
      p = 0.0                                                           
      l = 0.0                                                           &
      +(RC(:,170)*y(:,3 ))+(RC(:,176)*y(:,3 ))   

      tmp2(:) = 1.0+dts*l(:)

      CALL oneover_v(n_pnts, tmp2, tmp1)

      y(:,58) = (yp(:,58)+dts*p)*tmp1
!
        END DO ! End of iteration loop

      END DO  ! n_be_calls

      IF (lhook) CALL dr_hook('UKCA_DERIV_RAQ',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_DERIV_RAQ
