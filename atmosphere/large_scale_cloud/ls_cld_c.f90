! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Large-scale Cloud Scheme.
! Subroutine Interface:
! ======================================================================

!+ Large-scale Cloud Scheme Compression routine (Cloud points only).
! Subroutine Interface:
SUBROUTINE ls_cld_c(                                                    &
 p_f,rhcrit,qsl_f,qn_f,q_f,t_f,                                         &
 qcl_f,cf_f,grid_qc_f,bs_f,                                             &
 INDEX,points,rhc_row_length,rhc_rows,                                  &
 bl_levels,k, l_mixing_ratio)

  USE water_constants_mod,  ONLY: lc
  USE atmos_constants_mod,  ONLY: cp, r, repsilon
  USE yomhook,              ONLY: lhook, dr_hook
  USE parkind1,             ONLY: jprb, jpim
  USE atm_fields_bounds_mod,ONLY: qdims
  USE cloud_inputs_mod,     ONLY: l_eacf

  IMPLICIT NONE

! Purpose: Calculates liquid cloud water amounts and cloud amounts,
!          temperature and specific humidity from cloud-conserved and
!          other model variables. This is done for one model level.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No.29

!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
  INTEGER ::                                                            &
                        !, INTENT(IN)
    rhc_row_length,rhc_rows,                                            &
!       No. of gridpoints being processed.
    bl_levels,                                                          &
!       No. of boundary layer levels
    k,                                                                  &
!       Level no.
   points,                                                              &
!       No. of gridpoints with non-zero cloud
   index((1+qdims%j_end-qdims%j_start)*(1+qdims%i_end-qdims%i_start),2)
!       INDEX for  points with non-zero cloud from lowest model level.

  REAL ::                                                               &
                        !, INTENT(IN)
   rhcrit(rhc_row_length,rhc_rows),                                     &
!       Critical relative humidity.  See the paragraph incorporating
!       eqs P292.11 to P292.14.
   p_f(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       pressure (Pa).
   qsl_f(         qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       saturated humidity at temperature TL, and pressure P_F
   qn_f(          qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end)
!       Normalised super/subsaturation ( = QC/BS).

  LOGICAL ::                                                            &
                        !, INTENT(IN)
   l_mixing_ratio   !  Use mixing ratio formulation

  REAL ::                                                               &
                        !, INTENT(INOUT)
   q_f(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       On input : Vapour + liquid water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                   (kg water per kg air).
   t_f(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end)
!       On input : Liquid water temperature (TL) (K).
!       On output: Temperature at processed levels (K).

  REAL ::                                                               &
                        !, INTENT(OUT)
   qcl_f(         qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       Cloud liquid water content at processed levels (kg per kg air).
   cf_f(          qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       Liquid cloud fraction at processed levels.
   grid_qc_f(     qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       Super/subsaturation on processed levels. Input initially RMDI.
   bs_f(          qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end)
!       Value of bs at processed levels. Input initialized to RMDI.

!  Local parameters and other physical constants------------------------
  REAL :: alphl,lcrcp                  ! Derived parameters.
  PARAMETER (                                                           &
   alphl=repsilon*lc/r,                                                 &
                                        ! For liquid AlphaL calculation.
   lcrcp=lc/cp                                                          &
                                        ! Lat ht of condensation/Cp.
  )
  REAL :: wtn                          ! Weighting for ALPHAL iteration
  INTEGER ::                                                            &
   its                              ! Total number of iterations
  PARAMETER (its=5,wtn=0.75)

!  Local scalars--------------------------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
  REAL ::                                                               &
   al,                                                                  &
                         ! LOCAL AL (see equation P292.6).
   alphal,                                                              &
                         ! LOCAL ALPHAL (see equation P292.5).
   qn_adj,                                                              &
   rhcritx          ! scalar copy of RHCRIT(I,J)
  INTEGER ::                                                            &
   multrhc          ! Zero if (rhc_row_length*rhc_rows) le 1, else 1

!  (b) Others.
  INTEGER ::   i,ii,ij,n   ! Loop counters:I,II-horizontal field index.
!                                       : N - iteration counter.


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!  Local dynamic arrays-------------------------------------------------
!    7 blocks of real workspace are required.
  REAL ::                                                               &
   p(points),                                                           &
!       Pressure  (Pa).
   qs(points),                                                          &
!       Saturated spec humidity for temp T.
   qcn(points),                                                         &
!       Cloud water normalised with BS.
   t(points),                                                           &
!       temperature.
   q(points),                                                           &
!       specific humidity.
   bs(points),                                                          &
!       Sigmas*sqrt(6): sigmas the parametric standard deviation of
!       local cloud water content fluctuations.
   alphal_nm1(points)
!       ALPHAL at previous iteration.

!  External subroutine calls: ------------------------------------------
  EXTERNAL qsat_wat

!- End of Header


  IF (lhook) CALL dr_hook('LS_CLD_C',zhook_in,zhook_handle)

! ==Main Block==--------------------------------------------------------
! Operate on INDEXed points with non-zero cloud fraction.
! ----------------------------------------------------------------------
  IF ( (rhc_row_length * rhc_rows)  >   1) THEN
    multrhc = 1
  ELSE
    multrhc = 0
  END IF

!        RHCRITX = RHCRIT(1,1)
! Points_do1:
!CDIR NODEP
  DO i = 1, points
    ii = INDEX(i,1)
    ij = INDEX(i,2)
    IF ( multrhc== 1) THEN
      rhcritx = rhcrit(ii,ij)
    ELSE
      rhcritx = rhcrit(1,1)
    END IF
    p(i)  = p_f(ii,ij)
    qcn(i)= qn_f(ii,ij)
! ----------------------------------------------------------------------
! 1. Calculate ALPHAL (eq P292.5) and AL (P292.6).
!    CAUTION: T_F acts as TL (input value) until update in final section
!    CAUTION: Q_F acts as QW (input value) until update in final section
! ----------------------------------------------------------------------

    alphal = alphl * qsl_f(ii,ij) / (t_f(ii,ij) * t_f(ii,ij)) !P292.5
    al = 1.0 / (1.0 + (lcrcp * alphal))                    ! P292.6
    alphal_nm1(i) = alphal

! Rhcrit_if1:
    IF (rhcritx  <   1.) THEN
! ----------------------------------------------------------------------
! 2. Calculate BS (ie. sigma*sqrt(6), where sigma is
!    as in P292.14) and normalised cloud water QCN=qc/BS, using eqs
!    P292.15 & 16 if RHcrit < 1.
! N.B. QN (input) is initially in QCN
! N.B. QN does not depend on AL and so CF and QCN can be calculated
!      outside the iteration (which is performed in LS_CLD_C).
!      QN is > -1 for all points processed so CF > 0.
! ----------------------------------------------------------------------

      bs(i) = (1.0-rhcritx) * al * qsl_f(ii,ij)  ! P292.14
      IF (qcn(i)  <=  0.) THEN
        cf_f(ii,ij) = 0.5 * (1. + qcn(i)) * (1. + qcn(i))
        qcn(i)= (1. + qcn(i)) * (1. + qcn(i)) * (1. + qcn(i)) / 6.
      ELSE IF (qcn(i)  <   1.) THEN
        cf_f(ii,ij) = 1. - 0.5 * (1. - qcn(i)) * (1. - qcn(i))
        qcn(i)=qcn(i) + (1.-qcn(i)) * (1.-qcn(i)) * (1.-qcn(i))/6.
      ELSE ! QN  >=  1
        cf_f(ii,ij) = 1.
      END IF ! QCN_if

! ----------------------------------------------------------------------
! 3.b If necessary, modify cloud fraction using empirically adjusted
!     cloud fraction parametrization, but keep liquid content the same.
! ----------------------------------------------------------------------
      IF (l_eacf) THEN
            ! Adjust QN according to EACF parametrization

        IF (k <= bl_levels) THEN
          qn_adj=(qn_f(ii,ij)+0.184)/(1.-0.184)
        ELSE
          qn_adj=(qn_f(ii,ij)+0.0955)/(1.-0.0955)
        END IF

!         Calculate cloud fraction using adjusted QN
        IF (qn_adj  <=  0.) THEN
          cf_f(ii,ij) = 0.5 * (1.0 + qn_adj) * (1.0 + qn_adj)
        ELSE IF (qn_adj  <   1.) THEN
          cf_f(ii,ij) = 1. - 0.5 * (1.0 - qn_adj) * (1.0 - qn_adj)
        ELSE ! QN_ADJ  >=  1
          cf_f(ii,ij) = 1.
        END IF ! QN_ADJ_if

      END IF  ! l_eacf

    ELSE ! i.e. if RHcrit = 1
! ----------------------------------------------------------------------
! 3.a If RHcrit = 1., all points processed have QN > 0 and CF = 1.
! ----------------------------------------------------------------------
      bs(i) = al
      cf_f(ii,ij) = 1.
    END IF ! Rhcrit_if1

! ----------------------------------------------------------------------
! 3.1 Calculate 1st approx. to qc (store in QCL)
! ----------------------------------------------------------------------

    qcl_f(ii,ij) = qcn(i) * bs(i)

! ----------------------------------------------------------------------
! 3.2 Calculate 1st approx. specific humidity (total minus cloud water)
! ----------------------------------------------------------------------

    q(i) = q_f(ii,ij) - qcl_f(ii,ij)

! ----------------------------------------------------------------------
! 3.3 Calculate 1st approx. to temperature, adjusting for latent heating
! ----------------------------------------------------------------------

    t(i) = t_f(ii,ij) + lcrcp*qcl_f(ii,ij)
  END DO ! Points_do1

! ----------------------------------------------------------------------
! 4. Iteration to find better cloud water values.
! ----------------------------------------------------------------------
! Its_if:
  IF (its  >=  2) THEN
! Its_do:
    DO n = 2, its

! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qs,t,p,points,l_mixing_ratio)
! Points_do2:
      rhcritx = rhcrit(1,1)
      DO i = 1, points
        ii = INDEX(i,1)
        ij = INDEX(i,2)
        IF ( multrhc== 1) THEN
          rhcritx = rhcrit(ii,ij)
        ELSE
          rhcritx = rhcrit(1,1)
        END IF
! T_if:
        IF (t(i)  >   t_f(ii,ij)) THEN
!           NB. T > TL implies cloud fraction > 0.
          alphal = (qs(i) - qsl_f(ii,ij)) / (t(i) - t_f(ii,ij))
          alphal = wtn * alphal + (1.0 - wtn) * alphal_nm1(i)
          alphal_nm1(i) = alphal
          al = 1.0 / (1.0 + (lcrcp * alphal))
! Rhcrit_if2:
          IF (rhcritx  <   1.) THEN
            bs(i) = (1.0-rhcritx) * al * qsl_f(ii,ij)
!                                                             P292.14
          ELSE
            bs(i) = al
          END IF  ! Rhcrit_if2

! ----------------------------------------------------------------------
! 4.1 Calculate Nth approx. to qc (store in QCL).
! ----------------------------------------------------------------------

          qcl_f(ii,ij) = qcn(i) * bs(i)

! ----------------------------------------------------------------------
! 4.2 Calculate Nth approx. spec. humidity (total minus cloud water).
! ----------------------------------------------------------------------

          q(i) = q_f(ii,ij) - qcl_f(ii,ij)

! ----------------------------------------------------------------------
! 4.3 Calculate Nth approx. to temperature, adjusting for latent heating
! ----------------------------------------------------------------------

          t(i) = t_f(ii,ij) + lcrcp * qcl_f(ii,ij)

        END IF ! T_if
      END DO ! Points_do2
    END DO ! Its_do
  END IF ! Its_if

! ----------------------------------------------------------------------
! 5. Finally scatter back cloud point results to full field arrays.
!    CAUTION: T_F updated from TL (input) to T (output)
!    CAUTION: Q_F updated from QW (input) to Q (output)
! ----------------------------------------------------------------------

!DIR$ IVDEP
! Points_do3:
  DO i = 1, points
    ii = INDEX(i,1)
    ij = INDEX(i,2)
    q_f(ii,ij) = q(i)
    t_f(ii,ij) = t(i)
    grid_qc_f(ii,ij) = bs(i) * qn_f(ii,ij)
    bs_f(ii,ij) = bs(i)
  END DO ! Points_do3


  IF (lhook) CALL dr_hook('LS_CLD_C',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE ls_cld_c
! ======================================================================
