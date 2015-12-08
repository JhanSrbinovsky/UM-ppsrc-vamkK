! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Large-scale Cloud Scheme.
! Subroutine Interface:
SUBROUTINE ls_cld(                                                      &
!      Pressure related fields
 p_theta_levels, rhcrit,                                                &
!      Array dimensions
 levels, bl_levels,                                                     &
 rhc_row_length,rhc_rows,                                               &
!      From convection diagnosis (only used if A05_4A)
 ntml, cumulus, l_mixing_ratio,                                         &
!      Prognostic Fields
 t, cf, q, qcf, qcl,                                                    &
!      Liquid and frozen ice cloud fractions
 cfl, cff,                                                              &
 error)

  USE cv_run_mod,            ONLY: l_param_conv
  USE vectlib_mod,           ONLY: powr_v
  USE conversions_mod,       ONLY: pi
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: pdims, qdims, tdims
  USE cloud_inputs_mod,      ONLY: cloud_fraction_method,               &
   ice_fraction_method, l_eacf, overlap_ice_liquid, ctt_weight,         &
   t_weight, qsat_fixed, sub_cld

  USE c_cldsgs_mod, ONLY: qcfmin
  IMPLICIT NONE
!
! Purpose:
!   This subroutine calculates liquid and ice cloud fractional cover
!   for use with the enhanced precipitation microphysics scheme.

! Method:
!   Statistical cloud scheme separates input moisture into specific
!   humidity and cloud liquid water. Temperature calculated from liquid
!   water temperature. Cloud fractions calculated from statistical
!   relation between cloud fraction and cloud liquid/ice water content.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP No. 29


!  Global Variables:----------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
  INTEGER ::                                                            &
                        !, INTENT(IN)
   levels,                                                              &
!       No. of levels being processed.
   bl_levels,                                                           &
!       No. of boundary layer levels
   rhc_row_length,rhc_rows

  REAL ::                                                               &
                        !, INTENT(IN)
   qcf(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,levels),                    &
!       Cloud ice content at processed levels (kg water per kg air).
   p_theta_levels(pdims%i_start:pdims%i_end,                            & 
                  pdims%j_start:pdims%j_end,levels),                    &
!       pressure at all points (Pa).
   rhcrit(rhc_row_length,rhc_rows,levels)
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.

  INTEGER ::                                                            &
   ntml(          qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end)     
!       IN Height of diagnosed BL top

  LOGICAL ::                                                            &
   l_mixing_ratio  
!       IN true if using mixing ratios

  LOGICAL ::                                                            &
   cumulus(       qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end)  
!       IN Logical indicator of convection

  REAL ::                                                               &
                        !, INTENT(INOUT)
   q(             qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,levels),                    &
!       On input : Total water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels 
!                  (kg water per kg air).
   t(             tdims%i_start:tdims%i_end,                            & 
                  tdims%j_start:tdims%j_end,levels)
!       On input : Liquid/frozen water temperature (TL) (K).
!       On output: Temperature at processed levels (K).

  REAL ::                                                               &
                        !, INTENT(OUT)
   cf(            qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,levels),                    &
!       Cloud fraction at processed levels (decimal fraction).
   qcl(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,levels),                    &
!       Cloud liquid water content at processed levels (kg per kg air).
   cfl(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,levels),                    &
!       Liquid cloud fraction at processed levels (decimal fraction).
   cff(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,levels)
!       Frozen cloud fraction at processed levels (decimal fraction).

!     Error Status:
  INTEGER :: error     !, INTENT(OUT)  0 if OK; 1 if bad arguments.

!  Local parameters and other physical constants------------------------
  REAL :: rootwo       ! Sqrt(2.)
  REAL :: subgrid      ! Subgrid parameter in ice cloud calculation

!  Local scalars--------------------------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
  REAL ::                                                               &
   phiqcf,                                                              &
                        ! Arc-cosine term in Cloud ice fraction calc.
   cosqcf,                                                              &
                        ! Cosine term in Cloud ice fraction calc.
   overlap_max,                                                         &
                        ! Maximum possible overlap
   overlap_min,                                                         &
                        ! Minimum possible overlap
   overlap_random,                                                      &
                        ! Random overlap
   temp0,                                                               &
   temp1,                                                               &
   temp2,                                                               &
                        ! Temporaries for combinations of the
   qn_imp,                                                              &
   qn_adj
!                       ! overlap parameters

!  (b) Others.
  INTEGER :: k,i,j       ! Loop counters: K - vertical level index.
!                                        I,J - horizontal field indices.

  INTEGER :: qc_points,                                                 &
                        ! No. points with non-zero cloud
          multrhc   ! Zero if (rhc_row_length*rhc_rows) le 1, else 1

!  Local dynamic arrays-------------------------------------------------
!    6 blocks of real workspace are required.
  REAL ::                                                               &
   qcfrbs(        qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       qCF / bs
   qsl(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       Saturated specific humidity for temp TL or T.
   qsl_ctt(       qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       Saturated specific humidity wrt liquid at cloud top temperature
   qn(            qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end),                           &
!       Cloud water normalised with BS.
   grid_qc(       qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,levels),                    &
!       Gridbox mean saturation excess at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
   bs(            qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,levels),                    &
!       Maximum moisture fluctuation /6*sigma at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
   ctt(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,levels)
!       Ice cloud top temperature (K) - as coded it is really TL
  LOGICAL ::                                                            &
   lqc(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end)      
!       True for points with non-zero cloud
  INTEGER ::                                                            &
   index((1+qdims%j_end-qdims%j_start)*(1+qdims%i_end-qdims%i_start),2),&
!       Index for points with non-zero cloud
   llwic(         qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end)
!       Last Level With Ice Cloud
  REAL :: rhcritx              
!       Scalar copy of RHCRIT(I,J,K)

!       Variables for cache-blocking
  INTEGER            :: jj             !Block index
  INTEGER            :: kk

!       jblock is not a parameter. Value needs to be flexible at runtime.
  INTEGER            :: jblock

  INTEGER, PARAMETER :: kblock=4       !Block size

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!- End of Header

!Set value of blocking factor, jblock. May need to be changed on other 
! platforms or for OMP threads > 4
  jblock = 4

! ----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
! ----------------------------------------------------------------------
  error=0

  IF ( (rhc_row_length * rhc_rows)  >   1) THEN
    multrhc = 1
  ELSE
    multrhc = 0
  END IF

! ==Main Block==--------------------------------------------------------
! Subroutine structure :
! Loop round levels to be processed.
! ----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('LS_CLD',zhook_in,zhook_handle)

! Initialize cloud-top-temperature and last-level-with-ice-cloud arrays

! Levels_do1:
!Parameters used: kblock
!$OMP  PARALLEL DEFAULT(NONE)                                                   &
!$OMP& SHARED(llwic, ctt, levels, qcl, cfl, grid_qc, bs, p_theta_levels,        &
!$OMP& l_mixing_ratio, qcf, ice_fraction_method, ctt_weight, l_eacf,            &
!$OMP& rhc_row_length, rhc_rows, bl_levels, cloud_fraction_method, cf,          &
!$OMP& overlap_ice_liquid, cff, t_weight, sub_cld, q, t, jblock,                &
!$OMP& qsat_fixed, multrhc, cumulus, ntml, rhcrit, l_param_conv, qdims)             &
!$OMP& PRIVATE(k, j, i, rhcritx, qc_points, rootwo, subgrid, qsl, qsl_ctt,      &
!$OMP& phiqcf, cosqcf, qn_imp, qn_adj, overlap_max, overlap_min, overlap_random,&
!$OMP& temp0, temp1, temp2, qn, lqc, index, qcfrbs, kk, jj)

!Cache-blocking applied to loop over j. It is used here to allow
!OpenMP parallelism to be over j (looping over k is order-dependent),
!but still have j and k in the correct order for good cache use.
!If jblock=rows=(1+qdims%j_end-qdims%j_start), 
!then the inner "do j" loop does the most work, hence
!the loops over j and k are still the correct way around.
!$OMP DO SCHEDULE(DYNAMIC)
  DO jj = qdims%j_start, qdims%j_end, jblock
    DO j = jj, MIN((jj+jblock)-1,(1+qdims%j_end-qdims%j_start))
      DO i = qdims%i_start, qdims%i_end
        llwic(i,j)=0
      END DO
    END DO

    DO k = levels, 1, -1
      DO j = jj, MIN((jj+jblock)-1,(1+qdims%j_end-qdims%j_start))
        DO i = qdims%i_start, qdims%i_end

          IF (llwic(i,j)  /=  k+1) THEN
            ctt(i,j,k)=t(i,j,k)
          ELSE
            ctt(i,j,k)=ctt(i,j,k+1)
          END IF
          IF (qcf(i,j,k)  >   qcfmin) THEN
            llwic(i,j)=k
          END IF
        END DO
      END DO
    END DO
  END DO
!$OMP END DO

!Cache-blocking applied to loop over k.
!$OMP DO SCHEDULE(DYNAMIC)
  DO kk = levels, 1, -kblock
    DO k = kk, MAX((kk-kblock)+1,1), -1

! ----------------------------------------------------------------------
! 1. Calculate QSAT at liquid/ice water temperature, TL, and initialize
!    cloud water, sub-grid distribution and fraction arrays.
!    This requires a preliminary calculation of the pressure.
!    NB: On entry to the subroutine 'T' is TL and 'Q' is QW.
! ----------------------------------------------------------------------
!
!CDIR COLLAPSE
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          qcl(i,j,k) = 0.0
          cfl(i,j,k) = 0.0
          grid_qc(i,j,k) = rmdi
          bs(i,j,k) = rmdi
        END DO ! i
      END DO ! j

! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsl,t(1,1,k),p_theta_levels(1,1,k),             &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

      DO j = qdims%j_start, qdims%j_end

        DO i = qdims%i_start, qdims%i_end
          IF (multrhc==1) THEN
            rhcritx = rhcrit(i,j,k)
          ELSE
            rhcritx = rhcrit(1,1,k)
          END IF

! Omit CUMULUS points below (and including) NTML+1

          IF ( .NOT. l_param_conv .OR. (l_param_conv .AND.                      &
             (.NOT. cumulus(i,j) .OR. ( cumulus(i,j)                    &
             .AND. (k  >   ntml(i,j)+1) ))) ) THEN

! Rhcrit_if:
            IF (rhcritx  <   1.) THEN
! ----------------------------------------------------------------------
! 2. Calculate the quantity QN = QC/BS = (QW/QSL-1)/(1-RHcrit)
!    if RHcrit is less than 1
! ----------------------------------------------------------------------

              qn(i,j) = (q(i,j,k) / qsl(i,j) - 1.) /                    &
                        (1. - rhcritx)

! ----------------------------------------------------------------------
! 3. Set logical variable for cloud, LQC, for the case RHcrit < 1;
!    where QN > -1, i.e. qW/qSAT(TL,P) > RHcrit, there is cloud
! ----------------------------------------------------------------------

              lqc(i,j) = (qn(i,j)  >   -1.)
            ELSE
! ----------------------------------------------------------------------
! 2.a Calculate QN = QW - QSL if RHcrit equals 1
! ----------------------------------------------------------------------

              qn(i,j) = q(i,j,k) - qsl(i,j)

! ----------------------------------------------------------------------
! 3.a Set logical variable for cloud, LQC, for the case RHcrit = 1;
!     where QN > 0, i.e. qW > qSAT(TL,P), there is cloud
! ----------------------------------------------------------------------

              lqc(i,j) = (qn(i,j)  >   0.)
            END IF ! Rhcrit_if
          ELSE IF (l_param_conv) THEN
            lqc(i,j) = .FALSE.
          END IF  ! Test on CUMULUS and NTML for A05_4A only
        END DO ! i
      END DO ! j

! ----------------------------------------------------------------------
! 4. Form index of points where non-zero liquid cloud fraction
! ----------------------------------------------------------------------

      qc_points=0

      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          IF (lqc(i,j)) THEN
            qc_points = qc_points + 1
            INDEX(qc_points,1) = i
            INDEX(qc_points,2) = j
          END IF
        END DO ! i
      END DO ! j

! ----------------------------------------------------------------------
! 5. Call LS_CLD_C to calculate cloud water content, specific humidity,
!                  water cloud fraction and determine temperature.
! ----------------------------------------------------------------------
! Qc_points_if:
      IF (qc_points  >   0) THEN
! DEPENDS ON: ls_cld_c
        CALL ls_cld_c(p_theta_levels(1,1,k),rhcrit(1,1,k),qsl,qn,       &
        q(1,1,k),t(1,1,k),                                              &
                    qcl(1,1,k),cfl(1,1,k),grid_qc(1,1,k),bs(1,1,k),     &
        INDEX,qc_points,rhc_row_length,rhc_rows,                        &
        bl_levels,k, l_mixing_ratio)
      END IF ! Qc_points_if

! ----------------------------------------------------------------------
! 6. Calculate cloud fractions for ice clouds.
!    THIS IS STILL HIGHLY EXPERIMENTAL.
!    Begin by calculating Qsat_wat(T,P), at Temp. T, for estimate of bs.
! ----------------------------------------------------------------------

! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsl,t(1,1,k),p_theta_levels(1,1,k),             &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

      rootwo = SQRT(2.)

      IF (ice_fraction_method  ==  2) THEN
! Use cloud top temperature and a fixed qsat to give QCFRBS
! DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix(qsl_ctt,ctt(1,1,k),p_theta_levels(1,1,k),     &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

        CALL powr_v(                                                    &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           qsl_ctt,ctt_weight,qsl_ctt )
        CALL powr_v(                                                    &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           qsl,t_weight,qsl )

        subgrid = sub_cld ** (1.0-t_weight)                             &
                  / qsat_fixed ** (1.0-t_weight-ctt_weight)
      END IF ! ice_fraction_method eq 2

      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          IF ( multrhc== 1) THEN
            rhcritx = rhcrit(i,j,k)
          ELSE
            rhcritx = rhcrit(1,1,k)
          END IF
! ----------------------------------------------------------------------
! 6a Calculate qCF/bs.
! ----------------------------------------------------------------------
! Rhcrit_if2:
          IF (rhcritx  <   1.) THEN

            IF (ice_fraction_method  ==  1) THEN
              qcfrbs(i,j)=  qcf(i,j,k) / ((1.-rhcritx) * qsl(i,j))
            ELSE IF (ice_fraction_method  ==  2) THEN
              qcfrbs(i,j) = subgrid * qcf(i,j,k) / ((1.-rhcritx)        &
                       * qsl_ctt(i,j)*qsl(i,j))
            ELSE
! No ice cloud fraction method defined
            END IF ! ice_fraction_method

! ----------------------------------------------------------------------
! 6b Calculate frozen cloud fraction from frozen cloud water content.
! ----------------------------------------------------------------------
            IF (qcfrbs(i,j)  <=  0.) THEN
              cff(i,j,k) = 0.0
            ELSE IF (0.0<qcfrbs(i,j) .AND. (6.0*qcfrbs(i,j)) <= 1.0) THEN
              cff(i,j,k) = 0.5 * ((6. * qcfrbs(i,j))**(2./3.))
            ELSE IF (1.0<(6.*qcfrbs(i,j)) .AND. qcfrbs(i,j) < 1.0) THEN
              phiqcf = ACOS(rootwo * 0.75 * (1. - qcfrbs(i,j)))
              cosqcf = COS((phiqcf + (4. * pi)) / 3.)
              cff(i,j,k) = 1. - (4. * cosqcf * cosqcf)
            ELSE IF (qcfrbs(i,j)  >=  1.) THEN
              cff(i,j,k) = 1.
            END IF
            IF (l_eacf) THEN  ! Empirically adjusted cloud fraction
            ! Back out QN
              IF (0.0< qcfrbs(i,j) .AND. (6.0*qcfrbs(i,j)) <=  1.0) THEN
                qn_imp=SQRT(2.*cff(i,j,k))-1.
              ELSE IF (1.0<(6.0*qcfrbs(i,j)) .AND. qcfrbs(i,j)<1.0) THEN
                qn_imp=1.-SQRT((1.-cff(i,j,k))*2.)
              ELSE
                qn_imp = 1.
              END IF

            ! Modify QN with EACF relationship
              IF (k >  bl_levels) THEN
                qn_adj=(qn_imp+0.0955)/(1.-0.0955)
              ELSE
                qn_adj=(qn_imp+0.184)/(1.-0.184)
              END IF

            ! Recalculate ice cloud fraction with modified QN
              IF (qcfrbs(i,j)  <=  0.) THEN
                cff(i,j,k) = 0.0
              ELSE IF (qn_adj  <=  0.) THEN
                cff(i,j,k) = 0.5 * (1. + qn_adj) * (1. + qn_adj)
              ELSE IF (qn_adj  <   1.) THEN
                cff(i,j,k) = 1. - 0.5 * (1.-qn_adj) * (1.-qn_adj)
              ELSE
                cff(i,j,k) = 1.
              END IF

            END IF  ! L_eacf


          ELSE ! RHcrit = 1, set cloud fraction to 1 or 0

            IF (qcf(i,j,k)  >   0.0) THEN
              cff(i,j,k) = 1.0
            ELSE
              cff(i,j,k) = 0.0
            END IF

          END IF
        END DO ! i
      END DO ! j

! ----------------------------------------------------------------------
! 6c Calculate combined cloud fraction.
! ----------------------------------------------------------------------

      IF (cloud_fraction_method  ==  1) THEN

        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
!             Use minimum overlap condition
            cf(i,j,k) = MIN(cfl(i,j,k)+cff(i,j,k), 1.0)
          END DO
        END DO

      ELSE IF (cloud_fraction_method  ==  2) THEN

          !Hand-unrolled loop
        DO j = qdims%j_start, (1+qdims%j_end-qdims%j_start) -           &
                  MOD((1+qdims%j_end-qdims%j_start),2), 2

          DO i = qdims%i_start, qdims%i_end
! Calculate possible overlaps between ice and liquid in THIS layer
            overlap_max=MIN(cfl(i,j,k),cff(i,j,k))
            overlap_min=MAX(cfl(i,j,k)+cff(i,j,k)-1.0,0.0)
            overlap_random=cfl(i,j,k)*cff(i,j,k)
! Now use the specified degree of overlap to calculate the total
! cloud fraction (= cfice + cfliq - overlap)
            temp0=overlap_random
            temp1=0.5*(overlap_max-overlap_min)
            temp2=0.5*(overlap_max+overlap_min)-overlap_random
            cf(i,j,k)=cfl(i,j,k)+cff(i,j,k)                             &
                    -(temp0+temp1*overlap_ice_liquid                    &
                    +temp2*overlap_ice_liquid*overlap_ice_liquid)
! Check that the overlap wasnt negative
            cf(i,j,k)=MIN(cf(i,j,k),cfl(i,j,k)+cff(i,j,k))

          END DO

          DO i = qdims%i_start, qdims%i_end
! Calculate possible overlaps between ice and liquid in THIS layer
            overlap_max=MIN(cfl(i,j+1,k),cff(i,j+1,k))
            overlap_min=MAX(cfl(i,j+1,k)+cff(i,j+1,k)-1.0,0.0)
            overlap_random=cfl(i,j+1,k)*cff(i,j+1,k)
! Now use the specified degree of overlap to calculate the total
! cloud fraction (= cfice + cfliq - overlap)
            temp0=overlap_random
            temp1=0.5*(overlap_max-overlap_min)
            temp2=0.5*(overlap_max+overlap_min)-overlap_random
            cf(i,j+1,k)=cfl(i,j+1,k)+cff(i,j+1,k)                       &
                    -(temp0+temp1*overlap_ice_liquid                    &
                    +temp2*overlap_ice_liquid*overlap_ice_liquid)
! Check that the overlap wasnt negative
            cf(i,j+1,k)=MIN(cf(i,j+1,k),cfl(i,j+1,k)+cff(i,j+1,k))

          END DO
        END DO

          !Post-conditioning
        IF (MOD((1+qdims%j_end-qdims%j_start),2) == 1) THEN
          DO i = qdims%i_start, qdims%i_end
! Calculate possible overlaps between ice and liquid in THIS layer
            overlap_max=MIN( cfl(i,(1+qdims%j_end-qdims%j_start),k),    &
                             cff(i,(1+qdims%j_end-qdims%j_start),k) )
            overlap_min=MAX( cfl(i,(1+qdims%j_end-qdims%j_start),k) +   &
                             cff(i,(1+qdims%j_end-qdims%j_start),k)-1.0,&
                             0.0 )
            overlap_random=  cfl(i,(1+qdims%j_end-qdims%j_start),k) *   &
                             cff(i,(1+qdims%j_end-qdims%j_start),k)
! Now use the specified degree of overlap to calculate the total
! cloud fraction (= cfice + cfliq - overlap)
            temp0=overlap_random
            temp1=0.5*(overlap_max-overlap_min)
            temp2=0.5*(overlap_max+overlap_min)-overlap_random
            cf(i,(1+qdims%j_end-qdims%j_start),k)=                      &
              cfl(i,(1+qdims%j_end-qdims%j_start),k) +                  &
              cff(i,(1+qdims%j_end-qdims%j_start),k) -                  &
              (temp0+temp1*overlap_ice_liquid +                         &
              temp2*overlap_ice_liquid*overlap_ice_liquid)
! Check that the overlap wasnt negative
            cf(i,(1+qdims%j_end-qdims%j_start),k)=                      &
              MIN(cf(i,(1+qdims%j_end-qdims%j_start),k),                &
                 cfl(i,(1+qdims%j_end-qdims%j_start),k) +               &
                 cff(i,(1+qdims%j_end-qdims%j_start),k))

          END DO
        END IF !Post-conditioning

! CFF + CFL >= 1 implies CF = 1.
! Deal with this case separately to avoid roundoff issues:
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
            IF (cfl(i,j,k)+cff(i,j,k) >= 1.0) cf(i,j,k) = 1.0
          END DO
        END DO

      ELSE
! No total cloud fraction method defined


      END IF ! cloud_fraction_method

    END DO
  END DO ! Levels_do
!$OMP END DO

!$OMP END PARALLEL

9999 CONTINUE ! Error exit
  IF (lhook) CALL dr_hook('LS_CLD',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE ls_cld
