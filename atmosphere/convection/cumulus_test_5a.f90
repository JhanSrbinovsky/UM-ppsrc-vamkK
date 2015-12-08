! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Various tests to determine whether ascent is convective
!
MODULE cumulus_test_mod

IMPLICIT NONE

CONTAINS

! Subroutine Interface:
SUBROUTINE cumulus_test(npnts,mbl,wet_model_levels,                      &
                        ntml, k_plume, kcucheck,                         &
                        zhpar, land_frac, qw, cloud_fraction, z_theta,   &
                        z_lcl, nlcl, cumulus )

! ------------------------------------------------------------------------------
! Description:
!   Works out which unstable points are convective using top of parcel.
!   If parcel top is < 3000m performs a series of tests to check whether
!   the point is convective or stratocumulus. 
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
USE cv_diag_param_mod, ONLY:                                             &
    sc_cftol

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine arguments
 
INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,mbl                  & ! maximum number of boundary levels bl_levels-1
 ,wet_model_levels       ! Number of wet model levels 

INTEGER, INTENT(IN) :: &
  ntml(npnts)          & ! Top of parcel ascent
 ,k_plume(npnts)       & ! Starting model level for plume ascent
 ,kcucheck(npnts)        ! Starting model level for plume ascent

REAL, INTENT(IN) ::                      &
  zhpar(npnts)                           & ! height of parcel top (m)
 ,land_frac(npnts)                       & ! fraction of land
 ,qw(npnts,wet_model_levels)             & ! total water (kg/kg)
 ,cloud_fraction(npnts,wet_model_levels) & ! layer cloud fraction
 ,z_theta(npnts,wet_model_levels)        & ! height of model theta levels(m)      
 ,z_lcl(npnts)                             ! height of LCL (m)

INTEGER, INTENT(INOUT) :: &
  nlcl(npnts)                 ! Lifting condensation level

LOGICAL, INTENT(INOUT) :: & 
  cumulus(npnts)              ! true if convective point  


!-----------------------------------------------------------------------
! Local variables

INTEGER ::        & 
  ii,k              ! loop counter

REAL ::           &
  grad_cld        & ! SVL gradient in layer above LCL.
 ,grad_sub          ! SVL gradient in layer below LCL. 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CUMULUS_TEST',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------

DO ii=1, npnts

!-----------------------------------------------------------------------
!     Check lifting condensation level against height of parcel ascent.
!     If lifting condensation level lower than parcel ascent, and is
!     within bl_levels, then decide on type of cloudy layer.
!     Gradient tests and cumulus parametriztion require a minimum number 
!     of grid-levels between LCL and top of parcel ascent, otherwise 
!     define as stratocumulus.  To avoid resolution sensitivity, also 
!     require this depth to be more than 400m, i.e. cumulus clouds are 
!     deeper than 400m.  Note this depth is physically plausible and 
!     close to the two grid-level requirement in the original L38 
!     implementation.
!-----------------------------------------------------------------------
! Note the following test still uses MBL making it possibly resolution 
! dependent. This test is only likely to change results if the number of 
! boundary layer levels is reduced. This is unlikely to be done. 
!-----------------------------------------------------------------------

  IF ( ntml(ii)-nlcl(ii) >= 3 .AND. nlcl(ii)  <   MBL-1                 &
                              .AND. zhpar(ii)-z_lcl(ii) >= 400. ) THEN 

!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = zh
!     For cumulus top of mixed layer = Z_LCL
!     New test zhpar >= 3000m replaces (NTML >= MBL) i.e. parcel top 
!     above 3000m indicates convection.
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------

    IF (zhpar(ii) >= 3000.0) THEN
      
      cumulus(ii) = .TRUE.

      ! Added for low LCL levels 
      ! nlcl is not permitted to be less than 2 if Cu is diagnosed
      nlcl(ii) = MAX(2, nlcl(ii))


    ! Tops < 3000m and nlcl > kplume
    ELSE IF ( nlcl(ii)  >   k_plume(ii)) THEN   

! Current test is against a height of ~<2.5km
! This could be replaced by a scale height if a suitable method
! for determining a sensible height was possible from profile/cumulus
! depth information available in this routine

      IF (ntml(ii)  >   kcucheck(ii)                                    &
                  .AND. nlcl(ii)  <=  kcucheck(ii)-2) THEN

        grad_cld = ABS( QW(ii,kcucheck(ii)) - QW(ii,nlcl(ii)) )         &
                    /( z_theta(ii,kcucheck(ii)) - z_theta(ii,nlcl(ii)) )
      ELSE
        grad_cld = ABS( QW(ii,ntml(ii))     - QW(ii,nlcl(ii)) )         &
                    /( z_theta(ii,ntml(ii)) - z_theta(ii,nlcl(ii)) )
      END IF

      grad_sub   =  ABS( QW(ii,nlcl(ii)) - QW(ii,k_plume(ii)) )         &
                    /( z_theta(ii,nlcl(ii)) - z_theta(ii,k_plume(ii)) )

      IF (grad_cld  >   1.10*grad_sub) THEN
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level DOwn.
!     Note typical cumulus profiles are expected to have a fairly
!     uniform q profile from the surface to the cloud base and THEN a
!     decreasing profile of q above this in the cloud. Typical the
!     decreasing gradient from the cloud base to 2.5km will be the
!     order of > 1.10 the below cloud value.
!-----------------------------------------------------------------------

! test against a height ~ 2.5km

        IF (ntml(ii)  <=  kcucheck(ii)) THEN
                grad_cld = ABS( QW(ii,ntml(ii)-1) - QW(ii,nlcl(ii)) )   &
                      /( z_theta(ii,ntml(ii)-1) - z_theta(ii,nlcl(ii)) )
        END IF

        IF ( grad_cld  >   1.10*grad_sub) THEN
!-----------------------------------------------------------------------
!      Diagnose a cumulus layer
!-----------------------------------------------------------------------
          cumulus(ii) = .TRUE.
        END IF

      ELSE

! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
! and not yet been mixed (so could have been erroneously identIfied as
! well-mixed)

! First check using level below (recalculate grad_sub)

        IF (nlcl(ii) - k_plume(ii)  >=  2) THEN

          grad_sub = ABS( QW(ii,nlcl(ii)-1) - QW(ii,k_plume(ii)) )           &
                     /( z_theta(ii,nlcl(ii)-1) - z_theta(ii,k_plume(ii)) )

          IF ( grad_cld  >   1.10*grad_sub) THEN
            cumulus(ii) =.TRUE.
          END IF

        END IF

! If still diagnosing well-mixed, check using level above
! (recalculate grad_cld)

        IF (.NOT. cumulus(ii) ) THEN

          IF (ntml(ii)  >   kcucheck(ii)                                     &
                   .AND. nlcl(ii)  <=  kcucheck(ii)-2) THEN

            grad_cld = ABS( QW(ii,kcucheck(ii)) - QW(ii,nlcl(ii)+1) )        &
                       /( z_theta(ii,kcucheck(ii)) - z_theta(ii,nlcl(ii)+1) )
          ELSE
            grad_cld = ABS( QW(ii,ntml(ii)) - QW(ii,nlcl(ii)+1) )            &
                       /( z_theta(ii,ntml(ii)) - z_theta(ii,nlcl(ii)+1) )
          END IF

          IF ( grad_cld  >   1.10*grad_sub) THEN
            cumulus(ii) =.TRUE.
          END IF

        END IF   ! not cumulus
      END IF     ! test on cloud gradient
    END IF       ! test on cloud top height
  END IF         ! tests on nlcl

END DO           ! ii loop 

!-----------------------------------------------------------------------
!      Check that a cumulus layer has not been erroneously diagnosed in
!      a deep cloudy region
!      As the above checks are DOne on the total water rather than q it
!      is possible the conditions can be met in areas where the level of
!      prognostic qcl or qcf is high. The type of mistake is only
!      thought to occur over land. The wording of this comment implies
!      the problem occurred for situation where the gradient tests were 
!      done i.e. zhpar < 3000m.       
!-----------------------------------------------------------------------
! The following code still uses a MBL test. I think this test was designed
! with MBL ~ 3000m in mind and so should really be replaced with
!   zhpar(ii) < 3000.0 instead of ntml(ii)  <   MBL
! Changing this may alter results so I am not doing it yet. 
!-----------------------------------------------------------------------

DO ii=1, npnts
  K=nlcl(ii)

  IF ( (Land_frac(ii) > 0.0) .AND. cumulus(ii) .AND.              &
                                       ntml(ii)  <   MBL ) THEN
    DO WHILE( K  <=  ntml(ii) .AND. cloud_fraction(ii,K)          &
                                                  >=  sc_cftol )
      K = K + 1
    END DO
    IF (K  ==  ntml(ii)+1) cumulus(ii) = .FALSE.
  END IF
END DO       ! ii loop 

IF (lhook) CALL dr_hook('CUMULUS_TEST',zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
END SUBROUTINE cumulus_test

END MODULE cumulus_test_mod
