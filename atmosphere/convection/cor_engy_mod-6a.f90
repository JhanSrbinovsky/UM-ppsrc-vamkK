! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Ensure conservation of moist static energy

MODULE cor_engy_6a_mod

IMPLICIT NONE

! ------------------------------------------------------------------------------
!
! Description:
! Adjust the potential temperature increments to ensure the conservation
! of moist static energy.
!
! Method:
!  1. The current mass flux convection scheme was originally coded for a 
!    model using a vertical coordinate of hybrid pressure.
!    Column integrals were done by summing  using dp/g
!    At present a sum using this method best represents what the scheme 
!    will try to conserve.
!  2. The UM grid should now really be using summation using rho*r/a*r/a*dr
!     where a is the radius of the earth. There are problems with this as
!     the dynamics conserves moisture after interpolating it to the rho grid
!     in the vertical whereas the physics works on q, qcl, qcf on the theta
!     grid. Checking conservation of the theta grid will give a different 
!     answer to checking conservation on the rho grid. Both will be different
!     to using the pressure integrals. 
!
!  See UM Documentation paper No 27
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

CONTAINS

SUBROUTINE cor_engy_6a(np_field,npnts,nlev,index1, r2rho_th                  &
                    ,dr_across_th,exner_layer_centres,p_layer_boundaries     &
                    ,dqbydt,dqclbydt,dqcfbydt,rain,snow                      &
                    ,dthbydt)

USE cv_derived_constants_mod, ONLY: ra2, cv
USE earth_constants_mod, ONLY: g
USE water_constants_mod, ONLY: lc, lf
USE atmos_constants_mod, ONLY: cp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  np_field             & ! Length of data 
 ,npnts                & ! Number of convecting columns
 ,nlev                   ! Number of model layers


INTEGER, INTENT(IN) :: &
  index1(np_field)          ! index of points with convection

REAL,INTENT(IN) ::                   &
  r2rho_th(np_field,nlev)            & ! radius**2 density for theta lev (kg/m)
 ,dr_across_th(np_field,nlev)        & ! thickness of theta levels (m)
 ,exner_layer_centres(np_field,0:nlev) & ! Exner ratio
 ,p_layer_boundaries(np_field,0:nlev)    ! layer boundary (Pa)

REAL,INTENT(IN) :: dqclbydt(np_field,nlev) ! Increment to model cloud water
                                           ! (kg/kg/s)
REAL,INTENT(IN) :: dqcfbydt(np_field,nlev) ! Increment to model cloud ice
                                           ! (kg/kg/s)

REAL,INTENT(IN) ::      &
  rain(np_field)        & ! rain at surface (kg/m**2/s)
 ,snow(np_field)          ! Snow at surface (kg/m**2/s)


!----------------------------------------------------------------------
! Variables which are input but which are also updated in this routine
!----------------------------------------------------------------------

REAL, INTENT(INOUT) ::   &  
  dthbydt(np_field,nlev)   ! IN increments to model potential temperature
                           !    dur to convection (K/s)
                           ! OUT corrected increments to model potential 
                           !    temperature due to convection (K/s) 
REAL,INTENT(INOUT) :: dqbydt(np_field,nlev)  ! Increment to water vapour
                                             ! (kg/kg/s)

!---------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER ::        & 
  i,j,k             ! loop counters 

LOGICAL ::        &
  bposer(npnts)   & ! Mask for points in layer k at which increments to model
                    ! potential temperature due to convection are positive
 ,bcorr(npnts)      ! Mask for points at which enthalpy correction in
                    ! necessary

REAL :: qsum(npnts)    ! Mass weighted vertical summation of convective 
                       ! increments to humidity (kg/m2/s)
REAL :: qclsum(npnts)  ! Mass weighted vertical summation of convective 
                       ! increments to qcl (kg/m2/s)
REAL :: qcfsum(npnts)  ! Mass weighted vertical summation of convective 
                       ! increments to qcl (kg/m2/s)
REAL :: qspos(npnts)   ! Mass weighted vertical summation of convective 
                       ! positive increments to humidity (kg/m2/s)
REAL :: qsneg(npnts)   ! Mass weighted vertical summation of convective 
                       ! negative increments to humidity (kg/m2/s)
REAL :: qerror(npnts)  ! error in moisture conservation (kg/m2/s)
REAL :: qprecip(npnts) ! moisture lost from atmospheric column (kg/m2/s)
REAL :: dqscale(npnts) ! scaling factor used to correct the humidity
                       ! increments
REAL :: tspos(npnts)   ! Mass weighted vertical summation of convective 
                       ! positive increments to temperature (W/m2)
REAL :: tsneg(npnts)   ! Mass weighted vertical summation of convective 
                       ! negative increments to temperature (W/m2)
REAL :: terr(npnts)    ! The error in the vertical summation of energy 
                       ! (W/m2)
REAL :: dtscale(npnts) ! scaling factor used to correct the temperature
                       ! increments

! The var1 variables can be be used to check the conservation on the 
! pressure level grid. Useful for development but should be removed 
! when routine is actively used.
REAL :: qsum1(npnts)   ! Mass weighted vertical summation of convective 
                       ! increments to humidity (kg/m2/s)
REAL :: qclsum1(npnts) ! Mass weighted vertical summation of convective 
                       ! increments to qcl (kg/m2/s)
REAL :: qcfsum1(npnts) ! Mass weighted vertical summation of convective 
                       ! increments to qcl (kg/m2/s)
REAL :: qspos1(npnts)  ! Mass weighted vertical summation of convective 
                       ! positive increments to humidity (kg/m2/s)
REAL :: qsneg1(npnts)  ! Mass weighted vertical summation of convective 
                       ! negative increments to humidity (kg/m2/s)
REAL :: qerror1(npnts) ! error in moisture conservation (kg/m2/s)
REAL :: qprecip1(npnts)! moisture lost from atmospheric column (kg/m2/s)
REAL :: dqscale1(npnts)! scaling factor used to correct the humidity
                       ! increments
REAL :: tspos1(npnts)  ! Mass weighted vertical summation of convective 
                       ! positive increments to temperature (W/m2)
REAL :: tsneg1(npnts)  ! Mass weighted vertical summation of convective 
                       ! negative increments to temperature (W/m2)
REAL :: terr1(npnts)   ! The error in the vertical summation of energy 
                       ! (W/m2)
REAL :: dtscale1(npnts)! scaling factor used to correct the temperature
                       ! increments



REAL :: r2rhodr           ! (r/a)**2 rho dr 
REAL :: extempk           ! Exner ratio at the mid-point of layer k

REAL :: delpk             ! Difference in pressure across a layer (Pa)

INTEGER :: corr_moisture  ! Correct moisture conservation 
INTEGER :: corr_option    ! Option for the type of correction
                          ! 0: one-sided correction
                          ! 1: two-sided correction

INTEGER, PARAMETER :: onesidecorr = 0 ! Correct by apply a scaling to either 
                                      ! the positive or negative increments
INTEGER, PARAMETER :: bothsidecorr= 1 ! Correct by apply a scaling all the 
                                      ! increments
REAL, PARAMETER :: SmallNum = TINY(terr) !Small number


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('COR_ENGY_6A',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!  Check moisture conservation - not being used at present.
!  Assuming PC2 being used so increments to qcl & qcf.
!
!   total surface precip = sum { (dq+dqcl+qcf)*rho*(r/a)*(r/a)*dr}
! Or
!   total surface precip = sum { (dq+dqcl+qcf)*dp/g}
! ----------------------------------------------------------------------
corr_moisture = 1
corr_option   = 1
IF (corr_moisture == 1) THEN

  ! Initialise moisture variables.
  DO j=1,npnts
    qclsum(j)  = 0.0
    qcfsum(j)  = 0.0
    qspos(j)   = 0.0
    qsneg(j)   = 0.0
    dqscale(j) = 0.0
    qprecip(j) = rain(index1(j)) + snow(index1(j)) 
    bcorr(j)   = .FALSE.

    qclsum1(j)  = 0.0
    qcfsum1(j)  = 0.0
    qspos1(j)   = 0.0
    qsneg1(j)   = 0.0
    dqscale1(j) = 0.0
  END DO

  ! integrate change in column moisture
  DO k= 1,nlev
    DO j=1,npnts
      i = index1(j) 
  
      r2rhodr   = r2rho_th(i,k)*dr_across_th(i,k)*ra2
      qclsum(j) = qclsum(j) + dqclbydt(i,k)*r2rhodr  
      qcfsum(j) = qcfsum(j) + dqcfbydt(i,k)*r2rhodr  


!     Useful for check water conservation on pressure level grid
      delpk     = - (p_layer_boundaries(i,k) - p_layer_boundaries(i,k-1))
      qclsum1(j)= qclsum1(j) + dqclbydt(i,k)*delpk/g
      qcfsum1(j)= qcfsum1(j) + dqcfbydt(i,k)*delpk/g

      IF (dqbydt(i,k)  >   0.0) THEN
        qspos(j)  = qspos(j)  + dqbydt(i,k)*r2rhodr
        qspos1(j) = qspos1(j) + dqbydt(i,k)*delpk/g
      ELSE
        qsneg(j)  = qsneg(j)  + dqbydt(i,k)*r2rhodr
        qsneg1(j) = qsneg1(j) + dqbydt(i,k)*delpk/g
      END IF
    END DO
  END DO

  ! Calculate the absolute error
  DO j=1,npnts
    qerror(j)  = qprecip(j) + qspos(j)  + qsneg(j)  + qclsum(j)  + qcfsum(j)
    qerror1(j) = qprecip(j) + qspos1(j) + qsneg1(j) + qclsum1(j) + qcfsum1(j)
  END DO

  IF (corr_option == onesidecorr) THEN
    !Calculate the correction scaling
    DO j=1,npnts
      bposer(j) = qerror(j)  >   0.0
    
      IF (bposer(j) .AND. (qspos(j)  ==  0.0)) THEN
        bposer(j) = .FALSE.
      ELSE IF (.NOT.bposer(j) .AND. (qsneg(j)  ==  0.0)) THEN
        bposer(j) = .TRUE.
      END IF
    
      bcorr(j) = (ABS(qspos(j)) > SmallNum) &
            .OR. (ABS(qsneg(j)) > SmallNum)
    
      IF (bposer(j) .AND. bcorr(j)) THEN
        dqscale(j) = 1.0 - qerror(j)/qspos(j)
      ELSE IF (.NOT.bposer(j) .AND. bcorr(j)) THEN
        dqscale(j) = 1.0 - qerror(j)/qsneg(j)
      ELSE
        dqscale(j) = 1.0
      END IF
    END DO

      
    !Apply the correction scaling
    DO k=1,nlev
      DO j=1,npnts
        IF (bcorr(j) .AND. (( bposer(j) .AND. (dqbydt(index1(j),k) > 0.0))     &
            .OR. ( .NOT.bposer(j) .AND. (dqbydt(index1(j),k) < 0.0)))) THEN
           dqbydt(index1(j),k) = dqbydt(index1(j),k)*dqscale(j)
        END IF
      END DO  ! npnts
    END DO ! nlev

  ELSE IF (corr_option == bothsidecorr) THEN

    !Calculate the correction scaling
    DO j=1,npnts

      bcorr(j) = ABS(qerror(j)) > 10.*ABS(EPSILON(qspos)*(qspos(j)+qsneg(j)))
    
      IF (bcorr(j)) THEN
        dqscale(j) = 1.0 - qerror(j)/(qspos(j)+qsneg(j))
      ELSE
        dqscale(j) = 1.0
      END IF
    END DO

    DO j=1,npnts

      IF (ABS(qerror1(j)) > 10.*ABS(EPSILON(qspos1)*(qspos1(j)+qsneg1(j)))    &
          .AND. ((qspos1(j)+qsneg1(j))) > SmallNum) THEN
        dqscale1(j) = 1.0 - qerror1(j)/(qspos1(j)+qsneg1(j))
      ELSE
        dqscale1(j) = 1.0
      END IF
    END DO
      
    !Apply the correction scaling
    DO k=1,nlev
      DO j=1,npnts
        IF (bcorr(j)) THEN
           dqbydt(index1(j),k) = dqbydt(index1(j),k)*dqscale(j)
        END IF
      END DO  ! npnts
    END DO ! nlev
  END IF
END IF

! ----------------------------------------------------------------------
! Sum up mixing ratio and  +ve and -ve temperature increments
! ----------------------------------------------------------------------
! Moist static energy H = cpT + Lq 
! or                  H = cvT + gz + Lq      
!
!  T = th * exner
!
! For conservation
!   dH = Lf*snow = sum{ (cp dth*exner + Lc dq) * dp/g}
! Or
!   dH = Lf*snow = sum{ (cv dth*exner + Lc dq) *rho*(r/a)*(r/a)*dr} 
!
! Is the energy correct for PC2?
! Or should I check
!  dH = -lf*rain = sum{ (cv dth*exner + (Lc+lf)*dq+lf*dqcl)*rho*(r/a)*(r/a)*dr}
! Checking this makes the energy error more negative in most cases.
!
! Note prior to PC2 conservation of moisture would imply
!
!    rain+snow = sum{ (dq dp/g)}
!
! ----------------------------------------------------------------------

DO j=1,npnts
  qsum(j)     = 0.0
  qclsum(j)   = 0.0
  tspos(j)    = 0.0
  tsneg(j)    = 0.0
  dtscale(j)  = 0.0
  bcorr(j)   = .FALSE.
END DO

DO k=1,nlev
  DO j=1,npnts
    i = index1(j)

    r2rhodr   = r2rho_th(i,k)*dr_across_th(i,k)*ra2
    qsum(j)   = qsum(j)   + dqbydt(i,k)*r2rhodr  
    qclsum(j) = qclsum(j) + dqclbydt(i,k)*r2rhodr  

    IF (dthbydt(i,k)  >   0.0) THEN
      tspos(j)  = tspos(j)  + cp*dthbydt(i,k)*                              &
                                exner_layer_centres(i,k)*r2rhodr 
    ELSE
      tsneg(j)  = tsneg(j)  + cp*dthbydt(i,k)*                              &
                                exner_layer_centres(i,k)*r2rhodr 
    END IF

  END DO
END DO

! ----------------------------------------------------------------------
!   Calculate the error and apply the necessary correction
!
!   UM DOCUMENTATION PAPER 27
!
! Note method relies on conservation of moisture. This was approximately
! "true" if using the old pressure integrals but not for rho*r*r*dr.
! Given that moisture is not being conserved and the errors in this
! appear to be large in some cases then it is not surprising that the energy
! conservation by the scheme is very poor. Given the large errors
! it is not a very good idea to try to correct them as the calculations
! imply very large corrections are required in some columns.
! Until the scheme is altered to correctly take account of the geometry
! of the grid then it does not make sense to try to apply an energy 
! correction.
! ----------------------------------------------------------------------

! Calculate the absolute energy error
DO j=1,npnts
  terr(j) = (lc+lf)*qsum(j) +lf*qclsum(j) + lf*rain(index1(j))              &
                          + tspos(j) + tsneg(j)
END DO


IF (corr_option == onesidecorr) THEN
  !Calculate the correction scaling
  DO j=1,npnts
    bposer(j) = terr(j)  >   0.0
  
    IF (bposer(j) .AND. (tspos(j)  ==  0.0)) THEN
      bposer(j) = .FALSE.
    ELSE IF (.NOT.bposer(j) .AND. (tsneg(j)  ==  0.0)) THEN
      bposer(j) = .TRUE.
    END IF
  
    bcorr(j) = (ABS(tspos(j)) > SmallNum) .OR. (ABS(tsneg(j)) > SmallNum)
  
    IF (bposer(j) .AND. bcorr(j)) THEN
      dtscale(j) = 1.0 - terr(j)/tspos(j)
    ELSE IF (.NOT.bposer(j) .AND. bcorr(j)) THEN
      dtscale(j) = 1.0 - terr(j)/tsneg(j)
    ELSE
      dtscale(j) = 1.0
    END IF
  
  END DO

  
  !Apply the correction scaling
  DO k=1,nlev
    DO j=1,npnts
      IF (bcorr(j) .AND. (( bposer(j) .AND. (dthbydt(index1(j),k) > 0.0))     &
          .OR. ( .NOT.bposer(j) .AND. (dthbydt(index1(j),k) < 0.0)))) THEN
         dthbydt(index1(j),k) = dthbydt(index1(j),k)*dtscale(j)
      END IF
    END DO  ! npnts
  END DO ! nlev

ELSE IF (corr_option == bothsidecorr) THEN

  !Calculate the correction scaling
  DO j=1,npnts

    bcorr(j) = ABS(tspos(j)+tsneg(j)) > SmallNum
  
    IF (bcorr(j)) THEN
      dtscale(j) = 1.0 - terr(j)/(tspos(j)+tsneg(j))
    ELSE
      dtscale(j) = 1.0
    END IF
  
  END DO
  
  !Apply the correction scaling
  DO k=1,nlev
    DO j=1,npnts
      IF (bcorr(j)) THEN
         dthbydt(index1(j),k) = dthbydt(index1(j),k)*dtscale(j)
      END IF
    END DO  ! npnts
  END DO ! nlev

END IF

IF (lhook) CALL dr_hook('COR_ENGY_6A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE cor_engy_6a

END MODULE cor_engy_6a_mod
