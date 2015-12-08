! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Ensure conservation of moist static energy
!
SUBROUTINE cor_engy_5a(np_field,npnts,nlev,index1,r_theta, r_rho             &
                    ,r2rho_th, r2rho, dr_across_th, dr_across_rh             &
                    ,exner_layer_centres, theta, u, v                        &
                    ,dubydt, dvbydt, dqclbydt, dqcfbydt                      &
                    ,rain,snow,dqbydt,dthbydt)

USE cv_derived_constants_mod, ONLY: ra2, cv
USE cv_run_mod,               ONLY: l_mom
USE earth_constants_mod, ONLY: g, earth_radius
USE water_constants_mod, ONLY: lc, lf, tm
USE atmos_constants_mod, ONLY: cp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! NOTE this code is not currently recommended for use.
! It is still under development.
! Description:
! 1. Check moisture conservation and correct this by adjusting
!    the rain or snow if possible. If not, scale the negative dq increments.
!    Moisture is now conserved.
! 2. Evaluate the error in energy conservation by the convection scheme now
!    that the moisture is conserved.
! 3. Adjust the increments to potential temperature to give conservation.
! 4. Any columns where the required adjustment seems unreasonable will 
!    generate printout. Currently the positive/negative increments are only 
!    scaled by a factor between 0.4 and 1.6.
!
!  See UM Documentation paper No 27
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
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
  r_theta(np_field,0:nlev)           & ! radius of theta levels (m)
 ,r_rho(np_field,nlev)               & ! radius of rho levels (m)
 ,r2rho_th(np_field,nlev)            & ! radius**2 density for theta lev (kg/m)
 ,r2rho(np_field,nlev)               & ! radius**2 density for rho lev (kg/m)
 ,dr_across_th(np_field,nlev)        & ! thickness of theta levels (m)
 ,dr_across_rh(np_field,nlev)        & ! thickness of rho levels (m)
 ,exner_layer_centres(np_field,0:nlev) & ! Exner ratio
 ,theta(np_field,nlev)                 & ! theta level (K)
 ,u(np_field,nlev)                     & ! U wind (m/s)
 ,v(np_field,nlev)                       ! V wind (m/s)

REAL,INTENT(IN) ::            &
  dubydt(np_field,nlev)       & ! Increment to model u wind (m/s)
 ,dvbydt(np_field,nlev)       & ! Increment to model v wind (m/s)
 ,dqclbydt(np_field,nlev)     & ! Increment to model cloud water (kg/kg)
 ,dqcfbydt(np_field,nlev)       ! Increment to model cloud ice (kg/kg)

!----------------------------------------------------------------------
! Variables which are input but which are also updated in this routine
!----------------------------------------------------------------------
REAL,INTENT(INOUT) ::   &
  rain(np_field)        & ! rain at surface (kg/m**2/s)
 ,snow(np_field)          ! Snow at surface (kg/m**2/s)

REAL, INTENT(INOUT) ::   &  
  dqbydt(np_field,nlev)  & ! IN Increment to model water vapour (kg/kg)
                           ! OUT corrected increments 
 ,dthbydt(np_field,nlev)   ! IN increments to model potential temperature
                           !    dur to convection (K/s)
                           ! OUT corrected increments to model potential 
                           !    temperature due to convection (K/s) 
!-------------------------------------------------------------------------------
! May become output arrays in future
!-------------------------------------------------------------------------------
REAL ::                  &
  energy_err(np_field)   & ! Error in energy conservation (W/m2)
 ,ke_cmt_loss(np_field)    ! KE lost due to CMT (W/m2)

!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER ::        & 
  i,j,k,ip,ineg,ii  ! loop counters 

INTEGER ::        & 
  ienergy_method    ! Method for integrating energy 
                    ! 1 Use data on theta levels
                    ! 2 Use data on rho levels

INTEGER ::        & 
  indexn(npnts)   & ! index for negative points
 ,indexp(npnts)     ! index for positive points


LOGICAL ::        &
  bposer(npnts)   & ! .True.  - scale positive dtheta/dt increments
                    ! .false. - scale negative dtheta/dt increments
 ,bcorr(npnts)      ! Mask for points at which energy correction is
                    ! necessary

REAL ::         &
  qsum(npnts)   & ! Summation of increments to model mixing ratio due to
                  ! convection in the vertical, weighted according to the
                  ! mass of the layer (kg/m**2/s)
 ,qclsum(npnts) & ! Summation of increments to model qcl due to
                  ! convection in the vertical, weighted according to the
                  ! mass of the layer (kg/m**2/s)
 ,tspos(npnts)  & ! Summation of positive increments to model potential
                  ! temperature due to convection with height,weighted 
                  ! according to the mass of the layer (kg/m**2/s)
 ,tsneg(npnts)  & ! Summation of negative increments to model potential
                  ! temperature due to convection with height,weighted 
                  ! according to the mass of the layer (kg/m**2/s)
 ,qpos(npnts)   & ! Summation of positive increments to model q
                  ! due to convection with height,weighted 
                  ! according to the mass of the layer (kg/m**2/s)
 ,qneg(npnts)   & ! Summation of negative increments to model q
                  ! temperature due to convection with height,weighted 
                  ! according to the mass of the layer (kg/m**2/s)
 ,scalee(npnts) & ! Scaling factor used to correct dtheta/dt increments
 ,terr(npnts)     ! Summation of all increments to model potential
                  ! temperature due to convection with height,weighted 
                  ! according to the mass of the layer (kg/m**2/s)

REAL ::                     &
  dtbydt(npnts,nlev)        & ! dT
 ,dtbydt_rho(npnts,nlev)    & ! dT on rho levels
 ,dqbydt_rho(npnts,nlev)    & ! dq on rho levels
 ,dqclbydt_rho(npnts,nlev)  & ! dqcl on rho levels
 ,dqcfbydt_rho(npnts,nlev)    ! dqcf on rho levels

REAL ::            & 
  factorn(npnts)   & ! factor for negative points
 ,factorp(npnts)     ! factor for positive points

REAL ::              &
  qloss(np_field)    & ! moisture change in atmospheric column kg/kg/step
 ,qprecip(np_field)  & ! moisture lost from atmospheric column kg/kg/step
 ,qerror(np_field)     ! error in moisture conservation

REAL ::              &
  KE_loss(np_field)    ! Loss in KE due to convection J/m2/s

REAL ::              &
  r2rhodr              ! (r/a)**2 rho dr 

REAL ::              &
  t1                   ! level 1 T

REAL ::              &
  weight1, weight2, weight3, ww1,ww2   ! weights for interpolation

INTEGER ::         &
  icheck_moisture  & ! 1 Check moisture conservation 
 ,icorrect         & ! 1 apply an energy correction
 ,ike_cal          & ! 1 calculate KE change due to convection
 ,nfalse             ! Number of problem points
 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('COR_ENGY_5A',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! The energy integrals can be calculated using;
!
! Method 1 - tested and works
!  T and q on their original grid. This is consistent with the convection
! scheme. 
!
! Method 2
! Interpolating the increments to T and q to the rho levels. This is 
! consistent with the overall model energy correction but inconsistent
! with the convection scheme. (This version currently blewup in testing
! so is not available for use.)
! ----------------------------------------------------------------------
! Set up arrays
! output diagnostics
DO i=1,np_field
  energy_err(i)  = 0.0
  ke_cmt_loss(i) = 0.0
END DO

DO j=1,npnts
  qloss(j) = 0.0
  qpos(j)  = 0.0
  qneg(j)  = 0.0
  i= index1(j) 
  qprecip(j) = rain(i) + snow(i)  

  qsum(j)    = 0.0
  qclsum (j) = 0.0
  tspos(j) = 0.0
  tsneg(j) = 0.0
END DO

ienergy_method = 1  ! Only safe method at present

IF (ienergy_method == 1) THEN 

  ! ----------------------------------------------------------------------
  !  Check moisture conservation and correct error.
  !  Assuming PC2 being used so include increments to qcl & qcf.
  !
  !   total surface precip = sum { (dq+dqcl+qcf)*rho*(r/a)*(r/a)*dr}
  ! ----------------------------------------------------------------------
  ! integrate change in column moisture


  DO k= 1,nlev
    DO j=1,npnts
      i= index1(j) 
      qloss(j) = qloss(j)+ (dqbydt(i,k)+dqclbydt(i,k)+dqcfbydt(i,k))      &
                                 *r2rho_th(i,k)*dr_across_th(i,k)*ra2
      IF (dqbydt(i,k) > 0.0) THEN
        qpos(j) = qpos(j) + dqbydt(i,k)*r2rho_th(i,k)*dr_across_th(i,k)*ra2 
      END IF
      IF (dqbydt(i,k) < 0.0) THEN
        qneg(j) = qneg(j) + dqbydt(i,k)*r2rho_th(i,k)*dr_across_th(i,k)*ra2 
      END IF
    END DO
  END DO

  DO i=1,npnts
    qerror(i) = qprecip(i) + qloss(i)
  END DO


  ! Attempt to correct moisture budget
  ! Is there precipitation - if so attempt to increase /reduce to match
  ! increments?
  ineg=0
  ip=0
  DO i=1,npnts
    IF (qerror(i) < 0.0) THEN
    ! should be more precip, therefore increase precip
      IF (qprecip(i) == 0.0) THEN
      ! Need to decide whether precip is rain or snow
        t1 = theta(index1(i),1)*exner_layer_centres(index1(i),1)  
        IF(t1 < tm) THEN
          snow(index1(i)) = snow(index1(i)) - qerror(i)
        ELSE
          rain(index1(i)) = rain(index1(i)) - qerror(i)
        END IF 
      ELSE
        IF (rain(index1(i)) > 0.0) THEN
          rain(index1(i)) = rain(index1(i)) - qerror(i)
        ELSE
          snow(index1(i)) = snow(index1(i)) - qerror(i)
        END IF
      END IF
    ELSE  ! Error is positive
    ! should be less precip
      IF (qprecip(i) >= qerror(i)) THEN
        IF (rain(index1(i)) > 0.0) THEN
          rain(index1(i)) = rain(index1(i)) - qerror(i)
        ELSE
          snow(index1(i)) = snow(index1(i)) - qerror(i)
        END IF
      ELSE
        ! Existing precip less than required reduction  
        ! Zero precip if any and then decide what else to do?
        IF (qprecip(i) > 0.0 )THEN
          qerror(i) = qerror(i) - qprecip(i)
          IF (snow(index1(i)) > 0.0) THEN
            snow(index1(i)) = 0.0
          ELSE
            rain(index1(i)) = 0.0
          END IF
        END IF  
        ! Printing from a N48 UM8.1 run suggests the factorn is usually close
        ! to 1. so the increase is small.
        IF (qneg(i) < 0.0) THEN  
          ! negative increments present so can increase these to remove error
          ineg = ineg +1
          factorn(ineg) = 1.0 - qerror(i)/qneg(i) 
          indexn(ineg) = i            
        ELSE 
          ! Only positive increments so reduce these - odd must be removing 
          ! qcl or qcf? (This does appear to happen at present.)
          ip = ip +1
          factorp(ip) = 1.0 - qerror(i)/qpos(i) 
          indexp(ip) = i 
        END IF      
      END IF
    END IF
  END DO

  ! Correct negative dq 
  IF (ineg > 0) THEN
    DO k=1,nlev
      DO i = 1,ineg 
        ii=indexn(i)
        IF (dqbydt(index1(ii),k) < 0.0) THEN
          dqbydt(index1(ii),k) = dqbydt(index1(ii),k)*factorn(i)
        END IF
      END DO
    END DO
  END IF
  ! Correct positive dq 
  IF (ip > 0) THEN
    DO k=1,nlev
      DO i = 1,ip 
        ii=indexp(i)
        IF (dqbydt(index1(ii),k) > 0.0) THEN
          dqbydt(index1(ii),k) = dqbydt(index1(ii),k)*factorp(i)
        END IF
     END DO
    END DO
  END IF
 
  ! ----------------------------------------------------------------------
  ! Sum up mixing ratio and  +ve and -ve temperature increments
  ! ----------------------------------------------------------------------
  ! Moist energy E = cvT + gz + (Lc+Lf)q +Lfqcl (ignoring KE)
  !
  !  T = th * exner
  !
  ! For PC2 then 
  !
  ! dE = -lf*rain = sum{(cv*dth*exner+(Lc+lf)*dq+lf*dqcl)*rho*(r/a)*(r/a)*dr}
  !  Checking this makes the energy error more negative in most cases.
  ! ----------------------------------------------------------------------
  DO k=1,nlev
    DO i=1,npnts
      r2rhodr = r2rho_th(index1(i),k)*dr_across_th(index1(i),k)*ra2
      qsum(i)   = qsum(i)   + dqbydt(index1(i),k)  *r2rhodr  
      qclsum(i) = qclsum(i) + dqclbydt(index1(i),k)*r2rhodr  

      IF (dthbydt(index1(i),k)  >   0.0) THEN
        tspos(i) = tspos(i) + cv*dthbydt(index1(i),k)*                  &
                                exner_layer_centres(index1(i),k)*r2rhodr 
      ELSE
        tsneg(i) = tsneg(i) + cv*dthbydt(index1(i),k)*                  &
                                exner_layer_centres(index1(i),k)*r2rhodr 
      END IF
    END DO
  END DO

  ! ----------------------------------------------------------------------
  !  Calculate the error and apply the necessary energy correction
  !
  !   UM DOCUMENTATION PAPER 27
  !
  ! ----------------------------------------------------------------------

  DO i=1,npnts

    terr(i) = (lc+lf)*qsum(i) +lf*qclsum(i) + lf*rain(index1(i))             &
                          + tspos(i) + tsneg(i)

    ! scaling factor - safe value should get reset
    scalee(i) = 1.0    
  END DO

  nfalse=0
  DO i=1,npnts

    ! error is positive
    bposer(i) = terr(i)  >   0.0

    IF (bposer(i)) THEN
      ! Is it sensible to attempt to correct the positive increments by
      ! scaling them?
      IF ((tspos(i)  ==  0.0)) THEN
      ! No positive increments to correct
        bposer(i) = .FALSE.
      ELSE IF (terr(i) > tspos(i) .AND. tsneg(i) /= 0.0) THEN
      ! Correct negative increments instead. Not sure if this is a good idea.
      ! The alternative results in reversing the sign of the pos inc.
      ! Basicly the error in this case is too large to easily correct.
        bposer(i) = .FALSE.
      END IF
    ELSE    ! error is negative
      ! Is it sensible to attempt to correct the negative increments by
      ! scaling them?
      IF (tsneg(i)  ==  0.0) THEN
      ! No negative increments to correct
        bposer(i) = .TRUE.
      ELSE IF (terr(i) <  0.5*tsneg(i) .AND. tspos(i) /= 0.0) THEN
      ! Negative error too big to correct by reducing negative increments
      ! or correction of negative error will give a factor < 0.5
        bposer(i) = .TRUE.
      END IF
    END IF

    bcorr(i) = (tspos(i)  /=  0.0) .OR. (tsneg(i)  /=  0.0)

    IF (bposer(i) .AND. bcorr(i)) THEN
      scalee(i) = 1.0 - (terr(i)/tspos(i))
    ELSE IF (.NOT.bposer(i) .AND. bcorr(i)) THEN
      scalee(i) = 1.0 - (terr(i)/tsneg(i))
    END IF

    IF (.NOT.bcorr(i)) THEN
      ! Don't appear to be any of these false ascents with no theta inc
      nfalse = nfalse + 1
    ELSE  
      ! Stop correction if factors not reasonable 
      IF (scalee(i) > 1.6) THEN
        bcorr(i) = .FALSE.   ! more than 1.6 * orginal increments
        nfalse = nfalse + 1
      ELSE IF (scalee(i) <= 0.4 ) THEN
        bcorr(i) = .FALSE.   ! less than 0.4 * orginal increments
        nfalse = nfalse + 1
      END IF
    END IF
 
  END DO

  ! Scale positive or negative theta increments to correct energy
  DO k=1,nlev
    DO i=1,npnts
      IF (bcorr(i)) THEN   ! If to correct column

        IF (      ( bposer(i) .AND. (dthbydt(index1(i),k) > 0.0) )  &
        .OR. ( .NOT.bposer(i) .AND. (dthbydt(index1(i),k) < 0.0) ) )  THEN

          dthbydt(index1(i),k) = dthbydt(index1(i),k)*scalee(i)

        END IF

      END IF    ! test on correcting the column
    END DO  ! npnts
  END DO ! nlev


ELSE IF (ienergy_method == 2) THEN ! ienergy_method = 2

  ! ----------------------------------------------------------------------
  ! Method here for future development but unsafe at present
  ! ----------------------------------------------------------------------
  ! First interpolate dT and dq to rho levels after converting 
  ! dtheta to dT
  ! ----------------------------------------------------------------------

  DO k=1,nlev
    DO i=1,npnts
      dtbydt(i,k) = dthbydt(index1(i),k)*exner_layer_centres(index1(i),k)
    END DO  ! npnts
  END DO ! nlev
  ! level 1 on rho grid assumed to have same dt, dq as level 1 on theta grid
  DO i=1,npnts
    dtbydt_rho(i,1)   = dtbydt(i,1)
    dqbydt_rho(i,1)   = dqbydt(index1(i),1)
    dqclbydt_rho(i,1) = dqclbydt(index1(i),1)
    dqcfbydt_rho(i,1) = dqcfbydt(index1(i),1)
  END DO ! nlev
  DO k=2,nlev
    DO i=1,npnts
      j = index1(i)
      weight1 = r_theta(j,k) - r_rho(j,k)
      weight2 = r_rho(j,k) - r_theta(j,k-1)
      weight3 = r_theta(j,k) - r_theta(j,k-1)
      ww1 = weight1/weight3
      ww2 = weight2/weight3

      dtbydt_rho(i,k)   = ww2*dtbydt(i,k)   + ww1*dtbydt(i,k-1)
      dqbydt_rho(i,k)   = ww2*dqbydt(j,k)   + ww1*dqbydt(j,k-1)
      dqclbydt_rho(i,k) = ww2*dqclbydt(j,k) + ww1*dqclbydt(j,k-1)
      dqcfbydt_rho(i,k) = ww2*dqcfbydt(j,k) + ww1*dqcfbydt(j,k-1)
    END DO  ! npnts
  END DO ! nlev

  ! ----------------------------------------------------------------------
  !  Check moisture conservation and correct error.
  !  Assuming PC2 being used so increments to qcl & qcf.
  !
  !   total surface precip = sum { (dq+dqcl+qcf)*rho*(r/a)*(r/a)*dr}
  ! ----------------------------------------------------------------------
  ! integrate change in column moisture

  DO k= 1,nlev
    DO j=1,npnts
      i= index1(j) 
      r2rhodr = r2rho(i,k)*dr_across_rh(i,k)*ra2
      qloss(j) = qloss(j)+(dqbydt_rho(j,k)+dqclbydt_rho(j,k)+dqcfbydt_rho(j,k))&
                                 *r2rhodr
      IF (dqbydt_rho(j,k) > 0.0) THEN
        qpos(j) = qpos(j) + dqbydt_rho(j,k)*r2rhodr
      END IF
      IF (dqbydt_rho(j,k) < 0.0) THEN
        qneg(j) = qneg(j) + dqbydt_rho(j,k)*r2rhodr
      END IF
    END DO
  END DO

  DO i=1,npnts
    qerror(i) = qprecip(i) + qloss(i)
  END DO


  ! Attempt to correct moisture budget
  ! Is there precipitation - if so attempt to increase /reduce to match
  ! increments?
  ineg=0
  ip=0
  DO i=1,npnts
    IF (qerror(i) < 0.0) THEN
    ! should be more precip
      IF (qprecip(i) == 0.0) THEN
      ! Need to decide whether precip is rain or snow
        t1 = theta(index1(i),1)*exner_layer_centres(index1(i),1)  
        IF(t1 < tm) THEN
          snow(index1(i)) = snow(index1(i)) - qerror(i)
        ELSE
          rain(index1(i)) = rain(index1(i)) - qerror(i)
        END IF 
      ELSE
        IF (rain(index1(i)) > 0.0) THEN
          rain(index1(i)) = rain(index1(i)) - qerror(i)
        ELSE
          snow(index1(i)) = snow(index1(i)) - qerror(i)
        END IF
      END IF
    ELSE  ! Error is positive
    ! should be less precip
      IF (qprecip(i) >= qerror(i)) THEN
        IF (rain(index1(i)) > 0.0) THEN
          rain(index1(i)) = rain(index1(i)) - qerror(i)
        ELSE
          snow(index1(i)) = snow(index1(i)) - qerror(i)
        END IF
      ELSE
        ! Existing precip less than required reduction  
        ! Zero precip if any and then decide what else to do?
        IF (qprecip(i) > 0.0 )THEN
          qerror(i) = qerror(i) - qprecip(i)
          IF (snow(index1(i)) > 0.0) THEN
            snow(index1(i)) = 0.0
          ELSE
            rain(index1(i)) = 0.0
          END IF
        END IF  
        ! Printing from a N48 UM8.1 run suggests the factorn is usually close
        ! to 1. so the increase is small.
        IF (qneg(i) < 0.0) THEN  
          ! negative increments present so can increase these to remove error
          ineg = ineg +1
          factorn(ineg) = 1.0 - qerror(i)/qneg(i) 
          indexn(ineg) = i            
        ELSE 
          ! Only positive increments so reduce these - odd must be removing 
          ! qcl or qcf? (This does appear to happen at present.)
          ip = ip +1
          factorp(ip) = 1.0 - qerror(i)/qpos(i) 
          indexp(ip) = i 
        END IF      
      END IF
    END IF
  END DO

  ! Correct negative dq 
  IF (ineg > 0) THEN
    DO k=1,nlev
      DO i = 1,ineg 
        ii=indexn(i)
        IF (dqbydt_rho(ii,k) < 0.0) THEN
          dqbydt_rho(ii,k) = dqbydt_rho(ii,k)*factorn(i)
        END IF
      END DO
    END DO
    ! Need to interpolate dq correction back to theta grid
    DO i=1,ineg
      ii=indexn(i) 
      j = index1(ii)
      dqbydt(j,1)   = dqbydt_rho(ii,1) 
    END DO  ! npnts
    DO k=2,nlev
      DO i=1,ineg
        ii=indexn(i) 
        j = index1(ii)
        weight1 = r_rho(j,k) - r_theta(j,k-1)
        weight2 = r_theta(j,k-1) - r_rho(j,k-1)
        weight3 = r_rho(j,k) - r_rho(j,k-1)
        ww1 = weight1/weight3
        ww2 = weight2/weight3
        dqbydt(j,k)   = ww2*dqbydt_rho(ii,k)   + ww1*dqbydt_rho(ii,k-1)
      END DO  ! npnts
    END DO ! nlev

  END IF
  ! Correct positive dq 
  IF (ip > 0) THEN
    DO k=1,nlev
      DO i = 1,ip 
        ii=indexp(i)
        IF (dqbydt_rho(ii,k) > 0.0) THEN
          dqbydt_rho(ii,k) = dqbydt_rho(ii,k)*factorp(i)
        END IF
     END DO
    END DO
    ! Need to interpolate dq correction back to theta grid
    DO i=1,ip
      ii=indexp(i) 
      j = index1(ii)
      dqbydt(j,1)   = dqbydt_rho(ii,1) 
    END DO  ! npnts
    DO k=2,nlev
      DO i=1,ip
        ii=indexp(i) 
        j = index1(ii)
        weight1 = r_rho(j,k) - r_theta(j,k-1)
        weight2 = r_theta(j,k-1) - r_rho(j,k-1)
        weight3 = r_rho(j,k) - r_rho(j,k-1)
        ww1 = weight1/weight3
        ww2 = weight2/weight3
        dqbydt(j,k)   = ww2*dqbydt_rho(ii,k)   + ww1*dqbydt_rho(ii,k-1)
      END DO  ! npnts
    END DO ! nlev
  END IF

  ! ----------------------------------------------------------------------
  ! Sum up mixing ratio and  +ve and -ve temperature increments
  ! ----------------------------------------------------------------------
  ! Moist energy E = cvT + gz + (Lc+Lf)q +Lfqcl (ignoring KE)
  !
  !  T = th * exner
  !
  ! For PC2 then 
  !
  ! dE = -lf*rain = sum{(cv*dth*exner+(Lc+lf)*dq+lf*dqcl)*rho*(r/a)*(r/a)*dr}
  !  Checking this makes the energy error more negative in most cases.
  ! ----------------------------------------------------------------------
  DO k=1,nlev
    DO i=1,npnts
      j= index1(i) 
      r2rhodr = r2rho(j,k)*dr_across_rh(j,k)*ra2
      qsum(i) = qsum(i) + dqbydt_rho(i,k)*r2rhodr  
      qclsum(i) = qclsum(i) + dqclbydt_rho(i,k)*r2rhodr  

      IF (dtbydt_rho(i,k)  >   0.0) THEN
        tspos(i) = tspos(i) + cv*dtbydt_rho(i,k)*r2rhodr 
      ELSE
        tsneg(i) = tsneg(i) + cv*dtbydt_rho(i,k)*r2rhodr 
      END IF
    END DO
  END DO

  ! ----------------------------------------------------------------------
  !  Calculate the error and apply the necessary energy correction
  !
  !   UM DOCUMENTATION PAPER 27
  !
  ! ----------------------------------------------------------------------

  DO i=1,npnts
    terr(i) = (lc+lf)*qsum(i) +lf*qclsum(i) + lf*rain(index1(i))             &
                          + tspos(i) + tsneg(i)
    scalee(i) = 1.0
  END DO

  nfalse=0
  DO i=1,npnts

    bposer(i) = terr(i)  >   0.0

    IF (bposer(i)) THEN
      IF ((tspos(i)  ==  0.0)) THEN
      ! No positive increments to correct
        bposer(i) = .FALSE.
      ELSE IF (terr(i) > tspos(i) .AND. tsneg(i) /= 0.0) THEN
      ! Correct negative increments instead. Not sure if this is a good idea.
      ! The alternative results in reversing the sign of the pos inc.
      ! Basicly the error in this case is too large to easily correct.
        bposer(i) = .FALSE.
      END IF
    ELSE    ! error is negative
      IF (tsneg(i)  ==  0.0) THEN
      ! No negative increments to correct
        bposer(i) = .TRUE.
      ELSE IF (terr(i) <  0.5*tsneg(i) .AND. tspos(i) /= 0.0) THEN
      ! Negative error too big to correct by reducing negative increments
      ! or correction of negative error will give a factor < 0.5
        bposer(i) = .TRUE.
      END IF
    END IF

    bcorr(i) = (tspos(i)  /=  0.0) .OR. (tsneg(i)  /=  0.0)

    IF (bposer(i) .AND. bcorr(i)) THEN
      scalee(i) = 1.0 - (terr(i)/tspos(i))
    ELSE IF (.NOT.bposer(i) .AND. bcorr(i)) THEN
      scalee(i) = 1.0 - (terr(i)/tsneg(i))
    END IF

    IF (.NOT.bcorr(i)) THEN
      ! Don't appear to be any of these false ascents with no theta inc
      nfalse = nfalse + 1
    ELSE  
      ! Stop correction if factors not reasonable 
      IF (scalee(i) > 1.6) THEN
        bcorr(i) = .FALSE.   ! more than 1.6 * orginal increments
        nfalse = nfalse + 1
      ELSE IF (scalee(i) <= 0.4 ) THEN
        bcorr(i) = .FALSE.   ! less than 0.4 * orginal increments
        nfalse = nfalse + 1
      END IF
    END IF
 
  END DO


  ! Scale positive or negative T increments to correct energy
  DO k=1,nlev
    DO i=1,npnts
      IF (bcorr(i)) THEN   ! If to correct column

        IF  (       ( bposer(i) .AND. (dtbydt_rho(i,k)  >   0.0) )        &
          .OR. ( .NOT.bposer(i) .AND. (dtbydt_rho(i,k)  <   0.0) ) )  THEN

          dtbydt_rho(i,k) = dtbydt_rho(i,k)*scalee(i)

        END IF

      END IF    ! test on correcting the column
    END DO  ! npnts
  END DO ! nlev


  !----------------------------------------------------------------------- 
  ! Intepolate dT back to theta grid and convert back to dtheta 
  !----------------------------------------------------------------------- 
  ! level 1 on rho grid assumed to have same dt, dq as level 1 on theta grid
  DO i=1,npnts
    dtbydt(i,1)   = dtbydt_rho(i,1)
  END DO ! nlev
  DO k=2,nlev
    DO i=1,npnts
      j = index1(i)
      weight1 = r_rho(j,k) - r_theta(j,k-1)
      weight2 = r_theta(j,k-1) - r_rho(j,k-1)
      weight3 = r_rho(j,k) - r_rho(j,k-1)
      ww1 = weight1/weight3
      ww2 = weight2/weight3

      dtbydt(i,k)   = ww2*dtbydt_rho(i,k)   + ww1*dtbydt_rho(i,k-1)
    END DO  ! npnts
  END DO ! nlev

  DO k=1,nlev
    DO i=1,npnts
      dthbydt(index1(i),k) =dtbydt(i,k)/exner_layer_centres(index1(i),k)
    END DO  ! npnts
  END DO ! nlev


END IF  ! test on ienergy_method

!------------------------------------------------------------------------
! problem ascents - only prints if something not safe
!------------------------------------------------------------------------
IF (nfalse > 0) THEN
  DO i=1,npnts
    IF (.NOT.bcorr(i) ) THEN
      WRITE(6,'(a21,i6,a8,g26.18,a5,6g26.18)') 'Conv Energy problem: ',i,    &
              ' factor ',scalee(i),' err ',terr(i),tspos(i),tsneg(i),         &
                (lc+lf)*qsum(i),(lf)*qclsum(i),qerror(i)
      ! No appropriate format as if there is a problem want to know 
      ! variables to machine accuracy & nlev varies
      WRITE(6,*) ' dtheta ',(dthbydt(index1(i),k),k=1,nlev)
      WRITE(6,*) ' dq   ',(dqbydt(index1(i),k),k=1,nlev)
      WRITE(6,*) ' dqcl ',(dqclbydt(index1(i),k),k=1,nlev)
      WRITE(6,*) ' dqcf ',(dqcfbydt(index1(i),k),k=1,nlev)
    END IF
  END DO   
END IF

! Expand back to full array
DO i=1,npnts
  j=index1(i)
  energy_err(j) = terr(i)
END DO

!------------------------------------------------------------------------
! Calculate change in KE due to convection - values on rho grid
!------------------------------------------------------------------------
ike_cal = 1
IF (ike_cal == 1 .AND. l_mom) THEN
  DO i=1,npnts
    KE_loss(i) = 0.0
  END DO
  DO k=1,nlev
    DO i=1,npnts
      r2rhodr = r2rho(index1(i),k)*dr_across_rh(index1(i),k)*ra2
      KE_loss(i) = KE_loss(i) + 0.5*(2.0*dubydt(index1(i),k)*u(index1(i),k) &
                                   + 2.0*dvbydt(index1(i),k)*v(index1(i),k) &
                         + dubydt(index1(i),k)*dubydt(index1(i),k)          &
                         + dvbydt(index1(i),k)*dvbydt(index1(i),k) )*r2rhodr  
    END DO
  END DO


  ! Expand back to full array
  DO i=1,npnts
    j=index1(i)
    ke_cmt_loss(j) = KE_loss(i)
  END DO
 
END IF
!------------------------------------------------------------------------
IF (lhook) CALL dr_hook('COR_ENGY_5A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE cor_engy_5a
