! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Convection cloud Microphysics Scheme.
! Subroutine Interface:
SUBROUTINE cloud_w_4a5a (k, npnts, xpkp1, qclpkp1, qcfpkp1, prekp1, xsqkp1    &
               , blowst, flxkp1, xpk, qclpk, qcfpk, thekp1,qekp1, bwkp1       &
               , bland, qsekp1, bgmkp1, bterm, cca, iccb, icct, tcw, depth    &
               , ekp14, ekp34, delexkp1, cclwp, delpkp1, ccw, lcca, lcbase    &
               , lctop, l_shallow, l_q_interact, start_lev)

USE atmos_constants_mod, ONLY: cp, c_virtual
 
Use cv_run_mod, Only:                                                   &
    l_fix_udfactor, l_ccrad,                                            &
    ud_factor, mparwtr, qlmin, fac_qsat, ccw_for_precip_opt

USE cv_param_mod, ONLY:                                                 &
    critdsea, critdlnd, critdice

USE earth_constants_mod, ONLY: g

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Description : cloud microphysics routine
!
! Calculates preciptation produced in lifting parcel from layer k to k+1
!
! Calls CON_RAD to calculate parameters for radiation calculation
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
! Subroutine arguments (Order does not conform to coding standards as very old
!                       routine.)

INTEGER, INTENT(IN) :: &
  k                    & ! Present model layer number
 ,npnts                  ! Vector length

INTEGER, INTENT(IN) :: & 
  start_lev(npnts)       ! Level at which convection initiated 

LOGICAL, INTENT(IN) :: &
  blowst(npnts)        & ! Mask for those points at which stability is low 
                         ! enough for convection to occur
 ,bwkp1(npnts)         & ! Mask for whether condensate is liquid in layer k+1

 ,bland(npnts)         & ! Land/sea mask

 ,bgmkp1(npnts)        & ! Mask for parcels which are saturated in layer k+1

 ,bterm(npnts)         & ! Mask for parcels which terminate in layer k+1

 ,l_shallow(npnts)     & ! Mask for points where shallow convection is likely

 ,l_q_interact           ! Switch allows overwriting of parcel variables 
                         ! (will alter results).


REAL, INTENT(IN) ::     &
  xsqkp1(npnts)         & ! Excess parcel mixing ratio in layer k+1 (kg/kg)

 ,flxkp1(npnts)         & ! Parcel mass flux in layer k+1 (Pa/s)

 ,qclpk(npnts)          & ! Parcel liquid condensate mixing ratio in layer k
                          ! (kg/kg)
 ,qcfpk(npnts)          & ! Parcel frozen condensate mixing ratio in layer k
                          ! (kg/kg)
 ,thekp1(npnts)         & ! Potential temperature of cloud environment 
                          !    in layer k+1 (K)
 ,qekp1(npnts)          & ! mixing ratio of cloud environment in layer k+1
                          ! (kg/kg)
 ,qsekp1(npnts)         & ! Saturation mixing ratio of cloud environment
                          !    in layer k+1 (kg/kg)
 ,ekp14(npnts)          & ! entrainment rate at level k+1/4
                          ! multiplied by appropriate layer thickness
 ,ekp34(npnts)          & ! entrainment rate at level k+3/4
                          ! multiplied by appropriate layer thickness
 ,delexkp1(npnts)       & ! difference in exner ratio across layer k+1 (Pa)

 ,delpkp1(npnts)          ! Pressure difference across layer k+1


REAL, INTENT(INOUT) ::  &
  qclpkp1(npnts)        & ! IN Parcel liquid condensate mixing ratio
                          !    in layer k+1 before excess mixing
                          !    ratio water is added(kg/kg)
                          ! OUT Parcel liquid condensate mixing ratio
                          !     in layer k+1 (kg/kg)
 ,qcfpkp1(npnts)        & ! IN Parcel frozen condensate mixing ratio
                          !    in layer k+1 before excess mixing
                          !     ratio water is added(kg/kg)
                          ! OUT Parcel frozen condensate mixing ratio
                          !      in layer k+1 (kg/kg)
 ,xpk(npnts)            & ! IN Parcel cloud water in layer k (kg/kg)
                          ! OUT overwritten with qcl+qcf for layer k
                          ! if l_q_interact set to .true.
 ,tcw(npnts)            & ! IN  Total condensed water summed up to layer k
                          !     (kg/m**2/s)
                          ! OUT Total condensed water summed up to layer k
                          !     (kg/m**2/s)
 ,depth(npnts)          & ! IN  Depth of convective cloud to layer k (m)
                          ! OUT Depth of convective cloud to layer k+1 (m)
 ,cclwp(npnts)            ! IN  Condensed water path summed up to layer k
                          !     (kg/m**2)
                          ! OUT Condensed water path summed up to layer k+1
                          !     (kg/m**2)

INTEGER, INTENT(OUT) :: &
  iccb(npnts)           & ! Convective cloud base level

 ,icct(npnts)           & ! Convective cloud top level

 ,lcbase(npnts)         & ! Lowest convective cloud base level

 ,lctop(npnts)            ! Lowest convective cloud top level


REAL, INTENT(OUT) ::    &
  xpkp1(npnts)          & ! Parcel cloud water in layer k+1 (kg/kg)

 ,prekp1(npnts)         & ! Precipitation from parcel as it rises from layer
                          !  k to k+1 (kg/m**2/s)
 ,cca(npnts)            & ! Convective cloud amount (%)

 ,ccw(npnts)            & ! Convective cloud liquid water (g/kg) on model levels

 ,lcca(npnts)             ! Lowest convective cloud amount (%)



!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------

INTEGER :: I              ! Loop counter

LOGICAL, PARAMETER :: l_fallout =.TRUE.
!                              .true. is default. User DO NOT TOUCH!
!                              .false. kills precipitation
!                                      DIAGNOSTIC TEST CASE ONLY:
!                                     (prevent convection updating qT).

REAL ::          &
  dcrit          &  ! Critical depth at which Precipitation may form (m)

 ,mparmult       &  ! Factor used to multiply mparwtr, value between
                    ! 1. and ~3. being 1.0 for deeper clouds.   

 ,xmin           &  ! amount of cloud water retained by the parcel on 
                    ! Precipitation (kg/kg)

 ,epss           &  ! (1.0+ekp14)*(1.0+ekp34)

 ,ccw_ud(npnts)     ! cloud water for radiation


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CLOUD_W_4A5A',zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! Calculate parcel cloud condensate before Precipitation
!
! UM Documentation Paper 27
! Section (2B), Equation (13A or 13C)
!-------------------------------------------------------------------------------


DO i=1, npnts
  epss = (1.0+ekp14(i)) * (1.0+ekp34(i))
  xpkp1(i) = ( xpk(i)/epss ) + xsqkp1(i)

END DO

IF (l_q_interact) THEN

!  xpk(p1) are used for inputs only: need to work out how to use
!  qclpk and qcfpk instead. Meanwhile, overwrite xpk(p1).
 
 DO i=1, npnts
    IF (bwkp1(i)) THEN
      qclpkp1(i) = qclpkp1(i) + xsqkp1(i)
    ELSE
      qcfpkp1(i) = qcfpkp1(i) + xsqkp1(i)
    END IF
    xpk  (i) = qclpk  (i) + qcfpk  (i)
    xpkp1(i) = qclpkp1(i) + qcfpkp1(i)
  END DO 

END IF

!-------------------------------------------------------------------------------
! Store convective cloud liquid water before Precipitation
!-------------------------------------------------------------------------------
!
DO i=1,npnts
  ccw(i)     = xpkp1(i)
  prekp1(i)  = 0.0       ! initialise precipitation to zero
END DO

!-------------------------------------------------------------------------------
! Calculate cloud depth and assign crtical cloud depths
!
! UM Documentation Paper 27
! Section (8), Equation (34), (35)
!-------------------------------------------------------------------------------

DO I=1,npnts
  IF ( blowst(i) ) depth(i) = 0.0

! This could be improved and simplified by using actual model heights

  IF ( bgmkp1(i) )                                                         &
    depth(i) = depth(i) + ( cp * thekp1(i) * (1.0+c_virtual*qekp1(i)) *    &
                                                        delexkp1(i)/G )
END DO  

!-------------------------------------------------------------------------------
! Calculate Precipitation from layer k+1 and adjust cloud water
!
! UM Documentation Paper 27
! Section (8), Equation (36)  (ccw_for_precip_opt=0)
!-------------------------------------------------------------------------------

SELECT CASE (ccw_for_precip_opt)

  CASE (4)            ! xmin profile based on qsat with user defined 
                      ! minimum, maximum and qsat scaling.
    DO i=1,npnts
      
      xmin = MAX(MIN(mparwtr, fac_qsat*qsekp1(i)), qlmin)

      ! Precipitate if cloud water in the layer > xmin

      IF ( xpkp1(i) > xmin ) THEN

        prekp1(i) = (xpkp1(i) -xmin) * flxkp1(i) /g 

        IF (l_q_interact) THEN        ! PC2
        ! Update the parcel's liquid and frozen cloud condensate
          qclpkp1(i) = qclpkp1(i) * xmin / xpkp1(i)
          qcfpkp1(i) = qcfpkp1(i) * xmin / xpkp1(i)
        END IF 

        xpkp1(i) = xmin      ! cloud water in layer after precip

      END IF     ! test on whether to precipitate

      ! Note no application of ud_factor

      ccw_ud(i) = xpkp1(i)  

    END DO

  CASE (3)            ! Manoj's function for congestus
                                          
    Do i=1,npnts

      mparmult = 1.5 + 0.5*tanh((3000.0 -depth(i))/800.0)
      xmin = MIN(mparwtr*mparmult, fac_qsat*qsekp1(i)) 

      ! If a land point and liquid water in the layer 
      ! increase the minimum cloud water for Precipitation
      ! The reasons for this are an attempt to take some account of 
      ! more aerosols over land leading to more small cloud drops and
      ! therefore more cloud water before Precipitation.  

      IF (bwkp1(i) .and. bland(i)) xmin =xmin*2.0

      ! limit max value to 0.003 kg/kg
      xmin = MIN(xmin,0.003)          

      IF (l_q_interact) THEN   ! PC2
      ! Limit xmin
        xmin = max(xmin, qlmin)
      END IF

      ! Precipitate if cloud water in the layer > xmin

      IF ( xpkp1(i) > xmin ) THEN

        prekp1(i) = (xpkp1(i) -xmin) * flxkp1(i) /g 

        IF (l_q_interact) THEN        ! PC2
        ! Update the parcel's liquid and frozen cloud condensate
          qclpkp1(i) = qclpkp1(i) * xmin / xpkp1(i)
          qcfpkp1(i) = qcfpkp1(i) * xmin / xpkp1(i)
        END IF 

        xpkp1(i) = xmin      ! cloud water in layer after precip

      END IF     ! test on whether to precipitate

      ! Note no application of ud_factor

      ccw_ud(i) = xpkp1(i)  

    END DO

  CASE (2)            ! Manoj's function dependent on depth of cloud
                      ! Also removed extra test on l_fallout

    DO i=1,npnts

      mparmult = 2.0 + 1.0*tanh((1500.0 -depth(i))/1000.0)
      xmin = MIN(mparwtr*mparmult, fac_qsat*qsekp1(i)) 

      ! If a land point and liquid water in the layer 
      ! increase the minimum cloud water for Precipitation
      ! The reasons for this are an attempt to take some account of 
      ! more aerosols over land leading to more small cloud drops and
      ! therefore more cloud water before Precipitation.  

      IF (bwkp1(i) .and. bland(i)) xmin =xmin*2.0

      IF (l_q_interact) THEN   ! PC2
      ! Limit xmin
        xmin = max(xmin, qlmin)
      END IF

      ! Precipitate if cloud water in the layer > xmin

      IF ( xpkp1(i) > xmin ) THEN

        prekp1(i) = (xpkp1(i) -xmin) * flxkp1(i) /g 

        IF (l_q_interact) THEN        ! PC2
        ! Update the parcel's liquid and frozen cloud condensate
          qclpkp1(i) = qclpkp1(i) * xmin / xpkp1(i)
          qcfpkp1(i) = qcfpkp1(i) * xmin / xpkp1(i)
        END IF 

        xpkp1(i) = xmin      ! cloud water in layer after precip

      END IF     ! test on whether to precipitate

      ! Note no application of ud_factor

      ccw_ud(i) = xpkp1(i)  

    END DO

  CASE (1)            ! As old L_no_dcrit = .true. option
                      ! No test on a critical depth for Precipitation
                      ! Option in use for HadGEM2

    Do i=1,npnts
      xmin = MIN(mparwtr, fac_qsat*qsekp1(i)) 

      ! If a land point and liquid water in the layer 
      IF (bwkp1(i) .and. bland(i)) xmin =xmin*2.0

      IF (l_q_interact) THEN
      ! Limit xmin
        xmin = max(xmin, qlmin)
      END IF

      ! Precipitate if cloud water in the layer > xmin
! Original test with l_fallout left in but commented out so that could be
! reused for investigation
!      IF ( (xpkp1(i) > xmin) .and. l_fallout) THEN 
      IF ( (xpkp1(i) > xmin)) THEN

        prekp1(i) = (xpkp1(i) -xmin) * flxkp1(i) /g 

        IF (l_q_interact) THEN        ! PC2
        ! Update the parcel's liquid and frozen cloud condensate
          qclpkp1(i) = qclpkp1(i) * xmin / xpkp1(i)
          qcfpkp1(i) = qcfpkp1(i) * xmin / xpkp1(i)
        END IF 

        xpkp1(i) = xmin
        IF (l_ccrad .OR. l_fix_udfactor) THEN
          ccw_ud(i)= xpkp1(i)
        ELSE
          ccw_ud(i)= xpkp1(i)*ud_factor 
        END IF 
      ELSE     ! no Precipitation

        ! UD_FACTOR ought to be applied here too but, in order to 
        ! simplify the code, this is corrected by moving the point
        ! at which UD_FACTOR is applied to after convection 
        ! in CONV_CTL and passing UD_FACTOR of 1 to here
        ccw_ud(i) = xpkp1(i)  

      END IF     ! test on whether to precipitate

    END DO

  CASE DEFAULT       ! 0 - original convection code using a critical depth

    DO i=1,npnts

      ! Assign critical cloud depths
  
      IF (.not.bwkp1(i)) THEN    ! Ice water cloud in this layer 
        dcrit = critdice
      ELSE IF (bland(i)) THEN    ! liquid water over land
        dcrit = critdlnd
      ELSE                       ! liquid water over sea 
        dcrit = critdsea
      END IF
          
      xmin = MIN(mparwtr, fac_qsat*qsekp1(i)) 

      IF (l_q_interact) THEN
      ! Limit xmin
        xmin = max(xmin, qlmin)
      END IF

      ! precipitate if 
      ! Either  depth of cloud > dcrit OR convection is not shallow
      ! and    cloud water in the layer > xmin 
! orignal code with l_fallout left commented out so can be reused for tests
!      IF ( ( (depth(i) > dcrit).or.(.not.l_shallow(i)) )  &
!               .and. (xpkp1(i) > xmin) .and. l_fallout) THEN
      IF ( ( (depth(i) > dcrit).OR.(.NOT.l_shallow(i)) )  &
           .AND. (xpkp1(i) > xmin) ) THEN

        prekp1(i) = (xpkp1(i) -xmin) * flxkp1(i) /g 

        IF (l_q_interact) THEN        ! PC2
        ! Update the parcel's liquid and frozen cloud condensate
          qclpkp1(i) = qclpkp1(i) * xmin / xpkp1(i)
          qcfpkp1(i) = qcfpkp1(i) * xmin / xpkp1(i)
        END IF 
        xpkp1(i) = xmin

        IF (l_ccrad .OR. l_fix_udfactor) THEN
          ccw_ud(i)= xpkp1(i)
        ELSE
          ccw_ud(i)= xpkp1(i)*ud_factor 
        END IF 

      ELSE     ! no Precipitation

        ! UD_FACTOR ought to be applied here too but, in order to 
        ! simplify the code, this is corrected by moving the point
        ! at which UD_FACTOR is applied to after convection 
        ! in CONV_CTL and passing UD_FACTOR of 1 to here
        ccw_ud(i) = xpkp1(i)  

      END IF     ! test on whether to precipitate

    END DO

END SELECT      ! test on value of ccw_for_precip_opt 

!-------------------------------------------------------------------------------
! Calculate convective cloud base, convective cloud top, total condensed
! water/ice and convective cloud amount
!
! SUBROUTINE CON_RAD - now called after rainout has occurred
! UM Documentation paper 27
! Section (9)
!-------------------------------------------------------------------------------

! DEPENDS ON: con_rad_4a5a
CALL con_rad_4a5a(k,xpk,ccw_ud,flxkp1,bterm,cca,iccb,icct,start_lev,         &
                  tcw,ccw,cclwp,delpkp1,lcca,lcbase,lctop,npnts,l_q_interact)

!-------------------------------------------------------------------------------
! Store convective cloud liquid water after Precipitation
!-------------------------------------------------------------------------------

Do i=1,npnts
  ccw(i) = ccw_ud(i)
END DO

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('CLOUD_W_4A5A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE cloud_w_4a5a
