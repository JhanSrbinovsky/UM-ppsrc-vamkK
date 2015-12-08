! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the mass flux profile for deep CMT.
!
SUBROUTINE cmt_mass(np_field, nconv, nlevs, nterm, cu_term,                 &
                    kterm, cu_tend, n_0degc, nlcl, ntop,                    &
                    mb, p_0degc, plcl, ptop, phalf, p,                      &
                  ! Output arguments
                    mass_up,mass_dwn,visc)

USE cv_run_mod, ONLY:                                                       &
    deep_cmt_opt

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description : 
!  Calculates the mass flux profile for deep convection to be used in CMT
!  calculations. Uses the cloud-base mass flux from the plume scheme, but 
!  profile is not the same as used for the thermodynamic part of the 
!  convection scheme.
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
!------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  np_field             & ! Full field length
 ,nconv                & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nterm                & ! Number of points terminating
 ,cu_term(nterm)       & ! Indices for terminating points
 ,kterm(np_field)      & ! Terminating levels
 ,cu_tend(nterm)       & !
 ,n_0degc(nconv)       & ! Level corresponding to melting level
 ,nlcl(nconv)          & ! Lifting condensation level
 ,ntop(nconv)            ! Top level of convection

REAL, INTENT(IN)    :: & 
  mb(nconv)            & ! Cloud base mass flux
 ,p_0degc(nconv)       & ! Pressure of melting level (hPa)
 ,plcl(nconv)          & ! Pressure of LCL (hPa)
 ,ptop(nconv)          & ! Pressure at top of convection (hPa)
 ,phalf(nlevs,nconv)   & ! Pressure on model half levels (hPa)
 ,p(nlevs,nconv)         ! Pressure on model levels (hPa)

REAL, INTENT(OUT) ::    &
  mass_up(nlevs,nconv)  & ! Updraught mass flux profile (Pa/s)
 ,mass_dwn(nlevs,nconv) & ! Downdarught mass flux 
 ,visc(nlevs,nconv)       ! Viscosity


! Local variables

INTEGER ::        & 
  i,j,k,n           ! loop counters
INTEGER ::        & 
  error             ! error code

REAL ::           &
  beta_l          &
 ,beta_u          &
 ,a_m(nterm)      &
 ,alpha_up        &
 ,alpha_dwn       &
 ,delta_p         &
 ,shape_factor

CHARACTER (len=8), PARAMETER ::  routinename = 'cmt_mass'
CHARACTER (len=80) :: cmessage        ! error message

! Parameters (in future consider putting in a module).

REAL, PARAMETER ::  &
  a_max=1.5         & ! Ratio mass flux at melting level to cloud base flux
 ,a_top=0.2         & ! Ratio mass flux at top to cloud base flux
 ,a_dwn=0.0         & ! 
 ,p_ref=60000.0     & ! Used to scale ratios as p_0degC increase
 ,p_min=10000.0     & !
 ,top_press=15000.0 & !
 ,alpha_visc=0.30     !

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!------------------------------------------------------------------------
IF (lhook) CALL dr_hook('CMT_MASS',zhook_in,zhook_handle)
!------------------------------------------------------------------------
! Determine shape of mass flux profile. If the depth of convection below
! the zero degree level is large enough, and the depth above is as well,
! the mass flux profile is assumed to have a maximum at the zero degree
! level.

DO i=1,nterm
  j=cu_term(i)
  IF((p_0degc(j)-ptop(j)) >= p_min .AND. (plcl(j)-p_0degc(j)) >= p_min) THEN

! Mass flux profile has an elevated maximum, the value of the maximum
! mass flux increase with decreasing pressure of the zero degree level
! while it is below p_ref. Above p_ref the maximum value is fixed.

    IF(p_0degc(j) >= p_ref) THEN
      a_m(i)=1.0+(a_max-1.0)*(p_0degc(j)-(plcl(j)-p_min))/           &
                             (p_ref     -(plcl(j)-p_min))
    ELSE
      a_m(i)=a_max
    END IF
  ELSE
    a_m(i)=1.0
  END IF
END DO

DO i=1,nterm
  j=cu_term(i)
  n=cu_tend(i)
  mass_up(nlcl(j),j)=mb(j)
  mass_dwn(nlcl(j),j)=a_dwn*mb(j)
  beta_l=LOG(a_m(i))
  beta_u=LOG(a_m(i)/a_top)
  IF((p_0degc(j)-ptop(j)) >= p_min .AND. (plcl(j)-p_0degc(j)) >= p_min) THEN

    DO k=nlcl(j)+1,ntop(j)+1
      IF(phalf(k,j) >= p_0degc(j)) THEN
        mass_up(k,j)=a_m(i)*mb(j)*                                     &
                      EXP(-beta_l*((phalf(k,j)-p_0degc(j))/            &
                                   (plcl(j)-p_0degc(j)))**2)
        mass_dwn(k,j)=a_dwn*mb(j)

      ELSE
        mass_up(k,j)=a_m(i)*mb(j)*                                     &
                      EXP(-beta_u*((phalf(k,j)-p_0degc(j))/            &
                                   (ptop(j)-p_0degc(j)))**2)
        mass_dwn(k,j)=0.0
      END IF
    END DO
  ELSE
    DO k=nlcl(j)+1,ntop(j)+1

! Choose either the diagnosed top of convection or the level at which
! scheme detrains, whichever is higher (this prevents odd failures
! in the CMT scheme at mid-latitudes)

      IF(ptop(j) >= phalf(kterm(n)+1,j)) THEN
        mass_up(k,j)=a_m(i)*mb(j)*                                    &
                      EXP(-beta_u*((phalf(k,j)-plcl(j))/              &
                                   (phalf(kterm(n)+1,j)-plcl(j)))**2)
      ELSE
        mass_up(k,j)=a_m(i)*mb(j)*                                    &
                     EXP(-beta_u*((phalf(k,j)-plcl(j))/               &
                                  (ptop(j)-plcl(j)))**2)
      END IF
      mass_dwn(k,j)=0.0
    END DO
  END IF
END DO

! Calculate eddy viscosity (Note viscosity=0 at NLCL and  KTERM+2)

DO i=1,nterm
  j=cu_term(i)
  visc(nlcl(j),j)  =0.0
  visc(ntop(j)+2,j)=0.0
  delta_p=-(ptop(j)-plcl(j))

! alpha_up is a very crude representation of the non-dimensional profile
! of updraught vertical velocity. Tuning of the viscosity was done in the 
! SCM, comparing momentum fluxes with those derived from CRM simulation
! of periods during TOGA-COARE.
! Allows for possibility of a downdraught component, but set to zero at 
! present.

  IF (deep_cmt_opt == 0) THEN

    DO k=nlcl(j)+1,ntop(j)+1
      alpha_up=2.0*(phalf(k,j)-plcl(j))/(top_press-plcl(j))
      IF(p_0degc(j) < plcl(j) .AND. phalf(k,j) > p_0degc(j)) THEN
        alpha_dwn =2.0*SQRT((phalf(k,j)-p_0degc(j))/(plcl(j)-p_0degc(j)))
      ELSE
        alpha_dwn =0.0
      END IF
      visc(k,j)=(alpha_up*mass_up(k,j)+alpha_dwn*mass_dwn(k,j))*           &
                  alpha_visc*delta_p
    END DO

  ELSE IF (deep_cmt_opt == 1) THEN

! Added a shape factor to control gradient term
! Note as a_dwn=0.0 removed downward terms from calculation as waste of CPU

    DO k=nlcl(j)+1,ntop(j)+1
      alpha_up=2.0*(phalf(k,j)-plcl(j))/(top_press-plcl(j))
      shape_factor = 1.0-0.75*(phalf(k,j)-plcl(j))/(ptop(j)-plcl(j))
      visc(k,j)=alpha_up*mass_up(k,j)*alpha_visc*delta_p*shape_factor
    END DO

  ELSE IF (deep_cmt_opt == 5) THEN

! Added a quadratic shape factor

    DO k=nlcl(j)+1,ntop(j)+1
      alpha_up=2.0*(phalf(k,j)-plcl(j))/(top_press-plcl(j))
      shape_factor = (1.0-0.75*(phalf(k,j)-plcl(j))/(ptop(j)-plcl(j)))**2
      visc(k,j)=alpha_up*mass_up(k,j)*alpha_visc*delta_p*shape_factor
    END DO

  ELSE       ! unallowed option

    error = 1   
    WRITE(cmessage,'(a35,i3)') 'Unacceptable value of deep_cmt_opt ',      &
                deep_cmt_opt
    CALL ereport(routinename,error,cmessage)

  END IF

END DO

IF (lhook) CALL dr_hook('CMT_MASS',zhook_out,zhook_handle)

RETURN
END SUBROUTINE cmt_mass
