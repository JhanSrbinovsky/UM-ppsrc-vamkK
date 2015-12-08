! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE H_INT_CO----------------------------------------------
!
!    Purpose:  Calculates bi-linear horizontal interpolation
!              coefficients and gather indices for interpolating
!              between generalised latitude-longitude grids (eg
!              global, regional or rotated lat-lon grid) in which the
!              gridlength may vary with latitude and/or longitude. The
!              interpolation is carried out by subroutine
!              H_INT. Gather indices point to bottom left hand
!              corner and bottom right hand corner of each grid box on
!              source grid enclosing a target point. Two indices are
!              needed to cater for east-west (lambda direction) cyclic
!              boundaries when the source data is global. If a target po
!              falls outside the domain of the source data, one sided
!              differencing is used. The source latitude coordinates
!              must be supplied in decreasing order. The source long-
!              itude coordinates must be supplied in increasing order,
!              starting at any value, but not wrapping round. The
!              target points may be specified in any order.
!
!
!    Programming standard:
!             UMDP3 v8.2
!
!    System component: S121,S122
!
!    System task: S1
!
!    Documentation:
!              The interpolation formulae are described in
!              unified model on-line documentation paper S1.
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Grids
MODULE h_int_co_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE h_int_co                                               &
(index_b_l,index_b_r,weight_t_r,weight_b_r,weight_t_l,weight_b_l  &
,lambda_srce,phi_srce,lambda_targ,phi_targ                        &
,points_lambda_srce,points_phi_srce,points,cyclic)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN)  :: points_lambda_srce 
                          ! Number of lambda points on source grid
INTEGER, INTENT(IN)  :: points_phi_srce
                          ! Number of phi points on source grid
INTEGER, INTENT(IN)  :: points
                          ! Total number of points on target grid
                     
INTEGER, INTENT(OUT) :: index_b_l(points) 
                          ! Index of bottom lefthand corner of source gridbox
INTEGER, INTENT(OUT) :: index_b_r(points)  
                          ! Index of bottom righthand corner of source gridbox

REAL, INTENT(IN)     :: lambda_targ(points) 
                          !   Lambda coords of target grid in degrees
                          !   using same rotation as source grid
REAL, INTENT(IN)     :: phi_targ(points)    
                          !   Phi coords of target grid in degrees
                          !   using same rotation as source grid
REAL, INTENT(IN)     :: lambda_srce(points_lambda_srce)
                          !   Lambda coords of source grid in degrees                      
REAL, INTENT(IN)     :: phi_srce(points_phi_srce)
                          !   Phi coords of target grid in degree

REAL, INTENT(OUT)    :: weight_t_r(points)
                          !   Weight applied to value at top right
                          !   hand corner of source gridbox
REAL, INTENT(OUT)    :: weight_b_l(points)
                          !   Weight applied to value at bottom left
                          !   hand corner of source gridbox
REAL, INTENT(OUT)    :: weight_b_r(points)
                          !   Weight applied to value at bottom right
                          !   hand corner of source gridbox
REAL, INTENT(OUT)    :: weight_t_l(points)
                          !   Weight applied to value at top left
                          !   hand corner of source gridbox

LOGICAL,  INTENT(IN)  :: cyclic
                          !   =T, then source data is cyclic
                          !   =F, then source data is non-cyclic

! Local arrays:---------------------------------------------------------
REAL         ::  t_lambda(points) !Local value of target longitude

INTEGER      ::  ixp1(points)     !Longitudinal index plus 1
INTEGER      ::  ix(points)       !Longitudinal index
INTEGER      ::  iy(points)       !Latitudinal index

! new logical arrays used for vectorizing the set up of
! longitudinal and latitudinal indexes





!    Local variables:---------------------------------------------------
REAL         ::  a                !Longitudinal weight
REAL         ::  b                !Latitudinal weight

INTEGER      ::  i
INTEGER      ::  j                !Loop indices

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! ----------------------------------------------------------------------

! Variables for divide-and-conquer search of latitude/longitude arrays
INTEGER :: iupr   ! Uppermost-point
INTEGER :: ilwr   ! Lowermost-point 
INTEGER :: imid   ! Mid-point 

!     1. Initialise arrays

!         1.1 Scale target longitude so that it falls between
!             LAMBDA_SRCE(1) and LAMBDA_SRCE(1)+360

IF (lhook) CALL dr_hook('H_INT_CO',zhook_in,zhook_handle)

DO i=1,points
  t_lambda(i)=MOD(((lambda_targ(i)-lambda_srce(1))+720.),360.)  &
              +lambda_srce(1)
END DO

IF(cyclic)THEN
  DO i=1,points
    ix(i)=0
    iy(i)=1
  END DO
ELSE
  DO i=1,points
    ix(i)=1
    iy(i)=1
  END DO
END IF



!  2. Calculate lat and lon index of bottom left hand corner of
!     source grid box enclosing each target point.

! Longitude



!$OMP PARALLEL DO SCHEDULE(STATIC)                                &
!$OMP& SHARED(points,points_lambda_srce,t_lambda,lambda_srce)     &
!$OMP& SHARED(points_phi_srce,phi_srce,phi_targ,ix,iy)            &
!$OMP& PRIVATE(i,iupr,ilwr,imid) DEFAULT(NONE)
DO i=1,points

! Divide and conquer should more efficient than a brute-force loop over 
! all i

  ilwr = 0                        !first point(-1 in case we are cyclic)
  iupr = points_lambda_srce + 1   !last point (+1 in case we are cyclic) 
  DO WHILE ( (iupr-ilwr) > 1 )
    imid = (ilwr+iupr)/2
    IF ( lambda_srce(imid) > t_lambda(i) ) THEN
      iupr = imid
    END IF
    IF ( lambda_srce(imid) <= t_lambda(i) ) THEN
      ilwr = imid
    END IF
  ENDDO

  ix(i) = ilwr

  ilwr = 1                    !first point
  iupr = points_phi_srce      !last point
  DO WHILE ( (iupr-ilwr) > 1 )
    imid = (ilwr+iupr)/2
    IF ( phi_srce(imid) > phi_targ(i) ) THEN
      iupr = imid
    END IF
    IF ( phi_srce(imid) <= phi_targ(i) ) THEN
      ilwr = imid
    END IF
  ENDDO

  iy(i) = ilwr

END DO
!$OMP END PARALLEL DO



!   3. Correct 1-D indices for wrap around etc and then calculate
!      2-D indices of bottom left and bottom right hand corner
!      of each grid box.

IF(cyclic)THEN
!     3.1 Cyclic case

  DO i=1,points

! Set index for cyclic wrap around (not sure this is ever met since
! t_lambda always runs from lambda_srce(1) -> lambda_srce(points_lambda_srce)
! but we shall keep it here in case).
    IF(ix(i) <  1)THEN
      ix(i)=points_lambda_srce
      t_lambda(i)=t_lambda(i)+360.
    END IF

! Set index for one sided difference if target point to north or
! south of source area.
    iy(i)=MAX(iy(i),1)
    iy(i)=MIN(iy(i),points_phi_srce-1)

! 2-D indices
    index_b_l(i)=ix(i)+(iy(i)-1)*points_lambda_srce
    index_b_r(i)=index_b_l(i)+1

! Correct for cyclic boundaries if target point outside source grid.

    ixp1(i)=ix(i)+1
    IF(ix(i) == points_lambda_srce)THEN
      index_b_r(i)=index_b_r(i)-points_lambda_srce
      ixp1(i)=ixp1(i)-points_lambda_srce
    END IF

  END DO

ELSE

!     3.2 Non cyclic case
  DO i=1,points

! Check that the nearest source grid point is really the nearest.  Earlier we
! made sure t_lambda is between lambda_srce(1) and lambda_srce(1)+360.  If we
! are a LAM the first boundary might be nearer so make sure t_lambda reflects
! this.
    IF (ix(i) == points_lambda_srce) THEN
      IF (ABS(t_lambda(i)-lambda_srce(1)-360.) < &
          t_lambda(i)-lambda_srce(ix(i))) THEN
        t_lambda(i) = t_lambda(i) - 360.
        ix(i) = 1
      END IF
    END IF

! Set index for one sided difference if outside source area
    ix(i)=MAX(ix(i),1)
    ix(i)=MIN(ix(i),points_lambda_srce-1)
    IF (ix(i) <  1) THEN ! IX(I) < 1 if POINTS_LAMBDA_SRCE = 1
      ix(i)=1
    END IF

    ixp1(i)=ix(i)+1
    ixp1(i)=MIN(ixp1(i),points_lambda_srce)

! Set index for one sided difference if outside source area
    iy(i)=MAX(iy(i),1)
    iy(i)=MIN(iy(i),points_phi_srce-1)
    IF (iy(i) <  1) THEN ! IY(I) < 1 if POINTS_PHI_SRCE = 1
      iy(i)=1
    END IF


! 2-D indices
    index_b_l(i)=ix(i)  +(iy(i)-1)*points_lambda_srce
    index_b_r(i)=ixp1(i)+(iy(i)-1)*points_lambda_srce

  END DO

END IF

!  4. Compute interpolation weights

DO i=1,points

! Calculate basic weights (equation 2.2)
  a=(MOD(360.+lambda_srce(ixp1(i))-lambda_srce(ix(i)),360.))
  IF(a /= 0.)THEN
! If t_lambda - lambda_source is negative then just copy last value across.
    a=(MAX(t_lambda(i)-lambda_srce(ix(i)),0.0))/a
  ELSE
    a=0.
  END IF

! If we only have 1 row then we need to make sure we can cope.
  b=ABS(phi_srce(iy(i))-phi_srce(MIN(iy(i)+1,points_phi_srce)))
  IF (b /= 0.0) THEN
    b = MAX(phi_targ(i)-phi_srce(iy(i)),0.0)/b
  ELSE
    b = 0.
  END IF

! Calculate bi-linear interpolation weights as per equation (2.1)

  weight_t_r(i)=a*b
  weight_b_l(i)=(1.-a)*(1.-b)
  weight_t_l(i)=(1.-a)*b
  weight_b_r(i)=a*(1.-b)

END DO

IF (lhook) CALL dr_hook('H_INT_CO',zhook_out,zhook_handle)
RETURN
END SUBROUTINE h_int_co
END MODULE h_int_co_mod
