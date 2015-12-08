! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Convective cloud microphysics routine
!
MODULE con_rad_6a_mod

IMPLICIT NONE

! Description:
!   Calculates convective cloud top, base and amount
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CONTAINS

! Subroutine Interface:

      SUBROUTINE con_rad_6a (k, npnts, start_lev,                        &
                             ccwk, ccwkp1, tmp_ccwkp1, flxkp1, delpkp1,  &
                             l_q_interact, bterm,                        &
                             tcw, cclwp, lcca, lcbase, lctop,            &
                             iccb, icct, cca)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      
      USE earth_constants_mod, ONLY : g 
      
      IMPLICIT NONE



!----------------------------------------------------------------------
! Arguments with INTENT(IN)
!----------------------------------------------------------------------
INTEGER,INTENT(IN) :: k               ! Present model layer number
INTEGER,INTENT(IN) :: npnts           ! Vector length
INTEGER,INTENT(IN) :: start_lev(npnts)! Level at which convection initiated

REAL,INTENT(IN) :: ccwk(npnts)      ! Total condensate in level k (kg/kg)
REAL,INTENT(IN) :: ccwkp1(npnts)    ! Total condensate in level k+1 (kg/kg)
REAL,INTENT(IN) :: tmp_ccwkp1(npnts)! Total condensate in level k+1 (kg/kg)
                                    ! before precipitation
REAL,INTENT(IN) :: flxkp1(npnts)    ! Parcel mass flux in layer k+1 (Pa/s)
REAL,INTENT(IN) :: delpkp1(npnts)   ! pressure difference across layer k+1 (Pa)

LOGICAL,INTENT(IN) :: l_q_interact  ! True if PC2 is switched on
LOGICAL,INTENT(IN) :: bterm(npnts)  ! Mask for parcels which terminate 
                                    ! in layer k+1

! ----------------------------------------------------------------------
! Arguments with INTENT(IN/OUT)
! ----------------------------------------------------------------------
REAL,INTENT(INOUT) :: tcw(npnts)    ! Total condensed water (kg/m**2/s)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: cclwp(npnts)  ! Condensed water path (kg/m**2)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: lcca(npnts)   ! Lowest conv. cloud amount (%)

INTEGER,INTENT(INOUT) :: lcbase(npnts)! Lowest conv. cloud base level
INTEGER,INTENT(INOUT) :: lctop(npnts) ! Lowest conv. cloud top level

!----------------------------------------------------------------------
! Arguments with INTENT(OUT)
!----------------------------------------------------------------------
INTEGER,INTENT(OUT) :: iccb(npnts)  ! convective cloud base_level
INTEGER,INTENT(OUT) :: icct(npnts)  ! convective cloud top level

REAL,INTENT(OUT) :: cca(npnts)      ! convective cloud amount (%)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------
INTEGER :: i          ! Loop counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
!  Calculate cloud base and lowest cloud base
!  when cloud base set zero total condensed water
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CON_RAD_6A',zhook_in,zhook_handle)
      DO i=1, npnts
        IF ( ccwk(I)  <=  0.0 .AND. tmp_ccwkp1(I)  >   0.0 ) THEN
!       Assuming initial parcel condensate is zero
          iccb(i)   = k+1
          cclwp(i)  = 0.0

          IF ( lcbase(i)  ==  0 ) THEN
!         Lowest cloud base
            lcbase(i) = k+1
          END IF  ! If_lcbase_1

        ELSE IF ( l_q_interact .AND. k  ==  start_lev(i) ) THEN

          ! Do not know why this is done,
          ! but it lowers the ccb so that it is no longer consistent
          ! with ccw. Needs looking into.

          ! Non-zero initial parcel condensate (initialized to environment)
          iccb(i)  = k
          cclwp(i) = 0.0
          tcw(i)   = 0.0

          IF ( lcbase(i)  ==  0 ) THEN
!         Lowest cloud base
            lcbase(i) = k
          END IF  ! If_lcbase_2

        END IF

! ----------------------------------------------------------------------
!  calculate cloud top and lowest cloud top
! ----------------------------------------------------------------------
        IF (bterm(i)) THEN
          icct(i)  = k+1
        END IF

        IF (bterm(i) .AND.  lctop(i) == 0 ) THEN
          lctop(i) = k+1
        END IF

        IF ( flxkp1(i)  >   0.0) THEN

!---------------------------------------------------------------------
! Sum total condensed water per second - assumes that the initial
! convective layer is unsaturated.
! Uses the CCW before precipitation
!---------------------------------------------------------------------
          tcw(i)   = tcw(i)   + flxkp1(i) * tmp_ccwkp1(i) / g

!---------------------------------------------------------------------
! Sum conv condensed water path - assumes that the initial
! convective layer is unsaturated
! Uses the CCW after precipitation
!---------------------------------------------------------------------
          cclwp(i) = cclwp(i) + ccwkp1(i) * delpkp1(i) / g

        END IF

!---------------------------------------------------------------------
! calculate convective cloud amount if convection terminates in
! layer k and total condensed water path over a time step
!
! UM documentation paper p27
! section (9), equation (37)
!---------------------------------------------------------------------
        IF( bterm(i) .AND. tcw(i) >  0.0 ) THEN

          IF ( tcw(i)  <   2.002e-6 ) tcw(i) = 2.002E-6

          cca(i) = 0.7873 + 0.06 * LOG(tcw(i))
          IF (cca(i)  >   1.0) cca(i) = 1.0

          IF (lcca(i) <= 0.0) THEN
            lcca(i) = 0.7873 + 0.06 * LOG(tcw(i))
            IF (lcca(i)  >   1.0) lcca(i) = 1.0
          END IF

          tcw(i) = 0.0

        END IF
      END DO ! I loop over NPNTS

      IF (lhook) CALL dr_hook('CON_RAD_6A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE con_rad_6a
END MODULE con_rad_6a_mod
