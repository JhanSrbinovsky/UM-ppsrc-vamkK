! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Mixes the subcloud temperature and humidity increments
!
MODULE mix_ipert_6a_mod

IMPLICIT NONE

! Description:
! Mixes the convective increments from the initial parcel 
! perturbation in shallow/deep convection throughout the 
! boundary layer.
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

  SUBROUTINE mix_ipert_6a(npnts, nlev, nbl, ntml, p_layer_boundaries, &
                          exner_layer_centres, dthbydt, dqbydt, flx_init, &
                          thpert, qpert)

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

!----------------------------------------------------------------------
! Variables with intent in
!----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: npnts        ! Number of points

  INTEGER, INTENT(IN) :: nlev         ! Number of model levels

  INTEGER, INTENT(IN) :: nbl          ! in number of model layers
                                      ! potentially in the boundary layer

  INTEGER, INTENT(IN) :: ntml(npnts)  ! number of model layers for which
                                      ! increments are to be well-mixed.

  REAL, INTENT(IN)    :: flx_init(npnts) 
                                      ! the initial parcel mass flux
 
  REAL, INTENT(IN)    :: thpert(npnts)
                                      ! the initial parcel temperature 
                                      ! perturbation

  REAL, INTENT(IN)    :: qpert(npnts) ! in the initial parcel humidity
                                      ! perturbation                      

  REAL, INTENT(IN)    :: p_layer_boundaries(npnts,0:nlev)
                                      ! pressure at layer boundaries  

  REAL, INTENT(IN)    :: exner_layer_centres(npnts,0:nlev)
                                      ! exner pressure at layer centres

!----------------------------------------------------------------------
! variables with intent inout
!----------------------------------------------------------------------

  REAL, INTENT(INOUT) :: dthbydt(npnts,nlev)
                                      ! increment to potential
                                      ! temperature due to convection

  REAL, INTENT(INOUT) :: dqbydt(npnts,nlev)
                                      ! increment to mixing ratio
                                      ! due to convection

!----------------------------------------------------------------------
! Variables that are locally defined
!----------------------------------------------------------------------

  INTEGER :: i,k                      ! loop counters

  REAL    :: delpsum(npnts)           ! summation of model layer thicknesses
                                      ! with height. (pa)

  REAL    :: delpk(npnts,nbl)         ! difference in pressure across a
                                      ! layer (pa)

  REAL    :: delpexsum(npnts)         ! summation of model layer thicknesses
                                      ! multiplied by the exner pressure (pa).

  REAL    :: dthbydt_exdp(npnts)      ! increment to potential temperature
                                      ! due to intial perturbation at ntml
                                      ! multiplied by the layer thickness and 
                                      ! exner pressure

  REAL    :: dqbydt_dp(npnts)         ! increment to mixing ratio
                                      ! due to intial perturbation at ntml
                                      ! multiplied by the layer thickness

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Calculate the layer pressure thickness and sum up.
!----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('MIX_IPERT_6A',zhook_in,zhook_handle)

  DO i = 1, npnts
    delpk(i,1)   = p_layer_boundaries(i,0) - p_layer_boundaries(i,1)
    delpsum(i)   = delpk(i,1)
    delpexsum(i) = delpk(i,1) * exner_layer_centres(i,1)
  END DO

  DO k = 2, nbl
    DO i = 1, npnts
      IF (k  <=  ntml(i)) THEN
        delpk(i,k)   = p_layer_boundaries(i,k-1) - p_layer_boundaries(i,k)
        delpsum(i)   = delpsum(i)   + delpk(i,k)
        delpexsum(i) = delpexsum(i) + delpk(i,k) * exner_layer_centres(i,k)
      END IF
    END DO
  END DO

!----------------------------------------------------------------------
! Calculate the potential temperature and mixing ratio increments due 
! to the initial perturbation multiplied be the appropriate thickness
! nb. the delpk(ntml) in the numerator and denominator cancel
!----------------------------------------------------------------------

  DO i = 1, npnts
    dthbydt_exdp(i) = -flx_init(i) * thpert(i) * exner_layer_centres(i,ntml(i))
    dqbydt_dp(i)    = -flx_init(i) * qpert(i)
  END DO

!----------------------------------------------------------------------
! Mix the increments due to initial perturbation throughout the 
! sub-cloud layer.
!----------------------------------------------------------------------

  DO k = 1, nbl
    DO i = 1, npnts
      IF (k  <=  ntml(i)) THEN
        dthbydt(i,k) = dthbydt(i,k) + dthbydt_exdp(i) / delpexsum(i)
        dqbydt(i,k)  = dqbydt(i,k)  + dqbydt_dp(i)    / delpsum(i)
      END IF
    END DO
  END DO

  IF (lhook) CALL dr_hook('MIX_IPERT_6A',zhook_out,zhook_handle)
  RETURN

  END SUBROUTINE mix_ipert_6a

END MODULE mix_ipert_6a_mod
