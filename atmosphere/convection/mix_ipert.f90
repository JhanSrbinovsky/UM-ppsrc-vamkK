! -----------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! -----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
SUBROUTINE mix_ipert_4a5a(npnts, nreal, nlev, nbl, ntml,                  &
                          index1,p_layer_boundaries,                      &
                          exner_layer_centres, dthbydt, dqbydt, flx_init, &
                          thpert, qpert)

! -----------------------------------------------------------------------
! Purpose
! To well-mix convective increments from the initial parcel 
! perturbation in shallow/deep convection throughout the 
! boundary layer.
!
! -----------------------------------------------------------------------

  USE PrintStatus_mod
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE


!----------------------------------------------------------------------
! Variables with intent in
!----------------------------------------------------------------------

  INTEGER, INTENT(in) :: npnts        ! Number of points in arrays
  INTEGER, INTENT(in) :: nreal        ! Number of points to actually do

  INTEGER, INTENT(in) :: nlev         ! Number of model levels

  INTEGER, INTENT(in) :: nbl          ! in number of model layers
                                      ! potentially in the boundary layer

  INTEGER, INTENT(in) :: ntml(npnts)  ! number of model layers for which
                                      ! increments are to be well-mixed.
  INTEGER, INTENT(IN) :: index1(nreal) ! Location of points to be processed

  REAL, INTENT(in)    :: flx_init(npnts) 
                                      ! the initial parcel mass flux
 
  REAL, INTENT(in)    :: thpert(npnts)
                                      ! the initial parcel temperature 
                                      ! perturbation

  REAL, INTENT(in)    :: qpert(npnts) ! in the initial parcel humidity
                                      ! perturbation                      

  REAL, INTENT(in)    :: p_layer_boundaries(npnts,0:nlev)
                                      ! pressure at layer boundaries  

  REAL, INTENT(in)    :: exner_layer_centres(npnts,0:nlev)
                                      ! exner pressure at layer centres

!----------------------------------------------------------------------
! variables with intent inout
!----------------------------------------------------------------------

  REAL, INTENT(inout) :: dthbydt(npnts,nlev)
                                      ! increment to potential
                                      ! temperature due to convection

  REAL, INTENT(inout) :: dqbydt(npnts,nlev)
                                      ! increment to mixing ratio
                                      ! due to convection

!----------------------------------------------------------------------
! Variables that are locally defined
!----------------------------------------------------------------------

  INTEGER :: i,k                      ! loop counters

  REAL    :: delpsum(nreal)           ! summation of model layer thicknesses
                                      ! with height. (pa)

  REAL    :: delpk(nreal,nbl)         ! difference in pressure across a
                                      ! layer (pa)

  REAL    :: delpexsum(nreal)         ! summation of model layer thicknesses
                                      ! multiplied by the exner pressure (pa).

  REAL    :: dthbydt_exdp(nreal)      ! increment to potential temperature
                                      ! due to intial perturbation at ntml
                                      ! multiplied by the layer thickness and 
                                      ! exner pressure

  REAL    :: dqbydt_dp(nreal)         ! increment to mixing ratio
                                      ! due to intial perturbation at ntml
                                      ! multiplied by the layer thickness

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Calculate the layer pressure thickness and sum up.
!----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('MIX_IPERT_4A5A',zhook_in,zhook_handle)

  DO i = 1, nreal
    delpk(i,1)   = p_layer_boundaries(index1(i),0) -                         &
                                  p_layer_boundaries(index1(i),1)
    delpsum(i)   = delpk(i,1)
    delpexsum(i) = delpk(i,1) * exner_layer_centres(index1(i),1)
  END DO

  DO k = 2, nbl
    DO i = 1, nreal
      IF (k  <=  ntml(index1(i))) THEN
        delpk(i,k)   = p_layer_boundaries(index1(i),k-1) -                  &
                              p_layer_boundaries(index1(i),k)
        delpsum(i)   = delpsum(i)   + delpk(i,k)
        delpexsum(i) = delpexsum(i) + delpk(i,k)                            &
                                     * exner_layer_centres(index1(i),k)
      ENDIF
    ENDDO
  ENDDO

!----------------------------------------------------------------------
! Calculate the potential temperature and mixing ratio increments due 
! to the initial perturbation multiplied be the appropriate thickness
! nb. the delpk(ntml) in the numerator and denominator cancel
!----------------------------------------------------------------------

  DO i = 1, nreal
    dthbydt_exdp(i) = -flx_init(index1(i)) * thpert(index1(i))             &
                             * exner_layer_centres(index1(i),ntml(index1(i)))
    dqbydt_dp(i)    = -flx_init(index1(i)) * qpert(index1(i))
    IF (printstatus >= prstatus_normal) THEN
      IF (flx_init(index1(i)) <= 0.0) THEN
        WRITE(6,'(a19,g26.18)') 'flx_init(i) <= 0.0 ', flx_init(index1(i))
      END IF
    END IF
  END DO

!----------------------------------------------------------------------
! Mix the increments due to initial perturbation throughout the 
! sub-cloud layer.
!----------------------------------------------------------------------

  DO k = 1, nbl
    DO i = 1, nreal
      IF (k  <=  ntml(index1(i))) THEN
        dthbydt(index1(i),k) = dthbydt(index1(i),k) +                     &
                                           dthbydt_exdp(i) / delpexsum(i)
        dqbydt(index1(i),k)  = dqbydt(index1(i),k)  +                     &
                                            dqbydt_dp(i)    / delpsum(i)
      ENDIF
    ENDDO
  ENDDO

  IF (lhook) CALL dr_hook('MIX_IPERT_4A5A',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE mix_ipert_4a5a
