! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE backscatter_spectrum_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE backscatter_spectrum( gspect, nlim, n1, n2, tau,       &
                                timestep, tot_backscat, alpha)

 USE earth_constants_mod, ONLY: earth_radius

 USE conversions_mod, ONLY: pi

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim
 USE PrintStatus_mod
 IMPLICIT NONE
 
   ! This include contains the model PrintStatus

 INTEGER,INTENT(IN) :: nlim, n1, n2
 REAL, INTENT(IN)   :: tau, timestep, tot_backscat, alpha
 REAL, INTENT(OUT)  :: gspect(nlim)

!     Local variables
 REAL    :: gamman, p
 INTEGER :: n        ! Loop counter

!     Functions
 REAL :: chi         ! function

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('BACKSCATTER_SPECTRUM',zhook_in,zhook_handle)

! This power law is what Glenn Shutts observed in the CRM simulations
! (see Berner et al., 2009: J. Atmos. Sci, pp 603-626)
! Power of -1.54 based on practice at ECMWF
! e.g.: chi=n^^(2p+1); p = -1.27 => 2p+1 = -1.54
 p=-1.27

 gamman= 0.0
 DO n = n1, n2
! DEPENDS ON: chi
   gamman = gamman + (n+1)*(2*n+1)*chi(n)
 END DO
 gamman = gamman/alpha

 DO n = 1, nlim
   gspect(n) = 0.0
 END DO

! g(n) = b * n^p
! b = SQRT([4 * Pi * a^2 * tot_backscat]/[timestep * vz * GammaN]
! vz = variance of random numbers [-0.5; 0.5] = 1/12
!      tested by 1000 cases of random arrays of size = 1e9
! note: (4/vz) is pre-calculated = 48
 DO n = n1, n2
   gspect(n) = earth_radius * SQRT(48.*pi*tot_backscat/            &
              (timestep*gamman)) * n**p
 END DO

 IF  (printstatus  >  prstatus_normal) THEN
   WRITE(6,'("N (nlim=",I5," varies from N1 to N2 =",2I5)') nlim,n1,n2
   WRITE(6,'("GSPECT: wavenumber-dependent power-spectrum amplitude")')
   WRITE(6,'(12ES12.4)') gspect
 END IF
 IF (lhook) CALL dr_hook('BACKSCATTER_SPECTRUM',zhook_out,zhook_handle)
 RETURN

END SUBROUTINE backscatter_spectrum
END MODULE backscatter_spectrum_mod
