! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

 FUNCTION chi(n)

! This power law is what Glenn Shutts observed in the CRM simulations
! (see Berner et al., 2009: J. Atmos. Sci, pp 603-626)
! Power of -1.54 based on practice at ECMWF
! e.g.: chi=n^^(2p+1); p = -1.27 => 2p+1 = -1.54

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim

 IMPLICIT NONE

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

 INTEGER :: n
 REAL    :: chi

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('CHI',zhook_in,zhook_handle)

 chi= n**(-1.54)

 IF (lhook) CALL dr_hook('CHI',zhook_out,zhook_handle)
 RETURN

END FUNCTION chi
