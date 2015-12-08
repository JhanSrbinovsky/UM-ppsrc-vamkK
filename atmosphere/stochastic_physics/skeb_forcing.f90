! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE skeb_forcing_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE skeb_forcing( dpsidtc, dpsidts, n1, n2, nlim, dpsidt,        &
                         mu, nlat, timestep_number, ii)

! Uses cosine and sine Fourier coefficients that have been evolved in
! time using a Markov process to generate the Fourier series and
! convert this into grid space
! The Process uses some arcane functions ALF and EPS to create the
! Fourier series. This code is processed in parallel by distributing the
! latitude loop (in calling routine stph_skeb2) over multiple processors.
! Where possible, calculations on the first time-step are stored for
! future calls.

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE c_skeb2_mod, ONLY: nblock, nfacts
USE fourier_mod, ONLY: fourier
IMPLICIT NONE

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Factors nblock and nfacts needed for FFT arrays

 INTEGER :: timestep_number
 INTEGER :: n1,n2,nlim,lev
 REAL    :: dpsidtc(0:n2,n1:n2),dpsidts(0:n2,n1:n2)
 REAL    :: dpsidt(0:2*nlim+1)

 REAL    :: mu
 INTEGER :: i,m,m1,n,nn
 REAL    :: alf

! Latitude pointers
 INTEGER :: nlat,ii
 REAL,ALLOCATABLE,SAVE :: seqalf(:)

! N+2 equivalent to  2*nlim+2
 REAL    :: work((2*nlim+2)*nblock),trigs(2*nlim)
 INTEGER :: ifax(nfacts)

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('SKEB_FORCING',zhook_in,zhook_handle)

! Allocate variables (saved for subsequent calls)
 IF(.NOT.ALLOCATED(seqalf)) THEN
   ALLOCATE(seqalf(nlim*nlat*n2))
 END IF

 ifax=0
 dpsidt(:)= 0.0

 DO m1=0,2*nlim,2
   m= m1/2
   nn= MAX0(n1,m)
   DO n=nn,n2
     IF (timestep_number == 1) THEN
       ! DEPENDS ON:alf
       seqalf(ii)= alf(m,n,mu)
     END IF
     ! cosine coeffs.
     dpsidt(m1)=   dpsidt(m1)   + dpsidtc(m,n)*seqalf(ii)
     ! sine coeffs.
     dpsidt(m1+1)= dpsidt(m1+1) + dpsidts(m,n)*seqalf(ii)
     ii=ii+1
   END DO
 END DO

! Invert fourier series representation on each latitude
 CALL fourier( dpsidt, 2*nlim+2, trigs, ifax, 1, 2*nlim, 2*nlim,        &
              1, 1, work)
 IF (lhook) CALL dr_hook('SKEB_FORCING',zhook_out,zhook_handle)
 RETURN

END SUBROUTINE skeb_forcing

END MODULE skeb_forcing_mod
