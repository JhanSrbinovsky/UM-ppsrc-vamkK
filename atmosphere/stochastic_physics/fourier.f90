! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE fourier_mod

USE qpassm_mod, ONLY: qpassm
USE rpassm_mod, ONLY: rpassm
USE set99_mod,  ONLY: set99
IMPLICIT NONE

CONTAINS


 SUBROUTINE fourier( a, alen, trigs, ifax, inc, jump, n, lot, fsign,    &
                    work)

!     real transform of length n performed by removing redundant
!     operations from complex transform of length n

!     a is the array containing input & output data of size alen
!     work is an area of size (n+2)*min(lot,nblock)
!     trigs is a previously prepared list of trig function values
!     ifax is a previously prepared list of factors of n
!     inc is the increment within each data 'vector'
!         (e.g. inc=1 for consecutively stored data)
!     jump is the increment between the start of each data vector
!     n is the length of the data vectors
!     lot is the number of data vectors
!     fsign = +1 for transform from spectral to gridpoint
!           = -1 for transform from gridpoint to spectral

!     ordering of coefficients:
!         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
!         where b(0)=b(n/2)=0; (n+2) locations required

!     ordering of data:
!         x(0),x(1),x(2),...,x(n-1), 0 , 0 ; (n+2) locations required

!     n must be composed of factors 2,3 & 5 but does not have to be even

!     definition of transforms:
!     -------------------------

!     fsign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
!         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)

!     fsign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
!               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim

 USE c_skeb2_mod, ONLY: nblock, nfacts
 IMPLICIT NONE

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! nfacts and nblock values

 INTEGER, INTENT(IN) :: fsign, lot, n, inc, alen, jump
 INTEGER :: ifax(nfacts), nx, nfax, nblox, left, istart, nb, nvex, ia, i &
 ,         j, ii, jj, igo, la, k, ierr, ibase, jbase, ix, iz, ifac
 REAL    :: a(alen), work((n+2)*nblock), trigs(n)

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('FOURIER',zhook_in,zhook_handle)

IF (ifax(nfacts) /= n) CALL set99(trigs,ifax,n)
nfax=ifax(1)
nx=n+1
IF (mod(n,2) == 1) nx=n
nblox=1+(lot-1)/nblock
left=lot
IF (fsign == 1) THEN

!     fsign=+1, spectral to gridpoint transform
!     -----------------------------------------

      istart=1
      DO nb=1,nblox
         IF (left <= nblock) THEN
               nvex=left
            ELSE
               nvex=nblock
            END IF
         left=left-nvex
         ia=istart
         i=istart

         DO j=1,nvex
            a(i+inc)=0.5*a(i)
            i=i+jump
            END DO
         IF (mod(n,2) == 0) THEN
            i=istart+n*inc

            DO j=1,nvex
               a(i)=0.5*a(i)
               i=i+jump
               END DO
            END IF
         ia=istart+inc
         la=1
         igo=+1
         DO k=1,nfax
            ifac=ifax(k+1)
            ierr=-1
            IF (igo /= -1) THEN
              CALL rpassm(a(ia),a(ia+la*inc),work(1),                   &
                          work(ifac*la+1),trigs,inc,1,jump,nx,          &
                          nvex,n,ifac,la,ierr)
            ELSE
                CALL rpassm(work(1),work(la+1),a(ia),                   &
                            a(ia+ifac*la*inc),trigs,1,inc,nx,           &
                            jump,nvex,n,ifac,la,ierr)
            END IF
            IF (ierr /= 0) GO TO 500
            la=ifac*la
            igo=-igo
            ia=istart
            END DO

!     if necessary, copy results back to a
!     ------------------------------------

         IF (mod(nfax,2) == 1) THEN
            ibase=1
            jbase=ia
            DO jj=1,nvex
               i=ibase
               j=jbase
               DO ii=1,n
                  a(j)=work(i)
                  i=i+1
                  j=j+inc
                  END DO
               ibase=ibase+nx
               jbase=jbase+jump
               END DO
            END IF

!     fill in zeros at end
!     --------------------

         ix=istart+n*inc
         DO j=1,nvex
            a(ix)=0.0
            a(ix+inc)=0.0
            ix=ix+jump
            END DO

         istart=istart+nvex*jump
         END DO
         IF (lhook) CALL dr_hook('FOURIER',zhook_out,zhook_handle)
         RETURN

!     fsign=-1, gridpoint to spectral transform
!     -----------------------------------------

   ELSE
      istart=1
      DO nb=1,nblox
         IF (left <= nblock) THEN
               nvex=left
            ELSE
               nvex=nblock
            END IF
         left=left-nvex
         ia=istart
         la=n
         igo=+1

         DO k=1,nfax
            ifac=ifax(nfax+2-k)
            la=la/ifac
            ierr=-1
            IF (igo == 1) THEN
                  CALL qpassm(a(ia),a(ia+ifac*la*inc),work(1),          &
                              work(la+1),trigs,inc,1,jump,nx,           &
                              nvex,n,ifac,la,ierr)
               ELSE
                  CALL qpassm(work(1),work(ifac*la+1),a(ia),            &
                              a(ia+la*inc),trigs,1,inc,nx,              &
                              jump,nvex,n,ifac,la,ierr)
               END IF
            IF (ierr /= 0) GO TO 500
            igo=-igo
            ia=istart+inc
            END DO

!     if necessary, copy results back to a
!     ------------------------------------

         IF (mod(nfax,2) == 1) THEN
            ibase=1
            jbase=ia
            DO jj=1,nvex
               i=ibase
               j=jbase
               DO ii=1,n
                  a(j)=work(i)
                  i=i+1
                  j=j+inc
                  END DO
               ibase=ibase+nx
               jbase=jbase+jump
               END DO
            END IF

!     shift a(0) & fill in zero imag parts
!     ------------------------------------

         ix=istart

         DO j=1,nvex
            a(ix)=a(ix+inc)
            a(ix+inc)=0.0
            ix=ix+jump
            END DO
         IF (mod(n,2) == 0) THEN
            iz=istart+(n+1)*inc
            DO j=1,nvex
               a(iz)=0.0
               iz=iz+jump
               END DO
            END IF

         istart=istart+nvex*jump
         END DO
         IF (lhook) CALL dr_hook('FOURIER',zhook_out,zhook_handle)
         RETURN
   END IF

!     error messages
!     --------------

 500 CONTINUE
 IF( ierr == 1 ) THEN
      WRITE(6,'(A,I4,A)') 'vector length = ', nvex, ', greater than NBLOCK'
 ELSE IF( ierr == 2 ) THEN
      WRITE(6,'(A,I3,A)') 'factor = ', ifac, ' NOT catered for'
 ELSE
      WRITE(6,'(A,I3,A)') 'factor = ', ifac, ' ONLY catered for if la*ifac=n'
 END IF
 IF (lhook) CALL dr_hook('FOURIER',zhook_out,zhook_handle)
 RETURN

END SUBROUTINE fourier
END MODULE fourier_mod
