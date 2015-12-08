! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE rpassm_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE rpassm( a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot,      &
                   n, ifac, la,ierr )

!     performs one pass through data as part
!     of multiple real fft (fourier synthesis) routine

!     a is first real input vector
!         equivalence b(1) with a (la*inc1+1)
!     c is first real output vector
!         equivalence d(1) with c(ifac*la*inc2+1)
!     trigs is a precalculated list of sines & cosines
!     inc1 is the addressing increment for a
!     inc2 is the addressing increment for c
!     inc3 is the increment between input vectors a
!     inc4 is the increment between output vectors c
!     lot is the number of vectors
!     n is the length of the vectors
!     ifac is the current factor of n
!     la is the product of previous factors
!     ierr is an error indicator:
!              0 - pass completed without error
!              1 - lot greater than nblock
!              2 - ifac not catered for
!              3 - ifac only catered for if la=n/ifac



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

INTEGER :: n, lot, ifac, la, ierr, inc1, inc2, inc3, inc4
INTEGER :: m, iink, jink, ijump, kstop, ibad, ibase, igo                &
,          i, j, k, l, ijk, jbase, jump                                 &
,          ia, ib, ic, id, ie, if, ig, ih                               &
,          ja, jb, jc, jd, je, jf, jg, jh                               &
,          kb, kc, kd, ke, kf
REAL :: a(lot*(n+2)), b(lot*(n+2)), c(lot*(n+2)), d(lot*(n+2)), z       &
,       trigs(n), sin45, ssin36, ssin45, ssin60, ssin72, qqrt5          &
,       a0, a1, a2, a3, a4, a5, a6                                      &
,       b0, b1, b2, b3, b4, b5, b6                                      &
,       c1, c2, c3, c4, c5                                              &
,       s1, s2, s3, s4, s5                                              &
,       a10, a11, a20, a21, b10, b11, b20, b21

REAL :: sin36 = 0.587785252292473, sin72 = 0.951056516295154,           &
         qrt5 = 0.559016994374947, sin60 = 0.866025403784437

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('RPASSM',zhook_in,zhook_handle)

m=n/ifac
iink=la*inc1
jink=la*inc2
jump=(ifac-1)*jink
kstop=(n-ifac)/(2*ifac)

ibad=1
IF (lot > nblock) THEN
   ierr=ibad
   IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
   RETURN
   END IF
ibase=0
jbase=0
igo=ifac-1
IF (igo == 7) igo=6
ibad=2
ierr=ibad
IF (igo < 1.OR.igo > 6) THEN 
  IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
  RETURN
END IF
ibad=0
ierr=ibad
IF( igo == 1 ) THEN

!     coding for factor 2
!     -------------------

   ia=1
   ib=ia+(2*m-la)*inc1
   ja=1
   jb=ja+jink
   IF (la /= m) THEN
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)+a(ib+i)
            c(jb+j)=a(ia+i)-a(ib+i)
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ia=ia+iink
      iink=2*iink
      ib=ib-iink
      ibase=0
      jbase=jbase+jump
      jump=2*jump+jink
      IF (ia /= ib) THEN
         DO k=la,kstop,la
            kb=k+k
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            ibase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                 c(ja+j)=a(ia+i)+a(ib+i)
                 d(ja+j)=b(ia+i)-b(ib+i)
                 c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)+b(ib+i))
                 d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)+b(ib+i))
                 i=i+inc3
                 j=j+inc4
                 END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ia=ia+iink
            ib=ib-iink
            jbase=jbase+jump
            END DO
         IF (ia > ib) THEN 
           IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      ibase=0
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)
            c(jb+j)=-b(ia+i)
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
         RETURN
      END IF
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         c(ja+j)=2.0*(a(ia+i)+a(ib+i))
         c(jb+j)=2.0*(a(ia+i)-a(ib+i))
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 3
!     -------------------

ELSE IF (igo.eq.2) THEN
   ia=1
   ib=ia+(2*m-la)*inc1
   ic=ib
   ja=1
   jb=ja+jink
   jc=jb+jink
   IF (la /= m) THEN
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)+a(ib+i)
            c(jb+j)=(a(ia+i)-0.5*a(ib+i))-(sin60*(b(ib+i)))
            c(jc+j)=(a(ia+i)-0.5*a(ib+i))+(sin60*(b(ib+i)))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ia=ia+iink
      iink=2*iink
      ib=ib+iink
      ic=ic-iink
      jbase=jbase+jump
      jump=2*jump+jink
      IF (ia /= ic) THEN
         DO k=la,kstop,la
            kb=k+k
            kc=kb+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            ibase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                  c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                  d(ja+j)=b(ia+i)+(b(ib+i)-b(ic+i))
                  c(jb+j)=c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))           &
                                   -(sin60*(b(ib+i)+b(ic+i))))          &
                         -s1*((b(ia+i)-0.5*(b(ib+i)-b(ic+i)))           &
                                   +(sin60*(a(ib+i)-a(ic+i))))
                  d(jb+j)=s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))           &
                                   -(sin60*(b(ib+i)+b(ic+i))))          &
                         +c1*((b(ia+i)-0.5*(b(ib+i)-b(ic+i)))           &
                                   +(sin60*(a(ib+i)-a(ic+i))))
                  c(jc+j)=c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))           &
                                   +(sin60*(b(ib+i)+b(ic+i))))          &
                         -s2*((b(ia+i)-0.5*(b(ib+i)-b(ic+i)))           &
                                   -(sin60*(a(ib+i)-a(ic+i))))
                  d(jc+j)=s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))           &
                                   +(sin60*(b(ib+i)+b(ic+i))))          &
                         +c2*((b(ia+i)-0.5*(b(ib+i)-b(ic+i)))           &
                                   -(sin60*(a(ib+i)-a(ic+i))))
                  i=i+inc3
                  j=j+inc4
                  END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ia=ia+iink
            ib=ib+iink
            ic=ic-iink
            jbase=jbase+jump
            END DO
         IF (ia > ic) THEN 
           IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      ibase=0
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)+a(ib+i)
            c(jb+j)=(0.5*a(ia+i)-a(ib+i))-(sin60*b(ia+i))
            c(jc+j)=-(0.5*a(ia+i)-a(ib+i))-(sin60*b(ia+i))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
         RETURN
      END IF
   ssin60=2.0*sin60
   DO l=1,la
      i=ibase
      j=jbase

     DO ijk=1,lot
         c(ja+j)=2.0*(a(ia+i)+a(ib+i))
         c(jb+j)=(2.0*a(ia+i)-a(ib+i))-(ssin60*b(ib+i))
         c(jc+j)=(2.0*a(ia+i)-a(ib+i))+(ssin60*b(ib+i))
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 4
!     -------------------

ELSE IF ( igo == 3) THEN
   ia=1
   ib=ia+(2*m-la)*inc1
   ic=ib+2*m*inc1
   id=ib
   ja=1
   jb=ja+jink
   jc=jb+jink
   jd=jc+jink
   IF (la /= m) THEN
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=(a(ia+i)+a(ic+i))+a(ib+i)
            c(jb+j)=(a(ia+i)-a(ic+i))-b(ib+i)
            c(jc+j)=(a(ia+i)+a(ic+i))-a(ib+i)
            c(jd+j)=(a(ia+i)-a(ic+i))+b(ib+i)
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ia=ia+iink
      iink=2*iink
      ib=ib+iink
      ic=ic-iink
      id=id-iink
      jbase=jbase+jump
      jump=2*jump+jink
      IF (ib /= ic) THEN
         DO k=la,kstop,la
            kb=k+k
            kc=kb+kb
            kd=kc+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            c3=trigs(kd+1)
            s3=trigs(kd+2)
            ibase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                  c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
                  d(ja+j)=(b(ia+i)-b(ic+i))+(b(ib+i)-b(id+i))
                  c(jc+j)=c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))      &
                         -s2*((b(ia+i)-b(ic+i))-(b(ib+i)-b(id+i)))
                  d(jc+j)=s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))      &
                         +c2*((b(ia+i)-b(ic+i))-(b(ib+i)-b(id+i)))
                  c(jb+j)=c1*((a(ia+i)-a(ic+i))-(b(ib+i)+b(id+i)))      &
                         -s1*((b(ia+i)+b(ic+i))+(a(ib+i)-a(id+i)))
                  d(jb+j)=s1*((a(ia+i)-a(ic+i))-(b(ib+i)+b(id+i)))      &
                         +c1*((b(ia+i)+b(ic+i))+(a(ib+i)-a(id+i)))
                  c(jd+j)=c3*((a(ia+i)-a(ic+i))+(b(ib+i)+b(id+i)))      &
                         -s3*((b(ia+i)+b(ic+i))-(a(ib+i)-a(id+i)))
                  d(jd+j)=s3*((a(ia+i)-a(ic+i))+(b(ib+i)+b(id+i)))      &
                         +c3*((b(ia+i)+b(ic+i))-(a(ib+i)-a(id+i)))
                  i=i+inc3
                  j=j+inc4
                  END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ia=ia+iink
            ib=ib+iink
            ic=ic-iink
            id=id-iink
            jbase=jbase+jump
            END DO
         IF (ib > ic) THEN 
           IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      ibase=0
      sin45=SQRT(0.5)
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)+a(ib+i)
            c(jb+j)=sin45*((a(ia+i)-a(ib+i))-(b(ia+i)+b(ib+i)))
            c(jc+j)=b(ib+i)-b(ia+i)
            c(jd+j)=-sin45*((a(ia+i)-a(ib+i))+(b(ia+i)+b(ib+i)))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
         RETURN
      END IF
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         c(ja+j)=2.0*((a(ia+i)+a(ic+i))+a(ib+i))
         c(jb+j)=2.0*((a(ia+i)-a(ic+i))-b(ib+i))
         c(jc+j)=2.0*((a(ia+i)+a(ic+i))-a(ib+i))
         c(jd+j)=2.0*((a(ia+i)-a(ic+i))+b(ib+i))
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 5
!     -------------------

ELSE IF (igo == 4) THEN
   ia=1
   ib=ia+(2*m-la)*inc1
   ic=ib+2*m*inc1
   id=ic
   ie=ib
   ja=1
   jb=ja+jink
   jc=jb+jink
   jd=jc+jink
   je=jd+jink
   IF (la /= m) THEN
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
            c(jb+j)=((a(ia+i)-0.25*(a(ib+i)+a(ic+i)))                   &
                          +qrt5*(a(ib+i)-a(ic+i)))                      &
                          -(sin72*b(ib+i)+sin36*b(ic+i))
            c(jc+j)=((a(ia+i)-0.25*(a(ib+i)+a(ic+i)))                   &
                          -qrt5*(a(ib+i)-a(ic+i)))                      &
                          -(sin36*b(ib+i)-sin72*b(ic+i))
            c(jd+j)=((a(ia+i)-0.25*(a(ib+i)+a(ic+i)))                   &
                          -qrt5*(a(ib+i)-a(ic+i)))                      &
                          +(sin36*b(ib+i)-sin72*b(ic+i))
            c(je+j)=((a(ia+i)-0.25*(a(ib+i)+a(ic+i)))                   &
                          +qrt5*(a(ib+i)-a(ic+i)))                      &
                          +(sin72*b(ib+i)+sin36*b(ic+i))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ia=ia+iink
      iink=2*iink
      ib=ib+iink
      ic=ic+iink
      id=id-iink
      ie=ie-iink
      jbase=jbase+jump
      jump=2*jump+jink
      IF (ib /= id) THEN
         DO k=la,kstop,la
            kb=k+k
            kc=kb+kb
            kd=kc+kb
            ke=kd+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            c3=trigs(kd+1)
            s3=trigs(kd+2)
            c4=trigs(ke+1)
            s4=trigs(ke+2)
            ibase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                  a10=(a(ia+i)-0.25*((a(ib+i)+a(ie+i))                  &
                                      +(a(ic+i)+a(id+i))))              &
                          +qrt5*((a(ib+i)+a(ie+i))-(a(ic+i)             &
                          +a(id+i)))
                  a20=(a(ia+i)-0.25*((a(ib+i)+a(ie+i))                  &
                                      +(a(ic+i)+a(id+i))))              &
                          -qrt5*((a(ib+i)+a(ie+i))-(a(ic+i)             &
                          +a(id+i)))
                  b10=(b(ia+i)-0.25*((b(ib+i)-b(ie+i))                  &
                                      +(b(ic+i)-b(id+i))))              &
                          +qrt5*((b(ib+i)-b(ie+i))-(b(ic+i)             &
                          -b(id+i)))
                  b20=(b(ia+i)-0.25*((b(ib+i)-b(ie+i))                  &
                                      +(b(ic+i)-b(id+i))))              &
                          -qrt5*((b(ib+i)-b(ie+i))-(b(ic+i)             &
                          -b(id+i)))
                  a11=sin72*(b(ib+i)+b(ie+i))                           &
                          +sin36*(b(ic+i)+b(id+i))
                  a21=sin36*(b(ib+i)+b(ie+i))                           &
                          -sin72*(b(ic+i)+b(id+i))
                  b11=sin72*(a(ib+i)-a(ie+i))                           &
                          +sin36*(a(ic+i)-a(id+i))
                  b21=sin36*(a(ib+i)-a(ie+i))                           &
                          -sin72*(a(ic+i)-a(id+i))

                  c(ja+j)=a(ia+i)+((a(ib+i)+a(ie+i))+(a(ic+i)           &
                          +a(id+i)))
                  d(ja+j)=b(ia+i)+((b(ib+i)-b(ie+i))+(b(ic+i)           &
                          -b(id+i)))
                  c(jb+j)=c1*(a10-a11)-s1*(b10+b11)
                  d(jb+j)=s1*(a10-a11)+c1*(b10+b11)
                  c(je+j)=c4*(a10+a11)-s4*(b10-b11)
                  d(je+j)=s4*(a10+a11)+c4*(b10-b11)
                  c(jc+j)=c2*(a20-a21)-s2*(b20+b21)
                  d(jc+j)=s2*(a20-a21)+c2*(b20+b21)
                  c(jd+j)=c3*(a20+a21)-s3*(b20-b21)
                  d(jd+j)=s3*(a20+a21)+c3*(b20-b21)

                  i=i+inc3
                  j=j+inc4
                  END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ia=ia+iink
            ib=ib+iink
            ic=ic+iink
            id=id-iink
            ie=ie-iink
            jbase=jbase+jump
            END DO
         IF (ib > id) THEN 
           IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      ibase=0
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=(a(ia+i)+a(ib+i))+a(ic+i)
            c(jb+j)=(qrt5*(a(ia+i)-a(ib+i))                             &
                      +(0.25*(a(ia+i)+a(ib+i))-a(ic+i)))                &
                      -(sin36*b(ia+i)+sin72*b(ib+i))
            c(je+j)=-(qrt5*(a(ia+i)-a(ib+i))                            &
                      +(0.25*(a(ia+i)+a(ib+i))-a(ic+i)))                &
                      -(sin36*b(ia+i)+sin72*b(ib+i))
            c(jc+j)=(qrt5*(a(ia+i)-a(ib+i))                             &
                      -(0.25*(a(ia+i)+a(ib+i))-a(ic+i)))                &
                      -(sin72*b(ia+i)-sin36*b(ib+i))
            c(jd+j)=-(qrt5*(a(ia+i)-a(ib+i))                            &
                      -(0.25*(a(ia+i)+a(ib+i))-a(ic+i)))                &
                      -(sin72*b(ia+i)-sin36*b(ib+i))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
         RETURN
      END IF

   qqrt5=2.0*qrt5
   ssin36=2.0*sin36
   ssin72=2.0*sin72
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         c(ja+j)=2.0*(a(ia+i)+(a(ib+i)+a(ic+i)))
         c(jb+j)=(2.0*(a(ia+i)-0.25*(a(ib+i)+a(ic+i)))                  &
                     +qqrt5*(a(ib+i)-a(ic+i)))-(ssin72*b(ib+i)          &
                     +ssin36*b(ic+i))
         c(jc+j)=(2.0*(a(ia+i)-0.25*(a(ib+i)+a(ic+i)))                  &
                     -qqrt5*(a(ib+i)-a(ic+i)))-(ssin36*b(ib+i)          &
                     -ssin72*b(ic+i))
         c(jd+j)=(2.0*(a(ia+i)-0.25*(a(ib+i)+a(ic+i)))                  &
                     -qqrt5*(a(ib+i)-a(ic+i)))+(ssin36*b(ib+i)          &
                     -ssin72*b(ic+i))
         c(je+j)=(2.0*(a(ia+i)-0.25*(a(ib+i)+a(ic+i)))                  &
                     +qqrt5*(a(ib+i)-a(ic+i)))+(ssin72*b(ib+i)          &
                     +ssin36*b(ic+i))
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 6
!     -------------------

ELSE IF (igo == 5) THEN
   ia=1
   ib=ia+(2*m-la)*inc1
   ic=ib+2*m*inc1
   id=ic+2*m*inc1
   ie=ic
   if=ib
   ja=1
   jb=ja+jink
   jc=jb+jink
   jd=jc+jink
   je=jd+jink
   jf=je+jink
   IF (la /= m) THEN
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=(a(ia+i)+a(id+i))+(a(ib+i)+a(ic+i))
            c(jd+j)=(a(ia+i)-a(id+i))-(a(ib+i)-a(ic+i))
            c(jb+j)=((a(ia+i)-a(id+i))+0.5*(a(ib+i)-a(ic+i)))           &
                           -(sin60*(b(ib+i)+b(ic+i)))
            c(jf+j)=((a(ia+i)-a(id+i))+0.5*(a(ib+i)-a(ic+i)))           &
                           +(sin60*(b(ib+i)+b(ic+i)))
            c(jc+j)=((a(ia+i)+a(id+i))-0.5*(a(ib+i)+a(ic+i)))           &
                           -(sin60*(b(ib+i)-b(ic+i)))
            c(je+j)=((a(ia+i)+a(id+i))-0.5*(a(ib+i)+a(ic+i)))           &
                           +(sin60*(b(ib+i)-b(ic+i)))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ia=ia+iink
      iink=2*iink
      ib=ib+iink
      ic=ic+iink
      id=id-iink
      ie=ie-iink
      if=if-iink
      jbase=jbase+jump
      jump=2*jump+jink
      IF (ic /= id) THEN
         DO k=la,kstop,la
            kb=k+k
            kc=kb+kb
            kd=kc+kb
            ke=kd+kb
            kf=ke+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            c3=trigs(kd+1)
            s3=trigs(kd+2)
            c4=trigs(ke+1)
            s4=trigs(ke+2)
            c5=trigs(kf+1)
            s5=trigs(kf+2)
            ibase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                  a11= (a(ie+i)+a(ib+i))+(a(ic+i)+a(IF+i))
                  a20=(a(ia+i)+a(id+i))-0.5*a11
                  a21=sin60*((a(ie+i)+a(ib+i))-(a(ic+i)+a(IF+i)))
                  b11= (b(ib+i)-b(ie+i))+(b(ic+i)-b(IF+i))
                  b20=(b(ia+i)-b(id+i))-0.5*b11
                  b21=sin60*((b(ib+i)-b(ie+i))-(b(ic+i)-b(IF+i)))

                  c(ja+j)=(a(ia+i)+a(id+i))+a11
                  d(ja+j)=(b(ia+i)-b(id+i))+b11
                  c(jc+j)=c2*(a20-b21)-s2*(b20+a21)
                  d(jc+j)=s2*(a20-b21)+c2*(b20+a21)
                  c(je+j)=c4*(a20+b21)-s4*(b20-a21)
                  d(je+j)=s4*(a20+b21)+c4*(b20-a21)

                  a11=(a(ie+i)-a(ib+i))+(a(ic+i)-a(IF+i))
                  b11=(b(ie+i)+b(ib+i))-(b(ic+i)+b(IF+i))
                  a20=(a(ia+i)-a(id+i))-0.5*a11
                  a21=sin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(IF+i)))
                  b20=(b(ia+i)+b(id+i))+0.5*b11
                  b21=sin60*((b(ie+i)+b(ib+i))+(b(ic+i)+b(IF+i)))

                  c(jd+j)=c3*((a(ia+i)-a(id+i))+a11)                    &
                         -s3*((b(ia+i)+b(id+i))-b11)
                  d(jd+j)=s3*((a(ia+i)-a(id+i))+a11)                    &
                         +c3*((b(ia+i)+b(id+i))-b11)
                  c(jb+j)=c1*(a20-b21)-s1*(b20-a21)
                  d(jb+j)=s1*(a20-b21)+c1*(b20-a21)
                  c(jf+j)=c5*(a20+b21)-s5*(b20+a21)
                  d(jf+j)=s5*(a20+b21)+c5*(b20+a21)

                  i=i+inc3
                  j=j+inc4
                  END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ia=ia+iink
            ib=ib+iink
            ic=ic+iink
            id=id-iink
            ie=ie-iink
            if=if-iink
            jbase=jbase+jump
            END DO
         IF (ic > id) THEN
           IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      ibase=0
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ib+i)+(a(ia+i)+a(ic+i))
            c(jd+j)=b(ib+i)-(b(ia+i)+b(ic+i))
            c(jb+j)=(sin60*(a(ia+i)-a(ic+i)))                           &
                       -(0.5*(b(ia+i)+b(ic+i))+b(ib+i))
            c(jf+j)=-(sin60*(a(ia+i)-a(ic+i)))                          &
                       -(0.5*(b(ia+i)+b(ic+i))+b(ib+i))
            c(jc+j)=sin60*(b(ic+i)-b(ia+i))                             &
                       +(0.5*(a(ia+i)+a(ic+i))-a(ib+i))
            c(je+j)=sin60*(b(ic+i)-b(ia+i))                             &
                       -(0.5*(a(ia+i)+a(ic+i))-a(ib+i))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
         RETURN
      END IF
   ssin60=2.0*sin60
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         c(ja+j)=(2.0*(a(ia+i)+a(id+i)))+(2.0*(a(ib+i)+a(ic+i)))
         c(jd+j)=(2.0*(a(ia+i)-a(id+i)))-(2.0*(a(ib+i)-a(ic+i)))
         c(jb+j)=(2.0*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i)))              &
                  -(ssin60*(b(ib+i)+b(ic+i)))
         c(jf+j)=(2.0*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i)))              &
                  +(ssin60*(b(ib+i)+b(ic+i)))
         c(jc+j)=(2.0*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i)))              &
                  -(ssin60*(b(ib+i)-b(ic+i)))
         c(je+j)=(2.0*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i)))              &
                  +(ssin60*(b(ib+i)-b(ic+i)))
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 8
!     -------------------

ELSE
   ibad=3
   IF (la == m) THEN
      ia=1
      ib=ia+la*inc1
      ic=ib+2*la*inc1
      id=ic+2*la*inc1
      ie=id+2*la*inc1
      ja=1
      jb=ja+jink
      jc=jb+jink
      jd=jc+jink
      je=jd+jink
      jf=je+jink
      jg=jf+jink
      jh=jg+jink
      ssin45=SQRT(2.0)
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=2.0*(((a(ia+i)+a(ie+i))+a(ic+i))+(a(ib+i)           &
                    +a(id+i)))
            c(je+j)=2.0*(((a(ia+i)+a(ie+i))+a(ic+i))-(a(ib+i)           &
                    +a(id+i)))
            c(jc+j)=2.0*(((a(ia+i)+a(ie+i))-a(ic+i))-(b(ib+i)           &
                    -b(id+i)))
            c(jg+j)=2.0*(((a(ia+i)+a(ie+i))-a(ic+i))+(b(ib+i)           &
                    -b(id+i)))
            c(jb+j)=2.0*((a(ia+i)-a(ie+i))-b(ic+i))                     &
                    +ssin45*((a(ib+i)-a(id+i))-(b(ib+i)+b(id+i)))
            c(jf+j)=2.0*((a(ia+i)-a(ie+i))-b(ic+i))                     &
                    -ssin45*((a(ib+i)-a(id+i))-(b(ib+i)+b(id+i)))
            c(jd+j)=2.0*((a(ia+i)-a(ie+i))+b(ic+i))                     &
                    -ssin45*((a(ib+i)-a(id+i))+(b(ib+i)+b(id+i)))
            c(jh+j)=2.0*((a(ia+i)-a(ie+i))+b(ic+i))                     &
                    +ssin45*((a(ib+i)-a(id+i))+(b(ib+i)+b(id+i)))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ibad=0
      END IF
   END IF

!     final return
!     ------------

ierr=ibad
IF (lhook) CALL dr_hook('RPASSM',zhook_out,zhook_handle)
RETURN
END SUBROUTINE rpassm

END MODULE rpassm_mod
