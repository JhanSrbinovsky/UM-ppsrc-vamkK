! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE qpassm_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE qpassm( a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot,      &
                   n, ifac, la,ierr )

!     performs one pass through data as part
!     of multiple real fft (fourier analysis) routine

!     a is first real input vector
!         equivalence b(1) with a(ifac*la*inc1+1)
!     c is first real output vector
!         equivalence d(1) with c(la*inc2+1)
!     trigs is a precalculated list of sines & cosines
!     inc1 is the addressing increment for a
!     inc2 is the addressing increment for c
!     inc3 is the increment between input vectors a
!     inc4 is the increment between output vectors c
!     lot is the number of vectors
!     n is the length of the vectors
!     ifac is the current factor of n
!     la = n/(product of factors used so far)
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
INTEGER :: i, j, k, l, m, iink, jink, ijump, kstop, ibad, ibase, igo    &
,          ia, ib, ic, id, ie, if, ig, ih, ijk                          &
,          ja, jb, jc, jd, je, jf, jbase                                &
,          kb, kc, kd, ke, kf
REAL :: a(lot*(n+2)), b(lot*(n+2)), c(lot*(n+2)), d(lot*(n+2)), z       &
,       trigs(n), sin45, zsin36, zsin45, zsin60, zsin72, zqrt5          &
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

IF (lhook) CALL dr_hook('QPASSM',zhook_in,zhook_handle)

m=n/ifac
iink=la*inc1
jink=la*inc2
ijump=(ifac-1)*iink
kstop=(n-ifac)/(2*ifac)

ibad=1
ierr=ibad
IF (lot > nblock) THEN 
  IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
  RETURN
END IF
ibase=0
jbase=0
igo=ifac-1
IF (igo == 7) igo=6
ibad=2
ierr=ibad
IF (igo < 1.OR.igo > 6) THEN 
  IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
  RETURN
END IF
ibad=0
ierr=0
IF( igo == 1 ) THEN

!     coding for factor 2
!     -------------------

   ia=1
   ib=ia+iink
   ja=1
   jb=ja+(2*m-la)*inc2
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
      ja=ja+jink
      jink=2*jink
      jb=jb-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      IF (ja /= jb) THEN
         DO k=la,kstop,la
            kb=k+k
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            jbase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                  c(ja+j)=a(ia+i)+(c1*a(ib+i)+s1*b(ib+i))
                  c(jb+j)=a(ia+i)-(c1*a(ib+i)+s1*b(ib+i))
                  d(ja+j)=(c1*b(ib+i)-s1*a(ib+i))+b(ia+i)
                  d(jb+j)=(c1*b(ib+i)-s1*a(ib+i))-b(ia+i)
                  i=i+inc3
                  j=j+inc4
                  END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ibase=ibase+ijump
            ja=ja+jink
            jb=jb-jink
            END DO
         IF (ja > jb) THEN 
           IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      jbase=0
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)
            d(ja+j)=-a(ib+i)
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
         RETURN
      END IF
   z=1.0/FLOAT(n)
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         c(ja+j)=z*(a(ia+i)+a(ib+i))
         c(jb+j)=z*(a(ia+i)-a(ib+i))
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 3
!     -------------------

ELSE IF (igo == 2) THEN
   ia=1
   ib=ia+iink
   ic=ib+iink
   ja=1
   jb=ja+(2*m-la)*inc2
   jc=jb
   IF (la /= m) THEN
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
            c(jb+j)=a(ia+i)-0.5*(a(ib+i)+a(ic+i))
            d(jb+j)=sin60*(a(ic+i)-a(ib+i))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ja=ja+jink
      jink=2*jink
      jb=jb+jink
      jc=jc-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      IF (ja /= jc) THEN
         DO k=la,kstop,la
            kb=k+k
            kc=kb+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            jbase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                  a1=(c1*a(ib+i)+s1*b(ib+i))+(c2*a(ic+i)+s2*            &
                      b(ic+i))
                  b1=(c1*b(ib+i)-s1*a(ib+i))+(c2*b(ic+i)-s2*            &
                      a(ic+i))
                  a2=a(ia+i)-0.5*a1
                  b2=b(ia+i)-0.5*b1
                  a3=sin60*((c1*a(ib+i)+s1*b(ib+i))-(c2*                &
                             a(ic+i)                                    &
                            +s2*b(ic+i)))
                  b3=sin60*((c1*b(ib+i)-s1*a(ib+i))-(c2*b(ic+i)         &
                            -s2*a(ic+i)))
                  c(ja+j)=a(ia+i)+a1
                  d(ja+j)=b(ia+i)+b1
                  c(jb+j)=a2+b3
                  d(jb+j)=b2-a3
                  c(jc+j)=a2-b3
                  d(jc+j)=-(b2+a3)
                  i=i+inc3
                  j=j+inc4
                  END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ibase=ibase+ijump
            ja=ja+jink
            jb=jb+jink
            jc=jc-jink
            END DO
         IF (ja > jc) THEN 
           IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      jbase=0
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)+0.5*(a(ib+i)-a(ic+i))
            d(ja+j)=-sin60*(a(ib+i)+a(ic+i))
            c(jb+j)=a(ia+i)-(a(ib+i)-a(ic+i))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
         RETURN
      END IF
   z=1.0/FLOAT(n)
   zsin60=z*sin60
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         c(ja+j)=z*(a(ia+i)+(a(ib+i)+a(ic+i)))
         c(jb+j)=z*(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
         d(jb+j)=zsin60*(a(ic+i)-a(ib+i))
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 4
!     -------------------

ELSE IF (igo == 3) THEN
   ia=1
   ib=ia+iink
   ic=ib+iink
   id=ic+iink
   ja=1
   jb=ja+(2*m-la)*inc2
   jc=jb+2*m*inc2
   jd=jb
   IF (la /= m) THEN
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
            c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
            c(jb+j)=a(ia+i)-a(ic+i)
            d(jb+j)=a(id+i)-a(ib+i)
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ja=ja+jink
      jink=2*jink
      jb=jb+jink
      jc=jc-jink
      jd=jd-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      IF (jb /= jc) THEN
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
            jbase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                  a0=a(ia+i)+(c2*a(ic+i)+s2*b(ic+i))
                  a2=a(ia+i)-(c2*a(ic+i)+s2*b(ic+i))
                  a1=(c1*a(ib+i)+s1*b(ib+i))+(c3*a(id+i)+s3             &
                      *b(id+i))
                  a3=(c1*a(ib+i)+s1*b(ib+i))-(c3*a(id+i)+s3             &
                      *b(id+i))
                  b0=b(ia+i)+(c2*b(ic+i)-s2*a(ic+i))
                  b2=b(ia+i)-(c2*b(ic+i)-s2*a(ic+i))
                  b1=(c1*b(ib+i)-s1*a(ib+i))+(c3*b(id+i)-s3             &
                      *a(id+i))
                  b3=(c1*b(ib+i)-s1*a(ib+i))-(c3*b(id+i)-s3             &
                      *a(id+i))
                  c(ja+j)=a0+a1
                  c(jc+j)=a0-a1
                  d(ja+j)=b0+b1
                  d(jc+j)=b1-b0
                  c(jb+j)=a2+b3
                  c(jd+j)=a2-b3
                  d(jb+j)=b2-a3
                  d(jd+j)=-(b2+a3)
                  i=i+inc3
                  j=j+inc4
                  END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ibase=ibase+ijump
            ja=ja+jink
            jb=jb+jink
            jc=jc-jink
            jd=jd-jink
            END DO
         IF (jb > jc) THEN 
           IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      sin45=SQRT(0.5)
      jbase=0
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=a(ia+i)+sin45*(a(ib+i)-a(id+i))
            c(jb+j)=a(ia+i)-sin45*(a(ib+i)-a(id+i))
            d(ja+j)=-a(ic+i)-sin45*(a(ib+i)+a(id+i))
            d(jb+j)=a(ic+i)-sin45*(a(ib+i)+a(id+i))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
         RETURN
      END IF
   z=1.0/FLOAT(n)
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         c(ja+j)=z*((a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i)))
         c(jc+j)=z*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
         c(jb+j)=z*(a(ia+i)-a(ic+i))
         d(jb+j)=z*(a(id+i)-a(ib+i))
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 5
!     -------------------

ELSE IF (igo == 4) THEN
   ia=1
   ib=ia+iink
   ic=ib+iink
   id=ic+iink
   ie=id+iink
   ja=1
   jb=ja+(2*m-la)*inc2
   jc=jb+2*m*inc2
   jd=jc
   je=jb
   IF (la /= m) THEN
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            a1=a(ib+i)+a(ie+i)
            a3=a(ib+i)-a(ie+i)
            a2=a(ic+i)+a(id+i)
            a4=a(ic+i)-a(id+i)
            a5=a(ia+i)-0.25*(a1+a2)
            a6=qrt5*(a1-a2)
            c(ja+j)=a(ia+i)+(a1+a2)
            c(jb+j)=a5+a6
            c(jc+j)=a5-a6
            d(jb+j)=-sin72*a3-sin36*a4
            d(jc+j)=-sin36*a3+sin72*a4
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ja=ja+jink
      jink=2*jink
      jb=jb+jink
      jc=jc+jink
      jd=jd-jink
      je=je-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      IF (jb /= jd) THEN
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
            jbase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                  a1=(c1*a(ib+i)+s1*b(ib+i))+(c4*a(ie+i)+s4             &
                      *b(ie+i))
                  a3=(c1*a(ib+i)+s1*b(ib+i))-(c4*a(ie+i)+s4             &
                      *b(ie+i))
                  a2=(c2*a(ic+i)+s2*b(ic+i))+(c3*a(id+i)+s3             &
                      *b(id+i))
                  a4=(c2*a(ic+i)+s2*b(ic+i))-(c3*a(id+i)+s3             &
                      *b(id+i))
                  b1=(c1*b(ib+i)-s1*a(ib+i))+(c4*b(ie+i)-s4             &
                      *a(ie+i))
                  b3=(c1*b(ib+i)-s1*a(ib+i))-(c4*b(ie+i)-s4             &
                      *a(ie+i))
                  b2=(c2*b(ic+i)-s2*a(ic+i))+(c3*b(id+i)-s3             &
                      *a(id+i))
                  b4=(c2*b(ic+i)-s2*a(ic+i))-(c3*b(id+i)-s3             &
                      *a(id+i))
                  a5=a(ia+i)-0.25*(a1+a2)
                  a6=qrt5*(a1-a2)
                  b5=b(ia+i)-0.25*(b1+b2)
                  b6=qrt5*(b1-b2)
                  a10=a5+a6
                  a20=a5-a6
                  b10=b5+b6
                  b20=b5-b6
                  a11=sin72*b3+sin36*b4
                  a21=sin36*b3-sin72*b4
                  b11=sin72*a3+sin36*a4
                  b21=sin36*a3-sin72*a4
                  c(ja+j)=a(ia+i)+(a1+a2)
                  c(jb+j)=a10+a11
                  c(je+j)=a10-a11
                  c(jc+j)=a20+a21
                  c(jd+j)=a20-a21
                  d(ja+j)=b(ia+i)+(b1+b2)
                  d(jb+j)=b10-b11
                  d(je+j)=-(b10+b11)
                  d(jc+j)=b20-b21
                  d(jd+j)=-(b20+b21)
                  i=i+inc3
                  j=j+inc4
                  END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ibase=ibase+ijump
            ja=ja+jink
            jb=jb+jink
            jc=jc+jink
            jd=jd-jink
            je=je-jink
            END DO
         IF (jb > jd) THEN 
           IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      jbase=0
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            a1=a(ib+i)+a(ie+i)
            a3=a(ib+i)-a(ie+i)
            a2=a(ic+i)+a(id+i)
            a4=a(ic+i)-a(id+i)
            a5=a(ia+i)+0.25*(a3-a4)
            a6=qrt5*(a3+a4)
            c(ja+j)=a5+a6
            c(jb+j)=a5-a6
            c(jc+j)=a(ia+i)-(a3-a4)
            d(ja+j)=-sin36*a1-sin72*a2
            d(jb+j)=-sin72*a1+sin36*a2
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
         RETURN
      END IF
   z=1.0/FLOAT(n)
   zqrt5=z*qrt5
   zsin36=z*sin36
   zsin72=z*sin72
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         a1=a(ib+i)+a(ie+i)
         a3=a(ib+i)-a(ie+i)
         a2=a(ic+i)+a(id+i)
         a4=a(ic+i)-a(id+i)
         a5=z*(a(ia+i)-0.25*(a1+a2))
         a6=zqrt5*(a1-a2)
         c(ja+j)=z*(a(ia+i)+(a1+a2))
         c(jb+j)=a5+a6
         c(jc+j)=a5-a6
         d(jb+j)=-zsin72*a3-zsin36*a4
         d(jc+j)=-zsin36*a3+zsin72*a4
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 6
!     -------------------

ELSE IF (igo == 5) THEN
   ia=1
   ib=ia+iink
   ic=ib+iink
   id=ic+iink
   ie=id+iink
   if=ie+iink
   ja=1
   jb=ja+(2*m-la)*inc2
   jc=jb+2*m*inc2
   jd=jc+2*m*inc2
   je=jc
   jf=jb
   IF (la /= m) THEN
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            a11=(a(ic+i)+a(IF+i))+(a(ib+i)+a(ie+i))
            c(ja+j)=(a(ia+i)+a(id+i))+a11
            c(jc+j)=(a(ia+i)+a(id+i)-0.5*a11)
            d(jc+j)=sin60*((a(ic+i)+a(IF+i))-(a(ib+i)+a(ie+i)))
            a11=(a(ic+i)-a(IF+i))+(a(ie+i)-a(ib+i))
            c(jb+j)=(a(ia+i)-a(id+i))-0.5*a11
            d(jb+j)=sin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(IF+i)))
            c(jd+j)=(a(ia+i)-a(id+i))+a11
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
      ja=ja+jink
      jink=2*jink
      jb=jb+jink
      jc=jc+jink
      jd=jd-jink
      je=je-jink
      jf=jf-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      IF (jc /= jd) THEN
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
            jbase=0
            DO l=1,la
               i=ibase
               j=jbase

               DO ijk=1,lot
                  a1=c1*a(ib+i)+s1*b(ib+i)
                  b1=c1*b(ib+i)-s1*a(ib+i)
                  a2=c2*a(ic+i)+s2*b(ic+i)
                  b2=c2*b(ic+i)-s2*a(ic+i)
                  a3=c3*a(id+i)+s3*b(id+i)
                  b3=c3*b(id+i)-s3*a(id+i)
                  a4=c4*a(ie+i)+s4*b(ie+i)
                  b4=c4*b(ie+i)-s4*a(ie+i)
                  a5=c5*a(IF+i)+s5*b(IF+i)
                  b5=c5*b(IF+i)-s5*a(IF+i)
                  a11=(a2+a5)+(a1+a4)
                  a20=(a(ia+i)+a3)-0.5*a11
                  a21=sin60*((a2+a5)-(a1+a4))
                  b11=(b2+b5)+(b1+b4)
                  b20=(b(ia+i)+b3)-0.5*b11
                  b21=sin60*((b2+b5)-(b1+b4))
                  c(ja+j)=(a(ia+i)+a3)+a11
                  d(ja+j)=(b(ia+i)+b3)+b11
                  c(jc+j)=a20-b21
                  d(jc+j)=a21+b20
                  c(je+j)=a20+b21
                  d(je+j)=a21-b20
                  a11=(a2-a5)+(a4-a1)
                  a20=(a(ia+i)-a3)-0.5*a11
                  a21=sin60*((a4-a1)-(a2-a5))
                  b11=(b5-b2)-(b4-b1)
                  b20=(b3-b(ia+i))-0.5*b11
                  b21=sin60*((b5-b2)+(b4-b1))
                  c(jb+j)=a20-b21
                  d(jb+j)=a21-b20
                  c(jd+j)=a11+(a(ia+i)-a3)
                  d(jd+j)=b11+(b3-b(ia+i))
                  c(jf+j)=a20+b21
                  d(jf+j)=a21+b20
                  i=i+inc3
                  j=j+inc4
                  END DO
               ibase=ibase+inc1
               jbase=jbase+inc2
               END DO
            ibase=ibase+ijump
            ja=ja+jink
            jb=jb+jink
            jc=jc+jink
            jd=jd-jink
            je=je-jink
            jf=jf-jink
            END DO
         IF (jc.gt.jd) THEN 
           IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
           RETURN
           END IF
         END IF
      jbase=0
      DO l=1,la
         i=ibase
         j=jbase

         DO ijk=1,lot
            c(ja+j)=(a(ia+i)+0.5*(a(ic+i)-a(ie+i)))                     &
                          +sin60*(a(ib+i)-a(IF+i))
            d(ja+j)=-(a(id+i)+0.5*(a(ib+i)+a(IF+i)))                    &
                          -sin60*(a(ic+i)+a(ie+i))
            c(jb+j)=a(ia+i)-(a(ic+i)-a(ie+i))
            d(jb+j)=a(id+i)-(a(ib+i)+a(IF+i))
            c(jc+j)=(a(ia+i)+0.5*(a(ic+i)-a(ie+i)))                     &
                          -sin60*(a(ib+i)-a(IF+i))
            d(jc+j)=-(a(id+i)+0.5*(a(ib+i)+a(IF+i)))                    &
                          +sin60*(a(ic+i)+a(ie+i))
            i=i+inc3
            j=j+inc4
            END DO
         ibase=ibase+inc1
         jbase=jbase+inc2
         END DO
         IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
         RETURN
      END IF
   z=1.0/FLOAT(n)
   zsin60=z*sin60
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         a11=(a(ic+i)+a(IF+i))+(a(ib+i)+a(ie+i))
         c(ja+j)=z*((a(ia+i)+a(id+i))+a11)
         c(jc+j)=z*((a(ia+i)+a(id+i))-0.5*a11)
         d(jc+j)=zsin60*((a(ic+i)+a(IF+i))-(a(ib+i)+a(ie+i)))
         a11=(a(ic+i)-a(IF+i))+(a(ie+i)-a(ib+i))
         c(jb+j)=z*((a(ia+i)-a(id+i))-0.5*a11)
         d(jb+j)=zsin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(IF+i)))
         c(jd+j)=z*((a(ia+i)-a(id+i))+a11)
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
      IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
      RETURN

!     coding for factor 8
!     -------------------

ELSE
   ibad=3
   ierr=3
   IF (la /= m) THEN 
     IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
     RETURN
   END IF
   ia=1
   ib=ia+iink
   ic=ib+iink
   id=ic+iink
   ie=id+iink
   if=ie+iink
   ig=IF+iink
   ih=ig+iink
   ja=1
   jb=ja+la*inc2
   jc=jb+2*m*inc2
   jd=jc+2*m*inc2
   je=jd+2*m*inc2
   z=1.0/FLOAT(n)
   zsin45=z*SQRT(0.5)
   DO l=1,la
      i=ibase
      j=jbase

      DO ijk=1,lot
         c(ja+j)=z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))+              &
                   ((a(id+i)+a(ih+i))+(a(ib+i)+a(IF+i))))
         c(je+j)=z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))-              &
                   ((a(id+i)+a(ih+i))+(a(ib+i)+a(IF+i))))
         c(jc+j)=z*((a(ia+i)+a(ie+i))-(a(ic+i)+a(ig+i)))
         d(jc+j)=z*((a(id+i)+a(ih+i))-(a(ib+i)+a(IF+i)))
         c(jb+j)=z*(a(ia+i)-a(ie+i))                                    &
                  +zsin45*((a(ih+i)-a(id+i))-(a(IF+i)-a(ib+i)))
         c(jd+j)=z*(a(ia+i)-a(ie+i))                                    &
                  -zsin45*((a(ih+i)-a(id+i))-(a(IF+i)-a(ib+i)))
         d(jb+j)=zsin45*((a(ih+i)-a(id+i))+(a(IF+i)-a(ib+i)))           &
                      +z*(a(ig+i)-a(ic+i))
         d(jd+j)=zsin45*((a(ih+i)-a(id+i))+(a(IF+i)-a(ib+i)))           &
                  -z*(a(ig+i)-a(ic+i))
         i=i+inc3
         j=j+inc4
         END DO
      ibase=ibase+inc1
      jbase=jbase+inc2
      END DO
   END IF

!     final return
!     ------------

ibad=0
ierr=ibad
IF (lhook) CALL dr_hook('QPASSM',zhook_out,zhook_handle)
RETURN

END SUBROUTINE qpassm

END MODULE qpassm_mod
