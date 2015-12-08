! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer. 
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution. 
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Fast-j routine for calculating online photolysis rates
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!    This code is written to UMDP3 vn8.2 programming standards
!
! ######################################################################
!
      MODULE  FASTJ_MIE
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
      PUBLIC
      SAVE

      INTEGER, PARAMETER                      :: nl=2700
!                        Max no of levs after insertion of extra Mie levs

      INTEGER, PARAMETER                      :: nlevs_mie=2*nl
!                        No of levs in Mie grid: 2*(2*lpar+2+jaddto(1))+3

      INTEGER, PARAMETER                      :: m_gauss=4
!                        Number of Gauss points used (fixed at 4)

      INTEGER, PARAMETER                      :: m=1
!                        Solve eqn of R.T. only for first-order M=1

      INTEGER, PARAMETER                      :: n=m_gauss
!                        Number of quadrature points (fixed at 4)

      INTEGER, PARAMETER                      :: mfit=2*m_gauss
!                        Expansion of phase function (fixed at 8)

      INTEGER, PARAMETER                      :: nwb1=1
!                        First wlength bin to use - hardwired for vector'n      

      INTEGER, PARAMETER                      :: nwb2=7
!                        Last wlength bin to use - also use for dimensioning

      INTEGER                                 :: nwb3
      INTEGER                                 :: nd

      INTEGER, ALLOCATABLE                    :: jadsub(:)
      
      REAL                                    :: radius

      REAL,DIMENSION(m_gauss)                 :: wt,emu
      REAL,DIMENSION(m_gauss,2*m_gauss)       :: pm
      REAL, ALLOCATABLE                       :: zrefl(:)
      REAL, ALLOCATABLE                       :: zflux(:)
      REAL, ALLOCATABLE                       :: zu0(:)
      REAL, ALLOCATABLE                       :: v1(:,:)
      REAL, ALLOCATABLE                       :: pm0(:,:)
      REAL, ALLOCATABLE                       :: pomega(:,:,:)
      REAL, ALLOCATABLE                       :: ztau(:,:)
      REAL, ALLOCATABLE                       :: fz(:,:)
      REAL, ALLOCATABLE                       :: fj(:,:)
      REAL, ALLOCATABLE                       :: a(:,:,:)
      REAL, ALLOCATABLE                       :: c1(:,:,:)
      REAL, ALLOCATABLE                       :: h(:,:,:)
      REAL, ALLOCATABLE                       :: aa(:,:,:)
      REAL, ALLOCATABLE                       :: cc(:,:,:)
      REAL, ALLOCATABLE                       :: s(:,:,:)
      REAL, ALLOCATABLE                       :: rr(:,:,:)
      REAL, ALLOCATABLE                       :: b(:,:,:,:)
      REAL, ALLOCATABLE                       :: dd(:,:,:,:)
      REAL, ALLOCATABLE                       :: w(:,:,:,:)
      REAL, ALLOCATABLE                       :: u1(:,:,:,:)

      INTERFACE FASTJ_MIE_ALLOC_MEM
        MODULE PROCEDURE FASTJ_MIE_ALLOC_MEM
      END INTERFACE FASTJ_MIE_ALLOC_MEM

      INTERFACE FASTJ_MIE_DEALLOC_MEM
        MODULE PROCEDURE FASTJ_MIE_DEALLOC_MEM
      END INTERFACE FASTJ_MIE_DEALLOC_MEM

      CONTAINS

! ######################################################################
       SUBROUTINE FASTJ_MIE_ALLOC_MEM (kpcx_lo)
         USE        FASTJ_DATA,         ONLY: nc
         IMPLICIT   NONE

         INTEGER,INTENT(IN)                    :: kpcx_lo

         INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
         INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
         REAL(KIND=jprb)               :: zhook_handle

         IF (lhook) CALL dr_hook('FASTJ_MIE:FASTJ_MIE_ALLOC_MEM',zhook_in,zhook_handle)
         nwb3 = nwb2*kpcx_lo

         ALLOCATE(a(nwb3,m_gauss,nlevs_mie))
         ALLOCATE(b(m_gauss+1,m_gauss+1,nwb3,nlevs_mie))
         ALLOCATE(c1(nwb3,m_gauss,nlevs_mie))
         ALLOCATE(h(nwb3,m_gauss,nlevs_mie))
         ALLOCATE(aa(nwb3,m_gauss,m_gauss))
         ALLOCATE(cc(nwb3,m_gauss,m_gauss))
         ALLOCATE(s(nwb3,m_gauss,m_gauss))
         ALLOCATE(w(nwb3,m_gauss,m_gauss,2))
         ALLOCATE(u1(nwb3,m_gauss,m_gauss,2))
         ALLOCATE(v1(nwb3,m_gauss))
         ALLOCATE(pm0(kpcx_lo,2*m_gauss))
         ALLOCATE(pomega(nwb3,2*m_gauss,nlevs_mie))
         ALLOCATE(ztau(nwb3,nlevs_mie))
         ALLOCATE(fz(nwb3,nlevs_mie))
         ALLOCATE(fj(nwb3,nlevs_mie))
         ALLOCATE(dd(nwb3,m_gauss,m_gauss,nlevs_mie))
         ALLOCATE(rr(nwb3,m_gauss,nlevs_mie))
         ALLOCATE(zrefl(kpcx_lo))
         ALLOCATE(zflux(nwb3))
         ALLOCATE(zu0(kpcx_lo))
         ALLOCATE(jadsub(nc))

         IF (lhook) CALL dr_hook('FASTJ_MIE:FASTJ_MIE_ALLOC_MEM',zhook_out,zhook_handle)

         RETURN
       END SUBROUTINE FASTJ_MIE_ALLOC_MEM
       
! ######################################################################
       SUBROUTINE FASTJ_MIE_DEALLOC_MEM
  
         IMPLICIT   NONE

         INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
         INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
         REAL(KIND=jprb)               :: zhook_handle

         IF (lhook) CALL dr_hook('FASTJ_MIE:FASTJ_MIE_DEALLOC_MEM',zhook_in,zhook_handle)

         DEALLOCATE(a)
         DEALLOCATE(b)
         DEALLOCATE(c1)
         DEALLOCATE(h)
         DEALLOCATE(aa)
         DEALLOCATE(cc)
         DEALLOCATE(s)
         DEALLOCATE(w)
         DEALLOCATE(u1)
         DEALLOCATE(v1)
         DEALLOCATE(pm0)
         DEALLOCATE(pomega)
         DEALLOCATE(ztau)
         DEALLOCATE(fz)
         DEALLOCATE(fj)
         DEALLOCATE(dd)
         DEALLOCATE(rr)
         DEALLOCATE(zrefl)
         DEALLOCATE(zflux)
         DEALLOCATE(zu0)
         DEALLOCATE(jadsub)

         IF (lhook) CALL dr_hook('FASTJ_MIE:FASTJ_MIE_DEALLOC_MEM',zhook_out,zhook_handle)

         RETURN
       END SUBROUTINE FASTJ_MIE_DEALLOC_MEM

      END MODULE   FASTJ_MIE
