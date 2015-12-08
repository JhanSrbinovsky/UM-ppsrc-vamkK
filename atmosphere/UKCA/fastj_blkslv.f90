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
!
! ######################################################################
!
      SUBROUTINE FASTJ_BLKSLV

      USE     FASTJ_DATA
      USE     FASTJ_MIE
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      INTEGER :: i, j, k, id, KW

      REAL    ::  xsum

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
!     This routine solves the block tri-diagonal system:
!     A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('FASTJ_BLKSLV',zhook_in,zhook_handle)
! DEPENDS ON: fastj_gen
      CALL FASTJ_GEN

!-----------UPPER BOUNDARY ID=1

      DO KW=NWB1,NWB3
! DEPENDS ON: fastj_matin4
        CALL FASTJ_MATIN4 (B(1,1,KW,1))
      ENDDO
      
      DO I=1,N
        DO KW=NWB1,NWB3
          RR(KW,I,1) = 0.0d0
        ENDDO
        DO J=1,N
          DO KW=NWB1,NWB3
            xSUM = 0.0d0
            DO K=1,N
              xSUM = xSUM - B(I,K,KW,1)*CC(KW,K,J)
            ENDDO
            DD(KW,I,J,1) = xSUM
            RR(KW,I,1) = RR(KW,I,1) + B(I,J,KW,1)*H(KW,J,1)
          ENDDO
        ENDDO
      ENDDO

!----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1

      DO ID=2,ND-1
        IF (N == 4) THEN                     !kk Hardcode loop for n=4
          DO I=1,N
            DO KW=NWB1,NWB3
              B(I,1,KW,ID) = B(I,1,KW,ID) + A(KW,I,ID)*DD(KW,I,1,ID-1)
              B(I,2,KW,ID) = B(I,2,KW,ID) + A(KW,I,ID)*DD(KW,I,2,ID-1)
              B(I,3,KW,ID) = B(I,3,KW,ID) + A(KW,I,ID)*DD(KW,I,3,ID-1)
              B(I,4,KW,ID) = B(I,4,KW,ID) + A(KW,I,ID)*DD(KW,I,4,ID-1)
              H(KW,I,ID)   = H(KW,I,ID)   - A(KW,I,ID)*RR(KW,I,ID-1)
            ENDDO 
          ENDDO 
        ELSE 
          DO I=1,N
            DO KW=NWB1,NWB3
              DO J=1,N
              B(I,J,KW,ID) = B(I,J,KW,ID) + A(KW,I,ID)*DD(KW,I,J,ID-1)
              ENDDO
              H(KW,I,ID) = H(KW,I,ID) - A(KW,I,ID)*RR(KW,I,ID-1)
            ENDDO
          ENDDO
        ENDIF

        DO KW=NWB1,NWB3
!DEPENDS ON: fastj_matin4
          CALL FASTJ_MATIN4 (B(1,1,KW,ID))
        ENDDO

        IF (N == 4) THEN                     !kk Hardcode loop for n=4
          DO KW=NWB1,NWB3
!CDIR unroll=N
            DO I=1,N
!             j = 1
              RR(KW,I,ID) = B(I,1,KW,ID)*H(KW,1,ID)
              DD(KW,I,1,ID) = - B(I,1,KW,ID)*C1(KW,1,ID)
!             j = 2
              RR(KW,I,ID) = RR(KW,I,ID) + B(I,2,KW,ID)*H(KW,2,ID)
              DD(KW,I,2,ID) = - B(I,2,KW,ID)*C1(KW,2,ID)
!             j = 3
              RR(KW,I,ID) = RR(KW,I,ID) + B(I,3,KW,ID)*H(KW,3,ID)
              DD(KW,I,3,ID) = - B(I,3,KW,ID)*C1(KW,3,ID)
!             j = 4
              RR(KW,I,ID) = RR(KW,I,ID) + B(I,4,KW,ID)*H(KW,4,ID)
              DD(KW,I,4,ID) = - B(I,4,KW,ID)*C1(KW,4,ID)
            ENDDO
          ENDDO
        ELSE
          DO I=1,N
            DO KW=NWB1,NWB3
              RR(KW,I,ID) = 0.0d0
            ENDDO
            DO J=1,N
              DO KW=NWB1,NWB3
                RR(KW,I,ID) = RR(KW,I,ID) + B(I,J,KW,ID)*H(KW,J,ID)
                DD(KW,I,J,ID) = - B(I,J,KW,ID)*C1(KW,J,ID)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDDO                        ! End of ID loop

!---------FINAL DEPTH POINT: ND

       DO I=1,N
        DO J=1,N
          DO KW=NWB1,NWB3
            xSUM = 0.0d0
            DO K=1,N
              xSUM = xSUM + AA(KW,I,K)*DD(KW,K,J,ND-1)
            ENDDO
            B(I,J,KW,ND) = B(I,J,KW,ND) + xSUM
            H(KW,I,ND) = H(KW,I,ND) - AA(KW,I,J)*RR(KW,J,ND-1)
          ENDDO
        ENDDO
      ENDDO

      DO KW=NWB1,NWB3
!DEPENDS ON: fastj_matin4
        CALL FASTJ_MATIN4 (B(1,1,KW,ND))
      ENDDO
      DO I=1,N
        do KW=NWB1,NWB3
          RR(KW,I,ND) = 0.0d0
        enddo
        DO J=1,N
          do KW=NWB1,NWB3
            RR(KW,I,ND) = RR(KW,I,ND) + B(I,J,KW,ND)*H(KW,J,ND)
          enddo
        ENDDO
      ENDDO

!-----------BACK SOLUTION

      DO ID=ND-1,1,-1
        DO I=1,N
          DO KW=NWB1,NWB3
!CDIR unroll=N
            DO J=1,N
              RR(KW,I,ID) = RR(KW,I,ID) + DD(KW,I,J,ID)*RR(KW,J,ID+1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!----------MEAN J & H

      DO ID=1,ND
        DO KW=NWB1,NWB3
          FJ(KW,ID) = 0.0d0
        ENDDO
      ENDDO
      DO ID=1,ND,2
        DO I=1,N
          DO KW=NWB1,NWB3
            FJ(KW,ID) = FJ(KW,ID) + RR(KW,I,ID)*WT(I)
          ENDDO
        ENDDO
      ENDDO

      DO ID=2,ND,2
        DO I=1,N
          DO KW=NWB1,NWB3
            FJ(KW,ID) = FJ(KW,ID) + RR(KW,I,ID)*WT(I)*EMU(I)
          ENDDO
        ENDDO
      ENDDO

      IF (lhook) CALL dr_hook('FASTJ_BLKSLV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_BLKSLV
