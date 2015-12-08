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
!   Set aerosol/cloud types and define black carbon profile
!
!     MX       Number of different types of aerosol to be considered
!     MIEDX    Index of aerosol types in jv_spec.dat - hardwire in here
!     LAERO    Switch for aerosol climatology - use BC absorber if false
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
      SUBROUTINE FASTJ_SET_AER
      USE FASTJ_DATA
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE


      INTEGER            :: i          ! Loop variable

      CHARACTER (LEN=80) :: cmessage   ! Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     Initialise aerosol index

      IF (lhook) CALL dr_hook('FASTJ_SET_AER',zhook_in,zhook_handle)
      DO i = 1,mx
        miedx(i) = 0
      ENDDO

!     Select Aerosol/Cloud types to be used - define types here

      MIEDX(1) = 10    !  Water Cloud (Deirmenjian 8 micron)
      MIEDX(2) = 14    !  Irregular Ice Cloud (Mishchenko)
!     MIEDX(3) = 20    !  Dust - large particles
!     MIEDX(4) = 17    !  Dust - sub-micron
!     MIEDX(5) = 25    !  Sulphate
!     MIEDX(6) = 22    !  Black carbon
!     MIEDX(7) = 23    !  Organic Carbon
!     MIEDX(8) = 24    !  Sea-Salt

!     sulphur is used as 3. aerosol.
!     Additional aersol types are added in set_prof, please add here the
!     corresponding MIEDX values and change MX in mo_fastj_data.f  

      DO i = 3,mx
        miedx(i) =  25  !  Sulphate
      ENDDO

!     Convert from 550nm to 1000nm optical depths (dependent on aerosol type)
!     Necessary for aerosol climatology only

      zmiex(1) = 1.d0
      zmiex(2) = 1.d0
      DO i = 3,MX
        zmiex(i) = qaa(4,miedx(i))*4.d0/                            &
     &            (qaa(3,miedx(i))*3.d0+qaa(2,miedx(i)))
      ENDDO

!     Ensure all 'MX' types are valid selections

      DO i = 1,MX
        WRITE(6,*) miedx(i),titlea(miedx(i))
        IF (miedx(i) > naa .or. miedx(i) <= 0) THEN
          WRITE(6,*) MIEDX(i),NAA
          cmessage= 'MIEDX(i) is negative or less than naa'

          CALL EREPORT('FASTJ_SET_AER',naa,cmessage)
        ENDIF
      ENDDO

!     Approximate Black Carbon up to 10 km; surface 200 ng/m3  
!     (Liousse et al) Scale: 1 ng/m3 = 1.0d-15 g/cm3 
!     (1.0d-11 g/m2/cm as BREF is in cm))

!     Simple place-holder profile

      DO i = 1,51
        bref(i) = 10.d0*1.0d-11
        IF (i > 6) bref(i) = 0.d0
      ENDDO

!     Simple altitude scaling for aerosol climatology

      aref(1) = 0.5d0
      DO i = 2,5
        aref(i) = aref(i-1)*0.5d0
      ENDDO
      DO i = 6,51
        aref(i) = 0.d0
      ENDDO

      IF (lhook) CALL dr_hook('FASTJ_SET_AER',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_SET_AER
