! *****************************COPYRIGHT*******************************
! (c) [University of California] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
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
!       copyright noticed, this list of conditions and the following
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
!   Fast-jx is an updated routine for calculating online photolysis rates
!   This module contains routines that read in the phase factors etc from file
!   Based upon the fast_specs routine, though with large differences caused by
!   differences between fast-j and fast-jx
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
      MODULE FASTJX_SPECS

      USE    FASTJX_DATA
      USE PrintStatus_mod
      USE ereport_mod, ONLY : ereport
      USE yomhook,     ONLY: lhook, dr_hook
      USE parkind1,    ONLY: jprb, jpim
      IMPLICIT NONE

      !*********************************************
      ! STRINGS

      ! String containing description of data set
      CHARACTER(LEN=78) :: TITLE0

      ! String containing cloud/aerosol scattering
      CHARACTER(LEN=7)  ::  TITLAA(A_)

      ! String containing species being photolysed
      CHARACTER(LEN=7)  ::  TITLEJ(X_)

      ! Dummy strings containing duplicate species strings
      CHARACTER(LEN=7)  ::  TITLEJ2,TITLEJ3

      !**********************************************
      ! INTERFACES

      INTERFACE FASTJX_RD_JS
        MODULE PROCEDURE FASTJX_RD_JS
      END INTERFACE FASTJX_RD_JS

      INTERFACE FASTJX_RD_MIE
        MODULE PROCEDURE FASTJX_RD_MIE
      END INTERFACE FASTJX_RD_MIE

      INTERFACE FASTJX_RD_XXX
        MODULE PROCEDURE FASTJX_RD_XXX
      END INTERFACE FASTJX_RD_XXX

      PUBLIC FASTJX_RD_JS, FASTJX_RD_MIE, FASTJX_RD_XXX

      CONTAINS

! ######################################################################
      ! Reads in quantum yields and labels
      ! Copied from fastj_rd_js routine
      SUBROUTINE FASTJX_RD_JS

      USE UKCA_CHEM_DEFS_MOD,   ONLY: ratj_t, ratj_defs
      IMPLICIT NONE
!     Read in quantum yield 'jfacta' and fastj label 'jlabel'

!     jfacta    Quantum yield (or multiplication factor) for photolysis
!     jlabel    Reference label identifying appropriate J-value to use

      INTEGER :: i

      CHARACTER(LEN=10) :: adjusted_fname
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('FASTJX_SPECS.FASTJX_RD_JS',zhook_in,     &
                               zhook_handle)
      adjusted_fname = ' '

      DO i=1,jppj 
        jfacta(i)=ratj_defs(i)%jfacta/100.E0
        adjusted_fname=TRIM(ADJUSTL(ratj_defs(i)%fname))
        jlabel(i)=adjusted_fname(1:7)
        IF (PrintStatus >= PrStatus_Diag)                               &
                   WRITE(6,'(I6,E12.3,A12)') i,jfacta(i),jlabel(i)
      END DO

      IF (lhook) CALL dr_hook('FASTJX_SPECS.FASTJX_RD_JS',zhook_out,    &
                               zhook_handle)
      RETURN
      END SUBROUTINE FASTJX_RD_JS


! ######################################################################
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
!-----------------------------------------------------------------------  
      SUBROUTINE FASTJX_RD_MIE(nj1,namfil)
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nj1  ! Channel number for reading data file
      CHARACTER(LEN=*), INTENT(IN) :: namfil ! Name of scattering data file 
                                             ! (e.g., FJX_scat.dat)
      CHARACTER (LEN=70) :: cmessage        ! Contains string for error handling
      INTEGER :: errcode                    ! error code

      INTEGER  i, j, k                      ! Loop variables
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      !************************************
      ! End of Header
      IF (lhook) CALL dr_hook('FASTJX_SPECS.FASTJX_RD_MIE',zhook_in,    &
                               zhook_handle)

      ! Open data file containing aerosl/cloud data
      OPEN (nj1,FILE=NAMFIL,STATUS='OLD',FORM='FORMATTED')

      ! Read number of data types and title
      READ (nj1,'(i2,a78)') naa,title0

      ! If the number of data types exceeds maximum allowed exit with an error
      IF (naa > a_) then
        cmessage = 'Too many scattering data sets' 
        errcode = 100
        CALL EREPORT('FASTJX_RD_MIE',errcode,cmessage) 
      END IF

      ! Read Cloud layering variables
      READ (nj1,'(5x,I5,2F10.5)') jtaumx,atau,atau0
      IF (printstatus >= prstatus_oper) WRITE(6,'(a,2F9.5,I5)')         &
            ' atau/atau0/jmx',atau,atau0,jtaumx

      ! Read blank line
      READ (nj1,*)

      ! Loop over aerosol types
      DO J = 1,naa

        ! Read title, effective radius and density
        READ (NJ1,'(3x,a20,32x,f5.3,15x,f5.3)')  TITLAA(J),RAA(J),DAA(J)

        ! Loop over 5 wavelength bins
        DO k = 1,5     
          ! read wavelength, q, scattering albedo and phases
          READ (NJ1,'(f4.0,f7.4,f7.4,7f6.3,1x,f7.3,f8.4)') &
            waa(k,j),qaa(k,j),saa(k,j),(paa(i,k,j),i=2,8)
          ! set first phase to 0
          paa(1,k,j) = 1.E0

        END DO ! wavelengths
      END DO ! aerosols
 
      ! Close file
      CLOSE(NJ1)

      IF (printstatus >= prstatus_oper) THEN
      ! Output some information
        WRITE(6,'(A,9F8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):'   &
               ,(waa(k,1),k=1,5)

      ! Output file title
        WRITE(6,*) TITLE0

      ! Loop over aerosol types writing radius, density and q
        DO j=1,NAA
          WRITE(6,'(i3,1x,a8,7f8.3)') j,titlaa(j),raa(j),daa(j),        &
                                      (qaa(k,j),k=1,5)
        END DO
      END IF

      IF (lhook) CALL dr_hook('FASTJX_SPECS.FASTJX_RD_MIE',zhook_out,   &
                              zhook_handle)
      RETURN
      END SUBROUTINE FASTJX_RD_MIE

!########################################################################
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections. 
!
!>>>>NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 (parm_CTM.f) should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (parm_MIE.f).
!           if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
!
!-----------------------------------------------------------------------
      SUBROUTINE FASTJX_RD_XXX(nj1,namfil)

      IMPLICIT NONE

      INTEGER, INTENT(IN)       :: nj1 ! Channel number for reading data file
      CHARACTER(LEN=*), INTENT(IN)  :: namfil    ! Name of spectral data file 
                                                 ! (JX_spec.dat) 
      CHARACTER (LEN=70)        :: cmessage      ! String for error handling
      INTEGER                   :: errcode       ! errror code

      INTEGER ::  i, j, jj, k, iw                ! Loop variables
      INTEGER ::  nqqq                           ! No of cross sections read in
      INTEGER ::  nwww                           ! No of wavelength bins
      INTEGER ::  nqrd                           ! number of x-sections
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      !******************************
      ! EOH
      IF (lhook) CALL dr_hook('FASTJX_SPECS.FASTJX_RD_XXX',zhook_in,    &
                              zhook_handle)

! Initialise the temperatures to 0
      DO j = 1,x_
        DO k = 1,3
          tqq(k,j) = 0.E0
        ENDDO
      ENDDO

!----------spectral data----set for new format data J-ver8.3------------------
!         note that NJVAL = # J-values, but NQQQ (>NJVAL) = # Xsects read in 
!         for 2005a data, NJVAL = 62 (including a spare XXXX) and 
!              NQQQ = 64 so that 4 wavelength datasets read in for acetone
!         note NQQQ is not used outside this subroutine!
! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects 

      ! Open file containing cross sections
      OPEN (nj1,FILE=namfil,status='old',form='formatted')

      ! Read title of file
      READ (nj1,'(a)') title0

      ! Read number of photolysed species, number of x-sections & number of wavelength bins
      READ (NJ1,'(10x,5i5)') njval, nqrd, nwww

      ! set maximum and minimum wavelngth bins
      nw1 = 1
      nw2 = nwww

      ! Check that number of photolysed species and number of cross sections
      ! doesn't exceed maximum allowed. If either do then exit
      IF (njval > x_ .OR. nqrd > x_) THEN
        cmessage = 'Number of Cross Sections exceeds Maximum Allowed' 
        errcode = 100
        CALL EREPORT('FASTJX_RD_XXX',errcode,cmessage) 
      ENDIF
      
      !----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
      READ (NJ1,'(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            (wl(iw),iw=1,nwww)
      READ (NJ1,'(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            (fl(iw),iw=1,nwww)
      READ (NJ1,'(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            (qrayl(iw),iw=1,nwww)

      !---READ O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
      READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej(1),tqq(1,1), (qo2(iw,1),iw=1,nwww)
      READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej2,  tqq(2,1), (qo2(iw,2),iw=1,nwww)
      READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej3,  tqq(3,1), (qo2(iw,3),iw=1,nwww)

      READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej(2),tqq(1,2), (qo3(iw,1),iw=1,nwww)
      READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej2,  tqq(2,2), (qo3(iw,2),iw=1,nwww)
      READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej3,  tqq(3,2), (qo3(iw,3),iw=1,nwww)

      READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej(3),tqq(1,3), (q1d(iw,1),iw=1,nwww)
      READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej2,  tqq(2,3), (q1d(iw,2),iw=1,nwww)
      READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej3,  tqq(3,3), (q1d(iw,3),iw=1,nwww)

      ! WRITE information to unit 6
      IF (printstatus >= prstatus_oper) THEN
        DO j = 1,3
          WRITE(6,'(i6,a7,2e10.3)') j,titlej(j),(tqq(i,j),i=1,3)
        END DO
      END IF

      !---READ remaining species:  X-sections at 2 T_s
      JJ = 4
      DO J = 4,NQRD
 
        READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej(jj),tqq(1,jj),(qqq(iw,1,jj),iw=1,nwww)
        READ (NJ1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
            titlej2,  tqq(2,jj),(qqq(iw,2,jj),iw=1,nwww)

        !---include stratospheric J's (this also includes Cl and Br compounds!)
        IF (w_ == 18 .OR. titlej2(7:7) /= 'x') THEN
          IF (printstatus >= prstatus_oper)                               &
            WRITE(6,'(i6,a7,2e10.3)') jj,titlej(jj), (tqq(i,jj),i=1,2)
          JJ = JJ+1
        END IF

      END DO
      nqqq = jj-1
      njval = njval + (nqqq - nqrd)

      !---truncate number of wavelengths to DO troposphere-only
      IF (w_ /= wx_) THEN

        !---TROP-ONLY
        IF (w_ == 12) THEN
          IF (printstatus >= prstatus_oper) WRITE(6,'(a)')              &
            ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'
          nw2 = 12

          ! Remove first four wavelength bins from  total
          DO iw = 1,4
            wl(iw) = wl(iw+4)
            fl(iw) = fl(iw+4)
            qrayl(iw) = qrayl(iw+4)

            DO k = 1,3
              qo2(iw,k) = qo2(iw+4,k)
              qo3(iw,k) = qo3(iw+4,k)
              q1d(iw,k) = q1d(iw+4,k)
            END DO

            DO J = 4,nqqq
              qqq(iw,1,j) = qqq(iw+4,1,j)
              qqq(iw,2,j) = qqq(iw+4,2,j)
            END DO
          END DO

          ! Remove 9/10 wavelength bins from total
          DO iw = 5,12
            wl(iw) = wl(iw+6)
            fl(iw) = fl(iw+6)
            qrayl(iw) = qrayl(iw+6)

            DO k = 1,3
              qo2(iw,k) = qo2(iw+6,k)
              qo3(iw,k) = qo3(iw+6,k)
              q1d(iw,k) = q1d(iw+6,k)
            END DO
            DO j = 4,NQQQ
              qqq(iw,1,j) = qqq(iw+6,1,j)
              qqq(iw,2,j) = qqq(iw+6,2,j)
            END DO
          END DO

          !---TROP-QUICK  (must scale solar flux for W=5)
        ELSEIF (w_ == 8) THEN
          IF (printstatus >= prstatus_oper) WRITE(6,'(a)')              &
            ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'
          nw2 = 8

          DO iw = 1,1
            wl(iw) = wl(iw+4)
            fl(iw) = fl(iw+4)*2.E0
            qrayl(iw) = qrayl(iw+4)

            DO k = 1,3
              qo2(iw,k) = qo2(iw+4,k)
              qo3(iw,k) = qo3(iw+4,k)
              q1d(iw,k) = q1d(iw+4,k)
            END DO

            DO j = 4,nqqq
              qqq(iw,1,j) = qqq(iw+4,1,j)
              qqq(iw,2,j) = qqq(iw+4,2,j)
            END DO
          END DO

          DO iw = 2,8
            wl(iw) = wl(iw+10)
            fl(iw) = fl(iw+10)
            qrayl(iw) = qrayl(iw+10)

            DO k = 1,3
              qo2(iw,k) = qo2(iw+10,k)
              qo3(iw,k) = qo3(iw+10,k)
              q1d(iw,k) = q1d(iw+10,k)
            END DO

            DO j = 4,nqqq
              qqq(iw,1,j) = qqq(iw+10,1,j)
              qqq(iw,2,j) = qqq(iw+10,2,j)
            END DO
          END DO

        ELSE
           cmessage = 'Incorrect Number of Wavelength Bins, must be 8, 12 or 18' 
           errcode = 100
           CALL EREPORT('FASTJX_RD_XXX',errcode,cmessage) 
        END IF
      END IF

      ! Close file
      CLOSE(NJ1)
 
      !**************************************************************
      ! Map local indices to UKCA ones

      DO j = 1,NJVAL
        DO k = 1,jppj
          IF ( k == 1 .AND. printstatus >= prstatus_oper)               &
            WRITE(6,*) 'FASTJX Compare titles ', titlej(j), jlabel(k)
          IF (jlabel(k) == titlej(j)) jind(k) = j
        END DO 
      END DO

      DO k = 1,jppj

        IF (printstatus >= prstatus_oper) THEN
          WRITE(6,*) 'Comparing J-rate for photolysis reaction ',jlabel(k),' ?'
          WRITE(6,*) 'Using index ', jind(k), maxval(qqq(:,:,jind(k)))
        END IF

        IF (jfacta(k) < 1.E-20 .AND. printstatus >= prstatus_oper)       &
           WRITE(6,*) 'Not using photolysis reaction ',jlabel(k)

        IF (jind(k) == 0) THEN

          IF (jfacta(k) < 1.E-20) THEN
            jind(k)=1
          ELSE
            cmessage =                                                  &
             ' Which J-rate for photolysis reaction '//jlabel(k)//' ?'
            errcode = 1
            CALL EREPORT("FASTJX_RD_XXX",errcode,cmessage)
          END IF
        END IF
      END DO

      IF (lhook) CALL dr_hook('FASTJX_SPECS.FASTJX_RD_XXX',zhook_out,   &
                              zhook_handle)
      RETURN
      END SUBROUTINE FASTJX_RD_XXX

!#######################################################################

      END MODULE FASTJX_SPECS
