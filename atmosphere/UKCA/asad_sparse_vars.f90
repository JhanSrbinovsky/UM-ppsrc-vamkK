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
!    Module holding variables and routines for sparse algebra.
!    Part of the ASAD chemical solver. Contains the following
!    routines:
!      setup_spfuljac
!      spfuljac
!      splinslv2
!      spresolv2
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
!     Method
!     ------
!     Sparse algebra works in the same way as dense algebra, namely by LU
!     decomposition (Gaussian elimination) of the linear system. In contrast
!     to dense algebra, here we keep track of non-zero matrix elements
!     and hence cut out algebraic amnipulations whose result would be zero.
!     To make this efficient, species need to be reorder, such that those
!     with the fewest non-zero matrix elements (reactions) associated with
!     them, occur first in the list, and those with most (generally OH)
!     occur last.
!
! ######################################################################
!
      MODULE asad_sparse_vars

      USE ASAD_MOD
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE PrintStatus_mod
      USE Control_Max_Sizes
      IMPLICIT NONE

! TYPSIZE start
!   Description:
!     This file contains sizes needed for dynamic allocation of
!   main data arrays within the model. Sizes read in from the user
!   interface via NAMELISTs are passed by /COMMON/. Other control
!   sizes that are fundamental in the definition of data structures
!   are assigned by PARAMETER statements.
!
!   Declarations for the NLSIZES namelist are also held in the module
!   nlsizes_namelist_mod. That module is currently only used by the
!   reconfiguration, while the UM uses this include file.
!
! All sizes
! Not dependent on sub-model
! DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
! ATMOS START
! Main sizes of fields for each submodel
! Grid-related sizes for ATMOSPHERE submodel.
INTEGER :: ROW_LENGTH           ! No of points per local row
INTEGER :: global_ROW_LENGTH    ! Points per global row
INTEGER :: ROWS                 ! No of local (theta) rows
INTEGER :: global_ROWS          ! No of global (theta) rows
INTEGER :: MODEL_LEVELS         ! No of model levels
INTEGER :: LAND_FIELD           ! No of land points in field
INTEGER :: NTILES               ! No of land surface tiles
INTEGER :: NICE                 ! No. of sea ice thickness categories
INTEGER :: NICE_USE             ! No. of sea ice categories used fully
                                !  in surface exchange and radiation
                                !  (If nice>1 & nice_use=1, categories only 
                                !  partially used in surface exchange)

! Physics-related sizes for ATMOSPHERE submodel
INTEGER :: WET_LEVELS          ! No of moist-levels
INTEGER :: CLOUD_LEVELS        ! No of cloud-levels
INTEGER :: ST_LEVELS           ! No of soil temperature levels
INTEGER :: SM_LEVELS           ! No of soil moisture levels
INTEGER :: BL_LEVELS           ! No of boundary-layer-levels
INTEGER :: OZONE_LEVELS        ! No of ozone-levels
INTEGER :: TPPS_OZONE_LEVELS   ! No of tropopause-ozone-levels
INTEGER :: RIVER_ROWS          ! No of rows for river routing
INTEGER :: RIVER_ROW_LENGTH    ! Row length for river routing
! Dynamics-related sizes for ATMOSPHERE submodel

INTEGER :: TR_LEVELS            ! No of tracer-levels
INTEGER :: TR_VARS              ! No of passive tracers
INTEGER :: TR_LBC_VARS          ! No of tracers in lbcs 
INTEGER :: TR_UKCA              ! No of UKCA tracers
INTEGER :: TR_LBC_UKCA          ! No of UKCA tracer lbcs 

! For Small executables

! Grid related sizes for data structure
! Data structure sizes for ATMOSPHERE submodel
INTEGER :: A_PROG_LOOKUP     ! No of prognostic fields
INTEGER :: A_PROG_LEN        ! Total length of prog fields
INTEGER :: A_LEN_INTHD       ! Length of INTEGER header
INTEGER :: A_LEN_REALHD      ! Length of REAL header
INTEGER :: A_LEN2_LEVDEPC    ! No of LEVEL-dependent arrays
INTEGER :: A_LEN2_ROWDEPC    ! No of ROW-dependent arrays
INTEGER :: A_LEN2_COLDEPC    ! No of COLUMN-dependent arrays
INTEGER :: A_LEN2_FLDDEPC    ! No of FIELD arrays
INTEGER :: A_LEN_EXTCNST     ! No of EXTRA scalar constants
INTEGER :: A_LEN_CFI1        ! Length of compressed fld index 1
INTEGER :: A_LEN_CFI2        ! Length of compressed fld index 2
INTEGER :: A_LEN_CFI3        ! Length of compressed fld index 3
! atmos end

! Data structure sizes for ATMOSPHERE ANCILLARY file control
! routines
INTEGER :: NANCIL_LOOKUPSA  ! Max no of fields to be read

! Data structure sizes for ATMOSPHERE INTERFACE file control
! routines
INTEGER :: N_INTF_A          ! No of atmosphere interface areas
INTEGER :: MAX_INTF_MODEL_LEVELS ! Max no of model levs in all areas
INTEGER :: MAX_LBCROW_LENGTH ! Max no of lbc row length in all areas
INTEGER :: MAX_LBCROWS ! Max no of lbc rows in all areas

!  Data structure sizes for ATMOSPHERE BOUNDARY file control
! routines

! Sizes applicable to all configurations (DUMPS/FIELDSFILE)

INTEGER :: PP_LEN_INTHD   ! Length of PP file integer header
INTEGER :: PP_LEN_REALHD  ! Length of PP file real    header


      ! Grid related sizes for COUPLING between ATMOS and OCEAN
      ! submodels [For MPP, sizes are 'global' values over all
      ! PEs.Also needed for river routing]
      INTEGER:: AOCPL_IMT                ! Ocean rowlength
      INTEGER:: AOCPL_JMT                ! Ocean no. of rows
      INTEGER:: AOCPL_ROW_LENGTH         ! Atmos rowlength
      INTEGER:: AOCPL_P_ROWS             ! Atmos no. of p rows

      COMMON/SIZE_AOCPL/                                                &
        AOCPL_IMT, AOCPL_JMT, AOCPL_ROW_LENGTH, AOCPL_P_ROWS

! Other sizes passed from namelist into common blocks
! Any additions to this common block must be mirrored in nlsizes_namelist_mod.
COMMON/NLSIZES/                                                     &
    ROW_LENGTH,global_ROW_LENGTH,ROWS,global_ROWS,                  &
    LAND_FIELD,MODEL_LEVELS,WET_LEVELS,                             &
    NTILES, NICE, NICE_USE,                                         &
    CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
    OZONE_LEVELS,TPPS_OZONE_LEVELS,TR_VARS,TR_LBC_VARS,             &
    TR_UKCA,TR_LBC_UKCA,RIVER_ROWS,RIVER_ROW_LENGTH,                &
    A_PROG_LOOKUP,A_PROG_LEN,                                       &
    A_LEN_INTHD,A_LEN_REALHD,                                       &
    A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
    A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
    A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &    
    NANCIL_LOOKUPSA,                                                &    
    N_INTF_A, MAX_INTF_MODEL_LEVELS, MAX_LBCROW_LENGTH,             &
    MAX_LBCROWS, PP_LEN_INTHD,PP_LEN_REALHD

!-----------------------------------------------------------------
! data in STASHC#x member of the job library

! Data structure sizes for ATMOSPHERE submodel (config dependent)
INTEGER :: A_LEN2_LOOKUP   ! Total no of fields (incl diags)
INTEGER :: A_LEN_DATA      ! Total no of words of data
INTEGER :: A_LEN_D1        ! Total no of words in atmos D1

! Size of main data array for this configuration

INTEGER :: LEN_TOT             ! Length of D1 array
INTEGER :: N_OBJ_D1_MAX         ! No of objects in D1 array

COMMON/STSIZES/                                                     &
    A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                              &
    LEN_TOT,N_OBJ_D1_MAX
! global (ie. dump version) of *_LEN_DATA
INTEGER :: global_A_LEN_DATA

COMMON /MPP_STSIZES_extra/ global_A_LEN_DATA
! Sizes of Stash Auxillary Arrays and associated index arrays
! Initialised in UMINDEX and UMINDEX_A/O/W
INTEGER :: LEN_A_IXSTS
INTEGER :: LEN_A_SPSTS

COMMON /DSIZE_STS/                                                  &
    LEN_A_IXSTS, LEN_A_SPSTS
!     The number of land points is computed for each PE
!     before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field
      ! ----------------------------------------------------------------
      ! extra variables not passed through user interface

      ! fundamental data sizes :
      ! Fundamental parameter  sizes of data structure
      ! Sizes applicable to all configurations (HISTORY FILE)

      ! Length of history file in dump
      INTEGER, PARAMETER :: LEN_DUMPHIST = 0

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      ! Length of dump fixed header
      INTEGER, PARAMETER :: LEN_FIXHD = 256

      ! Size of a single LOOKUP header
      INTEGER, PARAMETER :: LEN1_LOOKUP  = 64
      INTEGER, PARAMETER :: MPP_LEN1_LOOKUP= 2

      ! Size of compressed LBC LOOKUP (only used internally and
      ! contains just the items which change between each set of LBCs
      INTEGER, PARAMETER :: LEN1_LBC_COMP_LOOKUP = 8

      ! Sizes applicable to all configurations (STASH)
      ! Moved to typstsz.h

      INTEGER:: INTF_LEN2_LEVDEPC !1st dim of interface out lev dep cons
      INTEGER:: INTF_LEN2_ROWDEPC !2nd dim of interface out Row dep cons
      INTEGER:: INTF_LEN2_COLDEPC !2nd dim of interface out Col dep cons
      
      COMMON /DSIZE/                                                    &
        INTF_LEN2_LEVDEPC,INTF_LEN2_ROWDEPC,INTF_LEN2_COLDEPC
      ! sub-model atmosphere   :
      ! Data structure sizes derived from grid size
      INTEGER:: A_LEN1_LEVDEPC ! IN: 1st dim of level  dep const
      INTEGER:: A_LEN1_ROWDEPC ! IN: 1st dim of row    dep const
      INTEGER:: A_LEN1_COLDEPC ! IN: 1st dim of column dep const
      INTEGER:: A_LEN1_FLDDEPC ! IN: 1st dim of field  dep const

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER :: INTF_LOOKUPSA        ! No of interface lookups.
      COMMON /DSIZE_A/                                                  &
        A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,A_LEN1_COLDEPC,    &
        INTF_LOOKUPSA

      ! sub-model atmosphere   : derived sizes
      ! derived from model grid/levels. Arakawa B-grid

                                  ! Size of fields on THETA grid:
      INTEGER :: THETA_FIELD_SIZE     ! IN: with no halos
      INTEGER :: THETA_OFF_SIZE       ! IN: with simple halos
      INTEGER :: THETA_HALO_SIZE      ! IN: with extended halos

                                  ! Size of fields on U grid:
      INTEGER :: U_FIELD_SIZE         ! IN: with no halos
      INTEGER :: U_OFF_SIZE           ! IN: with simple halos
      INTEGER :: U_HALO_SIZE          ! IN: with extended halos

                                  ! Size of fields on V grid
      INTEGER :: V_FIELD_SIZE         ! IN: with no halos
      INTEGER :: V_OFF_SIZE           ! IN: with simple halos
      INTEGER :: V_HALO_SIZE          ! IN: with extended halos

      INTEGER :: N_ROWS               ! IN: No of V-rows
      INTEGER :: N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/                                                  &
        N_ROWS,N_CCA_LEV,THETA_FIELD_SIZE,U_FIELD_SIZE,V_FIELD_SIZE,    &
        THETA_OFF_SIZE,U_OFF_SIZE,V_OFF_SIZE,                           &
        THETA_HALO_SIZE,U_HALO_SIZE,V_HALO_SIZE
      ! boundary updating      : derived values
      ! Variables describing the Atmosphere Lateral Boundary Conditions
      ! Local (per processor) information


! TYPSIZE end

      INTEGER, PARAMETER :: maxterms = 160  ! maximum number of nonzero
                                           ! terms for individual species
      INTEGER, PARAMETER :: maxfterms = 100    ! maximum number of terms
                                         ! involving fractional  products
      INTEGER :: total
      INTEGER :: total1    ! total number of nonzero entries in Jacobian

      INTEGER :: nposterms(spfjsize_max), nnegterms(spfjsize_max),      &
                 nfracterms(spfjsize_max)

      INTEGER :: posterms(spfjsize_max, maxterms)
      INTEGER :: negterms(spfjsize_max, maxterms)
      INTEGER :: fracterms(spfjsize_max, maxfterms)
      INTEGER :: base_tracer(spfjsize_max)

      REAL :: fraction(spfjsize_max, maxfterms)


      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SETUP_SPFULJAC(n_points)
!
!  Routine to set up pointer arrays for sparse full jacobian.
!  Note that this routine mirrors the "dense" fuljac routine
!  but only calculates where the nonsero entries would go
!  in a full Jacobian. It then calculates the pointer arrays
!  for a sparse representation of the Jacobian.

        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points

! Local variables

      INTEGER :: errcode                ! Variable passed to ereport


      INTEGER :: irj
      INTEGER :: ifamd
      INTEGER :: itrd
      INTEGER :: i
      INTEGER :: j
      INTEGER :: jc
      INTEGER :: itrcr
      INTEGER :: j3
      INTEGER :: jn
      INTEGER :: js
      INTEGER :: jl
      INTEGER :: i1
      INTEGER :: i2
      INTEGER :: is
      INTEGER :: kr
      INTEGER :: ikr
      INTEGER :: krj
      INTEGER :: ij1
      INTEGER :: itemp1
      INTEGER :: ik(jpmsp)
      INTEGER :: ij(jpmsp)
      INTEGER :: activity(jpctr)

      INTEGER, ALLOCATABLE :: p(:,:)    ! permutation matrix
      INTEGER, ALLOCATABLE :: temp(:,:) ! temporary index matrix
      INTEGER, ALLOCATABLE :: map(:,:)

      LOGICAL :: lj(jpmsp)

      CHARACTER :: ityped*2
      CHARACTER(LEN=72) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! index map for nonzero entries in Jacobian
!
!
      IF (lhook) CALL dr_hook('ASAD_SPARSE_VARS:SETUP_SPFULJAC',        &
                              zhook_in,zhook_handle)
      IF (.NOT. ALLOCATED(map)) ALLOCATE(map(jpctr, jpctr))
      map = 0
      DO i=1,jpctr
        map(i,i) = 1
      END DO
!
!
!     -----------------------------------------------------------------
!           2.  Calc. full Jacobian matrix.
!               ----- ---- -------- -------
!
      DO jc = 1, ntrf
        itrcr = nltrf(jc)
        DO j3 = 1,nmzjac(itrcr)
          irj = nzjac1(j3,itrcr)
!
          DO jn=1,jpmsp
            ij(jn)=njcoth(irj,jn)
            lj(jn)= (ij(jn) /= 0)
          END DO
          DO jn=1,jpmsp
            IF(lj(jn)) map(ij(jn),itrcr) = 1
          END DO
!
          IF(npdfr(irj,1) /= 0) THEN
            i1 = npdfr(irj,1)
            i2 = npdfr(irj,2)
            DO jn = i1, i2
              is = ntabpd(jn,1)
              map(is,itrcr) = 1
            END DO
          END IF
!
        END DO
      END DO
!
!  Go through the steady state additions to the Jacobian; currently
!  assume that only O(1D) and O(3P) are modelled, and that both are
!  required for the O3 loss rate.
!
      DO jc = 1, nstst
        DO j3 = 1,nmsjac(jc)
          irj = nsjac1(j3,jc)
          DO jn=1,jpmsp
            ij(jn)=njcoth(irj,jn)
            lj(jn)= (ij(jn) /= 0)
          END DO
          DO jn=1,jpmsp
            IF(lj(jn)) THEN
              map(ij(jn),ntro3 ) = 1
!              map(ij(jn),ntroh ) = 1
!              map(ij(jn),ntrho2) = 1
!              map(ij(jn),ntrno ) = 1
            END IF
          END DO
        END DO
      END DO
!
!     -----------------------------------------------------------------
!          4.  Add deposition terms to Jacobian diagonal.
!              --- ---------- ----- -- -------- ---------
!
      IF ( (ndepw  /=  0) .OR. (ndepd  /=  0) ) THEN
        DO js = 1, jpspec
          ifamd = moffam(js)
          itrd = madvtr(js)
          ityped = ctype(js)
          IF ( itrd /= 0 ) map(itrd,itrd) = 1
        END DO
      END IF
!
! produce pointer variables

      total = SUM(map)
      IF (mype == 0 .AND. printstatus >= prstatus_oper) WRITE(6,*)      &
        'TOTAL NUMBER OF NONZERO ENTRIES IN JACOBIAN: ',total

      is=0
      pointer2(:,:)=0
      DO i=1,jpctr
        DO j=1,jpctr
! calculate forward and backward pointers for sparse representation
          IF (map(i,j) == 1) THEN
            is = is + 1
! backward pointer
            pointer2(i,j) = is
          END IF
        END DO
      END DO

! Reorder species to minimize fill-in
      DO i=1,jpctr
        activity(i) = SUM(map(:,i)) + SUM(map(i,:))
        ro(i) = i
      END DO
      DO i=1,jpctr-1
        DO j=i+1,jpctr
          IF (activity(i) > activity(j)) THEN
! exchange i and j tracers if i is more active than j.
            itemp1 = ro(i)
            ro(i) = ro(j)
            ro(j) = itemp1
            itemp1 = activity(i)
            activity(i) = activity(j)
            activity(j) = itemp1
          END IF
        END DO
      END DO
      DO i=1,jpctr
        IF (mype == 0 .AND. printstatus >= prstatus_oper) THEN
          write(6,*) 'IN ASAD SETUP_SPFULJAC: ',&
               advt(ro(i)),' ',activity(i)
        END IF
      END DO

! reorganize pointer variable to account for varying fill-in
      IF (.NOT. ALLOCATED(p)) ALLOCATE(p(jpctr,jpctr))
      IF (.NOT. ALLOCATED(temp)) ALLOCATE(temp(jpctr,jpctr))

      p = 0
      DO i=1,jpctr
        p(i,ro(i)) = 1
      END DO

! Calculate P*A
      DO i=1,jpctr
        DO j=1,jpctr
          temp(i,j) = SUM(p(i,:) * pointer2(:,j))
        END DO
      END DO

! form P*A*P' -> A
      pointer(:,:) = 0.0
      DO i=1,jpctr
        DO j=1,jpctr
          pointer(i,j) = SUM(temp(i,:) * p(j,:))
        END DO
      END DO

      IF (ALLOCATED(p)) DEALLOCATE(p)
      IF (ALLOCATED(temp)) DEALLOCATE(temp)
      IF (ALLOCATED(map)) DEALLOCATE(map)

! calculate production and loss terms
      nposterms = 0
      nnegterms = 0
      nfracterms = 0
!     -----------------------------------------------------------------
!           2.  Calc. full Jacobian matrix.
!               ----- ---- -------- -------
!
      DO jc = 1, ntrf
        itrcr = nltrf(jc)
!
        DO j3 = 1,nmzjac(itrcr)
          irj = nzjac1(j3,itrcr)
!
          DO jn=1,jpmsp
            ij(jn)=njcoth(irj,jn)
            lj(jn)= (ij(jn) /= 0)
          END DO
          DO jn=1,2
            IF (lj(jn)) THEN
              i = pointer2(ij(jn),itrcr)
              nnegterms(i) = nnegterms(i) + 1
              IF (nnegterms(i) > maxterms) THEN
                errcode=1
                cmessage=' Increase maxterms'
                CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
              END IF
              negterms(i,nnegterms(i)) = irj
            END IF
          END DO
          IF (nfrpx(irj) == 0) THEN
            DO jn=3,jpmsp
              IF (lj(jn)) THEN
                i = pointer2(ij(jn),itrcr)
                nposterms(i) = nposterms(i) + 1
                IF (nposterms(i) > maxterms) THEN
                  errcode=2
                  cmessage=' Increase maxterms'
                  CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
                END IF
                posterms(i,nposterms(i)) = irj
              END IF
            END DO
          ELSE
            DO jn=3,jpmsp
              IF (lj(jn)) THEN
                i = pointer2(ij(jn),itrcr)
                nfracterms(i) = nfracterms(i) + 1
                IF (nfracterms(i) > maxfterms) THEN
                  errcode = 3
                  cmessage=' Increase maxfterms'
                  CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
                END IF
                IF (nfrpx(irj)+jn-3 > jpfrpx) THEN
                  cmessage = ' frpx array index > jpfrpx'
                  errcode = 3
                  write(6,*) cmessage
                  write(6,*) 'irj: ',irj,' nfrpx(irj): ',nfrpx(irj)
                  write(6,*) 'i: ',i,' jn: ',jn
                  CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
                END IF
                fracterms(i,nfracterms(i)) = irj
                fraction(i,nfracterms(i)) = frpx(nfrpx(irj)+jn-3)
              END IF
            END DO
          END IF
!
          IF (npdfr(irj,1) /= 0) THEN
            i1 = npdfr(irj,1)
            i2 = npdfr(irj,2)
            DO jn = i1, i2
              is = ntabpd(jn,1)
              i = pointer2(is,itrcr)
              nfracterms(i) = nfracterms(i) + 1
              IF (nfracterms(i) > maxfterms) THEN
                errcode=3
                cmessage=' Increase maxfterms'
                CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
              END IF
              fracterms(i,nfracterms(i)) = irj
              fraction(i,nfracterms(i)) = ztabpd(jn,1)
            END DO
          END IF
!
        END DO
      END DO

! calculate base_tracer (i.e., number of base tracer that goes with each
! matrix element
      DO i=1,jpctr
        DO j=1,jpctr
          i1 = pointer2(j,i)
          IF (i1 > 0) base_tracer(i1)=i
        END DO
      END DO

! calculate number of matrix elements after factorization. Make sure it is
! smaller than the set maximum.
      total1 = total
      pointer1 = pointer
      DO kr=1,jpctr
        DO i=kr+1,jpctr
          ikr = pointer1(i,kr)
          IF (ikr > 0) THEN
            DO j=kr+1,jpctr
              krj = pointer1(kr,j)
              IF (krj > 0) THEN
                ij1 = pointer1(i,j)
! Distinguish whether matrix element is zero or not. If not, proceed
! as in dense case. If it is, create new matrix element.
                IF (ij1 <= 0) THEN
                  total1 = total1 + 1
                  IF (total1 > spfjsize_max) THEN
                    errcode=total1
                    write(6,*) 'Total1 exceeded spfjsize_max: ',total1, &
                                spfjsize_max
                    cmessage='Increase spfjsize_max'
                    CALL EREPORT                                        &
                      ('SETUP_SPFULJAC',errcode,cmessage)
                  END IF
                  pointer1(i,j) = total1
                END IF
              END IF
            END DO
          END IF
        END DO
      END DO

      IF (mype == 0 .AND. printstatus >= prstatus_normal) THEN
        WRITE (6,*) 'Initial and final fill-in: ',total, total1
        WRITE (6,*) 'spfjsize_max is currently set to: ',spfjsize_max
        WRITE (6,*) 'Reduce spfjsize_max if much greater than total1 '
      END IF

      IF (lhook) CALL dr_hook('ASAD_SPARSE_VARS:SETUP_SPFULJAC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE setup_spfuljac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE spfuljac(n_points)
!
!  Routine to calculate the sparse full Jacobian
!  Filter for negatives before entering!
!
! 28/3/2006 Include derivative of SS species with respect to
! ozone (deriv) in calculation of Jacobian
!                                 Olaf Morgenstern
!
! 5/4/2006 Bug removed with indexing of dpd and dpw
!                                 Olaf Morgenstern
! 21/6/2006 Adapted from fuljac to calculate sparse
!           full Jacobian.
!                                 Olaf Morgenstern

        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points

! Local variables
      REAL :: deltt
      REAL :: fr

      CHARACTER (LEN=2) :: ityped

      LOGICAL, SAVE :: first = .TRUE.

      INTEGER :: p
      INTEGER :: irj
      INTEGER :: ifamd
      INTEGER :: itrd
      INTEGER :: i
      INTEGER :: j
      INTEGER :: jc
      INTEGER :: itrcr
      INTEGER :: j3
      INTEGER :: jn
      INTEGER :: js
      INTEGER :: jl
      INTEGER :: i1
      INTEGER :: i2
      INTEGER :: is
      INTEGER :: ij(jpmsp)
      INTEGER :: ik(jpmsp)

      LOGICAL :: lj(jpmsp)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
      IF (lhook) CALL dr_hook('ASAD_SPARSE_VARS:SPFULJAC',zhook_in,zhook_handle)
      deltt=1./cdt
!
      IF (first) THEN
! Determine number and positions of nonzero elements in sparse
! full Jacobian
        CALL setup_spfuljac(n_points)
        first = .FALSE.
      END IF

      spfj = 0.

! Calculate diagonal element of Jacobian
      DO i=1,jpctr
        spfj(:,pointer2(i,i))=-deltt*f(:,i)
      END DO
!
!
!     -----------------------------------------------------------------
!           2.  Calc. full Jacobian matrix.
!               ----- ---- -------- -------
!
! Sum up positive non-fractional, negative, and fractional terms
! Count down (I think that's faster...)
!
      DO p=1,total
        DO i=nposterms(p),1,-1
          spfj(:,p) = spfj(:,p) + prk(:,posterms(p,i))
        END DO
        DO i=nnegterms(p),1,-1
          spfj(:,p) = spfj(:,p) - prk(:,negterms(p,i))
        END DO
        DO i=nfracterms(p),1,-1
          spfj(:,p) = spfj(:,p) + fraction(p,i)*prk(:,fracterms(p,i))
        END DO
      END DO

!  Go through the steady state additions to the Jacobian; currently
!  assume that only O(1D) [and O(3P)] are modelled, and that [both] are
!  required for the O3 loss rate.
!
      DO jc = 1, nstst
        DO j3 = 1,nmsjac(jc)
          irj = nsjac1(j3,jc)
          DO jn=1,jpmsp
            ij(jn)=njcoth(irj,jn)
            lj(jn)=ij(jn) /= 0
          END DO
          DO jn=1,2
            IF (lj(jn)) THEN
              p = pointer2(ij(jn),ntro3)
              spfj(:,p) = spfj(:,p) -                                   &
                prk(:,irj)*deriv(:,jc,1) ! *prkrat(jc)
            END IF
          END DO
          DO jn=3,jpmsp
            IF (lj(jn)) THEN
              p = pointer2(ij(jn),ntro3)
              spfj(:,p) = spfj(:,p) +                                   &
                prk(:,irj)*deriv(:,jc,1) ! *prkrat(jc)
            END IF
          END DO
        END DO
      END DO
!
!     -----------------------------------------------------------------
!          4.  Add deposition terms to Jacobian diagonal.
!              --- ---------- ----- -- -------- ---------
!
      IF ( (ndepw  /=  0) .OR. (ndepd  /=  0) ) THEN
        DO js = 1, jpspec
          ifamd = moffam(js)
          itrd = madvtr(js)
          ityped = ctype(js)
!
          IF ( ifamd /= 0 ) THEN
            DO jl=1,n_points
              IF ((ityped == jpfm) .OR. ((ityped == jpif) .AND.         &
              linfam(jl,itrd))) THEN
                p = pointer2(ifamd,ifamd)
                spfj(jl,p)=spfj(jl,p)                                   &
                - nodd(js)*(dpd(jl,js)+dpw(jl,js))*y(jl,js)
              END IF
            END DO
          END IF
          IF ( itrd /= 0 ) THEN
            p = pointer2(itrd,itrd)
            spfj(:,p) = spfj(:,p) - (dpd(:,js)+dpw(:,js))*y(:,js)
          END IF
        END DO
      END IF
!
!     -------------------------------------------------------------
!          5.  Jacobian elements in final form
!              -------- -------- -- ----- ----
!
      DO p=1,total
        spfj(:,p)=spfj(:,p)/f(:,base_tracer(p))   ! filter f earlier!
      END DO
!
      IF (lhook) CALL dr_hook('ASAD_SPARSE_VARS:SPFULJAC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE spfuljac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE splinslv2(zb,zx,n_points,rafeps,rafbig)

        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE
!
!***** SUB -LINSLV- does simple Gaussian elimination
!*****     fj(N,N)*X(N) = B(N)  BY REDUCING THE A-MATRIX IN PLACE.
!***** THE ENTRY -RESOLV- ASSUMES THAT THE fj-MATRIX HAS BEEN PROPERLY
!*****     REDUCED AND JUST SOLVES FOR X(N).  THIS OPTION SAVES TIME
!*****     WHEN THE SYSTEM IS TO BE RESOLVED WITH A NEW B-VECTOR.
!
! 7/7/2006  Adapted to do sparse matrix algebra and species reordering
!           to improve throughput. For a typical application, this
!           reduces CPU time by at least a factor of 2.
!                                     Olaf Morgenstern
! 1/9/2006  Check diagonal elements of Jacobian whether they are close
!           to 0. Improves stability. Olaf Morgenstern
!
!
! Subroutine interface

      INTEGER, INTENT(IN) :: n_points
      REAL,    INTENT(IN) :: rafeps
      REAL,    INTENT(IN) :: rafbig

      REAL, INTENT(INOUT) :: zb(theta_field_size,jpctr)
      REAL, INTENT(OUT)   :: zx(theta_field_size,jpctr)

! Local variables
      INTEGER :: kr
      INTEGER :: jl
      INTEGER :: i
      INTEGER :: j
      INTEGER :: ikr
      INTEGER :: krj
      INTEGER :: ij

      REAL :: zb1(theta_field_size,jpctr)
      REAL :: zx1(theta_field_size,jpctr)
      REAL :: pivot(theta_field_size)
      REAL :: kfact(theta_field_size)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
! copy pointer variables into new representations; we do not
! want to change these.

      IF (lhook) CALL dr_hook('ASAD_SPARSE_VARS:SPLINSLV2',zhook_in,zhook_handle)
      spfj = MIN(MAX(spfj,-rafbig),rafbig)

      total1 = total
      pointer1 = pointer
      DO kr=1,jpctr
        pivot = spfj(:,pointer1(kr,kr))
        WHERE (ABS(pivot) > rafeps)
          pivot = 1./pivot
        ELSEWHERE
          pivot = rafbig
        ENDWHERE
!        PIVOT = 1./spfj(:,pointer1(kr,kr))
        DO i=kr+1,jpctr
          ikr = pointer1(i,kr)
          IF (ikr > 0) THEN
            KFACT = spfj(:,ikr)*pivot
            spfj(:,ikr) = KFACT
            DO j=kr+1,jpctr
              krj = pointer1(kr,j)
              IF (krj > 0) THEN
                ij = pointer1(i,j)
! Distinguish whether matrix element is zero or not. If not, proceed
! as in dense case. If it is, create new matrix element.
                IF (ij > 0) THEN
                  spfj(:,ij) = spfj(:,ij) - KFACT*spfj(:,krj)
                ELSE
                  total1 = total1 + 1
                  pointer1(i,j) = total1
                  spfj(:,total1) = -KFACT*spfj(:,krj)
                END IF
              END IF
            END DO
          END IF
        END DO
      END DO

! Filter spfj

      spfj = MIN(MAX(spfj,-rafbig),rafbig)

!
!c      call resolv2(zb,zx,n_points,rafeps)!!resolv2 manually inlined below
!
! Form P*zb (permutation of RHS)
      DO i=1,jpctr
        zb1(:,i) = zb(:,ro(i))
      END DO

      DO kr=1,jpctr-1
        DO i=kr+1,jpctr
          ikr = pointer1(i,kr)
          IF (ikr > 0)                                                  &
            zb1(:,i) = zb1(:,i) - spfj(:,ikr) * zb1(:,kr)
        END DO
      END DO
!
      DO kr=jpctr,1,-1
        zx1(:,kr) = zb1(:,kr)
        DO j=kr+1,jpctr
          krj = pointer1(kr,j)
          IF (krj > 0)                                                  &
            zx1(:,kr) = zx1(:,kr) - spfj(:,krj) * zx1(:,j)
        END DO
        pivot = spfj(:,pointer1(kr,kr))
        WHERE (ABS(pivot) < rafeps) pivot = rafeps
        zx1(:,kr) = MIN(MAX(zx1(:,kr),-rafbig),rafbig)/pivot
      END DO

! Form P' zx1 = zx (result of linear system)

      DO j=1,jpctr
        zx(:,ro(j)) = zx1(:,j)
      END DO

      IF (lhook) CALL dr_hook('ASAD_SPARSE_VARS:SPLINSLV2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE splinslv2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE spresolv2(zb,zx,n_points,rafeps)

!            ENTRY FOR BACK SOLUTION WITH DIFFERENT B-VALUE
! 7/7/2006 Modified to perform sparse matrix processing and species
!          reordering.               Olaf Morgenstern
! 17/8/2006 Bug removed              Olaf Morgenstern

        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: n_points
      REAL,    INTENT(IN) :: rafeps
      REAL, INTENT(INOUT) :: zb(theta_field_size,jpctr)
      REAL, INTENT(INOUT) :: zx(theta_field_size,jpctr)

! Local variables
      INTEGER :: kr
      INTEGER :: jl
      INTEGER :: i
      INTEGER :: j
      INTEGER :: krj
      INTEGER :: ikr

      REAL :: zb1(theta_field_size, jpctr)
      REAL :: zx1(theta_field_size, jpctr)
      REAL :: pivot(theta_field_size)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
! Form P*zb (permutation of RHS)
      IF (lhook) CALL dr_hook('ASAD_SPARSE_VARS:SPRESOLV2',zhook_in,zhook_handle)
      DO i=1,jpctr
        zb1(:,i) = zb(:,ro(i))
      END DO

      DO kr=1,jpctr-1
        DO i=kr+1,jpctr
          ikr = pointer1(i,kr)
          IF (ikr > 0) zb1(:,i) = zb1(:,i) - spfj(:,ikr)*zb1(:,kr)
        END DO
      END DO
!
      DO kr=jpctr,1,-1
        zx1(:,kr) = zb1(:,kr)
        DO j=kr+1,jpctr
          krj = pointer1(kr,j)
          IF (krj > 0) zx1(:,kr) = zx1(:,kr) - spfj(:,krj)*zx1(:,j)
        END DO
        pivot = spfj(:,pointer1(kr,kr))
        WHERE(ABS(pivot) < rafeps) pivot = rafeps
        zx1(:,kr) = zx1(:,kr) / pivot
      END DO

      DO i=1,jpctr
        zx(:,ro(i)) = zx1(:,i)
      END DO

      IF (lhook) CALL dr_hook('ASAD_SPARSE_VARS:SPRESOLV2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE spresolv2

      END MODULE asad_sparse_vars
