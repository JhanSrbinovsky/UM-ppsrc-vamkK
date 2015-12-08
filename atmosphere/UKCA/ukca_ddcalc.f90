! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Combine surface resistance rc with aerodynamic resistance ra and
!  quasi-laminar resistance rb to get overall dry deposition velocity.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_DDEPCTL.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!                     
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
      SUBROUTINE UKCA_DDCALC(row_length, rows, bl_levels, timestep,    &
        dzl, zbl, gsf, ra, rb, rc,                                     &
        nlev_in_bl, zdryrt)  

      USE ASAD_MOD,             ONLY: ndepd
      USE ukca_option_mod,      ONLY: jpdd
      USE nstypes,              ONLY: ntype, npft
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE Control_Max_Sizes
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: row_length, rows
      INTEGER, INTENT(IN) :: bl_levels

      REAL, INTENT(IN) :: timestep
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: zbl              
                                     ! boundary layer depth
      REAL, DIMENSION(row_length, rows, bl_levels), INTENT(IN) :: dzl  
                                     ! thickness of BL levels  
      REAL, DIMENSION(row_length,rows,ntype), INTENT(IN) :: gsf        
                                     ! surface heat flux
      REAL, DIMENSION(row_length,rows,ntype), INTENT(IN) :: ra         
                                     ! aerodynamic resistance
      REAL, DIMENSION(row_length,rows,jpdd), INTENT(IN) :: rb          
                                     ! quasi-laminar resistance
      REAL, DIMENSION(row_length,rows,ntype,jpdd), INTENT(IN) :: rc    
                                     ! surface resistance

      INTEGER, DIMENSION(row_length,rows), INTENT(OUT) :: nlev_in_bl   
                                     ! no of levs in BL
      REAL, DIMENSION(row_length,rows,jpdd), INTENT(OUT) :: zdryrt     
                                     ! dry deposition rate

      INTEGER :: i, j, k, l, n
      REAL :: dd
      REAL, DIMENSION(row_length,rows) :: layer_depth
      REAL, DIMENSION(row_length,rows,ntype,jpdd) :: vd    ! deposn velocity

      REAL :: r_nodep = 1.0e40

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     Set all arrays to zero

      IF (lhook) CALL dr_hook('UKCA_DDCALC',zhook_in,zhook_handle)
      DO j = 1, jpdd
        DO k = 1, rows
          DO i = 1, row_length
            zdryrt(i,k,j) = 0.0
            vd(i,k,1,j) = 0.0
            vd(i,k,2,j) = 0.0
            vd(i,k,3,j) = 0.0
            vd(i,k,4,j) = 0.0
            vd(i,k,5,j) = 0.0
            vd(i,k,6,j) = 0.0
            vd(i,k,7,j) = 0.0
            vd(i,k,8,j) = 0.0
            vd(i,k,9,j) = 0.0
          END DO
        END DO
      END DO

!     Calculate depth of highest model level completely contained
!     within the boundary layer

      DO k = 1, rows
        DO i = 1, row_length
          layer_depth(i,k) = dzl(i,k,1)
          nlev_in_bl(i,k) = 1
        END DO
      END DO

      DO l = 2, bl_levels
        DO k = 1, rows
          DO i = 1, row_length
            dd = layer_depth(i,k) + dzl(i,k,l)
            IF (dd < zbl(i,k)) THEN
              layer_depth(i,k) = dd
              nlev_in_bl(i,k) = l
            END IF
          END DO
        END DO
      END DO

!     Calculate overall dry deposition velocity [vd = 1/(ra + rb + rc)]
!     Do vegetated tiles first. Quasi-laminar resistance pre-multiplied by
!     ln[z0m/z0] = 2.0 for vegetated areas, or 1.0 for smooth surfaces
!     See Ganzeveld & Lelieveld, JGR 1995 Vol.100 No. D10 pp.20999-21012.

      DO k = 1, rows
        DO j = 1, ndepd
          DO n = 1, npft
            DO i = 1, row_length
              IF (rc(i,k,n,j) < r_nodep .AND. gsf(i,k,n) > 0.0) THEN
                vd(i,k,n,j) = 1.0 /                                    &
                  (ra(i,k,n) + 2.0*rb(i,k,j) + rc(i,k,n,j))
              END IF
            END DO
          END DO

!         Now do calculation for non-vegetated tiles

          DO n = npft+1, ntype
            DO i = 1, row_length
              IF (rc(i,k,n,j) < r_nodep .AND. gsf(i,k,n) > 0.0) THEN
                vd(i,k,n,j) = 1.0 /                                    &
                  (ra(i,k,n) + rb(i,k,j) + rc(i,k,n,j))
              END IF
            END DO
          END DO
        END DO
      END DO

!     VD() now contains dry deposition velocities for each tile 
!     in each grid sq. Calculate overall first-order loss rate 
!     over time "timestep" for each tile and sum over all tiles 
!     to obtain overall first-order loss rate zdryrt().

      DO n = 1, ntype
        DO j = 1, ndepd
          DO k = 1, rows
            DO i = 1, row_length
              IF (vd(i,k,n,j) > 0.0) THEN
                zdryrt(i,k,j) = zdryrt(i,k,j) + gsf(i,k,n) *           &
                  (1.0-EXP(-vd(i,k,n,j) * timestep / layer_depth(i,k)))
              END IF
            END DO
          END DO
        END DO
      END DO

!     ZDRYRT() contains loss rate over time "timestep".
!     Divide by timestep to get rate in s-1.

      DO j = 1, ndepd
        DO k = 1, rows
          DO i = 1, row_length
            zdryrt(i,k,j) = -LOG(1.0 - zdryrt(i,k,j)) / timestep
          END DO
        END DO
      END DO

      IF (lhook) CALL dr_hook('UKCA_DDCALC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_DDCALC
