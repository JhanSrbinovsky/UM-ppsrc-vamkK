! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculate the interpolation coefficients for spatial degradation

!  Purpose: To calculate the coefficients used to interpolate the
!           radiation quantities from neighbouring grid boxes to
!           a grid box where the radiation calculations were not
!           done on a radiation timestep.

!  Method:  In a chequer-board pattern, a box which has not done a
!           radiation calculation is surrounded on the north, south,
!           east and west by boxes which have calculated the radiation
!           quantities therefore a coefficient of 0.25 for these four
!           points would seem to be sufficient. However, the land-sea
!           contrast can sometimes be so sharp that we would not wish
!           to interpolate from one grid box to another if their
!           surface types do not match. Furthermore, at the boundary
!           of the whole model domain, halo points do not contain
!           meaningful values, so we do not use them.

!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Radiation Control

      SUBROUTINE coeffs_degrade(es_space_interp, land, first_row,       &
                                last_row, model_domain, mt_lam,         &
                                at_top_of_lpg, at_base_of_lpg,          &
                                at_left_of_lpg,at_right_of_lpg,         &
                                p_field, row_length, rows, offx, offy)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: p_field
!                              Full length of data on a PE
      INTEGER, INTENT(IN) :: row_length
!                              No of points per row
      INTEGER, INTENT(IN) :: first_row
!                              first row of non-polar data
      INTEGER, INTENT(IN) :: last_row
!                              last row of non-polar data
      INTEGER, INTENT(IN) :: model_domain

      INTEGER, INTENT(IN) :: mt_lam

      INTEGER, INTENT(IN) :: offx

      INTEGER, INTENT(IN) :: offy

      INTEGER, INTENT(IN) :: rows

      LOGICAL, INTENT(IN) :: land(1-offx:row_length+offx,               &
                                                   1-offy:rows+offy)
!            If true, then grid box is a land point
      LOGICAL, INTENT(IN) ::  at_top_of_lpg
!              Logical variable indicating if this PE at edge of LPG
      LOGICAL, INTENT(IN) ::  at_right_of_lpg
!              Logical variable indicating if this PE at edge of LPG
      LOGICAL, INTENT(IN) ::  at_base_of_lpg
!              Logical variable indicating if this PE at edge of LPG
      LOGICAL, INTENT(IN) ::  at_left_of_lpg
!              Logical variable indicating if this PE at edge of LPG

      REAL, INTENT(OUT) :: es_space_interp(4, row_length, rows)
!                The coefficients for radiation interpolation

!     Local variables

      INTEGER :: i              ! loop variable
      INTEGER :: j              ! loop variable
      INTEGER :: k              ! loop variable
      INTEGER :: row_length_tot ! Total row_length including halos

      LOGICAL :: surf_match(4, row_length, rows)
!                   True if central grid box and its neighbour
!                   have same surface type

      REAL :: fac1              ! Factors used for
      REAL :: fac2              ! computation of final
      REAL :: fac3              ! coefficients
      REAL :: fac4              !
      REAL :: num1              ! Temporary storage
      REAL :: tempnum           ! Temporary storage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! externals: none

      IF (lhook) CALL dr_hook('COEFFS_DEGRADE',zhook_in,zhook_handle)
      DO j=first_row,last_row
        DO i=1,row_length
          DO k=1,4
            es_space_interp(k,i,j)=0.
          END DO
        END DO
      END DO

!  Firstly, set all points assuming there are no boundary problems on
!  the whole domain.

      DO j=first_row,last_row
        DO i=1,row_length
          surf_match(1,i,j)=(land(i,j).EQV.land(i,j+1))
          surf_match(2,i,j)=(land(i,j).EQV.land(i+1,j))
          surf_match(3,i,j)=(land(i,j).EQV.land(i,j-1))
          surf_match(4,i,j)=(land(i,j).EQV.land(i-1,j))
        END DO
      END DO

!  Now the problem at the edges of the whole model domain are sorted
!  out. Note that global and non-global runs must be treated
!  differently - the latter type have a problem at the left and right
!  edges of the whole model domain, whereas global models only have
!  problems at the top and bottom of the whole domain.

      IF (at_top_of_lpg) THEN
        DO i=1, row_length
          surf_match(1,i,last_row)=.FALSE.
        END DO
      END IF

      IF (at_base_of_lpg) THEN
        DO i=1, row_length
          surf_match(3,i,first_row)=.FALSE.
        END DO
      END IF

      IF (model_domain  ==  mt_lam) THEN
        IF (at_right_of_lpg) THEN
          DO j=first_row,last_row
            surf_match(2,row_length,j)=.FALSE.
          END DO
        END IF

        IF (at_left_of_lpg) THEN
          DO j=first_row,last_row
            surf_match(4,1,j)=.FALSE.
          END DO
        END IF
      END IF

! Now calculate the coefficients for interpolation.

      DO j = first_row,last_row
        DO i = 1,row_length
          fac1=0.
          fac2=0.
          fac3=0.
          fac4=0.
          IF (surf_match(1,i,j)) THEN
            fac1=1.0
          END IF
          IF (surf_match(2,i,j)) THEN
            fac2=1.0
          END IF
          IF (surf_match(3,i,j)) THEN
            fac3=1.0
          END IF
          IF (surf_match(4,i,j)) THEN
            fac4=1.0
          END IF
          tempnum=fac1+fac2+fac3+fac4
          IF (tempnum <  0.1) THEN
            fac1=1.
            fac2=1.
            fac3=1.
            fac4=1.
          END IF
          num1=1./(fac1+fac2+fac3+fac4)
          es_space_interp(1,i,j)=fac1*num1
          es_space_interp(2,i,j)=fac2*num1
          es_space_interp(3,i,j)=fac3*num1
          es_space_interp(4,i,j)=fac4*num1
        END DO
      END DO

      IF (lhook) CALL dr_hook('COEFFS_DEGRADE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE coeffs_degrade
