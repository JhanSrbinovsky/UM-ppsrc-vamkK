! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine diffupper
!
! Purpose:
!          This routine used to evaluate diffusion at upper levels only.
!          Calculates conservative horizontal diffusion of a field
!          assuming that all model surfaces are at constant height.
!          This means we can cancel the delta_eta terms and the
!          r-terms can be absorbed into the diffusion coefficient.
! Method:
!          Is described in ;
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      SUBROUTINE eg_diffupper(                                          &
                           field, off_v,                                &
                           levels, active_levels,                       &
                           in_rows, cos_rows, row_length,               &
                           off_x, off_y, halo_i, halo_j,                &
                           sec_latitude, cos_latitude,                  &
                           j_start, j_stop,                             &
                           start_level, up_diff, max_upd_levels,        &
                           csxi2_on, csxi2_off, v_shift, L_u_ns_diff  )

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

! Arguments with INTENT IN. ie: Input variables.

      INTEGER, INTENT(In) ::                                            &
        row_length,                                                     &
                    ! number of point on a row.
        in_rows,                                                        &
                    ! number of rows for field
        cos_rows,                                                       &
                    ! number of rows for cos lat
        levels,                                                         &
                    ! number of levels in field
        max_upd_levels,                                                 &
                      ! max number of levels diffusion applied
        active_levels,                                                  &
                    ! number of levels for diffusion
        start_level,                                                    &
                    ! start level for diffusion
        j_start,                                                        &
                    ! NS-Loop pointer
        j_stop,                                                         &
                      ! NS-Loop pointer
        off_v,                                                          &
                    ! Offset for cos pointer (v=1, u,theta=0)
        halo_i,                                                         &
                    ! Size of halo in i direction.
        halo_j,                                                         &
                    ! Size of halo in j direction.
        off_x,                                                          &
                    ! Size of small halo in i
        off_y       ! Size of small halo in j.
 
      REAL, INTENT(In) ::  up_diff(max_upd_levels) 
                                     ! effective diffusion coefficient

      REAL, INTENT(In) ::                                               &
        sec_latitude (1-off_x:row_length+off_x, 1-off_y:in_rows+off_y), &
        cos_latitude (1-off_x:row_length+off_x, 1-off_y:cos_rows+off_y)

     Integer, Intent(In) :: v_shift

     Logical, Intent(In) :: L_u_ns_diff ! true if diffusing u

     Real,    Intent(In) :: csxi2_on(v_shift-halo_j:                    &
                                     in_rows-1+v_shift+halo_j),         &
                            csxi2_off(1-v_shift-halo_j:                 &
                                      cos_rows-v_shift+halo_j)

! Arguments with INTENT IN/OUT. ie: Input variables changed on Output.

      REAL, INTENT(InOut):: field(1-halo_i:row_length+halo_i,           &
                                  1-halo_j:in_rows+halo_j, levels)

! Local Variables.

      INTEGER                                                           &
        i, j, k, ka,                                                    &
                              ! Loop indices
        info

! Local arrays

     Real  :: csphi_on(v_shift-halo_j:in_rows-1+v_shift+halo_j),        &
              csphi_off(1-v_shift-halo_j:cos_rows-v_shift+halo_j),      &
              u_fac(v_shift-halo_j:in_rows-1+v_shift+halo_j)

      REAL                                                              &
        lambda_term( row_length, in_rows, active_levels ),              &
        phi_term( row_length, in_rows, active_levels ),                 &
        temp(1-off_x:row_length+off_x,1-off_y:in_rows+off_y )

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! No External Routines:

! ----------------------------------------------------------------------
! Section 1. Diffusion constant in wavenumber space
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('EG_DIFFUPPER',zhook_in,zhook_handle)

      DO k = 1, active_levels

        ka = k + start_level - 1

! ----------------------------------------------------------------------
! Section 1.1  Calculate lambda direction term.
! ----------------------------------------------------------------------
        DO j = j_start, j_stop
          DO i = 1-off_x, row_length
            temp(i,j) = ( field(i+1,j,ka) - field(i,j,ka) ) * up_diff(k)
          END DO
          DO i = 1, row_length
            lambda_term(i,j,k) = temp(i,j) - temp(i-1,j)
          END DO
        END DO

! ----------------------------------------------------------------------
! Section 1.2  Calculate phi direction term.
! ----------------------------------------------------------------------
        ! settings for diffusion of u
        If ( L_u_ns_diff ) Then
          Do j = v_shift-halo_j,in_rows-1+v_shift+halo_j
            csphi_on(j) = csxi2_on(j)*csxi2_on(j)
            u_fac(j) = 1.0/csxi2_on(j)
          EndDo
          Do j =1-v_shift-halo_j,cos_rows-v_shift+halo_j
            csphi_off(j) = csxi2_off(j)*csxi2_off(j)*csxi2_off(j)
          EndDo 
        Else
          Do j = v_shift-halo_j,in_rows-1+v_shift+halo_j
            csphi_on(j) = csxi2_on(j)
            u_fac(j) = 1.0
          EndDo
          Do j =1-v_shift-halo_j,cos_rows-v_shift+halo_j
            csphi_off(j) = csxi2_off(j)
          EndDo        
        End If


        DO j = j_start-1, j_stop
          DO i = 1, row_length
            temp(i,j) = ( field(i,j+1,ka)*u_fac(j+1)                   &
                        - field(i,j,ka)  *u_fac(j)  )                  &
                         *csphi_off(j)*up_diff(k)
          END DO
        END DO

        DO j = j_start, j_stop
          DO i = 1, row_length
            phi_term(i,j,k) = (temp(i,j) - temp(i,j-1)) *               &
                                    1.0/csphi_on(j+v_shift-1)
          END DO
        END DO

      END DO ! k = 1, active_levels

! ----------------------------------------------------------------------
! Section 2  Update fields
! ----------------------------------------------------------------------

      DO k = 1, active_levels
        ka = k + start_level - 1

        DO j =j_start, j_stop
          DO i = 1, row_length
            field(i,j,ka) = field(i,j,ka) +                             &
                               lambda_term(i,j,k) + phi_term(i,j,k)
          END DO
        END DO

      END DO ! k = 1, active_levels

      IF (lhook) CALL dr_hook('EG_DIFFUPPER',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eg_diffupper
