! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine IDL_Set_Ancil

      Subroutine IDL_Set_Ancil(                                         &
     &  row_length, rows, n_rows, model_levels                          &
     &, wet_model_levels, n_cca_levels                                  &
     &, halo_i,halo_j                                                   &
     &, cumulus, nbdsc, ntdsc                                           &
     &, cca, ccb, cct, cclwp                                            &
     &, tstar, land_sea_mask                                            &
     &, SW_incs, LW_incs                                                &
     &, t1_sd, q1_sd, zh                                                &
     &, area_cloud_fraction, bulk_cloud_fraction                        &
     &, cloud_fraction_liquid, cloud_fraction_frozen                    &
     &, ti, z0msea, ntml, u_0, v_0, u_0_p, v_0_p                        &
     &, theta, tstar_tile, ntiles, land_field, land_index               &
     &, L_code_test, rad_hr, theta_surface, bl_levels  )

! Purpose:
!             set idealised ancillaries
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, n_cca_levels                                                    &
     &, bl_levels                                                       &
     &, halo_i                                                          &
     &, halo_j

      Integer ntiles      ! Number of tiles in MOSESII
      Integer land_field  ! Number of land points
      Integer land_index(land_field)

      Logical                                                           &
     &  cumulus(row_length, rows)                                       &
                                  ! bl convection flag
     &, land_sea_mask(row_length, rows)                                 &
     &, L_code_test   ! User switch

      Real                                                              &
     &  SW_incs(row_length, rows, 0:model_levels+1)                     &
     &, LW_incs(row_length, rows, 0:model_levels)                       &
     &, t1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, q1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, ti(row_length, rows)                                            &
                              ! set equal to tstar
     &, zh(row_length, rows)                                            &
                             ! boundary layer height
     &, rad_hr(row_length, rows, bl_levels,2)                           &
                                                ! BL rad heating rates
     &, theta_surface    ! surface theta

      Real, Intent(InOut) :: theta(1-halo_i:row_length+halo_i,          &
     &                             1-halo_j:rows+halo_j, model_levels)

! Convection
      Real                                                              &
     &  cca (row_length, rows, n_cca_levels)                            &
     &, cclwp(row_length, rows) ! condensed water path (KG/M**2)

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)

      Integer                                                           &
     &  ntml(row_length, rows)                                          &
     &, nbdsc(row_length, rows)                                         &
     &, ntdsc(row_length, rows)

! Diagnostic variables

      Real                                                              &
     &  area_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, bulk_cloud_fraction(1-halo_i:row_length+halo_i,                 &
     &                      1-halo_j:rows+halo_j, wet_model_levels)     &
     &, cloud_fraction_liquid(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_model_levels)   &
     &, cloud_fraction_frozen(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_model_levels)

      Real                                                              &
     &  tstar(row_length, rows)                                         &
     &, tstar_tile(land_field,ntiles)

      Real                                                              &
     &  u_0(row_length, rows)                                           &
                                ! set to zero
     &, v_0(row_length, n_rows)                                         &
                                ! set to zero
     &, u_0_p(row_length, rows)                                         &
                                  ! set to zero
     &, v_0_p(row_length, rows) ! set to zero

      Real                                                              &
     &  z0msea (row_length, rows) ! veg/qrparm.veg.rough

! Locals
      Integer i,j,k
      Integer l1,l2
      Logical L_const_tsurface   ! Switch for constant surface temperat.
      Logical L_sea

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_SET_ANCIL',zhook_in,zhook_handle)
      do j = 1, rows
        do i = 1, row_length
          cumulus(i,j)=.false.
          nbdsc(i,j)=1
          ntdsc(i,j)=1
          t1_sd(i,j)=0.0
          q1_sd(i,j)=0.0
          zh(i,j)=100.0
          cclwp(i,j)=0.0
          ccb(i,j)=1
          cct(i,j)=1
          ti(i,j)=tstar(i,j)
          z0msea(i,j) = 1.0e-4
          ntml(i,j)=1
          u_0(i,j)=0.0
          u_0_p(i,j)=0.0
          v_0_p(i,j)=0.0
        enddo
      enddo

      ! Set surface temperature on land tiles and sea to theta_surface
      ! or level 1 theta (i.e. first layer is isentropic)

      ! At present, hardwire L_const_tsurface to .true.
      L_const_tsurface = .true.

      If (L_const_tsurface) Then

        Do j = 1, rows
          Do i = 1, row_length
            tstar(i,j) = theta_surface
          End Do
        End Do

        Do l2 = 1,ntiles
          Do l1 = 1, land_field
            j=(land_index(l1)-1)/row_length + 1
            i=land_index(l1) - (j-1)*row_length
            tstar_tile(l1,l2) = theta_surface
          End Do
        End Do

      Else ! Surface temperature = level 1 theta

        Do j = 1, rows
          Do i = 1, row_length
            tstar(i,j)=theta(i,j,1)
          End Do
        End Do

        Do l2 = 1,ntiles
          Do l1 = 1, land_field
            j=(land_index(l1)-1)/row_length + 1
            i=land_index(l1) - (j-1)*row_length
            tstar_tile(l1,l2) = theta(i,j,1)
          End Do
        End Do

      End If ! on L_const_tsurface

      do j = 1, n_rows
        do i = 1, row_length
          v_0(i,j)=0.0
        enddo
      enddo

      do k = 1, n_cca_levels
        do j = 1, rows
          do i = 1, row_length
            cca(i,j,k) = 0.0
          enddo
        enddo
      enddo
      do k = 1, wet_model_levels
        do j = 1, rows
          do i = 1 , row_length
            area_cloud_fraction(i,j,k) = 0.0
            bulk_cloud_fraction(i,j,k) = 0.0
            cloud_fraction_liquid(i,j,k) = 0.0
            cloud_fraction_frozen(i,j,k) = 0.0
          enddo
        enddo
      enddo

      do k = 0, model_levels
        do j = 1, rows
          do i = 1, row_length
            SW_incs(i,j,k) = 0.0
            LW_incs(i,j,k) = 0.0
          enddo
        enddo
      enddo
      do j = 1, rows
        do i = 1, row_length
          SW_incs(i,j,model_levels+1) = 0.0
        enddo
      enddo

      do k = 1, bl_levels
        do j = 1, rows
          do i = 1, row_length
            rad_hr(i,j,k,1) = 0.0
            rad_hr(i,j,k,2) = 0.0
          enddo
        enddo
      enddo

      IF (lhook) CALL dr_hook('IDL_SET_ANCIL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_Set_Ancil
