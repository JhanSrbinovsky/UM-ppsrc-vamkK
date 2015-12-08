! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************
!
! Subroutine Cloud_check.

      Subroutine Cloud_check(                                           &
     &                      j_start, j_end, i_start, i_end              &
     &,                     rows, row_length, wet_model_levels          &
     &,                     halo_i, halo_j                              &
     &,                     qcl, qcf                                    &
     &,                     cloud_fraction_liquid                       &
     &,                     cloud_fraction_frozen                       &
     &,                     area_cloud_fraction                         &
     &,                     bulk_cloud_fraction)

! Purpose:
!        Set appropriate cloud fractions to zero if qcl or qcf are zero
!
! Method:
!          Is described in ;
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: misc
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
              ! model dimensions
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, wet_model_levels                                                &
                             ! number of model levels
     &, halo_i                                                          &
                ! Size of halo in i direction
     &, halo_j                                                          &
                ! Size of halo in j direction
     &, j_start                                                         &
                       ! start j index for processor
     &, j_end                                                           &
                     ! end j index for processor
     &, i_start                                                         &
                       ! start i index for processor
     &, i_end        ! end i index for processor

! Cloud fields
      Real                                                              &
     &  area_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, bulk_cloud_fraction(1-halo_i:row_length+halo_i,                 &
     &                      1-halo_j:rows+halo_j, wet_model_levels)     &
     &, cloud_fraction_liquid(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_model_levels)   &
     &, cloud_fraction_frozen(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_model_levels)

! Primary Arrays
      Real                                                              &
     &  qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &          wet_model_levels)                                       &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &          wet_model_levels)

! Local Variables.

! scalars
      Integer                                                           &
     &  i, j, k           ! Loop indices

      Real                                                              &
     & cloud_tol

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 1.   Initialise set tolerance for testing cloud
!              when qcl and qcf have been updated
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CLOUD_CHECK',zhook_in,zhook_handle)
         cloud_tol = 1.0e-08

! ----------------------------------------------------------------------
! Section 2.   Check if qcl/qcf =0
!              If so set appropriate cloud fraction = 0
! ----------------------------------------------------------------------

         Do k = 1, wet_model_levels
          Do j = j_start, j_end
            Do i = i_start, i_end
             if (qcl(i,j,k) < cloud_tol) then
               cloud_fraction_liquid(i,j,k) = 0.0
             endif      !    qcl(i,j,k) <  cloud_tol
             if (qcf(i,j,k) < cloud_tol) then
               cloud_fraction_frozen(i,j,k) = 0.0
             endif     !    qcf(i,j,k) <  cloud_tol
             if (   cloud_fraction_liquid(i,j,k) <= 0.0                 &
     &        .and. cloud_fraction_frozen(i,j,k) <= 0.0 ) then
                   area_cloud_fraction(i,j,k) = 0.0
                   bulk_cloud_fraction(i,j,k) = 0.0
             endif  ! cloud_fraction = 0
            End Do
          End Do
         End Do

      IF (lhook) CALL dr_hook('CLOUD_CHECK',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Cloud_check
