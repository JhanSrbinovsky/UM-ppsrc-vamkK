! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine for the aerosol climatology for NWP
!
! Purpose:
!  Copy individual fields of the aerosol climatology for NWP into
!  a single array. Aerosol species that are not requested are not
!  copied. The array index corresponding to each requested component
!  are stored for latter access to the corresponding field.
!
! Method:
!   Straightforward.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Description of code:
!   FORTRAN 90
!- ---------------------------------------------------------------------
MODULE set_arcl_clim_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE set_arcl_clim(                                                      &
                                ! Array dimensions
  row_length,                                                                  &
  rows,                                                                        &
  model_levels,                                                                &
  n_arcl_compnts,                                                              &
                                ! Internal model switches
  l_use_arcl,                                                                  &
                                ! Climatologies from ancillary files:
                                !    biomass-burning
  arclbiom_fr, arclbiom_ag, arclbiom_ic,                                       &
                                !    black-carbon
  arclblck_fr, arclblck_ag,                                                    &
                                !    sea-salt
  arclsslt_fi, arclsslt_jt,                                                    &
                                !    sulphate
  arclsulp_ac, arclsulp_ak, arclsulp_di,                                       &
                                !    mineral dust
  arcldust_b1, arcldust_b2, arcldust_b3,                                       &
  arcldust_b4, arcldust_b5, arcldust_b6,                                       &
                                !    fossil-fuel organic carbon
  arclocff_fr, arclocff_ag, arclocff_ic,                                       &
                                !    delta aerosol
  arcldlta_dl,                                                                 &
                                ! Internal climatology array
  arcl,                                                                        &
                                ! Component array indices
  i_arcl_compnts                                                               &
  )

USE arcl_mod,        ONLY: npd_arcl_compnts, npd_arcl_species, ip_arcl_sulp,   &
                           ip_arcl_dust, ip_arcl_sslt, ip_arcl_blck,           &
                           ip_arcl_biom, ip_arcl_ocff, ip_arcl_dlta,           &
                           ip_arcl_sulp_ac, ip_arcl_sulp_ak, ip_arcl_sulp_di,  &
                           ip_arcl_dust_b1, ip_arcl_dust_b2, ip_arcl_dust_b3,  &
                           ip_arcl_dust_b4, ip_arcl_dust_b5, ip_arcl_dust_b6,  &
                           ip_arcl_sslt_fi, ip_arcl_sslt_jt, ip_arcl_blck_fr,  &
                           ip_arcl_blck_ag, ip_arcl_biom_fr, ip_arcl_biom_ag,  &
                           ip_arcl_biom_ic, ip_arcl_ocff_fr, ip_arcl_ocff_ag,  &
                           ip_arcl_ocff_ic, ip_arcl_dlta_dl

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!
! Arguments with intent(in)
!

!
! Array dimensions
!
INTEGER, INTENT(IN) ::                                                         &
  row_length,                                                                  &
  rows,                                                                        &
  model_levels,                                                                &
  n_arcl_compnts

!
! Internal model switches
!
LOGICAL, DIMENSION(npd_arcl_species), INTENT(IN) :: l_use_arcl

!
! Climatologies from ancillary files
!
!    Three components of biomass-burning
!
REAL, DIMENSION(row_length, rows, model_levels), INTENT(IN) ::                 &
  arclbiom_fr,                                                                 &
  arclbiom_ag,                                                                 &
  arclbiom_ic
!
!    Two components of black-carbon
!
REAL, DIMENSION(row_length, rows, model_levels), INTENT(IN) ::                 &
  arclblck_fr,                                                                 &
  arclblck_ag
!
!    Two components of sea-salt
!
REAL, DIMENSION(row_length, rows, model_levels), INTENT(IN) ::                 &
  arclsslt_fi,                                                                 &
  arclsslt_jt
!
!    Three components of sulphate
!
REAL, DIMENSION(row_length, rows, model_levels), INTENT(IN) ::                 &
  arclsulp_ac,                                                                 &
  arclsulp_ak,                                                                 &
  arclsulp_di
!
!    Six components of mineral dust
!
REAL, DIMENSION(row_length, rows, model_levels), INTENT(IN) ::                 &
  arcldust_b1,                                                                 &
  arcldust_b2,                                                                 &
  arcldust_b3,                                                                 &
  arcldust_b4,                                                                 &
  arcldust_b5,                                                                 &
  arcldust_b6

!
!    Three components of fossil-fuel organic carbon
!
REAL, DIMENSION(row_length, rows, model_levels), INTENT(IN) ::                 &
  arclocff_fr,                                                                 &
  arclocff_ag,                                                                 &
  arclocff_ic

!
!    One component of delta aerosol
!
REAL, DIMENSION(row_length, rows, model_levels), INTENT(IN) ::                 &
  arcldlta_dl

!
! Arguments with intent(out)
!

!
! Internal climatology array
!
REAL, DIMENSION(row_length, rows, model_levels, n_arcl_compnts),               &
  INTENT(OUT) :: arcl

!
! Component array indices
!
INTEGER, DIMENSION(npd_arcl_compnts), INTENT(OUT) :: i_arcl_compnts

!
! Local variables
!

INTEGER i_cmp

INTEGER i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('SET_ARCL_CLIM',zhook_in,zhook_handle)
i_cmp = 1 ! since this routine has been called, there is at least
! one component to process.

!
! For each requested species, copy the corresponding components
! into the gathering array.
! We also keep track of the index where component has been put.
! An index of -1 denotes components that are not used.
!
IF (l_use_arcl(ip_arcl_sulp)) THEN

  i_arcl_compnts(ip_arcl_sulp_ac) = i_cmp
  i_arcl_compnts(ip_arcl_sulp_ak) = i_cmp + 1
  i_arcl_compnts(ip_arcl_sulp_di) = i_cmp + 2

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        arcl(i, j, k, i_cmp  ) = arclsulp_ac(i, j, k)
        arcl(i, j, k, i_cmp+1) = arclsulp_ak(i, j, k)
        arcl(i, j, k, i_cmp+2) = arclsulp_di(i, j, k)

      END DO ! i
    END DO ! j
  END DO ! k

  i_cmp = i_cmp + 3

ELSE

  i_arcl_compnts(ip_arcl_sulp_ac) = -1
  i_arcl_compnts(ip_arcl_sulp_ak) = -1
  i_arcl_compnts(ip_arcl_sulp_di) = -1

END IF

IF (l_use_arcl(ip_arcl_dust)) THEN

  i_arcl_compnts(ip_arcl_dust_b1) = i_cmp
  i_arcl_compnts(ip_arcl_dust_b2) = i_cmp + 1
  i_arcl_compnts(ip_arcl_dust_b3) = i_cmp + 2
  i_arcl_compnts(ip_arcl_dust_b4) = i_cmp + 3
  i_arcl_compnts(ip_arcl_dust_b5) = i_cmp + 4
  i_arcl_compnts(ip_arcl_dust_b6) = i_cmp + 5

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        arcl(i, j, k, i_cmp  ) = arcldust_b1(i, j, k)
        arcl(i, j, k, i_cmp+1) = arcldust_b2(i, j, k)
        arcl(i, j, k, i_cmp+2) = arcldust_b3(i, j, k)
        arcl(i, j, k, i_cmp+3) = arcldust_b4(i, j, k)
        arcl(i, j, k, i_cmp+4) = arcldust_b5(i, j, k)
        arcl(i, j, k, i_cmp+5) = arcldust_b6(i, j, k)

      END DO ! i
    END DO ! j
  END DO ! k

  i_cmp = i_cmp + 6

ELSE

  i_arcl_compnts(ip_arcl_dust_b1) = -1
  i_arcl_compnts(ip_arcl_dust_b2) = -1
  i_arcl_compnts(ip_arcl_dust_b3) = -1
  i_arcl_compnts(ip_arcl_dust_b4) = -1
  i_arcl_compnts(ip_arcl_dust_b5) = -1
  i_arcl_compnts(ip_arcl_dust_b6) = -1

END IF

IF (l_use_arcl(ip_arcl_sslt)) THEN

  i_arcl_compnts(ip_arcl_sslt_fi) = i_cmp
  i_arcl_compnts(ip_arcl_sslt_jt) = i_cmp + 1

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        arcl(i, j, k, i_cmp  ) = arclsslt_fi(i, j, k)
        arcl(i, j, k, i_cmp+1) = arclsslt_jt(i, j, k)

      END DO ! i
    END DO ! j
  END DO ! k

  i_cmp = i_cmp + 2

ELSE

  i_arcl_compnts(ip_arcl_sslt_fi) = -1
  i_arcl_compnts(ip_arcl_sslt_jt) = -1

END IF

IF (l_use_arcl(ip_arcl_blck)) THEN

  i_arcl_compnts(ip_arcl_blck_fr) = i_cmp
  i_arcl_compnts(ip_arcl_blck_ag) = i_cmp + 1

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        arcl(i, j, k, i_cmp  ) = arclblck_fr(i, j, k)
        arcl(i, j, k, i_cmp+1) = arclblck_ag(i, j, k)

      END DO ! i
    END DO ! j
  END DO ! k

  i_cmp = i_cmp + 2

ELSE

  i_arcl_compnts(ip_arcl_blck_fr) = -1
  i_arcl_compnts(ip_arcl_blck_ag) = -1

END IF

IF (l_use_arcl(ip_arcl_biom)) THEN

  i_arcl_compnts(ip_arcl_biom_fr) = i_cmp
  i_arcl_compnts(ip_arcl_biom_ag) = i_cmp + 1
  i_arcl_compnts(ip_arcl_biom_ic) = i_cmp + 2

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        arcl(i, j, k, i_cmp  ) = arclbiom_fr(i, j, k)
        arcl(i, j, k, i_cmp+1) = arclbiom_ag(i, j, k)
        arcl(i, j, k, i_cmp+2) = arclbiom_ic(i, j, k)

      END DO ! i
    END DO ! j
  END DO ! k

  i_cmp = i_cmp + 3

ELSE

  i_arcl_compnts(ip_arcl_biom_fr) = -1
  i_arcl_compnts(ip_arcl_biom_ag) = -1
  i_arcl_compnts(ip_arcl_biom_ic) = -1

END IF

IF (l_use_arcl(ip_arcl_ocff)) THEN

  i_arcl_compnts(ip_arcl_ocff_fr) = i_cmp
  i_arcl_compnts(ip_arcl_ocff_ag) = i_cmp + 1
  i_arcl_compnts(ip_arcl_ocff_ic) = i_cmp + 2

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        arcl(i, j, k, i_cmp  ) = arclocff_fr(i, j, k)
        arcl(i, j, k, i_cmp+1) = arclocff_ag(i, j, k)
        arcl(i, j, k, i_cmp+2) = arclocff_ic(i, j, k)

      END DO ! i
    END DO ! j
  END DO ! k

  i_cmp = i_cmp + 3

ELSE

  i_arcl_compnts(ip_arcl_ocff_fr) = -1
  i_arcl_compnts(ip_arcl_ocff_ag) = -1
  i_arcl_compnts(ip_arcl_ocff_ic) = -1

END IF

IF (l_use_arcl(ip_arcl_dlta)) THEN

  i_arcl_compnts(ip_arcl_dlta_dl) = i_cmp

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        arcl(i, j, k, i_cmp  ) = arcldlta_dl(i, j, k)

      END DO ! i
    END DO ! j
  END DO ! k

  i_cmp = i_cmp + 1

ELSE

  i_arcl_compnts(ip_arcl_dlta_dl) = -1

END IF
IF (lhook) CALL dr_hook('SET_ARCL_CLIM',zhook_out,zhook_handle)
RETURN

END SUBROUTINE set_arcl_clim
END MODULE set_arcl_clim_mod
