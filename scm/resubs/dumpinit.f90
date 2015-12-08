! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE DumpInit
!   PURPOSE:- To initialise primary variables from RESDUMP read in
!             previously from tape

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

SUBROUTINE dumpinit                                                           &
  ( row_length, rows, nprimvars, land_points, nlevs, nwet, nbl_levs           &
  , nsoilt_levs, nsoilm_levs, ntrop, n_cca_lev, land_sea_mask, resdump, u, v  &
  , w, t, theta, q, qcl, qcf, layer_cloud, p, rho, t_deep_soil, smc           &
  , canopy_gb, snodep, tstar, zh, z0msea, cca, rccb, rcct, smcl )

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(In) 
!-----------------------------------------------------------------------------
  INTEGER, INTENT(In) :: &
    row_length           &! x dimension
  , rows                 &! y dimension
  , nlevs                &! Number of levels.
  , nwet                 &! Number of model levels in which Q is set.
  , nbl_levs             &! Number of Boundary layer levels
  , nsoilt_levs          &! Number of soil temperature levels
  , nsoilm_levs          &! Number of soil moisture levels
  , ntrop                &! Max number of levels in the troposphere
  , nprimvars            &! minimum number. of variables required to restart
  , land_points          &! Number of land_points
  , n_cca_lev             ! Number of levels of cca

  REAL, INTENT(In) ::                   &
    resdump(row_length,rows,nprimvars)   ! Contains restart dump

  LOGICAL, INTENT(In) ::                &
    land_sea_mask(row_length,rows)       ! True if land point

!-----------------------------------------------------------------------------
! Arguments with INTENT(Out) 
!-----------------------------------------------------------------------------
  REAL, INTENT(Out) ::                  &
    u(row_length,rows,nlevs)            &! Zonal wind (m/s^2)
  , v(row_length,rows,nlevs)            &! Meridional wind (m/s^2)
  , w(row_length,rows,0:nlevs)          &! vertical velocity      
  , p( row_length,rows,nlevs+1)         &
  , rho(row_length,rows,nlevs)          &
  , t(row_length,rows,nlevs)            &! Temperature at each level
  , theta(row_length,rows,nlevs)        &! Potential temp. (K)      
  , q(row_length,rows,nwet)             &! Specific humidity (kg/kg)
  , layer_cloud(row_length,rows,nwet)   &! Layer cloud amount (decimal)
  , qcl(row_length,rows,nwet)           &! Cloud water content (kg/kg)
  , qcf(row_length,rows,nwet)           &! Cloud ice content (kg/kg)
  , t_deep_soil(land_points,nsoilt_levs)&! Deep soil temperatures
  , smc(land_points)                    &! Soil moisture content (kg/m^2)
  , smcl(land_points,nsoilm_levs)       &! Soil moisture in layers (kg/m^2)
  , rccb(row_length,rows)               &! Convective cloud base
  , rcct(row_length,rows)               &! Convective cloud top
  , cca(row_length,rows,n_cca_lev)      &! Convective cloud amount
  , canopy_gb(land_points)              &! Canopy water content (kg/m^2)
  , snodep(row_length,rows)             &! Snow depth (kg/m^2)
  , tstar(row_length,rows)              &! Surface temperature (K)
  , z0msea(row_length,rows)             &! Sea surface roughness length 
  , zh(row_length,rows)                  ! Height above surface of to
                                         ! boundary layer (m)


!-----------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------------------------------

  INTEGER ::           &
    i,j,k              &! Loop counter
  , icount             &! Counter
  , land_cnt

!-----------------------------------------------------------------------------

  land_cnt = 0

  DO k=1, rows
    DO j=1, row_length

      IF (land_sea_mask(j,k)) THEN
        land_cnt = land_cnt + 1
      END IF

      DO i=1, nlevs
        u(j,k, i) = resdump(j,k, i)
      END DO

      icount = i
      DO i=icount, icount + nlevs-1
        v(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = i
      DO i=icount, icount + nlevs
        w(j,k, i-icount) = resdump(j,k, i)
      END DO

      icount = i
      DO i=icount, icount + nlevs-1
        t(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = i
      DO i=icount, icount + nlevs-1
        theta(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = i
      DO i=icount, icount + nwet-1
        q(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = i
      DO i=icount, icount + nwet-1
        qcl(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = i
      DO i=icount, icount + nwet-1
        qcf(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = i
      DO i=icount, icount + nwet-1
        layer_cloud(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = i
      IF (land_sea_mask(j,k)) THEN
        DO i=icount, icount + nsoilt_levs-1
          t_deep_soil(land_cnt, i-icount + 1) = resdump(j,k, i)
        END DO
        icount = i
      END IF

      DO i=icount, icount + nlevs+1
        p(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = i
      DO i=icount, icount + nlevs
        rho(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = i
      IF (land_sea_mask(j,k)) THEN
        smc(land_cnt) = resdump(j,k, icount)
        icount = icount + 1
        canopy_gb(land_cnt) = resdump(j,k, icount)
        icount = icount + 1
      END IF

      snodep(j,k) = resdump(j,k, icount)

      icount = icount + 1
      tstar(j,k) = resdump(j,k, icount)

      icount = icount + 1
      zh(j,k) = resdump(j,k, icount)

      icount = icount + 1
      z0msea(j,k) = resdump(j,k, icount)

      icount = icount + 1
      DO i=icount, icount + n_cca_lev
        cca(j,k, i-icount + 1) = resdump(j,k, i)
      END DO

      icount = icount + 1
      rccb(j,k) = resdump(j,k, icount)

      icount = icount + 1
      rcct(j,k) = resdump(j,k, icount)

      IF (land_sea_mask(j,k)) THEN
        DO i=1, nsoilm_levs
          smcl(land_cnt, i) = resdump(j,k, icount + i)
        END DO
      END IF

    END DO                     ! j
  END DO                     ! k

  RETURN
END SUBROUTINE dumpinit

