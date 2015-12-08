! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! To calculate the column integrated relative humidity
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
MODULE column_rh_mod


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  ! Model level heights from centre of Earth
  USE level_heights_mod, ONLY: &
    r_theta_levels             &  ! Radii on theta levels (m) 
   ,r_rho_levels                  ! Radii on rho levels (m)

  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE calc_column_rh(                 &
     npnts, wet_levels, nlcl                 &
     , row_length,rows,model_levels          &
     , index_i, index_j                      &
     , p, T, q, rho_only                     &
     , z_half, l_mixing_ratio, zmax          &
     , qsat, column_rh, column_rh_bl         &
     , column_q, column_q_bl )

    IMPLICIT NONE
    
    INTEGER, INTENT(in) ::     &
      npnts                    & ! Number of points 
    , wet_levels               & ! Number of levels
    , nlcl(npnts)                ! Level of LCL

    INTEGER, INTENT(in) ::     &
      row_length               & ! Local number of points on a row
    , rows                       ! Local number of rows in a theta field
    
    INTEGER, INTENT(in) ::     &
      model_levels               ! Number of model levels

    INTEGER, INTENT(in) ::     &
      index_i(npnts)           & ! Column number of unstable points
    , index_j(npnts)             ! Row number of unstable points
    
    REAL, INTENT(in) ::        &
       p(npnts, wet_levels)    & ! pressure (Pa)
       , T(npnts, wet_levels)  & ! temperature (K)
       , q(npnts, wet_levels)    ! water vapour (kg/kg)

    REAL, INTENT(in) ::        &  
       rho_only(row_length,rows,model_levels) 
                                 ! Density (kg/m3)

    REAL, INTENT(in) ::        &
       z_half(npnts, wet_levels) ! Height on half levels

    REAL, INTENT(in) ::        &
       zmax                      ! Top of column
    
    LOGICAL, INTENT(IN) ::     & 
       l_mixing_ratio            ! .TRUE. if input q a mixing ratio otherwise 
                                 ! assumes input q a specific humidity.

    REAL, INTENT(out) ::       &
       qsat(npnts, wet_levels)   ! q saturation value
       
    REAL, INTENT(out) ::       &
       column_rh(npnts)        & ! Column integrated RH (fraction)
       , column_rh_bl(npnts)     ! Column integrated RH up to LCL

    REAL, INTENT(out) ::       &
       column_q(npnts)         & ! Column integrated q 
       , column_q_bl(npnts)      ! Column integrated q up to LCL
    
    !------------------------------
    ! Local variables
    !------------------------------

    REAL :: column_qsat          ! column_qsat
    REAL :: temp_mass            ! column_qsat

    INTEGER :: k,ii,i,j          ! loop counters

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook( 'COLUMN_RH_MOD:CALC_COLUMN_RH'                   &
                           , zhook_in, zhook_handle )

    DO  k = 1,wet_levels

        ! DEPENDS ON: qsat_mix
        CALL qsat_mix(qsat(:,k),T(:,k),p(:,k),npnts,l_mixing_ratio)

    END DO

    DO ii=1,npnts
      i = index_i(ii)   
      j = index_j(ii) 

      column_q(ii) = 0.0
      column_qsat  = 0.0

      DO k = 2,wet_levels-1

        IF (z_half(ii,k) <= zmax) THEN

          temp_mass = rho_only(i,j,k)                                         &
                    *  r_theta_levels(i,j,k) * r_theta_levels(i,j,k)          & 
                    * (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))

          column_q(ii) = column_q(ii) + temp_mass * q(ii,k)
          column_qsat  = column_qsat  + temp_mass * qsat(ii,k)
          
          IF ( k == nlcl(ii) ) THEN
            column_rh_bl(ii) = column_q(ii)/column_qsat
          END IF
        END IF

      END DO

      column_rh(ii)=column_q(ii)/column_qsat

    END DO
    IF (lhook) CALL dr_hook( 'COLUMN_RH_MOD:CALC_COLUMN_RH'                   &
                           , zhook_out, zhook_handle )
    RETURN


  END SUBROUTINE calc_column_rh

END MODULE column_rh_mod
