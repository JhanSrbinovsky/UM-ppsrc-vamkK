! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! 
! Calculate the scaled divergence of the Eliassen-Palm flux from
! u, v, w and temperature on pressure levels.
! Residual mean meridional circulation, Meridional and Vertical
! components of Eliassen-Palm, Meridional heat flux and momentum 
! fluxes are returned also.
! 
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Climate Diagnostics
MODULE calc_div_ep_flux_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE Calc_Div_EP_Flux(                                            &
    rows, n_rows, row_length, global_row_length, sin_v_latitude,        &
    qvstarbar_p, qwstarbar_p, qepy_p, qepz_p, qdivep_p,                 &
    qzvptp_p, qzupvp_p,                                                 &
    u_press, u_p_levs, v_p_levs, w_p_levs, t_p_levs,                    &
    vstarbar_p_levs, wstarbar_p_levs,  Fy_p_levs, Fz_p_levs,            &
    divF_p_levs, zvptp_p_levs, zupvp_p_levs,                            &
    u_p, v_p, w_p, temperature,                                         &
    vstarbar, wstarbar, Fy, Fz, divF, zvptp, zupvp)

    USE dynamics_grid_mod, ONLY: l_vatpoles

    USE earth_constants_mod, ONLY: g, earth_radius, omega

    USE atmos_constants_mod, ONLY: kappa, p_zero, sclht

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    USE atm_fields_bounds_mod, ONLY:            &
                     udims, udims_s, vdims, vdims_s 

    USE UM_ParVars
    USE zonal_average_mod, ONLY: zonal_average
    IMPLICIT NONE

    INTEGER, Intent(in) :: rows, n_rows, row_length, global_row_length 

    ! No of levels for output of u_p, v_p, w_p and temperature
    INTEGER, intent(in) :: u_p_levs, v_p_levs, w_p_levs, t_p_levs


    ! No of levels for residual circulation
    INTEGER, Intent(IN) :: vstarbar_p_levs, wstarbar_p_levs,            &
                           Fy_p_levs, Fz_p_levs, divF_p_levs,           &
                           zvptp_p_levs, zupvp_p_levs

    ! Surface density at reference latitude
    REAL, Parameter  :: rho_surface = 1.212 

    REAL sin_v_latitude(vdims%i_start:vdims%i_end,                      &
                        vdims%j_start:vdims%j_end)   

    REAL zlogc(u_p_levs)  ! "log-pressure" coordinate

    REAL, Intent(in) ::                                                 &
         u_p(udims%i_start:udims%i_end,                                 &
             vdims%j_start:vdims%j_end,u_p_levs)                        &
                                           ! u at selected pressures
        ,v_p(udims%i_start:udims%i_end,                                 &
             vdims%j_start:vdims%j_end,v_p_levs)                        &
                                           ! v at selected pressures
        ,w_p(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,       &
             w_p_levs)                                                  
                                           ! w at selected pressures     

    REAL, Intent(in) :: u_press(u_p_levs) 

    LOGICAL, Intent(in) ::                                              &
         qvstarbar_p                                                    &
                 ! Flag for residual mean meridional circulation
        ,qwstarbar_p                                                    &
                 ! Flag for residual mean meridional circulation
        ,qepy_p                                                         &
                 ! Flag for Eliassen Palm flux (phi component)
        ,qepz_p                                                         &
                 ! Flag for Eliassen Palm flux (vertical component)
        ,qdivep_p                                                       &
                 ! Flag for divergence of Eliassen Palm flux
        ,qzvptp_p                                                       &
                 ! Flag for meridional heat flux
        ,qzupvp_p                                                       
                 ! Flag for meridional momentum flux

    INTEGER  start, end

    REAL     uv(udims%i_start:udims%i_end,                              &
                vdims%j_start:vdims%j_end, u_p_levs)                    &
           , theta_p(udims%i_start:udims%i_end,                         &
                     vdims%j_start:vdims%j_end, u_p_levs)               &
           , vT(udims%i_start:udims%i_end,                              &
                vdims%j_start:vdims%j_end, u_p_levs)                    &
           , vtheta(udims%i_start:udims%i_end,                          &
                    vdims%j_start:vdims%j_end, u_p_levs)                &      
           , temperature(udims%i_start:udims%i_end,                     &
                         vdims%j_start:vdims%j_end,  t_p_levs)          &
           , uw(udims%i_start:udims%i_end,                              &
                vdims%j_start:vdims%j_end, u_p_levs)                    &
           , zonal_u_halo(vdims_s%j_start:vdims_s%j_end, u_p_levs)      &
           , zonal_u(vdims%j_start:vdims%j_end, u_p_levs)               &
           , zonal_v(vdims%j_start:vdims%j_end, u_p_levs)               &
           , zonal_uv(vdims%j_start:vdims%j_end, u_p_levs)              &
           , zonal_theta(vdims%j_start:vdims%j_end, u_p_levs)           &
           , zonal_t(vdims%j_start:vdims%j_end, u_p_levs)               &
           , zonal_vT(vdims%j_start:vdims%j_end, u_p_levs)              &
           , zonal_vtheta(vdims%j_start:vdims%j_end, u_p_levs)          &
           , zonal_w(vdims%j_start:vdims%j_end, u_p_levs)               &
           , zonal_uw(vdims%j_start:vdims%j_end, u_p_levs)              &
           , zvpthp(vdims_s%j_start:vdims_s%j_end, u_p_levs)            &
           , zupwp(vdims%j_start:vdims%j_end, u_p_levs)                 &
           , dthdz(vdims_s%j_start:vdims_s%j_end, u_p_levs)             &
           , phi(vdims_s%j_start:vdims_s%j_end)                         &
           , rho0(u_p_levs)                                             &
           , dudz(vdims%j_start:vdims%j_end, u_p_levs)                  &
           , ducosdphi(vdims%j_start:vdims%j_end, u_p_levs)             &
           , divFz(vdims%j_start:vdims%j_end, u_p_levs)                 &
           , divFy(vdims%j_start:vdims%j_end, u_p_levs)                 &
           , Fy_halo(vdims_s%j_start:vdims_s%j_end, u_p_levs)

    REAL, Intent(out) ::                                                &    
         vstarbar(vdims%j_start:vdims%j_end, u_p_levs)                  &
        ,wstarbar(vdims%j_start:vdims%j_end, u_p_levs)                  &
        ,divF(vdims%j_start:vdims%j_end, u_p_levs)                      &
        ,Fy(vdims%j_start:vdims%j_end, u_p_levs)                        &
        ,Fz(vdims%j_start:vdims%j_end, u_p_levs)                        &
        ,zvptp(vdims%j_start:vdims%j_end, u_p_levs)                     &
        ,zupvp(vdims%j_start:vdims%j_end, u_p_levs)

    INTEGER :: x,y,z
    INTEGER :: top

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! STASH items 310 311: residual circulation on pressure surfaces/B grid
! ----------------------------------------------------------------------

    IF (lhook) CALL dr_hook('CALC_DIV_EP_FLUX',zhook_in,zhook_handle)

    IF(qvstarbar_p .OR. qwstarbar_p .OR. qepy_p .OR. qepz_p &
         .OR. qdivep_p .OR. qzvptp_p .OR. qzupvp_p) THEN
!
! Calculate the "log-pressure" coordinate z=-Hln(p/ps) and basic density
! from surface density at reference latitude.
!
      DO z=1,u_p_levs
        zlogc(z) = - sclht * ALOG(u_press(z)/(p_zero*0.01))
        rho0(z)      = rho_surface * u_press(z) / (p_zero*0.01)
      END DO

      DO z=1,u_p_levs
        DO y= vdims%j_start,vdims%j_end      ! latitude loop
          DO x=udims%i_start,udims%i_end     ! longitude loop
            theta_p(x,y,z) = temperature(x,y,z) *                     &
                                ((p_zero*0.01) / u_press(z))**kappa
          END DO
        END DO
      END DO

      DO z=1,u_p_levs
        DO y= vdims%j_start,vdims%j_end      ! latitude loop
          DO x= udims%i_start,udims%i_end    ! longitude loop
            uv(x,y,z)     = u_p(x,y,z) * v_p(x,y,z)
            vtheta(x,y,z) = v_p(x,y,z) * theta_p(x,y,z)
            vT(x,y,z)     = v_p(x,y,z) * temperature(x,y,z)
            uw(x,y,z)     = u_p(x,y,z) * w_p(x,y,z)
          END DO
        END DO
      END DO

      CALL zonal_average(w_p, zonal_w,udims%i_start,udims%i_end,      &
                                      vdims%j_start,vdims%j_end,      &
                                      global_row_length,              &
                                      w_p_levs, gc_proc_row_group)
      CALL zonal_average(u_p, zonal_u,udims%i_start,udims%i_end,      &
                                      vdims%j_start,vdims%j_end,      & 
                                      global_row_length, &
                                      u_p_levs, gc_proc_row_group)

      CALL zonal_average(v_p, zonal_v,udims%i_start,udims%i_end,      &
                                      vdims%j_start,vdims%j_end,      & 
                                      global_row_length, &
                                      u_p_levs, gc_proc_row_group)

      CALL zonal_average(uv, zonal_uv,udims%i_start,udims%i_end,      &
                                      vdims%j_start,vdims%j_end,      & 
                                      global_row_length, &
                                      u_p_levs, gc_proc_row_group)

      CALL zonal_average(theta_p, zonal_theta,                        &
                                      udims%i_start,udims%i_end,      &
                                      vdims%j_start,vdims%j_end,      & 
                                      global_row_length, &
                                      u_p_levs, gc_proc_row_group)

      CALL zonal_average(temperature, zonal_t,                        &
                                      udims%i_start,udims%i_end,      &
                                      vdims%j_start,vdims%j_end,      & 
                                      global_row_length, &
                                      u_p_levs, gc_proc_row_group)

      CALL zonal_average(vtheta, zonal_vtheta,                        &
                                      udims%i_start,udims%i_end,      &
                                      vdims%j_start,vdims%j_end,      & 
                                      global_row_length, &
                                      u_p_levs, gc_proc_row_group)

      CALL zonal_average(uw, zonal_uw,udims%i_start,udims%i_end,      &
                                      vdims%j_start,vdims%j_end,      & 
                                      global_row_length, &
                                      u_p_levs, gc_proc_row_group)

      CALL zonal_average(vT, zonal_vT,udims%i_start,udims%i_end,      &
                                      vdims%j_start,vdims%j_end,      & 
                                      global_row_length, &
                                      u_p_levs, gc_proc_row_group)

! Now work out the zonal average of the prime quantities needed
! for EP flux. eg {u'v'}={uv}-{u}{v}, where {}=zonal average.
! Calculate the meridional heat (v'T') and momentum (u'v') fluxes
!
      DO z=1,u_p_levs
        DO y=vdims%j_start,vdims%j_end

! Calculate meridional momentum flux. STASH item 316. 
         zupvp(y,z) =zonal_uv(y,z) - zonal_u(y,z) * zonal_v(y,z)

! Calculate meridional heat flux. STASH item 315.
         zvptp(y,z)=zonal_vT(y,z) - zonal_v(y,z) * zonal_T(y,z)
         zvpthp(y,z)=zonal_vtheta(y,z) - zonal_v(y,z) * zonal_theta(y,z)
         zupwp(y,z) =zonal_uw(y,z) - zonal_u(y,z) * zonal_w(y,z)
        END DO
      END DO

!
!  Calculate dtheta/dz
!
! For lowest level 
      DO y=vdims%j_start,vdims%j_end
        dthdz(y,1) = (zonal_theta(y,2) -                           &
                      zonal_theta(y,1)) /  (zlogc(2) - zlogc(1))
      END DO
      
! For highest level 
      DO y=vdims%j_start,vdims%j_end       
        dthdz(y,u_p_levs) = (zonal_theta(y,u_p_levs) -             &
                             zonal_theta(y,u_p_levs-1)) /          &
                             (zlogc(u_p_levs) - zlogc(u_p_levs-1))
      END DO
      
      DO z=2,u_p_levs-1
        DO y=vdims%j_start,vdims%j_end
          dthdz(y,z) = (zonal_theta(y,z+1) - zonal_theta(y,z-1)) / &
                             (zlogc(z+1)-zlogc(z-1))
        END DO
      END DO
    END IF

    IF(qvstarbar_p) THEN
!
!  Calculate vstarbar. STASH item 310
! ----------------------------------------------------------------------
!  v* = v_bar - (d/dz [rho0] / rho0 * v'th'/dthdz) - d/dz [ v'th'/dthdz ]
!
      top = u_p_levs
! For bottom level       
      DO y=vdims%j_start,vdims%j_end       
        vstarbar(y,1) = zonal_v(y,1) +                         &
                     ((zvpthp(y,1) / dthdz(y,1)) / sclht) -    &
                     (                                         &
                       ((zvpthp(y,2) / dthdz(y,2)) -           &
                        (zvpthp(y,1) / dthdz(y,1)) )/          &
                       (zlogc(2) - zlogc(1))                   &
                     )
      END DO
      
! For top level 
      DO y=vdims%j_start,vdims%j_end       
        vstarbar(y,top) = zonal_v(y,top) +                     &
                  ((zvpthp(y,top) / dthdz(y,top)) / sclht) -   &
                  (                                            &
                      ( (zvpthp(y,top) / dthdz(y,top) )-       &
                       ( zvpthp(y,top-1) / dthdz(y,top-1))) /  &
                       (zlogc(top) -zlogc(top-1))              &
                  )
      END DO
                  
      DO z=2,top-1
        DO y=vdims%j_start,vdims%j_end
              vstarbar(y,z) = zonal_v(y,z) +                         &
                              ((zvpthp(y,z) / dthdz(y,z)) / sclht) - &
                              (                                      &
                              (  (zvpthp(y,z+1) / dthdz(y,z+1) )-    &
                               (  zvpthp(y,z-1) / dthdz(y,z-1)))/    &
                                (zlogc(z+1) - zlogc(z-1))            &
                              )
        END DO
      END DO
    END IF ! on STASHflag

    IF(qwstarbar_p .OR. qepy_p .OR. qepz_p .OR. qdivep_p ) THEN
!
! Obtain latitude needed in calculation of wstarbar and Eliassen-Palm flux.
!
        DO y=vdims%j_start,vdims%j_end
          phi(y) = ASIN(sin_v_latitude(1,y))
        END DO
    END IF

    IF(qwstarbar_p) THEN
!
!  Calculate wstarbar. STASH item 311
! ----------------------------------------------------------------------
!  w* = w_bar - 1/(a cos{phi}) * d/dphi [ cos{phi} * v'th'/dthdz ]
!
! DEPENDS ON: swap_bounds
      CALL Swap_bounds(zvpthp,1,n_rows,u_p_levs,                    &
                       0,offy,fld_type_v,.TRUE.)
! DEPENDS ON: swap_bounds
      CALL Swap_bounds(dthdz,1,n_rows,u_p_levs,                     &
                       0,offy,fld_type_v,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL Swap_bounds(phi,1,n_rows,1,                              &
                       0,offy,fld_type_v,.FALSE.)

      DO z=u_p_levs,1,-1
        start = vdims%j_start
        end   = vdims%j_end

        IF (at_extremity(PSouth)) THEN
          start = vdims%j_start+1
          wstarbar(vdims%j_start,z) = zonal_w(vdims%j_start,z) +    &
                 (1. / (earth_radius * COS(phi(vdims%j_start)))) *  &
                   ( -SIN(phi(vdims%j_start))                       &
          * zvpthp(vdims%j_start,z) / dthdz(vdims%j_start,z) +      &
                     COS(phi(vdims%j_start)) *                      &
                     (                                              &
         (zvpthp(vdims%j_start+1,z) / dthdz(vdims%j_start+1,z)) -   &
         (zvpthp(vdims%j_start,z) / dthdz(vdims%j_start,z))         &
                     )                                              &
             / (phi(vdims%j_start+1) - phi(vdims%j_start))          &
                   )
        END IF

        IF (at_extremity(PNorth)) THEN
          end   = vdims%j_end-1
          wstarbar(vdims%j_end,z) = zonal_w(vdims%j_end,z) +        &
                   (1. / (earth_radius * COS(phi(vdims%j_end)))) *  &
                   (                                                &
           -SIN(phi(vdims%j_end)) * zvpthp(vdims%j_end,z)           &
                                      / dthdz(vdims%j_end,z) +      &
                     COS(phi(vdims%j_end)) *                        &
                     (                                              &
             (zvpthp(vdims%j_end,z) / dthdz(vdims%j_end,z)) -       &
             (zvpthp(vdims%j_end-1,z) / dthdz(vdims%j_end-1,z))     &
                     )                                              &
            / (phi(vdims%j_end) - phi(vdims%j_end-1))               &
                   )
        END IF

!CDIR NOUNROLL
        DO y=start,end
          wstarbar(y,z) = zonal_w(y,z) +                            &
                   (1. / (earth_radius * COS(phi(y)))) *            &
                   (                                                &
                     -SIN(phi(y)) * zvpthp(y,z) / dthdz(y,z) +      &
                     COS(phi(y)) *                                  &
                     (                                              &
                       (zvpthp(y+1,z) / dthdz(y+1,z)) -             &
                       (zvpthp(y-1,z) / dthdz(y-1,z))               &
                     )                                              &
                      / (phi(y+1) - phi(y-1))                       &
                   )
        END DO
      END DO
    END IF

    IF(qepy_p .OR. qdivep_p) THEN
!
!  Calculate dubar/dz if required
!
! bottom level 
      DO y=vdims%j_start,vdims%j_end 
        dudz(y,1) = (zonal_u(y,2) - zonal_u(y,1)) /                 &
                           (zlogc(2) - zlogc(1)) 
      END DO
                         
! top level 
      DO y=vdims%j_start,vdims%j_end 
        dudz(y,u_p_levs) = (zonal_u(y,u_p_levs) -                   &
                            zonal_u(y,u_p_levs-1)) /                &
                            (zlogc(u_p_levs) - zlogc(u_p_levs-1))
      END DO                      

      DO z=2,u_p_levs-1
        DO  y=vdims%j_start,vdims%j_end 
          dudz(y,z) = (zonal_u(y,z+1) - zonal_u(y,z-1)) /           &
                              (zlogc(z+1) - zlogc(z-1))
        END DO
      END DO
    END IF

    IF(qepz_p .OR. qdivep_p) THEN
!
!  Calculate d(ubar*cos(phi))/dphi * (1./radius*cos(phi))
!
        DO z=1,u_p_levs
          DO y=vdims%j_start,vdims%j_end 
            zonal_u_halo(y,z) = zonal_u(y,z)
          END DO
        END DO

! DEPENDS ON: swap_bounds
        CALL Swap_bounds(zonal_u_halo,1,n_rows,u_p_levs,            &
                         0,offy,fld_type_v,.TRUE.)

        DO z=1,u_p_levs
          start = vdims%j_start
          end   = vdims%j_end 
          IF (at_extremity(PSouth)) THEN
            start = vdims%j_start+1
            ducosdphi(vdims%j_start,z) = (                               &
    (zonal_u_halo(vdims%j_start,z) - zonal_u_halo(vdims%j_start+1,z)) /  &
     (phi(vdims%j_start) - phi(vdims%j_start+1)) -                       &
       zonal_u_halo(vdims%j_start,z) * TAN(phi(vdims%j_start))           &
                             ) / earth_radius
          END IF

          IF (at_extremity(PNorth)) THEN
            end = vdims%j_end-1
            ducosdphi(vdims%j_end,z) = (                                 &
        (zonal_u_halo(vdims%j_end-1,z) - zonal_u_halo(vdims%j_end,z)) /  &
         (phi(vdims%j_end-1) - phi(vdims%j_end)) -                       &
          zonal_u_halo(vdims%j_end,z) * TAN(phi(vdims%j_end))            &
                                  ) / earth_radius
          END IF
!CDIR NOUNROLL
          DO y=start,end
            ducosdphi(y,z) = (                                              &
                              (zonal_u_halo(y-1,z) - zonal_u_halo(y+1,z)) / &
                               (phi(y-1) - phi(y+1)) -                      &
                               zonal_u_halo(y,z) * TAN(phi(y))              &
                             ) / earth_radius
          END DO
        END DO
    END IF
! ----------------------------------------------------------------------
! STASH items 312 313: Eliassen-Palm flux (phi & z components) 
!                      on pressure surfaces/B grid
! ----------------------------------------------------------------------
!
! Calculate Eliassen-Palm flux
! ----------------------------------------------------------------------
!  F_phi = rho0 * a cos{phi} * ( dubardz * v'th'/dthdz - u'v' )
!  F_z   = rho0 * a cos{phi} * ( (2 Omega sin{phi} - 1/(a cos{phi}) * 
!               d/dphi [ubar cos{phi}] ) * v'th'/dthdz - u'w' )
!                                                                       
    IF(qepy_p .OR. qdivep_p) THEN
        DO z=1,u_p_levs
          DO y=vdims%j_start,vdims%j_end 
             Fy(y,z) = rho0(z) * earth_radius * COS(phi(y)) *       &
                        (                                           &
                          dudz(y,z) * zvpthp(y,z) / dthdz(y,z)      &
                          - zupvp(y,z)                              &
                        )
          END DO
        END DO
      END IF

      IF(qepz_p .OR. qdivep_p) THEN
        DO z=1,u_p_levs
          DO y=vdims%j_start,vdims%j_end 
             Fz(y,z) = rho0(z) * earth_radius * COS(phi(y)) *       &
                        (                                           &
                          (2. * omega * SIN(phi(y)) -               &
                           ducosdphi(y,z)) *                        &
                           (zvpthp(y,z) / dthdz(y,z)) - zupwp(y,z)  &
                        )
          END DO
        END DO
    END IF

    IF(qdivep_p) THEN
!
! Calculate Div F. STASH item 314
! ----------------------------------------------------------------------
! 1/(a cos{phi}) * d/dphi [ cos{phi} F_phi ] + d/dz [F_z] 
!
      DO z=1,u_p_levs
        DO y=vdims%j_start,vdims%j_end 
           Fy_halo(y,z) = Fy(y,z)
        END DO
      END DO

! DEPENDS ON: swap_bounds
      CALL Swap_bounds(Fy_halo,1,n_rows,u_p_levs,                   &
                       0,offy,fld_type_v,.TRUE.)

! bottom level 
      DO y=vdims%j_start,vdims%j_end 
        divFz(y,1) = (Fz(y,2) - Fz(y,1)) /                          &
                          (zlogc(2) - zlogc(1))
      END DO
                           
! top level 
      DO y=vdims%j_start,vdims%j_end 
        divFz(y,u_p_levs) = (Fz(y,u_p_levs) -                       &
                                    Fz(y,u_p_levs-1)) /             &
                              (zlogc(u_p_levs) - zlogc(u_p_levs-1))
      END DO                     
      
      DO z=2,u_p_levs-1
!CDIR NOUNROLL
        DO y=vdims%j_start,vdims%j_end 
            divFz(y,z) = (Fz(y,z+1) - Fz(y,z-1)) /                  &
                          (zlogc(z+1) - zlogc(z-1))
        END DO
      END DO

      DO z=1,u_p_levs
        start = vdims%j_start
        end   = vdims%j_end
        
        IF (at_extremity(PSouth)) THEN
          start = vdims%j_start+1
          divFy(vdims%j_start,z) = (                                &
      (Fy_halo(vdims%j_start,z) - Fy_halo(vdims%j_start+1,z)) /     &
       (phi(vdims%j_start) - phi(vdims%j_start+1)) -                &
              Fy_halo(vdims%j_start,z) * TAN(phi(vdims%j_start))    &
                       ) / earth_radius
        END IF

        IF (at_extremity(PNorth)) THEN
          end = vdims%j_end-1
          divFy(vdims%j_end,z) = (                                  &
             (Fy_halo(vdims%j_end-1,z) - Fy_halo(vdims%j_end,z)) /  &
                (phi(vdims%j_end-1) - phi(vdims%j_end)) -           &
                 Fy_halo(vdims%j_end,z) * TAN(phi(vdims%j_end))     &
                       ) / earth_radius
        END IF
!CDIR NOUNROLL
        DO y=start,end
            divFy(y,z) = (                                      &
                          (Fy_halo(y-1,z) - Fy_halo(y+1,z)) /   &
                           (phi(y-1) - phi(y+1)) -              &
                           Fy_halo(y,z) * TAN(phi(y))           &
                         ) / earth_radius
        END DO
!
! Calculate scaled Div F  
!----------------------------------------------------------------------
! Div.F / (ro0 * a cos{phi}
! Add code to cover for issues when cos(phi(y))=0 as at ENDGame poles. 
! set poles to zero to avoid divde by tiny values...

        start = vdims%j_start
        end   = vdims%j_end

        IF ( l_vatpoles ) THEN
        IF (at_extremity(PSouth)) THEN
          start = vdims%j_start+1
          divF(vdims%j_start,z) = 0.0          
        END IF

        IF (at_extremity(PNorth)) THEN
          end = vdims%j_end-1 
          divF(vdims%j_end,z) = 0.0
        END IF            
        END IF  ! vatpoles
!CDIR NOUNROLL
        DO y=start,end
        
          divF(y,z) = (divFy(y,z) + divFz(y,z)) /                &
                       (rho0(z) * earth_radius * COS(phi(y)))
        END DO
      END DO
    END IF

    IF (lhook) CALL dr_hook('CALC_DIV_EP_FLUX',zhook_out,zhook_handle)
    RETURN

    END SUBROUTINE Calc_Div_EP_Flux

END MODULE calc_div_ep_flux_mod
