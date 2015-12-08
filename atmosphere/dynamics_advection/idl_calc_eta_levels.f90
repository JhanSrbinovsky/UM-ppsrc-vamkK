! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine  IDL_Calc_eta_levels

      Subroutine IDL_Calc_eta_levels(                                   &
     &                      row_length, rows, model_levels              &
     &,                     me, n_proc                                  &
     &,                     first_constant_rho_level                    &
     &,                     eta_theta_levels, eta_rho_levels            &
!  Grid information
     &,                     grid_number, height_domain                  &
     &,                     first_theta_height, thin_theta_height       &
     &,                     big_layers, transit_layers, mod_layers      &
     &,                     big_factor, mag                             &
     &,                     L_code_test)

! Purpose:
!          Sets up vertical grid (eta_levels) for 3d LAM
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

      Integer                                                           &
     &  me                                                              &
                   ! My processor number
     &, n_proc     ! Total number of processors

      Logical                                                           &
     & L_code_test  ! User switch

      Integer                                                           &
     &  grid_number                                                     &
     &, big_layers, transit_layers, mod_layers

      Real                                                              &
     &  height_domain                                                   &
     &, first_theta_height                                              &
     &, thin_theta_height                                               &
     &, big_factor


      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, model_levels                                                    &
                         ! number of model levels
     &, first_constant_rho_level

      Real                                                              &
           ! vertical co-ordinate information
     &  eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

! local variables

      Real                                                              &
     &  delta_z                                                         &
     &, delta_eta                                                       &
     &, mag_value1, mag_value2, mag_value3                              &
     &, mag_factor1, mag_factor2, mag_factor3                           &
     &, top_mag_value                                                   &
     &, mag, mag3

      Integer                                                           &
     & k, fcrl, level

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


       Parameter( mag_value1 = 1.04 )  ! value for 38 levels
       Parameter( mag_value2 = 1.1001 )  ! value for 38 levels
       Parameter( mag_value3 = 1.1)  ! value for 38 levels
       Parameter( top_mag_value = 1.5 )  ! value for 38 levels
       Parameter( mag_factor1 = 0.01 )  ! value for 38 levels
       Parameter( mag_factor2 = 0.02 )  ! value for 38 levels
       Parameter( mag_factor3 = 0.02)  ! value for 38 levels

! Description: COMDECK containing vertical grid types
!  for use in idealised problems
!
      INTEGER, PARAMETER :: vert_regular=1
      INTEGER, PARAMETER :: vert_quadratic_theta=21
      INTEGER, PARAMETER :: vert_bi_quadratic=22
      INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
      INTEGER, PARAMETER :: vert_schar=3
      INTEGER, PARAMETER :: vert_dwd=4
      INTEGER, PARAMETER :: vert_stretch_plus_regular=5
      INTEGER, PARAMETER :: vert_quad_stretch_thin=6
      INTEGER, PARAMETER :: vert_regular_thin=7
      INTEGER, PARAMETER :: vert_geometric_theta=8
      INTEGER, PARAMETER :: vert_dump=10
      INTEGER, PARAMETER :: vert_idl_um_grid=11

! No External routines

! ----------------------------------------------------------------------
! Section 1.  Initialise Data fields.
!             Set reference vertical grid
!             Set eta_theta_levels
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_CALC_ETA_LEVELS',zhook_in,zhook_handle)

!   Change first_constant_rho_level to ensure levels are flat in
!   any stretching region
      if ( big_layers  >   0) then
        fcrl = model_levels - big_layers - transit_layers
        if(fcrl  <   first_constant_rho_level)then
          first_constant_rho_level = fcrl
        endif  !fcrl  <   first_constant_rho_level
      endif   !big_layers  >   0
!   fcrl shorthand for first_constant_rho_level
      fcrl = first_constant_rho_level

      eta_theta_levels(0) = 0.0

      if(grid_number  ==  vert_regular )then

!!!!!!!!!!      Regular grid start   !!!!!!!!!!!!!!!!!!!!!!!!!

        delta_eta = 1.0/Real(model_levels)
        eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to regular grid
        Do k=1, model_levels
          eta_theta_levels(k) = Real(k) * delta_eta
        end do
        if(me  ==  0)then
          Write(Unit=6,Fmt=*) ' Regular vertical grid selected'
          Write(Unit=6,Fmt=*) '   eta_theta_levels(top) = ',            &
     &                         eta_theta_levels(model_levels)
          Write(Unit=6,Fmt=*) '   height_domain =',height_domain,       &
     &                        ' metres'
          Write(Unit=6,Fmt=*) '   each layer thickness =',              &
     &                         height_domain*delta_eta
        Endif    !(me  ==  0)
! Set eta_rho_levels
        do k=1,model_levels
          eta_rho_levels(k) =                                           &
     &    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
        end do
!!!!!!!!!!      Regular grid end  !!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(grid_number  ==  vert_quadratic_theta )then

!!!!!!!!!!     Quadratic grid start  !!!!!!!!!!!!!!!!!!!!!!!!!

        delta_eta = 1.0/Real(model_levels)
        eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to quadratic grid
        Do k=1, model_levels
          eta_theta_levels(k) = (Real(k) * delta_eta)                   &
     &                        * (Real(k) * delta_eta)
        end do
        if(me  ==  0)then
          print*,'***** Quadratic grid for theta levels *****'
          print*,'eta_theta_levels(top) =',                             &
     &    eta_theta_levels(model_levels)
          print*,'height_domain =',height_domain,' metres'
        Endif      !(me  ==  0)
! Set eta_rho_levels
        do k=1,model_levels
          eta_rho_levels(k) =                                           &
     &    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
        end do
!!!!!!!!!!      Quadratic grid for theta levels end  !!!!!!!!

      elseif(grid_number  ==  vert_bi_quadratic)then

!!!!!!!!!!     Quadratic u & theta grids start  !!!!!!!!!!!!!!!!!
        delta_eta = 1.0/Real(model_levels)
        eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to quadratic grid
        Do k=1, model_levels
          eta_theta_levels(k) = (Real(k) * delta_eta)                   &
     &                        * (Real(k) * delta_eta)
        end do
        if(me  ==  0)then
          print*,'***** Quadratic grid for u & theta levels *****'
       print*,'eta_theta_levels(top) =',eta_theta_levels(model_levels)
          print*,'height_domain =',height_domain,' metres'
        Endif    !(me  ==  0)
!!!!!!!!!!      Quadratic grids for u & theta levels end  !!!!!!!!

      elseif(grid_number  ==  vert_quadratic_uvtheta)then

!!!!!!!!!!    Common Quadratic u & theta grid start  !!!!!!!!!!

        delta_eta = 1.0/Real(2*model_levels)
        eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to quadratic grid
        Do k=1, model_levels
          eta_rho_levels(k) = (Real(2*k-1) * delta_eta)                 &
     &                      * (Real(2*k-1) * delta_eta)
          eta_theta_levels(k) = (Real(2*k) * delta_eta)                 &
     &                      * (Real(2*k) * delta_eta)
        end do
        if(me  ==  0)then
          print*,'** Common Quadratic grid for u & theta levels **'
       print*,'eta_theta_levels(top) =',eta_theta_levels(model_levels)
          print*,'height_domain =',height_domain,' metres'
        Endif    !(me  ==  0)

!!!!!!!!  Common Quadratic grid for u & theta levels end  !!!!!!!!

      elseif(grid_number  ==  vert_schar)then

!!!!!!!!!!     Schar grid start  !!!!!!!!!!!!!!!!!!!!!!!!!

        Do k=0, model_levels
! These values for eta_theta_levels equate approximately to those used
! by Schar.
          eta_theta_levels(k) = Real(k)**1.2/                           &
     &                        real(model_levels)**1.2
        end do
        if(me  ==  0)then
          print*,'***** Schar grid selected *****'
       print*,'eta_theta_levels(top) =',eta_theta_levels(model_levels)
          print*,'height_domain =',height_domain,' metres'
        Endif    !(me  ==  0)

! Set eta_rho_levels
        do k=1,model_levels
          eta_rho_levels(k) =                                           &
     &    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
        end do
!!!!!!!!!!     Schar grid end  !!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(grid_number  ==  vert_dwd )then

!!!!!!!!!!      DWD stretched start   !!!!!!!!!!!!!!!!!!!!!!!!!

! These values for eta_theta_levels similar to DWD
! Set eta_theta_levels to actual height values
        eta_theta_levels(0) = 0.0
        eta_theta_levels(1) = 40.0
        eta_theta_levels(2) = 100.0
        delta_z = 100.0
        Do k = 3, model_levels
! These values for eta_theta_levels equate to stretched  grid
          eta_theta_levels(k) = eta_theta_levels(k-1) + delta_z
          delta_z = delta_z + 40.0
        end do
        height_domain = eta_theta_levels( model_levels)
! Now normalise eta_theta_levels
        Do k=1,model_levels
          eta_theta_levels(k) =                                         &
     &    eta_theta_levels(k)/eta_theta_levels( model_levels)
        end do
        if(me  ==  0)then
          print*,'***** DWD stretched grid selected *****'
       print*,'eta_theta_levels(top) =',eta_theta_levels(model_levels)
          print*,'height_domain =',height_domain,' metres'
        Endif    !(me  ==  0)
! Set eta_rho_levels
        do k=1,model_levels
          eta_rho_levels(k) =                                           &
     &    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
        end do
!!!!!!!!!!      DWD stretched end  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif(grid_number  ==  vert_stretch_plus_regular)then

!!!!!!!!!!     regular grid followed by stretched (big_layers)!!!!!

! mod_layers is the number of regular layers
        mod_layers=model_levels-transit_layers-big_layers
! mag is the magnification factor applied over transit_layers
        mag=big_factor**(1.0/real(max(1,transit_layers)))
! interval chosen to scale to 1 over domain
!  NB   mag * (1.0 - big_factor)/(1.0 - mag) +
!    if first transit layer magnified
        delta_eta = 1.0/ (   real(mod_layers) +                         &
     &             (1.0 - big_factor)/(1.0 - mag) +                     &
     &             big_factor * real(big_layers) )
! Make height domain for problem the full height domain
! For this grid the INPUT height domain is for regular levels
        height_domain = height_domain/(real(mod_layers)*delta_eta)
        eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to regular grid
        Do k=1, mod_layers
          eta_theta_levels(k) = Real(k) * delta_eta
        end do
        Do k = mod_layers + 1 ,mod_layers + transit_layers
! These values for eta_theta_levels equate to stretched  grid
          eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
          delta_eta = delta_eta * mag
        end do
! These values for eta_theta_levels equate to stretched  grid
        Do k = mod_layers + transit_layers + 1 , model_levels
          eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
        end do
        if(me  ==  0)then
          print*,'** regular grid followed by stretched grid selected'
         print*,'     regular grid =', height_domain*delta_eta,' metres'&
     &    ,' over ',mod_layers,' layers'
          print*,'     stretching over ',transit_layers,' layers',      &
     &    ' with magnification factor =',mag,' times'
          print*,'     thick ',big_factor,'* regular over ',            &
     &    big_layers,' layers'
          print*,'eta_theta_levels(top) =',                             &
     &    eta_theta_levels(model_levels)
          print*,'height_domain =',height_domain,' metres'
        Endif    !(me  ==  0)
! Set eta_rho_levels
        do k = 1, model_levels
          eta_rho_levels(k) =                                           &
     &    0.5*(eta_theta_levels(k-1) + eta_theta_levels(k) )
        end do

!!!!!!! end of regular grid followed by stretched (big_layers)!!!

      elseif(grid_number  ==  vert_quad_stretch_thin)then
        if(me  ==  0)then
          print*,'********grid option 6  ****'
          print*,'*quadratic grid from level 2 *'
          print*,'*        followed by expanding layers    *  '
          print*,'*  thin layer near surface    *  '
        Endif    !(me  ==  0)

! For this grid need to work with model_levels - 1 to generate
!  quadratic part
        delta_eta = 1.0/Real(model_levels - 1)
        eta_theta_levels(0) = 0.0
! These values for eta_theta_levels equate to quadratic grid
        Do k = 1, model_levels - 1
          eta_theta_levels(k) = (Real(k) * delta_eta)                   &
     &                      * (Real(k) * delta_eta)
        end do
! Use eta_rho_levels to hold ratio of succesive eta differences
        level = 0
        Do k = 1, model_levels - 2
         eta_rho_levels(k)=(eta_theta_levels(k+1)- eta_theta_levels(k)) &
     &                    /(eta_theta_levels(k)- eta_theta_levels(k-1))
          if (level == 0) then
            if (eta_rho_levels(k) <   mag_value1) level = k
          end if
        end do

        if (level  /=  0) then
          delta_eta = eta_theta_levels(level ) -                        &
     &                eta_theta_levels(level - 1)
          mag = mag_value1
          mag3 = mag_factor3
          Do k = level + 1 , model_levels - 1
            delta_eta = mag * delta_eta
            eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
            mag = mag + mag_factor1
            if(mag >  mag_value3)then
              mag = mag - mag_factor1 + mag3
              mag3 = 2.0 * mag3
              if(mag >  2.0)mag = 1.3
            else if(mag >  mag_value2)then
              mag = mag - mag_factor1 + mag_factor2
            end if
          end do
        end if  !  level  /=  0

! Re-normalise eta_theta_levels and set eta_rho_levels
! Store current top of quadratic grid in eta_theta_levels(model_levels)
        eta_theta_levels(model_levels) =                                &
     &                   eta_theta_levels(model_levels - 1)
! Work from top to allow level 1 values to be inserted
        do k= model_levels - 1, 2, -1
          eta_theta_levels(k) = eta_theta_levels(k-1) /                 &
     &                        eta_theta_levels(model_levels)
        end do
!!!!!  reset height domain to give r_theta = first_theta_height
!!!!!                                    at level 1
        height_domain = first_theta_height / eta_theta_levels(2)
        if(me  ==  0)then
          print*,'height_domain =',height_domain,' metres'
        Endif    !(me  ==  0)
! Re-normalise eta_theta_levels(model_levels)
        eta_theta_levels(model_levels) = 1.0
! Now add in thin bottom layer
        eta_theta_levels(1) = thin_theta_height / height_domain
        eta_rho_levels(1) = 0.5 * eta_theta_levels(1)
        do k= 2, model_levels
          eta_rho_levels(k) =                                           &
     &        0.5 * (eta_theta_levels(k-1) + eta_theta_levels(k) )
        end do

!!!!!! end of quadratic grid followed by expanding layers!!!

      elseif(grid_number  >   vert_dump )then
        if(me  ==  0)then
          print*,'********grid option not supported****'
        Endif    !(me  ==  0)

      endif  ! grid_number options

      IF (lhook) CALL dr_hook('IDL_CALC_ETA_LEVELS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_Calc_eta_levels

!  End Subroutine IDL_Calc_eta_levels

