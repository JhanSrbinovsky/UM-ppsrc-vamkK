! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Adds a temperature anomaly (bubble) to the temperature field

      SUBROUTINE IDL_initialise_bubble(                                 &
                            row_length, rows, halo_i, halo_j            &
      ,                     model_levels, wet_model_levels              &
      ,                     delta_lambda, delta_phi, base_phi           &
      ,                     r_theta_levels                              &
      ,                     p_zero, me, l_datastart                     &
      ,                     global_row_length, global_rows              &
      ,                     idl_max_num_bubbles, idl_bubble_option      &
      ,                     idl_bubble_max, idl_bubble_height           &
      ,                     idl_bubble_xoffset, idl_bubble_yoffset      &
      ,                     idl_bubble_width, idl_bubble_depth          &
      ,                     L_idl_bubble_saturate, L_trivial_trigs      &
      ,                     exner_theta_levels, theta, q                &
        )

      USE earth_constants_mod, ONLY: earth_radius
      USE atmos_constants_mod, ONLY: recip_kappa
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Description:
!   Adds a 3D temperature anomaly (bubble) to the potential temperature
!   field (warm or cold). There are a number of options for different
!   shaped bubbles. There is also the option to saturate the bubble.
!
! Method:
!   Sets up x,y,z functions containing the distance from the centre
!   of the bubble, then applies a function to generate the bubble
!   amplitude at all points, and adds this to the model field.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
! Subroutine arguments

      ! Array Dimensions
      Integer, Intent (In) :: row_length   ! no. points on a row
      Integer, Intent (In) :: rows         ! no. rows
      Integer, Intent (In) :: halo_i       ! Size of halo x-direction
      Integer, Intent (In) :: halo_j       ! Size of halo y-direction
      Integer, Intent (In) :: model_levels ! number of model levels
      Integer, Intent (In) :: wet_model_levels ! no. of wet levels

      ! Horizontal co-ordinate information
      Real, Intent (In) :: delta_lambda    ! grid spacing (x-direction)
      Real, Intent (In) :: delta_phi       ! grid spacing (y-direction)
      Real, Intent (In) :: base_phi        ! first latitude

      Real, Intent (In) :: r_theta_levels(1-halo_i:row_length+halo_i,   &
     &                           1-halo_j:rows+halo_j,0:model_levels)

       ! Physical Constants
       Real, Intent (In) :: p_zero         ! reference pressure

      ! Multi-processor array variables
      Integer, Intent (In) :: me                ! My processor number
      Integer, Intent (In) :: global_row_length ! global number of pts
      Integer, Intent (In) :: global_rows       ! global number of rows
      Integer, Intent (In) :: l_datastart(2)    ! Start row for pe


      ! Idealised options

      Integer, Intent (In) :: idl_max_num_bubbles
                            ! Maximum number of bubbles allowed
      Integer, Intent (In) :: idl_bubble_option(idl_max_num_bubbles)
                            ! Bubble option number

      Real, Intent (In)  :: idl_bubble_max(idl_max_num_bubbles)
                          ! Bubble maximum amplitude (K)
      Real, Intent (In)  :: idl_bubble_height(idl_max_num_bubbles)
                          ! Height of bubble centre (m)
      Real, Intent (In)  :: idl_bubble_xoffset(idl_max_num_bubbles)
               ! Bubble x-offset (normalised units: 0.5 = domain centre)
      Real, Intent (In)  :: idl_bubble_yoffset(idl_max_num_bubbles)
               ! Bubble y-offset (normalised units: 0.5 = domain centre)
      Real, Intent (In)  :: idl_bubble_width(idl_max_num_bubbles)
                          ! Horizontal width of bubble (m)
      Real, Intent (In)  :: idl_bubble_depth(idl_max_num_bubbles)
                          ! Horizontal depth of bubble (m)

      Logical, Intent (In) :: L_idl_bubble_saturate(idl_max_num_bubbles)
                            ! True if bubble saturated
      Logical, Intent (In) :: L_trivial_trigs


      ! Data arrays
      Real, Intent (In) ::                                              &
     &  exner_theta_levels(1-halo_i:row_length+halo_i                   &
     &,                    1-halo_j:rows+halo_j, model_levels)

      Real, Intent (InOut) ::                                           &
     &  theta(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &      model_levels)    ! Potential Temperature (K)

      Real, Intent (InOut) ::                                           &
     &  q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &    wet_model_levels)  ! Humidity (kg/kg)

! Local Variables

      Integer i, j, k   ! Loop indices
      Integer b         ! bubble number
      Integer length    ! length of a row+halos

      Real thbub(1-halo_i:row_length+halo_i) ! bubble potential temp.
      Real tbub(1-halo_i:row_length+halo_i)  ! bubble temperature
      Real qbub(1-halo_i:row_length+halo_i)  ! bubble specific humidity
      Real pbub(1-halo_i:row_length+halo_i)  ! bubble pressure
      Real tx,ty,gj,dist,hght
      Real bubblefn(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,   &
     &          model_levels)
      Real dx(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j)
      Real dy

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!- End of header
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_INITIALISE_BUBBLE',zhook_in,zhook_handle)
      If (me == 0) Then
        Write (6,*) ' '
        Write (6,*) ' BUBBLE PERTURBATIONS '
      End If

      ! Calculate x-direction gridlength (in metres)
      If (L_trivial_trigs) Then
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            dx(i,j) = delta_lambda*Earth_radius
          End Do
        End Do
      Else ! not trivial trigs
        Do j = 1-halo_j, rows+halo_j
          gj = l_datastart(2) + j - 1
          Do i = 1-halo_i, row_length+halo_i
            dx(i,j) = COS(base_phi + (gj-1) * delta_phi)                &
     &                    *delta_lambda*Earth_radius
          End Do
        End Do
      End If ! on trivial_trigs

      ! Calculate y-direction gridlength (in metres)
      dy = delta_phi*Earth_radius

      !-----------------------------------------------------------------
      !
      !                   Loop over bubbles
      !
      !-----------------------------------------------------------------
      Do b = 1,idl_max_num_bubbles

        If ( idl_bubble_option(b) > 0 ) Then

          If (me == 0) Then
            Write (6,*) ' '
            Write (6,Fmt='(A16,I2)') '  Bubble number ',b
          End If

          ! Initialise bubble function to 0
          bubblefn(:,:,:) = 0.0

          !-------------------------------------------------
          !
          ! Calculate a "Gaussian" bubble perturbation
          !  EXP(-(dist/scale)**2)
          !
          !-------------------------------------------------
          If (idl_bubble_option(b) == 1) Then

            If (me == 0) Then
             Write (6,*) ' Option 1: Gaussian bubble'
             Write (6,Fmt='(A31,F8.2,A2)')                              &
     &         '   Initiating bubble with centre at height ',           &
     &         idl_bubble_height(b),' m'
              Write (6,Fmt='(A25,F7.2,A22,F7.2,A2)')                    &
     &         '   with horizontal scale ',idl_bubble_width(b),         &
     &         ' m and vertical scale ',idl_bubble_depth(b),' m'
              Write (6,Fmt='(A44,F5.2,A2)')                             &
     &         '   Max potential temperature perturbation = ',          &
     &                             idl_bubble_max(b),' K'
            End If

            ! Calculate bubble 3D Gaussian function
            tx = l_datastart(1)-1-global_row_length                     &
     &                            *idl_bubble_xoffset(b)
            ty = l_datastart(2)-1-global_rows*idl_bubble_yoffset(b)
            Do k = 1, model_levels
              Do j = 1-halo_j, rows+halo_j
                Do i = 1-halo_i, row_length+halo_i
                  dist =                                                &
     &               ((tx + i - 0.5)*dx(i,j) / idl_bubble_width(b))**2  &
     &              +((ty + j - 0.5)*dy / idl_bubble_width(b))**2       &
     &              +((r_theta_levels(i,j,k)-Earth_radius               &
     &              - idl_bubble_height(b)) / idl_bubble_depth(b))**2
                  bubblefn(i,j,k) = EXP(-1.*dist)
                End Do
              End Do
            End Do

          End If  ! on (idl_bubble_option == 1)


          !---------------------------------------------
          !
          ! Calculate a "block" bubble perturbation
          !
          !---------------------------------------------
          If (idl_bubble_option(b) == 2) Then

            If (me == 0) Then
              Write (6,*) ' Option 3: Block bubble'
              Write (6,Fmt='(A31,F8.2,A2)')                             &
     &         '   Initiating bubble with centre at height ',           &
     &         idl_bubble_height(b),' m'
              Write (6,Fmt='(A25,F7.2,A22,F7.2,A2)')                    &
     &         '   with horizontal width ',idl_bubble_width(b),         &
     &         ' m and vertical depth ',idl_bubble_depth(b),' m'
              Write (6,Fmt='(A44,F5.2,A2)')                             &
     &         '   Potential temperature perturbation = ',              &
     &                             idl_bubble_max(b),' K'
            End If

            ! Calculate bubble 3D Gaussian function
            tx = l_datastart(1)-1-global_row_length                     &
     &                             *idl_bubble_xoffset(b)
            ty = l_datastart(2)-1-global_rows*idl_bubble_yoffset(b)
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                dist = SQRT(                                            &
     &               ((tx + i - 0.5)*dx(i,j) / idl_bubble_width(b))**2  &
     &              +((ty + j - 0.5)*dy / idl_bubble_width(b))**2 )
                hght = ABS(r_theta_levels(i,j,k)-Earth_radius           &
     &                 - idl_bubble_height(b))
                Do k = 1, model_levels
                  If (dist <= 0.5 .AND. hght <= idl_bubble_depth(b)/2.) &
     &             Then
                     bubblefn(i,j,k) = 1.0
                   End If
                End Do
              End Do
            End Do

          End If ! on (idl_bubble_option == 2)


          !---------------------------------------------------------
          !
          ! Add bubble temperature perturbation to the initial state
          !
          !---------------------------------------------------------

          Do k = 1, model_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                theta(i,j,k) = theta(i,j,k) +                           &
     &                         idl_bubble_max(b)*bubblefn(i,j,k)
              End Do
            End Do
          End Do


          !---------------------------------------------------------
          !
          ! Set bubble to saturation if required
          !
          !---------------------------------------------------------
          If (L_idl_bubble_saturate(b)) Then

            length = row_length + halo_i + halo_i
            Do k = 1, model_levels
              Do j = 1-halo_j, rows+halo_j
                Do i = 1-halo_i, row_length+halo_i
                  thbub(i) = idl_bubble_max(b)*bubblefn(i,j,k)
                  pbub(i) = p_zero*exner_theta_levels(i,j,k)            &
     &                             **recip_Kappa
                  tbub(i) = theta(i,j,k) * exner_theta_levels(i,j,k)
                End Do
! DEPENDS ON: qsat
                Call qsat( qbub(1-halo_i), tbub(1-halo_i),              &
     &                     pbub(1-halo_i), length)
                ! Set bubble to saturation
                Do i = 1-halo_i, row_length+halo_i
                  If (thbub(i)  >   0.1) q(i,j,k) = qbub(i)
                End Do
              End Do
            End Do

          End If ! on L_bubble_saturate


        End If ! on idl_bubble_option(b) > 0

      End Do ! on b, loop over bubbles


      IF (lhook) CALL dr_hook('IDL_INITIALISE_BUBBLE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_initialise_bubble

      !  End subroutine IDL_initialise_bubble
