! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Adds a random number field to t,q for idealised atmosphere runs
!
! Description:
!   Adds a random number field to the temperature and humidity fields
!   for idealised atmosphere runs dependent on options specified in the
!   idealised namelist.
!
! Method:
!   Uses RANDOM_NUMBER routine to generate a random number field from
!   0 to 1 for the global field on pe0, then scatters it across all
!   processors. Options to apply to just t or q or both from the lowest
!   model level up to a specified height.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection

      SUBROUTINE IDL_random_perturb(                                    &
     &                      row_length, rows, halo_i, halo_j            &
     &,                     model_levels, wet_model_levels              &
     &,                     eta_theta_levels, height_domain             &
     &,                     me, n_proc, all_proc_group                  &
     &,                     global_row_length, global_rows              &
     &,                     g_rows, g_row_length                        &
     &,                     L_perturb_t, perturb_magnitude_t            &
     &,                     L_perturb_q, perturb_magnitude_q            &
     &,                     L_perturb_correlate_tq                      &
     &,                     L_perturb_correlate_vert                    &
     &,                     L_perturb_correlate_time                    &
     &,                     perturb_type, perturb_height                &
     &,                     theta, q                                    &
     &                      )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      USE UM_Parvars, ONLY : g_datastart
IMPLICIT NONE

      ! Array Dimensions
      Integer, Intent (In) :: row_length       ! no. of points on a row
      Integer, Intent (In) :: rows             ! no. of rows
      Integer, Intent (In) :: halo_i           ! Size of halo in x
      Integer, Intent (In) :: halo_j           ! Size of halo in y
      Integer, Intent (In) :: model_levels     ! no. of model levels
      Integer, Intent (In) :: wet_model_levels ! no. of wet levels

      ! Grid information
      Real, Intent (In) :: eta_theta_levels(0:model_levels)
      Real, Intent (In) :: height_domain

      ! Multi-processor array variables
      Integer, Intent (In) :: me                ! My processor number
      Integer, Intent (In) :: n_proc            ! Total number of pes
      Integer, Intent (In) :: all_proc_group    ! Group identifier
      Integer, Intent (In) :: global_row_length ! global number of pts
      Integer, Intent (In) :: global_rows       ! global number of rows
      Integer, Intent (In) :: g_rows(0:n_proc-1)
      Integer, Intent (In) :: g_row_length(0:n_proc-1)


      ! Idealised perturbation settings
      Logical, Intent (In) :: L_perturb_t ! Switch for theta perturbs
      Logical, Intent (In) :: L_perturb_q ! Switch for moisture perturbs
      Logical, Intent (In) :: L_perturb_correlate_tq
      Logical, Intent (In) :: L_perturb_correlate_vert
      Logical, Intent (In) :: L_perturb_correlate_time

      Real, Intent (In) ::                                              &
     &  perturb_magnitude_t  ! Magnitude of theta perturbations (K)

      Real, Intent (In) ::                                              &
     &  perturb_magnitude_q  ! Magnitude of theta perturbations (kg/kg)

      Integer, Intent (In) ::                                           &
     &  perturb_type         ! Form of perturbation (eg. random)

      Real, Intent (InOut) ::                                           &
     &  perturb_height(2)    ! Min and max heights for perturbations


      ! Data arrays
      Real, Intent (InOut) ::                                           &
     &  q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &    wet_model_levels)
      Real, Intent (InOut) ::                                           &
     &  theta(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &      model_levels)

! Local Variables
      Integer                                                           &
     &  i, j, k, d, gi, gj                                              &
                            ! Loop indices
     &, levbot,levtop                                                   &
                            ! Bottom and top levels for perturbations
     &, ransize                                                         &
                            ! Random number seed size
     &, icode, info

      Integer, Dimension(:), Allocatable :: iseed

      Real                                                              &
     &  z_theta_levels(model_levels)                                    &
                                     !height of eta_theta_levels (m)
     &, maxran, minran      ! Max and min perturbations for diagnostic

      ! Random number generator variables
      Real                                                              &
     &  ran_pe0(1-halo_i:global_row_length+halo_i,                      &
     &          1-halo_j:global_rows+halo_j)                            &
     &, rbuf(1-halo_i:global_row_length+halo_i,                         &
     &       1-halo_j:global_rows+halo_j)                               &
     &, ran(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



!- End of header
! ----------------------------------------------------------------------

        IF (lhook) CALL dr_hook('IDL_RANDOM_PERTURB',zhook_in,zhook_handle)

        ! Pass in array, local dimensions, global dimensions, seed
        If (me  ==  0) Then
          Write(6,*) '----------------------'
          Write(6,*) ' Random Perturbations '
          Write(6,*) '----------------------'

          !-------------------------------------------------------------
          ! Create grid array of random numbers on pe0
          !-------------------------------------------------------------

          !Calling random number generator
          ! ransize is an output from routine RANDOM_SEED and contains
          ! the number of default type integers (N) that are needed to
          ! hold the value of the seed, which is an 8-byte variable.

          Call RANDOM_SEED(size=ransize)

          ! Initialise random number generator
          Allocate (iseed(ransize))

          ! Set up a random number seed for RANDOM_NUMBER call
          Call RANDOM_SEED(GET=iseed(1:ransize))

          ! Generate random numbers for global grid on processor 0
          Call RANDOM_NUMBER(ran_pe0)

          Deallocate (iseed)

          !-------------------------------------------------------------
          ! Apply any optional algorithms to random number field
          !-------------------------------------------------------------

          ! Gridscale random noise
          If (perturb_type  ==  1) Then

            ! Do nothing to random number field

! Option not yet coded
!          ! Define a 2D Gaussian BUMP and add it onto existing field
!          Else If (perturb_type  ==  2) Then
!
!            pert(ipwl,jpwl)=
!     &      pert(ipwl,jpwl)+((1.0/(GAUSSTDV*(2.0*Pi)**0.5)) *
!     &      RANUM(GAUSCENX,GAUSCENY) *
!     &      EXP(-(((i-GAUSCENX)**2.0)+((j-GAUSCENY)**2.0))
!     &      /(2.0*GAUSSTDV**2.0)))

          Else

            Write (6,*) '  Perturbation type not recognised'
            Write (6,*) '  perturb_type=',perturb_type
            Write (6,*) '  No perturbations applied.'

          End If

        End If   !(me  ==  0)

        !-------------------------------------------------------------
        ! Scatter random number array to all processors
        !-------------------------------------------------------------

! Could replace this section with SCATTER_FIELD
!        CALL SCATTER_FIELD(ran,ran_pe0,
!     &        row_length+2*halo_i,rows+2*halo_j,
!     &        global_row_length, global_rows,
!     &        fld_type_p, halo_type_extended,
!     &        0, gc_all_proc_group,
!     &        ICODE, CMessage)

        If (me == 0) Then
          ! Set local random number array on pe0
          Do j=1-halo_j,rows+halo_j
            Do i=1-halo_i, row_length+halo_i
              ran(i,j)=ran_pe0(i,j)
            End Do
          End Do

          ! Determine maximum and minimum random numbers for diagnostic
          maxran = 0.0
          minran = 0.0
          Do i = 1-halo_i,global_row_length+halo_i
            Do j = 1-halo_j,global_rows+halo_j
              If (ran_pe0(i,j)  >   maxran) maxran = ran_pe0(i,j)
              If (ran_pe0(i,j)  <   minran) minran = ran_pe0(i,j)
            End Do
          End Do

           Do d=1,n_proc-1
            Do j=1-halo_j,g_rows(d)+halo_j
              gj = g_datastart(2,d) + j - 1
              Do i=1-halo_i, g_row_length(d)+halo_i
                gi = g_datastart(1,d) + i - 1
                rbuf(i,j)=ran_pe0(gi,gj)
              End Do
            End Do

            ! Send (decompose) random number field to processors
            Call gc_rsend(100+d,                                        &
     &            (global_rows+halo_j*2)*(global_row_length+halo_i*2),  &
     &                  d, info, rbuf, rbuf)

          End Do  ! d
        End If  !me=0

        ! All processors receive random number field from pe0
        If (me  /=  0) Then
          Call gc_rrecv(100+me,                                         &
     &            (global_rows+halo_j*2)*(global_row_length+halo_i*2),  &
     &                    0, info, rbuf, rbuf)
          Do j=1-halo_j,rows+halo_j
            Do i=1-halo_i, row_length+halo_i
              ran(i,j)=rbuf(i,j)
            End Do
          End Do
        End If      ! me==0


! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &         ran, row_length, rows,                                   &
     &         1, halo_i, halo_j, fld_type_p, .false.)


        If (me  ==  0)                                                  &
     &    Write(Unit=6,Fmt='(A50,F6.3,A13,F6.3)')                       &
     &      '  Gridscale random numbers generated with minimum ',       &
     &       minran,' and maximum ',maxran


        !-------------------------------------------------------------
        ! Perturbations are applied to all levels within the height
        ! range specified in the namelist by perturb_height(1) and (2)
        ! perturb_height(1) is lowest height to apply perturbations
        ! perturb_height(2) is highest height to apply perturbations
        !-------------------------------------------------------------

        Do k = 1, model_levels
          z_theta_levels(k) = eta_theta_levels(k)*height_domain
        End Do

        ! If highest perturb_height is below first model level
        ! then reset to height of first model level
        ! (plus a bit to ensure IF test in following section
        !  is satisfied)
        If (perturb_height(2)  <   z_theta_levels(1)) Then
          perturb_height(2) = z_theta_levels(1)+0.00000001
        End If

        ! Work out bottom and top levels to apply perturbations
        levbot = 0
        levtop = 0
        Do k = 1, model_levels
          If ((z_theta_levels(k)  >=  perturb_height(1)) .AND.          &
     &        (z_theta_levels(k)  <=  perturb_height(2))) Then
            If (levbot  ==  0) levbot = k
            levtop = k
          End If
        End Do


        !-------------------------------------------------------------
        ! Apply grid scale random perturbations
        !-------------------------------------------------------------

        ! Generate an array of random numbers from 0 to 1
        ! Call IDL_Get_Random_Array

        ! Apply perturbations to potential temperature field
        If (L_perturb_t) Then
          Do k = levbot, levtop
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                theta(i,j,k) = theta(i,j,k)                             &
     &                  + 2.*perturb_magnitude_t*(ran(i,j)-0.5)
              End Do
            End Do
          End Do
          ! Write out information on applied perturbations
          If (me  ==  0)                                                &
     &      Write(6,*) 'Gridscale random temperature +-',               &
     &       perturb_magnitude_t,' perturbations'
        End If

        ! Generate a different array of random numbers from 0 to 1
        ! if theta and q perturbations are to be uncorrelated
        ! If (L_perturb_q .AND. .NOT. L_perturb_correlate_tq) Then
        !   Call IDL_Get_Random_Array
        ! End If

        ! Apply perturbations to water vapour field
        If (L_perturb_q) Then
          Do k = levbot, levtop
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                q(i,j,k) = q(i,j,k)                                     &
     &                   + 2.*perturb_magnitude_q*(ran(i,j)-0.5)
              End Do
            End Do
          End Do
          ! Write out information on applied perturbations
          If (me  ==  0)                                                &
     &      Write(6,*) 'Gridscale random moisture +-',                  &
     &       perturb_magnitude_q,' perturbations'
        End If



        !-------------------------------------------------------------
        ! Write out further information on applied perturbations
        !-------------------------------------------------------------

        If (me  ==  0) Then
          Write(6,*) '  applied from model level ',levbot,' to ',levtop
          Write(6,*) '  (Height: ',perturb_height(1),'m to ',           &
     &                     perturb_height(2),'m)'
          Write(6,*) '  Theta and moisture perturbations correlated'
          Write(6,*) '  Perturbations correlated in the vertical'
        End If


      IF (lhook) CALL dr_hook('IDL_RANDOM_PERTURB',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_random_perturb

      !  End subroutine IDL_random_perturb
