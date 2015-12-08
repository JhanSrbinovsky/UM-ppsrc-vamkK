! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine perturb_theta
      SUBROUTINE perturb_theta(                                         &
     &                         theta, row_length, rows, model_levels,   &
     &                         global_row_length, global_rows,          &
     &                         model_domain, at_extremity,              &
     &                         offx, offy, IntRand_Seed, l_datastart    &
     &                        )

! Purpose:
!          Perturb theta at the bit level using a random number  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
      Integer, Intent(In) :: ROW_LENGTH     ! No of points per local row
      Integer, Intent(In) :: ROWS           ! No of local (theta) rows
      Integer, Intent(In) :: MODEL_LEVELS   ! No of model levels
      Integer, Intent(In) :: Offx    ! standard halo size in East-West
      Integer, Intent(In) :: Offy    ! standard halo size in North-South
      Logical, Intent(In) :: At_extremity(4)
      Integer, Intent(In) :: model_domain
      Integer, Intent(In) :: IntRand_Seed   ! Seed for random numbers
      Integer, Intent(In) :: l_datastart(2)
      Integer, Intent(In) :: global_ROW_LENGTH    
      Integer, Intent(In) :: global_ROWs

      Real, Intent (InOut) ::                                           &
     &  theta(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &        model_levels) 

! Local Variables
      Integer :: i,j,k    ! loop variables
      Integer :: seed_len ! Length of random number seed vector
      Integer :: j_start,j_end
      Integer :: gi,gj

      Integer, Allocatable :: seed(:) ! Vector random number seed

      REAL :: RandomNumbers(global_row_length,global_rows,model_levels)
      REAL :: eps_machine
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



!----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('PERTURB_THETA',zhook_in,zhook_handle)
      if(model_domain==mt_global .and. at_extremity(Psouth))then
        j_start=2
      else 
        j_start=1
      endif
      if(model_domain==mt_global .and. at_extremity(PNorth))then
        j_end=rows-1
      else 
        j_end=rows
      endif

      write(6,*)'PERTURBING THETA FIELD'
      eps_machine=epsilon(1.0)
      write(6,*)' MACHINE epsilon: ', eps_machine
      write(6,*)' IntRand_Seed:    ', IntRand_Seed

      ! Initialise random number seed.
      ! (For some reason, the seeds n*4, n*4+1, n*4+2 and n*4+3 all
      ! give the same random numbers on the IBM. To get around this,
      ! we multiply the supplied seed by 4.)
      CALL RANDOM_SEED (SIZE = seed_len)
      ALLOCATE (seed(seed_len))
      seed(:) = IntRand_seed * 4
      CALL RANDOM_SEED (PUT = seed)
      DEALLOCATE (seed)

      ! Get random numbers in the range 0-1:
      CALL RANDOM_NUMBER (RandomNumbers)

      Do k=1,model_levels
        Do j=j_start,j_end
          gj=j+l_datastart(2)-1
          Do i=1, row_length
            gi=i+l_datastart(1)-1
            theta(i,j,k) = theta(i,j,k) +                               &
     &           theta(i,j,k) * (2.0*RandomNumbers(gi,gj,k) - 1.0) *    &
     &             eps_machine
          End Do
        End Do
      End Do

      IF (lhook) CALL dr_hook('PERTURB_THETA',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE perturb_theta
