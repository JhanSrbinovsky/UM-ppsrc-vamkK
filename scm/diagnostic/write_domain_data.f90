! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Writes file domain.dat for the SCM diagnostics

SUBROUTINE write_domain_data                                                  &
  ( row_length, rows, model_levels, z_top_of_model                            &
  , first_constant_r_rho_level, orog, eta_theta,eta_rho )

  IMPLICIT NONE

! Description:
!   Write file domain.dat which will contain information on the
!   size and dimensions of the model. This information will be
!   necessary for making sense of some of the information in the
!   diagnostic output data files.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran 90

  INTEGER :: row_length,rows       ! In Horizontal domain
  REAL    :: orog(row_length*rows) ! In Orography height
  INTEGER :: model_levels          ! In The total no. of model levels
  REAL    :: z_top_of_model        ! In The height of the "top of the
                                   !    atmosphere"

  INTEGER ::first_constant_r_rho_level
                                   ! In Lowest rho level that has
                                   !    constant height
  REAL ::                     &
    eta_theta(model_levels+1) & ! In The etas of the theta and
  , eta_rho(model_levels)       !    rho levels

  INTEGER :: i                    ! Counter
  INTEGER, PARAMETER :: unit = 10 ! Temporary unit no. to write to

  OPEN (unit=unit,file='domain.dat')

  ! write out the size of the model domain
  WRITE(unit,'(I4,1X,I4,1X,I4,1X,F9.2,1X,I4)')                      &
       row_length,rows,model_levels,                                &
       z_top_of_model,first_constant_r_rho_level

  ! write out the eta levels
  WRITE(unit,*)(eta_theta(i),i=1,model_levels+1)
  WRITE(unit,*)(eta_rho(i),i=1,model_levels)

  ! write out the orography data so that PV-wave programs can
  ! calculate heights from eta-levels
  WRITE(unit,*)(orog(i),i=1,row_length*rows)

  CLOSE(unit)
  RETURN

END SUBROUTINE write_domain_data

