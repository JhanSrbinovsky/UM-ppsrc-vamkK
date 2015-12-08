! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Check_sl_domain.

      Subroutine Check_sl_domain(                                       &
     &                        model_domain, depart_phi, depart_lambda,  &
     &                        row_length, rows_depart, model_levels,    &
     &                        domain_size_x, domain_size_y,             &
     &                        max_lambda, min_lambda,                   &
     &                        max_phi, min_phi)

! Purpose:
!          Checks data lies inside model domain.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE conversions_mod, ONLY: pi
      USE domain_params
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                     ! Dimension of u_adv array in i direction.
     &, rows_depart                                                     &
                     ! Dimension of depart arrays in j direction.
     &, model_levels                                                    &
                     ! Dimension of u_adv array in k direction.
     &, model_domain ! holds integer code for model domain

      Real                                                              &
     &  max_lambda                                                      &
     &, max_phi                                                         &
     &, min_lambda                                                      &
     &, min_phi

      Real :: domain_size_x
      Real :: domain_size_y

! Arguments with Intent IN/OUT.

      Real                                                              &
                      ! Departure point co-ordinates.
     &  depart_lambda (row_length, rows_depart, model_levels)           &
     &, depart_phi (row_length, rows_depart, model_levels)

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k      ! Loop indices

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines: None

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Check points lie inside model domain.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CHECK_SL_DOMAIN',zhook_in,zhook_handle)
      If (model_domain  ==  mt_Global) Then
! check horizontal
        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_phi(i,j,k)  <   min_phi ) Then
                depart_phi(i,j,k) = -pi -  depart_phi(i,j,k)
                depart_lambda(i,j,k) = Pi + depart_lambda(i,j,k)
              Else If (depart_phi(i,j,k)  >   max_phi ) Then
                depart_phi(i,j,k) = Pi -  depart_phi(i,j,k)
                depart_lambda(i,j,k) = Pi + depart_lambda(i,j,k)
              End If
              If (depart_lambda(i,j,k)  <   min_lambda ) Then
                depart_lambda(i,j,k) = Pi*2. + depart_lambda(i,j,k)
              Else If (depart_lambda(i,j,k)  >=  max_lambda ) Then
                depart_lambda(i,j,k) = depart_lambda(i,j,k) - Pi*2.
              End If
            End Do
          End Do
        End Do

      Else If (model_domain  ==  mt_LAM) Then

! check horizontal in the LAM at u points.

        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_phi(i,j,k)  <   min_phi ) Then
                depart_phi(i,j,k) = min_phi
              Else If (depart_phi(i,j,k) >  max_phi)Then
                depart_phi(i,j,k) = max_phi
              End If
              If (depart_lambda(i,j,k)  <   min_lambda ) Then
                depart_lambda(i,j,k) = min_lambda
              Else If(depart_lambda(i,j,k) >  max_lambda) Then
                depart_lambda(i,j,k) = max_lambda
              End If
            End Do
          End Do
        End Do

      Else If (model_domain  ==  mt_cyclic_LAM) Then

! check horizontal in the periodic in x LAM at u points.

        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_phi(i,j,k)  <   min_phi ) Then
                depart_phi(i,j,k) = min_phi
              Else If (depart_phi(i,j,k) >  max_phi)Then
                depart_phi(i,j,k) = max_phi
              End If
              If (depart_lambda(i,j,k)  <   min_lambda ) Then
                depart_lambda(i,j,k) = depart_lambda(i,j,k)             &
     &                               + domain_size_x             
              Else If(depart_lambda(i,j,k) >= max_lambda) Then
                depart_lambda(i,j,k) = depart_lambda(i,j,k)             &
     &                               - domain_size_x            
              End If
            End Do
          End Do
        End Do


      Else If (model_domain  ==  mt_bi_cyclic_lam ) Then
! check horizontal in the periodic in both x and y LAM at u points

        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_phi(i,j,k)  <   min_phi) then
                depart_phi(i,j,k) = depart_phi(i,j,k)                   &
     &                            + domain_size_y     
              Elseif (depart_phi(i,j,k)  >   max_phi) then
                depart_phi(i,j,k) = depart_phi(i,j,k)                   &
     &                            - domain_size_y     
              Endif
              If (depart_lambda(i,j,k)  <   min_lambda) then
                depart_lambda(i,j,k) = depart_lambda(i,j,k)             &
     &                               + domain_size_x          
              Elseif (depart_lambda(i,j,k)  >   max_lambda) then
                depart_lambda(i,j,k) = depart_lambda(i,j,k)             &
     &                               - domain_size_x          
              Endif
            Enddo
          Enddo
        Enddo

      End if

! End of routine.
      IF (lhook) CALL dr_hook('CHECK_SL_DOMAIN',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Check_sl_domain

