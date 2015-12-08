! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_Elliptic_Operator_2B
      SUBROUTINE GCR_Elliptic_Operator_2B(                              &
     &                                field, row_length, rows, n_rows,  &
     &                                model_levels, model_domain,       &
     &                                first_constant_r_rho_level,       &
     &                                first_constant_r_rho_level_m1,    &
     &                                eta_theta_levels,                 &
     &                                eta_rho_levels,                   &
     &                                FV_cos_theta_latitude,            &
     &                                lambda_p, phi_p, lambda_u, phi_v, &
     &                                dlambda_p, dphi_p,                &
     &                                dlambda_u, dphi_v,                &
     &                                recip_dlamp, recip_dphip,         &
     &                                recip_dlamu, recip_dphiv,         &
     &                                wt_lambda_p, wt_phi_p,            &
     &                                wt_lambda_u, wt_phi_v,            &
     &                                HM_Cxx1, HM_Cxx2, HM_Cyy1,        &
     &                                HM_Cyy2, HM_Czz, HM_Cz,           &
     &                                HM_Cxz, HM_Cyz,  HM_Cxp, HM_Cyp,  &
     &                                HM_Cxy1, HM_Cxy2,                 &
     &                                HM_Cyx1, HM_Cyx2,                 &
     &                                HM_C2n,HM_C3, HM_C4, HM_C5,       &
     &                                weight_upper, weight_lower,       &
     &                                L_of_field,                       &
     &                                offx, offy, halo_i, halo_j,       &
     &                                at_extremity, proc_row_group,     &
     &                                global_row_length,                &
     &                                i_start, i_stop, j_start, j_stop, &
     &                                j_begin, j_end,                   &
     &                                L_regular )

! Purpose:
!          Applies Elliptic Operator L to input field.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE global_2d_sums_mod, ONLY : global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v-rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, offx                                                            &
                         ! small halo size
     &, offy                                                            &
                         ! small halo size
     &, halo_i                                                          &
                         ! large halo size
     &, halo_j                                                          &
                         ! large halo size
     &, model_domain                                                    &
                         ! holds integer code for model domain
     &, first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, first_constant_r_rho_level_m1 ! value used to dimension
                                      ! arrays, max of (1 and
                                      ! first_constant_r_rho_level)
! input field
      Real                                                              &
     &  field(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)

! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

! Coefficients of the elliptic operator
      Real                                                              &
     & HM_Cxx1(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)  &
     &,HM_Cxx2(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)  &
     &,HM_Cxy1(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)  &
     &,HM_Cxy2(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)  &
     &,HM_Cyy1(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)  &
     &,HM_Cyy2(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)  &
     &,HM_Cyx1(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)  &
     &,HM_Cyx2(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)  &
     &,HM_Czz (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)  &
     &,HM_Cz (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)   &
     &,HM_C2n(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)   &
     &,HM_C3 (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)   &
     &,HM_C4 (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)   &
     &,HM_C5 (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)   &
     &, HM_Cxz (1-offx:row_length+offx, 1-offy:rows+offy,               &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyz (1-offx:row_length+offx, 1-offy:rows+offy,               &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxp (1-offx:row_length+offx, 1-offy:rows+offy,               &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyp (1-offx:row_length+offx, 1-offy:rows+offy,               &
     &         first_constant_r_rho_level_m1)

! Trigonometric functions
      Real                                                              &
     &  FV_cos_theta_latitude (1-offx:row_length+offx,1-offy:rows+offy)
                              ! finite volume cosine array

! vertical co-ordinate information
      Real                                                              &
     &  eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

!  VarRes horizontal co-ordinate spacing.
      Real                                                              &
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_p(1-halo_i:row_length+halo_i, 1-halo_j:  rows+halo_j)       &
     &, phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, dlambda_p(1-halo_i : row_length+halo_i)                         &
     &, dlambda_u(1-halo_i : row_length+halo_i)                         &
     &, dphi_p(1-halo_i:row_length+halo_i, 1-offy:rows+offy)            &
     &, dphi_v(1-halo_i:row_length+halo_i, 1-offy:n_rows+offy)          &
     &, recip_dlamp(1-halo_i : row_length + halo_i)                     &
     &, recip_dlamu(1-halo_i : row_length + halo_i)                     &
     &, recip_dphip(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
     &, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

!   parallel variables Local
      Integer info

      Integer                                                           &
     &  proc_row_group                                                  &
     &, global_row_length                                               &
     &, i_start, i_stop                                                 &
                                     ! loop bounds set in PE_Helmholtz
     &, j_start, j_stop                                                 &
                                     ! loop bounds set in PE_Helmholtz
     &, j_begin, j_end

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular               ! true for regular resolution


! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

! Arguments with Intent OUT. ie: variables Output only

      Real                                                              &
     &  L_of_field(row_length,rows,model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k      ! Loop indices

      Real                                                              &
     &  term_1                                                          &
                   ! a term in the equation
     &, term_2                                                          &
                   ! another term in the equation
     &, factor_1                                                        &
                   ! another bit of the equation
     &, factor_2                                                        &
                   ! and another bit of the equation
     &, factor_3                                                        &
                   ! and yet another bit of the equation
     &, interp_upper                                                    &
                     ! weights for vertical averaging in eta
     &, interp_lower                                                    &
                     ! weights for vertical averaging in eta
     &, term_lower                                                      &
     &, recip_g_row_len

! Local arrays

      Real                                                              &
     &  first_deriv_x(1-offx:row_length, 1-offy:rows+offy)              &
     &, x_term(1-offx:row_length, rows, model_levels)                   &
     &, first_deriv_y(1-offx:row_length+offx,1-offy: rows)              &
     &, y_term(row_length, 1-offy:rows, model_levels)                   &
     &, interpx(1-offx:row_length+offx,1-offy:rows+offy)                &
     &, interpy(row_length,1-offy:rows)                                 &
     &, term_upper(1-offx:row_length+offx,1-offy:rows+offy)             &
     &, l_s_poles(row_length,model_levels)                              &
     &, sum_s(model_levels)                                             &
     &, l_n_poles(row_length,model_levels)                              &
     &, sum_n(model_levels)                                             &
     &, factor_1_p(model_levels)

       Real                                                             &
                  !  VarRes
     &  work1( 1-offx:row_length+offx, 1-offy:rows+offy )               &
     &, work2( 1-offx:row_length+offx, 1-offy:rows+offy )

      integer :: ifirst
      real :: termj_lower, termj_upper
      logical :: L_first_level

! No External Routines:
       External  swap_bounds

       real Two_norm

       INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
       INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
       REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1.   Calculate Elliptic Operator applied to input field
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GCR_ELLIPTIC_OPERATOR_2B',zhook_in,zhook_handle)
      recip_g_row_len = 1. / global_row_length

! ----------------------------------------------------------------------
! Section 1.1   Calculate horizontal first derivatives.
! ----------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,2) DEFAULT(NONE)                     &
!$OMP& PRIVATE(k,i,j,first_deriv_x,first_deriv_y,interpx,               &
!$OMP&        termj_lower,termj_upper,work1,work2,L_first_level)        &
!$OMP& SHARED(L_of_field, field, x_term, y_term,model_levels,           &
!$OMP&        j_begin,j_end,i_start,i_stop,L_regular,j_start,j_stop,    &
!$OMP&        recip_dlamp,recip_dlamu,recip_dphip,recip_dphiv,          &
!$OMP&        wt_lambda_p,wt_lambda_u,wt_phi_p,wt_phi_v,                &
!$OMP&        first_constant_r_rho_level,HM_C2n,weight_upper,           &
!$OMP&        weight_lower,HM_Cxp,HM_Cyp,HM_Cxx2,                       &
!$OMP&        HM_Cxy2,HM_Cyx2,HM_Cyy2,HM_C4,HM_Cxx1,HM_Cyy1,            &
!$OMP&        HM_Cxy1,HM_Cyx1,FV_cos_theta_latitude) 
      Do k = 1, model_levels
        L_first_level = (k==1)

        If ( k < first_constant_r_rho_level ) Then

! first calculate first derivative with respect to x on constant r
! surface.

! calculate vertical derivative at all points.
! store in interpx to save re-calculation in y step.

            Do j = j_begin - 1, j_end + 1
              Do i = i_start-1, i_stop + 1
                termj_upper      = HM_C2n(i,j,k) *                       &
     &                               (field(i,j,k+1) - field(i,j,k))
               if(L_first_level)then
                interpx(i,j) = weight_upper(i,j,k) * termj_upper     
               else
                termj_lower      = HM_C2n(i,j,k-1) *                       &
     &                               (field(i,j,k) - field(i,j,k-1))
                interpx(i,j) = weight_upper(i,j,k) * termj_upper +  &
     &                         weight_lower(i,j,k) * termj_lower
               endif
              End Do
            End Do

! calculate first derivative on constant surface.
!!!   polar rows not included

          if( L_regular ) then

            Do j = j_begin - 1, j_end + 1
              Do i = i_start-1, i_stop
                first_deriv_x(i,j) = ( field(i+1,j,k) - field(i,j,k) ) -&
     &                                 HM_Cxp(i,j,k) * .5 *             &
     &                                ( interpx(i,j) + interpx(i+1,j) )
              EndDo
            EndDo

          else !  variable resolution

            Do j = j_begin - 1, j_end + 1
              Do i = i_start-1, i_stop
                first_deriv_x(i,j) = ( field(i+1,j,k) - field(i,j,k) ) *&
     &                                  recip_dlamp(i) - HM_Cxp(i,j,k) *&
     &                                (wt_lambda_p(i+1) * interpx(i,j) +&
     &                        (1.0 - wt_lambda_p(i+1)) * interpx(i+1,j))
              EndDo
            EndDo
          endif ! L_regular

! now calculate first derivative with respect to y on constant r
! surface.
          if( L_regular ) then

            Do j = j_begin - 1, j_end
              Do i = i_start-1, i_stop + 1
                first_deriv_y(i,j) = ( field(i,j+1,k) - field(i,j,k) ) -&
     &                                             0.5 * HM_Cyp(i,j,k) *&
     &                               ( interpx(i,j) + interpx(i,j+1) )
              EndDo
            EndDo

          else !  variable resolution

            Do j = j_begin - 1, j_end
              Do i = i_start-1, i_stop + 1
                first_deriv_y(i,j) = ( field(i,j+1,k) - field(i,j,k) ) *&
     &                                recip_dphip(i,j) - HM_Cyp(i,j,k) *&
     &                             ( wt_phi_p(i,j+1) * interpx(i,j) +   &
     &                       (1.0 - wt_phi_p(i,j+1)) * interpx(i,j+1) )
              EndDo
            EndDo

          endif ! L_regular

        Else  !   k >= first_constant_r_rho_level

! surface is constant.
! first calculate first derivative with respect to x on constant r
! surface.
          if( L_regular ) then

            Do j = j_begin - 1, j_end + 1
              Do i = i_start-1, i_stop
                first_deriv_x(i,j) = (field(i+1,j,k) - field(i,j,k))
              End do
            End do

          else !  variable resolution

            Do j = j_begin - 1, j_end + 1
              Do i = i_start-1, i_stop
                first_deriv_x(i,j) = (field(i+1,j,k) - field(i,j,k)) *  &
     &                                recip_dlamp(i)
              End do
            End do
          endif ! L_regular

! now calculate first derivative with respect to y on constant r
! surface.
          if( L_regular ) then

            Do j = j_begin - 1, j_end
              Do i = i_start-1, i_stop + 1
                first_deriv_y(i,j) = ( field(i,j+1,k) - field(i,j,k) )
              End do
            End do

          else !  variable resolution

            Do j = j_begin - 1, j_end
              Do i = i_start-1, i_stop + 1
                first_deriv_y(i,j) = ( field(i,j+1,k) - field(i,j,k) ) *&
     &                                 recip_dphip(i,j)
              End do
            End do

          endif ! L_regular

        End If  ! k < first_constant_r_rho_level

! ----------------------------------------------------------------------
! Section 1.2   Calculate the term to be differenced and averaged in x,
!               and the term to be differenced and averaged in y.
! ----------------------------------------------------------------------

        if( L_regular ) then
! x_term
          Do j = j_begin, j_end
            Do i = i_start - 1, i_stop
              x_term(i,j,k) = HM_Cxx2(i,j,k) * first_deriv_x(i,j) +     &
     &                                       .25 * HM_Cxy1(i,j,k) *     &
     &                    ( HM_Cxy2(i,j,k) * first_deriv_y(i,j) +       &
     &                      HM_Cxy2(i+1,j,k) * first_deriv_y(i+1,j) +   &
     &                      HM_Cxy2(i,j-1,k) * first_deriv_y(i,j-1) +   &
     &                      HM_Cxy2(i+1,j-1,k) * first_deriv_y(i+1,j-1))
            End Do
          End Do
! y_term
          Do j = j_begin-1, j_end
            Do i = i_start, i_stop
            y_term(i,j,k) = HM_Cyy2(i,j,k) * first_deriv_y(i,j) -       &
     &                                       .25 * HM_Cyx1(i,j,k) *     &
     &                    ( HM_Cyx2(i,j,k) * first_deriv_x(i,j) +       &
     &                      HM_Cyx2(i-1,j,k) * first_deriv_x(i-1,j) +   &
     &                      HM_Cyx2(i,j+1,k) * first_deriv_x(i,j+1) +   &
     &                      HM_Cyx2(i-1,j+1,k) * first_deriv_x(i-1,j+1))
            End Do
          End Do

        else !  variable resolution

! x_term
          Do j = j_begin, j_end
            Do i = i_start - 1, i_stop + 1
              work1(i,j) = wt_phi_v(i,j) *                              &
     &                       HM_Cxy2(i,j-1,k) * first_deriv_y(i,j-1) +  &
     &                       (1.0 - wt_phi_v(i,j)) *                    &
     &                       HM_Cxy2(i,j,k) * first_deriv_y(i,j)
            End Do
            Do i = i_start - 1, i_stop
              x_term(i,j,k) = HM_Cxx2(i,j,k) * first_deriv_x(i,j) +     &
     &                                             HM_Cxy1(i,j,k) *     &
     &                            ( wt_lambda_p(i+1) * work1(i,j) +     &
     &                      (1.0 - wt_lambda_p(i+1)) * work1(i+1,j) )
           End Do
          End Do
! y_term.
          Do j = j_begin-1, j_end + 1
            Do i = i_start, i_stop
              work1(i,j) =  wt_lambda_u(i) *                            &
     &                        HM_Cyx2(i-1,j,k) * first_deriv_x(i-1,j) + &
     &                        (1.0 - wt_lambda_u(i)) *                  &
     &                        HM_Cyx2(i,j,k)   * first_deriv_x(i,j)
            End Do
          End Do
          Do j = j_begin-1, j_end
            Do i = i_start, i_stop
            y_term(i,j,k) = HM_Cyy2(i,j,k) * first_deriv_y(i,j) -       &
     &                                           HM_Cyx1(i,j,k) *       &
     &                             ( wt_phi_p(i,j+1) * work1(i,j) +     &
     &                       (1.0 - wt_phi_p(i,j+1)) * work1(i,j+1) )
            End Do
          End Do

        endif ! L_regular

! ----------------------------------------------------------------------
! Section 2.1   Calculate constant terms of operator applied to field.
! ----------------------------------------------------------------------

        Do j = j_start, j_stop
          Do i = i_start, i_stop
            L_of_field(i,j,k) = - HM_C4(i,j,k) * field(i,j,k) *         &
     &                           FV_cos_theta_latitude(i,j)
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.2   Calculate lambda derivative terms of operator applied
!               to field, including d/dx(d/dy) terms.
! ----------------------------------------------------------------------

! Calculate Operator
        if( L_regular ) then

          Do j = j_begin, j_end
            Do i = i_start, i_stop
              L_of_field(i,j,k) = L_of_field(i,j,k) +                   &
     &                             ( x_term(i,j,k) * HM_Cxx1(i,j,k) -   &
     &                          x_term(i-1,j,k) * HM_Cxx1(i-1,j,k)) +   &
     &                              (y_term(i,j,k) * HM_Cyy1(i,j,k) -   &
     &                           y_term(i,j-1,k) * HM_Cyy1(i,j-1,k) )
            End Do
          End Do

        else !  variable resolution

          Do j = j_begin, j_end
            Do i = i_start-1, i_stop
              work1(i,j) = x_term(i,j,k) * HM_Cxx1(i,j,k)
            End Do
          End Do
          Do j = j_begin-1, j_end
            Do i = i_start, i_stop
              work2(i,j) = y_term(i,j,k) * HM_Cyy1(i,j,k)
            End Do
          End Do

          Do j = j_begin, j_end
            Do i = i_start, i_stop
              L_of_field(i,j,k) = L_of_field(i,j,k) +                   &
     &                                ( work1(i,j) - work1(i-1,j) ) *   &
     &                                         recip_dlamu(i-1) +       &
     &                                ( work2(i,j) - work2(i,j-1) ) *   &
     &                                         recip_dphiv(i,j-1)
            End Do
          End Do

        endif ! L_regular

      End Do ! k = 1, model_levels
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Section 2.3.2 Poles in Global Model
! ----------------------------------------------------------------------


! average the value and add on constant term, note any other polar
! point will do as all values are the same.

      If (model_domain == mt_global) Then

        If (at_extremity(PSouth)) Then
          If( L_regular ) then
          Do k = 1, model_levels
            Do i = 1, row_length
              l_s_poles(i,k) = y_term(i,1,k) * HM_Cyy1(i,1,k)
            End Do
          End Do
          Else !  variable resolution
          Do k = 1, model_levels
            Do i = 1, row_length
              l_s_poles(i,k) = y_term(i,1,k) * HM_Cyy1(i,1,k)           &
     &           * dlambda_p(i) * recip_dphiv(i,2)
            End Do
          End Do
          Endif ! L_regular

          CALL global_2d_sums(l_s_poles, row_length, 1, 0, 0,           &
                              model_levels, sum_s, proc_row_group)

          Do k = 1, model_levels
            L_of_field(1,1,k) = L_of_field(1,1,k) +                     &
     &                          sum_s(k) * recip_g_row_len
! Copy answer at one point to all others
            Do i =2,row_length
              L_of_field(i,1,k) = L_of_field(1,1,k)
            End Do
          End Do
        End If !  at_extremity(PSouth)

        If (at_extremity(PNorth)) Then
          If( L_regular ) then
          Do k = 1, model_levels
            Do i = 1, row_length
              l_n_poles(i,k) = y_term(i,rows-1,k) * HM_Cyy1(i,rows-1,k)
            End Do
          End Do
          Else !  variable resolution
          Do k = 1, model_levels
            Do i = 1, row_length
              l_n_poles(i,k) = y_term(i,rows-1,k) * HM_Cyy1(i,rows-1,k) &
     &           * dlambda_p(i) * recip_dphiv(i,n_rows)
            End Do
          End Do
          Endif ! L_regular

          CALL global_2d_sums(l_n_poles, row_length, 1, 0, 0,           &
                              model_levels, sum_n, proc_row_group)

          Do k = 1, model_levels
            L_of_field(1,rows,k) = L_of_field(1,rows,k)                 &
     &                             - sum_n(k) * recip_g_row_len
! Copy answer at one point to all others
            Do i =2,row_length
              L_of_field(i,rows,k) = L_of_field(1,rows,k)
            End Do
          End Do
        End If  ! at_extremity(PNorth)

      End If ! model_domain == mt_global

! ----------------------------------------------------------------------
! Section 2.4   Calculate vertical derivative terms of operator applied
!               to field.
! ----------------------------------------------------------------------

      If ( model_levels > 1 ) Then

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP& PRIVATE(k,i,j,factor_1,factor_2,factor_3,term_1,term_2)          &
!$OMP& SHARED(j_start, j_stop,field,l_of_field,eta_rho_levels,          &
!$OMP&        eta_theta_levels,fv_cos_theta_latitude,model_levels,      &
!$OMP&    weight_upper,weight_lower,hm_c3,hm_czz,hm_cz,i_start, i_stop)
        Do k = 1, model_levels

          if (k == 1) then

            factor_1 = 1. / (eta_theta_levels(k) -                      &
     &                       eta_theta_levels(k-1))
            factor_2 = 1. / (eta_rho_levels(k+1) - eta_rho_levels(k))

            Do j = j_start, j_stop
              Do i = i_start, i_stop
! Calculate upper derivative
                term_1 = ( field(i,j,k+1) - field(i,j,k) ) * factor_2
                L_of_field(i,j,k) = L_of_field(i,j,k) +                 &
     &                              FV_cos_theta_latitude(i,j) *        &
     &                              ( factor_1 * term_1 *               &
     &                               HM_Czz(i,j,k) + HM_C3(i,j,k) *     &
     &                               term_1 * HM_Cz(i,j,k) )
              End Do
            End Do

          else if (k == model_levels) then

            factor_1 = 1. / (eta_theta_levels(k) -                      &
     &                       eta_theta_levels(k-1))
            factor_3 = 1. / (eta_rho_levels(k) - eta_rho_levels(k-1))

            Do j = j_start, j_stop
              Do i = i_start, i_stop
! Calculate lower derivative.
                term_2 = ( field(i,j,k) - field(i,j,k-1) ) * factor_3
                L_of_field(i,j,k) = L_of_field(i,j,k) -                 &
     &                              FV_cos_theta_latitude(i,j) *        &
     &                       ( factor_1 * term_2 * HM_Czz(i,j,k-1) -    &
     &                      HM_C3(i,j,k) * term_2 * HM_Cz(i,j,k-1) *    &
     &                                       weight_lower(i,j,k) )
              End Do
            End Do

          else ! 1 < k < model_levels

            factor_1 = 1. / (eta_theta_levels(k) -                      &
     &                       eta_theta_levels(k-1))
            factor_2 = 1. / (eta_rho_levels(k+1) - eta_rho_levels(k))
            factor_3 = 1. / (eta_rho_levels(k) - eta_rho_levels(k-1))

            Do j = j_start, j_stop
              Do i = i_start, i_stop
                L_of_field(i,j,k) = L_of_field(i,j,k) +                 &
     &                              FV_cos_theta_latitude(i,j) *        &
     &                        ( field(i,j,k+1) - field(i,j,k) ) *       &
     &                    factor_2 * ( factor_1 * HM_Czz(i,j,k) +       &
     &                              HM_C3(i,j,k) * HM_Cz(i,j,k) *       &
     &                                      weight_upper(i,j,k) )   &
     &                           -  FV_cos_theta_latitude(i,j) *        &
     &                       ( field(i,j,k) - field(i,j,k-1) ) *        &
     &                    factor_3 * ( factor_1 * HM_Czz(i,j,k-1) -     &
     &                              HM_C3(i,j,k) * HM_Cz(i,j,k-1) *     &
     &                                      weight_lower(i,j,k) )
              End Do
            End Do

          endif !  k == 1

        End Do !  k = 1, model_levels
!$OMP END PARALLEL DO

      End If  ! model_levels > 1

! ----------------------------------------------------------------------
! Section 2.5   Calculate mixed vertical derivative terms of operator
!               applied to field.
! ----------------------------------------------------------------------

      If ( first_constant_r_rho_level > 1 ) Then

!$OMP PARALLEL DEFAULT(NONE)                                            &
!$OMP& PRIVATE(ifirst,k,j,i,interp_upper,interp_lower,interpx,interpy,  &
!$OMP&        term_lower,term_upper)                                    &
!$OMP& SHARED(first_constant_r_rho_level,eta_theta_levels,HM_C5,        &
!$OMP&     eta_rho_levels,j_begin,j_end,i_start,i_stop,wt_lambda_u,     &
!$OMP&     x_term,HM_Cxz,y_term,HM_Cyz,factor_1_p,dlambda_p,wt_phi_v,   &
!$OMP&     L_of_field,FV_cos_theta_latitude,model_domain,at_extremity,  &
!$OMP&     row_length,l_s_poles,l_n_poles,n_rows,L_regular)
        ifirst=0
!$OMP DO SCHEDULE(STATIC) 
        Do k = 1, first_constant_r_rho_level

          If (k == 1 ) Then
           ifirst=1

! Interpolate quantities to levels either side in vertical
            interp_upper = (eta_theta_levels(k) - eta_rho_levels(k)) /  &
     &                     (eta_rho_levels(k+1) - eta_rho_levels(k))
            interp_lower = 1.0 - interp_upper

            Do j = j_begin, j_end
              Do i = i_start - 1, i_stop
                interpx(i,j) = ( interp_upper * x_term(i,j,k+1) +       &
     &                           interp_lower * x_term(i,j,k) ) *       &
     &                                          HM_Cxz(i,j,k)
              End Do
            End Do

            Do j = j_begin-1, j_end
              Do i = i_start, i_stop
                interpy(i,j) = (interp_upper * y_term(i,j,k+1) +        &
     &                          interp_lower * y_term(i,j,k) ) *        &
     &                                         HM_Cyz(i,j,k)
              End Do
            End Do
! Calculate 1. / delta eta.
            factor_1_p(k) = 1. / (eta_theta_levels(k) -                 &
     &                            eta_theta_levels(k-1))

            if( L_regular) then

              Do j = j_begin, j_end
                Do i = i_start, i_stop
                  term_upper(i,j) = 0.5 *                               &
     &                               ( interpx(i,j) + interpx(i-1,j) +  &
     &                                 interpy(i,j) + interpy(i,j-1) ) *&
     &                                                  HM_C5(i,j,k)
                  L_of_field(i,j,k) = L_of_field(i,j,k) -               &
     &                                 term_upper(i,j) * factor_1_p(k) *&
     &                                      FV_cos_theta_latitude(i,j)
                End do
              End do

            else  !  variable resolution

              Do j = j_begin, j_end
                Do i = i_start, i_stop
                  term_upper(i,j) = ( wt_phi_v(i,j) * interpy(i,j-1) +  &
     &                            (1.0 - wt_phi_v(i,j)) * interpy(i,j) +&
     &                                wt_lambda_u(i)  * interpx(i-1,j) +&
     &                         (1.0 - wt_lambda_u(i)) * interpx(i,j) ) *&
     &                                                    HM_C5(i,j,k)
                  L_of_field(i,j,k) = L_of_field(i,j,k) -               &
     &                                 term_upper(i,j) * factor_1_p(k) *&
     &                                      FV_cos_theta_latitude(i,j)
                End do
              End do

            endif !  L_regular

! save terms to perform calculation at boundaries
! Boundaries assumed flat in limited area model so code omitted.
            If (model_domain == mt_global) Then
! south pole
              If(at_extremity(PSouth))then
                If( L_regular ) then
                  Do i= 1, row_length
                    l_s_poles(i,k) = interpy(i,1)
                  End Do
                Else !  variable resolution
                  Do i= 1, row_length
                    l_s_poles(i,k) = interpy(i,1) * dlambda_p(i)
                  End Do
                Endif ! L_regular
              End If
! north pole
              If(at_extremity(PNorth))then
                If( L_regular ) then
                  Do i= 1, row_length
                    l_n_poles(i,k) = interpy(i,n_rows)
                  End Do
                Else !  variable resolution
                  Do i= 1, row_length
                    l_n_poles(i,k) = interpy(i,n_rows) * dlambda_p(i)
                  End Do
                Endif ! L_regular
              End If

            End If  ! model_domain == mt_global

          Else If ( k == first_constant_r_rho_level ) Then

! Calculate 1. / delta eta.
            factor_1_p(k) = 1. / (eta_theta_levels(k) -                 &
     &                            eta_theta_levels(k-1))

            Do j = j_begin, j_end
              Do i = i_start, i_stop
                L_of_field(i,j,k) = L_of_field(i,j,k) +                 &
     &                              term_upper(i,j)                     &
     &                           * factor_1_p(k)                        &
     &                           * FV_cos_theta_latitude(i,j)
              End Do
            End Do

          Else ! 1 < k < first_constant_r_rho_level
            if(ifirst==0)then
             ifirst=1
         interp_upper = (eta_theta_levels(k-1) - eta_rho_levels(k-1)) /  &
     &                     (eta_rho_levels(k) - eta_rho_levels(k-1))
            interp_lower = 1.0 - interp_upper

            Do j = j_begin, j_end
              Do i = i_start - 1, i_stop
                interpx(i,j) = ( interp_upper * x_term(i,j,k) +       &
     &                           interp_lower * x_term(i,j,k-1) ) *       &
     &                                          HM_Cxz(i,j,k-1)
              End Do
            End Do

           Do j = j_begin-1, j_end
              Do i = i_start, i_stop
                interpy(i,j) = ( interp_upper * y_term(i,j,k) +       &
     &                           interp_lower * y_term(i,j,k-1) ) *       &
     &                                          HM_Cyz(i,j,k-1)
              End Do
            End Do

            if( L_regular ) then

              Do j = j_begin, j_end
                Do i = i_start, i_stop
                  term_upper(i,j) = 0.5 * (interpx(i,j) + interpx(i-1,j)&
     &                                 + interpy(i,j) + interpy(i,j-1)) &
     &                                   * HM_C5(i,j,k-1)
                End Do
              End Do

            else  !  variable resolution

              Do j = j_begin, j_end
                Do i = i_start, i_stop
                  term_upper(i,j) = (  wt_phi_v(i,j) * interpy(i,j-1) + &
     &                           (1.0 - wt_phi_v(i,j)) * interpy(i,j) + &
     &                                wt_lambda_u(i) * interpx(i-1,j) + &
     &                        (1.0 - wt_lambda_u(i)) * interpx(i,j) ) * &
     &                                                   HM_C5(i,j,k-1)
                End Do
              End Do

            endif !  L_regular
         End if  !ifirst

! Interpolate quantities to levels either side in vertical
            interp_upper = (eta_theta_levels(k) - eta_rho_levels(k)) /  &
     &                     (eta_rho_levels(k+1) - eta_rho_levels(k))
            interp_lower = 1.0 - interp_upper

            Do j = j_begin, j_end
              Do i = i_start - 1, i_stop
                interpx(i,j) = ( interp_upper * x_term(i,j,k+1) +       &
     &                           interp_lower * x_term(i,j,k) ) *       &
     &                                          HM_Cxz(i,j,k)
              End Do
            End Do

           Do j = j_begin-1, j_end
              Do i = i_start, i_stop
                interpy(i,j) = ( interp_upper * y_term(i,j,k+1) +       &
     &                           interp_lower * y_term(i,j,k) ) *       &
     &                                          HM_Cyz(i,j,k)
              End Do
            End Do

! Calculate 1. / delta eta.
            factor_1_p(k) = 1. / (eta_theta_levels(k) -                 &
     &                            eta_theta_levels(k-1))


            if( L_regular ) then

              Do j = j_begin, j_end
                Do i = i_start, i_stop
                  term_lower = term_upper(i,j)
                  term_upper(i,j) = 0.5 * (interpx(i,j) + interpx(i-1,j)&
     &                                 + interpy(i,j) + interpy(i,j-1)) &
     &                                   * HM_C5(i,j,k)
                  L_of_field(i,j,k) = L_of_field(i,j,k) -               &
     &                              ( term_upper(i,j) - term_lower ) *  &
     &                                                 factor_1_p(k) *  &
     &                                    FV_cos_theta_latitude(i,j)
                End Do
              End Do

            else  !  variable resolution

              Do j = j_begin, j_end
                Do i = i_start, i_stop
                  term_lower = term_upper(i,j)
                  term_upper(i,j) = (  wt_phi_v(i,j) * interpy(i,j-1) + &
     &                           (1.0 - wt_phi_v(i,j)) * interpy(i,j) + &
     &                                wt_lambda_u(i) * interpx(i-1,j) + &
     &                        (1.0 - wt_lambda_u(i)) * interpx(i,j) ) * &
     &                                                   HM_C5(i,j,k)
                  L_of_field(i,j,k) = L_of_field(i,j,k) -               &
     &                              ( term_upper(i,j) - term_lower ) *  &
     &                                                 factor_1_p(k) *  &
     &                                    FV_cos_theta_latitude(i,j)
                End Do
              End Do

            endif !  L_regular

! save terms to perform calculation at boundaries
! Boundaries assumed flat in limited area model so code omitted.
            If (model_domain == mt_global) Then
! south pole
              If(at_extremity(PSouth))then
                If( L_regular ) then
                  Do i= 1, row_length
                    l_s_poles(i,k) = interpy(i,1)
                  End Do
                Else !  variable resolution
                  Do i= 1, row_length
                    l_s_poles(i,k) = interpy(i,1) * dlambda_p(i)
                  End Do
                Endif ! L_regular
              End If
! north pole
              If(at_extremity(PNorth))then
                If( L_regular ) then
                  Do i= 1, row_length
                    l_n_poles(i,k) = interpy(i,n_rows)
                  End Do
                Else !  variable resolution
                  Do i= 1, row_length
                    l_n_poles(i,k) = interpy(i,n_rows) * dlambda_p(i)
                  End Do
                Endif ! L_regular
              End If

            End If  !  model_domain == mt_global

          End If !  k == 1

        End Do  !  k = 1, first_constant_r_rho_level
!$OMP END DO
!$OMP END PARALLEL

! average the value and add on constant term, note any other polar
! point will do as all values are the same.
        If ( model_domain == mt_global) Then

          If ( at_extremity(PSouth) ) Then
            CALL global_2d_sums(l_s_poles, row_length, 1, 0, 0,         &
                                first_constant_r_rho_level-1,           &
                                sum_s, proc_row_group)

            Do k = 1, first_constant_r_rho_level

              If (k == 1) Then
                term_upper(1,1) = 2.0 * sum_s(k) * HM_C5(1,1,k) *       &
     &                                            recip_g_row_len
                L_of_field(1,1,k) = L_of_field(1,1,k) - term_upper(1,1)*&
     &                                                   factor_1_p(k) *&
     &                                      FV_cos_theta_latitude(1,1)
                Do i = 2, row_length
                  L_of_field(i,1,k) = L_of_field(1,1,k)
                End Do

              Else If (k == first_constant_r_rho_level) Then
                L_of_field(1,1,k) = L_of_field(1,1,k) + term_upper(1,1)*&
     &                                                   factor_1_p(k) *&
     &                                      FV_cos_theta_latitude(1,1)
              Do i = 2, row_length
                L_of_field(i,1,k) = L_of_field(1,1,k)
              End Do

            Else
              term_lower = term_upper(1,1)
              term_upper(1,1) = 2.0 * sum_s(k) * HM_C5(1,1,k) *         &
     &                                        recip_g_row_len
              L_of_field(1,1,k) = L_of_field(1,1,k) -                   &
     &                           ( term_upper(1,1) - term_lower ) *     &
     &                                              factor_1_p(k) *     &
     &                                 FV_cos_theta_latitude(1,1)
              Do i = 2, row_length
                L_of_field(i,1,k) = L_of_field(1,1,k)
              End Do

              End If  ! k == 1

            End Do ! k = 1, first_constant_r_rho_level

          End If ! at_extremity(PSouth)

          If (at_extremity(PNorth)) Then

            CALL global_2d_sums(l_n_poles, row_length, 1, 0, 0,         &
                                first_constant_r_rho_level-1,           &
                                sum_n, proc_row_group)

            j = rows
            Do k = 1, first_constant_r_rho_level

              If (k == 1) Then
                term_upper(1,j) = 2.0 * sum_n(k) * HM_C5(1,j,k) *       &
     &                                            recip_g_row_len
                L_of_field(1,j,k) = L_of_field(1,j,k) - term_upper(1,j)*&
     &                                                   factor_1_p(k) *&
     &                                      FV_cos_theta_latitude(1,j)
                Do i = 2, row_length
                  L_of_field(i,j,k) = L_of_field(1,j,k)
                End Do

              Else If (k  ==  first_constant_r_rho_level) Then

                L_of_field(1,j,k) = L_of_field(1,j,k) + term_upper(1,j)*&
     &                                                   factor_1_p(k) *&
     &                                      FV_cos_theta_latitude(1,j)
                Do i = 2, row_length
                  L_of_field(i,j,k) = L_of_field(1,j,k)
                End Do

              Else

                term_lower = term_upper(1,j)
                term_upper(1,j) = 2.0 * sum_n(k) * HM_C5(1,j,k) *       &
     &                                            recip_g_row_len
                L_of_field(1,j,k) = L_of_field(1,j,k) -                 &
     &                              ( term_upper(1,j) - term_lower ) *  &
     &                                                 factor_1_p(k) *  &
     &                                    FV_cos_theta_latitude(1,j)
                Do i = 2, row_length
                  L_of_field(i,j,k) = L_of_field(1,j,k)
                End Do

              End If !  k == 1

            End Do ! k = 1, first_constant_r_rho_level

          End If ! at_extremity(PNorth)

        End If ! model_domain == mt_global

      End If ! first_constant_r_rho_level > 1

      IF (lhook) CALL dr_hook('GCR_ELLIPTIC_OPERATOR_2B',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_Elliptic_Operator_2B
