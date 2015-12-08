! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_Coefficient_2B

      SUBROUTINE GCR_Coefficient_2B(                                    &
     &                           first_term, second_term,               &
     &                           row_length, rows, model_levels,        &
     &                           model_domain, coefficient,             &
     &                           at_extremity, n_proc,                  &
     &                           gc_proc_col_group, gc_proc_row_group,  &
     &                           l_datastart,                           &
     &                           i_start, i_stop, j_begin, j_end        &
     &                           )

! Purpose:
!          Calculates a coefficient given by
!
!                              < first_term, second_term >
!          inner_product = -   __________________________
!
!                              < second_term, second_term >
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
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
     &, model_levels                                                    &
                         ! number of model levels.
     &, i_start                                                         &
                      ! loop bound set in PE_Helmholtz
     &, i_stop                                                          &
                      ! loop bound set in PE_Helmholtz
     &, j_begin                                                         &
                      ! loop bound set in PE_Helmholtz
     &, j_end        ! loop bound set in PE_Helmholtz

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  first_term(row_length,rows,model_levels)                        &
                                                  ! first term in inner
                                                  ! product definition
     &, second_term(row_length,rows,model_levels)
                                                  ! second term in inner
                                                  ! product definition

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

! Arguments with Intent OUT. ie: variables Output only

      Real                                                              &
     &  coefficient

! Local Variables.

      Integer                                                           &
     &  i, j, k

      Integer                                                           &
     &  l_datastart(3)

      Integer                                                           &
     &  istat                                                           &
     &, gc_proc_col_group                                               &
     &, gc_proc_row_group                                               &
     &, n_proc

! parallel Local arrays
      Real                                                              &
     &  inner_product(2)                                                &
     &, term(row_length,rows,2)                                         &
     &, sum_temp(rows,2)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid




      Integer jj,j1,j2,kk,len1,len2
 

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('GCR_COEFFICIENT_2B',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 1.   Calculate Error Norm.
! ----------------------------------------------------------------------

! top term is in term (i,j,1)
! bottom term is in term (i,j,2)
      Do i= 1, row_length
        Do j= 1, rows
          term(i,j,1)=0.0
          term(i,j,2)=0.0
        End Do
      End Do

      If (model_domain == mt_global .and. l_datastart(1) == 1) Then
! Global model only calculate inner product for one of the polar points.

        If ( at_extremity(PSouth) ) then
          Do k = 1, model_levels
            term(1,1,1)= term(1,1,1) +                                  &
     &                   first_term(1,1,k) * second_term(1,1,k)
            term(1,1,2)= term(1,1,2) +                                  &
     &                   second_term(1,1,k) * second_term(1,1,k)
          End Do
        End If
        if (at_extremity(PNorth) ) then
          Do k = 1, model_levels
            term(1,rows,1)= term(1,rows,1) + first_term(1,rows,k)       &
     &                      * second_term(1,rows,k)
            term(1,rows,2)= term(1,rows,2) + second_term(1,rows,k)      &
     &                      * second_term(1,rows,k)
          End Do
        End If

      EndIf ! model_domain == mt_global .and. l_datastart(1) == 1

! The next sets of code execute the following algorithm, albeit in two
! optimised ways, one for vector strip-mining and the other for
! cache-blocking with OpenMP
!
!      Do k = 1, model_levels
!        Do j = j_begin, j_end
!          Do i = i_start, i_stop
!            term(i,j,1) = term(i,j,1) + first_term(i,j,k)               &
!     &                                * second_term(i,j,k)
!            term(i,j,2) = term(i,j,2) + second_term(i,j,k)              &
!     &                                * second_term(i,j,k)
!          End Do
!        End Do
!      End Do 

! This version is best for IBM. Cache-Blocked by 4 and hand unrolled
! by 2 with some OpenMP added.
      len1=4
      len2=(j_end -j_begin)/len1+1
!
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP& PRIVATE(k,j,i,jj,j1,j2,kk) SHARED(len1,len2,j_begin,j_end ,      &
!$OMP& model_levels,i_start, i_stop,term,first_term,second_term)
      Do jj=1, len1
        j1=j_begin+(jj-1)*len2
        j2=min(j_begin+jj*len2-1, j_end )
        Do k = 1, model_levels-1,2
          Do j = j1, j2
            Do i = i_start, i_stop
              term(i,j,1) = term(i,j,1) + first_term(i,j,k)             &
     &                                  * second_term(i,j,k)            &
     &                                  + first_term(i,j,k+1)           &
     &                                  * second_term(i,j,k+1)

              term(i,j,2) = term(i,j,2) + second_term(i,j,k)            &
     &                                  * second_term(i,j,k)            &
     &                                  + second_term(i,j,k+1)          &
     &                                  * second_term(i,j,k+1)
            End Do
          End Do
        End Do ! k = 1, model_levels-1,2
        If (mod(model_levels,2) /= 0) Then
          kk = model_levels
          Do j = j1, j2
            Do i = i_start, i_stop
              term(i,j,1) = term(i,j,1) + first_term(i,j,kk)            &
     &                                  * second_term(i,j,kk)
              term(i,j,2) = term(i,j,2) + second_term(i,j,kk)           &
     &                                  * second_term(i,j,kk)
            End Do
          End Do
        End If
      End Do !jj
!$OMP END PARALLEL DO

      CALL global_2d_sums(term, row_length, rows, 0, 0, 2,              &
                          inner_product)

! Calculate coefficient

      coefficient = - inner_product(1) / inner_product(2)

      IF (lhook) CALL dr_hook('GCR_COEFFICIENT_2B',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_Coefficient_2B
