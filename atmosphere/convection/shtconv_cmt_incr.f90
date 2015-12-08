! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Shallow convection calculate increments to u and v
!

      SUBROUTINE SHTCONV_CMT_INCR(n_sh, nlevs,ntpar_max                 &
     &,                            r2rho, r2rho_th                      &
     &,                            dr_across_rh                         &
     &,                            UW,VW                                &
                                        ! output arguements
     &,                            dubydt,dvbydt)
!
! Purpose:
!  To calculate increments to U and V due to shallow convection.
!
!   Called by Shallow_conv (5A version)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
      integer, intent(in) ::                                            &
     &  n_sh                                                            &
                       ! total number of shallow convcetive points
     &, nlevs                                                           &
                       ! number of model levels
     &, ntpar_max      ! maximum cloud top level

      real, intent(in) ::                                               &
     &  r2rho(n_sh,nlevs)                                               &
                                 ! r2 rho on rho levels
     &, r2rho_th(n_sh,nlevs)                                            &
                                 ! r2 rho on theta levels
     &, dr_across_rh(n_sh,nlevs)                                        &
                                 ! dr across rho levels
     &, UW(n_sh,nlevs)                                                  &
                           ! U-Component of stress profile (M2
     &, VW(n_sh,nlevs)     ! V-Component of stress profile (M2
!
! Arguments with intent OUT:
!
      real, intent(out) ::                                              &
     &  dubydt(n_sh,nlevs+1)                                            &
                                  ! Tendency in U (MS-2)
     &, dvbydt(n_sh,nlevs+1)      ! Tendency in V (MS-2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

      Integer :: I,K       ! counters
!
      Real ::                                                           &
     & rhodz               ! r2 rho * dz

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
!-----------------------------------------------------------------------
! Input to this routine contains stress profile of uw & vw on levels
! from the surface to the top of the inversion.
!-----------------------------------------------------------------------
!
! Lowest uv level surface stress zero
!
      IF (lhook) CALL dr_hook('SHTCONV_CMT_INCR',zhook_in,zhook_handle)
      k=1
      Do i=1,n_sh
        rhodz  = r2rho(i,k)*dr_across_rh(i,k)
        dubydt(i,K) = -uw(i,k)*r2rho_th(i,k)/rhodz
        dvbydt(i,K) = -vw(i,k)*r2rho_th(i,k)/rhodz
      End do

!
! Mixed layer and Cloud layer increments
!
      Do K=2,ntpar_max+1    !
        Do i=1,n_sh

          rhodz  = r2rho(i,k)*dr_across_rh(i,k)
          dubydt(i,K) =                                                 &
     &        -(r2rho_th(i,k)*uw(i,k)-r2rho_th(i,k-1)*uw(i,K-1))/rhodz
          dvbydt(i,K) =                                                 &
     &        -(r2rho_th(i,k)*vw(i,k)-r2rho_th(i,k-1)*vw(i,K-1))/rhodz
        End Do
      End Do

      IF (lhook) CALL dr_hook('SHTCONV_CMT_INCR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SHTCONV_CMT_INCR
