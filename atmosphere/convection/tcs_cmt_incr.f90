! *****************************COPYRIGHT*******************************
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection

 SUBROUTINE TCS_CMT_INCR(n_npts, nlevs, ntop_max, kterm                &
                         ,r2rho, r2rho_th                              &
                         ,dr_across_rh                                 &
                         ,UW,VW                                        &
                                        ! output arguements
                                        !
                         ,dubydt,dvbydt)

! Description:
!
!  To calculate increments to U and V due to shallow convection.
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
   Integer, intent(in) :: &
     n_npts               & ! total number of convective columns
   , nlevs                & ! number of model levels                      
   , ntop_max             & ! maximum cloud top level
   , kterm(n_npts)          ! top of convection theta levels

   Real, intent(in) ::          &
     r2rho(n_npts,nlevs)        & ! r2 rho on rho levels (kg/m)
   , r2rho_th(n_npts,nlevs)     & ! r2 rho on theta levels (kg/m) 
   , dr_across_rh(n_npts,nlevs) & ! dr across rho levels (m)
   , UW(n_npts,nlevs)           & ! U-Component of stress profile (m2/s2)
   , VW(n_npts,nlevs)             ! V-Component of stress profile (m2/s2)

!
! Arguments with intent OUT:
!
   Real, intent(out) ::         &
     dubydt(n_npts,nlevs)       & ! Tendency in U (ms-2)
   , dvbydt(n_npts,nlevs)         ! Tendency in V (ms-2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

   Integer :: i,k      ! loop counters

   Real ::   rhodz     ! r2 rho * dz

   INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
   INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
   REAL(KIND=jprb)               :: zhook_handle


!
!-----------------------------------------------------------------------
! Input to this routine contains stress profile of uw & vw on levels
! from the surface to the top of the cloud.
!-----------------------------------------------------------------------
!
! Lowest uv level surface stress zero
!
   IF (lhook) CALL dr_hook('TCS_CMT_INCR',zhook_in,zhook_handle)
   k=1
   Do i=1,n_npts
     If (0 < kterm(i)) Then    ! Column has convection
       rhodz  = r2rho(i,k)*dr_across_rh(i,k)
       dubydt(i,K) = -uw(i,k)*r2rho_th(i,k)/rhodz
       dvbydt(i,K) = -vw(i,k)*r2rho_th(i,k)/rhodz
     End If    
   End do

!
! Mixed layer and Cloud layer increments 
! Note restricted to levels where increments apply as otherwise get 
! non-zero increments even though uw & vw = 0.0 due to numerics.
!
!CDIR NOUNROLL
   Do k=2,ntop_max+1    ! Should reduce cost where large stratosphere
!   Do k=2,nlevs   !
     Do i=1,n_npts
       If (k <= kterm(i)+1) Then   ! level where increments apply
         rhodz  = r2rho(i,k)*dr_across_rh(i,k)
         dubydt(i,k) =                                                 &
              -(r2rho_th(i,k)*uw(i,k)-r2rho_th(i,k-1)*uw(i,k-1))/rhodz
         dvbydt(i,k) =                                                 &
              -(r2rho_th(i,k)*vw(i,k)-r2rho_th(i,k-1)*vw(i,k-1))/rhodz
       End If    
     End Do
   End Do

   IF (lhook) CALL dr_hook('TCS_CMT_INCR',zhook_out,zhook_handle)
   RETURN

END SUBROUTINE TCS_CMT_INCR
