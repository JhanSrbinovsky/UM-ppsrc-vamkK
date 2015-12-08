! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  New Shallow convection scheme - calculation of d(w'u')/z
!
      SUBROUTINE SHTCONV_GRAD_STRESS(n_sh, nlev,max_cldlev              &
     &,                       nclev                                     &
     &,                       timestep, mb, wstar_up                    &
     &,                       zcld                                      &
     &,                       mf_cld, w_up, k_func                      &
     &,                       u,v                                       &
     &,                      r2rho,r2rho_th,dr_across_rh,dr_across_th   &
                                     ! output arguements
     &,                       uw_cld,vw_cld)
!
! Purpose:
!   Calculate the gradient component of the stress due to shallow
!   convection in the cloud
!
!   Called by shallow_conv.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
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
     &  nlev                                                            &
                 ! No. of model layers
     &, n_sh                                                            &
                 ! No. of shallow convection points
     &, nclev(n_sh)                                                     &
                     ! number of cloud levels
     &, max_cldlev   ! Maximum number of cloud levels

      real, intent(in) ::                                               &
     &  timestep                                                        &
                                 ! model timestep (s)
     &, zcld(n_sh)                                                      &
                                 ! cloud layer depth (m) for CMT
     &, mb(n_sh)                                                        &
                                 ! cloud base mass flux
     &, wstar_up(n_sh)                                                  &
                                 ! w* convective velocity scale (m/s)
     &, w_up(n_sh,nlev)                                                 &
                                 ! ensemble vertical velocity (m/s)
     &, mf_cld(n_sh,nlev)                                               &
                                 ! mass flux in (th lev) ensemble (m/s)
     &, k_func(n_sh,nlev)                                               &
                                 ! k function in cld
     &, u(n_sh,nlev)                                                    &
                                 ! U-component of mean wind (MS-1)
     &, v(n_sh,nlev)             ! V-component of mean wind (MS-1)
      real, intent(in) ::                                               &
     &  r2rho(n_sh,nlev)                                                &
     &, r2rho_th(n_sh,nlev)                                             &
     &, dr_across_th(n_sh,nlev)                                         &
     &, dr_across_rh(n_sh,nlev)

!
! Arguments with intent inout:
!
!
! Arguments with intent out:
!
      real, intent(out) ::                                              &
     &  uw_cld(n_sh,nlev)                                               &
                               ! U-component of stress
     &, vw_cld(n_sh,nlev)      ! V-component of stress

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
!
      Integer i,k    ! counters

      real                                                              &
     &  visc(n_sh,max_cldlev+1)                                         &
                                 ! viscosity profile (M2S-1)
     &, A(n_sh,max_cldlev+1)                                            &
                                 ! tridiagonal matrix elements
     &, B(n_sh,max_cldlev+1)                                            &
     &, C(n_sh,max_cldlev+1)                                            &
     &, ue_tp1(n_sh,max_cldlev+1)                                       &
                                   ! after timestep velocity vectors
     &, ve_tp1(n_sh,max_cldlev+1)                                       &
                                   ! for subsequent use
     &, dz,dz12
      Integer ::                                                        &
     &  nclevp1(n_sh)                                                   &
                           ! cloud level for uv cal
     &, max_cldlev1        ! max cloud levels plus 1

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      EXTERNAL TRIDIAG_ALL
!
!-----------------------------------------------------------------------
! Note 4A code assumed a 2 level inversion now go to 1 level
! ie drop top condition.
!-----------------------------------------------------------------------
!
!  wup and mass_flux input on required levels ?
!
! Calculate the eddy viscosity profile K(z/zcld)
!

      IF (lhook) CALL dr_hook('SHTCONV_GRAD_STRESS',zhook_in,zhook_handle)
      max_cldlev1 = max_cldlev+1

      Do k=1,max_cldlev1  ! from level above cloud base to just
        Do i=1,n_sh

          If (k <  nclev(i)) then         ! in cloud

            visc(i,k)=k_func(i,k)*mf_cld(i,k)                           &
     &                            *W_up(i,k)*zcld(i)/wstar_up(i)
          Else if(k == nclev(i)) then     ! inversion ?

            dz = dr_across_th(i,k)
            visc(i,k)=0.09*mb(i)*dz

          Else
            visc(i,k)=0.0
          Endif
        End Do
      End Do

!
! Calculate GRADIENT component of stress
!
!
! Use implicit timestepping
!

       k=1
       Do i=1,n_sh
             dz12 = dr_across_th(i,k)
             dz   = dr_across_rh(i,k)*r2rho(i,k)
             A(i,k)=0.0
             C(i,k)=-visc(i,k)*r2rho_th(i,k)                            &
     &                               *timestep/(dz*dz12)
             B(i,k)=1.0-A(i,k)-C(i,k)

           nclevp1(i) = nclev(i) + 1
       End do
       Do k=2,max_cldlev1
         Do i=1,n_sh
           dz = dr_across_rh(i,k)*r2rho(i,k)
           If (dz == 0.0) then
             write(6,*) ' dz = 0',dr_across_rh(i,k),r2rho(i,k),i,k
           endif
           If (k <= (nclev(i))) Then
             dz12 = dr_across_th(i,k-1)
           If (dz12 == 0.0) then
             write(6,*) ' dz12 = 0',dr_across_th(i,k-1),i,k
           endif
            A(i,k)=-visc(i,k-1)*r2rho_th(i,k-1)                         &
     &                               *timestep/(dz*dz12)
             dz12 = dr_across_th(i,k)
             C(i,k)=-visc(i,k)*r2rho_th(i,k)                            &
     &                               *timestep/(dz*dz12)
           Else if (k == (nclev(i)+1)) Then
             dz12 = dr_across_th(i,k-1)
             A(i,k)=-visc(i,k-1)*r2rho_th(i,k-1)                        &
     &                               *timestep/(dz*dz12)
             C(i,k)=0.0
           Else
! elements not required in calculation (zero)
             C(i,k) = 0.0
             A(i,k) = 0.0

           End if
           B(i,k)=1.0-A(i,k)-C(i,k)

         End do
       End do

!
! Calculate NEW timestep wind conponents using tridiag
!
! DEPENDS ON: tridiag_all
       CALL TRIDIAG_all(max_cldlev1,n_sh,nclevp1,A,B,C,u,ue_tp1)
! DEPENDS ON: tridiag_all
       CALL TRIDIAG_all(max_cldlev1,n_sh,nclevp1,A,B,C,v,ve_tp1)
!
! Calculate stress profile -Kdu/dz from latest u/v values
!

! Initialise uw,vw
      Do k=1,max_cldlev1
        Do i = 1,n_sh
           uw_cld(i,k) = 0.0
           vw_cld(i,k) = 0.0
        End Do
      End Do

      Do k=1,max_cldlev1
        Do i=1,n_sh
          If (k <= nclev(i)) then
            dz = dr_across_th(i,k)
            uw_cld(i,k)=-visc(i,k)*(ue_tp1(i,k+1)-ue_tp1(i,k))/dz
            vw_cld(i,k)=-visc(i,k)*(ve_tp1(i,k+1)-ve_tp1(i,k))/dz
          End if
        End Do
      End Do

!-----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('SHTCONV_GRAD_STRESS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SHTCONV_GRAD_STRESS
