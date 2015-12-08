! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Shallow convection calculate Cloud base stress
!
      SUBROUTINE SHTCONV_BASE_STRESS(n_sh, nlevs, ntml, ntpar,ntpar_max &
     &,                              timestep                           &
     &,                              mb,wsc_o_mb                        &
     &,                              UW0,VW0,du_cb,dv_cb                &
     &,                              rho_theta,zrho,ztheta              &
     &,                              flg_uw_shall,flg_vw_shall          &
                                           ! IN/OUT ARGUMENTS
     &,                              UW,VW                              &
                                          ! OUTPUT ARGUMENTS
     &,                              uw_shall,vw_shall)

!
! Purpose:
!  To calculate Cloud base stress for shallow cumulus.
!  Also complete calcuation ofstress profile.
!
!   Called by Shallow_conv (5A version)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
      Integer, intent(in) ::                                            &
     &  n_sh                                                            &
                    ! Total number of shallow points
     &, nlevs                                                           &
                    ! Number of model levels
     &, ntml(n_sh)                                                      &
                    ! levels of LCL
     &, ntpar(n_sh)                                                     &
                     ! levels of TOP ofcloud layer
     &,ntpar_max     ! max value of ntpar +1
!
      Logical, intent(in) ::                                            &
     & flg_uw_shall                                                     &
                        ! STASH FLAGS FOR SHALLOW
     &,flg_vw_shall     ! CONVECTION STRESS DIAGNOSTIC
!
      Real, intent(in) ::                                               &
     &  timestep                                                        &
                         ! MODEL timestep (S)
     &,  mb(n_sh)                                                       &
                         ! Cloud base mass flux (MS-1)
     &,  wsc_o_mb(n_sh)                                                 &
                         ! Cloud-layer velocity scale (MS-1)
     &,  UW0(n_sh)                                                      &
                         ! U-component of surface stress (M2S-2)
     &,  VW0(n_sh)                                                      &
                         ! V-component of surface stress (M2S-2)
     &,  du_cb(n_sh)                                                    &
                         ! dU across cloud base (m/s)
     &,  dv_cb(n_sh)                                                    &
                         ! dV across cloud base (m/s)
     &,  rho_theta(n_sh,nlevs)                                          &
                               ! Density model th levels (kgm-3)
     &,  zrho(n_sh,nlevs)                                               &
                               ! height of model rho levels (m)
     &,  ztheta(n_sh,nlevs)    ! height of model theta levels (m)

!
! Arguments with intent INOUT:
!
      Real, intent(inout) ::                                            &
     &  uw(n_sh,nlevs)                                                  &
                         ! U-component of STRESS PROFILE (M2S-2)
     &, vw(n_sh,nlevs)   ! V-component of STRESS PROFILE (M2S-2)
!
! Arguments with intent OUT:
!
      Real, intent(out) ::                                              &
     &  uw_shall(n_sh,nlevs)                                            &
                              ! STASH DIAGNOSTIC FOR U-COMP STREss
     &, vw_shall(n_sh,nlevs)  ! STASH DIAGNOSTIC FOR V-COMP STREss
!
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

      Integer I,k            ! COUNTERS
!
      Real                                                              &
     &  omg2_jump(n_sh)                                                 &
                                ! jump in Y component of vorticity
     &, omg1_jump(n_sh)                                                 &
                                ! jump in X-component of vorticity
     &, zlcl_cmt(n_sh)                                                  &
                                ! lcl for CMT    (m)
     &, fcmt,DZ                                                         &
     &, A,B,C                                                           &
                                ! coefficients
     &, T,DZ1                                                           &
     &, z_depth(n_sh)                                                   &
     &, expadt         ! exp (Adt)

!
      Real ::                                                           &
     &  beta  = 0.04                                                    &
                         ! Coefficient for CMT calculation
     &, delta = 2.3                                                     &
                         ! Coefficient for CMT calculation
     &, gamma = 1.63     ! Coefficient for CMT calculation

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



!-----------------------------------------------------------------------
!
! Calculate jumps in vorticity across cloud base. (This is done by
! assuming that during the time step du and dv vary as exp(-T/TAU).
! Needs to be done to avoid instability around cloud base.)
!
      IF (lhook) CALL dr_hook('SHTCONV_BASE_STRESS',zhook_in,zhook_handle)
      Do i=1,n_sh

        zlcl_cmt(i) = ztheta(i,ntml(i))   ! cloud base fo CMT
!
! dz (cloud base - theta level below cloud base (zlcl_cmt) )
!
        dz = zrho(i,ntml(i)+1)-zlcl_cmt(i)
!
! dz across cloud base
!
        dz1 = ztheta(i,ntml(i)+1) - zlcl_cmt(i)
!
! depth of cloud as defined for uv calculations (different from th &q)
! calculations ?
!
        z_depth(i) = ztheta(i,ntpar(i)) - zlcl_cmt(i)

!  f(z/zcld) = exp(-fcmt)
        fcmt=beta*wsc_o_mb(i)*(dz1/z_depth(i))
        B=(1.0/zlcl_cmt(i)-(exp(-fcmt)-1.0)/dz1)
! alpha in doc (but extra /dz factor)
        A=zlcl_cmt(i)*mb(i)*B/(delta*dz)
        expadt = exp(-A*timestep)

! beta in documentation

        C = ( B*(1.0-gamma/delta) - 1.0/zlcl_cmt(i) )*uw0(i)

        if (C == 0.0.and.du_cb(i) == 0.0) then
          omg2_jump=0.0
        else
          T=-ALOG( (C*(1.0-expadt)/A+du_cb(i)*(expadt-1.0))/            &
     &                     ((C-A*du_cb(i))*timestep))/A
          omg2_jump(i)=( C*(1.0-exp(-A*T))/A + du_cb(i)*exp(-A*T) )/dz
        endif

! beta in documentation
        C = ( B*(1.0-gamma/delta) - 1.0/zlcl_cmt(i) )*vw0(i)

        if (C == 0.0.and.dv_cb(i) == 0.0) then
          omg1_jump=0.0
        else
          T=-ALOG((C*(1.0-expadt)/A+dv_cb(i)*(expadt-1.0))/             &
     &                     ((C-A*dv_cb(i))*timestep))/A
          omg1_jump(I)=-(C*(1.0-exp(-A*T))/A+dv_cb(i)*exp(-A*T))/dz
        endif
      End Do
!
! Calculate the cloud-base stress components
! Equations 11 & 12 section 5.1
!

      Do I=1,n_sh
        uw(i,ntml(i))=zlcl_cmt(i)*(-mb(i)*omg2_jump(i)-                 &
     &                 gamma*UW0(i)/zlcl_cmt(i))/delta+UW0(i)
        vw(i,ntml(i))=zlcl_cmt(i)*(mb(i)*omg1_jump(i)-                  &
     &                 gamma*VW0(i)/zlcl_cmt(i))/delta+VW0(i)
      End Do

!
! Calculate non-gradient stress profile in cloud
! Altered numbering of uw arrays
!
! New form of Fcmt   where alpha >1?  (like fng term for thermo?)
! Fcmt(z/zlcd)  = exp(-alpha*(z/zcld))

      Do k=1,ntpar_max

        Do i=1,n_sh


          If(k >= (ntml(i)+1).and.k <= (ntpar(i)-1)) then
! F function

            fcmt=exp(-1.1*(ztheta(i,k)-zlcl_cmt(i))/z_depth(i))

! all cloud levels add non-gradient term to eddy viscosity term

            uw(i,k)=uw(i,k)+uw(i,ntml(i))*fcmt
            vw(i,k)=vw(i,k)+vw(i,ntml(i))*fcmt

          Else if (k <= (ntml(i)-1)) then

! Which is correct ?
!            uw(i,k) = uw(i,ntml(i))*ztheta(i,k)/zlcl_cmt(i)
!            vw(i,k) = vw(i,ntml(i))*ztheta(i,k)/zlcl_cmt(i)
! This if assuming uw on rho levels
            uw(i,k) = uw(i,ntml(i))*zrho(i,k)/zlcl_cmt(i)
            vw(i,k) = vw(i,ntml(i))*zrho(i,k)/zlcl_cmt(i)

          End if

        End Do

      End Do

!
! Copy stress to output arrays multiplying by density on theta levels
!
      IF(flg_uw_shall) then
        Do k=1,ntpar_max+1
          Do i=1,n_sh
            uw_shall(i,k)=uw(i,k)*rho_theta(i,k)
          End Do
        End Do
      End If
      IF(flg_vw_shall) then
        Do k=1,ntpar_max+1
          Do i=1,n_sh
            vw_shall(i,k)=vw(i,k)*rho_theta(i,k)
          End Do
        End Do
      End If


      IF (lhook) CALL dr_hook('SHTCONV_BASE_STRESS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SHTCONV_BASE_STRESS
