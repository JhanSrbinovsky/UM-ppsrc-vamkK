! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+  Calculate the type of convection
!
!

SUBROUTINE dts_conv_classify(n_dp,klclm,klclp,kfrm,kfrp,k_adm,k_adp        &
                            ,ktopm,ktopp,zlcl,zfr,h_ad,ztop                &
                            ,iconvclass)



USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
! Simple classification routine to identify where the freezing level
! and level of maximum buoyancy are relative to the lifting
! condensation level and the convective top
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------
! Note at present while this routine is called the information it 
! generates is not actually used.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) ::     &
  n_dp                     &  ! Number of deep points
 ,klclm(n_dp)              &  ! Model level below LCL 
 ,klclp(n_dp)              &  ! Model level above LCL 
 ,kfrm(n_dp)               &  ! Model level below freezing level
 ,kfrp(n_dp)               &  ! Model level above freezing level
 ,k_adm(n_dp)              &  ! Model level below max buoyancy
 ,k_adp(n_dp)              &  ! Model level above max buoyancy
 ,ktopm(n_dp)              &  ! Model level below cloud top
 ,ktopp(n_dp)                 ! Model level above cloud top
        
REAL, INTENT(IN) ::        &
  zlcl(n_dp)               &  ! Height of the lifting condensation level (m)
 ,zfr(n_dp)                &  ! Height of the freezing level (m)
 ,h_ad(n_dp)               &  ! Height of maximum buoyancy excess (m)
 ,ztop(n_dp)                  ! Height of cloud top (m)
        
INTEGER, INTENT(OUT) ::    &
  iconvclass(n_dp)            ! Number indicating classification of convection

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  i_dp

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
  IF (lhook) CALL dr_hook('DTS_CONV_CLASSIFY',zhook_in,zhook_handle)
  iconvclass(:) = 0

  DO i_dp=1,n_dp

    ! Basic requirement that the top is on at least one level higher than 
    ! the lcl
    IF(ktopm(i_dp) > klclm(i_dp)) THEN

      ! standard tropical convection: require that there be a rho level
      ! between zlcl and zfr and between zfr and ztop
      IF(kfrm(i_dp) > klclm(i_dp) .AND. kfrm(i_dp) < ktopm(i_dp)) THEN
        iconvclass(i_dp) = 10
        ! maximum buoy beneath freezing level
        IF(k_adm(i_dp) <= kfrm(i_dp)) THEN
          iconvclass(i_dp) = 11
        END IF
        ! maximum buoy beneath lcl
        IF(k_adm(i_dp) <= klclm(i_dp)) THEN
          iconvclass(i_dp) = 12
        END IF 
      END IF

      ! cold convection -- freezing level below the lcl
      IF(kfrm(i_dp) <= klclm(i_dp)) THEN
        iconvclass(i_dp) = 20
        ! max buoy beneath lcl
        IF(k_adm(i_dp) <= klclm(i_dp)) THEN 
          iconvclass(i_dp) = 22
        END IF
        ! max buoy beneath freezing level < lcl
        ! (different order from above because kfrm < klclm)
        IF(k_adm(i_dp) <= kfrm(i_dp)) THEN 
          iconvclass(i_dp) = 21
        END IF
      END IF

      ! warm convection -- freezing level above the top of the layer
      IF(kfrm(i_dp) >= ktopm(i_dp)) THEN 
        iconvclass(i_dp) = 30
        ! max buoy beneath lcl
        IF(k_adm(i_dp) <= klclm(i_dp)) THEN
          iconvclass(i_dp) = 32
        END IF
      END IF
    ELSE 

    ! top below lcl ie dry convection, a job for the boundary layer?

      iconvclass(i_dp) = -1

    END IF
          
  END DO
  IF (lhook) CALL dr_hook('DTS_CONV_CLASSIFY',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE dts_conv_classify
