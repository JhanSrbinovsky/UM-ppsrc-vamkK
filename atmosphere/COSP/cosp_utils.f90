! (c) British Crown Copyright 2008-2013, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!  Description:  Module that contains COSP utilities: precipitation mixing ratio
!                from fluxes, quality control and initialisation routines,
!                and error/warning report.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP

MODULE MOD_COSP_UTILS
  USE cosp_constants_mod
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE

  INTERFACE Z_TO_DBZ
    MODULE PROCEDURE Z_TO_DBZ_2D,Z_TO_DBZ_3D,Z_TO_DBZ_4D
  END INTERFACE

  INTERFACE COSP_CHECK_INPUT
    MODULE PROCEDURE COSP_CHECK_INPUT_1D,COSP_CHECK_INPUT_2D,COSP_CHECK_INPUT_3D
  END INTERFACE
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_PRECIP_MXRATIO --------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_PRECIP_MXRATIO(Npoints,Nlevels,Ncolumns,p,T,prec_frac,         &
                  prec_type,n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,      &
                  gamma2,gamma3,gamma4,flux,mxratio,reff)

    ! Input arguments, (IN)
    INTEGER,INTENT(IN) :: Npoints,Nlevels,Ncolumns
    REAL,INTENT(IN),DIMENSION(Npoints,Nlevels) :: p,T,flux
    REAL,INTENT(IN),DIMENSION(Npoints,Ncolumns,Nlevels) :: prec_frac
    REAL,INTENT(IN) :: n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2,    &
                       gamma3,gamma4,prec_type
    ! Input arguments, (OUT)
    REAL,INTENT(OUT),DIMENSION(Npoints,Ncolumns,Nlevels) :: mxratio
    REAL,INTENT(INOUT),DIMENSION(Npoints,Ncolumns,Nlevels) :: reff
    ! Local variables
    INTEGER :: i,j,k
    REAL :: sigma,one_over_xip1,xi,rho0,rho,lambda_x,gamma_4_3_2,delta

    mxratio = 0.0

    ! N_ax is used to control which hydrometeors need to be computed
    IF (n_ax >= 0.0) THEN
      xi      = d_x/(alpha_x + b_x - n_bx + 1.0)
      rho0    = 1.29
      sigma   = (gamma2/(gamma1*c_x))*(n_ax*a_x*gamma2)**xi
      one_over_xip1 = 1.0/(xi + 1.0)
      gamma_4_3_2 = 0.5*gamma4/gamma3
      delta = (alpha_x + b_x + d_x - n_bx + 1.0)

      DO k=1,Nlevels
        DO j=1,Ncolumns
          DO i=1,Npoints
            IF ((prec_frac(i,j,k)==prec_type).OR.(prec_frac(i,j,k)==3.)) THEN
              rho = p(i,k)/(287.05*T(i,k))
              mxratio(i,j,k)=(flux(i,k)*((rho/rho0)**g_x)*sigma)**one_over_xip1
              mxratio(i,j,k)=mxratio(i,j,k)/rho
              ! Compute effective radius
              IF ((reff(i,j,k) <= 0.0).AND.(flux(i,k) /= 0.0)) THEN
                lambda_x =                                                     &
                   (a_x*c_x*((rho0/rho)**g_x)*n_ax*gamma1/flux(i,k))**(1./delta)
                reff(i,j,k) = gamma_4_3_2/lambda_x
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
END SUBROUTINE COSP_PRECIP_MXRATIO


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE ZERO_INT -------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELEMENTAL SUBROUTINE ZERO_INT(x,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10, &
                                 y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, &
                                 y21,y22,y23,y24,y25,y26,y27,y28,y29,y30)

  INTEGER,INTENT(INOUT) :: x
  INTEGER,INTENT(INOUT),OPTIONAL :: y01,y02,y03,y04,y05,y06,y07,y08,y09,y10, &
                                    y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, &
                                    y21,y22,y23,y24,y25,y26,y27,y28,y29,y30
  x = 0
  IF (PRESENT(y01)) y01 = 0
  IF (PRESENT(y02)) y02 = 0
  IF (PRESENT(y03)) y03 = 0
  IF (PRESENT(y04)) y04 = 0
  IF (PRESENT(y05)) y05 = 0
  IF (PRESENT(y06)) y06 = 0
  IF (PRESENT(y07)) y07 = 0
  IF (PRESENT(y08)) y08 = 0
  IF (PRESENT(y09)) y09 = 0
  IF (PRESENT(y10)) y10 = 0
  IF (PRESENT(y11)) y11 = 0
  IF (PRESENT(y12)) y12 = 0
  IF (PRESENT(y13)) y13 = 0
  IF (PRESENT(y14)) y14 = 0
  IF (PRESENT(y15)) y15 = 0
  IF (PRESENT(y16)) y16 = 0
  IF (PRESENT(y17)) y17 = 0
  IF (PRESENT(y18)) y18 = 0
  IF (PRESENT(y19)) y19 = 0
  IF (PRESENT(y20)) y20 = 0
  IF (PRESENT(y21)) y21 = 0
  IF (PRESENT(y22)) y22 = 0
  IF (PRESENT(y23)) y23 = 0
  IF (PRESENT(y24)) y24 = 0
  IF (PRESENT(y25)) y25 = 0
  IF (PRESENT(y26)) y26 = 0
  IF (PRESENT(y27)) y27 = 0
  IF (PRESENT(y28)) y28 = 0
  IF (PRESENT(y29)) y29 = 0
  IF (PRESENT(y30)) y30 = 0
END SUBROUTINE  ZERO_INT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE ZERO_REAL ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELEMENTAL SUBROUTINE ZERO_REAL(x,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10, &
                                 y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, &
                                 y21,y22,y23,y24,y25,y26,y27,y28,y29,y30)

  REAL,INTENT(INOUT) :: x
  REAL,INTENT(INOUT),OPTIONAL :: y01,y02,y03,y04,y05,y06,y07,y08,y09,y10, &
                                 y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, &
                                 y21,y22,y23,y24,y25,y26,y27,y28,y29,y30
  x = 0.0
  IF (PRESENT(y01)) y01 = 0.0
  IF (PRESENT(y02)) y02 = 0.0
  IF (PRESENT(y03)) y03 = 0.0
  IF (PRESENT(y04)) y04 = 0.0
  IF (PRESENT(y05)) y05 = 0.0
  IF (PRESENT(y06)) y06 = 0.0
  IF (PRESENT(y07)) y07 = 0.0
  IF (PRESENT(y08)) y08 = 0.0
  IF (PRESENT(y09)) y09 = 0.0
  IF (PRESENT(y10)) y10 = 0.0
  IF (PRESENT(y11)) y11 = 0.0
  IF (PRESENT(y12)) y12 = 0.0
  IF (PRESENT(y13)) y13 = 0.0
  IF (PRESENT(y14)) y14 = 0.0
  IF (PRESENT(y15)) y15 = 0.0
  IF (PRESENT(y16)) y16 = 0.0
  IF (PRESENT(y17)) y17 = 0.0
  IF (PRESENT(y18)) y18 = 0.0
  IF (PRESENT(y19)) y19 = 0.0
  IF (PRESENT(y20)) y20 = 0.0
  IF (PRESENT(y21)) y21 = 0.0
  IF (PRESENT(y22)) y22 = 0.0
  IF (PRESENT(y23)) y23 = 0.0
  IF (PRESENT(y24)) y24 = 0.0
  IF (PRESENT(y25)) y25 = 0.0
  IF (PRESENT(y26)) y26 = 0.0
  IF (PRESENT(y27)) y27 = 0.0
  IF (PRESENT(y28)) y28 = 0.0
  IF (PRESENT(y29)) y29 = 0.0
  IF (PRESENT(y30)) y30 = 0.0
END SUBROUTINE  ZERO_REAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE Z_TO_DBZ_2D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Z_TO_DBZ_2D(mdi,z)
    REAL,INTENT(IN) :: mdi
    REAL,DIMENSION(:,:),INTENT(INOUT) :: z
    ! Reflectivity Z:
    ! Input in [m3]
    ! Output in dBZ, with Z in [mm6 m-3]

    ! 1.e18 to convert from [m3] to [mm6 m-3]
    z = 1.e18*z
    WHERE (z > 1.0e-6) ! Limit to -60 dBZ
      z = 10.0*LOG10(z)
    ELSEWHERE
      z = mdi
    END WHERE
  END SUBROUTINE Z_TO_DBZ_2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE Z_TO_DBZ_3D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Z_TO_DBZ_3D(mdi,z)
    REAL,INTENT(IN) :: mdi
    REAL,DIMENSION(:,:,:),INTENT(INOUT) :: z
    ! Reflectivity Z:
    ! Input in [m3]
    ! Output in dBZ, with Z in [mm6 m-3]

    ! 1.e18 to convert from [m3] to [mm6 m-3]
    z = 1.e18*z
    WHERE (z > 1.0e-6) ! Limit to -60 dBZ
      z = 10.0*LOG10(z)
    ELSEWHERE
      z = mdi
    END WHERE
  END SUBROUTINE Z_TO_DBZ_3D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE Z_TO_DBZ_4D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Z_TO_DBZ_4D(mdi,z)
    REAL,INTENT(IN) :: mdi
    REAL,DIMENSION(:,:,:,:),INTENT(INOUT) :: z
    ! Reflectivity Z:
    ! Input in [m3]
    ! Output in dBZ, with Z in [mm6 m-3]

    ! 1.e18 to convert from [m3] to [mm6 m-3]
    z = 1.e18*z
    WHERE (z > 1.0e-6) ! Limit to -60 dBZ
      z = 10.0*LOG10(z)
    ELSEWHERE
      z = mdi
    END WHERE
  END SUBROUTINE Z_TO_DBZ_4D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_1D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_1D(vname,x,min_val,max_val)
    CHARACTER(LEN=*) :: vname
    REAL,INTENT(INOUT) :: x(:)
    REAL,INTENT(IN),OPTIONAL :: min_val,max_val
    LOGICAL :: l_min,l_max
    REAL :: min_v,max_v
    CHARACTER(LEN=128) :: pro_name='COSP_CHECK_INPUT_1D'

    l_min=.FALSE.
    l_max=.FALSE.

    IF (PRESENT(min_val)) THEN
      min_v = min_val
      IF (ANY(x < min_val)) THEN
      l_min = .TRUE.
        WHERE (x < min_val)
          x = min_val
        END WHERE
      ENDIF
    ENDIF
    IF (PRESENT(max_val)) THEN
      max_v = max_val
      IF (ANY(x > max_val)) THEN
        l_max = .TRUE.
        WHERE (x > max_val)
          x = max_val
        END WHERE
      ENDIF
    ENDIF

    CALL COSP_CHECK_INPUT_PRINT(vname,l_min,l_max,min_v,max_v)
  END SUBROUTINE COSP_CHECK_INPUT_1D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_2D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_2D(vname,x,min_val,max_val)
    CHARACTER(LEN=*) :: vname
    REAL,INTENT(INOUT) :: x(:,:)
    REAL,INTENT(IN),OPTIONAL :: min_val,max_val
    LOGICAL :: l_min,l_max
    REAL :: min_v,max_v
    CHARACTER(LEN=128) :: pro_name='COSP_CHECK_INPUT_2D'

    l_min=.FALSE.
    l_max=.FALSE.

    IF (PRESENT(min_val)) THEN
      min_v = min_val
      IF (ANY(x < min_val)) THEN
      l_min = .TRUE.
        WHERE (x < min_val)
          x = min_val
        END WHERE
      ENDIF
    ENDIF
    IF (PRESENT(max_val)) THEN
      max_v = max_val
      IF (ANY(x > max_val)) THEN
        l_max = .TRUE.
        WHERE (x > max_val)
          x = max_val
        END WHERE
      ENDIF
    ENDIF

    CALL COSP_CHECK_INPUT_PRINT(vname,l_min,l_max,min_v,max_v)
  END SUBROUTINE COSP_CHECK_INPUT_2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_3D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_3D(vname,x,min_val,max_val)
    CHARACTER(LEN=*) :: vname
    REAL,INTENT(INOUT) :: x(:,:,:)
    REAL,INTENT(IN),OPTIONAL :: min_val,max_val
    LOGICAL :: l_min,l_max
    REAL :: min_v,max_v
    CHARACTER(LEN=128) :: pro_name='COSP_CHECK_INPUT_3D'
    CHARACTER(LEN=100) :: cmessage

    l_min=.FALSE.
    l_max=.FALSE.

    IF (PRESENT(min_val)) THEN
      min_v = min_val
      IF (ANY(x < min_val)) THEN
      l_min = .TRUE.
        WHERE (x < min_val)
          x = min_val
        END WHERE
      ENDIF
    ENDIF
    IF (PRESENT(max_val)) THEN
      max_v = max_val
      IF (ANY(x > max_val)) THEN
        l_max = .TRUE.
        WHERE (x > max_val)
          x = max_val
        END WHERE
      ENDIF
    ENDIF

    CALL COSP_CHECK_INPUT_PRINT(vname,l_min,l_max,min_v,max_v)
  END SUBROUTINE COSP_CHECK_INPUT_3D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE COSP_CHECK_INPUT_PRINT ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_PRINT(vname,l_min,l_max,min_val,max_val)
    CHARACTER(LEN=*) :: vname
    LOGICAL,INTENT(IN) :: l_min,l_max
    REAL,INTENT(IN) :: min_val,max_val
    CHARACTER(LEN=*),PARAMETER :: pro_name='COSP_CHECK_INPUT_PRINT'
    CHARACTER(LEN=100) :: cmessage
    integer :: errcode

    IF (l_min) THEN
      WRITE(cmessage,*) 'Minimum value of '//TRIM(vname)//' set to: ',min_val
      errcode = -9
      CALL COSP_EREPORT(pro_name,cmessage,errcode)
    ENDIF
    IF (l_max) THEN
      WRITE(cmessage,*) 'Maximum value of '//TRIM(vname)//' set to: ',max_val
      errcode = -9
      CALL COSP_EREPORT(pro_name,cmessage,errcode)
    ENDIF
  END SUBROUTINE COSP_CHECK_INPUT_PRINT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE COSP_EREPORT ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_EREPORT(routine_name,message,errcode)
    CHARACTER(LEN = *), INTENT(IN) :: routine_name
    CHARACTER(LEN = *), INTENT(IN) :: message
    INTEGER, INTENT(INOUT) :: errcode

!  Subroutine that handles errors and warnings.
!
!   ErrorStatus > 0 is an error
!   ErrorStatus < 0 is a warning
!   ErrorStatus = 0 is ok

    CALL ereport(routine_name,errcode,message)

  END SUBROUTINE COSP_EREPORT


END MODULE MOD_COSP_UTILS
