! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE T_INT_C-------------------------------------------------
!
!    Purpose:
!      Carries out linear interpolation in time between two fields at
!      times T1 and T2. If the missing data indicator is present at one
!      of the times, the value at the other time is used. The interpolat
!      is controlled by a field ZI. A prescribed value is inserted where
!      If ZI changes between 0 and non-zero in the period T1 - T2, then
!      the field is linearly interpolated between its value at the time
!      ZI is non-zero and the prescibed value at the time when ZI become
!      The fractional time at which ZI changes between 0 and non-zero in
!      period T1 - T2 must be provided as input.
!
!
!    Programming standard: UMDP3 v8.2
!
!    System component:S191
!
!    System task: S1
!
!    Documentation:  The interpolation formulae are described in
!              unified model on-line documentation paper S1.
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Grids
MODULE t_int_c_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE t_int_c(data_t1,t1,data_t2,t2,data_t3,t3,points        &
,frac_time,zi_t1,pres_value)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE t_int_mod, ONLY: t_int
IMPLICIT NONE

INTEGER, INTENT(IN)  ::  points          !IN No of points to be processed

REAL,    INTENT(IN)  :: data_t1(points)  !IN Data at T1
REAL,    INTENT(IN)  :: data_t2(points)  !IN Data at T2
REAL,    INTENT(IN)  :: zi_t1(points)    !IN Value of controlling fieled at T1
REAL,    INTENT(IN)  :: pres_value(points) 
                                         !IN Prescribed value of Data when ZI=0
REAL,    INTENT(IN)  :: frac_time(points)
                                         !IN Fractional time at which ZI changes betwee
                                         !zero and non-zero in this time range

REAL,    INTENT(IN)  :: t1               !IN Time of first data field
REAL,    INTENT(IN)  :: t2               !In Time of second data field
REAL,    INTENT(IN)  :: t3    !IN Time at which new field is required T1<=T3<=T2

REAL,    INTENT(OUT)  :: data_t3(points)  !OUT Data at T3




! Local arrays:---------------------------------------------------------
 REAL  ::  int_time(points)

! Local variables:------------------------------------------------------
REAL   :: alpha         !Fractional time
REAL   :: epsi          ! rounding error
REAL   :: alpha_plus    ! add rounding error to alpha
REAL   :: alpha_minus   ! alpha minus rounding error

INTEGER  ::  i     !Loop index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

IF (lhook) CALL dr_hook('T_INT_C',zhook_in,zhook_handle)

! set rounding error
epsi=1.0e-6

CALL t_int(data_t1,t1,data_t2,t2,data_t3,t3,points)
 
alpha=(t3-t1)/(t2-t1)
alpha_plus=alpha+epsi
alpha_minus=alpha-epsi

DO i=1,points

  IF(frac_time(i) /= rmdi)THEN
    IF(zi_t1(i) == 0.0)THEN
       IF(alpha_minus <  frac_time(i))THEN
         data_t3(i)=pres_value(i)
       ELSE
         int_time(i)=(alpha-frac_time(i))/(1.-frac_time(i))
         data_t3(i)=pres_value(i)*(1.-int_time(i))                &
                   +data_t2(i)*int_time(i)
       END IF
    ELSE
       IF(alpha_plus >  frac_time(i))THEN
         data_t3(i)=pres_value(i)
       ELSE
         int_time(i)=(frac_time(i)-alpha)/(frac_time(i))
         data_t3(i)=pres_value(i)*(1.-int_time(i))                &
                   +data_t1(i)*int_time(i)
       END IF
    END IF
  END IF

END DO

IF (lhook) CALL dr_hook('T_INT_C',zhook_out,zhook_handle)
RETURN
END SUBROUTINE t_int_c
END MODULE t_int_c_mod
