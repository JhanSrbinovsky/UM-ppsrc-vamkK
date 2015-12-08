! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generic routines for manipulating pp-fields within Fieldcalc

!=======================================================================

!=======================================================================
SUBROUTINE sstpert( factor, Fieldin, nens, fieldclim,   &  ! in
                    Fieldout, ErrorStatus )                ! inout

! Description:
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:                                                       &
  PP_Header_type,                                                       &
  PP_Field_type
USE Err_Mod, ONLY:                                                      &
  StatusOK,                                                             &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN)    :: fieldin
TYPE(PP_Field_type), INTENT(IN)    :: fieldclim(12)
                       ! Monthly DeltaSST climatology from ancil file
INTEGER, INTENT(IN) :: nens
                       ! Ensemble member number

TYPE(PP_Field_type), INTENT(INOUT) :: fieldout
INTEGER, INTENT(INOUT) :: ErrorStatus
REAL, INTENT(IN) :: factor
                       ! SST perturbation inflation factor (alpha)
INTEGER :: dt(8),                                                       &
                       ! date array
           iarg,                                                        &
                       ! used for random seed
           tam,                                                         &
                       ! used for random seed
           i_month1, i_month2
                       ! months adjacent to current date
INTEGER :: i,j
INTEGER, ALLOCATABLE ::                                                 &
           prevseed(:),                                                 &
                       ! random seed (current)
           iranseed(:) ! random seed (next)
REAL :: r_month1, r_month2
                       ! months adjacent to current date
REAL, ALLOCATABLE ::                                                    &
        rnum(:)
REAL, ALLOCATABLE ::                                                    &
        psif(:,:)

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "sstpert"
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

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

!-----------------------------------------------------
! Step 1: Setup random seed using field valid date
!-----------------------------------------------------
 dt(1) = fieldin % hdr % validyear
 dt(2) = fieldin % hdr % validmonth
 dt(3) = fieldin % hdr % validdate
 dt(4) = 0                ! Shift from UTC
 dt(5) = fieldin % hdr % validhour + 1
                          ! Set range 1 - 24
 dt(6) = fieldin % hdr % validmin
 dt(7) = nens             ! Ens mem 0-24 into second dimension
 dt(8) = nens + 100       ! Ens mem 100-124 into millisecond dim

 CALL random_seed(SIZE=tam)
 IF (.NOT.ALLOCATED(prevseed)) ALLOCATE(prevseed(tam))
 IF (.NOT.ALLOCATED(iranseed)) ALLOCATE(iranseed(tam))
 IF (.NOT.ALLOCATED(rnum)) ALLOCATE(rnum(tam))

! Use same method as in Atmos_Section35
 dt(1) = dt(1) - 100*INT(0.01*dt(1))
 iarg = (dt(3) - 32075 +                                            &
        1461*(dt(1) + 4800 + (dt(2) - 14)/12)/4 +                 &
        367*(dt(2) - 2 - (dt(2)-14)/12*12)/12 -                   &
        3*((dt(1)+4900+(dt(2)-14)/12)/100)/4)*1000 +              &
        dt(8)**2.86 + dt(7)**3.79 + dt(5)**5.12 +                 &
        dt(6)**3.24
! Generate initial random number set to use as power
 prevseed(:) = iarg
 CALL random_seed(PUT=prevseed(1:tam))
 CALL random_number(rnum)
! Range of seed from 0 to 2**31 (32-bit Int)
 iranseed(:)=iarg*rnum(:)
! Set final seed
 CALL random_seed(PUT=iranseed(1:tam))

 WRITE(6,*) 'SSTPERT SETUP: Sizeof Random Seed, Ensemble Member'
 WRITE(6,*) tam, nens
 WRITE(6,*) 'SSTPERT SETUP: model date = ',dt
 WRITE(6,*) 'SSTPERT SETUP: Random seed generated using model date'
 WRITE(6,*) iranseed

 DEALLOCATE(prevseed)
 DEALLOCATE(iranseed)
 DEALLOCATE(rnum)


! Initialise output field
 IF ( ASSOCIATED( fieldout % RData ) ) THEN
   DEALLOCATE( fieldout % RData )
 END IF
 fieldout % Hdr = fieldin % Hdr
 ALLOCATE( fieldout % RData(fieldout % Hdr % NumCols,                   &
                            fieldout % Hdr % NumRows) )

! Find the two adjacent months nearest AINITIAL date
!   and calculate ratio of each month
 IF ( dt(3) < 15 ) THEN
    i_month1 = dt(2) - 1
    r_month1 = (16 - dt(3))/30.
    r_month2 = 1. - r_month1
 ELSE
    i_month1 = dt(2)
    r_month2 = (dt(3) - 15)/30.
    r_month1 = 1. - r_month2
 ENDIF
 i_month2 = i_month1 + 1
 IF ( i_month2 > 12 ) i_month2 = i_month2 - 12
 IF ( i_month1 < 1 ) i_month1 = i_month1 + 12

! Generate Random field with specified power spectrum
 IF (.NOT.ALLOCATED( psif ))                                            &
 ALLOCATE( psif(fieldout % Hdr % NumCols, fieldout % Hdr % NumRows) )
! DEPENDS ON: sst_genpatt
 CALL sst_genpatt(fieldout % Hdr % NumCols,                             &
                  fieldout % Hdr % NumRows, psif)

! Find absolute sqrt of forcing pattern (power is in units K^2)
 DO j = 1, fieldout % Hdr % NumRows
   DO i = 1, fieldout % Hdr % NumCols
     psif(i,j) = SIGN( SQRT( ABS( psif(i,j) )), psif(i,j))
   END DO
 END DO

!  multiply by interpolated DSST field and add to surface temp
 DO j = 1, fieldout % Hdr % NumRows
   DO i = 1, fieldout % Hdr % NumCols
     fieldout % RData(i,j) = psif(i,j) * factor *                       &
                    (r_month1 * fieldclim(i_month1) % RData(i,j) +      &
                     r_month2 * fieldclim(i_month2) % RData(i,j))
   END DO
 END DO

 DEALLOCATE(psif)

9999 CONTINUE

END SUBROUTINE sstpert

