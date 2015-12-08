! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module dealing with some common/simple manipulations of the pp header

! Description:
!   pp files can have different pp headers to the corresponding lookup
!   elements in a fields file.  The module provides routines for converting
!   the field file lookup table elements into pp headers.

!   See the documentation for the individual subroutines for more details

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs

! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.

MODULE pp_header_manips

  USE lookup_addresses

  IMPLICIT NONE

  INTEGER, PARAMETER :: max_len_ilabel=45

  ! Module Common variables
  LOGICAL  :: loper   = .FALSE.
  LOGICAL  :: lzero   = .FALSE.
  PRIVATE  :: loper, lzero
CONTAINS

  LOGICAL FUNCTION is_gregorian(ilabel)
    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  ilabel(max_len_ilabel)
    is_gregorian = (MOD(ilabel(lbtim), 10) == 1)
    RETURN
  END FUNCTION is_gregorian
  
  LOGICAL FUNCTION valid_start_year(ilabel)
    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  ilabel(max_len_ilabel)
    valid_start_year = (ilabel(lbyr) >  0)
    RETURN
  END FUNCTION valid_start_year
  
  SUBROUTINE set_oper(oper)
    IMPLICIT NONE
    
! Description:
!   set whether headers will follow the operational convention

    LOGICAL, INTENT(IN)    :: oper
    loper = oper
    
  END SUBROUTINE set_oper
  
  SUBROUTINE set_zero(zero)
    IMPLICIT NONE
    
    LOGICAL, INTENT(IN)    :: zero
    lzero = zero
    
  END SUBROUTINE set_zero

  SUBROUTINE header_manip(ilabel)  ! change names

    IMPLICIT NONE

! Description:
!   coordinates the calling of any pp-header changes that
!   need to be made
    INTEGER ilabel(max_len_ilabel)
    
    IF (loper) THEN
      CALL operheader(ilabel)
    END IF
    
    IF (lzero) THEN
      CALL zero_da_fields(ilabel)
    END IF
    
  END SUBROUTINE header_manip

  SUBROUTINE operheader(ilabel)

! Description:
!   Operational pp headers use a different convention to the
!   UM fields file headers
!   This routine converts from UM convention to operational convention.
!   This means:
!     for time mean fields or accumulations the order of the
!     time elements are swapped
!     for Met08 screen temperature the Met08 code is adjusted for max
!     and min fields.

    USE ereport_mod, ONLY : ereport
    IMPLICIT NONE
  
    INTEGER, PARAMETER  :: err_not_gregorian = 1
    INTEGER, PARAMETER  :: err_not_valid_start = 2
  
    INTEGER, INTENT(INOUT)  :: ilabel(max_len_ilabel)
  
    INTEGER                 :: &          
        start_year,            & 
        start_month,           & 
        start_day,             & 
        start_hour,            & 
        start_minute,          & 
        start_day_number,      & 
        start_second,          & 
        start_time_days,       & 
        start_time_secs,       & 
        end_year,              & 
        end_month,             & 
        end_day,               & 
        end_hour,              & 
        end_minute,            & 
        end_day_number,        & 
        end_second,            & 
        end_time_days,         & 
        end_time_secs,         & 
        data_year,             & 
        data_month,            & 
        data_day,              & 
        data_hour,             & 
        data_minute,           & 
        data_second,           & 
        data_day_number,       & 
        data_time_days,        & 
        data_time_secs,        & 
        fcst_prd,              & 
        time_period_secs
  
    INTEGER     :: lbtyp_old
    LOGICAL     :: lcal360
    LOGICAL     :: l_fixtyp


  ! Initialise time period to negative. This will make sure it will only be 
  ! used for times which are calculated when lbtime /= 11
    time_period_secs = -1
  ! Initialise lcal360 to false rather than initialising it with SAVE in
  ! declaration.
    lcal360 = .FALSE.
  
    IF (.NOT. is_gregorian(ilabel)) THEN
      CALL ereport('operheader', err_not_gregorian, &
                   'operheader assumes gregorian calendar')
    END IF
  
    IF (.NOT. valid_start_year(ilabel)) THEN
      CALL ereport('operheader', err_not_valid_start, &
                   'pp header start year not greater than 0')
    END IF
  
! For either an accumulation or time mean (ie LBTIM /= 0) the start &
! end time are in a different order to the data and veri time for a
! snap shot type field. This anomaly has to be catered for operational
! use. Thus the PP package will not work properly on accum/time mn field
! for operational Fields files.

    IF(ilabel(lbtim) /= 11) THEN

!       re -calculate the data time from the end time  and fcst period
!     First calculate the no of seconds from day 0
! if initial data time is missing, cannot recalculate data time
! for accum/time mean fields, hence, the PP Package would fail.
! Therefore, if the anomaly arises and temporal data is missing,
! this subroutine will be aborted and data from this field will
! not be written to COS file.
      start_year=ilabel(lbyr)
      start_month=ilabel(lbmon)
      start_day=ilabel(lbdat)
      start_hour=ilabel(lbhr)
      start_minute=ilabel(lbmin)
!
      end_year=ilabel(lbyrd)
      end_month=ilabel(lbmond)
      end_day=ilabel(lbdatd)
      end_hour=ilabel(lbhrd)
      end_minute=ilabel(lbmind)
! Support old (lbrel=2) day_number in pp header
!     and new (lbrel=3) second
      IF (ilabel(lbrel) >= 3) THEN
        start_second     = ilabel(lbsec)
        start_day_number = 0
        end_second       = ilabel(lbsecd)
        end_day_number   = 0
      ELSE
        start_second     = 0
        start_day_number = ilabel(lbsec)
        end_second       = 0
        end_day_number   = ilabel(lbsecd)
      END IF
      fcst_prd=ilabel(lbft)
! DEPENDS ON: time2sec
      CALL time2sec (end_year,end_month,end_day,end_hour,             &
                     end_minute,end_second,0,0,                       &
                     end_time_days,end_time_secs,lcal360)
! DEPENDS ON: time2sec
      CALL time2sec (start_year,start_month,start_day,start_hour,     &
                     start_minute,start_second,0,0,                   &
                     start_time_days,start_time_secs,lcal360)

      time_period_secs = (end_time_days-start_time_days)*24*3600 +    &
                         (end_time_secs-start_time_secs)

!   Subtract forecast hours from end time in (days,seconds)

! DEPENDS ON: time_df
      CALL time_df(end_time_days,end_time_secs,0,-fcst_prd*3600,      &
                   data_time_days,data_time_secs)

!     Go back and re-calculate Year/Month/Day/Hour/Sec.
! DEPENDS ON: sec2time
      CALL sec2time(0,0,data_time_days,data_time_secs,                &
                    data_year,data_month,data_day,                    &
                    data_hour,data_minute,data_second,data_day_number,&
                    lcal360)
      ilabel(lbyrd)=data_year
      ilabel(lbmond)=data_month
      ilabel(lbdatd)=data_day
      ilabel(lbhrd)=data_hour
      ilabel(lbmind)=data_minute
      IF (ilabel(lbrel) >= 3) THEN
        ilabel(lbsecd)=data_second
      ELSE
        ilabel(lbsecd)=data_day_number
      END IF
      ilabel(lbyr)=end_year
      ilabel(lbmon)=end_month
      ilabel(lbdat)=end_day
      ilabel(lbhr)=end_hour
      ilabel(lbmin)=end_minute
      IF (ilabel(lbrel) >= 3) THEN
        ilabel(lbsec)=end_second
      ELSE
        ilabel(lbsec)=end_day_number
      END IF
    END IF

    l_fixtyp = .FALSE.
    lbtyp_old = ilabel(lbtyp)

!       fix to correct max/min temp M08 codes
    IF(ilabel(lbtyp) == 58) THEN
      l_fixtyp = .TRUE.
!         check lbproc for max or min
      IF(ilabel(lbproc) == 4096) ilabel(lbtyp)=157  ! MIN
      IF(ilabel(lbproc) == 8192) ilabel(lbtyp)=156  ! MAX
    ! Fix to correct large-scale rain rates.
    ELSE IF (ilabel(lbtyp) == 63) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours
      ! (and 6 hours)
      IF(ilabel(lbproc) == 128 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
        ilabel(lbtyp) = 537
      END IF
    ! Fix to correct convective rain rates
    ELSE IF (ilabel(lbtyp) == 64) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours
      ! (and 6 hours)
      IF(ilabel(lbproc) == 128 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
         ilabel(lbtyp) = 533
      END IF
    ! Fix to correct large-scale snow rates
    ELSE IF (ilabel(lbtyp) == 68) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours 
      ! (and 6 hours)
      IF(ilabel(lbproc) == 128 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
         ilabel(lbtyp) = 538
      END IF
    ! Fix to correct convective snow rates
    ELSE IF (ilabel(lbtyp) == 69) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours 
      ! (and 6 hours)
      IF(ilabel(lbproc) == 128 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
         ilabel(lbtyp) = 534
      END IF
    ! Fix to correct high cloud amounts
    ELSE IF (ilabel(lbtyp) == 80) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours 
      ! (and 6 hours)
      IF(ilabel(lbproc) == 128 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
         ilabel(lbtyp) = 539
      END IF
    ! Fix to correct med. cloud amounts.
    ELSE IF (ilabel(lbtyp) == 81) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours 
      ! (and 6 hours)
      IF(ilabel(lbproc) == 128 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
         ilabel(lbtyp) = 540
      END IF
    ! Fix to correct low cloud amounts.
    ELSE IF (ilabel(lbtyp) == 82) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours 
      ! (and 6 hours)
      IF(ilabel(lbproc) == 128 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
         ilabel(lbtyp) = 541
      END IF
    ! Fix to correct lowest conv. cloud amount.
    ELSE IF (ilabel(lbtyp) == 201) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours 
      ! (and 6 hours)
      IF(ilabel(lbproc) == 128 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
         ilabel(lbtyp) = 542
      END IF
    ! Fix to correct lowest conv. cloud base
    ELSE IF (ilabel(lbtyp) == 202) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours 
      ! (and 6 hours)
      IF(ilabel(lbproc) == 8192 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
         ilabel(lbtyp) = 543
      END IF
    ! Fix to correct lowest conv. cloud top
    ELSE IF (ilabel(lbtyp) == 203) THEN
      l_fixtyp = .TRUE.
      ! Default uses 1 hour we want a special version for 3 hours 
      ! (and 6 hours)
      IF(ilabel(lbproc) == 4096 .AND. &
         (time_period_secs == 3*3600 .OR. time_period_secs == 6*3600)) THEN
         ilabel(lbtyp) = 544
      END IF
    END IF
    IF (l_fixtyp) THEN
      WRITE(6,'(A,I5,A,I11)') &
        'fix to TYPE=',lbtyp_old,' proc=',ilabel(lbproc)
      WRITE(6,'(A,I5)')' type=',ilabel(lbtyp)
    END IF

  END SUBROUTINE operheader

  SUBROUTINE zero_da_fields(ilabel)

    IMPLICIT NONE
! Description:
!   Some fields in the pp headers are only applicable to direct access data sets.
!   This routine sets those fields to zero.

    INTEGER, INTENT(INOUT)  :: ilabel(max_len_ilabel)

    ilabel(lbnrec) = 0
    ilabel(lbegin) = 0
    
  END SUBROUTINE zero_da_fields
  
END MODULE pp_header_manips
