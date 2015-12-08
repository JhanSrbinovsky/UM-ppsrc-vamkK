! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Adds input time interval to input time
!

      subroutine calc_new_time(a_intf_freq_hr,                          &
                          a_intf_freq_mn,a_intf_freq_sc,                &
                          yy,mm,dd,hr,mn,ss,day_no,lcal360)

      implicit none

! Description:
!   Subroutine takes current time and an interval as input, and returns
!   the time with the interval added
!
! Method:
!   Calls time2sec to convert input time to days and seconds
!   Converts interval time into days and seconds
!   Calls sec2time to add the interval time in days and seconds
!   to the input time in days and seconds, returning the time in
!   yy,mm,dd,hr,mn,ss,day_no,lcal360 format
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

! Scalar arguments with intent (in)
      integer, intent(in) :: a_intf_freq_hr
      integer, intent(in) :: a_intf_freq_mn
      integer, intent(in) :: a_intf_freq_sc

! Scalar arguments with intent (ionout)
      integer, intent(inout) :: yy,mm,dd,hr,mn,ss,day_no

! Logical arguments with intent (inout)
      logical, intent(inout) :: lcal360

! Local variables
      integer :: basis_time_days
      integer :: basis_time_secs
      integer :: elapsed_days
      integer :: elapsed_secs
      integer :: previous_time_days
      integer :: previous_time_secs
      integer :: days
      integer :: interval_secs


      basis_time_days = 0
      basis_time_secs = 0

! End of header

! Call time2sec on the previous time to get a time in days and seconds
! relative ti the basis_time
! DEPENDS ON: time2sec
      call time2sec(yy,mm,dd,hr,mn,ss,                                  &
             basis_time_days,basis_time_secs,elapsed_days,elapsed_secs, &
             lcal360)


! Set previous days and seconds to value given by time2sec i.e. time
! from basis_time
      previous_time_days = elapsed_days
      previous_time_secs = elapsed_secs


! Calculate number of seconds between lbc times from variables in
! intfcnsta namelist
      interval_secs=(a_intf_freq_hr*3600)+(a_intf_freq_mn*60)+          &
                     a_intf_freq_sc

! Break the number of seconds down into days and seconds
      days=0
      do
        if(interval_secs > 86400)then
          interval_secs=interval_secs-86400
          days=days+1
        else
          exit
        endif
      enddo

! Call sec2time to add the time between the lbcs to the
! time since the basis time
! DEPENDS ON: sec2time
      call sec2time(days,interval_secs,previous_time_days,              &
                    previous_time_secs,                                 &
                    yy,mm,dd,hr,mn,ss,day_no,lcal360)

      return
      END SUBROUTINE calc_new_time
