! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE PR_INHDA---------------------------------------
!
!    Purpose: Prints out integer constants record and checks
!             validity of information.
!
!
!    Programming standard:  Unified Model Documentation Paper No 3
!
!    Documentation: Unified Model Documentation Paper No F3
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Dump I/O

SUBROUTINE pr_inhda                                               &
(inthd,len_inthd,row_length,p_rows                                &
,p_levels,q_levels,tr_levels                                      &
,st_levels,sm_levels,bl_levels                                    &
, tr_vars                                                        &
,icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER                                                           &
 len_inthd                                                        &
                !IN Length of integer header
,row_length                                                       &
                !IN No of points along a model row
,p_rows                                                           &
                !IN No of model rows
,p_levels                                                         &
                !IN No of model levels
,q_levels                                                         &
                !IN No of moist levels
,tr_levels                                                        &
                !IN No of tracer levels
,st_levels                                                        &
                !IN No of soil temperature levels
,sm_levels                                                        &
                !IN No of soil moisture levels
,bl_levels                                                        &
                !IN No of b. layer levels
,tr_vars        !IN No of tracer variables

INTEGER                                                           &
 inthd(len_inthd)                                                 &
                !IN Integer header
,icode          !OUT Return code; successful=0
                !                 error > 0

CHARACTER(LEN=80)                                                    &
 cmessage       !OUT Error message if ICODE > 0

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('PR_INHDA',zhook_in,zhook_handle)
icode=0
cmessage=' '

WRITE(6,'('' '')')
WRITE(6,'('' INTEGER CONSTANTS'')')
WRITE(6,'('' -----------------'')')
WRITE(6,'('' Number of timesteps since start of run -'',I7)')     &
inthd(1)
WRITE(6,'('' Meaning interval for mean fields (hrs) -'',I7)')     &
inthd(2)
WRITE(6,'('' Number of dumps used to generate mean  -'',I7)')     &
inthd(3)
WRITE(6,'(A,I7)')                                                 &
' No of hrs between neighbouring contiguous sections of means -'  &
,inthd(4)
WRITE(6,'(A,A,I7)')                                               &
' No of hrs between end of one contiguous section and',           &
' start of next -',inthd(5)

IF (row_length /= inthd(6) .OR. p_rows /= inthd(7)) THEN
  WRITE(6,'('' **FATAL ERROR** specifying model dimensions'')')
  WRITE(6,'('' Programmed dimensions ='',i7,'' x'',I7)')          &
        row_length,p_rows
  WRITE(6,'('' File dimensions ='',i7,'' x'',I7)')                &
        inthd(6),inthd(7)
  WRITE(6,'('' ***************'')')
  icode=4
  cmessage='PR_INHDA: Consistency check'
  IF (lhook) CALL dr_hook('PR_INHDA',zhook_out,zhook_handle)
  RETURN
ELSE
  WRITE(6,'('' Number of E-W x N-S points -'',i7,'' x'',I7)')     &
        inthd(6),inthd(7)
END IF

IF (p_levels /= inthd(8) .OR. q_levels /= inthd(9)) THEN
  WRITE(6,'('' **FATAL ERROR** specifying no of model levels'')')
  WRITE(6,'('' Programmed levels (wet) ='',i2,''('',i2,'')'')')   &
        p_levels,q_levels
  WRITE(6,'('' File levels (wet) ='',i2,''('',i2,'')'')')         &
        inthd(8),inthd(9)
  WRITE(6,'('' ***************'')')
  icode=4
  cmessage='PR_INHDA: Consistency check'
  IF (lhook) CALL dr_hook('PR_INHDA',zhook_out,zhook_handle)
  RETURN
ELSE
  WRITE(6,'('' Number of levels -'',I7)')                         &
        inthd(8)
  WRITE(6,'('' Number of wet levels -'',I7)')                     &
        inthd(9)
END IF
IF ((tr_levels /= inthd(12) .AND. inthd(14) /= 0)                 &
    .OR. st_levels /= inthd(10) .OR. sm_levels /= inthd(28)       &
    .OR. bl_levels /= inthd(13)) THEN
  WRITE(6,'('' **FATAL ERROR** specifying model level info'')')
  WRITE(6,'(A,A,4I3)')                                            &
  ' Programmed tracer, soil temperature,',                        &
  ' soil moisture,  b.l. levels =',                               &
          tr_levels,st_levels,sm_levels,bl_levels
  WRITE(6,'(A,A,4I3)')                                            &
  ' File tracer, soil temperature, soil moisture,' ,              &
  '  b.l. levels =',                                              &
          inthd(12),inthd(10),inthd(28),inthd(13)
  WRITE(6,'('' ***************'')')
  icode=4
  cmessage='PR_INHDA: Consistency check'
  IF (lhook) CALL dr_hook('PR_INHDA',zhook_out,zhook_handle)
  RETURN
ELSE
  WRITE(6,'('' Number of soil temperature levels -'',I7)')        &
        inthd(10)
  WRITE(6,'('' Number of soil moisture levels -'',I7)')           &
        inthd(28)
  WRITE(6,'('' Number of tracer levels -'',I7)')                  &
        inthd(12)
  WRITE(6,'('' Number of boundary layer levels -'',I7)')          &
        inthd(13)
END IF

IF (tr_vars /= inthd(14)) THEN
  WRITE(6,'('' **FATAL ERROR** specifying number of variables'')')
  WRITE (6,'('' Programmed tr_vars = '',I3)') tr_vars
  WRITE (6,'('' File       tr_vars = '',I3)') inthd(14)
  WRITE(6,'('' ***************'')')
  icode=4
  cmessage='PR_INHDA: Consistency check'
  IF (lhook) CALL dr_hook('PR_INHDA',zhook_out,zhook_handle)
  RETURN
ELSE
  WRITE(6,'('' Number of passive tracers advected -'',I7)')       &
        inthd(14)
END IF

IF (lhook) CALL dr_hook('PR_INHDA',zhook_out,zhook_handle)
RETURN
END SUBROUTINE pr_inhda
