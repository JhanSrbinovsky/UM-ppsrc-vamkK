! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!      Subroutine PRINTSUB
!
!     Single Column Unified Model routine to write out T and Q, and
!     increments DT and DQ due to statistical forcing.
!     Will print out any number of columns
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!=====================================================================

SUBROUTINE printsub (points, nlevs, nwet, title, istep, iday, t, q, dt, dq)

  IMPLICIT NONE

!---------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------

  INTEGER ::         &
    points           &! IN no of model columns.
  , nlevs            &! IN no of levels.
  , nwet              ! IN no of model levels in which Q is set

  CHARACTER(LEN=60) ::    &
    title             ! Heading

  INTEGER ::         &
    iday             &! Day number
  , istep             ! Timestep

  REAL ::            &
    q(points,nlevs)  &! Specific humidity (kg/kg)
  , t(points,nlevs)  &! Temperature (K)
  , dt(points,nlevs) &! Temperature increment(K)
  , dq(points,nlevs)  ! Specific humidity increment (kg/kg)

!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------

  CHARACTER(LEN=28) ::    &
    cfmt              ! Format statement for each row of variables Q T DT DQ

  CHARACTER(LEN=33) ::    &
    ctfmt             ! Format statement for title of each row

  INTEGER ::         &
    i                &! Write statement loop counter
  , l                &! Loop counter
  , element          &! Array element no.
  , lastrow          &! No. of elements in last row
  , nlevsrows        &! No. of rows and Do Loop counter
  , nlevscount       &
  , nwetrows         &
  , nwetcount

! Set format statements

  cfmt = '(''         '',  (1pe10.3,1x))'
  ctfmt = '(''0       '',  (3x,''Level'',i2,1x))'

! Loop over sites :

  DO l=1, points

! Heading

    WRITE (11,'(A32,I5,A11,I5,A60)')            &
      '0Day relative to winter solstice', iday, &
      '; timestep ', istep, title

! Write out variables T and DT for NLEVS, maximum of 10
! variables per row

    WRITE (11,'(A33,I3)') '0 Number of atmospheric levels = ', nlevs

! Calculate no. of rows and no. of elements in last row

    IF (MOD(nlevs,10) == 0) THEN
      nlevsrows = INT(nlevs/10)
      lastrow = 10
    ELSE
      nlevsrows = INT(nlevs/10) + 1
      lastrow = MOD(nlevs,10)
    END IF

    DO nlevscount=1, nlevsrows
      element = 10*(nlevscount-1)

      IF (nlevscount  <   nlevsrows) THEN

! Write out all complete rows ie of 10 variables per row

        WRITE(11,'(A16,i2,9(4x,"level",i2))')                             &
          '0          level', (element+i, i=1, 10)

        WRITE(11,'(A9,10(1pe10.3,1x)/,A9,10(1pe10.3,1x)/)')               &
          ' t k     ',(t(l,element+i), i=1,10),                           &
          ' dT K    ',(dt(l,element+i), i=1, 10)
      ELSE

! Write out last row. Use an internal format statement by creating a
! character string. This will enable a variable format to be created eg
! NF10.6 where N is the no. of elements in the last row which can be
! written into the format statement via an internal write statement.

        WRITE (ctfmt(13:14),'(i2)') lastrow
        WRITE (11,ctfmt) (element+i,i=1,lastrow)
        WRITE (cfmt(14:15),'(i2)') lastrow
        WRITE (cfmt(4:7),'(''T K '')')
        WRITE (11,cfmt) (t(l,i+element),i=1,lastrow)
        WRITE (cfmt(4:7),'(''dT K'')')
        WRITE (11,cfmt) (dt(l,i+element),i=1,lastrow)
      END IF
    END DO

! Write out variables Q and DQ for NWET, maximum of 10
! variables per row

    WRITE (11,'(A39,I3)')                                                 &
      '0 Number of moist atmospheric levels = ', nwet

! Calculate no. of rows and no. of elements in last row

    IF ( MOD(nwet,10)  ==  0) THEN
      nwetrows = INT(nwet/10)
      lastrow = 10
    ELSE
      nwetrows = INT(nwet/10) + 1
      lastrow = MOD(nwet,10)
    END IF

    DO nwetcount=1, nwetrows
      element = 10*(nwetcount-1)

      IF (nwetcount  <   nwetrows) THEN

! Write out all complete rows ie of 10 variables per row

        WRITE (11,'(A16,i2,9(4x,"level",i2))')                            &
          '0          level', (element+i, i=1, 10)

        WRITE (11,'(A9,10(1pe10.3,1x)/,A9,10(1pe10.3,1x)/)')              &
          ' Q Kg/Kg ', (q(l,element+i),  i=1, 10),                        &
          ' dQ Kg/Kg', (dq(l,element+i), i=1, 10)

      ELSE

! Write out last row. Use an internal format statement by creating a
! character string. This will enable a variable format to be created eg
! NF10.6 where N is the no. of elements in the last row which can be
! written into the format statement via an internal write statement.

        WRITE (ctfmt(13:14),'(i2)') lastrow
        WRITE (11,ctfmt) (element+i,i=1,lastrow)
        WRITE (cfmt(4:10),'(''Q Kg/Kg'')')
        WRITE (11,cfmt) (q(l,i+element),i=1,lastrow)
        WRITE (cfmt(4:11),'(''dQ Kg/Kg'')')
        WRITE (11,cfmt) (dq(l,i+element),i=1,lastrow)
      END IF

    END DO
  END DO                     ! l

  RETURN

END SUBROUTINE printsub

