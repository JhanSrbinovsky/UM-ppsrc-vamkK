! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE StatStep
! Purpose:-           To calculate statistical forcing required each
!                     timestep
! Programmer:-        J. LEAN - modified code from original SCM to
!                     meet UM standards + now includes interpolation
!                     of all forcing variables between subsequent days
!                     Last modified on 1/3/93 to include double-
!                     precision versions of NAG random number
!                     generator routines (ie G05) for compatability
!                     with workstation
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!=====================================================================

SUBROUTINE statstep                                                           &
  ( row_length, rows, nlevs, nwet, ntrop, deltan, px, py, daysteps            &
  , stepcount, dayno, tr, vnr, vpr, qr, wr, tbar, tsd, tdash, dbar, dsd       &
  , ddash, vnbar, vpbar, vnsd, wbar, wsd, ctbar, ctsd, at, cdbar, cdsd, ad    &
  , cvnbar, cvnsd, avn, cwbar, cwsd, aw, press, rpress, u, v, w, t, q         &
  , prinstat, u_inc, v_inc, w_inc, t_inc, q_inc, daycount, timestep, iv, ntab &
  , iy, idum )

  IMPLICIT NONE

!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------

  INTEGER ::           &
    row_length         &! In no of model rows/columns
  , rows               &
  , nlevs              &! In no of levels.
  , nwet               &! In no of model levels in which Q is set
  , ntrop              &! In Max number of levels in the troposphere
  , daycount           &! In Daynumber (1 represents 1st january)
  , dayno              &! In Daynumber relative to winter solstice
  , daysteps           &! In No. of timesteps in 1 day
  , stepcount          &! In Timestep counter
  , ntab               &! In Dimension of array used in random generator
  , iv(ntab), iy, idum  ! In state of number generator saved
                        !    on tape for continuation run

  LOGICAL ::           &
    prinstat            ! T if printout of increments due to statistical
                        !   forcing required each timestep

  REAL ::                          &
    ad(row_length,rows,nwet-1)     &! In     term a of equation 2.22 for dew
  , at(row_length,rows,nlevs-1)    &!        pt depression and temp.
  , avn(row_length,rows,nlevs-1)   &! In     term a of equation 2.22 for
  , aw(row_length,rows,ntrop-1)    &!        horiz. and vertical velocities
  , cdbar(row_length,rows,nwet)    &!        Mean and SD of random variable
  , cdsd(row_length,rows,nwet)     &!        for dew point depression
  , ctbar(row_length,rows,nlevs)   &! In     Mean and SD of random variable
  , ctsd(row_length,rows,nlevs)    &!        for temp.

  , cvnbar(row_length,rows,nlevs)  &! In    Mean and SD of random variable
  , cvnsd(row_length,rows,nlevs)   &!       for velocity VN
  , cwbar(row_length,rows,ntrop)   &
  , cwsd(row_length,rows,ntrop)    &! In    Mean and SD of random variable
                                    !       for vertical velocity
  , dbar(row_length,rows,nwet)     &! In    Mean and SD dewpoint
  , dsd(row_length,rows,nwet)      &!       depression at daycount days
                                    !       from winter solstice (K)
  , ddash(row_length,rows,nwet)    &! In    Dew pt. corrections
  , deltan(row_length,rows)        &! In    Radius of area (m)
  , dq1                            &! In    Spec. humidity differences
  , dq2                            &!       (kg/kg)
  , dt1                            &! In    Temp. differences (K)
  , dt2                            &
  , press(row_length,rows,nlevs)   &! In    Pressure coordinates (Pa)
  , px(row_length,rows,ntrop)      &! In    Reciprocal log functions for
  , py(row_length,rows,ntrop-1)    &!       calc. of vert. advection
                                    !       used in eqns 2.12 and 2.13
  , q(row_length,rows,nwet)        &! Out   Specific humidity (kg/kg)
  , qr(row_length,rows,nwet,2)     &! InOut Randomly sampled humidity
                                    !       (kg/kg)
  , rpress(row_length,rows,nlevs)  &! In    Reciprocal pressure for rho
                                    !       levels (1/HPa or 1/mb)
  , t(row_length,rows,nlevs)       &! Out   Temp (K)
  , tbar(row_length,rows,nlevs)    &! In    Mean and SD temperature at
                                    !       daycount days from
                                    !       winter solstice (K)
  , tdash(row_length,rows,nlevs)   &! In    Temp. corrections (K)
  , timestep                       &! In    model timestep (s)
  , tsd(row_length,rows,nlevs)     &! In    SD of temp. at daycount days
                                    !       from winter solstice (K)
  , tr(row_length,rows,nlevs,2)    &! InOut Randomly sampled temp. (K)
  , u(row_length,rows,nlevs)       &
  , v(row_length,rows,nlevs)       &! Out   Zonal and meridional winds (m/s)
  , vnbar(row_length,rows,nlevs)   &! In    Mean and SD velocity VN at
                                    !       daycount days from winter
                                    !       solstice (m/s)
  , vnr(row_length,rows,nlevs,2)   &! InOut Randomly sampled horizontal
                                    !       velocity (m/s)
  , vnsd(row_length,rows,nlevs)    &! In    Mean and SD velocity VN at
                                    !       daycount days from
                                    !       winter solstice (m/s)
  , vpbar(row_length,rows,nlevs)   &! In    Mean velocity VP at
                                    !       daycount days from
                                    !       winter solstice (m/s)
  , vpr(row_length,rows,nlevs,2)   &! InOut Randomly sampled horizontal
                                    !       velocity (m/s)
  , w(row_length,rows,0:nlevs)     &! Vertical velocity
  , wbar(row_length,rows,ntrop)    &
  , wsd(row_length,rows,ntrop)     &! In    Mean and SD vertical
                                    !       velocity at daycount days
                                    !       from winter solstice (mb/s)
  , wr(row_length,rows,ntrop,2)    &! InOut Randomly sampled vertical
                                    !       velocity (mb/s)
  , u_inc(row_length,rows,nlevs)   &! Out   Zonal and meridional wind
  , v_inc(row_length,rows,nlevs)   &!       increment due to large-scale
                                    !       horizontal and vertical
  , w_inc(row_length,rows,nlevs)   &! Out
  , t_inc(row_length,rows,nlevs)   &! Out   Temp increment due to
                                    !       large-scale horizontal and
                                    !       vertical advection
                                    !       (K/s/day)
  , q_inc(row_length,rows,nwet)     !       Specific humidity increment
                                    !       due to large-scale horizontal
                                    !       and vertical advection
                                    !       (kg/kg)/s/day

!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
  INTEGER ::                       &
    i, i1, j, j1, k, l              ! Loop counters

  REAL ::                          &
    cdd                            &! Randomly sampled variables for
  , ct                             &!  temp. and dew pt. depression
  , cvn                            &! Randomly sampled variables for
  , cw                             &!  horizontal and vertical velocity
  , d0                             &! Randomly sampled dew pt. depression for
                                    ! 1st level
  , dewpt(row_length,rows,nwet,2)  &! Dew-point
  , dr                             &! Randomly sampled dew pt. depression (K)
  , f1                             &! Used in calc of advection term
                                    ! in equation 2.12
  , n1, n2                         &! Constants for interpolation
  , qk(row_length,rows,nwet)       &! Factor to prevent Q becoming negative
  , qrint(row_length,rows,nwet)    &! Specific humidity (interpolated values)
                                    ! (kg/kg)
  , rdaysteps, rstepcount          &! Real values of timesteps in day and
                                    ! timestep counter
  , t0                             &! Randomly sampled temp. for 1st level (K)
  , trint(row_length,rows,nlevs)   &! Temp. (interpolated values) (K)
  , vnrint(row_length,rows,nlevs)  &! Horizontal velocities (linearly
  , vprint(row_length,rows,nlevs)  &! interpolated values) (m/s)
  , wrint(row_length,rows,ntrop)   &! Vertical velocity (linearly
                                    ! interpolated values) (mb/s)
  , g05dde                          ! Function that samples randomly from a
                                    ! Gaussian distribution with a given mean
                                    ! and SD



  IF (stepcount == 1) THEN

    IF (daycount == 1) THEN
      i1 = 1
      j1 = 2
    ELSE
      i1 = 2
      j1 = 2

! Save state of Random Number Generator for continuation STATS run done from
! tape to allow for the first day of a STATS run, when G05DDE is used twice
! as many times (to set up 2 profiles) and so the variables after forcing on
! a continuation run would be different from an unbroken run.
! - daycount ne 1.

! DEPENDS ON: g05cfe
      CALL g05cfe(idum,iv,iy)
      DO j=1, row_length
        DO l=1, rows

          DO k=1, nlevs
            tr(j,l,k,j1-1) = tr(j,l,k,j1)
            vnr(j,l,k,j1-1) = vnr(j,l,k,j1)
            vpr(j,l,k,j1-1) = vpr(j,l,k,j1)
          END DO

          DO k=1, nwet
            qr(j,l,k,j1-1) = qr(j,l,k,j1)
          END DO

          DO k=1, ntrop
            wr(j,l,k,j1-1) = wr(j,l,k,j1)
          END DO

        END DO
      END DO
    END IF

!---------------------------------------------------------------------
!       Create new profiles (2 for 1st day and 1 from then on)
!       G05DDE(A,B) returns a pseudo-random real number taken from
!       a normal (Gaussian) distribution with mean A and SD B.
!---------------------------------------------------------------------

    DO j=1, row_length
      DO l=1, rows
        DO k=i1, j1

! DEPENDS ON: g05dde
          t0 = g05dde(tbar(j,l,1),tsd(j,l,1))
          tr(j,l,1,k) = t0+tdash(j,l,1)
! DEPENDS ON: g05dde
          d0 = g05dde(dbar(j,l,1),dsd(j,l,1))
          dr = d0+ddash(j,l,1)
          dewpt(j,l,1,k) = tr(j,l,1,k)-dr
! DEPENDS ON: g05dde
          vnr(j,l,1,k) = g05dde(vnbar(j,l,1),vnsd(j,l,1))
! DEPENDS ON: g05dde
          wr(j,l,1,k) = g05dde(wbar(j,l,1),wsd(j,l,1))

          DO i=1, nlevs-1
! DEPENDS ON: g05dde
            ct = g05dde(ctbar(j,l,i),ctsd(j,l,i))
            t0 = at(j,l,i)*t0+ct
            tr(j,l,i+1,k) = t0+tdash(j,l,i+1)
! DEPENDS ON: g05dde
            cvn = g05dde(cvnbar(j,l,i),cvnsd(j,l,i))
            vnr(j,l,i+1,k) = avn(j,l,i)*vnr(j,l,i,k)+cvn
          END DO

          DO i=1, nwet-1
! DEPENDS ON: g05dde
            cdd = g05dde(cdbar(j,l,i),cdsd(j,l,i))
            d0 = ad(j,l,i)*d0+cdd
            dr = d0+ddash(j,l,i+1)
            dewpt(j,l,i+1,k) = tr(j,l,i+1,k)-dr
          END DO

          DO i=1, ntrop-1
! DEPENDS ON: g05dde
            cw = g05dde(cwbar(j,l,i),cwsd(j,l,i))
            wr(j,l,i+1,k) = aw(j,l,i)*wr(j,l,i,k)+cw
          END DO

          DO i=1, nlevs
! DEPENDS ON: g05dde
            vpr(j,l,i,k) = g05dde(vpbar(j,l,i),vnsd(j,l,i))
          END DO

! After the first profile is set up on the first day, save state of
! Random Number Generator for continuation STATS run done from tape to allow
! for the first day of a STATS run, when G05DDF is used twice as many times
! (to set up 2 profiles) and so the variables after forcing on a continuation
! run would be different from an unbroken run.
! - daycount EQ 1.

          IF (k == 1 .AND. l == 1) THEN
! DEPENDS ON: g05cfe
            CALL g05cfe(idum,iv,iy)
          END IF
        END DO
      END DO
    END DO

    DO k=i1, j1
! DEPENDS ON: qsat
      CALL qsat (qr(1,1,1,k), dewpt(1,1,1,k), press, (row_length*rows*nwet))
    END DO
  END IF                     !  stepcount=1

! Interpolate between 2 values

  rdaysteps = REAL(daysteps)
  rstepcount = REAL(stepcount)
  n1 = (rdaysteps-rstepcount+1.) / rdaysteps
  n2 = (rstepcount-1.0) / rdaysteps

  DO k=1, nlevs
    DO j=1, rows
      DO l=1, row_length
        trint(l,j,k)  = n1 * tr(l,j,k,1)  + n2 * tr(l,j,k,2)
        vnrint(l,j,k) = n1 * vnr(l,j,k,1) + n2 * vnr(l,j,k,2)
        vprint(l,j,k) = n1 * vpr(l,j,k,1) + n2 * vpr(l,j,k,2)
      END DO
    END DO
  END DO

  DO k=1, nwet
    DO j=1, rows
      DO l=1, row_length
        qrint(l,j,k) = n1 * qr(l,j,k,1) + n2 * qr(l,j,k,2)
      END DO
    END DO
  END DO

  DO k=1, ntrop
    DO j=1, rows
      DO l=1, row_length
        wrint(l,j,k) = n1 * wr(l,j,k,1) + n2 * wr(l,j,k,2)
      END DO
    END DO
  END DO

! Set U and V increments
  DO k=1, nlevs
    DO j=1, rows
      DO l=1, row_length
        u_inc(l,j,k) = - u(l,j,k) + vnrint(l,j,k)
        v_inc(l,j,k) = - v(l,j,k) + vprint(l,j,k)
      END DO
    END DO
  END DO

  DO k=1, ntrop
    DO j=1, rows
      DO l=1, row_length
        w_inc(l,j,k) = - w(l,j,k) + wrint(l,j,k)
      END DO
    END DO
  END DO

!---------------------------------------------------------------------
!     Add vertical advection increments to T and Q (eqns. 2 and 3)
!---------------------------------------------------------------------

  DO k=1, nlevs
    DO j=1, rows
      DO l=1, row_length
        t_inc(l,j,k) = 0.0
      END DO
    END DO
  END DO

  DO k=1, nwet
    DO j=1, rows
      DO l=1, row_length
        q_inc(l,j,k) = 0.0
      END DO
    END DO
  END DO

  DO j=1, rows
    DO l=1, row_length
      dt1 = t(l,j,2) - t(l,j,1)
      dq1 = q(l,j,2) - q(l,j,1)
      DO i=2, ntrop
        dt2 = t(l,j,i+1) - t(l,j,i)
        dq2 = q(l,j,i+1) - q(l,j,i)
        f1 = -wrint(l,j,i) * timestep * rpress(l,j,i)
        t_inc(l,j,i) = t_inc(l,j,i) + f1 * (dt2 * px(l,j,i) +       &
            dt1 * px(l,j,i-1)                                       &
          - (dt1 + dt2) * py(l,j,i-1) - .2856 * t(l,j,i))
        q_inc(l,j,i) = q_inc(l,j,i) + f1 * (dq2 * px(l,j,i) +       &
        dq1 * px(l,j,i-1)                                           &
          - (dq1+dq2) * py(l,j,i-1))
        dt1 = dt2
        dq1 = dq2
      END DO
    END DO
  END DO                     ! l

! Printout increments due to vert. advection if required
!
! These are commented out for now so that the call can be added in at
! a later date if so required.
!      If (prinstat)
!        Call PRINTSUB(
!     ! In
!        row_length, rows, nlevs, nwet,
!     !
!        ' Temperature/moisture profiles + incs due to vert advection '
!        ,stepcount , dayno, t, q, t_inc, q_inc)

!---------------------------------------------------------------------
!     Add horizontal increments to T and Q (eqns. 2 and 3)
!---------------------------------------------------------------------

  DO k=1, nlevs
    DO j=1, rows
      DO l=1, row_length
        t_inc(l,j,k) = t_inc(l,j,k) + timestep * ABS(vnrint(l,j,k)) &
          *    (trint(l,j,k)-t(l,j,k)) / deltan(l,j)
      END DO
    END DO
  END DO

! This section prevents the increment from allowing
! Q to become negative (lowest value Q can take is 1.E-6)

  DO k=1, nwet
    DO j=1, rows
      DO l=1, row_length
        q_inc(l,j,k) = q_inc(l,j,k) + timestep * ABS(vnrint(l,j,k)) &
          *          (qrint(l,j,k)-q(l,j,k))                        &
          /          deltan(l,j)
        qk(l,j,k) = -1 * q(l,j,k) + 1.e-6
        IF (q_inc(l,j,k)  <   qk(l,j,k)) THEN
          q_inc(l,j,k) = qk(l,j,k)
        END IF
      END DO
    END DO
  END DO

! Printout increments due to horiz. advection if required

  IF (prinstat) THEN
! DEPENDS ON: printsub
    CALL printsub                                                             &
      ( row_length*rows, nlevs, nwet                                          &
      , ' Temperature/moisture profiles + incs due to advection', stepcount   &
      , dayno, t(1,1:1,1:nlevs), q(1,1:1,1:nlevs), t_inc(1,1:1,1:nlevs)       &
      , q_inc(1,1:1,1:nlevs) )
  END IF

  RETURN
END SUBROUTINE statstep

