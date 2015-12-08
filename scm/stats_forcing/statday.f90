! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE StatDay
! Purpose:-           To calculate statistical forcing required each
!                     CLIM_CHANGE days
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!=====================================================================

SUBROUTINE statday                                                            &
  ( row_length, rows, nlevs, nwet, ntrop, atime, btime, dayno, deltan         &
  , daycount, tbara, tbarb, tsda, tsdb, dbara, dbarb, vnbara, vnbarb, vnsda   &
  , vnsdb, vpbara, vpbarb, wbara, wbarb, wsda, wsdb, alfada, alfadb, pa, pb   &
  , p, tgrada, tgradb, dgrada, dgradb, cort, cord, corvn, corw, tdash, ddash  &
  , ctbar, ctsd, at, cdbar, cdsd, ad, cvnbar, cvnsd, avn, cwbar, cwsd, aw     &
  , tbar, tsd, dbar, dsd, vnbar, vnsd, vpbar, wbar, wsd )

  IMPLICIT NONE

!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------

  INTEGER ::    &
    row_length  &! In dimensions of arrays
  , rows        &
  , nlevs       &! In Number of levels of the model.
  , nwet        &! In Number of humidity model levels
  , ntrop       &! In Max number of levels in the troposphere
  , dayno       &! In Dayno. relative to winter solstice
  , daycount     ! In Daynumber (ie 1 is 1st january)

  REAL ::                         &
    alfada(row_length,rows)       &! In  Amplitude and mean of seasonal
  , alfadb(row_length,rows)       &!     variation of tuning factor
  , at(row_length,rows,nlevs-1)   &! Out Variable a in eqn 2.22
  , ad(row_length,rows,nwet-1)    &!     used to calculate mean of
                                   !     random variable for
                                   !     temp. and dew point depression
  , atime, btime                  &! In  Constants for calculating
                                   !     annual cycle in eqn. 2.33
  , avn(row_length,rows,nlevs-1)  &! Out Variable a in eqn 2.22
  , aw(row_length,rows,ntrop-1)   &!     used to calculate mean of random
                                   !     variable for horiz. and vert
                                   !     velocity
  , cord(row_length,rows)         &! In  vertical correlation
                                   !     coeff. for dew pt. depress(0.9)
  , cort(row_length,rows)         &! In  vertical correlation
                                   !     coeff. for temp. (0.9)
  , corvn(row_length,rows)        &! In  vertical correlation
                                   !     coeff. for velocity VN (0.5)
  , corw(row_length,rows)         &! In  vertical correlation
                                   !     coeff. for velocity W (0.5)
  , cdbar(row_length,rows,nwet)   &! Out Mean and SD of random variable
  , cdsd(row_length,rows,nwet)    &!     for dew pt. depression (eqns
                                   !     2.22 and 2.23)
  , ctbar(row_length,rows,nlevs)  &! Out Mean and SD of random
  , ctsd(row_length,rows,nlevs)   &!     variable for temp. (eqns 2.22
                                   !     and 2.23)
  , cvnbar(row_length,rows,nlevs) &! Out Mean and SD of random variable
  , cvnsd(row_length,rows,nlevs)  &!     for velocity VN (eqns 2.22
                                   !     and 2.23)
  , cwbar(row_length,rows,ntrop)  &! Out Mean and SD of random variable
  , cwsd(row_length,rows,ntrop)   &!     for vertical velocity (eqns 2.22
                                   !     and 2.23)
  , dbar(row_length,rows,nwet)    &! Out Mean and SD dewpoint
  , dsd(row_length,rows,nwet)     &!     depression at daynumber relative
                                   !     to winter solstice (K)
  , dbara(row_length,rows,nwet)   &! In  Amplitude and mean of seasonal
  , dbarb(row_length,rows,nwet)   &!     variation of mean dew pt.
                                   !     depression (K)
  , ddash(row_length,rows,nwet)   &! Out Dew pt. corrections (K)
  , deltan(row_length,rows)       &! In  Radius of area (m)
  , dgrada(row_length,rows,nwet)  &! In  Amplitude and mean of seasonal
  , dgradb(row_length,rows,nwet)  &!     variation of dew pt. depression
                                   !     gradient (K/km)
  , p(row_length,rows,nlevs+1)    &! Out Pressure on rho levels (Pa)
  , pa(row_length,rows,nlevs+1)   &! Out Amplitude and mean of seasonal
  , pb(row_length,rows,nlevs+1)   &!     variation of pressure
  , tdash(row_length,rows,nlevs)  &!     Temp. correction (K)
  , tbar(row_length,rows,nlevs)   &! Out Mean and SD temperature at
  , tsd(row_length,rows,nlevs)    &!     daycount days from winter
                                   !     solstice (K)
  , tbara(row_length,rows,nlevs)  &! In  Amplitude and mean of
  , tbarb(row_length,rows,nlevs)  &!     seasonal variation of temp. (K)
  , tgrada(row_length,rows,nlevs) &! In  Amplitude and mean of seasonal
  , tgradb(row_length,rows,nlevs) &!     variation of temp. gradient
                                   !     (K/km)
  , tsda(row_length,rows,nlevs)   &! In  Amplitude and mean of seasonal
  , tsdb(row_length,rows,nlevs)   &!     variation of SD of temp. (K)
  , vnbar(row_length,rows,nlevs)  &! Out Mean and SD velocity VN at
  , vnsd(row_length,rows,nlevs)   &!     daycount days from
                                   !     winter solstice (m/s)
  , vpbar(row_length,rows,nlevs)  &! Out Mean  velocity VP at
                                   !     daycount days from
                                   !     winter solstice (m/s)
  , vnbara(row_length,rows,nlevs) &! In  Amplitude and mean of seasonal
  , vnbarb(row_length,rows,nlevs) &!     variation of velocity VN
                                   !     (m/s)
  , vnsda(row_length,rows,nlevs)  &! In  Amplitude and mean of seasonal
  , vnsdb(row_length,rows,nlevs)  &!     variation of SD of velocity VN
                                   !     (m/s)
  , vpbara(row_length,rows,nlevs) &! In  Amplitude and mean of seasonal
  , vpbarb(row_length,rows,nlevs) &!     variation of velocity VP
                                   !     (m/s)
  , wbar(row_length,rows,ntrop)   &! Out Mean and SD vertical
  , wsd(row_length,rows,ntrop)    &!     velocity at daycount days
                                   !     from winter solstice (mb/s)
  , wbara(row_length,rows,ntrop)  &! In  Amplitude and mean of seasonal
  , wbarb(row_length,rows,ntrop)  &!     variation of SD of vert. vel.
                                   !     (mb/s)
  , wsda(row_length,rows,ntrop)   &! In  Amplitude and mean of seasonal
  , wsdb(row_length,rows,ntrop)    !     variation of SD of vert. vel.
                                   !     (mb/s)

!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------

  REAL ::                         &
    alfad(row_length,rows)        &! Tuning factor at daycount days
                                   ! from winter solstice
  , daynew                        &! function to calculate SIN
                                   ! of argument (eqn 2.33)
  , dgrad                         &! dew pt. depression gradient
  , tgrad                         &! Temp. gradient
  , xt                             ! Argument of SIN distribution (eqn. 2.33)

  INTEGER ::                      &
    i, j, k                        ! Loop counter

! Calculate argument of SIN (in eqn. 12)
  IF (daycount == 1) THEN

! DEPENDS ON: daynew
    xt = daynew(atime, btime, dayno)

  ELSE

! DEPENDS ON: daynew
    xt = daynew(atime, btime, dayno+1)

  END IF

! Calculate sinusoidal distribution (eqn. 12)

! DEPENDS ON: xnew
  CALL xnew (tbar, tbara, tbarb, row_length, rows, nlevs, xt)

! DEPENDS ON: xnew
  CALL xnew (tsd, tsda, tsdb, row_length, rows, nlevs, xt)

! DEPENDS ON: xnew
  CALL xnew (dbar, dbara, dbarb, row_length, rows, nwet, xt)

! DEPENDS ON: xnew
  CALL xnew (vnbar, vnbara, vnbarb, row_length, rows, nlevs, xt)

! DEPENDS ON: xnew
  CALL xnew (vnsd, vnsda, vnsdb, row_length, rows, nlevs, xt)

! DEPENDS ON: xnew
  CALL xnew (vpbar, vpbara, vpbarb, row_length, rows, nlevs, xt)

! DEPENDS ON: xnew
  CALL xnew (wbar, wbara, wbarb, row_length, rows, ntrop, xt)

! DEPENDS ON: xnew
  CALL xnew (wsd, wsda, wsdb, row_length, rows, ntrop, xt)

! DEPENDS ON: xnew
  CALL xnew (alfad, alfada, alfadb, row_length, rows, 1, xt)

! DEPENDS ON: xnew
  CALL xnew (p, pa, pb, row_length, rows, 1, xt)

  DO i=1, row_length
    DO j=1, rows

      DO k=1, nlevs
        tgrad = tgrada(i,j,k)*xt + tgradb(i,j,k)
        IF (vnbar(i,j,k)  >   0.0) THEN
          tgrad = -tgrad
        END IF

! Calculate corrections for temp. and dew pt. depression

        tdash(i,j,k)=deltan(i,j)*tgrad*0.001
      END DO

      DO k=1, nwet
        dgrad = dgrada(i,j,k) * xt + dgradb(i,j,k)
        dsd(i,j,k) = alfad(i,j) * tsd(i,j,k)
        IF (vnbar(i,j,k)  >   0.0) THEN
          dgrad = -dgrad
        END IF

! Calculate corrections for temp. and dew pt. depression

        ddash(i,j,k) = deltan(i,j) * dgrad*0.001
      END DO

    END DO
  END DO


! Calculate mean and SD of random variable (eqns. 6 and 7)

! DEPENDS ON: acinit
  CALL acinit (tbar, tsd, at, ctbar, ctsd, cort, nlevs, row_length, rows)

! DEPENDS ON: acinit
  CALL acinit (dbar, dsd, ad, cdbar, cdsd, cord, nwet, row_length, rows)

! DEPENDS ON: acinit
  CALL acinit (vnbar, vnsd, avn, cvnbar, cvnsd, corvn, nlevs, row_length, rows)

! DEPENDS ON: acinit
  CALL acinit (wbar, wsd, aw, cwbar, cwsd, corw, ntrop, row_length, rows)

  RETURN
END SUBROUTINE statday

