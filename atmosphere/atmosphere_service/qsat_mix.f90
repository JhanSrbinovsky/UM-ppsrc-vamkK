! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Saturation Specific Humidity/mixing Scheme(Qsat):Vapour to Liquid/Ice
!
! Subroutine Interface:
      SUBROUTINE QSAT_mix (                                             &
!      Output field
     &  QmixS                                                           &
!      Input fields
     &, T, P                                                            &
!      Array dimensions
     &, NPNTS                                                           &
!      logical control
     &, lq_mix                                                          &
     &  )

      USE qsat_data, ONLY:                                              &
         repsilon,                                                       &
         one_minus_epsilon,                                             &
         T_low,                                                         &
         T_high,                                                        &
         delta_T,                                                       &
         zerodegc,                                                      &
         ES

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Purpose:
!   Returns a saturation specific humidity or mixing ratio given a
!   temperature and pressure using the saturation vapour pressure
!   calculated using the Goff-Gratch formulae, adopted by the WMO as
!   taken from Landolt-Bornstein, 1987 Numerical Data and Functional
!   Relationships in Science and Technolgy. Group V/vol 4B meteorology.
!   Phyiscal and Chemical properties or air, P35.
!
!   Values in the lookup table are over water above 0 deg C and over ice
!   below this temperature.
!
! Method:
!   Uses lookup tables to find eSAT, calculates qSAT directly from that.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Service
!
! Code Description:
!   Language: FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No.29
!
! Declarations:
!
!  Global Variables:----------------------------------------------------

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! arguments with intent in. ie: input variables.
!
      Integer, intent(in) ::                                            &
     &  npnts    ! Points (=horizontal dimensions) being processed by
                 !  qSAT scheme.
!
      Real, intent(in)  ::                                              &
     &  T(npnts)                                                        &
                      !  Temperature (K).
     &, P(npnts)      !  Pressure (Pa).

      logical, intent(in)  ::                                           &
     &  lq_mix      !  .true. return qsat as a mixing ratio
                    !  .false. return qsat as a specific humidity

!
! arguments with intent out
!
      Real, intent(out)   ::                                            &
     &  QmixS(npnts)  ! Output Saturation mixing ratio or saturation
                      ! specific humidity at temperature T and pressure
                      ! P (kg/kg).

!-----------------------------------------------------------------------
!  Local scalars
!-----------------------------------------------------------------------
      Integer :: itable, itable_p1    ! Work variables

      Real :: atable, atable_p1       ! Work variables

      Real :: fsubw, fsubw_p1
            ! FACTOR THAT CONVERTS FROM SAT VAPOUR PRESSURE IN A PURE
            ! WATER SYSTEM TO SAT VAPOUR PRESSURE IN AIR.

      Real :: TT, TT_P1

      Real, parameter :: R_delta_T = 1./delta_T

      Integer :: I, II

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('QSAT_MIX',zhook_in,zhook_handle)

!
! loop over points
!
      Do i = 1, npnts-1,2

!      Compute the factor that converts from sat vapour pressure in a
!      pure water system to sat vapour pressure in air, fsubw.
!      This formula is taken from equation A4.7 of Adrian Gill's book:
!      Atmosphere-Ocean Dynamics. Note that his formula works in terms
!      of pressure in mb and temperature in Celsius, so conversion of
!      units leads to the slightly different equation used here.
 
        fsubw    = 1.0 + 1.0E-8 * P(I)   * ( 4.5 +                      &
          6.0E-4*( T(I) - zerodegC )   * ( T(I) - zerodegC ) )
        fsubw_p1 = 1.0 + 1.0E-8 * P(I+1) * ( 4.5 +                      &
          6.0E-4*( T(I+1) - zerodegC ) * ( T(I+1) - zerodegC ) )

!      Use the lookup table to find saturated vapour pressure, and store
!      it in qmixs.
!
        TT = MAX(T_low,T(I))
        TT = MIN(T_high,TT)
        atable = (TT - T_low + delta_T) * R_delta_T
        itable = atable
        atable = atable - itable

        TT_p1 = MAX(T_LOW,T(I+1))
        TT_p1 = MIN(T_HIGH,TT_p1)
        ATABLE_p1 = (TT_p1 - T_LOW + DELTA_T) * R_DELTA_T
        ITABLE_p1 = ATABLE_p1
        ATABLE_p1 = ATABLE_p1 - ITABLE_p1

        QmixS(I)   = (1.0 - atable)   *ES(itable) + atable*ES(itable+1)
        QmixS(I+1) = (1.0 - atable_p1)*ES(itable_p1) +                 &
                                              atable_p1*ES(itable_p1+1)
 
!      Multiply by fsubw to convert to saturated vapour pressure in air
!      (equation A4.6 of Adrian Gill's book).
!
        QmixS(I)   = QmixS(I)   * fsubw
        QmixS(I+1) = QmixS(I+1) * fsubw_p1
!
!      Now form the accurate expression for qmixs, which is a rearranged
!      version of equation A4.3 of Gill's book.
!
!-----------------------------------------------------------------------
! For mixing ratio,  rsat = epsilon *e/(p-e)
! e - saturation vapour pressure
! Note applying the fix to qsat for specific humidity at low pressures
! is not possible, this implies mixing ratio qsat tends to infinity.
! If the pressure is very low then the mixing ratio value may become
! very large.
!-----------------------------------------------------------------------

       if (lq_mix) then

         QmixS(I)   = ( repsilon*QmixS(I) ) /                           &
                    ( MAX(P(I),  1.1*QmixS(I))   - QmixS(I) )
         QmixS(I+1) = ( repsilon*QmixS(I+1) ) /                         &
                    ( MAX(P(I+1),1.1*QmixS(I+1)) - QmixS(I+1) )

!-----------------------------------------------------------------------
! For specific humidity,   qsat = epsilon*e/(p-(1-epsilon)e)
!
! Note that at very low pressure we apply a fix, to prevent a
! singularity (qsat tends to 1. kg/kg).
!-----------------------------------------------------------------------
        else

          QmixS(I)   = ( repsilon*QmixS(I) ) /                          &
                ( MAX(P(I),  QmixS(I))   - one_minus_epsilon*QmixS(I) )
          QmixS(I+1) = ( repsilon*QmixS(I+1) ) /                        &
                ( MAX(P(I+1),QmixS(I+1)) - one_minus_epsilon*QmixS(I+1))

        endif  ! test on lq_mix
!
      End Do  ! Npnts_do_1

      II = I
      Do I = II, NPNTS
        fsubw = 1.0 + 1.0E-8*P(I)*( 4.5 +                               &
     &    6.0E-4*( T(I) - zerodegC )*( T(I) - zerodegC ) )
!
        TT = MAX(T_low,T(I))
        TT = MIN(T_high,TT)
        atable = (TT - T_low + delta_T) * R_delta_T
        itable = atable
        atable = atable - itable

        QmixS(I) = (1.0 - atable)*ES(itable) + atable*ES(itable+1)
        QmixS(I) = QmixS(I) * fsubw

       if (lq_mix) then
         QmixS(I) = ( repsilon*QmixS(I) ) /                             &
     &              ( MAX(P(I),1.1*QmixS(i)) - QmixS(I) )
       else
         QmixS(I) = ( repsilon*QmixS(I) ) /                             &
     &          ( MAX(P(I),QmixS(I)) - one_minus_epsilon*QmixS(I) )
       endif  ! test on lq_mix
      End Do  ! Npnts_do_1

      IF (lhook) CALL dr_hook('QSAT_MIX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE QSAT_mix
! ======================================================================
