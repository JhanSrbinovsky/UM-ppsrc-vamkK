! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Saturation Specific Humidity Scheme (Qsat_Wat): Vapour to Liquid.
! Subroutine Interface:
      SUBROUTINE qsat_wat_mix (                                         &
!      Output field
     &  QmixS                                                           &
!      Input fields
     &, T, P                                                            &
!      Array dimensions
     &, NPNTS                                                           &
!      logical control
     &, lq_mix                                                          &
     &  )
 
      USE qsat_wat_data, ONLY:                                          &
         repsilon,                                                      &
         one_minus_epsilon,                                             &
         T_low,                                                         &
         T_high,                                                        &
         delta_T,                                                       &
         zerodegc,                                                      &
         ES
      USE vectlib_mod, ONLY: oneover_v
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
 
! Purpose:
!   Returns a saturation specific humidity or mixing ratio given a
!   temperature and pressure using the saturation vapour pressure
!   calculated using the Goff-Gratch formulae, adopted by the WMO as
!   taken from Landolt-Bornstein, 1987 Numerical Data and Functional
!   Relationships in Science and Technolgy. Group V/vol 4B meteorology.
!   Phyiscal and Chemical properties or air, P35.
!
!   Values in the lookup table are over water above and below 0 deg C.
!
!   Note : For vapour pressure over water this formula is valid for
!   temperatures between 373K and 223K. The values for saturated vapour
!   over water in the lookup table below are out of the lower end of
!   this range. However it is standard WMO practice to use the formula
!   below its accepted range for use with the calculation of dew points
!   in the upper atmosphere
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
!
!  Subroutine Arguments:------------------------------------------------
!
! arguments with intent in. ie: input variables.
 
      Integer, intent(in) :: npnts
      ! Points (=horizontal dimensions) being processed by qSAT scheme.
 
      Real, intent(in)  :: T(npnts) !  Temperature (K).
      Real, intent(in)  :: P(npnts) !  Pressure (Pa).  

      logical, intent(in)  :: lq_mix
                    !  .true. return qsat as a mixing ratio
                    !  .false. return qsat as a specific humidity
 
! arguments with intent out
 
      Real, intent(out)   ::  QmixS(npnts)
             ! Output Saturation mixing ratio or saturation specific
             ! humidity at temperature T and pressure P (kg/kg).

!-----------------------------------------------------------------------
!  Local scalars
!-----------------------------------------------------------------------
 
      Integer :: itable    ! Work variables
 
      Real :: atable       ! Work variables
 
      Real :: fsubw, fsubw_p1    
            ! FACTOR THAT CONVERTS FROM SAT VAPOUR PRESSURE IN A PURE 
            ! WATER SYSTEM TO SAT VAPOUR PRESSURE IN AIR.

      Real :: TT
      
      Real :: vect_in(npnts), vect_out(npnts)
            ! temp array for calculating reciprocals
      Integer, parameter :: block_length = 400
      Integer :: end_pt 

      Real, parameter :: R_delta_T = 1./delta_T
 
      Integer :: I, II, J

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

 
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('QSAT_WAT_MIX',zhook_in,zhook_handle)
      
      
! loop over points
      Do j=1,npnts, block_length 

        end_pt = MIN(j+block_length-1,npnts)

        Do i = j, end_pt
          
          TT = MAX(T_low,T(I))
          TT = MIN(T_high,TT)
          atable = (TT - T_low + delta_T) * R_delta_T
          itable = atable
          atable = atable - itable
          
!      Use the lookup table to find saturated vapour pressure, and store
!      it in qmixs.
          
          QmixS(I)   = (1.0 - atable)    * ES(itable)    +              &
               atable*ES(itable+1)
        End Do

        If (lq_mix) then   

          Do i = j, end_pt
!      Compute the factor that converts from sat vapour pressure in a
!      pure water system to sat vapour pressure in air, fsubw.
!      This formula is taken from equation A4.7 of Adrian Gill's book:
!      Atmosphere-Ocean Dynamics. Note that his formula works in terms
!      of pressure in mb and temperature in Celsius, so conversion of
!      units leads to the slightly different equation used here.
 
            fsubw    = 1.0 + 1.0E-8*P(I)   * ( 4.5 +                    &
                 6.0E-4*( T(I)   - zerodegC ) * ( T(I) - zerodegC ) )
                
!      Multiply by fsubw to convert to saturated vapour pressure in air
!      (equation A4.6 of Adrian Gill's book).
 
            QmixS(I)   = QmixS(I)   * fsubw
          
!      Now form the accurate expression for qmixs, which is a rearranged
!      version of equation A4.3 of Gill's book.
 
!-----------------------------------------------------------------------
! For mixing ratio,  rsat = epsilon *e/(p-e)
! e - saturation vapour pressure
! Note applying the fix to qsat for specific humidity at low pressures
! is not possible, this implies mixing ratio qsat tends to infinity.
! If the pressure is very low then the mixing ratio value may become
! very large.
!-----------------------------------------------------------------------

            vect_in(I) = ( MAX(P(I),  1.1*QmixS(I))   - QmixS(I) )
          
          End Do

        Else
          
          Do i = j, end_pt
!      Compute the factor that converts from sat vapour pressure in a
!      pure water system to sat vapour pressure in air, fsubw.
!      This formula is taken from equation A4.7 of Adrian Gill's book:
!      Atmosphere-Ocean Dynamics. Note that his formula works in terms
!      of pressure in mb and temperature in Celsius, so conversion of
!      units leads to the slightly different equation used here.
            
            fsubw    = 1.0 + 1.0E-8*P(I)   * ( 4.5 +                    &
                 6.0E-4*( T(I)   - zerodegC ) * ( T(I) - zerodegC ) )
        
!      Multiply by fsubw to convert to saturated vapour pressure in air
!      (equation A4.6 of Adrian Gill's book).
 
            QmixS(I)   = QmixS(I)   * fsubw
        
!      Now form the accurate expression for qmixs, which is a rearranged
!      version of equation A4.3 of Gill's book.

!-----------------------------------------------------------------------
! For specific humidity,   qsat = epsilon*e/(p-(1-epsilon)e)
!
! Note that at very low pressure we apply a fix, to prevent a
! singularity (qsat tends to 1. kg/kg).
!-----------------------------------------------------------------------
            vect_in(I) = (MAX(P(I), QmixS(I))-one_minus_epsilon*QmixS(I))

          End Do    
        End If
      End Do
       
      CALL oneover_v(npnts, vect_in, vect_out)
      
      Do i=1, npnts
        QMixS(I) = repsilon*QmixS(I) * vect_out(I)
      End Do
      
      IF (lhook) CALL dr_hook('QSAT_WAT_MIX',zhook_out,zhook_handle)
      RETURN
    END SUBROUTINE qsat_wat_mix
! ======================================================================
