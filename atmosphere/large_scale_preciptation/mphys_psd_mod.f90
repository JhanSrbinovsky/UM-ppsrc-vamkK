! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE mphys_psd_mod

! Description:
! Holds particle-size distribution constants required by the large-scale
! precipitation scheme
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

  IMPLICIT NONE


!---------------------------------------------------------------------
! CALCULATION OF FALL SPEEDS
!---------------------------------------------------------------------

      ! Do we wish to calculate the ice fall velocities?
      ! TRUE if calculate speeds, FALSE if specify speeds

      LOGICAL, PARAMETER :: l_calcfall = .TRUE.

!---------------------------------------------------------------------
! RAIN PARAMETERS
!---------------------------------------------------------------------


      ! Drop size distribution for rain: N(D) =  N0 D^m exp(-lambda D)
      ! where N0 = X1R lambda^X2R  and m = X4R

!     REAL, PARAMETER :: x1r is set in the UMUI
!     REAL, PARAMETER :: x2r is set in the UMUI
!     REAL, PARAMETER :: x4r is set in the UMUI

      ! Fall speed diameter relationship for rain: vt(D) = CR D^DR

     REAL, PARAMETER :: cr = 386.8    
     REAL, PARAMETER :: dr = 0.67     

      ! Abel & Shipway (2007) fall speed diameter relationships for rain:
      !     vt(D) = C1R D^D1R exp(-h1r D) + C2R D^D2R exp(-h2r D)

     REAL, PARAMETER :: c1r = 4845.1
     REAL, PARAMETER :: d1r = 1.0
     REAL, PARAMETER :: h1r = 195.0
     REAL, PARAMETER :: c2r = -446.009
     REAL, PARAMETER :: d2r = 0.782127
     REAL, PARAMETER :: h2r = 4085.35

!---------------------------------------------------------------------
! ICE AGGREGATE PARAMETERS
!---------------------------------------------------------------------

      ! Particle size distribution for ice aggregates: 
      ! N(D) = N0 D^m exp(-lambda D)
      ! where N0 = X1I TCG lambda^X2I, 
      ! m = X4I and TCG = exp(- X3I T[deg C])

!     REAL, PARAMETER :: x1i is set in the UMUI

     REAL, PARAMETER :: x2i = 0.0
     REAL, PARAMETER :: x3i = 0.1222
     REAL, PARAMETER :: x4i = 0.0

      ! Mass diameter relationship for ice:  m(D) = AI D^BI
      ! These are set in the UMUI. 
      ! Recommended values for the generic particle size distribution are
      ! from Brown and Francis and are
      ! AI = AIC = 1.85E-2, BI = BIC = 1.9. If l_calcfall is changed from .true.
      ! then the generic psd values should be set below to ci = cic = 8.203
      ! and di = dic = 0.2888

      ! REAL, PARAMETER :: ai   )
      ! REAL, PARAMETER :: bi   ) Parameters are set in the UMUI
      ! REAL, PARAMETER :: aic  )
      ! REAL, PARAMETER :: bic  )

      ! The area diameter relationships are only used if
      ! L_CALCFALL = .TRUE.
      ! Area diameter relationship for ice:  Area(D) = RI D^SI

     REAL, PARAMETER :: ri = 0.131
     REAL, PARAMETER :: si = 1.88


      ! The Best/Reynolds relationships are only used if
      ! L_CALCFALL = .TRUE.
      ! Relationship between Best number and Reynolds number:
      ! Re(D)  = LSP_EI(C) Be^LSP_FI(C)
     
      ! REAL, PARAMETER :: lsp_ei  ) 
      ! REAL, PARAMETER :: lsp_fi  ) Set in the UMUI

      ! The fall speeds of ice particles are only used if
      ! L_CALCFALL = .FALSE.
      ! Fall speed diameter relationships for ice:
      ! vt(D) = CI D^DI

     REAL, PARAMETER :: ci0 = 14.3
     REAL, PARAMETER :: di0 = 0.416 

!---------------------------------------------------------------------
! ICE CRYSTAL PARAMETERS
!---------------------------------------------------------------------

      ! Particle size distribution for ice crystals: 

      ! N(D) = N0 D^m exp(-lambda D)
      ! where N0 = X1I TCG lambda^X2I, 
      ! m = X4I and TCG = exp(- X3I T[deg C])

!    REAL, PARAMETER :: x1ic is set in the UMUI

     REAL, PARAMETER :: x2ic = 0.0
  
     REAL, PARAMETER :: x3ic = 0.1222
     REAL, PARAMETER :: x4ic = 0.0

      ! The area diameter relationships are only used if
      ! L_CALCFALL = .TRUE.
      ! Area diameter relationship for ice:  Area(D) = RI D^SI


     REAL, PARAMETER :: ric = 0.131
     REAL, PARAMETER :: sic = 1.88

      ! The Best/Reynolds relationships are only used if
      ! L_CALCFALL = .TRUE.
      ! Relationship between Best number and Reynolds number:
      ! Re(D)  = LSP_EI(C) Be^LSP_FI(C)
     
      ! REAL, PARAMETER :: lsp_eic  ) 
      ! REAL, PARAMETER :: lsp_fic  ) Set in the UMUI
  
      ! The fall speeds of ice particles are only used if
      ! L_CALCFALL = .FALSE.
      ! Fall speed diameter relationships for ice:
      ! vt(D) = CI D^DI

     REAL, PARAMETER :: cic0 = 74.5
     REAL, PARAMETER :: dic0 = 0.640

!---------------------------------------------------------------------
! GRAUPEL PARAMETERS
!---------------------------------------------------------------------

      ! Drop size distribution for graupel: N(D) =  N0 D^m exp(-lambda D)
      ! where N0 = X1G lambda^X2G  and m = X4G

     REAL, PARAMETER :: x1g = 5.E25
     REAL, PARAMETER :: x2g = -4.0
     REAL, PARAMETER :: x4g = 2.5

      ! Mass diameter relationship for graupel:  m(D) = AG D^BG

     REAL, PARAMETER :: ag  = 261.8
     REAL, PARAMETER :: bg  = 3.0

      ! Fall speed diameter relationship for graupel: vt(D) = CG D^DG

     REAL, PARAMETER :: cg  = 253.0
     REAL, PARAMETER :: dg  = 0.734


END MODULE mphys_psd_mod
