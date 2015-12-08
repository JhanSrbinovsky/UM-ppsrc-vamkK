! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description: 
!   Calculate solar zenith angle for photolysis applications.
!   Adapted from UM SOLANG routine, but cos(SZA) is not 
!   assumed to be 0 at night.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE UKCA_SOLANG (sindec, t, dt, eqt, sinlat, longit,      &
                              k, cosz)        

      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE                                                    
    

      INTEGER, INTENT(IN) ::  k         ! Number of points

      REAL, INTENT(IN)    :: sindec     ! Sin(solar declination)
      REAL, INTENT(IN)    :: t          ! Start time (GMT)       
      REAL, INTENT(IN)    :: dt         ! timestep
      REAL, INTENT(IN)    :: eqt        ! Eqn of time
      
      REAL, DIMENSION(K), INTENT(IN)    :: sinlat  ! sin(latitude)
      REAL, DIMENSION(K), INTENT(IN)    :: longit  ! longitude

      REAL, DIMENSION(K), INTENT(OUT)   :: cosz    ! Mean Cos(sza)                                            

!     Local variables

      INTEGER          :: j         ! Loop counter over points            

      REAL, PARAMETER  :: s2r = pi/43200.0  ! sec-to-rads converter

      REAL             :: sinsin    ! Products of the sines and of the cosines    
      REAL             :: coscos    ! of solar declination and of latitude.       
      REAL             :: hat       ! Local hour angle at the start time.         
      REAL             :: difsin    ! A difference-of-sines intermediate value    
      REAL             :: trad      ! Start and length of timestep (T & DT)
      REAL             :: dtrad     ! converted to radians after midday GMT

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_SOLANG',zhook_in,zhook_handle)
      trad  =  t * s2r - pi                                                   
      dtrad = dt * s2r                                                      

      DO j = 1, k                                    
       sinsin  = sindec * sinlat(j)                                          
       coscos  = SQRT( (1.-sindec**2) * (1.-sinlat(j)**2) )                  
       hat     = longit(j) + trad + eqt                              
       difsin  = SIN(hat + dtrad) - SIN(hat)                 
       cosz(j) = difsin*coscos/dtrad + sinsin          
      END DO   
      
      IF (lhook) CALL dr_hook('UKCA_SOLANG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_SOLANG
