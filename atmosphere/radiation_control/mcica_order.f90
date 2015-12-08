! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to calculate the first subcol for each k-term, and no of
!  subcols per band

! Method:
!       Trivial

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
SUBROUTINE mcica_order(ierr                                             &
!                       General spectral properties
, i_first_band, i_last_band                                             &
!                       Gaseous absorption
, i_band_esft, index_absorb                                             &
!                       Dimensions of arrays
, nd_esft_term, nd_band, nd_species                                     &
  )

  USE rad_pcf
  USE mcica_mod
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE um_types
  IMPLICIT NONE



!                       Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
    ierr
!          Error flag

!                       Range of spectral bands:
INTEGER, INTENT(IN) ::                                                  &
    i_first_band                                                        &
!         First band
  , i_last_band
!         Last band

!                       Gaseous absorption
INTEGER, INTENT(IN) ::                                                  &
   i_band_esft(nd_band, nd_species)                                     &
!         Number of terms in band
  , index_absorb(nd_species, nd_band)
!         List of absorbers in bands

!                       Dimensions of arrays
INTEGER, INTENT(IN) ::                                                  &
    nd_band                                                             &
!         Size allocated for bands in spectral computation
  , nd_esft_term                                                        &
!         Size allocated for ESFT terms
  , nd_species
!         Size allocated for gaseous species

!                       Local variables.
INTEGER                                                                 &
    i                                                                   &
!         Loop variable
  , j                                                                   &
!         Loop variable
  , k                                                                   &
!         Loop variable
  , l                                                                   &
!         Loop variable
  , tot_subcol
!         Maximum number of sub-columns required

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MCICA_ORDER',zhook_in,zhook_handle)

  ALLOCATE(subcol_band(i_first_band:i_last_band))
  ALLOCATE(first_subcol_band(i_first_band:i_last_band+1))
  ALLOCATE(first_subcol_k(nd_band,nd_esft_term+1))

  tot_subcol=0
  subcol_band=0
  first_subcol_band(i_first_band)=1
  first_subcol_k=0

! These values may be read in from the spectral files eventually, but
! until that functionality is added, we set them here.
  DO i=i_first_band,i_last_band
    DO j=1, i_band_esft(i,index_absorb(1,i))
       subcol_band(i)=subcol_band(i)+subcol_k(i,j)
    END DO

    first_subcol_band(i+1)=first_subcol_band(i)+subcol_band(i)
    first_subcol_k(i,1)=first_subcol_band(i)

    DO j=1, i_band_esft(i,index_absorb(1,i))
      first_subcol_k(i,j+1)=first_subcol_k(i,j)+subcol_k(i,j)
    END DO

    tot_subcol=tot_subcol+subcol_band(i)
  END DO

  IF (lhook) CALL dr_hook('MCICA_ORDER',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE mcica_order
