! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Reads in the pre-computed variables for use in the
!+ interaction between UKCA aerosols and the radiation
!+ code.
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
SUBROUTINE ukca_radaer_read_precalc(                                    &
                                    ierr,                               &
                                    cmessage                            &
)

USE parkind1, ONLY: jpim, jprb 
USE yomhook,  ONLY: lhook, dr_hook
USE ukca_radaer_precalc
USE ukca_radaer_lut

USE spcrg3a_mod, ONLY: ip_solar, ip_infra_red
IMPLICIT NONE

!
! Arguments
!
!
! Error indicator
!
INTEGER ierr
!
! Error message
!
CHARACTER (len=80) :: cmessage

!
! Local variables
!
CHARACTER (len=80) :: filename
INTEGER ios


CHARACTER (len=9), DIMENSION(npd_ukca_spectrum), PARAMETER :: &
    spectrum_name = (/ 'shortwave', 'longwave ' /)

!
! 
!

INTEGER n_band

INTEGER aerosol_index

!
! Number of bands and AOD wavelengths in the spectral files.
! These must match the values read from the precalc file.
!
INTEGER, DIMENSION(npd_ukca_lut_spectrum) :: specf_n_band
INTEGER specf_n_aod_wavel

INTEGER n_x
INTEGER n_nr
INTEGER n_ni

INTEGER, PARAMETER :: funit = 165

CHARACTER (len=80) :: line ! buffer for neglected lines (headers,...)

!
! Loop indices
!
INTEGER i, &
        j, &
        k, &
        n

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_in, zhook_handle)

ierr = 0

!
! Get the number of wavebands as defined in the spectral files
!
! DEPENDS ON: ukca_get_specinfo
CALL ukca_get_specinfo(  &
  specf_n_band(1),       &
  specf_n_band(2),       &
  specf_n_aod_wavel      &
  )

!
! Get the name of the file to open and open it.
!
CALL get_file(funit, filename, 80, ios)
IF (ios /= 0) THEN
  ierr = -1
  cmessage = 'Error reading name of UKCA pre-computed file.'
  IF (lhook) CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_out, zhook_handle)
  RETURN
END IF

OPEN(unit=funit, file=filename, iostat=ios)
IF (ios /= 0) THEN
  ierr = -1
  cmessage = 'Error opening ' // filename
  IF (lhook) CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_out, zhook_handle)
  RETURN
END IF

!
! Read the formatted file.
!
!
! Number of integration points and weight of each point in the
! integration.
!
READ(funit, '(29x,i3)') precalc%n_integ_pts
IF (precalc%n_integ_pts > npd_integ_pts) THEN
  ierr = -1
  cmessage = 'Increase the maximum number of integration points.'
  CLOSE(funit)
  IF (lhook) CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_out, zhook_handle)
  RETURN
END IF

READ(funit, '(a)') line ! column header

DO n = 1, precalc%n_integ_pts
    
  READ(funit, '(4x,E12.5)') precalc%weight(n)
    
END DO ! n

!
! For each spectrum
!
DO i = 1, npd_ukca_spectrum
  
  READ(funit, '(a)') line ! spectrum header
  
  !
  ! Number of wavebands
  !
  READ(funit, '(20x,i3)') n_band
  IF (n_band /= specf_n_band(i)) THEN
    ierr = -1
    cmessage = 'UKCA file is not compatible with the ' // &
               spectrum_name(i) // ' spectral file.'
    CLOSE(funit)
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF
  
  IF (n_band > npd_ukca_band) THEN
    ierr = -1
    cmessage = &
      'Too many wavebands. Increase the maximum number of wavebands.'
    CLOSE(funit)
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF
  
  !
  ! Dimension of the UKCA LUTs for each aerosol mode.
  !
  DO j = 1, npd_ukca_lut_mode
  
    READ(funit, '(37x,3(i3,1x))') n_x, n_nr, n_ni
    IF ( n_x /= ukca_lut(j, i)%n_x  .or.   &
        n_nr /= ukca_lut(j, i)%n_nr .or.   &
        n_ni /= ukca_lut(j, i)%n_ni) THEN
      ierr = -1
      cmessage = &
        'UKCA file is not compatible with UKCA look-up tables.'
      CLOSE(funit)
      IF (lhook) THEN
        CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_out, zhook_handle)
      END IF
      RETURN
    END IF
    
  END DO ! j (UKCA aerosol mode)
  
  !
  ! Loop on all wavebands
  !
  DO j = 1, n_band
  
    READ(funit, '(a)') line ! Waveband header.
    READ(funit, '(a)') line ! column header
    
    !
    ! Wavelength at integration point and monochromatic solar
    ! irradiance or Planckian at this wavelength
    !
    DO n = 1, precalc%n_integ_pts
      
      READ(funit, '(6x,2(E12.5,1x))') &
        precalc%wavelength(n, j, i), precalc%irrad(n, j, i)
            
    END DO ! n (integration points)
    
    !
    ! Waveband-integrated irradiance or Planckian
    !
    READ(funit, '(27x,E12.5)') precalc%flux(j, i)
    
    !
    ! Complex refractive index of aerosols at integration points.
    !
    READ(funit, '(a)') line ! header
    
    DO k = 1, npd_ukca_maxcomptype
      
      READ(funit, '(24x,i2)') aerosol_index
      
      READ(funit, '(a)') line ! column header
      
      DO n = 1, precalc%n_integ_pts
        
        READ(funit, '(8x,2(E12.5,1x))') &
          precalc%realrefr(aerosol_index, n, j, i), &
          precalc%imagrefr(aerosol_index, n, j, i)
                
      END DO ! n
      
    END DO ! k
    
  END DO ! j (wavebands)
  
END DO ! i (spectrum)

!
! Number and values of wavelengths for the aerosol optical depth
!
READ(funit, '(a)') line ! header

READ(funit, '(22x,i3)') precalc%n_ukca_aod_wavel
IF (precalc%n_ukca_aod_wavel /= specf_n_aod_wavel) THEN
  ierr = -1
  cmessage = 'UKCA file is not compatible with the ' // &
             'longwave spectral file (AOD block).'
  CLOSE(funit)
  IF (lhook) CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_out, zhook_handle)
  RETURN
END IF
IF (precalc%n_ukca_aod_wavel > npd_ukca_aod_wavel) THEN
  ierr = -1
  cmessage = 'Increase npd_ukca_aod_wavel in UKCAPCALC.'
  CLOSE(funit)
  IF (lhook) CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_out, zhook_handle)
  RETURN
END IF

READ(funit, '(a)') line ! column header

DO i = 1, precalc%n_ukca_aod_wavel

  READ(funit, '(4x,E12.5)') precalc%aod_wavel(i)
  
END DO ! i

!
! Complex refractive index of aerosols at AOD wavelengths.
!
READ(funit, '(a)') line ! header

DO k = 1, npd_ukca_maxcomptype

  READ(funit, '(24x,i2)') aerosol_index
  
  READ(funit, '(a)') line ! column header
  
  DO n = 1, precalc%n_ukca_aod_wavel

    READ(funit, '(8x,2(E12.5,1x))') &
      precalc%aod_realrefr(aerosol_index, n), &
      precalc%aod_imagrefr(aerosol_index, n)
  
  END DO ! n

END DO ! k

CLOSE(funit)

IF (lhook) CALL dr_hook('UKCA_RADAER_READ_PRECALC', zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_read_precalc
