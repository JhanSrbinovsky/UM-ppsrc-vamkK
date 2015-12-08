! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Performs Source=8 field initialisation calculations

MODULE rcf_lsh_field_checks_mod

IMPLICIT NONE
!  Subroutine Rcf_Lsh_Field_Checks - Check LSH initialisation settings.

! Description:
!   Make consistency checks on the options set for initialising the large-scale
!   hydrology prognostics.

! Method:
!   If the fields for mean and standard deviation of topographic index are
!   not taken from the start dump or are interpolated from the start dump,
!   those fields which are calculated directly from them much be reset.
!   Check that variable fexp (stashcode 276) is set to be a global constant.
!   Check that none of the LSH variables fexp are set to zero of missing data.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE rcf_lsh_field_checks( data_source, fields_out, &
                                 fields_count_out, field_pos )

USE rcf_data_source_mod

USE rcf_field_type_mod, ONLY : &
    field_type

USE submodel_mod, ONLY :     &
    submodel_ident,          &
    atmos_im

USE ereport_mod, ONLY : &
    ereport

USE rcf_interp_weights_mod, ONLY : &
    h_int_active

USE UM_ParVars, ONLY :        &
    mype

USE PrintStatus_mod, ONLY :    &
    printstatus, prstatus_normal, ltimer

USE rcf_locate_mod, ONLY : &
    rcf_locate

USE rcf_stashcodes_mod, ONLY : &
    stashcode_ti_mean,         stashcode_ti_sig,          &
    stashcode_fexp,                                       &
    stashcode_gamtot,          stashcode_zw,              &
    stashcode_sthzw,           stashcode_fsat,            &
    stashcode_fwetl,                                      &
    stashcode_a_fsat,          stashcode_c_fsat,          &
    stashcode_a_fwet,          stashcode_c_fwet,          &
    stashcode_prog_sec

IMPLICIT NONE

! Arguments
TYPE( data_source_type ), POINTER    :: data_source( : )
TYPE( field_type ), POINTER          :: fields_out( : )
INTEGER, INTENT(IN)                  :: fields_count_out
INTEGER, INTENT(IN)                  :: field_pos

! Local variables
INTEGER                           :: errorstatus
INTEGER                           :: ti_mean_source
INTEGER                           :: ti_sig_source
INTEGER                           :: fexp_source
CHARACTER (LEN=*), PARAMETER      :: routinename='Rcf_Lsh_Field_Checks'
CHARACTER (LEN=80)                :: cmessage

INTEGER                           :: pos_ti_mean
INTEGER                           :: pos_ti_sig
INTEGER                           :: pos_fexp

EXTERNAL timer

! DEPENDS ON: timer
IF (ltimer) CALL timer( routinename, 3)


!-----------------------------------------------------------------------------
! Loop around source fields until a LSH field is found
!-----------------------------------------------------------------------------

!Initialise so it is clear whether they have been set in the loop below:
ti_mean_source=0
ti_sig_source=0
fexp_source=0

CALL rcf_locate ( stashcode_prog_sec, stashcode_ti_mean, fields_out,         &
          fields_count_out, pos_ti_mean,  .TRUE. )
CALL rcf_locate ( stashcode_prog_sec, stashcode_ti_sig, fields_out,          &
          fields_count_out, pos_ti_sig,  .TRUE. )
CALL rcf_locate ( stashcode_prog_sec, stashcode_fexp, fields_out,            &
          fields_count_out, pos_fexp,  .TRUE. )

IF (pos_ti_mean  /= 0 .AND. pos_ti_sig /= 0 .AND. pos_fexp /= 0 ) THEN

  ti_mean_source = data_source ( pos_ti_mean ) % source
  ti_sig_source = data_source ( pos_ti_sig ) % source
  fexp_source = data_source ( pos_fexp ) % source

ELSE

  errorstatus = 5
  cmessage = 'Could not find ti_mean and/or ti_sig and/or fexp'
  CALL ereport ( routinename, errorstatus, cmessage)

END IF


SELECT CASE( fields_out(field_pos) % stashmaster % model )

! Only atmos section:
CASE ( atmos_im )

  SELECT CASE( fields_out(field_pos) % stashmaster % section )

! Only section zero fields:
  CASE ( stashcode_prog_sec )

    SELECT CASE( fields_out(field_pos) % stashmaster % item )

! If the topographic fields or fexp are read in externally then
! recalculate the fields which are dependent on them:

    CASE( stashcode_gamtot,stashcode_a_fsat,stashcode_c_fsat,              &
          stashcode_a_fwet,stashcode_c_fwet )

      IF ( ti_mean_source /= input_dump .OR. ti_sig_source /= input_dump   &
         .OR. fexp_source /= input_dump ) THEN

        data_source( field_pos ) % source = field_calcs
        IF (printstatus >= prstatus_normal .AND. mype == 0) THEN

           WRITE(6,'(2a,i4,a)')                                            &
             'Mean &/or standard dev topographic index &/or fexp'          &
             ,' are not read from the input dump => LSH field stashno:'    &
             ,fields_out( field_pos ) % stashmaster % item                 &
             ,' is being recalculated in the reconfiguration'

        END IF

      END IF

! If the topographic fields are read in from the input dump but are then to
! be interpolated then must recalculate the fields dependent on them:

      IF (h_int_active .AND. &
          data_source( field_pos ) % source /= field_calcs) THEN

        data_source( field_pos ) % source = field_calcs
        IF (printstatus >= prstatus_normal .AND. mype == 0) THEN

          WRITE(6,'(2a,i4,3a)')                                            &
             'Mean and/or standard dev. topographic index are'             &
            ,' interpolated from the input dump => LSH stashno: '          &
            ,fields_out(field_pos) % stashmaster % item                    &
            ,' SOURCE is being reset'                                      &
            ,' to 8 and field is being ',&
            'recalculated in the reconfiguration'

        END IF

      END IF

      IF ( data_source( field_pos ) % source /= input_dump .AND.           &
        data_source( field_pos ) % source /= field_calcs ) THEN

        IF (mype == 0) THEN

        errorstatus = -10
        WRITE(cmessage,'(a,i4,2a)')'*ERROR Stashnumber: '                  &
            ,fields_out(field_pos) % stashmaster % item                    &
            ,' must be in the'                                             &
            ,' input dump or calculated in recon'
        CALL ereport ( routinename, errorstatus, cmessage)

        END IF

      END IF

    CASE( stashcode_fexp )

      IF ( data_source( field_pos ) % source /= set_to_const .AND.         &
           data_source( field_pos ) % source /= input_dump) THEN

        IF (mype == 0) THEN

          WRITE(6,'(a,i4,a)')'Stashnumber: ',stashcode_fexp                &
            ,' must be a global constant or read from the input dump.'

        END IF

        errorstatus = 10
        WRITE(cmessage,'(a,i4,a)') 'stashnumber ',stashcode_fexp           &
          ,' must be a global constant or read from the input dump.'
        CALL ereport ( routinename, errorstatus, cmessage)

      END IF

    CASE( stashcode_zw,stashcode_sthzw,stashcode_fsat,stashcode_fwetl )

      IF ( data_source( field_pos ) % source == set_to_zero .OR.           &
           data_source( field_pos ) % source == set_to_mdi ) THEN

        IF (mype == 0) THEN

          WRITE(6,'(a,i4)')'***ERROR: Stash: '                             &
            ,fields_out(field_pos) % stashmaster % item
          WRITE(6,'(a,4i4,a)')                                             &
            'It is not physically realistic to set stashnos:'              &
            ,stashcode_zw,stashcode_sthzw,stashcode_fsat,stashcode_fwetl   &
            ,' to zero or missing data.'
          WRITE(6,'(a)')'Read in from an input file or ancillary.'

        END IF

        errorstatus = 15
        WRITE(cmessage,'(i4,a)')fields_out(field_pos) % stashmaster % item &
          , ' should not be 0 or missing data'
        CALL ereport ( routinename, errorstatus, cmessage)

      END IF

    CASE DEFAULT
      errorstatus = 16
      WRITE(cmessage,'(i4,a)') fields_out(field_pos) % stashmaster % item  &
        , ' was not expected to be used to check LSH fields.'
      CALL ereport ( routinename, errorstatus, cmessage)
!--------------

    END SELECT

  END SELECT

END SELECT

! DEPENDS ON: timer
IF (ltimer) CALL timer( routinename, 4)

RETURN
END SUBROUTINE rcf_lsh_field_checks
END MODULE rcf_lsh_field_checks_mod
