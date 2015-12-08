! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine interface:
      subroutine dataw(rows,columns,fieldsize,nlevels,levn,len_extra,   &
     & fieldn,len1_lookup_all,lookup_all,fixhd,                         &
     & len_cfi, cfi1, cfi2, cfi3, fldsizelev,ftin1,ftout,               &
     & tracer_grid,add_wrap_pts,ibm_to_cray,compress,rmdi_input,wave,   &
     & lsmask, l_bit_32,                                                &
     & icode)

        USE um_types, ONLY: real32
        USE mask_compression, ONLY: compress_to_mask
        USE Decomp_DB, ONLY : decompose
        IMPLICIT NONE

!
! Description: This writes the data out using WRITFLD.
!               If compress oa_pack is used.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Subroutine arguments
!   Scalar arguments with intent(in):

      integer rows       ! number of rows in input pp field
      integer columns    ! number of columns in input pp field
      integer fieldsize  ! number of points in output anc. field
      integer nlevels    ! number of levels
      integer levn       ! current level number
      integer len_extra
      integer ftin1
      integer ftout
      integer fieldn
      integer len1_lookup_all
      integer icode          ! error status

      real rmdi_input

      logical tracer_grid
      logical add_wrap_pts
      logical ibm_to_cray
      logical compress
      logical l_bit_32

!   Array  arguments with intent(in):

      integer lookup_all(len1_lookup_all,*)
      integer fixhd(*)

      integer len_cfi(3)         ! dimensions of arrays
      integer cfi1(len_cfi(1))   !   compressed
      integer cfi2(len_cfi(2))   !   field index
      integer cfi3(len_cfi(3))   !   arrays
      integer fldsizelev(nlevels)   ! size of output field on each level

      logical lsmask(rows*columns)

! local arrays
      real(KIND=real32) datain(rows*columns)
      real data_field(rows*columns)
      real field_wrap(columns+2,rows)
      real field_to_write(fieldsize)
      real extra_data(len_extra+1)  ! space for extra data

! Local Scalars
      integer i
      integer j,istart,iend,ii
      integer field_type          ! 0 for tracers; 1 for velocities
      integer no_cmp   ! # of pts in full compressed field (all levels)
      integer no_rows_m   ! number of rows east-west on model grid
      integer n_sea_points

      logical LTimer      ! timer switch (set to false)
      logical cyclic_grid ! T => input field to OA_PACK has
                          !      overlap points
      logical wave ! creating wave dump

      CHARACTER(LEN=256) cmessage      ! error message


! Input arguments for decompose_smexe
      INTEGER                                                           &

     &  tot_levels        ! IN  :total number of levels


! Function & Subroutine calls:
      REAL P1,P2
      LOGICAL LNER
      LNER(P1,P2) = ((ABS(P1-P2))  <   (1.E-6*ABS(P1+P2)))

!- End of header


!L  1. Read data and do number format conversion if needed
! DEPENDS ON: readdata
      Call readdata( rows, columns, ftin1, ibm_to_cray, len_extra,      &
     &               l_bit_32, .FALSE., data_field, extra_data )


!L 1.1 Convert real missing data indicators
      if ( rmdi_input  /=  rmdi) then
        i=0
        do j = 1,rows*columns
          if ( LNER (data_field(j), rmdi_input) ) then
!         if ( rmdi_input  >   0.0 ) then
             data_field(j) = rmdi
             i=i+1
          end if
        end do
        if (i >  0) then
        write (6,*) i,' RMDI converted from ',rmdi_input,' to ',rmdi
        endif
      end if

!L 2. Add in wrap points when add_wrap_pts=t

      if (add_wrap_pts) then

        do j=1,rows
          do i=1,columns
            field_wrap(i,j)=data_field(i+(j-1)*columns)
          END DO
        END DO
        do j=1,rows
          field_wrap(columns+1,j)=field_wrap(1,j)
          field_wrap(columns+2,j)=field_wrap(2,j)
        END DO

!L 3. Pack data using compression indices when compress=t

       if (compress) then

         if(.not.wave) then

           if (tracer_grid) then
            field_type = 0
            no_rows_m = rows
           else
            field_type = 1
            no_rows_m = rows + 1
           end if

           no_cmp = 0
           do i = 1, nlevels   ! do not use levn in this loop
             no_cmp = no_cmp + fldsizelev(i)
           end do

           cyclic_grid = .TRUE.   ! input pp fields do not
                                   ! have wrap-points

           LTimer = .FALSE.
           icode  = 0

! DEPENDS ON: oa_pack
           call OA_PACK(icode, cmessage, LTimer,                        &
     &       no_rows_m, columns+2, nlevels, len_cfi(1), fieldsize,      &
     &       cfi1, cfi2, cfi3, no_cmp, rmdi,                            &
     &       levn, field_type, cyclic_grid, field_wrap,                 &
     &       field_to_write)


           if (icode  >   0) then
             write (6,*) 'error from OA_PACK:', cmessage
             go to 9999
           end if

         else       ! add_wrap .and. compress .and. wave

! compress using SLMASK for wave model - use SEA POINTS set to TRUE
! a value for n-SEA-points is returned from this subroutine

!!!!!!!!! This needs attention

            CALL compress_to_mask(data_field,field_to_write,lsmask,       &
                rows*columns,n_SEA_points)

            print*,'after to land points no_cmp is ',n_sea_points
            no_cmp=n_sea_points

         endif

        else      ! add_wrap .and. .not. compress

          do j=1,rows
            do i=1,columns+2
              field_to_write(i+(j-1)*(columns+2) ) =  field_wrap(i,j)
            END DO
          END DO

        endif

      else         ! .not. add_wrap

!L 3.1 Pack data using compression indices when compress=t

        if (compress) then

          if(.not. wave) then

           if (tracer_grid) then
            field_type = 0
            no_rows_m = rows
           else
            field_type = 1
            no_rows_m = rows + 1
           end if

           no_cmp = 0
           do i = 1, nlevels   ! do not use levn in this loop
             no_cmp = no_cmp + fldsizelev(i)
           end do

           cyclic_grid = .FALSE.   ! input pp fields do not
                                   ! have wrap-points

           LTimer = .FALSE.
           icode  = 0

! DEPENDS ON: oa_pack
           call OA_PACK(icode, cmessage, LTimer,                        &
     &       no_rows_m, columns, nlevels, len_cfi(1), fieldsize,        &
     &       cfi1, cfi2, cfi3, no_cmp, rmdi,                            &
     &       levn, field_type, cyclic_grid, data_field, field_to_write)


           if (icode  >   0) then
             write (6,*) 'error from OA_PACK:', cmessage
             go to 9999
           end if

          else         ! .not. add_wrap .and. compress .and. wave

! compress using SLMASK for wave model - use SEA POINTS set to TRUE
! a value for n-SEA-points is returned from this subroutine

            CALL compress_to_mask(data_field,field_to_write,lsmask,       &
                rows*columns,n_SEA_points)

            print*,'after to land points no_cmp is ',n_sea_points
            no_cmp=n_sea_points

          endif

         else         ! .not. add_wrap .and. .not. compress

          do j = 1, fieldsize
            field_to_write(j) =  data_field(j)
          end do

         endif

      end if

!L 5. Output data using WRITFLDS

!C TEMP print out data LSMASK for wave dump

      if(lookup_all(23,fieldn) == 38) then
        write(6,*) ' '
        print*,'before writing data array'
        istart=1
        iend=istart+columns-1
        do i=rows,1,-1
          print*, (field_to_write(ii),ii=istart,iend)
          istart=istart+columns
          iend=iend+columns
        enddo
      endif

      CALL decompose(columns, rows,0,0,1)

! DEPENDS ON: writflds
      CALL WRITFLDS(ftout,1,fieldn,lookup_all,                          &
     &              len1_lookup_all,field_to_write,fieldsize,           &
     &              fixhd,                                              &
     &              icode,cmessage )

      if (icode  >   0) then
        write (6,*) 'error from WRITFLDS:', cmessage
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE dataw
! Purpose: Works out the lookup tables for the dump/ancillary    *
!           file header from the pp fields                       *
!
! Subroutine interface:

!+ Skip namelists in f90 compiled UM code removing need for
!+ assign -f 77 g:sf  in script
!
! Subroutine Interface:
