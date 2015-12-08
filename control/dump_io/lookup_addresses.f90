! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Parameters for addressing elements in a UM file lookup table
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Dump I/O

MODULE lookup_addresses
  IMPLICIT NONE

  ! Validity time
  INTEGER, PARAMETER :: lbyr   =1   ! Year
  INTEGER, PARAMETER :: lbmon  =2   ! Month
  INTEGER, PARAMETER :: lbdat  =3   ! Day of month
  INTEGER, PARAMETER :: lbhr   =4   ! Hour
  INTEGER, PARAMETER :: lbmin  =5   ! Minute
  INTEGER, PARAMETER :: lbsec  =6   ! Second, from LBREL=3
  ! LBREL<3      Integer, Parameter :: LBDAY  =6   ! Day number

  ! Data time
  INTEGER, PARAMETER :: lbyrd  =7   ! Year
  INTEGER, PARAMETER :: lbmond =8   ! Month
  INTEGER, PARAMETER :: lbdatd =9   ! Day of month
  INTEGER, PARAMETER :: lbhrd  =10  ! Hour
  INTEGER, PARAMETER :: lbmind =11  ! Minute
  INTEGER, PARAMETER :: lbsecd =12  ! Second, from LBREL=3
  ! LBREL<3      Integer, Parameter :: LBDAYD =12  ! Day number

  INTEGER, PARAMETER :: lbtim  =13  ! Time indicator
  INTEGER, PARAMETER :: lbft   =14  ! Forcast period (hours)
  INTEGER, PARAMETER :: lblrec =15  ! Length of data record
  INTEGER, PARAMETER :: lbcode =16  ! Grid type code
  INTEGER, PARAMETER :: lbhem  =17  ! Hemisphere indicator
  INTEGER, PARAMETER :: lbrow  =18  ! Number of rows in grid
  INTEGER, PARAMETER :: lbnpt  =19  ! Number of points per row
  INTEGER, PARAMETER :: lbext  =20  ! Length of extra data
  INTEGER, PARAMETER :: lbpack =21  ! Packing method indicator
  INTEGER, PARAMETER :: lbrel  =22  ! Header release number
  INTEGER, PARAMETER :: lbfc   =23  ! Field code
  INTEGER, PARAMETER :: lbcfc  =24  ! Second field code
  INTEGER, PARAMETER :: lbproc =25  ! Processing code
  INTEGER, PARAMETER :: lbvc   =26  ! Vertical coordinate type
  INTEGER, PARAMETER :: lbrvc  =27  ! Coordinate type for reference level
  INTEGER, PARAMETER :: lbexp  =28  ! Experiment number
  INTEGER, PARAMETER :: lbegin =29  ! Start record
  INTEGER, PARAMETER :: lbnrec =30  ! No of records for direct access files
                                    ! - OR -
                                    ! Length on disk of record
                                    ! (in words) inclusive of any padding
  INTEGER, PARAMETER :: lbproj =31  ! Met-O-8 projection number
  INTEGER, PARAMETER :: lbtyp  =32  ! Met-O-8 field type
  INTEGER, PARAMETER :: lblev  =33  ! Met-O-8 level code
  INTEGER, PARAMETER :: lbrsvd1=34  ! Reserved for future PP-package use
  INTEGER, PARAMETER :: lbrsvd2=35  ! Reserved for future PP-package use
  INTEGER, PARAMETER :: lbrsvd3=36  ! Reserved for future PP-package use
  INTEGER, PARAMETER :: lbrsvd4=37  ! Used externally with PP-package as
                                    ! LBENS ensemble member number
  INTEGER, PARAMETER :: lbsrce =38  ! =1111 to indicate following apply to UM
  INTEGER, PARAMETER :: data_type=39! Indicator for real/int or timeseries
  INTEGER, PARAMETER :: naddr  =40  ! Start address in DATA_REAL or DATA_INT
  INTEGER, PARAMETER :: lbuser3=41  ! Free for user-defined function
  INTEGER, PARAMETER :: item_code=42!Stash code
  INTEGER, PARAMETER :: lbplev =43  ! Pseudo-level indicator (if defined)
  INTEGER, PARAMETER :: lbuser6=44  ! Free for user-defined function
  INTEGER, PARAMETER :: model_code=45 ! internal model identifier

  INTEGER, PARAMETER :: bulev  =46  ! Upper level boundary
  INTEGER, PARAMETER :: bhulev =47  ! Upper level boundary
  INTEGER, PARAMETER :: brsvd3 =48  ! Reserved for future PP-package use
  INTEGER, PARAMETER :: brsvd4 =49  ! Reserved for future PP-package use
  INTEGER, PARAMETER :: bdatum =50  ! Datum value
  INTEGER, PARAMETER :: bacc   =51  ! (Packed fields) Packing accuracy
  INTEGER, PARAMETER :: blev   =52  ! Level
  INTEGER, PARAMETER :: brlev  =53  ! Lower level boundary
  INTEGER, PARAMETER :: bhlev  =54  ! (Hybrid levels) A-level of value
  INTEGER, PARAMETER :: bhrlev =55  ! Lower level boundary
  INTEGER, PARAMETER :: bplat  =56  ! Real latitude of 'pseudo' N Pole
  INTEGER, PARAMETER :: bplon  =57  ! Real longitude of 'pseudo' N Pole
  INTEGER, PARAMETER :: bgor   =58  ! Grid orientation
  INTEGER, PARAMETER :: bzy    =59  ! Zeroth latitude
  INTEGER, PARAMETER :: bdy    =60  ! Latitude interval
  INTEGER, PARAMETER :: bzx    =61  ! Zeroth longitude
  INTEGER, PARAMETER :: bdx    =62  ! Longitude interval
  INTEGER, PARAMETER :: bmdi   =63  ! Missing data indicator
  INTEGER, PARAMETER :: bmks   =64  ! M,K,S scaling factor

  INTEGER, PARAMETER :: lbcc_lbyr   = 1  ! Year
  INTEGER, PARAMETER :: lbcc_lbmon  = 2  ! Month
  INTEGER, PARAMETER :: lbcc_lbdat  = 3  ! Day of the month
  INTEGER, PARAMETER :: lbcc_lbhr   = 4  ! Hour
  INTEGER, PARAMETER :: lbcc_lbmin  = 5  ! Minute
  INTEGER, PARAMETER :: lbcc_lbsec  = 6  ! Second, from LBREL=3
  ! LBREL<3      Integer, Parameter :: LBCC_LBDAY  = 6  ! Day number

  INTEGER, PARAMETER :: lbcc_lbegin = 7  ! Start record
  INTEGER, PARAMETER :: lbcc_naddr  = 8  ! Start address of DATA
  ! Mapping of MPP_LOOKUP; analogous to mapping in PP header

  INTEGER, PARAMETER :: p_naddr=1    ! Address on local PE
  INTEGER, PARAMETER :: p_lblrec=2   ! Local length of record

END MODULE lookup_addresses
