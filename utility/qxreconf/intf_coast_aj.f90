! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: INTF_COAST_AJ------------------------------------------
!LL
!LL  Purpose: To calculate a suitable value of NSEARCH for call to
!LL           SPIRAL_S and then call SPIRAL_S
!L
!LL  Reviewer:  Date of review:
!LL
!LL  Tested under compiler: cft77
!LL  Tested under OS version: UNICOS 7
!LL
!LL  Code version no: 1       Date: 17 November 1993
!LL
!LL  Modification History:
!LL
!LL  Model
!LL  Version    Date
!LL  4.0        26/07/95  Make local copy of land sea mask to work
!LL                       to correct error of mask not being reset
!LL                       after call to SPIRAL_S.
!LL  4.0        02/11/95  The Fortran for squaring an array needed to be
!LL       redefined for the DecAlpha because of lexcon. (N.Farnon)
!    4.5        30/07/97  Skip put of loop when nearest neighbour found.
!LL  5.1      22/6/00      Jump out of loop if max. search radius
!LL                        is reached. P. Selwood.
!LL
!LL Programming standard: UM Doc Paper 3, version
!LL
!LL Logucal component number:
!LL
!LL Project task: S1
!LL
!LL
!LL Documentation:
!LL   UMDP S1
!LL
!LL -------------------------------------------------------------
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration
      SUBROUTINE INTF_COAST_AJ                                          &
     &           (LAND_SEA_MASK,INDEX_UNRES,NO_POINT_UNRES,             &
     &           POINTS_PHI,POINTS_LAMBDA,DATA_FIELD,SEA_LAND,          &
     &           CYCLIC,MAXDIM)

      IMPLICIT NONE

!*L ARGUMENTS:---------------------------------------------------

      INTEGER                                                           &
     & POINTS_PHI                                                       &
                        !IN number of rows in grid
     &,POINTS_LAMBDA                                                    &
                        !IN number of columns in grid
     &,NO_POINT_UNRES                                                   &
                        !INOUT number of unresolved points
     &,LAND_SEA_MASK(POINTS_LAMBDA*POINTS_PHI)                          &
                        !IN land sea mask
     &,INDEX_UNRES(POINTS_LAMBDA*POINTS_PHI)                            &
                        !INOUT index to unresolved pts
     &,SEA_LAND         !IN =0 for sea field  =1/-1 for land field

      REAL                                                              &
     & DATA_FIELD(POINTS_LAMBDA*POINTS_PHI) !IN field

      LOGICAL                                                           &
     & CYCLIC           ! IN =T if data covers complete latitude circle

!*L LOCAL VARIABLES

      INTEGER                                                           &
     & I,J,JJ,JJJ,K,JK,NSCH                                             &
                            ! indices
     &,IPT,IROW,ICOL                                                    &
                      ! coordinate of unresolved point
     &,IPOINT,IUNRES                                                    &
                      ! do loop variable
     &,NPOINTS                                                          &
                      ! number of points in serach box
     &,MAXDIM                                                           &
                     ! largest dimension of field
     &,IR((1+2*MAXDIM)*(1+2*MAXDIM))                                    &
                         ! row numbers of points to serach
     &,IC((1+2*MAXDIM)*(1+2*MAXDIM))                                    &
                         ! col numbers of points to search
     &,IND_SEARCH((1+2*MAXDIM)*(1+2*MAXDIM))                            &
                         ! index to points to search
                     ! still unresolved after calling this subroutine
     &,ISUM_MASK                                                        &
                     ! number of surrounding points which have data
     &,ISEARCH                                                          &
                     ! largest dimension of field
     &,NSEARCH                                                          &
                     ! minimum search radius required.
     &,LAND_SEA_TEMP(POINTS_LAMBDA*POINTS_PHI) ! local copy of mask


!*L   EXTERNAL ROUTINES
      EXTERNAL SPIRAL_S

      NSEARCH=0

! Take local copy of land sea mask to work with
      DO IPOINT=1,POINTS_LAMBDA*POINTS_PHI
        LAND_SEA_TEMP(IPOINT)=LAND_SEA_MASK(IPOINT)
      ENDDO

! toggle land sea mask to exclude unresolved points from meaning process
      DO IUNRES=1,NO_POINT_UNRES
        IF(SEA_LAND == 0)THEN
          LAND_SEA_TEMP(INDEX_UNRES(IUNRES))=1
        ELSE
          LAND_SEA_TEMP(INDEX_UNRES(IUNRES))=0
        ENDIF
      ENDDO


! Loop around unresolved points
      DO IUNRES=1,NO_POINT_UNRES

! find unresolved point coordinate in terms of rows and cols
        IPT=INDEX_UNRES(IUNRES)
        IROW=INT(REAL(IPT)/REAL(POINTS_LAMBDA)+1)
        ICOL=IPT-(IROW-1)*POINTS_LAMBDA

        ISEARCH=MAXDIM
        DO I=1,MAXDIM

! Calculate number of points in search box
          NPOINTS=(1+2*I)**2 ! number of grid points in search box

! calculate surrounding points' coords in terms of rows and cols
          JJJ=1
          DO J=-I,I
            DO JJ=JJJ,JJJ+2*I
              IR(JJ)=IROW+J
            ENDDO
            JJJ=JJJ+1+2*I
          ENDDO

          JJJ=1+2*I
          JK=1
          DO J=-I,I
            DO JJ=0,2*I
              IC(JK+JJ*JJJ)=ICOL+J
            ENDDO
            JK=JK+1
          ENDDO

! Check that col and rows are in range of grid
          DO IPOINT=1,NPOINTS
            IF(IC(IPOINT) >  POINTS_LAMBDA) THEN
              IF(CYCLIC) THEN
                IC(IPOINT)=IC(IPOINT)-POINTS_LAMBDA
              ELSE
                IC(IPOINT)=POINTS_LAMBDA
              ENDIF
            ENDIF
            IF(IC(IPOINT) <  1) THEN
              IF(CYCLIC) THEN
                IC(IPOINT)=IC(IPOINT)+POINTS_LAMBDA
              ELSE
                IC(IPOINT)=1
              ENDIF
            ENDIF
            IF(IR(IPOINT) <  1) IR(IPOINT)=1
            IF(IR(IPOINT) >  POINTS_PHI) IR(IPOINT)=POINTS_PHI
          ENDDO

! Form index search array
          DO IPOINT=1,NPOINTS
            IND_SEARCH(IPOINT)=(IR(IPOINT)-1)*POINTS_LAMBDA+IC(IPOINT)
          ENDDO

! search for data around this point. If no data is found the point
! remains unresolved

          ISUM_MASK=0   ! number of points with data found

          DO IPOINT=1,NPOINTS
           IF(IABS(LAND_SEA_TEMP(IND_SEARCH(IPOINT))) == IABS(SEA_LAND) &
     &       .AND.DATA_FIELD(IND_SEARCH(IPOINT)) >= 0.0)THEN
             ISUM_MASK=ISUM_MASK+1
           ENDIF
         ENDDO

          IF(ISUM_MASK >  0)THEN
            ISEARCH=MIN0(ISEARCH,I)
            GOTO 100
          END IF

        ENDDO

 100    CONTINUE
        NSEARCH=MAX0(ISEARCH,NSEARCH)

! If the search radius reaches the maximum, we have no need to
! do testing for the rest of the unresolved points
        IF (NSEARCH >= MAXDIM) EXIT
      ENDDO


! DEPENDS ON: spiral_s
      CALL SPIRAL_S(LAND_SEA_TEMP,INDEX_UNRES,NO_POINT_UNRES,           &
     &                 POINTS_PHI,POINTS_LAMBDA,DATA_FIELD,NSEARCH,     &
     &                 SEA_LAND,CYCLIC)


      RETURN
      END SUBROUTINE INTF_COAST_AJ
