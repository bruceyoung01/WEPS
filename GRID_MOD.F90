!  $ID: GRID_MOD.F90 V01 01/13/2014 11:50 BRUCE EXP$
!
!******************************************************************************
!  MODULE GRID_MOD PROCESSES STUFF RELATED TO GRID BOX, LIKE GRID CENTER 
!  LATITUDE, GRID CENTER LONGITUDE, GRID CORNER LATITUDE, AND GRID CORNER 
!  LONGITUDE.
!
!  MODULE VARIABLES:
!  ============================================================================
!  (1 )
!
!
!  MODULE SUBROUTINES:
!  ============================================================================
!  (1 ) MK_WRF_GRIDS : CALCULATES AND WRITES OUT WRF GRID CENTER LATITUDE, 
!                      GRID CENTER LONGITUDE, GRID CORNER LATITUDE, 
!                      AND GRID CORNER LONGITUDE.
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (01/13/2014)
!******************************************************************************
!
      MODULE GRID_MOD

      ! USE OTHER MODULES
      USE MAP_UTILS

      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE
      TYPE(PROJ_INFO)  :: PROJ

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES
      ! AND ROUTINES FROM BEING SEEN OUTSIDE 'GRID_MOD.F90'
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES

      ! ... AND THESE ROUTINES
      PUBLIC :: MK_WRF_GRIDS

      ! MODULE VARIABLES

      !=================================================================
      ! MODULE ROUTINES -- FOLLOW BELOW THE "CONTAINS" STATEMENT
      !=================================================================

      CONTAINS

!------------------------------------------------------------------------------
!
!  $ID: MK_WRF_GRIDS V01 01/13/2014 13:15 BRUCE EXP$
!
!******************************************************************************
!     SUBTOURINE MK_WRF_GRIDS WRITES OUT GRID BOX LATITUDE/LONGITUDE INFO INTO 
!     TEXT FILES.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (01/13/2014)
!******************************************************************************
!
      SUBROUTINE MK_WRF_GRIDS (IDOM)

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : PROJ_CODE, DX, DY
      USE NAMELIST_ARRAY_MOD, ONLY : E_WE, E_SN
      USE NAMELIST_ARRAY_MOD, ONLY : CORNER_LAT, CORNER_LON
      USE NAMELIST_ARRAY_MOD, ONLY : STAND_LON
      USE NAMELIST_ARRAY_MOD, ONLY : TRUELAT1, TRUELAT2
      USE NAMELIST_ARRAY_MOD, ONLY : OUTPUTDIR
      USE PARAMETER_MOD,      ONLY : KNOWNI, KNOWNJ
      USE PARAMETER_MOD,      ONLY : LATINC, LONINC

      ! ARGUMENTS
      INTEGER                      :: IDOM

      ! LOCAL VARIABLES
      CHARACTER*255                :: FLNM
      CHARACTER*2                  :: CDOM
      INTEGER                      :: I, J, M4

      REAL                         :: REALI, REALJ
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1)    :: WLAT, WLON
      REAL, DIMENSION(E_WE(IDOM), E_SN(IDOM))        :: LAT,  LON
      REAL, DIMENSION(4, E_WE(IDOM)-1, E_SN(IDOM)-1) :: CORNER_WLAT
      REAL, DIMENSION(4, E_WE(IDOM)-1, E_SN(IDOM)-1) :: CORNER_WLON

      WRITE(CDOM, '(I2.2)') IDOM

      ! CALL MAP_SET TO SPECIFY THE WRF USED MAP, INCLUDING GRID BOX
      CALL MAP_SET                     &
            ( PROJ_LC, PROJ, CORNER_LAT(IDOM), CORNER_LON(IDOM),      &
              CORNER_LAT(IDOM), CORNER_LON(IDOM), KNOWNI, KNOWNJ,     &
              DX(IDOM), DY(IDOM), LATINC, LONINC, STAND_LON,          &
              TRUELAT1, TRUELAT2)

      ! CALL SUBROUTINE IJ_TO_LATLON TO GET THE LATITUDE AND LONGITUDE
      ! OF EACH GRID BOX IN WRFChem
      WLAT = 0.0
      WLON = 0.0

      DO J = 1, E_SN(IDOM) - 1
         DO I = 1, E_WE(IDOM) - 1
            REALI = REAL(I) + 0.5
            REALJ = REAL(J) + 0.5
            CALL IJ_TO_LATLON(PROJ, REALI, REALJ, &
                              WLAT(I, J), WLON(I, J))
         ENDDO
      ENDDO

      ! TO GET CORNER POINTS INFORTION
      LAT         = 0.0
      LON         = 0.0
      CORNER_WLAT = 0.0
      CORNER_WLON = 0.0

      DO J = 1, E_SN(IDOM)
         DO I = 1, E_WE(IDOM)
            REALI = REAL(I)
            REALJ = REAL(J)
            CALL IJ_TO_LATLON(PROJ, REALI, REALJ, &
                              LAT(I, J), LON(I, J))
            IF (LAT(I, J) .EQ. 0.0 .OR. &
                LON(I, J) .EQ. 0.0) THEN
               WRITE(6, *) 'I = ', I, 'J = ', J, &
                           'LAT = ', LAT(I, J),  &
                           'LON = ', LON(I, J)
            ENDIF
         ENDDO
      ENDDO

      DO J = 1, E_SN(IDOM) - 1
         DO I = 1, E_WE(IDOM) - 1
            CORNER_WLAT(1, I, J) = LAT(I, J)
            CORNER_WLON(1, I, J) = LON(I, J)
            CORNER_WLAT(2, I, J) = LAT(I + 1, J)
            CORNER_WLON(2, I, J) = LON(I + 1, J)
            CORNER_WLAT(3, I, J) = LAT(I + 1, J + 1)
            CORNER_WLON(3, I, J) = LON(I + 1, J + 1)
            CORNER_WLAT(4, I, J) = LAT(I, J + 1)
            CORNER_WLON(4, I, J) = LON(I, J + 1)
         ENDDO
      ENDDO

      ! WRITE OUT THE GRID LATITUDE AND LONGITUDE TO TEXT FILES

      ! WRITE GRID CENTER LATITUDE
      FLNM = TRIM(OUTPUTDIR)//'grids/lat_d'//CDOM//'.txt'
      OPEN (30, FILE = TRIM(FLNM))
      DO J = 1, E_SN(IDOM) - 1
         WRITE(30, 100) (WLAT(I, J), I = 1, E_WE(IDOM) - 1)
      ENDDO
      CLOSE(30)

      ! WRITE GRID CENTER LONGITUDE
      FLNM = TRIM(OUTPUTDIR)//'grids/lon_d'//CDOM//'.txt'
      OPEN (30, FILE = TRIM(FLNM))
      DO J = 1, E_SN(IDOM) - 1
         WRITE(30, 100) (WLON(I, J), I = 1, E_WE(IDOM) - 1)
      ENDDO
      CLOSE(30)

      ! WRITE GRID CORNER LATITUDE
      FLNM = TRIM(OUTPUTDIR)//'grids/corner_lat_d'//CDOM//'.txt'
      OPEN (30, FILE = TRIM(FLNM))
      DO J = 1, E_SN(IDOM) - 1
         DO I = 1, E_WE(IDOM) - 1
            WRITE(30, 101) (CORNER_WLAT(M4, I, J), M4 = 1, 4)
         ENDDO
      ENDDO
      CLOSE(30)


      ! WRITE GRID CORNER LONGITUDE
      FLNM = TRIM(OUTPUTDIR)//'grids/corner_lon_d'//CDOM//'.txt'
      OPEN (30, FILE = TRIM(FLNM))
      DO J = 1, E_SN(IDOM) - 1
         DO I = 1, E_WE(IDOM) - 1
            WRITE(30, 101) (CORNER_WLON(M4, I, J), M4 = 1, 4)
         ENDDO
      ENDDO
      CLOSE(30)

100   FORMAT(720F12.5)
101   FORMAT(4F12.5)

      END SUBROUTINE MK_WRF_GRIDS

      END MODULE GRID_MOD
