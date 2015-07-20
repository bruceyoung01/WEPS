!  $ID: FINN_MOD.F V01 07/04/2012 15:42 BRUCE EXP$
!
!******************************************************************************
!  MODULE FINN_MOD READS AND COMPUTES FINN EMISSION AND REWRITES IT
!  INTO THE FORMAT WHICH CAN BE READ BY convert_emiss.exe IN WRF-Chem.
!  
!  
!  MODULE VARIABLES:
!  ============================================================================
!  (1 )
!
!  MODULE ROUTINES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/04/2012)
!******************************************************************************
!
      MODULE FINN_MOD


      ! USE OTHER MODULES
      USE MAP_UTILS
      USE NEI_MOD,        ONLY : NHR
      USE NAMELIST_MOD,   ONLY : F_PM25
      USE GFED_MOD,       ONLY : VERTICAL_LEVEL_FIRE

      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE
      TYPE(PROJ_INFO) :: PROJ

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES 
      ! AND ROUTINES FROM BEING SEEN OUTSIDE "FINN_MOD.F"
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES
      PUBLIC   :: FINN3RD_OC, FINN3RD_BC, FINN3RD_PM25

      ! ... EXCEPT THESE ROUTINES
      PUBLIC   :: FINN
      PUBLIC   :: INIT_FINN_ARRAY
      PUBLIC   :: CLEANUP_FINN

      !=================================================================
      ! PARAMETERS
      !=================================================================


      REAL, ALLOCATABLE      :: FINN3RD_OC(:,:,:,:)
      REAL, ALLOCATABLE      :: FINN3RD_BC(:,:,:,:)
      REAL, ALLOCATABLE      :: FINN3RD_PM25(:,:,:,:)



      ! ... EXCEPT THESE VARIABLES ...

      CONTAINS


!------------------------------------------------------------------------------
!
!  $ID: FINN V01 07/23/2012 10:31 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE FINN READS AND REWRITES FINN FIRE EMISSION.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) 
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/23/2012)
!******************************************************************************
!
      SUBROUTINE FINN (IDOM, IYEAR, IMONTH, IJD)

      !REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD,   ONLY : START_HOUR, END_HOUR
      USE NAMELIST_ARRAY_MOD,   ONLY : PROJ_CODE, DX, DY
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD,   ONLY : CORNER_LAT, CORNER_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : STAND_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : TRUELAT1, TRUELAT2
      USE NAMELIST_ARRAY_MOD,   ONLY : INJ_HEIGHT
      USE NAMELIST_ARRAY_MOD,   ONLY : FINNDIR
      USE PARAMETER_MOD,        ONLY : KNOWNI, KNOWNJ
      USE PARAMETER_MOD,        ONLY : LATINC, LONINC
      USE NEI_MOD,              ONLY : ZFA

      ! ARGUMENTS
      INTEGER                :: IDOM, IYEAR, IMONTH, IJD

      ! PARAMETER
      INTEGER, PARAMETER     :: MAXL     = 500000
      INTEGER, PARAMETER     :: NMOZART  = 40
      INTEGER, PARAMETER     :: NSAPRC   = 41
      REAL,    PARAMETER     :: DAYTOSED = 3600.

      ! LOCAL VARIABLES
      CHARACTER (LEN = 255)  :: CYEAR, CMONTH, DIR, INPTF, HEADER
      INTEGER                :: X1, Y1
      INTEGER                :: IBOT, ITOP
      REAL                   :: X, Y
      REAL*8                 :: ZDIF, FRC

      ! FINN FILE VARIABLES
      CHARACTER (LEN = 10)   :: TMPTIME, TMPHR, CHR
      CHARACTER (LEN = 255), DIMENSION(MAXL) :: FINNTIME
      INTEGER                :: ILINE, NL, TMPJD, TMPGENVEG
      INTEGER                :: I, IHR, IPF
      INTEGER, DIMENSION(MAXL) :: FINNJD
      REAL                   :: TMPLATI, TMPLONGI, TMPAREA
      REAL                   :: TMPCO2, TMPCO, TMPH2
      REAL                   :: TMPNO, TMPNO2, TMPSO2, TMPNH3, TMPCH4
      REAL                   :: TMPNMOC, TMPBIGALD, TMPBIGALK, TMPBIGENE
      REAL                   :: TMPC10H16, TMPC2H4, TMPC2H5OH, TMPC2H6
      REAL                   :: TMPC3H6, TMPC3H8, TMPCH2O, TMPCH3CHO
      REAL                   :: TMPCH3COCH3, TMPCH3COCHO, TMPCH3COOH
      REAL                   :: TMPCH3OH, TMPCRESOL, TMPGLYALD, TMPHYAC
      REAL                   :: TMPISOP, TMPMACR, TMPMEK, TMPMVK, TMPHCN
      REAL                   :: TMPCH3CN, TMPTOLUENE, TMPPM25, TMPOC
      REAL                   :: TMPBC, TMPPM10, TMPHCOOH, TMPC2H2

      REAL, DIMENSION(MAXL)  :: FINNLAT, FINNLON, FINNAREA
      REAL, DIMENSION(MAXL, NMOZART) :: FLUX
      REAL :: FRACT
      INTEGER :: ICC


      ! CALL MAP_SET TO SPECIFY THE WRF USED MAP, INCLUDING GRID BOX
      CALL MAP_SET   &
             (proj_code=PROJ_LC,proj=PROJ,lat1=CORNER_LAT(IDOM),lon1=CORNER_LON(IDOM),    &
              knowni=KNOWNI,knownj=KNOWNJ,dx=DX(IDOM),stdlon=STAND_LON,truelat1=TRUELAT1, &
              truelat2=TRUELAT2)

      ! INITIALIZE FINN3RD
      FINN3RD_OC = 0.0
      FINN3RD_BC = 0.0
      FINN3RD_PM25 = 0.0

      ! SPECIFY DIRECTORY FOR FINN FIRE EMISSION
      WRITE(CYEAR,  110) IYEAR
      WRITE(CMONTH, 120) IMONTH
110   FORMAT(I4.4)
120   FORMAT(I2.2)
      DIR=FINNDIR(1: LEN_TRIM(FINNDIR))//CYEAR(1: LEN_TRIM(CYEAR))//'/'
      INPTF = DIR(1: LEN_TRIM(DIR))//'GLOB_MOZ4_'//                     &
              CYEAR(1: LEN_TRIM(CYEAR))//CMONTH(1: LEN_TRIM(CMONTH))//  &
              '.txt'

      ! OPEN FINN FIRE EMISSION FILE
      OPEN(21, FILE = trim(INPTF), STATUS = 'OLD')

      ! READ HEADER
      READ(21, *) HEADER

      ! DO LINE LOOP TO READ ALL THE DATA INTO DIFFERENT ARRAY
      DO ILINE = 1, MAXL
         READ(21, *, END = 100)TMPJD, TMPTIME, TMPGENVEG, TMPLATI,         &
                               TMPLONGI, TMPAREA,TMPCO2, TMPCO, TMPH2,     &
                               TMPNO, TMPNO2, TMPSO2, TMPNH3, TMPCH4,      &
                               TMPNMOC, TMPBIGALD, TMPBIGALK, TMPBIGENE,   &
                               TMPC10H16, TMPC2H4, TMPC2H5OH, TMPC2H6,     &
                               TMPC3H6, TMPC3H8, TMPCH2O, TMPCH3CHO,       &
                               TMPCH3COCH3, TMPCH3COCHO, TMPCH3COOH,       &
                               TMPCH3OH, TMPCRESOL, TMPGLYALD, TMPHYAC,    &
                               TMPISOP, TMPMACR, TMPMEK, TMPMVK, TMPHCN,   &
                               TMPCH3CN, TMPTOLUENE, TMPPM25, TMPOC,       &
                               TMPBC, TMPPM10, TMPHCOOH, TMPC2H2
         FINNJD(ILINE)   = TMPJD
         FINNTIME(ILINE) = TMPTIME
         FINNLAT(ILINE)  = TMPLATI
         FINNLON(ILINE)  = TMPLONGI
         FINNAREA(ILINE) = TMPAREA
         FLUX(ILINE, :)  =(/TMPCO2, TMPCO, TMPH2,                      &
                            TMPNO, TMPNO2, TMPSO2, TMPNH3, TMPCH4,     &
                            TMPNMOC, TMPBIGALD, TMPBIGALK, TMPBIGENE,  &
                            TMPC10H16, TMPC2H4, TMPC2H5OH, TMPC2H6,    &
                            TMPC3H6, TMPC3H8, TMPCH2O, TMPCH3CHO,      &
                            TMPCH3COCH3, TMPCH3COCHO, TMPCH3COOH,      &
                            TMPCH3OH, TMPCRESOL, TMPGLYALD, TMPHYAC,   &
                            TMPISOP, TMPMACR, TMPMEK, TMPMVK, TMPHCN,  &
                            TMPCH3CN, TMPTOLUENE, TMPPM25, TMPOC,      &
                            TMPBC, TMPPM10, TMPHCOOH, TMPC2H2/)
      ENDDO
100   CONTINUE
      CLOSE(21)

      ! CALL SUBROUTINE CERTICAL_LEVEL TO FIND K LEVEL OF 
      ! FLAMBE EMISSION
      CALL VERTICAL_LEVEL_FIRE(INJ_HEIGHT, IBOT, ITOP, ZDIF)

      ! GET THE REAL LINES OF THE FILE
      NL = ILINE - 1
      FRACT = 10.**9/(DX(IDOM)*DY(IDOM))/DAYTOSED
      DO IPF = 1, NL
         ICC = 0
         IF (FINNJD(IPF) == IJD) THEN
            IF (LEN_TRIM(FINNTIME(IPF)) .EQ. 3) THEN
               CHR   = ADJUSTL(FINNTIME(IPF))
               TMPHR = CHR(1:1)
               READ (TMPHR, '(I1.1)') IHR
               ICC = 1
            ELSE IF (LEN_TRIM(FINNTIME(IPF)) .EQ. 4) THEN
               CHR   = ADJUSTL(FINNTIME(IPF))
               TMPHR = CHR(1:2)
               READ (TMPHR, '(I2.2)') IHR
               ICC = 1
            ENDIF
         ELSE IF (FINNJD(IPF) == IJD+1) THEN
          ! for ihr = 24
            IF (LEN_TRIM(FINNTIME(IPF)) .LT. 3) THEN
               IHR = 24
               ICC = 1
            ENDIF
         ENDIF

         IF (ICC == 0) CYCLE

         CALL LATLON_TO_IJ(PROJ,FINNLAT(IPF),FINNLON(IPF),X,Y)

! for each value has been put at mid-point, e.g, at 1.5, 2.5 etc,
! instead of 1.0,2.0,...
! zf
!         IF (X .GE. 1 .AND. X .LE. E_WE(IDOM)-1 .AND.    &
!             Y .GE. 1 .AND. Y .LE. E_SN(IDOM)-1 )THEN
!             IF(X-INT(X) .GT. 0.5) X1 = INT(X+0.5)
!             IF(X-INT(X) .LE. 0.5) X1 = INT(X)
!             IF(Y-INT(Y) .GT. 0.5) Y1 = INT(Y+0.5)
!             IF(Y-INT(Y) .LE. 0.5) Y1 = INT(Y)
         IF (X .GE. 1 .AND. X .LT. E_WE(IDOM) .AND.     &
             Y .GE. 1 .AND. Y .LT. E_SN(IDOM) )THEN
            X1 = INT(X)
            Y1 = INT(Y)

            DO I = IBOT,ITOP
               FRC = (ZFA(I+1)-ZFA(I))/ZDIF

               ! FLUX UNIT IS kg/day for each FINNAREA, to ug/m2/sec
               IF (F_PM25) THEN
               FINN3RD_PM25(X1,I,Y1,IHR)=FINN3RD_PM25(X1,I,Y1,IHR) +  &
                                       FLUX(IPF, 35)*FRC*FRACT
               ELSE
               FINN3RD_OC(X1,I,Y1,IHR)=FINN3RD_OC(X1,I,Y1,IHR) +  &
                                       FLUX(IPF, 36)*FRC*FRACT

               FINN3RD_BC(X1,I,Y1,IHR)=FINN3RD_BC(X1,I,Y1,IHR) +  &
                                       FLUX(IPF, 37)*FRC*FRACT
               ENDIF
            ENDDO ! END OF I LOOP
         ENDIF
      ENDDO ! END OF IPF LOOP

130   FORMAT(I1.1)
140   FORMAT(I2.2)

      END SUBROUTINE FINN

!------------------------------------------------------------------------------
!  $ID: INIT_FINN_ARRAY V01 07/03/2012 09:39 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_FINN_ARRAY INITIALIZES FINN RELATED ARRAY.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) FINN3RD_OC    (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) FINN3RD_BC    (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) FINN3RD_PM25  (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/10/2012)
!******************************************************************************
!
   SUBROUTINE INIT_FINN_ARRAY(IDOM)

      ! REFERENCES TO F90 MODULES
      USE ERROR_MOD,            ONLY : ALLOC_ERR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT

      ! LOCAL VARIABLES
      INTEGER      :: AS, IDOM

      !========================================================================
      ! FINN OUTPUT VARIABELS
      !========================================================================

      ! FINN3RD_OC
      ALLOCATE(FINN3RD_OC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR),STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('FINN3RD_OC')
      FINN3RD_OC = 0D0

      ! FINN3RD_BC
      ALLOCATE(FINN3RD_BC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR),STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('FINN3RD_BC')
      FINN3RD_BC = 0D0

      ! FINN3RD_PM25
      ALLOCATE(FINN3RD_PM25(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR),STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('FINN3RD_PM25')
      FINN3RD_PM25 = 0D0

      ! RETURN TO THE CALLING ROUTINE
   END SUBROUTINE INIT_FINN_ARRAY


!------------------------------------------------------------------------------
!
!  $ID: CLEANUP_FINN V01 06/18/2012 08:48 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE CLEANUP_FINN DEALLOCATES ALL MODULE ARRAYS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) FINN3RD_OC (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) FINN3RD_BC (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) FINN3RD_PM25 (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE REFERING TO GEOS-CHEM FORMAT. (07/10/2012)
!******************************************************************************
!
   SUBROUTINE CLEANUP_FINN

      !=================================================================
      ! CLEANUP_FINN BEGINS HERE!
      !=================================================================
      IF ( ALLOCATED( FINN3RD_OC) ) DEALLOCATE( FINN3RD_OC   )
      IF ( ALLOCATED( FINN3RD_BC) ) DEALLOCATE( FINN3RD_BC   )
      IF ( ALLOCATED( FINN3RD_PM25) ) DEALLOCATE( FINN3RD_PM25   )

   END SUBROUTINE CLEANUP_FINN

END MODULE FINN_MOD
