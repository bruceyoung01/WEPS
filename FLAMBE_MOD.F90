! $ID: FLAMBE_MOD.F V01 06/15/2012 15:38 BRUCE EXP$
!
!******************************************************************************
!  MODULE FLAMBE_MOD READS FLAMBE FIRE EMISSION.
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
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/15/2012)
!******************************************************************************
!
      MODULE FLAMBE_MOD

      ! USE OTHER MODULES
      USE MAP_UTILS
      USE NAMELIST_MOD, ONLY : F_PM25
      USE GFED_MOD,     ONLY : FE, ER, SE
      USE GFED_MOD,     ONLY : VERTICAL_LEVEL_FIRE
      USE GFED_MOD,     ONLY : READ_NC_2D, READ_NC_3D

      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE
      TYPE(PROJ_INFO) :: PROJ

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES
      ! AND ROUTINES FROM BEING SEEN OUTSIDE 'FLAMBE_MOD.F'
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES
      PUBLIC :: FLAMBE3RD_PM25,FLAMBE3RD_OC,FLAMBE3RD_BC 

      ! ... AND THESE ROUTINES
      PUBLIC :: FLAMBE
      PUBLIC :: INIT_FLAMBE_ARRAY
      PUBLIC :: CLEANUP_FLAMBE

      ! MODULE VARIABLES
      REAL, ALLOCATABLE      :: FLAMBE3RD_PM25(:,:,:,:)
      REAL, ALLOCATABLE      :: FLAMBE3RD_OC  (:,:,:,:)
      REAL, ALLOCATABLE      :: FLAMBE3RD_BC  (:,:,:,:)

      INTEGER, PARAMETER     :: NHR = 24

      CONTAINS

!------------------------------------------------------------------------------
!
!  $ID: FLAMBE V01 07/04/2012 16:00 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE FLAMBE READS AND REWRITES FLAMBE FIRE EMISSION.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) 
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/04/2012)
!******************************************************************************
!
      SUBROUTINE FLAMBE (IDOM, IYEAR, IMONTH, IDAY)

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD,   ONLY : START_HOUR, END_HOUR
      USE NAMELIST_ARRAY_MOD,   ONLY : PROJ_CODE, DX, DY
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD,   ONLY : CORNER_LAT, CORNER_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : STAND_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : TRUELAT1, TRUELAT2
      USE NAMELIST_ARRAY_MOD,   ONLY : INJ_HEIGHT
      USE NAMELIST_ARRAY_MOD,   ONLY : FLAMBEDIR, ERDIR
      USE PARAMETER_MOD,        ONLY : KNOWNI, KNOWNJ
      USE PARAMETER_MOD,        ONLY : LATINC, LONINC
      USE NEI_MOD,              ONLY : ZFA

      ! ARGUMENTS
      INTEGER                 :: IDOM, IYEAR, IMONTH, IDAY

      ! PARAMETER
      INTEGER, PARAMETER      :: MAXL    = 50000
      REAL*8,  PARAMETER      :: HRTOSED = 3600.

      ! LOCAL VARIABLES
      CHARACTER (LEN = 255)   :: CTIME, DIR, SUF, INPTF, INPTF_ER_BC, &
                                 INPTF_ER_OC, VAR2D_ER
      CHARACTER*4             :: CYEAR
      CHARACTER*2             :: CDOM,CMONTH
      INTEGER                 :: I, ILINE, NL, IPF, IHR, HHMM
      INTEGER                 :: X1, Y1
      INTEGER                 :: IBOT, ITOP
      REAL                    :: X, Y
      REAL*8                  :: ZDIF, FRC
      REAL                    :: FRACT1,FRACT2

      ! VARIABLES FOR FLAMBE EMISSION BEFORE 2003      
      CHARACTER (LEN = 255)   :: NOUSE
      REAL                    :: TMPLAT1, TMPLON1, TMPLAT2, TMPLON2
      REAL                    :: TMPAREA, TMPFLUX
      INTEGER                 :: SAT, TMPPOSE

      ! VARIABLES FOR FLAMBE EMISSION AFTER 2004
      INTEGER                 :: INMIN, INMID, INMAX, HALF
      INTEGER*8               :: STAMP
      REAL                    :: CARBON
      REAL                    :: TMP1, TMP2, TMP3

      REAL, DIMENSION(MAXL)   :: FLAT, FLON, FAREA, FLUX
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1)  :: ER_BC, ER_OC
      INTEGER                 :: IHR1,IDAY1,IMONTH1,IYEAR1
      real, dimension(12) :: ndays ! the total day numbers
      data ndays /31,28,31,30,31,30,31,31,30,31,30,31/

      if (mod(iyear,4) == 0) then
         ndays(2) = 29
      else
         ndays(2) = 28
      endif

      ! CALL MAP_SET TO SPECIFY THE WRF USED MAP, INCLUDING GRID BOX
      CALL MAP_SET                                                   &
           (proj_code=PROJ_LC,proj=PROJ,lat1=CORNER_LAT(IDOM),       &
            lon1=CORNER_LON(IDOM), knowni=KNOWNI,knownj=KNOWNJ,      &
            dx=DX(IDOM),stdlon=STAND_LON,truelat1=TRUELAT1,          &
            truelat2=TRUELAT2)

      ! INITIALIZE FLAMBE3RD
      FLAMBE3RD_PM25 = 0.0
      FLAMBE3RD_OC = 0.0
      FLAMBE3RD_BC = 0.0
 
      ER_BC = 0.0; ER_OC = 0.0

      FRACT1 = 10.**9/(DX(IDOM)*DY(IDOM))/HRTOSED
      FRACT2 = 10.**6/(DX(IDOM)*DY(IDOM))/HRTOSED

      WRITE(CDOM  , '(I2.2)') IDOM
      WRITE(CYEAR , '(I4.4)') IYEAR
      WRITE(CMONTH, '(I2.2)') IMONTH

      ! CALL SUBROUTINE CERTICAL_LEVEL TO FIND K LEVEL OF 
      ! FLAMBE EMISSION
      CALL VERTICAL_LEVEL_FIRE(INJ_HEIGHT, IBOT, ITOP, ZDIF)

      IF (.not. F_PM25) THEN
         ! if f_PM25 = .false., then reading the emission ratios to
         ! PM25
         INPTF_ER_BC = ERDIR(1: LEN_TRIM(ERDIR))//'BC_2_PM2p5/ER'// &
                       CYEAR//CMONTH//'_d'//CDOM//'.nc'
         INPTF_ER_OC = ERDIR(1: LEN_TRIM(ERDIR))//'OC_2_PM2p5/ER'// &
                       CYEAR//CMONTH//'_d'//CDOM//'.nc'

         CALL READ_NC_2D(inptf_er_bc, SE, e_we(idom)-1, e_sn(idom)-1, &
                         ER_BC)
         CALL READ_NC_2D(inptf_er_oc, SE, e_we(idom)-1, e_sn(idom)-1, &
                         ER_OC)
      ENDIF

      ! DO HOUR LOOP
      DO IHR = START_HOUR(IDOM), END_HOUR(IDOM)
         HHMM    = IHR-1
         IDAY1   = IDAY
         IMONTH1 = IMONTH
         IYEAR1  = IYEAR
         IHR1    = HHMM
         IF (HHMM == 0) then
            IHR1 = 24
            IDAY1 = IDAY+1
            IF (IDAY1 > NDAYS(IMONTH)) THEN
               IDAY1 = 01
               IMONTH1 = IMONTH + 1
               IF (IMONTH1 > 12) THEN
                  IMONTH1 = 1
                  IYEAR1 = IYEAR + 1
               ENDIF
            ENDIF
         ENDIF

         HHMM = INT(HHMM*100)
         WRITE(CTIME, 110) IYEAR1, IMONTH1, IDAY1, HHMM
110      FORMAT(I4.4I2.2I2.2I4.4)

         !=================================================================
         ! READS FLAMBE EMISSION BEFORE 2003
         !=================================================================

         IF (IYEAR .LE. 2003) THEN
            DIR = FLAMBEDIR(1: LEN_TRIM(FLAMBEDIR))//'good'// &
                  CTIME(1: 4)//'/smoke_goes_'
            INPTF = DIR(1: LEN_TRIM(DIR))//CTIME(1: LEN_TRIM(CTIME))
            OPEN(1, FILE = INPTF, STATUS = 'OLD')
            DO ILINE = 1, MAXL
               TMPLAT1 = 0.0
               TMPLON1 = 0.0
               TMPLAT2 = 0.0
               TMPLON2 = 0.0
               SAT     = 0
               TMPAREA = 0.0
               TMPFLUX = 0.0
               TMPPOSE = 0
               READ(1, *, END = 100) TMPLAT1, TMPLON1, TMPLAT2, &
                                     TMPLON2, SAT, TMPAREA,     &
                                     TMPFLUX, TMPPOSE, NOUSE
               FLAT(ILINE)  = TMPLAT1
               FLON(ILINE)  = TMPLON1
               FAREA(ILINE) = TMPAREA
               FLUX(ILINE)  = TMPFLUX
            ENDDO
100         CONTINUE
            CLOSE(1)

         ELSE
         !=================================================================
         ! READS FLAMBE EMISSION AFTER 2004
         !=================================================================

            SUF = '.dat'
            DIR = FLAMBEDIR(1: LEN_TRIM(FLAMBEDIR))//CTIME(1: 4)// &
                  '/'//CTIME(1: 6)//'/flambe_arctas_'
            INPTF = DIR(1: LEN_TRIM(DIR))// &
                    CTIME(1: LEN_TRIM(CTIME))//SUF(1: LEN_TRIM(SUF))
            OPEN(1, FILE = INPTF, STATUS = 'OLD')
            DO ILINE = 1, MAXL
               STAMP   = 0
               TMPLON1 = 0.0
               TMPLAT1 = 0.0
               INMIN   = 0
               INMID   = 0
               INMAX   = 0
               SAT     = 0
               TMPAREA = 0.0
               TMPFLUX = 0.0
               HALF    = 0
               CARBON  = 0.0
               TMP1    = 0.0
               TMP2    = 0.0
               TMP3    = 0.0
               READ(1, *, END = 200) STAMP, TMPLON1, TMPLAT1, INMIN, &
                                     INMID, INMAX, SAT, TMPAREA,     &
                                     TMPFLUX, HALF, CARBON, TMP1,    &
                                     TMP2, TMP3
               FLAT(ILINE)  = TMPLAT1
               FLON(ILINE)  = TMPLON1
               FAREA(ILINE) = TMPAREA
               FLUX(ILINE)  = TMPFLUX
            ENDDO
200         CONTINUE
            CLOSE(1)
         ENDIF

         ! REAL # OF LINES IN THE FILES
         NL = ILINE - 1
         ! CALCULATE THE SMOKE EMISSION RATE OVER PER WRF/CHEM GRID BOX
         DO IPF = 1, NL
            CALL LATLON_TO_IJ(PROJ,FLAT(IPF),FLON(IPF),X,Y)
! for each value has been put at mid-point, e.g, at 1.5, 2.5 etc,
! instead of 1.0,2.0,...
! zf
!           IF(X .GE. 1 .AND. X .LE. E_WE(IDOM)-1 .AND.     &
!              Y .GE. 1 .AND. Y .LE. E_SN(IDOM)-1 )THEN
!           IF(X-INT(X) .GT. 0.5) X1 = INT(X+0.5)
!           IF(X-INT(X) .LE. 0.5) X1 = INT(X)
!           IF(Y-INT(Y) .GT. 0.5) Y1 = INT(Y+0.5)
!           IF(Y-INT(Y) .LE. 0.5) Y1 = INT(Y)
            IF (X .GE. 1 .AND. X .LT. E_WE(IDOM) .AND.     &
                Y .GE. 1 .AND. Y .LT. E_SN(IDOM) )THEN
               X1 = INT(X)
               Y1 = INT(Y)

               DO I = IBOT,ITOP
                  FRC = (ZFA(I+1)-ZFA(I))/ZDIF
                  ! THE FLUX UNIT IS DIFFERENT FROM < 2003 AND > 2003
                  ! < 2003 FLUX UNIT IS kg/m2/hr
                  ! > 2003 FLUX UNIT IS g/m2/hr
                  IF (IYEAR .LE.2003)THEN
                     IF (F_PM25) THEN
                        FLAMBE3RD_PM25(X1,I,Y1,IHR1) = &
                        FLAMBE3RD_PM25(X1,I,Y1,IHR1) + &
                        FLUX(IPF)*FAREA(IPF)*FRC*FRACT1
                     ELSE
                        FLAMBE3RD_OC(X1,I,Y1,IHR1) = &
                        FLAMBE3RD_OC(X1,I,Y1,IHR1) + &
                        FLUX(IPF)*ER_OC(X1,Y1)*FAREA(IPF)*FRC*FRACT1
                        FLAMBE3RD_BC(X1,I,Y1,IHR1) = &
                        FLAMBE3RD_BC(X1,I,Y1,IHR1) + &
                        FLUX(IPF)*ER_BC(X1,Y1)*FAREA(IPF)*FRC*FRACT1
                     ENDIF
                  ELSE
                     IF (F_PM25) THEN
                        FLAMBE3RD_PM25(X1,I,Y1,IHR1) = &
                        FLAMBE3RD_PM25(X1,I,Y1,IHR1) + &
                        FLUX(IPF)*FAREA(IPF)*FRC*FRACT2
                     ELSE
                        FLAMBE3RD_OC(X1,I,Y1,IHR1) = &
                        FLAMBE3RD_OC(X1,I,Y1,IHR1) + &
                        FLUX(IPF)*ER_OC(X1,Y1)*FAREA(IPF)*FRC*FRACT2
                        FLAMBE3RD_BC(X1,I,Y1,IHR1) = &
                        FLAMBE3RD_BC(X1,I,Y1,IHR1) + &
                        FLUX(IPF)*ER_BC(X1,Y1)*FAREA(IPF)*FRC*FRACT2
                     ENDIF
                  ENDIF
               ENDDO ! END OF I LOOP
            ENDIF
         ENDDO ! END OF IPF LOOP
      ENDDO ! END OF IHR LOOP

      ! END OF SUBROUTINE
      END SUBROUTINE FLAMBE


!------------------------------------------------------------------------------
!  $ID: INIT_FLAMBE_ARRAY V01 07/04/2012 17:48 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_FLAMBE_ARRAY INITIALIZES FLAMBE RELATED ARRAY.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) FLAMBE3RD_*    (REAL*8) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/04/2012)
!******************************************************************************
!
      SUBROUTINE INIT_FLAMBE_ARRAY(IDOM)

      ! REFERENCES TO F90 MODULES
      USE ERROR_MOD,            ONLY : ALLOC_ERR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT

      ! LOCAL VARIABLES
      INTEGER      :: AS, IDOM

      !========================================================================
      ! FLAMBE OUTPUT VARIABELS
      !========================================================================

      ! FLAMBE3RD_PM25
      ALLOCATE(FLAMBE3RD_PM25(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR),&
               STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('FLAMBE3RD_PM25')
      FLAMBE3RD_PM25 = 0D0

      ! FLAMBE3RD_BC
      ALLOCATE(FLAMBE3RD_BC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR),  &
               STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('FLAMBE3RD_BC')
      FLAMBE3RD_BC = 0D0

      ! FLAMBE3RD_OC
      ALLOCATE(FLAMBE3RD_OC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR),  &
               STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('FLAMBE3RD_OC')
      FLAMBE3RD_OC = 0D0

      ! RETURN TO THE CALLING ROUTINE
      END SUBROUTINE INIT_FLAMBE_ARRAY

!------------------------------------------------------------------------------
!
!  $ID: CLEANUP_FLAMBE V01 07/04/2012 17:50 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE CLEANUP_FLAMBE DEALLOCATES ALL MODULE ARRAYS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) FLAMBE3RD_*  (REAL*8) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE REFERING TO GEOS-CHEM FORMAT. (07/04/2012)
!******************************************************************************
!
      SUBROUTINE CLEANUP_FLAMBE

      !=================================================================
      ! CLEANUP_FLAMBE BEGINS HERE!
      !=================================================================
      IF ( ALLOCATED( FLAMBE3RD_PM25 ) ) DEALLOCATE( FLAMBE3RD_PM25 )
      IF ( ALLOCATED( FLAMBE3RD_OC   ) ) DEALLOCATE( FLAMBE3RD_OC   )
      IF ( ALLOCATED( FLAMBE3RD_BC   ) ) DEALLOCATE( FLAMBE3RD_BC   )

      END SUBROUTINE CLEANUP_FLAMBE

      ! END OF MODULE
      END MODULE FLAMBE_MOD
