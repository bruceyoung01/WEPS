!  $ID: WEPS.F V01 06/28/2012 10:47 BRUCE EXP$
!
!******************************************************************************
!
!  MODEL NAME     : WEPS
!  MODEL VERSION  : V01
!  UPDATED TIME   : 06/28/2012 10:49 CST
!  ============================================================================
!  PROGRAM WEPS IS THE MAIN PROGRAM OF WEPS, WHICH HAS THE CAPACITY AS
!  FOLLOWS:
!  (1 ) SET UP DOMAIN, SAME AS WPS
!  (2 ) READ BACKGROUND EMISSIONS AND FIRE EMISSIONS
!  (3 ) GRID EMISSIONS INTO WPS DOMAIN
!  (4 ) SPECIFY DIFFERENT EMISSION SPECIES INTO RADMS SPECIES
!
!  MODULES REFERENCED BY WEPS.F
!  ============================================================================
!  (1 ) READ_NAMELIST           : NAMELIST_MOD.F
!  (2 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE, REFERING TO RICHARD'S ART-RTM
!       MODEL. (06/28/2012)
!  (2 ) ADD QFED. (09/30/2013)
!******************************************************************************

      PROGRAM WEPS

      ! REFERENCED F90 MODULES
      USE TOOL_KIT_MOD,         ONLY : ANNUAL_JD,      JD_TO_MMDD
      USE NAMELIST_ARRAY_MOD,   ONLY : START_YEAR,     END_YEAR
      USE NAMELIST_ARRAY_MOD,   ONLY : START_MONTH,    END_MONTH
      USE NAMELIST_ARRAY_MOD,   ONLY : START_DAY,      END_DAY
      USE NAMELIST_ARRAY_MOD,   ONLY : MAX_N,          MAX_DOM
      USE NAMELIST_ARRAY_MOD,   ONLY : CORNER_LON,     CORNER_LAT
      USE NAMELIST_ARRAY_MOD,   ONLY : LNEI,           LINTEX
      USE NAMELIST_ARRAY_MOD,   ONLY : LFLAMBE,        LFINN
      USE NAMELIST_ARRAY_MOD,   ONLY : LGBBEP,         LGFED, LGFAS
      USE NAMELIST_ARRAY_MOD,   ONLY : LSEVIRI,        LQFED
      USE NAMELIST_MOD,         ONLY : READ_NAMELIST
      USE GRID_MOD,             ONLY : MK_WRF_GRIDS
      USE NEI_MOD,              ONLY : INIT_NEI_ARRAY
      USE NEI_MOD,              ONLY : CLEANUP_NEI
      USE NEI_MOD,              ONLY : NEI
      USE INTEX_MOD,            ONLY : INIT_INTEX_ARRAY
      USE INTEX_MOD,            ONLY : INTEX
      USE INTEX_MOD,            ONLY : CLEANUP_INTEX
      USE FLAMBE_MOD,           ONLY : INIT_FLAMBE_ARRAY
      USE FLAMBE_MOD,           ONLY : CLEANUP_FLAMBE
      USE FLAMBE_MOD,           ONLY : FLAMBE
      USE FINN_MOD,             ONLY : INIT_FINN_ARRAY
      USE FINN_MOD,             ONLY : FINN
      USE FINN_MOD,             ONLY : CLEANUP_FINN
      USE GBBEP_MOD,            ONLY : INIT_GBBEP_ARRAY
      USE GBBEP_MOD,            ONLY : CLEANUP_GBBEP
      USE GBBEP_MOD,            ONLY : GBBEP
      USE GFED_MOD,             ONLY : INIT_GFED_ARRAY
      USE GFED_MOD,             ONLY : GFED
      USE GFED_MOD,             ONLY : CLEANUP_GFED
      USE GFAS_MOD,             ONLY : INIT_GFAS_ARRAY
      USE GFAS_MOD,             ONLY : GFAS
      USE GFAS_MOD,             ONLY : CLEANUP_GFAS
      USE QFED_MOD,             ONLY : INIT_QFED_ARRAY
      USE QFED_MOD,             ONLY : QFED
      USE QFED_MOD,             ONLY : CLEANUP_QFED
      USE SEVIRI_MOD,           ONLY : INIT_SEVIRI_ARRAY
      USE SEVIRI_MOD,           ONLY : SEVIRI
      USE SEVIRI_MOD,           ONLY : CLEANUP_SEVIRI
      USE WRITE_MOD,            ONLY : WRITE_RADM2

      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE

      ! PARAMETERS
      INTEGER, PARAMETER          :: STARTT1 = 1
      INTEGER, PARAMETER          :: ENDT1   = 12
      INTEGER, PARAMETER          :: STARTT2 = 13
      INTEGER, PARAMETER          :: ENDT2   = 24

      ! LOCAL VARIABLES
      CHARACTER (LEN = 255)       :: SYSTEM_DATE, SYSTEM_TIME
      LOGICAL                     :: ILOGIC
      INTEGER                     :: IDOM
      INTEGER                     :: IJD
      INTEGER                     :: IYEAR, IMONTH, IDAY
      INTEGER                     :: PRE_MONTH
      INTEGER, DIMENSION(MAX_N)   :: START_JD, END_JD

      !=================================================================
      !  WEPS BEGINS HERE!
      !=================================================================

      ! READ VARIABLES FROM namelist.weps
      CALL READ_NAMELIST

      DO IDOM = 1, MAX_DOM

       !================================================================
       ! ECHO SYSTEM TIME
       !================================================================

       CALL DATE_AND_TIME(SYSTEM_DATE, SYSTEM_TIME)
       WRITE(6, 90  ) 'CURRENT TIME : ', SYSTEM_DATE, '-',     &
                      SYSTEM_TIME(1:2), ':', SYSTEM_TIME(3:4), &
                      ':', SYSTEM_TIME(5:6)
       WRITE(6, 100 ) REPEAT('-', 79)
       WRITE(6, 110 ) 'Now Processing Domain : ', IDOM
90     FORMAT(/, A, A8, A, A2, A, A2, A, A2/)
100    FORMAT(/, A  )
110    FORMAT(A, I2 )

       !================================================================
       ! INITIALIZATIONS...
       !================================================================

       ! LOGICAL FOR READING MONTHLY AND 3-HOURLY DATA OR NOT
       ! IF .TRUE., READ MONTHLY AND 3-HOURLY DATA, OTHERWISE, NOT.
       ILOGIC = .TRUE.

       ! INITIALIZE NEI ARRAYS
       IF (LNEI(IDOM)) THEN
        CALL INIT_NEI_ARRAY(IDOM)
       ENDIF

       ! INITIALIZE INTEX ARRAYS
       IF (LINTEX(IDOM)) THEN
        CALL INIT_INTEX_ARRAY(IDOM)
       ENDIF

       ! INITIALIZE FLAMBE ARRAYS
       IF (LFLAMBE(IDOM)) THEN
        CALL INIT_FLAMBE_ARRAY(IDOM)
       ENDIF

       ! INITIALIZE FINN ARRAYS
       IF (LFINN(IDOM)) THEN
        CALL INIT_FINN_ARRAY(IDOM)
       ENDIF

       ! INITIALIZE GBBEP ARRAYS
       IF (LGBBEP(IDOM)) THEN
        CALL INIT_GBBEP_ARRAY(IDOM)
       ENDIF

       ! INITIALIZE GFED ARRAYS
       IF (LGFED(IDOM)) THEN
        CALL INIT_GFED_ARRAY(IDOM)
       ENDIF

       ! INITIALIZE SEVIRI ARRAYS
       IF (LSEVIRI(IDOM)) THEN
        CALL INIT_SEVIRI_ARRAY(IDOM)
       ENDIF

       ! INITIALIZE GFAS ARRAYS
       IF (LGFAS(IDOM)) THEN
        CALL INIT_GFAS_ARRAY(IDOM)
       ENDIF

       ! INITIALIZE QFED ARRAYS
       IF (LQFED(IDOM)) THEN
        CALL INIT_QFED_ARRAY(IDOM)
       ENDIF

       !================================================================
       ! GET RELEVANT GEOGRAPHY INFO
       !================================================================

       ! CALCULATE WRF GRID INFO
       CALL MK_WRF_GRIDS (IDOM)

       !================================================================
       ! BACKGROUND EMISSIONS PROCESSING BEGINS HERE!
       !================================================================       

       ! NEI_2005
       IF (LNEI(IDOM)) THEN
        CALL NEI(IDOM)
       ENDIF

       ! INTEX
       IF (LINTEX(IDOM)) THEN
        CALL INTEX(IDOM)
       ENDIF

       !================================================================
       ! FIRE EMISSION PROCESSING BEGINS HERE!
       !================================================================

       PRE_MONTH = 0

       ! CONVERT CALENDAR DATE TO JULIAN DATE
       START_JD(IDOM) = ANNUAL_JD(START_YEAR(IDOM), START_MONTH(IDOM), &
                            START_DAY(IDOM))
       END_JD(IDOM)   = ANNUAL_JD(END_YEAR(IDOM),   END_MONTH(IDOM),   &
                            END_DAY(IDOM))
       
       DO IYEAR = START_YEAR(IDOM), END_YEAR(IDOM)
        DO IJD = START_JD(IDOM), END_JD(IDOM)


         ! CONVERT JULIAN DATE TO CALENDAR DATE
         CALL JD_TO_MMDD(IJD, IYEAR, IMONTH, IDAY)

         IF (PRE_MONTH .NE. IMONTH) THEN
            ILOGIC = .TRUE.
         ELSE
            ILOGIC = .FALSE.
         ENDIF
         PRE_MONTH = IMONTH

         ! GFED FIRE EMISSION
         IF (LGFED(IDOM)) THEN
            CALL GFED(ILOGIC, IDOM, IYEAR, IMONTH, IDAY)
         ENDIF

         !==============================================================
         ! USING GFED 3-HOURLY EMISSION DATA TO CALCULATE 3-HOURLY 
         ! EMISSION RATIO (BC/PM2.5, OC/PM2.5, BC/TPM, OC/TPM, AND 
         ! PM2.5/TPM)
         !==============================================================

         ! FLAMBE FIRE EMISSION
         IF (LFLAMBE(IDOM)) THEN
            CALL FLAMBE(IDOM, IYEAR, IMONTH, IDAY)
         ENDIF

         ! FINN FIRE EMISSION
         IF (LFINN(IDOM)) THEN
            CALL FINN(IDOM, IYEAR, IMONTH, IJD)
         ENDIF

         ! GBBEP FIRE EMISSION
         IF (LGBBEP(IDOM)) THEN
            CALL GBBEP(IDOM, IYEAR, IJD)
         ENDIF

         ! SEVIRI FIRE EMISSION
         IF (LSEVIRI(IDOM)) THEN
            CALL SEVIRI(IDOM, IYEAR, IMONTH, IDAY)
         ENDIF

         ! GFAS FIRE EMISSION
         IF (LGFAS(IDOM)) THEN
            CALL GFAS(IDOM, IYEAR, IMONTH, IDAY)
         ENDIF

         ! QFED FIRE EMISSION
         IF (LQFED(IDOM)) THEN
            CALL QFED(IDOM, IYEAR, IMONTH, IDAY)
         ENDIF

         ! WRITE BACKGROUND AND FIRE EMISSION INTO WRF-Chem FORMAT
         CALL WRITE_RADM2(IDOM, IYEAR, IMONTH, IDAY, STARTT1, ENDT1)
         CALL WRITE_RADM2(IDOM, IYEAR, IMONTH, IDAY, STARTT2, ENDT2)
        ENDDO ! IJD
       ENDDO ! IYEAR

       !================================================================
       ! CLEANUP BEGINS HERE!
       !================================================================

       ! CLEANUP_NEI
       CALL CLEANUP_NEI

       ! CLEANUP_INTEX
       CALL CLEANUP_INTEX

       ! CLEANUP_FLAMBE
       CALL CLEANUP_FLAMBE

       ! CLEANUP_FINN
       CALL CLEANUP_FINN

       ! CLEANUP GBBEP
       CALL CLEANUP_GBBEP

       ! CLEANUP_GFED
       CALL CLEANUP_GFED

       ! CLEANUP_SEVIRI
       CALL CLEANUP_SEVIRI

       ! CLEANUP_GFAS
       CALL CLEANUP_GFAS

       ! CLEANUP_QFED
       CALL CLEANUP_QFED

      ENDDO ! IDOM

      !================================================================
      ! ECHO SYSTEM TIME
      !================================================================

      CALL DATE_AND_TIME(SYSTEM_DATE, SYSTEM_TIME)
      WRITE(6, 90  ) 'CURRENT TIME : ', SYSTEM_DATE, '-',      &
                     SYSTEM_TIME(1:2), ':', SYSTEM_TIME(3:4),  &
                     ':', SYSTEM_TIME(5:6)

      WRITE(6, *) 'SUCCESSFUL COMPLETION' 
      ! END OF THE MAIN PROGRAM
      END PROGRAM WEPS
