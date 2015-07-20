! $ID: FIRE_BG_WRF_v01.F V01 06/09/2012 16:40 BRUCE EXP$
!
!******************************************************************************
!  PROGRAM FIRE_BG_WRF_v01 READS FIRE EMISSIONS, INCLUDING FLAMBE FIRE
!  EMISSION, NCAR FINN FIRE EMISSION, NOAA GBBEP FIRE EMISSION, AND GFED
!  FIRE EMISSION, AND BACKGROUND EMISSIONS, INCLUDING EPA NEI 2005 BACKGROUND 
!  EMISSION, FINALLY REWRITE THEM INTO RADM2 EMISSION FORMAT, WHICH CAN BE
!  READ BY WRF-CHEM.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/09/2012)
!  (2 ) UNIT OF PM2.5 EMISSION IS ug/m3 m/s, FOR OTHER SPECIES, CAN REFER TO 
!       convert_emiss.exe OUTPUT FILE. (06/16/2012)

!******************************************************************************

      PROGRAM FIRE_BG_WRF_v01

      USE map_utils
      USE TOOL_KIT_MOD, ONLY : JD_to_MMDD
      USE LOGICAL_MOD
      USE NEI_MOD
      USE FLAMBE_MOD
      USE FINN_MOD
      USE GBBEP_MOD
      USE GFED_MOD

!  FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE

!**********************************************************************
!  USER NEEDS TO SPECIFY VARIABLES
!  DOMAIN    : WRF-Chem DOMAIN ORDER NUMBER
!  PROJ_CODE : PROJECTION CODE:LAMBERT CONFOMERL=1
!  STARTYEAR : START YEAR OF CALCULATION
!  ENDYEAR   : END YEAR OF CALCULATION
!  STARTJD   : START JULIAN DAY OF CALCULATION
!  ENDJD     : END JULIAN DAY OF CALCULATION
!  IX        : NUMBER OF GRID POINTS OF X DIRECTION
!  JX        : NUMBER OF GRID POINTS OF Y DIRECTION
!  HEIGHT    : HEIGHT THAT THE SMOKE CAN RISE(INJECTION HEIGHT)
!  LAT1      : LATITUDE OF THE SOUTHWEST CORNER OF WRF MOTHER DOMAIN
!  LON1      : LONGITUDE OF THE SOUTHWEST CORNER OF WRF MOTHER DOMAIN
!  DX        : GRID INCREMENT OF WEST-EAST DIRECTION(M)
!  DY        : GRID INCREMENT OF SOUTH-NORTH DIRECTION(M)
!  STDLON    : STANDARD LONGITUDE
!  TRUELAT1  : TRUE LATITUDE
!  TRUELAT2  : TRUE LATITUDE
!  NSPECIES  : NUMBER OF SPECIES WHICH FIRE EMISSION IS ADDED ON
!  SPECIES   : SPECIES WHICH FIRE EMISSION IS ADDED ON
!  NEIDIR    : NEI 2005 EMISSION DIRECTORY
!  FLAMBEDIR : FLAMBE FIRE EMISSION DIRECTORY
!  FINNDIR   : NCAR FINN FIRE EMISSION DIRECTORY
!  GBBEPDIR  : NOAA GBBEP FIRE EMISSION DIRECTORY
!  GFEDDIR   : GFED FIRE EMISSION DIRECTORY
!  ==============================START=================================

      USE NAMELIST_MOD
      CHARACTER (LEN = 10)  :: DOMAIN    = 'd03'
      INTEGER               :: PROJ_CODE =  1
      INTEGER               :: STARTYEAR =  2011
      INTEGER               :: ENDYEAR   =  2011
      INTEGER               :: STARTJD   =  152
      INTEGER               :: ENDJD     =  166
      INTEGER               :: STARTHR   =  01
      INTEGER               :: ENDHR     =  24
      INTEGER               :: IX        =  112
      INTEGER               :: JX        =  85
      INTEGER               :: KX        =  19
      REAL                  :: HEIGHT    =  650.0
      REAL                  :: LAT1      =  39.48536
      REAL                  :: LON1      = -99.43076
      REAL                  :: DX        =  4.E3
      REAL                  :: DY        =  4.E3
      REAL                  :: STDLON    = -102.452
      REAL                  :: TRUELAT1  =  35.028
      REAL                  :: TRUELAT2  =  35.028
      INTEGER               :: NSPECIES  = 2
      INTEGER, DIMENSION(NSPECIES) :: SPECIES = (/27, 29/)
      CHARACTER (LEN = 120) :: NEIDIR    = &
      '/work/jwang7/bruceyoung01/parallel/data/NEI_2005/'
      CHARACTER (LEN = 120) :: FLAMBEDIR = &
      '/work/jwang7/bruceyoung01/parallel/data/smoke/FLAMBE/'
      CHARACTER (LEN = 120) :: FINNDIR   = &
      '/work/jwang7/bruceyoung01/parallel/data/smoke/FINN/'
      CHARACTER (LEN = 120) :: GBBEPDIR  = &
      '/work/jwang7/bruceyoung01/parallel/data/smoke/GBBEP/'
      CHARACTER (LEN = 120) :: GFEDDIR   = &
      '/work/jwang7/bruceyoung01/parallel/data/smoke/GFED/'
!  ==============================END===================================
!
!  DECLARE VATIABLES FOR WRF-CHEM GRID BOX
!  ********************************************************************
!  KNOWNI      : X GRID INDEX OF THE SOUTH-WEST POINT(ORIGIN POINT)
!  KNOWNJ      : Y GRID INDEX OF THE SOUTH-WEST POINT(ORIGIN POINT)
!  ********************************************************************
      REAL                  :: KNOWNI    =  1.
      REAL                  :: KNOWNJ    =  1.
      REAL                  :: LATINC    = -999.
      REAL                  :: LONINC    = -999.

!  DECLARE VARIABLES FOR SUBROUTINE OUTPUT_FILENAME
      CHARACTER(LEN = 120)  :: CDAY
      CHARACTER(LEN = 120)  :: SFILE1 = 'wrfem_00to12z'
      CHARACTER(LEN = 120)  :: SFILE2 = 'wrfem_12to24z'
      INTEGER               :: IYEAR, IJD, IMONTH, IDAY

!  DECLARE VARIABLES FOR SUBROUTINE WRITE_RADM2
!  ********************************************************************
!  OUTPUTFILE1 : NAME OF OUTPUT FILE 00 TO 12
!  OUTPUTFILE2 : NAME OF OUTPUT FILE 12 TO 24
!  STARTT1     : START HOUR FOR OUTPUT FILE 1
!  ENDT1       : END HOUR FOR OUTPUT FILE 1
!  STARTT2     : START HOUR FOR OUTPUT FILE 2
!  ENDT2       : END HOUR FOR OUTPUT FILE 2
!  ********************************************************************
      CHARACTER(LEN = 120)  :: OUTPUTFILE1, OUTPUTFILE2
      INTEGER               :: NRADM     = 30
      INTEGER               :: TIME      = 24
      INTEGER               :: STARTT1   = 1 
      INTEGER               :: ENDT1     = 12
      INTEGER               :: STARTT2   = 13
      INTEGER               :: ENDT2     = 24
      TYPE(PROJ_INFO)       :: PROJ

!  WRF_NEI_POINT : NEI 2005 EMISSION
      REAL, DIMENSION(IX2, KX, JX2, NRADM, TIME) :: WRF_NEI_POINT
      REAL, DIMENSION(IX2, KX, JX2, TIME)        :: WRF_FLAMBE

!  CALL SUBROUTINE NEI_POINT TO READ NEI 2005 POINT EMISSION
      CALL NEI_POINT &
           (PROJ_CODE, IX, JX, KX, KP, LAT1, LON1, DX, DY, KNOWNI, &
            KNOWNJ, LATINC, LONINC, STDLON, TRUELAT1, TRUELAT2,    &
            NEIDIR,                                                &
            WRF_NEI_POINT)

!  DO THE JULIAN DAY LOOP TO READ FIRE EMISSION IN SPECIFIC DAY
      DO IYEAR = STARTYEAR, ENDYEAR
       DO IJD = STARTJD, ENDJD

!  CALL SUBROUTINE JD_TO_MMDD WHICH IS USED TO CONVERT 
!  JULIAN DAY TO CALENDAR DAY
        CALL JD_TO_MMDD(IJD,IYEAR,IMONTH,IDAY)

!  CALL SUBROUTINE OUTPUT_FILENAME TO SPECIFY OUTPUT FILE NAME
        CALL OUTPUT_FILENAME                                  &
            (IYEAR, IMONTH, IDAY, DOMAIN, SFILE1, OUTPUTFILE1)
        CALL OUTPUT_FILENAME                                  &
            (IYEAR, IMONTH, IDAY, DOMAIN, SFILE2, OUTPUTFILE2)

!  CALL SUBROUTINE WRITE_SEPA_V01 TO WRITE SPECIFIED EMISSION INTO
!  WRFCHEM REQUIRED FORMAT
        CALL WRITE_RADM2                                      &
            (OUTPUTFILE1, STARTT1, ENDT1, TIME, IX2, KX, JX2, &
             NSPECIES, SPECIES, NRADM, ENAME,                 &
             WRF_NEI_POINT, WRF_FLAMBE)
        CALL WRITE_RADM2                                      &
            (OUTPUTFILE2, STARTT2, ENDT2, TIME, IX2, KX, JX2, &
             NSPECIES, SPECIES, NRADM, ENAME,                 &
             WRF_NEI_POINT, WRF_FLAMBE)
       END DO ! IJD
      END DO ! IYEAR

      END PROGRAM FIRE_BG_WRF_v01
