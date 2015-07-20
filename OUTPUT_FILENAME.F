! $ID: OUTPUT_FILENAME.F V01 06/16/2012 17:24 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE OUTPUT_FILENAME WRITES OUTPUT FILE NAME DEPENDING ON YEAR,
!  MONTH, DAY, AND SPECIFIC STRING.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) YEAR      (INTEGER)  : EMISSION YEAR                              [---]
!  (2 ) MONTH     (INTEGER)  : EMISSION MONTH                             [---]
!  (3 ) DAY       (INTEGER)  : EMISSION DAY                               [---]
!  (4 ) DOMAIN    (CHARACTER): WRF DOMAIN ORDER                           [---]
!  (5 ) SFILE     (CHARACTER): TWO 12-H EMISSION DATA FILE NAME STRING    [---]
!  (6 ) OUTPUTFILE(CHARACTER): OUTPUT FILE NAME USED FOR WRF-CHEM         [---]
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/16/2012)
!******************************************************************************

      SUBROUTINE OUTPUT_FILENAME                   &
                 (YEAR, MONTH, DAY, DOMAIN, SFILE, &
                  OUTPUTFILE)

      INTEGER              :: YEAR, MONTH, DAY
      CHARACTER(LEN = 120) :: DOMAIN, CDAY, SFILE, OUTPUTFILE
!  DEFINE THE CHARACTER OF FILE NAME
      WRITE(CDAY,'(I4.4A1I2.2A1I2.2)') YEAR,'-',MONTH,'-',DAY
      PRINT*, 'CDAY : ', CDAY
      OUTPUTFILE = CDAY(1: LEN_TRIM(CDAY))//SFILE(1: LEN_TRIM(SFILE))&
                   //'_'//DOMAIN(1: LEN_TRIM(DOMAIN))

      END SUBROUTINE OUTPUT_FILENAME
