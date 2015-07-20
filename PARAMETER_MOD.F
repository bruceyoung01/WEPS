! $ID: PARAMETER_MOD.F V01 07/01/2012 09:37 BRUCE EXP$
!
!******************************************************************************
!  MODULE PARAMETER_MOD INCLUDES PARAMETERS USED IN WEPS PROGRAM.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) KNOWNI      (REAL*8) : X-LOCATION OF KNOWN LAT/LON               [---]
!  (2 ) KNOWNJ      (REAL*8) : Y-LOCATION OF KNOWN LAT/LON               [---]
!  (3 ) LATINC      (REAL*8) : LATITUDE INCREMENT FOR CYLINDRICAL LAT/LON[---]
!  (4 ) LONINC      (REAL*8) : LONGITUDE INCREMENT FOR CYLINDRICAL
!                              LAT/LON, ALSO THE LON INCREMENT FOR
!                              GAUSSIAN GRID                             [---]
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/01/2012)
!******************************************************************************
!
      MODULE PARAMETER_MOD

      ! MAKE EVERYTHING PUBLIC
      PUBLIC

      ! PARAMETER FOR MAP_SET
      REAL, PARAMETER      :: KNOWNI =  1.0
      REAL, PARAMETER      :: KNOWNJ =  1.0
      REAL, PARAMETER      :: LATINC = -999.0
      REAL, PARAMETER      :: LONINC = -999.0

      END MODULE PARAMETER_MOD
