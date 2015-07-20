! $ID: GEOMETRY_MOD.F V01 08/20/2012 09:34 BRUCE EXP$
!
!******************************************************************************
!  MODULE GEOMETRY_MOD INCLUDES SUBROUTINES RELATED TO GEOMETRY CALCULATION
!
!  SUBROUTINES:
!  ============================================================================
!  (1 )DEG_TO_DIS   :CALCULATE GREAT CIRCLE DISTANCE ON THE EARTH, 
!                               DEPENDING ON LATITUDE AND LONGITUDE.
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (08/20/2012)
!******************************************************************************
!
      MODULE GEOMETRY_MOD

      ! USE OTHER MODULES

      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES
      ! AND ROUTINES FROM BEING SEEN OUTSIDE 'GEOMETRY_MOD.F'
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES

      ! ... AND THESE FUNCTIONS AND ROURINES
      PUBLIC :: DEG_TO_DIS

      ! MODULE VARIABLES

      CONTAINS

!------------------------------------------------------------------------------
!
!  $ID: DEG_TO_DIS V01 08/20/2012 10:32 BRUCE EXP$
!
!******************************************************************************
!  FUNCTION DEG_TO_DIS CALCULATES GREAT CIRCLE DISTANCE ON THE 
!  EARTH BETWEEN TWO DIFFERENT POINTS.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (08/20/2012)
!******************************************************************************
!
      REAL FUNCTION DEG_TO_DIS(LATITUDE1, LONGITUDE1, LATITUDE2, 
     &                         LONGITUDE2)

      ! USE OTHER MODULES
      USE CONSTANTS_MOD,     ONLY : EARTH_RADIUS_M

      REAL            :: LATITUDE1, LONGITUDE1
      REAL            :: LATITUDE2, LONGITUDE2
      REAL            :: DELTALAT,  DELTALON
      REAL            :: A, C

      DELTALAT = LATITUDE2  - LATITUDE1
      DELTALON = LONGITUDE2 - LONGITUDE1
      A = (SIN(DELTALAT/2))**2+COS(LATITUDE1)*COS(LATITUDE2)*
     &    (SIN(DELTALON/2))**2
      C = 2*ASIN(MIN(1.0, SQRT(A)))
      DEG_TO_DIS = EARTH_RADIUS_M*C

      RETURN
      END FUNCTION DEG_TO_DIS


      ! END OF MODULE
      END MODULE GEOMETRY_MOD
