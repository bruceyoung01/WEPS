! $ID: LOGICAL_MOD.F V01 06/09/2012 22:40 BRUCE EXP$
!
!******************************************************************************
!  MODULE LOGICAL_MOD CONTAINS ALL OF THE LOGICAL SWITCHES USED BY 
!  PROGRAM FIRE_BG_WRF_v01.F
!
!  MODULE VARIABLES:
!  ============================================================================
!  (1 ) LNEIPOINT  (LOGICAL) :ON/OFF SWITCH FOR NEI 2005 POINT EMISSION
!  (2 ) LNEIAREA   (LOGICAL) :ON/OFF SWITCH FOR NEI 2005 AREA EMISSION
!  (3 ) LFLAMBE    (LOGICAL) :ON/OFF SWITCH FOR FLAMBE FIRE EMISSION
!  (4 ) LFINN      (LOGICAL) :ON/OFF SWITCH FOR NCAR FINN FIRE EMISSION
!  (5 ) LGBBEP     (LOGICAL) :ON/OFF SWITCH FOR NOAA GBBEP FIRE EMISSION
!  (6 ) LGFED      (LOGICAL) :ON/OFF SWITCH FOR GFED FIRE EMISSION
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/09/2012)
!******************************************************************************

      MODULE LOGICAL_MOD

      IMPLICIT NONE

!  ====================================================================
!  DEFINE MODULE VARIABLES
!  ====================================================================

!  NEI 2005
      LOGICAL :: LNEIPOINT
      LOGICAL :: LNEIAREA

!  FLAMBE
      LOGICAL :: LFLAMBE

!  NCAR FINN
      LOGICAL :: LFINN

!  NOAA GBBEP
      LOGICAL :: LGBBEP

!  GFED
      LOGICAL :: LGFED


!  END OF MODULE
      END MODULE LOGICAL_MOD
