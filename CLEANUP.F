!  $ID: CLEANUP.F V01 06/26/2012 16:22 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE CLEANUP DEALLOCATES THE MEMORY ASSIGNED TO DYNAMIC ALLOCATABLE 
!  ARRAYS.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RTM MODEL BY BRUCE. (06/26/2012)
!******************************************************************************
!

      SUBROUTINE CLEANUP

      ! REFERENCES TO F90 MODULES
      USE NEI_MOD,       ONLY : CLEANUP_NEI
      USE INTEX_MOD,     ONLY : CLEANUP_INTEX
      USE FLAMBE_MOD,    ONLY : CLEANUP_FLAMBE
      USE FINN_MOD,      ONLY : CLEANUP_FINN
      USE GBBEP_MOD,     ONLY : CLEANUP_GBBEP
      USE GFED_MOD,      ONLY : CLEANUP_GFED

      IMPLICIT NONE

      !=================================================================
      ! CLEANUP BEGINS HERE!
      !=================================================================

      CALL CLEANUP_NEI
      CALL CLEANUP_INTEX
      CALL CLEANUP_FLAMBE
      CALL CLEANUP_FINN
      CALL CLEANUP_GBBEP
      CALL CLEANUP_GFED

      !=================================================================
      ! PRINT TO SCREEN
      !=================================================================
      WRITE(6, '(/A)') REPEAT('*', 79)
      WRITE(6, 100   )
100   FORMAT('-CLEANUP: Deallocating all the retrieval arrays now...')

      ! RETURN TO THE CALLING PROGRAM
      END SUBROUTINE CLEANUP
