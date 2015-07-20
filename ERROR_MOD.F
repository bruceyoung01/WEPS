! $ID: ERROR_MOD.F V01 06/18/2012 10:09 BRUCE EXP$
!
!******************************************************************************
!  MODULE ERROR_MOD CONTAINS ERROR CHECKING ROUTINES.
!
!  ROUTINE LIST:
!  ============================================================================
!  (1 ) ERROR_STOP : PRINT ERROR MESSAGE AND STOP
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE REFERING TO GEOS-CHEM FORMAT. (06/18/2012)
!******************************************************************************
!

      ! INTERFACE:
      MODULE ERROR_MOD

      ! USES:
      IMPLICIT NONE

      PRIVATE

      ! PUBLIC MEMBER FUNCTIONS:
      PUBLIC  :: ALLOC_ERR
      PUBLIC  :: ERROR_STOP

      ! THIS MODULE ROUTINES STARTS HERE ...

      CONTAINS

!------------------------------------------------------------------------------
!
!  $ID: ALLOC_ERR V01 06/30/2012 11:38 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE ALLOC_ERR PRINTS AN ERROR MESSAGE IF THERE IS NOT ENOUGH
!  MEMORY TO ALLOCATE A PARTICULAR ALLOCATABLE ARRAY
!
!  ARGUMENTS AS INPUT:
!  ============================================================================
!  (1 ) ARRAYNAME    (CHARACTER) : ALLOCATING ARRAY NAME
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RTM MODEL. (06/30/2012)
!******************************************************************************

      SUBROUTINE ALLOC_ERR(ARRAYNAME, AS)

      ! ARGUMENTS
      CHARACTER (LEN = *), INTENT(IN) :: ARRAYNAME
      INTEGER, OPTIONAL, INTENT(IN)   :: AS

      ! LOCAL VARIABLES
      CHARACTER (LEN = 255)           :: IFORT_ERRMSG, MSG
      CHARACTER (LEN = 255)           :: ERRMSG

      ! DEFINE ERROR MESSAGE
      ERRMSG = 'Allocation error in array: ' //TRIM(ARRAYNAME)

      IF (PRESENT(AS)) THEN
       ERRMSG = TRIM(ERRMSG)//' :: '//TRIM(MSG)
      ENDIF

      ! PRINT ERROR MESSAGE, DEALLOCATE MEMORY, AND STOP THE RUN
      CALL ERROR_STOP (ERRMSG, 'ALLOC_ERR.F')

      END SUBROUTINE ALLOC_ERR


!------------------------------------------------------------------------------

!  $ID: ERROR_STOP V01 06/26/2012 16:10 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE ERROR_STOP PRINTS AN ERROR MESSAGE AND STOP.
!
!  ARGUMENTS AS INPUT:
!  ============================================================================
!  (1 ) MESSAGE    (CHARACTER) : ERROR MESSAGE TO BE PRINTED
!  (2 ) LOCATION   (CHARACTER) : LOCATION WHERE ERROR_STOP IS CALLED FROM
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM GEOS-CHEM MODEL BY BRUCE. (06/26/2012)
!******************************************************************************
!
      SUBROUTINE ERROR_STOP(MESSAGE, LOCATION)

      ! ARGUMENTS
      CHARACTER (LEN = *), INTENT(IN) :: MESSAGE, LOCATION

      !=================================================================
      ! ERROR_STOP BEGINS HERE!
      !=================================================================

      ! WRITE MSG
      WRITE(6, '(A)') REPEAT('=', 79)
      WRITE(6, '(A)') 'PROGRAMING ERROR : '//TRIM(MESSAGE)
      WRITE(6, '(A)') 'STOP AT          : '//TRIM(LOCATION)
      WRITE(6, '(A)') REPEAT('=', 79)

      ! DEALLOCATE MEMORY AND STOP THE RUN
      WRITE(6, '(/A)') 'Deallocate memory and stop the run.'
!     CALL CLEANUP
      CALL EXIT(99999)

      ! RETURN TO CALLING PROGRAM
      END SUBROUTINE ERROR_STOP

!------------------------------------------------------------------------------

      END MODULE ERROR_MOD
