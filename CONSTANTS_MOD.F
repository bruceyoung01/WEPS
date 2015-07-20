!  $ID: CONSTANTS_MOD.F V01 06/28/2012 08:47 BRUCE EXP$
!
!******************************************************************************
!  MODULE CONSTANTS_MOD SETS UP PARAMETERS OF PHYSICS, MATHEMATICS, CHEMISTRY, 
!  AND SO FORTH. THEY ARE USED IN OTHER MODULES
!
!  VARIABLES:
!  ============================================================================
!  MATHEMATICS:
!  (1 ) PI       (REAL*8)  : CIRCUMFERENCE RATIO/ARCHIMEDES' CONSTANT/
!                            LUDOLPH'S NUMBER
!  (2 ) ZEROI    (INTEGER) : INTEGER TYPE ZERO
!  (3 ) ZEROR    (REAL*8)  : REAL TYPE ZERO
!  (4 ) ONEI     (INTEGER) : INTEGER TYPE ONE
!  (5 ) ONER     (REAL*8)  : REAL TYPE ONE
!  (6 ) IMAGU    (COMPLEX) : IMAGINARY UNIT
!  (7 ) E        (REAL*8)  : NAPIER'S CONSTANT/EULER'S NUMBER/BASE OF 
!                            NATURAL LOGARITHM
!  (8 ) GAMMAE   (REAL*8)  : EULER-MASCHERONI CONSTANT
!  (9 ) PHI      (REAL*8)  : GOLDEN RATIO
!  (10) RHO      (REAL*8)  : PLASTIC CONSTANT
!  (11) BETAS    (REAL*8)  : EMBREE-TREFETHEN CONSTANT
!  (12) DELTA    (REAL*8)  : FEIGENBAUM CONSTANT
!  (13) ALPHA    (REAL*8)  : FEIGENBAUM CONSTANT
!  
!  PHYSICS:
!  (1 ) EC       (REAL*8)  : ELEMENTARY CHARGE             [C]
!  (2 ) G        (REAL*8)  : GRAVITATIONAL CONSTANT
!  [m^3kg^-1s^-2]
!  (3 ) C        (REAL*8)  : SPEED OF LIGHT IN VACUUM      [m/s]
!  (4 ) H        (REAL*8)  : PLANCK'S CONSTANT             [Js]
!  (5 ) HBAR     (REAL*8)  : DIRAC'S CONSTANT              [Js]
!  (6 ) MUB      (REAL*8)  : BOHR MAGNETON                 [Am^2]
!  (7 ) ALPHA0   (REAL*8)  : BOHR RADIUS                   [A]
!  (8 ) SIGMA    (REAL*8)  : STEFAN-BOLTZMANN'S CONSTANT   [Wm^2K^-4]
!  (9 ) KW       (REAL*8)  : WIEN'S CONSTANT               [mK]
!  (10) R        (REAL*8)  : MOLAR GASCONSTANT             [J/mol]
!  (11) NA       (REAL*8)  : AVOGADRO'S CONSTANT           [mol^-1]
!  (12) K        (REAL*8)  : BOLTZMANN'S CONSTANT          [J/K]
!  (13) D0       (REAL*8)  : DIAMETER OF THE SUN           [m]
!  (14) M0       (REAL*8)  : MASS OF THE SUN               [kg]
!  (15) T0       (REAL*8)  : ROTATIONAL PERIOD OF THE SUN  [days]
!  (16) RA       (REAL*8)  : RADIUS OF EARTH               [m]
!  (17) MA       (REAL*8)  : MASS OF EARTH                 [kg]
!  (18) TA       (REAL*8)  : ROTATIONAL PERIOD OF EARTH    [hours]
!  (19) LJ       (REAL*8)  : LIGHT YEAR                    [m]
!  (20) RTOD     (REAL*8)  : RADIAN TO DEGREE              [DEGREE/RAD]
!  (21) DTOR     (REAL*8)  : DEGREE TO RADIAN              [RAD/DEGREE]
!  (22) TAUP     (REAL*8)  : DRY LAPSE RATE                [K/km]
!
!  ATMOSPHERIC SCIENCES:
!  (1 ) EARTH_RADIUS_M  (REAL*8)  : EARTH MEAN RADIUS      [m]
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/28/2012)
!******************************************************************************


      MODULE CONSTANTS_MOD

      ! MAKE EVERYTHING PUBLIC
      PUBLIC

!-----------------------------MATHEMATICAL CONSTANTS---------------------------
      INTEGER       :: ZEROI  = 0
      REAL*8        :: ZEROR  = 0.0
      INTEGER       :: ONEI   = 1
      REAL*8        :: ONER   = 1.0
      REAL*8        :: E      = 2.718281828459
      REAL*8        :: GAMMAE = 0.57721566490153286060651209008240243
      REAL*8        :: PHI    = 1.61803398874989484820458683436563811
      REAL*8        :: RHO    = 1.32471795724474602596090885447809734
      REAL*8        :: BETAS  = 0.70258
      REAL*8        :: DELTA  = 4.66920160910299067185320382046620161
      REAL*8        :: ALPHA  = 2.50290787509589282228390287321821578
!-----------------------------PHYSICAL CONSTANTS-------------------------------
      REAL*8        :: EC     = 1.60217733*10**(-19)
      REAL*8        :: G      = 6.67259*10**(-11)
      REAL*8        :: C      = 2.99792458*10**(8)
      REAL*8        :: H      = 6.6260755*10**(-34)
      REAL*8        :: HBAR   = 1.0545727*10**(-34)
      REAL*8        :: MUB    = 9.2741*10**(-24)
      REAL*8        :: ALPHA0 = 0.52918
      REAL*8        :: SIGMA  = 5.67032*10**(-8)
      REAL*8        :: KW     = 2.8978*10**(-3)
      REAL*8        :: R      = 8.31441
      REAL*8        :: NA     = 6.0221367*10**(23)
      REAL*8        :: K      = 1.380658*10**(-23)
      REAL*8        :: D0     = 1.392*10**(9)
      REAL*8        :: M0     = 1.989*10**(30)
      REAL*8        :: T0     = 25.38
      REAL*8        :: RA     = 6.378*10**(6)
      REAL*8        :: MA     = 5.976*10**(24)
      REAL*8        :: TA     = 23.96
      REAL*8        :: LJ     = 9.4605*10**(15)
      REAL*8        :: RTOD   = 57.29577951
      REAL*8        :: DTOR   = 0.01745329252
      REAL*8        :: TAUP   = 9.8

!-----------------------------ATMOSPHERIC SCIENCES-----------------------------
      REAL, PARAMETER :: PI          = 3.141592653589793
      REAL, PARAMETER :: DEG_PER_RAD = 180./PI
      REAL, PARAMETER :: RAD_PER_DEG = PI/180.
 
      ! FOR WPS MAP
      ! ANGULAR ROTATION RATE OF THE EARTH
      REAL, PARAMETER :: OMEGA_E = 7.292e-5

      ! MEAN EARTH RADIUS IN M.  THE VALUE BELOW IS CONSISTENT
      ! WITH NCEP'S ROUTINES AND GRIDS.
      REAL, PARAMETER :: A_WGS84  = 6378137.
      REAL, PARAMETER :: B_WGS84  = 6356752.314
      REAL, PARAMETER :: RE_WGS84 = A_WGS84
      REAL, PARAMETER :: E_WGS84  = 0.081819192
      REAL, PARAMETER :: A_NAD83  = 6378137.
      REAL, PARAMETER :: RE_NAD83 = A_NAD83
      REAL, PARAMETER :: E_NAD83  = 0.0818187034
      REAL, PARAMETER :: EARTH_RADIUS_M = 6370000.! SAME AS MM5 SYSTEM
      REAL, PARAMETER :: EARTH_CIRC_M = 2.*PI*EARTH_RADIUS_M

      END MODULE CONSTANTS_MOD
