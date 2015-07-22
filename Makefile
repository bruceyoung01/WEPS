#  $ID: Makefile V01 06/28/2012 11:09 BRUCE EXP$
#
#******************************************************************************
#  Makefile IS USED TO COMPILE WEPS.
#
#  NOTES:
#  ============================================================================
#  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/28/2012)
#  =====================bruces-iMac COMPATIBLE MODIFICATION START==============
#  (2 ) REMOVE READING ENVIRONMENT VARIABLES FROM configure.wps BY BRUCE. 
#       (11/18/2013)
#  (3 ) MODIFY FORTRAN COMPILER FROM IFORT TO GFORTRAN (gfortran-mp-4.8), WHICH 
#       IS INSTALLED FROM MACPORTS. (11/18/2013)
#  (4 ) IN ORDER TO BE CONSISTENT WITH GFORTRAN (gfortran-mp-4.8), NETCDF IS 
#       USED THE ONE COMPILED FROM GFORTRAN (gfortran-mp-4.8) (/opt/local), 
#       INSTEAD OF /usr/local. BRUCE (11/18/2013)
#  (5 ) TO BE CONSISTENT WITH F77 AND F90, THE FOLLOWING MODULES HAVE BEEN 
#       MODIFIED FROM .F TO .F90:
#       WEPS.F90,       TOOL_KIT_MOD.F90, NEI_MOD.F90,   INTEX_MOD.F90,
#       FLAMBE_MOD.F90, FINN_MOD.F90,     GBBEP_MOD.F90, GFED_MOD.F90, 
#       SEVIRI_MOD.F90, GFAS_MOD.F90,     QFED_MOD.F90,  WRITE_MOD.F90
#  =====================bruces-iMac COMPATIBLE MODIFICATION END================
#******************************************************************************

#  NAME OF COMPILER
NETCDF      = /opt/local
FFLAGS      = -c -g -I${NETCDF}/include -L${NETCDF}/lib -lnetcdf -lnetcdff -fno-range-check
FC          = gfortran-mp-4.8 $(FFLAGS)
FCC         = gfortran-mp-4.8 -free $(FFLAGS)
LINK.f      = gfortran-mp-4.8 -I${NETCDF}/include -L${NETCDF}/lib -lnetcdf -lnetcdff -fno-range-check

#  USED LIBRARIES
LIBS        =

#  NAME OF EXECUTABLE FILE
APPLICATION =  'WEPS.EXE'

#  OBJECT MODULES
OBJ_MAIN    = WEPS.o
OBJ_MODULES = PARALLEL_MODULE.o         \
              TOOL_KIT_MOD.o            \
              PARAMETER_MOD.o           \
              CONSTANTS_MOD.o           \
              GEOMETRY_MOD.o            \
              MODULE_DEBUG.o            \
              ERROR_MOD.o               \
              CHARPAK_MOD.o             \
              MISC_DEFINITIONS_MODULE.o \
              MODULE_MAP_UTILS.o        \
              NAMELIST_ARRAY_MOD.o      \
              NAMELIST_MOD.o            \
              GRID_MOD.o                \
              NEI_MOD.o                 \
              INTEX_MOD.o               \
              GFED_MOD.o                \
              FLAMBE_MOD.o              \
              FINN_MOD.o                \
              GBBEP_MOD.o               \
              SEVIRI_MOD.o              \
              GFAS_MOD.o                \
              QFED_MOD.o                \
              WRITE_MOD.o
#             CLEANUP.o
              
#  EXECUTABLES
#all: $(APPLICATION)

$(APPLICATION):   $(OBJ_MODULES)  \
                  $(OBJ_MAIN)
	$(LINK.f) $(OBJ_MODULES)  \
	          $(OBJ_MAIN)     \
                  $(LIBS) -o $(APPLICATION)

#  COMPILATIONS
WEPS.o				: WEPS.F90
PARALLEL_MODULE.o		: PARALLEL_MODULE.F90
TOOL_KIT_MOD.o			: TOOL_KIT_MOD.F90
PARAMETER_MOD.o			: PARAMETER_MOD.F
CONSTANTS_MOD.o			: CONSTANTS_MOD.F
GEOMETRY_MOD.o			: GEOMETRY_MOD.F
MODULE_DEBUG.o			: MODULE_DEBUG.F90
ERROR_MOD.o			: ERROR_MOD.F
CHARPAK_MOD.o			: CHARPAK_MOD.F
MISC_DEFINITIONS_MODULE.o	: MISC_DEFINITIONS_MODULE.F90
MODULE_MAP_UTILS.o		: MODULE_MAP_UTILS.F90
NAMELIST_ARRAY_MOD.o		: NAMELIST_ARRAY_MOD.F
NAMELIST_MOD.o			: NAMELIST_MOD.F
GRID_MOD.o                      : GRID_MOD.F90
NEI_MOD.o			: NEI_MOD.F90
INTEX_MOD.o			: INTEX_MOD.F90
GFED_MOD.o			: GFED_MOD.F90
FLAMBE_MOD.o			: FLAMBE_MOD.F90
FINN_MOD.o			: FINN_MOD.F90
GBBEP_MOD.o			: GBBEP_MOD.F90
SEVIRI_MOD.o			: SEVIRI_MOD.F90
GFAS_MOD.o			: GFAS_MOD.F90
QFED_MOD.o			: QFED_MOD.F90
WRITE_MOD.o			: WRITE_MOD.F90

clean:
	rm -f ./*.o *.mod *.EXE


.SUFFIXES: .f .F .f90 .F90
.f.o:			; $(FC) -c $*.f
.F.o:			; $(FC) -c $*.F
.f90.o:			; $(FCC) -c $*.f90
.F90.o:			; $(FCC) -c $*.F90
