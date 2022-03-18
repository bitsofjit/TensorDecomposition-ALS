# makefile
# Purpose: Create the executable for ALS
#
include makefile.inc

#MACPORTS = /opt/local

MACPORTS = /usr/local
LIBPATH = -L$(MACPORTS)/lib/

EXT_LIBS  = $(LIBPATH) -llapack -lblas

# On Linux can use: EXT_LIBS  = -lblas -llapack /usr/lib/libg2c.so.0
# (if DEFS is undefined in makefile.inc, then blas/lapack are not needed)

EXEC = ../bin/sepals

EXECTST = ../bin/quad

SRC   =  nrtype.f90 func_decom_als.f90

SRCTST   =  nrtype.f90 quad_test.f90

OBJS = ${SRC:.f90=.o}

OBJSTST = ${SRCTST:.f90=.o}

MOD  = ${SRC:.f90=.mod}

sepals: 
	@echo "Compiling $(G95) $(GDBFLAG) $(FFLAGS) $(SRC)"
	@$(G95) $(GDBFLAG) $(FFLAGS) $(SRC)
	@echo "Linking..."
	@$(G95) -o $(EXEC) $(OBJS) $(EXT_LIBS)

inttest:
	@echo "Compiling $(G95) $(GDBFLAG) $(FFLAGS) $(SRCTST)"
	@$(G95) $(GDBFLAG) $(FFLAGS) $(SRCTST)
	@echo "Linking..."
	@$(G95) -o $(EXECTST) $(OBJSTST) $(EXT_LIBS)

clean   :
	@echo "Removing *.o and ../bin/sepals"
	@$(RM) $(OBJS)
	@$(RM) $(EXEC)

