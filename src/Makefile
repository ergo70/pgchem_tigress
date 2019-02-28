# $PostgreSQL: pgsql/contrib/pgchem/Makefile
# Makefile for linux-x64
# Adjust the locations of OpenBabel according to your system
#
override with_llvm=no

OPENBABEL_SOURCE=openbabel-2.4.1
OPENBABEL_SOURCE_DIR=$(OPENBABEL_SOURCE)/src
OPENBABEL_INCLUDE_DIRS=-I$(OPENBABEL_SOURCE)/build/include/ -I$(OPENBABEL_SOURCE)/include/ -I$(OPENBABEL_SOURCE)/include/openbabel/ -I$(OPENBABEL_SOURCE)/build/ -I../../../src/backend/
OPENBABEL_LIB_DIRS=-L$(OPENBABEL_SOURCE)/build/lib/
MATCHMOL_LIB_DIR=barsoi/

MODULE_big = libpgchem
EXTENSION = pgchem_tigress
DATA = pgchem_tigress--4.0.sql pgchem_tigress--unpackaged--4.0.sql
PGFILEDESC = "pgchem_tigress - the chemoinformatics extension for PostgreSQL"

OBJS = reaction/reaction_op.o reaction/reaction_gist.o reaction/reaction_io.o molecule/molecule_op.o molecule/molecule_gist.o molecule/molecule_io.o obwrapper.o functions.o compression.o openbabel_addon/fingerpc8.o
SHLIB_LINK = -L${MATCHMOL_LIB_DIR} -lbarsoi ${OPENBABEL_LIB_DIRS} -lopenbabel -lz
PG_CPPFLAGS=-m64 -fPIC -O2 -finline-functions -funroll-loops -I. ${OPENBABEL_INCLUDE_DIRS} 

ifdef USE_PGXS
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
else
subdir = contrib/pgchem_tigress
top_builddir = ../../..
include $(top_builddir)/src/Makefile.global
include $(top_srcdir)/contrib/contrib-global.mk
endif
