##
## blockSQP -- Sequential quadratic programming for problems with
##             block-diagonal Hessian matrix.
## Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
##
## Licensed under the zlib license. See LICENSE for more details.
##

########################################################################
#                User configuration: qpOASES settings                  #
########################################################################

# Location of qpOASES header files and shared library
QPOASESDIR = ~/numerics/coin-qpOASES
QPOASESINCLUDE = $(QPOASESDIR)/include
QPOASESLIBDIR = $(QPOASESDIR)/bin

########################################################################
#                      End of user configuration                       #
########################################################################


OBJDIR = obj
INCLUDEDIR = include
SRCDIR = src
LIBDIR = lib
DEFS =

OPTIONS = -g -O0 -fPIC -fopenmp -I $(INCLUDEDIR) -I $(QPOASESINCLUDE) -D$(DEFS) -Wno-deprecated -Wno-write-strings -Wall
#OPTIONS = -g -O3 -fPIC -fopenmp -I $(INCLUDEDIR) -I $(QPOASESINCLUDE) -D$(DEFS) -Wno-deprecated -Wno-write-strings -Wall

OBJECTS = $(OBJDIR)/blocksqp_matrix.o \
		$(OBJDIR)/blocksqp_problemspec.o \
		$(OBJDIR)/blocksqp_general_purpose.o \
		$(OBJDIR)/blocksqp_glob.o \
		$(OBJDIR)/blocksqp_hess.o \
		$(OBJDIR)/blocksqp_iter.o \
		$(OBJDIR)/blocksqp_main.o \
		$(OBJDIR)/blocksqp_options.o \
		$(OBJDIR)/blocksqp_qp.o \
		$(OBJDIR)/blocksqp_restoration.o \
		$(OBJDIR)/blocksqp_stats.o

LIBS       = -L $(LIBDIR) -Xlinker -rpath -Xlinker $(LIBDIR) \
             -L $(QPOASESLIBDIR) -Xlinker -rpath -Xlinker $(QPOASESLIBDIR) \
             -L /usr/local/lib -Xlinker -rpath -Xlinker /usr/local/lib \
             -lblockSQP \
             -lqpOASES \
             -llapack

all: blockSQP

library: $(OBJECTS) | $(LIBDIR)
	g++ -shared -fopenmp -o $(LIBDIR)/libblockSQP.so $(OBJECTS)

blockSQP: library $(OBJDIR)/blocksqp_driver.o
	g++ -o blockSQP $(OBJDIR)/blocksqp_driver.o $(LIBS)

min: $(OBJECTS) | $(LIBDIR)
	g++ -shared -o $(LIBDIR)/libblockSQP_min.so $(OBJECTS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp makefile | $(OBJDIR)
	g++ -c $(OPTIONS) -o $@ $<

$(LIBDIR):
	mkdir $(LIBDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(LIBDIR) blockSQP
