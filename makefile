# Makefile for standalone version of blockSQP

QPOASESINCLUDE = $(QPOASESDIR)/include_aw
QPOASESLIBDIR = $(QPOASESDIR)/bin

OBJDIR = obj
INCLUDEDIR = include
SRCDIR = src
LIBDIR = lib
DEFS =

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
             -lqpOASES_aw \
             -llapack

OPTIONS = -g -O0 -fPIC -I $(INCLUDEDIR) -I $(QPOASESINCLUDE) $(DEFS) -Wno-deprecated -Wno-write-strings

all: blockSQP

library: $(OBJECTS) | $(LIBDIR)
	g++ -shared -o $(LIBDIR)/libblockSQP.so $(OBJECTS)

blockSQP: library $(OBJDIR)/blocksqp_driver.o
	g++ -o blockSQP $(OBJDIR)/blocksqp_driver.o $(LIBS)

min: $(OBJECTS) | $(LIBDIR)
	g++ -shared -o $(LIBDIR)/libblockSQP_min.so $(OBJECTS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	g++ -c $(OPTIONS) -o $@ $<

$(LIBDIR):
	mkdir $(LIBDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(LIBDIR) blockSQP
