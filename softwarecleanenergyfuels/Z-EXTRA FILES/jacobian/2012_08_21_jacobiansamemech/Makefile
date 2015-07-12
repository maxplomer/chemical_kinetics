# Makefile for Linux 64 by Ruiqin Shan

FC = ifort
CC = gcc
LN = gcc


MATROOT = /opt/matlab2009b
MATBITS = glnxa64
DLLEXT = mexa64


# matlab extern include files
MATINC=$(MATROOT)/extern/include
INCLUDES=-I$(MATINC)

# matlab lib files
MATLIB = $(MATROOT)/extern/lib/$(MATBITS)
LPATH =  $(MATROOT)/bin/$(MATBITS)
LIBS = -lmx -lmex -lmat -lm

#intel fortran libs
ILPATH = /opt/intel/Compiler/11.0/069/lib/intel64
ILIBS = -lifcore -limf -lirc -lsvml


# Flags for compilers
FFLAGS = -c -fPIC -vec-report0 -nologo
CFLAGS = -c -fPIC -DMATLAB_MEX_FILE $(INCLUDES)
LFLAGS = -shared -L$(LPATH) -L$(ILPATH) $(LIBS) $(ILIBS) -o $@.$(DLLEXT)
LFLAGS += -Wl,--version-script,$(MATLIB)/mexFunction.map
LFLAGS += -Wl,--no-undefined
LFLAGS += -Wl,-rpath,$(ILPATH)


.c.o:
	$(CC) $(CFLAGS) $<

.F.o:
	$(FC) $(FFLAGS) $<


ALL: jact

JACT_OBJS = jact.o DWDCT.o

jact: $(JACT_OBJS)

	$(LN) $(LFLAGS) $(JACT_OBJS)
clean:
	rm -f *.o *~ 
