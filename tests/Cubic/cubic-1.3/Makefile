###
### CONFIGURATION
###

# 1) Change this IF you used Configure --prefix <MY_PREFIX_DIRECTORY>
# PARIHOME=MY_PREFIX_DIRECTORY
# PARIHOME=$(GP_INSTALL_PREFIX)
PARIHOME=/home/cs-gb2khv2/Downloads/pari-2.11.2/GPDIR

# 2) Change this IF your pari library was not compiled with gmp (why ?), or
# if it has a non-standard name, e.g.
#LIB=-lpari-2-9-0 -lm
LIB=-lgmp -lpari -lm

###
### NO USER SERVICEABLE PARTS BELOW
###
VERSION=1.3
PARIDIR=$(PARIHOME)/include/pari
PARILIB=$(PARIHOME)/lib

INCLUDEDIR= -I$(PARIDIR)

DBGFLAG=-g -Wall -Wimplicit -Wmissing-prototypes -L$(PARILIB)
OPTFLAG=-O3 -Wall -DGCC_INLINE
PRFFLAG=-g -pg $(OPTFLAG)

CC=gcc
LD=gcc
LDFLAGS=-L$(PARILIB) -Wl,-rpath,$(PARILIB)
PROFLIB=/home/kb/pari/Olinux-i686.prf/libpari.a

OBJ=util.o ccubic.o rcubic.o
OBJNOPRINT=util-noprint.o ccubic.o rcubic.o
OBJDBG=util.dbg.o ccubic.dbg.o rcubic.dbg.o
OBJPRF=util.prf.o ccubic.prf.o rcubic.prf.o

all: cubic
dbg: cubic.dbg
prf: cubic.prf
noprint: cubic-noprint

tags::
	ctags *.[ch]

%.o: %.c cubic.h cubicinl.h
	$(CC) -DPRINT $(OPTFLAG) $(INCLUDEDIR) -o $@ -c $<
%-noprint.o: %.c cubic.h cubicinl.h
	$(CC) $(OPTFLAG) $(INCLUDEDIR) -o $@ -c $<
%.dbg.o: %.c cubic.h cubicinl.h
	$(CC) -DPRINT $(DBGFLAG) $(INCLUDEDIR) -o $@ -c $<
%.prf.o: %.c cubic.h cubicinl.h
	$(CC) $(PRFFLAG) $(INCLUDEDIR) -o $@ -c $<

all: cubic cubic-noprint

cubic: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIB)
cubic.dbg: $(OBJDBG)
	$(CC) $(LDFLAGS) $(DBGFLAG) -o $@ $(OBJDBG) $(LIB)
cubic.prf: $(OBJPRF)
	$(CC) $(PRFFLAG) -o $@ $(OBJPRF) $(PROFLIB) $(LIB)
cubic-noprint: $(OBJNOPRINT)
	$(LD) $(LDFLAGS) -o $@ $(OBJNOPRINT) $(LIB)

clean::
	-rm -f *.o cubic cubic.dbg cubic.prf cubic-noprint tags

distrib:: clean
	 cwd=`pwd`; cd /tmp; ln -s $$cwd cubic-$(VERSION); tar czf cubic-$(VERSION).tgz cubic-$(VERSION)/* --exclude=.git
	 mv /tmp/cubic-$(VERSION).tgz .

tags::
	ctags *.c *.h
