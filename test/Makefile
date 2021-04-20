vpath %.c ..
INCL = -I..
.c.o:
	$(CC) $(CFLAGS) $(INCL) -c -o ${@F} $<
#CFLAGS = -g -Wall -Wshadow -Wstrict-prototypes -DNDEBUG -O3
CFLAGS = -g -Wall -Wshadow -Wstrict-prototypes
CC = gcc
OBJ = likelihood.o chol.o amoeba.o misc.o io.o getcic.o zeroin.o symeig.o

XNINFOMAT= xninfomat.o ninfomat.o misc.o io.o getcic.o zeroin.o chol.o
xninfomat : $(XNINFOMAT)
	$(CC) $(CFLAGS) -o $@ $(XNINFOMAT) -lm

xagentdef: xagentdef.o io.o getcic.o zeroin.o misc.o
	$(CC) $(CFLAGS) -o $@ xagentdef.o io.o getcic.o zeroin.o misc.o -lm

xbonedef: xbonedef.o io.o getcic.o zeroin.o misc.o 
	$(CC) $(CFLAGS) -o $@ xbonedef.o io.o getcic.o zeroin.o misc.o -lm

xread_counts: xread_counts.o io.o getcic.o zeroin.o misc.o 
	$(CC) $(CFLAGS) -o $@ xread_counts.o io.o getcic.o zeroin.o misc.o -lm

XSYMEIG=xsymeig.o symeig.o misc.o
xsymeig: $(XSYMEIG)
	$(CC) $(CFLAGS) -o $@ $(XSYMEIG) -lm

XSYMEIG2=xsymeig2.o symeig.o misc.o io.o getcic.o zeroin.o
xsymeig2: $(XSYMEIG2)
	$(CC) $(CFLAGS) -o $@ $(XSYMEIG2) -lm

XLNL=xlnl.o likelihood.o misc.o io.o  zeroin.o getcic.o chol.o symeig.o 
xlnl : $(XLNL)
	$(CC) $(CFLAGS) -o $@ $(XLNL) -lm

XLUMP=xlump.o misc.o io.o
xlump : $(XLUMP)
	$(CC) $(CFLAGS) -o $@ xlump.o $(XLUMP) -lm

XXLUMP=xxlump.o misc.o io.o header.o likelihood.o unirand.o getcic.o \
       zeroin.o chol.o symeig.o
xxlump : $(XXLUMP)
	$(CC) $(CFLAGS) -o $@ $(XXLUMP) -lm

xxTax : xxTax.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ xxTax.o $(OBJ) -lm

clean : 
	rm -f *.o junk*
