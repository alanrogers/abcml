vpath %.c ../src
INCL := -I../src
.c.o:
	$(CC) $(CFLAGS) -c -o ${@F} $<
#CFLAGS = -g -Wall -Wshadow -Wstrict-prototypes -DNDEBUG -O3
CFLAGS = -g -Wall -Wshadow -Wstrict-prototypes $(INCL)
CC = gcc

XNINFOMAT := xninfomat.o ninfomat.o misc.o io.o getcic.o zeroin.o chol.o
xninfomat : $(XNINFOMAT)
	$(CC) $(CFLAGS) -o $@ $(XNINFOMAT) -lm

XAGENTDEF := xagentdef.o io.o getcic.o zeroin.o misc.o
xagentdef: $(XAGENTDEF)
	$(CC) $(CFLAGS) -o $@ $(XAGENTDEF) -lm

XBONEDEF := xbonedef.o io.o getcic.o zeroin.o misc.o 
xbonedef : $(XBONEDEF)
	$(CC) $(CFLAGS) -o $@ $(XBONEDEF) -lm

XREAD_COUNTS := xread_counts.o io.o getcic.o zeroin.o misc.o
xread_counts: $(XREAD_COUNTS)
	$(CC) $(CFLAGS) -o $@ $(XREAD_COUNTS) -lm

XSYMEIG := xsymeig.o symeig.o misc.o
xsymeig: $(XSYMEIG)
	$(CC) $(CFLAGS) -o $@ $(XSYMEIG) -lm

XSYMEIG2 := xsymeig2.o symeig.o misc.o io.o getcic.o zeroin.o
xsymeig2: $(XSYMEIG2)
	$(CC) $(CFLAGS) -o $@ $(XSYMEIG2) -lm

XLNL := xlnL.o likelihood.o misc.o io.o  zeroin.o getcic.o chol.o symeig.o 
xlnL : $(XLNL)
	$(CC) $(CFLAGS) -o $@ $(XLNL) -lm

# XLUMP doesn't work.
#XLUMP := xlump.o misc.o io.o
#xlump : $(XLUMP)
#	$(CC) $(CFLAGS) -o $@ xlump.o $(XLUMP) -lm

XXLUMP=xxlump.o misc.o io.o likelihood.o unirand.o getcic.o \
       zeroin.o chol.o symeig.o
xxlump : $(XXLUMP)
	$(CC) $(CFLAGS) -o $@ $(XXLUMP) -lm

# Doesn't work. Commented out for now. Details in xxTax.c
#XXTAX := xxTax.o likelihood.o chol.o amoeba.o misc.o io.o getcic.o zeroin.o \
#         symeig.o
#xxTax : $(XXTAX)
#	$(CC) $(CFLAGS) -o $@ $(XXTAX) -lm

clean : 
	rm -f *.o junk*
