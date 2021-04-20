#  where the executables should be installed:
BIN = $(HOME)/bin

# create string defining git version.
git_version := "$(shell git describe --dirty --always --tags)"

# define your compiler here:
CC = gcc

# define your copy command here:
COPY = cp

# define your remove (or delete) command here:
RM = rm -f

# suffix of object files?  (o for unix, obj for MSDOS):
O = o

# suffix for executable files? (empty for unix, .exe for dos):
E =

# command line arguments for C compiler
WARN1 = -Wall -Wshadow -Wpointer-arith \
  -Wcast-qual -Wcast-align -Wwrite-strings  \
  -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes \
  -Wmissing-declarations -Winline -pedantic 
WARN2 = -Wunused
OPT1 = -DNDEBUG
#OPT2 = -O
CFLAGS = -g $(WARN1) $(WARN2) $(OPT1) $(OPT2) \
 -DGIT_VERSION=\"$(git_version)\"

#
# There should be no need to make changes below this line
#
TARGETS = abcml$(E) abcsim$(E) tabcfg$(E) mau2cfg$(E) cplcfg$(E) gnaw$(E) \
  wgtcfg$(E)

#make all targets
all : $(TARGETS)
	touch all

#target wgtcfg: Weight the configurations in a .cfg file
WGTCFG = wgtcfg.$(O) getcic.$(O) io.$(O) header.$(O) misc.$(O) zeroin.$(O)
wgtcfg$(E) : $(WGTCFG)
	$(CC) $(CFLAGS) -o $@ $(WGTCFG) -lm

#target gnaw: subject an assemblage to attrition
GNAW = gnaw.$(O)  getcic.$(O) io.$(O) header.$(O) misc.$(O) zeroin.$(O) \
  unirand.$(O)
gnaw$(E) : $(GNAW)
	$(CC) $(CFLAGS) -o $@ $(GNAW) -lm

#target cplcfg: Calculate complement of .cfg file
CPLCFG = cplcfg.$(O) getcic.$(O) io.$(O) header.$(O) misc.$(O) zeroin.$(O)
cplcfg$(E) : $(CPLCFG)
	$(CC) $(CFLAGS) -o $@ $(CPLCFG) -lm

#target mau2cfg: convert .mau format to .cfg format
MAU2CFG = mau2cfg.$(O) getcic.$(O) misc.$(O) io.$(O) header.$(O) zeroin.$(O)
mau2cfg$(E) : $(MAU2CFG)
	$(CC) $(CFLAGS) -o $@ $(MAU2CFG) -lm

#target abcsim: simulation program
ABCSIM = abcsim.$(O) getcic.$(O) zeroin.$(O) misc.$(O) io.$(O) unirand.$(O) \
  header.$(O)
abcsim$(E) : $(ABCSIM) 
	$(CC) $(CFLAGS) -o $@ $(ABCSIM) -lm

#target abcml: estimation program
ABCML = abcml.$(O) io.$(O) getcic.$(O) zeroin.$(O) misc.$(O) likelihood.$(O) \
  prlike.$(O) chol.$(O) amoeba.$(O)  symeig.$(O) unirand.$(O) header.$(O) \
  ninfomat.$(O)
abcml$(E) : $(ABCML) 
	$(CC) $(CFLAGS) -o $@ $(ABCML) -lm

#target tabcfg: tabulate configurations
TABCFG = tabcfg.$(O) getcic.$(O) misc.$(O) io.$(O) header.$(O) zeroin.$(O)
tabcfg$(E) : $(TABCFG) 
	$(CC) $(CFLAGS) -o $@ $(TABCFG) -lm

install : all
	$(COPY) $(TARGETS) $(BIN)

clean : 
	$(RM) *.$(O) *~ all bone junk* core abcmlsrc.zip try dep $(TARGETS)

abcmlsrc.zip :
	zip -l abcmlsrc LICENSE Makefile.dos Makefile.unix *.c *.h -x try.*

depend :
	echo '#Automatically generated dependency info' > depend
	$(CC) -M *.c >> depend

include depend
