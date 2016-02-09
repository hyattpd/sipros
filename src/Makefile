CC=	g++
SRC=	alphabet.cpp \
	bitap.cpp \
	blosum.cpp \
	config.cpp \
	dbpeptide.cpp \
	denovo.cpp \
	directoryStructure.cpp \
	isotopologue.cpp \
	main.cpp \
	MassSpectralPeakCluster.cpp \
	MassSpectralPeak.cpp \
	ms2scan.cpp \
	outputtable.cpp \
	pattern.cpp \
	PeptideScorer.cpp \
	profrag.cpp \
	proNovoConfig.cpp \
	ptm.cpp \
	scanindex.cpp \
	scanmass.cpp \
	TandemMassSpectrum.cpp \
	text.cpp \
	tinystr.cpp \
	tinyxml.cpp \
	tinyxmlerror.cpp \
	tinyxmlparser.cpp \
	tokenvector.cpp 
PROC=	sipros
CFLAGS=	-Wall -O2
LIBS=	-lm
LDFLAGS=	$(LIBS)
OBJS=	${SRC:.cpp=.o}
all:	$(PROC)
clean:	$(OBJS)
	/bin/rm -rf $(OBJS) $(PROC)
$(PROC):	$(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS) $(CFLAGS)
alphabet.o:	alphabet.cpp
	$(CC) -c $(CFLAGS) alphabet.cpp
bitap.o:	bitap.cpp
	$(CC) -c $(CFLAGS) bitap.cpp
blosum.o:	blosum.cpp
	$(CC) -c $(CFLAGS) blosum.cpp
config.o:	config.cpp
	$(CC) -c $(CFLAGS) config.cpp
dbpeptide.o:	dbpeptide.cpp
	$(CC) -c $(CFLAGS) dbpeptide.cpp
denovo.o:	denovo.cpp
	$(CC) -c $(CFLAGS) denovo.cpp
directoryStructure.o:	directoryStructure.cpp
	$(CC) -c $(CFLAGS) directoryStructure.cpp
isotopologue.o:	isotopologue.cpp
	$(CC) -c $(CFLAGS) isotopologue.cpp
main.o:	main.cpp
	$(CC) -c $(CFLAGS) main.cpp
MassSpectralPeakCluster.o:	MassSpectralPeakCluster.cpp
	$(CC) -c $(CFLAGS) MassSpectralPeakCluster.cpp
MassSpectralPeak.o:	MassSpectralPeak.cpp
	$(CC) -c $(CFLAGS) MassSpectralPeak.cpp
ms2scan.o:	ms2scan.cpp
	$(CC) -c $(CFLAGS) ms2scan.cpp
outputtable.o:	outputtable.cpp
	$(CC) -c $(CFLAGS) outputtable.cpp
pattern.o:	pattern.cpp
	$(CC) -c $(CFLAGS) pattern.cpp
PeptideScorer.o:	PeptideScorer.cpp
	$(CC) -c $(CFLAGS) PeptideScorer.cpp
profrag.o:	profrag.cpp
	$(CC) -c $(CFLAGS) profrag.cpp
proNovoConfig.o:	proNovoConfig.cpp
	$(CC) -c $(CFLAGS) proNovoConfig.cpp
ptm.o:	ptm.cpp
	$(CC) -c $(CFLAGS) ptm.cpp
scanindex.o:	scanindex.cpp
	$(CC) -c $(CFLAGS) scanindex.cpp
scanmass.o:	scanmass.cpp
	$(CC) -c $(CFLAGS) scanmass.cpp
TandemMassSpectrum.o:	TandemMassSpectrum.cpp
	$(CC) -c $(CFLAGS) TandemMassSpectrum.cpp
text.o:	text.cpp
	$(CC) -c $(CFLAGS) text.cpp
tinystr.o:	tinystr.cpp
	$(CC) -c $(CFLAGS) tinystr.cpp
tinyxml.o:	tinyxml.cpp
	$(CC) -c $(CFLAGS) tinyxml.cpp
tinyxmlerror.o:	tinyxmlerror.cpp
	$(CC) -c $(CFLAGS) tinyxmlerror.cpp
tinyxmlparser.o:	tinyxmlparser.cpp
	$(CC) -c $(CFLAGS) tinyxmlparser.cpp
tokenvector.o:	tokenvector.cpp
	$(CC) -c $(CFLAGS) tokenvector.cpp
