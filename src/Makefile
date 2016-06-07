CC=gcc
INC=-I/opt/local/include/
L=-L/opt/local/lib/
LIBR=-lm -lsndfile -lconfig -lgsl  -lgslcblas
BIN=../bin/

# Uncomment for speed optimization
#optspeed=-O3

all:	Analysis Synthesis
AnalysisFunctions.o:	AnalysisFunctions.c AnalysisFunctions.h
	$(CC) -c AnalysisFunctions.c $(optspeed) $(L) $(INC)
Analysis.o:	Analysis.c AnalysisFunctions.h
	$(CC) -c Analysis.c $(optspeed) $(L) $(INC)
Analysis: AnalysisFunctions.o Analysis.o
	$(CC) -o Analysis AnalysisFunctions.o Analysis.o $(optspeed) $(L) $(INC) $(LIBR) 

SynthesisFunctions.o:	SynthesisFunctions.c SynthesisFunctions.h
	$(CC) -c SynthesisFunctions.c $(optspeed) $(L) $(INC)
Synthesis.o:	Synthesis.c SynthesisFunctions.h
	$(CC) -c Synthesis.c $(optspeed) $(L) $(INC)
Synthesis:	SynthesisFunctions.o Synthesis.o
	$(CC) -o Synthesis SynthesisFunctions.o Synthesis.o $(optspeed) $(L) $(INC) $(LIBR) 

install: all
	mkdir -p $(BIN)
	cp Analysis Synthesis $(BIN)

clean:
	rm -f *.o
	rm -f *~

distclean: clean
	rm -f Analysis Synthesis


