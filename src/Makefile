perga: pergaMain.o perga.o graph.o contig.o update.o  util.o hashPE.o \
	PEAssembly.o dtRowHash.o kmerseqBlock.o ridposBlock.o kmerBlock.o \
	reads.o lenStatistics.o svm.o \
	contiggraph.o candPath.o tandPath.o contigPath.o merge.o \
	scafPerga.o scafContigIndex.o scafMap.o scafLink.o scafOverlap.o \
	scafGap.o scafSeq.o
	gcc -g -o $@ $^ -lm  # $@:target; $^:dependences
	rm *.o
.c.o:  # .o files rely on the corresponding .c files
	gcc -g -c $< -lm

clean:
	rm *.o
	rm *~
	rm *.bak
