perga: srgaMain.o srga.o graph.o contig.o update.o  util.o hashPE.o \
	PEAssembly.o dtRowHash.o kmerseqBlock.o ridposBlock.o kmerBlock.o \
	reads.o align.o correction.o lenStatistics.o ref.o svm.o malign.o \
	msa.o
	gcc -g -o $@ $^ -lm  # $@:target; $^:dependences
	rm *.o
.c.o:  # .o files rely on the corresponding .c files
	gcc -g -c $< -lm

clean:
	rm *.o
	rm *~
	rm *.bak
