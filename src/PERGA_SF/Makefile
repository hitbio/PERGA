perga_sf: scafMain.o scafSRGA.o scafContig.o scafContigIndex.o scafMapping.o \
	scafReadList.o scafContigList.o scafContigLinking.o scafContigOverlap.o \
	scafGraph.o scafGapFilling.o scafUpdate.o scafPEAssembly.o scafContigMerge.o \
	util.o
	gcc -std=c99 -g -o $@ $^ -lm  # $@: targets; $^: dependences
	rm *.o
.c.o:  # .o relys on the corresponding .c.
	gcc -std=c99 -g -c $<

clean:
	rm *.o
	rm *~
	rm *.bak
