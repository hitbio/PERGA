#!/bin/sh
	mkdir bin
	cd src/PERGA
	make
	mv -f perga ../../bin
	cp -r model ../../bin
	cd ../../src/PERGA_SF
	make
	mv -f perga_sf ../../bin
	cd ../../
	
