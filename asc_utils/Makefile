# file: GNUMakefile
#

all: asc_utils_test

asc_utils_test: asc_utils_0.o asc_utils_1.o asc_utils_2.o asc_utils_test.o
	g++  asc_utils_0.o asc_utils_1.o asc_utils_2.o asc_utils_test.o -o asc_utils_test

asc_utils_test.o: asc_utils_test.cc asc_utils_0.cc asc_utils.h
	g++  -c asc_utils_test.cc -o asc_utils_test.o

asc_utils_0.o: asc_utils_0.cc asc_utils.h
	g++  -c asc_utils_0.cc -o asc_utils_0.o

asc_utils_1.o: asc_utils_1.cc asc_utils.h
	g++  -c asc_utils_1.cc -o asc_utils_1.o

asc_utils_2.o: asc_utils_2.cc asc_utils.h
	g++ -c asc_utils_2.cc -o asc_utils_2.o

clean:
	rm -f asc_utils_test.o asc_utils_0.o asc_utils_1.o asc_utils_2.o
	rm -f asc_utils_test

#
# end of file


