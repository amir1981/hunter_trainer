# file: GNUMakefile
#

all: spectrogram_test 

spectrogram_test: spectrogram_0.o spectrogram_1.o spectrogram_2.o kiss_fft.o spectrogram_test.o
	g++ spectrogram_0.o spectrogram_1.o spectrogram_2.o spectrogram_test.o kiss_fft.o -o spectrogram_test

spectrogram_test.o: spectrogram_test.cc spectrogram_0.cc spectrogram.h
	g++ -c spectrogram_test.cc -o spectrogram_test.o

spectrogram_0.o: spectrogram_0.cc spectrogram.h
	g++ -c spectrogram_0.cc -o spectrogram_0.o

spectrogram_1.o: spectrogram_1.cc spectrogram.h
	g++ -c spectrogram_1.cc -o spectrogram_1.o

spectrogram_2.o: spectrogram_2.cc  spectrogram.h 
	g++ -c spectrogram_2.cc   -o spectrogram_2.o

kiss_fft.o: kiss_fft130/kiss_fft.c kiss_fft130/kiss_fft.h
	g++  -c -lm kiss_fft130/kiss_fft.c  -o kiss_fft.o

clean:
	rm -f spectrogram_test.o spectrogram_0.o spectrogram_1.o spectrogram_2.o
	rm -f spectrogram_test 
	rm -f kiss_fft.o

#
# end of file


