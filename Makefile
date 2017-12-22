# file: GNUMakefile
#

all: asc

asc: asc_0.o asc_1.o asc_2.o asc.o spectrogram_0.o spectrogram_1.o spectrogram_2.o kiss_fft.o asc_model_0.o asc_model_1.o asc_model_2.o asc_utils_0.o asc_utils_1.o asc_utils_2.o
	g++  -g asc_0.o asc_1.o asc_2.o asc.o spectrogram_0.o spectrogram_1.o spectrogram_2.o kiss_fft.o asc_model_0.o asc_model_1.o asc_model_2.o asc_utils_0.o asc_utils_1.o asc_utils_2.o -o asc

asc.o: asc.cc asc_0.cc asc.h
	g++ -g -c asc.cc -o asc.o

asc_0.o: asc_0.cc asc.h
	g++ -g -c asc_0.cc -o asc_0.o

asc_1.o: asc_1.cc asc.h
	g++ -g -c asc_1.cc -o asc_1.o

asc_2.o: asc_2.cc asc.h
	g++ -g -c asc_2.cc -o asc_2.o
	
	
spectrogram_0.o: spectrogram/spectrogram_0.cc spectrogram/spectrogram.h
	g++ -g -c spectrogram/spectrogram_0.cc -o spectrogram_0.o

spectrogram_1.o: spectrogram/spectrogram_1.cc spectrogram/spectrogram.h
	g++ -g -c spectrogram/spectrogram_1.cc -o spectrogram_1.o

spectrogram_2.o: spectrogram/spectrogram_2.cc  spectrogram/spectrogram.h 
	g++ -g -c spectrogram/spectrogram_2.cc   -o spectrogram_2.o

kiss_fft.o: spectrogram/kiss_fft130/kiss_fft.c spectrogram/kiss_fft130/kiss_fft.h
	g++ -g  -c -lm spectrogram/kiss_fft130/kiss_fft.c  -o kiss_fft.o

asc_model_0.o: asc_model_0.cc asc_model.h
	g++ -g -c asc_model_0.cc -o asc_model_0.o

asc_model_1.o: asc_model_1.cc asc_model.h
	g++ -g -c asc_model_1.cc -o asc_model_1.o

asc_model_2.o: asc_model_2.cc asc_model.h
	g++ -g -c asc_model_2.cc -o asc_model_2.o

asc_utils_0.o:  asc_utils/asc_utils_0.cc  asc_utils/asc_utils.h
	g++  -g -c asc_utils/asc_utils_0.cc  -o  asc_utils_0.o

asc_utils_1.o:  asc_utils/asc_utils_1.cc  asc_utils/asc_utils.h
	g++ -g -c asc_utils/asc_utils_1.cc  -o  asc_utils_1.o
	
asc_utils_2.o:  asc_utils/asc_utils_2.cc  asc_utils/asc_utils.h
	g++  -g -c asc_utils/asc_utils_2.cc  -o  asc_utils_2.o
	
clean:
	rm -f *.o 
	rm -f asc

#
# end of file


