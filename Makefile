CC=g++
CFLAGS=-g -Wall -fPIC -O3
LIB=-litpp -lgsl -ljpeg -ltiff
DSPOBJS=distortion.o quantizer.o
AUOBJS=pqmf.o wav.o mpeg_binary.o mpeg_decoder.o mpeg_file.o mpeg_frame.o bit_stream.o
COMMOBJS=myldpc.o capacity.o mlc_msd.o convolutional_code.o length_adjuster.o turbo_code.o mymodulation.o
IMGOBJS=stbi_load.o stbi_write.o mybmp.o mytiff.o myentropy.o myjpeg.o
# UTL=myutl.o
INCLUDEPATH=./include/
SRCPATH=./SRC/


all: libmylib.so

libmylib.so: $(AUOBJS) $(IMGOBJS) $(COMMOBJS) $(DSPOBJS)
	$(CC) $(CFLAGS) -shared -o libmylib.so $(DSPOBJS) $(AUOBJS) $(IMGOBJS) $(COMMOBJS) $(LIB)

pqmf.o: $(INCLUDEPATH)pqmf.h $(SRCPATH)pqmf.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)pqmf.cpp

quantizer.o: $(INCLUDEPATH)quantizer.h $(SRCPATH)quantizer.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)quantizer.cpp

distortion.o: $(INCLUDEPATH)distortion.h $(SRCPATH)distortion.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)distortion.cpp

turbo_code.o: $(INCLUDEPATH)turbo_code.h $(SRCPATH)turbo_code.cpp $(INCLUDEPATH)mycomm.h
	$(CC) $(CFLAGS) -c $(SRCPATH)turbo_code.cpp

myldpc.o: $(INCLUDEPATH)myldpc.h $(SRCPATH)myldpc.cpp $(INCLUDEPATH)mycomm.h
	$(CC) $(CFLAGS) -c $(SRCPATH)myldpc.cpp -litpp

capacity.o: $(INCLUDEPATH)capacity.h $(SRCPATH)capacity.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)capacity.cpp -litpp

mlc_msd.o: $(INCLUDEPATH)mlc_msd.h $(SRCPATH)mlc_msd.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)mlc_msd.cpp -litpp

mymodulation.o: $(INCLUDEPATH)mymodulation.h $(SRCPATH)mymodulation.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)mymodulation.cpp -litpp

length_adjuster.o: $(INCLUDEPATH)length_adjuster.h $(SRCPATH)length_adjuster.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)length_adjuster.cpp -litpp

convolutional_code.o: $(INCLUDEPATH)convolutional_code.h $(SRCPATH)convolutional_code.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)convolutional_code.cpp -litpp

wav.o: $(INCLUDEPATH)wav.h $(SRCPATH)wav.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)wav.cpp

mpeg_binary.o: $(INCLUDEPATH)mpeg.h $(INCLUDEPATH)mpeg_binary.h $(SRCPATH)mpeg_binary.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)mpeg_binary.cpp

mpeg_decoder.o: $(INCLUDEPATH)bit_stream.h $(INCLUDEPATH)mpeg.h $(SRCPATH)mpeg_decoder.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)mpeg_decoder.cpp

mpeg_file.o: $(INCLUDEPATH)bit_stream.h $(INCLUDEPATH)mpeg.h $(SRCPATH)mpeg_file.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)mpeg_file.cpp

mpeg_frame.o: $(INCLUDEPATH)bit_stream.h $(INCLUDEPATH)mpeg.h $(SRCPATH)mpeg_frame.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)mpeg_frame.cpp

bit_stream.o: $(INCLUDEPATH)bit_stream.h $(SRCPATH)bit_stream.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)bit_stream.cpp

stbi_load.o: $(INCLUDEPATH)stbi_load.h $(SRCPATH)stbi_load.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)stbi_load.cpp

stbi_write.o: $(INCLUDEPATH)stbi_write.h $(SRCPATH)stbi_write.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)stbi_write.cpp

mybmp.o: $(INCLUDEPATH)mybmp.h $(SRCPATH)mybmp.cpp $(INCLUDEPATH)mymatrix.h
	$(CC) $(CFLAGS) -c $(SRCPATH)mybmp.cpp

mytiff.o: $(INCLUDEPATH)mytiff.h $(SRCPATH)mytiff.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)mytiff.cpp -ltiff

myjpeg.o: $(INCLUDEPATH)myjpeg.h $(SRCPATH)myjpeg.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)myjpeg.cpp -ljpeg -litpp

myentropy.o: $(INCLUDEPATH)myentropy.h $(SRCPATH)myentropy.cpp
	$(CC) $(CFLAGS) -c $(SRCPATH)myentropy.cpp -litpp

clean:
	rm $(AUOBJS) $(COMMOBJS) $(IMGOBJS) libmylib.so
