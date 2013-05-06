CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS = 
GLIBS += 
OBJECTS = trimSide.o 
HEADERS = globalConstants.h

ALL : trimSide.exe
	@echo "Listo!"

trimSide.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o trimSide.exe $(LIBS) $(GLIBS) $(CFLAGS)

trimSide.o : trimSide.cc $(HEADERS)
	$(CPP) -c trimSide.cc -o trimSide.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
