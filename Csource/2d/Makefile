CC=openmpicc
CFLAGS=-c -O3 -std=gnu99
LDFLAGS=-lfftw3_threads -lfftw3 -lhdf5 -lm
SOURCES=main_2d.c model_2d.c time_schemes.c operators.c fft_routines_2d.c hdf5_routines_2d.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=euler_2d
LIBPATTERNS = lib%.a
INC = -I/usr/local/hdf5-1.8.11/hdf5/include -I/usr/local/fftw/include/
LIB = -L/usr/local/hdf5-1.8.11/hdf5/lib/ -L/usr/local/fftw/lib/

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS) $(LIB)

.c.o: $(LDFLAGS)
	$(CC) $(CFLAGS) $< -o $@ $(INC) 

clean:
	rm *.o
