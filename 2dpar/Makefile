CC=openmpicc
CFLAGS=-c -O3
LDFLAGS=-lfftw3_mpi -lfftw3_threads -lfftw3 -lhdf5 -lm
SOURCES=main_2d_par.c model_2d_par.c time_schemes_par.c operators_par.c fft_routines_2d_par.c hdf5_routines_2d_par.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=euler_2d_par
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
