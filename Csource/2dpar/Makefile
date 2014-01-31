#CC=openmpicc
# Use this on Fionn
CC=mpiicc

CFLAGS=-c -O3
LDFLAGS=-lfftw3_mpi -lfftw3_threads -lfftw3 -lhdf5 -lm
SOURCES=main_2d_par.c model_2d_par.c time_schemes_par.c operators_par.c fft_routines_2d_par.c hdf5_routines_2d_par.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=euler_2d_par
LIBPATTERNS = lib%.a

#FFTW_INC = /usr/local/fftw/include
#FFTW_LIB = /usr/local/fftw/lib
#HDF5_INC = /usr/local/hdf5-1.8.11/hdf5/include
#HDF5_LIB = /usr/local/hdf5-1.8.11/hdf5/lib 

# Use this on Fionn
FFTW_INC = /ichec/home/users/viotti/fftw-3.3.3/include/
FFTW_LIB = /ichec/home/users/viotti/fftw-3.3.3/lib/
HDF5_INC = /ichec/home/users/viotti/hdf5-1.8.11/hdf5/include
HDF5_LIB = /ichec/home/users/viotti/hdf5-1.8.11/hdf5/lib/

INC = -I$(HDF5_INC) -I$(FFTW_INC)
LIB = -L$(HDF5_LIB) -L$(FFTW_LIB)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS) $(LIB)

.c.o: $(LDFLAGS)
	$(CC) $(CFLAGS) $< -o $@ $(INC) 

clean:
	rm *.o
