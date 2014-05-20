//
//  hdf5_routines_2d.h
//  euler
//
//  Created by Claudio Viotti on 6/24/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//


#ifndef euler_hdf5_routines_2d_h
#define euler_hdf5_routines_2d_h

#include "euler_2d_par.h"

#define DATASETNAME "ExtendibleArray"
#define RANK         2

void get_ic_2d(char* filename, double* u1, double* u2);

/* These functions simply wrap the basic HDF5 file manipulation routines. */
hid_t create_file_2d(char* filename);
hid_t open_file_2d(char* filename);
herr_t close_file_2d(hid_t file);

/* These functions write the output HDF5 data files. */
void write_header_2d(hid_t file, double time);
void write_field_2d(hid_t file, double* eta, double* phi);
void write_field_complex_2d(hid_t h5_file, const fftw_complex* heta, const fftw_complex* phi);
void write_extra_2d(hid_t file, double* array1, double* array2);

#endif
