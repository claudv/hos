//
//  hdf5_routines_2d.h
//  euler
//
//  Created by Claudio Viotti on 6/24/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include"euler_2d.h"

#ifndef euler_hdf5_routines_2d_h
#define euler_hdf5_routines_2d_h

#define DATASETNAME "ExtendibleArray"
#define RANK         2


#endif


void get_params(char* filename);
void get_ic_2d(char* filename, double* u1, double* u2);
hid_t create_file_2d(char* filename);
hid_t open_file_2d(char* filename);
herr_t close_file_2d(hid_t file);
void write_header_2d(hid_t file, double t);
void write_field_2d(hid_t file, double* eta, double* phi);
void write_extra_2d(hid_t file, double* Array1, double* Array2);