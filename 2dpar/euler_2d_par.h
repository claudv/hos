//
//  euler_2d.h
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//


#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<tgmath.h>
#include<fftw3-mpi.h>

#include "mpi.h"
#include "hdf5.h"

#define H5_DATATYPE "H5T_IEEE_F64LE"
#define SAVE_FILE_BUFSIZE 128
#define RUNTIME_DATA_BUFSIZE 32

#define FFT_TRANSPOSE 1

///////////////////////
/* Global variables. */


/* Parameters */
int             Nx;
int             Ny;
//ptrdiff_t       N;
int             runsubid;
double          Lx;
double          Ly;
double          Kx0;
double          Ky0;
double          g;


double          T;
double          dtsave;
int             saveflg;    // =1 -> basic output, >1 -> extra output

/* Dommermuth ramping */
int             rampflg;
double          Tramp;

/* MPI variables */
int         mpi_size;
int         mpi_rank;
MPI_Comm    comm;
MPI_Info    info;


/* Definitions global for debugging pourpose */
hid_t       savefileid;
hid_t       savefileid2;
char        savefile_buff[SAVE_FILE_BUFSIZE];
char        savefile2_buff[SAVE_FILE_BUFSIZE];
herr_t      status;

/* Coeffs for RK routines */
int             Nstages;
double*         a;
double*         bElem;
double*         c;
double*         cs;
double**        b;
fftw_complex    *htemp;
fftw_complex    *fun;


/* Working space for assembling rhs and dealiasing */
double*         eta_x;
double*         phi_x;

double*         temp1;
double*         temp2;

fftw_complex*   htemp1;
fftw_complex*   htemp2;


/* Variables used in computing vertical velocity */
int             NLevs;
double*         Coeff;
fftw_complex*   hetan;
fftw_complex*   hphin;
fftw_complex*   hwn;
fftw_complex*   hwM;
fftw_complex*   hwM2;
fftw_complex*   hw2M;
fftw_complex*   hw2M2;


/* Runtime datafile */
char        runtime_data_buff[RUNTIME_DATA_BUFSIZE];
FILE*       runtime_fid;
double       Ham;
double       Ham_glob;

