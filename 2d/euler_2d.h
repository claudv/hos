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
#include <fftw3.h>

#include "hdf5.h"

#define H5_DATATYPE "H5T_IEEE_F64LE"
#define SAVE_FILE_BUFSIZE 128

///////////////////////
/* Global variables. */


/* Parameters */
int             Nx;
int             Ny;
int             Nxl;
int             Nyl;
int             N;
int             runsubid;
double          Lx;
double          Ly;
double          Kx0;
double          Ky0;
double          g;


double          T;
double          dtsave;

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
fftw_complex    *funElem;
fftw_complex    *htemp;
fftw_complex    **fun;


/* Working space for assembling rhs and dealiasing */
double*         eta_x;
double*         phi_x;

double*         temp1;
double*         temp2;
double*         temp1_large;
double*         temp2_large;

fftw_complex*   htemp1;
fftw_complex*   htemp2;
fftw_complex*   htemp1_large;
fftw_complex*   htemp2_large;
