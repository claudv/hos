//
//  euler_2d.h
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#ifndef EULER_H
#define EULER_H
 

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<assert.h>
#include<tgmath.h>
#include<fftw3-mpi.h>

#include "mpi.h"
#include "hdf5.h"

#define H5_DATATYPE "H5T_IEEE_F64LE"
#define SAVE_FILE_BUFSIZE 128
#define RUNTIME_DATA_BUFSIZE 32

/*----------------------------------------------------------------------------------------------------------*/
/* Choose between the standard untransposed (FFT_TRANSPOSE = 0) or the faster transposed (FFT_TRANSPOSE = 1)*/
/* fftw output format. The Fourier coefficients are distributed differently across the parallel processes   */
/*                                                                                                          */
/*  REGULAR OUTPUT (row-major ordered in memory):                                                           */
/*                                                                                                          */
/*  0 1 ...                                  ... Ny/2   \                                                   */
/*  .                                               .    |                                                  */
/*  .                                                    |  process 0                                       */
/*  local_Nx - 1                                         |                                                  */
/*  .                                               .    |     | padding (internal fftw use only)           */
/*  .                                               .   /      |                                            */
/*  local_0_start                                       \                                                   */
/*  .                                                    |                                                  */
/*  .                                                    |  process 1                                       */
/*  local_0_start + local_Nx - 1                         |                                                  */
/*  .                                                    |     | padding (internal fftw use only)           */
/*  .                                                   /      |                                            */
/*                                                                                                          */
/*  .                                                          .                                            */
/*  .                                                          .                                            */
/*                                                                                                          */
/*  (Nproc-1)*local_0_start  ... Ny/2                   \                                                   */
/*  .                                               .    |                                                  */
/*  .                                               .    |  process Nproc-1                                 */
/*  (Nproc-1)*local_0_start + local_Nx - 1   ... Ny/2    |                                                  */
/*  .                                                    |     | padding (internal fftw use only)           */
/*  .                                                   /      |                                            */
/*                                                                                                          */
/*  TRANSPOSED OUTPUT (row-major ordered in memory):                                                        */
/*                                                                                                          */
/*  0 1 ...                                    ... Nx-1 \                                                   */
/*  .                                                 .  |                                                  */
/*  .                                                    |  process 0                                       */
/*  local_Nyhpo - 1                                      |                                                  */
/*  .                                                 .  |     | padding (internal fftw use only)           */
/*  .                                                 . /      |                                            */
/*  local_1_start                                     . \                                                   */
/*  .                                                    |                                                  */
/*  .                                                    |  process 1                                       */
/*  local_1_start + local_Nyhpo - 1                      |                                                  */
/*  .                                                    |     | padding (internal fftw use only)           */
/*  .                                                   /      |                                            */
/*                                                                                                          */
/*  .                                                          .                                            */
/*  .                                                          .                                            */
/*                                                                                                          */
/*  (Nproc-1)*local_1_start                           . \                                                   */
/*  .                                                 .  |                                                  */
/*  .                                                    |  process Nproc-1                                 */
/*  (Nproc-1)*local_1_start + local_Nyhpo - 1  ... Nx-1  |                                                  */
/*  .                                                    |   | padding (internal fftw use only)             */
/*  .                                                   /    |                                              */
#define FFT_TRANSPOSE 1


/*---------------------------------------------------------------------------------------------------------*/
/* Number of levels used in the HOS expansion.                                                             */
/* Particular cases are:                                                                                   */
/* Nlev=0 --> linear waves;                                                                                */
/* Nlev=2 --> Zakharov equations.                                                                          */
#define NLevs 2

typedef unsigned char flg_type;

/*--------------------------------------*/
/* Global variables --------------------*/

/* 
In the reader routine global variables are initialized,
therefore they are included as non-constants by defining
global_params_reader before header inclusion. They should 
be constant anywhere else.
THIS APPROACH COULD BE UNSAFE--UNDER TESTING 
*/
#ifdef global_params_reader
#define __TYPE__QUAL__ extern
#else
#define __TYPE__QUAL__ const
#endif

__TYPE__QUAL__ int          Nx;
__TYPE__QUAL__ int          Ny;
__TYPE__QUAL__ int          runsubid;
__TYPE__QUAL__ double       Lx;
__TYPE__QUAL__ double       Ly;
__TYPE__QUAL__ double       Kx0;
__TYPE__QUAL__ double       Ky0;
__TYPE__QUAL__ double       g;

__TYPE__QUAL__ double       T;
__TYPE__QUAL__ double       dtsave;
__TYPE__QUAL__ flg_type     saveflg;    // =1 -> basic output, >1 -> extended output

/* Dommermuth ramping -----------------*/
__TYPE__QUAL__ flg_type    rampflg;
__TYPE__QUAL__ double      Tramp;

/* Wind forcing -----------------------*/
__TYPE__QUAL__ flg_type    windflg;
__TYPE__QUAL__ double      Uwind_x;
__TYPE__QUAL__ double      Uwind_y;

/* Size of dealiased complex arrays ---*/
__TYPE__QUAL__ ptrdiff_t   mx;
__TYPE__QUAL__ ptrdiff_t   my;


fftw_complex*   hetan;
fftw_complex*   hphin;

/* MPI related variables --------------*/
int         mpi_size;
int         mpi_rank;
MPI_Comm    comm;
MPI_Info    info;

/* Global arrays for temporary storage */
double*         temp1;
double*         temp2;
fftw_complex*   htemp1;
fftw_complex*   htemp2;

/* Runtime datafile -----------------*/
char        runtime_data_buff[RUNTIME_DATA_BUFSIZE];
FILE*       runtime_fid;
double      Ham;
double      Ham_glob;


#endif /* EULER_H */
