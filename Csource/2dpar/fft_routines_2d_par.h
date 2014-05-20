//
//  fft_routines.h
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "euler_2d_par.h"


#ifndef euler_fft_routines_h
#define euler_fft_routines_h

/* FFT-related global variables */
double          *f;
fftw_complex    *hf;
fftw_plan       fftp;
fftw_plan       ifftp;

ptrdiff_t alloc_local;
ptrdiff_t local_Nx;
ptrdiff_t local_Nyhpo;
ptrdiff_t local_N;
ptrdiff_t local_0_start;
ptrdiff_t local_1_start;
ptrdiff_t i;
ptrdiff_t j;

ptrdiff_t fNx;
ptrdiff_t fNy;

void fft_2d(const double* u, fftw_complex* hu, fftw_plan plan);
void ifft_2d(const fftw_complex* hu, double* u, fftw_plan plan);
void fft_2d_large(double* u, fftw_complex* hu, fftw_plan plan);
void ifft_2d_large(fftw_complex* hu, double* u, fftw_plan plan);

#endif /*euler_fft_routines_h*/