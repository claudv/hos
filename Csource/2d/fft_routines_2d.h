//
//  fft_routines.h
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "euler_2d.h"
#include <fftw3.h>


#ifndef euler_fft_routines_h
#define euler_fft_routines_h
#endif

/* FFT related global variables */
double          *f;
fftw_complex    *hf;
double          *ff;
fftw_complex    *hff;
fftw_plan       fftp;
fftw_plan       ifftp;
fftw_plan       fftp_large;
fftw_plan       ifftp_large;


void fft_2d(double* u, fftw_complex* hu, fftw_plan plan);
void ifft_2d(fftw_complex* hu, double* u, fftw_plan plan);
void fft_2d_large(double* u, fftw_complex* hu, fftw_plan plan);
void ifft_2d_large(fftw_complex* hu, double* u, fftw_plan plan);