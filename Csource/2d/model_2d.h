//
//  model_2d.h
//  euler
//
//  Created by Claudio Viotti on 6/24/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//


#include "euler_2d.h"
#include "fft_routines_2d.h"
#include "operators.h"

#ifndef euler_model_h
#define euler_model_h



#endif

void rhs_test(fftw_complex* rhs, fftw_complex* u);
void rhs_hos(fftw_complex* rhs, fftw_complex* u, double t);
void rhs_hos_setup();

void Zvel(fftw_complex* hu,  fftw_complex* hwM, fftw_complex* hwM2, fftw_complex* hw2M, fftw_complex* hw2M2, double t);


int             NLevs;
double*         Coeff;
fftw_complex*   hetan;
fftw_complex*   hphin;
fftw_complex*   hwn;
fftw_complex*   hwM;
fftw_complex*   hwM2;
fftw_complex*   hw2M;
fftw_complex*   hw2M2;
