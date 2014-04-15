//
//  operators.h
//  euler
//
//  Created by Claudio Viotti on 7/10/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "fft_routines_2d_par.h"


#ifndef euler_operators_h
#define euler_operators_h



#endif


void Dx(fftw_complex* hu, fftw_complex* hu_x);
void Dy(fftw_complex* hu, fftw_complex* hu_y);
void Dz(const fftw_complex* hu, fftw_complex* hu_z);
void Mult(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hprod);
void Sum(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hsum);
void Dealias(fftw_complex* hu);
void Filter(fftw_complex* hu, double k_peak, double fb1, double fb2);

