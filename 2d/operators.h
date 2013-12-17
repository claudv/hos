//
//  operators.h
//  euler
//
//  Created by Claudio Viotti on 7/10/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "fft_routines_2d.h"

#ifndef euler_operators_h
#define euler_operators_h



#endif


void Dx(fftw_complex* hu, fftw_complex* hu_x);
void Dy(fftw_complex* hu, fftw_complex* hu_y);
void Dz(fftw_complex* hu, fftw_complex* hu_z);
void Extend(fftw_complex* hu, fftw_complex* hu_large);
void Shrink(fftw_complex* hu_large, fftw_complex* hu);
void MultDea(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hprod);
void Mult(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hprod);
void Sum(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hsum);


