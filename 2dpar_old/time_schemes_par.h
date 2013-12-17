//
//  time_schemes.h
//  euler
//
//  Created by Claudio Viotti on 3/9/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//


#include "euler_2d_par.h"
#include "fft_routines_2d_par.h"
#include "model_2d_par.h"



#ifndef euler_time_schemes_h
    #define euler_time_schemes_h
#endif


void sol_update_RK(fftw_complex* u_old,double* t,double dt,char* dtflag);
void Setup_TimeScheme(int scheme_flg);